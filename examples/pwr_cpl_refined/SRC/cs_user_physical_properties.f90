!-------------------------------------------------------------------------------

!                      code_saturne version 8.0
!                      ------------------------
! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
!> \file cs_user_physical_properties.f90
!>
!> \brief Definition of physical variable laws.
!>
!> See \ref physical_properties for examples.
!>

!===============================================================================
!> \brief Definition of physical variable laws.
!>
!> \section Warning
!>
!> It is \b forbidden to modify turbulent viscosity \c visct here
!> (a specific subroutine is dedicated to that: \ref usvist)
!>
!> - icp = 0 must <b> have been specified </b>
!>    in \ref usipsu if we wish to define a variable specific heat
!>    cpro_cp (otherwise: memory overwrite).
!>
!> - the kivisl field integer key (diffusivity_id)
!>    must <b> have been specified </b>
!>    in \ref usipsu if we wish to define a variable viscosity
!>    \c viscls.
!>
!>
!> \remarks
!>  - This routine is called at the beginning of each time step
!>    Thus, <b> AT THE FIRST TIME STEP </b> (non-restart case), the only
!>    values initialized before this call are those defined
!>      - in the GUI or  \ref usipsu (cs_user_parameters.f90)
!>             - density    (initialized at \c ro0)
!>             - viscosity  (initialized at \c viscl0)
!>      - in the GUI or \ref cs_user_initialization
!>             - calculation variables (initialized at 0 by defaut
!>             or to the value given in the GUI or in \ref cs_user_initialization)
!>
!>  - We may define here variation laws for cell properties, for:
!>     - density:                                    rom    kg/m3
!>     - density at boundary faces:                  romb   kg/m3)
!>     - molecular viscosity:                        cpro_viscl  kg/(m s)
!>     - specific heat:                              cpro_cp     J/(kg degrees)
!>     - diffusivities associated with scalars:      cpro_vscalt kg/(m s)
!>
!> \b Warning: if the scalar is the temperature, cpro_vscalt corresponds
!> to its conductivity (Lambda) in W/(m K)
!>
!>
!> The types of boundary faces at the previous time step are available
!>   (except at the first time step, where arrays \c itypfb and \c itrifb have
!>   not been initialized yet)
!>
!> It is recommended to keep only the minimum necessary in this file
!>   (i.e. remove all unused example code)
!>
!>
!> \section usphyv_cell_id Cells identification
!>
!> Cells may be identified using the \ref getcel subroutine.
!> The syntax of this subroutine is described in the
!> \ref cs_user_boundary_conditions subroutine,
!> but a more thorough description can be found in the user guide.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     mbrom         indicator of filling of romb array
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine usphyv &
 ( nvar   , nscal  ,                                              &
   mbrom  ,                                                       &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use field
use mesh
use lagran
use cs_c_bindings
use coupling_resolved

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mbrom

double precision dt(ncelet)

!===============================================================================
! Local variables

integer          ivart, iel, ifac
integer          iscal, ifcvsl
integer          cpro_index
double precision intp1, intp2, intp
double precision xvart
double precision factor
double precision maxlim, minlim

double precision, dimension(:), pointer :: coefap, coefbp

double precision, dimension(:), pointer :: cpro_rom, bfpro_rom
double precision, dimension(:), pointer :: cpro_viscl
double precision, dimension(:), pointer :: cpro_cp
double precision, dimension(:), pointer :: cpro_vscalt

double precision, dimension(:), pointer :: cvar_scalt
double precision, dimension(:), pointer :: temp

save maxlim, minlim

!===============================================================================

!Use uniformly distributed physical property table
if (itherm.eq.1) then
  ! Temperature, e.g. delta_t = 0.1C
  factor = 1.d0/increment
else if (itherm.eq.2) then
  ! Entalpy, e.g. delta_h = 1000J/kg 
  factor = 1.d0/increment
  call field_get_val_s_by_name('temperature', temp)
end if

if(physical_property_table) then
  if(ntcabs==ntpabs) call readin_properties_from_data_table(maxlim, minlim)
end if

if (iscalt.gt.0) then
  ivart = isca(iscalt)
  call field_get_val_s(ivarfl(ivart), cvar_scalt)
endif

!!Convert enthalpy to temperature
if (itherm.eq.2) then
  do iel = 1, ncel
    xvart = cvar_scalt(iel)
    if (xvart>maxlim) xvart = maxlim
    if (xvart<minlim) xvart = minlim
    cpro_index = int((xvart-user_prop(1,1))*factor)+1
    intp = user_prop(cpro_index+1,1)-user_prop(cpro_index,1)
    intp1 = (xvart-user_prop(cpro_index,1))/intp
    intp2 = 1.d0-intp1
    temp(iel) = user_prop(cpro_index, 2)*intp2&
              + user_prop(cpro_index+1, 2)*intp1
  end do
end if

!-------------------------------------------------------------------------------
!! Variable density
if (irovar.gt.0) then
  call field_get_val_s(icrom, cpro_rom)
  call field_get_val_s(ibrom, bfpro_rom)
else
  cpro_rom => NULL()
endif

if (irovar.gt.0 .and. irangp.le.0 .and. ntcabs.eq.ntpabs) then
  write(*,*)'User defined variable density is used'
end if

!! Stop if the density has not been defined as variable
!if (irovar.le.0) then
!  write(nfecra,3000) irovar
!  call csexit (1)
!endif

if (irovar.gt.0) then
  do iel = 1, ncel
    xvart = cvar_scalt(iel)
    if (xvart>maxlim) xvart = maxlim
    if (xvart<minlim) xvart = minlim
    cpro_index = int((xvart-user_prop(1,1))*factor)+1
    intp = user_prop(cpro_index+1,1)-user_prop(cpro_index,1)
    intp1 = (xvart-user_prop(cpro_index,1))/intp
    intp2 = 1.d0-intp1
    cpro_rom(iel) = user_prop(cpro_index, 3)*intp2&
                  + user_prop(cpro_index+1, 3)*intp1
  end do

!! Density at boundary faces
!! For coarse grid, boundary density should not be calculated using the following
!! method, as the wall function does not work.
  if(.false.) then
    call field_get_coefa_s(ivarfl(ivart), coefap)
    call field_get_coefb_s(ivarfl(ivart), coefbp)

    mbrom = 1

    do ifac = 1, nfabor
      iel = ifabor(ifac)
      xvart = coefap(ifac) + cvar_scalt(iel)*coefbp(ifac)
      if (xvart>maxlim) xvart = maxlim
      if (xvart<minlim) xvart = minlim
      cpro_index = int((xvart-user_prop(1,1))*factor)+1
      intp = user_prop(cpro_index+1,1)-user_prop(cpro_index,1)
      intp1 = (xvart-user_prop(cpro_index,1))/intp
      intp2 = 1.d0-intp1
      bfpro_rom(ifac) = user_prop(cpro_index, 3)*intp2 &
                      + user_prop(cpro_index+1, 3)*intp1
    enddo
  end if

end if

!-------------------------------------------------------------------------------
!! Variable molecular viscosity
if (ivivar.gt.0) then
  call field_get_val_s(iprpfl(iviscl), cpro_viscl)
else
  cpro_viscl => NULL()
endif

if (ivivar.gt.0 .and. irangp.le.0 .and. ntcabs.eq.ntpabs) then
  write(*,*)'User defined variable molecular viscosity is used'
end if

!! Stop if the volumetric viscosity has not been defined as variable
!if (ivivar.le.0) then
!  write(nfecra,2000) ivivar
!  call csexit (1)
!endif

if (ivivar.gt.0) then
  do iel = 1, ncel
    xvart = cvar_scalt(iel)
    cpro_index = int((xvart-user_prop(1,1))*factor)+1
    intp = user_prop(cpro_index+1,1)-user_prop(cpro_index,1)
    intp1 = (xvart-user_prop(cpro_index,1))/intp
    intp2 = 1.d0-intp1
    cpro_viscl(iel) = user_prop(cpro_index,5)*intp2 &
                    + user_prop(cpro_index+1,5)*intp1
  end do
end if

!-------------------------------------------------------------------------------
!! In case not solving energy
if (iscalt.le.0) return

!-------------------------------------------------------------------------------
!! Variable isobaric specific heat (heat capacity)
if (icp.ge.0) then
  call field_get_val_s(iprpfl(icp), cpro_cp)
else
  cpro_cp => NULL()
end if

if (icp.ge.0 .and. irangp.le.0 .and. ntcabs.eq.ntpabs) then
  write(*,*)'User defined variable specific heat is used'
end if

!! Stop if the isobaric specific heat cpro_cp has not been defined as variable
!if (icp.lt.0) then
!  write(nfecra,1000) icp
!  call csexit (1)
!endif

if (icp.ge.0) then
  do iel = 1, ncel
    xvart = cvar_scalt(iel)
    if (xvart>maxlim) xvart = maxlim
    if (xvart<minlim) xvart = minlim
    cpro_index = int((xvart-user_prop(1,1))*factor)+1
    intp = user_prop(cpro_index+1,1)-user_prop(cpro_index,1)
    intp1 = (xvart-user_prop(cpro_index,1))/intp
    intp2 = 1.d0-intp1
    cpro_cp(iel) = user_prop(cpro_index,4)*intp2 &
                 + user_prop(cpro_index+1,4)*intp1
  end do
end if

!-------------------------------------------------------------------------------
!! Variable molecular thermal conductivity

call field_get_key_int(ivarfl(isca(iscalt)), kivisl, ifcvsl)

if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_vscalt)
else
  cpro_vscalt => NULL()
endif

! Stop if the molecular thermal conductivity has not been defined as variable

if (ifcvsl.ge.0 .and. irangp.le.0 .and. ntcabs.eq.ntpabs) then
  write(*,*)'User defined variable molecular thermal conductivity is used'
end if
!if (ifcvsl.lt.0) then
!  write(nfecra,1010) iscalt
!  call csexit (1)
!endif

if (ifcvsl.ge.0) then
  ! if iscalt is not temperature
  if (itherm.ne.1) then
    !variable Cp
    if (icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)
    !constant Cp
    if (icp.le.0) then
      do iel = 1, ncel
        xvart = cvar_scalt(iel)
        cpro_cp(iel) = cp0
      enddo
    else
      ! --- If Cp is not uniform, we use cpro_cp above
      do iel = 1, ncel
        xvart = cvar_scalt(iel)
        if (xvart>maxlim) xvart = maxlim
        if (xvart<minlim) xvart = minlim
        cpro_index = int((xvart-user_prop(1,1))*factor)+1
        intp = user_prop(cpro_index+1,1)-user_prop(cpro_index,1)
        intp1 = (xvart-user_prop(cpro_index,1))/intp
        intp2 = 1.d0-intp1      
        cpro_vscalt(iel) = user_prop(cpro_index, 6)*intp2 &
                         + user_prop(cpro_index+1, 6)*intp1
        cpro_vscalt(iel) = cpro_vscalt(iel)/cpro_cp(iel)
      enddo
    endif
  ! if iscalt is temperature, the Cp division is not needed
  else
    do iel = 1, ncel
      xvart = cvar_scalt(iel)
      if (xvart>maxlim) xvart = maxlim
      if (xvart<minlim) xvart = minlim
      cpro_index = int((xvart-user_prop(1,1))*factor)+1
      intp = user_prop(cpro_index+1,1)-user_prop(cpro_index,1)
      intp1 = (xvart-user_prop(cpro_index,1))/intp
      intp2 = 1.d0-intp1      
      cpro_vscalt(iel) = user_prop(cpro_index, 6)*intp2 &
                       + user_prop(cpro_index+1, 6)*intp1
    enddo
  endif

end if

!--------
! Formats
!--------
 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@       in the GUI or in the user subroutine ''usipsu'', the',/, &
'@         isobaric specific heat is declared as a property',/,   &
'@         uniform in space: icp = ',i10   ,/,                    &
'@       in the user subroutine ''usphyv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the GUI input data or the',/, &
'@    user subroutines ''usipsu'' or ''usphyv''.',/,              &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@     For the scalar ',i10,/,                                    &
'@       in the GUI or in the user subroutine ''usipsu'', the',/, &
'@         thermal conductivity is declared as a property',/,     &
'@         uniform in space.',/,                                  &
'@       in the user subroutine ''usphyv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the GUI input data or the',/, &
'@    user subroutines ''usipsu'' or ''usphyv''.',/,              &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

 2000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@       in the user subroutine ''uscfx2'', the molecular',/,     &
'@         volumetric viscosity is declared as a property',/,     &
'@         uniform in space: iviscv = ',i10,/,                    &
'@       in the user subroutine ''usphyv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the user subroutines',/,      &
'@    ''uscfx2'' or ''usphyv''.',/,                               &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

 3000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@       in the user subroutine ''uscfx2'', the density',/,       &
'@         is declared as a property',/,                          &
'@         uniform in space: irovar = ',i10,/,                    &
'@       in the user subroutine ''usphyv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the user subroutines',/,      &
'@    ''uscfx2'' or ''usphyv''.',/,                               &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
!----
! End
!----

return
end subroutine usphyv

!===============================================================================
!> \brief Modify turbulent viscosity
!>
!> This subroutine is called at beginning of each time step
!> after the computation of the turbulent viscosity
!> (physical quantities have already been computed in \ref usphyv).
!>
!> Turbulent viscosity \f$ \mu_T \f$ (kg/(m s)) can be modified.
!>
!> A modification of the turbulent viscosity can lead to very
!> significant differences between solutions and even give wrong
!> results.
!>
!> This subroutine is therefore reserved to expert users.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icepdc        head loss cell numbering
!> \param[in]     icetsm        numbering of cells with mass source term
!> \param[in]     itypsm        kind of mass source for each variable
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        work array for head loss terms
!> \param[in]     smacel        values of variables related to mass source
!>                              term. If ivar=ipr, smacel=mass flux
!_______________________________________________________________________________

subroutine usvist &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)

!===============================================================================

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine usvist

!===============================================================================
