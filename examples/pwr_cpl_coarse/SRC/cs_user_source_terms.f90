!===============================================================================
! User source terms definition.
!
! 1) Momentum equation (coupled solver)
! 2) Species transport
! 3) Turbulence (Rij-epsilon)
!===============================================================================

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
!> \file cs_user_source_terms.f90
!>
!> \brief User subroutines for additional right-hand side source terms
!>
!> See \ref cs_user_source_terms and
!> \ref cs_user_source_terms-scalar_in_a_channel for examples.
!===============================================================================

!===============================================================================
!> \brief Additional right-hand side source terms for velocity components equation
!> (Navier-Stokes)
!>
!> \deprecated Use \ref cs_user_source_terms instead.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     ivar          index number of the current variable
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!>                              source terms or mass rate
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!_______________________________________________________________________________

subroutine ustsnv &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use mesh
use field
use SubChCFD

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          ivar

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(3,ncelet), crvimp(3,3,ncelet)

! Local variables

character*80    :: chaine, chaine1

integer         :: i, iel, ifac, ielt, nlfac, iflmas
integer         :: nlelt, cptfac

double precision :: tau_w, dpdzf, dpdzd, Kc, len_porous, ubulk

integer, allocatable, dimension(:) :: lstelt, lstfac
double precision, dimension(:), pointer ::  Re, cf, den, Dh, volume_porosi, &
                                            axial_resis_f, axial_resis_d, &
                                            porous_zone, imasfl, wall_shear_stress
double precision, dimension(:,:),pointer :: A_vel, vel

double precision :: flowint, flowref, surf, value, flowrate
double precision :: flown, flown1, dp, rho_ref
save flown, flown1, dp
!===============================================================================


!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(nfabor))

chaine = all_walls

!! Use subchannel velocity to calculate the frictional force
call field_get_val_v(ivarfl(iu),vel)
call field_get_val_v_by_name('AverageVelocity', A_vel)
call subchannel_velocity()

!! Calculate filtered velocity for the porous media zones
if(prs .and. iporo_cal==3) call filtered_velocity()

!!Calculate subchannel Reynolds number 
call field_get_val_s_by_name('AxialReynoldsNumber', Re)
call subchannel_axial_Re()

call field_get_val_s_by_name('AverageDensity', den)

!! Calculate the frictional factor
call field_get_val_s_by_name('AxialFrictionFactor', cf)
call subchannel_axial_f()

if(prs) then
  call field_get_val_s_by_name('HydraulicDiameter', Dh)
  call field_get_val_s_by_name('porosity', volume_porosi)
  call field_get_val_s_by_name('AxialResistanceFriction', axial_resis_f)
  call field_get_val_s_by_name('AxialResistanceDrag', axial_resis_d)
  call field_get_val_s_by_name('PorousZone', porous_zone)
!  call field_get_val_s_by_name('VolumePorosity', volume_porosi)
end if

!! If mometum equation is not solved, then return. Note, to specify the source 
!! term for the momentum equation, this routine is only called for ivar=iu.
if (ivar.ne.iu) return 

!!===========================
!!Source term calculation
!!===========================
call field_get_val_s_by_name('WallShearStress', wall_shear_stress)
call getfbr(chaine, nlelt, lstelt)

!! Impose the momentum source of wall frictional force
do ielt = 1, nlelt
  ifac = lstelt(ielt)
  iel = ifabor(ifac)

  tau_w = 0.25d0*cf(iel)*0.5d0*den(iel)*abs(A_vel(3,iel))*A_vel(3,iel)
  wall_shear_stress(iel) = abs(tau_w)
!  crvexp(3,iel) = -tau_w*surfbn(ifac)+crvexp(3,iel)

!! Implicit treatment may increase numerical stability
  crvimp(3,3,iel) = -tau_w*surfbn(ifac)/(vel(3,iel)+epzero)

!  if(prs.and.(1.d0-volume_porosi(iel))>1.d-12) then
  if(prs.and.porous_zone(iel)>0.5d0) then
    crvexp(3,iel) = tau_w*surfbn(ifac)+crvexp(3,iel)
  end if

end do

deallocate(lstelt)

if(prs) then
  call resistance_user()

  do iel=1,ncel
!!    if((1.d0-volume_porosi(iel))>1.d-12) then
    if(porous_zone(iel)>0.5d0) then
      crvexp(3,iel) = (axial_resis_f(iel)+axial_resis_d(iel)) &
                    * volume_porosi(iel)*volume(iel) &
                    + crvexp(3,iel)
    end if
  end do

end if

!!===========================
!!Periodicity (momentum)
!!===========================
if(.true.) return
!!Area of the cross section
!!surf = 6.49d-4 m2 (based on real geometry)
!!surf = 6.6059d-4 m2 (based on the coarse mesh)
!!volumetic flow rate: 60.7 m3/h
uref = 60.7d0*0.25d0/3600.d0/6.49d-4
rho_ref = den(1)
flowref = uref*rho_ref !Mass flux
flowint = flown !Mass flux of the previous time step
flown = 0.0d0
surf = 0.0d0
cptfac = 0
chaine1 = 'inlet'
allocate(lstfac(nfac))
call getfac(chaine1, nlfac, lstfac)
!! This gets the interger for a field associated with given key 
call field_get_key_int(ivarfl(iw),kimasf,iflmas)
!! Access mass flux 
call field_get_val_s(iflmas, imasfl) 

!! First loop to determine surface and current flow
do ielt = 1, nlfac
  ifac = lstfac(ielt)
  surf = surf + surfan(ifac) 
  flown = flown + abs(imasfl(ifac))
  cptfac = cptfac + 1
enddo

!! Incase of parralesim 
if (irangp.ge.0) then
  call parsom(surf)
  call parsom(flown)
  call parcpt(cptfac)
endif

flown = flown/surf

if (ntcabs .eq. 1) then
  dp = 0.0d0
  flown  = flowref
  flown1 = flown
else
  flown1 = flowint
endif

!! Check masss flow balances
dp = dp + (2.0d0*(flown-flowref)-(flown1-flowref))/(2.0d0*dt(1))

!!momentum source on all cells
do iel = 1, ncel
  if(prs) then
    value = -volume(iel)*dp*volume_porosi(iel)
  else
    value = -volume(iel)*dp
  end if
  crvexp(3,iel) = value+crvexp(3,iel)
enddo

if(irangp<=0.and.mod(ntcabs,10)==0)write(*,*) 'dp:', dp, 'Mass flux:',flown

if(mod(ntcabs,100)==0)then
  write(nfecra,*) '==FLOW RATE CONTROL=='
  write(nfecra,*) 'current time step number=',ntcabs
  write(nfecra,*) 'number of inlet faces = ', cptfac
  write(nfecra,*) 'inlet surface = ', surf
  write(nfecra,*) 'flow at time n at inlet = ', flown , 'flow forcing = ', value
  write(nfecra,*) 'pressure gradient = ',-dp
  write(nfecra,*) 'ref. flow', flowref
  write(nfecra,*) 'flow at time n-1 at inlet = ', flown1
  write(nfecra,*) 'mass flow rate = ', flowrate
end if

!--------
! Formats
!--------


!----
! End
!----

! Deallocate the temporary array
deallocate(lstfac)

return
end subroutine ustsnv

!===============================================================================

!===============================================================================
!> \brief    Additional right-hand side source terms for scalar equations (user
!>           scalars and specific physics scalars).
!>
!> \deprecated Use \ref cs_user_source_terms instead.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     iscal         index number of the current scalar
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!> \param[in]                   source terms or mass rate
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!______________________________________________________________________________!


subroutine ustssc &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field
use SubChCFD

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)

! Local variables
integer          iel

double precision, dimension(:), pointer :: imasfl

double precision, dimension(:,:), pointer :: vel, A_vel
double precision, dimension(:), pointer :: surf_rod, cvol, tfluxwall
integer, allocatable, dimension(:) :: lstelt
integer          nlelt, ielt
double precision ubulk, sumvol ,rho_ref, heat_flow
integer          iflmas, ifac, nn
character*80 :: chaine

!===============================================================================

!!===========================
!!Thermal periodicity
!!===========================
if(.true.) return

!! Impose heat sink for periodicity
if (iscal.ne.iscalt) return

call field_get_val_v(ivarfl(iu),vel)
call field_get_val_s_by_name('WallThermalFlux', tfluxwall)

allocate(lstelt(nfabor))
chaine = heated_rod_walls

call getfbr(chaine, nlelt, lstelt)
heat_flow = 0.d0
do ielt = 1, nlelt
  ifac = lstelt(ielt)
  heat_flow = heat_flow + surfbn(ifac)*tfluxwall(ifac)
end do

sumvol = 0.d0
ubulk = 1.d-15
do iel = 1, ncel
  sumvol = sumvol+volume(iel)
  ubulk = ubulk+ vel(3,iel)*volume(iel)
end do

!! Incase of parralesim sum all the totals from the different processors 
if (irangp.ge.0) then
  call parsom(heat_flow)
  call parsom(sumvol)
  call parsom(ubulk)
end if

!! Calculate heat sink
call field_get_val_s_by_name('SubchannelVolume', cvol)
call field_get_val_s_by_name('RodSurfaceArea', surf_rod)
call field_get_val_v_by_name('AverageVelocity', A_vel)
call subchannel_velocity()

do iel = 1, ncel
  crvimp(iel) = 0.d0
  crvexp(iel) = heat_flow*volume(iel)*vel(3,iel)/ubulk
end do
!--------
! Formats
!--------

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine ustssc

!===============================================================================

!===============================================================================
!> \brief    Additional right-hand side source terms for vectorial equations
!>           (user vectors and specific physics vectors).
!>
!> \deprecated Use \ref cs_user_source_terms instead.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss terms
!> \param[in]     ncesmp        number of cells with mass source terms
!> \param[in]     iscal         index number of the current scalar
!> \param[in]     icepdc        index number of cells with head loss terms
!> \param[in]     icetsm        index number of cells with mass source terms
!> \param[in]     itypsm        type of mass source term for each variable
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        head loss coefficient
!> \param[in]     smacel        value associated to each variable in the mass
!> \param[in]                   source terms or mass rate
!> \param[out]    crvexp        explicit part of the source term
!> \param[out]    crvimp        implicit part of the source term
!______________________________________________________________________________!


subroutine ustsvv &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   crvexp , crvimp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision crvexp(3,ncelet), crvimp(3,3,ncelet)

! Local variables

integer, allocatable, dimension(:) :: lstelt

!===============================================================================



!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))

!--------
! Formats
!--------

 1000 format(' User source terms for variable ',A8,/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine ustsvv
