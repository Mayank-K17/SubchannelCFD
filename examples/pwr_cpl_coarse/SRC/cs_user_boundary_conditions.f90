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
! Function:
! ---------

!> \file cs_user_boundary_conditions.f90
!>
!> \brief User subroutine which fills boundary conditions arrays
!> (\c icodcl, \c rcodcl) for solved variables.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[out]    icodcl        boundary condition code
!> \param[in]     itrifb        indirection for boundary faces ordering
!> \param[in,out] itypfb        boundary face types
!> \param[in,out] izfppp        boundary face zone number
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!_______________________________________________________________________________

subroutine cs_f_user_boundary_conditions &
 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     ,                                                       &
   rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use atchem
use atincl
use atsoil
use ctincl
use cs_fuel_incl
use mesh
use field
use turbomachinery
use iso_c_binding
use cs_c_bindings
use SubChCFD

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvar,3)

! Local variables

! INSERT_VARIABLE_DEFINITIONS_HERE

integer, allocatable, dimension(:) :: lstelt1, lstelt2
integer                            :: ifac, iel, ivar
integer                            :: ielt, nlelt1, nlelt2
integer                            :: symmtyp

double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cht

!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt1(nfabor))  ! temporary array for boundary faces selection
allocate(lstelt2(nfabor))

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

call field_get_val_v(ivarfl(iu),vel)
call field_get_val_s_by_name('CHTInternalCoupling', cht)

call getfbr(heated_rod_walls, nlelt1, lstelt1)
call getfbr(unheated_walls, nlelt2, lstelt2)

symmtyp = slip_type

do ielt = 1, nlelt1
  ifac = lstelt1(ielt)

  if(symmtyp==1) then

    itypfb(ifac) = isymet

  else if(symmtyp==2) then

    itypfb(ifac) = iparoi

    iel = ifabor(ifac)
    rcodcl(ifac, iu, 1) = vel(1,iel)
    rcodcl(ifac, iv, 1) = vel(2,iel)
    rcodcl(ifac, iw, 1) = vel(3,iel)

  else if(symmtyp==3) then

    itypfb(ifac) = iparoi

  else

!! Slip walls for velocities and homogeneous Neumann for turbulence variables
!  icodcl(ifac, iu) = 4
!  icodcl(ifac, iv) = 4
!  icodcl(ifac, iw) = 4
!  icodcl(ifac, ik) = 3
!  rcodcl(ifac,ik,1) = 0.d0
!  rcodcl(ifac,ik,2) = rinfin
!  rcodcl(ifac,ik,3) = 0.d0
!  icodcl(ifac, iep) = 3
!  rcodcl(ifac,iep,1) = 0.d0
!  rcodcl(ifac,iep,2) = rinfin
!  rcodcl(ifac,iep,3) = 0.d0

  end if

  !! Thermal boundaries on heated walls
  !! Set the heated & unheated sections
  if (nscal.gt.0) then
    icodcl(ifac, isca(iscalt)) = 3
    if(cht(ifac)>0.5d0) then
      rcodcl(ifac,isca(iscalt),3) = 0.d0
    else
      rcodcl(ifac,isca(iscalt),3) = -value_rw
      if(cdgfbo(3,ifac) < 0.03d0 .or. cdgfbo(3,ifac) > 0.63d0) then
        rcodcl(ifac,isca(iscalt),3) = 0.d0
      end if
    end if
  end if

enddo

!!--------------------------------------------------------
!! Set boundary conditions for the unheated walls
do ielt = 1, nlelt2
  ifac = lstelt2(ielt)

  if(symmtyp==1) then

    itypfb(ifac) = isymet

  else if(symmtyp==2) then

    itypfb(ifac) = iparoi

    iel = ifabor(ifac)
    rcodcl(ifac, iu, 1) = vel(1,iel)
    rcodcl(ifac, iv, 1) = vel(2,iel)
    rcodcl(ifac, iw, 1) = vel(3,iel)

  else if(symmtyp==3) then

    itypfb(ifac) = iparoi

  else

!! Slip walls for velocities and homogeneous Neumann for turbulence variables
!  icodcl(ifac, iu) = 4
!  icodcl(ifac, iv) = 4
!  icodcl(ifac, iw) = 4
!  icodcl(ifac, ik) = 3
!  rcodcl(ifac,ik,1) = 0.d0
!  rcodcl(ifac,ik,2) = rinfin
!  rcodcl(ifac,ik,3) = 0.d0
!  icodcl(ifac, iep) = 3
!  rcodcl(ifac,iep,1) = 0.d0
!  rcodcl(ifac,iep,2) = rinfin
!  rcodcl(ifac,iep,3) = 0.d0

  end if

enddo

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt1)  ! temporary array for boundary faces selection
deallocate(lstelt2)

return
end subroutine cs_f_user_boundary_conditions

!> \section examples Examples
!>   Several examples are provided
!>   \ref cs_user_boundary_conditions_examples "here".
