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
! Purpose:
! -------

!> \file cs_user_extra_operations.f90
!>
!> \brief This function is called at the end of each time step, and has a very
!>  general purpose
!>  (i.e. anything that does not have another dedicated user subroutine)
!>
!> See \ref cs_user_extra_operations_examples and
!> \ref cs_user_extra_operations-nusselt_calculation for examples.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine cs_f_user_extra_operations &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use lagran
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use field_operator
use turbomachinery
use cs_c_bindings
use SubChCFD

!===============================================================================

implicit none

interface

double precision function cs_timer_wtime()        &
  bind(C, name='cs_timer_wtime')
  use, intrinsic :: iso_c_binding
  implicit none
end function cs_timer_wtime

end interface

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

call subchannel_physical_properties()

if(prs .and. iporo_cal==3) call filtered_physical_properties()

call subchannel_Nu()

call subchannel_temperature()

call subchannel_Twall()

!! Output wall clock time/calculation time
calculation_time = cs_timer_wtime()-wall_clock_time
wall_clock_time = cs_timer_wtime()

if(irangp<=0) write(nfecra,*) 'wall clock time = ', wall_clock_time, &
                              'calculation time = ', calculation_time

!! Pressure reconstruction in a coupled simulation
if(cpl) then
  if(mod(ntcabs, fpre_rec)==0) &
    call pressure_reconstruction()
end if


!! Post-processing section
!!------------------------------------------------------------------------------
!! Calculate circumferential wall temperature of a 5x5 PWR bundle
if(.true.) call circumferential_twall ()

!! calculate the axial distribution of the circumferentially averaged 
!! wall temperature and some cross-sectional bulk quantities
if(.false.) call axial_avg_quantities ()

if(.false.) call energy_balance(nvar, nscal, dt)

!! Round Robin: extract subchannel velocities at z=235 mm
if(.false.) call round_robin_subchannel_vel()

!! In case of CHT, check Nusselt calculated based on simulation result
if(.false.) call subchannel_nusselt_check()

!! Extract the coolant channel wall temperature and thermal flux 
if(.false.) call circumferential_tflux_HTGR_GA()
!!------------------------------------------------------------------------------

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine cs_f_user_extra_operations
