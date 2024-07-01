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
use coupling_resolved

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

call subchannel_velocity()

call subchannel_axial_Re()

call subchannel_temperature()

!! Post-processing section
!!------------------------------------------------------------------------------
!! calculate the axial distribution of the circumferentially averaged 
!! wall temperature and some cross-sectional bulk quantities
if(.true.) call axial_avg_quantities ()

!! Calculate circumferential wall temperature, y plus and shear stress of a 5x5 
!! PWR bundle
if(.true.) call circumferential_quantities ()

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine cs_f_user_extra_operations
