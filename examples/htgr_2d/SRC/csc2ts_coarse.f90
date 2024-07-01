!-------------------------------------------------------------------------------

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

!> \file csc2ts.f90
!> \brief Code-code coupling with source terms.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     ncecpl        number of coupling
!> \param[in]     lcecpl
!> \param[in]     f_id          field index
!> \param[in]     f_dim         field dimension
!> \param[in]     rvcpce        distant variable array
!> \param[in,out] crvexp        explicit source term
!______________________________________________________________________________

subroutine csc2ts &
 ( ncecpl,  lcecpl , f_id, f_dim , rvcpce , rdpcpce, &
   crvexp,  crvimp)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use cplsat
use mesh
use field
use cs_c_bindings
use SubChCFD

!===============================================================================

implicit none

! Arguments

integer          ncecpl
integer          lcecpl(ncecpl)
integer          f_id, f_dim

double precision crvexp(f_dim, ncelet), crvimp(f_dim, f_dim, ncelet)
double precision rvcpce(f_dim, ncecpl), rdpcpce(3,ncecpl)

! Local variables

integer          numcpl
integer          isou
integer          ipt    , ielloc
double precision xdis   , xloc   , xtau   , rovtau
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer :: cvara_s
double precision, dimension(:,:), pointer :: cvara_v

integer :: i, iel

double precision, dimension(:), pointer :: volume_porosi, feedback_zone

double precision, allocatable, dimension(:,:) :: mom_imp
double precision, dimension(:,:), pointer :: mom_s, mom_cx, mom_hat, vol_f, &
                                             mom_pt
double precision, dimension(3) ::mom_c
double precision, dimension(3) :: itgtd_vol
double precision, dimension(3) :: xyzcenmin, xyzcenmax, &
                                  xyzfbmin, xyzfbmax
save xyzcenmin, xyzcenmax, xyzfbmin, xyzfbmax

double precision, pointer, dimension(:) :: pres, coefa_p, coefb_p
double precision, pointer, dimension(:,:) :: delta_p
integer :: inc, nswrgp, imligp, iwarnp
double precision :: epsrgp, climgp
double precision, allocatable, dimension(:,:) :: dfrcxt
double precision, allocatable, dimension(:,:) :: dp
type(var_cal_opt) :: vcopt_p

double precision, dimension(:), pointer :: mom_t, received_temp
double precision, allocatable, dimension(:) :: xcp
double precision, dimension(:), pointer :: cpro_cp
double precision :: mom_t_total, volume_total

!----------------------------------------------------------------------------------
if(.not.cpl) return

! For scalars
if (f_dim.eq.1) then
  if(.not.feedback_energy) return

  call field_get_val_prev_s(f_id, cvara_s)
  call field_get_val_s(icrom, crom)
  call field_get_val_s_by_name('EnergySourceCoupling', mom_t)
  call field_get_val_s_by_name('received_T', received_temp)
  if (icp.ge.0) call field_get_val_s(icp, cpro_cp)

  allocate(xcp(ncel))

! If it is a temperature
  if (itherm.eq.1) then
    if (icp.ge.0) then
      do iel = 1, ncel
        xcp(iel) = cpro_cp(iel)
      enddo
    else
      do iel = 1, ncel
        xcp(iel) = cp0
      enddo
    endif
  else
    do iel = 1, ncel
      xcp(iel) = 1.d0
    enddo
  endif

  if(iphydr==1) &
    call field_get_val_v_by_name('volume_forces', vol_f)

  xtau = 1.d0/lambda_energy*dtref
  !! activate energy source
  do ipt = 1, ncecpl

    ielloc = lcecpl(ipt)

    xdis = rvcpce(1,ipt)
    xloc = cvara_s(ielloc)
    received_temp(ielloc) = xdis

    mom_t(ielloc) = (xdis-xloc)*crom(ielloc)*xcp(ielloc)/xtau
    crvexp(1,ielloc) = crvexp(1,ielloc)+mom_t(ielloc)*cell_f_vol(ielloc)

  enddo

  mom_t_total = 0.d0
  volume_total = 0.d0
  do ipt = 1, ncecpl
    ielloc = lcecpl(ipt)
    mom_t_total = mom_t_total+mom_t(ielloc)*cell_f_vol(ielloc)
    volume_total = volume_total+cell_f_vol(ielloc)
  enddo

  if (irangp.ge.0) then
    call parsom(mom_t_total)
    call parsom(volume_total)
  endif

  if (irangp.le.0) write(*,*)mom_t_total

  mom_t_total = mom_t_total/volume_total

  do ipt = 1, ncecpl

    ielloc = lcecpl(ipt)
    crvexp(1,ielloc) = crvexp(1,ielloc)-mom_t_total*cell_f_vol(ielloc)

  enddo

  deallocate(xcp)

! For vectors and tensors
else

  call field_get_val_prev_v(f_id, cvara_v)
  call field_get_val_s(icrom, crom)
  call field_get_val_v_by_name('MomentumSource', mom_s)
  call field_get_val_v_by_name('MomSourceCX', mom_cx)
  call field_get_val_v_by_name('MomSourceSHat', mom_hat)
  call field_get_val_v_by_name('MomSourcePotential', mom_pt)
  call field_get_val_s_by_name('FeedbackFlag', feedback_zone)

  if(iphydr==1) &
    call field_get_val_v_by_name('volume_forces', vol_f)
  if(prs) &
    call field_get_val_s_by_name('porosity', volume_porosi)

  feedback_zone = 1.d0

  if(feedback_user) call feedback_zone_user()

  allocate (mom_imp(3,ncelet))

  xtau = 1.d0/lambda*dtref

  do ipt = 1, ncecpl

    ielloc = lcecpl(ipt)

    do isou = 1, 3
      xdis = rvcpce(isou,ipt)
      xloc = cvara_v(isou,ielloc)

      if(feedback .and. int(feedback_zone(ielloc))==1) then
        mom_s(isou,ielloc) = (xdis-xloc)*crom(ielloc)/xtau
        mom_imp(isou,ielloc) = crom(ielloc)/xtau
!!      Maybe better to move to the final stage for source term calculation
        if(prs)then
          mom_s(isou,ielloc) = mom_s(isou,ielloc)*volume_porosi(ielloc)
          mom_imp(isou,ielloc) = mom_imp(isou,ielloc)*volume_porosi(ielloc)
        end if
      else
        mom_s(isou,ielloc) = 0.d0
        mom_imp(isou,ielloc) = 0.d0
      end if

    enddo

  enddo

  select case (mom_crr)
  case (0)
    do i=1, 3
      do iel=1, ncel
        mom_hat(i,iel) = mom_s(i,iel)
      end do
    end do
  !! The following section may need revision when porosity is considered
  case (1)
    mom_hat = 0.d0
    mom_cx = 0.d0
    mom_c = 0.d0
    itgtd_vol = 0.d0
    do i=1, 3
      do iel=1, ncel
        mom_c(i) = mom_c(i)+mom_s(i,iel)*volume(iel)
        itgtd_vol(i) = itgtd_vol(i)+volume(iel)
      end do

      if (irangp.ge.0) then
        call parsom(mom_c(i))
        call parsom(itgtd_vol(i))
      endif

      mom_c(i) = mom_c(i)/itgtd_vol(i)

      do iel=1, ncel
        mom_cx(i,iel) = mom_c(i)
        mom_hat(i,iel) = mom_s(i,iel)-mom_cx(i,iel)
      end do

    end do

  case (2)
    mom_hat = 0.d0
    do numcpl=1,nbrcpl

    !! Calculate the coordinate ranges of the coupling region
    if(ntcabs==ntpabs+1) call coupling_xyz(numcpl, &
                              xyzcenmin, xyzcenmax, xyzfbmin, xyzfbmax)

    mom_cx = 0.d0
    mom_c = 0.d0
    itgtd_vol = 0.d0
    do i=1, 3
      do iel=1, ncel
        if(xyzcen(i,iel)>xyzfbmin(i) .and. xyzcen(i,iel)<xyzfbmax(i)) then
          mom_c(i) = mom_c(i)+mom_s(i,iel)*volume(iel)
          itgtd_vol(i) = itgtd_vol(i)+volume(iel)
        end if
      end do

      if (irangp.ge.0) then
        call parsom(mom_c(i))
        call parsom(itgtd_vol(i))
      endif

      mom_c(i) = mom_c(i)/itgtd_vol(i)
      do iel=1, ncel
        if(xyzcen(i,iel)>xyzfbmin(i) .and. xyzcen(i,iel)<xyzfbmax(i)) then
          mom_cx(i,iel) = mom_c(i)
          mom_hat(i,iel) = mom_hat(i,iel)+mom_s(i,iel)-mom_cx(i,iel)
        end if
      end do

    end do
    
    end do
      
  case (3)
  !!To be added

  end select

  if(ipre_rec==2.and.mod(ntcabs,fpre_rec)==0)then

    allocate(dp(3,ncelet))
    call field_get_val_s(ivarfl(ipr),pres)
    call field_get_coefa_s(ivarfl(ipr), coefa_p)
    call field_get_coefb_s(ivarfl(ipr), coefb_p)  

    call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)

    inc = 0

    call grdpor(inc)

    if (iphydr.eq.1.or.iifren.eq.1) inc = 1
    nswrgp = vcopt_p%nswrgr
    imligp = vcopt_p%imligr
    iwarnp = vcopt_p%iwarni
    epsrgp = vcopt_p%epsrgr
    climgp = vcopt_p%climgr

    allocate(dfrcxt(3,ncelet))
    dfrcxt = 0.d0

    !Gradient calculation for normal scalar field !TODO
    call gradient_s(ivarfl(ipr), imrgra, inc, nswrgp, imligp, &
                    iwarnp, epsrgp, climgp, pres, coefa_p, coefb_p, dp)

    !Gradient calculation for potential type scalar field
    !no longer available in version 8.0
    !call gradient_potential_s &
    ! (ivarfl(ipr)     , imrgra , inc    , iccocg , nswrgp , imligp , &
    !  iphydr , iwarnp , epsrgp , climgp , extrap ,                   &
    !  dfrcxt ,                                                       &
    !  pres   , coefa_p        , coefb_p        ,                     &
    !  dp  ) 

    do ipt = 1, ncecpl
      ielloc = lcecpl(ipt)
      if(feedback .and. int(feedback_zone(ielloc))==1) then
        do isou = 1, 3
          xdis = rdpcpce(isou,ipt)
          xloc = dp(isou,ielloc)
          mom_pt(isou,ielloc) = xloc-xdis
        enddo
      else
        mom_pt(isou,ielloc) = 0.d0
      end if
    enddo
      
    deallocate(dfrcxt, dp)

  end if

  do ipt = 1, ncecpl

    ielloc = lcecpl(ipt)
    
    do isou = 1, 3
      xdis = rvcpce(isou,ipt)
      xloc = cvara_v(isou,ielloc)
      crvexp(isou,ielloc) = crvexp(isou,ielloc) &
                          + mom_hat(isou,ielloc)*cell_f_vol(ielloc)
      if(imp_source)then
        crvexp(isou,ielloc) = crvexp(isou,ielloc) &
                            + mom_imp(isou,ielloc)*xloc*cell_f_vol(ielloc)
        crvimp(isou,isou,ielloc) = -mom_imp(isou,ielloc)*cell_f_vol(ielloc)
      end if

    enddo

  enddo

  deallocate (mom_imp)

  !--------
  !!check data passing 
  call field_get_val_v_by_name('delta_p', delta_p)

  do ipt = 1, ncecpl
    ielloc = lcecpl(ipt)
    do isou = 1, 3
      delta_p(isou,ielloc) = rdpcpce(isou,ipt)
    enddo
  enddo
  !--------

endif

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
