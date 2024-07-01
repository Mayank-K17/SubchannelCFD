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

!> \file csccel.f90
!> \brief Exchange of coupling variables between tow instances of code_saturne
!> thanks to cells.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     ivar          variable number
!> \param[out]    crvimp        working table for implicit part
!> \param[out]    crvexp        working table for explicit part
!______________________________________________________________________________

subroutine csccel &
 ( ivar   , f_dim,                                                       &
   crvexp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use cplsat
use field
use mesh
use cs_c_bindings
use coupling_resolved

!===============================================================================

implicit none

! Arguments

integer          ivar
integer          f_dim

double precision crvexp(3,ncelet)

! Local variables

integer          f_id
integer          numcpl
integer          ncesup , nfbsup
integer          ncecpl , nfbcpl , ncencp , nfbncp
integer          ncedis , nfbdis
integer          ncecpg , ncedig
integer          ityloc , ityvar

integer, allocatable, dimension(:) :: lcecpl , lfbcpl
integer, allocatable, dimension(:) :: locpts

double precision, allocatable, dimension(:,:) :: coopts , djppts , dofpts
double precision, allocatable, dimension(:) :: pndpts
double precision, allocatable, dimension(:,:) :: rvdis, rvcel

double precision, allocatable, dimension(:,:) :: rdpdis, rdpcel
double precision, allocatable, dimension(:,:) :: coefa_dp, dp
double precision, allocatable, dimension(:,:,:) :: coefb_dp
double precision, pointer, dimension(:) :: pres, coefa_p, coefb_p
double precision, pointer, dimension(:,:) :: delta_p
integer :: iel, ii
integer :: iccocg, inc, nswrgp, imligp, iwarnp
double precision :: epsrgp, climgp
double precision, allocatable, dimension(:,:) :: dfrcxt
type(var_cal_opt) :: vcopt_p

!===============================================================================

! cannot be called for non variables ...
if (ivar.le.0) then
  f_id = -1
  call csexit(1)
else
  f_id = ivarfl(ivar)
endif

!get the dimension of the variable
call field_get_dim(f_id, f_dim)

do numcpl = 1, nbrcpl

!===============================================================================
! 1.  DEFINITION DE CHAQUE COUPLAGE
!===============================================================================

  call nbecpl                                                     &
  !==========
 ( numcpl ,                                                       &
   ncesup , nfbsup ,                                              &
   ncecpl , nfbcpl , ncencp , nfbncp )

  ! Allocate temporary arrays for coupling information
  allocate(lcecpl(ncecpl))
  allocate(lfbcpl(nfbcpl))

  call lelcpl                                                     &
  !==========
 ( numcpl ,                                                       &
   ncecpl , nfbcpl ,                                              &
   lcecpl , lfbcpl )

  deallocate(lfbcpl)

!===============================================================================
! 2.  PREPARATION DES VARIABLES A ENVOYER SUR LES CELLULES
!===============================================================================

  ityvar = 1

! --- Informations géométriques de localisation

  call npdcpl(numcpl, ncedis, nfbdis)
  !==========

  ! Allocate temporary arrays for geometric quantities
  allocate(locpts(ncedis))
  allocate(coopts(3,ncedis), djppts(3,ncedis), dofpts(3,ncedis))
  allocate(pndpts(ncedis))

  ! Allocate temporary arrays for variables exchange
  allocate(rvdis(f_dim, ncedis))
  allocate(rvcel(f_dim, ncecpl))
  allocate(rdpdis(3, ncedis))
  allocate(rdpcel(3, ncecpl))
  rdpdis = 0.d0
  rdpcel = 0.d0

  call coocpl &
  !==========
( numcpl , ncedis , ityvar , &
  ityloc , locpts , coopts , &
  djppts , dofpts , pndpts )

  if (ityloc.eq.2) then
    write(nfecra,1000)
    call csexit(1)
    !==========
  endif

  ! We check that there exists some entities to be coupled at least on one
  ! rank
  ! otherwise there is no need to call gradients for instance
  ncecpg = ncecpl
  ncedig = ncedis
  if (irangp.ge.0) then
    call parcpt(ncecpg)
    call parcpt(ncedig)
  endif

  ! --- Prepare variables to be transfer

  if (ncedig.gt.0) then

    call cscpce &
  ( ncedis , f_id   , f_dim ,                                     &
    locpts ,                                                      &
    coopts , rvdis  )

    if(send_dp.and.mod(ntcabs,fsend_dp)==0.and.f_dim==3)then

      allocate(coefa_dp(3,ncelet))
      allocate(coefb_dp(3,3,ncelet))
      allocate(dp(3,ncelet))

      call field_get_val_s(ivarfl(ipr),pres)
      call field_get_coefa_s(ivarfl(ipr), coefa_p)
      call field_get_coefb_s(ivarfl(ipr), coefb_p)

      call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)

      iccocg = 1
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

      !!Gradient calculation for normal scalar field !TODO note changes of this
      !!subroutine between v7.0 and v8.0
      call gradient_s(ivarfl(ipr), imrgra, inc, nswrgp, imligp, iwarnp, &
                      epsrgp, climgp, pres, coefa_p, coefb_p, dp)

!--------
!!check data passing 
      call field_get_val_v_by_name('delta_p', delta_p)
      do ii=1,3
        do iel=1, ncelet
          delta_p(ii,iel)=dp(ii,iel)
        end do
      end do
!--------

      coefa_dp = 0.d0
      coefb_dp = 1.0d0

      call cscpce_dp                                                   &
       ( ncedis  , ipr,                                                &
         locpts  ,                                                     &
         dp      ,                                                     &
         coefa_dp, coefb_dp ,                                         &
         coopts  , rdpdis  )

      deallocate(coefa_dp)
      deallocate(coefb_dp)
      deallocate(dp)
      deallocate(dfrcxt)

    end if

  endif

  ! Free memory
  deallocate(locpts)
  deallocate(coopts, djppts, dofpts)
  deallocate(pndpts)

  ! Cet appel est symétrique, donc on teste sur NCEDIG et NCECPG
  ! (rien a envoyer, rien a recevoir)
  if (ncedig.gt.0.or.ncecpg.gt.0) then

    ! for vectorial exchange
    call varcpl &
  ( numcpl , ncedis , ncecpl , ityvar , f_dim  ,                  &
    rvdis  ,                                                      &
    rvcel  )

    if(send_dp.and.mod(ntcabs,fsend_dp)==0.and.f_dim==3) then
      call varcpl &
    ( numcpl , ncedis , ncecpl , ityvar , f_dim ,                   &
      rdpdis ,                                                      &
      rdpcel  )
    end if

  endif

  ! Free memory
  deallocate(rvdis)
  deallocate(rdpdis)

!===============================================================================
! 3.  TRADUCTION DU COUPLAGE EN TERME DE TERMES SOURCES
!===============================================================================

  if (ncecpg.gt.0) then

    call csc2ts(ncecpl, lcecpl, f_id, f_dim, rvcel, crvexp)

  endif

  ! Free memory
  deallocate(rvcel)
  deallocate(rdpcel)
  deallocate(lcecpl)

enddo
!     Fin de la boucle sur les couplages


!--------
! FORMATS
!--------
 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :                                             ',/,&
'@    =========                                               ',/,&
'@    LE COUPLAGE VIA LES FACES EN TANT QU''ELEMENTS          ',/,&
'@    SUPPORTS N''EST PAS ENCORE GERE PAR LE NOYAU.           ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
!----
! FIN
!----

return
end subroutine
