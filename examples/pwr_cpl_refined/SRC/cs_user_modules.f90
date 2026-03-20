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

!> \file cs_user_modules.f90
!>
!> \brief User-defined module: it allows to create any user array.
!>
!> See \ref cs_user_modules for examples.
!>
!> This file is compiled before all other user Fortran files.
!> To ensure this, it must not be renamed.
!>
!> The user may define an arbitrary number of modules here, even though
!> only one is defined in the example.
!
!> \cond DOXYGEN_SHOULD_SKIP_THIS

!-------------------------------------------------------------------------------

module coupling_resolved

implicit none

!=============================================================================== 
!! Parameters specification
!!------------------------------------------------------------------------------
!! Bundle type
! PWR: 1, AGR: 2
  integer,parameter :: bundle_type = 1
!! Subchannel mesh name
! PWR: 'PWR_subchannel_mesh.dat', AGR: 'AGR_subchannel_mesh.dat'
  character(len=30),parameter :: subchannel_mesh_name = 'PWR_subchannel_mesh.dat'
!! Boundary name
!        Pin wall |      Housing wall        | Guide tube wall | Inlet | Outlet
! PWR:  rod_walls | side_walls/housing_walls |    N/A          | inlet | outlet
! AGR:   RodWall  |       SleeveWall         |  GuideTubeWall  | Inlet | Outlet
  character(len=30),parameter :: pin =     'rod_walls',        &
                                 housing = 'housing_walls',       &
                                 guide   = 'N/A',              &
                                 inlet   = 'inlet',            &
                                 outlet  = 'outlet'

!! Geometrical information
!!------------------------------------------------------------------------------
!! Specific for AGR bundles
! voidage of the fuel pins
  double precision, parameter :: voidage = 0.729558
!! Subchannel parameters
! bundle sleeve radius (m)
  double precision,parameter :: bundle_R = 0.095d0
! coordinates of origin for the circular bundle
  double precision,parameter :: x_origin = 0.d0, &
                                y_origin = 0.d0, &
                                z_origin = 0.d0
! maximum number of nodes in subchannel mesh elements
  integer,parameter :: nCol_max = 8
  
!!------------------------------------------------------------------------------
!! Specific for PWR bundles
! Parameters for identification of the subchannel type
  double precision, parameter :: spacing_edge = 0.95d0
  double precision, parameter :: spacing_corner = 0.45d0
! Rod diameter
  double precision, parameter :: rod_diameter = 0.0095d0 !(m)
! Pitch-to-diameter ratio for interior subchannels
  double precision, parameter :: PbyD_interior = 1.326d0
! Pitch-to-diameter ratio for the edge subchannels
  double precision, parameter :: PbyD_edge = 1.326d0
! Pitch-to-diameter ratio for the corner subchannels
  double precision, parameter :: PbyD_corner = 1.326d0
! Heating or cooling
  logical, parameter :: heating = .true.

!!------------------------------------------------------------------------------
!! Mesh parameters  
! Scale the meshes if generated in other units than 'meters'
  logical :: subchannel_mesh_scale = .true.
  logical :: mesh_scale = .true.
  double precision, parameter :: subchannel_scaling_factor = 1.d-3
  double precision, parameter :: mesh_scaling_factor = 1.d-3
! Use user defined routines for non-ideal rod bundles
  logical :: user_subchannel_mesh = .false.
  logical :: user_frictional_factor = .false.
  logical :: user_nusselt_number = .false.

!!------------------------------------------------------------------------------
!! Physical parameters
! Heat flux on the heated rod walls (W/m2)
  double precision, parameter :: heat_flux = 6.76923d5
  double precision, parameter :: crr_Nu = 1.d0 !0.82d0
!! Use data table to define spatially variable physical properties
  logical :: physical_property_table = .false.
  character(len=50) :: name_data_table = 'water_properties_156bar.txt'
  !! Increment of the thermal variable in the property table (should be consistent
!! with the data table created. Data table should follow the formats shown in 
!! the examples)
  double precision,parameter :: increment = 0.1d0

!!------------------------------------------------------------------------------
!! Other parameters
!! Math constants  
  double precision, parameter :: ee = 2.71828183, pi = 3.14159265

!!------------------------------------------------------------------------------
!! Coupling module
!! Sending pressure gradient to coarse mesh baseline model
  logical :: send_dp
!! Frequency of sending pressure gradient to the coarse mesh baseline model   
  integer :: fsend_dp
!! Definition of boundary conditions of the coupling interfaces (useful for 
!! velocity, scalars and turbulence quantities, pressure is not included)
!! 0: default
!! 1: user defined
  integer :: interface_user = 0

!!------------------------------------------------------------------------------
!! Simulation statistics  
!! wall clock time
  double precision :: wall_clock_time = 0.
  double precision :: calculation_time = 0.

!===============================================================================  
!! Subchannel variables  
  integer :: n_cell
  integer, allocatable, dimension(:) :: cell_type
  double precision, allocatable, dimension(:) :: cell_height, cell_vol
  double precision, allocatable, dimension(:) :: x_channel, &
                                                 y_channel, &
                                                 z_channel
  double precision, allocatable, dimension(:) :: subchannel_surf,     &
                                                 subchannel_vol,      &
                                                 subchannel_rodsurf,  &
                                                 subchannel_wallsurf, &
                                                 PbyD
  double precision, allocatable, dimension(:,:) :: xx, yy, zz
  character(len=80) :: all_rod_walls,      &
                       heated_rod_walls,   &
                       unheated_walls,     &
                       all_walls 

  double precision, allocatable, dimension (:,:) :: user_prop
  integer :: n_rows

!===============================================================================

  contains
!===============================================================================
!! This routine is used to pass mesh scaling factor to cs_user_mesh.c
  function set_scaling_factor() &
    bind (C,name='set_scaling_factor')

    use, intrinsic :: iso_c_binding

    implicit none

    real(kind=c_double) :: set_scaling_factor

    if (mesh_scale) then
      set_scaling_factor = mesh_scaling_factor
    else
      set_scaling_factor = 1.d0
    end if

  end function set_scaling_factor

!===============================================================================
!! This routine is used to pass integer variable interface_user to .c routines
  function set_interface_user() bind (C,name='set_interface_user')
  
    use, intrinsic :: iso_c_binding
  
    implicit none
 
    integer(kind=c_int) :: set_interface_user
 
    set_interface_user = interface_user

  end function set_interface_user
  
!===============================================================================
!! This routine should be called after subchannel_mesh to allocate memories 
!! mainly for subchannel variables 
  subroutine allocate_var()

    implicit none

!! Arrays for subchannel wall surface area and volume (excluding solid)
    allocate(subchannel_surf(n_cell), subchannel_vol(n_cell), &
             subchannel_rodsurf(n_cell), subchannel_wallsurf(n_cell))
  
  end subroutine allocate_var


!===============================================================================
!! This routine should be called to dump the memory of common variables
!! after the whole of simulation is completed
  subroutine deallocate_var()

    implicit none

    deallocate(cell_height, cell_vol, cell_type)
    deallocate(subchannel_surf, subchannel_vol, subchannel_rodsurf, &
               subchannel_wallsurf, PbyD)
  
  end subroutine deallocate_var

!===============================================================================
!! This routine is used to bridge the subchannel mesh and computing mesh
  subroutine subchannel_mesh()

    implicit none

    select case (bundle_type)
    case (1)
      all_rod_walls = trim(pin)
      heated_rod_walls = trim(pin)
      unheated_walls = trim(housing)
      all_walls = trim(pin)//' or '//trim(housing)
      if (.not. user_subchannel_mesh) then
        call subchannel_mesh_PWR()
      else
        call subchannel_mesh_user()
      end if
    case (2)
      all_rod_walls = trim(pin)//' or '//trim(guide)
      heated_rod_walls = trim(pin)
      unheated_walls = trim(housing)//' or '//trim(guide)
      all_walls = trim(pin)//' or '//trim(guide)//' or '//trim(housing)
      if (.not. user_subchannel_mesh) then
        call subchannel_mesh_AGR()
      else
        call subchannel_mesh_user()
      end if
    end select

  end subroutine subchannel_mesh


!===============================================================================
!! This routine is to acquire the data dimensions of the subchannel mesh 
  subroutine subchannel_mesh_dimensions(mesh_file, n_node, n_entity)
    
    implicit none

!!  Aguments
    character(len=30), intent(IN) :: mesh_file
    integer, intent(OUT) :: n_node, n_entity

!!  Local variables
    character(1) :: dummy
    logical :: alive

!!  Check if mesh file exists
    inquire(file=mesh_file, exist=alive)

    if(alive) then
      open(100, file=mesh_file, status='old') 
    else
      write(*,*)"subchannel mesh doesn't exist. Program terminated"
      stop
    end if

!!  Read in number of nodes and 1D, 2D and 3D elements 
    read(100,*)n_node, n_entity

    close(100)

  end subroutine subchannel_mesh_dimensions


!===============================================================================  
!! This routine is to read in the data of the subchannel mesh
  subroutine subchannel_mesh_readin_pwr(mesh_file, n_node, n_entity, &
                                        x, y, z, index_cell, n_cell)
    
    implicit none

!!    Aguments
    character(len=30), intent(IN) :: mesh_file
    integer, intent(IN) :: n_node, n_entity
    double precision, dimension(:), intent (OUT) ::  x, y, z
    integer, dimension(:,:), intent (OUT) :: index_cell
    integer, intent (OUT) :: n_cell 

!!    Local variables
    character(1) :: dummy
    logical :: alive
    integer :: itemp, stat
    integer :: i

 !!  Check if mesh file exists
    inquire(file=mesh_file, exist=alive)

    if(alive) then
      open(100, file=mesh_file, status='old')
    else
      write(*,*)"subchannel mesh doesn't exist. Program terminated"
      stop
    end if

    read(100,*)

 !! Read in node cordinates of the subchannel mesh
    do i=1,n_node
      read(100,*)dummy, x(i), y(i), z(i)
    end do

 !! read in cell information
    n_cell = 0
    do
      read(100,*, iostat=stat)dummy, itemp
      if(stat/=0)exit
      if(mod(itemp,100)==8)then
        backspace(100)
        n_cell = n_cell+1
        read(100,*, iostat=stat)dummy, itemp, index_cell(n_cell, 1:8)
      end if
    end do

    close(100)

  end subroutine subchannel_mesh_readin_pwr


!===============================================================================
!! This routine is for standard PWR bundles
  subroutine subchannel_mesh_PWR()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: n_node, n_entity
    double precision, allocatable, dimension(:) ::  x, y, z
    integer, allocatable, dimension(:,:) :: index_cell
    integer :: i, iel
    double precision :: SubchannelID_err
    double precision, dimension(:), pointer :: subchannel_ID, &
                                               subchannel_type, &
                                               pitch2diameter

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)

    call subchannel_mesh_dimensions(subchannel_mesh_name, n_node, n_entity)

    allocate(x(n_node),y(n_node),z(n_node))
    allocate(index_cell(n_entity, 8))

    call subchannel_mesh_readin_pwr(subchannel_mesh_name, n_node, n_entity, &
                                    x, y, z, index_cell, n_cell)

!! Scale mesh from mm to m
    if(subchannel_mesh_scale)then
      x = x*subchannel_scaling_factor
      y = y*subchannel_scaling_factor
      z = z*subchannel_scaling_factor
    end if
           
!! Module variable, should be deallocated after the whole simulation
    allocate(xx(n_cell, 2), yy(n_cell, 2), zz(n_cell, 2))
    allocate(cell_height(n_cell))
    allocate(cell_vol(n_cell))
    allocate(cell_type(n_cell))
    allocate(PbyD(n_cell))
    allocate(x_channel(n_cell),y_channel(n_cell),z_channel(n_cell)) 
    do i=1, n_cell
      xx(i,1) = minval(x(index_cell(i,:)))
      xx(i,2) = maxval(x(index_cell(i,:)))
      yy(i,1) = minval(y(index_cell(i,:)))
      yy(i,2) = maxval(y(index_cell(i,:)))
      zz(i,1) = minval(z(index_cell(i,:)))
      zz(i,2) = maxval(z(index_cell(i,:)))
!! Calculate cell height and cell volume of the subchannel mesh
      cell_height(i) = abs(zz(i,2)-zz(i,1))
      cell_vol(i) = (xx(i,2)-xx(i,1))*(yy(i,2)-yy(i,1))*(zz(i,2)-zz(i,1))
!! Calculate cell center coordinates of the subchannel cell
      x_channel(i) = 0.5d0*(xx(i,2)+xx(i,1))
      y_channel(i) = 0.5d0*(yy(i,2)+yy(i,1))
      z_channel(i) = 0.5d0*(zz(i,2)+zz(i,1))
    end do

!! cell_type: 1 interior channel; 2 edge channel; 3 corner channel
    cell_type = 1
    where(cell_vol/cell_height<spacing_edge*maxval(cell_vol/cell_height)) &
      cell_type = 2
    where(cell_vol/cell_height<spacing_corner*maxval(cell_vol/cell_height)) &
      cell_type = 3

!! Set pitch-to-diameter ratio for the subchannels
    where(cell_type==1)PbyD = PbyD_interior
    where(cell_type==2)PbyD = PbyD_edge
    where(cell_type==3)PbyD = PbyD_corner

!! Find the subchannel index for each control volume
    do iel=1, ncel
      do i=1, n_cell
        if(xyzcen(1,iel)>=xx(i,1).and.xyzcen(1,iel)<=xx(i,2).and.&
           xyzcen(2,iel)>=yy(i,1).and.xyzcen(2,iel)<=yy(i,2).and.&
           xyzcen(3,iel)>=zz(i,1).and.xyzcen(3,iel)<=zz(i,2))then
          subchannel_ID(iel) = i
        end if
      end do
    end do

!! Check the subchannel index for each control volume
    SubchannelID_err = 0.d0
    do iel=1, ncel
      if(abs(subchannel_ID(iel))<1.d-12) SubchannelID_err = 1.d0
    end do

    if(irangp.ge.0)call parsom(SubchannelID_err)

    if(SubchannelID_err>0.99d0)then
      if(irangp<=0)then 
        write(*,*) 'Error: some cells have wrong subchannel IDs'
        write(*,*) 'Please check the subchannel mesh!'
      end if
      stop
    end if

    call field_get_val_s_by_name('PitchtoDiameterRatio', pitch2diameter)
    call field_get_val_s_by_name('SubchannelType', subchannel_type)
    
    do iel=1, ncel
      pitch2diameter(iel) = PbyD(int(subchannel_ID(iel)))
      subchannel_type(iel) = cell_type(int(subchannel_ID(iel)))
    end do

    deallocate(x,y,z)
    deallocate(index_cell)
    
  end subroutine subchannel_mesh_PWR


!===============================================================================  
!! This routine is to read in the data of the subchannel mesh
  subroutine subchannel_mesh_readin_agr(mesh_file, n_node, n_entity, &
                                        x, y, z, index_cell, n_cell)
    
    implicit none

!!    Aguments
    character(len=30), intent(IN) :: mesh_file
    integer, intent(IN) :: n_node, n_entity
    double precision, dimension(:), intent (OUT) ::  x, y, z
    integer, dimension(:,:), intent (OUT) :: index_cell
    integer, intent (OUT) :: n_cell 

!!    Local variables
    character(1) :: dummy
    Character(len=512) :: cLine
    logical :: alive
    integer :: itemp, stat
    integer :: i
    integer :: nCol

 !!  Check if mesh file exists
    inquire(file=mesh_file, exist=alive)

    if(alive) then
      open(100, file=mesh_file, status='old')
    else
      write(*,*)"subchannel mesh doesn't exist. Program terminated"
      stop
    end if

    read(100,*)

 !! Read in node cordinates of the subchannel mesh
    do i=1,n_node
      read(100,*)dummy, x(i), y(i), z(i)
    end do

 !! read in cell information
    n_cell = 0
    do
      read(100,*, iostat=stat)dummy, itemp
      if(stat/=0)exit
      if(itemp/100==3)then
        backspace(100)
        n_cell = n_cell+1
        read(100,'(a512)')cLine
        nCol = nColumns_irr(cLine)-2
        backspace(100)
        read(100,*, iostat=stat)dummy, itemp, index_cell(n_cell, 1:nCol)
      end if
    end do

    close(100)

  end subroutine subchannel_mesh_readin_agr


!===============================================================================
!! This routine is for standard AGR bundles
  subroutine subchannel_mesh_AGR()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: n_node, n_entity
    double precision, allocatable, dimension(:) ::  x, y, z
    double precision, allocatable, dimension(:,:) :: xx, yy, zz
    double precision, allocatable, dimension(:,:) :: x_face1, y_face1, &
                                                     z_face1, x_face2, &
                                                     y_face2, z_face2
    integer, allocatable, dimension(:,:) :: index_cell
    integer :: i, j, k1, k2, iel
    double precision, dimension(:), pointer :: subchannel_ID, &
                                               subchannel_type, &
                                               pitch2diameter
    logical :: within

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)

    call subchannel_mesh_dimensions(subchannel_mesh_name, n_node, n_entity)
    
    allocate(x(n_node),y(n_node),z(n_node))
    allocate(index_cell(n_entity, nCol_max))
    x = 0.
    y = 0.
    z = 0.
    index_cell = -1

    call subchannel_mesh_readin_agr(subchannel_mesh_name, n_node, n_entity, &
                                    x, y, z, index_cell, n_cell)

!! Scale mesh from mm to m
    if(subchannel_mesh_scale)then
      x = x*subchannel_scaling_factor
      y = y*subchannel_scaling_factor
      z = z*subchannel_scaling_factor
    end if

!! Check value for debug
!    do i=1,n_cell
!      write(*,*)index_cell(i,1:nCol_max)
!    end do

!! Calculate minimum and maximum coordinates value of the subchannel cell     
    allocate(xx(n_cell, nCol_max), &
             yy(n_cell, nCol_max), &
             zz(n_cell, nCol_max))
    allocate(x_face1(nCol_max/2,n_cell), &
             y_face1(nCol_max/2,n_cell), &
             z_face1(nCol_max/2,n_cell), &
             x_face2(nCol_max/2,n_cell), &
             y_face2(nCol_max/2,n_cell), &
             z_face2(nCol_max/2,n_cell))
    xx = 0.d0
    yy = 0.d0
    zz = 0.d0
    x_face1 = 0.d0
    x_face2 = 0.d0
    y_face1 = 0.d0
    y_face2 = 0.d0
    z_face1 = 0.d0
    z_face2 = 0.d0

!! Module variable, should be deallocated after the whole simulation
    allocate(cell_height(n_cell))
    allocate(cell_vol(n_cell))
    allocate(cell_type(n_cell))
    allocate(PbyD(n_cell))
    allocate(x_channel(n_cell),y_channel(n_cell),z_channel(n_cell))
    PbyD = 1.5

    do i=1, n_cell
      xx(i,:) = x(index_cell(i,:))
      yy(i,:) = y(index_cell(i,:))
      zz(i,:) = z(index_cell(i,:))

!! Calculate cell height of the subchannel mesh
      k1 = 1
      k2 = 1
      x_face1(k1,i) = xx(i,1)
      y_face1(k1,i) = yy(i,1)
      z_face1(k1,i) = zz(i,1)
      do j=2,6
        if(abs(zz(i,1)-zz(i,j))>1.e-10)then
          cell_height(i) = abs(zz(i,1)-zz(i,j))
          x_face2(k2,i) = xx(i,j)
          y_face2(k2,i) = yy(i,j)
          z_face2(k2,i) = zz(i,j)
          k2 = k2+1
        else
          k1 = k1+1
          x_face1(k1,i) = xx(i,j)
          y_face1(k1,i) = yy(i,j)
          z_face1(k1,i) = zz(i,j)
        end if
      end do

    end do

!! Calculate the cell volume of the subchannel mesh
  do i=1, n_cell
    call SubchannelVol(cell_vol(i))
  end do

!! Identify the cell type of the subchannel mesh
  do i=1, n_cell
    call SubchannelType(cell_type(i))
  end do

!! Find the subchannel indext for each control volume 
    do iel=1, ncel
      within = .false.
      do i=1, n_cell
        call find_suchannel_index(xyzcen(1,iel), xyzcen(2,iel), x_face1(:,i),&
                                  y_face1(:,i), index_cell(i,:), within)
        if((xyzcen(3,iel)-z_face1(1,i))*(xyzcen(3,iel)-z_face2(1,i))<0.&
           .and.within)&
          subchannel_ID(iel) = i
      end do
    end do

    call field_get_val_s_by_name('PitchtoDiameterRatio', pitch2diameter)
    call field_get_val_s_by_name('SubchannelType', subchannel_type)

    do iel=1, ncel
      pitch2diameter(iel) = PbyD(int(subchannel_ID(iel)))
      subchannel_type(iel) = cell_type(int(subchannel_ID(iel)))
    end do

    deallocate(x,y,z)
    deallocate(index_cell)
    deallocate(xx,yy,zz)
    deallocate(x_face1,y_face1,z_face1,x_face2,y_face2,z_face2)

  end subroutine subchannel_mesh_AGR


!===============================================================================  
!! This routine is used to acquire the number of columns in a irregular data 
!! file
  Integer function nColumns_irr(cStr)
    Character(len = *) , Intent(IN) :: cStr
    Integer :: i
    Logical :: bIsSeparator , bIsQuote
    nColumns_irr = 0
    bIsSeparator = .true.
    bIsQuote = .false.
    do i = 1,len_trim(cStr)
      select case(cStr(i:i))
      case('"' , "'")
        if (.not.bIsQuote) nColumns_irr = nColumns_irr + 1
        bIsQuote = .not.bIsQuote
        bIsSeparator = .false.
      case(" " , "," , char(9))
        if (.not.bIsQuote) then
          bIsSeparator = .true.
        end if
      case default
        if (bIsSeparator) then
          nColumns_irr = nColumns_irr + 1
        end if
        bIsSeparator = .false.
      end select
    end do

  end function nColumns_irr


!===============================================================================
!! This routine is used to acquire the number of columns in a regular data file 
  function nColumn_r(vars) result(numvars)

    character(len=*), intent(in) :: vars
    integer :: numvars
    character(len=len(vars)) :: tmpvars
    character(len=256) :: tmpvar
    tmpvars = trim(adjustl(vars))
    numvars = 0
    do while (len_trim(tmpvars)>0)
      read(tmpvars, *) tmpvar
      numvars = numvars + 1
      tmpvars = tmpvars(index(tmpvars, trim(tmpvar))+len_trim(tmpvar):)
    end do

  end function nColumn_r


!===============================================================================
!! This routine is used to calculate the volume of the subchannel mesh elements 
!! Need to be added (currently not necessary in the code) 
  subroutine SubchannelVol(cellVol)
    
    implicit none

    double precision, intent(OUT) :: cellVol

  end subroutine SubchannelVol


!===============================================================================
!! This routine is used to identify the type of the subchannel mesh elements 
!! Need to be added (This is important to calculate Cf)  
  subroutine SubchannelType(cellType)

    implicit none

    integer, intent(OUT) :: cellType

    cellType = 1

  end subroutine SubchannelType 


!===============================================================================
!! This routine is to judge if (x0,y0) is within an axial subchannel face cell
  subroutine find_suchannel_index(x0, y0, x, y, index_cell, within)

    implicit none

    double precision, dimension(:), intent(IN) :: x, y
    double precision, intent(IN) :: x0, y0
    integer, dimension(:), intent(IN) :: index_cell
    logical, intent(OUT) :: within

    double precision, dimension(:), allocatable :: same_side
    integer :: nNodes
    integer :: i,j,k

    within = .true.

    nNodes = 0
    do i=1, size(index_cell)
      if(index_cell(i)/=-1)nNodes = nNodes+1
    end do 
    nNodes = nNodes/2

    allocate(same_side(nNodes))

    do i=1, nNodes
      j = i
      k = i 
      if(i>nNodes-2)then
        k = i-nNodes
        j = i
      end if
      if(i>nNodes-1)then
        j = i-nNodes
      end if
      same_side(i) = (x0-x(j+1))*(y(i)-y(j+1))-(y0-y(j+1))*(x(i)-x(j+1))
      same_side(i) = same_side(i)*((x(k+2)-x(j+1)) &
                   * (y(i)-y(j+1))-(y(k+2)-y(j+1))*(x(i)-x(j+1)))
      if(edge_corr(x(i),y(i),x(j+1),y(j+1)))same_side(i)=1.d0
    end do

    do i=1,nNodes
      if(same_side(i)<0)within = .false.
    end do

  end subroutine find_suchannel_index


!===============================================================================
! This function is used for edge channel correction
  function edge_corr(x1, y1, x2, y2)

    implicit none

    double precision, intent(IN) :: x1, y1, x2, y2
    logical :: edge_corr

    double precision :: radius1, radius2

    edge_corr = .false.
    
    radius1 = sqrt((x1-x_origin)**2+(y1-y_origin)**2)
    radius2 = sqrt((x2-x_origin)**2+(y2-y_origin)**2)

    if(abs(radius1-bundle_R)<1.e-10 .and. abs(radius2-bundle_R)<1.e-10)&
      edge_corr = .true.

  end function edge_corr 

!===============================================================================
!! This routine is used to calculate the hydraulic diameter of each subchannel
!! Not accurate on very coarse computing mesh  
  subroutine hydraulic_diameter()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    character(len=80) :: chaine1, chaine2
    integer :: i, iel, ielt, ifac, nlelt1, nlelt2
    integer, allocatable, dimension(:) :: lstelt1, lstelt2
    double precision, dimension(:), pointer :: Dh, subchannel_ID
    double precision, dimension(:), pointer :: surf, cvol, surf_rod, surf_wall
    integer, allocatable, dimension(:) :: nf, nc
    
!! Check number of wall-adjacent cells in each suchannel cell  
    allocate(nf(n_cell), nc(n_cell))
    nf = 0
    nc = 0

!! Get all wall faces
    allocate(lstelt1(nfabor))
    allocate(lstelt2(nfabor))
    chaine1 = all_rod_walls
    chaine2 = housing
    call getfbr(chaine1, nlelt1, lstelt1)
    call getfbr(chaine2, nlelt2, lstelt2)

!! Field arrays used
    call field_get_val_s_by_name('HydraulicDiameter', Dh)
    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('TotalSurfaceArea', surf)
    call field_get_val_s_by_name('SubchannelVolume', cvol)
    call field_get_val_s_by_name('RodSurfaceArea', surf_rod)
    call field_get_val_s_by_name('HousingWallSurfaceArea', surf_wall)

!! Calculate the wall surface areas of each subchannel
    subchannel_surf = 0.d0
    subchannel_rodsurf = 0.d0
    subchannel_wallsurf = 0.d0
    do ielt=1, nlelt1
      ifac = lstelt1(ielt)
      iel = ifabor(ifac)
      do i=1, n_cell
        if(int(subchannel_ID(iel))==i)then
          subchannel_surf(i) = subchannel_surf(i)+surfbn(ifac)
          subchannel_rodsurf(i) = subchannel_rodsurf(i)+surfbn(ifac)
          nf(i) = nf(i)+1
        end if
      end do
    end do
    do ielt=1, nlelt2
      ifac = lstelt2(ielt)
      iel = ifabor(ifac)
      do i=1, n_cell
        if(int(subchannel_ID(iel))==i)then
          subchannel_surf(i) = subchannel_surf(i)+surfbn(ifac)
          subchannel_wallsurf(i) = subchannel_wallsurf(i)+surfbn(ifac)
          nf(i) = nf(i)+1
        end if
      end do
    end do

!! Calculate total volume of each subchannl (excluding solid)
    subchannel_vol = 0.d0
    do iel=1, ncel
      do i=1, n_cell
        if(int(subchannel_ID(iel))==i)then
          subchannel_vol(i) = subchannel_vol(i)+volume(iel)
          nc(i) = nc(i)+1
        end if
      end do
    end do

!! In case of parallel run
    if (irangp.ge.0) then
      call parrsm(n_cell, subchannel_surf)
      call parrsm(n_cell, subchannel_rodsurf)
      call parrsm(n_cell, subchannel_wallsurf)
      call parrsm(n_cell, subchannel_vol)
      call parism(n_cell, nc)
      call parism(n_cell, nf)
    endif

!! Calculate the hydraulic diameter for each subchannel 
!! Store total surface area and subchannel volume 
    do iel=1, ncel
      Dh(iel) = 4.d0*subchannel_vol(int(subchannel_ID(iel)))&
              / subchannel_surf(int(subchannel_ID(iel)))
      surf(iel) = subchannel_surf(int(subchannel_ID(iel)))
      surf_rod(iel) = subchannel_rodsurf(int(subchannel_ID(iel)))
      surf_wall(iel) = subchannel_wallsurf(int(subchannel_ID(iel)))
      cvol(iel) = subchannel_vol(int(subchannel_ID(iel)))
    end do

!! Deallocate temporary arrays
    deallocate(lstelt1, lstelt2)
    deallocate(nf, nc)
  
  end subroutine hydraulic_diameter

!===============================================================================
!! This routine is used to calculate the axial Reynolds number
  subroutine subchannel_axial_Re()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: i, ii, iel
    double precision, dimension(:), pointer :: subchannel_ID, Dh, axial_Re
    double precision, dimension(:), pointer :: A_vis, A_den
    double precision, dimension(:,:), pointer :: A_vel
    double precision, allocatable, dimension(:,:) :: subchannel_vel

    call field_get_val_v_by_name('AverageVelocity', A_vel)
    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('HydraulicDiameter', Dh)
    call field_get_val_s_by_name('AxialReynoldsNumber', axial_Re)
    call field_get_val_s_by_name('AverageDynamicViscosity', A_vis)
    call field_get_val_s_by_name('AverageDensity', A_den)

    do iel=1, ncel
      axial_Re(iel) = A_den(iel)*A_vel(3,iel)*Dh(iel)/A_vis(iel)
      axial_Re(iel) = abs(axial_Re(iel))
    end do

  end subroutine subchannel_axial_Re

!===============================================================================
!! This routine is used to calculate the axial frictional factor
  subroutine subchannel_axial_f()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: iel
    double precision, dimension(:), pointer :: axial_Re, axial_f, &
                                               subchannel_ID, A_Tbulk, &
                                               A_twall
    double precision :: temperature_ratio

    call field_get_val_s_by_name('AxialReynoldsNumber', axial_Re)
    call field_get_val_s_by_name('AxialFrictionFactor', axial_f)
    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('SubchannelBulkTemperature', A_Tbulk)
    call field_get_val_s_by_name('SubchannelWallTemperature', A_twall)

    do iel=1, ncel
      select case (bundle_type)
      case (1)
        if(axial_Re(iel)<1.d-30)then
          axial_f(iel) = 1.d30
        else
          axial_f(iel) = FRICTIONAL_FACTOR_PWR(subchannel_ID(iel), &
                                               axial_Re(iel), 1)
          if(user_frictional_factor) &
            axial_f(iel) = FRICTIONAL_FACTOR_USER(iel)
        end if
      case(2)
        temperature_ratio = A_twall(iel)/A_Tbulk(iel)
        if(axial_Re(iel)<1.d-30)then
          axial_f(iel) = 1.d30
        else
          axial_f(iel) = FRICTIONAL_FACTOR_AGR(voidage, temperature_ratio, &
                                       subchannel_ID(iel), axial_Re(iel), 1)
          if(user_frictional_factor) &
            axial_f(iel) = FRICTIONAL_FACTOR_USER(iel)
        end if
      end select
    end do

  end subroutine subchannel_axial_f

!===============================================================================
!! This is for PWR
  function FRICTIONAL_FACTOR_PWR(sub_ID, Re, direction)

    use numvar
    use field
    use mesh
    use parall
    use optcal

    implicit none

    double precision :: FRICTIONAL_FACTOR_PWR
!! Aguments
    double precision :: sub_ID, Re
    integer :: direction

!! Local variables 
    double precision, dimension(3) :: a, b1, b2
    double precision :: n, cf(3)

    integer :: isub

!! Cheng & Todreas (1986)
    if(Re<=2300.d0)then
      a  = (/35.55d0, 44.40d0, 58.83d0/)
      b1 = (/263.7d0, 256.7d0, 160.7d0/)
      b2 = (/-190.2d0, -267.6d0, -203.5d0/)
      n  = 1.d0
    else
      a  = (/0.1339d0, 0.1430d0, 0.1452d0/)
      b1 = (/0.09059d0, 0.04199d0, 0.02681d0/)
      b2 = (/-0.09926d0, -0.04428d0, -0.03411d0/)
      n  = 0.18d0
    end if

    isub = int(sub_ID)

    select case (direction)
    case(1)
      cf = a+b1*(PbyD(isub)-1.d0)+b2*(PbyD(isub)-1.d0)**2
      cf = cf/(Re**n+1.d-15)
!! Checking results
!      if(isub==2)write(*,*)cf(2), Re, PbyD(isub)
    case(2)

    case(3)

    end select

    select case (cell_type(isub))
    case(1)
      FRICTIONAL_FACTOR_PWR = cf(1)
    case(2)
      FRICTIONAL_FACTOR_PWR = cf(2)
    case(3)
      FRICTIONAL_FACTOR_PWR = cf(3)
    case default
      FRICTIONAL_FACTOR_PWR = 1.d30
    end select

  end function FRICTIONAL_FACTOR_PWR

!===============================================================================
!! This is for AGR
  function FRICTIONAL_FACTOR_AGR(voidage, temp_ratio, sub_ID, Re, direction)

    use numvar
    use field
    use mesh
    use parall

    implicit none

    double precision :: FRICTIONAL_FACTOR_AGR
!! Aguments
    double precision :: voidage, temp_ratio, sub_ID, Re
    integer :: direction

!! Local variables 
    double precision :: f0, f, cf
    double precision :: index3 = 0.55

    integer :: isub

    if (voidage<0.7d0 .or. voidage>=0.8d0)then
      if (irangp.le.0) print*,'Warning: voidage is out of range. Check the code'
    end if

!! EDF report: E Romero, TET/REP/0054/94 CEN/DFWG/P(92)4743/B/3, 1994
!! voidage = 0.7
    if(Re<100.d0)then
      if (irangp.le.0) print*, "Warning: cf calculation out of Re lower bound"
!      f0 = 35.6d0/Re
    else if(Re>=100.d0 .and. Re<2000.d0)then
      f0 = 35.6d0/Re  
    else if(Re>=2000.d0 .and. Re<3000.d0)then
      f0 = 0.44035d0*(Re**(-0.4221)) 
    else if(Re>=3000.d0 .and. Re<5.d3)then
      f0 = 0.146626d0*(Re**(-0.28475d0))
    else if(Re>=5.d3 .and. Re<2.d4)then
      f0 = 0.000181*(Re**0.50156d0)
    else if(Re>=2.d4 .and. Re<5.d4)then
      f0 = 0.004937*(Re**(0.16797))
    else if(Re>=5.d4 .and. Re<1.d6)then
      f0 = 0.027498*(Re**(0.009246))
    else
      if (irangp.le.0) print*, "Warning: cf calculation out of Re upper bound"
      f0 = 0.027498*(Re**(0.009246))
    end if

    isub = int(sub_ID)

    select case (direction)
    case(1)
      f = f0/(temp_ratio**index3)
    case(2)

    case(3)

    end select

!!  convert the Fanning frictional factor to Dacy frictional factor
    cf = 4.d0*f
    FRICTIONAL_FACTOR_AGR = cf

  end function FRICTIONAL_FACTOR_AGR

!===============================================================================
!!read in physical properties from data table
  subroutine readin_properties_from_data_table(maxlim, minlim)

    use parall

    implicit none

    !Arguments
    double precision, intent(out) :: maxlim, minlim

    !Local variables
    character(1) :: dummy
    integer :: i
    integer :: stat
    logical :: alive
    integer :: col = 6

    inquire(file=name_data_table, exist=alive)

    if(alive) then
      open(100,file=name_data_table)
    else
      if(irangp.le.0)&
        write(*,*) "Physical property table doesn't exist. Program terminated."
      stop
    end if

    do i=1,5
      read(100,*)
    end do

    n_rows = 0

    do
      read(100,*, iostat=stat) dummy
      if(stat/=0)exit
      n_rows = n_rows+1
    end do

    allocate(user_prop(n_rows,col))

    rewind(100)

    do i=1,5
      read(100,*)
    end do

    do i=1, n_rows
      read (100,*) user_prop(i,:)
    end do

    close(100)

    maxlim = user_prop(n_rows,1)
    minlim = user_prop(1,1)

  end subroutine readin_properties_from_data_table


!===============================================================================
!! This routine is used to calculate the average velocity in each subchannel  
  subroutine subchannel_velocity()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: i, ii, iel, j
    double precision, dimension(:), pointer :: subchannel_ID
    double precision, dimension(:,:), pointer :: A_vel, vel
    double precision, allocatable, dimension(:,:) :: subchannel_vel 
    
    call field_get_val_v(ivarfl(iu), vel)
    call field_get_val_v_by_name('AverageVelocity', A_vel)
    call field_get_val_s_by_name('SubchannelID', subchannel_ID) 
    
    allocate(subchannel_vel(n_cell,3))
    
    subchannel_vel = 0.d0
    do iel=1, ncel
      j = int(subchannel_ID(iel))
      do ii=1,3
        subchannel_vel(j,ii) = subchannel_vel(j,ii)+vel(ii,iel)*volume(iel)
      end do
    end do

    if (irangp.ge.0) then
      do ii=1,3
        call parrsm(n_cell, subchannel_vel(:,ii))
      end do
    endif

    do i=1,n_cell
      subchannel_vel(i,1) = subchannel_vel(i,1)/subchannel_vol(i)
      subchannel_vel(i,2) = subchannel_vel(i,2)/subchannel_vol(i)
      subchannel_vel(i,3) = subchannel_vel(i,3)/subchannel_vol(i)
    end do

    do iel=1, ncel
      j = int(subchannel_ID(iel))
      A_vel(1,iel) = subchannel_vel(j,1)
      A_vel(2,iel) = subchannel_vel(j,2)
      A_vel(3,iel) = subchannel_vel(j,3)
    end do
   
    deallocate(subchannel_vel)
    
  end subroutine subchannel_velocity

!===============================================================================
 !! This routine is used to initialize subchannel tempeatures
  subroutine subchannel_temperatures_init()

    use numvar
    use optcal
    use field
    use mesh
    use parall

    implicit none

    integer i,j,iel
    double precision, dimension(:), pointer :: subchannel_ID
    double precision, dimension(:), pointer :: A_temp, A_Tbulk, A_twall
    double precision, dimension(:), pointer :: tlocal
    double precision, allocatable, dimension(:) :: subchannel_temp 

    if (iscalt.le.0) return

    call field_get_val_s_by_name('SubchannelID', subchannel_ID) 
    call field_get_val_s_by_name('AverageTemperature', A_temp)
    call field_get_val_s_by_name('SubchannelBulkTemperature', A_Tbulk)
    call field_get_val_s_by_name('SubchannelWallTemperature', A_twall)
    call field_get_val_s(ivarfl(isca(iscalt)),tlocal)
    
    allocate(subchannel_temp(n_cell))

    subchannel_temp = 0.d0
    do iel=1, ncel
      j = int(subchannel_ID(iel))
      subchannel_temp(j) = subchannel_temp(j)+tlocal(iel)*volume(iel)
    end do

    if (irangp.ge.0) then
      call parrsm(n_cell, subchannel_temp)
    endif

    do i=1,n_cell
      subchannel_temp(i) = subchannel_temp(i)/subchannel_vol(i)
    end do

    do iel=1, ncel
      j = int(subchannel_ID(iel))
      A_temp(iel) = subchannel_temp(j)
      A_Tbulk(iel) = A_temp(iel)
      A_twall(iel) = A_temp(iel)
    end do

    deallocate(subchannel_temp)

  end subroutine subchannel_temperatures_init


!===============================================================================
!! This routine is used to calculate the average & bulk temprature of each 
!! subchannel
  subroutine subchannel_temperature()

    use numvar
    use optcal
    use field
    use mesh
    use parall

    implicit none

    integer :: i, ii, iel, j
    integer :: ivart
    double precision, dimension(:), pointer :: subchannel_ID
    double precision, dimension(:), pointer :: A_temp, A_Tbulk, temp, den
    double precision, dimension(:,:), pointer :: vel
    double precision, allocatable, dimension(:) :: subchannel_temp, &
                                                   subchannel_Tbulk, &
                                                   subchannel_mass_flux

    if (iscalt.gt.0) then
      ivart = isca(iscalt)
      call field_get_val_s(ivarfl(ivart), temp)
      call field_get_val_s_by_name('AverageTemperature', A_temp)
      call field_get_val_s_by_name('SubchannelBulkTemperature', A_Tbulk)
    else 
      return
    endif

    call field_get_val_v(ivarfl(iu), vel)
    call field_get_val_s(icrom, den) 
    call field_get_val_s_by_name('SubchannelID', subchannel_ID) 

    allocate(subchannel_temp(n_cell), subchannel_Tbulk(n_cell), &
             subchannel_mass_flux(n_cell))

    subchannel_temp = 0.d0
    subchannel_Tbulk = 0.d0
    subchannel_mass_flux = 0.d0
    do iel=1, ncel
      j = int(subchannel_ID(iel))
      subchannel_temp(j) = subchannel_temp(j)+temp(iel)*volume(iel)
      subchannel_Tbulk(j) = subchannel_Tbulk(j)+temp(iel)*volume(iel) &
                          * den(iel)*vel(3,iel)
      subchannel_mass_flux(j) = subchannel_mass_flux(j)+volume(iel) &
                              * den(iel)*vel(3,iel)
    end do

    if (irangp.ge.0) then
      call parrsm(n_cell, subchannel_temp)
      call parrsm(n_cell, subchannel_Tbulk)
      call parrsm(n_cell, subchannel_mass_flux)
    endif

    do i=1,n_cell
      subchannel_temp(i) = subchannel_temp(i)/subchannel_vol(i)
      if(abs(subchannel_mass_flux(i))<1.d-30)then
        subchannel_Tbulk(i) = subchannel_temp(i)
      else
        subchannel_Tbulk(i) = subchannel_Tbulk(i)/subchannel_mass_flux(i)
      end if
    end do

    do iel=1, ncel
      j = int(subchannel_ID(iel))
      A_temp(iel) = subchannel_temp(j)
      A_Tbulk(iel) = subchannel_Tbulk(j)
    end do

    deallocate(subchannel_temp, subchannel_Tbulk, subchannel_mass_flux)
 
  end subroutine subchannel_temperature

!===============================================================================
!! This routine is used to calculate the subchannel averaged physical properties   
  subroutine subchannel_physical_properties()
  
    use numvar
    use optcal
    use cstphy
    use field
    use mesh
    use parall

    implicit none

    integer :: i, ii, iel, j
    integer :: ifcvsl
    double precision, dimension(:), pointer :: subchannel_ID, A_den, A_vis, &
                                               A_cond, A_cp
    double precision, dimension(:), pointer :: den, dynamicvis, cond, cpro_cp
    double precision, allocatable, dimension(:) :: subchannel_den, &
                                                   subchannel_vis, &
                                                   subchannel_cond, &
                                                   subchannel_cp
    
    call field_get_val_s(iprpfl(iviscl), dynamicvis)
    call field_get_val_s(icrom, den)

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('AverageDensity', A_den)
    call field_get_val_s_by_name('AverageDynamicViscosity', A_vis)

    allocate(subchannel_den(n_cell), subchannel_vis(n_cell))

    subchannel_den = 0.d0
    subchannel_vis = 0.d0
    do iel=1, ncel
      j = int(subchannel_ID(iel))
      subchannel_den(j) = subchannel_den(j)+den(iel)*volume(iel)
      subchannel_vis(j) = subchannel_vis(j)+dynamicvis(iel)*volume(iel)
    end do

    if (irangp.ge.0) then
      call parrsm(n_cell, subchannel_den)
      call parrsm(n_cell, subchannel_vis)
    endif

    do i=1,n_cell
      subchannel_den(i) = subchannel_den(i)/subchannel_vol(i)
      subchannel_vis(i) = subchannel_vis(i)/subchannel_vol(i)
    end do

    do iel=1, ncel
      j = int(subchannel_ID(iel))
      A_den(iel) = subchannel_den(j)
      A_vis(iel) = subchannel_vis(j)
    end do

    deallocate(subchannel_den, subchannel_vis)

    if (iscalt.le.0) return

!! note: this is only for the case that thermal variable is tempeature
    call field_get_key_int(ivarfl(isca(iscalt)), kivisl, ifcvsl)
    if (ifcvsl.ge.0) then
      call field_get_val_s(ifcvsl, cond)
    else
      cond => NULL()
    endif

    if (icp.gt.0) then
      call field_get_val_s(iprpfl(icp), cpro_cp)
    else
      cpro_cp => NULL()
    end if

    call field_get_val_s_by_name('AverageThermalConductivity', A_cond)
    call field_get_val_s_by_name('AverageHeatCapacity', A_cp)

    allocate(subchannel_cond(n_cell), subchannel_cp(n_cell))

    subchannel_cond = 0.d0
    subchannel_cp = 0.d0
    do iel=1, ncel
      j = int(subchannel_ID(iel))
      if (ifcvsl.ge.0) then
        subchannel_cond(j) = subchannel_cond(j)+cond(iel)*volume(iel)
      else
        subchannel_cond(j) = subchannel_cond(j)+visls0(iscalt)*volume(iel)
      end if
      
      if (icp.gt.0) then
        subchannel_cp(j) = subchannel_cp(j)+cpro_cp(iel)*volume(iel)
      else
        subchannel_cp(j) = subchannel_cp(j)+cp0*volume(iel)
      end if
    end do

    if (irangp.ge.0) then
      call parrsm(n_cell, subchannel_cond)
      call parrsm(n_cell, subchannel_cp)
    endif

    do i=1,n_cell
      subchannel_cond(i) = subchannel_cond(i)/subchannel_vol(i)
      subchannel_cp(i) = subchannel_cp(i)/subchannel_vol(i)
    end do

    do iel=1, ncel
      j = int(subchannel_ID(iel))
      A_cond(iel) = subchannel_cond(j)
      A_cp(iel) = subchannel_cp(j)
    end do

    deallocate(subchannel_cond, subchannel_cp)


  end subroutine subchannel_physical_properties


!===============================================================================



!*******************************************************************************
!*                          User defined routines                              *
!*******************************************************************************

!===============================================================================
! The follwoing subroutines are used for user defined rod bundles
!===============================================================================

!===============================================================================
!! This routine is used to define subchannel mesh of non-ideal rod bundles
  subroutine subchannel_mesh_user()
  
    use numvar
    use field
    use mesh
    use parall
  
    implicit none
    
!!  Non-changable variables    
    integer :: n_node, n_entity
    double precision, allocatable, dimension(:) ::  x, y, z
    integer, allocatable, dimension(:,:) :: index_cell
    double precision, dimension(:), pointer :: subchannel_ID, &
                                               subchannel_type, &
                                               pitch2diameter

!!  User-defined variables
    integer :: i, iel
    double precision, allocatable, dimension(:,:) :: xx, yy, zz
    double precision :: pitch
    double precision, allocatable, dimension(:) :: balloon_diameter

!!  Non-changable lines
    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    
    call subchannel_mesh_dimensions(subchannel_mesh_name, n_node, n_entity)
    
    allocate(x(n_node),y(n_node),z(n_node))
    allocate(index_cell(n_entity, 8))
    
    call subchannel_mesh_readin_pwr(subchannel_mesh_name, n_node, n_entity, &
                                    x, y, z, index_cell, n_cell)

!! Module variables, should be deallocated after the whole simulation
    allocate(cell_height(n_cell))
    allocate(cell_vol(n_cell))
    allocate(cell_type(n_cell))
    allocate(PbyD(n_cell))
    allocate(x_channel(n_cell),y_channel(n_cell),z_channel(n_cell)) 
    
!! Scale mesh from mm to m
    if(subchannel_mesh_scale)then
      x = x*subchannel_scaling_factor
      y = y*subchannel_scaling_factor
      z = z*subchannel_scaling_factor
    end if

!-------------------------------------------------------------------------------
!! INSERT USER DEFIENED LINES FROM HERE!
!-------------------------------------------
!! This section is for a distort PWR fuel bundle (one column of rods shifted)
    if(.false.)then
    if(irangp.le.0) &
      write(*,*) 'Checking user defined routine for subchannel mesh ... '
!! Calculate minimum and maximum coordinates value of the subchannel cell
    allocate(xx(n_cell, 2), yy(n_cell, 2), zz(n_cell, 2))
    do i=1, n_cell
      xx(i,1) = minval(x(index_cell(i,:)))
      xx(i,2) = maxval(x(index_cell(i,:)))
      yy(i,1) = minval(y(index_cell(i,:)))
      yy(i,2) = maxval(y(index_cell(i,:)))
      zz(i,1) = minval(z(index_cell(i,:)))
      zz(i,2) = maxval(z(index_cell(i,:)))
      
!! Calculate cell height and cell volume of the subchannel mesh
      cell_height(i) = abs(zz(i,2)-zz(i,1))
      cell_vol(i) = (xx(i,2)-xx(i,1))*(yy(i,2)-yy(i,1))*(zz(i,2)-zz(i,1))
      
!! Calculate cell center coordinates of the subchannel cell
      x_channel(i) = 0.5d0*(xx(i,2)+xx(i,1))
      y_channel(i) = 0.5d0*(yy(i,2)+yy(i,1))
      z_channel(i) = 0.5d0*(zz(i,2)+zz(i,1))
    end do
    
!! cell_type: 1 interior channel; 2 edge channel; 3 corner channel
!! This simple method may be not good in other cases
    cell_type = 1
    where(cell_vol<0.7*maxval(cell_vol))cell_type = 2
    where(cell_vol<0.35*maxval(cell_vol))cell_type = 3
    
!! Set pitch-to-diameter ratio for the subchannels (for non-uniform rod 
!! distribution)
    where(cell_type==1)PbyD(:) = (xx(:,2)-xx(:,1)+yy(:,2)-yy(:,1))*0.5d0 &
                               / rod_diameter
    where(cell_type==2)PbyD(:) = (((min(xx(:,2)-xx(:,1),yy(:,2)-yy(:,1)) &
                               - 0.5d0*rod_diameter)+rod_diameter) &
                               / rod_diameter*2.d0+max(xx(:,2)-xx(:,1), &
                                 yy(:,2)-yy(:,1))/rod_diameter)/3.d0
    where(cell_type==3)PbyD(:) = ((xx(:,2)-xx(:,1)+yy(:,2)-yy(:,1))*0.5d0 &
                               - 0.5d0*rod_diameter+rod_diameter)/rod_diameter

!! Find the subchannel index for each control volume
    do iel=1, ncel
      do i=1, n_cell
        if(xyzcen(1,iel)>=xx(i,1).and.xyzcen(1,iel)<=xx(i,2).and.&
           xyzcen(2,iel)>=yy(i,1).and.xyzcen(2,iel)<=yy(i,2).and.&
           xyzcen(3,iel)>=zz(i,1).and.xyzcen(3,iel)<=zz(i,2))then
          subchannel_ID(iel) = i
        end if
      end do
    end do

    deallocate(xx,yy,zz)

    if(irangp.le.0) &
      write(*,*) 'User defined subchannel mesh completed.'

    end if
!-------------------------------------------
!! This section is for a locally ballooned PWR fuel bundle.
    if(.true.)then
    if(irangp.le.0) &
      write(*,*) 'Checking user defined routine for subchannel mesh ... '
!! Calculate minimum and maximum coordinates value of the subchannel cell
    allocate(xx(n_cell, 2), yy(n_cell, 2), zz(n_cell, 2))
    allocate(balloon_diameter(n_cell))
    do i=1, n_cell
      xx(i,1) = minval(x(index_cell(i,:)))
      xx(i,2) = maxval(x(index_cell(i,:)))
      yy(i,1) = minval(y(index_cell(i,:)))
      yy(i,2) = maxval(y(index_cell(i,:)))
      zz(i,1) = minval(z(index_cell(i,:)))
      zz(i,2) = maxval(z(index_cell(i,:)))
      
!! Calculate cell height and cell volume of the subchannel mesh
      cell_height(i) = abs(zz(i,2)-zz(i,1))
      cell_vol(i) = (xx(i,2)-xx(i,1))*(yy(i,2)-yy(i,1))*(zz(i,2)-zz(i,1))

!! Calculate cell center coordinates of the subchannel cell
      x_channel(i) = 0.5d0*(xx(i,2)+xx(i,1))
      y_channel(i) = 0.5d0*(yy(i,2)+yy(i,1))
      z_channel(i) = 0.5d0*(zz(i,2)+zz(i,1))
    end do

!! Calculate rod diameter for the ballooned section
    balloon_diameter = rod_diameter
    where(z_channel>0.3238d0 .and. z_channel<0.3492d0) &
      balloon_diameter = 2.d0*(7.07d0-5.d0)/25.4d0*(z_channel-0.3238d0) &
                       + rod_diameter
    where(z_channel>=0.3492d0 .and. z_channel<=0.3746d0) &
      balloon_diameter = 2.d0*7.07d-3
    where(z_channel>0.3746d0 .and. z_channel<0.4d0) &
      balloon_diameter = -2.d0*(7.07d0-5.d0)/25.4d0*(z_channel-0.4d0) &
                       + rod_diameter

!! cell_type: 1 interior channel; 2 edge channel; 3 corner channel
!! Irregular cell type: 4 center blocked; 5 side blocked; 6 corner blocked
    cell_type = 1
    where(cell_vol/cell_height<0.95*maxval(cell_vol/cell_height))cell_type = 2
    where(cell_vol/cell_height<0.65*maxval(cell_vol/cell_height))cell_type = 3
    pitch = 0.0137d0
    where(x_channel>0.d0 .and. x_channel<pitch .and. &
          y_channel>0.d0 .and. y_channel<pitch .and. &
          z_channel>0.3238d0 .and. z_channel<0.4d0) cell_type = 4
    where(x_channel>pitch .and. x_channel<2.d0*pitch .and. &
          y_channel>0.d0 .and. y_channel<pitch .and. &
          z_channel>0.3238d0 .and. z_channel<0.4d0) cell_type = 5
    where(x_channel>0.d0 .and. x_channel<pitch .and. &
          y_channel>pitch .and. y_channel<2.d0*pitch .and. &
          z_channel>0.3238d0 .and. z_channel<0.4d0) cell_type = 5
    where(x_channel>pitch .and. x_channel<2.d0*pitch .and. &
          y_channel>pitch .and. y_channel<2.d0*pitch .and. &
          z_channel>0.3238d0 .and. z_channel<0.4d0) cell_type = 6

!! Set pitch-to-diameter ratio for the subchannels 
    where(cell_type==1)PbyD = PbyD_interior
    where(cell_type==2)PbyD = PbyD_edge
    where(cell_type==3)PbyD = PbyD_corner
    where(cell_type==4)PbyD = pitch/balloon_diameter
    where(cell_type==5)PbyD = pitch/(balloon_diameter*0.25d0 &
                            + 0.5d0*(balloon_diameter+rod_diameter)*0.5d0 &
                            + rod_diameter*0.25d0)
    where(cell_type==6)PbyD = pitch/(rod_diameter*0.5d0 &
                            + 0.5d0*(balloon_diameter+rod_diameter)*0.5d0)

!! Find the subchannel index for each control volume
    do iel=1, ncel
      do i=1, n_cell
        if(xyzcen(1,iel)>=xx(i,1).and.xyzcen(1,iel)<=xx(i,2).and.&
           xyzcen(2,iel)>=yy(i,1).and.xyzcen(2,iel)<=yy(i,2).and.&
           xyzcen(3,iel)>=zz(i,1).and.xyzcen(3,iel)<=zz(i,2))then
          subchannel_ID(iel) = i
        end if
      end do
    end do

    deallocate(xx,yy,zz)
    deallocate(balloon_diameter)
    
    if(irangp.le.0) &
      write(*,*) 'User defined subchannel mesh completed.'

    end if
!------------------------------------------------------------------------------- 

    call field_get_val_s_by_name('PitchtoDiameterRatio', pitch2diameter)
    call field_get_val_s_by_name('SubchannelType', subchannel_type)
    
    do iel=1, ncel
      pitch2diameter(iel) = PbyD(int(subchannel_ID(iel)))
      subchannel_type(iel) = cell_type(int(subchannel_ID(iel)))
    end do 

    deallocate(x,y,z)
    deallocate(index_cell)

  end subroutine subchannel_mesh_user

!===============================================================================
!! This routine is used to define user frictional factor 
  function FRICTIONAL_FACTOR_USER(iel)

    use numvar
    use field
    use mesh
    use parall

    implicit none

    double precision :: FRICTIONAL_FACTOR_USER

!! Aguments
    integer :: iel

!! Local variables
    double precision, dimension(:), pointer :: axial_Re,         &
                                               subchannel_ID

    call field_get_val_s_by_name('AxialReynoldsNumber', axial_Re)
    call field_get_val_s_by_name('SubchannelID', subchannel_ID)

!!  User example 1: locally ballooned PWR bundle
    if(.true.) &
      FRICTIONAL_FACTOR_USER = friction_user_1(subchannel_ID(iel), &
                                               axial_Re(iel))
    
  end function FRICTIONAL_FACTOR_USER

!-------------------------------------------------------------------------------
!! This function is used for a locally ballooned PWR bundle
  function friction_user_1(sub_ID, Re)

    implicit none

    double precision :: friction_user_1

!! Aguments
    double precision :: sub_ID, Re

!! Local variables 
    double precision, dimension(3) :: a, b1, b2
    double precision :: n, cf(3)
    double precision :: cf_tube

    integer :: isub

  !! Cheng & Todreas (1986)
    if(Re<=2300.d0)then
      a  = (/35.55d0, 44.40d0, 58.83d0/)
      b1 = (/263.7d0, 256.7d0, 160.7d0/)
      b2 = (/-190.2d0, -267.6d0, -203.5d0/)
      n  = 1.d0
    else
      a  = (/0.1339d0, 0.1430d0, 0.1452d0/)
      b1 = (/0.09059d0, 0.04199d0, 0.02681d0/)
      b2 = (/-0.09926d0, -0.04428d0, -0.03411d0/)
      n  = 0.18d0
    end if

    isub = int(sub_ID)

    cf = a+b1*(PbyD(isub)-1.d0)+b2*(PbyD(isub)-1.d0)**2
    cf = cf/(Re**n+1.d-15)
!!  Checking results
!   if(isub==2)write(*,*)cf(2), Re, PbyD(isub)

    if(Re<=2300.d0)then
!!    Laminar frictional factor
      cf_tube = 64.d0/Re
    else if(Re>2300.d0 .and. Re<3.d4)then
!!    Blasius relation for turbulent frictional factor
      cf_tube = 0.316d0*Re**(-0.25d0)
    else if(Re>=3.d4 .and. Re<1.d6)then
!!    McAdams relation for turbulent frictional factor
      cf_tube = 0.184d0*Re**(-0.2d0)
    else
      cf_tube = 0.184d0*Re**(-0.2d0)
      write(*,*) 'Warning: Re number beyond the valid range for frictional & 
                  & factor calculation.'
    end if

    select case (cell_type(isub))
    case(1)
      friction_user_1 = cf(1)
    case(2)
      friction_user_1 = cf(2)
    case(3)
      friction_user_1 = cf(3)
    case(4)
!!    contact point for this type of subchannel
      if(PbyD(isub)>1.d0) then
        friction_user_1 = cf(1)
      else
        friction_user_1 = cf_tube
      end if
    case(5)
!!    contact point for this type of subchannel
      if(PbyD(isub)>1.156118d0) then
        friction_user_1 = cf(1)
      else
        friction_user_1 = cf(1)
      end if
    case(6)
!!    contact point for this type of subchannel
      if(PbyD(isub)>1.254005d0) then
        friction_user_1 = cf(1)
      else
        friction_user_1 = cf(1)
      end if
    case default
      friction_user_1 = 1.d30
    end select

  end function friction_user_1

!===============================================================================
!! This routine is used to implement user defined boudnary conditions for the
!! coupling interfaces (The current implementation works well for the 2-scale 
!! wall function)
  subroutine interface_user_boundary()

    use numvar
    use field
    use mesh
    use parall
    use optcal
    use cstphy

    implicit none

    integer :: i, ii, iel, ifac, j
    integer :: ilelt, nlelt
    integer, allocatable, dimension(:) :: lstelt
    integer :: ivar, f_id, ifcvsl
    double precision, dimension(:,:), pointer :: A_cplvel, cplvel
    double precision, dimension(:), pointer :: A_den, A_vis, subchannel_ID, &
                                               w_dist_solver,Dh,subchannel_vol,&
                                               cpro_viscl, subchannel_type
    double precision, dimension(:), pointer :: viscl, visct, crom, viscls, &
                                               cpro_cp
    double precision, dimension(:), pointer :: cplt, A_cplt, cplk, cplep, &
                                               cplomg, bld
    double precision, allocatable, dimension(:) :: w_dist
    double precision, allocatable, dimension(:,:) :: A_cplvel_base
    double precision, allocatable, dimension(:) :: A_cplt_base
    double precision, allocatable, dimension(:) :: wall_heat_flux
    integer, allocatable, dimension(:) :: iw_index
    double precision :: yplus, uplus, kplus, eplus, tplus
    double precision :: utau, tau_w, axial_Re, axial_f, reyfro, tstar
    double precision :: cpp, rkl, prdtl, visclc, visctc, turb_schmidt

    double precision, dimension(:), pointer :: w_dist_check, w_index_check, &
                                               wall_heat_flux_check, &
                                               tplus_check, tstar_check, &
                                               prdtl_check

    logical,parameter :: analytical_wall_distance = .true.
    save :: w_dist, iw_index, wall_heat_flux

    !! Get the target coupling interface
    allocate(lstelt(nfabor))
    call getfbr('interface_inlet', nlelt, lstelt)

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('SubchannelType', subchannel_type)
    call field_get_val_s_by_name('HydraulicDiameter', Dh)
    call field_get_val_s_by_name('AverageDensity', A_den)
    call field_get_val_s_by_name('AverageDynamicViscosity', A_vis)
    call field_get_val_s_by_name('SubchannelVolume', subchannel_vol)
    call field_get_val_v_by_name('InterfaceVelocity', cplvel)
    call field_get_val_v_by_name('InterfaceBulkVelocity', A_cplvel)
    call field_get_val_s_by_name('InterfaceKineticEnergy', cplk)
    call field_get_val_s_by_name('InterfaceEpsilon', cplep)
    call field_get_val_s_by_name('InterfaceOmega', cplomg)
    call field_get_val_s(iviscl, cpro_viscl)
    call field_get_val_s_by_name('BlendingFactor', bld)
    bld = 0.5d0

    call interface_subchannel_velocity()

    allocate(A_cplvel_base(3, nfabor))
    
    do ilelt=1,nlelt
      ifac = lstelt(ilelt)
      bld(ifac) = 1.0d0
      do ii=1,3
        A_cplvel_base(ii,ifac) = A_cplvel(ii,ifac)
      end do
    end do

    if(iscalt.gt.0) then
      allocate(A_cplt_base(nfabor))

      call field_get_val_s_by_name('InterfaceTemperature', cplt)
      call field_get_val_s_by_name('InterfaceBulkTemperature', A_cplt)

      call interface_subchannel_temperature()

      do ilelt=1,nlelt
        ifac = lstelt(ilelt)
        A_cplt_base(ifac) = A_cplt(ifac)
      end do

      ivar = isca(iscalt)
      f_id = ivarfl(ivar)

      call field_get_val_s(iviscl, viscl)
      call field_get_val_s(ivisct, visct)

      call field_get_key_int (f_id, kivisl, ifcvsl)

      if (ifcvsl .ge. 0) then
        call field_get_val_s(ifcvsl, viscls)
      endif

      call field_get_val_s(icrom, crom)

      if (icp.ge.0) then
        call field_get_val_s(icp, cpro_cp)
      endif

      call field_get_key_double(ivarfl(ivar), ksigmas, turb_schmidt)

    end if

    if(ntcabs==ntpabs+1) then
      allocate(w_dist(nfabor))
      allocate(iw_index(nfabor))
      allocate(wall_heat_flux(nfabor))

      call field_get_val_s_by_name('wall_distance', w_dist_solver)
      call field_get_val_s_by_name('WallDistanceCheck', w_dist_check)
      call field_get_val_s_by_name('WallIndexCheck', w_index_check)
      call field_get_val_s_by_name('WallHeatFluxCheck', wall_heat_flux_check)

      if(analytical_wall_distance)then

        call analytical_w_dist(w_dist, iw_index)

        do ilelt=1,nlelt
          ifac = lstelt(ilelt)
          if(iw_index(ifac)==1)then
            wall_heat_flux(ifac) = -8.68453d5
          else if(iw_index(ifac)==2)then
            wall_heat_flux(ifac) = -6.76923d5
          else if(iw_index(ifac)==3)then
            !wall_heat_flux(ifac) = 0.d0
            wall_heat_flux(ifac) = -6.76923d5
          end if
          w_dist_check(ifac) = w_dist(ifac)
          w_index_check(ifac) = real(iw_index(ifac))
          wall_heat_flux_check(ifac) = wall_heat_flux(ifac)
        end do

      else
        do ilelt=1,nlelt
          ifac = lstelt(ilelt)
          iel = ifabor(ifac)
          w_dist(ifac) = w_dist_solver(iel)
        end do
      end if

    end if

    call field_get_val_s_by_name('TplusCheck', tplus_check)
    call field_get_val_s_by_name('TstarCheck', tstar_check)
    call field_get_val_s_by_name('PrdtlCheck', prdtl_check)

    do ilelt=1,nlelt
      ifac = lstelt(ilelt)
      iel = ifabor(ifac)
      axial_Re = abs(A_den(iel)*A_cplvel(3,ifac)*Dh(iel)/A_vis(iel))
      axial_f = FRICTIONAL_FACTOR_PWR(subchannel_ID(iel), axial_Re, 1)
      tau_w = 0.25d0*axial_f*0.5d0*A_den(iel)*(A_cplvel(3,ifac)**2.)
      utau = sqrt(tau_w/A_den(iel))
      yplus = A_den(iel)*utau*w_dist(ifac)/A_vis(iel)
      reyfro = A_den(iel)*utau*Dh(iel)/2.d0/A_vis(iel)
      call Reichardt_law (yplus, reyfro, uplus, kplus, eplus)
      cplvel(3,ifac) = uplus*utau
      cplk(ifac) = 0.75d0*kplus*utau**2.
!      cplk(ifac) = kplus*utau**2
!      cplk(ifac) = tau_w/A_den(iel)/0.3d0+2.d0*w_dist(ifac)/Dh(iel) &
!                 * (0.002d0*(8.d0/7.d0*A_cplvel(3,iel))**2 &
!                 - tau_w/A_den(iel)/0.3d0)
      cplep(ifac) = A_den(iel)*eplus*utau**4./cpro_viscl(iel)
      cplomg(ifac) = cplep(ifac)/cmu/max(1.0d-10, cplk(ifac))

      if(iscalt.gt.0)then
        visclc = viscl(iel)
        visctc = visct(iel)
        cpp = 1.d0
        if (itherm.eq.1) then
          if (icp.ge.0) then
            cpp = cpro_cp(iel)
          else
            cpp = cp0
          endif
        endif

        if (ifcvsl.lt.0) then
          rkl = visls0(iscalt)
          prdtl = cpp*visclc/rkl
        else
          rkl = viscls(iel)
          prdtl = cpp*visclc/rkl
        endif

        call temperature_wallfun(yplus, prdtl, turb_schmidt, tplus)

        tstar = wall_heat_flux(ifac)/A_den(iel)/cpp/utau
        cplt(ifac) = tstar*tplus
        
        tplus_check(ifac) = tplus
        tstar_check(ifac) = tstar
        prdtl_check(ifac) = prdtl

      end if
    end do

!!  Calibrate the velocity and temperature based on A_cplvel and A_cplt
    call interface_subchannel_velocity()
    
    do ilelt=1, nlelt
      ifac = lstelt(ilelt)
      do ii=1,3
        cplvel(ii,ifac) = A_cplvel_base(ii,ifac)/(A_cplvel(ii,ifac)+1.d-30) &
                        * cplvel(ii,ifac)
      end do
    end do
    
    if(iscalt.gt.0) then

      call interface_subchannel_temperature()

      do ilelt=1, nlelt
        ifac = lstelt(ilelt)
        cplt(ifac) = A_cplt_base(ifac)-A_cplt(ifac)+cplt(ifac)
      end do
      
    end if

    deallocate(lstelt)
    deallocate(A_cplvel_base)
    if(iscalt.gt.0) deallocate(A_cplt_base)


  contains

  subroutine Reichardt_law (yplus, reyfro, uplus, kplus, eplus)

    implicit none

    double precision, intent(IN) :: yplus, reyfro
    double precision, intent(OUT) :: uplus, kplus, eplus
    double precision, parameter :: xkappa = 0.41d0

    uplus = log(1.d0+0.4d0*yplus)/xkappa &
          + 7.8d0*(1.d0-exp(-yplus/11.d0)-yplus/11.d0*exp(-0.33d0*yplus))

    kplus = 0.07d0*yplus*yplus*exp(-yplus/8.d0) &
          + (1.d0 - exp(-yplus/20.d0))*4.5d0/(1.d0 + 4.d0*yplus/reyfro)

!    kplus = 4.2d0*exp(-yplus/360.d0)

    eplus = (1.d0/xkappa)/(yplus**4+15.d0**4)**(0.25d0)

  end subroutine Reichardt_law

  subroutine temperature_wallfun(yplus, prl, prt, tplus)

    use numvar
    use field
    use mesh
    use parall
    use optcal
 
    implicit none

    double precision, intent(IN) :: yplus
    double precision, intent(IN) :: prt, prl
    double precision, intent(OUT) :: tplus

    double precision :: yp1, yp2, a2, beta2
    double precision, parameter :: xkappa = 0.41d0

    yp2 = xkappa*1000.d0/prt
    yp2 = sqrt(yp2)
    yp1 = (1000.d0/prl)**(1./3.)

    a2 = 15.d0*prl**(2./3.)
    beta2 = a2-500.d0/(yp2**2.)

    if (yplus >= yp1 .and. yplus < yp2) &
      tplus = a2-500.d0/(yplus*yplus)

    if (yplus >= yp2) &
      tplus = beta2+prt/xkappa*log(yplus/yp2)

  end subroutine temperature_wallfun

  subroutine analytical_w_dist(w_dist, iw_index)

    implicit none

    double precision, dimension(:) :: w_dist
    integer, dimension(:) :: iw_index

    double precision, dimension(4) :: x, y
    double precision :: xmin, ymin, xmax, ymax
    integer, dimension(4) :: iw_type
    double precision, dimension(4) :: distance
    integer :: node_index

    xmin = minval(xx(:,1))
    xmax = maxval(xx(:,2))
    ymin = minval(yy(:,1))
    ymax = maxval(yy(:,2))

    w_dist = 0.d0
    iw_index = 0
    do ilelt=1,nlelt
      ifac = lstelt(ilelt)
      iel = ifabor(ifac)
      j = int(subchannel_ID(iel))
      x(1) = xx(j,1)
      y(1) = yy(j,1)
      x(2) = xx(j,2)
      y(2) = yy(j,1)
      x(3) = xx(j,2)
      y(3) = yy(j,2)
      x(4) = xx(j,1)
      y(4) = yy(j,2)

      if(subchannel_type(iel)==1) then
        iw_type = 0
        distance(:) = sqrt((xyzcen(1,iel)-x(:))**2+(xyzcen(2,iel)-y(:))**2)
        w_dist(ifac) = minval(distance)-0.5d0*rod_diameter
        node_index = minloc(distance,1)
        if(xyzcen(1,iel)>0.01415d0 .and. xyzcen(2,iel)>0.01415d0 .and. &
           xyzcen(1,iel)<0.05195d0 .and. xyzcen(2,iel)<0.05195d0) then
          iw_index(ifac) = 1
        else
          iw_index(ifac) = 2
        end if

      else if (subchannel_type(iel)==2) then
        iw_type = 0
        do i=1,4
          if(abs(x(i)-xmin)<1.d-10 .or. abs(x(i)-xmax)<1.d-10) iw_type(i)=1
          if(abs(y(i)-ymin)<1.d-10 .or. abs(y(i)-ymax)<1.d-10) iw_type(i)=2
        end do
        do i=1,4
          if(iw_type(i)==0)then
            distance(i) = sqrt((xyzcen(1,iel)-x(i))**2+(xyzcen(2,iel)-y(i))**2)
            distance(i) = distance(i)-0.5d0*rod_diameter
          else if(iw_type(i)==1)then
            distance(i) = abs(xyzcen(1,iel)-x(i))
          else if(iw_type(i)==2)then
            distance(i) = abs(xyzcen(2,iel)-y(i))
          end if
        end do
        w_dist(ifac) = minval(distance)
        node_index = minloc(distance,1)
        if(iw_type(node_index)==0)then 
          iw_index(ifac) = 2
        else
          iw_index(ifac) = 3
        end if

      else if (subchannel_type(iel)==3) then
        iw_type = 0
        do i=1,4
          if(abs(x(i)-xmin)<1.d-10 .or. abs(x(i)-xmax)<1.d-10) iw_type(i)=1
          if(abs(y(i)-ymin)<1.d-10 .or. abs(y(i)-ymax)<1.d-10) iw_type(i)=2
          if(abs(x(i)-xmin)<1.d-10 .and. abs(y(i)-ymin)<1.d-10) iw_type(i)=3
          if(abs(x(i)-xmin)<1.d-10 .and. abs(y(i)-ymax)<1.d-10) iw_type(i)=3
          if(abs(x(i)-xmax)<1.d-10 .and. abs(y(i)-ymin)<1.d-10) iw_type(i)=3
          if(abs(x(i)-xmax)<1.d-10 .and. abs(y(i)-ymax)<1.d-10) iw_type(i)=3
        end do

        do i=1,4
          if(iw_type(i)==0)then
            distance(i) = sqrt((xyzcen(1,iel)-x(i))**2+(xyzcen(2,iel)-y(i))**2)
            distance(i) = distance(i)-0.5d0*rod_diameter
          else if(iw_type(i)==1)then
            distance(i) = abs(xyzcen(1,iel)-x(i))
          else if(iw_type(i)==2)then
            distance(i) = abs(xyzcen(2,iel)-y(i))
          else if(iw_type(i)==3)then
            distance(i) = sqrt((xyzcen(1,iel)-x(i))**2+(xyzcen(2,iel)-y(i))**2)
          end if
        end do
        w_dist(ifac) = minval(distance)
        node_index = minloc(distance,1)
        if(iw_type(node_index)==0)then 
          iw_index(ifac) = 2
        else
          iw_index(ifac) = 3
        end if

      end if

    end do

  end subroutine analytical_w_dist

  subroutine location_square(x,y,x1,y1,x2,y2,within)

    implicit none

    double precision, intent(IN) :: x,y
    double precision, intent(IN) :: x1,y1,x2,y2
    logical, intent(OUt) :: within

    double precision :: rx, ry, rx1, ry1
    logical :: condition1, condition2

    condition1 = .false.
    condition2 = .false.

    rx = x-x1
    rx1 = x2-x1
    ry = y-y1
    ry1 = y2-y1

    if(rx*rx1>0.d0 .and. abs(rx1)>abs(rx)) condition1 = .true.
    if(ry*ry1>0.d0 .and. abs(ry1)>abs(ry)) condition2 = .true.
    if(condition1 .and. condition2) then
      within = .true.
    else
      within = .false.
    end if

  end subroutine location_square

  end subroutine interface_user_boundary

!-------------------------------------------------------------------------------
!! This routine is called in interface_user_boundary() to calculate subchannel
!! bulk velocity at the coupling interfaces
  subroutine interface_subchannel_velocity()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: i, ii, iel, ifac, j
    integer :: ilelt, nlelt
    integer, allocatable, dimension(:) :: lstelt
    double precision, dimension(:), pointer :: subchannel_ID, subchannel_vol
    double precision, dimension(:,:), pointer :: A_cplvel, cplvel
    double precision, allocatable, dimension(:,:) :: subchannel_cplvel

    call field_get_val_v_by_name('InterfaceVelocity', cplvel)
    call field_get_val_v_by_name('InterfaceBulkVelocity', A_cplvel)
    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('SubchannelVolume', subchannel_vol) 

    !! Get all wall faces
    allocate(lstelt(nfabor))
    call getfbr('interface_inlet', nlelt, lstelt)
    allocate(subchannel_cplvel(n_cell,3))

    subchannel_cplvel = 0.d0
    do ilelt=1, nlelt
      ifac = lstelt(ilelt)
      iel = ifabor(ifac)
      j = int(subchannel_ID(iel))
      do ii=1,3
        subchannel_cplvel(j,ii) = subchannel_cplvel(j,ii) &
                                + cplvel(ii,ifac)*volume(iel)
      end do
    end do

    if (irangp.ge.0) then
      do ii=1,3
        call parrsm(n_cell, subchannel_cplvel(:,ii))
      end do
    endif

    do ilelt=1, nlelt
      ifac = lstelt(ilelt)
      iel = ifabor(ifac)
      j = int(subchannel_ID(iel))
      do ii=1,3
        A_cplvel(ii,ifac) = subchannel_cplvel(j,ii)/subchannel_vol(iel)
      end do
    end do

    deallocate(subchannel_cplvel)
    deallocate(lstelt)

  end subroutine interface_subchannel_velocity

!-------------------------------------------------------------------------------
!! This routine is called in interface_user_boundary() to calculate subchannel
!! bulk temperature at the coupling interfaces (This should be refined for
!! cases with significant variable density)
  subroutine interface_subchannel_temperature()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: i, iel, ifac, j
    integer :: ilelt, nlelt
    integer, allocatable, dimension(:) :: lstelt
    double precision, dimension(:), pointer :: subchannel_ID, subchannel_vol, &
                                               A_cplt, cplt
    double precision, dimension(:,:), pointer :: cplvel
    double precision, allocatable, dimension(:) :: subchannel_cplt,&
                                                   subchannel_cplt_w,&
                                                   subchannel_w

    call field_get_val_v_by_name('InterfaceVelocity', cplvel)
    call field_get_val_s_by_name('InterfaceTemperature', cplt)
    call field_get_val_s_by_name('InterfaceBulkTemperature', A_cplt)
    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('SubchannelVolume', subchannel_vol)

    !! Get all wall faces
    allocate(lstelt(nfabor))
    call getfbr('interface_inlet', nlelt, lstelt)
    allocate(subchannel_cplt(n_cell))
    allocate(subchannel_cplt_w(n_cell))
    allocate(subchannel_w(n_cell))

    subchannel_cplt = 0.d0
    subchannel_cplt_w = 0.d0
    subchannel_w = 0.d0
    do ilelt=1, nlelt
      ifac = lstelt(ilelt)
      iel = ifabor(ifac)
      j = int(subchannel_ID(iel))
      subchannel_cplt_w(j) = subchannel_cplt_w(j) &
                           + cplt(ifac)*cplvel(3,ifac)*volume(iel)
      subchannel_w(j) = subchannel_w(j) &
                      + cplvel(3,ifac)*volume(iel)
    end do

    if (irangp.ge.0) then
      call parrsm(n_cell, subchannel_cplt_w)
      call parrsm(n_cell, subchannel_w)
    endif

    do i=1, n_cell
      subchannel_cplt(i) = subchannel_cplt_w(i)/subchannel_w(i)
    end do

    do ilelt=1, nlelt
      ifac = lstelt(ilelt)
      iel = ifabor(ifac)
      j = int(subchannel_ID(iel))
      A_cplt(ifac) = subchannel_cplt(j)
    end do

    deallocate(subchannel_cplt, subchannel_cplt_w, subchannel_w)
    deallocate(lstelt)

  end subroutine interface_subchannel_temperature

!===============================================================================
! The follwoing subroutines are used for post-processings
!===============================================================================

!===============================================================================
!! This routine is used to extract rod wall tempeatures, y plus and shear stress
!! along the circumfrential direction at a specified axial location (The 
!! following routine is currently for the 5x5 PWR bundle and should be adjusted 
!! according to different data processing requirements)
!! ----------------------------
!! Rod IDs
!! | (11) (12) (13) (14) (15) | 
!! |                          |
!! | (21) (22) (23) (24) (25) |
!! |                          |
!! | (31) (32) (33) (34) (35) |
!! |                          |
!! | (41) (42) (43) (44) (45) |
!! |                          | 
!! | (51) (52) (53) (54) (55) |
!! ----------------------------
  subroutine circumferential_quantities ()
  
    use numvar
    use optcal
    use cstphy
    use field
    use mesh
    use parall  
    use entsor
       
    implicit none
    
    integer :: i, ii, j, jj, k, kk, iel, ilelt, ifac, nlelt
    integer, allocatable, dimension(:) :: lstelt 
    integer :: nm, mm, arr
    integer :: impout
    double precision, dimension(:,:,:), allocatable :: rod_center_xyz
    double precision, dimension(:,:), allocatable :: rod_radius
    double precision, dimension(:), allocatable :: angle, angle_g,           &
                                                   twall_para, twall_para_g, &
                                                   yplus_para, yplus_para_g, &
                                                   tplus_para, tplus_para_g, &
                                                   tstar_para, tstar_para_g, &
                                                   fwall_para, fwall_para_g, &
                                                   nusselt_para, nusselt_para_g
    double precision, dimension(:), pointer :: twall, yplusp, nusseltp, tplusp,&
                                               tstarp, subchannel_ID, cvar_scalt
    double precision, dimension(:,:), pointer :: fwall
    integer :: itwall, ifwall, iyplus, itplus, itstar, inusselt, ivart
    double precision :: x,x1,x2,y,y1,y2,radius,radius_cell,angle0,cell_thichness
    character(len=50) :: chaine, cfilename
    logical :: condition_a, condition_b
    integer :: nvariables
    double precision, dimension(:,:), allocatable :: array
    
    ! User specified parameters
    ! Geometrical parameters
    double precision, parameter :: pitch = 0.0126d0, &
                                   tol1 = 1.d-3, &
                                   tol2 = 1.d-3, &
                                   fuel_rod_diameter = 0.0095d0, &
                                   x11 = 0.00785d0, &
                                   y11 = 0.05825d0
    double precision, dimension(4) ::  axial_loc
    
    axial_loc(1) = -1.050d0
    axial_loc(2) = -1.015d0
    axial_loc(3) = -0.945d0
    axial_loc(4) = -0.875d0

    if (iscalt.le.0) return

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)


    if (iscalt.gt.0) then
      ivart = isca(iscalt)
      call field_get_val_s(ivarfl(ivart), cvar_scalt)
    endif
    
    call field_get_id_try('boundary_temperature', itwall)
    if (itwall.ge.0) call field_get_val_s (itwall, twall)
    
    call field_get_id_try('boundary_forces', ifwall)
    if (ifwall.ge.0) call field_get_val_v (ifwall, fwall)

!!  to enforce the field of yplus, tplus and tstar to exist, the relevant 
!!  options in the postprocessing section of GUI should be ticked or manually 
!!  created using cs_user_parameters.c (LocalNusseltNumber equivalent to
!!  boundary_layer_nusselt)
    call field_get_id_try('yplus', iyplus)
    call field_get_id_try('tplus', itplus)
    call field_get_id_try('tstar', itstar)

    if (iyplus.ge.0) call field_get_val_s (iyplus, yplusp)

    if (itstar.ge.0 .and. itplus.ge.0) then
      call field_get_val_s (itplus, tplusp)
      call field_get_val_s (itstar, tstarp)
    end if

    call field_get_id_try('LocalNusseltNumber', inusselt)
    call field_get_val_s (inusselt, nusseltp)

    allocate(rod_center_xyz(5,5,2), rod_radius(5,5))
!!  Modify the value 100 if the number of the circumferential divisions 
!!  are too large
    arr = 200
    allocate(angle(arr),angle_g(arr),twall_para(arr),twall_para_g(arr), &
             fwall_para(arr), fwall_para_g(arr), yplus_para(arr),       &
             yplus_para_g(arr), tplus_para(arr), tplus_para_g(arr),     &
             tstar_para(arr), tstar_para_g(arr), nusselt_para(arr),     &
             nusselt_para_g(arr))
    
    rod_radius = fuel_rod_diameter/2.d0

    rod_center_xyz(1,1,1) = x11
    rod_center_xyz(1,1,2) = y11
    do j=1,5
      if(j>1)rod_center_xyz(1,j,1)=rod_center_xyz(1,1,1)+(j-1)*pitch
      do i=1,5
        rod_center_xyz(i,j,1) = rod_center_xyz(1,j,1)
        rod_center_xyz(i,j,2) = rod_center_xyz(1,1,2)-(i-1)*pitch
      end do
    end do

    !! Get all wall faces
    allocate(lstelt(nfabor))
    chaine = heated_rod_walls
    call getfbr(chaine, nlelt, lstelt)
    call boundary_nusselt (nlelt, lstelt, nusseltp)

    do kk=1,4
!    do k=1,5
!    do j=1,5

!!  only calculate for the central rod
    do k=3,3
    do j=3,3
    
    nm = 1
    angle = 0.d0 
    do ii=0, nrangp-1
    mm = 1
    if(irangp==ii.or.irangp==-1)then
    do ilelt = 1, nlelt
      ifac = lstelt(ilelt)
      iel = ifabor(ifac)
      jj = subchannel_ID(iel)
      cell_thichness = cell_height(jj)
      x = rod_center_xyz(j,k,1)
      y = rod_center_xyz(j,k,2)
      radius = sqrt((cdgfbo(1,ifac)-x)**2+(cdgfbo(2,ifac)-y)**2)
      radius_cell = sqrt((xyzcen(1,iel)-x)**2+(xyzcen(2,iel)-y)**2)

      !! In the case that axial_loc coincides with grid lines
      condition_a = abs(cdgfbo(3,ifac)-axial_loc(kk)) &
                    < (0.5d0-tol1)*cell_thichness
      condition_b = abs(cdgfbo(3,ifac)-axial_loc(kk)) &
                    < (0.5d0+tol1)*cell_thichness&
                    .and. (cdgfbo(3,ifac)-axial_loc(kk)) > 0.d0
      condition_b = .not.(condition_a) .and. condition_b

      if(radius < rod_radius(j,k)*(1.d0+tol2) .and. &
         radius_cell > rod_radius(j,k) .and.        &
         (condition_a .or. condition_b)) then

        if((cdgfbo(2,ifac)-y) >= 0)then
          angle(mm) = acos((cdgfbo(1,ifac)-x)/radius)/pi*180.d0
        else
          angle(mm) = 360.d0-acos((cdgfbo(1,ifac)-x)/radius)/pi*180.d0 
        end if  
        if(itwall.ge.0) then
          twall_para(mm) = twall(ifac)
        else
          twall_para(mm) = 0.d0
        end if        
        if(ifwall.ge.0) then
          fwall_para(mm) = fwall(3,ifac)/surfbn(ifac)
        else
          fwall_para(mm) = 0.d0
        end if
        if(iyplus.ge.0) then
          yplus_para(mm) = yplusp(ifac)
        else
          yplus_para(mm) = 0.d0
        end if
        if (itstar.ge.0 .and. itplus.ge.0) then
          tplus_para(mm) = tplusp(ifac)
          tstar_para(mm) = tstarp(ifac)
        else
          tplus_para(mm) = 0.d0
          tstar_para(mm) = 0.d0
        end if
        if (itstar.ge.0 .and. itplus.ge.0 .and. inusselt.ge.0)then
          nusselt_para(mm) = nusseltp(ifac)
        else
          nusselt_para(mm) = 0.d0
        end if
        mm = mm+1
      end if
    enddo
!    write(*,*)'irangp',irangp,mm
    endif

    if(irangp.ge.0)call parall_bcast_i(ii, mm)

    do i=1,mm-1
      if(irangp.ge.0)then
        call parall_bcast_r(ii, angle(i))
        call parall_bcast_r(ii, twall_para(i))
        call parall_bcast_r(ii, fwall_para(i))
        call parall_bcast_r(ii, yplus_para(i))
        call parall_bcast_r(ii, tplus_para(i))
        call parall_bcast_r(ii, tstar_para(i))
        call parall_bcast_r(ii, nusselt_para(i))
      endif

      angle_g(nm) = angle(i)
      twall_para_g(nm) = twall_para(i)
      fwall_para_g(nm) = fwall_para(i)
      yplus_para_g(nm) = yplus_para(i)
      tplus_para_g(nm) = tplus_para(i)
      tstar_para_g(nm) = tstar_para(i)
      nusselt_para_g(nm) = nusselt_para(i)
      nm = nm+1

      if(irangp.ge.0)then
        call parall_bcast_i(ii, nm)
        call parall_bcast_r(ii, angle_g(nm))
        call parall_bcast_r(ii, twall_para_g(nm))
        call parall_bcast_r(ii, fwall_para_g(nm))
        call parall_bcast_r(ii, yplus_para_g(nm))
        call parall_bcast_r(ii, tplus_para_g(nm))
        call parall_bcast_r(ii, tstar_para_g(nm))
        call parall_bcast_r(ii, nusselt_para_g(nm))
      endif
    enddo
    
    enddo

    if (ntcabs.eq.ntmabs) then
    impout = impusr(1)
    if (irangp.le.0) then
      write(cfilename,'(a27,2i02.2,a4,I5,a4)')&
                                        'circumferential_quantities_',j,k,&
                                        '_z=',int(axial_loc(kk)*1.d3),'.dat'
      open(impout,file=Trim(cfilename))         
      write(impout,*)  &
           '   degree(o)      t_wall(C)      f_wall(N/m2)   yplus          &
           &tplus          tstar          nusselt_wall'
    endif
    endif

    nvariables = 7
    allocate (array(nm-1,nvariables))

    array(1:nm-1,1) = angle_g(1:nm-1)
    array(1:nm-1,2) = twall_para_g(1:nm-1)
    array(1:nm-1,3) = fwall_para_g(1:nm-1)
    array(1:nm-1,4) = yplus_para_g(1:nm-1)
    array(1:nm-1,5) = tplus_para_g(1:nm-1)
    array(1:nm-1,6) = tstar_para_g(1:nm-1)
    array(1:nm-1,7) = nusselt_para_g(1:nm-1)
    
    call sorting(array)

    angle_g(1:nm-1) = array(1:nm-1,1)
    twall_para_g(1:nm-1) = array(1:nm-1,2)
    fwall_para_g(1:nm-1) = array(1:nm-1,3)
    yplus_para_g(1:nm-1) = array(1:nm-1,4)
    tplus_para_g(1:nm-1) = array(1:nm-1,5)
    tstar_para_g(1:nm-1) = array(1:nm-1,6)
    nusselt_para_g(1:nm-1) = array(1:nm-1,7)

    deallocate(array)

    do i=1, nm-1
      if (ntcabs.eq.ntmabs) then
        if (irangp.le.0) write(impout,'(7Es15.5)') angle_g(i), twall_para_g(i),&
                        fwall_para_g(i), yplus_para_g(i), tplus_para_g(i),     &
                        tstar_para_g(i), nusselt_para_g(i)
      endif  
    end do
    
!! Serial version    
!    do ilelt = 1, nlelt
!      ifac = lstelt(ilelt)
!      x = rod_center_xyz(1,1,1)
!      y = rod_center_xyz(1,1,2)
!      radius = sqrt((cdgfbo(1,ifac)-x)**2+(cdgfbo(2,ifac)-y)**2)
!      if(radius < rod_radius(1,1)*(1.d0+tol2) &
!         .and. cdgfbo(3,ifac) > -0.5d0*thichness)then
!        if((cdgfbo(2,ifac)-y) >= 0)then
!          angle0 = acos((cdgfbo(1,ifac)-x)/radius)/pi*180.d0
!        else
!          angle0 = 360.d0-acos((cdgfbo(1,ifac)-x)/radius)/pi*180.d0 
!        end if  
!      if (irangp.le.0) write(impout,'(2f10.5)')angle0, twall(ifac)    
!      end if    
!     
!    enddo

    if (ntcabs.eq.ntmabs) then
      if (irangp.le.0) close(impout)  
    endif

    end do
    end do
    end do

    deallocate(rod_center_xyz, rod_radius)
    deallocate(angle, angle_g, twall_para, twall_para_g, fwall_para, &
               fwall_para_g, yplus_para, yplus_para_g, tplus_para,   &
               tplus_para_g, tstar_para, tstar_para_g, nusselt_para, &
               nusselt_para_g)
    deallocate(lstelt)

    contains

    subroutine sorting(matrix)

      implicit none

      double precision,dimension(:,:), allocatable :: matrix
      double precision,dimension(:), allocatable :: buf
      integer :: nrow, ncolumn, irow, krow
      integer :: i
       
      nrow = size(matrix(:,1))
      ncolumn = size(matrix(1,:))
      
      allocate (buf(ncolumn))
      
      do irow = 1, nrow
        krow = minloc(matrix(irow:nrow, 1), dim=1)+irow-1
        buf(:) = matrix(irow,:)
        matrix(irow,:) = matrix(krow,:)
        matrix(krow,:) = buf(:)
      enddo 

      deallocate(buf)

      return
      
    end subroutine sorting

  end subroutine circumferential_quantities


!===============================================================================
!! This routine is used to calculate the circumferentially averaged rod wall 
!! temperature and cross-sectional bulk quantities along the axial direction 
!! (The following routine is currently for the 5x5 PWR bundle and should be 
!! adjusted according to different data processing requirements)
!! ----------------------------
!! Rod IDs
!! | (11) (12) (13) (14) (15) | 
!! |                          |
!! | (21) (22) (23) (24) (25) |
!! |                          |
!! | (31) (32) (33) (34) (35) |
!! |                          |
!! | (41) (42) (43) (44) (45) |
!! |                          | 
!! | (51) (52) (53) (54) (55) |
!! ----------------------------
  subroutine axial_avg_quantities ()

    use numvar
    use optcal
    use cstphy
    use field
    use mesh
    use parall
    use entsor

    implicit none

    integer :: i, ii, j, jj, k, iel, ilelt, ifac, nlelt
    integer :: iel1, irang1, irangv
    integer :: npoint
    integer :: ivart
    integer, allocatable, dimension(:) :: lstelt
    integer, allocatable, dimension(:) :: n_faces
    integer :: impout
    double precision, allocatable, dimension(:) :: xabs
    double precision, allocatable, dimension(:) :: xt1, xt2, xt_bulk, xw, xp
    double precision, dimension(:,:,:), allocatable :: rod_center_xyz
    double precision, dimension(:,:), allocatable :: rod_radius
    double precision, dimension(:), pointer :: twall, subchannel_ID
    double precision, dimension(:), pointer :: cvar_p, cvar_scalt, cpro_rom
    double precision, dimension(:,:), pointer :: cvar_vel
    double precision :: x,x1,x2,y,y1,y2,radius,cell_thichness
    double precision :: xyz(3)
    double precision :: sum_surf, sum_surf_t
    double precision :: temp, sum_vol, sum_vol_t, sum_vol_w, sum_vol_p, &
                        sum_vol_rho_w, sum_vol_t_rho_w
    character(len=50) :: chaine, cfilename
    logical :: condition_a

    ! User specified parameters
    ! Geometrical parameters
    double precision, parameter :: pitch = 0.0126d0, &
                                   tol1 = 1.d-3, &
                                   tol2 = 1.d-3, &
                                   fuel_rod_diameter = 0.0095d0, &
                                   x11 = 0.00785d0, &
                                   y11 = 0.05825d0, &
                                   start_location = -1270.d-3, &
                                   end_location = -871.d-3, &
                                   increment = 1.d-3
                                   !increment = 2.d-3 !(mesh used 3 coarse)
    
    if (iscalt.le.0) return

!!  Set sampling axial locations
    iel1   = -999
    irang1 = -999
    npoint = int((end_location-start_location)/increment)
    allocate(xabs(npoint))
    do ii = 1, npoint
      xyz(1) = 0.004d0
      xyz(2) = 0.004d0
      xyz(3) = start_location+float(ii-1)/float(npoint-1)&
             * (end_location-start_location)

      call findpt(ncelet, ncel, xyzcen, xyz(1), xyz(2), xyz(3), iel, irangv)

      if ((iel.ne.iel1).or.(irangv.ne.irang1)) then
        iel1   = iel
        irang1 = irangv
        if (irangp.eq.irangv) xabs(ii) = xyzcen(3,iel)
        if (irangp.ge.0) then
          call parall_bcast_r(irangv, xabs(ii))
        endif
      endif

    enddo

    ivart = isca(iscalt)
    call field_get_val_s(ivarfl(ivart), cvar_scalt)
    call field_get_val_s(icrom, cpro_rom)
    call field_get_val_v(ivarfl(iu), cvar_vel)
    call field_get_val_s(ivarfl(ipr), cvar_p)

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    
    call field_get_val_s_by_name ('boundary_temperature',  twall)
    
    allocate(rod_center_xyz(5,5,2), rod_radius(5,5))

    rod_radius = fuel_rod_diameter/2.d0

    rod_center_xyz(1,1,1) = x11
    rod_center_xyz(1,1,2) = y11

    do j=1,5
      if(j>1)rod_center_xyz(1,j,1)=rod_center_xyz(1,1,1)+(j-1)*pitch
      do i=1,5
        rod_center_xyz(i,j,1) = rod_center_xyz(1,j,1)
        rod_center_xyz(i,j,2) = rod_center_xyz(1,1,2)-(i-1)*pitch
      end do
    end do

    !! Get all wall faces
    allocate(lstelt(nfabor))
    chaine = all_rod_walls
    call getfbr(chaine, nlelt, lstelt)

    allocate(xt1(npoint))
    allocate(xt2(npoint))
    allocate(xt_bulk(npoint))
    allocate(xw(npoint))
    allocate(xp(npoint))
    allocate(n_faces(npoint))

    do ii = 1, npoint
      sum_vol = 0.d0
      sum_vol_rho_w = 0.d0
      sum_vol_t = 0.d0
      sum_vol_t_rho_w = 0.d0
      sum_vol_w = 0.d0
      sum_vol_p = 0.d0
      n_faces(ii) = 0
      do iel=1, ncel
        temp = abs(xyzcen(3,iel)-xabs(ii))
        if(temp<1.e-5)then
          sum_vol = sum_vol+volume(iel)
          sum_vol_rho_w = sum_vol_rho_w+volume(iel)*cpro_rom(iel) &
                        * cvar_vel(3,iel)
          sum_vol_t = sum_vol_t+volume(iel)*cvar_scalt(iel)
          sum_vol_t_rho_w = sum_vol_t_rho_w+volume(iel)*cvar_scalt(iel)&
                          * cpro_rom(iel)*cvar_vel(3,iel)
          sum_vol_w = sum_vol_w+volume(iel)*cvar_vel(3,iel)
          sum_vol_p = sum_vol_p+volume(iel)*cvar_p(iel)
          n_faces(ii) = n_faces(ii)+1
        end if
      end do
      if (irangp.ge.0) then
        call parsom(sum_vol)
        call parsom(sum_vol_t)
        call parsom(sum_vol_w)
        call parsom(sum_vol_p)
        call parsom(sum_vol_rho_w)
        call parsom(sum_vol_t_rho_w)
        call parcpt(n_faces(ii))
      end if
      xt1(ii) = sum_vol_t/sum_vol
      xt_bulk(ii) = sum_vol_t_rho_w/sum_vol_rho_w
      xw(ii) = sum_vol_w/sum_vol
      xp(ii) = sum_vol_p/sum_vol
    end do

    if (ntcabs.eq.ntmabs) then
      impout = impusr(1)
      if (irangp.le.0) then
        write(cfilename,'(a20,a4)')'Axial_avg_quantities', '.dat'
        open(impout,file=Trim(cfilename))
        write(impout,*)  &
           'Axial_coordinate(m)  t_avg(C)  t_bulk(C)  w_bulk(m/s)  p_avg(Pa)&
           &    N_faces'
        do ii=1, npoint        
          write(impout,'(5Es15.5, I8)') xabs(ii), xt1(ii), xt_bulk(ii), &
                                        xw(ii), xp(ii), n_faces(ii)
        end do
        close(impout)
      endif
    endif

!    do k=1,5
!    do j=1,5

!    only calculate for the central rod
    do k=3,3
    do j=3,3

    do ii = 1, npoint
      sum_surf = 0.d0
      sum_surf_t = 0.d0
      n_faces(ii) = 0
      do ilelt = 1, nlelt
        ifac = lstelt(ilelt)
        iel = ifabor(ifac)
        jj = subchannel_ID(iel)
        cell_thichness = cell_height(jj)
        x = rod_center_xyz(j,k,1)
        y = rod_center_xyz(j,k,2)
        radius = sqrt((cdgfbo(1,ifac)-x)**2+(cdgfbo(2,ifac)-y)**2)
      
        if(radius < rod_radius(j,k)*(1.d0+tol2))then
          temp = abs(cdgfbo(3,ifac)-xabs(ii))
          if(temp<1.e-5)then
            sum_surf = sum_surf+surfbn(ifac)
            sum_surf_t = sum_surf_t+surfbn(ifac)*twall(ifac)
            n_faces(ii) = n_faces(ii)+1
          end if
        end if
      enddo
      if (irangp.ge.0) then
        call parsom(sum_surf)
        call parsom(sum_surf_t)
        call parcpt(n_faces(ii))
      end if
      xt2(ii) = sum_surf_t/sum_surf
    enddo

    if (ntcabs.eq.ntmabs) then
      impout = impusr(1)
      if (irangp.le.0) then
        write(cfilename,'(a16,2i02.2,a4)')'Axial_avg_twall_',j,k, '.dat'
        open(impout,file=Trim(cfilename))
        write(impout,*)  &
           'Axial_coordinate(m)  t_wall(C)  N_faces'
        do ii=1, npoint
          write(impout,'(2f10.5, I4)') xabs(ii), xt2(ii), n_faces(ii)
        end do
        close(impout)
      endif
    endif

    end do
    end do

    deallocate(rod_center_xyz, rod_radius)
    deallocate(lstelt)
    deallocate(xabs)
    deallocate(xt1)
    deallocate(xt2)
    deallocate(xt_bulk)
    deallocate(xw)
    deallocate(xp)
    deallocate(n_faces)

  end subroutine axial_avg_quantities


!===============================================================================
!! This routine is used to calculate the local Nusselt number of the rod walls
!! It is taken from post_boundary_nusselt in post_util.f90
  subroutine boundary_nusselt (nfbrps , lstfbr, bnussl)

    use paramx
    use pointe
    use entsor
    use cstnum
    use cstphy
    use optcal
    use numvar
    use parall
    use period
    use mesh
    use field
    use field_operator
    use cs_c_bindings

    implicit none

    ! Arguments

    integer, intent(in)                                :: nfbrps
    integer, dimension(nfbrps), intent(in)             :: lstfbr
    double precision, dimension(nfbrps), intent(out)   :: bnussl

    ! Local variables

    integer ::         inc, iprev
    integer ::         iel, ifac, iloc, ivar
    integer ::         ifcvsl, itplus, itstar

    double precision :: xvsl  , srfbn , heq
    double precision :: diipbx, diipby, diipbz
    double precision :: numer, denom

    double precision, allocatable, dimension(:) :: theipb
    double precision, allocatable, dimension(:,:) :: grad
    double precision, dimension(:), pointer :: cviscl
    double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
    double precision, dimension(:), pointer :: tplusp, tstarp
    double precision, dimension(:), pointer :: tscalp, hextp, hintp, dist_theipb
    character(len=50) :: chaine

    type(var_cal_opt) :: vcopt

    logical(c_bool), dimension(:), pointer ::  cpl_faces

    call field_get_id_try('tplus', itplus)
    call field_get_id_try('tstar', itstar)

    if (itstar.ge.0 .and. itplus.ge.0) then

      ivar = isca(iscalt)

      call field_get_val_prev_s(ivarfl(ivar), tscalp)

      call field_get_val_s(itplus, tplusp)
      call field_get_val_s(itstar, tstarp)

      ! Boundary condition pointers for diffusion

      call field_get_coefaf_s(ivarfl(ivar), cofafp)
      call field_get_coefbf_s(ivarfl(ivar), cofbfp)

      ! Boundary condition pointers for diffusion with coupling

      call field_get_hext(ivarfl(ivar), hextp)
      call field_get_hint(ivarfl(ivar), hintp)

      ! Compute variable values at boundary faces

      call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

      allocate(theipb(nfabor))
      theipb = 0.d0
  
      do iloc = 1, nfbrps
        ifac = lstfbr(iloc)
        iel = ifabor(ifac)
        theipb(ifac) = tscalp(iel)
      enddo

      ! Reconstructed fluxes
      if (vcopt%ircflu .gt. 0 .and. itbrrb.eq.1) then

        ! Compute gradient of temperature / enthalpy

        allocate(grad(3,ncelet))

        inc = 1
        iprev  = 1 ! previsous time step value

        call field_gradient_scalar &
          (ivarfl(ivar), iprev, inc, grad)

        ! Compute reconstructed temperature

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          iel = ifabor(ifac)

          diipbx = diipb(1,ifac)
          diipby = diipb(2,ifac)
          diipbz = diipb(3,ifac)

          theipb(ifac) = theipb(ifac)+diipbx*grad(1,iel)+diipby*grad(2,iel) &
                       + diipbz*grad(3,iel)
        enddo

        deallocate(grad)
      endif

      if (vcopt%icoupl.gt.0) then
        call field_get_coupled_faces(ivarfl(ivar), cpl_faces)
        allocate(dist_theipb(nfabor))
        call cs_ic_field_dist_data_by_face_id(ivarfl(ivar), 1, theipb, dist_theipb)
      endif

      ! Physical property pointers

      call field_get_key_int (ivarfl(ivar), kivisl, ifcvsl)
      if (ifcvsl .ge. 0) then
        call field_get_val_s(ifcvsl, cviscl)
      endif

      ! Boundary condition pointers for gradients and advection

      call field_get_coefa_s(ivarfl(ivar), coefap)
      call field_get_coefb_s(ivarfl(ivar), coefbp)

      ! Compute using reconstructed temperature value in boundary cells

      do iloc = 1, nfbrps
        ifac = lstfbr(iloc)
        iel = ifabor(ifac)

        if (ifcvsl.ge.0) then
          xvsl = cviscl(iel)
        else
          xvsl = visls0(iscalt)
        endif

        srfbn = surfbn(ifac)

        numer = (cofafp(ifac) + cofbfp(ifac)*theipb(ifac)) * distb(ifac)
        ! here numer = 0 if current face is coupled

        if (vcopt%icoupl.gt.0) then
          if (cpl_faces(ifac)) then
            heq = hextp(ifac) * hintp(ifac) / ((hextp(ifac) + hintp(ifac))*srfbn)
            numer = heq*(theipb(ifac)-dist_theipb(ifac)) * distb(ifac)
          endif
        endif

        denom = xvsl * tplusp(ifac)*tstarp(ifac)

        if (abs(denom).gt.1e-30) then
          bnussl(ifac) = numer / denom
        else
          bnussl(ifac) = 0.d0
        endif

      enddo

      if (vcopt%icoupl.gt.0) then
        deallocate(dist_theipb)
      endif

      deallocate(theipb)

    else ! default if not computable

      do iloc = 1, nfbrps
        bnussl(iloc) = -1.d0
      enddo

    endif

    return
  end subroutine boundary_nusselt

end module coupling_resolved

!-------------------------------------------------------------------------------

!> (DOXYGEN_SHOULD_SKIP_THIS) \endcond
