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

module SubChCFD

implicit none

!=============================================================================== 
!! Parameter specification
!!------------------------------------------------------------------------------
!! Bundle type
! PWR: 1, AGR: 2, prismatic HTGR (General Atomics design): 3
  integer,parameter :: bundle_type = 1
!! Subchannel mesh name
! PWR: 'PWR_subchannel_mesh.dat'
! AGR: 'AGR_subchannel_mesh.dat'
! prismatic HTGR (General Atomics design): 'GA_HTGR_subchannel_mesh.dat'
  character(len=80),parameter :: subchannel_mesh_name = 'PWR_subchannel_mesh.dat'
!! Boundary name
!        Pin wall |      Housing wall        | Guide tube wall | Inlet | Outlet
! PWR:  rod_walls | side_walls/housing_walls |    N/A          | inlet | outlet
! AGR:   RodWall  |       SleeveWall         |  GuideTubeWall  | Inlet | Outlet
! HTGR:channel_wall|            N/A           |    N/A          | Inlet | Outlet
  character(len=80),parameter :: pin =     'rod_walls',                        &
                                 pin_un =  'N/A',                              &
                                 housing = 'housing_walls',                    &
                                 guide   = 'N/A',                              &
                                 inlet   = 'inlet',                            &
                                 outlet  = 'outlet'

!! Geometrical information
!!------------------------------------------------------------------------------
!! For AGR rod bundles and HTGR (GA design) fuel assemblies
!! Sepcific for AGR rod bundles
! voidage of the fuel pins
  double precision, parameter :: voidage = 0.729558d0
!! Subchannel parameters
! bundle sleeve radius (m)
  double precision,parameter :: bundle_R = 0.95d0
! coordinates of origin for the circular bundle
  double precision,parameter :: x_origin = 0.d0, &
                                y_origin = 0.d0
!!For both
! maximum number of nodes in subchannel mesh elements
  integer,parameter :: nCol_max = 8

!!------------------------------------------------------------------------------
!! Specific for PWR bundles
! Parameters for identification of the subchannel type
  double precision, parameter :: spacing_edge   = 0.95d0
  double precision, parameter :: spacing_corner = 0.45d0
! Rod diameter
  double precision, parameter :: rod_diameter   = 0.0095d0 !(m)
! Pitch-to-diameter ratio for interior subchannels
  double precision, parameter :: PbyD_interior  = 1.326d0
! Pitch-to-diameter ratio for the edge subchannels
  double precision, parameter :: PbyD_edge      = 1.326d0
! Pitch-to-diameter ratio for the corner subchannels
  double precision, parameter :: PbyD_corner    = 1.326d0

!!------------------------------------------------------------------------------
!! Slip wall parameters
!! Method to decouple the first cell and the wall for velocity components
!  Set to symmety for all variables: 1 (useful for 0-equation turbulence models)
!  Explicit method for slip walls:   2 (useful for regular geometries)
!  Pseudo symmetry for velocities:   3 (useful for advanced turbulence models 
!                                       and wall heat transfer)
  integer, parameter :: slip_type = 3

!!------------------------------------------------------------------------------
!! Mesh parameters  
! Scale the meshes if generated in other units than 'meters'
  logical,parameter :: subchannel_mesh_scale = .true.,                         &
                       mesh_scale            = .true.
  double precision, parameter :: subchannel_scaling_factor = 1.d-3,            &
                                 mesh_scaling_factor       = 1.d-3
! Use user defined routines for non-ideal rod bundles
  logical, parameter :: user_subchannel_mesh   = .false.,                      &
                        user_frictional_factor = .false.,                      &
                        user_nusselt_number    = .false.

!!------------------------------------------------------------------------------
!! Physical parameters
! Thermal boundary type and value
! 1: specify the wall heat flux 
! 2: specify the wall temperature (still under test)
  integer, parameter :: thermal_type_rw  = 1,                                  &
                        thermal_type_rwu = 1,                                  &
                        thermal_type_hw  = 1,                                  &
                        thermal_type_o   = 1
  double precision, parameter :: value_rw  = 0.d0,                             &
                                 value_rwu = 0.d0,                             &
                                 value_hw  = 0.d0,                             &
                                 value_o   = 0.d0
! Heating or cooling (for the use in Dittus-Boelter Nusselt number correlation)
  logical, parameter :: heating = .true.
!!Nusselt number correction (not recommend to use, set to 1: no correction)
  double precision, parameter :: crr_Nu = 1.d0 !0.82d0
!! Use data table to define variable physical properties
  logical, parameter :: physical_property_table = .false.
  character(len=50),parameter :: name_data_table = 'water_enthalpy_based.txt'
!! Increment of the thermal variable in the property table (should be consistent
!! with the data table created. Data table should follow the formats shown in 
!! the examples)
  double precision, parameter :: increment = 1000.d0

!!------------------------------------------------------------------------------
!! Coupling module
!!----------------
!! General settings
!! Specify if using the coupling module
  logical, parameter :: cpl = .true.
!! Allow momentum feedback from the embedded sub-model to the SubChCFD model
!! (two-way coupling)
  logical, parameter :: feedback = .true.
!! Use user defined feedback zone
  logical, parameter :: feedback_user = .false.
!! Set convergence factor for momentum in two-way coupling
  double precision :: lambda = 0.1d0
!! Allow energy feedback from the embedded sub-model to the SubChCFD model
!! (two-way coupling)  
  logical, parameter :: feedback_energy = .false.
!! Set convergence factor for energy in two-way coupling
  double precision :: lambda_energy= 0.1d0
!! Use semi-implicit correction source term for momentum
  logical, parameter :: imp_source = .true.
!!----------------  
!! Source term correction type (integration method to ensure overall head losse
!!                              correct)
!! 0: no source term correction (default)
!! 1: constant CX covering the entire domain
!! 2: constant CX covering the coupling region
!! 3: varing CX (currently not available, computational expensive)
  integer, parameter :: mom_crr = 0
!!----------------  
!! Set pressure reconstruction method
!! 1: using mom_hat 
!! 2: using mom_pt (based on dp from the embedded sub-model)
  integer, parameter :: ipre_rec = 2
!! Set pressure reconstruction frequency
  integer, parameter :: fpre_rec = 100
!! Receive dp from the remote instance (embedded) for pressure reconstruction 
  logical, parameter :: receive_dp = .true.
!! frequency of receiving dp (Note: freceive_dp <= fpre_rec)
  integer, parameter :: freceive_dp = 100
!! Set reference pressure point
  double precision, dimension(3) :: xyz_ref = 0.d0  !(default)

!! interpolation scheme from the coarse mesh to the fine mesh for boundary 
!! condition definition of the latter (useful for velocity components and 
!! scalers, pressure and turbulence quantities are not included)
!! 0: 0th order projection
!! 1: 1st order projection
  integer, parameter :: interpol = 1
!!------------------------------------------------------------------------------
!! Porous module
!!----------------
!! Specify if using the porous module
  logical, parameter :: prs = .false.
!! Specify isotropic or anisotropic porosity model
!! 1: isotropic (only use volume porosity)
!! 2: anisotropic (use both volume and face porosity)
!! 3: porous model with improved numerics
  integer, parameter :: iprs = 1
!! Auxiliary mesh names
  character(len=30),parameter :: solid_mesh_volume  = 'solid_volume.dat'
  character(len=30),parameter :: solid_mesh_surface = 'solid_surface.dat'
  character(len=30),parameter :: filtering_mesh     = 'filtering_mesh.dat'
!! Auxiliary mesh scaling
  logical, parameter :: solid_mesh_scale_v = .false.,                           &
                        solid_mesh_scale_f = .false.
  double precision, parameter :: solid_scaling_factor_v = 1.d-3,               &
                                 solid_scaling_factor_f = 1.d-3
  logical,parameter :: filtering_mesh_scale = .false.
  double precision, parameter :: filtering_scaling_factor = 1.d-3
!! Porosity calculation method
!! sub-channel based:    1 (embedding method)
!! computing cell based: 2 (embedding method)
!! filtering mesh based: 3 (interfacing method)
  integer, parameter :: iporo_cal = 1
!! Porosity data file
  character(len=30),parameter :: v_porous = "volume_porosity.dat",             &
                                 f_porous = "face_porosity.dat"

!!------------------------------------------------------------------------------
!! Simulation statistics
!! wall clock time
  double precision :: wall_clock_time = 0.d0
  double precision :: calculation_time = 0.d0

!!------------------------------------------------------------------------------
!! Other parameters
!! Math constants  
  double precision, parameter :: ee = 2.71828183d0,                            &
                                 pi = 3.14159265d0

!===============================================================================  
!! Subchannel variables  
  integer :: n_cell
  integer, allocatable, dimension(:) :: cell_type
  double precision, allocatable, dimension(:) :: cell_height, cell_vol
  double precision, allocatable, dimension(:) :: x_channel,y_channel,z_channel,&
                                     x_channel_min,y_channel_min,z_channel_min,&
                                     x_channel_max,y_channel_max,z_channel_max
  double precision, allocatable, dimension(:) :: subchannel_surf,     &
                                                 subchannel_vol,      &
                                                 subchannel_rodsurf,  &
                                                 subchannel_wallsurf, &
                                                 PbyD
  double precision, allocatable, dimension(:) :: store1
  character(len=80) :: all_rod_walls,      &
                       heated_rod_walls,   &
                       unheated_walls,     &
                       all_walls 

  double precision, allocatable, dimension (:,:) :: user_prop
  integer :: n_rows

  double precision, allocatable, dimension(:,:) :: filter_cell_vol

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
!! This routine is used to pass logical variable cpl to .c routines
  function coupling() bind (C,name='coupling')

    use, intrinsic :: iso_c_binding

    implicit none

    logical(kind=c_bool) :: coupling

    coupling = cpl

  end function coupling


!===============================================================================
!! This routine is used to pass logical variable prs to .c routines
  function porosity() bind (C,name='porosity')

    use, intrinsic :: iso_c_binding

    implicit none

    logical(kind=c_bool) :: porosity

    porosity = prs

  end function porosity


!===============================================================================
!! This routine is used to pass integer variable iprs to .c routines
  function iporosity() bind (C,name='iporosity')

    use, intrinsic :: iso_c_binding

    implicit none

    integer(kind=c_int) :: iporosity

    iporosity = iprs

  end function iporosity


!===============================================================================
!! This routine is used to pass integer variable iporo_cal to .c routines
  function iporosity_cal() bind (C,name='iporosity_cal')

    use, intrinsic :: iso_c_binding

    implicit none

    integer(kind=c_int) :: iporosity_cal
 
    iporosity_cal = iporo_cal

  end function iporosity_cal


!===============================================================================
!! This routine is used to calculate porosities
  subroutine cal_porosity() bind (C,name='cal_porosity')

    use, intrinsic :: iso_c_binding

    implicit none

    call subchannel_mesh()
    call porosities()

  end subroutine cal_porosity


!===============================================================================
!! This routine should be called after subchannel_mesh to allocate memories 
!! mainly for subchannel variables 
  subroutine allocate_var()

    implicit none

!! Arrays for subchannel wall surface area and volume (excluding solid)
    allocate(subchannel_surf(n_cell), subchannel_vol(n_cell), &
             subchannel_rodsurf(n_cell), subchannel_wallsurf(n_cell))
             
!! Storage arrays for subchannel variables
    allocate(store1(n_cell))

  end subroutine allocate_var


!===============================================================================
!! This routine should be called to dump the memory of common variables
!! after the whole of simulation is completed
  subroutine deallocate_var()

    implicit none

    deallocate(cell_height, cell_vol, cell_type)
    deallocate(subchannel_surf, subchannel_vol, subchannel_rodsurf, &
               subchannel_wallsurf, PbyD)
    deallocate(store1)
  
  end subroutine deallocate_var

!===============================================================================
!! This routine is used to bridge the subchannel mesh and computing mesh
  subroutine subchannel_mesh()

    implicit none

    select case (bundle_type)
    case (1)
      if(trim(pin)=='N/A')then
        all_rod_walls = trim(pin_un)
      elseif(trim(pin_un)=='N/A')then
        all_rod_walls = trim(pin)
      else
        all_rod_walls = trim(pin)//' or '//trim(pin_un)
      end if
      heated_rod_walls = trim(pin)
      if (trim(pin_un)=='N/A')then
        unheated_walls = trim(housing)
      else
        unheated_walls = trim(housing)//' or '//trim(pin_un)
      end if
      all_walls = trim(all_rod_walls)//' or '//trim(housing)
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
     case(3)
      all_rod_walls = trim(pin)
      heated_rod_walls = trim(pin)
      unheated_walls = trim(housing)
      all_walls = trim(pin)//' or '//trim(housing)
      if (.not. user_subchannel_mesh) then
        call subchannel_mesh_HTGR_GA()
      else
        call subchannel_mesh_user()
      end if
    end select

    call allocate_var()

    call hydraulic_diameter()

  end subroutine subchannel_mesh


!===============================================================================
!! This routine is to acquire the data dimensions of the subchannel mesh 
  subroutine mesh_dimensions(mesh_file, n_node, n_entity)

    use parall

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
      if(irangp.le.0)&
        write(*,*) trim(mesh_file)//" doesn't exist. Program terminated."
      stop
    end if

!!  Read in number of nodes and 1D, 2D and 3D elements 
    read(100,*)n_node, n_entity

    close(100)

  end subroutine mesh_dimensions


!===============================================================================  
!! This routine is to read-in the data of the subchannel mesh
  subroutine subchannel_mesh_readin_pwr(mesh_file, n_node, n_entity, &
                                        x, y, z, index_cell, n_cell)
    use parall

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
      if(irangp.le.0)&
        write(*,*) trim(mesh_file)//" doesn't exist. Program terminated."
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
    double precision, allocatable, dimension(:,:) :: xx, yy, zz
    integer, allocatable, dimension(:,:) :: index_cell
    integer :: i, iel
    double precision :: SubchannelID_err
    double precision, dimension(:), pointer :: subchannel_ID, &
                                               subchannel_type, &
                                               pitch2diameter

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)

    call mesh_dimensions(subchannel_mesh_name, n_node, n_entity)

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

!! Calculate minimum and maximum coordinates value of the subchannel cell     
    allocate(xx(n_cell, 2), yy(n_cell, 2), zz(n_cell, 2))

!! Module variable, should be deallocated after the whole simulation
    allocate(cell_height(n_cell))
    allocate(cell_vol(n_cell))
    allocate(cell_type(n_cell))
    allocate(PbyD(n_cell))
    allocate(x_channel(n_cell),y_channel(n_cell),z_channel(n_cell))
    allocate(x_channel_min(n_cell),y_channel_min(n_cell),z_channel_min(n_cell))
    allocate(x_channel_max(n_cell),y_channel_max(n_cell),z_channel_max(n_cell))

    do i=1, n_cell
      xx(i,1) = minval(x(index_cell(i,:)))
      xx(i,2) = maxval(x(index_cell(i,:)))
      yy(i,1) = minval(y(index_cell(i,:)))
      yy(i,2) = maxval(y(index_cell(i,:)))
      zz(i,1) = minval(z(index_cell(i,:)))
      zz(i,2) = maxval(z(index_cell(i,:)))
!!  Stored for global use
      x_channel_min(i) = xx(i,1)
      y_channel_min(i) = yy(i,1)
      z_channel_min(i) = zz(i,1)
      x_channel_max(i) = xx(i,2)
      y_channel_max(i) = yy(i,2)
      z_channel_max(i) = zz(i,2)

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
    subchannel_ID = -1.d0
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
      if(subchannel_ID(iel)>0.d0)then
        pitch2diameter(iel) = PbyD(int(subchannel_ID(iel)))
        subchannel_type(iel) = cell_type(int(subchannel_ID(iel)))
      end if
    end do

!! Output for checking results
!    do i=1,n_cell
!      write(*,*)i, xx(i,1),xx(i,2),yy(i,1),yy(i,2),zz(i,1),zz(i,2)
!      write(*,*)i, x(i), y(i), z(i)
!      write(*,*) i, cell_height(i)
!      write(*,*) i, cell_vol(i), cell_type(i), PbyD(i)
!    end do

!    open(200,file='test1.txt')
!    do iel=1, ncel
!      write(200,*)iel, subchannel_ID(iel)
!    end do
!    close(200)

    deallocate(x,y,z)
    deallocate(index_cell)
    deallocate(xx,yy,zz)

  end subroutine subchannel_mesh_PWR


!===============================================================================  
!! This routine is to read in the data of the subchannel mesh
  subroutine subchannel_mesh_readin_agr(mesh_file, n_node, n_entity, &
                                        x, y, z, index_cell, n_cell)

    use parall

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
      if(irangp.le.0)&
        write(*,*) trim(mesh_file)//" doesn't exist. Program terminated."
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

    call mesh_dimensions(subchannel_mesh_name, n_node, n_entity)
    
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
    allocate(x_channel_min(n_cell),y_channel_min(n_cell),z_channel_min(n_cell))
    allocate(x_channel_max(n_cell),y_channel_max(n_cell),z_channel_max(n_cell))
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

!! Calculate cell center coordinates of the subchannel cells (to be added)

!! Calculate maximum and minimum xyz of the subchannel cells (to be added)

!! Find the subchannel indext for each control volume
    subchannel_ID = -1.d0
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
      if(subchannel_ID(iel)>0.d0)then
        pitch2diameter(iel) = PbyD(int(subchannel_ID(iel)))
        subchannel_type(iel) = cell_type(int(subchannel_ID(iel)))
      end if
    end do

!! Output for checking results
!    do i=1,n_cell
!      write(*,*)i, x_face1(1,i),y_face1(1,i),z_face1(1,i)
!      write(*,*)i, x_face1(2,i),y_face1(2,i),z_face1(2,i)
!      write(*,*)i, x_face1(3,i),y_face1(3,i),z_face1(3,i)
!      write(*,*)i, x_face2(1,i),y_face2(1,i),z_face2(1,i)
!      write(*,*)i, x_face2(2,i),y_face2(2,i),z_face2(2,i)
!      write(*,*)i, x_face2(3,i),y_face2(3,i),z_face2(3,i)
!    end do

!    open(200,file='test1.txt')
!    do iel=1, ncel 
!      write(200,*)iel, subchannel_ID(iel)
!    end do
!    close(200)

    deallocate(x,y,z)
    deallocate(index_cell)
    deallocate(xx,yy,zz)
    deallocate(x_face1,y_face1,z_face1,x_face2,y_face2,z_face2)

  end subroutine subchannel_mesh_AGR


!===============================================================================
!! This routine is for HTGR fuel (GA design) assemblies
  subroutine subchannel_mesh_HTGR_GA()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: n_node, n_entity
    double precision, allocatable, dimension(:) ::  x, y, z
    double precision, allocatable, dimension(:,:) :: xx, yy, zz
    double precision, allocatable, dimension(:,:) :: x_faceo1, y_faceo1, &
                                                     z_faceo1, x_faceo2, &
                                                     y_faceo2, z_faceo2, &
                                                     x_face1,  y_face1,  &
                                                     z_face1,  x_face2,  &
                                                     y_face2,  z_face2
    integer, allocatable, dimension(:,:) :: index_cello, index_cell
    integer :: i, j, k1, k2, iel, ilelt, nlelt, n_cello
    integer, allocatable, dimension(:) :: lstelt
    double precision, dimension(:), pointer :: subchannel_ID, &
                                               subchannel_type
    double precision, dimension(:), pointer :: cht
    character*80 :: chaine
    logical :: within, alive
    integer :: indicator

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('SubchannelType', subchannel_type)
    call field_get_val_s_by_name('CHTInternalCoupling', cht)
    cht = 1.d0

 !!  Check if subchannel ID data file exists
    inquire(file='subchannelID.dat', exist=alive)
    if(alive) then
      call io_subchannelID(1)
      return
    end if

!! Dealing with the fluid region only
    allocate(lstelt(ncel))
    chaine = 'channel'
    call getcel(chaine, nlelt, lstelt)

    call mesh_dimensions(subchannel_mesh_name, n_node, n_entity)

    allocate(x(n_node),y(n_node),z(n_node))
    allocate(index_cello(n_entity, nCol_max))
    allocate(index_cell(n_entity, nCol_max))
    x = 0.
    y = 0.
    z = 0.
    index_cello = -1
    index_cell = -1

    call subchannel_mesh_readin_agr(subchannel_mesh_name, n_node, n_entity, &
                                    x, y, z, index_cello, n_cello)

!! Scale mesh from mm to m
    if(subchannel_mesh_scale)then
      x = x*subchannel_scaling_factor
      y = y*subchannel_scaling_factor
      z = z*subchannel_scaling_factor
    end if

!! Calculate minimum and maximum coordinates value of the subchannel cell
    allocate(xx(n_cello, nCol_max), &
             yy(n_cello, nCol_max), &
             zz(n_cello, nCol_max))
    allocate(x_faceo1(nCol_max/2,n_cello), &
             y_faceo1(nCol_max/2,n_cello), &
             z_faceo1(nCol_max/2,n_cello), &
             x_faceo2(nCol_max/2,n_cello), &
             y_faceo2(nCol_max/2,n_cello), &
             z_faceo2(nCol_max/2,n_cello))
    allocate(x_face1(nCol_max/2,n_cello), &
             y_face1(nCol_max/2,n_cello), &
             z_face1(nCol_max/2,n_cello), &
             x_face2(nCol_max/2,n_cello), &
             y_face2(nCol_max/2,n_cello), &
             z_face2(nCol_max/2,n_cello))
    xx = 0.d0
    yy = 0.d0
    zz = 0.d0
    x_faceo1 = 0.d0
    x_faceo2 = 0.d0
    y_faceo1 = 0.d0
    y_faceo2 = 0.d0
    z_faceo1 = 0.d0
    z_faceo2 = 0.d0
    x_face1 = 0.d0
    x_face2 = 0.d0
    y_face1 = 0.d0
    y_face2 = 0.d0
    z_face1 = 0.d0
    z_face2 = 0.d0

    do i=1, n_cello
      xx(i,:) = x(index_cello(i,:))
      yy(i,:) = y(index_cello(i,:))
      zz(i,:) = z(index_cello(i,:))

      k1 = 1
      k2 = 1
      x_faceo1(k1,i) = xx(i,1)
      y_faceo1(k1,i) = yy(i,1)
      z_faceo1(k1,i) = zz(i,1)
      do j=2,6
        if(abs(zz(i,1)-zz(i,j))>1.e-10)then
          x_faceo2(k2,i) = xx(i,j)
          y_faceo2(k2,i) = yy(i,j)
          z_faceo2(k2,i) = zz(i,j)
          k2 = k2+1
        else
          k1 = k1+1
          x_faceo1(k1,i) = xx(i,j)
          y_faceo1(k1,i) = yy(i,j)
          z_faceo1(k1,i) = zz(i,j)
        end if
      end do
    end do

!! remove the subchannel cells that are not corresponding to fluid
    n_cell = 0
    do i=1, n_cello
      within = .false.
      indicator = 0
      do ilelt = 1, nlelt
        iel = lstelt(ilelt)
        call find_suchannel_index(xyzcen(1,iel), xyzcen(2,iel), x_faceo1(:,i),&
                                  y_faceo1(:,i), index_cello(i,:), within)
        if((xyzcen(3,iel)-z_faceo1(1,i))*(xyzcen(3,iel)-z_faceo2(1,i))<0.&
        .and.within) then
          indicator = 1
          exit
        end if
      end do
      if(irangp.ge.0) call parcmx(indicator)
      if(indicator==1)then
        n_cell = n_cell+1
        index_cell(n_cell,:) = index_cello(i,:)
        x_face1(:,n_cell) = x_faceo1(:,i)
        x_face2(:,n_cell) = x_faceo2(:,i)
        y_face1(:,n_cell) = y_faceo1(:,i)
        y_face2(:,n_cell) = y_faceo2(:,i)
        z_face1(:,n_cell) = z_faceo1(:,i)
        z_face2(:,n_cell) = z_faceo2(:,i)
      end if
    end do

!! Module variable, should be deallocated after the whole simulation
    allocate(cell_height(n_cell))
    allocate(x_channel(n_cell),y_channel(n_cell),z_channel(n_cell))
    allocate(x_channel_min(n_cell),y_channel_min(n_cell),z_channel_min(n_cell))
    allocate(x_channel_max(n_cell),y_channel_max(n_cell),z_channel_max(n_cell))

    do i=1, n_cell
      cell_height(i) = abs(z_face1(1,i)-z_face2(1,i))
    end do

!! Find the subchannel index for each control volume
    subchannel_ID = -1.d0
    subchannel_type = -1.d0
    do ilelt = 1, nlelt
      iel = lstelt(ilelt)
      subchannel_type(iel) = 1.d0
      within = .false.
      do i=1, n_cell
        call find_suchannel_index(xyzcen(1,iel), xyzcen(2,iel), x_face1(:,i),&
                                  y_face1(:,i), index_cell(i,:), within)
        if((xyzcen(3,iel)-z_face1(1,i))*(xyzcen(3,iel)-z_face2(1,i))<0.&
          .and.within)&
          subchannel_ID(iel) = i
      end do

!      do i=1, n_cello
!        call find_suchannel_index(xyzcen(1,iel), xyzcen(2,iel), x_faceo1(:,i),&
!                                  y_faceo1(:,i), index_cello(i,:), within)
!        if((xyzcen(3,iel)-z_faceo1(1,i))*(xyzcen(3,iel)-z_faceo2(1,i))<0.&
!          .and.within)&
!          subchannel_ID(iel) = i
!      end do

    end do

    call io_subchannelID(2)

    deallocate(x,y,z)
    deallocate(index_cell, index_cello)
    deallocate(xx,yy,zz)
    deallocate(x_faceo1,y_faceo1,z_faceo1,x_faceo2,y_faceo2,z_faceo2)
    deallocate(x_face1,y_face1,z_face1,x_face2,y_face2,z_face2)

  end subroutine subchannel_mesh_HTGR_GA


!===============================================================================
!! This routine is to read or write subchannel IDs
  subroutine io_subchannelID(mode)

    use numvar
    use field
    use mesh
    use parall

    implicit none

    interface
      subroutine global_cell_id(c_id, c_g_id)        &
        bind(C, name='global_cell_id')
        use, intrinsic :: iso_c_binding
        implicit none
        integer(kind=c_int) :: c_id
        integer(kind=c_int) :: c_g_id
      end subroutine global_cell_id
    end interface

    integer :: mode

    double precision, pointer, dimension(:) :: subchannel_ID
    integer :: i, ii, iel
    double precision, allocatable, dimension(:) :: subchannel_IDg
    integer :: c_g_id
    integer :: dummy
    integer :: aux, bcast

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    allocate(subchannel_IDg(ncelgb))

    if(mode==1) then
      open(100, file='subchannelID.dat', status='old') 
!!  read in subchannel IDs
      read(100, *)n_cell

      do i=1, ncelgb
        read(100, *) dummy, subchannel_IDg(i)
      end do
      close(100)

      do iel=1, ncel
        if (irangp.ge.0) then
          call global_cell_id(iel, c_g_id)
          subchannel_ID(iel) = subchannel_IDg(c_g_id)
        else
          subchannel_ID(iel) = subchannel_IDg(iel)
        end if
      end do

    else if(mode==2) then

      subchannel_IDg = 1.d12

      do iel=1, ncel
        if (irangp.ge.0) then
          call global_cell_id(iel, c_g_id)
          subchannel_IDg(c_g_id) = subchannel_ID(iel)
        else
          subchannel_IDg(iel) = subchannel_ID(iel)
        endif
      end do

      if (irangp.ge.0) then
        do i=1, ncelgb
          aux = 1
          call parfpt (aux, bcast, subchannel_IDg(i))
          call parall_bcast_r(bcast, subchannel_IDg(i))
        end do
      end if

      if(irangp.le.0)then
        open (100, file='subchannelID.dat', status='replace')
        write(100, *) n_cell
        do i=1, ncelgb
          write(100,*) i, subchannel_IDg(i)
        end do
        close(100)
      end if

    end if

    deallocate(subchannel_IDg)

  end subroutine io_subchannelID

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

    if(bundle_type/=2) return

    radius1 = sqrt((x1-x_origin)**2+(y1-y_origin)**2)
    radius2 = sqrt((x2-x_origin)**2+(y2-y_origin)**2)

    if(abs(radius1-bundle_R)<1.e-10 .and. abs(radius2-bundle_R)<1.e-10)&
      edge_corr = .true.

  end function edge_corr


!===============================================================================
!! This routine is to read-in data of the solid mesh for porosity calculation  
  subroutine porosities()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    logical :: alive, match
    integer :: i, iel
    integer :: itemp
    character(len=30) :: ctemp
    double precision, dimension(n_cell) :: sub_poro_zone

    double precision,  dimension(:), pointer :: porosity_v
    double precision,  dimension(:), pointer :: porous_zone
    double precision,  dimension(:), pointer :: subchannel_ID

    integer :: mode

    integer :: aux, bcast
    double precision :: temp

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('VolumePorosity', porosity_v)
    call field_get_val_s_by_name('PorousZone', porous_zone)

    if(irangp<=0)then
      if(iporo_cal==1)then
        write(*,*)'Porosity calculation method: subchannel based'
      else if(iporo_cal==2)then
        write(*,*)'Porosity calculation method: cell based'
      else if(iporo_cal==3)then
        write(*,*)'Porosity calculation method: filtering mesh based'
      else
        write(*,*)'Error set for porosity calculation method!'
        write(*,*)'Simulation is terminated.'
        stop
      end if
    end if

!! Check if porosity data file exists and read in porosities from file if so
    inquire(file=v_porous, exist=alive)

    if(alive) then
      open(100, file=v_porous, status='old')
      read(100, '(A30, I1)') ctemp, itemp
      close(100)

      if (itemp/=iporo_cal) then
        match = .false.
      else
        match = .true.
      end if

      if(match)then

        if(iporo_cal==3) call porous_filtering_mesh()

        if(irangp<=0)write(*,*)'reading volume porosities ...'
        mode = 1
        call io_porosity(porosity_v, mode)
        if(irangp<=0)write(*,*)'reading completed.'

      end if

    end if

!! Calculate volume porosities if data file does not exist
    if((.not.alive) .or. (.not.match)) then
      if(irangp<=0) then
        write(*,*)'porosity data file not available, porosity calculation &
                   & is required'
      end if
      call volume_porosity()
      if(irangp<=0)write(*,*)'porosity calculation completed.'
    end if

!! Mark the porous zone based on subchannel divisions (i.e. any region in a
!! subchannel has non unity porosity, the whole subchannel is marked as 
!! porous zone)
    sub_poro_zone = 0.d0
    porous_zone = 0.d0
    do iel=1, ncel
      if((1.d0-porosity_v(iel))>1.d-12.and.subchannel_ID(iel)>0.d0)then
        sub_poro_zone(int(subchannel_ID(iel)))=1.d0
      end if
    end do

!! In case subchannels sit across MPI ranks
    if (irangp.ge.0) then
      do i=1, n_cell
        aux = 1
        temp = -sub_poro_zone(i)
        call parfpt (aux, bcast, temp)
        call parall_bcast_r(bcast, sub_poro_zone(i))
      end do
    end if

    do iel=1, ncel
      if(sub_poro_zone(int(subchannel_ID(iel)))>0.d0 &
         .or.subchannel_ID(iel)<0.d0)then
         
        porous_zone(iel) = 1.d0
      end if
    end do

  end subroutine porosities


!===============================================================================
!! This routine is used to read auxiliary meshes with only one type of cells
  subroutine auxiliary_mesh_readin(mesh_file, n_node, n_entity, n_vertices,&
                               x, y, z, index_cell)

    use parall

    implicit none

!!  Aguments
    character(len=30), intent(IN) :: mesh_file
    integer, intent(IN) :: n_node, n_entity, n_vertices
    double precision, dimension(:), intent (OUT) ::  x, y, z
    integer, dimension(:,:), intent (OUT) :: index_cell

!!  Local variables
    character(1) :: dummy
    logical :: alive
    integer :: itemp, stat
    integer :: i

!!  Check if mesh file exists
    inquire(file=mesh_file, exist=alive)

    if(alive) then
      open(100, file=mesh_file, status='old')
    else
      if(irangp.le.0)&
        write(*,*) trim(mesh_file)//" doesn't exist. Program terminated."
      stop
    end if

    read(100,*)

!!  Read in node cordinates of the subchannel mesh
    do i=1,n_node
      read(100,*)dummy, x(i), y(i), z(i)
    end do

!!  Read in information of entities  
    do i=1, n_entity
      read(100,*, iostat=stat)dummy, itemp
      if(stat/=0)exit
      if(mod(itemp,100)==n_vertices)then
        backspace(100)
        read(100,*, iostat=stat)dummy, itemp, index_cell(i, 1:n_vertices)
      else
        if(irangp.le.0)&
          write(*,*) "Error reading "//trim(mesh_file)
        stop
      end if
    end do

    close(100)

  end subroutine auxiliary_mesh_readin


!===============================================================================
!! This routine is to calculate volume porosities
  subroutine volume_porosity()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: n_node, n_entity, n_vertices
    integer :: mode
    double precision, allocatable, dimension(:) ::  x, y, z
    integer, allocatable, dimension(:,:) :: index_cell

    integer :: i, ii, j, iel, iel1, iel2, ifac, n, n_fac, n_fbr
    double precision, allocatable, dimension(:,:) :: xx, yy, zz
    double precision, allocatable, dimension(:,:) :: cell_xyz
    double precision, allocatable, dimension(:) :: cell_vol

    double precision, dimension(3) :: vec_A, vec_B, vec_C, vec_AxB
    double precision, dimension(3) :: solid_min, solid_max
    double precision :: solid_vol

    integer, dimension(n_cell) :: index_sub
    double precision, dimension(n_cell) :: sub_poro

    integer, dimension(nfac) :: index_fac, flag_fac
    integer, dimension(nfabor) :: index_fbr
    double precision, dimension(ncel) :: solid_vol_cell
    double precision :: solid_vol_pyrmd

    integer :: n_node_f, n_entity_f, n_vertices_f
    double precision, allocatable, dimension(:) ::  x_f, y_f, z_f
    integer, allocatable, dimension(:,:) :: index_cell_f
    double precision, allocatable, dimension(:,:) :: x_filter, &
                                                     y_filter, &
                                                     z_filter
    double precision, allocatable, dimension(:) :: filter_vol
    double precision, allocatable, dimension(:) :: filter_poro

    double precision, dimension(:), pointer :: subchannel_ID
    double precision, dimension(:), pointer :: filter_ID
    double precision, dimension(:), pointer :: porosity_v

    logical :: cond1, cond2

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('VolumePorosity', porosity_v)

    call mesh_dimensions(solid_mesh_volume, n_node, n_entity)

!!  for pure tetrahedral solid mesh
    n_vertices = 4

    allocate(x(n_node),y(n_node),z(n_node))
    allocate(index_cell(n_entity, n_vertices))

    allocate(xx(n_entity,n_vertices))
    allocate(yy(n_entity,n_vertices))
    allocate(zz(n_entity,n_vertices))
    allocate(cell_xyz(n_entity,3))
    allocate(cell_vol(n_entity))

    call auxiliary_mesh_readin(solid_mesh_volume, n_node, n_entity, &
                               n_vertices, x, y, z, index_cell)

!!  Scale mesh from mm to m
    if(solid_mesh_scale_v)then
      x = x*solid_scaling_factor_v
      y = y*solid_scaling_factor_v
      z = z*solid_scaling_factor_v
    end if

!!  Calculate coordinates of the cell centers (solid mesh cells)
    do i=1, n_entity
      xx(i,:) = x(index_cell(i,:))
      yy(i,:) = y(index_cell(i,:))
      zz(i,:) = z(index_cell(i,:))
      cell_xyz(i,1) = sum(xx(i,:))/n_vertices
      cell_xyz(i,2) = sum(yy(i,:))/n_vertices
      cell_xyz(i,3) = sum(zz(i,:))/n_vertices
    end do

!!  Calculate cell volumes of the solid mesh cells
    do i=1, n_entity
      vec_A(1) = xx(i,2)-xx(i,1)
      vec_A(2) = yy(i,2)-yy(i,1)
      vec_A(3) = zz(i,2)-zz(i,1)
      vec_B(1) = xx(i,3)-xx(i,1)
      vec_B(2) = yy(i,3)-yy(i,1)
      vec_B(3) = zz(i,3)-zz(i,1)
      vec_C(1) = xx(i,4)-xx(i,1)
      vec_C(2) = yy(i,4)-yy(i,1)
      vec_C(3) = zz(i,4)-zz(i,1)
      vec_AxB(1) = vec_A(2)*vec_B(3)-vec_A(3)*vec_B(2)  
      vec_AxB(2) = vec_A(3)*vec_B(1)-vec_A(1)*vec_B(3)
      vec_AxB(3) = vec_A(1)*vec_B(2)-vec_A(2)*vec_B(1)
      cell_vol(i) = vec_AxB(1)*vec_C(1) &
                  + vec_AxB(2)*vec_C(2) &
                  + vec_AxB(3)*vec_C(3)
      cell_vol(i) = abs(1.d0/6.d0*cell_vol(i))
    end do

!! Prosity calculation method 1 (subchannel based)
!! Find the xyz range of the solid mesh
    solid_min(1) = minval(x)
    solid_min(2) = minval(y)
    solid_min(3) = minval(z)
    solid_max(1) = maxval(x)
    solid_max(2) = maxval(y)
    solid_max(3) = maxval(z)

!! Find subchannel indices to be looped over
    n = 0
    do ii=1, n_cell
      if(x_channel_min(ii) <= solid_max(1) .and. &
         y_channel_min(ii) <= solid_max(2) .and. &
         z_channel_min(ii) <= solid_max(3) .and. &
         x_channel_max(ii) >= solid_min(1) .and. &
         y_channel_max(ii) >= solid_min(2) .and. &
         z_channel_max(ii) >= solid_min(3)) then

        n = n+1
        index_sub(n) = ii
      end if
    end do

!! Find the subchannel index for each control volume (only for square lattice
!  rod bundles)
    sub_poro = 1.d0
    do ii=1, n
      solid_vol = 0.d0
      j = index_sub(ii)
      do i=1, n_entity
        if(cell_xyz(i,1) >= x_channel_min(j) .and. &
           cell_xyz(i,1) <= x_channel_max(j) .and. &
           cell_xyz(i,2) >= y_channel_min(j) .and. &
           cell_xyz(i,2) <= y_channel_max(j) .and. &
           cell_xyz(i,3) >= z_channel_min(j) .and. &
           cell_xyz(i,3) <= z_channel_max(j)) then

          solid_vol = solid_vol+cell_vol(i)
        end if
      end do
      sub_poro(j) = 1.d0-solid_vol/subchannel_vol(j)
    end do

    if(iporo_cal==1)then
      do iel=1, ncel
        if(subchannel_ID(iel)>0.d0) &
          porosity_v(iel) = sub_poro(int(subchannel_ID(iel)))
      end do

    end if

!! Prosity calculation method 2 (cell based)
!! Parallel run may be not as trustable as that of serial run depending on the 
!! quality of the solid mesh
    if(iporo_cal==2)then
!! Find face indices to be looped over
      flag_fac = 0
      n_fac = 0
      n_fbr = 0
      do ii=1, n
        j = index_sub(ii)
        do i=1, nfac
          cond1 = .false.
          cond2 = .false.
          iel1 = ifacel(1,i)
          iel2 = ifacel(2,i)
          if(subchannel_ID(iel1)==j) cond1 = .true.
          if(subchannel_ID(iel2)==j) cond2 = .true.
          if(cond1 .and. cond2) then
            n_fac = n_fac+1
            index_fac(n_fac) = i
          end if
          if(cond1 .or. cond2)then
            if(.not.(cond1 .and. cond2))then
              if(flag_fac(i)==0)then
                n_fac = n_fac+1
                flag_fac(i) = 1
                index_fac(n_fac) = i
              end if
            end if
          end if
        end do

        do i=1,nfabor
          iel = ifabor(i)
          if(subchannel_ID(iel)==j) then
            n_fbr = n_fbr+1
            index_fbr(n_fbr) = i
          end if  
        end do

      end do

      write(*,'(A10, I4, A6, I8, A6, I8)')'MPI rank: ', irangp, 'n_fac', &
                                            n_fac, 'n_fbr', n_fbr

      solid_vol_cell = 0
      do ii=1, n_fac
        ifac = index_fac(ii)
        iel1 = ifacel(1,ifac)
        iel2 = ifacel(2,ifac)

        if(iel1<=ncel)then
          call solid_vol_cell_interior(ifac, iel1, cell_xyz, cell_vol, &
                                       solid_vol_pyrmd)

          solid_vol_cell(iel1) = solid_vol_cell(iel1)+solid_vol_pyrmd
        end if

        if(iel2<=ncel)then
          call solid_vol_cell_interior(ifac, iel2, cell_xyz, cell_vol, &
                                       solid_vol_pyrmd)

          solid_vol_cell(iel2) = solid_vol_cell(iel2)+solid_vol_pyrmd
        end if

        write(*,'(A10, I4, A8, 3I8)')'MPI rank: ', irangp, 'i_fac', ii, &
                                      iel1, iel2

      end do

      do ii=1, n_fbr
        ifac = index_fbr(ii)
        iel = ifabor(ifac)

        if(iel<=ncel)then
          call solid_vol_cell_boundary(ifac, iel, cell_xyz, cell_vol, &
                                       solid_vol_pyrmd)

          solid_vol_cell(iel) = solid_vol_cell(iel)+solid_vol_pyrmd
        end if

        write(*,'(A10, I4, A8, 2I8)')'MPI rank: ', irangp, 'i_fbr', ii, iel

      end do

      do iel=1, ncel
        porosity_v(iel) = 1.d0-solid_vol_cell(iel)/volume(iel)
      end do

    end if

!! Porosity method calculation method 3 (local flitering mesh based)
    if(iporo_cal==3)then

      call field_get_val_s_by_name('FilterID', filter_ID)

      call mesh_dimensions(filtering_mesh, n_node_f, n_entity_f)

!!  for pure hexahedral filtering mesh
      n_vertices_f = 8

      allocate(x_f(n_node_f),y_f(n_node_f),z_f(n_node_f))
      allocate(index_cell_f(n_entity_f, n_vertices_f))

      allocate(filter_cell_vol(1,n_entity_f))

      allocate(x_filter(n_entity_f,2))
      allocate(y_filter(n_entity_f,2))
      allocate(z_filter(n_entity_f,2))

      allocate(filter_vol(n_entity_f))
      allocate(filter_poro(n_entity_f))

      call auxiliary_mesh_readin(filtering_mesh, n_node_f, n_entity_f, &
                                 n_vertices_f, x_f, y_f, z_f, index_cell_f)

!!  Scale mesh from mm to m
      if(filtering_mesh_scale)then
        x_f = x_f*filtering_scaling_factor
        y_f = y_f*filtering_scaling_factor
        z_f = z_f*filtering_scaling_factor
      end if

!!  Calculate volume of the filter blocks
      do i=1, n_entity_f
        x_filter(i,1) = minval(x_f(index_cell_f(i,:)))
        x_filter(i,2) = maxval(x_f(index_cell_f(i,:)))
        y_filter(i,1) = minval(y_f(index_cell_f(i,:)))
        y_filter(i,2) = maxval(y_f(index_cell_f(i,:)))
        z_filter(i,1) = minval(z_f(index_cell_f(i,:)))
        z_filter(i,2) = maxval(z_f(index_cell_f(i,:)))
        filter_vol(i) = (x_filter(i,2)-x_filter(i,1)) &
                      * (y_filter(i,2)-y_filter(i,1)) &
                      * (z_filter(i,2)-z_filter(i,1))
        filter_cell_vol(1,i) = filter_vol(i)
      end do

      filter_poro = 1.d0
      do ii=1, n_entity_f
        solid_vol = 0.d0
        do i=1, n_entity
          if(cell_xyz(i,1) >= x_filter(ii,1) .and. &
             cell_xyz(i,1) <= x_filter(ii,2) .and. &
             cell_xyz(i,2) >= y_filter(ii,1) .and. &
             cell_xyz(i,2) <= y_filter(ii,2) .and. &
             cell_xyz(i,3) >= z_filter(ii,1) .and. &
             cell_xyz(i,3) <= z_filter(ii,2)) then
           
            solid_vol = solid_vol+cell_vol(i)
          end if
        end do
        filter_poro(ii) = 1.d0-solid_vol/filter_vol(ii)
      end do

      filter_ID = -1.d0
      do iel=1, ncel
        do i=1, n_entity_f
          if(xyzcen(1,iel)>=x_filter(i,1) &
            .and.xyzcen(1,iel)<=x_filter(i,2) &
            .and.xyzcen(2,iel)>=y_filter(i,1) &
            .and.xyzcen(2,iel)<=y_filter(i,2) &
            .and.xyzcen(3,iel)>=z_filter(i,1) &
            .and.xyzcen(3,iel)<=z_filter(i,2))then
            
            filter_ID(iel) = i
          end if
        end do
      end do

      porosity_v = 1.d0
      do iel=1, ncel
        if(filter_ID(iel)>0.d0) &
          porosity_v(iel) = filter_poro(int(filter_ID(iel)))
      end do

      deallocate(x_f,y_f,z_f)
      deallocate(index_cell_f)

      deallocate(x_filter)
      deallocate(y_filter)
      deallocate(z_filter)
      
      deallocate(filter_vol)
      deallocate(filter_poro)

    end if

!! Write volume porosity field to file (in form of global cell indices)
    mode = 2
    call io_porosity(porosity_v, mode)


!!  Output for checking results 
!    if(irangp.le.0)write(*,*)n_node, n_entity
!    do i=1,n_entity
!      if(irangp.le.0 .and. i<=10)&
!        write(*,*)i, x(i), y(i), z(i)
!!        write(*,*)i, index_cell(i,:)
!    end do

!    open(200,file='test1.txt')
!    do iel=1, ncel 
!      write(200,*)iel, subchannel_ID(iel)
!    end do
!    close(200)

    deallocate(x,y,z)
    deallocate(index_cell)

    deallocate(xx,yy,zz)
    deallocate(cell_xyz)
    deallocate(cell_vol)

!! Testing subroutines (derive mesh quantities from C routines)
!    call face_vtx_test()
!    call userroutine()

!! Testing determinant calculation
!    call recursion_det()

!! Testing point_tetra function
!    call point_tetra_test()

  end subroutine volume_porosity


!===============================================================================
!! This routine is to read in the filtering meshes of the porous sub-models
  subroutine porous_filtering_mesh()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: i, iel
    
    integer :: n_node, n_entity, n_vertices
    integer, allocatable, dimension(:,:) :: index_cell
    
    double precision, allocatable, dimension(:) ::  x, y, z
    double precision, allocatable, dimension(:,:) :: x_filter, &
                                                     y_filter, &
                                                     z_filter
    double precision, allocatable, dimension(:) :: filter_vol

    double precision, dimension(:), pointer :: filter_ID

    call field_get_val_s_by_name('FilterID', filter_ID)

    call mesh_dimensions(filtering_mesh, n_node, n_entity)

!!  for pure hexahedral filtering mesh
    n_vertices = 8

    allocate(x(n_node),y(n_node),z(n_node))
    allocate(index_cell(n_entity, n_vertices))

    allocate(x_filter(n_entity,2))
    allocate(y_filter(n_entity,2))
    allocate(z_filter(n_entity,2))
    allocate(filter_vol(n_entity))

    if(.not.allocated(filter_cell_vol)) allocate(filter_cell_vol(1,n_entity))

    call auxiliary_mesh_readin(filtering_mesh, n_node, n_entity, &
                               n_vertices, x, y, z, index_cell)

!!  Scale mesh from mm to m
    if(filtering_mesh_scale)then
      x = x*filtering_scaling_factor
      y = y*filtering_scaling_factor
      z = z*filtering_scaling_factor
    end if

!!  Calculate volume of the filter blocks
    do i=1, n_entity
      x_filter(i,1) = minval(x(index_cell(i,:)))
      x_filter(i,2) = maxval(x(index_cell(i,:)))
      y_filter(i,1) = minval(y(index_cell(i,:)))
      y_filter(i,2) = maxval(y(index_cell(i,:)))
      z_filter(i,1) = minval(z(index_cell(i,:)))
      z_filter(i,2) = maxval(z(index_cell(i,:)))
      filter_vol(i) = (x_filter(i,2)-x_filter(i,1)) &
                    * (y_filter(i,2)-y_filter(i,1)) &
                    * (z_filter(i,2)-z_filter(i,1))
      filter_cell_vol(1,i) = filter_vol(i)
    end do

    filter_ID = -1.d0
    do iel=1, ncel
      do i=1, n_entity
        if(xyzcen(1,iel)>=x_filter(i,1) &
           .and.xyzcen(1,iel)<=x_filter(i,2) &
           .and.xyzcen(2,iel)>=y_filter(i,1) &
           .and.xyzcen(2,iel)<=y_filter(i,2) &
           .and.xyzcen(3,iel)>=z_filter(i,1) &
           .and.xyzcen(3,iel)<=z_filter(i,2))then
            
          filter_ID(iel) = i
        end if
      end do
    end do 

    deallocate(x,y,z)
    deallocate(index_cell)

    deallocate(x_filter)
    deallocate(y_filter)
    deallocate(z_filter)
    deallocate(filter_vol)

  end subroutine porous_filtering_mesh


!===============================================================================
!! This routine is to read or write volume porosities from file 
!! volume_porosity.dat
  subroutine io_porosity(porosity_v, mode)

    use numvar
    use field
    use mesh
    use parall

    implicit none

    interface
      subroutine global_cell_id(c_id, c_g_id)        &
        bind(C, name='global_cell_id')
        use, intrinsic :: iso_c_binding
        implicit none
        integer(kind=c_int) :: c_id
        integer(kind=c_int) :: c_g_id
      end subroutine global_cell_id
    end interface

    double precision, dimension(:) :: porosity_v
    integer :: mode

    integer :: i, ii, iel
    double precision, allocatable, dimension(:) :: porosity_v_g
    integer :: c_g_id
    integer :: dummy

    character(len=40) :: heading
    double precision :: const_poro
    double precision, dimension(:), pointer :: filter_ID

    integer :: aux, bcast

    allocate(porosity_v_g(ncelgb))

!!  read in porosities
    if(mode==1) then
      open (100, file=v_porous, status='old')

      read(100, *)
      read(100, *) heading
      if (trim(heading)=='constant'.and.iporo_cal==3)then
        write(*,*)'constant volume porosity used for Zone: 1'
        read(100, *) const_poro
        call field_get_val_s_by_name('FilterID', filter_ID)
        do iel=1, ncel
          if(filter_ID(iel)>0.d0)then
            porosity_v(iel) = const_poro
          else
            porosity_v(iel) = 1.d0
          end if
        end do
        return

      else
        backspace(100)
      end if

      do i=1, ncelgb
        read(100, *) dummy, porosity_v_g(i)
      end do
      close(100)

      do iel=1, ncel
        if (irangp.ge.0) then
          call global_cell_id(iel, c_g_id)
          porosity_v(iel) = porosity_v_g(c_g_id)
        else
          porosity_v(iel) = porosity_v_g(iel)
        end if
      end do

    end if

!!  write porosities to file
    if(mode==2) then

      porosity_v_g = 1.5d0

      do iel=1, ncel
        if (irangp.ge.0) then
          call global_cell_id(iel, c_g_id)
          porosity_v_g(c_g_id) = porosity_v(iel)
        else
          porosity_v_g(iel) = porosity_v(iel)
        endif

      end do

      if (irangp.ge.0) then

        do i=1, ncelgb
          aux = 1
          call parfpt (aux, bcast, porosity_v_g(i))
          call parall_bcast_r(bcast, porosity_v_g(i))

        end do
      end if

      if(irangp.le.0)then
        open (100, file=v_porous, status='replace')
        write(100, '(A30, I1)') 'Porosity calculation method:', iporo_cal
        do i=1, ncelgb
          write(100,*) i, porosity_v_g(i)
        end do
        close(100)
      end if
    end if

    deallocate(porosity_v_g)

  end subroutine io_porosity


!===============================================================================
!! This is a routine to test calling a C function(routine) from a Fortan routine
!! to access mesh quantities
  subroutine face_vtx_test()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    interface
      subroutine n_face_vtx(f_id, n_face_vertices)        &
        bind(C, name='n_face_vtx')
        use, intrinsic :: iso_c_binding
        implicit none
        integer(kind=c_int) :: f_id
        integer(kind=c_int) :: n_face_vertices
      end subroutine n_face_vtx
    end interface

    interface
      subroutine cord_face_vtx(f_id, n_face_vertices, cord_vtx)        &
        bind(C, name='cord_face_vtx')
        use, intrinsic :: iso_c_binding
        implicit none
        integer(kind=c_int) :: f_id
        integer(kind=c_int) :: n_face_vertices
        real(kind=c_double),dimension(*) :: cord_vtx
      end subroutine cord_face_vtx
    end interface

    integer :: nface_vtx
    integer :: face_id = 50
    double precision, allocatable, dimension(:,:) :: cord_vtx

    double precision, dimension(3) :: va, vb, vc
    double precision :: area

    call n_face_vtx(face_id, nface_vtx)
    write(*,*) nface_vtx
!    write(*,*) ipnfac(face_id+1) !counterpart of face_vtx_idx[fid] in fortran
                                  !note every fortran index = c index+1

    allocate (cord_vtx(nface_vtx,3))
    call cord_face_vtx(face_id, nface_vtx, cord_vtx)
!    write(*,*) cord_vtx(4,1)
!    write(*,*) nodfac(ipnfac(face_id+1)+1) !counterpart of face_vtx[s_id+1] 
                                            !in fortran note every fortran 
                                            !index = c index+1
    deallocate (cord_vtx)

  end subroutine face_vtx_test


!===============================================================================
!! This is a routine to test calling a C function(routine) from a Fortan routine
!! It should be noted that c arrays start from 0 whereas Fortran arrays start
!! from 1. So, converting a c array index to a fortran array index should add
!! the c array index by 1 to obtain the corresponding Fortran array index. 
!! In addition, a 2d index e.g. (i,j) in a c array corresponds to (j,i) in a
!! Fortran array  
  subroutine userroutine()

    use parall

    implicit none

    interface
      subroutine modifyarray(insv) BIND(C, name='modifyarray')
        use, intrinsic :: iso_c_binding
        real(c_double),dimension(*), intent(INOUT) :: insv
      end subroutine modifyarray
    end interface

    double precision, allocatable, dimension(:,:) :: insv
    integer :: i,j

    allocate (insv(10,3))
    insv = 2.d0

    write(*,*) insv(1,1)
    call modifyarray(insv)
    if (irangp.le.0) then
      do i=1,10
        do j=1,3
          write(*,*)i,j, insv(i,j)
        end do
      end do
    end if

    deallocate(insv)

  end subroutine userroutine


!===============================================================================
!! This routine to calculate the solid volume of a pyramid formed by the 
!! boundary face: ifac and the cell center of cell: iel
  subroutine solid_vol_cell_interior(ifac, iel, xyzcen_solid, cell_vol_solid, &
                                     vol_solid)

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: ifac
    integer :: iel
    double precision, dimension(:,:) :: xyzcen_solid
    double precision, dimension(:) :: cell_vol_solid
    double precision :: vol_solid

    integer :: i, ii
    integer :: s_id, e_id
    integer :: nface_vtx
    integer :: dim1
    double precision, allocatable, dimension(:,:) :: cord_vtx

    double precision, allocatable, dimension(:) :: area
    double precision :: total_area
    double precision :: dis0, dis1, dis2, dis3

    integer :: in_out

    s_id = ipnfac(ifac)
    e_id = ipnfac(ifac+1)
    nface_vtx = e_id-s_id

    allocate (cord_vtx(3,nface_vtx))

    do ii = 1, nface_vtx
      cord_vtx(:,ii) = xyznod(:,nodfac(s_id+ii-1))
    end do

!!  The vertices in a face are in good order
    allocate(area(nface_vtx-2))
    do ii=1, nface_vtx-2
      area(ii) = tri_area(cord_vtx(:,1), cord_vtx(:,ii+1), cord_vtx(:,ii+2))
    end do
    total_area = sum(area)
    if(abs(surfan(ifac)-total_area)>1.e-10) then
      write(*,*) 'Warning! the vertices of interior face: ',ifac, ' are &
                  &not in good order!'
      write(*,*) ifac, surfan(ifac), total_area
    end if
    deallocate(area)

    vol_solid = 0.d0
    dim1 = size(xyzcen_solid,1)
    do ii=1, nface_vtx-2
      do i=1, dim1
        in_out = -1

        dis0 = dis_points(xyzcen_solid(i,:), xyzcen(:,iel))
        dis1 = dis_points(cord_vtx(:,1), xyzcen(:,iel))
        dis2 = dis_points(cord_vtx(:,ii+1), xyzcen(:,iel))
        dis3 = dis_points(cord_vtx(:,ii+2), xyzcen(:,iel))

        if(dis0<=dis1 .or. dis0<=dis2 .or. dis0<=dis3) then
          in_out = point_tetra(xyzcen_solid(i,:), cord_vtx(:,1),   &
                               cord_vtx(:,ii+1), cord_vtx(:,ii+2), &
                               xyzcen(:,iel))
          if(in_out==1) vol_solid = vol_solid+cell_vol_solid(i)
          if(in_out==0) write(*,*)'Warning, abnormal solid cell detected!'
        end if

      end do
    end do

    deallocate(cord_vtx)

  end subroutine solid_vol_cell_interior


!===============================================================================
!! This routine to calculate the solid volume of a pyramid formed by the
!! interior face: ifac and the cell center of cell: iel
  subroutine solid_vol_cell_boundary(ifac, iel, xyzcen_solid, cell_vol_solid, &
                                     vol_solid)

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: ifac
    integer :: iel
    double precision, dimension(:,:) :: xyzcen_solid
    double precision, dimension(:) :: cell_vol_solid
    double precision :: vol_solid

    integer :: i, ii
    integer :: s_id, e_id
    integer :: nface_vtx
    integer :: dim1
    double precision, allocatable, dimension(:,:) :: cord_vtx

    double precision, allocatable, dimension(:) :: area
    double precision :: total_area
    double precision :: dis0, dis1, dis2, dis3

    integer :: in_out

    s_id = ipnfbr(ifac)
    e_id = ipnfbr(ifac+1)
    nface_vtx = e_id-s_id

    allocate (cord_vtx(3,nface_vtx))

    do ii = 1, nface_vtx
      cord_vtx(:,ii) = xyznod(:,nodfbr(s_id+ii-1))
    end do

!!  The vertices in a face are in good order
    allocate(area(nface_vtx-2))
    do ii=1, nface_vtx-2
      area(ii) = tri_area(cord_vtx(:,1), cord_vtx(:,ii+1), cord_vtx(:,ii+2))
    end do  
    total_area = sum(area)
    if(abs(surfbn(ifac)-total_area)>1.e-10) then
      write(*,*) 'Warning! the vertices of boundary face: ',ifac, ' are &
                  &not in good order!'
      write(*,*) ifac, surfbn(ifac), total_area
    end if
    deallocate(area)

    vol_solid = 0.d0
    dim1 = size(xyzcen_solid,1)
    do ii=1, nface_vtx-2
      do i=1, dim1
        in_out = -1

        dis0 = dis_points(xyzcen_solid(i,:), xyzcen(:,iel))
        dis1 = dis_points(cord_vtx(:,1), xyzcen(:,iel))
        dis2 = dis_points(cord_vtx(:,ii+1), xyzcen(:,iel))
        dis3 = dis_points(cord_vtx(:,ii+2), xyzcen(:,iel))

        if(dis0<=dis1 .or. dis0<=dis2 .or. dis0<=dis3) then
          in_out = point_tetra(xyzcen_solid(i,:), cord_vtx(:,1),   &
                               cord_vtx(:,ii+1), cord_vtx(:,ii+2), &
                               xyzcen(:,iel))
          if(in_out==1) vol_solid = vol_solid+cell_vol_solid(i)
          if(in_out==0) write(*,*)'Warning, abnormal solid cell detected!'
        end if 

      end do
    end do

    deallocate(cord_vtx)

  end subroutine solid_vol_cell_boundary


!===============================================================================
!! This subroutine is a simple test on point_tetra
  subroutine point_tetra_test()

    implicit none

    double precision, dimension(3) :: point, tet_a, tet_b, tet_c, tet_d
    data point /0.1d0, 0.1d0, -0.1d0/
    data tet_a /1.d0, 0.d0, 0.d0/
    data tet_b /0.d0, 1.d0, 0.d0/
    data tet_c /0.d0, 0.d0, 0.d0/
    data tet_d /0.d0, 0.d0, 1.d0/

    integer :: in_out

    in_out = point_tetra(point, tet_a, tet_b, tet_c, tet_d)

    if(in_out==1) then
      write(*,*)'Point is inside the tetrahedron.'
    else if(in_out==0) then
      write(*,*)'Point is on the surface of the tetrahedron.'
    else
      write(*,*)'Point is outside the tetrahedron.'
    end if

  end subroutine point_tetra_test


!===============================================================================
!! This function is to judge a point is inside/outside a tetrahedron
  function point_tetra(point, tet_a, tet_b, tet_c, tet_d)

    implicit none

    double precision, dimension(:) :: point, tet_a, tet_b, tet_c, tet_d
    integer :: point_tetra

    double precision, dimension(4,4) :: d0, d1, d2, d3, d4
    double precision :: det_d0, det_d1, det_d2, det_d3, det_d4

    d0 = 1.d0
    d1 = 1.d0
    d2 = 1.d0
    d3 = 1.d0
    d4 = 1.d0

    d0(1, 1:3) = tet_a(1:3)
    d0(2, 1:3) = tet_b(1:3)
    d0(3, 1:3) = tet_c(1:3)
    d0(4, 1:3) = tet_d(1:3)

    d1(1, 1:3) = point(1:3)
    d1(2, 1:3) = tet_b(1:3)
    d1(3, 1:3) = tet_c(1:3)
    d1(4, 1:3) = tet_d(1:3)

    d2(1, 1:3) = tet_a(1:3)
    d2(2, 1:3) = point(1:3)
    d2(3, 1:3) = tet_c(1:3)
    d2(4, 1:3) = tet_d(1:3)

    d3(1, 1:3) = tet_a(1:3)
    d3(2, 1:3) = tet_b(1:3)
    d3(3, 1:3) = point(1:3)
    d3(4, 1:3) = tet_d(1:3)

    d4(1, 1:3) = tet_a(1:3)
    d4(2, 1:3) = tet_b(1:3)
    d4(3, 1:3) = tet_c(1:3)
    d4(4, 1:3) = point(1:3)

    det_d0 = det_matrix(d0)
    det_d1 = det_matrix(d1)
    det_d2 = det_matrix(d2)
    det_d3 = det_matrix(d3)
    det_d4 = det_matrix(d4)

    if(det_d0*det_d1>0.d0 .and. det_d1*det_d2>0.d0 .and. det_d2*det_d3>0.d0 &
       .and. det_d3*det_d4>0.d0) then
      point_tetra = 1
    else if (det_d0*det_d1*det_d2*det_d3*det_d4==0.d0) then
      point_tetra = 0
    else
      point_tetra = -1
    end if

  end function point_tetra


!===============================================================================
!! This subroutine is a simple test on det_matrix
  subroutine recursion_det()

    implicit none

    double precision :: answer
    double precision, dimension(4,4) :: matrix
    data matrix /1,2,1,2,1,1,1,2,3,1,0,0,2,1,2,0/

    answer = det_matrix(matrix)

    write(*,*)'result=',answer

    return

  end subroutine recursion_det


!===============================================================================
!! This function is to calculate the determinant of a matrix
  recursive double precision function det_matrix(input)result(ans)

    implicit none

    integer::i,j,k
    integer::dim
    double precision res

    double precision, dimension(:,:) :: input
    double precision, allocatable, dimension(:,:) :: sub
    double precision, allocatable, dimension(:) :: multiply

    dim=size(input,1)

    allocate(multiply(dim))
    allocate(sub(dim-1,dim-1))

    ans=0.

    do j=1,dim
      multiply(j)=input(1,j)
      do i=1,dim-1
        do k=1,dim-1
          if (k<j)then
            sub(i,k)=input(i+1,k)
          else
            sub(i,k)=input(i+1,k+1) 
          end if
        end do
      end do

      if(dim==1) then
        ans=multiply(j)
        return
      end if

      ans=ans+multiply(j)*((-1)**(j+1))*det_matrix(sub)
    end do

    return

  end function det_matrix


!===============================================================================
!! This function is to calculate the surface area of a triangle by coorinates of 
!! the vertices
  function tri_area(va, vb, vc)

    implicit none

    double precision, dimension(3), intent(IN) :: va, vb, vc
    double precision :: tri_area

    double precision, dimension(3) :: vector_a, vector_b, vec_axb

    vector_a = vb-va
    vector_b = vc-va

    vec_axb(1) = vector_a(2)*vector_b(3)-vector_a(3)*vector_b(2)
    vec_axb(2) = vector_a(3)*vector_b(1)-vector_a(1)*vector_b(3)
    vec_axb(3) = vector_a(1)*vector_b(2)-vector_a(2)*vector_b(1)

    tri_area = 0.5d0*sqrt(vec_axb(1)**2+vec_axb(2)**2+vec_axb(3)**2)

  end function tri_area


!===============================================================================
!! This function is to calculate the distance of two points
  function dis_points(va, vb)

    implicit none

    double precision, dimension(3), intent(IN) :: va, vb
    double precision :: dis_points

    double precision, dimension(3) :: vector

    vector = vb-va

    dis_points= sqrt(vector(1)**2+vector(2)**2+vector(3)**2)

  end function dis_points

 
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
    double precision, allocatable, dimension(:) :: store

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
      if(subchannel_ID(iel)>0.d0)then
        Dh(iel) = 4.d0*subchannel_vol(int(subchannel_ID(iel)))&
                / subchannel_surf(int(subchannel_ID(iel)))
        surf(iel) = subchannel_surf(int(subchannel_ID(iel)))
        surf_rod(iel) = subchannel_rodsurf(int(subchannel_ID(iel)))
        surf_wall(iel) = subchannel_wallsurf(int(subchannel_ID(iel)))
        cvol(iel) = subchannel_vol(int(subchannel_ID(iel)))
      end if
    end do

!! Derive hydralic diameter from file (porous SubChCFD for 7x7 ballooned bundle)
!  allocate (store(n_cell))
!  open (400, file='hydraulic_diameter.dat')
!  do i=1, n_cell
!    read(400,*) store(i)
!  end do
!  do iel=1, ncel
!    Dh(iel) = store(int(subchannel_ID(iel)))
!  end do
!  close(400)
!  deallocate(store)

!! write hydraulic diameter to file (porous SubChCFD for 7x7 ballooned bundle)
!  open(300, file='hydraulic_diameter.dat')
!  do i=1,n_cell
!    write(300,*) 4.d0*subchannel_vol(i)/subchannel_surf(i)  
!  end do
!  close(300)


!! Output for checking results
!    write(*,*) n_cell
!    write(*,*) cell_height
!    open(200,file='test.txt')
!    do i=1, n_cell
!       write(200,*) i, subchannel_surf(i), nf(i)
!       write(200,*) i, subchannel_vol(i), nc(i)
!    end do
!    close(200)

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
      if(subchannel_ID(iel)>0.d0)then
        axial_Re(iel) = A_den(iel)*A_vel(3,iel)*Dh(iel)/A_vis(iel)
        axial_Re(iel) = abs(axial_Re(iel))
      end if
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
      if(subchannel_ID(iel)>0.d0)then
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
        case(3)
          if(axial_Re(iel)<1.d-30)then
            axial_f(iel) = 1.d30
          else
            axial_f(iel) = FRICTIONAL_FACTOR_HTGR_GA(subchannel_ID(iel), &
                                                     axial_Re(iel))
            if(user_frictional_factor) &
              axial_f(iel) = FRICTIONAL_FACTOR_USER(iel)
          end if
        end select
      end if
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

!!  Convert the Fanning frictional factor to Darcy frictional factor
    cf = 4.d0*f
    FRICTIONAL_FACTOR_AGR = cf

  end function FRICTIONAL_FACTOR_AGR


!===============================================================================
!! This is for HTGR (GA design)
  function FRICTIONAL_FACTOR_HTGR_GA(sub_ID, Re)

    use numvar
    use field
    use mesh
    use parall
    use optcal

    implicit none

    double precision :: FRICTIONAL_FACTOR_HTGR_GA
!! Aguments
    double precision :: sub_ID, Re

!! Local variables
    double precision :: cf

    integer :: isub

!! For fully developed flow in smooth circular channels
    if(Re<=2.3d3)then
      cf = 64.d0/Re
!! Blasius
    else if(Re>=2.3d3 .and. Re<3.d4) then
      cf = 0.316*(Re**(-0.25d0))
!! McAdams
    else if(Re>=3.d4 .and. Re<1.d6) then
      cf = 0.184*(Re**(-0.2d0))
!! To be added
    else
      cf = 0.184*(Re**(-0.2d0))
    end if

    isub = int(sub_ID)

    FRICTIONAL_FACTOR_HTGR_GA = cf

  end function FRICTIONAL_FACTOR_HTGR_GA


!===============================================================================
!! This routine is used to calculate suchannel Nusselt number
  subroutine subchannel_Nu()

    use numvar
    use optcal
    use field
    use mesh
    use parall

    implicit none

    integer :: iel
    double precision, dimension(:), pointer :: axial_Re, Nu, subchannel_ID, &
                                               Dh, A_Tbulk, A_twall
    double precision, dimension(:), pointer :: A_vis, A_cond, A_cp
    double precision :: temperature_ratio, Pr

    if (iscalt.le.0) return

    call field_get_val_s_by_name('AxialReynoldsNumber', axial_Re)
    call field_get_val_s_by_name('NusseltNumber', Nu)
    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('HydraulicDiameter', Dh)
    call field_get_val_s_by_name('SubchannelBulkTemperature', A_Tbulk)
    call field_get_val_s_by_name('SubchannelWallTemperature', A_twall)
    call field_get_val_s_by_name('AverageDynamicViscosity', A_vis)
    call field_get_val_s_by_name('AverageThermalConductivity', A_cond)
    call field_get_val_s_by_name('AverageHeatCapacity', A_cp)

    do iel=1, ncel
      if(subchannel_ID(iel)>0.d0)then
        select case (bundle_type)
        case (1)
          if(axial_Re(iel)<1.d-30)then
            Nu(iel) = 4.36
          else
            Pr = A_vis(iel)*A_cp(iel)/A_cond(iel)
            Nu(iel) = NUSSELT_NUMBER_PWR(subchannel_ID(iel), axial_Re(iel), Pr,&
                                         Dh(iel), heating, 1)*crr_Nu
            if(user_nusselt_number) &
              Nu(iel) = NUSSELT_NUMBER_USER(iel)
          end if
        case (2)
          temperature_ratio = A_twall(iel)/A_Tbulk(iel)
          if(axial_Re(iel)<1.d-30)then
            Nu(iel) = 4.36
          else
            Nu(iel) = NUSSELT_NUMBER_AGR(voidage, temperature_ratio, &
                                         subchannel_ID(iel), axial_Re(iel), 1) &
                                         *crr_Nu
            if(user_nusselt_number) &
              Nu(iel) = NUSSELT_NUMBER_USER(iel)
          end if
        case (3)
          if(axial_Re(iel)<1.d-30)then
            Nu(iel) = 4.36
          else
            Pr = A_vis(iel)*A_cp(iel)/A_cond(iel)
            Nu(iel) = NUSSELT_NUMBER_HTGR_GA(subchannel_ID(iel), axial_Re(iel),&
                                             Pr, heating)*crr_Nu
            if(user_nusselt_number) &
              Nu(iel) = NUSSELT_NUMBER_USER(iel)
          end if
        end select
      end if
    end do

  end subroutine subchannel_Nu


!===============================================================================
!! This is for PWR
  function NUSSELT_NUMBER_PWR(sub_ID, Re, Pr, Dh, heating, direction)

    use numvar
    use field
    use mesh
    use parall

    implicit none

    double precision :: NUSSELT_NUMBER_PWR
!! Aguments
    double precision :: sub_ID, Re, Pr, Dh
    logical :: heating
    integer :: direction

!! Local variables 
    double precision :: coef, exp_B, Nusselt, Nusselt0
    logical :: condition1, condition2, condition3, condition
    integer :: isub

    isub = int(sub_ID)

    if (heating)then
      Nusselt0 = 0.023d0*(Re**0.8d0)*(Pr**0.4d0)
    else
      Nusselt0 = 0.023d0*(Re**0.8d0)*(Pr**0.3d0)
    end if

    exp_B = -Dh/rod_diameter
    coef = 1.d0+0.912d0*(Re**(-0.1d0))*(Pr**0.4d0)*(1.d0-2.0043d0*(ee**exp_B))

    condition1 = Re>=3.d3 .and. Re<=1.d6
    condition2 = Pr>=0.66d0 .and. Pr<=5.d0
    condition3 = PbyD(isub)>=1.d0 .and. PbyD(isub)<=1.8d0
    condition = condition1 .and. condition2 .and. condition3

    Nusselt = coef*Nusselt0
    if(.not.condition) then
      !print*, 'Warning! Nusselt number calculation error: out of range at ',&
      !        'Subchannel ID ', isub
    end if
    NUSSELT_NUMBER_PWR = Nusselt

  end function NUSSELT_NUMBER_PWR


!===============================================================================
!! This is for AGR
  function NUSSELT_NUMBER_AGR(voidage, temp_ratio, sub_ID, Re, direction)

    use numvar
    use field
    use mesh
    use parall

    implicit none

    double precision :: NUSSELT_NUMBER_AGR
!! Aguments
    double precision :: voidage, temp_ratio, sub_ID, Re
    integer :: direction

!! Local variables 
    double precision, dimension(:), pointer :: A_vis, A_cond, A_cp
    double precision :: St0, St, Pr, Nusselt
    double precision :: index1 = 0.71, index2 = 0.5

    integer :: isub

    if (voidage<0.7d0 .or. voidage>=0.8d0)then
      if (irangp.le.0) print*,'Warning: voidage is out of range. Check the code'
    end if

!! EDF report: E Romero, TET/REP/0054/94 CEN/DFWG/P(92)4743/B/3, 1994
!! voidage = 0.7
    if(Re<100.d0)then
      if (irangp.le.0) print*, "Warning: st calculation out of Re lower bound"
      st0 = 17.6d0/Re
    else if(Re>=100.d0 .and. Re<2000.d0)then
      st0 = 17.6d0/Re
    else if(Re>=2000.d0 .and. Re<5000.d0)then
      st0 = 0.02664d0*(Re**(-0.14573d0)) 
    else if(Re>=5000.d0 .and. Re<1.d4)then
      st0 = 0.01064d0*(Re**(-0.03797d0))
    else if(Re>=1.d4 .and. Re<2.d4)then
      st0 = 0.001281d0*(Re**0.19188d0)
    else if(Re>=2.d4 .and. Re<5.d4)then
      st0 = 0.0101d0*(Re**(-0.01662d0))
    else if(Re>=5.d4 .and. Re<1.d6)then
      st0 = 0.0834d0*(Re**(-0.21179d0))       
    else
      if (irangp.le.0) print*, "Warning: cf calculation out of Re upper bound"
      st0 = 0.0834d0*(Re**(-0.21179d0)) 
    end if

    isub = int(sub_ID)

    select case (direction)
    case(1)
      st = st0/(temp_ratio**(index1*index2))
    case(2)

    case(3)

    end select

!!  calculate the Prandtl number
    call field_get_val_s_by_name('AverageDynamicViscosity', A_vis)
    call field_get_val_s_by_name('AverageThermalConductivity', A_cond)
    call field_get_val_s_by_name('AverageHeatCapacity', A_cp) 
    Pr = A_vis(isub)*A_cp(isub)/A_cond(isub) !(not correct, check!)

!!  convert the Stanton number into Nusselt number
    Nusselt = st*Re*Pr
    NUSSELT_NUMBER_AGR = Nusselt

  end function NUSSELT_NUMBER_AGR


!===============================================================================
!! This is for HTGR (GA design)
  function NUSSELT_NUMBER_HTGR_GA(sub_ID, Re, Pr, heating)

    use numvar
    use field
    use mesh
    use parall

    implicit none

    double precision :: NUSSELT_NUMBER_HTGR_GA
!! Aguments
    double precision :: sub_ID, Re, Pr
    logical :: heating

!! Local variables 
    double precision :: Nusselt
    logical :: condition1, condition2, condition3, condition
    integer :: isub

    isub = int(sub_ID)

    if(Re<=2300.d0)then
      Nusselt = 4.36d0
!! Dittus-Boelter correlation (developed flow, L/D>60)
!  For L/D<60, entrance effect should be considered.
!  0.023 based on Nuclear Systems I
!  0.021 based on Tak et al 2012
    else
      if (heating)then
  !      Nusselt = 0.023d0*(Re**0.8d0)*(Pr**0.4d0)
        Nusselt = 0.021d0*(Re**0.8d0)*(Pr**0.4d0)
      else
        Nusselt = 0.023d0*(Re**0.8d0)*(Pr**0.3d0)
      end if
    end if

!! It is better the correlation is used within the range of the below conditions
    condition1 = Re>=1.d4
    condition2 = Pr>=0.7d0 .and. Pr<=1.d2
    condition = condition1 .and. condition2

    if(.not.condition) then
      !print*, 'Warning! Nusselt number calculation error: out of range at ',&
      !        'Subchannel ID ', isub
    end if
    NUSSELT_NUMBER_HTGR_GA = Nusselt

  end function NUSSELT_NUMBER_HTGR_GA


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

    if(.not.allocated(user_prop)) allocate(user_prop(n_rows,col))

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
    double precision, dimension(:,:), pointer :: A_vel, B_vel, vel
    double precision, dimension(:), pointer :: den
    double precision, allocatable, dimension(:,:) :: subchannel_vel
    double precision, allocatable, dimension(:,:) :: subchannel_bulk_vel
    double precision, allocatable, dimension(:) :: subchannel_rho_A

    call field_get_val_v(ivarfl(iu), vel)
    call field_get_val_v_by_name('AverageVelocity', A_vel)
    call field_get_val_v_by_name('BulkVelocity', B_vel)
    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s(icrom, den)

    allocate(subchannel_vel(n_cell,3))
    allocate(subchannel_bulk_vel(n_cell,3))
    allocate(subchannel_rho_A(n_cell))

    subchannel_vel = 0.d0
    subchannel_bulk_vel = 0.d0
    subchannel_rho_A = 0.d0
    do iel=1, ncel
      j = int(subchannel_ID(iel))
      if(j>0)then
        do ii=1,3
          subchannel_vel(j,ii) = subchannel_vel(j,ii)+vel(ii,iel)*volume(iel)
          subchannel_bulk_vel(j,ii) = subchannel_bulk_vel(j,ii) &
                                    + vel(ii,iel)*den(iel)*volume(iel)
        end do
        subchannel_rho_A(j) = subchannel_rho_A(j)+den(iel)*volume(iel)
      end if
    end do

    if (irangp.ge.0) then
      do ii=1,3
        call parrsm(n_cell, subchannel_vel(:,ii))
        call parrsm(n_cell, subchannel_bulk_vel(:,ii))
      end do
      call parrsm(n_cell, subchannel_rho_A)
    endif

    do i=1,n_cell
      do ii=1,3
        subchannel_vel(i,ii) = subchannel_vel(i,ii)/subchannel_vol(i)
        subchannel_bulk_vel(i,ii) = subchannel_bulk_vel(i,ii) &
                                  / subchannel_rho_A(i)
      end do
    end do

    do iel=1, ncel
      j = int(subchannel_ID(iel))
      if(j>0)then
        do ii=1,3
          A_vel(ii,iel) = subchannel_vel(j,ii)
          B_vel(ii,iel) = subchannel_bulk_vel(j,ii)
        end do
      end if
    end do

!! One optimization has been done in this and the following subroutine
!! to reduce simulatons time. To locate modification, refer to the version
!! without optimization. 
   
!    if (ntcabs==100.and.irangp<=0)then
!      open (10, file='optimized.txt')
!      do j=1,n_cell
!        write(10,*)subchannel_vel(j,1),subchannel_vel(j,2),subchannel_vel(j,3)
!      end do
!      write(10,*)'The end of writing'
!    end if

    deallocate(subchannel_vel)
    deallocate(subchannel_bulk_vel)
    deallocate(subchannel_rho_A)

  end subroutine subchannel_velocity


!===============================================================================
 !! This routine is used to initialize subchannel tempeatures
  subroutine subchannel_temperatures_init()

    use numvar
    use optcal
    use field
    use mesh
    use parall
    use cstphy

    implicit none

    integer i,j,iel
    double precision, dimension(:), pointer :: subchannel_ID
    double precision, dimension(:), pointer :: A_temp, A_Tbulk, A_twall, A_tflux
    double precision, dimension(:), pointer :: tlocal, elocal
    double precision, allocatable, dimension(:) :: subchannel_temp

    if (iscalt.le.0) return

    call field_get_val_s_by_name('SubchannelID', subchannel_ID) 
    call field_get_val_s_by_name('AverageTemperature', A_temp)
    call field_get_val_s_by_name('SubchannelBulkTemperature', A_Tbulk)
    call field_get_val_s_by_name('SubchannelWallTemperature', A_twall)
    call field_get_val_s_by_name('SubchannelWallHeatFlux', A_tflux)
    if (itherm.ne.1) then
      call field_get_val_s_by_name('temperature', tlocal)
      tlocal = t0-273.15d0
    else
      call field_get_val_s(ivarfl(isca(iscalt)), tlocal)
    end if

    allocate(subchannel_temp(n_cell))

    subchannel_temp = 0.d0
    do iel=1, ncel
      j = int(subchannel_ID(iel))
      if(j>0) &
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
      if(j>0)then
        A_temp(iel) = subchannel_temp(j)
        A_Tbulk(iel) = A_temp(iel)
        A_twall(iel) = A_temp(iel)+5.d0
        A_tflux(iel) = 0.d0
      end if
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
    double precision, dimension(:), pointer :: A_temp, A_Tbulk, temp, enth, den
    double precision, dimension(:,:), pointer :: vel
    double precision, allocatable, dimension(:) :: subchannel_temp, &
                                                   subchannel_Tbulk, &
                                                   subchannel_mass_flux

    if (iscalt.gt.0) then
      ivart = isca(iscalt)

      if (itherm.ne.1) then
        call field_get_val_s_by_name('temperature', temp)
        call field_get_val_s(ivarfl(ivart), enth)
        do iel=1,ncel
          !call usthht (1, enth(iel), temp(iel)) !TODO
        end do
      else
        call field_get_val_s(ivarfl(ivart), temp)
      end if

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
      if(j>0)then
        subchannel_temp(j) = subchannel_temp(j)+temp(iel)*volume(iel)
        subchannel_Tbulk(j) = subchannel_Tbulk(j)+temp(iel)*volume(iel) &
                            * den(iel)*abs(vel(3,iel))
        subchannel_mass_flux(j) = subchannel_mass_flux(j)+volume(iel) &
                                * den(iel)*abs(vel(3,iel))
      end if
    end do

    if (irangp.ge.0) then
      call parrsm(n_cell, subchannel_temp)
      call parrsm(n_cell, subchannel_Tbulk)
      call parrsm(n_cell, subchannel_mass_flux)
    endif

    do i=1,n_cell
      subchannel_temp(i) = subchannel_temp(i)/subchannel_vol(i)
      if(abs(subchannel_mass_flux(i))<1.d-12)then
        subchannel_Tbulk(i) = subchannel_temp(i)
      else
        subchannel_Tbulk(i) = subchannel_Tbulk(i)/subchannel_mass_flux(i)
      end if
    end do

    do iel=1, ncel
      j = int(subchannel_ID(iel))
      if(j>0)then
        A_temp(iel) = subchannel_temp(j)
        A_Tbulk(iel) = subchannel_Tbulk(j)
      end if
    end do

!! Output for checking results
!    open(200,file='test.txt')
!    if(ntcabs==500.and.irangp<=0)then
!      do i=1, n_cell
!        write(200,*) i,subchannel_mass_flux(i),subchannel_Tbulk(i), &
!                     subchannel_temp(i)
!      end do
!    end if
!    close(200)

    deallocate(subchannel_temp, subchannel_Tbulk, subchannel_mass_flux)
 
  end subroutine subchannel_temperature


!===============================================================================
!! This routine is used to calculate the wall tempeatures
  subroutine subchannel_Twall()

    use numvar
    use optcal
    use field
    use mesh
    use parall
    use cs_c_bindings

    implicit none

    character(len=80) :: chaine
    integer :: i, j, iel, ielt, ifac, nlelt
    integer, allocatable, dimension(:) :: lstelt
    double precision, dimension(:), pointer :: subchannel_ID, A_twall, A_tflux,&
                                               Nusselt, tbulk, tlocal, Dh,     &
                                               cond, htc
    double precision, dimension(:), pointer :: twall0, twall1, tfluxwall
    double precision, dimension(:), allocatable :: subchannel_wall_Temp,       &
                                                   subchannel_wall_heat_flux,  &
                                                   subchannel_wall_area
    integer, allocatable, dimension(:) :: nf, nc

    if (iscalt.le.0) return

!! Check number of wall-adjacent cells in each suchannel cell  
    allocate(nf(n_cell), nc(n_cell))
    nf = 0
    nc = 0

    if (iscalt.le.0) return

!! Get all wall faces (only the SubChCFD domain)
    allocate(lstelt(nfabor))
    chaine = all_walls
    call getfbr(chaine, nlelt, lstelt)

!! Field arrays used
    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('HydraulicDiameter', Dh)
    call field_get_val_s_by_name('SubchannelWallTemperature', A_twall)
    call field_get_val_s_by_name('SubchannelWallHeatFlux', A_tflux)
    call field_get_val_s_by_name('NusseltNumber', Nusselt)
    call field_get_val_s_by_name('HTC', htc)
    call field_get_val_s_by_name('SubchannelBulkTemperature', tbulk)
    call field_get_val_s_by_name('WallTemperature0', twall0)
    call field_get_val_s_by_name('WallTemperature1', twall1)
    call field_get_val_s_by_name('WallThermalFlux', tfluxwall)
    call field_get_val_s_by_name('AverageThermalConductivity', cond)
    if (itherm.ne.1) then
      call field_get_val_s_by_name('temperature', tlocal)
    else
      call field_get_val_s(ivarfl(isca(iscalt)), tlocal)
    end if

    allocate(subchannel_wall_area(n_cell))
    allocate(subchannel_wall_Temp(n_cell))
    allocate(subchannel_wall_heat_flux(n_cell))

!! Disable solid cells in case of CHT
!    if (fluid_solid) call cs_porous_model_set_has_disable_flag(1)

!! Calculate the wall thermal fluxes of each boundary face
    call wall_thermal_flux (nlelt, lstelt, tfluxwall)

    subchannel_wall_Temp = 0.d0
    subchannel_wall_heat_flux = 0.d0
    subchannel_wall_area = 0.d0
    do ielt=1, nlelt
      ifac = lstelt(ielt)
      iel = ifabor(ifac)
      j = int(subchannel_ID(iel))
!      if(cell_is_active(iel)==0)cycle
      if(j<0)cycle
      htc(iel) = Nusselt(iel)*cond(iel)/Dh(iel)
      twall0(ifac) = tbulk(iel)-tfluxwall(ifac)/htc(iel)
      twall1(ifac) = tlocal(iel)-tfluxwall(ifac)/htc(iel)
      subchannel_wall_area(j) = subchannel_wall_area(j)+surfbn(ifac)
!!    Note here the average wall temperature is calculated based on twall1
      subchannel_wall_Temp(j) = subchannel_wall_Temp(j)+twall1(ifac) &
                              * surfbn(ifac)
      subchannel_wall_heat_flux(j) = subchannel_wall_heat_flux(j)    &
                                   + tfluxwall(ifac)*surfbn(ifac)
      nf(j) = nf(j)+1
    end do

!! Re-enable solid cells in case of CHT
!    if(fluid_solid) call cs_porous_model_set_has_disable_flag(0)

    if (irangp.ge.0) then
      call parrsm(n_cell, subchannel_wall_area)
      call parrsm(n_cell, subchannel_wall_Temp)
      call parrsm(n_cell, subchannel_wall_heat_flux)
      call parism(n_cell, nf)
    endif

    do i=1,n_cell
      subchannel_wall_Temp(i) = subchannel_wall_Temp(i)/subchannel_wall_area(i)
      subchannel_wall_heat_flux(i) = subchannel_wall_heat_flux(i)  &
                                   / subchannel_wall_area(i)
    end do

    do iel=1, ncel
      j = int(subchannel_ID(iel))
      if(j>0) then
        A_twall(iel) = subchannel_wall_Temp(j)
        A_tflux(iel) = subchannel_wall_heat_flux(j)
      end if
    end do

!! Output for checking results
!    open(200,file='test.txt')
!    if(ntcabs==500.and.irangp<=0)then
!      do i=1, n_cell
!        write(200,*) i,subchannel_wall_area(i),subchannel_wall_Temp(i), nf(i)
!      end do
!    end if
!    close(200)
                
!! Deallocate temporary arrays
    deallocate(lstelt)
    deallocate(subchannel_wall_area)
    deallocate(subchannel_wall_Temp)
    deallocate(subchannel_wall_heat_flux)
    deallocate(nc,nf)

  end subroutine subchannel_Twall


!===============================================================================
!! This routine is used to calculate the wall thermal flux
  subroutine wall_thermal_flux(nfbrps,lstfbr,bflux)

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
    integer, intent(in) :: nfbrps
    integer, dimension(:), intent(in) :: lstfbr
    double precision, dimension(:), intent(inout) :: bflux

    ! Local variables
    integer ::         inc, iprev
    integer ::         ifac, iloc, ivar
    integer ::         iel
    integer ::         iflmab

    double precision :: cpp, srfbn
    double precision :: flumab, diipbx, diipby, diipbz, tcel

    double precision, dimension(:), pointer :: cpro_cp

    double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
    double precision, allocatable, dimension(:,:) :: grad
    double precision, dimension(:), pointer :: bmasfl, tscalp

    double precision, dimension(:), pointer :: cht

    type(var_cal_opt) :: vcopt

    call field_get_val_s_by_name('CHTInternalCoupling', cht)

! Initialize variables to avoid compiler warnings
    if (iscalt.gt.0) then
      ivar = isca(iscalt)

! Boundary condition pointers for gradients and advection
      call field_get_coefa_s(ivarfl(ivar), coefap)
      call field_get_coefb_s(ivarfl(ivar), coefbp)

! Boundary condition pointers for diffusion
      call field_get_coefaf_s(ivarfl(ivar), cofafp)
      call field_get_coefbf_s(ivarfl(ivar), cofbfp)

! Pointers to fields and properties
      call field_get_val_s(ivarfl(ivar), tscalp)

      if (itherm.eq.1 .and. icp.ge.0) then
        call field_get_val_s(icp, cpro_cp)
      endif

      call field_get_key_int(ivarfl(ivar), kbmasf, iflmab)
      call field_get_val_s(iflmab, bmasfl)

! Compute variable values at boundary faces
      call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

! Reconstructed fluxes
      if (vcopt%ircflu .gt. 0 .and. itbrrb.eq.1) then
! Compute gradient of temperature / enthalpy
        allocate(grad(3,ncelet))

        inc = 1
        iprev = 1

        call field_gradient_scalar(ivarfl(ivar), iprev, inc, grad)

! Compute diffusive and convective flux using reconstructed temperature
        do iloc = 1, nfbrps

          ifac = lstfbr(iloc)
          iel = ifabor(ifac)

          if(cht(ifac)>0.5d0) cycle

          diipbx = diipb(1,ifac)
          diipby = diipb(2,ifac)
          diipbz = diipb(3,ifac)

          tcel = tscalp(iel)&
               + diipbx*grad(1,iel) + diipby*grad(2,iel) + diipbz*grad(3,iel)

          if (itherm.eq.1) then
            if (icp.ge.0) then
              cpp = cpro_cp(iel)
            else
              cpp = cp0
            endif
          else
            cpp = 1.d0
          endif

          srfbn = max(surfbn(ifac), epzero**2)
          flumab = bmasfl(ifac)

          bflux(ifac) = (cofafp(ifac) + cofbfp(ifac)*tcel)&
                      - flumab*cpp/srfbn*(coefap(ifac) + coefbp(ifac)*tcel)
        enddo

        deallocate(grad)

      else ! If flux is not reconstructed

! Compute diffusive and convective flux using non-reconstructed temperature

        do iloc = 1, nfbrps

          ifac = lstfbr(iloc)
          iel = ifabor(ifac)
          tcel = tscalp(iel)

          if(cht(ifac)>0.5d0) cycle

          if (itherm.eq.1) then
            if (icp.ge.0) then
              cpp = cpro_cp(iel)
            else
              cpp = cp0
            endif
          else
            cpp = 1.d0
          endif

          srfbn = max(surfbn(ifac), epzero**2)
          flumab = bmasfl(ifac)

          bflux(ifac) = (cofafp(ifac) + cofbfp(ifac)*tcel)&
                      - flumab*cpp/srfbn*(coefap(ifac) + coefbp(ifac)*tcel)

        enddo

      endif ! test on reconstruction

    else ! if thermal variable is not available

      do iloc = 1, nfbrps
        ifac = lstfbr(iloc)
        if(cht(ifac)>0.5d0) cycle
        bflux(iloc) = 0.d0
      enddo

    endif

    return
    
  end subroutine wall_thermal_flux


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
      if(j>0)then
        subchannel_den(j) = subchannel_den(j)+den(iel)*volume(iel)
        subchannel_vis(j) = subchannel_vis(j)+dynamicvis(iel)*volume(iel)
      end if
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
      if(j>0)then
        A_den(iel) = subchannel_den(j)
        A_vis(iel) = subchannel_vis(j)
      end if
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
      if(j>0)then
        if (itherm.ne.1) then
          if (ifcvsl.ge.0) then
            if (icp.gt.0) then
              subchannel_cond(j) = subchannel_cond(j)+cond(iel)*cpro_cp(iel) &
                                 * volume(iel)
            else
              subchannel_cond(j) = subchannel_cond(j)+cond(iel)*cp0*volume(iel)
            end if
          else
            subchannel_cond(j) = subchannel_cond(j)+visls0(iscalt)*volume(iel)
          end if
        else
          if (ifcvsl.ge.0) then
            subchannel_cond(j) = subchannel_cond(j)+cond(iel)*volume(iel)
          else
            subchannel_cond(j) = subchannel_cond(j)+visls0(iscalt)*volume(iel)
          end if
        end if
        
        if (icp.gt.0) then
          subchannel_cp(j) = subchannel_cp(j)+cpro_cp(iel)*volume(iel)
        else
          subchannel_cp(j) = subchannel_cp(j)+cp0*volume(iel)
        end if
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
      if(j>0)then
        A_cond(iel) = subchannel_cond(j)
        A_cp(iel) = subchannel_cp(j)
      end if
    end do

    deallocate(subchannel_cond, subchannel_cp)


  end subroutine subchannel_physical_properties


!===============================================================================
!! This routine is used to calculate the filtered velocity of each filter block
!! in the porous media zone 
  subroutine filtered_velocity()

    use numvar
    use field
    use mesh
    use parall

    implicit none

    integer :: i, ii, iel, j
    double precision, dimension(:), pointer :: filter_ID
    double precision, dimension(:,:), pointer :: F_vel, vel
    double precision, allocatable, dimension(:,:) :: filter_vel
    double precision, allocatable, dimension(:) :: filter_vol
    integer :: n_node, n_entity, n_vertices

    call field_get_val_v(ivarfl(iu), vel)
    call field_get_val_v_by_name('FilteredVelocity', F_vel)
    call field_get_val_s_by_name('FilterID', filter_ID) 

    call mesh_dimensions(filtering_mesh, n_node, n_entity)

    allocate(filter_vel(n_entity,3))
    allocate(filter_vol(n_entity))

    filter_vel = 0.d0
    do iel=1, ncel
      j = int(filter_ID(iel))
      if(j>0)then
        do ii=1,3
          filter_vel(j,ii) = filter_vel(j,ii)+vel(ii,iel)*volume(iel)
        end do
      end if
    end do

    if (irangp.ge.0) then
      call parrsm(n_entity, filter_vel(:,1))
      call parrsm(n_entity, filter_vel(:,2))
      call parrsm(n_entity, filter_vel(:,3))
    endif

    do i=1,n_entity
      filter_vol(i) = filter_cell_vol(1,i)
      filter_vel(i,1) = filter_vel(i,1)/filter_vol(i)
      filter_vel(i,2) = filter_vel(i,2)/filter_vol(i)
      filter_vel(i,3) = filter_vel(i,3)/filter_vol(i)
    end do

    do iel=1, ncel
      j = int(filter_ID(iel))
      if(j>0)then
        F_vel(1,iel) = filter_vel(j,1)
        F_vel(2,iel) = filter_vel(j,2)
        F_vel(3,iel) = filter_vel(j,3)
      end if
    end do

    deallocate(filter_vel)
    deallocate(filter_vol)

  end subroutine filtered_velocity


!===============================================================================
!! This routine is used to calculate the filtered physical properties for the 
!! porous media zones   
  subroutine filtered_physical_properties()

    use numvar
    use optcal
    use cstphy
    use field
    use mesh
    use parall

    implicit none

    integer :: i, ii, iel, j
    integer :: ifcvsl
    double precision, dimension(:), pointer :: filter_ID, F_den, F_vis, &
                                               F_cond, F_cp
    double precision, dimension(:), pointer :: den, dynamicvis, cond, cpro_cp
    double precision, allocatable, dimension(:) :: filter_den, &
                                                   filter_vis, &
                                                   filter_cond, &
                                                   filter_cp 
    double precision, allocatable, dimension(:) :: filter_vol
    integer :: n_node, n_entity

    call field_get_val_s(iprpfl(iviscl), dynamicvis)
    call field_get_val_s(icrom, den)

    call field_get_val_s_by_name('FilterID', filter_ID)
    call field_get_val_s_by_name('FilteredDensity', F_den)
    call field_get_val_s_by_name('FilteredDynamicViscosity', F_vis)

    call mesh_dimensions(filtering_mesh, n_node, n_entity)
    
    allocate(filter_vol(n_entity))
    allocate(filter_den(n_entity), filter_vis(n_entity))

    filter_den = 0.d0
    filter_vis = 0.d0
    do iel=1, ncel
      j = int(filter_ID(iel))
      if(j>0)then
        filter_den(j) = filter_den(j)+den(iel)*volume(iel)
        filter_vis(j) = filter_vis(j)+dynamicvis(iel)*volume(iel)
      end if
    end do

    if (irangp.ge.0) then
      call parrsm(n_entity, filter_den)
      call parrsm(n_entity, filter_vis)
    endif

    do i=1,n_entity
      filter_vol(i) = filter_cell_vol(1,i)
      filter_den(i) = filter_den(i)/filter_vol(i)
      filter_vis(i) = filter_vis(i)/filter_vol(i)
    end do

    do iel=1, ncel
      j = int(filter_ID(iel))
      if(j>0)then
        F_den(iel) = filter_den(j)
        F_vis(iel) = filter_vis(j)
      end if
    end do

    deallocate(filter_den, filter_vis)
    deallocate(filter_vol)

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

    call field_get_val_s_by_name('FilteredThermalConductivity', F_cond)
    call field_get_val_s_by_name('FilteredHeatCapacity', F_cp)

    allocate(filter_vol(n_entity))
    allocate(filter_cond(n_entity), filter_cp(n_entity))

    filter_cond = 0.d0
    filter_cp = 0.d0
    do iel=1, ncel
      j = int(filter_ID(iel))
      if(j>0)then

        if (itherm.ne.1) then
          if (ifcvsl.ge.0) then
            if (icp.gt.0) then
              filter_cond(j) = filter_cond(j)+cond(iel)*cpro_cp(iel)*volume(iel)
            else
              filter_cond(j) = filter_cond(j)+cond(iel)*cp0*volume(iel)
            end if
          else
            filter_cond(j) = filter_cond(j)+visls0(iscalt)*volume(iel)
          end if
        else
          if (ifcvsl.ge.0) then
            filter_cond(j) = filter_cond(j)+cond(iel)*volume(iel)
          else
            filter_cond(j) = filter_cond(j)+visls0(iscalt)*volume(iel)
          end if
        end if

        if (icp.gt.0) then
          filter_cp(j) = filter_cp(j)+cpro_cp(iel)*volume(iel)
        else
          filter_cp(j) = filter_cp(j)+cp0*volume(iel)
        end if
      end if
    end do

    if (irangp.ge.0) then
      call parrsm(n_entity, filter_cond)
      call parrsm(n_entity, filter_cp)
    endif

    do i=1,n_entity
      filter_vol(i) = filter_cell_vol(1,i)
      filter_cond(i) = filter_cond(i)/filter_vol(i)
      filter_cp(i) = filter_cp(i)/filter_vol(i)
    end do

    do iel=1, ncel
      j = int(filter_ID(iel))
      if(j>0)then
        F_cond(iel) = filter_cond(j)
        F_cp(iel) = filter_cp(j)
      end if
    end do

    deallocate(filter_cond, filter_cp)
    deallocate(filter_vol)

  end subroutine filtered_physical_properties


!===============================================================================
!! This routine is used to obtain xyz coordinate ranges of the coupling region
  subroutine coupling_xyz(numcpl, xyzcenmin, xyzcenmax, xyzfbmin, xyzfbmax) 

    use paramx
    use numvar
    use entsor
    use optcal
    use cstphy
    use cstnum
    use parall
    use cplsat
    use mesh

    implicit none

!!  Arguments
    integer, intent(IN) :: numcpl
    double precision, dimension(3), intent(OUT) :: xyzcenmin, xyzcenmax, &
                                                   xyzfbmin, xyzfbmax
                                                   
!!  Local variables
    integer :: i, ipt, ielloc
    integer :: ncesup, nfbsup
    integer :: ncecpl, nfbcpl, ncencp, nfbncp
    integer :: ncedis, nfbdis
    integer :: ityloc, ityvar

    integer, allocatable, dimension(:) :: lcecpl, lfbcpl, lcencp, lfbncp
    integer, allocatable, dimension(:) :: locpts

    double precision, allocatable, dimension(:,:) :: coopts , djppts , dofpts
    double precision, allocatable, dimension(:) :: pndpts

    double precision, allocatable, dimension(:,:) :: xyzcenloc

    call nbecpl(numcpl, ncesup, nfbsup, ncecpl, nfbcpl, ncencp, nfbncp)

    allocate(lcecpl(ncecpl))
    allocate(xyzcenloc(3,ncecpl))

    call lelcpl(numcpl, ncecpl, nfbcpl, lcecpl, lfbcpl)

    do ipt=1, ncecpl
      ielloc = lcecpl(ipt)
      xyzcenloc(:,ipt) = xyzcen(:,ielloc)
    end do

    do i=1, 3
      xyzcenmin(i) = minval(xyzcenloc(i,:))
      xyzcenmax(i) = maxval(xyzcenloc(i,:))
      if (irangp.ge.0) then
        call parmin(xyzcenmin(i))
        call parmax(xyzcenmax(i))
      end if
    end do

!    write(*,*) 'local cell min', xyzcenmin
!    write(*,*) 'local cell max', xyzcenmax

    ityvar = 2

    call npdcpl(numcpl, ncedis, nfbdis)

    ! Allocate temporary arrays for geometric quantities
    allocate(locpts(nfbdis))
    allocate(coopts(3,nfbdis), djppts(3,nfbdis), dofpts(3,nfbdis))
    allocate(pndpts(nfbdis))

    call coocpl(numcpl, nfbdis, ityvar, ityloc, locpts, coopts, &
                djppts, dofpts, pndpts)

    if (ityloc.eq.2) then
      write(nfecra,1000)
      call csexit(1)
    endif

    do i=1, 3
      xyzfbmin(i) = minval(coopts(i,:))
      xyzfbmax(i) = maxval(coopts(i,:))
      if (irangp.ge.0) then
        call parmin(xyzfbmin(i))
        call parmax(xyzfbmax(i))
      end if
    end do

!    write(*,*) 'distant face min', xyzfbmin
!    write(*,*) 'distant face max', xyzfbmax

    deallocate(lcecpl)
    deallocate(xyzcenloc)

    !--------
    ! Formats
    !--------
    1000 format(                                                     &
    '@                                                            ',/,&
    '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
    '@                                                            ',/,&
    '@ @@ ATTENTION :                                             ',/,&
    '@    =========                                               ',/,&
    '@    The couling via faces as supported elements is          ',/,&
    '@    not yet mamaged by the code.                            ',/,&
    '@                                                            ',/,&
    '@    Calculation cannot be performed.                        ',/,&
    '@                                                            ',/,&
    '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
    '@                                                            ',/)
    !----
    ! End
    !----

  end subroutine coupling_xyz


!===============================================================================
!! This routine is used to solve a steady poisson equation with a given RHS
  subroutine poisson &
   (fext   , dfext  ,                                              &
    phydr  , flumas , flumab ,                                     &
    coefap , coefbp ,                                              &
    cofafp , cofbfp ,                                              &
    viscf  , viscb  ,                                              &
    dam    , xam    ,                                              &
    dpvar  , smbr   )

    use atincl, only: iatmst
    use paramx
    use numvar
    use entsor
    use cstnum
    use optcal
    use parall
    use period
    use mesh
    use cs_c_bindings

    implicit none

    ! Arguments

    double precision fext(3,ncelet)
    double precision dfext(3,ncelet)
    double precision phydr(ncelet)
    double precision flumas(nfac), flumab(nfabor)
    double precision coefap(nfabor), coefbp(nfabor)
    double precision cofafp(nfabor), cofbfp(nfabor)
    double precision viscf(nfac), viscb(nfabor)
    double precision dam(ncelet), xam(nfac)
    double precision dpvar(ncelet)
    double precision smbr(ncelet)

    ! Local variables

    character(len=80) :: chaine
    integer          lchain
    integer          f_id,   inc   , init  , isym
    integer          iel   , ical
    integer          nswmpr
    integer          isweep, niterf
    integer          iphydp
    integer          nswrgp, imligp, iwgrp, iwarnp
    integer          idiffp, iconvp, ndircp
    integer          ibsize, iesize
    integer          imucpp, f_id0

    double precision residu, rnorm , rnrmf , rnrmdf
    double precision epsrgp, climgp, extrap, epsilp
    double precision precre, precab, thetap

    double precision rvoid(1)

    double precision, allocatable, dimension(:) :: rovsdt, div_dfext, viscce

    type(var_cal_opt) :: vcopt_u, vcopt_pr

    ! Allocate temporary arrays
    allocate(rovsdt(ncelet), div_dfext(ncelet), viscce(ncelet))

    call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
    call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_pr)

    f_id = -1
    chaine = 'hydrostatic_p'
    lchain = 16

    ! Symetric
    isym  = 1

    ! Matrix block size
    ibsize = 1
    iesize = 1

    f_id0 = -1

    do iel = 1, ncel
      rovsdt(iel) = 0.d0
    enddo

    do iel = 1, ncel
      viscce(iel) = 1.d0
    enddo

    call viscfa(imvisf, viscce, viscf, viscb)

    iconvp = 0
    idiffp = 1
    ndircp = 0
    thetap = 1.d0
    imucpp = 0

    call matrix (iconvp, idiffp, ndircp, isym, thetap, imucpp, coefbp, cofbfp, &
                 rovsdt, flumas, flumab, viscf, viscb, rvoid, dam, xam)

    init   = 1
    inc    = 0
    nswrgp = vcopt_pr%nswrgr
    imligp = vcopt_pr%imligr
    iwarnp = vcopt_pr%iwarni
    epsrgp = vcopt_pr%epsrgr
    climgp = vcopt_pr%climgr

    call projts (init, nswrgp, dfext, cofbfp, flumas, flumab, viscf, viscb, &
                 viscce, viscce, viscce)

    init = 1
    call divmas(init,flumas,flumab,div_dfext)
    rnorm = sqrt(cs_gdot(ncel,div_dfext,div_dfext))

    nswmpr = vcopt_pr%nswrsm

    do iel = 1,ncel
      phydr(iel) = 0.d0
      dpvar(iel) = 0.d0
      smbr(iel) = 0.d0
    enddo

    do isweep = 1, nswmpr
      do iel = 1, ncel
        smbr(iel) = - div_dfext(iel) - smbr(iel)
      enddo

      residu = sqrt(cs_gdot(ncel,smbr,smbr))

      if (residu.le.10.d0*vcopt_pr%epsrsm*rnorm) then
        goto 101
      endif

      do iel = 1, ncel
        dpvar(iel) = 0.d0
      enddo

      iwarnp = vcopt_pr%iwarni
      epsilp = vcopt_pr%epsilo
      ibsize = 1
      iesize = 1

      call sles_solve_native(f_id, chaine, isym, ibsize, iesize, dam, xam, &
                             epsilp, rnorm, niterf, residu, smbr, dpvar)

      if (isweep.eq.nswmpr) then
        do iel = 1, ncel
          phydr(iel) = phydr(iel) + dpvar(iel)
        enddo
      else
        do iel = 1, ncel
          phydr(iel) = phydr(iel) + dpvar(iel)
        enddo

        init = 1
        inc  = 1
        nswrgp = vcopt_pr%nswrgr
        imligp = vcopt_pr%imligr
        iwgrp  = vcopt_pr%iwgrec
        iwarnp = vcopt_pr%iwarni
        epsrgp = vcopt_pr%epsrgr
        climgp = vcopt_pr%climgr
        extrap = 0.d0
        iphydp = 1

        call itrgrp (f_id0, init, inc, imrgra, nswrgp, imligp, iphydp, iwgrp,  &
                     iwarnp, epsrgp, climgp, extrap, dfext, phydr, coefap, &
                     coefbp, cofafp, cofbfp, viscf, viscb, viscce, smbr)

      endif

    enddo

    101  continue

    deallocate(rovsdt, div_dfext, viscce)

    call sles_free_native(-1, chaine)

    return

  end subroutine poisson 


!===============================================================================
!! This routine is used to calculate the pressure field in case of coupling
  subroutine pressure_reconstruction()

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
    use cs_c_bindings

    implicit none

    integer :: iel, ifac
    integer :: iel1, irang1, irangv
    double precision hint, qimp
    double precision frcxt(3,ncelet), dfrcxt(3,ncelet)
    double precision viscf(nfac), viscb(ndimfb)
    double precision pressure_ref
    double precision, allocatable, dimension(:) :: dam, xam
    double precision, allocatable, dimension(:) :: dphi, rhs
    double precision, allocatable, dimension(:) :: coefap, coefbp
    double precision, allocatable, dimension(:) :: cofafp, cofbfp
    double precision, allocatable, dimension(:) :: iflux, bflux
    double precision, allocatable, dimension(:) :: res_phi

    double precision, dimension(:,:), pointer :: mom_hat, mom_pt
    double precision, dimension(:), pointer :: pressure_output, pressure_cal
    double precision, dimension(:,:), pointer :: grad_output
    integer :: f_id, imrgra_output, inc, nswrgp, imligp, iwarnp, iphydr_output,&
               iprev
    double precision :: epsrgp, climgp, extrap
    double precision, dimension(nfabor) :: coefap_output, coefbp_output

    type(var_cal_opt) :: vcopt_p, vcopt_u

    call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
    call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)

    allocate(dam(ncelet), xam(nfac))
    allocate(dphi(ncelet))
    allocate(rhs(ncelet))

    allocate(coefap(ndimfb), cofafp(ndimfb))
    allocate(coefbp(ndimfb), cofbfp(ndimfb))

    allocate(iflux(nfac), bflux(ndimfb))

    allocate(res_phi(ncelet))

    call field_get_val_v_by_name('MomSourceSHat', mom_hat)
    call field_get_val_v_by_name('MomSourcePotential', mom_pt)

    do ifac = 1, nfac
      iflux(ifac) = 0.d0
    enddo

    do ifac = 1, nfabor
      iel = ifabor(ifac)
      bflux(ifac) = 0.d0
      hint = 0.d0
      qimp = 0.d0

      call set_neumann_scalar (coefap(ifac), cofafp(ifac), coefbp(ifac), &
                               cofbfp(ifac), qimp, hint)
    enddo

    do iel = 1, ncel
      if(ipre_rec==1)then
        frcxt(1,iel)  = mom_hat(1,iel)
        frcxt(2,iel)  = mom_hat(2,iel)
        frcxt(3,iel)  = mom_hat(3,iel)
        dfrcxt(1,iel) = mom_hat(1,iel)
        dfrcxt(2,iel) = mom_hat(2,iel)
        dfrcxt(3,iel) = mom_hat(3,iel)
      elseif(ipre_rec==2) then
        frcxt(1,iel)  = mom_pt(1,iel)
        frcxt(2,iel)  = mom_pt(2,iel)
        frcxt(3,iel)  = mom_pt(3,iel)
        dfrcxt(1,iel) = mom_pt(1,iel)
        dfrcxt(2,iel) = mom_pt(2,iel)
        dfrcxt(3,iel) = mom_pt(3,iel)      
      end if
    enddo

    ! Parallelism and periodicity treatment
    if (irangp.ge.0.or.iperio.eq.1) then
      call synvin(frcxt)
      call synvin(dfrcxt)
    endif

    call poisson (frcxt, dfrcxt, res_phi, iflux, bflux, coefap, coefbp, &
                  cofafp, cofbfp, viscf, viscb, dam, xam, dphi, rhs)

    deallocate(coefap, cofafp)
    deallocate(coefbp, cofbfp)

    call field_get_id('PressureOutput', f_id)
    call field_get_val_s(f_id, pressure_output)
    !using name to access field
    !call field_get_val_s_by_name('PressureOutput', pressure_output)

    call field_get_val_s(ivarfl(ipr), pressure_cal)

    iel1   = -999
    irang1 = -999 
    call findpt(ncelet, ncel, xyzcen, xyz_ref(1), xyz_ref(2), xyz_ref(3),&
                iel, irangv)
    if ((iel.ne.iel1).or.(irangv.ne.irang1)) then
      iel1   = iel
      irang1 = irangv
      if (irangp.eq.irangv) pressure_ref = res_phi(iel)
      if (irangp.ge.0) then
        call parall_bcast_r(irangv, pressure_ref)
      endif
    endif

    do iel=1, ncel
      res_phi(iel) = res_phi(iel)-pressure_ref
    end do

    do iel=1, ncel
      pressure_output(iel) = pressure_cal(iel)-res_phi(iel)
    end do

    call field_get_val_v_by_name('PressureGradient', grad_output)

    f_id=-1
    imrgra_output=0
    inc=1
    nswrgp=2
    imligp=-1
    iwarnp=0
    epsrgp=1.d-6
    climgp=0.d0
    extrap=0.d0
    coefap_output=0.d0
    coefbp_output=1.d0

    iprev  = 0
    !Not work for user defined field properties !TODO check if this is still
    !the case in v8.0
    !call field_gradient_scalar(f_id, iprev, inc, grad_output)

    !Gradient calculation for normal scalar field !TODO note changes of this
    !subroutine between v6.0 and v8.0
    call gradient_s(f_id, imrgra, inc, nswrgp, imligp, iwarnp, epsrgp, climgp, &
                    pressure_output, coefap_output, coefbp_output, grad_output)

    !Gradient calculation for potential type scalar field !TODO The following
    !gradient calculation method should be used, however, it is no longer 
    !available in v8.0
    !iphydr_output = 1
    !call gradient_potential_s (f_id, imrgra_output, inc, iccocg, nswrgp, &
    !                           imligp, iphydr_output, iwarnp, epsrgp, climgp, &
    !                           extrap, dfrcxt, pressure_output, coefap_output,&
    !                           coefbp_output, grad_output)

  end subroutine pressure_reconstruction


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
    
    call mesh_dimensions(subchannel_mesh_name, n_node, n_entity)
    
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
    allocate(x_channel_min(n_cell),y_channel_min(n_cell),z_channel_min(n_cell))
    allocate(x_channel_max(n_cell),y_channel_max(n_cell),z_channel_max(n_cell))
        
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
!!  Stored for global use
      x_channel_min(i) = xx(i,1)
      y_channel_min(i) = yy(i,1)
      z_channel_min(i) = zz(i,1)
      x_channel_max(i) = xx(i,2)
      y_channel_max(i) = yy(i,2)
      z_channel_max(i) = zz(i,2)

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
    subchannel_ID = -1.d0
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
!!  Stored for global use
      x_channel_min(i) = xx(i,1)
      y_channel_min(i) = yy(i,1)
      z_channel_min(i) = zz(i,1)
      x_channel_max(i) = xx(i,2)
      y_channel_max(i) = yy(i,2)
      z_channel_max(i) = zz(i,2)
      
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
    subchannel_ID = -1.d0
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
      if(subchannel_ID(iel)>0.d0)then
        pitch2diameter(iel) = PbyD(int(subchannel_ID(iel)))
        subchannel_type(iel) = cell_type(int(subchannel_ID(iel)))
      end if
    end do 

    deallocate(x,y,z)
    deallocate(index_cell)

  end subroutine subchannel_mesh_user


!===============================================================================
!! This routine is used to define zones for feedback in two-way coupling
  subroutine feedback_zone_user()

    use numvar
    use field
    use mesh
    use parall

    integer :: iel
    double precision :: x1, y1, z1, x2, y2, z2
    double precision, pointer, dimension(:) :: feedback_zone

    call field_get_val_s_by_name('FeedbackFlag', feedback_zone)

    x1 = 0.d0
    y1 = 0.d0
    z1 = -1150.d-3
    x2 = 66.1d-3
    y2 = 66.1d-3
    z2 = -871.d-3

    feedback_zone = 0.d0

    do iel=1,ncel
    if(xyzcen(1,iel)>x1 .and. xyzcen(2,iel)>y1 .and. xyzcen(3,iel)>z1 .and. &
       xyzcen(1,iel)<x2 .and. xyzcen(2,iel)<y2 .and. xyzcen(3,iel)<z2) then
      feedback_zone(iel) = 1.d0
    end if
    end do

  end subroutine feedback_zone_user

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

!!  User example 2
!      FRICTIONAL_FACTOR_USER = friction_user_2(iel)

  end function FRICTIONAL_FACTOR_USER


!-------------------------------------------------------------------------------
!! User example 1
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

!-------------------------------------------------------------------------------
!! User example 2
!! This function is used for forced convection with property variation
  function friction_user_2(iel)

    use numvar
    use field
    use mesh
    use parall
    use optcal

    implicit none

    double precision :: friction_user_2

!! Aguments
    integer          :: iel

!! Local variables
    integer          :: icorr
    integer          :: i
    double precision :: cf, cf0
    double precision :: temp, resi
    double precision, dimension(:), pointer :: tbulk, A_twall, Dh
    double precision, dimension(:,:), pointer :: B_vel
    double precision :: rhob, visb, cpb, condb, Reb, Prb, rhow
    double precision :: intp, intp1, intp2
    integer :: cpro_index

!    call field_get_val_s_by_name('SubchannelBulkTemperature', tbulk)
    call field_get_val_s_by_name('AverageTemperature', tbulk)
    call field_get_val_s_by_name('SubchannelWallTemperature', A_twall)
    call field_get_val_s_by_name('HydraulicDiameter', Dh)
    call field_get_val_v_by_name('BulkVelocity', B_vel)

    !! Obtain subchannel bulk properties based on bulk temperature
    cpro_index = 1
    do while (tbulk(iel) >= user_prop(cpro_index,1))
      cpro_index = cpro_index + 1
    end do
    intp = user_prop(cpro_index,1)-user_prop(cpro_index-1,1)
    intp1 = (tbulk(iel)-user_prop(cpro_index-1,1))/intp
    intp2 = 1.d0-intp1

    !! Bulk density
    rhob = user_prop(cpro_index-1,3)*intp2+(user_prop(cpro_index,3))*intp1
    !! Bulk heat capacity
    cpb = user_prop(cpro_index-1,4)*intp2+(user_prop(cpro_index,4))*intp1
    !! Bulk viscosity
    visb = user_prop(cpro_index-1,5)*intp2+(user_prop(cpro_index,5))*intp1
    !! Bulk thermal conductivity
    condb = user_prop(cpro_index-1,6)*intp2+(user_prop(cpro_index,6))*intp1

    !! Calculate subchannel bulk Reynolds number
    Reb = rhob*abs(B_vel(3,iel))*Dh(iel)/visb
    !! Calculate subchannel bulk Prandtl number
    Prb = visb*cpb/condb

    !! Select correlation to be used for turbulent flow
    !! 1: Blasius and McAdams
    !! 2: Petukhov and Popov (1983)
    !! 3: Colebrook-White
    icorr = 2

    !! For fully developed flow in smooth circular channels
    if(Reb<=2.3d3)then
      cf = 64.d0/Reb
    end if

    if(Reb>2.3d3)then
      !! Blasius
      if(Reb<3.d4) then
        cf = 0.316d0*(Reb**(-0.25d0))
      !! McAdams
      else if(Reb>=3.d4 .and. Reb<1.d6) then
        cf = 0.184d0*(Reb**(-0.2d0))
      else
        cf = 0.184d0*(Reb**(-0.2d0))
      end if

      if(icorr==2) then
        cf = 1.d0/(1.82d0*log10(Reb/8.d0))**2
      end if

      if(icorr==3) then
        do i=1, 10000
          cf0 = cf
          temp = -2.*log10(2.51d0/Reb/sqrt(cf))
          cf = 1./temp**2
          resi = abs(cf-cf0)
          if (resi<1.d-5) exit
        end do
        if(resi>=1.d-5) write(*,*)"Warning: non-convergence for cf"
      end if

    end if

!! To introduce corrections for impacts of variable properties
!! Petukhov and Popov (1983)
    cpro_index = 1
    do while (A_twall(iel) >= user_prop(cpro_index,1))
      cpro_index = cpro_index + 1
    end do
    intp = user_prop(cpro_index,1)-user_prop(cpro_index-1,1)
    intp1 = (A_twall(iel)-user_prop(cpro_index-1,1))/intp
    intp2 = 1.d0-intp1

    !! Wall density
    rhow = user_prop(cpro_index-1,3)*intp2+(user_prop(cpro_index,3))*intp1

    cf = cf*((rhow/rhob)**0.4)

    friction_user_2 = cf

  end function friction_user_2


!===============================================================================
!! This routine is used to define user Nusselt number 
  function NUSSELT_NUMBER_USER(iel)

    use numvar
    use optcal
    use field
    use mesh
    use parall

    implicit none

    double precision :: NUSSELT_NUMBER_USER

!! Aguments
    integer :: iel

    if (iscalt.le.0) return

!!  User example 1
    if(.true.) &
      NUSSELT_NUMBER_USER = nusselt_user_1(iel)

!!  User example 2
!      NUSSELT_NUMBER_USER = nusselt_user_2(iel)

  end function NUSSELT_NUMBER_USER


!-------------------------------------------------------------------------------
!! User example 1: variable property correction for Nu in forced convections
!! Best suited for flow at supercritical pressures
!! Based on work of Krasnosechekov and Protopopov
  function nusselt_user_1(iel)

    use numvar
    use field
    use mesh
    use parall
    use optcal

    implicit none

    double precision :: nusselt_user_1

!! Aguments
    integer :: iel

!! Local variables
    double precision, dimension(:), pointer :: tbulk, A_twall, Dh
    double precision, dimension(:,:), pointer :: B_vel
    double precision :: Tpc, twall0, twall
    double precision :: rhob, enthb, visb, cpa, cpb, condb, Reb, Prb, &
                        rhow, enthw
    double precision :: htc, Nuf, n, tol
    double precision :: intp, intp1, intp2
    integer :: cpro_index

    if (xyzcen(3,iel) < 0.03d0 .or. xyzcen(3,iel) > 0.63d0) then
      nusselt_user_1 = 1.d0
      return
    end if

    Tpc = 377.51d0

    call field_get_val_s_by_name('SubchannelBulkTemperature', tbulk)
    call field_get_val_s_by_name('SubchannelWallTemperature', A_twall)
    call field_get_val_s_by_name('HydraulicDiameter', Dh)
    call field_get_val_v_by_name('BulkVelocity', B_vel)

    !! Obtain subchannel bulk properties based on bulk temperature
    cpro_index = 1
    do while (tbulk(iel) >= user_prop(cpro_index,2))
      cpro_index = cpro_index + 1
    end do
    intp = user_prop(cpro_index,2)-user_prop(cpro_index-1,2)
    intp1 = (tbulk(iel)-user_prop(cpro_index-1,2))/intp
    intp2 = 1.d0-intp1

    !! Bulk density
    rhob = user_prop(cpro_index-1,3)*intp2+(user_prop(cpro_index,3))*intp1
    !! Bulk enthalpy
    enthb = user_prop(cpro_index-1,1)*intp2+(user_prop(cpro_index,1))*intp1
    !! Bulk viscosity
    visb = user_prop(cpro_index-1,5)*intp2+(user_prop(cpro_index,5))*intp1
    !! Bulk heat capacity
    cpb = user_prop(cpro_index-1,4)*intp2+(user_prop(cpro_index,4))*intp1
    !! Bulk thermal conductivity
    condb = user_prop(cpro_index-1,6)*intp2+(user_prop(cpro_index,6))*intp1

    !! Calculate subchannel bulk Reynolds number
    Reb = rhob*B_vel(3,iel)*Dh(iel)/visb
    !! Calculate subchannel bulk Prandtl number
    Prb = visb*cpb/condb

    !! Calculate the Nusselt number
    twall0 = A_twall(iel)
    tol = 1.d-3
    !open(100,file='output_data.txt')
    do
      cpro_index = 1
      do while (twall0 >= user_prop(cpro_index,2))
        cpro_index = cpro_index + 1
      end do
      intp = user_prop(cpro_index,2)-user_prop(cpro_index-1,2)
      intp1 = (twall0-user_prop(cpro_index-1,2))/intp
      intp2 = 1.d0-intp1

      !! Wall density
      rhow = user_prop(cpro_index-1,3)*intp2+(user_prop(cpro_index,3))*intp1
      !! Wall enthalpy
      enthw = user_prop(cpro_index-1,1)*intp2+(user_prop(cpro_index,1))*intp1

      if ((tbulk(iel)<=twall0 .and. twall0<=Tpc) &
         .or. (tbulk(iel)>=1.2d0*Tpc .and. tbulk(iel)<=twall0)) then
        n = 0.4d0
      else if ((tbulk(iel)<=Tpc) .and. (Tpc<=twall0)) then
        n = 0.4d0+0.2d0*(twall0/Tpc-1.d0)
      else if ((Tpc<=tbulk(iel)) .and. (tbulk(iel)<=1.2d0*Tpc) &
              .and. (tbulk(iel)<=twall0)) then
        n = 0.4d0+0.2d0*(twall0/Tpc-1.d0)*(1.d0-5.d0*(tbulk(iel)/Tpc-1.d0))
      else
        if(irangp.le.0) write(*,*)'Warning: error in calculating Nusselt Number&
                                   & - cannot obtain parameter n'
      end if

      cpa = (enthw-enthb)/(twall0-tbulk(iel)+1.d-30)
      Nuf = 0.0183d0*(Reb**0.82)*(Prb**0.5)*((rhow/rhob)**0.3)*((cpa/cpb)**n)
      htc = Nuf*condb/Dh(iel)
      twall = tbulk(iel)+value_rw/(htc+1.d-30)
      if(abs(twall-twall0) < tol) exit
      twall0 = twall
      !! check data
      !if(iel==417)then
      !  if(irangp.le.0) then
      !    write(100,'(I8,9f12.2)') ntcabs, tbulk(iel), htc, Reb, Prb, rhow, &
      !                             rhob, cpa, cpb, twall
      !  end if
      !end if
    end do
    !close (100)

    nusselt_user_1 = Nuf

  end function nusselt_user_1

!-------------------------------------------------------------------------------
!! User example 2: variable property correction for Nu in forced convections
!! Best for subcritical conditions
!! Petukhov and Popov (1983)
  function nusselt_user_2(iel)

    use numvar
    use field
    use mesh
    use parall
    use optcal

    implicit none

    double precision :: nusselt_user_2

!! Aguments
    integer          :: iel

!! Local variables
    logical :: condition, condition1, condition2
    double precision, dimension(:), pointer :: tbulk, A_twall, A_tflux, Dh
    double precision, dimension(:,:), pointer :: B_vel
    double precision :: twall0, twall
    double precision :: rhob, visb, kvisb, cpb, condb, betab, Reb, Prb, rhow
    double precision :: htc, Nuf, cf, tol
    double precision :: Gr, Bo, Nur
    double precision :: intp, intp1, intp2
    integer :: cpro_index, bo_index

    double precision, dimension(53) :: tab_bo, tab_Nur

    double precision, parameter :: g = 9.81

    tab_bo = (/ 1.00000000e-07,1.15139540e-07,1.32571137e-07,1.52641797e-07,   &
                1.75751062e-07,2.02358965e-07,2.32995181e-07,2.68269580e-07,   &
                3.08884360e-07,3.55648031e-07,4.09491506e-07,4.71486636e-07,   &
                5.42867544e-07,6.25055193e-07,7.19685673e-07,8.28642773e-07,   &
                9.54095476e-07,1.09854114e-06,1.26485522e-06,1.45634848e-06,   &
                1.67683294e-06,1.93069773e-06,2.10000000e-06,2.13000000e-06,   &
                2.20000000e-06,2.22299648e-06,2.55954792e-06,2.94705170e-06,   &
                3.39322177e-06,3.90693994e-06,4.49843267e-06,5.17947468e-06,   &
                5.96362332e-06,6.86648845e-06,7.90604321e-06,9.10298178e-06,   &
                1.04811313e-05,1.20679264e-05,1.38949549e-05,1.59985872e-05,   &
                1.84206997e-05,2.12095089e-05,2.44205309e-05,2.81176870e-05,   &
                3.23745754e-05,3.72759372e-05,4.29193426e-05,4.94171336e-05,   &
                5.68986603e-05,6.55128557e-05,7.54312006e-05,8.68511374e-05,   &
                1.00000000e-04 /)

    tab_Nur = (/0.99416288,0.99326366,0.99222168,0.99101589,0.98961794,        &
                0.98799536,0.98610953,0.98391427,0.98135397,0.97836115,        &
                0.9748532 ,0.97072775,0.96585626,0.96007448,0.95316795,        &
                0.94484889,0.93471699,0.92218776,0.90634331,0.88559737,        &
                0.85666018,0.80958495,0.74395438,0.4719312 ,0.47844216,        &
                0.48054972,0.50978019,0.54023444,0.57190827,0.60479573,        &
                0.63889063,0.67418736,0.71068198,0.74837281,0.787261  ,        &
                0.82735146,0.86865282,0.91117783,0.95494356,0.99997139,        &
                1.04628673,1.09391984,1.14290416,1.19327764,1.24508184,        &
                1.29836195,1.35316657,1.40954824,1.46756163,1.52726537,        &
                1.58872091,1.6519927 ,1.7171481 /)

!    call field_get_val_s_by_name('SubchannelBulkTemperature', tbulk)
    call field_get_val_s_by_name('AverageTemperature', tbulk)
    call field_get_val_s_by_name('SubchannelWallTemperature', A_twall)
    call field_get_val_s_by_name('SubchannelWallHeatFlux', A_tflux)
    call field_get_val_s_by_name('HydraulicDiameter', Dh)
    call field_get_val_v_by_name('BulkVelocity', B_vel)

    !! Obtain subchannel bulk properties based on bulk temperature
    cpro_index = 1
    do while (tbulk(iel) >= user_prop(cpro_index,1))
      cpro_index = cpro_index + 1
    end do
    intp = user_prop(cpro_index,1)-user_prop(cpro_index-1,1)
    intp1 = (tbulk(iel)-user_prop(cpro_index-1,1))/intp
    intp2 = 1.d0-intp1

    !! Bulk density
    rhob = user_prop(cpro_index-1,3)*intp2+(user_prop(cpro_index,3))*intp1
    !! Bulk heat capacity
    cpb = user_prop(cpro_index-1,4)*intp2+(user_prop(cpro_index,4))*intp1
    !! Bulk viscosity
    visb = user_prop(cpro_index-1,5)*intp2+(user_prop(cpro_index,5))*intp1
    !! Bulk thermal conductivity
    condb = user_prop(cpro_index-1,6)*intp2+(user_prop(cpro_index,6))*intp1

    !! Bulk thermal expansion coefficient
    betab = user_prop(cpro_index-1,7)*intp2+(user_prop(cpro_index,7))*intp1
    kvisb = visb/rhob

    !! Calculate subchannel bulk Reynolds number
    Reb = rhob*abs(B_vel(3,iel))*Dh(iel)/visb
    !! Calculate subchannel bulk Prandtl number
    Prb = visb*cpb/condb

    !! calculate subchannel bulk Grashof number
    Gr = betab*g*abs(A_tflux(iel))*(Dh(iel)**4)/condb/(kvisb**2)
    !! calculate subchannel bulk Buoyancy parameter
    Bo = Gr/Reb**3.425/Prb**0.8

    if(Reb<=2300.d0) then
      Nuf = 4.36d0

    else
      !! Dittus-Boelter correlation (not used)
      if (heating)then
        Nuf = 0.021d0*(Reb**0.8d0)*(Prb**0.4d0)
      else
        Nuf = 0.023d0*(Reb**0.8d0)*(Prb**0.3d0)
      end if

      !! It is better the correlation is used within the range of the below conditions
      condition1 = Reb>=1.d4
      condition2 = Prb>=0.7d0 .and. Prb<=1.d2
      condition = condition1 .and. condition2

      if(.not.condition) then
        !print*, 'Warning! Nusselt number calculation error: out of range at ',&
        !        'Subchannel ID ', isub
      end if

      !! To introduce corrections to consider variable properties
      !! Calculate the Nusselt number
      twall0 = A_twall(iel)
      tol = 1.d-3
      do
        cpro_index = 1
        do while (twall0 >= user_prop(cpro_index,1))
          cpro_index = cpro_index + 1
        end do
        intp = user_prop(cpro_index,1)-user_prop(cpro_index-1,1)
        intp1 = (twall0-user_prop(cpro_index-1,1))/intp
        intp2 = 1.d0-intp1

        !! Wall density
        rhow = user_prop(cpro_index-1,3)*intp2+(user_prop(cpro_index,3))*intp1

        if (Reb<=1.d-30) then
          cf = 1.d30
        else if(Reb<2.3d3 .and. Reb>1.d-30) then
          cf = 64.d0/Reb
        else
          cf = 1.d0/(1.82d0*log10(Reb/8.d0))**2
        end if
        cf = cf*((rhow/rhob)**0.4)
        Nuf = cf/8.d0*Reb*Prb
        Nuf = Nuf/(1.d0+900.d0/Reb+12.7d0*((cf/8.d0)**0.5)*(Prb**(2.d0/3.d0)-1.d0))
        htc = Nuf*condb/Dh(iel)
        twall = tbulk(iel)+value_rw/(htc+1.d-30)
        if(abs(twall-twall0) < tol) exit
        twall0 = twall
      end do

    end if

    !!buoyancy correction
    Nur = 1.d0
!    Bo_index = 1
!    do while (Bo >= tab_Bo(Bo_index))
!      if (Bo_index>53) then
!        Bo = tab_Bo(53)
!        if(irangp.le.0) &
!          write(*,*) 'Warning: Buoyancy parameter exceeds upper limit'
!        exit
!      end if
!      Bo_index = Bo_index + 1
!    end do
!    intp = tab_Bo(Bo_index)-tab_Bo(Bo_index-1)
!    intp1 = (Bo-tab_Bo(Bo_index-1))/intp
!    intp2 = 1.d0-intp1
!    Nur = tab_Nur(Bo_index-1)*intp2+(tab_Nur(Bo_index))*intp1

    nusselt_user_2 = Nuf*Nur

  end function nusselt_user_2

!===============================================================================
!! This routine is used to calculate porous media resistance
  subroutine resistance_user()

    use numvar
    use field
    use mesh
    use cstphy
    use parall

    implicit none

    double precision, dimension(:), pointer :: axial_resis_f, axial_resis_d

    call field_get_val_s_by_name('AxialResistanceFriction', axial_resis_f)
    call field_get_val_s_by_name('AxialResistanceDrag', axial_resis_d)

    if(iporo_cal==1 .or. iporo_cal==2) then
      call round_robin_porous_I()
    else if(iporo_cal==3) then
!      call ballooned_7x7_porous_II_advanced()
!      call ballooned_7x7_porous_II_basic()
      call round_robin_porous_II_basic()
    end if

    contains
!-------------------------------------------------------------------------------
!! 7x7 ballooned rod bundle (porous method II advanced)
    subroutine ballooned_7x7_porous_II_advanced()

      implicit none

      integer :: i, j, iel

      double precision :: cf, Kc, len_porous, ubulk, Dh, dpdzf, dpdzd, Re, PbyD

      double precision, dimension(:), pointer :: filter_ID, F_den, F_mu
      double precision, dimension(:,:), pointer :: F_vel
      double precision :: a, b1, b2, n

      call field_get_val_s_by_name('FilterID', filter_ID)
      call field_get_val_s_by_name('FilteredDensity', F_den)
      call field_get_val_s_by_name('FilteredDynamicViscosity', F_mu)
      call field_get_val_v_by_name('FilteredVelocity', F_vel)

      Kc = 0.45d0
      len_porous = 0.0762d0
      ubulk = 1.74d0
      Dh = 0.0139d0
      PbyD = 1.37d0

      do iel=1,ncel
        j = filter_ID(iel)
        if(j>0.d0)then

          Re = F_den(iel)*F_vel(3,iel)*Dh/F_mu(iel)

         !! Cheng & Todreas (1986)
          if(Re<=2300.d0)then
            a  = 35.55d0
            b1 = 263.7d0
            b2 = -190.2d0
            n  = 1.d0
          else
            a  = 0.1339d0
            b1 = 0.09059d0
            b2 = -0.09926d0
            n  = 0.18d0
          end if

          cf = a+b1*(PbyD-1.d0)+b2*(PbyD-1.d0)**2
          cf = cf/(Re**n+1.d-15)

          dpdzf = 0.5d0*cf*F_den(iel)*abs(F_vel(3,iel))*F_vel(3,iel)/Dh
          dpdzd = 0.5d0*Kc*F_den(iel)*abs(ubulk)*ubulk/len_porous
          axial_resis_f(iel) = -dpdzf
          axial_resis_d(iel) = -dpdzd
        end if
      end do

    end subroutine ballooned_7x7_porous_II_advanced

!-------------------------------------------------------------------------------
!! 7x7 ballooned rod bundle (porous method II basic)
    subroutine ballooned_7x7_porous_II_basic()

      implicit none

      integer :: i, j, iel

      double precision :: cf, Kc, len_porous, ubulk, Dh, dpdzf, dpdzd, Re, PbyD

      double precision, dimension(:), pointer :: filter_ID, den, mu
      double precision, dimension(:,:), pointer :: vel
      double precision :: a, b1, b2, n

      call field_get_val_s_by_name('FilterID', filter_ID)

      call field_get_val_s(iprpfl(iviscl), mu)
      call field_get_val_s(icrom, den)
      call field_get_val_v(ivarfl(iu),vel)

      Kc = 0.45d0
      len_porous = 0.0762d0
      ubulk = 1.74d0
      Dh = 0.0139d0
      PbyD = 1.37d0

      do iel=1,ncel
        j = filter_ID(iel)
        if(j>0.d0)then

          Re = den(iel)*vel(3,iel)*Dh/mu(iel)

         !! Cheng & Todreas (1986)    
          if(Re<=2300.d0)then
            a  = 35.55d0
            b1 = 263.7d0
            b2 = -190.2d0
            n  = 1.d0
          else
            a  = 0.1339d0
            b1 = 0.09059d0
            b2 = -0.09926d0
            n  = 0.18d0
          end if

          cf = a+b1*(PbyD-1.d0)+b2*(PbyD-1.d0)**2
          cf = cf/(Re**n+1.d-15)

          dpdzf = 0.5d0*cf*den(iel)*abs(vel(3,iel))*vel(3,iel)/Dh
          dpdzd = 0.5d0*Kc*den(iel)*abs(ubulk)*ubulk/len_porous
          axial_resis_f(iel) = -dpdzf
          axial_resis_d(iel) = -dpdzd
        end if
      end do

    end subroutine ballooned_7x7_porous_II_basic

!-------------------------------------------------------------------------------
!! Round Robin 5x5 rod bundle (porous method II basic)
    subroutine round_robin_porous_II_basic()

      implicit none

      integer :: i, j, iel

      double precision :: cf, Kc, len_porous, ubulk, Dh, dpdzf, dpdzd, Re, PbyD

      double precision, dimension(:), pointer :: filter_ID, den, mu
      double precision, dimension(:,:), pointer :: vel
      double precision :: a, b1, b2, n

      call field_get_val_s_by_name('FilterID', filter_ID)
      
      call field_get_val_s(iprpfl(iviscl), mu)
      call field_get_val_s(icrom, den)
      call field_get_val_v(ivarfl(iu),vel)

      Kc = 0.429d0
      len_porous = 0.008d0
      ubulk = uref
      Dh = 0.01178d0 !(center channel)
      PbyD = 1.326d0

      do iel=1,ncel
        j = filter_ID(iel)
        if(j>0.d0)then

          Re = den(iel)*vel(3,iel)*Dh/mu(iel)

         !! Cheng & Todreas (1986)
          if(Re<=2300.d0)then
            a  = 35.55d0
            b1 = 263.7d0
            b2 = -190.2d0
            n  = 1.d0
          else
            a  = 0.1339d0
            b1 = 0.09059d0
            b2 = -0.09926d0
            n  = 0.18d0
          end if

          cf = a+b1*(PbyD-1.d0)+b2*(PbyD-1.d0)**2
          cf = cf/(Re**n+1.d-15)

          dpdzf = 0.5d0*cf*den(iel)*abs(vel(3,iel))*vel(3,iel)/Dh
          dpdzd = 0.5d0*Kc*den(iel)*abs(ubulk)*ubulk/len_porous
          axial_resis_f(iel) = -dpdzf
          axial_resis_d(iel) = -dpdzd
!          axial_resis_f(iel) = 0.d0
!          axial_resis_d(iel) = 0.d0
        end if
      end do

    end subroutine round_robin_porous_II_basic

!-------------------------------------------------------------------------------
!! Round Robin 5x5 rod bundle (porous method I)
    subroutine round_robin_porous_I()

      implicit none

      integer :: iel

      double precision :: Kc, len_porous, ubulk, dpdzf, dpdzd

      double precision, dimension(:), pointer :: den, cf, Dh
      double precision, dimension(:,:), pointer :: A_vel

      call field_get_val_s_by_name('AverageDensity', den)
      call field_get_val_v_by_name('AverageVelocity', A_vel)
      call field_get_val_s_by_name('AxialFrictionFactor', cf)
      call field_get_val_s_by_name('HydraulicDiameter', Dh)

      Kc = 0.429d0
      len_porous = 0.008d0

      !Exercise 1-1
      ubulk = uref
      !Exercise 2-1
      !ubulk = 5.5985d0

      do iel=1,ncel
        dpdzf = 0.5d0*cf(iel)*den(iel)*abs(A_vel(3,iel))*A_vel(3,iel)/Dh(iel)
        dpdzd = 0.5d0*Kc*den(iel)*abs(ubulk)*ubulk/len_porous
        axial_resis_f(iel) = -dpdzf
        axial_resis_d(iel) = -dpdzd
      end do

    end subroutine round_robin_porous_I
!-------------------------------------------------------------------------------

  end subroutine resistance_user


!===============================================================================
! The follwoing subroutines are used for post-processings
! Subroutines included here are called in cs_user_extra_operations.f90
!===============================================================================

!===============================================================================
!! This routine is used to extract rod wall temperatures along the circumfrential
!! direction at a specified axial location (The following routine is currently 
!! for the 5x5 PWR bundle and should be adjusted according to different data 
!! processing requirements)
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
  subroutine circumferential_twall ()
  
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
    double precision, dimension(:), allocatable :: angle, angle_g, &
                                                   twall_para,     &
                                                   twall_para_g
    double precision, dimension(:), pointer :: twall, subchannel_ID
    double precision :: x,x1,x2,y,y1,y2,radius,angle0,cell_thichness
    character(len=50) :: chaine, cfilename
    logical :: condition_a, condition_b
    integer :: nvariables
    double precision, dimension(:,:), allocatable :: array

    ! User specified parameters
    ! Local temperature (true) or subchannel bulk temperature (false) based
    logical :: local_based = .true.
    ! Geometrical parameters
    double precision, parameter :: pitch = 0.0126d0, &
                                   tol1 = 1.d-3, &
                                   tol2 = 1.d-3, &
                                   fuel_rod_diameter = 0.0095d0, &
                                   x11 = 0.00785d0, &
                                   y11 = 0.05825d0
    integer, parameter :: num_axial_loc = 1
    double precision, dimension(num_axial_loc) ::  axial_loc

    axial_loc(1) = -2.5d-4

    if (iscalt.le.0) return

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)

    if(local_based)then
      call field_get_val_s_by_name('WallTemperature1', twall)
    else
      call field_get_val_s_by_name('WallTemperature0', twall)
    end if

    allocate(rod_center_xyz(5,5,2), rod_radius(5,5))

    arr = 100
    allocate(angle(arr),angle_g(arr),twall_para(arr),twall_para_g(arr))

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

    do kk=1,num_axial_loc
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

      !! In the case that axial_loc coincides with grid lines
      condition_a = abs(cdgfbo(3,ifac)-axial_loc(kk)) < (0.5d0-tol1)*cell_thichness
      condition_b = abs(cdgfbo(3,ifac)-axial_loc(kk)) < (0.5d0+tol1)*cell_thichness&
                    .and. (cdgfbo(3,ifac)-axial_loc(kk)) > 0.d0
      condition_b = .not.(condition_a) .and. condition_b

      if(radius < rod_radius(j,k)*(1.d0+tol2) .and. &
          (condition_a .or. condition_b))then
        if((cdgfbo(2,ifac)-y) >= 0)then
          angle(mm) = acos((cdgfbo(1,ifac)-x)/radius)/pi*180.d0
        else
          angle(mm) = 360.d0-acos((cdgfbo(1,ifac)-x)/radius)/pi*180.d0
        end if
        twall_para(mm) = twall(ifac)
        mm = mm+1
      end if
    enddo
    endif

    if(irangp.ge.0)call parall_bcast_i(ii, mm)

    do i=1,mm-1
      if(irangp.ge.0)then
        call parall_bcast_r(ii, angle(i))
        call parall_bcast_r(ii, twall_para(i))
      endif

      angle_g(nm) = angle(i)
      twall_para_g(nm) = twall_para(i)
      nm = nm+1

      if(irangp.ge.0)then
        call parall_bcast_i(ii, nm)
        call parall_bcast_r(ii, angle_g(nm))
        call parall_bcast_r(ii, twall_para_g(nm))
      endif
    enddo

    enddo

    if (ntcabs.eq.ntmabs) then
    impout = impusr(1)
    if (irangp.le.0) then
      write(cfilename,'(a22,2i02.2,a3,I5,a4)')'circumferential_twall_',j,k,&
                                        '_z=',int(axial_loc(kk)*1.d3),'.dat'
      open(impout,file=Trim(cfilename))
      write(impout,*)  &
           'degree(o)  t_wall(C)'
    endif
    endif

    nvariables = 7
    allocate (array(nm-1,nvariables))


    array(1:nm-1,1) = angle_g(1:nm-1)
    array(1:nm-1,2) = twall_para_g(1:nm-1)

    call sorting(array)

    angle_g(1:nm-1) = array(1:nm-1,1)
    twall_para_g(1:nm-1) = array(1:nm-1,2)

    deallocate(array)

    do i=1, nm-1
      if (ntcabs.eq.ntmabs) then
        if (irangp.le.0) write(impout,'(2f10.5)') angle_g(i), twall_para_g(i)
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
!    enddo

    if (ntcabs.eq.ntmabs) then
      if (irangp.le.0) close(impout)
    endif

    end do
    end do
    end do

    deallocate(rod_center_xyz, rod_radius)
    deallocate(angle, angle_g, twall_para, twall_para_g)
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

  end subroutine circumferential_twall


!===============================================================================
!! This routine is used to extract coolant channel wall heat flux along the 
!! circumfrential direction at a specified axial location of a prismatic
!! HTGR (GA) fuel assembly.
  subroutine circumferential_tflux_HTGR_GA()
  
    use numvar
    use optcal
    use field
    use mesh
    use parall
    use cstphy
    use entsor
    use cs_c_bindings

    implicit none

    integer :: i, ii, j, jj, k, kk, iel, ilelt, ifac, nlelt
    integer, allocatable, dimension(:) :: lstelt
    integer :: nm, mm, arr
    integer :: impout
    double precision, dimension(:,:), allocatable :: channel_center_xyz
    double precision, dimension(:), allocatable :: channel_radius
    double precision, dimension(:), allocatable :: angle, angle_g, &
                                                   twall_para,     &
                                                   twall_para_g,   &
                                                   tflux_para,     &
                                                   tflux_para_g
    double precision, dimension(:), pointer :: twall, tflux, subchannel_ID
    double precision :: x,x1,x2,y,y1,y2,radius,angle0,cell_thichness
    character(len=50) :: chaine, cfilename
    logical :: condition_a, condition_b
    integer :: nvariables
    double precision, dimension(:,:), allocatable :: array

    ! User specified parameters
    ! Geometrical parameters
    double precision, parameter :: tol1 = 1.d-3, &
                                   tol2 = 1.d-3, &
                                   fuel_channel_diameter1 = 12.70d-3, &
                                   fuel_channel_diameter2 = 15.88d-3

    integer, parameter :: num_axial_loc = 1
    double precision, dimension(num_axial_loc) ::  axial_loc

    interface
      subroutine cal_boundary_thermal_flux()        &
        bind(C, name='cal_boundary_thermal_flux')
        use, intrinsic :: iso_c_binding
        implicit none
      end subroutine cal_boundary_thermal_flux
    end interface

    axial_loc(1) = 2.d-3

    if (iscalt.le.0) return

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)


    call field_get_val_s_by_name('WallThermalFlux', tflux)
    call field_get_val_s_by_name('WallTemperature0', twall)

    call cal_boundary_thermal_flux()

    allocate(channel_center_xyz(13,2), channel_radius(13))

    arr = 100
    allocate(angle(arr),angle_g(arr))
    allocate(twall_para(arr),twall_para_g(arr))
    allocate(tflux_para(arr),tflux_para_g(arr))

    channel_radius(1)    = fuel_channel_diameter1/2.d0
    channel_radius(2:13) = fuel_channel_diameter2/2.d0

    !Note: Orgin of the coordinate system is at the bottom left corner of the 
    !fuel block
    channel_center_xyz(1,1) = 103.923d-3
    channel_center_xyz(1,2) = 147.438d-3
    channel_center_xyz(2,1) = 103.923d-3
    channel_center_xyz(2,2) = 114.875d-3
    channel_center_xyz(3,1) = 75.723d-3
    channel_center_xyz(3,2) = 131.156d-3
    channel_center_xyz(4,1) = 103.923d-3
    channel_center_xyz(4,2) = 82.3123d-3
    channel_center_xyz(5,1) = 75.723d-3
    channel_center_xyz(5,2) = 98.5936d-3
    channel_center_xyz(6,1) = 103.923d-3
    channel_center_xyz(6,2) = 49.7497d-3
    channel_center_xyz(7,1) = 75.723d-3
    channel_center_xyz(7,2) = 66.031d-3
    channel_center_xyz(8,1) = 47.523d-3
    channel_center_xyz(8,2) = 82.3123d-3
    channel_center_xyz(9,1) = 103.923d-3
    channel_center_xyz(9,2) = 17.1871d-3
    channel_center_xyz(10,1) = 75.723d-3
    channel_center_xyz(10,2) = 33.4684d-3
    channel_center_xyz(11,1) = 47.523d-3
    channel_center_xyz(11,2) = 49.7497d-3
    channel_center_xyz(12,1) = 47.523d-3
    channel_center_xyz(12,2) = 17.1871d-3
    channel_center_xyz(13,1) = 19.323d-3
    channel_center_xyz(13,2) = 33.4684d-3

    !! Get all wall faces
    allocate(lstelt(nfabor))
    chaine = all_walls
    call getfbr(chaine, nlelt, lstelt)

!! Disable solid cells in case of CHT
    if (fluid_solid) call cs_porous_model_set_has_disable_flag(1)

    do kk=1,num_axial_loc
    do j=1,13

    nm = 1
    angle = 0.d0
    do ii=0, nrangp-1
    mm = 1
    if(irangp==ii.or.irangp==-1)then
    do ilelt = 1, nlelt
      ifac = lstelt(ilelt)
      iel = ifabor(ifac)
      if(cell_is_active(iel)==0)cycle
      jj = subchannel_ID(iel)
      cell_thichness = cell_height(jj)
      x = channel_center_xyz(j,1)
      y = channel_center_xyz(j,2)
      radius = sqrt((cdgfbo(1,ifac)-x)**2+(cdgfbo(2,ifac)-y)**2)

      !! In the case that axial_loc coincides with grid lines
      condition_a = abs(cdgfbo(3,ifac)-axial_loc(kk)) < (0.5d0-tol1)*cell_thichness
      condition_b = abs(cdgfbo(3,ifac)-axial_loc(kk)) < (0.5d0+tol1)*cell_thichness&
                    .and. (cdgfbo(3,ifac)-axial_loc(kk)) > 0.d0
      condition_b = .not.(condition_a) .and. condition_b

      if(radius < channel_radius(j)*(1.d0+tol2) .and. &
          (condition_a .or. condition_b))then
        if((cdgfbo(2,ifac)-y) >= 0)then
          angle(mm) = acos((cdgfbo(1,ifac)-x)/radius)/pi*180.d0
        else
          angle(mm) = 360.d0-acos((cdgfbo(1,ifac)-x)/radius)/pi*180.d0
        end if
        twall_para(mm) = twall(ifac)
        tflux_para(mm) = tflux(ifac)
        mm = mm+1
      end if
    enddo
    endif

    if(irangp.ge.0)call parall_bcast_i(ii, mm)

    do i=1,mm-1
      if(irangp.ge.0)then
        call parall_bcast_r(ii, angle(i))
        call parall_bcast_r(ii, twall_para(i))
        call parall_bcast_r(ii, tflux_para(i))
      endif

      angle_g(nm) = angle(i)
      twall_para_g(nm) = twall_para(i)
      tflux_para_g(nm) = tflux_para(i)
      nm = nm+1

      if(irangp.ge.0)then
        call parall_bcast_i(ii, nm)
        call parall_bcast_r(ii, angle_g(nm))
        call parall_bcast_r(ii, twall_para_g(nm))
        call parall_bcast_r(ii, tflux_para_g(nm))
      endif
    enddo

    enddo

    if (ntcabs.eq.ntmabs) then
    impout = impusr(1)
    if (irangp.le.0) then
      write(cfilename,'(a22,i02.2,a3,I1,a4)')'circumferential_tflux_',j,&
                                        '_z=',int(axial_loc(kk)*1.d3),'.dat'
      open(impout,file=Trim(cfilename))
      write(impout,*)  &
           'degree(o)  t_wall(C)   t_flux(J/m2s)'
    endif
    endif

    nvariables = 7
    allocate (array(nm-1,nvariables))


    array(1:nm-1,1) = angle_g(1:nm-1)
    array(1:nm-1,2) = twall_para_g(1:nm-1)
    array(1:nm-1,3) = tflux_para_g(1:nm-1)

    call sorting(array)

    angle_g(1:nm-1) = array(1:nm-1,1)
    twall_para_g(1:nm-1) = array(1:nm-1,2)
    tflux_para_g(1:nm-1) = array(1:nm-1,3)

    deallocate(array)

    do i=1, nm-1
      if (ntcabs.eq.ntmabs) then
        if (irangp.le.0) write(impout,'(2f10.4, ES15.4)') angle_g(i),      &
                                                          twall_para_g(i), &
                                                          tflux_para_g(i)
      endif
    end do

!! Serial version
!    do ilelt = 1, nlelt
!      ifac = lstelt(ilelt)
!      x = channel_center_xyz(1,1,1)
!      y = channel_center_xyz(1,1,2)
!      radius = sqrt((cdgfbo(1,ifac)-x)**2+(cdgfbo(2,ifac)-y)**2)
!      if(radius < channel_radius(1,1)*(1.d0+tol2) &
!         .and. cdgfbo(3,ifac) > -0.5d0*thichness)then
!        if((cdgfbo(2,ifac)-y) >= 0)then
!          angle0 = acos((cdgfbo(1,ifac)-x)/radius)/pi*180.d0
!        else
!          angle0 = 360.d0-acos((cdgfbo(1,ifac)-x)/radius)/pi*180.d0
!        end if  
!      if (irangp.le.0) write(impout,'(2f10.5)')angle0, tflux(ifac)
!      end if
!    enddo

    if (ntcabs.eq.ntmabs) then
      if (irangp.le.0) close(impout)
    endif

    end do
    end do

!! Re-enable solid cells in case of CHT
    if(fluid_solid) call cs_porous_model_set_has_disable_flag(0)

    deallocate(channel_center_xyz, channel_radius)
    deallocate(angle, angle_g)
    deallocate(twall_para, twall_para_g)
    deallocate(tflux_para, tflux_para_g)
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

  end subroutine circumferential_tflux_HTGR_GA


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
    ! Local temperature (true) or subchannel bulk temperature (false) based
    logical :: local_based = .true.
    ! Geometrical parameters
    double precision, parameter :: pitch = 0.0126d0, &
                                   tol1 = 1.d-3, &
                                   tol2 = 1.d-3, &
                                   fuel_rod_diameter = 0.0095d0, &
                                   x11 = 0.00785d0, &
                                   y11 = 0.05825d0, &
                                   start_location = -3950.d-3, &
                                   end_location = 100.d-3, &
                                   increment = 5.d-3

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
    if (itherm.ne.1) then
      call field_get_val_s_by_name('temperature', cvar_scalt)
    else
      call field_get_val_s(ivarfl(ivart), cvar_scalt)
    end if
    call field_get_val_s(icrom, cpro_rom)
    call field_get_val_v(ivarfl(iu), cvar_vel)
    call field_get_val_s(ivarfl(ipr), cvar_p)

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)

    if(local_based)then
      call field_get_val_s_by_name('WallTemperature1', twall)
    else
      call field_get_val_s_by_name('WallTemperature0', twall)
    end if

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
           'Axial_coordinate(m)  t_avg(C)  t_bulk(C)  w_bulk(m/s)   p_avg(Pa)&
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
!! This routine is specifically used for the Round Robin case to extract
!! subchannel velocities which are based on the subchannel indices in the 
!! report (z=235mm, eq. z=-44mm in the periodic simulation) 
  subroutine round_robin_subchannel_vel()

    use numvar
    use optcal
    use cstphy
    use field
    use mesh
    use parall
    use entsor

    implicit none

    integer :: i, j, iel
    double precision, dimension(36) :: sub_vel
    double precision, dimension(3,36) :: sub_coord
    double precision, dimension (3,n_cell) :: temp
    
    double precision, dimension(n_cell) :: subchannel_vel
    double precision, pointer, dimension(:,:) :: vel
    double precision, pointer, dimension(:) :: subchannel_ID
    
    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_v(ivarfl(iu), vel)

    do i=1, n_cell
      temp(1,i) = x_channel(i)
      temp(2,i) = y_channel(i)
      temp(3,i) = z_channel(i)
    end do

!!  mesh scheme 1
    sub_coord(1,22) = 0.0063d0
    sub_coord(2,22) = 0.0063d0
    sub_coord(3,22) = -0.043d0

    sub_coord(1,23) = 0.0063d0+0.0126d0
    sub_coord(2,23) = 0.0063d0
    sub_coord(3,23) = -0.043d0

    sub_coord(1,24) = 0.0063d0+0.022825d0
    sub_coord(2,24) = 0.0063d0
    sub_coord(3,24) = -0.043d0

    sub_coord(1,28) = sub_coord(1,22)
    sub_coord(2,28) = 0.0063d0+0.0126d0
    sub_coord(3,28) = -0.043d0

    sub_coord(1,29) = sub_coord(1,23)
    sub_coord(2,29) = 0.0063d0+0.0126d0
    sub_coord(3,29) = -0.043d0

    sub_coord(1,30) = sub_coord(1,24)
    sub_coord(2,30) = 0.0063d0+0.0126d0
    sub_coord(3,30) = -0.043d0

    sub_coord(1,34) = sub_coord(1,22)
    sub_coord(2,34) = 0.0063d0+0.022825d0
    sub_coord(3,34) = -0.043d0

    sub_coord(1,35) = sub_coord(1,23)
    sub_coord(2,35) = 0.0063d0+0.022825d0
    sub_coord(3,35) = -0.043d0

    sub_coord(1,36) = sub_coord(1,24)
    sub_coord(2,36) = 0.0063d0+0.022825d0
    sub_coord(3,36) = -0.043d0

    subchannel_vel = 0.d0
    do iel=1, ncel
      j = int(subchannel_ID(iel))
      if(j>0) &
        subchannel_vel(j) = subchannel_vel(j)+vel(3,iel)*volume(iel)
    end do

    if (irangp.ge.0) then
      call parrsm(n_cell, subchannel_vel)
    endif

    do i=1,n_cell
      subchannel_vel(i) = subchannel_vel(i)/subchannel_vol(i)
    end do

    do i=1, n_cell
      if(point_distance(temp(:,i),sub_coord(:,22))<1.d-10)&
        sub_vel(22) = subchannel_vel(i)
      if(point_distance(temp(:,i),sub_coord(:,23))<1.d-10)&
        sub_vel(23) = subchannel_vel(i)
      if(point_distance(temp(:,i),sub_coord(:,24))<1.d-10)&
        sub_vel(24) = subchannel_vel(i)
      if(point_distance(temp(:,i),sub_coord(:,28))<1.d-10)&
        sub_vel(28) = subchannel_vel(i)
      if(point_distance(temp(:,i),sub_coord(:,29))<1.d-10)&
        sub_vel(29) = subchannel_vel(i)
      if(point_distance(temp(:,i),sub_coord(:,30))<1.d-10)&
        sub_vel(30) = subchannel_vel(i)
      if(point_distance(temp(:,i),sub_coord(:,34))<1.d-10)&
        sub_vel(34) = subchannel_vel(i)
      if(point_distance(temp(:,i),sub_coord(:,35))<1.d-10)&
        sub_vel(35) = subchannel_vel(i)
      if(point_distance(temp(:,i),sub_coord(:,36))<1.d-10)&
        sub_vel(36) = subchannel_vel(i)
    end do

    sub_vel(19) = sub_vel(24)
    sub_vel(20) = sub_vel(23)
    sub_vel(21) = sub_vel(22)

    sub_vel(25) = sub_vel(30)
    sub_vel(26) = sub_vel(29)
    sub_vel(27) = sub_vel(28)

    sub_vel(31) = sub_vel(36)
    sub_vel(32) = sub_vel(35)
    sub_vel(33) = sub_vel(34)

    do i=1,6
      sub_vel(i) = sub_vel(i+30)
    end do

    do i=7,12
      sub_vel(i) = sub_vel(i+18)
    end do

    do i=13,18
      sub_vel(i) = sub_vel(i+6)
    end do

    if(irangp<=0) then
      open(100, file='Exported_data_subchannel_velocity_Z.dat')
        write(100,*)'SubchannelID Subchannel_velocity_Z'
        do i=1,36
          write(100,*)i, sub_vel(i)
        end do
      close(100)
    end if

    contains
    
    function point_distance(coord_a, coord_b)

      implicit none

      double precision, dimension(:) :: coord_a, coord_b
      double precision :: point_distance

      integer :: i

      point_distance = 0.d0
      do i=1,3
        point_distance = (coord_a(i)-coord_b(i))**2+point_distance
      end do
      point_distance = sqrt(point_distance)

    end function point_distance

  end subroutine round_robin_subchannel_vel

!===============================================================================
!! This routine is used to check energy balance of the overlapping region in the
!! 7x7 locally ballooned rod bundle
  subroutine energy_balance(nvar, nscal, dt)
  
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

    implicit none

    ! Arguments

    integer :: nvar, nscal

    double precision :: dt(ncelet)

    ! Local variables

    integer          iel    , ifac   , ivar
    integer          iel1   , iel2   , ieltsm
    integer          iortho
    integer          inc    , iprev
    integer          iflmas , iflmab , ipccp
    integer          iscal
    integer          ilelt  , nlelt

    double precision xvara  , xvar
    double precision xbilan , xbilvl , xbilpa
    double precision xbilen , xbilso , xbildv
    double precision xfluxf
    double precision heated_wall_area
    double precision flumab , ctb1, ctb2

    integer, allocatable, dimension(:) :: lstelt

    double precision, allocatable, dimension(:,:) :: grad
    double precision, allocatable, dimension(:) :: treco
    double precision, allocatable, dimension(:) :: xcp
    double precision, dimension(:), pointer :: imasfl, bmasfl
    double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
    double precision, dimension(:), pointer ::  crom
    double precision, dimension(:), pointer :: cpro_visct, cpro_cp
    double precision, dimension(:), pointer :: cvar_scal, cvara_scal

    if (iscalt.le.0) return

    allocate(lstelt(max(ncel,nfac,nfabor)))
    allocate(xcp(ncel))

    if (ntcabs==ntpabs+1 .or. mod(ntcabs, 100)==0) then
      ! xbilvl: volume contribution of unsteady terms
      ! xbildv: volume contribution due to to term in div(rho u)
      ! xbilpa: contribution from heated/cooled walls
      ! xbilen: contribution from inlets
      ! xbilso: contribution from outlets
      ! xbilan: total balance

      xbilvl = 0.d0
      xbildv = 0.d0
      xbilpa = 0.d0
      xbilen = 0.d0
      xbilso = 0.d0
      xbilan = 0.d0
      heated_wall_area = 0.d0

      iscal = iscalt         ! temperature scalar number
      ivar =  isca(iscal)    ! temperature variable number

      call field_get_val_s(ivarfl(ivar), cvar_scal)
      call field_get_val_prev_s(ivarfl(ivar), cvara_scal)

      ! Physical quantity numbers
      call field_get_val_s(icrom, crom)
      call field_get_val_s(ivisct, cpro_visct)

      ! Pointers to the mass fluxes
      call field_get_key_int(ivarfl(ivar), kimasf, iflmas)
      call field_get_key_int(ivarfl(ivar), kbmasf, iflmab)
      call field_get_val_s(iflmas, imasfl)
      call field_get_val_s(iflmab, bmasfl)

      ! The balance is in Joule, so store the specific heat when dealing
      ! with the Temperature.

      if (icp.ge.0) call field_get_val_s(icp, cpro_cp)

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

      ! Boundary condition pointers for gradients and advection

      call field_get_coefa_s(ivarfl(ivar), coefap)
      call field_get_coefb_s(ivarfl(ivar), coefbp)

      ! Boundary condition pointers for diffusion

      call field_get_coefaf_s(ivarfl(ivar), cofafp)
      call field_get_coefbf_s(ivarfl(ivar), cofbfp)

      ! --> Compute value reconstructed at I' for boundary faces

      allocate(treco(nfabor))

      iortho = 0

      ! --> General case (for non-orthogonal meshes)

      if (iortho.eq.0) then

        ! Allocate a work array for the gradient calculation
        allocate(grad(3,ncelet))

        ! --- Compute temperature gradient

        inc = 1
        iprev = 0

        call field_gradient_scalar (ivarfl(ivar), iprev, inc, grad)

        ! - Compute reconstructed value in boundary cells

        do ifac = 1, nfabor
          iel = ifabor(ifac)
          treco(ifac) =   cvar_scal(iel)           &
                      + diipb(1,ifac)*grad(1,iel)  &
                      + diipb(2,ifac)*grad(2,iel)  &
                      + diipb(3,ifac)*grad(3,iel)
        enddo

        ! Free memory
        deallocate(grad)

      ! --> Case of orthogonal meshes

      else

        ! Compute reconstructed value
        ! (here, we assign the non-reconstructed value)

        do ifac = 1, nfabor
          iel = ifabor(ifac)
          treco(ifac) = cvar_scal(iel)
        enddo

      endif

      ! 2.1 Compute the balance at time step n

      ! --> Balance on interior volumes

      ! If it is variable, the density 'rom' has been computed at the beginning
      ! of the time step using the temperature from the previous time step.

      do iel = 1, ncel
        xvara = cvara_scal(iel)
        xvar  = cvar_scal(iel)
        xbilvl = xbilvl                              &
               + volume(iel) * xcp(iel) * crom(iel)  &
               * (xvara - xvar) / dt(iel)
      enddo

      ! --> Balance on all faces (interior and boundary), for div(rho u)
      !     ------------------------------------------------------------

      do ifac = 1, nfac

        iel1 = ifacel(1,ifac)
        if (iel1.le.ncel) then
          ctb1 = imasfl(ifac)*xcp(iel1)*cvar_scal(iel1)
        else
          ctb1 = 0d0
        endif

        iel2 = ifacel(2,ifac)
        if (iel2.le.ncel) then
          ctb2 = imasfl(ifac)*xcp(iel2)*cvar_scal(iel2)
        else
          ctb2 = 0d0
        endif

        xbildv =  xbildv + (ctb1 - ctb2)
      enddo

      do ifac = 1, nfabor
        iel = ifabor(ifac)
        xbildv = xbildv+xcp(iel)*bmasfl(ifac)* cvar_scal(iel)
      enddo

      ! --> Balance on boundary faces

      ! - Compute the contribution from heated/cooled walls 

      call getfbr('rod_walls', nlelt, lstelt)

      do ilelt = 1, nlelt

        ifac = lstelt(ilelt)
        iel  = ifabor(ifac)   ! associated boundary cell
        
        heated_wall_area = heated_wall_area + surfbn(ifac)
        
        ! Physical variables

        flumab = bmasfl(ifac)

        ! Contribution to flux from the current face
        ! (diffusion and convection fluxes)

        xfluxf = - surfbn(ifac)*(cofafp(ifac)+cofbfp(ifac)*treco(ifac))      &
               - flumab*xcp(iel)*(coefap(ifac)+coefbp(ifac)*treco(ifac))

        xbilpa = xbilpa + xfluxf

      enddo

      ! Contribution from inlet

      call getfbr('inlet', nlelt, lstelt)

      do ilelt = 1, nlelt

        ifac = lstelt(ilelt)
        iel  = ifabor(ifac)   ! associated boundary cell

        ! Physical variables

        flumab = bmasfl(ifac)

        ! Contribution to flux from the current face
        ! (diffusion and convection fluxes)

        xfluxf = - surfbn(ifac)*(cofafp(ifac)+cofbfp(ifac)*treco(ifac))      &
               - flumab*xcp(iel)*(coefap(ifac)+coefbp(ifac)*treco(ifac))

        xbilen = xbilen + xfluxf

      enddo

      ! Contribution from outlet

      call getfbr('outlet', nlelt, lstelt)

      do ilelt = 1, nlelt

        ifac = lstelt(ilelt)
        iel  = ifabor(ifac)   ! associated boundary cell

        ! Physical variables

        flumab = bmasfl(ifac)

        ! Contribution to flux from the current face
        ! (diffusion and convection fluxes)

        xfluxf = - surfbn(ifac)*(cofafp(ifac)+cofbfp(ifac)*treco(ifac))      &
               - flumab*xcp(iel)*(coefap(ifac)+coefbp(ifac)*treco(ifac))

        xbilso = xbilso + xfluxf

      enddo

      ! Now the work array for the temperature can be freed
      deallocate(treco)

      ! Sum of values on all ranks (parallel calculations)

      if (irangp.ge.0) then
        call parsom(xbilvl)
        call parsom(xbildv)
        call parsom(xbilpa)
        call parsom(xbilen)
        call parsom(xbilso)
        call parsom(heated_wall_area)
      endif

      ! --> Total balance

      ! We add the different contributions calculated above.

      xbilan = xbilvl + xbildv + xbilpa + xbilen + xbilso

      ! 2.3 Write the balance at time step n

      write (nfecra, 2000)                                               &
        ntcabs, heated_wall_area, xbilvl, xbildv, xbilpa,                &
        xbilen, xbilso, xbilan

      if(ntcabs==ntpabs+1)then
        open(100,file='energy_balance.dat',status='replace')
        write (100, 2010)                                                &
        ntcabs, heated_wall_area, xbilvl, xbildv, xbilpa,                &
        xbilen, xbilso, xbilan       
        close(100)
      else
        open(100,file='energy_balance.dat',status='old', position='append')
        write(100, 2020)                                                 &
        ntcabs, heated_wall_area, xbilvl, xbildv, xbilpa,                &
        xbilen, xbilso, xbilan  
        close(100)
      end if

2000 format                                                              &
      (/,                                                                &
      3X,'** Thermal balance **', /,                                     &
      3X,'   ---------------', /,                                        &
      '---', '------',                                                   &
      '------------------------------------------------------------', /, &
      'bt ','  Iter',                                                    &
      '   heated_w_a     Volume     Divergence    h/c_Wall',             &
      '     Inlet     Outlet    Total', /,                               &
      'bt ', i6, 7e12.4, /,                                              &
      '---','------',                                                    &
      '------------------------------------------------------------')
2010 format                                                              &
      (/,                                                                &
      3X,'** Thermal balance **', /,                                     &
      '---', '------',                                                   &
      '------------------------------------------------------------', /, &
      '   Iter',                                                         &
      '   heated_w_a     Volume     Divergence    h/c_Wall',             &
      '     Inlet     Outlet    Total', /,                               &
      i6, 7e12.4)

2020 format (i6, 7e12.4)

    endif ! End of test on first time step

    ! Deallocate the temporary array
    deallocate(lstelt)
    deallocate(xcp)

  end subroutine energy_balance


!===============================================================================
!! This routine is used to calculate the subchannel Nusselt number based on the
!! simulated thermal flux and boundary temperature (useful for CHT cases)
  subroutine subchannel_nusselt_check()

    use numvar
    use optcal
    use field
    use mesh
    use parall
    use cs_c_bindings

    character(len=80) :: chaine
    integer :: i, j, iel, ielt, ifac, nlelt
    integer, allocatable, dimension(:) :: lstelt
    double precision, dimension(:), pointer :: subchannel_ID, tbulk, twall, dh,&
                                               acond, tfluxwall, nusselt_corr
    double precision, dimension(:), allocatable :: atwall, htc_calc,           &
                                                   nusselt_calc, atfluxwall
    double precision, dimension(:), allocatable :: subchannel_wall_flux,       &
                                                   subchannel_wall_temp,       &
                                                   subchannel_wall_area,       &
                                                   htc, htc_corr, nusselt, tb
    integer, allocatable, dimension(:) :: nf, nc

    interface
      subroutine cal_boundary_thermal_flux()        &
        bind(C, name='cal_boundary_thermal_flux')
        use, intrinsic :: iso_c_binding
        implicit none
      end subroutine cal_boundary_thermal_flux
    end interface

    if (iscalt.le.0) return

!! Check number of wall-adjacent cells in each suchannel cell
    allocate(nf(n_cell), nc(n_cell))
    nf = 0
    nc = 0

    if (iscalt.le.0) return

!! Get all wall faces (only the SubChCFD domain)
    allocate(lstelt(nfabor))
    chaine = all_walls
    call getfbr(chaine, nlelt, lstelt)

!! Field arrays used

    call field_get_val_s_by_name('SubchannelID', subchannel_ID)
    call field_get_val_s_by_name('HydraulicDiameter', dh)
    call field_get_val_s_by_name('AverageThermalConductivity', acond)
    call field_get_val_s_by_name('WallTemperature0', twall)
!    call field_get_val_s_by_name('boundary_temperature', twall)
    call field_get_val_s_by_name('SubchannelBulkTemperature', tbulk)
    call field_get_val_s_by_name('WallThermalFlux', tfluxwall)
    call field_get_val_s_by_name('NusseltNumber', nusselt_corr)

    allocate(atwall(ncel))
    allocate(atfluxwall(ncel))
    allocate(nusselt_calc(ncel))
    allocate(htc_calc(ncel))

    allocate(subchannel_wall_area(n_cell))
    allocate(subchannel_wall_temp(n_cell))
    allocate(subchannel_wall_flux(n_cell))
    allocate(htc(n_cell))
    allocate(htc_corr(n_cell))
    allocate(nusselt(n_cell))
    allocate(tb(n_cell))

!! Calculate the wall thermal fluxes of each boundary face
!! Note: the calculated boundary heat flux here is slightly different from that
!! calculated in clptur.f90. There is one time step difference between the two
!! calculations, so the difference is due to convergence?
!    call cal_boundary_thermal_flux()

!! Disable solid cells in case of CHT
    if (fluid_solid) call cs_porous_model_set_has_disable_flag(1)

    subchannel_wall_area = 0.d0
    subchannel_wall_temp = 0.d0
    subchannel_wall_flux = 0.d0

    do ielt=1, nlelt
      ifac = lstelt(ielt)
      iel = ifabor(ifac)
      if(cell_is_active(iel)==0)cycle
      j = int(subchannel_ID(iel))
      subchannel_wall_area(j) = subchannel_wall_area(j)+surfbn(ifac)
      subchannel_wall_temp(j) = subchannel_wall_temp(j)+twall(ifac) &
                              * surfbn(ifac)
      subchannel_wall_flux(j) = subchannel_wall_flux(j)+tfluxwall(ifac) &
                              * surfbn(ifac)
      nf(j) = nf(j)+1
    end do

!! Re-enable solid cells in case of CHT
    if(fluid_solid) call cs_porous_model_set_has_disable_flag(0)

    if (irangp.ge.0) then
      call parrsm(n_cell, subchannel_wall_area)
      call parrsm(n_cell, subchannel_wall_temp)
      call parrsm(n_cell, subchannel_wall_flux)
      call parism(n_cell, nf)
    endif

    do iel=1, ncel
      j = int(subchannel_ID(iel))
      if(j>0) then
        atwall(iel) = subchannel_wall_temp(j)/subchannel_wall_area(j)
        atfluxwall(iel) = subchannel_wall_flux(j)/subchannel_wall_area(j)
        htc_calc(iel) = atfluxwall(iel)/(tbulk(iel)-atwall(iel))
        nusselt_calc(iel) = htc_calc(iel)*dh(iel)/acond(iel)
      end if
    end do

    htc = -1.d10
    htc_corr = -1.d10
    nusselt = -1.d10
    tb = -1.d10
    do i=1,n_cell
      do iel=1,ncel
        j = int(subchannel_ID(iel))
        if(j>0 .and. i==j)then
          htc(i) = htc_calc(iel)
          htc_corr(i) = nusselt_corr(iel)*acond(iel)/dh(iel)
          nusselt(i) = nusselt_calc(iel)
          tb(i) = tbulk(iel)
          cycle
        end if
      end do
    end do

    if (irangp.ge.0) then
      do i=1, n_cell
        call parmax(htc(i))
        call parmax(htc_corr(i))
        call parmax(nusselt(i))
        call parmax(tb(i))
      end do
    endif

! Output for plotting
    if(ntcabs==ntmabs .and. irangp<=0) then
      open(200,file='nusselt_calc.dat')
      write(200,'(5A16)') 'Subchannel ID', 'bulk temp.', 'Nusselt No.', &
                          'HTC sim', 'HTC corr'
      do i=1, n_cell
        write(200,'(I16, 4f16.5)') i, tb(i), nusselt(i), htc(i), htc_corr(i)
      end do
      close(200)
    end if

!! Deallocate temporary arrays
    deallocate(lstelt)
    deallocate(atwall)
    deallocate(atfluxwall)
    deallocate(htc_calc)
    deallocate(nusselt_calc)
    deallocate(subchannel_wall_temp)
    deallocate(subchannel_wall_flux)
    deallocate(subchannel_wall_area)
    deallocate(htc, htc_corr, nusselt)
    deallocate(nc,nf)

  end subroutine subchannel_nusselt_check


end module SubChCFD

!-------------------------------------------------------------------------------

!> (DOXYGEN_SHOULD_SKIP_THIS) \endcond
