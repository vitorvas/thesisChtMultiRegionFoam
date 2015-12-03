! file src/general_module.f90
!
! Copyright 2009-2011 Dalton Harvie (daltonh@unimelb.edu.au)
! 
! This file is part of arb finite volume solver, referred to as `arb'.
! 
! arb is a software package designed to solve arbitrary partial
! differential equations on unstructured meshes using the finite volume
! method.  Primarily it consists of fortran source code, perl source
! code and shell scripts.  arb replies on certain third party software
! to run, most notably the computer algebra system maxima
! <http://maxima.sourceforge.net/> which is released under the GNU GPL.
! 
! The copyright of arb is held by Dalton Harvie.
! 
! arb is released under the GNU GPL.  arb is free software: you can
! redistribute it and/or modify it under the terms of the GNU General
! Public License (version 3) as published by the Free Software Foundation.
! You should have received a copy of the GNU General Public License
! along with arb (see file license/gpl.txt after unpacking).  If not,
! see <http://www.gnu.org/licenses/>.
! 
! arb is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
! 
! For full details of arb's license see the license directory.
! 
! The current homepage for the arb finite volume solver project is
! <http://www.chemeng.unimelb.edu.au/people/staff/daltonh/downloads/arb>.
!
!-------------------------------------------------------------------------
! module file for variables and routines common to most of the code
module general_module

! statements specifying different data types and parameters
implicit none

! allow all subroutines and functions to be public

integer, parameter :: totaldimensions=3 ! this is maximum number of dimensions, possibly hardcode in future

! a generic kernel for calculating averages and derivatives
type kernel_type
  character(len=4) :: centring ! whether input values are cell (i), face (j) or node (k) centred
  integer, dimension(:), allocatable :: ijk ! array of cell/face/node indicies that are used in this kernel
  double precision, dimension(:), allocatable :: v ! value of kernel in 1-to-1 correspondance with ijk
end type kernel_type

! this type specifies details of each node (vertex)
type node_type
  integer :: type ! integer specifying whether node is within the domain (1) or on a boundary (2)
  double precision, dimension(totaldimensions) :: x ! location of node
  integer, dimension(:), allocatable :: jface ! array storing j indicies of surrounding faces
  integer, dimension(:), allocatable :: icell ! array storing i indicies of surrounding cells
end type node_type
  
! this type specifies details of each cell face
type face_type
  integer :: type ! integer specifying whether face is within the domain (1) or on a boundary (2)
  double precision, dimension(totaldimensions) :: x ! location of face centre
  double precision :: area ! area of face, (length in 1D, zero in 0D)
  double precision :: dx ! distance between adjacent cell centres
  double precision, dimension(totaldimensions,totaldimensions) :: norm ! orthogonal orientation vectors for face: first norm(:,1) points normal to face in direction from cell icell(1) -> icell(2), second norm(:,2) points from node(1)->node(2) on face (2d and 3d) (first tangent) and third norm(:,3) is normal to both (second tangent in 3d)
  integer, dimension(:), allocatable :: icell ! array storing i indicies of the 2 adjacent cells
  integer, dimension(:), allocatable :: knode ! array storing k indicies of surrounding nodes
  integer, dimension(:), allocatable :: region_list ! list of regions that the face is a member of
  integer :: dimensions ! number of dimensions that cell is (0, 1 or 2)
  integer :: gtype ! gmsh element type for element geometry
  type(kernel_type), dimension(0:6) :: kernel ! kernel(m) = kernel for the average (m=0) or derivative in the m'th coordinate direction - m>=4 is for normals to the face for norm(:,m-3)
end type face_type

! this type specifies details of each cell
type cell_type
  integer :: type ! integer specifying whether cell is within the domain (1) or on a boundary (2)
  double precision, dimension(totaldimensions) :: x ! location of cell centre
  double precision :: vol ! volume (area in 2D, length in 1D)
  integer, dimension(:), allocatable :: knode ! array storing k indicies of surrounding nodes
  integer, dimension(:), allocatable :: jface ! array storing j indicies of surrounding faces
  integer, dimension(:), allocatable :: icell ! array storing i indicies of surrounding cells
  integer, dimension(:), allocatable :: region_list ! list of regions that the cell is a member of
  integer :: dimensions ! number of dimensions that cell is (1, 2 or 3)
  integer :: gtype ! gmsh element type for element geometry
  type(kernel_type), dimension(0:4) :: kernel ! kernel(m) = kernel for the average from faces (m=0) or derivative in the m'th coordinate direction, or average from nodes (m=4)
end type cell_type
  
! type for links between regions
type region_link_type
  character(len=1000) :: from_region ! name of the region which we are linking from
  character(len=1000) :: to_region ! name of the region which we are linking to
  integer :: from_region_number ! number of the region which we are linking from
  integer :: to_region_number ! number of the region which we are linking to
  character(len=4) :: from_centring ! centring of the region we are linking from
  character(len=4) :: to_centring ! centring of the region we are linking to
  integer, dimension(:), allocatable :: to_ns ! ns index within to_region from ns index of element within from_region
end type region_link_type

! type for regions
type region_type
  character(len=1000) :: name ! name of the region
  character(len=1000) :: location ! string specifying location of region
  character(len=4) :: centring ! whether cell or face centred
  integer :: dimensions ! maximum dimensions of the elements within the region
  integer, dimension(:), allocatable :: ij ! array of cell or face i or j indicies that are within this region - dimension of this is number of elements in region
  integer, dimension(:), allocatable :: ns ! array that specifies data number from i or j index
end type region_type

! data type for any functions that ultimately depend on field data
type funk_type
  double precision :: v ! value of function
  double precision, dimension(:), allocatable :: dv ! value of derivative, in 1:1 with pp
  integer, dimension(:), allocatable :: pp ! phi variable which derivative is taken wrt
  integer :: ndv ! number of derivative elements that currently contain valid data
end type funk_type

! meta data type for general variables
type var_type
  character(len=1000) :: name ! name of the variable
  character(len=1000) :: units ! character string of the units
  double precision :: multiplier ! multiplier appended to units when interacting with outside world
  double precision :: magnitude ! an order of magnitude estimate of the variable, calculated from initial conditions and only for phi variables
  character(len=4) :: centring
  character(len=100) :: type ! variable type: constant, transient, unknown, derived, equation, output, condition, local
  integer :: relstep ! relative timestep, with relstep=0 being the current step
  integer :: deriv ! whether derivative is to be calculated for this variable
  character(len=6) :: rank ! specifies whether this is a component of a scalar, vector or tensor compound
  character(len=1000) :: region ! name of the region in which it is applied
  integer :: region_number ! number of the region in which it is applied
  character(len=1000) :: compound_name ! name of the compound variable of which this scalar is a component
  integer :: compound_number ! number of the compound variable of which this scalar is a component
  integer, dimension(:), allocatable :: component ! ordered list of the compound variable of which this scalar is a component
  character(len=100), dimension(:), allocatable :: options ! array of options for this var
  type(funk_type), dimension(:), allocatable :: funk ! an array of the data, which depending on centring, may contain 1 (none), itotal (cell) or jtotal (face) elements
  integer :: someloop ! if this variable is not to be stored (ie, local is on) then this is the someloop number to be used instead (otherwise if local is off someloop = 0)
end type var_type

! data type for var_lists
type var_list_type
! character(len=4) :: centring
! character(len=100) :: type ! variable type
  integer, dimension(:), allocatable :: list
end type var_list_type

! meta data type for all compound variables
type compound_type
  character(len=1000) :: name ! name of the variable
  character(len=1000) :: units ! character string of the units
  double precision :: multiplier ! multiplier appended to units when interacting with outside world
  character(len=4) :: centring
  character(len=100) :: type ! variable type: constant, transient, unknown, derived, equation, output, condition, local
  character(len=6) :: rank ! specifies whether this is a component of a scalar, vector or tensor compound
  integer :: relstep ! relative timestep, with relstep=0 being the current step
  integer :: deriv ! whether derivative is to be calculated for this variable
  character(len=1000) :: region ! name of the region in which it is applied
  integer :: region_number ! number of the region in which it is applied
  integer, dimension(:), allocatable :: component ! ordered list of the compound variable of which this scalar is a component
  character(len=100), dimension(:), allocatable :: options ! array of options for this compound
end type compound_type

! mesh arrays
type(node_type), dimension(:), allocatable :: node ! array of nodes
type(face_type), dimension(:), allocatable :: face ! array of faces
type(cell_type), dimension(:), allocatable :: cell ! array of cells

! now define actual data arrays which are all saved
! each variable has a corresponding data type
type(var_type), dimension(:), allocatable :: var ! var(m) = general variable of type m, itself containing array of funk data for each var
type(compound_type), dimension(:), allocatable :: compound ! list of compound (scalar, vector and tensor) general variables

! regions
type(region_type), dimension(:), allocatable :: region ! array of regions

! region links
type(region_link_type), dimension(:), allocatable :: region_link ! array of region links

real :: last_cpu_time
integer :: idomain, iboundary, itotal ! number of domain (type 1) and boundary (type 2) cells, also total
integer :: jdomain, jboundary, jtotal ! number of domain (type 1) and boundary (type 2) faces, also total
integer :: kdomain, kboundary, ktotal ! number of domain (type 1) and boundary (type 2) nodes, also total
integer :: ptotal ! number of equations
integer :: relstepmax ! maximum relstep value for a transient simulation
integer :: timestep = 0, newtstep = 0 ! timestep and newtonstep indicies
integer :: maximum_dimensions = 0 ! maximum dimensions of any region used in the simulation
type(funk_type), dimension(:), allocatable :: someloop ! funk structure for local updating use to store someloop values and derivatives
double precision, dimension(:), allocatable :: equation_magnitude ! expresses the magnitude of each equation funk for use in the residual calculation
character(len=100) :: constants_file
type(var_list_type), dimension(:), allocatable :: var_list ! array of var_lists, according to type and centring
double precision, parameter :: pi = 4.d0*atan(1.d0)
character(len=100), parameter :: indexformat = 'i8', floatformat='g18.6', realformat='g15.8', stringformat='a18' ! formating used for output throughout program
! reals apparently have about 7 decimal places and width has to be d+7 (ifort)
integer, parameter :: fwarn = 11, fdetail = 12, foutput = 13, fgmsh = 14, fconstants = 15, fconverge = 16, foutputstep = 17 ! various file handles

! define some string parameter arrays
character(len=100), dimension(8), parameter :: var_types = [ "constant   ", "transient  ", "unknown    ", "derived    ", &
  "equation   ", "output     ", "condition  ", "local      " ]
character(len=100), dimension(3), parameter :: stepoutput_options = ["stepoutput        ", "stepoutputnoupdate", &
  "nostepoutput      "]
character(len=100), dimension(2), parameter :: output_options = ["output  ", "nooutput"]
character(len=100), dimension(5), parameter :: output_gmesh_options = ["output            ", "centringoutput    ", &
  "meshoutput        ", "centringmeshoutput", "nooutput          "]
character(len=100), dimension(5), parameter :: input_gmesh_options = ["input            ", "centringinput    ", "meshinput        ", &
  "centringmeshinput", "noinput          "]
character(len=100), dimension(3), parameter :: data_options = ["elementdata           ", "elementnodedata       ", &
  "elementnodelimiteddata"]
character(len=100), dimension(2), parameter :: input_options = ["input            ", "noinput          "]

! code version details
real, parameter :: version = 0.30 ! current version
real, parameter :: minimum_version = 0.25 ! minimum version constants.in file that will still work with this version
character(len=100), parameter :: versionname = "useful ulmar"

! the following are default values for various parameters, some of which may be overwritten by input file options
logical :: transient_simulation = .false. ! whether simulation is transient or steady-state
integer :: timestepmax = huge(timestepmax) ! maximum number of timesteps performed
integer :: timestepout = huge(timestepout) ! maximum number of timesteps between output
integer :: newtstepmax = 1000 ! maximum number of steps performed by newton proceedure
double precision :: newtrestol = 1.d-10 ! tolerance that indicates convergence of the newton proceedure
double precision, parameter :: limitertolerance = 1.d-10 ! tolerance used when calculating advection gradient limiting - set to small positive number
double precision, parameter :: limitercontgrad = 2.d0 ! factor that determines the gradient of the continuous advection limiter - set ~> 1.15 and ~< 2
double precision, parameter :: unknown_mag_limit = 1.d+6 ! ratio between unknown magnitude and specified order of unknown that signals an error 
double precision, parameter :: eps_dv = 1.d-20
character(len=100) :: linear_solver = "DEFAULT" ! type of linear solver used: DEFAULT will choose optimal solver available
character(len=100) :: output_step_file = "default" ! whether to print output.step file or not: default|on, newtstep, timestep, output, final, off
logical, parameter :: output_timings = .false. ! whether to time processes and output results to screen (see subroutine time_process)
logical, parameter :: kernel_details_file = .false. ! print out a text file (output/kernel_details.txt) with all the kernel details
logical, parameter :: mesh_details_file = .false. ! print out a text file (output/mesh_details.txt) with all the mesh details
logical, parameter :: region_details_file = .false. ! print out a text file (output/region_details.txt) with all the region details
logical, parameter :: link_details_file = .false. ! print out a text file (output/link_details.txt) with all the link details
logical, parameter :: convergence_details_file = .false. ! write some convergence debugging data to output/convergence_details.txt

!----------------------------------------------------------------------------

! interface to make general array resizing procedure that can handle all standard data types
interface resize_array
  module procedure resize_double_precision_array, resize_real_array, resize_integer_array, resize_character_array
end interface resize_array

interface push_array
  module procedure push_double_precision_array, push_real_array, push_integer_array, push_character_array
end interface push_array

interface allocatable_size
  module procedure allocatable_double_precision_size, allocatable_real_size, allocatable_integer_size, allocatable_character_size
end interface allocatable_size

interface copy_array
  module procedure copy_double_precision_array, copy_real_array, copy_integer_array, copy_character_array
end interface copy_array

contains

!----------------------------------------------------------------------------

subroutine remove_comments(textline)

character(len=*) :: textline
integer :: cut

cut=scan(textline,'#') ! find position of comment character if it is there
if (cut>0) textline=textline(1:cut-1)//repeat(' ',len(textline)-cut) ! remove comment if there
!if (cut>0) textline=textline(1:cut-1) ! remove comment if there

end subroutine remove_comments

!----------------------------------------------------------------------------

function var_number_from_name(name)

integer :: var_number_from_name, m
character(len=*), intent(in) :: name

var_number_from_name = 0 ! default is an error code

do m=1,ubound(var,1)
  if (trim(var(m)%name) == trim(name)) then
    var_number_from_name = m
    return
  end if
end do

end function var_number_from_name

!----------------------------------------------------------------------------

function region_number_from_name(name,centring,location,dimensions,creatable,existing)

! find and checks data regarding region_number, or alternatively sets new one if
!  creatable is present and set to true

integer :: region_number_from_name
character(len=*), intent(in) :: name
character(len=1000), intent(in), optional :: location
character(len=4), intent(in), optional :: centring
integer, intent(in), optional :: dimensions
logical, intent(in), optional :: creatable
logical, intent(out), optional :: existing
type(region_type) :: default_element
integer :: n

region_number_from_name = 0 ! default is an error code
if (present(existing)) existing = .false.

! see if region already exists
if (allocated(region)) then
  do n=1,ubound(region,1)
    if (trim(region(n)%name) == trim(name)) then
      region_number_from_name = n
      exit
    end if
  end do
end if

! check centring of existing region is consistent if it was previously set
if (region_number_from_name /= 0) then
  if (present(existing)) existing = .true.
  if (present(centring)) then
    if (trim(region(region_number_from_name)%centring) /= ''.and. &
        centring /= region(region_number_from_name)%centring) then
      write(*,*) 'ERROR: existing region = '//trim(name)//' has different centring than previously specified'
      region_number_from_name = 0 ! set number to indicate error
      return
    end if
  end if
  if (present(dimensions)) then
    if (dimensions >= 0 .and. dimensions /= region(region_number_from_name)%dimensions) then
      write(*,*) 'ERROR: existing region = '//trim(name)//' has different dimensions than previously specified'
      region_number_from_name = 0 ! set number to indicate error
      return
    end if
  end if
end if

! to create or alter region data logical creatable must be present and set to true
if (.not.present(creatable)) return 
if (.not.creatable) return 

! create new region
if (region_number_from_name == 0) then
  default_element%name = name
  if (present(location)) then
    default_element%location = location
  else
    default_element%location = ''
  end if
  if (present(centring)) then
    default_element%centring = centring
  else
    default_element%centring = ''
  end if
  if (present(dimensions)) then
    default_element%dimensions = dimensions
  else
    default_element%dimensions = -1 ! this value indicates that the dimensions are not known
  end if
  call resize_region_array(new_element=default_element,change=1)
  region_number_from_name = ubound(region,1)
else
! check existing region and update any info if presently empty
  if (present(centring).and.trim(region(region_number_from_name)%centring) == '') then
    write(*,'(a)') 'INFO: updating region = '//trim(name)//' centring to '//trim(centring)
    region(region_number_from_name)%centring = centring
  end if
  if (present(location).and.trim(region(region_number_from_name)%location) == '') then
    write(*,'(a)') 'INFO: updating region = '//trim(name)//' location to '//trim(location)
    region(region_number_from_name)%location = location
  end if
  if (present(dimensions).and.region(region_number_from_name)%dimensions == -1) then
    write(*,'(a,i1)') 'INFO: updating region = '//trim(name)//' dimensions to ',dimensions
    region(region_number_from_name)%dimensions = dimensions
  end if
end if

end function region_number_from_name

!----------------------------------------------------------------------------

subroutine resize_node_array(new_element,change,new_size)

type(node_type), optional, intent(in) :: new_element
integer, intent(in), optional :: change, new_size
type(node_type), dimension(:), allocatable :: node_old
integer :: old_size, min_size, new_size_l, change_l, k

change_l = 1
if (present(change)) change_l = change

old_size = 0
if (allocated(node)) old_size = ubound(node,1)

if (present(new_size)) then
  if (present(change)) stop 'ERROR: both change and new_size specified in resize_node_array'
  new_size_l = new_size
else
  new_size_l = old_size + change_l
end if
min_size = min(old_size,new_size_l)

! if any data is to be required in the future then save it now
if (min_size>0) then
  allocate(node_old(min_size))
  do k=1,min_size ! also create same allocatable structure within
    if (allocated(node(k)%jface)) allocate(node_old(k)%jface(ubound(node(k)%jface,1)))
    if (allocated(node(k)%icell)) allocate(node_old(k)%icell(ubound(node(k)%icell,1)))
  end do
  node_old=node(1:min_size)
end if

if (old_size>0) deallocate(node)

if (new_size_l>0) then
  allocate(node(new_size_l))
  if (min_size>0) then
    do k=1,min_size ! copy over allocatable structure within
      if (allocated(node_old(k)%jface)) allocate(node(k)%jface(ubound(node_old(k)%jface,1)))
      if (allocated(node_old(k)%icell)) allocate(node(k)%icell(ubound(node_old(k)%icell,1)))
    end do
    node(1:min_size)=node_old
    deallocate(node_old)
  end if
  if (min_size<new_size_l.and.present(new_element)) then
    do k=min_size+1,new_size_l ! copy over allocatable structure within
      if (allocated(new_element%jface)) allocate(node(k)%jface(ubound(new_element%jface,1)))
      if (allocated(new_element%icell)) allocate(node(k)%icell(ubound(new_element%icell,1)))
    end do
    node(min_size+1:new_size_l) = new_element ! initialise new values
  end if
end if

!write(*,*) 'at end of push_node_array with element = '//trim(print_node(node(new_size)))

end subroutine resize_node_array

!----------------------------------------------------------------------------

subroutine resize_face_array(new_element,change,new_size)

type(face_type), optional, intent(in) :: new_element
integer, intent(in), optional :: change, new_size
type(face_type), dimension(:), allocatable :: face_old
integer :: old_size, min_size, new_size_l, change_l, j, l

change_l = 1
if (present(change)) change_l = change

old_size = 0
if (allocated(face)) old_size = ubound(face,1)

if (present(new_size)) then
  if (present(change)) stop 'ERROR: both change and new_size specified in resize_face_array'
  new_size_l = new_size
else
  new_size_l = old_size + change_l
end if
min_size = min(old_size,new_size_l)

! if any data is to be required in the future then save it now
if (min_size>0) then
  allocate(face_old(min_size))
  do j=1,min_size ! also create same allocatable structure within
    if (allocated(face(j)%knode)) allocate(face_old(j)%knode(ubound(face(j)%knode,1)))
    if (allocated(face(j)%icell)) allocate(face_old(j)%icell(ubound(face(j)%icell,1)))
    do l = lbound(face(j)%kernel,1), ubound(face(j)%kernel,1)
      if (allocated(face(j)%kernel(l)%ijk)) allocate(face_old(j)%kernel(l)%ijk(ubound(face(j)%kernel(l)%ijk,1)))
      if (allocated(face(j)%kernel(l)%v)) allocate(face_old(j)%kernel(l)%v(ubound(face(j)%kernel(l)%v,1)))
    end do
  end do
  face_old=face(1:min_size)
end if

if (old_size>0) deallocate(face)

if (new_size_l>0) then
  allocate(face(new_size_l))
  if (min_size>0) then
    do j=1,min_size ! copy over allocatable structure within
      if (allocated(face_old(j)%knode)) allocate(face(j)%knode(ubound(face_old(j)%knode,1)))
      if (allocated(face_old(j)%icell)) allocate(face(j)%icell(ubound(face_old(j)%icell,1)))
      do l = lbound(face_old(j)%kernel,1), ubound(face_old(j)%kernel,1)
        if (allocated(face_old(j)%kernel(l)%ijk)) allocate(face(j)%kernel(l)%ijk(ubound(face_old(j)%kernel(l)%ijk,1)))
        if (allocated(face_old(j)%kernel(l)%v)) allocate(face(j)%kernel(l)%v(ubound(face_old(j)%kernel(l)%v,1)))
      end do
    end do
    face(1:min_size)=face_old
    deallocate(face_old)
  end if
  if (min_size<new_size_l.and.present(new_element)) then
    do j=min_size+1,new_size_l ! copy over allocatable structure within
      if (allocated(new_element%knode)) allocate(face(j)%knode(ubound(new_element%knode,1)))
      if (allocated(new_element%icell)) allocate(face(j)%icell(ubound(new_element%icell,1)))
      do l = lbound(new_element%kernel,1), ubound(new_element%kernel,1)
        if (allocated(new_element%kernel(l)%ijk)) allocate(face(j)%kernel(l)%ijk(ubound(new_element%kernel(l)%ijk,1)))
        if (allocated(new_element%kernel(l)%v)) allocate(face(j)%kernel(l)%v(ubound(new_element%kernel(l)%v,1)))
      end do
    end do
    face(min_size+1:new_size_l) = new_element ! initialise new values
  end if
end if

!write(*,*) 'at end of push_face_array with element = '//trim(print_face(face(new_size)))

end subroutine resize_face_array

!----------------------------------------------------------------------------

subroutine resize_cell_array(new_element,change,new_size)

type(cell_type), optional, intent(in) :: new_element
integer, intent(in), optional :: change, new_size
type(cell_type), dimension(:), allocatable :: cell_old
integer :: old_size, min_size, new_size_l, change_l, i, l

change_l = 1
if (present(change)) change_l = change

old_size = 0
if (allocated(cell)) old_size = ubound(cell,1)

if (present(new_size)) then
  if (present(change)) stop 'ERROR: both change and new_size specified in resize_cell_array'
  new_size_l = new_size
else
  new_size_l = old_size + change_l
end if
min_size = min(old_size,new_size_l)

! if any data is to be required in the future then save it now
if (min_size>0) then
  allocate(cell_old(min_size))
  do i=1,min_size ! also create same allocatable structure within
    if (allocated(cell(i)%knode)) allocate(cell_old(i)%knode(ubound(cell(i)%knode,1)))
    if (allocated(cell(i)%jface)) allocate(cell_old(i)%jface(ubound(cell(i)%jface,1)))
    if (allocated(cell(i)%icell)) allocate(cell_old(i)%icell(ubound(cell(i)%icell,1)))
    do l = lbound(cell(i)%kernel,1), ubound(cell(i)%kernel,1)
      if (allocated(cell(i)%kernel(l)%ijk)) allocate(cell_old(i)%kernel(l)%ijk(ubound(cell(i)%kernel(l)%ijk,1)))
      if (allocated(cell(i)%kernel(l)%v)) allocate(cell_old(i)%kernel(l)%v(ubound(cell(i)%kernel(l)%v,1)))
    end do
  end do
  cell_old=cell(1:min_size)
end if

if (old_size>0) deallocate(cell)

if (new_size_l>0) then
  allocate(cell(new_size_l))
  if (min_size>0) then
    do i=1,min_size ! copy over allocatable structure within
      if (allocated(cell_old(i)%knode)) allocate(cell(i)%knode(ubound(cell_old(i)%knode,1)))
      if (allocated(cell_old(i)%jface)) allocate(cell(i)%jface(ubound(cell_old(i)%jface,1)))
      if (allocated(cell_old(i)%icell)) allocate(cell(i)%icell(ubound(cell_old(i)%icell,1)))
      do l = lbound(cell_old(i)%kernel,1), ubound(cell_old(i)%kernel,1)
        if (allocated(cell_old(i)%kernel(l)%ijk)) allocate(cell(i)%kernel(l)%ijk(ubound(cell_old(i)%kernel(l)%ijk,1)))
        if (allocated(cell_old(i)%kernel(l)%v)) allocate(cell(i)%kernel(l)%v(ubound(cell_old(i)%kernel(l)%v,1)))
      end do
    end do
    cell(1:min_size)=cell_old
    deallocate(cell_old)
  end if
  if (min_size<new_size_l.and.present(new_element)) then
    do i=min_size+1,new_size_l ! copy over allocatable structure within
      if (allocated(new_element%knode)) allocate(cell(i)%knode(ubound(new_element%knode,1)))
      if (allocated(new_element%jface)) allocate(cell(i)%jface(ubound(new_element%jface,1)))
      if (allocated(new_element%icell)) allocate(cell(i)%icell(ubound(new_element%icell,1)))
      do l = lbound(new_element%kernel,1), ubound(new_element%kernel,1)
        if (allocated(new_element%kernel(l)%ijk)) allocate(cell(i)%kernel(l)%ijk(ubound(new_element%kernel(l)%ijk,1)))
        if (allocated(new_element%kernel(l)%v)) allocate(cell(i)%kernel(l)%v(ubound(new_element%kernel(l)%v,1)))
      end do
    end do
    cell(min_size+1:new_size_l) = new_element ! initialise new values
  end if
end if

end subroutine resize_cell_array

!----------------------------------------------------------------------------

subroutine resize_region_array(new_element,change,new_size)

type(region_type), optional, intent(in) :: new_element
integer, intent(in), optional :: change, new_size
type(region_type), dimension(:), allocatable :: region_old
integer :: old_size, min_size, new_size_l, change_l, i

change_l = 1
if (present(change)) change_l = change

old_size = 0
if (allocated(region)) old_size = ubound(region,1)

if (present(new_size)) then
  if (present(change)) stop 'ERROR: both change and new_size specified in resize_region_array'
  new_size_l = new_size
else
  new_size_l = old_size + change_l
end if
min_size = min(old_size,new_size_l)

! if any data is to be required in the future then save it now
if (min_size>0) then
  allocate(region_old(min_size))
  do i=1,min_size ! also create same allocatable structure within
    if (allocated(region(i)%ij)) allocate(region_old(i)%ij(ubound(region(i)%ij,1)))
  end do
  region_old=region(1:min_size)
end if

if (old_size>0) deallocate(region)

if (new_size_l>0) then
  allocate(region(new_size_l))
  if (min_size>0) then
    do i=1,min_size ! copy over allocatable structure within
      if (allocated(region_old(i)%ij)) allocate(region(i)%ij(ubound(region_old(i)%ij,1)))
    end do
    region(1:min_size)=region_old
    deallocate(region_old)
  end if
  if (min_size<new_size_l.and.present(new_element)) then
    do i=min_size+1,new_size_l ! copy over allocatable structure within
      if (allocated(new_element%ij)) allocate(region(i)%ij(ubound(new_element%ij,1)))
    end do
    region(min_size+1:new_size_l) = new_element ! initialise new values
  end if
end if

end subroutine resize_region_array

!----------------------------------------------------------------------------

subroutine resize_sparse_matrix(keep_data,aa,icn,irn,change,new_size)

! here we change the size of a sparse matrix (by change or too new_size) while 
!  possibly (by default) keeping the data

double precision, dimension(:), allocatable :: aa
integer, dimension(:), allocatable :: icn, irn
integer, optional, intent(in) :: change, new_size
logical, optional, intent(in) :: keep_data
logical :: keep_data_l

! by default keep the data in the arrays
keep_data_l = .true.
if (present(keep_data)) then
  if (.not.keep_data) keep_data_l = .false.
end if

if (present(change)) then
  call resize_double_precision_array(keep_data=keep_data_l,array=aa,change=change)
  call resize_integer_array(keep_data=keep_data_l,array=icn,change=change)
  call resize_integer_array(keep_data=keep_data_l,array=irn,change=change)
else if (present(new_size)) then
  call resize_double_precision_array(keep_data=keep_data_l,array=aa,new_size=new_size)
  call resize_integer_array(keep_data=keep_data_l,array=icn,new_size=new_size)
  call resize_integer_array(keep_data=keep_data_l,array=irn,new_size=new_size)
end if

end subroutine resize_sparse_matrix

!----------------------------------------------------------------------------

subroutine resize_double_precision_array(keep_data,array,change,new_size)

! here we change the size of an array (by change) while maintaining its data

double precision, dimension(:), allocatable :: array, array_store
integer, intent(in), optional :: change, new_size
integer :: change_l, old_size, new_size_l, min_size
logical, optional :: keep_data


change_l = 1
if (present(change)) change_l = change

old_size = 0
if (allocated(array)) old_size = ubound(array,1)

if (present(new_size)) then
  if (present(change)) stop 'ERROR: both change and new_size specified in resize_double_precision_array'
  new_size_l = new_size
else
  new_size_l = old_size + change_l
end if

min_size = min(old_size,new_size_l)

! if not keeping data (default is to keep data) then...
if (present(keep_data)) then
  if (.not.keep_data) min_size = 0
end if

! if any data is to be required in the future then save it now
if (min_size>0) then
  allocate(array_store(min_size))
  array_store = array(1:min_size) !bit pedantic - should be OK without indicies
end if

!if (old_size>0) deallocate(array)
if (allocated(array)) deallocate(array)

if (new_size_l>0) then
  allocate(array(new_size_l))
  if (min_size>0) then
    array(1:min_size) = array_store !bit pedantic - should be OK without indicies
    deallocate(array_store)
  end if
  if (min_size<new_size_l) array(min_size+1:new_size_l) = 0.d0 ! initialise new values
end if

end subroutine resize_double_precision_array

!----------------------------------------------------------------------------

subroutine resize_integer_array(keep_data,array,change,new_size)

! here we change the size of an array (by change or to new_size) while maintaining its data

integer, dimension(:), allocatable :: array, array_store
integer, intent(in), optional :: change, new_size
integer :: change_l, old_size, new_size_l, min_size
logical, optional :: keep_data

change_l = 1
if (present(change)) change_l = change

old_size = 0
if (allocated(array)) old_size = ubound(array,1)

if (present(new_size)) then
  if (present(change)) stop 'ERROR: both change and new_size specified in resize_integer_array'
  new_size_l = new_size
else
  new_size_l = old_size + change_l
end if

min_size = min(old_size,new_size_l)

! if not keeping data (default is to keep data) then...
if (present(keep_data)) then
  if (.not.keep_data) min_size = 0
end if

! if any data is to be required in the future then save it now
if (min_size>0) then
  allocate(array_store(min_size))
  array_store = array(1:min_size) !bit pedantic - should be OK without indicies
end if

!if (old_size>0) deallocate(array)
if (allocated(array)) deallocate(array)

if (new_size_l>0) then
  allocate(array(new_size_l))
  if (min_size>0) then
    array(1:min_size) = array_store !bit pedantic - should be OK without indicies
    deallocate(array_store)
  end if
  if (min_size<new_size_l) array(min_size+1:new_size_l) = 0 ! initialise new values
end if

end subroutine resize_integer_array

!----------------------------------------------------------------------------

subroutine resize_character_array(keep_data,array,change,new_size)

! here we change the size of an array (by change) while maintaining its data

character(len=*), dimension(:), allocatable :: array
character(len=len(array)), dimension(:), allocatable :: array_store
integer, intent(in), optional :: change, new_size
integer :: change_l, old_size, new_size_l, min_size
logical, optional :: keep_data

change_l = 1
if (present(change)) change_l = change

old_size = 0
if (allocated(array)) old_size = ubound(array,1)

if (present(new_size)) then
  if (present(change)) stop 'ERROR: both change and new_size specified in resize_character_array'
  new_size_l = new_size
else
  new_size_l = old_size + change_l
end if

min_size = min(old_size,new_size_l)

! if not keeping data (default is to keep data) then...
if (present(keep_data)) then
  if (.not.keep_data) min_size = 0
end if

! if any data is to be required in the future then save it now
if (min_size>0) then
  allocate(array_store(min_size))
  array_store = array(1:min_size) !bit pedantic - should be OK without indicies
end if

!if (old_size>0) deallocate(array)
if (allocated(array)) deallocate(array)

if (new_size_l>0) then
  allocate(array(new_size_l))
  if (min_size>0) then
    array(1:min_size) = array_store !bit pedantic - should be OK without indicies
    deallocate(array_store)
  end if
  if (min_size<new_size_l) array(min_size+1:new_size_l) = '' ! initialise new values
end if

end subroutine resize_character_array

!----------------------------------------------------------------------------

subroutine resize_real_array(keep_data,array,change,new_size)

! here we change the size of an array (by change) while maintaining its data

real, dimension(:), allocatable :: array, array_store
integer, intent(in), optional :: change, new_size
integer :: change_l, old_size, new_size_l, min_size
logical, optional :: keep_data

change_l = 1
if (present(change)) change_l = change

old_size = 0
if (allocated(array)) old_size = ubound(array,1)

if (present(new_size)) then
  if (present(change)) stop 'ERROR: both change and new_size specified in resize_real_array'
  new_size_l = new_size
else
  new_size_l = old_size + change_l
end if

min_size = min(old_size,new_size_l)

! if not keeping data (default is to keep data) then...
if (present(keep_data)) then
  if (.not.keep_data) min_size = 0
end if

! if any data is to be required in the future then save it now
if (min_size>0) then
  allocate(array_store(min_size))
  array_store = array(1:min_size) !bit pedantic - should be OK without indicies
end if

!if (old_size>0) deallocate(array)
if (allocated(array)) deallocate(array)

if (new_size_l>0) then
  allocate(array(new_size_l))
  if (min_size>0) then
    array(1:min_size) = array_store !bit pedantic - should be OK without indicies
    deallocate(array_store)
  end if
  if (min_size<new_size_l) array(min_size+1:new_size_l) = 0.e0 ! initialise new values
end if

end subroutine resize_real_array

!----------------------------------------------------------------------------

subroutine push_character_array(array,new_element,reverse)

! add element to an array, by default at the end, or with reverse to the start
character(len=*), dimension(:), allocatable :: array
character(len=*) :: new_element
logical, optional :: reverse
integer :: n

call resize_character_array(array=array,change=1)
array(ubound(array,1))=new_element

! if reverse is specified then add the element to the beginning of the array rather than the start
if (present(reverse)) then
  if (reverse) then
    do n = ubound(array,1)-1, 1, -1
      array(n+1) = array(n)
    end do
    array(1) = new_element
  end if
end if

end subroutine push_character_array

!----------------------------------------------------------------------------

subroutine unshift_character_array(array,new_element)

! add element to the start of an array
character(len=*), dimension(:), allocatable :: array
character(len=*) :: new_element

call push_character_array(array=array,new_element=new_element,reverse=.true.)

end subroutine unshift_character_array

!----------------------------------------------------------------------------

subroutine push_integer_array(array,new_element)

integer, dimension(:), allocatable :: array
integer :: new_element

call resize_integer_array(array=array,change=1)
array(ubound(array,1))=new_element

end subroutine push_integer_array

!----------------------------------------------------------------------------

subroutine push_double_precision_array(array,new_element)

double precision, dimension(:), allocatable :: array
double precision :: new_element

call resize_double_precision_array(array=array,change=1)
array(ubound(array,1))=new_element

end subroutine push_double_precision_array

!----------------------------------------------------------------------------

subroutine push_real_array(array,new_element)

real, dimension(:), allocatable :: array
real :: new_element

call resize_real_array(array=array,change=1)
array(ubound(array,1))=new_element

end subroutine push_real_array

!----------------------------------------------------------------------------

function allocatable_integer_size(array)

! find size of integer array that may or may not be allocated

integer, dimension(:), allocatable :: array
integer :: allocatable_integer_size

if (allocated(array)) then
  allocatable_integer_size = size(array)
else
  allocatable_integer_size = 0
end if

end function allocatable_integer_size

!----------------------------------------------------------------------------

function allocatable_double_precision_size(array)

! find size of double precision array that may or may not be allocated

double precision, dimension(:), allocatable :: array
integer :: allocatable_double_precision_size

if (allocated(array)) then
  allocatable_double_precision_size = size(array)
else
  allocatable_double_precision_size = 0
end if

end function allocatable_double_precision_size

!----------------------------------------------------------------------------

function allocatable_real_size(array)

! find size of real array that may or may not be allocated

real, dimension(:), allocatable :: array
integer :: allocatable_real_size

if (allocated(array)) then
  allocatable_real_size = size(array)
else
  allocatable_real_size = 0
end if

end function allocatable_real_size

!----------------------------------------------------------------------------

function allocatable_character_size(array)

! find size of character array that may or may not be allocated

character(len=*), dimension(:), allocatable :: array
integer :: allocatable_character_size

if (allocated(array)) then
  allocatable_character_size = size(array)
else
  allocatable_character_size = 0
end if

end function allocatable_character_size

!----------------------------------------------------------------------------

function print_node(element)

character(len=1000) :: print_node
type(node_type) :: element
integer :: n, error, njface, nicell
character(len=1000) :: formatline

njface = 0
if (allocated(element%jface)) njface = ubound(element%jface(:),1)
nicell = 0
if (allocated(element%icell)) nicell = ubound(element%icell(:),1)

formatline = '(a,i1,a'//repeat(',a,g12.5',totaldimensions)// &
  ',a'//repeat(',a,'//trim(indexformat),njface)// &
  ',a'//repeat(',a,'//trim(indexformat),nicell)//')'

write(print_node,fmt=formatline,iostat=error) 'type = ',element%type, &
  ': x =',(' ',element%x(n),n=1,totaldimensions), &
  ': jface =',(' ',element%jface(n),n=1,njface), &
  ': icell =',(' ',element%icell(n),n=1,nicell)

if (error /= 0) print_node = 'ERROR: problem in print_node: formatline = '//trim(formatline)

end function print_node

!----------------------------------------------------------------------------

function print_face(element)

character(len=10000) :: print_face
type(face_type) :: element
integer :: n, error, nknode, nicell, nregions, l
character(len=10000) :: formatline

nknode = 0
if (allocated(element%knode)) nknode = ubound(element%knode(:),1)
nicell = 0
if (allocated(element%icell)) nicell = ubound(element%icell(:),1)
nregions = 0
if (allocated(element%region_list)) nregions = ubound(element%region_list,1)

formatline = '(a,i1,a,i1,a,i2,a'//repeat(',a,g12.5',totaldimensions)// &
  ',a,g12.5,a,g12.5'// &
  repeat(',a,i1,a'//repeat(',a,g12.5',totaldimensions),totaldimensions)// &
  ',a'//repeat(',a,'//trim(indexformat),nknode)// &
  ',a'//repeat(',a,'//trim(indexformat),nicell)// &
  ',a'//repeat(',a,i3',nregions)//')'

write(print_face,fmt=formatline,iostat=error) 'type = ',element%type, &
  ': dimensions = ',element%dimensions, &
  ': gtype = ',element%gtype, &
  ': x =',(' ',element%x(n),n=1,totaldimensions), &
  ': area = ',element%area,': dx = ',element%dx, &
  (': norm(',l,') =',(' ',element%norm(n,l),n=1,totaldimensions),l=1,totaldimensions), &
  ': knode =',(' ',element%knode(n),n=1,nknode), &
  ': icell =',(' ',element%icell(n),n=1,nicell), &
  ': region_list =',(' ',element%region_list(n),n=1,nregions)

if (error /= 0) print_face = 'ERROR: problem in print_face: formatline = '//trim(formatline)

end function print_face

!----------------------------------------------------------------------------

function print_cell(element)

character(len=10000) :: print_cell
type(cell_type) :: element
integer :: n, njface, error, nknode, nicell, nregions
character(len=10000) :: formatline

nknode = 0
if (allocated(element%knode)) nknode = ubound(element%knode,1)
njface = 0
if (allocated(element%jface)) njface = ubound(element%jface,1)
nicell = 0
if (allocated(element%icell)) nicell = ubound(element%icell,1)
nregions = 0
if (allocated(element%region_list)) nregions = ubound(element%region_list,1)

formatline = '(a,i1,a,i1,a,i2,a'//repeat(',a,g12.5',totaldimensions)// &
  ',a,g12.5,a,i1'// &
  ',a'//repeat(',a,'//trim(indexformat),nknode)// &
  ',a'//repeat(',a,'//trim(indexformat),njface)// &
  ',a'//repeat(',a,'//trim(indexformat),nicell)// &
  ',a'//repeat(',a,i3',nregions)//')'

write(print_cell,fmt=formatline,iostat=error) 'type = ',element%type, &
  ': dimensions = ',element%dimensions, &
  ': gtype = ',element%gtype, &
  ': x =',(' ',element%x(n),n=1,totaldimensions), &
  ': vol = ', element%vol, ': njface = ',njface, &
  ': knode =',(' ',element%knode(n),n=1,nknode), &
  ': jface =',(' ',element%jface(n),n=1,njface), &
  ': icell =',(' ',element%icell(n),n=1,nicell), &
  ': region_list =',(' ',element%region_list(n),n=1,nregions)

if (error /= 0) print_cell = 'ERROR: problem in print_cell: formatline = '//trim(formatline)

end function print_cell

!----------------------------------------------------------------------------

subroutine add_to_dv(funke,deriv,funka)

! derivative info from dv and pp are added to general variable container
! now elements are placed in ascending pp order with no duplicates

double precision :: deriv ! factor by which dv data is multiplied as it is added
type(funk_type) :: funke ! existing funk container that will have funka dv data added to it
type(funk_type) :: funka ! container that has dv data to add to the funke
type(funk_type) :: funkt ! funk container which is used to assemble combined funk
integer na, ne, n
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'function add_to_dv'

if (debug) write(*,*) 'deriv = ',deriv
if (deriv == 0.d0) return
if (debug) write(*,*) 'funka%ndv = ',funka%ndv
if (funka%ndv == 0) return

! make sure funkt is large enough - base size on no duplication of elements
if (.not.allocated(funkt%pp)) then
  if (debug) write(*,*) 'allocating funkt with size ',funka%ndv+funke%ndv
  allocate(funkt%pp(funka%ndv+funke%ndv),funkt%dv(funka%ndv+funke%ndv))
else if (ubound(funkt%pp,1) < funka%ndv+funke%ndv) then
  deallocate(funkt%pp,funkt%dv)
  if (debug) write(*,*) 'reallocating funkt with size ',funka%ndv+funke%ndv
  allocate(funkt%pp(funka%ndv+funke%ndv),funkt%dv(funka%ndv+funke%ndv))
end if

if (debug) then
  write(*,*) 'adding funka: funka%ndv = ',funka%ndv
  if (funka%ndv > 0) then
    write(*,'(100(a,i8))') (' ',funka%pp(n),n=1,funka%ndv)
    write(*,'(100(a,g9.2))') (' ',deriv*funka%dv(n),n=1,funka%ndv)
  end if
  write(*,*) 'existing funke: funke%ndv = ',funke%ndv
  if (funke%ndv > 0) then
    write(*,'(100(a,i8))') (' ',funke%pp(n),n=1,funke%ndv)
    write(*,'(100(a,g9.2))') (' ',funke%dv(n),n=1,funke%ndv)
  end if
! if (allocated(funke%pp)) then
!   write(*,'(100(a,i8))') (' ',funke%pp(n),n=1,funke%ndv)
!   write(*,'(100(a,g9.2))') (' ',funke%dv(n),n=1,funke%ndv)
! else
!   write(*,*) 'existing funke not allocated'
! end if
end if

! cycle through elements of both funke and funka, adding them to funkt in ascending order

ne = 1 ! next read position in funke
na = 1 ! next read position in funka
funkt%ndv = 0 ! last write position in funkt

do
  if (na > funka%ndv) then ! we've come to the end of funka - add on remaining funke elements
    funkt%pp(funkt%ndv+1:funkt%ndv+1+funke%ndv-ne) = funke%pp(ne:funke%ndv)
    funkt%dv(funkt%ndv+1:funkt%ndv+1+funke%ndv-ne) = funke%dv(ne:funke%ndv)
    funkt%ndv = funkt%ndv+1+funke%ndv-ne
!   ne = funke%ndv + 1
    if (debug) write(*,*) 'adding on final funke elements: na = ',na
    exit
  else if (ne > funke%ndv) then ! we've come to the end of funke - add on remaining funka elements
    funkt%pp(funkt%ndv+1:funkt%ndv+1+funka%ndv-na) = funka%pp(na:funka%ndv)
    funkt%dv(funkt%ndv+1:funkt%ndv+1+funka%ndv-na) = deriv*funka%dv(na:funka%ndv)
    funkt%ndv = funkt%ndv+1+funka%ndv-na
!   na = funka%ndv + 1
    if (debug) write(*,*) 'adding on final funka elements: ne = ',ne
    exit
  else if (funka%pp(na) < funke%pp(ne)) then ! funka comes next
    funkt%ndv = funkt%ndv + 1
    funkt%pp(funkt%ndv:funkt%ndv) = funka%pp(na)
    funkt%dv(funkt%ndv:funkt%ndv) = deriv*funka%dv(na)
    na = na + 1
  else if (funke%pp(ne) < funka%pp(na)) then ! funke comes next
    funkt%ndv = funkt%ndv + 1
    funkt%pp(funkt%ndv:funkt%ndv) = funke%pp(ne)
    funkt%dv(funkt%ndv:funkt%ndv) = funke%dv(ne)
    ne = ne + 1
  else ! funke and funka have same pp element next
    funkt%ndv = funkt%ndv + 1
    funkt%pp(funkt%ndv:funkt%ndv) = funke%pp(ne)
    funkt%dv(funkt%ndv:funkt%ndv) = funke%dv(ne) + deriv*funka%dv(na)
    na = na + 1
    ne = ne + 1
    if (debug) write(*,*) 'duplicate elements: ne = ',ne,': na = ',na
  end if
end do
    
if (debug) write(*,*) 'all elements added to funkt: funkt%ndv = ',funkt%ndv

! now copy funtk back to funke
if (.not.allocated(funke%pp)) then
  if (debug) write(*,*) 'allocating funke with size ',funkt%ndv
  allocate(funke%pp(funkt%ndv),funke%dv(funkt%ndv))
else if (ubound(funke%pp,1) < funkt%ndv) then
  deallocate(funke%pp,funke%dv)
  if (debug) write(*,*) 'reallocating funke with size ',funkt%ndv
  allocate(funke%pp(funkt%ndv),funke%dv(funkt%ndv))
end if

if (.false.) then ! copy over wholeus boleus
  funke%ndv = funkt%ndv 
  funke%pp(1:funke%ndv) = funkt%pp(1:funke%ndv)
  funke%dv(1:funke%ndv) = funkt%dv(1:funke%ndv)
else ! copy over elements whose magnitude is larger than eps_dv
  funke%ndv = 0
  do n = 1, funkt%ndv
    if (abs(funkt%dv(n)) > eps_dv) then
      funke%ndv = funke%ndv + 1
      funke%pp(funke%ndv) = funkt%pp(n)
      funke%dv(funke%ndv) = funkt%dv(n)
    end if
  end do
end if

if (debug) then
  write(*,*) 'final funke: funke%ndv = ',funke%ndv
  write(*,'(100(a,i8))') (' ',funke%pp(n),n=1,funke%ndv)
  write(*,'(100(a,g9.2))') (' ',funke%dv(n),n=1,funke%ndv)
end if

if (debug) write(*,'(a/80(1h-))') 'function add_to_dv'

end subroutine add_to_dv

!----------------------------------------------------------------------------

subroutine multiply_dv(funke,deriv)

! derivatives dv in funke are multiplied by scalar deriv

double precision :: deriv ! factor by which dv data is multiplied
type(funk_type) :: funke ! existing funk container
integer :: n
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'function multiply_dv'

if (debug) then
  write(*,*) 'initial funke: funke%ndv = ',funke%ndv
  write(*,'(100(a,i8))') (' ',funke%pp(n),n=1,funke%ndv)
  write(*,'(100(a,g9.2))') (' ',funke%dv(n),n=1,funke%ndv)
  write(*,*) 'deriv = ',deriv
end if

if (allocated(funke%dv).and.funke%ndv > 0) funke%dv(1:funke%ndv) = funke%dv(1:funke%ndv)*deriv

if (debug) then
  write(*,*) 'final funke: funke%ndv = ',funke%ndv
  write(*,'(100(a,i8))') (' ',funke%pp(n),n=1,funke%ndv)
  write(*,'(100(a,g9.2))') (' ',funke%dv(n),n=1,funke%ndv)
end if

if (debug) write(*,'(a/80(1h-))') 'function multiply_dv'

end subroutine multiply_dv

!----------------------------------------------------------------------------

subroutine reset_funk(funkl)

! info in funk structure is reset

type(funk_type) :: funkl

funkl%v = 0.d0
funkl%ndv = 0
!if (allocated(funkl%dv)) funkl%dv(:) = 0.d0 ! not strictly necessary
!if (allocated(funkl%pp)) funkl%pp(:) = 0 ! not strictly necessary

end subroutine reset_funk

!----------------------------------------------------------------------------

function print_funk(funkl)

character(len=1000) :: print_funk, print_funk_saved
type(funk_type) :: funkl ! local funk variable
integer :: error, nn
character(len=1000) :: formatline

formatline = '(a,g12.5,a,i2)'
write(print_funk,fmt=formatline,iostat=error) 'v = ',funkl%v,': ndv = ',funkl%ndv
if (error /= 0) print_funk = print_funk//'ERROR: 1: problem in print_funk: formatline = '//trim(formatline)

print_funk_saved = print_funk
if (allocated(funkl%pp)) then
  formatline = '(a'//repeat(',a,i8',ubound(funkl%pp,1))//')'
  write(print_funk,fmt=formatline,iostat=error) trim(print_funk_saved)//': pp =',(' ',funkl%pp(nn),nn=1,ubound(funkl%pp,1))
  if (error /= 0) print_funk = trim(print_funk_saved)//'ERROR: 2: problem in print_funk: formatline = '//trim(formatline)
else
  formatline = '(a)'
  write(print_funk,fmt=formatline,iostat=error) trim(print_funk_saved)//': pp not allocated'
  if (error /= 0) print_funk = trim(print_funk_saved)//'ERROR: 3: problem in print_funk: formatline = '//trim(formatline)
end if

print_funk_saved = print_funk
if (allocated(funkl%dv)) then
  formatline = '(a'//repeat(',a,g12.5',ubound(funkl%dv,1))//')'
  write(print_funk,fmt=formatline,iostat=error) trim(print_funk_saved)//': dv =',(' ',funkl%dv(nn),nn=1,ubound(funkl%dv,1))
  if (error /= 0) print_funk = trim(print_funk_saved)//'ERROR: 4: problem in print_funk: formatline = '//trim(formatline)
else
  formatline = '(a)'
  write(print_funk,fmt=formatline,iostat=error) trim(print_funk_saved)//': dv not allocated'
  if (error /= 0) print_funk = trim(print_funk_saved)//'ERROR: 5: problem in print_funk: formatline = '//trim(formatline)
end if

end function print_funk

!----------------------------------------------------------------------------

function get_const(name)

character(len=*) :: name
integer :: m
double precision :: get_const

m = var_number_from_name(name)
if (m == 0.or.trim(var(m)%type) /= "constant") then
  write(*,*) 'ERROR: a constant variable with name '//trim(name)//' was not found in get_const'
  stop
end if
if (trim(var(m)%centring) /= "none") then
  write(*,*) 'ERROR: the constant variable '//trim(name)//' does not have none centring in get_const'
  stop
end if
get_const = var(m)%funk(1)%v

end function get_const

!----------------------------------------------------------------------------

function distance(x1,x2)

! this function finds distance between two points

!integer :: n
double precision, dimension(:) :: x1, x2
double precision :: distance

! write(*,*) 'in distance with: x1 = ',x1,': x2 = ',x2
! distance = 0.d0
! do n = 1, totaldimensions
!   distance = distance + (x1(n)-x2(n))**2
!   write(*,*) 'n = ',n,': distance = ',distance,' x1(n) = ',x1(n),': x2(n) = ',x2(n)
! end do
! distance = sqrt(distance)

distance = vector_magnitude(x1-x2) ! use fortran intrinsic

end function distance

!----------------------------------------------------------------------------

function cross_product(x1,x2)

! this function finds cross product between two vectors

double precision, dimension(:), intent(in) :: x1, x2
double precision, dimension(totaldimensions) :: cross_product

if (ubound(x1,1) /= totaldimensions .or. ubound(x2,1) /= totaldimensions) stop &
  'ERROR: x1 or x2 dimensions inappropriate in function cross_product'
cross_product(1) = x1(2) * x2(3) - x1(3) * x2(2)
cross_product(2) = x1(3) * x2(1) - x1(1) * x2(3)
cross_product(3) = x1(1) * x2(2) - x1(2) * x2(1)

end function cross_product

!----------------------------------------------------------------------------

function vector_magnitude(x)

! this function finds magnitude of a vector

double precision, dimension(:) :: x
double precision :: vector_magnitude

vector_magnitude = sqrt(dot_product(x,x))

end function vector_magnitude

!----------------------------------------------------------------------------

subroutine normalise_vector(x1)

! this function normalises the length of x1

double precision, dimension(:) :: x1
double precision :: length

length = vector_magnitude(x1)
if (length > 1.d-20) then
  x1 = x1/length
else
  x1 = 0.d0
  write(*,*) 'WARNING: vector of zero length found in normalise_vector'
end if
  
end subroutine normalise_vector

!-----------------------------------------------------------------

subroutine set_face(face_to_set,new_value)

! set the face in face_to_set to the new_value, including shape
type(face_type), intent(in) :: new_value
type(face_type), intent(inout) :: face_to_set
integer :: l

! make shape of allocatables in face_to_set equal to shape of new_value
! icell
call resize_integer_array(keep_data=.false.,array=face_to_set%icell,new_size=allocatable_size(new_value%icell))
! knode
call resize_integer_array(keep_data=.false.,array=face_to_set%knode,new_size=allocatable_size(new_value%knode))
! region_list
call resize_integer_array(keep_data=.false.,array=face_to_set%region_list,new_size=allocatable_size(new_value%region_list))

! copy kernels, including shape
do l = lbound(new_value%kernel,1), ubound(new_value%kernel,1)
  call copy_kernel(original=new_value%kernel(l),copy=face_to_set%kernel(l))
end do

! now set values (kernel values get set twice)
face_to_set = new_value

end subroutine set_face

!-----------------------------------------------------------------

subroutine set_cell(cell_to_set,new_value)

! set the cell in cell_to_set to the new_value, including shape
type(cell_type), intent(in) :: new_value
type(cell_type), intent(inout) :: cell_to_set
integer :: l

! make shape of allocatables in cell_to_set equal to shape of new_value
! knode
call resize_integer_array(keep_data=.false.,array=cell_to_set%knode,new_size=allocatable_size(new_value%knode))
! jface
call resize_integer_array(keep_data=.false.,array=cell_to_set%jface,new_size=allocatable_size(new_value%jface))
! icell
call resize_integer_array(keep_data=.false.,array=cell_to_set%icell,new_size=allocatable_size(new_value%icell))
! region_list
call resize_integer_array(keep_data=.false.,array=cell_to_set%region_list,new_size=allocatable_size(new_value%region_list))

! copy kernels, including shape
do l = lbound(new_value%kernel,1), ubound(new_value%kernel,1)
  call copy_kernel(original=new_value%kernel(l),copy=cell_to_set%kernel(l))
end do

! now set values (kernel values get set twice)
cell_to_set = new_value

end subroutine set_cell

!-----------------------------------------------------------------

function jface_from_knode_list(node_list)

integer :: j, jj, k, kk, jc, jjc
integer, dimension(:), allocatable :: node_list
integer :: jface_from_knode_list
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'function jface_from_knode_list'

if (debug) write(*,*) 'node_list = ',node_list

jface_from_knode_list = 0

if (.not.allocated(node(node_list(1))%jface)) return

face_loop: do jj = 1, ubound(node(node_list(1))%jface,1) ! loop through all faces that are associated with the first node
  j = node(node_list(1))%jface(jj)

  node_loop: do kk = 2, ubound(node_list,1)
    k = node_list(kk)

    if (.not.allocated(node(k)%jface)) return ! if these indices aren't even allocated, face must be new

    compare_face_loop: do jjc = 1, ubound(node(k)%jface,1)
      jc = node(k)%jface(jjc)

      if (jc == j) cycle node_loop ! face appears in both node lists so check next node

    end do compare_face_loop

    cycle face_loop ! here the original face was not found so move on to next one

  end do node_loop

  jface_from_knode_list = j ! j has been found in every node in the list - it must be the face

end do face_loop
  


! !   if (debug) write(*,*) 'kk = ',kk,': k = ',k
!     node_list_loop: do kkl = 1, ubound(node_list,1)
!       kl = node_list(kkl)
! !     if (debug) write(*,*) 'kkl = ',kkl,': kl = ',kl
!       if (k == kl) cycle node_loop ! if we find the node, move onto next node for face
! !     if (k == kl) then
! !       if (debug) write(*,*) 'matched face(j)%node (k) to node_list (kl): kk = ',kk,': k = ',k,': kkl = ',kkl,': kl = ',kl,': j = ',j
! !       cycle node_loop ! if we find the node, move onto next node for face
! !     end if
!     end do node_list_loop
! !   if (debug) write(*,*) 'did not find match for face%node node: kk = ',kk,': k = ',k,': j = ',j
!     cycle face_loop ! if we are here then the node is not present in this face, so move on to next face
!     stop 'ERROR: should not be here'
!   end do node_loop
!   jface_from_knode_list = j ! if we are here then all nodes matched - so face matches
!   exit face_loop
! end do face_loop

if (debug) then
  if (jface_from_knode_list == 0) then
    write(*,*) 'did not find face: jface_from_knode_list = ',jface_from_knode_list
  else
    write(*,*) 'found face: jface_from_knode_list = ',jface_from_knode_list,': face%knode = ',face(jface_from_knode_list)%knode
  end if
end if

if (debug) write(*,'(a/80(1h+))') 'function jface_from_knode_list'

end function jface_from_knode_list

! !-----------------------------------------------------------------

! function icell_from_jface_list(face_list)

! integer :: i, ii, k, kk, ic, iic
! integer, dimension(:), allocatable :: face_list
! integer :: icell_from_jface_list
! logical, parameter :: debug = .false.

! if (debug) write(*,'(80(1h+)/a)') 'function icell_from_jface_list'

! if (debug) write(*,*) 'face_list = ',face_list

! icell_from_jface_list = 0

! if (.not.allocated(face_list)) return

! cell_loop: do ii = 1, 2 ! loop through all cells that are associated with the first face
!   i = face(face_list(1))%icell(ii)
!   if (i == 0) cycle cell_loop

!   face_loop: do kk = 2, ubound(face_list,1)
!     k = face_list(kk)

!     compare_cell_loop: do iic = 1, 2
!       ic = face(k)%icell(iic)
!       
!       if (ic == i) cycle face_loop ! cell appears in both face lists so check next face

!     end do compare_cell_loop

!     cycle cell_loop ! here the original cell was not found so move on to next one

!   end do face_loop

!   icell_from_jface_list = i ! i has been found in every face in the list - it must be the cell

! end do cell_loop

! ! cell_loop: do i = 1, itotal
! ! 	face_loop: do jj = 1, ubound(cell(i)%jface,1)
! ! 		j = cell(i)%jface(jj)
! ! 		face_list_loop: do jjl = 1, ubound(face_list,1)
! ! 			jl = face_list(jjl)
! ! 			if (j == jl) cycle face_loop ! if we find the face, move onto next face for cell
! ! 		end do face_list_loop
! ! 		cycle cell_loop ! if we are here then the face is not present in this cell, so move on to next cell
! ! 	end do face_loop
! ! 	icell_from_jface_list = i ! if we are here then all faces matched - so cell matches
! ! 	exit cell_loop
! ! end do cell_loop

! if (debug) then
!   if (icell_from_jface_list == 0) then
!     write(*,*) 'did not find cell: icell_from_jface_list = ',icell_from_jface_list
!   else
!     write(*,*) 'found cell: icell_from_jface_list = ',icell_from_jface_list,': cell%jface = ',cell(icell_from_jface_list)%jface
!   end if
! end if

! if (debug) write(*,'(a/80(1h+))') 'function icell_from_jface_list'

! end function icell_from_jface_list

!-----------------------------------------------------------------

function icell_from_knode_list(node_list)

integer :: i, ii, k, kk, ic, iic
integer, dimension(:), allocatable :: node_list
integer :: icell_from_knode_list
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'function icell_from_knode_list'

if (debug) write(*,*) 'node_list = ',node_list

icell_from_knode_list = 0

! temp &&&
!return

if (.not.allocated(node(node_list(1))%icell)) return

cell_loop: do ii = 1, ubound(node(node_list(1))%icell,1) ! loop through all cells that are associated with the first node
  i = node(node_list(1))%icell(ii)

  node_loop: do kk = 2, ubound(node_list,1)
    k = node_list(kk)

    if (.not.allocated(node(k)%icell)) return ! if these indices aren't even allocated, cell must be new

    compare_cell_loop: do iic = 1, ubound(node(k)%icell,1)
      ic = node(k)%icell(iic)

      if (ic == i) cycle node_loop ! cell appears in both node lists so check next node

    end do compare_cell_loop

    cycle cell_loop ! here the original cell was not found so move on to next one

  end do node_loop

  icell_from_knode_list = i ! i has been found in every node in the list - it must be the cell

end do cell_loop

if (debug) then
  if (icell_from_knode_list == 0) then
    write(*,*) 'did not find cell: icell_from_knode_list = ',icell_from_knode_list
  else
    write(*,*) 'found cell: icell_from_knode_list = ',icell_from_knode_list,': cell%knode = ',cell(icell_from_knode_list)%knode
  end if
end if

if (debug) write(*,'(a/80(1h+))') 'function icell_from_knode_list'

end function icell_from_knode_list

!-----------------------------------------------------------------

subroutine find_2d_geometry(knode,area,norm,centre)

! find geometry info about a 2d convex and flat polygon
! convex means that any line drawn from one edge to another intersects only those two edges
! flat means that all points lie on the one plane (checked within this routine)
! norm is a unit normal to the surface ( normalised( (node(2)-node(1)) X (node(3)-node(2)) ) )
! area and centre are calculated by splitting surface into triangles, each with a vertex on the first node
! centre is the centroid, defined for our purposes as the point which is gives the same value
!  for a function as the average value of that function, assuming a linear dependence on space

integer, dimension(:), intent(in) :: knode ! ordered list of node indices that surround geometry
double precision, intent(out), optional :: area ! area of geometry
double precision, dimension(totaldimensions), intent(out), optional :: centre ! vectors specifying centre of surface
double precision, dimension(totaldimensions,totaldimensions), intent(out), optional :: norm ! vectors specifying normal to surface
double precision, dimension(totaldimensions) :: norms, normtriangle, centretriangle
integer :: kk
double precision :: aa, areatriangle
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h-)/a)') 'subroutine find_2d_geometry'

norms = cross_product( node(knode(2))%x-node(knode(1))%x , node(knode(3))%x-node(knode(2))%x )
call normalise_vector(norms)
if (present(norm)) then
  norm(:,1) = norms ! normal
  norm(:,2) = node(knode(2))%x-node(knode(1))%x ! first tangent, lying in plane of geometry from node 1->2
  call normalise_vector(norm(:,2))
  norm(:,3) = cross_product( norm(:,1), norm(:,2) ) ! second tangent, lying in plane of geometry
  call normalise_vector(norm(:,3))
end if

if (debug) write(*,*) 'finding area of 2d_geometry: knode = ',knode,': norms = ',norms

! new method which breaks surface into multiple triangles, each anchored to the first vertex

aa = 0.d0
if (present(centre)) centre = [ 0.d0, 0.d0, 0.d0 ]
do kk = 3, ubound(knode,1)
! so triangle has vertices knode(1), knode(kk-1) and knode(kk) in an anti-clockwise sense

! area of a triangle is half the area of a parallelgram formed between adjacent sides
  normtriangle = cross_product( node(knode(kk-1))%x-node(knode(1))%x , node(knode(kk))%x-node(knode(1))%x )
  areatriangle = vector_magnitude(normtriangle)/2.d0
  aa = aa + areatriangle

! while we're at it, check that surface is flat if it has more than 3 vertices
  if (kk > 2) then
    call normalise_vector(normtriangle)
    if (vector_magnitude(normtriangle-norms) > 1.d-10) then
      write(*,*) 'ERROR: a 2d geometry is not flat'
      write(*,*) 'normtriangle = ',normtriangle,': norms = ',norms,': knode = ',knode
      stop
    end if
  end if

! now centroid
  if (present(centre)) then
! centroid of any simplex (point, line, triangle, tetrahedron) given by sum of vertices divided by number of vertices (http://en.wikipedia.org/wiki/Centroid)
    centretriangle = (node(knode(1))%x + node(knode(kk-1))%x + node(knode(kk))%x)/3.d0
! centroid of compound form given by area weighted sum of centroid components (http://en.wikipedia.org/wiki/Centroid)
    centre = centre + areatriangle*centretriangle
  end if

end do

if (aa.lt.1.d-20) stop 'ERROR: zero or negative area in find_2d_geometry'

if (present(centre)) centre = centre/aa

if (present(area)) area = aa

! old method - does not give correct third (dummy) dimension for centroid
! area is calucated from gauss' theorem with F=x_1*delta_1
! centre is also calculated from gauss' theorem, with centre defined as point which is gives the same value
!  for a function as the average value of that function, assuming a linear dependence on space

!double precision, dimension(totaldimensions) :: norml, tangl, norms, normtriangle, centretriangle
!integer :: kk, kd, ku, kdim
!double precision :: ll, aa, normerror, areatriangle

! ! loop through faces (lines) surrounding cell in anti-clockwise direction
! aa = 0.d0
! if (present(centre)) centre = [ 0.d0, 0.d0, 0.d0 ]
! do kk = 1, ubound(knode,1)
!   kd = knode(kk)
!   if (kk == ubound(knode,1)) then
!     ku = knode(1)
!   else
!     ku = knode(kk+1)
!   end if
!   tangl = node(ku)%x-node(kd)%x ! tangent to edge
!   call normalise_vector(tangl)
!   norml = cross_product( tangl , norms ) ! outward facing normal to edge
!   call normalise_vector(norml)
!   ll = distance(node(ku)%x,node(kd)%x)
!   if (debug) write(*,*) 'incrementing area: tangl = ',tangl,': norml = ',norml,': ll = ',ll
! ! finally calculate volume using 2d gauss theorem
!   aa = aa + norml(1)*ll*( node(kd)%x(1) + tangl(1)*ll/2.d0 )
!   if (present(centre)) then
!     do kdim = 1, totaldimensions
!       centre(kdim) = centre(kdim) + (norml(kdim)/2.d0)*( node(kd)%x(kdim)**2*ll + &
!         node(kd)%x(kdim)*tangl(kdim)*ll**2 + tangl(kdim)**2*ll**3/3.d0 )
!     end do
!   end if
! end do

! if (aa.lt.1.d-20) stop 'ERROR: zero or negative area in find_2d_geometry'

! if (present(centre)) centre = centre/aa

! if (present(area)) area = aa

if (debug) write(*,'(a/80(1h-))') 'subroutine find_2d_geometry'

end subroutine find_2d_geometry

!-----------------------------------------------------------------

subroutine find_3d_geometry(jface,volume,centre)

! find geometry info about a 3d convex polyhedron
! it is assumed that each surface of the polyhedron is flat (not checked)
! convex means that any line drawn from one surface to another intersects only those two surfaces
! volume and centre are calculated by splitting surface into tetrahedra, each with a vertex on the first node of the first face
! centre is the centroid, defined for our purposes as the point which is gives the same value
!  for a function as the average value of that function, assuming a linear dependence on space

integer, dimension(:), intent(in) :: jface ! list of face indices that surround geometry
double precision, intent(out), optional :: volume ! volume of geometry
double precision, dimension(totaldimensions), intent(out), optional :: centre ! vector specifying centre of the polyhedron
double precision, dimension(totaldimensions) :: centretetrahedron
integer :: jj, kk
integer, dimension(4) :: knode
double precision :: vv, volumetetrahedron
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h-)/a)') 'subroutine find_3d_geometry'

if (debug) write(*,*) 'finding volume of 3d_geometry: jface = ',jface

! new method which breaks surface into multiple tetrahedra, each anchored to the first vertex of the first face

vv = 0.d0
if (present(centre)) centre = [ 0.d0, 0.d0, 0.d0 ]
knode(1) = face(jface(1))%knode(1) ! anchor point

jj_loop: do jj = 2, ubound(jface,1) ! cycle through faces that may possibly not have the anchor vertex as a member

! check that anchor node is not a member of this face
  do kk = 1, ubound(face(jface(jj))%knode,1)
    if (face(jface(jj))%knode(kk) == knode(1)) cycle jj_loop
  end do

! local anchor point on this face is the first node (a la find_2d_geometry)
  knode(2) = face(jface(jj))%knode(1)

! now cycle through the tetrahedron associated with this face
  do kk = 3, ubound(face(jface(jj))%knode,1)
! four nodes that define tetrahedron are now:
!   knode(1) = face(jface(1))%knode(1) - constant for all - set above
!   knode(2) = face(jface(jj))%knode(1) - constant for this face - set above
    knode(3) =  face(jface(jj))%knode(kk-1)
    knode(4) =  face(jface(jj))%knode(kk)
    
! according to http://en.wikipedia.org/wiki/Tetrahedron volume of this tet is given by
    volumetetrahedron = abs( dot_product( node(knode(2))%x-node(knode(1))%x , &
      cross_product( node(knode(3))%x-node(knode(1))%x , node(knode(4))%x-node(knode(1))%x ) ) )/6.d0
    vv = vv + volumetetrahedron

    if (present(centre)) then
! centroid of any simplex (point, line, triangle, tetrahedron) given by sum of vertices divided by number of vertices (http://en.wikipedia.org/wiki/Centroid)
      centretetrahedron = (node(knode(1))%x + node(knode(2))%x + node(knode(3))%x + node(knode(4))%x)/4.d0
! centroid of compound form given by volume weighted sum of centroid components (http://en.wikipedia.org/wiki/Centroid)
      centre = centre + volumetetrahedron*centretetrahedron
    end if

  end do

end do jj_loop

if (vv.lt.1.d-20) stop 'ERROR: zero or negative volume in find_3d_geometry'

if (present(centre)) centre = centre/vv

if (present(volume)) volume = vv

if (debug) write(*,'(a/80(1h-))') 'subroutine find_3d_geometry'

end subroutine find_3d_geometry

!----------------------------------------------------------------------------

subroutine add_jface_to_nodes(jface,knode)

! here we take a face (jface) and list of nodes (knode) and check that all the nodes have this face associated with them
integer, dimension(:), allocatable, intent(in) :: knode ! ordered list of node indices that surround geometry
integer, intent(in) :: jface
integer :: kk, k, j, jj

if (.not.allocated(knode)) return

kk_loop: do kk = 1, ubound(knode,1)
  k = knode(kk)
  if (allocated(node(k)%jface)) then
    do jj = 1, ubound(node(k)%jface,1)
      j = node(k)%jface(jj)
      if (j == jface) cycle kk_loop ! face is already in list, so move on to next node
    end do
  end if
  call push_integer_array(array=node(k)%jface,new_element=jface)
! call resize_integer_array(array=node(k)%jface,change=1)
! node(k)%jface(ubound(node(k)%jface,1)) = jface
end do kk_loop

end subroutine add_jface_to_nodes

!-----------------------------------------------------------------

subroutine add_icell_to_nodes(icell,knode)

! here we take a cell (icell) and list of nodes (knode) and check that all the nodes have this cell associated with them
integer, dimension(:), allocatable, intent(in) :: knode ! ordered list of node indices that surround geometry
integer, intent(in) :: icell
integer :: kk, k, i, ii

if (.not.allocated(knode)) return

kk_loop: do kk = 1, ubound(knode,1)
  k = knode(kk)
  if (allocated(node(k)%icell)) then
    do ii = 1, ubound(node(k)%icell,1)
      i = node(k)%icell(ii)
      if (i == icell) cycle kk_loop ! face is already in list, so move on to next node
    end do
  end if
  call push_integer_array(array=node(k)%icell,new_element=icell)
! call resize_integer_array(array=node(k)%icell,change=1)
! node(k)%icell(ubound(node(k)%icell,1)) = icell
end do kk_loop

end subroutine add_icell_to_nodes

!-----------------------------------------------------------------

function fortran_float(n)

! this is just the dble intrinsic, renamed so that it can pass through maxima

integer :: n
double precision :: fortran_float

fortran_float = dble(n)

end function fortran_float

!-----------------------------------------------------------------

function heaviside(x)

! this is the heaviside function

double precision :: x, heaviside

if (x > 0) then
  heaviside = 1.d0
else if (x < 0) then
  heaviside = 0.d0
else
  heaviside = 0.5d0
end if

end function heaviside

!-----------------------------------------------------------------

function divop(i,j)

! assemble coefficients required for divergence sum in domain equations

integer, intent(in) :: i,j
double precision :: divop
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'function divop'

divop = face(j)%area/cell(i)%vol
if (face(j)%icell(1) /= i) divop = -divop ! reverse sign if dot product of nface and ncell isn't 1

if (debug) write(*,'(a/80(1h-))') 'function divop'

end function divop

!-----------------------------------------------------------------

subroutine time_process(description)

character(len=*), intent(in), optional :: description
real :: this_cpu_time
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'subroutine time_process'

if (.not.output_timings) then
  if (debug) write(*,*) 'logical output_timings set to .false. - not doing timings'
  return
end if

call cpu_time(this_cpu_time)

if (present(description)) write(*,'(a,a,a,g11.4,a)') 'INFO: cpu time to complete ',trim(description), &
    ' routines = ',this_cpu_time-last_cpu_time,' s'

last_cpu_time = this_cpu_time

if (debug) write(*,'(a/80(1h-))') 'subroutine time_process'

end subroutine time_process

!-----------------------------------------------------------------

subroutine copy_kernel(original,copy)

type(kernel_type) :: original, copy
integer :: new_size
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'subroutine copy_kernel'

copy%centring = original%centring

new_size = 0
if (allocated(original%ijk)) new_size = ubound(original%ijk,1)

if (allocated(copy%ijk)) deallocate(copy%ijk)
if (allocated(copy%v)) deallocate(copy%v)

if (new_size == 0) return

allocate(copy%ijk(new_size),copy%v(new_size))
!copy = original ! this should work...
copy%ijk = original%ijk
copy%v = original%v

if (debug) write(*,'(a/80(1h-))') 'subroutine copy_kernel'

end subroutine copy_kernel

!-----------------------------------------------------------------

function location_in_list(array,element)

integer, dimension(:), allocatable :: array
integer :: element
integer :: location_in_list
integer :: i
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'function location_in_list'

location_in_list = 0

if (.not.allocated(array)) return

do i = 1, ubound(array,1)
  if (element == array(i)) then
    location_in_list = i
    exit
  end if
end do

if (debug) write(*,'(a/80(1h-))') 'function location_in_list'

end function location_in_list

!-----------------------------------------------------------------

function location_in_list_dummy(array,element)

! same as location in list but instead accepts a non-allocatable array

integer, dimension(:) :: array
integer :: element
integer :: location_in_list_dummy
integer :: i
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'function location_in_list_dummy'

location_in_list_dummy = 0

do i = 1, ubound(array,1)
  if (element == array(i)) then
    location_in_list_dummy = i
    exit
  end if
end do

if (debug) write(*,'(a/80(1h-))') 'function location_in_list_dummy'

end function location_in_list_dummy

!-----------------------------------------------------------------

function region_delta(ij,centring,name)

character(len=4) :: centring ! whether region is cell or face centred
character(len=*) :: name ! name of the region
double precision :: region_delta
integer :: ij, region_number
logical :: existing
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'function region_delta'

region_delta = 0.d0

region_number = region_number_from_name(name=name,existing=existing,centring=centring)
if (.not.existing) then
  write(*,*) 'ERROR: function region_delta references region '//trim(name)// &
    ' which does not exist'
  stop
end if
if (region_number == 0) then
  write(*,*) 'ERROR: problem in function region_delta with region '//trim(name)// &
    ':- most likely region has wrong centring'
  stop
end if

if (centring == 'cell') then
  if (location_in_list(array=cell(ij)%region_list,element=region_number) == 0) return
else if (centring == 'face') then
  if (location_in_list(array=face(ij)%region_list,element=region_number) == 0) return
else
  stop "ERROR: incorrect centring in region_delta"
end if

region_delta = 1.d0

if (debug) write(*,'(a/80(1h-))') 'function region_delta'

end function region_delta

!-----------------------------------------------------------------

function var_list_number(type,centring)

! little function to return number of var_list that corresponds to type and centring
character(len=*) :: centring ! whether region is cell or face centred
character(len=*) :: type ! type of var variable
integer :: var_list_number
integer :: ntype, ncentring, n
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'function var_list_number'

if (trim(type) == "all") then
  ntype = ubound(var_types,1) + 1
else
  ntype = 1
  do while (ntype <= ubound(var_types,1))
    if (trim(type) == trim(var_types(ntype))) exit
    ntype = ntype + 1
  end do
  if (ntype > ubound(var_types,1)) call error_stop('unknown type '//trim(type)//'in var_list_number')
end if

!if (trim(type) == "constant") then
!  ntype = 1
!else if (trim(type) == "transient") then
!  ntype = 2
!else if (trim(type) == "unknown") then
!  ntype = 3
!else if (trim(type) == "derived") then
!  ntype = 4
!else if (trim(type) == "equation") then
!  ntype = 5
!else if (trim(type) == "output") then
!  ntype = 6
!else if (trim(type) == "condition") then
!  ntype = 7
!else if (trim(type) == "local") then
!  ntype = 8
!else if (trim(type) == "all") then
!  ntype = 9
!else
!  stop "ERROR: unknown type in var_list_number"
!end if

if (trim(centring) == "cell") then
  ncentring = 1
else if (trim(centring) == "face") then
  ncentring = 2
else if (trim(centring) == "none") then
  ncentring = 3
else if (trim(centring) == "all") then
  ncentring = 4
else
  stop "ERROR: unknown centring in var_list_number"
end if

var_list_number = ntype + (ncentring-1)*(ubound(var_types,1)+1)

if (debug) write(*,'(a/80(1h-))') 'function var_list_number'

end function var_list_number

!-----------------------------------------------------------------

function var_list_lookup(type,centring)

! little function to return var_list as an array corresponding to type and centring
character(len=*) :: centring ! centring of var
character(len=*) :: type ! type of var variable
integer, dimension(:), allocatable :: var_list_lookup
integer :: var_list_length
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'function var_list_lookup'

var_list_length = allocatable_size(var_list(var_list_number(type=type,centring=centring))%list)

if (allocated(var_list_lookup)) deallocate(var_list_lookup)

allocate(var_list_lookup(var_list_length)) ! even if list is empty allocate with zero size, allowable in modern fortran

var_list_lookup = var_list(var_list_number(type=type,centring=centring))%list

if (debug) write(*,'(a/80(1h-))') 'function var_list_lookup'

end function var_list_lookup

!-----------------------------------------------------------------

function nsvar(m,ij,noerror)

! little function to lookup ns data number corresponding to location ij for var(m)
integer, intent(in) :: m
integer, intent(in), optional :: ij
logical, intent(in), optional :: noerror ! do not report an error even if out of range detected
logical :: noerror_l
integer :: nsvar

if (m > ubound(var,1) .or. m < 1) stop 'ERROR: variable index m out of range in nsvar'

if (var(m)%centring == "none") then
  nsvar = 1
else
  if (.not.present(ij)) call error_stop('ij is not present for a cell or face centred var in function nsvar: This '// &
    'means that an attempt is being made to reference variable '//trim(var(m)%name)//' at an erroreous location. '// &
    'Look at the use of this variable in the equations and check that its region and centring context is correct.')
  if (ij == 0) call error_stop('ij is equal to 0 in function nsvar: This means that an attempt is being made '// &
    'to reference variable '//trim(var(m)%name)//' at an erroreous location. '// &
    'Look at the use of this variable in the equations and check that its region and centring context is correct.')
  nsvar = region(var(m)%region_number)%ns(ij)
  if (nsvar == 0) then
    noerror_l = .false.
    if (present(noerror)) noerror_l = noerror
    if (noerror_l) then
      nsvar = 0
    else
      call error_stop('nsvar is equal to 0 in nsvar: This means that an attempt is being made '// &
        'to reference '//trim(var(m)%centring)//' centred variable '//trim(var(m)%name)//' outside of the region '// &
        trim(var(m)%region)//' in which it is defined. '// &
        'Look at the use of this variable in the equations and check that its region and centring context is correct.')
    end if
  end if
end if

end function nsvar

!-----------------------------------------------------------------

function ijvar(m,ns)

! little function to lookup ij index corresponding to data number ns for var(m)
integer :: m
integer :: ns
integer :: ijvar

if (m > ubound(var,1) .or. m < 1) stop 'ERROR: m out of range in function ijvar'

if (var(m)%centring == "none") then
  ijvar = 1
else
  ijvar = region(var(m)%region_number)%ij(ns)
end if

end function ijvar

!-----------------------------------------------------------------

subroutine reset_face(default_face)

! resets a face element
type(face_type) :: default_face
integer :: n
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'subroutine reset_face'

default_face%type = 0
default_face%x = 0.d0
default_face%area = 0.d0
default_face%dx = 0.d0
default_face%norm = 0.d0
if (allocated(default_face%icell)) deallocate(default_face%icell)
allocate(default_face%icell(2))
if (allocated(default_face%knode)) deallocate(default_face%knode)
if (allocated(default_face%region_list)) deallocate(default_face%region_list)
default_face%dimensions = 0
default_face%gtype = 0
!do n = 0, 2*totaldimensions
do n = lbound(default_face%kernel,1), ubound(default_face%kernel,1)
  default_face%kernel(n)%centring = ''
  if (allocated(default_face%kernel(n)%ijk)) deallocate(default_face%kernel(n)%ijk)
  if (allocated(default_face%kernel(n)%v)) deallocate(default_face%kernel(n)%v)
end do

if (debug) write(*,'(a/80(1h-))') 'subroutine reset_face'

end subroutine reset_face

!-----------------------------------------------------------------

subroutine reset_cell(default_cell)

! resets a cell element
type(cell_type) :: default_cell
integer :: n
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'subroutine reset_cell'

default_cell%type = 0
default_cell%x = 0.d0
default_cell%vol = 0.d0
if (allocated(default_cell%knode)) deallocate(default_cell%knode)
if (allocated(default_cell%jface)) deallocate(default_cell%jface)
if (allocated(default_cell%icell)) deallocate(default_cell%icell)
if (allocated(default_cell%region_list)) deallocate(default_cell%region_list)
default_cell%dimensions = 0
default_cell%gtype = 0
!do n = 0, totaldimensions
do n = lbound(default_cell%kernel,1), ubound(default_cell%kernel,1)
  default_cell%kernel(n)%centring = ''
  if (allocated(default_cell%kernel(n)%ijk)) deallocate(default_cell%kernel(n)%ijk)
  if (allocated(default_cell%kernel(n)%v)) deallocate(default_cell%kernel(n)%v)
end do

if (debug) write(*,'(a/80(1h-))') 'subroutine reset_cell'

end subroutine reset_cell

!----------------------------------------------------------------------------

subroutine copy_integer_array(original,copy)

! copy allocatable original array to copy array

integer, dimension(:), allocatable :: original, copy

call resize_integer_array(keep_data=.false.,array=copy,new_size=allocatable_integer_size(original))
if (allocated(original)) copy = original

end subroutine copy_integer_array

!-----------------------------------------------------------------

subroutine copy_double_precision_array(original,copy)

! copy allocatable original array to copy array

double precision, dimension(:), allocatable :: original, copy

call resize_double_precision_array(keep_data=.false.,array=copy,new_size=allocatable_double_precision_size(original))
if (allocated(original)) copy = original

end subroutine copy_double_precision_array

!-----------------------------------------------------------------

subroutine copy_character_array(original,copy)

! copy allocatable original array to copy array

character(len=*), dimension(:), allocatable :: original
character(len=len(original)), dimension(:), allocatable :: copy

call resize_character_array(keep_data=.false.,array=copy,new_size=allocatable_character_size(original))
if (allocated(original)) copy = original

end subroutine copy_character_array

!-----------------------------------------------------------------

subroutine copy_real_array(original,copy)

! copy allocatable original array to copy array

real, dimension(:), allocatable :: original, copy

call resize_real_array(keep_data=.false.,array=copy,new_size=allocatable_real_size(original))
if (allocated(original)) copy = original

end subroutine copy_real_array

!-----------------------------------------------------------------

function elements_in_common(array1,array2)

! this returns the number of common elements between 2 integer arrays
! could do fancier search by reordering/reducing array2 but would
! require the array to be copied

integer, dimension(:), allocatable :: array1, array2
integer :: elements_in_common
integer :: n, n2

elements_in_common = 0

if (.not.allocated(array1).or..not.allocated(array2)) return

do n = 1, ubound(array1,1)
  n2 = location_in_list(array=array2,element=array1(n))
  if (n2 /= 0) elements_in_common = elements_in_common + 1
end do

end function elements_in_common

!-----------------------------------------------------------------

subroutine memory_manage_dvs(type,action)

! here we deallocate or reallocate the dv arrays in the var%funks

character(len=*) :: type ! variable type this is applied to
character(len=*) :: action ! deallocate or reallocate
integer :: n, m, ns
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'subroutine memory_manage_dvs'

if (debug) write(*,*) 'doing action = '//trim(action)//' for type '//trim(type)

if (debug) pause

if (trim(action) == 'deallocate') then
  do n = 1, allocatable_size(var_list(var_list_number(centring="all",type=type))%list)
    m = var_list(var_list_number(centring="all",type=type))%list(n)
    if (debug) write(*,*) 'deallocating dv for variable '//trim(var(m)%name)
    do ns = 1, ubound(var(m)%funk,1)
      var(m)%funk(ns)%ndv = -abs(var(m)%funk(ns)%ndv) ! set ndv to be negative of previous value as a marker
      if (allocated(var(m)%funk(ns)%pp)) deallocate(var(m)%funk(ns)%pp)
      if (allocated(var(m)%funk(ns)%dv)) deallocate(var(m)%funk(ns)%dv)
    end do
  end do
else if (trim(action) == 'reallocate') then
  do n = 1, allocatable_size(var_list(var_list_number(centring="all",type=type))%list)
    m = var_list(var_list_number(centring="all",type=type))%list(n)
    if (debug) write(*,*) 'reallocating dv for variable '//trim(var(m)%name)
    do ns = 1, ubound(var(m)%funk,1)
      if (var(m)%funk(ns)%ndv < 0) then
        if (allocated(var(m)%funk(ns)%pp)) deallocate(var(m)%funk(ns)%pp) ! should not need this
        if (allocated(var(m)%funk(ns)%dv)) deallocate(var(m)%funk(ns)%dv) ! should not need this
        allocate(var(m)%funk(ns)%pp(-var(m)%funk(ns)%ndv))
        allocate(var(m)%funk(ns)%dv(-var(m)%funk(ns)%ndv))
        var(m)%funk(ns)%ndv = 0
      end if
    end do
  end do
else
  stop 'ERROR: unknown action in memory_manage_dvs'
end if

if (debug) pause

if (debug) write(*,'(a/80(1h-))') 'subroutine memory_manage_dvs'

end subroutine memory_manage_dvs

!-----------------------------------------------------------------

function scanstring(string,substring)

! this routine searches string for substring, and if found, returns
!  the index of the starting character

character(len=*) :: string
character(len=*) :: substring
integer :: scanstring
integer :: l, lstring, lsubstring

scanstring = 0
lstring = len(string)
lsubstring = len(substring)

do l = 1, lstring-lsubstring+1
  if (string(l:l+lsubstring-1) == substring) then
    scanstring = l
    return
  end if
end do

end function scanstring

!-----------------------------------------------------------------

function unknown_var_from_pp(pp)

! search for unknown variable that pp refers to

integer :: unknown_var_from_pp, nnvar, mm, pp

unknown_var_from_pp = 0
mm = 0

do nnvar = 1, allocatable_size(var_list(var_list_number(centring="all",type="unknown"))%list)
  mm = var_list(var_list_number(centring="all",type="unknown"))%list(nnvar)
  if (var(mm)%funk(1)%pp(1) <= pp .and. var(mm)%funk(ubound(var(mm)%funk,1))%pp(1) >= pp) then
    unknown_var_from_pp = mm
    exit
  end if
end do
      
end function unknown_var_from_pp

!-----------------------------------------------------------------

! fortran 2003 requires gfortran 4.6
!function number_is_valid(scalar)
!
!! check that a number is neither a NaN or infinite
!
!use ieee_arithmetic
!
!double precision, intent(in) :: scalar
!logical :: number_is_valid
!
!number_is_valid = .false.
!if (ieee_is_nan(scalar)) return
!if (.not.ieee_is_finite(scalar)) return
!number_is_valid = .true.
!      
!end function number_is_valid

function number_is_valid(scalar)

! check that a number is neither a NaN or infinite

double precision, intent(in) :: scalar
logical :: number_is_valid

number_is_valid = .false.
if (scalar /= scalar) return
if (abs(scalar) > huge(scalar)) return
number_is_valid = .true.
      
end function number_is_valid

!-----------------------------------------------------------------

function variable_location_string(m,ns)

! little routine to create a formated string describing a variable location

integer :: m, ns, i, j
character(len=1000) :: variable_location_string
character(len=1000) :: formatline

if (var(m)%centring == 'cell') then
  i = region(var(m)%region_number)%ij(ns)
  formatline = '(a,'//trim(dindexformat(i))//',a,3(1x,g10.4),a)'
  write(variable_location_string,fmt=formatline) 'cell i = ',i,' and x =',cell(i)%x, &
    ' ('//trim(var(m)%type)//' '//trim(var(m)%name)//')'
else if (var(m)%centring == 'face') then
  j = region(var(m)%region_number)%ij(ns)
  formatline = '(a,'//trim(dindexformat(j))//',a,3(1x,g10.4),a)'
  write(variable_location_string,fmt=formatline) 'face j = ',j,' and x =',face(j)%x, &
    ' ('//trim(var(m)%type)//' '//trim(var(m)%name)//')'
else
  write(variable_location_string,'(a)') 'nowhere ('//trim(var(m)%type)//' '//trim(var(m)%name)//' has none centring)'
end if

end function variable_location_string

!----------------------------------------------------------------------------

function dindexformat(i)

! this function gives the minimum integer format for a given number
! d is for dynamic

character(len=3) :: dindexformat
integer :: i, ipositive, ilength, icompare

if (i < 0) then
  ilength = 2
  ipositive = -i
else
  ilength = 1
  ipositive = i
end if

icompare = 10
do while (ipositive >= icompare)
  ilength = ilength + 1
  icompare = icompare*10
end do

write(dindexformat,'(i2)') ilength
dindexformat = 'i'//adjustl(dindexformat)

end function dindexformat

!-----------------------------------------------------------------

function basename(fullname)

! this function extracts a file's basename from a file's fullname

character(len=1000) :: basename
character(len=*) :: fullname
integer :: cutr, cutl

cutl=scan(fullname,'/',.true.) ! find rightmost occurance of directory separator
cutr=scan(fullname,'.',.true.) ! find rightmost occurance of . indicating extension
if (cutr <= cutl) cutr = len(fullname)+1
if (cutr - 1 - (cutl + 1) > 1000) stop "ERROR: file name too long in basename"
basename = trim(adjustl(fullname(cutl+1:cutr-1)))

end function basename

!-----------------------------------------------------------------

subroutine error_stop(comment)

! little subroutine to stop simulation and printout comment

character(len=*), optional :: comment
logical :: fconverge_opened, foutputstep_opened

inquire(unit=fconverge,opened=fconverge_opened)
inquire(unit=foutputstep,opened=foutputstep_opened)

if (present(comment)) then
  write(*,*) 'ERROR: '//trim(comment)
  if (fconverge_opened) write(fconverge,*) 'ERROR: '//trim(comment)
end if

if (fconverge_opened) close(fconverge)
if (foutputstep_opened) close(foutputstep)

stop

end subroutine error_stop

!-----------------------------------------------------------------

function check_option(array,possibilities)

character(len=*), dimension(:), allocatable, intent(in) :: array
character(len=*), dimension(:), intent(in) :: possibilities
character(len=100) :: check_option
integer :: n,m

check_option = ""

if (.not.allocated(array)) return

do n = 1, ubound(array,1) ! forward lookup
!do n = ubound(array,1), 1, -1 ! reverse lookup
  do m = 1, ubound(possibilities,1)
    if (trim(array(n)) == trim(possibilities(m))) then
      check_option = trim(array(n))
      return
    end if
  end do
end do

end function check_option

!-----------------------------------------------------------------

! ifort does not like passing nondimensioned nonsized character arrays as arguments, gives weird memory faults
!function check_stopfile(file_array)

! little file to check on the existence of some files

!character(len=*), dimension(:), intent(in) :: file_array
!logical :: check_stopfile
!integer :: n
!
!check_stopfile = .false.
!do n = 1, ubound(file_array,1)
!  inquire(file=trim(file_array(n)),exist=check_stopfile)
!  if (check_stopfile) return
!end do
    
!end function check_stopfile

function check_stopfile(file)

! little function to check on the existence of a stop or action file

character(len=*), intent(in) :: file
logical :: check_stopfile

inquire(file=trim(file),exist=check_stopfile)

end function check_stopfile

!-----------------------------------------------------------------

subroutine ring_bell

! little subroutine to ring the bell
! NB: this string is entered via cntrl-v cntrl-g in vi
write(*,*) ""

end subroutine ring_bell

!-----------------------------------------------------------------

end module general_module

!----------------------------------------------------------------------------
