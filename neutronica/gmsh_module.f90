! file src/gmsh_module.f90
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
module gmsh_module

! statements specifying different data types and parameters
implicit none

! these are the only subroutines that are accessible from outside this module
private
public setup_gmesh, finalise_gmesh, push_gmesh

type node_list_type
  integer, dimension(:), allocatable :: node
end type node_list_type
  
! this data structure holds info about how elements, faces and nodes are defined in gmsh
type gtype_list_type
  character(len=100) :: name
  integer :: dimensions ! 0 for point, 1 for line, 2 for plane, 3 for volume
  integer :: nnodes ! number of surrounding nodes
  integer :: nfaces ! number of surrounding faces
  type(node_list_type), dimension(:), allocatable :: face_nodes ! lists the node indicies of each of the nfaces faces
  logical :: supported ! suitable for arb?
end type gtype_list_type

! reference list of all the gtypes
type(gtype_list_type), public, dimension(:), allocatable :: gtype_list

! info about each of the gelements
type gelement_type
  integer :: icell ! corresponding icell, 0 if none
  integer :: jface ! corresponding jface, 0 if none
  integer, dimension(:), allocatable :: gregions ! list of gregions that it is a member of, or failing that, one entry of 0
end type gelement_type

! info about each of the gregions
type gregion_type
  integer :: region_number
  character(len=4) :: centring
end type gregion_type

! type for each gmsh that is read in or output
type gmesh_type
  character(len=1000) :: filename
  character(len=1000) :: basename
  character(len=100), dimension(:), allocatable :: options ! array of options for this gmesh
  integer :: dimensions
  integer :: ngnodes ! number of gnodes
  integer, dimension(:), allocatable :: knode_from_gnode ! 0 if undefined
  integer, dimension(:), allocatable :: gnode_from_knode ! 0 if undefined
  integer :: ngelements ! number of gelements, including repeats for multiple regions
  integer :: ngelements_face ! number of face gelements, including repeats for multiple regions
  integer :: ngelements_cell ! number of cell gelements, including repeats for multiple regions
  type(gelement_type), dimension(:), allocatable :: gelement ! list of gelements
  integer :: ngregions ! number of gregions
  type(gregion_type), dimension(:), allocatable :: gregion ! mapping from gregion number to arb region number, dimension gregion_max
end type gmesh_type

! gmsh meshes
type(gmesh_type), public, dimension(:), allocatable :: gmesh

!-----------------------------------------------------------------
contains

!-----------------------------------------------------------------

subroutine setup_gmesh

! here we allocate data to the reference gtype_list array

integer :: n
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'subroutine setup_gmesh'

allocate(gtype_list(31))

! default is that element is unsupported
do n=1,ubound(gtype_list,1)
  gtype_list(n)%name = "unknown"
  gtype_list(n)%supported = .false.
end do

! now set supported elements
! faces are definined in the context of the element having the same dimensions as the problem
! ie, as a cell
! if the element has a dimension that is less than the problem dimension, it is assumed to be a
! face (dimension-1)
! surface nodes are directed using the right-hand rule with a positive outward pointing normal

n = 15
gtype_list(n)%name = "1-node point"
gtype_list(n)%dimensions = 0
gtype_list(n)%nnodes = 1
gtype_list(n)%nfaces = 0
gtype_list(n)%supported = .true.

n = 1
gtype_list(n)%name = "2-node line"
gtype_list(n)%dimensions = 1
gtype_list(n)%nnodes = 2
gtype_list(n)%nfaces = 2
allocate(gtype_list(n)%face_nodes(gtype_list(n)%nfaces))
allocate(gtype_list(n)%face_nodes(1)%node(1))
gtype_list(n)%face_nodes(1)%node = [ 1 ]
allocate(gtype_list(n)%face_nodes(2)%node(1))
gtype_list(n)%face_nodes(2)%node = [ 2 ]
gtype_list(n)%supported = .true.

n = 2
gtype_list(n)%name = "3-node triangle"
gtype_list(n)%dimensions = 2
gtype_list(n)%nnodes = 3
gtype_list(n)%nfaces = 3
allocate(gtype_list(n)%face_nodes(gtype_list(n)%nfaces))
allocate(gtype_list(n)%face_nodes(1)%node(2))
gtype_list(n)%face_nodes(1)%node = [ 1, 2 ]
allocate(gtype_list(n)%face_nodes(2)%node(2))
gtype_list(n)%face_nodes(2)%node = [ 2, 3 ]
allocate(gtype_list(n)%face_nodes(3)%node(2))
gtype_list(n)%face_nodes(3)%node = [ 3, 1 ]
gtype_list(n)%supported = .true.

n = 3
gtype_list(n)%name = "4-node quadrangle"
gtype_list(n)%dimensions = 2
gtype_list(n)%nnodes = 4
gtype_list(n)%nfaces = 4
allocate(gtype_list(n)%face_nodes(gtype_list(n)%nfaces))
allocate(gtype_list(n)%face_nodes(1)%node(2))
gtype_list(n)%face_nodes(1)%node = [ 1, 2 ]
allocate(gtype_list(n)%face_nodes(2)%node(2))
gtype_list(n)%face_nodes(2)%node = [ 2, 3 ]
allocate(gtype_list(n)%face_nodes(3)%node(2))
gtype_list(n)%face_nodes(3)%node = [ 3, 4 ]
allocate(gtype_list(n)%face_nodes(4)%node(2))
gtype_list(n)%face_nodes(4)%node = [ 4, 1 ]
gtype_list(n)%supported = .true.

n = 4
gtype_list(n)%name = "4-node tetrahedron"
gtype_list(n)%dimensions = 3
gtype_list(n)%nnodes = 4
gtype_list(n)%nfaces = 4
allocate(gtype_list(n)%face_nodes(gtype_list(n)%nfaces))
allocate(gtype_list(n)%face_nodes(1)%node(3))
gtype_list(n)%face_nodes(1)%node = [ 1, 4, 3 ]
allocate(gtype_list(n)%face_nodes(2)%node(3))
gtype_list(n)%face_nodes(2)%node = [ 2, 3, 4 ]
allocate(gtype_list(n)%face_nodes(3)%node(3))
gtype_list(n)%face_nodes(3)%node = [ 1, 2, 4 ]
allocate(gtype_list(n)%face_nodes(4)%node(3))
gtype_list(n)%face_nodes(4)%node = [ 1, 3, 2 ]
gtype_list(n)%supported = .true.

n = 5
gtype_list(n)%name = "8-node hexahedron"
gtype_list(n)%dimensions = 3
gtype_list(n)%nnodes = 8
gtype_list(n)%nfaces = 6
allocate(gtype_list(n)%face_nodes(gtype_list(n)%nfaces))
allocate(gtype_list(n)%face_nodes(1)%node(4))
gtype_list(n)%face_nodes(1)%node = [ 1, 5, 8, 4 ]
allocate(gtype_list(n)%face_nodes(2)%node(4))
gtype_list(n)%face_nodes(2)%node = [ 5, 6, 7, 8 ]
allocate(gtype_list(n)%face_nodes(3)%node(4))
gtype_list(n)%face_nodes(3)%node = [ 2, 3, 7, 6 ]
allocate(gtype_list(n)%face_nodes(4)%node(4))
gtype_list(n)%face_nodes(4)%node = [ 1, 4, 3, 2 ]
allocate(gtype_list(n)%face_nodes(5)%node(4))
gtype_list(n)%face_nodes(5)%node = [ 1, 2, 6, 5 ]
allocate(gtype_list(n)%face_nodes(6)%node(4))
gtype_list(n)%face_nodes(6)%node = [ 3, 4, 8, 7 ]
gtype_list(n)%supported = .true.

n = 6
gtype_list(n)%name = "6-node prism"
gtype_list(n)%dimensions = 3
gtype_list(n)%nnodes = 6
gtype_list(n)%nfaces = 5
allocate(gtype_list(n)%face_nodes(gtype_list(n)%nfaces))
allocate(gtype_list(n)%face_nodes(1)%node(4))
gtype_list(n)%face_nodes(1)%node = [ 1, 2, 5, 4 ]
allocate(gtype_list(n)%face_nodes(2)%node(4))
gtype_list(n)%face_nodes(2)%node = [ 2, 3, 6, 5 ]
allocate(gtype_list(n)%face_nodes(3)%node(4))
gtype_list(n)%face_nodes(3)%node = [ 1, 4, 6, 3 ]
allocate(gtype_list(n)%face_nodes(4)%node(3))
gtype_list(n)%face_nodes(4)%node = [ 1, 3, 2 ]
allocate(gtype_list(n)%face_nodes(5)%node(3))
gtype_list(n)%face_nodes(5)%node = [ 4, 5, 6 ]
gtype_list(n)%supported = .true.

n = 7
gtype_list(n)%name = "5-node pyramid"
gtype_list(n)%dimensions = 3
gtype_list(n)%nnodes = 5
gtype_list(n)%nfaces = 5
allocate(gtype_list(n)%face_nodes(gtype_list(n)%nfaces))
allocate(gtype_list(n)%face_nodes(1)%node(3))
gtype_list(n)%face_nodes(1)%node = [ 1, 2, 5 ]
allocate(gtype_list(n)%face_nodes(2)%node(3))
gtype_list(n)%face_nodes(2)%node = [ 2, 3, 5 ]
allocate(gtype_list(n)%face_nodes(3)%node(3))
gtype_list(n)%face_nodes(3)%node = [ 3, 4, 5 ]
allocate(gtype_list(n)%face_nodes(4)%node(3))
gtype_list(n)%face_nodes(4)%node = [ 1, 5, 4 ]
allocate(gtype_list(n)%face_nodes(5)%node(4))
gtype_list(n)%face_nodes(5)%node = [ 1, 4, 3, 2 ]
gtype_list(n)%supported = .true.

if (debug) write(*,'(a/80(1h+))') 'subroutine setup_gmesh'

end subroutine setup_gmesh

!-----------------------------------------------------------------

subroutine finalise_gmesh

! here we calculate some reverse lookups and count elements

use general_module
integer :: gmesh_number, region_number, i, j, k, n, gregion, gnode, knode, gelement, iadd, jadd, &
  new_size, old_size, jj, o, nregion
integer, dimension(:), allocatable :: face_in_gregion
type(gelement_type), dimension(:), allocatable :: gelement_old
character(len=1000) :: formatline
logical, parameter :: debug = .false.
logical :: debug_sparse = .true. 
                  
if (debug) debug_sparse = .true.

if (debug_sparse) write(*,'(80(1h+)/a)') 'subroutine finalise_gmesh'

do gmesh_number = 0, ubound(gmesh,1)

  if (debug) write(*,*) 'gmesh_number = ',gmesh_number

!------------------------------
! setup output gmesh
  if (gmesh_number == 0) then

! setup output gmesh indices - completely overwrite any that were there previously
! knode and gnode have the same index
    if (allocated(gmesh(0)%knode_from_gnode)) deallocate(gmesh(0)%knode_from_gnode)
    if (ktotal > 0) then
      allocate(gmesh(0)%knode_from_gnode(ktotal))
      do k = 1, ktotal
        gmesh(0)%knode_from_gnode(k) = k
      end do
    end if
! gelements are arranged as: domain cells, boundary cells/faces, domain faces
    if (allocated(gmesh(0)%gelement)) deallocate(gmesh(0)%gelement)
    if (itotal+jdomain > 0) then
      allocate(gmesh(0)%gelement(itotal+jdomain))
      gmesh(0)%gelement(:)%icell = 0
      gmesh(0)%gelement(:)%jface = 0
      do i = 1, idomain
        gmesh(0)%gelement(i)%icell = i
      end do
      do i = idomain+1, itotal ! the same gmsh boundary element can refer to both the boundary cell and face
        gmesh(0)%gelement(i)%icell = i
        gmesh(0)%gelement(i)%jface = cell(i)%jface(1)
      end do
      n = itotal
      do j = 1, jtotal
        if (face(j)%type == 1) then
          n = n + 1
          if (n > itotal+jdomain) stop "ERROR: problem with setting up output gmesh in setup_mesh"
          gmesh(0)%gelement(n)%jface = j
        end if
      end do
    end if

! set output gmesh dimensions
    if (allocated(cell)) then
      gmesh(0)%dimensions = maxval(cell(:)%dimensions)
    else
      gmesh(0)%dimensions = 0
    end if

! include all non-SYSTEM regions in output gmesh
    if (allocated(gmesh(0)%gregion)) deallocate(gmesh(0)%gregion)
    nregion  = 0
    if (allocated(region)) then
      do region_number = 1, ubound(region,1)
        if (region(region_number)%location(1:6) /= 'SYSTEM') nregion = nregion + 1
      end do
      if (nregion > 0) then
        allocate(gmesh(0)%gregion(nregion))
        nregion  = 0
        do region_number = 1, ubound(region,1)
          if (region(region_number)%location(1:6) /= 'SYSTEM') then
            nregion = nregion + 1
            gmesh(0)%gregion(nregion)%region_number = region_number
            gmesh(0)%gregion(nregion)%centring = region(region_number)%centring
          end if
        end do
      end if
    end if

  else
!------------------------------
! operations for gmsh input gmeshes

! also include any faces that were not included originally but that border cells that are within this gmesh
    if (jtotal > 0) then
      allocate(face_in_gregion(jtotal))
      face_in_gregion = 0 ! indicates face is not in region
! first mark any faces that are already given gelements as 2
      if (allocated(gmesh(gmesh_number)%gelement)) then
        old_size = ubound(gmesh(gmesh_number)%gelement,1)
        do gelement = 1, old_size
          if (gmesh(gmesh_number)%gelement(gelement)%jface /= 0) face_in_gregion(gmesh(gmesh_number)%gelement(gelement)%jface) = 2
        end do
      else
        old_size = 0
      end if
! now loop through all cells in region, identifying faces that border these cells
      new_size = old_size
      do gelement = 1, old_size
        i = gmesh(gmesh_number)%gelement(gelement)%icell
        if (i /= 0) then
          do jj = 1, ubound(cell(i)%jface,1)
            j = cell(i)%jface(jj)
            if (face_in_gregion(j) == 0) then
              face_in_gregion(j) = 1 ! face is in region but has not been previously defined, so define it
              new_size = new_size + 1
            end if
          end do
        end if
      end do
! now expand the gmesh%gelement array to include the new faces
      if (new_size /= old_size) then
        if (old_size > 0) then
          allocate(gelement_old(old_size))
          gelement_old = gmesh(gmesh_number)%gelement
        end if
        deallocate(gmesh(gmesh_number)%gelement)
        allocate(gmesh(gmesh_number)%gelement(new_size))
        gmesh(gmesh_number)%gelement(:)%icell = 0
        gmesh(gmesh_number)%gelement(:)%jface = 0
        if (old_size > 0) then
          gmesh(gmesh_number)%gelement(1:old_size) = gelement_old ! allocate all the old data across
          deallocate(gelement_old)
        end if
        gelement = old_size
        do j = 1, jtotal
          if (face_in_gregion(j) == 1) then
            gelement = gelement + 1
            gmesh(gmesh_number)%gelement(gelement)%jface = j
          end if
        end do
        if (gelement /= new_size) stop "ERROR: problem in finalise_gmesh"
      end if
      deallocate(face_in_gregion)
    end if

  end if
!------------------------------
! operations for all gmeshes

! reverse lookup and counting of gnodes
  gmesh(gmesh_number)%ngnodes = 0
  if (ktotal > 0) then
    call resize_integer_array(array=gmesh(gmesh_number)%gnode_from_knode,new_size=ktotal)
    gmesh(gmesh_number)%gnode_from_knode = 0
    do gnode = 1, ubound(gmesh(gmesh_number)%knode_from_gnode,1)
      knode = gmesh(gmesh_number)%knode_from_gnode(gnode)
      if (knode /= 0) then
        gmesh(gmesh_number)%gnode_from_knode(knode) = gnode
        gmesh(gmesh_number)%ngnodes = gmesh(gmesh_number)%ngnodes + 1
      end if
    end do
  else if (allocated(gmesh(gmesh_number)%gnode_from_knode)) then
    deallocate(gmesh(gmesh_number)%gnode_from_knode)
  end if

! counting of gregions
  gmesh(gmesh_number)%ngregions = 0
  if (allocated(gmesh(gmesh_number)%gregion)) then
    do gregion = 1, ubound(gmesh(gmesh_number)%gregion,1)
      region_number = gmesh(gmesh_number)%gregion(gregion)%region_number
      if (region_number /= 0) gmesh(gmesh_number)%ngregions = gmesh(gmesh_number)%ngregions + 1
    end do
  end if

! counting of gelements and create arrays of gregions for each gelement
  gmesh(gmesh_number)%ngelements = 0
  gmesh(gmesh_number)%ngelements_face = 0
  gmesh(gmesh_number)%ngelements_cell = 0
  if (allocated(gmesh(gmesh_number)%gelement)) then
    do gelement = 1, ubound(gmesh(gmesh_number)%gelement,1)
      if (allocated(gmesh(gmesh_number)%gelement(gelement)%gregions)) deallocate(gmesh(gmesh_number)%gelement(gelement)%gregions)
      i = gmesh(gmesh_number)%gelement(gelement)%icell
      iadd = 0
      if (i /= 0.and.allocated(gmesh(gmesh_number)%gregion)) then
        do gregion = 1, ubound(gmesh(gmesh_number)%gregion,1)
  ! only include gelements that have the same dimension as the gregion (otherwise gmsh pickles)
          region_number = gmesh(gmesh_number)%gregion(gregion)%region_number
          if (region_number == 0) cycle
          if (cell(i)%dimensions /= region(region_number)%dimensions) cycle
  ! also only include non-boundary cells, otherwise there will be confusion when reading in file again
  ! NB: as per original gmsh file, boundary cells are created once the read in has taken place
  ! data can still be associated with them as there gelements for the boundary face elements are defined
          if (cell(i)%type /= 1) cycle
          if (location_in_list(array=cell(i)%region_list,element=region_number) /= 0) then
            iadd = iadd + 1
            call push_integer_array(gmesh(gmesh_number)%gelement(gelement)%gregions,new_element=gregion)
          end if
        end do
      end if
      j = gmesh(gmesh_number)%gelement(gelement)%jface
      jadd = 0
      if (j /= 0.and.allocated(gmesh(gmesh_number)%gregion)) then
        do gregion = 1, ubound(gmesh(gmesh_number)%gregion,1)
          region_number = gmesh(gmesh_number)%gregion(gregion)%region_number
          if (region_number == 0) cycle
          if (face(j)%dimensions /= region(region_number)%dimensions) cycle
          if (location_in_list(array=face(j)%region_list,element=region_number) /= 0) then
            jadd = jadd + 1
            call push_integer_array(gmesh(gmesh_number)%gelement(gelement)%gregions,new_element=gregion)
          end if
        end do
      end if
      if (iadd+jadd == 0) call push_integer_array(gmesh(gmesh_number)%gelement(gelement)%gregions,new_element=0)
      if (i /= 0) gmesh(gmesh_number)%ngelements_cell = gmesh(gmesh_number)%ngelements_cell + max(iadd+jadd,1)
      if (j /= 0) gmesh(gmesh_number)%ngelements_face = gmesh(gmesh_number)%ngelements_face + max(iadd+jadd,1)
      if (i /= 0 .or. j /= 0) gmesh(gmesh_number)%ngelements = gmesh(gmesh_number)%ngelements + max(iadd+jadd,1)
    end do
  end if

  if (debug) then
    write(87,*) 'GMESH: gmesh_number = ',gmesh_number,': basename = ',trim(gmesh(gmesh_number)%basename)
    write(87,*) 'ngelements = ',gmesh(gmesh_number)%ngelements
    write(87,*) 'ngelements_cell = ',gmesh(gmesh_number)%ngelements_cell
    write(87,*) 'ngelements_face = ',gmesh(gmesh_number)%ngelements_face
    if (allocated(gmesh(gmesh_number)%gelement)) then
      do gelement = 1, ubound(gmesh(gmesh_number)%gelement,1)
        write(87,*) 'gelement = ',gelement,': icell = ',gmesh(gmesh_number)%gelement(gelement)%icell, &
        ': jface = ',gmesh(gmesh_number)%gelement(gelement)%jface, &
        ': gregions =',(' ',gmesh(gmesh_number)%gelement(gelement)%gregions(gregion), &
        gregion=1,ubound(gmesh(gmesh_number)%gelement(gelement)%gregions,1))
      end do
    end if
  end if
    
end do

if (debug_sparse) then
  write(*,'(a)') 'INFO: gmeshes:'
  do gmesh_number = 0, ubound(gmesh,1) ! atleast mesh 0 must always be allocated
    formatline = '(a,'//trim(dindexformat(gmesh_number))// &
      ',a,'//trim(dindexformat(gmesh(gmesh_number)%ngnodes))// &
      ',a,'//trim(dindexformat(gmesh(gmesh_number)%ngelements))// &
      ',a,'//trim(dindexformat(gmesh(gmesh_number)%ngelements_cell))// &
      ',a,'//trim(dindexformat(gmesh(gmesh_number)%ngelements_face))// &
      ',a,'//trim(dindexformat(gmesh(gmesh_number)%ngregions))// &
      ',a'//repeat(',a',ubound(gmesh(gmesh_number)%options,1))//')'
    write(*,fmt=formatline) ' gmesh_number = ',gmesh_number,': basename = '//trim(gmesh(gmesh_number)%basename)//': ngnodes = ', &
      gmesh(gmesh_number)%ngnodes,': ngelements = ',gmesh(gmesh_number)%ngelements,': ngelements_cell = ', &
      gmesh(gmesh_number)%ngelements_cell,': ngelements_face = ',gmesh(gmesh_number)%ngelements_face,': ngregions = ', &
      gmesh(gmesh_number)%ngregions, &
      ': options (prioritised) =',(' '//trim(gmesh(gmesh_number)%options(o)),o=1,ubound(gmesh(gmesh_number)%options,1))
  end do
end if

if (debug_sparse) write(*,'(a/80(1h-))') 'subroutine finalise_gmesh'

end subroutine finalise_gmesh

!-----------------------------------------------------------------

subroutine push_gmesh(filename,gmesh_number)

! push new gmesh element onto gmesh array
! right now cannot handle allocated gregions_from_gelement

use general_module
character(len=*), intent(in) :: filename
integer, intent(out), optional :: gmesh_number
integer :: gmesh_number_local
type(gmesh_type), dimension(:), allocatable :: gmesh_tmp
integer :: n

if (allocated(gmesh)) then

! check if mesh is already defined
  do n = 0, ubound(gmesh,1)
    if (trim(gmesh(n)%basename) == trim(basename(filename))) then
      gmesh_number_local = n
      if (trim(gmesh(n)%filename) /= trim(filename)) then
        write(*,'(a)') 'ERROR: the same gmesh basename '//trim(gmesh(n)%basename)// &
          ' has been specified in multiple directories:'//'  this basename must be unique'
        stop
      end if
      if (present(gmesh_number)) gmesh_number = gmesh_number_local
      return
    end if
  end do

! if not already defined then increase gmesh array

  gmesh_number_local = ubound(gmesh,1)
  allocate(gmesh_tmp(0:gmesh_number_local))
  gmesh_tmp(0:gmesh_number_local) = gmesh(0:gmesh_number_local)
  deallocate(gmesh)
  gmesh_number_local = gmesh_number_local + 1
else
! output array should be the first one to be allocated, having index 0
  gmesh_number_local = 0
end if
allocate(gmesh(0:gmesh_number_local))
if (allocated(gmesh_tmp)) then
  gmesh(0:gmesh_number_local-1) = gmesh_tmp(0:gmesh_number_local-1)
  deallocate(gmesh_tmp)
end if

! set data for the new element
gmesh(gmesh_number_local)%filename = trim(filename)
gmesh(gmesh_number_local)%basename = trim(basename(filename))
! add default input and output options
if (gmesh_number_local == 0) then
  call push_character_array(array=gmesh(gmesh_number_local)%options,new_element='noinput')
  call push_character_array(array=gmesh(gmesh_number_local)%options,new_element='output')
else
  call push_character_array(array=gmesh(gmesh_number_local)%options,new_element='input')
  call push_character_array(array=gmesh(gmesh_number_local)%options,new_element='nooutput')
end if
if (present(gmesh_number)) gmesh_number = gmesh_number_local ! return gmesh_number if requested

end subroutine push_gmesh

!-----------------------------------------------------------------

end module gmsh_module

!-----------------------------------------------------------------
