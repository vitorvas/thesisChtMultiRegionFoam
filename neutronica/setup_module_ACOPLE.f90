! file src/setup_module.f90
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
module setup_module

implicit none

private
public setup ! only subroutine accessible from outside the module

! various setup related options

!-----------------------------------------------------------------
contains

!-----------------------------------------------------------------

subroutine setup

! this handles everything prior to a simulation starting

use general_module
use equations_module
use gmsh_module
use kernel_module
use output_module
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'subroutine setup'

if (debug) write(*,*) 'calling setup_gesh'

call setup_gmesh ! create reference gtype_list array which has details of the gmsh element types

if (debug) write(*,*) 'calling allocate_meta_arrays'

call allocate_meta_arrays ! allocate arrays that have meta data about variables

if (debug) write(*,*) 'calling read_constants_file'

call read_constants_file ! read constants file for numerical data, some region and mesh info

if (debug) write(*,*) 'calling setup_mesh'

call setup_mesh ! read in meshes

if (debug) write(*,*) 'calling setup_regions'

call setup_regions ! assign ij indices to each region not setup by mesh

if (debug) write(*,*) 'calling setup_kernels'

!!!! call setup_kernels ! create averaging and differencing kernels !!!! abernal, 13/2/14

if (debug) write(*,*) 'calling setup_vars'

call setup_vars ! do setup for each var variable

if (debug) write(*,*) 'calling finalise_gmesh'

call finalise_gmesh ! finalise the setup of all gmeshes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal. 13/06/2013
!call output_step(action="setup") ! setup the output_step file
!
!if (convergence_details_file) then
!
!  if (debug) write(*,*) 'calling setup_convergence_file'
!
!  call setup_convergence_file ! setup convergence file
!
!end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (debug) write(*,'(a/80(1h-))') 'subroutine setup'

end subroutine setup

!-----------------------------------------------------------------

subroutine setup_convergence_file

use general_module
character(len=100) :: filename
integer :: error

! open convergence output file if requested
filename = "output/convergence_details.txt"
open(fconverge,file=trim(filename),access='append',iostat=error)
if (error /= 0) call error_stop('problem opening file '//trim(filename))

end subroutine setup_convergence_file

!-----------------------------------------------------------------

subroutine read_constants_file

! here we read in all sorts of simulation data from the constants.in file, but
! not the actual constants values (which is done in read_constants)

use general_module
use gmsh_module
integer :: error, cut, n, m, gmesh_number
character(len=1000) :: textline, keyword, name, formatline, region_name, location, options
character(len=4) :: centring
character(len=1) :: delimiter
character(len=100) :: option
real :: versiontmp
logical :: existing
logical, parameter :: debug = .false.
                  
if (debug) write(*,'(80(1h+)/a)') 'subroutine read_constants_file'

! push default mesh (index 0) which will include everything
call push_gmesh(filename='output/output.msh')

! read in constants values from constants.in file

write(*,'(a)') "INFO: reading simulation information from constants file "//trim(constants_file)
open(unit=fconstants,file=trim(constants_file),status='old',iostat=error)
if (error /= 0) call error_stop('problem opening constants file '//trim(constants_file))

fileloop: do
  read(fconstants,'(a)',iostat=error) textline
  if (error /= 0) exit fileloop ! reached end of file
  call remove_comments(textline) ! remove comments from line
  cut=scan(textline,' ') ! split input line at first space
  if (cut<=1) cycle fileloop ! if there is space at first character or string is empty, cycle
  keyword=textline(1:cut-1)
  textline=adjustl(textline(cut:len(textline)))
! if (debug) write(*,'(a)') 'keyword = '//trim(keyword)//': textline = '//trim(textline)//''
  if (debug) write(*,'(a)') 'keyword = '//trim(keyword)

!---------------
! name
  if (trim(keyword) == 'FACE_REGION'.or.trim(keyword) == 'CELL_REGION') then ! read in variable names
    cut=scan(textline,'>') ! find end of name string
    if (textline(1:1) /= '<'.or.cut<2) call error_stop('region name in '//trim(keyword)// &
      ' in constants file incorrectly specified on line:'//trim(textline))
    name = textline(1:cut) ! name now includes <>
    textline=adjustl(textline(cut+1:len(textline)))
  end if

!---------------
! location

  if (trim(keyword) == 'FACE_REGION' .or. trim(keyword) == 'CELL_REGION') then
    read(textline,*,iostat=error) location ! with no format specification string in constants file should be quotted
    if (error /= 0) location = "" ! if no location is specified then we are just specifying region centring for gmsh region
    if (trim(keyword) == 'FACE_REGION') then
      centring = 'face'
    else
      centring = 'cell'
    end if
    m = region_number_from_name(name=name,location=location,centring=centring,existing=existing,creatable=.true.)
    if (m == 0.or.existing) call error_stop('allocation of region = '//trim(region_name)// &
      ' failed, most likely because the region has already been allocated previously in the constants file')
    formatline = '(a,'//trim(dindexformat(m))//',a)'
    if (trim(location) == "") then
      write(*,fmt=formatline) 'INFO: centring for region '//trim(region(m)%name)//' has been specified from constants file: '// &
        'region_number = ',m,': centring = '//trim(centring)
    else
      write(*,fmt=formatline) 'INFO: region '//trim(region(m)%name)//' has been created from constants file: region_number =  ',m, &
        ': centring = '//trim(centring)//': location = '//trim(region(m)%location)
    end if
  end if

!---------------
! mesh info

  if (trim(keyword) == 'READ_GMSH'.or.trim(keyword) == 'MSH_FILE') then
    if (trim(keyword) == 'READ_GMSH') write(*,*) 'WARNING: keyword READ_GMSH should be replaced by keyword MSH_FILE'
! extract gmsh filename, either delimited by or not
    if (textline(1:1) == '"' .or. textline(1:1) == "'") then
      delimiter = textline(1:1)
      textline = trim(textline(2:len(textline)))
      cut=scan(textline,delimiter) ! split input line at next delimiter
    else
      cut=scan(textline,' ') ! split input line at first space
    end if
    if (cut<=1) call error_stop('problem reading in name of gmsh mesh file from line '//textline) ! if there is space at first character or string is empty, then filename not found
    call push_gmesh(filename=textline(1:cut-1),gmesh_number=gmesh_number)
    textline=trim(adjustl(textline(cut+1:len(textline))))
! extract any options and add to the start of the options array
    do
      cut=scan(textline,',')
      if (cut == 1) then
        textline = trim(adjustl(textline(2:len(textline))))
      else if (cut > 1) then
        option = trim(adjustl(textline(1:cut-1)))
        if (option /= "") call unshift_character_array(array=gmesh(gmesh_number)%options,new_element=option)
        textline = trim(adjustl(textline(cut+1:len(textline))))
      else
        option = trim(adjustl(textline))
        if (option /= "") call unshift_character_array(array=gmesh(gmesh_number)%options,new_element=option)
        exit
      end if
    end do
    if (allocated(gmesh(gmesh_number)%options)) then
      write(options,'(100(a))') (' '//trim(gmesh(gmesh_number)%options(n)),n=1,ubound(gmesh(gmesh_number)%options,1))
    else
      options = " (none)"
    end if
    formatline = '(a,'//trim(dindexformat(gmesh_number))//',a)'
    write(*,fmt=formatline) 'INFO: gmesh created from constants file: gmesh_number = ',gmesh_number,': fullname = '// &
      trim(gmesh(gmesh_number)%filename)//': basename = '//trim(gmesh(gmesh_number)%basename)// &
      ': current (prioritised) options ='//trim(options)
  end if

!---------------
! dimensions

  if (trim(keyword) == 'DIMENSIONS') &
    write(*,'(a)') "WARNING: global dimensions statement now ignored - dimensions are specific to each gmesh element"

!---------------
! newtrestol

  if (trim(keyword) == 'NEWTRESTOL') then
    read(textline,*,iostat=error) newtrestol
    if (error /= 0) call error_stop('problem reading in newton loop tolerance from line '//textline)
    write(*,'(a,g14.6)') 'INFO: newtrestol = ',newtrestol
  end if

!---------------
! newtstepmax

  if (trim(keyword) == 'NEWTSTEPMAX') then
    read(textline,*,iostat=error) newtstepmax
    if (error /= 0) call error_stop('problem reading in maximum number of newton steps from line '//textline)
    formatline = '(a,'//trim(dindexformat(newtstepmax))//')'
    write(*,fmt=formatline) 'INFO: newtstepmax = ',newtstepmax
  end if

!---------------
! timestepmax

  if (trim(keyword) == 'TIMESTEPMAX') then
    read(textline,*,iostat=error) timestepmax
    if (error /= 0) call error_stop('problem reading in maximum number of time steps from line '//textline)
    formatline = '(a,'//trim(dindexformat(timestepmax))//')'
    write(*,fmt=formatline) 'INFO: timestepmax = ',timestepmax
  end if

!---------------
! timestepout

  if (trim(keyword) == 'TIMESTEPOUT') then
    read(textline,*,iostat=error) timestepout
    if (error /= 0) call error_stop('problem reading in number of time steps between output from line '//textline)
    formatline = '(a,'//trim(dindexformat(timestepout))//')'
    write(*,fmt=formatline) 'INFO: timestepout = ',timestepout
  end if

!---------------
! version

  if (trim(keyword) == 'VERSION') then
    read(textline,*,iostat=error) versiontmp
    if (abs(version - versiontmp) > 2.*tiny(1.e0)) then
      formatline = '(a,f4.2,a,f4.2,a)'
      if (versiontmp < minimum_version) then
        write(*,fmt=formatline) 'ERROR: unsafe version mismatch between constants.in (',versiontmp,') and the current (', &
          version,'): the language syntax has changed between these versions'
        stop
      else
        write(*,fmt=formatline) 'WARNING: (safe) version mismatch between constants.in (',versiontmp,') and the current (', &
          version,'): additional language features are now available'
      end if
    end if
  end if

!---------------
! linear solver type

  if (trim(keyword) == 'LINEAR_SOLVER') then
    read(textline,*,iostat=error) linear_solver
    if (error /= 0) call error_stop('problem reading in linear solver type from line '//textline)
    write(*,'(a)') 'INFO: linear_solver = '//trim(linear_solver)
  end if

!---------------

end do fileloop

close(fconstants)

if (debug) write(*,'(a/80(1h-))') 'subroutine read_constants_file'

end subroutine read_constants_file

!-----------------------------------------------------------------

subroutine setup_mesh

! here we do all sorts of mesh setting
! set types for all cells, faces and nodes
! check on icell associations for faces
! set i,j,k boundary and domain and numbers
! set x for faces and nodes

use general_module
use gmsh_module
use input_module
integer :: i, ii, j, jj, k, kk, ncells, gmesh_number, n, error
type(cell_type) :: default_cell
double precision, dimension(totaldimensions) :: tangc, normc
character(len=1000) :: formatline, filename
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'subroutine setup_mesh'

!-------------------------------------
! set initial empty element numbers
ktotal = 0
jtotal = 0
itotal = 0

!-------------------------------------
! read in the gmesh data
if (debug) write(*,*) 'reading in gmesh data'

call read_gmesh(contents='mesh')

!-------------------------------------
if (debug) write(*,*) 'doing generic mesh setup'

! run through cells checking icells and possibly setting icell(2) indicies
do i = 1, itotal
  cell(i)%type = 1 ! domain cell
  jj_loop: do jj = 1, ubound(cell(i)%jface,1)
    j = cell(i)%jface(jj)
    if (face(j)%icell(1) == 0) then
      write(*,*) 'ERROR: a face that is adjacent to a cell has no icell=1 attached: j = ',j,': face = ',trim(print_face(face(j)))
      stop
    end if
    if (face(j)%icell(1) == i) cycle jj_loop
    if (ubound(face(j)%icell,1) < 2) call push_integer_array(array=face(j)%icell,new_element=0)
    if (face(j)%icell(2) == i) cycle jj_loop
    if (face(j)%icell(2) /= 0) then
      write(*,*) 'ERROR: a face has more than 2 cells attached to it: j = ',j,': face = ',trim(print_face(face(j)))
      stop
    end if
    face(j)%icell(2) = i
  end do jj_loop
end do
      
!-------------------------------------
! run through nodes setting type to default
do k = 1, ktotal
  node(k)%type = 1
end do

!-------------------------------------
! now assign face type and create any boundary cells
idomain = itotal
jboundary = 0
jdomain = 0

if (debug) write(*,*) 'creating boundary cells'

do j = 1, jtotal
  if (face(j)%icell(1) == 0) then
    write(*,*) 'ERROR: a face has no icell=1 attached to it: j = ',j,': face = ',trim(print_face(face(j)))
    write(*,*) 'A possibility is that there is a problem with the dimensions in a .msh file or in a region definition.'
    write(*,*) 'Check that the dimensions in each msh file (physical entities) are consistent, and that the maximum dimension '// &
      'of each domain is listed (remember that each domain region must be given a physical entity name in gmsh).'
    write(*,*) 'Also check that any domains used in each msh file that have a dimension less than that of the maximum domain '// &
      'have their centrings explicity set in constants.in'
    stop
  end if
  if (ubound(face(j)%icell,1) < 2) call push_integer_array(array=face(j)%icell,new_element=0)
  if (face(j)%icell(2) == 0) then
    face(j)%type = 2 ! boundary face
    jboundary = jboundary + 1

! create new boundary cell, storing in default_cell
    call reset_cell(default_cell)
    default_cell%type = 2
    default_cell%dimensions = face(j)%dimensions
    default_cell%gtype = face(j)%gtype ! adopt gmsh element type - must be set
    call resize_integer_array(keep_data=.false.,array=default_cell%jface,new_size=1)
    default_cell%jface(1) = j
    call copy_integer_array(original=face(j)%knode,copy=default_cell%knode)

    itotal = itotal + 1
    if (ubound(cell,1) < itotal) call resize_cell_array(change=int(ktotal/10+1))
    call set_cell(cell_to_set=cell(itotal),new_value=default_cell)
    face(j)%icell(2) = itotal

!   if (debug) write(*,*) 'creating new boundary cell: i = ',itotal,': jface = ',default_cell%jface

! set node icell references for the new cell
    call add_icell_to_nodes(itotal,default_cell%knode)

! nodes associated with boundary faces are also boundary faces
    do kk = 1, ubound(face(j)%knode,1)
      k = face(j)%knode(kk)
      node(k)%type = 2
    end do

  else
    face(j)%type = 1 ! domain face
    jdomain = jdomain + 1
  end if

end do

iboundary = itotal - idomain

! resize cell arrays
call resize_cell_array(new_size=itotal)

! run through gmeshes checking that any boundary faces or cells have a corresponding boundary cell or face
do gmesh_number = 0, ubound(gmesh,1)
  if (allocated(gmesh(gmesh_number)%gelement)) then
    do n = 1, ubound(gmesh(gmesh_number)%gelement,1)
      i = gmesh(gmesh_number)%gelement(n)%icell
      j = gmesh(gmesh_number)%gelement(n)%jface
      if (j /= 0.and.i == 0) then
        if (face(j)%type == 2) gmesh(gmesh_number)%gelement(n)%icell = face(j)%icell(2)
      end if
      if (i /= 0.and.j == 0) then
        if (cell(i)%type == 2) gmesh(gmesh_number)%gelement(n)%jface = cell(i)%jface(1)
      end if
    end do
  end if
end do

!-------------------------------------
! run through nodes again counting types

kboundary = 0
do k = 1, ktotal
  if (node(k)%type == 2) kboundary = kboundary + 1
end do
kdomain = ktotal - kboundary

!-------------------------------------
! run through faces increasing icell

do j = 1, jtotal
! first two elements (the immediate cell neighbours) are already present

  if (.true.) then
! second elements are from the surrounding nodes
    kk_loop: do kk = 1, ubound(face(j)%knode,1)
      k = face(j)%knode(kk)
      ii_loop: do ii = 1, ubound(node(k)%icell,1)
        if (location_in_list(array=face(j)%icell,element=node(k)%icell(ii)) == 0) &
          call push_integer_array(array=face(j)%icell,new_element=node(k)%icell(ii))
      end do ii_loop
    end do kk_loop
  else
! new method based solely on expansion through number of degrees of separation
! if this works for kernels then should integrate it in kernel mask routine instead
    do n = 2, 6 ! loop through number of expansion stages
      ncells = ubound(face(j)%icell,1) ! save list of cells that are included at the start
      do ii = 1, ncells
        i = face(j)%icell(ii)
        do jj = 1, ubound(cell(i)%jface,1)
          if (location_in_list(array=face(j)%icell,element=face(cell(i)%jface(jj))%icell(1)) == 0) then
            call push_integer_array(array=face(j)%icell,new_element=face(cell(i)%jface(jj))%icell(1))
          else if (location_in_list(array=face(j)%icell,element=face(cell(i)%jface(jj))%icell(2)) == 0) then
            call push_integer_array(array=face(j)%icell,new_element=face(cell(i)%jface(jj))%icell(2))
          end if
        end do
      end do
    end do
  end if

end do

!-------------------------------------
! run through cells setting icell (surrounding cells)
! this list includes all cells that share a common node with the central cell,
! included in the order of:
! 1 = central cell: 2-ubound(cell%jface,1)+1 = cells sharing a face: rest are the
! cells that only share a node

do i = 1, itotal

! first element is central cell
  call resize_integer_array(keep_data=.false.,array=cell(i)%icell,new_size=1)
  cell(i)%icell(1) = i

! second elements are from the surrounding faces
  do jj = 1, ubound(cell(i)%jface,1)
    j = cell(i)%jface(jj)
    do ii = 1, 2
      if (face(j)%icell(ii) /= 0 .and. face(j)%icell(ii) /= i) then
        call push_integer_array(array=cell(i)%icell,new_element=face(j)%icell(ii))
      end if
    end do
  end do

! third elements are from the surrounding nodes
  do kk = 1, ubound(cell(i)%knode,1)
    k = cell(i)%knode(kk)
    ii2_loop: do ii = 1, ubound(node(k)%icell,1)
      if (location_in_list(array=cell(i)%icell,element=node(k)%icell(ii)) == 0) &
        call push_integer_array(array=cell(i)%icell,new_element=node(k)%icell(ii))
    end do ii2_loop
  end do

end do

!-------------------------------------
! run through faces and cells checking that number of components is consistent with dimensions

if (debug) write(*,*) 'checking mesh consistency'

do i = 1, itotal
  if (cell(i)%type == 2) then
    if (ubound(cell(i)%jface,1) /= 1 ) then
      write(*,*) 'ERROR: boundary cell does not have 1 face'
      write(*,'(a,i6,a)') 'cell: i = ',i,': ',trim(print_cell(cell(i)))
      stop
    end if
    if (cell(i)%dimensions == 2) then
      if (ubound(cell(i)%knode,1) < 3 ) then
        write(*,*) 'ERROR: 2d boundary cell has less than 3 nodes'
        write(*,'(a,i6,a)') 'cell: i = ',i,': ',trim(print_cell(cell(i)))
        stop
      end if
    else if (cell(i)%dimensions == 1) then
      if (ubound(cell(i)%knode,1) /= 2 ) then
        write(*,*) 'ERROR: 1d boundary cell does not have 2 nodes'
        write(*,'(a,i6,a)') 'cell: i = ',i,': ',trim(print_cell(cell(i)))
        stop
      end if
    else if (cell(i)%dimensions == 0) then
      if (ubound(cell(i)%knode,1) /= 1 ) then
        write(*,*) 'ERROR: 0d boundary cell does not have 1 node'
        write(*,'(a,i6,a)') 'cell: i = ',i,': ',trim(print_cell(cell(i)))
        stop
      end if
    else
      write(*,*) 'ERROR: boundary cell has incorrect dimensions/nodes'
      write(*,'(a,i6,a)') 'cell: i = ',i,': ',trim(print_cell(cell(i)))
      stop
    end if
  else
    if (cell(i)%dimensions == 3) then
      if (ubound(cell(i)%jface,1) < 4 .or. ubound(cell(i)%knode,1) < 4) then
        write(*,*) 'ERROR: 3d domain cell has less than 4 faces/nodes'
        write(*,'(a,i6,a)') 'cell: i = ',i,': ',trim(print_cell(cell(i)))
        stop
      end if
    else if (cell(i)%dimensions == 2) then
      if (ubound(cell(i)%jface,1) < 3 .or. ubound(cell(i)%knode,1) < 3) then
        write(*,*) 'ERROR: 2d domain cell has less than 3 faces/nodes'
        write(*,'(a,i6,a)') 'cell: i = ',i,': ',trim(print_cell(cell(i)))
        stop
      end if
    else if (cell(i)%dimensions == 1) then
      if (ubound(cell(i)%jface,1) /= 2  .or. ubound(cell(i)%knode,1) /= 2) then
        write(*,*) 'ERROR: 1d domain cell does not have 2 faces/nodes'
        write(*,'(a,i6,a)') 'cell: i = ',i,': ',trim(print_cell(cell(i)))
        stop
      end if
    else
      write(*,*) 'ERROR: domain cell has incorrect dimensions'
      write(*,'(a,i6,a)') 'cell: i = ',i,': ',trim(print_cell(cell(i)))
      stop
    end if
    do jj = 1, ubound(cell(i)%jface,1)
      j = cell(i)%jface(jj)
      if (face(j)%dimensions /= cell(i)%dimensions - 1) then
        write(*,*) 'ERROR: face has dimensions which are inconsistent with adjacent cell'
        write(*,'(a,i6,a)') 'cell: i = ',i,': ',trim(print_cell(cell(i)))
        write(*,'(a,i6,a)') 'face: j = ',j,': ',trim(print_face(face(j)))
        stop
      end if
    end do
  end if
end do

do j = 1, jtotal
  if (face(j)%dimensions == 2) then
    if (ubound(face(j)%knode,1) < 3 ) then
      write(*,*) 'ERROR: 2d face has less than 3 nodes'
      write(*,'(a,i6,a)') 'face: j = ',j,': ',trim(print_face(face(j)))
      stop
    end if
  else if (face(j)%dimensions == 1) then
    if (ubound(face(j)%knode,1) /= 2 ) then
      write(*,*) 'ERROR: 1d face does not have 2 nodes'
      write(*,'(a,i6,a)') 'face: j = ',j,': ',trim(print_face(face(j)))
      stop
    end if
  else if (face(j)%dimensions == 0) then
    if (ubound(face(j)%knode,1) /= 1 ) then
      write(*,*) 'ERROR: 0d face does not have 1 node'
      write(*,'(a,i6,a)') 'face: j = ',j,': ',trim(print_face(face(j)))
      stop
    end if
  else
    write(*,*) 'ERROR: face has incorrect dimensions/nodes'
    write(*,'(a,i6,a)') 'face: j = ',j,': ',trim(print_face(face(j)))
    stop
  end if
end do

!-------------------------------------
! run through cells calculating vol and x (centre) 

if (debug) write(*,*) 'calculating mesh cell geometries'

do i = 1, itotal
  if (cell(i)%dimensions == 3) then
    call find_3d_geometry(jface=cell(i)%jface,volume=cell(i)%vol,centre=cell(i)%x)
  else if (cell(i)%dimensions == 2) then
    call find_2d_geometry(knode=cell(i)%knode,area=cell(i)%vol,centre=cell(i)%x)
  else if (cell(i)%dimensions == 1) then
    cell(i)%vol = distance( node(cell(i)%knode(1))%x , node(cell(i)%knode(2))%x )
    cell(i)%x = ( node(cell(i)%knode(1))%x + node(cell(i)%knode(2))%x )/2.d0
  else if (cell(i)%dimensions == 0) then
    cell(i)%vol = 1.d0
    cell(i)%x = node(cell(i)%knode(1))%x
  end if
end do

!-------------------------------------
! run through faces calculating area, x (centre) and norm

if (debug) write(*,*) 'calculating mesh face geometries'

do j = 1, jtotal
  if (face(j)%dimensions == 2) then
    call find_2d_geometry(knode=face(j)%knode,area=face(j)%area,norm=face(j)%norm,centre=face(j)%x)
  else if (face(j)%dimensions == 1) then
    face(j)%area = distance( node(face(j)%knode(1))%x , node(face(j)%knode(2))%x )
    tangc = cell(face(j)%icell(2))%x - cell(face(j)%icell(1))%x ! vector that is in same plane as the adjacent cells
    normc = cross_product( tangc , node(face(j)%knode(2))%x - node(face(j)%knode(1))%x )
    face(j)%norm(:,1) = cross_product( node(face(j)%knode(2))%x - node(face(j)%knode(1))%x , normc )
    call normalise_vector(face(j)%norm(:,1))
    face(j)%norm(:,2) = node(face(j)%knode(2))%x - node(face(j)%knode(1))%x ! vector that is tangent to line
    call normalise_vector(face(j)%norm(:,2))
    face(j)%norm(:,3) = cross_product( face(j)%norm(:,1) , face(j)%norm(:,2) ) ! third vector which points out of solution plane
    call normalise_vector(face(j)%norm(:,3))
    face(j)%x = ( node(face(j)%knode(1))%x + node(face(j)%knode(2))%x )/2.d0
  else if (face(j)%dimensions == 0) then
    face(j)%area = 1.d0
    face(j)%norm(:,1) = cell(face(j)%icell(2))%x - cell(face(j)%icell(1))%x
    call normalise_vector(face(j)%norm(:,1))
! from orthogonal vectors select one which has smallest dot product with first
    face(j)%norm(:,2) = [1.d0, 0.d0, 0.d0]
    face(j)%norm(:,3) = [0.d0, 1.d0, 0.d0]
    if (abs(dot_product(face(j)%norm(:,2),face(j)%norm(:,1))) < abs(dot_product(face(j)%norm(:,3),face(j)%norm(:,1))) ) then
      face(j)%norm(:,3) = cross_product( face(j)%norm(:,1) , face(j)%norm(:,2) )
      face(j)%norm(:,2) = cross_product( face(j)%norm(:,1) , face(j)%norm(:,3) )
    else
      face(j)%norm(:,2) = cross_product( face(j)%norm(:,1) , face(j)%norm(:,3) )
      face(j)%norm(:,3) = cross_product( face(j)%norm(:,1) , face(j)%norm(:,2) )
    end if
    call normalise_vector(face(j)%norm(:,2))
    call normalise_vector(face(j)%norm(:,3))
    face(j)%x = node(face(j)%knode(1))%x
  end if
  face(j)%dx = distance( cell(face(j)%icell(1))%x, cell(face(j)%icell(2))%x ) 
end do

!------------------------------------------
! write out mesh details if requested

write(*,'(a)') 'INFO: number of mesh elements::'
formatline = '(3(a,'//trim(indexformat)//'))'
write(*,fmt=formatline) ' NODES: ktotal = ',ktotal,': kdomain = ',kdomain,': kboundary = ',kboundary
formatline = '(3(a,'//trim(indexformat)//'))'
write(*,fmt=formatline) ' FACES: jtotal = ',jtotal,': jdomain = ',jdomain,': jboundary = ',jboundary
formatline = '(3(a,'//trim(indexformat)//'))'
write(*,fmt=formatline) ' CELLS: itotal = ',itotal,': idomain = ',idomain,': iboundary = ',iboundary

if (mesh_details_file) then
  if (debug) write(*,*) 'writing mesh details to mesh_details.txt file'

  filename = "output/mesh_details.txt"
  open(fdetail,file=trim(filename),status='replace',iostat=error)
  if (error /= 0) call error_stop('problem opening file '//trim(filename))

  write(fdetail,'(a)') 'MESH DETAILS:'
  formatline = '(3(a,'//trim(indexformat)//'))'
  write(fdetail,fmt=formatline) 'NODES: ktotal = ',ktotal,': kdomain = ',kdomain,': kboundary = ',kboundary
  formatline = '(a,'//trim(indexformat)//',a,a)'
  do k = 1,ktotal
    write(fdetail,fmt=formatline) 'node: k = ',k,': ',trim(print_node(node(k)))
  end do
  formatline = '(3(a,'//trim(indexformat)//'))'
  write(fdetail,fmt=formatline) 'FACES: jtotal = ',jtotal,': jdomain = ',jdomain,': jboundary = ',jboundary
  formatline = '(a,'//trim(indexformat)//',a,a)'
  do j = 1,jtotal
    write(fdetail,fmt=formatline) 'face: j = ',j,': ',trim(print_face(face(j)))
  end do
  formatline = '(3(a,'//trim(indexformat)//'))'
  write(fdetail,fmt=formatline) 'CELLS: itotal = ',itotal,': idomain = ',idomain,': iboundary = ',iboundary
  formatline = '(a,'//trim(indexformat)//',a,a)'
  do i = 1,itotal
    write(fdetail,fmt=formatline) 'cell: i = ',i,': ',trim(print_cell(cell(i)))
  end do

  close(fdetail)
end if

if (debug) write(*,'(a/80(1h-))') 'subroutine setup_mesh'

end subroutine setup_mesh

!-----------------------------------------------------------------

subroutine setup_regions

! here we setup the regions by finding the i or j indices for each

use general_module
integer :: error, m, i, j, n, cut, nregion, ijregion, nsregion, ns, nscompound, ii, jj, ij, k, l, ijtotal
double precision :: tmp, tmpmax
double precision, dimension(totaldimensions) :: x, xmin, xmax ! a single location
character(len=1000) :: keyword, name, formatline, location, aregion, region_list, filename, geometry
character(len=4) :: centring
character(len=1) :: rsign
logical :: existing, debug_sparse = .false.!!!!.true.
logical, parameter :: debug = .false.

if (debug) debug_sparse = .true.
                  
if (debug) write(*,'(80(1h+)/a)') 'subroutine setup_regions'

! create system generated regions

location='SYSTEM'
centring='cell'

m = region_number_from_name(name='<all cells>',location=location,centring=centring,existing=existing,creatable=.true.)
if (existing) call error_stop("the system region = "//trim(name)// &
  " has already been set:- this name should not be used user-defined regions")
allocate(region(m)%ij(itotal))
do n = 1, itotal
  region(m)%ij(n) = n
end do
  
m = region_number_from_name(name='<domain>',location=location,centring=centring,existing=existing,creatable=.true.)
if (existing) call error_stop("the system region = "//trim(name)// &
  " has already been set:- this name should not be used user-defined regions")
allocate(region(m)%ij(idomain))
n = 0
do i = 1, itotal
  if (cell(i)%type == 1) then
    n = n + 1
    region(m)%ij(n) = i
  end if
end do
  
m = region_number_from_name(name='<boundary cells>',location=location,centring=centring,existing=existing,creatable=.true.)
if (existing) call error_stop("the system region = "//trim(name)// &
  " has already been set:- this name should not be used user-defined regions")
allocate(region(m)%ij(iboundary))
n = 0
do i = 1, itotal
  if (cell(i)%type == 2) then
    n = n + 1
    region(m)%ij(n) = i
  end if
end do
  
centring='face'

m = region_number_from_name(name='<all faces>',location=location,centring=centring,existing=existing,creatable=.true.)
if (existing) call error_stop("the system region = "//trim(name)// &
  " has already been set:- this name should not be used user-defined regions")
allocate(region(m)%ij(jtotal))
do n = 1, jtotal
  region(m)%ij(n) = n
end do
  
m = region_number_from_name(name='<domain faces>',location=location,centring=centring,existing=existing,creatable=.true.)
if (existing) call error_stop("the system region = "//trim(name)// &
  " has already been set:- this name should not be used user-defined regions")
allocate(region(m)%ij(jdomain))
n = 0
do j = 1, jtotal
  if (face(j)%type == 1) then
    n = n + 1
    region(m)%ij(n) = j
  end if
end do
  
m = region_number_from_name(name='<boundaries>',location=location,centring=centring,existing=existing,creatable=.true.)
if (existing) call error_stop("the system region = "//trim(name)// &
  " has already been set:- this name should not be used user-defined regions")
allocate(region(m)%ij(jboundary))
n = 0
do j = 1, jtotal
  if (face(j)%type == 2) then
    n = n + 1
    region(m)%ij(n) = j
  end if
end do
  
! now run through remaining regions doing any user-defined locating

do m=1,ubound(region,1)

! define the keyword which in this context means the type of region creation method
  keyword = 'UNKNOWN'
  if (region(m)%location(1:4) == "GMSH") then
    keyword = 'GMSH'
  else if (region(m)%location(1:2) == "AT") then
    keyword = 'AT'
  else if (region(m)%location(1:6) == "WITHIN") then
    keyword = 'WITHIN'
  else if (region(m)%location(1:11) == "BOUNDARY OF") then
    keyword = 'BOUNDARY OF'
  else if (region(m)%location(1:9) == "DOMAIN OF") then
    keyword = 'DOMAIN OF'
  else if (region(m)%location(1:15) == "ASSOCIATED WITH") then
    keyword = 'ASSOCIATED WITH'
  else if (region(m)%location(1:8) == "COMPOUND") then
    keyword = 'COMPOUND'
  end if

!---------------------
! a user defined region from the constants file that is a single point

  if (trim(keyword) == "AT") then

! if region is already defined then it is overwritten
    if (allocated(region(m)%ij)) then
      write(*,'(a/a)') &
        "NOTE: an AT region operator is acting on region "//trim(region(m)%name)// &
        " that already contains an element:", &
        " the previous element will be overwritten with the new"
      region(m)%ij(1) = 0
    else
      allocate(region(m)%ij(1))
    end if

    read(region(m)%location(3:1000),*,iostat=error) x
    if (error /= 0) call error_stop('AT location for region '//trim(region(m)%name)//' is not understood as a single point')

    tmpmax = 1.d+20
    if (region(m)%centring == "cell") then
      do i=1,itotal
        tmp = distance(x,cell(i)%x)
        if (tmp < tmpmax) then
          region(m)%ij(1) = i
          tmpmax = tmp
        end if
      end do
    else 
      do j=1,jtotal
        tmp = distance(x,face(j)%x)
        if (tmp < tmpmax) then
          region(m)%ij(1) = j
          tmpmax = tmp
        end if
      end do
    end if

!---------------------
! a user defined region from the constants file that is any elements within a geometry 

  else if (trim(keyword) == "WITHIN") then

! if region is already defined then it is overwritten
    if (allocated(region(m)%ij)) then
      write(*,'(a/a)') &
        "NOTE: an AT region operator is acting on region "//trim(region(m)%name)// &
        " that already contains an element:", &
        " the previous element will be overwritten with the new"
      deallocate(region(m)%ij)
    end if

! find the geometry type
    location = adjustl(region(m)%location(8:1000)//repeat(' ',7))
    cut = scan(location,' ')
    read(location(1:cut-1),*,iostat=error) geometry
    if (error /= 0) call error_stop('WITHIN geometry type for region '//trim(region(m)%name)//' is not understood')
    location = adjustl(location(cut+1:1000)//repeat(' ',cut))

    if (trim(geometry) == 'BOX') then
      if (debug) write(*,*) 'found WITHIN BOX geometry with location points: '//trim(location)

      read(location,*,iostat=error) xmin, xmax
      if (error /= 0) call error_stop('WITHIN BOX location for region '//trim(region(m)%name)// &
        ' is not understood as two three dimensional points')

! check that points are in min and max order and otherwise reorder
      do l = 1, 3
        if (xmin(l) > xmax(l)) then
          write(*,'(a,i1,a)') 'WARNING: dimension ',l,' of the points that define the BOX geometry in region '// &
            trim(region(m)%name)//' were incorrectly ordered: they should be in the order of the minimum coordinate values '// &
            'in each dimension, followed by the maximum coordinate values in each dimension.'
          write(*,'(2(a,i1,a,g14.6))') ' xmin(',l,') = ',xmin(l),' xmax(',l,') = ',xmax(l)
          write(*,*) 'These values will be swapped.'
          tmp = xmax(l)
          xmax(l) = xmin(l)
          xmin(l) = tmp
        end if
      end do

      if (region(m)%centring == "cell") then
        ijtotal = itotal
      else
        ijtotal = jtotal
      end if

      ij_loop: do ij=1,ijtotal
        if (region(m)%centring == "cell") then
          x = cell(ij)%x
        else
          x = face(ij)%x
        end if
        do l = 1, 3
          if ((x(l)-xmin(l))*(xmax(l)-x(l)) < 0.d0) cycle ij_loop
        end do
        call push_integer_array(array=region(m)%ij,new_element=ij)
      end do ij_loop

    else
      call error_stop('WITHIN geometry type of '//trim(geometry)//' for region '//trim(region(m)%name)//' is not understood')
    end if

    if (.not.allocated(region(m)%ij)) write(*,'(a)') 'WARNING: no elements were found within the WITHIN '//trim(geometry)// &
      ' region of '//trim(region(m)%name)

!---------------------
! a new region composed of the boundary to another region, or similar type of related domain

  else if (trim(keyword) == "BOUNDARY OF".or.trim(keyword) == "DOMAIN OF".or.trim(keyword) == "ASSOCIATED WITH") then

! note: this may mean that the region is already defined but defining twice won't hurt if the definition is the same
    if (allocated(region(m)%ij)) then
      write(*,'(a/a)') &
        "NOTE: a "//trim(keyword)//" region operator is acting on region "//trim(region(m)%name)// &
        " that already contains some elements:", &
        " the previous elements will be overwritten with the new"
      deallocate(region(m)%ij)
    end if

! check centring of requested region
    if (region(m)%centring /= 'face'.and.region(m)%centring /= 'cell') call &
      error_stop('incorrect centring for requested '//trim(keyword)//' region '//trim(region(m)%name))

! find constitutent region around which to find boundary
    aregion = adjustl(trim(region(m)%location(len(trim(keyword))+1:1000)))
    
    nregion = region_number_from_name(name=aregion,existing=existing,creatable=.false.)
! check that region exists and that centring is consistent
    if (.not.existing) call error_stop("region "//trim(aregion)//" which is specified in "//trim(keyword)//" operator for "// &
      trim(region(m)%name)//" is not found")
    if (nregion == 0) call error_stop("problem with region "//trim(aregion)//" in "//trim(keyword)//" region "// &
      trim(region(m)%name))
    if (region(nregion)%centring /= 'face'.and.region(nregion)%centring /= 'cell') &
      call error_stop('incorrect centring for '//trim(keyword)//' consitutent region '//trim(region(nregion)%name))

    do nsregion = 1,ubound(region(nregion)%ij,1)  ! loop through all ij indices in constituent region
      ijregion = region(nregion)%ij(nsregion) ! NB, ijregion and region(m)%ij may actually be i or j values depending on centring
      if (region(nregion)%centring == 'cell') then
        if (region(m)%centring == 'cell') then
! create cell region from cell region
          do ii = 1, ubound(cell(ijregion)%jface,1)+1 ! loop around cells that border cells
            i = cell(ijregion)%icell(ii)
            if (trim(keyword) == "BOUNDARY OF".and.cell(i)%type /= 2) cycle
            if (trim(keyword) == "DOMAIN OF".and.cell(i)%type /= 1) cycle
            if (location_in_list(array=region(m)%ij,element=i) == 0) call push_array(array=region(m)%ij,new_element=i)
          end do
        else
! create face region from cell region
          do jj = 1, ubound(cell(ijregion)%jface,1)
            j = cell(ijregion)%jface(jj)
            if (trim(keyword) == "BOUNDARY OF".and.face(j)%type /= 2) cycle
            if (trim(keyword) == "DOMAIN OF".and.face(j)%type /= 1) cycle
            if (location_in_list(array=region(m)%ij,element=j) == 0) call push_array(array=region(m)%ij,new_element=j)
          end do
        end if
      else
! create cell region from face region
! BOUNDARY OF will pick out boundary cells coincident with boundary faces
        if (region(m)%centring == 'cell') then
          do ii = 1, 2
            i = face(ijregion)%icell(ii)
            if (trim(keyword) == "BOUNDARY OF".and.cell(i)%type /= 2) cycle
            if (trim(keyword) == "DOMAIN OF".and.cell(i)%type /= 1) cycle
            if (location_in_list(array=region(m)%ij,element=i) == 0) call push_array(array=region(m)%ij,new_element=i)
          end do
        else
! create face region from face region
! BOUNDARY OF will pick out faces that are on the boundary
! ASSOCIATED WITH is nonsense - will just copy region
          j = ijregion
          if (trim(keyword) == "BOUNDARY OF".and.face(j)%type /= 2) cycle
          if (trim(keyword) == "DOMAIN OF".and.face(j)%type /= 1) cycle
          if (location_in_list(array=region(m)%ij,element=j) == 0) call push_array(array=region(m)%ij,new_element=j)
        end if
      end if
    end do

!---------------------
! a new region composed of a compound list of other regions

  else if (trim(keyword) == "COMPOUND") then

! note: this may mean that the region is already defined but defining twice won't hurt if the definition is the same
    if (allocated(region(m)%ij)) then
      write(*,'(a/a)') &
        "NOTE: a COMPOUND region operator is acting on region "//trim(region(m)%name)// &
        " that already contains some elements:", &
        " the previous elements will be overwritten with the new"
      deallocate(region(m)%ij)
    end if

    region_list = trim(region(m)%location(9:1000))

! loop through all regions
    do while (len_trim(region_list) /= 0) 

      region_list = adjustl(region_list(1:len(region_list))) ! remove leading spaces

      if (region_list(1:1) == "+" .or. region_list(1:1) == "-") then
        rsign = region_list(1:1)
        region_list = adjustl(region_list(2:len(region_list)))
      else
        rsign = "+"
      end if

      cut = scan(region_list,">")
      if (region_list(1:1) /= "<" .or. cut <=1 ) call &
        error_stop("format for equation region incorrect in COMPOUND operator list "//trim(region_list))

      aregion = region_list(1:cut)
      region_list = adjustl(region_list(cut+1:len(region_list)))

      nregion = region_number_from_name(name=aregion,centring=region(m)%centring,existing=existing,creatable=.false.)
! check that region exists and that centring is consistent
      if (.not.existing) call error_stop("region "//trim(aregion)//" which is part of COMPOUND region "// &
        trim(region(m)%name)//" is not found")
      if (nregion == 0) call error_stop("problem with region "//trim(aregion)//" in COMPOUND region "//trim(region(m)%name)// &
        ":- regions most likely have different centrings")

! loop through all existing ij indices in region
! - if indice exists and we are adding, ignore, otherwise do
! - if indice exists and we are subtracting, do, otherwise ignore

      if (.not.allocated(region(nregion)%ij)) call error_stop("region ij indices not allocated in COMPOUND operator for region "// &
        trim(region(nregion)%name))
        
      do nsregion = 1,ubound(region(nregion)%ij,1)  ! loop through all ij indices in constituent region
        ijregion = region(nregion)%ij(nsregion) ! NB, ijregion and region(m)%ij may actually be i or j values depending on centring
        
        nscompound = 0 ! this is the position of iregion in the compound region's ij indices
        if (allocated(region(m)%ij)) then
          do ns = 1,ubound(region(m)%ij,1)
            if (region(m)%ij(ns) == ijregion) then
              nscompound = ns ! region indice has been found in compound region's ij indices
              exit
            end if
          end do
        end if

        if (rsign == "+" .and. nscompound == 0) then
          call push_integer_array(array=region(m)%ij,new_element=ijregion) ! add region location to equation ij indices
        else if (rsign == "-" .and. nscompound /= 0) then
          if (nscompound /= ubound(region(m)%ij,1)) then ! unless it is the last element of indicies
            region(m)%ij(nscompound:ubound(region(m)%ij,1)-1) = & ! shift indicies one space to the left to remove reference
              region(m)%ij(nscompound+1:ubound(region(m)%ij,1))
          end if
          call resize_integer_array(array=region(m)%ij,change=-1) ! reduce ij array by 1
        end if

      end do

    end do

!---------------------
! this is a region that should have elements assigned to it from a mesh file

  else if (trim(keyword) == "GMSH") then

! if it doesn't then that indicates a problem
    if (.not.allocated(region(m)%ij)) write(*,'(a)') 'WARNING: the region '//trim(region(m)%name)// &
      ' which is supposed to be read in from a mesh file contains no elements: location = '//trim(region(m)%location)

!---------------------

  else if (.not.allocated(region(m)%ij)) then

! TODO: other fancier region specifications

    call error_stop('location for region '//trim(region(m)%name)//' is not understood: location = '//trim(region(m)%location))

  end if

! find ns indicies which give the data number corresponding to location i or j

  if (region(m)%centring == "cell") then
    allocate(region(m)%ns(itotal))
  else if (region(m)%centring == "face") then
    allocate(region(m)%ns(jtotal))
  else
    stop 'ERROR: a region has neither cell or face centring'
  end if

  region(m)%ns = 0 ! a zero indicates that this region does not include this i or j index
  if (allocatable_size(region(m)%ij) /= 0) then
    do ns = 1, ubound(region(m)%ij,1)
      region(m)%ns(region(m)%ij(ns)) = ns
    end do
  end if
    
  if (debug) then
    write(82,*) '# region = '//trim(region(m)%name)//': centring = '//region(m)%centring
    if (region(m)%centring == "cell") then
      write(82,*) '# i, ns'
    else
      write(82,*) '# j, ns'
    end if
    do ij = 1, ubound(region(m)%ns,1)
      write(82,*) ij, region(m)%ns(ij)
    end do
  end if

end do

!---------------------
! calculate maximum dimensions for each region
do m=1,ubound(region,1)
  region(m)%dimensions = 0
  if (allocatable_size(region(m)%ij) /= 0) then
    if (region(m)%centring == 'cell') then
      do ns = 1, ubound(region(m)%ij,1)
        region(m)%dimensions = max(region(m)%dimensions,cell(region(m)%ij(ns))%dimensions)
      end do
    else if (region(m)%centring == 'face') then
      do ns = 1, ubound(region(m)%ij,1)
        region(m)%dimensions = max(region(m)%dimensions,face(region(m)%ij(ns))%dimensions)
      end do
    end if
    maximum_dimensions = max(maximum_dimensions,region(m)%dimensions) ! set maximum number of dimensions of any region used in the simulation
  end if
end do
if (debug_sparse) write(*,'(a,i1)') 'INFO: the maximum number of dimensions of any region is ',maximum_dimensions

!---------------------
! now do reverse indicies - ie, list of regions which each cell is a member of

! cells
do i = 1, itotal
  do m = 1, ubound(region,1)
    if (allocatable_size(region(m)%ij) == 0) cycle
    if (region(m)%centring /= 'cell') cycle
    if (location_in_list(array=region(m)%ij,element=i) /= 0) call push_array(array=cell(i)%region_list,new_element=m)
  end do
end do

! faces
do j = 1, jtotal
  do m = 1, ubound(region,1)
    if (allocatable_size(region(m)%ij) == 0) cycle
    if (region(m)%centring /= 'face') cycle
    if (location_in_list(array=region(m)%ij,element=j) /= 0) call push_array(array=face(j)%region_list,new_element=m)
  end do
end do

!---------------------
! write out summary info about the regions
if (debug_sparse) write(*,'(a)') 'INFO: regions:'
do m=1,ubound(region,1)
  if (allocatable_size(region(m)%ij) == 0) then
    formatline = '(a,'//trim(dindexformat(m))//',a)'
    if (debug_sparse) write(*,fmt=formatline) ' region_number = ',m,': name = '//trim(region(m)%name)//': location = '// &
      trim(region(m)%location)//': centring = '//region(m)%centring//': contains no elements'
  else
    formatline = '(a,'//trim(dindexformat(m))//',a,'//trim(dindexformat(region(m)%dimensions))// &
      ',a,'//trim(dindexformat(region(m)%ij(1)))// &
      ',a,'//trim(dindexformat(ubound(region(m)%ij,1)))// &
      ',a,'//trim(dindexformat(region(m)%ij(ubound(region(m)%ij,1))))//')'
    if (debug_sparse) write(*,fmt=formatline) ' region_number = ',m,': name = '//trim(region(m)%name)//': location = '// &
      trim(region(m)%location)//': centring = '//region(m)%centring//': dimensions = ',region(m)%dimensions, &
      ': ij(1) = ',region(m)%ij(1),': ij(',ubound(region(m)%ij,1),') = ', &
      region(m)%ij(ubound(region(m)%ij,1))
  end if
end do
! and some warnings all the time if a region contains no elements
do m=1,ubound(region,1)
  if (allocatable_size(region(m)%ij) == 0) write(*,'(a)') 'WARNING: region '//trim(region(m)%name)//' contains no elements'
end do

!---------------------
! setup any region links
if (allocated(region_link)) then
  if (debug_sparse) write(*,'(a)') 'INFO: region_links:'
  do m = 1, ubound(region_link,1)
    call setup_region_link(m,debug_sparse)
  end do
end if

!---------------------
if (region_details_file) then
  if (debug_sparse) write(*,*) 'INFO: writing region details to region_details.txt file'

  filename = "output/region_details.txt"
  open(fdetail,file=trim(filename),status='replace',iostat=error)
  if (error /= 0) call error_stop('problem opening file '//trim(filename))

  write(fdetail,'(a)') 'MESH DETAILS:'
  formatline = '(3(a,'//trim(indexformat)//'))'
  write(fdetail,fmt=formatline) 'NODES: ktotal = ',ktotal,': kdomain = ',kdomain,': kboundary = ',kboundary
  formatline = '(a,'//trim(indexformat)//',a,a)'
  do k = 1,ktotal
    write(fdetail,fmt=formatline) 'node: k = ',k,': ',trim(print_node(node(k)))
  end do
  formatline = '(3(a,'//trim(indexformat)//'))'
  write(fdetail,fmt=formatline) 'FACES: jtotal = ',jtotal,': jdomain = ',jdomain,': jboundary = ',jboundary
  formatline = '(a,'//trim(indexformat)//',a,a)'
  do j = 1,jtotal
    write(fdetail,fmt=formatline) 'face: j = ',j,': ',trim(print_face(face(j)))
  end do
  formatline = '(3(a,'//trim(indexformat)//'))'
  write(fdetail,fmt=formatline) 'CELLS: itotal = ',itotal,': idomain = ',idomain,': iboundary = ',iboundary
  formatline = '(a,'//trim(indexformat)//',a,a)'
  do i = 1,itotal
    write(fdetail,fmt=formatline) 'cell: i = ',i,': ',trim(print_cell(cell(i)))
  end do

  close(fdetail)
end if
!---------------------

if (debug) write(*,'(a/80(1h-))') 'subroutine setup_regions'

end subroutine setup_regions

!-----------------------------------------------------------------

subroutine setup_region_link(m,debug_sparse)

! little subroutine to setup the link between a from_region and to_region
use general_module
integer :: m, nsf, nst, n_to, n_from, to_ij, from_ij, error
integer, dimension(:), allocatable :: from_ns ! temporary array of ns in from_region from ns in to_region
double precision :: maxdist, dist
double precision, dimension(totaldimensions) :: from_x, to_x ! single locations
character(len=1000) :: formatline, filename
logical :: existing, debug_sparse
logical, parameter :: debug = .false.
                  
if (debug) write(*,'(80(1h+)/a)') 'subroutine setup_region_link'

if (debug) write(*,'(a,i3)') ' finding region_link number ',m

! find to_region
region_link(m)%to_region_number = region_number_from_name(name=region_link(m)%to_region, &
  centring=region_link(m)%to_centring,existing=existing,creatable=.false.)
! check that region exists and that centring is consistent
if (.not.existing) then
  write(*,'(a)') "ERROR: "//trim(region_link(m)%to_centring)//" region "//trim(region_link(m)%to_region)// &
    " which is part of a region link function does not exist"
  stop
end if
if (region_link(m)%to_region_number == 0) then
  write(*,'(a)') "ERROR: problem with "//trim(region_link(m)%to_centring)//" region "//trim(region_link(m)%to_region)// &
    " which is part of a region link function:- region most likely has a centring which is inconsistent with its use"
  stop
end if

! find from_region
region_link(m)%from_region_number = region_number_from_name(name=region_link(m)%from_region, &
  centring=region_link(m)%from_centring,existing=existing,creatable=.false.)
! check that region exists and that centring is consistent
if (.not.existing) then
  write(*,'(a)') "ERROR: "//trim(region_link(m)%from_centring)//" region "//trim(region_link(m)%from_region)// &
    " which is part of a region link function does not exist"
  stop
end if
if (region_link(m)%from_region_number == 0) then
  write(*,'(a)') "ERROR: problem with "//trim(region_link(m)%from_centring)//" region "//trim(region_link(m)%to_region)// &
    " which is part of a region link function:- region most likely has a centring which is inconsistent with its use"
  stop
end if

! now create links - ns index to ns index
if (allocated(region_link(m)%to_ns)) deallocate(region_link(m)%to_ns)
allocate(region_link(m)%to_ns(ubound(region(region_link(m)%from_region_number)%ij,1)))
region_link(m)%to_ns = 0 ! default value if no link is found

! also create temporary storage of reverse indicies for checking purposes
allocate(from_ns(ubound(region(region_link(m)%to_region_number)%ij,1)))
from_ns = 0
n_from = 0
n_to = 0

do nsf = 1, ubound(region_link(m)%to_ns,1)
  maxdist = huge(1.d0)
  if (region_link(m)%from_centring == 'cell') then
    from_x = cell(region(region_link(m)%from_region_number)%ij(nsf))%x
  else
    from_x = face(region(region_link(m)%from_region_number)%ij(nsf))%x
  end if
  do nst = 1, ubound(region(region_link(m)%to_region_number)%ij,1)
    if (region_link(m)%to_centring == 'cell') then
      to_x = cell(region(region_link(m)%to_region_number)%ij(nst))%x
    else
      to_x = face(region(region_link(m)%to_region_number)%ij(nst))%x
    end if
    dist = distance(from_x,to_x)
    if (dist < maxdist) then
      region_link(m)%to_ns(nsf) = nst
      maxdist = dist
    end if
  end do
  if (region_link(m)%to_ns(nsf) /= 0) then
    n_from = n_from + 1 ! number of from elements that have links defined
    if (from_ns(region_link(m)%to_ns(nsf)) == 0) n_to = n_to + 1 ! number of to elements that have links defined
    from_ns(region_link(m)%to_ns(nsf)) = nsf
  end if
end do

! run some checks on linking

! check that each from_region cell has a valid link
if (n_to < ubound(region(region_link(m)%to_region_number)%ij,1)) then
  formatline = '(a,'//trim(dindexformat(ubound(region(region_link(m)%to_region_number)%ij,1)-n_to))//',a)'
  write(*,fmt=formatline) 'INFO: ',ubound(region(region_link(m)%to_region_number)%ij,1)-n_to,' elements in '// &
    trim(region_link(m)%to_centring)//' region '//trim(region_link(m)%to_region)//' do not have links from '// &
    trim(region_link(m)%from_centring)//' region '//trim(region_link(m)%from_region)
end if

! see if there is a unique to_region element for every from_region
if (n_from > n_to) then
  formatline = '(a,'//trim(dindexformat(n_from-n_to))//',a)'
  write(*,fmt=formatline) 'WARNING: there are ',n_from-n_to,' non-unique links to '//trim(region_link(m)%to_region)// &
    ' '//trim(region_link(m)%to_centring)//' elements coming from multiple '//trim(region_link(m)%from_region)// &
    ' '//trim(region_link(m)%from_centring)//' elements - i.e., there is not a ''one-to-one'''// &
    ' correspondance between the linked regions'
end if

! check that each from_region cell has a valid link
if (n_from < ubound(region_link(m)%to_ns,1)) then
  formatline = '(a,'//trim(dindexformat(ubound(region_link(m)%to_ns,1)-n_from))//',a)'
  write(*,fmt=formatline) 'WARNING: ',ubound(region_link(m)%to_ns,1)-n_from,' '//trim(region_link(m)%from_centring)// &
    ' elements in '//trim(region_link(m)%from_region)//' do not have links to '//trim(region_link(m)%to_centring)// &
    ' elements in '//trim(region_link(m)%to_region)
end if

if (debug_sparse) then
  formatline = '(a,'//trim(dindexformat(m))//',a,'//trim(dindexformat(n_from))//',a,'//trim(dindexformat(n_to))//',a)'
  write(*,fmt=formatline) " region_link ",m," from "//trim(region_link(m)%from_centring)//" region "// &
    trim(region_link(m)%from_region)//" to "//trim(region_link(m)%to_centring)//" region "// &
    trim(region_link(m)%to_region)//" contains ",n_from," forward and ",n_to," reverse links"
end if

if (link_details_file) then
  if (debug) write(*,*) 'writing region details to region_details.txt file'

  filename = "output/link_details.txt"
  open(fdetail,file=trim(filename),status='replace',iostat=error)
  if (error /= 0) call error_stop('problem opening file '//trim(filename))

  write(fdetail,'(a)') repeat('*',80)
  write(fdetail,'(a)') "region_link from "//trim(region_link(m)%from_region)//" to "//trim(region_link(m)%to_region)
  write(fdetail,'(/a)') 'FORWARD LINKS'
  do nsf = 1, ubound(region(region_link(m)%from_region_number)%ij,1)
    from_ij = region(region_link(m)%from_region_number)%ij(nsf)
    if (region_link(m)%to_ns(nsf) /= 0) then
      to_ij = region(region_link(m)%to_region_number)%ij(region_link(m)%to_ns(nsf))
    else
      to_ij = 0
    end if
    formatline = '(a,'//trim(dindexformat(nsf))//',a,'//trim(dindexformat(from_ij))//',a,'//trim(dindexformat(to_ij))//')'
    write(fdetail,fmt=formatline) '  region link ',nsf,' from ij = ',from_ij,' to ij = ',to_ij
  end do
  write(fdetail,'(/a)') 'REVERSE LINKS'
  do nst = 1, ubound(region(region_link(m)%to_region_number)%ij,1)
    to_ij = region(region_link(m)%to_region_number)%ij(nst)
    if (from_ns(nst) /= 0) then
      from_ij = region(region_link(m)%from_region_number)%ij(from_ns(nst))
    else
      from_ij = 0
    end if
    formatline = '(a,'//trim(dindexformat(nst))//',a,'//trim(dindexformat(to_ij))//',a,'//trim(dindexformat(from_ij))//')'
    write(fdetail,fmt=formatline) '  region link ',nst,' to ij = ',to_ij,' from ij = ',from_ij
  end do
  write(fdetail,'(/a)') repeat('*',80)

  close(fdetail)
end if

deallocate(from_ns)

if (debug) write(*,'(a/80(1h-))') 'subroutine setup_region_link'

end subroutine setup_region_link

!-----------------------------------------------------------------

subroutine read_constants

! here we read in values for the constants given in constants.in
! options are not handled here, but rather by setup_equations.pl

use general_module
integer :: error, cut, n, nn, ns, m, ij
character(len=1000) :: textline, keyword, name, formatline, region_name
character(len=4) :: centring
integer, dimension(:), allocatable :: region_list ! contains the last read in region_list by number
integer, dimension(:), allocatable :: region_number ! array of regions to use when reading in region_constants
double precision, dimension(:), allocatable :: region_value
logical :: existing
logical, parameter :: debug = .false.!!!!.true.
                  
if (debug) write(*,'(80(1h+)/a)') 'subroutine read_constants'

write(*,'(a)') "INFO: reading constant values from file "//trim(constants_file)
open(unit=fconstants,file=trim(constants_file),status='old',iostat=error)
if (error /= 0) call error_stop('problem opening constants file '//trim(constants_file))

fileloop: do
  read(fconstants,'(a)',iostat=error) textline
  if (error /= 0) exit fileloop ! reached end of file
  call remove_comments(textline) ! remove comments from line
  cut=scan(textline,' ') ! split input line at first space
  if (cut<=1) cycle fileloop ! if there is space at first character or string is empty, cycle
  keyword=textline(1:cut-1)
  textline=adjustl(textline(cut:len(textline)))
  if (keyword(1:5) == "CELL_".or.keyword(1:5) == "FACE_".or.keyword(1:5) == "NONE_") then
    if (keyword(1:5) == "CELL_") then
      centring = "cell"
    else if (keyword(1:5) == "FACE_") then
      centring = "face"
    else if (keyword(1:5) == "NONE_") then
      centring = "none"
    end if
    keyword=adjustl(keyword(6:len(keyword)))
  else
    centring = "none"
  end if
! if (debug) write(*,'(a)') 'INFO: in constants file found keyword = '//trim(keyword)//': with centring = '//centring
! if (debug) write(*,'(a)') 'textline = '//trim(textline)

  if (keyword == "REGION_CONSTANT" .and. centring == "none") then
    write(*,*) 'ERROR: problem in constants file: a REGION_CONSTANT must have either CELL or FACE centring'
    stop
  end if

!---------------
  if (trim(keyword) == 'REGION_LIST') then ! read in region names and create a local array (region_list) containing them
    if (allocated(region_list)) deallocate(region_list)
    do while (len_trim(textline)>0)
      cut=scan(textline,'>') ! split input line at region end delimiter
      if (textline(1:1) /= '<'.or.cut<2) then
        write(*,*) 'ERROR: a region name in '//trim(keyword)//' in constants file incorrectly specified on line:'
        write(*,*) trim(textline)
        stop
      end if
      region_name = textline(1:cut)
      textline=adjustl(textline(cut+1:len(textline)))
!write(*,*) 'region_name = '//trim(region_name)//': textline = '//trim(textline)
      m = region_number_from_name(name=region_name,centring=centring,creatable=.false.,existing=existing)
      if (.not.existing .or. m == 0) then
        write(*,*) 'ERROR: the region '//trim(region_name)//' which appears in a REGION_LIST does not exist or has ', &
          'inconsistent centring'
        stop
      end if
      call push_integer_array(array=region_list,new_element=m)
    end do
    if (debug) then
      write(*,'(a,i2,a)') 'INFO: a REGION_LIST has been found that contains ',ubound(region_list,1),' elements'
      do m=1,ubound(region_list,1)
        write(*,*) ' m = ',m,': region_list = ',region_list(m),': region_name = ', &
          trim(region(region_list(m))%name)
      end do
    end if
  end if
!---------------

!---------------
! name
  if (trim(keyword) == 'REGION_CONSTANT'.or.trim(keyword) == 'CONSTANT') then ! read in variable names
    cut=scan(textline,'>') ! find end of name string
    if (textline(1:1) /= '<'.or.cut<2) then
      write(*,*) 'ERROR: name in '//trim(keyword)//' in constants file incorrectly specified on line:'
      write(*,*) trim(textline)
      stop
    end if
    name = textline(1:cut) ! name includes <>
    textline=adjustl(textline(cut+1:len(textline)))

!---------------
! units

    cut=scan(textline,']') ! find end of units string
    if (textline(1:1) /= '['.or.cut<2) then
      write(*,*) 'WARNING: no units assigned to '//trim(keyword)//' '//trim(name)//' in constants file on line:'
      write(*,*) trim(textline)
    else
      textline=adjustl(textline(cut+1:len(textline)))
    end if

!---------------
! find constant var number

    m = var_number_from_name(name)
    if (m == 0.or.trim(var(m)%type) /= "constant") then
      write(*,*) 'ERROR: constant variable '//trim(name)//' was not found when reading constants file'
      stop
    end if
    if (trim(var(m)%centring) /= trim(centring)) then
      write(*,*) 'ERROR: constant variable '//trim(name)//' seems inconsistent when reading constants file'
      stop
    end if

!---------------

    if (allocated(region_value)) deallocate(region_value)
    if (allocated(region_number)) deallocate(region_number)

    if (trim(keyword) == 'REGION_CONSTANT') then
      if (.not.allocated(region_list)) then
        write(*,*) 'ERROR: a REGION_LIST statement must precede a REGION_CONSTANT statement in the constants file'
        stop
      end if
      allocate(region_number(ubound(region_list,1)))
      region_number = region_list
      allocate(region_value(ubound(region_list,1)))
    else if (trim(keyword) == 'CONSTANT') then
      allocate(region_value(1))
      allocate(region_number(1))
    end if

    read(textline,*,iostat=error) region_value
    if (error /= 0) call error_stop('problem reading in constant '//trim(name))
    region_value = region_value*var(m)%multiplier ! convert to specified units

!---------------
! location statements for constant

    if (trim(keyword) == 'CONSTANT' .and. trim(centring) /= "none") then
      cut=scanstring(textline,'ON ') ! find location sequence
      if (cut > 0) textline=adjustl(trim(textline(cut+2:len(textline))))
      cut=scan(textline,'>') ! find end of region name
      if (textline(1:1) == '<'.and.cut>1) then
        region_name = textline(1:cut)
      else if (centring == "cell") then
        region_name = "<all cells>"
      else if (centring == "face") then
        region_name = "<all faces>"
      end if
      region_number(1) = region_number_from_name(name=region_name,creatable=.false.,existing=existing)
      if (.not.existing .or. region_number(1) == 0) then
        write(*,*) 'ERROR: the region '//trim(region_name)//' which appears in a CONSTANT statement does not exist'
        stop
      end if
    end if

!---------------
! finally assign values to correct regions

    if (trim(centring) == "none") then
      var(m)%funk(1)%v = region_value(1)
    else
      do n = 1, ubound(region_number,1)
        if (trim(region(region_number(n))%centring) /= trim(var(m)%centring)) then
          write(*,*) 'ERROR: when reading the constants file, the centring of constant '//trim(var(m)%name)// &
            ' does not match that of region '//trim(region(region_number(n))%name)
          stop
        end if
        do nn = 1, ubound(region(region_number(n))%ij,1)
          ij = region(region_number(n))%ij(nn)
          ns = nsvar(m=m,ij=ij)
          if (ns /= 0) var(m)%funk(ns)%v = region_value(n)
        end do
      end do
    end if

    if (debug) then
      formatline = '(a,'//trim(floatformat)//')'
      if (centring == "none") then
        write(*,fmt=formatline) ' read in '//var(m)%centring//' '//trim(var(m)%type)//' '//trim(var(m)%name)// &
          ': value = ',var(m)%funk(1)%v
      else
        write(*,'(a)') ' read in '//var(m)%centring//' '//trim(var(m)%type)//' '//trim(var(m)%name)//': values per region:'
        do n = 1, ubound(region_number,1)
          write(*,fmt=formatline) '  '//trim(region(region_number(n))%name)//': value = ',region_value(n)
        end do
      end if
    end if
    
  end if
    
!---------------

end do fileloop

close(fconstants)

if (debug) write(*,'(a/80(1h-))') 'subroutine read_constants'

end subroutine read_constants

!----------------------------------------------------------------------------

subroutine setup_vars

! here we allocate array elements for the fields and functions and initialise fields

use general_module
use equations_module
integer :: m, ns, n, mc, pptotal, mtype, o
character(len=1000) :: formatline, component_list
logical :: existing, first
logical, parameter :: debug = .false.
logical :: debug_sparse = .false.!!!!.true. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 30/10/2013
integer xs_reading
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



if (debug) debug_sparse = .true.
                  
if (debug_sparse) write(*,'(80(1h+)/a)') 'subroutine setup_vars'

! allocate var funks, zero these funks and set region numbers
do m = 1, ubound(var,1)
  if (var(m)%centring /= "none") then
    var(m)%region_number = region_number_from_name(name=var(m)%region,centring=var(m)%centring,existing=existing)
! check that region exists and that the variable and region centring are consistent
    if (.not.existing) call error_stop('there is a problem with region '//trim(var(m)%region)//' which is associated with '// &
      trim(var(m)%type)//' '//trim(var(m)%name)//': the region does not exist')
    if (var(m)%region_number == 0) call error_stop('there is a problem with region '//trim(var(m)%region)// &
      ' which is associated with '//trim(var(m)%type)//' '//trim(var(m)%name)//': most probably the centrings are inconsistent')
    if (var(m)%type /= 'local') allocate(var(m)%funk(ubound(region(var(m)%region_number)%ij,1)))
  else
    var(m)%region_number = 0 ! dummy region number for none centred variables
    if (var(m)%type /= 'local') allocate(var(m)%funk(1))
  end if
  if (allocated(var(m)%funk)) then
    do ns = 1, ubound(var(m)%funk,1)
      call reset_funk(var(m)%funk(ns))
    end do
  end if
end do

! copy region numbers set in components to compound
do mc = 1, ubound(compound,1)
  first = .true.
  do n = 1, ubound(compound(mc)%component,1)
    m = compound(mc)%component(n)
    if (m /= 0) then
      if (first) then
        compound(mc)%region_number = var(m)%region_number
        first = .false.
      else
! also check whether the remaining valid components are from the same region
        if (compound(mc)%region_number /= var(m)%region_number) &
          call error_stop("components of compound "//trim(compound(mc)%name)//" are from inconsistent regions")
      end if
    end if
  end do
end do

! now setup var_lists
allocate(var_list(var_list_number(type="all",centring="all")))
do mtype = 1, ubound(var_types,1)
  do m = 1, ubound(var,1)
    if (trim(var(m)%type) /= trim(var_types(mtype))) cycle ! lists are created in order
    call push_array(array=var_list(var_list_number(type="all",centring="all"))%list,new_element=m)
    call push_array(array=var_list(var_list_number(type=var(m)%type,centring="all"))%list,new_element=m)
    call push_array(array=var_list(var_list_number(type="all",centring=var(m)%centring))%list,new_element=m)
    call push_array(array=var_list(var_list_number(type=var(m)%type,centring=var(m)%centring))%list,new_element=m)
  end do
end do
if (debug) then
  write(82,*) 'var_list details:'
  do m = 1, var_list_number(type="all",centring="all")
    write(82,*) 'var_list number = ',m,': allocatable_size = ',allocatable_size(var_list(m)%list)
    if (allocatable_size(var_list(m)%list) > 0) then
      do n = 1, allocatable_size(var_list(m)%list)
        write(82,*) '  contains '//trim(var(var_list(m)%list(n))%name)//' of centring = '// &
          trim(var(var_list(m)%list(n))%centring)//' and type = '//trim(var(var_list(m)%list(n))%type)
      end do
    end if
  end do
end if

! run through unknowns setting derivatives and calculating ptotal
ptotal = 0 ! this is the number of unknown variables (should equal pptotal)
do n = 1, allocatable_size(var_list(var_list_number(centring="all",type="unknown"))%list)
  m = var_list(var_list_number(centring="all",type="unknown"))%list(n)
  do ns = 1, ubound(var(m)%funk,1)
    ptotal = ptotal + 1
    call push_array(var(m)%funk(ns)%pp,ptotal)
    call push_array(var(m)%funk(ns)%dv,1.d0)
    var(m)%funk(ns)%ndv = 1
  end do
end do
  
! run through equations calculating pptotal
pptotal = 0 ! this is the number of equations (should equal ptotal)
do n = 1, allocatable_size(var_list(var_list_number(centring="all",type="equation"))%list)
  m = var_list(var_list_number(centring="all",type="equation"))%list(n)
  pptotal = pptotal + ubound(var(m)%funk,1)
end do

! count maximum relstep in the transients
relstepmax = 0
if (transient_simulation) then
  do n = 1, allocatable_size(var_list(var_list_number(centring="all",type="transient"))%list)
    m = var_list(var_list_number(centring="all",type="transient"))%list(n)
    relstepmax = max(relstepmax,var(m)%relstep)
  end do
end if

! add default output options - these come after user-set options so only take effect if user has not already set
! COMPONENTS
! output is off by default
! stepoutput is off by default
! input is off by default (values will be read in via compounds for unknowns by default)
! elementdata used by default
do n = 1, ubound(var,1)
  call push_character_array(array=var(n)%options,new_element='nooutput')
  call push_character_array(array=var(n)%options,new_element='nostepoutput')
  call push_character_array(array=var(n)%options,new_element='noinput')
  call push_character_array(array=var(n)%options,new_element='elementdata')
end do
! COMPOUNDS
! output - if var is an unknown or an output or a cell centred derived then on by default
!   for a transient simulation the last relstep data does not need to be output for restart
! stepoutput - on by default for none centred unknown, derived, transients and outputs for current relstep
!   note output variables that should only be updated when the msh files are written should have
!   stepoutputnoupdate or nosetupoutput set as options
! input is on for all unknowns
! elementdata used by default unless cell centred scalar, then elementnodedata
! elementnodedata cannot be used for face centred quantities (or none centred for that matter)
do n = 1, ubound(compound,1)
  if ((compound(n)%type == 'unknown'.or.compound(n)%type == 'output'.or. &
    (compound(n)%type == 'derived'.and.compound(n)%centring == 'cell').or. &
    compound(n)%type == 'transient').and.compound(n)%relstep < max(relstepmax,1)) then
    call push_character_array(array=compound(n)%options,new_element='output')
  else
    call push_character_array(array=compound(n)%options,new_element='nooutput')
  end if
  if ((compound(n)%type == 'unknown'.or.compound(n)%type == 'derived'.or. &
    compound(n)%type == 'transient'.or.compound(n)%type == 'output') &
    .and.compound(n)%centring == 'none'.and.compound(n)%relstep == 0) then
    call push_character_array(array=compound(n)%options,new_element='stepoutput')
  else
    call push_character_array(array=compound(n)%options,new_element='nostepoutput')
  end if
  if (compound(n)%type == 'unknown'.or.compound(n)%type == 'transient'.or. &
      compound(n)%type == 'constant') then
    call push_character_array(array=compound(n)%options,new_element='input')
  else
    call push_character_array(array=compound(n)%options,new_element='noinput')
  end if
  if (compound(n)%centring /= 'none') then
    if (compound(n)%rank == 'scalar'.and.compound(n)%centring == 'cell') then
      call push_character_array(array=compound(n)%options,new_element='elementnodedata')
    else
      call push_character_array(array=compound(n)%options,new_element='elementdata')
    end if
  end if
end do

! run through constants reading in file values and setting any equation determined ones

! ABERNAL: Instructions to choose betwen NEMTAB or CONSTANTS.IN to read the cross sections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 30/10/2013
!call read_constants
     open(unit=95,file='input',status='old')
     do n=1,11 ! Line 12 contains the information about cross section reading
        read(95,*)
     enddo
     read(95,*) xs_reading
     close (95)

     select case (xs_reading)
        case(1) ! file "constants.in"
           call read_constants
        case(2,3,4,5,6,7) ! NEMTAB files
           call xsnemtab
        case default
           print '(///)'
           print *,'ERROR READING THE LINE 12 OF FILE "input"'
           print '(///)'
           stop
     end select
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! ABERNAL: He visto que esta rutina no afecta al calculo y se consigue reducir
! tiempo de calculo si no se llama. Si hay algun problema con los resultados,
! hay que revisarla
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 14/06/2013
!call update_constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





! spray out some info about all the var variables
if (debug_sparse) then
  if (debug) then
    write(*,'(a)') '---------------------------------------------------------------'
    write(*,'(a)') 'INFO: component variables:'
    do mtype = 1, ubound(var_types,1)
      write(*,'(a/a)') '------------------------------------',trim(var_types(mtype))
      do n = 1, allocatable_size(var_list(var_list_number(centring="all",type=trim(var_types(mtype))))%list)
        m = var_list(var_list_number(centring="all",type=trim(var_types(mtype))))%list(n)
        if (var(m)%centring == 'none') then
          formatline = '(a,'// &
            trim(dindexformat(m))//',a,'// &
            trim(dindexformat(var(m)%relstep))//',a,i1,a,'// &
            trim(dindexformat(var(m)%someloop))//',a'// &
            repeat(',a',ubound(var(m)%options,1))//')'
          write(*,fmt=formatline) ' variable_number = ',m, &
            ': name = '//trim(var(m)%name)// &
            ': type = '//trim(var(m)%type)// &
            ': centring = '//var(m)%centring// &
            ': rank = '//trim(var(m)%rank)// &
            ': relstep = ',var(m)%relstep, &
            ': deriv = ',var(m)%deriv, &
            ': someloop = ',var(m)%someloop, &
            ': (prioritised) options =',(' '//trim(var(m)%options(o)),o=1,ubound(var(m)%options,1))
        else
          formatline = '(a,'// &
            trim(dindexformat(m))//',a,'// &
            trim(dindexformat(var(m)%region_number))//',a,'// &
            trim(dindexformat(var(m)%relstep))//',a,i1,a,'// &
            trim(dindexformat(var(m)%someloop))//',a,'// &
            trim(dindexformat(region(var(m)%region_number)%ij(1)))//',a,'// &
            trim(dindexformat(ubound(region(var(m)%region_number)%ij,1)))//',a,'// &
            trim(dindexformat(region(var(m)%region_number)%ij(ubound(region(var(m)%region_number)%ij,1))))// &
            ',a'//repeat(',a',ubound(var(m)%options,1))//')'
          write(*,fmt=formatline) ' variable_number = ',m, &
            ': name = '//trim(var(m)%name)// &
            ': type = '//trim(var(m)%type)// &
            ': centring = '//var(m)%centring// &
            ': rank = '//trim(var(m)%rank)// &
            ': relstep = ',var(m)%relstep, &
            ': deriv = ',var(m)%deriv, &
            ': someloop = ',var(m)%someloop, &
            ': region = '//trim(var(m)%region)// &
            ': region_number = ',var(m)%region_number, &
            ': ij(1) = ',region(var(m)%region_number)%ij(1), &
            ': ij(',ubound(region(var(m)%region_number)%ij,1),') = ', &
            region(var(m)%region_number)%ij(ubound(region(var(m)%region_number)%ij,1)), &
            ': (prioritised) options =',(' '//trim(var(m)%options(o)),o=1,ubound(var(m)%options,1))
        end if
      end do
    end do
  end if

  write(*,'(a)') '---------------------------------------------------------------'
  write(*,'(a)') 'INFO: compound variables:'
  do mtype = 1, ubound(var_types,1)
    write(*,'(a/a)') '------------------------------------',trim(var_types(mtype))
    do m = 1, ubound(compound,1)
      if (trim(compound(m)%type) /= trim(var_types(mtype))) cycle
      component_list = "["
      do o = 1, ubound(compound(m)%component,1)
        if (compound(m)%component(o) == 0) then
          component_list = trim(component_list)//' (empty)'
        else
          component_list = trim(component_list)//' '//trim(var(compound(m)%component(o))%name)
        end if
      end do
      component_list = trim(component_list)//' ]'
      if (compound(m)%centring == 'none') then
        formatline = '(a,'// &
          trim(dindexformat(m))//',a,'// &
          trim(dindexformat(compound(m)%relstep))//',a,i1,a'// &
          repeat(',a',ubound(compound(m)%options,1))//')'
        write(*,fmt=formatline) ' compound_number = ',m, &
          ': name = '//trim(compound(m)%name)// &
          ': type = '//trim(compound(m)%type)// &
          ': centring = '//compound(m)%centring// &
          ': rank = '//trim(compound(m)%rank)// &
          ': relstep = ',compound(m)%relstep, &
          ': deriv = ',compound(m)%deriv, &
          ': component_list ='//trim(component_list)// &
          ': (prioritised) options =',(' '//trim(compound(m)%options(o)),o=1,ubound(compound(m)%options,1))
      else
        formatline = '(a,'// &
          trim(dindexformat(m))//',a,'// &
          trim(dindexformat(compound(m)%relstep))//',a,i1,a,'// &
          trim(dindexformat(compound(m)%region_number))// &
          ',a'//repeat(',a',ubound(compound(m)%options,1))//')'
        write(*,fmt=formatline) ' compound_number = ',m, &
          ': name = '//trim(compound(m)%name)// &
          ': type = '//trim(compound(m)%type)// &
          ': centring = '//compound(m)%centring// &
          ': rank = '//trim(compound(m)%rank)// &
          ': relstep = ',compound(m)%relstep, &
          ': deriv = ',compound(m)%deriv, &
          ': region = '//trim(compound(m)%region)// &
          ': region_number = ',compound(m)%region_number, &
          ': component_list ='//trim(component_list)// &
          ': (prioritised) options =',(' '//trim(compound(m)%options(o)),o=1,ubound(compound(m)%options,1))
      end if
    end do
  end do
  write(*,'(a)') '---------------------------------------------------------------'
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 10/10/2013
!if (pptotal /= ptotal) then
!  write(*,*) 'ERROR: the total number of equations does not match the total number of unknown variables'
!  write(*,*) '  pptotal (number of equations) = ',pptotal
!  write(*,*) '  ptotal (number of unknown variables) = ',ptotal
!  stop
!else if (pptotal == 0) then
!  write(*,*) 'WARNING: no equations/unknowns are being solved for'
!end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! allocate equation_magnitude which is used in residual calculations
allocate(equation_magnitude(ptotal))
equation_magnitude = 0.d0

if (debug_sparse) write(*,'(a/80(1h-))') 'subroutine setup_vars'

end subroutine setup_vars

!-----------------------------------------------------------------

end module setup_module

!-----------------------------------------------------------------
