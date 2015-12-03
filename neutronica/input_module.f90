! file src/input_module.f90
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
module input_module

! statements specifying different data types and parameters
implicit none

! these are the only subroutines that are accessible from outside this module
private
public read_gmesh

!-----------------------------------------------------------------
contains

!-----------------------------------------------------------------

subroutine read_gmesh(contents,var_number)

! should only be called for constants, unknowns and initial transients so should never have local variables

use general_module
use gmsh_module
character(len=4) :: contents ! either mesh or data
integer, optional :: var_number ! if reading data then this is the var_number to read
integer :: gmesh_number
character(len=1000) :: formatline
character(len=100) :: input_option
logical :: success
logical, parameter :: debug = .false.

if (debug) write(*,'(80(1h+)/a)') 'subroutine read_gmesh'

if (.not.allocated(gmesh)) stop 'ERROR: no gmsh meshes have been specified'

success = .false.

do gmesh_number = 0, ubound(gmesh,1)
  formatline = '(a,'//trim(dindexformat(gmesh_number))//',a)'
  input_option = check_option(gmesh(gmesh_number)%options,input_gmesh_options) ! get input option for gmesh
  if (trim(input_option) == "noinput" .or. &
    (contents == 'data'.and.(trim(input_option) == "centringmeshinput".or.trim(input_option) == "meshinput"))) then
    if (gmesh_number /= 0.or.debug) &
      write(*,fmt=formatline) 'INFO: skipping '//trim(contents)//' read for gmesh file (',gmesh_number,') having basename '// &
      trim(gmesh(gmesh_number)%basename)//' based on gmesh input option '//trim(input_option)
  else if (contents == 'mesh') then
    write(*,fmt=formatline) 'INFO: performing mesh read for gmesh file (',gmesh_number,') having basename '// &
      trim(gmesh(gmesh_number)%basename)//' based on gmesh input option '//trim(input_option)
    success = .true.
    if (trim(input_option) == "centringinput" .or. trim(input_option) == "centringmeshinput") then
      write(*,'(a)') 'INFO: reading in cell element mesh:'
      call read_gmesh_file(gmesh_number=gmesh_number,file_centring='cell',contents=contents) ! read in cell mesh for gmesh_number
      write(*,'(a)') 'INFO: reading in face element mesh:'
      call read_gmesh_file(gmesh_number=gmesh_number,file_centring='face',contents=contents) ! read in face mesh for gmesh_number
    else
      call read_gmesh_file(gmesh_number=gmesh_number,file_centring='both',contents=contents) ! read in combined mesh for gmesh_number
    end if
  else if (contents == 'data') then
    if (.not.present(var_number)) call error_stop("read_gmesh called with data contents but no var_number")
    if (trim(check_option(var(var_number)%options,input_options)) == "noinput".and. &
      trim(check_option(compound(var(var_number)%compound_number)%options,input_options)) == "noinput") then
      if (debug) write(*,fmt=formatline) 'INFO: skipping data read for variable '//trim(var(var_number)%name)// &
        ' for gmesh file (',gmesh_number,') having basename '//trim(gmesh(gmesh_number)%basename)//' based on component option '// &
        trim(check_option(var(var_number)%options,input_options))//' and compound option '// &
        trim(check_option(compound(var(var_number)%compound_number)%options,input_options)) 
    else
      if (debug) write(*,fmt=formatline) 'INFO: initiating data read for variable '//trim(var(var_number)%name)// &
        ' for gmesh file (',gmesh_number,') having basename '//trim(gmesh(gmesh_number)%basename)// &
        ' based on gmesh input option '//trim(input_option)//', component option '// &
        trim(check_option(var(var_number)%options,input_options))//' and compound option '// &
        trim(check_option(compound(var(var_number)%compound_number)%options,input_options)) 
      if (trim(input_option) == "centringinput") then
        if (debug) write(*,'(a)') 'INFO: reading in '//var(var_number)%centring//' element data:'
        call read_gmesh_file(gmesh_number=gmesh_number,file_centring=var(var_number)%centring,contents=contents, &
          var_number=var_number) ! read in centred data for gmesh_number
      else
        call read_gmesh_file(gmesh_number=gmesh_number,file_centring='both',contents=contents,var_number=var_number) ! read in combined data for gmesh_number
      end if
    end if
  end if
end do

if (contents == 'mesh'.and..not.success) write(*,'(a)') "WARNING: no meshes have been read in"

if (debug) write(*,'(a/80(1h-))') 'subroutine read_gmesh'

end subroutine read_gmesh

!-----------------------------------------------------------------

subroutine read_gmesh_file(gmesh_number,file_centring,contents,var_number)

! this subroutine reads in a gmsh file, for contents = mesh or data
! specifically for contents=mesh it:
! - defines any regions not already defined (name, centring and location)
! - defines nodes (x, jface) and sets ktotal
! - defines cells (dimensions, jface, knode, gtype), adds to relevant region list and sets itotal
! - defines faces (dimensions, icell, knode, gtype), adds to relevant region list and sets jtotal

use general_module
use gmsh_module
integer :: gmesh_number, gmsh_bin, n, error
integer, optional :: var_number
real :: gmsh_version
character(len=1000) :: textline, filename
character(len=4) :: file_centring, contents
logical, parameter :: debug = .false. ! this is passed to called reading subroutines too

if (debug) write(*,'(80(1h+)/a)') 'subroutine read_gmesh_file'

filename = trim(gmesh(gmesh_number)%filename)
! for centringinput alter filename to suit either cell or face
if (file_centring == "cell" .or. file_centring == "face" .or. file_centring == "none") then
  n = scan(filename,'.',.true.)
  if (n == 0) call error_stop("problem constructing centring filename in subroutine read_gmesh_file: "//trim(filename))
  filename = filename(1:n)//trim(file_centring)//"."//filename(n+1:len(filename))
end if
if (debug) write(*,*) ' reading in filename = '//trim(filename)
open(unit=fgmsh,file=trim(filename),status='old',iostat=error)
if (error /= 0) call error_stop('problem opening mesh file '//trim(filename))

!---------------------------------------------
! check version and format of file
do 
  read(fgmsh,'(a)',iostat=error) textline
  if (error /= 0) call error_stop('gmsh file version not found for '//trim(filename))
  if (trim(textline) == "$MeshFormat") exit
end do
read(fgmsh,*,iostat=error) gmsh_version, gmsh_bin
if (error /= 0) call error_stop('problem reading in gmsh file version in '//trim(filename))
if (gmsh_version < 2.1) then
  write(*,'(a,f4.1,a)') 'ERROR: the gmsh file version (',gmsh_version,') is too old and is not supported (update to v2.1) in '// &
    trim(filename)
  stop
else if (gmsh_version > 2.2) then
  if (contents == 'mesh') write(*,'(a,f4.1,a)') 'WARNING: the gmsh mesh file version (',gmsh_version, &
    ') is not fully suported yet (but may work) in '//trim(filename)
end if
if (gmsh_bin /= 0) call error_stop('gmsh binary file type not supported in '//trim(filename))
if (debug) write(*,*) 'gmsh_version = ',gmsh_version
      
if (contents == 'mesh') then
!---------------------------------------------
! read in physical names and possibly create any required regions
! if reading the face file check that each physical region is defined

  if (file_centring /= "face") then
    call read_gmesh_regions(check=.false.,gmesh_number=gmesh_number,filename=filename,debug=debug) ! read in
    if (debug) write(*,*) 'gregions checked'
  else
    call read_gmesh_regions(check=.true.,gmesh_number=gmesh_number,filename=filename,debug=debug) ! check for consistency
    if (debug) write(*,*) 'gregions read in'
  end if

!---------------------------------------------
! read in nodes
! if reading the face file check that nodes are exactly the same as for the cell file

  if (file_centring /= "face") then
    call read_gmesh_nodes(check=.false.,gmesh_number=gmesh_number,filename=filename,debug=debug) ! read in
    if (debug) write(*,*) 'nodes checked'
  else
    call read_gmesh_nodes(check=.true.,gmesh_number=gmesh_number,filename=filename,debug=debug) ! check for consistency
    if (debug) write(*,*) 'nodes read in'
  end if

!---------------------------------------------
! read in elements, creating cells and faces as we go, if nodes have been read in previously
  if (allocated(gmesh(gmesh_number)%knode_from_gnode)) then
    call read_gmesh_elements(gmesh_number=gmesh_number,filename=filename,debug=debug)
    if (debug) write(*,*) 'elements read in'
  end if

!---------------------------------------------
else if (contents == 'data') then
! read in data
  if (.not.present(var_number)) call error_stop("read_gmesh_file called with data contents but no var_number")
  call read_gmesh_data(gmesh_number=gmesh_number,filename=filename,var_number=var_number,debug=debug)
  if (debug) write(*,*) 'data read in'

end if
!---------------------------------------------

close(fgmsh)

if (debug) write(*,'(a/80(1h-))') 'subroutine read_gmesh_file'

end subroutine read_gmesh_file

!-----------------------------------------------------------------

subroutine read_gmesh_data(gmesh_number,filename,var_number,debug)

use general_module
use gmsh_module
integer :: error, gmesh_number, ntags, var_number, compound_number, nrank, nelements, gelement, nelement, n, k, ij, ns, &
  region_number, nnodes
character(len=1000) :: textline, filename, formatline, name
character(len=100) :: data_type
character(len=4) :: centring
double precision, dimension(:,:), allocatable :: data_array
logical :: compoundl
logical :: debug

if (debug) write(*,'(a)') 'INFO: reading in data for variable '//trim(var(var_number)%name)//' from file '//trim(filename)//':'

rewind(fgmsh)

data_type = ''

main_loop: do 
  read(fgmsh,'(a)',iostat=error) textline
  if (error /= 0) exit
  if (data_type /= '') then
    if (trim(textline) == '$End'//trim(data_type)) data_type = ''
  else if (trim(textline) == "$Data" .or. trim(textline) == "$ElementData" .or. trim(textline) == "$ElementNodeData") then
    data_type = trim(textline(2:len(textline)))

! if no elements are defined then skip the data (really this is a file inconsistency though)
    if (.not.allocated(gmesh(gmesh_number)%gelement).and.trim(data_type) /= 'Data') then
      write(*,'(a)') 'WARNING: skipping orphaned '//trim(data_type)//' data in gmsh file '//trim(filename)// &
        ' for which there are no elements defined'
      cycle main_loop
    end if

! character data: ie, variable name
    read(fgmsh,*,iostat=error) ntags ! ntags is the number of character variables, the first of which should be the variable name
    if (error /= 0.or.ntags < 1) then
      write(*,*) 'WARNING: problem reading file for data type '//trim(data_type)//' in gmsh file '//trim(filename)
      cycle main_loop
    end if
    do n = 1, ntags
      read(fgmsh,*,iostat=error) textline
      if (error /= 0) then
        write(*,*) 'WARNING: problem reading file for data type '//trim(data_type)//' in gmsh file '//trim(filename)
        cycle main_loop
      end if
      if (n == 1) then
        name = var(var_number)%name
        centring = var(var_number)%centring
        region_number = var(var_number)%region_number
        if (trim(textline) == trim(var(var_number)%name)) then
          compoundl = .false.
          compound_number = 0
        else if (trim(textline) == trim(compound(var(var_number)%compound_number)%name)) then
          compoundl = .true.
          compound_number = var(var_number)%compound_number
        else
          cycle main_loop ! this data field does not contain the requested var so move on
        end if
      else
        write(*,*) 'WARNING: ignoring character string in '//trim(filename)//' in gmsh file'//trim(textline)
      end if
    end do
    if (debug) then
      if (compoundl) then
        write(*,*) 'found compound variable '//trim(name)//' for var '//trim(var(var_number)%name)//' in gmsh file '//trim(filename)
      else
        write(*,*) 'found var variable '//trim(name)//' in gmsh file '//trim(filename)
      end if
    end if

! perform some sanity checks on data_type
    if ((data_type == 'Data'.and.centring /= 'none').or.(data_type /= 'Data'.and.centring == 'none').or. &
      (data_type == 'ElementNodeData'.and.centring /= 'cell')) then
      write(*,*) 'WARNING: inconsistent variable centring for variable '//trim(name)//' in gmsh file '//trim(filename)
      cycle main_loop
    end if
  
! real data: empty
    read(fgmsh,*,iostat=error) ntags ! ntags is the number of real variables, should be 0
    if (error /= 0) then
      write(*,*) 'WARNING: problem reading file for variable '//trim(name)//' in gmsh file '//trim(filename)
      cycle main_loop
    end if
    do n = 1, ntags
      read(fgmsh,*,iostat=error) textline
    end do

! integer data: timestep, number of entries per elements, number of elements
    read(fgmsh,*,iostat=error) ntags ! ntags is the number of real variables, should be 3
    if (error /= 0.or.ntags < 3) then
      write(*,*) 'WARNING: problem reading file for variable '//trim(name)//' in gmsh file '//trim(filename)
      cycle main_loop
    end if
    read(fgmsh,*,iostat=error) timestep
    if (error /= 0) then
      write(*,*) 'WARNING: problem reading file for variable '//trim(name)//' in gmsh file '//trim(filename)
      cycle main_loop
    end if
    read(fgmsh,*,iostat=error) nrank
    if (error /= 0) then
      write(*,*) 'WARNING: problem reading file for variable '//trim(name)//' in gmsh file '//trim(filename)
      cycle main_loop
    end if
    read(fgmsh,*,iostat=error) nelements
    if (error /= 0) then
      write(*,*) 'WARNING: problem reading file for variable '//trim(name)//' in gmsh file '//trim(filename)
      cycle main_loop
    end if

! perform a sanity check on variable rank
    if (compoundl) then
      if (nrank /= ubound(compound(compound_number)%component,1)) then
        write(*,*) 'WARNING: rank of compound '//trim(name)//' does not match that in gmsh file '//trim(filename)
        cycle main_loop
      end if
    else
      if (nrank /= 1) then
        write(*,*) 'WARNING: rank of component variable '//trim(name)//' does not match that in gmsh file '//trim(filename)
        cycle main_loop
      end if
    end if

    if (debug) write(*,*) ' about to read in data for variable = '//trim(name)//': compoundl = ',compoundl,': var_number = ', &
      var_number,': compound_number = ',compound_number,': nrank = ',nrank,': nelements = ',nelements,': data_type = '// &
      trim(data_type)

! do a temporary allocation of data_array, assuming one element
    allocate(data_array(nrank,1))
    ns = 1

! now to read in actual data
! errors found from here onwards indicate a real problem, so terminate if any found
    do nelement = 1, nelements
      read(fgmsh,'(a)',iostat=error) textline
      if (error /= 0) call error_stop ('problem reading file for variable '//trim(name)//' in gmsh file '//trim(filename))
      if (trim(data_type) == 'ElementNodeData') then
        read(textline,*,iostat=error) gelement, ntags
      else
        read(textline,*,iostat=error) gelement
        ntags = 1
      end if
      if (error /= 0) call error_stop ('problem reading file for variable '//trim(name)//' in gmsh file '//trim(filename))
      if (gelement < 1.or.gelement > ubound(gmesh(gmesh_number)%gelement,1)) &
        call error_stop('gelement out of range for variable '//trim(name)//' in gmsh file '//trim(filename))
      if (ntags < 1) call error_stop('ntags incorrect for variable '//trim(name)//' in gmsh file '//trim(filename))

! find cell or face if it is not none centred and check that it is within variable's region
      if (centring /= 'none') then
        ij = 0
        if (centring == 'cell') then
          ij = gmesh(gmesh_number)%gelement(gelement)%icell
        else if  (centring == 'face') then
          ij = gmesh(gmesh_number)%gelement(gelement)%jface
        end if
        if (ij == 0.and.debug) write(*,*) 'ij = ',ij,': gelement = ',gelement,': textline = ',trim(textline)
        if (ij == 0) call error_stop(trim(centring)//' centred data references an invalid element for variable '//trim(name)// &
          ' in gmsh file '//trim(filename))
        ns = region(region_number)%ns(ij)
        if (ns == 0) cycle main_loop ! silently skip any data that is not within the variable's region
      end if

! find nnodes depending on data_type and for ElementNodeData reallocate data_array to match
      nnodes = 1 ! do not have to keep reallocating data_array
      if (trim(data_type) == 'ElementNodeData') then ! NB, this also implies cell centring
        nnodes = ubound(cell(ij)%knode,1)
        if (ntags /= nnodes) call error_stop('mismatch between number of nodes that element has and what file specifies for '// &
          'variable '//trim(name)//' in gmsh file '//trim(filename))
        if (ubound(data_array,2) /= nnodes) then
          deallocate(data_array)
          allocate(data_array(nrank,nnodes))
        end if
      end if

! now read in data
      if (trim(data_type) == 'ElementNodeData') then
        read(textline,*,iostat=error) gelement, ntags, ((data_array(n,k),n=1,nrank),k=1,nnodes)
      else
        read(textline,*,iostat=error) gelement, ((data_array(n,k),n=1,nrank),k=1,nnodes)
      end if
      if (error /= 0)  call error_stop('problem reading in data line for variable '//trim(name)// &
          ' in gmsh file '//trim(filename))

! finally use data to set variable value
      if (compoundl) then
        n = location_in_list(compound(compound_number)%component,var_number)
      else
        n = 1
      end if
      if (data_type == 'Data'.or.data_type == 'ElementData') then
        var(var_number)%funk(ns)%v = data_array(n,1)
      else
! have to use cell node kernel (4) to average data for elementnodedata
        var(var_number)%funk(ns)%v = 0.d0
        do k = 1, nnodes
          var(var_number)%funk(ns)%v = var(var_number)%funk(ns)%v + cell(ij)%kernel(4)%v(k)*data_array(n,k)
        end do
      end if

    end do
    deallocate(data_array)
    if (compoundl) then
      formatline = '(a,'//trim(dindexformat(nelements))//',a,'//trim(dindexformat(nrank))//',a)'
      write(*,fmt=formatline) 'INFO: read data for variable '//trim(var(var_number)%name)//' from compound variable '// &
        trim(compound(compound_number)%name)//': centring = '//trim(centring)//': data_type = '//trim(data_type)// &
        ': nelements = ',nelements,': nrank = ',nrank,': gmsh file = '//trim(filename)
    else
      formatline = '(a,'//trim(dindexformat(nelements))//',a)'
      write(*,fmt=formatline) 'INFO: read data for variable '//trim(var(var_number)%name)// &
        ' from component: centring = '//trim(centring)//': data_type = '//trim(data_type)//': nelements = ',nelements, &
        ': gmsh file = '//trim(filename)
    end if

  end if
end do main_loop

end subroutine read_gmesh_data

!-----------------------------------------------------------------

subroutine read_gmesh_elements(gmesh_number,filename,debug)

use general_module
use gmsh_module
integer :: error, gmesh_number
integer :: i, j, k, n, m, jj, nelements, nelements_tenth, gelement, gtype, ntags, region_number, nnode, nface, nfaces, &
  gelement_max, change, gelement_previous
integer, dimension(:), allocatable :: tags, local_gnodes, local_knodes
integer, dimension(:), allocatable :: icell_from_gelement, jface_from_gelement
type(cell_type) :: default_cell
type(face_type) :: default_face
character(len=4) :: centring
character(len=1000) :: textline, filename
logical :: debug

rewind(fgmsh)

allocate(default_face%icell(2)) ! allocate spaces for the two surrounding faces

do 
  read(fgmsh,'(a)',iostat=error) textline
! if (error /= 0) call error_stop('list of elements not found in gmsh file '//trim(filename))
  if (error /= 0) then
    write(*,'(a)') "WARNING: list of elements not found in gmsh file "//trim(filename)
    return
  end if
  if (trim(textline) == "$Elements") exit
end do
read(fgmsh,*,iostat=error) nelements
if (error /= 0) call error_stop('problem reading in number of elements in gmsh file '//trim(filename))
if (debug) write(*,*) 'number of elements: nelements = ',nelements

! allocate temporary gelement lookup arrays based on estimate of required size
allocate(icell_from_gelement(nelements))
allocate(jface_from_gelement(nelements))
gelement_max = 0 ! this records the highest gelement number
icell_from_gelement = 0
jface_from_gelement = 0

nelements_tenth = max(int(float(nelements)/10.),1)
if (mod(nelements,nelements_tenth) /= 0) nelements_tenth = nelements_tenth + 1

do n = 1, nelements

  if (mod(n,nelements_tenth) == 0 .or. n == nelements) write(*,'(a,i3,a)') &
    'INFO: reading in mesh:',min(int(float(n)*100./float(nelements)),100),'%' ! a real representing progress

  read(fgmsh,'(a)',iostat=error) textline
  if (error /= 0) call error_stop('problem reading elements in gmsh file '//trim(filename))
  read(textline,*,iostat=error) gelement, gtype, ntags
  if (error /= 0) call error_stop('problem reading initial element data in gmsh file '//trim(filename))
  if (gtype > ubound(gtype_list,1) .or. .not.gtype_list(gtype)%supported) stop &
    'ERROR: element type not supported in gmsh file'

! get tags and local_nodes arrays to the correct size ready for reading
  call resize_integer_array(array=tags,new_size=ntags)
  call resize_integer_array(array=local_gnodes,new_size=gtype_list(gtype)%nnodes)
  call resize_integer_array(array=local_knodes,new_size=gtype_list(gtype)%nnodes)

! re-read line now knowing number of tags and number of nodes
  read(textline,*,iostat=error) gelement, gtype, ntags, tags, local_gnodes
  if (error /= 0) call error_stop('problem reading element tags/nodes in gmsh file '//trim(filename))
! convert local gnode list to knode list
  do k = 1, ubound(local_gnodes,1)
    local_knodes(k) = gmesh(gmesh_number)%knode_from_gnode(local_gnodes(k))
  end do

  if (ubound(icell_from_gelement,1) < gelement) then
    change = max(gelement-ubound(icell_from_gelement,1),int(nelements/10)+1)
    call resize_integer_array(keep_data=.true.,array=icell_from_gelement,change=change)
    call resize_integer_array(keep_data=.true.,array=jface_from_gelement,change=change)
  end if
  gelement_max = max(gelement_max,gelement)

! determine based on region and dimension info whether element is a cell or face
! default is a cell if the element has the maximum dimensions, otherwise a face
! note, cells should always be a part of a gmsh physical entity originally, and hence always be associated with a region,
!  so if an element is not associated with a region then it must be a face
  if (gtype_list(gtype)%dimensions == gmesh(gmesh_number)%dimensions) then
    centring = 'cell'
  else
    centring = 'face'
  end if
! now check on the parent region, overwriting above defaults if necessary
  region_number = 0
  if (tags(1) /= 0.and.tags(1) <= ubound(gmesh(gmesh_number)%gregion,1)) then
    region_number = gmesh(gmesh_number)%gregion(tags(1))%region_number
    centring = gmesh(gmesh_number)%gregion(tags(1))%centring
  end if

  if (debug) write(81,*) 'NEW ELEMENT: gelement = ',gelement,': gtype = ',gtype,': name = ',trim(gtype_list(gtype)%name), &
    ': ntags = ',ntags,': tags = ',tags,': local_gnodes = ',local_gnodes,': local_knodes = ',local_knodes, &
    ': gtype_dimensions = ',gtype_list(gtype)%dimensions,': gmesh_dimensions = ',gmesh(gmesh_number)%dimensions, &
    ': centring = ',centring

!------------------
! import as a cell

  if (centring == 'cell') then

    nfaces = gtype_list(gtype)%nfaces
    if (nfaces == 0) call error_stop('ERROR: a cell has no faces in the gmsh file '//trim(filename))
    call reset_cell(default_cell)
    call resize_integer_array(array=default_cell%jface,new_size=nfaces)

! loop through face node lists
    do nface = 1, nfaces

      call reset_face(default_face)
      call resize_integer_array(array=default_face%knode,new_size=ubound(gtype_list(gtype)%face_nodes(nface)%node,1))

! assemble face_node list
      do nnode = 1, ubound(default_face%knode,1)
        default_face%knode(nnode) = local_knodes(gtype_list(gtype)%face_nodes(nface)%node(nnode))
      end do

      if (debug) write(81,*) '  nface = ',nface,': default_face%knode = ',default_face%knode

! find face number if this face has been previously defined
      j = jface_from_knode_list(default_face%knode)

! if face is not found then create new one
      if (j == 0) then
        jtotal = jtotal + 1
        if (ubound(face,1) < jtotal) call resize_face_array(change=int(ktotal/10)+1)
        j = jtotal
! find gtype for face that has correct dimension and number of nodes
        default_face%gtype = 0
        do m = 1, ubound(gtype_list,1)
          if (gtype_list(m)%dimensions == gtype_list(gtype)%dimensions - 1.and. &
            gtype_list(m)%nnodes == ubound(default_face%knode,1) ) then
            default_face%gtype = m
            exit
          end if
        end do
        if (default_face%gtype == 0) call error_stop('gtype not identified for generated face in gmsh file '//trim(filename))

! set node jface references for the new face
        call add_jface_to_nodes(j,default_face%knode)
          
        if (debug) write(81,*) 'creating new face: j = ',j,': knode = ',default_face%knode,': gtype = ',default_face%gtype
      else
        default_face%gtype = face(j)%gtype
        if (debug) write(81,*) 'overwriting existing face: j = ',j,': knode = ',default_face%knode
      end if

      default_face%icell = [ 0, 0 ] ! reset icell(1) later, scrap icell(2)

! set or overwrite the face data to ensure that this cell is the icell=1 one
      default_face%dimensions = gtype_list(default_face%gtype)%dimensions
      call set_face(face_to_set=face(j),new_value=default_face)

      default_cell%jface(nface) = j

    end do

! now check that cell doesn't exist already from knowledge of nodes
    i = icell_from_knode_list(local_knodes)

    if (i == 0) then
      default_cell%gtype = gtype
      default_cell%dimensions = gtype_list(gtype)%dimensions
      call copy_integer_array(original=local_knodes,copy=default_cell%knode)
      itotal = itotal + 1
      if (ubound(cell,1) < itotal) call resize_cell_array(change=int(ktotal/5)+1)
      i = itotal
      call set_cell(cell_to_set=cell(i),new_value=default_cell)
! set node icell references for the new cell
      call add_icell_to_nodes(i,default_cell%knode)
      if (debug) write(81,*) 'creating new cell: i = ',i,': jface = ',default_cell%jface
    end if

! add gmesh lookup reference
    icell_from_gelement(gelement)=i

! set icell(1) so that the surrounding face normals point in the direction of the right-hand rule for each set of face nodes
    do jj = 1, ubound(cell(i)%jface,1)
      j = cell(i)%jface(jj)
      face(j)%icell(1) = i
    end do

! add cell number to region list if it isn't there already
    if (region_number /= 0) then
      if (debug) write(81,*) 'cell i = ',i,': is being added to region = ',trim(region(region_number)%name), &
        ': cell dimensions = ',cell(i)%dimensions,': region dimensions = ',region(region_number)%dimensions
      if (region(region_number)%dimensions /= cell(i)%dimensions) call error_stop('a cell is being added to region '// &
        trim(region(region_number)%name)//' however the dimensions are inconsistent')
      if (region(region_number)%centring /= 'cell') call error_stop('a cell is being added to region '// &
        trim(region(region_number)%name)//' however the centring is inconsistent')
      if (location_in_list(array=region(region_number)%ij,element=i) == 0) &
        call push_integer_array(array=region(region_number)%ij,new_element=i)
    end if

!------------------
! import as a face

  else if (centring == 'face') then

! find face number if this face has been previously defined
    j = jface_from_knode_list(local_knodes)

! if face is not found then create new one
    if (j == 0) then
      call reset_face(default_face)
      call copy_integer_array(original=local_knodes,copy=default_face%knode)
      default_face%gtype = gtype
      default_face%dimensions = gtype_list(gtype)%dimensions
      default_face%icell = [ 0, 0 ]
      jtotal = jtotal + 1
      if (ubound(face,1) < jtotal) call resize_face_array(change=int(ktotal/5)+1)
      j = jtotal
      call set_face(face_to_set=face(j),new_value=default_face)
      if (debug) write(81,*) 'creating new face: j = ',j,': knode = ',default_face%knode
! set node jface references for the new face
      call add_jface_to_nodes(j,default_face%knode)
    end if

! add gmesh lookup reference
    jface_from_gelement(gelement)=j

! add face number to region list if it isn't there already
    if (region_number /= 0) then
      if (debug) write(81,*) 'face j = ',j,': is being added to region = ',trim(region(region_number)%name), &
        ': face dimensions = ',face(j)%dimensions,': region dimensions = ',region(region_number)%dimensions
      if (region(region_number)%dimensions /= face(j)%dimensions) call error_stop('a face is being added to region '// &
        trim(region(region_number)%name)//' however the dimensions are inconsistent')
      if (region(region_number)%centring /= 'face') call error_stop('a face is being added to region '// &
        trim(region(region_number)%name)//' however the centring is inconsistent')
      if (location_in_list(array=region(region_number)%ij,element=j) == 0) &
        call push_integer_array(array=region(region_number)%ij,new_element=j)
    end if

!------------------
  end if

end do

! minimise cell and face arrays sizes
if (ubound(cell,1) /= itotal) call resize_cell_array(new_size=itotal)
if (ubound(face,1) /= jtotal) call resize_face_array(new_size=jtotal)

! allocate real gelement lookup arrays
! first copy info already in gelement array (from previous cell mesh reading for instance) into local variables
if (allocated(gmesh(gmesh_number)%gelement)) then
  gelement_previous = ubound(gmesh(gmesh_number)%gelement,1)
  if (gelement_previous > gelement_max) then
    call resize_integer_array(keep_data=.true.,array=icell_from_gelement,new_size=gelement_previous)
    call resize_integer_array(keep_data=.true.,array=jface_from_gelement,new_size=gelement_previous)
  end if
  do gelement = 1, gelement_previous
    if (gmesh(gmesh_number)%gelement(gelement)%icell /= 0) then
      if (icell_from_gelement(gelement) == 0) then
        icell_from_gelement(gelement) = gmesh(gmesh_number)%gelement(gelement)%icell
      else if (icell_from_gelement(gelement) /= gmesh(gmesh_number)%gelement(gelement)%icell) then
        call error_stop("cell element mismatch between elements from centringinput for "//trim(filename))
      end if
    end if
    if (gmesh(gmesh_number)%gelement(gelement)%jface /= 0) then
      if (jface_from_gelement(gelement) == 0) then
        jface_from_gelement(gelement) = gmesh(gmesh_number)%gelement(gelement)%jface
      else if (jface_from_gelement(gelement) /= gmesh(gmesh_number)%gelement(gelement)%jface) then
        call error_stop("face element mismatch between elements from centringinput for "//trim(filename))
      end if
    end if
  end do
  gelement_max = max(gelement_max,gelement_previous)
  deallocate(gmesh(gmesh_number)%gelement)
end if

! now copy local arrays to global ones
allocate(gmesh(gmesh_number)%gelement(gelement_max))
do gelement = 1, gelement_max
  gmesh(gmesh_number)%gelement(gelement)%icell = icell_from_gelement(gelement)
  gmesh(gmesh_number)%gelement(gelement)%jface = jface_from_gelement(gelement)
end do
deallocate(icell_from_gelement,jface_from_gelement)

end subroutine read_gmesh_elements

!-----------------------------------------------------------------

subroutine read_gmesh_nodes(check,gmesh_number,filename,debug)

use general_module
use gmsh_module
integer :: error, k, nnodes, gnode_check, gmesh_number, gnode_max
integer, dimension(:), allocatable :: gnode
character(len=1000) :: textline, filename
double precision, dimension(totaldimensions) :: x_check
logical :: check, debug

rewind(fgmsh)
do 
  read(fgmsh,'(a)',iostat=error) textline
! if (error /= 0) call error_stop('list of nodes not found in gmsh file '//trim(filename))
  if (error /= 0) then
    write(*,'(a)') "WARNING: list of nodes not found in gmsh file "//trim(filename)
    return
  end if
  if (trim(textline) == "$Nodes") exit
end do
read(fgmsh,*,iostat=error) nnodes
if (error /= 0) call error_stop('problem reading in number of nodes in gmsh file '//trim(filename))
if (debug) write(*,*) 'number of nodes: nnodes = ',nnodes
if (.not.check) then
  call resize_node_array(change=nnodes)
  allocate(gnode(nnodes))
  do k = 1, nnodes
    read(fgmsh,*,iostat=error) gnode(k),node(k+ktotal)%x
    if (error /= 0) call error_stop('problem reading nodes in gmsh file '//trim(filename))
  end do
  ! now create lookup array
  allocate(gmesh(gmesh_number)%knode_from_gnode(maxval(gnode)))
  gmesh(gmesh_number)%knode_from_gnode = 0
  do k = 1, nnodes
    gmesh(gmesh_number)%knode_from_gnode(gnode(k)) = k + ktotal
  end do
  deallocate(gnode)
  ktotal = ktotal + nnodes
else
  gnode_max = 0
  do k = 1, nnodes
    read(fgmsh,*,iostat=error) gnode_check,x_check
    if (error /= 0) call error_stop('problem reading nodes in gmsh file '//trim(filename))
    if (distance(x_check,node(gmesh(gmesh_number)%knode_from_gnode(gnode_check))%x) > 1.d-10) call error_stop( &
      "node mismatch between msh files: cell and face msh files must come from the same simulation")
    gnode_max = max(gnode_max,gnode_check)
  end do
  if (ubound(gmesh(gmesh_number)%knode_from_gnode,1) /= gnode_max) call error_stop("number of nodes mismatched between msh "// &
    "files:  cell and face msh files must come from the same simulation")
end if

end subroutine read_gmesh_nodes

!-----------------------------------------------------------------

subroutine read_gmesh_regions(check,gmesh_number,filename,debug)

use general_module
use gmsh_module
integer :: error, ngregions, n, gmesh_number, region_number
integer, dimension(:), allocatable :: gregion_dimensions, gregion_number
character(len=1000), dimension(:), allocatable :: gregion_name
character(len=1000) :: formatline, textline, filename, location
character(len=4) :: centring
logical :: check, debug, existing

rewind(fgmsh)
do 
  read(fgmsh,'(a)',iostat=error) textline
! if (error /= 0) call error_stop('list of physical names not found in gmsh file '//trim(filename))
  if (error /= 0) then
    write(*,'(a)') "WARNING: list of physical names not found in gmsh file "//trim(filename)
    return
  end if
  if (trim(textline) == "$PhysicalNames") exit
end do
read(fgmsh,*,iostat=error) ngregions
if (error /= 0) call error_stop('problem reading in number of physical names in gmsh file '//trim(filename))
if (debug) write(*,*) 'number of physical names: ngregions = ',ngregions

! run through physical entities finding dimension and gregion_numbers
allocate(gregion_number(ngregions),gregion_dimensions(ngregions),gregion_name(ngregions))
do n=1,ngregions
  read(fgmsh,*,iostat=error) gregion_dimensions(n),gregion_number(n),gregion_name(n)
  if (debug) write(*,*) gregion_dimensions(n), gregion_number(n), trim(gregion_name(n))
  if (error /= 0) call error_stop('problem reading physical name in gmsh file '//trim(filename))
end do
if (debug) write(*,*) ' gmesh dimensions = ',maxval(gregion_dimensions)
if (check) then
  if (gmesh(gmesh_number)%dimensions /= maxval(gregion_dimensions)) call error_stop( &
    "dimension mismatch between msh files: cell and face msh files must come from the same simulation")
else
  gmesh(gmesh_number)%dimensions = maxval(gregion_dimensions)
end if

if (.not.check) then
  allocate(gmesh(gmesh_number)%gregion(maxval(gregion_number))) ! create a listing of the gmsh physical regions referenced by their gmsh number
  gmesh(gmesh_number)%gregion(:)%region_number = 0 ! default value of 0 indicates this region is not allocated to an arb region
  gmesh(gmesh_number)%gregion(:)%centring = ''
end if
do n=1,ngregions
! region existence possibilities
! 1) nothing has been set: centring is calculated from dimension related to the rest of the mesh file
! 2) centring has been set from constants.in but not location: this centring is to be applied to elements read in from gmsh file
! 3) centring and location has been set from previous GMSH read: check on consistency of centring and dimensions
! 4) centring and location has been set from constants.in: gmsh element data is to be ignored

  region_number = region_number_from_name(name=gregion_name(n),creatable=.false.)
  formatline = '(a,'//trim(dindexformat(gregion_number(n)))//',a)'
  write(location,fmt=formatline) 'GMSH physical entity number ',gregion_number(n),' from file '//trim(filename)

  if (region_number == 0) then
! region does not yet exist: determine centring and create new region
    if (check) call error_stop('region '//trim(gregion_name(n))//' in gmsh file '//trim(filename)// &
      ' is not properly defined on current read: all centring files with the same basename must come from the '// &
      'same simulation and have the same physical entities defined')
    if (debug) write(*,*) ' before read: region '//trim(gregion_name(n))//': not previously defined'
    if (gregion_dimensions(n) == gmesh(gmesh_number)%dimensions) then
      centring = "cell"
    else if (gregion_dimensions(n) <= gmesh(gmesh_number)%dimensions-1) then
      centring = "face"
    end if
    if (gregion_dimensions(n) < gmesh(gmesh_number)%dimensions-1) &
      write(*,'(a)') 'WARNING: region '//trim(gregion_name(n))//' in gmsh file '//trim(filename)// &
      ' assumed to have face centring based on dimensions relative to the other elements in this file'
! create new region and set gregion region_number at the same time
    gmesh(gmesh_number)%gregion(gregion_number(n))%region_number = region_number_from_name(name=gregion_name(n), &
      location=location,centring=centring,dimensions=gregion_dimensions(n),creatable=.true.)
    gmesh(gmesh_number)%gregion(gregion_number(n))%centring = centring
    region_number = gmesh(gmesh_number)%gregion(gregion_number(n))%region_number

  else if (trim(region(region_number)%centring) == "") then
! region is defined but not centring: this is an error
    call error_stop('region '//trim(gregion_name(n))//' in gmsh file '//trim(filename)//' is already allocated but'// &
      ' has no centring defined')

  else if (trim(region(region_number)%location(1:4)) == "GMSH") then
    if (check) then
      if (gmesh(gmesh_number)%gregion(gregion_number(n))%region_number /= region_number) call error_stop('region '// &
        trim(gregion_name(n))//' in gmsh file '//trim(filename)//' has inconsistent region_number on reread')
    else
      if (debug) write(*,*) ' before read: region '//trim(gregion_name(n))//': region number ',region_number, &
        ': centring = '//trim(region(region_number)%centring)//': location '//trim(region(region_number)%location)// &
        ': dimensions = ',region(region_number)%dimensions
      write(*,'(a)') 'WARNING: region '//trim(gregion_name(n))//' is contained in more than one gmsh file:'
      write(*,'(a)') ' original location: '//trim(region(region_number)%location)
      write(*,'(a)') ' current location: '//trim(location)
! if this gregion has more dimensions then the dimensions of the region expands
!     region(region_number)%dimensions = max(gregion_dimensions(n),region(region_number)%dimensions) 
! actually this should flag a problem as each element associated with a gmsh physical entity should have the same dimension
      if (region(region_number)%dimensions /= gregion_dimensions(n)) call error_stop('region '//trim(gregion_name(n))// &
        ' in gmsh file '//trim(filename)//' has dimensions that are inconsistent with a previous gmsh definition of this region')
! this allows data from multiple mesh files to be read into the same region
      gmesh(gmesh_number)%gregion(gregion_number(n))%region_number = region_number
      gmesh(gmesh_number)%gregion(gregion_number(n))%centring = region(region_number)%centring
    end if

  else if (trim(region(region_number)%location) == "") then
! only centring has been set: set location and dimensions
    if (check) call error_stop('region '//trim(gregion_name(n))//' in gmsh file '//trim(filename)// &
      ' is not properly defined on current read: all centring files with the same basename must come from the '// &
      'same simulation and have the same physical entities defined')
    if (debug) write(*,*) ' before read: region '//trim(gregion_name(n))//': region number ',region_number, &
      ': centring = '//trim(region(region_number)%centring)//': location '//trim(region(region_number)%location)// &
      ': dimensions = ',region(region_number)%dimensions
    if (region(region_number)%dimensions /= -1) call error_stop('region '//trim(gregion_name(n))//' in gmsh '// &
      'file '//trim(filename)//' that is about to be read already has dimensions set from some previous definition')
    region(region_number)%location = trim(location)
    region(region_number)%dimensions = gregion_dimensions(n)
    gmesh(gmesh_number)%gregion(gregion_number(n))%region_number = region_number
    gmesh(gmesh_number)%gregion(gregion_number(n))%centring = region(region_number)%centring

  else
! this is a region that is to be defined using user statements: gmsh definition will be ignored to allow reread
    if (.not.check) write(*,'(a)') "NOTE: region "//trim(gregion_name(n))//" defined in the gmsh file "//trim(filename)// &
      " conflicts with a region definition given in constants.in:  the region definition in the current gmsh file will be ignored"
! however centring is recorded so that the centring of cells associated with this gregion can be correctly identified
    gmesh(gmesh_number)%gregion(gregion_number(n))%centring = region(region_number)%centring

  end if

  if (debug) write(*,*) ' after read: region '//trim(gregion_name(n))//': region number ',region_number, &
    ': centring = '//trim(region(region_number)%centring)//': location '//trim(region(region_number)%location)// &
    ': dimensions = ',region(region_number)%dimensions

end do
deallocate(gregion_number,gregion_dimensions,gregion_name)

end subroutine read_gmesh_regions

!-----------------------------------------------------------------

end module input_module

!-----------------------------------------------------------------
