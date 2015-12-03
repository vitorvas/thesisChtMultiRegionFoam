! file src/kernel_module.f90
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
module kernel_module

implicit none

! only the setup_kernels subroutine is accessible from outside this module
private
public setup_kernels

! type for weight fluxing
type weight_flux_type
  integer :: ndonating ! the number of faces that are donating from this cell
  integer, dimension(:), allocatable :: donating ! the ii index of the cell that this face donates to
  double precision, dimension(:), allocatable :: proportion ! proportion of flux from cell ii in kernel that will be fluxed through this (cell(i)%face) face
end type weight_flux_type

! specifications for the kernels follow:
! both mls and optimisation methods:
character(len=100), parameter :: kernel_method = 'optimisation' ! polynomial basis kernels calculated via nonlinear optimisation
!character(len=100), parameter :: kernel_method = 'mls' ! polynomial basis kernels calculated via moving least squares technique
!character(len=100), parameter :: kernel_method = 'simple' ! simple over face differencing, only suitable for 1D problems.
logical, parameter :: cell_from_face_kernels = .false. ! calculate the cell centred derivatives by averaging the face centred derivatives to the cell centre.  This seems to work well, both from a computational efficiency and stability perspective.
logical, parameter :: face_relative_from_absolute = .true. ! calculate the face centred derivative kernels that are relative to the face normal from the already-calculated coordinate index kernels.  Again, seems to be a good idea  both from a computational efficiency and stability perspective.  Infact unless the kernel masks or weighting factors are different, the kernels calculated this way will be exactly the same as if calculated independently.
integer, parameter :: minimum_domain_separation = 2, minimum_boundary_separation = 2 ! these are the default separations for domain/boundary faces that specify the size of the kernels.  Large numbers produce large kernels.
integer, parameter :: maximum_domain_separation = 3, maximum_boundary_separation = 3 ! these are the maximum separations for domain/boundary faces that specify the size of the kernels.  These are the maximum sizes that can be used if adaptive_mask_size is on.
logical, parameter :: limit_mask_to_icell = .false. ! cells within kernel mask must come from icell - this is a locality constraint that all cells in mask must share at least one node with the relevant face
integer, parameter :: polynomial_order = 1 ! maximum order (power index) within polynomial basis functions
double precision, parameter :: kernel_h_multiplier = 1.0d0 ! multiple h_kernel by this
double precision, parameter :: weight_separation_amount = 0.1d0 ! for non-centred kernels this specifies the how the weights of each separation are related - the smaller the number, the tighter the kernel - use is slightly different for conservative and nonconservative weighting - for conservative_weighting this is the proportion of each cell weight that is fluxed to further away neighbours
logical, parameter :: remove_small_elements = .true. ! kernel values below small_element_minimum will be removed
double precision, parameter :: small_element_minimum = 1.d-10 ! minimum kernel size allowed (otherwise set to zero to save memory)
logical, parameter :: conservative_weighting = .false. ! use a conservation principle to calculate weights that are connected, otherwise use formula based on absolute separation
logical, parameter :: orientation_dependent_weights = .true. ! the kernel weights are difference depending on the direction of the kernel
logical, parameter :: radial_kernel_weighting = .true. ! the kernel weighting is also influenced by the radial distance
logical, parameter :: uniform_cell_averaging_kernels = .true. ! the face->cell and node->cell kernels are uniform - seems to happen anyway using a linear fit, but have not proven this generally

! mls options:
logical, parameter :: mls_check_minw = .true. ! check that the minw value is large enough
double precision, parameter :: mls_minimum_minw = 0.2d0 ! minimum value of SVD minw allowed for mask to be acceptable when using adaptive_mask_size

! optimisation options:
logical, parameter :: optimise_positise = .false. ! constrain kernel values in an attempt at getting positive kernels
logical, parameter :: optimise_positise_cautiously = .true. ! only allow change to kernel if it directly reduces negative index of kernel
double precision, parameter :: optimise_positise_cautiously_multiplier = 1.d0 ! when optimising cautiously, allow changes if negative index is less than this multiple of last negative index - a large number reproduces a non-cautious approach - 1.d0 is a cautious approach - less than one will be super cautious (ie, the kernel will not be changed much)
integer, parameter :: optimise_additional_elements = 0 ! setting this to zero will mean that each kernel has at least as many elements present as in the minimum separation level - a large negative number will give the smallest kernel, a positive number will increase the kernel size (really needs positise on for /= 0 on this variable - setting to a negative number is not advisable as kernel may not be structurally symmetric)
logical, parameter :: optimise_positise_sum_index = .true. ! base negative index (which measures kernel positivity) on sum of negative elements, rather than maximum negative element
double precision, parameter :: small_pp = 1.d-8 ! used as a cut-off for elements of the polynomial basis

! kernel usage notes:
! 1) conservative weighting does a better job for both mls and optimisation kernels as evidenced by lower errors, but there isn't much inbetween the conservative and nonconservative methods
! 2) zero_wayward_boundary_weights only needs to be set for mls, or optimisation with optimise_positise off
!-----------------------------------------------------------------
contains

!-----------------------------------------------------------------

subroutine setup_kernels

! here we calculate the face and cell centred kernels

use general_module
integer :: i, j, l, ii, jj, ij, ii2, i2, kk, ierror, n, kk2, l2
integer :: n_domain_kernels, n_domain_elements, n_boundary_kernels, n_boundary_elements, n_elements, &
  min_location, max_location, maximum_separation, minimum_separation, separation
double precision :: dx1, dx2, h_kernel, min_value, max_value, value, minw
double precision, dimension(:,:), allocatable :: r, norm, pp
double precision, dimension(:), allocatable :: kernel_error
integer, dimension(:), allocatable :: separation_index, separation_array
double precision, dimension(:,:), allocatable :: max_rel_face_kernel ! maximum of separation value / minimum of central values
integer, dimension(:,:), allocatable :: max_rel_face_jface
logical :: error
character(len=10000) :: formatline, filename
logical, parameter :: debug = .false.
logical :: debug_sparse = .false.

if (debug) debug_sparse = .true.

if (debug_sparse) write(*,'(80(1h+)/a)') 'subroutine setup_kernels'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 13/06/2013
!filename = "output/kernel_warnings.txt"
filename = "kernel_warnings.txt"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(fwarn,file=trim(filename),status='replace',iostat=ierror)
if (ierror /= 0) then
  write(*,*) 'ERROR: problem opening file ',filename
  stop
end if

!------------------------------------------
! run some checks on selected options
if (trim(kernel_method) == 'optimisation') then
  if (.not.face_relative_from_absolute) call error_stop('face_relative_from_absolute must be true for '//trim(kernel_method)// &
    ' kernels')
  if (.not.radial_kernel_weighting) write(*,'(a)') 'WARNING: inadvisable combination of optimisation kernel method without '// &
    'radial_kernel_weighting'
end if

!------------------------------------------
! run through cells and faces setting kernels
! face kernels:
! kernel(0) is the average based on surrounding cell values
! kernel(1-3) is the derivative in the 1-3rd direction based on surrounding cell values
! kernel(4-6) is the derivative in the face normal directions based on surrounding cell values
! cell kernels:
! kernel(0) is the average based on surrounding face values
! kernel(1-3) is the derivative in the 1-3rd direction based on surrounding cell values
! kernel(4) is the average based on surrounding node values

!------------------------------------------
! setting up face kernels

if (debug_sparse) write(*,'(a)') ' constructing face kernels using '//trim(kernel_method)//' method'

! mls and optimisation kernels
if (trim(kernel_method) /= 'simple') then

! zero separation level specific kernel maximums
  allocate(max_rel_face_jface(0:6,2:max(maximum_domain_separation,maximum_boundary_separation)))
  allocate(max_rel_face_kernel(0:6,2:max(maximum_domain_separation,maximum_boundary_separation)))
  max_rel_face_jface = 0
  max_rel_face_kernel = 0.d0

! temp &&&&
  do j = 1, jtotal
! do j = 1, 100 
! do j = 3954, 3954
! do j = 3174, 3174
! do j = 38165,38165
! do j = 13700,13700

    if (debug) write(83,*) '----------------------------'
    if (debug) write(83,*) 'FACE: j = ',j,': face type = ',face(j)%type,': face dimensions = ',face(j)%dimensions

! find h_kernel for this face which is independent of kernel mask and direction
! TODO: put in alternative h_kernel based on volume

    if (ubound(face(j)%knode,1) > 1) then
      h_kernel = 1.d+20
      do kk = 1, ubound(face(j)%knode,1)
        do kk2 = kk+1, ubound(face(j)%knode,1)
          h_kernel = min(h_kernel,distance(node(face(j)%knode(kk))%x , node(face(j)%knode(kk2))%x))
        end do
      end do
      h_kernel = h_kernel/2.d0
    else
      h_kernel = face(j)%dx
      if (face(j)%type == 1) h_kernel = h_kernel/2.d0 ! if not a boundary face need to divide this by 2
    end if
    h_kernel = kernel_h_multiplier*h_kernel
    if (debug) write(83,*) 'h_kernel = ',h_kernel

! set the (minimum) default separations
    if (face(j)%type == 2) then
      minimum_separation = minimum_boundary_separation
      maximum_separation = maximum_boundary_separation
    else 
      minimum_separation = minimum_domain_separation
      maximum_separation = maximum_domain_separation
    end if

! setup the kernel mask which is the same for all kernel directions

! include first two elements and assign their separations locally (specific to the face)
! make sure that first two elements are as per icell so that boundary values correctly applied
    face(j)%kernel(0)%centring = 'cell'
    call resize_integer_array(keep_data=.false.,array=face(j)%kernel(0)%ijk,new_size=2)
    face(j)%kernel(0)%ijk(1:2) = face(j)%icell(1:2)
    call resize_integer_array(keep_data=.false.,array=separation_index,new_size=1)
    separation_index(1) = 2 ! last index in kernel%ijk that has a cell with separation 1
    call resize_integer_array(keep_data=.false.,array=separation_array,new_size=2)
    separation_array = 1

! add elements to the kernel mask in increasing order of separation up to the maximum_separation, and update separation_index &
!  array accordingly
    call expand_kernel_mask(iarray=face(j)%icell,maximum_separation=maximum_separation,imask=face(j)%kernel(0)%ijk, &
      separation_index=separation_index,separation_array=separation_array)

! also size value array
    call resize_double_precision_array(keep_data=.false.,array=face(j)%kernel(0)%v,new_size=ubound(face(j)%kernel(0)%ijk,1))

! create r for all cells in the mask and scale it with h_kernel
    if (allocated(r)) deallocate(r)
    allocate(r(totaldimensions,ubound(face(j)%kernel(0)%ijk,1)))
    do ii = 1, ubound(face(j)%kernel(0)%ijk,1)
      r(:,ii) = cell(face(j)%kernel(0)%ijk(ii))%x - face(j)%x
    end do
    r = r/h_kernel

! and then convert to a consistent basis
! construct norm, find an orthogonal basis for r and convert r and the norm to this basis
    if (allocated(norm)) deallocate(norm)
    allocate(norm(totaldimensions,2*totaldimensions))
    norm = 0.d0
    norm(1,1) = 1.d0
    norm(2,2) = 1.d0
    norm(3,3) = 1.d0
    norm(:,4) = face(j)%norm(:,1)
    norm(:,5) = face(j)%norm(:,2)
    norm(:,6) = face(j)%norm(:,3)
    call construct_orthogonal_basis('face',r=r,norm=norm,error=error)
    if (error) call error_stop('unable to construct orthogonal basis vectors for face kernel')

! calculate polynomial basis pp tensor from list of r vectors
    call construct_polynomial_basis_tensor(r,polynomial_order,pp,error)
    if (error) call error_stop('unable to construct pp basis tensor for face kernel')

! check minw, enlarging the minimum_separation if required
    if (trim(kernel_method) == 'mls'.and.mls_check_minw) &
      call check_mask_minw(pp,separation_index,minimum_separation,minw)

! loop through all the directions required, doing face relative directions first

    face_direction_loop: do l = 6, 0, -1

      if (debug) write(83,*) 'START direction_loop: l = ',l

! copy and reset kernel
      if (l /= 0) call copy_kernel(original=face(j)%kernel(0),copy=face(j)%kernel(l))
      face(j)%kernel(l)%v = 0.d0

      if (l == 0.and.face(j)%type == 2) then
! for boundary cells averaging kernel don't do mls
        face(j)%kernel(0)%v(2) = 1.d0
        if (debug) write(83,*) 'boundary averaging kernel: type = ',face(j)%type

      else if (l >= 1.and.vector_magnitude(norm(:,l)) < 1.d-10) then
! if the norm is zero in this direction don't do either
        if (debug) then
          write(83,'(a)') 'norm component when expressed in basis is zero: skipping mls kernel construction'
          write(83,*) 'l = ',l,': norm(:,l) = ',norm(:,l),': vector_magnitude(norm(:,l)) = ',vector_magnitude(norm(:,l)) 
        end if

      else if (l >= 1.and.l <= 3.and.face_relative_from_absolute) then ! if doing directions in reverse
! construct these absolute kernels (1->3) from the relative ones calculated earlier (4->6)
        if (debug) write(83,*) 'constructing kernel from previous derivative kernels: l = ',l
        do ii = 1, ubound(face(j)%kernel(l)%ijk,1)
          do l2 = 4,6
            face(j)%kernel(l)%v(ii) = face(j)%kernel(l)%v(ii) + face(j)%norm(l,l2-3)*face(j)%kernel(l2)%v(ii)
          end do
        end do

!       else if (l >= 4.and.face_relative_from_absolute) then ! if doing directions forward
! construct these relative kernels (4->6) from the absolute ones calculated earlier (1->3)
!         if (debug) write(83,*) 'constructing kernel from previous derivative kernels: l = ',l
!         do ii = 1, ubound(face(j)%kernel(l)%ijk,1)
!           do l2 = 1,3
!             face(j)%kernel(l)%v(ii) = face(j)%kernel(l)%v(ii) + face(j)%norm(l2,l-3)*face(j)%kernel(l2)%v(ii)
!           end do
!         end do

      else
! create kernels via mls or optimisation method

        if (debug) write(83,*) 'calculating kernels ',l,' for face ',j,'j via '//trim(kernel_method)//' method'

        if (trim(kernel_method) == 'mls') then
          if (l == 0) then
            call mls_kernel(centring='face',ijk=j,l_kernel=l,rr=r,pp=pp,kernel=face(j)%kernel(l)%v, &
              local_polynomial_order=polynomial_order,minimum_separation=minimum_separation,separation_array=separation_array, &
              separation_index=separation_index,error=error)
          else
            call mls_kernel(centring='face',ijk=j,l_kernel=l,rr=r,norm=norm(:,l),pp=pp,kernel=face(j)%kernel(l)%v, &
              local_polynomial_order=polynomial_order,minimum_separation=minimum_separation,separation_array=separation_array, &
              separation_index=separation_index,error=error)
          end if
        else
          if (l == 0) then
            call optimisation_kernel(centring='face',ijk=j,l_kernel=l,rr=r,pp=pp,kernel=face(j)%kernel(l)%v, &
              local_polynomial_order=polynomial_order,minimum_separation=minimum_separation,separation_array=separation_array, &
              separation_index=separation_index,error=error)
          else if (l >= 4.and.l <= 6) then
            call optimisation_kernel(centring='face',ijk=j,l_kernel=l,rr=r,norm=norm(:,l),pp=pp,kernel=face(j)%kernel(l)%v, &
              local_polynomial_order=polynomial_order,minimum_separation=minimum_separation,separation_array=separation_array, &
              separation_index=separation_index,error=error)
          else
            call error_stop('optimisation kernel require face_relative_from_absolute to be true')
          end if
        end if
        if (error) call error_stop('error in calculating a face '//trim(kernel_method)//' kernel')

      end if
  
    end do face_direction_loop

! rescale derivative kernels
    do l = 1, 6
      face(j)%kernel(l)%v=face(j)%kernel(l)%v/h_kernel
    end do

! find maximum values for each separation level
    do l = 0, 6
      do separation = 2, ubound(separation_index,1)
        value = maxval(abs(face(j)%kernel(l)%v(separation_index(separation-1)+1:separation_index(separation))))/ &
          maxval(abs(face(j)%kernel(l)%v(1:2)))
        if (value > max_rel_face_kernel(l,separation)) then
          max_rel_face_kernel(l,separation) = value
          max_rel_face_jface(l,separation) = j
        end if
      end do
    end do
          
! print some debugging info about kernels

    if (debug_sparse) then
      formatline = '(a,'//trim(indexformat)//',a,i1,a,g9.2,a,i2)'
      write(83,fmt=formatline) 'END separation_loop: all kernels calculated for: j = ',j,'j: type = ',face(j)%type, &
        ': h_kernel = ',h_kernel,': minimum_separation = ',minimum_separation,': maximum_separation = ',maximum_separation
! print out details of all cells that are in the kernel
      do ii = 1, ubound(face(j)%kernel(0)%ijk,1)
        i = face(j)%kernel(0)%ijk(ii)
        formatline = '(a,i3,a,'//trim(dindexformat(i))//',a,i2,a,g9.2,a'//repeat(',1x,f7.3',7)//')'
        write(83,fmt=formatline) 'ii = ',ii,': i = ',i,': sep. = ',separation_array(ii),': rmag = ',vector_magnitude(r(:,ii)), &
          ': v = ',face(j)%kernel(0)%v(ii),(face(j)%kernel(l)%v(ii)*h_kernel,l=1,6)
      end do
      if (trim(kernel_method) == 'mls'.and.mls_check_minw) write(83,*) 'minw = ',minw

! temp &&&
!     if (.false..and.j == 12294) then
!       write(*,*) 'WARNING: writing out kernel debugging file for j = ',j
!       open(unit=84,file='kernel_debugging.msh')
!       write(84,'(a/a/a)') '$MeshFormat','2.2 0 8','$EndMeshFormat'
!       write(84,'(a,3(/a))') '$Nodes','1','1 0. 0. 0.','$EndNodes'
!       write(84,'(a/i2)') '$Elements',ubound(r,2)
!       do ii = 1, ubound(r,2)
!         write(84,'(i2,a)') ii,' 15 2 0 0 1'
!       end do
!       write(84,'(a)') '$EndElements'
!       if (ubound(r,1) /= 3) stop "only 3d vectors can be handles right now"

!       write(84,'(a,6(/a))') '$ElementData','1','"<r>"','0','3','0','3'
!       write(84,'(i2)') ubound(r,2)
!       do ii = 1, ubound(r,2)
!         write(84,'(i1,3(1x,f10.5))') ii,(real(r(l,ii)),l=1,3)
!       end do
!       write(84,'(a)') '$EndElementData'
!         
!       write(84,'(a,6(/a))') '$ElementData','1','"<facenorm>"','0','3','0','3'
!       write(84,'(i2)') 1
!       write(84,'(i1,3(1x,f10.5))') 1,(real(norm(l,4)),l=1,3)
!       write(84,'(a)') '$EndElementData'
!         
!       close(unit=84)
!     end if
    end if

  end do

! TODO: need to redo array bounds here to cope with separation=1
  if (.true.) then
    do l = 0, 6
      do separation = 2, ubound(separation_index,1)
        write(fwarn,'(a,i1,a,i1,a,g10.3,a,i8)') 'l = ',l,': separation = ',separation,': max_rel_face_kernel = ', &
          max_rel_face_kernel(l,separation),': j = ',max_rel_face_jface(l,separation)
      end do
    end do
  end if

  deallocate(max_rel_face_jface,max_rel_face_kernel)

!----------------------
! uber simple masks suitable for 1D applications only
! only 0 (average) and 4 (gradient in face direction) defined
else
  do j = 1, jtotal
    do l = 0, 6
      allocate(face(j)%kernel(l)%ijk(2),face(j)%kernel(l)%v(2))
      face(j)%kernel(l)%ijk = face(j)%icell(1:2)
      face(j)%kernel(l)%v = 0.d0
      if (l == 0.or.l == 4) then
        dx1 = abs(dot_product( cell(face(j)%icell(1))%x-face(j)%x , face(j)%norm(:,1) ))
        dx2 = abs(dot_product( cell(face(j)%icell(2))%x-face(j)%x , face(j)%norm(:,1) ))
        if (l == 0) then
          face(j)%kernel(0)%v(1) = dx2/(dx1+dx2)
          face(j)%kernel(0)%v(2) = dx1/(dx1+dx2)
        else if (l == 4) then
          face(j)%kernel(4)%v(1) = -1.d0/(dx1+dx2)
          face(j)%kernel(4)%v(2) = 1.d0/(dx1+dx2)
        end if
      end if
    end do
  end do
end if

if (allocated(r)) deallocate(r)
if (allocated(norm)) deallocate(norm)
if (allocated(separation_index)) deallocate(separation_index)
if (allocated(separation_array)) deallocate(separation_array)

!------------------------------------------
! setting up cell kernels

if (debug_sparse) write(*,'(a)') ' constructing cell kernels using '//trim(kernel_method)//' method'
if (trim(kernel_method) /= 'simple') then

  do i = 1, itotal

    if (debug_sparse) then
      write(83,*) '----------------------------'
      formatline = '(a,'//trim(indexformat)//',a,i1,a,i1)'
      write(83,fmt=formatline) 'CELL: i = ',i,'i: cell type = ',cell(i)%type,': cell dimensions = ',cell(i)%dimensions
      if (debug) write(83,*) 'surrounding faces: ',cell(i)%jface
    end if

! set h_kernel for all cell kernels based on maximum spacing between nodes on cell if cell has enough dimensions
! otherwise based on distance to neighbouring cell centroid
    if (cell(i)%dimensions > 0) then
      h_kernel = 1.d+20
      do kk = 1, ubound(cell(i)%knode,1)
        do kk2 = kk+1, ubound(cell(i)%knode,1)
          h_kernel = min(h_kernel,distance(node(cell(i)%knode(kk))%x , node(cell(i)%knode(kk2))%x))
        end do
      end do
      h_kernel = h_kernel/2.d0
    else ! this must be a boundary cell in a 1d domain, so use dx from corresponding boundary face 
      h_kernel = face(cell(i)%jface(1))%dx
      if (face(cell(i)%jface(1))%type /= 2) call error_stop('problem in setup_kernels')
    end if
    h_kernel = kernel_h_multiplier*h_kernel
    if (debug) write(83,*) 'h_kernel = ',h_kernel

    cell_direction_loop: do l = 0, 4

! set kernel centring
      if (l == 0) then
        cell(i)%kernel(l)%centring = 'face'
        call copy_integer_array(original=cell(i)%jface,copy=cell(i)%kernel(l)%ijk)
        allocate(cell(i)%kernel(l)%v(ubound(cell(i)%kernel(l)%ijk,1)))
        cell(i)%kernel(l)%v = 0.d0
      else if (l == 4) then
        cell(i)%kernel(l)%centring = 'node'
        call copy_integer_array(original=cell(i)%knode,copy=cell(i)%kernel(l)%ijk)
        allocate(cell(i)%kernel(l)%v(ubound(cell(i)%kernel(l)%ijk,1)))
        cell(i)%kernel(l)%v = 0.d0
      else
        cell(i)%kernel(l)%centring = 'cell'
      end if

      if (debug) then
        write(83,*)
        write(83,*) 'NEW KERNEL: i = ',i,': l = ',l,': kernel centring = '//trim(cell(i)%kernel(l)%centring)
      end if

!---------------
! boundary cells take on value from boundary face as set above as they are coincident
! also, in 1D boundary cells are coincident with boundary nodes also so take on that value
      if ((l == 0.or.(l == 4.and.cell(i)%dimensions == 0)).and.cell(i)%type == 2) then

        cell(i)%kernel(l)%v(1) = 1.d0 
        
!---------------
! boundary cells take on derivatives from boundary face as they are coincident
      else if (l >= 1.and.l <= 3.and.cell(i)%type == 2) then

        j = cell(i)%jface(1)
        call copy_kernel(original=face(j)%kernel(l),copy=cell(i)%kernel(l))
        cell(i)%kernel(l)%centring = 'cell' ! have to rewrite this
        cell(i)%kernel(l)%v = cell(i)%kernel(l)%v*h_kernel ! rescaling here so that all derivative kernels can be unscaled later
        if (debug) write(83,*) 'pulling face boundary kernel from j = ',j

!---------------
! the averaging kernels are just uniform, which seems (!) to be correct for the linear interpolation (which must be used here anyway due to number of points)
      else if (uniform_cell_averaging_kernels.and.l == 0.or.l == 4) then
        cell(i)%kernel(l)%v = 1.d0/ubound(cell(i)%kernel(l)%ijk,1)

!---------------
! construct derivative kernels from surrounding face derivative and cell averaging kernels
      else if (l >= 1.and.l <= 3.and.cell_from_face_kernels) then

! first create mask from surrounding face masks of the same derivative
        allocate(cell(i)%kernel(l)%ijk(1))
        cell(i)%kernel(l)%ijk(1) = i
        do jj = 1, ubound(cell(i)%kernel(0)%ijk,1) ! this is equivalent to cell(i)%jface
          j = cell(i)%kernel(0)%ijk(jj)
          do ii2 = 1, ubound(face(j)%kernel(l)%ijk,1)
            i2 = face(j)%kernel(l)%ijk(ii2)
            if (location_in_list(array=cell(i)%kernel(l)%ijk,element=i2) == 0) &
              call push_array(array=cell(i)%kernel(l)%ijk,new_element=i2) ! add element if it is not already on the list
          end do
        end do

! create value array and zero it
        allocate(cell(i)%kernel(l)%v(ubound(cell(i)%kernel(l)%ijk,1)))
        cell(i)%kernel(l)%v = 0.d0

! now use averaging kernel to create derivative entries
        do jj = 1, ubound(cell(i)%kernel(0)%ijk,1)
          j = cell(i)%kernel(0)%ijk(jj)
          do ii2 = 1, ubound(face(j)%kernel(l)%ijk,1)
            i2 = face(j)%kernel(l)%ijk(ii2)
            n = location_in_list(array=cell(i)%kernel(l)%ijk,element=i2)
            if (n == 0) stop 'ERROR: member of face kernel not found in cell mask kernel when summing face elements'
            cell(i)%kernel(l)%v(n) = cell(i)%kernel(l)%v(n) + cell(i)%kernel(0)%v(jj)*face(j)%kernel(l)%v(ii2)
          end do
        end do

        cell(i)%kernel(l)%v = cell(i)%kernel(l)%v*h_kernel ! rescaling here so that all derivative kernels can be unscaled later

!---------------
! otherwise use mls or optimisation method to create new kernels
      else 

        if (l == 0.or.l == 1.or.l == 4) then
! set the (minimum) default separations
          if (cell(i)%type == 2) then
            minimum_separation = minimum_boundary_separation
            maximum_separation = maximum_boundary_separation
          else 
            minimum_separation = minimum_domain_separation
            maximum_separation = maximum_domain_separation
          end if
        end if

! calculate separation arrays and for the derivatives create the mask
        if (l == 1) then
          call resize_integer_array(keep_data=.false.,array=separation_index,new_size=1)
          separation_index(1) = 1 ! last index in kernel%ijk that has a cell with separation 1
          call resize_integer_array(keep_data=.false.,array=separation_array,new_size=1)
          separation_array = 1
          call resize_integer_array(keep_data=.false.,array=cell(i)%kernel(l)%ijk,new_size=1)
          cell(i)%kernel(l)%ijk = i
! expand the separation_arrays to include all cells up to and including the maximum_separation
          call expand_kernel_mask(iarray=cell(i)%icell,maximum_separation=maximum_separation,imask=cell(i)%kernel(1)%ijk, &
            separation_index=separation_index,separation_array=separation_array)
        else if (l == 2.or.l == 3) then
          call copy_integer_array(original=cell(i)%kernel(1)%ijk,copy=cell(i)%kernel(l)%ijk)
! NB: separation_index and separation_array just get reused from the l = 1 case
        else ! (l == 0.or.l == 4)
          call resize_integer_array(keep_data=.false.,array=separation_index,new_size=1)
          separation_index(1) = ubound(cell(i)%kernel(l)%ijk,1) ! last index in kernel%ijk that has a cell with separation 1
          call resize_integer_array(keep_data=.false.,array=separation_array,new_size=separation_index(1))
          separation_array = 1
        end if

! allocate and zero value arrays
        if (l /= 0.and.l /= 4) then
          allocate(cell(i)%kernel(l)%v(ubound(cell(i)%kernel(l)%ijk,1)))
          cell(i)%kernel(l)%v = 0.d0
        end if

! create r, scale it with h_kernel, and then convert to a consistent basis
        if (l == 0.or.l == 1.or.l == 4) then
!         if (allocated(r)) deallocate(r)
          allocate(r(totaldimensions,ubound(cell(i)%kernel(l)%ijk,1)))
          do ii = 1, ubound(cell(i)%kernel(l)%ijk,1)
            if (cell(i)%kernel(l)%centring.eq.'cell') then
              r(:,ii) = cell(cell(i)%kernel(l)%ijk(ii))%x - cell(i)%x
            else if (cell(i)%kernel(l)%centring.eq.'face') then
              r(:,ii) = face(cell(i)%kernel(l)%ijk(ii))%x - cell(i)%x
            else if (cell(i)%kernel(l)%centring.eq.'node') then
              r(:,ii) = node(cell(i)%kernel(l)%ijk(ii))%x - cell(i)%x
            else
              call error_stop('ERROR: problem in setup_kernels with cell kernel centring')
            end if
          end do
          r = r/h_kernel
! construct norm if required and also convert to same basis
          if (l == 1) then ! only has to be calculated once
!           if (allocated(norm)) deallocate(norm)
            allocate(norm(totaldimensions,totaldimensions))
            norm = 0.d0
            norm(1,1) = 1.d0
            norm(2,2) = 1.d0
            norm(3,3) = 1.d0
            call construct_orthogonal_basis('cell',r=r,norm=norm,error=error)
          else
            call construct_orthogonal_basis('cell',r=r,error=error)
          end if
          if (error) call error_stop('ERROR: unable to construct orthogonal basis vectors for cell kernel')

! calculate polynomial basis pp tensor from list of r vectors
          if (l == 1) then
            call construct_polynomial_basis_tensor(r,polynomial_order,pp,error)
          else if (l == 0.or.l == 4) then
! the order of the averaging kernels is limited to be <= 1 as there are only a limited number of elements in the kernel mask
            call construct_polynomial_basis_tensor(r,local_polynomial_order=min(1,polynomial_order),pp=pp,error=error)
          end if
          if (error) call error_stop('unable to construct pp basis tensor for cell kernel')

! check minw, enlarging the minimum_separation if required
          if (trim(kernel_method) == 'mls'.and.mls_check_minw) &
            call check_mask_minw(pp,separation_index,minimum_separation,minw)

        end if

! use mls method to construct kernels
        if (l == 0.or.l == 4) then ! average from surrounding faces (l=0) or nodes (l=4)
          if (trim(kernel_method) == 'mls') then
            call mls_kernel(centring='cell',ijk=i,l_kernel=l,rr=r,pp=pp,kernel=cell(i)%kernel(l)%v, &
              local_polynomial_order=min(1,polynomial_order),minimum_separation=minimum_separation, &
              separation_array=separation_array,separation_index=separation_index,error=error)
          else
            call optimisation_kernel(centring='cell',ijk=i,l_kernel=l,rr=r,pp=pp,kernel=cell(i)%kernel(l)%v, &
              local_polynomial_order=min(1,polynomial_order),minimum_separation=minimum_separation, &
              separation_array=separation_array,separation_index=separation_index,error=error)
          end if
        else ! derivatives
          if (vector_magnitude(norm(:,l)) < 1.d-10) then
            if (debug) then
              write(83,'(a)') 'norm component when expressed in basis is zero: skipping mls kernel construction'
              write(83,*) 'l = ',l,': norm(:,l) = ',norm(:,l),': vector_magnitude(norm(:,l)) = ',vector_magnitude(norm(:,l)) 
            end if
          else if (trim(kernel_method) == 'mls') then
            call mls_kernel(centring='cell',ijk=i,l_kernel=l,rr=r,norm=norm(:,l),pp=pp,kernel=cell(i)%kernel(l)%v, &
              local_polynomial_order=polynomial_order,minimum_separation=minimum_separation,separation_array=separation_array, &
              separation_index=separation_index,error=error)
          else
            call optimisation_kernel(centring='cell',ijk=i,l_kernel=l,rr=r,norm=norm(:,l),pp=pp,kernel=cell(i)%kernel(l)%v, &
              local_polynomial_order=polynomial_order,minimum_separation=minimum_separation,separation_array=separation_array, &
              separation_index=separation_index,error=error)
          end if
        end if
        if (error) call error_stop('ERROR: problem calculating cell kernel')

      end if
!---------------
! print some debugging info about kernels

      if (debug_sparse) then
        if (allocated(separation_array)) then
          formatline = '(a,i1,a,100(a,f7.3,a,'//trim(indexformat)//',a,i2,a))'
          write(83,fmt=formatline) 'l = ',l,': cent. = '//trim(cell(i)%kernel(l)%centring)//': value(ijk,sep.) =', &
            (' ',cell(i)%kernel(l)%v(ii),' (',cell(i)%kernel(l)%ijk(ii),',',separation_array(ii),')', &
            ii=1,ubound(cell(i)%kernel(l)%ijk,1))
        else
          formatline = '(a,i1,a,100(a,f7.3,a,'//trim(indexformat)//',a))'
          write(83,fmt=formatline) 'l = ',l,': cent. = '//trim(cell(i)%kernel(l)%centring)//': value(ijk) =', &
            (' ',cell(i)%kernel(l)%v(ii),' (',cell(i)%kernel(l)%ijk(ii),')',ii=1,ubound(cell(i)%kernel(l)%ijk,1))
        end if
        if ((l <= 1.or.l == 4).and.trim(kernel_method) == 'mls'.and.mls_check_minw) write(83,*) 'minw = ',minw
      end if

      if (debug) then
        do ii = 1, ubound(cell(i)%kernel(l)%ijk,1)
          i2 = cell(i)%kernel(l)%ijk(ii)
          if (.not.allocated(r)) then
            formatline = '(a,i3,a,'//trim(indexformat)//',a,g9.2)'
            write(83,fmt=formatline) 'ii = ',ii,': i2 = ',i2,': v = ',cell(i)%kernel(l)%v(ii)
          else if (allocated(r).and..not.allocated(norm)) then
            formatline = '(a,i3,a,'//trim(indexformat)//',a'//repeat(',1x,f6.2',ubound(r,1))//',2(a,g9.2))'
            write(83,fmt=formatline) 'ii = ',ii,': i2 = ',i2,': r =',r(:,ii), &
              ': rmag = ',vector_magnitude(r(:,ii)),': v = ',cell(i)%kernel(l)%v(ii)
          else
            formatline = '(a,i3,a,'//trim(indexformat)//',a'//repeat(',1x,f6.2',ubound(r,1))//',3(a,g9.2))'
            write(83,fmt=formatline) 'ii = ',ii,': i2 = ',i2,': r =',r(:,ii), &
              ': rmag = ',vector_magnitude(r(:,ii)),': r.norm = ',dot_product(r(:,ii),norm(:,l)),': v = ',cell(i)%kernel(l)%v(ii)
          end if
        end do
      end if

      if (l >= 1.and. l <= 3) cell(i)%kernel(l)%v=cell(i)%kernel(l)%v/h_kernel

      if (allocated(r).and.(l == 0.or.l >= 3)) deallocate(r)
      if (allocated(norm).and.(l == 0.or.l >= 3)) deallocate(norm)
      if (allocated(separation_array).and.(l == 0.or.l >= 3)) deallocate(separation_array)
      if (allocated(separation_index).and.(l == 0.or.l >= 3)) deallocate(separation_index)

    end do cell_direction_loop

    if (debug_sparse) then
      formatline = '(a,'//trim(indexformat)//',a,i1)'
      write(83,fmt=formatline) 'END separation_loop: all kernels calculated for: i = ',i,'i: type = ',cell(i)%type
    end if

  end do

else
!----------------------
! uber simple masks based on number of elements
! only 0 (average from surrounding faces) and 4 (average from surrounding nodes) defined
! average from face kernel (l=0) which is inverse of number of faces in kernel
! average from node kernel (l=4) which is inverse of number of nodes in kernel
  do i=1,itotal
    do jj = 1, ubound(cell(i)%kernel(0)%ijk,1)
      cell(i)%kernel(0)%v(jj) = 1.d0/dble(ubound(cell(i)%kernel(0)%ijk,1))
    end do
    do kk = 1, ubound(cell(i)%kernel(4)%ijk,1)
      cell(i)%kernel(4)%v(kk) = 1.d0/dble(ubound(cell(i)%kernel(4)%ijk,1))
    end do
  end do
end if

if (allocated(r)) deallocate(r)
if (allocated(norm)) deallocate(norm)
if (allocated(separation_index)) deallocate(separation_index)
if (allocated(separation_array)) deallocate(separation_array)

!------------------------------------------
! perform efficiency measures

if (remove_small_elements) then

  if (debug_sparse) write(*,*) 'removing small kernel elements'

! count and remove zeroish kernel functions

!------------
! face

  do j = 1, jtotal
    do l = lbound(face(j)%kernel,1), ubound(face(j)%kernel,1)
      h_kernel = 1.d0 ! kernels should have a size of order 1/h_kernel
      if (l >= 1) h_kernel = face(j)%dx
      n_elements = 0
      do n = 1, ubound(face(j)%kernel(l)%ijk,1)
        if (abs(face(j)%kernel(l)%v(n)) >= small_element_minimum/h_kernel) then
          n_elements = n_elements + 1
          face(j)%kernel(l)%v(n_elements) = face(j)%kernel(l)%v(n)
          face(j)%kernel(l)%ijk(n_elements) = face(j)%kernel(l)%ijk(n)
        end if
      end do
! resize kernel arrays (keeping data), but keep 1 zero element to avoid problems with general handling later
!     call resize_integer_array(array=face(j)%kernel(l)%ijk,new_size=max(n_elements,1))
!     call resize_double_precision_array(array=face(j)%kernel(l)%v,new_size=max(n_elements,1))
!     if (n_elements == 0) then
!       face(j)%kernel(l)%ijk(1) = 0
!       face(j)%kernel(l)%v(1) = 0.d0
!     end if
! resize kernel arrays (keeping data), but keep 1 zero element to avoid problems with general handling later
      if (n_elements == 0) then
        deallocate(face(j)%kernel(l)%ijk)
        allocate(face(j)%kernel(l)%ijk(0))
        deallocate(face(j)%kernel(l)%v)
        allocate(face(j)%kernel(l)%v(0))
      else
        call resize_integer_array(keep_data=.true.,array=face(j)%kernel(l)%ijk,new_size=n_elements)
        call resize_double_precision_array(keep_data=.true.,array=face(j)%kernel(l)%v,new_size=n_elements)
      end if
    end do
  end do
      
!------------
! cell

  do i = 1, itotal
    do l = lbound(cell(i)%kernel,1), ubound(cell(i)%kernel,1)
      n_elements = 0
      h_kernel = 1.d0 ! kernels should have a size of order 1/h_kernel
      if (l >= 1.and. l <= 3) h_kernel = cell(i)%vol**(1.d0/dble(cell(i)%dimensions))
      do n = 1, ubound(cell(i)%kernel(l)%ijk,1)
        if (abs(cell(i)%kernel(l)%v(n)) >= small_element_minimum/h_kernel) then
          n_elements = n_elements + 1
          cell(i)%kernel(l)%v(n_elements) = cell(i)%kernel(l)%v(n)
          cell(i)%kernel(l)%ijk(n_elements) = cell(i)%kernel(l)%ijk(n)
        end if
      end do
! resize kernel arrays (keeping data), but keep 1 zero element to avoid problems with general handling later
!     call resize_array(array=cell(i)%kernel(l)%ijk,new_size=max(n_elements,1))
!     call resize_array(array=cell(i)%kernel(l)%v,new_size=max(n_elements,1))
!     if (n_elements == 0) then
!       cell(i)%kernel(l)%ijk(1) = 0
!       cell(i)%kernel(l)%v(1) = 0.d0
!     end if
      if (n_elements == 0) then
        deallocate(cell(i)%kernel(l)%ijk)
        allocate(cell(i)%kernel(l)%ijk(0))
        deallocate(cell(i)%kernel(l)%v)
        allocate(cell(i)%kernel(l)%v(0))
      else
        call resize_integer_array(keep_data=.true.,array=cell(i)%kernel(l)%ijk,new_size=n_elements)
        call resize_double_precision_array(keep_data=.true.,array=cell(i)%kernel(l)%v,new_size=n_elements)
      end if
    end do
  end do

end if

!---------------------------
! print warnings of any negative zeroth order kernels

if (debug_sparse) write(*,*) 'checking kernel consistencies and printing statistics'

formatline = '(a,'//trim(indexformat)//',a,g12.6,a,'//trim(indexformat)//')'
n = 0
min_value = 0.d0
min_location = 0
do i = 1, itotal
  value = minval(cell(i)%kernel(0)%v)
  if (value < 0.d0) then
    n = n + 1
    if (value < min_value) then
      min_value = value
      min_location = i
    end if
  end if
end do
if (n > 0) write(fwarn,fmt=formatline) &
  'WARNING: ',n,' negative cell averaging kernel elements detected: minimum kernel element = ',min_value,': at i = ',min_location

n = 0
min_value = 0.d0
min_location = 0
do j = 1, jtotal
  value = minval(face(j)%kernel(0)%v)
  if (value < 0.d0) then
    n = n + 1
    if (value < min_value) then
      min_value = value
      min_location = j
    end if
  end if
end do
if (n > 0) write(fwarn,fmt=formatline) &
  'WARNING: ',n,' negative face averaging kernel elements detected: minimum kernel element = ',min_value,': at j = ',min_location

! also run a check on the facenorm derivatives on the boundaries, looking for positive elements
n = 0
max_value = 0.d0
max_location = 0
do j = 1, jtotal
  if (face(j)%type /= 2) cycle
  value = maxval(face(j)%kernel(4)%v(3:ubound(face(j)%kernel(4)%v,1)))
  if (value > 0.d0) then
    n = n + 1
    if (value > max_value) then
      max_value = value
      max_location = j
    end if
  end if
end do
if (n > 0) write(fwarn,fmt=formatline) &
  'WARNING: ',n,' positive face normal derivative kernel elements detected: maximum kernel element = ',max_value,': at j = ',max_location

!---------------------------
! print warnings if not enough kernels have been constructed for averaging and in each direction

do i = 1, itotal
  formatline = '(a,'//trim(indexformat)//')'
  if (maxval(abs(cell(i)%kernel(0)%v)) < 1.d-10) write(fwarn,fmt=formatline) &
    'missing averaging kernel for cell: i = ',i
  formatline = '(a,'//trim(indexformat)//',a,i1,a,i1)'
  n = 0
  do l = 1, totaldimensions
    if (maxval(abs(cell(i)%kernel(l)%v)) >= 1.d-10) n = n + 1
  end do
  if (n < cell(i)%dimensions) write(fwarn,fmt=formatline) &
    'not enough derivative kernels have been constructed for cell: i = ',i, &
    ': number constructed = ',n,': cell dimensions = ',cell(i)%dimensions
end do
  
do j = 1, jtotal
  formatline = '(a,'//trim(indexformat)//')'
  if (maxval(abs(face(j)%kernel(0)%v)) < 1.d-10) write(fwarn,fmt=formatline) &
    'missing averaging kernel for face: j = ',j
  formatline = '(a,'//trim(indexformat)//',a,i1,a,i1)'
  n = 0
  do l = 1, 2*totaldimensions
    if (maxval(abs(face(j)%kernel(l)%v)) >= 1.d-10) n = n + 1
  end do
  if (n < 2*(face(j)%dimensions+1)) write(fwarn,fmt=formatline) &
    'not enough derivative kernels have been constructed for face: j = ',j, &
    ': number constructed = ',n,': face dimensions = ',face(j)%dimensions
end do
  
!---------------------------
! check and print warnings if kernels are inconsistent

! cell kernels
formatline = '(a,'//trim(indexformat)//',a,i1,a,i1,a,'//trim(floatformat)//')'
allocate(kernel_error(1:4))
do i = 1, itotal
  do l = 0, 4
    h_kernel = 1.d0 ! kernels should have a size of order 1/h_kernel
    if (l >= 1.and. l <= 3) h_kernel = cell(i)%vol**(1.d0/dble(cell(i)%dimensions))
    kernel_error = 0.d0
    if (maxval(abs(cell(i)%kernel(l)%v)) > 1.d-10) then
      if (l == 0) then ! averaging kernel l = 0, face centred
        kernel_error(1) = sum(cell(i)%kernel(l)%v) - 1.d0
        do n = 1, 3
          do ij = 1, ubound(cell(i)%kernel(l)%ijk,1)
            kernel_error(n+1) = kernel_error(n+1) + cell(i)%kernel(l)%v(ij)*face(cell(i)%kernel(l)%ijk(ij))%x(n)
          end do
          kernel_error(n+1) = kernel_error(n+1) - cell(i)%x(n)
        end do
      else if (l == 4) then ! averaging kernel l = 4, node centred
        kernel_error(1) = sum(cell(i)%kernel(l)%v) - 1.d0
        do n = 1, 3
          do ij = 1, ubound(cell(i)%kernel(l)%ijk,1)
            kernel_error(n+1) = kernel_error(n+1) + cell(i)%kernel(l)%v(ij)*node(cell(i)%kernel(l)%ijk(ij))%x(n)
          end do
          kernel_error(n+1) = kernel_error(n+1) - cell(i)%x(n)
        end do
      else ! derivative kernels l = 1,2,3, cell centred
        kernel_error(1) = sum(cell(i)%kernel(l)%v)
        do n = 1, 3
          do ij = 1, ubound(cell(i)%kernel(l)%ijk,1)
            kernel_error(n+1) = kernel_error(n+1) + cell(i)%kernel(l)%v(ij)*cell(cell(i)%kernel(l)%ijk(ij))%x(n)
          end do
          if (n == l) kernel_error(n+1) = kernel_error(n+1) - 1.d0
        end do
      end if
      do n = 1, 4
        if (abs(kernel_error(n)) > 1.d-8/h_kernel) write(fwarn,fmt=formatline) 'ERROR: for cell(',i,')%kernel(',l,') equation ', &
          n,' is inconsistent: kernel_error(n) = ',kernel_error(n)
      end do
    end if
  end do
end do
deallocate(kernel_error)
  
! face kernels
allocate(kernel_error(1:4))
do j = 1, jtotal
  do l = 0, 6
    h_kernel = 1.d0 ! kernels should have a size of order 1/h_kernel
    if (l >= 1) h_kernel = face(j)%dx
    kernel_error = 0.d0
    if (maxval(abs(face(j)%kernel(l)%v)) > 1.d-10) then
      if (l == 0) then ! averaging kernel l = 0, cell centred
        kernel_error(1) = sum(face(j)%kernel(l)%v) - 1.d0
        do n = 1, 3
          do ij = 1, ubound(face(j)%kernel(l)%ijk,1)
            kernel_error(n+1) = kernel_error(n+1) + face(j)%kernel(l)%v(ij)*cell(face(j)%kernel(l)%ijk(ij))%x(n)
          end do
          kernel_error(n+1) = kernel_error(n+1) - face(j)%x(n)
        end do
      else ! derivative kernels l = 1,2,3,4,5,6, cell centred
        kernel_error(1) = sum(face(j)%kernel(l)%v)
        do n = 1, 3
          do ij = 1, ubound(face(j)%kernel(l)%ijk,1)
            kernel_error(n+1) = kernel_error(n+1) + face(j)%kernel(l)%v(ij)*cell(face(j)%kernel(l)%ijk(ij))%x(n)
          end do
          if (l < 4) then
            if (n == l) kernel_error(n+1) = kernel_error(n+1) - 1.d0
          else
            kernel_error(n+1) = kernel_error(n+1) - face(j)%norm(n,l-3)
          end if
        end do
      end if
      do n = 1, 4
        if (abs(kernel_error(n)) > 1.d-8/h_kernel) write(fwarn,fmt=formatline) 'ERROR: for face(',j,')%kernel(',l,') equation ', &
          n,' is inconsistent: kernel_error(n) = ',kernel_error(n)
      end do
    end if
  end do
end do
deallocate(kernel_error)
  
!------------------------------------------
! calculate some kernel statistics

formatline = '(a,i1,a,'//trim(indexformat)//',a,f5.1)'

n_elements = 0
do l = 0, 6
  n_boundary_kernels = 0
  n_boundary_elements = 0
  n_domain_kernels = 0
  n_domain_elements = 0
  do j = 1, jtotal
    if (face(j)%type == 2) then
      n_boundary_kernels = n_boundary_kernels + 1
      n_boundary_elements = n_boundary_elements + ubound(face(j)%kernel(l)%ijk,1)
    else
      n_domain_kernels = n_domain_kernels + 1
      n_domain_elements = n_domain_elements + ubound(face(j)%kernel(l)%ijk,1)
    end if
    n_elements = n_elements + ubound(face(j)%kernel(l)%ijk,1)
  end do
  write(fwarn,fmt=formatline) 'INFO: for face kernel ',l, &
    ': number of face domain element kernels = ',n_domain_kernels, &
    ': average number of elements = ',dble(n_domain_elements)/dble(max(n_domain_kernels,1))
  write(fwarn,fmt=formatline) 'INFO: for face kernel ',l, &
    ': number of face boundary element kernels = ',n_boundary_kernels, &
    ': average number of elements = ',dble(n_boundary_elements)/dble(max(n_boundary_kernels,1))
end do

do l = 0, 4
  n_boundary_kernels = 0
  n_boundary_elements = 0
  n_domain_kernels = 0
  n_domain_elements = 0
  do i = 1, itotal
    if (cell(i)%type == 2) then
      n_boundary_kernels = n_boundary_kernels + 1
      n_boundary_elements = n_boundary_elements + ubound(cell(i)%kernel(l)%ijk,1)
    else
      n_domain_kernels = n_domain_kernels + 1
      n_domain_elements = n_domain_elements + ubound(cell(i)%kernel(l)%ijk,1)
    end if
    n_elements = n_elements + ubound(cell(i)%kernel(l)%ijk,1)
  end do
  write(fwarn,fmt=formatline) 'INFO: for cell kernel ',l, &
    ': number of cell domain element kernels = ',n_domain_kernels, &
    ': average number of elements = ',dble(n_domain_elements)/dble(max(n_domain_kernels,1))
  write(fwarn,fmt=formatline) 'INFO: for cell kernel ',l, &
    ': number of cell boundary element kernels = ',n_boundary_kernels, &
    ': average number of elements = ',dble(n_boundary_elements)/dble(max(n_boundary_kernels,1))
end do

!if (debug) write(*,'(a,i8)') 'INFO: total number of kernel elements = ',n_kernels
write(*,'(a,i8)') 'INFO: total number of kernel elements = ',n_elements
write(fwarn,'(a,i8)') 'INFO: total number of kernel elements = ',n_elements

close(fwarn)

!------------------------------------------
! write out kernel details if requested

if (kernel_details_file) then
  if (debug_sparse) write(*,*) 'writing kernel details to kernel_details.txt file'

  filename = "output/kernel_details.txt"
  open(fdetail,file=trim(filename),status='replace',iostat=ierror)
  if (ierror /= 0) then
    write(*,*) 'ERROR: problem opening file ',filename
    stop
  end if

  write(fdetail,'(a)') '-----------------------------------------------------------------------'
  write(fdetail,'(a)') 'CELL KERNEL DETAILS:'
  do i = 1, itotal
    do l = lbound(cell(i)%kernel,1), ubound(cell(i)%kernel,1)
      formatline = '(a,'//trim(indexformat)//',a,i1,a,i1,a,i1,a,a,a'//repeat(',a,'//trim(floatformat)//',a,'// &
        trim(indexformat)//',a',ubound(cell(i)%kernel(l)%ijk,1))//')'
      write(fdetail,fmt=formatline) ' i = ',i,'i: type = ',cell(i)%type,': dimension = ',cell(i)%dimensions,': l = ',l, &
        ': centring = ',cell(i)%kernel(l)%centring,': kernel v(ijk) =', &
        (' ',cell(i)%kernel(l)%v(n),'(',cell(i)%kernel(l)%ijk(n),')',n=1,ubound(cell(i)%kernel(l)%ijk,1))
    end do
  end do
  write(fdetail,'(a)') '-----------------------------------------------------------------------'
  write(fdetail,'(a)') 'FACE KERNEL DETAILS:'
  do j = 1, jtotal
    do l = lbound(face(j)%kernel,1), ubound(face(j)%kernel,1)
      formatline = '(a,'//trim(indexformat)//',a,i1,a,i1,a,i1,a,a,a'//repeat(',a,'//trim(floatformat)//',a,'// &
        trim(indexformat)//',a',ubound(face(j)%kernel(l)%ijk,1))//')'
      write(fdetail,fmt=formatline) ' j = ',j,'j: type = ',face(j)%type,': dimension = ',face(j)%dimensions,': l = ',l, &
        ': centring = ',face(j)%kernel(l)%centring,': kernel v(ijk) =', &
        (' ',face(j)%kernel(l)%v(n),'(',face(j)%kernel(l)%ijk(n),')',n=1,ubound(face(j)%kernel(l)%ijk,1))
    end do
  end do

  close(fdetail)
end if

!------------------------------------------

if (debug_sparse) write(*,'(a/80(1h-))') 'subroutine setup_kernels'

end subroutine setup_kernels

!-----------------------------------------------------------------

subroutine construct_orthogonal_basis(centring,r,norm,error)

! here we find orthogonal basis vectors for the space spanned by r, and put those vectors
!  and the optional norm array in terms of this new basis
! r and norm get resized within this routine

use general_module
use lapack_module
!!!!use numerical_recipes_module
character(len=*) :: centring
integer :: i, j, m, n, n_plus, nchange, basis_dimension, null_dimension, nnorm
double precision, dimension(:,:), allocatable, optional :: norm
double precision, dimension(:,:), allocatable :: r, a
double precision, dimension(:), allocatable :: w
double precision, dimension(:,:), allocatable :: u, v, r_basis
integer, dimension(:), allocatable :: basis_list, null_list
double precision, dimension(totaldimensions,totaldimensions) :: basis_vectors
logical :: error
double precision, parameter :: small_element = 1.d-7 ! had problems when this was 1.d-10 - think that lack of precision in .msh file may be an issue
logical, parameter :: numrec = .false.
logical, parameter :: debug = .false.

if (debug) write(83,'(80(1h+)/a)') 'subroutine construct_orthogonal_basis'

error = .false.
m = totaldimensions ! number of rows in r
n = ubound(r,2) ! number of vectors in r
if (ubound(r,1) /= totaldimensions) call error_stop('dimensions of r incorrect in construct_orthonormal_basis')
nnorm = 0
if (present(norm)) then
  nnorm = ubound(norm,2) ! number of vectors in norm
  if (ubound(norm,1) /= m) call error_stop('dimensions of norm incorrect in construct_orthonormal_basis')
end if
if (debug) then
  write(83,'(a)') 'centring = '//trim(centring)
  write(83,*) 'j:r'
  do j = 1, m
    write(83,'(i2,100(a,g9.2))') j,':',(r(j,i),' ',i=1,n)
  end do
  if (present(norm)) then
    write(83,*) 'j:norm'
    do j = 1, m
      write(83,'(i2,100(a,g9.2))') j,':',(norm(j,i),' ',i=1,nnorm)
    end do
  end if
end if

n_plus = max(m,n) ! increase n given to svd so that a minimum of m (3) basis vectors are returned
allocate(a(m,n_plus),u(m,m),w(n_plus),v(n_plus,n_plus))
a = 0.d0
a(:,1:n) = r ! copy the first n columns of r to a
! for svd dimensions should be: a(m,n_plus), u(m,m), w(n_plus), v(n_plus,n_plus)
if (numrec) then
  !!!!call numerical_recipes_singular_value_decomposition(a=a,u=u,w=w,v=v,error=error)
else
  call lapack_singular_value_decomposition(a=a,u=u,w=w,v=v,error=error)
end if
if (error) call error_stop('problem performing singular value decomposition in construct_orthogonal_basis')
if (debug) then
  write(83,*) 'done svd decomposition'
  write(83,*) 'j:u'
  do j = 1, m
    write(83,'(i2,100(a,g9.2))') j,':',(u(j,i),' ',i=1,m)
  end do
  write(83,*) 'j:w:v'
  do j = 1, n_plus
    write(83,'(i2,100(a,g9.2))') j,':',w(j),':',(v(j,i),' ',i=1,n_plus)
  end do
end if

do j = 1, m
  if (abs(w(j)) >= small_element) then
    call push_array(array=basis_list,new_element=j)
! else
  else if (vector_magnitude(u(:,j)) >= small_element ) then ! null_list vectors must be nonzero
    call push_array(array=null_list,new_element=j)
  end if
end do
! either of these could be zero
basis_dimension = allocatable_size(basis_list)
null_dimension = allocatable_size(null_list)

! check that a basis has been found that spans fully three dimensionsal test
if (basis_dimension+null_dimension /= totaldimensions) then
  write(*,*) 'ERROR: singular value decomposition was not able to find three orthonormal basis vectors'
  if (debug) write(83,*) 'ERROR: singular value decomposition was not able to find three orthonormal basis vectors'
  error = .true.
end if

! now assemble the basis vectors matrix
basis_vectors = 0.d0
do i = 1, basis_dimension
  basis_vectors(:,i) = u(:,basis_list(i))
end do
do i = basis_dimension+1, basis_dimension+null_dimension
  basis_vectors(:,i) = u(:,null_list(i-basis_dimension))
end do
  
if (debug) then
  write(83,*) 'basis_dimension = ',basis_dimension,': basis_list = ',basis_list
  write(83,*) 'null_dimension = ',null_dimension,': null_list = ',null_list
  write(83,*) 'j:basis_vectors'
  do j = 1, totaldimensions
    write(83,'(i2,100(a,g9.2))') j,':',(basis_vectors(j,i),' ',i=1,totaldimensions)
  end do
end if

! if the face vectors were brought in, use these as a basis instead
if (trim(centring) == 'face'.and.present(norm)) then
  if (nnorm /= 2*totaldimensions) call error_stop('face norm has wrong number of dimensions in construct_orthogonal_basis')
! need to check that the space spanned by the non-null svd basis vectors is the same as that spanned
!  by the same number of norm vectors
! only have to worry about 2 and 1 dimensional spaces - in the former check that the nullspace is equivalent, in the latter check the
!  the basis space is equivalent
  if (basis_dimension == 2.and.abs(sqrt(abs(dot_product(basis_vectors(:,3),norm(:,6))))-1.d0) > small_element) then
    call error_stop('space spanned by svd basis vectors and that spanned by face norm vectors is not the same in '// &
      'construct_orthogonal_basis')
  else if (basis_dimension == 1.and.abs(sqrt(abs(dot_product(basis_vectors(:,1),norm(:,4))))-1.d0) > small_element) then
    call error_stop('space spanned by svd basis vectors and that spanned by face norm vectors is not the same in '// &
      'construct_orthogonal_basis')
  end if
! TODO: maybe there are smarter solutions here
! TODO: if there is an error could leave basis as svd one - no, optimise kernels requires this to work
  basis_vectors(:,1:3) = norm(:,4:6)
else if (trim(centring) == 'cell'.and.present(norm)) then
  if (nnorm /= totaldimensions) call error_stop('cell norm has wrong number of dimensions in construct_orthogonal_basis')
  if (basis_dimension == 2) then
    nchange = maxloc(basis_vectors(:,3),dim=1)
    if (abs(sqrt(abs(dot_product(basis_vectors(:,3),norm(:,nchange))))-1.d0) > small_element) &
      call error_stop('space spanned by svd basis vectors is not representable by coordinate direction vectors in '// &
        'construct_orthogonal_basis')
! set basis vectors to the norm ones, swapping the null_direction one out to the back
    basis_vectors(:,1:3) = norm(:,1:3)
    basis_vectors(:,3) = norm(:,nchange)
    basis_vectors(:,nchange) = norm(:,3)
  else if (basis_dimension == 1) then
    nchange = maxloc(basis_vectors(:,1),dim=1)
    if (abs(sqrt(abs(dot_product(basis_vectors(:,1),norm(:,nchange))))-1.d0) > small_element) &
      call error_stop('space spanned by svd basis vectors is not representable by coordinate direction vectors in '// &
        'construct_orthogonal_basis')
! set basis vectors to the norm ones, swapping the null_direction one out to the back
    basis_vectors(:,1:3) = norm(:,1:3)
    basis_vectors(:,1) = norm(:,nchange)
    basis_vectors(:,nchange) = norm(:,1)
  else
    basis_vectors(:,1:3) = norm(:,1:3)
  end if
! put option here to be able to continue with mls kernels??
end if

if (debug) then
  write(83,*) 'j:basis_vectors after attempting to use the normal vectors'
  do j = 1, totaldimensions
    write(83,'(i2,100(a,g9.2))') j,':',(basis_vectors(j,i),' ',i=1,totaldimensions)
  end do
end if

deallocate(a,u,v,w)
if (allocated(null_list)) deallocate(null_list)
if (allocated(basis_list)) deallocate(basis_list)
if (error) return

! now to express all the r and norm vectors in terms of the new basis
allocate(r_basis(totaldimensions,n+nnorm))
r_basis(:,1:n) = r(:,1:n)
if (present(norm)) r_basis(:,n+1:n+nnorm) = norm(:,1:nnorm)
if (.true.) then
! noting that the basis_vectors matrix is orthogonal, its inverse is equal to its transpose
  r_basis = matmul(transpose(basis_vectors),r_basis)
else if (numrec) then
  !!!!call numerical_recipes_linear_solver(a=basis_vectors,b=r_basis,error=error)
else
  call lapack_linear_solver(a=basis_vectors,b=r_basis,error=error)
end if
if (error) then
  if (debug) write(83,*) 'ERROR: problem when converting vectors to new basis in construct_orthogonal_basis'
  deallocate(r_basis)
  return
end if

if (debug) then
  write(83,*) 'j:r_basis'
  do j = 1, totaldimensions
    write(83,'(i2,100(a,g9.2))') j,':',(r_basis(j,i),' ',i=1,n+nnorm)
  end do
end if

! check that basis_dimension is consistent with r_basis
do j = basis_dimension+1, totaldimensions
  do i = 1, n
    if (abs(r_basis(j,i)) >= small_element) call error_stop('ERROR: non-zero r_basis r component beyond basis_dimension')
  end do
  do i = n+1, n+nnorm
    if (abs(r_basis(j,i)) >= small_element) then
      if (debug) write(83,'(a,i1,a)') &
        'WARNING: normal vector ',i-n,' cannot be expressed in the r space basis vectors'
      r_basis(:,i) = 0.d0 ! zero all elements of this normal as a marker to calling routine
    end if
  end do
end do

! reshape norm and r to return
! keep number of rows of r greater than zero
deallocate(r)
allocate(r(max(basis_dimension,1),n))
r(:,:) = r_basis(1:max(basis_dimension,1),1:n)
if (present(norm)) then
  deallocate(norm)
  allocate(norm(max(basis_dimension,1),nnorm))
  norm(:,:) = r_basis(1:max(basis_dimension,1),n+1:n+nnorm)
end if
deallocate(r_basis)

if (debug) then
  write(83,*) 'j:r'
  do j = 1, max(basis_dimension,1)
    write(83,'(i2,100(a,g9.2))') j,':',(r(j,i),' ',i=1,n)
  end do
  if (present(norm)) then
    write(83,*) 'j:norm'
    do j = 1, max(basis_dimension,1)
      write(83,'(i2,100(a,g9.2))') j,':',(norm(j,i),' ',i=1,nnorm)
    end do
  end if
end if

if (debug) write(83,'(a/80(1h-))') 'subroutine construct_orthogonal_basis'

end subroutine construct_orthogonal_basis

!-----------------------------------------------------------------

subroutine mls_kernel(centring,ijk,l_kernel,rr,norm,pp,kernel,local_polynomial_order,minimum_separation,separation_array, &
  separation_index,error)

! here we calculate the zeroth and first order kernels using a mls method, adapted from
!  papers of L. Cueto-Felgueroso, 2006ish

use general_module
use lapack_module
!!!!use numerical_recipes_module

character(len=*), intent(in) :: centring ! face|cell|node|none that the kernel is associated with (not the location of the data)
integer, intent(in) :: ijk ! index of element
integer, intent(in) :: l_kernel ! (directional) number of the kernel
double precision, dimension(:,:), allocatable :: rr ! array of surrounding points relative to kernel centre
double precision, dimension(:), optional :: norm ! normal for first order direction
double precision, dimension(:), allocatable :: weight ! weighting function
integer, dimension(:), allocatable :: weight_importance ! weighting function importance (not used here)
integer, dimension(:), allocatable :: separation_array
integer, dimension(:), allocatable :: separation_index
double precision, dimension(:), allocatable :: kernel ! kernel to calculate
double precision, dimension(:,:), allocatable :: pp, mmm, mmm_inv, ww, identity
double precision, dimension(:), allocatable :: rr_centre
integer :: local_polynomial_order ! local polynomial_order
integer :: minimum_separation ! largest separation that will be considered in this kernel
integer :: spatial_dimensions, mm, nn, n, i, j, nn_total
logical :: error
double precision :: trace
logical, parameter :: numrec = .false.
logical, parameter :: debug = .false.

if (debug) write(83,'(80(1h+)/a)') 'subroutine mls_kernel'

kernel = 0.d0
error = .false.

! find spatial dimensions
!spatial_dimensions = dimensions
spatial_dimensions = ubound(rr,1)
if (present(norm)) then
  if (spatial_dimensions /= ubound(norm,1)) call error_stop('dimensions of rr and norm do not match in subroutine mls_kernel')
end if
! total number of kernel nodes, including those that may not be used as their separation_index is higher than minimum_separation
nn_total = ubound(rr,2)

! polynomial basis dimension from the pp tensor (number of rows in pp)
mm = ubound(pp,1)
if (ubound(pp,2) /= nn_total) call error_stop('dimensions of pp and rr do not match in mls_kernel')

! calculate the weights
if (trim(centring) == 'face') then
  call calc_weight(centring,ijk,l_kernel,rr,pp,minimum_separation,face(ijk)%kernel(l_kernel)%ijk,separation_array, &
    separation_index,weight,weight_importance)
else
  call calc_weight(centring,ijk,l_kernel,rr,pp,minimum_separation,cell(ijk)%kernel(l_kernel)%ijk,separation_array, &
    separation_index,weight,weight_importance)
end if
nn = ubound(weight,1) ! set number of elements that are in use, which corresponds to the size of weight

if (debug) then
  write(83,*) 'spatial_dimensions,mm,nn,nn_total'
  write(83,*) spatial_dimensions,mm,nn,nn_total
  if (present(norm)) then
    write(83,*) 'derivative: norm = ',norm
  else
    write(83,*) 'average: no norm'
  end if
end if

if (debug) then
  write(83,*) 'j:rr'
  do j = 1, ubound(rr,1)
    write(83,'(i2,100(a,g9.2))') j,':',(rr(j,i),' ',i=1,ubound(rr,2))
  end do
end if

! allocate various
allocate(mmm(mm,mm),mmm_inv(mm,mm),ww(nn,nn),identity(mm,mm))

if (debug) then
  write(83,*) 'j:pp'
  do j = 1, ubound(pp,1)
    write(83,'(i2,100(a,g9.2))') j,':',(pp(j,i),' ',i=1,ubound(pp,2))
  end do
end if

! form ww, which is the tensor diagonal of the weight vector
ww = 0.d0
do n = 1, nn
  ww(n,n) = weight(n)
end do

if (debug) then
  write(83,*) 'j:ww'
  do j = 1, ubound(ww,1)
    write(83,'(i2,100(a,g9.2))') j,':',(ww(j,i),' ',i=1,ubound(ww,2))
  end do
end if

! form mmm
mmm = matmul(matmul(pp(:,1:nn),ww),transpose(pp(:,1:nn)))

! invert mm using lu decompostion by inverting identity matrix column by column
! setup mm_inv as identity matrix
identity = 0.d0
do n = 1, mm
  identity(n,n) = 1.d0
end do
mmm_inv = identity
if (numrec) then
  !!!!call numerical_recipes_linear_solver(a=mmm,b=mmm_inv,error=error,overwritea=.false.)
else
  call lapack_linear_solver(a=mmm,b=mmm_inv,error=error,overwritea=.false.)
end if

! check accuracy of inversion - guarding against numerical singularity
!trace = norm2(matmul(mmm,mmm_inv)-identity) norm2 not supported by gfortran 4.5.1
trace = sqrt(sum((matmul(mmm,mmm_inv)-identity)**2))

if (debug) write(83,*) 'trace = ',trace
if (trace > 1.d-10) error = .true.

if (error) then
  if (debug) write(83,*) 'error occured when trying to invert mmm - returning'
  deallocate(mmm,mmm_inv,ww,identity)
  return
end if

if (debug) then
  write(83,*) 'j:mmm_inv'
  do j = 1, ubound(mmm_inv,1)
    write(83,'(i2,100(a,g9.2))') j,':',(mmm_inv(j,i),' ',i=1,ubound(mmm_inv,2))
  end do
end if

allocate(rr_centre(spatial_dimensions))
rr_centre = 0.d0
if (.not.present(norm)) then
  kernel(1:nn) = matmul(polynomial_basis_vector(rr=rr_centre,local_polynomial_order=local_polynomial_order), &
    matmul(mmm_inv,matmul(pp(:,1:nn),ww)))
else
  kernel(1:nn) = matmul(polynomial_basis_vector(rr=rr_centre,norm=norm,local_polynomial_order=local_polynomial_order), &
    matmul(mmm_inv,matmul(pp(:,1:nn),ww)))
end if

if (debug) then
  write(83,*) 'j:kernel'
  do j = 1, ubound(kernel,1)
    write(83,'(i2,100(a,g9.2))') j,':',kernel(j)
  end do
end if

deallocate(mmm,mmm_inv,ww,rr_centre,identity)

if (debug) write(83,'(a/80(1h-))') 'subroutine mls_kernel'

end subroutine mls_kernel

!-----------------------------------------------------------------

double precision function radial_kernel(rr,derivative_order)

! this is the kernel required for the mls kernel method
! this kernel is 1 at rr=0, 1/2 at rr=1, approaches 0 as rr approaches infinity

integer :: derivative_order
double precision :: rr
double precision :: pi = 4.d0*atan(1.d0)
!integer :: n = 1
!double precision :: a = 4.d0
!integer :: n = 2
!double precision :: a = 1.50d0
integer :: n = 3
! pre-optimisation default
double precision :: a = 1.20d0
!integer :: n = 4
!double precision :: a = 1.d0
!integer :: n = 8
!double precision :: a = 0.65d0
!double precision :: gamma = 1.4d0

radial_kernel = 0.d0

if (derivative_order == 0) then
!  radial_kernel = (cos(pi*rr/(1.d0+rr))+1.d0)/2.d0
!  radial_kernel = ((cos(pi*rr/(rr+1))+1)/2)**n
! if (.true.) then
!   radial_kernel = (cos(pi*a*rr/(a*rr+1))+1)**n/2**n
    radial_kernel = ((cos(pi*a*rr/(a*rr+1))+1)/2.d0)**n
! else
!   radial_kernel = exp(-(rr/gamma)**2)
! end if
else if (derivative_order == 1) then
!  radial_kernel = -pi*sin(pi*rr/(1.d0+rr))/(2.d0*(1.d0+rr)**2)
!  radial_kernel = -pi*n*(cos(pi*rr/(rr+1))+1)**(n-1)*sin(pi*rr/(rr+1))/(2**n*(rr+1)**2)
!  radial_kernel = -pi*a*n*(cos(pi*a*rr/(a*rr+1))+1)**(n-1)*sin(pi*a*rr/(a*rr+1))/(2**n*(a*rr+1)**2)
  radial_kernel = -pi*a*dble(n)*(cos(pi*a*rr/(a*rr+1))+1)**(n-1)*sin(pi*a*rr/(a*rr+1))/(2**n*(a*rr+1)**2)
end if

end function radial_kernel

!-----------------------------------------------------------------

function polynomial_basis_vector(rr,norm,local_polynomial_order)

! here we calculate the polynomial basis vector evaluated at rr
! if norm is present, calculate the derivative (gradient) of this vector in the direction of norm

use general_module
double precision, dimension(:) :: rr ! single location vector of unknown dimensions
double precision, dimension(:), optional :: norm ! normal for first order direction
double precision, dimension(:), allocatable :: polynomial_basis_vector ! result returned
integer :: local_polynomial_order ! local polynomial_order
integer :: spatial_dimensions, vector_dimension, n
logical, parameter :: debug = .false.

if (debug) write(83,'(80(1h+)/a)') 'function polynomial_basis_vector'

if (local_polynomial_order > 2) call error_stop &
  ('polynomial_basis_vector function can handle atmost 2nd order polynomials right now')
spatial_dimensions = ubound(rr,1)
if (present(norm)) then
  if (ubound(norm,1) /= spatial_dimensions) then
    write(*,*) 'norm and rr spatial dimensions don''t match in polynomial_basis_vector'
    stop
  end if
end if
if (local_polynomial_order == 2) then
  if (spatial_dimensions == 1) then
    vector_dimension = 3
  else if (spatial_dimensions == 2) then
    vector_dimension = 6
  else if (spatial_dimensions == 3) then
    vector_dimension = 10
  else
    stop 'ERROR: spatial dimensions incorrect in polynomial_basis_vector'
  end if
else
  vector_dimension = 1 + spatial_dimensions
end if

allocate(polynomial_basis_vector(vector_dimension))

if (.not.present(norm)) then
  polynomial_basis_vector(1) = 1.d0
  polynomial_basis_vector(2:spatial_dimensions+1) = rr
  if (local_polynomial_order == 2) then
    do n = 1, spatial_dimensions
      polynomial_basis_vector(n+spatial_dimensions+1) = rr(n)**2
    end do
    if (spatial_dimensions == 2) then
      polynomial_basis_vector(6) = rr(1)*rr(2)
    else if (spatial_dimensions == 3) then
      polynomial_basis_vector(8) = rr(1)*rr(2)
      polynomial_basis_vector(9) = rr(1)*rr(3)
      polynomial_basis_vector(10) = rr(2)*rr(3)
    end if
  end if
else
  polynomial_basis_vector(1) = 0.d0
  polynomial_basis_vector(2:spatial_dimensions+1) = norm
  if (local_polynomial_order == 2) then
    do n = 1, spatial_dimensions
      polynomial_basis_vector(n+spatial_dimensions+1) = 2.d0*rr(n)*norm(n)
    end do
    if (spatial_dimensions == 2) then
      polynomial_basis_vector(6) = rr(1)*norm(2)+rr(2)*norm(1)
    else if (spatial_dimensions == 3) then
      polynomial_basis_vector(8) = rr(1)*norm(2)+rr(2)*norm(1)
      polynomial_basis_vector(9) = rr(1)*norm(3)+rr(3)*norm(1)
      polynomial_basis_vector(10) = rr(2)*norm(3)+rr(3)*norm(2)
    end if
  end if
end if
  
if (debug) write(83,'(a/80(1h-))') 'function polynomial_basis_vector'

end function polynomial_basis_vector

!-----------------------------------------------------------------

subroutine expand_kernel_mask(iarray,maximum_separation,imask,separation_index,separation_array)

! here we find the kernel mask, and calculate the separation, defined as the number of faces you have to pass through
!  to get to that cell
! iarray is the list of neighbours that could possibly be included in the mask
! the imask array is ordered in increasing separation
! the separation information is stored in array separation_index, which specifies the last index for each separation
! separation_array specifies the separation of each element
! maximum_separation is the maximum separation of cells that should be included in the kernel

use general_module
integer, dimension(:), allocatable :: iarray ! the array of all surrounding cells, passed in as one of the cell or face arrays, not ordered
integer :: maximum_separation ! this is the maximum separation that we will consider
integer, dimension(:), allocatable :: imask ! this is the list of cells that is passed out, with a changed size
integer, dimension(:), allocatable :: separation_index ! this is the list of separation indicies is passed out, with a changed size
integer, dimension(:), allocatable :: separation_array ! this the separation indicies per element passed out, with a changed size
integer :: separation, i, ii, start_index, end_index, i2, ii2
logical, parameter :: debug = .false.

if (debug) write(83,'(80(1h+)/a)') 'subroutine expand_kernel_mask'

separation = ubound(separation_index,1) ! this is the previously highest defined separation

if (ubound(imask,1) /= separation_index(separation)) call error_stop( &
  'imask size does not match with final separation_index in calling arguments to expand_kernel_mask')
if (ubound(separation_array,1) /= separation_index(separation)) call error_stop( &
  'separation_array size does not match with final separation_index in calling arguments to expand_kernel_mask')

if (debug) then
  write(83,*) 'calling parameters: iarray = ',iarray
  write(83,*) 'calling parameters: maximum_separation = ',maximum_separation
  write(83,*) 'calling parameters: imask = ',imask
  write(83,*) 'calling parameters: separation_index = ',separation_index
  write(83,*) 'calling parameters: separation_array = ',separation_array
end if
  
if (separation == 1) then
  end_index = 0
else
  end_index = separation_index(separation-1)
end if

! loop through separations until we reach the requested maximum_separation
separation_loop: do while (separation < maximum_separation.and. &
  (limit_mask_to_icell.and.ubound(iarray,1) > ubound(imask,1).or..not.limit_mask_to_icell))

! define indicies of imask that bound last separation
  start_index = end_index + 1
  end_index = separation_index(separation)
  if (end_index < start_index) call error_stop('indicies out of order in expand_kernel_masks')

  separation = separation + 1
  call resize_integer_array(keep_data=.true.,array=separation_index,change=1)
  separation_index(separation) = separation_index(separation-1)

  if (debug) then
    write(83,'(3(a,i3))') 'separation = ',separation,': start_index = ',start_index,': end_index = ',end_index
    write(83,*) 'previous separation (separation = ',separation-1,') list = ',imask(start_index:end_index)
  end if

! loop through all cells which have the last separation looking for ones that have no separation and which are adjacent
  last_separation_loop: do ii = start_index, end_index
    i = imask(ii)

! loop through neighbours of this cell that share a face element
    neighbour_loop: do ii2 = 2, ubound(cell(i)%jface,1)+1
      i2 = cell(i)%icell(ii2)
      if (location_in_list(array=imask,element=i2) /= 0) cycle neighbour_loop ! cell is already in list
      if (limit_mask_to_icell.and.location_in_list(array=iarray,element=i2) == 0) cycle neighbour_loop ! cell is not in overall iarray list
! if we are here then cell i2 has the current separation
      if (debug) write(83,*) 'found new imask element = ',i2,' having separation = ',separation
      call push_integer_array(array=imask,new_element=i2)
      separation_index(separation) = ubound(imask,1)
    end do neighbour_loop

  end do last_separation_loop

! as cells are separated by common faces, if no cells are found at a given separation level then this is an error
  if (separation_index(separation) == separation_index(separation-1)) call error_stop( &
    'no elements are found at a separation level in expand_kernel_mask, indicating a mesh problem of some sort')
  
  if (debug) write(83,*) 'calculated separation (separation = ',separation,') list = ', &
    imask(separation_index(separation-1)+1:separation_index(separation))
  
end do separation_loop

! update separation_array
! TODO: could keep original data here but not much worth
! TODO: could put it in the main loop but would have to think about array size allocation - probably not worth it efficiency wise
call resize_integer_array(keep_data=.false.,array=separation_array,new_size=separation_index(ubound(separation_index,1)))
end_index = 0
do separation = 1, ubound(separation_index,1)
  start_index = end_index + 1
  end_index = separation_index(separation)
  separation_array(start_index:end_index) = separation
end do

if (debug) then
  write(83,*) 'final parameters: imask = ',imask
  write(83,*) 'final parameters: separation_index = ',separation_index
  write(83,*) 'final parameters: maximum_separation = ',maximum_separation
  write(83,*) 'final parameters: ubound(iarray,1) = ',ubound(iarray,1)
  write(83,*) 'final parameters: ubound(imask,1) = ',ubound(imask,1)
  write(83,*) 'final parameters: separation_array = ',separation_array
end if
  
if (debug) write(83,'(a/80(1h-))') 'subroutine expand_kernel_mask'

end subroutine expand_kernel_mask

!-----------------------------------------------------------------

subroutine calc_weight(centring,ijk,l_kernel,rr,pp,minimum_separation,imask,separation_array,separation_index, &
  weight,weight_importance)

! here we calculate the weight scaling factor for all cells in the kernel
! weight is either a true weighting (mls method) or a trial kernel value (optimisation method)
! also calculate weight_importance which is an integer measure of how much want the kernel to be weighted there
!  (0=high weighting, >0=lower weighting)

use general_module
character(len=*), intent(in) :: centring ! face|cell that the kernel is associated with (not the location of the data)
integer, intent(in) :: ijk ! index of element
integer, intent(in) :: l_kernel ! (directional) number of the kernel
double precision, dimension(:), allocatable :: weight ! the list of weights that is passed out
integer, dimension(:), allocatable :: weight_importance ! this specifies the importance (separation) index of each element
integer, dimension(:), allocatable :: imask ! list of cells that could be in the kernel
integer :: minimum_separation ! only elements with this separation will have non-zero weights for the mls method
integer, dimension(:), allocatable :: separation_index ! these are the final elements having that separation
integer, dimension(:), allocatable :: separation_array ! this is the separation indicies of each element
double precision, dimension(:,:), allocatable :: rr ! array of cell locations, first index is spatial dimension, second is point number
double precision, dimension(:,:), allocatable :: pp ! polynomial basis tensor, arranged as columns corresponding to each point in the kernel
integer :: l_kernel_local, nn_total, nn, separation, i, ii, start_index, end_index, i2, ii2, flux_to_cells
double precision :: weight_sum, rr_mag
logical, parameter :: debug = .false.

if (debug) write(83,'(80(1h+)/a)') 'subroutine calc_weight'

if (debug) then
  write(83,*) 'calling parameters: centring = '//trim(centring)
  write(83,*) 'calling parameters: ijk = ',ijk
  write(83,*) 'calling parameters: l_kernel = ',l_kernel
  write(83,*) 'calling parameters: minimum_separation = ',minimum_separation
  write(83,*) 'calling parameters: imask = ',imask
  write(83,*) 'calling parameters: separation_index = ',separation_index
  write(83,*) 'calling parameters: separation_array = ',separation_array
end if
  
! total number of elements
nn_total = ubound(rr,2)

! bound minimum separation by the maximum separation available
minimum_separation = min(minimum_separation,ubound(separation_index,1))

! calc nn, which is the number of elements that should be included in the kernel calc
! different for mls versus optimisation 
if (trim(kernel_method) == 'optimisation') then
  nn = nn_total ! all cells included at this stage
else
  nn = separation_index(min(minimum_separation,ubound(separation_index,1))) ! only those up to and including minimum separation
end if

! resize weight and zero value
call resize_double_precision_array(keep_data=.false.,array=weight,new_size=nn)
weight = 0.d0

! find weight_importance here
call copy_integer_array(original=separation_array,copy=weight_importance)
weight_importance = weight_importance - 1 ! 0 is the minimum in the weight index
! orientation_dependent_weights means that for face kernels and tangential derivatives, the cells adjacent to the face are given less weighting
if (orientation_dependent_weights.and.((trim(centring)=='face'.and.l_kernel >= 5).or. &
  (trim(centring)=='cell'.and.l_kernel >= 1.and.l_kernel <= 3))) then
  if (ubound(separation_index,1) >= 1) weight_importance(1:separation_index(1)) = 1
  if (ubound(separation_index,1) >= 2) weight_importance(separation_index(1)+1:separation_index(2)) = 0
end if

if (conservative_weighting.and.trim(centring)=='face'.and.(l_kernel == 0.or.l_kernel == 4)) then
! not recommended anymore

  weight(1:2) = 0.5d0 ! initialise first separation weights
  end_index = 0
  do separation = 1, min(ubound(separation_index,1),minimum_separation)-1 ! loop through separations that are to donate weights
    start_index = end_index+1
    end_index = separation_index(separation)
    do ii = start_index, end_index
      i = imask(ii) ! i is a cell that is about to possibly loose some material

! calculate/find the cells that are about to be fluxed-to from cell i
      flux_to_cells = 0
      do ii2 = separation_index(separation)+1, separation_index(separation+1) ! this is all cells having a separation 1 greater
        i2 = imask(ii2)
        if (location_in_list_dummy(array=cell(i)%icell(2:ubound(cell(i)%jface,1)+1),element=i2) /= 0) then ! we have found a neighbour
          flux_to_cells = flux_to_cells + 1 ! if so then it will receive some weight from i
        end if
      end do

      if (debug) write(83,*) 'for cell i = ',i,' having separation = ',separation,' found flux-to_cells = ',flux_to_cells

! if any cells will receive some from i then go ahead and do the fluxing
      if (flux_to_cells > 0) then

        do ii2 = separation_index(separation)+1, separation_index(separation+1) ! this is all cells having a separation 1 greater
          i2 = imask(ii2)
          if (location_in_list_dummy(array=cell(i)%icell(2:ubound(cell(i)%jface,1)+1),element=i2) /= 0) then ! we have found a neighbour
            if (debug) write(83,*) 'adding weight to cell i2 = ',i2,' with current weight = ',weight(ii2)
            weight(ii2) = weight(ii2) + weight(ii)*weight_separation_amount/dble(flux_to_cells)
          end if
        end do

        if (debug) write(83,*) 'removing weight from cell i = ',i,' with current weight = ',weight(ii)
        weight(ii) = weight(ii)*(1.d0 - weight_separation_amount)
      end if

    end do

  end do

else
! these are based on the weight_importance
! if orientation_dependent_weights is set then these are different for each l_kernel direction
! for the optimisation method, these weights are infact the trial kernels

  weight_sum = 0.d0
  do i = 1, nn

! l_kernel_local is the direction of the derivative, if in the range 1->3
    if (trim(centring) == 'face') then
      l_kernel_local = l_kernel - 3
    else 
      l_kernel_local = l_kernel
    end if
      
    rr_mag = vector_magnitude(rr(:,i))
    if (radial_kernel_weighting) then
      if (trim(kernel_method) == 'optimisation'.and.l_kernel_local >= 1.and.l_kernel_local <= 3) then
        weight(i) = (-radial_kernel(rr_mag,1)*rr(l_kernel_local,i)/max(rr_mag,1.d-10))* &
          (weight_separation_amount**weight_importance(i))
        weight_sum = weight_sum + weight(i)*rr(l_kernel_local,i)
      else
        weight(i) = radial_kernel(rr_mag,0)*(weight_separation_amount**weight_importance(i))
        weight_sum = weight_sum + weight(i)
      end if
    else
      weight(i) = weight_separation_amount**weight_importance(i)
      weight_sum = weight_sum + weight(i)
    end if

  end do
! special case normal boundary derivative - overwrite the centred trial kernel (which is zero anyway) to preserve straight sum
  if (radial_kernel_weighting.and.trim(kernel_method) == 'optimisation'.and.trim(centring) == 'face'.and. &
    l_kernel_local == 1) then
    if (face(ijk)%type == 2) weight(2) = -sum(weight)
  end if
  weight = weight/weight_sum

end if

if (debug) then
  write(83,*) 'final parameters: minimum_separation = ',minimum_separation
  write(83,*) 'final parameters: weight = ',weight
  write(83,*) 'final parameters: weight_importance = ',weight_importance
end if
  
if (debug) write(83,'(a/80(1h-))') 'subroutine calc_weight'

end subroutine calc_weight

!-----------------------------------------------------------------

subroutine check_mask_minw(pp,separation_index,minimum_separation,minw)

use general_module
use lapack_module
!!!!use numerical_recipes_module
double precision, dimension(:,:), allocatable :: pp ! array of polynomial basis vectors
integer, dimension(:), allocatable :: separation_index
integer :: minimum_separation
double precision :: minw ! minimum w found representing the singularity of the mask (ref chenoweth09, JCP, 228, p5592)
logical :: error ! whether there are any problems with this mask or not
character(len=1000) :: formatline
logical :: svd_error
double precision, dimension(:,:), allocatable :: bb, u, v
double precision, dimension(:), allocatable :: w
integer :: minimum_separation_in, i, j, mm, nn, nn_tot
logical :: numrec = .false.
logical, parameter :: debug = .false.
logical :: debug_sparse = .false.

if (debug) debug_sparse = .true.

if (debug_sparse) write(83,'(80(1h+)/a)') 'subroutine check_mask_minw'

error = .false.
minw = -1.d0 ! represents not set

mm = ubound(pp,1) ! this is the number of elements in each polynomial vector
nn_tot = ubound(pp,2) ! this is the total number of elements in the kernel
minimum_separation_in = min(minimum_separation,ubound(separation_index,1)) ! limit the minimum_separation to the separations available

if (debug_sparse) write(83,*) 'mm = ',mm,': nn_tot = ',nn_tot,': minimum_separation_in = ',minimum_separation_in

do minimum_separation = minimum_separation_in, ubound(separation_index,1)

  nn = separation_index(minimum_separation) ! find the number of elements in the kernel up to and including the minimum separation

  allocate(bb(nn,mm),u(nn,nn),w(mm),v(mm,mm))

  bb = transpose(pp(:,1:nn))

  if (debug) then
    write(83,*) 'before svd decomposition'
    write(83,*) 'j:bb'
    do j = 1, nn
      write(83,'(i2,100(a,g9.2))') j,':',(bb(j,i),' ',i=1,mm)
    end do
  end if

! for svd dimensions should be: a(m,n), u(m,m), w(n), v(n,n)
  if (numrec) then
    !!!!call numerical_recipes_singular_value_decomposition(a=bb,u=u,w=w,v=v,error=svd_error)
  else
    call lapack_singular_value_decomposition(a=bb,u=u,w=w,v=v,error=svd_error)
  end if
  if (svd_error) call error_stop('problem performing svd in check_kernel_mask')

  if (debug) then
    write(83,*) 'done svd decomposition'
    write(83,*) 'j:u'
    do j = 1, nn
      write(83,'(i2,100(a,g9.2))') j,':',(u(j,i),' ',i=1,mm)
    end do
    write(83,*) 'j:w:v'
    do j = 1, nn
      write(83,'(i2,100(a,g9.2))') j,':',w(j),':',(v(j,i),' ',i=1,nn)
    end do
  end if

  minw = minval(w)
  if (debug_sparse) write(83,*) 'minw = ',minw,': nn = ',nn,': minimum_separation = ',minimum_separation

  deallocate(bb,u,w,v)

  if (minw < mls_minimum_minw) then
    if (nn == nn_tot) then
      if (debug_sparse) then
        write(83,*) 'WARNING: no more elements available within maximum separation to further increase minw'
      end if
      exit
    end if
    if (debug_sparse) write(83,*) 'Adding additional elements to kernel to increase minw'
! temp &&&
!   if (.true.) write(83,*) 'Adding additional elements to kernel to increase minw'
  else
    exit
  end if

end do

if (debug_sparse) write(83,*) 'final minimum_separation = ',minimum_separation,': minw = ',minw

if (debug_sparse) write(83,'(a/80(1h-))') 'subroutine check_mask_minw'

end subroutine check_mask_minw

!-----------------------------------------------------------------

subroutine construct_polynomial_basis_tensor(r,local_polynomial_order,pp,error)

! subroutine to allocate and construct polynomial basis tensor

double precision, dimension(:,:), allocatable :: r
integer :: local_polynomial_order
double precision, dimension(:,:), allocatable :: pp
logical :: error
integer :: nn, mm, j, i, n
logical, parameter :: debug = .false.

if (debug) write(83,'(80(1h-)/a)') 'subroutine construct_polynomial_basis_tensor'

error = .false.
! the r tensor consists of nn columns of m-dimensional location vectors

! number of kernel nodes
nn = ubound(r,2)

! find polynomial basis dimension by sending a test vector to polynomial_basis_vector
mm = ubound(polynomial_basis_vector(rr=r(:,1),local_polynomial_order=local_polynomial_order),1) ! try first order polynomial basis used
if (nn < mm) then
  if (debug) write(83,*) 'not enough points (dim r) in kernel mask for the requested polynomial_order: mm = ',mm,': nn = ',nn
  error = .true.
  return
end if

if (debug) then
  write(83,*) 'mm,nn,local_polynomial_order'
  write(83,*) mm,nn,local_polynomial_order
end if

if (debug) then
  write(83,*) 'j:r'
  do j = 1, ubound(r,1)
    write(83,'(i2,100(a,g9.2))') j,':',(r(j,i),' ',i=1,ubound(r,2))
  end do
end if

! allocate pp if it doesn't have the correct size already
if (allocated(pp)) then
  if (ubound(pp,1) /= mm.or.ubound(pp,2) /= nn) deallocate(pp)
end if
if (.not.allocated(pp)) allocate(pp(mm,nn))

! form pp
do n = 1, nn
  pp(:,n) = polynomial_basis_vector(rr=r(:,n),local_polynomial_order=local_polynomial_order)
end do

if (debug) then
  write(83,*) 'j:pp'
  do j = 1, ubound(pp,1)
    write(83,'(i2,100(a,g9.2))') j,':',(pp(j,i),' ',i=1,ubound(pp,2))
  end do
end if

if (debug) write(83,'(a/80(1h-))') 'subroutine construct_polynomial_basis_tensor'

end subroutine construct_polynomial_basis_tensor

!-----------------------------------------------------------------

subroutine optimisation_kernel(centring,ijk,l_kernel,rr,norm,pp,kernel,local_polynomial_order,minimum_separation, &
  separation_array,separation_index,error)

! here we use an optimisation technique to calculate kernels which have a polynomial basis

use general_module
use lapack_module
!!!!use numerical_recipes_module
character(len=*), intent(in) :: centring ! face|cell|node|none that the kernel is associated with (not the location of the data)
integer, intent(in) :: ijk ! index of element
integer, intent(in) :: l_kernel ! (directional) number of the kernel
double precision, dimension(:,:), allocatable :: rr ! array of surrounding points relative to kernel centre
double precision, dimension(:), optional :: norm ! normal for first order direction
double precision, dimension(:,:), allocatable :: pp ! polynomial basis tensor, arranged as columns corresponding to each point in the kernel
double precision, dimension(:), allocatable :: kernel ! kernel to calculate
integer :: local_polynomial_order ! local polynomial_order
integer :: minimum_separation ! largest separation that will be considered in this kernel
integer, dimension(:), allocatable :: separation_array
integer, dimension(:), allocatable :: separation_index
double precision, dimension(:), allocatable :: weight ! weighting function
integer, dimension(:), allocatable :: weight_importance ! weighting function importance (not used here)
integer :: nn, mm, d, i, j, n_kernel_steps
double precision :: llnorm, l
double precision, dimension(:), allocatable :: y
double precision, dimension(:,:), allocatable :: ll, lll
integer, dimension(:), allocatable :: active
logical :: error, active_change, linear_error
character(len=10) :: phase ! tells the constraints algorithm what phase we are in
double precision, parameter :: llnorm_tol = 1.d-8 ! tolerance indicating linear solver has been successful
integer, parameter :: maximum_constraint_steps = 100
logical, parameter :: numrec = .false.
logical, parameter :: debug = .false.

error = .false.

if (debug) write(83,'(80(1h+)/a)') 'subroutine optimisation_kernel'

! find d, which is the element within the polynomial basis that this kernel refers to
if (present(norm)) then
! a present norm indicates that we are calculating a derivative in the direction corresponding to d
  d = maxloc(norm,dim=1) + 1
! check that the norm direction is in one of the coordinate directions
  if (abs(norm(d-1)-1.d0) > 1.d-10) call error_stop('the optimisation kernel routine cannot (yet) calculate '// &
    'derivatives that are not in the coordinate direction.  You could try: a) orientating your 1D or 2D mesh to be coincident '// &
    'with the coordinates; b) turn on cell_from_face_kernels (in kernel_module.f90); '// &
    'or c) generalise the optimisation routine and fix construct_orthogonal_basis error checks!')
else
  d = 1
end if
! now find some array dimensions
mm = ubound(pp,1) ! this is the number of elements in each polynomial basis vector, which is the number of rows in the pp tensor
if (d > mm) call error_stop('norm in optimisation_kernel inconsistent with polynomial basis')
nn = ubound(pp,2) ! this corresponds to the number of points in the kernel, which is the number of columns in the pp tensor
if (nn /= ubound(kernel,1)) call error_stop('the dimensions of the pp tensor and kernel vector do not match in optimisation_kernel')

! calculate the weights
if (trim(centring) == 'face') then
  call calc_weight(centring,ijk,l_kernel,rr,pp,minimum_separation,face(ijk)%kernel(l_kernel)%ijk,separation_array, &
    separation_index,weight,weight_importance)
else
  call calc_weight(centring,ijk,l_kernel,rr,pp,minimum_separation,cell(ijk)%kernel(l_kernel)%ijk,separation_array, &
    separation_index,weight,weight_importance)
end if

if (debug) then
  write(83,*) 'd = ',d,': mm = ',mm,': nn = ',nn
  write(83,'(a)') 'weight = '
  write(83,'(100(f8.2,a))') (weight(i),' ',i=1,ubound(weight,1))
  write(83,'(a)') 'weight_importance = '
  write(83,'(100(i8,a))') (weight_importance(i),' ',i=1,nn)
  write(83,'(a)') 'pp = '
  do j = 1, ubound(pp,1)
    write(83,'(100(f8.2,a))') (pp(j,i),' ',i=1,ubound(pp,2))
  end do
end if

! allocate arrays for the optimisation procedure
allocate(y(2*nn+mm))
allocate(ll(ubound(y,1),1))
allocate(lll(ubound(y,1),ubound(y,1)))
allocate(active(ubound(y,1)))

y = 1.d0 ! this sets all kernels and all lagrangian functions to 1

phase = 'first'
linear_error = .false.
call optimisation_kernel_constraints(pp,y,weight,separation_array,minimum_separation,d,nn,mm,active,active_change, &
  phase,linear_error,error)

if (error) call error_stop('unresolvable error during '//trim(phase)//' phase of optimisation_kernel_constraints')

call optimisation_kernel_update(pp,y,weight,weight_importance,separation_array,d,nn,mm,active,ll(:,1),lll,llnorm,l)

if (debug) write(83,'(2(a,g13.6))') 'INFO: entering constraint loop: initial llnorm = ',llnorm
n_kernel_steps = 0

constraint_loop: do

  n_kernel_steps = n_kernel_steps + 1

  lll = -lll

  if (numrec) then
    !!!!call numerical_recipes_linear_solver(a=lll,b=ll,error=linear_error)
  else
    call lapack_linear_solver(a=lll,b=ll,error=linear_error)
  end if

  if (debug.and.linear_error) write(83,*) 'linear error coming out of linear solver routine'

! check on validity of ll elements
  if (.not.linear_error) then
    do j = 1, ubound(ll,1)
      if (.not.number_is_valid(ll(j,1))) then
        linear_error = .true.
        exit
      end if
    end do
    if (debug.and.linear_error) write(83,*) 'linear error from NaN'
  end if

! also check on size of ll kernel elements
  if (.not.linear_error) then
    if (maxval(abs(ll(1:nn,1))) > 1.d+20) linear_error = .true.
    if (debug.and.linear_error) write(83,*) 'linear error from size'
  end if

! if all OK then update y and check on llnorm
  if (.not.linear_error) then
    y = y + ll(:,1)
    call optimisation_kernel_update(pp,y,weight,weight_importance,separation_array,d,nn,mm,active,ll(:,1),lll,llnorm,l)
    if (llnorm > llnorm_tol) linear_error = .true.
    if (debug.and.linear_error) write(83,*) 'linear error from large llnorm = ',llnorm
  end if

! run check on differential kernels
  if (.false.) then
    call optimisation_kernel_check(pp,y,weight,weight_importance,separation_array,d,nn,mm,active)
  end if

  call optimisation_kernel_constraints(pp,y,weight,separation_array,minimum_separation,d,nn,mm, &
    active,active_change,phase,linear_error,error)

  if (error) call error_stop('unresolvable error during '//trim(phase)//' phase of optimisation_kernel_constraints')

  if (.not.active_change) exit
  if (n_kernel_steps > maximum_constraint_steps) &
    call error_stop('constraint loop is out of control in optimisation_kernel_constraints')

  call optimisation_kernel_update(pp,y,weight,weight_importance,separation_array,d,nn,mm,active,ll(:,1),lll,llnorm,l)

end do constraint_loop

if (debug) write(83,*) ' CONSTRAINT LOOP CONVERGED'

kernel = y(1:nn)

deallocate(y,ll,lll,active,weight,weight_importance)

if (debug) write(83,'(a/80(1h-))') 'subroutine optimisation_kernel'

end subroutine optimisation_kernel

!-----------------------------------------------------------------

subroutine optimisation_kernel_update(pp,y,weight,weight_importance,separation_array,d,nn,mm,active,ll,lll,llnorm,l)

! here we calculate the f vector and lll tensor from the latest y vect
! this is all for solving the optimisation kernel problem

use general_module
integer :: d, nn, mm, i, j
integer, dimension(:), allocatable :: active, weight_importance
double precision, dimension(:) :: ll
double precision, dimension(:), allocatable :: y, weight
double precision, dimension(:,:), allocatable :: pp, lll
integer, dimension(:), allocatable :: separation_array
double precision :: llnorm, l, pp_sign, weight_factor, rr, weight_tmp, pp_tmp
logical :: debug = .false.

if (debug) write(83,'(80(1h+)/a)') 'subroutine optimisation_kernel_update'

l = 0.d0
ll = 0.d0
lll = 0.d0

!----------------
! this stuff due to weight maximisation
do i = 1, nn

! rr = vector_magnitude(pp(min(2,mm):min(4,mm),i))
! weight_factor = rr/weight(i)**2
! weight_factor = rr/weight(i) ! this is the one
! weight_factor = rr/(weight(i)*radial_kernel(rr,0))
! weight_factor = 1.d0/weight(i)
! weight_factor = rr
! weight_factor = 1.d0
! weight_factor = 1.d0/weight(i) - 1.d0/max(maxval(weight),0.5d0)
! weight_factor = 1.d0/weight(i) - 0.99d0/maxval(weight)
! weight_factor = weight_factor**2
    
! weight_factor = 10.d0**weight_importance(i)
  weight_factor = 1.d0/(weight_separation_amount**weight_importance(i))
! weight_factor = 1.d0/(max(abs(weight(i)),1.d-2))**2

  if (radial_kernel_weighting) then
    l = l + (y(i)-weight(i))**2*weight_factor
    ll(i) = ll(i) + 2*(y(i)-weight(i))*weight_factor
    lll(i,i) = lll(i,i) + 2*weight_factor
  else
    l = l + (y(i)*pp(d,i)-weight(i))**2*weight_factor
    ll(i) = ll(i) + 2*(y(i)*pp(d,i)-weight(i))*pp(d,i)*weight_factor
    lll(i,i) = lll(i,i) + 2*pp(d,i)**2*weight_factor
!   l = l + (y(i)*pp(d,i))**2*weight_factor
!   ll(i) = ll(i) + 2*(y(i)*pp(d,i))*pp(d,i)*weight_factor
!   lll(i,i) = lll(i,i) + 2*pp(d,i)**2*weight_factor
  end if

end do

!----------------
! this stuff due to the polynomial conditions
do j = 1, nn
  do i = nn+1, nn+mm
    ll(j) = ll(j) + y(i)*pp(i-nn,j)
    lll(j,i) = lll(j,i) + pp(i-nn,j)
    lll(i,j) = lll(i,j) + pp(i-nn,j)
  end do
end do

do i = nn+1, nn+mm
  do j = 1, nn
    l = l + y(i)*y(j)*pp(i-nn,j)
    ll(i) = ll(i) + y(j)*pp(i-nn,j)
  end do
  if (i-nn == d) then
    l = l - y(i)
    ll(i) = ll(i) - 1.d0
  end if
end do

!----------------
! this stuff due to the positivity conditions
do i = 1, nn
  if (active(i) == 1) then
!   pp_sign = sign(1.d0,pp(d,i))
    pp_sign = sign(max(abs(pp(d,i)),small_pp),pp(d,i)) ! this gives something that has the same sign as pp(d,i), but a magnitude that is >= small_pp
    l = l - y(i+mm+nn)*y(i)*pp_sign
    ll(i) = ll(i) - y(i+mm+nn)*pp_sign
    ll(i+nn+mm) = ll(i+nn+mm) - y(i)*pp_sign
    lll(i,i+mm+nn) = lll(i,i+mm+nn) - pp_sign
    lll(i+mm+nn,i) = lll(i+mm+nn,i) - pp_sign
  else
    ll(i+nn+mm) = ll(i+nn+mm) + y(i+nn+mm) - 1.d0
    lll(i+nn+mm,i+nn+mm) = lll(i+nn+mm,i+nn+mm) + 1.d0
  end if
end do

!----------------
  
llnorm = 0.d0
do j = 1, ubound(ll,1)
  llnorm = llnorm + ll(j)**2
end do
! normalise llnorm by maximum element size in ll now
llnorm = llnorm/maxval(abs(y))

if (debug) then
  write(83,*) 'j:active:k/e/i:k*pp/-/-:ll:lll'
  do j = 1, nn
    write(83,'(i2,100(a,g9.2))') j,':-:k',y(j),':',y(j)*pp(d,j),':',ll(j),':',(lll(j,i),' ',i=1,ubound(ll,1))
  end do
  do j = nn+1, nn+mm
    write(83,'(i2,100(a,g9.2))') j,':-:e',y(j),':---------:',ll(j),':',(lll(j,i),' ',i=1,ubound(ll,1))
  end do
  do j = mm+nn+1, mm+2*nn
    write(83,'(i2,a,i1,99(a,g9.2))') j,':',active(j-mm-nn),':i',y(j),':---------:',ll(j),':',(lll(j,i),' ',i=1,ubound(ll,1))
  end do
  write(83,*) 'llnorm = ',llnorm
  write(83,*) 'l = ',l
end if

if (debug) write(83,'(a/80(1h-))') 'subroutine optimisation_kernel_update'

end subroutine optimisation_kernel_update

!-----------------------------------------------------------------

subroutine optimisation_kernel_check(pp,y,weight,weight_importance,separation_array,d,nn,mm,active)

! here we check the analytical optimisaion kernel derivatives against differences

use general_module
integer :: d, nn, mm, i, j
double precision, dimension(:), allocatable :: ll, ll_difference, ll_u, ll_d
double precision, dimension(:), allocatable :: y, weight, y_difference
double precision, dimension(:,:), allocatable :: pp, lll, lll_difference
integer, dimension(:), allocatable :: active, weight_importance
integer, dimension(:), allocatable :: separation_array
double precision :: llnorm, l, l_u, l_d, ll_error, lll_error
integer :: jj
logical :: debug = .true.
double precision :: y_eps = 1.d-6 ! small difference increment to y

if (debug) write(83,'(80(1h+)/a)') 'subroutine optimisation_kernel_check'

jj = 2*nn+mm
allocate(ll(jj),lll(jj,jj))
allocate(y_difference(jj),ll_difference(jj),lll_difference(jj,jj))
allocate(ll_u(jj),ll_d(jj))

! calculate difference approximations
do j = 1, jj
  y_difference = y
  y_difference(j) = y_difference(j) + y_eps
  call optimisation_kernel_update(pp,y_difference,weight,weight_importance,separation_array,d,nn,mm,active,ll_u,lll, &
    llnorm,l_u)
  y_difference(j) = y_difference(j) - 2.d0*y_eps
  call optimisation_kernel_update(pp,y_difference,weight,weight_importance,separation_array,d,nn,mm,active,ll_d,lll, &
    llnorm,l_d)
  ll_difference(j) = (l_u-l_d)/(2.d0*y_eps)
  do i = 1, jj
    lll_difference(i,j) = (ll_u(i)-ll_d(i))/(2.d0*y_eps)
  end do
end do

! calculate differential result
call optimisation_kernel_update(pp,y,weight,weight_importance,separation_array,d,nn,mm,active,ll,lll,llnorm,l)

if (debug) then
  write(83,*) 'j:y:ll:lll: differential result'
  do j = 1, jj
    write(83,'(i2,100(a,f8.2))') j,':',y(j),':',ll(j),':',(lll(j,i),' ',i=1,jj)
  end do
  write(83,*) 'j:y:ll:lll: difference result'
  do j = 1, jj
    write(83,'(i2,100(a,f8.2))') j,':',y(j),':',ll_difference(j),':',(lll_difference(j,i),' ',i=1,jj)
  end do
  write(83,*) 'j:y:ll:lll: error'
  do j = 1, jj
    write(83,'(i2,100(a,f8.2))') j,':',y(j),':',ll(j)-ll_difference(j),':',(lll(j,i)-lll_difference(j,i),' ',i=1,jj)
  end do
  ll_error = 0.d0
  lll_error = 0.d0
  do j = 1, jj
    ll_error = ll_error + abs(ll(j)-ll_difference(j))
    do i = 1, jj
      lll_error = lll_error + abs(lll(j,i)-lll_difference(j,i))
    end do
  end do
  write(83,*) 'll_error = ',ll_error
  write(83,*) 'lll_error = ',lll_error
end if

deallocate(ll,lll,ll_difference,lll_difference,ll_u,ll_d)

if (debug) write(83,'(a/80(1h-))') 'subroutine optimisation_kernel_check'

end subroutine optimisation_kernel_check

!-----------------------------------------------------------------

subroutine optimisation_kernel_constraints(pp,y,weight,separation_array,minimum_separation,d,nn,mm,active, &
  active_change,phase,linear_error,error)

! here we calculate set any constraints on the kernel variables

use general_module
integer :: d, nn, mm, i, j, nfree, minimum_separation, answer, i_min, i_trial, i_negative
integer, dimension(:), allocatable :: active, active_tmp
integer, dimension(:), allocatable, save :: active_last, active_best
double precision, dimension(:), allocatable :: y, weight
double precision, dimension(:), allocatable, save :: y_last, y_best
double precision, dimension(:,:), allocatable :: pp
integer, dimension(:), allocatable :: separation_array
double precision :: pp_sign, pp_min, pp_max, lambda_index, lambda_index_min, negative_index, &
  local_negative_index, trial_max, trial_local, distant_negative_index, maximum_negative_index
double precision, save :: negative_index_last, negative_index_best
integer, save :: i_negative_last, minimum_kernel_size, i_negative_best
logical :: active_change ! indicates whether any changes have been made to the constraints
character(len=10) :: phase ! first|second|expand|contract|finish : this is the phase of development we are in
logical :: linear_error ! error logical from the linear system solution
logical :: error
logical :: debug = .false.

if (debug) write(83,'(80(1h+)/a)') 'subroutine optimisation_kernel_constraints'

active_change = .false.
error = .false.

if (debug) then
  write(83,*) 'entering:'
  write(83,*) 'phase = ',phase
  write(83,*) 'linear_error = ',linear_error
  write(83,*) 'negative_index_last = ',negative_index_last
  write(83,*) 'i_negative_last = ',i_negative_last
  write(83,*) 'minimum_kernel_size = ',minimum_kernel_size
end if

!----------------
! setup constraints

if (trim(phase) == 'first'.or.(trim(phase) == 'second'.and.linear_error)) then

  active_change = .true.
  active = 0 ! initialise all to zero 
  i_negative = -1 ! dummy value
  i_negative_last = -1 ! dummy value
  i_negative_best = -1 ! dummy value
  negative_index = -huge(1.d0) ! dummy value
  negative_index_last = negative_index
  negative_index_best = negative_index

! setup initial constraints
  do i = 1, nn
    if (separation_array(i) > minimum_separation) then
      active(i) = 1 ! start by putting constraints on all elements beyond the minimum_separation
    else if (.not.radial_kernel_weighting.and.abs(pp(d,i)) < small_pp) then
      if (vector_magnitude(pp(min(2,mm):mm,i)) < small_pp) then
        active(i) = 0 ! if polymer basis in this dimension is small and we ontop of the kernel centre then make constraint inactive
      else
        active(i) = 1 ! otherwise make the constraint active to set this kernel element to zero
      end if
    end if
  end do

! now check on the number of free variables: at least as many variables must be free as the minimum_kernel_size
  nfree = nn - sum(active) ! this is the number of kernel elements within the unconstrained kernel

! set the minimum number of kernel elements
  if (trim(phase) == 'first') then
! mm is the smallest possible kernel size that could satisfy the polynomial equations
! nfree is at this stage the number of kernel elements that result from all the applicable minimum_separation elements include in the kernel
! optimise_additional_elements is a factor to increase the kernel size
    minimum_kernel_size = max(mm,nfree+optimise_additional_elements)
  else
! on subsequent calls increase the kernel size based on what was used last time
    minimum_kernel_size = nn - sum(active_last) + 1 ! if there was an error then we need to increase the number of free variables
  end if

! in first instance try to remove some constraints above the minimum_separation
  if (nfree < minimum_kernel_size) then !
    do i = 1, nn
      if (separation_array(i) > minimum_separation.and.(abs(pp(d,i)) >= small_pp.or.radial_kernel_weighting)) then
        active(i) = 0
        nfree = nfree + 1
        if (nfree == minimum_kernel_size) exit
      end if
    end do
  end if

! if we still have a problem then move in from the outer elements setting constaints off irrespective of pp magnitude
  if (nfree < minimum_kernel_size) then
    do i = nn, 1, -1
      if (active(i) == 1) then
        active(i) = 0
        nfree = nfree + 1
        if (nfree == minimum_kernel_size) exit
      end if
    end do
  end if

! so by now if there is still a problem then not enough cells are available and we have an unresolvable problem
  if (nfree < minimum_kernel_size) then
    if (debug) write(83,*) 'ran out of elements in kernel mask without a successful linear solution being obtained: '// &
      'Either you should make the minimum kernel size larger, or something is amiss with the dimensions of the mesh.'
    error = .true.
  end if

  if (.not.error) then
    call copy_integer_array(original=active,copy=active_last)
    call copy_integer_array(original=active,copy=active_best)
    phase = 'second' ! this means that no kernel has been solved yet by the linear solver
  end if

!----------------
!else if (.false..or.d==1) then
else if (.false.) then
! this is a debugging tool for changing the constraints by hand

  write(*,*) 'introduce (+) or remove (-) new constraint?'
  read(*,*) answer
  if (answer > 0) then
    active(answer) = 1
    write(*,*) 'active(',answer,') = ',active(answer)
  else if (answer < 0) then
    active(-answer) = 0
    write(*,*) 'active(',-answer,') = ',active(-answer)
  end if
  active_change = .true.

!----------------
else if (.true..and.optimise_positise) then ! restrained kernel expansion

  if (trim(phase) == 'second') phase = 'expand' ! change phase to signify that one linear solution worked

! calculate latest statistics if a solution was found
  if (.not.linear_error) then
! loop through all free variables seeing if there are any negative pp*k elements, and find the location and value of the largest
    i_negative = 0
    maximum_negative_index = 0.d0
    negative_index = 0.d0
    do i = nn, 1, -1
      local_negative_index = y(i)*pp(d,i)
!     local_negative_index = y(i)*pp(d,i)/weight(i)
      if (active(i) == 0) then
        if (local_negative_index < maximum_negative_index) then ! this means that the element has current a negative pp*k which is largest in magnitude
          i_negative = i
          maximum_negative_index = local_negative_index
        end if
        if (local_negative_index < 0.d0) negative_index = negative_index + local_negative_index
      end if
    end do
  end if

  if (.not.optimise_positise_sum_index) negative_index = maximum_negative_index ! use maximum kernel value rather than sum

  if (linear_error.and.sum(active) > sum(active_last)) then
! this signifies that we have tried to increase the kernel size, and an error has occurred.  There is no way to resolve this, so exit with the currently-best kernel
! set the kernel to the best one we achieved and exit
    active = active_best
    y = y_best
    negative_index = negative_index_best
    i_negative = i_negative_best
    phase = 'finish'

  else if (linear_error.or.(trim(phase) == 'contract'.and.i_negative /= 0)) then
! we are here because there was an error in the linear solution after a contraction operation in the expand phase, or
!  during the contract phase a negative kernel element appeared
! eitherway, we need to rewind to the last saved kernel
    active = active_last
    y = y_last
    i_negative = i_negative_last
    negative_index = negative_index_last
! set minimum_kernel_size to current size to prevent the same step reoccurring
    minimum_kernel_size = nn - sum(active)
    if (trim(phase) == 'contract') phase = 'finish' ! the minimum_kernel_size would also achieve this, but make it explicit

  else if (trim(phase) /= 'contract'.and.optimise_positise_cautiously.and. &
    negative_index < optimise_positise_cautiously_multiplier*negative_index_last) then
! using the optimise_positise_cautiously option we only allow a change if it decreases the negative index
    active = active_best
    y = y_best
    i_negative = i_negative_best
    negative_index = negative_index_best
    phase = 'finish'

  else
! otherwise we are progressing, so save the current kernel
    call copy_integer_array(original=active,copy=active_last)
    call copy_double_precision_array(original=y,copy=y_last)
    i_negative_last = i_negative
    negative_index_last = negative_index
    if (negative_index >= negative_index_best) then
      call copy_integer_array(original=active,copy=active_best)
      call copy_double_precision_array(original=y,copy=y_best)
      i_negative_best = i_negative
      negative_index_best = negative_index
    end if
  end if
    
! now do any changes to decrease the negative index and/or kernel size
  if (trim(phase) /= 'finish') then

    if (nn-sum(active) > minimum_kernel_size.and.i_negative /= 0) then
! there is a negative kernel element and we have space to remove it
      active(i_negative) = 1
      active_change = .true.

    else if (i_negative /= 0) then
! there is still a negative kernel element that we will try to remove by adding negative lambda constrained cells, if possible

      lambda_index_min = 1.d0
      i_min = 0
      do i = 1, nn ! loop through all cells looking for a constrained element that has the largest negative lambda
!       if (separation_array(i) /= separation_array(max(i-1,1)).and.i_min > 0) exit ! we have found a candidate in the current separation level so use this
        if (separation_array(i) /= separation_array(max(i-1,1)).and.i_min > 0.and.separation_array(i) > minimum_separation) exit ! we have found a candidate in the current separation level so use this
        if (active(i) == 1.and.y(nn+mm+i) < 0.d0) then
          if (.not.radial_kernel_weighting.and.abs(pp(d,i)) < small_pp.and.vector_magnitude(pp(min(2,mm):mm,i)) > small_pp) cycle ! do not relax constraint on these elements that have same pp(d,i) as origin
          lambda_index = y(nn+mm+i) ! fiddle with this
!         lambda_index = y(nn+mm+i)/weight(i) ! fiddle with this
          if (lambda_index < lambda_index_min) then
            i_min = i
            lambda_index_min = lambda_index
          end if
        end if
      end do
      
      if (i_min > 0) then
! we have found a variable that could help the situation, so relax this constraint
        active(i_min) = 0
        active_change = .true.
      else
! otherwise no variables are left which could help the situation and there is nothing more that can be done so we never enter the contract phase
! set the kernel to the best one we achieved and exit
        active = active_best
        y = y_best
        negative_index = negative_index_best
        i_negative = i_negative_best
        phase = 'finish' ! not necessary to set this, but...
      end if

    else
! we are now in the contraction phase as i_negative == 0
      phase = 'contract'

    end if

  end if

end if
      
!if (debug.and.(.not.active_change.or.trim(phase)=='finish')) then
if (debug) then
  write(83,*) 'exiting:'
  write(83,*) 'phase = ',phase
  write(83,*) 'negative_index = ',negative_index
  write(83,*) 'negative_index_best = ',negative_index_best
  write(83,*) 'i_negative = ',i_negative
  write(83,*) 'minimum_kernel_size = ',minimum_kernel_size
  write(83,*) 'active_change = ',active_change
  write(83,*) 'i:active:active_last:active_best:sep:k:k*pp:lambda'
  do i = 1, nn
    write(83,'(i2,3(a,i1),a,i2,100(a,f12.2))') i,':',active(i),':',active_last(i),':',active_best(i), &
      ':',separation_array(i),':',y(i),':',y(i)*pp(d,i),':',y(i+mm+nn)
  end do
end if

if (debug) write(83,'(a/80(1h-))') 'subroutine optimisation_kernel_constraints'

end subroutine optimisation_kernel_constraints

!-----------------------------------------------------------------

end module kernel_module

!-----------------------------------------------------------------
