!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ALVARO BERNAL GARCIA   (abernal@iqn.upv.es)
! PhD Candidate
! ISIRYM , UPV (VALENCIA)
! 6/02/2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
module equations_module

implicit none
! allow all subroutines and functions to be public


double precision, allocatable, dimension(:,:) :: vol_poly_coeff
double precision, allocatable, dimension(:,:,:) :: surf_poly_coeff, grad_poly_coeff


contains

!----------------------------------------------------------------------------

subroutine allocate_meta_arrays

! here we allocate array elements of the meta data that corresponds to the input files

use general_module
logical, parameter :: debug = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! abernal
integer :: num_mat, i,mat
character*120 :: string_material
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  
if (debug) write(*,'(80(1h+)/a)') 'subroutine allocate_meta_arrays'

! allocate meta arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! abernal
open(unit=95,file='input',status='old')
read(95,*)
read(95,*) num_mat ! number of materials
close (95)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! allocate var and compound arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! abernal
allocate(var(7*num_mat)) ! (7*num_mat) constants {7 for each material}  
allocate(compound(7*num_mat))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! abernal
do mat=1,num_mat
    write (string_material,'(I120)') mat
    string_material=adjustl(string_material)
    !constant: <D_1>
    var(1+7*(mat-1))%name='<D_1_'//trim(string_material)//'>'
    var(1+7*(mat-1))%multiplier=1.d0
    var(1+7*(mat-1))%units='cm⁻¹'
    var(1+7*(mat-1))%centring='none'
    var(1+7*(mat-1))%type='constant'
    var(1+7*(mat-1))%rank='scalar'
    var(1+7*(mat-1))%relstep=0
    var(1+7*(mat-1))%deriv=0
    var(1+7*(mat-1))%region=''
    var(1+7*(mat-1))%compound_name='<D_1_'//trim(string_material)//'>'
    var(1+7*(mat-1))%compound_number=1+7*(mat-1)
    var(1+7*(mat-1))%someloop=0
    allocate(var(1+7*(mat-1))%component(1))
    var(1+7*(mat-1))%component=[ (1+7*(mat-1)) ]

    !constant: <D_2>
    var(2+7*(mat-1))%name='<D_2_'//trim(string_material)//'>'
    var(2+7*(mat-1))%multiplier=1.d0
    var(2+7*(mat-1))%units='cm⁻¹'
    var(2+7*(mat-1))%centring='none'
    var(2+7*(mat-1))%type='constant'
    var(2+7*(mat-1))%rank='scalar'
    var(2+7*(mat-1))%relstep=0
    var(2+7*(mat-1))%deriv=0
    var(2+7*(mat-1))%region=''
    var(2+7*(mat-1))%compound_name='<D_2_'//trim(string_material)//'>'
    var(2+7*(mat-1))%compound_number=(2+7*(mat-1))
    var(2+7*(mat-1))%someloop=0
    allocate(var(2+7*(mat-1))%component(1))
    var(2+7*(mat-1))%component=[ (2+7*(mat-1)) ]

    !constant: <sigma_a_1>
    var(3+7*(mat-1))%name='<sigma_a_1_'//trim(string_material)//'>'
    var(3+7*(mat-1))%multiplier=1.d0
    var(3+7*(mat-1))%units='cm⁻¹'
    var(3+7*(mat-1))%centring='none'
    var(3+7*(mat-1))%type='constant'
    var(3+7*(mat-1))%rank='scalar'
    var(3+7*(mat-1))%relstep=0
    var(3+7*(mat-1))%deriv=0
    var(3+7*(mat-1))%region=''
    var(3+7*(mat-1))%compound_name='<sigma_a_1_'//trim(string_material)//'>'
    var(3+7*(mat-1))%compound_number=(3+7*(mat-1))
    var(3+7*(mat-1))%someloop=0
    allocate(var(3+7*(mat-1))%component(1))
    var(3+7*(mat-1))%component=[ (3+7*(mat-1)) ]

    !constant: <nu_x_sigma_f_1>
    var(4+7*(mat-1))%name='<nu_x_sigma_f_1_'//trim(string_material)//'>'
    var(4+7*(mat-1))%multiplier=1.d0
    var(4+7*(mat-1))%units='cm⁻¹'
    var(4+7*(mat-1))%centring='none'
    var(4+7*(mat-1))%type='constant'
    var(4+7*(mat-1))%rank='scalar'
    var(4+7*(mat-1))%relstep=0
    var(4+7*(mat-1))%deriv=0
    var(4+7*(mat-1))%region=''
    var(4+7*(mat-1))%compound_name='<nu_x_sigma_f_1_'//trim(string_material)//'>'
    var(4+7*(mat-1))%compound_number=(4+7*(mat-1))
    var(4+7*(mat-1))%someloop=0
    allocate(var(4+7*(mat-1))%component(1))
    var(4+7*(mat-1))%component=[ (4+7*(mat-1)) ]

    !constant: <sigma_s_12>
    var(5+7*(mat-1))%name='<sigma_s_12_'//trim(string_material)//'>'
    var(5+7*(mat-1))%multiplier=1.d0
    var(5+7*(mat-1))%units='cm⁻¹'
    var(5+7*(mat-1))%centring='none'
    var(5+7*(mat-1))%type='constant'
    var(5+7*(mat-1))%rank='scalar'
    var(5+7*(mat-1))%relstep=0
    var(5+7*(mat-1))%deriv=0
    var(5+7*(mat-1))%region=''
    var(5+7*(mat-1))%compound_name='<sigma_s_12_'//trim(string_material)//'>'
    var(5+7*(mat-1))%compound_number=(5+7*(mat-1))
    var(5+7*(mat-1))%someloop=0
    allocate(var(5+7*(mat-1))%component(1))
    var(5+7*(mat-1))%component=[ (5+7*(mat-1)) ]

    !constant: <nu_x_sigma_f_2>
    var(6+7*(mat-1))%name='<nu_x_sigma_f_2_'//trim(string_material)//'>'
    var(6+7*(mat-1))%multiplier=1.d0
    var(6+7*(mat-1))%units='cm⁻¹'
    var(6+7*(mat-1))%centring='none'
    var(6+7*(mat-1))%type='constant'
    var(6+7*(mat-1))%rank='scalar'
    var(6+7*(mat-1))%relstep=0
    var(6+7*(mat-1))%deriv=0
    var(6+7*(mat-1))%region=''
    var(6+7*(mat-1))%compound_name='<nu_x_sigma_f_2_'//trim(string_material)//'>'
    var(6+7*(mat-1))%compound_number=(6+7*(mat-1))
    var(6+7*(mat-1))%someloop=0
    allocate(var(6+7*(mat-1))%component(1))
    var(6+7*(mat-1))%component=[ (6+7*(mat-1)) ]

    !constant: <sigma_a_2>
    var(7+7*(mat-1))%name='<sigma_a_2_'//trim(string_material)//'>'
    var(7+7*(mat-1))%multiplier=1.d0
    var(7+7*(mat-1))%units='cm⁻¹'
    var(7+7*(mat-1))%centring='none'
    var(7+7*(mat-1))%type='constant'
    var(7+7*(mat-1))%rank='scalar'
    var(7+7*(mat-1))%relstep=0
    var(7+7*(mat-1))%deriv=0
    var(7+7*(mat-1))%region=''
    var(7+7*(mat-1))%compound_name='<sigma_a_2_'//trim(string_material)//'>'
    var(7+7*(mat-1))%compound_number=(7+7*(mat-1))
    var(7+7*(mat-1))%someloop=0
    allocate(var(7+7*(mat-1))%component(1))
    var(7+7*(mat-1))%component=[ (7+7*(mat-1)) ]

 
    !compound: <D_1>
    compound(1+7*(mat-1))%name='<D_1_'//trim(string_material)//'>'
    compound(1+7*(mat-1))%units='cm'
    compound(1+7*(mat-1))%centring='none'
    compound(1+7*(mat-1))%type='constant'
    compound(1+7*(mat-1))%rank='scalar'
    compound(1+7*(mat-1))%relstep=0
    compound(1+7*(mat-1))%deriv=0
    compound(1+7*(mat-1))%region=''
    allocate(compound(1+7*(mat-1))%component(1))
    compound(1+7*(mat-1))%component=[ 1+7*(mat-1) ]

    !compound: <D_2>
    compound(2+7*(mat-1))%name='<D_2_'//trim(string_material)//'>'
    compound(2+7*(mat-1))%units='cm'
    compound(2+7*(mat-1))%centring='none'
    compound(2+7*(mat-1))%type='constant'
    compound(2+7*(mat-1))%rank='scalar'
    compound(2+7*(mat-1))%relstep=0
    compound(2+7*(mat-1))%deriv=0
    compound(2+7*(mat-1))%region=''
    allocate(compound(2+7*(mat-1))%component(1))
    compound(2+7*(mat-1))%component=[ 2+7*(mat-1) ]

    !compound: <sigma_a_1>
    compound(3+7*(mat-1))%name='<sigma_a_1_'//trim(string_material)//'>'
    compound(3+7*(mat-1))%units='cm?¹'
    compound(3+7*(mat-1))%centring='none'
    compound(3+7*(mat-1))%type='constant'
    compound(3+7*(mat-1))%rank='scalar'
    compound(3+7*(mat-1))%relstep=0
    compound(3+7*(mat-1))%deriv=0
    compound(3+7*(mat-1))%region=''
    allocate(compound(3+7*(mat-1))%component(1))
    compound(3+7*(mat-1))%component=[ 3+7*(mat-1) ]

    !compound: <nu_x_sigma_f_1>
    compound(4+7*(mat-1))%name='<nu_x_sigma_f_1_'//trim(string_material)//'>'
    compound(4+7*(mat-1))%units='cm?¹'
    compound(4+7*(mat-1))%centring='none'
    compound(4+7*(mat-1))%type='constant'
    compound(4+7*(mat-1))%rank='scalar'
    compound(4+7*(mat-1))%relstep=0
    compound(4+7*(mat-1))%deriv=0
    compound(4+7*(mat-1))%region=''
    allocate(compound(4+7*(mat-1))%component(1))
    compound(4+7*(mat-1))%component=[ 4+7*(mat-1) ]

    !compound: <sigma_s_12>
    compound(5+7*(mat-1))%name='<sigma_s_12_'//trim(string_material)//'>'
    compound(5+7*(mat-1))%units='cm?¹'
    compound(5+7*(mat-1))%centring='none'
    compound(5+7*(mat-1))%type='constant'
    compound(5+7*(mat-1))%rank='scalar'
    compound(5+7*(mat-1))%relstep=0
    compound(5+7*(mat-1))%deriv=0
    compound(5+7*(mat-1))%region=''
    allocate(compound(5+7*(mat-1))%component(1))
    compound(5+7*(mat-1))%component=[ 5+7*(mat-1) ]

    !compound: <nu_x_sigma_f_2>
    compound(6+7*(mat-1))%name='<nu_x_sigma_f_2_'//trim(string_material)//'>'
    compound(6+7*(mat-1))%units='cm?¹'
    compound(6+7*(mat-1))%centring='none'
    compound(6+7*(mat-1))%type='constant'
    compound(6+7*(mat-1))%rank='scalar'
    compound(6+7*(mat-1))%relstep=0
    compound(6+7*(mat-1))%deriv=0
    compound(6+7*(mat-1))%region=''
    allocate(compound(6+7*(mat-1))%component(1))
    compound(6+7*(mat-1))%component=[ 6+7*(mat-1) ]

    !compound: <sigma_a_2>
    compound(7+7*(mat-1))%name='<sigma_a_2_'//trim(string_material)//'>'
    compound(7+7*(mat-1))%units='cm?¹'
    compound(7+7*(mat-1))%centring='none'
    compound(7+7*(mat-1))%type='constant'
    compound(7+7*(mat-1))%rank='scalar'
    compound(7+7*(mat-1))%relstep=0
    compound(7+7*(mat-1))%deriv=0
    compound(7+7*(mat-1))%region=''
    allocate(compound(7+7*(mat-1))%component(1))
    compound(7+7*(mat-1))%component=[ 7+7*(mat-1) ]

enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


constants_file = "constants.in" ! filename that was used by setup_equations


! set transient_simulation logical
!<sub_string:transient_simulation>

if (debug) write(*,'(a/80(1h-))') 'subroutine allocate_meta_arrays'

end subroutine allocate_meta_arrays

!-----------------------------------------------------------------







subroutine update

! here the equations are discretized

use general_module

use materialscells_module ! Modulo hecho por abernal para saber los materiales de cada celda
use eigensolver_module ! Module developed by abernal. It is used to determine "trans_BC_virtcell" and vectors defining matrices



integer :: num_mat, boundary_conditions(6),mat, mat1, mat2
integer number_unknowns, bound_cond, i1, i2, i, j, cont, t, jcell, &
 & jcell1, jcell2, i12
integer, allocatable, dimension(:,:) :: pos_col_unknowns
double precision grad_term, div_term_1, div_term_2, value_bound_cond
double precision nu_x_sigma_f_1, nu_x_sigma_f_2, sigma_s_12, &
 & sigma_a_1, sigma_a_2, D_1, D_2, D_1_mat1, D_2_mat1, &
 & D_1_mat2, D_2_mat2
 double precision :: mult_adf1(2), mult_adf2(2)


! Determine the number of materials and the boundary conditions
open(unit=95,file='input',status='old')
read(95,*)
read(95,*) num_mat ! number of materials
read(95,*)
read(95,*) (boundary_conditions(i),i=1,6) ! boundary conditions of -X, +X, -Y, +Y, -Z, +Z
close (95)



! Instruction to determine the material in each cell
if(.not.(allocated(material_element))) call materials_cells !!!! 16/10/2014 , abernal,  acople



! Determine the number of unknowns
 number_unknowns=idomain+jdomain*2+jboundary
 

! Allocate vectors that will define the eigenvalue problem matrices
allocate (val_M11(number_unknowns))
allocate (val_M12(number_unknowns))
allocate (val_L21(number_unknowns))
allocate (col_M11_M12_L21(number_unknowns))
allocate (row_M11_M12_L21(idomain+1))
allocate (trans_row_M11_M12_L21(idomain))
allocate (row_L11_L22(number_unknowns+1))
allocate (col_L11_L22((idomain+jboundary+jdomain*2*2)*7)) ! 7 is the maximum number of polynomial coefficients for 3D cases
allocate (val_L11((idomain+jboundary+jdomain*2*2)*7))
allocate (val_L22((idomain+jboundary+jdomain*2*2)*7))



! Calculate the polynomial_coefficients
 call polynomial_coefficients


! Loop over the cells to determine the position of the unkown vector
allocate(pos_col_unknowns(idomain,7))
 cont=0
do i=1,idomain
   ! Loop over the polynomial terms
   do t=1,ubound(cell(i)%jface,1)+1
      cont=cont+1
      pos_col_unknowns(i,t)=cont
   enddo
enddo


! Loop over the cells to discretize the diffusion equation by using the FVM
num_col_L11_L22=0
num_col_M11_M12_L21=0
row_L11_L22(1)=1
row_M11_M12_L21(1)=1
 cont=0
do i=1,idomain
   mat=material_element(i) ! Determination of the material of this cell 
   ! Determination of cross sections
   nu_x_sigma_f_1=var(7*(mat-1)+4)%funk(1)%v
   nu_x_sigma_f_2=var(7*(mat-1)+6)%funk(1)%v
   sigma_s_12=var(7*(mat-1)+5)%funk(1)%v
   sigma_a_1=var(7*(mat-1)+3)%funk(1)%v
   sigma_a_2=var(7*(mat-1)+7)%funk(1)%v
   D_1=var(7*(mat-1)+1)%funk(1)%v
   D_2=var(7*(mat-1)+2)%funk(1)%v


   ! Loop over the polynomial terms
   do t=1,ubound(cell(i)%jface,1)+1
      num_col_L11_L22=num_col_L11_L22+1
      col_L11_L22(num_col_L11_L22)=pos_col_unknowns(i,t)

      num_col_M11_M12_L21=num_col_M11_M12_L21+1
      col_M11_M12_L21(num_col_M11_M12_L21)=pos_col_unknowns(i,t)

      val_M11(num_col_M11_M12_L21)=nu_x_sigma_f_1*vol_poly_coeff(i,t)
      val_M12(num_col_M11_M12_L21)=nu_x_sigma_f_2*vol_poly_coeff(i,t) 
      val_L21(num_col_M11_M12_L21)=-1*sigma_s_12*vol_poly_coeff(i,t)

      ! Loop over the faces to calculate the Divergence term for each polynomial term
      grad_term=0
      do jcell=1,ubound(cell(i)%jface,1)
         j=cell(i)%jface(jcell)
         grad_term=grad_term+grad_poly_coeff(i,jcell,t)*face(j)%area
      enddo
      div_term_1=(-1)*D_1/cell(i)%vol*grad_term ! divergence_1= -D1/vol*grad
      div_term_2=(-1)*D_2/cell(i)%vol*grad_term ! divergence_2= -D2/vol*grad

      val_L11(num_col_L11_L22)=(sigma_s_12+sigma_a_1)*vol_poly_coeff(i,t)+div_term_1 
      val_L22(num_col_L11_L22)=sigma_a_2*vol_poly_coeff(i,t)+div_term_2
   enddo
   cont=cont+1
   row_L11_L22(cont+1)=num_col_L11_L22+1
   row_M11_M12_L21(i+1)=num_col_M11_M12_L21+1
   trans_row_M11_M12_L21(i)=cont


   ! Loop over the faces surrounding this cell to impose the face equations: boundary conditions or inner face conditions
   do jcell=1,ubound(cell(i)%jface,1)
      j=cell(i)%jface(jcell)

      ! Determine if the face is a boundary face or an inner face
      select case(face(j)%type)

      case(1) ! Inner face

         ! Determine the adjacent cells of this face
         i1=face(j)%icell(1)
         i2=face(j)%icell(2)
         mat1=material_element(i1) ! Determination of the material of the first cell 
         mat2=material_element(i2) ! Determination of the material of the second cell
         D_1_mat1=var(7*(mat1-1)+1)%funk(1)%v
         D_2_mat1=var(7*(mat1-1)+2)%funk(1)%v
         D_1_mat2=var(7*(mat2-1)+1)%funk(1)%v
         D_2_mat2=var(7*(mat2-1)+2)%funk(1)%v 

         ! Determine if this cell is the first or the second adjacent cell
         if (i1==i) then
             i12=1
             ! Determine the surface number of the second cell containing this inner face
             jcell2=1
             do while (j/=cell(i2)%jface(jcell2))
                jcell2=jcell2+1
             enddo
             ! Determine the surface number of the first cell containing this inner face
             jcell1=jcell
         else
             i12=2
             ! Determine the surface number of the first cell containing this inner face
             jcell1=1
             do while (j/=cell(i1)%jface(jcell1))
                jcell1=jcell1+1
             enddo
             ! Determine the surface number of the second cell containing this inner face
             jcell2=jcell
         endif

         ! Impose the flux condition if this cell is the first adjacent cell.
         ! Impose the current condition if this cell is the second adjacent cell.
         select case(i12)

         case(1) ! Flux condition

            ! Determination of ADF if necessary
            call adf_cells(j,mat1,mat2,mult_adf1,mult_adf2)

            ! Loop over the polynomial terms of the first cell
            do t=1,ubound(cell(i1)%jface,1)+1
               num_col_L11_L22=num_col_L11_L22+1
               col_L11_L22(num_col_L11_L22)=pos_col_unknowns(i1,t)
               val_L11(num_col_L11_L22)=surf_poly_coeff(i1,jcell1,t)*mult_adf1(1)
               val_L22(num_col_L11_L22)=surf_poly_coeff(i1,jcell1,t)*mult_adf1(2)
            enddo
            ! Loop over the polynomial terms of the second cell
            do t=1,ubound(cell(i2)%jface,1)+1
               num_col_L11_L22=num_col_L11_L22+1
               col_L11_L22(num_col_L11_L22)=pos_col_unknowns(i2,t)
               val_L11(num_col_L11_L22)=-1*surf_poly_coeff(i2,jcell2,t)*mult_adf2(1)
               val_L22(num_col_L11_L22)=-1*surf_poly_coeff(i2,jcell2,t)*mult_adf2(2)
            enddo

         case(2) ! Current condition

            ! Loop over the polynomial terms of the first cell
            do t=1,ubound(cell(i1)%jface,1)+1
               num_col_L11_L22=num_col_L11_L22+1
               col_L11_L22(num_col_L11_L22)=pos_col_unknowns(i1,t)
               val_L11(num_col_L11_L22)=D_1_mat1*grad_poly_coeff(i1,jcell1,t)
               val_L22(num_col_L11_L22)=D_2_mat1*grad_poly_coeff(i1,jcell1,t)
            enddo
            ! Loop over the polynomial terms of the second cell
            do t=1,ubound(cell(i2)%jface,1)+1
               num_col_L11_L22=num_col_L11_L22+1
               col_L11_L22(num_col_L11_L22)=pos_col_unknowns(i2,t)
               val_L11(num_col_L11_L22)=D_1_mat2*grad_poly_coeff(i2,jcell2,t)
               val_L22(num_col_L11_L22)=D_2_mat2*grad_poly_coeff(i2,jcell2,t)
            enddo

         end select


      case(2) ! Boundary face

         ! Determine the boundary condition: 0=zero flux, 1=reflective flux
         if (region(face(j)%region_list(1))%name=='<-x>') then
            bound_cond=boundary_conditions(1)
         elseif (region(face(j)%region_list(1))%name=='<+x>') then
            bound_cond=boundary_conditions(2)
         elseif (region(face(j)%region_list(1))%name=='<-y>') then
            bound_cond=boundary_conditions(3)
         elseif (region(face(j)%region_list(1))%name=='<+y>') then
            bound_cond=boundary_conditions(4)
         elseif (region(face(j)%region_list(1))%name=='<-z>') then
            bound_cond=boundary_conditions(5)
         elseif (region(face(j)%region_list(1))%name=='<+z>') then
            bound_cond=boundary_conditions(6)
         endif

         ! Assign the equation corresponding to the boundary condition to this face
         select case(bound_cond)
         case(0) ! Zero flux
            ! Loop over the polynomial terms
            do t=1,ubound(cell(i)%jface,1)+1
               num_col_L11_L22=num_col_L11_L22+1
               col_L11_L22(num_col_L11_L22)=pos_col_unknowns(i,t)
               val_L11(num_col_L11_L22)=surf_poly_coeff(i,jcell,t)
               val_L22(num_col_L11_L22)=surf_poly_coeff(i,jcell,t)
            enddo
         case(1) ! Reflective flux
            ! Loop over the polynomial terms
            do t=1,ubound(cell(i)%jface,1)+1
               num_col_L11_L22=num_col_L11_L22+1
               col_L11_L22(num_col_L11_L22)=pos_col_unknowns(i,t)
               val_L11(num_col_L11_L22)=grad_poly_coeff(i,jcell,t)
               val_L22(num_col_L11_L22)=grad_poly_coeff(i,jcell,t)
            enddo
         case default ! ERROR
            write(*,*) 'ERROR: Some boundary condition is not 1 (reflective flux) or 0 (zero flux)'
            stop
         end select

      end select
      cont=cont+1
      row_L11_L22(cont+1)=num_col_L11_L22+1
   enddo
enddo



! Deallocate "pos_col_unknowns"
deallocate(pos_col_unknowns)


end subroutine update













! This subroutine calculates the constant multiplying each polynomial coeffcient for each cell
subroutine polynomial_coefficients

use general_module


 ! Local variables
 integer i,j, coef_poly(3,7) ! Maximum: 7 polynomial terms


 ! Allocate variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! abernal, 16/10/14,   acople
! allocate(vol_poly_coeff(idomain,7)) ! "idomain" is defined in "general_module"
! allocate(surf_poly_coeff(idomain,6,7)) ! Maximum: 6 faces, 7 polynomial terms
! allocate(grad_poly_coeff(idomain,6,7))
 if (.not.allocated(vol_poly_coeff)) allocate(vol_poly_coeff(idomain,7))
 if (.not.allocated(surf_poly_coeff)) allocate(surf_poly_coeff(idomain,6,7))
 if (.not.allocated(grad_poly_coeff)) allocate(grad_poly_coeff(idomain,6,7))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! Read polynomial coefficients fron the file "input"
 open(unit=95,file='input',status='old')
 do i=1,9
    read(95,*)
 enddo
 read(95,*) ((coef_poly(j,i),j=1,3),i=1,ubound(cell(1)%jface,1)+1) 
 close (95)

 ! Polynomial approximation: phi(x,y,z)=

 ! Loop over the cells
 do i=1,idomain
    select case (ubound(cell(i)%jface,1))
    case(3) ! Triangle
       ! Calculate volume averaged values multipliers
          call vol_poly_coef_triangle_orderN(i,coef_poly)
       ! Calculate surface and gradient averaged values multipliers
          call surf_and_grad_poly_coef_triangle_orderN(i,coef_poly)
    case(4) ! Tetraedron or Cuadrangle
       select case(cell(i)%dimensions)
       case(2) ! Cuadrangle
          ! Calculate volume averaged values multipliers
             call vol_poly_coef_cuadrangle_orderN(i,coef_poly)
          ! Calculate surface and gradient averaged values multipliers
             call surf_and_grad_poly_coef_cuadrangle_orderN(i,coef_poly)
       case(3) ! Tetraedron
          ! Calculate volume averaged values multipliers
             call vol_poly_coef_tetraedron_orderN(i,coef_poly)
          ! Calculate surface and gradient averaged values multipliers
             call surf_and_grad_poly_coef_tetraedron_orderN(i,coef_poly)
       end select
    case(6) ! Hexaedron
       ! Calculate volume averaged values multipliers
          call vol_poly_coef_hexaedron_orderN(i,coef_poly)
       ! Calculate surface and gradient averaged values multipliers
          call surf_and_grad_poly_coef_hexaedron_orderN(i,coef_poly)
    case default
       write(*,*) 'ERROR: The mesh contains an element which is not a tetraedron , a hexaedron, a triangle or a cuadrangle'
       stop
    end select

 enddo



end subroutine polynomial_coefficients















! This subroutine calculates the cell averaged valaues multiplying each polynomial coeffcient for a tetraedron cell
! by using the following polynomial approach: 
! phi(x,y,z)= a + b*x + c*y + d*z + e*x**2 + f*y**2 + g*z**2 + h*x*y + k*x*z + l*y*z + m*x**3 + n*y**3 + o*z**3 + 
!           + p*x**2*y + q*x**2*z + r*y**2*x + s*y**2*z + t*z**2*x + u*z**2*y + v*x*y*z
! But only with 5 polynomial terms
subroutine vol_poly_coef_tetraedron(c,vect_coef)

use general_module


 ! Global variables
 integer c, vect_coef(5)

 ! Local variables
 double precision Jx, Jy, Jz, Pxy, Pxz, Pyz, x1, x2, x3, x4, y1, y2, y3, y4, &
  & z1, z2, z3, z4, vol_poly_coef, x0(4,3)
 integer i, alpha, beta, gamma, delta, j

 ! Definition of: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
 x1=node(cell(c)%knode(1))%x(1)
 y1=node(cell(c)%knode(1))%x(2)
 z1=node(cell(c)%knode(1))%x(3)
 x2=node(cell(c)%knode(2))%x(1)
 y2=node(cell(c)%knode(2))%x(2)
 z2=node(cell(c)%knode(2))%x(3)
 x3=node(cell(c)%knode(3))%x(1)
 y3=node(cell(c)%knode(3))%x(2)
 z3=node(cell(c)%knode(3))%x(3)
 x4=node(cell(c)%knode(4))%x(1)
 y4=node(cell(c)%knode(4))%x(2)
 z4=node(cell(c)%knode(4))%x(3)



 ! Definition of tetrahedron coordinates
 do i=1,3
    do j=1,4
       x0(i,j)=node(cell(c)%knode(j))%x(i)
    enddo
 enddo



 ! The inertia tensor has been calculated by using the formulas exposed in paper:
 ! "Explicit Exact Formulas for the 3-D Tetrahedron Inertia Tensor in Terms of its Vertex Coordinates"

 Jx=6*cell(c)%vol/60.0*(y1**2 + y1*y2 + y2**2 + y1*y3 + y2*y3 + &
    & y3**2 + y1*y4 + y2*y4 + y3*y4 + y4**2 + z1**2 + z1*z2 + &
    & z2**2 + z1*z3 + z2*z3 + z3**2 + z1*z4 + z2*z4 + z3*z4 + z4**2)

 Jy=6*cell(c)%vol/60.0*(x1**2 + x1*x2 + x2**2 + x1*x3 + x2*x3 + &
    & x3**2 + x1*x4 + x2*x4 + x3*x4 + x4**2 + z1**2 + z1*z2 + &
    & z2**2 + z1*z3 + z2*z3 + z3**2 + z1*z4 + z2*z4 + z3*z4 + z4**2)

 Jz=6*cell(c)%vol/60.0*(x1**2 + x1*x2 + x2**2 + x1*x3 + x2*x3 + &
    & x3**2 + x1*x4 + x2*x4 + x3*x4 + x4**2 + y1**2 + y1*y2 + &
    & y2**2 + y1*y3 + y2*y3 + y3**2 + y1*y4 + y2*y4 + y3*y4 + y4**2)

 Pyz=6*cell(c)%vol/120.0*(2*y1*z1 + y2*z1 + y3*z1 + y4*z1 + y1*z2 + &
     & 2*y2*z2 + y3*z2 + y4*z2 + y1*z3 + y2*z3 + 2*y3*z3 + &
     & y4*z3 + y1*z4 + y2*z4 + y3*z4 + 2*y4*z4)

 Pxz=6*cell(c)%vol/120.0*(2*x1*z1 + x2*z1 + x3*z1 + x4*z1 + x1*z2 + &
     & 2*x2*z2 + x3*z2 + x4*z2 + x1*z3 + x2*z3 + 2*x3*z3 + &
     & x4*z3 + x1*z4 + x2*z4 + x3*z4 + 2*x4*z4)

 Pxy=6*cell(c)%vol/120.0*(2*x1*y1 + x2*y1 + x3*y1 + x4*y1 + x1*y2 + &
     & 2*x2*y2 + x3*y2 + x4*y2 + x1*y3 + x2*y3 + 2*x3*y3 + &
     & x4*y3 + x1*y4 + x2*y4 + x3*y4 + 2*x4*y4)



 ! phi(x,y,z)=a + b*x + c*y + d*z + e*x**2 + f*y**2 + g*z**2 + h*x*y + k*x*z + l*y*z + m*x*y*z
 do i=1,5
    select case (vect_coef(i))
    case (1) ! a : 1
       vol_poly_coef=1
    case (2) ! b*x : Xg  (coordinate X of geometry center)
       vol_poly_coef=cell(c)%x(1)
    case (3) ! c*y : Yg  (coordinate Y of geometry center)
       vol_poly_coef=cell(c)%x(2)
    case (4) ! d*z : Zg  (coordinate Z of geometry center)
       vol_poly_coef=cell(c)%x(3)
    case (5) ! e*x**2 : (Jy+Jz-Jx)/(2*volume) (Moment of Inertia of plane X / volume)
       vol_poly_coef=(Jy+Jz-Jx)/(2*cell(c)%vol)
    case (6) ! f*y**2 : (Jx+Jz-Jy)/(2*volume) (Moment of Inertia of plane Y / volume)
       vol_poly_coef=(Jx+Jz-Jy)/(2*cell(c)%vol)
    case (7) ! g*z**2 : (Jx+Jy-Jz)/(2*volume) (Moment of Inertia of plane Z / volume)
       vol_poly_coef=(Jx+Jy-Jz)/(2*cell(c)%vol)
    case (8) ! h*x*y : Pxy/volume (Product of Inertia of plane XY / volume)
       vol_poly_coef=Pxy/cell(c)%vol
    case (9) ! k*x*z : Pxz/volume (Product of Inertia of plane XZ / volume)
       vol_poly_coef=Pxz/cell(c)%vol
    case (10) ! l*y*z : Pyz/volume (Product of Inertia of plane YZ / volume)
       vol_poly_coef=Pyz/cell(c)%vol
    case(11) ! m*x**3
       alpha=1
       beta=1
       gamma=1
       call  integrate_order3_volume_tetraedron(c,alpha,beta,gamma,vol_poly_coef)
    case(12) ! n*y**3
       alpha=2
       beta=2
       gamma=2
       call  integrate_order3_volume_tetraedron(c,alpha,beta,gamma,vol_poly_coef)
    case(13) ! o*z**3
       alpha=3
       beta=3
       gamma=3
       call  integrate_order3_volume_tetraedron(c,alpha,beta,gamma,vol_poly_coef)
    case(14) ! p*x**2*y
       alpha=1
       beta=1
       gamma=2
       call  integrate_order3_volume_tetraedron(c,alpha,beta,gamma,vol_poly_coef)
    case(15) ! q*x**2*z
       alpha=1
       beta=1
       gamma=3
       call  integrate_order3_volume_tetraedron(c,alpha,beta,gamma,vol_poly_coef)
    case(16) ! r*y**2*x
       alpha=2
       beta=2
       gamma=1
       call  integrate_order3_volume_tetraedron(c,alpha,beta,gamma,vol_poly_coef)
    case(17) ! s*y**2*z
       alpha=2
       beta=2
       gamma=3
       call  integrate_order3_volume_tetraedron(c,alpha,beta,gamma,vol_poly_coef)
    case(18) ! t*z**2*x
       alpha=3
       beta=3
       gamma=1
       call  integrate_order3_volume_tetraedron(c,alpha,beta,gamma,vol_poly_coef)
    case(19) ! u*z**2*y
       alpha=3
       beta=3
       gamma=2
       call  integrate_order3_volume_tetraedron(c,alpha,beta,gamma,vol_poly_coef)
    case(20) ! v*x*y*z
       alpha=1
       beta=2
       gamma=3
       call  integrate_order3_volume_tetraedron(c,alpha,beta,gamma,vol_poly_coef)
    case(21) ! x**4
       alpha=1
       beta=1
       gamma=1
       delta=1
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(22) ! y**4
       alpha=2
       beta=2
       gamma=2
       delta=2
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(23) ! z**4
       alpha=3
       beta=3
       gamma=3
       delta=3
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(24) ! x**3*y
       alpha=1
       beta=1
       gamma=1
       delta=2
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(25) ! x**3*z
       alpha=1
       beta=1
       gamma=1
       delta=3
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(26) ! y**3*x
       alpha=2
       beta=2
       gamma=2
       delta=1
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(27) ! y**3*z
       alpha=2
       beta=2
       gamma=2
       delta=3
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(28) ! z**3*x
       alpha=3
       beta=3
       gamma=3
       delta=1
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(29) ! z**3*y
       alpha=3
       beta=3
       gamma=3
       delta=2
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(30) ! x**2*y**2
       alpha=1
       beta=1
       gamma=2
       delta=2
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(31) ! x**2*z**2
       alpha=1
       beta=1
       gamma=3
       delta=3
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(32) ! y**2*z**2
       alpha=2
       beta=2
       gamma=3
       delta=3
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(33) ! x**2*y*z
       alpha=1
       beta=1
       gamma=2
       delta=3
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(34) ! y**2*x*z
       alpha=2
       beta=2
       gamma=1
       delta=3
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    case(35) ! z**2*x*y
       alpha=3
       beta=3
       gamma=1
       delta=2
       call  integrate_order4_volume_tetraedron_coordinates(x0,alpha,beta,gamma,delta,vol_poly_coef)
    end select
    vol_poly_coeff(c,i)=vol_poly_coef
 enddo


end subroutine vol_poly_coef_tetraedron















! This subroutine calculates the cell averaged valaues multiplying each polynomial coeffcient for a hexaedron cell
! by using the following polynomial approach: 
! phi(x,y,z)= a + b*x + c*y + d*z + e*x**2 + f*y**2 + g*z**2 + h*x*y + k*x*z + l*y*z + m*x**3 + n*y**3 + o*z**3 + 
!           + p*x**2*y + q*x**2*z + r*y**2*x + s*y**2*z + t*z**2*x + u*z**2*y + v*x*y*z
! But only with 7 polynomial terms
subroutine vol_poly_coef_hexaedron(c,vect_coef)

use general_module


 ! Global variables
 integer c, vect_coef(7)

 ! Local variables
 double precision Jx, Jy, Jz, Pxy, Pxz, Pyz, x1, x2, x3, x4, y1, y2, y3, y4, &
  & z1, z2, z3, z4, vol_poly_coef, vol_tetraedron
 integer :: i, knode(4), nodes_tetraedron(5,4), k
 integer alpha, beta, gamma, delta



 ! The hexaedron is divided into 5 tetraedron and the inertia tensor is calculated by adding the inertia tensor of each tetraedron
 Jx=0
 Jy=0
 Jz=0
 Pxy=0
 Pxz=0
 Pyz=0
 nodes_tetraedron(1,:)=(/1,2,3,6/)
 nodes_tetraedron(2,:)=(/1,5,8,6/)
 nodes_tetraedron(3,:)=(/1,3,8,6/)
 nodes_tetraedron(4,:)=(/7,3,8,6/)
 nodes_tetraedron(5,:)=(/1,3,8,4/)
 do i=1,5 ! Loop over the 5 tetraedra
    ! Definition of: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
    do k=1,4
       knode(k)=cell(c)%knode(nodes_tetraedron(i,k))
    enddo
    x1=node(knode(1))%x(1)
    y1=node(knode(1))%x(2)
    z1=node(knode(1))%x(3)
    x2=node(knode(2))%x(1)
    y2=node(knode(2))%x(2)
    z2=node(knode(2))%x(3)
    x3=node(knode(3))%x(1)
    y3=node(knode(3))%x(2)
    z3=node(knode(3))%x(3)
    x4=node(knode(4))%x(1)
    y4=node(knode(4))%x(2)
    z4=node(knode(4))%x(3)
    vol_tetraedron=abs((x2-x1)*((y3-y1)*(z4-z1)-(y4-y1)*(z3-z1)) &
                   & -(y2-y1)*((x3-x1)*(z4-z1)-(x4-x1)*(z3-z1)) &
                   & +(z2-z1)*((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1)))/6.0

    ! The inertia tensor of each tetraedron has been calculated by using the formulas exposed in paper:
    ! "Explicit Exact Formulas for the 3-D Tetrahedron Inertia Tensor in Terms of its Vertex Coordinates"

    Jx=Jx+6*vol_tetraedron/60.0*(y1**2 + y1*y2 + y2**2 + y1*y3 + y2*y3 + &
       & y3**2 + y1*y4 + y2*y4 + y3*y4 + y4**2 + z1**2 + z1*z2 + &
       & z2**2 + z1*z3 + z2*z3 + z3**2 + z1*z4 + z2*z4 + z3*z4 + z4**2)

    Jy=Jy+6*vol_tetraedron/60.0*(x1**2 + x1*x2 + x2**2 + x1*x3 + x2*x3 + &
       & x3**2 + x1*x4 + x2*x4 + x3*x4 + x4**2 + z1**2 + z1*z2 + &
       & z2**2 + z1*z3 + z2*z3 + z3**2 + z1*z4 + z2*z4 + z3*z4 + z4**2)

    Jz=Jz+6*vol_tetraedron/60.0*(x1**2 + x1*x2 + x2**2 + x1*x3 + x2*x3 + &
       & x3**2 + x1*x4 + x2*x4 + x3*x4 + x4**2 + y1**2 + y1*y2 + &
       & y2**2 + y1*y3 + y2*y3 + y3**2 + y1*y4 + y2*y4 + y3*y4 + y4**2)

    Pyz=Pyz+6*vol_tetraedron/120.0*(2*y1*z1 + y2*z1 + y3*z1 + y4*z1 + y1*z2 + &
        & 2*y2*z2 + y3*z2 + y4*z2 + y1*z3 + y2*z3 + 2*y3*z3 + &
        & y4*z3 + y1*z4 + y2*z4 + y3*z4 + 2*y4*z4)

    Pxz=Pxz+6*vol_tetraedron/120.0*(2*x1*z1 + x2*z1 + x3*z1 + x4*z1 + x1*z2 + &
        & 2*x2*z2 + x3*z2 + x4*z2 + x1*z3 + x2*z3 + 2*x3*z3 + &
        & x4*z3 + x1*z4 + x2*z4 + x3*z4 + 2*x4*z4)

    Pxy=Pxy+6*vol_tetraedron/120.0*(2*x1*y1 + x2*y1 + x3*y1 + x4*y1 + x1*y2 + &
        & 2*x2*y2 + x3*y2 + x4*y2 + x1*y3 + x2*y3 + 2*x3*y3 + &
        & x4*y3 + x1*y4 + x2*y4 + x3*y4 + 2*x4*y4)
 enddo




 ! phi(x,y,z)=a + b*x + c*y + d*z + e*x**2 + f*y**2 + g*z**2 + h*x*y + k*x*z + l*y*z
 do i=1,7
    select case (vect_coef(i))
    case (1) ! a : 1
       vol_poly_coef=1
    case (2) ! b*x : Xg  (coordinate X of geometry center)
       vol_poly_coef=cell(c)%x(1)
    case (3) ! c*y : Yg  (coordinate Y of geometry center)
       vol_poly_coef=cell(c)%x(2)
    case (4) ! d*z : Zg  (coordinate Z of geometry center)
       vol_poly_coef=cell(c)%x(3)
    case (5) ! e*x**2 : (Jy+Jz-Jx)/(2*volume) (Moment of Inertia of plane X / volume)
       vol_poly_coef=(Jy+Jz-Jx)/(2*cell(c)%vol)
    case (6) ! f*y**2 : (Jx+Jz-Jy)/(2*volume) (Moment of Inertia of plane Y / volume)
       vol_poly_coef=(Jx+Jz-Jy)/(2*cell(c)%vol)
    case (7) ! g*z**2 : (Jx+Jy-Jz)/(2*volume) (Moment of Inertia of plane Z / volume)
       vol_poly_coef=(Jx+Jy-Jz)/(2*cell(c)%vol)
    case (8) ! h*x*y : Pxy/volume (Product of Inertia of plane XY / volume)
       vol_poly_coef=Pxy/cell(c)%vol
    case (9) ! k*x*z : Pxz/volume (Product of Inertia of plane XZ / volume)
       vol_poly_coef=Pxz/cell(c)%vol
    case (10) ! l*y*z : Pyz/volume (Product of Inertia of plane YZ / volume)
       vol_poly_coef=Pyz/cell(c)%vol
    case(11) ! m*x**3
       alpha=1
       beta=1
       gamma=1
       call  integrate_order3_volume_hexaedron(c,alpha,beta,gamma,vol_poly_coef)
    case(12) ! n*y**3
       alpha=2
       beta=2
       gamma=2
       call  integrate_order3_volume_hexaedron(c,alpha,beta,gamma,vol_poly_coef)
    case(13) ! o*z**3
       alpha=3
       beta=3
       gamma=3
       call  integrate_order3_volume_hexaedron(c,alpha,beta,gamma,vol_poly_coef)
    case(14) ! p*x**2*y
       alpha=1
       beta=1
       gamma=2
       call  integrate_order3_volume_hexaedron(c,alpha,beta,gamma,vol_poly_coef)
    case(15) ! q*x**2*z
       alpha=1
       beta=1
       gamma=3
       call  integrate_order3_volume_hexaedron(c,alpha,beta,gamma,vol_poly_coef)
    case(16) ! r*y**2*x
       alpha=2
       beta=2
       gamma=1
       call  integrate_order3_volume_hexaedron(c,alpha,beta,gamma,vol_poly_coef)
    case(17) ! s*y**2*z
       alpha=2
       beta=2
       gamma=3
       call  integrate_order3_volume_hexaedron(c,alpha,beta,gamma,vol_poly_coef)
    case(18) ! t*z**2*x
       alpha=3
       beta=3
       gamma=1
       call  integrate_order3_volume_hexaedron(c,alpha,beta,gamma,vol_poly_coef)
    case(19) ! u*z**2*y
       alpha=3
       beta=3
       gamma=2
       call  integrate_order3_volume_hexaedron(c,alpha,beta,gamma,vol_poly_coef)
    case(20) ! v*x*y*z
       alpha=1
       beta=2
       gamma=3
       call  integrate_order3_volume_hexaedron(c,alpha,beta,gamma,vol_poly_coef)
    case(21) ! x**4
       alpha=1
       beta=1
       gamma=1
       delta=1
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(22) ! y**4
       alpha=2
       beta=2
       gamma=2
       delta=2
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(23) ! z**4
       alpha=3
       beta=3
       gamma=3
       delta=3
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(24) ! x**3*y
       alpha=1
       beta=1
       gamma=1
       delta=2
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(25) ! x**3*z
       alpha=1
       beta=1
       gamma=1
       delta=3
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(26) ! y**3*x
       alpha=2
       beta=2
       gamma=2
       delta=1
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(27) ! y**3*z
       alpha=2
       beta=2
       gamma=2
       delta=3
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(28) ! z**3*x
       alpha=3
       beta=3
       gamma=3
       delta=1
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(29) ! z**3*y
       alpha=3
       beta=3
       gamma=3
       delta=2
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(30) ! x**2*y**2
       alpha=1
       beta=1
       gamma=2
       delta=2
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(31) ! x**2*z**2
       alpha=1
       beta=1
       gamma=3
       delta=3
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(32) ! y**2*z**2
       alpha=2
       beta=2
       gamma=3
       delta=3
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(33) ! x**2*y*z
       alpha=1
       beta=1
       gamma=2
       delta=3
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(34) ! y**2*x*z
       alpha=2
       beta=2
       gamma=1
       delta=3
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    case(35) ! z**2*x*y
       alpha=3
       beta=3
       gamma=1
       delta=2
       call  integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_poly_coef)
    end select
    vol_poly_coeff(c,i)=vol_poly_coef
 enddo


end subroutine vol_poly_coef_hexaedron







! This subroutine calculates the surface and gradient averaged valaues multiplying each polynomial 
! coeffcient for a tetraedron cell by using the following polynomial approach: 
! phi(x,y,z)= a + b*x + c*y + d*z + e*x**2 + f*y**2 + g*z**2 + h*x*y + k*x*z + l*y*z + m*x**3 + n*y**3 + o*z**3 + 
!           + p*x**2*y + q*x**2*z + r*y**2*x + s*y**2*z + t*z**2*x + u*z**2*y + v*x*y*z
! But only with 5 polynomial terms
subroutine surf_and_grad_poly_coef_tetraedron(c,vect_coef)

use general_module


 ! Global variables
 integer c, vect_coef(5)

 ! Local variables
 integer j, jcell, i, k
 double precision Xg,Yg,Zg, Jx, Jy,Jz, Pxy,Pxz, Pyz, &
 & coord_node1(3),coord_node2(3),coord_node3(3), surf_poly_coef, grad_poly_coef, &
 & ux, uy, uz, x0(3,3), surf1, surf2, surf3
 integer alpha, beta, gamma, delta
    

 ! Calculate surface averaged values multipliers
 do jcell=1,ubound(cell(c)%jface,1)
    j=cell(c)%jface(jcell)

    ! Calculation of centroid and inertia tensor for each surface
    Xg=face(j)%x(1)
    Yg=face(j)%x(2)
    Zg=face(j)%x(3)
    coord_node1=node(face(j)%knode(1))%x
    coord_node2=node(face(j)%knode(2))%x
    coord_node3=node(face(j)%knode(3))%x
    call inertia_tensor_triangle(coord_node1,coord_node2,coord_node3, &
         & Jx, Jy, Jz, Pxy, Pxz, Pyz)

    ! Determination of the correct sense of surface normal vector
    ux=face(j)%norm(1,1)
    uy=face(j)%norm(2,1)
    uz=face(j)%norm(3,1)
    if (c/=face(j)%icell(1)) then
       ux=(-1)*ux
       uy=(-1)*uy
       uz=(-1)*uz
    endif

    ! Definition of triangle coordinates in global coordinates
    do i=1,3
       do k=1,3
          x0(i,k)=node(face(j)%knode(k))%x(i)
       enddo
    enddo


    ! Determination of surface and gradient multipliers for each surface
    do i=1,5
       select case (vect_coef(i))
       case (1) ! a : 1
          surf_poly_coef=1
          grad_poly_coef=0
       case (2) ! b*x 
          surf_poly_coef=Xg
          grad_poly_coef=ux
       case (3) ! c*y 
          surf_poly_coef=Yg
          grad_poly_coef=uy
       case (4) ! d*z 
          surf_poly_coef=Zg
          grad_poly_coef=uz
       case (5) ! e*x**2 
          surf_poly_coef=(Jy+Jz-Jx)/(2*face(j)%area)
          grad_poly_coef=2*ux*Xg
       case (6) ! f*y**2 
          surf_poly_coef=(Jx+Jz-Jy)/(2*face(j)%area)
          grad_poly_coef=2*uy*Yg
       case (7) ! g*z**2 
          surf_poly_coef=(Jx+Jy-Jz)/(2*face(j)%area)
          grad_poly_coef=2*uz*Zg
       case (8) ! h*x*y 
          surf_poly_coef=Pxy/face(j)%area
          grad_poly_coef=ux*Yg+uy*Xg
       case (9) ! k*x*z 
          surf_poly_coef=Pxz/face(j)%area
          grad_poly_coef=ux*Zg+uz*Xg
       case (10) ! l*y*z 
          surf_poly_coef=Pyz/face(j)%area
          grad_poly_coef=uy*Zg+uz*Yg
       case(11) ! m*x**3
          alpha=1
          beta=1
          gamma=1
          call  integrate_order3_area_triangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=ux*3*(Jy+Jz-Jx)/(2*face(j)%area)
       case(12) ! n*y**3
          alpha=2
          beta=2
          gamma=2
          call  integrate_order3_area_triangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=uy*3*(Jx+Jz-Jy)/(2*face(j)%area)
       case(13) ! o*z**3
          alpha=3
          beta=3
          gamma=3
          call  integrate_order3_area_triangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=uz*3*(Jx+Jy-Jz)/(2*face(j)%area)
       case(14) ! p*x**2*y
          alpha=1
          beta=1
          gamma=2
          call  integrate_order3_area_triangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=ux*2*Pxy/face(j)%area+uy*(Jy+Jz-Jx)/(2*face(j)%area)
       case(15) ! q*x**2*z
          alpha=1
          beta=1
          gamma=3
          call  integrate_order3_area_triangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=ux*2*Pxz/face(j)%area+uz*(Jy+Jz-Jx)/(2*face(j)%area)
       case(16) ! r*y**2*x
          alpha=2
          beta=2
          gamma=1
          call  integrate_order3_area_triangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=ux*(Jx+Jz-Jy)/(2*face(j)%area)+uy*2*Pxy/face(j)%area
       case(17) ! s*y**2*z
          alpha=2
          beta=2
          gamma=3
          call  integrate_order3_area_triangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=uy*2*Pyz/face(j)%area+uz*(Jx+Jz-Jy)/(2*face(j)%area)
       case(18) ! t*z**2*x
          alpha=3
          beta=3
          gamma=1
          call  integrate_order3_area_triangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=ux*(Jx+Jy-Jz)/(2*face(j)%area)+uz*2*Pxz/face(j)%area
       case(19) ! u*z**2*y
          alpha=3
          beta=3
          gamma=2
          call  integrate_order3_area_triangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=uy*(Jx+Jy-Jz)/(2*face(j)%area)+uz*2*Pyz/face(j)%area
       case(20) ! v*x*y*z
          alpha=1
          beta=2
          gamma=3
          call  integrate_order3_area_triangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=ux*Pyz/face(j)%area+uy*Pxz/face(j)%area+uz*Pxy/face(j)%area
       case(21) ! x**4
          alpha=1
          beta=1
          gamma=1
          delta=1
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,alpha,beta,gamma,surf1)
          grad_poly_coef=ux*4*surf1
       case(22) ! y**4
          alpha=2
          beta=2
          gamma=2
          delta=2
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,alpha,beta,gamma,surf1)
          grad_poly_coef=uy*4*surf1
       case(23) ! z**4
          alpha=3
          beta=3
          gamma=3
          delta=3
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,alpha,beta,gamma,surf1)
          grad_poly_coef=uz*4*surf1
       case(24) ! x**3*y
          alpha=1
          beta=1
          gamma=1
          delta=2
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,1,1,2,surf1)
          call  integrate_order3_area_triangle_coordinates(j,x0,1,1,1,surf2)
          grad_poly_coef=ux*3*surf1+uy*surf2
       case(25) ! x**3*z
          alpha=1
          beta=1
          gamma=1
          delta=3
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,1,1,3,surf1)
          call  integrate_order3_area_triangle_coordinates(j,x0,1,1,1,surf2)
          grad_poly_coef=ux*3*surf1+uz*surf2
       case(26) ! y**3*x
          alpha=2
          beta=2
          gamma=2
          delta=1
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,2,2,1,surf1)
          call  integrate_order3_area_triangle_coordinates(j,x0,2,2,2,surf2)
          grad_poly_coef=uy*3*surf1+ux*surf2
       case(27) ! y**3*z
          alpha=2
          beta=2
          gamma=2
          delta=3
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,2,2,3,surf1)
          call  integrate_order3_area_triangle_coordinates(j,x0,2,2,2,surf2)
          grad_poly_coef=uy*3*surf1+uz*surf2
       case(28) ! z**3*x
          alpha=3
          beta=3
          gamma=3
          delta=1
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,3,3,1,surf1)
          call  integrate_order3_area_triangle_coordinates(j,x0,3,3,3,surf2)
          grad_poly_coef=uz*3*surf1+ux*surf2
       case(29) ! z**3*y
          alpha=3
          beta=3
          gamma=3
          delta=2
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,3,3,2,surf1)
          call  integrate_order3_area_triangle_coordinates(j,x0,3,3,3,surf2)
          grad_poly_coef=uz*3*surf1+uy*surf2
       case(30) ! x**2*y**2
          alpha=1
          beta=1
          gamma=2
          delta=2
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,1,2,2,surf1)
          call  integrate_order3_area_triangle_coordinates(j,x0,2,1,1,surf2)
          grad_poly_coef=ux*2*surf1+uy*2*surf2
       case(31) ! x**2*z**2
          alpha=1
          beta=1
          gamma=3
          delta=3
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,1,3,3,surf1)
          call  integrate_order3_area_triangle_coordinates(j,x0,3,1,1,surf2)
          grad_poly_coef=ux*2*surf1+uz*2*surf2
       case(32) ! y**2*z**2
          alpha=2
          beta=2
          gamma=3
          delta=3
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,2,3,3,surf1)
          call  integrate_order3_area_triangle_coordinates(j,x0,3,2,2,surf2)
          grad_poly_coef=uy*2*surf1+uz*2*surf2
       case(33) ! x**2*y*z
          alpha=1
          beta=1
          gamma=2
          delta=3
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,1,2,3,surf1)
          call  integrate_order3_area_triangle_coordinates(j,x0,1,1,3,surf2)
          call  integrate_order3_area_triangle_coordinates(j,x0,1,1,2,surf3)
          grad_poly_coef=ux*2*surf1+uy*surf2+uz*surf3
       case(34) ! y**2*x*z
          alpha=2
          beta=2
          gamma=1
          delta=3
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,1,2,3,surf1)
          call  integrate_order3_area_triangle_coordinates(j,x0,2,2,3,surf2)
          call  integrate_order3_area_triangle_coordinates(j,x0,2,2,1,surf3)
          grad_poly_coef=uy*2*surf1+ux*surf2+uz*surf3
       case(35) ! z**2*x*y
          alpha=3
          beta=3
          gamma=1
          delta=2
          call  integrate_order4_area_triangle_coordinates(j,x0,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_triangle_coordinates(j,x0,1,2,3,surf1)
          call  integrate_order3_area_triangle_coordinates(j,x0,3,3,2,surf2)
          call  integrate_order3_area_triangle_coordinates(j,x0,3,3,1,surf3)
          grad_poly_coef=uz*2*surf1+ux*surf2+uy*surf3
       end select
       surf_poly_coeff(c,jcell,i)=surf_poly_coef
       grad_poly_coeff(c,jcell,i)=grad_poly_coef
    enddo
 enddo


end subroutine surf_and_grad_poly_coef_tetraedron







! This subroutine calculates the surface and gradient averaged valaues multiplying each polynomial 
! coeffcient for a hexaedron cell by using the following polynomial approach: 
! phi(x,y,z)= a + b*x + c*y + d*z + e*x**2 + f*y**2 + g*z**2 + h*x*y + k*x*z + l*y*z + m*x**3 + n*y**3 + o*z**3 + 
!           + p*x**2*y + q*x**2*z + r*y**2*x + s*y**2*z + t*z**2*x + u*z**2*y + v*x*y*z
! But only with 7 polynomial terms
subroutine surf_and_grad_poly_coef_hexaedron(c,vect_coef)

use general_module


 ! Global variables
 integer c, vect_coef(7)

 ! Local variables
 integer j, jcell, i
 double precision Xg,Yg,Zg, Jx, Jy,Jz, Pxy,Pxz, Pyz, &
 & coord_node1(3),coord_node2(3),coord_node3(3), surf_poly_coef, grad_poly_coef, &
 & ux, uy, uz, coord_node4(3), surf1, surf2, surf3
 integer alpha, beta, gamma, delta
    

 ! Calculate surface averaged values multipliers
 do jcell=1,ubound(cell(c)%jface,1)
    j=cell(c)%jface(jcell)

    ! Calculation of centroid and inertia tensor for each surface
    Xg=face(j)%x(1)
    Yg=face(j)%x(2)
    Zg=face(j)%x(3)
    coord_node1=node(face(j)%knode(1))%x
    coord_node2=node(face(j)%knode(2))%x
    coord_node3=node(face(j)%knode(3))%x
    coord_node4=node(face(j)%knode(4))%x
    call inertia_tensor_cuadrangle(coord_node1,coord_node2,coord_node3, &
         & coord_node4, Jx, Jy, Jz, Pxy, Pxz, Pyz)

    ! Determination of the correct sense of surface normal vector
    ux=face(j)%norm(1,1)
    uy=face(j)%norm(2,1)
    uz=face(j)%norm(3,1)
    if (c/=face(j)%icell(1)) then
       ux=(-1)*ux
       uy=(-1)*uy
       uz=(-1)*uz
    endif

    ! Determination of surface and gradient multipliers for each surface
    do i=1,7
       select case (vect_coef(i))
       case (1) ! a : 1
          surf_poly_coef=1
          grad_poly_coef=0
       case (2) ! b*x 
          surf_poly_coef=Xg
          grad_poly_coef=ux
       case (3) ! c*y 
          surf_poly_coef=Yg
          grad_poly_coef=uy
       case (4) ! d*z 
          surf_poly_coef=Zg
          grad_poly_coef=uz
       case (5) ! e*x**2 
          surf_poly_coef=(Jy+Jz-Jx)/(2*face(j)%area)
          grad_poly_coef=2*ux*Xg
       case (6) ! f*y**2 
          surf_poly_coef=(Jx+Jz-Jy)/(2*face(j)%area)
          grad_poly_coef=2*uy*Yg
       case (7) ! g*z**2 
          surf_poly_coef=(Jx+Jy-Jz)/(2*face(j)%area)
          grad_poly_coef=2*uz*Zg
       case (8) ! h*x*y 
          surf_poly_coef=Pxy/face(j)%area
          grad_poly_coef=ux*Yg+uy*Xg
       case (9) ! k*x*z 
          surf_poly_coef=Pxz/face(j)%area
          grad_poly_coef=ux*Zg+uz*Xg
       case (10) ! l*y*z 
          surf_poly_coef=Pyz/face(j)%area
          grad_poly_coef=uy*Zg+uz*Yg
       case(11) ! m*x**3
          alpha=1
          beta=1
          gamma=1
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=ux*3*(Jy+Jz-Jx)/(2*face(j)%area)
       case(12) ! n*y**3
          alpha=2
          beta=2
          gamma=2
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=uy*3*(Jx+Jz-Jy)/(2*face(j)%area)
       case(13) ! o*z**3
          alpha=3
          beta=3
          gamma=3
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=uz*3*(Jx+Jy-Jz)/(2*face(j)%area)
       case(14) ! p*x**2*y
          alpha=1
          beta=1
          gamma=2
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=ux*2*Pxy/face(j)%area+uy*(Jy+Jz-Jx)/(2*face(j)%area)
       case(15) ! q*x**2*z
          alpha=1
          beta=1
          gamma=3
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=ux*2*Pxz/face(j)%area+uz*(Jy+Jz-Jx)/(2*face(j)%area)
       case(16) ! r*y**2*x
          alpha=2
          beta=2
          gamma=1
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=ux*(Jx+Jz-Jy)/(2*face(j)%area)+uy*2*Pxy/face(j)%area
       case(17) ! s*y**2*z
          alpha=2
          beta=2
          gamma=3
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=uy*2*Pyz/face(j)%area+uz*(Jx+Jz-Jy)/(2*face(j)%area)
       case(18) ! t*z**2*x
          alpha=3
          beta=3
          gamma=1
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=ux*(Jx+Jy-Jz)/(2*face(j)%area)+uz*2*Pxz/face(j)%area
       case(19) ! u*z**2*y
          alpha=3
          beta=3
          gamma=2
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=uy*(Jx+Jy-Jz)/(2*face(j)%area)+uz*2*Pyz/face(j)%area
       case(20) ! v*x*y*z
          alpha=1
          beta=2
          gamma=3
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf_poly_coef)
          grad_poly_coef=ux*Pyz/face(j)%area+uy*Pxz/face(j)%area+uz*Pxy/face(j)%area
       case(21) ! x**4
          alpha=1
          beta=1
          gamma=1
          delta=1
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf1)
          grad_poly_coef=ux*4*surf1
       case(22) ! y**4
          alpha=2
          beta=2
          gamma=2
          delta=2
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf1)
          grad_poly_coef=uy*4*surf1
       case(23) ! z**4
          alpha=3
          beta=3
          gamma=3
          delta=3
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,alpha,beta,gamma,surf1)
          grad_poly_coef=uz*4*surf1
       case(24) ! x**3*y
          alpha=1
          beta=1
          gamma=1
          delta=2
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,1,1,2,surf1)
          call  integrate_order3_area_cuadrangle(j,1,1,1,surf2)
          grad_poly_coef=ux*3*surf1+uy*surf2
       case(25) ! x**3*z
          alpha=1
          beta=1
          gamma=1
          delta=3
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,1,1,3,surf1)
          call  integrate_order3_area_cuadrangle(j,1,1,1,surf2)
          grad_poly_coef=ux*3*surf1+uz*surf2
       case(26) ! y**3*x
          alpha=2
          beta=2
          gamma=2
          delta=1
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,2,2,1,surf1)
          call  integrate_order3_area_cuadrangle(j,2,2,2,surf2)
          grad_poly_coef=uy*3*surf1+ux*surf2
       case(27) ! y**3*z
          alpha=2
          beta=2
          gamma=2
          delta=3
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,2,2,3,surf1)
          call  integrate_order3_area_cuadrangle(j,2,2,2,surf2)
          grad_poly_coef=uy*3*surf1+uz*surf2
       case(28) ! z**3*x
          alpha=3
          beta=3
          gamma=3
          delta=1
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,3,3,1,surf1)
          call  integrate_order3_area_cuadrangle(j,3,3,3,surf2)
          grad_poly_coef=uz*3*surf1+ux*surf2
       case(29) ! z**3*y
          alpha=3
          beta=3
          gamma=3
          delta=2
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,3,3,2,surf1)
          call  integrate_order3_area_cuadrangle(j,3,3,3,surf2)
          grad_poly_coef=uz*3*surf1+uy*surf2
       case(30) ! x**2*y**2
          alpha=1
          beta=1
          gamma=2
          delta=2
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,1,2,2,surf1)
          call  integrate_order3_area_cuadrangle(j,2,1,1,surf2)
          grad_poly_coef=ux*2*surf1+uy*2*surf2
       case(31) ! x**2*z**2
          alpha=1
          beta=1
          gamma=3
          delta=3
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,1,3,3,surf1)
          call  integrate_order3_area_cuadrangle(j,3,1,1,surf2)
          grad_poly_coef=ux*2*surf1+uz*2*surf2
       case(32) ! y**2*z**2
          alpha=2
          beta=2
          gamma=3
          delta=3
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,2,3,3,surf1)
          call  integrate_order3_area_cuadrangle(j,3,2,2,surf2)
          grad_poly_coef=uy*2*surf1+uz*2*surf2
       case(33) ! x**2*y*z
          alpha=1
          beta=1
          gamma=2
          delta=3
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,1,2,3,surf1)
          call  integrate_order3_area_cuadrangle(j,1,1,3,surf2)
          call  integrate_order3_area_cuadrangle(j,1,1,2,surf3)
          grad_poly_coef=ux*2*surf1+uy*surf2+uz*surf3
       case(34) ! y**2*x*z
          alpha=2
          beta=2
          gamma=1
          delta=3
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,1,2,3,surf1)
          call  integrate_order3_area_cuadrangle(j,2,2,3,surf2)
          call  integrate_order3_area_cuadrangle(j,2,2,1,surf3)
          grad_poly_coef=uy*2*surf1+ux*surf2+uz*surf3
       case(35) ! z**2*x*y
          alpha=3
          beta=3
          gamma=1
          delta=2
          call  integrate_order4_area_cuadrangle(j,alpha,beta,gamma,delta,surf_poly_coef)
          call  integrate_order3_area_cuadrangle(j,1,2,3,surf1)
          call  integrate_order3_area_cuadrangle(j,3,3,2,surf2)
          call  integrate_order3_area_cuadrangle(j,3,3,1,surf3)
          grad_poly_coef=uz*2*surf1+ux*surf2+uy*surf3
       end select
       surf_poly_coeff(c,jcell,i)=surf_poly_coef
       grad_poly_coeff(c,jcell,i)=grad_poly_coef
    enddo
 enddo


end subroutine surf_and_grad_poly_coef_hexaedron











! This subroutine calculates the inertia tensor of a triangle 
subroutine inertia_tensor_triangle(coord_node1,coord_node2,coord_node3, &
 & Jx_t, Jy_t, Jz_t, Pxy_t, Pxz_t, Pyz_t)

use general_module


 ! Global variables
 double precision, dimension(:), intent(in) :: coord_node1,coord_node2,coord_node3
 double precision, intent(out) :: Jx_t, Jy_t, Jz_t, Pxy_t, Pxz_t, Pyz_t

 ! Local variables
 double precision r12,r13, r23, a, h, b, Jx, Jy,Jz, Pxy,Pxz, Pyz, &
 & Pxy_r, Pxz_r, Pyz_r, Jx_r, Jy_r, Jz_r
 double precision, dimension(3) :: vector12, vector12y, vector12z, vector13 
 double precision, dimension(3) :: ux,uy,uz
 integer k,l
 double precision coord_centroid(3), normtriangle(3), area, rotation_matrix(3,3), &
   & inertia_tensor(3,3)

 ! Calculation of Product of Inertia and Moment of Inertia: Centroid and X axis is line joining nodes 1 and 2
 r12=distance(coord_node1,coord_node2)
 r13=distance(coord_node1,coord_node3)
 r23=distance(coord_node2,coord_node3)
 a=(r12**2+r13**2-r23**2)/(2*r12)
 h=sqrt(r13**2-a**2)
 b=r12
 Jx=b*h**3/36.0
 Jy=b*h/36.0*(b**2-a*(b-a))
 Jz=Jx+Jy
 Pxy=h**2*b*(2*a-b)/72.0
 Pxz=0
 Pyz=0
 inertia_tensor(1,:)=(/Jx,-1*Pxy,-1*Pxz/)
 inertia_tensor(2,:)=(/-1*Pxy,Jy,-1*Pyz/)
 inertia_tensor(3,:)=(/-1*Pxz,-1*Pyz,Jz/)
 


 ! Calculation of the rotation matrix
 ! Definition of the global axes vectors
 ux(:)=(/1.0 , 0.0 , 0.0/)
 uy(:)=(/0.0 , 1.0 , 0.0/)
 uz(:)=(/0.0 , 0.0 , 1.0/)
 ! Definition of the local axes vectors
 do k=1,3
    vector12(k)=coord_node2(k)-coord_node1(k)
    vector13(k)=coord_node3(k)-coord_node1(k)
 enddo
 vector12z=cross_product(vector12,vector13)
 vector12y=cross_product(vector12z,vector12)
 call normalise_vector(vector12)
 call normalise_vector(vector12y)
 call normalise_vector(vector12z)
 ! Calculation of the cosine of the angles between global and local axes
 rotation_matrix(1,1)=dot_product(ux,vector12)
 rotation_matrix(1,2)=dot_product(ux,vector12y)
 rotation_matrix(1,3)=dot_product(ux,vector12z)
 rotation_matrix(2,1)=dot_product(uy,vector12)
 rotation_matrix(2,2)=dot_product(uy,vector12y)
 rotation_matrix(2,3)=dot_product(uy,vector12z)
 rotation_matrix(3,1)=dot_product(uz,vector12)
 rotation_matrix(3,2)=dot_product(uz,vector12y)
 rotation_matrix(3,3)=dot_product(uz,vector12z)



 ! Calculation of Product of Inertia and Moment of Inertia: Rotation of axes
 Jx_r=0
 Jy_r=0
 Jz_r=0
 Pxy_r=0
 Pxz_r=0
 Pyz_r=0
 do k=1,3
    do l=1,3
       Jx_r=Jx_r+rotation_matrix(1,l)*inertia_tensor(l,k)*rotation_matrix(1,k)
       Jy_r=Jy_r+rotation_matrix(2,l)*inertia_tensor(l,k)*rotation_matrix(2,k)
       Jz_r=Jz_r+rotation_matrix(3,l)*inertia_tensor(l,k)*rotation_matrix(3,k)
       Pxy_r=Pxy_r+rotation_matrix(1,l)*inertia_tensor(l,k)*rotation_matrix(2,k)
       Pxz_r=Pxz_r+rotation_matrix(1,l)*inertia_tensor(l,k)*rotation_matrix(3,k)
       Pyz_r=Pyz_r+rotation_matrix(2,l)*inertia_tensor(l,k)*rotation_matrix(3,k)
    enddo
 enddo
 Pxy_r=-1*Pxy_r
 Pxz_r=-1*Pxz_r
 Pyz_r=-1*Pyz_r
 

 ! Calculation of Product of Inertia and Moment of Inertia: Translation of axes
 coord_centroid(1)=(coord_node1(1)+coord_node2(1)+coord_node3(1))/3.0
 coord_centroid(2)=(coord_node1(2)+coord_node2(2)+coord_node3(2))/3.0
 coord_centroid(3)=(coord_node1(3)+coord_node2(3)+coord_node3(3))/3.0
 normtriangle = cross_product(coord_node2-coord_node1,coord_node3-coord_node1)
 area = vector_magnitude(normtriangle)/2.0
 Jx_t=Jx_r+area*(coord_centroid(2)**2+coord_centroid(3)**2)
 Jy_t=Jy_r+area*(coord_centroid(1)**2+coord_centroid(3)**2)
 Jz_t=Jz_r+area*(coord_centroid(1)**2+coord_centroid(2)**2)
 Pxy_t=Pxy_r+area*coord_centroid(1)*coord_centroid(2)
 Pxz_t=Pxz_r+area*coord_centroid(1)*coord_centroid(3)
 Pyz_t=Pyz_r+area*coord_centroid(2)*coord_centroid(3)
 

end subroutine inertia_tensor_triangle











! This subroutine calculates the inertia tensor of a cuadrangle 
subroutine inertia_tensor_cuadrangle(coord_node1,coord_node2,coord_node3, &
         & coord_node4, Jx, Jy, Jz, Pxy, Pxz, Pyz)
 
use general_module


 ! Global variables
 double precision, dimension(:), intent(in) :: coord_node1,coord_node2,coord_node3, coord_node4
 double precision, intent(out) :: Jx, Jy, Jz, Pxy, Pxz, Pyz

 ! Local variables
 double precision Jx_1, Jy_1, Jz_1, Pxy_1, Pxz_1, Pyz_1, &
  & Jx_2, Jy_2, Jz_2, Pxy_2, Pxz_2, Pyz_2


 ! The inertia tensor will be calculated as a sum of 2 triangles
 ! Split the cuadrangle into 2 triangles: First one is composed of nodes 1, 2, 4 and the other is composed of
 ! nodes 2, 3, 4. 

 ! Determination of Moments of Intertia of the first triangle
 call inertia_tensor_triangle(coord_node1,coord_node2,coord_node4, &
 & Jx_1, Jy_1, Jz_1, Pxy_1, Pxz_1, Pyz_1)

 ! Determination of Moments of Intertia of the second triangle
 call inertia_tensor_triangle(coord_node2,coord_node3,coord_node4, &
 & Jx_2, Jy_2, Jz_2, Pxy_2, Pxz_2, Pyz_2)

 ! The inertia tensor will be calculated as a sum of 2 triangles
 Jx=Jx_1+Jx_2
 Jy=Jy_1+Jy_2
 Jz=Jz_1+Jz_2
 Pxy=Pxy_1+Pxy_2
 Pxz=Pxz_1+Pxz_2
 Pyz=Pyz_1+Pyz_2


end subroutine inertia_tensor_cuadrangle
















! This subroutine calculates the analytical integral of a third order polynomial term for a tetraedron cell
! This term is expressed as: alpha*beta*gamma;   The possible values for alpha,beta and gamma are: 1=x, 2=y, 3=z
subroutine integrate_order3_volume_tetraedron(c,alpha,beta,gamma,vol_int)

use general_module


 ! Global variables
 integer, intent(in) :: c, alpha, beta, gamma
 double precision, intent(out) :: vol_int

 ! Local variables
 double precision x(3,4), a1, a2, a3, &
 & integral_tau3, integral_tau2_x_tau1, integral_tau_x_tau_x_tau
 integer i,j,k


 ! Definition of tetrahedron coordinates
 do i=1,3
    do j=1,4
       x(i,j)=node(cell(c)%knode(j))%x(i)
    enddo
 enddo



 ! Calculation of: "integral_tau3", "integral_tau2_x_tau1",  "integral_tau_x_tau_x_tau"
 integral_tau3= 1.0/20.0
 integral_tau2_x_tau1= 1.0/60.0
 integral_tau_x_tau_x_tau= 1.0/120.0


 ! Calculation of "a1"
 a1=0.0
 do i=1,4
    a1=a1 + x(alpha,i)*x(beta,i)*x(gamma,i)
 enddo


 ! Calculation of "a2"
 a2=0.0
 do i=1,4
    do j=1,4
       if (j/=i) then
          a2=a2 + x(alpha,i)*x(beta,i)*x(gamma,j) + x(alpha,i)*x(gamma,i)*x(beta,j) + x(beta,i)*x(gamma,i)*x(alpha,j)
       endif
    enddo
 enddo


 ! Calculation of "a3"
 a3=0.0
 do i=1,4
    do j=1,4
       if (j/=i) then
          do k=1,4
             if ((k/=j).and.(k/=i)) then
                a3=a3 + x(alpha,i)*x(beta,j)*x(gamma,k)
             endif
          enddo
       endif
    enddo
 enddo  


 ! Calculation of the integral
 vol_int= a1*integral_tau3 + a2*integral_tau2_x_tau1 + a3*integral_tau_x_tau_x_tau

end subroutine integrate_order3_volume_tetraedron















! This subroutine calculates the analytical integral of a third order polynomial term for a triangle face
! This term is expressed as: alpha*beta*gamma;   The possible values for alpha,beta and gamma are: 1=x, 2=y, 3=z
subroutine integrate_order3_area_triangle(facej,alpha,beta,gamma,area_int)

use general_module


 ! Global variables
 integer, intent(in) :: facej,alpha,beta,gamma
 double precision, intent(out) :: area_int

 ! Local variables
 double precision x0(3,3), x(3), y(3), z
 double precision a1, a2, a3, a4, &
 & a5, a6, a7, a8, a9, a10, b1, b2, b3, c1, c2, c3, &
 & d1, d2, d3, e1, e2, e3, f1, f2, g1, g2, h1, h2, m ,n , &
 & p1, p2, p3, p4, p5, p6, p7
 double precision integral_tau3, integral_tau2_x_tau1, integral_tau_x_tau_x_tau, &
 & integral_tau2, integral_tau_x_tau, integral_tau
 double precision cos_x0_x, cos_x0_y, cos_y0_x, cos_y0_y, &
 & cos_z0_x, cos_z0_y, cos_x0_z, cos_y0_z, cos_z0_z, rot(3,3)
 double precision, dimension(3) :: ux, uy, uz
 integer i,j,k


 ! Definition of triangle coordinates in global coordinates
 do i=1,3
    do j=1,3
       x0(i,j)=node(face(facej)%knode(j))%x(i)
    enddo
 enddo



 ! Definition of the global axes vectors
 ux(:)=(/1.0 , 0.0 , 0.0/)
 uy(:)=(/0.0 , 1.0 , 0.0/)
 uz(:)=(/0.0 , 0.0 , 1.0/)


 ! Calculation of angles between local and global axes
 cos_x0_x=dot_product(ux,face(facej)%norm(:,2))
 cos_x0_y=dot_product(ux,face(facej)%norm(:,3))
 cos_x0_z=dot_product(ux,face(facej)%norm(:,1))
 cos_y0_x=dot_product(uy,face(facej)%norm(:,2))
 cos_y0_y=dot_product(uy,face(facej)%norm(:,3))
 cos_y0_z=dot_product(uy,face(facej)%norm(:,1))
 cos_z0_x=dot_product(uz,face(facej)%norm(:,2))
 cos_z0_y=dot_product(uz,face(facej)%norm(:,3))
 cos_z0_z=dot_product(uz,face(facej)%norm(:,1))


 ! Definition of the rotation matrix
 rot(1,1:3)=(/cos_x0_x, cos_x0_y, cos_x0_z/);
 rot(2,1:3)=(/cos_y0_x, cos_y0_y, cos_y0_z/);
 rot(3,1:3)=(/cos_z0_x, cos_z0_y, cos_z0_z/);


 ! Calculation of: "integral_tau3", "integral_tau2_x_tau1",  "integral_tau_x_tau_x_tau", 
 ! "integral_tau2", "integral_tau_x_tau", "integral_tau"
 integral_tau3= 1.0/10.0
 integral_tau2_x_tau1= 1.0/30.0
 integral_tau_x_tau_x_tau= 1.0/60.0
 integral_tau2=1.0/6.0
 integral_tau_x_tau=1.0/12.0
 integral_tau=1.0/3.0


 ! Calculation of "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10" 
 a1=rot(alpha,1)*rot(beta,1)*rot(gamma,1)
 a2=rot(alpha,2)*rot(beta,2)*rot(gamma,2)
 a3=rot(alpha,3)*rot(beta,3)*rot(gamma,3)
 a4=rot(alpha,1)*rot(beta,1)*rot(gamma,2) + rot(alpha,1)*rot(gamma,1)*rot(beta,2) + rot(beta,1)*rot(gamma,1)*rot(alpha,2)
 a5=rot(alpha,1)*rot(beta,1)*rot(gamma,3) + rot(alpha,1)*rot(gamma,1)*rot(beta,3) + rot(beta,1)*rot(gamma,1)*rot(alpha,3)
 a6=rot(alpha,2)*rot(beta,2)*rot(gamma,1) + rot(alpha,2)*rot(gamma,2)*rot(beta,1) + rot(beta,2)*rot(gamma,2)*rot(alpha,1)
 a7=rot(alpha,2)*rot(beta,2)*rot(gamma,3) + rot(alpha,2)*rot(gamma,2)*rot(beta,3) + rot(beta,2)*rot(gamma,2)*rot(alpha,3)
 a8=rot(alpha,3)*rot(beta,3)*rot(gamma,1) + rot(alpha,3)*rot(gamma,3)*rot(beta,1) + rot(beta,3)*rot(gamma,3)*rot(alpha,1)
 a9=rot(alpha,3)*rot(beta,3)*rot(gamma,2) + rot(alpha,3)*rot(gamma,3)*rot(beta,2) + rot(beta,3)*rot(gamma,3)*rot(alpha,2)
 a10=rot(alpha,1)*(rot(beta,2)*rot(gamma,3) + rot(gamma,2)*rot(beta,3)) + &
   & rot(beta,1)*(rot(alpha,2)*rot(gamma,3) + rot(gamma,2)*rot(alpha,3)) + &
   & rot(gamma,1)*(rot(alpha,2)*rot(beta,3) + rot(beta,2)*rot(alpha,3))

 
 ! Calculation of local coordinates of triangle coordinates
 do i=1,3
    x(i)=cos_x0_x*x0(1,i) + cos_y0_x*x0(2,i) + cos_z0_x*x0(3,i)
    y(i)=cos_x0_y*x0(1,i) + cos_y0_y*x0(2,i) + cos_z0_y*x0(3,i)
 enddo
 z=cos_x0_z*x0(1,1) + cos_y0_z*x0(2,1) + cos_z0_z*x0(3,1)


 ! Calculation of "b1", "c1", "d1", "e1", "f1", "g1", "h1", "m", "n"
 b1=0.0
 c1=0.0
 d1=0.0
 e1=0.0
 f1=0.0
 g1=0.0
 h1=0.0
 m=0.0
 n=0.0
 do i=1,3
    b1=b1 + x(i)**3
    c1=c1 + y(i)**3
    d1=d1 + x(i)**2*y(i)
    e1=e1 + x(i)*y(i)**2
    f1=f1 + x(i)**2
    g1=g1 + y(i)**2
    h1=h1 + x(i)*y(i)
    m=m + x(i)
    n=n + y(i)
 enddo


 ! Calculation of "b2", "c2", "d2", "e2", "f2", "g2", "h2"
 b2=0.0
 c2=0.0
 d2=0.0
 e2=0.0
 f2=0.0
 g2=0.0
 h2=0.0
 do i=1,3
    do j=1,3
       if (j/=i) then
          b2=b2 + 3*x(i)**2*x(j)
          c2=c2 + 3*y(i)**2*y(j)
          d2=d2 + x(i)**2*y(j) + 2*x(i)*y(i)*x(j)
          e2=e2 + y(i)**2*x(j) + 2*x(i)*y(i)*y(j)
          f2=f2 + x(i)*x(j)
          g2=g2 + y(i)*y(j)
          h2=h2 + x(i)*y(j)
       endif
    enddo
 enddo


 ! Calculation of "b3", "c3", "d3", "e3"
 b3=0.0
 c3=0.0
 d3=0.0
 e3=0.0
 do i=1,3
    do j=1,3
       if (j/=i) then
          do k=1,3
             if ((k/=j).and.(k/=i)) then
                b3=b3 + x(i)*x(j)*x(k)
                c3=c3 + y(i)*y(j)*y(k)
                d3=d3 + x(i)*x(j)*y(k)
                e3=e3 + y(i)*y(j)*x(k)
             endif
          enddo
       endif
    enddo
 enddo  


 
 ! Calculaltion of "p1", "p2", "p3", "p4", "p5", "p6", "p7"
 p1= a1*b1 + a2*c1 + a4*d1 + a6*e1
 p2= a1*b2 + a2*c2 + a4*d2 + a6*e2
 p3= a1*b3 + a2*c3 + a4*d3 + a6*e3 
 p4= a5*z*f1 + a7*z*g1 + a10*z*h1
 p5= a5*z*f2 + a7*z*g2 + a10*z*h2
 p6= a8*z**2*m + a9*z**2*n
 p7= a3*z**3


 ! Calculation of the integral
 area_int= p1*integral_tau3 + p2*integral_tau2_x_tau1 + p3*integral_tau_x_tau_x_tau + &
         & p4*integral_tau2 + p5*integral_tau_x_tau + p6*integral_tau + p7

end subroutine integrate_order3_area_triangle
















! This subroutine calculates the analytical integral of a third order polynomial term for a hexaedron cell
! This term is expressed as: alpha*beta*gamma;   The possible values for alpha,beta and gamma are: 1=x, 2=y, 3=z
subroutine integrate_order3_volume_hexaedron(c,alpha,beta,gamma,vol_int)

use general_module


 ! Global variables
 integer, intent(in) :: c, alpha, beta, gamma
 double precision, intent(out) :: vol_int

 ! Local variables
 double precision x(3,4), vol_tetraedron, aux_vol_int
 integer i,j,k, nodes_tetraedron(5,4), knode(4)


 ! The hexaedron is divided into 5 tetraedron and the integral is calculated by adding the integral of each tetraedron, 
 ! multiplied by its volume, and the sum is divided by the hexaedron's volume
 nodes_tetraedron(1,:)=(/1,2,3,6/)
 nodes_tetraedron(2,:)=(/1,5,8,6/)
 nodes_tetraedron(3,:)=(/1,3,8,6/)
 nodes_tetraedron(4,:)=(/7,3,8,6/)
 nodes_tetraedron(5,:)=(/1,3,8,4/)
 vol_int=0.0
 do i=1,5 ! Loop over the 5 tetraedra

    ! Definition of the coordinates of each tetraedron
    do k=1,4
       knode(k)=cell(c)%knode(nodes_tetraedron(i,k))
       do j=1,3
          x(j,k)=node(knode(k))%x(j)
       enddo
    enddo

    ! Caculation of the volume of this tetraedron
    vol_tetraedron=abs((x(1,2)-x(1,1))*((x(2,3)-x(2,1))*(x(3,4)-x(3,1))-(x(2,4)-x(2,1))*(x(3,3)-x(3,1))) &
                   & -(x(2,2)-x(2,1))*((x(1,3)-x(1,1))*(x(3,4)-x(3,1))-(x(1,4)-x(1,1))*(x(3,3)-x(3,1))) &
                   & +(x(3,2)-x(3,1))*((x(1,3)-x(1,1))*(x(2,4)-x(2,1))-(x(1,4)-x(1,1))*(x(2,3)-x(2,1))))/6.0

    ! Integrate the polynomial term in this tetraedron
    call integrate_order3_volume_tetraedron_coordinates(x,alpha,beta,gamma,aux_vol_int)

    ! Multiply the integral and the volume of this tetraedron and add it to "vol_int"
    vol_int=vol_int + aux_vol_int*vol_tetraedron
 enddo


 ! Divide the sum by the hexaedron volume
 vol_int=vol_int/cell(c)%vol


end subroutine integrate_order3_volume_hexaedron















! This subroutine calculates the analytical integral of a third order polynomial term for a tetraedron cell
! This term is expressed as: alpha*beta*gamma;   The possible values for alpha,beta and gamma are: 1=x, 2=y, 3=z
subroutine integrate_order3_volume_tetraedron_coordinates(x,alpha,beta,gamma,vol_int)

use general_module


 ! Global variables
 integer, intent(in) :: alpha, beta, gamma
 double precision, intent(in), dimension(3,4) :: x ! x(1,1:4): X coordinates of the 4 tetrahedron's nodes
                                                   ! x(2,1:4): Y coordinates of the 4 tetrahedron's nodes
                                                   ! x(3,1:4): Z coordinates of the 4 tetrahedron's nodes
 double precision, intent(out) :: vol_int

 ! Local variables
 double precision a1, a2, a3, &
 & integral_tau3, integral_tau2_x_tau1, integral_tau_x_tau_x_tau
 integer i,j,k



 ! Calculation of: "integral_tau3", "integral_tau2_x_tau1",  "integral_tau_x_tau_x_tau"
 integral_tau3= 1.0/20.0
 integral_tau2_x_tau1= 1.0/60.0
 integral_tau_x_tau_x_tau= 1.0/120.0


 ! Calculation of "a1"
 a1=0.0
 do i=1,4
    a1=a1 + x(alpha,i)*x(beta,i)*x(gamma,i)
 enddo


 ! Calculation of "a2"
 a2=0.0
 do i=1,4
    do j=1,4
       if (j/=i) then
          a2=a2 + x(alpha,i)*x(beta,i)*x(gamma,j) + x(alpha,i)*x(gamma,i)*x(beta,j) + x(beta,i)*x(gamma,i)*x(alpha,j)
       endif
    enddo
 enddo


 ! Calculation of "a3"
 a3=0.0
 do i=1,4
    do j=1,4
       if (j/=i) then
          do k=1,4
             if ((k/=j).and.(k/=i)) then
                a3=a3 + x(alpha,i)*x(beta,j)*x(gamma,k)
             endif
          enddo
       endif
    enddo
 enddo  


 ! Calculation of the integral
 vol_int= a1*integral_tau3 + a2*integral_tau2_x_tau1 + a3*integral_tau_x_tau_x_tau

end subroutine integrate_order3_volume_tetraedron_coordinates















! This subroutine calculates the analytical integral of a third order polynomial term for a cuadrangle face
! This term is expressed as: alpha*beta*gamma;   The possible values for alpha,beta and gamma are: 1=x, 2=y, 3=z
subroutine integrate_order3_area_cuadrangle(facej,alpha,beta,gamma,area_int)

use general_module


 ! Global variables
 integer, intent(in) :: facej,alpha,beta,gamma
 double precision, intent(out) :: area_int

 ! Local variables
 double precision, dimension(3,3) :: x0 ! x(1,1:3): X coordinates of the 3 triangle's nodes
                                                    ! x(2,1:3): Y coordinates of the 3 triangle's nodes
                                                    ! x(3,1:3): Z coordinates of the 3 triangle's nodes
 double precision area_int1, area_int2, area1, area2, normtriangle(3), v1(3), v2(3)
 integer i,j,k


 ! Split the cuadrangle into 2 triangles: First one is composed of nodes 1, 2, 4 and the other is composed of
 ! nodes 2, 3, 4. 

 ! The integral is calculated by adding the integral of each triangle, 
 ! multiplied by its area, and the sum is divided by the cuadrangles's area which is the sum of the 2 triangle areas


 ! Calculation of the first triangle
 ! Definition of triangle coordinates
 do i=1,3
    x0(i,1)=node(face(facej)%knode(1))%x(i)
    x0(i,2)=node(face(facej)%knode(2))%x(i)
    x0(i,3)=node(face(facej)%knode(4))%x(i)
 enddo
 ! Calculation of the area
 do i=1,3
    v1(i)=x0(i,2)-x0(i,1)
    v2(i)=x0(i,3)-x0(i,1)
 enddo
 normtriangle = cross_product(v1,v2)
 area1 = vector_magnitude(normtriangle)/2.0
 ! Integrate the polynomial term in this triangle
 call integrate_order3_area_triangle_coordinates(facej,x0,alpha,beta,gamma,area_int1)


 ! Calculation of the second triangle
 ! Definition of triangle coordinates
 do i=1,3
    x0(i,1)=node(face(facej)%knode(2))%x(i)
    x0(i,2)=node(face(facej)%knode(3))%x(i)
    x0(i,3)=node(face(facej)%knode(4))%x(i)
 enddo
 ! Calculation of the area
 do i=1,3
    v1(i)=x0(i,2)-x0(i,1)
    v2(i)=x0(i,3)-x0(i,1)
 enddo
 normtriangle = cross_product(v1,v2)
 area2 = vector_magnitude(normtriangle)/2.0
 ! Integrate the polynomial term in this triangle
 call integrate_order3_area_triangle_coordinates(facej,x0,alpha,beta,gamma,area_int2)


 ! Calculation on the integral in the cuadrangle
 area_int=(area_int1*area1 + area_int2*area2)/(area1+area2)


end subroutine integrate_order3_area_cuadrangle















! This subroutine calculates the analytical integral of a third order polynomial term for a triangle face
! This term is expressed as: alpha*beta*gamma;   The possible values for alpha,beta and gamma are: 1=x, 2=y, 3=z
subroutine integrate_order3_area_triangle_coordinates(facej,x0,alpha,beta,gamma,area_int)

use general_module


 ! Global variables
 integer, intent(in) :: facej,alpha,beta,gamma
 double precision, intent(in), dimension(3,3) :: x0 ! x(1,1:3): X coordinates of the 3 triangle's nodes
                                                    ! x(2,1:3): Y coordinates of the 3 triangle's nodes
                                                    ! x(3,1:3): Z coordinates of the 3 triangle's nodes
 double precision, intent(out) :: area_int

 ! Local variables
 double precision x(3), y(3), z
 double precision a1, a2, a3, a4, &
 & a5, a6, a7, a8, a9, a10, b1, b2, b3, c1, c2, c3, &
 & d1, d2, d3, e1, e2, e3, f1, f2, g1, g2, h1, h2, m ,n , &
 & p1, p2, p3, p4, p5, p6, p7
 double precision integral_tau3, integral_tau2_x_tau1, integral_tau_x_tau_x_tau, &
 & integral_tau2, integral_tau_x_tau, integral_tau
 double precision cos_x0_x, cos_x0_y, cos_y0_x, cos_y0_y, &
 & cos_z0_x, cos_z0_y, cos_x0_z, cos_y0_z, cos_z0_z, rot(3,3)
 double precision, dimension(3) :: ux, uy, uz
 integer i,j,k



 ! Definition of the global axes vectors
 ux(:)=(/1.0 , 0.0 , 0.0/)
 uy(:)=(/0.0 , 1.0 , 0.0/)
 uz(:)=(/0.0 , 0.0 , 1.0/)


 ! Calculation of angles between local and global axes
 cos_x0_x=dot_product(ux,face(facej)%norm(:,2))
 cos_x0_y=dot_product(ux,face(facej)%norm(:,3))
 cos_x0_z=dot_product(ux,face(facej)%norm(:,1))
 cos_y0_x=dot_product(uy,face(facej)%norm(:,2))
 cos_y0_y=dot_product(uy,face(facej)%norm(:,3))
 cos_y0_z=dot_product(uy,face(facej)%norm(:,1))
 cos_z0_x=dot_product(uz,face(facej)%norm(:,2))
 cos_z0_y=dot_product(uz,face(facej)%norm(:,3))
 cos_z0_z=dot_product(uz,face(facej)%norm(:,1))


 ! Definition of the rotation matrix
 rot(1,1:3)=(/cos_x0_x, cos_x0_y, cos_x0_z/);
 rot(2,1:3)=(/cos_y0_x, cos_y0_y, cos_y0_z/);
 rot(3,1:3)=(/cos_z0_x, cos_z0_y, cos_z0_z/);


 ! Calculation of: "integral_tau3", "integral_tau2_x_tau1",  "integral_tau_x_tau_x_tau", 
 ! "integral_tau2", "integral_tau_x_tau", "integral_tau"
 integral_tau3= 1.0/10.0
 integral_tau2_x_tau1= 1.0/30.0
 integral_tau_x_tau_x_tau= 1.0/60.0
 integral_tau2=1.0/6.0
 integral_tau_x_tau=1.0/12.0
 integral_tau=1.0/3.0


 ! Calculation of "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10" 
 a1=rot(alpha,1)*rot(beta,1)*rot(gamma,1)
 a2=rot(alpha,2)*rot(beta,2)*rot(gamma,2)
 a3=rot(alpha,3)*rot(beta,3)*rot(gamma,3)
 a4=rot(alpha,1)*rot(beta,1)*rot(gamma,2) + rot(alpha,1)*rot(gamma,1)*rot(beta,2) + rot(beta,1)*rot(gamma,1)*rot(alpha,2)
 a5=rot(alpha,1)*rot(beta,1)*rot(gamma,3) + rot(alpha,1)*rot(gamma,1)*rot(beta,3) + rot(beta,1)*rot(gamma,1)*rot(alpha,3)
 a6=rot(alpha,2)*rot(beta,2)*rot(gamma,1) + rot(alpha,2)*rot(gamma,2)*rot(beta,1) + rot(beta,2)*rot(gamma,2)*rot(alpha,1)
 a7=rot(alpha,2)*rot(beta,2)*rot(gamma,3) + rot(alpha,2)*rot(gamma,2)*rot(beta,3) + rot(beta,2)*rot(gamma,2)*rot(alpha,3)
 a8=rot(alpha,3)*rot(beta,3)*rot(gamma,1) + rot(alpha,3)*rot(gamma,3)*rot(beta,1) + rot(beta,3)*rot(gamma,3)*rot(alpha,1)
 a9=rot(alpha,3)*rot(beta,3)*rot(gamma,2) + rot(alpha,3)*rot(gamma,3)*rot(beta,2) + rot(beta,3)*rot(gamma,3)*rot(alpha,2)
 a10=rot(alpha,1)*(rot(beta,2)*rot(gamma,3) + rot(gamma,2)*rot(beta,3)) + &
   & rot(beta,1)*(rot(alpha,2)*rot(gamma,3) + rot(gamma,2)*rot(alpha,3)) + &
   & rot(gamma,1)*(rot(alpha,2)*rot(beta,3) + rot(beta,2)*rot(alpha,3))

 
 ! Calculation of local coordinates of triangle coordinates
 do i=1,3
    x(i)=cos_x0_x*x0(1,i) + cos_y0_x*x0(2,i) + cos_z0_x*x0(3,i)
    y(i)=cos_x0_y*x0(1,i) + cos_y0_y*x0(2,i) + cos_z0_y*x0(3,i)
 enddo
 z=cos_x0_z*x0(1,1) + cos_y0_z*x0(2,1) + cos_z0_z*x0(3,1)


 ! Calculation of "b1", "c1", "d1", "e1", "f1", "g1", "h1", "m", "n"
 b1=0.0
 c1=0.0
 d1=0.0
 e1=0.0
 f1=0.0
 g1=0.0
 h1=0.0
 m=0.0
 n=0.0
 do i=1,3
    b1=b1 + x(i)**3
    c1=c1 + y(i)**3
    d1=d1 + x(i)**2*y(i)
    e1=e1 + x(i)*y(i)**2
    f1=f1 + x(i)**2
    g1=g1 + y(i)**2
    h1=h1 + x(i)*y(i)
    m=m + x(i)
    n=n + y(i)
 enddo


 ! Calculation of "b2", "c2", "d2", "e2", "f2", "g2", "h2"
 b2=0.0
 c2=0.0
 d2=0.0
 e2=0.0
 f2=0.0
 g2=0.0
 h2=0.0
 do i=1,3
    do j=1,3
       if (j/=i) then
          b2=b2 + 3*x(i)**2*x(j)
          c2=c2 + 3*y(i)**2*y(j)
          d2=d2 + x(i)**2*y(j) + 2*x(i)*y(i)*x(j)
          e2=e2 + y(i)**2*x(j) + 2*x(i)*y(i)*y(j)
          f2=f2 + x(i)*x(j)
          g2=g2 + y(i)*y(j)
          h2=h2 + x(i)*y(j)
       endif
    enddo
 enddo


 ! Calculation of "b3", "c3", "d3", "e3"
 b3=0.0
 c3=0.0
 d3=0.0
 e3=0.0
 do i=1,3
    do j=1,3
       if (j/=i) then
          do k=1,3
             if ((k/=j).and.(k/=i)) then
                b3=b3 + x(i)*x(j)*x(k)
                c3=c3 + y(i)*y(j)*y(k)
                d3=d3 + x(i)*x(j)*y(k)
                e3=e3 + y(i)*y(j)*x(k)
             endif
          enddo
       endif
    enddo
 enddo  


 
 ! Calculaltion of "p1", "p2", "p3", "p4", "p5", "p6", "p7"
 p1= a1*b1 + a2*c1 + a4*d1 + a6*e1
 p2= a1*b2 + a2*c2 + a4*d2 + a6*e2
 p3= a1*b3 + a2*c3 + a4*d3 + a6*e3 
 p4= a5*z*f1 + a7*z*g1 + a10*z*h1
 p5= a5*z*f2 + a7*z*g2 + a10*z*h2
 p6= a8*z**2*m + a9*z**2*n
 p7= a3*z**3


 ! Calculation of the integral
 area_int= p1*integral_tau3 + p2*integral_tau2_x_tau1 + p3*integral_tau_x_tau_x_tau + &
         & p4*integral_tau2 + p5*integral_tau_x_tau + p6*integral_tau + p7

end subroutine integrate_order3_area_triangle_coordinates















! This subroutine calculates the analytical integral of a fourth order polynomial term for a tetraedron cell
! This term is expressed as: alpha*beta*gamma*delta;   The possible values for alpha,beta,gamma and delta are: 1=x, 2=y, 3=z
subroutine integrate_order4_volume_tetraedron(c,alpha,beta,gamma,delta,vol_int)

use general_module


 ! Global variables
 integer, intent(in) :: c, alpha, beta, gamma, delta
 double precision, intent(out) :: vol_int

 ! Local variables
 double precision a1, a2, a3, a4, a5, &
 & integral_tau4, integral_tau3_x_tau1, integral_tau2_x_tau2, &
 & integral_tau2_x_tau_x_tau ,integral_tau_x_tau_x_tau_x_tau
 integer i,j,k,l
 double precision, dimension(3,4) :: x ! x(1,1:4): X coordinates of the 4 tetrahedron's nodes
                                       ! x(2,1:4): Y coordinates of the 4 tetrahedron's nodes
                                       ! x(3,1:4): Z coordinates of the 4 tetrahedron's nodes




 ! Definition of tetrahedron coordinates
 do i=1,3
    do j=1,4
       x(i,j)=node(cell(c)%knode(j))%x(i)
    enddo
 enddo



 ! Calculation of integrals in tetrehedral coordinates
 integral_tau4= 1.0/35.0
 integral_tau3_x_tau1= 1.0/140.0
 integral_tau2_x_tau2= 1.0/210.0
 integral_tau2_x_tau_x_tau= 1.0/420.0
 integral_tau_x_tau_x_tau_x_tau= 1.0/840.0



 ! Calculation of "a1"
 a1=0.0
 do i=1,4
    a1=a1 + x(alpha,i)*x(beta,i)*x(gamma,i)*x(delta,i)
 enddo


 ! Calculation of "a2"
 a2=0.0
 do i=1,4
    do j=1,4
       if (j/=i) then
          a2=a2 + x(alpha,i)*x(beta,i)*x(gamma,i)*x(delta,j) + x(alpha,i)*x(beta,i)*x(gamma,j)*x(delta,i) &
            &   + x(alpha,i)*x(beta,j)*x(gamma,i)*x(delta,i) + x(alpha,j)*x(beta,i)*x(gamma,i)*x(delta,i)
       endif
    enddo
 enddo


 ! Calculation of "a3"
 a3=0.0
 do i=1,4
    do j=1,4
       if (j>i) then
          a3=a3 + x(alpha,i)*x(beta,i)*x(gamma,j)*x(delta,j) + x(alpha,i)*x(beta,j)*x(gamma,i)*x(delta,j) &
            &   + x(alpha,i)*x(beta,j)*x(gamma,j)*x(delta,i) + x(alpha,j)*x(beta,i)*x(gamma,i)*x(delta,j) &
            &   + x(alpha,j)*x(beta,i)*x(gamma,j)*x(delta,i) + x(alpha,j)*x(beta,j)*x(gamma,i)*x(delta,i)
       endif
    enddo
 enddo


 ! Calculation of "a4"
 a4=0.0
 do i=1,4
    do j=1,4
       if (j/=i) then
          do k=1,4
             if ((k>j).and.(k/=i)) then
                a4=a4 + x(alpha,i)*x(beta,i)*x(gamma,j)*x(delta,k) + x(alpha,i)*x(beta,i)*x(gamma,k)*x(delta,j) &
                  &   + x(alpha,i)*x(beta,j)*x(gamma,i)*x(delta,k) + x(alpha,i)*x(beta,j)*x(gamma,k)*x(delta,i) &
                  &   + x(alpha,i)*x(beta,k)*x(gamma,i)*x(delta,j) + x(alpha,i)*x(beta,k)*x(gamma,j)*x(delta,i) &
                  &   + x(alpha,j)*x(beta,i)*x(gamma,i)*x(delta,k) + x(alpha,j)*x(beta,i)*x(gamma,k)*x(delta,i) &
                  &   + x(alpha,j)*x(beta,k)*x(gamma,i)*x(delta,i) + x(alpha,k)*x(beta,i)*x(gamma,i)*x(delta,j) &
                  &   + x(alpha,k)*x(beta,i)*x(gamma,j)*x(delta,i) + x(alpha,k)*x(beta,j)*x(gamma,i)*x(delta,i)
             endif
          enddo
       endif
    enddo
 enddo


 ! Calculation of "a5"
 a5=0.0
 do i=1,4
    do j=1,4
       if (j/=i) then
          do k=1,4
             if ((k/=i).and.(k/=j)) then
                do l=1,4
                   if ((l/=k).and.(l/=j).and.(l/=i)) then
                      a5=a5 + x(alpha,i)*x(beta,j)*x(gamma,k)*x(delta,l)
                   endif
                enddo
             endif
          enddo
       endif
    enddo
 enddo  


 ! Calculation of the integral
 vol_int= a1*integral_tau4 + a2*integral_tau3_x_tau1 + a3*integral_tau2_x_tau2 + a4*integral_tau2_x_tau_x_tau + &
        & a5*integral_tau_x_tau_x_tau_x_tau

end subroutine integrate_order4_volume_tetraedron














! This subroutine calculates the analytical integral of a fourth order polynomial term for a tetraedron cell
! This term is expressed as: alpha*beta*gamma*delta;   The possible values for alpha,beta,gamma and delta are: 1=x, 2=y, 3=z
subroutine integrate_order4_volume_tetraedron_coordinates(x,alpha,beta,gamma,delta,vol_int)

use general_module


 ! Global variables
 integer, intent(in) :: alpha, beta, gamma, delta
 double precision, intent(in), dimension(3,4) :: x ! x(1,1:4): X coordinates of the 4 tetrahedron's nodes
                                                   ! x(2,1:4): Y coordinates of the 4 tetrahedron's nodes
                                                   ! x(3,1:4): Z coordinates of the 4 tetrahedron's nodes
 double precision, intent(out) :: vol_int

 ! Local variables
 double precision a1, a2, a3, a4, a5, &
 & integral_tau4, integral_tau3_x_tau1, integral_tau2_x_tau2, &
 & integral_tau2_x_tau_x_tau ,integral_tau_x_tau_x_tau_x_tau
 integer i,j,k,l



 ! Calculation of integrals in tetrehedral coordinates
 integral_tau4= 1.0/35.0
 integral_tau3_x_tau1= 1.0/140.0
 integral_tau2_x_tau2= 1.0/210.0
 integral_tau2_x_tau_x_tau= 1.0/420.0
 integral_tau_x_tau_x_tau_x_tau= 1.0/840.0



 ! Calculation of "a1"
 a1=0.0
 do i=1,4
    a1=a1 + x(alpha,i)*x(beta,i)*x(gamma,i)*x(delta,i)
 enddo


 ! Calculation of "a2"
 a2=0.0
 do i=1,4
    do j=1,4
       if (j/=i) then
          a2=a2 + x(alpha,i)*x(beta,i)*x(gamma,i)*x(delta,j) + x(alpha,i)*x(beta,i)*x(gamma,j)*x(delta,i) &
            &   + x(alpha,i)*x(beta,j)*x(gamma,i)*x(delta,i) + x(alpha,j)*x(beta,i)*x(gamma,i)*x(delta,i)
       endif
    enddo
 enddo


 ! Calculation of "a3"
 a3=0.0
 do i=1,4
    do j=1,4
       if (j>i) then
          a3=a3 + x(alpha,i)*x(beta,i)*x(gamma,j)*x(delta,j) + x(alpha,i)*x(beta,j)*x(gamma,i)*x(delta,j) &
            &   + x(alpha,i)*x(beta,j)*x(gamma,j)*x(delta,i) + x(alpha,j)*x(beta,i)*x(gamma,i)*x(delta,j) &
            &   + x(alpha,j)*x(beta,i)*x(gamma,j)*x(delta,i) + x(alpha,j)*x(beta,j)*x(gamma,i)*x(delta,i)
       endif
    enddo
 enddo


 ! Calculation of "a4"
 a4=0.0
 do i=1,4
    do j=1,4
       if (j/=i) then
          do k=1,4
             if ((k>j).and.(k/=i)) then
                a4=a4 + x(alpha,i)*x(beta,i)*x(gamma,j)*x(delta,k) + x(alpha,i)*x(beta,i)*x(gamma,k)*x(delta,j) &
                  &   + x(alpha,i)*x(beta,j)*x(gamma,i)*x(delta,k) + x(alpha,i)*x(beta,j)*x(gamma,k)*x(delta,i) &
                  &   + x(alpha,i)*x(beta,k)*x(gamma,i)*x(delta,j) + x(alpha,i)*x(beta,k)*x(gamma,j)*x(delta,i) &
                  &   + x(alpha,j)*x(beta,i)*x(gamma,i)*x(delta,k) + x(alpha,j)*x(beta,i)*x(gamma,k)*x(delta,i) &
                  &   + x(alpha,j)*x(beta,k)*x(gamma,i)*x(delta,i) + x(alpha,k)*x(beta,i)*x(gamma,i)*x(delta,j) &
                  &   + x(alpha,k)*x(beta,i)*x(gamma,j)*x(delta,i) + x(alpha,k)*x(beta,j)*x(gamma,i)*x(delta,i)
             endif
          enddo
       endif
    enddo
 enddo


 ! Calculation of "a5"
 a5=0.0
 do i=1,4
    do j=1,4
       if (j/=i) then
          do k=1,4
             if ((k/=i).and.(k/=j)) then
                do l=1,4
                   if ((l/=k).and.(l/=j).and.(l/=i)) then
                      a5=a5 + x(alpha,i)*x(beta,j)*x(gamma,k)*x(delta,l)
                   endif
                enddo
             endif
          enddo
       endif
    enddo
 enddo  


 ! Calculation of the integral
 vol_int= a1*integral_tau4 + a2*integral_tau3_x_tau1 + a3*integral_tau2_x_tau2 + a4*integral_tau2_x_tau_x_tau + &
        & a5*integral_tau_x_tau_x_tau_x_tau

end subroutine integrate_order4_volume_tetraedron_coordinates















! This subroutine calculates the analytical integral of a fourth order polynomial term for a triangle face
! This term is expressed as: alpha*beta*gamma*delta;   The possible values for alpha,beta,gamma and delta are: 1=x, 2=y, 3=z
subroutine integrate_order4_area_triangle_coordinates(facej,x0,alpha,beta,gamma,delta,area_int)

use general_module


 ! Global variables
 integer, intent(in) :: facej,alpha,beta,gamma,delta
 double precision, intent(in), dimension(3,3) :: x0 ! x(1,1:3): X coordinates of the 3 triangle's nodes
                                                    ! x(2,1:3): Y coordinates of the 3 triangle's nodes
                                                    ! x(3,1:3): Z coordinates of the 3 triangle's nodes
 double precision, intent(out) :: area_int

 ! Local variables
 double precision x(3), y(3), z
 double precision b(15), a(10,15), c(11)
 double precision integral_tau3, integral_tau2_x_tau1, integral_tau_x_tau_x_tau, &
 & integral_tau2, integral_tau_x_tau, integral_tau, &
 & integral_tau4, integral_tau3_x_tau1, integral_tau2_x_tau2, &
 & integral_tau2_x_tau_x_tau
 double precision cos_x0_x, cos_x0_y, cos_y0_x, cos_y0_y, &
 & cos_z0_x, cos_z0_y, cos_x0_z, cos_y0_z, cos_z0_z, rot(3,3)
 double precision, dimension(3) :: ux, uy, uz
 integer i,j,k, aux(4)



 ! Definition of the global axes vectors
 ux(:)=(/1.0 , 0.0 , 0.0/)
 uy(:)=(/0.0 , 1.0 , 0.0/)
 uz(:)=(/0.0 , 0.0 , 1.0/)


 ! Calculation of angles between local and global axes
 cos_x0_x=dot_product(ux,face(facej)%norm(:,2))
 cos_x0_y=dot_product(ux,face(facej)%norm(:,3))
 cos_x0_z=dot_product(ux,face(facej)%norm(:,1))
 cos_y0_x=dot_product(uy,face(facej)%norm(:,2))
 cos_y0_y=dot_product(uy,face(facej)%norm(:,3))
 cos_y0_z=dot_product(uy,face(facej)%norm(:,1))
 cos_z0_x=dot_product(uz,face(facej)%norm(:,2))
 cos_z0_y=dot_product(uz,face(facej)%norm(:,3))
 cos_z0_z=dot_product(uz,face(facej)%norm(:,1))


 ! Definition of the rotation matrix
 rot(1,1:3)=(/cos_x0_x, cos_x0_y, cos_x0_z/);
 rot(2,1:3)=(/cos_y0_x, cos_y0_y, cos_y0_z/);
 rot(3,1:3)=(/cos_z0_x, cos_z0_y, cos_z0_z/);


 ! Calculation of integrals in triangular coordinates
 integral_tau4= 1.0/15.0
 integral_tau3_x_tau1= 1.0/60.0
 integral_tau2_x_tau2= 1.0/90.0
 integral_tau2_x_tau_x_tau= 1.0/180.0
 integral_tau3= 1.0/10.0
 integral_tau2_x_tau1= 1.0/30.0
 integral_tau_x_tau_x_tau= 1.0/60.0
 integral_tau2=1.0/6.0
 integral_tau_x_tau=1.0/12.0
 integral_tau=1.0/3.0

 
 ! Calculation of local coordinates of triangle coordinates
 do i=1,3
    x(i)=cos_x0_x*x0(1,i) + cos_y0_x*x0(2,i) + cos_z0_x*x0(3,i)
    y(i)=cos_x0_y*x0(1,i) + cos_y0_y*x0(2,i) + cos_z0_y*x0(3,i)
 enddo
 z=cos_x0_z*x0(1,1) + cos_y0_z*x0(2,1) + cos_z0_z*x0(3,1)


 ! Calculation of "b": 
 ! x0(alpha)*x0(beta)*x0(gamma)*x0(delta)= b(1)*x**4 + b(2)*y**4 + b(3)*x**3*y + b(4)*y**3*x +
 ! + b(5)*x**2*y**2 + b(6)*x**3 + b(7)*y**3 + b(8)*x**2*y + b(9)*y**2*x + b(10)*x**2 + b(11)*y**2 + 
 ! + b(12)*x*y + b(13)*x + b(14)*y + b15
 b(1)=rot(alpha,1)*rot(beta,1)*rot(gamma,1)*rot(delta,1)
 b(2)=rot(alpha,2)*rot(beta,2)*rot(gamma,2)*rot(delta,2)
 b(15)=rot(alpha,3)*rot(beta,3)*rot(gamma,3)*rot(delta,3)*z**4
 b(3)=rot(alpha,1)*rot(beta,1)*rot(gamma,1)*rot(delta,2)+rot(alpha,1)*rot(beta,1)*rot(gamma,2)*rot(delta,1) &
  &  +rot(alpha,1)*rot(beta,2)*rot(gamma,1)*rot(delta,1)+rot(alpha,2)*rot(beta,1)*rot(gamma,1)*rot(delta,1) 
 b(4)=rot(alpha,2)*rot(beta,2)*rot(gamma,2)*rot(delta,1)+rot(alpha,2)*rot(beta,2)*rot(gamma,1)*rot(delta,2) &
  &  +rot(alpha,2)*rot(beta,1)*rot(gamma,2)*rot(delta,2)+rot(alpha,1)*rot(beta,2)*rot(gamma,2)*rot(delta,2) 
 b(6)=(rot(alpha,1)*rot(beta,1)*rot(gamma,1)*rot(delta,3)+rot(alpha,1)*rot(beta,1)*rot(gamma,3)*rot(delta,1) &
  &   +rot(alpha,1)*rot(beta,3)*rot(gamma,1)*rot(delta,1)+rot(alpha,3)*rot(beta,1)*rot(gamma,1)*rot(delta,1))*z 
 b(7)=(rot(alpha,2)*rot(beta,2)*rot(gamma,2)*rot(delta,3)+rot(alpha,2)*rot(beta,2)*rot(gamma,3)*rot(delta,2) &
  &   +rot(alpha,2)*rot(beta,3)*rot(gamma,2)*rot(delta,2)+rot(alpha,3)*rot(beta,2)*rot(gamma,2)*rot(delta,2))*z 
 b(13)=(rot(alpha,3)*rot(beta,3)*rot(gamma,3)*rot(delta,1)+rot(alpha,3)*rot(beta,3)*rot(gamma,1)*rot(delta,3) &
  &    +rot(alpha,3)*rot(beta,1)*rot(gamma,3)*rot(delta,3)+rot(alpha,1)*rot(beta,3)*rot(gamma,3)*rot(delta,3))*z**3 
 b(14)=(rot(alpha,3)*rot(beta,3)*rot(gamma,3)*rot(delta,2)+rot(alpha,3)*rot(beta,3)*rot(gamma,2)*rot(delta,3) &
  &    +rot(alpha,3)*rot(beta,2)*rot(gamma,3)*rot(delta,3)+rot(alpha,2)*rot(beta,3)*rot(gamma,3)*rot(delta,3))*z**3
 b(5)=rot(alpha,1)*rot(beta,1)*rot(gamma,2)*rot(delta,2)+rot(alpha,1)*rot(beta,2)*rot(gamma,1)*rot(delta,2) &
  &  +rot(alpha,1)*rot(beta,2)*rot(gamma,2)*rot(delta,1)+rot(alpha,2)*rot(beta,1)*rot(gamma,1)*rot(delta,2) & 
  &  +rot(alpha,2)*rot(beta,1)*rot(gamma,2)*rot(delta,1)+rot(alpha,2)*rot(beta,2)*rot(gamma,1)*rot(delta,1)
 b(10)=(rot(alpha,1)*rot(beta,1)*rot(gamma,3)*rot(delta,3)+rot(alpha,1)*rot(beta,3)*rot(gamma,1)*rot(delta,3) &
  &    +rot(alpha,1)*rot(beta,3)*rot(gamma,3)*rot(delta,1)+rot(alpha,3)*rot(beta,1)*rot(gamma,1)*rot(delta,3) & 
  &    +rot(alpha,3)*rot(beta,1)*rot(gamma,3)*rot(delta,1)+rot(alpha,3)*rot(beta,3)*rot(gamma,1)*rot(delta,1))*z**2
 b(11)=(rot(alpha,2)*rot(beta,2)*rot(gamma,3)*rot(delta,3)+rot(alpha,2)*rot(beta,3)*rot(gamma,2)*rot(delta,3) &
  &    +rot(alpha,2)*rot(beta,3)*rot(gamma,3)*rot(delta,2)+rot(alpha,3)*rot(beta,2)*rot(gamma,2)*rot(delta,3) & 
  &    +rot(alpha,3)*rot(beta,2)*rot(gamma,3)*rot(delta,2)+rot(alpha,3)*rot(beta,3)*rot(gamma,2)*rot(delta,2))*z**2
 b(8)=(rot(alpha,1)*rot(beta,1)*rot(gamma,2)*rot(delta,3)+rot(alpha,1)*rot(beta,1)*rot(gamma,3)*rot(delta,2) &
  &   +rot(alpha,1)*rot(beta,2)*rot(gamma,1)*rot(delta,3)+rot(alpha,1)*rot(beta,2)*rot(gamma,3)*rot(delta,1) & 
  &   +rot(alpha,1)*rot(beta,3)*rot(gamma,1)*rot(delta,2)+rot(alpha,1)*rot(beta,3)*rot(gamma,2)*rot(delta,1) & 
  &   +rot(alpha,2)*rot(beta,1)*rot(gamma,1)*rot(delta,3)+rot(alpha,2)*rot(beta,1)*rot(gamma,3)*rot(delta,1) & 
  &   +rot(alpha,2)*rot(beta,3)*rot(gamma,1)*rot(delta,1)+rot(alpha,3)*rot(beta,1)*rot(gamma,1)*rot(delta,2) & 
  &   +rot(alpha,3)*rot(beta,1)*rot(gamma,2)*rot(delta,1)+rot(alpha,3)*rot(beta,2)*rot(gamma,1)*rot(delta,1))*z
 b(9)=(rot(alpha,2)*rot(beta,2)*rot(gamma,1)*rot(delta,3)+rot(alpha,2)*rot(beta,2)*rot(gamma,3)*rot(delta,1) &
  &   +rot(alpha,2)*rot(beta,1)*rot(gamma,2)*rot(delta,3)+rot(alpha,2)*rot(beta,1)*rot(gamma,3)*rot(delta,2) & 
  &   +rot(alpha,2)*rot(beta,3)*rot(gamma,2)*rot(delta,1)+rot(alpha,2)*rot(beta,3)*rot(gamma,1)*rot(delta,2) & 
  &   +rot(alpha,1)*rot(beta,2)*rot(gamma,2)*rot(delta,3)+rot(alpha,1)*rot(beta,2)*rot(gamma,3)*rot(delta,2) & 
  &   +rot(alpha,1)*rot(beta,3)*rot(gamma,2)*rot(delta,2)+rot(alpha,3)*rot(beta,2)*rot(gamma,2)*rot(delta,1) & 
  &   +rot(alpha,3)*rot(beta,2)*rot(gamma,1)*rot(delta,2)+rot(alpha,3)*rot(beta,1)*rot(gamma,2)*rot(delta,2))*z
 b(12)=(rot(alpha,3)*rot(beta,3)*rot(gamma,2)*rot(delta,1)+rot(alpha,3)*rot(beta,3)*rot(gamma,1)*rot(delta,2) &
  &    +rot(alpha,3)*rot(beta,2)*rot(gamma,3)*rot(delta,1)+rot(alpha,3)*rot(beta,2)*rot(gamma,1)*rot(delta,3) & 
  &    +rot(alpha,3)*rot(beta,1)*rot(gamma,3)*rot(delta,2)+rot(alpha,3)*rot(beta,1)*rot(gamma,2)*rot(delta,3) & 
  &    +rot(alpha,2)*rot(beta,3)*rot(gamma,3)*rot(delta,1)+rot(alpha,2)*rot(beta,3)*rot(gamma,1)*rot(delta,3) & 
  &    +rot(alpha,2)*rot(beta,1)*rot(gamma,3)*rot(delta,3)+rot(alpha,1)*rot(beta,3)*rot(gamma,3)*rot(delta,2) & 
  &    +rot(alpha,1)*rot(beta,3)*rot(gamma,2)*rot(delta,3)+rot(alpha,1)*rot(beta,2)*rot(gamma,3)*rot(delta,3))*z**2




 ! Calculation of "a":
 ! x(alpha)*x(beta)*x(gamma)*x(delta)= a(1)*ti**4 + a(2)*ti**3*tj + a(3)*ti**2*tj**2 + a(4)*ti**2*tj*tk 
 ! x(alpha)*x(beta)*x(gamma)= a(5)*ti**3 + a(6)*ti**2*tj + a(7)*ti*tj*tk 
 ! x(alpha)*x(beta)= a(8)*ti**2 + a(9)*ti*tj 
 a=0.0
 do i=1,3
    a(1,1)=a(1,1)+x(i)**4
    a(1,2)=a(1,2)+y(i)**4
    a(1,3)=a(1,3)+x(i)**3*y(i)
    a(1,4)=a(1,4)+y(i)**3*x(i)
    a(1,5)=a(1,5)+x(i)**2*y(i)**2

    a(5,6)=a(5,6)+x(i)**3
    a(5,7)=a(5,7)+y(i)**3
    a(5,8)=a(5,8)+x(i)**2*y(i)
    a(5,9)=a(5,9)+y(i)**2*x(i)

    a(8,10)=a(8,10)+x(i)**2
    a(8,11)=a(8,11)+y(i)**2
    a(8,12)=a(8,12)+x(i)*y(i)

    a(10,13)=a(10,13)+x(i)
    a(10,14)=a(10,14)+y(i)
    do j=1,3
       if (j>i) then
          a(3,1)=a(3,1)+6*x(i)**2*x(j)**2
          a(3,2)=a(3,2)+6*y(i)**2*y(j)**2
          a(3,3)=a(3,3)+3*x(i)**2*x(j)*y(j)+3*x(j)**2*x(i)*y(i)
          a(3,4)=a(3,4)+3*y(i)**2*y(j)*x(j)+3*y(j)**2*y(i)*x(i)
          a(3,5)=a(3,5)+x(i)**2*y(j)**2+4*x(i)*x(j)*y(i)*y(j)+x(j)**2*y(i)**2
       endif
       if (j/=i) then
          a(2,1)=a(2,1)+4*x(i)**3*x(j)
          a(2,2)=a(2,2)+4*y(i)**3*y(j)
          a(2,3)=a(2,3)+x(i)**3*y(j)+3*x(i)**2*x(j)*y(i)
          a(2,4)=a(2,4)+y(i)**3*x(j)+3*y(i)**2*y(j)*x(i)
          a(2,5)=a(2,5)+2*x(i)**2*y(i)*y(j)+2*y(i)**2*x(i)*x(j)

          a(6,6)=a(6,6)+3*x(i)**2*x(j)
          a(6,7)=a(6,7)+3*y(i)**2*y(j)
          a(6,8)=a(6,8)+x(i)**2*y(j)+2*x(i)*x(j)*y(i)
          a(6,9)=a(6,9)+y(i)**2*x(j)+2*y(i)*y(j)*x(i)

          a(9,10)=a(9,10)+x(i)*x(j)
          a(9,11)=a(9,11)+y(i)*y(j)
          a(9,12)=a(9,12)+x(i)*y(j)
          do k=1,3
             if ((k/=i).and.(k/=j)) then
                if (k>j) then 
                   a(4,1)=a(4,1)+12*x(i)**2*x(j)*x(k)
                   a(4,2)=a(4,2)+12*y(i)**2*y(j)*y(k)
                   a(4,3)=a(4,3)+3*x(i)**2*x(j)*y(k)+3*x(i)**2*x(k)*y(j)+6*x(i)*x(j)*x(k)*y(i)
                   a(4,4)=a(4,4)+3*y(i)**2*y(j)*x(k)+3*y(i)**2*y(k)*x(j)+6*y(i)*y(j)*y(k)*x(i)
                   a(4,5)=a(4,5)+2*x(i)**2*y(j)*y(k)+4*x(i)*x(j)*y(i)*y(k)+4*x(i)*x(k)*y(i)*y(j)&
                         & +2*x(j)*x(k)*y(i)**2
                endif

                a(7,6)=a(7,6)+x(i)*x(j)*x(k)
                a(7,7)=a(7,7)+y(i)*y(j)*y(k)
                a(7,8)=a(7,8)+x(i)*x(j)*y(k)
                a(7,9)=a(7,9)+y(i)*y(j)*x(k)
             endif
          enddo
       endif
    enddo
 enddo



 ! Calculation of "c":
 ! x0(alpha)*x0(beta)*x0(gamma)*x0(delta)= c(1)*ti**4 + c(2)*ti**3*tj + c(3)*ti**2*tj**2 + c(4)*ti**2*tj*tk + 
 ! + c(5)*ti**3 + c(6)*ti**2*tj + c(7)*ti*tj*tk + c(8)*ti**2 + c(9)*ti*tj + c(10)*ti + c(11)
 c=0.0
 do j=1,4
    do i=1,5
       c(j)=c(j)+b(i)*a(j,i)
    enddo
 enddo
 do j=5,7
    do i=6,9
       c(j)=c(j)+b(i)*a(j,i)
    enddo
 enddo
 do j=8,9
    do i=10,12
       c(j)=c(j)+b(i)*a(j,i)
    enddo
 enddo
 c(10)=a(10,13)*b(13)+a(10,14)*b(14)
 c(11)=b(15)




 ! Calculation of the integral
 area_int= c(1)*integral_tau4 + c(2)*integral_tau3_x_tau1 + c(3)*integral_tau2_x_tau2 + c(4)*integral_tau2_x_tau_x_tau + &
         & c(5)*integral_tau3 + c(6)*integral_tau2_x_tau1 + c(7)*integral_tau_x_tau_x_tau + &
         & c(8)*integral_tau2 + c(9)*integral_tau_x_tau + c(10)*integral_tau + c(11)

end subroutine integrate_order4_area_triangle_coordinates
















! This subroutine calculates the analytical integral of a fourth order polynomial term for a hexaedron cell
! This term is expressed as: alpha*beta*gamma;   The possible values for alpha,beta,gamma and delta are: 1=x, 2=y, 3=z
subroutine integrate_order4_volume_hexaedron(c,alpha,beta,gamma,delta,vol_int)

use general_module


 ! Global variables
 integer, intent(in) :: c, alpha, beta, gamma, delta
 double precision, intent(out) :: vol_int

 ! Local variables
 double precision x(3,4), vol_tetraedron, aux_vol_int
 integer i,j,k, nodes_tetraedron(5,4), knode(4)


 ! The hexaedron is divided into 5 tetraedron and the integral is calculated by adding the integral of each tetraedron, 
 ! multiplied by its volume, and the sum is divided by the hexaedron's volume
 nodes_tetraedron(1,:)=(/1,2,3,6/)
 nodes_tetraedron(2,:)=(/1,5,8,6/)
 nodes_tetraedron(3,:)=(/1,3,8,6/)
 nodes_tetraedron(4,:)=(/7,3,8,6/)
 nodes_tetraedron(5,:)=(/1,3,8,4/)
 vol_int=0.0
 do i=1,5 ! Loop over the 5 tetraedra

    ! Definition of the coordinates of each tetraedron
    do k=1,4
       knode(k)=cell(c)%knode(nodes_tetraedron(i,k))
       do j=1,3
          x(j,k)=node(knode(k))%x(j)
       enddo
    enddo

    ! Caculation of the volume of this tetraedron
    vol_tetraedron=abs((x(1,2)-x(1,1))*((x(2,3)-x(2,1))*(x(3,4)-x(3,1))-(x(2,4)-x(2,1))*(x(3,3)-x(3,1))) &
                   & -(x(2,2)-x(2,1))*((x(1,3)-x(1,1))*(x(3,4)-x(3,1))-(x(1,4)-x(1,1))*(x(3,3)-x(3,1))) &
                   & +(x(3,2)-x(3,1))*((x(1,3)-x(1,1))*(x(2,4)-x(2,1))-(x(1,4)-x(1,1))*(x(2,3)-x(2,1))))/6.0

    ! Integrate the polynomial term in this tetraedron
    call integrate_order4_volume_tetraedron_coordinates(x,alpha,beta,gamma,delta,aux_vol_int)

    ! Multiply the integral and the volume of this tetraedron and add it to "vol_int"
    vol_int=vol_int + aux_vol_int*vol_tetraedron
 enddo


 ! Divide the sum by the hexaedron volume
 vol_int=vol_int/cell(c)%vol


end subroutine integrate_order4_volume_hexaedron















! This subroutine calculates the analytical integral of a fourth order polynomial term for a cuadrangle face
! This term is expressed as: alpha*beta*gamma;   The possible values for alpha,beta,gamma and delta are: 1=x, 2=y, 3=z
subroutine integrate_order4_area_cuadrangle(facej,alpha,beta,gamma,delta,area_int)

use general_module


 ! Global variables
 integer, intent(in) :: facej,alpha,beta,gamma, delta
 double precision, intent(out) :: area_int

 ! Local variables
 double precision, dimension(3,3) :: x0 ! x(1,1:3): X coordinates of the 3 triangle's nodes
                                                    ! x(2,1:3): Y coordinates of the 3 triangle's nodes
                                                    ! x(3,1:3): Z coordinates of the 3 triangle's nodes
 double precision area_int1, area_int2, area1, area2, normtriangle(3), v1(3), v2(3)
 integer i,j,k


 ! Split the cuadrangle into 2 triangles: First one is composed of nodes 1, 2, 4 and the other is composed of
 ! nodes 2, 3, 4. 

 ! The integral is calculated by adding the integral of each triangle, 
 ! multiplied by its area, and the sum is divided by the cuadrangles's area which is the sum of the 2 triangle areas


 ! Calculation of the first triangle
 ! Definition of triangle coordinates
 do i=1,3
    x0(i,1)=node(face(facej)%knode(1))%x(i)
    x0(i,2)=node(face(facej)%knode(2))%x(i)
    x0(i,3)=node(face(facej)%knode(4))%x(i)
 enddo
 ! Calculation of the area
 do i=1,3
    v1(i)=x0(i,2)-x0(i,1)
    v2(i)=x0(i,3)-x0(i,1)
 enddo
 normtriangle = cross_product(v1,v2)
 area1 = vector_magnitude(normtriangle)/2.0
 ! Integrate the polynomial term in this triangle
 call integrate_order4_area_triangle_coordinates(facej,x0,alpha,beta,gamma,delta,area_int1)


 ! Calculation of the second triangle
 ! Definition of triangle coordinates
 do i=1,3
    x0(i,1)=node(face(facej)%knode(2))%x(i)
    x0(i,2)=node(face(facej)%knode(3))%x(i)
    x0(i,3)=node(face(facej)%knode(4))%x(i)
 enddo
 ! Calculation of the area
 do i=1,3
    v1(i)=x0(i,2)-x0(i,1)
    v2(i)=x0(i,3)-x0(i,1)
 enddo
 normtriangle = cross_product(v1,v2)
 area2 = vector_magnitude(normtriangle)/2.0
 ! Integrate the polynomial term in this triangle
 call integrate_order4_area_triangle_coordinates(facej,x0,alpha,beta,gamma,delta,area_int2)


 ! Calculation on the integral in the cuadrangle
 area_int=(area_int1*area1 + area_int2*area2)/(area1+area2)


end subroutine integrate_order4_area_cuadrangle















! This subroutine calculates the volume averaged valaues of integral volume for a tetraedron cell of the following
! polynomial term expressed in tetrahedral coordinates: 
! p(t1,t2,t3,t4)= t1**n1 * t2**n2 * t3**n3 * t4**n4  
! ti are the tetrahedral coordinates and ni are the polynomial exponents
subroutine volume_integrate_tetraedral_coord(alpha,beta,gamma,delta,vol_int)

use general_module


 ! Global variables
 integer, intent(in) :: alpha,beta,gamma,delta
 double precision, intent(out) :: vol_int

 ! Local variables
 integer i, k, polynomial(4)
 double precision   p(5)


 ! Polynomial term
 polynomial(1)=alpha
 polynomial(2)=beta
 polynomial(3)=gamma
 polynomial(4)=delta

 ! Calculation of the factorial products to determine the volume integral 
 do k=1,4
    p(k)=1.0
    if (polynomial(k)/=0) then
       do i=1,polynomial(k)
          p(k)=p(k)*i
       enddo
    endif
 enddo
 k=5
 p(k)=1.0
 do i=1,sum(polynomial(:))+3
    p(k)=p(k)*i
 enddo    


 ! Calculate the volume integral
 vol_int=p(1)*p(2)*p(3)*p(4)*6.0/p(5)

 
end subroutine volume_integrate_tetraedral_coord















! This subroutine calculates the surface averaged valaues of surface integral for a triangular face of the following
! polynomial term expressed in triangular coordinates: 
! p(t1,t2,t3)= t1**n1 * t2**n2 * t3**n3  
! ti are the triangular coordinates and ni are the polynomial exponents
subroutine area_integrate_triangular_coord(alpha,beta,gamma,area_int)

use general_module


 ! Global variables
 integer, intent(in) :: alpha,beta,gamma
 double precision, intent(out) :: area_int

 ! Local variables
 integer i, k, polynomial(3)
 double precision   p(4)


 ! Polynomial term
 polynomial(1)=alpha
 polynomial(2)=beta
 polynomial(3)=gamma

 ! Calculation of the factorial products to determine the surface integral
 do k=1,3
    p(k)=1.0
    if (polynomial(k)/=0) then
       do i=1,polynomial(k)
          p(k)=p(k)*i
       enddo
    endif
 enddo
 k=4
 p(k)=1.0
 do i=1,sum(polynomial(:))+2
    p(k)=p(k)*i
 enddo    


 ! Calculate the surface integral
 area_int=p(1)*p(2)*p(3)*2.0/p(4)
 
end subroutine area_integrate_triangular_coord














! This subroutine calculates the analytical volume integral of a N order polynomial term for a tetraedron cell
! This term is expressed as: x**nx * y**ny * z**nz 
! nx, ny, nz are no negative integers 
subroutine integrate_orderN_volume_tetraedron(x,nx,ny,nz,vol_int)

use general_module


 ! Global variables
 integer, intent(in) :: nx,ny,nz
 double precision, intent(in), dimension(3,4) :: x ! x(1,1:4): X coordinates of the 4 tetrahedron's nodes
                                                   ! x(2,1:4): Y coordinates of the 4 tetrahedron's nodes
                                                   ! x(3,1:4): Z coordinates of the 4 tetrahedron's nodes
 double precision, intent(out) :: vol_int

 ! Local variables
 double precision prod, aux_vol_int
 integer n,m,j, alpha, beta, gamma, delta, cont
 integer, allocatable, dimension(:) :: i, xyz


 ! Determine the polynomial order (n)
 n=nx+ny+nz

 ! Calculation of the volume integral
 if (n>0) then
    vol_int=0.0
    allocate(i(n),xyz(n))
    i=1
    cont=0
    if (nx>0) then
       do j=1,nx
          cont=cont+1
          xyz(cont)=1
       enddo
    endif 
    if (ny>0) then
       do j=1,ny
          cont=cont+1
          xyz(cont)=2
       enddo
    endif 
    if (nz>0) then
       do j=1,nz
          cont=cont+1
          xyz(cont)=3
       enddo
    endif 

    do cont=1,4**n
       prod=1.0
       alpha=0
       beta=0
       gamma=0
       delta=0
       do m=1,n
          select case(i(m))
          case(1)
             alpha=alpha+1
          case(2)
             beta=beta+1
          case(3)
             gamma=gamma+1
          case(4)
             delta=delta+1
          end select
          prod=prod*x(xyz(m),i(m))
       enddo
       call volume_integrate_tetraedral_coord(alpha,beta,gamma,delta,aux_vol_int)
       vol_int=vol_int+prod*aux_vol_int
       i(n)=i(n)+1
       j=n
       do
          if((i(j)>4).and.(cont<4**n)) then
             i(j)=1
             i(j-1)=i(j-1)+1
             j=j-1
          else
             exit
          endif
       enddo
    enddo
    deallocate(i,xyz)
 else
    vol_int=1.0
 endif


end subroutine integrate_orderN_volume_tetraedron














! This subroutine calculates the analytical surface integral of a N order polynomial term for a triangular face
! This term is expressed as: x**nx * y**ny * z**nz 
! nx, ny, nz are no negative integers 
subroutine integrate_orderN_area_triangle(facej,x0,nx,ny,nz,area_int)

use general_module


 ! Global variables
 integer, intent(in) :: nx,ny,nz, facej
 double precision, intent(in), dimension(3,3) :: x0 ! x(1,1:3): X coordinates of the 3 triangle's nodes in global coordinates
                                                    ! x(2,1:3): Y coordinates of the 3 triangle's nodes in global coordinates
                                                    ! x(3,1:3): Z coordinates of the 3 triangle's nodes in global coordinates
 double precision, intent(out) :: area_int

 ! Local variables
 double precision prod, aux_area_int, z, rot(3,3), x(3), y(3), ux(3), uy(3), uz(3), &
 & cos_x0_x, cos_x0_y, cos_x0_z, cos_y0_x, cos_y0_y, cos_y0_z, &
 & cos_z0_x, cos_z0_y, cos_z0_z
 integer n,m,j, alpha, beta, gamma, cont
 integer, allocatable, dimension(:) :: i, x0y0z0





 ! Definition of the global axes vectors
 ux(:)=(/1.0 , 0.0 , 0.0/)
 uy(:)=(/0.0 , 1.0 , 0.0/)
 uz(:)=(/0.0 , 0.0 , 1.0/)

 ! Calculation of angles between local and global axes
 cos_x0_x=dot_product(ux,face(facej)%norm(:,2))
 cos_x0_y=dot_product(ux,face(facej)%norm(:,3))
 cos_x0_z=dot_product(ux,face(facej)%norm(:,1))
 cos_y0_x=dot_product(uy,face(facej)%norm(:,2))
 cos_y0_y=dot_product(uy,face(facej)%norm(:,3))
 cos_y0_z=dot_product(uy,face(facej)%norm(:,1))
 cos_z0_x=dot_product(uz,face(facej)%norm(:,2))
 cos_z0_y=dot_product(uz,face(facej)%norm(:,3))
 cos_z0_z=dot_product(uz,face(facej)%norm(:,1))

 ! Definition of the rotation matrix
 rot(1,1:3)=(/cos_x0_x, cos_x0_y, cos_x0_z/);
 rot(2,1:3)=(/cos_y0_x, cos_y0_y, cos_y0_z/);
 rot(3,1:3)=(/cos_z0_x, cos_z0_y, cos_z0_z/);

  ! Calculation of local coordinates of triangle coordinates
 do j=1,3
    x(j)=cos_x0_x*x0(1,j) + cos_y0_x*x0(2,j) + cos_z0_x*x0(3,j)
    y(j)=cos_x0_y*x0(1,j) + cos_y0_y*x0(2,j) + cos_z0_y*x0(3,j)
 enddo
 z=cos_x0_z*x0(1,1) + cos_y0_z*x0(2,1) + cos_z0_z*x0(3,1)




 ! Determine the polynomial order (n)
 n=nx+ny+nz

 ! Calculation of the surface integral
 if (n>0) then
    area_int=0.0
    allocate(i(n),x0y0z0(n))
    i=1
    cont=0
    if (nx>0) then
       do j=1,nx
          cont=cont+1
          x0y0z0(cont)=1
       enddo
    endif 
    if (ny>0) then
       do j=1,ny
          cont=cont+1
          x0y0z0(cont)=2
       enddo
    endif 
    if (nz>0) then
       do j=1,nz
          cont=cont+1
          x0y0z0(cont)=3
       enddo
    endif 

    do cont=1,4**n
       prod=1.0
       alpha=0
       beta=0
       gamma=0
       do m=1,n
          select case(i(m))
          case(1)
             alpha=alpha+1
             prod=prod*(rot(x0y0z0(m),1)*x(i(m))+rot(x0y0z0(m),2)*y(i(m)))
          case(2)
             beta=beta+1
             prod=prod*(rot(x0y0z0(m),1)*x(i(m))+rot(x0y0z0(m),2)*y(i(m)))
          case(3)
             gamma=gamma+1
             prod=prod*(rot(x0y0z0(m),1)*x(i(m))+rot(x0y0z0(m),2)*y(i(m)))
          case(4)
             prod=prod*rot(x0y0z0(m),3)*z
          end select
       enddo
       call area_integrate_triangular_coord(alpha,beta,gamma,aux_area_int)
       area_int=area_int+prod*aux_area_int
       i(n)=i(n)+1
       j=n
       do
          if((i(j)>4).and.(cont<4**n)) then
             i(j)=1
             i(j-1)=i(j-1)+1
             j=j-1
          else
             exit
          endif
       enddo
    enddo
    deallocate(i,x0y0z0)
 else
    area_int=1.0
 endif


end subroutine integrate_orderN_area_triangle














! This subroutine calculates the analytical surface integral of the gradient a N order polynomial term for a triangular face
! This term is expressed as: x**nx * y**ny * z**nz 
! nx, ny, nz are no negative integers 
subroutine integrate_orderN_area_triangle_grad(c,j,x0,nx,ny,nz,area_int)

use general_module


 ! Global variables
 integer, intent(in) :: nx,ny,nz, j, c
 double precision, intent(in), dimension(3,3) :: x0 ! x(1,1:3): X coordinates of the 3 triangle's nodes in global coordinates
                                                    ! x(2,1:3): Y coordinates of the 3 triangle's nodes in global coordinates
                                                    ! x(3,1:3): Z coordinates of the 3 triangle's nodes in global coordinates
 double precision, intent(out) :: area_int

 ! Local variables
 double precision ux, uy, uz, area_intx, area_inty, area_intz
 integer i



 ! Determination of the correct sense of surface normal vector
 ux=face(j)%norm(1,1)
 uy=face(j)%norm(2,1)
 uz=face(j)%norm(3,1)
 if (c/=face(j)%icell(1)) then
    ux=(-1)*ux
    uy=(-1)*uy
    uz=(-1)*uz
 endif

 ! Calculation of the gradient in X direction
 area_intx=0.0
 if ((ux/=0).and.(nx>0)) then
    call integrate_orderN_area_triangle(j,x0,nx-1,ny,nz,area_intx)
    area_intx=area_intx*nx*ux
 endif

 ! Calculation of the gradient in Y direction
 area_inty=0.0
 if ((uy/=0).and.(ny>0)) then
    call integrate_orderN_area_triangle(j,x0,nx,ny-1,nz,area_inty)
    area_inty=area_inty*ny*uy
 endif

 ! Calculation of the gradient in Z direction
 area_intz=0.0
 if ((uz/=0).and.(nz>0)) then
    call integrate_orderN_area_triangle(j,x0,nx,ny,nz-1,area_intz)
    area_intz=area_intz*nz*uz
 endif

 ! Addition of all the direcctions
 area_int=area_intx+area_inty+area_intz


end subroutine integrate_orderN_area_triangle_grad
















! This subroutine calculates the analytical volume integral of a N order polynomial term for a hexaedron cell
! This term is expressed as: x**nx * y**ny * z**nz 
! nx, ny, nz are no negative integers 
subroutine integrate_orderN_volume_hexaedron(c,nx,ny,nz,vol_int)

use general_module


 ! Global variables
 integer, intent(in) :: c, nx,ny,nz
 double precision, intent(out) :: vol_int

 ! Local variables
 double precision x(3,4), vol_tetraedron, aux_vol_int
 integer i,j,k, nodes_tetraedron(5,4), knode(4)


 ! The hexaedron is divided into 5 tetraedron and the integral is calculated by adding the integral of each tetraedron, 
 ! multiplied by its volume, and the sum is divided by the hexaedron's volume
 nodes_tetraedron(1,:)=(/1,2,3,6/)
 nodes_tetraedron(2,:)=(/1,5,8,6/)
 nodes_tetraedron(3,:)=(/1,3,8,6/)
 nodes_tetraedron(4,:)=(/7,3,8,6/)
 nodes_tetraedron(5,:)=(/1,3,8,4/)
 vol_int=0.0
 do i=1,5 ! Loop over the 5 tetraedra

    ! Definition of the relatives coordinates of each tetraedron with respect to the hexaedron center
    do k=1,4
       knode(k)=cell(c)%knode(nodes_tetraedron(i,k))
       do j=1,3
          x(j,k)=node(knode(k))%x(j)-cell(c)%x(j)
       enddo
    enddo

    ! Caculation of the volume of this tetraedron
    vol_tetraedron=abs((x(1,2)-x(1,1))*((x(2,3)-x(2,1))*(x(3,4)-x(3,1))-(x(2,4)-x(2,1))*(x(3,3)-x(3,1))) &
                   & -(x(2,2)-x(2,1))*((x(1,3)-x(1,1))*(x(3,4)-x(3,1))-(x(1,4)-x(1,1))*(x(3,3)-x(3,1))) &
                   & +(x(3,2)-x(3,1))*((x(1,3)-x(1,1))*(x(2,4)-x(2,1))-(x(1,4)-x(1,1))*(x(2,3)-x(2,1))))/6.0

    ! Integrate the polynomial term in this tetraedron
    call integrate_orderN_volume_tetraedron(x,nx,ny,nz,aux_vol_int)

    ! Multiply the integral and the volume of this tetraedron and add it to "vol_int"
    vol_int=vol_int + aux_vol_int*vol_tetraedron
 enddo


 ! Divide the sum by the hexaedron volume
 vol_int=vol_int/cell(c)%vol


end subroutine integrate_orderN_volume_hexaedron















! This subroutine calculates the analytical surface integral of a N order polynomial term and its gradient 
! for a cuadrangular face.
! This term is expressed as: x**nx * y**ny * z**nz 
! nx, ny, nz are no negative integers 
subroutine integrate_orderN_area_cuadrangle_and_grad(c,facej,nx,ny,nz,area_int,grad_int)

use general_module


 ! Global variables
 integer, intent(in) :: c,facej,nx,ny,nz
 double precision, intent(out) :: area_int, grad_int

 ! Local variables
 double precision, dimension(3,3) :: x0 ! x(1,1:3): X coordinates of the 3 triangle's nodes
                                                    ! x(2,1:3): Y coordinates of the 3 triangle's nodes
                                                    ! x(3,1:3): Z coordinates of the 3 triangle's nodes
 double precision area_int1, area_int2, area1, area2, normtriangle(3), v1(3), v2(3), grad_int1, grad_int2
 integer i,j,k


 ! Split the cuadrangle into 2 triangles: First one is composed of nodes 1, 2, 4 and the other is composed of
 ! nodes 2, 3, 4. 

 ! The integral is calculated by adding the integral of each triangle, 
 ! multiplied by its area, and the sum is divided by the cuadrangles's area which is the sum of the 2 triangle areas


 ! Calculation of the first triangle
 ! Definition of triangle relative coordinates with respecto to the hexaedron center
 do i=1,3
    x0(i,1)=node(face(facej)%knode(1))%x(i)-cell(c)%x(i)
    x0(i,2)=node(face(facej)%knode(2))%x(i)-cell(c)%x(i)
    x0(i,3)=node(face(facej)%knode(4))%x(i)-cell(c)%x(i)
 enddo
 ! Calculation of the area
 do i=1,3
    v1(i)=x0(i,2)-x0(i,1)
    v2(i)=x0(i,3)-x0(i,1)
 enddo
 normtriangle = cross_product(v1,v2)
 area1 = vector_magnitude(normtriangle)/2.0
 ! Integrate the polynomial term in this triangle
 call integrate_orderN_area_triangle(facej,x0,nx,ny,nz,area_int1)
 call integrate_orderN_area_triangle_grad(c,facej,x0,nx,ny,nz,grad_int1)


 ! Calculation of the second triangle
 ! Definition of triangle relative coordinates with respecto to the hexaedron center
 do i=1,3
    x0(i,1)=node(face(facej)%knode(2))%x(i)-cell(c)%x(i)
    x0(i,2)=node(face(facej)%knode(3))%x(i)-cell(c)%x(i)
    x0(i,3)=node(face(facej)%knode(4))%x(i)-cell(c)%x(i)
 enddo
 ! Calculation of the area
 do i=1,3
    v1(i)=x0(i,2)-x0(i,1)
    v2(i)=x0(i,3)-x0(i,1)
 enddo
 normtriangle = cross_product(v1,v2)
 area2 = vector_magnitude(normtriangle)/2.0
 ! Integrate the polynomial term in this triangle
 call integrate_orderN_area_triangle(facej,x0,nx,ny,nz,area_int2)
 call integrate_orderN_area_triangle_grad(c,facej,x0,nx,ny,nz,grad_int2)


 ! Calculation on the integral in the cuadrangle
 area_int=(area_int1*area1 + area_int2*area2)/(area1+area2)
 grad_int=(grad_int1*area1 + grad_int2*area2)/(area1+area2)


end subroutine integrate_orderN_area_cuadrangle_and_grad















! This subroutine calculates the cell averaged valaues multiplying each polynomial coeffcient for a tetraedron cell
! by using the following polynomial approach: 
! phi(x,y,z)= a1 * x**nx1 * y**ny1 * z**nz1 + a2 * x**nx2 * y**ny2 * z**nz2 + 
!           + a3 * x**nx3 * y**ny3 * z**nz3 + a4 * x**nx4 * y**ny4 * z**nz4 +
!           + a5 * x**nx5 * y**ny5 * z**nz5
! Only with 5 polynomial terms
subroutine vol_poly_coef_tetraedron_orderN(c,vect_coef)

use general_module


 ! Global variables
 integer c, vect_coef(3,5)

 ! Local variables
 double precision vol_poly_coef, x0(3,4)
 integer i, nx, ny, nz,j


 ! Definition of tetrahedron relative coordinates with respect to the tetrahedron center
 do i=1,3
    do j=1,4
       x0(i,j)=node(cell(c)%knode(j))%x(i)-cell(c)%x(i)
    enddo
 enddo


 ! Calculation on volume integrals
 do i=1,5
    nx=vect_coef(1,i)
    ny=vect_coef(2,i)
    nz=vect_coef(3,i)
    call integrate_orderN_volume_tetraedron(x0,nx,ny,nz,vol_poly_coef)
    vol_poly_coeff(c,i)=vol_poly_coef
 enddo


end subroutine vol_poly_coef_tetraedron_orderN















! This subroutine calculates the cell averaged valaues multiplying each polynomial coeffcient for a hexaedron cell
! by using the following polynomial approach: 
! phi(x,y,z)= a1 * x**nx1 * y**ny1 * z**nz1 + a2 * x**nx2 * y**ny2 * z**nz2 + 
!           + a3 * x**nx3 * y**ny3 * z**nz3 + a4 * x**nx4 * y**ny4 * z**nz4 +
!           + a5 * x**nx5 * y**ny5 * z**nz5 + a6 * x**nx6 * y**ny6 * z**nz6 +
!           + a7 * x**nx7 * y**ny7 * z**nz7
! Only with 7 polynomial terms
subroutine vol_poly_coef_hexaedron_orderN(c,vect_coef)

use general_module


 ! Global variables
 integer c, vect_coef(3,7)

 ! Local variables
 double precision vol_poly_coef
 integer i, nx, ny, nz



 ! Calculation on volume integrals
 do i=1,7
    nx=vect_coef(1,i)
    ny=vect_coef(2,i)
    nz=vect_coef(3,i)
    call integrate_orderN_volume_hexaedron(c,nx,ny,nz,vol_poly_coef)
    vol_poly_coeff(c,i)=vol_poly_coef
 enddo


end subroutine vol_poly_coef_hexaedron_orderN








! This subroutine calculates the surface and gradient averaged valaues multiplying each polynomial 
! coeffcient for a tetraedron cell by using the following polynomial approach: 
! phi(x,y,z)= a1 * x**nx1 * y**ny1 * z**nz1 + a2 * x**nx2 * y**ny2 * z**nz2 + 
!           + a3 * x**nx3 * y**ny3 * z**nz3 + a4 * x**nx4 * y**ny4 * z**nz4 +
!           + a5 * x**nx5 * y**ny5 * z**nz5
! Only with 5 polynomial terms
subroutine surf_and_grad_poly_coef_tetraedron_orderN(c,vect_coef)

use general_module


 ! Global variables
 integer c, vect_coef(3,5)

 ! Local variables
 integer j, jcell, i, k
 double precision surf_poly_coef, grad_poly_coef, x0(3,3)
 integer nx,ny,nz
    

 ! Calculate surface averaged values multipliers
 do jcell=1,ubound(cell(c)%jface,1)
    j=cell(c)%jface(jcell)

    ! Definition of triangle reative coordinates with respect to tetrahedron center
    do i=1,3
       do k=1,3
          x0(i,k)=node(face(j)%knode(k))%x(i)-cell(c)%x(i)
       enddo
    enddo

    ! Determination of surface and gradient multipliers for each surface
    do i=1,5
       nx=vect_coef(1,i)
       ny=vect_coef(2,i)
       nz=vect_coef(3,i)
       call integrate_orderN_area_triangle(j,x0,nx,ny,nz,surf_poly_coef)
       call integrate_orderN_area_triangle_grad(c,j,x0,nx,ny,nz,grad_poly_coef)
       surf_poly_coeff(c,jcell,i)=surf_poly_coef
       grad_poly_coeff(c,jcell,i)=grad_poly_coef
    enddo
 enddo


end subroutine surf_and_grad_poly_coef_tetraedron_orderN








! This subroutine calculates the surface and gradient averaged valaues multiplying each polynomial 
! coeffcient for a hexaedron cell by using the following polynomial approach: 
! phi(x,y,z)= a1 * x**nx1 * y**ny1 * z**nz1 + a2 * x**nx2 * y**ny2 * z**nz2 + 
!           + a3 * x**nx3 * y**ny3 * z**nz3 + a4 * x**nx4 * y**ny4 * z**nz4 +
!           + a5 * x**nx5 * y**ny5 * z**nz5 + a6 * x**nx6 * y**ny6 * z**nz6 +
!           + a7 * x**nx7 * y**ny7 * z**nz7
! Only with 7 polynomial terms
subroutine surf_and_grad_poly_coef_hexaedron_orderN(c,vect_coef)

use general_module


 ! Global variables
 integer c, vect_coef(3,7)

 ! Local variables
 integer j, jcell, i, k
 double precision surf_poly_coef, grad_poly_coef
 integer nx,ny,nz
    

 ! Calculate surface averaged values multipliers
 do jcell=1,ubound(cell(c)%jface,1)
    j=cell(c)%jface(jcell)


    ! Determination of surface and gradient multipliers for each surface
    do i=1,7
       nx=vect_coef(1,i)
       ny=vect_coef(2,i)
       nz=vect_coef(3,i)
       call integrate_orderN_area_cuadrangle_and_grad(c,j,nx,ny,nz,surf_poly_coef,grad_poly_coef)
       surf_poly_coeff(c,jcell,i)=surf_poly_coef
       grad_poly_coeff(c,jcell,i)=grad_poly_coef
    enddo
 enddo


end subroutine surf_and_grad_poly_coef_hexaedron_orderN















! This subroutine calculates the cell averaged valaues multiplying each polynomial coeffcient for a triangle cell
! by using the following polynomial approach: 
! phi(x,y,z)= a1 * x**nx1 * y**ny1 * z**nz1 + a2 * x**nx2 * y**ny2 * z**nz2 + 
!           + a3 * x**nx3 * y**ny3 * z**nz3 + a4 * x**nx4 * y**ny4 * z**nz4
! Only with 4 polynomial terms
subroutine vol_poly_coef_triangle_orderN(c,vect_coef)

use general_module


 ! Global variables
 integer c, vect_coef(3,4)

 ! Local variables
 double precision vol_poly_coef, x0(3,3)
 integer i, nx, ny, nz,j


 ! Definition of tetrahedron relative coordinates with respect to the tetrahedron center
 do i=1,3
    do j=1,3
       x0(i,j)=node(cell(c)%knode(j))%x(i)-cell(c)%x(i)
    enddo
 enddo


 ! Calculation on volume integrals
 do i=1,4
    nx=vect_coef(1,i)
    ny=vect_coef(2,i)
    nz=vect_coef(3,i)
    call integrate_orderN_volume_triangle(x0,nx,ny,nz,vol_poly_coef)
    vol_poly_coeff(c,i)=vol_poly_coef
 enddo


end subroutine vol_poly_coef_triangle_orderN














! This subroutine calculates the analytical volume integral of a N order polynomial term for a triangle cell
! This term is expressed as: x**nx * y**ny * z**nz 
! nx, ny, nz are no negative integers 
subroutine integrate_orderN_volume_triangle(x,nx,ny,nz,vol_int)

use general_module


 ! Global variables
 integer, intent(in) :: nx,ny,nz
 double precision, intent(in), dimension(3,3) :: x ! x(1,1:3): X coordinates of the 3 triangle's nodes
                                                   ! x(2,1:3): Y coordinates of the 3 triangle's nodes
                                                   ! x(3,1:3): Z coordinates of the 3 triangle's nodes
 double precision, intent(out) :: vol_int

 ! Local variables
 double precision prod, aux_vol_int
 integer n,m,j, alpha, beta, gamma, cont
 integer, allocatable, dimension(:) :: i, xyz


 ! Determine the polynomial order (n)
 n=nx+ny+nz

 ! Calculation of the volume integral
 if (n>0) then
    vol_int=0.0
    allocate(i(n),xyz(n))
    i=1
    cont=0
    if (nx>0) then
       do j=1,nx
          cont=cont+1
          xyz(cont)=1
       enddo
    endif 
    if (ny>0) then
       do j=1,ny
          cont=cont+1
          xyz(cont)=2
       enddo
    endif 
    if (nz>0) then
       do j=1,nz
          cont=cont+1
          xyz(cont)=3
       enddo
    endif 

    do cont=1,3**n
       prod=1.0
       alpha=0
       beta=0
       gamma=0
       do m=1,n
          select case(i(m))
          case(1)
             alpha=alpha+1
          case(2)
             beta=beta+1
          case(3)
             gamma=gamma+1
          end select
          prod=prod*x(xyz(m),i(m))
       enddo
       call area_integrate_triangular_coord(alpha,beta,gamma,aux_vol_int)
       vol_int=vol_int+prod*aux_vol_int
       i(n)=i(n)+1
       j=n
       do
          if((i(j)>3).and.(cont<3**n)) then
             i(j)=1
             i(j-1)=i(j-1)+1
             j=j-1
          else
             exit
          endif
       enddo
    enddo
    deallocate(i,xyz)
 else
    vol_int=1.0
 endif


end subroutine integrate_orderN_volume_triangle








! This subroutine calculates the surface and gradient averaged valaues multiplying each polynomial 
! coeffcient for a triangle cell by using the following polynomial approach: 
! phi(x,y,z)= a1 * x**nx1 * y**ny1 * z**nz1 + a2 * x**nx2 * y**ny2 * z**nz2 + 
!           + a3 * x**nx3 * y**ny3 * z**nz3 + a4 * x**nx4 * y**ny4 * z**nz4
! Only with 4 polynomial terms
subroutine surf_and_grad_poly_coef_triangle_orderN(c,vect_coef)

use general_module


 ! Global variables
 integer c, vect_coef(3,4)

 ! Local variables
 integer j, jcell, i, k
 double precision surf_poly_coef, grad_poly_coef, x0(3,2)
 integer nx,ny,nz
    

 ! Calculate surface averaged values multipliers
 do jcell=1,ubound(cell(c)%jface,1)
    j=cell(c)%jface(jcell)

    ! Definition of line relative coordinates with respect to triangle center
    do i=1,3
       do k=1,2
          x0(i,k)=node(face(j)%knode(k))%x(i)-cell(c)%x(i)
       enddo
    enddo

    ! Determination of surface and gradient multipliers for each surface
    do i=1,4
       nx=vect_coef(1,i)
       ny=vect_coef(2,i)
       nz=vect_coef(3,i)
       call integrate_orderN_area_line(j,x0,nx,ny,nz,surf_poly_coef)
       call integrate_orderN_area_line_grad(c,j,x0,nx,ny,nz,grad_poly_coef)
       surf_poly_coeff(c,jcell,i)=surf_poly_coef
       grad_poly_coeff(c,jcell,i)=grad_poly_coef
    enddo
 enddo


end subroutine surf_and_grad_poly_coef_triangle_orderN














! This subroutine calculates the analytical surface integral of a N order polynomial term for a line face
! This term is expressed as: x**nx * y**ny * z**nz 
! nx, ny, nz are no negative integers 
subroutine integrate_orderN_area_line(facej,x0,nx,ny,nz,area_int)

use general_module


 ! Global variables
 integer, intent(in) :: nx,ny,nz, facej
 double precision, intent(in), dimension(3,2) :: x0 ! x0(1,1:2): X coordinates of the 2 line's nodes in global coordinates
                                                    ! x0(2,1:2): Y coordinates of the 2 line's nodes in global coordinates
                                                    ! x0(3,1:2): Z coordinates of the 2 line's nodes in global coordinates
 double precision, intent(out) :: area_int

 ! Local variables
 double precision prod, aux_area_int, z, rot(3,3), x(2), y, ux(3), uy(3), uz(3), &
 & cos_x0_x, cos_x0_y, cos_x0_z, cos_y0_x, cos_y0_y, cos_y0_z, &
 & cos_z0_x, cos_z0_y, cos_z0_z
 integer n,m,j, alpha, beta, gamma, cont
 integer, allocatable, dimension(:) :: i, x0y0z0





 ! Definition of the global axes vectors
 ux(:)=(/1.0 , 0.0 , 0.0/)
 uy(:)=(/0.0 , 1.0 , 0.0/)
 uz(:)=(/0.0 , 0.0 , 1.0/)

 ! Calculation of angles between local and global axes
 cos_x0_x=dot_product(ux,face(facej)%norm(:,2))
 cos_x0_y=dot_product(ux,face(facej)%norm(:,3))
 cos_x0_z=dot_product(ux,face(facej)%norm(:,1))
 cos_y0_x=dot_product(uy,face(facej)%norm(:,2))
 cos_y0_y=dot_product(uy,face(facej)%norm(:,3))
 cos_y0_z=dot_product(uy,face(facej)%norm(:,1))
 cos_z0_x=dot_product(uz,face(facej)%norm(:,2))
 cos_z0_y=dot_product(uz,face(facej)%norm(:,3))
 cos_z0_z=dot_product(uz,face(facej)%norm(:,1))

 ! Definition of the rotation matrix
 rot(1,1:3)=(/cos_x0_x, cos_x0_y, cos_x0_z/);
 rot(2,1:3)=(/cos_y0_x, cos_y0_y, cos_y0_z/);
 rot(3,1:3)=(/cos_z0_x, cos_z0_y, cos_z0_z/);

  ! Calculation of local coordinates of line coordinates
 do j=1,2
    x(j)=cos_x0_x*x0(1,j) + cos_y0_x*x0(2,j) + cos_z0_x*x0(3,j)
 enddo
 y=cos_x0_y*x0(1,1) + cos_y0_y*x0(2,1) + cos_z0_y*x0(3,1)
 z=cos_x0_z*x0(1,1) + cos_y0_z*x0(2,1) + cos_z0_z*x0(3,1)




 ! Determine the polynomial order (n)
 n=nx+ny+nz

 ! Calculation of the surface integral
 if (n>0) then
    area_int=0.0
    allocate(i(n),x0y0z0(n))
    i=1
    cont=0
    if (nx>0) then
       do j=1,nx
          cont=cont+1
          x0y0z0(cont)=1
       enddo
    endif 
    if (ny>0) then
       do j=1,ny
          cont=cont+1
          x0y0z0(cont)=2
       enddo
    endif 
    if (nz>0) then
       do j=1,nz
          cont=cont+1
          x0y0z0(cont)=3
       enddo
    endif 

    do cont=1,3**n
       prod=1.0
       alpha=0
       beta=0
       do m=1,n
          select case(i(m))
          case(1)
             alpha=alpha+1
             prod=prod*rot(x0y0z0(m),1)*x(i(m))
          case(2)
             beta=beta+1
             prod=prod*rot(x0y0z0(m),1)*x(i(m))
          case(3)
             prod=prod*rot(x0y0z0(m),3)*z
          end select
       enddo
       call area_integrate_linear_coord(alpha,beta,aux_area_int)
       area_int=area_int+prod*aux_area_int
       i(n)=i(n)+1
       j=n
       do
          if((i(j)>3).and.(cont<3**n)) then
             i(j)=1
             i(j-1)=i(j-1)+1
             j=j-1
          else
             exit
          endif
       enddo
    enddo
    deallocate(i,x0y0z0)
 else
    area_int=1.0
 endif


end subroutine integrate_orderN_area_line














! This subroutine calculates the analytical surface integral of the gradient a N order polynomial term for a triangular face
! This term is expressed as: x**nx * y**ny * z**nz 
! nx, ny, nz are no negative integers 
subroutine integrate_orderN_area_line_grad(c,j,x0,nx,ny,nz,area_int)

use general_module


 ! Global variables
 integer, intent(in) :: nx,ny,nz, j, c
 double precision, intent(in), dimension(3,2) :: x0 ! x0(1,1:2): X coordinates of the 2 line's nodes in global coordinates
                                                    ! x0(2,1:2): Y coordinates of the 2 line's nodes in global coordinates
                                                    ! x0(3,1:2): Z coordinates of the 2 line's nodes in global coordinates
 double precision, intent(out) :: area_int

 ! Local variables
 double precision ux, uy, uz, area_intx, area_inty, area_intz
 integer i



 ! Determination of the correct sense of surface normal vector
 ux=face(j)%norm(1,1)
 uy=face(j)%norm(2,1)
 uz=face(j)%norm(3,1)
 if (c/=face(j)%icell(1)) then
    ux=(-1)*ux
    uy=(-1)*uy
    uz=(-1)*uz
 endif

 ! Calculation of the gradient in X direction
 area_intx=0.0
 if ((ux/=0).and.(nx>0)) then
    call integrate_orderN_area_line(j,x0,nx-1,ny,nz,area_intx)
    area_intx=area_intx*nx*ux
 endif

 ! Calculation of the gradient in Y direction
 area_inty=0.0
 if ((uy/=0).and.(ny>0)) then
    call integrate_orderN_area_line(j,x0,nx,ny-1,nz,area_inty)
    area_inty=area_inty*ny*uy
 endif

 ! Calculation of the gradient in Z direction
 area_intz=0.0
 if ((uz/=0).and.(nz>0)) then
    call integrate_orderN_area_line(j,x0,nx,ny,nz-1,area_intz)
    area_intz=area_intz*nz*uz
 endif

 ! Addition of all the direcctions
 area_int=area_intx+area_inty+area_intz


end subroutine integrate_orderN_area_line_grad















! This subroutine calculates the surface averaged valaues of surface integral for a linear face of the following
! polynomial term expressed in triangular coordinates: 
! p(t1,t2)= t1**n1 * t2**n2
! ti are the triangular coordinates and ni are the polynomial exponents
subroutine area_integrate_linear_coord(alpha,beta,area_int)

use general_module


 ! Global variables
 integer, intent(in) :: alpha,beta
 double precision, intent(out) :: area_int

 ! Local variables
 integer i, k, polynomial(2)
 double precision   p(3)


 ! Polynomial term
 polynomial(1)=alpha
 polynomial(2)=beta

 ! Calculation of the factorial products to determine the surface integral
 do k=1,2
    p(k)=1.0
    if (polynomial(k)/=0) then
       do i=1,polynomial(k)
          p(k)=p(k)*i
       enddo
    endif
 enddo
 k=3
 p(k)=1.0
 do i=1,sum(polynomial(:))+1
    p(k)=p(k)*i
 enddo    


 ! Calculate the surface integral
 area_int=p(1)*p(2)*1.0/p(3)
 
end subroutine area_integrate_linear_coord















! This subroutine calculates the cell averaged valaues multiplying each polynomial coeffcient for a cuadrangle cell
! by using the following polynomial approach: 
! phi(x,y,z)= a1 * x**nx1 * y**ny1 * z**nz1 + a2 * x**nx2 * y**ny2 * z**nz2 + 
!           + a3 * x**nx3 * y**ny3 * z**nz3 + a4 * x**nx4 * y**ny4 * z**nz4 +
!           + a5 * x**nx5 * y**ny5 * z**nz5 
! Only with 5 polynomial terms
subroutine vol_poly_coef_cuadrangle_orderN(c,vect_coef)

use general_module


 ! Global variables
 integer c, vect_coef(3,5)

 ! Local variables
 double precision vol_poly_coef
 integer i, nx, ny, nz



 ! Calculation on volume integrals
 do i=1,5
    nx=vect_coef(1,i)
    ny=vect_coef(2,i)
    nz=vect_coef(3,i)
    call integrate_orderN_volume_cuadrangle(c,nx,ny,nz,vol_poly_coef)
    vol_poly_coeff(c,i)=vol_poly_coef
 enddo


end subroutine vol_poly_coef_cuadrangle_orderN
















! This subroutine calculates the analytical volume integral of a N order polynomial term for a cuadrangle cell
! This term is expressed as: x**nx * y**ny * z**nz 
! nx, ny, nz are no negative integers 
subroutine integrate_orderN_volume_cuadrangle(c,nx,ny,nz,vol_int)

use general_module


 ! Global variables
 integer, intent(in) :: c, nx,ny,nz
 double precision, intent(out) :: vol_int

 ! Local variables
 double precision x(3,3), vol_triangle, aux_vol_int, vectnorm_triangle(3)
 integer i,j,k, nodes_triangle(2,3), knode(3)


 ! The cuadrangle is divided into 2 triangles and the integral is calculated by adding the integral of each triangle, 
 ! multiplied by its area, and the sum is divided by the cuadrangle's area
 nodes_triangle(1,:)=(/1,2,3/)
 nodes_triangle(2,:)=(/1,3,4/)
 vol_int=0.0
 do i=1,2 ! Loop over the 2 triangles

    ! Definition of the relatives coordinates of each triangle with respect to the cuarangle center
    do k=1,3
       knode(k)=cell(c)%knode(nodes_triangle(i,k))
       do j=1,3
          x(j,k)=node(knode(k))%x(j)-cell(c)%x(j)
       enddo
    enddo

    ! Caculation of the area of this triangle
    vectnorm_triangle(1)=(x(2,2)-x(2,1))*(x(3,3)-x(3,2))-(x(3,2)-x(3,1))*(x(2,3)-x(2,2))
    vectnorm_triangle(2)=(x(1,2)-x(1,1))*(x(3,3)-x(3,2))-(x(3,2)-x(3,1))*(x(1,3)-x(1,2))
    vectnorm_triangle(3)=(x(1,2)-x(1,1))*(x(2,3)-x(2,2))-(x(2,2)-x(2,1))*(x(1,3)-x(1,2))
    vol_triangle=sqrt(vectnorm_triangle(1)**2+vectnorm_triangle(2)**2+vectnorm_triangle(3)**2)/2.0

    ! Integrate the polynomial term in this triangle
    call integrate_orderN_volume_triangle(x,nx,ny,nz,aux_vol_int)

    ! Multiply the integral and the volume of this tetraedron and add it to "vol_int"
    vol_int=vol_int + aux_vol_int*vol_triangle
 enddo


 ! Divide the sum by the cuadrangle area
 vol_int=vol_int/cell(c)%vol


end subroutine integrate_orderN_volume_cuadrangle








! This subroutine calculates the surface and gradient averaged valaues multiplying each polynomial 
! coeffcient for a cuadrangle cell by using the following polynomial approach: 
! phi(x,y,z)= a1 * x**nx1 * y**ny1 * z**nz1 + a2 * x**nx2 * y**ny2 * z**nz2 + 
!           + a3 * x**nx3 * y**ny3 * z**nz3 + a4 * x**nx4 * y**ny4 * z**nz4 +
!           + a5 * x**nx5 * y**ny5 * z**nz5 
! Only with 5 polynomial terms
subroutine surf_and_grad_poly_coef_cuadrangle_orderN(c,vect_coef)

use general_module


 ! Global variables
 integer c, vect_coef(3,5)

 ! Local variables
 integer j, jcell, i, k
 double precision surf_poly_coef, grad_poly_coef, x0(3,2)
 integer nx,ny,nz
    

 ! Calculate surface averaged values multipliers
 do jcell=1,ubound(cell(c)%jface,1)
    j=cell(c)%jface(jcell)

    ! Definition of line relative coordinates with respect to cuadrangle center
    do i=1,3
       do k=1,2
          x0(i,k)=node(face(j)%knode(k))%x(i)-cell(c)%x(i)
       enddo
    enddo

    ! Determination of surface and gradient multipliers for each surface
    do i=1,5
       nx=vect_coef(1,i)
       ny=vect_coef(2,i)
       nz=vect_coef(3,i)
       call integrate_orderN_area_line(j,x0,nx,ny,nz,surf_poly_coef)
       call integrate_orderN_area_line_grad(c,j,x0,nx,ny,nz,grad_poly_coef)
       surf_poly_coeff(c,jcell,i)=surf_poly_coef
       grad_poly_coeff(c,jcell,i)=grad_poly_coef
    enddo
 enddo


end subroutine surf_and_grad_poly_coef_cuadrangle_orderN

















! This subroutine determines the ADFs for the adjoint cells of face "facej". The values
! calculated are "mult_adf1" and "mult_adf2" for cell 1 and 2 respectively
subroutine adf_cells(facej,mat1,mat2,mult_adf1,mult_adf2)

 use general_module
 use variables_xs_module

 ! Global variables
 integer, intent (in) :: facej, mat1, mat2
 double precision, intent(out) :: mult_adf1(2), mult_adf2(2)


 ! Local variables
 double precision :: uj(3)



 ! Define the norm vector of face "facej"
 uj=face(facej)%norm(1:3,1)  



 ! Apply ADF if it is specified in file "input"
 if (rotadf) then

    ! ADF is applied if the adjoint cells have different materials
    if (mat1==mat2) then

       ! ADF is applied in -X, +X,-Y, +Y
       if (abs(abs(uj(3))-1)>=1e-2) then

          ! -X, +X
          if (abs(abs(uj(1))-1)<=1e-2) then

             ! West ADF for cell 1 and east ADF for cell 2
             if (uj(1)<0) then
                mult_adf1(1)=mat_adf(4,1,mat1)
                mult_adf1(2)=mat_adf(4,2,mat1)
                mult_adf2(1)=mat_adf(2,1,mat2)
                mult_adf2(2)=mat_adf(2,2,mat2)

             ! East ADF for cell 1 and west ADF for cell 2
             else
                mult_adf1(1)=mat_adf(2,1,mat1)
                mult_adf1(2)=mat_adf(2,2,mat1)
                mult_adf2(1)=mat_adf(4,1,mat2)
                mult_adf2(2)=mat_adf(4,2,mat2)
             endif

          ! -Y, +Y
          elseif (abs(abs(uj(2))-1)<=1e-2) then

             ! South ADF for cell 1 and north ADF for cell 2
             if (uj(2)<0) then
                mult_adf1(1)=mat_adf(3,1,mat1)
                mult_adf1(2)=mat_adf(3,2,mat1)
                mult_adf2(1)=mat_adf(1,1,mat2)
                mult_adf2(2)=mat_adf(1,2,mat2)

             ! North ADF for cell 1 and south ADF for cell 2
             else
                mult_adf1(1)=mat_adf(1,1,mat1)
                mult_adf1(2)=mat_adf(1,2,mat1)
                mult_adf2(1)=mat_adf(3,1,mat2)
                mult_adf2(2)=mat_adf(3,2,mat2)
             endif

          else
             write(*,*) 'ERROR: ADF applied to a surface which is not -X, +X, -Y or +Y'
             stop
          endif
       else
          mult_adf1=1
          mult_adf2=1   
       endif
    else
       mult_adf1=1
       mult_adf2=1
    endif

 else
    mult_adf1=1
    mult_adf2=1
 endif

end subroutine adf_cells







!-----------------------------------------------------------------

end module equations_module

!----------------------------------------------------------------------------
