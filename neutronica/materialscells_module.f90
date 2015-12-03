module materialscells_module

implicit none
public materials_cells


integer :: number_cells
integer, allocatable, dimension (:) :: material_element



contains



subroutine materials_cells

use general_module
!use equations_module
use gmsh_module
integer :: gmesh_number ! output data concerned with this gmesh
character(len=4) :: centring ! whether to output only cell, face, none or all centred elements (default all)
integer :: error, i, j, k, l, m, region_number, kk, nnodes, n, mc, nrank, ns, &
  gregion, gnode, gelement, ngelements, ij, mvar, jj, number_points(5), &
  & num_ent_cells, type_cells(5), num_mat, num_region, mat 

double precision, dimension(:), allocatable :: cellvaluel
double precision, dimension(:,:), allocatable :: cellgradl
character(len=1000) :: filename, linkname, formatline, system_command
character(len=100) :: data_option,aux_number
integer, dimension(:), allocatable :: select_cells, select_faces, &
& num_region_material
logical, parameter :: select_elements = .false.
logical, parameter :: debug = .false.

gmesh_number=0


! Determine the number of materials
open(unit=95,file='input',status='old')
read(95,*)
read(95,*) num_mat ! number of materials
close (95)


! Determine "number_cells"
  if (gmesh(gmesh_number)%ngelements > 0) then
    num_ent_cells=0
    number_cells=0
    do gelement = 1, ubound(gmesh(gmesh_number)%gelement,1)
      i = gmesh(gmesh_number)%gelement(gelement)%icell 
      if (i /= 0.and.(gmesh(gmesh_number)%gelement(gelement)%JFACE == 0)) then
        num_ent_cells = num_ent_cells+number_points(cell(i)%gtype)+1
        number_cells=number_cells+1
      end if
    end do
  end if

! determine which regions are <material_> 
  allocate(num_region_material(num_mat))
  if (gmesh(gmesh_number)%ngregions > 0) then
    do gregion = 1, ubound(gmesh(gmesh_number)%gregion,1)
      region_number = gmesh(gmesh_number)%gregion(gregion)%region_number
      if (len_trim(region(region_number)%name)>=11) then
          if (region(region_number)%name(1:10) == '<material_') then
             i=scan(region(region_number)%name,'>') !!!! abernal, 20/11/2013
             read(region(region_number)%name(11:i-1),*) j !read(region(region_number)%name(11:10+i),*) j !!!! abernal, 20/11/2013
             num_region_material(j)=region_number
          endif
      endif
    end do
  end if

! Get the material distribution for each element (cell)
  allocate (material_element(number_cells)) 
  k=0
  if (gmesh(gmesh_number)%ngelements > 0) then
    do gelement = 1, ubound(gmesh(gmesh_number)%gelement,1)
      i = gmesh(gmesh_number)%gelement(gelement)%icell 
      if (i /= 0.and.(gmesh(gmesh_number)%gelement(gelement)%JFACE == 0)) then
          !
        do gregion = 1, ubound(gmesh(gmesh_number)%gelement(gelement)%gregions,1)
           j=0
           num_region=gmesh(gmesh_number)%gelement(gelement)%gregions(gregion)
           do
              j=j+1
              if (j<=num_mat) then
                  if (num_region==num_region_material(j)) then
                      k=k+1
                      material_element(k)=j
                      exit
                  endif
              else
                 exit
              endif
           enddo
        end do
      end if
    end do
  end if
  deallocate(num_region_material)


end subroutine materials_cells


end module materialscells_module 



