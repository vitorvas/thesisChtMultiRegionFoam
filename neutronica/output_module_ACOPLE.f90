! file src/output_module.f90
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
module output_module

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                        abernal, 4/06/2013
public write_vtk, calculo_potencia
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!-----------------------------------------------------------------
contains





subroutine write_vtk

! dumps the data to a file fit for printing

use general_module
use equations_module
use gmsh_module
use eigensolver_module
integer :: gmesh_number ! output data concerned with this gmesh
character(len=4) :: centring ! whether to output only cell, face, none or all centred elements (default all)
integer :: error, i, j, k, l, m, region_number, kk, nnodes, n, mc, nrank, ns, &
  gregion, gnode, gelement, ngelements, ij, mvar, jj, number_points(5), &
  & num_ent_cells, type_cells(5), number_cells, num_mat, num_region
double precision, dimension(:), allocatable :: cellvaluel
double precision, dimension(:,:), allocatable :: cellgradl
character(len=1000) :: filename, linkname, formatline, system_command
character(len=100) :: data_option,aux_number
integer, dimension(:), allocatable :: select_cells, select_faces, &
& num_region_material
logical, parameter :: select_elements = .false.
logical, parameter :: debug = .false.

gmesh_number=0
                  
if (debug) write(*,'(80(1h+)/a)') 'subroutine write_gmesh'

  
!----------------------------------------------------
! assemble filename based on options
filename = trim(gmesh(gmesh_number)%basename)
filename = trim(filename)//'.vtr'


! open file
open(foutput,file=trim(filename),status='replace',iostat=error)
if (error /= 0) call error_stop('problem opening file '//trim(filename))

! write vtk format intro stuff
write(foutput,'(a/a/a/a)') '# vtk DataFile Version 2.0', &
       &'output in VTK format','ASCII', 'DATASET UNSTRUCTURED_GRID'


! list nodes 
  if (gmesh(gmesh_number)%ngnodes > 0) then
    write(foutput,*) 'POINTS  ',gmesh(gmesh_number)%ngnodes,'  float'
    do gnode = 1, ubound(gmesh(gmesh_number)%knode_from_gnode,1)
      k = gmesh(gmesh_number)%knode_from_gnode(gnode)
      if (k /= 0) write(foutput,*) node(k)%x(1),' ',node(k)%x(2),' ',node(k)%x(3)
    end do
  end if

! list the nodes composing each cell
  number_points(1)=2
  number_points(2)=3
  number_points(3)=4
  number_points(4)=4
  number_points(5)=8
  if (gmesh(gmesh_number)%ngelements > 0) then
    num_ent_cells=0
    number_cells=0
    ! loop to determine "num_ent_cells"
    do gelement = 1, ubound(gmesh(gmesh_number)%gelement,1)
      i = gmesh(gmesh_number)%gelement(gelement)%icell 
      if (i /= 0.and.(gmesh(gmesh_number)%gelement(gelement)%JFACE == 0)) then
        num_ent_cells = num_ent_cells+number_points(cell(i)%gtype)+1
        number_cells=number_cells+1
      end if
    end do
    !
    write(foutput,*) 'CELLS   ', number_cells,'  ',num_ent_cells
    do gelement = 1, ubound(gmesh(gmesh_number)%gelement,1)
      i = gmesh(gmesh_number)%gelement(gelement)%icell 
      if (i /= 0.and.(gmesh(gmesh_number)%gelement(gelement)%JFACE == 0)) then
          write(foutput,*) number_points(cell(i)%gtype),' ', & ! set elementary entity equal to 
            (gmesh(gmesh_number)%gnode_from_knode(cell(i)%knode(l))-1,l=1,ubound(cell(i)%knode,1))
      end if
    end do
  end if

! list the types of each cell
  type_cells(1)=3
  type_cells(2)=5
  type_cells(3)=9
  type_cells(4)=10
  type_cells(5)=12
  if (gmesh(gmesh_number)%ngelements > 0) then
    write(foutput,*) 'CELL_TYPES   ', number_cells
    do gelement = 1, ubound(gmesh(gmesh_number)%gelement,1)
      i = gmesh(gmesh_number)%gelement(gelement)%icell 
      if (i /= 0.and.(gmesh(gmesh_number)%gelement(gelement)%JFACE == 0)) then
          write(foutput,*) type_cells(cell(i)%gtype)
      end if
    end do
  end if


! Determine the number of materials from the file "input"
open(unit=95,file='input',status='old')
read(95,*)
read(95,*) num_mat ! number of materials
close (95)


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


! write the material distribution
  if (gmesh(gmesh_number)%ngelements > 0) then
    write(foutput,*) 'CELL_DATA   ', number_cells
    write(foutput,*) 'SCALARS materiales int 1'
    write(foutput,*) 'LOOKUP_TABLE default'
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
                      write(foutput,*) j
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

 ! Ahora se escriben los flujos y potencias
 if (num_autovalores>=1) then
    write(foutput,*)
    write(foutput,*) 'FIELD fieldData ', 3*num_autovalores
    do i=1,num_autovalores
       write(aux_number,*) i
       write(foutput,*) 'flujo_1_autovalor_',trim(adjustl(aux_number)),&
                        &' 1 ',idomain,' float'
       do j=1,idomain
          write(foutput,*) autovectores_cel(j,i)
       enddo
       write(foutput,*) 'flujo_2_autovalor_',trim(adjustl(aux_number)),&
                        &' 1 ',idomain,' float'
       do j=1+idomain,2*idomain
          write(foutput,*) autovectores_cel(j,i)
       enddo
       write(foutput,*) 'potencia_autovalor_',trim(adjustl(aux_number)),&
                        &' 1 ',idomain,' float'
       do j=1,idomain
          write(foutput,*) potencia(j,i)
       enddo
    enddo
 endif

  close(foutput)

end subroutine write_vtk











subroutine calculo_potencia(number_cells_power)


use general_module
use eigensolver_module
use materialscells_module

integer, intent(out) :: number_cells_power
integer :: i, j, mat
double precision nu_x_sigma_f_1, nu_x_sigma_f_2


! Calculate the power: P = nu_x_sigma_f_1 * flux_1 + nu_x_sigma_f_2 * flux_2
 number_cells_power=idomain ! Número de celdas que tienen fisión
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! abernal, 16/10/14,   acople
! allocate (potencia(idomain,num_autovalores))
 if (.not.allocated(potencia)) allocate (potencia(idomain,num_autovalores))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do j=1,idomain
    mat=material_element(j)
    do i=1,num_autovalores
       nu_x_sigma_f_1 = var(4+7*(mat-1))%funk(1)%v
       ! flux_1 = autovectores_cel(j,i)
       nu_x_sigma_f_2 = var(6+7*(mat-1))%funk(1)%v
       ! flux_2 = autovectores_cel(j+idomain,i)
       potencia(j,i)=nu_x_sigma_f_1*autovectores_cel(j,i) + &
       & nu_x_sigma_f_2*autovectores_cel(j+idomain,i)
    enddo
    ! Se descuentan las celdas que no tienen fisión para normalizar después
    if ((nu_x_sigma_f_2==0).and.(nu_x_sigma_f_1==0)) number_cells_power=number_cells_power-1
 enddo


end subroutine calculo_potencia








subroutine nodes_results(postprocess,converged,power)  !!!! abernal, 15/10/2014


use general_module
use equations_module
use eigensolver_module
use materialscells_module

 ! Declaration of global variables
 integer postprocess, converged
 double precision, intent (out), dimension(:) :: power  !!!! abernal, 15/10/2014

 ! Declaration of local variables
 integer i, j, k, aux, nx, ny, nz, num_rc_reactor_nodes, &
 & num_materials, mat, reflector, i1, j1, k1
 integer, allocatable :: radmap(:,:), numbering_radmap(:,:,:), rc_reactor_nodes(:,:)
 double precision, allocatable :: vol_mat(:), p(:,:), f1(:,:), f2(:,:), &
 & sum_vol(:), sum_p_x_vol(:), assembly_power(:,:,:), axial_power(:,:), &
 & axial_vol(:), assembly_vol(:,:)




  ! Determine the number of materials from file "input"
  open(unit=13,file='input',status='old')
  read(13,*)
  read(13,*) num_materials



  ! Obtain other post-processing options from file "input"
  do i=3,23
     read(13,*)
  enddo
  read(13,*) nx, ny, nz
  read(13,*) 


  ! Obtain the radial map from file "input"
  allocate(radmap(ny,nx))
  do j=1,ny
     read(13,*) (radmap(j,i),i=1,nx)
  enddo
  close(13)



  ! Determine the row and colum of the reactor nodes
  aux=0
  reflector=0
  allocate(numbering_radmap(ny,nx,nz),rc_reactor_nodes(num_materials,3))
  do k=1,nz
     do j=1,ny
        do i=1,nx
           if (radmap(j,i)>0) then
              aux=aux+1
              rc_reactor_nodes(aux,1)=j
              rc_reactor_nodes(aux,2)=i
              rc_reactor_nodes(aux,3)=k
              numbering_radmap(j,i,k)=aux
           else
              numbering_radmap(j,i,k)=0
           end if
           if ((reflector==0).and.(radmap(j,i)==1)) reflector=1 ! signal to determine if there is reflector
        end do
     enddo
  enddo
  num_rc_reactor_nodes=aux



  ! Obtain the average results for each material (node)
  ! Initailize the post-processed variables
  allocate(vol_mat(num_materials))
  allocate(p(num_materials,num_autovalores))
  allocate(f1(num_materials,num_autovalores))
  allocate(f2(num_materials,num_autovalores))
  vol_mat=0.0
  p=0.0
  f1=0.0
  f2=0.0
  ! Loop pver the cells to calculate the nodes results
  do i=1,idomain
     mat=material_element(i)
     do j=1,num_autovalores
        p(mat,j)=p(mat,j)+potencia(i,j)*cell(i)%vol
        f1(mat,j)=f1(mat,j)+autovectores_cel(i,j)*cell(i)%vol
        f2(mat,j)=f2(mat,j)+autovectores_cel(i+idomain,j)*cell(i)%vol
     enddo
     vol_mat(mat)=vol_mat(mat)+cell(i)%vol
  enddo
  do mat=1,num_materials
     do j=1,num_autovalores
        p(mat,j)=p(mat,j)/vol_mat(mat)
        f1(mat,j)=f1(mat,j)/vol_mat(mat)
        f2(mat,j)=f2(mat,j)/vol_mat(mat)
     enddo
  enddo
  ! Renormalize the results to accomplish mean power is equal to 1. Same format as PARCS
  allocate(sum_p_x_vol(num_autovalores))
  sum_p_x_vol=0.0
  allocate (sum_vol(1))
  sum_vol(1)=0.0
  i=1
  do j=1,num_materials
     j1=rc_reactor_nodes(j,1)
     i1=rc_reactor_nodes(j,2)
     k1=rc_reactor_nodes(j,3)
     if ((radmap(j1,i1)==2).and.((reflector==0).or.((k1/=1).and.(k1/=nz)))) then
        sum_p_x_vol(i)=sum_p_x_vol(i)+p(j,i)*vol_mat(j)
        sum_vol(1)=sum_vol(1)+vol_mat(j)
     endif
  enddo
  if (num_autovalores>1) then
     do i=2,num_autovalores
        do j=1,num_materials
           j1=rc_reactor_nodes(j,1)
           i1=rc_reactor_nodes(j,2)
           k1=rc_reactor_nodes(j,3)
           if ((radmap(j1,i1)==2).and.((reflector==0).or.((k1/=1).and.(k1/=nz)))) then
              sum_p_x_vol(i)=sum_p_x_vol(i)+abs(p(j,i))*vol_mat(j)
           endif
        enddo
     enddo
  endif
  do i=1,num_autovalores
     do j=1,num_materials
        f1(j,i)=f1(j,i)/(sum_p_x_vol(i)/sum_vol(1))
        f2(j,i)=f2(j,i)/(sum_p_x_vol(i)/sum_vol(1))
        p(j,i)=p(j,i)/(sum_p_x_vol(i)/sum_vol(1))
     enddo
  enddo




  ! Write the nodes results in files "nodes_power.txt", "nodes_flux1.txt", "nodes_flux2.txt"
  if (converged==1) then
  open(unit=14,file='nodes_power.txt')
  open(unit=15,file='nodes_flux1.txt')
  open(unit=16,file='nodes_flux2.txt')
  do aux=1,num_autovalores
     write(14,*) 'EIGENVALUE ',aux
     write(15,*) 'EIGENVALUE ',aux
     write(16,*) 'EIGENVALUE ',aux
     do k=1,nz
        write(14,*) 
        write(15,*) 
        write(16,*) 
        write(14,*) 'PLANE Z=',k
        write(15,*) 'PLANE Z=',k
        write(16,*) 'PLANE Z=',k
        do j=1,ny
           write(14,*) 
           write(15,*) 
           write(16,*) 
           do i=1,nx
              if (numbering_radmap(j,i,k)/=0) then
                 mat=numbering_radmap(j,i,k)
                 write(14,fmt='(f10.4)',advance='no') p(mat,aux)
                 write(15,fmt='(f10.4)',advance='no') f1(mat,aux)
                 write(16,fmt='(f10.4)',advance='no') f2(mat,aux)
              else
                 write(14,fmt='(A10)',advance='no') '          '
                 write(15,fmt='(A10)',advance='no') '          ' 
                 write(16,fmt='(A10)',advance='no') '          ' 
              endif
           enddo
        enddo
        write(14,*) 
        write(15,*) 
        write(16,*) 
     enddo
     write(14,*) 
     write(15,*) 
     write(16,*)
     write(14,*) 
     write(15,*) 
     write(16,*)  
  enddo
  close(14)
  close(15)
  close(16)
  endif



  ! Assembly and axial power distributions
  if (postprocess==2) then
     ! Allocate power variables
     allocate(axial_power(nz,num_autovalores),assembly_power(ny,nx,num_autovalores))

     ! Axial power distribution
     allocate(axial_vol(nz))
     do k=1,nz
        do j=1,num_autovalores
           axial_power(k,j)=0.0
           if (j==1) axial_vol(k)=0.0
           do i=1+(k-1)*num_rc_reactor_nodes/nz,num_rc_reactor_nodes/nz*k
              j1=rc_reactor_nodes(i,1)
              i1=rc_reactor_nodes(i,2)
              if ((radmap(j1,i1)==2).and.((reflector==0).or.((k/=1).and.(k/=nz)))) then
                 axial_power(k,j)=axial_power(k,j)+p(i,j)*vol_mat(i)
                 if (j==1) axial_vol(k)=axial_vol(k)+vol_mat(i)
              endif
           enddo
           if (axial_vol(k)/=0) axial_power(k,j)=axial_power(k,j)/axial_vol(k)
        enddo
     enddo

     ! Normalize the axial power distribution. Same format as PARCS
     j=1
     sum_p_x_vol(j)=0.0
     sum_vol(1)=0.0
     do k=1,nz
        if (axial_vol(k)/=0) then
           sum_p_x_vol(j)=sum_p_x_vol(j)+axial_power(k,j)*axial_vol(k)
           sum_vol(1)=sum_vol(1)+axial_vol(k)
        endif
     enddo
     if (num_autovalores>1) then
        do j=2,num_autovalores
           sum_p_x_vol(j)=0.0
           do k=1,nz
              if (axial_vol(k)/=0) then
                 sum_p_x_vol(j)=sum_p_x_vol(j)+abs(axial_power(k,j))*axial_vol(k)
              endif
           enddo
        enddo
     endif
     do j=1,num_autovalores
        do k=1,nz
           axial_power(k,j)=axial_power(k,j)/(sum_p_x_vol(j)/sum_vol(1))
        enddo
     enddo

     ! Assembly power distribution
     allocate(assembly_vol(ny,nx))
     do j=1,ny
        do i=1,nx
           if (radmap(j,i)>0) then
              do aux=1,num_autovalores
                 assembly_power(j,i,aux)=0.0
                 if (aux==1) assembly_vol(j,i)=0.0
                 do k=1,nz
                    mat=numbering_radmap(j,i,k)
                    if ((radmap(j,i)==2).and.((reflector==0).or.((k/=1).and.(k/=nz)))) then
                       assembly_power(j,i,aux)=assembly_power(j,i,aux)+p(mat,aux)*vol_mat(mat)
                       if (aux==1) assembly_vol(j,i)=assembly_vol(j,i)+vol_mat(mat)
                    endif
                 enddo
                 if (assembly_vol(j,i)/=0) assembly_power(j,i,aux)=assembly_power(j,i,aux)/assembly_vol(j,i)
              enddo
           endif
        enddo
     enddo

     ! Normalize the assembly power distribution. Same format as PARCS
     aux=1
     sum_p_x_vol(aux)=0.0
     sum_vol(1)=0.0
     do j=1,ny
        do i=1,nx
           if (radmap(j,i)>0) then
              if (assembly_vol(j,i)/=0) then
                 sum_p_x_vol(aux)=sum_p_x_vol(aux)+assembly_power(j,i,aux)*assembly_vol(j,i)
                 sum_vol(1)=sum_vol(1)+assembly_vol(j,i)
              endif
           endif
        enddo
     enddo
     if (num_autovalores>1) then
        do aux=2,num_autovalores
           sum_p_x_vol(aux)=0.0
           do j=1,ny
              do i=1,nx
                 if (radmap(j,i)>0) then
                    if (assembly_vol(j,i)/=0) then
                       sum_p_x_vol(aux)=sum_p_x_vol(aux)+abs(assembly_power(j,i,aux))*assembly_vol(j,i)
                    endif
                 endif
              enddo
           enddo
        enddo
     endif
     do aux=1,num_autovalores
        do j=1,ny
           do i=1,nx
              if (radmap(j,i)>0) then
                 assembly_power(j,i,aux)=assembly_power(j,i,aux)/(sum_p_x_vol(aux)/sum_vol(1))
              endif
           enddo
        enddo
     enddo

     ! Write the axial and assembly power distributions in "power_distributions.txt"
     if (converged==1) then
     open(unit=17,file='power_distributions.txt')
     do aux=1,num_autovalores
        write(17,*) 'EIGENVALUE ',aux
        write(17,*)
        write(17,*) 'AXIAL POWER DISTRIBUTION'
        write(17,*)
        write(17,'(A23)') '  PLANE Z      POWER  '
        do k=1,nz
           write(17,'(i9,a2,f10.4)') (nz+1-k),'  ',axial_power(nz+1-k,aux)
        enddo
        write(17,*)
        write(17,*) 'ASSEMBLY POWER DISTRIBUTION'
        write(17,*)
        do j=1,ny
           do i=1,nx
              if (radmap(j,i)>0) then
                 write(17,fmt='(f10.4)',advance='no') assembly_power(j,i,aux)
              else
                 write(17,fmt='(A10)',advance='no') '          '
              endif
           enddo
           write(17,*)
        enddo
        write(17,*)
        write(17,*)
     enddo
     close(17)
     endif  
     deallocate(sum_vol,axial_power,assembly_power,axial_vol,assembly_vol,sum_p_x_vol)   
  endif



  ! Store the power in vector "power"  !!!! abernal, 15/10/2014
  power(1:num_materials)=p(1:num_materials,1)




  ! Deallocate variables
  deallocate(numbering_radmap,radmap,vol_mat,p,f1,f2)   


end subroutine nodes_results











subroutine calculo_potencia_adj(number_cells_power)


use general_module
use eigensolver_module
use materialscells_module

integer, intent(out) :: number_cells_power
integer :: i, j, mat
double precision nu_x_sigma_f_1, nu_x_sigma_f_2


! Calculate the power: P = nu_x_sigma_f_1 * flux_1 + nu_x_sigma_f_2 * flux_2
 number_cells_power=idomain ! Número de celdas que tienen fisión
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! abernal, 16/10/14,   acople
! allocate (potencia_adj(idomain,num_autovalores))
 if (.not.allocated(potencia)) allocate (potencia_adj(idomain,num_autovalores))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do j=1,idomain
    mat=material_element(j)
    do i=1,num_autovalores
       nu_x_sigma_f_1 = var(4+7*(mat-1))%funk(1)%v
       ! flux_1 = autovectores_adj_cel(j,i)
       nu_x_sigma_f_2 = var(6+7*(mat-1))%funk(1)%v
       ! flux_2 = autovectores_adj_cel(j+idomain,i)
       potencia_adj(j,i)=nu_x_sigma_f_1*autovectores_adj_cel(j,i) + &
       & nu_x_sigma_f_2*autovectores_adj_cel(j+idomain,i)
    enddo
    ! Se descuentan las celdas que no tienen fisión para normalizar después
    if ((nu_x_sigma_f_2==0).and.(nu_x_sigma_f_1==0)) number_cells_power=number_cells_power-1
 enddo


end subroutine calculo_potencia_adj








subroutine nodes_results_adj(postprocess)


use general_module
use equations_module
use eigensolver_module
use materialscells_module

 ! Declaration of global variables
 integer postprocess

 ! Declaration of local variables
 integer i, j, k, aux, nx, ny, nz, num_rc_reactor_nodes, &
 & num_materials, mat, reflector, i1, j1, k1
 integer, allocatable :: radmap(:,:), numbering_radmap(:,:,:), rc_reactor_nodes(:,:)
 double precision, allocatable :: vol_mat(:), p(:,:), f1(:,:), f2(:,:), &
 & sum_vol(:), sum_p_x_vol(:), assembly_power(:,:,:), axial_power(:,:), &
 & axial_vol(:), assembly_vol(:,:)




  ! Determine the number of materials from file "input"
  open(unit=13,file='input',status='old')
  read(13,*)
  read(13,*) num_materials



  ! Obtain other post-processing options from file "input"
  do i=3,23
     read(13,*)
  enddo
  read(13,*) nx, ny, nz
  read(13,*) 


  ! Obtain the radial map from file "input"
  allocate(radmap(ny,nx))
  do j=1,ny
     read(13,*) (radmap(j,i),i=1,nx)
  enddo
  close(13)



  ! Determine the row and colum of the reactor nodes
  aux=0
  reflector=0
  allocate(numbering_radmap(ny,nx,nz),rc_reactor_nodes(num_materials,3))
  do k=1,nz
     do j=1,ny
        do i=1,nx
           if (radmap(j,i)>0) then
              aux=aux+1
              rc_reactor_nodes(aux,1)=j
              rc_reactor_nodes(aux,2)=i
              rc_reactor_nodes(aux,3)=k
              numbering_radmap(j,i,k)=aux
           else
              numbering_radmap(j,i,k)=0
           end if
           if ((reflector==0).and.(radmap(j,i)==1)) reflector=1 ! signal to determine if there is reflector
        end do
     enddo
  enddo
  num_rc_reactor_nodes=aux



  ! Obtain the average results for each material (node)
  ! Initailize the post-processed variables
  allocate(vol_mat(num_materials))
  allocate(p(num_materials,num_autovalores))
  allocate(f1(num_materials,num_autovalores))
  allocate(f2(num_materials,num_autovalores))
  vol_mat=0.0
  p=0.0
  f1=0.0
  f2=0.0
  ! Loop pver the cells to calculate the nodes results
  do i=1,idomain
     mat=material_element(i)
     do j=1,num_autovalores
        p(mat,j)=p(mat,j)+potencia_adj(i,j)*cell(i)%vol
        f1(mat,j)=f1(mat,j)+autovectores_adj_cel(i,j)*cell(i)%vol
        f2(mat,j)=f2(mat,j)+autovectores_adj_cel(i+idomain,j)*cell(i)%vol
     enddo
     vol_mat(mat)=vol_mat(mat)+cell(i)%vol
  enddo
  do mat=1,num_materials
     do j=1,num_autovalores
        p(mat,j)=p(mat,j)/vol_mat(mat)
        f1(mat,j)=f1(mat,j)/vol_mat(mat)
        f2(mat,j)=f2(mat,j)/vol_mat(mat)
     enddo
  enddo
  ! Renormalize the results to accomplish mean power is equal to 1. Same format as PARCS
  allocate(sum_p_x_vol(num_autovalores))
  sum_p_x_vol=0.0
  allocate (sum_vol(1))
  sum_vol(1)=0.0
  i=1
  do j=1,num_materials
     j1=rc_reactor_nodes(j,1)
     i1=rc_reactor_nodes(j,2)
     k1=rc_reactor_nodes(j,3)
     if ((radmap(j1,i1)==2).and.((reflector==0).or.((k1/=1).and.(k1/=nz)))) then
        sum_p_x_vol(i)=sum_p_x_vol(i)+p(j,i)*vol_mat(j)
        sum_vol(1)=sum_vol(1)+vol_mat(j)
     endif
  enddo
  if (num_autovalores>1) then
     do i=2,num_autovalores
        do j=1,num_materials
           j1=rc_reactor_nodes(j,1)
           i1=rc_reactor_nodes(j,2)
           k1=rc_reactor_nodes(j,3)
           if ((radmap(j1,i1)==2).and.((reflector==0).or.((k1/=1).and.(k1/=nz)))) then
              sum_p_x_vol(i)=sum_p_x_vol(i)+abs(p(j,i))*vol_mat(j)
           endif
        enddo
     enddo
  endif
  do i=1,num_autovalores
     do j=1,num_materials
        f1(j,i)=f1(j,i)/(sum_p_x_vol(i)/sum_vol(1))
        f2(j,i)=f2(j,i)/(sum_p_x_vol(i)/sum_vol(1))
        p(j,i)=p(j,i)/(sum_p_x_vol(i)/sum_vol(1))
     enddo
  enddo




  ! Write the nodes results in files "nodes_power_adj.txt", "nodes_flux1_adj.txt", "nodes_flux2_adj.txt"
  open(unit=14,file='nodes_power_adj.txt')
  open(unit=15,file='nodes_flux1_adj.txt')
  open(unit=16,file='nodes_flux2_adj.txt')
  do aux=1,num_autovalores
     write(14,*) 'EIGENVALUE ',aux
     write(15,*) 'EIGENVALUE ',aux
     write(16,*) 'EIGENVALUE ',aux
     do k=1,nz
        write(14,*) 
        write(15,*) 
        write(16,*) 
        write(14,*) 'PLANE Z=',k
        write(15,*) 'PLANE Z=',k
        write(16,*) 'PLANE Z=',k
        do j=1,ny
           write(14,*) 
           write(15,*) 
           write(16,*) 
           do i=1,nx
              if (numbering_radmap(j,i,k)/=0) then
                 mat=numbering_radmap(j,i,k)
                 write(14,fmt='(f10.4)',advance='no') p(mat,aux)
                 write(15,fmt='(f10.4)',advance='no') f1(mat,aux)
                 write(16,fmt='(f10.4)',advance='no') f2(mat,aux)
              else
                 write(14,fmt='(A10)',advance='no') '          '
                 write(15,fmt='(A10)',advance='no') '          ' 
                 write(16,fmt='(A10)',advance='no') '          ' 
              endif
           enddo
        enddo
        write(14,*) 
        write(15,*) 
        write(16,*) 
     enddo
     write(14,*) 
     write(15,*) 
     write(16,*)
     write(14,*) 
     write(15,*) 
     write(16,*)  
  enddo
  close(14)
  close(15)
  close(16)



  ! Assembly and axial power distributions
  if (postprocess==2) then
     ! Allocate power variables
     allocate(axial_power(nz,num_autovalores),assembly_power(ny,nx,num_autovalores))

     ! Axial power distribution
     allocate(axial_vol(nz))
     do k=1,nz
        do j=1,num_autovalores
           axial_power(k,j)=0.0
           if (j==1) axial_vol(k)=0.0
           do i=1+(k-1)*num_rc_reactor_nodes/nz,num_rc_reactor_nodes/nz*k
              j1=rc_reactor_nodes(i,1)
              i1=rc_reactor_nodes(i,2)
              if ((radmap(j1,i1)==2).and.((reflector==0).or.((k/=1).and.(k/=nz)))) then
                 axial_power(k,j)=axial_power(k,j)+p(i,j)*vol_mat(i)
                 if (j==1) axial_vol(k)=axial_vol(k)+vol_mat(i)
              endif
           enddo
           if (axial_vol(k)/=0) axial_power(k,j)=axial_power(k,j)/axial_vol(k)
        enddo
     enddo

     ! Normalize the axial power distribution. Same format as PARCS
     j=1
     sum_p_x_vol(j)=0.0
     sum_vol(1)=0.0
     do k=1,nz
        if (axial_vol(k)/=0) then
           sum_p_x_vol(j)=sum_p_x_vol(j)+axial_power(k,j)*axial_vol(k)
           sum_vol(1)=sum_vol(1)+axial_vol(k)
        endif
     enddo
     if (num_autovalores>1) then
        do j=2,num_autovalores
           sum_p_x_vol(j)=0.0
           do k=1,nz
              if (axial_vol(k)/=0) then
                 sum_p_x_vol(j)=sum_p_x_vol(j)+abs(axial_power(k,j))*axial_vol(k)
              endif
           enddo
        enddo
     endif
     do j=1,num_autovalores
        do k=1,nz
           axial_power(k,j)=axial_power(k,j)/(sum_p_x_vol(j)/sum_vol(1))
        enddo
     enddo

     ! Assembly power distribution
     allocate(assembly_vol(ny,nx))
     do j=1,ny
        do i=1,nx
           if (radmap(j,i)>0) then
              do aux=1,num_autovalores
                 assembly_power(j,i,aux)=0.0
                 if (aux==1) assembly_vol(j,i)=0.0
                 do k=1,nz
                    mat=numbering_radmap(j,i,k)
                    if ((radmap(j,i)==2).and.((reflector==0).or.((k/=1).and.(k/=nz)))) then
                       assembly_power(j,i,aux)=assembly_power(j,i,aux)+p(mat,aux)*vol_mat(mat)
                       if (aux==1) assembly_vol(j,i)=assembly_vol(j,i)+vol_mat(mat)
                    endif
                 enddo
                 if (assembly_vol(j,i)/=0) assembly_power(j,i,aux)=assembly_power(j,i,aux)/assembly_vol(j,i)
              enddo
           endif
        enddo
     enddo

     ! Normalize the assembly power distribution. Same format as PARCS
     aux=1
     sum_p_x_vol(aux)=0.0
     sum_vol(1)=0.0
     do j=1,ny
        do i=1,nx
           if (radmap(j,i)>0) then
              if (assembly_vol(j,i)/=0) then
                 sum_p_x_vol(aux)=sum_p_x_vol(aux)+assembly_power(j,i,aux)*assembly_vol(j,i)
                 sum_vol(1)=sum_vol(1)+assembly_vol(j,i)
              endif
           endif
        enddo
     enddo
     if (num_autovalores>1) then
        do aux=2,num_autovalores
           sum_p_x_vol(aux)=0.0
           do j=1,ny
              do i=1,nx
                 if (radmap(j,i)>0) then
                    if (assembly_vol(j,i)/=0) then
                       sum_p_x_vol(aux)=sum_p_x_vol(aux)+abs(assembly_power(j,i,aux))*assembly_vol(j,i)
                    endif
                 endif
              enddo
           enddo
        enddo
     endif
     do aux=1,num_autovalores
        do j=1,ny
           do i=1,nx
              if (radmap(j,i)>0) then
                 assembly_power(j,i,aux)=assembly_power(j,i,aux)/(sum_p_x_vol(aux)/sum_vol(1))
              endif
           enddo
        enddo
     enddo

     ! Write the axial and assembly power distributions in "power_distributions_adj.txt"
     open(unit=17,file='power_distributions_adj.txt')
     do aux=1,num_autovalores
        write(17,*) 'EIGENVALUE ',aux
        write(17,*)
        write(17,*) 'AXIAL POWER DISTRIBUTION'
        write(17,*)
        write(17,'(A23)') '  PLANE Z      POWER  '
        do k=1,nz
           write(17,'(i9,a2,f10.4)') (nz+1-k),'  ',axial_power(nz+1-k,aux)
        enddo
        write(17,*)
        write(17,*) 'ASSEMBLY POWER DISTRIBUTION'
        write(17,*)
        do j=1,ny
           do i=1,nx
              if (radmap(j,i)>0) then
                 write(17,fmt='(f10.4)',advance='no') assembly_power(j,i,aux)
              else
                 write(17,fmt='(A10)',advance='no') '          '
              endif
           enddo
           write(17,*)
        enddo
        write(17,*)
        write(17,*)
     enddo
     close(17)  
     deallocate(sum_vol,axial_power,assembly_power,axial_vol,assembly_vol,sum_p_x_vol)   
  endif

  ! Deallocate variables
  deallocate(numbering_radmap,radmap,vol_mat,p,f1,f2)   






end subroutine nodes_results_adj










!-----------------------------------------------------------------

end module output_module

!-----------------------------------------------------------------