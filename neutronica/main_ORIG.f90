! ALVARO BERNAL GARCIA (abernal@iqn.upv.es)
! PhD Candidate
! ISISRYM, UNIVERSITAT POLITECNICA DE VALENCIA
! SPAIN
! 
!
!-------------------------------------------------------------------------
subroutine voldif(vector_valores) !!!!!program voldif ! abernal, 21/07/2014                 


use general_module
use setup_module
use equations_module
use output_module
use eigensolver_module


implicit none

real :: vector_valores(10)


integer nev ! number of eigenvalues to compute
integer ncv ! number of column vector to be used by the solution algorithm (largest dimension of the working subspace)
integer mpd ! maximum projected dimension used in the eigenproblem
integer default_param ! integer to activate the eigensolver: 0=No , 1=Yes , 2=Yes and default parameters
integer i,j, t, cont, tadj
double precision, dimension(:), allocatable :: sum_potencia_x_vol, sum_vol ! used to normalize the power
double precision, dimension(:), allocatable :: sum_potencia_adj_x_vol ! used to normalize the adjoint power
integer number_cells_power, postprocess


write(*,*) 'COMPROBACION LINKADO DE LIBRERIAS'
write(*,*) (vector_valores(i),i=1,10)


write(*,*)
write(*,*)
write(*,*) 'program VOLDIF 3D '
write(*,*)
write(*,*)
write(*,*) 'calculation of geometry properties: program arb'
write(*,*) 'eigensolver: SLEPc library'
write(*,*)
write(*,*)
write(*,*) '---------------------------------------------------------'
write(*,*)
write(*,*)




call time_process

call setup ! sets up variable metadata, reads in values, allocates arrays, creates mesh, initialises fields etc
call time_process(description='setup')


call update ! discretise equations
call time_process(description='equations discretization')






!!!! The eigensolver parameters are read from the file "input"
 open(unit=13,file='input',status='old')
 do i=1,5
    read(13,*)
 enddo
 read(13,*) default_param
!If the option default is assigned:
 if ((default_param==2).or.(default_param==4)) then
    nev=5
    ncv=2*nev
    mpd=ncv
 else
    read(13,*) 
    read(13,*) nev,ncv,mpd
 endif
 close(13)
 ! Determination of adjoint calculation
 select case(default_param)
 case(3)
    adjoint=1
    default_param=1
 case(4)
    adjoint=1
    default_param=2
 case default
    adjoint=0
 end select


if ((default_param==1).OR.(default_param==2)) then
  write(*,*)
  write(*,*) 'Eigenvalue Problem, calling eigenvalue solver...'
  call eigen

  call time_process(description='eigensolver')

  if (num_autovalores>=1) then
!       Se escriben los autovalores y los errores en el archivo
!       "autovalores.txt"
       open(unit=88,file='autovalores.txt')
       write(88,*) 'Autovalores:'
       write(88,*) (autovalores(i),i=1,num_autovalores)
       write(88,*)
       write(88,*) 'Errores:'
       write(88,*) (errores_autopares(i),i=1,num_autovalores)
       close(88)


! Calculo de los valores promediados en cada celda del flujo
       allocate(autovectores_cel(2*idomain,num_autovalores))
       do j=1,num_autovalores
          cont=0
          do i=1,idomain
             autovectores_cel(i,j)=0
             autovectores_cel(i+idomain,j)=0
             do t=1,ubound(cell(i)%jface,1)+1
                ! Flujo rapido
                cont=t+(ubound(cell(i)%jface,1)+1)*(i-1)
                autovectores_cel(i,j)=autovectores_cel(i,j)+vol_poly_coeff(i,t)*autovectores(cont,j)
                ! Flujo termico
                cont=t+(ubound(cell(i)%jface,1)+1)*(i+idomain-1)
                autovectores_cel(i+idomain,j)=autovectores_cel(i+idomain,j)+vol_poly_coeff(i,t)*autovectores(cont,j)
             enddo
          enddo
       enddo



! Calculo de los valores promediados en cada celda del flujo adjunto
       if (adjoint==1) then
          allocate(autovectores_adj_cel(2*idomain,num_autovalores))
          do j=1,num_autovalores
             cont=0
             do i=1,idomain
                autovectores_adj_cel(i,j)=0
                autovectores_adj_cel(i+idomain,j)=0
                do t=1,ubound(cell(i)%jface,1)+1
                   tadj=ubound(cell(i)%jface,1)+2-t
                   ! Flujo rapido
                   cont=t+(ubound(cell(i)%jface,1)+1)*(i-1)
                   if (t==1) autovectores_adj_cel(i,j)=autovectores_adj(cont,j)!!!!autovectores_adj_cel(i,j)=autovectores_adj_cel(i,j)+vol_poly_coeff(i,tadj)*autovectores_adj(cont,j)
                   ! Flujo termico
                   cont=t+(ubound(cell(i)%jface,1)+1)*(i+idomain-1)
                   if (t==1) autovectores_adj_cel(i+idomain,j)=autovectores_adj(cont,j)!!!!autovectores_adj_cel(i+idomain,j)=autovectores_adj_cel(i+idomain,j)+vol_poly_coeff(i,tadj)*autovectores_adj(cont,j)
                enddo
             enddo
          enddo
       endif
    


! Normalización de los flujos para que la potencia media sea igual 1
! Se calcula la potencia media pero con valores absolutos, y se divide la potencia y los flujos por la
! potencia media. Potencia media= sum(potencia*volumen)/sum(volumen) ,
! para los volumenes de las celdas en que la potencia no sea 0      
       call calculo_potencia(number_cells_power)  ! Calcula la potencia para todos los autovalores sin normalizar
       allocate(sum_potencia_x_vol(num_autovalores))
       sum_potencia_x_vol=0
       if (adjoint==1) then
          call calculo_potencia_adj(number_cells_power)
          allocate(sum_potencia_adj_x_vol(num_autovalores))
          sum_potencia_adj_x_vol=0
       endif
       allocate (sum_vol(num_autovalores))
       sum_vol=0
       i=1
          do j=1,idomain
             sum_potencia_x_vol(i)=sum_potencia_x_vol(i)+potencia(j,i)*cell(j)%vol ! Obtener la suma de la potencia de las celdas para todos los autovectores
             if (potencia(j,i) /= 0) sum_vol(i)=sum_vol(i)+cell(j)%vol 
          enddo
       if (adjoint==1) then
          do j=1,idomain
             sum_potencia_adj_x_vol(i)=sum_potencia_adj_x_vol(i)+potencia_adj(j,i)*cell(j)%vol ! Obtener la suma de la potencia de las celdas para todos los autovectores
          enddo
       endif
       if (num_autovalores>1) then
          do i=2,num_autovalores
             do j=1,idomain
                sum_potencia_x_vol(i)=sum_potencia_x_vol(i)+abs(potencia(j,i))*cell(j)%vol ! Obtener la suma de la potencia en valor absoluto de las celdas para todos los autovectores
                if (potencia(j,i) /= 0) sum_vol(i)=sum_vol(i)+cell(j)%vol 
             enddo
          enddo
          if (adjoint==1) then
             do i=2,num_autovalores
                do j=1,idomain
                   sum_potencia_adj_x_vol(i)=sum_potencia_adj_x_vol(i)+abs(potencia_adj(j,i))*cell(j)%vol ! Obtener la suma de la potencia de las celdas para todos los autovectores
                enddo
             enddo
          endif
       endif
       ! Multiplicar todos los autovectores y potencias por: 1/(sum_potencia_x_vol/sum_vol)
       do i=1,num_autovalores
          do j=1,idomain
             autovectores_cel(j,i)=autovectores_cel(j,i)/(sum_potencia_x_vol(i)/sum_vol(i))
             autovectores_cel(j+idomain,i)=autovectores_cel(j+idomain,i)/(sum_potencia_x_vol(i)/sum_vol(i))
             potencia(j,i)=potencia(j,i)/(sum_potencia_x_vol(i)/sum_vol(i))
          enddo
       enddo
       if (adjoint==1) then
          do i=1,num_autovalores
             do j=1,idomain
                autovectores_adj_cel(j,i)=autovectores_adj_cel(j,i)/(sum_potencia_adj_x_vol(i)/sum_vol(i))
                autovectores_adj_cel(j+idomain,i)=autovectores_adj_cel(j+idomain,i)/(sum_potencia_adj_x_vol(i)/sum_vol(i))
                potencia_adj(j,i)=potencia_adj(j,i)/(sum_potencia_adj_x_vol(i)/sum_vol(i))
             enddo
          enddo
       endif                                       

!      Se escriben los autovectores en el archivo "autovectores.txt"
!      Cada columna corresponde a cada autovector
       open(unit=89,file='autovectores.txt')
       do j=1,2*idomain
          write(89,*) (autovectores_cel(j,i),i=1,num_autovalores)
       enddo
       close(89)

!      Se escriben los autovectores del adjunto en el archivo "autovectores_adj.txt"
!      Cada columna corresponde a cada autovector
       if (adjoint==1) then
          open(unit=89,file='autovectores_adj.txt')
          do j=1,2*idomain
             write(89,*) (autovectores_adj_cel(j,i),i=1,num_autovalores)
          enddo
          close(89)
       endif 


!      Se libera la memoria de "sum_potencia_x_vol" y "sum_vol"
       deallocate(sum_potencia_x_vol,sum_vol)
       if (adjoint==1) then
          deallocate(sum_potencia_adj_x_vol)
       endif 

       call write_vtk

       ! Post-processing: Nodes results and axial and radial distributions
       open(unit=13,file='input',status='old')
       do i=1,21
          read(13,*)
       enddo
       read(13,*) postprocess
       close(13)
       if (postprocess/=0) call nodes_results(postprocess)
       if ((postprocess/=0).and.(adjoint==1)) call nodes_results_adj(postprocess)


       call time_process(description='post-processing: normalization and writing in VTK format')

       write(*,'(a)') "SUCCESS: the simulation finished gracefully"
  else
     write(*,*) 'ERROR: problem in eigensolver since no eigenvalues converged'
  endif
else
   write(*,*) 'ERROR: eigenvalue was not calculated. Check line 6 of file "input"'
endif




end subroutine voldif !!!!end program voldif ! abernal, 23/07/2014

!-----------------------------------------------------------------
