!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2011, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!     
!  SLEPc is free software: you can redistribute it and/or modify it under  the
!  terms of version 3 of the GNU Lesser General Public License as published by
!  the Free Software Foundation.
!
!  SLEPc  is  distributed in the hope that it will be useful, but WITHOUT  ANY 
!  WARRANTY;  without even the implied warranty of MERCHANTABILITY or  FITNESS 
!  FOR  A  PARTICULAR PURPOSE. See the GNU Lesser General Public  License  for 
!  more details.
!
!  You  should have received a copy of the GNU Lesser General  Public  License
!  along with SLEPc. If not, see <http://www.gnu.org/licenses/>.
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! ---------------------------------------------------------------------- 
!
module eigensolver_module


#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>
#include <finclude/petscvec.h90>

      double precision, allocatable :: autovalores(:), errores_autopares(:)
      double precision, allocatable :: autovectores(:,:), potencia(:,:), autovectores_adj(:,:), &
      & potencia_adj(:,:)
      double precision, dimension(:,:), allocatable :: autovectores_cel, autovectores_adj_cel
      integer :: num_autovalores, num_incognitas, adjoint



double precision, dimension(:), allocatable :: val_M11, val_M12, val_L21, & 
& val_L11, val_L22 
integer, dimension(:), allocatable :: row_L11_L22, col_L11_L22, row_M11_M12_L21, col_M11_M12_L21, &
 & trans_row_M11_M12_L21
! example: row_L11_L22=[1, 10, 15] 
!          Matrices L11 and L22 have 2 rows. The first row rangs from position 1 of vector "col_L11_L22"
!          to position 9 of vector "col_L11_L22". The second row rangs from position 10 to position 14
!          Vector "col_L11_L22" contains the positions of the eigenvector, that is, the column of the matrix.
!          Vector "val_L11" contains the value which will multiply the corresponding position of the eigenvector.
integer num_col_L11_L22, num_col_M11_M12_L21 ! length of array "col_L11_L22" and array "col_M11_M12_L21"




! PETSc variables
      Mat     L11, L22, L21, M11, M12, mx, L11t, L22t
      KSP     kL11, kL22, kL11t, kL22t
      Vec     w, w_adj



contains










subroutine eigen
implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/slepcsys.h>
#include <finclude/slepceps.h>
#include <finclude/petscvec.h90>

      integer i,j,ini_col,fin_col,default_param, number_unknowns,aux_i,cont
      integer, allocatable, dimension(:) :: rows_mx
      double precision, allocatable, dimension(:,:) :: orthonormal_basis



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!  Variables:

      Mat            Mshell, Mshell_adj
      EPS            eps, eps_adj
      EPSType        tname,method
      EPSWhich       which
      PetscReal      tol, error_eigenpairs
      PetscScalar    kr, ki
      PetscInt       n, npre
      PetscInt       nev, maxit, its,nconv, ncv, mpd, nconv_adj    
      PetscInt       col
      PetscMPIInt    rank
      PetscErrorCode ierr
      PetscBool      flg
      PetscScalar    value
      ST             st
      KSP            ksp, ksp_phi2
      PC             pc, pc_phi2
      Vec            xr, xi, xr2, b, xpi, bi, ui, x1_adji, x2_adji, uii
      PetscScalar, pointer :: auxiliar(:)
      MatFactorInfo  info(11) !!!! alv
      IS             rperm, cperm !!!! alv



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
      


      n=ubound(row_L11_L22,1)-1



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Compute the operator matrix that defines the eigensystem, Mshell*X = k*X
!                                         Mshell=L11_inv*(M11-M12*L22_inv*L21)
!                                                          
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      call MatCreate(PETSC_COMM_WORLD,L11,ierr)
      call MatSetSizes(L11,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
      call MatSetFromOptions(L11,ierr)
      npre=min(14,n)
      call MatSeqAIJSetPreallocation(L11,npre,PETSC_NULL_INTEGER,ierr)
      call MatSetUp(L11,ierr)

      call MatCreate(PETSC_COMM_WORLD,M11,ierr)
      call MatSetSizes(M11,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
      call MatSetFromOptions(M11,ierr)
      npre=min(7,n)
      call MatSeqAIJSetPreallocation(M11,npre,PETSC_NULL_INTEGER,ierr)
      call MatSetUp(M11,ierr)

      call MatCreate(PETSC_COMM_WORLD,M12,ierr)
      call MatSetSizes(M12,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
      call MatSetFromOptions(M12,ierr)
      npre=min(7,n)
      call MatSeqAIJSetPreallocation(M12,npre,PETSC_NULL_INTEGER,ierr)
      call MatSetUp(M12,ierr)

      call MatCreate(PETSC_COMM_WORLD,L21,ierr)
      call MatSetSizes(L21,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
      call MatSetFromOptions(L21,ierr)
      npre=min(7,n)
      call MatSeqAIJSetPreallocation(L21,npre,PETSC_NULL_INTEGER,ierr)
      call MatSetUp(L21,ierr)

      call MatCreate(PETSC_COMM_WORLD,L22,ierr)
      call MatSetSizes(L22,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
      call MatSetFromOptions(L22,ierr)
      npre=min(14,n)
      call MatSeqAIJSetPreallocation(L22,npre,PETSC_NULL_INTEGER,ierr)
      call MatSetUp(L22,ierr)








      ! Loop over the rows to define the matrices entries 
       cont=1
       do i=1,n
          aux_i=trans_row_M11_M12_L21(cont)      
          ini_col=row_L11_L22(i)
          fin_col=row_L11_L22(i+1)-1
          do j=ini_col,fin_col
             col=col_L11_L22(j)

             value=val_L11(j)
             call MatSetValues(L11,1,i-1,1,col-1,value,INSERT_VALUES,ierr)

             value=val_L22(j)
             call MatSetValues(L22,1,i-1,1,col-1,value,INSERT_VALUES,ierr)
          enddo
          if (i==aux_i) then
             ini_col=row_M11_M12_L21(cont)
             fin_col=row_M11_M12_L21(cont+1)-1
             do j=ini_col,fin_col
                col=col_M11_M12_L21(j)

                value=val_L21(j)
                call MatSetValues(L21,1,i-1,1,col-1,value,INSERT_VALUES,ierr)

                value=val_M11(j)
                call MatSetValues(M11,1,i-1,1,col-1,value,INSERT_VALUES,ierr)

                value=val_M12(j)
                call MatSetValues(M12,1,i-1,1,col-1,value,INSERT_VALUES,ierr)
             enddo
             cont=cont+1
             if (cont>ubound(trans_row_M11_M12_L21,1)) cont=1
          endif
       enddo







      call MatAssemblyBegin(L11,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(L11,MAT_FINAL_ASSEMBLY,ierr)

      call MatAssemblyBegin(M11,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(M11,MAT_FINAL_ASSEMBLY,ierr)

      call MatAssemblyBegin(M12,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(M12,MAT_FINAL_ASSEMBLY,ierr)

      call MatAssemblyBegin(L21,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(L21,MAT_FINAL_ASSEMBLY,ierr)

      call MatAssemblyBegin(L22,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(L22,MAT_FINAL_ASSEMBLY,ierr)






      ! Deallocate vectors which define the matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! abernal, 16/10/14,   acople
!      deallocate(val_M11,val_M12,val_L21,val_L11,val_L22,row_L11_L22, col_L11_L22)
      deallocate(val_M11,val_M12,val_L21,val_L11,val_L22,row_L11_L22,&
       & col_L11_L22,row_M11_M12_L21,col_M11_M12_L21,trans_row_M11_M12_L21)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Build KSP objects
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     call KSPcreate(PETSC_COMM_WORLD,kL11,ierr)
     call KSPSetOperators(kL11,L11,L11,DIFFERENT_NONZERO_PATTERN,ierr)
     call KSPSetType(kL11,KSPPREONLY,ierr)
     call KSPGetPC(kL11,pc,ierr)
     call PCSetType(pc,PCLU,ierr)
     !call KSPSetFromOptions(kL11,ierr)
      call PCFactorSetMatSolverPackage (pc,MATSOLVERMUMPS,ierr) 

     call KSPcreate(PETSC_COMM_WORLD,kL22,ierr)
     call KSPSetOperators(kL22,L22,L22,DIFFERENT_NONZERO_PATTERN,ierr)
     call KSPSetType(kL22,KSPPREONLY,ierr)
     call KSPGetPC(kL22,pc,ierr)
     call PCSetType(pc,PCLU,ierr)
     !call KSPSetFromOptions(kL22,ierr)
      call PCFactorSetMatSolverPackage (pc,MATSOLVERMUMPS,ierr) 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Create the eigensolver and display info
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!     ** Create eigensolver context
      call EPSCreate(PETSC_COMM_WORLD,eps,ierr)


!     Define Matrix shell: Mshell
      call MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,PETSC_NULL_INTEGER,Mshell,ierr)
      call MatShellSetOperation(Mshell,MATOP_MULT,obtain_matrix,ierr) ! "obtain_matrix" is the subroutine which calculates Mshell*phi_1

      call MatGetVecs(Mshell,w,PETSC_NULL_OBJECT,ierr)

!     ** Set operators
      call EPSSetOperators(eps,Mshell,PETSC_NULL_OBJECT,ierr) ! Standar eigenproblem
      call EPSSetProblemType(eps,EPS_NHEP,ierr) 

!     ** Set solver parameters
!     First the parameters are read from the file "input"
      open(unit=13,file='input',status='old')
      do i=1,5
          read(13,*)
      enddo
      read(13,*) default_param
!     If the option default is assigned:
      if (default_param==2) then
         nev=5
         ncv=2*nev
         mpd=ncv
      else
         read(13,*) 
         read(13,*) nev,ncv,mpd
      endif
      close(13)
!     Then the parameters are set
      call EPSSetDimensions (eps,nev,ncv,mpd,ierr) 
      which=EPS_LARGEST_MAGNITUDE
      call EPSSetWhichEigenpairs (eps,which,ierr)
      method=EPSKRYLOVSCHUR
      call EPSSetType (eps,method,ierr)
      call EPSSetFromOptions(eps,ierr)

!     ** Pasar las instrucciones de entrada por pantalla para que funcione bien SLEPC
      !value=1
      !call EPSSetTarget (eps,value,ierr)
      !call EPSGetST (eps,st,ierr)
      !call STSetType (st,STSINVERT,ierr)
      !call STGetKSP (st,ksp,ierr)
      !call KSPSetType (ksp,KSPPREONLY,ierr)
      !call KSPGetPC (ksp,pc,ierr)
      !call PCSetType (pc,PCLU,ierr)
      !call PCFactorSetMatSolverPackage (pc,MATSOLVERMUMPS,ierr) 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Solve the eigensystem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      call EPSSolve(eps,ierr) 
      call EPSGetIterationNumber(eps,its,ierr)
      if (rank .eq. 0) then
        write(*,*)
        write(*,*) 'Generalized Eigenvalue Problem'
        write(*,110) its
      endif
 110  format (/' Number of iterations of the method:',I4)
  
!     ** Optional: Get some information from the solver and display it
      call EPSGetType(eps,tname,ierr)
      if (rank .eq. 0) then
        write(*,120) tname
      endif
 120  format (' Solution method: ',A)
      call EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,                 &
     &                      PETSC_NULL_INTEGER,ierr)
      if (rank .eq. 0) then
        write(*,130) nev
      endif
 130  format (' Number of requested eigenvalues:',I2)
      call EPSGetTolerances(eps,tol,maxit,ierr)
      if (rank .eq. 0) then
        write(*,140) tol, maxit
      endif
 140  format (' Stopping condition: tol=',1P,E10.4,', maxit=',I4)

!     ** Obtain the number of converged eigenpairs
      call EPSGetConverged(eps,nconv,ierr)
      if (rank .eq. 0) then
        write(*,150) nconv
      endif
 150  format (' Number of converged eigenpairs:',I2/)






!      ** Get converged eigenpairs:
      num_incognitas=2*n
      num_autovalores=0
      if (nconv.gt.0) then
        num_autovalores=nconv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! abernal, 16/10/14,   acople
!        allocate (autovalores(nconv))
!        allocate (errores_autopares(nconv))
!        allocate (autovectores(num_incognitas,nconv))
        if (.not.allocated(autovalores)) allocate (autovalores(nconv))
        if (.not.allocated(errores_autopares)) allocate (errores_autopares(nconv))
        if (.not.allocated(autovectores)) allocate (autovectores(num_incognitas,nconv))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        call EPSPrintSolution(eps,PETSC_NULL_OBJECT,ierr)

        call VecDuplicate(w,xr,ierr)
        call VecDuplicate(w,xi,ierr)
        call VecDuplicate(w,xr2,ierr)
        call VecDuplicate(w,b,ierr)

        do i=0,nconv-1
          call EPSGetEigenpair(eps,i,kr,ki,xr,xi,ierr)
          autovalores(i+1)=PetscRealPart(kr)
          call VecGetArrayF90(xr,auxiliar,ierr)
          autovectores(1:num_incognitas/2,i+1)=auxiliar
          call VecRestoreArrayF90(xr,auxiliar,ierr)

          ! Calculate thermal flux (xr2): -L22*phi_2=L21*xr=b -->  {phi_2=-xr2} ---> L22*xr2=b  --> Solved by means of KSP
          call MatMult(L21,xr,b,ierr)
          call KSPSolve(kL22,b,xr2,ierr) !  L22*xr2=b

          call VecGetArrayF90(xr2,auxiliar,ierr)
          autovectores(1+num_incognitas/2:num_incognitas,i+1)=-1*auxiliar ! phi_2=-xr2
          call VecRestoreArrayF90(xr2,auxiliar,ierr)

!          ** Compute the relative error associated to each eigenpair
          call EPSComputeRelativeError(eps,i,error_eigenpairs,ierr)
          errores_autopares(i+1)=error_eigenpairs
        enddo







        ! ADJOINT CALCULATION
        if ((adjoint==1).and.(nconv.gt.0)) then

           ! Create the matrix whose columns are the fast flux for each eigenvalue (i.e. the eigenvectors) : matrix "mx"
           call MatCreate(PETSC_COMM_WORLD,mx,ierr)
           call MatSetSizes(mx,PETSC_DECIDE,PETSC_DECIDE,n,nconv,ierr)
           call MatSetFromOptions(mx,ierr)
           npre=min(nconv,n)
           call MatSeqAIJSetPreallocation(mx,npre,PETSC_NULL_INTEGER,ierr)
           call MatSetUp(mx,ierr)

           ! Loop over the eigenvalues to create "mx", which is an orthonormal basis of the eigenvectors
           ! Each column of "mx" will be each vector of this basis
           ! The orthonormalization is achieved by using the Gram-Schmidt process
           allocate(rows_mx(n),orthonormal_basis(n,nconv))
           do i=1,n
              rows_mx(i)=i-1
           enddo
           ! First eigenvector
           i=0
           ! Get the first eigenvector
           call EPSGetEigenpair(eps,i,kr,ki,xr,xi,ierr)
           ! Normalize the first eigenvector
           call VecNormalize(xr,PETSC_NULL_REAL,ierr)
           ! Add the first eigenvector to the first column of "mx"
           call VecGetArrayF90(xr,auxiliar,ierr)
           orthonormal_basis(:,i+1)=auxiliar
           call MatSetValues(mx,n,rows_mx,1,i,auxiliar,INSERT_VALUES,ierr)
           call VecRestoreArrayF90(xr,auxiliar,ierr)
           ! Continue the process for the rest eigenvectors
           if (nconv>1) then
              do i=1,nconv-1
                 ! Get the eigenvector "i"
                 call EPSGetEigenpair(eps,i,kr,ki,xr,xi,ierr)
                 ! Orthogonalize this eigenvector by using the Gram-Schmidt process
                 do j=0,i-1
                    ! Get the previous eigenvector orthonormalized : "j"
                    call VecSetValues(xr2,n,rows_mx,orthonormal_basis(:,j+1),INSERT_VALUES,ierr)
                    ! Calculate the dot product of the eigenvectors "i" and "j"
                    call VecDot(xr2,xr,value,ierr)
                    ! Calculate the sum: Vi=Vi-dot_product(Uj,Vi)*Uj
                    value=-1*value
                    call VecAXPY(xr,value,xr2,ierr)
                 enddo
                 ! Normalize the eigenvector "i"
                 call VecNormalize(xr,PETSC_NULL_REAL,ierr)
                 ! Add the eigenvector "i" to its corresponding column of the matrix
                 call VecGetArrayF90(xr,auxiliar,ierr)
                 orthonormal_basis(:,i+1)=auxiliar
                 call MatSetValues(mx,n,rows_mx,1,i,orthonormal_basis(:,i+1),INSERT_VALUES,ierr)
                 call VecRestoreArrayF90(xr,auxiliar,ierr)
              enddo
           endif
           deallocate(orthonormal_basis)
           ! Instructions to complete the creation of the matrix "mx"
           call MatAssemblyBegin(mx,MAT_FINAL_ASSEMBLY,ierr)
           call MatAssemblyEnd(mx,MAT_FINAL_ASSEMBLY,ierr)
           call VecAssemblyBegin(xr2,ierr)
           call VecAssemblyEnd(xr2,ierr)

           ! Create KSP context for L11' and L22'
           call MatTranspose(L11,MAT_INITIAL_MATRIX,L11t,ierr)
           call MatTranspose(L22,MAT_INITIAL_MATRIX,L22t,ierr)
           call KSPcreate(PETSC_COMM_WORLD,kL11t,ierr)
           call KSPSetOperators(kL11t,L11t,L11t,DIFFERENT_NONZERO_PATTERN,ierr)
           call KSPSetType(kL11t,KSPPREONLY,ierr)
           call KSPGetPC(kL11t,pc,ierr)
           call PCSetType(pc,PCLU,ierr)
           call PCFactorSetMatSolverPackage (pc,MATSOLVERMUMPS,ierr) 
           call KSPcreate(PETSC_COMM_WORLD,kL22t,ierr)
           call KSPSetOperators(kL22t,L22t,L22t,DIFFERENT_NONZERO_PATTERN,ierr)
           call KSPSetType(kL22t,KSPPREONLY,ierr)
           call KSPGetPC(kL22t,pc,ierr)
           call PCSetType(pc,PCLU,ierr)
           call PCFactorSetMatSolverPackage (pc,MATSOLVERMUMPS,ierr) 

           ! Create matrix "Mshell_adj": Mshell_adj=mx'*Mshell'*mx
           call MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,nconv,nconv,PETSC_NULL_INTEGER,Mshell_adj,ierr)
           call MatShellSetOperation(Mshell_adj,MATOP_MULT,obtain_matrix_adj,ierr) ! "obtain_matrix_adj" is the subroutine which calculates Mshell_adj*vector
           call MatGetVecs(Mshell_adj,w_adj,PETSC_NULL_OBJECT,ierr)

           ! Solve the following several eigenvalue problem: Mshell_adj*u=k*u  ,for "nconv" eigenvalues
           call EPSCreate(PETSC_COMM_WORLD,eps_adj,ierr)
           call EPSSetOperators(eps_adj,Mshell_adj,PETSC_NULL_OBJECT,ierr) ! Standar eigenproblem
           call EPSSetProblemType(eps_adj,EPS_NHEP,ierr) 
           call EPSSetDimensions (eps_adj,nconv,2*nconv,2*nconv,ierr) 
           which=EPS_LARGEST_MAGNITUDE
           call EPSSetWhichEigenpairs (eps_adj,which,ierr)
           method=EPSKRYLOVSCHUR
           call EPSSetType (eps_adj,method,ierr)
           call EPSSetFromOptions(eps_adj,ierr)
           call EPSSolve(eps_adj,ierr) 
           call EPSGetConverged(eps_adj,nconv_adj,ierr)
           print '(// A21)',' ADJOINT EIGENVALUES:'
           call EPSPrintSolution(eps_adj,PETSC_NULL_OBJECT,ierr) !!!! BORRAR

           ! Loop over the eigenvevalues calculated in the previous step, to calculate the adjoint fluxes:
           ! 1st : eigenvectors of   Mshell': xpi=mx*ui
           ! 2nd : L11'*x1_adji=xpi  --> solve with KSP
           ! 3rd : bi = 1/eigenvalue_i * M12' * x1_adji
           ! 4th : L22'*x2_adji=bi --> solve with KSP           call VecGetArrayF90(xr2,auxiliar,ierr)

           ! Define vectors and matrixes needed
           call VecDuplicate(w_adj,ui,ierr)
           call VecDuplicate(w_adj,uii,ierr)
           call VecDuplicate(xr,xpi,ierr)
           call VecDuplicate(xr,x1_adji,ierr)
           call VecDuplicate(xr,bi,ierr)
           call VecDuplicate(xr,x2_adji,ierr)
           ! Allocate "autovectores_adj"
           allocate (autovectores_adj(num_incognitas,nconv))
           ! Loop over the eigenvalues
           do i=0,nconv-1
              ! Get "ui"
              call EPSGetEigenpair(eps_adj,i,kr,ki,ui,uii,ierr)
              ! Calculate "xpi"
              call MatMult(mx,ui,xpi,ierr)
              ! Solve with KSP: L11'*x1_adji=xpi
              call KSPSolve(kL11t,xpi,x1_adji,ierr) 
              ! Calculate "bi" : bi = 1/eigenvalue_i * M12' * x1_adji
              call MatMultTranspose(M12,x1_adji,bi,ierr)
              value=1.0/autovalores(i+1)
              call VecScale(bi,value,ierr)
              ! Solve with KSP: L22'*x2_adji=bi
              call KSPSolve(kL22t,bi,x2_adji,ierr) 
              ! Get the fortran array of the adjoint fast flux
              call VecGetArrayF90(x1_adji,auxiliar,ierr)
              autovectores_adj(1:num_incognitas/2,i+1)=auxiliar 
              call VecRestoreArrayF90(x1_adji,auxiliar,ierr)
              ! Get the fortran array of the adjoint thermal flux
              call VecGetArrayF90(x2_adji,auxiliar,ierr)
              autovectores_adj(1+num_incognitas/2:num_incognitas,i+1)=auxiliar 
              call VecRestoreArrayF90(x2_adji,auxiliar,ierr)
           enddo
           ! Free work space  
           call EPSDestroy(eps_adj,ierr)
           call MatDestroy(mx,ierr)
           call VecDestroy(ui,ierr)
           call VecDestroy(xpi,ierr)
           call VecDestroy(x1_adji,ierr)
           call VecDestroy(x2_adji,ierr)
           call VecDestroy(bi,ierr)
           call KSPDestroy(kL11t,ierr)
           call KSPDestroy(kL22t,ierr)
           call MatDestroy(Mshell_adj,ierr)
        endif 
      endif     



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!     ** Free work space
      call EPSDestroy(eps,ierr)
      call KSPDestroy(kL11,ierr)
      call KSPDestroy(kL22,ierr)
      call MatDestroy(L11,ierr)
      call MatDestroy(M11,ierr)
      call MatDestroy(M12,ierr)
      call MatDestroy(L21,ierr)
      call MatDestroy(L22,ierr)
      call MatDestroy(Mshell,ierr)
      call VecDestroy(xr,ierr)
      call VecDestroy(xi,ierr)
      call VecDestroy(xr2,ierr)
      call VecDestroy(b,ierr)
      call VecDestroy(w,ierr)

!      call SlepcFinalize(ierr)

      
 end subroutine eigen









  
! This subroutine calculates the shell matrix needed to resolve the standard eigenvalue problem:
! Mshell*vector=[L11_inv*(M11-M12*L22_inv*L21)]*vector
 subroutine obtain_matrix(Mshell,initial,final,ierr)

implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscvec.h90>
#include <finclude/petscviewer.h>

     Mat            Mshell
     PetscErrorCode ierr
     Vec            initial, final
     PetscScalar    a

     call MatMult(L21,initial,final,ierr)
     call KSPSolve(kL22,final,w,ierr)
     call MatMult(M12,w,final,ierr)
     call MatMult(M11,initial,w,ierr)
     a=-1
     call VecAXPY(w,a,final,ierr)
     call KSPSolve(kL11,w,final,ierr)

 end subroutine obtain_matrix









  
! This subroutine calculates the shell matrix needed to resolve the standard eigenvalue problem:
! Mshell_adj*vector=[mx'*Mshell'*mx)]*vector= mx'*[(M11'-L21'*L22_inv'*M12')*L11_inv']*mx*vector
 subroutine obtain_matrix_adj(Mshell_adj,initial,final,ierr)

implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscvec.h90>
#include <finclude/petscviewer.h>

     Mat            Mshell_adj
     PetscErrorCode ierr
     Vec            initial, aux, aux2, aux3, final
     PetscScalar    a

     ! Define vectors
     call VecDuplicate(w,aux,ierr)
     call VecDuplicate(w,aux2,ierr)
     call VecDuplicate(w,aux3,ierr)

     ! Calculation of Mshell_adj*vector
     call MatMult(mx,initial,aux,ierr) 
     call KSPSolve(kL11t,aux,aux2,ierr)
     call MatMultTranspose(M11,aux2,aux3,ierr)
     call MatMultTranspose(M12,aux2,aux,ierr)
     call KSPSolve(kL22t,aux,aux2,ierr)
     call MatMultTranspose(L21,aux2,aux,ierr)
     a=-1
     call VecAYPX(aux,a,aux3,ierr)
     call MatMultTranspose(mx,aux,final,ierr)

 end subroutine obtain_matrix_adj







end module eigensolver_module

