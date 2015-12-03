!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ALVARO BERNAL GARCIA
! Phd Candidate
! ISIRYM, UNIVERSIDAD POLITECNICA DE VALENCIA
! VALENCIA, SPAIN
! 30/10/2013
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! NODES NUMBERING IN .geo FILE. EXAMPLE
!
!
! Reference axis:               Axial plane example:          Plane numbers: 2
!
!  y                                   * *
!   |                                * * * *
!   |__ __ x                         * * * *
!                                      * *
!
!
!
! Nodes numbering:
!
!      PLANE Z=1             PLANE Z=2
!    
!        1   2                13  14
!    3   4   5   6        15  16  17  18
!    7   8   9  10        19  20  21  22  
!       11  12                23  24 
!
!
! The node's name has to be: <material_#>        #: Node number






! This subroutine reads cross sections from files NEMTAB and assigns them to var(m),
! where m is the variable number corresponding to the diffusion equations constants.
 
 subroutine xsnemtab

     use variables_xs_module ! This module contains NEMTAB cross sections variables
     use general_module ! This module contains diffusion equations coefficients variables

     ! Variables declaration
     integer i, mat, num_mat, l, k, remainder, quotient
      

     ! Read the number of materials "num_mat" and the number of compositions 
     ! in NEMTAB files "ncmpur13" . These values are read from file "input"
     open(unit=95,file='input',status='old')
     read(95,*)
     read(95,*) num_mat ! number of materials
     do i=3,19 ! Line 20 contains "ncmpur13"
        read(95,*)
     enddo
     !"ncmpur13" and "ncmpr13" are defined in "variables_xs_module"
     read(95,*) ncmpur13 ! number of compositions in NEMTAB files
     ncmpr13=ncmpur13 ! 
     close (95)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  abernal, 18/10/14
 if (thcoupling==0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! This subroutine reads NEMTAB files
     call pbttinput


     ! This subroutine reads GEOM file, containing nodes geometry data
     call read_GEOM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  abernal, 18/10/14
 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! This subroutine determines the cross sections for a given fuel temperature, 
     !moderator density, boron concentration and control rod distribution.
     call xsecfb13 


     ! Initialize "mat_adf" if ADF are used
     if (rotadf) allocate(mat_adf(4,2,num_mat)) ! 4 ADF (n,e,s,w) and 2 energy groups



     ! Loop through all nodes to assign its cross sections. 
     ! (It is important to use the nodes numbering explained at the beginning)
     do mat=1,num_mat
        
        ! Calculation of auxiliar values to determine "l" and "k"
        quotient=mat/number_nodes_pz ! "number_nodes_pz" is defined in "variables_xs_module"
        remainder=mat-quotient*number_nodes_pz

        ! k: plane number
        ! l: node number on each plane
        if (remainder==0) then
           k=quotient
           l=number_nodes_pz
        else
           k=quotient+1
           l=remainder
        endif
       
        ! "var" is defined in "general_module"
        ! "xsd","xsa", "xss" and "xsnf" are defined in "variables_xs_module"
        var(1+7*(mat-1))%funk(1)%v=xsd(1,l,k) ! D_1
        var(2+7*(mat-1))%funk(1)%v=xsd(2,l,k) ! D_2
        var(3+7*(mat-1))%funk(1)%v=xsa(1,l,k) ! sigma_a_1
        var(4+7*(mat-1))%funk(1)%v=xsnf(1,l,k) ! nu_x_sigma_f_1
        var(5+7*(mat-1))%funk(1)%v=xss(l,k) ! sigma_s_12
        var(6+7*(mat-1))%funk(1)%v=xsnf(2,l,k) ! nu_x_sigma_f_2
        var(7+7*(mat-1))%funk(1)%v=xsa(2,l,k) ! sigma_a_2

        ! Determination of ADF for each material if necessary
        if (rotadf) then
           do i=1,2
              mat_adf(1,i,mat)=adfn(i,l,k)
              mat_adf(2,i,mat)=adfe(i,l,k)
              mat_adf(3,i,mat)=adfs(i,l,k)
              mat_adf(4,i,mat)=adfw(i,l,k)
           enddo
        endif
     enddo


     return
 end subroutine xsnemtab
