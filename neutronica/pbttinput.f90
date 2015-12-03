  subroutine pbttinput
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 20/9/2013
      use variables_xs_module
      implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! Reads in the PBTT x-section data both rodded and unrodded
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 24/9/2013
!      DIMENSION tb(48)
      double precision aa
      double precision, allocatable, dimension(:) :: tb
      integer i,j,k,l, ixsecur, ixsecr, ixsecur2, ixsecr2, ixsecur3, ixsecr3, &
      & icomp, nnn, icr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CHARACTER*80 dum
!
!
!
      PR_exist=.true.
!
      ixsecur=21
      ixsecr=22
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!
!!!!   Claubia - UFMG rmiro - UPV                                        !!!!
!
      ibor=0
      inquire (FILE='PPM_ranges', EXIST=PR_exist)
      if (PR_exist.eqv..true.) then
          ibor=1
!!!!   Claubia - UFMG rmiro - UPV      04/03/2010                        !!!!
         OPEN(UNIT=27, FILE='PPM_ranges', STATUS='OLD')
	   read(27,*)
	   read(27,*)ppml1
	   read(27,*)ppml2
	   read(27,*)ppml3
	   CLOSE(27)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ixsecur2=23
         ixsecr2=24
         OPEN(UNIT=ixsecur2, FILE='nemtab2', STATUS='OLD')
         OPEN(UNIT=ixsecr2, FILE='nemtabr2', STATUS='OLD')
!!!!
         ixsecur3=25
         ixsecr3=26
         OPEN(UNIT=ixsecur3, FILE='nemtab3', STATUS='OLD')
         OPEN(UNIT=ixsecr3, FILE='nemtabr3', STATUS='OLD')
      end if 
!!!!                                                                     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      rmiro - UPV                                                    !!!!
      OPEN(UNIT=ixsecur, FILE='nemtab', STATUS='OLD')
      OPEN(UNIT=ixsecr, FILE='nemtabr', STATUS='OLD')

!      OPEN(UNIT=ixsecur, FILE='PBTT_xsec_unrodded', STATUS='OLD')
!      OPEN(UNIT=ixsecr, FILE='PBTT_xsec_rodded', STATUS='OLD')
!!!!                                                                     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read the total number of compositions in x-section set
! Read the unrodded set first and rodded set next
!
! ncmpur13 is the number of unrodded compositions
!
!    1st do loop  reads these 5 lines
! *
! *  NEM-Cross Section Table Input                                                
! *                                                                               
! *    T Fuel        Density        Boron ppm.    T Mod.        Void              
!        6             6              0             0             0    
      do i=1,4
        READ(ixsecur,'(a80)')dum
      enddo

      READ(ixsecur,*)ntfentry,ndmentry
	ntentry=ntfentry+ndmentry
	ntdentry=(ntfentry*ndmentry)+(ntfentry+ndmentry)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 24/9/2013
   call ini_variables_pbttinput ! It initializes the variables that will be used
   allocate(tb(ntdentry))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! loop over all compositions
      do i = 1, ncmpur13
        READ(ixsecur,'(a80)')dum ! *
        READ(ixsecur,'(a80)')dum ! **** x-section set # ncmpur
        READ(ixsecur,*)icomp     ! read icomp
        do j=1,2
          READ(ixsecur,'(a80)')dum ! *
          READ(ixsecur,'(a80)')dum ! Group No. j 
          do k=1,5
            READ(ixsecur,'(a80)')dum ! *
            READ(ixsecur,'(a80)')dum ! ********
            READ(ixsecur,'(a80)')dum ! *
            if (j.eq.1 .and. k.eq.5) READ(ixsecur,'(a80)')dum ! *** from group 1 to 2
            READ(ixsecur,*) (tb(l), l=1,ntdentry)
            do l=1,ntdentry
              xspbttur(l,k,j,icomp) = tb(l)
            enddo
          enddo
          if (j.eq.2) then           ! read the 6th xs set of group 2: the Xe micro xs
            READ(ixsecur,'(a80)')dum ! *        
            READ(ixsecur,'(a80)')dum ! ********
            READ(ixsecur,'(a80)')dum ! *	    
            READ(ixsecur,*) (tb(l), l=1,ntdentry)
            do l=1,ntdentry
!             units for micro xs, needs 1.0e-24 to convert to Parcs units, jcr
              xspbttur(l,6,j,icomp) = tb(l) * 1.0e-24
            enddo
          endif
          do k=1,2                   ! read ADF's : pbtt final specs, jcr
            READ(ixsecur,'(a80)')dum ! *
            READ(ixsecur,'(a80)')dum ! ********
            READ(ixsecur,'(a80)')dum ! *	    
            READ(ixsecur,*) (tb(l), l=1,ntdentry)
            do l=1,ntdentry
              adfpbur(l,k,j,icomp) = tb(l)
            enddo
          enddo
! detector data
! k=1 Detector Flux Ratio Table, k=2 Detector Microscopic X-Section Table
          do k=1,2
             READ(ixsecur,'(a80)')dum ! *
             READ(ixsecur,'(a80)')dum ! ********
             READ(ixsecur,'(a80)')dum ! *	 
             READ(ixsecur,*) (tb(l), l=1,ntdentry)
             do l=1,ntdentry
                detpbur(l,k,j,icomp) = tb(l)
             enddo
           enddo
        enddo
        READ(ixsecur,'(a80)')dum ! *
        READ(ixsecur,'(a80)')dum ! ***** Effective Delayed Neutron Yield 
        READ(ixsecur,'(a80)')dum ! *
        READ(ixsecur,*)(betapbur(j,i),j=1,6)
        READ(ixsecur,'(a80)')dum ! *
        READ(ixsecur,'(a80)')dum ! ***** Decay Constants
        READ(ixsecur,'(a80)')dum ! *
        READ(ixsecur,*)(lambdapbur(j,i),j=1,6)
        READ(ixsecur,'(a80)')dum ! *
        READ(ixsecur,'(a80)')dum ! ***** inv velocities
        READ(ixsecur,'(a80)')dum ! *
        READ(ixsecur,*)(rvelpbur(j,i),j=1,2)
        READ(ixsecur,'(a80)')dum ! *
      enddo
      close(ixsecur)
!
! Rodded set next
!
! ncmpr13 is the number of rodded compositions
!
      do i=1,ncmpur13
        icrcomp(i)=0
      enddo
      do i=1,5
        READ(ixsecr,'(a80)')dum
      enddo
      do i = 1, ncmpr13
        READ(ixsecr,'(a80)')dum
        READ(ixsecr,'(a80)')dum
        READ(ixsecr,*)icr
        do j=1,2
          READ(ixsecr,'(a80)')dum
          READ(ixsecr,'(a80)')dum
          do k=1,5
            READ(ixsecr,'(a80)')dum
            READ(ixsecr,'(a80)')dum
            READ(ixsecr,'(a80)')dum
            if (j.eq.1 .and. k.eq.5) READ(ixsecr,'(a80)')dum
!
! icr are the compositions that are rodded eg 1-15, 34-48, ...
!
            icrcomp(icr)=i
            READ(ixsecr,*) (tb(l), l=1,ntdentry)
!
            do l=1,ntdentry
              xspbttr(l,k,j,i) = tb(l)
            enddo
          enddo
          if ( j.eq.2) then         ! read the 6th xs in group 2
            READ(ixsecr,'(a80)')dum ! *        
            READ(ixsecr,'(a80)')dum ! ********
            READ(ixsecr,'(a80)')dum ! *	    
            READ(ixsecr,*) (tb(l), l=1,ntdentry)
            do l=1,ntdentry
              xspbttr(l,6,j,i) = tb(l) * 1.0e-24
            enddo
          endif
          do k=1,2                  ! read ADF's : pbtt final specs, jcr
            READ(ixsecr,'(a80)')dum ! *
            READ(ixsecr,'(a80)')dum ! ********
            READ(ixsecr,'(a80)')dum ! *	    
            READ(ixsecr,*) (tb(l), l=1,ntdentry)
            do l=1,ntdentry
              adfpbr(l,k,j,i) = tb(l)
            enddo
          enddo
! detector data
! k=1 Detector Flux Ratio Table, k=2 Detector Microscopic X-Section Table
          do k=1,2
             READ(ixsecr,'(a80)')dum ! *
             READ(ixsecr,'(a80)')dum ! ********
             READ(ixsecr,'(a80)')dum ! *	 
             READ(ixsecr,*) (tb(l), l=1,ntdentry)
             do l=1,ntdentry
                detpbr(l,k,j,i) = tb(l)
             enddo
          enddo
        enddo
        READ(ixsecr,'(a80)')dum ! *
        READ(ixsecr,'(a80)')dum ! ***** Effective Delayed Neutron Yield 
        READ(ixsecr,'(a80)')dum ! *
        READ(ixsecr,*)(betapbr(j,i),j=1,6)
        READ(ixsecr,'(a80)')dum ! *
        READ(ixsecr,'(a80)')dum ! ***** Decay Constants
        READ(ixsecr,'(a80)')dum ! *
        READ(ixsecr,*)(lambdapbr(j,i),j=1,6)
        READ(ixsecr,'(a80)')dum ! *
        READ(ixsecr,'(a80)')dum ! ***** inv velocities
        READ(ixsecr,'(a80)')dum ! *
        READ(ixsecr,*)(rvelpbr(j,i),j=1,2)
        READ(ixsecr,'(a80)')dum ! *
      enddo
      close(ixsecr)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Claubia - UFMG rmiro - UPV                                        !!!!
      if(ibor.eq.1)then
!
!     reads nemtab2 and nemtabr2
!
      do i=1,4
        READ(ixsecur2,'(a80)')dum
      enddo
      READ(ixsecur2,*)nnn,nnn
!	ntentry=ntfentry+ndmentry
!	ntdentry=(ntfentry*ndmentry)+(ntfentry+ndmentry)

! loop over all compositions
      do i = 1, ncmpur13
        READ(ixsecur2,'(a80)')dum ! *
        READ(ixsecur2,'(a80)')dum ! **** x-section set # ncmpur
        READ(ixsecur2,*)icomp     ! read icomp
        do j=1,2
          READ(ixsecur2,'(a80)')dum ! *
          READ(ixsecur2,'(a80)')dum ! Group No. j 
          do k=1,5
            READ(ixsecur2,'(a80)')dum ! *
            READ(ixsecur2,'(a80)')dum ! ********
            READ(ixsecur2,'(a80)')dum ! *
            if (j.eq.1 .and. k.eq.5) READ(ixsecur2,'(a80)')dum ! *** from group 1 to 2
           READ(ixsecur2,*) (tb(l), l=1,ntdentry)
            do l=1,ntdentry
              xspbttur2(l,k,j,icomp) = tb(l)
            enddo
          enddo
          if (j.eq.2) then           ! read the 6th xs set of group 2: the Xe micro xs
            READ(ixsecur2,'(a80)')dum ! *        
            READ(ixsecur2,'(a80)')dum ! ********
            READ(ixsecur2,'(a80)')dum ! *	    
            READ(ixsecur2,*) (tb(l), l=1,ntdentry)
            do l=1,ntdentry
!             units for micro xs, needs 1.0e-24 to convert to Parcs units, jcr
              xspbttur2(l,6,j,icomp) = tb(l) * 1.0e-24
            enddo

          endif
          do k=1,2                   ! read ADF's : pbtt final specs, jcr
            READ(ixsecur2,'(a80)')dum ! *
            READ(ixsecur2,'(a80)')dum ! ********
            READ(ixsecur2,'(a80)')dum ! *	    
            READ(ixsecur2,*) (tb(l), l=1,ntdentry)
!            do l=1,ntdentry
!              adfpbur(l,k,j,icomp) = tb(l)
!            enddo
          enddo
! detector data
! k=1 Detector Flux Ratio Table, k=2 Detector Microscopic X-Section Table
          do k=1,2
             READ(ixsecur2,'(a80)')dum ! *
             READ(ixsecur2,'(a80)')dum ! ********
             READ(ixsecur2,'(a80)')dum ! *	 
             READ(ixsecur2,*) (tb(l), l=1,ntdentry)
!             do l=1,ntdentry
!                detpbur(l,k,j,icomp) = tb(l)
!             enddo
           enddo
        enddo
        READ(ixsecur2,'(a80)')dum ! *
        READ(ixsecur2,'(a80)')dum ! ***** Effective Delayed Neutron Yield 
        READ(ixsecur2,'(a80)')dum ! *
        READ(ixsecur2,*)(aa,j=1,6)
        READ(ixsecur2,'(a80)')dum ! *
        READ(ixsecur2,'(a80)')dum ! ***** Decay Constants
        READ(ixsecur2,'(a80)')dum ! *
        READ(ixsecur2,*)(aa,j=1,6)
        READ(ixsecur2,'(a80)')dum ! *
        READ(ixsecur2,'(a80)')dum ! ***** inv velocities
        READ(ixsecur2,'(a80)')dum ! *
        READ(ixsecur2,*)(aa,j=1,2)
        READ(ixsecur2,'(a80)')dum ! *
      enddo
      close(ixsecur2)
!
! Rodded set next
!
! ncmpr13 is the number of rodded compositions
!
      do i=1,ncmpur13
        icrcomp(i)=0
      enddo
      do i=1,5
        READ(ixsecr2,'(a80)')dum
      enddo
      do i = 1, ncmpr13
        READ(ixsecr2,'(a80)')dum
        READ(ixsecr2,'(a80)')dum
        READ(ixsecr2,*)icr
        do j=1,2
          READ(ixsecr2,'(a80)')dum
          READ(ixsecr2,'(a80)')dum
          do k=1,5
            READ(ixsecr2,'(a80)')dum
            READ(ixsecr2,'(a80)')dum
            READ(ixsecr2,'(a80)')dum
            if (j.eq.1 .and. k.eq.5) READ(ixsecr2,'(a80)')dum
!
! icr are the compositions that are rodded eg 1-15, 34-48, ...
!
            icrcomp(icr)=i
            READ(ixsecr2,*) (tb(l), l=1,ntdentry)
!
            do l=1,ntdentry
              xspbttr2(l,k,j,i) = tb(l)
            enddo
          enddo
          if ( j.eq.2) then         ! read the 6th xs in group 2
            READ(ixsecr2,'(a80)')dum ! *        
            READ(ixsecr2,'(a80)')dum ! ********
            READ(ixsecr2,'(a80)')dum ! *	    
            READ(ixsecr2,*) (tb(l), l=1,ntdentry)
            do l=1,ntdentry
              xspbttr2(l,6,j,i) = tb(l) * 1.0e-24
            enddo

          endif
          do k=1,2                  ! read ADF's : pbtt final specs, jcr
            READ(ixsecr2,'(a80)')dum ! *
            READ(ixsecr2,'(a80)')dum ! ********
            READ(ixsecr2,'(a80)')dum ! *	    
            READ(ixsecr2,*) (tb(l), l=1,ntdentry)
!            do l=1,ntdentry
!              adfpbr(l,k,j,i) = tb(l)
!            enddo
          enddo
! detector data
! k=1 Detector Flux Ratio Table, k=2 Detector Microscopic X-Section Table
          do k=1,2
             READ(ixsecr2,'(a80)')dum ! *
             READ(ixsecr2,'(a80)')dum ! ********
             READ(ixsecr2,'(a80)')dum ! *	 
             READ(ixsecr2,*) (tb(l), l=1,ntdentry)
!             do l=1,ntdentry
!                detpbr(l,k,j,i) = tb(l)
!             enddo
          enddo
        enddo
        READ(ixsecr2,'(a80)')dum ! *
        READ(ixsecr2,'(a80)')dum ! ***** Effective Delayed Neutron Yield 
        READ(ixsecr2,'(a80)')dum ! *
        READ(ixsecr2,*)(aa,j=1,6)
        READ(ixsecr2,'(a80)')dum ! *
        READ(ixsecr2,'(a80)')dum ! ***** Decay Constants
        READ(ixsecr2,'(a80)')dum ! *
        READ(ixsecr2,*)(aa,j=1,6)
        READ(ixsecr2,'(a80)')dum ! *
        READ(ixsecr2,'(a80)')dum ! ***** inv velocities
        READ(ixsecr2,'(a80)')dum ! *
        READ(ixsecr2,*)(aa,j=1,2)
        READ(ixsecr2,'(a80)')dum ! *
      enddo
      close(ixsecr2)
!
!     reads nemtab3 and nemtabr3
!
      do i=1,4
        READ(ixsecur3,'(a80)')dum
      enddo
      READ(ixsecur3,*)nnn,nnn
!	ntentry=ntfentry+ndmentry
!	ntdentry=(ntfentry*ndmentry)+(ntfentry+ndmentry)

! loop over all compositions
      do i = 1, ncmpur13
        READ(ixsecur3,'(a80)')dum ! *
        READ(ixsecur3,'(a80)')dum ! **** x-section set # ncmpur
        READ(ixsecur3,*)icomp     ! read icomp
        do j=1,2
          READ(ixsecur3,'(a80)')dum ! *
          READ(ixsecur3,'(a80)')dum ! Group No. j 
          do k=1,5
            READ(ixsecur3,'(a80)')dum ! *
            READ(ixsecur3,'(a80)')dum ! ********
            READ(ixsecur3,'(a80)')dum ! *
            if (j.eq.1 .and. k.eq.5) READ(ixsecur3,'(a80)')dum ! *** from group 1 to 2
           READ(ixsecur3,*) (tb(l), l=1,ntdentry)
            do l=1,ntdentry
              xspbttur3(l,k,j,icomp) = tb(l)
            enddo
          enddo
          if (j.eq.2) then           ! read the 6th xs set of group 2: the Xe micro xs
            READ(ixsecur3,'(a80)')dum ! *        
            READ(ixsecur3,'(a80)')dum ! ********
            READ(ixsecur3,'(a80)')dum ! *	    
            READ(ixsecur3,*) (tb(l), l=1,ntdentry)
            do l=1,ntdentry
!             units for micro xs, needs 1.0e-24 to convert to Parcs units, jcr
              xspbttur3(l,6,j,icomp) = tb(l) * 1.0e-24
            enddo

          endif
          do k=1,2                   ! read ADF's : pbtt final specs, jcr
            READ(ixsecur3,'(a80)')dum ! *
            READ(ixsecur3,'(a80)')dum ! ********
            READ(ixsecur3,'(a80)')dum ! *	    
            READ(ixsecur3,*) (tb(l), l=1,ntdentry)
!            do l=1,ntdentry
!              adfpbur(l,k,j,icomp) = tb(l)
!            enddo
          enddo
! detector data
! k=1 Detector Flux Ratio Table, k=2 Detector Microscopic X-Section Table
          do k=1,2
             READ(ixsecur3,'(a80)')dum ! *
             READ(ixsecur3,'(a80)')dum ! ********
             READ(ixsecur3,'(a80)')dum ! *	 
             READ(ixsecur3,*) (tb(l), l=1,ntdentry)
!             do l=1,ntdentry
!                detpbur(l,k,j,icomp) = tb(l)
!             enddo
           enddo
        enddo
        READ(ixsecur3,'(a80)')dum ! *
        READ(ixsecur3,'(a80)')dum ! ***** Effective Delayed Neutron Yield 
        READ(ixsecur3,'(a80)')dum ! *
        READ(ixsecur3,*)(aa,j=1,6)
        READ(ixsecur3,'(a80)')dum ! *
        READ(ixsecur3,'(a80)')dum ! ***** Decay Constants
        READ(ixsecur3,'(a80)')dum ! *
        READ(ixsecur3,*)(aa,j=1,6)
        READ(ixsecur3,'(a80)')dum ! *
        READ(ixsecur3,'(a80)')dum ! ***** inv velocities
        READ(ixsecur3,'(a80)')dum ! *
        READ(ixsecur3,*)(aa,j=1,2)
        READ(ixsecur3,'(a80)')dum ! *
      enddo
      close(ixsecur3)
!
! Rodded set next
!
! ncmpr13 is the number of rodded compositions
!
      do i=1,ncmpur13
        icrcomp(i)=0
      enddo
      do i=1,5
        READ(ixsecr3,'(a80)')dum
      enddo
      do i = 1, ncmpr13
        READ(ixsecr3,'(a80)')dum
        READ(ixsecr3,'(a80)')dum
        READ(ixsecr3,*)icr
        do j=1,2
          READ(ixsecr3,'(a80)')dum
          READ(ixsecr3,'(a80)')dum
          do k=1,5
            READ(ixsecr3,'(a80)')dum
            READ(ixsecr3,'(a80)')dum
            READ(ixsecr3,'(a80)')dum
            if (j.eq.1 .and. k.eq.5) READ(ixsecr3,'(a80)')dum
!
! icr are the compositions that are rodded eg 1-15, 34-48, ...
!
            icrcomp(icr)=i
            READ(ixsecr3,*) (tb(l), l=1,ntdentry)
!
            do l=1,ntdentry
              xspbttr3(l,k,j,i) = tb(l)
            enddo
          enddo
          if ( j.eq.2) then         ! read the 6th xs in group 2
            READ(ixsecr3,'(a80)')dum ! *        
            READ(ixsecr3,'(a80)')dum ! ********
            READ(ixsecr3,'(a80)')dum ! *	    
            READ(ixsecr3,*) (tb(l), l=1,ntdentry)
            do l=1,ntdentry
              xspbttr3(l,6,j,i) = tb(l) * 1.0e-24
            enddo

          endif
          do k=1,2                  ! read ADF's : pbtt final specs, jcr
            READ(ixsecr3,'(a80)')dum ! *
            READ(ixsecr3,'(a80)')dum ! ********
            READ(ixsecr3,'(a80)')dum ! *	    
            READ(ixsecr3,*) (tb(l), l=1,ntdentry)
!            do l=1,ntdentry
!              adfpbr(l,k,j,i) = tb(l)
!            enddo
          enddo
! detector data
! k=1 Detector Flux Ratio Table, k=2 Detector Microscopic X-Section Table
          do k=1,2
             READ(ixsecr3,'(a80)')dum ! *
             READ(ixsecr3,'(a80)')dum ! ********
             READ(ixsecr3,'(a80)')dum ! *	 
             READ(ixsecr3,*) (tb(l), l=1,ntdentry)
!             do l=1,ntdentry
!                detpbr(l,k,j,i) = tb(l)
!             enddo
          enddo
        enddo
        READ(ixsecr3,'(a80)')dum ! *
        READ(ixsecr3,'(a80)')dum ! ***** Effective Delayed Neutron Yield 
        READ(ixsecr3,'(a80)')dum ! *
        READ(ixsecr3,*)(aa,j=1,6)
        READ(ixsecr3,'(a80)')dum ! *
        READ(ixsecr3,'(a80)')dum ! ***** Decay Constants
        READ(ixsecr3,'(a80)')dum ! *
        READ(ixsecr3,*)(aa,j=1,6)
        READ(ixsecr3,'(a80)')dum ! *
        READ(ixsecr3,'(a80)')dum ! ***** inv velocities
        READ(ixsecr3,'(a80)')dum ! *
        READ(ixsecr3,*)(aa,j=1,2)
        READ(ixsecr3,'(a80)')dum ! *
      enddo
      close(ixsecr3)
!
      end if  ! (ibor)
!!!!                                                                     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 23/10/2013
      deallocate(tb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
 end subroutine pbttinput

