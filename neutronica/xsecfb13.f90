 subroutine xsecfb13
!
! calculate nodal cross sections at given t/h condition per PBTT
! Benchmark spec
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 14/10/2013
    use variables_xs_module
!    real, allocatable :: DED (:,:,:) !!!! abernal
    double precision tlap, tfmax, txsec, Tfunif, Dmunif, Tmunif, tval, dval
    integer k, kth, lth, l, ipr, la, icomp, ii, xs_reading, err, dim_zcr, icr
    double precision, dimension (:,:), allocatable :: dcool, tdopl, rambda, &
    & crbdens, tcool
    logical fdbk, rffdbk
    double precision, dimension(:), allocatable :: z_cr
    double precision, dimension(2):: rodfrac  ! 2: number of groups
    integer reactor_type ! 0:PWR, 1:BWR
    integer num_nodos_comb
    integer,  dimension (:), allocatable :: nodos_comb


    !txsec=0 !!!! abernal
    timeedt=-1.0
    eps=0.001

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 18/10/2014
!    call ini_variables_pbttxsec ! To allocate the variables that will store the cross sections
    if (thcoupling==0) call ini_variables_pbttxsec ! To allocate the variables that will store the cross sections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 18/10/2014
!    rotadf=.false. !!!! abernal, 4/7/2014
    if (thcoupling==0) rotadf=.false. !!!! abernal, 4/7/2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 18/10/2014
    if (thcoupling==1) fdbk=.true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! Instruction to determine if the thermalhydraulic variables will be read from
! the files "DENS" AND "TFUS" OR FROM "input"  !!!! abernal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 14/10/2013
    open(unit=8,file='input',status='OLD')
    do ii=1,11  ! The line 12 of file "input" contains the integer ("xs_reading") to determine the XS reading
       read(8,*)
    end do
    read(8,*,iostat=err) xs_reading
    ! close(8) !!!! abernal, 17/10/2014

    if (err==0) then
       select case (xs_reading) ! 2:3D thermalhydraulic data, 3:uniform thermalhydraulic data
       case (2,4)
          fdbk=.false.
          rffdbk=.true.
          if (xs_reading==4) rotadf=.true. !!!! abernal, 4/7/2014
       case (3,5)
          fdbk=.false.
          rffdbk=.false.
          if (xs_reading==5) rotadf=.true. !!!! abernal, 4/7/2014
          read(8,*)
          read(8,*,iostat=err) Dmunif, Tfunif, Tmunif
          if (err/=0) then
             print '(///)'
             print *,'ERROR READING THE LINE 14 OF FILE "input"'
             print '(///)'
             stop
          end if
       case (6,7) !!!! abernal, 30/10/2014
          fdbk=.true. !!!! abernal, 30/10/2014
          rffdbk=.true. !!!! abernal, 30/10/2014
          if (xs_reading==7) rotadf=.true. !!!! abernal, 30/10/2014
       end select
    else
       print '(///)'
       print *,'ERROR READING THE LINE 12 OF FILE "input"'
       print '(///)'
       stop
    end if
    close(8) !!!! abernal, 17/10/2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




! Instructions to determine Xenon Concentration !!!! abernal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 12/11/2013
 if (thcoupling==0) then !!!!  abernal, ISIRYM, 18/10/2014
    open(unit=9,file='3D_Xenon_Number_Density.txt',status='OLD',iostat=err)
    rnxe=0.0 ! Toda la matriz se iguala a 0
    if (err==0) then ! Se ha encontrado el fichero
       ! Primero hay que saber que nodos del plano radial son combustible
       num_nodos_comb=0
       allocate(nodos_comb(number_nodes_pz))
       do l=1,number_nodes_pz
          if (materials_nodes(2,l)<ncmpur13-2) then ! Las 3 ultimas composiciones son las del reflector
             num_nodos_comb=num_nodos_comb+1
             nodos_comb(num_nodos_comb)=l
          end if
       end do
       ! Se lee la concentracion del Xenon
       do l=1,num_nodos_comb
          ii=nodos_comb(l)
          read(9,*) lth, (rnxe(ii,k),k=2,number_zplane-1)
          do k=2,number_zplane-1
             rnxe(ii,k)=rnxe(ii,k)*(1e+24)
          end do
       end do
    end if
    close(9)
 endif  !!!!  abernal, ISIRYM, 18/10/2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
!
!!!!
!!!!      rmiro - UPV     / abernal - UPV                                !!!!
!    allocate(DED(0:(number_nodes_pz+1),0:number_zplane,3)) !!!! abernal
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    !tlap=dclock() !!!! abernal
!
    tfmax=0
!



! Determination of the axial position of the control rods
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 14/10/2013
    dim_zcr=maxval(control_rod_map)
    allocate(z_cr(dim_zcr))
    open(unit=8,file='input',status='OLD')
    do ii=1,15  ! The line 16 of file "input" contains the control rods axial position
       read(8,*)
    end do
    read(8,*,iostat=err) (z_cr(l),l=1,dim_zcr) ! Read the number of steps of each control rod bank
    if (err/=0) then
       print '(///)'
       print *,'ERROR READING THE LINE 16 OF FILE "input"'
       print '(///)'
       stop
    end if
    close(8)
    do ii=1,dim_zcr
       z_cr(ii)=fully_inserted_z+step_size*z_cr(ii) ! Transform the number of steps into axial length
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





! Calculation of the control rod volume fraction: crbdens(:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 30/10/2013
    ! First, the rector type is read from file "input"
    open(unit=8,file='input',status='OLD')
    do ii=1,17  ! The line 18 of file "input" contains the reactor type: PWR or BWR
       read(8,*)
    end do
    read(8,*,iostat=err) reactor_type ! Read the number of steps of each control rod bank
    if (err/=0) then
       print '(///)'
       print *,'ERROR READING THE LINE 16 OF FILE "input"'
       print '(///)'
       stop
    end if
    close(8)
    ! Loop through all the nodes to determine the control rod density
    allocate(crbdens(number_nodes_pz,number_zplane))
	do k=1,number_zplane
	   do l=1,number_nodes_pz
	      if ((control_rod_map(l)>0).and.(k>1).and.(k<number_zplane)) then
	         ii=control_rod_map(l)
	         select case(reactor_type)
	         case (0) ! PWR
	            if (z_cr(ii)<axial_mesh_nodes(k+1)) then
	               crbdens(l,k)=(axial_mesh_nodes(k+1)-z_cr(ii))/&
	               &(axial_mesh_nodes(k+1)-axial_mesh_nodes(k))
	               if (crbdens(l,k)>1) crbdens(l,k)=1.0
	            else
	               crbdens(l,k)=0
	            end if
             case (1) ! BWR
	            if (z_cr(ii)>axial_mesh_nodes(k)) then
	               crbdens(l,k)=(z_cr(ii)-axial_mesh_nodes(k))/&
	               &(axial_mesh_nodes(k+1)-axial_mesh_nodes(k))
	               if (crbdens(l,k)>1) crbdens(l,k)=1.0
	            else
	               crbdens(l,k)=0
	            end if
	         case default
                print '(///)'
                print *,'ERROR READING THE LINE 18 OF FILE "input". 0:PWR OR 1:BWR'
                print '(///)'
                stop
	         end select
	      else
	         crbdens(l,k)=0
	      end if
	   end do
	end do
	deallocate(z_cr) ! deallocate(z_cr,axial_mesh_nodes) !!!! abernal, ISIRYM, 22/10/2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      rmiro - UPV / Claubia - UFMG 16/03/2010  / abernal - UPV 10/10/2013
!!!!                                               !!!!
	if (.not.fdbk .and. rffdbk) then
       if (timeedt.lt.0.0) then
	      open(unit=969,file='DENS',status='OLD')
	      open(unit=970,file='TFUS',status='OLD')
          if (ibor.eq.1) then
		     open(unit=972,file='BORO',status='OLD')
	      end if

	      read(969,*)
	      read(970,*)
	      if (ibor.eq.1) then
	         read(972,*)
	         allocate (ppml(number_nodes_pz,number_zplane)) !!!! abernal
	      end if
	      allocate (dcool(number_zplane,number_nodes_pz),tdopl(number_zplane,number_nodes_pz)) !!!! abernal
	      do kth=1,number_zplane
	         read(969,*)
             read(970,*)
             read(969,*)(dcool(kth,lth),lth=1,number_nodes_pz)
             read(970,*)(tdopl(kth,lth),lth=1,number_nodes_pz)
	         if (ibor.eq.1) then
	            read(972,*)
	            read(972,*)(ppml(lth,kth),lth=1,number_nodes_pz)
             end if
		     do lth=1,number_nodes_pz
		        tdopl(kth,lth)=sqrt(tdopl(kth,lth))
	         end do
          end do
	   end if
       if (timeedt.lt.0.0) then
	      rewind(969)
	      rewind(970)
	      if (ibor.eq.1) then
	         rewind(972)
	      end if
       end if
       !if (time.ge.tend) then !!!! abernal
	      close(969)
	      close(970)
	      if (ibor.eq.1) then
	         close(972)
		  end if
	   !end if !!!! abernal --> if (time.ge.tend)
    end if

   if (.not.fdbk .and. .not.rffdbk) allocate(tcool(number_zplane,number_nodes_pz)) !!!! abernal

!

      allocate (rambda(6,ncmpur13)) !!!! abernal

      do k=1,number_zplane
         ipr=k !ipr=izpr(k)  !!!! abernal
         kth=k
         do l=1,number_nodes_pz
            ! DED(l,k,1) = fracdc  !!!! abernal
            ! DED(l,k,2) = fracdvb  !!!! abernal
            ! DED(l,k,3) = fracdwr  !!!! abernal
            la=l ! la=ltola(l) !!!! abernal
            ! lfa=ltolfa(l) !!!! abernal
            lth=l ! lth=ltochan(l) !!!! abernal
            icomp=materials_nodes(ipr,la) ! icomp=iprcomp(la,ipr) !!!! abernal
! rmmiller PBTT+
!           we assume there is little difference between the rodded and unrodded lambda's
!
            do ii=1,6
!             rambda(ii,icomp)=lambdapbur(ii,1)
               rambda(ii,icomp)=lambdapbur(ii,icomp)
            enddo
! rmmiller PBTT-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 28/10/2013
!             icr=icrcomp(icomp)
            if (crbdens(la,k)==0) then
               icr=0
            else
               icr=icomp
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(icr.eq.0) then
               rodfrac(1)=0.0
               rodfrac(2)=0.0
            else
               rodfrac(1)=crbdens(la,k)
               rodfrac(2)=crbdens(la,k)
            endif
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      rmiro - UPV                                                    !!!!
	        if (.not.fdbk .and. rffdbk) then
               tval = tdopl(kth,l)*tdopl(kth,l)
               dval = dcool(kth,l)
               tfmax = max(tfmax,tval)
	        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      rmiro - UPV                                                    !!!!
	        if (.not.fdbk .and. .not.rffdbk) then
               tval = Tfunif*Tfunif
	           dval = Dmunif
	           tfmax = max(tfmax,tval)
	           tdopl(kth,l) = Tfunif
	           dcool(kth,l) = Dmunif
	           tcool(kth,l) = Tmunif
	        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  abernal, ISIRYM, 18/10/2014
            if (fdbk) then !if (thcoupling==1) then !!!! abernal, 30/10/2014
               tval=temp_coupling((k-1)*number_nodes_pz+l)
               dval=dens_coupling((k-1)*number_nodes_pz+l)
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            call pbttxsec(k,l,kth,lth,icomp,icr,rodfrac,tval,dval)
            xsd(1,l,k)=1/(xstr(1,l,k)*3)
            xsd(2,l,k)=1/(xstr(2,l,k)*3)
            xst(1,l,k)=xsa(1,l,k)+xss(l,k)
            xst(2,l,k)=xsa(2,l,k)

         enddo
      enddo

!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      rmiro - UPV                                                    !!!!
      !if (time.gt.0 .and. timeedt.eq.time) then !!!! abernal
	     ! deallocate(DED) !!!! abernal
	     !return !!!! abernal
	  !end if !!!! abernal
!!!!      if(time.gt.0 .and. timeedt.eq.time) return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

!
      !txsec=txsec+dclock()-tlap !!!! abernal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      rmiro - UPV                                                    !!!!
      ! deallocate(DED) !!!! abernal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
 end subroutine xsecfb13

