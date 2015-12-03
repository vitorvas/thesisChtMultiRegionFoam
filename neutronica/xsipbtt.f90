 subroutine xsipbtt(kth,lth,tf,dm,xstab,adftab,dettab,xsout,adfout &
     &,                  df,dxs)
!
! find linearly interploated cross section for the given fuel temperature (tf)
! and moderator density condition (dm)
!
!  - initiate mask condition search starting from the mask index stored during
!    the in the previous call to this routine for the t/h node specified by
!    kth (axial level index) and lth (radial channel number)
!
! I/O's
!  - xstab:  tabularized cross section for each composition
!  - adftab: tabularized adf's for each composition
!  - dettab: tabularized detector reponses for each composition
!  - xsout:  two-group cross sections for 5 reaction types
!            stored common
!  - adf
!  - det
!  - itabindx(idim,kth,lth): the starting mask index
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 22/10/2013
      use variables_xs_module

      integer kth,lth
      double precision tf,dm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      logical, parameter :: linearint=.true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      rmiro - UPV     /   abernal - ISIRYM                           !!!!
!!!!
!      integer ntfentry,ndmentry,ntentry,ibor !!!! abernal 21/10/2013
!!!!   Claubia - UFMG rmiro - UPV      04/03/2010                        !!!!
!      double precision ppml1,ppml2,ppml3 !!!! abernal 21/10/2013
!	common/simtab/ntfentry,ndmentry,ntentry,ibor,ppml1,ppml2,ppml3 !!!! abernal 21/10/2013
!!!!
!      parameter (ntabdim=2) !table dimension !!!! abernal 22/10/2013
!!!!  parameter (ntabdim=2,ntfentry=6,ndmentry=6) !table dimension
!!!!  parameter (ntentry=ntfentry+ndmentry)
!!!!  parameter (ntyb=ntentry+1,ntx2b=ntfentry+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 22/10/2013
!      dimension tabval(ntentry),xvals(ntabdim)
!      dimension indxl(ntabdim),indxr(ntabdim),nvalbase(ntabdim)
!      dimension xswfl(ntabdim),xswfr(ntabdim),nentries(ntabdim)
!!      real xstab, xsout, adftab, dettab, adfout, df, dxs ! jcr !orig
!      Double Precision xstab, xsout, adftab, dettab, adfout, df, dxs ! jcr !Tony's
!      dimension xstab(48,6,ng)  ! input xsec set
!      dimension adftab(48,2,ng) ! input adf set
!      dimension dettab(48,2,ng) ! input det set
!      dimension xsout(6,ng)     ! output xsec
!      dimension adfout(4,ng)    ! output xsec
!      dimension df(ng)          ! output xsec
!      dimension dxs(ng)         ! output xsec
      double precision tabval(ntentry)
      double precision xswfl(2), xswfr(2) ! 2:ntabdim --> table dimension
      double precision xstab(ntdentry,6,2), adftab(ntdentry,2,2), dettab(ntdentry,2,2) ! (:,:,2) --> 2=ng --> energy groups
      double precision xsout(6,2), adfout(4,2)  ! (:,2) --> 2=ng --> energy groups
      double precision df(2), dxs(2) ! 2=ng --> energy groups
      double precision xvals(2) ! 2:ntabdim --> table dimension
      integer indxl(2), indxr(2), nvalbase(2), nentries(2), itabindxpb(2) ! 2:ntabdim --> table dimension
      double precision xsb, xst_2

      integer ng, ntabdim, ntyb, ntx2b, idim, ileft, iright, nbase, nentry, aux_while, &
      & locb, loct, locbl, locbr, loctl, loctr, i, ir, m
      double precision xval, big

      ntabdim=2
      ng=2
      big=1.0e30
      itabindxpb(:)=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      rmiro - UPV                                                    !!!!
!!!!  data nentries/ntfentry,ndmentry/
!!!!  data nvalbase/0,ntfentry/
      nentries(1)=ntfentry
      nentries(2)=ndmentry
      nvalbase(1)=0
      nvalbase(2)=ntfentry
      ntyb=ntentry+1
      ntx2b=ntfentry+1
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!prevent extrapolation                                     !tk 2.1m3
      if (tf.gt.xstab(ntfentry,1,1)) then
         tf=xstab(ntfentry,1,1)
      elseif (tf.lt.xstab(1,1,1)) then
         tf=xstab(1,1,1)
      endif
      if (dm.gt.xstab(ntfentry+ndmentry,1,1)) then
         dm=xstab(ntfentry+ndmentry,1,1)
      elseif (dm.lt.xstab(ntfentry+1,1,1)) then
         dm=xstab(ntfentry+1,1,1)
      endif
      if(linearint)then
!
! find the appropriate positions in the variable mask
         xvals(1)=tf
         xvals(2)=dm
         do idim=1,ntabdim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 22/10/201
!            ileft=itabindxpb(idim,kth,lth)
            ileft=itabindxpb(idim)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            iright=ileft+1
            nbase=nvalbase(idim)
            xval=xvals(idim)
            nentry=nentries(idim)
            do i=1,nentry
               tabval(i)=xstab(i+nbase,1,1)
            enddo
            tabval(nentry+1)=big
!
            if(xval.lt.tabval(ileft)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 22/10/2013
!               do i=ileft-1,1,-1  !decrease index
               aux_while=0
               i=ileft
               do while (aux_while==0)  !decrease index
                  i=i-1
                  if (i==1) aux_while=2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  if(xval.ge.tabval(i)) then
                     ileft=i
                     iright=ileft+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 22/10/2013
!                     go to 10
                     aux_while=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  endif
               enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 22/10/2013
!!                           xval lower than allowed minimum  !!!! IT WAS COMMENTED
!               ileft=1
!!               iright=ileft  !!!! IT WAS COMMENTED
!               iright=ileft+1  !allow extrapolation
!            endif
!            if(xval.gt.tabval(iright)) then
!               do i=iright+1,nentry  !increase index
!                  if(xval.le.tabval(i)) then
!                     ileft=i-1
!                     iright=ileft+1
!                     go to 10
!                  endif
!               enddo
!!                           xval greater than allowed maximum  !!!! IT WAS COMMENTED
!!               ileft=nentry  !!!! IT WAS COMMENTED
!!               iright=ileft  !!!! IT WAS COMMENTED
!               ileft=nentry-1 !allow extrapolation
!               iright=ileft+1
!            endif
!   10       continue
               if (aux_while==2) then
                  ileft=1
                  iright=ileft+1  !allow extrapolation
               endif
            else
               aux_while=0
            endif
            if (aux_while/=1) then
               if(xval.gt.tabval(iright)) then
                  aux_while=0
                  i=iright
                  do while (aux_while==0)
                     i=i+1
                     if (i==nentry) aux_while=2
                     if(xval.le.tabval(i)) then
                        ileft=i-1
                        iright=ileft+1
                        aux_while=1
                     endif
                  enddo
                  if (aux_while==2) then
                     ileft=nentry-1 !allow extrapolation
                     iright=ileft+1
                  endif
               endif
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            indxl(idim)=ileft
            indxr(idim)=iright
!
            if(ileft.ne.iright) then
               xswfl(idim)=(tabval(iright)-xval)                           &
     &                 /(tabval(iright)-tabval(ileft))
            else
               xswfl(idim)=1
            endif
            xswfr(idim)=1-xswfl(idim)
!
! update table index
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 22/10/201
!            itabindxpb(idim,kth,lth)=ileft
            itabindxpb(idim)=ileft
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         enddo
!
! determine the location at the table
         locb=(indxl(2)-1)*ntfentry+ntentry
         loct=(indxr(2)-1)*ntfentry+ntentry
         locbl=locb+indxl(1)
         locbr=locb+indxr(1)
         loctl=loct+indxl(1)
         loctr=loct+indxr(1)
!
         do m=1,ng
            do ir=1,5
               xsb=xswfl(1)*xstab(locbl,ir,m)+xswfr(1)*xstab(locbr,ir,m)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 22/10/201
!               xst=xswfl(1)*xstab(loctl,ir,m)+xswfr(1)*xstab(loctr,ir,m)
!               xsout(ir,m)=xsb*xswfl(2)+xst*xswfr(2)
               xst_2=xswfl(1)*xstab(loctl,ir,m)+xswfr(1)*xstab(loctr,ir,m)
               xsout(ir,m)=xsb*xswfl(2)+xst_2*xswfr(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            enddo
         enddo
!
         xsb=xswfl(1)*xstab(locbl,6,2)+xswfr(1)*xstab(locbr,6,2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 22/10/201
!         xst=xswfl(1)*xstab(loctl,6,2)+xswfr(1)*xstab(loctr,6,2)
!         xsout(6,2)=xsb*xswfl(2)+xst*xswfr(2)
         xst_2=xswfl(1)*xstab(loctl,6,2)+xswfr(1)*xstab(loctr,6,2)
         xsout(6,2)=xsb*xswfl(2)+xst_2*xswfr(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! adf's, jcr
! pbtt input:1=w, 2=s // Parcs notation: 1=n, 2=e, 3=s, 4=w
! in pbtt, n=w and e=s
         do m=1,ng
            do ir=1,2
               xsb=xswfl(1)*adftab(locbl,ir,m)+xswfr(1)*adftab(locbr,ir,m)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 22/10/201
!               xst=xswfl(1)*adftab(loctl,ir,m)+xswfr(1)*adftab(loctr,ir,m)
!               adfout(ir,m)=xsb*xswfl(2)+xst*xswfr(2)
               xst_2=xswfl(1)*adftab(loctl,ir,m)+xswfr(1)*adftab(loctr,ir,m)
               adfout(ir,m)=xsb*xswfl(2)+xst_2*xswfr(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            enddo
	        adfout(3,m)=adfout(2,m)
            adfout(4,m)=adfout(1,m)
         enddo
!
         ir=1
	     do m=1,ng
	        xsb=xswfl(1)*dettab(locbl,ir,m)+xswfr(1)*dettab(locbr,ir,m)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 22/10/201
!            xst=xswfl(1)*dettab(loctl,ir,m)+xswfr(1)*dettab(loctr,ir,m)
!            df(m)=xsb*xswfl(2)+xst*xswfr(2)
            xst_2=xswfl(1)*dettab(loctl,ir,m)+xswfr(1)*dettab(loctr,ir,m)
            df(m)=xsb*xswfl(2)+xst_2*xswfr(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         enddo
!
! detector xs, jcr
!
         ir=2
	     do m=1,ng
	        xsb=xswfl(1)*dettab(locbl,ir,m)+xswfr(1)*dettab(locbr,ir,m)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                       abernal, 22/10/201
!            xst=xswfl(1)*dettab(loctl,ir,m)+xswfr(1)*dettab(loctr,ir,m)
!            dxs(m)=xsb*xswfl(2)+xst*xswfr(2)
            xst_2=xswfl(1)*dettab(loctl,ir,m)+xswfr(1)*dettab(loctr,ir,m)
            dxs(m)=xsb*xswfl(2)+xst_2*xswfr(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         enddo
      else
        
         do m=1,ng
            do ir=1,5
               call qtab2d(ntfentry,ndmentry,xstab(1,ir,m), &
     &              xstab(ntx2b,ir,m),xstab(ntyb,ir,m), &
     &              tf,dm,xsout(ir,m))
            enddo
         enddo
!
         ir=6
         m=2
         call qtab2d(ntfentry,ndmentry,xstab(1,ir,m), &
     &              xstab(ntx2b,ir,m),xstab(ntyb,ir,m), &
     &              tf,dm,xsout(ir,m))
        
!
! adf's, jcr
! pbtt input:1=w, 2=s // Parcs notation: 1=n, 2=e, 3=s, 4=w
! in pbtt, n=w and e=s
         do m=1,ng
            do ir=1,2
               call qtab2d(ntfentry,ndmentry,adftab(1,ir,m), &
     &              adftab(ntx2b,ir,m),adftab(ntyb,ir,m), &
     &              tf,dm,adfout(ir,m))
            enddo
	        adfout(3,m)=adfout(2,m)
            adfout(4,m)=adfout(1,m)
         enddo
!
         ir=1
	     do m=1,ng
            call qtab2d(ntfentry,ndmentry,dettab(1,ir,m), &
     &              dettab(ntx2b,ir,m),dettab(ntyb,ir,m), &
     &              tf,dm,df(m))
         enddo
!
! detector xs, jcr
!
         ir=2
	     do m=1,ng
            call qtab2d(ntfentry,ndmentry,dettab(1,ir,m), &
     &              dettab(ntx2b,ir,m),dettab(ntyb,ir,m), &
     &              tf,dm,dxs(m))
         enddo
      endif
	  return
 end subroutine xsipbtt
