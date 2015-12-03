 subroutine pbttxsec(k,l,kth,lth,icomp,icr,rodfrac,tval,dval)
! This subroutine incorporates all PBTT feedback mechanisms
! into x-sections
!
!
!   no xenon absorption in group 1 for the PBTT benchmark
!   no samarium computations in PBTT benchmark
!   xsxea(1,l,k)=0.   |
!   xssma(1,l,k)=0.   |> already set to zero in initxsec.F
!   xssma(2,l,k)=0.   |
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                               abernal, ISIRYM, 21/10/2013
      use variables_xs_module

      integer k, kth, lth, l, icomp, icr
      double precision tval, dval, rodfrac(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!
!
! PBTT header
!
! added for test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                               abernal, ISIRYM, 21/10/2013
!      logical premier
!      save premier
!      data premier/.true./
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                               abernal, ISIRYM, 21/10/2013
!      dimension rodfrac(ng)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER it(4)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  rmiro - UPV        /     abernal - ISIRYM  (21/10/2013)           !!!!
!!!!  Double Precision tb(48), x(4), kconv(2), fxs, xs(6,ng), xsr(6,ng) !Tony's
!      Double Precision tb(48), x(4), kconv(2), fxs, xs(6,ng), xsr(6,ng) !rmiro
!     &,xsb(6,ng), xsrb(6,ng) !rmiro
      Double Precision x(4), kconv(2), fxs, xs(6,2), xsr(6,2) & !!!! abernal
     &,xsb(6,2), xsrb(6,2)  !!!! abernal
      double precision, allocatable, dimension(:) :: tb
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Double Precision adf(4,ng) , df(ng), dxs(ng)                      !Tony's
       Double Precision adf(4,2) , df(2), dxs(2) !!!! abernal
     ! data eps/0.001/ !!!! abernal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Claubia - UFMG rmiro - UPV                                       !!!!
!!!!
      double precision xsaux

!
!!!!                                                                    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                               abernal, ISIRYM, 21/10/2013
      integer irot, m, mn, ng
      double precision unrodfrac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                               abernal, ISIRYM, 21/10/2013
      ng=2 ! number of energy groups
      allocate(tb(ntdentry))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
      irot=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                               abernal, ISIRYM, 21/10/2013
!      if(rotadf) irot=lrotfa(l)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
      kconv(1)=.3213e-10
      kconv(2)=.3206e-10
!

      if (icr.eq.0 .or. rodfrac(1).ne.1.0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Claubia - UFMG / rmiro - UPV                                     !!!!
	     if (ibor.eq.0) then
            call xsipbtt(kth,lth,tval,dval,xspbttur(1,1,1,icomp) &
     &,               adfpbur(1,1,1,icomp),detpbur(1,1,1,icomp),xs &
     &,               adf,df,dxs)                                ! jcr
         end if
!!!!
         if (ibor.eq.1)then
	        if (ppml(l,k).lt.ppml1 .or. ppml(l,k).gt.ppml3) then
	           write(*,*)'Boron Concentration out of range',l,k,ppml(l,k)
	           STOP
		    end if
            if (ppml(l,k).ge.ppml1 .and. ppml(l,k).le.ppml2) then
               call xsipbtt(kth,lth,tval,dval,xspbttur(1,1,1,icomp) &
     &,             adfpbur(1,1,1,icomp),detpbur(1,1,1,icomp),xs &
     &,             adf,df,dxs)                                  ! jcr
	           call xsipbtt(kth,lth,tval,dval,xspbttur2(1,1,1,icomp) &
     &,             adfpbur(1,1,1,icomp),detpbur(1,1,1,icomp),xsb &
     &,             adf,df,dxs)                                  ! jcr
               do m=1,ng
                  do mn=1,6
                     xsaux=0.0
                     xsaux= (((ppml2-ppml(l,k))*xs(mn,m))/(ppml2-ppml1))+ &
     &               (((ppml(l,k)-ppml1)*xsb(mn,m))/(ppml2-ppml1))
                     xs(mn,m)=xsaux
                  end do
               end do
            end if
!!!!
            if (ppml(l,k).gt.ppml2 .and. ppml(l,k).le.ppml3) then
	           call xsipbtt(kth,lth,tval,dval,xspbttur2(1,1,1,icomp) &
     &,             adfpbur(1,1,1,icomp),detpbur(1,1,1,icomp),xs &
     &,             adf,df,dxs)                                  ! jcr
               call xsipbtt(kth,lth,tval,dval,xspbttur3(1,1,1,icomp) &
     &,             adfpbur(1,1,1,icomp),detpbur(1,1,1,icomp),xsb &
     &,             adf,df,dxs)                                  ! jcr
               do m=1,ng
                  do mn=1,6
                     xsaux=0.0
                     xsaux= (((ppml3-ppml(l,k))*xs(mn,m))/(ppml3-ppml2))+ &
     &               (((ppml(l,k)-ppml2)*xsb(mn,m))/(ppml3-ppml2))
                     xs(mn,m)=xsaux
                  end do
               end do
            end if
         endif
!!!!                                                                    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do m=1,ng
            unrodfrac=1-rodfrac(m)
            xstr(m,l,k) = (1/(3*xs(1,m)))*unrodfrac
            xsa(m,l,k)  = xs(2,m)*unrodfrac
            xsnf(m,l,k) = xs(4,m)*unrodfrac
            xskf(m,l,k) = kconv(m)*xs(3,m)*unrodfrac
            xsf(m,l,k)=  xs(3,m)*unrodfrac

!if rotadf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                               abernal, ISIRYM, 21/10/2013
!            if(rotadf) then
!               adfn(m,l,k) = adf(mod(irot  ,4)+1,m)*unrodfrac
!               adfe(m,l,k) = adf(mod(irot+1,4)+1,m)*unrodfrac
!               adfs(m,l,k) = adf(mod(irot+2,4)+1,m)*unrodfrac
!               adfw(m,l,k) = adf(mod(irot+3,4)+1,m)*unrodfrac
!            else
!              adfn(m,l,k) = 1.
!              adfe(m,l,k) = 1.
!              adfw(m,l,k) = 1.
!              adfs(m,l,k) = 1.
!            endif
!            IF (detector) THEN
!               detflx(m,l,k) =  df(m)*(1-rodfrac(2))
!               detxs(m,l,k)  = dxs(m)*df(m)*(1-rodfrac(2))  !v25r1m0
!            END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         enddo
         xss(l,k) = xs(5,1)*(1-rodfrac(1))
!
! HE COMENTADO LA PARTE DEL XENON. SI SE QUIERE INCORPORAR, SE NECESITA LA
! VARIABLE "rnxe(l,k)", QUE ES "XENON NUMBER DENSITY"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                               abernal, ISIRYM, 21/10/2013
!#ifdef PBTT_XE_ON
!         xsxea(2,l,k)=xs(6,2) * (1-rodfrac(2) )
!!        Sigma_a = Sigma_a - Sigma_a,XE + N_Xe * sigma_a,Xe , jcr  !ESTA LINEA TIENE QUE ESTAR COMENTADA
!         xsa(2,l,k)  = xsa(2,l,k) - xs(5,2)*(1-rodfrac(2)) +
!     &                       rnxe(l,k)*xsxea(2,l,k)
!#endif
         xsxea(2,l,k)=xs(6,2) * (1-rodfrac(2) )
         !Sigma_a = Sigma_a - Sigma_a,XE + N_Xe * sigma_a,Xe , jcr  !ESTA LINEA TIENE QUE ESTAR COMENTADA
         xsa(2,l,k)  = xsa(2,l,k) - xs(5,2)*(1-rodfrac(2)) + &
     &                       rnxe(l,k)*xsxea(2,l,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
         if(rodfrac(2).eq.0.0) return
!
      else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Claubia - UFMG rmiro - UPV                                       !!!!
	     if (ibor.eq.0) then
            call xsipbtt(kth,lth,tval,dval,xspbttr(1,1,1,icr) &
     &,                 adfpbr(1,1,1,icr),detpbr(1,1,1,icr),xs  &
     &,                 adf,df,dxs)                              ! jcr
         end if
!!!!
         if (ibor.eq.1)then
	        if (ppml(l,k).lt.ppml1 .or. ppml(l,k).gt.ppml3) then
	           write(*,*)'Boron Concentration out of range',l,k,ppml(l,k)
	           STOP
		    end if
            if (ppml(l,k).ge.ppml1 .and. ppml(l,k).le.ppml2) then
               call xsipbtt(kth,lth,tval,dval,xspbttr(1,1,1,icr) &
     &,             adfpbr(1,1,1,icr),detpbr(1,1,1,icr),xs  &
     &,             adf,df,dxs)                                  ! jcr
	           call xsipbtt(kth,lth,tval,dval,xspbttr2(1,1,1,icr) &
     &,             adfpbr(1,1,1,icr),detpbr(1,1,1,icr),xsb &
     &,             adf,df,dxs)                                  ! jcr
               do m=1,ng
                  do mn=1,6
                     xsaux=0.0
                     xsaux= (((ppml2-ppml(l,k))*xs(mn,m))/(ppml2-ppml1))+ &
     &               (((ppml(l,k)-ppml1)*xsb(mn,m))/(ppml2-ppml1))
                     xs(mn,m)=xsaux
                  end do
               end do
            end if
!!!!
            if (ppml(l,k).gt.ppml2 .and. ppml(l,k).le.ppml3) then
	           call xsipbtt(kth,lth,tval,dval,xspbttr2(1,1,1,icr) &
     &,             adfpbr(1,1,1,icr),detpbr(1,1,1,icr),xs &
     &,             adf,df,dxs)                                  ! jcr
               call xsipbtt(kth,lth,tval,dval,xspbttr3(1,1,1,icr) &
     &,             adfpbr(1,1,1,icr),detpbr(1,1,1,icr),xsb &
     &,             adf,df,dxs)                                  ! jcr
               do m=1,ng
                  do mn=1,6
                     xsaux=0.0
                     xsaux= (((ppml3-ppml(l,k))*xs(mn,m))/(ppml3-ppml2))+ &
     &               (((ppml(l,k)-ppml2)*xsb(mn,m))/(ppml3-ppml2))
                     xs(mn,m)=xsaux
                  end do
               end do
            end if
         endif
!!!!                                                                    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
         do m=1,ng
            xstr(m,l,k) = 1/(3*xs(1,m))
            xsa(m,l,k)  = xs(2,m)
            xsnf(m,l,k) = xs(4,m)
            xskf(m,l,k) = kconv(m)*xs(3,m)
            xsf(m,l,k)=  xs(3,m)
!if rotadf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                               abernal, ISIRYM, 21/10/2013
!            if(rotadf) then
!               adfn(m,l,k) = adf(mod(irot  ,4)+1,m)
!               adfe(m,l,k) = adf(mod(irot+1,4)+1,m)
!               adfs(m,l,k) = adf(mod(irot+2,4)+1,m)
!               adfw(m,l,k) = adf(mod(irot+3,4)+1,m)
!            else
!              adfn(m,l,k) = 1.
!              adfe(m,l,k) = 1.
!              adfw(m,l,k) = 1.
!              adfs(m,l,k) = 1.
!            endif
!            if((l.eq.42).and.(k.eq.4)) then
!               continue
!            endif
!            IF (detector) THEN
!               detflx(m,l,k) = df(m)
!               detxs(m,l,k)  = dxs(m)*df(m)   !v25r1m0
!            END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         enddo
         xss(l,k) = xs(5,1)
!

! HE COMENTADO LA PARTE DEL XENON. SI SE QUIERE INCORPORAR, SE NECESITA LA
! VARIABLE "rnxe(l,k)", QUE ES "XENON NUMBER DENSITY"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                               abernal, ISIRYM, 21/10/2013
!#ifdef PBTT_XE_ON
!         xsxea(2,l,k)=xs(6,2)
!!        Sigma_a = Sigma_a - Sigma_a,XE + N_Xe * sigma_a,Xe , jcr  !ESTA LINEA TIENE QUE ESTAR COMENTADA
!         xsa(2,l,k)  = xsa(2,l,k) - xs(5,2) + rnxe(l,k)*xsxea(2,l,k)
!#endif
         xsxea(2,l,k)=xs(6,2)
         !Sigma_a = Sigma_a - Sigma_a,XE + N_Xe * sigma_a,Xe , jcr  !ESTA LINEA TIENE QUE ESTAR COMENTADA
         xsa(2,l,k)  = xsa(2,l,k) - xs(5,2) + rnxe(l,k)*xsxea(2,l,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
         return
      endif
!
! this part adds the unrodded and rodded xsections together for PRN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Claubia - UFMG rmiro - UPV                                       !!!!
	  if (ibor.eq.0) then
         call xsipbtt(kth,lth,tval,dval,xspbttr(1,1,1,icr) &
     &,              adfpbr(1,1,1,icr),detpbr(1,1,1,icr),xsr &
     &,              adf,df,dxs)                                 ! jcr
      end if
!!!!
      if (ibor.eq.1)then
	     if (ppml(l,k).lt.ppml1 .or. ppml(l,k).gt.ppml3) then
	        write(*,*)'Boron Concentration out of range',l,k,ppml(l,k)
	        STOP
	     end if
         if (ppml(l,k).ge.ppml1 .and. ppml(l,k).le.ppml2) then
            call xsipbtt(kth,lth,tval,dval,xspbttr(1,1,1,icr) &
     &,          adfpbr(1,1,1,icr),detpbr(1,1,1,icr),xsr &
     &,          adf,df,dxs)                                     ! jcr
	        call xsipbtt(kth,lth,tval,dval,xspbttr2(1,1,1,icr) &
     &,          adfpbr(1,1,1,icr),detpbr(1,1,1,icr),xsrb &
     &,          adf,df,dxs)                                     ! jcr
            do m=1,ng
               do mn=1,6
                  xsaux=0.0
                  xsaux= (((ppml2-ppml(l,k))*xsr(mn,m))/(ppml2-ppml1))+ &
     &            (((ppml(l,k)-ppml1)*xsrb(mn,m))/(ppml2-ppml1))
                  xsr(mn,m)=xsaux
               end do
            end do
         end if
!!!!
         if (ppml(l,k).gt.ppml2 .and. ppml(l,k).le.ppml3) then
	        call xsipbtt(kth,lth,tval,dval,xspbttr2(1,1,1,icr) &
     &,          adfpbr(1,1,1,icr),detpbr(1,1,1,icr),xsr &
     &,          adf,df,dxs)                                     ! jcr
            call xsipbtt(kth,lth,tval,dval,xspbttr3(1,1,1,icr) &
     &,          adfpbr(1,1,1,icr),detpbr(1,1,1,icr),xsrb &
     &,          adf,df,dxs)                                     ! jcr
            do m=1,ng
               do mn=1,6
                  xsaux=0.0
                  xsaux= (((ppml3-ppml(l,k))*xsr(mn,m))/(ppml3-ppml2))+ &
     &            (((ppml(l,k)-ppml2)*xsrb(mn,m))/(ppml3-ppml2))
                  xsr(mn,m)=xsaux
               end do
            end do
         end if
      endif
!!!!                                                                    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do m=1,ng
         xstr(m,l,k) = xstr(m,l,k) + (1/(3*xsr(1,m)))*rodfrac(m)
         xsa(m,l,k)  = xsa(m,l,k)  + xsr(2,m)*rodfrac(m)
         xsnf(m,l,k) = xsnf(m,l,k) + xsr(4,m)*rodfrac(m)
         xskf(m,l,k) = xskf(m,l,k) + kconv(m)*xsr(3,m)*rodfrac(m)
         xsf(m,l,k)=  xsf(m,l,k) + xsr(3,m)*rodfrac(m)

!if rotadf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                               abernal, ISIRYM, 21/10/2013
!         if(rotadf) then
!            adfn(m,l,k) = adfn(m,l,k) + adf(mod(irot  ,4)+1,m)*         &
!     &                    rodfrac(m)
!            adfe(m,l,k) = adfe(m,l,k) + adf(mod(irot+1,4)+1,m)*         &
!     &                    rodfrac(m)
!            adfs(m,l,k) = adfs(m,l,k) + adf(mod(irot+2,4)+1,m)*         &
!     &                    rodfrac(m)
!            adfw(m,l,k) = adfw(m,l,k) + adf(mod(irot+3,4)+1,m)*         &
!     &                    rodfrac(m)
!         else
!            adfn(m,l,k) = 1.
!            adfe(m,l,k) = 1.
!            adfw(m,l,k) = 1.
!            adfs(m,l,k) = 1.
!         endif
!         if((l.eq.42).and.(k.eq.4)) then
!            continue
!         endif
!         IF (detector) THEN
!            detflx(m,l,k) = detflx(m,l,k) +  df(m)*rodfrac(2)
!            detxs(m,l,k)  = detxs(m,l,k)  + dxs(m)*df(m)*rodfrac(2)  !v25r1m0
!         END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
      xss(l,k) = xss(l,k)+ xsr(5,1)*rodfrac(1)
!

! HE COMENTADO LA PARTE DEL XENON. SI SE QUIERE INCORPORAR, SE NECESITA LA
! VARIABLE "rnxe(l,k)", QUE ES "XENON NUMBER DENSITY"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                               abernal, ISIRYM, 21/10/2013
!#ifdef PBTT_XE_ON
!      xsxea(2,l,k)=xs(6,2)*rodfrac(2)
!!     Sigma_a = Sigma_a - Sigma_a,XE + N_Xe * sigma_a,Xe , jcr  !ESTA LINEA TIENE QUE ESTAR COMENTADA
!      xsa(2,l,k)  = xsa(2,l,k) - xs(5,2)*rodfrac(2) +
!     &              rnxe(l,k)*xsxea(2,l,k)
!#endif
      xsxea(2,l,k)=xs(6,2)*rodfrac(2)
      !Sigma_a = Sigma_a - Sigma_a,XE + N_Xe * sigma_a,Xe , jcr  !ESTA LINEA TIENE QUE ESTAR COMENTADA
      xsa(2,l,k)  = xsa(2,l,k) - xs(5,2)*rodfrac(2) + &
     &              rnxe(l,k)*xsxea(2,l,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                abernal, ISIRYM, 23/10/2013
      deallocate(tb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      return
 end subroutine pbttxsec

