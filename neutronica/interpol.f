!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Justin Watson
!      10-4-97
!     SUBROUTINE interpol(n, xval, xin, yin, f) 
!
!     This program is designed to read data from tables and Then call 
!     the subroutine "interp" which will interpolate the specified 
!     points from the table.
!
!     INPUT:
!     n                 =     Number of data points.
!     xval              =     The value for which the interpolation is done for.
!     xin               =     The strictly increasing or decreasing sequence of 
!                             knots.
!     yin               =     The prescribed function values at the knots.
!
!     OUTPUT:
!     f                 =     The interpolated value.
!
!     Called by:  table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE interpol(n, xval, xin, yin, f)
      IMPLICIT NONE
      integer i,n
      double precision xval,xin(n),yin(n),f
      DOUBLE PRECISION bot(100),cot(100),dot(100),eot(100),fot(100),q,p
!
!     Set up some intiial values.
!
      i = 0
!
!     quinat is the subroutine that calculates the coefficients between 
!     each set of data points.
!
      CALL quinat(n,xin,yin,bot,cot,dot,eot,fot)
!
!     Determine which coefficients should be used to do the interpolation 
!     and then use those coefficients to do the interpolation.
!
      DO 10 i=1,n
        IF ((xval .GE. xin(i)).AND.(xval .LE. xin(i+1))) THEN
          IF (i .EQ. n-1) THEN
            q = xin(i+1) - xval
            f = ((((-fot(i)*q + eot(i+1))*q - dot(i+1))*q 
     +          + cot(i+1))*q - bot(i+1))*q + yin(i+1)
          ELSE
            p = xval-xin(i)
            f = ((((fot(i)*p + eot(i))*p + dot(i))*p + cot(i))*p
     +          +bot(i))*p  + yin(i)
          END IF
         END IF
10    CONTINUE
      RETURN
      END
