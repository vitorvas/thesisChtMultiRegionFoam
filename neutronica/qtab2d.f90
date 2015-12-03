      SUBROUTINE qtab2d(n1,n2,tx1,tx2,ty,x1,x2,f)                                    !
      IMPLICIT none
      INTEGER n1,n2, i,j 
      double precision tx1(n1),tx2(n2),ty(n1,n2),x1,x2,f,y(n2)
      DO j=1, n2
         CALL interpol(n1,x1,tx1,ty(1,j),y(j))
      ENDDO
      CALL interpol(n2,x2,tx2,y,f)
      RETURN
      END
