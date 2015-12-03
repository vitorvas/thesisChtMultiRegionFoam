C     ALGORITHM 600, COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.9, NO. 2,
C     JUN., 1983, P. 258-259.
      SUBROUTINE QUINAT(N, X, Y, B, C, D, E, F)                         QUI   10
C
C	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	IMPLICIT NONE
      INTEGER N
C      REAL X(N), Y(N), B(N), C(N), D(N), E(N), F(N)
      DOUBLE PRECISION X(N), Y(N), B(N), C(N), D(N), E(N), F(N)
C
C
C
C     QUINAT COMPUTES THE COEFFICIENTS OF A QUINTIC NATURAL QUINTIC SPLI
C     S(X) WITH KNOTS X(I) INTERPOLATING THERE TO GIVEN FUNCTION VALUES:
C               S(X(I)) = Y(I)  FOR I = 1,2, ..., N.
C     IN EACH INTERVAL (X(I),X(I+1)) THE SPLINE FUNCTION S(XX) IS A
C     POLYNOMIAL OF FIFTH DEGREE:
C     S(XX) = ((((F(I)*P+E(I))*P+D(I))*P+C(I))*P+B(I))*P+Y(I)    (*)
C           = ((((-F(I)*Q+E(I+1))*Q-D(I+1))*Q+C(I+1))*Q-B(I+1))*Q+Y(I+1)
C     WHERE  P = XX - X(I)  AND  Q = X(I+1) - XX.
C     (NOTE THE FIRST SUBSCRIPT IN THE SECOND EXPRESSION.)
C     THE DIFFERENT POLYNOMIALS ARE PIECED TOGETHER SO THAT S(X) AND
C     ITS DERIVATIVES UP TO S"" ARE CONTINUOUS.
C
C        INPUT:
C
C     N          NUMBER OF DATA POINTS, (AT LEAST THREE, I.E. N > 2)
C     X(1:N)     THE STRICTLY INCREASING OR DECREASING SEQUENCE OF
C                KNOTS.  THE SPACING MUST BE SUCH THAT THE FIFTH POWER
C                OF X(I+1) - X(I) CAN BE FORMED WITHOUT OVERFLOW OR
C                UNDERFLOW OF EXPONENTS.
C     Y(1:N)     THE PRESCRIBED FUNCTION VALUES AT THE KNOTS.
C
C        OUTPUT:
C
C     B,C,D,E,F  THE COMPUTED SPLINE COEFFICIENTS AS IN (*).
C         (1:N)  SPECIFICALLY
C                B(I) = S'(X(I)), C(I) = S"(X(I))/2, D(I) = S"'(X(I))/6,
C                E(I) = S""(X(I))/24,  F(I) = S""'(X(I))/120.
C                F(N) IS NEITHER USED NOR ALTERED.  THE FIVE ARRAYS
C                B,C,D,E,F MUST ALWAYS BE DISTINCT.
C
C        OPTION:
C
C     IT IS POSSIBLE TO SPECIFY VALUES FOR THE FIRST AND SECOND
C     DERIVATIVES OF THE SPLINE FUNCTION AT ARBITRARILY MANY KNOTS.
C     THIS IS DONE BY RELAXING THE REQUIREMENT THAT THE SEQUENCE OF
C     KNOTS BE STRICTLY INCREASING OR DECREASING.  SPECIFICALLY:
C
C     IF X(J) = X(J+1) THEN S(X(J)) = Y(J) AND S'(X(J)) = Y(J+1),
C     IF X(J) = X(J+1) = X(J+2) THEN IN ADDITION S"(X(J)) = Y(J+2).
C
C     NOTE THAT S""(X) IS DISCONTINUOUS AT A DOUBLE KNOT AND, IN
C     ADDITION, S"'(X) IS DISCONTINUOUS AT A TRIPLE KNOT.  THE
C     SUBROUTINE ASSIGNS Y(I) TO Y(I+1) IN THESE CASES AND ALSO TO
C     Y(I+2) AT A TRIPLE KNOT.  THE REPRESENTATION (*) REMAINS
C     VALID IN EACH OPEN INTERVAL (X(I),X(I+1)).  AT A DOUBLE KNOT,
C     X(J) = X(J+1), THE OUTPUT COEFFICIENTS HAVE THE FOLLOWING VALUES:
C       Y(J) = S(X(J))          = Y(J+1)
C       B(J) = S'(X(J))         = B(J+1)
C       C(J) = S"(X(J))/2       = C(J+1)
C       D(J) = S"'(X(J))/6      = D(J+1)
C       E(J) = S""(X(J)-0)/24     E(J+1) = S""(X(J)+0)/24
C       F(J) = S""'(X(J)-0)/120   F(J+1) = S""'(X(J)+0)/120
C     AT A TRIPLE KNOT, X(J) = X(J+1) = X(J+2), THE OUTPUT
C     COEFFICIENTS HAVE THE FOLLOWING VALUES:
C       Y(J) = S(X(J))         = Y(J+1)    = Y(J+2)
C       B(J) = S'(X(J))        = B(J+1)    = B(J+2)
C       C(J) = S"(X(J))/2      = C(J+1)    = C(J+2)
C       D(J) = S"'((X(J)-0)/6    D(J+1) = 0  D(J+2) = S"'(X(J)+0)/6
C       E(J) = S""(X(J)-0)/24    E(J+1) = 0  E(J+2) = S""(X(J)+0)/24
C       F(J) = S""'(X(J)-0)/120  F(J+1) = 0  F(J+2) = S""'(X(J)+0)/120
C
      INTEGER I, M
C      REAL B1, P, PQ, PQQR, PR, P2, P3, Q, QR, Q2, Q3, R, 
C     +                 R2, S, T, U, V
      DOUBLE PRECISION B1, P, PQ, PQQR, PR, P2, P3, Q, QR, Q2, Q3, R, 
     +                 R2, S, T, U, V
C
      IF (N.LE.2) GO TO 190
C
C     COEFFICIENTS OF A POSITIVE DEFINITE, PENTADIAGONAL MATRIX,
C     STORED IN D,E,F FROM 2 TO N-2.
C
      M = N - 2
      Q = X(2) - X(1)
      R = X(3) - X(2)
      Q2 = Q*Q
      R2 = R*R
      QR = Q + R
      D(1) = 0.
      E(1) = 0.
      D(2) = 0.
      IF (Q.NE.0.) D(2) = 6.*Q*Q2/(QR*QR)
C
      IF (M.LT.2) GO TO 40
      DO 30 I=2,M
        P = Q
        Q = R
        R = X(I+2) - X(I+1)
        P2 = Q2
        Q2 = R2
        R2 = R*R
        PQ = QR
        QR = Q + R
        IF (Q) 20, 10, 20
   10   D(I+1) = 0.
        E(I) = 0.
        F(I-1) = 0.
        GO TO 30
   20   Q3 = Q2*Q
        PR = P*R
        PQQR = PQ*QR
        D(I+1) = 6.*Q3/(QR*QR)
        D(I) = D(I) + (Q+Q)*(15.*PR*PR+(P+R)*Q*(20.*PR+7.*Q2)+Q2*(8.*
     *   (P2+R2)+21.*PR+Q2+Q2))/(PQQR*PQQR)
        D(I-1) = D(I-1) + 6.*Q3/(PQ*PQ)
        E(I) = Q2*(P*QR+3.*PQ*(QR+R+R))/(PQQR*QR)
        E(I-1) = E(I-1) + Q2*(R*PQ+3.*QR*(PQ+P+P))/(PQQR*PQ)
        F(I-1) = Q3/PQQR
   30 CONTINUE
C
   40 IF (R.NE.0.) D(M) = D(M) + 6.*R*R2/(QR*QR)
C
C     FIRST AND SECOND ORDER DIVIDED DIFFERENCES OF THE GIVEN FUNCTION
C     VALUES, STORED IN B FROM 2 TO N AND IN C FROM 3 TO N
C     RESPECTIVELY. CARE IS TAKEN OF DOUBLE AND TRIPLE KNOTS.
C
      DO 60 I=2,N
        IF (X(I).NE.X(I-1)) GO TO 50
        B(I) = Y(I)
        Y(I) = Y(I-1)
        GO TO 60
   50   B(I) = (Y(I)-Y(I-1))/(X(I)-X(I-1))
   60 CONTINUE
      DO 80 I=3,N
        IF (X(I).NE.X(I-2)) GO TO 70
        C(I) = B(I)*0.5
        B(I) = B(I-1)
        GO TO 80
   70   C(I) = (B(I)-B(I-1))/(X(I)-X(I-2))
   80 CONTINUE
C
C     SOLVE THE LINEAR SYSTEM WITH C(I+2) - C(I+1) AS RIGHT-HAND SIDE.
C
      IF (M.LT.2) GO TO 100
      P = 0.
      C(1) = 0.
      E(M) = 0.
      F(1) = 0.
      F(M-1) = 0.
      F(M) = 0.
      C(2) = C(4) - C(3)
      D(2) = 1./D(2)
C
      IF (M.LT.3) GO TO 100
      DO 90 I=3,M
        Q = D(I-1)*E(I-1)
        D(I) = 1./(D(I)-P*F(I-2)-Q*E(I-1))
        E(I) = E(I) - Q*F(I-1)
        C(I) = C(I+2) - C(I+1) - P*C(I-2) - Q*C(I-1)
        P = D(I-1)*F(I-1)
   90 CONTINUE
C
  100 I = N - 1
      C(N-1) = 0.
      C(N) = 0.
      IF (N.LT.4) GO TO 120
      DO 110 M=4,N
C        I = N-2, ..., 2
        I = I - 1
        C(I) = (C(I)-E(I)*C(I+1)-F(I)*C(I+2))*D(I)
  110 CONTINUE
C
C     INTEGRATE THE THIRD DERIVATIVE OF S(X).
C
  120 M = N - 1
      Q = X(2) - X(1)
      R = X(3) - X(2)
      B1 = B(2)
      Q3 = Q*Q*Q
      QR = Q + R
      IF (QR) 140, 130, 140
  130 V = 0.
      T = 0.
      GO TO 150
  140 V = C(2)/QR
      T = V
  150 F(1) = 0.
      IF (Q.NE.0.) F(1) = V/Q
      DO 180 I=2,M
        P = Q
        Q = R
        R = 0.
        IF (I.NE.M) R = X(I+2) - X(I+1)
        P3 = Q3
        Q3 = Q*Q*Q
        PQ = QR
        QR = Q + R
        S = T
        T = 0.
        IF (QR.NE.0.) T = (C(I+1)-C(I))/QR
        U = V
        V = T - S
        IF (PQ) 170, 160, 170
  160   C(I) = C(I-1)
        D(I) = 0.
        E(I) = 0.
        F(I) = 0.
        GO TO 180
  170   F(I) = F(I-1)
        IF (Q.NE.0.) F(I) = V/Q
        E(I) = 5.*S
        D(I) = 10.*(C(I)-Q*S)
        C(I) = D(I)*(P-Q) + (B(I+1)-B(I)+(U-E(I))*P3-(V+E(I))*Q3)/PQ
        B(I) = (P*(B(I+1)-V*Q3)+Q*(B(I)-U*P3))/PQ -
     *   P*Q*(D(I)+E(I)*(Q-P))
  180 CONTINUE
C
C     END POINTS X(1) AND X(N).
C
      P = X(2) - X(1)
      S = F(1)*P*P*P
      E(1) = 0.
      D(1) = 0.
      C(1) = C(2) - 10.*S
      B(1) = B1 - (C(1)+S)*P
C
      Q = X(N) - X(N-1)
      T = F(N-1)*Q*Q*Q
      E(N) = 0.
      D(N) = 0.
      C(N) = C(N-1) + 10.*T
      B(N) = B(N) + (C(N)-T)*Q
  190 RETURN
      END
