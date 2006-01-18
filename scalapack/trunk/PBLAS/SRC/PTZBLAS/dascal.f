      SUBROUTINE DASCAL( N, ALPHA, X, INCX )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DASCAL performs the following operation:
*
*                   x := abs( alpha ) * abs( x ),
*
*  where alpha is a scalar and x is an n vector.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          On entry, N specifies the length of the vector x. N  must  be
*          at least zero.
*
*  ALPHA   (input) DOUBLE PRECISION
*          On entry, ALPHA specifies the scalar alpha.
*
*  X       (input/output) DOUBLE PRECISION array of dimension at least
*          ( 1 + ( n - 1 )*abs( INCX ) ). Before entry,  the incremented
*          array  X  must  contain the vector x. On exit, entries of the
*          incremented array X are mutiplied by alpha in absolute value.
*
*  INCX    (input) INTEGER
*          On entry, INCX specifies the increment for the elements of X.
*          INCX must not be zero.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, IX, M, MP1
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MOD
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = 1
      ELSE IF( INCX.EQ.0 ) THEN
         INFO = 4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DASCAL', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.LE.0 )
     $   RETURN
*
*     Form  x := abs( alpha ) * abs( x )
*
      IF( INCX.EQ.1 )
     $   GO TO 40
*
*     code for increments not equal to 1
*
*     Set up the start point in  X.
*
      IF( INCX.GT.0 ) THEN
         IX = 1
      ELSE
         IX = 1 - ( N - 1 ) * INCX
      END IF
*
      IF( ALPHA.EQ.ZERO ) THEN
         DO 10 I = 1, N
           X( IX ) = ZERO
           IX = IX + INCX
   10    CONTINUE
      ELSE IF( ALPHA.EQ.ONE ) THEN
         DO 20 I = 1, N
           X( IX ) = ABS( X( IX ) )
           IX = IX + INCX
   20    CONTINUE
      ELSE
         DO 30 I = 1, N
           X( IX ) = ABS( ALPHA * X( IX ) )
           IX = IX + INCX
   30    CONTINUE
      END IF
*
      RETURN
*
*     code for increment equal to 1
*
*     clean-up loop
*
   40 M = MOD( N, 4 )
*
      IF( M.EQ.0 )
     $   GO TO 80
*
      IF( ALPHA.EQ.ZERO ) THEN
         DO 50 I = 1, M
           X( I ) = ZERO
   50    CONTINUE
      ELSE IF( ALPHA.EQ.ONE ) THEN
         DO 60 I = 1, M
           X( I ) = ABS( X( I ) )
   60    CONTINUE
      ELSE
         DO 70 I = 1, M
           X( I ) = ABS( ALPHA * X( I ) )
   70    CONTINUE
      END IF
*
      IF( N.LT.4 )
     $   RETURN
*
   80 MP1 = M + 1
*
      IF( ALPHA.EQ.ZERO ) THEN
         DO 90 I = MP1, N, 4
            X( I     ) = ZERO
            X( I + 1 ) = ZERO
            X( I + 2 ) = ZERO
            X( I + 3 ) = ZERO
   90    CONTINUE
      ELSE IF( ALPHA.EQ.ONE ) THEN
         DO 100 I = MP1, N, 4
            X( I     ) = ABS( X( I     ) )
            X( I + 1 ) = ABS( X( I + 1 ) )
            X( I + 2 ) = ABS( X( I + 2 ) )
            X( I + 3 ) = ABS( X( I + 3 ) )
  100    CONTINUE
      ELSE
         DO 110 I = MP1, N, 4
            X( I     ) = ABS( ALPHA * X( I     ) )
            X( I + 1 ) = ABS( ALPHA * X( I + 1 ) )
            X( I + 2 ) = ABS( ALPHA * X( I + 2 ) )
            X( I + 3 ) = ABS( ALPHA * X( I + 3 ) )
  110    CONTINUE
      END IF
*
      RETURN
*
*     End of DASCAL
*
      END
