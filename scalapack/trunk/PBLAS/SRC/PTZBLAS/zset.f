      SUBROUTINE ZSET( N, ALPHA, X, INCX )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*16         ALPHA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZSET sets the entries of an n vector x to the scalar alpha.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          On entry, N specifies the length of the vector x. N  must  be
*          at least zero.
*
*  ALPHA   (input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  X       (input/output) COMPLEX*16 array of dimension at least
*          ( 1 + ( n - 1 )*abs( INCX ) ). Before entry,  the incremented
*          array  X  must  contain the vector x. On exit, entries of the
*          incremented array X are set to alpha.
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
*     .. Local Scalars ..
      INTEGER            I, INFO, IX, M, MP1
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
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
         CALL XERBLA( 'ZSET', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.LE.0 )
     $   RETURN
*
*     Form  x := alpha
*
      IF( INCX.EQ.1 )
     $   GO TO 20
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
      DO 10 I = 1, N
        X( IX ) = ALPHA
        IX = IX + INCX
   10 CONTINUE
*
      RETURN
*
*     code for increment equal to 1
*
*     clean-up loop
*
   20 M = MOD( N, 4 )
*
      IF( M.EQ.0 )
     $   GO TO 40
*
      DO 30 I = 1, M
        X( I ) = ALPHA
   30 CONTINUE
      IF( N.LT.4 )
     $   RETURN
*
   40 MP1 = M + 1
      DO 50 I = MP1, N, 4
         X( I     ) = ALPHA
         X( I + 1 ) = ALPHA
         X( I + 2 ) = ALPHA
         X( I + 3 ) = ALPHA
   50 CONTINUE
*
      RETURN
*
*     End of ZSET
*
      END
