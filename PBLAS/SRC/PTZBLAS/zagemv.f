      SUBROUTINE ZAGEMV( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y,
     $                   INCY )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        TRANS
      INTEGER            INCX, INCY, LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Y( * )
      COMPLEX*16         A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZAGEMV performs one of the matrix-vector operations
*
*     y := abs( alpha )*abs( A )*abs( x )+ abs( beta*y ),
*
*     or
*
*     y := abs( alpha )*abs( A' )*abs( x ) + abs( beta*y ),
*
*     or
*
*     y := abs( alpha )*abs( conjg( A' ) )*abs( x ) + abs( beta*y ),
*
*  where  alpha  and  beta  are real scalars, y is a real vector, x is a
*  vector and A is an m by n matrix.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          On entry,  TRANS  specifies the  operation to be performed as
*          follows:
*
*             TRANS = 'N' or 'n':
*                y := abs( alpha )*abs( A )*abs( x )+ abs( beta*y )
*
*             TRANS = 'T' or 't':
*                y := abs( alpha )*abs( A' )*abs( x ) + abs( beta*y )
*
*             TRANS = 'C' or 'c':
*                y := abs( alpha )*abs( conjg( A' ) )*abs( x ) +
*                     abs( beta*y )
*
*  M       (input) INTEGER
*          On entry, M  specifies the number of rows of the matrix  A. M
*          must be at least zero.
*
*  N       (input) INTEGER
*          On entry, N  specifies the number of columns of the matrix A.
*          N must be at least zero.
*
*  ALPHA   (input) DOUBLE PRECISION
*          On entry, ALPHA specifies the real scalar alpha.
*
*  A       (input) COMPLEX*16 array of dimension ( LDA, n ).
*          On entry, A  is an array of dimension ( LDA, N ). The leading
*          m by n part of the array  A  must contain the matrix of coef-
*          ficients.
*
*  LDA     (input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least max( 1, M ).
*
*  X       (input) COMPLEX*16 array of dimension at least
*          ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' and  at
*          least ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.  Before entry,
*          the incremented array X must contain the vector x.
*
*  INCX    (input) INTEGER
*          On entry, INCX specifies the increment for the elements of X.
*          INCX must not be zero.
*
*  BETA    (input) DOUBLE PRECISION
*          On entry,  BETA  specifies the real scalar beta. When BETA is
*          supplied as zero then Y need not be set on input.
*
*  Y       (input/output) DOUBLE PRECISION array of dimension at least
*          ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' and  at
*          least ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.  Before  entry
*          with BETA non-zero, the incremented array  Y must contain the
*          vector y. On exit, the incremented array  Y is overwritten by
*          the updated vector y.
*
*  INCY    (input) INTEGER
*          On entry, INCY specifies the increment for the elements of Y.
*          INCY must not be zero.
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
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
      DOUBLE PRECISION   ABSX, TALPHA, TEMP
      COMPLEX*16         ZDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( .NOT.LSAME( TRANS, 'N' ) .AND.
     $    .NOT.LSAME( TRANS, 'T' ) .AND.
     $    .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = 1
      ELSE IF( M.LT.0 ) THEN
         INFO = 2
      ELSE IF( N.LT.0 ) THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 ) THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 ) THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZAGEMV', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 ) THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 ) THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := abs( beta*y ).
*
      IF( INCY.EQ.1 ) THEN
         IF( BETA.EQ.ZERO ) THEN
            DO 10, I = 1, LENY
               Y( I ) = ZERO
   10       CONTINUE
         ELSE IF( BETA.EQ.ONE ) THEN
            DO 20, I = 1, LENY
               Y( I ) = ABS( Y( I ) )
   20       CONTINUE
         ELSE
            DO 30, I = 1, LENY
               Y( I ) = ABS( BETA * Y( I ) )
   30       CONTINUE
         END IF
      ELSE
         IY = KY
         IF( BETA.EQ.ZERO ) THEN
            DO 40, I = 1, LENY
               Y( IY ) = ZERO
               IY      = IY + INCY
   40       CONTINUE
         ELSE IF( BETA.EQ.ONE ) THEN
            DO 50, I = 1, LENY
               Y( IY ) = ABS( Y( IY ) )
               IY      = IY + INCY
   50       CONTINUE
         ELSE
            DO 60, I = 1, LENY
               Y( IY ) = ABS( BETA * Y( IY ) )
               IY      = IY + INCY
   60       CONTINUE
         END IF
      END IF
*
      IF( ALPHA.EQ.ZERO )
     $   RETURN
*
      TALPHA = ABS( ALPHA )
*
      IF( LSAME( TRANS, 'N' ) ) THEN
*
*        Form  y := abs( alpha ) * abs( A ) * abs( x ) + y.
*
         JX = KX
         IF( INCY.EQ.1 ) THEN
            DO 80, J = 1, N
               ABSX = CABS1( X( JX ) )
               IF( ABSX.NE.ZERO ) THEN
                  TEMP = TALPHA * ABSX
                  DO 70, I = 1, M
                     Y( I ) = Y( I ) + TEMP * CABS1( A( I, J ) )
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         ELSE
            DO 100, J = 1, N
               ABSX = CABS1( X( JX ) )
               IF( ABSX.NE.ZERO ) THEN
                  TEMP = TALPHA * ABSX
                  IY   = KY
                  DO 90, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP * CABS1( A( I, J ) )
                     IY      = IY      + INCY
   90             CONTINUE
               END IF
               JX = JX + INCX
  100       CONTINUE
         END IF
*
      ELSE
*
*        Form  y := abs( alpha ) * abs( A' ) * abs( x ) + y.
*
         JY = KY
         IF( INCX.EQ.1 ) THEN
            DO 120, J = 1, N
               TEMP = ZERO
               DO 110, I = 1, M
                  TEMP = TEMP + CABS1( A( I, J ) ) * CABS1( X( I ) )
  110          CONTINUE
               Y( JY ) = Y( JY ) + TALPHA * TEMP
               JY      = JY      + INCY
  120       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 130, I = 1, M
                  TEMP = TEMP + CABS1( A( I, J ) ) * CABS1( X( IX ) )
                  IX   = IX   + INCX
  130          CONTINUE
               Y( JY ) = Y( JY ) + TALPHA * TEMP
               JY      = JY      + INCY
  140       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZAGEMV
*
      END
