      SUBROUTINE CAHEMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y,
     $                   INCY )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO
      INTEGER            INCX, INCY, LDA, N
      REAL               ALPHA, BETA
*     ..
*     .. Array Arguments ..
      REAL               Y( * )
      COMPLEX            A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  CAHEMV performs the following matrix-vector operation
*
*     y := abs( alpha )*abs( A )*abs( x )+ abs( beta*y ),
*
*  where  alpha  and  beta  are real scalars, y is a real vector, x is a
*  vector and A is an n by n Hermitian matrix.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          On entry, UPLO  specifies whether the upper or lower triangu-
*          lar part of the array A is to be referenced as follows:
*
*             UPLO = 'U' or 'u'   Only the upper triangular part of A is
*                                 to be referenced.
*             UPLO = 'L' or 'l'   Only the lower triangular part of A is
*                                 to be referenced.
*
*  N       (input) INTEGER
*          On entry, N specifies the order of the matrix A. N must be at
*          least zero.
*
*  ALPHA   (input) REAL
*          On entry, ALPHA specifies the real scalar alpha.
*
*  A       (input) COMPLEX array
*          On entry,  A  is an array of dimension (LDA,N).  Before entry
*          with UPLO = 'U' or 'u', the leading n by n part of the  array
*          A must contain the upper triangular part of the Hermitian ma-
*          trix and the strictly lower triangular part of A is not refe-
*          renced.  When  UPLO = 'L' or 'l',  the leading n by n part of
*          the array  A  must  contain  the lower triangular part of the
*          Hermitian matrix and the strictly upper trapezoidal part of A
*          is not referenced.
*          Note that the  imaginary parts  of the local entries  corres-
*          ponding to the offdiagonal elements of A need not  be set and
*          assumed to be zero.
*
*  LDA     (input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least max( 1, N ).
*
*  X       (input) COMPLEX array of dimension at least
*          ( 1 + ( n - 1 )*abs( INCX ) ).  Before entry, the incremented
*          array X must contain the vector x.
*
*  INCX    (input) INTEGER
*          On entry, INCX specifies the increment for the elements of X.
*          INCX must not be zero.
*
*  BETA    (input) REAL
*          On entry,  BETA  specifies the real scalar beta. When BETA is
*          supplied as zero then Y need not be set on input.
*
*  Y       (input/output) REAL array of dimension at least
*          ( 1 + ( n - 1 )*abs( INCY ) ).  Before  entry with  BETA non-
*          zero, the  incremented array  Y must contain the vector y. On
*          exit, the  incremented array  Y is overwritten by the updated
*          vector y.
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
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
      REAL               TALPHA, TEMP0, TEMP1, TEMP2
      COMPLEX            ZDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CONJG, MAX, REAL
*     ..
*     .. Statement Functions ..
      REAL               CABS1
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CAHEMV', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF( INCX.GT.0 ) THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 ) * INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 ) * INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
*     First form  y := abs( beta * y ).
*
      IF( BETA.NE.ONE ) THEN
         IF( INCY.EQ.1 ) THEN
            IF( BETA.EQ.ZERO ) THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = ABS( BETA * Y( I ) )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO ) THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = ABS( BETA * Y( IY ) )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
*
      IF( ALPHA.EQ.ZERO )
     $   RETURN
*
      TALPHA = ABS( ALPHA )
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Form  y  when A is stored in upper triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) ) THEN
            DO 60, J = 1, N
               TEMP1 = TALPHA * CABS1( X( J ) )
               TEMP2 = ZERO
               DO 50, I = 1, J - 1
                  TEMP0  = CABS1( A( I, J ) )
                  Y( I ) = Y( I ) + TEMP1 * TEMP0
                  TEMP2  = TEMP2  + TEMP0 * CABS1( X( I ) )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1 * ABS( REAL( A( J, J ) ) ) +
     $                  ALPHA * TEMP2
*
   60       CONTINUE
*
         ELSE
*
            JX = KX
            JY = KY
*
            DO 80, J = 1, N
               TEMP1 = TALPHA * CABS1( X( JX ) )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
*
               DO 70, I = 1, J - 1
                  TEMP0   = CABS1( A( I, J ) )
                  Y( IY ) = Y( IY ) + TEMP1 * TEMP0
                  TEMP2   = TEMP2   + TEMP0 * CABS1( X( IX ) )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1 * ABS( REAL( A( J, J ) ) ) +
     $                   ALPHA * TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
*
   80       CONTINUE
*
         END IF
*
      ELSE
*
*        Form  y  when A is stored in lower triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) ) THEN
*
            DO 100, J = 1, N
*
               TEMP1 = TALPHA * CABS1( X( J ) )
               TEMP2  = ZERO
               Y( J ) = Y( J ) + TEMP1 * ABS( REAL( A( J, J ) ) )
*
               DO 90, I = J + 1, N
                  TEMP0  = CABS1( A( I, J ) )
                  Y( I ) = Y( I ) + TEMP1 * TEMP0
                  TEMP2  = TEMP2  + TEMP0 * CABS1( X( I ) )
*
   90          CONTINUE
*
               Y( J ) = Y( J ) + ALPHA * TEMP2
*
  100       CONTINUE
*
         ELSE
*
            JX = KX
            JY = KY
*
            DO 120, J = 1, N
               TEMP1 = TALPHA * CABS1( X( JX ) )
               TEMP2   = ZERO
               Y( JY ) = Y( JY ) + TEMP1 * ABS( REAL( A( J, J ) ) )
               IX      = JX
               IY      = JY
*
               DO 110, I = J + 1, N
*
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  TEMP0   = CABS1( A( I, J ) )
                  Y( IY ) = Y( IY ) + TEMP1 * TEMP0
                  TEMP2   = TEMP2   + TEMP0 * CABS1( X( IX ) )
*
  110          CONTINUE
*
               Y( JY ) = Y( JY ) + ALPHA * TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
*
  120       CONTINUE
*
         END IF
*
      END IF
*
      RETURN
*
*     End of CAHEMV
*
      END
