      SUBROUTINE CATRMV( UPLO, TRANS, DIAG, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        DIAG, TRANS, UPLO
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
*  CATRMV performs one of the matrix-vector operations
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
*  vector and A is an n by n unit or non-unit, upper or lower triangular
*  matrix.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          On entry,  UPLO  specifies  whether the matrix is an upper or
*          lower triangular matrix as follows:
*
*             UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*             UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*  TRANS   (input) CHARACTER*1
*          On entry,  TRANS  specifies  the operation to be performed as
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
*  DIAG    (input) CHARACTER*1
*          On entry, DIAG  specifies whether or not A is unit triangular
*          as follows:
*
*             DIAG = 'U' or 'u'  A is assumed to be unit triangular.
*
*             DIAG = 'N' or 'n'  A is not assumed to be unit triangular.
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
*          A must contain the upper triangular part of the matrix A  and
*          the  strictly  lower  triangular part of A is not referenced.
*          When  UPLO = 'L' or 'l', the leading n by n part of the array
*          A  must contain the lower triangular part of the matrix A and
*          the  strictly  upper trapezoidal part of A is not referenced.
*          Note that when  DIAG = 'U' or 'u', the diagonal elements of A
*          are not referenced either, but are assumed to be unity.
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
      LOGICAL            NOUNIT
      REAL               ABSX, TALPHA, TEMP
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
      INTRINSIC          ABS, AIMAG, MAX, REAL
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
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 7
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 9
      ELSE IF( INCY.EQ.0 ) THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CATRMV', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
      NOUNIT = LSAME( DIAG , 'N' )
*
*     Set up the start points in  X  and  Y.
*
      IF( INCX.GT.0 ) THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 ) * INCX
      END IF
      IF( INCY.GT.0 ) THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 ) * INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := abs( beta*y ).
*
      IF( INCY.EQ.1 ) THEN
         IF( BETA.EQ.ZERO ) THEN
            DO 10, I = 1, N
               Y( I ) = ZERO
   10       CONTINUE
         ELSE IF( BETA.EQ.ONE ) THEN
            DO 20, I = 1, N
               Y( I ) = ABS( Y( I ) )
   20       CONTINUE
         ELSE
            DO 30, I = 1, N
               Y( I ) = ABS( BETA * Y( I ) )
   30       CONTINUE
         END IF
      ELSE
         IY = KY
         IF( BETA.EQ.ZERO ) THEN
            DO 40, I = 1, N
               Y( IY ) = ZERO
               IY      = IY + INCY
   40       CONTINUE
         ELSE IF( BETA.EQ.ONE ) THEN
            DO 50, I = 1, N
               Y( IY ) = ABS( Y( IY ) )
               IY      = IY + INCY
   50       CONTINUE
         ELSE
            DO 60, I = 1, N
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
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := abs( alpha ) * abs( A ) * abs( x ) + y.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            JX = KX
            IF( INCY.EQ.1 ) THEN
               DO 80, J = 1, N
                  ABSX = CABS1( X( JX ) )
                  IF( ABSX.NE.ZERO ) THEN
                     TEMP = TALPHA * ABSX
                     DO 70, I = 1, J - 1
                        Y( I ) = Y( I ) + TEMP * CABS1( A( I, J ) )
   70                CONTINUE
*
                     IF( NOUNIT ) THEN
                        Y( J ) = Y( J ) + TEMP * CABS1( A( J, J ) )
                     ELSE
                        Y( J ) = Y( J ) + TEMP
                     END IF
                  END IF
                  JX = JX + INCX
   80          CONTINUE
*
            ELSE
*
               DO 100, J = 1, N
                  ABSX = CABS1( X( JX ) )
                  IF( ABSX.NE.ZERO ) THEN
                     TEMP = TALPHA * ABSX
                     IY   = KY
                     DO 90, I = 1, J - 1
                        Y( IY ) = Y( IY ) + TEMP * CABS1( A( I, J ) )
                        IY      = IY      + INCY
   90                CONTINUE
*
                     IF( NOUNIT ) THEN
                        Y( IY ) = Y( IY ) + TEMP * CABS1( A( J, J ) )
                     ELSE
                        Y( IY ) = Y( IY ) + TEMP
                     END IF
                  END IF
                  JX = JX + INCX
  100          CONTINUE
*
            END IF
*
         ELSE
*
            JX = KX
            IF( INCY.EQ.1 ) THEN
               DO 120, J = 1, N
                  ABSX = CABS1( X( JX ) )
                  IF( ABSX.NE.ZERO ) THEN
*
                     TEMP = TALPHA * ABSX
*
                     IF( NOUNIT ) THEN
                        Y( J ) = Y( J ) + TEMP * CABS1( A( J, J ) )
                     ELSE
                        Y( J ) = Y( J ) + TEMP
                     END IF
*
                     DO 110, I = J + 1, N
                        Y( I ) = Y( I ) + TEMP * CABS1( A( I, J ) )
  110                CONTINUE
                  END IF
                  JX = JX + INCX
  120          CONTINUE
*
            ELSE
*
               DO 140, J = 1, N
                  ABSX = CABS1( X( JX ) )
                  IF( ABSX.NE.ZERO ) THEN
                     TEMP = TALPHA * ABSX
                     IY   = KY + ( J - 1 ) * INCY
*
                     IF( NOUNIT ) THEN
                        Y( IY ) = Y( IY ) + TEMP * CABS1( A( J, J ) )
                     ELSE
                        Y( IY ) = Y( IY ) + TEMP
                     END IF
*
                     DO 130, I = J + 1, N
                        IY      = IY      + INCY
                        Y( IY ) = Y( IY ) + TEMP * CABS1( A( I, J ) )
  130                CONTINUE
                  END IF
                  JX = JX + INCX
  140          CONTINUE
*
            END IF
*
         END IF
*
      ELSE
*
*        Form  y := abs( alpha ) * abs( A' ) * abs( x ) + y.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            JY = KY
            IF( INCX.EQ.1 ) THEN
               DO 160, J = 1, N
*
                  TEMP = ZERO
*
                  DO 150, I = 1, J - 1
                     TEMP = TEMP + CABS1( A( I, J ) ) * CABS1( X( I ) )
  150             CONTINUE
*
                  IF( NOUNIT ) THEN
                     TEMP = TEMP + CABS1( A( J, J ) ) * CABS1( X( J ) )
                  ELSE
                     TEMP = TEMP + CABS1( X( J ) )
                  END IF
*
                  Y( JY ) = Y( JY ) + TALPHA * TEMP
                  JY      = JY      + INCY
*
  160          CONTINUE
*
            ELSE
*
               DO 180, J = 1, N
                  TEMP = ZERO
                  IX   = KX
                  DO 170, I = 1, J - 1
                     TEMP = TEMP + CABS1( A( I, J ) ) * CABS1( X( IX ) )
                     IX   = IX   + INCX
  170             CONTINUE
*
                  IF( NOUNIT ) THEN
                     TEMP = TEMP + CABS1( A( J, J ) ) * CABS1( X( IX ) )
                  ELSE
                     TEMP = TEMP + CABS1( X( IX ) )
                  END IF
*
                  Y( JY ) = Y( JY ) + TALPHA * TEMP
                  JY      = JY      + INCY
*
  180          CONTINUE
*
            END IF
*
         ELSE
*
            JY = KY
            IF( INCX.EQ.1 ) THEN
*
               DO 200, J = 1, N
*
                  IF( NOUNIT ) THEN
                     TEMP = CABS1( A( J, J ) ) * CABS1( X( J ) )
                  ELSE
                     TEMP = CABS1( X( J ) )
                  END IF
*
                  DO 190, I = J + 1, N
                     TEMP = TEMP + CABS1( A( I, J ) ) * CABS1( X( I ) )
  190             CONTINUE
*
                  Y( JY ) = Y( JY ) + TALPHA * TEMP
                  JY      = JY      + INCY
*
  200          CONTINUE
*
            ELSE
*
               DO 220, J = 1, N
*
                  IX   = KX + ( J - 1 ) * INCX
*
                  IF( NOUNIT ) THEN
                     TEMP = CABS1( A( J, J ) ) * CABS1( X( IX ) )
                  ELSE
                     TEMP = CABS1( X( IX ) )
                  END IF
*
                  DO 210, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + CABS1( A( I, J ) ) * CABS1( X( IX ) )
  210             CONTINUE
                  Y( JY ) = Y( JY ) + TALPHA * TEMP
                  JY      = JY      + INCY
  220          CONTINUE
            END IF
         END IF
*
      END IF
*
      RETURN
*
*     End of CATRMV
*
      END
