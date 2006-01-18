      SUBROUTINE PBDVECADD( ICONTXT, MODE, N, ALPHA, X, INCX, BETA, Y,
     $                      INCY )
*
*  -- PB-BLAS routine (version 2.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory.
*     April 28, 1996
*
*     .. Scalar Arguments ..
      CHARACTER*1           MODE
      INTEGER               ICONTXT, INCX, INCY, N
      DOUBLE PRECISION      ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION      X( * ), Y( * )
*
*     ..
*
*  Purpose
*  =======
*
*  PBDVECADD performs a vector X to be added to Y
*    Y := alpha*op(X) + beta*Y,
*  where alpha and beta are scalars, and X and Y are n vectors,
*  and op(X) = X**H if MODE = 'C',
*
*  Arguments
*  =========
*
*  ICONTXT (input) INTEGER
*          ICONTXT is the BLACS mechanism for partitioning communication
*          space.  A defining property of a context is that a message in
*          a context cannot be sent or received in another context.  The
*          BLACS context includes the definition of a grid, and each
*          process' coordinates in it.
*
*  MODE   (input) CHARACTER*1
*          Specifies the transposed, or conjugate transposed vector X
*          to be added to the vector Y
*          = 'C': Conjugate vector X is added for complex data set.
*                 Y = alpha * X**H + beta * Y
*          ELSE : Vector X is added.  Y = alpha*X + beta*Y
*                 if MODE = 'V', BLAS routine may be used.
*
*  N       (input) INTEGER
*          The number of elements of the vectors X and Y to be added.
*          N >= 0.
*
*  ALPHA   (input) DOUBLE PRECISION
*          ALPHA specifies the scalar alpha.
*
*  X       (input) DOUBLE PRECISION array of DIMENSION at least
*          ( 1 + ( N - 1 )*abs( INCX ) )
*          The incremented array X must contain the vector X.
*
*  INCX    (input) INTEGER
*          INCX specifies the increment for the elements of X.
*          INCX <> 0.
*
*  BETA    (input) DOUBLE PRECISION
*          BETA specifies the scalar beta.
*
*  Y       (input/output) DOUBLE PRECISION array of DIMENSION at least
*          ( 1 + ( N - 1 )*abs( INCY ) )
*          On entry with BETA non-zero, the incremented array Y must
*          contain the vector Y.
*          On exit, Y is overwritten by the updated vector Y.
*
*  INCY  - (input) INTEGER
*          INCY specifies the increment for the elements of Y.
*          INCY <> 0.
*
*  =====================================================================
*
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO,          ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0)
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX, IY
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DCOPY, DAXPY
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 .OR. ( ALPHA.EQ.ZERO .AND. BETA.EQ.ONE ) ) RETURN
*
      IF( ALPHA.EQ.ZERO ) THEN
         IF( BETA.EQ.ZERO ) THEN
            IF( INCY.EQ.1 ) THEN
               DO 10 I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               IY = 1
               DO 20 I = 1, N
                  Y( IY ) = ZERO
                  IY = IY + INCY
   20          CONTINUE
            END IF
*
         ELSE
            IF( LSAME( MODE, 'V' ) ) THEN
               CALL DSCAL( N, BETA, Y, INCY )
            ELSE IF( INCY.EQ.1 ) THEN
               DO 30 I = 1, N
                  Y( I ) = BETA * Y( I )
   30          CONTINUE
            ELSE
               IY = 1
               DO 40 I = 1, N
                  Y( IY ) = BETA * Y( IY )
                  IY = IY + INCY
   40          CONTINUE
            END IF
         END IF
*
      ELSE
         IF( ALPHA.EQ.ONE ) THEN
            IF( BETA.EQ.ZERO ) THEN
               IF( LSAME( MODE, 'V' ) ) THEN
                  CALL DCOPY( N, X, INCX, Y, INCY )
               ELSE IF( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
                  DO 50 I = 1, N
                     Y( I ) = X( I )
   50             CONTINUE
               ELSE
                  IX = 1
                  IY = 1
                  DO 60 I = 1, N
                     Y( IY ) = X( IX )
                     IX = IX + INCX
                     IY = IY + INCY
   60             CONTINUE
               END IF
*
            ELSE IF( BETA.EQ.ONE ) THEN
               IF( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
                  DO 70 I = 1, N
                     Y( I ) = X( I ) + Y( I )
   70             CONTINUE
               ELSE
                  IX = 1
                  IY = 1
                  DO 80 I = 1, N
                     Y( IY ) = X( IX ) + Y( IY )
                     IX = IX + INCX
                     IY = IY + INCY
   80             CONTINUE
               END IF
*
            ELSE
               IF( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
                  DO 90 I = 1, N
                     Y( I ) = X( I ) + BETA * Y( I )
   90             CONTINUE
               ELSE
                  IX = 1
                  IY = 1
                  DO 100 I = 1, N
                     Y( IY ) = X( IX ) + BETA * Y( IY )
                     IX = IX + INCX
                     IY = IY + INCY
  100             CONTINUE
               END IF
            END IF
*
         ELSE
            IF( BETA.EQ.ZERO ) THEN
               IF( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
                  DO 110 I = 1, N
                     Y( I ) = ALPHA * X( I )
  110             CONTINUE
               ELSE
                  IX = 1
                  IY = 1
                  DO 120 I = 1, N
                     Y( IY ) = X( IX )
                     IX = IX + INCX
                     IY = IY + INCY
  120             CONTINUE
               END IF
*
            ELSE IF( BETA.EQ.ONE ) THEN
               IF( LSAME( MODE, 'V' ) ) THEN
                  CALL DAXPY( N, ALPHA, X, INCX, Y, INCY )
               ELSE IF( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
                  DO 130 I = 1, N
                     Y( I ) = ALPHA * X( I ) + Y( I )
  130             CONTINUE
               ELSE
                  IX = 1
                  IY = 1
                  DO 140 I = 1, N
                     Y( IY ) = ALPHA * X( IX ) + Y( IY )
                     IX = IX + INCX
                     IY = IY + INCY
  140             CONTINUE
               END IF
*
            ELSE
               IF( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
                  DO 150 I = 1, N
                     Y( I ) = ALPHA * X( I ) + BETA * Y( I )
  150             CONTINUE
               ELSE
                  IX = 1
                  IY = 1
                  DO 160 I = 1, N
                     Y( IY ) = ALPHA * X( IX ) + BETA * Y( IY )
                     IX = IX + INCX
                     IY = IY + INCY
  160             CONTINUE
               END IF
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of PBDVECADD
*
      END
