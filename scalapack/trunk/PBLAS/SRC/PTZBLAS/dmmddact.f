      SUBROUTINE DMMDDACT( M, N, ALPHA, A, LDA, BETA, B, LDB )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, M, N
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DMMDDACT performs the following operation:
*
*     A := alpha * A + beta * B',
*
*  where alpha, beta are scalars; A is an m by n matrix and B is an n by
*  m matrix.
*
*  Arguments
*  =========
*
*  M       (local input) INTEGER
*          On entry, M  specifies the number of rows of A and the number
*          of columns of B. M must be at least zero.
*
*  N       (local input) INTEGER
*          On entry, N  specifies the number of rows of B and the number
*          of columns of A. N must be at least zero.
*
*  ALPHA   (local input) DOUBLE PRECISION
*          On entry,  ALPHA  specifies the scalar alpha. When  ALPHA  is
*          supplied as zero then the local entries of the array  A  need
*          not be set on input.
*
*  A       (local input/local output) DOUBLE PRECISION array
*          On entry, A is an array of dimension ( LDA, N ). On exit, the
*          leading n by m part of B has been added into the leading m by
*          n part of A.
*
*  LDA     (local input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least max( 1, M ).
*
*  BETA    (local input) DOUBLE PRECISION
*          On entry,  BETA  specifies the scalar beta. When BETA is sup-
*          plied as zero then the local entries of the array B need  not
*          be set on input.
*
*  B       (local input) DOUBLE PRECISION array
*          On entry, B is an array of dimension ( LDB, M ).
*
*  LDB     (local input) INTEGER
*          On entry, LDB specifies the leading dimension of the array B.
*          LDB must be at least max( 1, N ).
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
      INTEGER            I, J
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DSCAL
*     ..
*     .. Executable Statements ..
*
      IF( M.GE.N ) THEN
         IF( BETA.EQ.ONE ) THEN
            IF( ALPHA.EQ.ZERO ) THEN
               DO 20 J = 1, N
                  CALL DCOPY( M, B( J, 1 ), LDB, A( 1, J ), 1 )
*                 DO 10 I = 1, M
*                    A( I, J ) = B( J, I )
*  10             CONTINUE
   20          CONTINUE
            ELSE IF( ALPHA.NE.ONE ) THEN
               DO 40 J = 1, N
                  DO 30 I = 1, M
                     A( I, J ) = B( J, I ) + ALPHA * A( I, J )
   30             CONTINUE
   40          CONTINUE
            ELSE
               DO 60 J = 1, N
                  CALL DAXPY( M, ONE, B( J, 1 ), LDB, A( 1, J ), 1 )
*                 DO 50 I = 1, M
*                    A( I, J ) = B( J, I ) + A( I, J )
*  50             CONTINUE
   60          CONTINUE
            END IF
         ELSE IF( BETA.NE.ZERO ) THEN
            IF( ALPHA.EQ.ZERO ) THEN
               DO 80 J = 1, N
                  DO 70 I = 1, M
                     A( I, J ) = BETA * B( J, I )
   70             CONTINUE
   80          CONTINUE
            ELSE IF( ALPHA.NE.ONE ) THEN
               DO 100 J = 1, N
                  DO 90 I = 1, M
                     A( I, J ) = BETA * B( J, I ) + ALPHA * A( I, J )
   90             CONTINUE
  100          CONTINUE
            ELSE
               DO 120 J = 1, N
                  CALL DAXPY( M, BETA, B( J, 1 ), LDB, A( 1, J ), 1 )
*                 DO 110 I = 1, M
*                    A( I, J ) = BETA * B( J, I ) + A( I, J )
* 110             CONTINUE
  120          CONTINUE
            END IF
         ELSE
            IF( ALPHA.EQ.ZERO ) THEN
               DO 140 J = 1, N
                  DO 130 I = 1, M
                     A( I, J ) = ZERO
  130             CONTINUE
  140          CONTINUE
            ELSE IF( ALPHA.NE.ONE ) THEN
               DO 160 J = 1, N
                  CALL DSCAL( M, ALPHA, A( 1, J ), 1 )
*                 DO 150 I = 1, M
*                    A( I, J ) = ALPHA * A( I, J )
* 150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( BETA.EQ.ONE ) THEN
            IF( ALPHA.EQ.ZERO ) THEN
               DO 180 J = 1, M
                  CALL DCOPY( N, B( 1, J ), 1, A( J, 1 ), LDA )
*                 DO 170 I = 1, N
*                    A( J, I ) = B( I, J )
* 170             CONTINUE
  180          CONTINUE
            ELSE IF( ALPHA.NE.ONE ) THEN
               DO 200 J = 1, M
                  DO 190 I = 1, N
                     A( J, I ) = B( I, J ) + ALPHA * A( J, I )
  190             CONTINUE
  200          CONTINUE
            ELSE
               DO 220 J = 1, M
                  CALL DAXPY( N, ONE, B( 1, J ), 1, A( J, 1 ), LDA )
*                 DO 210 I = 1, N
*                    A( J, I ) = B( I, J ) + A( J, I )
* 210             CONTINUE
  220          CONTINUE
            END IF
         ELSE IF( BETA.NE.ZERO ) THEN
            IF( ALPHA.EQ.ZERO ) THEN
               DO 240 J = 1, M
                  DO 230 I = 1, N
                     A( J, I ) = BETA * B( I, J )
  230             CONTINUE
  240          CONTINUE
            ELSE IF( ALPHA.NE.ONE ) THEN
               DO 260 J = 1, M
                  DO 250 I = 1, N
                     A( J, I ) = BETA * B( I, J ) + ALPHA * A( J, I )
  250             CONTINUE
  260          CONTINUE
            ELSE
               DO 280 J = 1, M
                  CALL DAXPY( N, BETA, B( 1, J ), 1, A( J, 1 ), LDA )
*                 DO 270 I = 1, N
*                    A( J, I ) = BETA * B( I, J ) + A( J, I )
* 270             CONTINUE
  280          CONTINUE
            END IF
         ELSE
            IF( ALPHA.EQ.ZERO ) THEN
               DO 300 J = 1, N
                  DO 290 I = 1, M
                     A( I, J ) = ZERO
  290             CONTINUE
  300          CONTINUE
            ELSE IF( ALPHA.NE.ONE ) THEN
               DO 320 J = 1, N
                  CALL DSCAL( M, ALPHA, A( 1, J ), 1 )
*                 DO 310 I = 1, M
*                    A( I, J ) = ALPHA * A( I, J )
* 310             CONTINUE
  320          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DMMDDACT
*
      END
