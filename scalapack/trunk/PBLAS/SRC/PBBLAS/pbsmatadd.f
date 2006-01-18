      SUBROUTINE PBSMATADD( ICONTXT, MODE, M, N, ALPHA, A, LDA, BETA, B,
     $                      LDB )
*
*  -- PB-BLAS routine (version 2.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory.
*     April 28, 1996
*
*     .. Scalar Arguments ..
      CHARACTER*1        MODE
      INTEGER            ICONTXT, LDA, LDB, M, N
      REAL               ALPHA, BETA
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  PBSMATADD performs the matrix add operation B := alpha*A + beta*B,
*  where alpha and beta are scalars, and A and B are m-by-n
*  upper/lower trapezoidal matrices, or rectangular matrices.
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
*  MODE    (input) CHARACTER*1
*          Specifies the part of the matrix A, or (conjugate) transposed
*          matrix A to be added to the matrix B,
*          = 'U':  Upper triangular part
*                  up(B) = alpha*up(A) + beta*up(B)
*          = 'L':  Lower triangular part
*                  lo(B) = alpha*lo(A) + beta*lo(B)
*          = 'T':  Transposed matrix A
*                  B = alpha*A**T + beta*B
*          = 'C':  Conjugate transposed matrix A
*                  B = alpha*A**H + beta*B
*          Otherwise: B = alpha*A + beta*B
*              if M = LDA = LDB: use one BLAS loop
*              if MODE = 'V' :   columnwise copy using BLAS if possible
*              else :            use double loops
*
*  M       (input) INTEGER
*          M specifies the number of columns of the matrix A if
*          MODE != 'T'/'C', and it specifies the number of rows of the
*          matrix A otherwise.  It also specifies the number of rows of
*          the matrix B. M >= 0.
*
*  N       (input) INTEGER
*          N specifies the number of rows of the matrix A if
*          MODE != 'T'/'C', and it specifies the number of columns of
*          the matrix A otherwise. It also specifies the number of
*          columns of the matrix B. N >= 0.
*
*  ALPHA   (input) REAL
*          ALPHA specifies the scalar alpha.
*
*  A       (input) REAL array, dimension (LDA,N)
*          The m by n matrix A if MODE != 'T'/'C'.
*          If MODE = 'U', only the upper triangle or trapezoid is
*          accessed; if MODE = 'L', only the lower triangle or
*          trapezoid is accessed.  Otherwise all m-by-n data matrix
*          is accessed.
*          And the n by m matrix A if MODE = 'T'/'C'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M) if
*          MODE != 'T'/'C'.  And LDA >= max(1,N) if MODE = 'T'/'C'.
*
*  BETA    (input) REAL
*          BETA specifies the scalar beta.
*
*  B       (input) REAL array, dimension (LDB,N)
*          On exit, B = alpha*A + beta*B
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO,          ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0)
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL, SCOPY, SAXPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( M.LE.0 .OR. N.LE.0 .OR. ( ALPHA.EQ.ZERO.AND.BETA.EQ.ONE ) )
     $   RETURN
*
*     A is upper triangular or upper trapezoidal,
*
      IF( LSAME( MODE, 'U' ) ) THEN
         IF( ALPHA.EQ.ZERO ) THEN
            IF( BETA.EQ.ZERO ) THEN
               DO 20 J = 1, N
                  DO 10 I = 1, MIN( J, M )
                     B( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40 J = 1, N
                  DO 30 I = 1, MIN( J, M )
                     B( I, J ) = BETA * B( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
*
         ELSE IF( ALPHA.EQ.ONE ) THEN
            IF( BETA.EQ.ZERO ) THEN
               DO 60 J = 1, N
                  DO 50 I = 1, MIN( J, M )
                     B( I, J ) = A( I, J )
   50             CONTINUE
   60          CONTINUE
            ELSE IF( BETA.EQ.ONE ) THEN
               DO 80 J = 1, N
                  DO 70 I = 1, MIN( J, M )
                     B( I, J ) = A( I, J ) + B( I, J )
   70             CONTINUE
   80          CONTINUE
            ELSE
               DO 100 J = 1, N
                  DO 90 I = 1, MIN( J, M )
                     B( I, J ) = A( I, J ) + BETA * B( I, J )
   90             CONTINUE
  100          CONTINUE
            END IF
*
         ELSE
            IF( BETA.EQ.ZERO ) THEN
               DO 120 J = 1, N
                  DO 110 I = 1, MIN( J, M )
                     B( I, J ) = ALPHA * A( I, J )
  110             CONTINUE
  120          CONTINUE
            ELSE IF( BETA.EQ.ONE ) THEN
               DO 140 J = 1, N
                  DO 130 I = 1, MIN( J, M )
                     B( I, J ) = ALPHA * A( I, J ) + B( I, J )
  130             CONTINUE
  140          CONTINUE
            ELSE
               DO 160 J = 1, N
                  DO 150 I = 1, MIN( J, M )
                     B( I, J ) = ALPHA * A( I, J ) + BETA * B( I, J )
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
*
*     A is lower triangular or upper trapezoidal,
*
      ELSE IF( LSAME( MODE, 'L' ) ) THEN
         IF( ALPHA.EQ.ZERO ) THEN
            IF( BETA.EQ.ZERO ) THEN
               DO 180 J = 1, N
                  DO 170 I = J, M
                     B( I, J ) = ZERO
  170             CONTINUE
  180          CONTINUE
            ELSE
               DO 200 J = 1, N
                  DO 190 I = J, M
                     B( I, J ) = BETA * B( I, J )
  190             CONTINUE
  200          CONTINUE
            END IF
*
         ELSE IF( ALPHA.EQ.ONE ) THEN
            IF( BETA.EQ.ZERO ) THEN
               DO 220 J = 1, N
                  DO 210 I = J, M
                     B( I, J ) = A( I, J )
  210             CONTINUE
  220          CONTINUE
            ELSE IF( BETA.EQ.ONE ) THEN
               DO 240 J = 1, N
                  DO 230 I = J, M
                     B( I, J ) = A( I, J ) + B( I, J )
  230             CONTINUE
  240          CONTINUE
            ELSE
               DO 260 J = 1, N
                  DO 250 I = J, M
                     B( I, J ) = A( I, J ) + BETA * B( I, J )
  250             CONTINUE
  260          CONTINUE
            END IF
*
         ELSE
            IF( BETA.EQ.ZERO ) THEN
               DO 280 J = 1, N
                  DO 270 I = J, M
                     B( I, J ) = ALPHA * A( I, J )
  270             CONTINUE
  280          CONTINUE
            ELSE IF( BETA.EQ.ONE ) THEN
               DO 300 J = 1, N
                  DO 290 I = J, M
                     B( I, J ) = ALPHA * A( I, J ) + B( I, J )
  290             CONTINUE
  300          CONTINUE
            ELSE
               DO 320 J = 1, N
                  DO 310 I = J, M
                     B( I, J ) = ALPHA * A( I, J ) + BETA * B( I, J )
  310             CONTINUE
  320          CONTINUE
            END IF
         END IF
*
*     If MODE = 'Transpose'/'Conjugate'
*
      ELSE IF( LSAME( MODE, 'T' ) .OR. LSAME( MODE, 'C' ) ) THEN
         IF( ALPHA.EQ.ZERO ) THEN
            IF( BETA.EQ.ZERO ) THEN
               DO 340 J = 1, N
                  DO 330 I = 1, M
                     B( I, J ) = ZERO
  330             CONTINUE
  340          CONTINUE
            ELSE
               DO 360 J = 1, N
                  DO 350 I = 1, M
                     B( I, J ) = BETA * B( I, J )
  350             CONTINUE
  360          CONTINUE
            END IF
*
         ELSE IF( ALPHA.EQ.ONE ) THEN
            IF( BETA.EQ.ZERO ) THEN
               DO 380 J = 1, N
                  DO 370 I = 1, M
                     B( I, J ) = A( J, I )
  370             CONTINUE
  380          CONTINUE
            ELSE IF( BETA.EQ.ONE ) THEN
               DO 400 J = 1, N
                  DO 390 I = 1, M
                     B( I, J ) = A( J, I ) + B( I, J )
  390             CONTINUE
  400          CONTINUE
            ELSE
               DO 420 J = 1, N
                  DO 410 I = 1, M
                     B( I, J ) = A( J, I ) + BETA * B( I, J )
  410             CONTINUE
  420          CONTINUE
            END IF
*
         ELSE
            IF( BETA.EQ.ZERO ) THEN
               DO 440 J = 1, N
                  DO 430 I = 1, M
                     B( I, J ) = ALPHA * A( J, I )
  430             CONTINUE
  440          CONTINUE
            ELSE IF( BETA.EQ.ONE ) THEN
               DO 460 J = 1, N
                  DO 450 I = 1, M
                     B( I, J ) = ALPHA * A( J, I ) + B( I, J )
  450             CONTINUE
  460          CONTINUE
            ELSE
               DO 480 J = 1, N
                  DO 470 I = 1, M
                     B( I, J ) = ALPHA * A( J, I ) + BETA * B( I, J )
  470             CONTINUE
  480          CONTINUE
            END IF
         END IF
*
*     Other cases (for genral matrix additions)
*
      ELSE
         IF( ALPHA.EQ.ZERO ) THEN
            IF( BETA.EQ.ZERO ) THEN
               DO 500 J = 1, N
                  DO 490 I = 1, M
                     B( I, J ) = ZERO
  490             CONTINUE
  500          CONTINUE
*
            ELSE
               IF( M.EQ.LDB ) THEN
                  CALL SSCAL( M*N, BETA, B( 1, 1 ), 1 )
               ELSE IF( LSAME( MODE, 'V' ) ) THEN
                  DO 510 J = 1, N
                     CALL SSCAL( M, BETA, B( 1, J ), 1 )
  510             CONTINUE
               ELSE
                  DO 530 J = 1, N
                     DO 520 I = 1, M
                        B( I, J ) = BETA * B( I, J )
  520                CONTINUE
  530             CONTINUE
               END IF
            END IF
*
         ELSE IF( ALPHA.EQ.ONE ) THEN
            IF( BETA.EQ.ZERO ) THEN
               IF( M.EQ.LDA .AND. M.EQ.LDB ) THEN
                  CALL SCOPY( M*N, A( 1, 1 ), 1, B( 1, 1 ), 1 )
               ELSE IF( LSAME( MODE, 'V' ) ) THEN
                  DO 540 J = 1, N
                     CALL SCOPY( M, A( 1, J ), 1, B( 1, J ), 1 )
  540             CONTINUE
               ELSE
                  DO 560 J = 1, N
                     DO 550 I = 1, M
                        B( I, J ) = A( I, J )
  550                CONTINUE
  560             CONTINUE
               END IF
*
            ELSE IF( BETA.EQ.ONE ) THEN
               DO 580 J = 1, N
                  DO 570 I = 1, M
                     B( I, J ) = A( I, J ) + B( I, J )
  570             CONTINUE
  580          CONTINUE
*
            ELSE
               DO 600 J = 1, N
                  DO 590 I = 1, M
                     B( I, J ) = A( I, J ) + BETA * B( I, J )
  590             CONTINUE
  600          CONTINUE
            END IF
*
         ELSE
            IF( BETA.EQ.ZERO ) THEN
               DO 620 J = 1, N
                  DO 610 I = 1, M
                     B( I, J ) = ALPHA * A( I, J )
  610             CONTINUE
  620          CONTINUE
*
            ELSE IF( BETA.EQ.ONE ) THEN
               IF( M.EQ.LDA .AND. M.EQ.LDB ) THEN
                  CALL SAXPY( M*N, ALPHA, A( 1, 1 ), 1, B( 1, 1 ), 1 )
               ELSE IF( LSAME( MODE, 'V' ) ) THEN
                  DO 630 J = 1, N
                     CALL SAXPY( M, ALPHA, A( 1, J ), 1, B( 1, J ), 1 )
  630             CONTINUE
               ELSE
                  DO 650 J = 1, N
                     DO 640 I = 1, M
                        B( I, J ) = ALPHA * A( I, J ) + B( I, J )
  640                CONTINUE
  650             CONTINUE
               END IF
*
            ELSE
               DO 670 J = 1, N
                  DO 660 I = 1, M
                     B( I, J ) = ALPHA * A( I, J ) + BETA * B( I, J )
  660             CONTINUE
  670          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of PBSMATADD
*
      END
