      SUBROUTINE SMATADD( M, N, ALPHA, A, LDA, BETA, C, LDC )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDC, M, N
      REAL               ALPHA, BETA
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  SMATADD performs the following local matrix-matrix operation
*
*               C := alpha * A + beta * C,
*
*  where alpha and beta are scalars, and A and C are m by n arrays.
*
*  Arguments
*  =========
*
*  M       (local input) INTEGER
*          The number of rows of the array A. M >= 0.
*
*  N       (local input) INTEGER
*          The number of columns of the array A. N >= 0.
*
*  ALPHA   (local input) REAL
*          The scalar ALPHA.
*
*  A       (local input) REAL
*          Array, dimension (LDA,*), the array A.
*
*  LDA     (local input) INTEGER
*          The leading dimension of the array A, LDA >= MAX(1, M)
*
*  BETA    (local input) REAL
*          The scalar BETA.
*
*  C       (local input/local output) REAL
*          Array, dimension (LDC,*), the array C.
*
*  LDC     (local input) INTEGER
*          The leading dimension of the array C, LDC >= MAX(1, M)
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE  = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      IF( (M.EQ.0).OR.(N.EQ.0).OR.((ALPHA.EQ.ZERO).AND.(BETA.EQ.ONE)) )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( BETA.EQ.ZERO ) THEN
            IF( ALPHA.EQ.ZERO ) THEN
               DO 10 I = 1, M
                  C( I, 1 ) = ZERO
   10          CONTINUE
            ELSE
               DO 20  I = 1, M
                  C( I, 1 ) = ALPHA*A( I, 1 )
   20          CONTINUE
            END IF
         ELSE
            IF( ALPHA.EQ.ONE ) THEN
               IF( BETA.EQ.ONE ) THEN
                  DO 30 I = 1, M
                     C( I, 1 ) = A( I, 1 ) + C( I, 1 )
   30             CONTINUE
               ELSE
                  DO 40 I = 1, M
                     C( I, 1 ) = A( I, 1 ) + BETA*C( I, 1 )
   40             CONTINUE
               END IF
            ELSE IF( BETA.EQ.ONE ) THEN
               DO 50 I = 1, M
                  C( I, 1 ) = ALPHA*A( I, 1 ) + C( I, 1 )
   50          CONTINUE
            ELSE
               DO 60 I = 1, M
                  C( I, 1 ) = ALPHA*A( I, 1 ) + BETA*C( I, 1 )
   60          CONTINUE
            END IF
         END IF
      ELSE
         IF( BETA.EQ.ZERO ) THEN
            IF( ALPHA.EQ.ZERO ) THEN
               DO 80 J = 1, N
                  DO 70 I = 1, M
                     C( I, J ) = ZERO
   70             CONTINUE
   80          CONTINUE
            ELSE
               DO 100 J = 1, N
                  DO 90 I = 1, M
                     C( I, J ) = ALPHA * A( I, J )
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
            IF( ALPHA.EQ.ONE ) THEN
               IF( BETA.EQ.ONE ) THEN
                  DO 120 J = 1, N
                     DO 110 I = 1, M
                        C( I, J ) = A( I, J ) + C( I, J )
  110                CONTINUE
  120             CONTINUE
               ELSE
                  DO 140 J = 1, N
                     DO 130 I = 1, M
                        C( I, J ) = A( I, J ) + BETA * C( I, J )
  130                CONTINUE
  140             CONTINUE
               END IF
            ELSE IF( BETA.EQ.ONE ) THEN
               DO 160 J = 1, N
                  DO 150 I = 1, M
                     C( I, J ) = C( I, J ) + ALPHA * A( I, J )
  150             CONTINUE
  160          CONTINUE
            ELSE
               DO 180 J = 1, N
                  DO 170 I = 1, M
                     C( I, J ) = ALPHA * A( I, J ) + BETA * C( I, J )
  170             CONTINUE
  180          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of SMATADD
*
      END
