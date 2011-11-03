      SUBROUTINE BSLAAPP( ISIDE, M, N, NB, A, LDA, NITRAF, ITRAF,
     $                    DTRAF, WORK )
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            ISIDE, LDA, M, N, NB, NITRAF
*     ..
*     .. Array Arguments ..
      INTEGER            ITRAF( * )
      REAL               A( LDA, * ), DTRAF( * ), WORK( * )
*
*
*  Purpose
*  =======
*
*  BSLAAPP computes
*
*          B = Q**T * A       or       B = A * Q,
*
*  where A is an M-by-N matrix and Q is an orthogonal matrix represented
*  by the parameters in the arrays ITRAF and DTRAF as described in
*  BSTREXC.
*
*  This is an auxiliary routine called by BDTRSEN.
*
*  Arguments
*  =========
*
*  ISIDE   (input) INTEGER
*          Specifies whether Q multiplies A from the left or right as
*          follows:
*          = 0: compute B = Q**T * A;
*          = 1: compute B = A * Q.
*
*  M       (input) INTEGER
*          The number of rows of A.
*
*  N       (input) INTEGER
*          The number of columns of A.
*
*  NB      (input) INTEGER
*          If ISIDE = 0, the Q is applied block column-wise to the rows
*          of A and NB specifies the maximal width of the block columns.
*          If ISIDE = 1, this variable is not referenced.
*
*  A       (input/output) REAL             array, dimension (LDA,N)
*          On entry, the matrix A.
*          On exit, A is overwritten by B.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  NITRAF  (input) INTEGER
*          Length of the array ITRAF. NITRAF >= 0.
*
*  ITRAF   (input) INTEGER array, length NITRAF
*          List of parameters for representing the transformation
*          matrix Q, see BSTREXC.
*
*  DTRAF   (output) REAL             array, length k, where
*          List of parameters for representing the transformation
*          matrix Q, see BSTREXC.
*
*  WORK    (workspace) REAL             array, dimension (N)
*
*  =====================================================================
*

*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IT, J, NNB, PD
      REAL               TAU
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARFX, SROT
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( ISIDE.EQ.0 ) THEN
*
*        Apply Q from left.
*
         DO 20 J = 1, N, NB
            PD = 1
            NNB = MIN( NB, N - J + 1 )
            DO 10 I = 1, NITRAF
               IT = ITRAF(I)
               IF( IT.LE.M ) THEN
*
*                 Apply Givens rotation.
*
                  CALL SROT( NNB, A(IT,J), LDA, A(IT+1,J), LDA,
     $                       DTRAF(PD), DTRAF(PD+1) )
                  PD = PD + 2
               ELSE IF( IT.LE.2*M ) THEN
*
*                 Apply Householder reflector of first kind.
*
                  TAU = DTRAF(PD)
                  DTRAF(PD) = ONE
                  CALL SLARFX( 'Left', 3, NNB, DTRAF(PD), TAU,
     $                         A(IT-M,J), LDA, WORK )
                  DTRAF(PD) = TAU
                  PD = PD + 3
               ELSE
*
*                 Apply Householder reflector of second kind.
*
                  TAU = DTRAF(PD+2)
                  DTRAF(PD+2) = ONE
                  CALL SLARFX( 'Left', 3, NNB, DTRAF(PD), TAU,
     $                         A(IT-2*M,J), LDA, WORK )
                  DTRAF(PD+2) = TAU
                  PD = PD + 3
               END IF
   10       CONTINUE
   20    CONTINUE
      ELSE
         PD = 1
         DO 30 I = 1, NITRAF
            IT = ITRAF(I)
            IF( IT.LE.N ) THEN
*
*              Apply Givens rotation.
*
               CALL SROT( M, A(1,IT), 1, A(1,IT+1), 1, DTRAF(PD),
     $                    DTRAF(PD+1) )
               PD = PD + 2
            ELSE IF( IT.LE.2*N ) THEN
*
*              Apply Householder reflector of first kind.
*
               TAU = DTRAF(PD)
               DTRAF(PD) = ONE
               CALL SLARFX( 'Right', M, 3, DTRAF(PD), TAU, A(1,IT-N),
     $                      LDA, WORK )
               DTRAF(PD) = TAU
               PD = PD + 3
            ELSE
*
*              Apply Householder reflector of second kind.
*
               TAU = DTRAF(PD+2)
               DTRAF(PD+2) = ONE
               CALL SLARFX( 'Right', M, 3, DTRAF(PD), TAU, A(1,IT-2*N),
     $                      LDA, WORK )
               DTRAF(PD+2) = TAU
               PD = PD + 3
            END IF
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of BSLAAPP
*
      END
