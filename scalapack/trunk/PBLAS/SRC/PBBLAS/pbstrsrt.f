      SUBROUTINE PBSTRSRT( ICONTXT, ADIST, M, N, NB, A, LDA, BETA, B,
     $                     LDB, LCMP, LCMQ, NINT )
*
*  -- PB-BLAS routine (version 2.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory.
*     April 28, 1996
*
*     .. Scalar Arguments ..
      CHARACTER*1        ADIST
      INTEGER            ICONTXT, LCMP, LCMQ, LDA, LDB, M, N, NB, NINT
      REAL               BETA
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  PBSTRSRT forms   T <== A + beta * T, where T is a sorted
*  condensed block row (or column) from a block column (or row) of A
*  with sorting index ISRT
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Variables ..
      INTEGER            JA, JB, K, KK, NJUMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           PBSMATADD
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL
      EXTERNAL           ICEIL, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( ADIST, 'R' ) ) THEN
         NJUMP = NB * LCMQ
         DO 20 K = 0, LCMQ-1
            JA = NINT * MOD( K*LCMP, LCMQ ) + 1
            JB = K * NB + 1
*
            DO 10 KK = 1, ICEIL( NINT, NB )
               IF( N.LT.JB ) GO TO 20
               CALL PBSMATADD( ICONTXT, 'G', M, MIN( N-JB+1, NB ), ONE,
     $                         A(1, JA), LDA, BETA, B(1, JB), LDB )
               JA = JA + NB
               JB = JB + NJUMP
   10       CONTINUE
   20    CONTINUE
*
*     if( LSAME( ADIST, 'C') ) then
*
      ELSE
         NJUMP = NB * LCMP
         DO 40 K = 0, LCMP-1
            JA = 1
            JB = K * NB + 1
*
            DO 30 KK = 1, ICEIL( NINT, NB )
               IF( M.LT.JB ) GO TO 40
               CALL PBSMATADD( ICONTXT, 'G', MIN( M-JB+1, NB ), N, ONE,
     $                         A(JA, N*MOD(K*LCMQ,LCMP)+1), LDA, BETA,
     $                         B(JB, 1), LDB )
               JA = JA + NB
               JB = JB + NJUMP
   30       CONTINUE
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of PBSTRSRT
*
      END
