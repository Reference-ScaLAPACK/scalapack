      SUBROUTINE PBSTRGET( ICONTXT, ADIST, M, N, MNB, A, LDA, MCROW,
     $                     MCCOL, IGD, MYROW, MYCOL, NPROW, NPCOL )
*
*  -- PB-BLAS routine (version 2.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory.
*     April 28, 1996
*
*     .. Scalar Arguments ..
      CHARACTER*1        ADIST
      INTEGER            ICONTXT, IGD, LDA, M, MCCOL, MCROW, MNB, MYCOL,
     $                   MYROW, N, NPCOL, NPROW
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  PBSTRGET forms a row block of A from scattered row subblocks if
*  ADIST = 'R', or forms a column block of A from scattered column
*  subblocks,  if ADIST = 'C'.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, TWO
      PARAMETER          ( ONE = 1.0E+0, TWO = 2.0E+0 )
*     ..
*     .. Local Variables ..
      INTEGER            KINT, KINT2, KLEN, KMOD, KPPOS, NLEN, NNUM,
     $                   NTLEN
      REAL               TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, NUMROC
      EXTERNAL           LSAME,  ICEIL, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGERV2D, SGESD2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*
*     if A is a row block, it needs to communicate columnwise.
*
      IF( LSAME( ADIST, 'R' ) ) THEN
         KPPOS = MOD( NPROW+MYROW-MCROW, NPROW )
         IF( MOD( KPPOS, IGD ).EQ.0 ) THEN
            KINT = IGD
            NLEN = N
            NNUM = MIN( NPROW/IGD, MNB-MCCOL )
            TEMP = REAL( NNUM )
            NTLEN = N * NNUM
            NNUM = IGD * NNUM
            IF( KPPOS.GE.NNUM ) GO TO 30
            KPPOS = MOD( KPPOS, NPROW )
*
   10       CONTINUE
            IF( TEMP.GT.ONE ) THEN
               KINT2 = 2 * KINT
               KMOD = MOD( KPPOS, KINT2 )
*
               IF( KMOD.EQ.0 ) THEN
                  IF( KPPOS+KINT.LT.NNUM ) THEN
                     KLEN = NTLEN - (KPPOS/KINT2)*(KINT2/IGD)*N
                     KLEN = MIN( KLEN-NLEN, NLEN )
                     CALL SGERV2D( ICONTXT, M, KLEN, A(1,NLEN+1), LDA,
     $                             MOD(MYROW+KINT, NPROW), MYCOL )
                     NLEN = NLEN + KLEN
                  END IF
               ELSE
                  CALL SGESD2D( ICONTXT, M, NLEN, A, LDA,
     $                          MOD(NPROW+MYROW-KINT, NPROW), MYCOL )
                  GO TO 30
               END IF
*
               KINT = KINT2
               TEMP = TEMP / TWO
               GO TO 10
            END IF
         END IF
*
*     if A is a column block, it needs to communicate rowwise.
*
      ELSE IF( LSAME( ADIST, 'C' ) ) THEN
*
         KPPOS = MOD( NPCOL+MYCOL-MCCOL, NPCOL )
         IF( MOD( KPPOS, IGD ).EQ.0 ) THEN
            KINT = IGD
            NLEN = N
            NNUM = MIN( NPCOL/IGD, MNB-MCROW )
            TEMP = REAL( NNUM )
            NTLEN = N * NNUM
            NNUM = IGD * NNUM
            IF( KPPOS.GE.NNUM ) GO TO 30
            KPPOS = MOD( KPPOS, NPCOL )
*
   20       CONTINUE
            IF( TEMP.GT.ONE ) THEN
               KINT2 = 2 * KINT
               KMOD = MOD( KPPOS, KINT2 )
*
               IF( KMOD.EQ.0 ) THEN
                  IF( KPPOS+KINT.LT.NNUM ) THEN
                     KLEN = NTLEN - (KPPOS/KINT2)*(KINT2/IGD)*N
                     KLEN = MIN( KLEN-NLEN, NLEN )
                     CALL SGERV2D( ICONTXT, M, KLEN, A(1,NLEN+1), LDA,
     $                             MYROW, MOD(MYCOL+KINT, NPCOL) )
                     NLEN = NLEN + KLEN
                  END IF
               ELSE
                  CALL SGESD2D( ICONTXT, M, NLEN, A, LDA, MYROW,
     $                          MOD(NPCOL+MYCOL-KINT, NPCOL) )
                  GO TO 30
               END IF
*
               KINT = KINT2
               TEMP = TEMP / TWO
               GO TO 20
            END IF
         END IF
      END IF
*
   30 CONTINUE
*
      RETURN
*
*     End of PBSTRGET
*
      END
