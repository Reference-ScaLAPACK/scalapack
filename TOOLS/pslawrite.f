      SUBROUTINE PSLAWRITE( FILNAM, M, N, A, IA, JA, DESCA, IRWRIT,
     $                      ICWRIT, WORK )
*
*  -- ScaLAPACK tools routine (version 1.8) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     written by Antoine Petitet, August 1995 (petitet@cs.utk.edu)
*     adapted by Julie Langou, April 2007 (julie@cs.utk.edu)
*
*     .. Scalar Arguments ..
      INTEGER            IA, ICWRIT, IRWRIT, JA, M, N
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      FILNAM
      INTEGER            DESCA( * )
      REAL               A( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLAWRITE writes to a file named FILNAMa distributed matrix sub( A )
*  denoting A(IA:IA+M-1,JA:JA+N-1). The local pieces are sent to and
*  written by the process of coordinates (IRWWRITE, ICWRIT).
*
*  WORK must be of size >= MB_ = DESCA( MB_ ).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NOUT
      PARAMETER          ( NOUT = 13 )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            H, I, IACOL, IAROW, IB, ICTXT, ICURCOL,
     $                   ICURROW, II, IIA, IN, J, JB, JJ, JJA, JN, K,
     $                   LDA, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_GRIDINFO, INFOG2L,
     $                   SGERV2D, SGESD2D
*     ..
*     .. External Functions ..
      INTEGER            ICEIL
      EXTERNAL           ICEIL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( MYROW.EQ.IRWRIT .AND. MYCOL.EQ.ICWRIT ) THEN
         OPEN( NOUT, FILE=FILNAM, STATUS='UNKNOWN' )
         WRITE( NOUT, FMT = * ) M, N
      END IF
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $              IIA, JJA, IAROW, IACOL )
      ICURROW = IAROW
      ICURCOL = IACOL
      II = IIA
      JJ = JJA
      LDA = DESCA( LLD_ )
*
*     Handle the first block of column separately
*
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      JB = JN-JA+1
      DO 60 H = 0, JB-1
         IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+M-1 )
         IB = IN-IA+1
         IF( ICURROW.EQ.IRWRIT .AND. ICURCOL.EQ.ICWRIT ) THEN
            IF( MYROW.EQ.IRWRIT .AND. MYCOL.EQ.ICWRIT ) THEN
               DO 10 K = 0, IB-1
                  WRITE( NOUT, FMT = 9999 ) A( II+K+(JJ+H-1)*LDA )
   10          CONTINUE
            END IF
         ELSE
            IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
               CALL SGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ), LDA,
     $                       IRWRIT, ICWRIT )
            ELSE IF( MYROW.EQ.IRWRIT .AND. MYCOL.EQ.ICWRIT ) THEN
               CALL SGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                       ICURROW, ICURCOL )
               DO 20 K = 1, IB
                  WRITE( NOUT, FMT = 9999 ) WORK( K )
   20          CONTINUE
            END IF
         END IF
         IF( MYROW.EQ.ICURROW )
     $      II = II + IB
         ICURROW = MOD( ICURROW+1, NPROW )
         CALL BLACS_BARRIER( ICTXT, 'All' )
*
*        Loop over remaining block of rows
*
         DO 50 I = IN+1, IA+M-1, DESCA( MB_ )
            IB = MIN( DESCA( MB_ ), IA+M-I )
            IF( ICURROW.EQ.IRWRIT .AND. ICURCOL.EQ.ICWRIT ) THEN
               IF( MYROW.EQ.IRWRIT .AND. MYCOL.EQ.ICWRIT ) THEN
                  DO 30 K = 0, IB-1
                     WRITE( NOUT, FMT = 9999 ) A( II+K+(JJ+H-1)*LDA )
   30             CONTINUE
               END IF
            ELSE
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                  CALL SGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                          LDA, IRWRIT, ICWRIT )
               ELSE IF( MYROW.EQ.IRWRIT .AND. MYCOL.EQ.ICWRIT ) THEN
                  CALL SGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                          ICURROW, ICURCOL )
                  DO 40 K = 1, IB
                     WRITE( NOUT, FMT = 9999 ) WORK( K )
   40             CONTINUE
               END IF
            END IF
            IF( MYROW.EQ.ICURROW )
     $         II = II + IB
            ICURROW = MOD( ICURROW+1, NPROW )
            CALL BLACS_BARRIER( ICTXT, 'All' )
   50    CONTINUE
*
        II = IIA
        ICURROW = IAROW
   60 CONTINUE
*
      IF( MYCOL.EQ.ICURCOL )
     $   JJ = JJ + JB
      ICURCOL = MOD( ICURCOL+1, NPCOL )
      CALL BLACS_BARRIER( ICTXT, 'All' )
*
*     Loop over remaining column blocks
*
      DO 130 J = JN+1, JA+N-1, DESCA( NB_ )
         JB = MIN(  DESCA( NB_ ), JA+N-J )
         DO 120 H = 0, JB-1
            IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+M-1 )
            IB = IN-IA+1
            IF( ICURROW.EQ.IRWRIT .AND. ICURCOL.EQ.ICWRIT ) THEN
               IF( MYROW.EQ.IRWRIT .AND. MYCOL.EQ.ICWRIT ) THEN
                  DO 70 K = 0, IB-1
                     WRITE( NOUT, FMT = 9999 ) A( II+K+(JJ+H-1)*LDA )
   70             CONTINUE
               END IF
            ELSE
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                  CALL SGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                          LDA, IRWRIT, ICWRIT )
               ELSE IF( MYROW.EQ.IRWRIT .AND. MYCOL.EQ.ICWRIT ) THEN
                  CALL SGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                          ICURROW, ICURCOL )
                  DO 80 K = 1, IB
                     WRITE( NOUT, FMT = 9999 ) WORK( K )
   80             CONTINUE
               END IF
            END IF
            IF( MYROW.EQ.ICURROW )
     $         II = II + IB
            ICURROW = MOD( ICURROW+1, NPROW )
            CALL BLACS_BARRIER( ICTXT, 'All' )
*
*           Loop over remaining block of rows
*
            DO 110 I = IN+1, IA+M-1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IA+M-I )
               IF( ICURROW.EQ.IRWRIT .AND. ICURCOL.EQ.ICWRIT ) THEN
                  IF( MYROW.EQ.IRWRIT .AND. MYCOL.EQ.ICWRIT ) THEN
                     DO 90 K = 0, IB-1
                        WRITE( NOUT, FMT = 9999 ) A( II+K+(JJ+H-1)*LDA )
   90                CONTINUE
                  END IF
               ELSE
                  IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                     CALL SGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                             LDA, IRWRIT, ICWRIT )
                   ELSE IF( MYROW.EQ.IRWRIT .AND. MYCOL.EQ.ICWRIT ) THEN
                     CALL SGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                             ICURROW, ICURCOL )
                     DO 100 K = 1, IB
                        WRITE( NOUT, FMT = 9999 ) WORK( K )
  100                CONTINUE
                  END IF
               END IF
               IF( MYROW.EQ.ICURROW )
     $            II = II + IB
               ICURROW = MOD( ICURROW+1, NPROW )
               CALL BLACS_BARRIER( ICTXT, 'All' )
  110       CONTINUE
*
            II = IIA
            ICURROW = IAROW
  120    CONTINUE
*
         IF( MYCOL.EQ.ICURCOL )
     $      JJ = JJ + JB
         ICURCOL = MOD( ICURCOL+1, NPCOL )
         CALL BLACS_BARRIER( ICTXT, 'All' )
*
  130 CONTINUE
*
      IF( MYROW.EQ.IRWRIT .AND. MYCOL.EQ.ICWRIT ) THEN
         CLOSE( NOUT )
      END IF
*
 9999 FORMAT( E15.8 )
*
      RETURN
*
*     End of PSLAWRITE
*
      END
