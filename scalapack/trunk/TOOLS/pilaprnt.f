      SUBROUTINE PILAPRNT( M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT,
     $                     CMATNM, NOUT, WORK )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, ICPRNT, IRPRNT, JA, M, N, NOUT
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      CMATNM
      INTEGER            DESCA( * )
      INTEGER            A( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PILAPRNT prints to the standard output a distributed matrix sub( A )
*  denoting A(IA:IA+M-1,JA:JA+N-1). The local pieces are sent and
*  printed by the process of coordinates (IRPRNT, ICPRNT).
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input) @(typec) pointer into the local memory to a
*          local array of dimension (LLD_A, LOCc(JA+N-1) ) containing
*          the local pieces of the distributed matrix sub( A ).
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  IRPRNT  (global input) INTEGER
*          The row index of the printing process.
*
*  ICPRNT  (global input) INTEGER
*          The column index of the printing process.
*
*  CMATNM  (global input) CHARACTER*(*)
*          Identifier of the distributed matrix to be printed.
*
*  NOUT    (global input) INTEGER
*          The unit number for output file. NOUT = 6, ouput to screen,
*          NOUT = 0, output to stderr.
*
*  WORK    (local workspace) @(typec)
*          Working array of minimum size equal to MB_A.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
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
     $                   IGERV2D, IGESD2D
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
         IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
            IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
               DO 10 K = 0, IB-1
                  WRITE( NOUT, FMT = 9999 )
     $                   CMATNM, IA+K, JA+H, A( II+K+(JJ+H-1)*LDA )
   10          CONTINUE
            END IF
         ELSE
            IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
               CALL IGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ), LDA,
     $                       IRPRNT, ICPRNT )
            ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
               CALL IGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                       ICURROW, ICURCOL )
               DO 20 K = 1, IB
                  WRITE( NOUT, FMT = 9999 )
     $                   CMATNM, IA+K-1, JA+H, WORK( K )
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
            IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
               IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  DO 30 K = 0, IB-1
                     WRITE( NOUT, FMT = 9999 )
     $                      CMATNM, I+K, JA+H, A( II+K+(JJ+H-1)*LDA )
   30             CONTINUE
               END IF
            ELSE
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                  CALL IGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                          LDA, IRPRNT, ICPRNT )
               ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  CALL IGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                          ICURROW, ICURCOL )
                  DO 40 K = 1, IB
                     WRITE( NOUT, FMT = 9999 )
     $                      CMATNM, I+K-1, JA+H, WORK( K )
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
            IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
               IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  DO 70 K = 0, IB-1
                     WRITE( NOUT, FMT = 9999 )
     $                      CMATNM, IA+K, J+H, A( II+K+(JJ+H-1)*LDA )
   70             CONTINUE
               END IF
            ELSE
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                  CALL IGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                          LDA, IRPRNT, ICPRNT )
               ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  CALL IGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                          ICURROW, ICURCOL )
                  DO 80 K = 1, IB
                     WRITE( NOUT, FMT = 9999 )
     $                      CMATNM, IA+K-1, J+H, WORK( K )
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
               IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
                  IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                     DO 90 K = 0, IB-1
                        WRITE( NOUT, FMT = 9999 )
     $                         CMATNM, I+K, J+H, A( II+K+(JJ+H-1)*LDA )
   90                CONTINUE
                  END IF
               ELSE
                  IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                     CALL IGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                             LDA, IRPRNT, ICPRNT )
                   ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                     CALL IGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                             ICURROW, ICURCOL )
                     DO 100 K = 1, IB
                        WRITE( NOUT, FMT = 9999 )
     $                         CMATNM, I+K-1, J+H, WORK( K )
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
 9999 FORMAT(A,'(',I6,',',I6,')=',I8)
*
      RETURN
*
*     End of PILAPRNT
*
      END
