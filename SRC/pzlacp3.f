      SUBROUTINE PZLACP3( M, I, A, DESCA, B, LDB, II, JJ, REV )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     July 31, 2001
*
*     .. Scalar Arguments ..
      INTEGER            I, II, JJ, LDB, M, REV
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         A( * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  PZLACP3 is an auxiliary routine that copies from a global parallel
*    array into a local replicated array or vise versa.  Notice that
*    the entire submatrix that is copied gets placed on one node or
*    more.  The receiving node can be specified precisely, or all nodes
*    can receive, or just one row or column of nodes.
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
*          M is the order of the square submatrix that is copied.
*          M >= 0.
*          Unchanged on exit
*
*  I       (global input) INTEGER
*          A(I,I) is the global location that the copying starts from.
*          Unchanged on exit.
*
*  A       (global input/output) COMPLEX*16 array, dimension
*          (DESCA(LLD_),*)
*          On entry, the parallel matrix to be copied into or from.
*          On exit, if REV=1, the copied data.
*          Unchanged on exit if REV=0.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  B       (local input/output) COMPLEX*16 array of size (LDB,M)
*          If REV=0, this is the global portion of the array
*             A(I:I+M-1,I:I+M-1).
*          If REV=1, this is the unchanged on exit.
*
*  LDB     (local input) INTEGER
*          The leading dimension of B.
*
*  II      (global input) INTEGER
*          By using REV 0 & 1, data can be sent out and returned again.
*          If REV=0, then II is destination row index for the node(s)
*             receiving the replicated B.
*             If II>=0,JJ>=0, then node (II,JJ) receives the data
*             If II=-1,JJ>=0, then all rows in column JJ receive the
*                             data
*             If II>=0,JJ=-1, then all cols in row II receive the data
*             If II=-1,JJ=-1, then all nodes receive the data
*          If REV<>0, then II is the source row index for the node(s)
*             sending the replicated B.
*
*  JJ      (global input) INTEGER
*          Similar description as II above
*
*  REV     (global input) INTEGER
*          Use REV = 0 to send global A into locally replicated B
*             (on node (II,JJ)).
*          Use REV <> 0 to send locally replicated B from node (II,JJ)
*             to its owner (which changes depending on its location in
*             A) into the global A.
*
*  Further Details
*  ===============
*
*  Implemented by:  M. Fahey, May 28, 1999
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            COL, CONTXT, HBL, ICOL1, ICOL2, IDI, IDJ, IFIN,
     $                   III, IROW1, IROW2, ISTOP, ISTOPI, ISTOPJ, ITMP,
     $                   JJJ, LDA, MYCOL, MYROW, NPCOL, NPROW, ROW
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG1L, ZGEBR2D, ZGEBS2D,
     $                   ZGERV2D, ZGESD2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
      IF( M.LE.0 )
     $   RETURN
*
      HBL = DESCA( MB_ )
      CONTXT = DESCA( CTXT_ )
      LDA = DESCA( LLD_ )
*
      CALL BLACS_GRIDINFO( CONTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( REV.EQ.0 ) THEN
         DO 20 IDI = 1, M
            DO 10 IDJ = 1, M
               B( IDI, IDJ ) = ZERO
   10       CONTINUE
   20    CONTINUE
      END IF
*
      IFIN = I + M - 1
*
      IF( MOD( I+HBL, HBL ).NE.0 ) THEN
         ISTOP = MIN( I+HBL-MOD( I+HBL, HBL ), IFIN )
      ELSE
         ISTOP = I
      END IF
      IDJ = I
      ISTOPJ = ISTOP
      IF( IDJ.LE.IFIN ) THEN
   30    CONTINUE
         IDI = I
         ISTOPI = ISTOP
         IF( IDI.LE.IFIN ) THEN
   40       CONTINUE
            ROW = MOD( ( IDI-1 ) / HBL, NPROW )
            COL = MOD( ( IDJ-1 ) / HBL, NPCOL )
            CALL INFOG1L( IDI, HBL, NPROW, ROW, 0, IROW1, ITMP )
            IROW2 = NUMROC( ISTOPI, HBL, ROW, 0, NPROW )
            CALL INFOG1L( IDJ, HBL, NPCOL, COL, 0, ICOL1, ITMP )
            ICOL2 = NUMROC( ISTOPJ, HBL, COL, 0, NPCOL )
            IF( ( MYROW.EQ.ROW ) .AND. ( MYCOL.EQ.COL ) ) THEN
               IF( ( II.EQ.-1 ) .AND. ( JJ.EQ.-1 ) ) THEN
*
*                 Send the message to everyone
*
                  IF( REV.EQ.0 ) THEN
                     CALL ZGEBS2D( CONTXT, 'All', ' ', IROW2-IROW1+1,
     $                             ICOL2-ICOL1+1, A( ( ICOL1-1 )*LDA+
     $                             IROW1 ), LDA )
                  END IF
               END IF
               IF( ( II.EQ.-1 ) .AND. ( JJ.NE.-1 ) ) THEN
*
*                 Send the message to Column MYCOL which better be JJ
*
                  IF( REV.EQ.0 ) THEN
                     CALL ZGEBS2D( CONTXT, 'Col', ' ', IROW2-IROW1+1,
     $                             ICOL2-ICOL1+1, A( ( ICOL1-1 )*LDA+
     $                             IROW1 ), LDA )
                  END IF
               END IF
               IF( ( II.NE.-1 ) .AND. ( JJ.EQ.-1 ) ) THEN
*
*                 Send the message to Row MYROW which better be II
*
                  IF( REV.EQ.0 ) THEN
                     CALL ZGEBS2D( CONTXT, 'Row', ' ', IROW2-IROW1+1,
     $                             ICOL2-ICOL1+1, A( ( ICOL1-1 )*LDA+
     $                             IROW1 ), LDA )
                  END IF
               END IF
               IF( ( II.NE.-1 ) .AND. ( JJ.NE.-1 ) .AND.
     $             ( ( MYROW.NE.II ) .OR. ( MYCOL.NE.JJ ) ) ) THEN
*
*                 Recv/Send the message to (II,JJ)
*
                  IF( REV.EQ.0 ) THEN
                     CALL ZGESD2D( CONTXT, IROW2-IROW1+1, ICOL2-ICOL1+1,
     $                             A( ( ICOL1-1 )*LDA+IROW1 ), LDA, II,
     $                             JJ )
                  ELSE
                     CALL ZGERV2D( CONTXT, IROW2-IROW1+1, ICOL2-ICOL1+1,
     $                             B( IDI-I+1, IDJ-I+1 ), LDB, II, JJ )
                  END IF
               END IF
               IF( REV.EQ.0 ) THEN
                  DO 60 JJJ = ICOL1, ICOL2
                     DO 50 III = IROW1, IROW2
                        B( IDI+III-IROW1+1-I, IDJ+JJJ-ICOL1+1-I )
     $                     = A( ( JJJ-1 )*LDA+III )
   50                CONTINUE
   60             CONTINUE
               ELSE
                  DO 80 JJJ = ICOL1, ICOL2
                     DO 70 III = IROW1, IROW2
                        A( ( JJJ-1 )*LDA+III ) = B( IDI+III-IROW1+1-I,
     $                     IDJ+JJJ-ICOL1+1-I )
   70                CONTINUE
   80             CONTINUE
               END IF
            ELSE
               IF( ( II.EQ.-1 ) .AND. ( JJ.EQ.-1 ) ) THEN
                  IF( REV.EQ.0 ) THEN
                     CALL ZGEBR2D( CONTXT, 'All', ' ', IROW2-IROW1+1,
     $                             ICOL2-ICOL1+1, B( IDI-I+1, IDJ-I+1 ),
     $                             LDB, ROW, COL )
                  END IF
               END IF
               IF( ( II.EQ.-1 ) .AND. ( JJ.EQ.MYCOL ) ) THEN
                  IF( REV.EQ.0 ) THEN
                     CALL ZGEBR2D( CONTXT, 'Col', ' ', IROW2-IROW1+1,
     $                             ICOL2-ICOL1+1, B( IDI-I+1, IDJ-I+1 ),
     $                             LDB, ROW, COL )
                  END IF
               END IF
               IF( ( II.EQ.MYROW ) .AND. ( JJ.EQ.-1 ) ) THEN
                  IF( REV.EQ.0 ) THEN
                     CALL ZGEBR2D( CONTXT, 'Row', ' ', IROW2-IROW1+1,
     $                             ICOL2-ICOL1+1, B( IDI-I+1, IDJ-I+1 ),
     $                             LDB, ROW, COL )
                  END IF
               END IF
               IF( ( II.EQ.MYROW ) .AND. ( JJ.EQ.MYCOL ) ) THEN
                  IF( REV.EQ.0 ) THEN
                     CALL ZGERV2D( CONTXT, IROW2-IROW1+1, ICOL2-ICOL1+1,
     $                             B( IDI-I+1, IDJ-I+1 ), LDB, ROW,
     $                             COL )
                  ELSE
                     CALL ZGESD2D( CONTXT, IROW2-IROW1+1, ICOL2-ICOL1+1,
     $                             B( IDI-I+1, IDJ-I+1 ), LDB, ROW,
     $                             COL )
*                    CALL ZGESD2D(CONTXT, IROW2-IROW1+1, ICOL2-ICOL1+1,
*    $                            A((ICOL1-1)*LDA+IROW1),LDA, ROW, COL)
                  END IF
               END IF
            END IF
            IDI = ISTOPI + 1
            ISTOPI = MIN( ISTOPI+HBL, IFIN )
            IF( IDI.LE.IFIN )
     $         GO TO 40
         END IF
         IDJ = ISTOPJ + 1
         ISTOPJ = MIN( ISTOPJ+HBL, IFIN )
         IF( IDJ.LE.IFIN )
     $      GO TO 30
      END IF
      RETURN
*
*     End of PZLACP3
*
      END
