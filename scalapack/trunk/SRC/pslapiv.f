      SUBROUTINE PSLAPIV( DIREC, ROWCOL, PIVROC, M, N, A, IA, JA,
     $                    DESCA, IPIV, IP, JP, DESCIP, IWORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*     .. Scalar Arguments ..
      CHARACTER*1        DIREC, PIVROC, ROWCOL
      INTEGER            IA, IP, JA, JP, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCIP( * ), IPIV( * ), IWORK( * )
      REAL               A( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLAPIV applies either P (permutation matrix indicated by IPIV)
*  or inv( P ) to a general M-by-N distributed matrix
*  sub( A ) = A(IA:IA+M-1,JA:JA+N-1), resulting in row or column
*  pivoting. The pivot vector may be distributed across a process row
*  or a column. The pivot vector should be aligned with the distributed
*  matrix A. This routine will transpose the pivot vector if necessary.
*  For example if the row pivots should be applied to the columns of
*  sub( A ), pass ROWCOL='C' and PIVROC='C'.
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
*  Restrictions
*  ============
*
*  IPIV must always be a distributed vector (not a matrix).  Thus:
*  IF( ROWPIV .EQ. 'C' ) THEN
*     JP must be 1
*  ELSE
*     IP must be 1
*  END IF
*
*  The following restrictions apply when IPIV must be transposed:
*  IF( ROWPIV.EQ.'C' .AND. PIVROC.EQ.'C') THEN
*      DESCIP(MB_) must equal DESCA(NB_)
*  ELSE IF( ROWPIV.EQ.'R" .AND. PIVROC.EQ.'R') THEN
*      DESCIP(NB_) must equal DESCA(MB_)
*  END IF
*
*  Arguments
*  =========
*
*  DIREC   (global input) CHARACTER*1
*          Specifies in which order the permutation is applied:
*            = 'F' (Forward) Applies pivots Forward from top of matrix.
*                  Computes P*sub( A ).
*            = 'B' (Backward) Applies pivots Backward from bottom of
*                  matrix. Computes inv( P )*sub( A ).
*
*  ROWCOL  (global input) CHARACTER*1
*          Specifies if the rows or columns are to be permuted:
*             = 'R' Rows will be permuted,
*             = 'C' Columns will be permuted.
*
*  PIVROC  (global input) CHARACTER*1
*          Specifies whether IPIV is distributed over a process row
*          or column:
*          = 'R' IPIV distributed over a process row
*          = 'C' IPIV distributed over a process column
*
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of
*          rows of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) REAL pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          distributed submatrix sub( A ) to which the row or column
*          interchanges will be applied. On exit, the local pieces
*          of the permuted distributed submatrix.
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
*  IPIV    (local input) INTEGER array, dimension (LIPIV) where LIPIV is
*          when ROWCOL='R' or 'r':
*             >= LOCr( IA+M-1 ) + MB_A      if PIVROC='C' or 'c',
*             >= LOCc( M + MOD(JP-1,NB_P) ) if PIVROC='R' or 'r', and,
*          when ROWCOL='C' or 'c':
*             >= LOCr( N + MOD(IP-1,MB_P) ) if PIVROC='C' or 'c',
*             >= LOCc( JA+N-1 ) + NB_A      if PIVROC='R' or 'r'.
*          This array contains the pivoting information. IPIV(i) is the
*          global row (column), local row (column) i was swapped with.
*          When ROWCOL='R' or 'r' and PIVROC='C' or 'c', or ROWCOL='C'
*          or 'c' and PIVROC='R' or 'r', the last piece of this array of
*          size MB_A (resp. NB_A) is used as workspace. In those cases,
*          this array is tied to the distributed matrix A.
*
*  IP      (global input) INTEGER
*          The row index in the global array P indicating the first
*          row of sub( P ).
*
*  JP      (global input) INTEGER
*          The column index in the global array P indicating the
*          first column of sub( P ).
*
*  DESCIP  (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed vector IPIV.
*
*  IWORK   (local workspace) INTEGER array, dimension (LDW)
*          where LDW is equal to the workspace necessary for
*          transposition, and the storage of the tranposed IPIV:
*
*          Let LCM be the least common multiple of NPROW and NPCOL.
*          IF( ROWCOL.EQ.'R' .AND. PIVROC.EQ.'R' ) THEN
*             IF( NPROW.EQ.NPCOL ) THEN
*                LDW = LOCr( N_P + MOD(JP-1, NB_P) ) + NB_P
*             ELSE
*                LDW = LOCr( N_P + MOD(JP-1, NB_P) ) +
*                      NB_P * CEIL( CEIL(LOCc(N_P)/NB_P) / (LCM/NPCOL) )
*             END IF
*          ELSE IF( ROWCOL.EQ.'C' .AND. PIVROC.EQ.'C' ) THEN
*             IF( NPROW.EQ.NPCOL ) THEN
*                LDW = LOCc( M_P + MOD(IP-1, MB_P) ) + MB_P
*             ELSE
*                LDW = LOCc( M_P + MOD(IP-1, MB_P) ) +
*                      MB_P * CEIL( CEIL(LOCr(M_P)/MB_P) / (LCM/NPROW) )
*             END IF
*          ELSE
*             IWORK is not referenced.
*          END IF
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
      LOGICAL            ROWPVT
      INTEGER            I, ICTXT, ICURCOL, ICURROW, IIP, ITMP, IPT,
     $                   JJP, JPT, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            DESCPT( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGEBR2D, IGEBS2D,
     $                   INFOG2L, PICOL2ROW, PIROW2COL, PSLAPV2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC, INDXG2P
      EXTERNAL           LSAME, NUMROC, INDXG2P
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      ROWPVT = LSAME( ROWCOL, 'R' )
*
*     If we're pivoting the rows of sub( A )
*
      IF( ROWPVT ) THEN
         IF( M.LE.1 .OR. N.LT.1 )
     $      RETURN
*
*        If the pivot vector is already distributed correctly
*
         IF( LSAME( PIVROC, 'C' ) ) THEN
            CALL PSLAPV2( DIREC, ROWCOL, M, N, A, IA, JA, DESCA, IPIV,
     $                    IP, JP, DESCIP )
*
*        Otherwise, we must redistribute IPIV to match PSLAPV2
*
         ELSE
*
*           Take IPIV distributed over row 0, and store it in
*           iwork, distributed over column 0
*
            IPT = MOD( JP-1, DESCA(MB_) )
            DESCPT(M_) = M + IPT + NPROW*DESCA(MB_)
            DESCPT(N_) = 1
            DESCPT(MB_) = DESCA(MB_)
            DESCPT(NB_) = 1
            DESCPT(RSRC_) = INDXG2P( IA, DESCA(MB_), IA, DESCA(RSRC_),
     $                               NPROW )
            DESCPT(CSRC_) = MYCOL
            DESCPT(CTXT_) = ICTXT
            DESCPT(LLD_) = NUMROC( DESCPT(M_), DESCPT(MB_), MYROW,
     $                             DESCPT(RSRC_), NPROW )
            ITMP = NUMROC( DESCIP(N_), DESCIP(NB_), MYCOL,
     $                     DESCIP(CSRC_), NPCOL )
            CALL INFOG2L( IP, JP-IPT, DESCIP, NPROW, NPCOL, MYROW,
     $                    MYCOL, IIP, JJP, ICURROW, ICURCOL )
            CALL PIROW2COL( ICTXT, M+IPT, 1, DESCIP(NB_), IPIV(JJP),
     $                      ITMP, IWORK, DESCPT(LLD_), 0, ICURCOL,
     $                      DESCPT(RSRC_),
     $                      MYCOL, IWORK(DESCPT(LLD_)-DESCPT(MB_)+1) )
*
*           Send column-distributed pivots to all columns
*
            ITMP = DESCPT(LLD_) - DESCPT(MB_)
            IF( MYCOL.EQ.0 ) THEN
               CALL IGEBS2D( ICTXT, 'Row', ' ', ITMP, 1, IWORK, ITMP )
            ELSE
               CALL IGEBR2D( ICTXT, 'Row', ' ', ITMP, 1, IWORK, ITMP,
     $                       MYROW, 0 )
            END IF
*
*           Adjust pivots so they are relative to the start of IWORK,
*           not IPIV
*
            IPT = IPT + 1
            DO 10 I = 1, ITMP
               IWORK(I) = IWORK(I) - JP + IPT
   10       CONTINUE
            CALL PSLAPV2( DIREC, ROWCOL, M, N, A, IA, JA, DESCA, IWORK,
     $                    IPT, 1, DESCPT )
         END IF
*
*     Otherwise, we're pivoting the columns of sub( A )
*
      ELSE
         IF( M.LT.1 .OR. N.LE.1 )
     $      RETURN
*
*        If the pivot vector is already distributed correctly
*
         IF( LSAME( PIVROC, 'R' ) ) THEN
            CALL PSLAPV2( DIREC, ROWCOL, M, N, A, IA, JA, DESCA, IPIV,
     $                    IP, JP, DESCIP )
*
*        Otherwise, we must redistribute IPIV to match PSLAPV2
*
         ELSE
*
*           Take IPIV distributed over column 0, and store it in
*           iwork, distributed over row 0
*
            JPT = MOD( IP-1, DESCA(NB_) )
            DESCPT(M_) = 1
            DESCPT(N_) = N + JPT + NPCOL*DESCA(NB_)
            DESCPT(MB_) = 1
            DESCPT(NB_) = DESCA(NB_)
            DESCPT(RSRC_) = MYROW
            DESCPT(CSRC_) = INDXG2P( JA, DESCA(NB_), JA, DESCA(CSRC_),
     $                               NPCOL )
            DESCPT(CTXT_) = ICTXT
            DESCPT(LLD_) = 1
            CALL INFOG2L( IP-JPT, JP, DESCIP, NPROW, NPCOL, MYROW,
     $                    MYCOL, IIP, JJP, ICURROW, ICURCOL )
            ITMP = NUMROC( N+JPT, DESCPT(NB_), MYCOL, DESCPT(CSRC_),
     $                     NPCOL )
            CALL PICOL2ROW( ICTXT, N+JPT, 1, DESCIP(MB_), IPIV(IIP),
     $                      DESCIP(LLD_), IWORK, MAX(1, ITMP), ICURROW,
     $                      0, 0, DESCPT(CSRC_), IWORK(ITMP+1) )
*
*           Send row-distributed pivots to all rows
*
            IF( MYROW.EQ.0 ) THEN
               CALL IGEBS2D( ICTXT, 'Column', ' ', ITMP, 1, IWORK,
     $                       ITMP )
            ELSE
               CALL IGEBR2D( ICTXT, 'Column', ' ', ITMP, 1, IWORK,
     $                       ITMP, 0, MYCOL )
            END IF
*
*           Adjust pivots so they are relative to the start of IWORK,
*           not IPIV
*
            JPT = JPT + 1
            DO 20 I = 1, ITMP
               IWORK(I) = IWORK(I) - IP + JPT
   20       CONTINUE
            CALL PSLAPV2( DIREC, ROWCOL, M, N, A, IA, JA, DESCA, IWORK,
     $                    1, JPT, DESCPT )
         END IF
      END IF
*
      RETURN
*
*     End of PSLAPIV
*
      END
