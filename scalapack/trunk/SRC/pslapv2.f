      SUBROUTINE PSLAPV2( DIREC, ROWCOL, M, N, A, IA, JA, DESCA, IPIV,
     $                    IP, JP, DESCIP )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          DIREC, ROWCOL
      INTEGER            IA, IP, JA, JP, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCIP( * ), IPIV( * )
      REAL               A( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLAPV2 applies either P (permutation matrix indicated by IPIV)
*  or inv( P ) to a M-by-N distributed matrix sub( A ) denoting
*  A(IA:IA+M-1,JA:JA+N-1), resulting in row or column pivoting.  The
*  pivot vector should be aligned with the distributed matrix A.  For
*  pivoting the rows of sub( A ), IPIV should be distributed along a
*  process column and replicated over all process rows.  Similarly,
*  IPIV should be distributed along a process row and replicated over
*  all process columns for column pivoting.
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
*  DIREC   (global input) CHARACTER
*          Specifies in which order the permutation is applied:
*            = 'F' (Forward) Applies pivots Forward from top of matrix.
*                  Computes P * sub( A );
*            = 'B' (Backward) Applies pivots Backward from bottom of
*                  matrix. Computes inv( P ) * sub( A ).
*
*  ROWCOL  (global input) CHARACTER
*          Specifies if the rows or columns are to be permuted:
*            = 'R' Rows will be permuted,
*            = 'C' Columns will be permuted.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) REAL pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this local array contains the local pieces of the
*          distributed matrix sub( A ) to which the row or columns
*          interchanges will be applied. On exit, this array contains
*          the local pieces of the permuted distributed matrix.
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
*  IPIV    (input) INTEGER array, dimension >= LOCr(M_A)+MB_A if
*          ROWCOL = 'R', LOCc(N_A)+NB_A otherwise. It contains
*          the pivoting information. IPIV(i) is the global row (column),
*          local row (column) i was swapped with.  The last piece of the
*          array of size MB_A (resp. NB_A) is used as workspace. IPIV is
*          tied to the distributed matrix A.
*
*  IP      (global input) INTEGER
*          IPIV's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JP      (global input) INTEGER
*          IPIV's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCIP  (global and local input) INTEGER array of dimension 8
*          The array descriptor for the distributed matrix IPIV.
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
      LOGICAL            FORWRD, ROWPVT
      INTEGER            I, IB, ICTXT, ICURCOL, ICURROW, IIP, IP1, ITMP,
     $                   IPVWRK, J, JB, JJP, JP1, K, MA, MBA, MYCOL,
     $                   MYROW, NBA, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGEBS2D, IGEBR2D, INFOG2L,
     $                   PSSWAP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, NUMROC
      EXTERNAL           ICEIL, LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
      ROWPVT = LSAME( ROWCOL, 'R' )
      IF( ROWPVT ) THEN
         IF( M.LE.1 .OR. N.LT.1 )
     $      RETURN
      ELSE
         IF( M.LT.1 .OR. N.LE.1 )
     $      RETURN
      END IF
      FORWRD = LSAME( DIREC, 'F' )
*
*
*     Get grid and matrix parameters
*
      MA    = DESCA( M_ )
      MBA   = DESCA( MB_ )
      NBA   = DESCA( NB_ )
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     If I'm applying pivots from beginning to end (e.g., repeating
*     pivoting done earlier).  Thus this section computes P * sub( A ).
*
      IF( FORWRD ) THEN
         CALL INFOG2L( IP, JP, DESCIP, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIP, JJP, ICURROW, ICURCOL )
*
*        If I'm pivoting the rows of sub( A )
*
         IF( ROWPVT ) THEN
            IPVWRK = NUMROC( DESCIP( M_ ), DESCIP( MB_ ), MYROW,
     $                       DESCIP( RSRC_ ), NPROW ) + 1 -
     $                       DESCIP( MB_ )
*
*          Loop over rows of sub( A )
*
            I = IA
            IB = MIN( M, ICEIL( IA, MBA ) * MBA - IA + 1 )
   10       CONTINUE
*
*              Find local pointer into IPIV, and broadcast this block's
*              pivot information to everyone in process column
*
               IF( MYROW.EQ.ICURROW ) THEN
                  CALL IGEBS2D( ICTXT, 'Columnwise', ' ', IB, 1,
     $                          IPIV( IIP ), IB )
                  ITMP = IIP
                  IIP = IIP + IB
               ELSE
                  ITMP = IPVWRK
                  CALL IGEBR2D( ICTXT, 'Columnwise', ' ', IB, 1,
     $                          IPIV( ITMP ), IB, ICURROW, MYCOL )
               END IF
*
*              Pivot the block of rows
*
               DO 20 K = I, I+IB-1
                  IP1 = IPIV( ITMP ) - IP + IA
                  IF( IP1.NE.K )
     $               CALL PSSWAP( N, A, K, JA, DESCA, MA, A, IP1, JA,
     $                            DESCA, MA )
                  ITMP = ITMP + 1
   20          CONTINUE
*
*              Go on to next row of processes, increment row counter,
*              and figure number of rows to pivot next
*
               ICURROW = MOD( ICURROW+1, NPROW )
               I = I + IB
               IB = MIN( MBA, M-I+IA )
            IF( IB .GT. 0 ) GOTO 10
*
*        If I am pivoting the columns of sub( A )
*
         ELSE
            IPVWRK = NUMROC( DESCIP( N_ ), DESCIP( NB_ ), MYCOL,
     $                       DESCIP( CSRC_ ), NPCOL ) + 1 -
     $                       DESCIP( NB_ )
*
*          Loop over columns of sub( A )
*
            J = JA
            JB = MIN( N, ICEIL( JA, NBA ) * NBA - JA + 1 )
   30       CONTINUE
*
*              Find local pointer into IPIV, and broadcast this block's
*              pivot information to everyone in process row
*
               IF( MYCOL.EQ.ICURCOL ) THEN
                  CALL IGEBS2D( ICTXT, 'Rowwise', ' ', JB, 1,
     $                          IPIV( JJP ), JB )
                  ITMP = JJP
                  JJP = JJP + JB
               ELSE
                  ITMP = IPVWRK
                  CALL IGEBR2D( ICTXT, 'Rowwise', ' ', JB, 1,
     $                          IPIV( ITMP ), JB, MYROW, ICURCOL )
               END IF
*
*              Pivot the block of columns
*
               DO 40 K = J, J+JB-1
                  JP1 = IPIV( ITMP ) - JP + JA
                  IF( JP1.NE.K )
     $               CALL PSSWAP( M, A, IA, K, DESCA, 1, A, IA, JP1,
     $                            DESCA, 1 )
                  ITMP = ITMP + 1
   40          CONTINUE
*
*              Go on to next column of processes, increment column
*              counter, and figure number of columns to pivot next
*
               ICURCOL = MOD( ICURCOL+1, NPCOL )
               J = J + JB
               JB = MIN( NBA, N-J+JA )
            IF( JB .GT. 0 ) GOTO 30
         END IF
*
*     If I want to apply pivots in reverse order, i.e. reversing
*     pivoting done earlier.  Thus this section computes
*     inv( P ) * sub( A ).
*
      ELSE
*
*        If I'm pivoting the rows of sub( A )
*
         IF( ROWPVT ) THEN
            CALL INFOG2L( IP+M-1, JP, DESCIP, NPROW, NPCOL, MYROW,
     $                    MYCOL, IIP, JJP, ICURROW, ICURCOL )
*
            IPVWRK = NUMROC( DESCIP( M_ ), DESCIP( MB_ ), MYROW,
     $                       DESCIP( RSRC_ ), NPROW ) + 1 -
     $                       DESCIP( MB_ )
*
*           If I'm not in the current process row, my IIP points out
*           past end of pivot vector (since I don't own a piece of the
*           last row). Adjust IIP so it points at last pivot entry.
*
            IF( MYROW.NE.ICURROW ) IIP = IIP - 1
*
*           Loop over rows in reverse order, starting at last row
*
            I = IA + M - 1
            IB = MOD( I, MBA )
            IF( IB .EQ. 0 ) IB = MBA
            IB = MIN( IB, M )
   50       CONTINUE
*
*              Find local pointer into IPIV, and broadcast this block's
*              pivot information to everyone in process column
*
               IF( MYROW.EQ.ICURROW ) THEN
                  ITMP = IIP
                  IIP = IIP - IB
                  CALL IGEBS2D( ICTXT, 'Columnwise', ' ', IB, 1,
     $                          IPIV( IIP+1 ), IB )
               ELSE
                  CALL IGEBR2D( ICTXT, 'Columnwise', ' ', IB, 1,
     $                          IPIV( IPVWRK ), IB, ICURROW, MYCOL )
                  ITMP = IPVWRK + IB - 1
               END IF
*
*              Pivot the block of rows
*
               DO 60 K = I, I-IB+1, -1
                  IP1 = IPIV( ITMP ) - IP + IA
                  IF( IP1.NE.K )
     $               CALL PSSWAP( N, A, K, JA, DESCA, MA, A, IP1, JA,
     $                            DESCA, MA )
                  ITMP = ITMP - 1
   60          CONTINUE
*
*              Go to previous row of processes, decrement row counter,
*              and figure number of rows to be pivoted next
*
               ICURROW = MOD( NPROW+ICURROW-1, NPROW )
               I = I - IB
               IB = MIN( MBA, I-IA+1 )
            IF( IB .GT. 0 ) GOTO 50
*
*        Otherwise, I'm pivoting the columns of sub( A )
*
         ELSE
            CALL INFOG2L( IP, JP+N-1, DESCIP, NPROW, NPCOL, MYROW,
     $                    MYCOL, IIP, JJP, ICURROW, ICURCOL )
            IPVWRK = NUMROC( DESCIP( N_ ), DESCIP( NB_ ), MYCOL,
     $                       DESCIP( CSRC_ ), NPCOL ) + 1 -
     $                       DESCIP( NB_ )
*
*           If I'm not in the current process column, my JJP points out
*           past end of pivot vector (since I don't own a piece of the
*           last column). Adjust JJP so it points at last pivot entry.
*
            IF( MYCOL.NE.ICURCOL ) JJP = JJP - 1
*
*          Loop over columns in reverse order starting at last column
*
            J = JA + N - 1
            JB = MOD( J, NBA )
            IF( JB .EQ. 0 ) JB = NBA
            JB = MIN( JB, N )
   70       CONTINUE
*
*              Find local pointer into IPIV, and broadcast this block's
*              pivot information to everyone in process row
*
               IF( MYCOL.EQ.ICURCOL ) THEN
                  ITMP = JJP
                  JJP = JJP - JB
                  CALL IGEBS2D( ICTXT, 'Rowwise', ' ', JB, 1,
     $                          IPIV( JJP+1 ), JB )
               ELSE
                  CALL IGEBR2D( ICTXT, 'Rowwise', ' ', JB, 1,
     $                          IPIV( IPVWRK ), JB, MYROW, ICURCOL )
                  ITMP = IPVWRK + JB - 1
               END IF
*
*              Pivot a block of columns
*
               DO 80 K = J, J-JB+1, -1
                  JP1 = IPIV( ITMP ) - JP + JA
                  IF( JP1.NE.K )
     $               CALL PSSWAP( M, A, IA, K, DESCA, 1, A, IA, JP1,
     $                            DESCA, 1 )
                  ITMP = ITMP - 1
   80          CONTINUE
*
*              Go to previous row of processes, decrement row counter,
*              and figure number of rows to be pivoted next
*
               ICURCOL = MOD( NPCOL+ICURCOL-1, NPCOL )
               J = J - JB
               JB = MIN( NBA, J-JA+1 )
            IF( JB .GT. 0 ) GOTO 70
         END IF
*
      END IF
*
      RETURN
*
*     End PSLAPV2
*
      END
