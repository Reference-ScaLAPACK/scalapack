      SUBROUTINE PDGETRRV( M, N, A, IA, JA, DESCA, IPIV, WORK )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 28, 2001
*
*     .. Scalar Arguments ..
      INTEGER            IA, JA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), IPIV( * )
      DOUBLE PRECISION   A( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDGETRRV reforms sub( A ) = A(IA:IA+M-1,JA:JA+N-1) from the
*  triangular matrices L and U returned by PDGETRF.  It multiplies
*  an upper triangular matrix stored in the upper triangle of sub( A )
*  times the unit lower triangular matrix stored in the lower triangle.
*  To accomplish this, the routine basically performs the PDGETRF
*  routine in reverse.
*
*  It computes L*U first, and then apply P: P*L*U => sub( A ). In the
*  J-th loop, the block column (or column panel), which has the lower
*  triangular unit matrix L is multiplied with the block row (or row
*  panel), which contains the upper triangular matrix U.
*
*                     ( L1 )             ( 0  0 )   ( L1*U1   L1*U2   )
*   A` = L * U + A` = (    ) * (U1 U2) + (      ) = (                 )
*                     ( L2 )             ( 0  A`)   ( L2*U1  L2*U2+A` )
*
*  where L1 is a lower unit triangular matrix and U1 is an upper
*  triangular matrix.
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
*          The number of rows to be operated on, i.e. the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, the local pieces of the distributed matrix sub( A )
*          contains the the factors L and U from the factorization
*          sub( A ) = P*L*U; the unit diagonal elements of L are not
*          stored. On exit, the original distributed matrix sub( A )
*          is restored.
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
*  IPIV    (local input) INTEGER array, dimension ( LOCr(M_A)+MB_A )
*          This array contains the pivoting information.
*          IPIV(i) -> The global row local row i was swapped with.
*          This array is tied to the distributed matrix A.
*
*  WORK    (local workspace) DOUBLE PRECISION array of dimension (LWORK)
*          LWORK >= MpA0 * NB_A + NqA0 * MB_A, where
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          MpA0 = NUMROC( M+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*          NqA0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),
*
*          WORK is used to store a block of columns of L, and a block of
*          rows of U. INDXG2P and NUMROC are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            IACOL, IAROW, ICTXT, IL, IPL, IPU, IROFF, J,
     $                   JB, JL, JN, MN, MP, MYCOL, MYROW, NPCOL, NPROW
*     .. Local Arrays ..
      INTEGER            DESCIP( DLEN_ ), DESCL( DLEN_ ),
     $                   DESCU( DLEN_ ), IDUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, PDGEMM, PDLACPY,
     $                   PDLAPIV, PDLASET, PB_TOPGET, PB_TOPSET
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, INDXG2P, NUMROC
      EXTERNAL           ICEIL, INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters.
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IROFF = MOD( IA-1, DESCA( MB_ ) )
      IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ), NPROW )
      MP = NUMROC( M+IROFF, DESCA( MB_ ), MYROW, IAROW, NPROW )
      IPL = 1
      IPU = IPL + MP * DESCA( NB_ )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'S-ring' )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
*
*     Define array descriptors for L and U
*
      MN = MIN( M, N )
      IL = MAX( ( ( IA+MN-2 ) / DESCA( MB_ ) ) * DESCA( MB_ ) + 1, IA )
      JL = MAX( ( ( JA+MN-2 ) / DESCA( NB_ ) ) * DESCA( NB_ ) + 1, JA )
      JN = MIN( ICEIL( JA, DESCA( NB_ ) )*DESCA( NB_ ), JA+MN-1 )
      IAROW = INDXG2P( IL, DESCA( MB_ ), MYROW, DESCA( RSRC_ ), NPROW )
      IACOL = INDXG2P( JL, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ), NPCOL )
*
      CALL DESCSET( DESCL, IA+M-IL, DESCA( NB_ ), DESCA( MB_ ),
     $              DESCA( NB_ ), IAROW, IACOL, ICTXT, MAX( 1, MP ) )
*
      CALL DESCSET( DESCU, DESCA( MB_ ), JA+N-JL, DESCA( MB_ ),
     $              DESCA( NB_ ), IAROW, IACOL, ICTXT, DESCA( MB_ ) )
*
      CALL DESCSET( DESCIP, DESCA( M_ ) + DESCA( MB_ )*NPROW, 1,
     $              DESCA( MB_ ), 1, DESCA( RSRC_ ), MYCOL, ICTXT,
     $              NUMROC( DESCA( M_ ), DESCA( MB_ ), MYROW,
     $                      DESCA( RSRC_ ), NPROW ) + DESCA( MB_ ) )
*
*
      DO 10 J = JL, JN+1, -DESCA( NB_ )
*
         JB = MIN( JA+MN-J, DESCA( NB_ ) )
*
*        Copy unit lower triangular part of sub( A ) into WORK
*
         CALL PDLACPY( 'Lower', M-IL+IA, JB, A, IL, J, DESCA,
     $                 WORK( IPL ), 1, 1, DESCL )
         CALL PDLASET( 'Upper', M-IL+IA, JB, ZERO, ONE, WORK( IPL ),
     $                 1, 1, DESCL )
*
*        Copy upper triangular part of sub( A ) into WORK(IPU)
*
         CALL PDLACPY( 'Upper', JB, JA+N-J, A, IL, J, DESCA,
     $                 WORK( IPU ), 1, 1, DESCU )
         CALL PDLASET( 'Lower', JB-1, JA+N-J, ZERO, ZERO,
     $                 WORK( IPU ), 2, 1, DESCU )
*
*        Zero the strict lower triangular piece of the current block.
*
         CALL PDLASET( 'Lower', IA+M-IL-1, JB, ZERO, ZERO, A, IL+1, J,
     $                 DESCA )
*
*        Zero the upper triangular piece of the current block.
*
         CALL PDLASET( 'Upper', JB, JA+N-J, ZERO, ZERO, A, IL, J,
     $                 DESCA )
*
*        Update the matrix sub( A ).
*
         CALL PDGEMM( 'No transpose', 'No transpose', IA+M-IL,
     $                JA+N-J, JB, ONE, WORK( IPL ), 1, 1, DESCL,
     $                WORK( IPU ), 1, 1, DESCU, ONE, A, IL, J, DESCA )
*
         IL = IL - DESCA( MB_ )
         DESCL( M_ ) = DESCL( M_ ) + DESCL( MB_ )
         DESCL( RSRC_ ) = MOD( DESCL( RSRC_ ) + NPROW - 1, NPROW )
         DESCL( CSRC_ ) = MOD( DESCL( CSRC_ ) + NPCOL - 1, NPCOL )
         DESCU( N_ ) = DESCU( N_ ) + DESCU( NB_ )
         DESCU( RSRC_ ) = DESCL( RSRC_ )
         DESCU( CSRC_ ) = DESCL( CSRC_ )
*
   10 CONTINUE
*
*     Handle first block separately
*
      JB = MIN( JN-JA+1, DESCA( NB_ ) )
*
*     Copy unit lower triangular part of sub( A ) into WORK
*
      CALL PDLACPY( 'Lower', M, JB, A, IA, JA, DESCA, WORK( IPL ),
     $              1, 1, DESCL )
      CALL PDLASET( 'Upper', M, JB, ZERO, ONE, WORK( IPL ), 1, 1,
     $              DESCL )
*
*     Copy upper triangular part of sub( A ) into WORK(IPU)
*
      CALL PDLACPY( 'Upper', JB, N, A, IA, JA, DESCA, WORK( IPU ), 1,
     $              1, DESCU )
      CALL PDLASET( 'Lower', JB-1, N, ZERO, ZERO, WORK( IPU ), 2, 1,
     $              DESCU )
*
*     Zero the strict lower triangular piece of the current block.
*
      CALL PDLASET( 'Lower', M-1, JB, ZERO, ZERO, A, IA+1, JA, DESCA )
*
*     Zero the upper triangular piece of the current block.
*
      CALL PDLASET( 'Upper', JB, N, ZERO, ZERO, A, IA, JA, DESCA )
*
*     Update the matrix sub( A ).
*
      CALL PDGEMM( 'No transpose', 'No transpose', M, N, JB, ONE,
     $             WORK( IPL ), 1, 1, DESCL, WORK( IPU ), 1, 1,
     $             DESCU, ONE, A, IA, JA, DESCA )
*
*     Apply pivots so that sub( A ) = P*L*U
*
      CALL PDLAPIV( 'Backward', 'Row', 'Col', MIN( M, N ), N, A, IA, JA,
     $              DESCA, IPIV, IA, 1, DESCIP, IDUM )
*
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      RETURN
*
*     End of PDGETRRV
*
      END
