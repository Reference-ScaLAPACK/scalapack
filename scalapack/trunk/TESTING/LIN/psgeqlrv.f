      SUBROUTINE PSGEQLRV( M, N, A, IA, JA, DESCA, TAU, WORK )
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
      INTEGER            DESCA( * )
      REAL               A( * ),  TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSGEQLRV computes sub( A ) = A(IA:IA+M-1,JA:JA+N-1) from L, Q
*  computed by PSGEQLF.
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
*  A       (local input/local output) REAL pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, sub( A ) contains the the factors L and Q computed
*          by PSGEQLF. On exit, the original matrix is restored.
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
*  TAU     (local input) REAL, array, dimension LOCc(N_A).
*          This array contains the scalar factors TAU of the elementary
*          reflectors computed by PSGEQLF. TAU is tied to the dis-
*          tributed matrix A.
*
*  WORK    (local workspace) REAL array, dimension (LWORK)
*          LWORK = NB_A * ( 2*Mp0 + Nq0 + NB_A ), where
*          Mp0   = NUMROC( M+IROFF, MB_A, MYROW, IAROW, NPROW ) * NB_A,
*          Nq0   = NUMROC( N+ICOFF, NB_A, MYCOL, IACOL, NPCOL ) * MB_A,
*          IROFF = MOD( IA-1, MB_A ), ICOFF = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
*                           NPROW ),
*          IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
*                           NPCOL ),
*          and NUMROC, INDXG2P are ScaLAPACK tool functions;
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
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            IACOL, IAROW, ICTXT, IIA, IPT, IPV, IPW, IROFF,
     $                   IV, J, JB, JJA, JN, K, MP, MYCOL, MYROW, NPCOL,
     $                   NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            DESCV( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, INFOG2L, PSLACPY,
     $                   PSLARFB, PSLARFT, PSLASET, PB_TOPGET,
     $                   PB_TOPSET
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, NUMROC
      EXTERNAL           ICEIL, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      K = MIN( M, N )
      JN = MIN( ICEIL( JA+N-K, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
*
      IROFF = MOD( IA-1, DESCA( MB_ ) )
      CALL INFOG2L( IA, JA+N-K, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $              IIA, JJA, IAROW, IACOL )
      MP = NUMROC( M+IROFF, DESCA( MB_ ), MYROW, IAROW, NPROW )
      IPV = 1
      IPT = IPV + MP * DESCA( NB_ )
      IPW = IPT + DESCA( NB_ ) * DESCA( NB_ )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'I-ring' )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
*
      CALL DESCSET( DESCV, M+IROFF, DESCA( NB_ ), DESCA( MB_ ),
     $              DESCA( NB_ ), IAROW, IACOL, ICTXT, MAX( 1, MP ) )
*
*     Handle first block separately
*
      IV = 1 + M - K + IROFF
      JB = JN - JA - N + K + 1
*
*     Compute upper triangular matrix T
*
      CALL PSLARFT( 'Backward', 'Columnwise', M-N+JN-JA+1, JB, A, IA,
     $              JA+N-K, DESCA, TAU, WORK( IPT ), WORK( IPW ) )
*
*     Copy Householder vectors into workspace
*
      CALL PSLACPY( 'All', M-N+JN-JA+1, JB, A, IA, JA+N-K, DESCA,
     $              WORK( IPV ), IROFF+1, 1, DESCV )
      CALL PSLASET( 'Lower', JB, JB, ZERO, ONE, WORK( IPV ), IV,
     $              1, DESCV )
*
*     Zeoes the strict upper triangular part of A to get block
*     row of L
*
      CALL PSLASET( 'All', M-K, JB, ZERO, ZERO, A, IA, JA+N-K,
     $              DESCA )
      CALL PSLASET( 'Upper', JB, JB-1, ZERO, ZERO, A, IA+M-K,
     $              JA+N-K+1, DESCA )
*
*     Apply block Householder transformation
*
      CALL PSLARFB( 'Left', 'No transpose', 'Backward', 'Columnwise',
     $              M-N+JN-JA+1, JN-JA+1, JB, WORK( IPV ), IROFF+1, 1,
     $              DESCV, WORK( IPT ), A, IA, JA, DESCA, WORK( IPW ) )
*
      DESCV( CSRC_ ) = MOD( DESCV( CSRC_ ) + 1, NPCOL )
*
*     Loop over the remaining column blocks
*
      DO 10 J = JN+1, JA+N-1, DESCA( NB_ )
         JB = MIN( JA+N-J, DESCA( NB_ ) )
         IV = 1 + M - N + J - JA + IROFF
*
*        Compute upper triangular matrix T
*
         CALL PSLARFT( 'Backward', 'Columnwise', M-N+J+JB-JA, JB, A, IA,
     $                 J, DESCA, TAU, WORK( IPT ), WORK( IPW ) )
*
*        Copy Householder vectors into workspace
*
         CALL PSLACPY( 'All', M-N+J+JB-JA, JB, A, IA, J, DESCA,
     $                 WORK( IPV ), IROFF+1, 1, DESCV )
         CALL PSLASET( 'Lower', JB, JB, ZERO, ONE, WORK( IPV ), IV,
     $                 1, DESCV )
*
*        Zeoes the strict upper triangular part of sub( A ) to get
*        block row of L
*
         CALL PSLASET( 'All', M-N+J-JA, JB, ZERO, ZERO, A, IA, J,
     $                 DESCA )
         CALL PSLASET( 'Upper', JB, JB-1, ZERO, ZERO, A, IA+M-N+J-JA,
     $                 J+1, DESCA )
*
*        Apply block Householder transformation
*
         CALL PSLARFB( 'Left', 'No transpose', 'Backward', 'Columnwise',
     $                 M-N+J+JB-JA, J+JB-JA, JB, WORK( IPV ), IROFF+1,
     $                 1, DESCV, WORK( IPT ), A, IA, JA, DESCA,
     $                 WORK( IPW ) )
*
         DESCV( CSRC_ ) = MOD( DESCV( CSRC_ ) + 1, NPCOL )
*
   10 CONTINUE
*
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      RETURN
*
*     End of PSGEQLRV
*
      END
