      SUBROUTINE PCGERQRV( M, N, A, IA, JA, DESCA, TAU, WORK )
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
      COMPLEX            A( * ),  TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PCGERQRV computes sub( A ) = A(IA:IA+M-1,JA:JA+N-1) from R, Q
*  computed by PCGERQF.
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
*  A       (local input/local output) COMPLEX pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, sub( A ) contains the the factors R and Q computed
*          by PCGERQF. On exit, the original matrix is restored.
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
*  TAU     (local input) COMPLEX, array, dimension LOCr(M_A).
*          This array contains the scalar factors TAU of the elementary
*          reflectors computed by PCGERQF. TAU is tied to the dis-
*          tributed matrix A.
*
*  WORK    (local workspace) COMPLEX array, dimension (LWORK)
*          LWORK = MB_A * ( Mp0 + 2*Nq0 + MB_A ), where
*          Mp0   = NUMROC( M+IROFF, MB_A, MYROW, IAROW, NPROW ) * NB_A,
*          Nq0   = NUMROC( N+ICOFF, NB_A, MYCOL, IACOL, NPCOL ) * MB_A,
*          IROFF = MOD( IA-1, MB_A ), ICOFF = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
*                           NPROW ),
*          IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
*                           NPCOL ),
*          and NUMROC, INDXG2P are ScaLAPACK tool functions.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE  = ( 1.0E+0, 0.0E+0 ),
     $                     ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            I, IACOL, IAROW, IB, ICOFF, ICTXT, IIA, IN,
     $                   IPT, IPV, IPW, JJA, JV, K, MYCOL, MYROW, NPCOL,
     $                   NPROW, NQ
*     ..
*     .. Local Arrays ..
      INTEGER            DESCV( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, INFOG2L, PCLACPY,
     $                   PCLARFB, PCLARFT, PCLASET, PB_TOPGET,
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
      IN = MIN( ICEIL( IA+M-K, DESCA( MB_ ) ) * DESCA( MB_ ), IA+M-1 )
*
      ICOFF = MOD( JA-1, DESCA( NB_ ) )
      CALL INFOG2L( IA+M-K, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $              IIA, JJA, IAROW, IACOL )
      NQ = NUMROC( N+ICOFF, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IPV = 1
      IPT = IPV + NQ * DESCA( MB_ )
      IPW = IPT + DESCA( MB_ ) * DESCA( MB_ )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ' ' )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', 'I-ring' )
*
      CALL DESCSET( DESCV, DESCA( MB_), N + ICOFF, DESCA( MB_ ),
     $              DESCA( NB_ ), IAROW, IACOL, ICTXT, DESCA( MB_ ) )
*
*     Handle first block separately
*
      IB = IN - IA - M + K + 1
      JV = 1 + N - K + ICOFF
*
*     Compute upper triangular matrix T
*
      CALL PCLARFT( 'Backward', 'Rowwise', N-M+IN-IA+1, IB, A, IA+M-K,
     $              JA, DESCA, TAU, WORK( IPT ), WORK( IPW ) )
*
*     Copy Householder vectors into workspace
*
      CALL PCLACPY( 'All', IB, N-M+IN-IA+1, A, IA+M-K, JA, DESCA,
     $              WORK( IPV ), 1, ICOFF+1, DESCV )
      CALL PCLASET( 'Upper', IB, IB, ZERO, ONE, WORK( IPV ), 1, JV,
     $              DESCV )
*
*     Zeoes the strict lower triangular part of sub( A ) to get block
*     column of R
*
      CALL PCLASET( 'All', IB, N-K, ZERO, ZERO, A, IA+M-K, JA,
     $              DESCA )
      CALL PCLASET( 'Lower', IB-1, IB, ZERO, ZERO, A, IA+M-K+1,
     $              JA+N-K, DESCA )
*
*     Apply block Householder transformation
*
      CALL PCLARFB( 'Right', 'Conjugate transpose', 'Backward',
     $              'Rowwise', IN-IA+1, N-M+IN-IA+1, IB, WORK( IPV ), 1,
     $              ICOFF+1, DESCV, WORK( IPT ), A, IA, JA, DESCA,
     $              WORK( IPW ) )
*
      DESCV( RSRC_ ) = MOD( DESCV( RSRC_ ) + 1, NPROW )
*
*     Loop over the remaining row blocks
*
      DO 10 I = IN+1, IA+M-1, DESCA( MB_ )
         IB = MIN( IA+M-I, DESCA( MB_ ) )
         JV = 1 + N - M + I - IA + ICOFF
*
*        Compute upper triangular matrix T
*
         CALL PCLARFT( 'Backward', 'Rowwise', N-M+I+IB-IA, IB, A, I, JA,
     $                 DESCA, TAU, WORK( IPT ), WORK( IPW ) )
*
*        Copy Householder vectors into workspace
*
         CALL PCLACPY( 'All', IB, N-M+I+IB-IA, A, I, JA, DESCA,
     $                 WORK( IPV ), 1, ICOFF+1, DESCV )
         CALL PCLASET( 'Upper', IB, IB, ZERO, ONE, WORK( IPV ), 1, JV,
     $                 DESCV )
*
*        Zeoes the strict Lower triangular part of sub( A ) to get
*        block column of R
*
         CALL PCLASET( 'All', IB, N-M+I-IA, ZERO, ZERO, A, I, JA,
     $                 DESCA )
         CALL PCLASET( 'Lower', IB-1, IB, ZERO, ZERO, A, I+1,
     $                 JA+N-M+I-IA, DESCA )
*
*        Apply block Householder transformation
*
         CALL PCLARFB( 'Right', 'Conjugate transpose', 'Backward',
     $                 'Rowwise', I+IB-IA, N-M+I+IB-IA, IB, WORK( IPV ),
     $                 1, ICOFF+1, DESCV, WORK( IPT ), A, IA, JA, DESCA,
     $                 WORK( IPW ) )
*
         DESCV( RSRC_ ) = MOD( DESCV( RSRC_ ) + 1, NPROW )
*
   10 CONTINUE
*
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      RETURN
*
*     End of PCGERQRV
*
      END
