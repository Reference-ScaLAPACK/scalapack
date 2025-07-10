      SUBROUTINE PSGELQRV( M, N, A, IA, JA, DESCA, TAU, WORK )
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
      REAL               A( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSGELQRV computes sub( A ) = A(IA:IA+M-1,JA:JA+N-1) from L, Q
*  computed by PSGELQF.
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
*          by PSGELQF. On exit, the original matrix is restored.
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
*  TAU     (local input) REAL, array, dimension
*          LOCr(IA+MIN(M,N)-1). This array contains the scalar factors
*          TAU of the elementary reflectors computed by PSGELQF. TAU
*          is tied to the distributed matrix A.
*
*  WORK    (local workspace) REAL array, dimension
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
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            I, IACOL, IAROW, IB, ICOFF, ICTXT, IIA, IL, IN,
     $                   IPT, IPV, IPW, J, JJA, JV, K, MYCOL, MYROW,
     $                   NPCOL, NPROW, NQ
*     ..
*     .. Local Arrays ..
      INTEGER            DESCV( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, INFOG2L, PSLACPY,
     $                   PSLARFB, PSLARFT, PSLASET
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
      IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+K-1 )
      IL = MAX( ( (IA+K-2) / DESCA( MB_ ) ) * DESCA( MB_ ) + 1, IA )
*
      ICOFF = MOD( JA-1, DESCA( NB_ ) )
      CALL INFOG2L( IL, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      NQ = NUMROC( N+ICOFF, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IPV = 1
      IPT = IPV + NQ * DESCA( MB_ )
      IPW = IPT + DESCA( MB_ ) * DESCA( MB_ )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ' ' )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', 'D-ring' )
*
      CALL DESCSET( DESCV, DESCA( MB_ ), N + ICOFF, DESCA( MB_ ),
     $              DESCA( NB_ ), IAROW, IACOL, ICTXT, DESCA( MB_ ) )
*
      DO 10 I = IL, IN+1, -DESCA( MB_ )
         IB = MIN( IA+K-I, DESCA( MB_ ) )
         J  = JA + I - IA
         JV = 1  + I - IA + ICOFF
*
*        Compute upper triangular matrix T
*
         CALL PSLARFT( 'Forward', 'Rowwise', N-J+JA, IB, A, I, J, DESCA,
     $                 TAU, WORK( IPT ), WORK( IPW ) )
*
*        Copy Householder vectors into workspace
*
         CALL PSLACPY( 'Upper', IB, N-J+JA, A, I, J, DESCA, WORK( IPV ),
     $                 1, JV, DESCV )
         CALL PSLASET( 'Lower', IB, N-J+JA, ZERO, ONE, WORK( IPV ), 1,
     $                 JV, DESCV )
*
*        Zeroes the strict upper triangular part of sub( A ) to get
*        block column of L
*
         CALL PSLASET( 'Upper', IB, N-J+JA-1, ZERO, ZERO, A, I, J+1,
     $                 DESCA )
*
*        Apply block Householder transformation
*
         CALL PSLARFB( 'Right', 'Transpose', 'Forward', 'Rowwise',
     $                 M-I+IA, N-J+JA, IB, WORK( IPV ), 1, JV, DESCV,
     $                 WORK( IPT ), A, I, J, DESCA, WORK( IPW ) )
*
         DESCV( RSRC_ ) = MOD( DESCV( RSRC_ ) + NPROW - 1, NPROW )
*
   10 CONTINUE
*
*     Handle first block separately
*
      IB = IN - IA + 1
*
*     Compute upper triangular matrix T
*
      CALL PSLARFT( 'Forward', 'Rowwise', N, IB, A, IA, JA, DESCA, TAU,
     $              WORK( IPT ), WORK( IPW ) )
*
*     Copy Householder vectors into workspace
*
      CALL PSLACPY( 'Upper', IB, N, A, IA, JA, DESCA, WORK( IPV ), 1,
     $              ICOFF+1, DESCV )
      CALL PSLASET( 'Lower', IB, N, ZERO, ONE, WORK, 1, ICOFF+1, DESCV )
*
*     Zeroes the strict upper triangular part of sub( A ) to get
*     block column of L
*
      CALL PSLASET( 'Upper', IB, N-1, ZERO, ZERO, A, IA, JA+1, DESCA )
*
*     Apply block Householder transformation
*
      CALL PSLARFB( 'Right', 'Transpose', 'Forward', 'Rowwise', M, N,
     $              IB, WORK( IPV ), 1, ICOFF+1, DESCV, WORK( IPT ), A,
     $              IA, JA, DESCA, WORK( IPW ) )
*
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      RETURN
*
*     End of PSGELQRV
*
      END
