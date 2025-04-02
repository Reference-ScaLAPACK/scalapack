      SUBROUTINE PZTZRZRV( M, N, A, IA, JA, DESCA, TAU, WORK )
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
      COMPLEX*16         A( * ),  TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZTZRZRV computes sub( A ) = A(IA:IA+M-1,JA:JA+N-1) from T, Z
*  computed by PZTZRZF.
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
*          columns of the distributed submatrix sub( A ). N >= M >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, sub( A ) contains the the factors T and Z computed
*          by PZTZRZF. On exit, the original matrix is restored.
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
*  TAU     (local input) COMPLEX*16, array, dimension LOCr(M_A).
*          This array contains the scalar factors TAU of the elementary
*          reflectors computed by PZTZRZF. TAU is tied to the dis-
*          tributed matrix A.
*
*  WORK    (local workspace) COMPLEX*16 array, dimension (LWORK)
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
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            I, IACOL, IAROW, IB, ICOFF, ICTXT, IIA, IN,
     $                   IPT, IPV, IPW, JJA, JM1, JV, L, MYCOL, MYROW,
     $                   NPCOL, NPROW, NQ
*     ..
*     .. Local Arrays ..
      INTEGER            DESCV( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, INFOG2L, PB_TOPGET,
     $                   PB_TOPSET, PZLACPY, PZLARZB, PZLARZT,
     $                   PZLASET
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
*     Quick return if possible
*
      IF( N.LT.M )
     $   RETURN
*
      L = N - M
      JM1 = JA + MIN( M+1, N ) - 1
      IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+M-1 )
      ICOFF = MOD( JA-1, DESCA( NB_ ) )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      NQ = NUMROC( N+ICOFF, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IPV = 1
      IPT = IPV + NQ * DESCA( MB_ )
      IPW = IPT + DESCA( MB_ ) * DESCA( MB_ )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ' ' )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', 'I-ring' )
*
      CALL DESCSET( DESCV, DESCA( MB_ ), N + ICOFF, DESCA( MB_ ),
     $              DESCA( NB_ ), IAROW, IACOL, ICTXT, DESCA( MB_ ) )
*
*     Handle first block separately
*
      IB = IN - IA + 1
      JV = ICOFF + JM1 - JA + 1
*
*     Compute upper triangular matrix T
*
      CALL PZLARZT( 'Backward', 'Rowwise', L, IB, A, IA, JM1, DESCA,
     $              TAU, WORK( IPT ), WORK( IPW ) )
*
*     Copy Householder vectors into workspace
*
      CALL PZLACPY( 'All', IB, L, A, IA, JM1, DESCA, WORK( IPV ), 1,
     $              JV, DESCV )
*
*     Save temporarily strict lower part of A(IA:IA+IB-1,JA:JA+IB-1)
*
      CALL PZLACPY( 'Lower', IB-1, IB-1, A, IA+1, JA, DESCA,
     $              WORK( IPV ), 1, ICOFF+1, DESCV )
*
*     Zeroes the row panel of sub( A ) to get T(IA:IN,JA:JA+N-1)
*
      CALL PZLASET( 'All', IB, L, ZERO, ZERO, A, IA, JM1, DESCA )
      CALL PZLASET( 'Lower', IB-1, IB-1, ZERO, ZERO, A, IA+1, JA,
     $              DESCA )
*
*     Apply block Householder transformation
*
      CALL PZLARZB( 'Right', 'Conjugate transpose', 'Backward',
     $              'Rowwise', IN-IA+1, N, IB, L, WORK( IPV ), 1, JV,
     $              DESCV, WORK( IPT ), A, IA, JA, DESCA, WORK( IPW ) )
*
*     Restore strict lower part of A( IA:IA+IB-1, JA:JA+N-1 )
*
      CALL PZLACPY( 'Lower', IB-1, IB-1, WORK( IPV ), 1, ICOFF+1, DESCV,
     $              A, IA+1, JA, DESCA )
*
      DESCV( RSRC_ ) = MOD( DESCV( RSRC_ ) + 1, NPROW )
*
*     Loop over the remaining row blocks
*
      DO 10 I = IN+1, IA+M-1, DESCA( MB_ )
         IB = MIN( IA+M-I, DESCA( MB_ ) )
*
*        Compute upper triangular matrix T
*
         CALL PZLARZT( 'Backward', 'Rowwise', L, IB, A, I, JM1, DESCA,
     $                 TAU, WORK( IPT ), WORK( IPW ) )
*
*        Copy Householder vectors into workspace
*
         CALL PZLACPY( 'All', IB, L, A, I, JM1, DESCA, WORK( IPV ), 1,
     $                 JV, DESCV )
*
*        Save temporarily strict lower part of A(I:I+IB-1,J:J+IB-1 )
*
         CALL PZLACPY( 'Lower', IB-1, IB-1, A, I+1, JA+I-IA, DESCA,
     $                 WORK( IPV ), 1, ICOFF+1+I-IA, DESCV )
*
*        Zeoes the row panel of sub( A ) to get T(IA:I-1,JA+I-IA:JA+N-1)
*
         CALL PZLASET( 'All', IB, L, ZERO, ZERO, A, I, JM1, DESCA )
         CALL PZLASET( 'Lower', IB-1, IB-1, ZERO, ZERO, A, I+1, JA+I-IA,
     $                 DESCA )
*
*        Apply block Householder transformation
*
         CALL PZLARZB( 'Right', 'Conjugate transpose', 'Backward',
     $                 'Rowwise', I+IB-IA, N-I+IA, IB, L, WORK( IPV ),
     $                 1, JV, DESCV, WORK( IPT ), A, IA, JA+I-IA, DESCA,
     $                 WORK( IPW ) )
*
         CALL PZLACPY( 'Lower', IB-1, IB-1, WORK( IPV ), 1,
     $                 ICOFF+1+I-IA, DESCV, A, I+1, JA+I-IA, DESCA )
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
*     End of PZTZRZRV
*
      END
