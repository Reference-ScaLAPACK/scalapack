      SUBROUTINE PDGEHDRV( N, ILO, IHI, A, IA, JA, DESCA, TAU, WORK )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 28, 2001
*
*     .. Scalar Arguments ..
      INTEGER              IA, IHI, ILO, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER              DESCA( * )
      DOUBLE PRECISION     A( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDGEHDRV computes sub( A ) = A(IA:IA+N-1,JA:JA+N-1) from the
*  orthogonal matrix Q, the Hessenberg matrix, and the array TAU
*  returned by PDGEHRD:
*                       sub( A ) := Q * H * Q'
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  ILO     (global input) INTEGER
*  IHI     (global input) INTEGER
*          It is assumed that sub( A ) is already upper triangular in
*          rows and columns 1:ILO-1 and IHI+1:N. If N > 0,
*          1 <= ILO <= IHI <= N; otherwise set ILO = 1, IHI = N.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the N-by-N
*          general distributed matrix sub( A ) reduced to Hessenberg
*          form by PDGEHRD. The upper triangle and the first sub-
*          diagonal of sub( A ) contain the upper Hessenberg matrix H,
*          and the elements below the first subdiagonal, with the array
*          TAU, represent the orthogonal matrix Q as a product of
*          elementary reflectors. On exit, the original distributed
*          N-by-N matrix sub( A ) is recovered.
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
*  TAU     (local input) DOUBLE PRECISION array, dimension LOCc(JA+N-2)
*          The scalar factors of the elementary reflectors returned by
*          PDGEHRD. TAU is tied to the distributed matrix A.
*
*  WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK).
*          LWORK >= NB*NB + NB*IHLP + MAX[ NB*( IHLP+INLQ ),
*                   NB*( IHLQ + MAX[ IHIP,
*                   IHLP+NUMROC( NUMROC( IHI-ILO+LOFF+1, NB, 0, 0,
*                   NPCOL ), NB, 0, 0, LCMQ ) ] ) ]
*
*          where NB = MB_A = NB_A,
*          LCM is the least common multiple of NPROW and NPCOL,
*          LCM = ILCM( NPROW, NPCOL ), LCMQ = LCM / NPCOL,
*
*          IROFFA = MOD( IA-1, NB ),
*          IAROW = INDXG2P( IA, NB, MYROW, RSRC_A, NPROW ),
*          IHIP = NUMROC( IHI+IROFFA, NB, MYROW, IAROW, NPROW ),
*
*          ILROW = INDXG2P( IA+ILO-1, NB, MYROW, RSRC_A, NPROW ),
*          ILCOL = INDXG2P( JA+ILO-1, NB, MYCOL, CSRC_A, NPCOL ),
*          IHLP = NUMROC( IHI-ILO+IROFFA+1, NB, MYROW, ILROW, NPROW ),
*          IHLQ = NUMROC( IHI-ILO+IROFFA+1, NB, MYCOL, ILCOL, NPCOL ),
*          INLQ = NUMROC( N-ILO+IROFFA+1, NB, MYCOL, ILCOL, NPCOL ).
*
*          ILCM, INDXG2P and NUMROC are ScaLAPACK tool functions;
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
      DOUBLE PRECISION     ZERO
      PARAMETER            ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER              I, IACOL, IAROW, ICTXT, IHLP, II, IOFF, IPT,
     $                     IPV, IPW, IV, J, JB, JJ, JL, K, MYCOL, MYROW,
     $                     NB, NPCOL, NPROW
*     ..
*     .. Local Arrays ..
      INTEGER              DESCV( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER              INDXG2P, NUMROC
      EXTERNAL             INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL             BLACS_GRIDINFO, DESCSET, INFOG2L, PDLARFB,
     $                     PDLARFT, PDLACPY, PDLASET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC            MAX, MIN, MOD
*     ..
*     .. Executable statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Quick return if possible
*
      IF( IHI-ILO.LE.0 )
     $   RETURN
*
      NB = DESCA( MB_ )
      IOFF = MOD( IA+ILO-2, NB )
      CALL INFOG2L( IA+ILO-1, JA+ILO-1, DESCA, NPROW, NPCOL, MYROW,
     $              MYCOL, II, JJ, IAROW, IACOL )
      IHLP = NUMROC( IHI-ILO+IOFF+1, NB, MYROW, IAROW, NPROW )
*
      IPT = 1
      IPV = IPT + NB * NB
      IPW = IPV + IHLP * NB
      JL = MAX( ( ( JA+IHI-2 ) / NB ) * NB + 1, JA + ILO - 1 )
      CALL DESCSET( DESCV, IHI-ILO+IOFF+1, NB, NB, NB, IAROW,
     $              INDXG2P( JL, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $              NPCOL ), ICTXT, MAX( 1, IHLP ) )
*
      DO 10 J = JL, ILO+JA+NB-IOFF-1, -NB
         JB = MIN( JA+IHI-J-1, NB )
         I  = IA + J - JA
         K  = I - IA + 1
         IV = K - ILO + IOFF + 1
*
*        Compute upper triangular matrix T from TAU.
*
         CALL PDLARFT( 'Forward', 'Columnwise', IHI-K, JB, A, I+1, J,
     $                 DESCA, TAU, WORK( IPT ), WORK( IPW ) )
*
*        Copy Householder vectors into workspace.
*
         CALL PDLACPY( 'All', IHI-K, JB, A, I+1, J, DESCA, WORK( IPV ),
     $                 IV+1, 1, DESCV )
*
*        Zero out the strict lower triangular part of A.
*
         CALL PDLASET( 'Lower', IHI-K-1, JB, ZERO, ZERO, A, I+2, J,
     $                 DESCA )
*
*        Apply block Householder transformation from Left.
*
         CALL PDLARFB( 'Left', 'No transpose', 'Forward', 'Columnwise',
     $                 IHI-K, N-K+1, JB, WORK( IPV ), IV+1, 1, DESCV,
     $                 WORK( IPT ), A, I+1, J, DESCA, WORK( IPW ) )
*
*        Apply block Householder transformation from Right.
*
         CALL PDLARFB( 'Right', 'Transpose', 'Forward', 'Columnwise',
     $                 IHI, IHI-K, JB, WORK( IPV ), IV+1, 1, DESCV,
     $                 WORK( IPT ), A, IA, J+1, DESCA, WORK( IPW ) )
*
         DESCV( CSRC_ ) = MOD( DESCV( CSRC_ ) + NPCOL - 1, NPCOL )
*
   10 CONTINUE
*
*     Handle the first block separately
*
      IV = IOFF + 1
      I = IA + ILO - 1
      J = JA + ILO - 1
      JB = MIN( NB-IOFF, JA+IHI-J-1 )
*
*     Compute upper triangular matrix T from TAU.
*
      CALL PDLARFT( 'Forward', 'Columnwise', IHI-ILO, JB, A, I+1, J,
     $              DESCA, TAU, WORK( IPT ), WORK( IPW ) )
*
*     Copy Householder vectors into workspace.
*
      CALL PDLACPY( 'All', IHI-ILO, JB, A, I+1, J, DESCA, WORK( IPV ),
     $              IV+1, 1, DESCV )
*
*     Zero out the strict lower triangular part of A.
*
      IF( IHI-ILO.GT.0 )
     $   CALL PDLASET( 'Lower', IHI-ILO-1, JB, ZERO, ZERO, A, I+2, J,
     $                 DESCA )
*
*     Apply block Householder transformation from Left.
*
      CALL PDLARFB( 'Left', 'No transpose', 'Forward', 'Columnwise',
     $              IHI-ILO, N-ILO+1, JB, WORK( IPV ), IV+1, 1, DESCV,
     $              WORK( IPT ), A, I+1, J, DESCA, WORK( IPW ) )
*
*     Apply block Householder transformation from Right.
*
      CALL PDLARFB( 'Right', 'Transpose', 'Forward', 'Columnwise', IHI,
     $              IHI-ILO, JB, WORK( IPV ), IV+1, 1, DESCV,
     $              WORK( IPT ), A, IA, J+1, DESCA, WORK( IPW ) )
*
      RETURN
*
*     End of PDGEHDRV
*
      END
