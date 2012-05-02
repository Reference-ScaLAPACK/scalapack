      SUBROUTINE PCGBTRS( TRANS, N, BWL, BWU, NRHS, A, JA, DESCA, IPIV,
     $                    B, IB, DESCB, AF, LAF, WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            BWU, BWL, IB, INFO, JA, LAF, LWORK, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), IPIV(*)
      COMPLEX            A( * ), AF( * ), B( * ), WORK( * )
*     ..
*
*
*  Purpose
*  =======
*
*  PCGBTRS solves a system of linear equations
*
*            A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
*                                    or
*            A(1:N, JA:JA+N-1)' * X = B(IB:IB+N-1, 1:NRHS)
*
*  where A(1:N, JA:JA+N-1) is the matrix used to produce the factors
*  stored in A(1:N,JA:JA+N-1) and AF by PCGBTRF.
*  A(1:N, JA:JA+N-1) is an N-by-N complex
*  banded distributed
*  matrix with bandwidth BWL, BWU.
*
*  Routine PCGBTRF MUST be called first.
*
*  =====================================================================
*
*  Arguments
*  =========
*
*
*  TRANS   (global input) CHARACTER
*          = 'N':  Solve with A(1:N, JA:JA+N-1);
*          = 'C':  Solve with conjugate_transpose( A(1:N, JA:JA+N-1) );
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix A(1:N, JA:JA+N-1). N >= 0.
*
*  BWL     (global input) INTEGER
*          Number of subdiagonals. 0 <= BWL <= N-1
*
*  BWU     (global input) INTEGER
*          Number of superdiagonals. 0 <= BWU <= N-1
*
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the distributed submatrix B(IB:IB+N-1, 1:NRHS).
*          NRHS >= 0.
*
*  A       (local input/local output) COMPLEX pointer into
*          local memory to an array with first dimension
*          LLD_A >=(2*bwl+2*bwu+1) (stored in DESCA).
*          On entry, this array contains the local pieces of the
*          N-by-N unsymmetric banded distributed Cholesky factor L or
*          L^T A(1:N, JA:JA+N-1).
*          This local portion is stored in the packed banded format
*            used in LAPACK. Please see the Notes below and the
*            ScaLAPACK manual for more detail on the format of
*            distributed matrices.
*
*  JA      (global input) INTEGER
*          The index in the global array A that points to the start of
*          the matrix to be operated on (which may be either all of A
*          or a submatrix of A).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN.
*          if 1D type (DTYPE_A=501), DLEN >= 7;
*          if 2D type (DTYPE_A=1), DLEN >= 9 .
*          The array descriptor for the distributed matrix A.
*          Contains information of mapping of A to memory. Please
*          see NOTES below for full description and options.
*
*  IPIV    (local output) INTEGER array, dimension >= DESCA( NB ).
*          Pivot indices for local factorizations.
*          Users *should not* alter the contents between
*          factorization and solve.
*
*  B       (local input/local output) COMPLEX pointer into
*          local memory to an array of local lead dimension lld_b>=NB.
*          On entry, this array contains the
*          the local pieces of the right hand sides
*          B(IB:IB+N-1, 1:NRHS).
*          On exit, this contains the local piece of the solutions
*          distributed matrix X.
*
*  IB      (global input) INTEGER
*          The row index in the global array B that points to the first
*          row of the matrix to be operated on (which may be either
*          all of B or a submatrix of B).
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN.
*          if 1D type (DTYPE_B=502), DLEN >=7;
*          if 2D type (DTYPE_B=1), DLEN >= 9.
*          The array descriptor for the distributed matrix B.
*          Contains information of mapping of B to memory. Please
*          see NOTES below for full description and options.
*
*  AF      (local output) COMPLEX array, dimension LAF.
*          Auxiliary Fillin Space.
*          Fillin is created during the factorization routine
*          PCGBTRF and this is stored in AF. If a linear system
*          is to be solved using PCGBTRS after the factorization
*          routine, AF *must not be altered* after the factorization.
*
*  LAF     (local input) INTEGER
*          Size of user-input Auxiliary Fillin space AF. Must be >=
*          (NB+bwu)*(bwl+bwu)+6*(bwl+bwu)*(bwl+2*bwu)
*          If LAF is not large enough, an error code will be returned
*          and the minimum acceptable size will be returned in AF( 1 )
*
*  WORK    (local workspace/local output)
*          COMPLEX temporary workspace. This space may
*          be overwritten in between calls to routines. WORK must be
*          the size given in LWORK.
*          On exit, WORK( 1 ) contains the minimal LWORK.
*
*  LWORK   (local input or global input) INTEGER
*          Size of user-input workspace WORK.
*          If LWORK is too small, the minimal acceptable size will be
*          returned in WORK(1) and an error code is returned. LWORK>=
*          NRHS*(NB+2*bwl+4*bwu)
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  =====================================================================
*
*
*  Restrictions
*  ============
*
*  The following are restrictions on the input parameters. Some of these
*    are temporary and will be removed in future releases, while others
*    may reflect fundamental technical limitations.
*
*    Non-cyclic restriction: VERY IMPORTANT!
*      P*NB>= mod(JA-1,NB)+N.
*      The mapping for matrices must be blocked, reflecting the nature
*      of the divide and conquer algorithm as a task-parallel algorithm.
*      This formula in words is: no processor may have more than one
*      chunk of the matrix.
*
*    Blocksize cannot be too small:
*      If the matrix spans more than one processor, the following
*      restriction on NB, the size of each block on each processor,
*      must hold:
*      NB >= (BWL+BWU)+1
*      The bulk of parallel computation is done on the matrix of size
*      O(NB) on each processor. If this is too small, divide and conquer
*      is a poor choice of algorithm.
*
*    Submatrix reference:
*      JA = IB
*      Alignment restriction that prevents unnecessary communication.
*
*
*  =====================================================================
*
*
*  Notes
*  =====
*
*  If the factorization routine and the solve routine are to be called
*    separately (to solve various sets of righthand sides using the same
*    coefficient matrix), the auxiliary space AF *must not be altered*
*    between calls to the factorization routine and the solve routine.
*
*  The best algorithm for solving banded and tridiagonal linear systems
*    depends on a variety of parameters, especially the bandwidth.
*    Currently, only algorithms designed for the case N/P >> bw are
*    implemented. These go by many names, including Divide and Conquer,
*    Partitioning, domain decomposition-type, etc.
*
*  Algorithm description: Divide and Conquer
*
*    The Divide and Conqer algorithm assumes the matrix is narrowly
*      banded compared with the number of equations. In this situation,
*      it is best to distribute the input matrix A one-dimensionally,
*      with columns atomic and rows divided amongst the processes.
*      The basic algorithm divides the banded matrix up into
*      P pieces with one stored on each processor,
*      and then proceeds in 2 phases for the factorization or 3 for the
*      solution of a linear system.
*      1) Local Phase:
*         The individual pieces are factored independently and in
*         parallel. These factors are applied to the matrix creating
*         fillin, which is stored in a non-inspectable way in auxiliary
*         space AF. Mathematically, this is equivalent to reordering
*         the matrix A as P A P^T and then factoring the principal
*         leading submatrix of size equal to the sum of the sizes of
*         the matrices factored on each processor. The factors of
*         these submatrices overwrite the corresponding parts of A
*         in memory.
*      2) Reduced System Phase:
*         A small (max(bwl,bwu)* (P-1)) system is formed representing
*         interaction of the larger blocks, and is stored (as are its
*         factors) in the space AF. A parallel Block Cyclic Reduction
*         algorithm is used. For a linear system, a parallel front solve
*         followed by an analagous backsolve, both using the structure
*         of the factored matrix, are performed.
*      3) Backsubsitution Phase:
*         For a linear system, a local backsubstitution is performed on
*         each processor in parallel.
*
*
*  Descriptors
*  ===========
*
*  Descriptors now have *types* and differ from ScaLAPACK 1.0.
*
*  Note: banded codes can use either the old two dimensional
*    or new one-dimensional descriptors, though the processor grid in
*    both cases *must be one-dimensional*. We describe both types below.
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
*
*  One-dimensional descriptors:
*
*  One-dimensional descriptors are a new addition to ScaLAPACK since
*    version 1.0. They simplify and shorten the descriptor for 1D
*    arrays.
*
*  Since ScaLAPACK supports two-dimensional arrays as the fundamental
*    object, we allow 1D arrays to be distributed either over the
*    first dimension of the array (as if the grid were P-by-1) or the
*    2nd dimension (as if the grid were 1-by-P). This choice is
*    indicated by the descriptor type (501 or 502)
*    as described below.
*
*    IMPORTANT NOTE: the actual BLACS grid represented by the
*    CTXT entry in the descriptor may be *either*  P-by-1 or 1-by-P
*    irrespective of which one-dimensional descriptor type
*    (501 or 502) is input.
*    This routine will interpret the grid properly either way.
*    ScaLAPACK routines *do not support intercontext operations* so that
*    the grid passed to a single ScaLAPACK routine *must be the same*
*    for all array descriptors passed to that routine.
*
*    NOTE: In all cases where 1D descriptors are used, 2D descriptors
*    may also be used, since a one-dimensional array is a special case
*    of a two-dimensional array with one dimension of size unity.
*    The two-dimensional array used in this case *must* be of the
*    proper orientation:
*      If the appropriate one-dimensional descriptor is DTYPEA=501
*      (1 by P type), then the two dimensional descriptor must
*      have a CTXT value that refers to a 1 by P BLACS grid;
*      If the appropriate one-dimensional descriptor is DTYPEA=502
*      (P by 1 type), then the two dimensional descriptor must
*      have a CTXT value that refers to a P by 1 BLACS grid.
*
*
*  Summary of allowed descriptors, types, and BLACS grids:
*  DTYPE           501         502         1         1
*  BLACS grid      1xP or Px1  1xP or Px1  1xP       Px1
*  -----------------------------------------------------
*  A               OK          NO          OK        NO
*  B               NO          OK          NO        OK
*
*  Note that a consequence of this chart is that it is not possible
*    for *both* DTYPE_A and DTYPE_B to be 2D_type(1), as these lead
*    to opposite requirements for the orientation of the BLACS grid,
*    and as noted before, the *same* BLACS context must be used in
*    all descriptors in a single ScaLAPACK subroutine call.
*
*  Let A be a generic term for any 1D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN  EXPLANATION
*  --------------- ---------- ------------------------------------------
*  DTYPE_A(global) DESCA( 1 ) The descriptor type. For 1D grids,
*                                TYPE_A = 501: 1-by-P grid.
*                                TYPE_A = 502: P-by-1 grid.
*  CTXT_A (global) DESCA( 2 ) The BLACS context handle, indicating
*                                the BLACS process grid A is distribu-
*                                ted over. The context itself is glo-
*                                bal, but the handle (the integer
*                                value) may vary.
*  N_A    (global) DESCA( 3 ) The size of the array dimension being
*                                distributed.
*  NB_A   (global) DESCA( 4 ) The blocking factor used to distribute
*                                the distributed dimension of the array.
*  SRC_A  (global) DESCA( 5 ) The process row or column over which the
*                                first row or column of the array
*                                is distributed.
*  LLD_A  (local)  DESCA( 6 ) The leading dimension of the local array
*                                storing the local blocks of the distri-
*                                buted array A. Minimum value of LLD_A
*                                depends on TYPE_A.
*                                TYPE_A = 501: LLD_A >=
*                                   size of undistributed dimension, 1.
*                                TYPE_A = 502: LLD_A >=NB_A, 1.
*  Reserved        DESCA( 7 ) Reserved for future use.
*
*
*
*  =====================================================================
*
*  Implemented for ScaLAPACK by:
*     Andrew J. Cleary, Livermore National Lab and University of Tenn.,
*     and Marbwus Hegland, Australian Natonal University. Feb., 1997.
*  Based on code written by    : Peter Arbenz, ETH Zurich, 1996.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0 )
      PARAMETER          ( ZERO = 0.0E+0 )
      COMPLEX            CONE, CZERO
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
      INTEGER            INT_ONE
      PARAMETER          ( INT_ONE = 1 )
      INTEGER            DESCMULT, BIGNUM
      PARAMETER          ( DESCMULT = 100, BIGNUM = DESCMULT*DESCMULT )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            APTR, BBPTR, BM, BMN, BN, BNN, BW, CSRC,
     $                   FIRST_PROC, ICTXT, ICTXT_NEW, ICTXT_SAVE,
     $                   IDUM2, IDUM3, J, JA_NEW, L, LBWL, LBWU, LDBB,
     $                   LDW, LLDA, LLDB, LM, LMJ, LN, LPTR, MYCOL,
     $                   MYROW, NB, NEICOL, NP, NPACT, NPCOL, NPROW,
     $                   NPSTR, NP_SAVE, ODD_SIZE, PART_OFFSET,
     $                   RECOVERY_VAL, RETURN_CODE, STORE_M_B,
     $                   STORE_N_A, WORK_SIZE_MIN, WPTR
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA_1XP( 7 ), DESCB_PX1( 7 ),
     $                   PARAM_CHECK( 17, 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESC_CONVERT, GLOBCHK, PXERBLA,
     $                   RESHAPE
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      EXTERNAL           LSAME
      EXTERNAL           NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, MOD
*     ..
*     .. Executable Statements ..
*
*
*     Test the input parameters
*
      INFO = 0
*
*     Convert descriptor into standard form for easy access to
*        parameters, check that grid is of right shape.
*
      DESCA_1XP( 1 ) = 501
      DESCB_PX1( 1 ) = 502
*
      CALL DESC_CONVERT( DESCA, DESCA_1XP, RETURN_CODE )
*
      IF( RETURN_CODE .NE. 0) THEN
         INFO = -( 8*100 + 2 )
      ENDIF
*
      CALL DESC_CONVERT( DESCB, DESCB_PX1, RETURN_CODE )
*
      IF( RETURN_CODE .NE. 0) THEN
         INFO = -( 11*100 + 2 )
      ENDIF
*
*     Consistency checks for DESCA and DESCB.
*
*     Context must be the same
      IF( DESCA_1XP( 2 ) .NE. DESCB_PX1( 2 ) ) THEN
         INFO = -( 11*100 + 2 )
      ENDIF
*
*        These are alignment restrictions that may or may not be removed
*        in future releases. -Andy Cleary, April 14, 1996.
*
*     Block sizes must be the same
      IF( DESCA_1XP( 4 ) .NE. DESCB_PX1( 4 ) ) THEN
         INFO = -( 11*100 + 4 )
      ENDIF
*
*     Source processor must be the same
*
      IF( DESCA_1XP( 5 ) .NE. DESCB_PX1( 5 ) ) THEN
         INFO = -( 11*100 + 5 )
      ENDIF
*
*     Get values out of descriptor for use in code.
*
      ICTXT = DESCA_1XP( 2 )
      CSRC = DESCA_1XP( 5 )
      NB = DESCA_1XP( 4 )
      LLDA = DESCA_1XP( 6 )
      STORE_N_A = DESCA_1XP( 3 )
      LLDB = DESCB_PX1( 6 )
      STORE_M_B = DESCB_PX1( 3 )
*
*     Get grid parameters
*
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NP = NPROW * NPCOL
*
*
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         IDUM2 = ICHAR( 'N' )
      ELSE IF ( LSAME( TRANS, 'C' ) ) THEN
         IDUM2 = ICHAR( 'C' )
      ELSE
         INFO = -1
      END IF
*
      IF( LWORK .LT. -1) THEN
         INFO = -16
      ELSE IF ( LWORK .EQ. -1 ) THEN
         IDUM3 = -1
      ELSE
         IDUM3 = 1
      ENDIF
*
      IF( N .LT. 0 ) THEN
         INFO = -2
      ENDIF
*
      IF( N+JA-1 .GT. STORE_N_A ) THEN
         INFO = -( 8*100 + 6 )
      ENDIF
*
      IF(( BWL .GT. N-1 ) .OR.
     $   ( BWL .LT. 0 ) ) THEN
         INFO = -3
      ENDIF
*
      IF(( BWU .GT. N-1 ) .OR.
     $   ( BWU .LT. 0 ) ) THEN
         INFO = -4
      ENDIF
*
      IF( LLDA .LT. (2*BWL+2*BWU+1) ) THEN
         INFO = -( 8*100 + 6 )
      ENDIF
*
      IF( NB .LE. 0 ) THEN
         INFO = -( 8*100 + 4 )
      ENDIF
*
      BW = BWU+BWL
*
      IF( N+IB-1 .GT. STORE_M_B ) THEN
         INFO = -( 11*100 + 3 )
      ENDIF
*
      IF( LLDB .LT. NB ) THEN
         INFO = -( 11*100 + 6 )
      ENDIF
*
      IF( NRHS .LT. 0 ) THEN
         INFO = -5
      ENDIF
*
*     Current alignment restriction
*
      IF( JA .NE. IB) THEN
         INFO = -7
      ENDIF
*
*     Argument checking that is specific to Divide & Conquer routine
*
      IF( NPROW .NE. 1 ) THEN
         INFO = -( 8*100+2 )
      ENDIF
*
      IF( N .GT. NP*NB-MOD( JA-1, NB )) THEN
         INFO = -( 2 )
         CALL PXERBLA( ICTXT,
     $      'PCGBTRS, D&C alg.: only 1 block per proc',
     $      -INFO )
         RETURN
      ENDIF
*
      IF((JA+N-1.GT.NB) .AND. ( NB.LT.(BWL+BWU+1) )) THEN
         INFO = -( 8*100+4 )
         CALL PXERBLA( ICTXT,
     $      'PCGBTRS, D&C alg.: NB too small',
     $      -INFO )
         RETURN
      ENDIF
*
*
*     Check worksize
*
      WORK_SIZE_MIN = NRHS*(NB+2*BWL+4*BWU)
*
      WORK( 1 ) = WORK_SIZE_MIN
*
      IF( LWORK .LT. WORK_SIZE_MIN ) THEN
         IF( LWORK .NE. -1 ) THEN
         INFO = -16
         CALL PXERBLA( ICTXT,
     $      'PCGBTRS: worksize error ',
     $      -INFO )
         ENDIF
         RETURN
      ENDIF
*
*     Pack params and positions into arrays for global consistency check
*
      PARAM_CHECK( 17, 1 ) = DESCB(5)
      PARAM_CHECK( 16, 1 ) = DESCB(4)
      PARAM_CHECK( 15, 1 ) = DESCB(3)
      PARAM_CHECK( 14, 1 ) = DESCB(2)
      PARAM_CHECK( 13, 1 ) = DESCB(1)
      PARAM_CHECK( 12, 1 ) = IB
      PARAM_CHECK( 11, 1 ) = DESCA(5)
      PARAM_CHECK( 10, 1 ) = DESCA(4)
      PARAM_CHECK(  9, 1 ) = DESCA(3)
      PARAM_CHECK(  8, 1 ) = DESCA(1)
      PARAM_CHECK(  7, 1 ) = JA
      PARAM_CHECK(  6, 1 ) = NRHS
      PARAM_CHECK(  5, 1 ) = BWU
      PARAM_CHECK(  4, 1 ) = BWL
      PARAM_CHECK(  3, 1 ) = N
      PARAM_CHECK(  2, 1 ) = IDUM3
      PARAM_CHECK(  1, 1 ) = IDUM2
*
      PARAM_CHECK( 17, 2 ) = 1105
      PARAM_CHECK( 16, 2 ) = 1104
      PARAM_CHECK( 15, 2 ) = 1103
      PARAM_CHECK( 14, 2 ) = 1102
      PARAM_CHECK( 13, 2 ) = 1101
      PARAM_CHECK( 12, 2 ) = 10
      PARAM_CHECK( 11, 2 ) = 805
      PARAM_CHECK( 10, 2 ) = 804
      PARAM_CHECK(  9, 2 ) = 803
      PARAM_CHECK(  8, 2 ) = 801
      PARAM_CHECK(  7, 2 ) = 7
      PARAM_CHECK(  6, 2 ) = 5
      PARAM_CHECK(  5, 2 ) = 4
      PARAM_CHECK(  4, 2 ) = 3
      PARAM_CHECK(  3, 2 ) = 2
      PARAM_CHECK(  2, 2 ) = 16
      PARAM_CHECK(  1, 2 ) = 1
*
*     Want to find errors with MIN( ), so if no error, set it to a big
*     number. If there already is an error, multiply by the the
*     descriptor multiplier.
*
      IF( INFO.GE.0 ) THEN
         INFO = BIGNUM
      ELSE IF( INFO.LT.-DESCMULT ) THEN
         INFO = -INFO
      ELSE
         INFO = -INFO * DESCMULT
      END IF
*
*     Check consistency across processors
*
      CALL GLOBCHK( ICTXT, 17, PARAM_CHECK, 17,
     $              PARAM_CHECK( 1, 3 ), INFO )
*
*     Prepare output: set info = 0 if no error, and divide by DESCMULT
*     if error is not in a descriptor entry.
*
      IF( INFO.EQ.BIGNUM ) THEN
         INFO = 0
      ELSE IF( MOD( INFO, DESCMULT ) .EQ. 0 ) THEN
         INFO = -INFO / DESCMULT
      ELSE
         INFO = -INFO
      END IF
*
      IF( INFO.LT.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCGBTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( NRHS.EQ.0 )
     $   RETURN
*
*
*     Adjust addressing into matrix space to properly get into
*        the beginning part of the relevant data
*
      PART_OFFSET = NB*( (JA-1)/(NPCOL*NB) )
*
      IF ( (MYCOL-CSRC) .LT. (JA-PART_OFFSET-1)/NB ) THEN
         PART_OFFSET = PART_OFFSET + NB
      ENDIF
*
      IF ( MYCOL .LT. CSRC ) THEN
         PART_OFFSET = PART_OFFSET - NB
      ENDIF
*
*     Form a new BLACS grid (the "standard form" grid) with only procs
*        holding part of the matrix, of size 1xNP where NP is adjusted,
*        starting at csrc=0, with JA modified to reflect dropped procs.
*
*     First processor to hold part of the matrix:
*
      FIRST_PROC = MOD( ( JA-1 )/NB+CSRC, NPCOL )
*
*     Calculate new JA one while dropping off unused processors.
*
      JA_NEW = MOD( JA-1, NB ) + 1
*
*     Save and compute new value of NP
*
      NP_SAVE = NP
      NP = ( JA_NEW+N-2 )/NB + 1
*
*     Call utility routine that forms "standard-form" grid
*
      CALL RESHAPE( ICTXT, INT_ONE, ICTXT_NEW, INT_ONE,
     $              FIRST_PROC, INT_ONE, NP )
*
*     Use new context from standard grid as context.
*
      ICTXT_SAVE = ICTXT
      ICTXT = ICTXT_NEW
      DESCA_1XP( 2 ) = ICTXT_NEW
      DESCB_PX1( 2 ) = ICTXT_NEW
*
*     Get information about new grid.
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Drop out processors that do not have part of the matrix.
*
      IF( MYROW .LT. 0 ) THEN
         GOTO 1234
      ENDIF
*
*
*
*     Begin main code
*
*     Move data into workspace - communicate/copy (overlap)
*
      IF (MYCOL .LT. NPCOL-1) THEN
         CALL CGESD2D( ICTXT, BWU, NRHS, B(NB-BWU+1), LLDB,
     $        0, MYCOL + 1)
      ENDIF
*
      IF (MYCOL .LT. NPCOL-1) THEN
         LM = NB-BWU
      ELSE
         LM = NB
      ENDIF
*
      IF (MYCOL .GT. 0) THEN
         WPTR = BWU+1
      ELSE
         WPTR = 1
      ENDIF
*
      LDW = NB+BWU + 2*BW+BWU
*
      CALL CLAMOV( 'G', LM, NRHS, B(1), LLDB, WORK( WPTR ), LDW )
*
*     Zero out rest of work
*
      DO 1501 J=1, NRHS
        DO 1502 L=WPTR+LM, LDW
          WORK( (J-1)*LDW+L ) = CZERO
 1502   CONTINUE
 1501 CONTINUE
*
      IF (MYCOL .GT. 0) THEN
         CALL CGERV2D( ICTXT, BWU, NRHS, WORK(1), LDW,
     $        0, MYCOL-1)
      ENDIF
*
********************************************************************
*       PHASE 1: Local computation phase -- Solve L*X = B
********************************************************************
*
*     Size of main (or odd) partition in each processor
*
      ODD_SIZE = NUMROC( N, NB, MYCOL, 0, NPCOL )
*
      IF (MYCOL .NE. 0) THEN
         LBWL = BW
         LBWU = 0
         APTR = 1
      ELSE
         LBWL = BWL
         LBWU = BWU
         APTR = 1+BWU
      ENDIF
*
      IF (MYCOL .NE. NPCOL-1) THEN
         LM = NB - LBWU
         LN = NB - BW
      ELSE IF (MYCOL .NE. 0) THEN
         LM = ODD_SIZE + BWU
         LN = MAX(ODD_SIZE-BW,0)
      ELSE
         LM = N
         LN = MAX( N-BW, 0 )
      ENDIF
*
      DO 21 J = 1, LN
*
         LMJ = MIN(LBWL,LM-J)
         L = IPIV( J )
*
         IF( L.NE.J ) THEN
            CALL CSWAP(NRHS, WORK(L), LDW, WORK(J), LDW)
         ENDIF
*
         LPTR = BW+1 + (J-1)*LLDA + APTR
*
         CALL CGERU(LMJ,NRHS,-CONE, A(LPTR),1, WORK(J),LDW,
     $           WORK(J+1),LDW)
*
   21 CONTINUE
*
********************************************************************
*       PHASE 2: Global computation phase -- Solve L*X = B
********************************************************************
*
*     Define the initial dimensions of the diagonal blocks
*     The offdiagonal blocks (for MYCOL > 0) are of size BM by BW
*
      IF (MYCOL .NE. NPCOL-1) THEN
         BM = BW - LBWU
         BN = BW
      ELSE
         BM = MIN(BW,ODD_SIZE) + BWU
         BN = MIN(BW,ODD_SIZE)
      ENDIF
*
*     Pointer to first element of block bidiagonal matrix in AF
*     Leading dimension of block bidiagonal system
*
      BBPTR = (NB+BWU)*BW + 1
      LDBB   = 2*BW + BWU
*
      IF (NPCOL.EQ.1) THEN
*
*        In this case the loop over the levels will not be
*        performed.
         CALL CGETRS( 'N', N-LN, NRHS, AF(BBPTR+BW*LDBB), LDBB,
     $        IPIV(LN+1), WORK( LN+1 ), LDW, INFO)
*
      ENDIF
*
* Loop over levels ...
*
*     The two integers NPACT (nu. of active processors) and NPSTR
*     (stride between active processors) is used to control the
*     loop.
*
      NPACT = NPCOL
      NPSTR = 1
*
*     Begin loop over levels
  200 IF (NPACT .LE. 1) GOTO 300
*
*     Test if processor is active
          IF (MOD(MYCOL,NPSTR) .EQ. 0) THEN
*
*   Send/Receive blocks
*
             IF (MOD(MYCOL,2*NPSTR) .EQ. 0) THEN
*
                NEICOL = MYCOL + NPSTR
*
                IF (NEICOL/NPSTR .LE. NPACT-1) THEN
*
                   IF (NEICOL/NPSTR .LT. NPACT-1) THEN
                      BMN = BW
                   ELSE
                      BMN = MIN(BW,NUMROC(N, NB, NEICOL, 0, NPCOL))+BWU
                   ENDIF
*
                   CALL CGESD2D( ICTXT, BM, NRHS,
     $                  WORK(LN+1), LDW, 0, NEICOL )
*
                   IF( NPACT .NE. 2 )THEN
*
*                     Receive answers back from partner processor
*
                      CALL CGERV2D(ICTXT, BM+BMN-BW, NRHS,
     $                   WORK( LN+1 ), LDW, 0, NEICOL )
*
                      BM = BM+BMN-BW
*
                   ENDIF
*
                ENDIF
*
             ELSE
*
                NEICOL = MYCOL - NPSTR
*
                IF (NEICOL .EQ. 0) THEN
                   BMN = BW - BWU
                ELSE
                   BMN = BW
                ENDIF
*
                CALL CLAMOV( 'G', BM, NRHS, WORK(LN+1), LDW,
     $               WORK(NB+BWU+BMN+1), LDW )
*
                CALL CGERV2D( ICTXT, BMN, NRHS, WORK( NB+BWU+1 ),
     $                  LDW, 0, NEICOL )
*
*               and do the permutations and eliminations
*
                IF (NPACT .NE. 2) THEN
*
*                  Solve locally for BW variables
*
                   CALL CLASWP( NRHS, WORK(NB+BWU+1), LDW, 1, BW,
     $                  IPIV(LN+1), 1)
*
                   CALL CTRSM('L','L','N','U', BW, NRHS, CONE,
     $                  AF(BBPTR+BW*LDBB), LDBB, WORK(NB+BWU+1), LDW)
*
*                  Use soln just calculated to update RHS
*
                   CALL CGEMM( 'N', 'N', BM+BMN-BW, NRHS, BW,
     $                -CONE, AF(BBPTR+BW*LDBB+BW), LDBB,
     $                WORK(NB+BWU+1), LDW,
     $                CONE, WORK(NB+BWU+1+BW), LDW )
*
*                  Give answers back to partner processor
*
                   CALL CGESD2D( ICTXT, BM+BMN-BW, NRHS,
     $                WORK(NB+BWU+1+BW), LDW, 0, NEICOL )
*
                ELSE
*
*                  Finish up calculations for final level
*
                   CALL CLASWP( NRHS, WORK(NB+BWU+1), LDW, 1, BM+BMN,
     $                  IPIV(LN+1), 1)
*
                   CALL CTRSM('L','L','N','U', BM+BMN, NRHS, CONE,
     $                  AF(BBPTR+BW*LDBB), LDBB, WORK(NB+BWU+1), LDW)
                ENDIF
*
             ENDIF
*
             NPACT = (NPACT + 1)/2
             NPSTR = NPSTR * 2
             GOTO 200
*
         ENDIF
*
  300 CONTINUE
*
*
**************************************
*     BACKSOLVE
********************************************************************
*       PHASE 2: Global computation phase -- Solve U*Y = X
********************************************************************
*
      IF (NPCOL.EQ.1) THEN
*
*        In this case the loop over the levels will not be
*        performed.
*        In fact, the backsolve portion was done in the call to
*          CGETRS in the frontsolve.
*
      ENDIF
*
*     Compute variable needed to reverse loop structure in
*        reduced system.
*
      RECOVERY_VAL = NPACT*NPSTR - NPCOL
*
*     Loop over levels
*      Terminal values of NPACT and NPSTR from frontsolve are used
*
 2200 IF( NPACT .GE. NPCOL ) GOTO 2300
*
         NPSTR = NPSTR/2
*
         NPACT = NPACT*2
*
*        Have to adjust npact for non-power-of-2
*
         NPACT = NPACT-MOD( (RECOVERY_VAL/NPSTR), 2 )
*
*        Find size of submatrix in this proc at this level
*
         IF( MYCOL/NPSTR .LT. NPACT-1 ) THEN
            BN = BW
         ELSE
            BN = MIN(BW, NUMROC(N, NB, NPCOL-1, 0, NPCOL) )
         ENDIF
*
*        If this processor is even in this level...
*
         IF( MOD( MYCOL, 2*NPSTR ) .EQ. 0 ) THEN
*
            NEICOL = MYCOL+NPSTR
*
            IF( NEICOL/NPSTR .LE. NPACT-1 ) THEN
*
               IF( NEICOL/NPSTR .LT. NPACT-1 ) THEN
                  BMN = BW
                  BNN = BW
               ELSE
                  BMN = MIN(BW,NUMROC(N, NB, NEICOL, 0, NPCOL))+BWU
                  BNN = MIN(BW, NUMROC(N, NB, NEICOL, 0, NPCOL) )
               ENDIF
*
               IF( NPACT .GT. 2 ) THEN
*
                  CALL CGESD2D( ICTXT, 2*BW, NRHS, WORK( LN+1 ),
     $                  LDW, 0, NEICOL )
*
                  CALL CGERV2D( ICTXT, BW, NRHS, WORK( LN+1 ),
     $                  LDW, 0, NEICOL )
*
               ELSE
*
                  CALL CGERV2D( ICTXT, BW, NRHS, WORK( LN+1 ),
     $                  LDW, 0, NEICOL )
*
               ENDIF
*
            ENDIF
*
         ELSE
*           This processor is odd on this level
*
            NEICOL = MYCOL - NPSTR
*
            IF (NEICOL .EQ. 0) THEN
               BMN = BW - BWU
            ELSE
               BMN = BW
            ENDIF
*
            IF( NEICOL .LT. NPCOL-1 ) THEN
               BNN = BW
            ELSE
               BNN = MIN(BW, NUMROC(N, NB, NEICOL, 0, NPCOL) )
            ENDIF
*
            IF( NPACT .GT. 2 ) THEN
*
*              Move RHS to make room for received solutions
*
               CALL CLAMOV( 'G', BW, NRHS, WORK(NB+BWU+1),
     $               LDW, WORK(NB+BWU+BW+1), LDW )
*
               CALL CGERV2D( ICTXT, 2*BW, NRHS, WORK( LN+1 ),
     $                  LDW, 0, NEICOL )
*
               CALL CGEMM( 'N', 'N', BW, NRHS, BN,
     $                -CONE, AF(BBPTR), LDBB,
     $                WORK(LN+1), LDW,
     $                CONE, WORK(NB+BWU+BW+1), LDW )
*
*
               IF( MYCOL .GT. NPSTR ) THEN
*
                  CALL CGEMM( 'N', 'N', BW, NRHS, BW,
     $                -CONE, AF(BBPTR+2*BW*LDBB), LDBB,
     $                WORK(LN+BW+1), LDW,
     $                CONE, WORK(NB+BWU+BW+1), LDW )
*
               ENDIF
*
               CALL CTRSM('L','U','N','N', BW, NRHS, CONE,
     $                AF(BBPTR+BW*LDBB), LDBB, WORK(NB+BWU+BW+1), LDW)
*
*              Send new solution to neighbor
*
               CALL CGESD2D( ICTXT, BW, NRHS,
     $                WORK( NB+BWU+BW+1 ), LDW, 0, NEICOL )
*
*              Copy new solution into expected place
*
               CALL CLAMOV( 'G', BW, NRHS, WORK(NB+BWU+1+BW),
     $               LDW, WORK(LN+BW+1), LDW )
*
            ELSE
*
*              Solve with local diagonal block
*
               CALL CTRSM( 'L','U','N','N', BN+BNN, NRHS, CONE,
     $                  AF(BBPTR+BW*LDBB), LDBB, WORK(NB+BWU+1), LDW)
*
*              Send new solution to neighbor
*
               CALL CGESD2D( ICTXT, BW, NRHS,
     $                WORK(NB+BWU+1), LDW, 0, NEICOL )
*
*              Shift solutions into expected positions
*
               CALL CLAMOV( 'G', BNN+BN-BW, NRHS, WORK(NB+BWU+1+BW),
     $               LDW, WORK(LN+1), LDW )
*
*
               IF( (NB+BWU+1) .NE. (LN+1+BW) ) THEN
*
*                 Copy one row at a time since spaces may overlap
*
                  DO 1064 J=1, BW
                     CALL CCOPY( NRHS, WORK(NB+BWU+J), LDW,
     $                                      WORK(LN+BW+J), LDW )
 1064             CONTINUE
*
               ENDIF
*
            ENDIF
*
         ENDIF
*
      GOTO 2200
*
 2300 CONTINUE
*     End of loop over levels
*
********************************************************************
*       PHASE 1: (Almost) Local computation phase -- Solve U*Y = X
********************************************************************
*
*     Reset BM to value it had before reduced system frontsolve...
*
      IF (MYCOL .NE. NPCOL-1) THEN
         BM = BW - LBWU
      ELSE
         BM = MIN(BW,ODD_SIZE) + BWU
      ENDIF
*
*     First metastep is to account for the fillin blocks AF
*
      IF( MYCOL .LT. NPCOL-1 ) THEN
*
         CALL CGESD2D( ICTXT, BW, NRHS, WORK( NB-BW+1 ),
     $                  LDW, 0, MYCOL+1 )
*
      ENDIF
*
      IF( MYCOL .GT. 0 ) THEN
*
         CALL CGERV2D( ICTXT, BW, NRHS, WORK( NB+BWU+1 ),
     $                  LDW, 0, MYCOL-1 )
*
*        Modify local right hand sides with received rhs's
*
         CALL CGEMM( 'N', 'N', LM-BM, NRHS, BW, -CONE,
     $           AF( 1 ), LM, WORK( NB+BWU+1 ), LDW, CONE,
     $           WORK( 1 ), LDW )
*
      ENDIF
*
      DO 2021 J = LN, 1, -1
*
         LMJ = MIN( BW, ODD_SIZE-1 )
*
         LPTR = BW-1+J*LLDA+APTR
*
*        In the following, the TRANS=T option is used to reverse
*           the order of multiplication, not as a true transpose
*
         CALL CGEMV( 'T', LMJ, NRHS, -CONE, WORK( J+1), LDW,
     $           A( LPTR ), LLDA-1, CONE, WORK( J ), LDW )
*
*        Divide by diagonal element
*
         CALL CSCAL( NRHS, CONE/A( LPTR-LLDA+1 ),
     $               WORK( J ), LDW )
 2021 CONTINUE
*
*
*
      CALL CLAMOV( 'G', ODD_SIZE, NRHS, WORK( 1 ), LDW,
     $             B( 1 ), LLDB )
*
*     Free BLACS space used to hold standard-form grid.
*
      ICTXT = ICTXT_SAVE
      IF( ICTXT .NE. ICTXT_NEW ) THEN
         CALL BLACS_GRIDEXIT( ICTXT_NEW )
      ENDIF
*
 1234 CONTINUE
*
*     Restore saved input parameters
*
      NP = NP_SAVE
*
*     Output worksize
*
      WORK( 1 ) = WORK_SIZE_MIN
*
      RETURN
*
*     End of PCGBTRS
*
      END
