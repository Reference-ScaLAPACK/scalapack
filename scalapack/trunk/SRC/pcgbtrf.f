      SUBROUTINE PCGBTRF( N, BWL, BWU, A, JA, DESCA, IPIV, AF, LAF,
     $                    WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      INTEGER            BWL, BWU, INFO, JA, LAF, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), IPIV( * )
      COMPLEX            A( * ), AF( * ), WORK( * )
*     ..
*
*
*  Purpose
*  =======
*
*  PCGBTRF computes a LU factorization
*  of an N-by-N complex banded
*  distributed matrix
*  with bandwidth BWL, BWU: A(1:N, JA:JA+N-1).
*  Reordering is used to increase parallelism in the factorization.
*  This reordering results in factors that are DIFFERENT from those
*  produced by equivalent sequential codes. These factors cannot
*  be used directly by users; however, they can be used in
*  subsequent calls to PCGBTRS to solve linear systems.
*
*  The factorization has the form
*
*          P A(1:N, JA:JA+N-1) Q = L U
*
*  where U is a banded upper triangular matrix and L is banded
*  lower triangular, and P and Q are permutation matrices.
*  The matrix Q represents reordering of columns
*  for parallelism's sake, while P represents
*  reordering of rows for numerical stability using
*  classic partial pivoting.
*
*  =====================================================================
*
*  Arguments
*  =========
*
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
*  A       (local input/local output) COMPLEX pointer into
*          local memory to an array with first dimension
*          LLD_A >=(2*bwl+2*bwu+1) (stored in DESCA).
*          On entry, this array contains the local pieces of the
*          N-by-N unsymmetric banded distributed matrix
*          A(1:N, JA:JA+N-1) to be factored.
*          This local portion is stored in the packed banded format
*            used in LAPACK. Please see the Notes below and the
*            ScaLAPACK manual for more detail on the format of
*            distributed matrices.
*          On exit, this array contains information containing details
*            of the factorization.
*          Note that permutations are performed on the matrix, so that
*            the factors returned are different from those returned
*            by LAPACK.
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
*          1
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = K<=NPROCS, the submatrix stored on processor
*                INFO and factored locally was not
*                nonsingular,  and
*                the factorization was not completed.
*                If INFO = K>NPROCS, the submatrix stored on processor
*                INFO-NPROCS representing interactions with other
*                processors was not
*                nonsingular,
*                and the factorization was not completed.
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
*     and Markus Hegland, Australian Natonal University. Feb., 1997.
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
      INTEGER            APTR, BBPTR, BIPTR, BM, BM1, BM2, BMN, BN, BW,
     $                   CSRC, DBPTR, FIRST_PROC, I, ICTXT, ICTXT_NEW,
     $                   ICTXT_SAVE, IDUM3, J, JA_NEW, JPTR, L, LAF_MIN,
     $                   LBWL, LBWU, LDB, LDBB, LLDA, LM, LMJ, LN, LNJ,
     $                   LPTR, MYCOL, MYROW, MY_NUM_COLS, NB, NEICOL,
     $                   NP, NPACT, NPCOL, NPROW, NPSTR, NP_SAVE, NRHS,
     $                   ODD_N, ODD_SIZE, ODPTR, OFST, PART_OFFSET,
     $                   PART_SIZE, RETURN_CODE, STORE_N_A,
     $                   WORK_SIZE_MIN
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA_1XP( 7 ), PARAM_CHECK( 9, 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDEXIT, BLACS_GRIDINFO, CAXPY, CGEMM,
     $                   CGERV2D, CGESD2D, CLAMOV, CLATCPY, CPBTRF,
     $                   CPOTRF, CSYRK, CTBTRS, CTRMM, CTRRV2D, CTRSD2D,
     $                   CTRSM, CTRTRS, DESC_CONVERT, GLOBCHK, PXERBLA,
     $                   RESHAPE
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      EXTERNAL           LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, MIN, MOD
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
*
      CALL DESC_CONVERT( DESCA, DESCA_1XP, RETURN_CODE )
*
      IF( RETURN_CODE .NE. 0) THEN
         INFO = -( 6*100 + 2 )
      ENDIF
*
*     Get values out of descriptor for use in code.
*
      ICTXT = DESCA_1XP( 2 )
      CSRC = DESCA_1XP( 5 )
      NB = DESCA_1XP( 4 )
      LLDA = DESCA_1XP( 6 )
      STORE_N_A = DESCA_1XP( 3 )
*
*     Get grid parameters
*
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NP = NPROW * NPCOL
*
*
*
      IF( LWORK .LT. -1) THEN
         INFO = -11
      ELSE IF ( LWORK .EQ. -1 ) THEN
         IDUM3 = -1
      ELSE
         IDUM3 = 1
      ENDIF
*
      IF( N .LT. 0 ) THEN
         INFO = -1
      ENDIF
*
      IF( N+JA-1 .GT. STORE_N_A ) THEN
         INFO = -( 6*100 + 6 )
      ENDIF
*
      IF(( BWL .GT. N-1 ) .OR.
     $   ( BWL .LT. 0 ) ) THEN
         INFO = -2
      ENDIF
*
      IF(( BWU .GT. N-1 ) .OR.
     $   ( BWU .LT. 0 ) ) THEN
         INFO = -3
      ENDIF
*
      IF( LLDA .LT. (2*BWL+2*BWU+1) ) THEN
         INFO = -( 6*100 + 6 )
      ENDIF
*
      IF( NB .LE. 0 ) THEN
         INFO = -( 6*100 + 4 )
      ENDIF
*
      BW = BWU+BWL
*
*     Argument checking that is specific to Divide & Conquer routine
*
      IF( NPROW .NE. 1 ) THEN
         INFO = -( 6*100+2 )
      ENDIF
*
      IF( N .GT. NP*NB-MOD( JA-1, NB )) THEN
         INFO = -( 1 )
         CALL PXERBLA( ICTXT,
     $      'PCGBTRF, D&C alg.: only 1 block per proc',
     $      -INFO )
         RETURN
      ENDIF
*
      IF((JA+N-1.GT.NB) .AND. ( NB.LT.(BWL+BWU+1) )) THEN
         INFO = -( 6*100+4 )
         CALL PXERBLA( ICTXT,
     $      'PCGBTRF, D&C alg.: NB too small',
     $      -INFO )
         RETURN
      ENDIF
*
*
*     Check auxiliary storage size
*
      LAF_MIN = (NB+BWU)*(BWL+BWU)+6*(BWL+BWU)*(BWL+2*BWU)
*
      IF( LAF .LT. LAF_MIN ) THEN
         INFO = -9
*        put minimum value of laf into AF( 1 )
         AF( 1 ) = LAF_MIN
         CALL PXERBLA( ICTXT,
     $      'PCGBTRF: auxiliary storage error ',
     $      -INFO )
         RETURN
      ENDIF
*
*     Check worksize
*
      WORK_SIZE_MIN = 1
*
      WORK( 1 ) = WORK_SIZE_MIN
*
      IF( LWORK .LT. WORK_SIZE_MIN ) THEN
         IF( LWORK .NE. -1 ) THEN
         INFO = -11
*        put minimum value of work into work( 1 )
         WORK( 1 ) = WORK_SIZE_MIN
         CALL PXERBLA( ICTXT,
     $      'PCGBTRF: worksize error ',
     $      -INFO )
         ENDIF
         RETURN
      ENDIF
*
*     Pack params and positions into arrays for global consistency check
*
      PARAM_CHECK(  9, 1 ) = DESCA(5)
      PARAM_CHECK(  8, 1 ) = DESCA(4)
      PARAM_CHECK(  7, 1 ) = DESCA(3)
      PARAM_CHECK(  6, 1 ) = DESCA(1)
      PARAM_CHECK(  5, 1 ) = JA
      PARAM_CHECK(  4, 1 ) = BWU
      PARAM_CHECK(  3, 1 ) = BWL
      PARAM_CHECK(  2, 1 ) = N
      PARAM_CHECK(  1, 1 ) = IDUM3
*
      PARAM_CHECK(  9, 2 ) = 605
      PARAM_CHECK(  8, 2 ) = 604
      PARAM_CHECK(  7, 2 ) = 603
      PARAM_CHECK(  6, 2 ) = 601
      PARAM_CHECK(  5, 2 ) = 5
      PARAM_CHECK(  4, 2 ) = 3
      PARAM_CHECK(  3, 2 ) = 2
      PARAM_CHECK(  2, 2 ) = 1
      PARAM_CHECK(  1, 2 ) = 11
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
      CALL GLOBCHK( ICTXT, 9, PARAM_CHECK, 9,
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
         CALL PXERBLA( ICTXT, 'PCGBTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
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
*     ********************************
*     Values reused throughout routine
*
*     User-input value of partition size
*
      PART_SIZE = NB
*
*     Number of columns in each processor
*
      MY_NUM_COLS = NUMROC( N, PART_SIZE, MYCOL, 0, NPCOL )
*
*     Offset in columns to beginning of main partition in each proc
*
      IF ( MYCOL .EQ. 0 ) THEN
        PART_OFFSET = PART_OFFSET+MOD( JA_NEW-1, PART_SIZE )
        MY_NUM_COLS = MY_NUM_COLS - MOD(JA_NEW-1, PART_SIZE )
      ENDIF
*
*     Offset in elements
*
      OFST = PART_OFFSET*LLDA
*
*     Size of main (or odd) partition in each processor
*
      ODD_SIZE = NUMROC( N, PART_SIZE, MYCOL, 0, NPCOL )
*
*
*     Zero out space for fillin
*
      DO 10 I = 1, LAF_MIN
         AF( I ) = CZERO
   10 CONTINUE
*
      DO 9 J = 1, ODD_SIZE
         DO 8 I = 1, BW
            A( I+(J-1)*LLDA ) = CZERO
    8    CONTINUE
    9 CONTINUE
*
*     Begin main code
*
********************************************************************
*     PHASE 1: Local computation phase.
********************************************************************
*
*
*     Transfer triangle B_i of local matrix to next processor
*     for fillin. Overlap the send with the factorization of A_i.
*
      IF (MYCOL .LE. NPCOL-2) THEN
*
*     The last processor does not need to send anything.
*     BIPTR = location of triangle B_i in memory
         BIPTR = (NB-BW)*LLDA + 2*BW+1
*
         CALL CTRSD2D( ICTXT, 'U', 'N',
     $       MIN( BW, BWU+NUMROC( N, NB, MYCOL+1, 0, NPCOL ) ),
     $       BW, A(BIPTR), LLDA-1, 0, MYCOL+1)
*
      ENDIF
*
*     Factor main partition P_i A_i = L_i U_i on each processor
*
*     LBWL, LBWU: lower and upper bandwidth of local solver
*     Note that for MYCOL > 0 one has lower triangular blocks!
*     LM is the number of rows which is usually NB except for
*     MYCOL = 0 where it is BWU less and MYCOL=NPCOL-1 where it
*     is NR+BWU where NR is the number of columns on the last processor
*     Finally APTR is the pointer to the first element of A. As LAPACK
*     has a slightly different matrix format than Scalapack the pointer
*     has to be adjusted on processor MYCOL=0.
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
      IF (LN .GT. 0) THEN
*
         CALL CGBTRF(LM,LN, LBWL,LBWU, A(APTR),LLDA, IPIV, INFO)
*
         IF( INFO.NE.0 ) THEN
            INFO = INFO + NB*MYCOL
            GO TO 90
         END IF
*
         NRHS = BW
         LDB = LLDA-1
*
*     Update the last BW columns of A_i (code modified from CGBTRS)
*
*     Only the eliminations of unknowns > LN-BW have an effect on
*     the last BW columns. Loop over them...
*
            DO 23 J = MAX(LN-BW+1,1), LN
*
               LMJ = MIN( LBWL, LM-J )
               LNJ = MIN( BW, J+BW-LN+APTR-1 )
*
               L = IPIV( J )
*
               JPTR = J-(LN+1)+2*BW+1-LBWL + LN*LLDA
*
               IF( L.NE.J ) THEN
*
*        Element (L,LN+1) is swapped with element (J,LN+1) etc
*        Furthermore, the elements in the same row are LDB=LLDA-1 apart
*        The complicated formulas are to cope with the banded
*          data format:
*
                  LPTR = L-(LN+1)+2*BW+1-LBWL + LN*LLDA
*
                  CALL CSWAP( LNJ, A(LPTR),LDB, A(JPTR), LDB )
*
               ENDIF
*
*              LPTR is the pointer to the beginning of the
*                 coefficients of L
*
               LPTR = BW+1+APTR + (J-1)*LLDA
*
               CALL CGERU(LMJ,LNJ,-CONE, A(LPTR),1, A(JPTR),LDB,
     $           A(JPTR+1),LDB)
   23       CONTINUE
*
      ENDIF
*
*      Compute spike fill-in, L_i F_i = P_i B_{i-1}
*
*      Receive triangle B_{i-1} from previous processor
*
      IF (MYCOL .GT. 0) THEN
*
         CALL CTRRV2D( ICTXT, 'U', 'N', MIN(BW, LM), BW, AF( 1 ),
     $        LM, 0, MYCOL-1)
*
*
*      Permutation and forward elimination (triang. solve)
*
         DO 24 J = 1, LN
*
            LMJ = MIN( LBWL, LM-J )
            L = IPIV( J )
*
            IF( L .NE. J ) THEN
*
                CALL CSWAP(NRHS,  AF(L), LM, AF(J), LM )
            ENDIF
*
            LPTR = BW+1+APTR + (J-1)*LLDA
*
            CALL CGERU( LMJ,NRHS, -CONE, A(LPTR),1,
     $           AF(J), LM, AF(J+1), LM)
*
   24   CONTINUE
*
      ENDIF
*
   90 CONTINUE
*
********************************************************************
*     PHASE 2: Formation and factorization of Reduced System.
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
*     Copy from A and AF into block bidiagonal matrix (tail of AF)
*
*     DBPTR = Pointer to diagonal blocks in A
      DBPTR = BW+1 + LBWU + LN*LLDA
*
      CALL CLAMOV('G',BM,BN, A(DBPTR),LLDA-1,
     $     AF(BBPTR + BW*LDBB),LDBB)
*
*     Zero out any junk entries that were copied
*
      DO 870 J=1, BM
         DO 880 I=J+LBWL, BM-1
            AF( BBPTR+BW*LDBB+(J-1)*LDBB+I ) = CZERO
  880    CONTINUE
  870 CONTINUE
*
      IF (MYCOL .NE. 0) THEN
*
*        ODPTR = Pointer to offdiagonal blocks in A
*
         ODPTR = LM-BM+1
         CALL CLAMOV('G',BM,BW, AF(ODPTR),LM,
     $        AF(BBPTR +2*BW*LDBB),LDBB)
      ENDIF
*
      IF (NPCOL.EQ.1) THEN
*
*        In this case the loop over the levels will not be
*        performed.
         CALL CGETRF( N-LN, N-LN, AF(BBPTR+BW*LDBB), LDBB,
     $        IPIV(LN+1), INFO)
*
      ENDIF
*
*     Loop over levels ... only occurs if npcol > 1
*
*     The two integers NPACT (nu. of active processors) and NPSTR
*     (stride between active processors) are used to control the
*     loop.
*
        NPACT = NPCOL
        NPSTR = 1
*
*       Begin loop over levels
*
  200   IF (NPACT .LE. 1) GOTO 300
*
*         Test if processor is active
*
          IF (MOD(MYCOL,NPSTR) .EQ. 0) THEN
*
*   Send/Receive blocks
*
*
             IF (MOD(MYCOL,2*NPSTR) .EQ. 0) THEN
*
*            This node will potentially do more work later
*
                NEICOL = MYCOL + NPSTR
*
                IF (NEICOL/NPSTR .LT. NPACT-1) THEN
                   BMN = BW
                ELSE IF (NEICOL/NPSTR .EQ. NPACT-1) THEN
                   ODD_N = NUMROC(N, NB, NPCOL-1, 0, NPCOL)
                   BMN = MIN(BW,ODD_N) + BWU
                ELSE
*
*                  Last processor skips to next level
                   GOTO 250
                ENDIF
*
*               BM1 = M for 1st block on proc pair, BM2 2nd block
*
                BM1 = BM
                BM2 = BMN
*
                IF (NEICOL/NPSTR .LE. NPACT-1 )THEN
*
                   CALL CGESD2D( ICTXT, BM, 2*BW, AF(BBPTR+BW*LDBB),
     $                  LDBB, 0, NEICOL )
*
                   CALL CGERV2D( ICTXT, BMN, 2*BW, AF(BBPTR+BM),
     $                  LDBB, 0, NEICOL)
*
                   IF( NPACT .EQ. 2 ) THEN
*
*                     Copy diagonal block to align whole system
*
                      CALL CLAMOV( 'G', BMN, BW, AF( BBPTR+BM ),
     $                  LDBB, AF( BBPTR+2*BW*LDBB+BM ), LDBB )
                   ENDIF
*
                ENDIF
*
             ELSE
*
*               This node stops work after this stage -- an extra copy
*               is required to make the odd and even frontal matrices
*               look identical
*
                NEICOL = MYCOL - NPSTR
*
                IF (NEICOL .EQ. 0) THEN
                   BMN = BW - BWU
                ELSE
                   BMN = BW
                ENDIF
*
                BM1 = BMN
                BM2 = BM
*
                CALL CGESD2D( ICTXT, BM, 2*BW, AF(BBPTR+BW*LDBB),
     $               LDBB, 0, NEICOL )
*
                CALL CLAMOV('G',BM, 2*BW, AF(BBPTR+BW*LDBB),LDBB,
     $               AF(BBPTR+BMN),LDBB)
*
                DO 31 J=BBPTR+2*BW*LDBB, BBPTR+3*BW*LDBB-1, LDBB
                   DO 32 I=0, LDBB-1
                      AF(I+J) = CZERO
   32              CONTINUE
   31           CONTINUE
*
                CALL CGERV2D( ICTXT, BMN, 2*BW, AF(BBPTR+BW*LDBB),
     $               LDBB, 0, NEICOL)
*
                IF( NPACT .EQ. 2 ) THEN
*
*                  Copy diagonal block to align whole system
*
                   CALL CLAMOV( 'G', BM, BW, AF( BBPTR+BMN ),
     $               LDBB, AF( BBPTR+2*BW*LDBB+BMN ), LDBB )
                ENDIF
*
             ENDIF
*
*            LU factorization with partial pivoting
*
             IF (NPACT .NE. 2) THEN
*
                CALL CGETRF(BM+BMN, BW, AF(BBPTR+BW*LDBB), LDBB,
     $               IPIV(LN+1), INFO)
*
*   Backsolve left side
*
                DO 301 J=BBPTR,BBPTR+BW*LDBB-1, LDBB
                   DO 302 I=0, BM1-1
                      AF(I+J) = CZERO
  302              CONTINUE
  301           CONTINUE
*
                CALL CLASWP(BW, AF(BBPTR), LDBB, 1, BW,
     $               IPIV(LN+1), 1)
*
                CALL CTRSM('L','L','N','U', BW, BW, CONE,
     $                AF(BBPTR+BW*LDBB), LDBB, AF(BBPTR), LDBB)
*
*               Use partial factors to update remainder
*
                CALL CGEMM( 'N', 'N', BM+BMN-BW, BW, BW,
     $                -CONE, AF(BBPTR+BW*LDBB+BW), LDBB,
     $                AF( BBPTR ), LDBB, CONE,
     $                AF( BBPTR+BW ), LDBB )
*
*   Backsolve right side
*
                NRHS = BW
*
                CALL CLASWP(NRHS, AF(BBPTR+2*BW*LDBB), LDBB, 1, BW,
     $               IPIV(LN+1), 1)
*
                CALL CTRSM('L','L','N','U', BW, NRHS, CONE,
     $               AF(BBPTR+BW*LDBB), LDBB, AF(BBPTR+2*BW*LDBB), LDBB)
*
*               Use partial factors to update remainder
*
                CALL CGEMM( 'N', 'N', BM+BMN-BW, NRHS, BW,
     $                -CONE, AF(BBPTR+BW*LDBB+BW), LDBB,
     $                AF( BBPTR+2*BW*LDBB ), LDBB, CONE,
     $                AF( BBPTR+2*BW*LDBB+BW ), LDBB )
*
*
*     Test if processor is active in next round
*
                IF (MOD(MYCOL,2*NPSTR) .EQ. 0) THEN
*
*                  Reset BM
*
                   BM = BM1+BM2-BW
*
*                  Local copying in the block bidiagonal area
*
*
                   CALL CLAMOV('G',BM,BW,
     $                  AF(BBPTR+BW),
     $                  LDBB, AF(BBPTR+BW*LDBB), LDBB)
                   CALL CLAMOV('G',BM,BW,
     $                  AF(BBPTR+2*BW*LDBB+BW),
     $                  LDBB, AF(BBPTR+2*BW*LDBB), LDBB)
*
*                  Zero out space that held original copy
*
                   DO 1020 J=0, BW-1
                      DO 1021 I=0, BM-1
                         AF(BBPTR+2*BW*LDBB+BW+J*LDBB+I) = CZERO
 1021                 CONTINUE
 1020              CONTINUE
*
                ENDIF
*
             ELSE
*
*               Factor the final 2 by 2 block matrix
*
                CALL CGETRF(BM+BMN,BM+BMN, AF(BBPTR+BW*LDBB), LDBB,
     $               IPIV(LN+1), INFO)
             ENDIF
*
          ENDIF
*
*        Last processor in an odd-sized NPACT skips to here
*
  250    CONTINUE
*
         NPACT = (NPACT + 1)/2
         NPSTR = NPSTR * 2
         GOTO 200
*
  300 CONTINUE
*     End loop over levels
*
 1000 CONTINUE
*     If error was found in Phase 1, processors jump here.
*
*     Free BLACS space used to hold standard-form grid.
*
      ICTXT = ICTXT_SAVE
      IF( ICTXT.NE.ICTXT_NEW ) THEN
         CALL BLACS_GRIDEXIT( ICTXT_NEW )
      END IF
*
 1234 CONTINUE
*     If this processor did not hold part of the grid it
*        jumps here.
*
*     Restore saved input parameters
*
      ICTXT = ICTXT_SAVE
      NP = NP_SAVE
*
*     Output worksize
*
      WORK( 1 ) = WORK_SIZE_MIN
*
*         Make INFO consistent across processors
*
          CALL IGAMX2D( ICTXT, 'A', ' ', 1, 1, INFO, 1, INFO, INFO,
     $              -1, 0, 0 )
*
          IF( MYCOL.EQ.0 ) THEN
             CALL IGEBS2D( ICTXT, 'A', ' ', 1, 1, INFO, 1 )
          ELSE
             CALL IGEBR2D( ICTXT, 'A', ' ', 1, 1, INFO, 1, 0, 0 )
          ENDIF
*
*
      RETURN
*
*     End of PCGBTRF
*
      END
*
