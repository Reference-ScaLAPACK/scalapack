      SUBROUTINE PDPBTRF( UPLO, N, BW, A, JA, DESCA, AF, LAF, WORK,
     $                    LWORK, INFO )
*
*  -- ScaLAPACK routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            BW, INFO, JA, LAF, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), AF( * ), WORK( * )
*     ..
*
*
*  Purpose
*  =======
*
*  PDPBTRF computes a Cholesky factorization
*  of an N-by-N real banded
*  symmetric positive definite distributed matrix
*  with bandwidth BW: A(1:N, JA:JA+N-1).
*  Reordering is used to increase parallelism in the factorization.
*  This reordering results in factors that are DIFFERENT from those
*  produced by equivalent sequential codes. These factors cannot
*  be used directly by users; however, they can be used in
*  subsequent calls to PDPBTRS to solve linear systems.
*
*  The factorization has the form
*
*          P A(1:N, JA:JA+N-1) P^T = U' U ,  if UPLO = 'U', or
*
*          P A(1:N, JA:JA+N-1) P^T = L L', if UPLO = 'L'
*
*  where U is a banded upper triangular matrix and L is banded
*  lower triangular, and P is a permutation matrix.
*
*  =====================================================================
*
*  Arguments
*  =========
*
*  UPLO    (global input) CHARACTER
*          = 'U':  Upper triangle of A(1:N, JA:JA+N-1) is stored;
*          = 'L':  Lower triangle of A(1:N, JA:JA+N-1) is stored.
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix A(1:N, JA:JA+N-1). N >= 0.
*
*  BW      (global input) INTEGER
*          Number of subdiagonals in L or U. 0 <= BW <= N-1
*
*  A       (local input/local output) DOUBLE PRECISION pointer into
*          local memory to an array with first dimension
*          LLD_A >=(bw+1) (stored in DESCA).
*          On entry, this array contains the local pieces of the
*          N-by-N symmetric banded distributed matrix
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
*  AF      (local output) DOUBLE PRECISION array, dimension LAF.
*          Auxiliary Fillin Space.
*          Fillin is created during the factorization routine
*          PDPBTRF and this is stored in AF. If a linear system
*          is to be solved using PDPBTRS after the factorization
*          routine, AF *must not be altered* after the factorization.
*
*  LAF     (local input) INTEGER
*          Size of user-input Auxiliary Fillin space AF. Must be >=
*          (NB+2*bw)*bw
*          If LAF is not large enough, an error code will be returned
*          and the minimum acceptable size will be returned in AF( 1 )
*
*  WORK    (local workspace/local output)
*          DOUBLE PRECISION temporary workspace. This space may
*          be overwritten in between calls to routines. WORK must be
*          the size given in LWORK.
*          On exit, WORK( 1 ) contains the minimal LWORK.
*
*  LWORK   (local input or global input) INTEGER
*          Size of user-input workspace WORK.
*          If LWORK is too small, the minimal acceptable size will be
*          returned in WORK(1) and an error code is returned. LWORK>=
*          bw*bw
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = K<=NPROCS, the submatrix stored on processor
*                INFO and factored locally was not
*                positive definite,  and
*                the factorization was not completed.
*                If INFO = K>NPROCS, the submatrix stored on processor
*                INFO-NPROCS representing interactions with other
*                processors was not
*                positive definite,
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
*      NB >= 2*BW
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
*         A small (BW* (P-1)) system is formed representing
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
*  Code Developer: Andrew J. Cleary, University of Tennessee.
*    Current address: Lawrence Livermore National Labs.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
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
      INTEGER            COMM_PROC, CSRC, FIRST_PROC, I, ICTXT,
     $                   ICTXT_NEW, ICTXT_SAVE, IDUM1, IDUM3, JA_NEW,
     $                   LAF_MIN, LEVEL_DIST, LLDA, MBW2, MYCOL, MYROW,
     $                   MY_NUM_COLS, NB, NEXT_TRI_SIZE_M,
     $                   NEXT_TRI_SIZE_N, NP, NPCOL, NPROW, NP_SAVE,
     $                   ODD_SIZE, OFST, PART_OFFSET, PART_SIZE,
     $                   PREV_TRI_SIZE_M, PREV_TRI_SIZE_N, RETURN_CODE,
     $                   STORE_N_A, WORK_SIZE_MIN
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA_1XP( 7 ), PARAM_CHECK( 9, 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDEXIT, BLACS_GRIDINFO, DAXPY,
     $                   DESC_CONVERT, DGEMM, DGERV2D, DGESD2D, DLAMOV,
     $                   DLATCPY, DPBTRF, DPOTRF, DSYRK, DTBTRS, DTRMM,
     $                   DTRRV2D, DTRSD2D, DTRSM, DTRTRS, GLOBCHK,
     $                   IGAMX2D, IGEBR2D, IGEBS2D, PXERBLA, RESHAPE
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
      IF( RETURN_CODE.NE.0 ) THEN
         INFO = -( 6*100+2 )
      END IF
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
*     Pre-calculate bw^2
*
      MBW2 = BW*BW
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NP = NPROW*NPCOL
*
*
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         IDUM1 = ICHAR( 'U' )
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         IDUM1 = ICHAR( 'L' )
      ELSE
         INFO = -1
      END IF
*
      IF( LWORK.LT.-1 ) THEN
         INFO = -10
      ELSE IF( LWORK.EQ.-1 ) THEN
         IDUM3 = -1
      ELSE
         IDUM3 = 1
      END IF
*
      IF( N.LT.0 ) THEN
         INFO = -2
      END IF
*
      IF( N+JA-1.GT.STORE_N_A ) THEN
         INFO = -( 6*100+6 )
      END IF
*
      IF( ( BW.GT.N-1 ) .OR. ( BW.LT.0 ) ) THEN
         INFO = -3
      END IF
*
      IF( LLDA.LT.( BW+1 ) ) THEN
         INFO = -( 6*100+6 )
      END IF
*
      IF( NB.LE.0 ) THEN
         INFO = -( 6*100+4 )
      END IF
*
*     Argument checking that is specific to Divide & Conquer routine
*
      IF( NPROW.NE.1 ) THEN
         INFO = -( 6*100+2 )
      END IF
*
      IF( N.GT.NP*NB-MOD( JA-1, NB ) ) THEN
         INFO = -( 2 )
         CALL PXERBLA( ICTXT, 'PDPBTRF, D&C alg.: only 1 block per proc'
     $                 , -INFO )
         RETURN
      END IF
*
      IF( ( JA+N-1.GT.NB ) .AND. ( NB.LT.2*BW ) ) THEN
         INFO = -( 6*100+4 )
         CALL PXERBLA( ICTXT, 'PDPBTRF, D&C alg.: NB too small', -INFO )
         RETURN
      END IF
*
*
*     Check auxiliary storage size
*
      LAF_MIN = ( NB+2*BW )*BW
*
      IF( LAF.LT.LAF_MIN ) THEN
         INFO = -8
*        put minimum value of laf into AF( 1 )
         AF( 1 ) = LAF_MIN
         CALL PXERBLA( ICTXT, 'PDPBTRF: auxiliary storage error ',
     $                 -INFO )
         RETURN
      END IF
*
*     Check worksize
*
      WORK_SIZE_MIN = BW*BW
*
      WORK( 1 ) = WORK_SIZE_MIN
*
      IF( LWORK.LT.WORK_SIZE_MIN ) THEN
         IF( LWORK.NE.-1 ) THEN
            INFO = -10
            CALL PXERBLA( ICTXT, 'PDPBTRF: worksize error ', -INFO )
         END IF
         RETURN
      END IF
*
*     Pack params and positions into arrays for global consistency check
*
      PARAM_CHECK( 9, 1 ) = DESCA( 5 )
      PARAM_CHECK( 8, 1 ) = DESCA( 4 )
      PARAM_CHECK( 7, 1 ) = DESCA( 3 )
      PARAM_CHECK( 6, 1 ) = DESCA( 1 )
      PARAM_CHECK( 5, 1 ) = JA
      PARAM_CHECK( 4, 1 ) = BW
      PARAM_CHECK( 3, 1 ) = N
      PARAM_CHECK( 2, 1 ) = IDUM3
      PARAM_CHECK( 1, 1 ) = IDUM1
*
      PARAM_CHECK( 9, 2 ) = 605
      PARAM_CHECK( 8, 2 ) = 604
      PARAM_CHECK( 7, 2 ) = 603
      PARAM_CHECK( 6, 2 ) = 601
      PARAM_CHECK( 5, 2 ) = 5
      PARAM_CHECK( 4, 2 ) = 3
      PARAM_CHECK( 3, 2 ) = 2
      PARAM_CHECK( 2, 2 ) = 10
      PARAM_CHECK( 1, 2 ) = 1
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
         INFO = -INFO*DESCMULT
      END IF
*
*     Check consistency across processors
*
      CALL GLOBCHK( ICTXT, 9, PARAM_CHECK, 9, PARAM_CHECK( 1, 3 ),
     $              INFO )
*
*     Prepare output: set info = 0 if no error, and divide by DESCMULT
*     if error is not in a descriptor entry.
*
      IF( INFO.EQ.BIGNUM ) THEN
         INFO = 0
      ELSE IF( MOD( INFO, DESCMULT ).EQ.0 ) THEN
         INFO = -INFO / DESCMULT
      ELSE
         INFO = -INFO
      END IF
*
      IF( INFO.LT.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDPBTRF', -INFO )
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
      PART_OFFSET = NB*( ( JA-1 ) / ( NPCOL*NB ) )
*
      IF( ( MYCOL-CSRC ).LT.( JA-PART_OFFSET-1 ) / NB ) THEN
         PART_OFFSET = PART_OFFSET + NB
      END IF
*
      IF( MYCOL.LT.CSRC ) THEN
         PART_OFFSET = PART_OFFSET - NB
      END IF
*
*     Form a new BLACS grid (the "standard form" grid) with only procs
*        holding part of the matrix, of size 1xNP where NP is adjusted,
*        starting at csrc=0, with JA modified to reflect dropped procs.
*
*     First processor to hold part of the matrix:
*
      FIRST_PROC = MOD( ( JA-1 ) / NB+CSRC, NPCOL )
*
*     Calculate new JA one while dropping off unused processors.
*
      JA_NEW = MOD( JA-1, NB ) + 1
*
*     Save and compute new value of NP
*
      NP_SAVE = NP
      NP = ( JA_NEW+N-2 ) / NB + 1
*
*     Call utility routine that forms "standard-form" grid
*
      CALL RESHAPE( ICTXT, INT_ONE, ICTXT_NEW, INT_ONE, FIRST_PROC,
     $              INT_ONE, NP )
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
      IF( MYROW.LT.0 ) THEN
         GO TO 120
      END IF
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
      IF( MYCOL.EQ.0 ) THEN
         PART_OFFSET = PART_OFFSET + MOD( JA_NEW-1, PART_SIZE )
         MY_NUM_COLS = MY_NUM_COLS - MOD( JA_NEW-1, PART_SIZE )
      END IF
*
*     Offset in elements
*
      OFST = PART_OFFSET*LLDA
*
*     Size of main (or odd) partition in each processor
*
      ODD_SIZE = MY_NUM_COLS
      IF( MYCOL.LT.NP-1 ) THEN
         ODD_SIZE = ODD_SIZE - BW
      END IF
*
*
*       Zero out space for fillin
*
      DO 10 I = 1, LAF_MIN
         AF( I ) = ZERO
   10 CONTINUE
*
*       Zero out space for work
*
      DO 20 I = 1, WORK_SIZE_MIN
         WORK( I ) = ZERO
   20 CONTINUE
*
*     Begin main code
*
      IF( LSAME( UPLO, 'L' ) ) THEN
*
********************************************************************
*       PHASE 1: Local computation phase.
********************************************************************
*
*
*       Sizes of the extra triangles communicated bewtween processors
*
         IF( MYCOL.GT.0 ) THEN
            PREV_TRI_SIZE_M = MIN( BW, NUMROC( N, PART_SIZE, MYCOL, 0,
     $                        NPCOL ) )
            PREV_TRI_SIZE_N = MIN( BW, NUMROC( N, PART_SIZE, MYCOL-1, 0,
     $                        NPCOL ) )
         END IF
*
         IF( MYCOL.LT.NPCOL-1 ) THEN
            NEXT_TRI_SIZE_M = MIN( BW, NUMROC( N, PART_SIZE, MYCOL+1, 0,
     $                        NPCOL ) )
            NEXT_TRI_SIZE_N = MIN( BW, NUMROC( N, PART_SIZE, MYCOL, 0,
     $                        NPCOL ) )
         END IF
*
         IF( MYCOL.LT.NP-1 ) THEN
*         Transfer last triangle D_i of local matrix to next processor
*         which needs it to calculate fillin due to factorization of
*         its main (odd) block A_i.
*         Overlap the send with the factorization of A_i.
*
            CALL DTRSD2D( ICTXT, 'U', 'N', NEXT_TRI_SIZE_M,
     $                    NEXT_TRI_SIZE_N, A( OFST+ODD_SIZE*LLDA+( BW+
     $                    1 ) ), LLDA-1, 0, MYCOL+1 )
*
         END IF
*
*
*       Factor main partition A_i = L_i {L_i}^T in each processor
*
         CALL DPBTRF( UPLO, ODD_SIZE, BW, A( OFST+1 ), LLDA, INFO )
*
         IF( INFO.NE.0 ) THEN
            INFO = MYCOL + 1
            GO TO 30
         END IF
*
*
         IF( MYCOL.LT.NP-1 ) THEN
*         Apply factorization to odd-even connection block B_i
*
*         transpose the connection block in preparation.
*
            CALL DLATCPY( 'U', BW, BW, A( ( OFST+( BW+1 )+( ODD_SIZE-
     $                    BW )*LLDA ) ), LLDA-1,
     $                    AF( ODD_SIZE*BW+2*MBW2+1+BW-BW ), BW )
*
*         Perform the triangular system solve {L_i}{{B'}_i}^T = {B_i}^T
*
            CALL DTRTRS( 'L', 'N', 'N', BW, BW,
     $                   A( OFST+1+( ODD_SIZE-BW )*LLDA ), LLDA-1,
     $                   AF( ODD_SIZE*BW+2*MBW2+1 ), BW, INFO )
*
*
*         transpose resulting block to its location
*           in main storage.
*
            CALL DLATCPY( 'L', BW, BW, AF( ODD_SIZE*BW+2*MBW2+1+BW-BW ),
     $                    BW, A( ( OFST+( BW+1 )+( ODD_SIZE-BW )*
     $                    LLDA ) ), LLDA-1 )
*
*
*         Compute contribution to diagonal block(s) of reduced system.
*          {C'}_i = {C_i}-{{B'}_i}{{B'}_i}^T
*
*         The following method uses more flops than necessary but
*           does not necessitate the writing of a new BLAS routine.
*
*
            CALL DSYRK( UPLO, 'T', BW, BW, -ONE,
     $                  AF( ODD_SIZE*BW+2*MBW2+1 ), BW, ONE,
     $                  A( OFST+1+ODD_SIZE*LLDA ), LLDA-1 )
*
         END IF
*       End of "if ( MYCOL .lt. NP-1 )..." loop
*
*
   30    CONTINUE
*       If the processor could not locally factor, it jumps here.
*
         IF( MYCOL.NE.0 ) THEN
*         Discard temporary matrix stored beginning in
*           AF( (odd_size+2*bw)*bw+1 ) and use for
*           off_diagonal block of reduced system.
*
*         Receive previously transmitted matrix section, which forms
*         the right-hand-side for the triangular solve that calculates
*         the "spike" fillin.
*
*
            CALL DTRRV2D( ICTXT, 'U', 'N', PREV_TRI_SIZE_M,
     $                    PREV_TRI_SIZE_N, AF( 1 ), ODD_SIZE, 0,
     $                    MYCOL-1 )
*
            IF( INFO.EQ.0 ) THEN
*
*         Calculate the "spike" fillin, ${L_i} {{G}_i}^T = {D_i}$ .
*
               CALL DTBTRS( 'L', 'N', 'N', ODD_SIZE, BW, BW,
     $                      A( OFST+1 ), LLDA, AF( 1 ), ODD_SIZE, INFO )
*
*
*         Calculate the update block for previous proc, E_i = G_i{G_i}^T
*
               CALL DSYRK( 'L', 'T', BW, ODD_SIZE, -ONE, AF( 1 ),
     $                     ODD_SIZE, ZERO, AF( 1+( ODD_SIZE+2*BW )*BW ),
     $                     BW )
*
*
*         Initiate send of E_i to previous processor to overlap
*           with next computation.
*
               CALL DGESD2D( ICTXT, BW, BW, AF( ODD_SIZE*BW+2*MBW2+1 ),
     $                       BW, 0, MYCOL-1 )
*
*
               IF( MYCOL.LT.NP-1 ) THEN
*
*           Calculate off-diagonal block(s) of reduced system.
*           Note: for ease of use in solution of reduced system, store
*           L's off-diagonal block in transpose form.
*           {F_i}^T =  {H_i}{{B'}_i}^T
*
*           Copy matrix H_i (the last bw cols of G_i) to AF storage
*             as per requirements of BLAS routine DTRMM.
*             Since we have G_i^T stored, transpose
*             H_i^T to H_i.
*
                  CALL DLATCPY( 'N', BW, BW, AF( ODD_SIZE-BW+1 ),
     $                          ODD_SIZE, AF( ( ODD_SIZE )*BW+1 ), BW )
*
                  CALL DTRMM( 'R', 'U', 'T', 'N', BW, BW, -ONE,
     $                        A( ( OFST+( BW+1 )+( ODD_SIZE-BW )*
     $                        LLDA ) ), LLDA-1, AF( ( ODD_SIZE )*BW+1 ),
     $                        BW )
*
*
               END IF
*
            END IF
*       End of "if ( MYCOL .ne. 0 )..."
*
         END IF
*       End of "if (info.eq.0) then"
*
*
*       Check to make sure no processors have found errors
*
         CALL IGAMX2D( ICTXT, 'A', ' ', 1, 1, INFO, 1, INFO, INFO, -1,
     $                 0, 0 )
*
         IF( MYCOL.EQ.0 ) THEN
            CALL IGEBS2D( ICTXT, 'A', ' ', 1, 1, INFO, 1 )
         ELSE
            CALL IGEBR2D( ICTXT, 'A', ' ', 1, 1, INFO, 1, 0, 0 )
         END IF
*
         IF( INFO.NE.0 ) THEN
            GO TO 110
         END IF
*       No errors found, continue
*
*
********************************************************************
*       PHASE 2: Formation and factorization of Reduced System.
********************************************************************
*
*       Gather up local sections of reduced system
*
*
*     The last processor does not participate in the factorization of
*       the reduced system, having sent its E_i already.
         IF( MYCOL.EQ.NPCOL-1 ) THEN
            GO TO 60
         END IF
*
*       Initiate send of off-diag block(s) to overlap with next part.
*       Off-diagonal block needed on neighboring processor to start
*       algorithm.
*
         IF( ( MOD( MYCOL+1, 2 ).EQ.0 ) .AND. ( MYCOL.GT.0 ) ) THEN
*
            CALL DGESD2D( ICTXT, BW, BW, AF( ODD_SIZE*BW+1 ), BW, 0,
     $                    MYCOL-1 )
*
         END IF
*
*       Copy last diagonal block into AF storage for subsequent
*         operations.
*
         CALL DLAMOV( 'N', BW, BW, A( OFST+ODD_SIZE*LLDA+1 ), LLDA-1,
     $                AF( ODD_SIZE*BW+MBW2+1 ), BW )
*
*       Receive cont. to diagonal block that is stored on this proc.
*
         IF( MYCOL.LT.NPCOL-1 ) THEN
*
            CALL DGERV2D( ICTXT, BW, BW, AF( ODD_SIZE*BW+2*MBW2+1 ), BW,
     $                    0, MYCOL+1 )
*
*          Add contribution to diagonal block
*
            CALL DAXPY( MBW2, ONE, AF( ODD_SIZE*BW+2*MBW2+1 ), 1,
     $                  AF( ODD_SIZE*BW+MBW2+1 ), 1 )
*
         END IF
*
*
*       *************************************
*       Modification Loop
*
*       The distance for sending and receiving for each level starts
*         at 1 for the first level.
         LEVEL_DIST = 1
*
*       Do until this proc is needed to modify other procs' equations
*
   40    CONTINUE
         IF( MOD( ( MYCOL+1 ) / LEVEL_DIST, 2 ).NE.0 )
     $      GO TO 50
*
*         Receive and add contribution to diagonal block from the left
*
         IF( MYCOL-LEVEL_DIST.GE.0 ) THEN
            CALL DGERV2D( ICTXT, BW, BW, WORK( 1 ), BW, 0,
     $                    MYCOL-LEVEL_DIST )
*
            CALL DAXPY( MBW2, ONE, WORK( 1 ), 1,
     $                  AF( ODD_SIZE*BW+MBW2+1 ), 1 )
*
         END IF
*
*         Receive and add contribution to diagonal block from the right
*
         IF( MYCOL+LEVEL_DIST.LT.NPCOL-1 ) THEN
            CALL DGERV2D( ICTXT, BW, BW, WORK( 1 ), BW, 0,
     $                    MYCOL+LEVEL_DIST )
*
            CALL DAXPY( MBW2, ONE, WORK( 1 ), 1,
     $                  AF( ODD_SIZE*BW+MBW2+1 ), 1 )
*
         END IF
*
         LEVEL_DIST = LEVEL_DIST*2
*
         GO TO 40
   50    CONTINUE
*       [End of GOTO Loop]
*
*
*       *********************************
*       Calculate and use this proc's blocks to modify other procs'...
*
*       Factor diagonal block
*
         CALL DPOTRF( 'L', BW, AF( ODD_SIZE*BW+MBW2+1 ), BW, INFO )
*
         IF( INFO.NE.0 ) THEN
            INFO = NPCOL + MYCOL
         END IF
*
*       ****************************************************************
*       Receive offdiagonal block from processor to right.
*         If this is the first group of processors, the receive comes
*         from a different processor than otherwise.
*
         IF( LEVEL_DIST.EQ.1 ) THEN
            COMM_PROC = MYCOL + 1
*
*           Move block into place that it will be expected to be for
*             calcs.
*
            CALL DLAMOV( 'N', BW, BW, AF( ODD_SIZE*BW+1 ), BW,
     $                   AF( ODD_SIZE*BW+2*MBW2+1 ), BW )
*
         ELSE
            COMM_PROC = MYCOL + LEVEL_DIST / 2
         END IF
*
         IF( MYCOL / LEVEL_DIST.LE.( NPCOL-1 ) / LEVEL_DIST-2 ) THEN
*
            CALL DGERV2D( ICTXT, BW, BW, AF( ODD_SIZE*BW+1 ), BW, 0,
     $                    COMM_PROC )
*
            IF( INFO.EQ.0 ) THEN
*
*
*         Modify upper off_diagonal block with diagonal block
*
*
               CALL DTRSM( 'L', 'L', 'N', 'N', BW, BW, ONE,
     $                     AF( ODD_SIZE*BW+MBW2+1 ), BW,
     $                     AF( ODD_SIZE*BW+1 ), BW )
*
            END IF
*         End of "if ( info.eq.0 ) then"
*
*         Calculate contribution from this block to next diagonal block
*
            CALL DSYRK( 'L', 'T', BW, BW, -ONE, AF( ( ODD_SIZE )*BW+1 ),
     $                  BW, ZERO, WORK( 1 ), BW )
*
*         Send contribution to diagonal block's owning processor.
*
            CALL DGESD2D( ICTXT, BW, BW, WORK( 1 ), BW, 0,
     $                    MYCOL+LEVEL_DIST )
*
         END IF
*       End of "if( mycol/level_dist .le. (npcol-1)/level_dist-2 )..."
*
*
*       ****************************************************************
*       Receive off_diagonal block from left and use to finish with this
*         processor.
*
         IF( ( MYCOL / LEVEL_DIST.GT.0 ) .AND.
     $       ( MYCOL / LEVEL_DIST.LE.( NPCOL-1 ) / LEVEL_DIST-1 ) ) THEN
*
            IF( LEVEL_DIST.GT.1 ) THEN
*
*           Receive offdiagonal block(s) from proc level_dist/2 to the
*           left
*
               CALL DGERV2D( ICTXT, BW, BW, AF( ODD_SIZE*BW+2*MBW2+1 ),
     $                       BW, 0, MYCOL-LEVEL_DIST / 2 )
*
            END IF
*
*
            IF( INFO.EQ.0 ) THEN
*
*         Use diagonal block(s) to modify this offdiagonal block
*
               CALL DTRSM( 'R', 'L', 'T', 'N', BW, BW, ONE,
     $                     AF( ODD_SIZE*BW+MBW2+1 ), BW,
     $                     AF( ODD_SIZE*BW+2*MBW2+1 ), BW )
*
            END IF
*         End of "if( info.eq.0 ) then"
*
*         Use offdiag block(s) to calculate modification to diag block
*           of processor to the left
*
            CALL DSYRK( 'L', 'N', BW, BW, -ONE,
     $                  AF( ( ODD_SIZE+2*BW )*BW+1 ), BW, ZERO,
     $                  WORK( 1 ), BW )
*
*         Send contribution to diagonal block's owning processor.
*
            CALL DGESD2D( ICTXT, BW, BW, WORK( 1 ), BW, 0,
     $                    MYCOL-LEVEL_DIST )
*
*         *******************************************************
*
            IF( MYCOL / LEVEL_DIST.LE.( NPCOL-1 ) / LEVEL_DIST-2 ) THEN
*
*           Decide which processor offdiagonal block(s) goes to
*
               IF( ( MOD( MYCOL / ( 2*LEVEL_DIST ), 2 ) ).EQ.0 ) THEN
                  COMM_PROC = MYCOL + LEVEL_DIST
               ELSE
                  COMM_PROC = MYCOL - LEVEL_DIST
               END IF
*
*           Use offdiagonal blocks to calculate offdiag
*             block to send to neighboring processor. Depending
*             on circumstances, may need to transpose the matrix.
*
               CALL DGEMM( 'N', 'N', BW, BW, BW, -ONE,
     $                     AF( ODD_SIZE*BW+2*MBW2+1 ), BW,
     $                     AF( ODD_SIZE*BW+1 ), BW, ZERO, WORK( 1 ),
     $                     BW )
*
*           Send contribution to offdiagonal block's owning processor.
*
               CALL DGESD2D( ICTXT, BW, BW, WORK( 1 ), BW, 0,
     $                       COMM_PROC )
*
            END IF
*
         END IF
*       End of "if( mycol/level_dist.le. (npcol-1)/level_dist -1 )..."
*
   60    CONTINUE
*
      ELSE
*
*       CASE UPLO = 'U'
*
********************************************************************
*       PHASE 1: Local computation phase.
********************************************************************
*
*
*       Sizes of the extra triangles communicated bewtween processors
*
         IF( MYCOL.GT.0 ) THEN
            PREV_TRI_SIZE_M = MIN( BW, NUMROC( N, PART_SIZE, MYCOL, 0,
     $                        NPCOL ) )
            PREV_TRI_SIZE_N = MIN( BW, NUMROC( N, PART_SIZE, MYCOL-1, 0,
     $                        NPCOL ) )
         END IF
*
         IF( MYCOL.LT.NPCOL-1 ) THEN
            NEXT_TRI_SIZE_M = MIN( BW, NUMROC( N, PART_SIZE, MYCOL+1, 0,
     $                        NPCOL ) )
            NEXT_TRI_SIZE_N = MIN( BW, NUMROC( N, PART_SIZE, MYCOL, 0,
     $                        NPCOL ) )
         END IF
*
*
*
*       Factor main partition A_i^T = U_i {U_i}^T in each processor
*
         CALL DPBTRF( UPLO, ODD_SIZE, BW, A( OFST+1 ), LLDA, INFO )
*
         IF( INFO.NE.0 ) THEN
            INFO = MYCOL + 1
            GO TO 70
         END IF
*
*
         IF( MYCOL.LT.NP-1 ) THEN
*         Apply factorization to odd-even connection block B_i
*
*         Move the connection block in preparation.
*
            CALL DLAMOV( 'L', BW, BW, A( ( OFST+1+ODD_SIZE*LLDA ) ),
     $                   LLDA-1, AF( ODD_SIZE*BW+2*MBW2+1+BW-BW ), BW )
*
*
*         Perform the triangular solve {L_i}{{B'}_i}^T = {B_i}^T
*
            CALL DTRTRS( 'U', 'T', 'N', BW, BW,
     $                   A( OFST+BW+1+( ODD_SIZE-BW )*LLDA ), LLDA-1,
     $                   AF( ODD_SIZE*BW+2*MBW2+1 ), BW, INFO )
*
*         Move the resulting block back to its location in main storage.
*
            CALL DLAMOV( 'L', BW, BW, AF( ODD_SIZE*BW+2*MBW2+1+BW-BW ),
     $                   BW, A( ( OFST+1+ODD_SIZE*LLDA ) ), LLDA-1 )
*
*
*         Compute contribution to diagonal block(s) of reduced system.
*          {C'}_i^T = {C_i}^T-{{B'}_i}^T{{B'}_i}
*
*         The following method uses more flops than necessary but
*           does not necessitate the writing of a new BLAS routine.
*
*
            CALL DSYRK( UPLO, 'T', BW, BW, -ONE,
     $                  AF( ODD_SIZE*BW+2*MBW2+1 ), BW, ONE,
     $                  A( OFST+BW+1+ODD_SIZE*LLDA ), LLDA-1 )
*
         END IF
*       End of "if ( MYCOL .lt. NP-1 )..." loop
*
*
   70    CONTINUE
*       If the processor could not locally factor, it jumps here.
*
         IF( MYCOL.NE.0 ) THEN
*         Discard temporary matrix stored beginning in
*           AF( (odd_size+2*bw)*bw+1 ) and use for
*           off_diagonal block of reduced system.
*
*         Calculate the "spike" fillin, ${L_i} {{G}_i}^T = {D_i}$ .
*
*
*         Copy D block into AF storage for solve.
*
            CALL DLATCPY( 'L', PREV_TRI_SIZE_N, PREV_TRI_SIZE_M,
     $                    A( OFST+1 ), LLDA-1, AF( 1 ), ODD_SIZE )
*
            IF( INFO.EQ.0 ) THEN
*
               CALL DTBTRS( 'U', 'T', 'N', ODD_SIZE, BW, BW,
     $                      A( OFST+1 ), LLDA, AF( 1 ), ODD_SIZE, INFO )
*
*
*         Calculate the update block for previous proc, E_i = G_i{G_i}^T
*
               CALL DSYRK( 'L', 'T', BW, ODD_SIZE, -ONE, AF( 1 ),
     $                     ODD_SIZE, ZERO, AF( 1+( ODD_SIZE+2*BW )*BW ),
     $                     BW )
*
*
*         Initiate send of E_i to previous processor to overlap
*           with next computation.
*
               CALL DGESD2D( ICTXT, BW, BW, AF( ODD_SIZE*BW+2*MBW2+1 ),
     $                       BW, 0, MYCOL-1 )
*
*
               IF( MYCOL.LT.NP-1 ) THEN
*
*           Calculate off-diagonal block(s) of reduced system.
*           Note: for ease of use in solution of reduced system, store
*           L's off-diagonal block in transpose form.
*           {F_i}^T =  {H_i}{{B'}_i}^T
*
*           Copy matrix H_i (the last bw cols of G_i) to AF storage
*             as per requirements of BLAS routine DTRMM.
*             Since we have G_i^T stored, transpose
*             H_i^T to H_i.
*
                  CALL DLATCPY( 'N', BW, BW, AF( ODD_SIZE-BW+1 ),
     $                          ODD_SIZE, AF( ( ODD_SIZE )*BW+1 ), BW )
*
                  CALL DTRMM( 'R', 'L', 'N', 'N', BW, BW, -ONE,
     $                        A( ( OFST+1+ODD_SIZE*LLDA ) ), LLDA-1,
     $                        AF( ( ODD_SIZE )*BW+1 ), BW )
*
               END IF
*
            END IF
*       End of "if ( MYCOL .ne. 0 )..."
*
         END IF
*       End of "if (info.eq.0) then"
*
*
*       Check to make sure no processors have found errors
*
         CALL IGAMX2D( ICTXT, 'A', ' ', 1, 1, INFO, 1, INFO, INFO, -1,
     $                 0, 0 )
*
         IF( MYCOL.EQ.0 ) THEN
            CALL IGEBS2D( ICTXT, 'A', ' ', 1, 1, INFO, 1 )
         ELSE
            CALL IGEBR2D( ICTXT, 'A', ' ', 1, 1, INFO, 1, 0, 0 )
         END IF
*
         IF( INFO.NE.0 ) THEN
            GO TO 110
         END IF
*       No errors found, continue
*
*
********************************************************************
*       PHASE 2: Formation and factorization of Reduced System.
********************************************************************
*
*       Gather up local sections of reduced system
*
*
*     The last processor does not participate in the factorization of
*       the reduced system, having sent its E_i already.
         IF( MYCOL.EQ.NPCOL-1 ) THEN
            GO TO 100
         END IF
*
*       Initiate send of off-diag block(s) to overlap with next part.
*       Off-diagonal block needed on neighboring processor to start
*       algorithm.
*
         IF( ( MOD( MYCOL+1, 2 ).EQ.0 ) .AND. ( MYCOL.GT.0 ) ) THEN
*
            CALL DGESD2D( ICTXT, BW, BW, AF( ODD_SIZE*BW+1 ), BW, 0,
     $                    MYCOL-1 )
*
         END IF
*
*       Transpose last diagonal block into AF storage for subsequent
*         operations.
*
         CALL DLATCPY( 'U', BW, BW, A( OFST+ODD_SIZE*LLDA+1+BW ),
     $                 LLDA-1, AF( ODD_SIZE*BW+MBW2+1 ), BW )
*
*       Receive cont. to diagonal block that is stored on this proc.
*
         IF( MYCOL.LT.NPCOL-1 ) THEN
*
            CALL DGERV2D( ICTXT, BW, BW, AF( ODD_SIZE*BW+2*MBW2+1 ), BW,
     $                    0, MYCOL+1 )
*
*          Add contribution to diagonal block
*
            CALL DAXPY( MBW2, ONE, AF( ODD_SIZE*BW+2*MBW2+1 ), 1,
     $                  AF( ODD_SIZE*BW+MBW2+1 ), 1 )
*
         END IF
*
*
*       *************************************
*       Modification Loop
*
*       The distance for sending and receiving for each level starts
*         at 1 for the first level.
         LEVEL_DIST = 1
*
*       Do until this proc is needed to modify other procs' equations
*
   80    CONTINUE
         IF( MOD( ( MYCOL+1 ) / LEVEL_DIST, 2 ).NE.0 )
     $      GO TO 90
*
*         Receive and add contribution to diagonal block from the left
*
         IF( MYCOL-LEVEL_DIST.GE.0 ) THEN
            CALL DGERV2D( ICTXT, BW, BW, WORK( 1 ), BW, 0,
     $                    MYCOL-LEVEL_DIST )
*
            CALL DAXPY( MBW2, ONE, WORK( 1 ), 1,
     $                  AF( ODD_SIZE*BW+MBW2+1 ), 1 )
*
         END IF
*
*         Receive and add contribution to diagonal block from the right
*
         IF( MYCOL+LEVEL_DIST.LT.NPCOL-1 ) THEN
            CALL DGERV2D( ICTXT, BW, BW, WORK( 1 ), BW, 0,
     $                    MYCOL+LEVEL_DIST )
*
            CALL DAXPY( MBW2, ONE, WORK( 1 ), 1,
     $                  AF( ODD_SIZE*BW+MBW2+1 ), 1 )
*
         END IF
*
         LEVEL_DIST = LEVEL_DIST*2
*
         GO TO 80
   90    CONTINUE
*       [End of GOTO Loop]
*
*
*       *********************************
*       Calculate and use this proc's blocks to modify other procs'...
*
*       Factor diagonal block
*
         CALL DPOTRF( 'L', BW, AF( ODD_SIZE*BW+MBW2+1 ), BW, INFO )
*
         IF( INFO.NE.0 ) THEN
            INFO = NPCOL + MYCOL
         END IF
*
*       ****************************************************************
*       Receive offdiagonal block from processor to right.
*         If this is the first group of processors, the receive comes
*         from a different processor than otherwise.
*
         IF( LEVEL_DIST.EQ.1 ) THEN
            COMM_PROC = MYCOL + 1
*
*           Move block into place that it will be expected to be for
*             calcs.
*
            CALL DLAMOV( 'N', BW, BW, AF( ODD_SIZE*BW+1 ), BW,
     $                   AF( ODD_SIZE*BW+2*MBW2+1 ), BW )
*
         ELSE
            COMM_PROC = MYCOL + LEVEL_DIST / 2
         END IF
*
         IF( MYCOL / LEVEL_DIST.LE.( NPCOL-1 ) / LEVEL_DIST-2 ) THEN
*
            CALL DGERV2D( ICTXT, BW, BW, AF( ODD_SIZE*BW+1 ), BW, 0,
     $                    COMM_PROC )
*
            IF( INFO.EQ.0 ) THEN
*
*
*         Modify upper off_diagonal block with diagonal block
*
*
               CALL DTRSM( 'L', 'L', 'N', 'N', BW, BW, ONE,
     $                     AF( ODD_SIZE*BW+MBW2+1 ), BW,
     $                     AF( ODD_SIZE*BW+1 ), BW )
*
            END IF
*         End of "if ( info.eq.0 ) then"
*
*         Calculate contribution from this block to next diagonal block
*
            CALL DSYRK( 'L', 'T', BW, BW, -ONE, AF( ( ODD_SIZE )*BW+1 ),
     $                  BW, ZERO, WORK( 1 ), BW )
*
*         Send contribution to diagonal block's owning processor.
*
            CALL DGESD2D( ICTXT, BW, BW, WORK( 1 ), BW, 0,
     $                    MYCOL+LEVEL_DIST )
*
         END IF
*       End of "if( mycol/level_dist .le. (npcol-1)/level_dist-2 )..."
*
*
*       ****************************************************************
*       Receive off_diagonal block from left and use to finish with this
*         processor.
*
         IF( ( MYCOL / LEVEL_DIST.GT.0 ) .AND.
     $       ( MYCOL / LEVEL_DIST.LE.( NPCOL-1 ) / LEVEL_DIST-1 ) ) THEN
*
            IF( LEVEL_DIST.GT.1 ) THEN
*
*           Receive offdiagonal block(s) from proc level_dist/2 to the
*           left
*
               CALL DGERV2D( ICTXT, BW, BW, AF( ODD_SIZE*BW+2*MBW2+1 ),
     $                       BW, 0, MYCOL-LEVEL_DIST / 2 )
*
            END IF
*
*
            IF( INFO.EQ.0 ) THEN
*
*         Use diagonal block(s) to modify this offdiagonal block
*
               CALL DTRSM( 'R', 'L', 'T', 'N', BW, BW, ONE,
     $                     AF( ODD_SIZE*BW+MBW2+1 ), BW,
     $                     AF( ODD_SIZE*BW+2*MBW2+1 ), BW )
*
            END IF
*         End of "if( info.eq.0 ) then"
*
*         Use offdiag block(s) to calculate modification to diag block
*           of processor to the left
*
            CALL DSYRK( 'L', 'N', BW, BW, -ONE,
     $                  AF( ( ODD_SIZE+2*BW )*BW+1 ), BW, ZERO,
     $                  WORK( 1 ), BW )
*
*         Send contribution to diagonal block's owning processor.
*
            CALL DGESD2D( ICTXT, BW, BW, WORK( 1 ), BW, 0,
     $                    MYCOL-LEVEL_DIST )
*
*         *******************************************************
*
            IF( MYCOL / LEVEL_DIST.LE.( NPCOL-1 ) / LEVEL_DIST-2 ) THEN
*
*           Decide which processor offdiagonal block(s) goes to
*
               IF( ( MOD( MYCOL / ( 2*LEVEL_DIST ), 2 ) ).EQ.0 ) THEN
                  COMM_PROC = MYCOL + LEVEL_DIST
               ELSE
                  COMM_PROC = MYCOL - LEVEL_DIST
               END IF
*
*           Use offdiagonal blocks to calculate offdiag
*             block to send to neighboring processor. Depending
*             on circumstances, may need to transpose the matrix.
*
               CALL DGEMM( 'N', 'N', BW, BW, BW, -ONE,
     $                     AF( ODD_SIZE*BW+2*MBW2+1 ), BW,
     $                     AF( ODD_SIZE*BW+1 ), BW, ZERO, WORK( 1 ),
     $                     BW )
*
*           Send contribution to offdiagonal block's owning processor.
*
               CALL DGESD2D( ICTXT, BW, BW, WORK( 1 ), BW, 0,
     $                       COMM_PROC )
*
            END IF
*
         END IF
*       End of "if( mycol/level_dist.le. (npcol-1)/level_dist -1 )..."
*
  100    CONTINUE
*
      END IF
*
  110 CONTINUE
*
*
*     Free BLACS space used to hold standard-form grid.
*
      IF( ICTXT_SAVE.NE.ICTXT_NEW ) THEN
         CALL BLACS_GRIDEXIT( ICTXT_NEW )
      END IF
*
  120 CONTINUE
*
*     Restore saved input parameters
*
      ICTXT = ICTXT_SAVE
      NP = NP_SAVE
*
*     Output minimum worksize
*
      WORK( 1 ) = WORK_SIZE_MIN
*
*         Make INFO consistent across processors
*
      CALL IGAMX2D( ICTXT, 'A', ' ', 1, 1, INFO, 1, INFO, INFO, -1, 0,
     $              0 )
*
      IF( MYCOL.EQ.0 ) THEN
         CALL IGEBS2D( ICTXT, 'A', ' ', 1, 1, INFO, 1 )
      ELSE
         CALL IGEBR2D( ICTXT, 'A', ' ', 1, 1, INFO, 1, 0, 0 )
      END IF
*
*
      RETURN
*
*     End of PDPBTRF
*
      END
