      SUBROUTINE PCPBTRS( UPLO, N, BW, NRHS, A, JA, DESCA, B, IB, DESCB,
     $                    AF, LAF, WORK, LWORK, INFO )
*
*
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     August 7, 2001 
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            BW, IB, INFO, JA, LAF, LWORK, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX            A( * ), AF( * ), B( * ), WORK( * )
*     ..
*
*
*  Purpose
*  =======
*
*  PCPBTRS solves a system of linear equations
*
*            A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
*
*  where A(1:N, JA:JA+N-1) is the matrix used to produce the factors
*  stored in A(1:N,JA:JA+N-1) and AF by PCPBTRF.
*  A(1:N, JA:JA+N-1) is an N-by-N complex
*  banded symmetric positive definite distributed
*  matrix with bandwidth BW.
*  Depending on the value of UPLO, A stores either U or L in the equn
*  A(1:N, JA:JA+N-1) = U'*U or L*L' as computed by PCPBTRF.
*
*  Routine PCPBTRF MUST be called first.
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
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the distributed submatrix B(IB:IB+N-1, 1:NRHS).
*          NRHS >= 0.
*
*  A       (local input/local output) COMPLEX pointer into
*          local memory to an array with first dimension
*          LLD_A >=(bw+1) (stored in DESCA).
*          On entry, this array contains the local pieces of the
*          N-by-N symmetric banded distributed Cholesky factor L or
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
*          PCPBTRF and this is stored in AF. If a linear system
*          is to be solved using PCPBTRS after the factorization
*          routine, AF *must not be altered* after the factorization.
*
*  LAF     (local input) INTEGER
*          Size of user-input Auxiliary Fillin space AF. Must be >=
*          (NB+2*bw)*bw
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
*          (bw*NRHS)
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
*  Code Developer: Andrew J. Cleary, University of Tennessee.
*    Current address: Lawrence Livermore National Labs.
*  This version released: August, 2001.
*
*  =====================================================================
*
*     ..
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
      PARAMETER          (DESCMULT = 100, BIGNUM = DESCMULT * DESCMULT)
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            CSRC, FIRST_PROC, ICTXT, ICTXT_NEW, ICTXT_SAVE,
     $                   IDUM1, IDUM3, JA_NEW, LLDA, LLDB, MYCOL, MYROW,
     $                   NB, NP, NPCOL, NPROW, NP_SAVE, PART_OFFSET,
     $                   RETURN_CODE, STORE_M_B, STORE_N_A,
     $                   WORK_SIZE_MIN
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA_1XP( 7 ), DESCB_PX1( 7 ),
     $                   PARAM_CHECK( 16, 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESC_CONVERT, GLOBCHK,
     $                   PCPBTRSV, PXERBLA, RESHAPE
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
      DESCB_PX1( 1 ) = 502
*
      CALL DESC_CONVERT( DESCA, DESCA_1XP, RETURN_CODE )
*
      IF( RETURN_CODE .NE. 0) THEN
         INFO = -( 7*100 + 2 )
      ENDIF
*
      CALL DESC_CONVERT( DESCB, DESCB_PX1, RETURN_CODE )
*
      IF( RETURN_CODE .NE. 0) THEN
         INFO = -( 10*100 + 2 )
      ENDIF
*
*     Consistency checks for DESCA and DESCB.
*
*     Context must be the same
      IF( DESCA_1XP( 2 ) .NE. DESCB_PX1( 2 ) ) THEN
         INFO = -( 10*100 + 2 )
      ENDIF
*
*        These are alignment restrictions that may or may not be removed
*        in future releases. -Andy Cleary, April 14, 1996.
*
*     Block sizes must be the same
      IF( DESCA_1XP( 4 ) .NE. DESCB_PX1( 4 ) ) THEN
         INFO = -( 10*100 + 4 )
      ENDIF
*
*     Source processor must be the same
*
      IF( DESCA_1XP( 5 ) .NE. DESCB_PX1( 5 ) ) THEN
         INFO = -( 10*100 + 5 )
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
      IF( LSAME( UPLO, 'U' ) ) THEN
         IDUM1 = ICHAR( 'U' )
      ELSE IF ( LSAME( UPLO, 'L' ) ) THEN
         IDUM1 = ICHAR( 'L' )
      ELSE
         INFO = -1
      END IF
*
      IF( LWORK .LT. -1) THEN
         INFO = -14
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
         INFO = -( 7*100 + 6 )
      ENDIF
*
      IF(( BW .GT. N-1 ) .OR.
     $   ( BW .LT. 0 ) ) THEN
         INFO = -3
      ENDIF
*
      IF( LLDA .LT. (BW+1) ) THEN
         INFO = -( 7*100 + 6 )
      ENDIF
*
      IF( NB .LE. 0 ) THEN
         INFO = -( 7*100 + 4 )
      ENDIF
*
      IF( N+IB-1 .GT. STORE_M_B ) THEN
         INFO = -( 10*100 + 3 )
      ENDIF
*
      IF( LLDB .LT. NB ) THEN
         INFO = -( 10*100 + 6 )
      ENDIF
*
      IF( NRHS .LT. 0 ) THEN
         INFO = -3
      ENDIF
*
*     Current alignment restriction
*
      IF( JA .NE. IB) THEN
         INFO = -6
      ENDIF
*
*     Argument checking that is specific to Divide & Conquer routine
*
      IF( NPROW .NE. 1 ) THEN
         INFO = -( 7*100+2 )
      ENDIF
*
      IF( N .GT. NP*NB-MOD( JA-1, NB )) THEN
         INFO = -( 2 )
         CALL PXERBLA( ICTXT,
     $      'PCPBTRS, D&C alg.: only 1 block per proc',
     $      -INFO )
         RETURN
      ENDIF
*
      IF((JA+N-1.GT.NB) .AND. ( NB.LT.2*BW )) THEN
         INFO = -( 7*100+4 )
         CALL PXERBLA( ICTXT,
     $      'PCPBTRS, D&C alg.: NB too small',
     $      -INFO )
         RETURN
      ENDIF
*
*
      WORK_SIZE_MIN =
     $           (BW*NRHS)
*
      WORK( 1 ) = WORK_SIZE_MIN
*
      IF( LWORK .LT. WORK_SIZE_MIN ) THEN
         IF( LWORK .NE. -1 ) THEN
         INFO = -14
         CALL PXERBLA( ICTXT,
     $      'PCPBTRS: worksize error',
     $      -INFO )
         ENDIF
         RETURN
      ENDIF
*
*     Pack params and positions into arrays for global consistency check
*
      PARAM_CHECK( 16, 1 ) = DESCB(5)
      PARAM_CHECK( 15, 1 ) = DESCB(4)
      PARAM_CHECK( 14, 1 ) = DESCB(3)
      PARAM_CHECK( 13, 1 ) = DESCB(2)
      PARAM_CHECK( 12, 1 ) = DESCB(1)
      PARAM_CHECK( 11, 1 ) = IB
      PARAM_CHECK( 10, 1 ) = DESCA(5)
      PARAM_CHECK(  9, 1 ) = DESCA(4)
      PARAM_CHECK(  8, 1 ) = DESCA(3)
      PARAM_CHECK(  7, 1 ) = DESCA(1)
      PARAM_CHECK(  6, 1 ) = JA
      PARAM_CHECK(  5, 1 ) = NRHS
      PARAM_CHECK(  4, 1 ) = BW
      PARAM_CHECK(  3, 1 ) = N
      PARAM_CHECK(  2, 1 ) = IDUM3
      PARAM_CHECK(  1, 1 ) = IDUM1
*
      PARAM_CHECK( 16, 2 ) = 1005
      PARAM_CHECK( 15, 2 ) = 1004
      PARAM_CHECK( 14, 2 ) = 1003
      PARAM_CHECK( 13, 2 ) = 1002
      PARAM_CHECK( 12, 2 ) = 1001
      PARAM_CHECK( 11, 2 ) = 9
      PARAM_CHECK( 10, 2 ) = 705
      PARAM_CHECK(  9, 2 ) = 704
      PARAM_CHECK(  8, 2 ) = 703
      PARAM_CHECK(  7, 2 ) = 701
      PARAM_CHECK(  6, 2 ) = 6
      PARAM_CHECK(  5, 2 ) = 4
      PARAM_CHECK(  4, 2 ) = 3
      PARAM_CHECK(  3, 2 ) = 2
      PARAM_CHECK(  2, 2 ) = 14
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
      CALL GLOBCHK( ICTXT, 16, PARAM_CHECK, 16,
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
         CALL PXERBLA( ICTXT, 'PCPBTRS', -INFO )
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
      INFO = 0
*
*     Call frontsolve routine
*
      IF( LSAME( UPLO, 'L' ) ) THEN
*
         CALL PCPBTRSV( 'L', 'N', N, BW, NRHS, A( PART_OFFSET+1 ),
     $                  JA_NEW, DESCA_1XP, B, IB, DESCB_PX1, AF, LAF,
     $                  WORK, LWORK, INFO )
*
      ELSE
*
         CALL PCPBTRSV( 'U', 'C', N, BW, NRHS, A( PART_OFFSET+1 ),
     $                  JA_NEW, DESCA_1XP, B, IB, DESCB_PX1, AF, LAF,
     $                  WORK, LWORK, INFO )
*
      ENDIF
*
*     Call backsolve routine
*
      IF( LSAME( UPLO, 'L' ) ) THEN
*
         CALL PCPBTRSV( 'L', 'C', N, BW, NRHS, A( PART_OFFSET+1 ),
     $                  JA_NEW, DESCA_1XP, B, IB, DESCB_PX1, AF, LAF,
     $                  WORK, LWORK, INFO )
*
      ELSE
*
         CALL PCPBTRSV( 'U', 'N', N, BW, NRHS, A( PART_OFFSET+1 ),
     $                  JA_NEW, DESCA_1XP, B, IB, DESCB_PX1, AF, LAF,
     $                  WORK, LWORK, INFO )
*
      ENDIF
 1000 CONTINUE
*
*
*     Free BLACS space used to hold standard-form grid.
*
      IF( ICTXT_SAVE .NE. ICTXT_NEW ) THEN
         CALL BLACS_GRIDEXIT( ICTXT_NEW )
      ENDIF
*
 1234 CONTINUE
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
*
      RETURN
*
*     End of PCPBTRS
*
      END
