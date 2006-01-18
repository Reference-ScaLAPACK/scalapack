      SUBROUTINE PSPTTRF( N, D, E, JA, DESCA, AF, LAF, WORK, LWORK,
     $                    INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 3, 2000
*
*     .. Scalar Arguments ..
      INTEGER            INFO, JA, LAF, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      REAL               AF( * ), D( * ), E( * ), WORK( * )
*     ..
*
*
*  Purpose
*  =======
*
*  PSPTTRF computes a Cholesky factorization
*  of an N-by-N real tridiagonal
*  symmetric positive definite distributed matrix
*  A(1:N, JA:JA+N-1).
*  Reordering is used to increase parallelism in the factorization.
*  This reordering results in factors that are DIFFERENT from those
*  produced by equivalent sequential codes. These factors cannot
*  be used directly by users; however, they can be used in
*  subsequent calls to PSPTTRS to solve linear systems.
*
*  The factorization has the form
*
*          P A(1:N, JA:JA+N-1) P^T = U' D U  or
*
*          P A(1:N, JA:JA+N-1) P^T = L D L',
*
*  where U is a tridiagonal upper triangular matrix and L is tridiagonal
*  lower triangular, and P is a permutation matrix.
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
*  D       (local input/local output) REAL pointer to local
*          part of global vector storing the main diagonal of the
*          matrix.
*          On exit, this array contains information containing the
*            factors of the matrix.
*          Must be of size >= DESCA( NB_ ).
*
*  E       (local input/local output) REAL pointer to local
*          part of global vector storing the upper diagonal of the
*          matrix. Globally, DU(n) is not referenced, and DU must be
*          aligned with D.
*          On exit, this array contains information containing the
*            factors of the matrix.
*          Must be of size >= DESCA( NB_ ).
*
*  JA      (global input) INTEGER
*          The index in the global array A that points to the start of
*          the matrix to be operated on (which may be either all of A
*          or a submatrix of A).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN.
*          if 1D type (DTYPE_A=501 or 502), DLEN >= 7;
*          if 2D type (DTYPE_A=1), DLEN >= 9.
*          The array descriptor for the distributed matrix A.
*          Contains information of mapping of A to memory. Please
*          see NOTES below for full description and options.
*
*  AF      (local output) REAL array, dimension LAF.
*          Auxiliary Fillin Space.
*          Fillin is created during the factorization routine
*          PSPTTRF and this is stored in AF. If a linear system
*          is to be solved using PSPTTRS after the factorization
*          routine, AF *must not be altered* after the factorization.
*
*  LAF     (local input) INTEGER
*          Size of user-input Auxiliary Fillin space AF. Must be >=
*          (NB+2)
*          If LAF is not large enough, an error code will be returned
*          and the minimum acceptable size will be returned in AF( 1 )
*
*  WORK    (local workspace/local output)
*          REAL temporary workspace. This space may
*          be overwritten in between calls to routines. WORK must be
*          the size given in LWORK.
*          On exit, WORK( 1 ) contains the minimal LWORK.
*
*  LWORK   (local input or global input) INTEGER
*          Size of user-input workspace WORK.
*          If LWORK is too small, the minimal acceptable size will be
*          returned in WORK(1) and an error code is returned. LWORK>=
*          8*NPCOL
*
*  INFO    (local output) INTEGER
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
*      NB >= 2
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
*    For tridiagonal matrices, it is obvious: N/P >> bw(=1), and so D&C
*    algorithms are the appropriate choice.
*
*  Algorithm description: Divide and Conquer
*
*    The Divide and Conqer algorithm assumes the matrix is narrowly
*      banded compared with the number of equations. In this situation,
*      it is best to distribute the input matrix A one-dimensionally,
*      with columns atomic and rows divided amongst the processes.
*      The basic algorithm divides the tridiagonal matrix up into
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
*         A small ((P-1)) system is formed representing
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
*  Note: tridiagonal codes can use either the old two dimensional
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
*    However, for tridiagonal matrices, since the objects being
*    distributed are the individual vectors storing the diagonals, we
*    have adopted the convention that both the P-by-1 descriptor and
*    the 1-by-P descriptor are allowed and are equivalent for
*    tridiagonal matrices. Thus, for tridiagonal matrices,
*    DTYPE_A = 501 or 502 can be used interchangeably
*    without any other change.
*  We require that the distributed vectors storing the diagonals of a
*    tridiagonal matrix be aligned with each other. Because of this, a
*    single descriptor, DESCA, serves to describe the distribution of
*    of all diagonals simultaneously.
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
*  A               OK          OK          OK        NO
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
*  Ignored         DESCA( 6 ) Ignored for tridiagonal matrices.
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
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
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
     $                   ICTXT_NEW, ICTXT_SAVE, IDUM3, INT_TEMP, JA_NEW,
     $                   LAF_MIN, LEVEL_DIST, LLDA, MYCOL, MYROW,
     $                   MY_NUM_COLS, NB, NP, NPCOL, NPROW, NP_SAVE,
     $                   ODD_SIZE, PART_OFFSET, PART_SIZE, RETURN_CODE,
     $                   STORE_N_A, TEMP, WORK_SIZE_MIN
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA_1XP( 7 ), PARAM_CHECK( 7, 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDEXIT, BLACS_GRIDINFO, DESC_CONVERT,
     $                   GLOBCHK, IGAMX2D, IGEBR2D, IGEBS2D, PXERBLA,
     $                   RESHAPE, SGERV2D, SGESD2D, SPTTRF, SPTTRSV,
     $                   STRRV2D, STRSD2D
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD, REAL
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
      TEMP = DESCA( DTYPE_ )
      IF( TEMP.EQ.502 ) THEN
*        Temporarily set the descriptor type to 1xP type
         DESCA( DTYPE_ ) = 501
      END IF
*
      CALL DESC_CONVERT( DESCA, DESCA_1XP, RETURN_CODE )
*
      DESCA( DTYPE_ ) = TEMP
*
      IF( RETURN_CODE.NE.0 ) THEN
         INFO = -( 5*100+2 )
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
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NP = NPROW*NPCOL
*
*
*
      IF( LWORK.LT.-1 ) THEN
         INFO = -9
      ELSE IF( LWORK.EQ.-1 ) THEN
         IDUM3 = -1
      ELSE
         IDUM3 = 1
      END IF
*
      IF( N.LT.0 ) THEN
         INFO = -1
      END IF
*
      IF( N+JA-1.GT.STORE_N_A ) THEN
         INFO = -( 5*100+6 )
      END IF
*
*     Argument checking that is specific to Divide & Conquer routine
*
      IF( NPROW.NE.1 ) THEN
         INFO = -( 5*100+2 )
      END IF
*
      IF( N.GT.NP*NB-MOD( JA-1, NB ) ) THEN
         INFO = -( 1 )
         CALL PXERBLA( ICTXT, 'PSPTTRF, D&C alg.: only 1 block per proc'
     $                 , -INFO )
         RETURN
      END IF
*
      IF( ( JA+N-1.GT.NB ) .AND. ( NB.LT.2*INT_ONE ) ) THEN
         INFO = -( 5*100+4 )
         CALL PXERBLA( ICTXT, 'PSPTTRF, D&C alg.: NB too small', -INFO )
         RETURN
      END IF
*
*
*     Check auxiliary storage size
*
      LAF_MIN = ( 12*NPCOL+3*NB )
*
      IF( LAF.LT.LAF_MIN ) THEN
         INFO = -7
*        put minimum value of laf into AF( 1 )
         AF( 1 ) = LAF_MIN
         CALL PXERBLA( ICTXT, 'PSPTTRF: auxiliary storage error ',
     $                 -INFO )
         RETURN
      END IF
*
*     Check worksize
*
      WORK_SIZE_MIN = 8*NPCOL
*
      WORK( 1 ) = WORK_SIZE_MIN
*
      IF( LWORK.LT.WORK_SIZE_MIN ) THEN
         IF( LWORK.NE.-1 ) THEN
            INFO = -9
            CALL PXERBLA( ICTXT, 'PSPTTRF: worksize error ', -INFO )
         END IF
         RETURN
      END IF
*
*     Pack params and positions into arrays for global consistency check
*
      PARAM_CHECK( 7, 1 ) = DESCA( 5 )
      PARAM_CHECK( 6, 1 ) = DESCA( 4 )
      PARAM_CHECK( 5, 1 ) = DESCA( 3 )
      PARAM_CHECK( 4, 1 ) = DESCA( 1 )
      PARAM_CHECK( 3, 1 ) = JA
      PARAM_CHECK( 2, 1 ) = N
      PARAM_CHECK( 1, 1 ) = IDUM3
*
      PARAM_CHECK( 7, 2 ) = 505
      PARAM_CHECK( 6, 2 ) = 504
      PARAM_CHECK( 5, 2 ) = 503
      PARAM_CHECK( 4, 2 ) = 501
      PARAM_CHECK( 3, 2 ) = 4
      PARAM_CHECK( 2, 2 ) = 1
      PARAM_CHECK( 1, 2 ) = 9
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
      CALL GLOBCHK( ICTXT, 7, PARAM_CHECK, 7, PARAM_CHECK( 1, 3 ),
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
         CALL PXERBLA( ICTXT, 'PSPTTRF', -INFO )
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
         GO TO 90
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
*     Size of main (or odd) partition in each processor
*
      ODD_SIZE = MY_NUM_COLS
      IF( MYCOL.LT.NP-1 ) THEN
         ODD_SIZE = ODD_SIZE - INT_ONE
      END IF
*
*
*       Zero out space for fillin
*
      DO 10 I = 1, LAF_MIN
         AF( I ) = ZERO
   10 CONTINUE
*
*     Begin main code
*
*
********************************************************************
*       PHASE 1: Local computation phase.
********************************************************************
*
*
      IF( MYCOL.LT.NP-1 ) THEN
*         Transfer last triangle D_i of local matrix to next processor
*         which needs it to calculate fillin due to factorization of
*         its main (odd) block A_i.
*         Overlap the send with the factorization of A_i.
*
         CALL STRSD2D( ICTXT, 'U', 'N', 1, 1,
     $                 E( PART_OFFSET+ODD_SIZE+1 ), LLDA-1, 0, MYCOL+1 )
*
      END IF
*
*
*       Factor main partition A_i = L_i {L_i}^T in each processor
*       Or A_i = {U_i}^T {U_i} if E is the upper superdiagonal
*
      CALL SPTTRF( ODD_SIZE, D( PART_OFFSET+1 ), E( PART_OFFSET+1 ),
     $             INFO )
*
      IF( INFO.NE.0 ) THEN
         INFO = MYCOL + 1
         GO TO 20
      END IF
*
*
      IF( MYCOL.LT.NP-1 ) THEN
*         Apply factorization to odd-even connection block B_i
*
*
*         Perform the triangular system solve {L_i}{{B'}_i}^T = {B_i}^T
*           by dividing B_i by diagonal element
*
         E( PART_OFFSET+ODD_SIZE ) = E( PART_OFFSET+ODD_SIZE ) /
     $                               D( PART_OFFSET+ODD_SIZE )
*
*
*
*         Compute contribution to diagonal block(s) of reduced system.
*          {C'}_i = {C_i}-{{B'}_i}{{B'}_i}^T
*
         D( PART_OFFSET+ODD_SIZE+1 ) = D( PART_OFFSET+ODD_SIZE+1 ) -
     $                                 D( PART_OFFSET+ODD_SIZE )*
     $                                 ( E( PART_OFFSET+ODD_SIZE )*
     $                                 ( E( PART_OFFSET+ODD_SIZE ) ) )
*
      END IF
*       End of "if ( MYCOL .lt. NP-1 )..." loop
*
*
   20 CONTINUE
*       If the processor could not locally factor, it jumps here.
*
      IF( MYCOL.NE.0 ) THEN
*
*         Receive previously transmitted matrix section, which forms
*         the right-hand-side for the triangular solve that calculates
*         the "spike" fillin.
*
*
         CALL STRRV2D( ICTXT, 'U', 'N', 1, 1, AF( 1 ), ODD_SIZE, 0,
     $                 MYCOL-1 )
*
         IF( INFO.EQ.0 ) THEN
*
*         Calculate the "spike" fillin, ${L_i} {{G}_i}^T = {D_i}$ .
*
            CALL SPTTRSV( 'N', ODD_SIZE, INT_ONE, D( PART_OFFSET+1 ),
     $                    E( PART_OFFSET+1 ), AF( 1 ), ODD_SIZE, INFO )
*
*         Divide by D
*
            DO 30 I = 1, ODD_SIZE
               AF( I ) = AF( I ) / D( PART_OFFSET+I )
   30       CONTINUE
*
*
*         Calculate the update block for previous proc, E_i = G_i{G_i}^T
*
*
*         Since there is no element-by-element vector multiplication in
*           the BLAS, this loop must be hardwired in without a BLAS call
*
            INT_TEMP = ODD_SIZE*INT_ONE + 2 + 1
            AF( INT_TEMP ) = 0
*
            DO 40 I = 1, ODD_SIZE
               AF( INT_TEMP ) = AF( INT_TEMP ) -
     $                          D( PART_OFFSET+I )*( AF( I )*
     $                          ( AF( I ) ) )
   40       CONTINUE
*
*
*         Initiate send of E_i to previous processor to overlap
*           with next computation.
*
            CALL SGESD2D( ICTXT, INT_ONE, INT_ONE, AF( ODD_SIZE+3 ),
     $                    INT_ONE, 0, MYCOL-1 )
*
*
            IF( MYCOL.LT.NP-1 ) THEN
*
*           Calculate off-diagonal block(s) of reduced system.
*           Note: for ease of use in solution of reduced system, store
*           L's off-diagonal block in transpose form.
*           {F_i}^T =  {H_i}{{B'}_i}^T
*
               AF( ODD_SIZE+1 ) = -D( PART_OFFSET+ODD_SIZE )*
     $                            ( E( PART_OFFSET+ODD_SIZE )*
     $                            AF( ODD_SIZE ) )
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
      CALL IGAMX2D( ICTXT, 'A', ' ', 1, 1, INFO, 1, INFO, INFO, -1, 0,
     $              0 )
*
      IF( MYCOL.EQ.0 ) THEN
         CALL IGEBS2D( ICTXT, 'A', ' ', 1, 1, INFO, 1 )
      ELSE
         CALL IGEBR2D( ICTXT, 'A', ' ', 1, 1, INFO, 1, 0, 0 )
      END IF
*
      IF( INFO.NE.0 ) THEN
         GO TO 80
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
         GO TO 70
      END IF
*
*       Initiate send of off-diag block(s) to overlap with next part.
*       Off-diagonal block needed on neighboring processor to start
*       algorithm.
*
      IF( ( MOD( MYCOL+1, 2 ).EQ.0 ) .AND. ( MYCOL.GT.0 ) ) THEN
*
         CALL SGESD2D( ICTXT, INT_ONE, INT_ONE, AF( ODD_SIZE+1 ),
     $                 INT_ONE, 0, MYCOL-1 )
*
      END IF
*
*       Copy last diagonal block into AF storage for subsequent
*         operations.
*
      AF( ODD_SIZE+2 ) = REAL( D( PART_OFFSET+ODD_SIZE+1 ) )
*
*       Receive cont. to diagonal block that is stored on this proc.
*
      IF( MYCOL.LT.NPCOL-1 ) THEN
*
         CALL SGERV2D( ICTXT, INT_ONE, INT_ONE, AF( ODD_SIZE+2+1 ),
     $                 INT_ONE, 0, MYCOL+1 )
*
*          Add contribution to diagonal block
*
         AF( ODD_SIZE+2 ) = AF( ODD_SIZE+2 ) + AF( ODD_SIZE+3 )
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
   50 CONTINUE
      IF( MOD( ( MYCOL+1 ) / LEVEL_DIST, 2 ).NE.0 )
     $   GO TO 60
*
*         Receive and add contribution to diagonal block from the left
*
      IF( MYCOL-LEVEL_DIST.GE.0 ) THEN
         CALL SGERV2D( ICTXT, INT_ONE, INT_ONE, WORK( 1 ), INT_ONE, 0,
     $                 MYCOL-LEVEL_DIST )
*
         AF( ODD_SIZE+2 ) = AF( ODD_SIZE+2 ) + WORK( 1 )
*
      END IF
*
*         Receive and add contribution to diagonal block from the right
*
      IF( MYCOL+LEVEL_DIST.LT.NPCOL-1 ) THEN
         CALL SGERV2D( ICTXT, INT_ONE, INT_ONE, WORK( 1 ), INT_ONE, 0,
     $                 MYCOL+LEVEL_DIST )
*
         AF( ODD_SIZE+2 ) = AF( ODD_SIZE+2 ) + WORK( 1 )
*
      END IF
*
      LEVEL_DIST = LEVEL_DIST*2
*
      GO TO 50
   60 CONTINUE
*       [End of GOTO Loop]
*
*
*       *********************************
*       Calculate and use this proc's blocks to modify other procs'...
      IF( AF( ODD_SIZE+2 ).EQ.ZERO ) THEN
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
         AF( ODD_SIZE+3 ) = AF( ODD_SIZE+1 )
*
      ELSE
         COMM_PROC = MYCOL + LEVEL_DIST / 2
      END IF
*
      IF( MYCOL / LEVEL_DIST.LE.( NPCOL-1 ) / LEVEL_DIST-2 ) THEN
*
         CALL SGERV2D( ICTXT, INT_ONE, INT_ONE, AF( ODD_SIZE+1 ),
     $                 INT_ONE, 0, COMM_PROC )
*
         IF( INFO.EQ.0 ) THEN
*
*
*         Modify upper off_diagonal block with diagonal block
*
*
            AF( ODD_SIZE+1 ) = AF( ODD_SIZE+1 ) / AF( ODD_SIZE+2 )
*
         END IF
*         End of "if ( info.eq.0 ) then"
*
*         Calculate contribution from this block to next diagonal block
*
         WORK( 1 ) = -ONE*AF( ODD_SIZE+1 )*AF( ODD_SIZE+2 )*
     $               ( AF( ODD_SIZE+1 ) )
*
*         Send contribution to diagonal block's owning processor.
*
         CALL SGESD2D( ICTXT, INT_ONE, INT_ONE, WORK( 1 ), INT_ONE, 0,
     $                 MYCOL+LEVEL_DIST )
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
     $    ( MYCOL / LEVEL_DIST.LE.( NPCOL-1 ) / LEVEL_DIST-1 ) ) THEN
*
         IF( LEVEL_DIST.GT.1 ) THEN
*
*           Receive offdiagonal block(s) from proc level_dist/2 to the
*           left
*
            CALL SGERV2D( ICTXT, INT_ONE, INT_ONE, AF( ODD_SIZE+2+1 ),
     $                    INT_ONE, 0, MYCOL-LEVEL_DIST / 2 )
*
         END IF
*
*
         IF( INFO.EQ.0 ) THEN
*
*         Use diagonal block(s) to modify this offdiagonal block
*
            AF( ODD_SIZE+3 ) = ( AF( ODD_SIZE+3 ) ) / AF( ODD_SIZE+2 )
*
         END IF
*         End of "if( info.eq.0 ) then"
*
*         Use offdiag block(s) to calculate modification to diag block
*           of processor to the left
*
         WORK( 1 ) = -ONE*AF( ODD_SIZE+3 )*AF( ODD_SIZE+2 )*
     $               ( AF( ODD_SIZE+3 ) )
*
*         Send contribution to diagonal block's owning processor.
*
         CALL SGESD2D( ICTXT, INT_ONE, INT_ONE, WORK( 1 ), INT_ONE, 0,
     $                 MYCOL-LEVEL_DIST )
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
            WORK( 1 ) = -ONE*AF( ODD_SIZE+3 )*AF( ODD_SIZE+2 )*
     $                  AF( ODD_SIZE+1 )
*
*           Send contribution to offdiagonal block's owning processor.
*
            CALL SGESD2D( ICTXT, INT_ONE, INT_ONE, WORK( 1 ), INT_ONE,
     $                    0, COMM_PROC )
*
         END IF
*
      END IF
*       End of "if( mycol/level_dist.le. (npcol-1)/level_dist -1 )..."
*
   70 CONTINUE
*
*
   80 CONTINUE
*
*
*     Free BLACS space used to hold standard-form grid.
*
      IF( ICTXT_SAVE.NE.ICTXT_NEW ) THEN
         CALL BLACS_GRIDEXIT( ICTXT_NEW )
      END IF
*
   90 CONTINUE
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
*     End of PSPTTRF
*
      END
