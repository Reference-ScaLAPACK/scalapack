      SUBROUTINE PSPTSV( N, NRHS, D, E, JA, DESCA, B, IB, DESCB, WORK,
     $                   LWORK, INFO )
*
*
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IB, INFO, JA, LWORK, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      REAL               B( * ), D( * ), E( * ), WORK( * )
*     ..
*
*
*  Purpose
*  =======
*
*  PSPTSV solves a system of linear equations
*
*                      A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
*
*  where A(1:N, JA:JA+N-1) is an N-by-N real
*  tridiagonal symmetric positive definite distributed
*  matrix.
*
*  Cholesky factorization is used to factor a reordering of
*  the matrix into L L'.
*
*  See PSPTTRF and PSPTTRS for details.
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
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the distributed submatrix B(IB:IB+N-1, 1:NRHS).
*          NRHS >= 0.
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
*  B       (local input/local output) REAL pointer into
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
*          (12*NPCOL + 3*NB)
*          +max((10+2*min(100,NRHS))*NPCOL+4*NRHS, 8*NPCOL)
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
*  Ignored         DESCA( 6 ) Ignored for tridiagonal matrices.
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
      INTEGER            ICTXT, MYCOL, MYROW, NB, NPCOL, NPROW,
     $                   WS_FACTOR
*     ..
*     .. External Subroutines ..
      EXTERNAL           PSPTTRF, PSPTTRS, PXERBLA
*     ..
*     .. Executable Statements ..
*
*     Note: to avoid duplication, most error checking is not performed
*           in this routine and is left to routines
*           PSPTTRF and PSPTTRS.
*
*     Begin main code
*
      INFO = 0
*
*     Get block size to calculate workspace requirements
*
      IF( DESCA( DTYPE_ ) .EQ. BLOCK_CYCLIC_2D ) THEN
         NB = DESCA( NB_ )
         ICTXT = DESCA( CTXT_ )
      ELSEIF( DESCA( DTYPE_ ) .EQ. 501 ) THEN
         NB = DESCA( 4 )
         ICTXT = DESCA( 2 )
      ELSEIF( DESCA( DTYPE_ ) .EQ. 502 ) THEN
         NB = DESCA( 4 )
         ICTXT = DESCA( 2 )
      ELSE
         INFO = -( 5*100 + DTYPE_ )
         CALL PXERBLA( ICTXT,
     $      'PSPTSV',
     $      -INFO )
         RETURN
      ENDIF
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*
*     Size needed for AF in factorization
*
      WS_FACTOR = (12*NPCOL + 3*NB)
*
*     Factor the matrix
*
      CALL PSPTTRF( N, D, E, JA, DESCA, WORK, MIN( LWORK, WS_FACTOR ),
     $              WORK( 1+WS_FACTOR ), LWORK-WS_FACTOR, INFO )
*
*     Check info for error conditions
*
      IF( INFO.NE.0 ) THEN
         IF( INFO .LT. 0 ) THEN
            CALL PXERBLA( ICTXT, 'PSPTSV', -INFO )
         ENDIF
         RETURN
      END IF
*
*     Solve the system using the factorization
*
      CALL PSPTTRS( N, NRHS, D, E, JA, DESCA, B, IB, DESCB, WORK,
     $              MIN( LWORK, WS_FACTOR ), WORK( 1+WS_FACTOR),
     $              LWORK-WS_FACTOR, INFO )
*
*     Check info for error conditions
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PSPTSV', -INFO )
         RETURN
      END IF
*
      RETURN
*
*     End of PSPTSV
*
      END
