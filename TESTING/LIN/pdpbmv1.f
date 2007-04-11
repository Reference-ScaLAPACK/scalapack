      SUBROUTINE PDPBDCMV( LDBW, BW, UPLO, N, A, JA, DESCA, NRHS, B, IB,
     $                     DESCB, X, WORK, LWORK, INFO )
*
*
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            BW, IB, INFO, JA, LDBW, LWORK, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * ), WORK( * ), X( * )
*     ..
*
*
*  Purpose
*  =======
*
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
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0 )
      PARAMETER          ( ZERO = 0.0D+0 )
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
      INTEGER            CSRC, DL_N_M, DL_N_N, DL_P_M, DL_P_N,
     $                   FIRST_PROC, I, ICTXT, ICTXT_NEW, ICTXT_SAVE,
     $                   IDUM1, IDUM3, J, JA_NEW, LLDA, LLDB, MYCOL,
     $                   MYROW, MY_NUM_COLS, NB, NP, NPCOL, NPROW,
     $                   NP_SAVE, ODD_SIZE, OFST, PART_OFFSET,
     $                   PART_SIZE, STORE_M_B, STORE_N_A
      INTEGER NUMROC_SIZE
*     ..
*     .. Local Arrays ..
      INTEGER            PARAM_CHECK( 16, 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PXERBLA, RESHAPE
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
      ICTXT = DESCA( CTXT_ )
      CSRC = DESCA( CSRC_ )
      NB = DESCA( NB_ )
      LLDA = DESCA( LLD_ )
      STORE_N_A = DESCA( N_ )
      LLDB = DESCB( LLD_ )
      STORE_M_B = DESCB( M_ )
*
*
*     Pre-calculate bw^2
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
*     Argument checking that is specific to Divide & Conquer routine
*
      IF( NPROW .NE. 1 ) THEN
         INFO = -( 7*100+2 )
      ENDIF
*
      IF( N .GT. NP*NB-MOD( JA-1, NB )) THEN
         INFO = -( 2 )
         CALL PXERBLA( ICTXT,
     $      'PDPBDCMV, D&C alg.: only 1 block per proc',
     $      -INFO )
         RETURN
      ENDIF
*
      IF((JA+N-1.GT.NB) .AND. ( NB.LT.2*BW )) THEN
         INFO = -( 7*100+4 )
         CALL PXERBLA( ICTXT,
     $      'PDPBDCMV, D&C alg.: NB too small',
     $      -INFO )
         RETURN
      ENDIF
*
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
         CALL PXERBLA( ICTXT, 'PDPBDCMV', -INFO )
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
      ODD_SIZE = MY_NUM_COLS
      IF ( MYCOL .LT. NP-1 ) THEN
         ODD_SIZE = ODD_SIZE - BW
      ENDIF
*
*
*
*       Zero out solution to use to accumulate answer
*
        NUMROC_SIZE =
     $    NUMROC( N, PART_SIZE, MYCOL, 0, NPCOL)
*
        DO 2279 J=1,NRHS
          DO 4502 I=1,NUMROC_SIZE
            X( (J-1)*LLDB + I ) = ZERO
 4502     CONTINUE
 2279   CONTINUE
*
        DO 5642 I=1, (BW+2)*BW
          WORK( I ) = ZERO
 5642   CONTINUE
*
*     Begin main code
*
*
**************************************
*
      IF ( LSAME( UPLO, 'L' ) ) THEN
*
*       Sizes of the extra triangles communicated bewtween processors
*
        IF( MYCOL .GT. 0 ) THEN
*
          DL_P_M= MIN( BW,
     $          NUMROC( N, PART_SIZE, MYCOL, 0, NPCOL ) )
          DL_P_N= MIN( BW,
     $          NUMROC( N, PART_SIZE, MYCOL-1, 0, NPCOL ) )
        ENDIF
*
        IF( MYCOL .LT. NPCOL-1 ) THEN
*
          DL_N_M= MIN( BW,
     $          NUMROC( N, PART_SIZE, MYCOL+1, 0, NPCOL ) )
          DL_N_N= MIN( BW,
     $          NUMROC( N, PART_SIZE, MYCOL, 0, NPCOL ) )
        ENDIF
*
*
        IF( MYCOL .LT. NPCOL-1 ) THEN
*         ...must send triangle in upper half of matrix to right
*
*         Transpose triangle in preparation for sending
*
          CALL DLATCPY( 'U', BW, BW,
     $          A( LLDA*( NUMROC_SIZE-BW )+1+BW ),
     $          LLDA-1, WORK( 1 ), BW )
*
*         Send the triangle to neighboring processor to right
*
          CALL DTRSD2D(ICTXT, 'L', 'N',
     $                  BW, BW,
     $                  WORK( 1 ),
     $                  BW, MYROW, MYCOL+1 )
*
        ENDIF
*
*       Use main partition in each processor to multiply locally
*
        CALL DSBMV( 'L', NUMROC_SIZE, BW, ONE, A( OFST+1 ), LLDA,
     $              B(PART_OFFSET+1), 1, ZERO, X( PART_OFFSET+1 ), 1 )
*
*
*
        IF ( MYCOL .LT. NPCOL-1 ) THEN
*
*         Do the multiplication of the triangle in the lower half
*
          CALL DCOPY( DL_N_N,
     $                  B( NUMROC_SIZE-DL_N_N+1 ),
     $                  1, WORK( BW*BW+1+BW-DL_N_N ), 1 )
*
         CALL DTRMV( 'U', 'N', 'N', BW,
     $            A( LLDA*( NUMROC_SIZE-BW )+1+BW ), LLDA-1,
     $            WORK( BW*BW+1 ), 1)
*
*        Zero out extraneous elements caused by TRMV if any
*
         IF( DL_N_M .GT. DL_N_N ) THEN
        DO 10  I = DL_N_M-DL_N_N, DL_N_M
                WORK( BW*BW+I ) = 0
   10   CONTINUE
         ENDIF
*
*         Send the result to the neighbor
*
          CALL DGESD2D( ICTXT, BW, 1,
     $       WORK( BW*BW+1 ), BW, MYROW, MYCOL+1 )
*
        ENDIF
*
        IF ( MYCOL .GT. 0 ) THEN
*
        DO 20  I=1, BW*( BW+2 )
          WORK( I ) = ZERO
   20   CONTINUE
*
*         Do the multiplication of the triangle in the upper half
*
*         Copy vector to workspace
*
          CALL DCOPY( DL_P_M, B( 1 ), 1,
     $                  WORK( BW*BW+1 ), 1)
*
*         Receive the triangle prior to multiplying by it.
*
          CALL DTRRV2D(ICTXT, 'L', 'N',
     $                  BW, BW,
     $                  WORK( 1 ), BW, MYROW, MYCOL-1 )
*
          CALL DTRMV(
     $     'L',
     $     'N',
     $     'N', BW,
     $     WORK( 1 ), BW,
     $     WORK( BW*BW+1 ), 1 )
*
*         Zero out extraneous results from TRMV if any
*
          IF( DL_P_M .GT. DL_P_N ) THEN
        DO 30  I=1, DL_P_M-DL_P_N
              WORK( BW*BW+I ) = 0
   30   CONTINUE
          ENDIF
*
*         Send result back
*
          CALL DGESD2D( ICTXT, BW, 1, WORK(BW*BW+1 ),
     $                   BW, MYROW, MYCOL-1 )
*
*         Receive vector result from neighboring processor
*
          CALL DGERV2D( ICTXT, BW, 1, WORK( BW*BW+1 ),
     $                    BW, MYROW, MYCOL-1 )
*
*         Do addition of received vector
*
          CALL DAXPY( BW, ONE,
     $                  WORK( BW*BW+1 ), 1,
     $                  X( 1 ), 1 )
*
        ENDIF
*
*
*
         IF( MYCOL .LT. NPCOL-1 ) THEN
*
*          Receive returned result
*
           CALL DGERV2D( ICTXT, BW, 1, WORK( BW*BW+1 ),
     $                    BW, MYROW, MYCOL+1 )
*
*          Do addition of received vector
*
           CALL DAXPY( BW, ONE,
     $                  WORK( BW*BW+1 ), 1,
     $                  X( NUMROC_SIZE-BW+1 ), 1)
*
         ENDIF
*
*
      ENDIF
*
*     End of LSAME if
*
**************************************
*
      IF ( LSAME( UPLO, 'U' ) ) THEN
*
*       Sizes of the extra triangles communicated bewtween processors
*
        IF( MYCOL .GT. 0 ) THEN
*
          DL_P_M= MIN( BW,
     $          NUMROC( N, PART_SIZE, MYCOL, 0, NPCOL ) )
          DL_P_N= MIN( BW,
     $          NUMROC( N, PART_SIZE, MYCOL-1, 0, NPCOL ) )
        ENDIF
*
        IF( MYCOL .LT. NPCOL-1 ) THEN
*
          DL_N_M= MIN( BW,
     $          NUMROC( N, PART_SIZE, MYCOL+1, 0, NPCOL ) )
          DL_N_N= MIN( BW,
     $          NUMROC( N, PART_SIZE, MYCOL, 0, NPCOL ) )
        ENDIF
*
*
        IF( MYCOL .GT. 0 ) THEN
*         ...must send triangle in lower half of matrix to left
*
*         Transpose triangle in preparation for sending
*
          CALL DLATCPY( 'L', BW, BW, A( OFST+1 ),
     $          LLDA-1, WORK( 1 ), BW )
*
*         Send the triangle to neighboring processor to left
*
          CALL DTRSD2D(ICTXT, 'U', 'N',
     $                  BW, BW,
     $                  WORK( 1 ),
     $                  BW, MYROW, MYCOL-1 )
*
        ENDIF
*
*       Use main partition in each processor to multiply locally
*
        CALL DSBMV( 'U', NUMROC_SIZE, BW, ONE, A( OFST+1 ), LLDA,
     $              B(PART_OFFSET+1), 1, ZERO, X( PART_OFFSET+1 ), 1 )
*
*
*
        IF ( MYCOL .LT. NPCOL-1 ) THEN
*
*         Do the multiplication of the triangle in the lower half
*
          CALL DCOPY( DL_N_N,
     $                  B( NUMROC_SIZE-DL_N_N+1 ),
     $                  1, WORK( BW*BW+1+BW-DL_N_N ), 1 )
*
*         Receive the triangle prior to multiplying by it.
*
          CALL DTRRV2D(ICTXT, 'U', 'N',
     $                  BW, BW,
     $                  WORK( 1 ), BW, MYROW, MYCOL+1 )
*
         CALL DTRMV( 'U', 'N', 'N', BW,
     $            WORK( 1 ), BW,
     $            WORK( BW*BW+1 ), 1)
*
*        Zero out extraneous elements caused by TRMV if any
*
         IF( DL_N_M .GT. DL_N_N ) THEN
        DO 40  I = DL_N_M-DL_N_N, DL_N_M
                WORK( BW*BW+I ) = 0
   40   CONTINUE
         ENDIF
*
*         Send the result to the neighbor
*
          CALL DGESD2D( ICTXT, BW, 1,
     $       WORK( BW*BW+1 ), BW, MYROW, MYCOL+1 )
*
        ENDIF
*
        IF ( MYCOL .GT. 0 ) THEN
*
        DO 50  I=1, BW*( BW+2 )
          WORK( I ) = ZERO
   50   CONTINUE
*
*         Do the multiplication of the triangle in the upper half
*
*         Copy vector to workspace
*
          CALL DCOPY( DL_P_M, B( 1 ), 1,
     $                  WORK( BW*BW+1 ), 1)
*
          CALL DTRMV(
     $     'L',
     $     'N',
     $     'N', BW,
     $     A( 1 ), LLDA-1,
     $     WORK( BW*BW+1 ), 1 )
*
*         Zero out extraneous results from TRMV if any
*
          IF( DL_P_M .GT. DL_P_N ) THEN
        DO 60  I=1, DL_P_M-DL_P_N
              WORK( BW*BW+I ) = 0
   60   CONTINUE
          ENDIF
*
*         Send result back
*
          CALL DGESD2D( ICTXT, BW, 1, WORK(BW*BW+1 ),
     $                   BW, MYROW, MYCOL-1 )
*
*         Receive vector result from neighboring processor
*
          CALL DGERV2D( ICTXT, BW, 1, WORK( BW*BW+1 ),
     $                    BW, MYROW, MYCOL-1 )
*
*         Do addition of received vector
*
          CALL DAXPY( BW, ONE,
     $                  WORK( BW*BW+1 ), 1,
     $                  X( 1 ), 1 )
*
        ENDIF
*
*
*
         IF( MYCOL .LT. NPCOL-1 ) THEN
*
*          Receive returned result
*
           CALL DGERV2D( ICTXT, BW, 1, WORK( BW*BW+1 ),
     $                    BW, MYROW, MYCOL+1 )
*
*          Do addition of received vector
*
           CALL DAXPY( BW, ONE,
     $                  WORK( BW*BW+1 ), 1,
     $                  X( NUMROC_SIZE-BW+1 ), 1)
*
         ENDIF
*
*
      ENDIF
*
*     End of LSAME if
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
*
      RETURN
*
*     End of PDBsBMV1
*
      END
