      SUBROUTINE PDGESV( N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB,
     $                   DESCB, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     Jan 30, 2006
*
*     .. Scalar Arguments ..
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      DOUBLE PRECISION   A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PDGESV computes the solution to a real system of linear equations
*
*                        sub( A ) * X = sub( B ),
*
*  where sub( A ) = A(IA:IA+N-1,JA:JA+N-1) is an N-by-N distributed
*  matrix and X and sub( B ) = B(IB:IB+N-1,JB:JB+NRHS-1) are N-by-NRHS
*  distributed matrices.
*
*  The LU decomposition with partial pivoting and row interchanges is
*  used to factor sub( A ) as sub( A ) = P * L * U, where P is a permu-
*  tation matrix, L is unit lower triangular, and U is upper triangular.
*  L and U are stored in sub( A ). The factored form of sub( A ) is then
*  used to solve the system of equations sub( A ) * X = sub( B ).
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
*  This routine requires square block decomposition ( MB_A = NB_A ).
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the distributed submatrix sub( B ). NRHS >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, the local pieces of the N-by-N distributed matrix
*          sub( A ) to be factored. On exit, this array contains the
*          local pieces of the factors L and U from the factorization
*          sub( A ) = P*L*U; the unit diagonal elements of L are not
*          stored.
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
*  IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
*          This array contains the pivoting information.
*          IPIV(i) -> The global row local row i was swapped with.
*          This array is tied to the distributed matrix A.
*
*  B       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension
*          (LLD_B,LOCc(JB+NRHS-1)).  On entry, the right hand side
*          distributed matrix sub( B ). On exit, if INFO = 0, sub( B )
*          is overwritten by the solution distributed matrix X.
*
*  IB      (global input) INTEGER
*          The row index in the global array B indicating the first
*          row of sub( B ).
*
*  JB      (global input) INTEGER
*          The column index in the global array B indicating the
*          first column of sub( B ).
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = K, U(IA+K-1,JA+K-1) is exactly zero.
*                The factorization has been completed, but the factor U
*                is exactly singular, so the solution could not be
*                computed.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            IAROW, IBROW, ICOFFA, ICTXT, IROFFA, IROFFB,
     $                   MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK2MAT, PDGETRF,
     $                   PDGETRS, PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P
      EXTERNAL           INDXG2P
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(600+CTXT_)
      ELSE
         CALL CHK1MAT( N, 1, N, 1, IA, JA, DESCA, 6, INFO )
         CALL CHK1MAT( N, 1, NRHS, 2, IB, JB, DESCB, 11, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $                       NPROW )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            IF( IROFFA.NE.0 ) THEN
               INFO = -4
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(600+NB_)
            ELSE IF( IBROW.NE.IAROW .OR. ICOFFA.NE.IROFFB ) THEN
               INFO = -9
            ELSE IF( DESCB( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(1100+NB_)
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -(1100+CTXT_)
            END IF
         END IF
         CALL PCHK2MAT( N, 1, N, 1, IA, JA, DESCA, 6, N, 1, NRHS, 2,
     $                  IB, JB, DESCB, 11, 0, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDGESV', -INFO )
         RETURN
      END IF
*
*     Compute the LU factorization of sub( A ).
*
      CALL PDGETRF( N, N, A, IA, JA, DESCA, IPIV, INFO )
*
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system sub( A ) * X = sub( B ), overwriting sub( B )
*        with X.
*
         CALL PDGETRS( 'No transpose', N, NRHS, A, IA, JA, DESCA, IPIV,
     $                 B, IB, JB, DESCB, INFO )
*
      END IF
*
      RETURN
*
*     End of PDGESV
*
      END
