      SUBROUTINE PSPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     $                    DESCB, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      REAL               A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PSPOTRS solves a system of linear equations
*
*                      sub( A ) * X = sub( B )
*          A(IA:IA+N-1,JA:JA+N-1)*X = B(IB:IB+N-1,JB:JB+NRHS-1)
*
*  where sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1) and is a N-by-N
*  symmetric positive definite distributed matrix using the Cholesky
*  factorization sub( A ) = U**T*U or L*L**T computed by PSPOTRF.
*  sub( B ) denotes the distributed matrix B(IB:IB+N-1,JB:JB+NRHS-1).
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
*  UPLO    (global input) CHARACTER
*          = 'U':  Upper triangle of sub( A ) is stored;
*          = 'L':  Lower triangle of sub( A ) is stored.
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the distributed submatrix sub( B ).  NRHS >= 0.
*
*  A       (local input) REAL pointer into local memory to
*          an array of dimension (LLD_A, LOCc(JA+N-1)). On entry, this
*          array contains the factors L or U from the Cholesky facto-
*          rization sub( A ) = L*L**T or U**T*U, as computed by PSPOTRF.
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
*  B       (local input/local output) REAL pointer into the
*          local memory to an array of local dimension
*          (LLD_B,LOCc(JB+NRHS-1)).  On entry, this array contains the
*          the local pieces of the right hand sides sub( B ).
*          On exit, this array contains the local pieces of the solution
*          distributed matrix X.
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
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            IAROW, IBROW, ICTXT, IROFFA, IROFFB, ICOFFA,
     $                   MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK2MAT, PSTRSM,
     $                   PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P
      EXTERNAL           INDXG2P, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters.
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters.
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(700+CTXT_)
      ELSE
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 2, NRHS, 3, IB, JB, DESCB, 11, INFO )
         UPPER = LSAME( UPLO, 'U' )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $                       NPROW )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IF ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(700+NB_)
            ELSE IF( IROFFB.NE.0 .OR. IBROW.NE.IAROW ) THEN
               INFO = -9
            ELSE IF( DESCB( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(1100+NB_)
            END IF
         END IF
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 1
         CALL PCHK2MAT( N, 2, N, 2, IA, JA, DESCA, 7, N, 2, NRHS,
     $                  3, IB, JB, DESCB, 11, 1, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PSPOTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Solve sub( A ) * X = sub( B ) where sub( A ) = U'*U.
*
*        Solve U'*X = sub( B ), overwriting sub( B ) with X.
*
         CALL PSTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $              ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
*
*        Solve U*X = sub( B ), overwriting sub( B ) with X.
*
         CALL PSTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $              NRHS, ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
      ELSE
*
*        Solve sub( A ) *X = sub( B ) where sub( A ) = L*L'.
*
*        Solve L*X = sub( B ), overwriting sub( B ) with X.
*
         CALL PSTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', N,
     $              NRHS, ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
*
*        Solve L'*X = sub( B ), overwriting sub( B ) with X.
*
         CALL PSTRSM( 'Left', 'Lower', 'Transpose', 'Non-unit', N, NRHS,
     $              ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
      END IF
*
      RETURN
*
*     End of PSPOTRS
*
      END
