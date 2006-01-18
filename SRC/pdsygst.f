*
*
      SUBROUTINE PDSYGST( IBTYPE, UPLO, N, A, IA, JA, DESCA, B, IB, JB,
     $                    DESCB, SCALE, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, IB, IBTYPE, INFO, JA, JB, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PDSYGST reduces a real symmetric-definite generalized eigenproblem
*  to standard form.
*
*  In the following sub( A ) denotes A( IA:IA+N-1, JA:JA+N-1 ) and
*  sub( B ) denotes B( IB:IB+N-1, JB:JB+N-1 ).
*
*  If IBTYPE = 1, the problem is sub( A )*x = lambda*sub( B )*x,
*  and sub( A ) is overwritten by inv(U**T)*sub( A )*inv(U) or
*  inv(L)*sub( A )*inv(L**T)
*
*  If IBTYPE = 2 or 3, the problem is sub( A )*sub( B )*x = lambda*x or
*  sub( B )*sub( A )*x = lambda*x, and sub( A ) is overwritten by
*  U*sub( A )*U**T or L**T*sub( A )*L.
*
*  sub( B ) must have been previously factorized as U**T*U or L*L**T by
*  PDPOTRF.
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
*  Arguments
*  =========
*
*  IBTYPE   (global input) INTEGER
*          = 1: compute inv(U**T)*sub( A )*inv(U) or
*               inv(L)*sub( A )*inv(L**T);
*          = 2 or 3: compute U*sub( A )*U**T or L**T*sub( A )*L.
*
*  UPLO    (global input) CHARACTER
*          = 'U':  Upper triangle of sub( A ) is stored and sub( B ) is
*                  factored as U**T*U;
*          = 'L':  Lower triangle of sub( A ) is stored and sub( B ) is
*                  factored as L*L**T.
*
*  N       (global input) INTEGER
*          The order of the matrices sub( A ) and sub( B ).  N >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          N-by-N symmetric distributed matrix sub( A ). If UPLO = 'U',
*          the leading N-by-N upper triangular part of sub( A ) contains
*          the upper triangular part of the matrix, and its strictly
*          lower triangular part is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of sub( A ) contains
*          the lower triangular part of the matrix, and its strictly
*          upper triangular part is not referenced.
*
*          On exit, if INFO = 0, the transformed matrix, stored in the
*          same format as sub( A ).
*
*  IA      (global input) INTEGER
*          A's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  B       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_B, LOCc(JB+N-1)). On entry,
*          this array contains the local pieces of the triangular factor
*          from the Cholesky factorization of sub( B ), as returned by
*          PDPOTRF.
*
*  IB      (global input) INTEGER
*          B's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JB      (global input) INTEGER
*          B's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
*
*  SCALE   (global output) DOUBLE PRECISION
*          Amount by which the eigenvalues should be scaled to
*          compensate for the scaling performed in this routine.
*          At present, SCALE is always returned as 1.0, it is
*          returned here to allow for future enhancement.
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
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D+0, HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            IACOL, IAROW, IBCOL, IBROW, ICOFFA, ICOFFB,
     $                   ICTXT, IROFFA, IROFFB, K, KB, MYCOL, MYROW, NB,
     $                   NPCOL, NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 2 ), IDUM2( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK2MAT, PDSYGS2,
     $                   PDSYMM, PDSYR2K, PDTRMM, PDTRSM, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, MIN, MOD
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, INDXG2P
      EXTERNAL           LSAME, ICEIL, INDXG2P
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Get grid parameters
*
      SCALE = ONE
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 700+CTXT_ )
      ELSE
         UPPER = LSAME( UPLO, 'U' )
         CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 3, N, 3, IB, JB, DESCB, 11, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $              NPROW )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $              NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $              NPCOL )
            IBCOL = INDXG2P( JB, DESCB( NB_ ), MYCOL, DESCB( CSRC_ ),
     $              NPCOL )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            ICOFFB = MOD( JB-1, DESCB( NB_ ) )
            IF( IBTYPE.LT.1 .OR. IBTYPE.GT.3 ) THEN
               INFO = -1
            ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -2
            ELSE IF( N.LT.0 ) THEN
               INFO = -3
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            ELSE IF( IROFFB.NE.0 .OR. IBROW.NE.IAROW ) THEN
               INFO = -9
            ELSE IF( ICOFFB.NE.0 .OR. IBCOL.NE.IACOL ) THEN
               INFO = -10
            ELSE IF( DESCB( MB_ ).NE.DESCA( MB_ ) ) THEN
               INFO = -( 1100+MB_ )
            ELSE IF( DESCB( NB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 1100+NB_ )
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -( 1100+CTXT_ )
            END IF
         END IF
         IDUM1( 1 ) = IBTYPE
         IDUM2( 1 ) = 1
         IF( UPPER ) THEN
            IDUM1( 2 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'L' )
         END IF
         IDUM2( 2 ) = 2
         CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 7, N, 3, N, 3, IB,
     $                  JB, DESCB, 11, 2, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDSYGST', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( IBTYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
*
*           Compute inv(U')*sub( A )*inv(U)
*
            K = 1
            NB = DESCA( NB_ )
            KB = MIN( ICEIL( JA, NB )*NB, JA+N-1 ) - JA + 1
*
   10       CONTINUE
*
*           Update the upper triangle of A(ia+k-1:ia+n-1,ja+k-1:ja+n-1)
*
            CALL PDSYGS2( IBTYPE, UPLO, KB, A, IA+K-1, JA+K-1, DESCA, B,
     $                    IB+K-1, IB+K-1, DESCB, INFO )
            IF( K+KB.LE.N ) THEN
               CALL PDTRSM( 'Left', UPLO, 'Transpose', 'Non-unit', KB,
     $                      N-K-KB+1, ONE, B, IB+K-1, JB+K-1, DESCB, A,
     $                      IA+K-1, JA+K+KB-1, DESCA )
               CALL PDSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF, A,
     $                      IA+K-1, JA+K-1, DESCA, B, IB+K-1, JB+K+KB-1,
     $                      DESCB, ONE, A, IA+K-1, JA+K+KB-1, DESCA )
               CALL PDSYR2K( UPLO, 'Transpose', N-K-KB+1, KB, -ONE, A,
     $                       IA+K-1, JA+K+KB-1, DESCA, B, IB+K-1,
     $                       JB+K+KB-1, DESCB, ONE, A, IA+K+KB-1,
     $                       JA+K+KB-1, DESCA )
               CALL PDSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF, A,
     $                      IA+K-1, JA+K-1, DESCA, B, IB+K-1, JB+K+KB-1,
     $                      DESCB, ONE, A, IA+K-1, JA+K+KB-1, DESCA )
               CALL PDTRSM( 'Right', UPLO, 'No transpose', 'Non-unit',
     $                      KB, N-K-KB+1, ONE, B, IB+K+KB-1, JB+K+KB-1,
     $                      DESCB, A, IA+K-1, JA+K+KB-1, DESCA )
            END IF
            K = K + KB
            KB = MIN( N-K+1, NB )
*
            IF( K.LE.N )
     $         GO TO 10
*
         ELSE
*
*           Compute inv(L)*sub( A )*inv(L')
*
            K = 1
            NB = DESCA( MB_ )
            KB = MIN( ICEIL( IA, NB )*NB, IA+N-1 ) - IA + 1
*
   20       CONTINUE
*
*           Update the lower triangle of A(ia+k-1:ia+n-1,ja+k-1:ja+n-1)
*
            CALL PDSYGS2( IBTYPE, UPLO, KB, A, IA+K-1, JA+K-1, DESCA, B,
     $                    IB+K-1, JB+K-1, DESCB, INFO )
            IF( K+KB.LE.N ) THEN
               CALL PDTRSM( 'Right', UPLO, 'Transpose', 'Non-unit',
     $                      N-K-KB+1, KB, ONE, B, IB+K-1, JB+K-1, DESCB,
     $                      A, IA+K+KB-1, JA+K-1, DESCA )
               CALL PDSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF, A,
     $                      IA+K-1, JA+K-1, DESCA, B, IB+K+KB-1, JB+K-1,
     $                      DESCB, ONE, A, IA+K+KB-1, JA+K-1, DESCA )
               CALL PDSYR2K( UPLO, 'No transpose', N-K-KB+1, KB, -ONE,
     $                       A, IA+K+KB-1, JA+K-1, DESCA, B, IB+K+KB-1,
     $                       JB+K-1, DESCB, ONE, A, IA+K+KB-1,
     $                       JA+K+KB-1, DESCA )
               CALL PDSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF, A,
     $                      IA+K-1, JA+K-1, DESCA, B, IB+K+KB-1, JB+K-1,
     $                      DESCB, ONE, A, IA+K+KB-1, JA+K-1, DESCA )
               CALL PDTRSM( 'Left', UPLO, 'No transpose', 'Non-unit',
     $                      N-K-KB+1, KB, ONE, B, IB+K+KB-1, JB+K+KB-1,
     $                      DESCB, A, IA+K+KB-1, JA+K-1, DESCA )
            END IF
            K = K + KB
            KB = MIN( N-K+1, NB )
*
            IF( K.LE.N )
     $         GO TO 20
*
         END IF
*
      ELSE
*
         IF( UPPER ) THEN
*
*           Compute U*sub( A )*U'
*
            K = 1
            NB = DESCA( NB_ )
            KB = MIN( ICEIL( JA, NB )*NB, JA+N-1 ) - JA + 1
*
   30       CONTINUE
*
*           Update the upper triangle of A(ia:ia+k+kb-2,ja:ja+k+kb-2)
*
            CALL PDTRMM( 'Left', UPLO, 'No transpose', 'Non-unit', K-1,
     $                   KB, ONE, B, IB, JB, DESCB, A, IA, JA+K-1,
     $                   DESCA )
            CALL PDSYMM( 'Right', UPLO, K-1, KB, HALF, A, IA+K-1,
     $                   JA+K-1, DESCA, B, IB, JB+K-1, DESCB, ONE, A,
     $                   IA, JA+K-1, DESCA )
            CALL PDSYR2K( UPLO, 'No transpose', K-1, KB, ONE, A, IA,
     $                    JA+K-1, DESCA, B, IB, JB+K-1, DESCB, ONE, A,
     $                    IA, JA, DESCA )
            CALL PDSYMM( 'Right', UPLO, K-1, KB, HALF, A, IA+K-1,
     $                   JA+K-1, DESCA, B, IB, JB+K-1, DESCB, ONE, A,
     $                   IA, JA+K-1, DESCA )
            CALL PDTRMM( 'Right', UPLO, 'Transpose', 'Non-unit', K-1,
     $                   KB, ONE, B, IB+K-1, JB+K-1, DESCB, A, IA,
     $                   JA+K-1, DESCA )
            CALL PDSYGS2( IBTYPE, UPLO, KB, A, IA+K-1, JA+K-1, DESCA, B,
     $                    IB+K-1, JB+K-1, DESCB, INFO )
*
            K = K + KB
            KB = MIN( N-K+1, NB )
*
            IF( K.LE.N )
     $         GO TO 30
*
         ELSE
*
*           Compute L'*sub( A )*L
*
            K = 1
            NB = DESCA( MB_ )
            KB = MIN( ICEIL( IA, NB )*NB, IA+N-1 ) - IA + 1
*
   40       CONTINUE
*
*           Update the lower triangle of A(ia:ia+k+kb-2,ja:ja+k+kb-2)
*
            CALL PDTRMM( 'Right', UPLO, 'No transpose', 'Non-unit', KB,
     $                   K-1, ONE, B, IB, JB, DESCB, A, IA+K-1, JA,
     $                   DESCA )
            CALL PDSYMM( 'Left', UPLO, KB, K-1, HALF, A, IA+K-1, JA+K-1,
     $                   DESCA, B, IB+K-1, JB, DESCB, ONE, A, IA+K-1,
     $                   JA, DESCA )
            CALL PDSYR2K( UPLO, 'Transpose', K-1, KB, ONE, A, IA+K-1,
     $                    JA, DESCA, B, IB+K-1, JB, DESCB, ONE, A, IA,
     $                    JA, DESCA )
            CALL PDSYMM( 'Left', UPLO, KB, K-1, HALF, A, IA+K-1, JA+K-1,
     $                   DESCA, B, IB+K-1, JB, DESCB, ONE, A, IA+K-1,
     $                   JA, DESCA )
            CALL PDTRMM( 'Left', UPLO, 'Transpose', 'Non-unit', KB, K-1,
     $                   ONE, B, IB+K-1, JB+K-1, DESCB, A, IA+K-1, JA,
     $                   DESCA )
            CALL PDSYGS2( IBTYPE, UPLO, KB, A, IA+K-1, JA+K-1, DESCA, B,
     $                    IB+K-1, JB+K-1, DESCB, INFO )
*
            K = K + KB
            KB = MIN( N-K+1, NB )
*
            IF( K.LE.N )
     $         GO TO 40
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PDSYGST
*
      END
