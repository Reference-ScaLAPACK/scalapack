      SUBROUTINE PZHENGST( IBTYPE, UPLO, N, A, IA, JA, DESCA, B, IB, JB,
     $                     DESCB, SCALE, WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     October 15, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, IB, IBTYPE, INFO, JA, JB, LWORK, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX*16         A( * ), B( * ), WORK( * )
*     ..
*
*  Purpose
*
*  =======
*
*  PZHENGST reduces a complex Hermitian-definite generalized
*  eigenproblem to standard form.
*
*  PZHENGST performs the same function as PZHEGST, but is based on
*  rank 2K updates, which are faster and more scalable than
*  triangular solves (the basis of PZHENGST).
*
*  PZHENGST calls PZHEGST when UPLO='U', hence PZHENGST provides
*  improved performance only when UPLO='L', IBTYPE=1.
*
*  PZHENGST also calls PZHEGST when insufficient workspace is
*  provided,  hence PZHENGST provides improved
*  performance only when LWORK >= 2 * NP0 * NB + NQ0 * NB + NB * NB
*
*  In the following sub( A ) denotes A( IA:IA+N-1, JA:JA+N-1 ) and
*  sub( B ) denotes B( IB:IB+N-1, JB:JB+N-1 ).
*
*  If IBTYPE = 1, the problem is sub( A )*x = lambda*sub( B )*x,
*  and sub( A ) is overwritten by inv(U**H)*sub( A )*inv(U) or
*  inv(L)*sub( A )*inv(L**H)
*
*  If IBTYPE = 2 or 3, the problem is sub( A )*sub( B )*x = lambda*x or
*  sub( B )*sub( A )*x = lambda*x, and sub( A ) is overwritten by
*  U*sub( A )*U**H or L**H*sub( A )*L.
*
*  sub( B ) must have been previously factorized as U**H*U or L*L**H by
*  PZPOTRF.
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
*          = 1: compute inv(U**H)*sub( A )*inv(U) or
*               inv(L)*sub( A )*inv(L**H);
*          = 2 or 3: compute U*sub( A )*U**H or L**H*sub( A )*L.
*
*  UPLO    (global input) CHARACTER
*          = 'U':  Upper triangle of sub( A ) is stored and sub( B ) is
*                  factored as U**H*U;
*          = 'L':  Lower triangle of sub( A ) is stored and sub( B ) is
*                  factored as L*L**H.
*
*  N       (global input) INTEGER
*          The order of the matrices sub( A ) and sub( B ).  N >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          N-by-N Hermitian distributed matrix sub( A ). If UPLO = 'U',
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
*  B       (local input) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_B, LOCc(JB+N-1)). On entry,
*          this array contains the local pieces of the triangular factor
*          from the Cholesky factorization of sub( B ), as returned by
*          PZPOTRF.
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
*  WORK    (local workspace/local output) COMPLEX*16 array,
*                                                  dimension (LWORK)
*          On exit, WORK( 1 ) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= MAX( NB * ( NP0 +1 ), 3 * NB )
*
*          When IBTYPE = 1 and UPLO = 'L', PZHENGST provides improved
*          performance when LWORK >= 2 * NP0 * NB + NQ0 * NB + NB * NB
*
*          where NB = MB_A = NB_A,
*          NP0 = NUMROC( N, NB, 0, 0, NPROW ),
*          NQ0 = NUMROC( N, NB, 0, 0, NPROW ),
*
*          NUMROC ia a ScaLAPACK tool functions
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the
*          optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
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
*
*     .. Parameters ..
      COMPLEX*16         ONEHALF, ONE, MONE
      DOUBLE PRECISION   RONE
      PARAMETER          ( ONEHALF = ( 0.5D0, 0.0D0 ),
     $                   ONE = ( 1.0D0, 0.0D0 ),
     $                   MONE = ( -1.0D0, 0.0D0 ), RONE = 1.0D0 )
      INTEGER            DLEN_, CTXT_, MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( DLEN_ = 9, CTXT_ = 2, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IACOL, IAROW, IBCOL, IBROW, ICOFFA, ICOFFB,
     $                   ICTXT, INDAA, INDG, INDR, INDRT, IROFFA,
     $                   IROFFB, J, K, KB, LWMIN, LWOPT, MYCOL, MYROW,
     $                   NB, NP0, NPCOL, NPK, NPROW, NQ0, POSTK
*     ..
*     .. Local Arrays ..
      INTEGER            DESCAA( DLEN_ ), DESCG( DLEN_ ),
     $                   DESCR( DLEN_ ), DESCRT( DLEN_ ), IDUM1( 2 ),
     $                   IDUM2( 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           LSAME, INDXG2P, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, PCHK2MAT,
     $                   PXERBLA, PZGEMM, PZHEGST, PZHEMM, PZHER2K,
     $                   PZLACPY, PZTRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DCONJG, ICHAR, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      SCALE = 1.0D0
*
      NB = DESCA( MB_ )
*
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
            NP0 = NUMROC( N, NB, 0, 0, NPROW )
            NQ0 = NUMROC( N, NB, 0, 0, NPCOL )
            LWMIN = MAX( NB*( NP0+1 ), 3*NB )
            IF( IBTYPE.EQ.1 .AND. .NOT.UPPER ) THEN
               LWOPT = 2*NP0*NB + NQ0*NB + NB*NB
            ELSE
               LWOPT = LWMIN
            END IF
            WORK( 1 ) = DCMPLX( DBLE( LWOPT ) )
            LQUERY = ( LWORK.EQ.-1 )
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
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -13
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
         CALL PXERBLA( ICTXT, 'PZHENGST', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*
      IF( IBTYPE.NE.1 .OR. UPPER .OR. LWORK.LT.LWOPT ) THEN
         CALL PZHEGST( IBTYPE, UPLO, N, A, IA, JA, DESCA, B, IB, JB,
     $                 DESCB, SCALE, INFO )
         RETURN
      END IF
*
      CALL DESCSET( DESCG, N, NB, NB, NB, IAROW, IACOL, ICTXT, NP0 )
      CALL DESCSET( DESCR, N, NB, NB, NB, IAROW, IACOL, ICTXT, NP0 )
      CALL DESCSET( DESCRT, NB, N, NB, NB, IAROW, IACOL, ICTXT, NB )
      CALL DESCSET( DESCAA, NB, NB, NB, NB, IAROW, IACOL, ICTXT, NB )
*
      INDG = 1
      INDR = INDG + DESCG( LLD_ )*NB
      INDAA = INDR + DESCR( LLD_ )*NB
      INDRT = INDAA + DESCAA( LLD_ )*NB
*
      DO 30 K = 1, N, NB
*
         KB = MIN( N-K+1, NB )
         POSTK = K + KB
         NPK = N - POSTK + 1
*
*
         CALL PZLACPY( 'A', N-POSTK+1, KB, B, POSTK+IB-1, K+JB-1, DESCB,
     $                 WORK( INDG ), POSTK, 1, DESCG )
         CALL PZLACPY( 'A', N-POSTK+1, KB, A, POSTK+IA-1, K+JA-1, DESCA,
     $                 WORK( INDR ), POSTK, 1, DESCR )
         CALL PZLACPY( 'A', KB, K-1, A, K+IA-1, JA, DESCA,
     $                 WORK( INDRT ), 1, 1, DESCRT )
*
         CALL PZLACPY( 'L', KB, KB, A, K+IA-1, K+JA-1, DESCA,
     $                 WORK( INDR ), K, 1, DESCR )
         CALL PZTRSM( 'Right', 'L', 'N', 'N', NPK, KB, MONE, B, K+IB-1,
     $                K+JB-1, DESCB, WORK( INDG ), POSTK, 1, DESCG )
*
         CALL PZHEMM( 'Right', 'L', NPK, KB, ONEHALF, A, K+IA-1, K+JA-1,
     $                DESCA, WORK( INDG ), POSTK, 1, DESCG, ONE,
     $                WORK( INDR ), POSTK, 1, DESCR )
*
         CALL PZHER2K( 'Lower', 'No T', NPK, KB, ONE, WORK( INDG ),
     $                 POSTK, 1, DESCG, WORK( INDR ), POSTK, 1, DESCR,
     $                 RONE, A, POSTK+IA-1, POSTK+JA-1, DESCA )
*
         CALL PZGEMM( 'No T', 'No Conj', NPK, K-1, KB, ONE,
     $                WORK( INDG ), POSTK, 1, DESCG, WORK( INDRT ), 1,
     $                1, DESCRT, ONE, A, POSTK+IA-1, JA, DESCA )
*
         CALL PZHEMM( 'Right', 'L', NPK, KB, ONE, WORK( INDR ), K, 1,
     $                DESCR, WORK( INDG ), POSTK, 1, DESCG, ONE, A,
     $                POSTK+IA-1, K+JA-1, DESCA )
*
         CALL PZTRSM( 'Left', 'Lower', 'No Conj', 'Non-unit', KB, K-1,
     $                ONE, B, K+IB-1, K+JB-1, DESCB, A, K+IA-1, JA,
     $                DESCA )
*
         CALL PZLACPY( 'L', KB, KB, A, K+IA-1, K+JA-1, DESCA,
     $                 WORK( INDAA ), 1, 1, DESCAA )
*
         IF( MYROW.EQ.DESCAA( RSRC_ ) .AND. MYCOL.EQ.DESCAA( CSRC_ ) )
     $        THEN
            DO 20 I = 1, KB
               DO 10 J = 1, I
                  WORK( INDAA+J-1+( I-1 )*DESCAA( LLD_ ) )
     $               = DCONJG( WORK( INDAA+I-1+( J-1 )*
     $               DESCAA( LLD_ ) ) )
   10          CONTINUE
   20       CONTINUE
         END IF
*
         CALL PZTRSM( 'Left', 'Lower', 'No Conj', 'Non-unit', KB, KB,
     $                ONE, B, K+IB-1, K+JB-1, DESCB, WORK( INDAA ), 1,
     $                1, DESCAA )
*
         CALL PZTRSM( 'Right', 'Lower', 'Conj', 'Non-unit', KB, KB, ONE,
     $                B, K+IB-1, K+JB-1, DESCB, WORK( INDAA ), 1, 1,
     $                DESCAA )
*
         CALL PZLACPY( 'L', KB, KB, WORK( INDAA ), 1, 1, DESCAA, A,
     $                 K+IA-1, K+JA-1, DESCA )
*
         CALL PZTRSM( 'Right', 'Lower', 'Conj', 'Non-unit', NPK, KB,
     $                ONE, B, K+IB-1, K+JB-1, DESCB, A, POSTK+IA-1,
     $                K+JA-1, DESCA )
*
         DESCR( CSRC_ ) = MOD( DESCR( CSRC_ )+1, NPCOL )
         DESCG( CSRC_ ) = MOD( DESCG( CSRC_ )+1, NPCOL )
         DESCRT( RSRC_ ) = MOD( DESCRT( RSRC_ )+1, NPROW )
         DESCAA( RSRC_ ) = MOD( DESCAA( RSRC_ )+1, NPROW )
         DESCAA( CSRC_ ) = MOD( DESCAA( CSRC_ )+1, NPCOL )
   30 CONTINUE
*
      WORK( 1 ) = DCMPLX( DBLE( LWOPT ) )
*
      RETURN
      END
