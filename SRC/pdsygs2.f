*
*
      SUBROUTINE PDSYGS2( IBTYPE, UPLO, N, A, IA, JA, DESCA, B, IB, JB,
     $                    DESCB, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, IB, IBTYPE, INFO, JA, JB, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PDSYGS2 reduces a real symmetric-definite generalized eigenproblem
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
     $                   ICTXT, IIA, IIB, IOFFA, IOFFB, IROFFA, IROFFB,
     $                   JJA, JJB, K, LDA, LDB, MYCOL, MYROW, NPCOL,
     $                   NPROW
      DOUBLE PRECISION   AKK, BKK, CT
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GRIDINFO, CHK1MAT, DAXPY,
     $                   DSCAL, DSYR2, DTRMV, DTRSV, INFOG2L, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P
      EXTERNAL           LSAME, INDXG2P
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters.
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
            ELSE IF( N+ICOFFA.GT.DESCA( NB_ ) ) THEN
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
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDSYGS2', -INFO )
         CALL BLACS_EXIT( ICTXT )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. ( MYROW.NE.IAROW .OR. MYCOL.NE.IACOL ) )
     $   RETURN
*
*     Compute local information
*
      LDA = DESCA( LLD_ )
      LDB = DESCB( LLD_ )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIB, JJB,
     $              IBROW, IBCOL )
*
      IF( IBTYPE.EQ.1 ) THEN
*
         IF( UPPER ) THEN
*
            IOFFA = IIA + JJA*LDA
            IOFFB = IIB + JJB*LDB
*
*           Compute inv(U')*sub( A )*inv(U)
*
            DO 10 K = 1, N
*
*              Update the upper triangle of
*              A(ia+k-1:ia+n-a,ia+k-1:ia+n-1)
*
               AKK = A( IOFFA-LDA )
               BKK = B( IOFFB-LDB )
               AKK = AKK / BKK**2
               A( IOFFA-LDA ) = AKK
               IF( K.LT.N ) THEN
                  CALL DSCAL( N-K, ONE / BKK, A( IOFFA ), LDA )
                  CT = -HALF*AKK
                  CALL DAXPY( N-K, CT, B( IOFFB ), LDB, A( IOFFA ),
     $                        LDA )
                  CALL DSYR2( UPLO, N-K, -ONE, A( IOFFA ), LDA,
     $                        B( IOFFB ), LDB, A( IOFFA+1 ), LDA )
                  CALL DAXPY( N-K, CT, B( IOFFB ), LDB, A( IOFFA ),
     $                        LDA )
                  CALL DTRSV( UPLO, 'Transpose', 'Non-unit', N-K,
     $                        B( IOFFB+1 ), LDB, A( IOFFA ), LDA )
               END IF
*
*              A( IOFFA ) -> A( K, K+1 )
*              B( IOFFB ) -> B( K, K+1 )
*
               IOFFA = IOFFA + LDA + 1
               IOFFB = IOFFB + LDB + 1
*
   10       CONTINUE
*
         ELSE
*
            IOFFA = IIA + 1 + ( JJA-1 )*LDA
            IOFFB = IIB + 1 + ( JJB-1 )*LDB
*
*           Compute inv(L)*sub( A )*inv(L')
*
            DO 20 K = 1, N
*
*              Update the lower triangle of
*              A(ia+k-1:ia+n-a,ia+k-1:ia+n-1)
*
               AKK = A( IOFFA-1 )
               BKK = B( IOFFB-1 )
               AKK = AKK / BKK**2
               A( IOFFA-1 ) = AKK
*
               IF( K.LT.N ) THEN
                  CALL DSCAL( N-K, ONE / BKK, A( IOFFA ), 1 )
                  CT = -HALF*AKK
                  CALL DAXPY( N-K, CT, B( IOFFB ), 1, A( IOFFA ), 1 )
                  CALL DSYR2( UPLO, N-K, -ONE, A( IOFFA ), 1,
     $                        B( IOFFB ), 1, A( IOFFA+LDA ), LDA )
                  CALL DAXPY( N-K, CT, B( IOFFB ), 1, A( IOFFA ), 1 )
                  CALL DTRSV( UPLO, 'No transpose', 'Non-unit', N-K,
     $                        B( IOFFB+LDB ), LDB, A( IOFFA ), 1 )
               END IF
*
*              A( IOFFA ) -> A( K+1, K )
*              B( IOFFB ) -> B( K+1, K )
*
               IOFFA = IOFFA + LDA + 1
               IOFFB = IOFFB + LDB + 1
*
   20       CONTINUE
*
         END IF
*
      ELSE
*
         IF( UPPER ) THEN
*
            IOFFA = IIA + ( JJA-1 )*LDA
            IOFFB = IIB + ( JJB-1 )*LDB
*
*           Compute U*sub( A )*U'
*
            DO 30 K = 1, N
*
*              Update the upper triangle of A(ia:ia+k-1,ja:ja+k-1)
*
               AKK = A( IOFFA+K-1 )
               BKK = B( IOFFB+K-1 )
               CALL DTRMV( UPLO, 'No transpose', 'Non-unit', K-1,
     $                     B( IIB+( JJB-1 )*LDB ), LDB, A( IOFFA ), 1 )
               CT = HALF*AKK
               CALL DAXPY( K-1, CT, B( IOFFB ), 1, A( IOFFA ), 1 )
               CALL DSYR2( UPLO, K-1, ONE, A( IOFFA ), 1, B( IOFFB ), 1,
     $                     A( IIA+( JJA-1 )*LDA ), LDA )
               CALL DAXPY( K-1, CT, B( IOFFB ), 1, A( IOFFA ), 1 )
               CALL DSCAL( K-1, BKK, A( IOFFA ), 1 )
               A( IOFFA+K-1 ) = AKK*BKK**2
*
*              A( IOFFA ) -> A( 1, K )
*              B( IOFFB ) -> B( 1, K )
*
               IOFFA = IOFFA + LDA
               IOFFB = IOFFB + LDB
*
   30       CONTINUE
*
         ELSE
*
            IOFFA = IIA + ( JJA-1 )*LDA
            IOFFB = IIB + ( JJB-1 )*LDB
*
*           Compute L'*sub( A )*L
*
            DO 40 K = 1, N
*
*              Update the lower triangle of A(ia:ia+k-1,ja:ja+k-1)
*
               AKK = A( IOFFA+( K-1 )*LDA )
               BKK = B( IOFFB+( K-1 )*LDB )
               CALL DTRMV( UPLO, 'Transpose', 'Non-unit', K-1,
     $                     B( IIB+( JJB-1 )*LDB ), LDB, A( IOFFA ),
     $                     LDA )
               CT = HALF*AKK
               CALL DAXPY( K-1, CT, B( IOFFB ), LDB, A( IOFFA ), LDA )
               CALL DSYR2( UPLO, K-1, ONE, A( IOFFA ), LDA, B( IOFFB ),
     $                     LDB, A( IIA+( JJA-1 )*LDA ), LDA )
               CALL DAXPY( K-1, CT, B( IOFFB ), LDB, A( IOFFA ), LDA )
               CALL DSCAL( K-1, BKK, A( IOFFA ), LDA )
               A( IOFFA+( K-1 )*LDA ) = AKK*BKK**2
*
*              A( IOFFA ) -> A( K, 1 )
*              B( IOFFB ) -> B( K, 1 )
*
               IOFFA = IOFFA + 1
               IOFFB = IOFFB + 1
*
   40       CONTINUE
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PDSYGS2
*
      END
