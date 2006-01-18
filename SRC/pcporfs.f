      SUBROUTINE PCPORFS( UPLO, N, NRHS, A, IA, JA, DESCA, AF, IAF, JAF,
     $                    DESCAF, B, IB, JB, DESCB, X, IX, JX, DESCX,
     $                    FERR, BERR, WORK, LWORK, RWORK, LRWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, IAF, IB, INFO, IX, JA, JAF, JB, JX,
     $                   LRWORK, LWORK, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCAF( * ), DESCB( * ),
     $                   DESCX( * )
      COMPLEX            A( * ), AF( * ), B( * ), WORK( * ), X( * )
      REAL               BERR( * ), FERR( * ), RWORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PCPORFS improves the computed solution to a system of linear
*  equations when the coefficient matrix is Hermitian positive definite
*  and provides error bounds and backward error estimates for the
*  solutions.
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
*  In the following comments, sub( A ), sub( X ) and sub( B ) denote
*  respectively A(IA:IA+N-1,JA:JA+N-1), X(IX:IX+N-1,JX:JX+NRHS-1) and
*  B(IB:IB+N-1,JB:JB+NRHS-1).
*
*  Arguments
*  =========
*
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix sub( A ) is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The order of the matrix sub( A ).  N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrices sub( B ) and sub( X ).  NRHS >= 0.
*
*  A       (local input) COMPLEX pointer into the local
*          memory to an array of local dimension (LLD_A,LOCc(JA+N-1) ).
*          This array contains the local pieces of the N-by-N Hermitian
*          distributed matrix sub( A ) to be factored.
*          If UPLO = 'U', the leading N-by-N upper triangular part of
*          sub( A ) contains the upper triangular part of the matrix,
*          and its strictly lower triangular part is not referenced.
*          If UPLO = 'L', the leading N-by-N lower triangular part of
*          sub( A ) contains the lower triangular part of the distribu-
*          ted matrix, and its strictly upper triangular part is not
*          referenced.
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
*  AF      (local input) COMPLEX pointer into the local memory
*          to an array of local dimension (LLD_AF,LOCc(JA+N-1)).
*          On entry, this array contains the factors L or U from the
*          Cholesky factorization sub( A ) = L*L**H or U**H*U, as
*          computed by PCPOTRF.
*
*  IAF     (global input) INTEGER
*          The row index in the global array AF indicating the first
*          row of sub( AF ).
*
*  JAF     (global input) INTEGER
*          The column index in the global array AF indicating the
*          first column of sub( AF ).
*
*  DESCAF  (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix AF.
*
*  B       (local input) COMPLEX pointer into the local memory
*          to an array of local dimension (LLD_B, LOCc(JB+NRHS-1) ).
*          On entry, this array contains the the local pieces of the
*          right hand sides sub( B ).
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
*  X       (local input) COMPLEX pointer into the local memory
*          to an array of local dimension (LLD_X, LOCc(JX+NRHS-1) ).
*          On entry, this array contains the the local pieces of the
*          solution vectors sub( X ). On exit, it contains the
*          improved solution vectors.
*
*  IX      (global input) INTEGER
*          The row index in the global array X indicating the first
*          row of sub( X ).
*
*  JX      (global input) INTEGER
*          The column index in the global array X indicating the
*          first column of sub( X ).
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*
*  FERR    (local output) REAL array of local dimension
*          LOCc(JB+NRHS-1).
*          The estimated forward error bound for each solution vector
*          of sub( X ).  If XTRUE is the true solution corresponding
*          to sub( X ), FERR is an estimated upper bound for the
*          magnitude of the largest element in (sub( X ) - XTRUE)
*          divided by the magnitude of the largest element in sub( X ).
*          The estimate is as reliable as the estimate for RCOND, and
*          is almost always a slight overestimate of the true error.
*          This array is tied to the distributed matrix X.
*
*  BERR    (local output) REAL array of local dimension
*          LOCc(JB+NRHS-1). The componentwise relative backward
*          error of each solution vector (i.e., the smallest re-
*          lative change in any entry of sub( A ) or sub( B )
*          that makes sub( X ) an exact solution).
*          This array is tied to the distributed matrix X.
*
*  WORK    (local workspace/local output) COMPLEX array,
*                                                   dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= 2*LOCr( N + MOD( IA-1, MB_A ) )
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  RWORK   (local workspace/local output) REAL array,
*                                                    dimension (LRWORK)
*          On exit, RWORK(1) returns the minimal and optimal LRWORK.
*
*  LRWORK  (local or global input) INTEGER
*          The dimension of the array RWORK.
*          LRWORK is local input and must be at least
*          LRWORK >= LOCr( N + MOD( IB-1, MB_B ) ).
*
*          If LRWORK = -1, then LRWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  Internal Parameters
*  ===================
*
*  ITMAX is the maximum number of steps of iterative refinement.
*
*  Notes
*  =====
*
*  This routine temporarily returns when N <= 1.
*
*  The distributed submatrices op( A ) and op( AF ) (respectively
*  sub( X ) and sub( B ) ) should be distributed the same way on the
*  same processes. These conditions ensure that sub( A ) and sub( AF )
*  (resp. sub( X ) and sub( B ) ) are "perfectly" aligned.
*
*  Moreover, this routine requires the distributed submatrices sub( A ),
*  sub( AF ), sub( X ), and sub( B ) to be aligned on a block boundary,
*  i.e., if f(x,y) = MOD( x-1, y ):
*  f( IA, DESCA( MB_ ) ) = f( JA, DESCA( NB_ ) ) = 0,
*  f( IAF, DESCAF( MB_ ) ) = f( JAF, DESCAF( NB_ ) ) = 0,
*  f( IB, DESCB( MB_ ) ) = f( JB, DESCB( NB_ ) ) = 0, and
*  f( IX, DESCX( MB_ ) ) = f( JX, DESCX( NB_ ) ) = 0.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      REAL               ZERO, RONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0E+0, RONE = 1.0E+0, TWO = 2.0E+0,
     $                     THREE = 3.0E+0 )
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            COUNT, IACOL, IAFCOL, IAFROW, IAROW, IXBCOL,
     $                   IXBROW, IXCOL, IXROW, ICOFFA, ICOFFAF, ICOFFB,
     $                   ICOFFX, ICTXT, ICURCOL, IDUM, II, IIXB, IIW,
     $                   IOFFXB, IPB, IPR, IPV, IROFFA, IROFFAF, IROFFB,
     $                   IROFFX, IW, J, JBRHS, JJ, JJFBE, JJXB, JN, JW,
     $                   K, KASE, LDXB, LRWMIN, LWMIN, MYCOL, MYRHS,
     $                   MYROW, NP, NP0, NPCOL, NPMOD, NPROW, NZ
      REAL               EPS, EST, LSTRES, S, SAFE1, SAFE2, SAFMIN
      COMPLEX            ZDUM
*     ..
*     .. Local Arrays ..
      INTEGER            DESCW( DLEN_ ), IDUM1( 5 ), IDUM2( 5 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, INDXG2P, NUMROC
      REAL               PSLAMCH
      EXTERNAL           ICEIL, INDXG2P, LSAME, NUMROC, PSLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CGEBR2D, CGEBS2D, CHK1MAT,
     $                   DESCSET, INFOG2L, PCAHEMV, PCAXPY, PCHK2MAT,
     $                   PCCOPY, PCHEMV, PCPOTRS, PCLACON,
     $                   PXERBLA, SGAMX2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, ICHAR, MAX, MIN, MOD, REAL
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
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
         INFO = -(700+CTXT_)
      ELSE
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 2, N, 2, IAF, JAF, DESCAF, 11, INFO )
         CALL CHK1MAT( N, 2, NRHS, 3, IB, JB, DESCB, 15, INFO )
         CALL CHK1MAT( N, 2, NRHS, 3, IX, JX, DESCX, 19, INFO )
         IF( INFO.EQ.0 ) THEN
            UPPER = LSAME( UPLO, 'U' )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFAF = MOD( IAF-1, DESCAF( MB_ ) )
            ICOFFAF = MOD( JAF-1, DESCAF( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            ICOFFB = MOD( JB-1, DESCB( NB_ ) )
            IROFFX = MOD( IX-1, DESCX( MB_ ) )
            ICOFFX = MOD( JX-1, DESCX( NB_ ) )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IAFCOL = INDXG2P( JAF, DESCAF( NB_ ), MYCOL,
     $                        DESCAF( CSRC_ ), NPCOL )
            IAFROW = INDXG2P( IAF, DESCAF( MB_ ), MYROW,
     $                        DESCAF( RSRC_ ), NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIXB, JJXB, IXBROW, IXBCOL )
            IXROW = INDXG2P( IX, DESCX( MB_ ), MYROW, DESCX( RSRC_ ),
     $                       NPROW )
            IXCOL = INDXG2P( JX, DESCX( NB_ ), MYCOL, DESCX( CSRC_ ),
     $                       NPCOL )
            NPMOD = NUMROC( N+IROFFA, DESCA( MB_ ), MYROW, IAROW,
     $                      NPROW )
            LWMIN = 2 * NPMOD
            LRWMIN = NPMOD
            WORK( 1 ) = CMPLX( REAL( LWMIN ) )
            RWORK( 1 ) = REAL( LRWMIN )
            LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
*
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( N.LT.0 ) THEN
               INFO = -2
            ELSE IF( NRHS.LT.0 ) THEN
               INFO = -3
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 700 + NB_ )
            ELSE IF( DESCA( MB_ ).NE.DESCAF( MB_ ) ) THEN
               INFO = -( 1100 + MB_ )
            ELSE IF( IROFFAF.NE.0 .OR. IAROW.NE.IAFROW ) THEN
               INFO = -9
            ELSE IF( DESCA( NB_ ).NE.DESCAF( NB_ ) ) THEN
               INFO = -( 1100 + NB_ )
            ELSE IF( ICOFFAF.NE.0 .OR. IACOL.NE.IAFCOL ) THEN
               INFO = -10
            ELSE IF( ICTXT.NE.DESCAF( CTXT_ ) ) THEN
               INFO = -( 1100 + CTXT_ )
            ELSE IF( IROFFA.NE.IROFFB .OR. IAROW.NE.IXBROW ) THEN
               INFO = -13
            ELSE IF( DESCA( MB_ ).NE.DESCB( MB_ ) ) THEN
               INFO = -( 1500 + MB_ )
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -( 1500 + CTXT_ )
            ELSE IF( DESCB( MB_ ).NE.DESCX( MB_ ) ) THEN
               INFO = -( 1900 + MB_ )
            ELSE IF( IROFFX.NE.0 .OR. IXBROW.NE.IXROW ) THEN
               INFO = -17
            ELSE IF( DESCB( NB_ ).NE.DESCX( NB_ ) ) THEN
               INFO = -( 1900 + NB_ )
            ELSE IF( ICOFFB.NE.ICOFFX .OR. IXBCOL.NE.IXCOL ) THEN
               INFO = -18
            ELSE IF( ICTXT.NE.DESCX( CTXT_ ) ) THEN
               INFO = -( 1900 + CTXT_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -23
            ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -25
            END IF
         END IF
*
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 1
         IDUM1( 2 ) = N
         IDUM2( 2 ) = 2
         IDUM1( 3 ) = NRHS
         IDUM2( 3 ) = 3
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 4 ) = -1
         ELSE
            IDUM1( 4 ) = 1
         END IF
         IDUM2( 4 ) = 23
         IF( LRWORK.EQ.-1 ) THEN
            IDUM1( 5 ) = -1
         ELSE
            IDUM1( 5 ) = 1
         END IF
         IDUM2( 5 ) = 25
         CALL PCHK2MAT( N, 2, N, 2, IA, JA, DESCA, 7, N, 2, N, 2, IAF,
     $                  JAF, DESCAF, 11, 0, IDUM1, IDUM2, INFO )
         CALL PCHK2MAT( N, 2, NRHS, 3, IB, JB, DESCB, 15, N, 2, NRHS, 3,
     $                  IX, JX, DESCX, 19, 5, IDUM1, IDUM2, INFO )
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCPORFS', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      JJFBE = JJXB
      MYRHS = NUMROC( JB+NRHS-1, DESCB( NB_ ), MYCOL, DESCB( CSRC_ ),
     $                NPCOL )
*
*     Quick return if possible
*
      IF( N.LE.1 .OR. NRHS.EQ.0 ) THEN
         DO 10 JJ = JJFBE, MYRHS
            FERR( JJ ) = ZERO
            BERR( JJ ) = ZERO
   10    CONTINUE
         RETURN
      END IF
*
      NP0 = NUMROC( N+IROFFB, DESCB( MB_ ), MYROW, IXBROW, NPROW )
      CALL DESCSET( DESCW, N+IROFFB, 1, DESCA( MB_ ), 1, IXBROW, IXBCOL,
     $              ICTXT, MAX( 1, NP0 ) )
      IPB = 1
      IPR = 1
      IPV = IPR + NP0
      IF( MYROW.EQ.IXBROW ) THEN
         IIW = 1 + IROFFB
         NP = NP0 - IROFFB
      ELSE
         IIW = 1
         NP = NP0
      END IF
      IW = 1 + IROFFB
      JW = 1
      LDXB = DESCB( LLD_ )
      IOFFXB = ( JJXB-1 )*LDXB
*
*     NZ = 1 + maximum number of nonzero entries in each row of sub( A )
*
      NZ = N + 1
      EPS = PSLAMCH( ICTXT, 'Epsilon' )
      SAFMIN = PSLAMCH( ICTXT, 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS
      JN = MIN( ICEIL( JB, DESCB( NB_ ) ) * DESCB( NB_ ), JB+NRHS-1 )
*
*     Handle first block separately
*
      JBRHS = JN - JB + 1
      DO 100 K = 0, JBRHS-1
*
         COUNT = 1
         LSTRES = THREE
   20    CONTINUE
*
*        Loop until stopping criterion is satisfied.
*
*        Compute residual R = sub(B) - op(sub(A)) * sub(X)
*
         CALL PCCOPY( N, B, IB, JB+K, DESCB, 1, WORK( IPR ), IW, JW,
     $                DESCW, 1 )
         CALL PCHEMV( UPLO, N, -ONE, A, IA, JA, DESCA, X, IX, JX+K,
     $                DESCX, 1, ONE, WORK( IPR ), IW, JW, DESCW, 1 )
*
*        Compute componentwise relative backward error from formula
*
*        max(i) ( abs(R(i))/(abs(sub(A))*abs(sub(X))+abs(sub(B)) )(i) )
*
*        where abs(Z) is the componentwise absolute value of the
*        matrix or vector Z.  If the i-th component of the
*        denominator is less than SAFE2, then SAFE1 is added to
*        the i-th components of the numerator and denominator
*        before dividing.
*
         IF( MYCOL.EQ.IXBCOL ) THEN
            IF( NP.GT.0 ) THEN
               DO 30 II = IIXB, IIXB + NP - 1
                  RWORK( IIW+II-IIXB ) = CABS1( B( II+IOFFXB ) )
   30          CONTINUE
            END IF
         END IF
*
         CALL PCAHEMV( UPLO, N, RONE, A, IA, JA, DESCA, X, IX, JX+K,
     $                 DESCX, 1, RONE, RWORK( IPB ), IW, JW, DESCW, 1 )
*
         S = ZERO
         IF( MYCOL.EQ.IXBCOL ) THEN
            IF( NP.GT.0 ) THEN
               DO 40 II = IIW-1, IIW+NP-2
                  IF( RWORK( IPB+II ).GT.SAFE2 ) THEN
                     S = MAX( S, CABS1( WORK( IPR+II ) ) /
     $                           RWORK( IPB+II ) )
                  ELSE
                     S = MAX( S, ( CABS1( WORK( IPR+II ) )+SAFE1 ) /
     $                           ( RWORK( IPB+II )+SAFE1 ) )
                  END IF
   40          CONTINUE
            END IF
         END IF
*
         CALL SGAMX2D( ICTXT, 'All', ' ', 1, 1, S, 1, IDUM, IDUM, 1,
     $                      -1, MYCOL )
         IF( MYCOL.EQ.IXBCOL )
     $      BERR( JJFBE ) = S
*
*        Test stopping criterion. Continue iterating if
*         1) The residual BERR(J) is larger than machine epsilon, and
*         2) BERR(J) decreased by at least a factor of 2 during the
*            last iteration, and
*         3) At most ITMAX iterations tried.
*
         IF( S.GT.EPS .AND. TWO*S.LE.LSTRES .AND. COUNT.LE.ITMAX ) THEN
*
*           Update solution and try again.
*
            CALL PCPOTRS( UPLO, N, 1, AF, IAF, JAF, DESCAF,
     $                    WORK( IPR ), IW, JW, DESCW, INFO )
            CALL PCAXPY( N, ONE, WORK( IPR ), IW, JW, DESCW, 1, X, IX,
     $                   JX+K, DESCX, 1 )
            LSTRES = S
            COUNT = COUNT + 1
            GO TO 20
         END IF
*
*        Bound error from formula
*
*        norm(sub(X) - XTRUE) / norm(sub(X)) .le. FERR =
*        norm( abs(inv(sub(A)))*
*            ( abs(R) +
*        NZ*EPS*( abs(sub(A))*abs(sub(X))+abs(sub(B)) ))) / norm(sub(X))
*
*        where
*          norm(Z) is the magnitude of the largest component of Z
*          inv(sub(A)) is the inverse of sub(A)
*          abs(Z) is the componentwise absolute value of the matrix
*          or vector Z
*          NZ is the maximum number of nonzeros in any row of sub(A),
*          plus 1
*          EPS is machine epsilon
*
*        The i-th component of
*               abs(R)+NZ*EPS*(abs(sub(A))*abs(sub(X))+abs(sub(B)))
*        is incremented by SAFE1 if the i-th component of
*        abs(sub(A))*abs(sub(X)) + abs(sub(B)) is less than SAFE2.
*
*        Use PCLACON to estimate the infinity-norm of the matrix
*        inv(sub(A)) * diag(W), where
*        W = abs(R) + NZ*EPS*( abs(sub(A))*abs(sub(X))+abs(sub(B)))))
*
         IF( MYCOL.EQ.IXBCOL ) THEN
            IF( NP.GT.0 ) THEN
               DO 50 II = IIW-1, IIW+NP-2
                  IF( RWORK( IPB+II ).GT.SAFE2 ) THEN
                     RWORK( IPB+II ) = CABS1( WORK( IPR+II ) ) +
     $                                 NZ*EPS*RWORK( IPB+II )
                  ELSE
                     RWORK( IPB+II ) = CABS1( WORK( IPR+II ) ) +
     $                                 NZ*EPS*RWORK( IPB+II ) + SAFE1
                  END IF
   50          CONTINUE
            END IF
         END IF
*
         KASE = 0
   60    CONTINUE
         IF( MYCOL.EQ.IXBCOL ) THEN
            CALL CGEBS2D( ICTXT, 'Rowwise', ' ', NP, 1, WORK( IPR ),
     $                    DESCW( LLD_ ) )
         ELSE
            CALL CGEBR2D( ICTXT, 'Rowwise', ' ', NP, 1, WORK( IPR ),
     $                    DESCW( LLD_ ), MYROW, IXBCOL )
         END IF
         DESCW( CSRC_ ) = MYCOL
         CALL PCLACON( N, WORK( IPV ), IW, JW, DESCW, WORK( IPR ),
     $                 IW, JW, DESCW, EST, KASE )
         DESCW( CSRC_ ) = IXBCOL
*
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
*
*              Multiply by diag(W)*inv(sub(A)').
*
               CALL PCPOTRS( UPLO, N, 1, AF, IAF, JAF, DESCAF,
     $                       WORK( IPR ), IW, JW, DESCW, INFO )
*
               IF( MYCOL.EQ.IXBCOL ) THEN
                  IF( NP.GT.0 ) THEN
                     DO 70 II = IIW-1, IIW+NP-2
                        WORK( IPR+II ) = RWORK( IPB+II )*WORK( IPR+II )
   70                CONTINUE
                  END IF
               END IF
            ELSE
*
*              Multiply by inv(sub(A))*diag(W).
*
               IF( MYCOL.EQ.IXBCOL ) THEN
                  IF( NP.GT.0 ) THEN
                     DO 80 II = IIW-1, IIW+NP-2
                        WORK( IPR+II ) = RWORK( IPB+II )*WORK( IPR+II )
   80                CONTINUE
                  END IF
               END IF
*
               CALL PCPOTRS( UPLO, N, 1, AF, IAF, JAF, DESCAF,
     $                       WORK( IPR ), IW, JW, DESCW, INFO )
            END IF
            GO TO 60
         END IF
*
*           Normalize error.
*
         LSTRES = ZERO
         IF( MYCOL.EQ.IXBCOL ) THEN
            IF( NP.GT.0 ) THEN
               DO 90 II = IIXB, IIXB+NP-1
                  LSTRES = MAX( LSTRES, CABS1( X( IOFFXB+II ) ) )
   90          CONTINUE
            END IF
            CALL SGAMX2D( ICTXT, 'Column', ' ', 1, 1, LSTRES, 1, IDUM,
     $                    IDUM, 1, -1, MYCOL )
            IF( LSTRES.NE.ZERO )
     $         FERR( JJFBE ) = EST / LSTRES
*
            JJXB = JJXB + 1
            JJFBE = JJFBE + 1
            IOFFXB = IOFFXB + LDXB
*
         END IF
*
  100 CONTINUE
*
      ICURCOL = MOD( IXBCOL+1, NPCOL )
*
*     Do for each right hand side
*
      DO 200 J = JN+1, JB+NRHS-1, DESCB( NB_ )
         JBRHS = MIN( JB+NRHS-J, DESCB( NB_ ) )
         DESCW( CSRC_ ) = ICURCOL
*
         DO 190 K = 0, JBRHS-1
*
            COUNT = 1
            LSTRES = THREE
  110       CONTINUE
*
*           Loop until stopping criterion is satisfied.
*
*           Compute residual R = sub( B ) - sub( A )*sub( X ).
*
            CALL PCCOPY( N, B, IB, J+K, DESCB, 1, WORK( IPR ), IW, JW,
     $                   DESCW, 1 )
            CALL PCHEMV( UPLO, N, -ONE, A, IA, JA, DESCA, X, IX, J+K,
     $                  DESCX, 1, ONE, WORK( IPR ), IW, JW, DESCW, 1 )
*
*           Compute componentwise relative backward error from formula
*
*           max(i) ( abs(R(i)) /
*                    ( abs(sub(A))*abs(sub(X)) + abs(sub(B)) )(i) )
*
*           where abs(Z) is the componentwise absolute value of the
*           matrix or vector Z.  If the i-th component of the
*           denominator is less than SAFE2, then SAFE1 is added to the
*           i-th components of the numerator and denominator before
*           dividing.
*
            IF( MYCOL.EQ.ICURCOL ) THEN
               IF( NP.GT.0 ) THEN
                  DO 120 II = IIXB, IIXB+NP-1
                     RWORK( IIW+II-IIXB ) = CABS1( B( II+IOFFXB ) )
  120             CONTINUE
               END IF
            END IF
*
            CALL PCAHEMV( UPLO, N, RONE, A, IA, JA, DESCA, X, IX, J+K,
     $                    DESCX, 1, RONE, RWORK( IPB ), IW, JW, DESCW,
     $                    1 )
*
            S = ZERO
            IF( MYCOL.EQ.ICURCOL ) THEN
               IF( NP.GT.0 )THEN
                  DO 130 II = IIW-1, IIW+NP-2
                     IF( RWORK( IPB+II ).GT.SAFE2 ) THEN
                        S = MAX( S, CABS1( WORK( IPR+II ) ) /
     $                              RWORK( IPB+II ) )
                     ELSE
                        S = MAX( S, ( CABS1( WORK( IPR+II ) )+SAFE1 ) /
     $                              ( RWORK( IPB+II )+SAFE1 ) )
                     END IF
  130             CONTINUE
               END IF
            END IF
*
            CALL SGAMX2D( ICTXT, 'All', ' ', 1, 1, S, 1, IDUM, IDUM, 1,
     $                    -1, MYCOL )
            IF( MYCOL.EQ.ICURCOL )
     $         BERR( JJFBE ) = S
*
*           Test stopping criterion. Continue iterating if
*             1) The residual BERR(J+K) is larger than machine epsilon,
*                and
*             2) BERR(J+K) decreased by at least a factor of 2 during
*                the last iteration, and
*             3) At most ITMAX iterations tried.
*
            IF( S.GT.EPS .AND. TWO*S.LE.LSTRES .AND.
     $          COUNT.LE.ITMAX ) THEN
*
*              Update solution and try again.
*
               CALL PCPOTRS( UPLO, N, 1, AF, IAF, JAF, DESCAF,
     $                       WORK( IPR ), IW, JW, DESCW, INFO )
               CALL PCAXPY( N, ONE, WORK( IPR ), IW, JW, DESCW, 1, X,
     $                      IX, J+K, DESCX, 1 )
               LSTRES = S
               COUNT = COUNT + 1
               GO TO 110
            END IF
*
*           Bound error from formula
*
*           norm(sub(X) - XTRUE) / norm(sub(X)) .le. FERR =
*           norm( abs(inv(sub(A)))*
*               ( abs(R) + NZ*EPS*(
*                 abs(sub(A))*abs(sub(X))+abs(sub(B)) )))/norm(sub(X))
*
*           where
*             norm(Z) is the magnitude of the largest component of Z
*             inv(sub(A)) is the inverse of sub(A)
*             abs(Z) is the componentwise absolute value of the matrix
*                or vector Z
*             NZ is the maximum number of nonzeros in any row of sub(A),
*                plus 1
*             EPS is machine epsilon
*
*           The i-th component of abs(R)+NZ*EPS*(abs(sub(A))*abs(sub(X))
*           +abs(sub(B))) is incremented by SAFE1 if the i-th component
*           of abs(sub(A))*abs(sub(X)) + abs(sub(B)) is less than SAFE2.
*
*           Use PCLACON to estimate the infinity-norm of the matrix
*           inv(sub(A)) * diag(W), where
*           W = abs(R) + NZ*EPS*( abs(sub(A))*abs(sub(X))+abs(sub(B)))))
*
            IF( MYCOL.EQ.ICURCOL ) THEN
               IF( NP.GT.0 ) THEN
                  DO 140 II = IIW-1, IIW+NP-2
                     IF( RWORK( IPB+II ).GT.SAFE2 ) THEN
                        RWORK( IPB+II ) = CABS1( WORK( IPR+II ) ) +
     $                                    NZ*EPS*RWORK( IPB+II )
                     ELSE
                        RWORK( IPB+II ) = CABS1( WORK( IPR+II ) ) +
     $                                    NZ*EPS*RWORK( IPB+II ) + SAFE1
                    END IF
  140             CONTINUE
               END IF
            END IF
*
            KASE = 0
  150       CONTINUE
            IF( MYCOL.EQ.ICURCOL ) THEN
               CALL CGEBS2D( ICTXT, 'Rowwise', ' ', NP, 1, WORK( IPR ),
     $                       DESCW( LLD_ ) )
            ELSE
               CALL CGEBR2D( ICTXT, 'Rowwise', ' ', NP, 1, WORK( IPR ),
     $                       DESCW( LLD_ ), MYROW, ICURCOL )
            END IF
            DESCW( CSRC_ ) = MYCOL
            CALL PCLACON( N, WORK( IPV ), IW, JW, DESCW, WORK( IPR ),
     $                    IW, JW, DESCW, EST, KASE )
            DESCW( CSRC_ ) = ICURCOL
*
            IF( KASE.NE.0 ) THEN
               IF( KASE.EQ.1 ) THEN
*
*                 Multiply by diag(W)*inv(sub(A)').
*
                  CALL PCPOTRS( UPLO, N, 1, AF, IAF, JAF, DESCAF,
     $                         WORK( IPR ), IW, JW, DESCW, INFO )
*
                  IF( MYCOL.EQ.ICURCOL ) THEN
                     IF( NP.GT.0 ) THEN
                        DO 160 II = IIW-1, IIW+NP-2
                           WORK( IPR+II ) = RWORK( IPB+II )*
     $                                      WORK( IPR+II )
  160                   CONTINUE
                     END IF
                  END IF
               ELSE
*
*                 Multiply by inv(sub(A))*diag(W).
*
                  IF( MYCOL.EQ.ICURCOL ) THEN
                     IF( NP.GT.0 ) THEN
                        DO 170 II = IIW-1, IIW+NP-2
                           WORK( IPR+II ) = RWORK( IPB+II )*
     $                                      WORK( IPR+II )
  170                   CONTINUE
                     END IF
                  END IF
*
                  CALL PCPOTRS( UPLO, N, 1, AF, IAF, JAF, DESCAF,
     $                          WORK( IPR ), IW, JW, DESCW, INFO )
               END IF
               GO TO 150
            END IF
*
*           Normalize error.
*
            LSTRES = ZERO
            IF( MYCOL.EQ.ICURCOL ) THEN
               IF( NP.GT.0 ) THEN
                  DO 180 II = IIXB, IIXB+NP-1
                     LSTRES = MAX( LSTRES, CABS1( X( IOFFXB+II ) ) )
  180             CONTINUE
               END IF
               CALL SGAMX2D( ICTXT, 'Column', ' ', 1, 1, LSTRES, 1,
     $                       IDUM, IDUM, 1, -1, MYCOL )
               IF( LSTRES.NE.ZERO )
     $            FERR( JJFBE ) = EST / LSTRES
*
               JJXB = JJXB + 1
               JJFBE = JJFBE + 1
               IOFFXB = IOFFXB + LDXB
*
            END IF
*
  190    CONTINUE
*
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  200 CONTINUE
*
      WORK( 1 ) = CMPLX( REAL( LWMIN ) )
      RWORK( 1 ) = REAL( LRWMIN )
*
      RETURN
*
*     End of PCPORFS
*
      END
