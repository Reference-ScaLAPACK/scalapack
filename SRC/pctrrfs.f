      SUBROUTINE PCTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     $                    B, IB, JB, DESCB, X, IX, JX, DESCX, FERR,
     $                    BERR, WORK, LWORK, RWORK, LRWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, IA, IB, IX, JA, JB, JX, LRWORK, LWORK,
     $                   N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCX( * )
      REAL               BERR( * ), FERR( * ), RWORK( * )
      COMPLEX            A( * ), B( * ), WORK( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  PCTRRFS provides error bounds and backward error estimates for the
*  solution to a system of linear equations with a triangular
*  coefficient matrix.
*
*  The solution matrix X must be computed by PCTRTRS or some other
*  means before entering this routine.  PCTRRFS does not do iterative
*  refinement because doing so cannot improve the backward error.
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
*          = 'U':  sub( A ) is upper triangular;
*          = 'L':  sub( A ) is lower triangular.
*
*  TRANS   (global input) CHARACTER*1
*          Specifies the form of the system of equations.
*          = 'N': sub( A ) * sub( X ) = sub( B )          (No transpose)
*          = 'T': sub( A )**T * sub( X ) = sub( B )          (Transpose)
*          = 'C': sub( A )**H * sub( X ) = sub( B )
*                                                  (Conjugate transpose)
*
*  DIAG    (global input) CHARACTER*1
*          = 'N':  sub( A ) is non-unit triangular;
*          = 'U':  sub( A ) is unit triangular.
*
*  N       (global input) INTEGER
*          The order of the matrix sub( A ).  N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrices sub( B ) and sub( X ).  NRHS >= 0.
*
*  A       (local input) COMPLEX pointer into the local memory
*          to an array of local dimension (LLD_A,LOCc(JA+N-1) ). This
*          array contains the local pieces of the original triangular
*          distributed matrix sub( A ).
*          If UPLO = 'U', the leading N-by-N upper triangular part of
*          sub( A ) contains the upper triangular part of the matrix,
*          and its strictly lower triangular part is not referenced.
*          If UPLO = 'L', the leading N-by-N lower triangular part of
*          sub( A ) contains the lower triangular part of the distribu-
*          ted matrix, and its strictly upper triangular part is not
*          referenced.
*          If DIAG = 'U', the diagonal elements of sub( A ) are also
*          not referenced and are assumed to be 1.
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
*          solution vectors sub( X ).
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
*          LOCc(JB+NRHS-1). The estimated forward error bounds for
*          each solution vector of sub( X ).  If XTRUE is the true
*          solution, FERR bounds the magnitude of the largest entry
*          in (sub( X ) - XTRUE) divided by the magnitude of the
*          largest entry in sub( X ).  The estimate is as reliable as
*          the estimate for RCOND, and is almost always a slight
*          overestimate of the true error.
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
*          LWORK >= 2*LOCr( N + MOD( IA-1, MB_A ) ).
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  RWORK   (local workspace/local output) REAL array,
*                                                 dimension (LRWORK)
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
*  Notes
*  =====
*
*  This routine temporarily returns when N <= 1.
*
*  The distributed submatrices sub( X ) and sub( B ) should be
*  distributed the same way on the same processes.  These conditions
*  ensure that sub( X ) and sub( B ) are "perfectly" aligned.
*
*  Moreover, this routine requires the distributed submatrices sub( A ),
*  sub( X ), and sub( B ) to be aligned on a block boundary,
*  i.e., if f(x,y) = MOD( x-1, y ):
*  f( IA, DESCA( MB_ ) ) = f( JA, DESCA( NB_ ) ) = 0,
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
      REAL               ZERO, RONE
      PARAMETER          ( ZERO = 0.0E+0, RONE = 1.0E+0 )
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, NOTRAN, NOUNIT, UPPER
      CHARACTER          TRANSN, TRANST
      INTEGER            IAROW, IXBCOL, IXBROW, IXCOL, IXROW, ICOFFA,
     $                   ICOFFB, ICOFFX, ICTXT, ICURCOL, IDUM, II, IIXB,
     $                   IIW, IOFFXB, IPB, IPR, IPV, IROFFA, IROFFB,
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
     $                   DESCSET, INFOG2L, PCATRMV, PCAXPY, PCHK1MAT,
     $                   PCHK2MAT, PCCOPY, PCLACON, PCTRMV,
     $                   PCTRSV, PXERBLA, SGAMX2D
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
         INFO = -( 900+CTXT_ )
      ELSE
         CALL CHK1MAT( N, 4, N, 4, IA, JA, DESCA, 9, INFO )
         CALL CHK1MAT( N, 4, NRHS, 5, IB, JB, DESCB, 13, INFO )
         CALL CHK1MAT( N, 4, NRHS, 5, IX, JX, DESCX, 17, INFO )
         IF( INFO.EQ.0 ) THEN
            UPPER = LSAME( UPLO, 'U' )
            NOTRAN = LSAME( TRANS, 'N' )
            NOUNIT = LSAME( DIAG, 'N' )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            ICOFFB = MOD( JB-1, DESCB( NB_ ) )
            IROFFX = MOD( IX-1, DESCX( MB_ ) )
            ICOFFX = MOD( JX-1, DESCX( NB_ ) )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIXB, JJXB, IXBROW, IXBCOL )
            IXROW = INDXG2P( IX, DESCX( MB_ ), MYROW, DESCX( RSRC_ ),
     $                       NPROW )
            IXCOL = INDXG2P( JX, DESCX( NB_ ), MYCOL, DESCX( CSRC_ ),
     $                       NPCOL )
            NPMOD = NUMROC( N+IROFFA, DESCA( MB_ ), MYROW, IAROW,
     $                      NPROW )
            LWMIN = 2*NPMOD
            WORK( 1 ) = REAL( LWMIN )
            LRWMIN = NPMOD
            RWORK( 1 ) = REAL( LRWMIN )
            LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
*
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND.
     $         .NOT.LSAME( TRANS, 'C' ) ) THEN
               INFO = -2
            ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
               INFO = -3
            ELSE IF( N.LT.0 ) THEN
               INFO = -4
            ELSE IF( NRHS.LT.0 ) THEN
               INFO = -5
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -7
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -8
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 900+NB_ )
            ELSE IF( IROFFA.NE.IROFFB .OR. IAROW.NE.IXBROW ) THEN
               INFO = -11
            ELSE IF( DESCA( MB_ ).NE.DESCB( MB_ ) ) THEN
               INFO = -( 1300+MB_ )
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -( 1300+CTXT_ )
            ELSE IF( IROFFX.NE.0 .OR. IXBROW.NE.IXROW ) THEN
               INFO = -15
            ELSE IF( ICOFFB.NE.ICOFFX .OR. IXBCOL.NE.IXCOL ) THEN
               INFO = -16
            ELSE IF( DESCB( MB_ ).NE.DESCX( MB_ ) ) THEN
               INFO = -( 1700+MB_ )
            ELSE IF( DESCB( NB_ ).NE.DESCX( NB_ ) ) THEN
               INFO = -( 1700+NB_ )
            ELSE IF( ICTXT.NE.DESCX( CTXT_ ) ) THEN
               INFO = -( 1700+CTXT_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -21
            ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -23
            END IF
         END IF
*
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 1
         IF( NOTRAN ) THEN
            IDUM1( 2 ) = ICHAR( 'N' )
         ELSE IF( LSAME( TRANS, 'T' ) ) THEN
            IDUM1( 2 ) = ICHAR( 'T' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'C' )
         END IF
         IDUM2( 2 ) = 2
         IF( NOUNIT ) THEN
            IDUM1( 3 ) = ICHAR( 'N' )
         ELSE
            IDUM1( 3 ) = ICHAR( 'U' )
         END IF
         IDUM2( 3 ) = 3
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 4 ) = -1
         ELSE
            IDUM1( 4 ) = 1
         END IF
         IDUM2( 4 ) = 21
         IF( LRWORK.EQ.-1 ) THEN
            IDUM1( 5 ) = -1
         ELSE
            IDUM1( 5 ) = 1
         END IF
         IDUM2( 5 ) = 23
         CALL PCHK1MAT( N, 4, N, 4, IA, JA, DESCA, 9, 0, IDUM1, IDUM2,
     $                  INFO )
         CALL PCHK2MAT( N, 4, NRHS, 5, IB, JB, DESCB, 13, N, 4, NRHS, 5,
     $                  IX, JX, DESCX, 17, 5, IDUM1, IDUM2, INFO )
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCTRRFS', -INFO )
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
      IF( NOTRAN ) THEN
         TRANSN = 'N'
         TRANST = 'C'
      ELSE
         TRANSN = 'C'
         TRANST = 'N'
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
*     NZ = maximum number of nonzero entries in each row of A, plus 1
*
      NZ = N + 1
      EPS = PSLAMCH( ICTXT, 'Epsilon' )
      SAFMIN = PSLAMCH( ICTXT, 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS
      JN = MIN( ICEIL( JB, DESCB( NB_ ) )*DESCB( NB_ ), JB+NRHS-1 )
*
*     Handle first block separately
*
      JBRHS = JN - JB + 1
      DO 90 K = 0, JBRHS - 1
*
*        Compute residual R = B - op(A) * X,
*        where op(A) = A or A', depending on TRANS.
*
         CALL PCCOPY( N, X, IX, JX+K, DESCX, 1, WORK( IPR ), IW, JW,
     $                DESCW, 1 )
         CALL PCTRMV( UPLO, TRANS, DIAG, N, A, IA, JA, DESCA,
     $                WORK( IPR ), IW, JW, DESCW, 1 )
         CALL PCAXPY( N, -ONE, B, IB, JB+K, DESCB, 1, WORK( IPR ), IW,
     $                JW, DESCW, 1 )
*
*        Compute componentwise relative backward error from formula
*
*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
*
*        where abs(Z) is the componentwise absolute value of the matrix
*        or vector Z.  If the i-th component of the denominator is less
*        than SAFE2, then SAFE1 is added to the i-th components of the
*        numerator and denominator before dividing.
*
         IF( MYCOL.EQ.IXBCOL ) THEN
            IF( NP.GT.0 ) THEN
               DO 20 II = IIXB, IIXB + NP - 1
                  RWORK( IIW+II-IIXB ) = CABS1( B( II+IOFFXB ) )
   20          CONTINUE
            END IF
         END IF
*
         CALL PCATRMV( UPLO, TRANS, DIAG, N, RONE, A, IA, JA, DESCA, X,
     $                 IX, JX+K, DESCX, 1, RONE, RWORK( IPB ), IW, JW,
     $                 DESCW, 1 )
*
         S = ZERO
         IF( MYCOL.EQ.IXBCOL ) THEN
            IF( NP.GT.0 ) THEN
               DO 30 II = IIW - 1, IIW + NP - 2
                  IF( RWORK( IPB+II ).GT.SAFE2 ) THEN
                     S = MAX( S, CABS1( WORK( IPR+II ) ) /
     $                           RWORK( IPB+II ) )
                  ELSE
                     S = MAX( S, ( CABS1( WORK( IPR+II ) )+SAFE1 ) /
     $                           ( RWORK( IPB+II )+SAFE1 ) )
                  END IF
   30          CONTINUE
            END IF
         END IF
*
         CALL SGAMX2D( ICTXT, 'All', ' ', 1, 1, S, 1, IDUM, IDUM, 1,
     $                 -1, MYCOL )
         IF( MYCOL.EQ.IXBCOL )
     $      BERR( JJFBE ) = S
*
*        Bound error from formula
*
*        norm(X - XTRUE) / norm(X) .le. FERR =
*        norm( abs(inv(op(A)))*
*           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)
*
*        where
*          norm(Z) is the magnitude of the largest component of Z
*          inv(op(A)) is the inverse of op(A)
*          abs(Z) is the componentwise absolute value of the matrix or
*             vector Z
*          NZ is the maximum number of nonzeros in any row of A, plus 1
*          EPS is machine epsilon
*
*        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
*        is incremented by SAFE1 if the i-th component of
*        abs(op(A))*abs(X) + abs(B) is less than SAFE2.
*
*        Use PCLACON to estimate the infinity-norm of the matrix
*           inv(op(A)) * diag(W),
*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
*
         IF( MYCOL.EQ.IXBCOL ) THEN
            IF( NP.GT.0 ) THEN
               DO 40 II = IIW - 1, IIW + NP - 2
                  IF( RWORK( IPB+II ).GT.SAFE2 ) THEN
                     RWORK( IPB+II ) = CABS1( WORK( IPR+II ) ) +
     $                                 NZ*EPS*RWORK( IPB+II )
                  ELSE
                     RWORK( IPB+II ) = CABS1( WORK( IPR+II ) ) +
     $                                 NZ*EPS*RWORK( IPB+II ) + SAFE1
                  END IF
   40          CONTINUE
            END IF
         END IF
*
         KASE = 0
   50    CONTINUE
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
*              Multiply by diag(W)*inv(op(A)').
*
               CALL PCTRSV( UPLO, TRANST, DIAG, N, A, IA, JA, DESCA,
     $                      WORK( IPR ), IW, JW, DESCW, 1 )
               IF( MYCOL.EQ.IXBCOL ) THEN
                  IF( NP.GT.0 ) THEN
                     DO 60 II = IIW - 1, IIW + NP - 2
                        WORK( IPR+II ) = RWORK( IPB+II )*WORK( IPR+II )
   60                CONTINUE
                  END IF
               END IF
            ELSE
*
*              Multiply by inv(op(A))*diag(W).
*
               IF( MYCOL.EQ.IXBCOL ) THEN
                  IF( NP.GT.0 ) THEN
                     DO 70 II = IIW - 1, IIW + NP - 2
                        WORK( IPR+II ) = RWORK( IPB+II )*WORK( IPR+II )
   70                CONTINUE
                  END IF
               END IF
               CALL PCTRSV( UPLO, TRANSN, DIAG, N, A, IA, JA, DESCA,
     $                      WORK( IPR ), IW, JW, DESCW, 1 )
            END IF
            GO TO 50
         END IF
*
*        Normalize error.
*
         LSTRES = ZERO
         IF( MYCOL.EQ.IXBCOL ) THEN
            IF( NP.GT.0 ) THEN
               DO 80 II = IIXB, IIXB + NP - 1
                  LSTRES = MAX( LSTRES, CABS1( X( IOFFXB+II ) ) )
   80          CONTINUE
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
   90 CONTINUE
*
      ICURCOL = MOD( IXBCOL+1, NPCOL )
*
*     Do for each right hand side
*
      DO 180 J = JN + 1, JB + NRHS - 1, DESCB( NB_ )
         JBRHS = MIN( JB+NRHS-J, DESCB( NB_ ) )
         DESCW( CSRC_ ) = ICURCOL
*
         DO 170 K = 0, JBRHS - 1
*
*           Compute residual R = B - op(A) * X,
*           where op(A) = A or A', depending on TRANS.
*
            CALL PCCOPY( N, X, IX, J+K, DESCX, 1, WORK( IPR ), IW, JW,
     $                   DESCW, 1 )
            CALL PCTRMV( UPLO, TRANS, DIAG, N, A, IA, JA, DESCA,
     $                   WORK( IPR ), IW, JW, DESCW, 1 )
            CALL PCAXPY( N, -ONE, B, IB, J+K, DESCB, 1, WORK( IPR ),
     $                   IW, JW, DESCW, 1 )
*
*           Compute componentwise relative backward error from formula
*
*           max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
*
*           where abs(Z) is the componentwise absolute value of the
*           matrix or vector Z.  If the i-th component of the
*           denominator is less than SAFE2, then SAFE1 is added to the
*           i-th components of the numerator and denominator before
*           dividing.
*
            IF( MYCOL.EQ.IXBCOL ) THEN
               IF( NP.GT.0 ) THEN
                  DO 100 II = IIXB, IIXB + NP - 1
                     RWORK( IIW+II-IIXB ) = CABS1( B( II+IOFFXB ) )
  100             CONTINUE
               END IF
            END IF
*
            CALL PCATRMV( UPLO, TRANS, DIAG, N, RONE, A, IA, JA, DESCA,
     $                    X, IX, J+K, DESCX, 1, RONE, RWORK( IPB ), IW,
     $                    JW, DESCW, 1 )
*
            S = ZERO
            IF( MYCOL.EQ.IXBCOL ) THEN
               IF( NP.GT.0 ) THEN
                  DO 110 II = IIW - 1, IIW + NP - 2
                     IF( RWORK( IPB+II ).GT.SAFE2 ) THEN
                        S = MAX( S, CABS1( WORK( IPR+II ) ) /
     $                              RWORK( IPB+II ) )
                     ELSE
                        S = MAX( S, ( CABS1( WORK( IPR+II ) )+SAFE1 ) /
     $                              ( RWORK( IPB+II )+SAFE1 ) )
                     END IF
  110             CONTINUE
               END IF
            END IF
*
            CALL SGAMX2D( ICTXT, 'All', ' ', 1, 1, S, 1, IDUM, IDUM, 1,
     $                    -1, MYCOL )
            IF( MYCOL.EQ.IXBCOL )
     $         BERR( JJFBE ) = S
*
*           Bound error from formula
*
*           norm(X - XTRUE) / norm(X) .le. FERR =
*           norm( abs(inv(op(A)))*
*              ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))/norm(X)
*
*           where
*             norm(Z) is the magnitude of the largest component of Z
*             inv(op(A)) is the inverse of op(A)
*             abs(Z) is the componentwise absolute value of the matrix
*                or vector Z
*             NZ is the maximum number of nonzeros in any row of A,
*                plus 1
*             EPS is machine epsilon
*
*           The i-th component of
*              abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
*           is incremented by SAFE1 if the i-th component of
*           abs(op(A))*abs(X) + abs(B) is less than SAFE2.
*
*           Use PCLACON to estimate the infinity-norm of the matrix
*              inv(op(A)) * diag(W),
*           where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
*
            IF( MYCOL.EQ.IXBCOL ) THEN
               IF( NP.GT.0 ) THEN
                  DO 120 II = IIW - 1, IIW + NP - 2
                     IF( RWORK( IPB+II ).GT.SAFE2 ) THEN
                        RWORK( IPB+II ) = CABS1( WORK( IPR+II ) ) +
     $                                    NZ*EPS*RWORK( IPB+II )
                     ELSE
                        RWORK( IPB+II ) = CABS1( WORK( IPR+II ) ) +
     $                                    NZ*EPS*RWORK( IPB+II ) + SAFE1
                     END IF
  120             CONTINUE
               END IF
            END IF
*
            KASE = 0
  130       CONTINUE
            IF( MYCOL.EQ.IXBCOL ) THEN
               CALL CGEBS2D( ICTXT, 'Rowwise', ' ', NP, 1, WORK( IPR ),
     $                       DESCW( LLD_ ) )
            ELSE
               CALL CGEBR2D( ICTXT, 'Rowwise', ' ', NP, 1, WORK( IPR ),
     $                       DESCW( LLD_ ), MYROW, IXBCOL )
            END IF
            DESCW( CSRC_ ) = MYCOL
            CALL PCLACON( N, WORK( IPV ), IW, JW, DESCW, WORK( IPR ),
     $                    IW, JW, DESCW, EST, KASE )
            DESCW( CSRC_ ) = IXBCOL
*
            IF( KASE.NE.0 ) THEN
               IF( KASE.EQ.1 ) THEN
*
*                 Multiply by diag(W)*inv(op(A)').
*
                  CALL PCTRSV( UPLO, TRANST, DIAG, N, A, IA, JA, DESCA,
     $                         WORK( IPR ), IW, JW, DESCW, 1 )
                  IF( MYCOL.EQ.IXBCOL ) THEN
                     IF( NP.GT.0 ) THEN
                        DO 140 II = IIW - 1, IIW + NP - 2
                           WORK( IPR+II ) = RWORK( IPB+II )*
     $                                      WORK( IPR+II )
  140                   CONTINUE
                     END IF
                  END IF
               ELSE
*
*                 Multiply by inv(op(A))*diag(W).
*
                  IF( MYCOL.EQ.IXBCOL ) THEN
                     IF( NP.GT.0 ) THEN
                        DO 150 II = IIW - 1, IIW + NP - 2
                           WORK( IPR+II ) = RWORK( IPB+II )*
     $                                      WORK( IPR+II )
  150                   CONTINUE
                     END IF
                  END IF
                  CALL PCTRSV( UPLO, TRANSN, DIAG, N, A, IA, JA, DESCA,
     $                         WORK( IPR ), IW, JW, DESCW, 1 )
               END IF
               GO TO 130
            END IF
*
*           Normalize error.
*
            LSTRES = ZERO
            IF( MYCOL.EQ.IXBCOL ) THEN
               IF( NP.GT.0 ) THEN
                  DO 160 II = IIXB, IIXB + NP - 1
                     LSTRES = MAX( LSTRES, CABS1( X( IOFFXB+II ) ) )
  160             CONTINUE
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
  170    CONTINUE
*
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  180 CONTINUE
*
      WORK( 1 ) = CMPLX( REAL( LWMIN ) )
      RWORK( 1 ) = REAL( LRWMIN )
*
      RETURN
*
*     End of PCTRRFS
*
      END
