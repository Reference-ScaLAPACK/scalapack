      SUBROUTINE PSLAED2( ICTXT, K, N, N1, NB, D, DROW, DCOL, Q, LDQ,
     $                    RHO, Z, W, DLAMDA, Q2, LDQ2, QBUF, CTOT, PSM,
     $                    NPCOL, INDX, INDXC, INDXP, INDCOL, COLTYP, NN,
     $                    NN1, NN2, IB1, IB2 )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      INTEGER            DCOL, DROW, IB1, IB2, ICTXT, K, LDQ, LDQ2, N,
     $                   N1, NB, NN, NN1, NN2, NPCOL
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            COLTYP( * ), CTOT( 0: NPCOL-1, 4 ),
     $                   INDCOL( N ), INDX( * ), INDXC( * ), INDXP( * ),
     $                   PSM( 0: NPCOL-1, 4 )
      REAL               D( * ), DLAMDA( * ), Q( LDQ, * ),
     $                   Q2( LDQ2, * ), QBUF( * ), W( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLAED2 sorts the two sets of eigenvalues together into a single
*  sorted set.  Then it tries to deflate the size of the problem.
*  There are two ways in which deflation can occur:  when two or more
*  eigenvalues are close together or if there is a tiny entry in the
*  Z vector.  For each such occurrence the order of the related secular
*  equation problem is reduced by one.
*
*  Arguments
*  =========
*
*  ICTXT  (global input) INTEGER
*         The BLACS context handle, indicating the global context of
*         the operation on the matrix. The context itself is global.
*
*  K      (output) INTEGER
*         The number of non-deflated eigenvalues, and the order of the
*         related secular equation. 0 <= K <=N.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  N1     (input) INTEGER
*         The location of the last eigenvalue in the leading sub-matrix.
*         min(1,N) < N1 < N.
*
*  NB      (global input) INTEGER
*          The blocking factor used to distribute the columns of the
*          matrix. NB >= 1.
*
*  D      (input/output) REAL array, dimension (N)
*         On entry, D contains the eigenvalues of the two submatrices to
*         be combined.
*         On exit, D contains the trailing (N-K) updated eigenvalues
*         (those which were deflated) sorted into increasing order.
*
*  DROW   (global input) INTEGER
*          The process row over which the first row of the matrix D is
*          distributed. 0 <= DROW < NPROW.
*
*  DCOL   (global input) INTEGER
*          The process column over which the first column of the
*          matrix D is distributed. 0 <= DCOL < NPCOL.
*
*  Q      (input/output) REAL array, dimension (LDQ, N)
*         On entry, Q contains the eigenvectors of two submatrices in
*         the two square blocks with corners at (1,1), (N1,N1)
*         and (N1+1, N1+1), (N,N).
*         On exit, Q contains the trailing (N-K) updated eigenvectors
*         (those which were deflated) in its last N-K columns.
*
*  LDQ    (input) INTEGER
*         The leading dimension of the array Q.  LDQ >= max(1,NQ).
*
*  RHO    (global input/output) REAL
*         On entry, the off-diagonal element associated with the rank-1
*         cut which originally split the two submatrices which are now
*         being recombined.
*         On exit, RHO has been modified to the value required by
*         PSLAED3.
*
*  Z      (global input) REAL array, dimension (N)
*         On entry, Z contains the updating vector (the last
*         row of the first sub-eigenvector matrix and the first row of
*         the second sub-eigenvector matrix).
*         On exit, the contents of Z have been destroyed by the updating
*         process.
*
*  DLAMDA (global output) REAL array, dimension (N)
*         A copy of the first K eigenvalues which will be used by
*         SLAED3 to form the secular equation.
*
*  W      (global output) REAL array, dimension (N)
*         The first k values of the final deflation-altered z-vector
*         which will be passed to SLAED3.
*
*  Q2     (output) REAL array, dimension (LDQ2, NQ)
*         A copy of the first K eigenvectors which will be used by
*
*  LDQ2    (input) INTEGER
*         The leading dimension of the array Q2.
*
*  QBUF   (workspace) REAL array, dimension 3*N
*
*  CTOT   (workspace) INTEGER array, dimension( NPCOL, 4)
*
*  PSM    (workspace) INTEGER array, dimension( NPCOL, 4)
*
*  NPCOL   (global input) INTEGER
*          The total number of columns over which the distributed
*           submatrix is distributed.
*
*  INDX   (workspace) INTEGER array, dimension (N)
*         The permutation used to sort the contents of DLAMDA into
*         ascending order.
*
*  INDXC  (output) INTEGER array, dimension (N)
*         The permutation used to arrange the columns of the deflated
*         Q matrix into three groups:  the first group contains non-zero
*         elements only at and above N1, the second contains
*         non-zero elements only below N1, and the third is dense.
*
*  INDXP  (workspace) INTEGER array, dimension (N)
*         The permutation used to place deflated values of D at the end
*         of the array.  INDXP(1:K) points to the nondeflated D-values
*         and INDXP(K+1:N) points to the deflated eigenvalues.
*
*  INDCOL (workspace) INTEGER array, dimension (N)
*
*  COLTYP (workspace/output) INTEGER array, dimension (N)
*         During execution, a label which will indicate which of the
*         following types a column in the Q2 matrix is:
*         1 : non-zero in the upper half only;
*         2 : dense;
*         3 : non-zero in the lower half only;
*         4 : deflated.
*
*  NN     (global output) INTEGER, the order of matrix U, (PSLAED1).
*  NN1    (global output) INTEGER, the order of matrix Q1, (PSLAED1).
*  NN2    (global output) INTEGER, the order of matrix Q2, (PSLAED1).
*  IB1    (global output) INTEGER, pointeur on Q1, (PSLAED1).
*  IB2    (global output) INTEGER, pointeur on Q2, (PSLAED1).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0E0, ZERO = 0.0E0, ONE = 1.0E0,
     $                   TWO = 2.0E0, EIGHT = 8.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            COL, CT, I, IAM, IE1, IE2, IMAX, INFO, J, JJQ2,
     $                   JJS, JMAX, JS, K2, MYCOL, MYROW, N1P1, N2, NJ,
     $                   NJCOL, NJJ, NP, NPROCS, NPROW, PJ, PJCOL, PJJ
      REAL               C, EPS, S, T, TAU, TOL
*     ..
*     .. External Functions ..
      INTEGER            INDXG2L, INDXL2G, ISAMAX, NUMROC
      REAL               PSLAMCH, SLAPY2
      EXTERNAL           INDXG2L, INDXL2G, ISAMAX, NUMROC, PSLAMCH,
     $                   SLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, BLACS_PINFO, INFOG1L, SCOPY,
     $                   SGERV2D, SGESD2D, SLAPST, SROT, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD, SQRT
*     ..
*     .. External Functions ..
*     ..
*     .. Local Arrays ..
      INTEGER            PTT( 4 )
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      CALL BLACS_PINFO( IAM, NPROCS )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NP = NUMROC( N, NB, MYROW, DROW, NPROW )
*
      N2 = N - N1
      N1P1 = N1 + 1
*
      IF( RHO.LT.ZERO ) THEN
         CALL SSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
*
*     Normalize z so that norm(z) = 1.  Since z is the concatenation of
*     two normalized vectors, norm2(z) = sqrt(2).
*
      T = ONE / SQRT( TWO )
      CALL SSCAL( N, T, Z, 1 )
*
*     RHO = ABS( norm(z)**2 * RHO )
*
      RHO = ABS( TWO*RHO )
*
*     Calculate the allowable deflation tolerance
*
      IMAX = ISAMAX( N, Z, 1 )
      JMAX = ISAMAX( N, D, 1 )
      EPS = PSLAMCH( ICTXT, 'Epsilon' )
      TOL = EIGHT*EPS*MAX( ABS( D( JMAX ) ), ABS( Z( IMAX ) ) )
*
*     If the rank-1 modifier is small enough, no more needs to be done
*     except to reorganize Q so that its columns correspond with the
*     elements in D.
*
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         GO TO 220
      END IF
*
*     If there are multiple eigenvalues then the problem deflates.  Here
*     the number of equal eigenvalues are found.  As each equal
*     eigenvalue is found, an elementary reflector is computed to rotate
*     the corresponding eigensubspace so that the corresponding
*     components of Z are zero in this new basis.
*
*
      CALL SLAPST( 'I', N, D, INDX, INFO )
*
      DO 10 I = 1, N1
         COLTYP( I ) = 1
   10 CONTINUE
      DO 20 I = N1P1, N
         COLTYP( I ) = 3
   20 CONTINUE
      COL = DCOL
      DO 40 I = 1, N, NB
         DO 30 J = 0, NB - 1
            IF( I+J.LE.N )
     $         INDCOL( I+J ) = COL
   30    CONTINUE
         COL = MOD( COL+1, NPCOL )
   40 CONTINUE
*
      K = 0
      K2 = N + 1
      DO 50 J = 1, N
         NJ = INDX( J )
         IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
*
*           Deflate due to small z component.
*
            K2 = K2 - 1
            COLTYP( NJ ) = 4
            INDXP( K2 ) = NJ
            IF( J.EQ.N )
     $         GO TO 80
         ELSE
            PJ = NJ
            GO TO 60
         END IF
   50 CONTINUE
   60 CONTINUE
      J = J + 1
      NJ = INDX( J )
      IF( J.GT.N )
     $   GO TO 80
      IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
*
*        Deflate due to small z component.
*
         K2 = K2 - 1
         COLTYP( NJ ) = 4
         INDXP( K2 ) = NJ
      ELSE
*
*        Check if eigenvalues are close enough to allow deflation.
*
         S = Z( PJ )
         C = Z( NJ )
*
*        Find sqrt(a**2+b**2) without overflow or
*        destructive underflow.
*
         TAU = SLAPY2( C, S )
         T = D( NJ ) - D( PJ )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
*
*           Deflation is possible.
*
            Z( NJ ) = TAU
            Z( PJ ) = ZERO
            IF( COLTYP( NJ ).NE.COLTYP( PJ ) )
     $         COLTYP( NJ ) = 2
            COLTYP( PJ ) = 4
            CALL INFOG1L( NJ, NB, NPCOL, MYCOL, DCOL, NJJ, NJCOL )
            CALL INFOG1L( PJ, NB, NPCOL, MYCOL, DCOL, PJJ, PJCOL )
            IF( INDCOL( PJ ).EQ.INDCOL( NJ ) .AND. MYCOL.EQ.NJCOL ) THEN
               CALL SROT( NP, Q( 1, PJJ ), 1, Q( 1, NJJ ), 1, C, S )
            ELSE IF( MYCOL.EQ.PJCOL ) THEN
               CALL SGESD2D( ICTXT, NP, 1, Q( 1, PJJ ), NP, MYROW,
     $                       NJCOL )
               CALL SGERV2D( ICTXT, NP, 1, QBUF, NP, MYROW, NJCOL )
               CALL SROT( NP, Q( 1, PJJ ), 1, QBUF, 1, C, S )
            ELSE IF( MYCOL.EQ.NJCOL ) THEN
               CALL SGESD2D( ICTXT, NP, 1, Q( 1, NJJ ), NP, MYROW,
     $                       PJCOL )
               CALL SGERV2D( ICTXT, NP, 1, QBUF, NP, MYROW, PJCOL )
               CALL SROT( NP, QBUF, 1, Q( 1, NJJ ), 1, C, S )
            END IF
            T = D( PJ )*C**2 + D( NJ )*S**2
            D( NJ ) = D( PJ )*S**2 + D( NJ )*C**2
            D( PJ ) = T
            K2 = K2 - 1
            I = 1
   70       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( PJ ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = PJ
                  I = I + 1
                  GO TO 70
               ELSE
                  INDXP( K2+I-1 ) = PJ
               END IF
            ELSE
               INDXP( K2+I-1 ) = PJ
            END IF
            PJ = NJ
         ELSE
            K = K + 1
            DLAMDA( K ) = D( PJ )
            W( K ) = Z( PJ )
            INDXP( K ) = PJ
            PJ = NJ
         END IF
      END IF
      GO TO 60
   80 CONTINUE
*
*     Record the last eigenvalue.
*
      K = K + 1
      DLAMDA( K ) = D( PJ )
      W( K ) = Z( PJ )
      INDXP( K ) = PJ
*
*     Count up the total number of the various types of columns, then
*     form a permutation which positions the four column types into
*     four uniform groups (although one or more of these groups may be
*     empty).
*
      DO 100 J = 1, 4
         DO 90 I = 0, NPCOL - 1
            CTOT( I, J ) = 0
   90    CONTINUE
         PTT( J ) = 0
  100 CONTINUE
      DO 110 J = 1, N
         CT = COLTYP( J )
         COL = INDCOL( J )
         CTOT( COL, CT ) = CTOT( COL, CT ) + 1
  110 CONTINUE
*
*     PSM(*) = Position in SubMatrix (of types 1 through 4)
*
      DO 120 COL = 0, NPCOL - 1
         PSM( COL, 1 ) = 1
         PSM( COL, 2 ) = 1 + CTOT( COL, 1 )
         PSM( COL, 3 ) = PSM( COL, 2 ) + CTOT( COL, 2 )
         PSM( COL, 4 ) = PSM( COL, 3 ) + CTOT( COL, 3 )
  120 CONTINUE
      PTT( 1 ) = 1
      DO 140 I = 2, 4
         CT = 0
         DO 130 J = 0, NPCOL - 1
            CT = CT + CTOT( J, I-1 )
  130    CONTINUE
         PTT( I ) = PTT( I-1 ) + CT
  140 CONTINUE
*
*     Fill out the INDXC array so that the permutation which it induces
*     will place all type-1 columns first, all type-2 columns next,
*     then all type-3's, and finally all type-4's.
*
      DO 150 J = 1, N
         JS = INDXP( J )
         COL = INDCOL( JS )
         CT = COLTYP( JS )
         I = INDXL2G( PSM( COL, CT ), NB, COL, DCOL, NPCOL )
         INDX( J ) = I
         INDXC( PTT( CT ) ) = I
         PSM( COL, CT ) = PSM( COL, CT ) + 1
         PTT( CT ) = PTT( CT ) + 1
  150 CONTINUE
      DO 160 J = 1, N
         JS = INDXP( J )
         JJS = INDXG2L( JS, NB, J, J, NPCOL )
         COL = INDCOL( JS )
         IF( COL.EQ.MYCOL ) THEN
            I = INDX( J )
            JJQ2 = INDXG2L( I, NB, J, J, NPCOL )
            CALL SCOPY( NP, Q( 1, JJS ), 1, Q2( 1, JJQ2 ), 1 )
         END IF
  160 CONTINUE
*
*
*     The deflated eigenvalues and their corresponding vectors go back
*     into the last N - K slots of D and Q respectively.
*
      CALL SCOPY( N, D, 1, Z, 1 )
      DO 170 J = K + 1, N
         JS = INDXP( J )
         I = INDX( J )
         D( I ) = Z( JS )
  170 CONTINUE
*
      PTT( 1 ) = 1
      DO 190 I = 2, 4
         CT = 0
         DO 180 J = 0, NPCOL - 1
            CT = CT + CTOT( J, I-1 )
  180    CONTINUE
         PTT( I ) = PTT( I-1 ) + CT
  190 CONTINUE
*
*
      IB1 = INDXC( 1 )
      IE1 = IB1
      IB2 = INDXC( PTT( 2 ) )
      IE2 = IB2
      DO 200 I = 2, PTT( 3 ) - 1
         IB1 = MIN( IB1, INDXC( I ) )
         IE1 = MAX( IE1, INDXC( I ) )
  200 CONTINUE
      DO 210 I = PTT( 2 ), PTT( 4 ) - 1
         IB2 = MIN( IB2, INDXC( I ) )
         IE2 = MAX( IE2, INDXC( I ) )
  210 CONTINUE
      NN1 = IE1 - IB1 + 1
      NN2 = IE2 - IB2 + 1
      NN = MAX( IE1, IE2 ) - MIN( IB1, IB2 ) + 1
  220 CONTINUE
      RETURN
*
*     End of PSLAED2
*
      END
