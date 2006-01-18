      SUBROUTINE CSTEQR2( COMPZ, N, D, E, Z, LDZ, NR, WORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     November 15, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N, NR
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * ), WORK( * )
      COMPLEX            Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CSTEQR2 is a modified version of LAPACK routine CSTEQR.
*  CSTEQR2 computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the implicit QL or QR method.
*  CSTEQR2 is modified from CSTEQR to allow each ScaLAPACK process
*  running CSTEQR2 to perform updates on a distributed matrix Q.
*  Proper usage of CSTEQR2 can be gleaned from
*  examination of ScaLAPACK's *  PCHEEV.
*  CSTEQR2 incorporates changes attributed to Greg Henry.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'I':  Compute eigenvalues and eigenvectors of the
*                  tridiagonal matrix.  Z must be initialized to the
*                  identity matrix by PCLASET or CLASET prior
*                  to entering this subroutine.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) REAL array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  Z       (local input/local output) COMPLEX array, global
*          dimension (N, N), local dimension (LDZ, NR).
*          On entry, if  COMPZ = 'V', then Z contains the orthogonal
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original symmetric matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          eigenvectors are desired, then  LDZ >= max(1,N).
*
*  NR      (input) INTEGER
*          NR = MAX(1, NUMROC( N, NB, MYPROW, 0, NPROCS ) ).
*          If COMPZ = 'N', then NR is not referenced.
*
*  WORK    (workspace) REAL array, dimension (max(1,2*N-2))
*          If COMPZ = 'N', then WORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm has failed to find all the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero; on exit, D
*                and E contain the elements of a symmetric tridiagonal
*                matrix which is orthogonally similar to the original
*                matrix.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO, THREE, HALF
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0,
     $                   THREE = 3.0E0, HALF = 0.5E0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E0, 1.0E0 ) )
      INTEGER            MAXIT, NMAXLOOK
      PARAMETER          ( MAXIT = 30, NMAXLOOK = 15 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ILAST, ISCALE, J, JTOT, K, L,
     $                   L1, LEND, LENDM1, LENDP1, LENDSV, LM1, LSV, M,
     $                   MM, MM1, NLOOK, NM1, NMAXIT
      REAL               ANORM, B, C, EPS, EPS2, F, G, GP, OLDEL, OLDGP,
     $                   OLDRP, P, R, RP, RT1, RT2, S, SAFMAX, SAFMIN,
     $                   SSFMAX, SSFMIN, TST, TST1
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANST, SLAPY2
      EXTERNAL           LSAME, SLAMCH, SLANST, SLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLASR, CSWAP, SLAEV2, SLARTG, SLASCL, SSTERF,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      ILAST = 0
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSEIF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 1
      ELSE
         ICOMPZ = -1
      ENDIF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSEIF( N.LT.0 ) THEN
         INFO = -2
      ELSEIF( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, NR ) ) THEN
         INFO = -6
      ENDIF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSTEQR2', -INFO )
         RETURN
      ENDIF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     If eigenvectors aren't not desired, this is faster
*
      IF( ICOMPZ.EQ.0 ) THEN
         CALL SSTERF( N, D, E, INFO )
         RETURN
      ENDIF
*
      IF( N.EQ.1 ) THEN
         Z( 1, 1 ) = CONE
         RETURN
      ENDIF
*
*     Determine the unit roundoff and over/underflow thresholds.
*
      EPS = SLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = SLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues and eigenvectors of the tridiagonal
*     matrix.
*
      NMAXIT = N*MAXIT
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GOTO 220
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GOTO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GOTO 30
            ENDIF
   20    CONTINUE
      ENDIF
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GOTO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = SLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GOTO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSEIF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      ENDIF
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      ENDIF
*
      IF( LEND.GT.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GOTO 60
   50       CONTINUE
         ENDIF
*
         M = LEND
*
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GOTO 110
*
*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L+1 ) THEN
            CALL SLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
            WORK( L ) = C
            WORK( N-1+L ) = S
            CALL CLASR( 'R', 'V', 'B', NR, 2, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GOTO 40
            GOTO 200
         ENDIF
*
         IF( JTOT.EQ.NMAXIT )
     $      GOTO 200
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = SLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
*
         IF( ICOMPZ.EQ.0 ) THEN
*           Do not do a lookahead!
            GOTO 90
         ENDIF
*
         OLDEL = ABS( E( L ) )
         GP = G
         RP = R
         TST = ABS( E( L ) )**2
         TST = TST / ( ( EPS2*ABS( D( L ) ) )*ABS( D( L+1 ) )+SAFMIN )
*
         NLOOK = 1
         IF( ( TST.GT.ONE ) .AND. ( NLOOK.LE.NMAXLOOK ) ) THEN
   70       CONTINUE
*
*           This is the lookahead loop, going until we have
*           convergence or too many steps have been taken.
*
            S = ONE
            C = ONE
            P = ZERO
            MM1 = M - 1
            DO 80 I = MM1, L, -1
               F = S*E( I )
               B = C*E( I )
               CALL SLARTG( GP, F, C, S, RP )
               GP = D( I+1 ) - P
               RP = ( D( I )-GP )*S + TWO*C*B
               P = S*RP
               IF( I.NE.L )
     $            GP = C*RP - B
   80       CONTINUE
            OLDGP = GP
            OLDRP = RP
*           Find GP & RP for the next iteration
            IF( ABS( C*OLDRP-B ).GT.SAFMIN ) THEN
               GP = ( ( OLDGP+P )-( D( L )-P ) ) / ( TWO*( C*OLDRP-B ) )
            ELSE
*
*           Goto put in by G. Henry to fix ALPHA problem
*
               GOTO 90
*              GP = ( ( OLDGP+P )-( D( L )-P ) ) /
*    $              ( TWO*( C*OLDRP-B )+SAFMIN )
            ENDIF
            RP = SLAPY2( GP, ONE )
            GP = D( M ) - ( D( L )-P ) +
     $           ( ( C*OLDRP-B ) / ( GP+SIGN( RP, GP ) ) )
            TST1 = TST
            TST = ABS( C*OLDRP-B )**2
            TST = TST / ( ( EPS2*ABS( D( L )-P ) )*ABS( OLDGP+P )+
     $            SAFMIN )
*           Make sure that we are making progress
            IF( ABS( C*OLDRP-B ).GT.0.9E0*OLDEL ) THEN
               IF( ABS( C*OLDRP-B ).GT.OLDEL ) THEN
                  GP = G
                  RP = R
               ENDIF
               TST = HALF
            ELSE
               OLDEL = ABS( C*OLDRP-B )
            ENDIF
            NLOOK = NLOOK + 1
            IF( ( TST.GT.ONE ) .AND. ( NLOOK.LE.NMAXLOOK ) )
     $         GOTO 70
         ENDIF
*
         IF( ( TST.LE.ONE ) .AND. ( TST.NE.HALF ) .AND.
     $       ( ABS( P ).LT.EPS*ABS( D( L ) ) ) .AND.
     $       ( ILAST.EQ.L ) .AND. ( ABS( E( L ) )**2.LE.10000.0E0*
     $       ( ( EPS2*ABS( D( L ) ) )*ABS( D( L+1 ) )+SAFMIN ) ) ) THEN
*
*           Skip the current step: the subdiagonal info is just noise.
*
            M = L
            E( M ) = ZERO
            P = D( L )
            JTOT = JTOT - 1
            GOTO 110
         ENDIF
         G = GP
         R = RP
*
*        Lookahead over
*
   90    CONTINUE
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         MM1 = M - 1
         DO 100 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL SLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            WORK( I ) = C
            WORK( N-1+I ) = -S
*
  100    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         MM = M - L + 1
         CALL CLASR( 'R', 'V', 'B', NR, MM, WORK( L ), WORK( N-1+L ),
     $               Z( 1, L ), LDZ )
*
         D( L ) = D( L ) - P
         E( L ) = G
         ILAST = L
         GOTO 40
*
*        Eigenvalue found.
*
  110    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GOTO 40
         GOTO 200
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
  120    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 130 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GOTO 140
  130       CONTINUE
         ENDIF
*
         M = LEND
*
  140    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GOTO 190
*
*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L-1 ) THEN
            CALL SLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
            WORK( M ) = C
            WORK( N-1+M ) = S
            CALL CLASR( 'R', 'V', 'F', NR, 2, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, L-1 ), LDZ )
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GOTO 120
            GOTO 200
         ENDIF
*
         IF( JTOT.EQ.NMAXIT )
     $      GOTO 200
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = SLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
*
         IF( ICOMPZ.EQ.0 ) THEN
*           Do not do a lookahead!
            GOTO 170
         ENDIF
*
         OLDEL = ABS( E( L-1 ) )
         GP = G
         RP = R
         TST = ABS( E( L-1 ) )**2
         TST = TST / ( ( EPS2*ABS( D( L ) ) )*ABS( D( L-1 ) )+SAFMIN )
         NLOOK = 1
         IF( ( TST.GT.ONE ) .AND. ( NLOOK.LE.NMAXLOOK ) ) THEN
  150       CONTINUE
*
*           This is the lookahead loop, going until we have
*           convergence or too many steps have been taken.
*
            S = ONE
            C = ONE
            P = ZERO
*
*        Inner loop
*
            LM1 = L - 1
            DO 160 I = M, LM1
               F = S*E( I )
               B = C*E( I )
               CALL SLARTG( GP, F, C, S, RP )
               GP = D( I ) - P
               RP = ( D( I+1 )-GP )*S + TWO*C*B
               P = S*RP
               IF( I.LT.LM1 )
     $            GP = C*RP - B
  160       CONTINUE
            OLDGP = GP
            OLDRP = RP
*           Find GP & RP for the next iteration
            IF( ABS( C*OLDRP-B ).GT.SAFMIN ) THEN
               GP = ( ( OLDGP+P )-( D( L )-P ) ) / ( TWO*( C*OLDRP-B ) )
            ELSE
*
*           Goto put in by G. Henry to fix ALPHA problem
*
               GOTO 170
*              GP = ( ( OLDGP+P )-( D( L )-P ) ) /
*    $              ( TWO*( C*OLDRP-B )+SAFMIN )
            ENDIF
            RP = SLAPY2( GP, ONE )
            GP = D( M ) - ( D( L )-P ) +
     $           ( ( C*OLDRP-B ) / ( GP+SIGN( RP, GP ) ) )
            TST1 = TST
            TST = ABS( ( C*OLDRP-B ) )**2
            TST = TST / ( ( EPS2*ABS( D( L )-P ) )*ABS( OLDGP+P )+
     $            SAFMIN )
*           Make sure that we are making progress
            IF( ABS( C*OLDRP-B ).GT.0.9E0*OLDEL ) THEN
               IF( ABS( C*OLDRP-B ).GT.OLDEL ) THEN
                  GP = G
                  RP = R
               ENDIF
               TST = HALF
            ELSE
               OLDEL = ABS( C*OLDRP-B )
            ENDIF
            NLOOK = NLOOK + 1
            IF( ( TST.GT.ONE ) .AND. ( NLOOK.LE.NMAXLOOK ) )
     $         GOTO 150
         ENDIF
         IF( ( TST.LE.ONE ) .AND. ( TST.NE.HALF ) .AND.
     $       ( ABS( P ).LT.EPS*ABS( D( L ) ) ) .AND.
     $       ( ILAST.EQ.L ) .AND. ( ABS( E( L-1 ) )**2.LE.10000.0E0*
     $       ( ( EPS2*ABS( D( L-1 ) ) )*ABS( D( L ) )+SAFMIN ) ) ) THEN
*
*           Skip the current step: the subdiagonal info is just noise.
*
            M = L
            E( M-1 ) = ZERO
            P = D( L )
            JTOT = JTOT - 1
            GOTO 190
         ENDIF
*
         G = GP
         R = RP
*
*        Lookahead over
*
  170    CONTINUE
*
         S = ONE
         C = ONE
         P = ZERO
         DO 180 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL SLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            WORK( I ) = C
            WORK( N-1+I ) = S
*
  180    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         MM = L - M + 1
         CALL CLASR( 'R', 'V', 'F', NR, MM, WORK( M ), WORK( N-1+M ),
     $               Z( 1, M ), LDZ )
*
         D( L ) = D( L ) - P
         E( LM1 ) = G
         ILAST = L
         GOTO 120
*
*        Eigenvalue found.
*
  190    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GOTO 120
         GOTO 200
*
      ENDIF
*
*     Undo scaling if necessary
*
  200 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL SLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL SLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSEIF( ISCALE.EQ.2 ) THEN
         CALL SLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL SLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ENDIF
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.LT.NMAXIT )
     $   GOTO 10
      DO 210 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  210 CONTINUE
      GOTO 250
*
*     Order eigenvalues and eigenvectors.
*
  220 CONTINUE
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
      DO 240 II = 2, N
         I = II - 1
         K = I
         P = D( I )
         DO 230 J = II, N
            IF( D( J ).LT.P ) THEN
               K = J
               P = D( J )
            ENDIF
  230    CONTINUE
         IF( K.NE.I ) THEN
            D( K ) = D( I )
            D( I ) = P
            CALL CSWAP( NR, Z( 1, I ), 1, Z( 1, K ), 1 )
         ENDIF
  240 CONTINUE
*
  250 CONTINUE
*     WRITE( *, FMT = * )'JTOT', JTOT
      RETURN
*
*     End of SSTEQR2
*
      END
