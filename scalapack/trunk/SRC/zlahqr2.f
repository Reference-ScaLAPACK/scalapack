      SUBROUTINE ZLAHQR2( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
     $                    IHIZ, Z, LDZ, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 22, 2000
*
*     .. Scalar Arguments ..
      LOGICAL            WANTT, WANTZ
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLAHQR2 is an auxiliary routine called by ZHSEQR to update the
*    eigenvalues and Schur decomposition already computed by ZHSEQR, by
*    dealing with the Hessenberg submatrix in rows and columns ILO to IHI.
*  This version of ZLAHQR (not the standard LAPACK version) uses a
*    double-shift algorithm (like LAPACK's DLAHQR).
*  Unlike the standard LAPACK convention, this does not assume the
*    subdiagonal is real, nor does it work to preserve this quality if
*    given.
*
*  Arguments
*  =========
*
*  WANTT   (input) LOGICAL
*          = .TRUE. : the full Schur form T is required;
*          = .FALSE.: only eigenvalues are required.
*
*  WANTZ   (input) LOGICAL
*          = .TRUE. : the matrix of Schur vectors Z is required;
*          = .FALSE.: Schur vectors are not required.
*
*  N       (input) INTEGER
*          The order of the matrix H.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          It is assumed that H is already upper triangular in rows and
*          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).
*          ZLAHQR works primarily with the Hessenberg submatrix in rows
*          and columns ILO to IHI, but applies transformations to all of
*          H if WANTT is .TRUE..
*          1 <= ILO <= max(1,IHI); IHI <= N.
*
*  H       (input/output) COMPLEX*16 array, dimension (LDH,N)
*          On entry, the upper Hessenberg matrix H.
*          On exit, if WANTT is .TRUE., H is upper triangular in rows
*          and columns ILO:IHI.  If WANTT is .FALSE., the contents of H
*          are unspecified on exit.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H. LDH >= max(1,N).
*
*  W       (output) COMPLEX*16 array, dimension (N)
*          The computed eigenvalues ILO to IHI are stored in the
*          corresponding elements of W. If WANTT is .TRUE., the
*          eigenvalues are stored in the same order as on the diagonal
*          of the Schur form returned in H, with W(i) = H(i,i).
*
*  ILOZ    (input) INTEGER
*  IHIZ    (input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if WANTZ is .TRUE..
*          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
*
*  Z       (input/output) COMPLEX*16 array, dimension (LDZ,N)
*          If WANTZ is .TRUE., on entry Z must contain the current
*          matrix Z of transformations, and on exit Z has been updated;
*          transformations are applied only to the submatrix
*          Z(ILOZ:IHIZ,ILO:IHI).  If WANTZ is .FALSE., Z is not
*          referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z. LDZ >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          > 0: if INFO = i, ZLAHQR failed to compute all the
*               eigenvalues ILO to IHI in a total of 30*(IHI-ILO+1)
*               iterations; elements i+1:ihi of W contain those
*               eigenvalues which have been successfully computed.
*
*  Further Details
*  ===============
*
*  Modified by Mark R. Fahey, June, 2000
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   RZERO, RONE
      PARAMETER          ( RZERO = 0.0D+0, RONE = 1.0D+0 )
      DOUBLE PRECISION   DAT1, DAT2
      PARAMETER          ( DAT1 = 0.75D+0, DAT2 = -0.4375D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I1, I2, ITN, ITS, J, K, L, M, NH, NR, NZ
      DOUBLE PRECISION   CS, OVFL, S, SMLNUM, TST1, ULP, UNFL
      COMPLEX*16         CDUM, H00, H10, H11, H12, H21, H22, H33, H33S,
     $                   H43H34, H44, H44S, SN, SUM, T1, T2, T3, V1, V2,
     $                   V3
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 1 )
      COMPLEX*16         V( 3 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANHS
      EXTERNAL           DLAMCH, ZLANHS
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, ZCOPY, ZLANV2, ZLARFG, ZROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( ILO.EQ.IHI ) THEN
         W( ILO ) = H( ILO, ILO )
         RETURN
      END IF
*
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
*
*     Set machine-dependent constants for the stopping criterion.
*     If norm(H) <= sqrt(OVFL), overflow should not occur.
*
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = RONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( NH / ULP )
*
*     I1 and I2 are the indices of the first row and last column of H
*     to which transformations must be applied. If eigenvalues only are
*     being computed, I1 and I2 are set inside the main loop.
*
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
*
*     ITN is the total number of QR iterations allowed.
*
      ITN = 30*NH
*
*     The main loop begins here. I is the loop index and decreases from
*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
*     with the active submatrix in rows and columns L to I.
*     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
*     H(L,L-1) is negligible so that the matrix splits.
*
      I = IHI
   10 CONTINUE
      L = ILO
      IF( I.LT.ILO )
     $   GO TO 150
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 or 2 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      DO 130 ITS = 0, ITN
*
*        Look for a single small subdiagonal element.
*
         DO 20 K = I, L + 1, -1
            TST1 = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) )
            IF( TST1.EQ.RZERO )
     $         TST1 = ZLANHS( '1', I-L+1, H( L, L ), LDH, RWORK )
            IF( CABS1( H( K, K-1 ) ).LE.MAX( ULP*TST1, SMLNUM ) )
     $         GO TO 30
   20    CONTINUE
   30    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
*
*           H(L,L-1) is negligible
*
            H( L, L-1 ) = ZERO
         END IF
*
*        Exit from loop if a submatrix of order 1 or 2 has split off.
*
         IF( L.GE.I-1 )
     $      GO TO 140
*
*        Now the active submatrix is in rows and columns L to I. If
*        eigenvalues only are being computed, only the active submatrix
*        need be transformed.
*
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
*
         IF( ITS.EQ.10 .OR. ITS.EQ.20 ) THEN
*
*           Exceptional shift.
*
*            S = ABS( DBLE( H( I,I-1 ) ) ) + ABS( DBLE( H( I-1,I-2 ) ) )
            S = CABS1( H( I, I-1 ) ) + CABS1( H( I-1, I-2 ) )
            H44 = DAT1*S
            H33 = H44
            H43H34 = DAT2*S*S
         ELSE
*
*           Prepare to use Wilkinson's shift.
*
            H44 = H( I, I )
            H33 = H( I-1, I-1 )
            H43H34 = H( I, I-1 )*H( I-1, I )
         END IF
*
*        Look for two consecutive small subdiagonal elements.
*
         DO 40 M = I - 2, L, -1
*
*           Determine the effect of starting the double-shift QR
*           iteration at row M, and see if this would make H(M,M-1)
*           negligible.
*
            H11 = H( M, M )
            H22 = H( M+1, M+1 )
            H21 = H( M+1, M )
            H12 = H( M, M+1 )
            H44S = H44 - H11
            H33S = H33 - H11
            V1 = ( H33S*H44S-H43H34 ) / H21 + H12
            V2 = H22 - H11 - H33S - H44S
            V3 = H( M+2, M+1 )
            S = CABS1( V1 ) + CABS1( V2 ) + ABS( V3 )
            V1 = V1 / S
            V2 = V2 / S
            V3 = V3 / S
            V( 1 ) = V1
            V( 2 ) = V2
            V( 3 ) = V3
            IF( M.EQ.L )
     $         GO TO 50
            H00 = H( M-1, M-1 )
            H10 = H( M, M-1 )
            TST1 = CABS1( V1 )*( CABS1( H00 )+CABS1( H11 )+
     $             CABS1( H22 ) )
            IF( CABS1( H10 )*( CABS1( V2 )+CABS1( V3 ) ).LE.ULP*TST1 )
     $         GO TO 50
   40    CONTINUE
   50    CONTINUE
*
*        Double-shift QR step
*
         DO 120 K = M, I - 1
*
*           The first iteration of this loop determines a reflection G
*           from the vector V and applies it from left and right to H,
*           thus creating a nonzero bulge below the subdiagonal.
*
*           Each subsequent iteration determines a reflection G to
*           restore the Hessenberg form in the (K-1)th column, and thus
*           chases the bulge one step toward the bottom of the active
*           submatrix.  NR is the order of G
*
            NR = MIN( 3, I-K+1 )
            IF( K.GT.M )
     $         CALL ZCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL ZLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
               IF( K.LT.I-1 )
     $            H( K+2, K-1 ) = ZERO
            ELSE IF( M.GT.L ) THEN
*              The real double-shift code uses H( K, K-1 ) = -H( K, K-1 )
*              instead of the following.
               H( K, K-1 ) = H( K, K-1 ) - DCONJG( T1 )*H( K, K-1 )
            END IF
            V2 = V( 2 )
            T2 = T1*V2
            IF( NR.EQ.3 ) THEN
               V3 = V( 3 )
               T3 = T1*V3
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
               DO 60 J = K, I2
                  SUM = DCONJG( T1 )*H( K, J ) +
     $                  DCONJG( T2 )*H( K+1, J ) +
     $                  DCONJG( T3 )*H( K+2, J )
                  H( K, J ) = H( K, J ) - SUM
                  H( K+1, J ) = H( K+1, J ) - SUM*V2
                  H( K+2, J ) = H( K+2, J ) - SUM*V3
   60          CONTINUE
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to min(K+3,I).
*
               DO 70 J = I1, MIN( K+3, I )
                  SUM = T1*H( J, K ) + T2*H( J, K+1 ) + T3*H( J, K+2 )
                  H( J, K ) = H( J, K ) - SUM
                  H( J, K+1 ) = H( J, K+1 ) - SUM*DCONJG( V2 )
                  H( J, K+2 ) = H( J, K+2 ) - SUM*DCONJG( V3 )
   70          CONTINUE
*
               IF( WANTZ ) THEN
*
*              Accumulate transformations in the matrix Z
*
                  DO 80 J = ILOZ, IHIZ
                     SUM = T1*Z( J, K ) + T2*Z( J, K+1 ) +
     $                     T3*Z( J, K+2 )
                     Z( J, K ) = Z( J, K ) - SUM
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*DCONJG( V2 )
                     Z( J, K+2 ) = Z( J, K+2 ) - SUM*DCONJG( V3 )
   80             CONTINUE
               END IF
            ELSE IF( NR.EQ.2 ) THEN
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
               DO 90 J = K, I2
                  SUM = DCONJG( T1 )*H( K, J ) +
     $                  DCONJG( T2 )*H( K+1, J )
                  H( K, J ) = H( K, J ) - SUM
                  H( K+1, J ) = H( K+1, J ) - SUM*V2
   90          CONTINUE
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to min(K+2,I).
*
               DO 100 J = I1, MIN( K+2, I )
                  SUM = T1*H( J, K ) + T2*H( J, K+1 )
                  H( J, K ) = H( J, K ) - SUM
                  H( J, K+1 ) = H( J, K+1 ) - SUM*DCONJG( V2 )
  100          CONTINUE
*
               IF( WANTZ ) THEN
*
*                 Accumulate transformations in the matrix Z
*
                  DO 110 J = ILOZ, IHIZ
                     SUM = T1*Z( J, K ) + T2*Z( J, K+1 )
                     Z( J, K ) = Z( J, K ) - SUM
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*DCONJG( V2 )
  110             CONTINUE
               END IF
            END IF
*
*           Since at the start of the QR step we have for M > L
*              H( K, K-1 ) = H( K, K-1 ) - DCONJG( T1 )*H( K, K-1 )
*           then we don't need to do the following
*           IF( K.EQ.M .AND. M.GT.L ) THEN
*              If the QR step was started at row M > L because two
*              consecutive small subdiagonals were found, then H(M,M-1)
*              must also be updated by a factor of (1-T1).
*              TEMP = ONE - T1
*              H( m, m-1 ) = H( m, m-1 )*DCONJG( TEMP )
*           END IF
  120    CONTINUE
*
*        Ensure that H(I,I-1) is real.
*
  130 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      INFO = I
      RETURN
*
  140 CONTINUE
*
      IF( L.EQ.I ) THEN
*
*        H(I,I-1) is negligible: one eigenvalue has converged.
*
         W( I ) = H( I, I )
*
      ELSE IF( L.EQ.I-1 ) THEN
*
*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
*
*        Transform the 2-by-2 submatrix to standard Schur form,
*        and compute and store the eigenvalues.
*
         CALL ZLANV2( H( I-1, I-1 ), H( I-1, I ), H( I, I-1 ),
     $                H( I, I ), W( I-1 ), W( I ), CS, SN )
*
         IF( WANTT ) THEN
*
*           Apply the transformation to the rest of H.
*
            IF( I2.GT.I )
     $         CALL ZROT( I2-I, H( I-1, I+1 ), LDH, H( I, I+1 ), LDH,
     $                    CS, SN )
            CALL ZROT( I-I1-1, H( I1, I-1 ), 1, H( I1, I ), 1, CS,
     $                 DCONJG( SN ) )
         END IF
         IF( WANTZ ) THEN
*
*           Apply the transformation to Z.
*
            CALL ZROT( NZ, Z( ILOZ, I-1 ), 1, Z( ILOZ, I ), 1, CS,
     $                 DCONJG( SN ) )
         END IF
*
      END IF
*
*     Decrement number of remaining iterations, and return to start of
*     the main loop with new value of I.
*
      ITN = ITN - ITS
      I = L - 1
      GO TO 10
*
  150 CONTINUE
      RETURN
*
*     End of ZLAHQR2
*
      END
