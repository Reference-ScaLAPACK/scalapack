      SUBROUTINE DLAR1VA(N, B1, BN, LAMBDA, D, L, LD, LLD, 
     $           PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ, MINGMA, 
     $           R, ISUPPZ, NRMINV, RESID, RQCORR, WORK )
*
      IMPLICIT NONE
*
*  -- ScaLAPACK computational routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ of Colorado Denver
*     July 4, 2010
*
*     .. Scalar Arguments ..
      LOGICAL            WANTNC
      INTEGER   B1, BN, N, NEGCNT, R
      DOUBLE PRECISION   GAPTOL, LAMBDA, MINGMA, NRMINV, PIVMIN, RESID,
     $                   RQCORR, ZTZ
*     ..
*     .. Array Arguments ..
      INTEGER            ISUPPZ( * )
      DOUBLE PRECISION   D( * ), L( * ), LD( * ), LLD( * ),
     $                  WORK( * )
      DOUBLE PRECISION Z( * )
*
*  Purpose
*  =======
*
*  DLAR1VA computes the (scaled) r-th column of the inverse of
*  the sumbmatrix in rows B1 through BN of the tridiagonal matrix
*  L D L^T - sigma I. When sigma is close to an eigenvalue, the
*  computed vector is an accurate eigenvector. Usually, r corresponds
*  to the index where the eigenvector is largest in magnitude.
*  The following steps accomplish this computation :
*  (a) Stationary qd transform,  L D L^T - sigma I = L(+) D(+) L(+)^T,
*  (b) Progressive qd transform, L D L^T - sigma I = U(-) D(-) U(-)^T,
*  (c) Computation of the diagonal elements of the inverse of
*      L D L^T - sigma I by combining the above transforms, and choosing
*      r as the index where the diagonal of the inverse is (one of the)
*      largest in magnitude.
*  (d) Computation of the (scaled) r-th column of the inverse using the
*      twisted factorization obtained by combining the top part of the
*      the stationary and the bottom part of the progressive transform.
*
*  Arguments
*  =========
*
*  N        (input) INTEGER
*           The order of the matrix L D L^T.
*
*  B1       (input) INTEGER
*           First index of the submatrix of L D L^T.
*
*  BN       (input) INTEGER
*           Last index of the submatrix of L D L^T.
*
*  LAMBDA    (input) DOUBLE PRECISION
*           The shift. In order to compute an accurate eigenvector,
*           LAMBDA should be a good approximation to an eigenvalue
*           of L D L^T.
*
*  L        (input) DOUBLE PRECISION array, dimension (N-1)
*           The (n-1) subdiagonal elements of the unit bidiagonal matrix
*           L, in elements 1 to N-1.
*
*  D        (input) DOUBLE PRECISION array, dimension (N)
*           The n diagonal elements of the diagonal matrix D.
*
*  LD       (input) DOUBLE PRECISION array, dimension (N-1)
*           The n-1 elements L(i)*D(i).
*
*  LLD      (input) DOUBLE PRECISION array, dimension (N-1)
*           The n-1 elements L(i)*L(i)*D(i).
*
*  PIVMIN   (input) DOUBLE PRECISION
*           The minimum pivot in the Sturm sequence.
*           
*  GAPTOL   (input) DOUBLE PRECISION
*           Tolerance that indicates when eigenvector entries are negligible
*           w.r.t. their contribution to the residual.
*
*  Z        (input/output) DOUBLE PRECISION array, dimension (N)
*           On input, all entries of Z must be set to 0.
*           On output, Z contains the (scaled) r-th column of the
*           inverse. The scaling is such that Z(R) equals 1.
*
*  WANTNC   (input) LOGICAL
*           Specifies whether NEGCNT has to be computed.
*
*  NEGCNT   (output) INTEGER
*           If WANTNC is .TRUE. then NEGCNT = the number of pivots < pivmin 
*           in the  matrix factorization L D L^T, and NEGCNT = -1 otherwise.
*
*  ZTZ      (output) DOUBLE PRECISION
*           The square of the 2-norm of Z.
*
*  MINGMA   (output) DOUBLE PRECISION
*           The reciprocal of the largest (in magnitude) diagonal
*           element of the inverse of L D L^T - sigma I.
*
*  R        (input/output) INTEGER
*           The twist index for the twisted factorization used to
*           compute Z.
*           On input, 0 <= R <= N. If R is input as 0, R is set to
*           the index where (L D L^T - sigma I)^{-1} is largest
*           in magnitude. If 1 <= R <= N, R is unchanged.
*           On output, R contains the twist index used to compute Z.
*           Ideally, R designates the position of the maximum entry in the
*           eigenvector.
*
*  ISUPPZ   (output) INTEGER array, dimension (2)
*           The support of the vector in Z, i.e., the vector Z is
*           nonzero only in elements ISUPPZ(1) through ISUPPZ( 2 ).
*
*  NRMINV   (output) DOUBLE PRECISION
*           NRMINV = 1/SQRT( ZTZ )
*
*  RESID    (output) DOUBLE PRECISION
*           The residual of the FP vector.
*           RESID = ABS( MINGMA )/SQRT( ZTZ )
*
*  RQCORR   (output) DOUBLE PRECISION
*           The Rayleigh Quotient correction to LAMBDA.
*           RQCORR = MINGMA*TMP
*
*  WORK     (workspace) DOUBLE PRECISION array, dimension (4*N)
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Beresford Parlett, University of California, Berkeley, USA
*     Jim Demmel, University of California, Berkeley, USA
*     Inderjit Dhillon, University of Texas, Austin, USA
*     Osni Marques, LBNL/NERSC, USA
*     Christof Voemel, University of California, Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLKLEN
      PARAMETER          ( BLKLEN = 16 )
       DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )

*     ..
*     .. Local Scalars ..
      LOGICAL            SAWNAN1, SAWNAN2
      INTEGER            BI, I, INDLPL, INDP, INDS, INDUMN, NB, NEG1,
     $                   NEG2, NX, R1, R2, TO
      DOUBLE PRECISION            ABSZCUR, ABSZPREV, DMINUS, DPLUS, EPS,
     $                            S, TMP, ZPREV
*     ..
*     .. External Functions ..
      LOGICAL DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DISNAN, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, DBLE
*     ..
*     .. Executable Statements ..
*
      EPS = DLAMCH( 'Precision' )

      
      IF( R.EQ.0 ) THEN
         R1 = B1
         R2 = BN
      ELSE
         R1 = R
         R2 = R
      END IF

*     Storage for LPLUS
      INDLPL = 0
*     Storage for UMINUS
      INDUMN = N
      INDS = 2*N + 1
      INDP = 3*N + 1

      IF( B1.EQ.1 ) THEN
         WORK( INDS ) = ZERO
      ELSE
         WORK( INDS+B1-1 ) = LLD( B1-1 )
      END IF

*
*     Compute the stationary transform (using the differential form)
*     until the index R2.
*
      SAWNAN1 = .FALSE.
      NEG1 = 0
      S = WORK( INDS+B1-1 ) - LAMBDA
      DO 50 I = B1, R1 - 1
         DPLUS = D( I ) + S
         WORK( INDLPL+I ) = LD( I ) / DPLUS
         IF(DPLUS.LT.ZERO) NEG1 = NEG1 + 1
         WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
         S = WORK( INDS+I ) - LAMBDA
 50   CONTINUE
      SAWNAN1 = DISNAN( S )
      IF( SAWNAN1 ) GOTO 60     
      DO 51 I = R1, R2 - 1
         DPLUS = D( I ) + S
         WORK( INDLPL+I ) = LD( I ) / DPLUS
         WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
         S = WORK( INDS+I ) - LAMBDA
 51   CONTINUE
      SAWNAN1 = DISNAN( S )
*
 60   CONTINUE
      IF( SAWNAN1 ) THEN
*        Runs a slower version of the above loop if a NaN is detected
         NEG1 = 0
         S = WORK( INDS+B1-1 ) - LAMBDA
         DO 70 I = B1, R1 - 1
            DPLUS = D( I ) + S
            IF(ABS(DPLUS).LT.PIVMIN) DPLUS = -PIVMIN
            WORK( INDLPL+I ) = LD( I ) / DPLUS
            IF(DPLUS.LT.ZERO) NEG1 = NEG1 + 1
            WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
            IF( WORK( INDLPL+I ).EQ.ZERO )
     $                      WORK( INDS+I ) = LLD( I )
            S = WORK( INDS+I ) - LAMBDA
 70      CONTINUE
         DO 71 I = R1, R2 - 1
            DPLUS = D( I ) + S
            IF(ABS(DPLUS).LT.PIVMIN) DPLUS = -PIVMIN
            WORK( INDLPL+I ) = LD( I ) / DPLUS
            WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
            IF( WORK( INDLPL+I ).EQ.ZERO ) 
     $                      WORK( INDS+I ) = LLD( I )
            S = WORK( INDS+I ) - LAMBDA
 71      CONTINUE
      END IF
*
*     Compute the progressive transform (using the differential form)
*     until the index R1
*
      SAWNAN2 = .FALSE.
      NEG2 = 0
      WORK( INDP+BN-1 ) = D( BN ) - LAMBDA
      DO 80 I = BN - 1, R1, -1
         DMINUS = LLD( I ) + WORK( INDP+I )
         TMP = D( I ) / DMINUS
         IF(DMINUS.LT.ZERO) NEG2 = NEG2 + 1
         WORK( INDUMN+I ) = L( I )*TMP
         WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - LAMBDA
 80   CONTINUE
      TMP = WORK( INDP+R1-1 )
      SAWNAN2 = DISNAN( TMP )	
      IF( SAWNAN2 ) THEN
*        Runs a slower version of the above loop if a NaN is detected
         NEG2 = 0
         DO 100 I = BN-1, R1, -1
            DMINUS = LLD( I ) + WORK( INDP+I )
            IF(ABS(DMINUS).LT.PIVMIN) DMINUS = -PIVMIN
            TMP = D( I ) / DMINUS
            IF(DMINUS.LT.ZERO) NEG2 = NEG2 + 1
            WORK( INDUMN+I ) = L( I )*TMP
            WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - LAMBDA
            IF( TMP.EQ.ZERO ) 
     $          WORK( INDP+I-1 ) = D( I ) - LAMBDA
 100     CONTINUE
      END IF
*
*     Find the index (from R1 to R2) of the largest (in magnitude)
*     diagonal element of the inverse
*
      MINGMA = WORK( INDS+R1-1 ) + WORK( INDP+R1-1 )
      IF( MINGMA.LT.ZERO ) NEG1 = NEG1 + 1
      IF( WANTNC ) THEN
         NEGCNT = NEG1 + NEG2
      ELSE
         NEGCNT = -1
      ENDIF
      IF( ABS(MINGMA).EQ.ZERO )
     $   MINGMA = EPS*WORK( INDS+R1-1 )
      R = R1
      DO 110 I = R1, R2 - 1
         TMP = WORK( INDS+I ) + WORK( INDP+I )
         IF( TMP.EQ.ZERO )
     $      TMP = EPS*WORK( INDS+I )
         IF( ABS( TMP ).LE.ABS( MINGMA ) ) THEN
            MINGMA = TMP
            R = I + 1
         END IF
 110  CONTINUE
*
*     Compute the FP vector: solve N^T v = e_r
*
      ISUPPZ( 1 ) = B1
      ISUPPZ( 2 ) = BN
      Z( R ) = ONE
      ZTZ = ONE
*
*     Compute the FP vector upwards from R
*
      NB = INT((R-B1)/BLKLEN)
      NX = R-NB*BLKLEN
      IF( .NOT.SAWNAN1 ) THEN
         DO 210 BI = R-1, NX, -BLKLEN
            TO = BI-BLKLEN+1
            DO 205 I = BI, TO, -1
               Z( I ) = -( WORK(INDLPL+I)*Z(I+1) )
               ZTZ = ZTZ + Z( I )*Z( I )
 205        CONTINUE
            IF( ABS(Z(TO)).LT.EPS .AND. 
     $        ABS(Z(TO+1)).LT.EPS ) THEN
               ISUPPZ(1) = TO
               GOTO 220
	    ENDIF
 210     CONTINUE
         DO 215 I = NX-1, B1, -1
            Z( I ) = -( WORK(INDLPL+I)*Z(I+1) )
            ZTZ = ZTZ + Z( I )*Z( I )
 215     CONTINUE
 220     CONTINUE
      ELSE
*        Run slower loop if NaN occurred.
         DO 230 BI = R-1, NX, -BLKLEN
            TO = BI-BLKLEN+1
            DO 225 I = BI, TO, -1
               IF( Z( I+1 ).EQ.ZERO ) THEN
                  Z( I ) = -( LD( I+1 ) / LD( I ) )*Z( I+2 )
               ELSE
                  Z( I ) = -( WORK( INDLPL+I )*Z( I+1 ) )
               END IF
               ZTZ = ZTZ + Z( I )*Z( I )
 225        CONTINUE
            IF( ABS(Z(TO)).LT.EPS .AND. 
     $        ABS(Z(TO+1)).LT.EPS ) THEN
               ISUPPZ(1) = TO
               GOTO 240
	    ENDIF
 230     CONTINUE
         DO 235 I = NX-1, B1, -1
            IF( Z( I+1 ).EQ.ZERO ) THEN
               Z( I ) = -( LD( I+1 ) / LD( I ) )*Z( I+2 )
            ELSE
               Z( I ) = -( WORK( INDLPL+I )*Z( I+1 ) )
            END IF
            ZTZ = ZTZ + Z( I )*Z( I )
 235     CONTINUE
 240     CONTINUE
      ENDIF
      DO 245 I= B1, (ISUPPZ(1)-1)
         Z(I) = ZERO
 245  CONTINUE
      
*     Compute the FP vector downwards from R in blocks of size BLKLEN
      IF( .NOT.SAWNAN2 ) THEN
         DO 260 BI = R+1, BN, BLKLEN
            TO = BI+BLKLEN-1
            IF ( TO.LE.BN ) THEN
               DO 250 I = BI, TO
                  Z(I) = -(WORK(INDUMN+I-1)*Z(I-1))
                  ZTZ = ZTZ + Z( I )*Z( I )
 250           CONTINUE   
               IF( ABS(Z(TO)).LE.EPS .AND. 
     $             ABS(Z(TO-1)).LE.EPS ) THEN
                  ISUPPZ(2) = TO
                  GOTO 265
	       ENDIF
            ELSE
               DO 255 I = BI, BN
                  Z(I) = -(WORK(INDUMN+I-1)*Z(I-1))
                  ZTZ = ZTZ + Z( I )*Z( I )
 255           CONTINUE   
            ENDIF
 260     CONTINUE
 265     CONTINUE
      ELSE
*        Run slower loop if NaN occurred.
         DO 280 BI = R+1, BN, BLKLEN
            TO = BI+BLKLEN-1
            IF ( TO.LE.BN ) THEN
               DO 270 I = BI, TO
                  ZPREV = Z(I-1)
                  ABSZPREV = ABS(ZPREV)
                  IF( ZPREV.NE.ZERO ) THEN
                     Z(I)= -(WORK(INDUMN+I-1)*ZPREV)
                  ELSE
                     Z(I)= -(LD(I-2)/LD(I-1))*Z(I-2)
                  END IF
                  ABSZCUR = ABS(Z(I))
                  ZTZ = ZTZ + ABSZCUR**2
 270           CONTINUE
               IF( ABSZCUR.LT.EPS .AND. 
     $             ABSZPREV.LT.EPS ) THEN
                  ISUPPZ(2) = I
                  GOTO 285
	       ENDIF
            ELSE
               DO 275 I = BI, BN
                  ZPREV = Z(I-1)
                  ABSZPREV = ABS(ZPREV)
                  IF( ZPREV.NE.ZERO ) THEN
                     Z(I)= -(WORK(INDUMN+I-1)*ZPREV)
                  ELSE
                     Z(I)= -(LD(I-2)/LD(I-1))*Z(I-2)
                  END IF
                  ABSZCUR = ABS(Z(I))
                  ZTZ = ZTZ + ABSZCUR**2
 275           CONTINUE
            ENDIF
 280     CONTINUE
 285     CONTINUE
      END IF
      DO 290 I= ISUPPZ(2)+1,BN
         Z(I) = ZERO
 290  CONTINUE
*
*     Compute quantities for convergence test
*     
      TMP = ONE / ZTZ
      NRMINV = SQRT( TMP )
      RESID = ABS( MINGMA )*NRMINV
      RQCORR = MINGMA*TMP
*
      RETURN
*
*     End of DLAR1VA
*
      END
