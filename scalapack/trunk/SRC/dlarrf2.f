      SUBROUTINE DLARRF2( N, D, L, LD, CLSTRT, CLEND, 
     $                   CLMID1, CLMID2, W, WGAP, WERR, TRYMID,
     $                   SPDIAM, CLGAPL, CLGAPR, PIVMIN, SIGMA,
     $                   DPLUS, LPLUS, WORK, INFO )
*
*  -- ScaLAPACK computational routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ of Colorado Denver
*     July 4, 2010
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            CLSTRT, CLEND, CLMID1, CLMID2, INFO, N
      DOUBLE PRECISION   CLGAPL, CLGAPR, PIVMIN, SIGMA, SPDIAM
      LOGICAL TRYMID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DPLUS( * ), L( * ), LD( * ), 
     $          LPLUS( * ), W( * ), WGAP( * ), WERR( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  Given the initial representation L D L^T and its cluster of close
*  eigenvalues (in a relative measure), W( CLSTRT ), W( CLSTRT+1 ), ...
*  W( CLEND ), DLARRF2 finds a new relatively robust representation
*  L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the
*  eigenvalues of L(+) D(+) L(+)^T is relatively isolated.
*
*  This is an enhanced version of DLARRF that also tries shifts in
*  the middle of the cluster, should there be a large gap, in order to
*  break large clusters into at least two pieces.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix (subblock, if the matrix splitted).
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The N diagonal elements of the diagonal matrix D.
*
*  L       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (N-1) subdiagonal elements of the unit bidiagonal
*          matrix L.
*
*  LD      (input) DOUBLE PRECISION array, dimension (N-1)
*          The (N-1) elements L(i)*D(i).
*
*  CLSTRT  (input) INTEGER
*          The index of the first eigenvalue in the cluster.
*
*  CLEND   (input) INTEGER
*          The index of the last eigenvalue in the cluster.
*
*  CLMID1,2(input) INTEGER
*          The index of a middle eigenvalue pair with large gap
*
*  W       (input) DOUBLE PRECISION array, dimension >=  (CLEND-CLSTRT+1)
*          The eigenvalue APPROXIMATIONS of L D L^T in ascending order.
*          W( CLSTRT ) through W( CLEND ) form the cluster of relatively
*          close eigenalues.
*
*  WGAP    (input/output) DOUBLE PRECISION array, dimension >=  (CLEND-CLSTRT+1)
*          The separation from the right neighbor eigenvalue in W.
*
*  WERR    (input) DOUBLE PRECISION array, dimension >=  (CLEND-CLSTRT+1)
*          WERR contain the semiwidth of the uncertainty
*          interval of the corresponding eigenvalue APPROXIMATION in W
*
*  SPDIAM (input) estimate of the spectral diameter obtained from the
*          Gerschgorin intervals
*
*  CLGAPL, CLGAPR (input) absolute gap on each end of the cluster.
*          Set by the calling routine to protect against shifts too close
*          to eigenvalues outside the cluster.
*
*  PIVMIN  (input) DOUBLE PRECISION
*          The minimum pivot allowed in the sturm sequence.
*
*  SIGMA   (output) DOUBLE PRECISION
*          The shift used to form L(+) D(+) L(+)^T.
*
*  DPLUS   (output) DOUBLE PRECISION array, dimension (N)
*          The N diagonal elements of the diagonal matrix D(+).
*
*  LPLUS   (output) DOUBLE PRECISION array, dimension (N-1)
*          The first (N-1) elements of LPLUS contain the subdiagonal
*          elements of the unit bidiagonal matrix L(+).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
*          Workspace.
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
      DOUBLE PRECISION   FOUR, MAXGROWTH1, MAXGROWTH2, ONE, QUART, TWO
      PARAMETER          ( ONE = 1.0D0, TWO = 2.0D0,
     $                     FOUR = 4.0D0, QUART = 0.25D0,
     $                     MAXGROWTH1 = 8.D0,
     $                     MAXGROWTH2 = 8.D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL   DORRR1, NOFAIL, SAWNAN1, SAWNAN2, TRYRRR1
      INTEGER      BI,I,J,KTRY,KTRYMAX,SLEFT,SRIGHT,SMID,SHIFT
      PARAMETER   ( KTRYMAX = 1, SMID =0, SLEFT = 1, SRIGHT = 2 )

*     DSTQDS loops will be blocked to detect NaNs earlier if they occur
      INTEGER BLKLEN
      PARAMETER ( BLKLEN = 512 )


      DOUBLE PRECISION   AVGAP, BESTSHIFT, CLWDTH, EPS, FACT, FAIL,
     $                   FAIL2, GROWTHBOUND, LDELTA, LDMAX, LEASTGROWTH,
     $                   LSIGMA, MAX1, MAX2, MINGAP, MSIGMA1, MSIGMA2,
     $                   OLDP, PROD, RDELTA, RDMAX, RRR1, RRR2, RSIGMA,
     $                   S, TMP, ZNM2
*     ..
*     .. External Functions ..
      LOGICAL DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DISNAN, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      FACT = DBLE(2**KTRYMAX)
      EPS = DLAMCH( 'Precision' )
      SHIFT = 0
      
*     Decide whether the code should accept the best among all 
*     representations despite large element growth or signal INFO=1
      NOFAIL = .TRUE.
*

*     Compute the average gap length of the cluster
      CLWDTH = ABS(W(CLEND)-W(CLSTRT)) + WERR(CLEND) + WERR(CLSTRT)
      AVGAP = CLWDTH / DBLE(CLEND-CLSTRT)
      MINGAP = MIN(CLGAPL, CLGAPR)

*     Initial values for shifts to both ends of cluster
      LSIGMA = MIN(W( CLSTRT ),W( CLEND )) - WERR( CLSTRT )
      RSIGMA = MAX(W( CLSTRT ),W( CLEND )) + WERR( CLEND )
      MSIGMA1 = W( CLMID1 ) + WERR( CLMID1 )
      MSIGMA2 = W( CLMID2 ) - WERR( CLMID2 )

*     Use a small fudge to make sure that we really shift to the outside
      LSIGMA = LSIGMA - ABS(LSIGMA)* FOUR * EPS
      RSIGMA = RSIGMA + ABS(RSIGMA)* FOUR * EPS

*     Compute upper bounds for how much to back off the initial shifts
      LDMAX = QUART * MINGAP + TWO * PIVMIN
      RDMAX = QUART * MINGAP + TWO * PIVMIN
	
      LDELTA = MAX(AVGAP,WGAP( CLSTRT ))/FACT
      RDELTA = MAX(AVGAP,WGAP( CLEND-1 ))/FACT
*
*     Initialize the record of the best representation found
*
      S = DLAMCH( 'S' )
      LEASTGROWTH = ONE / S 
      FAIL = DBLE(N-1)*MINGAP/(SPDIAM*EPS)
      FAIL2 = DBLE(N-1)*MINGAP/(SPDIAM*SQRT(EPS))
      GROWTHBOUND = MAXGROWTH1*SPDIAM

*
*     Set default best shift
*
      BESTSHIFT = LSIGMA


      IF(.NOT.TRYMID) GOTO 4
*
*     Try shifts in the middle
*     
      SHIFT = SMID

      DO 3 J=1,2
         SAWNAN1 = .FALSE.
         IF(J.EQ.1) THEN
*           Try left middle point
            SIGMA = MSIGMA1
         ELSE
*           Try left middle point
            SIGMA = MSIGMA2
	 ENDIF   
 
         S = -SIGMA
         DPLUS( 1 ) = D( 1 ) + S
         MAX1 = ABS( DPLUS( 1 ) )
         DO 2 BI = 1, N-1, BLKLEN
            DO 1 I = BI, MIN( BI+BLKLEN-1, N-1)
               LPLUS( I ) = LD( I ) / DPLUS( I )
               S = S*LPLUS( I )*L( I ) - SIGMA
               DPLUS( I+1 ) = D( I+1 ) + S
               MAX1 = MAX( MAX1,ABS(DPLUS(I+1)) )
 1          CONTINUE
            SAWNAN1=SAWNAN1 .OR. DISNAN(MAX1)
            IF (SAWNAN1) GOTO 3
 2       CONTINUE

         IF( .NOT.SAWNAN1 ) THEN
            IF( MAX1.LE.GROWTHBOUND ) THEN
               GOTO 100
            ELSE IF( MAX1.LE.LEASTGROWTH ) THEN           
               LEASTGROWTH = MAX1
               BESTSHIFT = SIGMA
            ENDIF
         ENDIF
 3    CONTINUE


 4    CONTINUE
*
*     Shifts in the middle not tried or not succeeded
*     Find best shift on the outside of the cluster
*
*     while (KTRY <= KTRYMAX)
      KTRY = 0 
*
*
*
 5    CONTINUE

*     Compute element growth when shifting to both ends of the cluster
*     accept shift if there is no element growth at one of the two ends

*     Left end
      SAWNAN1 = .FALSE.
      S = -LSIGMA
      DPLUS( 1 ) = D( 1 ) + S
      MAX1 = ABS( DPLUS( 1 ) )
      DO 12 BI = 1, N-1, BLKLEN
         DO 11 I = BI, MIN( BI+BLKLEN-1, N-1)
            LPLUS( I ) = LD( I ) / DPLUS( I )
            S = S*LPLUS( I )*L( I ) - LSIGMA
            DPLUS( I+1 ) = D( I+1 ) + S
            MAX1 = MAX( MAX1,ABS(DPLUS(I+1)) )
 11      CONTINUE
         SAWNAN1=SAWNAN1 .OR. DISNAN(MAX1)
         IF (SAWNAN1) GOTO 13
 12   CONTINUE
      IF( .NOT.SAWNAN1 ) THEN
         IF( MAX1.LE.GROWTHBOUND ) THEN
            SIGMA = LSIGMA
            SHIFT = SLEFT
            GOTO 100
         ELSE IF( MAX1.LE.LEASTGROWTH ) THEN           
            LEASTGROWTH = MAX1
            BESTSHIFT = LSIGMA
         ENDIF
      ENDIF
 13   CONTINUE

*     Right end      
      SAWNAN2 = .FALSE.
      S = -RSIGMA
      WORK( 1 ) = D( 1 ) + S
      MAX2 = ABS( WORK( 1 ) )
      DO 22 BI = 1, N-1, BLKLEN
         DO 21 I = BI, MIN( BI+BLKLEN-1, N-1)
            WORK( N+I ) = LD( I ) / WORK( I )
            S = S*WORK( N+I )*L( I ) - RSIGMA
            WORK( I+1 ) = D( I+1 ) + S
            MAX2 = MAX( MAX2,ABS(WORK(I+1)) )
 21      CONTINUE
         SAWNAN2=SAWNAN2 .OR. DISNAN(MAX2)
         IF (SAWNAN2) GOTO 23
 22   CONTINUE
      IF( .NOT.SAWNAN2 ) THEN
         IF( MAX2.LE.GROWTHBOUND ) THEN
            SIGMA = RSIGMA
	    SHIFT = SRIGHT
            GOTO 100
         ELSE IF( MAX2.LE.LEASTGROWTH ) THEN           
            LEASTGROWTH = MAX2
            BESTSHIFT = RSIGMA
         ENDIF
      ENDIF
 23   CONTINUE

*     If we are at this point, both shifts led to too much element growth

 50   CONTINUE

      IF (KTRY.LT.KTRYMAX) THEN
*        If we are here, both shifts failed also the RRR test.
*        Back off to the outside      
         LSIGMA = MAX( LSIGMA - LDELTA, 
     $     LSIGMA - LDMAX)
         RSIGMA = MIN( RSIGMA + RDELTA, 
     $     RSIGMA + RDMAX )
         LDELTA = TWO * LDELTA      
         RDELTA = TWO * RDELTA
*        Ensure that we do not back off too much of the initial shifts
         LDELTA = MIN(LDMAX,LDELTA)
         RDELTA = MIN(RDMAX,RDELTA)
         KTRY = KTRY + 1
         GOTO 5
      ELSE     
*        None of the representations investigated satisfied our
*        criteria. Take the best one we found.
         IF((LEASTGROWTH.LT.FAIL).OR.NOFAIL) THEN
            LSIGMA = BESTSHIFT
            SAWNAN1 = .FALSE.
            S = -LSIGMA
            DPLUS( 1 ) = D( 1 ) + S
            DO 6 I = 1, N - 1
               LPLUS( I ) = LD( I ) / DPLUS( I )
               S = S*LPLUS( I )*L( I ) - LSIGMA
               DPLUS( I+1 ) = D( I+1 ) + S
               IF(ABS(DPLUS(I+1)).LT.PIVMIN) THEN
                  DPLUS(I+1) = -PIVMIN
               ENDIF
 6          CONTINUE
            SIGMA = LSIGMA
    	    SHIFT = SLEFT
            GOTO 100
         ELSE
            INFO = 1
            RETURN
         ENDIF
      END IF           

 100  CONTINUE
      IF (SHIFT.EQ.SLEFT .OR. SHIFT.EQ.SMID ) THEN
      ELSEIF (SHIFT.EQ.SRIGHT) THEN
*        store new L and D back into DPLUS, LPLUS
         CALL DCOPY( N, WORK, 1, DPLUS, 1 )
         CALL DCOPY( N-1, WORK(N+1), 1, LPLUS, 1 )
      ENDIF

      RETURN
*
*     End of DLARRF2
*
      END
