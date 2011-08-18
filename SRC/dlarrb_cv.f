      SUBROUTINE DLARRB_CV( N, D, LLD, IFIRST, ILAST, RTOL1,
     $                   RTOL2, OFFSET, W, WGAP, WERR, WORK, IWORK,
     $                   PIVMIN, SPDIAM, TWIST, INFO )
*
*  -- LAPACK computational routine (version *TBA*) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     July 4, 2010
*
*     .. Scalar Arguments ..
      INTEGER            IFIRST, ILAST, INFO, N, OFFSET, TWIST
      DOUBLE PRECISION   PIVMIN, RTOL1, RTOL2, SPDIAM
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), LLD( * ), W( * ),
     $                   WERR( * ), WGAP( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  Given the relatively robust representation(RRR) L D L^T, DLARRB
*  does "limited" bisection to refine the eigenvalues of L D L^T,
*  W( IFIRST-OFFSET ) thru' W( ILAST-OFFSET ), to more accuracy. Initial
*  guesses for these eigenvalues are input in W, the corresponding estimate
*  of the error in these guesses and their gaps are input in WERR
*  and WGAP, respectively. During bisection, intervals
*  [left, right] are maintained by storing their mid-points and
*  semi-widths in the arrays W and WERR respectively.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The N diagonal elements of the diagonal matrix D.
*
*  LLD     (input) DOUBLE PRECISION array, dimension (N-1)
*          The (N-1) elements L(i)*L(i)*D(i).
*
*  IFIRST  (input) INTEGER
*          The index of the first eigenvalue to be computed.
*
*  ILAST   (input) INTEGER
*          The index of the last eigenvalue to be computed.
*
*  RTOL1   (input) DOUBLE PRECISION
*  RTOL2   (input) DOUBLE PRECISION
*          Tolerance for the convergence of the bisection intervals.
*          An interval [LEFT,RIGHT] has converged if
*          RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )
*          where GAP is the (estimated) distance to the nearest
*          eigenvalue.
*
*  OFFSET  (input) INTEGER
*          Offset for the arrays W, WGAP and WERR, i.e., the IFIRST-OFFSET
*          thru' ILAST-OFFSET elements of these arrays are to be used.
*
*  W       (input/output) DOUBLE PRECISION array, dimension (N)
*          On input, W( IFIRST-OFFSET ) thru' W( ILAST-OFFSET ) are
*          estimates of the eigenvalues of L D L^T indexed IFIRST thru' ILAST.
*          On output, these estimates are refined.
*
*  WGAP    (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On input, the (estimated) gaps between consecutive
*          eigenvalues of L D L^T, i.e., WGAP(I-OFFSET) is the gap between
*          eigenvalues I and I+1. Note that if IFIRST.EQ.ILAST
*          then WGAP(IFIRST-OFFSET) must be set to ZERO.
*          On output, these gaps are refined.
*
*  WERR    (input/output) DOUBLE PRECISION array, dimension (N)
*          On input, WERR( IFIRST-OFFSET ) thru' WERR( ILAST-OFFSET ) are
*          the errors in the estimates of the corresponding elements in W.
*          On output, these errors are refined.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
*          Workspace.
*
*  IWORK   (workspace) INTEGER array, dimension (2*N)
*          Workspace.
*
*  PIVMIN  (input) DOUBLE PRECISION
*          The minimum pivot in the sturm sequence.
*
*  SPDIAM  (input) DOUBLE PRECISION
*          The spectral diameter of the matrix.
*
*  TWIST   (input) INTEGER
*          The twist index for the twisted factorization that is used
*          for the negcount. 
*          TWIST = N: Compute negcount from L D L^T - LAMBDA I = L+ D+ L+^T
*          TWIST = 1: Compute negcount from L D L^T - LAMBDA I = U- D- U-^T
*          TWIST = R: Compute negcount from L D L^T - LAMBDA I = N(r) D(r) N(r)
*
*  INFO    (output) INTEGER
*          Error flag.
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
      DOUBLE PRECISION   ZERO, TWO, HALF
      PARAMETER        ( ZERO = 0.0D0, TWO = 2.0D0,
     $                   HALF = 0.5D0 )
      INTEGER   MAXITR
*     ..
*     .. Local Scalars ..
      INTEGER            I, I1, II, IP, ITER, K, NEGCNT, NEXT, NINT,
     $                   OLNINT, PREV, R
      DOUBLE PRECISION   BACK, CVRGD, GAP, LEFT, LGAP, MID, MNWDTH,
     $                   RGAP, RIGHT, SAVGAP, TMP, WIDTH
      LOGICAL   PARANOID
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN
      DOUBLE PRECISION   DLAMCH
      INTEGER            DLANEG_CV
      EXTERNAL           DISNAN, DLAMCH, DLANEG_CV
*
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Turn on paranoid check for rounding errors
*     invalidating uncertainty intervals of eigenvalues
*
      PARANOID = .TRUE.
*     
      MAXITR = INT( ( LOG( SPDIAM+PIVMIN )-LOG( PIVMIN ) ) /
     $           LOG( TWO ) ) + 2
      MNWDTH = TWO * PIVMIN
*
      R = TWIST
      IF((R.LT.1).OR.(R.GT.N)) R = N
*
*     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
*     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while
*     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 )
*     for an unconverged interval is set to the index of the next unconverged
*     interval, and is -1 or 0 for a converged interval. Thus a linked
*     list of unconverged intervals is set up.
*
      I1 = IFIRST
*     The number of unconverged intervals 
      NINT = 0
*     The last unconverged interval found
      PREV = 0
     
      RGAP = WGAP( I1-OFFSET )
      DO 75 I = I1, ILAST
         K = 2*I
         II = I - OFFSET
         LEFT = W( II ) - WERR( II )
         RIGHT = W( II ) + WERR( II )

*        Check that interval ends are not close to underflow
*        It will mess up the qds algorithms
*        This typically happens when the random perturbation fails
*        to separate eigenvalues that agree up to the underflow threshold
         IF ((ABS(LEFT).LE.16*PIVMIN).OR.(ABS(RIGHT).LE.16*PIVMIN))
     $      THEN
            INFO = -1
            RETURN
         ENDIF
 
         LGAP = RGAP
         RGAP = WGAP( II )
         GAP = MIN( LGAP, RGAP )

         IF( PARANOID ) THEN
*        Make sure that [LEFT,RIGHT] contains the desired eigenvalue
*        Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT 
*	 
*        Do while( NEGCNT(LEFT).GT.I-1 )
*	 
         BACK = WERR( II )
 20      CONTINUE
         NEGCNT = DLANEG_CV( N, D, LLD, LEFT, PIVMIN, R )
         IF( NEGCNT.GT.I-1 ) THEN
            LEFT = LEFT - BACK
            BACK = TWO*BACK
            GO TO 20
         END IF
*
*        Do while( NEGCNT(RIGHT).LT.I )
*        Compute negcount from dstqds facto L+D+L+^T = L D L^T - RIGHT 
*	 
         BACK = WERR( II )
 50      CONTINUE

         NEGCNT = DLANEG_CV( N, D, LLD, RIGHT, PIVMIN, R )
         IF( NEGCNT.LT.I ) THEN
             RIGHT = RIGHT + BACK
             BACK = TWO*BACK
             GO TO 50
         END IF
         ENDIF
*
         WIDTH = HALF*ABS( LEFT - RIGHT )
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )
         CVRGD = MAX(RTOL1*GAP,RTOL2*TMP)
         IF( WIDTH.LE.CVRGD .OR. WIDTH.LE.MNWDTH ) THEN
*           This interval has already converged and does not need refinement.
*           (Note that the gaps might change through refining the 
*            eigenvalues, however, they can only get bigger.)
*           Remove it from the list.
            IWORK( K-1 ) = -1
*           Make sure that I1 always points to the first unconverged interval
            IF((I.EQ.I1).AND.(I.LT.ILAST)) I1 = I + 1
            IF((PREV.GE.I1).AND.(I.LE.ILAST)) IWORK( 2*PREV-1 ) = I + 1
         ELSE
*           unconverged interval found
            PREV = I
            NINT = NINT + 1
            IWORK( K-1 ) = I + 1
            IWORK( K ) = NEGCNT
         END IF
         WORK( K-1 ) = LEFT
         WORK( K ) = RIGHT
 75   CONTINUE

*
*     Do while( NINT.GT.0 ), i.e. there are still unconverged intervals
*     and while (ITER.LT.MAXITR)
*
      ITER = 0 
 80   CONTINUE
      PREV = I1 - 1
      I = I1
      OLNINT = NINT

      DO 100 IP = 1, OLNINT
         K = 2*I
         II = I - OFFSET
         RGAP = WGAP( II )
         LGAP = RGAP
         IF(II.GT.1) LGAP = WGAP( II-1 ) 
         GAP = MIN( LGAP, RGAP )
         NEXT = IWORK( K-1 )
         LEFT = WORK( K-1 )
         RIGHT = WORK( K )
         MID = HALF*( LEFT + RIGHT ) 
*        semiwidth of interval
         WIDTH = RIGHT - MID
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )
         CVRGD = MAX(RTOL1*GAP,RTOL2*TMP)
         IF( ( WIDTH.LE.CVRGD ) .OR. ( WIDTH.LE.MNWDTH ).OR.
     $       ( ITER.EQ.MAXITR ) )THEN
*           reduce number of unconverged intervals
            NINT = NINT - 1
*           Mark interval as converged. 
            IWORK( K-1 ) = 0
            IF( I1.EQ.I ) THEN
               I1 = NEXT
            ELSE
*              Prev holds the last unconverged interval previously examined
               IF(PREV.GE.I1) IWORK( 2*PREV-1 ) = NEXT
            END IF
            I = NEXT
            GO TO 100
         END IF
         PREV = I

*
*        Perform one bisection step
*
         NEGCNT = DLANEG_CV( N, D, LLD, MID, PIVMIN, R )
         IF( NEGCNT.LE.I-1 ) THEN
            WORK( K-1 ) = MID
         ELSE
            WORK( K ) = MID
         END IF
         I = NEXT
 100  CONTINUE
      ITER = ITER + 1
*     do another loop if there are still unconverged intervals
*     However, in the last iteration, all intervals are accepted
*     since this is the best we can do.
      IF( ( NINT.GT.0 ).AND.(ITER.LE.MAXITR) ) GO TO 80
*
*
*     At this point, all the intervals have converged
*
*     save this gap to restore it after the loop
      SAVGAP = WGAP( ILAST-OFFSET )
*
      LEFT = WORK( 2*IFIRST-1 )
      DO 110 I = IFIRST, ILAST
         K = 2*I
         II = I - OFFSET
*        RIGHT is the right boundary of this current interval
         RIGHT = WORK( K ) 
*        All intervals marked by '0' have been refined.
         IF( IWORK( K-1 ).EQ.0 ) THEN
            W( II ) = HALF*( LEFT+RIGHT )
            WERR( II ) = RIGHT - W( II )
         END IF
*        Left is the boundary of the next interval
         LEFT = WORK( K +1 ) 
         WGAP( II ) = MAX( ZERO, LEFT - RIGHT )
 110  CONTINUE
*     restore the last gap which was overwritten by garbage
      WGAP( ILAST-OFFSET ) = SAVGAP

      RETURN
*
*     End of DLARRB
*
      END
*
*
*
      FUNCTION DLANEG_CV( N, D, LLD, SIGMA, PIVMIN, R )
      INTEGER DLANEG_CV
*
*     .. Scalar Arguments ..
      INTEGER            N, R
      DOUBLE PRECISION   PIVMIN, SIGMA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), LLD( * )
*
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, NEG1, NEG2, NEGCNT
      DOUBLE PRECISION   DMINUS, DPLUS, GAMMA, P, S, T, TMP
      LOGICAL SAWNAN
*     ..
*     .. External Functions ..
      LOGICAL DISNAN
      EXTERNAL DISNAN
      
      NEGCNT = 0
      
*     I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
      NEG1 = 0
      S = ZERO
      DO 21 J = 1, R - 1
         T = S - SIGMA
         DPLUS = D( J ) + T
         S = T*LLD( J ) / DPLUS 
         IF( DPLUS.LT.ZERO ) NEG1 = NEG1 + 1
 21   CONTINUE
      SAWNAN = DISNAN( S )
*     Run a slower version of the above loop if a NaN is detected
      IF( SAWNAN ) THEN
         NEG1 = 0
         S = ZERO
         T = -SIGMA
         DO 22 J = 1, R - 1
            DPLUS = D( J ) + T
            IF(ABS(DPLUS).LT.PIVMIN) DPLUS = -PIVMIN
            TMP = LLD( J ) / DPLUS
            IF( DPLUS.LT.ZERO ) NEG1 = NEG1 + 1
            S = T*TMP
            IF( TMP.EQ.ZERO ) S = LLD( J )
            T = S - SIGMA
 22      CONTINUE
      END IF
      NEGCNT = NEGCNT + NEG1
*     
*     II) lower part: L D L^T - SIGMA I = U- D- U-^T
      NEG2 = 0
      P = D( N ) - SIGMA
      DO 23 J = N - 1, R, -1
         DMINUS = LLD( J ) + P
         P = P*D( J )/DMINUS - SIGMA
         IF( DMINUS.LT.ZERO ) NEG2 = NEG2 + 1
 23   CONTINUE
      SAWNAN = DISNAN( P )
      IF( SAWNAN ) THEN
         NEG2 = 0
         P = D( N ) - SIGMA
         DO 24 J = N - 1, R, -1
            DMINUS = LLD( J ) + P
            IF(ABS(DMINUS).LT.PIVMIN) DMINUS = -PIVMIN
            TMP = D( J ) / DMINUS
            IF( DMINUS.LT.ZERO ) NEG2 = NEG2 + 1
            P = P*TMP - SIGMA
            IF( TMP.EQ.ZERO ) P = D( J ) - SIGMA
 24      CONTINUE
      END IF
      NEGCNT = NEGCNT + NEG2
*     
*     III) Twist index
      GAMMA = S + P
      IF( GAMMA.LT.ZERO ) NEGCNT = NEGCNT+1

      DLANEG_CV = NEGCNT
      END
