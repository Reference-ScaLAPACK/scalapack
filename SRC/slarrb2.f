      SUBROUTINE SLARRB2( N, D, LLD, IFIRST, ILAST, RTOL1,
     $                   RTOL2, OFFSET, W, WGAP, WERR, WORK, IWORK,
     $                   PIVMIN, LGPVMN, LGSPDM, TWIST, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     July 4, 2010
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            IFIRST, ILAST, INFO, N, OFFSET, TWIST
      REAL               LGPVMN, LGSPDM, PIVMIN, 
     $                   RTOL1, RTOL2
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               D( * ), LLD( * ), W( * ),
     $                   WERR( * ), WGAP( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  Given the relatively robust representation(RRR) L D L^T, SLARRB2
*  does "limited" bisection to refine the eigenvalues of L D L^T,
*  W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial
*  guesses for these eigenvalues are input in W, the corresponding estimate
*  of the error in these guesses and their gaps are input in WERR
*  and WGAP, respectively. During bisection, intervals
*  [left, right] are maintained by storing their mid-points and
*  semi-widths in the arrays W and WERR respectively.
*
*  NOTE: 
*  There are very few minor differences between SLARRB from LAPACK
*  and this current subroutine SLARRB2.
*  The most important reason for creating this nearly identical copy
*  is profiling: in the ScaLAPACK MRRR algorithm, eigenvalue computation 
*  using SLARRB2 is used for refinement in the construction of 
*  the representation tree, as opposed to the initial computation of the
*  eigenvalues for the root RRR which uses SLARRB. When profiling,
*  this allows an easy quantification of refinement work vs. computing
*  eigenvalues of the root.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  D       (input) REAL             array, dimension (N)
*          The N diagonal elements of the diagonal matrix D.
*
*  LLD     (input) REAL             array, dimension (N-1)
*          The (N-1) elements L(i)*L(i)*D(i).
*
*  IFIRST  (input) INTEGER
*          The index of the first eigenvalue to be computed.
*
*  ILAST   (input) INTEGER
*          The index of the last eigenvalue to be computed.
*
*  RTOL1   (input) REAL            
*  RTOL2   (input) REAL            
*          Tolerance for the convergence of the bisection intervals.
*          An interval [LEFT,RIGHT] has converged if
*          RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )
*          where GAP is the (estimated) distance to the nearest
*          eigenvalue.
*
*  OFFSET  (input) INTEGER
*          Offset for the arrays W, WGAP and WERR, i.e., the IFIRST-OFFSET
*          through ILAST-OFFSET elements of these arrays are to be used.
*
*  W       (input/output) REAL             array, dimension (N)
*          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are
*          estimates of the eigenvalues of L D L^T indexed IFIRST through ILAST.
*          On output, these estimates are refined.
*
*  WGAP    (input/output) REAL             array, dimension (N-1)
*          On input, the (estimated) gaps between consecutive
*          eigenvalues of L D L^T, i.e., WGAP(I-OFFSET) is the gap between
*          eigenvalues I and I+1. Note that if IFIRST.EQ.ILAST
*          then WGAP(IFIRST-OFFSET) must be set to ZERO.
*          On output, these gaps are refined.
*
*  WERR    (input/output) REAL             array, dimension (N)
*          On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are
*          the errors in the estimates of the corresponding elements in W.
*          On output, these errors are refined.
*
*  WORK    (workspace) REAL             array, dimension (4*N)
*          Workspace.
*
*  IWORK   (workspace) INTEGER array, dimension (2*N)
*          Workspace.
*
*  PIVMIN  (input) REAL             
*          The minimum pivot in the sturm sequence.
*
*  LGPVMN  (input) REAL            
*          Logarithm of PIVMIN, precomputed.
*
*  LGSPDM  (input) REAL             
*          Logarithm of the spectral diameter, precomputed.
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
*     .. Parameters ..
      REAL               ZERO, TWO, HALF
      PARAMETER        ( ZERO = 0.0E0, TWO = 2.0E0,
     $                   HALF = 0.5E0 )
      INTEGER   MAXITR
*     ..
*     .. Local Scalars ..
      INTEGER            I, I1, II, INDLLD, IP, ITER, J, K, NEGCNT,
     $                   NEXT, NINT, OLNINT, PREV, R
      REAL               BACK, CVRGD, GAP, LEFT, LGAP, MID, MNWDTH,
     $                   RGAP, RIGHT, SAVGAP, TMP, WIDTH
      LOGICAL   PARANOID
*     ..
*     .. External Functions ..
      LOGICAL            SISNAN
      REAL               SLAMCH
      INTEGER            SLANEG2A
      EXTERNAL           SISNAN, SLAMCH, 
     $                   SLANEG2A
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
      MAXITR = INT( ( LGSPDM - LGPVMN ) / LOG( TWO ) ) + 2
      MNWDTH = TWO * PIVMIN
*
      R = TWIST
*
      INDLLD = 2*N     
      DO 5 J = 1, N-1 
         I=2*J
         WORK(INDLLD+I-1) = D(J)
         WORK(INDLLD+I) = LLD(J)
  5   CONTINUE
      WORK(INDLLD+2*N-1) = D(N)
*
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
         LGAP = RGAP
         RGAP = WGAP( II )
         GAP = MIN( LGAP, RGAP )

         IF((ABS(LEFT).LE.16*PIVMIN).OR.(ABS(RIGHT).LE.16*PIVMIN))
     $      THEN
            INFO = -1
            RETURN
         ENDIF

         IF( PARANOID ) THEN
*        Make sure that [LEFT,RIGHT] contains the desired eigenvalue
*        Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT 
*	 
*        Do while( NEGCNT(LEFT).GT.I-1 )
*	 
         BACK = WERR( II )
 20      CONTINUE
         NEGCNT = SLANEG2A( N, WORK(INDLLD+1), LEFT, PIVMIN, R )
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
         NEGCNT = SLANEG2A( N, WORK(INDLLD+1),RIGHT, PIVMIN, R )

         IF( NEGCNT.LT.I ) THEN
             RIGHT = RIGHT + BACK
             BACK = TWO*BACK
             GO TO 50
         END IF
         ENDIF

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
         NEGCNT = SLANEG2A( N, WORK(INDLLD+1), MID, PIVMIN, R )
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
*     End of SLARRB2
*
      END
*
*
*
      FUNCTION SLANEG2( N, D, LLD, SIGMA, PIVMIN, R )
*
      IMPLICIT NONE
*
      INTEGER SLANEG2
*
*     .. Scalar Arguments ..
      INTEGER            N, R
      REAL               PIVMIN, SIGMA
*     ..
*     .. Array Arguments ..
      REAL               D( * ), LLD( * )
*
      REAL               ZERO
      PARAMETER        ( ZERO = 0.0E0 )

      INTEGER BLKLEN
      PARAMETER ( BLKLEN = 2048 )
*     ..
*     .. Local Scalars ..
      INTEGER            BJ, J, NEG1, NEG2, NEGCNT, TO
      REAL               DMINUS, DPLUS, GAMMA, P, S, T, TMP, XSAV
      LOGICAL SAWNAN
*     ..
*     .. External Functions ..
      LOGICAL SISNAN
      EXTERNAL SISNAN
      
      NEGCNT = 0
*      
*     I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
*     run dstqds block-wise to avoid excessive work when NaNs occur 
*
      S = ZERO
      DO 210 BJ = 1, R-1, BLKLEN
         NEG1 = 0
         XSAV = S
         TO = BJ+BLKLEN-1 
         IF ( TO.LE.R-1 ) THEN
            DO 21 J = BJ, TO
               T = S - SIGMA
               DPLUS = D( J ) + T
               IF( DPLUS.LT.ZERO ) NEG1=NEG1 + 1
               S = T*LLD( J ) / DPLUS 
 21         CONTINUE
         ELSE
            DO 22 J = BJ, R-1
               T = S - SIGMA
               DPLUS = D( J ) + T
               IF( DPLUS.LT.ZERO ) NEG1=NEG1 + 1
               S = T*LLD( J ) / DPLUS 
 22         CONTINUE
         ENDIF
         SAWNAN = SISNAN( S )
*         
         IF( SAWNAN ) THEN
            NEG1 = 0
            S = XSAV
            TO = BJ+BLKLEN-1 
            IF ( TO.LE.R-1 ) THEN
               DO 23 J = BJ, TO
                  T = S - SIGMA
                  DPLUS = D( J ) + T
                  IF(ABS(DPLUS).LT.PIVMIN) 
     $               DPLUS = -PIVMIN
                  TMP = LLD( J ) / DPLUS
                  IF( DPLUS.LT.ZERO ) 
     $               NEG1 = NEG1 + 1
                  S = T*TMP
                  IF( TMP.EQ.ZERO ) S = LLD( J )
 23            CONTINUE
            ELSE
               DO 24 J = BJ, R-1
                  T = S - SIGMA
                  DPLUS = D( J ) + T
                  IF(ABS(DPLUS).LT.PIVMIN) 
     $               DPLUS = -PIVMIN
                  TMP = LLD( J ) / DPLUS
                  IF( DPLUS.LT.ZERO ) NEG1=NEG1+1
                  S = T*TMP
                  IF( TMP.EQ.ZERO ) S = LLD( J )
 24            CONTINUE
            ENDIF
         END IF
         NEGCNT = NEGCNT + NEG1
 210  CONTINUE
*
*     II) lower part: L D L^T - SIGMA I = U- D- U-^T
*     
      P = D( N ) - SIGMA
      DO 230 BJ = N-1, R, -BLKLEN
         NEG2 = 0
         XSAV = P
         TO = BJ-BLKLEN+1
         IF ( TO.GE.R ) THEN
            DO 25 J = BJ, TO, -1
               DMINUS = LLD( J ) + P
               IF( DMINUS.LT.ZERO ) NEG2=NEG2+1
               TMP = P / DMINUS
               P = TMP * D( J ) - SIGMA
 25         CONTINUE
         ELSE
            DO 26 J = BJ, R, -1
               DMINUS = LLD( J ) + P
               IF( DMINUS.LT.ZERO ) NEG2=NEG2+1
               TMP = P / DMINUS
               P = TMP * D( J ) - SIGMA
 26         CONTINUE
         ENDIF
         SAWNAN = SISNAN( P )
*
         IF( SAWNAN ) THEN
            NEG2 = 0
            P = XSAV
            TO = BJ-BLKLEN+1
            IF ( TO.GE.R ) THEN
               DO 27 J = BJ, TO, -1
                  DMINUS = LLD( J ) + P
                  IF(ABS(DMINUS).LT.PIVMIN) 
     $               DMINUS = -PIVMIN
                  TMP = D( J ) / DMINUS
                  IF( DMINUS.LT.ZERO ) 
     $               NEG2 = NEG2 + 1
                  P = P*TMP - SIGMA
                  IF( TMP.EQ.ZERO ) 
     $               P = D( J ) - SIGMA
 27            CONTINUE
            ELSE
               DO 28 J = BJ, R, -1
                  DMINUS = LLD( J ) + P
                  IF(ABS(DMINUS).LT.PIVMIN) 
     $               DMINUS = -PIVMIN
                  TMP = D( J ) / DMINUS
                  IF( DMINUS.LT.ZERO ) 
     $               NEG2 = NEG2 + 1
                  P = P*TMP - SIGMA
                  IF( TMP.EQ.ZERO ) 
     $               P = D( J ) - SIGMA
 28            CONTINUE
            ENDIF
         END IF
         NEGCNT = NEGCNT + NEG2
 230  CONTINUE
*     
*     III) Twist index
*
      GAMMA = S + P
      IF( GAMMA.LT.ZERO ) NEGCNT = NEGCNT+1

      SLANEG2 = NEGCNT
      END
*
*
*
      FUNCTION SLANEG2A( N, DLLD, SIGMA, PIVMIN, R )
*
      IMPLICIT NONE
*
      INTEGER SLANEG2A
*
*     .. Scalar Arguments ..
      INTEGER            N, R
      REAL               PIVMIN, SIGMA
*     ..
*     .. Array Arguments ..
      REAL               DLLD( * )
*
      REAL               ZERO
      PARAMETER        ( ZERO = 0.0E0 )

      INTEGER BLKLEN
      PARAMETER ( BLKLEN = 512 )
*
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT
*     ..
*     .. Local Scalars ..
      INTEGER            BJ, I, J, NB, NEG1, NEG2, NEGCNT, NX
      REAL               DMINUS, DPLUS, GAMMA, P, S, T, TMP, XSAV
      LOGICAL SAWNAN
*     ..
*     .. External Functions ..
      LOGICAL SISNAN
      EXTERNAL SISNAN
      
      NEGCNT = 0
*      
*     I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
*     run dstqds block-wise to avoid excessive work when NaNs occur, 
*     first in chunks of size BLKLEN and then the remainder
*
      NB = INT((R-1)/BLKLEN)
      NX = NB*BLKLEN
      S = ZERO      
      DO 210 BJ = 1, NX, BLKLEN
         NEG1 = 0
         XSAV = S
         DO 21 J = BJ, BJ+BLKLEN-1 
            I = 2*J
            T = S - SIGMA
            DPLUS = DLLD( I-1 ) + T
            IF( DPLUS.LT.ZERO ) NEG1=NEG1 + 1
            S = T*DLLD( I ) / DPLUS 
 21      CONTINUE
         SAWNAN = SISNAN( S )
*         
         IF( SAWNAN ) THEN
            NEG1 = 0
            S = XSAV
            DO 23 J = BJ, BJ+BLKLEN-1 
               I = 2*J
               T = S - SIGMA
               DPLUS = DLLD( I-1 ) + T
               IF(ABS(DPLUS).LT.PIVMIN) 
     $            DPLUS = -PIVMIN
               TMP = DLLD( I ) / DPLUS
               IF( DPLUS.LT.ZERO ) 
     $            NEG1 = NEG1 + 1
               S = T*TMP
               IF( TMP.EQ.ZERO ) S = DLLD( I )
 23         CONTINUE
         END IF
         NEGCNT = NEGCNT + NEG1
 210  CONTINUE
*
      NEG1 = 0
      XSAV = S
      DO 22 J = NX+1, R-1
         I = 2*J
         T = S - SIGMA
         DPLUS = DLLD( I-1 ) + T
         IF( DPLUS.LT.ZERO ) NEG1=NEG1 + 1
         S = T*DLLD( I ) / DPLUS 
 22   CONTINUE
      SAWNAN = SISNAN( S )
*         
      IF( SAWNAN ) THEN
         NEG1 = 0
         S = XSAV
         DO 24 J = NX+1, R-1
            I = 2*J
            T = S - SIGMA
            DPLUS = DLLD( I-1 ) + T
            IF(ABS(DPLUS).LT.PIVMIN) 
     $         DPLUS = -PIVMIN
            TMP = DLLD( I ) / DPLUS
            IF( DPLUS.LT.ZERO ) NEG1=NEG1+1
            S = T*TMP
            IF( TMP.EQ.ZERO ) S = DLLD( I )
 24      CONTINUE
      ENDIF
      NEGCNT = NEGCNT + NEG1
*
*     II) lower part: L D L^T - SIGMA I = U- D- U-^T
*     
      NB = INT((N-R)/BLKLEN)
      NX = N-NB*BLKLEN
      P = DLLD( 2*N-1 ) - SIGMA
      DO 230 BJ = N-1, NX, -BLKLEN
         NEG2 = 0
         XSAV = P
         DO 25 J = BJ, BJ-BLKLEN+1, -1
            I = 2*J
            DMINUS = DLLD( I ) + P
            IF( DMINUS.LT.ZERO ) NEG2=NEG2+1
            TMP = P / DMINUS
            P = TMP * DLLD( I-1 ) - SIGMA
 25      CONTINUE
         SAWNAN = SISNAN( P )
*
         IF( SAWNAN ) THEN
            NEG2 = 0
            P = XSAV
            DO 27 J = BJ, BJ-BLKLEN+1, -1
               I = 2*J
               DMINUS = DLLD( I ) + P
               IF(ABS(DMINUS).LT.PIVMIN) 
     $            DMINUS = -PIVMIN
               TMP = DLLD( I-1 ) / DMINUS
               IF( DMINUS.LT.ZERO ) 
     $            NEG2 = NEG2 + 1
               P = P*TMP - SIGMA
               IF( TMP.EQ.ZERO ) 
     $            P = DLLD( I-1 ) - SIGMA
 27         CONTINUE
         END IF
         NEGCNT = NEGCNT + NEG2
 230  CONTINUE

      NEG2 = 0
      XSAV = P
      DO 26 J = NX-1, R, -1
         I = 2*J
         DMINUS = DLLD( I ) + P
         IF( DMINUS.LT.ZERO ) NEG2=NEG2+1
         TMP = P / DMINUS
         P = TMP * DLLD( I-1 ) - SIGMA
 26   CONTINUE
      SAWNAN = SISNAN( P )
*
      IF( SAWNAN ) THEN
         NEG2 = 0
         P = XSAV
         DO 28 J = NX-1, R, -1
            I = 2*J
            DMINUS = DLLD( I ) + P
            IF(ABS(DMINUS).LT.PIVMIN) 
     $         DMINUS = -PIVMIN
            TMP = DLLD( I-1 ) / DMINUS
            IF( DMINUS.LT.ZERO ) 
     $         NEG2 = NEG2 + 1
            P = P*TMP - SIGMA
            IF( TMP.EQ.ZERO ) 
     $         P = DLLD( I-1 ) - SIGMA
 28      CONTINUE
      END IF
      NEGCNT = NEGCNT + NEG2
*     
*     III) Twist index
*
      GAMMA = S + P
      IF( GAMMA.LT.ZERO ) NEGCNT = NEGCNT+1

      SLANEG2A = NEGCNT
      END

