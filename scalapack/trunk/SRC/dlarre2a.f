      SUBROUTINE DLARRE2A( RANGE, N, VL, VU, IL, IU, D, E, E2,
     $                    RTOL1, RTOL2, SPLTOL, NSPLIT, ISPLIT, 
     $                    M, DOL, DOU, NEEDIL, NEEDIU,
     $                    W, WERR, WGAP, IBLOCK, INDEXW, GERS, 
     $                    SDIAM, PIVMIN, WORK, IWORK, MINRGP, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ of Colorado Denver
*     July 4, 2010
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          RANGE
      INTEGER            DOL, DOU, IL, INFO, IU, M, N, NSPLIT,
     $                   NEEDIL, NEEDIU
      DOUBLE PRECISION   MINRGP, PIVMIN, RTOL1, RTOL2, SPLTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * ),
     $                   INDEXW( * )
      DOUBLE PRECISION   D( * ), E( * ), E2( * ), GERS( * ), 
     $                   SDIAM( * ), W( * ),WERR( * ), 
     $                   WGAP( * ), WORK( * )
*
*  Purpose
*  =======
*
*  To find the desired eigenvalues of a given real symmetric
*  tridiagonal matrix T, DLARRE2 sets any "small" off-diagonal
*  elements to zero, and for each unreduced block T_i, it finds
*  (a) a suitable shift at one end of the block's spectrum,
*  (b) the base representation, T_i - sigma_i I = L_i D_i L_i^T, and
*  (c) eigenvalues of each L_i D_i L_i^T.
*
*  NOTE:
*  The algorithm obtains a crude picture of all the wanted eigenvalues
*  (as selected by RANGE). However, to reduce work and improve scalability,
*  only the eigenvalues DOL to DOU are refined. Furthermore, if the matrix 
*  splits into blocks, RRRs for blocks that do not contain eigenvalues
*  from DOL to DOU are skipped.
*  The DQDS algorithm (subroutine DLASQ2) is not used, unlike in
*  the sequential case. Instead, eigenvalues are computed in parallel to some 
*  figures using bisection.

*
*  Arguments
*  =========
*
*  RANGE   (input) CHARACTER
*          = 'A': ("All")   all eigenvalues will be found.
*          = 'V': ("Value") all eigenvalues in the half-open interval
*                           (VL, VU] will be found.
*          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
*                           entire matrix) will be found.
*
*  N       (input) INTEGER
*          The order of the matrix. N > 0.
*
*  VL      (input/output) DOUBLE PRECISION
*  VU      (input/output) DOUBLE PRECISION
*          If RANGE='V', the lower and upper bounds for the eigenvalues.
*          Eigenvalues less than or equal to VL, or greater than VU,
*          will not be returned.  VL < VU.
*          If RANGE='I' or ='A', DLARRE2A computes bounds on the desired 
*          part of the spectrum.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the N diagonal elements of the tridiagonal
*          matrix T.
*          On exit, the N diagonal elements of the diagonal
*          matrices D_i.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the first (N-1) entries contain the subdiagonal
*          elements of the tridiagonal matrix T; E(N) need not be set.
*          On exit, E contains the subdiagonal elements of the unit
*          bidiagonal matrices L_i. The entries E( ISPLIT( I ) ),
*          1 <= I <= NSPLIT, contain the base points sigma_i on output.
*
*  E2      (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the first (N-1) entries contain the SQUARES of the 
*          subdiagonal elements of the tridiagonal matrix T; 
*          E2(N) need not be set.
*          On exit, the entries E2( ISPLIT( I ) ),
*          1 <= I <= NSPLIT, have been set to zero
*
*  RTOL1   (input) DOUBLE PRECISION
*  RTOL2   (input) DOUBLE PRECISION
*           Parameters for bisection.
*           An interval [LEFT,RIGHT] has converged if
*           RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )
*
*  SPLTOL (input) DOUBLE PRECISION
*          The threshold for splitting.
*
*  NSPLIT  (output) INTEGER
*          The number of blocks T splits into. 1 <= NSPLIT <= N.
*
*  ISPLIT  (output) INTEGER array, dimension (N)
*          The splitting points, at which T breaks up into blocks.
*          The first block consists of rows/columns 1 to ISPLIT(1),
*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
*          etc., and the NSPLIT-th consists of rows/columns
*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
*
*  M       (output) INTEGER
*          The total number of eigenvalues (of all L_i D_i L_i^T)
*          found.
*
*  DOL     (input) INTEGER
*  DOU     (input) INTEGER
*          If the user wants to work on only a selected part of the 
*          representation tree, he can specify an index range DOL:DOU.
*          Otherwise, the setting DOL=1, DOU=N should be applied. 
*          Note that DOL and DOU refer to the order in which the eigenvalues 
*          are stored in W. 
*
*  NEEDIL  (output) INTEGER
*  NEEDIU  (output) INTEGER
*          The indices of the leftmost and rightmost eigenvalues
*          of the root node RRR which are
*          needed to accurately compute the relevant part of the 
*          representation tree.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          The first M elements contain the eigenvalues. The
*          eigenvalues of each of the blocks, L_i D_i L_i^T, are
*          sorted in ascending order ( DLARRE2A may use the
*          remaining N-M elements as workspace).
*          Note that immediately after exiting this routine, only 
*          the eigenvalues from position DOL:DOU in W are 
*          reliable on this processor
*          because the eigenvalue computation is done in parallel.
*
*  WERR    (output) DOUBLE PRECISION array, dimension (N)
*          The error bound on the corresponding eigenvalue in W.
*          Note that immediately after exiting this routine, only 
*          the uncertainties from position DOL:DOU in WERR are
*          reliable on this processor
*          because the eigenvalue computation is done in parallel.
*
*  WGAP    (output) DOUBLE PRECISION array, dimension (N)
*          The separation from the right neighbor eigenvalue in W.
*          The gap is only with respect to the eigenvalues of the same block
*          as each block has its own representation tree.
*          Exception: at the right end of a block we store the left gap
*          Note that immediately after exiting this routine, only 
*          the gaps from position DOL:DOU in WGAP are
*          reliable on this processor
*          because the eigenvalue computation is done in parallel.
*
*  IBLOCK  (output) INTEGER array, dimension (N)
*          The indices of the blocks (submatrices) associated with the
*          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue
*          W(i) belongs to the first block from the top, =2 if W(i)
*          belongs to the second block, etc.
*
*  INDEXW  (output) INTEGER array, dimension (N)
*          The indices of the eigenvalues within each block (submatrix);
*          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the
*          i-th eigenvalue W(i) is the 10-th eigenvalue in block 2
*
*  GERS    (output) DOUBLE PRECISION array, dimension (2*N)
*          The N Gerschgorin intervals (the i-th Gerschgorin interval
*          is (GERS(2*i-1), GERS(2*i)).
*
*  PIVMIN  (output) DOUBLE PRECISION
*          The minimum pivot in the sturm sequence for T.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (6*N)
*          Workspace.
*
*  IWORK   (workspace) INTEGER array, dimension (5*N)
*          Workspace.
*
*  MINRGP  (input) DOUBLE PRECISION
*          The minimum relativ gap threshold to decide whether an eigenvalue
*          or a cluster boundary is reached.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          > 0:  A problem occured in DLARRE2A.
*          < 0:  One of the called subroutines signaled an internal problem. 
*                Needs inspection of the corresponding parameter IINFO
*                for further information.
*
*          =-1:  Problem in DLARRD2. 
*          = 2:  No base representation could be found in MAXTRY iterations.
*                Increasing MAXTRY and recompilation might be a remedy.
*          =-3:  Problem in DLARRB2 when computing the refined root 
*                representation
*          =-4:  Problem in DLARRB2 when preforming bisection on the 
*                desired part of the spectrum.
*          = -9  Problem: M < DOU-DOL+1, that is the code found fewer
*                eigenvalues than it was supposed to
*
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   FAC, FOUR, FOURTH, FUDGE, HALF, HNDRD,
     $                   MAXGROWTH, ONE, PERT, TWO, ZERO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, 
     $                     TWO = 2.0D0, FOUR=4.0D0,
     $                     HNDRD = 100.0D0,
     $                     PERT = 8.0D0,
     $                     HALF = ONE/TWO, FOURTH = ONE/FOUR, FAC= HALF,
     $                     MAXGROWTH = 64.0D0, FUDGE = 2.0D0 )
      INTEGER            MAXTRY
      PARAMETER          ( MAXTRY = 6 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOREP, RNDPRT, USEDQD
      INTEGER            CNT, CNT1, CNT2, I, IBEGIN, IDUM, IEND, IINFO,
     $                   IN, INDL, INDU, IRANGE, J, JBLK, MB, MM,
     $                   MYINDL, MYINDU, MYWBEG, MYWEND, WBEGIN, WEND
      DOUBLE PRECISION   AVGAP, BSRTOL, CLWDTH, DMAX, DPIVOT, EABS,
     $                   EMAX, EOLD, EPS, GL, GU, ISLEFT, ISRGHT,
     $                   LGPVMN, LGSPDM, RTL, S1, S2, SAFMIN, SGNDEF,
     $                   SIGMA, SPDIAM, TAU, TMP, TMP1, TMP2


*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION            DLAMCH
      EXTERNAL           DLAMCH, LSAME

*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLARNV, DLARRA, DLARRB2,
     $                   DLARRC, DLARRD2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN

*     ..
*     .. Executable Statements ..
*

      INFO = 0

*     Dis-/Enable a small random perturbation of the root representation
      RNDPRT = .TRUE.
*
*     Decode RANGE
*
      IF( LSAME( RANGE, 'A' ) ) THEN
         IRANGE = 1
      ELSE IF( LSAME( RANGE, 'V' ) ) THEN
         IRANGE = 2
      ELSE IF( LSAME( RANGE, 'I' ) ) THEN
         IRANGE = 3
      END IF

      M = 0

*     Get machine constants
      SAFMIN = DLAMCH( 'S' )
      EPS = DLAMCH( 'P' )

*     Set parameters
      RTL = SQRT(EPS)
      BSRTOL = 1.0D-1

*     Treat case of 1x1 matrix for quick return
      IF( N.EQ.1 ) THEN
         IF( (IRANGE.EQ.1).OR.
     $       ((IRANGE.EQ.2).AND.(D(1).GT.VL).AND.(D(1).LE.VU)).OR.
     $       ((IRANGE.EQ.3).AND.(IL.EQ.1).AND.(IU.EQ.1)) ) THEN
            M = 1
            W(1) = D(1)
*           The computation error of the eigenvalue is zero
            WERR(1) = ZERO
            WGAP(1) = ZERO
            IBLOCK( 1 ) = 1
            INDEXW( 1 ) = 1
            GERS(1) = D( 1 ) 
            GERS(2) = D( 1 ) 
         ENDIF       
*        store the shift for the initial RRR, which is zero in this case 
         E(1) = ZERO
         RETURN
      END IF

*     General case: tridiagonal matrix of order > 1

*     Init WERR, WGAP. 
      DO 1 I =1,N
         WERR(I) = ZERO
 1    CONTINUE
      DO 2 I =1,N
         WGAP(I) = ZERO
 2    CONTINUE

*     Compute Gerschgorin intervals and spectral diameter.
*     Compute maximum off-diagonal entry and pivmin.
      GL = D(1)
      GU = D(1)
      EOLD = ZERO
      EMAX = ZERO
      E(N) = ZERO
      DO 5 I = 1,N
         EABS = ABS( E(I) )
         IF( EABS .GE. EMAX ) THEN
            EMAX = EABS
         END IF
         TMP = EABS + EOLD
         EOLD  = EABS
         TMP1 = D(I) - TMP
         TMP2 = D(I) + TMP
         GL = MIN( GL, TMP1 )
         GU = MAX( GU, TMP2 )
         GERS( 2*I-1) = TMP1
         GERS( 2*I ) = TMP2
 5    CONTINUE
*     The minimum pivot allowed in the sturm sequence for T
      PIVMIN = SAFMIN * MAX( ONE, EMAX**2 )      
*     Compute spectral diameter. The Gerschgorin bounds give an
*     estimate that is wrong by at most a factor of SQRT(2)
      SPDIAM = GU - GL

*     Compute splitting points
      CALL DLARRA( N, D, E, E2, SPLTOL, SPDIAM, 
     $                    NSPLIT, ISPLIT, IINFO )

      IF( IRANGE.EQ.1 ) THEN
*        Set interval [VL,VU] that contains all eigenvalues 
         VL = GL
         VU = GU
      ENDIF

*     We call DLARRD2 to find crude approximations to the eigenvalues
*     in the desired range. In case IRANGE = 3, we also obtain the
*     interval (VL,VU] that contains all the wanted eigenvalues.
*     An interval [LEFT,RIGHT] has converged if
*     RIGHT-LEFT.LT.RTOL*MAX(ABS(LEFT),ABS(RIGHT))
*     DLARRD2 needs a WORK of size 4*N, IWORK of size 3*N
      CALL DLARRD2( RANGE, 'B', N, VL, VU, IL, IU, GERS, 
     $                 BSRTOL, D, E, E2, PIVMIN, NSPLIT, ISPLIT, 
     $                 MM, W, WERR, VL, VU, IBLOCK, INDEXW, 
     $                 WORK, IWORK, DOL, DOU, IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = -1
         RETURN
      ENDIF       
*     Make sure that the entries M+1 to N in W, WERR, IBLOCK, INDEXW are 0
      DO 14 I = MM+1,N
         W( I ) = ZERO
         WERR( I ) = ZERO
         IBLOCK( I ) = 0
         INDEXW( I ) = 0
 14   CONTINUE


***
*     Loop over unreduced blocks
      IBEGIN = 1
      WBEGIN = 1
      DO 170 JBLK = 1, NSPLIT
         IEND = ISPLIT( JBLK )
         IN = IEND - IBEGIN + 1

*        1 X 1 block
         IF( IN.EQ.1 ) THEN
            IF( (IRANGE.EQ.1).OR.( (IRANGE.EQ.2).AND.
     $         ( D( IBEGIN ).GT.VL ).AND.( D( IBEGIN ).LE.VU ) )
     $        .OR. ( (IRANGE.EQ.3).AND.(IBLOCK(WBEGIN).EQ.JBLK))
     $        ) THEN
               M = M + 1
               W( M ) = D( IBEGIN )
               WERR(M) = ZERO
*              The gap for a single block doesn't matter for the later 
*              algorithm and is assigned an arbitrary large value
               WGAP(M) = ZERO
               IBLOCK( M ) = JBLK
               INDEXW( M ) = 1
               WBEGIN = WBEGIN + 1
            ENDIF
*           E( IEND ) holds the shift for the initial RRR
            E( IEND ) = ZERO
            IBEGIN = IEND + 1
            GO TO 170
         END IF
*
*        Blocks of size larger than 1x1
*
*        E( IEND ) will hold the shift for the initial RRR, for now set it =0
         E( IEND ) = ZERO

         IF( ( IRANGE.EQ.1 ) .OR.
     $       ((IRANGE.EQ.3).AND.(IL.EQ.1.AND.IU.EQ.N)) ) THEN
*           MB =  number of eigenvalues to compute
            MB = IN
            WEND = WBEGIN + MB - 1
            INDL = 1
            INDU = IN
         ELSE
*           Count the number of eigenvalues in the current block.
            MB = 0
            DO 20 I = WBEGIN,MM
               IF( IBLOCK(I).EQ.JBLK ) THEN
                  MB = MB+1
               ELSE
                  GOTO 21
               ENDIF 
 20         CONTINUE
 21         CONTINUE

            IF( MB.EQ.0) THEN
*              No eigenvalue in the current block lies in the desired range
*              E( IEND ) holds the shift for the initial RRR
               E( IEND ) = ZERO
               IBEGIN = IEND + 1
               GO TO 170
            ENDIF
*
            WEND = WBEGIN + MB - 1
*           Find local index of the first and last desired evalue.
            INDL = INDEXW(WBEGIN)
            INDU = INDEXW( WEND )
	 ENDIF
*        
         IF( (WEND.LT.DOL).OR.(WBEGIN.GT.DOU) ) THEN
*           if this subblock contains no desired eigenvalues,
*           skip the computation of this representation tree
            IBEGIN = IEND + 1
            WBEGIN = WEND + 1
            M = M + MB
            GO TO 170
         END IF
*
         IF(.NOT. ( IRANGE.EQ.1 ) ) THEN

*           At this point, the sequential code decides
*	    whether dqds or bisection is more efficient.
*           Note: in the parallel code, we do not use dqds.
*           However, we do not change the shift strategy
*           if USEDQD is TRUE, then the same shift is used as for
*           the sequential code when it uses dqds.
*	    
            USEDQD = ( MB .GT. FAC*IN )
*	    
*           Calculate gaps for the current block
*           In later stages, when representations for individual 
*           eigenvalues are different, we use SIGMA = E( IEND ).
            SIGMA = ZERO
            DO 30 I = WBEGIN, WEND - 1
               WGAP( I ) = MAX( ZERO, 
     $                     W(I+1)-WERR(I+1) - (W(I)+WERR(I)) )
 30         CONTINUE
            WGAP( WEND ) = MAX( ZERO, 
     $                  VU - SIGMA - (W( WEND )+WERR( WEND )))
         ENDIF

*
*        Find local outer bounds GL,GU for the block
         GL = D(IBEGIN)
         GU = D(IBEGIN)
         DO 15 I = IBEGIN , IEND
            GL = MIN( GERS( 2*I-1 ), GL )
            GU = MAX( GERS( 2*I ), GU )
 15      CONTINUE
         SPDIAM = GU - GL
*        Save local spectral diameter for later use
         SDIAM(JBLK) = SPDIAM

*        Find approximations to the extremal eigenvalues of the block
         CALL DLARRK( IN, 1, GL, GU, D(IBEGIN),
     $               E2(IBEGIN), PIVMIN, RTL, TMP, TMP1, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = -1
            RETURN
         ENDIF       
         ISLEFT = MAX(GL, TMP - TMP1
     $            - HNDRD * EPS* ABS(TMP - TMP1))
         CALL DLARRK( IN, IN, GL, GU, D(IBEGIN),
     $               E2(IBEGIN), PIVMIN, RTL, TMP, TMP1, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = -1
            RETURN
         ENDIF       
         ISRGHT = MIN(GU, TMP + TMP1
     $                 + HNDRD * EPS * ABS(TMP + TMP1))
         IF( ( IRANGE.EQ.1 ).OR.USEDQD ) THEN
*           Case of DQDS shift
*           Improve the estimate of the spectral diameter
            SPDIAM = ISRGHT - ISLEFT
         ELSE
*           Case of bisection
*           Find approximations to the wanted extremal eigenvalues
            ISLEFT = MAX(GL, W(WBEGIN) - WERR(WBEGIN) 
     $                  - HNDRD * EPS*ABS(W(WBEGIN)- WERR(WBEGIN) ))
            ISRGHT = MIN(GU,W(WEND) + WERR(WEND)
     $                  + HNDRD * EPS * ABS(W(WEND)+ WERR(WEND)))
	 ENDIF


*        Decide whether the base representation for the current block
*        L_JBLK D_JBLK L_JBLK^T = T_JBLK - sigma_JBLK I
*        should be on the left or the right end of the current block.
*        The strategy is to shift to the end which is "more populated"
         IF( IRANGE.EQ.1 ) THEN
*           If all the eigenvalues have to be computed, we use dqd            
            USEDQD = .TRUE.
*           INDL is the local index of the first eigenvalue to compute
            INDL = 1
            INDU = IN
*           MB =  number of eigenvalues to compute
            MB = IN
            WEND = WBEGIN + MB - 1
*           Define 1/4 and 3/4 points of the spectrum
            S1 = ISLEFT + FOURTH * SPDIAM
	    S2 = ISRGHT - FOURTH * SPDIAM
         ELSE        
*           DLARRD2 has computed IBLOCK and INDEXW for each eigenvalue 
*           approximation. 
*           choose sigma
            IF( USEDQD ) THEN
               S1 = ISLEFT + FOURTH * SPDIAM
	       S2 = ISRGHT - FOURTH * SPDIAM
            ELSE
               TMP = MIN(ISRGHT,VU) -  MAX(ISLEFT,VL)
               S1 =  MAX(ISLEFT,VL) + FOURTH * TMP
               S2 =  MIN(ISRGHT,VU) - FOURTH * TMP
            ENDIF
         ENDIF       

*        Compute the negcount at the 1/4 and 3/4 points
         IF(MB.GT.2) THEN
	    CALL DLARRC( 'T', IN, S1, S2, D(IBEGIN), 
     $                    E(IBEGIN), PIVMIN, CNT, CNT1, CNT2, IINFO)
         ENDIF

	 IF(MB.LE.2) THEN
            SIGMA = GL	 
            SGNDEF = ONE
         ELSEIF( CNT1 - INDL .GE. INDU - CNT2 ) THEN
            IF( IRANGE.EQ.1 ) THEN
               SIGMA = MAX(ISLEFT,GL)
            ELSEIF( USEDQD ) THEN
*              use Gerschgorin bound as shift to get pos def matrix
               SIGMA = ISLEFT
            ELSE
*              use approximation of the first desired eigenvalue of the
*              block as shift
               SIGMA = MAX(ISLEFT,VL)
            ENDIF
            SGNDEF = ONE
         ELSE
            IF( IRANGE.EQ.1 ) THEN
               SIGMA = MIN(ISRGHT,GU)
            ELSEIF( USEDQD ) THEN
*              use Gerschgorin bound as shift to get neg def matrix
*              for dqds                  
               SIGMA = ISRGHT
            ELSE
*              use approximation of the first desired eigenvalue of the
*              block as shift
               SIGMA = MIN(ISRGHT,VU)
            ENDIF
            SGNDEF = -ONE
         ENDIF

 
*        An initial SIGMA has been chosen that will be used for computing
*        T - SIGMA I = L D L^T
*        Define the increment TAU of the shift in case the initial shift 
*        needs to be refined to obtain a factorization with not too much 
*        element growth.
         IF( USEDQD ) THEN
            TAU = SPDIAM*EPS*N + TWO*PIVMIN
            TAU = MAX(TAU,EPS*ABS(SIGMA))
         ELSE
            IF(MB.GT.1) THEN        
               CLWDTH = W(WEND) + WERR(WEND) - W(WBEGIN) - WERR(WBEGIN)
               AVGAP = ABS(CLWDTH / DBLE(WEND-WBEGIN))
               IF( SGNDEF.EQ.ONE ) THEN
                  TAU = HALF*MAX(WGAP(WBEGIN),AVGAP)
                  TAU = MAX(TAU,WERR(WBEGIN))
               ELSE
                  TAU = HALF*MAX(WGAP(WEND-1),AVGAP)
                  TAU = MAX(TAU,WERR(WEND))
               ENDIF
	    ELSE
               TAU = WERR(WBEGIN)
	    ENDIF
         ENDIF
*
         DO 80 IDUM = 1, MAXTRY
*           Compute L D L^T factorization of tridiagonal matrix T - sigma I. 
*           Store D in WORK(1:IN), L in WORK(IN+1:2*IN), and reciprocals of 
*           pivots in WORK(2*IN+1:3*IN)
            DPIVOT = D( IBEGIN ) - SIGMA
            WORK( 1 ) = DPIVOT
            DMAX = ABS( WORK(1) )
            J = IBEGIN
            DO 70 I = 1, IN - 1
               WORK( 2*IN+I ) = ONE / WORK( I )
               TMP = E( J )*WORK( 2*IN+I )
               WORK( IN+I ) = TMP
               DPIVOT = ( D( J+1 )-SIGMA ) - TMP*E( J )
               WORK( I+1 ) = DPIVOT
               DMAX = MAX( DMAX, ABS(DPIVOT) )
               J = J + 1
 70         CONTINUE
*           check for element growth
            IF( DMAX .GT. MAXGROWTH*SPDIAM ) THEN
               NOREP = .TRUE.
	    ELSE
               NOREP = .FALSE.
            ENDIF
	    IF(NOREP) THEN
*              Note that in the case of IRANGE=1, we use the Gerschgorin
*              shift which makes the matrix definite. So we should end up
*              here really only in the case of IRANGE = 2,3                
               IF( IDUM.EQ.MAXTRY-1 ) THEN
                  IF( SGNDEF.EQ.ONE ) THEN 
*                    The fudged Gerschgorin shift should succeed
                     SIGMA = 
     $                    GL - FUDGE*SPDIAM*EPS*N - FUDGE*TWO*PIVMIN
                  ELSE
                     SIGMA = 
     $                    GU + FUDGE*SPDIAM*EPS*N + FUDGE*TWO*PIVMIN
                  END IF
               ELSE
                  SIGMA = SIGMA - SGNDEF * TAU 
                  TAU = TWO * TAU
               END IF
            ELSE    
*              an initial RRR is found 
               GO TO 83 
            END IF
 80      CONTINUE
*        if the program reaches this point, no base representation could be 
*        found in MAXTRY iterations.
         INFO = 2
         RETURN

 83      CONTINUE
*        At this point, we have found an initial base representation
*        T - SIGMA I = L D L^T with not too much element growth.
*        Store the shift.
         E( IEND ) = SIGMA
*        Store D and L.         
         CALL DCOPY( IN, WORK, 1, D( IBEGIN ), 1 )
         CALL DCOPY( IN-1, WORK( IN+1 ), 1, E( IBEGIN ), 1 )

	
         IF(RNDPRT .AND. MB.GT.1 ) THEN
*
*           Perturb each entry of the base representation by a small 
*           (but random) relative amount to overcome difficulties with 
*           glued matrices.
*
            DO 122 I = 1, 4
               ISEED( I ) = 1
 122        CONTINUE

            CALL DLARNV(2, ISEED, 2*IN-1, WORK(1))
            DO 125 I = 1,IN-1
               D(IBEGIN+I-1) = D(IBEGIN+I-1)*(ONE+EPS*PERT*WORK(2*I-1))
               E(IBEGIN+I-1) = E(IBEGIN+I-1)*(ONE+EPS*PERT*WORK(2*I))
 125        CONTINUE
            D(IEND) = D(IEND)*(ONE+EPS*PERT*WORK(2*IN-1))
*
         ENDIF
*
*        Compute the required eigenvalues of L D L' by bisection
*        Shift the eigenvalue approximations
*        according to the shift of their representation. 
         DO 134 J=WBEGIN,WEND
            W(J) = W(J) - SIGMA
            WERR(J) = WERR(J) + ABS(W(J)) * EPS
 134     CONTINUE
*        call DLARRB2 to reduce eigenvalue error of the approximations
*        from DLARRD2
         DO 135 I = IBEGIN, IEND-1
            WORK( I ) = D( I ) * E( I )**2
 135     CONTINUE
*        use bisection to find EV from INDL to INDU
         INDL = INDEXW( WBEGIN )
         INDU = INDEXW( WEND )
*
*        Indicate that the current block contains eigenvalues that
*        are potentially needed later.
*
         NEEDIL = MIN(NEEDIL,WBEGIN)
         NEEDIU = MAX(NEEDIU,WEND)
*
*        For the parallel distributed case, only compute
*        those eigenvalues that have to be computed as indicated by DOL, DOU
*
         MYWBEG = MAX(WBEGIN,DOL) 
         MYWEND = MIN(WEND,DOU)
*
         IF(MYWBEG.GT.WBEGIN) THEN
*           This is the leftmost block containing wanted eigenvalues
*           as well as unwanted ones. To save on communication,
*           check if NEEDIL can be increased even further:
*           on the left end, only the eigenvalues of the cluster
*           including MYWBEG are needed
            DO 136 I = WBEGIN, MYWBEG-1
               IF ( WGAP(I).GE.MINRGP*ABS(W(I)) ) THEN
                  NEEDIL = MAX(I+1,NEEDIL)
               ENDIF
 136        CONTINUE
         ENDIF
         IF(MYWEND.LT.WEND) THEN
*           This is the rightmost block containing wanted eigenvalues
*           as well as unwanted ones. To save on communication,
*           Check if NEEDIU can be decreased even further.
            DO 137 I = MYWEND,WEND-1
               IF ( WGAP(I).GE.MINRGP*ABS(W(I)) ) THEN
                  NEEDIU = MIN(I,NEEDIU)
                  GOTO 138
               ENDIF
 137        CONTINUE
 138        CONTINUE
         ENDIF
*
*        Only compute eigenvalues from MYINDL to MYINDU
*        instead of INDL to INDU
*
         MYINDL = INDEXW( MYWBEG )
         MYINDU = INDEXW( MYWEND )
*
         LGPVMN = LOG( PIVMIN )
         LGSPDM = LOG( SPDIAM + PIVMIN )

         CALL DLARRB2(IN, D(IBEGIN), WORK(IBEGIN),
     $               MYINDL, MYINDU, RTOL1, RTOL2, MYINDL-1,
     $               W(MYWBEG), WGAP(MYWBEG), WERR(MYWBEG),
     $               WORK( 2*N+1 ), IWORK, PIVMIN, 
     $               LGPVMN, LGSPDM, IN, IINFO )
         IF( IINFO .NE. 0 ) THEN
            INFO = -4
            RETURN
         END IF
*        DLARRB2 computes all gaps correctly except for the last one
*        Record distance to VU/GU
         WGAP( WEND ) = MAX( ZERO, 
     $           ( VU-SIGMA ) - ( W( WEND ) + WERR( WEND ) ) )
         DO 140 I = INDL, INDU
            M = M + 1
            IBLOCK(M) = JBLK
            INDEXW(M) = I 
 140     CONTINUE
*
*        proceed with next block
         IBEGIN = IEND + 1
         WBEGIN = WEND + 1
 170  CONTINUE
*
      IF (M.LT.DOU-DOL+1) THEN
         INFO = -9
      ENDIF


      RETURN
*     
*     end of DLARRE2A
*
      END
