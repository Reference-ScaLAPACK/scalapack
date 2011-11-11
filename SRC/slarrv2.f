      SUBROUTINE SLARRV2( N, VL, VU, D, L, PIVMIN,
     $                   ISPLIT, M, DOL, DOU, NEEDIL, NEEDIU,
     $                   MINRGP, RTOL1, RTOL2, W, WERR, WGAP,
     $                   IBLOCK, INDEXW, GERS, SDIAM, 
     $                   Z, LDZ, ISUPPZ,
     $                   WORK, IWORK, VSTART, FINISH, 
     $                   MAXCLS, NDEPTH, PARITY, ZOFFSET, INFO )

*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     July 4, 2010
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            DOL, DOU, INFO, LDZ, M, N, MAXCLS,
     $                   NDEPTH, NEEDIL, NEEDIU, PARITY, ZOFFSET
      REAL               MINRGP, PIVMIN, RTOL1, RTOL2, VL, VU
      LOGICAL VSTART, FINISH 
*     ..
*     .. Array Arguments ..
      INTEGER            IBLOCK( * ), INDEXW( * ), ISPLIT( * ),
     $                   ISUPPZ( * ), IWORK( * )
      REAL               D( * ), GERS( * ), L( * ), SDIAM( * ), 
     $                   W( * ), WERR( * ),
     $                   WGAP( * ), WORK( * )
      REAL              Z( LDZ, * )
*
*  Purpose
*  =======
*
*  SLARRV2 computes the eigenvectors of the tridiagonal matrix
*  T = L D L^T given L, D and APPROXIMATIONS to the eigenvalues of L D L^T.
*  The input eigenvalues should have been computed by SLARRE2A
*  or by precious calls to SLARRV2.
*
*  The major difference between the parallel and the sequential construction
*  of the representation tree is that in the parallel case, not all eigenvalues
*  of a given cluster might be computed locally. Other processors might "own"
*  and refine part of an eigenvalue cluster. This is crucial for scalability. 
*  Thus there might be communication necessary before the current level of the 
*  representation tree can be parsed. 
*
*  Please note:
*  1. The calling sequence has two additional INTEGER parameters, 
*     DOL and DOU, that should satisfy M>=DOU>=DOL>=1. 
*     These parameters are only relevant for the case JOBZ = 'V'.
*     SLARRV2  ONLY computes the eigenVECTORS 
*     corresponding to eigenvalues DOL through DOU in W. (That is,
*     instead of computing the eigenvectors belonging to W(1) 
*     through W(M), only the eigenvectors belonging to eigenvalues
*     W(DOL) through W(DOU) are computed. In this case, only the
*     eigenvalues DOL:DOU are guaranteed to be accurately refined
*     to all figures by Rayleigh-Quotient iteration.
*
*  2. The additional arguments VSTART, FINISH, NDEPTH, PARITY, ZOFFSET 
*     are included as a thread-safe implementation equivalent to SAVE variables.
*     These variables store details about the local representation tree which is
*     computed layerwise. For scalability reasons, eigenvalues belonging to the 
*     locally relevant representation tree might be computed on other processors.
*     These need to be communicated before the inspection of the RRRs can proceed
*     on any given layer.           
*     Note that only when the variable FINISH is true, the computation has ended
*     All eigenpairs between DOL and DOU have been computed. M is set = DOU - DOL + 1.
*
*  3. SLARRV2 needs more workspace in Z than the sequential SLARRV. 
*     It is used to store the conformal embedding of the local representation tree.  
* 
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  VL      (input) REAL            
*  VU      (input) REAL            
*          Lower and upper bounds of the interval that contains the desired
*          eigenvalues. VL < VU. Needed to compute gaps on the left or right
*          end of the extremal eigenvalues in the desired RANGE.
*          VU is currently not used but kept as parameter in case needed.
*
*  D       (input/output) REAL             array, dimension (N)
*          On entry, the N diagonal elements of the diagonal matrix D.
*          On exit, D is overwritten.
*
*  L       (input/output) REAL             array, dimension (N)
*          On entry, the (N-1) subdiagonal elements of the unit
*          bidiagonal matrix L are in elements 1 to N-1 of L 
*          (if the matrix is not splitted.) At the end of each block
*          is stored the corresponding shift as given by SLARRE.
*          On exit, L is overwritten.
*
*  PIVMIN  (in) DOUBLE PRECISION
*          The minimum pivot allowed in the sturm sequence.
*
*  ISPLIT  (input) INTEGER array, dimension (N)
*          The splitting points, at which T breaks up into blocks.
*          The first block consists of rows/columns 1 to
*          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
*          through ISPLIT( 2 ), etc.
*
*  M       (input) INTEGER
*          The total number of input eigenvalues.  0 <= M <= N.
*
*  DOL     (input) INTEGER
*  DOU     (input) INTEGER
*          If the user wants to compute only selected eigenvectors from all
*          the eigenvalues supplied, he can specify an index range DOL:DOU.
*          Or else the setting DOL=1, DOU=M should be applied. 
*          Note that DOL and DOU refer to the order in which the eigenvalues 
*          are stored in W. 
*          If the user wants to compute only selected eigenpairs, then
*          the columns DOL-1 to DOU+1 of the eigenvector space Z contain the
*          computed eigenvectors. All other columns of Z are set to zero.
*          If DOL > 1, then Z(:,DOL-1-ZOFFSET) is used.
*          If DOU < M, then Z(:,DOU+1-ZOFFSET) is used.
*
*
*  NEEDIL  (input/output) INTEGER
*  NEEDIU  (input/output) INTEGER
*          Describe which are the left and right outermost eigenvalues 
*          that still need to be included in the computation. These indices
*          indicate whether eigenvalues from other processors are needed to
*          correctly compute the conformally embedded representation tree.
*          When DOL<=NEEDIL<=NEEDIU<=DOU, all required eigenvalues are local
*          to the processor and no communication is required to compute its
*          part of the representation tree.
*
*  MINRGP  (input) REAL            
*          The minimum relativ gap threshold to decide whether an eigenvalue
*          or a cluster boundary is reached.
*
*  RTOL1   (input) REAL            
*  RTOL2   (input) REAL            
*           Parameters for bisection.
*           An interval [LEFT,RIGHT] has converged if
*           RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )
*
*  W       (input/output) REAL             array, dimension (N)
*          The first M elements of W contain the APPROXIMATE eigenvalues for
*          which eigenvectors are to be computed. The eigenvalues
*          should be grouped by split-off block and ordered from
*          smallest to largest within the block. (The output array
*          W from SSTEGR2A is expected here.) Furthermore, they are with
*          respect to the shift of the corresponding root representation
*          for their block. On exit, 
*          W holds those UNshifted eigenvalues
*          for which eigenvectors have already been computed.
*
*  WERR    (input/output) REAL             array, dimension (N)
*          The first M elements contain the semiwidth of the uncertainty
*          interval of the corresponding eigenvalue in W
*
*  WGAP    (input/output) REAL             array, dimension (N)
*          The separation from the right neighbor eigenvalue in W.
*
*  IBLOCK  (input) INTEGER array, dimension (N)
*          The indices of the blocks (submatrices) associated with the
*          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue
*          W(i) belongs to the first block from the top, =2 if W(i)
*          belongs to the second block, etc.
*
*  INDEXW  (input) INTEGER array, dimension (N)
*          The indices of the eigenvalues within each block (submatrix);
*          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the
*          i-th eigenvalue W(i) is the 10-th eigenvalue in the second block.
*
*  GERS    (input) REAL             array, dimension (2*N)
*          The N Gerschgorin intervals (the i-th Gerschgorin interval
*          is (GERS(2*i-1), GERS(2*i)). The Gerschgorin intervals should
*          be computed from the original UNshifted matrix.
*          Currently NOT used but kept as parameter in case it becomes
*          needed in the future.
*
*  SDIAM   (input) REAL             array, dimension (N)
*          The spectral diameters for all unreduced blocks.
*
*  Z       (output) REAL             array, dimension (LDZ, max(1,M) )
*          If INFO = 0, the first M columns of Z contain the
*          orthonormal eigenvectors of the matrix T
*          corresponding to the input eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          In the distributed version, only a subset of columns
*          is accessed, see DOL,DOU and ZOFFSET.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )
*          The support of the eigenvectors in Z, i.e., the indices
*          indicating the nonzero elements in Z. The I-th eigenvector
*          is nonzero only in elements ISUPPZ( 2*I-1 ) through
*          ISUPPZ( 2*I ).
*
*  WORK    (workspace) REAL             array, dimension (12*N)
*
*  IWORK   (workspace) INTEGER array, dimension (7*N)
*
*  VSTART  (input/output) LOGICAL 
*          .TRUE. on initialization, set to .FALSE. afterwards.
*
*  FINISH  (input/output) LOGICAL 
*          A flag that indicates whether all eigenpairs have been computed.
*
*  MAXCLS  (input/output) INTEGER
*          The largest cluster worked on by this processor in the 
*          representation tree.
*
*  NDEPTH  (input/output) INTEGER
*          The current depth of the representation tree. Set to
*          zero on initial pass, changed when the deeper levels of
*          the representation tree are generated. 
*
*  PARITY  (input/output) INTEGER
*          An internal parameter needed for the storage of the
*          clusters on the current level of the representation tree.
*
*  ZOFFSET (input) INTEGER
*          Offset for storing the eigenpairs when Z is distributed
*          in 1D-cyclic fashion.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*
*          > 0:  A problem occured in SLARRV2.
*          < 0:  One of the called subroutines signaled an internal problem. 
*                Needs inspection of the corresponding parameter IINFO
*                for further information.
*
*          =-1:  Problem in SLARRB2 when refining a child's eigenvalues.
*          =-2:  Problem in SLARRF2 when computing the RRR of a child.
*                When a child is inside a tight cluster, it can be difficult
*                to find an RRR. A partial remedy from the user's point of
*                view is to make the parameter MINRGP smaller and recompile.
*                However, as the orthogonality of the computed vectors is 
*                proportional to 1/MINRGP, the user should be aware that 
*                he might be trading in precision when he decreases MINRGP.
*          =-3:  Problem in SLARRB2 when refining a single eigenvalue
*                after the Rayleigh correction was rejected.
*          = 5:  The Rayleigh Quotient Iteration failed to converge to 
*                full accuracy in MAXITR steps.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXITR, USE30, USE31, USE32A, USE32B
      PARAMETER          ( MAXITR = 10, USE30=30, USE31=31, 
     $                     USE32A=3210, USE32B = 3211 )
      REAL               ZERO, ONE, TWO, THREE, FOUR, HALF
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, 
     $                     TWO = 2.0E0, THREE = 3.0E0,
     $                     FOUR = 4.0E0, HALF = 0.5E0)
*     ..
*     .. Local Arrays ..
      INTEGER            SPLACE( 4 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DELREF, ESKIP, NEEDBS, ONLYLC, STP2II, TRYMID,
     $                   TRYRQC, USEDBS, USEDRQ
      INTEGER            I, IBEGIN, IEND, II, IINCLS, IINDC1, IINDC2,
     $                   IINDWK, IINFO, IM, IN, INDEIG, INDLD, INDLLD,
     $                   INDWRK, ISUPMN, ISUPMX, ITER, ITMP1, ITWIST, J,
     $                   JBLK, K, KK, MINIWSIZE, MINWSIZE, MYWFST,
     $                   MYWLST, NCLUS, NEGCNT, NEWCLS, NEWFST, NEWFTT,
     $                   NEWLST, NEWSIZ, OFFSET, OLDCLS, OLDFST, OLDIEN,
     $                   OLDLST, OLDNCL, P, Q, VRTREE, WBEGIN, WEND,
     $                   WINDEX, WINDMN, WINDPL, ZFROM, ZINDEX, ZTO,
     $                   ZUSEDL, ZUSEDU, ZUSEDW
      REAL               AVGAP, BSTRES, BSTW, ENUFGP, EPS, FUDGE, GAP,
     $                   GAPTOL, LAMBDA, LEFT, LGAP, LGPVMN, LGSPDM,
     $                   LOG_IN, MGAP, MINGMA, MYERR, NRMINV, NXTERR,
     $                   ORTOL, RESID, RGAP, RIGHT, RLTL30, RQCORR,
     $                   RQTOL, SAVEGP, SGNDEF, SIGMA, SPDIAM, SSIGMA,
     $                   TAU, TMP, TOL, ZTZ
*     ..
*     .. External Functions ..
      REAL              SLAMCH
      REAL               SDOT, SNRM2
      EXTERNAL           SDOT, SLAMCH, SNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY, SLAR1VA, SLARRB2,
     $                   SLARRF2, SLASET, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS, REAL, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*     ..


      INFO = 0
*     The first N entries of WORK are reserved for the eigenvalues
      INDLD = N+1
      INDLLD= 2*N+1
      INDWRK= 3*N+1
      MINWSIZE = 12 * N

*     IWORK(IINCLS+JBLK) holds the number of clusters on the current level 
*     of the reptree for block JBLK  
      IINCLS = 0
*     IWORK(IINDC1+1:IINC2+N) are used to store the clusters of the current
*     layer and the one above.
      IINDC1 = N
      IINDC2 = 2*N
      IINDWK = 3*N + 1
      MINIWSIZE = 7 * N

      EPS = SLAMCH( 'Precision' )
      RQTOL = TWO * EPS

      TRYRQC = .TRUE.
*     Decide which representation tree criterion to use
*     USE30 = Lapack 3.0 criterion
*     USE31 = LAPACK 3.1 criterion
*     USE32A = two criteria, determines singletons with USE31, and groups with avgap.
*     USE32B = two criteria, determines singletons with USE31, and groups with USE30.
      VRTREE = USE32A
*
      LGPVMN = LOG( PIVMIN )


      IF(VSTART) THEN
*      
*        PREPROCESSING, DONE ONLY IN THE FIRST CALL
*
         VSTART = .FALSE.   
*
         MAXCLS = 1

*        Set delayed eigenvalue refinement
*        In order to enable more parallelism, refinement
*        must be done immediately and cannot be delayed until
*        the next representation tree level.
         DELREF = .FALSE.

         DO 1 I= 1,MINWSIZE
            WORK( I ) = ZERO 
 1       CONTINUE

         DO 2 I= 1,MINIWSIZE
            IWORK( I ) = 0
 2       CONTINUE

         ZUSEDL = 1
         IF(DOL.GT.1) THEN
*           Set lower bound for use of Z
            ZUSEDL = DOL-1
         ENDIF
         ZUSEDU = M
         IF(DOU.LT.M) THEN
*           Set lower bound for use of Z
            ZUSEDU = DOU+1
         ENDIF
*        The width of the part of Z that is used
         ZUSEDW = ZUSEDU - ZUSEDL + 1
*
         CALL SLASET( 'Full', N, ZUSEDW, ZERO, ZERO, 
     $                    Z(1,(ZUSEDL-ZOFFSET)), LDZ )

*        Initialize NDEPTH, the current depth of the representation tree
         NDEPTH = 0
*        Initialize parity 
         PARITY = 1

*        Go through blocks, initialize data structures
         IBEGIN = 1
         WBEGIN = 1
         DO 10 JBLK = 1, IBLOCK( M )
            IEND = ISPLIT( JBLK )
            SIGMA = L( IEND )
            WEND = WBEGIN - 1
 3          CONTINUE
            IF( WEND.LT.M ) THEN
               IF( IBLOCK( WEND+1 ).EQ.JBLK ) THEN
                  WEND = WEND + 1
                  GO TO 3
               END IF
            END IF
            IF( WEND.LT.WBEGIN ) THEN
               IWORK( IINCLS + JBLK ) = 0
               IBEGIN = IEND + 1
               GO TO 10
            ELSEIF( (WEND.LT.DOL).OR.(WBEGIN.GT.DOU) ) THEN
               IWORK( IINCLS + JBLK ) = 0
               IBEGIN = IEND + 1
               WBEGIN = WEND + 1
               GO TO 10
            END IF
*           The number of eigenvalues in the current block
            IM = WEND - WBEGIN + 1
*           This is for a 1x1 block
            IF( IBEGIN.EQ.IEND ) THEN
               IWORK( IINCLS + JBLK ) = 0
               Z( IBEGIN, (WBEGIN-ZOFFSET) ) = ONE
               ISUPPZ( 2*WBEGIN-1 ) = IBEGIN
               ISUPPZ( 2*WBEGIN ) = IBEGIN
               W( WBEGIN ) = W( WBEGIN ) + SIGMA
               WORK( WBEGIN ) = W( WBEGIN )
               IBEGIN = IEND + 1
               WBEGIN = WBEGIN + 1
               GO TO 10
            END IF
            CALL SCOPY( IM, W( WBEGIN ), 1, 
     &                WORK( WBEGIN ), 1 )	 
*           We store in W the eigenvalue approximations w.r.t. the original
*           matrix T.
            DO 5 I=1,IM
               W(WBEGIN+I-1) = W(WBEGIN+I-1)+SIGMA
 5          CONTINUE

*           Initialize cluster counter for this block
            IWORK( IINCLS + JBLK ) = 1
            IWORK( IINDC1+IBEGIN ) = 1
            IWORK( IINDC1+IBEGIN+1 ) = IM
*
            IBEGIN = IEND + 1
            WBEGIN = WEND + 1
10       CONTINUE
*
      ENDIF 

*     Init NEEDIL and NEEDIU
      NEEDIL = DOU
      NEEDIU = DOL      

*     Here starts the main loop
*     Only one pass through the loop is done until no collaboration
*     with other processors is needed. 
 40   CONTINUE

      PARITY = 1 - PARITY

*     For each block, build next level of representation tree
*     if there are still remaining clusters 
      IBEGIN = 1
      WBEGIN = 1
      DO 170 JBLK = 1, IBLOCK( M )
         IEND = ISPLIT( JBLK )
         SIGMA = L( IEND )
*        Find the eigenvectors of the submatrix indexed IBEGIN
*        through IEND.
         IF(M.EQ.N) THEN
*           all eigenpairs are computed
            WEND = IEND
         ELSE
*           count how many wanted eigenpairs are in this block
            WEND = WBEGIN - 1
 15         CONTINUE
            IF( WEND.LT.M ) THEN
               IF( IBLOCK( WEND+1 ).EQ.JBLK ) THEN
                  WEND = WEND + 1
                  GO TO 15
               END IF
            END IF
         ENDIF

         OLDNCL = IWORK( IINCLS + JBLK )
         IF( OLDNCL.EQ.0 ) THEN
            IBEGIN = IEND + 1
            WBEGIN = WEND + 1
            GO TO 170
         END IF
*        OLDIEN is the last index of the previous block
         OLDIEN = IBEGIN - 1
*        Calculate the size of the current block
         IN = IEND - IBEGIN + 1
*        The number of eigenvalues in the current block
         IM = WEND - WBEGIN + 1

*        Find local spectral diameter of the block
         SPDIAM = SDIAM(JBLK)
         LGSPDM = LOG( SPDIAM + PIVMIN )
*        Compute ORTOL parameter, similar to SSTEIN
         ORTOL = SPDIAM*1.0E-3
*        Compute average gap
         AVGAP = SPDIAM/REAL(IN-1)
*        Compute the minimum of average gap and ORTOL parameter 
*        This can used as a lower bound for acceptable separation 
*        between eigenvalues 
         ENUFGP = MIN(ORTOL,AVGAP)

*        Any 1x1 block has been treated before

*        loop while( OLDNCLS.GT.0 )
*        generate the next representation tree level for the current block
         IF( OLDNCL.GT.0 ) THEN
*           This is a crude protection against infinitely deep trees
            IF( NDEPTH.GT.M ) THEN
               INFO = -2
               RETURN
            ENDIF
*           breadth first processing of the current level of the representation
*           tree: OLDNCL = number of clusters on current level
*           NCLUS is the number of clusters for the next level of the reptree
*           reset NCLUS to count the number of child clusters 
            NCLUS = 0
*
            LOG_IN = LOG(REAL(IN))
*
            RLTL30 = MIN( 1.0E-2, ONE / REAL( IN ) )
*
            IF( PARITY.EQ.0 ) THEN
               OLDCLS = IINDC1+IBEGIN-1
               NEWCLS = IINDC2+IBEGIN-1
            ELSE
               OLDCLS = IINDC2+IBEGIN-1
               NEWCLS = IINDC1+IBEGIN-1
            END IF
*           Process the clusters on the current level
            DO 150 I = 1, OLDNCL
               J = OLDCLS + 2*I
*              OLDFST, OLDLST = first, last index of current cluster.
*                               cluster indices start with 1 and are relative
*                               to WBEGIN when accessing W, WGAP, WERR, Z
               OLDFST = IWORK( J-1 )
               OLDLST = IWORK( J )
               IF( NDEPTH.GT.0 ) THEN
*                 Retrieve relatively robust representation (RRR) of cluster
*                 that has been computed at the previous level
*                 The RRR is stored in Z and overwritten once the eigenvectors
*                 have been computed or when the cluster is refined 

                  IF((DOL.EQ.1).AND.(DOU.EQ.M)) THEN
*                    Get representation from location of the leftmost evalue
*                    of the cluster
                     J = WBEGIN + OLDFST - 1
                  ELSE
                     IF(WBEGIN+OLDFST-1.LT.DOL) THEN
*                       Get representation from the left end of Z array 
                        J = DOL - 1
                     ELSEIF(WBEGIN+OLDFST-1.GT.DOU) THEN
*                       Get representation from the right end of Z array 
                        J = DOU
                     ELSE
                        J = WBEGIN + OLDFST - 1
                     ENDIF
                  ENDIF
                  CALL SCOPY( IN, Z( IBEGIN, (J-ZOFFSET) ), 
     $               1, D( IBEGIN ), 1 )
                  CALL SCOPY( IN-1, Z( IBEGIN, (J+1-ZOFFSET) ), 
     $               1, L( IBEGIN ),1 )
                  SIGMA = Z( IEND, (J+1-ZOFFSET) )
*                 Set the corresponding entries in Z to zero
                  CALL SLASET( 'Full', IN, 2, ZERO, ZERO,
     $                         Z( IBEGIN, (J-ZOFFSET) ), LDZ )
               END IF

*              Compute DL and DLL of current RRR
               DO 50 J = IBEGIN, IEND-1
                  TMP = D( J )*L( J )
                  WORK( INDLD-1+J ) = TMP
                  WORK( INDLLD-1+J ) = TMP*L( J )
   50          CONTINUE

               IF( NDEPTH.GT.0 .AND. DELREF ) THEN
*                 P and Q are index of the first and last eigenvalue to compute
*                 within the current block
                  P = INDEXW( WBEGIN-1+OLDFST )
                  Q = INDEXW( WBEGIN-1+OLDLST )
*                 Offset for the arrays WORK, WGAP and WERR, i.e., th P-OFFSET
*                 thru' Q-OFFSET elements of these arrays are to be used.
C                  OFFSET = P-OLDFST
                  OFFSET = INDEXW( WBEGIN ) - 1
*                 perform limited bisection (if necessary) to get approximate 
*                 eigenvalues to the precision needed.
                  CALL SLARRB2( IN, D( IBEGIN ), 
     $                         WORK(INDLLD+IBEGIN-1),
     $                         P, Q, RTOL1, RTOL2, OFFSET, 
     $                         WORK(WBEGIN),WGAP(WBEGIN),WERR(WBEGIN),
     $                         WORK( INDWRK ), IWORK( IINDWK ),
     $                         PIVMIN, LGPVMN, LGSPDM, IN, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     INFO = -1
                     RETURN
                  ENDIF       
*                 We also recompute the extremal gaps. W holds all eigenvalues
*                 of the unshifted matrix and must be used for computation
*                 of WGAP, the entries of WORK might stem from RRRs with 
*                 different shifts. The gaps from WBEGIN-1+OLDFST to
*                 WBEGIN-1+OLDLST are correctly computed in SLARRB2.
*                 However, we only allow the gaps to become greater since 
*                 this is what should happen when we decrease WERR
                  IF( OLDFST.GT.1) THEN
                     WGAP( WBEGIN+OLDFST-2 ) = 
     $             MAX(WGAP(WBEGIN+OLDFST-2),
     $                 W(WBEGIN+OLDFST-1)-WERR(WBEGIN+OLDFST-1) 
     $                 - W(WBEGIN+OLDFST-2)-WERR(WBEGIN+OLDFST-2) )
                  ENDIF
                  IF( WBEGIN + OLDLST -1 .LT. WEND ) THEN
                     WGAP( WBEGIN+OLDLST-1 ) = 
     $               MAX(WGAP(WBEGIN+OLDLST-1), 
     $                   W(WBEGIN+OLDLST)-WERR(WBEGIN+OLDLST) 
     $                   - W(WBEGIN+OLDLST-1)-WERR(WBEGIN+OLDLST-1) )
                  ENDIF
*                 Each time the eigenvalues in WORK get refined, we store
*                 the newly found approximation with all shifts applied in W
                  DO 53 J=OLDFST,OLDLST
                     W(WBEGIN+J-1) = WORK(WBEGIN+J-1)+SIGMA
 53               CONTINUE
               ELSEIF( (NDEPTH.EQ.0) .OR. (.NOT.DELREF) ) THEN 
*                 Some of the eigenvalues might have been computed on
*                 other processors                  
*                 Recompute gaps for this cluster 
*                 (all eigenvalues have the same
*                 representation, i.e. the same shift, so this is easy)
                  DO 54 J = OLDFST, OLDLST-1
                     MYERR = WERR(WBEGIN + J - 1) 
                     NXTERR = WERR(WBEGIN + J )
                     WGAP(WBEGIN+J-1) = MAX(WGAP(WBEGIN+J-1),
     $                    (   WORK(WBEGIN+J) - NXTERR ) 
     $                  - ( WORK(WBEGIN+J-1) + MYERR )
     $                                     )
 54               CONTINUE
               END IF
*
*              Process the current node.
*
               NEWFST = OLDFST
               DO 140 J = OLDFST, OLDLST
                  IF( J.EQ.OLDLST ) THEN
*                    we are at the right end of the cluster, this is also the
*                    boundary of the child cluster                    
                     NEWLST = J
                  ELSE 
                     IF (VRTREE.EQ.USE30) THEN
                        IF(WGAP( WBEGIN + J -1).GE.
     $                     RLTL30 * ABS(WORK(WBEGIN + J -1)) ) THEN
*                          the right relgap is big enough by the Lapack 3.0 criterion
                           NEWLST = J
                        ELSE
*                          inside a child cluster, the relative gap is not
*                          big enough.
                           GOTO 140
			ENDIF
                     ELSE IF (VRTREE.EQ.USE31) THEN
                        IF ( WGAP( WBEGIN + J -1).GE.
     $                      MINRGP* ABS( WORK(WBEGIN + J -1) ) ) THEN
*                          the right relgap is big enough by the Lapack 3.1 criterion
*                          (NEWFST,..,NEWLST) is well separated from the following 
                           NEWLST = J
                        ELSE
*                          inside a child cluster, the relative gap is not
*                          big enough.
                           GOTO 140
			ENDIF
                     ELSE IF (VRTREE.EQ.USE32A) THEN
                        IF( (J.EQ.OLDFST).AND.( WGAP(WBEGIN+J-1).GE.
     $                      MINRGP* ABS(WORK(WBEGIN+J-1)) ) ) THEN
*                          the right relgap is big enough by the Lapack 3.1 criterion
*                          Found a singleton
                           NEWLST = J
                        ELSE IF( (J.GT.OLDFST).AND.(J.EQ.NEWFST).AND.
     $                           ( WGAP(WBEGIN+J-2).GE.
     $                             MINRGP* ABS(WORK(WBEGIN+J-1)) ).AND. 
     $                           ( WGAP(WBEGIN+J-1).GE.
     $                             MINRGP* ABS(WORK(WBEGIN+J-1)) ) 
     $                     ) THEN
*                          Found a singleton
                           NEWLST = J
                        ELSE IF( (J.GT.NEWFST).AND.WGAP(WBEGIN+J-1).GE.
     $                     (MINRGP*ABS(WORK(WBEGIN+J-1)) ) ) 
     $                     THEN
*                          the right relgap is big enough by the Lapack 3.1 criterion
                           NEWLST = J
                        ELSE IF((J.GT.NEWFST).AND.(J+1.LT.OLDLST).AND.
     $                     (WGAP(WBEGIN+J-1).GE.ENUFGP))
     $                     THEN
*                          the right gap is bigger than ENUFGP
*                          Care needs to be taken with this criterion to make
*                          sure it does not create a remaining `false' singleton
                           NEWLST = J
                        ELSE
*                          inside a child cluster, the relative gap is not
*                          big enough.
                           GOTO 140
			ENDIF
                     ELSE IF (VRTREE.EQ.USE32B) THEN
                        IF( (J.EQ.OLDFST).AND.( WGAP(WBEGIN+J-1).GE.
     $                      MINRGP* ABS(WORK(WBEGIN+J-1)) ) ) THEN
*                          the right relgap is big enough by the Lapack 3.1 criterion
*                          Found a singleton
                           NEWLST = J
                        ELSE IF( (J.GT.OLDFST).AND.(J.EQ.NEWFST).AND.
     $                           ( WGAP(WBEGIN+J-2).GE.
     $                             MINRGP* ABS(WORK(WBEGIN+J-1)) ).AND. 
     $                           ( WGAP(WBEGIN+J-1).GE.
     $                             MINRGP* ABS(WORK(WBEGIN+J-1)) ) 
     $                     ) THEN
*                          Found a singleton
                           NEWLST = J
                        ELSE IF( (J.GT.NEWFST).AND.WGAP(WBEGIN+J-1).GE.
     $                     (MINRGP*ABS(WORK(WBEGIN+J-1)) ) ) 
     $                     THEN
*                          the right relgap is big enough by the Lapack 3.1 criterion
                           NEWLST = J
                        ELSE IF((J.GT.NEWFST).AND.(J+1.LT.OLDLST).AND.
     $                     (WGAP( WBEGIN + J -1).GE.
     $                     RLTL30 * ABS(WORK(WBEGIN + J -1)) ))
     $                     THEN
*                          the right relgap is big enough by the Lapack 3.0 criterion
*                          Care needs to be taken with this criterion to make
*                          sure it does not create a remaining `false' singleton
                           NEWLST = J
                        ELSE
*                          inside a child cluster, the relative gap is not
*                          big enough.
                           GOTO 140
			ENDIF
                     END IF
                  END IF

*                 Compute size of child cluster found
                  NEWSIZ = NEWLST - NEWFST + 1
                  MAXCLS = MAX( NEWSIZ, MAXCLS )

*                 NEWFTT is the place in Z where the new RRR or the computed
*                 eigenvector is to be stored
                  IF((DOL.EQ.1).AND.(DOU.EQ.M)) THEN
*                    Store representation at location of the leftmost evalue
*                    of the cluster
                     NEWFTT = WBEGIN + NEWFST - 1
                  ELSE
                     IF(WBEGIN+NEWFST-1.LT.DOL) THEN
*                       Store representation at the left end of Z array 
                        NEWFTT = DOL - 1
                     ELSEIF(WBEGIN+NEWFST-1.GT.DOU) THEN
*                       Store representation at the right end of Z array 
                        NEWFTT = DOU
                     ELSE
                        NEWFTT = WBEGIN + NEWFST - 1
                     ENDIF
                  ENDIF
*                 FOR 1D-DISTRIBUTED Z, COMPUTE NEWFTT SHIFTED BY ZOFFSET
                  NEWFTT = NEWFTT - ZOFFSET

                  IF( NEWSIZ.GT.1) THEN
*
*                    Current child is not a singleton but a cluster.
*
*
                     IF((WBEGIN+NEWLST-1.LT.DOL).OR.
     $                  (WBEGIN+NEWFST-1.GT.DOU)) THEN
*                       if the cluster contains no desired eigenvalues
*                       skip the computation of that branch of the rep. tree
                        GOTO 139
                     ENDIF

*                    Compute left and right cluster gap.
*
                     IF( NEWFST.EQ.1 ) THEN
                        LGAP = MAX( ZERO, 
     $                       W(WBEGIN)-WERR(WBEGIN) - VL )
                     ELSE
                        LGAP = WGAP( WBEGIN+NEWFST-2 )
                     ENDIF
                     RGAP = WGAP( WBEGIN+NEWLST-1 )
*
*                    For larger clusters, record the largest gap observed 
*                    somewhere near the middle of the cluster as a possible 
*                    alternative position for a shift when TRYMID is TRUE
*		     
                     MGAP = ZERO
                     IF(NEWSIZ.GE.50) THEN
                        KK = NEWFST
                        DO 545 K =NEWFST+NEWSIZ/3,NEWLST-NEWSIZ/3
		           IF(MGAP.LT.WGAP( WBEGIN+K-1 )) THEN
		              KK = K
		              MGAP = WGAP( WBEGIN+K-1 )
                           ENDIF
 545	                CONTINUE
                     ENDIF
		     
*
*                    Record the left- and right-most eigenvalues needed
*                    for the next level of the representation tree
                     NEEDIL = MIN(NEEDIL,WBEGIN+NEWFST-1)
                     NEEDIU = MAX(NEEDIU,WBEGIN+NEWLST-1)

*
*                    Check if middle gap is large enough to shift there
*
                     GAP = MIN(LGAP,RGAP)
		     TRYMID = (MGAP.GT.GAP)

		     SPLACE(1) = NEWFST
		     SPLACE(2) = NEWLST
		     IF(TRYMID) THEN
		        SPLACE(3) = KK
                        SPLACE(4) = KK+1
		     ELSE
		        SPLACE(3) = NEWFST
		        SPLACE(4) = NEWLST
		     ENDIF
*
*                    Compute left- and rightmost eigenvalue of child
*                    to high precision in order to shift as close
*                    as possible and obtain as large relative gaps
*                    as possible
*

                     DO 55 K =1,4
                        P = INDEXW( WBEGIN-1+SPLACE(K) )
                        OFFSET = INDEXW( WBEGIN ) - 1
                        CALL SLARRB2( IN, D(IBEGIN), 
     $                       WORK( INDLLD+IBEGIN-1 ),P,P,
     $                       RQTOL, RQTOL, OFFSET, 
     $                       WORK(WBEGIN),WGAP(WBEGIN),
     $                       WERR(WBEGIN),WORK( INDWRK ), 
     $                       IWORK( IINDWK ), 
     $                       PIVMIN, LGPVMN, LGSPDM, IN, IINFO )
 55                  CONTINUE
*
*                    Compute RRR of child cluster.
*                    Note that the new RRR is stored in Z                     
*
C                    SLARRF2 needs LWORK = 2*N
                     CALL SLARRF2( IN, D( IBEGIN ), L( IBEGIN ),
     $                         WORK(INDLD+IBEGIN-1), 
     $                         SPLACE(1), SPLACE(2), 
     $                         SPLACE(3), SPLACE(4), WORK(WBEGIN),
     $                         WGAP(WBEGIN), WERR(WBEGIN), TRYMID,
     $                         SPDIAM, LGAP, RGAP, PIVMIN, TAU, 
     $                         Z( IBEGIN, NEWFTT ),
     $                         Z( IBEGIN, NEWFTT+1 ),
     $                         WORK( INDWRK ), IINFO )
                     IF( IINFO.EQ.0 ) THEN
*                       a new RRR for the cluster was found by SLARRF2
*                       update shift and store it         
                        SSIGMA = SIGMA + TAU
                        Z( IEND, NEWFTT+1 ) = SSIGMA
*                       WORK() are the midpoints and WERR() the semi-width
*                       Note that the entries in W are unchanged.
                        DO 116 K = NEWFST, NEWLST
                           FUDGE = 
     $                          THREE*EPS*ABS(WORK(WBEGIN+K-1))
                           WORK( WBEGIN + K - 1 ) = 
     $                          WORK( WBEGIN + K - 1) - TAU
                           FUDGE = FUDGE + 
     $                          FOUR*EPS*ABS(WORK(WBEGIN+K-1))
*                          Fudge errors
                           WERR( WBEGIN + K - 1 ) =
     $                          WERR( WBEGIN + K - 1 ) + FUDGE
 116                    CONTINUE

                        NCLUS = NCLUS + 1
                        K = NEWCLS + 2*NCLUS
                        IWORK( K-1 ) = NEWFST
                        IWORK( K ) = NEWLST
*
                        IF(.NOT.DELREF) THEN
                           ONLYLC = .TRUE.
*
                           IF(ONLYLC) THEN
                              MYWFST = MAX(WBEGIN-1+NEWFST,DOL-1)
                              MYWLST = MIN(WBEGIN-1+NEWLST,DOU+1)
                           ELSE
                              MYWFST = WBEGIN-1+NEWFST
                              MYWLST = WBEGIN-1+NEWLST 
                           ENDIF

*                          Compute LLD of new RRR
                           DO 5000 K = IBEGIN, IEND-1
                              WORK( INDWRK-1+K ) = 
     $                        Z(K,NEWFTT)*
     $                       (Z(K,NEWFTT+1)**2)
 5000                      CONTINUE
*                          P and Q are index of the first and last 
*                          eigenvalue to compute within the new cluster
                           P = INDEXW( MYWFST )
                           Q = INDEXW( MYWLST )
*                          Offset for the arrays WORK, WGAP and WERR
                           OFFSET = INDEXW( WBEGIN ) - 1
*                          perform limited bisection (if necessary) to get approximate 
*                          eigenvalues to the precision needed.
                           CALL SLARRB2( IN, 
     $                         Z(IBEGIN, NEWFTT ),
     $                         WORK(INDWRK+IBEGIN-1),
     $                         P, Q, RTOL1, RTOL2, OFFSET, 
     $                         WORK(WBEGIN),WGAP(WBEGIN),WERR(WBEGIN),
     $                         WORK( INDWRK+N ), IWORK( IINDWK ),
     $                         PIVMIN, LGPVMN, LGSPDM, IN, IINFO )
                           IF( IINFO.NE.0 ) THEN
                              INFO = -1
                              RETURN
                           ENDIF       
*                          Each time the eigenvalues in WORK get refined, we store
*                          the newly found approximation with all shifts applied in W
                           DO 5003 K=NEWFST,NEWLST
                              W(WBEGIN+K-1) = WORK(WBEGIN+K-1)+SSIGMA
 5003                      CONTINUE
                        ENDIF
*
                     ELSE    
                        INFO = -2
                        RETURN
                     ENDIF      
	          ELSE
*
*                    Compute eigenvector of singleton
*
                     ITER = 0
*                    
                     TOL = FOUR * LOG_IN * EPS
*
                     K = NEWFST
                     WINDEX = WBEGIN + K - 1
                     ZINDEX = WINDEX - ZOFFSET
                     WINDMN = MAX(WINDEX - 1,1)
                     WINDPL = MIN(WINDEX + 1,M)
                     LAMBDA = WORK( WINDEX )
*                    Check if eigenvector computation is to be skipped
                     IF((WINDEX.LT.DOL).OR.
     $                  (WINDEX.GT.DOU)) THEN
                        ESKIP = .TRUE.
                        GOTO 125
                     ELSE
                        ESKIP = .FALSE.
                     ENDIF
                     LEFT = WORK( WINDEX ) - WERR( WINDEX )
                     RIGHT = WORK( WINDEX ) + WERR( WINDEX )
                     INDEIG = INDEXW( WINDEX )
                     IF( K .EQ. 1) THEN
                        LGAP = EPS*MAX(ABS(LEFT),ABS(RIGHT))
                     ELSE
                        LGAP = WGAP(WINDMN)
                     ENDIF
                     IF( K .EQ. IM) THEN
                        RGAP = EPS*MAX(ABS(LEFT),ABS(RIGHT))
                     ELSE
                        RGAP = WGAP(WINDEX)
                     ENDIF
                     GAP = MIN( LGAP, RGAP )
                     IF(( K .EQ. 1).OR.(K .EQ. IM)) THEN
                        GAPTOL = ZERO
                     ELSE
                        GAPTOL = GAP * EPS
                     ENDIF
                     ISUPMN = IN
                     ISUPMX = 1
*                    Update WGAP so that it holds the minimum gap 
*                    to the left or the right. This is crucial in the
*                    case where bisection is used to ensure that the
*                    eigenvalue is refined up to the required precision.
*                    The correct value is restored afterwards.
                     SAVEGP = WGAP(WINDEX)
                     WGAP(WINDEX) = GAP
*                    We want to use the Rayleigh Quotient Correction
*                    as often as possible since it converges quadratically
*                    when we are close enough to the desired eigenvalue.
*                    However, the Rayleigh Quotient can have the wrong sign
*                    and lead us away from the desired eigenvalue. In this
*                    case, the best we can do is to use bisection.
                     USEDBS = .FALSE.
                     USEDRQ = .FALSE.
*                    Bisection is initially turned off unless it is forced
                     NEEDBS =  .NOT.TRYRQC 
*                    Reset ITWIST
                     ITWIST = 0
 120                 CONTINUE
*                    Check if bisection should be used to refine eigenvalue
                     IF(NEEDBS) THEN
*                       Take the bisection as new iterate
                        USEDBS = .TRUE.
*                       Temporary copy of twist index needed
                        ITMP1 = ITWIST
                        OFFSET = INDEXW( WBEGIN ) - 1
                        CALL SLARRB2( IN, D(IBEGIN), 
     $                       WORK(INDLLD+IBEGIN-1),INDEIG,INDEIG,
     $                       ZERO, TWO*EPS, OFFSET, 
     $                       WORK(WBEGIN),WGAP(WBEGIN),
     $                       WERR(WBEGIN),WORK( INDWRK ), 
     $                       IWORK( IINDWK ), 
     $                       PIVMIN, LGPVMN, LGSPDM, ITMP1, IINFO )
                        IF( IINFO.NE.0 ) THEN
                           INFO = -3
                           RETURN
                        ENDIF       
                        LAMBDA = WORK( WINDEX )
*                       Reset twist index from inaccurate LAMBDA to
*                       force computation of true MINGMA  
                        ITWIST = 0
                     ENDIF
*                    Given LAMBDA, compute the eigenvector.
                     CALL SLAR1VA( IN, 1, IN, LAMBDA, D(IBEGIN),
     $                    L( IBEGIN ), WORK(INDLD+IBEGIN-1), 
     $                    WORK(INDLLD+IBEGIN-1),
     $                    PIVMIN, GAPTOL, Z( IBEGIN, ZINDEX),
     $                    .NOT.USEDBS, NEGCNT, ZTZ, MINGMA,
     $                    ITWIST, ISUPPZ( 2*WINDEX-1 ),
     $                    NRMINV, RESID, RQCORR, WORK( INDWRK ) )
                     IF(ITER .EQ. 0) THEN
                        BSTRES = RESID
                        BSTW = LAMBDA
                     ELSEIF(RESID.LT.BSTRES) THEN
                        BSTRES = RESID
                        BSTW = LAMBDA
                     ENDIF
                     ISUPMN = MIN(ISUPMN,ISUPPZ( 2*WINDEX-1 ))
                     ISUPMX = MAX(ISUPMX,ISUPPZ( 2*WINDEX ))
                     ITER = ITER + 1
*		     
*                    Convergence test for Rayleigh-Quotient iteration
*                    (omitted when Bisection has been used)
*
                     IF( RESID.GT.TOL*GAP .AND. ABS( RQCORR ).GT.
     $                    RQTOL*ABS( LAMBDA ) .AND. .NOT. USEDBS) 
     $                    THEN
*                       We need to check that the RQCORR update doesn't
*                       move the eigenvalue away from the desired one and
*                       towards a neighbor. -> protection with bisection
                        IF(INDEIG.LE.NEGCNT) THEN
*                          The wanted eigenvalue lies to the left
                           SGNDEF = -ONE
                        ELSE
*                          The wanted eigenvalue lies to the right
                           SGNDEF = ONE
                        ENDIF
*                       We only use the RQCORR if it improves the
*                       the iterate reasonably.
                        IF( ( RQCORR*SGNDEF.GE.ZERO )
     $                       .AND.( LAMBDA + RQCORR.LE. RIGHT)
     $                       .AND.( LAMBDA + RQCORR.GE. LEFT)
     $                       ) THEN
                           USEDRQ = .TRUE.
*                          Store new midpoint of bisection interval in WORK
                           IF(SGNDEF.EQ.ONE) THEN
*                             The current LAMBDA is on the left of the true
*                             eigenvalue
                              LEFT = LAMBDA
                           ELSE   
*                             The current LAMBDA is on the right of the true
*                             eigenvalue
                              RIGHT = LAMBDA
                           ENDIF  
                           WORK( WINDEX ) = 
     $                       HALF * (RIGHT + LEFT)
*                          Take RQCORR since it has the correct sign and
*                          improves the iterate reasonably
                           LAMBDA = LAMBDA + RQCORR
*                          Update width of error interval
                           WERR( WINDEX ) =                
     $                             HALF * (RIGHT-LEFT)
                        ELSE
                           NEEDBS = .TRUE.
                        ENDIF
                        IF(RIGHT-LEFT.LT.RQTOL*ABS(LAMBDA)) THEN
*                             The eigenvalue is computed to bisection accuracy
*                             compute eigenvector and stop
                           USEDBS = .TRUE.
                           GOTO 120
                        ELSEIF( ITER.LT.MAXITR ) THEN
                           GOTO 120
                        ELSEIF( ITER.EQ.MAXITR ) THEN
                           NEEDBS = .TRUE.
                           GOTO 120
                        ELSE
                           INFO = 5
                           RETURN
                        END IF
                     ELSE 
                        STP2II = .FALSE.
                     	IF(USEDRQ .AND. USEDBS .AND. 
     $                     BSTRES.LE.RESID) THEN
                           LAMBDA = BSTW
                           STP2II = .TRUE.
                        ENDIF
                        IF (STP2II) THEN
                           CALL SLAR1VA( IN, 1, IN, LAMBDA,
     $                          D( IBEGIN ), L( IBEGIN ), 
     $                          WORK(INDLD+IBEGIN-1), 
     $                          WORK(INDLLD+IBEGIN-1),
     $                          PIVMIN, GAPTOL, 
     $                          Z( IBEGIN, ZINDEX ),
     $                          .NOT.USEDBS, NEGCNT, ZTZ, MINGMA,
     $                          ITWIST, 
     $                          ISUPPZ( 2*WINDEX-1 ),
     $                          NRMINV, RESID, RQCORR, WORK( INDWRK ) )
                        ENDIF
                        WORK( WINDEX ) = LAMBDA
                     END IF
*
*                    Compute FP-vector support w.r.t. whole matrix
*
                     ISUPPZ( 2*WINDEX-1 ) = ISUPPZ( 2*WINDEX-1 )+OLDIEN
                     ISUPPZ( 2*WINDEX ) = ISUPPZ( 2*WINDEX )+OLDIEN
                     ZFROM = ISUPPZ( 2*WINDEX-1 )
                     ZTO = ISUPPZ( 2*WINDEX )
                     ISUPMN = ISUPMN + OLDIEN
                     ISUPMX = ISUPMX + OLDIEN
*                    Ensure vector is ok if support in the RQI has changed
                     IF(ISUPMN.LT.ZFROM) THEN
                        DO 122 II = ISUPMN,ZFROM-1
                           Z( II, ZINDEX ) = ZERO
 122                    CONTINUE
                     ENDIF   
                     IF(ISUPMX.GT.ZTO) THEN
                        DO 123 II = ZTO+1,ISUPMX
                           Z( II, ZINDEX ) = ZERO
 123                    CONTINUE
                     ENDIF   
                     CALL SSCAL( ZTO-ZFROM+1, NRMINV,
     $                       Z( ZFROM, ZINDEX ), 1 )
 125                 CONTINUE
*                    Update W 
                     W( WINDEX ) = LAMBDA+SIGMA
*                    Recompute the gaps on the left and right
*                    But only allow them to become larger and not
*                    smaller (which can only happen through "bad"
*                    cancellation and doesn't reflect the theory
*                    where the initial gaps are underestimated due
*                    to WERR being too crude.)
                     IF(.NOT.ESKIP) THEN
                        IF( K.GT.1) THEN
                           WGAP( WINDMN ) = MAX( WGAP(WINDMN),
     $                          W(WINDEX)-WERR(WINDEX) 
     $                          - W(WINDMN)-WERR(WINDMN) )
                        ENDIF
                        IF( WINDEX.LT.WEND ) THEN
                           WGAP( WINDEX ) = MAX( SAVEGP, 
     $                          W( WINDPL )-WERR( WINDPL ) 
     $                          - W( WINDEX )-WERR( WINDEX) )
                        ENDIF
                     ENDIF
                  ENDIF
*                 here ends the code for the current child
*
 139              CONTINUE
*                 Proceed to any remaining child nodes
                  NEWFST = J + 1
 140           CONTINUE
 150        CONTINUE
*           Store number of clusters             
            IWORK( IINCLS + JBLK ) = NCLUS
*
         END IF
         IBEGIN = IEND + 1
         WBEGIN = WEND + 1
 170  CONTINUE
*
*     Check if everything is done: no clusters left for 
*     this processor in any block
*
      FINISH = .TRUE.
      DO 180 JBLK = 1, IBLOCK( M )      
         FINISH = FINISH .AND. (IWORK(IINCLS + JBLK).EQ.0)
 180  CONTINUE

      IF(.NOT.FINISH) THEN
         NDEPTH = NDEPTH + 1
         IF((NEEDIL.GE.DOL).AND.(NEEDIU.LE.DOU)) THEN
*           Once this processor's part of the 
*           representation tree consists exclusively of eigenvalues
*           between DOL and DOU, it can work independently from all 
*           others.
            GOTO 40
         ENDIF
      ENDIF
*

      RETURN
*
*     End of SLARRV2
*
      END
