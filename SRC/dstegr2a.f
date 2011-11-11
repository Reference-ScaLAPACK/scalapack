      SUBROUTINE DSTEGR2A( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
     $                   M, W, Z, LDZ, NZC, WORK, LWORK, IWORK,
     $                   LIWORK, DOL, DOU, NEEDIL, NEEDIU,
     $                   INDERR, NSPLIT, PIVMIN, SCALE, WL, WU,
     $                   INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     July 4, 2010
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE
      INTEGER            DOL, DOU, IL, INDERR, INFO, IU, LDZ, LIWORK,
     $                   LWORK, M, N, NEEDIL, NEEDIU, NSPLIT, NZC
      DOUBLE PRECISION PIVMIN, SCALE, VL, VU, WL, WU

*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
      DOUBLE PRECISION   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSTEGR2A computes selected eigenvalues and initial representations.
*  needed for eigenvector computations in DSTEGR2B. It is invoked in the 
*  ScaLAPACK MRRR driver PDSYEVR and the corresponding Hermitian
*  version when both eigenvalues and eigenvectors are computed in parallel.
*  on multiple processors. For this case, DSTEGR2A implements the FIRST 
*  part of the MRRR algorithm, parallel eigenvalue computation and finding
*  the root RRR. At the end of DSTEGR2A,
*  other processors might have a part of the spectrum that is needed to
*  continue the computation locally. Once this eigenvalue information has
*  been received by the processor, the computation can then proceed by calling 
*  the SECOND part of the parallel MRRR algorithm, DSTEGR2B.
*
*  Please note:
*  1. The calling sequence has two additional INTEGER parameters, 
*     (compared to LAPACK's DSTEGR), these are
*     DOL and DOU and should satisfy M>=DOU>=DOL>=1. 
*     These parameters are only relevant for the case JOBZ = 'V'.
*
*     Globally invoked over all processors, DSTEGR2A computes 
*     ALL the eigenVALUES specified by RANGE. 
*     RANGE= 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in (VL,VU] will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*     DSTEGR2A LOCALLY only computes the eigenvalues 
*     corresponding to eigenvalues DOL through DOU in W. (That is,
*     instead of computing the eigenvectors belonging to W(1) 
*     through W(M), only the eigenvectors belonging to eigenvalues
*     W(DOL) through W(DOU) are computed. In this case, only the
*     eigenvalues DOL:DOU are guaranteed to be fully accurate.
*
*  2. M is NOT the number of eigenvalues specified by RANGE, but it is 
*     M = DOU - DOL + 1. Instead, M refers to the number of eigenvalues computed on 
*     this processor.
*
*  3. While no eigenvectors are computed in DSTEGR2A itself (this is
*     done later in DSTEGR2B), the interface
*     If JOBZ = 'V' then, depending on RANGE and DOL, DOU, DSTEGR2A 
*     might need more workspace in Z then the original DSTEGR. 
*     In particular, the arrays W and Z might not contain all the wanted eigenpairs
*     locally, instead this information is distributed over other 
*     processors.
*  
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the N diagonal elements of the tridiagonal matrix
*          T. On exit, D is overwritten.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the (N-1) subdiagonal elements of the tridiagonal
*          matrix T in elements 1 to N-1 of E. E(N) need not be set on
*          input, but is used internally as workspace.
*          On exit, E is overwritten.
*
*  VL      (input) DOUBLE PRECISION
*  VU      (input) DOUBLE PRECISION
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  M       (output) INTEGER
*          Globally summed over all processors, M equals 
*          the total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*          The local output equals M = DOU - DOL + 1.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          The first M elements contain approximations to the selected 
*          eigenvalues in ascending order. Note that immediately after 
*          exiting this routine, only the eigenvalues from
*          position DOL:DOU are to reliable on this processor
*          because the eigenvalue computation is done in parallel.          
*          The other entries outside DOL:DOU are very crude preliminary
*          approximations. Other processors hold reliable information on 
*          these other parts of the W array. 
*          This information is communicated in the ScaLAPACK driver.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) )
*          DSTEGR2A does not compute eigenvectors, this is done 
*          in DSTEGR2B. The argument Z as well as all related
*          other arguments only appear to keep the interface consistent
*          and to signal to the user that this subroutine is meant to 
*          be used when eigenvectors are computed.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', then LDZ >= max(1,N).
*
*  NZC     (input) INTEGER
*          The number of eigenvectors to be held in the array Z.  
*          If RANGE = 'A', then NZC >= max(1,N).
*          If RANGE = 'V', then NZC >= the number of eigenvalues in (VL,VU].
*          If RANGE = 'I', then NZC >= IU-IL+1.
*          If NZC = -1, then a workspace query is assumed; the
*          routine calculates the number of columns of the array Z that
*          are needed to hold the eigenvectors. 
*          This value is returned as the first entry of the Z array, and
*          no error message related to NZC is issued.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal
*          (and minimal) LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,18*N)
*          if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'.
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued.
*
*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.  LIWORK >= max(1,10*N)
*          if the eigenvectors are desired, and LIWORK >= max(1,8*N)
*          if only the eigenvalues are to be computed.
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued.
*
*  DOL     (input) INTEGER
*  DOU     (input) INTEGER
*          From all the eigenvalues W(1:M), only eigenvalues
*          W(DOL:DOU) are computed.
*
*  NEEDIL  (output) INTEGER
*  NEEDIU  (output) INTEGER
*          The indices of the leftmost and rightmost eigenvalues
*          needed to accurately compute the relevant part of the 
*          representation tree. This information can be used to 
*          find out which processors have the relevant eigenvalue
*          information needed so that it can be communicated.
*
*  INDERR  (output) INTEGER
*          INDERR points to the place in the work space where 
*          the eigenvalue uncertainties (errors) are stored.
*
*  NSPLIT  (output) INTEGER
*          The number of blocks T splits into. 1 <= NSPLIT <= N.
*
*  PIVMIN  (output) DOUBLE PRECISION
*          The minimum pivot in the sturm sequence for T.
*
*  SCALE   (output) DOUBLE PRECISION 
*          The scaling factor for the tridiagonal T.
*
*  WL      (output) DOUBLE PRECISION
*  WU      (output) DOUBLE PRECISION
*          The interval (WL, WU] contains all the wanted eigenvalues.         
*          It is either given by the user or computed in DLARRE2A.
*
*  INFO    (output) INTEGER
*          On exit, INFO
*          = 0:  successful exit
*          other:if INFO = -i, the i-th argument had an illegal value
*                if INFO = 10X, internal error in DLARRE2A,
*                Here, the digit X = ABS( IINFO ) < 10, where IINFO is 
*                the nonzero error code returned by DLARRE2A.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, FOUR, MINRGP
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0,
     $                     FOUR = 4.0D0,
     $                     MINRGP = 1.0D-3 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LQUERY, VALEIG, WANTZ, ZQUERY
      INTEGER            IIL, IINDBL, IINDW, IINDWK, IINFO, IINSPL, IIU,
     $                   INDE2, INDGP, INDGRS, INDSDM, INDWRK, ITMP,
     $                   ITMP2, J, LIWMIN, LWMIN, NZCMIN
      DOUBLE PRECISION   BIGNUM, EPS, RMAX, RMIN, RTOL1, RTOL2, SAFMIN,
     $                   SMLNUM, THRESH, TNRM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, DLAMCH, DLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARRC, DLARRE2A, DSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
*
      LQUERY = ( ( LWORK.EQ.-1 ).OR.( LIWORK.EQ.-1 ) )
      ZQUERY = ( NZC.EQ.-1 )

*     DSTEGR2A needs WORK of size 6*N, IWORK of size 3*N.
*     In addition, DLARRE2A needs WORK of size 6*N, IWORK of size 5*N.
*     Furthermore, DLARRV2 needs WORK of size 12*N, IWORK of size 7*N.
*     Workspace is kept consistent with DSTEGR2B even though 
*     DLARRV2 is not called here.
      IF( WANTZ ) THEN
         LWMIN = 18*N
         LIWMIN = 10*N
      ELSE
*        need less workspace if only the eigenvalues are wanted         
         LWMIN = 12*N
         LIWMIN = 8*N
      ENDIF

      WL = ZERO
      WU = ZERO
      IIL = 0
      IIU = 0

      IF( VALEIG ) THEN
*        We do not reference VL, VU in the cases RANGE = 'I','A'
*        The interval (WL, WU] contains all the wanted eigenvalues.         
*        It is either given by the user or computed in DLARRE2A.
         WL = VL
         WU = VU
      ELSEIF( INDEIG ) THEN
*        We do not reference IL, IU in the cases RANGE = 'V','A'
         IIL = IL
         IIU = IU
      ENDIF  
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( VALEIG .AND. N.GT.0 .AND. WU.LE.WL ) THEN
         INFO = -7
      ELSE IF( INDEIG .AND. ( IIL.LT.1 .OR. IIL.GT.N ) ) THEN
         INFO = -8
      ELSE IF( INDEIG .AND. ( IIU.LT.IIL .OR. IIU.GT.N ) ) THEN
         INFO = -9
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -13
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -17
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -19
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
*
         IF( WANTZ .AND. ALLEIG ) THEN
            NZCMIN = N
            IIL = 1
            IIU = N
         ELSE IF( WANTZ .AND. VALEIG ) THEN
            CALL DLARRC( 'T', N, VL, VU, D, E, SAFMIN, 
     $                            NZCMIN, ITMP, ITMP2, INFO )
            IIL = ITMP+1
            IIU = ITMP2
         ELSE IF( WANTZ .AND. INDEIG ) THEN
            NZCMIN = IIU-IIL+1
         ELSE 
*           WANTZ .EQ. FALSE.   
            NZCMIN = 0
         ENDIF  
         IF( ZQUERY .AND. INFO.EQ.0 ) THEN
            Z( 1,1 ) = NZCMIN
         ELSE IF( NZC.LT.NZCMIN .AND. .NOT.ZQUERY ) THEN
            INFO = -14
         END IF
      END IF

      IF ( WANTZ ) THEN
         IF ( DOL.LT.1 .OR. DOL.GT.NZCMIN ) THEN 
            INFO = -20
         ENDIF
         IF ( DOU.LT.1 .OR. DOU.GT.NZCMIN .OR. DOU.LT.DOL) THEN 
            INFO = -21
         ENDIF
      ENDIF

      IF( INFO.NE.0 ) THEN
*
C         Disable sequential error handler
C         for parallel case
C         CALL XERBLA( 'DSTEGR2A', -INFO )
*
         RETURN
      ELSE IF( LQUERY .OR. ZQUERY ) THEN
         RETURN
      END IF

*     Initialize NEEDIL and NEEDIU, these values are changed in DLARRE2A
      NEEDIL = DOU
      NEEDIU = DOL
*
*     Quick return if possible
*
      M = 0
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = D( 1 )
         ELSE
            IF( WL.LT.D( 1 ) .AND. WU.GE.D( 1 ) ) THEN
               M = 1
               W( 1 ) = D( 1 )
            END IF
         END IF
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
      INDGRS = 1
      INDERR = 2*N + 1
      INDGP = 3*N + 1
      INDSDM = 4*N + 1
      INDE2 = 5*N + 1
      INDWRK = 6*N + 1
*
      IINSPL = 1
      IINDBL = N + 1
      IINDW = 2*N + 1
      IINDWK = 3*N + 1
*
*     Scale matrix to allowable range, if necessary.
*
      SCALE = ONE
      TNRM = DLANST( 'M', N, D, E )
      IF( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) THEN
         SCALE = RMIN / TNRM
      ELSE IF( TNRM.GT.RMAX ) THEN
         SCALE = RMAX / TNRM
      END IF
      IF( SCALE.NE.ONE ) THEN
         CALL DSCAL( N, SCALE, D, 1 )
         CALL DSCAL( N-1, SCALE, E, 1 )
         TNRM = TNRM*SCALE
         IF( VALEIG ) THEN
*           If eigenvalues in interval have to be found, 
*           scale (WL, WU] accordingly
            WL = WL*SCALE
            WU = WU*SCALE
         ENDIF
      END IF
*
*     Compute the desired eigenvalues of the tridiagonal after splitting
*     into smaller subblocks if the corresponding off-diagonal elements
*     are small
*     THRESH is the splitting parameter for DLARRA in DLARRE2A      
*     A negative THRESH forces the old splitting criterion based on the
*     size of the off-diagonal.
      THRESH = -EPS
      IINFO = 0

*     Store the squares of the offdiagonal values of T
      DO 5 J = 1, N-1
         WORK( INDE2+J-1 ) = E(J)**2
 5    CONTINUE

*     Set the tolerance parameters for bisection
      IF( .NOT.WANTZ ) THEN
*        DLARRE2A computes the eigenvalues to full precision.   
         RTOL1 = FOUR * EPS
         RTOL2 = FOUR * EPS
      ELSE   
*        DLARRE2A computes the eigenvalues to less than full precision.
*        DLARRV2 will refine the eigenvalue approximations, and we can
*        need less accurate initial bisection in DLARRE2A.
         RTOL1 = FOUR*SQRT(EPS)
         RTOL2 = MAX( SQRT(EPS)*5.0D-3, FOUR * EPS )
      ENDIF
      CALL DLARRE2A( RANGE, N, WL, WU, IIL, IIU, D, E, 
     $             WORK(INDE2), RTOL1, RTOL2, THRESH, NSPLIT, 
     $             IWORK( IINSPL ), M, DOL, DOU, NEEDIL, NEEDIU,
     $             W, WORK( INDERR ),
     $             WORK( INDGP ), IWORK( IINDBL ),
     $             IWORK( IINDW ), WORK( INDGRS ), 
     $             WORK( INDSDM ), PIVMIN,
     $             WORK( INDWRK ), IWORK( IINDWK ), 
     $             MINRGP, IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = 100 + ABS( IINFO )
         RETURN
      END IF
*     Note that if RANGE .NE. 'V', DLARRE2A computes bounds on the desired
*     part of the spectrum. All desired eigenvalues are contained in
*     (WL,WU]


      RETURN
*
*     End of DSTEGR2A
*
      END
