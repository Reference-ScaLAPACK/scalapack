      SUBROUTINE DSTEGR2( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
     $                   M, W, Z, LDZ, NZC, ISUPPZ, WORK, LWORK, IWORK,
     $                   LIWORK, DOL, DOU, ZOFFSET, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     July 4, 2010
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE
      INTEGER            DOL, DOU, IL, INFO, IU, 
     $                   LDZ, NZC, LIWORK, LWORK, M, N, ZOFFSET
      DOUBLE PRECISION VL, VU

*     ..
*     .. Array Arguments ..
      INTEGER            ISUPPZ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
      DOUBLE PRECISION   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSTEGR2 computes selected eigenvalues and, optionally, eigenvectors
*  of a real symmetric tridiagonal matrix T. It is invoked in the 
*  ScaLAPACK MRRR driver PDSYEVR and the corresponding Hermitian
*  version either when only eigenvalues are to be computed, or when only
*  a single processor is used (the sequential-like case).
*
*  DSTEGR2 has been adapted from LAPACK's DSTEGR. Please note the
*  following crucial changes.
*
*  1. The calling sequence has two additional INTEGER parameters, 
*     DOL and DOU, that should satisfy M>=DOU>=DOL>=1. 
*     DSTEGR2  ONLY computes the eigenpairs
*     corresponding to eigenvalues DOL through DOU in W. (That is,
*     instead of computing the eigenpairs belonging to W(1) 
*     through W(M), only the eigenvectors belonging to eigenvalues
*     W(DOL) through W(DOU) are computed. In this case, only the
*     eigenvalues DOL:DOU are guaranteed to be fully accurate.
*
*  2. M is NOT the number of eigenvalues specified by RANGE, but is 
*     M = DOU - DOL + 1. This concerns the case where only eigenvalues
*     are computed, but on more than one processor. Thus, in this case
*     M refers to the number of eigenvalues computed on this processor.
*  
*  3. The arrays W and Z might not contain all the wanted eigenpairs
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
*          The first M elements contain the selected eigenvalues in
*          ascending order. Note that immediately after exiting this  
*          routine, only the eigenvalues from
*          position DOL:DOU are to reliable on this processor
*          because the eigenvalue computation is done in parallel.          
*          Other processors will hold reliable information on other
*          parts of the W array. This information is communicated in
*          the ScaLAPACK driver.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) )
*          If JOBZ = 'V', and if INFO = 0, then the first M columns of Z
*          contain some of the orthonormal eigenvectors of the matrix T
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and can be computed with a workspace
*          query by setting NZC = -1, see below.
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
*  ISUPPZ  (output) INTEGER ARRAY, dimension ( 2*max(1,M) )
*          The support of the eigenvectors in Z, i.e., the indices
*          indicating the nonzero elements in Z. The i-th computed eigenvector
*          is nonzero only in elements ISUPPZ( 2*i-1 ) through
*          ISUPPZ( 2*i ). This is relevant in the case when the matrix 
*          is split. ISUPPZ is only set if N>2.
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
*          From the eigenvalues W(1:M), only eigenvectors 
*          Z(:,DOL) to Z(:,DOU) are computed.
*          If DOL > 1, then Z(:,DOL-1-ZOFFSET) is used and overwritten.
*          If DOU < M, then Z(:,DOU+1-ZOFFSET) is used and overwritten.
*
*  ZOFFSET (input) INTEGER
*          Offset for storing the eigenpairs when Z is distributed
*          in 1D-cyclic fashion
*
*  INFO    (output) INTEGER
*          On exit, INFO
*          = 0:  successful exit
*          other:if INFO = -i, the i-th argument had an illegal value
*                if INFO = 10X, internal error in DLARRE2,
*                if INFO = 20X, internal error in DLARRV.
*                Here, the digit X = ABS( IINFO ) < 10, where IINFO is 
*                the nonzero error code returned by DLARRE2 or 
*                DLARRV, respectively.
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
      INTEGER            I, IIL, IINDBL, IINDW, IINDWK, IINFO, IINSPL,
     $                   IIU, INDE2, INDERR, INDGP, INDGRS, INDWRK,
     $                   ITMP, ITMP2, J, JJ, LIWMIN, LWMIN, NSPLIT,
     $                   NZCMIN
      DOUBLE PRECISION   BIGNUM, EPS, PIVMIN, RMAX, RMIN, RTOL1, RTOL2,
     $                   SAFMIN, SCALE, SMLNUM, THRESH, TMP, TNRM, WL,
     $                   WU
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, DLAMCH, DLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAE2, DLAEV2, DLARRC, DLARRE2,
     $                   DLARRV, DLASRT, DSCAL, DSWAP
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

*     DSTEGR2 needs WORK of size 6*N, IWORK of size 3*N.
*     In addition, DLARRE2 needs WORK of size 6*N, IWORK of size 5*N.
*     Furthermore, DLARRV needs WORK of size 12*N, IWORK of size 7*N.
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
*        It is either given by the user or computed in DLARRE2.
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
C         CALL XERBLA( 'DSTEGR2', -INFO )
*
         RETURN
      ELSE IF( LQUERY .OR. ZQUERY ) THEN
         RETURN
      END IF
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
*     THRESH is the splitting parameter for DLARRE2      
*     A negative THRESH forces the old splitting criterion based on the
*     size of the off-diagonal. A positive THRESH switches to splitting
*     which preserves relative accuracy. 
*
      IINFO = -1
*     Set the splitting criterion
      IF (IINFO.EQ.0) THEN
         THRESH = EPS
      ELSE
         THRESH = -EPS
      ENDIF
*
*     Store the squares of the offdiagonal values of T
      DO 5 J = 1, N-1
         WORK( INDE2+J-1 ) = E(J)**2
 5    CONTINUE

*     Set the tolerance parameters for bisection
      IF( .NOT.WANTZ ) THEN
*        DLARRE2 computes the eigenvalues to full precision.   
         RTOL1 = FOUR * EPS
         RTOL2 = FOUR * EPS
      ELSE   
*        DLARRE2 computes the eigenvalues to less than full precision.
*        DLARRV will refine the eigenvalue approximations, and we can
*        need less accurate initial bisection in DLARRE2.
*        Note: these settings do only affect the subset case and DLARRE2
         RTOL1 = SQRT(EPS)
         RTOL2 = MAX( SQRT(EPS)*5.0D-3, FOUR * EPS )
      ENDIF
      CALL DLARRE2( RANGE, N, WL, WU, IIL, IIU, D, E, 
     $             WORK(INDE2), RTOL1, RTOL2, THRESH, NSPLIT, 
     $             IWORK( IINSPL ), M, DOL, DOU,
     $             W, WORK( INDERR ),
     $             WORK( INDGP ), IWORK( IINDBL ),
     $             IWORK( IINDW ), WORK( INDGRS ), PIVMIN,
     $             WORK( INDWRK ), IWORK( IINDWK ), IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = 100 + ABS( IINFO )
         RETURN
      END IF
*     Note that if RANGE .NE. 'V', DLARRE2 computes bounds on the desired
*     part of the spectrum. All desired eigenvalues are contained in
*     (WL,WU]


      IF( WANTZ ) THEN
*
*        Compute the desired eigenvectors corresponding to the computed
*        eigenvalues
*
         CALL DLARRV( N, WL, WU, D, E,
     $                PIVMIN, IWORK( IINSPL ), M, 
     $                DOL, DOU, MINRGP, RTOL1, RTOL2, 
     $                W, WORK( INDERR ), WORK( INDGP ), IWORK( IINDBL ),
     $                IWORK( IINDW ), WORK( INDGRS ), Z, LDZ,
     $                ISUPPZ, WORK( INDWRK ), IWORK( IINDWK ), IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 200 + ABS( IINFO )
            RETURN
         END IF
      ELSE
*        DLARRE2 computes eigenvalues of the (shifted) root representation
*        DLARRV returns the eigenvalues of the unshifted matrix.
*        However, if the eigenvectors are not desired by the user, we need
*        to apply the corresponding shifts from DLARRE2 to obtain the 
*        eigenvalues of the original matrix. 
         DO 20 J = 1, M
            ITMP = IWORK( IINDBL+J-1 )
            W( J ) = W( J ) + E( IWORK( IINSPL+ITMP-1 ) )
 20      CONTINUE
      END IF
*

*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( SCALE.NE.ONE ) THEN
         CALL DSCAL( M, ONE / SCALE, W, 1 )
      END IF
*
*     Correct M if needed 
*
      IF ( WANTZ ) THEN
         IF( DOL.NE.1 .OR. DOU.NE.M ) THEN
            M = DOU - DOL +1
         ENDIF
      ENDIF
*
*     If eigenvalues are not in increasing order, then sort them, 
*     possibly along with eigenvectors.
*
      IF( NSPLIT.GT.1 ) THEN
         IF( .NOT. WANTZ ) THEN
            CALL DLASRT( 'I', DOU - DOL +1, W(DOL), IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = 3
               RETURN
            END IF
         ELSE
            DO 60 J = DOL, DOU - 1
               I = 0
               TMP = W( J )
               DO 50 JJ = J + 1, M
                  IF( W( JJ ).LT.TMP ) THEN
                     I = JJ
                     TMP = W( JJ )
                  END IF
 50            CONTINUE
               IF( I.NE.0 ) THEN
                  W( I ) = W( J )
                  W( J ) = TMP
                  IF( WANTZ ) THEN
                     CALL DSWAP( N, Z( 1, I-ZOFFSET ), 
     $                                 1, Z( 1, J-ZOFFSET ), 1 )
                     ITMP = ISUPPZ( 2*I-1 )
                     ISUPPZ( 2*I-1 ) = ISUPPZ( 2*J-1 )
                     ISUPPZ( 2*J-1 ) = ITMP
                     ITMP = ISUPPZ( 2*I )
                     ISUPPZ( 2*I ) = ISUPPZ( 2*J )
                     ISUPPZ( 2*J ) = ITMP
                  END IF
               END IF
 60         CONTINUE
         END IF
      ENDIF
*
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of DSTEGR2
*
      END
