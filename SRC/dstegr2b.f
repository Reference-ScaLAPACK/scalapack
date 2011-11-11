      SUBROUTINE DSTEGR2B( JOBZ, N, D, E, 
     $                   M, W, Z, LDZ, NZC, ISUPPZ, WORK, LWORK, IWORK,
     $                   LIWORK, DOL, DOU, NEEDIL, NEEDIU,
     $                   INDWLC, PIVMIN, SCALE, WL, WU,
     $                   VSTART, FINISH, MAXCLS,
     $                   NDEPTH, PARITY, ZOFFSET, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     July 4, 2010
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ
      INTEGER            DOL, DOU, INDWLC, INFO, LDZ, LIWORK, LWORK, M,
     $                   MAXCLS, N, NDEPTH, NEEDIL, NEEDIU, NZC, PARITY,
     $                   ZOFFSET

      DOUBLE PRECISION PIVMIN, SCALE, WL, WU
      LOGICAL VSTART, FINISH

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
*  DSTEGR2B should only be called after a call to DSTEGR2A.
*  From eigenvalues and initial representations computed by DSTEGR2A,
*  DSTEGR2B computes the selected eigenvalues and eigenvectors
*  of the real symmetric tridiagonal matrix in parallel 
*  on multiple processors. It is potentially invoked multiple times
*  on a given processor because the locally relevant representation tree 
*  might depend on spectral information that is "owned" by other processors
*  and might need to be communicated. 
* 
*  Please note:
*  1. The calling sequence has two additional INTEGER parameters, 
*     DOL and DOU, that should satisfy M>=DOU>=DOL>=1. 
*     These parameters are only relevant for the case JOBZ = 'V'.
*     DSTEGR2B  ONLY computes the eigenVECTORS 
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
*  3. DSTEGR2B needs more workspace in Z than the sequential DSTEGR. 
*     It is used to store the conformal embedding of the local representation tree.  
*  
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
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
*  M       (input) INTEGER
*          The total number of eigenvalues found
*          in DSTEGR2A.  0 <= M <= N.
*
*  W       (input) DOUBLE PRECISION array, dimension (N)
*          The first M elements contain approximations to the selected 
*          eigenvalues in ascending order. Note that only the eigenvalues from
*          the locally relevant part of the representation tree, that is
*          all the clusters that include eigenvalues from DOL:DOU, are reliable 
*          on this processor. (It does not need to know about any others anyway.)
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M) )
*          If JOBZ = 'V', and if INFO = 0, then 
*          a subset of the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix T
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          See DOL, DOU for more information.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', then LDZ >= max(1,N).
*
*  NZC     (input) INTEGER
*          The number of eigenvectors to be held in the array Z.  
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
*  NEEDIL  (input/output) INTEGER 
*  NEEDIU  (input/output) INTEGER
*          Describes which are the left and right outermost eigenvalues 
*          still to be computed. Initially computed by DLARRE2A,
*          modified in the course of the algorithm.
*
*  INDWLC  (output) DOUBLE PRECISION
*          Pointer into the workspace, location where the local
*          eigenvalue representations are stored. ("Local eigenvalues"
*          are those relative to the individual shifts of the RRRs.)
*
*  PIVMIN  (input) DOUBLE PRECISION
*          The minimum pivot in the sturm sequence for T.
*
*  SCALE   (input) DOUBLE PRECISION 
*          The scaling factor for T. Used for unscaling the eigenvalues
*          at the very end of the algorithm.
*
*  WL      (input) DOUBLE PRECISION
*  WU      (input) DOUBLE PRECISION
*          The interval (WL, WU] contains all the wanted eigenvalues.         
*
*  VSTART  (input/output) LOGICAL 
*          .TRUE. on initialization, set to .FALSE. afterwards.
*
*  FINISH  (input/output) LOGICAL
*          indicates whether all eigenpairs have been computed
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
*          in 1D-cyclic fashion
*
*  INFO    (output) INTEGER
*          On exit, INFO
*          = 0:  successful exit
*          other:if INFO = -i, the i-th argument had an illegal value
*                if INFO = 20X, internal error in DLARRV2.
*                Here, the digit X = ABS( IINFO ) < 10, where IINFO is 
*                the nonzero error code returned by DLARRV2.
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, FOUR, MINRGP
      PARAMETER          ( ONE = 1.0D0,
     $                     FOUR = 4.0D0,
     $                     MINRGP = 1.0D-3 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, WANTZ, ZQUERY
      INTEGER            IINDBL, IINDW, IINDWK, IINFO, IINSPL, INDERR,
     $                   INDGP, INDGRS, INDSDM, INDWRK, ITMP, J, LIWMIN,
     $                   LWMIN
      DOUBLE PRECISION   EPS, RTOL1, RTOL2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, DLAMCH, DLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARRV2, DSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
*
      LQUERY = ( ( LWORK.EQ.-1 ).OR.( LIWORK.EQ.-1 ) )
      ZQUERY = ( NZC.EQ.-1 )

*     DSTEGR2B needs WORK of size 6*N, IWORK of size 3*N.
*     In addition, DLARRE2A needed WORK of size 6*N, IWORK of size 5*N.
*     Workspace is kept consistent even though DLARRE2A is not called here.
*     Furthermore, DLARRV2 needs WORK of size 12*N, IWORK of size 7*N.
      IF( WANTZ ) THEN
         LWMIN = 18*N
         LIWMIN = 10*N
      ELSE
*        need less workspace if only the eigenvalues are wanted         
         LWMIN = 12*N
         LIWMIN = 8*N
      ENDIF
*
      INFO = 0
*
*     Get machine constants.
*
      EPS = DLAMCH( 'Precision' )
*
      IF( (N.EQ.0).OR.(N.EQ.1) ) THEN 
         FINISH = .TRUE.       
         RETURN
      ENDIF

      IF(ZQUERY.OR.LQUERY)
     $   RETURN
*
      INDGRS = 1
      INDERR = 2*N + 1
      INDGP = 3*N + 1
      INDSDM = 4*N + 1
      INDWRK = 6*N + 1
      INDWLC = INDWRK
*
      IINSPL = 1
      IINDBL = N + 1
      IINDW = 2*N + 1
      IINDWK = 3*N + 1

*     Set the tolerance parameters for bisection
      RTOL1 = FOUR*SQRT(EPS)
      RTOL2 = MAX( SQRT(EPS)*5.0D-3, FOUR * EPS )


      IF( WANTZ ) THEN
*
*        Compute the desired eigenvectors corresponding to the computed
*        eigenvalues
*
         CALL DLARRV2( N, WL, WU, D, E,
     $                PIVMIN, IWORK( IINSPL ), M, 
     $                DOL, DOU, NEEDIL, NEEDIU, MINRGP, RTOL1, RTOL2, 
     $                W, WORK( INDERR ), WORK( INDGP ), IWORK( IINDBL ),
     $                IWORK( IINDW ), WORK( INDGRS ), 
     $                WORK( INDSDM ), Z, LDZ,
     $                ISUPPZ, WORK( INDWRK ), IWORK( IINDWK ), 
     $                VSTART, FINISH, 
     $                MAXCLS, NDEPTH, PARITY, ZOFFSET, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 200 + ABS( IINFO )
            RETURN
         END IF
*
      ELSE
*        DLARRE2A computed eigenvalues of the (shifted) root representation
*        DLARRV2 returns the eigenvalues of the unshifted matrix.
*        However, if the eigenvectors are not desired by the user, we need
*        to apply the corresponding shifts from DLARRE2A to obtain the 
*        eigenvalues of the original matrix. 
         DO 30 J = 1, M
            ITMP = IWORK( IINDBL+J-1 )
            W( J ) = W( J ) + E( IWORK( IINSPL+ITMP-1 ) )
 30      CONTINUE
*
         FINISH = .TRUE.
*
      END IF
*

      IF(FINISH) THEN        
*        All eigenpairs have been computed       

*
*        If matrix was scaled, then rescale eigenvalues appropriately.
*
         IF( SCALE.NE.ONE ) THEN
            CALL DSCAL( M, ONE / SCALE, W, 1 )
         END IF
*
*        Correct M if needed 
*
         IF ( WANTZ ) THEN
            IF( DOL.NE.1 .OR. DOU.NE.M ) THEN
               M = DOU - DOL +1
            ENDIF
         ENDIF
*
*        No sorting of eigenpairs is done here, done later in the
*        calling subroutine
*
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      ENDIF

      RETURN
*
*     End of DSTEGR2B
*
      END
