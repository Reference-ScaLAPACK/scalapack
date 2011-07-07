      SUBROUTINE PSSYEVX( JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, VL,
     $                    VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, Z, IZ,
     $                    JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, IFAIL,
     $                    ICLUSTR, GAP, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001 
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LWORK, M,
     $                   N, NZ
      REAL               ABSTOL, ORFAC, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * ), ICLUSTR( * ),
     $                   IFAIL( * ), IWORK( * )
      REAL               A( * ), GAP( * ), W( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PSSYEVX computes selected eigenvalues and, optionally, eigenvectors
*  of a real symmetric matrix A by calling the recommended sequence
*  of ScaLAPACK routines.  Eigenvalues/vectors can be selected by
*  specifying a range of values or a range of indices for the desired
*  eigenvalues.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  PSSYEVX assumes IEEE 754 standard compliant arithmetic.  To port
*  to a system which does not have IEEE 754 arithmetic, modify
*  the appropriate SLmake.inc file to include the compiler switch
*  -DNO_IEEE.  This switch only affects the compilation of pslaiect.c.
*
*  Arguments
*  =========
*
*     NP = the number of rows local to a given process.
*     NQ = the number of columns local to a given process.
*
*  JOBZ    (global input) CHARACTER*1
*          Specifies whether or not to compute the eigenvectors:
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (global input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the interval [VL,VU] will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns of the matrix A.  N >= 0.
*
*  A       (local input/workspace) block cyclic REAL array,
*          global dimension (N, N),
*          local dimension ( LLD_A, LOCc(JA+N-1) )
*
*          On entry, the symmetric matrix A.  If UPLO = 'U', only the
*          upper triangular part of A is used to define the elements of
*          the symmetric matrix.  If UPLO = 'L', only the lower
*          triangular part of A is used to define the elements of the
*          symmetric matrix.
*
*          On exit, the lower triangle (if UPLO='L') or the upper
*          triangle (if UPLO='U') of A, including the diagonal, is
*          destroyed.
*
*  IA      (global input) INTEGER
*          A's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*          If DESCA( CTXT_ ) is incorrect, PSSYEVX cannot guarantee
*          correct error reporting.
*
*  VL      (global input) REAL
*          If RANGE='V', the lower bound of the interval to be searched
*          for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
*
*  VU      (global input) REAL
*          If RANGE='V', the upper bound of the interval to be searched
*          for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (global input) INTEGER
*          If RANGE='I', the index (from smallest to largest) of the
*          smallest eigenvalue to be returned.  IL >= 1.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  IU      (global input) INTEGER
*          If RANGE='I', the index (from smallest to largest) of the
*          largest eigenvalue to be returned.  min(IL,N) <= IU <= N.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (global input) REAL
*          If JOBZ='V', setting ABSTOL to PSLAMCH( CONTEXT, 'U') yields
*          the most orthogonal eigenvectors.
*
*          The absolute error tolerance for the eigenvalues.
*          An approximate eigenvalue is accepted as converged
*          when it is determined to lie in an interval [a,b]
*          of width less than or equal to
*
*                  ABSTOL + EPS *   max( |a|,|b| ) ,
*
*          where EPS is the machine precision.  If ABSTOL is less than
*          or equal to zero, then EPS*norm(T) will be used in its place,
*          where norm(T) is the 1-norm of the tridiagonal matrix
*          obtained by reducing A to tridiagonal form.
*
*          Eigenvalues will be computed most accurately when ABSTOL is
*          set to twice the underflow threshold 2*PSLAMCH('S') not zero.
*          If this routine returns with ((MOD(INFO,2).NE.0) .OR.
*          (MOD(INFO/8,2).NE.0)), indicating that some eigenvalues or
*          eigenvectors did not converge, try setting ABSTOL to
*          2*PSLAMCH('S').
*
*          See "Computing Small Singular Values of Bidiagonal Matrices
*          with Guaranteed High Relative Accuracy," by Demmel and
*          Kahan, LAPACK Working Note #3.
*
*          See "On the correctness of Parallel Bisection in Floating
*          Point" by Demmel, Dhillon and Ren, LAPACK Working Note #70
*
*  M       (global output) INTEGER
*          Total number of eigenvalues found.  0 <= M <= N.
*
*  NZ      (global output) INTEGER
*          Total number of eigenvectors computed.  0 <= NZ <= M.
*          The number of columns of Z that are filled.
*          If JOBZ .NE. 'V', NZ is not referenced.
*          If JOBZ .EQ. 'V', NZ = M unless the user supplies
*          insufficient space and PSSYEVX is not able to detect this
*          before beginning computation.  To get all the eigenvectors
*          requested, the user must supply both sufficient
*          space to hold the eigenvectors in Z (M .LE. DESCZ(N_))
*          and sufficient workspace to compute them.  (See LWORK below.)
*          PSSYEVX is always able to detect insufficient space without
*          computation unless RANGE .EQ. 'V'.
*
*  W       (global output) REAL array, dimension (N)
*          On normal exit, the first M entries contain the selected
*          eigenvalues in ascending order.
*
*  ORFAC   (global input) REAL
*          Specifies which eigenvectors should be reorthogonalized.
*          Eigenvectors that correspond to eigenvalues which are within
*          tol=ORFAC*norm(A) of each other are to be reorthogonalized.
*          However, if the workspace is insufficient (see LWORK),
*          tol may be decreased until all eigenvectors to be
*          reorthogonalized can be stored in one process.
*          No reorthogonalization will be done if ORFAC equals zero.
*          A default value of 10^-3 is used if ORFAC is negative.
*          ORFAC should be identical on all processes.
*
*  Z       (local output) REAL array,
*          global dimension (N, N),
*          local dimension ( LLD_Z, LOCc(JZ+N-1) )
*          If JOBZ = 'V', then on normal exit the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix
*          corresponding to the selected eigenvalues.  If an eigenvector
*          fails to converge, then that column of Z contains the latest
*          approximation to the eigenvector, and the index of the
*          eigenvector is returned in IFAIL.
*          If JOBZ = 'N', then Z is not referenced.
*
*  IZ      (global input) INTEGER
*          Z's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JZ      (global input) INTEGER
*          Z's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*          DESCZ( CTXT_ ) must equal DESCA( CTXT_ )
*
*  WORK    (local workspace/output) REAL array,
*          dimension max(3,LWORK)
*          On return, WORK(1) contains the optimal amount of
*          workspace required for efficient execution.
*          if JOBZ='N' WORK(1) = optimal amount of workspace
*             required to compute eigenvalues efficiently
*          if JOBZ='V' WORK(1) = optimal amount of workspace
*             required to compute eigenvalues and eigenvectors
*             efficiently with no guarantee on orthogonality.
*             If RANGE='V', it is assumed that all eigenvectors
*             may be required.
*
*  LWORK   (local input) INTEGER
*          Size of WORK
*          See below for definitions of variables used to define LWORK.
*          If no eigenvectors are requested (JOBZ = 'N') then
*             LWORK >= 5 * N + MAX( 5 * NN, NB * ( NP0 + 1 ) )
*          If eigenvectors are requested (JOBZ = 'V' ) then
*             the amount of workspace required to guarantee that all
*             eigenvectors are computed is:
*             LWORK >= 5*N + MAX( 5*NN, NP0 * MQ0 + 2 * NB * NB ) +
*               ICEIL( NEIG, NPROW*NPCOL)*NN
*
*             The computed eigenvectors may not be orthogonal if the
*             minimal workspace is supplied and ORFAC is too small.
*             If you want to guarantee orthogonality (at the cost
*             of potentially poor performance) you should add
*             the following to LWORK:
*                (CLUSTERSIZE-1)*N
*             where CLUSTERSIZE is the number of eigenvalues in the
*             largest cluster, where a cluster is defined as a set of
*             close eigenvalues: { W(K),...,W(K+CLUSTERSIZE-1) |
*                                  W(J+1) <= W(J) + ORFAC*2*norm(A) }
*          Variable definitions:
*             NEIG = number of eigenvectors requested
*             NB = DESCA( MB_ ) = DESCA( NB_ ) =
*                  DESCZ( MB_ ) = DESCZ( NB_ )
*             NN = MAX( N, NB, 2 )
*             DESCA( RSRC_ ) = DESCA( NB_ ) = DESCZ( RSRC_ ) =
*                              DESCZ( CSRC_ ) = 0
*             NP0 = NUMROC( NN, NB, 0, 0, NPROW )
*             MQ0 = NUMROC( MAX( NEIG, NB, 2 ), NB, 0, 0, NPCOL )
*             ICEIL( X, Y ) is a ScaLAPACK function returning
*             ceiling(X/Y)
*
*          When LWORK is too small:
*             If LWORK is too small to guarantee orthogonality,
*             PSSYEVX attempts to maintain orthogonality in
*             the clusters with the smallest
*             spacing between the eigenvalues.
*             If LWORK is too small to compute all the eigenvectors
*             requested, no computation is performed and INFO=-23
*             is returned.  Note that when RANGE='V', PSSYEVX does
*             not know how many eigenvectors are requested until
*             the eigenvalues are computed.  Therefore, when RANGE='V'
*             and as long as LWORK is large enough to allow PSSYEVX to
*             compute the eigenvalues, PSSYEVX will compute the
*             eigenvalues and as many eigenvectors as it can.
*
*          Relationship between workspace, orthogonality & performance:
*             Greater performance can be achieved if adequate workspace
*             is provided.  On the other hand, in some situations,
*             performance can decrease as the workspace provided
*             increases above the workspace amount shown below:
*
*             For optimal performance, greater workspace may be
*             needed, i.e.
*                LWORK >=  MAX( LWORK, 5*N + NSYTRD_LWOPT )
*                Where:
*                  LWORK, as defined previously, depends upon the number
*                     of eigenvectors requested, and
*                  NSYTRD_LWOPT = N + 2*( ANB+1 )*( 4*NPS+2 ) +
*                    ( NPS + 3 ) *  NPS
*
*                  ANB = PJLAENV( DESCA( CTXT_), 3, 'PSSYTTRD', 'L',
*                     0, 0, 0, 0)
*                  SQNPC = INT( SQRT( DBLE( NPROW * NPCOL ) ) )
*                  NPS = MAX( NUMROC( N, 1, 0, 0, SQNPC ), 2*ANB )
*
*                  NUMROC is a ScaLAPACK tool functions;
*                  PJLAENV is a ScaLAPACK envionmental inquiry function
*                  MYROW, MYCOL, NPROW and NPCOL can be determined by
*                    calling the subroutine BLACS_GRIDINFO.
*
*                For large N, no extra workspace is needed, however the
*                biggest boost in performance comes for small N, so it
*                is wise to provide the extra workspace (typically less
*                than a Megabyte per process).
*
*             If CLUSTERSIZE >= N/SQRT(NPROW*NPCOL), then providing
*             enough space to compute all the eigenvectors
*             orthogonally will cause serious degradation in
*             performance. In the limit (i.e. CLUSTERSIZE = N-1)
*             PSSTEIN will perform no better than SSTEIN on 1
*             processor.
*             For CLUSTERSIZE = N/SQRT(NPROW*NPCOL) reorthogonalizing
*             all eigenvectors will increase the total execution time
*             by a factor of 2 or more.
*             For CLUSTERSIZE > N/SQRT(NPROW*NPCOL) execution time will
*             grow as the square of the cluster size, all other factors
*             remaining equal and assuming enough workspace.  Less
*             workspace means less reorthogonalization but faster
*             execution.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the size
*          required for optimal performance for all work arrays. Each of
*          these values is returned in the first entry of the
*          corresponding work arrays, and no error message is issued by
*          PXERBLA.
*
*  IWORK   (local workspace) INTEGER array
*          On return, IWORK(1) contains the amount of integer workspace
*          required.
*
*  LIWORK  (local input) INTEGER
*          size of IWORK
*          LIWORK >= 6 * NNP
*          Where:
*            NNP = MAX( N, NPROW*NPCOL + 1, 4 )
*          If LIWORK = -1, then LIWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  IFAIL   (global output) INTEGER array, dimension (N)
*          If JOBZ = 'V', then on normal exit, the first M elements of
*          IFAIL are zero.  If (MOD(INFO,2).NE.0) on exit, then
*          IFAIL contains the
*          indices of the eigenvectors that failed to converge.
*          If JOBZ = 'N', then IFAIL is not referenced.
*
*  ICLUSTR (global output) integer array, dimension (2*NPROW*NPCOL)
*          This array contains indices of eigenvectors corresponding to
*          a cluster of eigenvalues that could not be reorthogonalized
*          due to insufficient workspace (see LWORK, ORFAC and INFO).
*          Eigenvectors corresponding to clusters of eigenvalues indexed
*          ICLUSTR(2*I-1) to ICLUSTR(2*I), could not be
*          reorthogonalized due to lack of workspace. Hence the
*          eigenvectors corresponding to these clusters may not be
*          orthogonal.  ICLUSTR() is a zero terminated array.
*          (ICLUSTR(2*K).NE.0 .AND. ICLUSTR(2*K+1).EQ.0) if and only if
*          K is the number of clusters
*          ICLUSTR is not referenced if JOBZ = 'N'
*
*  GAP     (global output) REAL array,
*             dimension (NPROW*NPCOL)
*          This array contains the gap between eigenvalues whose
*          eigenvectors could not be reorthogonalized. The output
*          values in this array correspond to the clusters indicated
*          by the array ICLUSTR. As a result, the dot product between
*          eigenvectors correspoding to the I^th cluster may be as high
*          as ( C * n ) / GAP(I) where C is a small constant.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  if (MOD(INFO,2).NE.0), then one or more eigenvectors
*                  failed to converge.  Their indices are stored
*                  in IFAIL.  Ensure ABSTOL=2.0*PSLAMCH( 'U' )
*                  Send e-mail to scalapack@cs.utk.edu
*                if (MOD(INFO/2,2).NE.0),then eigenvectors corresponding
*                  to one or more clusters of eigenvalues could not be
*                  reorthogonalized because of insufficient workspace.
*                  The indices of the clusters are stored in the array
*                  ICLUSTR.
*                if (MOD(INFO/4,2).NE.0), then space limit prevented
*                  PSSYEVX from computing all of the eigenvectors
*                  between VL and VU.  The number of eigenvectors
*                  computed is returned in NZ.
*                if (MOD(INFO/8,2).NE.0), then PSSTEBZ failed to compute
*                  eigenvalues.  Ensure ABSTOL=2.0*PSLAMCH( 'U' )
*                  Send e-mail to scalapack@cs.utk.edu
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrices A(IA:*, JA:*) and C(IC:IC+M-1,JC:JC+N-1)
*  must verify some alignment properties, namely the following
*  expressions should be true:
*
*  ( MB_A.EQ.NB_A.EQ.MB_Z .AND. IROFFA.EQ.IROFFZ .AND. IROFFA.EQ.0 .AND.
*    IAROW.EQ.IZROW )
*  where
*  IROFFA = MOD( IA-1, MB_A ) and ICOFFA = MOD( JA-1, NB_A ).
*
*  =====================================================================
*
*  Differences between PSSYEVX and SSYEVX
*  ======================================
*
*  A, LDA -> A, IA, JA, DESCA
*  Z, LDZ -> Z, IZ, JZ, DESCZ
*  WORKSPACE needs are larger for PSSYEVX.
*  LIWORK parameter added
*
*  ORFAC, ICLUSTER() and GAP() parameters added
*  meaning of INFO is changed
*
*  Functional differences:
*  PSSYEVX does not promise orthogonality for eigenvectors associated
*  with tighly clustered eigenvalues.
*  PSSYEVX does not reorthogonalize eigenvectors
*  that are on different processes. The extent of reorthogonalization
*  is controlled by the input parameter LWORK.
*
*  Version 1.4 limitations:
*     DESCA(MB_) = DESCA(NB_)
*     DESCA(M_) = DESCZ(M_)
*     DESCA(N_) = DESCZ(N_)
*     DESCA(MB_) = DESCZ(MB_)
*     DESCA(NB_) = DESCZ(NB_)
*     DESCA(RSRC_) = DESCZ(RSRC_)
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ZERO, ONE, TEN, FIVE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TEN = 10.0E+0,
     $                   FIVE = 5.0E+0 )
      INTEGER            IERREIN, IERRCLS, IERRSPC, IERREBZ
      PARAMETER          ( IERREIN = 1, IERRCLS = 2, IERRSPC = 4,
     $                   IERREBZ = 8 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LOWER, LQUERY, QUICKRETURN,
     $                   VALEIG, WANTZ
      CHARACTER          ORDER
      INTEGER            ANB, CSRC_A, I, IAROW, ICOFFA, ICTXT, IINFO,
     $                   INDD, INDD2, INDE, INDE2, INDIBL, INDISP,
     $                   INDTAU, INDWORK, IROFFA, IROFFZ, ISCALE,
     $                   ISIZESTEBZ, ISIZESTEIN, IZROW, LALLWORK,
     $                   LIWMIN, LLWORK, LWMIN, LWOPT, MAXEIGS, MB_A,
     $                   MQ0, MYCOL, MYROW, NB, NB_A, NEIG, NN, NNP,
     $                   NP0, NPCOL, NPROCS, NPROW, NPS, NSPLIT,
     $                   NSYTRD_LWOPT, NZZ, OFFSET, RSRC_A, RSRC_Z,
     $                   SIZEORMTR, SIZESTEIN, SIZESYEVX, SQNPC
      REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN,
     $                   SIGMA, SMLNUM, VLL, VUU
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 4 ), IDUM2( 4 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, INDXG2P, NUMROC, PJLAENV
      REAL               PSLAMCH, PSLANSY
      EXTERNAL           LSAME, ICEIL, INDXG2P, NUMROC, PJLAENV,
     $                   PSLAMCH, PSLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, IGAMN2D, PCHK1MAT,
     $                   PCHK2MAT, PSELGET, PSLARED1D, PSLASCL, PSORMTR,
     $                   PSSTEBZ, PSSTEIN, PSSYNTRD, PXERBLA, SGEBR2D,
     $                   SGEBS2D, SLASRT, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, ICHAR, INT, MAX, MIN, MOD, REAL,
     $                   SQRT
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      QUICKRETURN = ( N.EQ.0 )
*
*     Test the input arguments.
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
*
      WANTZ = LSAME( JOBZ, 'V' )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 800+CTXT_ )
      ELSE IF( WANTZ ) THEN
         IF( ICTXT.NE.DESCZ( CTXT_ ) ) THEN
            INFO = -( 2100+CTXT_ )
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 4, N, 4, IA, JA, DESCA, 8, INFO )
         IF( WANTZ )
     $      CALL CHK1MAT( N, 4, N, 4, IZ, JZ, DESCZ, 21, INFO )
*
         IF( INFO.EQ.0 ) THEN
*
*     Get machine constants.
*
            SAFMIN = PSLAMCH( ICTXT, 'Safe minimum' )
            EPS = PSLAMCH( ICTXT, 'Precision' )
            SMLNUM = SAFMIN / EPS
            BIGNUM = ONE / SMLNUM
            RMIN = SQRT( SMLNUM )
            RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
            NPROCS = NPROW*NPCOL
            LOWER = LSAME( UPLO, 'L' )
            ALLEIG = LSAME( RANGE, 'A' )
            VALEIG = LSAME( RANGE, 'V' )
            INDEIG = LSAME( RANGE, 'I' )
*
*     Set up pointers into the WORK array
*
            INDTAU = 1
            INDE = INDTAU + N
            INDD = INDE + N
            INDD2 = INDD + N
            INDE2 = INDD2 + N
            INDWORK = INDE2 + N
            LLWORK = LWORK - INDWORK + 1
*
*     Set up pointers into the IWORK array
*
            ISIZESTEIN = 3*N + NPROCS + 1
            ISIZESTEBZ = MAX( 4*N, 14, NPROCS )
            INDIBL = ( MAX( ISIZESTEIN, ISIZESTEBZ ) ) + 1
            INDISP = INDIBL + N
*
*     Compute the total amount of space needed
*
            LQUERY = .FALSE.
            IF( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
     $         LQUERY = .TRUE.
*
            NNP = MAX( N, NPROCS+1, 4 )
            LIWMIN = 6*NNP
*
            NPROCS = NPROW*NPCOL
            NB_A = DESCA( NB_ )
            MB_A = DESCA( MB_ )
            NB = NB_A
            NN = MAX( N, NB, 2 )
*
            RSRC_A = DESCA( RSRC_ )
            CSRC_A = DESCA( CSRC_ )
            IROFFA = MOD( IA-1, MB_A )
            ICOFFA = MOD( JA-1, NB_A )
            IAROW = INDXG2P( 1, NB_A, MYROW, RSRC_A, NPROW )
            NP0 = NUMROC( N+IROFFA, NB, 0, 0, NPROW )
            MQ0 = NUMROC( N+ICOFFA, NB, 0, 0, NPCOL )
            IF( WANTZ ) THEN
               RSRC_Z = DESCZ( RSRC_ )
               IROFFZ = MOD( IZ-1, MB_A )
               IZROW = INDXG2P( 1, NB_A, MYROW, RSRC_Z, NPROW )
            ELSE
               IROFFZ = 0
               IZROW = 0
            END IF
*
            IF( ( .NOT.WANTZ ) .OR. ( VALEIG .AND. ( .NOT.LQUERY ) ) )
     $           THEN
               LWMIN = 5*N + MAX( 5*NN, NB*( NP0+1 ) )
               IF( WANTZ ) THEN
                  MQ0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPCOL )
                  LWOPT = 5*N + MAX( 5*NN, NP0*MQ0+2*NB*NB )
               ELSE
                  LWOPT = LWMIN
               END IF
               NEIG = 0
            ELSE
               IF( ALLEIG .OR. VALEIG ) THEN
                  NEIG = N
               ELSE IF( INDEIG ) THEN
                  NEIG = IU - IL + 1
               END IF
               MQ0 = NUMROC( MAX( NEIG, NB, 2 ), NB, 0, 0, NPCOL )
               LWMIN = 5*N + MAX( 5*NN, NP0*MQ0+2*NB*NB ) +
     $                 ICEIL( NEIG, NPROW*NPCOL )*NN
               LWOPT = LWMIN
*
            END IF
*
*           Compute how much workspace is needed to use the
*           new TRD code
*
            ANB = PJLAENV( ICTXT, 3, 'PSSYTTRD', 'L', 0, 0, 0, 0 )
            SQNPC = INT( SQRT( DBLE( NPROW*NPCOL ) ) )
            NPS = MAX( NUMROC( N, 1, 0, 0, SQNPC ), 2*ANB )
            NSYTRD_LWOPT = 2*( ANB+1 )*( 4*NPS+2 ) + ( NPS+4 )*NPS
            LWOPT = MAX( LWOPT, 5*N+NSYTRD_LWOPT )
*
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
               WORK( 1 ) = ABSTOL
               IF( VALEIG ) THEN
                  WORK( 2 ) = VL
                  WORK( 3 ) = VU
               ELSE
                  WORK( 2 ) = ZERO
                  WORK( 3 ) = ZERO
               END IF
               CALL SGEBS2D( ICTXT, 'ALL', ' ', 3, 1, WORK, 3 )
            ELSE
               CALL SGEBR2D( ICTXT, 'ALL', ' ', 3, 1, WORK, 3, 0, 0 )
            END IF
            IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
               INFO = -1
            ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
               INFO = -2
            ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
               INFO = -3
            ELSE IF( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) THEN
               INFO = -10
            ELSE IF( INDEIG .AND. ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) )
     $                THEN
               INFO = -11
            ELSE IF( INDEIG .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) )
     $                THEN
               INFO = -12
            ELSE IF( LWORK.LT.LWMIN .AND. LWORK.NE.-1 ) THEN
               INFO = -23
            ELSE IF( LIWORK.LT.LIWMIN .AND. LIWORK.NE.-1 ) THEN
               INFO = -25
            ELSE IF( VALEIG .AND. ( ABS( WORK( 2 )-VL ).GT.FIVE*EPS*
     $               ABS( VL ) ) ) THEN
               INFO = -9
            ELSE IF( VALEIG .AND. ( ABS( WORK( 3 )-VU ).GT.FIVE*EPS*
     $               ABS( VU ) ) ) THEN
               INFO = -10
            ELSE IF( ABS( WORK( 1 )-ABSTOL ).GT.FIVE*EPS*ABS( ABSTOL ) )
     $                THEN
               INFO = -13
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 800+NB_ )
            END IF
            IF( WANTZ ) THEN
               IF( IROFFA.NE.IROFFZ ) THEN
                  INFO = -19
               ELSE IF( IAROW.NE.IZROW ) THEN
                  INFO = -19
               ELSE IF( DESCA( M_ ).NE.DESCZ( M_ ) ) THEN
                  INFO = -( 2100+M_ )
               ELSE IF( DESCA( N_ ).NE.DESCZ( N_ ) ) THEN
                  INFO = -( 2100+N_ )
               ELSE IF( DESCA( MB_ ).NE.DESCZ( MB_ ) ) THEN
                  INFO = -( 2100+MB_ )
               ELSE IF( DESCA( NB_ ).NE.DESCZ( NB_ ) ) THEN
                  INFO = -( 2100+NB_ )
               ELSE IF( DESCA( RSRC_ ).NE.DESCZ( RSRC_ ) ) THEN
                  INFO = -( 2100+RSRC_ )
               ELSE IF( DESCA( CSRC_ ).NE.DESCZ( CSRC_ ) ) THEN
                  INFO = -( 2100+CSRC_ )
               ELSE IF( ICTXT.NE.DESCZ( CTXT_ ) ) THEN
                  INFO = -( 2100+CTXT_ )
               END IF
            END IF
         END IF
         IF( WANTZ ) THEN
            IDUM1( 1 ) = ICHAR( 'V' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'N' )
         END IF
         IDUM2( 1 ) = 1
         IF( LOWER ) THEN
            IDUM1( 2 ) = ICHAR( 'L' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'U' )
         END IF
         IDUM2( 2 ) = 2
         IF( ALLEIG ) THEN
            IDUM1( 3 ) = ICHAR( 'A' )
         ELSE IF( INDEIG ) THEN
            IDUM1( 3 ) = ICHAR( 'I' )
         ELSE
            IDUM1( 3 ) = ICHAR( 'V' )
         END IF
         IDUM2( 3 ) = 3
         IF( LQUERY ) THEN
            IDUM1( 4 ) = -1
         ELSE
            IDUM1( 4 ) = 1
         END IF
         IDUM2( 4 ) = 4
         IF( WANTZ ) THEN
            CALL PCHK2MAT( N, 4, N, 4, IA, JA, DESCA, 8, N, 4, N, 4, IZ,
     $                     JZ, DESCZ, 21, 4, IDUM1, IDUM2, INFO )
         ELSE
            CALL PCHK1MAT( N, 4, N, 4, IA, JA, DESCA, 8, 4, IDUM1,
     $                     IDUM2, INFO )
         END IF
         WORK( 1 ) = REAL( LWOPT )
         IWORK( 1 ) = LIWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PSSYEVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( QUICKRETURN ) THEN
         IF( WANTZ ) THEN
            NZ = 0
            ICLUSTR( 1 ) = 0
         END IF
         M = 0
         WORK( 1 ) = REAL( LWOPT )
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF
*
*     Scale matrix to allowable range, if necessary.
*
      ABSTLL = ABSTOL
      ISCALE = 0
      IF( VALEIG ) THEN
         VLL = VL
         VUU = VU
      ELSE
         VLL = ZERO
         VUU = ZERO
      END IF
*
      ANRM = PSLANSY( 'M', UPLO, N, A, IA, JA, DESCA, WORK( INDWORK ) )
*
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
         ANRM = ANRM*SIGMA
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
         ANRM = ANRM*SIGMA
      END IF
*
      IF( ISCALE.EQ.1 ) THEN
         CALL PSLASCL( UPLO, ONE, SIGMA, N, N, A, IA, JA, DESCA, IINFO )
         IF( ABSTOL.GT.0 )
     $      ABSTLL = ABSTOL*SIGMA
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
            IF( VUU.EQ.VLL ) THEN
               VUU = VUU + 2*MAX( ABS( VUU )*EPS, SAFMIN )
            END IF
         END IF
      END IF
*
*     Call PSSYNTRD to reduce symmetric matrix to tridiagonal form.
*
      LALLWORK = LLWORK
*
      CALL PSSYNTRD( UPLO, N, A, IA, JA, DESCA, WORK( INDD ),
     $               WORK( INDE ), WORK( INDTAU ), WORK( INDWORK ),
     $               LLWORK, IINFO )
*
*
*     Copy the values of D, E to all processes
*
*     Here PxLARED1D is used to redistribute the tridiagonal matrix.
*     PxLARED1D, however, doesn't yet work with arbritary matrix
*     distributions so we have PxELGET as a backup.
*
      OFFSET = 0
      IF( IA.EQ.1 .AND. JA.EQ.1 .AND. RSRC_A.EQ.0 .AND. CSRC_A.EQ.0 )
     $     THEN
         CALL PSLARED1D( N, IA, JA, DESCA, WORK( INDD ), WORK( INDD2 ),
     $                   WORK( INDWORK ), LLWORK )
*
         CALL PSLARED1D( N, IA, JA, DESCA, WORK( INDE ), WORK( INDE2 ),
     $                   WORK( INDWORK ), LLWORK )
         IF( .NOT.LOWER )
     $      OFFSET = 1
      ELSE
         DO 10 I = 1, N
            CALL PSELGET( 'A', ' ', WORK( INDD2+I-1 ), A, I+IA-1,
     $                    I+JA-1, DESCA )
   10    CONTINUE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 I = 1, N - 1
               CALL PSELGET( 'A', ' ', WORK( INDE2+I-1 ), A, I+IA-1,
     $                       I+JA, DESCA )
   20       CONTINUE
         ELSE
            DO 30 I = 1, N - 1
               CALL PSELGET( 'A', ' ', WORK( INDE2+I-1 ), A, I+IA,
     $                       I+JA-1, DESCA )
   30       CONTINUE
         END IF
      END IF
*
*     Call PSSTEBZ and, if eigenvectors are desired, PSSTEIN.
*
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
*
      CALL PSSTEBZ( ICTXT, RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL,
     $              WORK( INDD2 ), WORK( INDE2+OFFSET ), M, NSPLIT, W,
     $              IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWORK ),
     $              LLWORK, IWORK( 1 ), ISIZESTEBZ, IINFO )
*
*
*     IF PSSTEBZ fails, the error propogates to INFO, but
*     we do not propogate the eigenvalue(s) which failed because:
*     1)  This should never happen if the user specifies
*         ABSTOL = 2 * PSLAMCH( 'U' )
*     2)  PSSTEIN will confirm/deny whether the eigenvalues are
*         close enough.
*
      IF( IINFO.NE.0 ) THEN
         INFO = INFO + IERREBZ
         DO 40 I = 1, M
            IWORK( INDIBL+I-1 ) = ABS( IWORK( INDIBL+I-1 ) )
   40    CONTINUE
      END IF
      IF( WANTZ ) THEN
*
         IF( VALEIG ) THEN
*
*           Compute the maximum number of eigenvalues that we can
*           compute in the
*           workspace that we have, and that we can store in Z.
*
*           Loop through the possibilities looking for the largest
*           NZ that we can feed to PSSTEIN and PSORMTR
*
*           Since all processes must end up with the same value
*           of NZ, we first compute the minimum of LALLWORK
*
            CALL IGAMN2D( ICTXT, 'A', ' ', 1, 1, LALLWORK, 1, 1, 1, -1,
     $                    -1, -1 )
*
            MAXEIGS = DESCZ( N_ )
*
            DO 50 NZ = MIN( MAXEIGS, M ), 0, -1
               MQ0 = NUMROC( NZ, NB, 0, 0, NPCOL )
               SIZESTEIN = ICEIL( NZ, NPROCS )*N + MAX( 5*N, NP0*MQ0 )
               SIZEORMTR = MAX( ( NB*( NB-1 ) ) / 2, ( MQ0+NP0 )*NB ) +
     $                     NB*NB
*
               SIZESYEVX = MAX( SIZESTEIN, SIZEORMTR )
               IF( SIZESYEVX.LE.LALLWORK )
     $            GO TO 60
   50       CONTINUE
   60       CONTINUE
         ELSE
            NZ = M
         END IF
         NZ = MAX( NZ, 0 )
         IF( NZ.NE.M ) THEN
            INFO = INFO + IERRSPC
*
            DO 70 I = 1, M
               IFAIL( I ) = 0
   70       CONTINUE
*
*     The following code handles a rare special case
*       - NZ .NE. M means that we don't have enough room to store
*         all the vectors.
*       - NSPLIT .GT. 1 means that the matrix split
*     In this case, we cannot simply take the first NZ eigenvalues
*     because PSSTEBZ sorts the eigenvalues by block when
*     a split occurs.  So, we have to make another call to
*     PSSTEBZ with a new upper limit - VUU.
*
            IF( NSPLIT.GT.1 ) THEN
               CALL SLASRT( 'I', M, W, IINFO )
               NZZ = 0
               IF( NZ.GT.0 ) THEN
*
                  VUU = W( NZ ) - TEN*( EPS*ANRM+SAFMIN )
                  IF( VLL.GE.VUU ) THEN
                     NZZ = 0
                  ELSE
                     CALL PSSTEBZ( ICTXT, RANGE, ORDER, N, VLL, VUU, IL,
     $                             IU, ABSTLL, WORK( INDD2 ),
     $                             WORK( INDE2+OFFSET ), NZZ, NSPLIT, W,
     $                             IWORK( INDIBL ), IWORK( INDISP ),
     $                             WORK( INDWORK ), LLWORK, IWORK( 1 ),
     $                             ISIZESTEBZ, IINFO )
                  END IF
*
                  IF( MOD( INFO / IERREBZ, 1 ).EQ.0 ) THEN
                     IF( NZZ.GT.NZ .OR. IINFO.NE.0 ) THEN
                        INFO = INFO + IERREBZ
                     END IF
                  END IF
               END IF
               NZ = MIN( NZ, NZZ )
*
            END IF
         END IF
         CALL PSSTEIN( N, WORK( INDD2 ), WORK( INDE2+OFFSET ), NZ, W,
     $                 IWORK( INDIBL ), IWORK( INDISP ), ORFAC, Z, IZ,
     $                 JZ, DESCZ, WORK( INDWORK ), LALLWORK, IWORK( 1 ),
     $                 ISIZESTEIN, IFAIL, ICLUSTR, GAP, IINFO )
*
         IF( IINFO.GE.NZ+1 )
     $      INFO = INFO + IERRCLS
         IF( MOD( IINFO, NZ+1 ).NE.0 )
     $      INFO = INFO + IERREIN
*
*     Z = Q * Z
*
*
         IF( NZ.GT.0 ) THEN
            CALL PSORMTR( 'L', UPLO, 'N', N, NZ, A, IA, JA, DESCA,
     $                    WORK( INDTAU ), Z, IZ, JZ, DESCZ,
     $                    WORK( INDWORK ), LLWORK, IINFO )
         END IF
*
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         CALL SSCAL( M, ONE / SIGMA, W, 1 )
      END IF
*
      WORK( 1 ) = REAL( LWOPT )
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of PSSYEVX
*
      END
