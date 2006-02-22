      SUBROUTINE PZHEGVX( IBTYPE, JOBZ, RANGE, UPLO, N, A, IA, JA,
     $                    DESCA, B, IB, JB, DESCB, VL, VU, IL, IU,
     $                    ABSTOL, M, NZ, W, ORFAC, Z, IZ, JZ, DESCZ,
     $                    WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK,
     $                    IFAIL, ICLUSTR, GAP, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     October 15, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IA, IB, IBTYPE, IL, INFO, IU, IZ, JA, JB, JZ,
     $                   LIWORK, LRWORK, LWORK, M, N, NZ
      DOUBLE PRECISION   ABSTOL, ORFAC, VL, VU
*     ..
*     .. Array Arguments ..
*
      INTEGER            DESCA( * ), DESCB( * ), DESCZ( * ),
     $                   ICLUSTR( * ), IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   GAP( * ), RWORK( * ), W( * )
      COMPLEX*16         A( * ), B( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*
*  =======
*
*  PZHEGVX computes all the eigenvalues, and optionally,
*  the eigenvectors
*  of a complex generalized Hermitian-definite eigenproblem, of the form
*  sub( A )*x=(lambda)*sub( B )*x,  sub( A )*sub( B )x=(lambda)*x,  or
*  sub( B )*sub( A )*x=(lambda)*x.
*  Here sub( A ) denoting A( IA:IA+N-1, JA:JA+N-1 ) is assumed to be
*  Hermitian, and sub( B ) denoting B( IB:IB+N-1, JB:JB+N-1 ) is assumed
*  to be Hermitian positive definite.
*
*  Notes
*  =====
*
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
*
*  Arguments
*  =========
*
*  IBTYPE   (global input) INTEGER
*          Specifies the problem type to be solved:
*          = 1:  sub( A )*x = (lambda)*sub( B )*x
*          = 2:  sub( A )*sub( B )*x = (lambda)*x
*          = 3:  sub( B )*sub( A )*x = (lambda)*x
*
*  JOBZ    (global input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (global input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the interval [VL,VU] will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  UPLO    (global input) CHARACTER*1
*          = 'U':  Upper triangles of sub( A ) and sub( B ) are stored;
*          = 'L':  Lower triangles of sub( A ) and sub( B ) are stored.
*
*  N       (global input) INTEGER
*          The order of the matrices sub( A ) and sub( B ).  N >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          N-by-N Hermitian distributed matrix sub( A ). If UPLO = 'U',
*          the leading N-by-N upper triangular part of sub( A ) contains
*          the upper triangular part of the matrix.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of sub( A ) contains
*          the lower triangular part of the matrix.
*
*          On exit, if JOBZ = 'V', then if INFO = 0, sub( A ) contains
*          the distributed matrix Z of eigenvectors.  The eigenvectors
*          are normalized as follows:
*          if IBTYPE = 1 or 2, Z**H*sub( B )*Z = I;
*          if IBTYPE = 3, Z**H*inv( sub( B ) )*Z = I.
*          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
*          or the lower triangle (if UPLO='L') of sub( A ), including
*          the diagonal, is destroyed.
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*          If DESCA( CTXT_ ) is incorrect, PZHEGVX cannot guarantee
*          correct error reporting.
*
*  B       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_B, LOCc(JB+N-1)).
*          On entry, this array contains the local pieces of the
*          N-by-N Hermitian distributed matrix sub( B ). If UPLO = 'U',
*          the leading N-by-N upper triangular part of sub( B ) contains
*          the upper triangular part of the matrix.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of sub( B ) contains
*          the lower triangular part of the matrix.
*
*          On exit, if INFO <= N, the part of sub( B ) containing the
*          matrix is overwritten by the triangular factor U or L from
*          the Cholesky factorization sub( B ) = U**H*U or
*          sub( B ) = L*L**H.
*
*  IB      (global input) INTEGER
*          The row index in the global array B indicating the first
*          row of sub( B ).
*
*  JB      (global input) INTEGER
*          The column index in the global array B indicating the
*          first column of sub( B ).
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
*          DESCB( CTXT_ ) must equal DESCA( CTXT_ )
*
*  VL      (global input) DOUBLE PRECISION
*          If RANGE='V', the lower bound of the interval to be searched
*          for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
*
*  VU      (global input) DOUBLE PRECISION
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
*  ABSTOL  (global input) DOUBLE PRECISION
*          If JOBZ='V', setting ABSTOL to PDLAMCH( CONTEXT, 'U') yields
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
*          set to twice the underflow threshold 2*PDLAMCH('S') not zero.
*          If this routine returns with ((MOD(INFO,2).NE.0) .OR.
*          (MOD(INFO/8,2).NE.0)), indicating that some eigenvalues or
*          eigenvectors did not converge, try setting ABSTOL to
*          2*PDLAMCH('S').
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
*          insufficient space and PZHEGVX is not able to detect this
*          before beginning computation.  To get all the eigenvectors
*          requested, the user must supply both sufficient
*          space to hold the eigenvectors in Z (M .LE. DESCZ(N_))
*          and sufficient workspace to compute them.  (See LWORK below.)
*          PZHEGVX is always able to detect insufficient space without
*          computation unless RANGE .EQ. 'V'.
*
*  W       (global output) DOUBLE PRECISION array, dimension (N)
*          On normal exit, the first M entries contain the selected
*          eigenvalues in ascending order.
*
*  ORFAC   (global input) DOUBLE PRECISION
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
*  Z       (local output) COMPLEX*16 array,
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
*          The row index in the global array Z indicating the first
*          row of sub( Z ).
*
*  JZ      (global input) INTEGER
*          The column index in the global array Z indicating the
*          first column of sub( Z ).
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*          DESCZ( CTXT_ ) must equal DESCA( CTXT_ )
*
*  WORK    (local workspace/output) COMPLEX*16 array,
*             dimension (LWORK)
*          WORK(1) returns the optimal workspace.
*
*  LWORK   (local input) INTEGER
*          Size of WORK array.  If only eigenvalues are requested:
*            LWORK >= N + MAX( NB * ( NP0 + 1 ), 3 )
*          If eigenvectors are requested:
*            LWORK >= N + ( NP0 + MQ0 + NB ) * NB
*          with NQ0 = NUMROC( NN, NB, 0, 0, NPCOL ).
*
*          For optimal performance, greater workspace is needed, i.e.
*            LWORK >= MAX( LWORK, N + NHETRD_LWOPT,
*              NHEGST_LWOPT )
*          Where LWORK is as defined above, and
*            NHETRD_LWORK = 2*( ANB+1 )*( 4*NPS+2 ) +
*              ( NPS + 1 ) * NPS
*            NHEGST_LWOPT =  2*NP0*NB + NQ0*NB + NB*NB
*
*            NB = DESCA( MB_ )
*            NP0 = NUMROC( N, NB, 0, 0, NPROW )
*            NQ0 = NUMROC( N, NB, 0, 0, NPCOL )
*            ICTXT = DESCA( CTXT_ )
*            ANB = PJLAENV( ICTXT, 3, 'PZHETTRD', 'L', 0, 0, 0, 0 )
*            SQNPC = SQRT( DBLE( NPROW * NPCOL ) )
*            NPS = MAX( NUMROC( N, 1, 0, 0, SQNPC ), 2*ANB )
*
*            NUMROC is a ScaLAPACK tool functions;
*            PJLAENV is a ScaLAPACK envionmental inquiry function
*            MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*            the subroutine BLACS_GRIDINFO.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the optimal
*          size for all work arrays. Each of these values is returned
*          in the first entry of the correspondingwork array, and no
*          error message is issued by PXERBLA.
*
*  RWORK   (local workspace/output) DOUBLE PRECISION array,
*             dimension max(3,LRWORK)
*          On return, RWORK(1) contains the amount of workspace
*             required for optimal efficiency
*          if JOBZ='N' RWORK(1) = optimal amount of workspace
*             required to compute eigenvalues efficiently
*          if JOBZ='V' RWORK(1) = optimal amount of workspace
*             required to compute eigenvalues and eigenvectors
*             efficiently with no guarantee on orthogonality.
*             If RANGE='V', it is assumed that all eigenvectors
*             may be required when computing optimal workspace.
*
*  LRWORK  (local input) INTEGER
*          Size of RWORK
*          See below for definitions of variables used to define LRWORK.
*          If no eigenvectors are requested (JOBZ = 'N') then
*             LRWORK >= 5 * NN + 4 * N
*          If eigenvectors are requested (JOBZ = 'V' ) then
*             the amount of workspace required to guarantee that all
*             eigenvectors are computed is:
*             LRWORK >= 4*N + MAX( 5*NN, NP0 * MQ0 ) +
*               ICEIL( NEIG, NPROW*NPCOL)*NN
*
*             The computed eigenvectors may not be orthogonal if the
*             minimal workspace is supplied and ORFAC is too small.
*             If you want to guarantee orthogonality (at the cost
*             of potentially poor performance) you should add
*             the following to LRWORK:
*                (CLUSTERSIZE-1)*N
*             where CLUSTERSIZE is the number of eigenvalues in the
*             largest cluster, where a cluster is defined as a set of
*             close eigenvalues: { W(K),...,W(K+CLUSTERSIZE-1) |
*                                  W(J+1) <= W(J) + ORFAC*2*norm(A) }
*          Variable definitions:
*             NEIG = number of eigenvectors requested
*             NB = DESCA( MB_ ) = DESCA( NB_ ) = DESCZ( MB_ ) =
*                  DESCZ( NB_ )
*             NN = MAX( N, NB, 2 )
*             DESCA( RSRC_ ) = DESCA( NB_ ) = DESCZ( RSRC_ ) =
*                              DESCZ( CSRC_ ) = 0
*             NP0 = NUMROC( NN, NB, 0, 0, NPROW )
*             MQ0 = NUMROC( MAX( NEIG, NB, 2 ), NB, 0, 0, NPCOL )
*             ICEIL( X, Y ) is a ScaLAPACK function returning
*             ceiling(X/Y)
*
*          When LRWORK is too small:
*             If LRWORK is too small to guarantee orthogonality,
*             PZHEGVX attempts to maintain orthogonality in
*             the clusters with the smallest
*             spacing between the eigenvalues.
*             If LRWORK is too small to compute all the eigenvectors
*             requested, no computation is performed and INFO=-25
*             is returned.  Note that when RANGE='V', PZHEGVX does
*             not know how many eigenvectors are requested until
*             the eigenvalues are computed.  Therefore, when RANGE='V'
*             and as long as LRWORK is large enough to allow PZHEGVX to
*             compute the eigenvalues, PZHEGVX will compute the
*             eigenvalues and as many eigenvectors as it can.
*
*          Relationship between workspace, orthogonality & performance:
*             If CLUSTERSIZE >= N/SQRT(NPROW*NPCOL), then providing
*             enough space to compute all the eigenvectors
*             orthogonally will cause serious degradation in
*             performance. In the limit (i.e. CLUSTERSIZE = N-1)
*             PZSTEIN will perform no better than ZSTEIN on 1 processor.
*             For CLUSTERSIZE = N/SQRT(NPROW*NPCOL) reorthogonalizing
*             all eigenvectors will increase the total execution time
*             by a factor of 2 or more.
*             For CLUSTERSIZE > N/SQRT(NPROW*NPCOL) execution time will
*             grow as the square of the cluster size, all other factors
*             remaining equal and assuming enough workspace.  Less
*             workspace means less reorthogonalization but faster
*             execution.
*
*          If LRWORK = -1, then LRWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
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
*
*          If LIWORK = -1, then LIWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  IFAIL   (output) INTEGER array, dimension (N)
*          IFAIL provides additional information when INFO .NE. 0
*          If (MOD(INFO/16,2).NE.0) then IFAIL(1) indicates the order of
*          the smallest minor which is not positive definite.
*          If (MOD(INFO,2).NE.0) on exit, then IFAIL contains the
*          indices of the eigenvectors that failed to converge.
*
*          If neither of the above error conditions hold and JOBZ = 'V',
*          then the first M elements of IFAIL are set to zero.
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
*  GAP     (global output) DOUBLE PRECISION array,
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
*                  in IFAIL.  Send e-mail to scalapack@cs.utk.edu
*                if (MOD(INFO/2,2).NE.0),then eigenvectors corresponding
*                  to one or more clusters of eigenvalues could not be
*                  reorthogonalized because of insufficient workspace.
*                  The indices of the clusters are stored in the array
*                  ICLUSTR.
*                if (MOD(INFO/4,2).NE.0), then space limit prevented
*                  PZHEGVX from computing all of the eigenvectors
*                  between VL and VU.  The number of eigenvectors
*                  computed is returned in NZ.
*                if (MOD(INFO/8,2).NE.0), then PZSTEBZ failed to
*                  compute eigenvalues.
*                  Send e-mail to scalapack@cs.utk.edu
*                if (MOD(INFO/16,2).NE.0), then B was not positive
*                  definite.  IFAIL(1) indicates the order of
*                  the smallest minor which is not positive definite.
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrices A(IA:*, JA:*), C(IC:IC+M-1,JC:JC+N-1),
*  and B( IB:IB+N-1, JB:JB+N-1 ) must verify some alignment properties,
*  namely the following expressions should be true:
*
*     DESCA(MB_) = DESCA(NB_)
*     IA = IB = IZ
*     JA = IB = JZ
*     DESCA(M_) = DESCB(M_) =DESCZ(M_)
*     DESCA(N_) = DESCB(N_)= DESCZ(N_)
*     DESCA(MB_) = DESCB(MB_) = DESCZ(MB_)
*     DESCA(NB_) = DESCB(NB_) = DESCZ(NB_)
*     DESCA(RSRC_) = DESCB(RSRC_) = DESCZ(RSRC_)
*     DESCA(CSRC_) = DESCB(CSRC_) = DESCZ(CSRC_)
*     MOD( IA-1, DESCA( MB_ ) ) = 0
*     MOD( JA-1, DESCA( NB_ ) ) = 0
*     MOD( IB-1, DESCB( MB_ ) ) = 0
*     MOD( JB-1, DESCB( NB_ ) ) = 0
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      COMPLEX*16         ONE
      PARAMETER          ( ONE = 1.0D+0 )
      DOUBLE PRECISION   FIVE, ZERO
      PARAMETER          ( FIVE = 5.0D+0, ZERO = 0.0D+0 )
      INTEGER            IERRNPD
      PARAMETER          ( IERRNPD = 16 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LQUERY, UPPER, VALEIG, WANTZ
      CHARACTER          TRANS
      INTEGER            ANB, IACOL, IAROW, IBCOL, IBROW, ICOFFA,
     $                   ICOFFB, ICTXT, IROFFA, IROFFB, LIWMIN, LRWMIN,
     $                   LRWOPT, LWMIN, LWOPT, MQ0, MYCOL, MYROW, NB,
     $                   NEIG, NHEGST_LWOPT, NHETRD_LWOPT, NN, NP0,
     $                   NPCOL, NPROW, NPS, NQ0, SQNPC
      DOUBLE PRECISION   EPS, SCALE
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 5 ), IDUM2( 5 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, INDXG2P, NUMROC, PJLAENV
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           LSAME, ICEIL, INDXG2P, NUMROC, PJLAENV, PDLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DGEBR2D, DGEBS2D,
     $                   DSCAL, PCHK1MAT, PCHK2MAT, PXERBLA, PZHEEVX,
     $                   PZHENGST, PZPOTRF, PZTRMM, PZTRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, ICHAR, INT, MAX, MIN, MOD,
     $                   SQRT
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 900+CTXT_ )
      ELSE IF( DESCA( CTXT_ ).NE.DESCB( CTXT_ ) ) THEN
         INFO = -( 1300+CTXT_ )
      ELSE IF( DESCA( CTXT_ ).NE.DESCZ( CTXT_ ) ) THEN
         INFO = -( 2600+CTXT_ )
      ELSE
*
*     Get machine constants.
*
         EPS = PDLAMCH( DESCA( CTXT_ ), 'Precision' )
*
         WANTZ = LSAME( JOBZ, 'V' )
         UPPER = LSAME( UPLO, 'U' )
         ALLEIG = LSAME( RANGE, 'A' )
         VALEIG = LSAME( RANGE, 'V' )
         INDEIG = LSAME( RANGE, 'I' )
         CALL CHK1MAT( N, 4, N, 4, IA, JA, DESCA, 9, INFO )
         CALL CHK1MAT( N, 4, N, 4, IB, JB, DESCB, 13, INFO )
         CALL CHK1MAT( N, 4, N, 4, IZ, JZ, DESCZ, 26, INFO )
         IF( INFO.EQ.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
               RWORK( 1 ) = ABSTOL
               IF( VALEIG ) THEN
                  RWORK( 2 ) = VL
                  RWORK( 3 ) = VU
               ELSE
                  RWORK( 2 ) = ZERO
                  RWORK( 3 ) = ZERO
               END IF
               CALL DGEBS2D( DESCA( CTXT_ ), 'ALL', ' ', 3, 1, RWORK,
     $                       3 )
            ELSE
               CALL DGEBR2D( DESCA( CTXT_ ), 'ALL', ' ', 3, 1, RWORK, 3,
     $                       0, 0 )
            END IF
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $              NPROW )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $              NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $              NPCOL )
            IBCOL = INDXG2P( JB, DESCB( NB_ ), MYCOL, DESCB( CSRC_ ),
     $              NPCOL )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            ICOFFB = MOD( JB-1, DESCB( NB_ ) )
*
*     Compute the total amount of space needed
*
            LQUERY = .FALSE.
            IF( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
     $         LQUERY = .TRUE.
*
            LIWMIN = 6*MAX( N, ( NPROW*NPCOL )+1, 4 )
*
            NB = DESCA( MB_ )
            NN = MAX( N, NB, 2 )
            NP0 = NUMROC( NN, NB, 0, 0, NPROW )
*
            IF( ( .NOT.WANTZ ) .OR. ( VALEIG .AND. ( .NOT.LQUERY ) ) )
     $           THEN
               LWMIN = N + MAX( NB*( NP0+1 ), 3 )
               LWOPT = LWMIN
               LRWMIN = 5*NN + 4*N
               IF( WANTZ ) THEN
                  MQ0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPCOL )
                  LRWOPT = 4*N + MAX( 5*NN, NP0*MQ0 )
               ELSE
                  LRWOPT = LRWMIN
               END IF
               NEIG = 0
            ELSE
               IF( ALLEIG .OR. VALEIG ) THEN
                  NEIG = N
               ELSE IF( INDEIG ) THEN
                  NEIG = IU - IL + 1
               END IF
               MQ0 = NUMROC( MAX( NEIG, NB, 2 ), NB, 0, 0, NPCOL )
               LWMIN = N + ( NP0+MQ0+NB )*NB
               LWOPT = LWMIN
               LRWMIN = 4*N + MAX( 5*NN, NP0*MQ0 ) +
     $                  ICEIL( NEIG, NPROW*NPCOL )*NN
               LRWOPT = LRWMIN
*
            END IF
*
*           Compute how much workspace is needed to use the
*           new TRD and GST algorithms
*
            ANB = PJLAENV( ICTXT, 3, 'PZHETTRD', 'L', 0, 0, 0, 0 )
            SQNPC = INT( SQRT( DBLE( NPROW*NPCOL ) ) )
            NPS = MAX( NUMROC( N, 1, 0, 0, SQNPC ), 2*ANB )
            NHETRD_LWOPT = 2*( ANB+1 )*( 4*NPS+2 ) + ( NPS+4 )*NPS
            NB = DESCA( MB_ )
            NP0 = NUMROC( N, NB, 0, 0, NPROW )
            NQ0 = NUMROC( N, NB, 0, 0, NPCOL )
            NHEGST_LWOPT = 2*NP0*NB + NQ0*NB + NB*NB
            LWOPT = MAX( LWOPT, N+NHETRD_LWOPT, NHEGST_LWOPT )
*
*       Version 1.0 Limitations
*
            IF( IBTYPE.LT.1 .OR. IBTYPE.GT.3 ) THEN
               INFO = -1
            ELSE IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
               INFO = -2
            ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
               INFO = -3
            ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -4
            ELSE IF( N.LT.0 ) THEN
               INFO = -5
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -7
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -8
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 900+NB_ )
            ELSE IF( DESCA( M_ ).NE.DESCB( M_ ) ) THEN
               INFO = -( 1300+M_ )
            ELSE IF( DESCA( N_ ).NE.DESCB( N_ ) ) THEN
               INFO = -( 1300+N_ )
            ELSE IF( DESCA( MB_ ).NE.DESCB( MB_ ) ) THEN
               INFO = -( 1300+MB_ )
            ELSE IF( DESCA( NB_ ).NE.DESCB( NB_ ) ) THEN
               INFO = -( 1300+NB_ )
            ELSE IF( DESCA( RSRC_ ).NE.DESCB( RSRC_ ) ) THEN
               INFO = -( 1300+RSRC_ )
            ELSE IF( DESCA( CSRC_ ).NE.DESCB( CSRC_ ) ) THEN
               INFO = -( 1300+CSRC_ )
            ELSE IF( DESCA( CTXT_ ).NE.DESCB( CTXT_ ) ) THEN
               INFO = -( 1300+CTXT_ )
            ELSE IF( DESCA( M_ ).NE.DESCZ( M_ ) ) THEN
               INFO = -( 2200+M_ )
            ELSE IF( DESCA( N_ ).NE.DESCZ( N_ ) ) THEN
               INFO = -( 2200+N_ )
            ELSE IF( DESCA( MB_ ).NE.DESCZ( MB_ ) ) THEN
               INFO = -( 2200+MB_ )
            ELSE IF( DESCA( NB_ ).NE.DESCZ( NB_ ) ) THEN
               INFO = -( 2200+NB_ )
            ELSE IF( DESCA( RSRC_ ).NE.DESCZ( RSRC_ ) ) THEN
               INFO = -( 2200+RSRC_ )
            ELSE IF( DESCA( CSRC_ ).NE.DESCZ( CSRC_ ) ) THEN
               INFO = -( 2200+CSRC_ )
            ELSE IF( DESCA( CTXT_ ).NE.DESCZ( CTXT_ ) ) THEN
               INFO = -( 2200+CTXT_ )
            ELSE IF( IROFFB.NE.0 .OR. IBROW.NE.IAROW ) THEN
               INFO = -11
            ELSE IF( ICOFFB.NE.0 .OR. IBCOL.NE.IACOL ) THEN
               INFO = -12
            ELSE IF( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) THEN
               INFO = -15
            ELSE IF( INDEIG .AND. ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) )
     $                THEN
               INFO = -16
            ELSE IF( INDEIG .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) )
     $                THEN
               INFO = -17
            ELSE IF( VALEIG .AND. ( ABS( RWORK( 2 )-VL ).GT.FIVE*EPS*
     $               ABS( VL ) ) ) THEN
               INFO = -14
            ELSE IF( VALEIG .AND. ( ABS( RWORK( 3 )-VU ).GT.FIVE*EPS*
     $               ABS( VU ) ) ) THEN
               INFO = -15
            ELSE IF( ABS( RWORK( 1 )-ABSTOL ).GT.FIVE*EPS*
     $               ABS( ABSTOL ) ) THEN
               INFO = -18
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -28
            ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -30
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -32
            END IF
         END IF
         IDUM1( 1 ) = IBTYPE
         IDUM2( 1 ) = 1
         IF( WANTZ ) THEN
            IDUM1( 2 ) = ICHAR( 'V' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'N' )
         END IF
         IDUM2( 2 ) = 2
         IF( UPPER ) THEN
            IDUM1( 3 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 3 ) = ICHAR( 'L' )
         END IF
         IDUM2( 3 ) = 3
         IF( ALLEIG ) THEN
            IDUM1( 4 ) = ICHAR( 'A' )
         ELSE IF( INDEIG ) THEN
            IDUM1( 4 ) = ICHAR( 'I' )
         ELSE
            IDUM1( 4 ) = ICHAR( 'V' )
         END IF
         IDUM2( 4 ) = 4
         IF( LQUERY ) THEN
            IDUM1( 5 ) = -1
         ELSE
            IDUM1( 5 ) = 1
         END IF
         IDUM2( 5 ) = 5
         CALL PCHK2MAT( N, 4, N, 4, IA, JA, DESCA, 9, N, 4, N, 4, IB,
     $                  JB, DESCB, 13, 5, IDUM1, IDUM2, INFO )
         CALL PCHK1MAT( N, 4, N, 4, IZ, JZ, DESCZ, 26, 0, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IWORK( 1 ) = LIWMIN
      WORK( 1 ) = DCMPLX( DBLE( LWOPT ) )
      RWORK( 1 ) = DBLE( LRWOPT )
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZHEGVX ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Form a Cholesky factorization of sub( B ).
*
      CALL PZPOTRF( UPLO, N, B, IB, JB, DESCB, INFO )
      IF( INFO.NE.0 ) THEN
         IWORK( 1 ) = LIWMIN
         WORK( 1 ) = DCMPLX( DBLE( LWOPT ) )
         RWORK( 1 ) = DBLE( LRWOPT )
         IFAIL( 1 ) = INFO
         INFO = IERRNPD
         RETURN
      END IF
*
*     Transform problem to standard eigenvalue problem and solve.
*
      CALL PZHENGST( IBTYPE, UPLO, N, A, IA, JA, DESCA, B, IB, JB,
     $               DESCB, SCALE, WORK, LWORK, INFO )
      CALL PZHEEVX( JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, VL, VU, IL,
     $              IU, ABSTOL, M, NZ, W, ORFAC, Z, IZ, JZ, DESCZ, WORK,
     $              LWORK, RWORK, LRWORK, IWORK, LIWORK, IFAIL, ICLUSTR,
     $              GAP, INFO )
*
      IF( WANTZ ) THEN
*
*        Backtransform eigenvectors to the original problem.
*
         NEIG = M
         IF( IBTYPE.EQ.1 .OR. IBTYPE.EQ.2 ) THEN
*
*           For sub( A )*x=(lambda)*sub( B )*x and
*           sub( A )*sub( B )*x=(lambda)*x; backtransform eigenvectors:
*           x = inv(L)'*y or inv(U)*y
*
            IF( UPPER ) THEN
               TRANS = 'N'
            ELSE
               TRANS = 'C'
            END IF
*
            CALL PZTRSM( 'Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE,
     $                   B, IB, JB, DESCB, Z, IZ, JZ, DESCZ )
*
         ELSE IF( IBTYPE.EQ.3 ) THEN
*
*           For sub( B )*sub( A )*x=(lambda)*x;
*           backtransform eigenvectors: x = L*y or U'*y
*
            IF( UPPER ) THEN
               TRANS = 'C'
            ELSE
               TRANS = 'N'
            END IF
*
            CALL PZTRMM( 'Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE,
     $                   B, IB, JB, DESCB, Z, IZ, JZ, DESCZ )
         END IF
      END IF
*
      IF( SCALE.NE.ONE ) THEN
         CALL DSCAL( N, SCALE, W, 1 )
      END IF
*
      IWORK( 1 ) = LIWMIN
      WORK( 1 ) = DCMPLX( DBLE( LWOPT ) )
      RWORK( 1 ) = DBLE( LRWOPT )
      RETURN
*
*     End of PZHEGVX
*
      END
