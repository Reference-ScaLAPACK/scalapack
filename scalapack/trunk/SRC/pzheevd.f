      SUBROUTINE PZHEEVD( JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ,
     $                    DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK,
     $                    LIWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 25, 2002
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            IA, INFO, IZ, JA, JZ, LIWORK, LRWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * ), IWORK( * )
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( * ), WORK( * ), Z( * )
*     
*
*  Purpose
*  =======
*
*  PZHEEVD computes all the eigenvalues and eigenvectors of a Hermitian
*  matrix A by using a divide and conquer algorithm.
*
*  Arguments
*  =========
*
*     NP = the number of rows local to a given process.
*     NQ = the number of columns local to a given process.
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;    (NOT IMPLEMENTED YET)
*          = 'V':  Compute eigenvalues and eigenvectors.
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
*  A       (local input/workspace) block cyclic COMPLEX*16 array,
*          global dimension (N, N), local dimension ( LLD_A,
*          LOCc(JA+N-1) )
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
*          If DESCA( CTXT_ ) is incorrect, PZHEEV cannot guarantee
*          correct error reporting.
*
*  W       (global output) DOUBLE PRECISION array, dimension (N)
*          If INFO=0, the eigenvalues in ascending order.
*
*  Z       (local output) COMPLEX*16 array,
*          global dimension (N, N),
*          local dimension ( LLD_Z, LOCc(JZ+N-1) )
*          Z contains the orthonormal eigenvectors of the matrix A.
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
*  WORK    (local workspace/output) COMPLEX*16 array,
*          dimension (LWORK)
*          On output, WORK(1) returns the workspace needed for the
*          computation.
*
*  LWORK   (local input) INTEGER
*          If eigenvectors are requested:
*            LWORK = N + ( NP0 + MQ0 + NB ) * NB,
*          with  NP0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPROW )
*                MQ0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPCOL )
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine calculates the size for all
*          work arrays. Each of these values is returned in the first
*          entry of the corresponding work array, and no error message
*          is issued by PXERBLA.
*
*  RWORK   (local workspace/output) DOUBLE PRECISION array,
*          dimension (LRWORK)
*          On output RWORK(1) returns the real workspace needed to
*          guarantee completion.  If the input parameters are incorrect,
*          RWORK(1) may also be incorrect.
*
*  LRWORK  (local input) INTEGER
*          Size of RWORK array.
*          LRWORK >= 1 + 9*N + 3*NP*NQ,
*          NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
*          NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
*
*  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
*          On output IWORK(1) returns the integer workspace needed.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          LIWORK = 7*N + 8*NPCOL + 2
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = 1 through N, the i(th) eigenvalue did not
*                converge in PDLAED3.
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrices sub( A ), sub( Z ) must verify
*  some alignment properties, namely the following expression
*  should be true:
*  ( MB_A.EQ.NB_A.EQ.MB_Z.EQ.NB_Z .AND. IROFFA.EQ.ICOFFA .AND.
*    IROFFA.EQ.0 .AND.IROFFA.EQ.IROFFZ. AND. IAROW.EQ.IZROW)
*    with IROFFA = MOD( IA-1, MB_A )
*     and ICOFFA = MOD( JA-1, NB_A ).
*
*  Further Details
*  ======= =======
*
*  Contributed by Francoise Tisseur, University of Manchester.
*
*  Reference:  F. Tisseur and J. Dongarra, "A Parallel Divide and
*              Conquer Algorithm for the Symmetric Eigenvalue Problem
*              on Distributed Memory Architectures",
*              SIAM J. Sci. Comput., 6:20 (1999), pp. 2223--2236.
*              (see also LAPACK Working Note 132)
*                http://www.netlib.org/lapack/lawns/lawn132.ps
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION               ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY
      INTEGER            CSRC_A, I, IACOL, IAROW, ICOFFA, IINFO, IIZ,
     $                   INDD, INDE, INDE2, INDRWORK, INDTAU, INDWORK,
     $                   INDZ, IPR, IPZ, IROFFA, IROFFZ, ISCALE, IZCOL,
     $                   IZROW, J, JJZ, LDR, LDZ, LIWMIN, LLRWORK,
     $                   LLWORK, LRWMIN, LWMIN, MB_A, MYCOL, MYROW, NB,
     $                   NB_A, NN, NP0, NPCOL, NPROW, NQ, NQ0, OFFSET,
     $                   RSRC_A
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER            DESCRZ( 9 ), IDUM1( 2 ), IDUM2( 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2L, INDXG2P, NUMROC
      DOUBLE PRECISION   PZLANHE, PDLAMCH
      EXTERNAL           LSAME, INDXG2L, INDXG2P, NUMROC, PZLANHE,
     $                   PDLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCINIT, INFOG2L,
     $                   PZELGET, PZHETRD, PCHK2MAT, PZLASCL, PZLASET,
     $                   PZUNMTR, PDLARED1D, PDLASET, PDSTEDC, PXERBLA,
     $                   DSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, ICHAR, MAX, MIN, MOD, DBLE, SQRT
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      INFO = 0
*
*     Quick return
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Test the input arguments.
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 700+CTXT_ )
      ELSE 
         CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 3, N, 3, IZ, JZ, DESCZ, 12, INFO )
         IF( INFO.EQ.0 ) THEN
            LOWER = LSAME( UPLO, 'L' )
            NB_A = DESCA( NB_ )
            MB_A = DESCA( MB_ )
            NB = NB_A
            RSRC_A = DESCA( RSRC_ )
            CSRC_A = DESCA( CSRC_ )
            IROFFA = MOD( IA-1, MB_A )
            ICOFFA = MOD( JA-1, NB_A )
            IAROW = INDXG2P( IA, NB_A, MYROW, RSRC_A, NPROW )
            IACOL = INDXG2P( JA, MB_A, MYCOL, CSRC_A, NPCOL )
            NP0 = NUMROC( N, NB, MYROW, IAROW, NPROW )
            NQ0 = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
            IROFFZ = MOD( IZ-1, MB_A )
            CALL INFOG2L( IZ, JZ, DESCZ, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIZ, JJZ, IZROW, IZCOL )
            LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
*
*           Compute the total amount of space needed
*
            NN = MAX( N, NB, 2 )
            NQ = NUMROC( NN, NB, 0, 0, NPCOL )
            LWMIN = N + ( NP0+NQ+NB )*NB
            LRWMIN = 1 + 9*N + 3*NP0*NQ0
            LIWMIN = 7*N + 8*NPCOL + 2
            WORK( 1 ) = DCMPLX( LWMIN )
            RWORK( 1 ) = DBLE( LRWMIN )
            IWORK( 1 ) = LIWMIN
            IF( .NOT.LSAME( JOBZ, 'V' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
               INFO = -2
            ELSE IF( LWORK.LT.LWMIN .AND. LWORK.NE.-1 ) THEN
               INFO = -14
            ELSE IF( LRWORK.LT.LRWMIN .AND. LRWORK.NE.-1 ) THEN
               INFO = -16
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -4
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            ELSE IF( IROFFA.NE.IROFFZ ) THEN
               INFO = -10
            ELSE IF( IAROW.NE.IZROW ) THEN
               INFO = -10
            ELSE IF( DESCA( M_ ).NE.DESCZ( M_ ) ) THEN
               INFO = -( 1200+M_ )
            ELSE IF( DESCA( N_ ).NE.DESCZ( N_ ) ) THEN
               INFO = -( 1200+N_ )
            ELSE IF( DESCA( MB_ ).NE.DESCZ( MB_ ) ) THEN
               INFO = -( 1200+MB_ )
            ELSE IF( DESCA( NB_ ).NE.DESCZ( NB_ ) ) THEN
               INFO = -( 1200+NB_ )
            ELSE IF( DESCA( RSRC_ ).NE.DESCZ( RSRC_ ) ) THEN
               INFO = -( 1200+RSRC_ )
            ELSE IF( DESCA( CTXT_ ).NE.DESCZ( CTXT_ ) ) THEN
               INFO = -( 1200+CTXT_ )
            END IF
         END IF
         IF( LOWER ) THEN
            IDUM1( 1 ) = ICHAR( 'L' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'U' )
         END IF
         IDUM2( 1 ) = 2
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 14
         CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 7, N, 3, N, 3, IZ,
     $                  JZ, DESCZ, 11, 2, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PZHEEVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = PDLAMCH( DESCA( CTXT_ ), 'Safe minimum' )
      EPS = PDLAMCH( DESCA( CTXT_ ), 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
*     Set up pointers into the WORK array
*
      INDTAU = 1
      INDWORK = INDTAU + N
      LLWORK = LWORK - INDWORK + 1
*
*     Set up pointers into the RWORK array
*
      INDE = 1
      INDD = INDE + N
      INDE2 = INDD + N
      INDRWORK = INDE2 + N
      LLRWORK = LRWORK - INDRWORK + 1
*
*     Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
*
      ANRM = PZLANHE( 'M', UPLO, N, A, IA, JA, DESCA,
     $       RWORK( INDRWORK ) )
*
*
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
*
      IF( ISCALE.EQ.1 ) THEN
         CALL PZLASCL( UPLO, ONE, SIGMA, N, N, A, IA, JA, DESCA, IINFO )
      END IF
*
*     Reduce Hermitian matrix to tridiagonal form.
*
      CALL PZHETRD( UPLO, N, A, IA, JA, DESCA, RWORK( INDD ),
     $              RWORK( INDE2 ), WORK( INDTAU ), WORK( INDWORK ),
     $              LLWORK, IINFO )
*
*     Copy the values of D, E to all processes
*
*     Here PxLARED1D is used to redistribute the tridiagonal matrix.
*     PxLARED1D, however, doesn't yet workMx Mawith arbritary matrix
*     distributions so we have PxELGET as a backup.
*
      OFFSET = 0
      IF( IA.EQ.1 .AND. JA.EQ.1 .AND. RSRC_A.EQ.0 .AND. CSRC_A.EQ.0 )
     $     THEN
         CALL PDLARED1D( N, IA, JA, DESCA, RWORK( INDD ), W,
     $                   RWORK( INDRWORK ), LLRWORK )
*
         CALL PDLARED1D( N, IA, JA, DESCA, RWORK( INDE2 ),
     $                   RWORK( INDE ), RWORK( INDRWORK ), LLRWORK )
         IF( .NOT.LOWER )
     $      OFFSET = 1
      ELSE
         DO 10 I = 1, N
            CALL PZELGET( 'A', ' ', WORK( INDWORK ), A, I+IA-1, I+JA-1,
     $                    DESCA )
            W( I ) = DBLE( WORK( INDWORK ) )
   10    CONTINUE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 I = 1, N - 1
               CALL PZELGET( 'A', ' ', WORK( INDWORK ), A, I+IA-1, I+JA,
     $                       DESCA )
               RWORK( INDE+I-1 ) = DBLE( WORK( INDWORK ) )
   20       CONTINUE
         ELSE
            DO 30 I = 1, N - 1
               CALL PZELGET( 'A', ' ', WORK( INDWORK ), A, I+IA, I+JA-1,
     $                       DESCA )
               RWORK( INDE+I-1 ) = DBLE( WORK( INDWORK ) )
   30       CONTINUE
         END IF
      END IF
*
*     Call PDSTEDC to compute eigenvalues and eigenvectors.
*
      INDZ = INDE + N
      INDRWORK = INDZ + NP0*NQ0
      LLRWORK = LRWORK - INDRWORK + 1
      LDR = MAX( 1, NP0 )
      CALL DESCINIT( DESCRZ, DESCZ( M_ ), DESCZ( N_ ), DESCZ( MB_ ),
     $               DESCZ( NB_ ), DESCZ( RSRC_ ), DESCZ( CSRC_ ),
     $               DESCZ( CTXT_ ), LDR, INFO )
      CALL PZLASET( 'Full', N, N, CZERO, CONE, Z, IZ, JZ, DESCZ )
      CALL PDLASET( 'Full', N, N, ZERO, ONE, RWORK( INDZ ), 1, 1,
     $              DESCRZ )
      CALL PDSTEDC( 'I', N, W, RWORK( INDE+OFFSET ), RWORK( INDZ ), IZ,
     $              JZ, DESCRZ, RWORK( INDRWORK ), LLRWORK, IWORK,
     $              LIWORK, IINFO )
*
      LDZ = DESCZ( LLD_ )
      LDR = DESCRZ( LLD_ )
      IIZ = INDXG2L( IZ, NB, MYROW, MYROW, NPROW )
      JJZ = INDXG2L( JZ, NB, MYCOL, MYCOL, NPCOL )
      IPZ = IIZ + ( JJZ-1 )*LDZ
      IPR = INDZ - 1 + IIZ + ( JJZ-1 )*LDR
      DO 50 J = 0, NQ0 - 1
         DO 40 I = 0, NP0 - 1
            Z( IPZ+I+J*LDZ ) = RWORK( IPR+I+J*LDR )
   40    CONTINUE
   50 CONTINUE
*
*     Z = Q * Z
*
      CALL PZUNMTR( 'L', UPLO, 'N', N, N, A, IA, JA, DESCA,
     $              WORK( INDTAU ), Z, IZ, JZ, DESCZ, WORK( INDWORK ),
     $              LLWORK, IINFO )
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         CALL DSCAL( N, ONE / SIGMA, W, 1 )
      END IF
*
      WORK( 1 ) = DCMPLX( LWMIN )
      RWORK( 1 ) = DBLE( LRWMIN )
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of PZHEEVD
*
      END
