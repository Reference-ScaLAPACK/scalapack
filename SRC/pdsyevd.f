      SUBROUTINE PDSYEVD( JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ,
     $                    DESCZ, WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 14, 2000
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            IA, INFO, IZ, JA, JZ, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * ), IWORK( * )
      DOUBLE PRECISION   A( * ), W( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PDSYEVD computes  all the eigenvalues and eigenvectors
*  of a real symmetric matrix A by calling the recommended sequence
*  of ScaLAPACK routines.
*
*  In its present form, PDSYEVD assumes a homogeneous system and makes
*  no checks for consistency of the eigenvalues or eigenvectors across
*  the different processes.  Because of this, it is possible that a
*  heterogeneous system may return incorrect results without any error
*  messages.
*
*  Arguments
*  =========
*
*     NP = the number of rows local to a given process.
*     NQ = the number of columns local to a given process.
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;     (NOT IMPLEMENTED YET)
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/workspace) block cyclic DOUBLE PRECISION array,
*          global dimension (N, N), local dimension ( LLD_A,
*          LOCc(JA+N-1) )
*          On entry, the symmetric matrix A.  If UPLO = 'U', only the
*          upper triangular part of A is used to define the elements of
*          the symmetric matrix.  If UPLO = 'L', only the lower
*          triangular part of A is used to define the elements of the
*          symmetric matrix.
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
*
*  W       (global output) DOUBLE PRECISION array, dimension (N)
*          If INFO=0, the eigenvalues in ascending order.
*
*  Z       (local output) DOUBLE PRECISION array,
*          global dimension (N, N),
*          local dimension ( LLD_Z, LOCc(JZ+N-1) )
*          Z contains the orthonormal eigenvectors
*          of the symmetric matrix A.
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
*  WORK    (local workspace/output) DOUBLE PRECISION array,
*          dimension (LWORK)
*          On output, WORK(1) returns the workspace required.
*
*  LWORK   (local input) INTEGER
*          LWORK >= MAX( 1+6*N+2*NP*NQ, TRILWMIN ) + 2*N
*          TRILWMIN = 3*N + MAX( NB*( NP+1 ), 3*NB )
*          NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
*          NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
*
*          If LWORK = -1, the LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          size for the WORK array.  The required workspace is returned
*          as the first element of WORK and no error message is issued
*          by PXERBLA.
*
*  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
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
*          > 0:  The algorithm failed to compute the INFO/(N+1) th
*                eigenvalue while working on the submatrix lying in
*                global rows and columns mod(INFO,N+1).
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
*
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            IACOL, IAROW, ICOFFA, ICOFFZ, ICTXT, IINFO,
     $                   INDD, INDE, INDE2, INDTAU, INDWORK, INDWORK2,
     $                   IROFFA, IROFFZ, ISCALE, LIWMIN, LLWORK,
     $                   LLWORK2, LWMIN, MYCOL, MYROW, NB, NP, NPCOL,
     $                   NPROW, NQ, OFFSET, TRILWMIN
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. Local Arrays ..
*     ..
      INTEGER            IDUM1( 2 ), IDUM2( 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      DOUBLE PRECISION   PDLAMCH, PDLANSY
      EXTERNAL           LSAME, INDXG2P, NUMROC, PDLAMCH, PDLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DSCAL, PCHK1MAT,
     $                   PDLARED1D, PDLASCL, PDLASET, PDORMTR, PDSTEDC,
     $                   PDSYTRD, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, ICHAR, MAX, MIN, MOD, SQRT
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Quick return
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Test the input arguments.
*
      ICTXT = DESCZ( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 600+CTXT_ )
      ELSE
         CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 3, N, 3, IZ, JZ, DESCZ, 12, INFO )
         IF( INFO.EQ.0 ) THEN
            UPPER = LSAME( UPLO, 'U' )
            NB = DESCA( NB_ )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFZ = MOD( IZ-1, DESCZ( MB_ ) )
            ICOFFZ = MOD( JZ-1, DESCZ( NB_ ) )
            IAROW = INDXG2P( IA, NB, MYROW, DESCA( RSRC_ ), NPROW )
            IACOL = INDXG2P( JA, NB, MYCOL, DESCA( CSRC_ ), NPCOL )
            NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
            NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
*
            LQUERY = ( LWORK.EQ.-1 )
            TRILWMIN = 3*N + MAX( NB*( NP+1 ), 3*NB )
            LWMIN = MAX( 1+6*N+2*NP*NQ, TRILWMIN ) + 2*N
            LIWMIN = 7*N + 8*NPCOL + 2
            WORK( 1 ) = DBLE( LWMIN )
            IWORK( 1 ) = LIWMIN
            IF( .NOT.LSAME( JOBZ, 'V' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -2
            ELSE IF( IROFFA.NE.ICOFFA .OR. ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( IROFFA.NE.IROFFZ .OR. ICOFFA.NE.ICOFFZ ) THEN
               INFO = -10
            ELSE IF( DESCA( M_ ).NE.DESCZ( M_ ) ) THEN
               INFO = -( 1200+M_ )
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            ELSE IF( DESCZ( MB_ ).NE.DESCZ( NB_ ) ) THEN
               INFO = -( 1200+NB_ )
            ELSE IF( DESCA( MB_ ).NE.DESCZ( MB_ ) ) THEN
               INFO = -( 1200+MB_ )
            ELSE IF( DESCA( CTXT_ ).NE.DESCZ( CTXT_ ) ) THEN
               INFO = -( 1200+CTXT_ )
            ELSE IF( DESCA( RSRC_ ).NE.DESCZ( RSRC_ ) ) THEN
               INFO = -( 1200+RSRC_ )
            ELSE IF( DESCA( CSRC_ ).NE.DESCZ( CSRC_ ) ) THEN
               INFO = -( 1200+CSRC_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -14
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -16
            END IF
         END IF
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 2
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 14
         CALL PCHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, 2, IDUM1, IDUM2,
     $                  INFO )
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDSYEVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Set up pointers into the WORK array
*
      INDTAU = 1
      INDE = INDTAU + N
      INDD = INDE + N
      INDE2 = INDD + N
      INDWORK = INDE2 + N
      LLWORK = LWORK - INDWORK + 1
      INDWORK2 = INDD
      LLWORK2 = LWORK - INDWORK2 + 1
*
*     Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
      SAFMIN = PDLAMCH( DESCA( CTXT_ ), 'Safe minimum' )
      EPS = PDLAMCH( DESCA( CTXT_ ), 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
      ANRM = PDLANSY( 'M', UPLO, N, A, IA, JA, DESCA, WORK( INDWORK ) )
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
         CALL PDLASCL( UPLO, ONE, SIGMA, N, N, A, IA, JA, DESCA, IINFO )
      END IF
*
*     Reduce symmetric matrix to tridiagonal form.
*
*
      CALL PDSYTRD( UPLO, N, A, IA, JA, DESCA, WORK( INDD ),
     $              WORK( INDE2 ), WORK( INDTAU ), WORK( INDWORK ),
     $              LLWORK, IINFO )
*
*     Copy the values of D, E to all processes.
*
      CALL PDLARED1D( N, IA, JA, DESCA, WORK( INDD ), W,
     $                WORK( INDWORK ), LLWORK )
*
      CALL PDLARED1D( N, IA, JA, DESCA, WORK( INDE2 ), WORK( INDE ),
     $                WORK( INDWORK ), LLWORK )
*
      CALL PDLASET( 'Full', N, N, ZERO, ONE, Z, 1, 1, DESCZ )
*
      IF( UPPER ) THEN
         OFFSET = 1
      ELSE
         OFFSET = 0
      END IF
      CALL PDSTEDC( 'I', N, W, WORK( INDE+OFFSET ), Z, IZ, JZ, DESCZ,
     $              WORK( INDWORK2 ), LLWORK2, IWORK, LIWORK, INFO )
*
      CALL PDORMTR( 'L', UPLO, 'N', N, N, A, IA, JA, DESCA,
     $              WORK( INDTAU ), Z, IZ, JZ, DESCZ, WORK( INDWORK2 ),
     $              LLWORK2, IINFO )
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         CALL DSCAL( N, ONE / SIGMA, W, 1 )
      END IF
*
      RETURN
*
*     End of PDSYEVD
*
      END
