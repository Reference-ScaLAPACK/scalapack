*
*
      SUBROUTINE PSGSEPTST( DESCA, UPLO, N, MATTYPE, IBTYPE, SUBTESTS,
     $                      THRESH, ORDER, ABSTOL, ISEED, A, COPYA, B,
     $                      COPYB, Z, LDA, WIN, WNEW, IFAIL, ICLUSTR,
     $                      GAP, IPREPAD, IPOSTPAD, WORK, LWORK, IWORK,
     $                      LIWORK, NOUT, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          SUBTESTS, UPLO
      INTEGER            IBTYPE, INFO, IPOSTPAD, IPREPAD, LDA, LIWORK,
     $                   LWORK, MATTYPE, N, NOUT, ORDER
      REAL               ABSTOL, THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), ICLUSTR( * ), IFAIL( * ),
     $                   ISEED( 4 ), IWORK( * )
      REAL               A( LDA, * ), B( LDA, * ), COPYA( LDA, * ),
     $                   COPYB( LDA, * ), GAP( * ), WIN( * ), WNEW( * ),
     $                   WORK( * ), Z( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  PSGSEPTST builds a random matrix A, and a well conditioned
*  matrix B, runs PSSYGVX() to compute the eigenvalues
*  and eigenvectors and then calls PSSYGVCHK to compute
*  the residual.
*
*  The random matrix built depends upon the following parameters:
*     N, NB, ISEED, ORDER
*
*  Arguments
*  =========
*
*     NP = the number of rows local to a given process.
*     NQ = the number of columns local to a given process.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_
*          The array descriptor for the distributed matrices
*
*  UPLO     (global input) CHARACTER*1
*           Specifies whether the upper or lower triangular part of the
*           symmetric matrix A is stored:
*           = 'U':  Upper triangular
*           = 'L':  Lower triangular
*
*  N        (global input) INTEGER
*           Size of the matrix to be tested.  (global size)
*
*  MATTYPE  (global input) INTEGER
*           Matrix type
*  Currently, the list of possible types is:
*
*  (1)  The zero matrix.
*  (2)  The identity matrix.
*
*  (3)  A diagonal matrix with evenly spaced entries
*       1, ..., ULP  and random signs.
*       (ULP = (first number larger than 1) - 1 )
*  (4)  A diagonal matrix with geometrically spaced entries
*       1, ..., ULP  and random signs.
*  (5)  A diagonal matrix with "clustered" entries 1, ULP, ..., ULP
*       and random signs.
*
*  (6)  Same as (4), but multiplied by SQRT( overflow threshold )
*  (7)  Same as (4), but multiplied by SQRT( underflow threshold )
*
*  (8)  A matrix of the form  U' D U, where U is orthogonal and
*       D has evenly spaced entries 1, ..., ULP with random signs
*       on the diagonal.
*
*  (9)  A matrix of the form  U' D U, where U is orthogonal and
*       D has geometrically spaced entries 1, ..., ULP with random
*       signs on the diagonal.
*
*  (10) A matrix of the form  U' D U, where U is orthogonal and
*       D has "clustered" entries 1, ULP,..., ULP with random
*       signs on the diagonal.
*
*  (11) Same as (8), but multiplied by SQRT( overflow threshold )
*  (12) Same as (8), but multiplied by SQRT( underflow threshold )
*
*  (13) symmetric matrix with random entries chosen from (-1,1).
*  (14) Same as (13), but multiplied by SQRT( overflow threshold )
*  (15) Same as (13), but multiplied by SQRT( underflow threshold )
*  (16) Same as (8), but diagonal elements are all positive.
*  (17) Same as (9), but diagonal elements are all positive.
*  (18) Same as (10), but diagonal elements are all positive.
*  (19) Same as (16), but multiplied by SQRT( overflow threshold )
*  (20) Same as (16), but multiplied by SQRT( underflow threshold )
*  (21) A tridiagonal matrix that is a direct sum of smaller diagonally
*       dominant submatrices. Each unreduced submatrix has geometrically
*       spaced diagonal entries 1, ..., ULP.
*  (22) A matrix of the form  U' D U, where U is orthogonal and
*       D has ceil(lgN) "clusters" at 0,1,2,...,ceil(lgN)-1. The
*       size of the cluster at the value I is 2^I.
*
*  IBTYPE   (global input) INTEGER
*          Specifies the problem type to be solved:
*          = 1:  sub( A )*x = (lambda)*sub( B )*x
*          = 2:  sub( A )*sub( B )*x = (lambda)*x
*          = 3:  sub( B )*sub( A )*x = (lambda)*x
*
*
*  SUBTESTS (global input) CHARACTER*1
*           'Y' - Perform subset tests
*           'N' - Do not perform subset tests
*
*  THRESH   (global input) REAL
*          A test will count as "failed" if the "error", computed as
*          described below, exceeds THRESH.  Note that the error
*          is scaled to be O(1), so THRESH should be a reasonably
*          small multiple of 1, e.g., 10 or 100.  In particular,
*          it should not depend on the precision (single vs. double)
*          or the size of the matrix.  It must be at least zero.
*
*  ORDER    (global input) INTEGER
*           Number of reflectors used in test matrix creation.
*           If ORDER is large, it will
*           take more time to create the test matrices but they will
*           be closer to random.
*           ORDER .lt. N not implemented
*
*  ABSTOL   (global input) REAL
*           The absolute tolerance for the eigenvalues. An
*           eigenvalue is considered to be located if it has
*           been determined to lie in an interval whose width
*           is "abstol" or less. If "abstol" is less than or equal
*           to zero, then ulp*|T| will be used, where |T| is
*           the 1-norm of the matrix. If eigenvectors are
*           desired later by inverse iteration ("PSSTEIN"),
*           "abstol" MUST NOT be bigger than ulp*|T|.
*
*           For the purposes of this test, ABSTOL=0.0 is fine.
*           THis test does not test for high relative accuracy.
*
*  ISEED   (global input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  A       (local workspace) REAL array, dim (N*N)
*          global dimension (N, N), local dimension (LDA, NQ)
*            A is distributed in a block cyclic manner over both rows
*            and columns.  The actual location of a particular element
*            in A is controlled by the values of NPROW, NPCOL, and NB.
*          The test matrix, which is then modified by PSSYGVX
*
*  COPYA   (local workspace) REAL array, dim (N, N)
*          COPYA is used to hold an identical copy of the array A
*          identical in both form and content to A
*
*  B       (local workspace) REAL array, dim (N*N)
*          global dimension (N, N), local dimension (LDA, NQ)
*            A is distributed in a block cyclic manner over both rows
*            and columns.
*          The B test matrix, which is then modified by PSSYGVX
*
*  COPYB   (local workspace) REAL array, dim (N, N)
*          COPYB is used to hold an identical copy of the array B
*          identical in both form and content to B
*
*  Z       (local workspace) REAL array, dim (N*N)
*            Z is distributed in the same manner as A
*          Z is used as workspace by the test routines
*          PSGSEPCHK
*
*  W       (local workspace) REAL array, dimension (N)
*          On normal exit from PSSYGVX, the first M entries
*          contain the selected eigenvalues in ascending order.
*
*  IFAIL   (global workspace) INTEGER array, dimension (N)
*
*  WORK    (local workspace) REAL array, dimension (LWORK)
*
*  LWORK   (local input) INTEGER
*          The length of the array WORK.  LWORK >= SIZETST as
*          returned by PSLASIZEGSEP
*
*  IWORK   (local workspace) INTEGER array, dimension (LIWORK)
*
*  LIWORK  (local input) INTEGER
*          The length of the array IWORK.  LIWORK >= ISIZETST as
*          returned by PSLASIZEGSEP
*
*  NOUT   (local input) INTEGER
*         The unit number for output file.  Only used on node 0.
*         NOUT = 6, output to screen,
*         NOUT = 0, output to stderr.
*         NOUT = 13, output to file, divide thresh by 10.0
*         NOUT = 14, output to file, divide thresh by 20.0
*         (This hack allows us to test more stringently internally
*         so that when errors on found on other computers they will
*         be serious enough to warrant our attention.)
*
*  INFO (global output) INTEGER
*         -3       This process is not involved
*         0        Test succeeded (passed |AQ -QL| and |QT*Q - I| tests)
*         1        At least one test failed
*         2        Residual test were not performed, thresh <= 0.0
*         3        Test was skipped because of inadequate memory space
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ZERO, ONE, TEN, HALF
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TEN = 10.0E+0,
     $                   HALF = 0.5E+0 )
      REAL               PADVAL
      PARAMETER          ( PADVAL = 19.25E+0 )
      INTEGER            MAXTYP
      PARAMETER          ( MAXTYP = 22 )
*     ..
*
*     .. Local Scalars ..
      LOGICAL            WKNOWN
      CHARACTER          JOBZ, RANGE
      CHARACTER*14       PASSED
      INTEGER            CONTEXT, I, IAM, IINFO, IL, IMODE, IN, INDD,
     $                   INDWORK, ISIZESUBTST, ISIZESYEVX, ISIZETST,
     $                   ITYPE, IU, J, LLWORK, LSYEVXSIZE, MAXSIZE,
     $                   MYCOL, MYROW, NB, NGEN, NLOC, NNODES, NP,
     $                   NPCOL, NPROW, NQ, RES, SIZECHK, SIZEMQRLEFT,
     $                   SIZEMQRRIGHT, SIZEQRF, SIZEQTQ, SIZESUBTST,
     $                   SIZESYEVX, SIZETMS, SIZETST, VALSIZE, VECSIZE
      REAL               ANINV, ANORM, COND, MAXQTQNRM, MAXTSTNRM, OVFL,
     $                   QTQNRM, RTOVFL, RTUNFL, TEMP1, TSTNRM, ULP,
     $                   ULPINV, UNFL, VL, VU
*     ..
*     .. Local Arrays ..
      INTEGER            ISEEDIN( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ),
     $                   KTYPE( MAXTYP )
      DOUBLE PRECISION   CTIME( 10 ), WTIME( 10 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      REAL               PSLAMCH, SLARAN
      EXTERNAL           LSAME, NUMROC, PSLAMCH, SLARAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, BLACS_PINFO, IGAMX2D, IGEBR2D,
     $                   IGEBS2D, PSCHEKPAD, PSELSET, PSFILLPAD,
     $                   PSGSEPSUBTST, PSLASET, PSLASIZEGSEP,
     $                   PSLASIZESYEVX, PSLATMS, PSMATGEN, SLABAD,
     $                   SLASRT, SLATMS, SLCOMBINE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, MAX, MIN, MOD, REAL, SQRT
*     ..
*     .. Data statements ..
      DATA               KTYPE / 1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8,
     $                   8, 8, 9, 9, 9, 9, 9, 10, 11 /
      DATA               KMAGN / 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1,
     $                   2, 3, 1, 1, 1, 2, 3, 1, 1 /
      DATA               KMODE / 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0,
     $                   0, 0, 4, 3, 1, 4, 4, 3, 0 /
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      INFO = 0
      PASSED = 'PASSED      '
      CONTEXT = DESCA( CTXT_ )
      NB = DESCA( NB_ )
*
      CALL BLACS_PINFO( IAM, NNODES )
      CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
*
*     Make sure that we have enough memory
*
*
      CALL PSLASIZEGSEP( DESCA, IPREPAD, IPOSTPAD, SIZEMQRLEFT,
     $                   SIZEMQRRIGHT, SIZEQRF, SIZETMS, SIZEQTQ,
     $                   SIZECHK, SIZESYEVX, ISIZESYEVX, SIZESUBTST,
     $                   ISIZESUBTST, SIZETST, ISIZETST )
*
      IF( LWORK.LT.SIZETST ) THEN
         INFO = 3
      END IF
*
      CALL IGAMX2D( CONTEXT, 'a', ' ', 1, 1, INFO, 1, 1, 1, -1, -1, 0 )
*
      IF( INFO.EQ.0 ) THEN
*
         INDD = 1
         INDWORK = INDD + N
         LLWORK = LWORK - INDWORK + 1
*
         ULP = PSLAMCH( CONTEXT, 'P' )
         ULPINV = ONE / ULP
         UNFL = PSLAMCH( CONTEXT, 'Safe min' )
         OVFL = ONE / UNFL
         CALL SLABAD( UNFL, OVFL )
         RTUNFL = SQRT( UNFL )
         RTOVFL = SQRT( OVFL )
         ANINV = ONE / REAL( MAX( 1, N ) )
*
*     This ensures that everyone starts out with the same seed.
*
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            CALL IGEBS2D( CONTEXT, 'a', ' ', 4, 1, ISEED, 4 )
         ELSE
            CALL IGEBR2D( CONTEXT, 'a', ' ', 4, 1, ISEED, 4, 0, 0 )
         END IF
         ISEEDIN( 1 ) = ISEED( 1 )
         ISEEDIN( 2 ) = ISEED( 2 )
         ISEEDIN( 3 ) = ISEED( 3 )
         ISEEDIN( 4 ) = ISEED( 4 )
*
*     Compute the matrix A
*
*     Control parameters:
*
*     KMAGN  KMODE        KTYPE
*     =1  O(1)   clustered 1  zero
*     =2  large  clustered 2  identity
*     =3  small  exponential  (none)
*     =4         arithmetic   diagonal, (w/ eigenvalues)
*     =5         random log   symmetric, w/ eigenvalues
*     =6         random       (none)
*     =7                      random diagonal
*     =8                      random symmetric
*     =9                      positive definite
*     =10                     block diagonal with tridiagonal blocks
*     =11                     Geometrically sized clusters.
*
         ITYPE = KTYPE( MATTYPE )
         IMODE = KMODE( MATTYPE )
*
*     Compute norm
*
         GO TO ( 10, 20, 30 )KMAGN( MATTYPE )
*
   10    CONTINUE
         ANORM = ONE
         GO TO 40
*
   20    CONTINUE
         ANORM = ( RTOVFL*ULP )*ANINV
         GO TO 40
*
   30    CONTINUE
         ANORM = RTUNFL*N*ULPINV
         GO TO 40
*
   40    CONTINUE
         IF( MATTYPE.LE.15 ) THEN
            COND = ULPINV
         ELSE
            COND = ULPINV*ANINV / TEN
         END IF
*
*     Special Matrices
*
*     Zero
*
*
         IF( ITYPE.EQ.1 ) THEN
*
*     Zero Matrix
*
            DO 50 I = 1, N
               WORK( INDD+I-1 ) = ZERO
   50       CONTINUE
            CALL PSLASET( 'All', N, N, ZERO, ZERO, COPYA, 1, 1, DESCA )
            WKNOWN = .TRUE.
*
         ELSE IF( ITYPE.EQ.2 ) THEN
*
*     Identity Matrix
*
            DO 60 I = 1, N
               WORK( INDD+I-1 ) = ONE
   60       CONTINUE
            CALL PSLASET( 'All', N, N, ZERO, ONE, COPYA, 1, 1, DESCA )
            WKNOWN = .TRUE.
*
         ELSE IF( ITYPE.EQ.4 ) THEN
*
*     Diagonal Matrix, [Eigen]values Specified
*
            CALL PSFILLPAD( DESCA( CTXT_ ), SIZETMS, 1, WORK( INDWORK ),
     $                      SIZETMS, IPREPAD, IPOSTPAD, PADVAL+1.0E+0 )
*
            CALL PSLATMS( N, N, 'S', ISEED, 'S', WORK( INDD ), IMODE,
     $                    COND, ANORM, 0, 0, 'N', COPYA, 1, 1, DESCA,
     $                    ORDER, WORK( INDWORK+IPREPAD ), SIZETMS,
     $                    IINFO )
            WKNOWN = .TRUE.
*
            CALL PSCHEKPAD( DESCA( CTXT_ ), 'PSLATMS1-WORK', SIZETMS, 1,
     $                      WORK( INDWORK ), SIZETMS, IPREPAD, IPOSTPAD,
     $                      PADVAL+1.0E+0 )
*
         ELSE IF( ITYPE.EQ.5 ) THEN
*
*     symmetric, eigenvalues specified
*
            CALL PSFILLPAD( DESCA( CTXT_ ), SIZETMS, 1, WORK( INDWORK ),
     $                      SIZETMS, IPREPAD, IPOSTPAD, PADVAL+2.0E+0 )
*
            CALL PSLATMS( N, N, 'S', ISEED, 'S', WORK( INDD ), IMODE,
     $                    COND, ANORM, N, N, 'N', COPYA, 1, 1, DESCA,
     $                    ORDER, WORK( INDWORK+IPREPAD ), SIZETMS,
     $                    IINFO )
*
            CALL PSCHEKPAD( DESCA( CTXT_ ), 'PSLATMS2-WORK', SIZETMS, 1,
     $                      WORK( INDWORK ), SIZETMS, IPREPAD, IPOSTPAD,
     $                      PADVAL+2.0E+0 )
*
            WKNOWN = .TRUE.
*
         ELSE IF( ITYPE.EQ.8 ) THEN
*
*     symmetric, random eigenvalues
*
            NP = NUMROC( N, DESCA( MB_ ), MYROW, 0, NPROW )
            NQ = NUMROC( N, DESCA( NB_ ), MYCOL, 0, NPCOL )
            CALL PSMATGEN( DESCA( CTXT_ ), 'S', 'N', N, N, DESCA( MB_ ),
     $                     DESCA( NB_ ), COPYA, DESCA( LLD_ ),
     $                     DESCA( RSRC_ ), DESCA( CSRC_ ), ISEED( 1 ),
     $                     0, NP, 0, NQ, MYROW, MYCOL, NPROW, NPCOL )
            INFO = 0
            WKNOWN = .FALSE.
*
         ELSE IF( ITYPE.EQ.9 ) THEN
*
*     Positive definite, eigenvalues specified.
*
*
            CALL PSFILLPAD( DESCA( CTXT_ ), SIZETMS, 1, WORK( INDWORK ),
     $                      SIZETMS, IPREPAD, IPOSTPAD, PADVAL+3.0E+0 )
*
            CALL PSLATMS( N, N, 'S', ISEED, 'S', WORK( INDD ), IMODE,
     $                    COND, ANORM, N, N, 'N', COPYA, 1, 1, DESCA,
     $                    ORDER, WORK( INDWORK+IPREPAD ), SIZETMS,
     $                    IINFO )
*
            WKNOWN = .TRUE.
*
            CALL PSCHEKPAD( DESCA( CTXT_ ), 'PSLATMS3-WORK', SIZETMS, 1,
     $                      WORK( INDWORK ), SIZETMS, IPREPAD, IPOSTPAD,
     $                      PADVAL+3.0E+0 )
*
         ELSE IF( ITYPE.EQ.10 ) THEN
*
*     Block diagonal matrix with each block being a positive
*     definite tridiagonal submatrix.
*
            CALL PSLASET( 'All', N, N, ZERO, ZERO, COPYA, 1, 1, DESCA )
            NP = NUMROC( N, DESCA( MB_ ), 0, 0, NPROW )
            NQ = NUMROC( N, DESCA( NB_ ), 0, 0, NPCOL )
            NLOC = MIN( NP, NQ )
            NGEN = 0
   70       CONTINUE
*
            IF( NGEN.LT.N ) THEN
               IN = MIN( 1+INT( SLARAN( ISEED )*REAL( NLOC ) ), N-NGEN )
*
               CALL SLATMS( IN, IN, 'S', ISEED, 'P', WORK( INDD ),
     $                      IMODE, COND, ANORM, 1, 1, 'N', A, LDA,
     $                      WORK( INDWORK ), IINFO )
*
               DO 80 I = 2, IN
                  TEMP1 = ABS( A( I-1, I ) ) /
     $                    SQRT( ABS( A( I-1, I-1 )*A( I, I ) ) )
                  IF( TEMP1.GT.HALF ) THEN
                     A( I-1, I ) = HALF*SQRT( ABS( A( I-1, I-1 )*A( I,
     $                             I ) ) )
                     A( I, I-1 ) = A( I-1, I )
                  END IF
   80          CONTINUE
               CALL PSELSET( COPYA, NGEN+1, NGEN+1, DESCA, A( 1, 1 ) )
               DO 90 I = 2, IN
                  CALL PSELSET( COPYA, NGEN+I, NGEN+I, DESCA,
     $                          A( I, I ) )
                  CALL PSELSET( COPYA, NGEN+I-1, NGEN+I, DESCA,
     $                          A( I-1, I ) )
                  CALL PSELSET( COPYA, NGEN+I, NGEN+I-1, DESCA,
     $                          A( I, I-1 ) )
   90          CONTINUE
               NGEN = NGEN + IN
               GO TO 70
            END IF
            WKNOWN = .FALSE.
*
         ELSE IF( ITYPE.EQ.11 ) THEN
*
*     Geometrically sized clusters.  Eigenvalues:  0,1,1,2,2,2,2, ...
*
            NGEN = 0
            J = 1
            TEMP1 = ZERO
  100       CONTINUE
            IF( NGEN.LT.N ) THEN
               IN = MIN( J, N-NGEN )
               DO 110 I = 0, IN - 1
                  WORK( INDD+NGEN+I ) = TEMP1
  110          CONTINUE
               TEMP1 = TEMP1 + ONE
               J = 2*J
               NGEN = NGEN + IN
               GO TO 100
            END IF
*
*
            CALL PSFILLPAD( DESCA( CTXT_ ), SIZETMS, 1, WORK( INDWORK ),
     $                      SIZETMS, IPREPAD, IPOSTPAD, PADVAL+4.0E+0 )
*
            CALL PSLATMS( N, N, 'S', ISEED, 'S', WORK( INDD ), IMODE,
     $                    COND, ANORM, 0, 0, 'N', COPYA, 1, 1, DESCA,
     $                    ORDER, WORK( INDWORK+IPREPAD ), SIZETMS,
     $                    IINFO )
*
            CALL PSCHEKPAD( DESCA( CTXT_ ), 'PSLATMS4-WORK', SIZETMS, 1,
     $                      WORK( INDWORK ), SIZETMS, IPREPAD, IPOSTPAD,
     $                      PADVAL+4.0E+0 )
*
*
*     WKNOWN ... NOT SET, GUESS A DEFAULT
*
            WKNOWN = .TRUE.
         ELSE
            IINFO = 1
         END IF
*
         IF( WKNOWN )
     $      CALL SLASRT( 'I', N, WORK( INDD ), IINFO )
*
*    Create the B matrix
*
         CALL PSFILLPAD( DESCA( CTXT_ ), SIZETMS, 1, WORK( INDWORK ),
     $                   SIZETMS, IPREPAD, IPOSTPAD, PADVAL+3.3E+0 )
*
         ANORM = ONE
*
*           Update ISEED so that {SLAGSY creates a different Q
*
         ISEED( 4 ) = MOD( ISEED( 4 )+257, 4096 )
         ISEED( 3 ) = MOD( ISEED( 3 )+192, 4096 )
         ISEED( 2 ) = MOD( ISEED( 2 )+35, 4096 )
         ISEED( 1 ) = MOD( ISEED( 1 )+128, 4096 )
         CALL PSLATMS( N, N, 'S', ISEED, 'P', WORK( INDD ), 3, TEN,
     $                 ANORM, N, N, 'N', COPYB, 1, 1, DESCA, ORDER,
     $                 WORK( INDWORK+IPREPAD ), SIZETMS, IINFO )
*
         CALL PSCHEKPAD( DESCA( CTXT_ ), 'PSLATMS5-WORK', SIZETMS, 1,
     $                   WORK( INDWORK ), SIZETMS, IPREPAD, IPOSTPAD,
     $                   PADVAL+3.3E+0 )
*
*
*      These values aren't actually used, but they make ftncheck happy.
*
         IL = -1
         IU = -2
         VL = ONE
         VU = -ONE
*
         CALL PSLASIZESYEVX( WKNOWN, 'A', N, DESCA, VL, VU, IL, IU,
     $                       ISEED, WORK( INDD ), MAXSIZE, VECSIZE,
     $                       VALSIZE )
*
         LSYEVXSIZE = MIN( MAXSIZE, LWORK )
         WKNOWN = .FALSE.
*
         CALL PSGSEPSUBTST( WKNOWN, IBTYPE, 'v', 'a', UPLO, N, VL, VU,
     $                      IL, IU, THRESH, ABSTOL, A, COPYA, B, COPYB,
     $                      Z, 1, 1, DESCA, WORK( INDD ), WIN, IFAIL,
     $                      ICLUSTR, GAP, IPREPAD, IPOSTPAD,
     $                      WORK( INDWORK ), LLWORK, LSYEVXSIZE, IWORK,
     $                      ISIZESYEVX, RES, TSTNRM, QTQNRM, NOUT )
*
*
*
         MAXTSTNRM = TSTNRM
         MAXQTQNRM = QTQNRM
*
         IF( THRESH.LE.ZERO ) THEN
            PASSED = 'SKIPPED       '
            INFO = 2
         ELSE IF( RES.NE.0 ) THEN
            PASSED = 'FAILED        '
            INFO = 1
         END IF
      END IF
*
      IF( THRESH.GT.ZERO .AND. LSAME( SUBTESTS, 'Y' ) ) THEN
*
*        Subtest 1:  JOBZ = 'V', RANGE = 'A', minimum memory
*
         IF( INFO.EQ.0 ) THEN
*
            JOBZ = 'V'
            RANGE = 'A'
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
*
            LSYEVXSIZE = VECSIZE
*
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               PASSED = 'FAILED stest 1'
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               INFO = 1
            END IF
         END IF
*
*        Subtest 2:  JOBZ = 'V', RANGE = 'A', random memory
*
         IF( INFO.EQ.0 ) THEN
            JOBZ = 'V'
            RANGE = 'A'
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
*
            LSYEVXSIZE = VECSIZE + INT( SLARAN( ISEED )*
     $                   REAL( MAXSIZE-VECSIZE ) )
*
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               PASSED = 'FAILED stest 2'
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               INFO = 1
            END IF
         END IF
*
*        Subtest 3:  JOBZ = 'N', RANGE = 'A', minimum memory
*
         IF( INFO.EQ.0 ) THEN
*
            JOBZ = 'N'
            RANGE = 'A'
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
*
            LSYEVXSIZE = VALSIZE
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               PASSED = 'FAILED stest 3'
               INFO = 1
            END IF
         END IF
*
*        Subtest 4:  JOBZ = 'N', RANGE = 'I', minimum memory
*
         IF( INFO.EQ.0 ) THEN
*
            IL = -1
            IU = -1
            JOBZ = 'N'
            RANGE = 'I'
*
*           We use PSLASIZESYEVX to choose IL and IU for us.
*
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
*
            LSYEVXSIZE = VALSIZE
*
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               PASSED = 'FAILED stest 4'
               INFO = 1
            END IF
         END IF
*
*        Subtest 5:  JOBZ = 'V', RANGE = 'I', maximum memory
*
         IF( INFO.EQ.0 ) THEN
*
            IL = -1
            IU = -1
            JOBZ = 'V'
            RANGE = 'I'
*
*           We use PSLASIZESYEVX to choose IL and IU for us.
*
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
*
            LSYEVXSIZE = MAXSIZE
*
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               PASSED = 'FAILED stest 5'
               INFO = 1
            END IF
         END IF
*
*        Subtest 6:  JOBZ = 'V', RANGE = 'I', minimum memory
*
         IF( INFO.EQ.0 ) THEN
            IL = -1
            IU = -1
            JOBZ = 'V'
            RANGE = 'I'
*
*           We use PSLASIZESYEVX to choose IL and IU for us.
*
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
*
            LSYEVXSIZE = VECSIZE
*
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               PASSED = 'FAILED stest 6'
               INFO = 1
            END IF
         END IF
*
*        Subtest 7:  JOBZ = 'V', RANGE = 'I', random memory
*
         IF( INFO.EQ.0 ) THEN
            IL = -1
            IU = -1
            JOBZ = 'V'
            RANGE = 'I'
*
*           We use PSLASIZESYEVX to choose IL and IU for us.
*
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
            LSYEVXSIZE = VECSIZE + INT( SLARAN( ISEED )*
     $                   REAL( MAXSIZE-VECSIZE ) )
*
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               PASSED = 'FAILED stest 7'
               INFO = 1
            END IF
         END IF
*
*        Subtest 8:  JOBZ = 'N', RANGE = 'V', minimum memory
*
         IF( INFO.EQ.0 ) THEN
            VL = ONE
            VU = -ONE
            JOBZ = 'N'
            RANGE = 'V'
*
*           We use PSLASIZESYEVX to choose VL and VU for us.
*
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
*
            LSYEVXSIZE = VALSIZE
*
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               PASSED = 'FAILED stest 8'
               INFO = 1
            END IF
         END IF
*
*        Subtest 9:  JOBZ = 'V', RANGE = 'V', maximum memory
*
         IF( INFO.EQ.0 ) THEN
            VL = ONE
            VU = -ONE
            JOBZ = 'V'
            RANGE = 'V'
*
*           We use PSLASIZESYEVX to choose VL and VU for us.
*
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
*
            LSYEVXSIZE = MAXSIZE
*
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               PASSED = 'FAILED stest 9'
               INFO = 1
            END IF
         END IF
*
*        Subtest 10:  JOBZ = 'V', RANGE = 'V',
*                     minimum memory required for eigenvectors
*
         IF( INFO.EQ.0 ) THEN
            VL = ONE
            VU = -ONE
            JOBZ = 'V'
            RANGE = 'V'
*
*           We use PSLASIZESYEVX to choose VL and VU for us.
*
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
*
            LSYEVXSIZE = VECSIZE
*
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               PASSED = 'FAILED stest10'
               INFO = 1
            END IF
         END IF
*
*        Subtest 11:  JOBZ = 'V', RANGE = 'V',
*                     random memory (enough for all eigenvectors
*                     but not enough to guarantee orthogonality
*
         IF( INFO.EQ.0 ) THEN
            VL = ONE
            VU = -ONE
            JOBZ = 'V'
            RANGE = 'V'
*
*           We use PSLASIZESYEVX to choose VL and VU for us.
*
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
*
            LSYEVXSIZE = VECSIZE + INT( SLARAN( ISEED )*
     $                   REAL( MAXSIZE-VECSIZE ) )
*
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               PASSED = 'FAILED stest11'
               INFO = 1
            END IF
         END IF
*
*        Subtest 12:  JOBZ = 'V', RANGE = 'V',
*                    miniimum memory required for eigenvalues only
*
         IF( INFO.EQ.0 ) THEN
            VL = ONE
            VU = -ONE
            JOBZ = 'V'
            RANGE = 'V'
*
*           We use PSLASIZESYEVX to choose VL and VU for us.
*
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
*
            LSYEVXSIZE = VALSIZE
*
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               PASSED = 'FAILED stest12'
               INFO = 1
            END IF
         END IF
*
*        Subtest 13:  JOBZ = 'V', RANGE = 'V',
*                     random memory (more than minimum required
*                     for eigenvalues, less than required for vectors)
*
         IF( INFO.EQ.0 ) THEN
            VL = ONE
            VU = -ONE
            JOBZ = 'V'
            RANGE = 'V'
*
*           We use PSLASIZESYEVX to choose VL and VU for us.
*
            CALL PSLASIZESYEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
*
            LSYEVXSIZE = VALSIZE + INT( SLARAN( ISEED )*
     $                   REAL( VECSIZE-VALSIZE ) )
*
            CALL PSGSEPSUBTST( .TRUE., IBTYPE, JOBZ, RANGE, UPLO, N, VL,
     $                         VU, IL, IU, THRESH, ABSTOL, A, COPYA, B,
     $                         COPYB, Z, 1, 1, DESCA, WIN( 1+IPREPAD ),
     $                         WNEW, IFAIL, ICLUSTR, GAP, IPREPAD,
     $                         IPOSTPAD, WORK( INDWORK ), LLWORK,
     $                         LSYEVXSIZE, IWORK, ISIZESYEVX, RES,
     $                         TSTNRM, QTQNRM, NOUT )
*
            IF( RES.NE.0 ) THEN
               MAXTSTNRM = MAX( TSTNRM, MAXTSTNRM )
               MAXQTQNRM = MAX( QTQNRM, MAXQTQNRM )
               PASSED = 'FAILED stest13'
               INFO = 1
            END IF
         END IF
      END IF
*
*
*
      CALL IGAMX2D( CONTEXT, 'All', ' ', 1, 1, INFO, 1, -1, -1, -1, -1,
     $              -1 )
*
      IF( INFO.EQ.1 ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9994 )'C  '
            WRITE( NOUT, FMT = 9993 )ISEEDIN( 1 )
            WRITE( NOUT, FMT = 9992 )ISEEDIN( 2 )
            WRITE( NOUT, FMT = 9991 )ISEEDIN( 3 )
            WRITE( NOUT, FMT = 9990 )ISEEDIN( 4 )
            IF( LSAME( UPLO, 'L' ) ) THEN
               WRITE( NOUT, FMT = 9994 )'      UPLO= ''L'' '
            ELSE
               WRITE( NOUT, FMT = 9994 )'      UPLO= ''U'' '
            END IF
            IF( LSAME( SUBTESTS, 'Y' ) ) THEN
               WRITE( NOUT, FMT = 9994 )'      SUBTESTS= ''Y'' '
            ELSE
               WRITE( NOUT, FMT = 9994 )'      SUBTESTS= ''N'' '
            END IF
            WRITE( NOUT, FMT = 9989 )N
            WRITE( NOUT, FMT = 9988 )NPROW
            WRITE( NOUT, FMT = 9987 )NPCOL
            WRITE( NOUT, FMT = 9986 )NB
            WRITE( NOUT, FMT = 9985 )MATTYPE
            WRITE( NOUT, FMT = 9984 )IBTYPE
            WRITE( NOUT, FMT = 9982 )ABSTOL
            WRITE( NOUT, FMT = 9981 )THRESH
            WRITE( NOUT, FMT = 9994 )'C  '
         END IF
      END IF
*
      CALL SLCOMBINE( CONTEXT, 'All', '>', 'W', 6, 1, WTIME )
      CALL SLCOMBINE( CONTEXT, 'All', '>', 'C', 6, 1, CTIME )
      IF( IAM.EQ.0 ) THEN
         IF( INFO.EQ.0 .OR. INFO.EQ.1 ) THEN
            IF( WTIME( 1 ).GE.0.0 ) THEN
               WRITE( NOUT, FMT = 9999 )N, NB, NPROW, NPCOL, MATTYPE,
     $            IBTYPE, SUBTESTS, WTIME( 1 ), CTIME( 1 ), MAXTSTNRM,
     $            PASSED
            ELSE
               WRITE( NOUT, FMT = 9998 )N, NB, NPROW, NPCOL, MATTYPE,
     $            IBTYPE, SUBTESTS, CTIME( 1 ), MAXTSTNRM, PASSED
            END IF
         ELSE IF( INFO.EQ.2 ) THEN
            IF( WTIME( 1 ).GE.0.0 ) THEN
               WRITE( NOUT, FMT = 9997 )N, NB, NPROW, NPCOL, MATTYPE,
     $            IBTYPE, SUBTESTS, WTIME( 1 ), CTIME( 1 )
            ELSE
               WRITE( NOUT, FMT = 9996 )N, NB, NPROW, NPCOL, MATTYPE,
     $            IBTYPE, SUBTESTS, CTIME( 1 )
            END IF
         ELSE IF( INFO.EQ.3 ) THEN
            WRITE( NOUT, FMT = 9995 )N, NB, NPROW, NPCOL, MATTYPE,
     $         IBTYPE, SUBTESTS
         END IF
      END IF
*
  120 CONTINUE
*
      RETURN
 9999 FORMAT( 1X, I5, 1X, I3, 1X, I3, 1X, I3, 1X, I3, 3X, I3, 4X, A1,
     $      1X, F8.2, 1X, F8.2, 1X, G9.2, 1X, A14 )
 9998 FORMAT( 1X, I5, 1X, I3, 1X, I3, 1X, I3, 1X, I3, 3X, I3, 4X, A1,
     $      1X, 8X, 1X, F8.2, 1X, G9.2, A14 )
 9997 FORMAT( 1X, I5, 1X, I3, 1X, I3, 1X, I3, 1X, I3, 3X, I3, 4X, A1,
     $      1X, F8.2, 1X, F8.2, 11X, 'Bypassed' )
 9996 FORMAT( 1X, I5, 1X, I3, 1X, I3, 1X, I3, 1X, I3, 3X, I3, 4X, A1,
     $      1X, 8X, 1X, F8.2, 11X, 'Bypassed' )
 9995 FORMAT( 1X, I5, 1X, I3, 1X, I3, 1X, I3, 1X, I3, 3X, I3, 4X, A1,
     $      22X, 'Bad MEMORY parameters' )
 9994 FORMAT( A )
 9993 FORMAT( '      ISEED( 1 ) =', I8 )
 9992 FORMAT( '      ISEED( 2 ) =', I8 )
 9991 FORMAT( '      ISEED( 3 ) =', I8 )
 9990 FORMAT( '      ISEED( 4 ) =', I8 )
 9989 FORMAT( '      N=', I8 )
 9988 FORMAT( '      NPROW=', I8 )
 9987 FORMAT( '      NPCOL=', I8 )
 9986 FORMAT( '      NB=', I8 )
 9985 FORMAT( '      MATTYPE=', I8 )
 9984 FORMAT( '      IBTYPE=', I8 )
 9983 FORMAT( '      SUBTESTS=', A1 )
 9982 FORMAT( '      ABSTOL=', D16.6 )
 9981 FORMAT( '      THRESH=', D16.6 )
 9980 FORMAT( ' Increase TOTMEM in PSGSEPDRIVER' )
*
*     End of PSGSEPTST
*
      END
