*
*
      SUBROUTINE PZSEPSUBTST( WKNOWN, JOBZ, RANGE, UPLO, N, VL, VU, IL,
     $                        IU, THRESH, ABSTOL, A, COPYA, Z, IA, JA,
     $                        DESCA, WIN, WNEW, IFAIL, ICLUSTR, GAP,
     $                        IPREPAD, IPOSTPAD, WORK, LWORK, RWORK,
     $                        LRWORK, LWORK1, IWORK, LIWORK, RESULT,
     $                        TSTNRM, QTQNRM, NOUT )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*     .. Scalar Arguments ..
      LOGICAL            WKNOWN
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IA, IL, IPOSTPAD, IPREPAD, IU, JA, LIWORK,
     $                   LRWORK, LWORK, LWORK1, N, NOUT, RESULT
      DOUBLE PRECISION   ABSTOL, QTQNRM, THRESH, TSTNRM, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), ICLUSTR( * ), IFAIL( * ),
     $                   IWORK( * )
      DOUBLE PRECISION   GAP( * ), RWORK( * ), WIN( * ), WNEW( * )
      COMPLEX*16         A( * ), COPYA( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PZSEPSUBTST calls PZHEEVX and then tests the output of
*  PZHEEVX
*  If JOBZ = 'V' then the following two tests are performed:
*     |AQ -QL| / (abstol + eps * norm(A) ) < N*THRESH
*     |QT * Q - I| / eps * norm(A) < N*THRESH
*  If WKNOWN then
*     we check to make sure that the eigenvalues match expectations
*     i.e. |WIN - WNEW(1+IPREPAD)| / (eps * |WIN|) < THRESH
*     where WIN is the array of eigenvalues as computed by
*     PZHEEVX when eigenvectors are requested
*
*  Arguments
*  =========
*
*     NP = the number of rows local to a given process.
*     NQ = the number of columns local to a given process.
*
*  WKNOWN  (global input) INTEGER
*          .FALSE.:  WIN does not contain the eigenvalues
*          .TRUE.:   WIN does contain the eigenvalues
*
*  JOBZ    (global input) CHARACTER*1
*          Specifies whether or not to compute the eigenvectors:
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvalues and eigenvectors.
*          Must be 'V' on first call to PZSEPSUBTST
*
*  RANGE   (global input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the interval [VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*          Must be 'A' on first call to PZSEPSUBTST
*
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          Size of the matrix to be tested.  (global size)
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
*  THRESH  (global input) DOUBLE PRECISION
*          A test will count as "failed" if the "error", computed as
*          described below, exceeds THRESH.  Note that the error
*          is scaled to be O(1), so THRESH should be a reasonably
*          small multiple of 1, e.g., 10 or 100.  In particular,
*          it should not depend on the precision (single vs. double)
*          or the size of the matrix.  It must be at least zero.
*
*  ABSTOL  (global input) DOUBLE PRECISION
*          The absolute tolerance for the eigenvalues. An
*          eigenvalue is considered to be located if it has
*          been determined to lie in an interval whose width
*          is "abstol" or less. If "abstol" is less than or equal
*          to zero, then ulp*|T| will be used, where |T| is
*          the 1-norm of the matrix. If eigenvectors are
*          desired later by inverse iteration ("PZSTEIN"),
*          "abstol" MUST NOT be bigger than ulp*|T|.
*
*  A       (local workspace) COMPLEX*16 array
*          global dimension (N, N), local dimension (DESCA(DLEN_), NQ)
*            A is distributed in a block cyclic manner over both rows
*            and columns.
*            See PZHEEVX for a description of block cyclic layout.
*          The test matrix, which is then modified by PZHEEVX
*          A has already been padded front and back, use A(1+IPREPAD)
*
*  COPYA   (local input) COMPLEX*16 array, dimension(N*N)
*          COPYA holds a copy of the original matrix A
*          identical in both form and content to A
*
*  Z       (local workspace) COMPLEX*16 array, dim (N*N)
*          Z is distributed in the same manner as A
*          Z contains the eigenvector matrix
*          Z is used as workspace by the test routines
*          PZSEPCHK and PZSEPQTQ.
*          Z has already been padded front and back, use Z(1+IPREPAD)
*
*  IA      (global input) INTEGER
*          On entry, IA specifies the global row index of the submatrix
*          of the global matrix A, COPYA and Z to operate on.
*
*  JA      (global input) INTEGER
*          On entry, IA specifies the global column index of the submat
*          of the global matrix A, COPYA and Z to operate on.
*
*  DESCA   (global/local input) INTEGER array of dimension 8
*          The array descriptor for the matrix A, COPYA and Z.
*
*  WIN     (global input) DOUBLE PRECISION array, dimension (N)
*          If .not. WKNOWN, WIN is ignored on input
*          Otherwise, WIN() is taken as the standard by which the
*          eigenvalues are to be compared against.
*
*  WNEW    (global workspace)  DOUBLE PRECISION array, dimension (N)
*          The eigenvalues as copmuted by this call to PZHEEVX
*          If JOBZ <> 'V' or RANGE <> 'A' these eigenvalues are
*          compared against those in WIN().
*          WNEW has already been padded front and back,
*          use WNEW(1+IPREPAD)
*
*  IFAIL   (global output) INTEGER array, dimension (N)
*          If JOBZ = 'V', then on normal exit, the first M elements of
*          IFAIL are zero.  If INFO > 0 on exit, then IFAIL contains the
*          indices of the eigenvectors that failed to converge.
*          If JOBZ = 'N', then IFAIL is not referenced.
*          IFAIL has already been padded front and back,
*          use IFAIL(1+IPREPAD)
*
*  ICLUSTR (global workspace) integer array, dimension (2*NPROW*NPCOL)
*
*  GAP     (global workspace) DOUBLE PRECISION array,
*          dimension (NPROW*NPCOL)
*
*  WORK    (local workspace) COMPLEX*16 array, dimension (LWORK)
*          WORK has already been padded front and back,
*          use WORK(1+IPREPAD)
*
*  LWORK   (local input) INTEGER
*          The actual length of the array WORK after padding.
*
*  RWORK   (local workspace) DOUBLE PRECISION array, dimension (LRWORK)
*          RWORK has already been padded front and back,
*          use RWORK(1+IPREPAD)
*
*  LRWORK   (local input) INTEGER
*          The actual length of the array RWORK after padding.
*
*  LWORK1  (local input) INTEGER
*          The amount of real workspace to pass to PZHEEVX
*
*  IWORK   (local workspace) INTEGER array, dimension (LIWORK)
*          IWORK has already been padded front and back,
*          use IWORK(1+IPREPAD)
*
*  LIWORK  (local input) INTEGER
*          The length of the array IWORK after padding.
*
*  RESULT  (global output) INTEGER
*          The result of this call to PZHEEVX
*          RESULT = -3   =>  This process did not participate
*          RESULT = 0    =>  All tests passed
*          RESULT = 1    =>  ONe or more tests failed
*
*  TSTNRM  (global output) DOUBLE PRECISION
*          |AQ- QL| / (ABSTOL+EPS*|A|)*N
*
*  QTQNRM  (global output) DOUBLE PRECISION
*          |QTQ -I| / N*EPS
*
*     .. Parameters ..
*
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   PADVAL, FIVE, NEGONE
      PARAMETER          ( PADVAL = 13.5285D+0, FIVE = 5.0D+0,
     $                   NEGONE = -1.0D+0 )
      COMPLEX*16         ZPADVAL
      PARAMETER          ( ZPADVAL = ( 13.989D+0, 1.93D+0 ) )
      INTEGER            IPADVAL
      PARAMETER          ( IPADVAL = 927 )
*     ..
*     .. Local Scalars ..
      LOGICAL            MISSLARGEST, MISSSMALLEST
      INTEGER            I, IAM, INDIWRK, INFO, ISIZEHEEVX, ISIZESUBTST,
     $                   ISIZETST, J, M, MAXEIGS, MAXIL, MAXIU, MAXSIZE,
     $                   MINIL, MQ, MYCOL, MYIL, MYROW, NCLUSTERS, NP,
     $                   NPCOL, NPROW, NQ, NZ, OLDIL, OLDIU, OLDNZ, RES,
     $                   RSIZECHK, RSIZEHEEVX, RSIZEQTQ, RSIZESUBTST,
     $                   RSIZETST, SIZEHEEVX, SIZEMQRLEFT, SIZEMQRRIGHT,
     $                   SIZEQRF, SIZESUBTST, SIZETMS, SIZETST, VALSIZE,
     $                   VECSIZE,  SIZEHEEVD, RSIZEHEEVD, ISIZEHEEVD
      DOUBLE PRECISION   EPS, EPSNORMA, ERROR, MAXERROR, MAXVU,
     $                   MINERROR, MINVL, NORMWIN, OLDVL, OLDVU, ORFAC,
     $                   SAFMIN
*     ..
*     .. Local Arrays ..
      INTEGER            DESCZ( DLEN_ ), DSEED( 4 ), ITMP( 2 )
*     ..
*     .. External Functions ..
*
      LOGICAL            LSAME
      INTEGER            NUMROC
      DOUBLE PRECISION   PDLAMCH, PZLANHE
      EXTERNAL           LSAME, NUMROC, PDLAMCH, PZLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCINIT, DGAMN2D, DGAMX2D,
     $                   IGAMN2D, IGAMX2D, PDCHEKPAD, PDFILLPAD,
     $                   PICHEKPAD, PIFILLPAD, PZCHEKPAD, PZELSET,
     $                   PZFILLPAD, PZHEEVX, PZLASIZEHEEVX, PZLASIZESEP,
     $                   PZSEPCHK, PZSEPQTQ, SLBOOT, SLTIMER, ZLACPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
      CALL PZLASIZESEP( DESCA, IPREPAD, IPOSTPAD, SIZEMQRLEFT,
     $                  SIZEMQRRIGHT, SIZEQRF, SIZETMS, RSIZEQTQ,
     $                  RSIZECHK, SIZEHEEVX, RSIZEHEEVX, ISIZEHEEVX,
     $                  SIZEHEEVD, RSIZEHEEVD, ISIZEHEEVD,
     $                  SIZESUBTST, RSIZESUBTST, ISIZESUBTST, SIZETST,
     $                  RSIZETST, ISIZETST )
*
      TSTNRM = NEGONE
      QTQNRM = NEGONE
      EPS = PDLAMCH( DESCA( CTXT_ ), 'Eps' )
      SAFMIN = PDLAMCH( DESCA( CTXT_ ), 'Safe min' )
*
      NORMWIN = SAFMIN / EPS
      IF( N.GE.1 )
     $   NORMWIN = MAX( ABS( WIN( 1 ) ), ABS( WIN( N ) ), NORMWIN )
*
*     Make sure that we aren't using information from previous calls
*
      NZ = -13
      OLDNZ = NZ
      OLDIL = IL
      OLDIU = IU
      OLDVL = VL
      OLDVU = VU
*
      DO 10 I = 1, LWORK1, 1
         RWORK( I+IPREPAD ) = 14.3D+0
   10 CONTINUE
      DO 20 I = 1, LIWORK, 1
         IWORK( I+IPREPAD ) = 14
   20 CONTINUE
      DO 30 I = 1, LWORK, 1
         WORK( I+IPREPAD ) = ( 15.63D+0, 1.1D+0 )
   30 CONTINUE
*
      DO 40 I = 1, N
         WNEW( I+IPREPAD ) = 3.14159D+0
   40 CONTINUE
*
      ICLUSTR( 1+IPREPAD ) = 139
*
      IF( LSAME( JOBZ, 'N' ) ) THEN
         MAXEIGS = 0
      ELSE
         IF( LSAME( RANGE, 'A' ) ) THEN
            MAXEIGS = N
         ELSE IF( LSAME( RANGE, 'I' ) ) THEN
            MAXEIGS = IU - IL + 1
         ELSE
            MINVL = VL - NORMWIN*FIVE*EPS - ABSTOL
            MAXVU = VU + NORMWIN*FIVE*EPS + ABSTOL
            MINIL = 1
            MAXIU = 0
            DO 50 I = 1, N
               IF( WIN( I ).LT.MINVL )
     $            MINIL = MINIL + 1
               IF( WIN( I ).LE.MAXVU )
     $            MAXIU = MAXIU + 1
   50       CONTINUE
*
            MAXEIGS = MAXIU - MINIL + 1
         END IF
      END IF
*
*
      CALL DESCINIT( DESCZ, DESCA( M_ ), DESCA( N_ ), DESCA( MB_ ),
     $               DESCA( NB_ ), DESCA( RSRC_ ), DESCA( CSRC_ ),
     $               DESCA( CTXT_ ), DESCA( LLD_ ), INFO )
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      INDIWRK = 1 + IPREPAD + NPROW*NPCOL + 1
*
      IAM = 1
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $   IAM = 0
*
*     If this process is not involved in this test, bail out now
*
      RESULT = -3
      IF( MYROW.GE.NPROW .OR. MYROW.LT.0 )
     $   GO TO 160
      RESULT = 0
*
*
*     DSEED is not used in this call to PZLASIZEHEEVX, the
*     following line just makes ftnchek happy.
*
      DSEED( 1 ) = 1
*
      CALL PZLASIZEHEEVX( WKNOWN, RANGE, N, DESCA, VL, VU, IL, IU,
     $                    DSEED, WIN, MAXSIZE, VECSIZE, VALSIZE )
*
      NP = NUMROC( N, DESCA( MB_ ), MYROW, 0, NPROW )
      NQ = NUMROC( N, DESCA( NB_ ), MYCOL, 0, NPCOL )
      MQ = NUMROC( MAXEIGS, DESCA( NB_ ), MYCOL, 0, NPCOL )
*
      CALL ZLACPY( 'A', NP, NQ, COPYA, DESCA( LLD_ ), A( 1+IPREPAD ),
     $             DESCA( LLD_ ) )
*
      CALL PZFILLPAD( DESCA( CTXT_ ), NP, NQ, A, DESCA( LLD_ ), IPREPAD,
     $                IPOSTPAD, ZPADVAL )
*
      CALL PZFILLPAD( DESCZ( CTXT_ ), NP, MQ, Z, DESCZ( LLD_ ), IPREPAD,
     $                IPOSTPAD, ZPADVAL+1.0D+0 )
*
      CALL PDFILLPAD( DESCA( CTXT_ ), N, 1, WNEW, N, IPREPAD, IPOSTPAD,
     $                PADVAL+2.0D+0 )
*
      CALL PDFILLPAD( DESCA( CTXT_ ), NPROW*NPCOL, 1, GAP, NPROW*NPCOL,
     $                IPREPAD, IPOSTPAD, PADVAL+3.0D+0 )
*
      CALL PDFILLPAD( DESCA( CTXT_ ), LWORK1, 1, RWORK, LWORK1, IPREPAD,
     $                IPOSTPAD, PADVAL+4.0D+0 )
*
      CALL PIFILLPAD( DESCA( CTXT_ ), LIWORK, 1, IWORK, LIWORK, IPREPAD,
     $                IPOSTPAD, IPADVAL )
*
      CALL PIFILLPAD( DESCA( CTXT_ ), N, 1, IFAIL, N, IPREPAD, IPOSTPAD,
     $                IPADVAL )
*
      CALL PIFILLPAD( DESCA( CTXT_ ), 2*NPROW*NPCOL, 1, ICLUSTR,
     $                2*NPROW*NPCOL, IPREPAD, IPOSTPAD, IPADVAL )
*
      CALL PZFILLPAD( DESCA( CTXT_ ), LWORK, 1, WORK, LWORK, IPREPAD,
     $                IPOSTPAD, ZPADVAL+4.1D+0 )
*
*     Make sure that PZHEEVX does not cheat (i.e. use answers
*     already computed.)
*
      DO 70 I = 1, N, 1
         DO 60 J = 1, MAXEIGS, 1
            CALL PZELSET( Z( 1+IPREPAD ), I, J, DESCA,
     $                    ( 13.0D+0, 1.34D+0 ) )
   60    CONTINUE
   70 CONTINUE
*
      ORFAC = -1.0D+0
*
      CALL SLBOOT
      CALL SLTIMER( 1 )
      CALL SLTIMER( 6 )
      CALL PZHEEVX( JOBZ, RANGE, UPLO, N, A( 1+IPREPAD ), IA, JA, DESCA,
     $              VL, VU, IL, IU, ABSTOL, M, NZ, WNEW( 1+IPREPAD ),
     $              ORFAC, Z( 1+IPREPAD ), IA, JA, DESCA,
     $              WORK( 1+IPREPAD ), SIZEHEEVX, RWORK( 1+IPREPAD ),
     $              LWORK1, IWORK( 1+IPREPAD ), LIWORK,
     $              IFAIL( 1+IPREPAD ), ICLUSTR( 1+IPREPAD ),
     $              GAP( 1+IPREPAD ), INFO )
      CALL SLTIMER( 6 )
      CALL SLTIMER( 1 )
*
      IF( THRESH.LE.0 ) THEN
         RESULT = 0
      ELSE
         CALL PZCHEKPAD( DESCA( CTXT_ ), 'PZHEEVX-A', NP, NQ, A,
     $                   DESCA( LLD_ ), IPREPAD, IPOSTPAD, ZPADVAL )
*
         CALL PZCHEKPAD( DESCZ( CTXT_ ), 'PZHEEVX-Z', NP, MQ, Z,
     $                   DESCZ( LLD_ ), IPREPAD, IPOSTPAD,
     $                   ZPADVAL+1.0D+0 )
*
         CALL PDCHEKPAD( DESCA( CTXT_ ), 'PZHEEVX-WNEW', N, 1, WNEW, N,
     $                   IPREPAD, IPOSTPAD, PADVAL+2.0D+0 )
*
         CALL PDCHEKPAD( DESCA( CTXT_ ), 'PZHEEVX-GAP', NPROW*NPCOL, 1,
     $                   GAP, NPROW*NPCOL, IPREPAD, IPOSTPAD,
     $                   PADVAL+3.0D+0 )
*
         CALL PDCHEKPAD( DESCA( CTXT_ ), 'PZHEEVX-rWORK', LWORK1, 1,
     $                   RWORK, LWORK1, IPREPAD, IPOSTPAD,
     $                   PADVAL+4.0D+0 )
*
         CALL PZCHEKPAD( DESCA( CTXT_ ), 'PZHEEVX-WORK', LWORK, 1, WORK,
     $                   LWORK, IPREPAD, IPOSTPAD, ZPADVAL+4.1D+0 )
*
         CALL PICHEKPAD( DESCA( CTXT_ ), 'PZHEEVX-IWORK', LIWORK, 1,
     $                   IWORK, LIWORK, IPREPAD, IPOSTPAD, IPADVAL )
*
         CALL PICHEKPAD( DESCA( CTXT_ ), 'PZHEEVX-IFAIL', N, 1, IFAIL,
     $                   N, IPREPAD, IPOSTPAD, IPADVAL )
*
         CALL PICHEKPAD( DESCA( CTXT_ ), 'PZHEEVX-ICLUSTR',
     $                   2*NPROW*NPCOL, 1, ICLUSTR, 2*NPROW*NPCOL,
     $                   IPREPAD, IPOSTPAD, IPADVAL )
*
*
*     Since we now know the spectrum, we can potentially reduce MAXSIZE.
*
         IF( LSAME( RANGE, 'A' ) ) THEN
            CALL PZLASIZEHEEVX( .TRUE., RANGE, N, DESCA, VL, VU, IL, IU,
     $                          DSEED, WNEW( 1+IPREPAD ), MAXSIZE,
     $                          VECSIZE, VALSIZE )
         END IF
*
*
*     Check INFO
*
*
*     Make sure that all processes return the same value of INFO
*
         ITMP( 1 ) = INFO
         ITMP( 2 ) = INFO
*
         CALL IGAMN2D( DESCA( CTXT_ ), 'a', ' ', 1, 1, ITMP, 1, 1, 1,
     $                 -1, -1, 0 )
         CALL IGAMX2D( DESCA( CTXT_ ), 'a', ' ', 1, 1, ITMP( 2 ), 1, 1,
     $                 1, -1, -1, 0 )
*
*
         IF( ITMP( 1 ).NE.ITMP( 2 ) ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = * )
     $         'Different processes return different INFO'
            RESULT = 1
         ELSE IF( MOD( INFO, 2 ).EQ.1 .OR. INFO.GT.7 .OR. INFO.LT.0 )
     $             THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9999 )INFO
            RESULT = 1
         ELSE IF( MOD( INFO / 2, 2 ).EQ.1 .AND. LWORK1.GE.MAXSIZE ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9996 )INFO
            RESULT = 1
         ELSE IF( MOD( INFO / 4, 2 ).EQ.1 .AND. LWORK1.GE.VECSIZE ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9996 )INFO
            RESULT = 1
         END IF
*
*
         IF( LSAME( JOBZ, 'V' ) .AND. ( ICLUSTR( 1+IPREPAD ).NE.
     $       0 ) .AND. ( MOD( INFO / 2, 2 ).NE.1 ) ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9995 )
            RESULT = 1
         END IF
*
*     Check M
*
         IF( ( M.LT.0 ) .OR. ( M.GT.N ) ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9994 )
            RESULT = 1
         ELSE IF( LSAME( RANGE, 'A' ) .AND. ( M.NE.N ) ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9993 )
            RESULT = 1
         ELSE IF( LSAME( RANGE, 'I' ) .AND. ( M.NE.IU-IL+1 ) ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9992 )
            RESULT = 1
         ELSE IF( LSAME( JOBZ, 'V' ) .AND.
     $            ( .NOT.( LSAME( RANGE, 'V' ) ) ) .AND. ( M.NE.NZ ) )
     $             THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9991 )
            RESULT = 1
         END IF
*
*        Check NZ
*
         IF( LSAME( JOBZ, 'V' ) ) THEN
            IF( LSAME( RANGE, 'V' ) ) THEN
               IF( NZ.GT.M ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9990 )
                  RESULT = 1
               END IF
               IF( NZ.LT.M .AND. MOD( INFO / 4, 2 ).NE.1 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9989 )
                  RESULT = 1
               END IF
            ELSE
               IF( NZ.NE.M ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9988 )
                  RESULT = 1
               END IF
            END IF
         END IF
         IF( RESULT.EQ.0 ) THEN
*
*     Make sure that all processes return the same # of eigenvalues
*
            ITMP( 1 ) = M
            ITMP( 2 ) = M
*
            CALL IGAMN2D( DESCA( CTXT_ ), 'a', ' ', 1, 1, ITMP, 1, 1, 1,
     $                    -1, -1, 0 )
            CALL IGAMX2D( DESCA( CTXT_ ), 'a', ' ', 1, 1, ITMP( 2 ), 1,
     $                    1, 1, -1, -1, 0 )
*
            IF( ITMP( 1 ).NE.ITMP( 2 ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9987 )
               RESULT = 1
            ELSE
*
*     Make sure that different processes return the same eigenvalues
*
               DO 80 I = 1, M
                  RWORK( I ) = WNEW( I+IPREPAD )
                  RWORK( I+M ) = WNEW( I+IPREPAD )
   80          CONTINUE
*
               CALL DGAMN2D( DESCA( CTXT_ ), 'a', ' ', M, 1, RWORK, M,
     $                       1, 1, -1, -1, 0 )
               CALL DGAMX2D( DESCA( CTXT_ ), 'a', ' ', M, 1,
     $                       RWORK( 1+M ), M, 1, 1, -1, -1, 0 )
*
               DO 90 I = 1, M
*
                  IF( RESULT.EQ.0 .AND. ( ABS( RWORK( I )-RWORK( M+
     $                I ) ).GT.FIVE*EPS*ABS( RWORK( I ) ) ) ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9986 )
                     RESULT = 1
                  END IF
   90          CONTINUE
            END IF
         END IF
*
*     Make sure that all processes return the same # of clusters
*
         IF( LSAME( JOBZ, 'V' ) ) THEN
            NCLUSTERS = 0
            DO 100 I = 0, NPROW*NPCOL - 1
               IF( ICLUSTR( 1+IPREPAD+2*I ).EQ.0 )
     $            GO TO 110
               NCLUSTERS = NCLUSTERS + 1
  100       CONTINUE
  110       CONTINUE
            ITMP( 1 ) = NCLUSTERS
            ITMP( 2 ) = NCLUSTERS
*
            CALL IGAMN2D( DESCA( CTXT_ ), 'a', ' ', 1, 1, ITMP, 1, 1, 1,
     $                    -1, -1, 0 )
            CALL IGAMX2D( DESCA( CTXT_ ), 'a', ' ', 1, 1, ITMP( 2 ), 1,
     $                    1, 1, -1, -1, 0 )
*
            IF( ITMP( 1 ).NE.ITMP( 2 ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9985 )
               RESULT = 1
            ELSE
*
*     Make sure that different processes return the same clusters
*
               DO 120 I = 1, NCLUSTERS
                  IWORK( INDIWRK+I ) = ICLUSTR( I+IPREPAD )
                  IWORK( INDIWRK+I+NCLUSTERS ) = ICLUSTR( I+IPREPAD )
  120          CONTINUE
               CALL IGAMN2D( DESCA( CTXT_ ), 'a', ' ', NCLUSTERS*2+1, 1,
     $                       IWORK( INDIWRK+1 ), NCLUSTERS*2+1, 1, 1,
     $                       -1, -1, 0 )
               CALL IGAMX2D( DESCA( CTXT_ ), 'a', ' ', NCLUSTERS*2+1, 1,
     $                       IWORK( INDIWRK+1+NCLUSTERS ),
     $                       NCLUSTERS*2+1, 1, 1, -1, -1, 0 )
*
*
               DO 130 I = 1, NCLUSTERS
                  IF( RESULT.EQ.0 .AND. IWORK( INDIWRK+I ).NE.
     $                IWORK( INDIWRK+NCLUSTERS+I ) ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9984 )
                     RESULT = 1
                  END IF
  130          CONTINUE
*
               IF( ICLUSTR( 1+IPREPAD+NCLUSTERS*2 ).NE.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9983 )
                  RESULT = 1
               END IF
            END IF
         END IF
*
*
         CALL IGAMX2D( DESCA( CTXT_ ), 'a', ' ', 1, 1, RESULT, 1, 1, 1,
     $                 -1, -1, 0 )
         IF( RESULT.NE.0 )
     $      GO TO 160
*
*     Compute eps * norm(A)
*
         IF( N.EQ.0 ) THEN
            EPSNORMA = EPS
         ELSE
            EPSNORMA = PZLANHE( 'I', UPLO, N, COPYA, IA, JA, DESCA,
     $                 RWORK )*EPS
         END IF
*
*     Note that a couple key variables get redefined in PZSEPCHK
*     as described by this table:
*
*     PZSEPTST name         PZSEPCHK name
*     -------------         -------------
*     COPYA                 A
*     Z                     Q
*     A                     C
*
*
         IF( LSAME( JOBZ, 'V' ) ) THEN
*
*     Perform the |AQ - QE| test
*
            CALL PDFILLPAD( DESCA( CTXT_ ), RSIZECHK, 1, RWORK,
     $                      RSIZECHK, IPREPAD, IPOSTPAD, 4.3D+0 )
*
            CALL PZSEPCHK( N, NZ, COPYA, IA, JA, DESCA,
     $                     MAX( ABSTOL+EPSNORMA, SAFMIN ), THRESH,
     $                     Z( 1+IPREPAD ), IA, JA, DESCZ,
     $                     A( 1+IPREPAD ), IA, JA, DESCA,
     $                     WNEW( 1+IPREPAD ), RWORK( 1+IPREPAD ),
     $                     RSIZECHK, TSTNRM, RES )
*
            CALL PDCHEKPAD( DESCA( CTXT_ ), 'PZSEPCHK-rWORK', RSIZECHK,
     $                      1, RWORK, RSIZECHK, IPREPAD, IPOSTPAD,
     $                      4.3D+0 )
*
            IF( RES.NE.0 )
     $         RESULT = 1
*
*     Perform the |QTQ - I| test
*
            CALL PDFILLPAD( DESCA( CTXT_ ), RSIZEQTQ, 1, RWORK,
     $                      RSIZEQTQ, IPREPAD, IPOSTPAD, 4.3D+0 )
*
*
            CALL PZSEPQTQ( N, NZ, THRESH, Z( 1+IPREPAD ), IA, JA, DESCZ,
     $                     A( 1+IPREPAD ), IA, JA, DESCA,
     $                     IWORK( 1+IPREPAD+1 ), ICLUSTR( 1+IPREPAD ),
     $                     GAP( 1+IPREPAD ), RWORK( IPREPAD+1 ),
     $                     RSIZEQTQ, QTQNRM, INFO, RES )
*
            CALL PDCHEKPAD( DESCA( CTXT_ ), 'PZSEPQTQ-rWORK', RSIZEQTQ,
     $                      1, RWORK, RSIZEQTQ, IPREPAD, IPOSTPAD,
     $                      4.3D+0 )
*
            IF( RES.NE.0 )
     $         RESULT = 1
*
            IF( INFO.NE.0 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9998 )INFO
               RESULT = 1
            END IF
         END IF
*
*     Check to make sure that we have the right eigenvalues
*
         IF( WKNOWN ) THEN
*
*        Set up MYIL if necessary
*
            MYIL = IL
*
            IF( LSAME( RANGE, 'V' ) ) THEN
               MYIL = 1
               MINIL = 1
               MAXIL = N - M + 1
            ELSE
               IF( LSAME( RANGE, 'A' ) ) THEN
                  MYIL = 1
               END IF
               MINIL = MYIL
               MAXIL = MYIL
            END IF
*
*     Find the largest difference between the computed
*     and expected eigenvalues
*
            MINERROR = NORMWIN
*
            DO 150 MYIL = MINIL, MAXIL
               MAXERROR = 0
*
*     Make sure that we aren't skipping any important eigenvalues
*
               MISSSMALLEST = .TRUE.
               IF( .NOT.LSAME( RANGE, 'V' ) .OR. ( MYIL.EQ.1 ) )
     $            MISSSMALLEST = .FALSE.
               IF( MISSSMALLEST .AND. ( WIN( MYIL-1 ).LT.VL+NORMWIN*
     $             FIVE*THRESH*EPS ) )MISSSMALLEST = .FALSE.
               MISSLARGEST = .TRUE.
               IF( .NOT.LSAME( RANGE, 'V' ) .OR. ( MYIL.EQ.MAXIL ) )
     $            MISSLARGEST = .FALSE.
               IF( MISSLARGEST .AND. ( WIN( MYIL+M ).GT.VU-NORMWIN*FIVE*
     $             THRESH*EPS ) )MISSLARGEST = .FALSE.
               IF( .NOT.MISSSMALLEST ) THEN
                  IF( .NOT.MISSLARGEST ) THEN
*
*     Make sure that the eigenvalues that we report are OK
*
                     DO 140 I = 1, M
                        ERROR = ABS( WIN( I+MYIL-1 )-WNEW( I+IPREPAD ) )
                        MAXERROR = MAX( MAXERROR, ERROR )
  140                CONTINUE
*
                     MINERROR = MIN( MAXERROR, MINERROR )
                  END IF
               END IF
  150       CONTINUE
*
*
*           If JOBZ = 'V' and RANGE='A', we might be comparing
*           against our estimate of what the eigenvalues ought to
*           be, rather than comparing against what PxHEEVX computed
*           last time around, so we have to be more generous.
*
            IF( LSAME( JOBZ, 'V' ) .AND. LSAME( RANGE, 'A' ) ) THEN
               IF( MINERROR.GT.NORMWIN*FIVE*FIVE*THRESH*EPS ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 )MINERROR, NORMWIN
                  RESULT = 1
               END IF
            ELSE
               IF( MINERROR.GT.NORMWIN*FIVE*THRESH*EPS ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 )MINERROR, NORMWIN
                  RESULT = 1
               END IF
            END IF
         END IF
*
*
*     Make sure that the IL, IU, VL and VU were not altered
*
         IF( IL.NE.OLDIL .OR. IU.NE.OLDIU .OR. VL.NE.OLDVL .OR. VU.NE.
     $       OLDVU ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9982 )
            RESULT = 1
         END IF
*
         IF( LSAME( JOBZ, 'N' ) .AND. ( NZ.NE.OLDNZ ) ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9981 )
            RESULT = 1
         END IF
*
      END IF
*
*     All processes should report the same result
*
      CALL IGAMX2D( DESCA( CTXT_ ), 'a', ' ', 1, 1, RESULT, 1, 1, 1, -1,
     $              -1, 0 )
*
  160 CONTINUE
*
*
      RETURN
*
 9999 FORMAT( 'PZHEEVX returned INFO=', I7 )
 9998 FORMAT( 'PZSEPQTQ returned INFO=', I7 )
 9997 FORMAT( 'PZSEPSUBTST minerror =', D11.2, ' normwin=', D11.2 )
 9996 FORMAT( 'PZHEEVX returned INFO=', I7,
     $      ' despite adequate workspace' )
 9995 FORMAT( 'ICLUSTR(1).NE.0 but mod(INFO/2,2).NE.1' )
 9994 FORMAT( 'M not in the range 0 to N' )
 9993 FORMAT( 'M not equal to N' )
 9992 FORMAT( 'M not equal to IU-IL+1' )
 9991 FORMAT( 'M not equal to NZ' )
 9990 FORMAT( 'NZ > M' )
 9989 FORMAT( 'NZ < M' )
 9988 FORMAT( 'NZ not equal to M' )
 9987 FORMAT( 'Different processes return different values for M' )
 9986 FORMAT( 'Different processes return different eigenvalues' )
 9985 FORMAT( 'Different processes return ',
     $      'different numbers of clusters' )
 9984 FORMAT( 'Different processes return different clusters' )
 9983 FORMAT( 'ICLUSTR not zero terminated' )
 9982 FORMAT( 'IL, IU, VL or VU altered by PZHEEVX' )
 9981 FORMAT( 'NZ altered by PZHEEVX with JOBZ=N' )
*
*     End of PZSEPSUBTST
*
      END
