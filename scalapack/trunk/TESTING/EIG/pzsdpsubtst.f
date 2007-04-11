      SUBROUTINE PZSDPSUBTST( WKNOWN, UPLO, N, THRESH, ABSTOL, A, COPYA,
     $                        Z, IA, JA, DESCA, WIN, WNEW, IPREPAD,
     $                        IPOSTPAD, WORK, LWORK, RWORK, LRWORK,
     $                        LWORK1, IWORK, LIWORK, RESULT, TSTNRM,
     $                        QTQNRM, NOUT )
*
*  -- ScaLAPACK testing routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     February 28, 2000
*
*     .. Scalar Arguments ..
      LOGICAL            WKNOWN
      CHARACTER          UPLO
      INTEGER            IA, IPOSTPAD, IPREPAD, JA, LIWORK, LRWORK,
     $                   LWORK, LWORK1, N, NOUT, RESULT
      DOUBLE PRECISION   ABSTOL, QTQNRM, THRESH, TSTNRM
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), IWORK( * )
      DOUBLE PRECISION   RWORK( * ), WIN( * ), WNEW( * )
      COMPLEX*16         A( * ), COPYA( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PZSDPSUBTST calls PZHEEVD and then tests the output of
*  PZHEEVD
*  If JOBZ = 'V' then the following two tests are performed:
*     |AQ -QL| / (abstol + eps * norm(A) ) < N*THRESH
*     |QT * Q - I| / eps * norm(A) < N*THRESH
*  If WKNOWN then
*     we check to make sure that the eigenvalues match expectations
*     i.e. |WIN - WNEW(1+IPREPAD)| / (eps * |WIN|) < THRESH
*     where WIN is the array of eigenvalues as computed by
*     PZHEEVD when eigenvectors are requested
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
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          Size of the matrix to be tested.  (global size)
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
*            See PZHEEVD for a description of block cyclic layout.
*          The test matrix, which is then modified by PZHEEVD
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
*          The eigenvalues as copmuted by this call to PZHEEVD
*          If JOBZ <> 'V' or RANGE <> 'A' these eigenvalues are
*          compared against those in WIN().
*          WNEW has already been padded front and back,
*          use WNEW(1+IPREPAD)
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
*          The amount of DOUBLE PRECISION workspace to pass to PZHEEVD
*
*  IWORK   (local workspace) INTEGER array, dimension (LIWORK)
*          IWORK has already been padded front and back,
*          use IWORK(1+IPREPAD)
*
*  LIWORK  (local input) INTEGER
*          The length of the array IWORK after padding.
*
*  RESULT  (global output) INTEGER
*          The result of this call to PZHEEVD
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
      DOUBLE PRECISION               PADVAL, FIVE, NEGONE
      PARAMETER          ( PADVAL = 13.5285D+0, FIVE = 5.0D+0,
     $                   NEGONE = -1.0D+0 )
      COMPLEX*16          CPADVAL
      PARAMETER          ( CPADVAL = ( 13.989D+0, 1.93D+0 ) )
      INTEGER            IPADVAL
      PARAMETER          ( IPADVAL = 927 )
      COMPLEX*16          CZERO, CONE, CNEGONE
      PARAMETER          ( CZERO = 0.0D+0, CONE = 1.0D+0,
     $                   CNEGONE = -1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IAM, INFO, ISIZEHEEVD, ISIZEHEEVX,
     $                   ISIZESUBTST, ISIZETST, MYCOL, MYROW, NP, NPCOL,
     $                   NPROW, NQ, RES, RSIZECHK, RSIZEHEEVD,
     $                   RSIZEHEEVX, RSIZEQTQ, RSIZESUBTST, RSIZETST,
     $                   SIZEHEEVD, SIZEHEEVX, SIZEMQRLEFT,
     $                   SIZEMQRRIGHT, SIZEQRF, SIZESUBTST, SIZETMS,
     $                   SIZETST
      DOUBLE PRECISION   EPS, EPSNORMA, ERROR, MAXERROR, MINERROR, NORM,
     $                   NORMWIN, SAFMIN, ULP
*     ..
*     .. Local Arrays ..
      INTEGER            ITMP( 2 )
*     ..
*     .. External Functions ..
*
      INTEGER            NUMROC
      DOUBLE PRECISION   PZLANGE, PZLANHE, PDLAMCH
      EXTERNAL           NUMROC, PZLANGE, PZLANHE, PDLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, ZLACPY, IGAMN2D, IGAMX2D,
     $                   PZCHEKPAD, PZFILLPAD, PZGEMM, PZHEEVD, PZLASET,
     $                   PZLASIZESEP, PZSEPCHK, PICHEKPAD, PIFILLPAD,
     $                   PDCHEKPAD, PDFILLPAD, SLBOOT, SLTIMER
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, DBLE
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      CALL PZLASIZESEP( DESCA, IPREPAD, IPOSTPAD, SIZEMQRLEFT,
     $                  SIZEMQRRIGHT, SIZEQRF, SIZETMS, RSIZEQTQ,
     $                  RSIZECHK, SIZEHEEVX, RSIZEHEEVX, ISIZEHEEVX,
     $                  SIZEHEEVD, RSIZEHEEVD, ISIZEHEEVD, SIZESUBTST,
     $                  RSIZESUBTST, ISIZESUBTST, SIZETST, RSIZETST,
     $                  ISIZETST )
*
      TSTNRM = NEGONE
      QTQNRM = NEGONE
      EPS = PDLAMCH( DESCA( CTXT_ ), 'Eps' )
      SAFMIN = PDLAMCH( DESCA( CTXT_ ), 'Safe min' )
*
      NORMWIN = SAFMIN / EPS
      IF( N.GE.1 )
     $   NORMWIN = MAX( ABS( WIN( 1+IPREPAD ) ),
     $             ABS( WIN( N+IPREPAD ) ), NORMWIN )
*
      DO 10 I = 1, LWORK1, 1
         RWORK( I+IPREPAD ) = 14.3D+0
   10 CONTINUE
      DO 20 I = 1, LIWORK, 1
         IWORK( I ) = 14
   20 CONTINUE
      DO 30 I = 1, LWORK, 1
         WORK( I+IPREPAD ) = ( 15.63D+0, 1.1D+0 )
   30 CONTINUE
*
      DO 40 I = 1, N
         WNEW( I+IPREPAD ) = 3.14159D+0
   40 CONTINUE
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
      IAM = 1
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $   IAM = 0
*
*     If this process is not involved in this test, bail out now
*
      RESULT = -3
      IF( MYROW.GE.NPROW .OR. MYROW.LT.0 )
     $   GO TO 60
      RESULT = 0
*
      NP = NUMROC( N, DESCA( MB_ ), MYROW, 0, NPROW )
      NQ = NUMROC( N, DESCA( NB_ ), MYCOL, 0, NPCOL )
*
      CALL ZLACPY( 'A', NP, NQ, COPYA, DESCA( LLD_ ), A( 1+IPREPAD ),
     $             DESCA( LLD_ ) )
*
      CALL PZFILLPAD( DESCA( CTXT_ ), NP, NQ, A, DESCA( LLD_ ), IPREPAD,
     $                IPOSTPAD, CPADVAL )
*
      CALL PZFILLPAD( DESCA( CTXT_ ), NP, NQ, Z, DESCA( LLD_ ), IPREPAD,
     $                IPOSTPAD, CPADVAL+1.0D+0 )
*
      CALL PDFILLPAD( DESCA( CTXT_ ), N, 1, WNEW, N, IPREPAD, IPOSTPAD,
     $                PADVAL+2.0D+0 )
*
      CALL PDFILLPAD( DESCA( CTXT_ ), LWORK1, 1, RWORK, LWORK1, IPREPAD,
     $                IPOSTPAD, PADVAL+4.0D+0 )
*
      CALL PIFILLPAD( DESCA( CTXT_ ), LIWORK, 1, IWORK, LIWORK, IPREPAD,
     $                IPOSTPAD, IPADVAL )
*
      CALL PZFILLPAD( DESCA( CTXT_ ), LWORK, 1, WORK, LWORK, IPREPAD,
     $                IPOSTPAD, CPADVAL+4.1D+0 )
*
      CALL SLBOOT
      CALL SLTIMER( 1 )
      CALL SLTIMER( 6 )
*
      CALL PZHEEVD( 'V', UPLO, N, A( 1+IPREPAD ), IA, JA, DESCA,
     $              WNEW( 1+IPREPAD ), Z( 1+IPREPAD ), IA, JA, DESCA,
     $              WORK( 1+IPREPAD ), SIZEHEEVD, RWORK( 1+IPREPAD ),
     $              LWORK1, IWORK( 1+IPREPAD ), LIWORK, INFO )
      CALL SLTIMER( 6 )
      CALL SLTIMER( 1 )
*
      IF( THRESH.LE.0 ) THEN
         RESULT = 0
      ELSE
         CALL PZCHEKPAD( DESCA( CTXT_ ), 'PZHEEVD-A', NP, NQ, A,
     $                   DESCA( LLD_ ), IPREPAD, IPOSTPAD, CPADVAL )
*
         CALL PZCHEKPAD( DESCA( CTXT_ ), 'PZHEEVD-Z', NP, NQ, Z,
     $                   DESCA( LLD_ ), IPREPAD, IPOSTPAD,
     $                   CPADVAL+1.0D+0 )
*
         CALL PDCHEKPAD( DESCA( CTXT_ ), 'PZHEEVD-WNEW', N, 1, WNEW, N,
     $                   IPREPAD, IPOSTPAD, PADVAL+2.0D+0 )
*
         CALL PDCHEKPAD( DESCA( CTXT_ ), 'PZHEEVD-rWORK', LWORK1, 1,
     $                   RWORK, LWORK1, IPREPAD, IPOSTPAD,
     $                   PADVAL+4.0D+0 )
*
         CALL PZCHEKPAD( DESCA( CTXT_ ), 'PZHEEVD-WORK', LWORK, 1, WORK,
     $                   LWORK, IPREPAD, IPOSTPAD, CPADVAL+4.1D+0 )
*
         CALL PICHEKPAD( DESCA( CTXT_ ), 'PZHEEVD-IWORK', LIWORK, 1,
     $                   IWORK, LIWORK, IPREPAD, IPOSTPAD, IPADVAL )
*
*     Check INFO
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
         ELSE IF( INFO.NE.0 ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9996 )INFO
            RESULT = 1
         END IF
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
*     Perform the |AQ - QE| test
*
         CALL PDFILLPAD( DESCA( CTXT_ ), RSIZECHK, 1, RWORK, RSIZECHK,
     $                   IPREPAD, IPOSTPAD, 4.3D+0 )
*
         CALL PZSEPCHK( N, N, COPYA, IA, JA, DESCA,
     $                  MAX( ABSTOL+EPSNORMA, SAFMIN ), THRESH,
     $                  Z( 1+IPREPAD ), IA, JA, DESCA, A( 1+IPREPAD ),
     $                  IA, JA, DESCA, WNEW( 1+IPREPAD ),
     $                  RWORK( 1+IPREPAD ), RSIZECHK, TSTNRM, RES )
*
         CALL PDCHEKPAD( DESCA( CTXT_ ), 'PZSDPCHK-rWORK', RSIZECHK, 1,
     $                   RWORK, RSIZECHK, IPREPAD, IPOSTPAD, 4.3D+0 )
*
         IF( RES.NE.0 ) THEN
            RESULT = 1
            WRITE( NOUT, FMT = 9995 )
         END IF
*
*     Perform the |QTQ - I| test
*
         CALL PDFILLPAD( DESCA( CTXT_ ), RSIZEQTQ, 1, RWORK, RSIZEQTQ,
     $                   IPREPAD, IPOSTPAD, 4.3D+0 )
*
*
         RES = 0
         ULP = PDLAMCH( DESCA( CTXT_ ), 'P' )
         CALL PZLASET( 'A', N, N, CZERO, CONE, A( 1+IPREPAD ), IA, JA,
     $                 DESCA )
         CALL PZGEMM( 'Conjugate transpose', 'N', N, N, N, CNEGONE,
     $                Z( 1+IPREPAD ), IA, JA, DESCA, Z( 1+IPREPAD ), IA,
     $                JA, DESCA, CONE, A( 1+IPREPAD ), IA, JA, DESCA )
         NORM = PZLANGE( '1', N, N, A( 1+IPREPAD ), IA, JA, DESCA,
     $          WORK( 1+IPREPAD ) )
         QTQNRM = NORM / ( DBLE( MAX( N, 1 ) )*ULP )
         IF( QTQNRM.GT.THRESH ) THEN
            RES = 1
         END IF
         CALL PDCHEKPAD( DESCA( CTXT_ ), 'PZSEPQTQ-rWORK', RSIZEQTQ, 1,
     $                   RWORK, RSIZEQTQ, IPREPAD, IPOSTPAD, 4.3D+0 )
*
         IF( RES.NE.0 ) THEN
            RESULT = 1
            WRITE( NOUT, FMT = 9994 )
         END IF
*
         IF( INFO.NE.0 ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )INFO
            RESULT = 1
         END IF
*
         IF( INFO.NE.0 ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )INFO
            RESULT = 1
         END IF
      END IF
*
*     Check to make sure that we have the right eigenvalues
*
      IF( WKNOWN .AND. N.GT.0 ) THEN
*
*     Find the largest difference between the computed
*     and expected eigenvalues
*
         MINERROR = NORMWIN
         MAXERROR = 0.0D+00
*
         DO 50 I = 1, N
            ERROR = ABS( WIN( I+IPREPAD )-WNEW( I+IPREPAD ) )
            MAXERROR = MAX( MAXERROR, ERROR )
   50    CONTINUE
         MINERROR = MIN( MAXERROR, MINERROR )
*
         IF( MINERROR.GT.NORMWIN*FIVE*THRESH*EPS ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 )MINERROR, NORMWIN
            RESULT = 1
         END IF
      END IF
*
*
*     All processes should report the same result
*
      CALL IGAMX2D( DESCA( CTXT_ ), 'a', ' ', 1, 1, RESULT, 1, 1, 1, -1,
     $              -1, 0 )
*
   60 CONTINUE
*
      RETURN
*
 9999 FORMAT( 'PZHEEVD returned INFO=', I7 )
 9998 FORMAT( 'PZSEPQTQ returned INFO=', I7 )
 9997 FORMAT( 'PZSDPSUBTST minerror =', D11.2, ' normwin=', D11.2 )
 9996 FORMAT( 'PZHEEVD returned INFO=', I7,
     $      ' despite adequate workspace' )
 9995 FORMAT( 'PZHEEVD failed the |AQ -QE| test' )
 9994 FORMAT( 'PZHEEVD failed the |QTQ -I| test' )
*
*     End of PZSDPSUBTST
*
      END
