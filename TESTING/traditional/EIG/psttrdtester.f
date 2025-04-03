      SUBROUTINE PSTTRDTESTER( IAM, NPROCS, CHECK, NOUT, THRESH, NVAL,
     $                         NMAT, MEM, TOTMEM, KPASS, KFAIL, KSKIP )
*
*  -- ScaLAPACK test routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     February 24, 2000
*
*     .. Scalar Arguments ..
      LOGICAL            CHECK
      INTEGER            IAM, KFAIL, KPASS, KSKIP, NMAT, NOUT, NPROCS,
     $                   TOTMEM
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            NVAL( * )
      REAL               MEM( * )
*     ..
*
*     Purpose
*     =======
*
*     PSTTRDTESTER tests PSSYTTRD
*
*     Arguments
*     =========
*
*     IAM     (local input) INTEGER
*     The local process number
*
*     NPROCS  (global input) INTEGER
*     The number of processors
*
*     CHECK   (global input) LOGICAL
*     Specifies whether the user wants to check the answer
*
*     NOUT    (local input) INTEGER
*     File descriptor
*
*     THRESH  (global input) REAL
*     Acceptable error threshold
*
*     NVAL    (global input) INTEGER array dimension NMAT
*     The matrix sizes to test
*
*     NMAT    (global input) INTEGER
*     The number of matrix sizes to test
*
*     MEM     (local input) REAL array dimension MEMSIZ
*     Where:
*       MEMSIZ = TOTMEM / REALSZ
*
*     TOTMEM  (global input) INTEGER
*     Number of bytes in MEM
*
*     KPASS   (local input/output) INTEGER
*     The number of tests which passed.  Only relevant on
*     processor 0.
*
*     KFAIL   (local input/output) INTEGER
*     The number of tests which failed.  Only relevant on
*     processor 0.
*
*     KSKIP   (local input/output) INTEGER
*     The number of tests which were skipped.  Only relevant on
*     processor 0.
*
*     ================================================================
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            REALSZ
      REAL               PADVAL
      PARAMETER          ( REALSZ = 4, PADVAL = -9923.0E+0 )
      INTEGER            TIMETESTS
      PARAMETER          ( TIMETESTS = 11 )
      INTEGER            TESTS
      PARAMETER          ( TESTS = 8 )
      INTEGER            MINTIMEN
      PARAMETER          ( MINTIMEN = 8 )
*     ..
*     .. Local Scalars ..
      LOGICAL            TIME
      CHARACTER          UPLO
      CHARACTER*6        PASSED
      INTEGER            DUMMY, IASEED, ICTXT, IMIDPAD, INFO, IPA, IPD,
     $                   IPE, IPOSTPAD, IPREPAD, IPT, IPW, ITEMP, J, K,
     $                   LCM, LWMIN, MAXTESTS, MEMSIZ, MYCOL, MYROW, N,
     $                   NB, NDIAG, NGRIDS, NN, NOFFD, NP, NPCOL, NPROW,
     $                   NPS, NQ, SPLITSTIMED, WORKSIZ, WORKTRD
      REAL               ANORM, FRESID
      DOUBLE PRECISION   NOPS, TMFLOPS
*     ..
*     .. Local Arrays ..
      INTEGER            ANBTEST( TESTS ), ANBTIME( TIMETESTS ),
     $                   BALTEST( TESTS ), BALTIME( TIMETESTS ),
     $                   DESCA( DLEN_ ), DESCD( DLEN_ ), IERR( 1 ),
     $                   INTERTEST( TESTS ), INTERTIME( TIMETESTS ),
     $                   PNBTEST( TESTS ), PNBTIME( TIMETESTS ),
     $                   TWOGEMMTEST( TESTS ), TWOGEMMTIME( TIMETESTS )
      DOUBLE PRECISION   CTIME( 100 ), WTIME( 100 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, DESCINIT,
     $                   IGEBR2D, IGEBS2D, IGSUM2D, PSCHEKPAD,
     $                   PSFILLPAD, PSLAFCHK, PSLATRAN, PSMATGEN,
     $                   PSSYTDRV, PSSYTTRD, SLBOOT, SLCOMBINE, SLTIMER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, ILCM, NUMROC, PJLAENV
      REAL               PSLANSY
      EXTERNAL           LSAME, ICEIL, ILCM, NUMROC, PJLAENV, PSLANSY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX, REAL, SQRT
*     ..
*
*     .. Scalars in Common ..
      INTEGER            ANB, BALANCED, BCKBLOCK, GSTBLOCK, INTERLEAVE,
     $                   LLTBLOCK, MINSZ, PNB, TIMEINTERNALS, TIMING,
     $                   TRSBLOCK, TWOGEMMS
*     ..
*     .. Common blocks ..
      COMMON             / BLOCKSIZES / GSTBLOCK, LLTBLOCK, BCKBLOCK,
     $                   TRSBLOCK
      COMMON             / MINSIZE / MINSZ
      COMMON             / PJLAENVTIMING / TIMING
      COMMON             / TAILOREDOPTS / PNB, ANB, INTERLEAVE,
     $                   BALANCED, TWOGEMMS
      COMMON             / TIMECONTROL / TIMEINTERNALS
*     ..
*     .. Data statements ..
      DATA               BALTIME / 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0 /
      DATA               INTERTIME / 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1 /
      DATA               TWOGEMMTIME / 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 /
      DATA               ANBTIME / 16, 16, 16, 16, 16, 8, 8, 32, 32, 16,
     $                   16 /
      DATA               PNBTIME / 32, 32, 32, 32, 32, 32, 32, 32, 32,
     $                   16, 64 /
      DATA               BALTEST / 0, 0, 0, 0, 1, 1, 1, 1 /
      DATA               INTERTEST / 0, 0, 1, 1, 0, 0, 1, 1 /
      DATA               TWOGEMMTEST / 0, 1, 0, 1, 0, 1, 0, 1 /
      DATA               ANBTEST / 1, 2, 3, 16, 1, 2, 3, 16 /
      DATA               PNBTEST / 1, 16, 8, 1, 16, 8, 1, 16 /
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*
      IASEED = 100
      SPLITSTIMED = 0
      NB = 1
      UPLO = 'L'
      MEMSIZ = TOTMEM / REALSZ
*
*     Print headings
*
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9995 )
         WRITE( NOUT, FMT = 9994 )
         WRITE( NOUT, FMT = 9993 )
         WRITE( NOUT, FMT = * )
      END IF
*
*     Loop over different process grids
*
      NGRIDS = INT( SQRT( REAL( NPROCS ) ) )
*
      DO 30 NN = 1, NGRIDS
*
         NPROW = NN
         NPCOL = NN
         IERR( 1 ) = 0
*
*        Define process grid
*
         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
         CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*        Go to bottom of loop if this case doesn't use my process
*
         IF( MYROW.GE.NPROW .OR. MYCOL.GE.NPCOL )
     $      GO TO 30
*
         DO 20 J = 1, NMAT
*
            N = NVAL( J )
*
*           Make sure matrix information is correct
*
            IERR( 1 ) = 0
            IF( N.LT.1 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9999 )'MATRIX', 'N', N
               IERR( 1 ) = 1
            END IF
*
*           Make sure no one had error
*
            CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
            IF( IERR( 1 ).GT.0 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 )'matrix'
               KSKIP = KSKIP + 1
               GO TO 20
            END IF
*
*           Loop over different blocking sizes
*
            IF( N.GT.MINTIMEN ) THEN
*
*              For timing tests, we perform one or two extra tests.
*              Both of these extra tests are performed with the
*              default values for the performance tuning parameters.
*              The second extra test (which is only performed if
*              split times are non-zero) is performed with timeinternals
*              set to 1 (which forces barrier syncs between many
*              phases of the computation).
*
               TIME = .TRUE.
               MAXTESTS = TIMETESTS + 2
            ELSE
               TIME = .FALSE.
               MAXTESTS = TESTS
            END IF
*
*
            DO 10 K = 1, MAXTESTS
               TIMEINTERNALS = 0
               IF( TIME ) THEN
                  IF( K.GE.MAXTESTS-1 ) THEN
*
*                    For the last two timings, we let pjlaenv set
*                    the execution path values.  These dummy
*                    initializations aren't really necessary,
*                    but they illustrate the fact that these values are
*                    set in xpjlaenv.  The dummy call to pjlaenv
*                    has the side effect of setting ANB.
*
                     MINSZ = -13
                     BALANCED = -13
                     INTERLEAVE = -13
                     TWOGEMMS = -13
                     ANB = -13
                     PNB = -13
                     TIMING = 1
                     DUMMY = PJLAENV( ICTXT, 3, 'PSSYTTRD', 'L', 0, 0,
     $                       0, 0 )
                     IF( K.EQ.MAXTESTS )
     $                  TIMEINTERNALS = 1
                  ELSE
                     TIMING = 0
                     MINSZ = 1
                     BALANCED = BALTIME( K )
                     INTERLEAVE = INTERTIME( K )
                     TWOGEMMS = TWOGEMMTIME( K )
                     ANB = ANBTIME( K )
                     PNB = PNBTIME( K )
                  END IF
               ELSE
                  TIMING = 0
                  MINSZ = 1
                  BALANCED = BALTEST( K )
                  INTERLEAVE = INTERTEST( K )
                  TWOGEMMS = TWOGEMMTEST( K )
                  ANB = ANBTEST( K )
                  PNB = PNBTEST( K )
               END IF
*
*              Skip the last test (with timeinternals = 1) if
*              PSSYTTRD is not collecting the split times.
*
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  CALL IGEBS2D( ICTXT, 'All', ' ', 1, 1, SPLITSTIMED,
     $                          1 )
               ELSE
                  CALL IGEBR2D( ICTXT, 'All', ' ', 1, 1, SPLITSTIMED, 1,
     $                          0, 0 )
               END IF
*
*
               IF( SPLITSTIMED.EQ.0 .AND. K.EQ.MAXTESTS )
     $            GO TO 10
*
*              The following hack tests to make sure that PNB need not
*              be the same on all processes.  (Provided that PNB is set
*              to 1 in the TRD.dat file.)
*
               IF( PNB.EQ.1 )
     $            PNB = 1 + IAM
*
*              Padding constants
*
               NP = NUMROC( N, NB, MYROW, 0, NPROW )
               NQ = NUMROC( N, NB, MYCOL, 0, NPCOL )
               IF( CHECK ) THEN
                  IPREPAD = MAX( NB, NP )
                  IMIDPAD = NB
                  IPOSTPAD = MAX( NB, NQ )
               ELSE
                  IPREPAD = 0
                  IMIDPAD = 0
                  IPOSTPAD = 0
               END IF
*
*              Initialize the array descriptor for the matrix A
*
*
               CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, ICTXT,
     $                        MAX( 1, NP )+IMIDPAD, IERR( 1 ) )
*
               CALL DESCINIT( DESCD, 1, N, NB, NB, 0, 0, ICTXT, 1,
     $                        INFO )
*
*              Check all processes for an error
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
               IF( IERR( 1 ).LT.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 )'descriptor'
                  KSKIP = KSKIP + 1
                  GO TO 10
               END IF
*
*              Assign pointers into MEM for SCALAPACK arrays, A is
*              allocated starting at position MEM( IPREPAD+1 )
*
               NDIAG = NQ
               IF( LSAME( UPLO, 'U' ) ) THEN
                  NOFFD = NQ
               ELSE
                  NOFFD = NUMROC( N-1, NB, MYCOL, 0, NPCOL )
               END IF
*
               IPA = IPREPAD + 1
               IPD = IPA + DESCA( LLD_ )*NQ + IPOSTPAD + IPREPAD
               IPE = IPD + NDIAG + IPOSTPAD + IPREPAD
               IPT = IPE + NOFFD + IPOSTPAD + IPREPAD
               IPW = IPT + NQ + IPOSTPAD + IPREPAD
*
*              Calculate the amount of workspace required for the
*              reduction
*
               NPS = MAX( NUMROC( N, 1, 0, 0, NPROW ), 2*ANB )
               LWMIN = 2*( ANB+1 )*( 4*NPS+2 ) + NPS
*
               WORKTRD = LWMIN + IPOSTPAD
               WORKSIZ = WORKTRD
*
*              Figure the amount of workspace required by the check
*
               IF( CHECK ) THEN
                  ITEMP = 2*NQ + NP
                  IF( NPROW.NE.NPCOL ) THEN
                     LCM = ILCM( NPROW, NPCOL )
                     ITEMP = NB*ICEIL( ICEIL( NP, NB ), LCM / NPROW ) +
     $                       ITEMP
                  END IF
                  ITEMP = MAX( ITEMP, 2*( NB+NP )*NB )
                  WORKSIZ = MAX( LWMIN, ITEMP ) + IPOSTPAD
               END IF
*
*              Check for adequate memory for problem size
*
               IERR( 1 ) = 0
               IF( IPW+WORKSIZ.GT.MEMSIZ ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9996 )'Tridiagonal reduction',
     $               ( IPW+WORKSIZ )*REALSZ
                  IERR( 1 ) = 1
               END IF
*
*              Check all processes for an error
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
               IF( IERR( 1 ).GT.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 )'MEMORY'
                  KSKIP = KSKIP + 1
                  GO TO 10
               END IF
*
*
*
*              Generate the matrix A
*
               CALL PSMATGEN( ICTXT, 'Hemm', 'N', DESCA( M_ ),
     $                        DESCA( N_ ), DESCA( MB_ ), DESCA( NB_ ),
     $                        MEM( IPA ), DESCA( LLD_ ), DESCA( RSRC_ ),
     $                        DESCA( CSRC_ ), IASEED, 0, NP, 0, NQ,
     $                        MYROW, MYCOL, NPROW, NPCOL )
*
*
*              Need Infinity-norm of A for checking
*
               IF( CHECK ) THEN
                  CALL PSFILLPAD( ICTXT, NP, NQ, MEM( IPA-IPREPAD ),
     $                            DESCA( LLD_ ), IPREPAD, IPOSTPAD,
     $                            PADVAL )
                  CALL PSFILLPAD( ICTXT, NDIAG, 1, MEM( IPD-IPREPAD ),
     $                            NDIAG, IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSFILLPAD( ICTXT, NOFFD, 1, MEM( IPE-IPREPAD ),
     $                            NOFFD, IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSFILLPAD( ICTXT, NQ, 1, MEM( IPT-IPREPAD ), NQ,
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  ANORM = PSLANSY( 'I', UPLO, N, MEM( IPA ), 1, 1,
     $                    DESCA, MEM( IPW ) )
                  CALL PSCHEKPAD( ICTXT, 'PSLANSY', NP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSLANSY', WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSFILLPAD( ICTXT, WORKTRD-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKTRD-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
               END IF
*
               CALL SLBOOT
               CALL BLACS_BARRIER( ICTXT, 'All' )
               CALL SLTIMER( 1 )
*
*              Reduce to symmetric tridiagonal form
*
               CALL PSSYTTRD( UPLO, N, MEM( IPA ), 1, 1, DESCA,
     $                        MEM( IPD ), MEM( IPE ), MEM( IPT ),
     $                        MEM( IPW ), LWMIN, INFO )
*
               CALL SLTIMER( 1 )
*
               IF( CHECK ) THEN
*
*                 Check for memory overwrite
*
                  CALL PSCHEKPAD( ICTXT, 'PSSYTTRD', NP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTTRD', NDIAG, 1,
     $                            MEM( IPD-IPREPAD ), NDIAG, IPREPAD,
     $                            IPOSTPAD, PADVAL )
*
                  CALL PSCHEKPAD( ICTXT, 'PSSYTTRDc', NOFFD, 1,
     $                            MEM( IPE-IPREPAD ), NOFFD, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTTRDd', NQ, 1,
     $                            MEM( IPT-IPREPAD ), NQ, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTTRDe', WORKTRD-IPOSTPAD,
     $                            1, MEM( IPW-IPREPAD ),
     $                            WORKTRD-IPOSTPAD, IPREPAD, IPOSTPAD,
     $                            PADVAL )
                  CALL PSFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
*
*                 Compute fctres = ||A - QTQ'|| / (||A|| * N * eps)
*
                  CALL PSSYTDRV( UPLO, N, MEM( IPA ), 1, 1, DESCA,
     $                           MEM( IPD ), MEM( IPE ), MEM( IPT ),
     $                           MEM( IPW ), IERR( 1 ) )
*
*                 TTRD does not preserve the upper triangular part of A.
*                 The following call to PSLATRAN means that we only
*                 check the lower triangular part of A - QTQ'
*
                  CALL PSLATRAN( N, 1, MEM( IPA ), 1, 1, DESCA,
     $                           MEM( IPW ) )
                  CALL PSLAFCHK( 'Hemm', 'No', N, N, MEM( IPA ), 1, 1,
     $                           DESCA, IASEED, ANORM, FRESID,
     $                           MEM( IPW ) )
*
*                 Check for memory overwrite
*
                  CALL PSCHEKPAD( ICTXT, 'PSSYTDRVf', NP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTDRVg', NDIAG, 1,
     $                            MEM( IPD-IPREPAD ), NDIAG, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTDRVh', NOFFD, 1,
     $                            MEM( IPE-IPREPAD ), NOFFD, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTDRVi', WORKSIZ-IPOSTPAD,
     $                            1, MEM( IPW-IPREPAD ),
     $                            WORKSIZ-IPOSTPAD, IPREPAD, IPOSTPAD,
     $                            PADVAL )
*
*                 Test residual and detect NaN result
*
                  IF( FRESID.LE.THRESH .AND. FRESID-FRESID.EQ.
     $                0.0E+0 .AND. IERR( 1 ).EQ.0 ) THEN
                     KPASS = KPASS + 1
                     PASSED = 'PASSED'
                  ELSE
                     IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $                  WRITE( NOUT, FMT = 9991 )FRESID
                     KFAIL = KFAIL + 1
                     PASSED = 'FAILED'
*
*
                  END IF
*
*
                  IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 .AND. IERR( 1 ).NE.0 )
     $               WRITE( NOUT, FMT = * )'D or E copies incorrect ...'
               ELSE
*
*                 Don't perform the checking, only the timing operation
*
                  KPASS = KPASS + 1
                  FRESID = FRESID - FRESID
                  PASSED = 'BYPASS'
               END IF
*
*              Gather maximum of all CPU and WALL clock timings
*
               CALL SLCOMBINE( ICTXT, 'All', '>', 'W', 50, 1, WTIME )
               CALL SLCOMBINE( ICTXT, 'All', '>', 'C', 50, 1, CTIME )
*
*              Print results
*
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
*
*                 TRD requires 16/3 N^3 floating point operations
*
                  NOPS = DBLE( N )
                  NOPS = ( 16.0D+0 / 3.0D+0 )*NOPS**3
                  NOPS = NOPS / 1.0D+6
*
*                 Print WALL time
*
                  IF( WTIME( 1 ).GT.0.0D+0 ) THEN
                     TMFLOPS = NOPS / WTIME( 1 )
                  ELSE
                     TMFLOPS = 0.0D+0
                  END IF
                  IF( WTIME( 1 ).GE.0.0D+0 )
     $               WRITE( NOUT, FMT = 9992 )'WALL', N, INTERLEAVE,
     $               TWOGEMMS, BALANCED, ANB, PNB, NPROW*NPCOL,
     $               WTIME( 1 ), TMFLOPS, FRESID, PASSED
*
*                 Print CPU time
*
                  IF( CTIME( 1 ).GT.0.0D+0 ) THEN
                     TMFLOPS = NOPS / CTIME( 1 )
                  ELSE
                     TMFLOPS = 0.0D+0
                  END IF
                  IF( CTIME( 1 ).GE.0.0D+0 )
     $               WRITE( NOUT, FMT = 9992 )'CPU ', N, INTERLEAVE,
     $               TWOGEMMS, BALANCED, ANB, PNB, NPROW*NPCOL,
     $               CTIME( 1 ), TMFLOPS, FRESID, PASSED
*
*
*                 If split times were collected (in PSSYttrd.f), print
*                 them out.
*
                  IF( WTIME( 13 )+WTIME( 15 )+WTIME( 16 ).GT.0.0D+0 .OR.
     $                CTIME( 13 )+CTIME( 15 )+CTIME( 16 ).GT.0.0D+0 )
     $                 THEN
                     SPLITSTIMED = 1
                  END IF
                  IF( SPLITSTIMED.EQ.1 ) THEN
                     WRITE( NOUT, FMT = 9990 )WTIME( 10 ), WTIME( 11 ),
     $                  WTIME( 12 ), WTIME( 13 ), WTIME( 14 ),
     $                  WTIME( 15 )
                     WRITE( NOUT, FMT = 9989 )WTIME( 16 ), WTIME( 17 ),
     $                  WTIME( 18 ), WTIME( 19 ), WTIME( 20 ),
     $                  WTIME( 21 )
*
                     WRITE( NOUT, FMT = 9988 )CTIME( 10 ), CTIME( 11 ),
     $                  CTIME( 12 ), CTIME( 13 ), CTIME( 14 ),
     $                  CTIME( 15 )
                     WRITE( NOUT, FMT = 9987 )CTIME( 16 ), CTIME( 17 ),
     $                  CTIME( 18 ), CTIME( 19 ), CTIME( 20 ),
     $                  CTIME( 21 )
                     WRITE( NOUT, FMT = 9986 )N, NPROW*NPCOL, PNB, ANB,
     $                  INTERLEAVE, BALANCED, TWOGEMMS, TIMEINTERNALS
                  END IF
               END IF
   10       CONTINUE
   20    CONTINUE
*
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            IF( SPLITSTIMED.EQ.1 ) THEN
               WRITE( NOUT, FMT = 9985 )
               WRITE( NOUT, FMT = 9984 )
               WRITE( NOUT, FMT = 9983 )
               WRITE( NOUT, FMT = 9982 )
               WRITE( NOUT, FMT = 9981 )
               WRITE( NOUT, FMT = 9980 )
               WRITE( NOUT, FMT = 9979 )
               WRITE( NOUT, FMT = 9978 )
               WRITE( NOUT, FMT = 9977 )
               WRITE( NOUT, FMT = 9976 )
               WRITE( NOUT, FMT = 9975 )
               WRITE( NOUT, FMT = 9974 )
               WRITE( NOUT, FMT = 9973 )
            END IF
         END IF
*
*
         CALL BLACS_GRIDEXIT( ICTXT )
   30 CONTINUE
      RETURN
*
 9999 FORMAT( 'ILLEGAL ', A6, ': ', A5, ' = ', I3,
     $      '; It should be at least 1' )
 9998 FORMAT( 'ILLEGAL GRID: nprow*npcol = ', I4, '. It can be at most',
     $      I4 )
 9997 FORMAT( 'Bad ', A6, ' parameters: going on to next test case.' )
 9996 FORMAT( 'Unable to perform ', A, ': need TOTMEM of at least',
     $      I11 )
*
 9995 FORMAT( 'PSSYTTRD, tailored reduction to tridiagonal form, test.'
     $       )
 9994 FORMAT( 'TIME N     int 2gm bal anb pnb prcs TRD Time ',
     $      '     MFLOPS Residual  CHECK' )
 9993 FORMAT( '---- ----  --- --- --- --- --- ---- -------- ',
     $      '----------- -------- ------' )
 9992 FORMAT( A4, 1X, I5, 1X, I3, 1X, I3, 1X, I3, 1X, I3, 1X, I3, 1X,
     $      I5, 1X, F9.2, 1X, F11.2, 1X, F8.2, 1X, A6 )
 9991 FORMAT( '||A - Q*T*Q''|| / (||A|| * N * eps) = ', G25.7 )
 9990 FORMAT( 'wsplit1=[wsplit1;', F9.2, 1X, F9.2, 1X, F9.2, 1X, F9.2,
     $      1X, F9.2, 1X, F9.2, ' ];' )
 9989 FORMAT( 'wsplit2=[wsplit2;', F9.2, 1X, F9.2, 1X, F9.2, 1X, F9.2,
     $      1X, F9.2, 1X, F9.2, ' ];' )
 9988 FORMAT( 'csplit1=[csplit1;', F9.2, 1X, F9.2, 1X, F9.2, 1X, F9.2,
     $      1X, F9.2, 1X, F9.2, ' ];' )
 9987 FORMAT( 'csplit2=[csplit2;', F9.2, 1X, F9.2, 1X, F9.2, 1X, F9.2,
     $      1X, F9.2, 1X, F9.2, ' ];' )
 9986 FORMAT( 'size_opts=[size_opts;', I4, 1X, I4, 1X, I4, 1X, I4, 1X,
     $      I4, 1X, I4, 1X, I4, 1X, I4, 1X, ' ];' )
 9985 FORMAT( 'N=1; NPROCS=2; PNB=3; ANB=4; INTERLEAVE=5; BALANCED=6;',
     $      ' TWOGEMMS=7; TIMEINTERNALS=8;' )
 9984 FORMAT( 'S1_OVERHEAD = 1; % Should be mainly cost of barrier' )
 9983 FORMAT( 'S1_BARRIER = 2; % Cost of barrier' )
 9982 FORMAT( 'S1_UPDCURCOL = 3; % Update the current column' )
 9981 FORMAT( 'S1_HOUSE = 4; % Compute the householder vector' )
 9980 FORMAT( 'S1_SPREAD = 5; % Spread across' )
 9979 FORMAT( 'S1_TRANSPOSE = 6; % Transpose' )
 9978 FORMAT( 'S2_UPDCURBLK = 1; % Update the current block column' )
 9977 FORMAT( 'S2_TRMVT = 2; % TRMVT v = A * h; vt = ht * A'' ' )
 9976 FORMAT( 'S2_UPD_V = 3; %  v = v + V * HT * h + H * VT * h ' )
 9975 FORMAT( 'S2_TRANS_SUM = 4; %  v = v + vt'' ' )
 9974 FORMAT( 'S2_DOT = 5; %  c = v'' * h ' )
 9973 FORMAT( 'S2_R2K = 6; %  A = A - v * h'' - h * v'' ' )
*
*
*     End of PSTTRDTESTER
*
      END
