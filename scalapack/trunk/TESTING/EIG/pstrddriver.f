      PROGRAM PSTRDDRIVER
*
*  -- ScaLAPACK testing driver (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     October 15, 1999
*
*  Purpose
*  ========
*
*  PSTRDDRIVER is the main test program for the REAL
*  SCALAPACK TRD (symmetric tridiagonal reduction) routines.
*
*  The program must be driven by a short data file.  An annotated
*  example of a data file can be obtained by deleting the first 3
*  characters from the following 13 lines:
*  'ScaLAPACK TRD computation input file'
*  'PVM machine'
*  'TRD.out'       output file name
*  6               device out
*  'L'             define Lower or Upper
*  3               number of problems sizes
*  5 31 201        values of N
*  3               number of NB's
*  2 3 5           values of NB
*  7               number of process grids (ordered pairs of P & Q)
*  1 2 1 4 2 3 8   values of P
*  1 2 4 1 3 2 1   values of Q
*  1.0             threshold
*
*  Internal Parameters
*  ===================
*
*  TOTMEM   INTEGER, default = 2000000
*           TOTMEM is a machine-specific parameter indicating the
*           maximum amount of available memory in bytes.
*           The user should customize TOTMEM to his platform.  Remember
*           to leave room in memory for the operating system, the BLACS
*           buffer, etc.  For example, on a system with 8 MB of memory
*           per process (e.g., one processor on an Intel iPSC/860), the
*           parameters we use are TOTMEM=6200000 (leaving 1.8 MB for OS,
*           code, BLACS buffer, etc).  However, for PVM, we usually set
*           TOTMEM = 2000000.  Some experimenting with the maximum value
*           of TOTMEM may be required.
*
*  INTGSZ   INTEGER, default = 4 bytes.
*  REALSZ   INTEGER, default = 4 bytes.
*           INTGSZ and REALSZ indicate the length in bytes on the
*           given platform for an integer and a single precision real.
*  MEM      REAL array, dimension ( TOTMEM / REALSZ )
*
*           All arrays used by SCALAPACK routines are allocated from
*           this array and referenced by pointers.  The integer IPA,
*           for example, is a pointer to the starting element of MEM for
*           the matrix A.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            REALSZ, TOTMEM, MEMSIZ, NTESTS
      REAL               PADVAL
      PARAMETER          ( REALSZ = 4, TOTMEM = 2000000,
     $                   MEMSIZ = TOTMEM / REALSZ, NTESTS = 20,
     $                   PADVAL = -9923.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            CHECK
      CHARACTER          UPLO
      CHARACTER*6        PASSED
      CHARACTER*80       OUTFILE
      INTEGER            I, IAM, IASEED, ICTXT, IMIDPAD, INFO, IPA, IPD,
     $                   IPE, IPOSTPAD, IPREPAD, IPT, IPW, ITEMP, J, K,
     $                   KFAIL, KPASS, KSKIP, KTESTS, LCM, LWORK, MYCOL,
     $                   MYROW, N, NB, NDIAG, NGRIDS, NMAT, NNB, NOFFD,
     $                   NOUT, NP, NPCOL, NPROCS, NPROW, NQ, WORKSIZ,
     $                   WORKTRD
      REAL               ANORM, FRESID, THRESH
      DOUBLE PRECISION   NOPS, TMFLOPS
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), IERR( 1 ), NBVAL( NTESTS ),
     $                   NVAL( NTESTS ), PVAL( NTESTS ), QVAL( NTESTS )
      REAL               MEM( MEMSIZ )
      DOUBLE PRECISION   CTIME( 1 ), WTIME( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_EXIT, BLACS_GET,
     $                   BLACS_GRIDEXIT, BLACS_GRIDINFO, BLACS_GRIDINIT,
     $                   BLACS_PINFO, DESCINIT, IGSUM2D, PSCHEKPAD,
     $                   PSFILLPAD, PSLAFCHK, PSMATGEN, PSSYTDRV,
     $                   PSSYTRD, PSTRDINFO, PSTTRDTESTER, SLBOOT,
     $                   SLCOMBINE, SLTIMER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, ILCM, NUMROC
      REAL               PSLANSY
      EXTERNAL           LSAME, ICEIL, ILCM, NUMROC, PSLANSY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Data statements ..
      DATA               KTESTS, KPASS, KFAIL, KSKIP / 4*0 /
*     ..
*     .. Executable Statements ..
*
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )STOP
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
      IASEED = 100
      CALL PSTRDINFO( OUTFILE, NOUT, UPLO, NMAT, NVAL, NTESTS, NNB,
     $                NBVAL, NTESTS, NGRIDS, PVAL, NTESTS, QVAL, NTESTS,
     $                THRESH, MEM, IAM, NPROCS )
      CHECK = ( THRESH.GE.0.0E+0 )
*
*     Print headings
*
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9995 )
         WRITE( NOUT, FMT = 9994 )
         WRITE( NOUT, FMT = * )
      END IF
*
*     Loop over different process grids
*
      DO 30 I = 1, NGRIDS
*
         NPROW = PVAL( I )
         NPCOL = QVAL( I )
*
*        Make sure grid information is correct
*
         IERR( 1 ) = 0
         IF( NPROW.LT.1 ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9999 )'GRID', 'nprow', NPROW
            IERR( 1 ) = 1
         ELSE IF( NPCOL.LT.1 ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9999 )'GRID', 'npcol', NPCOL
            IERR( 1 ) = 1
         ELSE IF( NPROW*NPCOL.GT.NPROCS ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )NPROW*NPCOL, NPROCS
            IERR( 1 ) = 1
         END IF
*
         IF( IERR( 1 ).GT.0 ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 )'grid'
            KSKIP = KSKIP + 1
            GO TO 30
         END IF
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
            DO 10 K = 1, NNB
*
               NB = NBVAL( K )
*
*              Make sure nb is legal
*
               IERR( 1 ) = 0
               IF( NB.LT.1 ) THEN
                  IERR( 1 ) = 1
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9999 )'NB', 'NB', NB
               END IF
*
*              Check all processes for an error
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
               IF( IERR( 1 ).GT.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 )'NB'
                  KSKIP = KSKIP + 1
                  GO TO 10
               END IF
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
               CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, ICTXT,
     $                        MAX( 1, NP )+IMIDPAD, IERR( 1 ) )
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
               LWORK = MAX( NB*( NP+1 ), 3*NB )
               WORKTRD = LWORK + IPOSTPAD
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
                  WORKSIZ = MAX( LWORK, ITEMP ) + IPOSTPAD
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
*              Generate the matrix A
*
               CALL PSMATGEN( ICTXT, 'Symm', 'N', DESCA( M_ ),
     $                        DESCA( N_ ), DESCA( MB_ ), DESCA( NB_ ),
     $                        MEM( IPA ), DESCA( LLD_ ), DESCA( RSRC_ ),
     $                        DESCA( CSRC_ ), IASEED, 0, NP, 0, NQ,
     $                        MYROW, MYCOL, NPROW, NPCOL )
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
               CALL PSSYTRD( UPLO, N, MEM( IPA ), 1, 1, DESCA,
     $                       MEM( IPD ), MEM( IPE ), MEM( IPT ),
     $                       MEM( IPW ), LWORK, INFO )
*
               CALL SLTIMER( 1 )
*
               IF( CHECK ) THEN
*
*                 Check for memory overwrite
*
                  CALL PSCHEKPAD( ICTXT, 'PSSYTRD', NP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTRD', NDIAG, 1,
     $                            MEM( IPD-IPREPAD ), NDIAG, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTRD', NOFFD, 1,
     $                            MEM( IPE-IPREPAD ), NOFFD, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTRD', NQ, 1,
     $                            MEM( IPT-IPREPAD ), NQ, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTRD', WORKTRD-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKTRD-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
*
*                 Compute fctres = ||A - QTQ'|| / (||A|| * N * eps)
*
                  CALL PSSYTDRV( UPLO, N, MEM( IPA ), 1, 1, DESCA,
     $                           MEM( IPD ), MEM( IPE ), MEM( IPT ),
     $                           MEM( IPW ), IERR( 1 ) )
                  CALL PSLAFCHK( 'Symm', 'No', N, N, MEM( IPA ), 1, 1,
     $                           DESCA, IASEED, ANORM, FRESID,
     $                           MEM( IPW ) )
*
*                 Check for memory overwrite
*
                  CALL PSCHEKPAD( ICTXT, 'PSSYTDRV', NP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTDRV', NDIAG, 1,
     $                            MEM( IPD-IPREPAD ), NDIAG, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTDRV', NOFFD, 1,
     $                            MEM( IPE-IPREPAD ), NOFFD, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSSYTDRV', WORKSIZ-IPOSTPAD,
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
     $                  WRITE( NOUT, FMT = 9986 )FRESID
                     KFAIL = KFAIL + 1
                     PASSED = 'FAILED'
                  END IF
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
               CALL SLCOMBINE( ICTXT, 'All', '>', 'W', 1, 1, WTIME )
               CALL SLCOMBINE( ICTXT, 'All', '>', 'C', 1, 1, CTIME )
*
*              Print results
*
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
*
*                 TRD requires 4/3 N^3 floating point operations
*
                  NOPS = DBLE( N )
*
                  NOPS = ( 4.0D+0 / 3.0D+0 )*NOPS**3
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
     $               WRITE( NOUT, FMT = 9993 )'WALL', UPLO, N, NB,
     $               NPROW, NPCOL, WTIME( 1 ), TMFLOPS, FRESID, PASSED
*
*                 Print CPU time
*
                  IF( CTIME( 1 ).GT.0.0D+0 ) THEN
                     TMFLOPS = NOPS / CTIME( 1 )
                  ELSE
                     TMFLOPS = 0.0D+0
                  END IF
                  IF( CTIME( 1 ).GE.0.0D+0 )
     $               WRITE( NOUT, FMT = 9993 )'CPU ', UPLO, N, NB,
     $               NPROW, NPCOL, CTIME( 1 ), TMFLOPS, FRESID, PASSED
               END IF
   10       CONTINUE
   20    CONTINUE
*
         CALL BLACS_GRIDEXIT( ICTXT )
   30 CONTINUE
*
      CALL PSTTRDTESTER( IAM, NPROCS, CHECK, NOUT, THRESH, NVAL, NMAT,
     $                   MEM, TOTMEM, KPASS, KFAIL, KSKIP )
*
*     Print ending messages and close output file
*
      IF( IAM.EQ.0 ) THEN
         KTESTS = KPASS + KFAIL + KSKIP
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9992 )KTESTS
         IF( CHECK ) THEN
            WRITE( NOUT, FMT = 9991 )KPASS
            WRITE( NOUT, FMT = 9989 )KFAIL
         ELSE
            WRITE( NOUT, FMT = 9990 )KPASS
         END IF
         WRITE( NOUT, FMT = 9988 )KSKIP
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9987 )
         IF( NOUT.NE.6 .AND. NOUT.NE.0 )
     $      CLOSE ( NOUT )
      END IF
*
      CALL BLACS_EXIT( 0 )
*
 9999 FORMAT( 'ILLEGAL ', A6, ': ', A5, ' = ', I3,
     $      '; It should be at least 1' )
 9998 FORMAT( 'ILLEGAL GRID: nprow*npcol = ', I4, '. It can be at most',
     $      I4 )
 9997 FORMAT( 'Bad ', A6, ' parameters: going on to next test case.' )
 9996 FORMAT( 'Unable to perform ', A, ': need TOTMEM of at least',
     $      I11 )
 9995 FORMAT( 'TIME UPLO      N  NB     P     Q  TRD Time ',
     $      '     MFLOPS Residual  CHECK' )
 9994 FORMAT( '---- ---- ------ --- ----- ----- --------- ',
     $      '----------- -------- ------' )
 9993 FORMAT( A4, 1X, A4, 1X, I6, 1X, I3, 1X, I5, 1X, I5, 1X, F9.2, 1X,
     $      F11.2, 1X, F8.2, 1X, A6 )
 9992 FORMAT( 'Finished', I4, ' tests, with the following results:' )
 9991 FORMAT( I5, ' tests completed and passed residual checks.' )
 9990 FORMAT( I5, ' tests completed without checking.' )
 9989 FORMAT( I5, ' tests completed and failed residual checks.' )
 9988 FORMAT( I5, ' tests skipped because of illegal input values.' )
 9987 FORMAT( 'END OF TESTS.' )
 9986 FORMAT( '||A - Q*T*Q''|| / (||A|| * N * eps) = ', G25.7 )
*
      STOP
*
*     End of PSTRDDRIVER
*
      END
