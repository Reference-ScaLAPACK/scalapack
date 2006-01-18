      PROGRAM PSHRDDRIVER
*
*  -- ScaLAPACK testing driver (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 13, 2000
*
*  Purpose
*  =======
*
*  PSHRDDRIVER is the main test program for the REAL
*  ScaLAPACK HRD (Hessenberg Reduction) routines.
*
*  The program must be driven by a short data file.  An annotated
*  example of a data file can be obtained by deleting the first 3
*  characters from the following 14 lines:
*  'ScaLAPACK HRD input file'
*  'PVM machine'
*  'HRD.out'            output file name (if any)
*  6                    device out
*  2                    number of problems sizes
*  100 101              values of N
*  2   1                values of ILO
*  99  101              values of IHI
*  3                    number of NB's
*  2 3 5                values of NB
*  7                    number of process grids (ordered pairs of P & Q)
*  1 2 1 4 2 3 8        values of P
*  1 2 4 1 3 2 1        values of Q
*  3.0                  threshold
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
      INTEGER            MEMSIZ, NTESTS, REALSZ, TOTMEM
      REAL               PADVAL
      PARAMETER          ( REALSZ = 4, TOTMEM = 2000000,
     $                     MEMSIZ = TOTMEM / REALSZ, NTESTS = 20,
     $                     PADVAL = -9923.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            CHECK
      CHARACTER*6        PASSED
      CHARACTER*80       OUTFILE
      INTEGER            I, IAM, IASEED, ICTXT, IHI, IHIP, IHLP, IHLQ,
     $                   ILCOL, ILO, ILROW, INFO, INLQ, IMIDPAD, IPA,
     $                   IPT, IPW, IPOSTPAD, IPREPAD, ITEMP, J, K,
     $                   KFAIL, KPASS, KSKIP, KTESTS, LCM, LCMQ, LOFF,
     $                   LWORK, MYCOL, MYROW, N, NB, NGRIDS, NMAT, NNB,
     $                   NPROCS, NOUT, NP, NPCOL, NPROW, NQ, WORKHRD,
     $                   WORKSIZ
      REAL               ANORM, FRESID, THRESH
      DOUBLE PRECISION   NOPS, TMFLOPS
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), IERR( 1 ), NBVAL( NTESTS ),
     $                   NVAL( NTESTS ), NVHI( NTESTS ), NVLO( NTESTS ),
     $                   PVAL( NTESTS ), QVAL( NTESTS )
      REAL               MEM( MEMSIZ )
      DOUBLE PRECISION   CTIME( 1 ), WTIME( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_EXIT, BLACS_GET,
     $                   BLACS_GRIDEXIT, BLACS_GRIDINIT, BLACS_GRIDINFO,
     $                   DESCINIT, IGSUM2D, BLACS_PINFO, PSFILLPAD,
     $                   PSLAFCHK, PSGEHDRV, PSGEHRD,
     $                   PSHRDINFO, PSMATGEN, SLBOOT,
     $                   SLCOMBINE, SLTIMER
*     ..
*     .. External Functions ..
      INTEGER            ILCM, INDXG2P, NUMROC
      REAL               PSLANGE
      EXTERNAL           ILCM, INDXG2P, NUMROC, PSLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Data statements ..
      DATA               KTESTS, KPASS, KFAIL, KSKIP / 4*0 /
*     ..
*     .. Executable Statements ..
*
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
      IASEED = 100
      CALL PSHRDINFO( OUTFILE, NOUT, NMAT, NVAL, NVLO, NVHI, NTESTS,
     $                NNB, NBVAL, NTESTS, NGRIDS, PVAL, NTESTS, QVAL,
     $                NTESTS, THRESH, MEM, IAM, NPROCS )
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
     $         WRITE( NOUT, FMT = 9999 ) 'GRID', 'nprow', NPROW
            IERR( 1 ) = 1
         ELSE IF( NPCOL.LT.1 ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9999 ) 'GRID', 'npcol', NPCOL
            IERR( 1 ) = 1
         ELSE IF( NPROW*NPCOL.GT.NPROCS ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )NPROW*NPCOL, NPROCS
            IERR( 1 ) = 1
         END IF
*
         IF( IERR( 1 ).GT.0 ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 ) 'grid'
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
     $      GOTO 30
*
         DO 20 J = 1, NMAT
*
            N = NVAL( J )
            ILO = NVLO( J )
            IHI = NVHI( J )
*
*           Make sure matrix information is correct
*
            IERR( 1 ) = 0
            IF( N.LT.1 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9999 ) 'MATRIX', 'N', N
               IERR( 1 ) = 1
            END IF
*
*           Check all processes for an error
*
            CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
            IF( IERR( 1 ).GT.0 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'matrix'
               KSKIP = KSKIP + 1
               GO TO 20
            END IF
*
            DO 10 K = 1, NNB
               NB = NBVAL( K )
*
*              Make sure nb is legal
*
               IERR( 1 ) = 0
               IF( NB.LT.1 ) THEN
                  IERR( 1 ) = 1
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9999 ) 'NB', 'NB', NB
               END IF
*
*              Check all processes for an error
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
               IF( IERR( 1 ).GT.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 ) 'NB'
                  KSKIP = KSKIP + 1
                  GO TO 10
               END IF
*
               NP = NUMROC( N, NB, MYROW, 0, NPROW )
               NQ = NUMROC( N, NB, MYCOL, 0, NPCOL )
               IF( CHECK ) THEN
                  IPREPAD  = MAX( NB, NP )
                  IMIDPAD  = NB
                  IPOSTPAD = MAX( NB, NQ )
               ELSE
                  IPREPAD  = 0
                  IMIDPAD  = 0
                  IPOSTPAD = 0
               END IF
*
*              Initialize the array descriptor for the matrix A
*
               CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, ICTXT,
     $                        MAX( 1, NP ) + IMIDPAD, INFO )
*
*              Check all processes for an error
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
               IF( IERR( 1 ).LT.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 ) 'descriptor'
                  KSKIP = KSKIP + 1
                  GO TO 10
               END IF
*
*              Assign pointers into MEM for SCALAPACK arrays, A is
*              allocated starting at position MEM( IPREPAD+1 )
*
               IPA = IPREPAD + 1
               IPT = IPA + DESCA( LLD_ )*NQ + IPOSTPAD + IPREPAD
               IPW = IPT + NQ + IPOSTPAD + IPREPAD
*
*              Calculate the amount of workspace required for the
*              reduction
*
               IHIP = NUMROC( IHI, NB, MYROW, DESCA( RSRC_ ), NPROW )
               LOFF = MOD( ILO-1, NB )
               ILROW = INDXG2P( ILO, NB, MYROW, DESCA( RSRC_ ), NPROW )
               ILCOL = INDXG2P( ILO, NB, MYCOL, DESCA( CSRC_ ), NPCOL )
               IHLP = NUMROC( IHI-ILO+LOFF+1, NB, MYROW, ILROW, NPROW )
               INLQ = NUMROC( N-ILO+LOFF+1, NB, MYCOL, ILCOL, NPCOL )
               LWORK = NB*( NB + MAX( IHIP+1, IHLP+INLQ ) )
               WORKHRD = LWORK + IPOSTPAD
               WORKSIZ = WORKHRD
*
*              Figure the amount of workspace required by the check
*
               IF( CHECK ) THEN
                  LCM = ILCM( NPROW, NPCOL )
                  LCMQ = LCM / NPCOL
                  IHLQ = NUMROC( IHI-ILO+LOFF+1, NB, MYCOL, ILCOL,
     $                           NPCOL )
                  ITEMP = NB*MAX( IHLP+INLQ, IHLQ+MAX( IHIP,
     $                    IHLP+NUMROC( NUMROC( IHI-ILO+LOFF+1, NB, 0, 0,
     $                    NPCOL ), NB, 0, 0, LCMQ ) ) )
                  WORKSIZ = MAX( NB*NB + NB*IHLP + ITEMP, NB * NP ) +
     $                           IPOSTPAD
               END IF
*
*              Check for adequate memory for problem size
*
               IERR( 1 ) = 0
               IF( IPW+WORKSIZ.GT.MEMSIZ ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9996 ) 'Hessenberg reduction',
     $                      ( IPW+WORKSIZ )*REALSZ
                  IERR( 1 ) = 1
               END IF
*
*              Check all processes for an error
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
               IF( IERR( 1 ).GT.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 ) 'MEMORY'
                  KSKIP = KSKIP + 1
                  GO TO 10
               END IF
*
*              Generate A
*
               CALL PSMATGEN( ICTXT, 'No', 'No', DESCA( M_ ),
     $                        DESCA( N_ ), DESCA( MB_ ), DESCA( NB_ ),
     $                        MEM( IPA ), DESCA( LLD_ ), DESCA( RSRC_ ),
     $                        DESCA( CSRC_ ),
     $                        IASEED, 0, NP, 0, NQ, MYROW, MYCOL,
     $                        NPROW, NPCOL )
*
*              Need Infinity-norm of A for checking
*
               IF( CHECK ) THEN
                  CALL PSFILLPAD( ICTXT, NP, NQ, MEM( IPA-IPREPAD ),
     $                            DESCA( LLD_ ), IPREPAD, IPOSTPAD,
     $                            PADVAL )
                  CALL PSFILLPAD( ICTXT, NQ, 1, MEM( IPT-IPREPAD ),
     $                            NQ, IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  ANORM = PSLANGE( 'I', N, N, MEM( IPA ), 1, 1, DESCA,
     $                             MEM( IPW ) )
                  CALL PSCHEKPAD( ICTXT, 'PSLANGE', NP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSLANGE',
     $                            WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSFILLPAD( ICTXT, WORKHRD-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKHRD-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
               END IF
*
               CALL SLBOOT()
               CALL BLACS_BARRIER( ICTXT, 'All' )
               CALL SLTIMER( 1 )
*
*              Reduce Hessenberg form
*
               CALL PSGEHRD( N, ILO, IHI, MEM( IPA ), 1, 1, DESCA,
     $                       MEM( IPT ), MEM( IPW ), LWORK, INFO )
               CALL SLTIMER( 1 )
*
               IF( CHECK ) THEN
*
*                 Check for memory overwrite
*
                  CALL PSCHEKPAD( ICTXT, 'PSGEHRD', NP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSGEHRD', NQ, 1,
     $                            MEM( IPT-IPREPAD ), NQ, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSGEHRD', WORKHRD-IPOSTPAD,
     $                            1, MEM( IPW-IPREPAD ),
     $                            WORKHRD-IPOSTPAD, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PSFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
*
*                 Compute fctres = ||A - Q H Q'|| / (||A||*N*eps)
*
                  CALL PSGEHDRV( N, ILO, IHI, MEM( IPA ), 1, 1, DESCA,
     $                           MEM( IPT ), MEM( IPW ) )
                  CALL PSLAFCHK( 'No', 'No', N, N, MEM( IPA ), 1, 1,
     $                           DESCA, IASEED, ANORM, FRESID,
     $                           MEM( IPW ) )
*
*                 Check for memory overwrite
*
                  CALL PSCHEKPAD( ICTXT, 'PSGEHDRV', NP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSGEHDRV', NQ, 1,
     $                            MEM( IPT-IPREPAD ), NQ, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSGEHDRV',
     $                            WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
*
*                 Test residual and detect NaN result
*
                  IF( FRESID.LE.THRESH .AND. FRESID-FRESID.EQ.0.0E+0 )
     $                THEN
                     KPASS = KPASS + 1
                     PASSED = 'PASSED'
                  ELSE
                     IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $                  WRITE( NOUT, FMT = 9986 ) FRESID
                     KFAIL = KFAIL + 1
                     PASSED = 'FAILED'
                  END IF
               ELSE
*
*                 Don't perform the checking, only the timing operation
*
                  KPASS = KPASS + 1
                  FRESID = FRESID - FRESID
                  PASSED = 'BYPASS'
               END IF
*
*              Gather max. of all CPU and WALL clock timings
*
               CALL SLCOMBINE( ICTXT, 'All', '>', 'W', 1, 1, WTIME )
               CALL SLCOMBINE( ICTXT, 'All', '>', 'C', 1, 1, CTIME )
*
*              Print results
*
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
*
*                 HRD requires 10/3 * N^3 floating point ops. (flops)
*                 more precisely,
*                 HRD requires 4/3*(IHI-ILO)^3 + 2*IHI*(IHI-ILO)^2 flops
*
                  NOPS = DBLE( IHI-ILO )
                  NOPS = NOPS * NOPS *
     $                   ( 2.0D0*DBLE( IHI ) + (4.0D0/3.0D0)*NOPS )
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
     $               WRITE( NOUT, FMT = 9993 ) 'WALL',  N, ILO, IHI, NB,
     $                      NPROW, NPCOL, WTIME( 1 ), TMFLOPS, FRESID,
     $                      PASSED
*
*                 Print CPU time
*
                  IF( CTIME( 1 ).GT.0.0D+0 ) THEN
                     TMFLOPS = NOPS / CTIME( 1 )
                  ELSE
                     TMFLOPS = 0.0D+0
                  END IF
                  IF( CTIME( 1 ).GE.0.0D+0 )
     $               WRITE( NOUT, FMT = 9993 ) 'CPU ', N, ILO, IHI, NB,
     $                      NPROW, NPCOL, CTIME( 1 ), TMFLOPS, FRESID,
     $                      PASSED
               END IF
   10       CONTINUE
   20    CONTINUE
*
         CALL BLACS_GRIDEXIT( ICTXT )
   30 CONTINUE
*
*     Print ending messages and close output file
*
      IF( IAM.EQ.0 ) THEN
         KTESTS = KPASS + KFAIL + KSKIP
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9992 ) KTESTS
         IF( CHECK ) THEN
            WRITE( NOUT, FMT = 9991 ) KPASS
            WRITE( NOUT, FMT = 9989 ) KFAIL
         ELSE
            WRITE( NOUT, FMT = 9990 ) KPASS
         END IF
         WRITE( NOUT, FMT = 9988 ) KSKIP
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9987 )
         IF( NOUT.NE.6 .AND. NOUT.NE.0 )
     $      CLOSE( NOUT )
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
 9995 FORMAT( 'TIME      N    ILO    IHI  NB     P     Q  HRD Time ',
     $        '     MFLOPS Residual  CHECK' )
 9994 FORMAT( '---- ------ ------ ------ --- ----- ----- --------- ',
     $        '----------- -------- ------' )
 9993 FORMAT( A4, 1X, I6, 1X, I6, 1X, I6, 1X, I3, 1X, I5, 1X, I5, 1X,
     $        F9.2, 1X, F11.2, 1X, F8.2, 1X, A6 )
 9992 FORMAT( 'Finished', I4, ' tests, with the following results:' )
 9991 FORMAT( I5, ' tests completed and passed residual checks.' )
 9990 FORMAT( I5, ' tests completed without checking.' )
 9989 FORMAT( I5, ' tests completed and failed residual checks.' )
 9988 FORMAT( I5, ' tests skipped because of illegal input values.' )
 9987 FORMAT( 'END OF TESTS.' )
 9986 FORMAT( '||A - Q*H*Q''|| / (||A|| * N * eps) = ', G25.7 )
*
      STOP
*
*     End of PSHRDDRIVER
*
      END
