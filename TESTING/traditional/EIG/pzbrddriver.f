      PROGRAM PZBRDDRIVER
*
*  -- ScaLAPACK testing driver (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 13, 2000
*
*  Purpose
*  =======
*
*  PZBRDDRIVER is the main test program for the COMPLEX*16
*  ScaLAPACK BRD (bidiagonal reduction) routines.
*
*  The program must be driven by a short data file.  An annotated
*  example of a data file can be obtained by deleting the first 3
*  characters from the following 13 lines:
*  'ScaLAPACK BRD computation input file'
*  'PVM machine'
*  'BRD.out'       output file name
*  6               device out
*  3               number of problems sizes
*  16 20 18        values of M
*  16 18 20        values of N
*  3               number of NB's
*  2 3 5           values of NB
*  7               number of process grids (ordered pairs of P & Q)
*  1 2 1 4 2 3 8   values of P
*  1 2 4 1 3 2 1   values of Q
*  1.0             threshold
*
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
*  ZPLXSZ   INTEGER, default = 16 bytes.
*           INTGSZ and ZPLXSZ indicate the length in bytes on the
*           given platform for an integer and a double precision
*           complex.
*  MEM      COMPLEX*16 array, dimension ( TOTMEM / ZPLXSZ )
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
      INTEGER            MEMSIZ, NTESTS, TOTMEM, ZPLXSZ, DBLESZ
      COMPLEX*16         PADVAL
      PARAMETER          ( TOTMEM = 2000000, ZPLXSZ = 16, DBLESZ = 8,
     $                     MEMSIZ = TOTMEM / ZPLXSZ, NTESTS = 20,
     $                     PADVAL = ( -9923.0D+0, -9923.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            CHECK
      CHARACTER*6        PASSED
      CHARACTER*80       OUTFILE
      INTEGER            I, IAM, IASEED, ICTXT, IMIDPAD, INFO, IPA, IPD,
     $                   IPE, IPOSTPAD, IPREPAD, IPTP, IPTQ, IPW, J, K,
     $                   KFAIL, KPASS, KSKIP, KTESTS, LWORK, M, MAXMN,
     $                   MINMN, MNP, MNQ, MP, MYCOL, MYROW, N, NB,
     $                   NDIAG, NGRIDS, NMAT, NNB, NOFFD, NOUT, NPCOL,
     $                   NPROCS, NPROW, NQ, WORKBRD, WORKSIZ
      REAL               THRESH
      DOUBLE PRECISION   ANORM, FRESID, NOPS, TMFLOPS
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), IERR( 1 ), NBVAL( NTESTS ),
     $                   MVAL( NTESTS ), NVAL( NTESTS ),
     $                   PVAL( NTESTS ), QVAL( NTESTS )
      DOUBLE PRECISION   CTIME( 1 ), WTIME( 1 )
      COMPLEX*16         MEM( MEMSIZ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_EXIT, BLACS_GET,
     $                   BLACS_GRIDEXIT, BLACS_GRIDINFO, BLACS_GRIDINIT,
     $                   BLACS_PINFO, DESCINIT, IGSUM2D, PZCHEKPAD,
     $                   PZBRDINFO, PZFILLPAD, PZLAFCHK,
     $                   PZMATGEN, PZGEBDRV, PZGEBRD, SLBOOT,
     $                   SLCOMBINE, SLTIMER
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, NUMROC
      DOUBLE PRECISION   PZLANGE
      EXTERNAL           ICEIL, NUMROC, PZLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
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
      CALL PZBRDINFO( OUTFILE, NOUT, NMAT, MVAL, NTESTS, NVAL, NTESTS,
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
     $         WRITE( NOUT, FMT = 9998 ) NPROW*NPCOL, NPROCS
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
         IF( MYROW.GE.NPROW .OR. MYCOL.GE.NPCOL )
     $      GO TO 30
*
*        Go to bottom of loop if this case doesn't use my process
*
         DO 20 J = 1, NMAT
*
            M = MVAL( J )
            N = NVAL( J )
*
*           Make sure matrix information is correct
*
            IERR( 1 ) = 0
            IF( M.LT.1 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9999 ) 'MATRIX', 'M', M
               IERR( 1 ) = 1
            ELSE IF( N.LT.1 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9999 ) 'MATRIX', 'N', N
               IERR( 1 ) = 1
            END IF
*
*           Make sure no one had error
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
*              Padding constants
*
               MP = NUMROC( M, NB, MYROW, 0, NPROW )
               NQ = NUMROC( N, NB, MYCOL, 0, NPCOL )
               MNP = NUMROC( MIN( M, N ), NB, MYROW, 0, NPROW )
               MNQ = NUMROC( MIN( M, N ), NB, MYCOL, 0, NPCOL )
               IF( CHECK ) THEN
                  IPREPAD  = MAX( NB, MP )
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
               CALL DESCINIT( DESCA, M, N, NB, NB, 0, 0, ICTXT,
     $                        MAX( 1, MP )+IMIDPAD, IERR( 1 ) )
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
               IF( M.GE.N ) THEN
                  NDIAG = MNQ
                  NOFFD = MNP
                  NDIAG = ICEIL( DBLESZ*NDIAG, ZPLXSZ )
                  NOFFD = ICEIL( DBLESZ*NOFFD, ZPLXSZ )
               ELSE
                  NDIAG = MNP
                  NOFFD = NUMROC( MIN( M, N )-1, NB, MYCOL, 0, NPCOL )
                  NDIAG = ICEIL( DBLESZ*NDIAG, ZPLXSZ )
                  NOFFD = ICEIL( DBLESZ*NOFFD, ZPLXSZ )
               END IF
*
               IPA  = IPREPAD + 1
               IPD  = IPA + DESCA( LLD_ )*NQ + IPOSTPAD + IPREPAD
               IPE  = IPD + NDIAG + IPOSTPAD + IPREPAD
               IPTQ = IPE + NOFFD + IPOSTPAD + IPREPAD
               IPTP = IPTQ + MNQ + IPOSTPAD + IPREPAD
               IPW  = IPTP + MNP + IPOSTPAD + IPREPAD
*
*              Calculate the amount of workspace required for the
*              reduction
*
               LWORK = NB*( MP+NQ+1 ) + NQ
               WORKBRD = LWORK + IPOSTPAD
               WORKSIZ = WORKBRD
*
*              Figure the amount of workspace required by the check
*
               IF( CHECK ) THEN
                  WORKSIZ = MAX( LWORK, 2*NB*( MP+NQ+NB ) ) + IPOSTPAD
               END IF
*
*              Check for adequate memory for problem size
*
               IERR( 1 ) = 0
               IF( IPW+WORKSIZ.GT.MEMSIZ ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9996 ) 'Bidiagonal reduction',
     $                      ( IPW+WORKSIZ )*ZPLXSZ
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
*              Generate the matrix A
*
               CALL PZMATGEN( ICTXT, 'No', 'No', DESCA( M_ ),
     $                        DESCA( N_ ), DESCA( MB_ ), DESCA( NB_ ),
     $                        MEM( IPA ), DESCA( LLD_ ), DESCA( RSRC_ ),
     $                        DESCA( CSRC_ ), IASEED, 0, MP, 0, NQ,
     $                        MYROW, MYCOL, NPROW, NPCOL )
*
*              Need Infinity-norm of A for checking
*
               IF( CHECK ) THEN
                  CALL PZFILLPAD( ICTXT, MP, NQ, MEM( IPA-IPREPAD ),
     $                            DESCA( LLD_ ), IPREPAD, IPOSTPAD,
     $                            PADVAL )
                  CALL PZFILLPAD( ICTXT, NDIAG, 1, MEM( IPD-IPREPAD ),
     $                            NDIAG, IPREPAD, IPOSTPAD, PADVAL )
                  CALL PZFILLPAD( ICTXT, NOFFD, 1, MEM( IPE-IPREPAD ),
     $                            NOFFD, IPREPAD, IPOSTPAD, PADVAL )
                  CALL PZFILLPAD( ICTXT, MNQ, 1, MEM( IPTQ-IPREPAD ),
     $                            MNQ, IPREPAD, IPOSTPAD, PADVAL )
                  CALL PZFILLPAD( ICTXT, MNP, 1, MEM( IPTP-IPREPAD ),
     $                            MNP, IPREPAD, IPOSTPAD, PADVAL )
                  CALL PZFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  ANORM = PZLANGE( 'I', M, N, MEM( IPA ), 1, 1, DESCA,
     $                             MEM( IPW ) )
                  CALL PZCHEKPAD( ICTXT, 'PZLANGE', MP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PZCHEKPAD( ICTXT, 'PZLANGE', WORKSIZ-IPOSTPAD,
     $                            1, MEM( IPW-IPREPAD ),
     $                            WORKSIZ-IPOSTPAD, IPREPAD, IPOSTPAD,
     $                            PADVAL )
                  CALL PZFILLPAD( ICTXT, WORKBRD-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKBRD-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
               END IF
*
               CALL SLBOOT()
               CALL BLACS_BARRIER( ICTXT, 'All' )
               CALL SLTIMER( 1 )
*
*              Reduce to bidiagonal form
*
               CALL PZGEBRD( M, N, MEM( IPA ), 1, 1, DESCA, MEM( IPD ),
     $                       MEM( IPE ), MEM( IPTQ ), MEM( IPTP ),
     $                       MEM( IPW ), LWORK, INFO )
*
               CALL SLTIMER( 1 )
*
               IF( CHECK ) THEN
*
*                 Check for memory overwrite
*
                  CALL PZCHEKPAD( ICTXT, 'PZGEBRD', MP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PZCHEKPAD( ICTXT, 'PZGEBRD', NDIAG, 1,
     $                            MEM( IPD-IPREPAD ), NDIAG, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PZCHEKPAD( ICTXT, 'PZGEBRD', NOFFD, 1,
     $                            MEM( IPE-IPREPAD ), NOFFD, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PZCHEKPAD( ICTXT, 'PZGEBRD', MNQ, 1,
     $                            MEM( IPTQ-IPREPAD ), MNQ, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PZCHEKPAD( ICTXT, 'PZGEBRD', MNP, 1,
     $                            MEM( IPTP-IPREPAD ), MNP, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PZCHEKPAD( ICTXT, 'PZGEBRD', WORKBRD-IPOSTPAD,
     $                            1, MEM( IPW-IPREPAD ),
     $                            WORKBRD-IPOSTPAD, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PZFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
*
*                 Compute fctres = ||A-Q*B*P|| / (||A|| * N * eps)
*
                  CALL PZGEBDRV( M, N, MEM( IPA ), 1, 1, DESCA,
     $                           MEM( IPD ), MEM( IPE ), MEM( IPTQ ),
     $                           MEM( IPTP ), MEM( IPW ), IERR( 1 ) )
                  CALL PZLAFCHK( 'No', 'No', M, N, MEM( IPA ), 1, 1,
     $                           DESCA, IASEED, ANORM, FRESID,
     $                           MEM( IPW ) )
*
*                 Check for memory overwrite
*
                  CALL PZCHEKPAD( ICTXT, 'PZGEBDRV', MP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PZCHEKPAD( ICTXT, 'PZGEBDRV', NDIAG, 1,
     $                            MEM( IPD-IPREPAD ), NDIAG, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PZCHEKPAD( ICTXT, 'PZGEBDRV', NOFFD, 1,
     $                            MEM( IPE-IPREPAD ), NOFFD, IPREPAD,
     $                            IPOSTPAD, PADVAL )
                  CALL PZCHEKPAD( ICTXT, 'PZGEBDRV', WORKSIZ-IPOSTPAD,
     $                           1, MEM( IPW-IPREPAD ),
     $                           WORKSIZ-IPOSTPAD, IPREPAD,
     $                           IPOSTPAD, PADVAL )
*
*                 Test residual and detect NaN result
*
                  IF( FRESID.LE.THRESH .AND. FRESID-FRESID.EQ.0.0D+0
     $                .AND. IERR( 1 ).EQ.0 ) THEN
                     KPASS = KPASS + 1
                     PASSED = 'PASSED'
                  ELSE
                     IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $                  WRITE( NOUT, FMT = 9986 ) FRESID
*
                     KFAIL = KFAIL + 1
                     PASSED = 'FAILED'
                  END IF
*
                  IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 .AND. IERR( 1 ).NE.0 )
     $               WRITE( NOUT, FMT = * )
     $                     'D or E copies incorrect ...'
               ELSE
*
*                 Don't perform the checking, only the timing operation
*
                  KPASS = KPASS + 1
                  FRESID = FRESID - FRESID
                  PASSED = 'BYPASS'
*
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
*                 BRD requires 32/3 N^3 floating point operations
*
                  MAXMN = MAX( M, N )
                  MINMN = MIN( M, N )
                  NOPS = 16.0D+0 * DBLE( MINMN ) * DBLE( MINMN ) *
     $                   ( DBLE( MAXMN ) - DBLE( MINMN ) / 3.0D+0 )
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
     $               WRITE( NOUT, FMT = 9993 ) 'WALL', M, N, NB, NPROW,
     $                      NPCOL, WTIME( 1 ), TMFLOPS, FRESID, PASSED
*
*                 Print CPU time
*
                  IF( CTIME( 1 ).GT.0.0D+0 ) THEN
                     TMFLOPS = NOPS / CTIME( 1 )
                  ELSE
                     TMFLOPS = 0.0D+0
                  END IF
                  IF( CTIME( 1 ).GE.0.0D+0 )
     $               WRITE( NOUT, FMT = 9993 ) 'CPU ', M, N, NB, NPROW,
     $                      NPCOL, CTIME( 1 ), TMFLOPS, FRESID, PASSED
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
         IF( NOUT.NE.6 .AND. NOUT.NE.0 ) CLOSE ( NOUT )
      END IF
*
      CALL BLACS_EXIT( 0 )
*
 9999 FORMAT( 'ILLEGAL ', A6, ': ', A5, ' = ', I3,
     $        '; It should be at least 1' )
 9998 FORMAT( 'ILLEGAL GRID: nprow*npcol = ', I4, '. It can be at most',
     $        I4 )
 9997 FORMAT( 'Bad ', A6, ' parameters: going on to next test case.' )
 9996 FORMAT( 'Unable to perform ', A, ': need TOTMEM of at least',
     $      I11 )
 9995 FORMAT( 'TIME      M      N  NB     P     Q  BRD Time ',
     $        '     MFLOPS Residual  CHECK' )
 9994 FORMAT( '---- ------ ------ --- ----- ----- --------- ',
     $        '----------- -------- ------' )
 9993 FORMAT( A4, 1X, I6, 1X, I6, 1X, I3, 1X, I5, 1X, I5, 1X, F9.2, 1X,
     $        F11.2, 1X, F8.2, 1X, A6 )
 9992 FORMAT( 'Finished', I4, ' tests, with the following results:' )
 9991 FORMAT( I5, ' tests completed and passed residual checks.' )
 9990 FORMAT( I5, ' tests completed without checking.' )
 9989 FORMAT( I5, ' tests completed and failed residual checks.' )
 9988 FORMAT( I5, ' tests skipped because of illegal input values.' )
 9987 FORMAT( 'END OF TESTS.' )
 9986 FORMAT( '||A - Q*B*P|| / (||A|| * N * eps) = ', G25.7 )
*
      STOP
*
*     End of PZBRDDRIVER
*
      END
