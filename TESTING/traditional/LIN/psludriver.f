      PROGRAM PSLUDRIVER
*
*  -- ScaLAPACK testing driver (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*  Purpose
*  ========
*
*  PSLUDRIVER is the main test program for the REAL
*  SCALAPACK LU routines.  This test driver performs an LU factorization
*  and solve. If the input matrix is non-square, only the factorization
*  is performed.  Condition estimation and iterative refinement are
*  optionally performed.
*
*  The program must be driven by a short data file.  An annotated
*  example of a data file can be obtained by deleting the first 3
*  characters from the following 18 lines:
*  'SCALAPACK, Version 2.0,  LU factorization input file'
*  'Intel iPSC/860 hypercube, gamma model.'
*  'LU.out'             output file name (if any)
*  6                    device out
*  1                    number of problems sizes
*  31 201               values of M
*  31 201               values of N
*  1                    number of NB's
*  2                    values of NB
*  1                    number of NRHS's
*  1                    values of NRHS
*  1                    number of NBRHS's
*  1                    values of NBRHS
*  1                    number of process grids (ordered pairs of P & Q)
*  2 1 4 2 3 8          values of P
*  2 4 1 3 2 1          values of Q
*  1.0                  threshold
*  T                    (T or F) Test Cond. Est. and Iter. Ref. Routines
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
      INTEGER            INTGSZ, MEMSIZ, NTESTS, REALSZ, TOTMEM
      REAL               PADVAL, ZERO
      PARAMETER          ( INTGSZ = 4, REALSZ = 4, TOTMEM = 2000000,
     $                     MEMSIZ = TOTMEM / REALSZ, NTESTS = 20,
     $                     PADVAL = -9923.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            CHECK, EST
      CHARACTER*6        PASSED
      CHARACTER*80       OUTFILE
      INTEGER            HH, I, IAM, IASEED, IBSEED, ICTXT, IMIDPAD,
     $                   INFO, IPA, IPA0, IPB, IPB0, IPBERR, IPFERR,
     $                   IPOSTPAD, IPPIV, IPREPAD, IPW, IPW2, J, K,
     $                   KFAIL, KK, KPASS, KSKIP, KTESTS, LCM, LCMQ,
     $                   LIPIV, LIWORK, LWORK, LW2, M, MAXMN,
     $                   MINMN, MP, MYCOL, MYRHS, MYROW, N, NB, NBRHS,
     $                   NGRIDS, NMAT, NNB, NNBR, NNR, NOUT, NP, NPCOL,
     $                   NPROCS, NPROW, NQ, NRHS, WORKSIZ
      REAL               ANORM, ANORM1, FRESID, RCOND, SRESID, SRESID2,
     $                   THRESH
      DOUBLE PRECISION   NOPS, TMFLOPS
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), IERR( 1 ),
     $                   MVAL( NTESTS ), NBRVAL( NTESTS ),
     $                   NBVAL( NTESTS ), NRVAL( NTESTS ),
     $                   NVAL( NTESTS ), PVAL( NTESTS ),
     $                   QVAL( NTESTS )
      REAL               MEM( MEMSIZ )
      DOUBLE PRECISION   CTIME( 2 ), WTIME( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_EXIT, BLACS_GET,
     $                   BLACS_GRIDEXIT, BLACS_GRIDINFO, BLACS_GRIDINIT,
     $                   BLACS_PINFO, DESCINIT, IGSUM2D, PSCHEKPAD,
     $                   PSFILLPAD, PSGECON, PSGERFS,
     $                   PSGETRF, PSGETRRV, PSGETRS,
     $                   PSLAFCHK, PSLASCHK, PSLUINFO,
     $                   PSMATGEN, SLBOOT, SLCOMBINE, SLTIMER
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, ILCM, NUMROC
      REAL               PSLANGE
      EXTERNAL           ICEIL, ILCM, NUMROC, PSLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Data Statements ..
      DATA               KFAIL, KPASS, KSKIP, KTESTS / 4*0 /
*     ..
*     .. Executable Statements ..
*
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
      IASEED = 100
      IBSEED = 200
      CALL PSLUINFO( OUTFILE, NOUT, NMAT, MVAL, NVAL, NTESTS, NNB,
     $               NBVAL, NTESTS, NNR, NRVAL, NTESTS, NNBR, NBRVAL,
     $               NTESTS, NGRIDS, PVAL, NTESTS, QVAL, NTESTS, THRESH,
     $               EST, MEM, IAM, NPROCS )
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
      DO 50 I = 1, NGRIDS
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
            GO TO 50
         END IF
*
*        Define process grid
*
         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
         CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*        Go to bottom of process grid loop if this case doesn't use my
*        process
*
         IF( MYROW.GE.NPROW .OR. MYCOL.GE.NPCOL )
     $      GO TO 50
*
         DO 40 J = 1, NMAT
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
*           Check all processes for an error
*
            CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
            IF( IERR( 1 ).GT.0 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'matrix'
               KSKIP = KSKIP + 1
               GO TO 40
            END IF
*
            DO 30 K = 1, NNB
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
                  GO TO 30
               END IF
*
*              Padding constants
*
               MP = NUMROC( M, NB, MYROW, 0, NPROW )
               NP = NUMROC( N, NB, MYROW, 0, NPROW )
               NQ = NUMROC( N, NB, MYCOL, 0, NPCOL )
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
*              Check all processes for an error
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
               IF( IERR( 1 ).LT.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 ) 'descriptor'
                  KSKIP = KSKIP + 1
                  GO TO 30
               END IF
*
*              Assign pointers into MEM for SCALAPACK arrays, A is
*              allocated starting at position MEM( IPREPAD+1 )
*
               IPA = IPREPAD+1
               IF( EST .AND. M.EQ.N ) THEN
                  IPA0 = IPA + DESCA( LLD_ )*NQ + IPOSTPAD + IPREPAD
                  IPPIV = IPA0 + DESCA( LLD_ )*NQ + IPOSTPAD + IPREPAD
               ELSE
                  IPPIV = IPA + DESCA( LLD_ )*NQ + IPOSTPAD + IPREPAD
               END IF
               LIPIV = ICEIL( INTGSZ*( MP+NB ), REALSZ )
               IPW = IPPIV + LIPIV + IPOSTPAD + IPREPAD
*
               IF( CHECK ) THEN
*
*                 Calculate the amount of workspace required by the
*                 checking routines PSLANGE, PSGETRRV, and
*                 PSLAFCHK
*
                  WORKSIZ = MAX( 2, NQ )
*
                  WORKSIZ = MAX( WORKSIZ, MP*DESCA( NB_ )+
     $                      NQ*DESCA( MB_ ) )
*
                  WORKSIZ = MAX( WORKSIZ, MP * DESCA( NB_ ) )
*
                  WORKSIZ = WORKSIZ + IPOSTPAD
*
               ELSE
*
                  WORKSIZ = IPOSTPAD
*
               END IF
*
*              Check for adequate memory for problem size
*
               IERR( 1 ) = 0
               IF( IPW+WORKSIZ.GT.MEMSIZ ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9996 ) 'factorization',
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
                  GO TO 30
               END IF
*
*              Generate matrix A of Ax = b
*
               CALL PSMATGEN( ICTXT, 'No transpose', 'No transpose',
     $                        DESCA( M_ ), DESCA( N_ ), DESCA( MB_ ),
     $                        DESCA( NB_ ), MEM( IPA ), DESCA( LLD_ ),
     $                        DESCA( RSRC_ ), DESCA( CSRC_ ), IASEED, 0,
     $                        MP, 0, NQ, MYROW, MYCOL, NPROW, NPCOL )
*
*              Calculate inf-norm of A for residual error-checking
*
               IF( CHECK ) THEN
                  CALL PSFILLPAD( ICTXT, MP, NQ, MEM( IPA-IPREPAD ),
     $                            DESCA( LLD_ ), IPREPAD, IPOSTPAD,
     $                            PADVAL )
                  CALL PSFILLPAD( ICTXT, LIPIV, 1, MEM( IPPIV-IPREPAD ),
     $                            LIPIV, IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  ANORM = PSLANGE( 'I', M, N, MEM( IPA ), 1, 1, DESCA,
     $                             MEM( IPW ) )
                  ANORM1 = PSLANGE( '1', M, N, MEM( IPA ), 1, 1, DESCA,
     $                             MEM( IPW ) )
                  CALL PSCHEKPAD( ICTXT, 'PSLANGE', MP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSLANGE', WORKSIZ-IPOSTPAD,
     $                            1, MEM( IPW-IPREPAD ),
     $                            WORKSIZ-IPOSTPAD, IPREPAD, IPOSTPAD,
     $                            PADVAL )
               END IF
*
               IF( EST .AND. M.EQ.N ) THEN
                  CALL PSMATGEN( ICTXT, 'No transpose', 'No transpose',
     $                           DESCA( M_ ), DESCA( N_ ), DESCA( MB_ ),
     $                           DESCA( NB_ ), MEM( IPA0 ),
     $                           DESCA( LLD_ ), DESCA( RSRC_ ),
     $                           DESCA( CSRC_ ), IASEED, 0, MP, 0, NQ,
     $                           MYROW, MYCOL, NPROW, NPCOL )
                  IF( CHECK )
     $               CALL PSFILLPAD( ICTXT, MP, NQ, MEM( IPA0-IPREPAD ),
     $                               DESCA( LLD_ ), IPREPAD, IPOSTPAD,
     $                               PADVAL )
               END IF
*
               CALL SLBOOT()
               CALL BLACS_BARRIER( ICTXT, 'All' )
               CALL SLTIMER( 1 )
*
*              Perform LU factorization
*
               CALL PSGETRF( M, N, MEM( IPA ), 1, 1, DESCA,
     $                       MEM( IPPIV ), INFO )
*
               CALL SLTIMER( 1 )
*
               IF( INFO.NE.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = * ) 'PSGETRF INFO=', INFO
                  KFAIL = KFAIL + 1
                  RCOND = ZERO
                  GO TO 30
               END IF
*
               IF( CHECK ) THEN
*
*                 Check for memory overwrite in LU factorization
*
                  CALL PSCHEKPAD( ICTXT, 'PSGETRF', MP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSGETRF', LIPIV, 1,
     $                            MEM( IPPIV-IPREPAD ), LIPIV, IPREPAD,
     $                            IPOSTPAD, PADVAL )
               END IF
*
               IF( M.NE.N ) THEN
*
*                 For non-square matrices, factorization only
*
                  NRHS = 0
                  NBRHS = 0
*
                  IF( CHECK ) THEN
*
*                    Compute FRESID = ||A - P*L*U|| / (||A|| * N * eps)
*
                     CALL PSGETRRV( M, N, MEM( IPA ), 1, 1, DESCA,
     $                              MEM( IPPIV ), MEM( IPW ) )
                     CALL PSLAFCHK( 'No', 'No', M, N, MEM( IPA ), 1, 1,
     $                              DESCA, IASEED, ANORM, FRESID,
     $                              MEM( IPW ) )
*
*                    Check for memory overwrite
*
                     CALL PSCHEKPAD( ICTXT, 'PSGETRRV', MP, NQ,
     $                               MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     CALL PSCHEKPAD( ICTXT, 'PSGETRRV', LIPIV, 1,
     $                               MEM( IPPIV-IPREPAD ), LIPIV,
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     CALL PSCHEKPAD( ICTXT, 'PSGETRRV',
     $                               WORKSIZ-IPOSTPAD, 1,
     $                               MEM( IPW-IPREPAD ),
     $                               WORKSIZ-IPOSTPAD, IPREPAD,
     $                               IPOSTPAD, PADVAL )
*
*                    Test residual and detect NaN result
*
                     IF( ( FRESID.LE.THRESH          ) .AND.
     $                   ( (FRESID-FRESID).EQ.0.0E+0 ) ) THEN
                        KPASS = KPASS + 1
                        PASSED = 'PASSED'
                     ELSE
                        KFAIL = KFAIL + 1
                        PASSED = 'FAILED'
                        IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $                     WRITE( NOUT, FMT = 9986 ) FRESID
                     END IF
*
                  ELSE
*
*                    Don't perform the checking, only timing
*
                     KPASS = KPASS + 1
                     FRESID = FRESID - FRESID
                     PASSED = 'BYPASS'
*
                  END IF
*
*                 Gather maximum of all CPU and WALL clock timings
*
                  CALL SLCOMBINE( ICTXT, 'All', '>', 'W', 1, 1,
     $                            WTIME )
                  CALL SLCOMBINE( ICTXT, 'All', '>', 'C', 1, 1,
     $                            CTIME )
*
*                 Print results
*
                  IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
*
                     MAXMN = MAX( M, N )
                     MINMN = MIN( M, N )
*
*                    M N^2 - 1/3 N^3 - 1/2 N^2 flops for LU
*                    factorization when M >= N
*
                     NOPS = DBLE( MAXMN )*( DBLE( MINMN )**2 ) -
     $                      (1.0D+0 / 3.0D+0)*( DBLE( MINMN )**3 ) -
     $                      (1.0D+0 / 2.0D+0)*( DBLE( MINMN )**2 )
*
*                    Calculate total megaflops -- factorization only,
*                    -- for WALL and CPU time, and print output
*
*                    Print WALL time if machine supports it
*
                     IF( WTIME( 1 ).GT.0.0D+0 ) THEN
                        TMFLOPS = NOPS / ( WTIME( 1 ) * 1.0D+6 )
                     ELSE
                        TMFLOPS = 0.0D+0
                     END IF
*
                     WTIME( 2 ) = 0.0D+0
                     IF( WTIME( 1 ).GE.0.0D+0 )
     $                  WRITE( NOUT, FMT = 9993 ) 'WALL', M, N, NB,
     $                         NRHS, NBRHS, NPROW, NPCOL, WTIME( 1 ),
     $                         WTIME( 2 ), TMFLOPS, PASSED
*
*                    Print CPU time if machine supports it
*
                     IF( CTIME( 1 ).GT.0.0D+0 ) THEN
                        TMFLOPS = NOPS / ( CTIME( 1 ) * 1.0D+6 )
                     ELSE
                        TMFLOPS = 0.0D+0
                     END IF
*
                     CTIME( 2 ) = 0.0D+0
                     IF( CTIME( 1 ).GE.0.0D+0 )
     $                  WRITE( NOUT, FMT = 9993 ) 'CPU ', M, N, NB,
     $                         NRHS, NBRHS, NPROW, NPCOL, CTIME( 1 ),
     $                         CTIME( 2 ), TMFLOPS, PASSED
                  END IF
*
               ELSE
*
*                 If M = N
*
                  IF( EST ) THEN
*
*                    Calculate workspace required for PSGECON
*
                     LWORK = MAX( 1, 2*NP ) + MAX( 1, 2*NQ ) +
     $                       MAX( 2, DESCA( NB_ )*
     $                       MAX( 1, ICEIL( NPROW-1, NPCOL ) ),
     $                       NQ + DESCA( NB_ )*
     $                       MAX( 1, ICEIL( NPCOL-1, NPROW ) ) )
                     IPW2  = IPW + LWORK + IPOSTPAD + IPREPAD
                     LIWORK = MAX( 1, NP )
                     LW2 = ICEIL( LIWORK*INTGSZ, REALSZ ) + IPOSTPAD
*
                     IERR( 1 ) = 0
                     IF( IPW2+LW2.GT.MEMSIZ ) THEN
                        IF( IAM.EQ.0 )
     $                     WRITE( NOUT, FMT = 9996 )'cond est',
     $                     ( IPW2+LW2 )*REALSZ
                        IERR( 1 ) = 1
                     END IF
*
*                    Check all processes for an error
*
                     CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                             -1, 0 )
*
                     IF( IERR( 1 ).GT.0 ) THEN
                        IF( IAM.EQ.0 )
     $                     WRITE( NOUT, FMT = 9997 ) 'MEMORY'
                        KSKIP = KSKIP + 1
                        GO TO 30
                     END IF
*
                     IF( CHECK ) THEN
                        CALL PSFILLPAD( ICTXT, LWORK, 1,
     $                                  MEM( IPW-IPREPAD ), LWORK,
     $                                  IPREPAD, IPOSTPAD, PADVAL )
                        CALL PSFILLPAD( ICTXT, LW2-IPOSTPAD, 1,
     $                                  MEM( IPW2-IPREPAD ),
     $                                  LW2-IPOSTPAD, IPREPAD,
     $                                  IPOSTPAD, PADVAL )
                     END IF
*
*                    Compute condition number of the matrix
*
                     CALL PSGECON( '1', N, MEM( IPA ), 1, 1, DESCA,
     $                             ANORM1, RCOND, MEM( IPW ), LWORK,
     $                             MEM( IPW2 ), LIWORK, INFO )
*
                     IF( CHECK ) THEN
                        CALL PSCHEKPAD( ICTXT, 'PSGECON', NP, NQ,
     $                                  MEM( IPA-IPREPAD ),
     $                                  DESCA( LLD_ ), IPREPAD,
     $                                  IPOSTPAD, PADVAL )
                        CALL PSCHEKPAD( ICTXT, 'PSGECON', LWORK, 1,
     $                                  MEM( IPW-IPREPAD ), LWORK,
     $                                  IPREPAD, IPOSTPAD, PADVAL )
                        CALL PSCHEKPAD( ICTXT, 'PSGECON',
     $                                  LW2-IPOSTPAD, 1,
     $                                  MEM( IPW2-IPREPAD ),
     $                                  LW2-IPOSTPAD, IPREPAD,
     $                                  IPOSTPAD, PADVAL )
                     END IF
                  END IF
*
*                 Loop over the different values for NRHS
*
                  DO 20 HH = 1, NNR
*
                     NRHS = NRVAL( HH )
*
                     DO 10 KK = 1, NNBR
*
                        NBRHS = NBRVAL( KK )
*
*                       Initialize Array Descriptor for rhs
*
                        CALL DESCINIT( DESCB, N, NRHS, NB, NBRHS, 0, 0,
     $                                 ICTXT, MAX( 1, NP )+IMIDPAD,
     $                                 IERR( 1 ) )
*
*                       Check all processes for an error
*
                        CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                                -1, 0 )
*
                        IF( IERR( 1 ).LT.0 ) THEN
                           IF( IAM.EQ.0 )
     $                        WRITE( NOUT, FMT = 9997 ) 'descriptor'
                           KSKIP = KSKIP + 1
                           GO TO 10
                        END IF
*
*                       move IPW to allow room for RHS
*
                        MYRHS = NUMROC( DESCB( N_ ), DESCB( NB_ ),
     $                                  MYCOL, DESCB( CSRC_ ), NPCOL )
                        IPB = IPW
*
                        IF( EST ) THEN
                           IPB0 = IPB + DESCB( LLD_ )*MYRHS + IPOSTPAD +
     $                            IPREPAD
                           IPFERR = IPB0 + DESCB( LLD_ )*MYRHS +
     $                              IPOSTPAD + IPREPAD
                           IPBERR = MYRHS + IPFERR + IPOSTPAD + IPREPAD
                           IPW = MYRHS + IPBERR + IPOSTPAD + IPREPAD
                        ELSE
                           IPW = IPB + DESCB( LLD_ )*MYRHS + IPOSTPAD +
     $                           IPREPAD
                        END IF
*
*                       Set worksiz: routines requiring most workspace
*                       is PSLASCHK
*
                        IF( CHECK ) THEN
                           LCM = ILCM( NPROW, NPCOL )
                           LCMQ = LCM / NPCOL
                           WORKSIZ = MAX( WORKSIZ-IPOSTPAD,
     $                       NQ * NBRHS + NP * NBRHS +
     $                       MAX( MAX( NQ*NB, 2*NBRHS ),
     $                       NBRHS * NUMROC( NUMROC(N,NB,0,0,NPCOL),NB,
     $                       0,0,LCMQ ) ) )
                           WORKSIZ = IPOSTPAD + WORKSIZ
                        ELSE
                           WORKSIZ = IPOSTPAD
                        END IF
*
                        IERR( 1 ) = 0
                        IF( IPW+WORKSIZ.GT.MEMSIZ ) THEN
                           IF( IAM.EQ.0 )
     $                        WRITE( NOUT, FMT = 9996 )'solve',
     $                        ( IPW+WORKSIZ )*REALSZ
                           IERR( 1 ) = 1
                        END IF
*
*                       Check all processes for an error
*
                        CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                                -1, 0 )
*
                        IF( IERR( 1 ).GT.0 ) THEN
                           IF( IAM.EQ.0 )
     $                        WRITE( NOUT, FMT = 9997 ) 'MEMORY'
                           KSKIP = KSKIP + 1
                           GO TO 10
                        END IF
*
*                       Generate RHS
*
                        CALL PSMATGEN( ICTXT, 'No', 'No', DESCB( M_ ),
     $                                 DESCB( N_ ), DESCB( MB_ ),
     $                                 DESCB( NB_ ), MEM( IPB ),
     $                                 DESCB( LLD_ ), DESCB( RSRC_ ),
     $                                 DESCB( CSRC_ ), IBSEED, 0, NP, 0,
     $                                 MYRHS, MYROW, MYCOL, NPROW,
     $                                 NPCOL )
*
                        IF( CHECK )
     $                     CALL PSFILLPAD( ICTXT, NP, MYRHS,
     $                                     MEM( IPB-IPREPAD ),
     $                                     DESCB( LLD_ ), IPREPAD,
     $                                     IPOSTPAD, PADVAL )
*
                        IF( EST ) THEN
                           CALL PSMATGEN( ICTXT, 'No', 'No',
     $                                    DESCB( M_ ), DESCB( N_ ),
     $                                    DESCB( MB_ ), DESCB( NB_ ),
     $                                    MEM( IPB0 ), DESCB( LLD_ ),
     $                                    DESCB( RSRC_ ),
     $                                    DESCB( CSRC_ ), IBSEED, 0, NP,
     $                                    0, MYRHS, MYROW, MYCOL, NPROW,
     $                                    NPCOL )
                           IF( CHECK ) THEN
                              CALL PSFILLPAD( ICTXT, NP, MYRHS,
     $                                        MEM( IPB0-IPREPAD ),
     $                                        DESCB( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSFILLPAD( ICTXT, 1, MYRHS,
     $                                        MEM( IPFERR-IPREPAD ), 1,
     $                                        IPREPAD, IPOSTPAD,
     $                                        PADVAL )
                              CALL PSFILLPAD( ICTXT, 1, MYRHS,
     $                                        MEM( IPBERR-IPREPAD ), 1,
     $                                        IPREPAD, IPOSTPAD,
     $                                        PADVAL )
                           END IF
                        END IF
*
                        CALL BLACS_BARRIER( ICTXT, 'All' )
                        CALL SLTIMER( 2 )
*
*                       Solve linear sytem via LU factorization
*
                        CALL PSGETRS( 'No', N, NRHS, MEM( IPA ), 1, 1,
     $                                DESCA, MEM( IPPIV ), MEM( IPB ),
     $                                1, 1, DESCB, INFO )
*
                        CALL SLTIMER( 2 )
*
                        IF( CHECK ) THEN
*
*                          check for memory overwrite
*
                           CALL PSCHEKPAD( ICTXT, 'PSGETRS', NP, NQ,
     $                                     MEM( IPA-IPREPAD ),
     $                                     DESCA( LLD_ ), IPREPAD,
     $                                     IPOSTPAD, PADVAL )
                           CALL PSCHEKPAD( ICTXT, 'PSGETRS', LIPIV, 1,
     $                                     MEM( IPPIV-IPREPAD ), LIPIV,
     $                                     IPREPAD, IPOSTPAD, PADVAL )
                           CALL PSCHEKPAD( ICTXT, 'PSGETRS', NP,
     $                                     MYRHS, MEM( IPB-IPREPAD ),
     $                                     DESCB( LLD_ ), IPREPAD,
     $                                     IPOSTPAD, PADVAL )
*
                           CALL PSFILLPAD( ICTXT, WORKSIZ-IPOSTPAD,
     $                                     1, MEM( IPW-IPREPAD ),
     $                                     WORKSIZ-IPOSTPAD, IPREPAD,
     $                                     IPOSTPAD, PADVAL )
*
*                          check the solution to rhs
*
                           CALL PSLASCHK( 'No', 'N', N, NRHS,
     $                                    MEM( IPB ), 1, 1, DESCB,
     $                                    IASEED, 1, 1, DESCA, IBSEED,
     $                                    ANORM, SRESID, MEM( IPW ) )
*
                           IF( IAM.EQ.0 .AND. SRESID.GT.THRESH )
     $                        WRITE( NOUT, FMT = 9985 ) SRESID
*
*                          check for memory overwrite
*
                           CALL PSCHEKPAD( ICTXT, 'PSLASCHK', NP,
     $                                     MYRHS, MEM( IPB-IPREPAD ),
     $                                     DESCB( LLD_ ), IPREPAD,
     $                                     IPOSTPAD, PADVAL )
                           CALL PSCHEKPAD( ICTXT, 'PSLASCHK',
     $                                     WORKSIZ-IPOSTPAD, 1,
     $                                     MEM( IPW-IPREPAD ),
     $                                     WORKSIZ-IPOSTPAD,
     $                                     IPREPAD, IPOSTPAD, PADVAL )
*
*                          The second test is a NaN trap
*
                           IF( SRESID.LE.THRESH .AND.
     $                         ( SRESID-SRESID ).EQ.0.0E+0 ) THEN
                              KPASS = KPASS + 1
                              PASSED = 'PASSED'
                           ELSE
                              KFAIL = KFAIL + 1
                              PASSED = 'FAILED'
                           END IF
                        ELSE
                           KPASS = KPASS + 1
                           SRESID = SRESID - SRESID
                           PASSED = 'BYPASS'
                        END IF
*
                        IF( EST ) THEN
*
*                          Calculate workspace required for PSGERFS
*
                           LWORK = MAX( 1, 3*NP )
                           IPW2  = IPW + LWORK + IPOSTPAD + IPREPAD
                           LIWORK = MAX( 1, NP )
                           LW2 = ICEIL( LIWORK*INTGSZ, REALSZ ) +
     $                           IPOSTPAD
*
                           IERR( 1 ) = 0
                           IF( IPW2+LW2.GT.MEMSIZ ) THEN
                              IF( IAM.EQ.0 )
     $                           WRITE( NOUT, FMT = 9996 )
     $                           'iter ref', ( IPW2+LW2 )*REALSZ
                              IERR( 1 ) = 1
                           END IF
*
*                          Check all processes for an error
*
                           CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1,
     $                                   IERR, 1, -1, 0 )
*
                           IF( IERR( 1 ).GT.0 ) THEN
                              IF( IAM.EQ.0 )
     $                           WRITE( NOUT, FMT = 9997 )
     $                           'MEMORY'
                              KSKIP = KSKIP + 1
                              GO TO 10
                           END IF
*
                           IF( CHECK ) THEN
                              CALL PSFILLPAD( ICTXT, LWORK, 1,
     $                                        MEM( IPW-IPREPAD ),
     $                                        LWORK, IPREPAD, IPOSTPAD,
     $                                        PADVAL )
                              CALL PSFILLPAD( ICTXT, LW2-IPOSTPAD, 1,
     $                                        MEM( IPW2-IPREPAD ),
     $                                        LW2-IPOSTPAD, IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                           END IF
*
*                          Use iterative refinement to improve the
*                          computed solution
*
                           CALL PSGERFS( 'No', N, NRHS, MEM( IPA0 ), 1,
     $                                   1, DESCA, MEM( IPA ), 1, 1,
     $                                   DESCA, MEM( IPPIV ),
     $                                   MEM( IPB0 ), 1, 1, DESCB,
     $                                   MEM( IPB ), 1, 1, DESCB,
     $                                   MEM( IPFERR ), MEM( IPBERR ),
     $                                   MEM( IPW ), LWORK, MEM( IPW2 ),
     $                                   LIWORK, INFO )
*
                           IF( CHECK ) THEN
                              CALL PSCHEKPAD( ICTXT, 'PSGERFS', NP,
     $                                        NQ, MEM( IPA0-IPREPAD ),
     $                                        DESCA( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSGERFS', NP,
     $                                        NQ, MEM( IPA-IPREPAD ),
     $                                        DESCA( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSGERFS', LIPIV,
     $                                        1, MEM( IPPIV-IPREPAD ),
     $                                        LIPIV, IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSGERFS', NP,
     $                                        MYRHS, MEM( IPB-IPREPAD ),
     $                                        DESCB( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSGERFS', NP,
     $                                        MYRHS,
     $                                        MEM( IPB0-IPREPAD ),
     $                                        DESCB( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSGERFS', 1,
     $                                        MYRHS,
     $                                        MEM( IPFERR-IPREPAD ), 1,
     $                                        IPREPAD, IPOSTPAD,
     $                                        PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSGERFS', 1,
     $                                        MYRHS,
     $                                        MEM( IPBERR-IPREPAD ), 1,
     $                                        IPREPAD, IPOSTPAD,
     $                                        PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSGERFS', LWORK,
     $                                        1, MEM( IPW-IPREPAD ),
     $                                        LWORK, IPREPAD, IPOSTPAD,
     $                                        PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSGERFS',
     $                                        LW2-IPOSTPAD, 1,
     $                                        MEM( IPW2-IPREPAD ),
     $                                        LW2-IPOSTPAD, IPREPAD,
     $                                        IPOSTPAD, PADVAL )
*
                              CALL PSFILLPAD( ICTXT, WORKSIZ-IPOSTPAD,
     $                                        1, MEM( IPW-IPREPAD ),
     $                                        WORKSIZ-IPOSTPAD, IPREPAD,
     $                                        IPOSTPAD, PADVAL )
*
*                             check the solution to rhs
*
                              CALL PSLASCHK( 'No', 'N', N, NRHS,
     $                                       MEM( IPB ), 1, 1, DESCB,
     $                                       IASEED, 1, 1, DESCA,
     $                                       IBSEED, ANORM, SRESID2,
     $                                       MEM( IPW ) )
*
                              IF( IAM.EQ.0 .AND. SRESID2.GT.THRESH )
     $                           WRITE( NOUT, FMT = 9985 ) SRESID2
*
*                             check for memory overwrite
*
                              CALL PSCHEKPAD( ICTXT, 'PSLASCHK', NP,
     $                                        MYRHS, MEM( IPB-IPREPAD ),
     $                                        DESCB( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSLASCHK',
     $                                        WORKSIZ-IPOSTPAD, 1,
     $                                        MEM( IPW-IPREPAD ),
     $                                        WORKSIZ-IPOSTPAD, IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                           END IF
                        END IF
*
*                       Gather max. of all CPU and WALL clock timings
*
                        CALL SLCOMBINE( ICTXT, 'All', '>', 'W', 2, 1,
     $                                  WTIME )
                        CALL SLCOMBINE( ICTXT, 'All', '>', 'C', 2, 1,
     $                                  CTIME )
*
*                       Print results
*
                        IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
*
*                          2/3 N^3 - 1/2 N^2 flops for LU factorization
*
                           NOPS = (2.0D+0/3.0D+0)*( DBLE(N)**3 ) -
     $                            (1.0D+0/2.0D+0)*( DBLE(N)**2 )
*
*                          nrhs * 2 N^2 flops for LU solve.
*
                           NOPS = NOPS + 2.0D+0*(DBLE(N)**2)*DBLE(NRHS)
*
*                          Calculate total megaflops -- factorization
*                          and solve -- for WALL and CPU time, and print
*                          output
*
*                          Print WALL time if machine supports it
*
                           IF( WTIME( 1 ) + WTIME( 2 ) .GT. 0.0D+0 )
     $                        THEN
                              TMFLOPS = NOPS /
     $                            ( ( WTIME( 1 )+WTIME( 2 ) ) * 1.0D+6 )
                           ELSE
                              TMFLOPS = 0.0D+0
                           END IF
*
*                          Print WALL time if supported
*
                           IF( WTIME( 2 ).GE.0.0D+0 )
     $                        WRITE( NOUT, FMT = 9993 ) 'WALL', M, N,
     $                               NB, NRHS, NBRHS, NPROW, NPCOL,
     $                               WTIME( 1 ), WTIME( 2 ), TMFLOPS,
     $                               PASSED
*
*                          Print CPU time if supported
*
                           IF( CTIME( 1 )+CTIME( 2 ).GT.0.0D+0 )
     $                        THEN
                              TMFLOPS = NOPS /
     $                            ( ( CTIME( 1 )+CTIME( 2 ) ) * 1.0D+6 )
                           ELSE
                              TMFLOPS = 0.0D+0
                           END IF
*
                           IF( CTIME( 2 ).GE.0.0D+0 )
     $                        WRITE( NOUT, FMT = 9993 ) 'CPU ', M, N,
     $                               NB, NRHS, NBRHS, NPROW, NPCOL,
     $                               CTIME( 1 ), CTIME( 2 ), TMFLOPS,
     $                               PASSED
                        END IF
   10                CONTINUE
   20             CONTINUE
*
                  IF( CHECK.AND.( SRESID.GT.THRESH ) ) THEN
*
*                    Compute fresid = ||A - P*L*U|| / (||A|| * N * eps)
*
                     CALL PSGETRRV( M, N, MEM( IPA ), 1, 1, DESCA,
     $                              MEM( IPPIV ), MEM( IPW ) )
                     CALL PSLAFCHK( 'No', 'No', M, N, MEM( IPA ), 1,
     $                              1, DESCA, IASEED, ANORM, FRESID,
     $                              MEM( IPW ) )
*
*                    Check for memory overwrite
*
                     CALL PSCHEKPAD( ICTXT, 'PSGETRRV', NP, NQ,
     $                               MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     CALL PSCHEKPAD( ICTXT, 'PSGETRRV', LIPIV,
     $                               1, MEM( IPPIV-IPREPAD ), LIPIV,
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     CALL PSCHEKPAD( ICTXT, 'PSGETRRV',
     $                               WORKSIZ-IPOSTPAD, 1,
     $                               MEM( IPW-IPREPAD ),
     $                               WORKSIZ-IPOSTPAD, IPREPAD,
     $                               IPOSTPAD, PADVAL )
*
                     IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $                  WRITE( NOUT, FMT = 9986 ) FRESID
                  END IF
               END IF
   30       CONTINUE
   40    CONTINUE
         CALL BLACS_GRIDEXIT( ICTXT )
   50 CONTINUE
*
*     Print ending messages and close output file
*
   60 CONTINUE
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
     $        '; It should be at least 1' )
 9998 FORMAT( 'ILLEGAL GRID: nprow*npcol = ', I4, '. It can be at most',
     $        I4 )
 9997 FORMAT( 'Bad ', A6, ' parameters: going on to next test case.' )
 9996 FORMAT( 'Unable to perform ', A, ': need TOTMEM of at least',
     $        I11 )
 9995 FORMAT( 'TIME     M     N  NB NRHS NBRHS    P    Q  LU Time ',
     $        'Sol Time  MFLOPS  CHECK' )
 9994 FORMAT( '---- ----- ----- --- ---- ----- ---- ---- -------- ',
     $        '-------- -------- ------' )
 9993 FORMAT( A4, 1X, I5, 1X, I5, 1X, I3, 1X, I5, 1X, I4, 1X, I4, 1X,
     $        I4, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, A6 )
 9992 FORMAT( 'Finished ', I6, ' tests, with the following results:' )
 9991 FORMAT( I5, ' tests completed and passed residual checks.' )
 9990 FORMAT( I5, ' tests completed without checking.' )
 9989 FORMAT( I5, ' tests completed and failed residual checks.' )
 9988 FORMAT( I5, ' tests skipped because of illegal input values.' )
 9987 FORMAT( 'END OF TESTS.' )
 9986 FORMAT( '||A - P*L*U|| / (||A|| * N * eps) = ', G25.7 )
 9985 FORMAT( '||Ax-b||/(||x||*||A||*eps*N) ', F25.7 )
*
      STOP
*
*     End of PSLUDRIVER
*
      END
