      PROGRAM PSLLTDRIVER
*
*  -- ScaLAPACK testing driver (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*  Purpose
*  =======
*
*  PSLLTDRIVER is the main test program for the REAL
*  ScaLAPACK Cholesky routines.  This test driver performs an
*  A = L*L**T or A = U**T*U factorization and solve, and optionally
*  performs condition estimation and iterative refinement.
*
*  The program must be driven by a short data file.  An annotated
*  example of a data file can be obtained by deleting the first 3
*  characters from the following 18 lines:
*  'ScaLAPACK LLt factorization input file'
*  'Intel iPSC/860 hypercube, gamma model.'
*  'LLT.out'            output file name (if any)
*  6                    device out
*  'U'                  define Lower or Upper
*  1                    number of problems sizes
*  31 100 200           values of N
*  1                    number of NB's
*  2 10 24              values of NB
*  1                    number of NRHS's
*  1                    values of NRHS
*  1                    Number of NBRHS's
*  1                    values of NBRHS
*  1                    number of process grids (ordered pairs of P & Q)
*  2                    values of P
*  2                    values of Q
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
      CHARACTER          UPLO
      CHARACTER*6        PASSED
      CHARACTER*80       OUTFILE
      INTEGER            HH, I, IAM, IASEED, IBSEED, ICTXT, IMIDPAD,
     $                   INFO, IPA, IPA0, IPB, IPB0, IPBERR, IPFERR,
     $                   IPREPAD, IPOSTPAD, IPW, IPW2, ITEMP, J, K,
     $                   KFAIL, KK, KPASS, KSKIP, KTESTS, LCM, LCMQ,
     $                   LIWORK, LWORK, LW2, MYCOL, MYRHS, MYROW, N, NB,
     $                   NBRHS, NGRIDS, NMAT, NNB, NNBR, NNR, NOUT, NP,
     $                   NPCOL, NPROCS, NPROW, NQ, NRHS, WORKSIZ
      REAL               ANORM, ANORM1, FRESID, RCOND, SRESID, SRESID2,
     $                   THRESH
      DOUBLE PRECISION   NOPS, TMFLOPS
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), IERR( 1 ),
     $                   NBRVAL( NTESTS ), NBVAL( NTESTS ),
     $                   NRVAL( NTESTS ), NVAL( NTESTS ),
     $                   PVAL( NTESTS ), QVAL( NTESTS )
      REAL               MEM( MEMSIZ )
      DOUBLE PRECISION   CTIME( 2 ), WTIME( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_EXIT, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, DESCINIT,
     $                   IGSUM2D, BLACS_PINFO, PSCHEKPAD, PSFILLPAD,
     $                   PSLAFCHK, PSLASCHK, PSLLTINFO,
     $                   PSMATGEN, PSPOCON, PSPORFS,
     $                   PSPOTRF, PSPOTRRV, PSPOTRS, SLBOOT,
     $                   SLCOMBINE, SLTIMER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, ILCM, NUMROC
      REAL               PSLANSY
      EXTERNAL           ICEIL, ILCM, LSAME, NUMROC, PSLANSY
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
      CALL PSLLTINFO( OUTFILE, NOUT, UPLO, NMAT, NVAL, NTESTS, NNB,
     $                NBVAL, NTESTS, NNR, NRVAL, NTESTS, NNBR, NBRVAL,
     $                NTESTS, NGRIDS, PVAL, NTESTS, QVAL, NTESTS,
     $                THRESH, EST, MEM, IAM, NPROCS )
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
            N = NVAL( J )
*
*           Make sure matrix information is correct
*
            IERR( 1 ) = 0
            IF( N.LT.1 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9999 ) 'MATRIX', 'N', N
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
     $                        MAX( 1, NP )+IMIDPAD, IERR( 1 ) )
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
               IF( EST ) THEN
                  IPA0 = IPA + DESCA( LLD_ )*NQ + IPOSTPAD + IPREPAD
                  IPW = IPA0 + DESCA( LLD_ )*NQ + IPOSTPAD + IPREPAD
               ELSE
                  IPW = IPA + DESCA( LLD_ )*NQ + IPOSTPAD + IPREPAD
               END IF
*
*
               IF( CHECK ) THEN
*
*                 Calculate the amount of workspace required by
*                 the checking routines PSLAFCHK, PSPOTRRV, and
*                 PSLANSY
*
                  WORKSIZ = NP * DESCA( NB_ )
*
                  WORKSIZ = MAX( WORKSIZ, DESCA( MB_ ) * DESCA( NB_ ) )
*
                  LCM = ILCM( NPROW, NPCOL )
                  ITEMP = MAX( 2, 2 * NQ ) + NP
                  IF( NPROW.NE.NPCOL ) THEN
                     ITEMP = ITEMP +
     $                       NB * ICEIL( ICEIL( NP, NB ), LCM / NPROW )
                  END IF
                  WORKSIZ = MAX( WORKSIZ, ITEMP )
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
     $               WRITE( NOUT, FMT = 9997 ) 'MEMORY'
                  KSKIP = KSKIP + 1
                  GO TO 30
               END IF
*
*              Generate a symmetric positive definite matrix A
*
               CALL PSMATGEN( ICTXT, 'Symm', 'Diag', DESCA( M_ ),
     $                        DESCA( N_ ), DESCA( MB_ ), DESCA( NB_ ),
     $                        MEM( IPA ), DESCA( LLD_ ), DESCA( RSRC_ ),
     $                        DESCA( CSRC_ ), IASEED, 0, NP, 0, NQ,
     $                        MYROW, MYCOL, NPROW, NPCOL )
*
*              Calculate inf-norm of A for residual error-checking
*
               IF( CHECK ) THEN
                  CALL PSFILLPAD( ICTXT, NP, NQ, MEM( IPA-IPREPAD ),
     $                             DESCA( LLD_ ), IPREPAD, IPOSTPAD,
     $                             PADVAL )
                  CALL PSFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                             MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                             IPREPAD, IPOSTPAD, PADVAL )
                  ANORM = PSLANSY( 'I', UPLO, N, MEM( IPA ), 1, 1,
     $                             DESCA, MEM( IPW ) )
                  ANORM1 = PSLANSY( '1', UPLO, N, MEM( IPA ), 1, 1,
     $                             DESCA, MEM( IPW ) )
                  CALL PSCHEKPAD( ICTXT, 'PSLANSY', NP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSLANSY',
     $                            WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
               END IF
*
               IF( EST ) THEN
                  CALL PSMATGEN( ICTXT, 'Symm', 'Diag', DESCA( M_ ),
     $                           DESCA( N_ ), DESCA( MB_ ),
     $                           DESCA( NB_ ), MEM( IPA0 ),
     $                           DESCA( LLD_ ), DESCA( RSRC_ ),
     $                           DESCA( CSRC_ ), IASEED, 0, NP, 0, NQ,
     $                           MYROW, MYCOL, NPROW, NPCOL )
                  IF( CHECK )
     $               CALL PSFILLPAD( ICTXT, NP, NQ,
     $                               MEM( IPA0-IPREPAD ), DESCA( LLD_ ),
     $                               IPREPAD, IPOSTPAD, PADVAL )
               END IF
*
               CALL SLBOOT()
               CALL BLACS_BARRIER( ICTXT, 'All' )
*
*              Perform LLt factorization
*
               CALL SLTIMER( 1 )
*
               CALL PSPOTRF( UPLO, N, MEM( IPA ), 1, 1, DESCA, INFO )
*
               CALL SLTIMER( 1 )
*
               IF( INFO.NE.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = * ) 'PSPOTRF INFO=', INFO
                  KFAIL = KFAIL + 1
                  RCOND = ZERO
                  GO TO 60
               END IF
*
               IF( CHECK ) THEN
*
*                 Check for memory overwrite in LLt factorization
*
                  CALL PSCHEKPAD( ICTXT, 'PSPOTRF', NP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
               END IF
*
               IF( EST ) THEN
*
*                 Calculate workspace required for PSPOCON
*
                  LWORK = MAX( 1, 2*NP ) + MAX( 1, 2*NQ ) +
     $                    MAX( 2, DESCA( NB_ )*
     $                    MAX( 1, ICEIL( NPROW-1, NPCOL ) ),
     $                    NQ + DESCA( NB_ )*
     $                    MAX( 1, ICEIL( NPCOL-1, NPROW ) ) )
                  IPW2  = IPW + LWORK + IPOSTPAD + IPREPAD
                  LIWORK = MAX( 1, NP )
                  LW2 = ICEIL( LIWORK*INTGSZ, REALSZ ) + IPOSTPAD
*
                  IERR( 1 ) = 0
                  IF( IPW2+LW2.GT.MEMSIZ ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9996 )'cond est',
     $                  ( IPW2+LW2 )*REALSZ
                     IERR( 1 ) = 1
                  END IF
*
*                 Check all processes for an error
*
                  CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                          -1, 0 )
*
                  IF( IERR( 1 ).GT.0 ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9997 ) 'MEMORY'
                     KSKIP = KSKIP + 1
                     GO TO 60
                  END IF
*
                  IF( CHECK ) THEN
                     CALL PSFILLPAD( ICTXT, LWORK, 1,
     $                               MEM( IPW-IPREPAD ), LWORK,
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     CALL PSFILLPAD( ICTXT, LW2-IPOSTPAD, 1,
     $                               MEM( IPW2-IPREPAD ),
     $                               LW2-IPOSTPAD, IPREPAD,
     $                               IPOSTPAD, PADVAL )
                  END IF
*
*                 Compute condition number of the matrix
*
                  CALL PSPOCON( UPLO, N, MEM( IPA ), 1, 1, DESCA,
     $                          ANORM1, RCOND, MEM( IPW ), LWORK,
     $                          MEM( IPW2 ), LIWORK, INFO )
*
                  IF( CHECK ) THEN
                     CALL PSCHEKPAD( ICTXT, 'PSPOCON', NP, NQ,
     $                               MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     CALL PSCHEKPAD( ICTXT, 'PSPOCON',
     $                               LWORK, 1, MEM( IPW-IPREPAD ),
     $                               LWORK, IPREPAD, IPOSTPAD,
     $                               PADVAL )
                     CALL PSCHEKPAD( ICTXT, 'PSPOCON',
     $                               LW2-IPOSTPAD, 1,
     $                               MEM( IPW2-IPREPAD ), LW2-IPOSTPAD,
     $                               IPREPAD, IPOSTPAD, PADVAL )
                  END IF
               END IF
*
*              Loop over the different values for NRHS
*
               DO 20 HH = 1, NNR
*
                  NRHS = NRVAL( HH )
*
                  DO 10 KK = 1, NNBR
*
                     NBRHS = NBRVAL( KK )
*
*                    Initialize Array Descriptor for rhs
*
                     CALL DESCINIT( DESCB, N, NRHS, NB, NBRHS, 0, 0,
     $                              ICTXT, MAX( 1, NP )+IMIDPAD,
     $                              IERR( 1 ) )
*
*                    move IPW to allow room for RHS
*
                     MYRHS = NUMROC( DESCB( N_ ), DESCB( NB_ ), MYCOL,
     $                               DESCB( CSRC_ ), NPCOL )
                     IPB = IPW
*
                     IF( EST ) THEN
                        IPB0 = IPB +  DESCB( LLD_ )*MYRHS + IPOSTPAD +
     $                           IPREPAD
                        IPFERR = IPB0 +  DESCB( LLD_ )*MYRHS + IPOSTPAD
     $                           + IPREPAD
                        IPBERR = MYRHS + IPFERR + IPOSTPAD + IPREPAD
                        IPW = MYRHS + IPBERR + IPOSTPAD + IPREPAD
                     ELSE
                        IPW = IPB +  DESCB( LLD_ )*MYRHS + IPOSTPAD +
     $                        IPREPAD
                     END IF
*
                     IF( CHECK ) THEN
*
*                       Calculate the amount of workspace required by
*                       the checking routines PSLASCHK
*
                        LCMQ = LCM / NPCOL
                        WORKSIZ = MAX( WORKSIZ-IPOSTPAD,
     $                    NQ * NBRHS + NP * NBRHS +
     $                    MAX( MAX( NQ*NB, 2*NBRHS ),
     $                    NBRHS * NUMROC( NUMROC(N,NB,0,0,NPCOL),NB,
     $                    0,0,LCMQ ) ) )
                        WORKSIZ = IPOSTPAD + WORKSIZ
                     ELSE
                        WORKSIZ = IPOSTPAD
                     END IF
*
                     IERR( 1 ) = 0
                     IF( IPW+WORKSIZ.GT.MEMSIZ ) THEN
                        IF( IAM.EQ.0 )
     $                     WRITE( NOUT, FMT = 9996 )'solve',
     $                            ( IPW+WORKSIZ )*REALSZ
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
                        GO TO 10
                     END IF
*
*                    Generate RHS
*
                     CALL PSMATGEN( ICTXT, 'No', 'No', DESCB( M_ ),
     $                              DESCB( N_ ), DESCB( MB_ ),
     $                              DESCB( NB_ ), MEM( IPB ),
     $                              DESCB( LLD_ ), DESCB( RSRC_ ),
     $                              DESCB( CSRC_ ), IBSEED, 0, NP, 0,
     $                              MYRHS, MYROW, MYCOL, NPROW, NPCOL )
*
                     IF( CHECK )
     $                  CALL PSFILLPAD( ICTXT, NP, MYRHS,
     $                                  MEM( IPB-IPREPAD ),
     $                                  DESCB( LLD_ ),
     $                                  IPREPAD, IPOSTPAD, PADVAL )
*
                     IF( EST ) THEN
                        CALL PSMATGEN( ICTXT, 'No', 'No', DESCB( M_ ),
     $                                 DESCB( N_ ), DESCB( MB_ ),
     $                                 DESCB( NB_ ), MEM( IPB0 ),
     $                                 DESCB( LLD_ ), DESCB( RSRC_ ),
     $                                 DESCB( CSRC_ ), IBSEED, 0, NP, 0,
     $                                 MYRHS, MYROW, MYCOL, NPROW,
     $                                 NPCOL )
*
                        IF( CHECK ) THEN
                           CALL PSFILLPAD( ICTXT, NP, MYRHS,
     $                                     MEM( IPB0-IPREPAD ),
     $                                     DESCB( LLD_ ), IPREPAD,
     $                                     IPOSTPAD, PADVAL )
                           CALL PSFILLPAD( ICTXT, 1, MYRHS,
     $                                     MEM( IPFERR-IPREPAD ), 1,
     $                                     IPREPAD, IPOSTPAD,
     $                                     PADVAL )
                           CALL PSFILLPAD( ICTXT, 1, MYRHS,
     $                                     MEM( IPBERR-IPREPAD ), 1,
     $                                     IPREPAD, IPOSTPAD,
     $                                     PADVAL )
                        END IF
                     END IF
*
                     CALL BLACS_BARRIER( ICTXT, 'All' )
                     CALL SLTIMER( 2 )
*
*                    Solve linear system via Cholesky factorization
*
                     CALL PSPOTRS( UPLO, N, NRHS, MEM( IPA ), 1, 1,
     $                             DESCA, MEM( IPB ), 1, 1, DESCB,
     $                             INFO )
*
                     CALL SLTIMER( 2 )
*
                     IF( CHECK ) THEN
*
*                       check for memory overwrite
*
                        CALL PSCHEKPAD( ICTXT, 'PSPOTRS', NP, NQ,
     $                                  MEM( IPA-IPREPAD ),
     $                                  DESCA( LLD_ ),
     $                                  IPREPAD, IPOSTPAD, PADVAL )
                        CALL PSCHEKPAD( ICTXT, 'PSPOTRS', NP,
     $                                  MYRHS, MEM( IPB-IPREPAD ),
     $                                  DESCB( LLD_ ), IPREPAD,
     $                                  IPOSTPAD, PADVAL )
*
                        CALL PSFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                                  MEM( IPW-IPREPAD ),
     $                                  WORKSIZ-IPOSTPAD, IPREPAD,
     $                                  IPOSTPAD, PADVAL )
*
*                       check the solution to rhs
*
                        CALL PSLASCHK( 'Symm', 'Diag', N, NRHS,
     $                                 MEM( IPB ), 1, 1, DESCB,
     $                                 IASEED, 1, 1, DESCA, IBSEED,
     $                                 ANORM, SRESID, MEM( IPW ) )
*
                        IF( IAM.EQ.0 .AND. SRESID.GT.THRESH )
     $                        WRITE( NOUT, FMT = 9985 ) SRESID
*
*                       check for memory overwrite
*
                        CALL PSCHEKPAD( ICTXT, 'PSLASCHK', NP,
     $                                  MYRHS, MEM( IPB-IPREPAD ),
     $                                  DESCB( LLD_ ), IPREPAD,
     $                                  IPOSTPAD, PADVAL )
                        CALL PSCHEKPAD( ICTXT, 'PSLASCHK',
     $                                  WORKSIZ-IPOSTPAD, 1,
     $                                  MEM( IPW-IPREPAD ),
     $                                  WORKSIZ-IPOSTPAD, IPREPAD,
     $                                  IPOSTPAD, PADVAL )
*
*                       The second test is a NaN trap
*
                        IF( ( SRESID.LE.THRESH          ).AND.
     $                      ( (SRESID-SRESID).EQ.0.0E+0 ) ) THEN
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
*                       Calculate workspace required for PSPORFS
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
                           CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR,
     $                                   1, -1, 0 )
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
                              CALL PSFILLPAD( ICTXT, LW2-IPOSTPAD,
     $                                        1, MEM( IPW2-IPREPAD ),
     $                                        LW2-IPOSTPAD,
     $                                        IPREPAD, IPOSTPAD,
     $                                        PADVAL )
                           END IF
*
*                          Use iterative refinement to improve the
*                          computed solution
*
                           CALL PSPORFS( UPLO, N, NRHS, MEM( IPA0 ),
     $                                   1, 1, DESCA, MEM( IPA ), 1, 1,
     $                                   DESCA, MEM( IPB0 ), 1, 1,
     $                                   DESCB, MEM( IPB ), 1, 1, DESCB,
     $                                   MEM( IPFERR ), MEM( IPBERR ),
     $                                   MEM( IPW ), LWORK, MEM( IPW2 ),
     $                                   LIWORK, INFO )
*
*                          check for memory overwrite
*
                           IF( CHECK ) THEN
                              CALL PSCHEKPAD( ICTXT, 'PSPORFS', NP,
     $                                        NQ, MEM( IPA0-IPREPAD ),
     $                                        DESCA( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSPORFS', NP,
     $                                        NQ, MEM( IPA-IPREPAD ),
     $                                        DESCA( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSPORFS', NP,
     $                                        MYRHS, MEM( IPB-IPREPAD ),
     $                                        DESCB( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSPORFS', NP,
     $                                        MYRHS,
     $                                        MEM( IPB0-IPREPAD ),
     $                                        DESCB( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSPORFS', 1,
     $                                        MYRHS,
     $                                        MEM( IPFERR-IPREPAD ), 1,
     $                                        IPREPAD, IPOSTPAD,
     $                                        PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSPORFS', 1,
     $                                        MYRHS,
     $                                        MEM( IPBERR-IPREPAD ), 1,
     $                                        IPREPAD, IPOSTPAD,
     $                                        PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSPORFS', LWORK,
     $                                        1, MEM( IPW-IPREPAD ),
     $                                        LWORK, IPREPAD, IPOSTPAD,
     $                                        PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSPORFS',
     $                                        LW2-IPOSTPAD, 1,
     $                                        MEM( IPW2-IPREPAD ),
     $                                        LW2-IPOSTPAD,
     $                                        IPREPAD, IPOSTPAD,
     $                                        PADVAL )
*
                              CALL PSFILLPAD( ICTXT, WORKSIZ-IPOSTPAD,
     $                                        1, MEM( IPW-IPREPAD ),
     $                                        WORKSIZ-IPOSTPAD, IPREPAD,
     $                                        IPOSTPAD, PADVAL )
*
*                             check the solution to rhs
*
                              CALL PSLASCHK( 'Symm', 'Diag', N, NRHS,
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
     $                                        WORKSIZ-IPOSTPAD,
     $                                        IPREPAD, IPOSTPAD,
     $                                        PADVAL )
                        END IF
                     END IF
*
*                    Gather maximum of all CPU and WALL clock timings
*
                     CALL SLCOMBINE( ICTXT, 'All', '>', 'W', 2, 1,
     $                               WTIME )
                     CALL SLCOMBINE( ICTXT, 'All', '>', 'C', 2, 1,
     $                               CTIME )
*
*                    Print results
*
                     IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
*
*                       1/3 N^3 + 1/2 N^2 flops for LLt factorization
*
                        NOPS = (DBLE(N)**3)/3.0D+0 +
     $                         (DBLE(N)**2)/2.0D+0
*
*                       nrhs * 2 N^2 flops for LLt solve.
*
                        NOPS = NOPS + 2.0D+0*(DBLE(N)**2)*DBLE(NRHS)
*
*                       Calculate total megaflops -- factorization and
*                       solve -- for WALL and CPU time, and print output
*
*                       Print WALL time if machine supports it
*
                        IF( WTIME( 1 ) + WTIME( 2 ) .GT. 0.0D+0 ) THEN
                           TMFLOPS = NOPS /
     $                            ( ( WTIME( 1 )+WTIME( 2 ) ) * 1.0D+6 )
                        ELSE
                           TMFLOPS = 0.0D+0
                        END IF
*
                        IF( WTIME( 2 ).GE.0.0D+0 )
     $                     WRITE( NOUT, FMT = 9993 ) 'WALL', UPLO, N,
     $                            NB, NRHS, NBRHS, NPROW, NPCOL,
     $                            WTIME( 1 ), WTIME( 2 ), TMFLOPS,
     $                            PASSED
*
*                       Print CPU time if machine supports it
*
                        IF( CTIME( 1 )+CTIME( 2 ).GT.0.0D+0 ) THEN
                           TMFLOPS = NOPS /
     $                            ( ( CTIME( 1 )+CTIME( 2 ) ) * 1.0D+6 )
                        ELSE
                           TMFLOPS = 0.0D+0
                        END IF
*
                        IF( CTIME( 2 ).GE.0.0D+0 )
     $                     WRITE( NOUT, FMT = 9993 ) 'CPU ', UPLO, N,
     $                            NB, NRHS, NBRHS, NPROW, NPCOL,
     $                            CTIME( 1 ), CTIME( 2 ), TMFLOPS,
     $                            PASSED
*
                     END IF
   10             CONTINUE
   20          CONTINUE
*
               IF( CHECK .AND. SRESID.GT.THRESH ) THEN
*
*                 Compute FRESID = ||A - LL'|| / (||A|| * N * eps)
*
                  CALL PSPOTRRV( UPLO, N, MEM( IPA ), 1, 1, DESCA,
     $                           MEM( IPW ) )
                  CALL PSLAFCHK( 'Symm', 'Diag', N, N, MEM( IPA ), 1, 1,
     $                           DESCA, IASEED, ANORM, FRESID,
     $                           MEM( IPW ) )
*
*                 Check for memory overwrite
*
                  CALL PSCHEKPAD( ICTXT, 'PSPOTRRV', NP, NQ,
     $                            MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                            IPREPAD, IPOSTPAD, PADVAL )
                  CALL PSCHEKPAD( ICTXT, 'PSGETRRV',
     $                            WORKSIZ-IPOSTPAD, 1,
     $                            MEM( IPW-IPREPAD ), WORKSIZ-IPOSTPAD,
     $                            IPREPAD, IPOSTPAD, PADVAL )
*
                  IF( IAM.EQ.0 ) THEN
                     IF( LSAME( UPLO, 'L' ) ) THEN
                        WRITE( NOUT, FMT = 9986 ) 'L*L''', FRESID
                     ELSE
                        WRITE( NOUT, FMT = 9986 ) 'U''*U', FRESID
                     END IF
                  END IF
               END IF
*
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
     $      CLOSE ( NOUT )
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
 9995 FORMAT( 'TIME UPLO     N  NB NRHS NBRHS    P    Q LLt Time ',
     $        'Slv Time   MFLOPS CHECK' )
 9994 FORMAT( '---- ---- ----- --- ---- ----- ---- ---- -------- ',
     $        '-------- -------- ------' )
 9993 FORMAT( A4, 4X, A1, 1X, I5, 1X, I3, 1X, I4, 1X, I5, 1X, I4, 1X,
     $        I4, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, A6 )
 9992 FORMAT( 'Finished ', I6, ' tests, with the following results:' )
 9991 FORMAT( I5, ' tests completed and passed residual checks.' )
 9990 FORMAT( I5, ' tests completed without checking.' )
 9989 FORMAT( I5, ' tests completed and failed residual checks.' )
 9988 FORMAT( I5, ' tests skipped because of illegal input values.' )
 9987 FORMAT( 'END OF TESTS.' )
 9986 FORMAT( '||A - ', A4, '|| / (||A|| * N * eps) = ', G25.7 )
 9985 FORMAT( '||Ax-b||/(||x||*||A||*eps*N) ', F25.7 )
*
      STOP
*
*     End of PSLLTDRIVER
*
      END
