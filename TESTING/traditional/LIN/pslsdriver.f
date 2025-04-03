      PROGRAM PSLSDRIVER
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 28, 2001
*
*  Purpose
*  =======
*
*  PSLSDRIVER is the main test program for the REAL
*  SCALAPACK (full rank) Least Squares routines. This test driver solves
*  full-rank least square problems.
*
*  The program must be driven by a short data file.  An annotated
*  example of a data file can be obtained by deleting the first 3
*  characters from the following 17 lines:
*  'ScaLapack LS solve input file'
*  'Intel iPSC/860 hypercube, gamma model.'
*  'LS.out'                    output file name (if any)
*  6                           device out
*  4                           number of problems sizes
*  55 17 31 201                values of M
*  5 71 31 201                 values of N
*  3                           number of NB's
*  2 3 5                       values of NB
*  3                           number of NRHS's
*  2 3 5                       values of NRHS
*  2                           number of NBRHS's
*  1 2                         values of NBRHS
*  7                           number of process grids (ordered P & Q)
*  1 2 1 4 2 3 8               values of P
*  7 2 4 1 3 2 1               values of Q
*  3.0                         threshold
*
*  Internal Parameters
*  ===================
*
*  TOTMEM   INTEGER, default = 6200000.
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
*  INTGSZ   INTEGER, default = 4 bytes.
*  REALSZ   INTEGER, default = 4 bytes.
*           INTGSZ and REALSZ indicate the length in bytes on the
*           given platform for an integer and a single precision real.
*  MEM      REAL array, dimension ( TOTMEM / REALSZ )
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
      REAL               ONE, ZERO
      PARAMETER          ( REALSZ = 4, TOTMEM = 2000000,
     $                     MEMSIZ = TOTMEM / REALSZ, NTESTS = 20,
     $                     PADVAL = -9923.0E+0 )
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            CHECK, TPSD
      CHARACTER          TRANS
      CHARACTER*6        PASSED
      CHARACTER*80       OUTFILE
      INTEGER            HH, I, IAM, IASEED, IBSEED, ICTXT, II, IMIDPAD,
     $                   INFO, IPA, IPB, IPOSTPAD, IPREPAD, IPW, IPW2,
     $                   IPX, ISCALE, ITRAN, ITYPE, J, JJ, K, KFAIL, KK,
     $                   KPASS, KSKIP, KTESTS, LCM, LCMP, LTAU, LWF,
     $                   LWORK, LWS, M, MNP, MNRHSP, MP, MQ, MYCOL,
     $                   MYROW, N, NB, NBRHS, NCOLS, NGRIDS, NMAT, NNB,
     $                   NNBR, NNR, NNRHSQ, NOUT, NP, NPCOL, NPROCS,
     $                   NPROW, NROWS, NQ, NRHS, NRHSP, NRHSQ, WORKSIZ
      REAL               ANORM, BNORM, SRESID, THRESH
      DOUBLE PRECISION   ADDFAC, ADDS, MULFAC, MULTS, NOPS, TMFLOPS
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), DESCW( LLD_ ),
     $                   DESCX( DLEN_ ), IERR( 2 ), MVAL( NTESTS ),
     $                   NBRVAL( NTESTS ), NBVAL( NTESTS ),
     $                   NRVAL( NTESTS ), NVAL( NTESTS ),
     $                   PVAL( NTESTS ), QVAL( NTESTS )
      REAL               MEM( MEMSIZ ), RESULT( 2 )
      DOUBLE PRECISION   CTIME( 1 ), WTIME( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_EXIT, BLACS_GET,
     $                   BLACS_GRIDEXIT, BLACS_GRIDINFO, BLACS_GRIDINIT,
     $                   BLACS_PINFO, DESCINIT, IGSUM2D, PSCHEKPAD,
     $                   PSFILLPAD, PSGELS, PSGEMM, PSLACPY,
     $                   PSLSINFO, PSMATGEN, PSNRM2, PSSCAL,
     $                   PSQRT13, PSQRT16, SLBOOT, SLCOMBINE, SLTIMER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, ILCM, NUMROC
      REAL               PSLANGE, PSQRT14, PSQRT17
      EXTERNAL           ICEIL, ILCM, LSAME, NUMROC, PSLANGE,
     $                   PSQRT14, PSQRT17
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Data Statements ..
      DATA               KTESTS, KPASS, KFAIL, KSKIP / 4*0 /
*     ..
*     .. Executable Statements ..
*
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
*
      IASEED = 100
      IBSEED = 200
      CALL PSLSINFO( OUTFILE, NOUT, NMAT, MVAL, NTESTS, NVAL,
     $               NTESTS, NNB, NBVAL, NTESTS, NNR, NRVAL, NTESTS,
     $               NNBR, NBRVAL, NTESTS, NGRIDS, PVAL, NTESTS, QVAL,
     $               NTESTS, THRESH, MEM, IAM, NPROCS )
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
      DO 90 I = 1, NGRIDS
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
            GO TO 90
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
         IF( ( MYROW.GE.NPROW ).OR.( MYCOL.GE.NPCOL ) )
     $      GO TO 90
*
         DO 80 J = 1, NMAT
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
               GO TO 80
            END IF
*
*           Loop over different blocking sizes
*
            DO 70 K = 1, NNB
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
                  GO TO 70
               END IF
*
*              Padding constants
*
               MP = NUMROC( M, NB, MYROW, 0, NPROW )
               MQ = NUMROC( M, NB, MYCOL, 0, NPCOL )
               NP = NUMROC( N, NB, MYROW, 0, NPROW )
               MNP = MAX( MP, NP )
               NQ = NUMROC( N, NB, MYCOL, 0, NPCOL )
*
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
     $                        MAX( 1, MP ) + IMIDPAD, IERR( 1 ) )
*
*              Check all processes for an error
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
               IF( IERR( 1 ).LT.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 ) 'descriptor'
                  KSKIP = KSKIP + 1
                  GO TO 70
               END IF
*
               DO 60 ISCALE = 1, 3
*
                  ITYPE = ISCALE
*
*                 Assign pointers into MEM for SCALAPACK arrays, A is
*                 allocated starting at position MEM( IPREPAD+1 )
*
                  IPA = IPREPAD + 1
                  IPX = IPA + DESCA( LLD_ )*NQ + IPOSTPAD + IPREPAD
                  IPW = IPX
*
                  WORKSIZ = NQ + IPOSTPAD
*
*                 Check for adequate memory for problem size
*
                  IERR( 1 ) = 0
                  IF( ( IPW+WORKSIZ ).GT.MEMSIZ ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9996 ) 'MEMORY',
     $                         ( IPX+WORKSIZ )*REALSZ
                     IERR( 1 ) = 1
                  END IF
*
*                 Check all processes for an error
*
                  CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                          0 )
*
                  IF( IERR( 1 ).GT.0 ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9997 ) 'MEMORY'
                     KSKIP = KSKIP + 1
                     GO TO 70
                  END IF
*
                  IF( CHECK ) THEN
                     CALL PSFILLPAD( ICTXT, MP, NQ, MEM( IPA-IPREPAD ),
     $                               DESCA( LLD_ ), IPREPAD, IPOSTPAD,
     $                               PADVAL )
                     CALL PSFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                               MEM( IPW-IPREPAD ),
     $                               WORKSIZ-IPOSTPAD, IPREPAD,
     $                               IPOSTPAD, PADVAL )
                  END IF
*
*                 Generate the matrix A and calculate its 1-norm
*
                  CALL PSQRT13( ISCALE, M, N, MEM( IPA ), 1, 1,
     $                          DESCA, ANORM, IASEED, MEM( IPW ) )
*
                  IF( CHECK ) THEN
                     CALL PSCHEKPAD( ICTXT, 'PSQRT13', MP, NQ,
     $                               MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     CALL PSCHEKPAD( ICTXT, 'PSQRT13',
     $                               WORKSIZ-IPOSTPAD, 1,
     $                               MEM( IPW-IPREPAD ),
     $                               WORKSIZ-IPOSTPAD, IPREPAD,
     $                               IPOSTPAD, PADVAL )
                  END IF
*
                  DO 50 ITRAN = 1, 2
*
                     IF( ITRAN.EQ.1 ) THEN
                        NROWS = M
                        NCOLS = N
                        TRANS = 'N'
                        TPSD  = .FALSE.
                     ELSE
                        NROWS = N
                        NCOLS = M
                        TRANS = 'T'
                        TPSD  = .TRUE.
                     END IF
*
*                    Loop over the different values for NRHS
*
                     DO 40 HH =  1, NNR
*
                        NRHS = NRVAL( HH )
*
                        DO 30 KK = 1, NNBR
*
                           NBRHS = NBRVAL( KK )
*
                           NRHSP = NUMROC( NRHS, NBRHS, MYROW, 0,
     $                                     NPROW )
                           NRHSQ = NUMROC( NRHS, NBRHS, MYCOL, 0,
     $                                     NPCOL )
*
*                          Define Array descriptor for rhs MAX(M,N)xNRHS
*
                           CALL DESCINIT( DESCX, MAX( M, N ), NRHS, NB,
     $                                    NBRHS, 0, 0, ICTXT,
     $                                    MAX( 1, MNP ) + IMIDPAD,
     $                                    IERR( 1 ) )
                           IF( TPSD ) THEN
                              CALL DESCINIT( DESCW, M, NRHS, NB, NBRHS,
     $                                       0, 0, ICTXT, MAX( 1, MP ) +
     $                                       IMIDPAD, IERR( 2 ) )
                           ELSE
                              CALL DESCINIT( DESCW, N, NRHS, NB, NBRHS,
     $                                       0, 0, ICTXT, MAX( 1, NP ) +
     $                                       IMIDPAD, IERR( 2 ) )
                           END IF
*
*                          Check all processes for an error
*
                           CALL IGSUM2D( ICTXT, 'All', ' ', 2, 1, IERR,
     $                                   2, -1, 0 )
*
                           IF( IERR( 1 ).LT.0 .OR. IERR( 2 ).LT.0 ) THEN
                              IF( IAM.EQ.0 )
     $                           WRITE( NOUT, FMT = 9997 ) 'descriptor'
                              KSKIP = KSKIP + 1
                              GO TO 30
                           END IF
*
*                          Check for enough memory
*
                           IPX = IPA + DESCA( LLD_ )*NQ + IPOSTPAD +
     $                           IPREPAD
                           IPW = IPX + DESCX( LLD_ )*NRHSQ + IPOSTPAD +
     $                           IPREPAD
                           WORKSIZ = DESCW( LLD_ )*NRHSQ + IPOSTPAD
*
                           IERR( 1 ) = 0
                           IF( IPW+WORKSIZ.GT.MEMSIZ ) THEN
                              IF( IAM.EQ.0 )
     $                           WRITE( NOUT, FMT = 9996 ) 'Generation',
     $                                  ( IPW+WORKSIZ )*REALSZ
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
     $                           WRITE( NOUT, FMT = 9997 ) 'MEMORY'
                              KSKIP = KSKIP + 1
                              GO TO 30
                           END IF
*
*                          Generate RHS
*
                           IF( TPSD ) THEN
                              CALL PSMATGEN( ICTXT, 'No', 'No',
     $                                       DESCW( M_ ), DESCW( N_ ),
     $                                       DESCW( MB_ ), DESCW( NB_ ),
     $                                       MEM( IPW ), DESCW( LLD_ ),
     $                                       DESCW( RSRC_ ),
     $                                       DESCW( CSRC_ ), IBSEED, 0,
     $                                       MP, 0, NRHSQ, MYROW, MYCOL,
     $                                       NPROW, NPCOL )
                           ELSE
                              CALL PSMATGEN( ICTXT, 'No', 'No',
     $                                       DESCW( M_ ), DESCW( N_ ),
     $                                       DESCW( MB_ ), DESCW( NB_ ),
     $                                       MEM( IPW ), DESCW( LLD_ ),
     $                                       DESCW( RSRC_ ),
     $                                       DESCW( CSRC_ ), IBSEED, 0,
     $                                       NP, 0, NRHSQ, MYROW, MYCOL,
     $                                       NPROW, NPCOL )
                           END IF
*
                           IF( CHECK ) THEN
                              CALL PSFILLPAD( ICTXT, MNP, NRHSQ,
     $                                        MEM( IPX-IPREPAD ),
     $                                        DESCX( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              IF( TPSD ) THEN
                                 CALL PSFILLPAD( ICTXT, MP, NRHSQ,
     $                                           MEM( IPW-IPREPAD ),
     $                                           DESCW( LLD_ ), IPREPAD,
     $                                           IPOSTPAD, PADVAL )
                              ELSE
                                 CALL PSFILLPAD( ICTXT, NP, NRHSQ,
     $                                           MEM( IPW-IPREPAD ),
     $                                           DESCW( LLD_ ), IPREPAD,
     $                                           IPOSTPAD, PADVAL )
                              END IF
                           END IF
*
                           DO 10 JJ = 1, NRHS
                              CALL PSNRM2( NCOLS, BNORM, MEM( IPW ), 1,
     $                                     JJ, DESCW, 1 )
                              IF( BNORM.GT.ZERO )
     $                           CALL PSSCAL( NCOLS, ONE / BNORM,
     $                                        MEM( IPW ), 1, JJ, DESCW,
     $                                        1 )
   10                      CONTINUE
*
                           CALL PSGEMM( TRANS, 'N', NROWS, NRHS, NCOLS,
     $                                  ONE, MEM( IPA ), 1, 1, DESCA,
     $                                  MEM( IPW ), 1, 1, DESCW, ZERO,
     $                                  MEM( IPX ), 1, 1, DESCX )
*
                           IF( CHECK ) THEN
*
*                             check for memory overwrite
*
                              CALL PSCHEKPAD( ICTXT, 'Generation', MP,
     $                                        NQ, MEM( IPA-IPREPAD ),
     $                                        DESCA( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'Generation', MNP,
     $                                        NRHSQ, MEM( IPX-IPREPAD ),
     $                                        DESCX( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              IF( TPSD ) THEN
                                 CALL PSCHEKPAD( ICTXT, 'Generation',
     $                                           MP, NRHSQ,
     $                                           MEM( IPW-IPREPAD ),
     $                                           DESCW( LLD_ ), IPREPAD,
     $                                           IPOSTPAD, PADVAL )
                              ELSE
                                 CALL PSCHEKPAD( ICTXT, 'Generation',
     $                                           NP, NRHSQ,
     $                                           MEM( IPW-IPREPAD ),
     $                                           DESCW( LLD_ ), IPREPAD,
     $                                           IPOSTPAD, PADVAL )
                              END IF
*
*                             Allocate space for copy of rhs
*
                              IPB = IPW
*
                              IF( TPSD ) THEN
                                 CALL DESCINIT( DESCB, N, NRHS, NB,
     $                                     NBRHS, 0, 0, ICTXT,
     $                                     MAX( 1, NP ) + IMIDPAD,
     $                                     IERR( 1 ) )
                              ELSE
                                 CALL DESCINIT( DESCB, M, NRHS, NB,
     $                                     NBRHS, 0, 0, ICTXT,
     $                                     MAX( 1, MP ) + IMIDPAD,
     $                                     IERR( 1 ) )
                              END IF
*
*                             Check all processes for an error
*
                              CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1,
     $                                      IERR, 1, -1, 0 )
*
                              IF( IERR( 1 ).LT.0 ) THEN
                                 IF( IAM.EQ.0 )
     $                              WRITE( NOUT, FMT = 9997 )
     $                                     'descriptor'
                                 KSKIP = KSKIP + 1
                                 GO TO 30
                              END IF
*
                              IPW = IPB + DESCB( LLD_ )*NRHSQ +
     $                              IPOSTPAD + IPREPAD
*
                           END IF
*
*                          Calculate the amount of workspace for PSGELS
*
                           IF( M.GE.N ) THEN
                              LTAU = NUMROC( MIN(M,N), NB, MYCOL, 0,
     $                                       NPCOL )
                              LWF  = NB * ( MP + NQ + NB )
                              LWS  = MAX( ( NB*( NB - 1 ) ) / 2,
     $                                    ( MP + NRHSQ ) * NB ) + NB*NB
                           ELSE
                              LCM = ILCM( NPROW, NPCOL )
                              LCMP = LCM / NPROW
                              LTAU = NUMROC( MIN(M,N), NB, MYROW, 0,
     $                                       NPROW )
                              LWF  = NB * ( MP + NQ + NB )
                              LWS  = MAX( ( NB*( NB - 1 ) ) / 2, ( NP +
     $                               MAX( NQ + NUMROC( NUMROC( N, NB, 0,
     $                               0, NPROW ), NB, 0, 0, LCMP ),
     $                               NRHSQ ) ) * NB ) + NB*NB
                           END IF
*
                           LWORK = LTAU + MAX( LWF, LWS )
                           WORKSIZ = LWORK + IPOSTPAD
*
*                          Check for adequate memory for problem size
*
                           IERR( 1 ) = 0
                           IF( IPW+WORKSIZ.GT.MEMSIZ ) THEN
                              IF( IAM.EQ.0 )
     $                           WRITE( NOUT, FMT = 9996 ) 'solve',
     $                                  ( IPW+WORKSIZ )*REALSZ
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
     $                           WRITE( NOUT, FMT = 9997 ) 'MEMORY'
                              KSKIP = KSKIP + 1
                              GO TO 30
                           END IF
*
                           IF( CHECK ) THEN
*
*                             Make the copy of the right hand side
*
                              CALL PSLACPY( 'All', NROWS, NRHS,
     $                                      MEM( IPX ), 1, 1, DESCX,
     $                                      MEM( IPB ), 1, 1, DESCB )
*
                              IF( TPSD ) THEN
                                 CALL PSFILLPAD( ICTXT, NP, NRHSQ,
     $                                           MEM( IPB-IPREPAD ),
     $                                           DESCB( LLD_ ), IPREPAD,
     $                                           IPOSTPAD, PADVAL )
                              ELSE
                                 CALL PSFILLPAD( ICTXT, MP, NRHSQ,
     $                                           MEM( IPB-IPREPAD ),
     $                                           DESCB( LLD_ ), IPREPAD,
     $                                           IPOSTPAD, PADVAL )
                              END IF
                              CALL PSFILLPAD( ICTXT, LWORK, 1,
     $                                        MEM( IPW-IPREPAD ),
     $                                        LWORK, IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                           END IF
*
                           CALL SLBOOT( )
                           CALL BLACS_BARRIER( ICTXT, 'All' )
                           CALL SLTIMER( 1 )
*
*                          Solve the LS or overdetermined system
*
                           CALL PSGELS( TRANS, M, N, NRHS, MEM( IPA ),
     $                                  1, 1, DESCA, MEM( IPX ), 1, 1,
     $                                  DESCX, MEM( IPW ), LWORK, INFO )
*
                           CALL SLTIMER( 1 )
*
                           IF( CHECK ) THEN
*
*                             check for memory overwrite
*
                              CALL PSCHEKPAD( ICTXT, 'PSGELS', MP,
     $                                        NQ, MEM( IPA-IPREPAD ),
     $                                        DESCA( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSGELS', MNP,
     $                                        NRHSQ, MEM( IPX-IPREPAD ),
     $                                        DESCX( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSGELS', LWORK,
     $                                        1, MEM( IPW-IPREPAD ),
     $                                        LWORK, IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                           END IF
*
*                          Regenerate A in place for testing and next
*                          iteration
*
                           CALL PSQRT13( ISCALE, M, N, MEM( IPA ), 1, 1,
     $                                   DESCA, ANORM, IASEED,
     $                                   MEM( IPW ) )
*
*                          check the solution to rhs
*
                           IF( CHECK ) THEN
*
*                             Am I going to call PSQRT17 ?
*
                              IF( ( M.GE.N .AND. ( .NOT.TPSD ) ) .OR.
     $                            ( M.LT.N .AND. TPSD ) ) THEN
*
*                                Call PSQRT17 first, A, X, and B remain
*                                unchanged.  Solving LS system
*
*                                Check amount of memory for PSQRT17
*
                                 IF( TPSD ) THEN
                                    WORKSIZ = NP*NRHSQ + NRHSP*MQ
                                    IPW2 = IPW + WORKSIZ
                                    WORKSIZ = WORKSIZ +
     $                                     MAX( NQ, MAX( MQ, NRHSQ ) ) +
     $                                     IPOSTPAD
                                 ELSE
                                    WORKSIZ = MP*NRHSQ + NRHSP*NQ
                                    IPW2 = IPW + WORKSIZ
                                    WORKSIZ = WORKSIZ +
     $                                     MAX( NQ, NRHSQ ) +
     $                                     IPOSTPAD
                                 END IF
*
*                                Check for adequate memory for problem
*                                size
*
                                 IERR( 1 ) = 0
                                 IF( ( IPW+WORKSIZ ).GT.MEMSIZ ) THEN
                                    IF( IAM.EQ.0 )
     $                                 WRITE( NOUT, FMT = 9996 )
     $                                  'MEMORY', ( IPW+WORKSIZ )*REALSZ
                                   IERR( 1 ) = 1
                                 END IF
*
*                                Check all processes for an error
*
                                 CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1,
     $                                         IERR, 1, -1, 0 )
*
                                 IF( IERR( 1 ).GT.0 ) THEN
                                    IF( IAM.EQ.0 )
     $                                 WRITE( NOUT, FMT = 9997 )
     $                                        'MEMORY'
                                    KSKIP = KSKIP + 1
                                    GO TO 30
                                 END IF
*
                                 CALL PSFILLPAD( ICTXT,
     $                                           WORKSIZ-IPOSTPAD, 1,
     $                                           MEM( IPW-IPREPAD ),
     $                                           WORKSIZ-IPOSTPAD,
     $                                           IPREPAD, IPOSTPAD,
     $                                           PADVAL )
*
                                 RESULT( 2 ) = PSQRT17( TRANS, 1, M, N,
     $                                                  NRHS,
     $                                                  MEM( IPA ),
     $                                                  1, 1, DESCA,
     $                                                  MEM( IPX ), 1,
     $                                                  1, DESCX,
     $                                                  MEM( IPB ),
     $                                                  1, 1, DESCB,
     $                                                  MEM( IPW ),
     $                                                  MEM( IPW2 ) )
                                 SRESID = RESULT( 2 )
*
                                 CALL PSCHEKPAD( ICTXT, 'PSQRT17',
     $                                           MP, NQ,
     $                                           MEM( IPA-IPREPAD ),
     $                                           DESCA( LLD_ ),
     $                                           IPREPAD, IPOSTPAD,
     $                                           PADVAL )
                                 CALL PSCHEKPAD( ICTXT, 'PSQRT17',
     $                                           MNP, NRHSQ,
     $                                           MEM( IPX-IPREPAD ),
     $                                           DESCX( LLD_ ), IPREPAD,
     $                                           IPOSTPAD, PADVAL )
                                 IF( TPSD ) THEN
                                    CALL PSCHEKPAD( ICTXT, 'PSQRT17',
     $                                              NP, NRHSQ,
     $                                              MEM( IPB-IPREPAD ),
     $                                              DESCB( LLD_ ),
     $                                              IPREPAD, IPOSTPAD,
     $                                              PADVAL )
                                 ELSE
                                    CALL PSCHEKPAD( ICTXT, 'PSQRT17',
     $                                              MP, NRHSQ,
     $                                              MEM( IPB-IPREPAD ),
     $                                              DESCB( LLD_ ),
     $                                              IPREPAD, IPOSTPAD,
     $                                              PADVAL )
                                 END IF
                                 CALL PSCHEKPAD( ICTXT, 'PSQRT17',
     $                                           WORKSIZ-IPOSTPAD, 1,
     $                                           MEM( IPW-IPREPAD ),
     $                                           WORKSIZ-IPOSTPAD,
     $                                           IPREPAD, IPOSTPAD,
     $                                           PADVAL )
                              END IF
*
*                             Call PSQRT16, B will be destroyed.
*
                              IF( TPSD ) THEN
                                 WORKSIZ = MP + IPOSTPAD
                              ELSE
                                 WORKSIZ = NQ + IPOSTPAD
                              END IF
*
*                             Check for adequate memory for problem size
*
                              IERR( 1 ) = 0
                              IF( ( IPW+WORKSIZ ).GT.MEMSIZ ) THEN
                                 IF( IAM.EQ.0 )
     $                              WRITE( NOUT, FMT = 9996 ) 'MEMORY',
     $                                    ( IPW+WORKSIZ )*REALSZ
                                IERR( 1 ) = 1
                              END IF
*
*                             Check all processes for an error
*
                              CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1,
     $                                 IERR, 1, -1, 0 )
*
                              IF( IERR( 1 ).GT.0 ) THEN
                                 IF( IAM.EQ.0 )
     $                              WRITE( NOUT, FMT = 9997 ) 'MEMORY'
                                 KSKIP = KSKIP + 1
                                 GO TO 30
                              END IF
*
                              CALL PSFILLPAD( ICTXT,
     $                                        WORKSIZ-IPOSTPAD, 1,
     $                                        MEM( IPW-IPREPAD ),
     $                                        WORKSIZ-IPOSTPAD,
     $                                        IPREPAD, IPOSTPAD,
     $                                        PADVAL )
*
                              CALL PSQRT16( TRANS, M, N, NRHS,
     $                                      MEM( IPA ), 1, 1, DESCA,
     $                                      MEM( IPX ), 1, 1, DESCX,
     $                                      MEM( IPB ), 1, 1, DESCB,
     $                                      MEM( IPW ), RESULT( 1 ) )
*
                              CALL PSCHEKPAD( ICTXT, 'PSQRT16',
     $                                        MP, NQ,
     $                                        MEM( IPA-IPREPAD ),
     $                                        DESCA( LLD_ ),
     $                                        IPREPAD, IPOSTPAD,
     $                                        PADVAL )
                              CALL PSCHEKPAD( ICTXT, 'PSQRT16',
     $                                        MNP, NRHSQ,
     $                                        MEM( IPX-IPREPAD ),
     $                                        DESCX( LLD_ ), IPREPAD,
     $                                        IPOSTPAD, PADVAL )
                              IF( TPSD ) THEN
                                 CALL PSCHEKPAD( ICTXT, 'PSQRT16',
     $                                           NP, NRHSQ,
     $                                           MEM( IPB-IPREPAD ),
     $                                           DESCB( LLD_ ),
     $                                           IPREPAD, IPOSTPAD,
     $                                           PADVAL )
                              ELSE
                                 CALL PSCHEKPAD( ICTXT, 'PSQRT16',
     $                                           MP, NRHSQ,
     $                                           MEM( IPB-IPREPAD ),
     $                                           DESCB( LLD_ ),
     $                                           IPREPAD, IPOSTPAD,
     $                                           PADVAL )
                              END IF
                              CALL PSCHEKPAD( ICTXT, 'PSQRT16',
     $                                        WORKSIZ-IPOSTPAD, 1,
     $                                        MEM( IPW-IPREPAD ),
     $                                        WORKSIZ-IPOSTPAD,
     $                                        IPREPAD, IPOSTPAD,
     $                                        PADVAL )
*
*                             Call PSQRT14
*
                              IF( ( M.GE.N .AND. TPSD ) .OR.
     $                            ( M.LT.N .AND. ( .NOT.TPSD ) ) ) THEN
*
                                 IPW = IPB
*
                                 IF( TPSD ) THEN
*
                                    NNRHSQ = NUMROC( N+NRHS, NB, MYCOL,
     $                                               0, NPCOL )
                                    LTAU = NUMROC( MIN( M, N+NRHS ), NB,
     $                                             MYCOL, 0, NPCOL )
                                    LWF = NB * ( NB + MP + NNRHSQ )
                                    WORKSIZ = MP * NNRHSQ + LTAU + LWF +
     $                                        IPOSTPAD
*
                                 ELSE
*
                                    MNRHSP = NUMROC( M+NRHS, NB, MYROW,
     $                                               0, NPROW )
                                    LTAU = NUMROC( MIN( M+NRHS, N ), NB,
     $                                             MYROW, 0, NPROW )
                                    LWF = NB * ( NB + MNRHSP + NQ )
                                    WORKSIZ = MNRHSP * NQ + LTAU + LWF +
     $                                        IPOSTPAD
*
                                 END IF
*
*                                Check for adequate memory for problem
*                                size
*
                                 IERR( 1 ) = 0
                                 IF( ( IPW+WORKSIZ ).GT.MEMSIZ ) THEN
                                    IF( IAM.EQ.0 )
     $                                 WRITE( NOUT, FMT = 9996 )
     $                                 'MEMORY', ( IPW+WORKSIZ )*REALSZ
                                    IERR( 1 ) = 1
                                 END IF
*
*                                Check all processes for an error
*
                                 CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1,
     $                                         IERR, 1, -1, 0 )
*
                                 IF( IERR( 1 ).GT.0 ) THEN
                                    IF( IAM.EQ.0 )
     $                                 WRITE( NOUT, FMT = 9997 )
     $                                       'MEMORY'
                                    KSKIP = KSKIP + 1
                                    GO TO 30
                                 END IF
*
                                 CALL PSFILLPAD( ICTXT,
     $                                           WORKSIZ-IPOSTPAD, 1,
     $                                           MEM( IPW-IPREPAD ),
     $                                           WORKSIZ-IPOSTPAD,
     $                                           IPREPAD, IPOSTPAD,
     $                                           PADVAL )
*
*                                Solve underdetermined system
*
                                 RESULT( 2 ) = PSQRT14( TRANS, M, N,
     $                                                  NRHS,
     $                                                  MEM( IPA ), 1,
     $                                                  1, DESCA,
     $                                                  MEM( IPX ),
     $                                                  1, 1, DESCX,
     $                                                  MEM( IPW ) )
                                 SRESID = RESULT( 2 )
*
                                 CALL PSCHEKPAD( ICTXT, 'PSQRT14',
     $                                           MP, NQ,
     $                                           MEM( IPA-IPREPAD ),
     $                                           DESCA( LLD_ ),
     $                                           IPREPAD, IPOSTPAD,
     $                                           PADVAL )
                                 CALL PSCHEKPAD( ICTXT, 'PSQRT14',
     $                                           MNP, NRHSQ,
     $                                           MEM( IPX-IPREPAD ),
     $                                           DESCX( LLD_ ), IPREPAD,
     $                                           IPOSTPAD, PADVAL )
                                 CALL PSCHEKPAD( ICTXT, 'PSQRT14',
     $                                           WORKSIZ-IPOSTPAD, 1,
     $                                           MEM( IPW-IPREPAD ),
     $                                           WORKSIZ-IPOSTPAD,
     $                                           IPREPAD, IPOSTPAD,
     $                                           PADVAL )
                              END IF
*
*                             Print information about the tests that
*                             did not pass the threshold.
*
                              PASSED = 'PASSED'
                              DO 20 II = 1, 2
                                 IF( ( RESULT( II ).GE.THRESH ) .AND.
     $                             ( RESULT( II )-RESULT( II ).EQ.0.0E+0
     $                              ) ) THEN
                                    IF( IAM.EQ.0 )
     $                                 WRITE( NOUT, FMT = 9986 )TRANS,
     $                                 M, N, NRHS, NB, ITYPE, II,
     $                                 RESULT( II )
                                    KFAIL = KFAIL + 1
                                    PASSED = 'FAILED'
                                 ELSE
                                    KPASS = KPASS + 1
                                 END IF
   20                         CONTINUE
*
                           ELSE
*
*                             By-pass the solve check
*
                              KPASS = KPASS + 1
                              SRESID = SRESID - SRESID
                              PASSED = 'BYPASS'
*
                           END IF
*
*                          Gather maximum of all CPU and WALL clock
*                          timings
*
                           CALL SLCOMBINE( ICTXT, 'All', '>', 'W', 1, 1,
     $                                     WTIME )
                           CALL SLCOMBINE( ICTXT, 'All', '>', 'C', 1, 1,
     $                                     CTIME )
*
*                          Print results
*
                           IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                              ADDFAC = 1
                              MULFAC = 1
                              IF( M.GE.N ) THEN
*
*                                NOPS = SOPLA( 'SGEQRF', M, N, 0, 0,
*                                NB ) + SOPLA( 'SORMQR', M, NRHS, N,
*                                0, NB )
*
                                 MULTS = N*( ( ( 23.D0 / 6.D0 )+M+N /
     $                                   2.D0 )+ N*( M-N / 3.D0 ) ) +
     $                                   N*NRHS*( 2.D0*M+2.D0-N )
                                 ADDS = N*( ( 5.D0 / 6.D0 )+N*
     $                                  ( 1.D0 / 2.D0+( M-N / 3.D0 ) ) )
     $                                  + N*NRHS*( 2.D0*M+1.D0-N )
                              ELSE
*
*                                NOPS = SOPLA( 'SGELQF', M, N, 0, 0,
*                                       NB ) + SOPLA( 'SORMLQ', M,
*                                       NRHS, N, 0, NB )
*
                                 MULTS = M*( ( ( 29.D0 / 6.D0 )+2.D0*N-M
     $                                   / 2.D0 )+M*( N-M / 3.D0 ) )
     $                                   + N*NRHS*( 2.D0*M+2.D0-N )
                                 ADDS = M*( ( 5.D0 / 6.D0 )+M / 2.D0+M*
     $                                  ( N-M / 3.D0 ) )
     $                                  + N*NRHS*( 2.D0*M+1.D0-N )
                              END IF
                              NOPS = ADDFAC*ADDS + MULFAC*MULTS
*
*                             Calculate total megaflops, for WALL and
*                             CPU time, and print output
*
*                             Print WALL time if machine supports it
*
                              IF( WTIME( 1 ).GT.0.0D+0 ) THEN
                                 TMFLOPS = NOPS / ( WTIME( 1 )*1.0D+6 )
                              ELSE
                                 TMFLOPS = 0.0D+0
                              END IF
*
                              IF( WTIME( 1 ).GE.0.0D+0 )
     $                           WRITE( NOUT, FMT = 9993 )
     $                                  'WALL', TRANS, M, N, NB, NRHS,
     $                                  NBRHS, NPROW, NPCOL, WTIME( 1 ),
     $                                  TMFLOPS, PASSED
*
*                             Print CPU time if machine supports it
*
                              IF( CTIME( 1 ).GT.0.0D+0 ) THEN
                                 TMFLOPS = NOPS / ( CTIME( 1 )*1.0D+6 )
                              ELSE
                                 TMFLOPS = 0.0D+0
                              END IF
*
                              IF( CTIME( 1 ).GE.0.0D+0 )
     $                           WRITE( NOUT, FMT = 9993 )
     $                                  'CPU ', TRANS, M, N, NB, NRHS,
     $                                  NBRHS, NPROW, NPCOL, CTIME( 1 ),
     $                                  TMFLOPS, PASSED
                           END IF
   30                   CONTINUE
   40                CONTINUE
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
         CALL BLACS_GRIDEXIT( ICTXT )
   90 CONTINUE
*
*     Print out ending messages and close output file
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
 9995 FORMAT( 'Time TRANS      M     N   NB  NRHS NBRHS     P     Q ',
     $        'LS Time     MFLOPS  CHECK' )
 9994 FORMAT( '---- ----- ------ ------ --- ----- ----- ----- ----- ',
     $        '--------- -------- ------' )
 9993 FORMAT( A4, 3X, A1, 3X, I6, 1X, I6, 1X, I3, 1X, I5, 1X, I5, 1X,
     $        I5, 1X, I5, 1X, F9.2, 1X, F8.2, 1X, A6 )
 9992 FORMAT( 'Finished', I6, ' tests, with the following results:' )
 9991 FORMAT( I5, ' tests completed and passed residual checks.' )
 9990 FORMAT( I5, ' tests completed without checking.' )
 9989 FORMAT( I5, ' tests completed and failed residual checks.' )
 9988 FORMAT( I5, ' tests skipped because of illegal input values.' )
 9987 FORMAT( 'END OF TESTS.' )
 9986 FORMAT( ' TRANS=''', A1, ''', M=', I5, ', N=', I5, ', NRHS=', I4,
     $      ', NB=', I4, ', type', I2, ', test(', I2, ')=', G12.5 )
*
      STOP
*
*     End of PSLSDRIVER
*
      END
