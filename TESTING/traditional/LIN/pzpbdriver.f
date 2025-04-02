      PROGRAM PZPBDRIVER
*
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*  Purpose
*  =======
*
*  PZPBDRIVER is a test program for the
*  ScaLAPACK Band Cholesky routines corresponding to the options
*  indicated by ZPB.  This test driver performs an
*  A = L*L**H factorization
*  and solves a linear system with the factors for 1 or more RHS.
*
*  The program must be driven by a short data file.
*  Here's an example file:
*'ScaLAPACK, Version 1.2, banded linear systems input file'
*'PVM.'
*''                              output file name (if any)
*6                               device out
*'L'                             define Lower or Upper
*9                               number of problem sizes
*1 5 17 28 37 121 200 1023 2048 3073     values of N
*6                               number of bandwidths
*1 2 4 10 31 64                values of BW
*1                               number of NB's
*-1 3 4 5                        values of NB (-1 for automatic choice)
*1                               number of NRHS's (must be 1)
*8                               values of NRHS
*1                               number of NBRHS's (ignored)
*1                               values of NBRHS (ignored)
*6                               number of process grids
*1 2 3 4 5 7 8 15 26 47 64       values of "Number of Process Columns"
*3.0                             threshold
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
*
*  INTGSZ   INTEGER, default = 4 bytes.
*  ZPLXSZ   INTEGER, default = 16 bytes.
*           INTGSZ and ZPLXSZ indicate the length in bytes on the
*           given platform for an integer and a double precision
*           complex.
*  MEM      COMPLEX*16 array, dimension ( TOTMEM/ZPLXSZ )
*           All arrays used by ScaLAPACK routines are allocated from
*           this array and referenced by pointers.  The integer IPB,
*           for example, is a pointer to the starting element of MEM for
*           the solution vector(s) B.
*
*  =====================================================================
*
*  Code Developer: Andrew J. Cleary, University of Tennessee.
*    Current address: Lawrence Livermore National Labs.
*  This version released: August, 2001.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            TOTMEM
      PARAMETER          ( TOTMEM = 3000000 )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*
      DOUBLE PRECISION   ZERO
      INTEGER            MEMSIZ, NTESTS, ZPLXSZ
      COMPLEX*16         PADVAL
      PARAMETER          ( ZPLXSZ = 16,
     $                     MEMSIZ = TOTMEM / ZPLXSZ, NTESTS = 20,
     $                     PADVAL = ( -9923.0D+0, -9923.0D+0 ),
     $                     ZERO = 0.0D+0 )
      INTEGER            INT_ONE
      PARAMETER          ( INT_ONE = 1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            CHECK
      CHARACTER          UPLO
      CHARACTER*6        PASSED
      CHARACTER*80       OUTFILE
      INTEGER            BW, BW_NUM, FILLIN_SIZE, FREE_PTR, H, HH, I,
     $                   IAM, IASEED, IBSEED, ICTXT, ICTXTB, IERR_TEMP,
     $                   IMIDPAD, INFO, IPA, IPB, IPOSTPAD, IPREPAD,
     $                   IPW, IPW_SIZE, IPW_SOLVE, IPW_SOLVE_SIZE,
     $                   IP_DRIVER_W, IP_FILLIN, J, K, KFAIL, KPASS,
     $                   KSKIP, KTESTS, MYCOL, MYRHS_SIZE, MYROW, N, NB,
     $                   NBW, NGRIDS, NMAT, NNB, NNBR, NNR, NOUT, NP,
     $                   NPCOL, NPROCS, NPROCS_REAL, NPROW, NQ, NRHS,
     $                   N_FIRST, N_LAST, WORKSIZ
      REAL               THRESH
            DOUBLE PRECISION    ANORM, NOPS, NOPS2, SRESID, TMFLOPS,
     $                          TMFLOPS2
*     ..
*     .. Local Arrays ..
      INTEGER            BWVAL( NTESTS ), DESCA( 7 ), DESCA2D( DLEN_ ),
     $                   DESCB( 7 ), DESCB2D( DLEN_ ), IERR( 1 ),
     $                   NBRVAL( NTESTS ), NBVAL( NTESTS ),
     $                   NRVAL( NTESTS ), NVAL( NTESTS ),
     $                   PVAL( NTESTS ), QVAL( NTESTS )
      DOUBLE PRECISION   CTIME( 2 ), WTIME( 2 )
      COMPLEX*16         MEM( MEMSIZ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_EXIT, BLACS_GET,
     $                   BLACS_GRIDEXIT, BLACS_GRIDINFO, BLACS_GRIDINIT,
     $                   BLACS_PINFO, DESCINIT, IGSUM2D, PZBMATGEN,
     $                   PZCHEKPAD, PZFILLPAD, PZMATGEN, PZPBINFO,
     $                   PZPBLASCHK, PZPBTRF, PZPBTRS, SLBOOT,
     $                   SLCOMBINE, SLTIMER
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      LOGICAL            LSAME
      DOUBLE PRECISION   PZLANGE
      EXTERNAL           LSAME, NUMROC, PZLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, MOD
*     ..
*     .. Data Statements ..
      DATA               KFAIL, KPASS, KSKIP, KTESTS / 4*0 /
*     ..
*
*
*
*     .. Executable Statements ..
*
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
      IASEED = 100
      IBSEED = 200
*
      CALL PZPBINFO( OUTFILE, NOUT, UPLO, NMAT, NVAL, NTESTS, NBW,
     $               BWVAL, NTESTS, NNB, NBVAL, NTESTS, NNR, NRVAL,
     $               NTESTS, NNBR, NBRVAL, NTESTS, NGRIDS, PVAL, NTESTS,
     $               QVAL, NTESTS, THRESH, MEM, IAM, NPROCS )
*
      CHECK = ( THRESH.GE.0.0D+0 )
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
      DO 60 I = 1, NGRIDS
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
*
*
*        Define transpose process grid
*
         CALL BLACS_GET( -1, 0, ICTXTB )
         CALL BLACS_GRIDINIT( ICTXTB, 'Column-major', NPCOL, NPROW )
*
*        Go to bottom of process grid loop if this case doesn't use my
*        process
*
         CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
         IF( MYROW.LT.0 .OR. MYCOL.LT.0 ) THEN
            GO TO 50
         ENDIF
*
         DO 40 J = 1, NMAT
*
           IERR( 1 ) = 0
*
           N = NVAL( J )
*
*          Make sure matrix information is correct
*
           IF( N.LT.1 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9999 ) 'MATRIX', 'N', N
               IERR( 1 ) = 1
           END IF
*
*          Check all processes for an error
*
           CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                    -1, 0 )
*
           IF( IERR( 1 ).GT.0 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'size'
               KSKIP = KSKIP + 1
               GO TO 40
           END IF
*
*
           DO 45 BW_NUM = 1, NBW
*
             IERR( 1 ) = 0
*
             BW = BWVAL( BW_NUM )
             IF( BW.LT.0 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9999 ) 'Band', 'bw', BW
               IERR( 1 ) = 1
             END IF
*
             IF( BW.GT.N-1 ) THEN
               IERR( 1 ) = 1
             END IF
*
*            Check all processes for an error
*
             CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                    -1, 0 )
*
             IF( IERR( 1 ).GT.0 ) THEN
               KSKIP = KSKIP + 1
               GO TO 45
             END IF
*
             DO 30 K = 1, NNB
*
               IERR( 1 ) = 0
*
               NB = NBVAL( K )
               IF( NB.LT.0 ) THEN
                  NB =( (N-(NPCOL-1)*BW-1)/NPCOL + 1 )
     $               + BW
                  NB = MAX( NB, 2*BW )
                  NB = MIN( N, NB )
               END IF
*
*              Make sure NB is legal
*
               IERR( 1 ) = 0
               IF( NB.LT.MIN( 2*BW, N ) ) THEN
                  IERR( 1 ) = 1
               ENDIF
*
*              Check all processes for an error
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1,
     $                       -1, 0 )
*
               IF( IERR( 1 ).GT.0 ) THEN
                  KSKIP = KSKIP + 1
                  GO TO 30
               END IF
*
*              Padding constants
*
               NP = NUMROC( (BW+1), (BW+1),
     $                      MYROW, 0, NPROW )
               NQ = NUMROC( N, NB, MYCOL, 0, NPCOL )
*
               IF( CHECK ) THEN
                  IPREPAD  = ((BW+1)+10)
                  IMIDPAD  = 10
                  IPOSTPAD = ((BW+1)+10)
               ELSE
                  IPREPAD  = 0
                  IMIDPAD  = 0
                  IPOSTPAD = 0
               END IF
*
*              Initialize the array descriptor for the matrix A
*
               CALL DESCINIT( DESCA2D, (BW+1), N,
     $                       (BW+1), NB, 0, 0,
     $                       ICTXT,((BW+1)+10), IERR( 1 ) )
*
*              Convert this to 1D descriptor
*
               DESCA( 1 ) = 501
               DESCA( 3 ) = N
               DESCA( 4 ) = NB
               DESCA( 5 ) = 0
               DESCA( 2 ) = ICTXT
               DESCA( 6 ) = ((BW+1)+10)
               DESCA( 7 ) = 0
*
               IERR_TEMP = IERR( 1 )
               IERR( 1 ) = 0
               IERR( 1 ) = MIN( IERR( 1 ), IERR_TEMP )
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
               FREE_PTR = 1
               IPB = 0
*
*              Save room for prepadding
               FREE_PTR = FREE_PTR + IPREPAD
*
               IPA = FREE_PTR
               FREE_PTR = FREE_PTR + DESCA2D( LLD_ )*
     $                     DESCA2D( NB_ )
     $                     + IPOSTPAD
*
*              Add memory for fillin
*              Fillin space needs to store:
*                Fillin spike:
*                Contribution to previous proc's diagonal block of
*                  reduced system:
*                Off-diagonal block of reduced system:
*                Diagonal block of reduced system:
*
               FILLIN_SIZE =
     $            (NB+2*BW)*BW
*
*              Claim memory for fillin
*
               FREE_PTR = FREE_PTR + IPREPAD
               IP_FILLIN = FREE_PTR
               FREE_PTR = FREE_PTR + FILLIN_SIZE
*
*              Workspace needed by computational routines:
*
               IPW_SIZE = 0
*
*              factorization:
*
               IPW_SIZE = BW*BW
*
*              Claim memory for IPW
*
               IPW = FREE_PTR
               FREE_PTR = FREE_PTR + IPW_SIZE
*
*              Check for adequate memory for problem size
*
               IERR( 1 ) = 0
               IF( FREE_PTR.GT.MEMSIZ ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9996 )
     $               'divide and conquer factorization',
     $               (FREE_PTR )*ZPLXSZ
                  IERR( 1 ) = 1
               END IF
*
*              Check all processes for an error
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR,
     $                       1, -1, 0 )
*
               IF( IERR( 1 ).GT.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 ) 'MEMORY'
                  KSKIP = KSKIP + 1
                  GO TO 30
               END IF
*
*              Worksize needed for LAPRNT
               WORKSIZ = MAX( ((BW+1)+10), NB )
*
               IF( CHECK ) THEN
*
*                 Calculate the amount of workspace required by
*                 the checking routines.
*
*                 PZLANGE
                  WORKSIZ = MAX( WORKSIZ, DESCA2D( NB_ ) )
*
*                 PZPBLASCHK
                  WORKSIZ = MAX( WORKSIZ,
     $                   MAX(5,MAX(BW*(BW+2),NB))+2*NB )
               END IF
*
               FREE_PTR = FREE_PTR + IPREPAD
               IP_DRIVER_W = FREE_PTR
               FREE_PTR = FREE_PTR + WORKSIZ + IPOSTPAD
*
*
*              Check for adequate memory for problem size
*
               IERR( 1 ) = 0
               IF( FREE_PTR.GT.MEMSIZ ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9996 ) 'factorization',
     $               ( FREE_PTR )*ZPLXSZ
                  IERR( 1 ) = 1
               END IF
*
*              Check all processes for an error
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR,
     $                       1, -1, 0 )
*
               IF( IERR( 1 ).GT.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 ) 'MEMORY'
                  KSKIP = KSKIP + 1
                  GO TO 30
               END IF
*
               CALL PZBMATGEN( ICTXT, UPLO, 'B', BW, BW, N, (BW+1), NB,
     $                         MEM( IPA ), ((BW+1)+10), 0, 0, IASEED,
     $                         MYROW, MYCOL, NPROW, NPCOL )
*
               CALL PZFILLPAD( ICTXT, NP, NQ, MEM( IPA-IPREPAD ),
     $                          ((BW+1)+10), IPREPAD, IPOSTPAD,
     $                          PADVAL )
*
               CALL PZFILLPAD( ICTXT, WORKSIZ, 1,
     $                          MEM( IP_DRIVER_W-IPREPAD ), WORKSIZ,
     $                          IPREPAD, IPOSTPAD, PADVAL )
*
*              Calculate norm of A for residual error-checking
*
               IF( CHECK ) THEN
*
                  ANORM = PZLANGE( '1', (BW+1),
     $                            N, MEM( IPA ), 1, 1,
     $                            DESCA2D, MEM( IP_DRIVER_W ) )
                  CALL PZCHEKPAD( ICTXT, 'PZLANGE', NP, NQ,
     $                         MEM( IPA-IPREPAD ), ((BW+1)+10),
     $                         IPREPAD, IPOSTPAD, PADVAL )
                  CALL PZCHEKPAD( ICTXT, 'PZLANGE',
     $                            WORKSIZ, 1,
     $                            MEM( IP_DRIVER_W-IPREPAD ), WORKSIZ,
     $                            IPREPAD, IPOSTPAD, PADVAL )
               END IF
*
*
               CALL SLBOOT()
               CALL BLACS_BARRIER( ICTXT, 'All' )
*
*              Perform factorization
*
               CALL SLTIMER( 1 )
*
               CALL PZPBTRF( UPLO, N, BW, MEM( IPA ), 1, DESCA,
     $                       MEM( IP_FILLIN ), FILLIN_SIZE, MEM( IPW ),
     $                       IPW_SIZE, INFO )
*
               CALL SLTIMER( 1 )
*
               IF( INFO.NE.0 ) THEN
                  IF( IAM.EQ.0 ) THEN
                    WRITE( NOUT, FMT = * ) 'PZPBTRF INFO=', INFO
                  ENDIF
                  KFAIL = KFAIL + 1
                  GO TO 30
               END IF
*
               IF( CHECK ) THEN
*
*                 Check for memory overwrite in factorization
*
                  CALL PZCHEKPAD( ICTXT, 'PZPBTRF', NP,
     $                    NQ, MEM( IPA-IPREPAD ), ((BW+1)+10),
     $                    IPREPAD, IPOSTPAD, PADVAL )
               END IF
*
*
*              Loop over the different values for NRHS
*
               DO 20 HH = 1, NNR
*
                  IERR( 1 ) = 0
*
                  NRHS = NRVAL( HH )
*
*                    Initialize Array Descriptor for rhs
*
                     CALL DESCINIT( DESCB2D, N, NRHS, NB, 1, 0, 0,
     $                             ICTXTB, NB+10, IERR( 1 ) )
*
*                    Convert this to 1D descriptor
*
                     DESCB( 1 ) = 502
                     DESCB( 3 ) = N
                     DESCB( 4 ) = NB
                     DESCB( 5 ) = 0
                     DESCB( 2 ) = ICTXT
                     DESCB( 6 ) = DESCB2D( LLD_ )
                     DESCB( 7 ) = 0
*
*                    reset free_ptr to reuse space for right hand sides
*
                     IF( IPB .GT. 0 ) THEN
                       FREE_PTR = IPB
                     ENDIF
*
                     FREE_PTR = FREE_PTR + IPREPAD
                     IPB = FREE_PTR
                     FREE_PTR = FREE_PTR + NRHS*DESCB2D( LLD_ )
     $                          + IPOSTPAD
*
*                    Allocate workspace for workspace in TRS routine:
*
                     IPW_SOLVE_SIZE = (BW*NRHS)
*
                     IPW_SOLVE = FREE_PTR
                     FREE_PTR = FREE_PTR + IPW_SOLVE_SIZE
*
                     IERR( 1 ) = 0
                     IF( FREE_PTR.GT.MEMSIZ ) THEN
                        IF( IAM.EQ.0 )
     $                     WRITE( NOUT, FMT = 9996 )'solve',
     $                            ( FREE_PTR )*ZPLXSZ
                        IERR( 1 ) = 1
                     END IF
*
*                    Check all processes for an error
*
                     CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1,
     $                             IERR, 1, -1, 0 )
*
                     IF( IERR( 1 ).GT.0 ) THEN
                        IF( IAM.EQ.0 )
     $                     WRITE( NOUT, FMT = 9997 ) 'MEMORY'
                        KSKIP = KSKIP + 1
                        GO TO 15
                     END IF
*
                     MYRHS_SIZE = NUMROC( N, NB, MYCOL, 0, NPCOL )
*
*                    Generate RHS
*
                     CALL PZMATGEN(ICTXTB, 'No', 'No',
     $                        DESCB2D( M_ ), DESCB2D( N_ ),
     $                        DESCB2D( MB_ ), DESCB2D( NB_ ),
     $                        MEM( IPB ),
     $                        DESCB2D( LLD_ ), DESCB2D( RSRC_ ),
     $                        DESCB2D( CSRC_ ),
     $                        IBSEED, 0, MYRHS_SIZE, 0, NRHS, MYCOL,
     $                        MYROW, NPCOL, NPROW )
*
                     IF( CHECK ) THEN
                        CALL PZFILLPAD( ICTXTB, NB, NRHS,
     $                                  MEM( IPB-IPREPAD ),
     $                                  DESCB2D( LLD_ ),
     $                                  IPREPAD, IPOSTPAD,
     $                                  PADVAL )
                        CALL PZFILLPAD( ICTXT, WORKSIZ, 1,
     $                                  MEM( IP_DRIVER_W-IPREPAD ),
     $                                  WORKSIZ, IPREPAD,
     $                                  IPOSTPAD, PADVAL )
                     END IF
*
*
                     CALL BLACS_BARRIER( ICTXT, 'All')
                     CALL SLTIMER( 2 )
*
*                    Solve linear system via factorization
*
                     CALL PZPBTRS( UPLO, N, BW, NRHS, MEM( IPA ), 1,
     $                             DESCA, MEM( IPB ), 1, DESCB,
     $                             MEM( IP_FILLIN ), FILLIN_SIZE,
     $                             MEM( IPW_SOLVE ), IPW_SOLVE_SIZE,
     $                             INFO )
*
                     CALL SLTIMER( 2 )
*
                     IF( INFO.NE.0 ) THEN
                       IF( IAM.EQ.0 )
     $  WRITE( NOUT, FMT = * ) 'PZPBTRS INFO=', INFO
                       KFAIL = KFAIL + 1
                       PASSED = 'FAILED'
                       GO TO 20
                     END IF
*
                     IF( CHECK ) THEN
*
*                       check for memory overwrite
*
                        CALL PZCHEKPAD( ICTXT, 'PZPBTRS-work',
     $                                  WORKSIZ, 1,
     $                                  MEM( IP_DRIVER_W-IPREPAD ),
     $                                  WORKSIZ, IPREPAD,
     $                                  IPOSTPAD, PADVAL )
*
*                       check the solution to rhs
*
                        SRESID = ZERO
*
                        CALL PZPBLASCHK( 'H', UPLO, N, BW, BW, NRHS,
     $                              MEM( IPB ), 1, 1, DESCB2D,
     $                              IASEED, MEM( IPA ), 1, 1, DESCA2D,
     $                              IBSEED, ANORM, SRESID,
     $                              MEM( IP_DRIVER_W ), WORKSIZ )
*
                        IF( IAM.EQ.0 ) THEN
                           IF( SRESID.GT.THRESH )
     $                        WRITE( NOUT, FMT = 9985 ) SRESID
                        END IF
*
*                       The second test is a NaN trap
*
                        IF( ( SRESID.LE.THRESH          ).AND.
     $                      ( (SRESID-SRESID).EQ.0.0D+0 ) ) THEN
                           KPASS = KPASS + 1
                           PASSED = 'PASSED'
                        ELSE
                           KFAIL = KFAIL + 1
                           PASSED = 'FAILED'
                        END IF
*
                     END IF
*
   15                CONTINUE
*                    Skipped tests jump to here to print out "SKIPPED"
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
                        NOPS = 0
                        NOPS2 = 0
*
                        N_FIRST = NB
                        NPROCS_REAL = ( N-1 )/NB + 1
                        N_LAST = MOD( N-1, NB ) + 1
*
*
                        NOPS = NOPS + DBLE(BW)*( -2.D0 / 3.D0+DBLE(BW)*
     $                        ( -1.D0+DBLE(BW)*( -1.D0 / 3.D0 ) ) ) +
     $                        DBLE(N)*( 1.D0+DBLE(BW)*( 3.D0 /
     $                        2.D0+DBLE(BW)*( 1.D0 / 2.D0 ) ) )
                        NOPS = NOPS + DBLE(BW)*( -1.D0 / 6.D0+DBLE(BW)
     $                        *( -1.D0 /2.D0+DBLE(BW)
     $                        *( -1.D0 / 3.D0 ) ) ) +
     $                        DBLE(N)*( DBLE(BW) /
     $                        2.D0*( 1.D0+DBLE(BW) ) )
*
                        NOPS = NOPS +
     $                         DBLE(NRHS)*( ( 2*DBLE(N)-DBLE(BW) )*
     $                         ( DBLE(BW)+1.D0 ) )+ DBLE(NRHS)*
     $                         ( DBLE(BW)*( 2*DBLE(N)-
     $                         ( DBLE(BW)+1.D0 ) ) )
*
*
*                       Second calc to represent actual hardware speed
*
*                     NB bw^2  flops for LLt factorization in 1st proc
*
                      NOPS2 = ( (DBLE(N_FIRST))* DBLE(BW)**2  )
*
                      IF ( NPROCS_REAL .GT. 1) THEN
*                       4 NB bw^2  flops for LLt factorization and
*                         spike calc in last processor
*
                        NOPS2 = NOPS2 +
     $                          4*( (DBLE(N_LAST)*DBLE(BW)**2) )
                      ENDIF
*
                      IF ( NPROCS_REAL .GT. 2) THEN
*                       4 NB bw^2  flops for LLt factorization and
*                         spike calc in other processors
*
                        NOPS2 = NOPS2 + (NPROCS_REAL-2)*
     $                          4*( (DBLE(NB)*DBLE(BW)**2) )
                      ENDIF
*
*                     Reduced system
*
                      NOPS2 = NOPS2 +
     $                  ( NPROCS_REAL-1 ) * ( BW*BW*BW/3 )
                      IF( NPROCS_REAL .GT. 1 ) THEN
                        NOPS2 = NOPS2 +
     $                     ( NPROCS_REAL-2 ) * ( 2 * BW*BW*BW )
                      ENDIF
*
*
*                     nrhs * 4 n_first*bw flops for LLt solve in proc 1.
*
                      NOPS2 = NOPS2 +
     $                    ( 4.0D+0*(DBLE(N_FIRST)*DBLE(BW))*DBLE(NRHS) )
*
                      IF ( NPROCS_REAL .GT. 1 ) THEN
*
*                     2*nrhs*4 n_last*bw flops for LLt solve in last.
*
                        NOPS2 = NOPS2 +
     $                  2*( 4.0D+0*(DBLE(N_LAST)*DBLE(BW))*DBLE(NRHS) )
                      ENDIF
*
                      IF ( NPROCS_REAL .GT. 2 ) THEN
*
*                     2 * nrhs * 4 NB*bw flops for LLt solve in others.
*
                        NOPS2 = NOPS2 +
     $                    ( NPROCS_REAL-2)*2*
     $                    ( 4.0D+0*(DBLE(NB)*DBLE(BW))*DBLE(NRHS) )
                      ENDIF
*
*                     Reduced system
*
                      NOPS2 = NOPS2 +
     $                  NRHS*( NPROCS_REAL-1 ) * ( BW*BW )
                      IF( NPROCS_REAL .GT. 1 ) THEN
                        NOPS2 = NOPS2 +
     $                   NRHS*( NPROCS_REAL-2 ) * ( 3 * BW*BW )
                      ENDIF
*
*
*                     Multiply by 4 to get complex count
*
                      NOPS2 = NOPS2 * DBLE(4)
*
*                       Calculate total megaflops - factorization and/or
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
                        IF( WTIME( 1 )+WTIME( 2 ).GT.0.0D+0 ) THEN
                           TMFLOPS2 = NOPS2 /
     $                            ( ( WTIME( 1 )+WTIME( 2 ) ) * 1.0D+6 )
                        ELSE
                           TMFLOPS2 = 0.0D+0
                        END IF
*
                        IF( WTIME( 2 ).GE.0.0D+0 )
     $                     WRITE( NOUT, FMT = 9993 ) 'WALL', UPLO,
     $                            N,
     $                            BW,
     $                            NB, NRHS, NPROW, NPCOL,
     $                            WTIME( 1 ), WTIME( 2 ), TMFLOPS,
     $                            TMFLOPS2, PASSED
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
                        IF( CTIME( 1 )+CTIME( 2 ).GT.0.0D+0 ) THEN
                           TMFLOPS2 = NOPS2 /
     $                            ( ( CTIME( 1 )+CTIME( 2 ) ) * 1.0D+6 )
                        ELSE
                           TMFLOPS2 = 0.0D+0
                        END IF
*
                        IF( CTIME( 2 ).GE.0.0D+0 )
     $                     WRITE( NOUT, FMT = 9993 ) 'CPU ', UPLO,
     $                            N,
     $                            BW,
     $                            NB, NRHS, NPROW, NPCOL,
     $                            CTIME( 1 ), CTIME( 2 ), TMFLOPS,
     $                            TMFLOPS2, PASSED
*
                     END IF
   20          CONTINUE
*
*
   30       CONTINUE
*           NNB loop
*
   45      CONTINUE
*          BW[] loop
*
   40   CONTINUE
*       NMAT loop
*
        CALL BLACS_GRIDEXIT( ICTXT )
        CALL BLACS_GRIDEXIT( ICTXTB )
*
   50   CONTINUE
*       NGRIDS DROPOUT
   60 CONTINUE
*     NGRIDS loop
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
 9995 FORMAT( 'TIME UL      N  BW   NB  NRHS  P    Q L*U Time ',
     $        'Slv Time   MFLOPS   MFLOP2  CHECK' )
 9994 FORMAT( '---- -- ------ --- ---- ----- -- ---- -------- ',
     $        '--------   ------   ------ ------' )
 9993 FORMAT( A4, 2X, A1, 1X, I6, 1X, I3, 1X, I4, 1X,
     $        I5, 1X, I2, 1X,
     $        I4, 1X, F8.3, F9.4, F9.2, F9.2, 1X, A6 )
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
*     End of PZPBTRS_DRIVER
*
      END
*
