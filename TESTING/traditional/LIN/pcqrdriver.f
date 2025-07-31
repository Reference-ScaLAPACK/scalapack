      PROGRAM PCQRDRIVER
*
*  -- ScaLAPACK testing driver (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 28, 2001
*
*  Purpose
*  =======
*
*  PCQRDRIVER is the main test program for the COMPLEX
*  SCALAPACK QR factorization routines. This test driver performs a QR
*  QL, LQ, RQ, QP (QR factorization with column pivoting) or TZ
*  (complete unitary factorization) factorization and checks the
*  results.
*
*  The program must be driven by a short data file.  An annotated
*  example of a data file can be obtained by deleting the first 3
*  characters from the following 16 lines:
*  'ScaLAPACK QR factorizations input file'
*  'PVM machine'
*  'QR.out'                      output file name (if any)
*  6                             device out
*  6                             number of factorizations
*  'QR' 'QL' 'LQ' 'RQ' 'QP' 'TZ' factorization: QR, QL, LQ, RQ, QP, TZ
*  4                             number of problems sizes
*  55 17 31 201                  values of M
*  5 71 31 201                   values of N
*  3                             number of MB's and NB's
*  4 3 5                         values of MB
*  4 7 3                         values of NB
*  7                             number of process grids (ordered P & Q)
*  1 2 1 4 2 3 8                 values of P
*  7 2 4 1 3 2 1                 values of Q
*  1.0                           threshold
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
*  CPLXSZ   INTEGER, default = 8 bytes.
*           INTGSZ, REALSZ and CPLXSZ indicate the length in bytes on
*           the given platform for an integer, a single precision real
*           and a single precision complex.
*  MEM      COMPLEX array, dimension ( TOTMEM / CPLXSZ )
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
      INTEGER            CPLXSZ, INTGSZ, MEMSIZ, NTESTS, REALSZ, TOTMEM
      COMPLEX            PADVAL
      PARAMETER          ( CPLXSZ = 8, INTGSZ = 4, REALSZ = 4,
     $                     TOTMEM = 2000000, MEMSIZ = TOTMEM / CPLXSZ,
     $                     NTESTS = 20,
     $                     PADVAL = ( -9923.0E+0, -9923.0E+0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        FACT
      CHARACTER*6        PASSED
      CHARACTER*7        ROUT
      CHARACTER*8        ROUTCHK
      CHARACTER*80       OUTFILE
      LOGICAL            CHECK
      INTEGER            I, IAM, IASEED, ICTXT, IMIDPAD, INFO, IPA,
     $                   IPOSTPAD, IPPIV, IPREPAD, IPTAU, IPRW, IPW, J,
     $                   K, KFAIL, KPASS, KSKIP, KTESTS, L, LIPIV,
     $                   LRWORK, LTAU, LWORK, M, MAXMN, MB, MINMN, MNP,
     $                   MNQ, MP, MYCOL, MYROW, N, NB, NFACT, NGRIDS,
     $                   NMAT, NNB, NOUT, NPCOL, NPROCS, NPROW, NQ,
     $                   WORKFCT, WORKRFCT, WORKSIZ
      REAL               ANORM, FRESID, THRESH
      DOUBLE PRECISION   NOPS, TMFLOPS
*     ..
*     .. Arrays ..
      CHARACTER*2        FACTOR( NTESTS )
      INTEGER            DESCA( DLEN_ ), IERR( 1 ), MBVAL( NTESTS ),
     $                   MVAL( NTESTS ), NBVAL( NTESTS ),
     $                   NVAL( NTESTS ), PVAL( NTESTS ), QVAL( NTESTS )
      DOUBLE PRECISION   CTIME( 1 ), WTIME( 1 )
      COMPLEX            MEM( MEMSIZ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_EXIT, BLACS_GET,
     $                   BLACS_GRIDEXIT, BLACS_GRIDINFO, BLACS_GRIDINIT,
     $                   BLACS_PINFO, DESCINIT, IGSUM2D, PCCHEKPAD,
     $                   PCFILLPAD, PCGELQF, PCGELQRV,
     $                   PCGEQLF, PCGEQLRV, PCGEQPF,
     $                   PCQPPIV, PCGEQRF, PCGEQRRV,
     $                   PCGERQF, PCGERQRV, PCTZRZRV,
     $                   PCMATGEN, PCLAFCHK, PCQRINFO,
     $                   PCTZRZF, SLBOOT, SLCOMBINE, SLTIMER
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      INTEGER            ICEIL, NUMROC
      REAL               PCLANGE
      EXTERNAL           ICEIL, LSAMEN, NUMROC, PCLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Data Statements ..
      DATA               KTESTS, KPASS, KFAIL, KSKIP /4*0/
*     ..
*     .. Executable Statements ..
*
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
      IASEED = 100
      CALL PCQRINFO( OUTFILE, NOUT, NFACT, FACTOR, NTESTS, NMAT, MVAL,
     $               NTESTS, NVAL, NTESTS, NNB, MBVAL, NTESTS, NBVAL,
     $               NTESTS, NGRIDS, PVAL, NTESTS, QVAL, NTESTS,
     $               THRESH, MEM, IAM, NPROCS )
      CHECK = ( THRESH.GE.0.0E+0 )
*
*     Loop over the different factorization types
*
      DO 40 I = 1, NFACT
*
         FACT = FACTOR( I )
*
*        Print headings
*
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = * )
            IF( LSAMEN( 2, FACT, 'QR' ) ) THEN
               ROUT = 'PCGEQRF'
               ROUTCHK = 'PCGEQRRV'
               WRITE( NOUT, FMT = 9986 )
     $                'QR factorization tests.'
            ELSE IF( LSAMEN( 2, FACT, 'QL' ) ) THEN
               ROUT = 'PCGEQLF'
               ROUTCHK = 'PCGEQLRV'
               WRITE( NOUT, FMT = 9986 )
     $                'QL factorization tests.'
            ELSE IF( LSAMEN( 2, FACT, 'LQ' ) ) THEN
               ROUT = 'PCGELQF'
               ROUTCHK = 'PCGELQRV'
               WRITE( NOUT, FMT = 9986 )
     $                'LQ factorization tests.'
            ELSE IF( LSAMEN( 2, FACT, 'RQ' ) ) THEN
               ROUT = 'PCGERQF'
               ROUTCHK = 'PCGERQRV'
               WRITE( NOUT, FMT = 9986 )
     $                'RQ factorization tests.'
            ELSE IF( LSAMEN( 2, FACT, 'QP' ) ) THEN
               ROUT = 'PCGEQPF'
               ROUTCHK = 'PCGEQRRV'
               WRITE( NOUT, FMT = 9986 )
     $                'QR factorization with column pivoting tests.'
            ELSE IF( LSAMEN( 2, FACT, 'TZ' ) ) THEN
               ROUT = 'PCTZRZF'
               ROUTCHK = 'PCTZRZRV'
               WRITE( NOUT, FMT = 9986 )
     $                'Complete unitary factorization tests.'
            END IF
            WRITE( NOUT, FMT = * )
            WRITE( NOUT, FMT = 9995 )
            WRITE( NOUT, FMT = 9994 )
            WRITE( NOUT, FMT = * )
         END IF
*
*        Loop over different process grids
*
         DO 30 J = 1, NGRIDS
*
            NPROW = PVAL( J )
            NPCOL = QVAL( J )
*
*           Make sure grid information is correct
*
            IERR( 1 ) = 0
            IF( NPROW.LT.1 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9999 ) 'GRID', 'nprow', NPROW
               IERR( 1 ) = 1
            ELSE IF( NPCOL.LT.1 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9999 ) 'GRID', 'npcol', NPCOL
               IERR( 1 ) = 1
            ELSE IF( NPROW*NPCOL.GT.NPROCS ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9998 ) NPROW*NPCOL, NPROCS
               IERR( 1 ) = 1
            END IF
*
            IF( IERR( 1 ).GT.0 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'grid'
               KSKIP = KSKIP + 1
               GO TO 30
            END IF
*
*           Define process grid
*
            CALL BLACS_GET( -1, 0, ICTXT )
            CALL BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
            CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*           Go to bottom of loop if this case doesn't use my process
*
            IF( MYROW.GE.NPROW .OR. MYCOL.GE.NPCOL )
     $         GO TO 30
*
            DO 20 K = 1, NMAT
*
               M = MVAL( K )
               N = NVAL( K )
*
*              Make sure matrix information is correct
*
               IERR(1) = 0
               IF( M.LT.1 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9999 ) 'MATRIX', 'M', M
                  IERR( 1 ) = 1
               ELSE IF( N.LT.1 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9999 ) 'MATRIX', 'N', N
                  IERR( 1 ) = 1
               END IF
*
*              Make sure no one had error
*
               CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
               IF( IERR( 1 ).GT.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9997 ) 'matrix'
                  KSKIP = KSKIP + 1
                  GO TO 20
               END IF
*
*              Loop over different blocking sizes
*
               DO 10 L = 1, NNB
*
                  MB = MBVAL( L )
                  NB = NBVAL( L )
*
*                 Make sure mb is legal
*
                  IERR( 1 ) = 0
                  IF( MB.LT.1 ) THEN
                     IERR( 1 ) = 1
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9999 ) 'MB', 'MB', MB
                  END IF
*
*                 Check all processes for an error
*
                  CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                          0 )
*
                  IF( IERR( 1 ).GT.0 ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9997 ) 'MB'
                     KSKIP = KSKIP + 1
                     GO TO 10
                  END IF
*
*                 Make sure nb is legal
*
                  IERR( 1 ) = 0
                  IF( NB.LT.1 ) THEN
                     IERR( 1 ) = 1
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9999 ) 'NB', 'NB', NB
                  END IF
*
*                 Check all processes for an error
*
                  CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                          0 )
*
                  IF( IERR( 1 ).GT.0 ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9997 ) 'NB'
                     KSKIP = KSKIP + 1
                     GO TO 10
                  END IF
*
*                 Padding constants
*
                  MP  = NUMROC( M, MB, MYROW, 0, NPROW )
                  NQ  = NUMROC( N, NB, MYCOL, 0, NPCOL )
                  MNP = NUMROC( MIN( M, N ), MB, MYROW, 0, NPROW )
                  MNQ = NUMROC( MIN( M, N ), NB, MYCOL, 0, NPCOL )
                  IF( CHECK ) THEN
                     IPREPAD  = MAX( MB, MP )
                     IMIDPAD  = NB
                     IPOSTPAD = MAX( NB, NQ )
                  ELSE
                     IPREPAD  = 0
                     IMIDPAD  = 0
                     IPOSTPAD = 0
                  END IF
*
*                 Initialize the array descriptor for the matrix A
*
                  CALL DESCINIT( DESCA, M, N, MB, NB, 0, 0, ICTXT,
     $                           MAX( 1, MP ) + IMIDPAD, IERR( 1 ) )
*
*                 Check all processes for an error
*
                  CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1,
     $                          0 )
*
                  IF( IERR( 1 ).LT.0 ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9997 ) 'descriptor'
                     KSKIP = KSKIP + 1
                     GO TO 10
                  END IF
*
*                 Assign pointers into MEM for ScaLAPACK arrays, A is
*                 allocated starting at position MEM( IPREPAD+1 )
*
                  IPA   = IPREPAD+1
                  IPTAU = IPA + DESCA( LLD_ ) * NQ + IPOSTPAD + IPREPAD
*
                  IF( LSAMEN( 2, FACT, 'QR' ) ) THEN
*
                     LTAU = MNQ
                     IPW  = IPTAU + LTAU + IPOSTPAD + IPREPAD
*
*                    Figure the amount of workspace required by the QR
*                    factorization
*
                     LWORK = DESCA( NB_ ) * ( MP + NQ + DESCA( NB_ ) )
                     WORKFCT = LWORK + IPOSTPAD
                     WORKSIZ = WORKFCT
*
                     IF( CHECK ) THEN
*
*                       Figure the amount of workspace required by the
*                       checking routines PCLAFCHK, PCGEQRRV and
*                       PCLANGE
*
                        WORKSIZ = LWORK + MP*DESCA( NB_ ) + IPOSTPAD
*
                     END IF
*
                  ELSE IF( LSAMEN( 2, FACT, 'QL' ) ) THEN
*
                     LTAU = NQ
                     IPW = IPTAU + LTAU + IPOSTPAD + IPREPAD
*
*                    Figure the amount of workspace required by the QL
*                    factorization
*
                     LWORK = DESCA( NB_ ) * ( MP + NQ + DESCA( NB_ ) )
                     WORKFCT = LWORK + IPOSTPAD
                     WORKSIZ = WORKFCT
*
                     IF( CHECK ) THEN
*
*                       Figure the amount of workspace required by the
*                       checking routines PCLAFCHK, PCGEQLRV and
*                       PCLANGE
*
                        WORKSIZ = LWORK + MP*DESCA( NB_ ) + IPOSTPAD
*
                     END IF
*
                  ELSE IF( LSAMEN( 2, FACT, 'LQ' ) ) THEN
*
                     LTAU = MNP
                     IPW = IPTAU + LTAU + IPOSTPAD + IPREPAD
*
*                    Figure the amount of workspace required by the LQ
*                    factorization
*
                     LWORK = DESCA( MB_ ) * ( MP + NQ + DESCA( MB_ ) )
                     WORKFCT = LWORK + IPOSTPAD
                     WORKSIZ = WORKFCT
*
                     IF( CHECK ) THEN
*
*                       Figure the amount of workspace required by the
*                       checking routines PCLAFCHK, PCGELQRV and
*                       PCLANGE
*
                        WORKSIZ = LWORK +
     $                            MAX( MP*DESCA( NB_ ), NQ*DESCA( MB_ )
     $                            ) + IPOSTPAD
*
                     END IF
*
                  ELSE IF( LSAMEN( 2, FACT, 'RQ' ) ) THEN
*
                     LTAU = MP
                     IPW = IPTAU + LTAU + IPOSTPAD + IPREPAD
*
*                    Figure the amount of workspace required by the QR
*                    factorization
*
                     LWORK = DESCA( MB_ ) * ( MP + NQ + DESCA( MB_ ) )
                     WORKFCT = LWORK + IPOSTPAD
                     WORKSIZ = WORKFCT
*
                     IF( CHECK ) THEN
*
*                       Figure the amount of workspace required by the
*                       checking routines PCLAFCHK, PCGERQRV and
*                       PCLANGE
*
                        WORKSIZ = LWORK +
     $                            MAX( MP*DESCA( NB_ ), NQ*DESCA( MB_ )
     $                            ) + IPOSTPAD
*
                     END IF
*
                  ELSE IF( LSAMEN( 2, FACT, 'QP' ) ) THEN
*
                     LTAU = MNQ
                     IPPIV = IPTAU + LTAU + IPOSTPAD + IPREPAD
                     LIPIV = ICEIL( INTGSZ*NQ, CPLXSZ )
                     IPW = IPPIV + LIPIV + IPOSTPAD + IPREPAD
*
*                    Figure the amount of workspace required by the
*                    factorization i.e from IPW on.
*
                     LWORK = MAX( 3, MP + MAX( 1, NQ ) )
                     WORKFCT = LWORK + IPOSTPAD
                     LRWORK = MAX( 1, 2 * NQ )
                     WORKRFCT = ICEIL( LRWORK*REALSZ, CPLXSZ ) +
     $                          IPOSTPAD
                     IPRW = IPW + WORKFCT + IPREPAD
                     WORKSIZ = WORKFCT + IPREPAD + WORKRFCT
*
                     IF( CHECK ) THEN
*
*                       Figure the amount of workspace required by the
*                       checking routines PCLAFCHK, PCGEQRRV,
*                       PCLANGE.
*
                        WORKSIZ = MAX( WORKSIZ - IPOSTPAD,
     $                    DESCA( NB_ )*( 2*MP + NQ + DESCA( NB_ ) ) ) +
     $                    IPOSTPAD
                     END IF
*
                  ELSE IF( LSAMEN( 2, FACT, 'TZ' ) ) THEN
*
                     LTAU = MP
                     IPW = IPTAU + LTAU + IPOSTPAD + IPREPAD
*
*                    Figure the amount of workspace required by the TZ
*                    factorization
*
                     LWORK = DESCA( MB_ ) * ( MP + NQ + DESCA( MB_ ) )
                     WORKFCT = LWORK + IPOSTPAD
                     WORKSIZ = WORKFCT
*
                     IF( CHECK ) THEN
*
*                       Figure the amount of workspace required by the
*                       checking routines PCLAFCHK, PCTZRZRV and
*                       PCLANGE
*
                        WORKSIZ = LWORK +
     $                            MAX( MP*DESCA( NB_ ), NQ*DESCA( MB_ )
     $                            ) + IPOSTPAD
*
                     END IF
*
                  END IF
*
*                 Check for adequate memory for problem size
*
                  IERR( 1 ) = 0
                  IF( IPW+WORKSIZ.GT.MEMSIZ ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9996 )
     $                         FACT // ' factorization',
     $                         ( IPW+WORKSIZ )*CPLXSZ
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
                     GO TO 10
                  END IF
*
*                 Generate the matrix A
*
                  CALL PCMATGEN( ICTXT, 'N', 'N', DESCA( M_ ),
     $                           DESCA( N_ ), DESCA( MB_ ),
     $                           DESCA( NB_ ), MEM( IPA ),
     $                           DESCA( LLD_ ), DESCA( RSRC_ ),
     $                           DESCA( CSRC_ ), IASEED, 0, MP, 0, NQ,
     $                           MYROW, MYCOL, NPROW, NPCOL )
*
*                 Need the Infinity of A for checking
*
                  IF( CHECK ) THEN
                     CALL PCFILLPAD( ICTXT, MP, NQ, MEM( IPA-IPREPAD ),
     $                               DESCA( LLD_ ), IPREPAD, IPOSTPAD,
     $                               PADVAL )
                     IF( LSAMEN( 2, FACT, 'QP' ) ) THEN
                        CALL PCFILLPAD( ICTXT, LIPIV, 1,
     $                                  MEM( IPPIV-IPREPAD ), LIPIV,
     $                                  IPREPAD, IPOSTPAD, PADVAL )
                     END IF
                     CALL PCFILLPAD( ICTXT, LTAU, 1,
     $                               MEM( IPTAU-IPREPAD ), LTAU,
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     CALL PCFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                               MEM( IPW-IPREPAD ),
     $                               WORKSIZ-IPOSTPAD,
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     ANORM = PCLANGE( 'I', M, N, MEM( IPA ), 1, 1,
     $                                DESCA, MEM( IPW ) )
                     CALL PCCHEKPAD( ICTXT, 'PCLANGE', MP, NQ,
     $                               MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     CALL PCCHEKPAD( ICTXT, 'PCLANGE',
     $                               WORKSIZ-IPOSTPAD, 1,
     $                               MEM( IPW-IPREPAD ),
     $                               WORKSIZ-IPOSTPAD, IPREPAD,
     $                               IPOSTPAD, PADVAL )
                     IF( LSAMEN( 2, FACT, 'QP' ) ) THEN
                        CALL PCFILLPAD( ICTXT, WORKRFCT-IPOSTPAD, 1,
     $                                  MEM( IPRW-IPREPAD ),
     $                                  WORKRFCT-IPOSTPAD,
     $                                  IPREPAD, IPOSTPAD, PADVAL )
                     END IF
                     CALL PCFILLPAD( ICTXT, WORKFCT-IPOSTPAD, 1,
     $                               MEM( IPW-IPREPAD ),
     $                               WORKFCT-IPOSTPAD,
     $                               IPREPAD, IPOSTPAD, PADVAL )
                  END IF
*
                  CALL SLBOOT()
                  CALL BLACS_BARRIER( ICTXT, 'All' )
*
*                 Perform QR factorizations
*
                  IF( LSAMEN( 2, FACT, 'QR' ) ) THEN
                     CALL SLTIMER( 1 )
                     CALL PCGEQRF( M, N, MEM( IPA ), 1, 1, DESCA,
     $                             MEM( IPTAU ), MEM( IPW ), LWORK,
     $                             INFO )
                     CALL SLTIMER( 1 )
                  ELSE IF( LSAMEN( 2, FACT, 'QL' ) ) THEN
                     CALL SLTIMER( 1 )
                     CALL PCGEQLF( M, N, MEM( IPA ), 1, 1, DESCA,
     $                             MEM( IPTAU ), MEM( IPW ), LWORK,
     $                             INFO )
                     CALL SLTIMER( 1 )
                  ELSE IF( LSAMEN( 2, FACT, 'LQ' ) ) THEN
                     CALL SLTIMER( 1 )
                     CALL PCGELQF( M, N, MEM( IPA ), 1, 1, DESCA,
     $                             MEM( IPTAU ), MEM( IPW ), LWORK,
     $                             INFO )
                     CALL SLTIMER( 1 )
                  ELSE IF( LSAMEN( 2, FACT, 'RQ' ) ) THEN
                     CALL SLTIMER( 1 )
                     CALL PCGERQF( M, N, MEM( IPA ), 1, 1, DESCA,
     $                             MEM( IPTAU ), MEM( IPW ), LWORK,
     $                             INFO )
                     CALL SLTIMER( 1 )
                  ELSE IF( LSAMEN( 2, FACT, 'QP' ) ) THEN
                     CALL SLTIMER( 1 )
                     CALL PCGEQPF( M, N, MEM( IPA ), 1, 1, DESCA,
     $                             MEM( IPPIV ), MEM( IPTAU ),
     $                             MEM( IPW ), LWORK, MEM( IPRW ),
     $                             LRWORK, INFO )
                     CALL SLTIMER( 1 )
                  ELSE IF( LSAMEN( 2, FACT, 'TZ' ) ) THEN
                     CALL SLTIMER( 1 )
                     IF( N.GE.M )
     $                  CALL PCTZRZF( M, N, MEM( IPA ), 1, 1, DESCA,
     $                                MEM( IPTAU ), MEM( IPW ), LWORK,
     $                                INFO )
                     CALL SLTIMER( 1 )
                  END IF
*
                  IF( CHECK ) THEN
*
*                    Check for memory overwrite in factorization
*
                     CALL PCCHEKPAD( ICTXT, ROUT, MP, NQ,
     $                               MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     CALL PCCHEKPAD( ICTXT, ROUT, LTAU, 1,
     $                               MEM( IPTAU-IPREPAD ), LTAU,
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     IF( LSAMEN( 2, FACT, 'QP' ) ) THEN
                        CALL PCCHEKPAD( ICTXT, ROUT, LIPIV, 1,
     $                                  MEM( IPPIV-IPREPAD ), LIPIV,
     $                                  IPREPAD, IPOSTPAD, PADVAL )
                        CALL PCCHEKPAD( ICTXT, ROUT, WORKRFCT-IPOSTPAD,
     $                                  1, MEM( IPRW-IPREPAD ),
     $                                  WORKRFCT-IPOSTPAD,
     $                                  IPREPAD, IPOSTPAD, PADVAL )
                     END IF
                     CALL PCCHEKPAD( ICTXT, ROUT, WORKFCT-IPOSTPAD, 1,
     $                               MEM( IPW-IPREPAD ),
     $                               WORKFCT-IPOSTPAD, IPREPAD,
     $                               IPOSTPAD, PADVAL )
                     CALL PCFILLPAD( ICTXT, WORKSIZ-IPOSTPAD, 1,
     $                               MEM( IPW-IPREPAD ),
     $                               WORKSIZ-IPOSTPAD,
     $                               IPREPAD, IPOSTPAD, PADVAL )
*
                     IF( LSAMEN( 2, FACT, 'QR' ) ) THEN
*
*                       Compute residual = ||A-Q*R|| / (||A||*N*eps)
*
                        CALL PCGEQRRV( M, N, MEM( IPA ), 1, 1, DESCA,
     $                                 MEM( IPTAU ), MEM( IPW ) )
                        CALL PCLAFCHK( 'No', 'No', M, N, MEM( IPA ), 1,
     $                              1, DESCA, IASEED, ANORM, FRESID,
     $                              MEM( IPW ) )
                     ELSE IF( LSAMEN( 2, FACT, 'QL' ) ) THEN
*
*                       Compute residual = ||A-Q*L|| / (||A||*N*eps)
*
                        CALL PCGEQLRV( M, N, MEM( IPA ), 1, 1, DESCA,
     $                                 MEM( IPTAU ), MEM( IPW ) )
                        CALL PCLAFCHK( 'No', 'No', M, N, MEM( IPA ), 1,
     $                              1, DESCA, IASEED, ANORM, FRESID,
     $                              MEM( IPW ) )
                     ELSE IF( LSAMEN( 2, FACT, 'LQ' ) ) THEN
*
*                       Compute residual = ||A-L*Q|| / (||A||*N*eps)
*
                        CALL PCGELQRV( M, N, MEM( IPA ), 1, 1, DESCA,
     $                                 MEM( IPTAU ), MEM( IPW ) )
                        CALL PCLAFCHK( 'No', 'No', M, N, MEM( IPA ), 1,
     $                              1, DESCA, IASEED, ANORM, FRESID,
     $                              MEM( IPW ) )
                     ELSE IF( LSAMEN( 2, FACT, 'RQ' ) ) THEN
*
*                       Compute residual = ||A-R*Q|| / (||A||*N*eps)
*
                        CALL PCGERQRV( M, N, MEM( IPA ), 1, 1, DESCA,
     $                                 MEM( IPTAU ), MEM( IPW ) )
                        CALL PCLAFCHK( 'No', 'No', M, N, MEM( IPA ), 1,
     $                              1, DESCA, IASEED, ANORM, FRESID,
     $                              MEM( IPW ) )
                     ELSE IF( LSAMEN( 2, FACT, 'QP' ) ) THEN
*
*                       Compute residual = ||AP-Q*R|| / (||A||*N*eps)
*
                        CALL PCGEQRRV( M, N, MEM( IPA ), 1, 1, DESCA,
     $                                 MEM( IPTAU ), MEM( IPW ) )
                     ELSE IF( LSAMEN( 2, FACT, 'TZ' ) ) THEN
*
*                       Compute residual = ||A-T*Z|| / (||A||*N*eps)
*
                        IF( N.GE.M ) THEN
                           CALL PCTZRZRV( M, N, MEM( IPA ), 1, 1, DESCA,
     $                                    MEM( IPTAU ), MEM( IPW ) )
                        END IF
                        CALL PCLAFCHK( 'No', 'No', M, N, MEM( IPA ), 1,
     $                                 1, DESCA, IASEED, ANORM, FRESID,
     $                                 MEM( IPW ) )
                     END IF
*
*                    Check for memory overwrite
*
                     CALL PCCHEKPAD( ICTXT, ROUTCHK, MP, NQ,
     $                               MEM( IPA-IPREPAD ), DESCA( LLD_ ),
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     CALL PCCHEKPAD( ICTXT, ROUTCHK, LTAU, 1,
     $                               MEM( IPTAU-IPREPAD ), LTAU,
     $                               IPREPAD, IPOSTPAD, PADVAL )
                     CALL PCCHEKPAD( ICTXT, ROUTCHK, WORKSIZ-IPOSTPAD,
     $                               1, MEM( IPW-IPREPAD ),
     $                               WORKSIZ-IPOSTPAD, IPREPAD,
     $                               IPOSTPAD, PADVAL )
*
                     IF( LSAMEN( 2, FACT, 'QP' ) ) THEN
*
                        CALL PCQPPIV( M, N, MEM( IPA ), 1, 1, DESCA,
     $                                MEM( IPPIV ) )
*
*                       Check for memory overwrite
*
                        CALL PCCHEKPAD( ICTXT, 'PCQPPIV', MP, NQ,
     $                                  MEM( IPA-IPREPAD ),
     $                                  DESCA( LLD_ ),
     $                                  IPREPAD, IPOSTPAD, PADVAL )
                        CALL PCCHEKPAD( ICTXT, 'PCQPPIV', LIPIV, 1,
     $                                  MEM( IPPIV-IPREPAD ), LIPIV,
     $                                  IPREPAD, IPOSTPAD, PADVAL )
*
                        CALL PCLAFCHK( 'No', 'No', M, N, MEM( IPA ), 1,
     $                                 1, DESCA, IASEED, ANORM, FRESID,
     $                                 MEM( IPW ) )
*
*                       Check for memory overwrite
*
                        CALL PCCHEKPAD( ICTXT, 'PCLAFCHK', MP, NQ,
     $                                  MEM( IPA-IPREPAD ),
     $                                  DESCA( LLD_ ),
     $                                  IPREPAD, IPOSTPAD, PADVAL )
                        CALL PCCHEKPAD( ICTXT, 'PCLAFCHK',
     $                                  WORKSIZ-IPOSTPAD, 1,
     $                                  MEM( IPW-IPREPAD ),
     $                                  WORKSIZ-IPOSTPAD, IPREPAD,
     $                                  IPOSTPAD, PADVAL )
                     END IF
*
*                    Test residual and detect NaN result
*
                     IF( LSAMEN( 2, FACT, 'TZ' ) .AND. N.LT.M ) THEN
                        KSKIP = KSKIP + 1
                        PASSED = 'BYPASS'
                     ELSE
                        IF( FRESID.LE.THRESH .AND.
     $                      (FRESID-FRESID).EQ.0.0E+0 ) THEN
                           KPASS = KPASS + 1
                           PASSED = 'PASSED'
                        ELSE
                           KFAIL = KFAIL + 1
                           PASSED = 'FAILED'
                        END IF
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
                  CALL SLCOMBINE( ICTXT, 'All', '>', 'W', 1, 1, WTIME )
                  CALL SLCOMBINE( ICTXT, 'All', '>', 'C', 1, 1, CTIME )
*
*                 Print results
*
                  IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
*
                     MINMN = MIN( M, N )
                     MAXMN = MAX( M, N )
*
                     IF( LSAMEN( 2, FACT, 'TZ' ) ) THEN
                        IF( M.GE.N ) THEN
                           NOPS = 0.0D+0
                        ELSE
*
*                          9 ( M^2 N - M^3 ) + 13 M N - M^2 for
*                          complete unitary factorization (M <= N).
*
                           NOPS = 9.0D+0 * (
     $                            DBLE( N )*( DBLE( M )**2 ) -
     $                            DBLE( M )**3 ) +
     $                            13.0D+0*DBLE( N )*DBLE( M ) -
     $                            DBLE( M )**2
                        END IF
*
                     ELSE
*
*                       8 M N^2 - 8/3 N^2 + 6 M N + 8 N^2 for QR type
*                       factorization when M >= N.
*
                        NOPS = 8.0D+0 * ( DBLE( MINMN )**2 ) *
     $                     ( DBLE( MAXMN )-DBLE( MINMN ) / 3.0D+0 ) +
     $                     ( 6.0D+0 * DBLE( MAXMN ) +
     $                       8.0D+0 * DBLE( MINMN ) ) *
     $                     DBLE( MINMN )
                     END IF
*
*                    Print WALL time
*
                     IF( WTIME( 1 ).GT.0.0D+0 ) THEN
                        TMFLOPS = NOPS / ( WTIME( 1 ) * 1.0D+6 )
                     ELSE
                        TMFLOPS = 0.0D+0
                     END IF
                     IF( WTIME( 1 ).GE.0.0D+0 )
     $                  WRITE( NOUT, FMT = 9993 ) 'WALL', M, N, MB, NB,
     $                         NPROW, NPCOL, WTIME( 1 ), TMFLOPS,
     $                         PASSED, FRESID
*
*                    Print CPU time
*
                     IF( CTIME( 1 ).GT.0.0D+0 ) THEN
                        TMFLOPS = NOPS / ( CTIME( 1 ) * 1.0D+6 )
                     ELSE
                        TMFLOPS = 0.0D+0
                     END IF
                     IF( CTIME( 1 ).GE.0.0D+0 )
     $                  WRITE( NOUT, FMT = 9993 ) 'CPU ', M, N, MB, NB,
     $                         NPROW, NPCOL, CTIME( 1 ), TMFLOPS,
     $                         PASSED, FRESID
*
                  END IF
*
   10          CONTINUE
*
   20       CONTINUE
*
            CALL BLACS_GRIDEXIT( ICTXT )
*
   30    CONTINUE
*
   40 CONTINUE
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
 9995 FORMAT( 'TIME      M      N  MB  NB     P     Q Fact Time ',
     $        '     MFLOPS  CHECK  Residual' )
 9994 FORMAT( '---- ------ ------ --- --- ----- ----- --------- ',
     $        '----------- ------  --------' )
 9993 FORMAT( A4, 1X, I6, 1X, I6, 1X, I3, 1X, I3, 1X, I5, 1X, I5, 1X,
     $        F9.2, 1X, F11.2, 1X, A6, 2X, G8.1 )
 9992 FORMAT( 'Finished ', I6, ' tests, with the following results:' )
 9991 FORMAT( I5, ' tests completed and passed residual checks.' )
 9990 FORMAT( I5, ' tests completed without checking.' )
 9989 FORMAT( I5, ' tests completed and failed residual checks.' )
 9988 FORMAT( I5, ' tests skipped because of illegal input values.' )
 9987 FORMAT( 'END OF TESTS.' )
 9986 FORMAT( A )
*
      STOP
*
*     End of PCQRDRIVER
*
      END
*
      SUBROUTINE PCQPPIV( M, N, A, IA, JA, DESCA, IPIV )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, JA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), IPIV( * )
      COMPLEX            A( * )
*     ..
*
*  Purpose
*  =======
*
*  PCQPPIV applies to sub( A ) = A(IA:IA+M-1,JA:JA+N-1) the pivots
*  returned by PCGEQPF in reverse order for checking purposes.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) COMPLEX pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, the local pieces of the M-by-N distributed matrix
*          sub( A ) which is to be permuted. On exit, the local pieces
*          of the distributed permuted submatrix sub( A ) * Inv( P ).
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  IPIV    (local input) INTEGER array, dimension LOCc(JA+N-1).
*          On exit, if IPIV(I) = K, the local i-th column of sub( A )*P
*          was the global K-th column of sub( A ). IPIV is tied to the
*          distributed matrix A.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            IACOL, ICOFFA, ICTXT, IITMP, IPVT, IPCOL,
     $                   IPROW, ITMP, J, JJ, JJA, KK, MYCOL, MYROW,
     $                   NPCOL, NPROW, NQ
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGEBR2D, IGEBS2D, IGERV2D,
     $                   IGESD2D, IGAMN2D, INFOG1L, PCSWAP
*     ..
*     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      EXTERNAL           INDXL2G, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      CALL INFOG1L( JA, DESCA( NB_ ), NPCOL, MYCOL, DESCA( CSRC_ ), JJA,
     $              IACOL )
      ICOFFA = MOD( JA-1, DESCA( NB_ ) )
      NQ = NUMROC( N+ICOFFA, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IF( MYCOL.EQ.IACOL )
     $   NQ = NQ - ICOFFA
*
      DO 20 J = JA, JA+N-2
*
         IPVT = JA+N-1
         ITMP = JA+N
*
*        Find first the local minimum candidate for pivoting
*
         CALL INFOG1L( J, DESCA( NB_ ), NPCOL, MYCOL, DESCA( CSRC_ ),
     $                 JJ, IACOL )
         DO 10 KK = JJ, JJA+NQ-1
            IF( IPIV( KK ).LT.IPVT )THEN
               IITMP = KK
               IPVT = IPIV( KK )
            END IF
   10    CONTINUE
*
*        Find the global minimum pivot
*
         CALL IGAMN2D( ICTXT, 'Rowwise', ' ', 1, 1, IPVT, 1, IPROW,
     $                 IPCOL, 1, -1, MYCOL )
*
*        Broadcast the corresponding index to the other process columns
*
         IF( MYCOL.EQ.IPCOL ) THEN
            ITMP = INDXL2G( IITMP, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                      NPCOL )
            CALL IGEBS2D( ICTXT, 'Rowwise', ' ', 1, 1, ITMP, 1 )
            IF( IPCOL.NE.IACOL ) THEN
               CALL IGERV2D( ICTXT, 1, 1, IPIV( IITMP ), 1, MYROW,
     $                       IACOL )
            ELSE
               IF( MYCOL.EQ.IACOL )
     $            IPIV( IITMP ) = IPIV( JJ )
            END IF
         ELSE
            CALL IGEBR2D( ICTXT, 'Rowwise', ' ', 1, 1, ITMP, 1, MYROW,
     $                    IPCOL )
            IF( MYCOL.EQ.IACOL .AND. IPCOL.NE.IACOL )
     $         CALL IGESD2D( ICTXT, 1, 1, IPIV( JJ ), 1, MYROW, IPCOL )
         END IF
*
*        Swap the columns of A
*
         CALL PCSWAP( M, A, IA, ITMP, DESCA, 1, A, IA, J, DESCA, 1 )
*
   20 CONTINUE
*
*     End of PCQPPIV
*
      END
