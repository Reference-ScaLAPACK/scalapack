      BLOCK DATA
      INTEGER NSUBS
      PARAMETER (NSUBS = 7)
      CHARACTER*7        SNAMES( NSUBS )
      COMMON             /SNAMEC/SNAMES
      DATA               SNAMES/'PDGEMV ', 'PDSYMV ', 'PDTRMV ',
     $                   'PDTRSV ', 'PDGER  ', 'PDSYR  ',
     $                   'PDSYR2 '/
      END BLOCK DATA
      
      PROGRAM PDBLA2TST
*
*  -- PBLAS testing driver (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*  Purpose
*  =======
*
*  PDBLA2TST is the main testing program for the PBLAS Level 2 routines.
*
*  The program must be driven by a short data file.  An  annotated exam-
*  ple of a data file can be obtained by deleting the first 3 characters
*  from the following 60 lines:
*  'Level 2 PBLAS, Testing input file'
*  'Intel iPSC/860 hypercube, gamma model.'
*  'PDBLAS2TST.SUMM'     output file name (if any)
*  6       device out
*  F       logical flag, T to stop on failures
*  F       logical flag, T to test error exits
*  0       verbosity, 0 for pass/fail, 1-3 for matrix dump on errors
*  10      the leading dimension gap
*  16.0    threshold value of test ratio
*  10              value of the logical computational blocksize NB
*  1               number of process grids (ordered pairs of P & Q)
*  2 2 1 4 2 3 8   values of P
*  2 2 4 1 3 2 1   values of Q
*  1.0D0           value of ALPHA
*  1.0D0           value of BETA
*  2               number of tests problems
*  'U' 'L'         values of UPLO
*  'N' 'T'         values of TRANS
*  'N' 'U'         values of DIAG
*  3  4            values of M
*  3  4            values of N
*  6 10            values of M_A
*  6 10            values of N_A
*  2  5            values of IMB_A
*  2  5            values of INB_A
*  2  5            values of MB_A
*  2  5            values of NB_A
*  0  1            values of RSRC_A
*  0  0            values of CSRC_A
*  1  1            values of IA
*  1  1            values of JA
*  6 10            values of M_X
*  6 10            values of N_X
*  2  5            values of IMB_X
*  2  5            values of INB_X
*  2  5            values of MB_X
*  2  5            values of NB_X
*  0  1            values of RSRC_X
*  0  0            values of CSRC_X
*  1  1            values of IX
*  1  1            values of JX
*  1  1            values of INCX
*  6 10            values of M_Y
*  6 10            values of N_Y
*  2  5            values of IMB_Y
*  2  5            values of INB_Y
*  2  5            values of MB_Y
*  2  5            values of NB_Y
*  0  1            values of RSRC_Y
*  0  0            values of CSRC_Y
*  1  1            values of IY
*  1  1            values of JY
*  6  1            values of INCY
*  PDGEMV  T  put F for no test in the same column
*  PDSYMV  T  put F for no test in the same column
*  PDTRMV  T  put F for no test in the same column
*  PDTRSV  T  put F for no test in the same column
*  PDGER   T  put F for no test in the same column
*  PDSYR   T  put F for no test in the same column
*  PDSYR2  T  put F for no test in the same column
*
*  Internal Parameters
*  ===================
*
*  TOTMEM  INTEGER
*          TOTMEM  is  a machine-specific parameter indicating the maxi-
*          mum  amount  of  available  memory per  process in bytes. The
*          user  should  customize TOTMEM to his  platform.  Remember to
*          leave  room  in  memory  for the  operating system, the BLACS
*          buffer, etc.  For  example,  on  a system with 8 MB of memory
*          per process (e.g., one processor  on an Intel iPSC/860),  the
*          parameters we use are TOTMEM=6200000  (leaving 1.8 MB for OS,
*          code, BLACS buffer, etc).  However,  for PVM,  we usually set
*          TOTMEM = 2000000.  Some experimenting  with the maximum value
*          of TOTMEM may be required. By default, TOTMEM is 2000000.
*
*  DBLESZ  INTEGER
*          DBLESZ  indicates  the  length in bytes on the given platform
*          for  a  double  precision  real. By default, DBLESZ is set to
*          eight.
*
*  MEM     DOUBLE PRECISION array
*          MEM is an array of dimension TOTMEM / DBLESZ.
*          All arrays used by SCALAPACK routines are allocated from this
*          array MEM and referenced by pointers. The  integer  IPA,  for
*          example, is a pointer to the starting element of MEM for  the
*          matrix A.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXTESTS, MAXGRIDS, GAPMUL, DBLESZ, TOTMEM,
     $                   MEMSIZ, NSUBS
      DOUBLE PRECISION   ONE, PADVAL, ZERO, ROGUE
      PARAMETER          ( MAXTESTS = 20, MAXGRIDS = 20, GAPMUL = 10,
     $                   DBLESZ = 8, TOTMEM = 2000000,
     $                   MEMSIZ = TOTMEM / DBLESZ, ZERO = 0.0D+0,
     $                   ONE = 1.0D+0, PADVAL = -9923.0D+0,
     $                   NSUBS = 7, ROGUE = -1.0D+10 )
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ERRFLG, SOF, TEE
      CHARACTER*1        AFORM, DIAG, DIAGDO, TRANS, UPLO
      INTEGER            CSRCA, CSRCX, CSRCY, I, IA, IAM, IASEED, ICTXT,
     $                   IGAP, IMBA, IMBX, IMBY, IMIDA, IMIDX, IMIDY,
     $                   INBA, INBX, INBY, INCX, INCY, IPA, IPG, IPMATA,
     $                   IPMATX, IPMATY, IPOSTA, IPOSTX, IPOSTY, IPREA,
     $                   IPREX, IPREY, IPX, IPY, IVERB, IX, IXSEED, IY,
     $                   IYSEED, J, JA, JX, JY, K, LDA, LDX, LDY, M, MA,
     $                   MBA, MBX, MBY, MEMREQD, MPA, MPX, MPY, MX, MY,
     $                   MYCOL, MYROW, N, NA, NBA, NBX, NBY, NCOLA,
     $                   NGRIDS, NLX, NLY, NOUT, NPCOL, NPROCS, NPROW,
     $                   NQA, NQX, NQY, NROWA, NTESTS, NX, NY, OFFD,
     $                   RSRCA, RSRCX, RSRCY, TSKIP, TSTCNT
      REAL               THRESH
      DOUBLE PRECISION   ALPHA, BETA, SCALE
*     ..
*     .. Local Arrays ..
      LOGICAL            LTEST( NSUBS ), YCHECK( NSUBS )
      CHARACTER*1        DIAGVAL( MAXTESTS ), TRANVAL( MAXTESTS ),
     $                   UPLOVAL( MAXTESTS )
      CHARACTER*80       OUTFILE
      INTEGER            CSCAVAL( MAXTESTS ), CSCXVAL( MAXTESTS ),
     $                   CSCYVAL( MAXTESTS ), DESCA( DLEN_ ),
     $                   DESCAR( DLEN_ ), DESCX( DLEN_ ),
     $                   DESCXR( DLEN_ ), DESCY( DLEN_ ),
     $                   DESCYR( DLEN_ ), IAVAL( MAXTESTS ), IERR( 6 ),
     $                   IMBAVAL( MAXTESTS ), IMBXVAL( MAXTESTS ),
     $                   IMBYVAL( MAXTESTS ), INBAVAL( MAXTESTS ),
     $                   INBXVAL( MAXTESTS ), INBYVAL( MAXTESTS ),
     $                   INCXVAL( MAXTESTS ), INCYVAL( MAXTESTS ),
     $                   IXVAL( MAXTESTS ), IYVAL( MAXTESTS ),
     $                   JAVAL( MAXTESTS ), JXVAL( MAXTESTS ),
     $                   JYVAL( MAXTESTS )
      INTEGER            KFAIL( NSUBS ), KPASS( NSUBS ), KSKIP( NSUBS ),
     $                   KTESTS( NSUBS ), MAVAL( MAXTESTS ),
     $                   MBAVAL( MAXTESTS ), MBXVAL( MAXTESTS ),
     $                   MBYVAL( MAXTESTS ), MVAL( MAXTESTS ),
     $                   MXVAL( MAXTESTS ), MYVAL( MAXTESTS ),
     $                   NAVAL( MAXTESTS ), NBAVAL( MAXTESTS ),
     $                   NBXVAL( MAXTESTS ), NBYVAL( MAXTESTS ),
     $                   NVAL( MAXTESTS ), NXVAL( MAXTESTS ),
     $                   NYVAL( MAXTESTS ), PVAL( MAXTESTS ),
     $                   QVAL( MAXTESTS ), RSCAVAL( MAXTESTS ),
     $                   RSCXVAL( MAXTESTS ), RSCYVAL( MAXTESTS )
      DOUBLE PRECISION   MEM( MEMSIZ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $                   IGSUM2D, PB_DCHEKPAD, PB_DESCSET2, PB_DFILLPAD,
     $                   PB_DLASCAL, PB_DLASET, PB_PDLAPRNT,
     $                   PDBLA2TSTINFO, PDBLAS2TSTCHK, PDBLAS2TSTCHKE,
     $                   PDCHKARG2, PDCHKVOUT, PDGEMV, PDGER, PDLAGEN,
     $                   PDLASCAL, PDLASET, PDMPRNT, PDSYMV, PDSYR,
     $                   PDSYR2, PDTRMV, PDTRSV, PDVPRNT, PMDESCCHK,
     $                   PMDIMCHK, PVDESCCHK, PVDIMCHK
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MOD
*     ..
*     .. Common Blocks ..
      CHARACTER*7        SNAMES( NSUBS )
      LOGICAL            ABRTFLG
      INTEGER            INFO, NBLOG
      COMMON             /SNAMEC/SNAMES
      COMMON             /INFOC/INFO, NBLOG
      COMMON             /PBERRORC/NOUT, ABRTFLG
*     ..
*     .. Data Statements ..
      DATA               YCHECK/.TRUE., .TRUE., .FALSE., .FALSE.,
     $                   .TRUE., .FALSE., .TRUE./
*     ..
*     .. Executable Statements ..
*
*     Initialization
*
*     Set flag so that the PBLAS error handler won't abort on errors, so
*     that the tester will detect unsupported operations.
*
      ABRTFLG = .FALSE.
*
*     So far no error, will become true as soon as one error is found.
*
      ERRFLG = .FALSE.
*
*     Test counters
*
      TSKIP  = 0
      TSTCNT = 0
*
*     Seeds for random matrix generations.
*
      IASEED = 100
      IXSEED = 200
      IYSEED = 300
*
*     So far no tests have been performed.
*
      DO 10 I = 1, NSUBS
         KPASS( I )  = 0
         KSKIP( I )  = 0
         KFAIL( I )  = 0
         KTESTS( I ) = 0
   10 CONTINUE
*
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
      CALL PDBLA2TSTINFO( OUTFILE, NOUT, NTESTS, DIAGVAL, TRANVAL,
     $                    UPLOVAL, MVAL, NVAL, MAVAL, NAVAL, IMBAVAL,
     $                    MBAVAL, INBAVAL, NBAVAL, RSCAVAL, CSCAVAL,
     $                    IAVAL, JAVAL, MXVAL, NXVAL, IMBXVAL, MBXVAL,
     $                    INBXVAL, NBXVAL, RSCXVAL, CSCXVAL, IXVAL,
     $                    JXVAL, INCXVAL, MYVAL, NYVAL, IMBYVAL,
     $                    MBYVAL, INBYVAL, NBYVAL, RSCYVAL, CSCYVAL,
     $                    IYVAL, JYVAL, INCYVAL, MAXTESTS, NGRIDS,
     $                    PVAL, MAXGRIDS, QVAL, MAXGRIDS, NBLOG, LTEST,
     $                    SOF, TEE, IAM, IGAP, IVERB, NPROCS, THRESH,
     $                    ALPHA, BETA, MEM )
*
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9975 )
         WRITE( NOUT, FMT = * )
      END IF
*
*     If TEE is set then Test Error Exits of routines.
*
      IF( TEE )
     $   CALL PDBLAS2TSTCHKE( LTEST, NOUT, NPROCS )
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
     $         WRITE( NOUT, FMT = 9999 ) 'GRID SIZE', 'NPROW', NPROW
            IERR( 1 ) = 1
         ELSE IF( NPCOL.LT.1 ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9999 ) 'GRID SIZE', 'NPCOL', NPCOL
            IERR( 1 ) = 1
         ELSE IF( NPROW*NPCOL.GT.NPROCS ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 ) NPROW*NPCOL, NPROCS
            IERR( 1 ) = 1
         END IF
*
         IF( IERR( 1 ).GT.0 ) THEN
            IF( IAM.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 ) 'GRID'
            TSKIP = TSKIP + 1
            GO TO 60
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
     $      GO TO 60
*
*        Loop over number of tests
*
         DO 50 J = 1, NTESTS
*
*           Get the test parameters
*
            DIAG  = DIAGVAL( J )
            TRANS = TRANVAL( J )
            UPLO  = UPLOVAL( J )
*
            M     = MVAL( J )
            N     = NVAL( J )
*
            MA    = MAVAL( J )
            NA    = NAVAL( J )
            IMBA  = IMBAVAL( J )
            INBA  = INBAVAL( J )
            MBA   = MBAVAL( J )
            NBA   = NBAVAL( J )
            RSRCA = RSCAVAL( J )
            CSRCA = CSCAVAL( J )
            IA    = IAVAL( J )
            JA    = JAVAL( J )
*
            MX    = MXVAL( J )
            NX    = NXVAL( J )
            IMBX  = IMBXVAL( J )
            INBX  = INBXVAL( J )
            MBX   = MBXVAL( J )
            NBX   = NBXVAL( J )
            RSRCX = RSCXVAL( J )
            CSRCX = CSCXVAL( J )
            IX    = IXVAL( J )
            JX    = JXVAL( J )
            INCX  = INCXVAL( J )
*
            MY    = MYVAL( J )
            NY    = NYVAL( J )
            IMBY  = IMBYVAL( J )
            INBY  = INBYVAL( J )
            MBY   = MBYVAL( J )
            NBY   = NBYVAL( J )
            RSRCY = RSCYVAL( J )
            CSRCY = CSCYVAL( J )
            IY    = IYVAL( J )
            JY    = JYVAL( J )
            INCY  = INCYVAL( J )
*
            IF( IAM.EQ.0 ) THEN
               TSTCNT = TSTCNT + 1
               WRITE( NOUT, FMT = * )
               WRITE( NOUT, FMT = 9996 ) TSTCNT, NPROW, NPCOL
               WRITE( NOUT, FMT = * )
*
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9994 )
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9993 ) M, N, UPLO, TRANS, DIAG
*
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9992 )
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9991 ) IA, JA, MA, NA, IMBA, INBA,
     $                                   MBA, NBA, RSRCA, CSRCA
*
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9990 )
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9989 ) IX, JX, MX, NX, IMBX, INBX,
     $                                   MBX, NBX, RSRCX, CSRCX, INCX
*
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9988 )
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9989 ) IY, JY, MY, NY, IMBY, INBY,
     $                                   MBY, NBY, RSRCY, CSRCY, INCY
*
               WRITE( NOUT, FMT = 9995 )
*
            END IF
*
*           Check the validity of the input test parameters
*
            IF( .NOT.LSAME( UPLO, 'U' ).AND.
     $          .NOT.LSAME( UPLO, 'L' ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'UPLO'
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
            IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $          .NOT.LSAME( TRANS, 'T' ).AND.
     $          .NOT.LSAME( TRANS, 'C' ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'TRANS'
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
            IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' ) )THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) TRANS
               WRITE( NOUT, FMT = 9997 ) 'DIAG'
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
*           Check and initialize the matrix descriptors
*
            CALL PMDESCCHK( ICTXT, NOUT, 'A', DESCA,
     $                      BLOCK_CYCLIC_2D_INB, MA, NA, IMBA, INBA,
     $                      MBA, NBA, RSRCA, CSRCA, MPA, NQA, IPREA,
     $                      IMIDA, IPOSTA, IGAP, GAPMUL, IERR( 1 ) )
            CALL PVDESCCHK( ICTXT, NOUT, 'X', DESCX,
     $                      BLOCK_CYCLIC_2D_INB, MX, NX, IMBX, INBX,
     $                      MBX, NBX, RSRCX, CSRCX, INCX, MPX, NQX,
     $                      IPREX, IMIDX, IPOSTX, IGAP, GAPMUL,
     $                      IERR( 2 ) )
            CALL PVDESCCHK( ICTXT, NOUT, 'Y', DESCY,
     $                      BLOCK_CYCLIC_2D_INB, MY, NY, IMBY, INBY,
     $                      MBY, NBY, RSRCY, CSRCY, INCY, MPY, NQY,
     $                      IPREY, IMIDY, IPOSTY, IGAP, GAPMUL,
     $                      IERR( 3 ) )
*
            IF( IERR( 1 ).GT.0 .OR. IERR( 2 ).GT.0 .OR.
     $          IERR( 3 ).GT.0 ) THEN
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
            LDA = MAX( 1, MA )
            LDX = MAX( 1, MX )
            LDY = MAX( 1, MY )
*
*           Assign pointers into MEM for matrices corresponding to
*           the distributed matrices A, X and Y.
*
            IPA = IPREA + 1
            IPX = IPA + DESCA( LLD_ )*NQA + IPOSTA + IPREX
            IPY = IPX + DESCX( LLD_ )*NQX + IPOSTX + IPREY
            IPMATA = IPY + DESCY( LLD_ )*NQY + IPOSTY
            IPMATX = IPMATA + MA*NA
            IPMATY = IPMATX + MX*NX
            IPG = IPMATY + MAX( MX*NX, MY*NY )
*
*           Check if sufficient memory.
*           Requirement = mem for local part of parallel matrices +
*                         mem for whole matrices for comp. check +
*                         mem for recving comp. check error vals.
*
            MEMREQD = IPG + MAX( M, N ) - 1 +
     $                MAX( MAX( IMBA, MBA ),
     $                     MAX( MAX( IMBX, MBX ),
     $                          MAX( IMBY, MBY ) ) )
            IERR( 1 ) = 0
            IF( MEMREQD.GT.MEMSIZ ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9986 ) MEMREQD*DBLESZ
               IERR( 1 ) = 1
            END IF
*
*           Check all processes for an error
*
            CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
            IF( IERR( 1 ).GT.0 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9987 )
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
*           Loop over all PBLAS 2 routines
*
            DO 30 K = 1, NSUBS
*
*              Continue only if this subroutine has to be tested.
*
               IF( .NOT.LTEST( K ) )
     $            GO TO 30
*
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = 9985 ) SNAMES( K )
               END IF
*
*              Define the size of the operands
*
               IF( K.EQ.1 ) THEN
                  NROWA = M
                  NCOLA = N
                  IF( LSAME( TRANS, 'N' ) ) THEN
                     NLX = N
                     NLY = M
                  ELSE
                     NLX = M
                     NLY = N
                  END IF
               ELSE IF( K.EQ.5 ) THEN
                  NROWA = M
                  NCOLA = N
                  NLX = M
                  NLY = N
               ELSE
                  NROWA = N
                  NCOLA = N
                  NLX = N
                  NLY = N
               END IF
*
*              Check the validity of the operand sizes
*
               CALL PMDIMCHK( ICTXT, NOUT, NROWA, NCOLA, 'A', IA, JA,
     $                        DESCA, IERR( 1 ) )
               CALL PVDIMCHK( ICTXT, NOUT, NLX, 'X', IX, JX, DESCX,
     $                        INCX, IERR( 2 ) )
               CALL PVDIMCHK( ICTXT, NOUT, NLY, 'Y', IY, JY, DESCY,
     $                        INCY, IERR( 3 ) )
*
               IF( IERR( 1 ).NE.0 .OR. IERR( 2 ).NE.0 .OR.
     $             IERR( 3 ).NE.0 ) THEN
                  KSKIP( K ) = KSKIP( K ) + 1
                  GO TO 30
               END IF
*
*              Generate distributed matrices A, X and Y
*
               IF( K.EQ.2 .OR. K.EQ.6 .OR. K.EQ.7 ) THEN
                  AFORM  = 'S'
                  DIAGDO = 'N'
                  OFFD   = IA - JA
               ELSE IF( ( K.EQ.4 ).AND.( LSAME( DIAG, 'N' ) ) ) THEN
                  AFORM  = 'N'
                  DIAGDO = 'D'
                  OFFD   = IA - JA
               ELSE
                  AFORM  = 'N'
                  DIAGDO = 'N'
                  OFFD   = 0
               END IF
*
               CALL PDLAGEN( .FALSE., AFORM, DIAGDO, OFFD, MA, NA,
     $                       1, 1, DESCA, IASEED, MEM( IPA ),
     $                       DESCA( LLD_ ) )
               CALL PDLAGEN( .FALSE., 'None', 'No diag', 0, MX, NX, 1,
     $                       1, DESCX, IXSEED, MEM( IPX ),
     $                       DESCX( LLD_ ) )
               IF( YCHECK( K ) )
     $            CALL PDLAGEN( .FALSE., 'None', 'No diag', 0, MY, NY,
     $                          1, 1, DESCY, IYSEED, MEM( IPY ),
     $                          DESCY( LLD_ ) )
*
*              Generate entire matrices on each process.
*
               CALL PB_DESCSET2( DESCAR, MA, NA, IMBA, INBA, MBA, NBA,
     $                           -1, -1, ICTXT, MAX( 1, MA ) )
               CALL PDLAGEN( .FALSE., AFORM, DIAGDO, OFFD, MA, NA,
     $                       1, 1, DESCAR, IASEED, MEM( IPMATA ),
     $                       DESCAR( LLD_ ) )
               CALL PB_DESCSET2( DESCXR, MX, NX, IMBX, INBX, MBX, NBX,
     $                           -1, -1, ICTXT, MAX( 1, MX ) )
               CALL PDLAGEN( .FALSE., 'None', 'No diag', 0, MX, NX, 1,
     $                       1, DESCXR, IXSEED, MEM( IPMATX ),
     $                       DESCXR( LLD_ ) )
               IF( YCHECK( K ) ) THEN
*
                  CALL PB_DESCSET2( DESCYR, MY, NY, IMBY, INBY, MBY,
     $                              NBY, -1, -1, ICTXT, MAX( 1, MY ) )
                  CALL PDLAGEN( .FALSE., 'None', 'No diag', 0, MY, NY,
     $                          1, 1, DESCYR, IYSEED, MEM( IPMATY ),
     $                          DESCYR( LLD_ ) )
*
               ELSE
*
*                 If Y is not needed, generate a copy of X instead
*
                  CALL PB_DESCSET2( DESCYR, MX, NX, IMBX, INBX, MBX,
     $                              NBX, -1, -1, ICTXT, MAX( 1, MX ) )
                  CALL PDLAGEN( .FALSE., 'None', 'No diag', 0, MX, NX,
     $                          1, 1, DESCYR, IXSEED, MEM( IPMATY ),
     $                          DESCYR( LLD_ ) )
*
               END IF
*
*              Zero non referenced part of the matrices A
*
               IF( ( K.EQ.2 .OR. K.EQ.6 .OR. K.EQ.7 ).AND.
     $             ( MAX( NROWA, NCOLA ).GT.1 ) ) THEN
*
*                 The distributed matrix A is symmetric
*
                  IF( LSAME( UPLO, 'L' ) ) THEN
*
*                    Zeros the strict upper triangular part of A.
*
                     CALL PDLASET( 'Upper', NROWA-1, NCOLA-1, ROGUE,
     $                             ROGUE, MEM( IPA ), IA, JA+1, DESCA )
                     IF( K.NE.2 ) THEN
                        CALL PB_DLASET( 'Upper', NROWA-1, NCOLA-1, 0,
     $                                  ROGUE, ROGUE,
     $                                  MEM( IPMATA+IA-1+JA*LDA ), LDA )
                     END IF
*
                  ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
*                    Zeros the strict lower triangular part of A.
*
                     CALL PDLASET( 'Lower', NROWA-1, NCOLA-1, ROGUE,
     $                             ROGUE, MEM( IPA ), IA+1, JA, DESCA )
                     IF( K.NE.2 ) THEN
                        CALL PB_DLASET( 'Lower', NROWA-1, NCOLA-1, 0,
     $                                  ROGUE, ROGUE,
     $                                  MEM( IPMATA+IA+(JA-1)*LDA ),
     $                                  LDA )
                     END IF
*
                  END IF
*
               ELSE IF( K.EQ.3 .OR. K.EQ.4 ) THEN
*
                  IF( LSAME( UPLO, 'L' ) ) THEN
*
*                    The distributed matrix A is lower triangular
*
                     IF( LSAME( DIAG, 'N' ) ) THEN
*
                        IF( MAX( NROWA, NCOLA ).GT.1 ) THEN
                           CALL PDLASET( 'Upper', NROWA-1, NCOLA-1,
     $                                   ROGUE, ROGUE, MEM( IPA ), IA,
     $                                   JA+1, DESCA )
                           CALL PB_DLASET( 'Upper', NROWA-1, NCOLA-1, 0,
     $                                     ZERO, ZERO,
     $                                     MEM( IPMATA+IA-1+JA*LDA ),
     $                                     LDA )
                        END IF
*
                     ELSE
*
                        CALL PDLASET( 'Upper', NROWA, NCOLA, ROGUE, ONE,
     $                                MEM( IPA ), IA, JA, DESCA )
                        CALL PB_DLASET( 'Upper', NROWA, NCOLA, 0, ZERO,
     $                                  ONE,
     $                                  MEM( IPMATA+IA-1+(JA-1)*LDA ),
     $                                  LDA )
                        IF( ( K.EQ.4 ).AND.
     $                      ( MAX( NROWA, NCOLA ).GT.1 ) ) THEN
                           SCALE = ONE / DBLE( MAX( NROWA, NCOLA ) )
                           CALL PDLASCAL( 'Lower', NROWA-1, NCOLA-1,
     $                                    SCALE, MEM( IPA ), IA+1, JA,
     $                                    DESCA )
                           CALL PB_DLASCAL( 'Lower', NROWA-1, NCOLA-1,
     $                                  0, SCALE,
     $                                  MEM( IPMATA+IA+(JA-1)*LDA ),
     $                                  LDA )
                        END IF
*
                     END IF
*
                  ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
*                    The distributed matrix A is upper triangular
*
                     IF( LSAME( DIAG, 'N' ) ) THEN
*
                        IF( MAX( NROWA, NCOLA ).GT.1 ) THEN
                           CALL PDLASET( 'Lower', NROWA-1, NCOLA-1,
     $                                   ROGUE, ROGUE, MEM( IPA ), IA+1,
     $                                   JA, DESCA )
                           CALL PB_DLASET( 'Lower', NROWA-1, NCOLA-1, 0,
     $                                     ZERO, ZERO,
     $                                     MEM( IPMATA+IA+(JA-1)*LDA ),
     $                                     LDA )
                        END IF
*
                     ELSE
*
                        CALL PDLASET( 'Lower', NROWA, NCOLA, ROGUE, ONE,
     $                                MEM( IPA ), IA, JA, DESCA )
                        CALL PB_DLASET( 'Lower', NROWA, NCOLA, 0, ZERO,
     $                                  ONE,
     $                                  MEM( IPMATA+IA-1+(JA-1)*LDA ),
     $                                  LDA )
                        IF( ( K.EQ.4 ).AND.
     $                      ( MAX( NROWA, NCOLA ).GT.1 ) ) THEN
                           SCALE = ONE / DBLE( MAX( NROWA, NCOLA ) )
                           CALL PDLASCAL( 'Upper', NROWA-1, NCOLA-1,
     $                                    SCALE, MEM( IPA ), IA, JA+1,
     $                                    DESCA )
                           CALL PB_DLASCAL( 'Upper', NROWA-1, NCOLA-1,
     $                                  0, SCALE,
     $                                  MEM( IPMATA+IA-1+JA*LDA ), LDA )
                        END IF
*
                     END IF
*
                  END IF
*
               END IF
*
*              Pad the guard zones of A, X and Y
*
               CALL PB_DFILLPAD( ICTXT, MPA, NQA, MEM( IPA-IPREA ),
     $                           DESCA( LLD_ ), IPREA, IPOSTA, PADVAL )
*
               CALL PB_DFILLPAD( ICTXT, MPX, NQX, MEM( IPX-IPREX ),
     $                           DESCX( LLD_ ), IPREX, IPOSTX, PADVAL )
*
               IF( YCHECK( K ) ) THEN
                  CALL PB_DFILLPAD( ICTXT, MPY, NQY, MEM( IPY-IPREY ),
     $                              DESCY( LLD_ ), IPREY, IPOSTY,
     $                              PADVAL )
               END IF
*
*              Initialize the check for INPUT-only arguments.
*
               INFO = 0
               CALL PDCHKARG2( ICTXT, NOUT, SNAMES( K ), UPLO, TRANS,
     $                         DIAG, M, N, ALPHA, IA, JA, DESCA, IX,
     $                         JX, DESCX, INCX, BETA, IY, JY, DESCY,
     $                         INCY, INFO )
*
*              Print initial parallel data if IVERB >= 2.
*
               IF( IVERB.EQ.2 ) THEN
                  CALL PB_PDLAPRNT( NROWA, NCOLA, MEM( IPA ), IA, JA,
     $                              DESCA, 0, 0, 'PARALLEL_INITIAL_A',
     $                              NOUT, MEM( IPG ) )
               ELSE IF( IVERB.GE.3 ) THEN
                  CALL PB_PDLAPRNT( MA, NA, MEM( IPA ), 1, 1, DESCA, 0,
     $                              0, 'PARALLEL_INITIAL_A', NOUT,
     $                              MEM( IPG ) )
               END IF
*
               IF( IVERB.EQ.2 ) THEN
                  IF( INCX.EQ.DESCX( M_ ) ) THEN
                     CALL PB_PDLAPRNT( 1, NLX, MEM( IPX ), IX, JX,
     $                                 DESCX, 0, 0,
     $                                 'PARALLEL_INITIAL_X', NOUT,
     $                                 MEM( IPG ) )
                  ELSE
                     CALL PB_PDLAPRNT( NLX, 1, MEM( IPX ), IX, JX,
     $                                 DESCX, 0, 0,
     $                                 'PARALLEL_INITIAL_X', NOUT,
     $                                 MEM( IPG ) )
                  END IF
               ELSE IF( IVERB.GE.3 ) THEN
                  CALL PB_PDLAPRNT( MX, NX, MEM( IPX ), 1, 1, DESCX, 0,
     $                              0, 'PARALLEL_INITIAL_X', NOUT,
     $                              MEM( IPG ) )
               END IF
*
               IF( YCHECK( K ) ) THEN
                  IF( IVERB.EQ.2 ) THEN
                     IF( INCY.EQ.DESCY( M_ ) ) THEN
                        CALL PB_PDLAPRNT( 1, NLY, MEM( IPY ), IY, JY,
     $                                    DESCY, 0, 0,
     $                                    'PARALLEL_INITIAL_Y', NOUT,
     $                                    MEM( IPG ) )
                     ELSE
                        CALL PB_PDLAPRNT( NLY, 1, MEM( IPY ), IY, JY,
     $                                    DESCY, 0, 0,
     $                                    'PARALLEL_INITIAL_Y', NOUT,
     $                                    MEM( IPG ) )
                     END IF
                  ELSE IF( IVERB.GE.3 ) THEN
                     CALL PB_PDLAPRNT( MY, NY, MEM( IPY ), 1, 1, DESCY,
     $                                 0, 0, 'PARALLEL_INITIAL_Y', NOUT,
     $                                 MEM( IPG ) )
                  END IF
               END IF
*
*              Call the Level 2 PBLAS routine
*
               INFO = 0
               IF( K.EQ.1 ) THEN
*
*                 Test PDGEMV
*
                  CALL PDGEMV( TRANS, M, N, ALPHA, MEM( IPA ), IA, JA,
     $                         DESCA, MEM( IPX ), IX, JX, DESCX, INCX,
     $                         BETA, MEM( IPY ), IY, JY, DESCY, INCY )
*
               ELSE IF( K.EQ.2 ) THEN
*
*                 Test PDSYMV
*
                  CALL PDSYMV( UPLO, N, ALPHA, MEM( IPA ), IA, JA,
     $                         DESCA, MEM( IPX ), IX, JX, DESCX, INCX,
     $                         BETA, MEM( IPY ), IY, JY, DESCY, INCY )
*
               ELSE IF( K.EQ.3 ) THEN
*
*                 Test PDTRMV
*
                  CALL PDTRMV( UPLO, TRANS, DIAG, N, MEM( IPA ), IA, JA,
     $                         DESCA, MEM( IPX ), IX, JX, DESCX, INCX )
*
               ELSE IF( K.EQ.4 ) THEN
*
*                 Test PDTRSV
*
                  CALL PDTRSV( UPLO, TRANS, DIAG, N, MEM( IPA ), IA, JA,
     $                         DESCA, MEM( IPX ), IX, JX, DESCX, INCX )
*
               ELSE IF( K.EQ.5 ) THEN
*
*                 Test PDGER
*
                  CALL PDGER( M, N, ALPHA, MEM( IPX ), IX, JX, DESCX,
     $                        INCX, MEM( IPY ), IY, JY, DESCY, INCY,
     $                        MEM( IPA ), IA, JA, DESCA )
*
               ELSE IF( K.EQ.6 ) THEN
*
*                 Test PDSYR
*
                  CALL PDSYR( UPLO, N, ALPHA, MEM( IPX ), IX, JX, DESCX,
     $                         INCX, MEM( IPA ), IA, JA, DESCA )
*
               ELSE IF( K.EQ.7 ) THEN
*
*                 Test PDSYR2
*
                  CALL PDSYR2( UPLO, N, ALPHA, MEM( IPX ), IX, JX,
     $                         DESCX, INCX, MEM( IPY ), IY, JY, DESCY,
     $                         INCY, MEM( IPA ), IA, JA, DESCA )
*
               END IF
*
*              Check if the operation has been performed.
*
               IF( INFO.NE.0 ) THEN
                  KSKIP( K ) = KSKIP( K ) + 1
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9974 ) INFO
                  GO TO 30
               END IF
*
*              Check padding
*
               CALL PB_DCHEKPAD( ICTXT, SNAMES( K ), MPA, NQA,
     $                           MEM( IPA-IPREA ), DESCA( LLD_ ), IPREA,
     $                           IPOSTA, PADVAL )
*
               CALL PB_DCHEKPAD( ICTXT, SNAMES( K ), MPX, NQX,
     $                           MEM( IPX-IPREX ), DESCX( LLD_ ), IPREX,
     $                           IPOSTX, PADVAL )
*
               IF( YCHECK( K ) ) THEN
                  CALL PB_DCHEKPAD( ICTXT, SNAMES( K ), MPY, NQY,
     $                              MEM( IPY-IPREY ), DESCY( LLD_ ),
     $                              IPREY, IPOSTY, PADVAL )
               END IF
*
*              Check the computations
*
               CALL PDBLAS2TSTCHK( ICTXT, NOUT, K, UPLO, TRANS, DIAG, M,
     $                             N, ALPHA, MEM( IPMATA ), MEM( IPA ),
     $                             IA, JA, DESCA, MEM( IPMATX ),
     $                             MEM( IPX ), IX, JX, DESCX, INCX,
     $                             BETA, MEM( IPMATY ), MEM( IPY ), IY,
     $                             JY, DESCY, INCY, THRESH, ROGUE,
     $                             MEM( IPG ), INFO )
               IF( MOD( INFO, 2 ).EQ.1 ) THEN
                  IERR( 1 ) = 1
               ELSE IF( MOD( INFO / 2, 2 ).EQ.1 ) THEN
                  IERR( 2 ) = 1
               ELSE IF( MOD( INFO / 4, 2 ).EQ.1 ) THEN
                  IERR( 3 ) = 1
               ELSE IF( INFO.NE.0 ) THEN
                  IERR( 1 ) = 1
                  IERR( 2 ) = 1
                  IERR( 3 ) = 1
               END IF
*
*              Check input-only scalar arguments
*
               INFO = 1
               CALL PDCHKARG2( ICTXT, NOUT, SNAMES( K ), UPLO, TRANS,
     $                         DIAG, M, N, ALPHA, IA, JA, DESCA, IX,
     $                         JX, DESCX, INCX, BETA, IY, JY, DESCY,
     $                         INCY, INFO )
*
*              Check input-only array arguments
*
               CALL PDCHKMOUT( NROWA, NCOLA, MEM( IPMATA ), MEM( IPA ),
     $                         IA, JA, DESCA, IERR( 4 ) )
               CALL PDCHKVOUT( NLX, MEM( IPMATX ), MEM( IPX ), IX, JX,
     $                         DESCX, INCX, IERR( 5 ) )
*
               IF( IERR( 4 ).NE.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9982 ) 'PARALLEL_A',
     $                                         SNAMES( K )
               END IF
*
               IF( IERR( 5 ).NE.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9982 ) 'PARALLEL_X',
     $                                         SNAMES( K )
               END IF
*
               IF( YCHECK( K ) ) THEN
                  CALL PDCHKVOUT( NLY, MEM( IPMATY ), MEM( IPY ), IY,
     $                            JY, DESCY, INCY, IERR( 6 ) )
                  IF( IERR( 6 ).NE.0 ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9982 ) 'PARALLEL_Y',
     $                                            SNAMES( K )
                  END IF
               END IF
*
*              Only node 0 prints computational test result
*
               IF( INFO.NE.0 .OR. IERR( 1 ).NE.0 .OR.
     $             IERR( 2 ).NE.0 .OR. IERR( 3 ).NE.0 .OR.
     $             IERR( 4 ).NE.0 .OR. IERR( 5 ).NE.0 .OR.
     $             IERR( 6 ).NE.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9984 ) SNAMES( K )
                  KFAIL( K ) = KFAIL( K ) + 1
                  ERRFLG = .TRUE.
               ELSE
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9983 ) SNAMES( K )
                  KPASS( K ) = KPASS( K ) + 1
               END IF
*
*              Dump matrix if IVERB >= 1 and error.
*
               IF( IVERB.GE.1 .AND. ERRFLG ) THEN
                  IF( IERR( 4 ).NE.0 .OR. IVERB.GE.3 ) THEN
                     CALL PDMPRNT( ICTXT, NOUT, MA, NA, MEM( IPMATA ),
     $                             LDA, 0, 0, 'SERIAL_A' )
                     CALL PB_PDLAPRNT( MA, NA, MEM( IPA ), 1, 1, DESCA,
     $                                 0, 0, 'PARALLEL_A', NOUT,
     $                                MEM( IPMATA ) )
                  ELSE IF( IERR( 1 ).NE.0 ) THEN
                     IF( ( NROWA.GT.0 ).AND.( NCOLA.GT.0 ) )
     $                  CALL PDMPRNT( ICTXT, NOUT, NROWA, NCOLA,
     $                                MEM( IPMATA+IA-1+(JA-1)*LDA ),
     $                                LDA, 0, 0, 'SERIAL_A' )
                     CALL PB_PDLAPRNT( NROWA, NCOLA, MEM( IPA ), IA, JA,
     $                                 DESCA, 0, 0, 'PARALLEL_A',
     $                                 NOUT, MEM( IPMATA ) )
                  END IF
                  IF( IERR( 5 ).NE.0 .OR. IVERB.GE.3 ) THEN
                     CALL PDMPRNT( ICTXT, NOUT, MX, NX, MEM( IPMATX ),
     $                             LDX, 0, 0, 'SERIAL_X' )
                     CALL PB_PDLAPRNT( MX, NX, MEM( IPX ), 1, 1, DESCX,
     $                                 0, 0, 'PARALLEL_X', NOUT,
     $                                 MEM( IPMATX ) )
                  ELSE IF( IERR( 2 ).NE.0 ) THEN
                     IF( NLX.GT.0 )
     $                  CALL PDVPRNT( ICTXT, NOUT, NLX,
     $                                MEM( IPMATX+IX-1+(JX-1)*LDX ),
     $                                INCX, 0, 0, 'SERIAL_X' )
                     IF( INCX.EQ.DESCX( M_ ) ) THEN
                        CALL PB_PDLAPRNT( 1, NLX, MEM( IPX ), IX, JX,
     $                                    DESCX, 0, 0, 'PARALLEL_X',
     $                                    NOUT, MEM( IPMATX ) )
                     ELSE
                        CALL PB_PDLAPRNT( NLX, 1, MEM( IPX ), IX, JX,
     $                                    DESCX, 0, 0, 'PARALLEL_X',
     $                                    NOUT, MEM( IPMATX ) )
                     END IF
                  END IF
                  IF( YCHECK( K ) ) THEN
                     IF( IERR( 6 ).NE.0 .OR. IVERB.GE.3 ) THEN
                        CALL PDMPRNT( ICTXT, NOUT, MY, NY,
     $                                MEM( IPMATY ), LDY, 0, 0,
     $                                'SERIAL_Y' )
                        CALL PB_PDLAPRNT( MY, NY, MEM( IPY ), 1, 1,
     $                                    DESCY, 0, 0, 'PARALLEL_Y',
     $                                    NOUT, MEM( IPMATX ) )
                     ELSE IF( IERR( 3 ).NE.0 ) THEN
                        IF( NLY.GT.0 )
     $                     CALL PDVPRNT( ICTXT, NOUT, NLY,
     $                                   MEM( IPMATY+IY-1+(JY-1)*LDY ),
     $                                   INCY, 0, 0, 'SERIAL_Y' )
                        IF( INCY.EQ.DESCY( M_ ) ) THEN
                           CALL PB_PDLAPRNT( 1, NLY, MEM( IPY ), IY, JY,
     $                                       DESCY, 0, 0, 'PARALLEL_Y',
     $                                       NOUT, MEM( IPMATX ) )
                        ELSE
                           CALL PB_PDLAPRNT( NLY, 1, MEM( IPY ), IY, JY,
     $                                       DESCY, 0, 0, 'PARALLEL_Y',
     $                                       NOUT, MEM( IPMATX ) )
                        END IF
                     END IF
                  END IF
               END IF
*
*              Leave if error and "Stop On Failure"
*
               IF( SOF.AND.ERRFLG )
     $            GO TO 70
*
   30       CONTINUE
*
   40       IF( IAM.EQ.0 ) THEN
               WRITE( NOUT, FMT = * )
               WRITE( NOUT, FMT = 9981 ) J
            END IF
*
   50   CONTINUE
*
        CALL BLACS_GRIDEXIT( ICTXT )
*
   60 CONTINUE
*
*     Come here, if error and "Stop On Failure"
*
   70 CONTINUE
*
*     Before printing out final stats, add TSKIP to all skips
*
      DO 80 I = 1, NSUBS
         IF( LTEST( I ) ) THEN
            KSKIP( I ) = KSKIP( I ) + TSKIP
            KTESTS( I ) = KSKIP( I ) + KFAIL( I ) + KPASS( I )
         END IF
   80 CONTINUE
*
*     Print results
*
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9977 )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9979 )
         WRITE( NOUT, FMT = 9978 )
*
         DO 90 I = 1, NSUBS
            WRITE( NOUT, FMT = 9980 ) '|', SNAMES( I ), KTESTS( I ),
     $                                KPASS( I ), KFAIL( I ), KSKIP( I )
   90    CONTINUE
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9976 )
         WRITE( NOUT, FMT = * )
*
      END IF
*
      CALL BLACS_EXIT( 0 )
*
 9999 FORMAT( 'ILLEGAL ', A, ': ', A, ' = ', I10,
     $        ' should be at least 1' )
 9998 FORMAT( 'ILLEGAL GRID: NPROW*NPCOL = ', I4,
     $        '. It can be at most', I4 )
 9997 FORMAT( 'Bad ', A, ' parameters: going on to next test case.' )
 9996 FORMAT( 2X, 'Test number ', I4 , ' started on a ', I6, ' x ',
     $        I6, ' process grid.' )
 9995 FORMAT( 2X, '   ------------------------------------------------',
     $        '--------------------------' )
 9994 FORMAT( 2X, '        M      N       UPLO       TRANS       DIAG' )
 9993 FORMAT( 5X,I6,1X,I6,9X,A1,11X,A1,10X,A1 )
 9992 FORMAT( 2X, '       IA     JA     MA     NA   IMBA   INBA',
     $        '    MBA    NBA RSRCA CSRCA' )
 9991 FORMAT( 5X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,I6,
     $        1X,I5,1X,I5 )
 9990 FORMAT( 2X, '       IX     JX     MX     NX   IMBX   INBX',
     $        '    MBX    NBX RSRCX CSRCX   INCX' )
 9989 FORMAT( 5X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,I6,
     $        1X,I5,1X,I5,1X,I6 )
 9988 FORMAT( 2X, '       IY     JY     MY     NY   IMBY   INBY',
     $        '    MBY    NBY RSRCY CSRCY   INCY' )
 9987 FORMAT( 'Not enough memory for this test: going on to',
     $        ' next test case.' )
 9986 FORMAT( 'Not enough memory. Need: ', I12 )
 9985 FORMAT( 2X, '   Tested Subroutine: ', A )
 9984 FORMAT( 2X, '   ***** Computational check: ', A, '       ',
     $        ' FAILED ',' *****' )
 9983 FORMAT( 2X, '   ***** Computational check: ', A, '       ',
     $        ' PASSED ',' *****' )
 9982 FORMAT( 2X, '   ***** ERROR ***** Matrix operand ', A,
     $        ' modified by ', A, ' *****' )
 9981 FORMAT( 2X, 'Test number ', I4, ' completed.' )
 9980 FORMAT( 2X,A1,2X,A7,8X,I4,6X,I4,5X,I4,4X,I4 )
 9979 FORMAT( 2X, '   SUBROUTINE  TOTAL TESTS  PASSED   FAILED  ',
     $        'SKIPPED' )
 9978 FORMAT( 2X, '   ----------  -----------  ------   ------  ',
     $        '-------' )
 9977 FORMAT( 2X, 'Testing Summary')
 9976 FORMAT( 2X, 'End of Tests.' )
 9975 FORMAT( 2X, 'Tests started.' )
 9974 FORMAT( 2X, '   ***** Operation not supported, error code: ',
     $        I5, ' *****' )
*
      STOP
*
*     End of PDBLA2TST
*
      END
      SUBROUTINE PDBLA2TSTINFO( SUMMRY, NOUT, NMAT, DIAGVAL, TRANVAL,
     $                          UPLOVAL, MVAL, NVAL, MAVAL, NAVAL,
     $                          IMBAVAL, MBAVAL, INBAVAL, NBAVAL,
     $                          RSCAVAL, CSCAVAL, IAVAL, JAVAL,
     $                          MXVAL, NXVAL, IMBXVAL, MBXVAL,
     $                          INBXVAL, NBXVAL, RSCXVAL, CSCXVAL,
     $                          IXVAL, JXVAL, INCXVAL, MYVAL, NYVAL,
     $                          IMBYVAL, MBYVAL, INBYVAL, NBYVAL,
     $                          RSCYVAL, CSCYVAL, IYVAL, JYVAL,
     $                          INCYVAL, LDVAL, NGRIDS, PVAL, LDPVAL,
     $                          QVAL, LDQVAL, NBLOG, LTEST, SOF, TEE,
     $                          IAM, IGAP, IVERB, NPROCS, THRESH, ALPHA,
     $                          BETA, WORK )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      LOGICAL            SOF, TEE
      INTEGER            IAM, IGAP, IVERB, LDPVAL, LDQVAL, LDVAL, NBLOG,
     $                   NGRIDS, NMAT, NOUT, NPROCS
      REAL               THRESH
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      CHARACTER*( * )    SUMMRY
      CHARACTER*1        DIAGVAL( LDVAL ), TRANVAL( LDVAL ),
     $                   UPLOVAL( LDVAL )
      LOGICAL            LTEST( * )
      INTEGER            CSCAVAL( LDVAL ), CSCXVAL( LDVAL ),
     $                   CSCYVAL( LDVAL ), IAVAL( LDVAL ),
     $                   IMBAVAL( LDVAL ), IMBXVAL( LDVAL ),
     $                   IMBYVAL( LDVAL ), INBAVAL( LDVAL ),
     $                   INBXVAL( LDVAL ), INBYVAL( LDVAL ),
     $                   INCXVAL( LDVAL ), INCYVAL( LDVAL ),
     $                   IXVAL( LDVAL ), IYVAL( LDVAL ), JAVAL( LDVAL ),
     $                   JXVAL( LDVAL ), JYVAL( LDVAL ), MAVAL( LDVAL ),
     $                   MBAVAL( LDVAL ), MBXVAL( LDVAL ),
     $                   MBYVAL( LDVAL ), MVAL( LDVAL ), MXVAL( LDVAL ),
     $                   MYVAL( LDVAL ), NAVAL( LDVAL ),
     $                   NBAVAL( LDVAL ), NBXVAL( LDVAL ),
     $                   NBYVAL( LDVAL ), NVAL( LDVAL ), NXVAL( LDVAL ),
     $                   NYVAL( LDVAL ), PVAL( LDPVAL ), QVAL( LDQVAL ),
     $                   RSCAVAL( LDVAL ), RSCXVAL( LDVAL ),
     $                   RSCYVAL( LDVAL ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDBLA2TSTINFO  get the needed startup information for testing various
*  Level 2 PBLAS routines, and transmits it to all processes.
*
*  Notes
*  =====
*
*  For packing the information we assumed that the length in bytes of an
*  integer is equal to the length in bytes of a real single precision.
*
*  Arguments
*  =========
*
*  SUMMRY  (global output) CHARACTER*(*)
*          On  exit,  SUMMRY  is  the  name of output (summary) file (if
*          any). SUMMRY is only defined for process 0.
*
*  NOUT    (global output) INTEGER
*          On exit, NOUT  specifies the unit number for the output file.
*          When NOUT is 6, output to screen,  when  NOUT is 0, output to
*          stderr. NOUT is only defined for process 0.
*
*  NMAT    (global output) INTEGER
*          On exit,  NMAT  specifies the number of different test cases.
*
*  DIAGVAL (global output) CHARACTER array
*          On entry,  DIAGVAL  is  an array of dimension LDVAL. On exit,
*          this array contains the values of DIAG to run the code with.
*
*  TRANVAL (global output) CHARACTER array
*          On entry, TRANVAL  is an array of dimension LDVAL.  On  exit,
*          this array contains  the  values  of  TRANS  to  run the code
*          with.
*
*  UPLOVAL (global output) CHARACTER array
*          On entry, UPLOVAL  is an array of dimension LDVAL.  On  exit,
*          this array contains the values of UPLO to run the code with.
*
*  MVAL    (global output) INTEGER array
*          On entry, MVAL is an array of dimension LDVAL.  On exit, this
*          array contains the values of M to run the code with.
*
*  NVAL    (global output) INTEGER array
*          On entry, NVAL is an array of dimension LDVAL.  On exit, this
*          array contains the values of N to run the code with.
*
*  MAVAL   (global output) INTEGER array
*          On entry, MAVAL is an array of dimension LDVAL. On exit, this
*          array  contains  the values  of  DESCA( M_ )  to run the code
*          with.
*
*  NAVAL   (global output) INTEGER array
*          On entry, NAVAL is an array of dimension LDVAL. On exit, this
*          array  contains  the values  of  DESCA( N_ )  to run the code
*          with.
*
*  IMBAVAL (global output) INTEGER array
*          On entry,  IMBAVAL  is an array of  dimension LDVAL. On exit,
*          this  array  contains  the values of DESCA( IMB_ ) to run the
*          code with.
*
*  MBAVAL  (global output) INTEGER array
*          On entry,  MBAVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains  the values of DESCA( MB_ ) to  run the
*          code with.
*
*  INBAVAL (global output) INTEGER array
*          On entry,  INBAVAL  is an array of  dimension LDVAL. On exit,
*          this  array  contains  the values of DESCA( INB_ ) to run the
*          code with.
*
*  NBAVAL  (global output) INTEGER array
*          On entry,  NBAVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains  the values of DESCA( NB_ ) to  run the
*          code with.
*
*  RSCAVAL (global output) INTEGER array
*          On entry, RSCAVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains the values of DESCA( RSRC_ ) to run the
*          code with.
*
*  CSCAVAL (global output) INTEGER array
*          On entry, CSCAVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains the values of DESCA( CSRC_ ) to run the
*          code with.
*
*  IAVAL   (global output) INTEGER array
*          On entry, IAVAL is an array of dimension LDVAL. On exit, this
*          array  contains the values of IA to run the code with.
*
*  JAVAL   (global output) INTEGER array
*          On entry, JAVAL is an array of dimension LDVAL. On exit, this
*          array  contains the values of JA to run the code with.
*
*  MXVAL   (global output) INTEGER array
*          On entry, MXVAL is an array of dimension LDVAL. On exit, this
*          array  contains  the values  of  DESCX( M_ )  to run the code
*          with.
*
*  NXVAL   (global output) INTEGER array
*          On entry, NXVAL is an array of dimension LDVAL. On exit, this
*          array  contains  the values  of  DESCX( N_ )  to run the code
*          with.
*
*  IMBXVAL (global output) INTEGER array
*          On entry,  IMBXVAL  is an array of  dimension LDVAL. On exit,
*          this  array  contains  the values of DESCX( IMB_ ) to run the
*          code with.
*
*  MBXVAL  (global output) INTEGER array
*          On entry,  MBXVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains  the values of DESCX( MB_ ) to  run the
*          code with.
*
*  INBXVAL (global output) INTEGER array
*          On entry,  INBXVAL  is an array of  dimension LDVAL. On exit,
*          this  array  contains  the values of DESCX( INB_ ) to run the
*          code with.
*
*  NBXVAL  (global output) INTEGER array
*          On entry,  NBXVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains  the values of DESCX( NB_ ) to  run the
*          code with.
*
*  RSCXVAL (global output) INTEGER array
*          On entry, RSCXVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains the values of DESCX( RSRC_ ) to run the
*          code with.
*
*  CSCXVAL (global output) INTEGER array
*          On entry, CSCXVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains the values of DESCX( CSRC_ ) to run the
*          code with.
*
*  IXVAL   (global output) INTEGER array
*          On entry, IXVAL is an array of dimension LDVAL. On exit, this
*          array  contains the values of IX to run the code with.
*
*  JXVAL   (global output) INTEGER array
*          On entry, JXVAL is an array of dimension LDVAL. On exit, this
*          array  contains the values of JX to run the code with.
*
*  INCXVAL (global output) INTEGER array
*          On entry,  INCXVAL  is  an array of dimension LDVAL. On exit,
*          this array  contains the values of INCX to run the code with.
*
*  MYVAL   (global output) INTEGER array
*          On entry, MYVAL is an array of dimension LDVAL. On exit, this
*          array  contains  the values  of  DESCY( M_ )  to run the code
*          with.
*
*  NYVAL   (global output) INTEGER array
*          On entry, NYVAL is an array of dimension LDVAL. On exit, this
*          array  contains  the values  of  DESCY( N_ )  to run the code
*          with.
*
*  IMBYVAL (global output) INTEGER array
*          On entry,  IMBYVAL  is an array of  dimension LDVAL. On exit,
*          this  array  contains  the values of DESCY( IMB_ ) to run the
*          code with.
*
*  MBYVAL  (global output) INTEGER array
*          On entry,  MBYVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains  the values of DESCY( MB_ ) to  run the
*          code with.
*
*  INBYVAL (global output) INTEGER array
*          On entry,  INBYVAL  is an array of  dimension LDVAL. On exit,
*          this  array  contains  the values of DESCY( INB_ ) to run the
*          code with.
*
*  NBYVAL  (global output) INTEGER array
*          On entry,  NBYVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains  the values of DESCY( NB_ ) to  run the
*          code with.
*
*  RSCYVAL (global output) INTEGER array
*          On entry, RSCYVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains the values of DESCY( RSRC_ ) to run the
*          code with.
*
*  CSCYVAL (global output) INTEGER array
*          On entry, CSCYVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains the values of DESCY( CSRC_ ) to run the
*          code with.
*
*  IYVAL   (global output) INTEGER array
*          On entry, IYVAL is an array of dimension LDVAL. On exit, this
*          array  contains the values of IY to run the code with.
*
*  JYVAL   (global output) INTEGER array
*          On entry, JYVAL is an array of dimension LDVAL. On exit, this
*          array  contains the values of JY to run the code with.
*
*  INCYVAL (global output) INTEGER array
*          On entry,  INCYVAL  is  an array of dimension LDVAL. On exit,
*          this array  contains the values of INCY to run the code with.
*
*  LDVAL   (global input) INTEGER
*          On entry, LDVAL specifies the maximum number of different va-
*          lues that can be used for DIAG, TRANS, UPLO, M, N,  DESCA(:),
*          IA,  JA,  DESCX(:),  IX, JX, INCX, DESCY(:), IY, JY and INCY.
*          This is also the maximum number of test cases.
*
*  NGRIDS  (global output) INTEGER
*          On exit, NGRIDS specifies the number of different values that
*          can be used for P and Q.
*
*  PVAL    (global output) INTEGER array
*          On entry, PVAL is an array of dimension LDPVAL. On exit, this
*          array contains the values of P to run the code with.
*
*  LDPVAL  (global input) INTEGER
*          On entry,  LDPVAL  specifies  the maximum number of different
*          values that can be used for P.
*
*  QVAL    (global output) INTEGER array
*          On entry, QVAL is an array of dimension LDQVAL. On exit, this
*          array contains the values of Q to run the code with.
*
*  LDQVAL  (global input) INTEGER
*          On entry,  LDQVAL  specifies  the maximum number of different
*          values that can be used for Q.
*
*  NBLOG   (global output) INTEGER
*          On exit, NBLOG specifies the logical computational block size
*          to run the tests with. NBLOG must be at least one.
*
*  LTEST   (global output) LOGICAL array
*          On entry,  LTEST  is an array of dimension at least seven. On
*          exit, if LTEST( i ) is .TRUE., the i-th Level 2 PBLAS routine
*          will be tested.  See  the  input file for the ordering of the
*          routines.
*
*  SOF     (global output) LOGICAL
*          On exit, if SOF is .TRUE., the tester will  stop on the first
*          detected failure. Otherwise, it won't.
*
*  TEE     (global output) LOGICAL
*          On exit, if TEE is .TRUE., the tester will  perform the error
*          exit tests. These tests won't be performed otherwise.
*
*  IAM     (local input) INTEGER
*          On entry,  IAM  specifies the number of the process executing
*          this routine.
*
*  IGAP    (global output) INTEGER
*          On exit, IGAP  specifies the user-specified gap used for pad-
*          ding. IGAP must be at least zero.
*
*  IVERB   (global output) INTEGER
*          On exit,  IVERB  specifies  the output verbosity level: 0 for
*          pass/fail, 1, 2 or 3 for matrix dump on errors.
*
*  NPROCS  (global input) INTEGER
*          On entry, NPROCS specifies the total number of processes.
*
*  THRESH  (global output) REAL
*          On exit,  THRESH  specifies the threshhold value for the test
*          ratio.
*
*  ALPHA   (global output) DOUBLE PRECISION
*          On exit, ALPHA specifies the value of alpha to be used in all
*          the test cases.
*
*  BETA    (global output) DOUBLE PRECISION
*          On exit, BETA  specifies the value of beta  to be used in all
*          the test cases.
*
*  WORK    (local workspace) INTEGER array
*          On   entry,   WORK   is   an  array  of  dimension  at  least
*          MAX( 3, 2*NGRIDS+37*NMAT+NSUBS+4 )  with  NSUBS  equal  to 7.
*          This array is used to pack all output arrays in order to send
*          the information in one message.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NIN, NSUBS
      PARAMETER          ( NIN = 11, NSUBS = 7 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LTESTT
      INTEGER            I, ICTXT, J
      DOUBLE PRECISION   EPS
*     ..
*     .. Local Arrays ..
      CHARACTER*7        SNAMET
      CHARACTER*79       USRINFO
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINIT, BLACS_SETUP, DGEBR2D, DGEBS2D,
     $                   ICOPY, IGEBR2D, IGEBS2D, SGEBR2D, SGEBS2D
*ype real dble cplx zplx
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, MAX, MIN
*     ..
*     .. Common Blocks ..
      CHARACTER*7        SNAMES( NSUBS )
      COMMON             /SNAMEC/SNAMES
*     ..
*     .. Executable Statements ..
*
*     Process 0 reads the input data, broadcasts to other processes and
*     writes needed information to NOUT
*
      IF( IAM.EQ.0 ) THEN
*
*        Open file and skip data file header
*
         OPEN( NIN, FILE='PDBLAS2TST.dat', STATUS='OLD' )
         READ( NIN, FMT = * ) SUMMRY
         SUMMRY = ' '
*
*        Read in user-supplied info about machine type, compiler, etc.
*
         READ( NIN, FMT = 9999 ) USRINFO
*
*        Read name and unit number for summary output file
*
         READ( NIN, FMT = * ) SUMMRY
         READ( NIN, FMT = * ) NOUT
         IF( NOUT.NE.0 .AND. NOUT.NE.6 )
     $      OPEN( NOUT, FILE = SUMMRY, STATUS = 'UNKNOWN' )
*
*        Read and check the parameter values for the tests.
*
*        Read the flag that indicates if Stop on Failure
*
         READ( NIN, FMT = * ) SOF
*
*        Read the flag that indicates if Test Error Exits
*
         READ( NIN, FMT = * ) TEE
*
*        Read the verbosity level
*
         READ( NIN, FMT = * ) IVERB
         IF( IVERB.LT.0 .OR. IVERB.GT.3 )
     $      IVERB = 0
*
*        Read the leading dimension gap
*
         READ( NIN, FMT = * ) IGAP
         IF( IGAP.LT.0 )
     $      IGAP = 0
*
*        Read the threshold value for test ratio
*
         READ( NIN, FMT = * ) THRESH
         IF( THRESH.LT.0.0 )
     $      THRESH = 16.0
*
*        Get logical computational block size
*
         READ( NIN, FMT = * ) NBLOG
         IF( NBLOG.LT.1 )
     $      NBLOG = 32
*
*        Get number of grids
*
         READ( NIN, FMT = * ) NGRIDS
         IF( NGRIDS.LT.1 .OR. NGRIDS.GT.LDPVAL ) THEN
            WRITE( NOUT, FMT = 9998 ) 'Grids', LDPVAL
            GO TO 120
         ELSE IF( NGRIDS.GT.LDQVAL ) THEN
            WRITE( NOUT, FMT = 9998 ) 'Grids', LDQVAL
            GO TO 120
         END IF
*
*        Get values of P and Q
*
         READ( NIN, FMT = * ) ( PVAL( I ), I = 1, NGRIDS )
         READ( NIN, FMT = * ) ( QVAL( I ), I = 1, NGRIDS )
*
*        Read ALPHA, BETA
*
         READ( NIN, FMT = * ) ALPHA
         READ( NIN, FMT = * ) BETA
*
*        Read number of tests.
*
         READ( NIN, FMT = * ) NMAT
         IF( NMAT.LT.1 .OR. NMAT.GT.LDVAL ) THEN
            WRITE( NOUT, FMT = 9998 ) 'Tests', LDVAL
            GO TO 120
         ENDIF
*
*        Read in input data into arrays.
*
         READ( NIN, FMT = * ) ( UPLOVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( TRANVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( DIAGVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( MVAL( I ),     I = 1, NMAT )
         READ( NIN, FMT = * ) ( NVAL( I ),     I = 1, NMAT )
         READ( NIN, FMT = * ) ( MAVAL( I ),    I = 1, NMAT )
         READ( NIN, FMT = * ) ( NAVAL( I ),    I = 1, NMAT )
         READ( NIN, FMT = * ) ( IMBAVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( INBAVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( MBAVAL( I ),   I = 1, NMAT )
         READ( NIN, FMT = * ) ( NBAVAL( I ),   I = 1, NMAT )
         READ( NIN, FMT = * ) ( RSCAVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( CSCAVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( IAVAL( I ),    I = 1, NMAT )
         READ( NIN, FMT = * ) ( JAVAL( I ),    I = 1, NMAT )
         READ( NIN, FMT = * ) ( MXVAL( I ),    I = 1, NMAT )
         READ( NIN, FMT = * ) ( NXVAL( I ),    I = 1, NMAT )
         READ( NIN, FMT = * ) ( IMBXVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( INBXVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( MBXVAL( I ),   I = 1, NMAT )
         READ( NIN, FMT = * ) ( NBXVAL( I ),   I = 1, NMAT )
         READ( NIN, FMT = * ) ( RSCXVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( CSCXVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( IXVAL( I ),    I = 1, NMAT )
         READ( NIN, FMT = * ) ( JXVAL( I ),    I = 1, NMAT )
         READ( NIN, FMT = * ) ( INCXVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( MYVAL( I ),    I = 1, NMAT )
         READ( NIN, FMT = * ) ( NYVAL( I ),    I = 1, NMAT )
         READ( NIN, FMT = * ) ( IMBYVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( INBYVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( MBYVAL( I ),   I = 1, NMAT )
         READ( NIN, FMT = * ) ( NBYVAL( I ),   I = 1, NMAT )
         READ( NIN, FMT = * ) ( RSCYVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( CSCYVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( IYVAL( I ),    I = 1, NMAT )
         READ( NIN, FMT = * ) ( JYVAL( I ),    I = 1, NMAT )
         READ( NIN, FMT = * ) ( INCYVAL( I ),  I = 1, NMAT )
*
*        Read names of subroutines and flags which indicate
*        whether they are to be tested.
*
         DO 10 I = 1, NSUBS
            LTEST( I ) = .FALSE.
   10    CONTINUE
   20    CONTINUE
         READ( NIN, FMT = 9996, END = 50 ) SNAMET, LTESTT
         DO 30 I = 1, NSUBS
            IF( SNAMET.EQ.SNAMES( I ) )
     $         GO TO 40
   30    CONTINUE
*
         WRITE( NOUT, FMT = 9995 )SNAMET
         GO TO 120
*
   40    CONTINUE
         LTEST( I ) = LTESTT
         GO TO 20
*
   50    CONTINUE
*
*        Close input file
*
         CLOSE ( NIN )
*
*        For pvm only: if virtual machine not set up, allocate it and
*        spawn the correct number of processes.
*
         IF( NPROCS.LT.1 ) THEN
            NPROCS = 0
            DO 60 I = 1, NGRIDS
               NPROCS = MAX( NPROCS, PVAL( I )*QVAL( I ) )
   60       CONTINUE
            CALL BLACS_SETUP( IAM, NPROCS )
         END IF
*
*        Temporarily define blacs grid to include all processes so
*        information can be broadcast to all processes
*
         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS )
*
*        Compute machine epsilon
*
         EPS = PDLAMCH( ICTXT, 'eps' )
*
*        Pack information arrays and broadcast
*
         CALL SGEBS2D( ICTXT, 'All', ' ', 1, 1, THRESH, 1 )
         CALL DGEBS2D( ICTXT, 'All', ' ', 1, 1, ALPHA, 1 )
         CALL DGEBS2D( ICTXT, 'All', ' ', 1, 1, BETA, 1 )
*
         WORK( 1 ) = NGRIDS
         WORK( 2 ) = NMAT
         WORK( 3 ) = NBLOG
         CALL IGEBS2D( ICTXT, 'All', ' ', 3, 1, WORK, 3 )
*
         I = 1
         IF( SOF ) THEN
            WORK( I ) = 1
         ELSE
            WORK( I ) = 0
         END IF
         I = I + 1
         IF( TEE ) THEN
            WORK( I ) = 1
         ELSE
            WORK( I ) = 0
         END IF
         I = I + 1
         WORK( I ) = IVERB
         I = I + 1
         WORK( I ) = IGAP
         I = I + 1
         DO 70 J = 1, NMAT
            WORK( I )   = ICHAR( DIAGVAL( J ) )
            WORK( I+1 ) = ICHAR( TRANVAL( J ) )
            WORK( I+2 ) = ICHAR( UPLOVAL( J ) )
            I = I + 3
   70    CONTINUE
         CALL ICOPY( NGRIDS, PVAL,     1, WORK( I ), 1 )
         I = I + NGRIDS
         CALL ICOPY( NGRIDS, QVAL,     1, WORK( I ), 1 )
         I = I + NGRIDS
         CALL ICOPY( NMAT,   MVAL,     1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   NVAL,     1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   MAVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   NAVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   IMBAVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   INBAVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   MBAVAL,   1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   NBAVAL,   1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   RSCAVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   CSCAVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   IAVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   JAVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   MXVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   NXVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   IMBXVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   INBXVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   MBXVAL,   1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   NBXVAL,   1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   RSCXVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   CSCXVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   IXVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   JXVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   INCXVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   MYVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   NYVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   IMBYVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   INBYVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   MBYVAL,   1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   NBYVAL,   1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   RSCYVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   CSCYVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   IYVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   JYVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   INCYVAL,  1, WORK( I ), 1 )
         I = I + NMAT
*
         DO 80 J = 1, NSUBS
            IF( LTEST( J ) ) THEN
               WORK( I ) = 1
            ELSE
               WORK( I ) = 0
            END IF
            I = I + 1
   80    CONTINUE
         I = I - 1
         CALL IGEBS2D( ICTXT, 'All', ' ', I, 1, WORK, I )
*
*        regurgitate input
*
         WRITE( NOUT, FMT = 9999 ) 'Level 2 PBLAS testing program.'
         WRITE( NOUT, FMT = 9999 ) USRINFO
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'Tests of the real double precision '//
     $               'Level 2 PBLAS'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9993 ) NMAT
         WRITE( NOUT, FMT = 9979 ) NBLOG
         WRITE( NOUT, FMT = 9992 ) NGRIDS
         WRITE( NOUT, FMT = 9990 )
     $               'P', ( PVAL(I), I = 1, MIN(NGRIDS, 5) )
         IF( NGRIDS.GT.5 )
     $      WRITE( NOUT, FMT = 9991 ) ( PVAL(I), I = 6,
     $                                  MIN( 10, NGRIDS ) )
         IF( NGRIDS.GT.10 )
     $      WRITE( NOUT, FMT = 9991 ) ( PVAL(I), I = 11,
     $                                  MIN( 15, NGRIDS ) )
         IF( NGRIDS.GT.15 )
     $      WRITE( NOUT, FMT = 9991 ) ( PVAL(I), I = 16, NGRIDS )
         WRITE( NOUT, FMT = 9990 )
     $               'Q', ( QVAL(I), I = 1, MIN(NGRIDS, 5) )
         IF( NGRIDS.GT.5 )
     $      WRITE( NOUT, FMT = 9991 ) ( QVAL(I), I = 6,
     $                                  MIN( 10, NGRIDS ) )
         IF( NGRIDS.GT.10 )
     $      WRITE( NOUT, FMT = 9991 ) ( QVAL(I), I = 11,
     $                                  MIN( 15, NGRIDS ) )
         IF( NGRIDS.GT.15 )
     $      WRITE( NOUT, FMT = 9991 ) ( QVAL(I), I = 16, NGRIDS )
         WRITE( NOUT, FMT = 9988 ) SOF
         WRITE( NOUT, FMT = 9987 ) TEE
         WRITE( NOUT, FMT = 9983 ) IGAP
         WRITE( NOUT, FMT = 9986 ) IVERB
         WRITE( NOUT, FMT = 9980 ) THRESH
         WRITE( NOUT, FMT = 9982 ) ALPHA
         WRITE( NOUT, FMT = 9981 ) BETA
         IF( LTEST( 1 ) ) THEN
            WRITE( NOUT, FMT = 9985 ) SNAMES( 1 ), ' ... Yes'
         ELSE
            WRITE( NOUT, FMT = 9985 ) SNAMES( 1 ), ' ... No '
         END IF
         DO 90 I = 2, NSUBS
            IF( LTEST( I ) ) THEN
               WRITE( NOUT, FMT = 9984 ) SNAMES( I ), ' ... Yes'
            ELSE
               WRITE( NOUT, FMT = 9984 ) SNAMES( I ), ' ... No '
            END IF
   90    CONTINUE
         WRITE( NOUT, FMT = 9994 ) EPS
         WRITE( NOUT, FMT = * )
*
      ELSE
*
*        If in pvm, must participate setting up virtual machine
*
         IF( NPROCS.LT.1 )
     $      CALL BLACS_SETUP( IAM, NPROCS )
*
*        Temporarily define blacs grid to include all processes so
*        information can be broadcast to all processes
*
         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS )
*
*        Compute machine epsilon
*
         EPS = PDLAMCH( ICTXT, 'eps' )
*
         CALL SGEBR2D( ICTXT, 'All', ' ', 1, 1, THRESH, 1, 0, 0 )
         CALL DGEBR2D( ICTXT, 'All', ' ', 1, 1, ALPHA, 1, 0, 0 )
         CALL DGEBR2D( ICTXT, 'All', ' ', 1, 1, BETA, 1, 0, 0 )
*
         CALL IGEBR2D( ICTXT, 'All', ' ', 3, 1, WORK, 3, 0, 0 )
         NGRIDS = WORK( 1 )
         NMAT   = WORK( 2 )
         NBLOG  = WORK( 3 )
*
         I = 2*NGRIDS + 37*NMAT + NSUBS + 4
         CALL IGEBR2D( ICTXT, 'All', ' ', I, 1, WORK, I, 0, 0 )
*
         I = 1
         IF( WORK( I ).EQ.1 ) THEN
            SOF = .TRUE.
         ELSE
            SOF = .FALSE.
         END IF
         I = I + 1
         IF( WORK( I ).EQ.1 ) THEN
            TEE = .TRUE.
         ELSE
            TEE = .FALSE.
         END IF
         I = I + 1
         IVERB = WORK( I )
         I = I + 1
         IGAP = WORK( I )
         I = I + 1
         DO 100 J = 1, NMAT
            DIAGVAL( J )  = CHAR( WORK( I ) )
            TRANVAL( J )  = CHAR( WORK( I+1 ) )
            UPLOVAL( J )  = CHAR( WORK( I+2 ) )
            I = I + 3
  100    CONTINUE
         CALL ICOPY( NGRIDS, WORK( I ), 1, PVAL,     1 )
         I = I + NGRIDS
         CALL ICOPY( NGRIDS, WORK( I ), 1, QVAL,     1 )
         I = I + NGRIDS
         CALL ICOPY( NMAT,   WORK( I ), 1, MVAL,     1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, NVAL,     1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, MAVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, NAVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, IMBAVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, INBAVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, MBAVAL,   1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, NBAVAL,   1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, RSCAVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, CSCAVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, IAVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, JAVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, MXVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, NXVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, IMBXVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, INBXVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, MBXVAL,   1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, NBXVAL,   1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, RSCXVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, CSCXVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, IXVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, JXVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, INCXVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, MYVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, NYVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, IMBYVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, INBYVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, MBYVAL,   1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, NBYVAL,   1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, RSCYVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, CSCYVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, IYVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, JYVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, INCYVAL,  1 )
         I = I + NMAT
*
         DO 110 J = 1, NSUBS
            IF( WORK( I ).EQ.1 ) THEN
               LTEST( J ) = .TRUE.
            ELSE
               LTEST( J ) = .FALSE.
            END IF
            I = I + 1
  110    CONTINUE
*
      END IF
*
      CALL BLACS_GRIDEXIT( ICTXT )
*
      RETURN
*
  120 WRITE( NOUT, FMT = 9997 )
      CLOSE( NIN )
      IF( NOUT.NE.6 .AND. NOUT.NE.0 )
     $   CLOSE( NOUT )
      CALL BLACS_ABORT( ICTXT, 1 )
*
      STOP
*
 9999 FORMAT( A )
 9998 FORMAT( ' Number of values of ',5A, ' is less than 1 or greater ',
     $        'than ', I2 )
 9997 FORMAT( ' Illegal input in file ',40A,'.  Aborting run.' )
 9996 FORMAT( A7, L2 )
 9995 FORMAT( '  Subprogram name ', A7, ' not recognized',
     $        /' ******* TESTS ABANDONED *******' )
 9994 FORMAT( 2X, 'Relative machine precision (eps) is taken to be ',
     $        E18.6 )
 9993 FORMAT( 2X, 'Number of Tests           : ', I6 )
 9992 FORMAT( 2X, 'Number of process grids   : ', I6 )
 9991 FORMAT( 2X, '                          : ', 5I6 )
 9990 FORMAT( 2X, A1, '                         : ', 5I6 )
 9988 FORMAT( 2X, 'Stop on failure flag      : ', L6 )
 9987 FORMAT( 2X, 'Test for error exits flag : ', L6 )
 9986 FORMAT( 2X, 'Verbosity level           : ', I6 )
 9985 FORMAT( 2X, 'Routines to be tested     :      ', A, A8 )
 9984 FORMAT( 2X, '                                 ', A, A8 )
 9983 FORMAT( 2X, 'Leading dimension gap     : ', I6 )
 9982 FORMAT( 2X, 'Alpha                     : ', G16.6 )
 9981 FORMAT( 2X, 'Beta                      : ', G16.6 )
 9980 FORMAT( 2X, 'Threshold value           : ', G16.6 )
 9979 FORMAT( 2X, 'Logical block size        : ', I6 )
*
*     End of PDBLA2TSTINFO
*
      END
      SUBROUTINE PDBLAS2TSTCHKE( LTEST, INOUT, NPROCS )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            INOUT, NPROCS
*     ..
*     .. Array Arguments ..
      LOGICAL            LTEST( * )
*     ..
*
*  Purpose
*  =======
*
*  PDBLAS2TSTCHKE tests the error exits of the Level 2 PBLAS.
*
*  Arguments
*  =========
*
*  LTEST   (global input) LOGICAL array
*          On entry, LTEST is an array of dimension at least 7 (NSUBS).
*             If LTEST( 1 ) is .TRUE., PDGEMV will be tested;
*             If LTEST( 2 ) is .TRUE., PDSYMV will be tested;
*             If LTEST( 3 ) is .TRUE., PDTRMV will be tested;
*             If LTEST( 4 ) is .TRUE., PDTRSV will be tested;
*             If LTEST( 5 ) is .TRUE., PDGER  will be tested;
*             If LTEST( 6 ) is .TRUE., PDSYR  will be tested;
*             If LTEST( 7 ) is .TRUE., PDSYR2 will be tested;
*
*  INOUT   (global input) INTEGER
*          On entry,  INOUT  specifies  the unit number for output file.
*          When INOUT is 6, output to screen,  when INOUT = 0, output to
*          stderr. INOUT is only defined in process 0.
*
*  NPROCS  (global input) INTEGER
*          On entry, NPROCS specifies the total number of processes cal-
*          ling this routine.
*
*  Calling sequence encodings
*  ==========================
*
*  code Formal argument list                                Examples
*
*  11   (n,      v1,v2)                                     _SWAP, _COPY
*  12   (n,s1,   v1   )                                     _SCAL, _SCAL
*  13   (n,s1,   v1,v2)                                     _AXPY, _DOT_
*  14   (n,s1,i1,v1   )                                     _AMAX
*  15   (n,u1,   v1   )                                     _ASUM, _NRM2
*
*  21   (     trans,     m,n,s1,m1,v1,s2,v2)                _GEMV
*  22   (uplo,             n,s1,m1,v1,s2,v2)                _SYMV, _HEMV
*  23   (uplo,trans,diag,  n,   m1,v1      )                _TRMV, _TRSV
*  24   (                m,n,s1,v1,v2,m1)                   _GER_
*  25   (uplo,             n,s1,v1,   m1)                   _SYR
*  26   (uplo,             n,u1,v1,   m1)                   _HER
*  27   (uplo,             n,s1,v1,v2,m1)                   _SYR2, _HER2
*
*  31   (          transa,transb,     m,n,k,s1,m1,m2,s2,m3) _GEMM
*  32   (side,uplo,                   m,n,  s1,m1,m2,s2,m3) _SYMM, _HEMM
*  33   (     uplo,trans,               n,k,s1,m1,   s2,m3) _SYRK
*  34   (     uplo,trans,               n,k,u1,m1,   u2,m3) _HERK
*  35   (     uplo,trans,               n,k,s1,m1,m2,s2,m3) _SYR2K
*  36   (     uplo,trans,               n,k,s1,m1,m2,u2,m3) _HER2K
*  37   (                             m,n,  s1,m1,   s2,m3) _TRAN_
*  38   (side,uplo,transa,       diag,m,n,  s1,m1,m2      ) _TRMM, _TRSM
*  39   (          trans,             m,n,  s1,m1,   s2,m3) _GEADD
*  40   (     uplo,trans,             m,n,  s1,m1,   s2,m3) _TRADD
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 7 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ABRTSAV
      INTEGER            I, ICTXT, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            SCODE( NSUBS )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GET, BLACS_GRIDEXIT, BLACS_GRIDINFO,
     $                   BLACS_GRIDINIT, PDDIMEE, PDGEMV, PDGER,
     $                   PDMATEE, PDOPTEE, PDSYMV, PDSYR, PDSYR2,
     $                   PDTRMV, PDTRSV, PDVECEE
*     ..
*     .. Common Blocks ..
      LOGICAL            ABRTFLG
      INTEGER            NOUT
      CHARACTER*7        SNAMES( NSUBS )
      COMMON             /SNAMEC/SNAMES
      COMMON             /PBERRORC/NOUT, ABRTFLG
*     ..
*     .. Data Statements ..
      DATA               SCODE/21, 22, 23, 23, 24, 25, 27/
*     ..
*     .. Executable Statements ..
*
*     Temporarily define blacs grid to include all processes so
*     information can be broadcast to all processes.
*
      CALL BLACS_GET( -1, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Set ABRTFLG to FALSE so that the PBLAS error handler won't abort
*     on errors during these tests and set the output device unit for
*     it.
*
      ABRTSAV = ABRTFLG
      ABRTFLG = .FALSE.
      NOUT    = INOUT
*
*     Test PDGEMV
*
      I = 1
      IF( LTEST( I ) ) THEN
         CALL PDOPTEE( ICTXT, NOUT, PDGEMV, SCODE( I ), SNAMES( I ) )
         CALL PDDIMEE( ICTXT, NOUT, PDGEMV, SCODE( I ), SNAMES( I ) )
         CALL PDMATEE( ICTXT, NOUT, PDGEMV, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDGEMV, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDSYMV
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDOPTEE( ICTXT, NOUT, PDSYMV, SCODE( I ), SNAMES( I ) )
         CALL PDDIMEE( ICTXT, NOUT, PDSYMV, SCODE( I ), SNAMES( I ) )
         CALL PDMATEE( ICTXT, NOUT, PDSYMV, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDSYMV, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDTRMV
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDOPTEE( ICTXT, NOUT, PDTRMV, SCODE( I ), SNAMES( I ) )
         CALL PDDIMEE( ICTXT, NOUT, PDTRMV, SCODE( I ), SNAMES( I ) )
         CALL PDMATEE( ICTXT, NOUT, PDTRMV, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDTRMV, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDTRSV
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDOPTEE( ICTXT, NOUT, PDTRSV, SCODE( I ), SNAMES( I ) )
         CALL PDDIMEE( ICTXT, NOUT, PDTRSV, SCODE( I ), SNAMES( I ) )
         CALL PDMATEE( ICTXT, NOUT, PDTRSV, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDTRSV, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDGER
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDDIMEE( ICTXT, NOUT, PDGER, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDGER, SCODE( I ), SNAMES( I ) )
         CALL PDMATEE( ICTXT, NOUT, PDGER, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDSYR
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDOPTEE( ICTXT, NOUT, PDSYR, SCODE( I ), SNAMES( I ) )
         CALL PDDIMEE( ICTXT, NOUT, PDSYR, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDSYR, SCODE( I ), SNAMES( I ) )
         CALL PDMATEE( ICTXT, NOUT, PDSYR, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDSYR2
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDOPTEE( ICTXT, NOUT, PDSYR2, SCODE( I ), SNAMES( I ) )
         CALL PDDIMEE( ICTXT, NOUT, PDSYR2, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDSYR2, SCODE( I ), SNAMES( I ) )
         CALL PDMATEE( ICTXT, NOUT, PDSYR2, SCODE( I ), SNAMES( I ) )
      END IF
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $   WRITE( NOUT, FMT = 9999 )
*
      CALL BLACS_GRIDEXIT( ICTXT )
*
*     Reset ABRTFLG to the value it had before calling this routine
*
      ABRTFLG = ABRTSAV
*
 9999 FORMAT( 2X, 'Error-exit tests completed.' )
*
      RETURN
*
*     End of PDBLAS2TSTCHKE
*
      END
      SUBROUTINE PDCHKARG2( ICTXT, NOUT, SNAME, UPLO, TRANS, DIAG, M,
     $                      N, ALPHA, IA, JA, DESCA, IX, JX, DESCX,
     $                      INCX, BETA, IY, JY, DESCY, INCY, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        DIAG, TRANS, UPLO
      INTEGER            IA, ICTXT, INCX, INCY, INFO, IX, IY, JA, JX,
     $                   JY, M, N, NOUT
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      SNAME
      INTEGER            DESCA( * ), DESCX( * ), DESCY( * )
*     ..
*
*  Purpose
*  =======
*
*  PDCHKARG2 checks the input-only arguments of the Level 2 PBLAS.  When
*  INFO = 0, this routine makes a copy of its arguments (which are INPUT
*  only arguments to PBLAS routines). Otherwise, it verifies the  values
*  of these arguments against the saved copies.
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  NOUT    (global input) INTEGER
*          On entry, NOUT specifies the unit number for the output file.
*          When NOUT is 6, output to screen,  when  NOUT is 0, output to
*          stderr. NOUT is only defined for process 0.
*
*  SNAME   (global input) CHARACTER*(*)
*          On entry, SNAME specifies the subroutine  name  calling  this
*          subprogram.
*
*  UPLO    (global input) CHARACTER*1
*          On entry, UPLO specifies the UPLO option in the Level 2 PBLAS
*          operation.
*
*  TRANS   (global input) CHARACTER*1
*          On entry,  TRANS  specifies  the TRANS  option in the Level 2
*          PBLAS operation.
*
*  DIAG    (global input) CHARACTER*1
*          On entry, DIAG specifies the DIAG option in the Level 2 PBLAS
*          operation.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies  the  dimension of the submatrix ope-
*          rands.
*
*  N       (global input) INTEGER
*          On entry,  N  specifies  the  dimension of the submatrix ope-
*          rands.
*
*  ALPHA   (global input) DOUBLE PRECISION
*          On entry, ALPHA specifies the scalar alpha.
*
*  IA      (global input) INTEGER
*          On entry, IA  specifies A's global row index, which points to
*          the beginning of the submatrix sub( A ).
*
*  JA      (global input) INTEGER
*          On entry, JA  specifies A's global column index, which points
*          to the beginning of the submatrix sub( A ).
*
*  DESCA   (global and local input) INTEGER array
*          On entry, DESCA  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix A.
*
*  IX      (global input) INTEGER
*          On entry, IX  specifies X's global row index, which points to
*          the beginning of the submatrix sub( X ).
*
*  JX      (global input) INTEGER
*          On entry, JX  specifies X's global column index, which points
*          to the beginning of the submatrix sub( X ).
*
*  DESCX   (global and local input) INTEGER array
*          On entry, DESCX  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix X.
*
*  INCX    (global input) INTEGER
*          On entry,  INCX   specifies  the  global  increment  for  the
*          elements of  X.  Only two values of  INCX   are  supported in
*          this version, namely 1 and M_X. INCX  must not be zero.
*
*  BETA    (global input) DOUBLE PRECISION
*          On entry, BETA specifies the scalar beta.
*
*  IY      (global input) INTEGER
*          On entry, IY  specifies Y's global row index, which points to
*          the beginning of the submatrix sub( Y ).
*
*  JY      (global input) INTEGER
*          On entry, JY  specifies Y's global column index, which points
*          to the beginning of the submatrix sub( Y ).
*
*  DESCY   (global and local input) INTEGER array
*          On entry, DESCY  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix Y.
*
*  INCY    (global input) INTEGER
*          On entry,  INCY   specifies  the  global  increment  for  the
*          elements of  Y.  Only two values of  INCY   are  supported in
*          this version, namely 1 and M_Y. INCY  must not be zero.
*
*  INFO    (global input/global output) INTEGER
*          When INFO = 0 on entry, the values of the arguments which are
*          INPUT only arguments to a PBLAS routine are copied into  sta-
*          tic variables and INFO is unchanged on exit.  Otherwise,  the
*          values  of  the  arguments are compared against the saved co-
*          pies. In case no error has been found INFO is zero on return,
*          otherwise it is non zero.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      CHARACTER*1        DIAGREF, TRANSREF, UPLOREF
      INTEGER            I, IAREF, INCXREF, INCYREF, IXREF, IYREF,
     $                   JAREF, JXREF, JYREF, MREF, MYCOL, MYROW, NPCOL,
     $                   NPROW, NREF
      DOUBLE PRECISION   ALPHAREF, BETAREF
*     ..
*     .. Local Arrays ..
      CHARACTER*15       ARGNAME
      INTEGER            DESCAREF( DLEN_ ), DESCXREF( DLEN_ ),
     $                   DESCYREF( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGSUM2D
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Save Statements ..
      SAVE
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Check if first call. If yes, then save.
*
      IF( INFO.EQ.0 ) THEN
*
         DIAGREF = DIAG
         TRANSREF = TRANS
         UPLOREF = UPLO
         MREF = M
         NREF = N
         ALPHAREF = ALPHA
         IAREF = IA
         JAREF = JA
         DO 10 I = 1, DLEN_
            DESCAREF( I ) = DESCA( I )
   10    CONTINUE
         IXREF = IX
         JXREF = JX
         DO 20 I = 1, DLEN_
            DESCXREF( I ) = DESCX( I )
   20    CONTINUE
         INCXREF = INCX
         BETAREF = BETA
         IYREF = IY
         JYREF = JY
         DO 30 I = 1, DLEN_
            DESCYREF( I ) = DESCY( I )
   30    CONTINUE
         INCYREF = INCY
*
      ELSE
*
*        Test saved args. Return with first mismatch.
*
         ARGNAME = ' '
         IF( .NOT. LSAME( DIAG, DIAGREF ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DIAG'
         ELSE IF( .NOT. LSAME( TRANS, TRANSREF ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'TRANS'
         ELSE IF( .NOT. LSAME( UPLO, UPLOREF ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'UPLO'
         ELSE IF( M.NE.MREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'M'
         ELSE IF( N.NE.NREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'N'
         ELSE IF( ALPHA.NE.ALPHAREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'ALPHA'
         ELSE IF( IA.NE.IAREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'IA'
         ELSE IF( JA.NE.JAREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'JA'
         ELSE IF( DESCA( DTYPE_ ).NE.DESCAREF( DTYPE_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCA( DTYPE_ )'
         ELSE IF( DESCA( M_ ).NE.DESCAREF( M_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCA( M_ )'
         ELSE IF( DESCA( N_ ).NE.DESCAREF( N_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCA( N_ )'
         ELSE IF( DESCA( IMB_ ).NE.DESCAREF( IMB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCA( IMB_ )'
         ELSE IF( DESCA( INB_ ).NE.DESCAREF( INB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCA( INB_ )'
         ELSE IF( DESCA( MB_ ).NE.DESCAREF( MB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCA( MB_ )'
         ELSE IF( DESCA( NB_ ).NE.DESCAREF( NB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCA( NB_ )'
         ELSE IF( DESCA( RSRC_ ).NE.DESCAREF( RSRC_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCA( RSRC_ )'
         ELSE IF( DESCA( CSRC_ ).NE.DESCAREF( CSRC_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCA( CSRC_ )'
         ELSE IF( DESCA( CTXT_ ).NE.DESCAREF( CTXT_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCA( CTXT_ )'
         ELSE IF( DESCA( LLD_ ).NE.DESCAREF( LLD_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCA( LLD_ )'
         ELSE IF( IX.NE.IXREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'IX'
         ELSE IF( JX.NE.JXREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'JX'
         ELSE IF( DESCX( DTYPE_ ).NE.DESCXREF( DTYPE_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCX( DTYPE_ )'
         ELSE IF( DESCX( M_ ).NE.DESCXREF( M_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCX( M_ )'
         ELSE IF( DESCX( N_ ).NE.DESCXREF( N_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCX( N_ )'
         ELSE IF( DESCX( IMB_ ).NE.DESCXREF( IMB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCX( IMB_ )'
         ELSE IF( DESCX( INB_ ).NE.DESCXREF( INB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCX( INB_ )'
         ELSE IF( DESCX( MB_ ).NE.DESCXREF( MB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCX( MB_ )'
         ELSE IF( DESCX( NB_ ).NE.DESCXREF( NB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCX( NB_ )'
         ELSE IF( DESCX( RSRC_ ).NE.DESCXREF( RSRC_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCX( RSRC_ )'
         ELSE IF( DESCX( CSRC_ ).NE.DESCXREF( CSRC_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCX( CSRC_ )'
         ELSE IF( DESCX( CTXT_ ).NE.DESCXREF( CTXT_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCX( CTXT_ )'
         ELSE IF( DESCX( LLD_ ).NE.DESCXREF( LLD_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCX( LLD_ )'
         ELSE IF( INCX.NE.INCXREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'INCX'
         ELSE IF( BETA.NE.BETAREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'BETA'
         ELSE IF( IY.NE.IYREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'IY'
         ELSE IF( JY.NE.JYREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'JY'
         ELSE IF( DESCY( DTYPE_ ).NE.DESCYREF( DTYPE_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCY( DTYPE_ )'
         ELSE IF( DESCY( M_ ).NE.DESCYREF( M_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCY( M_ )'
         ELSE IF( DESCY( N_ ).NE.DESCYREF( N_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCY( N_ )'
         ELSE IF( DESCY( IMB_ ).NE.DESCYREF( IMB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCY( IMB_ )'
         ELSE IF( DESCY( INB_ ).NE.DESCYREF( INB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCY( INB_ )'
         ELSE IF( DESCY( MB_ ).NE.DESCYREF( MB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCY( MB_ )'
         ELSE IF( DESCY( NB_ ).NE.DESCYREF( NB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCY( NB_ )'
         ELSE IF( DESCY( RSRC_ ).NE.DESCYREF( RSRC_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCY( RSRC_ )'
         ELSE IF( DESCY( CSRC_ ).NE.DESCYREF( CSRC_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCY( CSRC_ )'
         ELSE IF( DESCY( CTXT_ ).NE.DESCYREF( CTXT_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCY( CTXT_ )'
         ELSE IF( DESCY( LLD_ ).NE.DESCYREF( LLD_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCY( LLD_ )'
         ELSE IF( INCY.NE.INCYREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'INCY'
         ELSE
            INFO = 0
         END IF
*
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, 0 )
*
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
*
            IF( INFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 ) ARGNAME, SNAME
            ELSE
               WRITE( NOUT, FMT = 9998 ) SNAME
            END IF
*
         END IF
*
      END IF
*
 9999 FORMAT( 2X, '   ***** Input-only parameter check: ', A,
     $        ' FAILED  changed ', A, ' *****' )
 9998 FORMAT( 2X, '   ***** Input-only parameter check: ', A,
     $        ' PASSED  *****' )
*
      RETURN
*
*     End of PDCHKARG2
*
      END
      SUBROUTINE PDBLAS2TSTCHK( ICTXT, NOUT, NROUT, UPLO, TRANS, DIAG,
     $                          M, N, ALPHA, A, PA, IA, JA, DESCA, X,
     $                          PX, IX, JX, DESCX, INCX, BETA, Y, PY,
     $                          IY, JY, DESCY, INCY, THRESH, ROGUE,
     $                          WORK, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        DIAG, TRANS, UPLO
      INTEGER            IA, ICTXT, INCX, INCY, INFO, IX, IY, JA, JX,
     $                   JY, M, N, NOUT, NROUT
      REAL THRESH
      DOUBLE PRECISION   ALPHA, BETA, ROGUE
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCX( * ), DESCY( * )
      DOUBLE PRECISION   A( * ), PA( * ), PX( * ), PY( * ), WORK( * ),
     $                   X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  PDBLAS2TSTCHK performs the computational tests of the Level 2 PBLAS.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information  required to
*  establish the  mapping  between a  matrix entry and its corresponding
*  process and memory location.
*
*  In  the  following  comments,   the character _  should  be  read  as
*  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN       EXPLANATION
*  ---------------- --------------- ------------------------------------
*  DTYPE_A (global) DESCA( DTYPE_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT_  ) The BLACS context handle, indicating
*                                   the NPROW x NPCOL BLACS process grid
*                                   A  is distributed over.  The context
*                                   itself  is  global,  but  the handle
*                                   (the integer value) may vary.
*  M_A     (global) DESCA( M_     ) The  number of rows in the distribu-
*                                   ted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N_     ) The number of columns in the distri-
*                                   buted matrix A, N_A >= 0.
*  IMB_A   (global) DESCA( IMB_   ) The number of rows of the upper left
*                                   block of the matrix A, IMB_A > 0.
*  INB_A   (global) DESCA( INB_   ) The  number  of columns of the upper
*                                   left   block   of   the   matrix  A,
*                                   INB_A > 0.
*  MB_A    (global) DESCA( MB_    ) The blocking factor used to  distri-
*                                   bute the last  M_A-IMB_A rows of  A,
*                                   MB_A > 0.
*  NB_A    (global) DESCA( NB_    ) The blocking factor used to  distri-
*                                   bute the last  N_A-INB_A  columns of
*                                   A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC_  ) The process row over which the first
*                                   row of the matrix  A is distributed,
*                                   NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC_  ) The  process  column  over which the
*                                   first  column of  A  is distributed.
*                                   NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD_   ) The  leading  dimension of the local
*                                   array  storing  the  local blocks of
*                                   the distributed matrix A,
*                                   IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                   ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_NUMROC:
*  Lr( IA, K ) = PB_NUMROC( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_NUMROC( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  NOUT    (global input) INTEGER
*          On entry, NOUT specifies the unit number for the output file.
*          When NOUT is 6, output to screen,  when  NOUT is 0, output to
*          stderr. NOUT is only defined for process 0.
*
*  NROUT   (global input) INTEGER
*          On entry,  NROUT  specifies  which  routine will be tested as
*          follows:
*             If NROUT = 1,      PDGEMV will be tested;
*             else if NROUT = 2, PDSYMV will be tested;
*             else if NROUT = 3, PDTRMV will be tested;
*             else if NROUT = 4, PDTRSV will be tested;
*             else if NROUT = 5, PDGER  will be tested;
*             else if NROUT = 6, PDSYR  will be tested;
*             else if NROUT = 7, PDSYR2 will be tested;
*
*  UPLO    (global input) CHARACTER*1
*          On entry,  UPLO  specifies  if the upper or lower part of the
*          matrix operand is to be referenced.
*
*  TRANS   (global input) CHARACTER*1
*          On entry, TRANS  specifies if the matrix  operand  A is to be
*          transposed.
*
*  DIAG    (global input) CHARACTER*1
*          On entry, DIAG specifies if the  triangular matrix operand is
*          unit or non-unit.
*
*  M       (global input) INTEGER
*          On entry, M specifies the number of rows of A.
*
*  N       (global input) INTEGER
*          On entry, N specifies the number of columns of A.
*
*  ALPHA   (global input) DOUBLE PRECISION
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (local input/local output) DOUBLE PRECISION array
*          On entry, A is an array of  dimension  (DESCA( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PA.
*
*  PA      (local input) DOUBLE PRECISION array
*          On entry, PA is an array of dimension (DESCA( LLD_ ),*). This
*          array contains the local entries of the matrix PA.
*
*  IA      (global input) INTEGER
*          On entry, IA  specifies A's global row index, which points to
*          the beginning of the submatrix sub( A ).
*
*  JA      (global input) INTEGER
*          On entry, JA  specifies A's global column index, which points
*          to the beginning of the submatrix sub( A ).
*
*  DESCA   (global and local input) INTEGER array
*          On entry, DESCA  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix A.
*
*  X       (local input/local output) DOUBLE PRECISION array
*          On entry, X is an array of  dimension  (DESCX( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PX.
*
*  PX      (local input) DOUBLE PRECISION array
*          On entry, PX is an array of dimension (DESCX( LLD_ ),*). This
*          array contains the local entries of the matrix PX.
*
*  IX      (global input) INTEGER
*          On entry, IX  specifies X's global row index, which points to
*          the beginning of the submatrix sub( X ).
*
*  JX      (global input) INTEGER
*          On entry, JX  specifies X's global column index, which points
*          to the beginning of the submatrix sub( X ).
*
*  DESCX   (global and local input) INTEGER array
*          On entry, DESCX  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix X.
*
*  INCX    (global input) INTEGER
*          On entry,  INCX   specifies  the  global  increment  for  the
*          elements of  X.  Only two values of  INCX   are  supported in
*          this version, namely 1 and M_X. INCX  must not be zero.
*
*  BETA    (global input) DOUBLE PRECISION
*          On entry, BETA specifies the scalar beta.
*
*  Y       (local input/local output) DOUBLE PRECISION array
*          On entry, Y is an array of  dimension  (DESCY( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PY.
*
*  PY      (local input) DOUBLE PRECISION array
*          On entry, PY is an array of dimension (DESCY( LLD_ ),*). This
*          array contains the local entries of the matrix PY.
*
*  IY      (global input) INTEGER
*          On entry, IY  specifies Y's global row index, which points to
*          the beginning of the submatrix sub( Y ).
*
*  JY      (global input) INTEGER
*          On entry, JY  specifies Y's global column index, which points
*          to the beginning of the submatrix sub( Y ).
*
*  DESCY   (global and local input) INTEGER array
*          On entry, DESCY  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix Y.
*
*  INCY    (global input) INTEGER
*          On entry,  INCY   specifies  the  global  increment  for  the
*          elements of  Y.  Only two values of  INCY   are  supported in
*          this version, namely 1 and M_Y. INCY  must not be zero.
*
*  THRESH  (global input) REAL
*          On entry, THRESH is the threshold value for the test ratio.
*
*  ROGUE   (global input) DOUBLE PRECISION
*          On entry,  ROGUE  specifies  the  constant  used  to  pad the
*          non-referenced part of triangular or symmetric matrices.
*
*  WORK    (workspace) DOUBLE PRECISION array
*          On entry, WORK  is an array of dimension LWORK where LWORK is
*          at least  MAX( M, N ). This array is used to store the compu-
*          ted gauges (see PDMVCH).
*
*  INFO    (global output) INTEGER
*          On exit, if INFO = 0,  no  error  has  been  found, otherwise
*          if( MOD( INFO,   2 ) = 1 ) then an error on A has been found,
*          if( MOD( INFO/2, 2 ) = 1 ) then an error on X has been found,
*          if( MOD( INFO/4, 2 ) = 1 ) then an error on Y has been found.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, MYCOL, MYROW, NPCOL, NPROW
      DOUBLE PRECISION   ERR
*     ..
*     .. Local Arrays ..
      INTEGER            IERR( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DTRSV, PB_DLASET, PDCHKMIN,
     $                   PDCHKVIN, PDMVCH, PDTRMV, PDVMCH, PDVMCH2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( ( M.LE.0 ).OR.( N.LE.0 ) )
     $   RETURN
*
*     Start the operations
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      DO 10 I = 1, 3
         IERR( I ) = 0
   10 CONTINUE
*
      IF( NROUT.EQ.1 ) THEN
*
*        Test PDGEMV
*
*        Check the resulting vector Y
*
         CALL PDMVCH( ICTXT, TRANS, M, N, ALPHA, A, IA, JA, DESCA, X,
     $                IX, JX, DESCX, INCX, BETA, Y, PY, IY, JY, DESCY,
     $                INCY, WORK, ERR, IERR( 3 ) )
*
         IF( IERR( 3 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9996 ) ERR
         END IF
*
*        Check the input-only arguments
*
         CALL PDCHKMIN( ERR, M, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
         IF( LSAME( TRANS, 'N' ) ) THEN
            CALL PDCHKVIN( ERR, N, X, PX, IX, JX, DESCX, INCX,
     $                     IERR( 2 ) )
         ELSE
            CALL PDCHKVIN( ERR, M, X, PX, IX, JX, DESCX, INCX,
     $                     IERR( 2 ) )
         END IF
*
      ELSE IF( NROUT.EQ.2 ) THEN
*
*        Test PDSYMV
*
*        Check the resulting vector Y
*
         CALL PDMVCH( ICTXT, 'No transpose', N, N, ALPHA, A, IA, JA,
     $                DESCA, X, IX, JX, DESCX, INCX, BETA, Y, PY, IY,
     $                JY, DESCY, INCY, WORK, ERR, IERR( 3 ) )
*
         IF( IERR( 3 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9996 ) ERR
         END IF
*
*        Check the input-only arguments
*
         IF( LSAME( UPLO, 'L' ) ) THEN
            CALL PB_DLASET( 'Upper', N-1, N-1, 0, ROGUE, ROGUE,
     $                      A( IA+JA*DESCA( M_ ) ), DESCA( M_ ) )
         ELSE
            CALL PB_DLASET( 'Lower', N-1, N-1, 0, ROGUE, ROGUE,
     $                      A( IA+1+(JA-1)*DESCA( M_ ) ), DESCA( M_ ) )
         END IF
         CALL PDCHKMIN( ERR, N, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
         CALL PDCHKVIN( ERR, N, X, PX, IX, JX, DESCX, INCX, IERR( 2 ) )
*
      ELSE IF( NROUT.EQ.3 ) THEN
*
*        Test PDTRMV
*
*        Check the resulting vector X
*
         CALL PDMVCH( ICTXT, TRANS, N, N, ONE, A, IA, JA, DESCA, Y, IX,
     $                JX, DESCX, INCX, ZERO, X, PX, IX, JX, DESCX, INCX,
     $                WORK, ERR, IERR( 2 ) )
*
         IF( IERR( 2 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9996 ) ERR
         END IF
*
*        Check the input-only arguments
*
         IF( LSAME( UPLO, 'L' ) ) THEN
            IF( LSAME( DIAG, 'N' ) ) THEN
               CALL PB_DLASET( 'Upper', N-1, N-1, 0, ROGUE, ROGUE,
     $                         A( IA+JA*DESCA( M_ ) ), DESCA( M_ ) )
            ELSE
               CALL PB_DLASET( 'Upper', N, N, 0, ROGUE, ONE,
     $                         A( IA+(JA-1)*DESCA( M_ ) ), DESCA( M_ ) )
            END IF
         ELSE
            IF( LSAME( DIAG, 'N' ) ) THEN
               CALL PB_DLASET( 'Lower', N-1, N-1, 0, ROGUE, ROGUE,
     $                         A( IA+1+(JA-1)*DESCA( M_ ) ),
     $                         DESCA( M_ ) )
            ELSE
               CALL PB_DLASET( 'Lower', N, N, 0, ROGUE, ONE,
     $                         A( IA+(JA-1)*DESCA( M_ ) ), DESCA( M_ ) )
            END IF
         END IF
         CALL PDCHKMIN( ERR, N, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
*
      ELSE IF( NROUT.EQ.4 ) THEN
*
*        Test PDTRSV
*
*        Check the resulting vector X
*
         CALL DTRSV( UPLO, TRANS, DIAG, N, A( IA+(JA-1)*DESCA( M_ ) ),
     $               DESCA( M_ ), X( IX+(JX-1)*DESCX( M_ ) ), INCX )
         CALL PDTRMV( UPLO, TRANS, DIAG, N, PA, IA, JA, DESCA, PX, IX,
     $                JX, DESCX, INCX )
         CALL PDMVCH( ICTXT, TRANS, N, N, ONE, A, IA, JA, DESCA, X, IX,
     $                JX, DESCX, INCX, ZERO, Y, PX, IX, JX, DESCX, INCX,
     $                WORK, ERR, IERR( 2 ) )
*
         IF( IERR( 2 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9996 ) ERR
         END IF
*
*        Check the input-only arguments
*
         IF( LSAME( UPLO, 'L' ) ) THEN
            IF( LSAME( DIAG, 'N' ) ) THEN
               CALL PB_DLASET( 'Upper', N-1, N-1, 0, ROGUE, ROGUE,
     $                         A( IA+JA*DESCA( M_ ) ), DESCA( M_ ) )
            ELSE
               CALL PB_DLASET( 'Upper', N, N, 0, ROGUE, ONE,
     $                         A( IA+(JA-1)*DESCA( M_ ) ), DESCA( M_ ) )
            END IF
         ELSE
            IF( LSAME( DIAG, 'N' ) ) THEN
               CALL PB_DLASET( 'Lower', N-1, N-1, 0, ROGUE, ROGUE,
     $                         A( IA+1+(JA-1)*DESCA( M_ ) ),
     $                         DESCA( M_ ) )
            ELSE
               CALL PB_DLASET( 'Lower', N, N, 0, ROGUE, ONE,
     $                         A( IA+(JA-1)*DESCA( M_ ) ), DESCA( M_ ) )
            END IF
         END IF
         CALL PDCHKMIN( ERR, N, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
*
      ELSE IF( NROUT.EQ.5 ) THEN
*
*        Test PDGER
*
*        Check the resulting matrix A
*
         CALL PDVMCH( ICTXT, 'Ge', M, N, ALPHA, X, IX, JX, DESCX,
     $                INCX, Y, IY, JY, DESCY, INCY, A, PA, IA, JA,
     $                DESCA, WORK, ERR, IERR( 1 ) )
         IF( IERR( 1 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9996 ) ERR
         END IF
*
*        Check the input-only arguments
*
         CALL PDCHKVIN( ERR, M, X, PX, IX, JX, DESCX, INCX, IERR( 2 ) )
         CALL PDCHKVIN( ERR, N, Y, PY, IY, JY, DESCY, INCY, IERR( 3 ) )
*
      ELSE IF( NROUT.EQ.6 ) THEN
*
*        Test PDSYR
*
*        Check the resulting matrix A
*
         CALL PDVMCH( ICTXT, UPLO, N, N, ALPHA, X, IX, JX, DESCX,
     $                INCX, X, IX, JX, DESCX, INCX, A, PA, IA, JA,
     $                DESCA, WORK, ERR, IERR( 1 ) )
         IF( IERR( 1 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9996 ) ERR
         END IF
*
*        Check the input-only arguments
*
         CALL PDCHKVIN( ERR, N, X, PX, IX, JX, DESCX, INCX, IERR( 2 ) )
*
      ELSE IF( NROUT.EQ.7 ) THEN
*
*        Test PDSYR2
*
*        Check the resulting matrix A
*
         CALL PDVMCH2( ICTXT, UPLO, N, N, ALPHA, X, IX, JX, DESCX, INCX,
     $                 Y, IY, JY, DESCY, INCY, A, PA, IA, JA, DESCA,
     $                 WORK, ERR, IERR( 1 ) )
         IF( IERR( 1 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9996 ) ERR
         END IF
*
*        Check the input-only arguments
*
         CALL PDCHKVIN( ERR, N, X, PX, IX, JX, DESCX, INCX, IERR( 2 ) )
         CALL PDCHKVIN( ERR, N, Y, PY, IY, JY, DESCY, INCY, IERR( 3 ) )
*
      END IF
*
      IF( IERR( 1 ).NE.0 ) THEN
         INFO = INFO + 1
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9999 ) 'A'
      END IF
*
      IF( IERR( 2 ).NE.0 ) THEN
         INFO = INFO + 2
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9998 ) 'X'
      END IF
*
      IF( IERR( 3 ).NE.0 ) THEN
         INFO = INFO + 4
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9998 ) 'Y'
      END IF
*
 9999 FORMAT( 2X, '   ***** ERROR: Matrix operand ', A,
     $        ' is incorrect.' )
 9998 FORMAT( 2X, '   ***** ERROR: Vector operand ', A,
     $        ' is incorrect.' )
 9997 FORMAT( 2X, '   ***** FATAL ERROR - Computed result is less ',
     $        'than half accurate *****' )
 9996 FORMAT( 2X, '   ***** Test completed with maximum test ratio: ',
     $        F11.5, ' SUSPECT *****' )
*
      RETURN
*
*     End of PDBLAS2TSTCHK
*
      END
