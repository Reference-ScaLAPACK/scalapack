      BLOCK DATA
      INTEGER NSUBS
      PARAMETER (NSUBS = 11)
      CHARACTER*7        SNAMES( NSUBS )
      COMMON             /SNAMEC/SNAMES
      DATA               SNAMES/'PZGEMM ', 'PZSYMM ', 'PZHEMM ',
     $                   'PZSYRK ', 'PZHERK ', 'PZSYR2K',
     $                   'PZHER2K', 'PZTRMM ', 'PZTRSM ',
     $                   'PZGEADD', 'PZTRADD'/
      END BLOCK DATA
                   
      PROGRAM PZBLA3TST
*
*  -- PBLAS testing driver (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*  Purpose
*  =======
*
*  PZBLA3TST is the main testing program for the Level 3 PBLAS routines.
*
*  The program must be driven by a short data file.  An  annotated exam-
*  ple of a data file can be obtained by deleting the first 3 characters
*
*  from the following 64 lines:
*  'Level 3 PBLAS, Testing input file'
*  'Intel iPSC/860 hypercube, gamma model.'
*  'PZBLAS3TST.SUMM'     output file name (if any)
*  6     device out
*  F     logical flag, T to stop on failures
*  F     logical flag, T to test error exits
*  0     verbosity, 0 for pass/fail, 1-3 for matrix dump on errors
*  10    the leading dimension gap
*  16.0  threshold value of test ratio
*  10              value of the logical computational blocksize NB
*  1               number of process grids (ordered pairs of P & Q)
*  2 2 1 4 2 3 8   values of P
*  2 2 4 1 3 2 1   values of Q
*  (1.0D0, 0.0D0)  value of ALPHA
*  (1.0D0, 0.0D0)  value of BETA
*  2               number of tests problems
*  'N' 'U'         values of DIAG
*  'L' 'R'         values of SIDE
*  'N' 'T'         values of TRANSA
*  'N' 'T'         values of TRANSB
*  'U' 'L'         values of UPLO
*  3  4            values of M
*  3  4            values of N
*  3  4            values of K
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
*  6 10            values of M_B
*  6 10            values of N_B
*  2  5            values of IMB_B
*  2  5            values of INB_B
*  2  5            values of MB_B
*  2  5            values of NB_B
*  0  1            values of RSRC_B
*  0  0            values of CSRC_B
*  1  1            values of IB
*  1  1            values of JB
*  6 10            values of M_C
*  6 10            values of N_C
*  2  5            values of IMB_C
*  2  5            values of INB_C
*  2  5            values of MB_C
*  2  5            values of NB_C
*  0  1            values of RSRC_C
*  0  0            values of CSRC_C
*  1  1            values of IC
*  1  1            values of JC
*  PZGEMM  T  put F for no test in the same column
*  PZSYMM  T  put F for no test in the same column
*  PZHEMM  T  put F for no test in the same column
*  PZSYRK  T  put F for no test in the same column
*  PZHERK  T  put F for no test in the same column
*  PZSYR2K T  put F for no test in the same column
*  PZHER2K T  put F for no test in the same column
*  PZTRMM  T  put F for no test in the same column
*  PZTRSM  T  put F for no test in the same column
*  PZGEADD T  put F for no test in the same column
*  PZTRADD T  put F for no test in the same column
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
*  ZPLXSZ  INTEGER
*          DBLESZ  and  ZPLXSZ indicate the length in bytes on the given
*          platform  for a double  precision real and a double precision
*          complex. By default, DBLESZ is set to eight and ZPLXSZ is set
*          to sixteen.
*
*  MEM     COMPLEX*16 array
*          MEM is an array of dimension TOTMEM / ZPLXSZ.
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
      INTEGER            MAXTESTS, MAXGRIDS, GAPMUL, ZPLXSZ, TOTMEM,
     $                   MEMSIZ, NSUBS, DBLESZ
      COMPLEX*16         ONE, PADVAL, ZERO, ROGUE
      PARAMETER          ( MAXTESTS = 20, MAXGRIDS = 20, GAPMUL = 10,
     $                   ZPLXSZ = 16, TOTMEM = 2000000,
     $                   MEMSIZ = TOTMEM / ZPLXSZ, DBLESZ = 8,
     $                   PADVAL = ( -9923.0D+0, -9923.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ROGUE = ( -1.0D+10, 1.0D+10 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ), NSUBS = 11 )
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
      CHARACTER*1        ADIAGDO, AFORM, CFORM, DIAG, SIDE, TRANSA,
     $                   TRANSB, UPLO
      INTEGER            CSRCA, CSRCB, CSRCC, I, IA, IAM, IASEED, IB,
     $                   IBSEED, IC, ICSEED, ICTXT, IGAP, IMBA, IMBB,
     $                   IMBC, IMIDA, IMIDB, IMIDC, INBA, INBB, INBC,
     $                   IPA, IPB, IPC, IPG, IPMATA, IPMATB, IPMATC,
     $                   IPOSTA, IPOSTB, IPOSTC, IPREA, IPREB, IPREC,
     $                   IPW, IVERB, J, JA, JB, JC, K, L, LDA, LDB, LDC,
     $                   M, MA, MB, MBA, MBB, MBC, MC, MEMREQD, MPA,
     $                   MPB, MPC, MYCOL, MYROW, N, NA, NB, NBA, NBB,
     $                   NBC, NC, NCOLA, NCOLB, NCOLC, NGRIDS, NOUT,
     $                   NPCOL, NPROCS, NPROW, NQA, NQB, NQC, NROWA,
     $                   NROWB, NROWC, NTESTS, OFFDA, OFFDC, RSRCA,
     $                   RSRCB, RSRCC, TSKIP, TSTCNT
      REAL               THRESH
      COMPLEX*16         ALPHA, BETA, SCALE
*     ..
*     .. Local Arrays ..
      LOGICAL            BCHECK( NSUBS ), CCHECK( NSUBS ),
     $                   LTEST( NSUBS )
      CHARACTER*1        DIAGVAL( MAXTESTS ), SIDEVAL( MAXTESTS ),
     $                   TRNAVAL( MAXTESTS ), TRNBVAL( MAXTESTS ),
     $                   UPLOVAL( MAXTESTS )
      CHARACTER*80       OUTFILE
      INTEGER            CSCAVAL( MAXTESTS ), CSCBVAL( MAXTESTS ),
     $                   CSCCVAL( MAXTESTS ), DESCA( DLEN_ ),
     $                   DESCAR( DLEN_ ), DESCB( DLEN_ ),
     $                   DESCBR( DLEN_ ), DESCC( DLEN_ ),
     $                   DESCCR( DLEN_ ), IAVAL( MAXTESTS ),
     $                   IBVAL( MAXTESTS ), ICVAL( MAXTESTS ),
     $                   IERR( 6 ), IMBAVAL( MAXTESTS ),
     $                   IMBBVAL( MAXTESTS ), IMBCVAL( MAXTESTS ),
     $                   INBAVAL( MAXTESTS ), INBBVAL( MAXTESTS ),
     $                   INBCVAL( MAXTESTS ), JAVAL( MAXTESTS ),
     $                   JBVAL( MAXTESTS ), JCVAL( MAXTESTS )
      INTEGER            KFAIL( NSUBS ), KPASS( NSUBS ), KSKIP( NSUBS ),
     $                   KTESTS( NSUBS ), KVAL( MAXTESTS ),
     $                   MAVAL( MAXTESTS ), MBAVAL( MAXTESTS ),
     $                   MBBVAL( MAXTESTS ), MBCVAL( MAXTESTS ),
     $                   MBVAL( MAXTESTS ), MCVAL( MAXTESTS ),
     $                   MVAL( MAXTESTS ), NAVAL( MAXTESTS ),
     $                   NBAVAL( MAXTESTS ), NBBVAL( MAXTESTS ),
     $                   NBCVAL( MAXTESTS ), NBVAL( MAXTESTS ),
     $                   NCVAL( MAXTESTS ), NVAL( MAXTESTS ),
     $                   PVAL( MAXTESTS ), QVAL( MAXTESTS ),
     $                   RSCAVAL( MAXTESTS ), RSCBVAL( MAXTESTS ),
     $                   RSCCVAL( MAXTESTS )
      COMPLEX*16         MEM( MEMSIZ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $                   IGSUM2D, PB_DESCSET2, PB_PZLAPRNT, PB_ZCHEKPAD,
     $                   PB_ZFILLPAD, PB_ZLASCAL, PB_ZLASET, PMDESCCHK,
     $                   PMDIMCHK, PZBLA3TSTINFO, PZBLAS3TSTCHK,
     $                   PZBLAS3TSTCHKE, PZCHKARG3, PZCHKMOUT, PZGEADD,
     $                   PZGEMM, PZHEMM, PZHER2K, PZHERK, PZIPSET,
     $                   PZLAGEN, PZLASCAL, PZLASET, PZMPRNT, PZSYMM,
     $                   PZSYR2K, PZSYRK, PZTRADD, PZTRMM, PZTRSM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            PB_FCEIL
      EXTERNAL           PB_FCEIL, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, MAX, MOD, REAL
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
      DATA               BCHECK/.TRUE., .TRUE., .TRUE., .FALSE.,
     $                   .FALSE., .TRUE., .TRUE., .TRUE., .TRUE.,
     $                   .FALSE., .FALSE./
      DATA               CCHECK/.TRUE., .TRUE., .TRUE., .TRUE., .TRUE.,
     $                   .TRUE., .TRUE., .FALSE., .FALSE., .TRUE.,
     $                   .TRUE./
*     ..
*     .. Executable Statements ..
*
*     Initialization
*
*     Set flag so that the PBLAS error handler won't abort on errors,
*     so that the tester will detect unsupported operations.
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
      IBSEED = 200
      ICSEED = 300
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
      CALL PZBLA3TSTINFO( OUTFILE, NOUT, NTESTS, DIAGVAL, SIDEVAL,
     $                    TRNAVAL, TRNBVAL, UPLOVAL, MVAL, NVAL,
     $                    KVAL, MAVAL, NAVAL, IMBAVAL, MBAVAL,
     $                    INBAVAL, NBAVAL, RSCAVAL, CSCAVAL, IAVAL,
     $                    JAVAL, MBVAL, NBVAL, IMBBVAL, MBBVAL,
     $                    INBBVAL, NBBVAL, RSCBVAL, CSCBVAL, IBVAL,
     $                    JBVAL, MCVAL, NCVAL, IMBCVAL, MBCVAL,
     $                    INBCVAL, NBCVAL, RSCCVAL, CSCCVAL, ICVAL,
     $                    JCVAL, MAXTESTS, NGRIDS, PVAL, MAXGRIDS,
     $                    QVAL, MAXGRIDS, NBLOG, LTEST, SOF, TEE, IAM,
     $                    IGAP, IVERB, NPROCS, THRESH, ALPHA, BETA,
     $                    MEM )
*
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9976 )
         WRITE( NOUT, FMT = * )
      END IF
*
*     If TEE is set then Test Error Exits of routines.
*
      IF( TEE )
     $   CALL PZBLAS3TSTCHKE( LTEST, NOUT, NPROCS )
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
            DIAG   = DIAGVAL( J )
            SIDE   = SIDEVAL( J )
            TRANSA = TRNAVAL( J )
            TRANSB = TRNBVAL( J )
            UPLO   = UPLOVAL( J )
*
            M      = MVAL( J )
            N      = NVAL( J )
            K      = KVAL( J )
*
            MA    = MAVAL( J )
            NA    = NAVAL( J )
            IMBA  = IMBAVAL( J )
            MBA   = MBAVAL( J )
            INBA  = INBAVAL( J )
            NBA   = NBAVAL( J )
            RSRCA = RSCAVAL( J )
            CSRCA = CSCAVAL( J )
            IA    = IAVAL( J )
            JA    = JAVAL( J )
*
            MB    = MBVAL( J )
            NB    = NBVAL( J )
            IMBB  = IMBBVAL( J )
            MBB   = MBBVAL( J )
            INBB  = INBBVAL( J )
            NBB   = NBBVAL( J )
            RSRCB = RSCBVAL( J )
            CSRCB = CSCBVAL( J )
            IB    = IBVAL( J )
            JB    = JBVAL( J )
*
            MC    = MCVAL( J )
            NC    = NCVAL( J )
            IMBC  = IMBCVAL( J )
            MBC   = MBCVAL( J )
            INBC  = INBCVAL( J )
            NBC   = NBCVAL( J )
            RSRCC = RSCCVAL( J )
            CSRCC = CSCCVAL( J )
            IC    = ICVAL( J )
            JC    = JCVAL( J )
*
            IF( IAM.EQ.0 ) THEN
*
               TSTCNT = TSTCNT + 1
*
               WRITE( NOUT, FMT = * )
               WRITE( NOUT, FMT = 9996 ) TSTCNT, NPROW, NPCOL
               WRITE( NOUT, FMT = * )
*
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9994 )
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9993 ) M, N, K, SIDE, UPLO, TRANSA,
     $                                   TRANSB, DIAG
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
               WRITE( NOUT, FMT = 9991 ) IB, JB, MB, NB, IMBB, INBB,
     $                                   MBB, NBB, RSRCB, CSRCB
*
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9989 )
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9991 ) IC, JC, MC, NC, IMBC, INBC,
     $                                   MBC, NBC, RSRCC, CSRCC
*
               WRITE( NOUT, FMT = 9995 )
*
            END IF
*
*           Check the validity of the input test parameters
*
            IF( .NOT.LSAME( SIDE, 'L' ).AND.
     $          .NOT.LSAME( SIDE, 'R' ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'SIDE'
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
            IF( .NOT.LSAME( UPLO, 'U' ).AND.
     $          .NOT.LSAME( UPLO, 'L' ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'UPLO'
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
            IF( .NOT.LSAME( TRANSA, 'N' ).AND.
     $          .NOT.LSAME( TRANSA, 'T' ).AND.
     $          .NOT.LSAME( TRANSA, 'C' ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'TRANSA'
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
            IF( .NOT.LSAME( TRANSB, 'N' ).AND.
     $          .NOT.LSAME( TRANSB, 'T' ).AND.
     $          .NOT.LSAME( TRANSB, 'C' ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'TRANSB'
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
            IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $          .NOT.LSAME( DIAG , 'N' ) )THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'DIAG'
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
*
            CALL PMDESCCHK( ICTXT, NOUT, 'B', DESCB,
     $                      BLOCK_CYCLIC_2D_INB, MB, NB, IMBB, INBB,
     $                      MBB, NBB, RSRCB, CSRCB, MPB, NQB, IPREB,
     $                      IMIDB, IPOSTB, IGAP, GAPMUL, IERR( 2 ) )
*
            CALL PMDESCCHK( ICTXT, NOUT, 'C', DESCC,
     $                      BLOCK_CYCLIC_2D_INB, MC, NC, IMBC, INBC,
     $                      MBC, NBC, RSRCC, CSRCC, MPC, NQC, IPREC,
     $                      IMIDC, IPOSTC, IGAP, GAPMUL, IERR( 3 ) )
*
            IF( IERR( 1 ).GT.0 .OR. IERR( 2 ).GT.0 .OR.
     $          IERR( 3 ).GT.0 ) THEN
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
            LDA = MAX( 1, MA )
            LDB = MAX( 1, MB )
            LDC = MAX( 1, MC )
*
*           Assign pointers into MEM for matrices corresponding to
*           the distributed matrices A, X and Y.
*
            IPA = IPREA + 1
            IPB = IPA + DESCA( LLD_ )*NQA + IPOSTA + IPREB
            IPC = IPB + DESCB( LLD_ )*NQB + IPOSTB + IPREC
            IPMATA = IPC + DESCC( LLD_ )*NQC + IPOSTC
            IPMATB = IPMATA + MA*NA
            IPMATC = IPMATB + MB*NB
            IPG = IPMATC + MAX( MB*NB, MC*NC )
*
*           Check if sufficient memory.
*           Requirement = mem for local part of parallel matrices +
*                         mem for whole matrices for comp. check +
*                         mem for recving comp. check error vals.
*
            IPW = IPG + MAX( MAX( MAX( IMBA, MBA ),
     $                            MAX( IMBB, MBB ) ),
     $                       MAX( IMBC, MBC ) ) + MAX( M, MAX( N, K ) )
            MEMREQD = IPW + PB_FCEIL( REAL( MAX( M, MAX( N, K ) ) ) *
     $                REAL( DBLESZ ), REAL( ZPLXSZ ) ) - 1
            IERR( 1 ) = 0
            IF( MEMREQD.GT.MEMSIZ ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9987 ) MEMREQD*ZPLXSZ
               IERR( 1 ) = 1
            END IF
*
*           Check all processes for an error
*
            CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
            IF( IERR( 1 ).GT.0 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9988 )
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
*           Loop over all PBLAS 3 routines
*
            DO 30 L = 1, NSUBS
*
*              Continue only if this subroutine has to be tested.
*
               IF( .NOT.LTEST( L ) )
     $            GO TO 30
*
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = 9986 ) SNAMES( L )
               END IF
*
*              Define the size of the operands
*
               IF( L.EQ.1 ) THEN
*
*                 PZGEMM
*
                  NROWC = M
                  NCOLC = N
                  IF( LSAME( TRANSA, 'N' ) ) THEN
                     NROWA = M
                     NCOLA = K
                  ELSE
                     NROWA = K
                     NCOLA = M
                  END IF
                  IF( LSAME( TRANSB, 'N' ) ) THEN
                     NROWB = K
                     NCOLB = N
                  ELSE
                     NROWB = N
                     NCOLB = K
                  END IF
*
               ELSE IF( L.EQ.2 .OR. L.EQ.3 ) THEN
*
*                 PZSYMM, PZHEMM
*
                  NROWC = M
                  NCOLC = N
                  NROWB = M
                  NCOLB = N
                  IF( LSAME( SIDE, 'L' ) ) THEN
                     NROWA = M
                     NCOLA = M
                  ELSE
                     NROWA = N
                     NCOLA = N
                  END IF
*
               ELSE IF( L.EQ.4 .OR. L.EQ.5 ) THEN
*
*                 PZSYRK, PZHERK
*
                  NROWC = N
                  NCOLC = N
                  IF( LSAME( TRANSA, 'N' ) ) THEN
                     NROWA = N
                     NCOLA = K
                  ELSE
                     NROWA = K
                     NCOLA = N
                  END IF
                  NROWB = 0
                  NCOLB = 0
*
               ELSE IF( L.EQ.6 .OR. L.EQ.7 ) THEN
*
*                 PZSYR2K, PZHER2K
*
                  NROWC = N
                  NCOLC = N
                  IF( LSAME( TRANSA, 'N' ) ) THEN
                     NROWA = N
                     NCOLA = K
                     NROWB = N
                     NCOLB = K
                  ELSE
                     NROWA = K
                     NCOLA = N
                     NROWB = K
                     NCOLB = N
                  END IF
*
               ELSE IF( L.EQ.8 .OR. L.EQ.9 ) THEN
                  NROWB = M
                  NCOLB = N
                  IF( LSAME( SIDE, 'L' ) ) THEN
                     NROWA = M
                     NCOLA = M
                  ELSE
                     NROWA = N
                     NCOLA = N
                  END IF
                  NROWC = 0
                  NCOLC = 0
*
               ELSE IF( L.EQ.10 .OR. L.EQ.11 ) THEN
*
*                 PZGEADD, PZTRADD
*
                  IF( LSAME( TRANSA, 'N' ) ) THEN
                     NROWA = M
                     NCOLA = N
                  ELSE
                     NROWA = N
                     NCOLA = M
                  END IF
                  NROWC = M
                  NCOLC = N
                  NROWB = 0
                  NCOLB = 0
*
               END IF
*
*              Check the validity of the operand sizes
*
               CALL PMDIMCHK( ICTXT, NOUT, NROWA, NCOLA, 'A', IA, JA,
     $                        DESCA, IERR( 1 ) )
               CALL PMDIMCHK( ICTXT, NOUT, NROWB, NCOLB, 'B', IB, JB,
     $                        DESCB, IERR( 2 ) )
               CALL PMDIMCHK( ICTXT, NOUT, NROWC, NCOLC, 'C', IC, JC,
     $                        DESCC, IERR( 3 ) )
*
               IF( IERR( 1 ).NE.0 .OR. IERR( 2 ).NE.0 .OR.
     $             IERR( 3 ).NE.0 ) THEN
                  KSKIP( L ) = KSKIP( L ) + 1
                  GO TO 30
               END IF
*
*              Check special values of TRANSA for symmetric and
*              hermitian rank-k and rank-2k updates.
*
               IF( L.EQ.4 .OR. L.EQ.6 ) THEN
                  IF( .NOT.LSAME( TRANSA, 'N' ).AND.
     $                .NOT.LSAME( TRANSA, 'T' ) ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9975 ) 'TRANSA'
                     KSKIP( L ) = KSKIP( L ) + 1
                     GO TO 30
                  END IF
               ELSE IF( L.EQ.5 .OR. L.EQ.7 ) THEN
                  IF( .NOT.LSAME( TRANSA, 'N' ).AND.
     $                .NOT.LSAME( TRANSA, 'C' ) ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9975 ) 'TRANSA'
                     KSKIP( L ) = KSKIP( L ) + 1
                     GO TO 30
                  END IF
               END IF
*
*              Generate distributed matrices A, B and C
*
               IF( L.EQ.2 ) THEN
*
*                 PZSYMM
*
                  AFORM   = 'S'
                  ADIAGDO = 'N'
                  OFFDA   = IA - JA
                  CFORM   = 'N'
                  OFFDC   = 0
*
               ELSE IF( L.EQ.3 ) THEN
*
*                 PZHEMM
*
                  AFORM   = 'H'
                  ADIAGDO = 'N'
                  OFFDA   = IA - JA
                  CFORM   = 'N'
                  OFFDC   = 0
*
               ELSE IF( L.EQ.4 .OR. L.EQ.6 ) THEN
*
*                 PZSYRK, PZSYR2K
*
                  AFORM   = 'N'
                  ADIAGDO = 'N'
                  OFFDA   = 0
                  CFORM   = 'S'
                  OFFDC   = IC - JC
*
               ELSE IF( L.EQ.5 .OR. L.EQ.7 ) THEN
*
*                 PZHERK, PZHER2K
*
                  AFORM   = 'N'
                  ADIAGDO = 'N'
                  OFFDA   = 0
                  CFORM   = 'H'
                  OFFDC   = IC - JC
*
               ELSE IF( ( L.EQ.9 ).AND.( LSAME( DIAG, 'N' ) ) ) THEN
*
*                 PZTRSM
*
                  AFORM   = 'N'
                  ADIAGDO = 'D'
                  OFFDA   = IA - JA
                  CFORM   = 'N'
                  OFFDC   = 0
*
               ELSE
*
*                 Default values
*
                  AFORM   = 'N'
                  ADIAGDO = 'N'
                  OFFDA   = 0
                  CFORM   = 'N'
                  OFFDC   = 0
*
               END IF
*
               CALL PZLAGEN( .FALSE., AFORM, ADIAGDO, OFFDA, MA, NA,
     $                       1, 1, DESCA, IASEED, MEM( IPA ),
     $                       DESCA( LLD_ ) )
*
               IF( BCHECK( L ) )
     $            CALL PZLAGEN( .FALSE., 'None', 'No diag', 0, MB, NB,
     $                          1, 1, DESCB, IBSEED, MEM( IPB ),
     $                          DESCB( LLD_ ) )
*
               IF( CCHECK( L ) )
     $            CALL PZLAGEN( .FALSE., CFORM, 'No diag', OFFDC, MC,
     $                          NC, 1, 1, DESCC, ICSEED, MEM( IPC ),
     $                          DESCC( LLD_ ) )
*
*              Generate entire matrices on each process.
*
               CALL PB_DESCSET2( DESCAR, MA, NA, IMBA, INBA, MBA, NBA,
     $                           -1, -1, ICTXT, MAX( 1, MA ) )
               CALL PZLAGEN( .FALSE., AFORM, ADIAGDO, OFFDA, MA, NA,
     $                       1, 1, DESCAR, IASEED, MEM( IPMATA ),
     $                       DESCAR( LLD_ ) )
*
               IF( BCHECK( L ) ) THEN
                  CALL PB_DESCSET2( DESCBR, MB, NB, IMBB, INBB, MBB,
     $                              NBB, -1, -1, ICTXT, MAX( 1, MB ) )
                  CALL PZLAGEN( .FALSE., 'None', 'No diag', 0, MB, NB,
     $                          1, 1, DESCBR, IBSEED, MEM( IPMATB ),
     $                          DESCBR( LLD_ ) )
               END IF
*
               IF( CCHECK( L ) ) THEN
*
                  CALL PB_DESCSET2( DESCCR, MC, NC, IMBC, INBC, MBC,
     $                              NBC, -1, -1, ICTXT, MAX( 1, MC ) )
                  CALL PZLAGEN( .FALSE., CFORM, 'No diag', OFFDC, MC,
     $                          NC, 1, 1, DESCCR, ICSEED, MEM( IPMATC ),
     $                          DESCCR( LLD_ ) )
*
               ELSE
*
*                 If C is not needed, generate a copy of B instead
*
                  CALL PB_DESCSET2( DESCCR, MB, NB, IMBB, INBB, MBB,
     $                              NBB, -1, -1, ICTXT, MAX( 1, MB ) )
                  CALL PZLAGEN( .FALSE., 'None', 'No diag', 0, MB, NB,
     $                          1, 1, DESCCR, IBSEED, MEM( IPMATC ),
     $                          DESCCR( LLD_ ) )
*
               END IF
*
*              Zero non referenced part of the matrices A, B, C
*
               IF( ( ( L.EQ.2 ).OR. ( L.EQ.3 ) ).AND.
     $             ( MAX( NROWA, NCOLA ).GT.1 ) ) THEN
*
*                 The distributed matrix A is symmetric or Hermitian
*
                  IF( LSAME( UPLO, 'L' ) ) THEN
*
*                    Zeros the strict upper triangular part of A.
*
                     CALL PZLASET( 'Upper', NROWA-1, NCOLA-1, ROGUE,
     $                             ROGUE, MEM( IPA ), IA, JA+1, DESCA )
*
                  ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
*                    Zeros the strict lower triangular part of A.
*
                     CALL PZLASET( 'Lower', NROWA-1, NCOLA-1, ROGUE,
     $                             ROGUE, MEM( IPA ), IA+1, JA, DESCA )
*
                  END IF
*
               ELSE IF( ( ( L.EQ.4 ).OR.( L.EQ.5 ).OR.( L.EQ.6 ).OR.
     $                    ( L.EQ.7 ) ).AND.
     $                  ( MAX( NROWC, NCOLC ).GT.1 ) ) THEN
*
*                 The distributed matrix C is symmetric or Hermitian
*
                  IF( LSAME( UPLO, 'L' ) ) THEN
*
*                    Zeros the strict upper triangular part of C.
*
                     IF( MAX( NROWC, NCOLC ).GT.1 ) THEN
                        CALL PZLASET( 'Upper', NROWC-1, NCOLC-1, ROGUE,
     $                                ROGUE, MEM( IPC ), IC, JC+1,
     $                                DESCC )
                        CALL PB_ZLASET( 'Upper', NROWC-1, NCOLC-1, 0,
     $                                  ROGUE, ROGUE,
     $                                  MEM( IPMATC+IC-1+JC*LDC ), LDC )
                     END IF
*
                  ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
*                    Zeros the strict lower triangular part of C.
*
                     IF( MAX( NROWC, NCOLC ).GT.1 ) THEN
                        CALL PZLASET( 'Lower', NROWC-1, NCOLC-1, ROGUE,
     $                                ROGUE, MEM( IPC ), IC+1, JC,
     $                                DESCC )
                        CALL PB_ZLASET( 'Lower', NROWC-1, NCOLC-1, 0,
     $                                  ROGUE, ROGUE,
     $                                  MEM( IPMATC+IC+(JC-1)*LDC ),
     $                                  LDC )
                     END IF
*
                  END IF
*
               ELSE IF( L.EQ.8 .OR. L.EQ.9 ) THEN
*
                  IF( LSAME( UPLO, 'L' ) ) THEN
*
*                    The distributed matrix A is lower triangular
*
                     IF( LSAME( DIAG, 'N' ) ) THEN
*
                        IF( MAX( NROWA, NCOLA ).GT.1 ) THEN
                           CALL PZLASET( 'Upper', NROWA-1, NCOLA-1,
     $                                   ROGUE, ROGUE, MEM( IPA ), IA,
     $                                   JA+1, DESCA )
                           CALL PB_ZLASET( 'Upper', NROWA-1, NCOLA-1, 0,
     $                                     ZERO, ZERO,
     $                                     MEM( IPMATA+IA-1+JA*LDA ),
     $                                     LDA )
                        END IF
*
                     ELSE
*
                        CALL PZLASET( 'Upper', NROWA, NCOLA, ROGUE, ONE,
     $                                MEM( IPA ), IA, JA, DESCA )
                        CALL PB_ZLASET( 'Upper', NROWA, NCOLA, 0, ZERO,
     $                                  ONE,
     $                                  MEM( IPMATA+IA-1+(JA-1)*LDA ),
     $                                  LDA )
                        IF( ( L.EQ.9 ).AND.
     $                      ( MAX( NROWA, NCOLA ).GT.1 ) ) THEN
                           SCALE = ONE /
     $                             DCMPLX( DBLE( MAX( NROWA, NCOLA ) ) )
                           CALL PZLASCAL( 'Lower', NROWA-1, NCOLA-1,
     $                                    SCALE, MEM( IPA ), IA+1, JA,
     $                                    DESCA )
                           CALL PB_ZLASCAL( 'Lower', NROWA-1, NCOLA-1,
     $                                  0, SCALE,
     $                                  MEM( IPMATA+IA+(JA-1)*LDA ),
     $                                  LDA )
                        END IF
                     END IF
*
                  ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
*                    The distributed matrix A is upper triangular
*
                     IF( LSAME( DIAG, 'N' ) ) THEN
*
                        IF( MAX( NROWA, NCOLA ).GT.1 ) THEN
                           CALL PZLASET( 'Lower', NROWA-1, NCOLA-1,
     $                                   ROGUE, ROGUE, MEM( IPA ), IA+1,
     $                                   JA, DESCA )
                           CALL PB_ZLASET( 'Lower', NROWA-1, NCOLA-1, 0,
     $                                     ZERO, ZERO,
     $                                     MEM( IPMATA+IA+(JA-1)*LDA ),
     $                                     LDA )
                        END IF
*
                     ELSE
*
                        CALL PZLASET( 'Lower', NROWA, NCOLA, ROGUE, ONE,
     $                                MEM( IPA ), IA, JA, DESCA )
                        CALL PB_ZLASET( 'Lower', NROWA, NCOLA, 0, ZERO,
     $                                  ONE,
     $                                  MEM( IPMATA+IA-1+(JA-1)*LDA ),
     $                                  LDA )
                        IF( ( L.EQ.9 ).AND.
     $                      ( MAX( NROWA, NCOLA ).GT.1 ) ) THEN
                           SCALE = ONE /
     $                             DCMPLX( DBLE( MAX( NROWA, NCOLA ) ) )
                           CALL PZLASCAL( 'Upper', NROWA-1, NCOLA-1,
     $                                    SCALE, MEM( IPA ), IA, JA+1,
     $                                    DESCA )
                           CALL PB_ZLASCAL( 'Upper', NROWA-1, NCOLA-1,
     $                                  0, SCALE,
     $                                  MEM( IPMATA+IA-1+JA*LDA ), LDA )
                       END IF
*
                     END IF
*
                  END IF
*
               ELSE IF( L.EQ.11 ) THEN
*
                  IF( LSAME( UPLO, 'L' ) ) THEN
*
*                    The distributed matrix C is lower triangular
*
                     IF( MAX( NROWC, NCOLC ).GT.1 ) THEN
                        CALL PZLASET( 'Upper', NROWC-1, NCOLC-1,
     $                                ROGUE, ROGUE, MEM( IPC ), IC,
     $                                JC+1, DESCC )
                        CALL PB_ZLASET( 'Upper', NROWC-1, NCOLC-1, 0,
     $                                  ROGUE, ROGUE,
     $                                  MEM( IPMATC+IC-1+JC*LDC ), LDC )
                     END IF
*
                  ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
*                    The distributed matrix C is upper triangular
*
                     IF( MAX( NROWC, NCOLC ).GT.1 ) THEN
                        CALL PZLASET( 'Lower', NROWC-1, NCOLC-1,
     $                                ROGUE, ROGUE, MEM( IPC ), IC+1,
     $                                JC, DESCC )
                        CALL PB_ZLASET( 'Lower', NROWC-1, NCOLC-1, 0,
     $                                  ROGUE, ROGUE,
     $                                  MEM( IPMATC+IC+(JC-1)*LDC ),
     $                                  LDC )
                     END IF
*
                  END IF
*
               END IF
*
*              Pad the guard zones of A, B and C
*
               CALL PB_ZFILLPAD( ICTXT, MPA, NQA, MEM( IPA-IPREA ),
     $                           DESCA( LLD_ ), IPREA, IPOSTA, PADVAL )
*
               IF( BCHECK( L ) ) THEN
                  CALL PB_ZFILLPAD( ICTXT, MPB, NQB, MEM( IPB-IPREB ),
     $                              DESCB( LLD_ ), IPREB, IPOSTB,
     $                              PADVAL )
               END IF
*
               IF( CCHECK( L ) ) THEN
                  CALL PB_ZFILLPAD( ICTXT, MPC, NQC, MEM( IPC-IPREC ),
     $                              DESCC( LLD_ ), IPREC, IPOSTC,
     $                              PADVAL )
               END IF
*
*              Initialize the check for INPUT-only arguments.
*
               INFO = 0
               CALL PZCHKARG3( ICTXT, NOUT, SNAMES( L ), SIDE, UPLO,
     $                         TRANSA, TRANSB, DIAG, M, N, K, ALPHA, IA,
     $                         JA, DESCA, IB, JB, DESCB, BETA, IC, JC,
     $                         DESCC, INFO )
*
*              Print initial parallel data if IVERB >= 2.
*
               IF( IVERB.EQ.2 ) THEN
                  CALL PB_PZLAPRNT( NROWA, NCOLA, MEM( IPA ), IA, JA,
     $                              DESCA, 0, 0,
     $                              'PARALLEL_INITIAL_A', NOUT,
     $                              MEM( IPW ) )
               ELSE IF( IVERB.GE.3 ) THEN
                  CALL PB_PZLAPRNT( MA, NA, MEM( IPA ), 1, 1, DESCA,
     $                              0, 0, 'PARALLEL_INITIAL_A', NOUT,
     $                              MEM( IPW ) )
               END IF
*
               IF( BCHECK( L ) ) THEN
                  IF( IVERB.EQ.2 ) THEN
                     CALL PB_PZLAPRNT( NROWB, NCOLB, MEM( IPB ), IB, JB,
     $                                 DESCB, 0, 0,
     $                                 'PARALLEL_INITIAL_B', NOUT,
     $                                 MEM( IPW ) )
                  ELSE IF( IVERB.GE.3 ) THEN
                     CALL PB_PZLAPRNT( MB, NB, MEM( IPB ), 1, 1, DESCB,
     $                                 0, 0, 'PARALLEL_INITIAL_B', NOUT,
     $                                 MEM( IPW ) )
                  END IF
               END IF
*
               IF( CCHECK( L ) ) THEN
                  IF( IVERB.EQ.2 ) THEN
                     CALL PB_PZLAPRNT( NROWC, NCOLC, MEM( IPC ), IC, JC,
     $                                 DESCC, 0, 0,
     $                                 'PARALLEL_INITIAL_C', NOUT,
     $                                 MEM( IPW ) )
                  ELSE IF( IVERB.GE.3 ) THEN
                     CALL PB_PZLAPRNT( MC, NC, MEM( IPC ), 1, 1, DESCC,
     $                                 0, 0, 'PARALLEL_INITIAL_C', NOUT,
     $                                 MEM( IPW ) )
                  END IF
               END IF
*
*              Call the Level 3 PBLAS routine
*
               INFO = 0
               IF( L.EQ.1 ) THEN
*
*                 Test PZGEMM
*
                  CALL PZGEMM( TRANSA, TRANSB, M, N, K, ALPHA,
     $                         MEM( IPA ), IA, JA, DESCA, MEM( IPB ),
     $                         IB, JB, DESCB, BETA, MEM( IPC ), IC, JC,
     $                         DESCC )
*
               ELSE IF( L.EQ.2 ) THEN
*
*                 Test PZSYMM
*
                  CALL PZSYMM( SIDE, UPLO, M, N, ALPHA, MEM( IPA ), IA,
     $                         JA, DESCA, MEM( IPB ), IB, JB, DESCB,
     $                         BETA, MEM( IPC ), IC, JC, DESCC )
*
               ELSE IF( L.EQ.3 ) THEN
*
*                 Test PZHEMM
*
                  CALL PZIPSET( 'Bignum', NROWA, MEM( IPA ), IA, JA,
     $                          DESCA )
*
                  CALL PZHEMM( SIDE, UPLO, M, N, ALPHA, MEM( IPA ), IA,
     $                         JA, DESCA, MEM( IPB ), IB, JB, DESCB,
     $                         BETA, MEM( IPC ), IC, JC, DESCC )
*
                  CALL PZIPSET( 'Zero', NROWA, MEM( IPA ), IA, JA,
     $                          DESCA )
*
               ELSE IF( L.EQ.4 ) THEN
*
*                 Test PZSYRK
*
                  CALL PZSYRK( UPLO, TRANSA, N, K, ALPHA, MEM( IPA ),
     $                         IA, JA, DESCA, BETA, MEM( IPC ), IC, JC,
     $                         DESCC )
*
               ELSE IF( L.EQ.5 ) THEN
*
*                 Test PZHERK
*
                  IF( ( ( DCMPLX( DBLE( ALPHA ) ).NE.ZERO ).AND.
     $                  ( K.NE.0 ) ).OR.
     $                ( DCMPLX( DBLE( BETA  ) ).NE.ONE  ) )
     $               CALL PZIPSET( 'Bignum', N, MEM( IPC ), IC, JC,
     $                             DESCC )
*
                  CALL PZHERK( UPLO, TRANSA, N, K, DBLE( ALPHA ),
     $                         MEM( IPA ), IA, JA, DESCA, DBLE( BETA ),
     $                         MEM( IPC ), IC, JC, DESCC )
*
               ELSE IF( L.EQ.6 ) THEN
*
*                 Test PZSYR2K
*
                  CALL PZSYR2K( UPLO, TRANSA, N, K, ALPHA, MEM( IPA ),
     $                          IA, JA, DESCA, MEM( IPB ), IB, JB,
     $                          DESCB, BETA, MEM( IPC ), IC, JC,
     $                          DESCC )
*
               ELSE IF( L.EQ.7 ) THEN
*
*                 Test PZHER2K
*
                  IF( ( ( ALPHA.NE.ZERO ).AND.( K.NE.0 ) ).OR.
     $                ( DCMPLX( DBLE( BETA ) ).NE.ONE  ) )
     $               CALL PZIPSET( 'Bignum', N, MEM( IPC ), IC, JC,
     $                             DESCC )
*
                  CALL PZHER2K( UPLO, TRANSA, N, K, ALPHA, MEM( IPA ),
     $                          IA, JA, DESCA, MEM( IPB ), IB, JB,
     $                          DESCB, DBLE( BETA ), MEM( IPC ), IC, JC,
     $                          DESCC )
*
               ELSE IF( L.EQ.8 ) THEN
*
*                 Test PZTRMM
*
                  CALL PZTRMM( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,
     $                         MEM( IPA ), IA, JA, DESCA, MEM( IPB ),
     $                         IB, JB, DESCB )
*
               ELSE IF( L.EQ.9 ) THEN
*
*                 Test PZTRSM
*
                  CALL PZTRSM( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,
     $                         MEM( IPA ), IA, JA, DESCA, MEM( IPB ),
     $                         IB, JB, DESCB )
*
*
               ELSE IF( L.EQ.10 ) THEN
*
*                 Test PZGEADD
*
                  CALL PZGEADD( TRANSA, M, N, ALPHA, MEM( IPA ), IA, JA,
     $                          DESCA, BETA, MEM( IPC ), IC, JC, DESCC )
*
               ELSE IF( L.EQ.11 ) THEN
*
*                 Test PZTRADD
*
                  CALL PZTRADD( UPLO, TRANSA, M, N, ALPHA, MEM( IPA ),
     $                          IA, JA, DESCA, BETA, MEM( IPC ), IC, JC,
     $                          DESCC )
*
               END IF
*
*              Check if the operation has been performed.
*
               IF( INFO.NE.0 ) THEN
                  KSKIP( L ) = KSKIP( L ) + 1
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9974 ) INFO
                  GO TO 30
               END IF
*
*              Check padding
*
               CALL PB_ZCHEKPAD( ICTXT, SNAMES( L ), MPA, NQA,
     $                           MEM( IPA-IPREA ), DESCA( LLD_ ),
     $                           IPREA, IPOSTA, PADVAL )
*
               IF( BCHECK( L ) ) THEN
                  CALL PB_ZCHEKPAD( ICTXT, SNAMES( L ), MPB, NQB,
     $                              MEM( IPB-IPREB ), DESCB( LLD_ ),
     $                              IPREB, IPOSTB, PADVAL )
               END IF
*
               IF( CCHECK( L ) ) THEN
                  CALL PB_ZCHEKPAD( ICTXT, SNAMES( L ), MPC, NQC,
     $                              MEM( IPC-IPREC ), DESCC( LLD_ ),
     $                              IPREC, IPOSTC, PADVAL )
               END IF
*
*              Check the computations
*
               CALL PZBLAS3TSTCHK( ICTXT, NOUT, L, SIDE, UPLO, TRANSA,
     $                             TRANSB, DIAG, M, N, K, ALPHA,
     $                             MEM( IPMATA ), MEM( IPA ), IA, JA,
     $                             DESCA, MEM( IPMATB ), MEM( IPB ),
     $                             IB, JB, DESCB, BETA, MEM( IPMATC ),
     $                             MEM( IPC ), IC, JC, DESCC, THRESH,
     $                             ROGUE, MEM( IPG ), MEM( IPW ), INFO )
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
               CALL PZCHKARG3( ICTXT, NOUT, SNAMES( L ), SIDE, UPLO,
     $                         TRANSA, TRANSB, DIAG, M, N, K, ALPHA, IA,
     $                         JA, DESCA, IB, JB, DESCB, BETA, IC, JC,
     $                         DESCC, INFO )
*
*              Check input-only array arguments
*
               CALL PZCHKMOUT( NROWA, NCOLA, MEM( IPMATA ),
     $                         MEM( IPA ), IA, JA, DESCA, IERR( 4 ) )
               IF( IERR( 4 ).NE.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9983 ) 'PARALLEL_A',
     $                                            SNAMES( L )
               END IF
*
               IF( BCHECK( L ) ) THEN
                  CALL PZCHKMOUT( NROWB, NCOLB, MEM( IPMATB ),
     $                            MEM( IPB ), IB, JB, DESCB, IERR( 5 ) )
                  IF( IERR( 5 ).NE.0 ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9983 ) 'PARALLEL_B',
     $                                            SNAMES( L )
                  END IF
               END IF
*
               IF( CCHECK( L ) ) THEN
                  CALL PZCHKMOUT( NROWC, NCOLC, MEM( IPMATC ),
     $                            MEM( IPC ), IC, JC, DESCC, IERR( 6 ) )
                  IF( IERR( 6 ).NE.0 ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9983 ) 'PARALLEL_C',
     $                                            SNAMES( L )
                  END IF
               END IF
*
*              Only node 0 prints computational test result
*
               IF( INFO.NE.0 .OR. IERR( 1 ).NE.0 .OR.
     $             IERR( 2 ).NE.0 .OR. IERR( 3 ).NE.0 .OR.
     $             IERR( 4 ).NE.0 .OR. IERR( 5 ).NE.0 .OR.
     $             IERR( 6 ).NE.0 ) THEN
                  KFAIL( L ) = KFAIL( L ) + 1
                  ERRFLG = .TRUE.
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9985 ) SNAMES( L )
               ELSE
                  KPASS( L ) = KPASS( L ) + 1
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9984 ) SNAMES( L )
               END IF
*
*              Dump matrix if IVERB >= 1 and error.
*
               IF( IVERB.GE.1 .AND. ERRFLG ) THEN
                  IF( IERR( 4 ).NE.0 .OR. IVERB.GE.3 ) THEN
                     CALL PZMPRNT( ICTXT, NOUT, MA, NA, MEM( IPMATA ),
     $                             LDA, 0, 0, 'SERIAL_A' )
                     CALL PB_PZLAPRNT( MA, NA, MEM( IPA ), 1, 1, DESCA,
     $                                 0, 0, 'PARALLEL_A', NOUT,
     $                                 MEM( IPMATA ) )
                  ELSE IF( IERR( 1 ).NE.0 ) THEN
                     IF( ( NROWA.GT.0 ).AND.( NCOLA.GT.0 ) )
     $                  CALL PZMPRNT( ICTXT, NOUT, NROWA, NCOLA,
     $                                MEM( IPMATA+IA-1+(JA-1)*LDA ),
     $                                LDA, 0, 0, 'SERIAL_A' )
                     CALL PB_PZLAPRNT( NROWA, NCOLA, MEM( IPA ), IA, JA,
     $                                 DESCA, 0, 0, 'PARALLEL_A', NOUT,
     $                                 MEM( IPMATA ) )
                  END IF
                  IF( BCHECK( L ) ) THEN
                     IF( IERR( 5 ).NE.0 .OR. IVERB.GE.3 ) THEN
                        CALL PZMPRNT( ICTXT, NOUT, MB, NB,
     $                                MEM( IPMATB ), LDB, 0, 0,
     $                                'SERIAL_B' )
                        CALL PB_PZLAPRNT( MB, NB, MEM( IPB ), 1, 1,
     $                                    DESCB, 0, 0, 'PARALLEL_B',
     $                                    NOUT, MEM( IPMATB ) )
                     ELSE IF( IERR( 2 ).NE.0 ) THEN
                        IF( ( NROWB.GT.0 ).AND.( NCOLB.GT.0 ) )
     $                     CALL PZMPRNT( ICTXT, NOUT, NROWB, NCOLB,
     $                                   MEM( IPMATB+IB-1+(JB-1)*LDB ),
     $                                   LDB, 0, 0, 'SERIAL_B' )
                        CALL PB_PZLAPRNT( NROWB, NCOLB, MEM( IPB ), IB,
     $                                    JB, DESCB, 0, 0, 'PARALLEL_B',
     $                                    NOUT, MEM( IPMATB ) )
                     END IF
                  END IF
                  IF( CCHECK( L ) ) THEN
                     IF( IERR( 6 ).NE.0 .OR. IVERB.GE.3 ) THEN
                        CALL PZMPRNT( ICTXT, NOUT, MC, NC,
     $                                MEM( IPMATC ), LDC, 0, 0,
     $                                'SERIAL_C' )
                        CALL PB_PZLAPRNT( MC, NC, MEM( IPC ), 1, 1,
     $                                    DESCC, 0, 0, 'PARALLEL_C',
     $                                    NOUT, MEM( IPMATC ) )
                     ELSE IF( IERR( 3 ).NE.0 ) THEN
                        IF( ( NROWB.GT.0 ).AND.( NCOLB.GT.0 ) )
     $                     CALL PZMPRNT( ICTXT, NOUT, NROWC, NCOLC,
     $                                   MEM( IPMATC+IC-1+(JC-1)*LDC ),
     $                                   LDC, 0, 0, 'SERIAL_C' )
                        CALL PB_PZLAPRNT( NROWC, NCOLC, MEM( IPC ), IC,
     $                                    JC, DESCC, 0, 0, 'PARALLEL_C',
     $                                    NOUT, MEM( IPMATC ) )
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
               WRITE( NOUT, FMT = 9982 ) J
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
         WRITE( NOUT, FMT = 9978 )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9980 )
         WRITE( NOUT, FMT = 9979 )
*
         DO 90 I = 1, NSUBS
            WRITE( NOUT, FMT = 9981 ) '|', SNAMES( I ), KTESTS( I ),
     $                                KPASS( I ), KFAIL( I ), KSKIP( I )
   90    CONTINUE
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9977 )
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
     $        '-------------------' )
 9994 FORMAT( 2X, '        M      N      K    SIDE  UPLO  TRANSA  ',
     $        'TRANSB  DIAG' )
 9993 FORMAT( 5X,I6,1X,I6,1X,I6,6X,A1,5X,A1,7X,A1,7X,A1,5X,A1 )
 9992 FORMAT( 2X, '       IA     JA     MA     NA   IMBA   INBA',
     $        '    MBA    NBA RSRCA CSRCA' )
 9991 FORMAT( 5X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,I6,
     $        1X,I5,1X,I5 )
 9990 FORMAT( 2X, '       IB     JB     MB     NB   IMBB   INBB',
     $        '    MBB    NBB RSRCB CSRCB' )
 9989 FORMAT( 2X, '       IC     JC     MC     NC   IMBC   INBC',
     $        '    MBC    NBC RSRCC CSRCC' )
 9988 FORMAT( 'Not enough memory for this test: going on to',
     $        ' next test case.' )
 9987 FORMAT( 'Not enough memory. Need: ', I12 )
 9986 FORMAT( 2X, '   Tested Subroutine: ', A )
 9985 FORMAT( 2X, '   ***** Computational check: ', A, '       ',
     $        ' FAILED ',' *****' )
 9984 FORMAT( 2X, '   ***** Computational check: ', A, '       ',
     $        ' PASSED ',' *****' )
 9983 FORMAT( 2X, '   ***** ERROR ***** Matrix operand ', A,
     $        ' modified by ', A, ' *****' )
 9982 FORMAT( 2X, 'Test number ', I4, ' completed.' )
 9981 FORMAT( 2X,A1,2X,A7,8X,I4,6X,I4,5X,I4,4X,I4 )
 9980 FORMAT( 2X, '   SUBROUTINE  TOTAL TESTS  PASSED   FAILED  ',
     $        'SKIPPED' )
 9979 FORMAT( 2X, '   ----------  -----------  ------   ------  ',
     $        '-------' )
 9978 FORMAT( 2X, 'Testing Summary')
 9977 FORMAT( 2X, 'End of Tests.' )
 9976 FORMAT( 2X, 'Tests started.' )
 9975 FORMAT( 2X, '   ***** ', A, ' has an incorrect value:     ',
     $            ' BYPASS  *****' )
 9974 FORMAT( 2X, '   ***** Operation not supported, error code: ',
     $        I5, ' *****' )
*
      STOP
*
*     End of PZBLA3TST
*
      END
      SUBROUTINE PZBLA3TSTINFO( SUMMRY, NOUT, NMAT, DIAGVAL, SIDEVAL,
     $                          TRNAVAL, TRNBVAL, UPLOVAL, MVAL,
     $                          NVAL, KVAL, MAVAL, NAVAL, IMBAVAL,
     $                          MBAVAL, INBAVAL, NBAVAL, RSCAVAL,
     $                          CSCAVAL, IAVAL, JAVAL, MBVAL, NBVAL,
     $                          IMBBVAL, MBBVAL, INBBVAL, NBBVAL,
     $                          RSCBVAL, CSCBVAL, IBVAL, JBVAL,
     $                          MCVAL, NCVAL, IMBCVAL, MBCVAL,
     $                          INBCVAL, NBCVAL, RSCCVAL, CSCCVAL,
     $                          ICVAL, JCVAL, LDVAL, NGRIDS, PVAL,
     $                          LDPVAL, QVAL, LDQVAL, NBLOG, LTEST, SOF,
     $                          TEE, IAM, IGAP, IVERB, NPROCS, THRESH,
     $                          ALPHA, BETA, WORK )
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
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      CHARACTER*( * )    SUMMRY
      CHARACTER*1        DIAGVAL( LDVAL ), SIDEVAL( LDVAL ),
     $                   TRNAVAL( LDVAL ), TRNBVAL( LDVAL ),
     $                   UPLOVAL( LDVAL )
      LOGICAL            LTEST( * )
      INTEGER            CSCAVAL( LDVAL ), CSCBVAL( LDVAL ),
     $                   CSCCVAL( LDVAL ), IAVAL( LDVAL ),
     $                   IBVAL( LDVAL ), ICVAL( LDVAL ),
     $                   IMBAVAL( LDVAL ), IMBBVAL( LDVAL ),
     $                   IMBCVAL( LDVAL ), INBAVAL( LDVAL ),
     $                   INBBVAL( LDVAL ), INBCVAL( LDVAL ),
     $                   JAVAL( LDVAL ), JBVAL( LDVAL ), JCVAL( LDVAL ),
     $                   KVAL( LDVAL ), MAVAL( LDVAL ), MBAVAL( LDVAL ),
     $                   MBBVAL( LDVAL ), MBCVAL( LDVAL ),
     $                   MBVAL( LDVAL ), MCVAL( LDVAL ), MVAL( LDVAL ),
     $                   NAVAL( LDVAL ), NBAVAL( LDVAL ),
     $                   NBBVAL( LDVAL ), NBCVAL( LDVAL ),
     $                   NBVAL( LDVAL ), NCVAL( LDVAL ), NVAL( LDVAL ),
     $                   PVAL( LDPVAL ), QVAL( LDQVAL ),
     $                   RSCAVAL( LDVAL ), RSCBVAL( LDVAL ),
     $                   RSCCVAL( LDVAL ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZBLA3TSTINFO  get the needed startup information for testing various
*  Level 3 PBLAS routines, and transmits it to all processes.
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
*  SIDEVAL (global output) CHARACTER array
*          On entry,  SIDEVAL  is  an array of dimension LDVAL. On exit,
*          this array contains the values of SIDE to run the code with.
*
*  TRNAVAL (global output) CHARACTER array
*          On entry, TRNAVAL  is an array of dimension LDVAL.  On  exit,
*          this array contains  the  values  of  TRANSA  to run the code
*          with.
*
*  TRNBVAL (global output) CHARACTER array
*          On entry, TRNBVAL  is an array of dimension LDVAL.  On  exit,
*          this array contains  the  values  of  TRANSB  to run the code
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
*  KVAL    (global output) INTEGER array
*          On entry, KVAL is an array of dimension LDVAL.  On exit, this
*          array contains the values of K to run the code with.
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
*  MBVAL   (global output) INTEGER array
*          On entry, MBVAL is an array of dimension LDVAL. On exit, this
*          array  contains  the values  of  DESCB( M_ )  to run the code
*          with.
*
*  NBVAL   (global output) INTEGER array
*          On entry, NBVAL is an array of dimension LDVAL. On exit, this
*          array  contains  the values  of  DESCB( N_ )  to run the code
*          with.
*
*  IMBBVAL (global output) INTEGER array
*          On entry,  IMBBVAL  is an array of  dimension LDVAL. On exit,
*          this  array  contains  the values of DESCB( IMB_ ) to run the
*          code with.
*
*  MBBVAL  (global output) INTEGER array
*          On entry,  MBBVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains  the values of DESCB( MB_ ) to  run the
*          code with.
*
*  INBBVAL (global output) INTEGER array
*          On entry,  INBBVAL  is an array of  dimension LDVAL. On exit,
*          this  array  contains  the values of DESCB( INB_ ) to run the
*          code with.
*
*  NBBVAL  (global output) INTEGER array
*          On entry,  NBBVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains  the values of DESCB( NB_ ) to  run the
*          code with.
*
*  RSCBVAL (global output) INTEGER array
*          On entry, RSCBVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains the values of DESCB( RSRC_ ) to run the
*          code with.
*
*  CSCBVAL (global output) INTEGER array
*          On entry, CSCBVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains the values of DESCB( CSRC_ ) to run the
*          code with.
*
*  IBVAL   (global output) INTEGER array
*          On entry, IBVAL is an array of dimension LDVAL. On exit, this
*          array  contains the values of IB to run the code with.
*
*  JBVAL   (global output) INTEGER array
*          On entry, JBVAL is an array of dimension LDVAL. On exit, this
*          array  contains the values of JB to run the code with.
*
*  MCVAL   (global output) INTEGER array
*          On entry, MCVAL is an array of dimension LDVAL. On exit, this
*          array  contains  the values  of  DESCC( M_ )  to run the code
*          with.
*
*  NCVAL   (global output) INTEGER array
*          On entry, NCVAL is an array of dimension LDVAL. On exit, this
*          array  contains  the values  of  DESCC( N_ )  to run the code
*          with.
*
*  IMBCVAL (global output) INTEGER array
*          On entry,  IMBCVAL  is an array of  dimension LDVAL. On exit,
*          this  array  contains  the values of DESCC( IMB_ ) to run the
*          code with.
*
*  MBCVAL  (global output) INTEGER array
*          On entry,  MBCVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains  the values of DESCC( MB_ ) to  run the
*          code with.
*
*  INBCVAL (global output) INTEGER array
*          On entry,  INBCVAL  is an array of  dimension LDVAL. On exit,
*          this  array  contains  the values of DESCC( INB_ ) to run the
*          code with.
*
*  NBCVAL  (global output) INTEGER array
*          On entry,  NBCVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains  the values of DESCC( NB_ ) to  run the
*          code with.
*
*  RSCCVAL (global output) INTEGER array
*          On entry, RSCCVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains the values of DESCC( RSRC_ ) to run the
*          code with.
*
*  CSCCVAL (global output) INTEGER array
*          On entry, CSCCVAL  is an array of  dimension  LDVAL. On exit,
*          this  array  contains the values of DESCC( CSRC_ ) to run the
*          code with.
*
*  ICVAL   (global output) INTEGER array
*          On entry, ICVAL is an array of dimension LDVAL. On exit, this
*          array  contains the values of IC to run the code with.
*
*  JCVAL   (global output) INTEGER array
*          On entry, JCVAL is an array of dimension LDVAL. On exit, this
*          array  contains the values of JC to run the code with.
*
*  LDVAL   (global input) INTEGER
*          On entry, LDVAL specifies the maximum number of different va-
*          lues  that  can be used for DIAG, SIDE, TRANSA, TRANSB, UPLO,
*          M,  N,  K,  DESCA(:), IA, JA, DESCB(:), IB, JB, DESCC(:), IC,
*          JC. This is also the maximum number of test cases.
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
*          On entry, LTEST  is an array of dimension at least eleven. On
*          exit, if LTEST( i ) is .TRUE., the i-th Level 3 PBLAS routine
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
*  ALPHA   (global output) COMPLEX*16
*          On exit, ALPHA specifies the value of alpha to be used in all
*          the test cases.
*
*  BETA    (global output) COMPLEX*16
*          On exit, BETA  specifies the value of beta  to be used in all
*          the test cases.
*
*  WORK    (local workspace) INTEGER array
*          On   entry,   WORK   is   an  array  of  dimension  at  least
*          MAX( 3, 2*NGRIDS+38*NMAT+NSUBS+4 )  with  NSUBS  equal to 11.
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
      PARAMETER          ( NIN = 11, NSUBS = 11 )
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
     $                   BLACS_GRIDINIT, BLACS_SETUP, ICOPY, IGEBR2D,
     $                   IGEBS2D, SGEBR2D, SGEBS2D, ZGEBR2D, ZGEBS2D
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
         OPEN( NIN, FILE='PZBLAS3TST.dat', STATUS='OLD' )
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
         READ( NIN, FMT = * ) ( DIAGVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( SIDEVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( TRNAVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( TRNBVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( UPLOVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( MVAL   ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( NVAL   ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( KVAL   ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( MAVAL  ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( NAVAL  ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( IMBAVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( INBAVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( MBAVAL ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( NBAVAL ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( RSCAVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( CSCAVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( IAVAL  ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( JAVAL  ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( MBVAL  ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( NBVAL  ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( IMBBVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( INBBVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( MBBVAL ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( NBBVAL ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( RSCBVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( CSCBVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( IBVAL  ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( JBVAL  ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( MCVAL  ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( NCVAL  ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( IMBCVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( INBCVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( MBCVAL ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( NBCVAL ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( RSCCVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( CSCCVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( ICVAL  ( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( JCVAL  ( I ), I = 1, NMAT )
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
         CALL ZGEBS2D( ICTXT, 'All', ' ', 1, 1, ALPHA, 1 )
         CALL ZGEBS2D( ICTXT, 'All', ' ', 1, 1, BETA,  1 )
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
            WORK( I   ) = ICHAR( DIAGVAL( J ) )
            WORK( I+1 ) = ICHAR( SIDEVAL( J ) )
            WORK( I+2 ) = ICHAR( TRNAVAL( J ) )
            WORK( I+3 ) = ICHAR( TRNBVAL( J ) )
            WORK( I+4 ) = ICHAR( UPLOVAL( J ) )
            I = I + 5
   70    CONTINUE
         CALL ICOPY( NGRIDS, PVAL,     1, WORK( I ), 1 )
         I = I + NGRIDS
         CALL ICOPY( NGRIDS, QVAL,     1, WORK( I ), 1 )
         I = I + NGRIDS
         CALL ICOPY( NMAT,   MVAL,     1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   NVAL,     1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   KVAL,     1, WORK( I ), 1 )
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
         CALL ICOPY( NMAT,   MBVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   NBVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   IMBBVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   INBBVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   MBBVAL,   1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   NBBVAL,   1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   RSCBVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   CSCBVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   IBVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   JBVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   MCVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   NCVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   IMBCVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   INBCVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   MBCVAL,   1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   NBCVAL,   1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   RSCCVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   CSCCVAL,  1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   ICVAL,    1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   JCVAL,    1, WORK( I ), 1 )
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
         WRITE( NOUT, FMT = 9999 ) 'Level 3 PBLAS testing program.'
         WRITE( NOUT, FMT = 9999 ) USRINFO
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'Tests of the complex double precision '//
     $               'Level 3 PBLAS'
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
         CALL ZGEBR2D( ICTXT, 'All', ' ', 1, 1, ALPHA, 1, 0, 0 )
         CALL ZGEBR2D( ICTXT, 'All', ' ', 1, 1, BETA,  1, 0, 0 )
*
         CALL IGEBR2D( ICTXT, 'All', ' ', 3, 1, WORK, 3, 0, 0 )
         NGRIDS = WORK( 1 )
         NMAT   = WORK( 2 )
         NBLOG  = WORK( 3 )
*
         I = 2*NGRIDS + 38*NMAT + NSUBS + 4
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
            DIAGVAL( J ) = CHAR( WORK( I   ) )
            SIDEVAL( J ) = CHAR( WORK( I+1 ) )
            TRNAVAL( J ) = CHAR( WORK( I+2 ) )
            TRNBVAL( J ) = CHAR( WORK( I+3 ) )
            UPLOVAL( J ) = CHAR( WORK( I+4 ) )
            I = I + 5
  100    CONTINUE
         CALL ICOPY( NGRIDS, WORK( I ), 1, PVAL,     1 )
         I = I + NGRIDS
         CALL ICOPY( NGRIDS, WORK( I ), 1, QVAL,     1 )
         I = I + NGRIDS
         CALL ICOPY( NMAT,   WORK( I ), 1, MVAL,     1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, NVAL,     1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, KVAL,     1 )
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
         CALL ICOPY( NMAT,   WORK( I ), 1, MBVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, NBVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, IMBBVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, INBBVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, MBBVAL,   1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, NBBVAL,   1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, RSCBVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, CSCBVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, IBVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, JBVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, MCVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, NCVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, IMBCVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, INBCVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, MBCVAL,   1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, NBCVAL,   1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, RSCCVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, CSCCVAL,  1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, ICVAL,    1 )
         I = I + NMAT
         CALL ICOPY( NMAT,   WORK( I ), 1, JCVAL,    1 )
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
 9982 FORMAT( 2X, 'Alpha                     :      (', G16.6,
     $        ',', G16.6, ')' )
 9981 FORMAT( 2X, 'Beta                      :      (', G16.6,
     $        ',', G16.6, ')' )
 9980 FORMAT( 2X, 'Threshold value           : ', G16.6 )
 9979 FORMAT( 2X, 'Logical block size        : ', I6 )
*
*     End of PZBLA3TSTINFO
*
      END
      SUBROUTINE PZBLAS3TSTCHKE( LTEST, INOUT, NPROCS )
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
*  PZBLAS3TSTCHKE tests the error exits of the Level 3 PBLAS.
*
*  Arguments
*  =========
*
*  LTEST   (global input) LOGICAL array
*          On entry, LTEST is an array of dimension at least 11 (NSUBS).
*             If LTEST(  1 ) is .TRUE., PZGEMM  will be tested;
*             If LTEST(  2 ) is .TRUE., PZSYMM  will be tested;
*             If LTEST(  3 ) is .TRUE., PZHEMM  will be tested;
*             If LTEST(  4 ) is .TRUE., PZSYRK  will be tested;
*             If LTEST(  5 ) is .TRUE., PZHERK  will be tested;
*             If LTEST(  6 ) is .TRUE., PZSYR2K will be tested;
*             If LTEST(  7 ) is .TRUE., PZHER2K will be tested;
*             If LTEST(  8 ) is .TRUE., PZTRMM  will be tested;
*             If LTEST(  9 ) is .TRUE., PZTRSM  will be tested;
*             If LTEST( 10 ) is .TRUE., PZGEADD will be tested;
*             If LTEST( 11 ) is .TRUE., PZTRADD will be tested;
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
      PARAMETER          ( NSUBS = 11 )
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
     $                   BLACS_GRIDINIT, PZDIMEE, PZGEADD, PZGEMM,
     $                   PZHEMM, PZHER2K, PZHERK, PZMATEE, PZOPTEE,
     $                   PZSYMM, PZSYR2K, PZSYRK, PZTRADD, PZTRMM,
     $                   PZTRSM
*     ..
*     .. Common Blocks ..
      LOGICAL            ABRTFLG
      INTEGER            NOUT
      CHARACTER*7        SNAMES( NSUBS )
      COMMON             /SNAMEC/SNAMES
      COMMON             /PBERRORC/NOUT, ABRTFLG
*     ..
*     .. Data Statements ..
      DATA               SCODE/31, 32, 32, 33, 34, 35, 36, 38, 38, 39,
     $                   40/
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
*     Test PZGEMM
*
      I = 1
      IF( LTEST( I ) ) THEN
         CALL PZOPTEE( ICTXT, NOUT, PZGEMM, SCODE( I ), SNAMES( I ) )
         CALL PZDIMEE( ICTXT, NOUT, PZGEMM, SCODE( I ), SNAMES( I ) )
         CALL PZMATEE( ICTXT, NOUT, PZGEMM, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PZSYMM
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PZOPTEE( ICTXT, NOUT, PZSYMM, SCODE( I ), SNAMES( I ) )
         CALL PZDIMEE( ICTXT, NOUT, PZSYMM, SCODE( I ), SNAMES( I ) )
         CALL PZMATEE( ICTXT, NOUT, PZSYMM, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PZHEMM
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PZOPTEE( ICTXT, NOUT, PZHEMM, SCODE( I ), SNAMES( I ) )
         CALL PZDIMEE( ICTXT, NOUT, PZHEMM, SCODE( I ), SNAMES( I ) )
         CALL PZMATEE( ICTXT, NOUT, PZHEMM, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PZSYRK
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PZOPTEE( ICTXT, NOUT, PZSYRK, SCODE( I ), SNAMES( I ) )
         CALL PZDIMEE( ICTXT, NOUT, PZSYRK, SCODE( I ), SNAMES( I ) )
         CALL PZMATEE( ICTXT, NOUT, PZSYRK, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PZHERK
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PZOPTEE( ICTXT, NOUT, PZHERK, SCODE( I ), SNAMES( I ) )
         CALL PZDIMEE( ICTXT, NOUT, PZHERK, SCODE( I ), SNAMES( I ) )
         CALL PZMATEE( ICTXT, NOUT, PZHERK, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PZSYR2K
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PZOPTEE( ICTXT, NOUT, PZSYR2K, SCODE( I ), SNAMES( I ) )
         CALL PZDIMEE( ICTXT, NOUT, PZSYR2K, SCODE( I ), SNAMES( I ) )
         CALL PZMATEE( ICTXT, NOUT, PZSYR2K, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PZHER2K
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PZOPTEE( ICTXT, NOUT, PZHER2K, SCODE( I ), SNAMES( I ) )
         CALL PZDIMEE( ICTXT, NOUT, PZHER2K, SCODE( I ), SNAMES( I ) )
         CALL PZMATEE( ICTXT, NOUT, PZHER2K, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PZTRMM
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PZOPTEE( ICTXT, NOUT, PZTRMM, SCODE( I ), SNAMES( I ) )
         CALL PZDIMEE( ICTXT, NOUT, PZTRMM, SCODE( I ), SNAMES( I ) )
         CALL PZMATEE( ICTXT, NOUT, PZTRMM, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PZTRSM
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PZOPTEE( ICTXT, NOUT, PZTRSM, SCODE( I ), SNAMES( I ) )
         CALL PZDIMEE( ICTXT, NOUT, PZTRSM, SCODE( I ), SNAMES( I ) )
         CALL PZMATEE( ICTXT, NOUT, PZTRSM, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PZGEADD
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PZOPTEE( ICTXT, NOUT, PZGEADD, SCODE( I ), SNAMES( I ) )
         CALL PZDIMEE( ICTXT, NOUT, PZGEADD, SCODE( I ), SNAMES( I ) )
         CALL PZMATEE( ICTXT, NOUT, PZGEADD, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PZTRADD
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PZOPTEE( ICTXT, NOUT, PZTRADD, SCODE( I ), SNAMES( I ) )
         CALL PZDIMEE( ICTXT, NOUT, PZTRADD, SCODE( I ), SNAMES( I ) )
         CALL PZMATEE( ICTXT, NOUT, PZTRADD, SCODE( I ), SNAMES( I ) )
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
*     End of PZBLAS3TSTCHKE
*
      END
      SUBROUTINE PZCHKARG3( ICTXT, NOUT, SNAME, SIDE, UPLO, TRANSA,
     $                      TRANSB, DIAG, M, N, K, ALPHA, IA, JA,
     $                      DESCA, IB, JB, DESCB, BETA, IC, JC, DESCC,
     $                      INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        DIAG, SIDE, TRANSA, TRANSB, UPLO
      INTEGER            IA, IB, IC, ICTXT, INFO, JA, JB, JC, K, M, N,
     $                   NOUT
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      CHARACTER*7        SNAME
      INTEGER            DESCA( * ), DESCB( * ), DESCC( * )
*     ..
*
*  Purpose
*  =======
*
*  PZCHKARG3 checks the input-only arguments of the Level 3 PBLAS.  When
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
*  SIDE    (global input) CHARACTER*1
*          On entry, SIDE specifies the SIDE option in the Level 3 PBLAS
*          operation.
*
*  UPLO    (global input) CHARACTER*1
*          On entry, UPLO specifies the UPLO option in the Level 3 PBLAS
*          operation.
*
*  TRANSA  (global input) CHARACTER*1
*          On entry,  TRANSA  specifies the TRANSA option in the Level 3
*          PBLAS operation.
*
*  TRANSB  (global input) CHARACTER*1
*          On entry,  TRANSB  specifies the TRANSB option in the Level 3
*          PBLAS operation.
*
*  DIAG    (global input) CHARACTER*1
*          On entry, DIAG specifies the DIAG option in the Level 3 PBLAS
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
*  K       (global input) INTEGER
*          On entry,  K  specifies  the  dimension of the submatrix ope-
*          rands.
*
*  ALPHA   (global input) COMPLEX*16
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
*  IB      (global input) INTEGER
*          On entry, IB  specifies B's global row index, which points to
*          the beginning of the submatrix sub( B ).
*
*  JB      (global input) INTEGER
*          On entry, JB  specifies B's global column index, which points
*          to the beginning of the submatrix sub( B ).
*
*  DESCB   (global and local input) INTEGER array
*          On entry, DESCB  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix B.
*
*  BETA    (global input) COMPLEX*16
*          On entry, BETA specifies the scalar beta.
*
*  IC      (global input) INTEGER
*          On entry, IC  specifies C's global row index, which points to
*          the beginning of the submatrix sub( C ).
*
*  JC      (global input) INTEGER
*          On entry, JC  specifies C's global column index, which points
*          to the beginning of the submatrix sub( C ).
*
*  DESCC   (global and local input) INTEGER array
*          On entry, DESCC  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix C.
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
      CHARACTER*1        DIAGREF, SIDEREF, TRANSAREF, TRANSBREF, UPLOREF
      INTEGER            I, IAREF, IBREF, ICREF, JAREF, JBREF, JCREF,
     $                   KREF, MREF, MYCOL, MYROW, NPCOL, NPROW, NREF
      COMPLEX*16         ALPHAREF, BETAREF
*     ..
*     .. Local Arrays ..
      CHARACTER*15       ARGNAME
      INTEGER            DESCAREF( DLEN_ ), DESCBREF( DLEN_ ),
     $                   DESCCREF( DLEN_ )
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
         DIAGREF   = DIAG
         SIDEREF   = SIDE
         TRANSAREF = TRANSA
         TRANSBREF = TRANSB
         UPLOREF   = UPLO
         MREF      = M
         NREF      = N
         KREF      = K
         ALPHAREF  = ALPHA
         IAREF     = IA
         JAREF     = JA
         DO 10 I = 1, DLEN_
            DESCAREF( I ) = DESCA( I )
   10    CONTINUE
         IBREF     = IB
         JBREF     = JB
         DO 20 I = 1, DLEN_
            DESCBREF( I ) = DESCB( I )
   20    CONTINUE
         BETAREF   = BETA
         ICREF     = IC
         JCREF     = JC
         DO 30 I = 1, DLEN_
            DESCCREF( I ) = DESCC( I )
   30    CONTINUE
*
      ELSE
*
*        Test saved args. Return with first mismatch.
*
         ARGNAME = ' '
         IF( .NOT. LSAME( DIAG, DIAGREF ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DIAG'
         ELSE IF( .NOT. LSAME( SIDE, SIDEREF ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'SIDE'
         ELSE IF( .NOT. LSAME( TRANSA, TRANSAREF ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'TRANSA'
         ELSE IF( .NOT. LSAME( TRANSB, TRANSBREF ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'TRANSB'
         ELSE IF( .NOT. LSAME( UPLO, UPLOREF ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'UPLO'
         ELSE IF( M.NE.MREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'M'
         ELSE IF( N.NE.NREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'N'
         ELSE IF( K.NE.KREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'K'
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
         ELSE IF( IB.NE.IBREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'IB'
         ELSE IF( JB.NE.JBREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'JB'
         ELSE IF( DESCB( DTYPE_ ).NE.DESCBREF( DTYPE_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCB( DTYPE_ )'
         ELSE IF( DESCB( M_ ).NE.DESCBREF( M_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCB( M_ )'
         ELSE IF( DESCB( N_ ).NE.DESCBREF( N_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCB( N_ )'
         ELSE IF( DESCB( IMB_ ).NE.DESCBREF( IMB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCB( IMB_ )'
         ELSE IF( DESCB( INB_ ).NE.DESCBREF( INB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCB( INB_ )'
         ELSE IF( DESCB( MB_ ).NE.DESCBREF( MB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCB( MB_ )'
         ELSE IF( DESCB( NB_ ).NE.DESCBREF( NB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCB( NB_ )'
         ELSE IF( DESCB( RSRC_ ).NE.DESCBREF( RSRC_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCB( RSRC_ )'
         ELSE IF( DESCB( CSRC_ ).NE.DESCBREF( CSRC_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCB( CSRC_ )'
         ELSE IF( DESCB( CTXT_ ).NE.DESCBREF( CTXT_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCB( CTXT_ )'
         ELSE IF( DESCB( LLD_ ).NE.DESCBREF( LLD_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCB( LLD_ )'
         ELSE IF( BETA.NE.BETAREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'BETA'
         ELSE IF( IC.NE.ICREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'IC'
         ELSE IF( JC.NE.JCREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'JC'
         ELSE IF( DESCC( DTYPE_ ).NE.DESCCREF( DTYPE_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCC( DTYPE_ )'
         ELSE IF( DESCC( M_ ).NE.DESCCREF( M_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCC( M_ )'
         ELSE IF( DESCC( N_ ).NE.DESCCREF( N_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCC( N_ )'
         ELSE IF( DESCC( IMB_ ).NE.DESCCREF( IMB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCC( IMB_ )'
         ELSE IF( DESCC( INB_ ).NE.DESCCREF( INB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCC( INB_ )'
         ELSE IF( DESCC( MB_ ).NE.DESCCREF( MB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCC( MB_ )'
         ELSE IF( DESCC( NB_ ).NE.DESCCREF( NB_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCC( NB_ )'
         ELSE IF( DESCC( RSRC_ ).NE.DESCCREF( RSRC_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCC( RSRC_ )'
         ELSE IF( DESCC( CSRC_ ).NE.DESCCREF( CSRC_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCC( CSRC_ )'
         ELSE IF( DESCC( CTXT_ ).NE.DESCCREF( CTXT_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCC( CTXT_ )'
         ELSE IF( DESCC( LLD_ ).NE.DESCCREF( LLD_ ) ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'DESCC( LLD_ )'
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
*     End of PZCHKARG3
*
      END
      SUBROUTINE PZBLAS3TSTCHK( ICTXT, NOUT, NROUT, SIDE, UPLO, TRANSA,
     $                          TRANSB, DIAG, M, N, K, ALPHA, A, PA, IA,
     $                          JA, DESCA, B, PB, IB, JB, DESCB, BETA,
     $                          C, PC, IC, JC, DESCC, THRESH, ROGUE,
     $                          WORK, RWORK, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        DIAG, SIDE, TRANSA, TRANSB, UPLO
      INTEGER            IA, IB, IC, ICTXT, INFO, JA, JB, JC, K, M, N,
     $                   NOUT, NROUT
      REAL THRESH
      COMPLEX*16         ALPHA, BETA, ROGUE
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCC( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( * ), B( * ), C( * ), PA( * ), PB( * ),
     $                   PC( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZBLAS3TSTCHK performs the computational tests of the Level 3 PBLAS.
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
*             If NROUT =  1, PZGEMM  will be tested;
*             else if NROUT =  2, PZSYMM  will be tested;
*             else if NROUT =  3, PZHEMM  will be tested;
*             else if NROUT =  4, PZSYRK  will be tested;
*             else if NROUT =  5, PZHERK  will be tested;
*             else if NROUT =  6, PZSYR2K will be tested;
*             else if NROUT =  7, PZHER2K will be tested;
*             else if NROUT =  8, PZTRMM  will be tested;
*             else if NROUT =  9, PZTRSM  will be tested;
*             else if NROUT = 10, PZGEADD will be tested;
*             else if NROUT = 11, PZTRADD will be tested;
*
*  SIDE    (global input) CHARACTER*1
*          On entry, SIDE specifies if the multiplication should be per-
*          formed from the left or the right.
*
*  UPLO    (global input) CHARACTER*1
*          On entry,  UPLO  specifies  if the upper or lower part of the
*          matrix operand is to be referenced.
*
*  TRANSA  (global input) CHARACTER*1
*          On entry, TRANSA specifies if the matrix  operand  A is to be
*          transposed.
*
*  TRANSB  (global input) CHARACTER*1
*          On entry, TRANSB specifies if the matrix  operand  B is to be
*          transposed.
*
*  DIAG    (global input) CHARACTER*1
*          On entry, DIAG specifies if the  triangular matrix operand is
*          unit or non-unit.
*
*  M       (global input) INTEGER
*          On entry, M specifies the number of rows of C.
*
*  N       (global input) INTEGER
*          On entry, N specifies the number of columns of C.
*
*  K       (global input) INTEGER
*          On entry, K specifies the number of columns (resp. rows) of A
*          when  TRANSA = 'N'  (resp. TRANSA <> 'N')  in PxGEMM, PxSYRK,
*          PxSYR2K, PxHERK and PxHER2K.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (local input/local output) COMPLEX*16 array
*          On entry, A is an array of  dimension  (DESCA( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PA.
*
*  PA      (local input) COMPLEX*16 array
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
*  B       (local input/local output) COMPLEX*16 array
*          On entry, B is an array of  dimension  (DESCB( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PB.
*
*  PB      (local input) COMPLEX*16 array
*          On entry, PB is an array of dimension (DESCB( LLD_ ),*). This
*          array contains the local entries of the matrix PB.
*
*  IB      (global input) INTEGER
*          On entry, IB  specifies B's global row index, which points to
*          the beginning of the submatrix sub( B ).
*
*  JB      (global input) INTEGER
*          On entry, JB  specifies B's global column index, which points
*          to the beginning of the submatrix sub( B ).
*
*  DESCB   (global and local input) INTEGER array
*          On entry, DESCB  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix B.
*
*  BETA    (global input) COMPLEX*16
*          On entry, BETA specifies the scalar beta.
*
*  C       (local input/local output) COMPLEX*16 array
*          On entry, C is an array of  dimension  (DESCC( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PC.
*
*  PC      (local input) COMPLEX*16 array
*          On entry, PC is an array of dimension (DESCC( LLD_ ),*). This
*          array contains the local pieces of the matrix PC.
*
*  IC      (global input) INTEGER
*          On entry, IC  specifies C's global row index, which points to
*          the beginning of the submatrix sub( C ).
*
*  JC      (global input) INTEGER
*          On entry, JC  specifies C's global column index, which points
*          to the beginning of the submatrix sub( C ).
*
*  DESCC   (global and local input) INTEGER array
*          On entry, DESCC  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix C.
*
*  THRESH  (global input) REAL
*          On entry, THRESH is the threshold value for the test ratio.
*
*  ROGUE   (global input) COMPLEX*16
*          On entry,  ROGUE  specifies  the  constant  used  to  pad the
*          non-referenced part of triangular, symmetric or Hermitian ma-
*          trices.
*
*  WORK    (workspace) COMPLEX*16 array
*          On entry, WORK  is an array of dimension LWORK where LWORK is
*          at least MAX( M, MAX( N, K ) ).  This  array is used to store
*          a copy of a column of C (see PZMMCH).
*
*  RWORK   (workspace) DOUBLE PRECISION array
*          On entry, RWORK  is an array of dimension LRWORK where LRWORK
*          is at least MAX( M, MAX( N, K ) ). This array is used to sto-
*          re the computed gauges (see PZMMCH).
*
*  INFO    (global output) INTEGER
*          On exit, if INFO = 0, no  error  has  been  found,  otherwise
*          if( MOD( INFO,   2 ) = 1 ) then an error on A has been found,
*          if( MOD( INFO/2, 2 ) = 1 ) then an error on B has been found,
*          if( MOD( INFO/4, 2 ) = 1 ) then an error on C has been found.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   RZERO
      PARAMETER          ( RZERO = 0.0D+0 )
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                     ZERO = ( 0.0D+0, 0.0D+0 ) )
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
      COMPLEX*16         ALPHA1, BETA1
*     ..
*     .. Local Arrays ..
      INTEGER            IERR( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PB_ZLASET, PZCHKMIN, PZMMCH,
     $                   PZMMCH1, PZMMCH2, PZMMCH3, PZTRMM, ZTRSM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX
*     ..
*     .. Executable Statements ..
*
      INFO    = 0
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
*        Test PZGEMM
*
*        Check the resulting matrix C
*
         CALL PZMMCH( ICTXT, TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA,
     $                DESCA, B, IB, JB, DESCB, BETA, C, PC, IC, JC,
     $                DESCC, WORK, RWORK, ERR, IERR( 3 ) )
*
         IF( IERR( 3 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 ) ERR
         END IF
*
*        Check the input-only arguments
*
         IF( LSAME( TRANSA, 'N' ) ) THEN
            CALL PZCHKMIN( ERR, M, K, A, PA, IA, JA, DESCA, IERR( 1 ) )
         ELSE
            CALL PZCHKMIN( ERR, K, M, A, PA, IA, JA, DESCA, IERR( 1 ) )
         END IF
         IF( LSAME( TRANSB, 'N' ) ) THEN
            CALL PZCHKMIN( ERR, K, N, B, PB, IB, JB, DESCB, IERR( 2 ) )
         ELSE
            CALL PZCHKMIN( ERR, N, K, B, PB, IB, JB, DESCB, IERR( 2 ) )
         END IF
*
      ELSE IF( NROUT.EQ.2 ) THEN
*
*        Test PZSYMM
*
*        Check the resulting matrix C
*
         IF( LSAME( SIDE, 'L' ) ) THEN
            CALL PZMMCH( ICTXT, 'No transpose', 'No transpose', M, N, M,
     $                   ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB,
     $                   BETA, C, PC, IC, JC, DESCC, WORK, RWORK, ERR,
     $                   IERR( 3 ) )
         ELSE
            CALL PZMMCH( ICTXT, 'No transpose', 'No transpose', M, N, N,
     $                   ALPHA, B, IB, JB, DESCB, A, IA, JA, DESCA,
     $                   BETA, C, PC, IC, JC, DESCC, WORK, RWORK, ERR,
     $                   IERR( 3 ) )
         END IF
*
         IF( IERR( 3 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 ) ERR
         END IF
*
*        Check the input-only arguments
*
         IF( LSAME( UPLO, 'L' ) ) THEN
            IF( LSAME( SIDE, 'L' ) ) THEN
               CALL PB_ZLASET( 'Upper', M-1, M-1, 0, ROGUE, ROGUE,
     $                         A( IA+JA*DESCA( M_ ) ), DESCA( M_ ) )
            ELSE
               CALL PB_ZLASET( 'Upper', N-1, N-1, 0, ROGUE, ROGUE,
     $                         A( IA+JA*DESCA( M_ ) ), DESCA( M_ ) )
            END IF
         ELSE
            IF( LSAME( SIDE, 'L' ) ) THEN
               CALL PB_ZLASET( 'Lower', M-1, M-1, 0, ROGUE, ROGUE,
     $                         A( IA+1+(JA-1)*DESCA( M_ ) ),
     $                         DESCA( M_ ) )
            ELSE
               CALL PB_ZLASET( 'Lower', N-1, N-1, 0, ROGUE, ROGUE,
     $                         A( IA+1+(JA-1)*DESCA( M_ ) ),
     $                         DESCA( M_ ) )
            END IF
         END IF
*
         IF( LSAME( SIDE, 'L' ) ) THEN
            CALL PZCHKMIN( ERR, M, M, A, PA, IA, JA, DESCA, IERR( 1 ) )
         ELSE
            CALL PZCHKMIN( ERR, N, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
         END IF
         CALL PZCHKMIN( ERR, M, N, B, PB, IB, JB, DESCB, IERR( 2 ) )
*
      ELSE IF( NROUT.EQ.3 ) THEN
*
*        Test PZHEMM
*
*        Check the resulting matrix C
*
         IF( LSAME( SIDE, 'L' ) ) THEN
            CALL PZMMCH( ICTXT, 'No transpose', 'No transpose', M, N, M,
     $                   ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB,
     $                   BETA, C, PC, IC, JC, DESCC, WORK, RWORK, ERR,
     $                   IERR( 3 ) )
         ELSE
            CALL PZMMCH( ICTXT, 'No transpose', 'No transpose', M, N, N,
     $                   ALPHA, B, IB, JB, DESCB, A, IA, JA, DESCA,
     $                   BETA, C, PC, IC, JC, DESCC, WORK, RWORK, ERR,
     $                   IERR( 3 ) )
         END IF
*
         IF( IERR( 3 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 ) ERR
         END IF
*
*        Check the input-only arguments
*
         IF( LSAME( UPLO, 'L' ) ) THEN
            IF( LSAME( SIDE, 'L' ) ) THEN
               CALL PB_ZLASET( 'Upper', M-1, M-1, 0, ROGUE, ROGUE,
     $                         A( IA+JA*DESCA( M_ ) ), DESCA( M_ ) )
            ELSE
               CALL PB_ZLASET( 'Upper', N-1, N-1, 0, ROGUE, ROGUE,
     $                         A( IA+JA*DESCA( M_ ) ), DESCA( M_ ) )
            END IF
         ELSE
            IF( LSAME( SIDE, 'L' ) ) THEN
               CALL PB_ZLASET( 'Lower', M-1, M-1, 0, ROGUE, ROGUE,
     $                         A( IA+1+(JA-1)*DESCA( M_ ) ),
     $                         DESCA( M_ ) )
            ELSE
               CALL PB_ZLASET( 'Lower', N-1, N-1, 0, ROGUE, ROGUE,
     $                         A( IA+1+(JA-1)*DESCA( M_ ) ),
     $                         DESCA( M_ ) )
            END IF
         END IF
*
         IF( LSAME( SIDE, 'L' ) ) THEN
            CALL PZCHKMIN( ERR, M, M, A, PA, IA, JA, DESCA, IERR( 1 ) )
         ELSE
            CALL PZCHKMIN( ERR, N, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
         END IF
         CALL PZCHKMIN( ERR, M, N, B, PB, IB, JB, DESCB, IERR( 2 ) )
*
      ELSE IF( NROUT.EQ.4 ) THEN
*
*        Test PZSYRK
*
*        Check the resulting matrix C
*
         IF( LSAME( TRANSA, 'N' ) ) THEN
            CALL PZMMCH1( ICTXT, UPLO, 'No transpose', N, K, ALPHA, A,
     $                    IA, JA, DESCA, BETA, C, PC, IC, JC, DESCC,
     $                    WORK, RWORK, ERR, IERR( 3 ) )
         ELSE
            CALL PZMMCH1( ICTXT, UPLO, 'Transpose', N, K, ALPHA, A, IA,
     $                    JA, DESCA, BETA, C, PC, IC, JC, DESCC, WORK,
     $                    RWORK, ERR, IERR( 3 ) )
         END IF
*
         IF( IERR( 3 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 ) ERR
         END IF
*
*        Check the input-only arguments
*
         IF( LSAME( TRANSA, 'N' ) ) THEN
            CALL PZCHKMIN( ERR, N, K, A, PA, IA, JA, DESCA, IERR( 1 ) )
         ELSE
            CALL PZCHKMIN( ERR, K, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
         END IF
*
      ELSE IF( NROUT.EQ.5 ) THEN
*
*        Test PZHERK
*
*        Check the resulting matrix C
*
         BETA1  = DCMPLX( DBLE( BETA ), RZERO )
         ALPHA1 = DCMPLX( DBLE( ALPHA ), RZERO )
         IF( LSAME( TRANSA, 'N' ) ) THEN
            CALL PZMMCH1( ICTXT, UPLO, 'Hermitian', N, K, ALPHA1, A, IA,
     $                    JA, DESCA, BETA1, C, PC, IC, JC, DESCC, WORK,
     $                    RWORK, ERR, IERR( 3 ) )
         ELSE
            CALL PZMMCH1( ICTXT, UPLO, 'Conjugate transpose', N, K,
     $                    ALPHA1, A, IA, JA, DESCA, BETA1, C, PC, IC,
     $                    JC, DESCC, WORK, RWORK, ERR, IERR( 3 ) )
         END IF
*
         IF( IERR( 3 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 ) ERR
         END IF
*
*        Check the input-only arguments
*
         IF( LSAME( TRANSA, 'N' ) ) THEN
            CALL PZCHKMIN( ERR, N, K, A, PA, IA, JA, DESCA, IERR( 1 ) )
         ELSE
            CALL PZCHKMIN( ERR, K, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
         END IF
*
      ELSE IF( NROUT.EQ.6 ) THEN
*
*        Test PZSYR2K
*
*        Check the resulting matrix C
*
         IF( LSAME( TRANSA, 'N' ) ) THEN
            CALL PZMMCH2( ICTXT, UPLO, 'No transpose', N, K, ALPHA, A,
     $                    IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, PC,
     $                    IC, JC, DESCC, WORK, RWORK, ERR, IERR( 3 ) )
         ELSE
            CALL PZMMCH2( ICTXT, UPLO, 'Transpose', N, K, ALPHA, A,
     $                    IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, PC,
     $                    IC, JC, DESCC, WORK, RWORK, ERR,
     $                    IERR( 3 ) )
         END IF
*
         IF( IERR( 3 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 ) ERR
         END IF
*
*        Check the input-only arguments
*
         IF( LSAME( TRANSA, 'N' ) ) THEN
            CALL PZCHKMIN( ERR, N, K, A, PA, IA, JA, DESCA, IERR( 1 ) )
            CALL PZCHKMIN( ERR, N, K, B, PB, IB, JB, DESCB, IERR( 2 ) )
         ELSE
            CALL PZCHKMIN( ERR, K, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
            CALL PZCHKMIN( ERR, K, N, B, PB, IB, JB, DESCB, IERR( 2 ) )
         END IF
*
      ELSE IF( NROUT.EQ.7 ) THEN
*
*        Test PZHER2K
*
*        Check the resulting matrix C
*
         BETA1 = DCMPLX( DBLE( BETA ), RZERO )
         IF( LSAME( TRANSA, 'N' ) ) THEN
            CALL PZMMCH2( ICTXT, UPLO, 'Hermitian', N, K, ALPHA, A, IA,
     $                    JA, DESCA, B, IB, JB, DESCB, BETA1, C, PC, IC,
     $                    JC, DESCC, WORK, RWORK, ERR, IERR( 3 ) )
         ELSE
            CALL PZMMCH2( ICTXT, UPLO, 'Conjugate transpose', N, K,
     $                    ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB,
     $                    BETA1, C, PC, IC, JC, DESCC, WORK, RWORK, ERR,
     $                    IERR( 3 ) )
         END IF
*
         IF( IERR( 3 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 ) ERR
         END IF
*
*        Check the input-only arguments
*
         IF( LSAME( TRANSA, 'N' ) ) THEN
            CALL PZCHKMIN( ERR, N, K, A, PA, IA, JA, DESCA, IERR( 1 ) )
            CALL PZCHKMIN( ERR, N, K, B, PB, IB, JB, DESCB, IERR( 2 ) )
         ELSE
            CALL PZCHKMIN( ERR, K, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
            CALL PZCHKMIN( ERR, K, N, B, PB, IB, JB, DESCB, IERR( 2 ) )
         END IF
*
      ELSE IF( NROUT.EQ.8 ) THEN
*
*        Test PZTRMM
*
*        Check the resulting matrix B
*
         IF( LSAME( SIDE, 'L' ) ) THEN
            CALL PZMMCH( ICTXT, TRANSA, 'No transpose', M, N, M,
     $                   ALPHA, A, IA, JA, DESCA, C, IB, JB, DESCB,
     $                   ZERO, B, PB, IB, JB, DESCB, WORK, RWORK, ERR,
     $                   IERR( 2 ) )
         ELSE
            CALL PZMMCH( ICTXT, 'No transpose', TRANSA, M, N, N,
     $                   ALPHA, C, IB, JB, DESCB, A, IA, JA, DESCA,
     $                   ZERO, B, PB, IB, JB, DESCB, WORK, RWORK, ERR,
     $                   IERR( 2 ) )
         END IF
*
         IF( IERR( 2 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 ) ERR
         END IF
*
*        Check the input-only arguments
*
         IF( LSAME( SIDE, 'L' ) ) THEN
            IF( LSAME( UPLO, 'L' ) ) THEN
               IF( LSAME( DIAG, 'N' ) ) THEN
                  CALL PB_ZLASET( 'Upper', M-1, M-1, 0, ROGUE, ROGUE,
     $                            A( IA+JA*DESCA( M_ ) ), DESCA( M_ ) )
               ELSE
                  CALL PB_ZLASET( 'Upper', M, M, 0, ROGUE, ONE,
     $                            A( IA+(JA-1)*DESCA( M_ ) ),
     $                            DESCA( M_ ) )
               END IF
            ELSE
               IF( LSAME( DIAG, 'N' ) ) THEN
                  CALL PB_ZLASET( 'Lower', M-1, M-1, 0, ROGUE, ROGUE,
     $                            A( IA+1+(JA-1)*DESCA( M_ ) ),
     $                            DESCA( M_ ) )
               ELSE
                  CALL PB_ZLASET( 'Lower', M, M, 0, ROGUE, ONE,
     $                            A( IA+(JA-1)*DESCA( M_ ) ),
     $                            DESCA( M_ ) )
               END IF
            END IF
            CALL PZCHKMIN( ERR, M, M, A, PA, IA, JA, DESCA, IERR( 1 ) )
         ELSE
            IF( LSAME( UPLO, 'L' ) ) THEN
               IF( LSAME( DIAG, 'N' ) ) THEN
                  CALL PB_ZLASET( 'Upper', N-1, N-1, 0, ROGUE, ROGUE,
     $                            A( IA+JA*DESCA( M_ ) ), DESCA( M_ ) )
               ELSE
                  CALL PB_ZLASET( 'Upper', N, N, 0, ROGUE, ONE,
     $                            A( IA+(JA-1)*DESCA( M_ ) ),
     $                            DESCA( M_ ) )
               END IF
            ELSE
               IF( LSAME( DIAG, 'N' ) ) THEN
                  CALL PB_ZLASET( 'Lower', N-1, N-1, 0, ROGUE, ROGUE,
     $                            A( IA+1+(JA-1)*DESCA( M_ ) ),
     $                            DESCA( M_ ) )
               ELSE
                  CALL PB_ZLASET( 'Lower', N, N, 0, ROGUE, ONE,
     $                            A( IA+(JA-1)*DESCA( M_ ) ),
     $                            DESCA( M_ ) )
               END IF
            END IF
            CALL PZCHKMIN( ERR, N, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
         END IF
*
      ELSE IF( NROUT.EQ.9 ) THEN
*
*        Test PZTRSM
*
*        Check the resulting matrix B
*
         CALL ZTRSM( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,
     $               A( IA+(JA-1)*DESCA( M_ ) ), DESCA( M_ ),
     $               B( IB+(JB-1)*DESCB( M_ ) ), DESCB( M_ ) )
         CALL PZTRMM( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, PA, IA, JA,
     $                DESCA, PB, IB, JB, DESCB )
         IF( LSAME( SIDE, 'L' ) ) THEN
            CALL PZMMCH( ICTXT, TRANSA, 'No transpose', M, N, M, ALPHA,
     $                   A, IA, JA, DESCA, B, IB, JB, DESCB, ZERO, C,
     $                   PB, IB, JB, DESCB, WORK, RWORK, ERR,
     $                   IERR( 2 ) )
         ELSE
            CALL PZMMCH( ICTXT, 'No transpose', TRANSA, M, N, N, ALPHA,
     $                   B, IB, JB, DESCB, A, IA, JA, DESCA, ZERO, C,
     $                   PB, IB, JB, DESCB, WORK, RWORK, ERR,
     $                   IERR( 2 ) )
         END IF
*
         IF( IERR( 2 ).NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9998 )
         ELSE IF( ERR.GT.DBLE( THRESH ) ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $         WRITE( NOUT, FMT = 9997 ) ERR
         END IF
*
*        Check the input-only arguments
*
         IF( LSAME( SIDE, 'L' ) ) THEN
            IF( LSAME( UPLO, 'L' ) ) THEN
               IF( LSAME( DIAG, 'N' ) ) THEN
                  CALL PB_ZLASET( 'Upper', M-1, M-1, 0, ROGUE, ROGUE,
     $                            A( IA+JA*DESCA( M_ ) ), DESCA( M_ ) )
               ELSE
                  CALL PB_ZLASET( 'Upper', M, M, 0, ROGUE, ONE,
     $                            A( IA+(JA-1)*DESCA( M_ ) ),
     $                            DESCA( M_ ) )
               END IF
            ELSE
               IF( LSAME( DIAG, 'N' ) ) THEN
                  CALL PB_ZLASET( 'Lower', M-1, M-1, 0, ROGUE, ROGUE,
     $                            A( IA+1+(JA-1)*DESCA( M_ ) ),
     $                            DESCA( M_ ) )
               ELSE
                  CALL PB_ZLASET( 'Lower', M, M, 0, ROGUE, ONE,
     $                            A( IA+(JA-1)*DESCA( M_ ) ),
     $                            DESCA( M_ ) )
               END IF
            END IF
            CALL PZCHKMIN( ERR, M, M, A, PA, IA, JA, DESCA, IERR( 1 ) )
         ELSE
            IF( LSAME( UPLO, 'L' ) ) THEN
               IF( LSAME( DIAG, 'N' ) ) THEN
                  CALL PB_ZLASET( 'Upper', N-1, N-1, 0, ROGUE, ROGUE,
     $                            A( IA+JA*DESCA( M_ ) ), DESCA( M_ ) )
               ELSE
                  CALL PB_ZLASET( 'Upper', N, N, 0, ROGUE, ONE,
     $                            A( IA+(JA-1)*DESCA( M_ ) ),
     $                            DESCA( M_ ) )
               END IF
            ELSE
               IF( LSAME( DIAG, 'N' ) ) THEN
                  CALL PB_ZLASET( 'Lower', N-1, N-1, 0, ROGUE, ROGUE,
     $                            A( IA+1+(JA-1)*DESCA( M_ ) ),
     $                            DESCA( M_ ) )
               ELSE
                  CALL PB_ZLASET( 'Lower', N, N, 0, ROGUE, ONE,
     $                            A( IA+(JA-1)*DESCA( M_ ) ),
     $                            DESCA( M_ ) )
               END IF
            END IF
            CALL PZCHKMIN( ERR, N, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
         END IF
      ELSE IF( NROUT.EQ.10 ) THEN
*
*        Test PZGEADD
*
*        Check the resulting matrix C
*
         CALL PZMMCH3( 'All', TRANSA, M, N, ALPHA, A, IA, JA, DESCA,
     $                 BETA, C, PC, IC, JC, DESCC, ERR, IERR( 3 ) )
*
*        Check the input-only arguments
*
         IF( LSAME( TRANSA, 'N' ) ) THEN
            CALL PZCHKMIN( ERR, M, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
         ELSE
            CALL PZCHKMIN( ERR, N, M, A, PA, IA, JA, DESCA, IERR( 1 ) )
         END IF
*
      ELSE IF( NROUT.EQ.11 ) THEN
*
*        Test PZTRADD
*
*        Check the resulting matrix C
*
         CALL PZMMCH3( UPLO, TRANSA, M, N, ALPHA, A, IA, JA, DESCA,
     $                 BETA, C, PC, IC, JC, DESCC, ERR, IERR( 3 ) )
*
*        Check the input-only arguments
*
         IF( LSAME( TRANSA, 'N' ) ) THEN
            CALL PZCHKMIN( ERR, M, N, A, PA, IA, JA, DESCA, IERR( 1 ) )
         ELSE
            CALL PZCHKMIN( ERR, N, M, A, PA, IA, JA, DESCA, IERR( 1 ) )
         END IF
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
     $      WRITE( NOUT, FMT = 9999 ) 'B'
      END IF
*
      IF( IERR( 3 ).NE.0 ) THEN
         INFO = INFO + 4
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9999 ) 'C'
      END IF
*
 9999 FORMAT( 2X, '   ***** ERROR: Matrix operand ', A,
     $        ' is incorrect.' )
 9998 FORMAT( 2X, '   ***** FATAL ERROR - Computed result is less ',
     $        'than half accurate *****' )
 9997 FORMAT( 2X, '   ***** Test completed with maximum test ratio: ',
     $        F11.5, ' SUSPECT *****' )
*
      RETURN
*
*     End of PZBLAS3TSTCHK
*
      END
