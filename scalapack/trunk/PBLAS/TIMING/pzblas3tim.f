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
      
      PROGRAM PZBLA3TIM
*
*  -- PBLAS timing driver (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*  Purpose
*  =======
*
*  PZBLA3TIM  is the main timing program for the Level 3 PBLAS routines.
*
*  The program must be driven by a short data file.  An  annotated exam-
*  ple of a data file can be obtained by deleting the first 3 characters
*  from the following 59 lines:
*  'Level 3 PBLAS, Timing input file'
*  'Intel iPSC/860 hypercube, gamma model.'
*  'PZBLAS3TIM.SUMM'     output file name (if any)
*  6     device out
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
      INTEGER            MAXTESTS, MAXGRIDS, ZPLXSZ, TOTMEM, MEMSIZ,
     $                   NSUBS
      COMPLEX*16         ONE
      PARAMETER          ( MAXTESTS = 20, MAXGRIDS = 20, ZPLXSZ = 16,
     $                   ONE = ( 1.0D+0, 0.0D+0 ), TOTMEM = 2000000,
     $                   NSUBS = 11, MEMSIZ = TOTMEM / ZPLXSZ )
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      CHARACTER*1        ADIAGDO, AFORM, CFORM, DIAG, SIDE, TRANSA,
     $                   TRANSB, UPLO
      INTEGER            CSRCA, CSRCB, CSRCC, I, IA, IAM, IASEED, IB,
     $                   IBSEED, IC, ICSEED, ICTXT, IMBA, IMBB, IMBC,
     $                   IMIDA, IMIDB, IMIDC, INBA, INBB, INBC, IPA,
     $                   IPB, IPC, IPOSTA, IPOSTB, IPOSTC, IPREA, IPREB,
     $                   IPREC, J, JA, JB, JC, K, L, M, MA, MB, MBA,
     $                   MBB, MBC, MC, MEMREQD, MPA, MPB, MPC, MYCOL,
     $                   MYROW, N, NA, NB, NBA, NBB, NBC, NC, NCOLA,
     $                   NCOLB, NCOLC, NGRIDS, NOUT, NPCOL, NPROCS,
     $                   NPROW, NQA, NQB, NQC, NROWA, NROWB, NROWC,
     $                   NTESTS, OFFDA, OFFDC, RSRCA, RSRCB, RSRCC
      DOUBLE PRECISION   CFLOPS, NOPS, WFLOPS
      COMPLEX*16         ALPHA, BETA, SCALE
*     ..
*     .. Local Arrays ..
      LOGICAL            LTEST( NSUBS ), BCHECK( NSUBS ),
     $                   CCHECK( NSUBS )
      CHARACTER*1        DIAGVAL( MAXTESTS ), SIDEVAL( MAXTESTS ),
     $                   TRNAVAL( MAXTESTS ), TRNBVAL( MAXTESTS ),
     $                   UPLOVAL( MAXTESTS )
      CHARACTER*80       OUTFILE
      INTEGER            CSCAVAL( MAXTESTS ), CSCBVAL( MAXTESTS ),
     $                   CSCCVAL( MAXTESTS ), DESCA( DLEN_ ),
     $                   DESCB( DLEN_ ), DESCC( DLEN_ ),
     $                   IAVAL( MAXTESTS ), IBVAL( MAXTESTS ),
     $                   ICVAL( MAXTESTS ), IERR( 3 ),
     $                   IMBAVAL( MAXTESTS ), IMBBVAL( MAXTESTS ),
     $                   IMBCVAL( MAXTESTS ), INBAVAL( MAXTESTS ),
     $                   INBBVAL( MAXTESTS ), INBCVAL( MAXTESTS ),
     $                   JAVAL( MAXTESTS ), JBVAL( MAXTESTS ),
     $                   JCVAL( MAXTESTS ), KVAL( MAXTESTS ),
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
      DOUBLE PRECISION   CTIME( 1 ), WTIME( 1 )
      COMPLEX*16         MEM( MEMSIZ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_EXIT, BLACS_GET,
     $                   BLACS_GRIDEXIT, BLACS_GRIDINFO, BLACS_GRIDINIT,
     $                   BLACS_PINFO, IGSUM2D, PB_BOOT, PB_COMBINE,
     $                   PB_TIMER, PMDESCCHK, PMDIMCHK, PZBLA3TIMINFO,
     $                   PZGEADD, PZGEMM, PZHEMM, PZHER2K, PZHERK,
     $                   PZLAGEN, PZLASCAL, PZSYMM, PZSYR2K, PZSYRK,
     $                   PZTRADD, PZTRMM, PZTRSM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PDOPBL3
      EXTERNAL           LSAME, PDOPBL3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX
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
*     Set flag so that the PBLAS error handler won't abort on errors, so
*     that the tester will detect unsupported operations.
*
      ABRTFLG = .FALSE.
*
*     Seeds for random matrix generations.
*
      IASEED = 100
      IBSEED = 200
      ICSEED = 300
*
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
      CALL PZBLA3TIMINFO( OUTFILE, NOUT, NTESTS, DIAGVAL, SIDEVAL,
     $                    TRNAVAL, TRNBVAL, UPLOVAL, MVAL, NVAL,
     $                    KVAL, MAVAL, NAVAL, IMBAVAL, MBAVAL,
     $                    INBAVAL, NBAVAL, RSCAVAL, CSCAVAL, IAVAL,
     $                    JAVAL, MBVAL, NBVAL, IMBBVAL, MBBVAL,
     $                    INBBVAL, NBBVAL, RSCBVAL, CSCBVAL, IBVAL,
     $                    JBVAL, MCVAL, NCVAL, IMBCVAL, MBCVAL,
     $                    INBCVAL, NBCVAL, RSCCVAL, CSCCVAL, ICVAL,
     $                    JCVAL, MAXTESTS, NGRIDS, PVAL, MAXGRIDS,
     $                    QVAL, MAXGRIDS, NBLOG, LTEST, IAM, NPROCS,
     $                    ALPHA, BETA, MEM )
*
      IF( IAM.EQ.0 )
     $   WRITE( NOUT, FMT = 9984 )
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
               WRITE( NOUT, FMT = * )
               WRITE( NOUT, FMT = 9996 ) J, NPROW, NPCOL
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
               WRITE( NOUT, FMT = 9980 )
*
            END IF
*
*           Check the validity of the input test parameters
*
            IF( .NOT.LSAME( SIDE, 'L' ).AND.
     $          .NOT.LSAME( SIDE, 'R' ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'SIDE'
               GO TO 40
            END IF
*
            IF( .NOT.LSAME( UPLO, 'U' ).AND.
     $          .NOT.LSAME( UPLO, 'L' ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'UPLO'
               GO TO 40
            END IF
*
            IF( .NOT.LSAME( TRANSA, 'N' ).AND.
     $          .NOT.LSAME( TRANSA, 'T' ).AND.
     $          .NOT.LSAME( TRANSA, 'C' ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'TRANSA'
               GO TO 40
            END IF
*
            IF( .NOT.LSAME( TRANSB, 'N' ).AND.
     $          .NOT.LSAME( TRANSB, 'T' ).AND.
     $          .NOT.LSAME( TRANSB, 'C' ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'TRANSB'
               GO TO 40
            END IF
*
            IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $          .NOT.LSAME( DIAG , 'N' ) )THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'DIAG'
               GO TO 40
            END IF
*
*           Check and initialize the matrix descriptors
*
            CALL PMDESCCHK( ICTXT, NOUT, 'A', DESCA,
     $                      BLOCK_CYCLIC_2D_INB, MA, NA, IMBA, INBA,
     $                      MBA, NBA, RSRCA, CSRCA, MPA, NQA, IPREA,
     $                      IMIDA, IPOSTA, 0, 0, IERR( 1 ) )
*
            CALL PMDESCCHK( ICTXT, NOUT, 'B', DESCB,
     $                      BLOCK_CYCLIC_2D_INB, MB, NB, IMBB, INBB,
     $                      MBB, NBB, RSRCB, CSRCB, MPB, NQB, IPREB,
     $                      IMIDB, IPOSTB, 0, 0, IERR( 2 ) )
*
            CALL PMDESCCHK( ICTXT, NOUT, 'C', DESCC,
     $                      BLOCK_CYCLIC_2D_INB, MC, NC, IMBC, INBC,
     $                      MBC, NBC, RSRCC, CSRCC, MPC, NQC, IPREC,
     $                      IMIDC, IPOSTC, 0, 0, IERR( 3 ) )
*
            IF( IERR( 1 ).GT.0 .OR. IERR( 2 ).GT.0 .OR.
     $          IERR( 3 ).GT.0 ) THEN
               GO TO 40
            END IF
*
*           Assign pointers into MEM for matrices corresponding to
*           the distributed matrices A, X and Y.
*
            IPA = IPREA + 1
            IPB = IPA + DESCA( LLD_ )*NQA
            IPC = IPB + DESCB( LLD_ )*NQB
*
*           Check if sufficient memory.
*
            MEMREQD = IPC + DESCC( LLD_ )*NQC - 1
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
               ELSE IF( L.EQ.8 .OR. L.EQ.9 ) THEN
*
*                 PZTRMM, PZTRSM
*
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
     $                  WRITE( NOUT, FMT = 9983 ) SNAMES( L ), 'TRANSA'
                     GO TO 30
                  END IF
               ELSE IF( L.EQ.5 .OR. L.EQ.7 ) THEN
                  IF( .NOT.LSAME( TRANSA, 'N' ).AND.
     $                .NOT.LSAME( TRANSA, 'C' ) ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9983 ) SNAMES( L ), 'TRANSA'
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
                  AFORM = 'N'
                  ADIAGDO = 'N'
                  OFFDA = 0
                  CFORM = 'H'
                  OFFDC = IC - JC
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
               IF( ( L.EQ.9 ).AND.( .NOT.( LSAME( DIAG, 'N' ) ) ).AND.
     $             ( MAX( NROWA, NCOLA ).GT.1 ) ) THEN
                  SCALE = ONE / DCMPLX( DBLE( MAX( NROWA, NCOLA ) ) )
                  IF( LSAME( UPLO, 'L' ) ) THEN
                     CALL PZLASCAL( 'Lower', NROWA-1, NCOLA-1, SCALE,
     $                              MEM( IPA ), IA+1, JA, DESCA )
                  ELSE
                     CALL PZLASCAL( 'Upper', NROWA-1, NCOLA-1, SCALE,
     $                              MEM( IPA ), IA, JA+1, DESCA )
                  END IF
*
               END IF
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
               INFO = 0
               CALL PB_BOOT()
               CALL BLACS_BARRIER( ICTXT, 'All' )
*
*              Call the Level 3 PBLAS routine
*
               IF( L.EQ.1 ) THEN
*
*                 Test PZGEMM
*
                  NOPS = PDOPBL3( SNAMES( L ), M, N, K )
*
                  CALL PB_TIMER( 1 )
                  CALL PZGEMM( TRANSA, TRANSB, M, N, K, ALPHA,
     $                         MEM( IPA ), IA, JA, DESCA, MEM( IPB ),
     $                         IB, JB, DESCB, BETA, MEM( IPC ), IC, JC,
     $                         DESCC )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( L.EQ.2 ) THEN
*
*                 Test PZSYMM
*
                  IF( LSAME( SIDE, 'L' ) ) THEN
                     NOPS = PDOPBL3( SNAMES( L ), M, N, 0 )
                  ELSE
                     NOPS = PDOPBL3( SNAMES( L ), M, N, 1 )
                  END IF
*
                  CALL PB_TIMER( 1 )
                  CALL PZSYMM( SIDE, UPLO, M, N, ALPHA, MEM( IPA ), IA,
     $                         JA, DESCA, MEM( IPB ), IB, JB, DESCB,
     $                         BETA, MEM( IPC ), IC, JC, DESCC )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( L.EQ.3 ) THEN
*
*                 Test PZHEMM
*
                  IF( LSAME( SIDE, 'L' ) ) THEN
                     NOPS = PDOPBL3( SNAMES( L ), M, N, 0 )
                  ELSE
                     NOPS = PDOPBL3( SNAMES( L ), M, N, 1 )
                  END IF
*
                  CALL PB_TIMER( 1 )
                  CALL PZHEMM( SIDE, UPLO, M, N, ALPHA, MEM( IPA ), IA,
     $                         JA, DESCA, MEM( IPB ), IB, JB, DESCB,
     $                         BETA, MEM( IPC ), IC, JC, DESCC )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( L.EQ.4 ) THEN
*
*                 Test PZSYRK
*
                  NOPS = PDOPBL3( SNAMES( L ), N, N, K )
*
                  CALL PB_TIMER( 1 )
                  CALL PZSYRK( UPLO, TRANSA, N, K, ALPHA, MEM( IPA ),
     $                         IA, JA, DESCA, BETA, MEM( IPC ), IC, JC,
     $                         DESCC )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( L.EQ.5 ) THEN
*
*                 Test PZHERK
*
                  NOPS = PDOPBL3( SNAMES( L ), N, N, K )
*
                  CALL PB_TIMER( 1 )
                  CALL PZHERK( UPLO, TRANSA, N, K, DBLE( ALPHA ),
     $                         MEM( IPA ), IA, JA, DESCA, DBLE( BETA ),
     $                         MEM( IPC ), IC, JC, DESCC )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( L.EQ.6 ) THEN
*
*                 Test PZSYR2K
*
                  NOPS = PDOPBL3( SNAMES( L ), N, N, K )
*
                  CALL PB_TIMER( 1 )
                  CALL PZSYR2K( UPLO, TRANSA, N, K, ALPHA, MEM( IPA ),
     $                          IA, JA, DESCA, MEM( IPB ), IB, JB,
     $                          DESCB, BETA, MEM( IPC ), IC, JC,
     $                          DESCC )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( L.EQ.7 ) THEN
*
*                 Test PZHER2K
*
                  NOPS = PDOPBL3( SNAMES( L ), N, N, K )
*
                  CALL PB_TIMER( 1 )
                  CALL PZHER2K( UPLO, TRANSA, N, K, ALPHA, MEM( IPA ),
     $                          IA, JA, DESCA, MEM( IPB ), IB, JB,
     $                          DESCB, DBLE( BETA ), MEM( IPC ), IC, JC,
     $                          DESCC )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( L.EQ.8 ) THEN
*
*                 Test PZTRMM
*
                  IF( LSAME( SIDE, 'L' ) ) THEN
                     NOPS = PDOPBL3( SNAMES( L ), M, N, 0 )
                  ELSE
                     NOPS = PDOPBL3( SNAMES( L ), M, N, 1 )
                  END IF
*
                  CALL PB_TIMER( 1 )
                  CALL PZTRMM( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,
     $                         MEM( IPA ), IA, JA, DESCA, MEM( IPB ),
     $                         IB, JB, DESCB )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( L.EQ.9 ) THEN
*
*                 Test PZTRSM
*
                  IF( LSAME( SIDE, 'L' ) ) THEN
                     NOPS = PDOPBL3( SNAMES( L ), M, N, 0 )
                  ELSE
                     NOPS = PDOPBL3( SNAMES( L ), M, N, 1 )
                  END IF
*
                  CALL PB_TIMER( 1 )
                  CALL PZTRSM( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,
     $                         MEM( IPA ), IA, JA, DESCA, MEM( IPB ),
     $                         IB, JB, DESCB )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( L.EQ.10 ) THEN
*
*                 Test PZGEADD
*
                  NOPS = PDOPBL3( SNAMES( L ), M, N, M )
*
                  CALL PB_TIMER( 1 )
                  CALL PZGEADD( TRANSA, M, N, ALPHA, MEM( IPA ), IA, JA,
     $                          DESCA, BETA, MEM( IPC ), IC, JC, DESCC )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( L.EQ.11 ) THEN
*
*                 Test PZTRADD
*
                  IF( LSAME( UPLO, 'U' ) ) THEN
                     NOPS = PDOPBL3( SNAMES( L ), M, N, 0 )
                  ELSE
                     NOPS = PDOPBL3( SNAMES( L ), M, N, 1 )
                  END IF
*
                  CALL PB_TIMER( 1 )
                  CALL PZTRADD( UPLO, TRANSA, M, N, ALPHA, MEM( IPA ),
     $                          IA, JA, DESCA, BETA, MEM( IPC ), IC, JC,
     $                          DESCC )
                  CALL PB_TIMER( 1 )
*
               END IF
*
*              Check if the operation has been performed.
*
               IF( INFO.NE.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9982 ) INFO
                  GO TO 30
               END IF
*
               CALL PB_COMBINE( ICTXT, 'All', '>', 'W', 1, 1, WTIME )
               CALL PB_COMBINE( ICTXT, 'All', '>', 'C', 1, 1, CTIME )
*
*              Only node 0 prints timing test result
*
               IF( IAM.EQ.0 ) THEN
*
*                 Print WALL time if machine supports it
*
                  IF( WTIME( 1 ).GT.0.0D+0 ) THEN
                     WFLOPS = NOPS / ( WTIME( 1 ) * 1.0D+6 )
                  ELSE
                     WFLOPS = 0.0D+0
                  END IF
*
*                 Print CPU time if machine supports it
*
                  IF( CTIME( 1 ).GT.0.0D+0 ) THEN
                     CFLOPS = NOPS / ( CTIME( 1 ) * 1.0D+6 )
                  ELSE
                     CFLOPS = 0.0D+0
                  END IF
*
                  WRITE( NOUT, FMT = 9981 ) SNAMES( L ), WTIME( 1 ),
     $                                      WFLOPS, CTIME( 1 ), CFLOPS
*
               END IF
*
   30       CONTINUE
*
   40       IF( IAM.EQ.0 ) THEN
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = * )
               WRITE( NOUT, FMT = 9986 ) J
            END IF
*
   50   CONTINUE
*
        CALL BLACS_GRIDEXIT( ICTXT )
*
   60 CONTINUE
*
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9985 )
         WRITE( NOUT, FMT = * )
      END IF
*
      CALL BLACS_EXIT( 0 )
*
 9999 FORMAT( 'ILLEGAL ', A, ': ', A, ' = ', I10,
     $        ' should be at least 1' )
 9998 FORMAT( 'ILLEGAL GRID: NPROW*NPCOL = ', I4,
     $        '. It can be at most', I4 )
 9997 FORMAT( 'Bad ', A, ' parameters: going on to next test case.' )
 9996 FORMAT( 2X, 'Test number ', I2 , ' started on a ', I4, ' x ',
     $        I4, ' process grid.' )
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
 9986 FORMAT( 2X, 'Test number ', I2, ' completed.' )
 9985 FORMAT( 2X, 'End of Tests.' )
 9984 FORMAT( 2X, 'Tests started.' )
 9983 FORMAT( 5X, A, '     ***** ', A, ' has an incorrect value:     ',
     $            ' BYPASS  *****' )
 9982 FORMAT( 2X, '   ***** Operation not supported, error code: ',
     $        I5, ' *****' )
 9981 FORMAT( 2X, '|  ', A, 2X, F13.3, 2X, F13.3, 2X, F13.3, 2X, F13.3 )
 9980 FORMAT( 2X, '            WALL time (s)    WALL Mflops ',
     $        '  CPU time (s)     CPU Mflops' )
*
      STOP
*
*     End of PZBLA3TIM
*
      END
      SUBROUTINE PZBLA3TIMINFO( SUMMRY, NOUT, NMAT, DIAGVAL, SIDEVAL,
     $                          TRNAVAL, TRNBVAL, UPLOVAL, MVAL,
     $                          NVAL, KVAL, MAVAL, NAVAL, IMBAVAL,
     $                          MBAVAL, INBAVAL, NBAVAL, RSCAVAL,
     $                          CSCAVAL, IAVAL, JAVAL, MBVAL, NBVAL,
     $                          IMBBVAL, MBBVAL, INBBVAL, NBBVAL,
     $                          RSCBVAL, CSCBVAL, IBVAL, JBVAL,
     $                          MCVAL, NCVAL, IMBCVAL, MBCVAL,
     $                          INBCVAL, NBCVAL, RSCCVAL, CSCCVAL,
     $                          ICVAL, JCVAL, LDVAL, NGRIDS, PVAL,
     $                          LDPVAL, QVAL, LDQVAL, NBLOG, LTEST,
     $                          IAM, NPROCS, ALPHA, BETA, WORK )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            IAM, LDPVAL, LDQVAL, LDVAL, NBLOG, NGRIDS,
     $                   NMAT, NOUT, NPROCS
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
*  PZBLA3TIMINFO  get  the needed startup information for timing various
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
*  IAM     (local input) INTEGER
*          On entry,  IAM  specifies the number of the process executing
*          this routine.
*
*  NPROCS  (global input) INTEGER
*          On entry, NPROCS specifies the total number of processes.
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
*          MAX( 3, 2*NGRIDS+38*NMAT+NSUBS ) with NSUBS = 11. This  array
*          is  used  to  pack all output arrays in order to send info in
*          one message.
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
         OPEN( NIN, FILE='PZBLAS3TIM.dat', STATUS='OLD' )
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
*        Pack information arrays and broadcast
*
         CALL ZGEBS2D( ICTXT, 'All', ' ', 1, 1, ALPHA, 1 )
         CALL ZGEBS2D( ICTXT, 'All', ' ', 1, 1, BETA,  1 )
*
         WORK( 1 ) = NGRIDS
         WORK( 2 ) = NMAT
         WORK( 3 ) = NBLOG
         CALL IGEBS2D( ICTXT, 'All', ' ', 3, 1, WORK, 3 )
*
         I = 1
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
         WRITE( NOUT, FMT = 9999 )
     $               'Level 3 PBLAS timing program.'
         WRITE( NOUT, FMT = 9999 ) USRINFO
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'Tests of the complex double precision '//
     $               'Level 3 PBLAS'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9992 ) NMAT
         WRITE( NOUT, FMT = 9986 ) NBLOG
         WRITE( NOUT, FMT = 9991 ) NGRIDS
         WRITE( NOUT, FMT = 9989 )
     $               'P', ( PVAL(I), I = 1, MIN(NGRIDS, 5) )
         IF( NGRIDS.GT.5 )
     $      WRITE( NOUT, FMT = 9990 ) ( PVAL(I), I = 6,
     $                                  MIN( 10, NGRIDS ) )
         IF( NGRIDS.GT.10 )
     $      WRITE( NOUT, FMT = 9990 ) ( PVAL(I), I = 11,
     $                                  MIN( 15, NGRIDS ) )
         IF( NGRIDS.GT.15 )
     $      WRITE( NOUT, FMT = 9990 ) ( PVAL(I), I = 16, NGRIDS )
         WRITE( NOUT, FMT = 9989 )
     $               'Q', ( QVAL(I), I = 1, MIN(NGRIDS, 5) )
         IF( NGRIDS.GT.5 )
     $      WRITE( NOUT, FMT = 9990 ) ( QVAL(I), I = 6,
     $                                  MIN( 10, NGRIDS ) )
         IF( NGRIDS.GT.10 )
     $      WRITE( NOUT, FMT = 9990 ) ( QVAL(I), I = 11,
     $                                  MIN( 15, NGRIDS ) )
         IF( NGRIDS.GT.15 )
     $      WRITE( NOUT, FMT = 9990 ) ( QVAL(I), I = 16, NGRIDS )
         WRITE( NOUT, FMT = 9994 ) ALPHA
         WRITE( NOUT, FMT = 9993 ) BETA
         IF( LTEST( 1 ) ) THEN
            WRITE( NOUT, FMT = 9988 ) SNAMES( 1 ), ' ... Yes'
         ELSE
            WRITE( NOUT, FMT = 9988 ) SNAMES( 1 ), ' ... No '
         END IF
         DO 90 I = 2, NSUBS
            IF( LTEST( I ) ) THEN
               WRITE( NOUT, FMT = 9987 ) SNAMES( I ), ' ... Yes'
            ELSE
               WRITE( NOUT, FMT = 9987 ) SNAMES( I ), ' ... No '
            END IF
   90    CONTINUE
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
         CALL ZGEBR2D( ICTXT, 'All', ' ', 1, 1, ALPHA, 1, 0, 0 )
         CALL ZGEBR2D( ICTXT, 'All', ' ', 1, 1, BETA,  1, 0, 0 )
*
         CALL IGEBR2D( ICTXT, 'All', ' ', 3, 1, WORK, 3, 0, 0 )
         NGRIDS = WORK( 1 )
         NMAT   = WORK( 2 )
         NBLOG  = WORK( 3 )
*
         I = 2*NGRIDS + 38*NMAT + NSUBS
         CALL IGEBR2D( ICTXT, 'All', ' ', I, 1, WORK, I, 0, 0 )
*
         I = 1
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
 9994 FORMAT( 2X, 'Alpha                     :      (', G16.6,
     $        ',', G16.6, ')' )
 9993 FORMAT( 2X, 'Beta                      :      (', G16.6,
     $        ',', G16.6, ')' )
 9992 FORMAT( 2X, 'Number of Tests           : ', I6 )
 9991 FORMAT( 2X, 'Number of process grids   : ', I6 )
 9990 FORMAT( 2X, '                          : ', 5I6 )
 9989 FORMAT( 2X, A1, '                         : ', 5I6 )
 9988 FORMAT( 2X, 'Routines to be tested     :      ', A, A8 )
 9987 FORMAT( 2X, '                                 ', A, A8 )
 9986 FORMAT( 2X, 'Logical block size        : ', I6 )
*
*     End of PZBLA3TIMINFO
*
      END
