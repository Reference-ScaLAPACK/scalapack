      BLOCK DATA
      INTEGER NSUBS
      PARAMETER (NSUBS = 8)
      CHARACTER*7        SNAMES( 8 )
      COMMON             /SNAMEC/SNAMES
      DATA               SNAMES/'PCGEMV ', 'PCHEMV ', 'PCTRMV ',
     $                   'PCTRSV ', 'PCGERU ', 'PCGERC ',
     $                   'PCHER  ', 'PCHER2 '/
      END BLOCK DATA
      
      PROGRAM PCBLA2TIM
*
*  -- PBLAS timing driver (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*  Purpose
*  =======
*
*  PCBLA2TIM  is the main timing program for the Level 2 PBLAS routines.
*
*  The program must be driven by a short data file.  An  annotated exam-
*  ple of a data file can be obtained by deleting the first 3 characters
*  from the following 56 lines:
*  'Level 2 PBLAS, Timing input file'
*  'Intel iPSC/860 hypercube, gamma model.'
*  'PCBLAS2TIM.SUMM'   output file name (if any)
*  6       device out
*  10              value of the logical computational blocksize NB
*  1               number of process grids (ordered pairs of P & Q)
*  2 2 1 4 2 3 8   values of P
*  2 2 4 1 3 2 1   values of Q
*  (1.0E0, 0.0E0)  value of ALPHA
*  (1.0E0, 0.0E0)  value of BETA
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
*  PCGEMV  T  put F for no test in the same column
*  PCHEMV  T  put F for no test in the same column
*  PCTRMV  T  put F for no test in the same column
*  PCTRSV  T  put F for no test in the same column
*  PCGERU  T  put F for no test in the same column
*  PCGERC  T  put F for no test in the same column
*  PCHER   T  put F for no test in the same column
*  PCHER2  T  put F for no test in the same column
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
*  REALSZ  INTEGER
*  CPLXSZ  INTEGER
*          REALSZ  and  CPLXSZ indicate the length in bytes on the given
*          platform  for a  single precision real and a single precision
*          complex. By default,  REALSZ is set to four and CPLXSZ is set
*          to eight.
*
*  MEM     COMPLEX array
*          MEM is an array of dimension TOTMEM / CPLXSZ.
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
      INTEGER            MAXTESTS, MAXGRIDS, CPLXSZ, TOTMEM, MEMSIZ,
     $                   NSUBS
      COMPLEX            ONE
      PARAMETER          ( MAXTESTS = 20, MAXGRIDS = 20, CPLXSZ = 8,
     $                   ONE = ( 1.0E+0, 0.0E+0 ), TOTMEM = 2000000,
     $                   NSUBS = 8, MEMSIZ = TOTMEM / CPLXSZ )
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      CHARACTER*1        AFORM, DIAG, DIAGDO, TRANS, UPLO
      INTEGER            CSRCA, CSRCX, CSRCY, I, IA, IAM, IASEED, ICTXT,
     $                   IMBA, IMBX, IMBY, IMIDA, IMIDX, IMIDY, INBA,
     $                   INBX, INBY, INCX, INCY, IPA, IPOSTA, IPOSTX,
     $                   IPOSTY, IPREA, IPREX, IPREY, IPX, IPY, IX,
     $                   IXSEED, IY, IYSEED, J, JA, JX, JY, K, M, MA,
     $                   MBA, MBX, MBY, MEMREQD, MPA, MPX, MPY, MX, MY,
     $                   MYCOL, MYROW, N, NA, NBA, NBX, NBY, NCOLA,
     $                   NGRIDS, NLX, NLY, NOUT, NPCOL, NPROCS, NPROW,
     $                   NQA, NQX, NQY, NROWA, NTESTS, NX, NY, OFFD,
     $                   RSRCA, RSRCX, RSRCY
      DOUBLE PRECISION   CFLOPS, NOPS, WFLOPS
      COMPLEX            ALPHA, BETA, SCALE
*     ..
*     .. Local Arrays ..
      LOGICAL            LTEST( NSUBS ), YCHECK( NSUBS )
      CHARACTER*1        DIAGVAL( MAXTESTS ), TRANVAL( MAXTESTS ),
     $                   UPLOVAL( MAXTESTS )
      CHARACTER*80       OUTFILE
      INTEGER            CSCAVAL( MAXTESTS ), CSCXVAL( MAXTESTS ),
     $                   CSCYVAL( MAXTESTS ), DESCA( DLEN_ ),
     $                   DESCX( DLEN_ ), DESCY( DLEN_ ),
     $                   IAVAL( MAXTESTS ), IERR( 3 ),
     $                   IMBAVAL( MAXTESTS ), IMBXVAL( MAXTESTS ),
     $                   IMBYVAL( MAXTESTS ), INBAVAL( MAXTESTS ),
     $                   INBXVAL( MAXTESTS ), INBYVAL( MAXTESTS ),
     $                   INCXVAL( MAXTESTS ), INCYVAL( MAXTESTS ),
     $                   IXVAL( MAXTESTS ), IYVAL( MAXTESTS ),
     $                   JAVAL( MAXTESTS ), JXVAL( MAXTESTS ),
     $                   JYVAL( MAXTESTS ), MAVAL( MAXTESTS ),
     $                   MBAVAL( MAXTESTS ), MBXVAL( MAXTESTS ),
     $                   MBYVAL( MAXTESTS ), MVAL( MAXTESTS ),
     $                   MXVAL( MAXTESTS ), MYVAL( MAXTESTS ),
     $                   NAVAL( MAXTESTS ), NBAVAL( MAXTESTS ),
     $                   NBXVAL( MAXTESTS ), NBYVAL( MAXTESTS ),
     $                   NVAL( MAXTESTS ), NXVAL( MAXTESTS ),
     $                   NYVAL( MAXTESTS ), PVAL( MAXTESTS ),
     $                   QVAL( MAXTESTS ), RSCAVAL( MAXTESTS ),
     $                   RSCXVAL( MAXTESTS ), RSCYVAL( MAXTESTS )
      DOUBLE PRECISION   CTIME( 1 ), WTIME( 1 )
      COMPLEX            MEM( MEMSIZ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_EXIT, BLACS_GET,
     $                   BLACS_GRIDEXIT, BLACS_GRIDINFO, BLACS_GRIDINIT,
     $                   BLACS_PINFO, IGSUM2D, PB_BOOT, PB_COMBINE,
     $                   PB_TIMER, PCBLA2TIMINFO, PCGEMV, PCGERC,
     $                   PCGERU, PCHEMV, PCHER, PCHER2, PCLAGEN,
     $                   PCLASCAL, PCTRMV, PCTRSV, PMDESCCHK, PMDIMCHK,
     $                   PVDESCCHK, PVDIMCHK
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PDOPBL2
      EXTERNAL           LSAME, PDOPBL2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, DBLE, MAX, REAL
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
     $                   .TRUE., .TRUE., .FALSE., .TRUE./
*     ..
*     .. Executable Statements ..
*
*     Initialization
*
*     Set flag so that the PBLAS error handler won't abort on errors, so
*     that the tester will detect unsupported operations.
*
      ABRTFLG = .TRUE.
*
*     Seeds for random matrix generations.
*
      IASEED = 100
      IXSEED = 200
      IYSEED = 300
*
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
      CALL PCBLA2TIMINFO( OUTFILE, NOUT, NTESTS, DIAGVAL, TRANVAL,
     $                    UPLOVAL, MVAL, NVAL, MAVAL, NAVAL, IMBAVAL,
     $                    MBAVAL, INBAVAL, NBAVAL, RSCAVAL, CSCAVAL,
     $                    IAVAL, JAVAL, MXVAL, NXVAL, IMBXVAL, MBXVAL,
     $                    INBXVAL, NBXVAL, RSCXVAL, CSCXVAL, IXVAL,
     $                    JXVAL, INCXVAL, MYVAL, NYVAL, IMBYVAL,
     $                    MBYVAL, INBYVAL, NBYVAL, RSCYVAL,
     $                    CSCYVAL, IYVAL, JYVAL, INCYVAL, MAXTESTS,
     $                    NGRIDS, PVAL, MAXGRIDS, QVAL, MAXGRIDS,
     $                    NBLOG, LTEST, IAM, NPROCS, ALPHA, BETA, MEM )
*
      IF( IAM.EQ.0 )
     $   WRITE( NOUT, FMT = 9983 )
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
            MBA   = MBAVAL( J )
            INBA  = INBAVAL( J )
            NBA   = NBAVAL( J )
            RSRCA = RSCAVAL( J )
            CSRCA = CSCAVAL( J )
            IA    = IAVAL( J )
            JA    = JAVAL( J )
*
            MX    = MXVAL( J )
            NX    = NXVAL( J )
            IMBX  = IMBXVAL( J )
            MBX   = MBXVAL( J )
            INBX  = INBXVAL( J )
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
            MBY   = MBYVAL( J )
            INBY  = INBYVAL( J )
            NBY   = NBYVAL( J )
            RSRCY = RSCYVAL( J )
            CSRCY = CSCYVAL( J )
            IY    = IYVAL( J )
            JY    = JYVAL( J )
            INCY  = INCYVAL( J )
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
               WRITE( NOUT, FMT = 9980 )
*
            END IF
*
*           Check the validity of the input test parameters
*
            IF( .NOT.LSAME( UPLO, 'U' ).AND.
     $          .NOT.LSAME( UPLO, 'L' ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'UPLO'
               GO TO 40
            END IF
*
            IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $          .NOT.LSAME( TRANS, 'T' ).AND.
     $          .NOT.LSAME( TRANS, 'C' ) ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) 'TRANS'
               GO TO 40
            END IF
*
            IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' ) )THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9997 ) TRANS
               WRITE( NOUT, FMT = 9997 ) 'DIAG'
               GO TO 40
            END IF
*
*           Check and initialize the matrix descriptors
*
            CALL PMDESCCHK( ICTXT, NOUT, 'A', DESCA,
     $                      BLOCK_CYCLIC_2D_INB, MA, NA, IMBA, INBA,
     $                      MBA, NBA, RSRCA, CSRCA, MPA, NQA, IPREA,
     $                      IMIDA, IPOSTA, 0, 0, IERR( 1 ) )
            CALL PVDESCCHK( ICTXT, NOUT, 'X', DESCX,
     $                      BLOCK_CYCLIC_2D_INB, MX, NX, IMBX, INBX,
     $                      MBX, NBX, RSRCX, CSRCX, INCX, MPX, NQX,
     $                      IPREX, IMIDX, IPOSTX, 0, 0, IERR( 2 ) )
            CALL PVDESCCHK( ICTXT, NOUT, 'Y', DESCY,
     $                      BLOCK_CYCLIC_2D_INB, MY, NY, IMBY, INBY,
     $                      MBY, NBY, RSRCY, CSRCY, INCY, MPY, NQY,
     $                      IPREY, IMIDY, IPOSTY, 0, 0, IERR( 3 ) )
*
            IF( IERR( 1 ).GT.0 .OR. IERR( 2 ).GT.0 .OR.
     $          IERR( 3 ).GT.0 ) THEN
               GO TO 40
            END IF
*
*           Assign pointers into MEM for matrices corresponding to
*           the distributed matrices A, X and Y.
*
            IPA = 1
            IPX = IPA + DESCA( LLD_ ) * NQA
            IPY = IPX + DESCX( LLD_ ) * NQX
*
*           Check if sufficient memory.
*
            MEMREQD = IPY + DESCY( LLD_ ) * NQY - 1
            IERR( 1 ) = 0
            IF( MEMREQD.GT.MEMSIZ ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9986 ) MEMREQD*CPLXSZ
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
               ELSE IF( K.EQ.5 .OR. K.EQ.6 ) THEN
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
                  GO TO 30
               END IF
*
*              Generate distributed matrices A, X and Y
*
               IF( K.EQ.2 .OR. K.EQ.7 .OR. K.EQ.8 ) THEN
                  AFORM  = 'H'
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
               CALL PCLAGEN( .FALSE., AFORM, DIAGDO, OFFD, MA, NA,
     $                       1, 1, DESCA, IASEED, MEM( IPA ),
     $                       DESCA( LLD_ ) )
               CALL PCLAGEN( .FALSE., 'None', 'No diag', 0, MX, NX,
     $                       1, 1, DESCX, IXSEED, MEM( IPX ),
     $                       DESCX( LLD_ ) )
               IF( YCHECK( K ) )
     $            CALL PCLAGEN( .FALSE., 'None', 'No diag', 0, MY,
     $                          NY, 1, 1, DESCY, IYSEED, MEM( IPY ),
     $                          DESCY( LLD_ ) )
*
               IF( ( K.EQ.4 ).AND.( .NOT.( LSAME( DIAG, 'N' ) ) ).AND.
     $             ( MAX( NROWA, NCOLA ).GT.1 ) ) THEN
                  SCALE = ONE / CMPLX( REAL( MAX( NROWA, NCOLA ) ) )
                  IF( LSAME( UPLO, 'L' ) ) THEN
                     CALL PCLASCAL( 'Lower', NROWA-1, NCOLA-1, SCALE,
     $                              MEM( IPA ), IA+1, JA, DESCA )
                  ELSE
                     CALL PCLASCAL( 'Upper', NROWA-1, NCOLA-1, SCALE,
     $                              MEM( IPA ), IA, JA+1, DESCA )
                  END IF
               END IF
*
               INFO = 0
               CALL PB_BOOT()
               CALL BLACS_BARRIER( ICTXT, 'All' )
*
*              Call the Level 2 PBLAS routine
*
               IF( K.EQ.1 ) THEN
*
*                 Test PCGEMV
*
                  CALL PB_TIMER( 1 )
                  CALL PCGEMV( TRANS, M, N, ALPHA, MEM( IPA ), IA, JA,
     $                         DESCA, MEM( IPX ), IX, JX, DESCX, INCX,
     $                         BETA, MEM( IPY ), IY, JY, DESCY, INCY )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( K.EQ.2 ) THEN
*
*                 Test PCHEMV
*
                  CALL PB_TIMER( 1 )
                  CALL PCHEMV( UPLO, N, ALPHA, MEM( IPA ), IA, JA,
     $                         DESCA, MEM( IPX ), IX, JX, DESCX, INCX,
     $                         BETA, MEM( IPY ), IY, JY, DESCY, INCY )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( K.EQ.3 ) THEN
*
*                 Test PCTRMV
*
                  CALL PB_TIMER( 1 )
                  CALL PCTRMV( UPLO, TRANS, DIAG, N, MEM( IPA ), IA, JA,
     $                         DESCA, MEM( IPX ), IX, JX, DESCX, INCX )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( K.EQ.4 ) THEN
*
*                 Test PCTRSV
*
                  CALL PB_TIMER( 1 )
                  CALL PCTRSV( UPLO, TRANS, DIAG, N, MEM( IPA ), IA, JA,
     $                         DESCA, MEM( IPX ), IX, JX, DESCX, INCX )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( K.EQ.5 ) THEN
*
*                 Test PCGERU
*
                  CALL PB_TIMER( 1 )
                  CALL PCGERU( M, N, ALPHA, MEM( IPX ), IX, JX, DESCX,
     $                         INCX, MEM( IPY ), IY, JY, DESCY, INCY,
     $                         MEM( IPA ), IA, JA, DESCA )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( K.EQ.6 ) THEN
*
*                 Test PCGERC
*
                  CALL PB_TIMER( 1 )
                  CALL PCGERC( M, N, ALPHA, MEM( IPX ), IX, JX, DESCX,
     $                         INCX, MEM( IPY ), IY, JY, DESCY, INCY,
     $                         MEM( IPA ), IA, JA, DESCA )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( K.EQ.7 ) THEN
*
*                 Test PCHER
*
                  CALL PB_TIMER( 1 )
                  CALL PCHER( UPLO, N, REAL( ALPHA ), MEM( IPX ), IX,
     $                        JX, DESCX, INCX, MEM( IPA ), IA, JA,
     $                        DESCA )
                  CALL PB_TIMER( 1 )
*
               ELSE IF( K.EQ.8 ) THEN
*
*                 Test PCHER2
*
                  CALL PB_TIMER( 1 )
                  CALL PCHER2( UPLO, N, ALPHA, MEM( IPX ), IX, JX,
     $                         DESCX, INCX, MEM( IPY ), IY, JY, DESCY,
     $                         INCY, MEM( IPA ), IA, JA, DESCA )
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
*                 Calculate total flops
*
                  NOPS = PDOPBL2( SNAMES( K ), NROWA, NCOLA, 0, 0 )
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
                  WRITE( NOUT, FMT = 9981 ) SNAMES( K ), WTIME( 1 ),
     $                                      WFLOPS, CTIME( 1 ), CFLOPS
*
               END IF
*
   30       CONTINUE
*
   40       IF( IAM.EQ.0 ) THEN
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = * )
               WRITE( NOUT, FMT = 9985 ) J
            END IF
*
   50   CONTINUE
*
        CALL BLACS_GRIDEXIT( ICTXT )
*
   60 CONTINUE
*
*     Print results
*
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9984 )
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
 9985 FORMAT( 2X, 'Test number ', I2, ' completed.' )
 9984 FORMAT( 2X, 'End of Tests.' )
 9983 FORMAT( 2X, 'Tests started.' )
 9982 FORMAT( 2X, '   ***** Operation not supported, error code: ',
     $        I5, ' *****' )
 9981 FORMAT( 2X, '|  ', A, 2X, F13.3, 2X, F13.3, 2X, F13.3, 2X, F13.3 )
 9980 FORMAT( 2X, '            WALL time (s)    WALL Mflops ',
     $        '  CPU time (s)     CPU Mflops' )
*
      STOP
*
*     End of PCBLA2TIM
*
      END
      SUBROUTINE PCBLA2TIMINFO( SUMMRY, NOUT, NMAT, DIAGVAL, TRANVAL,
     $                          UPLOVAL, MVAL, NVAL, MAVAL, NAVAL,
     $                          IMBAVAL, MBAVAL, INBAVAL, NBAVAL,
     $                          RSCAVAL, CSCAVAL, IAVAL, JAVAL,
     $                          MXVAL, NXVAL, IMBXVAL, MBXVAL,
     $                          INBXVAL, NBXVAL, RSCXVAL, CSCXVAL,
     $                          IXVAL, JXVAL, INCXVAL, MYVAL, NYVAL,
     $                          IMBYVAL, MBYVAL, INBYVAL, NBYVAL,
     $                          RSCYVAL, CSCYVAL, IYVAL, JYVAL,
     $                          INCYVAL, LDVAL, NGRIDS, PVAL, LDPVAL,
     $                          QVAL, LDQVAL, NBLOG, LTEST, IAM, NPROCS,
     $                          ALPHA, BETA, WORK )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            IAM, LDPVAL, LDQVAL, LDVAL, NBLOG, NGRIDS,
     $                   NMAT, NOUT, NPROCS
      COMPLEX            ALPHA, BETA
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
*  PCBLA2TIMINFO  get  the needed startup information for timing various
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
*          On entry, LTEST  is an array of dimension at least eight.  On
*          exit, if LTEST( i ) is .TRUE., the i-th Level 2 PBLAS routine
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
*  ALPHA   (global output) COMPLEX
*          On exit, ALPHA specifies the value of alpha to be used in all
*          the test cases.
*
*  BETA    (global output) COMPLEX
*          On exit, BETA  specifies the value of beta  to be used in all
*          the test cases.
*
*  WORK    (local workspace) INTEGER array
*          On   entry,   WORK   is   an  array  of  dimension  at  least
*          MAX( 3, 2*NGRIDS+37*NMAT+NSUBS ) with NSUBS = 8.  This  array
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
      PARAMETER          ( NIN = 11, NSUBS = 8 )
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
     $                   BLACS_GRIDINIT, BLACS_SETUP, CGEBR2D, CGEBS2D,
     $                   ICOPY, IGEBR2D, IGEBS2D, SGEBR2D, SGEBS2D
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
         OPEN( NIN, FILE='PCBLAS2TIM.dat', STATUS='OLD' )
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
         END IF
*
*        Read in input data into arrays.
*
         READ( NIN, FMT = * ) ( UPLOVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( TRANVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( DIAGVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( MVAL   ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( NVAL   ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( MAVAL  ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( NAVAL  ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( IMBAVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( INBAVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( MBAVAL ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( NBAVAL ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( RSCAVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( CSCAVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( IAVAL  ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( JAVAL  ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( MXVAL  ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( NXVAL  ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( IMBXVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( INBXVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( MBXVAL ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( NBXVAL ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( RSCXVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( CSCXVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( IXVAL  ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( JXVAL  ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( INCXVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( MYVAL  ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( NYVAL  ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( IMBYVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( INBYVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( MBYVAL ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( NBYVAL ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( RSCYVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( CSCYVAL( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( IYVAL  ( I ),  I = 1, NMAT )
         READ( NIN, FMT = * ) ( JYVAL  ( I ),  I = 1, NMAT )
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
*        Pack information arrays and broadcast
*
         CALL CGEBS2D( ICTXT, 'All', ' ', 1, 1, ALPHA, 1 )
         CALL CGEBS2D( ICTXT, 'All', ' ', 1, 1, BETA, 1 )
*
         WORK( 1 ) = NGRIDS
         WORK( 2 ) = NMAT
         WORK( 3 ) = NBLOG
         CALL IGEBS2D( ICTXT, 'All', ' ', 3, 1, WORK, 3 )
*
         I = 1
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
         WRITE( NOUT, FMT = 9999 )
     $               'Level 2 PBLAS timing program.'
         WRITE( NOUT, FMT = 9999 ) USRINFO
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'Tests of the complex single precision '//
     $               'Level 2 PBLAS'
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
         DO 90 I = 1, NSUBS
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
         CALL CGEBR2D( ICTXT, 'All', ' ', 1, 1, ALPHA, 1, 0, 0 )
         CALL CGEBR2D( ICTXT, 'All', ' ', 1, 1, BETA, 1, 0, 0 )
*
         CALL IGEBR2D( ICTXT, 'All', ' ', 3, 1, WORK, 3, 0, 0 )
         NGRIDS = WORK( 1 )
         NMAT   = WORK( 2 )
         NBLOG  = WORK( 3 )
*
         I = 2*NGRIDS + 37*NMAT + NSUBS
         CALL IGEBR2D( ICTXT, 'All', ' ', I, 1, WORK, I, 0, 0 )
*
         I = 1
         DO 100 J = 1, NMAT
            DIAGVAL( J ) = CHAR( WORK( I ) )
            TRANVAL( J ) = CHAR( WORK( I+1 ) )
            UPLOVAL( J ) = CHAR( WORK( I+2 ) )
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
*     End of PCBLA2TIMINFO
*
      END
