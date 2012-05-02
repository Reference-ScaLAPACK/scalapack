      BLOCK DATA
      INTEGER NSUBS
      PARAMETER (NSUBS = 8)
      CHARACTER*7        SNAMES( NSUBS )
      COMMON             /SNAMEC/SNAMES
      DATA               SNAMES/'PDSWAP ', 'PDSCAL ', 'PDCOPY ',
     $                   'PDAXPY ', 'PDDOT  ', 'PDNRM2 ',
     $                   'PDASUM ', 'PDAMAX '/
      END BLOCK DATA

      PROGRAM PDBLA1TST
*
*  -- PBLAS testing driver (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*  Purpose
*  =======
*
*  PDBLA1TST is the main testing program for the PBLAS Level 1 routines.
*
*  The program must be driven by a short data file.  An  annotated exam-
*  ple of a data file can be obtained by deleting the first 3 characters
*  from the following 44 lines:
*  'Level 1 PBLAS, Testing input file'
*  'Intel iPSC/860 hypercube, gamma model.'
*  'PDBLAS1TST.SUMM'            output file name (if any)
*  6       device out
*  F       logical flag, T to stop on failures
*  F       logical flag, T to test error exits
*  0       verbosity, 0 for pass/fail, 1-3 for matrix dump on errors
*  10      the leading dimension gap
*  1       number of process grids (ordered pairs of P & Q)
*  2 2 1 4 2 3 8        values of P
*  2 2 4 1 3 2 1        values of Q
*  1.0D0                value of ALPHA
*  2                    number of tests problems
*  3  4                 values of N
*  6 10                 values of M_X
*  6 10                 values of N_X
*  2  5                 values of IMB_X
*  2  5                 values of INB_X
*  2  5                 values of MB_X
*  2  5                 values of NB_X
*  0  1                 values of RSRC_X
*  0  0                 values of CSRC_X
*  1  1                 values of IX
*  1  1                 values of JX
*  1  1                 values of INCX
*  6 10                 values of M_Y
*  6 10                 values of N_Y
*  2  5                 values of IMB_Y
*  2  5                 values of INB_Y
*  2  5                 values of MB_Y
*  2  5                 values of NB_Y
*  0  1                 values of RSRC_Y
*  0  0                 values of CSRC_Y
*  1  1                 values of IY
*  1  1                 values of JY
*  6  1                 values of INCY
*  PDSWAP  T            put F for no test in the same column
*  PDSCAL  T            put F for no test in the same column
*  PDCOPY  T            put F for no test in the same column
*  PDAXPY  T            put F for no test in the same column
*  PDDOT   T            put F for no test in the same column
*  PDNRM2  T            put F for no test in the same column
*  PDASUM  T            put F for no test in the same column
*  PDAMAX  T            put F for no test in the same column
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
      DOUBLE PRECISION   PADVAL, ZERO
      PARAMETER          ( MAXTESTS = 20, MAXGRIDS = 20, GAPMUL = 10,
     $                   DBLESZ = 8, TOTMEM = 2000000,
     $                   MEMSIZ = TOTMEM / DBLESZ, ZERO = 0.0D+0,
     $                   PADVAL = -9923.0D+0, NSUBS = 8 )
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
      INTEGER            CSRCX, CSRCY, I, IAM, ICTXT, IGAP, IMBX, IMBY,
     $                   IMIDX, IMIDY, INBX, INBY, INCX, INCY, IPMATX,
     $                   IPMATY, IPOSTX, IPOSTY, IPREX, IPREY, IPW, IPX,
     $                   IPY, IVERB, IX, IXSEED, IY, IYSEED, J, JX, JY,
     $                   K, LDX, LDY, MBX, MBY, MEMREQD, MPX, MPY, MX,
     $                   MY, MYCOL, MYROW, N, NBX, NBY, NGRIDS, NOUT,
     $                   NPCOL, NPROCS, NPROW, NQX, NQY, NTESTS, NX, NY,
     $                   PISCLR, RSRCX, RSRCY, TSKIP, TSTCNT
      DOUBLE PRECISION   ALPHA, PSCLR, PUSCLR
*     ..
*     .. Local Arrays ..
      CHARACTER*80       OUTFILE
      LOGICAL            LTEST( NSUBS ), YCHECK( NSUBS )
      INTEGER            CSCXVAL( MAXTESTS ), CSCYVAL( MAXTESTS ),
     $                   DESCX( DLEN_ ), DESCXR( DLEN_ ),
     $                   DESCY( DLEN_ ), DESCYR( DLEN_ ), IERR( 4 ),
     $                   IMBXVAL( MAXTESTS ), IMBYVAL( MAXTESTS ),
     $                   INBXVAL( MAXTESTS ), INBYVAL( MAXTESTS ),
     $                   INCXVAL( MAXTESTS ), INCYVAL( MAXTESTS ),
     $                   IXVAL( MAXTESTS ), IYVAL( MAXTESTS ),
     $                   JXVAL( MAXTESTS ), JYVAL( MAXTESTS ),
     $                   KFAIL( NSUBS ), KPASS( NSUBS ), KSKIP( NSUBS ),
     $                   KTESTS( NSUBS ), MBXVAL( MAXTESTS ),
     $                   MBYVAL( MAXTESTS ), MXVAL( MAXTESTS ),
     $                   MYVAL( MAXTESTS ), NBXVAL( MAXTESTS ),
     $                   NBYVAL( MAXTESTS ), NVAL( MAXTESTS ),
     $                   NXVAL( MAXTESTS ), NYVAL( MAXTESTS ),
     $                   PVAL( MAXTESTS ), QVAL( MAXTESTS ),
     $                   RSCXVAL( MAXTESTS ), RSCYVAL( MAXTESTS )
      DOUBLE PRECISION   MEM( MEMSIZ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $                   IGSUM2D, PB_DCHEKPAD, PB_DESCSET2, PB_DFILLPAD,
     $                   PB_PDLAPRNT, PDAMAX, PDASUM, PDAXPY,
     $                   PDBLA1TSTINFO, PDBLAS1TSTCHK, PDBLAS1TSTCHKE,
     $                   PDCHKARG1, PDCHKVOUT, PDCOPY, PDDOT, PDLAGEN,
     $                   PDMPRNT, PDNRM2, PDSCAL, PDSWAP, PDVPRNT,
     $                   PVDESCCHK, PVDIMCHK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MOD
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
      DATA               YCHECK/.TRUE., .FALSE., .TRUE., .TRUE., .TRUE.,
     $                   .FALSE., .FALSE., .FALSE./
*     ..
*     .. Executable Statements ..
*
*     Initialization
*
*     Set flag so that the PBLAS error handler will abort on errors.
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
      IXSEED = 100
      IYSEED = 200
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
      CALL PDBLA1TSTINFO( OUTFILE, NOUT, NTESTS, NVAL, MXVAL, NXVAL,
     $                    IMBXVAL, MBXVAL, INBXVAL, NBXVAL, RSCXVAL,
     $                    CSCXVAL, IXVAL, JXVAL, INCXVAL, MYVAL,
     $                    NYVAL, IMBYVAL, MBYVAL, INBYVAL, NBYVAL,
     $                    RSCYVAL, CSCYVAL, IYVAL, JYVAL, INCYVAL,
     $                    MAXTESTS, NGRIDS, PVAL, MAXGRIDS, QVAL,
     $                    MAXGRIDS, LTEST, SOF, TEE, IAM, IGAP, IVERB,
     $                    NPROCS, ALPHA, MEM )
*
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9979 )
         WRITE( NOUT, FMT = * )
      END IF
*
*     If TEE is set then Test Error Exits of routines.
*
      IF( TEE )
     $   CALL PDBLAS1TSTCHKE( LTEST, NOUT, NPROCS )
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
            N     = NVAL( J )
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
               TSTCNT = TSTCNT + 1
               WRITE( NOUT, FMT = * )
               WRITE( NOUT, FMT = 9996 ) TSTCNT, NPROW, NPCOL
               WRITE( NOUT, FMT = * )
*
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9994 )
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9993 ) N, IX, JX, MX, NX, IMBX, INBX,
     $                                   MBX, NBX, RSRCX, CSRCX, INCX
*
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9992 )
               WRITE( NOUT, FMT = 9995 )
               WRITE( NOUT, FMT = 9993 ) N, IY, JY, MY, NY, IMBY, INBY,
     $                                   MBY, NBY, RSRCY, CSRCY, INCY
               WRITE( NOUT, FMT = 9995 )
            END IF
*
*           Check the validity of the input and initialize DESC_
*
            CALL PVDESCCHK( ICTXT, NOUT, 'X', DESCX,
     $                      BLOCK_CYCLIC_2D_INB, MX, NX, IMBX, INBX,
     $                      MBX, NBX, RSRCX, CSRCX, INCX, MPX, NQX,
     $                      IPREX, IMIDX, IPOSTX, IGAP, GAPMUL,
     $                      IERR( 1 ) )
            CALL PVDESCCHK( ICTXT, NOUT, 'Y', DESCY,
     $                      BLOCK_CYCLIC_2D_INB, MY, NY, IMBY, INBY,
     $                      MBY, NBY, RSRCY, CSRCY, INCY, MPY, NQY,
     $                      IPREY, IMIDY, IPOSTY, IGAP, GAPMUL,
     $                      IERR( 2 ) )
*
            IF( IERR( 1 ).GT.0 .OR. IERR( 2 ).GT.0 ) THEN
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
            LDX = MAX( 1, MX )
            LDY = MAX( 1, MY )
*
*           Assign pointers into MEM for matrices corresponding to
*           vectors X and Y. Ex: IPX starts at position MEM( IPREX+1 ).
*
            IPX    = IPREX + 1
            IPY    = IPX + DESCX( LLD_ ) * NQX + IPOSTX + IPREY
            IPMATX = IPY + DESCY( LLD_ ) * NQY + IPOSTY
            IPMATY = IPMATX + MX * NX
            IPW    = IPMATY + MY * NY
*
*           Check if sufficient memory.
*           Requirement = mem for local part of parallel matrices +
*                         mem for whole matrices for comp. check +
*                         mem for recving comp. check error vals.
*
            MEMREQD = IPW - 1 +
     $                MAX( MAX( IMBX, MBX ), MAX( IMBY, MBY ) )
            IERR( 1 ) = 0
            IF( MEMREQD.GT.MEMSIZ ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9990 ) MEMREQD*DBLESZ
               IERR( 1 ) = 1
            END IF
*
*           Check all processes for an error
*
            CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, IERR, 1, -1, 0 )
*
            IF( IERR( 1 ).GT.0 ) THEN
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9991 )
               TSKIP = TSKIP + 1
               GO TO 40
            END IF
*
*           Loop over all PBLAS 1 routines
*
            DO 30 K = 1, NSUBS
*
*              Continue only if this sub has to be tested.
*
               IF( .NOT.LTEST( K ) )
     $            GO TO 30
*
               IF( IAM.EQ.0 ) THEN
                  WRITE( NOUT, FMT = * )
                  WRITE( NOUT, FMT = 9989 ) SNAMES( K )
               END IF
*
*              Check the validity of the operand sizes
*
               CALL PVDIMCHK( ICTXT, NOUT, N, 'X', IX, JX, DESCX, INCX,
     $                        IERR( 1 ) )
               CALL PVDIMCHK( ICTXT, NOUT, N, 'Y', IY, JY, DESCY, INCY,
     $                        IERR( 2 ) )
*
               IF( IERR( 1 ).NE.0 .OR. IERR( 2 ).NE.0 ) THEN
                  KSKIP( K ) = KSKIP( K ) + 1
                  GO TO 30
               END IF
*
*              Generate distributed matrices X and Y
*
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
               CALL PB_DESCSET2( DESCXR, MX, NX, IMBX, INBX, MBX, NBX,
     $                           -1, -1, ICTXT, MAX( 1, MX ) )
               CALL PDLAGEN( .FALSE., 'None', 'No diag', 0, MX, NX, 1,
     $                       1, DESCXR, IXSEED, MEM( IPMATX ),
     $                       DESCXR( LLD_ ) )
               IF( YCHECK( K ) ) THEN
                  CALL PB_DESCSET2( DESCYR, MY, NY, IMBY, INBY, MBY,
     $                              NBY, -1, -1, ICTXT, MAX( 1, MY ) )
                  CALL PDLAGEN( .FALSE., 'None', 'No diag', 0, MY, NY,
     $                          1, 1, DESCYR, IYSEED, MEM( IPMATY ),
     $                          DESCYR( LLD_ ) )
               END IF
*
*              Pad the guard zones of X, and Y
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
*              Initialize the check for INPUT only args.
*
               INFO = 0
               CALL PDCHKARG1( ICTXT, NOUT, SNAMES( K ), N, ALPHA, IX,
     $                         JX, DESCX, INCX, IY, JY, DESCY, INCY,
     $                         INFO )
*
               INFO = 0
               PSCLR  = ZERO
               PUSCLR = ZERO
               PISCLR = 0
*
*              Print initial parallel data if IVERB >= 2.
*
               IF( IVERB.EQ.2 ) THEN
                  IF( INCX.EQ.DESCX( M_ ) ) THEN
                     CALL PB_PDLAPRNT( 1, N, MEM( IPX ), IX, JX, DESCX,
     $                                 0, 0, 'PARALLEL_INITIAL_X', NOUT,
     $                                 MEM( IPW ) )
                  ELSE
                     CALL PB_PDLAPRNT( N, 1, MEM( IPX ), IX, JX, DESCX,
     $                                 0, 0, 'PARALLEL_INITIAL_X', NOUT,
     $                                 MEM( IPW ) )
                  END IF
                  IF( YCHECK( K ) ) THEN
                     IF( INCY.EQ.DESCY( M_ ) ) THEN
                        CALL PB_PDLAPRNT( 1, N, MEM( IPY ), IY, JY,
     $                                    DESCY, 0, 0,
     $                                    'PARALLEL_INITIAL_Y', NOUT,
     $                                    MEM( IPW ) )
                     ELSE
                        CALL PB_PDLAPRNT( N, 1, MEM( IPY ), IY, JY,
     $                                    DESCY, 0, 0,
     $                                    'PARALLEL_INITIAL_Y', NOUT,
     $                                    MEM( IPW ) )
                     END IF
                  END IF
               ELSE IF( IVERB.GE.3 ) THEN
                  CALL PB_PDLAPRNT( MX, NX, MEM( IPX ), 1, 1, DESCX, 0,
     $                              0, 'PARALLEL_INITIAL_X', NOUT,
     $                              MEM( IPW ) )
                  IF( YCHECK( K ) )
     $               CALL PB_PDLAPRNT( MY, NY, MEM( IPY ), 1, 1, DESCY,
     $                                 0, 0, 'PARALLEL_INITIAL_Y', NOUT,
     $                                 MEM( IPW ) )
               END IF
*
*              Call the PBLAS routine
*
               IF( K.EQ.1 ) THEN
*
*                 Test PDSWAP
*
                  CALL PDSWAP( N, MEM( IPX ), IX, JX, DESCX, INCX,
     $                         MEM( IPY ), IY, JY, DESCY, INCY )
*
               ELSE IF( K.EQ.2 ) THEN
*
*                 Test PDSCAL
*
                  PSCLR = ALPHA
                  CALL PDSCAL( N, ALPHA, MEM( IPX ), IX, JX, DESCX,
     $                         INCX )
*
               ELSE IF( K.EQ.3 ) THEN
*
*                 Test PDCOPY
*
                  CALL PDCOPY( N, MEM( IPX ), IX, JX, DESCX, INCX,
     $                         MEM( IPY ), IY, JY, DESCY, INCY )
*
               ELSE IF( K.EQ.4 ) THEN
*
*                 Test PDAXPY
*
                  PSCLR = ALPHA
                  CALL PDAXPY( N, ALPHA, MEM( IPX ), IX, JX, DESCX,
     $                         INCX, MEM( IPY ), IY, JY, DESCY, INCY )
*
               ELSE IF( K.EQ.5 ) THEN
*
*                 Test PDDOT
*
                  CALL PDDOT( N, PSCLR, MEM( IPX ), IX, JX, DESCX, INCX,
     $                        MEM( IPY ), IY, JY, DESCY, INCY )
*
               ELSE IF( K.EQ.6 ) THEN
*
*                 Test PDNRM2
*
                  CALL PDNRM2( N, PUSCLR, MEM( IPX ), IX, JX, DESCX,
     $                         INCX )
*
               ELSE IF( K.EQ.7 ) THEN
*
*                 Test PDASUM
*
                  CALL PDASUM( N, PUSCLR, MEM( IPX ), IX, JX, DESCX,
     $                         INCX )
*
               ELSE IF( K.EQ.8 ) THEN
*
                  CALL PDAMAX( N, PSCLR, PISCLR, MEM( IPX ), IX, JX,
     $                         DESCX, INCX )
*
               END IF
*
*              Check if the operation has been performed.
*
               IF( INFO.NE.0 ) THEN
                  KSKIP( K ) = KSKIP( K ) + 1
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9978 ) INFO
                  GO TO 30
               END IF
*
*              Check the computations
*
               CALL PDBLAS1TSTCHK( ICTXT, NOUT, K, N, PSCLR, PUSCLR,
     $                             PISCLR, MEM( IPMATX ), MEM( IPX ),
     $                             IX, JX, DESCX, INCX, MEM( IPMATY ),
     $                             MEM( IPY ), IY, JY, DESCY, INCY,
     $                             INFO )
               IF( MOD( INFO, 2 ).EQ.1 ) THEN
                  IERR( 1 ) = 1
               ELSE IF( MOD( INFO / 2, 2 ).EQ.1 ) THEN
                  IERR( 2 ) = 1
               ELSE IF( INFO.NE.0 ) THEN
                  IERR( 1 ) = 1
                  IERR( 2 ) = 1
               END IF
*
*              Check padding
*
               CALL PB_DCHEKPAD( ICTXT, SNAMES( K ), MPX, NQX,
     $                           MEM( IPX-IPREX ), DESCX( LLD_ ),
     $                           IPREX, IPOSTX, PADVAL )
               IF( YCHECK( K ) ) THEN
                  CALL PB_DCHEKPAD( ICTXT, SNAMES( K ), MPY, NQY,
     $                              MEM( IPY-IPREY ), DESCY( LLD_ ),
     $                              IPREY, IPOSTY, PADVAL )
               END IF
*
*              Check input-only scalar arguments
*
               INFO = 1
               CALL PDCHKARG1( ICTXT, NOUT, SNAMES( K ), N, ALPHA, IX,
     $                         JX, DESCX, INCX, IY, JY, DESCY, INCY,
     $                         INFO )
*
*              Check input-only array arguments
*
               CALL PDCHKVOUT( N, MEM( IPMATX ), MEM( IPX ), IX, JX,
     $                         DESCX, INCX, IERR( 3 ) )
*
               IF( IERR( 3 ).NE.0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9986 ) 'PARALLEL_X', SNAMES( K )
               END IF
*
               IF( YCHECK( K ) ) THEN
                  CALL PDCHKVOUT( N, MEM( IPMATY ), MEM( IPY ), IY, JY,
     $                            DESCY, INCY, IERR( 4 ) )
                  IF( IERR( 4 ).NE.0 ) THEN
                     IF( IAM.EQ.0 )
     $                  WRITE( NOUT, FMT = 9986 ) 'PARALLEL_Y',
     $                                       SNAMES( K )
                  END IF
               END IF
*
*              Only node 0 prints computational test result
*
               IF( INFO.NE.0 .OR. IERR( 1 ).NE.0 .OR.
     $             IERR( 2 ).NE.0 .OR. IERR( 3 ).NE.0 .OR.
     $             IERR( 4 ).NE. 0 ) THEN
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9988 ) SNAMES( K )
                  KFAIL( K ) = KFAIL( K ) + 1
                  ERRFLG = .TRUE.
               ELSE
                  IF( IAM.EQ.0 )
     $               WRITE( NOUT, FMT = 9987 ) SNAMES( K )
                  KPASS( K ) = KPASS( K ) + 1
               END IF
*
*              Dump matrix if IVERB >= 1 and error.
*
               IF( IVERB.GE.1 .AND. ERRFLG ) THEN
                  IF( IERR( 3 ).NE.0 .OR. IVERB.GE.3 ) THEN
                     CALL PDMPRNT( ICTXT, NOUT, MX, NX, MEM( IPMATX ),
     $                             LDX, 0, 0, 'SERIAL_X' )
                     CALL PB_PDLAPRNT( MX, NX, MEM( IPX ), 1, 1, DESCX,
     $                                 0, 0, 'PARALLEL_X', NOUT,
     $                                 MEM( IPMATX ) )
                  ELSE IF( IERR( 1 ).NE.0 ) THEN
                     IF( N.GT.0 )
     $                  CALL PDVPRNT( ICTXT, NOUT, N,
     $                                MEM( IPMATX+IX-1+(JX-1)*LDX ),
     $                                INCX, 0, 0, 'SERIAL_X' )
                     IF( INCX.EQ.DESCX( M_ ) ) THEN
                        CALL PB_PDLAPRNT( 1, N, MEM( IPX ), IX, JX,
     $                                    DESCX, 0, 0, 'PARALLEL_X',
     $                                    NOUT, MEM( IPMATX ) )
                     ELSE
                        CALL PB_PDLAPRNT( N, 1, MEM( IPX ), IX, JX,
     $                                    DESCX, 0, 0, 'PARALLEL_X',
     $                                    NOUT, MEM( IPMATX ) )
                     END IF
                  END IF
                  IF( YCHECK( K ) ) THEN
                     IF( IERR( 4 ).NE.0 .OR. IVERB.GE.3 ) THEN
                        CALL PDMPRNT( ICTXT, NOUT, MY, NY,
     $                                MEM( IPMATY ), LDY, 0, 0,
     $                                'SERIAL_Y' )
                        CALL PB_PDLAPRNT( MY, NY, MEM( IPY ), 1, 1,
     $                                    DESCY, 0, 0, 'PARALLEL_Y',
     $                                    NOUT, MEM( IPMATX ) )
                     ELSE IF( IERR( 2 ).NE.0 ) THEN
                        IF( N.GT.0 )
     $                     CALL PDVPRNT( ICTXT, NOUT, N,
     $                                   MEM( IPMATY+IY-1+(JY-1)*LDY ),
     $                                   INCY, 0, 0, 'SERIAL_Y' )
                        IF( INCY.EQ.DESCY( M_ ) ) THEN
                           CALL PB_PDLAPRNT( 1, N, MEM( IPY ), IY, JY,
     $                                       DESCY, 0, 0, 'PARALLEL_Y',
     $                                       NOUT, MEM( IPMATX ) )
                        ELSE
                           CALL PB_PDLAPRNT( N, 1, MEM( IPY ), IY, JY,
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
               WRITE( NOUT, FMT = 9985 ) J
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
         WRITE( NOUT, FMT = 9981 )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9983 )
         WRITE( NOUT, FMT = 9982 )
*
         DO 90 I = 1, NSUBS
            WRITE( NOUT, FMT = 9984 ) '|', SNAMES( I ), KTESTS( I ),
     $                                KPASS( I ), KFAIL( I ), KSKIP( I )
   90    CONTINUE
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9980 )
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
 9995 FORMAT( 2X, '---------------------------------------------------',
     $        '--------------------------' )
 9994 FORMAT( 2X, '     N     IX     JX     MX     NX  IMBX  INBX',
     $        '   MBX   NBX RSRCX CSRCX   INCX' )
 9993 FORMAT( 2X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,I5,1X,I5,1X,I5,1X,I5,1X,
     $        I5,1X,I5,1X,I6 )
 9992 FORMAT( 2X, '     N     IY     JY     MY     NY  IMBY  INBY',
     $        '   MBY   NBY RSRCY CSRCY   INCY' )
 9991 FORMAT( 'Not enough memory for this test: going on to',
     $        ' next test case.' )
 9990 FORMAT( 'Not enough memory. Need: ', I12 )
 9989 FORMAT( 2X, '   Tested Subroutine: ', A )
 9988 FORMAT( 2X, '   ***** Computational check: ', A, '       ',
     $        ' FAILED ',' *****' )
 9987 FORMAT( 2X, '   ***** Computational check: ', A, '       ',
     $        ' PASSED ',' *****' )
 9986 FORMAT( 2X, '   ***** ERROR ***** Matrix operand ', A,
     $        ' modified by ', A, ' *****' )
 9985 FORMAT( 2X, 'Test number ', I4, ' completed.' )
 9984 FORMAT( 2X,A1,2X,A7,8X,I4,6X,I4,5X,I4,4X,I4 )
 9983 FORMAT( 2X, '   SUBROUTINE  TOTAL TESTS  PASSED   FAILED  ',
     $        'SKIPPED' )
 9982 FORMAT( 2X, '   ----------  -----------  ------   ------  ',
     $        '-------' )
 9981 FORMAT( 2X, 'Testing Summary')
 9980 FORMAT( 2X, 'End of Tests.' )
 9979 FORMAT( 2X, 'Tests started.' )
 9978 FORMAT( 2X, '   ***** Operation not supported, error code: ',
     $        I5, ' *****' )
*
      STOP
*
*     End of PDBLA1TST
*
      END
      SUBROUTINE PDBLA1TSTINFO( SUMMRY, NOUT, NMAT, NVAL, MXVAL,
     $                          NXVAL, IMBXVAL, MBXVAL, INBXVAL,
     $                          NBXVAL, RSCXVAL, CSCXVAL, IXVAL,
     $                          JXVAL, INCXVAL, MYVAL, NYVAL, IMBYVAL,
     $                          MBYVAL, INBYVAL, NBYVAL, RSCYVAL,
     $                          CSCYVAL, IYVAL, JYVAL, INCYVAL,
     $                          LDVAL, NGRIDS, PVAL, LDPVAL, QVAL,
     $                          LDQVAL, LTEST, SOF, TEE, IAM, IGAP,
     $                          IVERB, NPROCS, ALPHA, WORK )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      LOGICAL            SOF, TEE
      INTEGER            IAM, IGAP, IVERB, LDPVAL, LDQVAL, LDVAL,
     $                   NGRIDS, NMAT, NOUT, NPROCS
      DOUBLE PRECISION   ALPHA
*     ..
*     .. Array Arguments ..
      CHARACTER*( * )    SUMMRY
      LOGICAL            LTEST( * )
      INTEGER            CSCXVAL( LDVAL ), CSCYVAL( LDVAL ),
     $                   IMBXVAL( LDVAL ), IMBYVAL( LDVAL ),
     $                   INBXVAL( LDVAL ), INBYVAL( LDVAL ),
     $                   INCXVAL( LDVAL ), INCYVAL( LDVAL ),
     $                   IXVAL( LDVAL ), IYVAL( LDVAL ), JXVAL( LDVAL ),
     $                   JYVAL( LDVAL ), MBXVAL( LDVAL ),
     $                   MBYVAL( LDVAL ), MXVAL( LDVAL ),
     $                   MYVAL( LDVAL ), NBXVAL( LDVAL ),
     $                   NBYVAL( LDVAL ), NVAL( LDVAL ), NXVAL( LDVAL ),
     $                   NYVAL( LDVAL ), PVAL( LDPVAL ), QVAL( LDQVAL ),
     $                   RSCXVAL( LDVAL ), RSCYVAL( LDVAL ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDBLA1TSTINFO  get the needed startup information for testing various
*  Level 1 PBLAS routines, and transmits it to all processes.
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
*  NVAL    (global output) INTEGER array
*          On entry, NVAL is an array of dimension LDVAL.  On exit, this
*          array contains the values of N to run the code with.
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
*          lues that can be used for  DESCX(:),  IX, JX, INCX, DESCY(:),
*          IY,  JY  and  INCY.  This  is also the maximum number of test
*          cases.
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
*  LTEST   (global output) LOGICAL array
*          On entry,  LTEST  is an array of dimension at least eight. On
*          exit, if LTEST( i ) is .TRUE., the i-th Level 1 PBLAS routine
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
*  ALPHA   (global output) DOUBLE PRECISION
*          On exit, ALPHA specifies the value of alpha to be used in all
*          the test cases.
*
*  WORK    (local workspace) INTEGER array
*          On   entry,   WORK   is   an  array  of  dimension  at  least
*          MAX( 2, 2*NGRIDS+23*NMAT+NSUBS+4 )  with  NSUBS  equal  to 8.
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
      PARAMETER          ( NIN = 11, NSUBS = 8 )
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
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
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
         OPEN( NIN, FILE='PDBLAS1TST.dat', STATUS='OLD' )
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
*        Get number of grids
*
         READ( NIN, FMT = * ) NGRIDS
         IF( NGRIDS.LT.1 .OR. NGRIDS.GT.LDPVAL ) THEN
            WRITE( NOUT, FMT = 9998 ) 'Grids', LDPVAL
            GO TO 100
         ELSE IF( NGRIDS.GT.LDQVAL ) THEN
            WRITE( NOUT, FMT = 9998 ) 'Grids', LDQVAL
            GO TO 100
         END IF
*
*        Get values of P and Q
*
         READ( NIN, FMT = * ) ( PVAL( I ), I = 1, NGRIDS )
         READ( NIN, FMT = * ) ( QVAL( I ), I = 1, NGRIDS )
*
*        Read ALPHA
*
         READ( NIN, FMT = * ) ALPHA
*
*        Read number of tests.
*
         READ( NIN, FMT = * ) NMAT
         IF( NMAT.LT.1 .OR. NMAT.GT.LDVAL ) THEN
            WRITE( NOUT, FMT = 9998 ) 'Tests', LDVAL
            GO TO 100
         END IF
*
*        Read in input data into arrays.
*
         READ( NIN, FMT = * ) ( NVAL( I ),     I = 1, NMAT )
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
         GO TO 100
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
         CALL DGEBS2D( ICTXT, 'All', ' ', 1, 1, ALPHA, 1 )
*
         WORK( 1 ) = NGRIDS
         WORK( 2 ) = NMAT
         CALL IGEBS2D( ICTXT, 'All', ' ', 2, 1, WORK, 2 )
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
         CALL ICOPY( NGRIDS, PVAL,     1, WORK( I ), 1 )
         I = I + NGRIDS
         CALL ICOPY( NGRIDS, QVAL,     1, WORK( I ), 1 )
         I = I + NGRIDS
         CALL ICOPY( NMAT,   NVAL,     1, WORK( I ), 1 )
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
         DO 70 J = 1, NSUBS
            IF( LTEST( J ) ) THEN
               WORK( I ) = 1
            ELSE
               WORK( I ) = 0
            END IF
            I = I + 1
   70    CONTINUE
         I = I - 1
         CALL IGEBS2D( ICTXT, 'All', ' ', I, 1, WORK, I )
*
*        regurgitate input
*
         WRITE( NOUT, FMT = 9999 ) 'Level 1 PBLAS testing program.'
         WRITE( NOUT, FMT = 9999 ) USRINFO
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'Tests of the real double precision '//
     $               'Level 1 PBLAS'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'The following parameter values will be used:'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9993 ) NMAT
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
         WRITE( NOUT, FMT = 9982 ) ALPHA
         IF( LTEST( 1 ) ) THEN
            WRITE( NOUT, FMT = 9985 ) SNAMES( 1 ), ' ... Yes'
         ELSE
            WRITE( NOUT, FMT = 9985 ) SNAMES( 1 ), ' ... No '
         END IF
         DO 80 I = 2, NSUBS
            IF( LTEST( I ) ) THEN
               WRITE( NOUT, FMT = 9984 ) SNAMES( I ), ' ... Yes'
            ELSE
               WRITE( NOUT, FMT = 9984 ) SNAMES( I ), ' ... No '
            END IF
   80    CONTINUE
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
         CALL DGEBR2D( ICTXT, 'All', ' ', 1, 1, ALPHA, 1, 0, 0 )
*
         CALL IGEBR2D( ICTXT, 'All', ' ', 2, 1, WORK, 2, 0, 0 )
         NGRIDS = WORK( 1 )
         NMAT   = WORK( 2 )
*
         I = 2*NGRIDS + 23*NMAT + NSUBS + 4
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
         CALL ICOPY( NGRIDS, WORK( I ), 1, PVAL,     1 )
         I = I + NGRIDS
         CALL ICOPY( NGRIDS, WORK( I ), 1, QVAL,     1 )
         I = I + NGRIDS
         CALL ICOPY( NMAT,   WORK( I ), 1, NVAL,     1 )
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
         DO 90 J = 1, NSUBS
            IF( WORK( I ).EQ.1 ) THEN
               LTEST( J ) = .TRUE.
            ELSE
               LTEST( J ) = .FALSE.
            END IF
            I = I + 1
   90    CONTINUE
*
      END IF
*
      CALL BLACS_GRIDEXIT( ICTXT )
*
      RETURN
*
  100 WRITE( NOUT, FMT = 9997 )
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
*
*     End of PDBLA1TSTINFO
*
      END
      SUBROUTINE PDBLAS1TSTCHKE( LTEST, INOUT, NPROCS )
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
*  PDBLAS1TSTCHKE tests the error exits of the Level 1 PBLAS.
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
*  LTEST   (global input) LOGICAL array
*          On entry, LTEST is an array of dimension at least 8 (NSUBS).
*             If LTEST( 1 ) is .TRUE., PDSWAP will be tested;
*             If LTEST( 2 ) is .TRUE., PDSCAL will be tested;
*             If LTEST( 3 ) is .TRUE., PDCOPY will be tested;
*             If LTEST( 4 ) is .TRUE., PDAXPY will be tested;
*             If LTEST( 5 ) is .TRUE., PDDOT  will be tested;
*             If LTEST( 6 ) is .TRUE., PDNRM2 will be tested;
*             If LTEST( 7 ) is .TRUE., PDASUM will be tested;
*             If LTEST( 8 ) is .TRUE., PDAMAX will be tested.
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
      PARAMETER          ( NSUBS = 8 )
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
     $                   BLACS_GRIDINIT, PDAMAX, PDASUM, PDAXPY, PDCOPY,
     $                   PDDIMEE, PDDOT, PDNRM2, PDSCAL, PDSWAP,
     $                   PDVECEE
*     ..
*     .. Common Blocks ..
      LOGICAL            ABRTFLG
      INTEGER            NOUT
      CHARACTER*7        SNAMES( NSUBS )
      COMMON             /SNAMEC/SNAMES
      COMMON             /PBERRORC/NOUT, ABRTFLG
*     ..
*     .. Data Statements ..
      DATA               SCODE/11, 12, 11, 13, 13, 15, 15, 14/
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
*     Test PDSWAP
*
      I = 1
      IF( LTEST( I ) ) THEN
         CALL PDDIMEE( ICTXT, NOUT, PDSWAP, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDSWAP, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDSCAL
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDDIMEE( ICTXT, NOUT, PDSCAL, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDSCAL, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDCOPY
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDDIMEE( ICTXT, NOUT, PDCOPY, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDCOPY, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDAXPY
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDDIMEE( ICTXT, NOUT, PDAXPY, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDAXPY, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDDOT
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDDIMEE( ICTXT, NOUT, PDDOT, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDDOT, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDNRM2
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDDIMEE( ICTXT, NOUT, PDNRM2, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDNRM2, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDASUM
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDDIMEE( ICTXT, NOUT, PDASUM, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDASUM, SCODE( I ), SNAMES( I ) )
      END IF
*
*     Test PDAMAX
*
      I = I + 1
      IF( LTEST( I ) ) THEN
         CALL PDDIMEE( ICTXT, NOUT, PDAMAX, SCODE( I ), SNAMES( I ) )
         CALL PDVECEE( ICTXT, NOUT, PDAMAX, SCODE( I ), SNAMES( I ) )
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
*     End of PDBLAS1TSTCHKE
*
      END
      SUBROUTINE PDCHKARG1( ICTXT, NOUT, SNAME, N, ALPHA, IX, JX,
     $                      DESCX, INCX, IY, JY, DESCY, INCY, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, INCX, INCY, INFO, IX, IY, JX, JY, N,
     $                   NOUT
      DOUBLE PRECISION   ALPHA
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      SNAME
      INTEGER            DESCX( * ), DESCY( * )
*     ..
*
*  Purpose
*  =======
*
*  PDCHKARG1 checks the input-only arguments of the Level 1 PBLAS.  When
*  INFO = 0, this routine makes a copy of its arguments (which are INPUT
*  only arguments to PBLAS routines). Otherwise, it verifies the  values
*  of these arguments against the saved copies.
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
*  SNAME   (global input) CHARACTER*(*)
*          On entry, SNAME specifies the subroutine  name  calling  this
*          subprogram.
*
*  N       (global input) INTEGER
*          On entry, N specifies the length of the subvector operands.
*
*  ALPHA   (global input) DOUBLE PRECISION
*          On entry, ALPHA specifies the scalar alpha.
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
      INTEGER            I, INCXREF, INCYREF, IXREF, IYREF, JXREF,
     $                   JYREF, MYCOL, MYROW, NPCOL, NPROW, NREF
      DOUBLE PRECISION   ALPHAREF
*     ..
*     .. Local Arrays ..
      CHARACTER*15       ARGNAME
      INTEGER            DESCXREF( DLEN_ ), DESCYREF( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGSUM2D
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
         NREF = N
         IXREF = IX
         JXREF = JX
         DO 10 I = 1, DLEN_
            DESCXREF( I ) = DESCX( I )
   10    CONTINUE
         INCXREF = INCX
         IYREF = IY
         JYREF = JY
         DO 20 I = 1, DLEN_
            DESCYREF( I ) = DESCY( I )
   20    CONTINUE
         INCYREF = INCY
         ALPHAREF = ALPHA
*
      ELSE
*
*        Test saved args. Return with first mismatch.
*
         ARGNAME = ' '
         IF( N.NE.NREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'N'
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
         ELSE IF( ALPHA.NE.ALPHAREF ) THEN
            WRITE( ARGNAME, FMT = '(A)' ) 'ALPHA'
         ELSE
            INFO = 0
         END IF
*
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, 0 )
*
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
*
            IF( INFO.GT.0 ) THEN
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
*     End of PDCHKARG1
*
      END
      LOGICAL FUNCTION PISINSCOPE( ICTXT, N, IX, JX, DESCX, INCX )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, INCX, IX, JX, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCX( * )
*     ..
*
*  Purpose
*  =======
*
*  PISINSCOPE returns  .TRUE.  if the calling process is in the scope of
*  sub( X ) = X( IX+(JX-1)*DESCX(M_)+(i-1)*INCX ) and  .FALSE.  if it is
*  not.  This  routine is used to determine which processes should check
*  the answer returned by some Level 1 PBLAS routines.
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
*  N       (global input) INTEGER
*          The length of the subvector sub( X ).
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
      LOGICAL            COLREP, ROWREP
      INTEGER            IIX, IXCOL, IXROW, JJX, MYCOL, MYROW, NPCOL,
     $                   NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PB_INFOG2L
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      CALL PB_INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIX, JJX, IXROW, IXCOL )
      ROWREP = ( IXROW.EQ.-1 )
      COLREP = ( IXCOL.EQ.-1 )
*
      IF( DESCX( M_ ).EQ.1 .AND. N.EQ.1 ) THEN
*
*        This is the special case, find process owner of IX, JX, and
*        only this process is the scope.
*
         PISINSCOPE = ( ( IXROW.EQ.MYROW .OR. ROWREP ) .AND.
     $                   ( IXCOL.EQ.MYCOL .OR. COLREP ) )
*
      ELSE
*
         IF( INCX.EQ.DESCX( M_ ) ) THEN
*
*           row vector
*
            PISINSCOPE = ( MYROW.EQ.IXROW .OR. ROWREP )
*
         ELSE
*
*           column vector
*
            PISINSCOPE = ( MYCOL.EQ.IXCOL .OR. COLREP )
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PISINSCOPE
*
      END
      SUBROUTINE PDBLAS1TSTCHK( ICTXT, NOUT, NROUT, N, PSCLR, PUSCLR,
     $                          PISCLR, X, PX, IX, JX, DESCX, INCX, Y,
     $                          PY, IY, JY, DESCY, INCY, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, INCX, INCY, INFO, IX, IY, JX, JY, N,
     $                   NOUT, NROUT, PISCLR
      DOUBLE PRECISION   PSCLR, PUSCLR
*     ..
*     .. Array Arguments ..
      INTEGER            DESCX( * ), DESCY( * )
      DOUBLE PRECISION   PX( * ), PY( * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  PDBLAS1TSTCHK performs the computational tests of the Level 1 PBLAS.
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
*             If NROUT = 1,      PDSWAP will be tested;
*             else if NROUT = 2, PDSCAL will be tested;
*             else if NROUT = 3, PDCOPY will be tested;
*             else if NROUT = 4, PDAXPY will be tested;
*             else if NROUT = 5, PDDOT  will be tested;
*             else if NROUT = 6, PDNRM2 will be tested;
*             else if NROUT = 7, PDASUM will be tested;
*             else if NROUT = 8, PDAMAX will be tested.
*
*  N       (global input) INTEGER
*          On entry, N specifies the length of the subvector operands.
*
*  PSCLR   (global input) DOUBLE PRECISION
*          On entry, depending on the value of  NROUT,  PSCLR  specifies
*          the scalar ALPHA, or the output scalar returned by the PBLAS,
*          i.e., the dot product, the 2-norm,  the  absolute sum  or the
*          value of AMAX.
*
*  PUSCLR  (global input) DOUBLE PRECISION
*          On entry, PUSCLR specifies the real part of the  scalar ALPHA
*          used  by  the  real  scaling, the 2-norm, or the absolute sum
*          routines.  PUSCLR  is  not  used in the real versions of this
*          routine.
*
*  PISCLR  (global input) DOUBLE PRECISION
*          On entry, PISCLR  specifies the value of the global index re-
*          turned by PDAMAX, otherwise PISCLR is not used.
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
*  INFO    (global output) INTEGER
*          On exit, if INFO = 0,  no  error  has  been  found, otherwise
*          if( MOD( INFO,   2 ) = 1 ) then an error on X has been found,
*          if( MOD( INFO/2, 2 ) = 1 ) then an error on Y has been found.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, INXSCOPE, INYSCOPE, ROWREP
      INTEGER            I, IB, ICURCOL, ICURROW, IDUMM, IIX, IIY, IN,
     $                   IOFFX, IOFFY, ISCLR, IXCOL, IXROW, IYCOL,
     $                   IYROW, J, JB, JJX, JJY, JN, KK, LDX, LDY,
     $                   MYCOL, MYROW, NPCOL, NPROW
      DOUBLE PRECISION   ERR, ERRMAX, PREC, SCLR, USCLR
*     ..
*     .. Local Arrays ..
      INTEGER            IERR( 6 )
      CHARACTER*5        ARGIN1, ARGIN2, ARGOUT1, ARGOUT2
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DCOPY, DSWAP, IGAMX2D,
     $                   PB_INFOG2L, PDCHKVIN, PDERRASUM, PDERRAXPY,
     $                   PDERRDOT, PDERRNRM2, PDERRSCAL
*     ..
*     .. External Functions ..
      LOGICAL            PISINSCOPE
      INTEGER            IDAMAX
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           IDAMAX, PDLAMCH, PISINSCOPE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      INFO    = 0
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      ARGIN1  = '     '
      ARGIN2  = '     '
      ARGOUT1 = '     '
      ARGOUT2 = '     '
      DO 10 I = 1, 6
         IERR( I ) = 0
   10 CONTINUE
*
      PREC = PDLAMCH( ICTXT, 'precision' )
*
      IF( NROUT.EQ.1 ) THEN
*
*        Test PDSWAP
*
         IOFFX = IX + ( JX - 1 ) * DESCX( M_ )
         IOFFY = IY + ( JY - 1 ) * DESCY( M_ )
         CALL DSWAP( N, X( IOFFX ), INCX, Y( IOFFY ), INCY )
         CALL PDCHKVIN( ERRMAX, N, X, PX, IX, JX, DESCX, INCX,
     $                    IERR( 1 ) )
         CALL PDCHKVIN( ERRMAX, N, Y, PY, IY, JY, DESCY, INCY,
     $                    IERR( 2 ) )
*
      ELSE IF( NROUT.EQ.2 ) THEN
*
*        Test PDSCAL
*
         LDX   = DESCX( LLD_ )
         IOFFX = IX + ( JX - 1 ) * DESCX( M_ )
         CALL PB_INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIX, JJX, IXROW, IXCOL )
         ICURROW = IXROW
         ICURCOL = IXCOL
         ROWREP = ( IXROW.EQ.-1 )
         COLREP = ( IXCOL.EQ.-1 )
*
         IF( INCX.EQ.DESCX( M_ ) ) THEN
*
*           sub( X ) is a row vector
*
            JB = DESCX( INB_ ) - JX + 1
            IF( JB.LE.0 )
     $         JB = ( (-JB ) / DESCX( NB_ ) + 1 ) * DESCX( NB_ ) + JB
            JB = MIN( JB, N )
            JN = JX + JB - 1
*
            DO 20 J = JX, JN
*
               CALL PDERRSCAL( ERR, PSCLR, X( IOFFX ), PREC )
*
               IF( ( MYROW.EQ.ICURROW .OR. ROWREP ) .AND.
     $             ( MYCOL.EQ.ICURCOL .OR. COLREP ) ) THEN
                  IF( ABS( PX( IIX+(JJX-1)*LDX ) - X( IOFFX ) ).GT.
     $                ERR )
     $             IERR( 1 ) = 1
                  JJX = JJX + 1
               END IF
*
               IOFFX = IOFFX + INCX
*
   20       CONTINUE
*
            ICURCOL = MOD( ICURCOL+1, NPCOL )
*
            DO 40 J = JN+1, JX+N-1, DESCX( NB_ )
               JB = MIN( JX+N-J, DESCX( NB_ ) )
*
               DO 30 KK = 0, JB-1
*
                  CALL PDERRSCAL( ERR, PSCLR, X( IOFFX ), PREC )
*
                  IF( ( MYROW.EQ.ICURROW .OR. ROWREP ) .AND.
     $                ( MYCOL.EQ.ICURCOL .OR. COLREP ) ) THEN
                     IF( ABS( PX( IIX+(JJX-1)*LDX ) - X( IOFFX ) ).GT.
     $                   ERR )
     $                  IERR( 1 ) = 1
                     JJX = JJX + 1
                  END IF
*
                  IOFFX = IOFFX + INCX
*
   30          CONTINUE
*
               ICURCOL = MOD( ICURCOL+1, NPCOL )
*
   40       CONTINUE
*
         ELSE
*
*           sub( X ) is a column vector
*
            IB = DESCX( IMB_ ) - IX + 1
            IF( IB.LE.0 )
     $         IB = ( (-IB ) / DESCX( MB_ ) + 1 ) * DESCX( MB_ ) + IB
            IB = MIN( IB, N )
            IN = IX + IB - 1
*
            DO 50 I = IX, IN
*
               CALL PDERRSCAL( ERR, PSCLR, X( IOFFX ), PREC )
*
               IF( ( MYROW.EQ.ICURROW .OR. ROWREP ) .AND.
     $             ( MYCOL.EQ.ICURCOL .OR. COLREP ) ) THEN
                  IF( ABS( PX( IIX+(JJX-1)*LDX ) - X( IOFFX ) ).GT.
     $                ERR )
     $               IERR( 1 ) = 1
                  IIX = IIX + 1
               END IF
*
               IOFFX = IOFFX + INCX
*
   50       CONTINUE
*
            ICURROW = MOD( ICURROW+1, NPROW )
*
            DO 70 I = IN+1, IX+N-1, DESCX( MB_ )
               IB = MIN( IX+N-I, DESCX( MB_ ) )
*
               DO 60 KK = 0, IB-1
*
                  CALL PDERRSCAL( ERR, PSCLR, X( IOFFX ), PREC )
*
                  IF( ( MYROW.EQ.ICURROW .OR. ROWREP ) .AND.
     $                ( MYCOL.EQ.ICURCOL .OR. COLREP ) ) THEN
                     IF( ABS( PX( IIX+(JJX-1)*LDX ) - X( IOFFX ) ).GT.
     $                   ERR )
     $                  IERR( 1 ) = 1
                     IIX = IIX + 1
                  END IF
*
                  IOFFX = IOFFX + INCX
   60          CONTINUE
*
               ICURROW = MOD( ICURROW+1, NPROW )
*
   70       CONTINUE
*
         END IF
*
      ELSE IF( NROUT.EQ.3 ) THEN
*
*        Test PDCOPY
*
         IOFFX = IX + ( JX - 1 ) * DESCX( M_ )
         IOFFY = IY + ( JY - 1 ) * DESCY( M_ )
         CALL DCOPY( N, X( IOFFX ), INCX, Y( IOFFY ), INCY )
         CALL PDCHKVIN( ERRMAX, N, X, PX, IX, JX, DESCX, INCX,
     $                  IERR( 1 ) )
         CALL PDCHKVIN( ERRMAX, N, Y, PY, IY, JY, DESCY, INCY,
     $                  IERR( 2 ) )
*
      ELSE IF( NROUT.EQ.4 ) THEN
*
*        Test PDAXPY
*
         CALL PDCHKVIN( ERRMAX, N, X, PX, IX, JX, DESCX, INCX,
     $                  IERR( 1 ) )
         LDY = DESCY( LLD_ )
         IOFFX = IX + ( JX - 1 ) * DESCX( M_ )
         IOFFY = IY + ( JY - 1 ) * DESCY( M_ )
         CALL PB_INFOG2L( IY, JY, DESCY, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIY, JJY, IYROW, IYCOL )
         ICURROW = IYROW
         ICURCOL = IYCOL
         ROWREP  = ( IYROW.EQ.-1 )
         COLREP  = ( IYCOL.EQ.-1 )
*
         IF( INCY.EQ.DESCY( M_ ) ) THEN
*
*           sub( Y ) is a row vector
*
            JB = DESCY( INB_ ) - JY + 1
            IF( JB.LE.0 )
     $         JB = ( (-JB ) / DESCY( NB_ ) + 1 ) * DESCY( NB_ ) + JB
            JB = MIN( JB, N )
            JN = JY + JB - 1
*
            DO 140 J = JY, JN
*
               CALL PDERRAXPY( ERR, PSCLR, X( IOFFX ), Y( IOFFY ),
     $                         PREC )
*
               IF( ( MYROW.EQ.ICURROW .OR. ROWREP ) .AND.
     $             ( MYCOL.EQ.ICURCOL .OR. COLREP ) ) THEN
                  IF( ABS( PY( IIY+(JJY-1)*LDY ) - Y( IOFFY ) ).GT.
     $                ERR ) THEN
                     IERR( 2 ) = 1
                  END IF
                  JJY = JJY + 1
               END IF
*
               IOFFX = IOFFX + INCX
               IOFFY = IOFFY + INCY
*
  140       CONTINUE
*
            ICURCOL = MOD( ICURCOL+1, NPCOL )
*
            DO 160 J = JN+1, JY+N-1, DESCY( NB_ )
               JB = MIN( JY+N-J, DESCY( NB_ ) )
*
               DO 150 KK = 0, JB-1
*
                  CALL PDERRAXPY( ERR, PSCLR, X( IOFFX ), Y( IOFFY ),
     $                            PREC )
*
                  IF( ( MYROW.EQ.ICURROW .OR. ROWREP ) .AND.
     $                ( MYCOL.EQ.ICURCOL .OR. COLREP ) ) THEN
                     IF( ABS( PY( IIY+(JJY-1)*LDY ) - Y( IOFFY ) ).GT.
     $                   ERR ) THEN
                        IERR( 2 ) = 1
                     END IF
                     JJY = JJY + 1
                  END IF
*
                  IOFFX = IOFFX + INCX
                  IOFFY = IOFFY + INCY
*
  150          CONTINUE
*
               ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  160       CONTINUE
*
         ELSE
*
*           sub( Y ) is a column vector
*
            IB = DESCY( IMB_ ) - IY + 1
            IF( IB.LE.0 )
     $         IB = ( (-IB ) / DESCY( MB_ ) + 1 ) * DESCY( MB_ ) + IB
            IB = MIN( IB, N )
            IN = IY + IB - 1
*
            DO 170 I = IY, IN
*
               CALL PDERRAXPY( ERR, PSCLR, X( IOFFX ), Y( IOFFY ),
     $                         PREC )
*
               IF( ( MYROW.EQ.ICURROW .OR. ROWREP ) .AND.
     $             ( MYCOL.EQ.ICURCOL .OR. COLREP ) ) THEN
                  IF( ABS( PY( IIY+(JJY-1)*LDY ) - Y( IOFFY ) ).GT.
     $                ERR ) THEN
                     IERR( 2 ) = 1
                  END IF
                  IIY = IIY + 1
               END IF
*
               IOFFX = IOFFX + INCX
               IOFFY = IOFFY + INCY
*
  170       CONTINUE
*
            ICURROW = MOD( ICURROW+1, NPROW )
*
            DO 190 I = IN+1, IY+N-1, DESCY( MB_ )
               IB = MIN( IY+N-I, DESCY( MB_ ) )
*
               DO 180 KK = 0, IB-1
*
                  CALL PDERRAXPY( ERR, PSCLR, X( IOFFX ), Y( IOFFY ),
     $                            PREC )
*
                  IF( ( MYROW.EQ.ICURROW .OR. ROWREP ) .AND.
     $                ( MYCOL.EQ.ICURCOL .OR. COLREP ) ) THEN
                     IF( ABS( PY( IIY+(JJY-1)*LDY ) - Y( IOFFY ) ).GT.
     $                   ERR ) THEN
                        IERR( 2 ) = 1
                     END IF
                     IIY = IIY + 1
                  END IF
*
                  IOFFX = IOFFX + INCX
                  IOFFY = IOFFY + INCY
*
  180          CONTINUE
*
               ICURROW = MOD( ICURROW+1, NPROW )
*
  190       CONTINUE
*
         END IF
*
      ELSE IF( NROUT.EQ.5 ) THEN
*
*        Test PDDOT
*
         CALL PDCHKVIN( ERRMAX, N, X, PX, IX, JX, DESCX, INCX,
     $                  IERR( 1 ) )
         CALL PDCHKVIN( ERRMAX, N, Y, PY, IY, JY, DESCY, INCY,
     $                  IERR( 2 ) )
         IOFFX = IX + ( JX - 1 ) * DESCX( M_ )
         IOFFY = IY + ( JY - 1 ) * DESCY( M_ )
         CALL PDERRDOT( ERR, N, SCLR, X( IOFFX ), INCX, Y( IOFFY ),
     $                  INCY, PREC )
         INXSCOPE = PISINSCOPE( ICTXT, N, IX, JX, DESCX, INCX )
         INYSCOPE = PISINSCOPE( ICTXT, N, IY, JY, DESCY, INCY )
         IF( INXSCOPE.OR.INYSCOPE ) THEN
            IF( ABS( PSCLR - SCLR ).GT.ERR ) THEN
               IERR( 3 ) = 1
               WRITE( ARGIN1, FMT = '(A)' ) 'DOT'
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 9998 ) ARGIN1
                  WRITE( NOUT, FMT = 9996 ) SCLR, PSCLR
               END IF
            END IF
         ELSE
            SCLR = ZERO
            IF( PSCLR.NE.SCLR ) THEN
               IERR( 4 ) = 1
               WRITE( ARGOUT1, FMT = '(A)' ) 'DOT'
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 9997 ) ARGOUT1
                  WRITE( NOUT, FMT = 9996 ) SCLR, PSCLR
               END IF
            END IF
         END IF
*
      ELSE IF( NROUT.EQ.6 ) THEN
*
*        Test PDNRM2
*
         CALL PDCHKVIN( ERRMAX, N, X, PX, IX, JX, DESCX, INCX,
     $                  IERR( 1 ) )
         IOFFX = IX + ( JX - 1 ) * DESCX( M_ )
         CALL PDERRNRM2( ERR, N, USCLR, X( IOFFX ), INCX, PREC )
         IF( PISINSCOPE( ICTXT, N, IX, JX, DESCX, INCX ) ) THEN
            IF( ABS( PUSCLR - USCLR ).GT.ERR ) THEN
               IERR( 3 ) = 1
               WRITE( ARGIN1, FMT = '(A)' ) 'NRM2'
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 9998 ) ARGIN1
                  WRITE( NOUT, FMT = 9996 ) USCLR, PUSCLR
               END IF
            END IF
         ELSE
            USCLR = ZERO
            IF( PUSCLR.NE.USCLR ) THEN
               IERR( 4 ) = 1
               WRITE( ARGOUT1, FMT = '(A)' ) 'NRM2'
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 9997 ) ARGOUT1
                  WRITE( NOUT, FMT = 9996 ) USCLR, PUSCLR
               END IF
            END IF
         END IF
*
      ELSE IF( NROUT.EQ.7 ) THEN
*
*        Test PDASUM
*
         CALL PDCHKVIN( ERRMAX, N, X, PX, IX, JX, DESCX, INCX,
     $                  IERR( 1 ) )
         IOFFX = IX + ( JX - 1 ) * DESCX( M_ )
         CALL PDERRASUM( ERR, N, USCLR, X( IOFFX ), INCX, PREC )
         IF( PISINSCOPE( ICTXT, N, IX, JX, DESCX, INCX ) ) THEN
            IF( ABS( PUSCLR - USCLR ) .GT. ERR ) THEN
               IERR( 3 ) = 1
               WRITE( ARGIN1, FMT = '(A)' ) 'ASUM'
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 9998 ) ARGIN1
                  WRITE( NOUT, FMT = 9996 ) USCLR, PUSCLR
               END IF
            END IF
         ELSE
            USCLR = ZERO
            IF( PUSCLR.NE.USCLR ) THEN
               IERR( 4 ) = 1
               WRITE( ARGOUT1, FMT = '(A)' ) 'ASUM'
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 9997 ) ARGOUT1
                  WRITE( NOUT, FMT = 9996 ) USCLR, PUSCLR
               END IF
            END IF
         END IF
*
      ELSE IF( NROUT.EQ.8 ) THEN
*
*        Test PDAMAX
*
         CALL PDCHKVIN( ERRMAX, N, X, PX, IX, JX, DESCX, INCX,
     $                  IERR( 1 ) )
         IOFFX = IX + ( JX - 1 ) * DESCX( M_ )
         IF( PISINSCOPE( ICTXT, N, IX, JX, DESCX, INCX ) ) THEN
            ISCLR = IDAMAX( N, X( IOFFX ), INCX )
            IF( N.LT.1 ) THEN
               SCLR = ZERO
            ELSE IF( ( INCX.EQ.1 ).AND.( DESCX( M_ ).EQ.1 ).AND.
     $               ( N.EQ.1 ) ) THEN
               ISCLR = JX
               SCLR = X( IOFFX )
            ELSE IF( INCX.EQ.DESCX( M_ ) ) THEN
               ISCLR = JX + ISCLR - 1
               SCLR = X( IX + ( ISCLR - 1 ) * DESCX( M_ ) )
            ELSE
               ISCLR = IX + ISCLR - 1
               SCLR = X( ISCLR + ( JX - 1 ) * DESCX( M_ ) )
            END IF
*
            IF( PSCLR.NE.SCLR ) THEN
               IERR( 3 ) = 1
               WRITE( ARGIN1, FMT = '(A)' ) 'AMAX'
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 9998 ) ARGIN1
                  WRITE( NOUT, FMT = 9996 ) SCLR, PSCLR
               END IF
            END IF
*
            IF( PISCLR.NE.ISCLR ) THEN
               IERR( 5 ) = 1
               WRITE( ARGIN2, FMT = '(A)' ) 'INDX'
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 9998 ) ARGIN2
                  WRITE( NOUT, FMT = 9995 ) ISCLR, PISCLR
               END IF
            END IF
         ELSE
            ISCLR = 0
            SCLR  = ZERO
            IF( PSCLR.NE.SCLR ) THEN
               IERR( 4 ) = 1
               WRITE( ARGOUT1, FMT = '(A)' ) 'AMAX'
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 9997 ) ARGOUT1
                  WRITE( NOUT, FMT = 9996 ) SCLR, PSCLR
               END IF
            END IF
            IF( PISCLR.NE.ISCLR ) THEN
               IERR( 6 ) = 1
               WRITE( ARGOUT2, FMT = '(A)' ) 'INDX'
               IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                  WRITE( NOUT, FMT = 9997 ) ARGOUT2
                  WRITE( NOUT, FMT = 9995 ) ISCLR, PISCLR
               END IF
            END IF
         END IF
*
      END IF
*
*     Find IERR across all processes
*
      CALL IGAMX2D( ICTXT, 'All', ' ', 6, 1, IERR, 6, IDUMM, IDUMM, -1,
     $              -1, 0 )
*
*     Encode the errors found in INFO
*
      IF( IERR( 1 ).NE.0 ) THEN
         INFO = INFO + 1
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9999 ) 'X'
      END IF
*
      IF( IERR( 2 ).NE.0 ) THEN
         INFO = INFO + 2
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9999 ) 'Y'
      END IF
*
      IF( IERR( 3 ).NE.0 )
     $   INFO = INFO + 4
*
      IF( IERR( 4 ).NE.0 )
     $   INFO = INFO + 8
*
      IF( IERR( 5 ).NE.0 )
     $   INFO = INFO + 16
*
      IF( IERR( 6 ).NE.0 )
     $   INFO = INFO + 32
*
 9999 FORMAT( 2X, '   ***** ERROR: Vector operand ', A,
     $        ' is incorrect.' )
 9998 FORMAT( 2X, '   ***** ERROR: Output scalar result ', A,
     $        ' in scope is incorrect.' )
 9997 FORMAT( 2X, '   ***** ERROR: Output scalar result ', A,
     $        ' out of scope is incorrect.' )
 9996 FORMAT( 2X, '   ***** Expected value is: ', D30.18, /2X,
     $        '         Obtained value is: ', D30.18 )
 9995 FORMAT( 2X, '   ***** Expected value is: ', I6, /2X,
     $        '         Obtained value is: ', I6 )
*
      RETURN
*
*     End of PDBLAS1TSTCHK
*
      END
      SUBROUTINE PDERRDOT( ERRBND, N, SCLR, X, INCX, Y, INCY, PREC )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   ERRBND, PREC, SCLR
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  PDERRDOT  serially  computes  the  dot product X**T * Y and returns a
*  scaled relative acceptable error bound on the result.
*
*  Notes
*  =====
*
*  If dot1 = SCLR and  dot2 are two different computed results, and dot1
*  is being assumed to be correct, we require
*
*     abs( dot1 - dot2 ) <= ERRBND = ERRFACT * abs( dot1 ),
*
*  where ERRFACT is computed as the maximum of the positive and negative
*  partial  sums  multiplied  by  a constant proportional to the machine
*  precision.
*
*  Arguments
*  =========
*
*  ERRBND  (global output) DOUBLE PRECISION
*          On exit, ERRBND  specifies the scaled relative acceptable er-
*          ror bound.
*
*  N       (global input) INTEGER
*          On entry, N specifies the length of the vector operands.
*
*  SCLR    (global output) DOUBLE PRECISION
*          On exit,  SCLR  specifies  the dot product of the two vectors
*          X and Y.
*
*  X       (global input) DOUBLE PRECISION array
*          On   entry,   X   is   an   array   of   dimension  at  least
*          ( 1 + ( n - 1 )*abs( INCX ) ).  Before  entry,  the incremen-
*          ted array X must contain the vector x.
*
*  INCX    (global input) INTEGER.
*          On entry, INCX specifies the increment for the elements of X.
*          INCX must not be zero.
*
*  Y       (global input) DOUBLE PRECISION array
*          On   entry,   Y   is   an   array   of   dimension  at  least
*          ( 1 + ( n - 1 )*abs( INCY ) ).  Before  entry,  the incremen-
*          ted array Y must contain the vector y.
*
*  INCY    (global input) INTEGER.
*          On entry, INCY specifies the increment for the elements of Y.
*          INCY must not be zero.
*
*  PREC    (global input) DOUBLE PRECISION
*          On entry, PREC specifies the machine precision.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO, ZERO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0,
     $                   ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX, IY
      DOUBLE PRECISION   ADDBND, FACT, SUMNEG, SUMPOS, TMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      IX = 1
      IY = 1
      SCLR = ZERO
      SUMPOS = ZERO
      SUMNEG = ZERO
      FACT = TWO * ( ONE + PREC )
      ADDBND = TWO * TWO * TWO * PREC
*
      DO 10 I = 1, N
         TMP = X( IX ) * Y( IY )
         SCLR = SCLR + TMP
         IF( TMP.GE.ZERO ) THEN
            SUMPOS = SUMPOS + TMP * FACT
         ELSE
            SUMNEG = SUMNEG - TMP * FACT
         END IF
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
*
      ERRBND = ADDBND * MAX( SUMPOS, SUMNEG )
*
      RETURN
*
*     End of PDERRDOT
*
      END
      SUBROUTINE PDERRNRM2( ERRBND, N, USCLR, X, INCX, PREC )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ERRBND, PREC, USCLR
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  PDERRNRM2  serially  computes  the  2-norm the vector X and returns a
*  scaled relative acceptable error bound on the result.
*
*  Notes
*  =====
*
*  If  norm1 = SCLR  and  norm2  are two different computed results, and
*  norm1 being assumed to be correct, we require
*
*     abs( norm1 - norm2 ) <= ERRBND = ERRFACT * abs( norm1 ),
*
*  where ERRFACT is computed as the maximum of the positive and negative
*  partial  sums  multiplied  by  a constant proportional to the machine
*  precision.
*
*  Arguments
*  =========
*
*  ERRBND  (global output) DOUBLE PRECISION
*          On exit, ERRBND  specifies the scaled relative acceptable er-
*          ror bound.
*
*  N       (global input) INTEGER
*          On entry, N specifies the length of the vector operand.
*
*  USCLR   (global output) DOUBLE PRECISION
*          On exit, USCLR specifies the 2-norm of the vector X.
*
*  X       (global input) DOUBLE PRECISION array
*          On   entry,   X   is   an   array   of   dimension  at  least
*          ( 1 + ( n - 1 )*abs( INCX ) ).  Before  entry,  the incremen-
*          ted array X must contain the vector x.
*
*  INCX    (global input) INTEGER.
*          On entry, INCX specifies the increment for the elements of X.
*          INCX must not be zero.
*
*  PREC    (global input) DOUBLE PRECISION
*          On entry, PREC specifies the machine precision.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO, ZERO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0,
     $                   ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI, ADDBND, FACT, SCALE, SSQ, SUMSCA, SUMSSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      USCLR = ZERO
      SUMSSQ = ONE
      SUMSCA = ZERO
      ADDBND = TWO * TWO * TWO * PREC
      FACT = ONE + TWO * ( ( ONE + PREC )**3 - ONE )
*
      SCALE = ZERO
      SSQ = ONE
      DO 10 IX = 1, 1 + ( N - 1 )*INCX, INCX
         IF( X( IX ).NE.ZERO ) THEN
            ABSXI = ABS( X( IX ) )
            IF( SCALE.LT.ABSXI )THEN
               SUMSSQ = ONE + ( SSQ*( SCALE/ABSXI )**2 ) * FACT
               ERRBND = ADDBND * SUMSSQ
               SUMSSQ = SUMSSQ + ERRBND
               SSQ    = ONE + SSQ*( SCALE/ABSXI )**2
               SUMSCA = ABSXI
               SCALE  = ABSXI
            ELSE
               SUMSSQ = SSQ + ( ( ABSXI/SCALE )**2 ) * FACT
               ERRBND = ADDBND * SUMSSQ
               SUMSSQ = SUMSSQ + ERRBND
               SSQ    = SSQ + ( ABSXI/SCALE )**2
            END IF
         END IF
   10 CONTINUE
*
      USCLR = SCALE * SQRT( SSQ )
*
*     Error on square root
*
      ERRBND = SQRT( SUMSSQ ) * ( ONE + TWO * ( 1.00001D+0 * PREC ) )
*
      ERRBND = ( SUMSCA * ERRBND ) - USCLR
*
      RETURN
*
*     End of PDERRNRM2
*
      END
      SUBROUTINE PDERRASUM( ERRBND, N, USCLR, X, INCX, PREC )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ERRBND, PREC, USCLR
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  PDERRASUM  serially computes the sum of absolute values of the vector
*  X and returns a scaled relative acceptable error bound on the result.
*
*  Arguments
*  =========
*
*  ERRBND  (global output) DOUBLE PRECISION
*          On exit, ERRBND  specifies a scaled relative acceptable error
*          bound. In this case the error bound is just the absolute  sum
*          multiplied  by  a constant proportional to the machine preci-
*          sion.
*
*  N       (global input) INTEGER
*          On entry, N specifies the length of the vector operand.
*
*  USCLR   (global output) DOUBLE PRECISION
*          On exit, USCLR  specifies  the  sum of absolute values of the
*          vector X.
*
*  X       (global input) DOUBLE PRECISION array
*          On   entry,   X   is   an   array   of   dimension  at  least
*          ( 1 + ( n - 1 )*abs( INCX ) ).  Before  entry,  the incremen-
*          ted array X must contain the vector x.
*
*  INCX    (global input) INTEGER.
*          On entry, INCX specifies the increment for the elements of X.
*          INCX must not be zero.
*
*  PREC    (global input) DOUBLE PRECISION
*          On entry, PREC specifies the machine precision.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   TWO, ZERO
      PARAMETER          ( TWO = 2.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ADDBND
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IX = 1
      USCLR = ZERO
      ADDBND = TWO * TWO * TWO * PREC
*
      DO 10 IX = 1, 1 + ( N - 1 )*INCX, INCX
         USCLR = USCLR + ABS( X( IX ) )
   10 CONTINUE
*
      ERRBND = ADDBND * USCLR
*
      RETURN
*
*     End of PDERRASUM
*
      END
      SUBROUTINE PDERRSCAL( ERRBND, PSCLR, X, PREC )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ERRBND, PREC, PSCLR, X
*     ..
*
*  Purpose
*  =======
*
*  PDERRSCAL serially computes the product PSCLR * X and returns a sca-
*  led relative acceptable error bound on the result.
*
*  Notes
*  =====
*
*  If s1 = PSCLR*X and  s2 are two different computed results, and s1 is
*  being assumed to be correct, we require
*
*        abs( s1 - s2 ) <= ERRBND = ERRFACT * abs( s1 ),
*
*  where ERRFACT is computed as two times the machine precision.
*
*  Arguments
*  =========
*
*  ERRBND  (global output) DOUBLE PRECISION
*          On exit, ERRBND  specifies the scaled relative acceptable er-
*          ror bound.
*
*  PSCLR   (global input) DOUBLE PRECISION
*          On entry, PSCLR specifies the scale factor.
*
*  X       (global input/global output) DOUBLE PRECISION
*          On entry, X  specifies the scalar to be scaled. On exit, X is
*          the scaled entry.
*
*  PREC    (global input) DOUBLE PRECISION
*          On entry, PREC specifies the machine precision.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D+0 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      X = PSCLR * X
*
      ERRBND = ( TWO * PREC ) * ABS( X )
*
      RETURN
*
*     End of PDERRSCAL
*
      END
      SUBROUTINE PDERRAXPY( ERRBND, PSCLR, X, Y, PREC )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ERRBND, PREC, PSCLR, X, Y
*     ..
*
*  Purpose
*  =======
*
*  PDERRAXPY  serially computes Y := Y + PSCLR * X and returns a scaled
*  relative acceptable error bound on the result.
*
*  Arguments
*  =========
*
*  ERRBND  (global output) DOUBLE PRECISION
*          On exit, ERRBND  specifies the scaled relative acceptable er-
*          ror bound.
*
*  PSCLR   (global input) DOUBLE PRECISION
*          On entry, PSCLR specifies the scale factor.
*
*  X       (global input) DOUBLE PRECISION
*          On entry, X  specifies the scalar to be scaled.
*
*  Y       (global input/global output) DOUBLE PRECISION
*          On entry, Y specifies the scalar to be added. On exit, Y con-
*          tains the resulting scalar.
*
*  PREC    (global input) DOUBLE PRECISION
*          On entry, PREC specifies the machine precision.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO, ZERO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0,
     $                   ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   ADDBND, FACT, SUMPOS, SUMNEG, TMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      SUMPOS = ZERO
      SUMNEG = ZERO
      FACT = ONE + TWO * PREC
      ADDBND = TWO * TWO * TWO * PREC
*
      TMP = PSCLR * X
      IF( TMP.GE.ZERO ) THEN
         SUMPOS = SUMPOS + TMP * FACT
      ELSE
         SUMNEG = SUMNEG - TMP * FACT
      END IF
*
      TMP = Y
      IF( TMP.GE.ZERO ) THEN
         SUMPOS = SUMPOS + TMP
      ELSE
         SUMNEG = SUMNEG - TMP
      END IF
*
      Y = Y + ( PSCLR * X )
*
      ERRBND = ADDBND * MAX( SUMPOS, SUMNEG )
*
      RETURN
*
*     End of PDERRAXPY
*
      END
