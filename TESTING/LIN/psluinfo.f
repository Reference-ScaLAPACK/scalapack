      SUBROUTINE PSLUINFO( SUMMRY, NOUT, NMAT, MVAL, NVAL, LDNVAL, NNB,
     $                     NBVAL, LDNBVAL, NNR, NRVAL, LDNRVAL, NNBR,
     $                     NBRVAL, LDNBRVAL, NGRIDS, PVAL, LDPVAL, QVAL,
     $                     LDQVAL, THRESH, EST, WORK, IAM, NPROCS )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      LOGICAL            EST
      CHARACTER*( * )    SUMMRY
      INTEGER            IAM, LDNBRVAL, LDNBVAL, LDNRVAL, LDNVAL,
     $                   LDPVAL, LDQVAL, NGRIDS, NMAT, NNB, NNBR,
     $                   NPROCS, NNR, NOUT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            MVAL( LDNVAL ), NBRVAL( LDNBRVAL ),
     $                   NBVAL( LDNBVAL ), NRVAL( LDNRVAL ),
     $                   NVAL( LDNVAL ), PVAL( LDPVAL ), QVAL( LDQVAL ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLUINFO gets needed startup information for LU factorization
*  and transmits it to all processes.
*
*  Arguments
*  =========
*
*  SUMMRY   (global output) CHARACTER*(*)
*           Name of output (summary) file (if any). Only defined for
*           process 0.
*
*  NOUT     (global output) INTEGER
*           The unit number for output file. NOUT = 6, ouput to screen,
*           NOUT = 0, output to stderr.  Only defined for process 0.
*
*  NMAT     (global output) INTEGER
*           The number of different values that can be used for M and N.
*
*  MVAL     (global output) INTEGER array, dimension (LDNVAL)
*           The values of M (number of rows in matrix) to run the code
*           with.
*
*  NVAL     (global output) INTEGER array, dimension (LDNVAL)
*           The values of N (number of columns in matrix) to run the
*           code with.
*
*  LDNVAL   (global input) INTEGER
*           The maximum number of different values that can be used for
*           M and N, LDNVAL > =  NMAT.
*
*  NNB      (global output) INTEGER
*           The number of different values that can be used for NB.
*
*  NBVAL    (global output) INTEGER array, dimension (LDNBVAL)
*           The values of NB (blocksize) to run the code with.
*
*  LDNBVAL  (global input) INTEGER
*           The maximum number of different values that can be used for
*           NB, LDNBVAL >= NNB.
*
*  NNR      (global output) INTEGER
*           The number of different values that can be used for NRHS.
*
*  NRVAL    (global output) INTEGER array, dimension(LDNRVAL)
*           The values of NRHS (# of Right Hand Sides) to run the code
*           with.
*
*  LDNRVAL  (global input) INTEGER
*           The maximum number of different values that can be used for
*           NRHS, LDNRVAL >= NNR.
*
*  NNBR     (global output) INTEGER
*           The number of different values that can be used for NBRHS.
*
*  NBRVAL   (global output) INTEGER array, dimension (LDNBRVAL)
*           The values of NBRHS (RHS blocksize) to run the code with.
*
*  LDNBRVAL (global input) INTEGER
*           The maximum number of different values that can be used for
*           NBRHS, LDNBRVAL >= NBRVAL.
*
*  NGRIDS   (global output) INTEGER
*           The number of different values that can be used for P & Q.
*
*  PVAL     (global output) INTEGER array, dimension (LDPVAL)
*           The values of P (number of process rows) to run the code
*           with.
*
*  LDPVAL   (global input) INTEGER
*           The maximum number of different values that can be used for
*           P, LDPVAL >= NGRIDS.
*
*  QVAL     (global output) INTEGER array, dimension (LDQVAL)
*           The values of Q (number of process columns) to run the code
*           with.
*
*  LDQVAL   (global input) INTEGER
*           The maximum number of different values that can be used for
*           Q, LDQVAL >= NGRIDS.
*
*  THRESH   (global output) REAL
*           Indicates what error checks shall be run and printed out:
*           < 0 : Perform no error checking
*           > 0 : report all residuals greater than THRESH, perform
*                 factor check only if solve check fails
*
*  EST      (global output) LOGICAL
*           Flag indicating if condition estimation and iterative
*           refinement routines are to be exercised.
*
*  WORK     (local workspace) INTEGER array of dimension >=
*           MAX( 6, 2*LDNVAL+LDNBVAL+LDNRVAL+LDNBRVAL+LDPVAL+LDQVAL )
*           Used to pack all input arrays in order to send info in one
*           message.
*
*  IAM      (local input) INTEGER
*           My process number.
*
*  NPROCS   (global input) INTEGER
*           The total number of processes.
*
* ======================================================================
*
* Note: For packing the information we assumed that the length in bytes
* ===== of an integer is equal to the length in bytes of a real single
*       precision.
*
* ======================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            NIN
      PARAMETER          ( NIN = 11 )
*     ..
*     .. Local Scalars ..
      CHARACTER*79       USRINFO
      INTEGER            I, ICTXT
      REAL               EPS
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINIT, BLACS_SETUP, ICOPY, IGEBR2D,
     $                   IGEBS2D, SGEBR2D, SGEBS2D
*     ..
*     .. External Functions ..
      REAL               PSLAMCH
      EXTERNAL           PSLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
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
         OPEN( NIN, FILE='LU.dat', STATUS='OLD' )
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
*        Get number of matrices and their dimensions
*
         READ( NIN, FMT = * ) NMAT
         IF( NMAT.LT.1 .OR. NMAT.GT.LDNVAL ) THEN
            WRITE( NOUT, FMT = 9994 ) 'N', LDNVAL
            GO TO 20
         END IF
         READ( NIN, FMT = * ) ( MVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( NVAL( I ), I = 1, NMAT )
*
*        Get values of NB
*
         READ( NIN, FMT = * ) NNB
         IF( NNB.LT.1 .OR. NNB.GT.LDNBVAL ) THEN
            WRITE( NOUT, FMT = 9994 ) 'NB', LDNBVAL
            GO TO 20
         END IF
         READ( NIN, FMT = * ) ( NBVAL( I ), I = 1, NNB )
*
*        Get values of NRHS
*
         READ( NIN, FMT = * ) NNR
         IF( NNR.LT.1 .OR. NNR.GT.LDNRVAL ) THEN
            WRITE( NOUT, FMT = 9994 ) 'NRHS', LDNRVAL
            GO TO 20
         END IF
         READ( NIN, FMT = * ) ( NRVAL( I ), I = 1, NNR )
*
*        Get values of NBRHS
*
         READ( NIN, FMT = * ) NNBR
         IF( NNBR.LT.1 .OR. NNBR.GT.LDNBRVAL ) THEN
            WRITE( NOUT, FMT = 9994 ) 'NBRHS', LDNBRVAL
            GO TO 20
         END IF
         READ( NIN, FMT = * ) ( NBRVAL( I ), I = 1, NNBR )
*
*        Get number of grids
*
         READ( NIN, FMT = * ) NGRIDS
         IF( NGRIDS.LT.1 .OR. NGRIDS.GT.LDPVAL ) THEN
            WRITE( NOUT, FMT = 9994 ) 'Grids', LDPVAL
            GO TO 20
         ELSE IF( NGRIDS.GT.LDQVAL ) THEN
            WRITE( NOUT, FMT = 9994 ) 'Grids', LDQVAL
            GO TO 20
         END IF
*
*        Get values of P and Q
*
         READ( NIN, FMT = * ) ( PVAL( I ), I = 1, NGRIDS )
         READ( NIN, FMT = * ) ( QVAL( I ), I = 1, NGRIDS )
*
*        Get level of checking
*
         READ( NIN, FMT = * ) THRESH
*
*        Read the flag that indicates whether to test the condition
*        estimation and iterative refinement routines.
*
         READ( NIN, FMT = * ) EST
*
*        Close input file
*
         CLOSE( NIN )
*
*        For pvm only: if virtual machine not set up, allocate it and
*        spawn the correct number of processes.
*
         IF( NPROCS.LT.1 ) THEN
            NPROCS = 0
            DO 10 I = 1, NGRIDS
               NPROCS = MAX( NPROCS, PVAL( I )*QVAL( I ) )
   10       CONTINUE
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
         EPS = PSLAMCH( ICTXT, 'eps' )
*
*        Pack information arrays and broadcast
*
         CALL SGEBS2D( ICTXT, 'All', ' ', 1, 1, THRESH, 1 )
*
         WORK( 1 ) = NMAT
         WORK( 2 ) = NNB
         WORK( 3 ) = NNR
         WORK( 4 ) = NNBR
         WORK( 5 ) = NGRIDS
         IF( EST ) THEN
            WORK( 6 ) = 1
         ELSE
            WORK( 6 ) = 0
         END IF
         CALL IGEBS2D( ICTXT, 'All', ' ', 6, 1, WORK, 6 )
*
         I = 1
         CALL ICOPY( NMAT, MVAL, 1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT, NVAL, 1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NNB, NBVAL, 1, WORK( I ), 1 )
         I = I + NNB
         CALL ICOPY( NNR, NRVAL, 1, WORK( I ), 1 )
         I = I + NNR
         CALL ICOPY( NNBR, NBRVAL, 1, WORK( I ), 1 )
         I = I + NNBR
         CALL ICOPY( NGRIDS, PVAL, 1, WORK( I ), 1 )
         I = I + NGRIDS
         CALL ICOPY( NGRIDS, QVAL, 1, WORK( I ), 1 )
         I = I + NGRIDS - 1
         CALL IGEBS2D( ICTXT, 'All', ' ', I, 1, WORK, I )
*
*        regurgitate input
*
         WRITE( NOUT, FMT = 9999 )
     $               'ScaLAPACK Ax=b by LU factorization.'
         WRITE( NOUT, FMT = 9999 ) USRINFO
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'Tests of the parallel '//
     $               'real single precision LU factorization '//
     $               'and solve.'
         WRITE( NOUT, FMT = 9999 )
     $               'The following scaled residual '//
     $               'checks will be computed:'
         WRITE( NOUT, FMT = 9999 )
     $               ' Solve residual         = ||Ax - b|| / '//
     $               '(||x|| * ||A|| * eps * N)'
         WRITE( NOUT, FMT = 9999 )
     $               ' Factorization residual = ||A - LU|| / '//
     $               '(||A|| * eps * N)'
         WRITE( NOUT, FMT = 9999 )
     $               'The matrix A is randomly '//
     $               'generated for each test.'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'An explanation of the input/output '//
     $               'parameters follows:'
         WRITE( NOUT, FMT = 9999 )
     $               'TIME    : Indicates whether WALL or '//
     $               'CPU time was used.'
*
         WRITE( NOUT, FMT = 9999 )
     $               'M       : The number of rows in the '//
     $               'matrix A.'
         WRITE( NOUT, FMT = 9999 )
     $               'N       : The number of columns in the '//
     $               'matrix A.'
         WRITE( NOUT, FMT = 9999 )
     $               'NB      : The size of the square blocks the'//
     $               ' matrix A is split into.'
         WRITE( NOUT, FMT = 9999 )
     $               'NRHS    : The total number of RHS to solve'//
     $               ' for.'
         WRITE( NOUT, FMT = 9999 )
     $               'NBRHS   : The number of RHS to be put on '//
     $               'a column of processes before going'
         WRITE( NOUT, FMT = 9999 )
     $               '          on to the next column of processes.'
         WRITE( NOUT, FMT = 9999 )
     $               'P       : The number of process rows.'
         WRITE( NOUT, FMT = 9999 )
     $               'Q       : The number of process columns.'
         WRITE( NOUT, FMT = 9999 )
     $               'THRESH  : If a residual value is less than'//
     $               ' THRESH, CHECK is flagged as PASSED'
         WRITE( NOUT, FMT = 9999 )
     $               'LU time : Time in seconds to factor the'//
     $               ' matrix'
         WRITE( NOUT, FMT = 9999 )
     $               'Sol Time: Time in seconds to solve the'//
     $               ' system.'
         WRITE( NOUT, FMT = 9999 )
     $               'MFLOPS  : Rate of execution for factor '//
     $               'and solve.'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'The following parameter values will be used:'
         WRITE( NOUT, FMT = 9996 )
     $               'M    ', ( MVAL(I), I = 1, MIN(NMAT, 10) )
         IF( NMAT.GT.10 )
     $      WRITE( NOUT, FMT = 9997 ) ( MVAL(I), I = 11, NMAT )
         WRITE( NOUT, FMT = 9996 )
     $               'N    ', ( NVAL(I), I = 1, MIN(NMAT, 10) )
         IF( NMAT.GT.10 )
     $      WRITE( NOUT, FMT = 9997 ) ( NVAL(I), I = 11, NMAT )
         WRITE( NOUT, FMT = 9996 )
     $               'NB   ', ( NBVAL(I), I = 1, MIN(NNB, 10) )
         IF( NNB.GT.10 )
     $      WRITE( NOUT, FMT = 9997 ) ( NBVAL(I), I = 11, NNB )
         WRITE( NOUT, FMT = 9996 )
     $               'NRHS ', ( NRVAL(I), I = 1, MIN(NNR, 10) )
         IF( NNR.GT.10 )
     $      WRITE( NOUT, FMT = 9997 ) ( NRVAL(I), I = 11, NNR )
         WRITE( NOUT, FMT = 9996 )
     $               'NBRHS', ( NBRVAL(I), I = 1, MIN(NNBR, 10) )
         IF( NNBR.GT.10 )
     $      WRITE( NOUT, FMT = 9997 ) ( NBRVAL(I), I = 11, NNBR )
         WRITE( NOUT, FMT = 9996 )
     $               'P    ', ( PVAL(I), I = 1, MIN(NGRIDS, 10) )
         IF( NGRIDS.GT.10 )
     $      WRITE( NOUT, FMT = 9997) ( PVAL(I), I = 11, NGRIDS )
         WRITE( NOUT, FMT = 9996 )
     $               'Q    ', ( QVAL(I), I = 1, MIN(NGRIDS, 10) )
         IF( NGRIDS.GT.10 )
     $      WRITE( NOUT, FMT = 9997 ) ( QVAL(I), I = 11, NGRIDS )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9995 ) EPS
         WRITE( NOUT, FMT = 9998 ) THRESH
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
         EPS = PSLAMCH( ICTXT, 'eps' )
*
         CALL SGEBR2D( ICTXT, 'All', ' ', 1, 1, THRESH, 1, 0, 0 )
         CALL IGEBR2D( ICTXT, 'All', ' ', 6, 1, WORK, 6, 0, 0 )
         NMAT   = WORK( 1 )
         NNB    = WORK( 2 )
         NNR    = WORK( 3 )
         NNBR   = WORK( 4 )
         NGRIDS = WORK( 5 )
         IF( WORK( 6 ).EQ.1 ) THEN
            EST = .TRUE.
         ELSE
            EST = .FALSE.
         END IF
*
         I = 2*NMAT + NNB + NNR + NNBR + 2*NGRIDS
         CALL IGEBR2D( ICTXT, 'All', ' ', I, 1, WORK, I, 0, 0 )
         I = 1
         CALL ICOPY( NMAT, WORK( I ), 1, MVAL, 1 )
         I = I + NMAT
         CALL ICOPY( NMAT, WORK( I ), 1, NVAL, 1 )
         I = I + NMAT
         CALL ICOPY( NNB, WORK( I ), 1, NBVAL, 1 )
         I = I + NNB
         CALL ICOPY( NNR, WORK( I ), 1, NRVAL, 1 )
         I = I + NNR
         CALL ICOPY( NNBR, WORK( I ), 1, NBRVAL, 1 )
         I = I + NNBR
         CALL ICOPY( NGRIDS, WORK( I ), 1, PVAL, 1 )
         I = I + NGRIDS
         CALL ICOPY( NGRIDS, WORK( I ), 1, QVAL, 1 )
*
      END IF
*
      CALL BLACS_GRIDEXIT( ICTXT )
*
      RETURN
*
   20 WRITE( NOUT, FMT = 9993 )
      CLOSE( NIN )
      IF( NOUT.NE.6 .AND. NOUT.NE.0 )
     $   CLOSE( NOUT )
      CALL BLACS_ABORT( ICTXT, 1 )
*
      STOP
*
 9999 FORMAT( A )
 9998 FORMAT( 'Routines pass computational tests if scaled residual ',
     $        'is less than ', G12.5 )
 9997 FORMAT( '                 ', 10I6 )
 9996 FORMAT( 2X, A5, '   :        ', 10I6 )
 9995 FORMAT( 'Relative machine precision (eps) is taken to be ',
     $        E18.6 )
 9994 FORMAT( ' Number of values of ',5A, ' is less than 1 or greater ',
     $        'than ', I2 )
 9993 FORMAT( ' Illegal input in file ',40A,'.  Aborting run.' )
*
*     End of PSLUINFO
*
      END
