      SUBROUTINE PZQRINFO( SUMMRY, NOUT, NFACT, FACTOR, LDFACT, NMAT,
     $                     MVAL, LDMVAL, NVAL, LDNVAL, NNB, MBVAL,
     $                     LDMBVAL, NBVAL, LDNBVAL, NGRIDS, PVAL,
     $                     LDPVAL, QVAL, LDQVAL, THRESH, WORK, IAM,
     $                     NPROCS )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IAM, LDFACT, LDMBVAL, LDMVAL, LDNBVAL, LDNVAL,
     $                   LDPVAL, LDQVAL, NFACT, NGRIDS, NMAT, NNB,
     $                   NPROCS, NOUT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      CHARACTER*2        FACTOR( LDFACT )
      CHARACTER*(*)      SUMMRY
      INTEGER            MBVAL( LDMBVAL ), MVAL( LDMVAL ),
     $                   NBVAL( LDNBVAL ), NVAL( LDNVAL ),
     $                   PVAL( LDPVAL ), QVAL( LDQVAL ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZQRINFO gets needed startup information for the QR factoriza-
*  tion routines and transmits it to all processes.
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
*  NFACT    (global output) INTEGER
*           The number of different factorization types to be tested.
*
*  FACTOR   (global output) CHARACTER*2 array of dimension of LDFACT,
*           The factorization types to be tested:
*           if FACTOR(i) = 'QR' then QR factorization,
*           if FACTOR(i) = 'QL' then QL factorization,
*           if FACTOR(i) = 'LQ' then LQ factorization,
*           if FACTOR(i) = 'RQ' then RQ factorization,
*           if FACTOR(i) = 'QP' then QR factorization with column
*                                    pivoting.
*           if FACTOR(i) = 'TZ' then complete unitary factorization.
*
*  LDFACT   (global input) INTEGER
*           The maximum number of different factorization types to be
*           tested.  LDFACT >=  NFACT.
*
*  NMAT     (global output) INTEGER
*           The number of different values that can be used for N.
*
*  MVAL     (global output) INTEGER array of dimension (LDNVAL), the
*           values of M (number of rows in matrix) to run the code
*           with.
*
*  LDMVAL   (global input) INTEGER
*           The maximum number of different values that can be used for
*           M, LDNVAL > =  NMAT.
*
*  NVAL     (global output) INTEGER array of dimension (LDNVAL), the
*           values of N (number of columns in matrix) to run the code
*           with.
*
*  LDNVAL   (global input) INTEGER
*           The maximum number of different values that can be used for
*           N, LDNVAL > =  NMAT.
*
*  NNB      (global output) INTEGER
*           The number of different values that can be used for MB and
*           NB.
*
*  MBVAL    (global output) INTEGER array of dimension (LDMBVAL), the
*           values of MB (row blocksize) to run the code with.
*
*  LDMBVAL  (global input) INTEGER
*           The maximum number of different values that can be used for
*           MB, LDMBVAL >= NNB.
*
*  NBVAL    (global output) INTEGER array of dimension (LDNBVAL), the
*           values of NB (column blocksize) to run the code with.
*
*  LDNBVAL  (global input) INTEGER
*           The maximum number of different values that can be used for
*           NB, LDNBVAL >= NNB.
*
*  NGRIDS   (global output) INTEGER
*           The number of different values that can be used for P & Q.
*
*  PVAL     (global output) INTEGER array of dimension (LDPVAL), the
*           values of P (number of process rows) to run the code with.
*
*  LDPVAL   (global input) INTEGER
*           The maximum number of different values that can be used for
*           P, LDPVAL >= NGRIDS.
*
*  QVAL     (global output) INTEGER array of dimension (LDQVAL), the
*           values of Q (number of process columns) to run the code
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
*  WORK     (local workspace) INTEGER array of dimension >=
*           MAX( 4, LDFACT+LDMVAL+LDNVAL+LDMBVAL+LDNBVAL+LDPVAL+LDQVAL )
*           used to pack all input arrays in order to send info in one
*           message.
*
*  IAM      (local input) INTEGER
*           My process number.
*
*  NPROCS   (global input) INTEGER
*           The total number of processes.
*
*  Note
*  ====
*
*  For packing the information we assumed that the length in bytes of an
*  integer is equal to the length in bytes of a real single precision.
*
*  =====================================================================
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
      INTEGER            I, ICTXT, K
      DOUBLE PRECISION   EPS
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINIT, BLACS_SETUP, ICOPY, IGEBR2D,
     $                   IGEBS2D, SGEBR2D, SGEBS2D
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           LSAMEN, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
         OPEN( NIN, FILE='QR.dat', STATUS='OLD' )
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
*        Get the matrix types to be tested
*
         READ( NIN, FMT = * ) NFACT
         IF( NFACT.LT.1 .OR. NFACT.GT.LDFACT ) THEN
            WRITE( NOUT, FMT = 9994 ) 'nb of factorization', LDFACT
            GO TO 40
         END IF
         READ( NIN, FMT = * ) ( FACTOR( I ), I = 1, NFACT )
*
*        Get number of matrices and their dimensions
*
         READ( NIN, FMT = * ) NMAT
         IF( NMAT.LT.1 .OR. NMAT.GT.LDNVAL ) THEN
            WRITE( NOUT, FMT = 9994 ) 'N', LDNVAL
            GO TO 40
         ELSE IF( NMAT.GT.LDMVAL ) THEN
            WRITE( NOUT, FMT = 9994 ) 'M', LDMVAL
            GO TO 40
         END IF
         READ( NIN, FMT = * ) ( MVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( NVAL( I ), I = 1, NMAT )
*
*        Get values of NB
*
         READ( NIN, FMT = * ) NNB
         IF( NNB.LT.1 .OR. NNB.GT.LDMBVAL ) THEN
            WRITE( NOUT, FMT = 9994 ) 'MB', LDMBVAL
            GO TO 40
         ELSE IF( NNB.GT.LDNBVAL ) THEN
            WRITE( NOUT, FMT = 9994 ) 'NB', LDNBVAL
            GO TO 40
         END IF
         READ( NIN, FMT = * ) ( MBVAL( I ), I = 1, NNB )
         READ( NIN, FMT = * ) ( NBVAL( I ), I = 1, NNB )
*
*        Get number of grids
*
         READ( NIN, FMT = * ) NGRIDS
         IF( NGRIDS.LT.1 .OR. NGRIDS.GT.LDPVAL ) THEN
            WRITE( NOUT, FMT = 9994 ) 'Grids', LDPVAL
            GO TO 40
         ELSE IF( NGRIDS.GT.LDQVAL ) THEN
            WRITE( NOUT, FMT = 9994 ) 'Grids', LDQVAL
            GO TO 40
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
               NPROCS = MAX( NPROCS, PVAL( I ) * QVAL( I ) )
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
         EPS = PDLAMCH( ICTXT, 'eps' )
*
*        Pack information arrays and broadcast
*
         CALL SGEBS2D( ICTXT, 'All', ' ', 1, 1, THRESH, 1 )
         WORK( 1 ) = NMAT
         WORK( 2 ) = NNB
         WORK( 3 ) = NGRIDS
         WORK( 4 ) = NFACT
         CALL IGEBS2D( ICTXT, 'All', ' ', 4, 1, WORK, 4 )
*
         I = 1
         DO 20 K = 1, NFACT
            IF( LSAMEN( 2, FACTOR( K ), 'QR' ) ) THEN
               WORK( I ) = 1
               I = I + 1
            ELSE IF( LSAMEN( 2, FACTOR( K ), 'QL' ) ) THEN
               WORK( I ) = 2
               I = I + 1
            ELSE IF( LSAMEN( 2, FACTOR( K ), 'LQ' ) ) THEN
               WORK( I ) = 3
               I = I + 1
            ELSE IF( LSAMEN( 2, FACTOR( K ), 'RQ' ) ) THEN
               WORK( I ) = 4
               I = I + 1
            ELSE IF( LSAMEN( 2, FACTOR( K ), 'QP' ) ) THEN
               WORK( I ) = 5
               I = I + 1
            ELSE IF( LSAMEN( 2, FACTOR( K ), 'TZ' ) ) THEN
               WORK( I ) = 6
               I = I + 1
            END IF
   20    CONTINUE
*
         CALL ICOPY( NMAT, MVAL, 1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT, NVAL, 1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NNB, MBVAL, 1, WORK( I ), 1 )
         I = I + NNB
         CALL ICOPY( NNB, NBVAL, 1, WORK( I ), 1 )
         I = I + NNB
         CALL ICOPY( NGRIDS, PVAL, 1, WORK( I ), 1 )
         I = I + NGRIDS
         CALL ICOPY( NGRIDS, QVAL, 1, WORK( I ), 1 )
         I = I + NGRIDS - 1
         CALL IGEBS2D( ICTXT, 'All', ' ', I, 1, WORK, I )
*
*        regurgitate input
*
         WRITE( NOUT, FMT = 9999 )
     $               'ScaLAPACK QR factorizations routines.'
         WRITE( NOUT, FMT = 9999 ) USRINFO
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'Tests of the parallel '//
     $               'complex single precision QR factorizations '//
     $               'routines.'
         WRITE( NOUT, FMT = 9999 )
     $               'The following scaled residual '//
     $               'checks will be computed:'
         WRITE( NOUT, FMT = 9999 )
     $               ' || A - QR  || / (|| A || * eps * N) and/or'
         WRITE( NOUT, FMT = 9999 )
     $               ' || A - QL  || / (|| A || * eps * N) and/or'
         WRITE( NOUT, FMT = 9999 )
     $               ' || A - LQ  || / (|| A || * eps * N) and/or'
         WRITE( NOUT, FMT = 9999 )
     $               ' || A - RQ  || / (|| A || * eps * N) and/or'
         WRITE( NOUT, FMT = 9999 )
     $               ' || A - QRP || / (|| A || * eps * N) and/or'
         WRITE( NOUT, FMT = 9999 )
     $               ' || A - TZ  || / (|| A || * eps * N)'
         WRITE( NOUT, FMT = 9999 )
     $               'The matrix A is randomly '//
     $               'generated for each test.'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'An explanation of the input/output '//
     $               'parameters follows:'
         WRITE( NOUT, FMT = 9999 )
     $               'TIME     : Indicates whether WALL or '//
     $               'CPU time was used.'
*
         WRITE( NOUT, FMT = 9999 )
     $               'M        : The number of rows in the '//
     $               'matrix A.'
         WRITE( NOUT, FMT = 9999 )
     $               'N        : The number of columns in the '//
     $               'matrix A.'
         WRITE( NOUT, FMT = 9999 )
     $               'MB       : The row blocksize of the blocks'//
     $               ' the matrix A is split into.'
         WRITE( NOUT, FMT = 9999 )
     $               'NB       : The column blocksize of the blocks'//
     $               ' the matrix A is split into.'
         WRITE( NOUT, FMT = 9999 )
     $               'P        : The number of process rows.'
         WRITE( NOUT, FMT = 9999 )
     $               'Q        : The number of process columns.'
         WRITE( NOUT, FMT = 9999 )
     $               'THRESH   : If a residual value is less than'//
     $               ' THRESH, CHECK is flagged as PASSED'
         WRITE( NOUT, FMT = 9999 )
         WRITE( NOUT, FMT = 9999 )
     $               'Fact Time: Time in seconds to factor the'//
     $               ' matrix.'
         WRITE( NOUT, FMT = 9999 )
     $               'MFLOPS   : Execution rate of the '//
     $               'factorization.'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'The following parameter values will be used:'
         WRITE( NOUT, FMT = 9996 )
     $               'M    ', ( MVAL( I ), I = 1, MIN( NMAT, 10 ) )
         IF( NMAT.GT.10 )
     $      WRITE( NOUT, FMT = 9997 ) ( MVAL( I ), I = 11, NMAT )
         WRITE( NOUT, FMT = 9996 )
     $               'N    ', ( NVAL( I ), I = 1, MIN( NMAT, 10 ) )
         IF( NMAT.GT.10 )
     $      WRITE( NOUT, FMT = 9997 ) ( NVAL( I ), I = 11, NMAT )
         WRITE( NOUT, FMT = 9996 )
     $               'MB   ', ( MBVAL( I ), I = 1, MIN( NNB, 10 ) )
         IF( NNB.GT.10 )
     $      WRITE( NOUT, FMT = 9997 ) ( MBVAL( I ), I = 11, NNB )
         WRITE( NOUT, FMT = 9996 )
     $               'NB   ', ( NBVAL( I ), I = 1, MIN( NNB, 10 ) )
         IF( NNB.GT.10 )
     $      WRITE( NOUT, FMT = 9997 ) ( NBVAL( I ), I = 11, NNB )
         WRITE( NOUT, FMT = 9996 )
     $               'P    ', ( PVAL( I ), I = 1, MIN( NGRIDS, 10 ) )
         IF( NGRIDS.GT.10 )
     $      WRITE( NOUT, FMT = 9997) ( PVAL( I ), I = 11, NGRIDS )
         WRITE( NOUT, FMT = 9996 )
     $               'Q    ', ( QVAL( I ), I = 1, MIN( NGRIDS, 10 ) )
         IF( NGRIDS.GT.10 )
     $      WRITE( NOUT, FMT = 9997 ) ( QVAL( I ), I = 11, NGRIDS )
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
*        all processes have needed startup information
*
         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS )
*
*        Compute machine epsilon
*
         EPS = PDLAMCH( ICTXT, 'eps' )
*
         CALL SGEBR2D( ICTXT, 'All', ' ', 1, 1, THRESH, 1, 0, 0 )
         CALL IGEBR2D( ICTXT, 'All', ' ', 4, 1, WORK, 4, 0, 0 )
         NMAT   = WORK( 1 )
         NNB    = WORK( 2 )
         NGRIDS = WORK( 3 )
         NFACT  = WORK( 4 )
*
         I = NFACT + 2*NMAT + 2*NNB + 2*NGRIDS
         CALL IGEBR2D( ICTXT, 'All', ' ', I, 1, WORK, I, 0, 0 )
*
         DO 30 K = 1, NFACT
            IF( WORK( K ).EQ.1 ) THEN
               FACTOR( K ) = 'QR'
            ELSE IF( WORK( K ).EQ.2 ) THEN
               FACTOR( K ) = 'QL'
            ELSE IF( WORK( K ).EQ.3 ) THEN
               FACTOR( K ) = 'LQ'
            ELSE IF( WORK( K ).EQ.4 ) THEN
               FACTOR( K ) = 'RQ'
            ELSE IF( WORK( K ).EQ.5 ) THEN
               FACTOR( K ) = 'QP'
            ELSE IF( WORK( K ).EQ.6 ) THEN
               FACTOR( K ) = 'TZ'
            END IF
   30    CONTINUE
*
         I = NFACT + 1
         CALL ICOPY( NMAT, WORK( I ), 1, MVAL, 1 )
         I = I + NMAT
         CALL ICOPY( NMAT, WORK( I ), 1, NVAL, 1 )
         I = I + NMAT
         CALL ICOPY( NNB, WORK( I ), 1, MBVAL, 1 )
         I = I + NNB
         CALL ICOPY( NNB, WORK( I ), 1, NBVAL, 1 )
         I = I + NNB
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
   40 WRITE( NOUT, FMT = 9993 )
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
*     End of PZQRINFO
*
      END
