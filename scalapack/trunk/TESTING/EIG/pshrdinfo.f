      SUBROUTINE PSHRDINFO( SUMMRY, NOUT, NMAT, NVAL, NVLO, NVHI,
     $                      LDNVAL, NNB, NBVAL, LDNBVAL, NGRIDS, PVAL,
     $                      LDPVAL, QVAL, LDQVAL, THRESH, WORK, IAM,
     $                      NPROCS )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER               IAM, LDNBVAL, LDNVAL, LDPVAL, LDQVAL,
     $                      NGRIDS, NMAT, NNB, NOUT, NPROCS
      REAL                  THRESH
*     ..
*     .. Array Arguments ..
      CHARACTER*( * )       SUMMRY
      INTEGER               NBVAL( LDNBVAL ), NVAL( LDNVAL ),
     $                      NVHI( LDNVAL ), NVLO( LDNVAL ),
     $                      PVAL( LDPVAL ), QVAL( LDQVAL ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSHRDINFO get the needed startup information for the Hessenberg
*  reduction tests and transmits it to all processes.
*
*  Arguments
*  =========
*
*  SUMMRY   (global output) CHARACTER*(*)
*           Name of output (summary) file (if any). Only defined for
*           process 0.
*
*  NOUT     (global output) INTEGER
*           The unit number for output file. NOUT = 6, output to screen,
*           NOUT = 0, output to stderr.  Only defined for process 0.
*
*  NMAT     (global output) INTEGER
*           The number of different values that can be used for
*           N, IHI & ILO.
*
*  NVAL     (global output) INTEGER array, dimension (LDNVAL)
*           The values of N (number of rows & columns in matrix).
*
*  NVLO     (global output) INTEGER array, dimension (LDNVAL)
*           The values of ILO.
*
*  NVHI     (global output) INTEGER array, dimension (LDNVAL)
*           The values of IHI.
*
*  LDNVAL   (global input) INTEGER
*           The maximum number of different values that can be used for
*           N, ILO and IHI. LDNVAL >=  NMAT.
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
*           = 0 : Perform no error checking
*           > 0 : report all residuals greater than THRESH.
*
*  WORK     (local workspace) INTEGER array, dimension >=
*           3*LDNVAL+LDNBVAL+2*LDPVAL.  Used to pack all input arrays
*           in order to send info in one message.
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
         OPEN( UNIT = NIN, FILE = 'HRD.dat', STATUS = 'OLD' )
         READ( NIN, FMT = * )SUMMRY
         SUMMRY = ' '
*
*        Read in user-supplied info about machine type, compiler, etc.
*
         READ( NIN, FMT = * ) USRINFO
*
*        Read name and unit number for summary output file
*
         READ( NIN, FMT = * ) SUMMRY
         READ( NIN, FMT = * ) NOUT
         IF( NOUT.NE.0 .AND. NOUT.NE.6 )
     $      OPEN( UNIT = NOUT, FILE = SUMMRY, STATUS = 'UNKNOWN' )
*
*        Read and check the parameter values for the tests.
*
*        Get number of matrices
*
         READ( NIN, FMT = * ) NMAT
         IF( NMAT.LT.1. .OR. NMAT.GT.LDNVAL ) THEN
            WRITE( NOUT, FMT = 9997 ) 'N', LDNVAL
            GO TO 20
         END IF
*
*        Get values of N, ILO, IHI
*
         READ( NIN, FMT = * ) ( NVAL( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( NVLO( I ), I = 1, NMAT )
         READ( NIN, FMT = * ) ( NVHI( I ), I = 1, NMAT )
*
*        Get values of NB
*
         READ( NIN, FMT = * ) NNB
         IF( NNB.LT.1 .OR. NNB.GT.LDNBVAL ) THEN
            WRITE( NOUT, FMT = 9997 ) 'NB', LDNBVAL
            GO TO 20
         END IF
         READ( NIN, FMT = * ) ( NBVAL( I ), I = 1, NNB )
*
*        Get number of grids
*
         READ( NIN, FMT = * ) NGRIDS
         IF( NGRIDS.LT.1 .OR. NGRIDS.GT.LDPVAL ) THEN
            WRITE( NOUT, FMT = 9997 ) 'Grids', LDPVAL
            GO TO 20
         ELSE IF( NGRIDS.GT.LDQVAL ) THEN
            WRITE( NOUT, FMT = 9997 ) 'Grids', LDQVAL
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
         WORK( 3 ) = NGRIDS
         CALL IGEBS2D( ICTXT, 'All', ' ', 1, 3, WORK, 1 )
*
         I = 1
         CALL ICOPY( NMAT, NVAL, 1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT, NVLO, 1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NMAT, NVHI, 1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NNB, NBVAL, 1, WORK( I ), 1 )
         I = I + NNB
         CALL ICOPY( NGRIDS, PVAL, 1, WORK( I ), 1 )
         I = I + NGRIDS
         CALL ICOPY( NGRIDS, QVAL, 1, WORK( I ), 1 )
         I = I + NGRIDS -1
         CALL IGEBS2D( ICTXT, 'All', ' ', 1, I, WORK, 1 )
*
*        regurgitate input
*
         WRITE( NOUT, FMT = 9999 )
     $      'ScaLAPACK Reduction routine to Hessenberg form.'
         WRITE( NOUT, FMT = 9999 ) USRINFO
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'Tests of the parallel '//
     $               'real single precision Hessenberg '
         WRITE( NOUT, FMT = 9999 ) 'reduction routines.'
         WRITE( NOUT, FMT = 9999 )
     $               'The following scaled residual '//
     $               'checks will be computed:'
         WRITE( NOUT, FMT = 9999 )
     $               ' ||A - Q H Q''|| / (||A|| * eps * N)'
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
         WRITE( NOUT, FMT = 9999 )
     $               'N        : The number of rows and columns '//
     $               'of the matrix A.'
         WRITE( NOUT, FMT = 9999 )
     $               'NB       : The size of the square blocks'//
     $               ' the matrix A is split into.'
         WRITE( NOUT, FMT = 9999 )
     $      '         on to the next column of processes.'
         WRITE( NOUT, FMT = 9999 )
     $               'P        : The number of process rows.'
         WRITE( NOUT, FMT = 9999 )
     $               'Q        : The number of process columns.'
         WRITE( NOUT, FMT = 9999 )
     $               'HRD time : Time in seconds to compute HRD '
         WRITE( NOUT, FMT = 9999 )
     $               'MFLOPS   : Rate of execution for HRD ' //
     $               'reduction.'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'The following parameter values will be used:'
         WRITE( NOUT, FMT = 9995 )
     $               'N    ', ( NVAL( I ), I = 1, MIN( NMAT, 10 ) )
         IF( NMAT.GT.10 )
     $      WRITE( NOUT, FMT = 9994 ) ( NVAL( I ), I = 11, NMAT )
         WRITE( NOUT, FMT = 9995 )
     $               'ILO    ', ( NVLO( I ), I = 1, MIN( NMAT, 10 ) )
         IF( NMAT.GT.10 )
     $      WRITE( NOUT, FMT = 9994 ) ( NVLO( I ), I = 11, NMAT )
         WRITE( NOUT, FMT = 9995 )
     $               'IHI    ', ( NVHI( I ), I = 1, MIN( NMAT, 10 ) )
         IF( NMAT.GT.10 )
     $      WRITE( NOUT, FMT = 9994 ) ( NVHI( I ), I = 11, NMAT )
         WRITE( NOUT, FMT = 9995 )
     $               'NB   ', ( NBVAL( I ), I = 1, MIN( NNB, 10 ) )
         IF( NNB.GT.10 )
     $      WRITE( NOUT, FMT = 9994 ) ( NBVAL( I ), I = 11, NNB )
         WRITE( NOUT, FMT = 9995 )
     $               'P    ', ( PVAL( I ), I = 1, MIN( NGRIDS, 10 ) )
         IF( NGRIDS.GT.10 )
     $      WRITE( NOUT, FMT = 9994 ) ( PVAL( I ), I = 11, NGRIDS )
         WRITE( NOUT, FMT = 9995 )
     $               'Q    ', ( QVAL( I ), I = 1, MIN( NGRIDS, 10 ) )
         IF( NGRIDS.GT.10 )
     $      WRITE( NOUT, FMT = 9994 ) ( QVAL( I ), I = 11, NGRIDS )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9996 ) EPS
         WRITE( NOUT, FMT = 9993 ) THRESH
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
         EPS = PSLAMCH( ICTXT, 'eps' )
*
         CALL SGEBR2D( ICTXT, 'All', ' ', 1, 1, THRESH, 1, 0, 0 )
         CALL IGEBR2D( ICTXT, 'All', ' ', 1, 3, WORK, 1, 0, 0 )
         NMAT   = WORK( 1 )
         NNB    = WORK( 2 )
         NGRIDS = WORK( 3 )
*
         I = 3*NMAT + NNB + 2*NGRIDS
         CALL IGEBR2D( ICTXT, 'All', ' ', 1, I, WORK, 1, 0, 0 )
*
         I = 1
         CALL ICOPY( NMAT, WORK( I ), 1, NVAL, 1 )
         I = I + NMAT
         CALL ICOPY( NMAT, WORK( I ), 1, NVLO, 1 )
         I = I + NMAT
         CALL ICOPY( NMAT, WORK( I ), 1, NVHI, 1 )
         I = I + NMAT
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
   20 CONTINUE
      WRITE( NOUT, FMT = 9998 )
      CLOSE( NIN )
      IF( NOUT.NE.6 .AND. NOUT.NE.0 )
     $   CLOSE( NOUT )
      CALL BLACS_ABORT( ICTXT, 1 )
*
      STOP
*
 9999 FORMAT( A )
 9998 FORMAT( ' ILLEGAL INPUT IN FILE ', 40A, '.  ABORTING RUN.' )
 9997 FORMAT( ' NUMBER OF VALUES OF ', 5A,
     $      ' IS LESS THAN 1 OR GREATER ', 'THAN ', I2 )
 9996 FORMAT( 'Relative machine precision (eps) is taken to be ',
     $      E18.6 )
 9995 FORMAT( 2X, A5, ':        ', 10I6 )
 9994 FORMAT( '              ', 10I6 )
 9993 FORMAT( 'Routines pass computational tests if scaled residual is',
     $      ' less than ', G14.7 )
*
*     End of PSHRDINFO
*
      END
