      PROGRAM PSSVDDRIVER
*
*  -- ScaLAPACK testing driver (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997      
*
*  Purpose
*  ========
*
*  Parallel Real singular value decomposition test driver.
*
*  INPUT:
*  =====
*  This routine tests PDGESVD, the parallel singular value
*  decomposition  solver. We would like to cover possible combinations
*  of: matrix size, process configuration (nprow and npcol), block 
*  size (nb), matrix type, and workspace available.
*
*  Current format of the input file SVD.dat lists the following:
*        device out
*        Threshold
*        number of matrices
*        number of rows for every matrix
*        number of columns for every matrix
*        number of process configurations (P, Q, NB)
*        values of P (NPROW) for every configuration
*        values of Q (NPCOL) for every configuration
*        values of NB for every configuration.
*  Here threshold is an integer constant with a value between 1 and
*  100, which meaning is explained in comments to PSSVDTST.
*
*  WHAT IT DOES:
*  ============
*  PSVDDRIVER checks floating-point arithmetic and parameters
*  provided by the user in initialization file SVD.dat. It reads and
*  broadcasts to all process parameters required to run actual testing 
*  code PSVDTST. In case all tests are successful it tells you so. For
*  the actual "meat" of the tests see comments to PSVDTST.
*
*=======================================================================
*
*     .. Local Scalars ..
      CHARACTER*80       SUMMARY
      INTEGER            CONTEXT, ERR, I, IAM, J, K, LWORK, MAXNODES,
     $                   NMATSIZES, NOUT, NPCONFIGS, NPROCS
      REAL               THRESH
*     ..
*     .. Parameters ..
      INTEGER            MAXSETSIZE, NIN, DBLSIZ, TOTMEM, MEMSIZ
      PARAMETER          ( MAXSETSIZE = 50, NIN = 11, DBLSIZ = 8,
     $                   TOTMEM = 2000000, MEMSIZ = TOTMEM / DBLSIZ )
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 ), MM( MAXSETSIZE ),
     $                   NBS( MAXSETSIZE ), NN( MAXSETSIZE ),
     $                   NPCOLS( MAXSETSIZE ), NPROWS( MAXSETSIZE ),
     $                   RESULT( 9 )
      REAL               WORK( MEMSIZ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINIT, BLACS_PINFO, BLACS_SETUP,
     $                   SGEBR2D, SGEBS2D, IGEBR2D, IGEBS2D, PSSVDTST
*     ..
*     .. Executable Statements ..
*
*     Get starting information.
*
      CALL BLACS_PINFO( IAM, NPROCS )
*
*     Open file and skip data header; read output device.
*
      IF( IAM.EQ.0 ) THEN
         OPEN( UNIT = NIN, FILE = 'SVD.dat', STATUS = 'OLD' )
         READ( NIN, FMT = * )SUMMARY
         READ( NIN, FMT = * )NOUT
         READ( NIN, FMT = * )MAXNODES
      END IF
*
      IF( NPROCS.LT.1 ) THEN
         CALL BLACS_SETUP( IAM, MAXNODES )
         NPROCS = MAXNODES
      END IF
*
      CALL BLACS_GET( -1, 0, CONTEXT )
      CALL BLACS_GRIDINIT( CONTEXT, 'R', 1, NPROCS )
*
*     Initialize variables, arrays, and grids.
*
      ERR = 0
      NMATSIZES = 0
      NPCONFIGS = 0
      LWORK = MEMSIZ
      ISEED( 1 ) = 139
      ISEED( 2 ) = 1139
      ISEED( 3 ) = 2139
      ISEED( 4 ) = 3139
*
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9992 )
         WRITE( NOUT, FMT = 9991 )
         WRITE( NOUT, FMT = 9990 )
         WRITE( NOUT, FMT = 9989 )
         WRITE( NOUT, FMT = 9988 )
         WRITE( NOUT, FMT = 9987 )
         WRITE( NOUT, FMT = 9986 )
         WRITE( NOUT, FMT = 9985 )
         WRITE( NOUT, FMT = 9984 )
         WRITE( NOUT, FMT = 9983 )
         WRITE( NOUT, FMT = 9982 )
         WRITE( NOUT, FMT = 9981 )
         WRITE( NOUT, FMT = 9980 )
         WRITE( NOUT, FMT = 9979 )
         WRITE( NOUT, FMT = 9978 )
         WRITE( NOUT, FMT = 9977 )
         WRITE( NOUT, FMT = 9976 )
         WRITE( NOUT, FMT = 9975 )
         WRITE( NOUT, FMT = 9974 )
         WRITE( NOUT, FMT = 9973 )
         WRITE( NOUT, FMT = 9972 )
         WRITE( NOUT, FMT = 9971 )
         WRITE( NOUT, FMT = 9970 )
         WRITE( NOUT, FMT = 9969 )
         WRITE( NOUT, FMT = 9968 )
         WRITE( NOUT, FMT = 9967 )
         WRITE( NOUT, FMT = 9966 )
         WRITE( NOUT, FMT = 9965 )
      END IF
*
*     Process 0 reads values in input file and broadcasts them to 
*     all other processes.
*
   10 CONTINUE
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )SUMMARY
         READ( NIN, FMT = * )SUMMARY
         READ( NIN, FMT = * )THRESH
         WRITE( NOUT, FMT = 9965 )SUMMARY
         CALL SGEBS2D( CONTEXT, 'All', ' ', 1, 1, THRESH, 1 )
      ELSE
         CALL SGEBR2D( CONTEXT, 'All', ' ', 1, 1, THRESH, 1, 0, 0 )
      END IF
      IF( THRESH.EQ.-1 ) THEN
         GO TO 80
      END IF
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )NMATSIZES
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, 1, NMATSIZES, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, 1, NMATSIZES, 1, 0, 0 )
      END IF
*     Deal with error
      IF( NMATSIZES.LT.1 .OR. NMATSIZES.GT.MAXSETSIZE ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9999 )'Matrix size', NMATSIZES, 1,
     $         MAXSETSIZE
         END IF
         ERR = -1
         GO TO 80
      END IF
*
*     Read array of MATSIZES.
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )( MM( I ), I = 1, NMATSIZES )
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, NMATSIZES, MM, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, NMATSIZES, MM, 1, 0, 0 )
      END IF
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )( NN( I ), I = 1, NMATSIZES )
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, NMATSIZES, NN, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, NMATSIZES, NN, 1, 0, 0 )
      END IF
*
*     Read and broadcast NPCONFIGS.
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )NPCONFIGS
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, 1, NPCONFIGS, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, 1, NPCONFIGS, 1, 0, 0 )
      END IF
*     Deal with error
      IF( NPCONFIGS.LT.1 .OR. NPCONFIGS.GT.MAXSETSIZE ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9999 )'# proc configs', NPCONFIGS, 1,
     $         MAXSETSIZE
         END IF
         ERR = -1
         GO TO 80
      END IF
*
*     Read and broadcast array of NPROWS.
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )( NPROWS( I ), I = 1, NPCONFIGS )
*
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, NPCONFIGS, NPROWS, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, NPCONFIGS, NPROWS, 1, 0,
     $                 0 )
      END IF
*     Deal with error
      DO 20 I = 1, NPCONFIGS
         IF( NPROWS( I ).LE.0 )
     $      ERR = -1
   20 CONTINUE
      IF( ERR.EQ.-1 ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9997 )' NPROW'
         END IF
         GO TO 80
      END IF
*
*     Read and broadcast array of NPCOLS.
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )( NPCOLS( I ), I = 1, NPCONFIGS )
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, NPCONFIGS, NPCOLS, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, NPCONFIGS, NPCOLS, 1, 0,
     $                 0 )
      END IF
*
*     Deal with error.
*
      DO 30 I = 1, NPCONFIGS
         IF( NPCOLS( I ).LE.0 )
     $      ERR = -1
   30 CONTINUE
      IF( ERR.EQ.-1 ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9997 )' NPCOL'
         END IF
         GO TO 80
      END IF
*
*     Read and broadcast array of NBs.
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )( NBS( I ), I = 1, NPCONFIGS )
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, NPCONFIGS, NBS, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, NPCONFIGS, NBS, 1, 0, 0 )
      END IF
*     Deal with error
      DO 40 I = 1, NPCONFIGS
         IF( NBS( I ).LE.0 )
     $      ERR = -1
   40 CONTINUE
      IF( ERR.EQ.-1 ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9997 )' NB'
         END IF
         GO TO 80
      END IF
*
      DO 70 J = 1, NMATSIZES
         DO 60 I = 1, NPCONFIGS
*
            DO 50 K = 1, 9
               RESULT( K ) = 0
   50       CONTINUE
            CALL PSSVDTST( MM( J ), NN( J ), NPROWS( I ), NPCOLS( I ),
     $                     NBS( I ), ISEED, THRESH, WORK, RESULT, LWORK,
     $                     NOUT )
*
   60    CONTINUE
   70 CONTINUE
*
      GO TO 10
*
   80 CONTINUE
      IF( IAM.EQ.0 ) THEN
         CLOSE ( NIN )
         CLOSE ( NOUT )
      END IF
*
      CALL BLACS_GRIDEXIT( CONTEXT )
*
      CALL BLACS_EXIT( 0 )
      STOP
*
*     End of PSSVDDRIVER
*
 9999 FORMAT( A20, ' is:', I5, ' must be between:', I5, ' and', I5 )
 9998 FORMAT( A )
 9997 FORMAT( A20, ' must be positive' )
 9996 FORMAT( A )
 9995 FORMAT( 'M = ', I5, ' N = ', I5, ' NPOW = ', I5, 'NPCOL = ', I5,
     $      ' NB = ', I5 )
*
 9994 FORMAT( 'Test #', I5, 'for this configuration has failed' )
 9993 FORMAT( 'All test passed for this configuration' )
 9992 FORMAT( ' ' )
 9991 FORMAT( 'Running tests of the parallel singular value ',
     $      'decomposition  routine:  PSGESVD' )
 9990 FORMAT( 'The following scaled residual checks will be',
     $      'computed:' )
 9989 FORMAT( ' || A - U*diag(S)*VT ||/( ||A||*max(M,N)*ulp )' )
 9988 FORMAT( ' || I - UT*U ||/( M*ulp )' )
 9987 FORMAT( ' || I - VT*V ||/( N*ulp )' )
 9986 FORMAT( ' ' )
 9985 FORMAT( 'An explanation of the input/output parameters',
     $      ' follows:' )
 9984 FORMAT( 'RESULT   : passed; or an indication of which',
     $      ' jobtype test failed' )
 9983 FORMAT( 'M        : The number of rows of the matrix A.' )
 9982 FORMAT( 'N        : The number of columns of the matrix A.' )
 9981 FORMAT( 'P        : The number of process rows.' )
 9980 FORMAT( 'Q        : The number of process columns.' )
 9979 FORMAT( 'NB       : The size of the square blocks the',
     $      ' matrix A is split into.' )
 9978 FORMAT( 'THRESH   : If a residual value is less than ',
     $      ' THRESH, RESULT is flagged as PASSED.' )
 9977 FORMAT( 'MTYPE    : matrix type (see pssvdtst.f).' )
 9976 FORMAT( 'CHK      : || A - U*diag(S)*VT ||/( ||A||',
     $      '*max(M,N)*ulp )' )
 9975 FORMAT( 'MTM      : maximum of two values:',/,
     $      '           || I - UT*U ||/( M*ulp ) and',
     $      '  || I - VT*V ||/( N*ulp )' )
 9974 FORMAT( 'DELTA    : maximum of three values:',/,
     $      '           || U - UC ||/( M*ulp*THRESH ),' )
 9973 FORMAT( '           || VT - VTC ||/( N*ulp*THRESH ), and' )
 9972 FORMAT( '           || S - SC || / ( SIZE*ulp*|S|*THRESH ), ' )
 9971 FORMAT( '           where UC, VTC, SC are singular vectors ',
     $      'and values' )
 9970 FORMAT( '           for JOBTYPE.NE.1 (see pdsvdcmp.f) ' )
 9969 FORMAT( 'HET      : P if heterogeneity was detected by PDGESVD' )
 9968 FORMAT( '           T if detected by the PDSVSTST, N if',
     $      ' undetected' )
 9967 FORMAT( ' ' )
 9966 FORMAT( 'RESULT      WALL       CPU     M     N   P   Q',
     $        '   NB MTYPE   CHK   MTM DELTA  HET' )
 9965 FORMAT( A )
      END
