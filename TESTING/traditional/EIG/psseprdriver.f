      PROGRAM PSSEPRDRIVER
*
*     Parallel REAL             symmetric eigenproblem test driver for PSSYEVR
*
      IMPLICIT NONE
*
*     The user should modify TOTMEM to indicate the maximum amount of
*     memory in bytes her system has.  Remember to leave room in memory
*     for operating system, the BLACS buffer, etc.  REALSZ
*     indicates the length in bytes on the given platform for a number,
*     real for SINGLE/DOUBLE PRECISION, and complex for COMPLEX/COMPLEX*16.
*     For example, on a standard system, the length of a
*     REAL is 4, and an integer takes up 4 bytes. Some playing around
*     to discover what the maximum value you can set MEMSIZ to may be
*     required.
*     All arrays used by factorization and solve are allocated out of
*     big array called MEM.
*
*     TESTS PERFORMED
*     ===============
*
*     This routine performs tests for combinations of:  matrix size, process 
*     configuration (nprow and npcol), block size (nb), 
*     matrix type, range of eigenvalue (all, by value, by index), 
*     and upper vs. lower storage.
*
*     It returns an error message when heterogeneity is detected.
*
*     The input file allows multiple requests where each one is 
*     of the following sets:
*       matrix sizes:                     n
*       process configuration triples:  nprow, npcol, nb
*       matrix types:
*       eigenvalue requests:              all, by value, by position
*       storage (upper vs. lower):        uplo
*
*     TERMS:
*       Request - means a set of tests, which is the cross product of
*       a set of specifications from the input file.
*       Test - one element in the cross product, i.e. a specific input
*       size and type, process configuration, etc.
*
*     .. Parameters ..
*
      INTEGER            TOTMEM, REALSZ, NIN
      PARAMETER          ( TOTMEM = 100000000, REALSZ = 4, NIN = 11 )
      INTEGER            MEMSIZ
      PARAMETER          ( MEMSIZ = TOTMEM / REALSZ )
*     ..
*     .. Local Scalars ..
      CHARACTER          HETERO
      CHARACTER*80       SUMMRY, USRINFO
      INTEGER            CONTEXT, IAM, INFO, ISIEEE, MAXNODES, NNOCHECK,
     $                   NOUT, NPASSED, NPROCS, NSKIPPED, NTESTS
*     ..
*     .. Local Arrays ..
*
      INTEGER            ISEED( 4 )
      REAL               MEM( MEMSIZ )
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
*
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT, 
     $                   BLACS_GRIDINIT, BLACS_PINFO, BLACS_SETUP, 
     $                   IGAMN2D, PSLACHKIEEE, PSLASNBT, PSSEPRREQ 
*     ..
*     .. Executable Statements ..
*
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
*
*
      IF( IAM.EQ.0 ) THEN
*
*        Open file and skip data file header
*
         OPEN( UNIT = NIN, FILE = 'SEPR.dat', STATUS = 'OLD' )
         READ( NIN, FMT = * )SUMMRY
         SUMMRY = ' '
*
*        Read in user-supplied info about machine type, compiler, etc.
*
         READ( NIN, FMT = 9999 )USRINFO
*
*        Read name and unit number for summary output file
*
         READ( NIN, FMT = * )SUMMRY
         READ( NIN, FMT = * )NOUT
         IF( NOUT.NE.0 .AND. NOUT.NE.6 )
     $      OPEN( NOUT, FILE = SUMMRY, STATUS = 'UNKNOWN' )
         READ( NIN, FMT = * )MAXNODES
         READ( NIN, FMT = * )HETERO
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
      CALL PSLASNBT( ISIEEE )
*
      CALL IGAMN2D( CONTEXT, 'a', ' ', 1, 1, ISIEEE, 1, 1, 1, -1, -1,
     $              0 )
*
      IF( ( ISIEEE.NE.0 ) ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9997 )
            WRITE( NOUT, FMT = 9996 )
            WRITE( NOUT, FMT = 9995 )
         END IF
*
         CALL PSLACHKIEEE( ISIEEE, SLAMCH( 'O' ), SLAMCH( 'U' ) )
*
         CALL IGAMN2D( CONTEXT, 'a', ' ', 1, 1, ISIEEE, 1, 1, 1, -1, -1,
     $                 0 )
*
         IF( ISIEEE.EQ.0 ) THEN
            GO TO 20
         END IF
*
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9986 )
         END IF
*
      END IF
*
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9999 )
     $      'Test ScaLAPACK symmetric eigendecomposition routine.'
         WRITE( NOUT, FMT = 9999 )USRINFO
         WRITE( NOUT, FMT = 9999 )' '
         WRITE( NOUT, FMT = 9999 )'Running tests of the parallel ' //
     $      'symmetric eigenvalue routine:  PSSYEVR.'
         WRITE( NOUT, FMT = 9999 )'The following scaled residual ' //
     $      'checks will be computed:'
         WRITE( NOUT, FMT = 9999 )' ||AQ - QL|| ' //
     $      '/ ((abstol + ||A|| * eps) * N)'
         WRITE( NOUT, FMT = 9999 )' ||Q^T*Q - I|| ' // '/ (N * eps)'
         WRITE( NOUT, FMT = 9999 )
         WRITE( NOUT, FMT = 9999 )'An explanation of the ' //
     $      'input/output parameters follows:'
         WRITE( NOUT, FMT = 9999 )'RESULT   : passed; or ' //
     $      'an indication of which eigen request test failed'
         WRITE( NOUT, FMT = 9999 )
     $      'N        : The number of rows and columns ' //
     $      'of the matrix A.'
         WRITE( NOUT, FMT = 9999 )
     $      'P        : The number of process rows.'
         WRITE( NOUT, FMT = 9999 )
     $      'Q        : The number of process columns.'
         WRITE( NOUT, FMT = 9999 )
     $      'NB       : The size of the square blocks' //
     $      ' the matrix A is split into.'
         WRITE( NOUT, FMT = 9999 )
     $      'THRESH   : If a residual value is less ' //
     $      'than THRESH, RESULT = PASSED.'
         WRITE( NOUT, FMT = 9999 )
     $      'TYP      : matrix type (see PSSEPRTST).'
         WRITE( NOUT, FMT = 9999 )'SUB      : Subtests (Y/N).'
         WRITE( NOUT, FMT = 9999 )'WALL     : Wallclock time.'
         WRITE( NOUT, FMT = 9999 )'CPU      : CPU time.'
         WRITE( NOUT, FMT = 9999 )'CHK      : ||AQ - QL|| ' //
     $      '/ ((abstol + ||A|| * eps) * N)'
         WRITE( NOUT, FMT = 9999 )'QTQ      : ||Q^T*Q - I||/ (N * eps)'
         WRITE( NOUT, FMT = 9999 )
     $      '         : when the adjusted QTQ norm exceeds THRESH',
     $      '           it is printed,'
         WRITE( NOUT, FMT = 9999 )
     $      '           otherwise the true QTQ norm is printed.'
         WRITE( NOUT, FMT = 9999 )
     $      '         : If more than one test is done, CHK and QTQ ' 
         WRITE( NOUT, FMT = 9999 )
     $      '           are the max over all eigentests performed.'
         WRITE( NOUT, FMT = 9999 )
     $      'TEST     : EVR - testing PSSYEVR'
         WRITE( NOUT, FMT = 9999 )' '
      END IF
*
      NTESTS = 0
      NPASSED = 0
      NSKIPPED = 0
      NNOCHECK = 0
*
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9979 )
         WRITE( NOUT, FMT = 9978 )
      END IF
*
   10 CONTINUE
*
      ISEED( 1 ) = 139
      ISEED( 2 ) = 1139
      ISEED( 3 ) = 2139
      ISEED( 4 ) = 3139
*
      CALL PSSEPRREQ( HETERO, NIN, MEM, MEMSIZ, NOUT, ISEED, NTESTS,
     $               NSKIPPED, NNOCHECK, NPASSED, INFO )
      IF( INFO.EQ.0 )
     $   GO TO 10
*
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9985 )NTESTS
         WRITE( NOUT, FMT = 9984 )NPASSED
         WRITE( NOUT, FMT = 9983 )NNOCHECK
         WRITE( NOUT, FMT = 9982 )NSKIPPED
         WRITE( NOUT, FMT = 9981 )NTESTS - NPASSED - NSKIPPED -
     $      NNOCHECK
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9980 )
      END IF
*
*     Uncomment this line on SUN systems to avoid the useless print out
*
c      CALL IEEE_FLAGS( 'clear', 'exception', 'underflow', ' ')
*
   20 CONTINUE
      IF( IAM.EQ.0 ) THEN
         CLOSE ( NIN )
         IF( NOUT.NE.6 .AND. NOUT.NE.0 )
     $      CLOSE ( NOUT )
      END IF
*
      CALL BLACS_GRIDEXIT( CONTEXT )
*
      CALL BLACS_EXIT( 0 )
      STOP
*
 9999 FORMAT( A )
 9997 FORMAT( 'Check if overflow is handled in ieee default manner.' )
 9996 FORMAT( 'If this is the last output you see, you should assume')
 9995 FORMAT( 'that overflow caused a floating point exception.' )
*
 9986 FORMAT( 'Test ok. The system appears to handle ieee overflow.' )
*
 9985 FORMAT( 'Finished ', I6, ' tests, with the following results:' )
 9984 FORMAT( I5, ' tests completed and passed residual checks.' )
 9983 FORMAT( I5, ' tests completed without checking.' )
 9982 FORMAT( I5, ' tests skipped for lack of memory.' )
 9981 FORMAT( I5, ' tests completed and failed.' )
 9980 FORMAT( 'END OF TESTS.' )
 9979 FORMAT( '     N  NB   P   Q TYP SUB   WALL      CPU  ',
     $      '    CHK       QTQ    CHECK    TEST' )
 9978 FORMAT( ' ----- --- --- --- --- --- -------- --------',
     $      ' --------- --------- -----    ----' )
*
*     End of PSSEPRDRIVER
*
      END



