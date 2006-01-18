*
*
      PROGRAM PDSEPDRIVER
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     Parallel DOUBLE PRECISION symmetric eigenproblem test driver
*
*     The user should modify TOTMEM to indicate the maximum amount of
*     memory in bytes her system has.  Remember to leave room in memory
*     for operating system, the BLACS buffer, etc.  INTSIZ and DBLSIZ
*     indicate the length in bytes on the given platform for an integer
*     and a double precision real.
*     For example, on our system with 8 MB of memory, TOTMEM=6500000
*     (leaves 1.5 MB for OS, code, BLACS buffer, etc), the length of a
*     DOUBLE is 8, and an integer takes up 4 bytes.  Some playing around
*     to discover what the maximum value you can set MEMSIZ to may be
*     required.
*     All arrays used by factorization and solve are allocated out of
*     big array called MEM.
*
*     The full tester requires approximately (5 n + 5 n^2/p + slop)
*     DOUBLE PRECISION words and 6*n integer words.
*     So, TOTMEM should be set to at least 1.1 * 8 * (5n + 5n^2/p)
*
*     WHAT WE TEST
*     ============
*
*     This routine tests PDSYEVX, the expert driver for the parallel
*     symmetric eigenvalue problem, PDSYEV and PDSYEVD.  We would like 
*     to cover all possible combinations of:  matrix size, process 
*     configuration (nprow and npcol), block size (nb), 
*     matrix type (??), range of eigenvalue (all, by value, 
*     by position), sorting options, and upper vs. lower storage.
*
*     As PDSYEV returns an error message when heterogeneity is detected,
*     the PDSYEV tests can be suppressed by changing the appropiate
*     entry in the input file.
*
*     We intend to provide two types of test input files, an
*     installation test and a thorough test.
*
*     We also intend that the reports be meaningful.  Our input file
*     will allow multiple requests where each request is a cross product
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
      INTEGER            TOTMEM, DBLESZ, NIN
      PARAMETER          ( TOTMEM = 2000000, DBLESZ = 8, NIN = 11 )
      INTEGER            MEMSIZ
      PARAMETER          ( MEMSIZ = TOTMEM / DBLESZ )
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
      DOUBLE PRECISION   MEM( MEMSIZ )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
*
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT, 
     $                   BLACS_GRIDINIT, BLACS_PINFO, BLACS_SETUP, 
     $                   IGAMN2D, PDLACHKIEEE, PDLASNBT, PDSEPREQ 
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
         OPEN( UNIT = NIN, FILE = 'SEP.dat', STATUS = 'OLD' )
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
      CALL PDLASNBT( ISIEEE )
*
      CALL IGAMN2D( CONTEXT, 'a', ' ', 1, 1, ISIEEE, 1, 1, 1, -1, -1,
     $              0 )
*
      IF( ( ISIEEE.NE.0 ) ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9998 )
            WRITE( NOUT, FMT = 9997 )
            WRITE( NOUT, FMT = 9996 )
            WRITE( NOUT, FMT = 9995 )
            WRITE( NOUT, FMT = 9994 )
            WRITE( NOUT, FMT = 9993 )
            WRITE( NOUT, FMT = 9992 )
            WRITE( NOUT, FMT = 9991 )
            WRITE( NOUT, FMT = 9990 )
         END IF
*
         CALL PDLACHKIEEE( ISIEEE, DLAMCH( 'O' ), DLAMCH( 'U' ) )
*
         CALL IGAMN2D( CONTEXT, 'a', ' ', 1, 1, ISIEEE, 1, 1, 1, -1, -1,
     $                 0 )
*
         IF( ISIEEE.EQ.0 ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT, FMT = 9989 )
               WRITE( NOUT, FMT = 9988 )
               WRITE( NOUT, FMT = 9987 )
            END IF
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
     $      'SCALAPACK symmetric Eigendecomposition routines.'
         WRITE( NOUT, FMT = 9999 )USRINFO
         WRITE( NOUT, FMT = 9999 )' '
         WRITE( NOUT, FMT = 9999 )'Running tests of the parallel ' //
     $      'symmetric eigenvalue routine:  PDSYEVX & '//
     $       ' PDSYEV & PDSYEVD.'
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
     $      'than THRESH, RESULT is flagged as PASSED.'
         WRITE( NOUT, FMT = 9999 )
     $      '         : the QTQ norm is allowed to exceed THRESH' //
     $      ' for those eigenvectors'
         WRITE( NOUT, FMT = 9999 )'         :  which could not be ' //
     $      'reorthogonalized for lack of workspace.'
         WRITE( NOUT, FMT = 9999 )
     $      'TYP      : matrix type (see PDSEPtst.f).'
         WRITE( NOUT, FMT = 9999 )'SUB      : Subtests ' //
     $      '(see PDSEPtst).f'
         WRITE( NOUT, FMT = 9999 )'CHK      : ||AQ - QL|| ' //
     $      '/ ((abstol + ||A|| * eps) * N)'
         WRITE( NOUT, FMT = 9999 )'QTQ      : ||Q^T*Q - I||/ (N * eps)'
         WRITE( NOUT, FMT = 9999 )
     $      '         : when the adjusted QTQ exceeds THRESH',
     $      ' the adjusted QTQ norm is printed'
         WRITE( NOUT, FMT = 9999 )
     $      '         : otherwise the true QTQ norm is printed'
         WRITE( NOUT, FMT = 9999 )
     $      '           If NT>1, CHK and QTQ are the max over all ' //
     $      'eigen request tests'
         WRITE( NOUT, FMT = 9999 )
     $      'TEST     : EVX - testing PDSYEVX, EV - testing PDSYEV, '//
     $        'EVD - testing PDSYEVD'
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
      CALL PDSEPREQ( HETERO, NIN, MEM, MEMSIZ, NOUT, ISEED, NTESTS,
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
*
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
*
 9999 FORMAT( A )
 9998 FORMAT( ' I am about to check to make sure that overflow' )
 9997 FORMAT( ' is handled in the ieee default manner.  If this' )
 9996 FORMAT( ' is the last output you see, you should assume' )
 9995 FORMAT( ' that overflow caused a floating point exception.' )
 9994 FORMAT( ' In that case, we recommend that you add -DNO_IEEE' )
 9993 FORMAT( ' to the CDEFS line in SLmake.inc.' )
 9992 FORMAT( ' Alternatively, you could set CDEFS in SLmake.inc ' )
 9991 FORMAT( ' to enable the default ieee behaviour, However, this' )
 9990 FORMAT( ' may result in good or very bad performance.' )
 9989 FORMAT( ' Either signed zeroes or signed infinities ' )
 9988 FORMAT( ' work incorrectly or your system.  Change your' )
 9987 FORMAT( ' SLmake.inc as suggested above.' )
*
 9986 FORMAT( ' Your system appears to handle ieee overflow.' )
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
*     End of PDSEPDRIVER
*
      END
