*
*
      PROGRAM PSRPTSEPTST
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*     Repeat parallel symmetric eigenproblem test
*     .. Parameters ..
*
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            MAXN, LWORK, LIWORK
      PARAMETER          ( MAXN = 200, LWORK = 500000,
     $                   LIWORK = 6*MAXN+4 )
*     ..
*     .. Local Scalars ..
      CHARACTER          HETERO, SUBTESTS, UPLO
      INTEGER            CONTEXT, IAM, INFO, IPOSTPAD, IPREPAD, LDA,
     $                   MATTYPE, N, NB, NPCOL, NPROCS, NPROW
      REAL               ABSTOL, THRESH
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), ICLUSTR( MAXN ), IFAIL( MAXN ),
     $                   ISEED( 4 ), IWORK( LIWORK )
      REAL               A( MAXN*MAXN ), COPYA( MAXN*MAXN ), 
     $                   GAP( MAXN ), WIN( MAXN ), WNEW( MAXN ), 
     $                   WORK( LWORK ), Z( MAXN*MAXN )
*     ..
*
*
*     .. External Subroutines ..
*
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDINIT,
     $                   BLACS_PINFO, BLACS_SETUP, DESCINIT, PSSEPTST
*     ..
*     .. Executable Statements ..
*
      IPREPAD = 3
      IPOSTPAD = 3
      LDA = MAXN
*
*     Set HETERO to 'Y' if you want to turn off the PxSYEV tests
*
      HETERO = 'N'
*
*     These lines should be replaced by the output from pxSEPdriver
*
*
      ISEED( 1 ) = 2312
      ISEED( 2 ) = 3709
      ISEED( 3 ) = 666
      ISEED( 4 ) = 3371
      UPLO = 'U'
      SUBTESTS = 'Y'
      N = 33
      NPROW = 2
      NPCOL = 2
      NB = 4
      MATTYPE = 9
*  note: the printout often makes a mess of ABSTOL
      ABSTOL = 0.1175494351E-37
      THRESH = .350000E+01
*
      CALL BLACS_PINFO( IAM, NPROCS )
      IF( NPROCS.LT.1 ) THEN
*
         NPROCS = NPROW*NPCOL
         CALL BLACS_SETUP( IAM, NPROCS )
      END IF
      CALL BLACS_GET( -1, 0, CONTEXT )
      CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
*
      CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
*
      CALL PSSEPTST( DESCA, UPLO, N, MATTYPE, SUBTESTS, THRESH, N,
     $               ABSTOL, ISEED, A, COPYA, Z, LDA, WIN, WNEW, IFAIL,
     $               ICLUSTR, GAP, IPREPAD, IPOSTPAD, WORK,
     $               LWORK-IPREPAD-IPOSTPAD, IWORK,
     $               LIWORK-IPREPAD-IPOSTPAD, HETERO, 6, INFO )
*
*
*
*     Uncomment this line on SUN systems to avoid the useless print out
*
*      CALL IEEE_FLAGS( 'clear', 'exception', 'underflow', '')
*
*
*
*
      CALL BLACS_EXIT( 0 )
      STOP
*
*
*
*     End of PSRPTSEPTST
*
      END
