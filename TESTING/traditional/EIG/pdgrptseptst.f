*
*
      PROGRAM PDRPTGSEPTST
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*     Repeat generalized parallel symmetric eigenproblem test
*     .. Parameters ..
*
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            MAXN, LWORK, LIWORK
      PARAMETER          ( MAXN = 200, LWORK = 500000,
     $                   LIWORK = 6*MAXN+4 )
*     ..
*     .. Local Scalars ..
      CHARACTER          SUBTESTS, UPLO
      INTEGER            CONTEXT, IAM, IBTYPE, INFO, IPOSTPAD, IPREPAD,
     $                   LDA, MATTYPE, N, NB, NPCOL, NPROCS, NPROW
      DOUBLE PRECISION   ABSTOL, THRESH
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), ICLUSTR( MAXN ), IFAIL( MAXN ),
     $                   ISEED( 4 ), IWORK( LIWORK )
      DOUBLE PRECISION   A( MAXN*MAXN ), B( MAXN, MAXN ),
     $                   COPYA( MAXN*MAXN ), COPYB( MAXN, MAXN ),
     $                   GAP( MAXN ), WIN( MAXN ), WNEW( MAXN ),
     $                   WORK( LWORK ), Z( MAXN*MAXN )
*     ..
*
*
*     .. External Subroutines ..
*
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDINIT,
     $                   BLACS_PINFO, BLACS_SETUP, DESCINIT, PDGSEPTST
*     ..
*     .. Executable Statements ..
*
      IPREPAD = 3
      IPOSTPAD = 3
      LDA = MAXN
*
*     These lines should be replaced by the output from pxGSEPdriver
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
      IBTYPE = 1
*  note: the printout often makes a mess of ABSTOL
      ABSTOL = 0.1175494351D-37
      THRESH = .350000D+01
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
      CALL PDGSEPTST( DESCA, UPLO, N, MATTYPE, IBTYPE, SUBTESTS, THRESH,
     $                N, ABSTOL, ISEED, A, COPYA, B, COPYB, Z, LDA, WIN,
     $                WNEW, IFAIL, ICLUSTR, GAP, IPREPAD, IPOSTPAD,
     $                WORK, LWORK-IPREPAD-IPOSTPAD, IWORK,
     $                LIWORK-IPREPAD-IPOSTPAD, 6, INFO )
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
*     End of PDRPTGSEPTST
*
      END
