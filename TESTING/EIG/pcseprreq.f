      SUBROUTINE PCSEPRREQ( HETERO, NIN, MEM, MEMSIZE, NOUT, ISEED,
     $                     NTESTS, NSKIPPED, NNOCHECK, NPASSED, INFO )
*
*  -- ScaLAPACK routine (@(MODE)version *TBA*) --
*     University of California, Berkeley and
*     University of Tennessee, Knoxville. 
*     October 21, 2006
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          HETERO
      INTEGER            INFO, MEMSIZE, NIN, NNOCHECK, NOUT, NPASSED,
     $                   NSKIPPED, NTESTS
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      COMPLEX            MEM( MEMSIZE )     
*
*  Purpose
*  =======
*
*  PCSEPRREQ performs one request from the input file 'SEPR.dat'
*  A request is the cross product of the specifications in the
*  input file. It prints one line per test.
*
*  Arguments
*  =========
*
*  NIN      (local input) INTEGER
*           The unit number for the input file 'SEPR.dat'
*
*  MEM      (local input ) COMPLEX          ARRAY, dimension MEMSIZE
*           Array encompassing the available single precision memory
*
*  MEMSIZE  (local input)  INTEGER
*           Size of MEM array
*
*  NOUT     (local input) INTEGER
*           The unit number for output file.
*           NOUT = 6, output to screen,
*           NOUT = 0, output to stderr.
*           NOUT = 13, output to file, divide thresh by 10
*           NOUT = 14, output to file, divide thresh by 20
*           Only used on node 0.
*           NOUT = 13, 14 allow the threshold to be tighter for our
*           internal testing which means that when a user reports
*           a threshold error, it is more likely to be significant.
*
*  ISEED    (global input/output) INTEGER array, dimension 4
*           Random number generator seed
*
*  NTESTS   (global input/output) INTEGER
*           NTESTS = NTESTS + tests requested
*
*  NSKIPPED (global input/output) INTEGER
*           NSKIPPED = NSKIPPED + tests skipped
*
*  NNOCHECK (global input/output) INTEGER
*           NNOCHECK = NNOCHECK + tests completed but not checked
*
*  NPASSED  (global input/output) INTEGER
*           NPASSED = NPASSED + tests which passed all checks
*
*  INFO     (global output) INTEGER
*           0 = test request ran
*          -1 = end of file
*          -2 = incorrect .dat file
*
*     .. Parameters ..
*
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            REALSZ, INTGSZ
      PARAMETER          ( REALSZ = 4, INTGSZ = 4 )
      INTEGER            KMPXSZ
      PARAMETER          ( KMPXSZ = 8 )
      INTEGER            MAXSETSIZE
      PARAMETER          ( MAXSETSIZE = 50 )
*     ..
*     .. Local Scalars ..
      CHARACTER          SUBTESTS
      INTEGER            CONTEXT, IAM, IMIDPAD, INITCON, IPOSTPAD,
     $                   IPREPAD, ISIZESUBTST, ISIZEEVR, ISIZETST,
     $                   LDA, LLWORK, MATSIZE, MATTYPE, MYCOL, MYROW, N,
     $                   NB, NMATSIZES, NMATTYPES, NNODES, NP, NPCOL,
     $                   NPCONFIGS, NPROW, NQ, NUPLOS, ORDER, PCONFIG,
     $                   PTRA, PTRCOPYA, PTRGAP, PTRICLUS, PTRIFAIL,
     $                   PTRIWRK, PTRW, PTRW2, PTRWORK, PTRZ, RES,
     $                   SIZECHK, SIZEMQRLEFT, SIZEMQRRIGHT, SIZEQRF,
     $                   SIZEQTQ, SIZESUBTST, SIZEEVR,
     $                   SIZETMS, SIZETST, UPLO
      INTEGER            PTRRWORK, RSIZEEVR, RSIZESUBTST, RSIZETST
*
      REAL               ABSTOL, THRESH
*     ..
*     .. Local Arrays ..
      CHARACTER          UPLOS( 2 )
      INTEGER            DESCA( DLEN_ ), MATSIZES( MAXSETSIZE ),
     $                   MATTYPES( MAXSETSIZE ), NBS( MAXSETSIZE ),
     $                   NPCOLS( MAXSETSIZE ), NPROWS( MAXSETSIZE )
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, NUMROC
      EXTERNAL           ICEIL, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GET, BLACS_GRIDEXIT, 
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO, 
     $                   DESCINIT, PCLASIZESEPR, PSSEPINFO, PCSEPRTST
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_PINFO( IAM, NNODES )
      CALL BLACS_GET( -1, 0, INITCON )
      CALL BLACS_GRIDINIT( INITCON, 'R', 1, NNODES )
*
      CALL PSSEPINFO( INITCON, IAM, NIN, NOUT, MAXSETSIZE, NMATSIZES,
     $                MATSIZES, NUPLOS, UPLOS, NPCONFIGS, NPROWS,
     $                NPCOLS, NBS, NMATTYPES, MATTYPES, 22, SUBTESTS,
     $                THRESH, ORDER, ABSTOL, INFO )
*
      CALL BLACS_GRIDEXIT( INITCON )
*
      IF( INFO.EQ.0 ) THEN
*
         DO 40 MATSIZE = 1, NMATSIZES
*
            DO 30 PCONFIG = 1, NPCONFIGS
*
               DO 20 MATTYPE = 1, NMATTYPES
*
                  DO 10 UPLO = 1, NUPLOS
*
                     N = MATSIZES( MATSIZE )
                     ORDER = N
*
                     NPROW = NPROWS( PCONFIG )
                     NPCOL = NPCOLS( PCONFIG )
                     NB = NBS( PCONFIG )
*
                     NP = NUMROC( N, NB, 0, 0, NPROW )
                     NQ = NUMROC( N, NB, 0, 0, NPCOL )
                     IPREPAD = MAX( NB, NP )
                     IMIDPAD = NB
                     IPOSTPAD = MAX( NB, NQ )
*
                     LDA = MAX( NP, 1 ) + IMIDPAD
*
                     CALL BLACS_GET( -1, 0, CONTEXT )
                     CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
                     CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW,
     $                                    MYCOL )
*
                     IF( MYROW.GE.0 ) THEN
                        CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0,
     $                                 CONTEXT, LDA, INFO )
                        CALL PCLASIZESEPR( DESCA, IPREPAD, IPOSTPAD,
     $                                     SIZEMQRLEFT, SIZEMQRRIGHT,
     $                                     SIZEQRF, SIZETMS, SIZEQTQ,
     $                                     SIZECHK, SIZEEVR, RSIZEEVR,
     $                                     ISIZEEVR, SIZESUBTST, 
     $                                     RSIZESUBTST, ISIZESUBTST,
     $                                     SIZETST, RSIZETST, ISIZETST )
*
                        PTRA = 1
                        PTRZ = PTRA + LDA*NQ + IPREPAD + IPOSTPAD
                        PTRCOPYA = PTRZ + LDA*NQ + IPREPAD + IPOSTPAD
                        PTRW = PTRCOPYA + LDA*NQ + IPREPAD + IPOSTPAD
                        PTRW2 = PTRW + ICEIL( MAX( N, 1 )+IPREPAD+
     $                          IPOSTPAD, KMPXSZ / REALSZ )
                        PTRWORK = PTRW2 + ICEIL( MAX( N, 1 )+IPREPAD+
     $                            IPOSTPAD, KMPXSZ / REALSZ )
                        PTRGAP = PTRWORK + SIZETST + IPREPAD + IPOSTPAD
                        PTRIFAIL = PTRGAP + ICEIL( NPROW*NPCOL+IPREPAD+
     $                             IPOSTPAD, KMPXSZ / REALSZ )
                        PTRICLUS = PTRIFAIL + ICEIL( N+IPREPAD+IPOSTPAD,
     $                             KMPXSZ / INTGSZ )
                        PTRIWRK = PTRICLUS + ICEIL( 2*NPROW*NPCOL+
     $                            IPREPAD+IPOSTPAD, KMPXSZ / INTGSZ )
                        PTRRWORK = PTRIWRK + ICEIL( ISIZETST+IPREPAD+
     $                             IPOSTPAD, KMPXSZ / INTGSZ )
                        LLWORK = ( MEMSIZE-PTRRWORK+1 )*KMPXSZ / REALSZ

                        NTESTS = NTESTS + 1
                        IF( LLWORK.LT.RSIZETST ) THEN
                           NSKIPPED = NSKIPPED + 1
                        ELSE
                           CALL PCSEPRTST( DESCA, UPLOS( UPLO ), N,
     $                                    MATTYPES( MATTYPE ), SUBTESTS,
     $                                    THRESH, N, ABSTOL, ISEED,
     $                                    MEM( PTRA ), MEM( PTRCOPYA ),
     $                                    MEM( PTRZ ), LDA, MEM( PTRW ),
     $                                    MEM( PTRW2 ), MEM( PTRIFAIL ),
     $                                    MEM( PTRICLUS ),
     $                                    MEM( PTRGAP ), IPREPAD,
     $                                    IPOSTPAD, MEM( PTRWORK ),
     $                                    SIZETST, MEM( PTRRWORK ),
     $                                    LLWORK, MEM( PTRIWRK ),
     $                                    ISIZETST, HETERO, NOUT, RES )
*
                           IF( RES.EQ.0 ) THEN
                              NPASSED = NPASSED + 1
                           ELSE IF( RES.EQ.2 ) THEN
                              NNOCHECK = NNOCHECK + 1
                           ELSE IF( RES.EQ.3 ) THEN
                              NSKIPPED = NSKIPPED + 1
                              WRITE( NOUT, FMT = * )' PCSEPRREQ failed'
                              CALL BLACS_ABORT( CONTEXT, -1 )
                           END IF
                        END IF
                        CALL BLACS_GRIDEXIT( CONTEXT )
                     END IF
   10             CONTINUE
   20          CONTINUE
   30       CONTINUE
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of PCSEPRREQ
*
      END
