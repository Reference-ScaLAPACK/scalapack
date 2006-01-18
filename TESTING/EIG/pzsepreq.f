*
*
      SUBROUTINE PZSEPREQ( NIN, MEM, MEMSIZE, NOUT, ISEED, NTESTS,
     $                     NSKIPPED, NNOCHECK, NPASSED, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INFO, MEMSIZE, NIN, NNOCHECK, NOUT, NPASSED,
     $                   NSKIPPED, NTESTS
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      COMPLEX*16         MEM( MEMSIZE )
*     ..
*
*  Purpose
*  =======
*
*  PZSEPREQ performs one request from the input file 'SEP.dat'
*  A request is the cross product of the specifications in the
*  input file.  PZSEPREQ prints one line per test.
*
*  Arguments
*  =========
*
*  NIN      (local input) INTEGER
*           The unit number for the input file 'SEP.dat'
*
*  MEM      (local input) COMPLEX*16         ARRAY, dimension MEMSIZE
*           Array encompassing the available single precision memory
*
*  MEMSIZE  (local input)  INTEGER
*           Size of MEM array
*
*  NOUT   (local input) INTEGER
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
*            Random number generator seed
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
*
*
*  INFO     (global output) INTEGER
*           0 = test request ran
*           -1 = end of file
*           -2 = incorrect .dat file
*
*
*     .. Parameters ..
*
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            ZPLXSZ, INTGSZ
      PARAMETER          ( ZPLXSZ = 16, INTGSZ = 4 )
      INTEGER            DBLESZ
      PARAMETER          ( DBLESZ = 8 )
      INTEGER            MAXSETSIZE
      PARAMETER          ( MAXSETSIZE = 50 )
*     ..
*     .. Local Scalars ..
      CHARACTER          SUBTESTS
      INTEGER            CONTEXT, IAM, IMIDPAD, INITCON, IPOSTPAD,
     $                   IPREPAD, ISIZEHEEVX, ISIZESUBTST, ISIZETST,
     $                   LDA, LLRWORK, MATSIZE, MATTYPE, MYCOL, MYROW,
     $                   N, NB, NMATSIZES, NMATTYPES, NNODES, NP, NPCOL,
     $                   NPCONFIGS, NPROW, NQ, NUPLOS, ORDER, PCONFIG,
     $                   PTRA, PTRCOPYA, PTRGAP, PTRICLUS, PTRIFAIL,
     $                   PTRIWRK, PTRRWORK, PTRW, PTRW2, PTRWORK, PTRZ,
     $                   RES, RSIZECHK, RSIZEHEEVX, RSIZEQTQ,
     $                   RSIZESUBTST, RSIZETST, SIZEHEEVX, SIZEMQRLEFT,
     $                   SIZEMQRRIGHT, SIZEQRF, SIZESUBTST, SIZETMS,
     $                   SIZETST, UPLO, SIZEHEEVD, RSIZEHEEVD,
     $                   ISIZEHEEVD
      DOUBLE PRECISION   ABSTOL, THRESH
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
     $                   DESCINIT, PDSEPINFO, PZLASIZESEP, PZSEPTST
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      CALL BLACS_PINFO( IAM, NNODES )
      CALL BLACS_GET( -1, 0, INITCON )
      CALL BLACS_GRIDINIT( INITCON, 'R', 1, NNODES )
*
      CALL PDSEPINFO( INITCON, IAM, NIN, NOUT, MAXSETSIZE, NMATSIZES,
     $                MATSIZES, NUPLOS, UPLOS, NPCONFIGS, NPROWS,
     $                NPCOLS, NBS, NMATTYPES, MATTYPES, 22, SUBTESTS,
     $                THRESH, ORDER, ABSTOL, INFO )
*
      CALL BLACS_GRIDEXIT( INITCON )
*
      IF( INFO.EQ.0 ) THEN
*
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
                     IF( MYROW.GE.0 ) THEN
                        CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0,
     $                                 CONTEXT, LDA, INFO )
                        CALL PZLASIZESEP( DESCA, IPREPAD, IPOSTPAD,
     $                                    SIZEMQRLEFT, SIZEMQRRIGHT,
     $                                    SIZEQRF, SIZETMS, RSIZEQTQ,
     $                                    RSIZECHK, SIZEHEEVX,
     $                                    RSIZEHEEVX, ISIZEHEEVX,
     $                                    SIZEHEEVD, RSIZEHEEVD,
     $                                    ISIZEHEEVD,
     $                                    SIZESUBTST, RSIZESUBTST,
     $                                    ISIZESUBTST, SIZETST,
     $                                    RSIZETST, ISIZETST )
*
                        PTRA = 1
                        PTRZ = PTRA + LDA*NQ + IPREPAD + IPOSTPAD
                        PTRCOPYA = PTRZ + LDA*NQ + IPREPAD + IPOSTPAD
                        PTRW = PTRCOPYA + LDA*NQ + IPREPAD + IPOSTPAD
                        PTRW2 = PTRW + ICEIL( MAX( N, 1 )+IPREPAD+
     $                          IPOSTPAD, ZPLXSZ / DBLESZ )
                        PTRWORK = PTRW2 + ICEIL( MAX( N, 1 )+IPREPAD+
     $                            IPOSTPAD, ZPLXSZ / DBLESZ )
                        PTRGAP = PTRWORK + SIZETST + IPREPAD + IPOSTPAD
                        PTRIFAIL = PTRGAP + ICEIL( NPROW*NPCOL+IPREPAD+
     $                             IPOSTPAD, ZPLXSZ / DBLESZ )
                        PTRICLUS = PTRIFAIL + ICEIL( N+IPREPAD+IPOSTPAD,
     $                             ZPLXSZ / INTGSZ )
                        PTRIWRK = PTRICLUS + ICEIL( 2*NPROW*NPCOL+
     $                            IPREPAD+IPOSTPAD, ZPLXSZ / INTGSZ )
                        PTRRWORK = PTRIWRK + ICEIL( ISIZETST+IPREPAD+
     $                             IPOSTPAD, ZPLXSZ / INTGSZ )
                        LLRWORK = ( MEMSIZE-PTRRWORK+1 )*ZPLXSZ / DBLESZ
C                       LLRWORK = ( MEMSIZE-PTRRWORK-IPREPAD-IPOSTPAD+1)
C    $                       *ZPLXSZ / DBLESZ
                        NTESTS = NTESTS + 1
                        IF( LLRWORK.LT.RSIZETST ) THEN
                           NSKIPPED = NSKIPPED + 1
                        ELSE
                           CALL PZSEPTST( DESCA, UPLOS( UPLO ), N,
     $                                    MATTYPES( MATTYPE ), SUBTESTS,
     $                                    THRESH, N, ABSTOL, ISEED,
     $                                    MEM( PTRA ), MEM( PTRCOPYA ),
     $                                    MEM( PTRZ ), LDA, MEM( PTRW ),
     $                                    MEM( PTRW2 ), MEM( PTRIFAIL ),
     $                                    MEM( PTRICLUS ),
     $                                    MEM( PTRGAP ), IPREPAD,
     $                                    IPOSTPAD, MEM( PTRWORK ),
     $                                    SIZETST, MEM( PTRRWORK ),
     $                                    LLRWORK, MEM( PTRIWRK ),
     $                                    ISIZETST, NOUT, RES )
*
                           IF( RES.EQ.0 ) THEN
                              NPASSED = NPASSED + 1
                           ELSE IF( RES.EQ.2 ) THEN
                              NNOCHECK = NNOCHECK + 1
                           ELSE IF( RES.EQ.3 ) THEN
                              NSKIPPED = NSKIPPED + 1
                              WRITE( NOUT, FMT=* )' PZSEPREQ failed'
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
*
      RETURN
*
*     End of PZDSEPREQ
*
      END
