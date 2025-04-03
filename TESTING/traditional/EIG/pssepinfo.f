*
*
      SUBROUTINE PSSEPINFO( CONTEXT, IAM, NIN, NOUT, MAXSETSIZE,
     $                      NMATSIZES, MATSIZES, NUPLOS, UPLOS,
     $                      NPCONFIGS, NPROWS, NPCOLS, NBS, NMATTYPES,
     $                      MATTYPES, MAXTYPE, SUBTESTS, THRESH, ORDER,
     $                      ABSTOL, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*  Purpose
*  =======
*
*  PSSEPINFO reads the input test data file (INFILE), copies the
*  information therein to all processes and returns this information
*  in the corresponding parameters.
*
*  Arguments
*  =========
*
*  CONTEXT    (global input) INTEGER
*             BLACS Context
*
*  IAM        (local input) INTEGER
*             process number.
*             IAM.EQ.0 on the proceesor that performs I/O
*
*  NIN        (global input) INTEGER
*             The unit number of the input file.
*
*  NOUT       (global output) INTEGER
*             The unit number for output file.
*             if NOUT = 6, ouput to screen,
*             if NOUT = 0, output to stderr
*             Only defined for process 0.
*
*  MAXSETSIZE (global output) INTEGER
*             Maximum set size.  Size of the following arrays:
*             MATSIZES, MATTYPES, NBS, NPCOLS, NPROWS
*
*  NMATSIZES  (global output) INTEGER
*             Number of matrix sizes to test
*
*  MATSIZES   (global output) INTEGER array dimension MAXSETSIZE
*             Matrix sizes to test
*
*  NUPLOS     (global output) INTEGER
*             Number of UPLO values to test
*
*  UPLOS      (global output) CHARACTER*1 array dimension 2
*             Values of UPLO to test
*
*  NPCONFIGS  (global output) INTEGER
*             Number of process configuratins (NPROW, NPCOL, NB)
*
*  NPROWS     (global output) INTEGER array dimension MAXSETSIZE
*             Values of NPROW to test
*
*  NPCOLS     (global output) INTEGER array dimension MAXSETSIZE
*             Values of NPCOL to test
*
*  NBS       (global output) INTEGER array dimension MAXSETSIZE
*             Values of NB to test
*
*  NMATTYPES  (global output) INTEGER
*             Number of matrix types to test
*
*  MATTYPES   (global output) INTEGER array dimension MAXSETSIZE
*             Matrix types to test
*             Refer to PSSEPTST for a complete description of the
*             supported matrix types.
*
*  MAXTYPE    (global input) INTEGER
*             Maximum allowed matrix type
*
*  SUBTESTS  (global output) CHARACTER
*             'N' = Do not perform subtests
*             'Y' = Perfrom subtests
*
*
*  THRESH   (global output) @(tupc)
*          A test will count as "failed" if the "error", computed as
*          described below, exceeds THRESH.  Note that the error
*          is scaled to be O(1), so THRESH should be a reasonably
*          small multiple of 1, e.g., 10 or 100.  In particular,
*          it should not depend on the precision (single vs. double)
*          or the size of the matrix.  It must be at least zero.
*          ( THRESH is set to 1/10 of the value defined in the .dat
*          file when NOUT = 13.  THRESH is set to 1/20 of the value
*          defined in the .dat file when NOUT = 14. This allows us
*          to specify more stringent criteria for our internal testing )
*
*  ORDER    (global output) INTEGER
*           Number of reflectors used in test matrix creation.
*           If ORDER is large, it will
*           take more time to create the test matrices but they will
*           be closer to random.
*           ORDER .lt. N not implemented
*
*  ABSTOL   (global output) REAL
*           The absolute tolerance for the eigenvalues. An
*           eigenvalue is considered to be located if it has
*           been determined to lie in an interval whose width
*           is "abstol" or less. If "abstol" is less than or equal
*           to zero, then ulp*|T| will be used, where |T| is
*           the 1-norm of the matrix. If eigenvectors are
*           desired later by inverse iteration ("PSSTEIN"),
*           "abstol" MUST NOT be bigger than ulp*|T|.
*
*           If ( ABSTOL .EQ. 0 in SEP.dat, it is set to
*           2.0 * PSLAMCH( 'u' ) in this routine.
*
*  INFO   (global output) INTEGER
*           0 = normal return
*           -1 = end of file
*           -2 = incorrrect data specification
*
*     .. Scalar Arguments ..
      CHARACTER          SUBTESTS
      INTEGER            CONTEXT, IAM, INFO, MAXSETSIZE, MAXTYPE, NIN,
     $                   NMATSIZES, NMATTYPES, NOUT, NPCONFIGS, NUPLOS,
     $                   ORDER
      REAL               ABSTOL, THRESH
*     ..
*     .. Array Arguments ..
      CHARACTER          UPLOS( 2 )
      INTEGER            MATSIZES( MAXSETSIZE ), MATTYPES( MAXSETSIZE ),
     $                   NBS( MAXSETSIZE ), NPCOLS( MAXSETSIZE ),
     $                   NPROWS( MAXSETSIZE )
*     ..
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               TWO, TEN, TWENTY
      PARAMETER          ( TWO = 2.0E0, TEN = 10.0E0, TWENTY = 20.0E0 )
*     ..
*     .. Local Scalars ..
      CHARACTER*80       TESTSUMMRY
      INTEGER            I, ISUBTESTS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               PSLAMCH
      EXTERNAL           LSAME, PSLAMCH
*     ..
*
*     .. External Subroutines ..
      EXTERNAL           IGEBR2D, IGEBS2D, SGEBR2D, SGEBS2D
*     ..
*
*     .. Local Arrays ..
      INTEGER            IUPLOS( 2 )
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = 9997 )TESTSUMMRY
         TESTSUMMRY = ' '
         READ( NIN, FMT = 9997 )TESTSUMMRY
         WRITE( NOUT, FMT = 9997 )TESTSUMMRY
      END IF
*
*     assign a default
      INFO = 0
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )NMATSIZES
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, 1, NMATSIZES, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, 1, NMATSIZES, 1, 0, 0 )
      END IF
      IF( NMATSIZES.EQ.-1 ) THEN
         INFO = -1
         GO TO 70
      END IF
      IF( NMATSIZES.LT.1 .OR. NMATSIZES.GT.MAXSETSIZE ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9999 )'Matrix size', NMATSIZES, 1,
     $         MAXSETSIZE
         END IF
         INFO = -2
         GO TO 70
      END IF
*
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )( MATSIZES( I ), I = 1, NMATSIZES )
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, NMATSIZES, MATSIZES, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, NMATSIZES, MATSIZES, 1,
     $                 0, 0 )
      END IF
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )NUPLOS
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, 1, NUPLOS, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, 1, NUPLOS, 1, 0, 0 )
      END IF
      IF( NUPLOS.LT.1 .OR. NUPLOS.GT.2 ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9999 )'# of UPLOs', NUPLOS, 1, 2
         END IF
         INFO = -2
         GO TO 70
      END IF
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )( UPLOS( I ), I = 1, NUPLOS )
         DO 10 I = 1, NUPLOS
            IF( LSAME( UPLOS( I ), 'L' ) ) THEN
               IUPLOS( I ) = 1
            ELSE
               IUPLOS( I ) = 2
            END IF
   10    CONTINUE
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, NUPLOS, IUPLOS, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, NUPLOS, IUPLOS, 1, 0, 0 )
      END IF
      DO 20 I = 1, NUPLOS
         IF( IUPLOS( I ).EQ.1 ) THEN
            UPLOS( I ) = 'L'
         ELSE
            UPLOS( I ) = 'U'
         END IF
   20 CONTINUE
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )NPCONFIGS
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, 1, NPCONFIGS, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, 1, NPCONFIGS, 1, 0, 0 )
      END IF
      IF( NPCONFIGS.LT.1 .OR. NPCONFIGS.GT.MAXSETSIZE ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9999 )'# proc configs', NPCONFIGS, 1,
     $         MAXSETSIZE
         END IF
         INFO = -2
         GO TO 70
      END IF
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )( NPROWS( I ), I = 1, NPCONFIGS )
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, NPCONFIGS, NPROWS, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, NPCONFIGS, NPROWS, 1, 0,
     $                 0 )
      END IF
      DO 30 I = 1, NPCONFIGS
         IF( NPROWS( I ).LE.0 )
     $      INFO = -2
   30 CONTINUE
      IF( INFO.EQ.-2 ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' NPROW'
         END IF
         GO TO 70
      END IF
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )( NPCOLS( I ), I = 1, NPCONFIGS )
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, NPCONFIGS, NPCOLS, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, NPCONFIGS, NPCOLS, 1, 0,
     $                 0 )
      END IF
      DO 40 I = 1, NPCONFIGS
         IF( NPCOLS( I ).LE.0 )
     $      INFO = -2
   40 CONTINUE
      IF( INFO.EQ.-2 ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' NPCOL'
         END IF
         GO TO 70
      END IF
*
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )( NBS( I ), I = 1, NPCONFIGS )
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, NPCONFIGS, NBS, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, NPCONFIGS, NBS, 1, 0, 0 )
      END IF
      DO 50 I = 1, NPCONFIGS
         IF( NBS( I ).LE.0 )
     $      INFO = -2
   50 CONTINUE
      IF( INFO.EQ.-2 ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9996 )' NB'
         END IF
         GO TO 70
      END IF
*
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )NMATTYPES
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, 1, NMATTYPES, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, 1, NMATTYPES, 1, 0, 0 )
      END IF
      IF( NMATTYPES.LT.1 .OR. NMATTYPES.GT.MAXSETSIZE ) THEN
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9999 )'matrix types', NMATTYPES, 1,
     $         MAXSETSIZE
         END IF
         INFO = -2
         GO TO 70
      END IF
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )( MATTYPES( I ), I = 1, NMATTYPES )
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, NMATTYPES, MATTYPES, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, NMATTYPES, MATTYPES, 1,
     $                 0, 0 )
      END IF
*
      DO 60 I = 1, NMATTYPES
         IF( MATTYPES( I ).LT.1 .OR. MATTYPES( I ).GT.MAXTYPE ) THEN
            IF( IAM.EQ.0 ) THEN
               WRITE( NOUT, FMT = 9999 )'matrix type', MATTYPES( I ),
     $            1, MAXTYPE
            END IF
            MATTYPES( I ) = 1
         END IF
   60 CONTINUE
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )SUBTESTS
         IF( LSAME( SUBTESTS, 'Y' ) ) THEN
            ISUBTESTS = 2
         ELSE
            ISUBTESTS = 1
         END IF
         CALL IGEBS2D( CONTEXT, 'All', ' ', 1, 1, ISUBTESTS, 1 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'All', ' ', 1, 1, ISUBTESTS, 1, 0, 0 )
      END IF
      IF( ISUBTESTS.EQ.2 ) THEN
         SUBTESTS = 'Y'
      ELSE
         SUBTESTS = 'N'
      END IF
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )THRESH
         IF( NOUT.EQ.13 )
     $      THRESH = THRESH / TEN
         IF( NOUT.EQ.14 )
     $      THRESH = THRESH / TWENTY
         CALL SGEBS2D( CONTEXT, 'All', ' ', 1, 1, THRESH, 1 )
      ELSE
         CALL SGEBR2D( CONTEXT, 'All', ' ', 1, 1, THRESH, 1, 0, 0 )
      END IF
*
      ORDER = 0
*
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * )ABSTOL
         CALL SGEBS2D( CONTEXT, 'All', ' ', 1, 1, ABSTOL, 1 )
      ELSE
         CALL SGEBR2D( CONTEXT, 'All', ' ', 1, 1, ABSTOL, 1, 0, 0 )
      END IF
      IF( ABSTOL.LT.0 )
     $   ABSTOL = TWO*PSLAMCH( CONTEXT, 'U' )
*
      INFO = 0
*
   70 CONTINUE
      RETURN
*
 9999 FORMAT( A20, ' is:', I5, ' must be between:', I5, ' and', I5 )
 9998 FORMAT( A20, ' is:', I5, ' must be:', I5, ' or', I5 )
 9997 FORMAT( A )
 9996 FORMAT( A20, ' must be positive' )
*
*     End of PSSEPINFO
*
      END
