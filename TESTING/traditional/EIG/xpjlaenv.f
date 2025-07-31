      INTEGER          FUNCTION PJLAENV( ICTXT, ISPEC, NAME, OPTS, N1,
     $                 N2, N3, N4 )
*
*  -- ScaLAPACK test routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 2, 2000
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ICTXT, ISPEC, N1, N2, N3, N4
*     ..
*
*  xpjlaenv.f versus pjlaenv.f
*  ===========================
*
*  xpjlaenv.f is used during testing to allow the timer/tester to
*  control pjlaenv's return values by setting common variables.
*  xpjlaenv.f guarantees that the return value is the same as the
*  corresponding value in common.  xpjlaenv.f either reads values
*  from common and uses them as return values or it writes the
*  return value to common.  Either way, xpjlaenv.f's return
*  value and the correpsonding value in common will always match.
*
*  When the common variable "TIMING" is set, the other common
*  variables are set to the values returned by xpjlaenv.f, else
*  xpjlaenv.f returns the values as set in common.
*
*  Purpose
*
*  =======
*
*  PJLAENV is called from the ScaLAPACK symmetric and Hermitian
*  tailored eigen-routines to choose
*  problem-dependent parameters for the local environment.  See ISPEC
*  for a description of the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (global input) INTEGER
*          Specifies the parameter to be returned as the value of
*          PJLAENV.
*          = 1: the data layout blocksize;
*          = 2: the panel blocking factor;
*          = 3: the algorithmic blocking factor;
*          = 4: execution path control;
*          = 5: maximum size for direct call to the LAPACK routine
*
*  NAME    (global input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (global input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (global input) INTEGER
*  N2      (global input) INTEGER
*  N3      (global input) INTEGER
*  N4      (global input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
*          At present, only N1 is used, and it (N1) is used only for
*          'TTRD'
*
* (PJLAENV) (global or local output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if PJLAENV = -k, the k-th argument had an illegal
*          value.
*
*          Most parameters set via a call to PJLAENV must be identical
*          on all processors and hence PJLAENV will return the same
*          value to all procesors (i.e. global output).  However some,
*          in particular, the panel blocking factor can be different
*          on each processor and hence PJLAENV can return different
*          values on different processors (i.e. local output).
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling PJLAENV from
*  the ScaLAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by PJLAENV is checked for validity
*      in the calling subroutine.  For example, PJLAENV is used to
*      retrieve the optimal blocksize for STRTRI as follows:
*
*      NB = PJLAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  PJLAENV is patterned after ILAENV and keeps the same interface in
*  anticipation of future needs, even though PJLAENV is only sparsely
*  used at present in ScaLAPACK.  Most ScaLAPACK codes use the input
*  data layout blocking factor as the algorithmic blocking factor -
*  hence there is no need or opportunity to set the algorithmic or
*  data decomposition blocking factor.
*
*  pXYYtevx.f and pXYYtgvx.f and pXYYttrd.f are the only codes which
*  call PJLAENV in this release.  pXYYtevx.f and pXYYtgvx.f redistribute
*  the data to the best data layout for each transformation.  pXYYttrd.f
*  uses a data layout blocking factor of 1 and a
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      LOGICAL            CNAME, GLOBAL, SNAME, TIME
      CHARACTER          C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*8        SUBNAM
      INTEGER            I, IC, IDUMM, IZ, MSZ, NB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR
*     ..
*
*
*     .. Scalars in Common ..
      INTEGER            ANB, BALANCED, BCKBLOCK, GSTBLOCK, INTERLEAVE,
     $                   LLTBLOCK, MINSZ, PNB, TIMING, TRSBLOCK,
     $                   TWOGEMMS
*     ..
*     .. External Subroutines ..
      EXTERNAL           IGAMX2D
*     ..
*     .. Common blocks ..
      COMMON             / BLOCKSIZES / GSTBLOCK, LLTBLOCK, BCKBLOCK,
     $                   TRSBLOCK
      COMMON             / MINSIZE / MINSZ
      COMMON             / PJLAENVTIMING / TIMING
      COMMON             / TAILOREDOPTS / PNB, ANB, INTERLEAVE,
     $                   BALANCED, TWOGEMMS
*     ..
*     .. Executable Statements ..
*
      TIME = ( TIMING.EQ.1 )
*
*
      GO TO ( 10, 10, 10, 10, 10 )ISPEC
*
*     Invalid value for ISPEC
*
      PJLAENV = -1
      RETURN
*
   10 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      PJLAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.100 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $             I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 2: 2 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 3: 4 )
      C3 = SUBNAM( 5: 7 )
      C4 = C3( 2: 3 )
*
*     This is to keep ftnchek happy
*
      IF( ( N2+N3+N4 )*0.NE.0 ) THEN
         C4 = OPTS
         C3 = C4
      END IF
*
      GO TO ( 50, 60, 70, 80, 90 )ISPEC
*
   50 CONTINUE
*
*     ISPEC = 1:  data layout block size
*     (global - all processes must use the same value)
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'SY' .OR. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'LLT' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
            IF( TIME ) THEN
               LLTBLOCK = NB
            ELSE
               NB = LLTBLOCK
               IF( NB.LE.0 ) THEN
                  PRINT *, 'xpjlaenv.f ERROR common variable LLTBLOCK',
     $               ' may be unitialized'
c                 CALL EXIT( 13 )
                  STOP
               END IF
            END IF
         ELSE IF( C3.EQ.'TTR' ) THEN
            IF( SNAME ) THEN
               NB = 1
            ELSE
               NB = 1
            END IF
         ELSE IF( C3.EQ.'GST' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
            IF( TIME ) THEN
               GSTBLOCK = NB
            ELSE
               NB = GSTBLOCK
               IF( NB.LE.0 ) THEN
                  PRINT *, 'xpjlaenv.f ERROR common variable GSTBLOCK',
     $               ' may be unitialized'
c                 CALL EXIT( 13 )
                  STOP
               END IF
            END IF
         ELSE IF( C3.EQ.'BCK' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
            IF( TIME ) THEN
               BCKBLOCK = NB
            ELSE
               NB = BCKBLOCK
               IF( NB.LE.0 ) THEN
                  PRINT *, 'xpjlaenv.f ERROR common variable BCKBLOCK',
     $               ' may be unitialized'
c                 CALL EXIT( 13 )
                  STOP
               END IF
            END IF
         ELSE IF( C3.EQ.'TRS' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
            IF( TIME ) THEN
               TRSBLOCK = NB
            ELSE
               NB = TRSBLOCK
               IF( NB.LE.0 ) THEN
                  PRINT *, 'xpjlaenv.f ERROR common variable TRSBLOCK',
     $               ' may be unitialized'
c                 CALL EXIT( 13 )
                  STOP 
               END IF
            END IF
         END IF
      END IF
*
*
      PJLAENV = NB
      GLOBAL = .TRUE.
      GO TO 100
*
   60 CONTINUE
*
*     ISPEC = 2:  panel blocking factor (Used only in PxyyTTRD)
*     (local - different processes may use different values)
*
      NB = 16
      IF( C2.EQ.'SY' .OR. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TTR' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         END IF
      END IF
      IF( TIME ) THEN
         PNB = NB
      ELSE
         NB = PNB
         IF( NB.LE.0 ) THEN
            PRINT *, 'xpjlaenv.f ERROR common variable PNB',
     $         ' may be unitialized'
c           CALL EXIT( 13 )
            STOP
         END IF
      END IF
      PJLAENV = NB
      GLOBAL = .FALSE.
      GO TO 100
*
*
   70 CONTINUE
*
*     ISPEC = 3:  algorithmic blocking factor (Used only in PxyyTTRD)
*     (global - all processes must use the same value)
*
      NB = 16
      NB = 1
      IF( C2.EQ.'SY' .OR. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TTR' ) THEN
            IF( SNAME ) THEN
               NB = 16
            ELSE
               NB = 16
            END IF
         END IF
      END IF
      IF( TIME ) THEN
         ANB = NB
      ELSE
         NB = ANB
         IF( NB.LE.0 ) THEN
            PRINT *, 'xpjlaenv.f ERROR common variable ANB',
     $         ' may be unitialized'
c           CALL EXIT( 13 )
            STOP
         END IF
      END IF
      PJLAENV = NB
      GLOBAL = .TRUE.
      GO TO 100
*
   80 CONTINUE
*
*     ISPEC = 4:  Execution path options (Used only in PxyyTTRD)
*     (global - all processes must use the same value)
*
      PJLAENV = -4
      IF( C2.EQ.'SY' .OR. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TTR' ) THEN
*           V and H interleaved (default is not interleaved)
            IF( N1.EQ.1 ) THEN
               PJLAENV = 1
               IF( TIME ) THEN
                  INTERLEAVE = PJLAENV
               ELSE
                  PJLAENV = INTERLEAVE
               END IF
            END IF
*
*           Two ZGEMMs (default is one ZGEMM)
            IF( N1.EQ.2 ) THEN
               PJLAENV = 0
               IF( TIME ) THEN
                  TWOGEMMS = PJLAENV
               ELSE
                  PJLAENV = TWOGEMMS
               END IF
            END IF
*           Balanced Update (default is minimum communication update)
            IF( N1.EQ.3 ) THEN
               PJLAENV = 0
               IF( TIME ) THEN
                  BALANCED = PJLAENV
               ELSE
                  PJLAENV = BALANCED
               END IF
            END IF
         END IF
      END IF
      GLOBAL = .TRUE.
      GO TO 100
*
   90 CONTINUE
*
*     ISPEC = 5:  Minimum size to justify call to parallel code
*     (global - all processes must use the same value)
*
      MSZ = 0
      IF( C2.EQ.'SY' .OR. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TTR' ) THEN
            IF( SNAME ) THEN
               MSZ = 100
            ELSE
               MSZ = 100
            END IF
         END IF
      END IF
      IF( TIME ) THEN
         MINSZ = MSZ
      ELSE
         MSZ = MINSZ
      END IF
      PJLAENV = MSZ
      GLOBAL = .TRUE.
      GO TO 100
*
  100 CONTINUE
*
      IF( GLOBAL ) THEN
         IDUMM = 0
         CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, PJLAENV, 1, IDUMM,
     $                 IDUMM, -1, -1, IDUMM )
      END IF
*
*
*
      RETURN
*
*     End of PJLAENV
*
      END
