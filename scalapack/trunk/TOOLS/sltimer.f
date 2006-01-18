      SUBROUTINE SLBOOT()
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*  Purpose
*  =======
*
*  SLBOOT (re)sets all timers to 0, and enables SLtimer.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTIMER
      PARAMETER          ( NTIMER = 64 )
      DOUBLE PRECISION   STARTFLAG, ZERO
      PARAMETER          ( STARTFLAG = -5.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. Common Blocks ..
      LOGICAL            DISABLED
      DOUBLE PRECISION   CPUSEC( NTIMER ), CPUSTART( NTIMER ),
     $                   WALLSEC( NTIMER ), WALLSTART( NTIMER )
      COMMON /SLTIMER00/ CPUSEC, WALLSEC, CPUSTART, WALLSTART, DISABLED
*     ..
*     .. Executable Statements ..
*
      DISABLED = .FALSE.
      DO 10 I = 1, NTIMER
         CPUSEC( I )  = ZERO
         WALLSEC( I ) = ZERO
         CPUSTART( I )  = STARTFLAG
         WALLSTART( I ) = STARTFLAG
   10 CONTINUE
*
      RETURN
*
*     End of SLBOOT
*
      END
*
      SUBROUTINE SLTIMER( I )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            I
*     ..
*
*  Purpose
*  =======
*
*  SLtimer provides a "stopwatch" functionality cpu/wall timer
*  (in seconds).  Up to 64 separate timers can be functioning at once.
*  The first call starts the timer, and the second stops it.  This
*  routine can be disenabled, so that calls to the timer are ignored.
*  This feature can be used to make sure certain sections of code do
*  not affect timings, even if they call routines which have SLtimer
*  calls in them.
*
*  Arguments
*  =========
*
*  I       (global input) INTEGER
*          The timer to stop/start.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTIMER
      PARAMETER          ( NTIMER = 64 )
      DOUBLE PRECISION   STARTFLAG
      PARAMETER          ( STARTFLAG = -5.0D+0 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DCPUTIME00, DWALLTIME00
      EXTERNAL           DCPUTIME00, DWALLTIME00
*     ..
*     .. Common Blocks ..
      LOGICAL            DISABLED
      DOUBLE PRECISION   CPUSEC( NTIMER ), CPUSTART( NTIMER ),
     $                   WALLSEC( NTIMER ), WALLSTART( NTIMER )
      COMMON /SLTIMER00/ CPUSEC, WALLSEC, CPUSTART, WALLSTART, DISABLED
*     ..
*     .. Executable Statements ..
*
*     If timing disabled, return
*
      IF( DISABLED )
     $   RETURN
*
      IF( WALLSTART( I ).EQ.STARTFLAG ) THEN
*
*        If timer has not been started, start it
*
         WALLSTART( I ) = DWALLTIME00()
         CPUSTART( I )  = DCPUTIME00()
*
      ELSE
*
*        Stop timer and add interval to count
*
         CPUSEC( I ) = CPUSEC( I ) + DCPUTIME00() - CPUSTART( I )
         WALLSEC( I ) = WALLSEC( I ) + DWALLTIME00() - WALLSTART( I )
         WALLSTART( I ) = STARTFLAG
*
      END IF
*
      RETURN
*
*     End of SLTIMER
*
      END
*
      SUBROUTINE SLENABLE()
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*  Purpose
*  =======
*
*  SLENABLE sets it so calls to SLtimer are not ignored.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTIMER
      PARAMETER          ( NTIMER = 64 )
*     ..
*     .. Common Blocks ..
      LOGICAL            DISABLED
      DOUBLE PRECISION   CPUSEC( NTIMER ), CPUSTART( NTIMER ),
     $                   WALLSEC( NTIMER ), WALLSTART( NTIMER )
      COMMON /SLTIMER00/ CPUSEC, WALLSEC, CPUSTART, WALLSTART, DISABLED
*     ..
*     .. Executable Statements ..
*
      DISABLED = .FALSE.
*
      RETURN
*
      END
*
      SUBROUTINE SLDISABLE()
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*  Purpose
*  =======
*
*  SLDISABLE sets it so calls to SLTIMER are ignored.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTIMER
      PARAMETER          ( NTIMER = 64 )
*     ..
*     .. Common Blocks ..
      LOGICAL            DISABLED
      DOUBLE PRECISION   CPUSEC( NTIMER ), CPUSTART( NTIMER ),
     $                   WALLSEC( NTIMER ), WALLSTART( NTIMER )
      COMMON /SLTIMER00/ CPUSEC, WALLSEC, CPUSTART, WALLSTART, DISABLED
*     ..
*     .. Executable Statements ..
*
      DISABLED = .TRUE.
*
      RETURN
*
*     End of SLDISABLE
*
      END
*
      DOUBLE PRECISION FUNCTION SLINQUIRE( TIMETYPE, I )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER*1        TIMETYPE
      INTEGER            I
*     ..
*
*  Purpose
*  =======
*
*  SLINQUIRE returns wall or cpu time that has accumulated in timer I.
*
*  Arguments
*  =========
*
*  TIMETYPE (global input) CHARACTER
*           Controls what time will be returned:
*           = 'W': wall clock time is returned,
*           = 'C': CPU time is returned (default).
*
*  I        (global input) INTEGER
*           The timer to return.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTIMER
      PARAMETER          ( NTIMER = 64 )
      DOUBLE PRECISION   ERRFLAG
      PARAMETER          ( ERRFLAG = -1.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   TIME
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DCPUTIME00, DWALLTIME00
      EXTERNAL           DCPUTIME00, DWALLTIME00, LSAME
*     ..
*     .. Common Blocks ..
      LOGICAL            DISABLED
      DOUBLE PRECISION   CPUSEC( NTIMER ), CPUSTART( NTIMER ),
     $                   WALLSEC( NTIMER ), WALLSTART( NTIMER )
      COMMON /SLTIMER00/ CPUSEC, WALLSEC, CPUSTART, WALLSTART, DISABLED
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( TIMETYPE, 'W' ) ) THEN
*
*        If walltime not available on this machine, return -1 flag
*
         IF( DWALLTIME00().EQ.ERRFLAG ) THEN
            TIME = ERRFLAG
         ELSE
            TIME = WALLSEC( I )
         END IF
      ELSE
         IF( DCPUTIME00().EQ.ERRFLAG ) THEN
            TIME = ERRFLAG
         ELSE
            TIME = CPUSEC( I )
         END IF
      END IF
*
      SLINQUIRE = TIME
*
      RETURN
*
*     End of SLINQUIRE
*
      END
*
      SUBROUTINE SLCOMBINE( ICTXT, SCOPE, OP, TIMETYPE, N, IBEG,
     $                      TIMES )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          OP, SCOPE, TIMETYPE
      INTEGER            IBEG, ICTXT, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   TIMES( N )
*     ..
*
*  Purpose
*  =======
*
*  SLCOMBINE takes the timing information stored on a scope of processes
*  and combines them into the user's TIMES array.
*
*  Arguments
*  =========
*
*  ICTXT    (local input) INTEGER
*           The BLACS context handle.
*
*  SCOPE    (global input) CHARACTER
*           Controls what processes in grid participate in combine.
*           Options are 'Rowwise', 'Columnwise', or 'All'.
*
*  OP       (global input) CHARACTER
*           Controls what combine should be done:
*           = '>': get maximal time on any process (default),
*           = '<': get minimal time on any process,
*           = '+': get sum of times across processes.
*
*  TIMETYPE (global input) CHARACTER
*           Controls what time will be returned in TIMES:
*           = 'W': wall clock time,
*           = 'C': CPU time (default).
*
*  N        (global input) INTEGER
*           The number of timers to combine.
*
*  IBEG     (global input) INTEGER
*           The first timer to be combined.
*
*  TIMES    (global output) DOUBLE PRECISION array, dimension (N)
*           The requested timing information is returned in this array.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTIMER
      PARAMETER          ( NTIMER = 64 )
      DOUBLE PRECISION   ERRFLAG
      PARAMETER          ( ERRFLAG = -1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            TMPDIS
      INTEGER            I
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGAMX2D, DGAMN2D, DGSUM2D
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DCPUTIME00, DWALLTIME00
      EXTERNAL           DCPUTIME00, DWALLTIME00, LSAME
*     ..
*     .. Common Blocks ..
      LOGICAL            DISABLED
      DOUBLE PRECISION   CPUSEC( NTIMER ), CPUSTART( NTIMER ),
     $                   WALLSEC( NTIMER ), WALLSTART( NTIMER )
      COMMON /SLTIMER00/ CPUSEC, WALLSEC, CPUSTART, WALLSTART, DISABLED
*     ..
*     .. Executable Statements ..
*
*     Disable timer for combine operation
*
      TMPDIS = DISABLED
      DISABLED = .TRUE.
*
*     Copy timer information into user's times array
*
      IF( LSAME( TIMETYPE, 'W' ) ) THEN
*
*        If walltime not available on this machine, fill in times
*        with -1 flag, and return
*
         IF( DWALLTIME00().EQ.ERRFLAG ) THEN
            DO 10 I = 1, N
               TIMES( I ) = ERRFLAG
   10       CONTINUE
            RETURN
         ELSE
            DO 20 I = 1, N
               TIMES( I ) = WALLSEC( IBEG + I - 1 )
   20       CONTINUE
         END IF
      ELSE
         IF( DCPUTIME00().EQ.ERRFLAG ) THEN
            DO 30 I = 1, N
               TIMES( I ) = ERRFLAG
   30       CONTINUE
            RETURN
         ELSE
            DO 40 I = 1, N
               TIMES( I ) = CPUSEC( IBEG + I - 1 )
   40       CONTINUE
         END IF
      ENDIF
*
*     Combine all nodes' information, restore disabled, and return
*
      IF( OP.EQ.'>' ) THEN
         CALL DGAMX2D( ICTXT, SCOPE, ' ', N, 1, TIMES, N, -1, -1,
     $                 -1, -1, 0 )
      ELSE IF( OP.EQ.'<' ) THEN
         CALL DGAMN2D( ICTXT, SCOPE, ' ', N, 1, TIMES, N, -1, -1,
     $                 -1, -1, 0 )
      ELSE IF( OP.EQ.'+' ) THEN
         CALL DGSUM2D( ICTXT, SCOPE, ' ', N, 1, TIMES, N, -1, 0 )
      ELSE
         CALL DGAMX2D( ICTXT, SCOPE, ' ', N, 1, TIMES, N, -1, -1,
     $                 -1, -1, 0 )
      END IF
*
      DISABLED = TMPDIS
*
      RETURN
*
*     End of SLCOMBINE
*
      END
