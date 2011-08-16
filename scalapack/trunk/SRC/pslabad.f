      SUBROUTINE PSLABAD( ICTXT, SMALL, LARGE )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT
      REAL               LARGE, SMALL
*     ..
*
*  Purpose
*  =======
*
*  PSLABAD takes as input the values computed by PSLAMCH for underflow
*  and overflow, and returns the square root of each of these values if
*  the log of LARGE is sufficiently large.  This subroutine is intended
*  to identify machines with a large exponent range, such as the Crays,
*  and redefine the underflow and overflow limits to be the square roots
*  of the values computed by PSLAMCH.  This subroutine is needed because
*  PSLAMCH does not compensate for poor arithmetic in the upper half of
*  the exponent range, as is found on a Cray.
*
*  In addition, this routine performs a global minimization and maximi-
*  zation on these values, to support heterogeneous computing networks.
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle in which the computation takes
*          place.
*
*  SMALL   (local input/local output) REAL
*          On entry, the underflow threshold as computed by PSLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of SMALL, otherwise unchanged.
*
*  LARGE   (local input/local output) REAL
*          On entry, the overflow threshold as computed by PSLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of LARGE, otherwise unchanged.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            IDUMM
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGAMN2D, SGAMX2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
*     ..
*     .. Executable Statements ..
*
*     If it looks like we're on a Cray, take the square root of
*     SMALL and LARGE to avoid overflow and underflow problems.
*
      IF( LOG10( LARGE ).GT.2000. ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
      IDUMM = 0
*
      CALL SGAMX2D( ICTXT, 'All', ' ', 1, 1, SMALL, 1, IDUMM,
     $              IDUMM, -1, -1, IDUMM )
      CALL SGAMN2D( ICTXT, 'All', ' ', 1, 1, LARGE, 1, IDUMM,
     $              IDUMM, -1, -1, IDUMM )
*
      RETURN
*
*     End of PSLABAD
*
      END
