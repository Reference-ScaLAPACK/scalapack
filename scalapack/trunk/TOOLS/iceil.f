      INTEGER FUNCTION ICEIL( INUM, IDENOM )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IDENOM, INUM
*     ..
*
*  Purpose
*  =======
*
*  ICEIL returns the ceiling of the division of two integers.
*
*  Arguments
*  =========
*
*  INUM     (local input) INTEGER
*           The numerator,
*
*  IDENOM   (local input) INTEGER
*           and the denominator of the fraction to be evaluated.
*
*  =====================================================================
*
*     .. Executable Statements ..
*
      ICEIL = (INUM+IDENOM-1) / IDENOM
*
      RETURN
*
*     End of ICEIL
*
      END
