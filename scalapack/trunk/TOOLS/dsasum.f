      DOUBLE PRECISION   FUNCTION DSASUM( N, X, INCX )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
*     ..
*     .. Array Arguments ..
      REAL               X( * )
*     ..
*
*  Purpose
*  =======
*
*  DSASUM is a simple FORTRAN wrapper around the BLAS function SASUM
*  returning the result as a double allowing it to be callable by C
*  programs.
*
*  =====================================================================
*
*     .. External Functions ..
      REAL               SASUM
      EXTERNAL           SASUM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
      DSASUM = DBLE( SASUM( N, X, INCX ) )
*
      RETURN
*
*     End of DSASUM
*
      END
