      DOUBLE PRECISION   FUNCTION DSCASUM( N, X, INCX )
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
      COMPLEX            X( * )
*     ..
*
*  Purpose
*  =======
*
*  DSCASUM is a simple FORTRAN wrapper around the BLAS function SCASUM
*  returning the result as a double allowing it to be callable by C
*  programs.
*
*  =====================================================================
*
*     .. External Functions ..
      REAL               SCASUM
      EXTERNAL           SCASUM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
      DSCASUM = DBLE( SCASUM( N, X, INCX ) )
*
      RETURN
*
*     End of DSCASUM
*
      END
