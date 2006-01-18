      SUBROUTINE SSDOT( N, DOT, X, INCX, Y, INCY )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      REAL               DOT
*     ..
*     .. Array Arguments ..
      REAL               X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  SSDOT is a simple FORTRAN wrapper around the BLAS function
*  SDOT returning the result in the parameter list instead.
*
*  =====================================================================
*
*     .. External Functions ..
      REAL               SDOT
      EXTERNAL           SDOT
*     ..
*     .. Executable Statements ..
*
      DOT = SDOT( N, X, INCX, Y, INCY )
*
      RETURN
*
*     End of SSDOT
*
      END
