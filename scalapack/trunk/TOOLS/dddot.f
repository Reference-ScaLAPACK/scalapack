      SUBROUTINE DDDOT( N, DOT, X, INCX, Y, INCY )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   DOT
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DDDOT is a simple FORTRAN wrapper around the BLAS function
*  DDOT returning the result in the parameter list instead.
*
*  =====================================================================
*
*     .. External Functions ..
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT
*     ..
*     .. Executable Statements ..
*
      DOT = DDOT( N, X, INCX, Y, INCY )
*
      RETURN
*
*     End of DDDOT
*
      END
