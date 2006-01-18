      SUBROUTINE CCDOTU( N, DOTU, X, INCX, Y, INCY )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      COMPLEX            DOTU
*     ..
*     .. Array Arguments ..
      COMPLEX            X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  CCDOTU is a simple FORTRAN wrapper around the BLAS function
*  CDOTU returning the result in the parameter list instead.
*
*  =====================================================================
*
*     .. External Functions ..
      COMPLEX            CDOTU
      EXTERNAL           CDOTU
*     ..
*     .. Executable Statements ..
*
      DOTU = CDOTU( N, X, INCX, Y, INCY )
*
      RETURN
*
*     End of CCDOTU
*
      END
