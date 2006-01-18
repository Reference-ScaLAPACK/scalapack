      SUBROUTINE ZZDOTU( N, DOTU, X, INCX, Y, INCY )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      COMPLEX*16         DOTU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZZDOTU is a simple FORTRAN wrapper around the BLAS function
*  ZDOTU returning the result in the parameter list instead.
*
*  =====================================================================
*
*     .. External Functions ..
      COMPLEX*16         ZDOTU
      EXTERNAL           ZDOTU
*     ..
*     .. Executable Statements ..
*
      DOTU = ZDOTU( N, X, INCX, Y, INCY )
*
      RETURN
*
*     End of ZZDOTU
*
      END
