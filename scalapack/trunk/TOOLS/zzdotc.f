      SUBROUTINE ZZDOTC( N, DOTC, X, INCX, Y, INCY )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      COMPLEX*16         DOTC
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZZDOTC is a simple FORTRAN wrapper around the BLAS function
*  ZDOTC returning the result in the parameter list instead.
*
*  =====================================================================
*
*     .. External Functions ..
      COMPLEX*16         ZDOTC
      EXTERNAL           ZDOTC
*     ..
*     .. Executable Statements ..
*
      DOTC = ZDOTC( N, X, INCX, Y, INCY )
*
      RETURN
*
*     End of ZZDOTC
*
      END
