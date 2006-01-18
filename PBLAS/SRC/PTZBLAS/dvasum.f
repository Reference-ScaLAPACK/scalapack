      SUBROUTINE DVASUM( N, ASUM, X, INCX )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ASUM
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DVASUM  returns the sum of absolute values of the entries of a vector
*  x.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          On entry, N specifies the length of the vector x. N  must  be
*          at least zero.
*
*  ASUM    (output) DOUBLE PRECISION
*          On exit, ASUM specifies the sum of absolute values.
*
*  X       (input) DOUBLE PRECISION array of dimension at least
*          ( 1 + ( n - 1 )*abs( INCX ) ). Before entry,  the incremented
*          array X must contain the vector x.
*
*  INCX    (input) INTEGER
*          On entry, INCX specifies the increment for the elements of X.
*          INCX must not be zero.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. External Functions ..
      DOUBLE PRECISION   DASUM
      EXTERNAL           DASUM
*     ..
*     .. Executable Statements ..
*
      ASUM = DASUM( N, X, INCX )
*
      RETURN
*
*     End of DVASUM
*
      END
