      SUBROUTINE DVVDOT( N, DOT, X, INCX, Y, INCY )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
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
*  DVVDOT computes the following dot product:
*
*     dot = dot + x**T * y,
*
*  where x and y are n vectors.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          On entry, N  specifies the length of the  vectors x and y.  N
*          must be at least zero.
*
*  DOT     (input/output) DOUBLE PRECISION
*          On exit, DOT is updated with the dot product of the vectors x
*          and y.
*
*  X       (input) DOUBLE PRECISION array of dimension at least
*          ( 1 + ( n - 1 )*abs( INCX ) ). Before entry,  the incremented
*          array X must contain the vector x.
*
*  INCX    (input) INTEGER
*          On entry, INCX specifies the increment for the elements of X.
*          INCX must not be zero.
*
*  Y       (input) DOUBLE PRECISION array of dimension at least
*          ( 1 + ( n - 1 )*abs( INCY ) ). Before entry,  the incremented
*          array Y must contain the vector y.
*
*  INCY    (input) INTEGER
*          On entry, INCY specifies the increment for the elements of Y.
*          INCY must not be zero.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. External Functions ..
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT
*     ..
*     .. Executable Statements ..
*
      DOT = DOT + DDOT( N, X, INCX, Y, INCY )
*
      RETURN
*
*     End of DVVDOT
*
      END
