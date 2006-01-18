      SUBROUTINE DASQRTB( A, B, C )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C
*     ..
*
*  Purpose
*  =======
*
*  DASQRTB computes c := a * sqrt( b ) where a, b and c are scalars.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*          On entry, A specifies the scalar a.
*
*  B       (input) DOUBLE PRECISION
*          On entry, B specifies the scalar b.
*
*  C       (output) DOUBLE PRECISION
*          On entry, C specifies the scalar c. On exit, c is overwritten
*          by the product of a and the square root of b.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          SQRT
*     ..
*     .. Executable Statements ..
*
      C = A * SQRT( B )
*
      RETURN
*
*     End of DASQRTB
*
      END
