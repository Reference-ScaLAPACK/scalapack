      SUBROUTINE SASQRTB( A, B, C )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      REAL               A, B, C
*     ..
*
*  Purpose
*  =======
*
*  SASQRTB computes c := a * sqrt( b ) where a, b and c are scalars.
*
*  Arguments
*  =========
*
*  A       (input) REAL
*          On entry, A specifies the scalar a.
*
*  B       (input) REAL
*          On entry, B specifies the scalar b.
*
*  C       (output) REAL
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
*     End of SASQRTB
*
      END
