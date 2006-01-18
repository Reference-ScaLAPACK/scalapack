      DOUBLE PRECISION   FUNCTION DSNRM2( N, X, INCX )
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
      REAL               X( * )
*     ..
*
*  Purpose
*  =======
*
*  DSNRM2 is a simple FORTRAN wrapper around the BLAS function SNRM2
*  returning the result as a double allowing it to be callable by C
*  programs.
*
*  =====================================================================
*
*     .. External Functions ..
      REAL               SNRM2
      EXTERNAL           SNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
      DSNRM2 = DBLE( SNRM2( N, X, INCX ) )
*
      RETURN
*
*     End of DSNRM2
*
      END
