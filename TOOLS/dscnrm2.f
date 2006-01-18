      DOUBLE PRECISION   FUNCTION DSCNRM2( N, X, INCX )
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
      COMPLEX            X( * )
*     ..
*
*  Purpose
*  =======
*
*  DSCNRM2 is a simple FORTRAN wrapper around the BLAS function SCNRM2
*  returning the result as a double allowing it to be callable by C
*  programs.
*
*  =====================================================================
*
*     .. External Functions ..
      REAL               SCNRM2
      EXTERNAL           SCNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
      DSCNRM2 = DBLE( SCNRM2( N, X, INCX ) )
*
      RETURN
*
*     End of DSCNRM2
*
      END
