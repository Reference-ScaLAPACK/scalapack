      FUNCTION DISANAN( DIN1, DIN2 )
      IMPLICIT NONE
      LOGICAL DISANAN
*
*  -- LAPACK auxiliary routine (version *TBA*) --
*     July 4, 2010
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DIN1, DIN2
*     ..
*
*  Purpose
*  =======
*
*  DISANAN checks for NaNs by comparing its two arguments for
*  inequality.  NaN is the only floating-point value where NaN.ne.NaN
*  returns .TRUE.  To check for NaNs, pass the same variable as both
*  arguments.
*
*  The function is only called with the two input arguments being
*  identical. The Fortran compiler however must assume that the two 
*  arguments are not the same variable, and the test will not be 
*  optimized away. Warning! Interprocedural or whole-program optimization 
*  might defeat this test. This has not been oberved so far though.  
*
*  =====================================================================
*
*  .. Executable Statements ..

      DISANAN = (DIN1.NE.DIN2)
*
      END
