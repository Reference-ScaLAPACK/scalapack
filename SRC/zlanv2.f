      SUBROUTINE ZLANV2( A, B, C, D, RT1, RT2, CS, SN )
*
*  -- ScaLAPACK routine (version 1.7) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     May 28, 1999
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS
      COMPLEX*16         A, B, C, D, RT1, RT2, SN
*     ..
*
*  Purpose
*  =======
*
*  ZLANV2 computes the Schur factorization of a complex 2-by-2
*  nonhermitian matrix in standard form:
*
*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
*       [ C  D ]   [ SN  CS ] [  0  DD ] [-SN  CS ]
*
*  Arguments
*  =========
*
*  A       (input/output) COMPLEX*16
*  B       (input/output) COMPLEX*16
*  C       (input/output) COMPLEX*16
*  D       (input/output) COMPLEX*16
*          On entry, the elements of the input matrix.
*          On exit, they are overwritten by the elements of the
*          standardised Schur form.
*
*  RT1     (output) COMPLEX*16
*  RT2     (output) COMPLEX*16
*          The two eigenvalues.
*
*  CS      (output) DOUBLE PRECISION
*  SN      (output) COMPLEX*16
*          Parameters of the rotation matrix.
*
*  Further Details
*  ===============
*
*  Implemented by Mark R. Fahey, May 28, 1999
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   RZERO, HALF, RONE
      PARAMETER          ( RZERO = 0.0D+0, HALF = 0.5D+0,
     $                   RONE = 1.0D+0 )
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         AA, BB, DD, T, TEMP, TEMP2, U, X, Y
*     ..
*     .. External Functions ..
      COMPLEX*16         ZLADIV
      EXTERNAL           ZLADIV
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLARTG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DCONJG, DIMAG, SQRT
*     ..
*     .. Executable Statements ..
*
*     Initialize CS and SN
*
      CS = RONE
      SN = ZERO
*
      IF( C.EQ.ZERO ) THEN
         GO TO 10
*
      ELSE IF( B.EQ.ZERO ) THEN
*
*        Swap rows and columns
*
         CS = RZERO
         SN = ONE
         TEMP = D
         D = A
         A = TEMP
         B = -C
         C = ZERO
         GO TO 10
      ELSE IF( ( A-D ).EQ.ZERO ) THEN
         TEMP = SQRT( B*C )
         A = A + TEMP
         D = D - TEMP
         IF( ( B+C ).EQ.ZERO ) THEN
            CS = SQRT( HALF )
            SN = DCMPLX( RZERO, RONE )*CS
         ELSE
            TEMP = SQRT( B+C )
            TEMP2 = ZLADIV( SQRT( B ), TEMP )
            CS = DBLE( TEMP2 )
            SN = ZLADIV( SQRT( C ), TEMP )
         END IF
         B = B - C
         C = ZERO
         GO TO 10
      ELSE
*
*        Compute eigenvalue closest to D
*
         T = D
         U = B*C
         X = HALF*( A-T )
         Y = SQRT( X*X+U )
         IF( DBLE( X )*DBLE( Y )+DIMAG( X )*DIMAG( Y ).LT.RZERO )
     $      Y = -Y
         T = T - ZLADIV( U, ( X+Y ) )
*
*        Do one QR step with exact shift T - resulting 2 x 2 in
*        triangular form.
*
         CALL ZLARTG( A-T, C, CS, SN, AA )
*
         D = D - T
         BB = CS*B + SN*D
         DD = -DCONJG( SN )*B + CS*D
*
         A = AA*CS + BB*DCONJG( SN ) + T
         B = -AA*SN + BB*CS
         C = ZERO
         D = T
*
      END IF
*
   10 CONTINUE
*
*     Store eigenvalues in RT1 and RT2.
*
      RT1 = A
      RT2 = D
      RETURN
*
*     End of ZLANV2
*
      END
