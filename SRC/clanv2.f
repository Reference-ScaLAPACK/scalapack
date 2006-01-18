      SUBROUTINE CLANV2( A, B, C, D, RT1, RT2, CS, SN )
*
*  -- ScaLAPACK routine (version 1.7) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     May 28, 1999
*
*     .. Scalar Arguments ..
      REAL               CS
      COMPLEX            A, B, C, D, RT1, RT2, SN
*     ..
*
*  Purpose
*  =======
*
*  CLANV2 computes the Schur factorization of a complex 2-by-2
*  nonhermitian matrix in standard form:
*
*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
*       [ C  D ]   [ SN  CS ] [  0  DD ] [-SN  CS ]
*
*  Arguments
*  =========
*
*  A       (input/output) COMPLEX
*  B       (input/output) COMPLEX
*  C       (input/output) COMPLEX
*  D       (input/output) COMPLEX
*          On entry, the elements of the input matrix.
*          On exit, they are overwritten by the elements of the
*          standardised Schur form.
*
*  RT1     (output) COMPLEX
*  RT2     (output) COMPLEX
*          The two eigenvalues.
*
*  CS      (output) REAL
*  SN      (output) COMPLEX
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
      REAL               RZERO, HALF, RONE
      PARAMETER          ( RZERO = 0.0E+0, HALF = 0.5E+0,
     $                   RONE = 1.0E+0 )
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ),
     $                   ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX            AA, BB, DD, T, TEMP, TEMP2, U, X, Y
*     ..
*     .. External Functions ..
      COMPLEX            CLADIV
      EXTERNAL           CLADIV
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARTG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, CMPLX, CONJG, AIMAG, SQRT
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
            SN = CMPLX( RZERO, RONE )*CS
         ELSE
            TEMP = SQRT( B+C )
            TEMP2 = CLADIV( SQRT( B ), TEMP )
            CS = REAL( TEMP2 )
            SN = CLADIV( SQRT( C ), TEMP )
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
         IF( REAL( X )*REAL( Y )+AIMAG( X )*AIMAG( Y ).LT.RZERO )
     $      Y = -Y
         T = T - CLADIV( U, ( X+Y ) )
*
*        Do one QR step with exact shift T - resulting 2 x 2 in
*        triangular form.
*
         CALL CLARTG( A-T, C, CS, SN, AA )
*
         D = D - T
         BB = CS*B + SN*D
         DD = -CONJG( SN )*B + CS*D
*
         A = AA*CS + BB*CONJG( SN ) + T
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
*     End of CLANV2
*
      END
