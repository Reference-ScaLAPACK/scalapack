      SUBROUTINE BSLAEXC( N, T, LDT, J1, N1, N2, ITRAF, DTRAF, WORK,
     $                    INFO )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      INTEGER            INFO, J1, LDT, N, N1, N2
*     ..
*     .. Array Arguments ..
      INTEGER            ITRAF( * )
      REAL               DTRAF( * ), T( LDT, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  BSLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
*  an upper quasi-triangular matrix T by an orthogonal similarity
*  transformation.
*
*  In contrast to the LAPACK routine DLAEXC, the orthogonal
*  transformation matrix Q is not explicitly constructed but
*  represented by paramaters contained in the arrays ITRAF and DTRAF,
*  see the description of BSTREXC for more details.
*
*  T must be in Schur canonical form, that is, block upper triangular
*  with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
*  has its diagonal elemnts equal and its off-diagonal elements of
*  opposite sign.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input/output) REAL             array, dimension (LDT,N)
*          On entry, the upper quasi-triangular matrix T, in Schur
*          canonical form.
*          On exit, the updated matrix T, again in Schur canonical form.
*
*  LDT     (input)  INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  J1      (input) INTEGER
*          The index of the first row of the first block T11.
*
*  N1      (input) INTEGER
*          The order of the first block T11. N1 = 0, 1 or 2.
*
*  N2      (input) INTEGER
*          The order of the second block T22. N2 = 0, 1 or 2.
*
*  ITRAF   (output) INTEGER array, length k, where
*             k = 1, if N1+N2 = 2;
*             k = 2, if N1+N2 = 3;
*             k = 4, if N1+N2 = 4.
*          List of parameters for representing the transformation
*          matrix Q, see BSTREXC.
*
*  DTRAF   (output) REAL             array, length k, where
*             k =  2, if N1+N2 = 2;
*             k =  5, if N1+N2 = 3;
*             k = 10, if N1+N2 = 4.
*          List of parameters for representing the transformation
*          matrix Q, see BSTREXC.
*
*  WORK    (workspace) REAL             array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          = 1: the transformed matrix T would be too far from Schur
*               form; the blocks are not swapped and T and Q are
*               unchanged.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      REAL               TEN
      PARAMETER          ( TEN = 10.0 )
      INTEGER            LDD, LDX
      PARAMETER          ( LDD = 4, LDX = 2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IERR, J2, J3, J4, K, LD, LI, ND
      REAL               CS, DNORM, EPS, SCALE, SMLNUM, SN, T11, T22,
     $                   T33, TAU, TAU1, TAU2, TEMP, THRESH, WI1, WI2,
     $                   WR1, WR2, XNORM
*     ..
*     .. Local Arrays ..
      REAL               D( LDD, 4 ), X( LDX, 2 )
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAMOV, SLANV2, SLARFG, SLARFX, SLARTG, SLASY2,
     $                   SROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. N1.EQ.0 .OR. N2.EQ.0 )
     $   RETURN
      IF( J1+N1.GT.N )
     $   RETURN
*
      J2 = J1 + 1
      J3 = J1 + 2
      J4 = J1 + 3
*
      IF( N1.EQ.1 .AND. N2.EQ.1 ) THEN
*
*        Swap two 1-by-1 blocks.
*
         T11 = T( J1, J1 )
         T22 = T( J2, J2 )
*
*        Determine the transformation to perform the interchange.
*
         CALL SLARTG( T( J1, J2 ), T22-T11, CS, SN, TEMP )
*
*        Apply transformation to the matrix T.
*
         IF( J3.LE.N )
     $      CALL SROT( N-J1-1, T( J1, J3 ), LDT, T( J2, J3 ), LDT, CS,
     $                 SN )
         CALL SROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
*
         T( J1, J1 ) = T22
         T( J2, J2 ) = T11
*
         ITRAF( 1 ) = J1
         DTRAF( 1 ) = CS
         DTRAF( 2 ) = SN
*
      ELSE
*
*        Swapping involves at least one 2-by-2 block.
*
*        Copy the diagonal block of order N1+N2 to the local array D
*        and compute its norm.
*
         ND = N1 + N2
         CALL SLAMOV( 'Full', ND, ND, T( J1, J1 ), LDT, D, LDD )
         DNORM = SLANGE( 'Max', ND, ND, D, LDD, WORK )
*
*        Compute machine-dependent threshold for test for accepting
*        swap.
*
         EPS = SLAMCH( 'P' )
         SMLNUM = SLAMCH( 'S' ) / EPS
         THRESH = MAX( TEN*EPS*DNORM, SMLNUM )
*
*        Solve T11*X - X*T22 = scale*T12 for X.
*
         CALL SLASY2( .FALSE., .FALSE., -1, N1, N2, D, LDD,
     $                D( N1+1, N1+1 ), LDD, D( 1, N1+1 ), LDD, SCALE, X,
     $                LDX, XNORM, IERR )
*
*        Swap the adjacent diagonal blocks.
*
         K = N1 + N1 + N2 - 3
         GO TO ( 10, 20, 30 )K
*
   10    CONTINUE
*
*        N1 = 1, N2 = 2: generate elementary reflector H so that:
*
*        ( scale, X11, X12 ) H = ( 0, 0, * )
*
         DTRAF( 1 ) = SCALE
         DTRAF( 2 ) = X( 1, 1 )
         DTRAF( 3 ) = X( 1, 2 )
         CALL SLARFG( 3, DTRAF( 3 ), DTRAF, 1, TAU )
         DTRAF( 3 ) = ONE
         T11 = T( J1, J1 )
*
*        Perform swap provisionally on diagonal block in D.
*
         CALL SLARFX( 'Left', 3, 3, DTRAF, TAU, D, LDD, WORK )
         CALL SLARFX( 'Right', 3, 3, DTRAF, TAU, D, LDD, WORK )
*
*        Test whether to reject swap.
*
         IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 3,
     $       3 )-T11 ) ).GT.THRESH )GO TO 50
*
*        Accept swap: apply transformation to the entire matrix T.
*
         CALL SLARFX( 'Left', 3, N-J1+1, DTRAF, TAU, T( J1, J1 ), LDT,
     $                WORK )
         CALL SLARFX( 'Right', J2, 3, DTRAF, TAU, T( 1, J1 ), LDT,
     $                WORK )
*
         T( J3, J1 ) = ZERO
         T( J3, J2 ) = ZERO
         T( J3, J3 ) = T11
*
         ITRAF( 1 ) = 2*N + J1
         LI = 2
         DTRAF( 3 ) = TAU
         LD = 4
         GO TO 40
*
   20    CONTINUE
*
*        N1 = 2, N2 = 1: generate elementary reflector H so that:
*
*        H (  -X11 ) = ( * )
*          (  -X21 ) = ( 0 )
*          ( scale ) = ( 0 )
*
         DTRAF( 1 ) = -X( 1, 1 )
         DTRAF( 2 ) = -X( 2, 1 )
         DTRAF( 3 ) = SCALE
         CALL SLARFG( 3, DTRAF( 1 ), DTRAF( 2 ), 1, TAU )
         DTRAF( 1 ) = ONE
         T33 = T( J3, J3 )
*
*        Perform swap provisionally on diagonal block in D.
*
         CALL SLARFX( 'Left', 3, 3, DTRAF, TAU, D, LDD, WORK )
         CALL SLARFX( 'Right', 3, 3, DTRAF, TAU, D, LDD, WORK )
*
*        Test whether to reject swap.
*
         IF( MAX( ABS( D( 2, 1 ) ), ABS( D( 3, 1 ) ), ABS( D( 1,
     $       1 )-T33 ) ).GT.THRESH )GO TO 50
*
*        Accept swap: apply transformation to the entire matrix T.
*
         CALL SLARFX( 'Right', J3, 3, DTRAF, TAU, T( 1, J1 ), LDT,
     $                WORK )
         CALL SLARFX( 'Left', 3, N-J1, DTRAF, TAU, T( J1, J2 ), LDT,
     $                WORK )
*
         T( J1, J1 ) = T33
         T( J2, J1 ) = ZERO
         T( J3, J1 ) = ZERO
*
         ITRAF( 1 ) = N + J1
         LI = 2
         DTRAF( 1 ) = TAU
         LD = 4
         GO TO 40
*
   30    CONTINUE
*
*        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
*        that:
*
*        H(2) H(1) (  -X11  -X12 ) = (  *  * )
*                  (  -X21  -X22 )   (  0  * )
*                  ( scale    0  )   (  0  0 )
*                  (    0  scale )   (  0  0 )
*
         DTRAF( 1 ) = -X( 1, 1 )
         DTRAF( 2 ) = -X( 2, 1 )
         DTRAF( 3 ) = SCALE
         CALL SLARFG( 3, DTRAF( 1 ), DTRAF( 2 ), 1, TAU1 )
         DTRAF( 1 ) = ONE
*
         TEMP = -TAU1*( X( 1, 2 )+DTRAF( 2 )*X( 2, 2 ) )
         DTRAF( 4 ) = -TEMP*DTRAF( 2 ) - X( 2, 2 )
         DTRAF( 5 ) = -TEMP*DTRAF( 3 )
         DTRAF( 6 ) = SCALE
         CALL SLARFG( 3, DTRAF( 4 ), DTRAF( 5 ), 1, TAU2 )
         DTRAF( 4 ) = ONE
*
*        Perform swap provisionally on diagonal block in D.
*
         CALL SLARFX( 'Left', 3, 4, DTRAF, TAU1, D, LDD, WORK )
         CALL SLARFX( 'Right', 4, 3, DTRAF, TAU1, D, LDD, WORK )
         CALL SLARFX( 'Left', 3, 4, DTRAF( 4 ), TAU2, D( 2, 1 ), LDD,
     $                WORK )
         CALL SLARFX( 'Right', 4, 3, DTRAF( 4 ), TAU2, D( 1, 2 ), LDD,
     $                WORK )
*
*        Test whether to reject swap.
*
         IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 4, 1 ) ),
     $       ABS( D( 4, 2 ) ) ).GT.THRESH )GO TO 50
*
*        Accept swap: apply transformation to the entire matrix T.
*
         CALL SLARFX( 'Left', 3, N-J1+1, DTRAF, TAU1, T( J1, J1 ), LDT,
     $                WORK )
         CALL SLARFX( 'Right', J4, 3, DTRAF, TAU1, T( 1, J1 ), LDT,
     $                WORK )
         CALL SLARFX( 'Left', 3, N-J1+1, DTRAF( 4 ), TAU2, T( J2, J1 ),
     $                LDT, WORK )
         CALL SLARFX( 'Right', J4, 3, DTRAF( 4 ), TAU2, T( 1, J2 ), LDT,
     $                WORK )
*
         T( J3, J1 ) = ZERO
         T( J3, J2 ) = ZERO
         T( J4, J1 ) = ZERO
         T( J4, J2 ) = ZERO
*
         ITRAF( 1 ) = N + J1
         ITRAF( 2 ) = N + J2
         LI = 3
         DTRAF( 1 ) = TAU1
         DTRAF( 4 ) = TAU2
         LD = 7
         GO TO 40
*
   40    CONTINUE
*
         IF( N2.EQ.2 ) THEN
*
*           Standardize new 2-by-2 block T11
*
            CALL SLANV2( T( J1, J1 ), T( J1, J2 ), T( J2, J1 ),
     $                   T( J2, J2 ), WR1, WI1, WR2, WI2, CS, SN )
            CALL SROT( N-J1-1, T( J1, J1+2 ), LDT, T( J2, J1+2 ), LDT,
     $                 CS, SN )
            CALL SROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
            ITRAF( LI ) = J1
            LI = LI + 1
            DTRAF( LD ) = CS
            DTRAF( LD+1 ) = SN
            LD = LD + 2
         END IF
*
         IF( N1.EQ.2 ) THEN
*
*           Standardize new 2-by-2 block T22
*
            J3 = J1 + N2
            J4 = J3 + 1
            CALL SLANV2( T( J3, J3 ), T( J3, J4 ), T( J4, J3 ),
     $                   T( J4, J4 ), WR1, WI1, WR2, WI2, CS, SN )
            IF( J3+2.LE.N )
     $         CALL SROT( N-J3-1, T( J3, J3+2 ), LDT, T( J4, J3+2 ),
     $                    LDT, CS, SN )
            CALL SROT( J3-1, T( 1, J3 ), 1, T( 1, J4 ), 1, CS, SN )
            ITRAF( LI ) = J3
            DTRAF( LD ) = CS
            DTRAF( LD+1 ) = SN
         END IF
*
      END IF
      RETURN
*
*     Exit with INFO = 1 if swap was rejected.
*
   50 CONTINUE
      INFO = 1
      RETURN
*
*     End of BSLAEXC
*
      END
