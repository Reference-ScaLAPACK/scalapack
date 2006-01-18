      SUBROUTINE DLAMSH ( S, LDS, NBULGE, JBLK, H, LDH, N, ULP )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            LDS, NBULGE, JBLK, LDH, N
      DOUBLE PRECISION   ULP
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   S(LDS,*), H(LDH,*)
*     ..
*
*  Purpose
*  =======
*
*  DLAMSH sends multiple shifts through a small (single node) matrix to
*     see how consecutive small subdiagonal elements are modified by 
*     subsequent shifts in an effort to maximize the number of bulges 
*     that can be sent through.
*  DLAMSH should only be called when there are multiple shifts/bulges 
*     (NBULGE > 1) and the first shift is starting in the middle of an
*     unreduced Hessenberg matrix because of two or more consecutive small
*     subdiagonal elements.
*
*  Arguments
*  =========
*
*  S       (local input/output) DOUBLE PRECISION array, (LDS,*)
*          On entry, the matrix of shifts.  Only the 2x2 diagonal of S is
*             referenced.  It is assumed that S has JBLK double shifts
*             (size 2).
*          On exit, the data is rearranged in the best order for 
*             applying.
*
*  LDS     (local input) INTEGER
*          On entry, the leading dimension of S.  Unchanged on exit.
*              1 < NBULGE <= JBLK <= LDS/2
*
*  NBULGE  (local input/output) INTEGER
*          On entry, the number of bulges to send through H ( >1 ).
*              NBULGE should be less than the maximum determined (JBLK).
*              1 < NBULGE <= JBLK <= LDS/2
*          On exit, the maximum number of bulges that can be sent 
*              through.
*
*  JBLK    (local input) INTEGER
*          On entry, the number of shifts determined for S.
*          Unchanged on exit.
*
*  H       (local input/output) DOUBLE PRECISION array (LDH,N)
*          On entry, the local matrix to apply the shifts on.
*              H should be aligned so that the starting row is 2.
*          On exit, the data is destroyed.
*
*  LDS     (local input) INTEGER
*          On entry, the leading dimension of S.  Unchanged on exit.
*
*  N       (local input) INTEGER
*          On entry, the size of H.  If all the bulges are expected to
*              go through, N should be at least 4*NBULGE+2.
*              Otherwise, NBULGE may be reduced by this routine.
*
*  ULP     (local input) DOUBLE PRECISION
*          On entry, machine precision
*          Unchanged on exit.
*
*  Implemented by:  G. Henry, May 1, 1997
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO, TEN
      PARAMETER ( ZERO = 0.0D+0, TEN = 10.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER          K, IBULGE, M, NR, J, IVAL, I
      DOUBLE PRECISION H44, H33, H43H34, H11, H22, H21, H12, H44S,
     $                 H33S, V1, V2, V3, H00, H10, TST1, T1, T2, T3,
     $                 SUM, S1, DVAL
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION V(3)
*     ..
*     .. External Subroutines ..
      EXTERNAL         DLARFG, DCOPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        MAX, ABS
*     ..
*     .. Executable Statements ..
*
      M = 2
      DO 10 IBULGE = 1, NBULGE
         H44 = S(2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+2)
         H33 = S(2*JBLK-2*IBULGE+1,2*JBLK-2*IBULGE+1)
         H43H34 = S(2*JBLK-2*IBULGE+1,2*JBLK-2*IBULGE+2)*
     $            S(2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+1)
         H11 = H( M, M )
         H22 = H( M+1, M+1 )
         H21 = H( M+1, M )
         H12 = H( M, M+1 )
         H44S = H44 - H11
         H33S = H33 - H11
         V1 = ( H33S*H44S-H43H34 ) / H21 + H12
         V2 = H22 - H11 - H33S - H44S
         V3 = H( M+2, M+1 )
         S1 = ABS( V1 ) + ABS( V2 ) + ABS( V3 )
         V1 = V1 / S1
         V2 = V2 / S1
         V3 = V3 / S1
         V( 1 ) = V1
         V( 2 ) = V2
         V( 3 ) = V3
         H00 = H( M-1, M-1 )
         H10 = H( M, M-1 )
         TST1 = ABS( V1 )*( ABS( H00 )+ABS( H11 )+ABS( H22 ) )
         IF( ABS( H10 )*( ABS( V2 )+ABS( V3 ) ).GT.ULP*TST1 ) THEN
*           Find minimum
            DVAL = (ABS(H10)*(ABS(V2)+ABS(V3))) / (ULP*TST1)   
            IVAL = IBULGE
            DO 15 I = IBULGE+1, NBULGE
               H44 = S(2*JBLK-2*I+2, 2*JBLK-2*I+2)
               H33 = S(2*JBLK-2*I+1,2*JBLK-2*I+1)
               H43H34 = S(2*JBLK-2*I+1,2*JBLK-2*I+2)*
     $                  S(2*JBLK-2*I+2, 2*JBLK-2*I+1)
               H11 = H( M, M )
               H22 = H( M+1, M+1 )
               H21 = H( M+1, M )
               H12 = H( M, M+1 )
               H44S = H44 - H11
               H33S = H33 - H11
               V1 = ( H33S*H44S-H43H34 ) / H21 + H12
               V2 = H22 - H11 - H33S - H44S
               V3 = H( M+2, M+1 )
               S1 = ABS( V1 ) + ABS( V2 ) + ABS( V3 )
               V1 = V1 / S1
               V2 = V2 / S1
               V3 = V3 / S1
               V( 1 ) = V1
               V( 2 ) = V2
               V( 3 ) = V3
               H00 = H( M-1, M-1 )
               H10 = H( M, M-1 )
               TST1 = ABS( V1 )*( ABS( H00 )+ABS( H11 )+ABS( H22 ) )
               IF ( (DVAL.GT.(ABS(H10)*(ABS(V2)+ABS(V3)))/(ULP*TST1))
     $             .AND. ( DVAL .GT. 1.D0 ) ) THEN
                  DVAL = (ABS(H10)*(ABS(V2)+ABS(V3))) / (ULP*TST1)   
                  IVAL = I
               END IF
  15        CONTINUE
            IF ( (DVAL .LT. TEN) .AND. (IVAL .NE. IBULGE) ) THEN
               H44 = S(2*JBLK-2*IVAL+2, 2*JBLK-2*IVAL+2)
               H33 = S(2*JBLK-2*IVAL+1,2*JBLK-2*IVAL+1)
               H43H34 = S(2*JBLK-2*IVAL+1,2*JBLK-2*IVAL+2)
               H10 =    S(2*JBLK-2*IVAL+2, 2*JBLK-2*IVAL+1)
               S(2*JBLK-2*IVAL+2,2*JBLK-2*IVAL+2) = 
     $              S(2*JBLK-2*IBULGE+2,2*JBLK-2*IBULGE+2)
               S(2*JBLK-2*IVAL+1,2*JBLK-2*IVAL+1) =
     $              S(2*JBLK-2*IBULGE+1,2*JBLK-2*IBULGE+1)
               S(2*JBLK-2*IVAL+1,2*JBLK-2*IVAL+2) =
     $              S(2*JBLK-2*IBULGE+1,2*JBLK-2*IBULGE+2) 
               S(2*JBLK-2*IVAL+2, 2*JBLK-2*IVAL+1) =
     $              S(2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+1) 
               S(2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+2) = H44
               S(2*JBLK-2*IBULGE+1,2*JBLK-2*IBULGE+1) = H33
               S(2*JBLK-2*IBULGE+1,2*JBLK-2*IBULGE+2) = H43H34
               S(2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+1) = H10
            END IF
            H44 = S(2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+2)
            H33 = S(2*JBLK-2*IBULGE+1,2*JBLK-2*IBULGE+1)
            H43H34 = S(2*JBLK-2*IBULGE+1,2*JBLK-2*IBULGE+2)*
     $               S(2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+1)
            H11 = H( M, M )
            H22 = H( M+1, M+1 )
            H21 = H( M+1, M )
            H12 = H( M, M+1 )
            H44S = H44 - H11
            H33S = H33 - H11
            V1 = ( H33S*H44S-H43H34 ) / H21 + H12
            V2 = H22 - H11 - H33S - H44S
            V3 = H( M+2, M+1 )
            S1 = ABS( V1 ) + ABS( V2 ) + ABS( V3 )
            V1 = V1 / S1
            V2 = V2 / S1
            V3 = V3 / S1
            V( 1 ) = V1
            V( 2 ) = V2
            V( 3 ) = V3
            H00 = H( M-1, M-1 )
            H10 = H( M, M-1 )
            TST1 = ABS( V1 )*( ABS( H00 )+ABS( H11 )+ABS( H22 ) )
         END IF
         IF( ABS( H10 )*( ABS( V2 )+ABS( V3 ) ).GT.TEN*ULP*TST1 ) THEN
*           IBULGE better not be 1 here or we have a bug!
            NBULGE = MAX(IBULGE -1,1)
            RETURN
         END IF
         DO 120 K = M, N - 1
            NR = MIN( 3, N-K+1 )
            IF( K.GT.M )
     $         CALL DCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL DLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
               IF( K.LT.N-1 )
     $            H( K+2, K-1 ) = ZERO
            ELSE
               H( K, K-1 ) = -H( K, K-1 )
            END IF
            V2 = V( 2 )
            T2 = T1*V2
            IF( NR.EQ.3 ) THEN
               V3 = V( 3 )
               T3 = T1*V3
               DO 60 J = K, N
                  SUM = H( K, J ) + V2*H( K+1, J ) + V3*H( K+2, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
                  H( K+2, J ) = H( K+2, J ) - SUM*T3
   60          CONTINUE
               DO 70 J = 1, MIN( K+3, N )
                  SUM = H( J, K ) + V2*H( J, K+1 ) + V3*H( J, K+2 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
                  H( J, K+2 ) = H( J, K+2 ) - SUM*T3
   70          CONTINUE
            END IF
  120    CONTINUE
   10 CONTINUE
* 
      RETURN
      END
