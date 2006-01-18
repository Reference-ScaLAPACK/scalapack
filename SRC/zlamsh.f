      SUBROUTINE ZLAMSH( S, LDS, NBULGE, JBLK, H, LDH, N, ULP )
*
*  -- ScaLAPACK routine (version 1.7) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     May 28, 1999
*
*     .. Scalar Arguments ..
      INTEGER            JBLK, LDH, LDS, N, NBULGE
      DOUBLE PRECISION   ULP
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), S( LDS, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLAMSH sends multiple shifts through a small (single node) matrix to
*     see how consecutive small subdiagonal elements are modified by
*     subsequent shifts in an effort to maximize the number of bulges
*     that can be sent through.
*  ZLAMSH should only be called when there are multiple shifts/bulges
*     (NBULGE > 1) and the first shift is starting in the middle of an
*     unreduced Hessenberg matrix because of two or more consecutive
*     small subdiagonal elements.
*
*  Arguments
*  =========
*
*  S       (local input/output) COMPLEX*16 array, ( LDS,* )
*          On entry, the matrix of shifts.  Only the 2x2 diagonal of S
*          is referenced.  It is assumed that S has JBLK double shifts
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
*  H       (local input/output) COMPLEX*16 array ( LDH,N )
*          On entry, the local matrix to apply the shifts on.
*              H should be aligned so that the starting row is 2.
*          On exit, the data is destroyed.
*
*  LDH     (local input) INTEGER
*          On entry, the leading dimension of H.  Unchanged on exit.
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
*  Further Details
*  ===============
*
*  Implemented by:  M. Fahey, May 28, 1999
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   RONE, TEN
      PARAMETER          ( RONE = 1.0D+0, TEN = 10.0D+0 )
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IBULGE, IVAL, J, K, M, NR
      DOUBLE PRECISION   DVAL, S1, TST1
      COMPLEX*16         CDUM, H00, H10, H11, H12, H21, H22, H33, H33S,
     $                   H43H34, H44, H44S, SUM, T1, T2, T3, V1, V2, V3
*     ..
*     .. Local Arrays ..
      COMPLEX*16         V( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZLARFG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
      M = 2
      DO 50 IBULGE = 1, NBULGE
         H44 = S( 2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+2 )
         H33 = S( 2*JBLK-2*IBULGE+1, 2*JBLK-2*IBULGE+1 )
         H43H34 = S( 2*JBLK-2*IBULGE+1, 2*JBLK-2*IBULGE+2 )*
     $            S( 2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+1 )
         H11 = H( M, M )
         H22 = H( M+1, M+1 )
         H21 = H( M+1, M )
         H12 = H( M, M+1 )
         H44S = H44 - H11
         H33S = H33 - H11
         V1 = ( H33S*H44S-H43H34 ) / H21 + H12
         V2 = H22 - H11 - H33S - H44S
         V3 = H( M+2, M+1 )
         S1 = CABS1( V1 ) + CABS1( V2 ) + CABS1( V3 )
         V1 = V1 / S1
         V2 = V2 / S1
         V3 = V3 / S1
         V( 1 ) = V1
         V( 2 ) = V2
         V( 3 ) = V3
         H00 = H( M-1, M-1 )
         H10 = H( M, M-1 )
         TST1 = CABS1( V1 )*( CABS1( H00 )+CABS1( H11 )+CABS1( H22 ) )
         IF( CABS1( H10 )*( CABS1( V2 )+CABS1( V3 ) ).GT.ULP*TST1 ) THEN
*           Find minimum
            DVAL = ( CABS1( H10 )*( CABS1( V2 )+CABS1( V3 ) ) ) /
     $             ( ULP*TST1 )
            IVAL = IBULGE
            DO 10 I = IBULGE + 1, NBULGE
               H44 = S( 2*JBLK-2*I+2, 2*JBLK-2*I+2 )
               H33 = S( 2*JBLK-2*I+1, 2*JBLK-2*I+1 )
               H43H34 = S( 2*JBLK-2*I+1, 2*JBLK-2*I+2 )*
     $                  S( 2*JBLK-2*I+2, 2*JBLK-2*I+1 )
               H11 = H( M, M )
               H22 = H( M+1, M+1 )
               H21 = H( M+1, M )
               H12 = H( M, M+1 )
               H44S = H44 - H11
               H33S = H33 - H11
               V1 = ( H33S*H44S-H43H34 ) / H21 + H12
               V2 = H22 - H11 - H33S - H44S
               V3 = H( M+2, M+1 )
               S1 = CABS1( V1 ) + CABS1( V2 ) + CABS1( V3 )
               V1 = V1 / S1
               V2 = V2 / S1
               V3 = V3 / S1
               V( 1 ) = V1
               V( 2 ) = V2
               V( 3 ) = V3
               H00 = H( M-1, M-1 )
               H10 = H( M, M-1 )
               TST1 = CABS1( V1 )*( CABS1( H00 )+CABS1( H11 )+
     $                CABS1( H22 ) )
               IF( ( DVAL.GT.( CABS1( H10 )*( CABS1( V2 )+
     $             CABS1( V3 ) ) ) / ( ULP*TST1 ) ) .AND.
     $             ( DVAL.GT.RONE ) ) THEN
                  DVAL = ( CABS1( H10 )*( CABS1( V2 )+CABS1( V3 ) ) ) /
     $                   ( ULP*TST1 )
                  IVAL = I
               END IF
   10       CONTINUE
            IF( ( DVAL.LT.TEN ) .AND. ( IVAL.NE.IBULGE ) ) THEN
               H44 = S( 2*JBLK-2*IVAL+2, 2*JBLK-2*IVAL+2 )
               H33 = S( 2*JBLK-2*IVAL+1, 2*JBLK-2*IVAL+1 )
               H43H34 = S( 2*JBLK-2*IVAL+1, 2*JBLK-2*IVAL+2 )
               H10 = S( 2*JBLK-2*IVAL+2, 2*JBLK-2*IVAL+1 )
               S( 2*JBLK-2*IVAL+2, 2*JBLK-2*IVAL+2 ) = S( 2*JBLK-2*
     $            IBULGE+2, 2*JBLK-2*IBULGE+2 )
               S( 2*JBLK-2*IVAL+1, 2*JBLK-2*IVAL+1 ) = S( 2*JBLK-2*
     $            IBULGE+1, 2*JBLK-2*IBULGE+1 )
               S( 2*JBLK-2*IVAL+1, 2*JBLK-2*IVAL+2 ) = S( 2*JBLK-2*
     $            IBULGE+1, 2*JBLK-2*IBULGE+2 )
               S( 2*JBLK-2*IVAL+2, 2*JBLK-2*IVAL+1 ) = S( 2*JBLK-2*
     $            IBULGE+2, 2*JBLK-2*IBULGE+1 )
               S( 2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+2 ) = H44
               S( 2*JBLK-2*IBULGE+1, 2*JBLK-2*IBULGE+1 ) = H33
               S( 2*JBLK-2*IBULGE+1, 2*JBLK-2*IBULGE+2 ) = H43H34
               S( 2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+1 ) = H10
            END IF
            H44 = S( 2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+2 )
            H33 = S( 2*JBLK-2*IBULGE+1, 2*JBLK-2*IBULGE+1 )
            H43H34 = S( 2*JBLK-2*IBULGE+1, 2*JBLK-2*IBULGE+2 )*
     $               S( 2*JBLK-2*IBULGE+2, 2*JBLK-2*IBULGE+1 )
            H11 = H( M, M )
            H22 = H( M+1, M+1 )
            H21 = H( M+1, M )
            H12 = H( M, M+1 )
            H44S = H44 - H11
            H33S = H33 - H11
            V1 = ( H33S*H44S-H43H34 ) / H21 + H12
            V2 = H22 - H11 - H33S - H44S
            V3 = H( M+2, M+1 )
            S1 = CABS1( V1 ) + CABS1( V2 ) + CABS1( V3 )
            V1 = V1 / S1
            V2 = V2 / S1
            V3 = V3 / S1
            V( 1 ) = V1
            V( 2 ) = V2
            V( 3 ) = V3
            H00 = H( M-1, M-1 )
            H10 = H( M, M-1 )
            TST1 = CABS1( V1 )*( CABS1( H00 )+CABS1( H11 )+
     $             CABS1( H22 ) )
         END IF
         IF( CABS1( H10 )*( CABS1( V2 )+CABS1( V3 ) ).GT.TEN*ULP*TST1 )
     $        THEN
*           IBULGE better not be 1 here or we have a bug!
            NBULGE = MAX( IBULGE-1, 1 )
            RETURN
         END IF
         DO 40 K = M, N - 1
            NR = MIN( 3, N-K+1 )
            IF( K.GT.M )
     $         CALL ZCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL ZLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
               IF( K.LT.N-1 )
     $            H( K+2, K-1 ) = ZERO
            ELSE
*              H(m,m-1) must be updated,
*
               H( K, K-1 ) = H( K, K-1 ) - DCONJG( T1 )*H( K, K-1 )
            END IF
            V2 = V( 2 )
            T2 = T1*V2
            IF( NR.EQ.3 ) THEN
               V3 = V( 3 )
               T3 = T1*V3
               DO 20 J = K, N
                  SUM = DCONJG( T1 )*H( K, J ) +
     $                  DCONJG( T2 )*H( K+1, J ) +
     $                  DCONJG( T3 )*H( K+2, J )
                  H( K, J ) = H( K, J ) - SUM
                  H( K+1, J ) = H( K+1, J ) - SUM*V2
                  H( K+2, J ) = H( K+2, J ) - SUM*V3
   20          CONTINUE
               DO 30 J = 1, MIN( K+3, N )
                  SUM = T1*H( J, K ) + T2*H( J, K+1 ) + T3*H( J, K+2 )
                  H( J, K ) = H( J, K ) - SUM
                  H( J, K+1 ) = H( J, K+1 ) - SUM*DCONJG( V2 )
                  H( J, K+2 ) = H( J, K+2 ) - SUM*DCONJG( V3 )
   30          CONTINUE
            END IF
   40    CONTINUE
   50 CONTINUE
*
      RETURN
*
*     End of ZLAMSH
*
      END
