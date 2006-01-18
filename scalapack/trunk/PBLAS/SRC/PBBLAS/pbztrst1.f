      SUBROUTINE PBZTRST1( ICONTXT, XDIST, N, NB, NZ, X, INCX, BETA, Y,
     $                     INCY, LCMP, LCMQ, NINT )
*
*  -- PB-BLAS routine (version 2.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory.
*     April 28, 1996
*
*     .. Scalar Arguments ..
      CHARACTER*1        XDIST
      INTEGER            ICONTXT, INCX, INCY, LCMP, LCMQ, N, NB, NINT,
     $                   NZ
      COMPLEX*16         BETA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  PBZTRST1 forms  y <== x + beta * y, where y is a sorted
*  condensed row (or column) vector from a column (or row) vector of x.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE  = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Variables ..
      INTEGER            ITER, IX, IY, K, KK, KZ, NJUMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           PBZVECADD
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL
      EXTERNAL           ICEIL, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX, MOD
*     ..
*     .. Executable Statements ..
*
      ITER = ICEIL( NINT,  NB )
      KZ   = NZ
*
      IF( LSAME( XDIST, 'R' ) ) THEN
         NJUMP = NB * LCMQ
*
         DO 20 KK = 0, LCMQ-1
            IX = NINT * MOD( KK*LCMP, LCMQ )
            IY = MAX( 0, NB*KK-NZ )
            IF( N.LT.IY ) GO TO 50
*
            IF( ITER.GT.1 ) THEN
               CALL PBZVECADD( ICONTXT, 'G', NB-KZ, ONE, X(IX*INCX+1),
     $                         INCX, BETA, Y(IY*INCY+1), INCY )
               IX = IX + NB - KZ
               IY = IY + NJUMP - KZ
               KZ = 0
*
               DO 10 K = 2, ITER-1
                  CALL PBZVECADD( ICONTXT, 'G', NB, ONE, X(IX*INCX+1),
     $                            INCX, BETA, Y(IY*INCY+1), INCY )
                  IX = IX + NB
                  IY = IY + NJUMP
   10          CONTINUE
            END IF
*
            CALL PBZVECADD( ICONTXT, 'G', MIN(NB-KZ,N-IY), ONE,
     $                      X(IX*INCX+1), INCX, BETA, Y(IY*INCY+1),
     $                      INCY )
            KZ = 0
   20    CONTINUE
*
*     if( LSAME( XDIST, 'C' ) ) then
*
      ELSE
         NJUMP = NB * LCMP
*
         DO 40 KK = 0, LCMP-1
            IX = NINT * MOD( KK*LCMQ, LCMP )
            IY = MAX( 0, NB*KK-NZ )
            IF( N.LT.IY ) GO TO 50
*
            IF( ITER.GT.1 ) THEN
               CALL PBZVECADD( ICONTXT, 'G', NB-KZ, ONE, X(IX*INCX+1),
     $                         INCX, BETA, Y(IY*INCY+1), INCY )
               IX = IX + NB - KZ
               IY = IY + NJUMP - KZ
               KZ = 0
*
               DO 30 K = 2, ITER-1
                  CALL PBZVECADD( ICONTXT, 'G', NB, ONE, X(IX*INCX+1),
     $                            INCX, BETA, Y(IY*INCY+1), INCY )
                  IX = IX + NB
                  IY = IY + NJUMP
   30          CONTINUE
            END IF
*
            CALL PBZVECADD( ICONTXT, 'G', MIN(NB-KZ,N-IY), ONE,
     $                      X(IX*INCX+1), INCX, BETA, Y(IY*INCY+1),
     $                      INCY )
            KZ = 0
   40    CONTINUE
      END IF
*
   50 CONTINUE
*
      RETURN
*
*     End of PBZTRST1
*
      END
