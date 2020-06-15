      SUBROUTINE ZZDOTC( N, DOTC, X, INCX, Y, INCY )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      COMPLEX*16         DOTC
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZZDOTC is a simple FORTRAN wrapper around the BLAS function
*  ZDOTC returning the result in the parameter list instead.
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX*16         ZTEMP
      INTEGER            I,IX,IY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG
*     ..
*     .. Executable Statements ..
*
      ZTEMP = (0.0d0,0.0d0)
      DOTC = (0.0d0,0.0d0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO i = 1,N
            ZTEMP = ZTEMP + DCONJG(X(I))*Y(I)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            ZTEMP = ZTEMP + DCONJG(X(IX))*Y(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      DOTC = ZTEMP
*
      RETURN
*
*     End of ZZDOTC
*
      END
