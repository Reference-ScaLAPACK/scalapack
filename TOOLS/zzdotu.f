      SUBROUTINE ZZDOTU( N, DOTU, X, INCX, Y, INCY )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      COMPLEX*16         DOTU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZZDOTU is a simple FORTRAN wrapper around the BLAS function
*  ZDOTU returning the result in the parameter list instead.
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX*16         ZTEMP
      INTEGER            I,IX,IY
*     ..
*     .. Executable Statements ..
*
      ZTEMP = (0.0d0,0.0d0)
      DOTU = (0.0d0,0.0d0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO i = 1,N
            ZTEMP = ZTEMP + X(I)*Y(I)
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
            ZTEMP = ZTEMP + X(IX)*Y(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      DOTU = ZTEMP
*
      RETURN
*
*     End of ZZDOTU
*
      END
