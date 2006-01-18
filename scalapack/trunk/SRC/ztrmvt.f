      SUBROUTINE ZTRMVT( UPLO, N, T, LDT, X, INCX, Y, INCY, W, INCW, Z,
     $                   INCZ )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 13, 2000
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INCW, INCX, INCY, INCZ, LDT, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         T( LDT, * ), W( * ), X( * ), Y( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  ZTRMVT  performs the matrix-vector operations
*
*     x := conjg( T' ) *y, and w := T *z,
*
*  where x is an n element vector and  T is an n by n
*  upper or lower triangular matrix.
*
*  Arguments 
*  =========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  T      - COMPLEX*16 array of DIMENSION ( LDT, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array T must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           T is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array T must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           T is not referenced.
*
*  LDT    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - COMPLEX*16 array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           On exit, X = T' * y
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16 array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.  Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  W      - COMPLEX*16 array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCW ) ).
*           On exit, W = T * z
*
*  INCW   - INTEGER.
*           On entry, INCW specifies the increment for the elements of
*           W. INCW must not be zero.
*           Unchanged on exit.
*
*  Z      - COMPLEX*16 array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCZ ) ).
*           Before entry, the incremented array Z must contain the n
*           element vector z.  Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entrz, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*
*     .. Local Scalars ..
      INTEGER            INFO
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = 1
      ELSE IF( N.LT.0 ) THEN
         INFO = 2
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = 4
      ELSE IF( INCW.EQ.0 ) THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 ) THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 ) THEN
         INFO = 10
      ELSE IF( INCZ.EQ.0 ) THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRMVT', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*
*
      IF( INCX.NE.1 .OR. INCY.NE.1 .OR. INCW.NE.1 .OR. INCZ.NE.1 .OR.
     $    .TRUE. ) THEN
         CALL ZCOPY( N, Y, INCY, X, INCX )
         CALL ZTRMV( UPLO, 'C', 'N', N, T, LDT, X, INCX )
         CALL ZCOPY( N, Z, INCZ, W, INCW )
         CALL ZTRMV( UPLO, 'N', 'N', N, T, LDT, W, INCW )
         RETURN
      END IF
*
      RETURN
*
*     End of ZTRMVT.
*
      END
