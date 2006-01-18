      SUBROUTINE ILACPY( UPLO, M, N, A, LDA, B, LDB )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  ILACPY copies all or part of a local array A to another array B.
*
*  Arguments
*  =========
*
*  UPLO    (local input) CHARACTER*1
*          Specifies the part of the array A to be copied to B.
*          = 'U':      Upper triangular part
*          = 'L':      Lower triangular part
*          Otherwise:  All of the array A
*
*  M       (local input) INTEGER
*          The number of rows of the array A.  M >= 0.
*
*  N       (local input) INTEGER
*          The number of columns of the array A.  N >= 0.
*
*  A       (local input) INTEGER
*          Array, dimension (LDA,N), the m by n array A.
*          If UPLO = 'U', only the upper triangle or trapezoid is
*          accessed; if UPLO = 'L', only the lower triangle or trapezoid
*          is accessed.
*
*  LDA     (local input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (local output) INTEGER
*          Array, dimension (LDB,N), on exit, B = A in the locations
*          specified by UPLO.
*
*  LDB     (local input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
*
*     End of ILACPY
*
      END
