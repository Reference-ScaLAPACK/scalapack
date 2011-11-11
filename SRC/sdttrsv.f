      SUBROUTINE SDTTRSV( UPLO, TRANS, N, NRHS, DL, D, DU,
     $                   B, LDB, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*
*     Written by Andrew J. Cleary, University of Tennessee.
*     August, 1996.
*     Modified from SGTTRS:
*  -- LAPACK routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, TRANS
      INTEGER            INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      REAL               B( LDB, * ), D( * ), DL( * ), DU( * )
*     ..
*
*  Purpose
*  =======
*
*  SDTTRSV solves one of the systems of equations
*     L * X = B,  L**T * X = B,  or  L**H * X = B,
*     U * X = B,  U**T * X = B,  or  U**H * X = B,
*  with factors of the tridiagonal matrix A from the LU factorization
*  computed by SDTTRF.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether to solve with L or U.
*
*  TRANS   (input) CHARACTER
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B     (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  DL      (input) COMPLEX array, dimension (N-1)
*          The (n-1) multipliers that define the matrix L from the
*          LU factorization of A.
*
*  D       (input) COMPLEX array, dimension (N)
*          The n diagonal elements of the upper triangular matrix U from
*          the LU factorization of A.
*
*  DU      (input) COMPLEX array, dimension (N-1)
*          The (n-1) elements of the first superdiagonal of U.
*
*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, B is overwritten by the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LOWER, NOTRAN
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SDTTRSV', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
      IF( LOWER ) THEN
*        Solve L*X = B, overwriting B with X.
*
         DO 35 J = 1, NRHS
*
*           Solve L*x = b.
*
            DO 10 I = 1, N - 1
                  B( I+1, J ) = B( I+1, J ) - DL( I )*B( I, J )
   10       CONTINUE
   35    CONTINUE
*
      ELSE
*        Solve U*x = b.
*
         DO 30 J = 1, NRHS
            B( N, J ) = B( N, J ) / D( N )
            IF( N.GT.1 )
     $         B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) /
     $                       D( N-1 )
            DO 20 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J ) ) / D( I )
   20       CONTINUE
   30    CONTINUE
*
      ENDIF
*
      ELSE
*
         IF( .NOT. LOWER ) THEN
*        Solve U**T * X = B, overwriting B with X.
*
         DO 65 J = 1, NRHS
*
*           Solve U**T * x = b.
*
            B( 1, J ) = B( 1, J ) / D( 1 )
            IF( N.GT.1 )
     $         B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 )
            DO 40 I = 3, N
               B( I, J ) = ( B( I, J )-DU( I-1 )*B( I-1, J ) ) / D( I )
   40       CONTINUE
   65    CONTINUE
*
         ELSE
*
*        Solve L**T * X = B, overwriting B with X.
         DO 60 J = 1, NRHS
*
*           Solve L**T * x = b.
*
            DO 50 I = N - 1, 1, -1
                  B( I, J ) = B( I, J ) - DL( I )*B( I+1, J )
   50       CONTINUE
   60    CONTINUE
         ENDIF
      END IF
*
*     End of SDTTRSV
*
      END
