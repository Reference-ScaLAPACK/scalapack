      SUBROUTINE ZPTTRSV( UPLO, TRANS, N, NRHS, D, E, B, LDB,
     $                        INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*
*     Written by Andrew J. Cleary, University of Tennessee.
*     November, 1996.
*     Modified from ZPTTRS:
*  -- LAPACK routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, TRANS
      INTEGER            INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
      COMPLEX*16         B( LDB, * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  ZPTTRSV solves one of the triangular systems
*     L * X = B, or  L**H * X = B,
*     U * X = B, or  U**H * X = B,
*  where L or U is the Cholesky factor of a Hermitian positive
*  definite tridiagonal matrix A such that
*  A = U**H*D*U or A = L*D*L**H (computed by ZPTTRF).
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the superdiagonal or the subdiagonal
*          of the tridiagonal matrix A is stored and the form of the
*          factorization:
*          = 'U':  E is the superdiagonal of U, and A = U'*D*U;
*          = 'L':  E is the subdiagonal of L, and A = L*D*L'.
*          (The two forms are equivalent if A is real.)
*
*  TRANS   (input) CHARACTER
*          Specifies the form of the system of equations:
*          = 'N':  L * X = B     (No transpose)
*          = 'N':  L * X = B     (No transpose)
*          = 'C':  U**H * X = B  (Conjugate transpose)
*          = 'C':  L**H * X = B  (Conjugate transpose)
*
*  N       (input) INTEGER
*          The order of the tridiagonal matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  D       (input) REAL array, dimension (N)
*          The n diagonal elements of the diagonal matrix D from the
*          factorization computed by ZPTTRF.
*
*  E       (input) COMPLEX array, dimension (N-1)
*          The (n-1) off-diagonal elements of the unit bidiagonal
*          factor U or L from the factorization computed by ZPTTRF
*          (see UPLO).
*
*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
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
      LOGICAL            NOTRAN, UPPER
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
      INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPTTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
      IF( .NOT.NOTRAN ) THEN
*
         DO 30 J = 1, NRHS
*
*           Solve U**T (or H) * x = b.
*
            DO 10 I = 2, N
               B( I, J ) = B( I, J ) - B( I-1, J )*DCONJG( E( I-1 ) )
   10       CONTINUE
   30    CONTINUE
*
      ELSE
*
         DO 35 J = 1, NRHS
*
*           Solve U * x = b.
*
            DO 20 I = N - 1, 1, -1
               B( I, J ) = B( I, J ) - B( I+1, J )*E( I )
   20       CONTINUE
   35    CONTINUE
      ENDIF
*
      ELSE
*
      IF( NOTRAN ) THEN
*
         DO 60 J = 1, NRHS
*
*           Solve L * x = b.
*
            DO 40 I = 2, N
               B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
   40       CONTINUE
   60    CONTINUE
*
      ELSE
*
         DO 65 J = 1, NRHS
*
*           Solve L**H * x = b.
*
            DO 50 I = N - 1, 1, -1
               B( I, J ) = B( I, J ) -
     $                     B( I+1, J )*DCONJG( E( I ) )
   50       CONTINUE
   65    CONTINUE
      ENDIF
*
      END IF
*
      RETURN
*
*     End of ZPTTRS
*
      END
