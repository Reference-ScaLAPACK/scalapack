      SUBROUTINE SPTTRSV( TRANS, N, NRHS, D, E, B, LDB,
     $                        INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*
*     Written by Andrew J. Cleary, University of Tennessee.
*     November, 1996.
*     Modified from SPTTRS:
*  -- LAPACK routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      REAL               D( * )
      REAL               B( LDB, * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  SPTTRSV solves one of the triangular systems
*     L**T* X = B, or  L * X = B,
*  where L is the Cholesky factor of a Hermitian positive
*  definite tridiagonal matrix A such that
*  A = L*D*L**H (computed by SPTTRF).
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER
*          Specifies the form of the system of equations:
*          = 'N':  L * X = B     (No transpose)
*          = 'T':  L**T * X = B  (Transpose)
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
*          factorization computed by SPTTRF.
*
*  E       (input) COMPLEX array, dimension (N-1)
*          The (n-1) off-diagonal elements of the unit bidiagonal
*          factor U or L from the factorization computed by SPTTRF
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
      LOGICAL            NOTRAN
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
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPTTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
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
     $                     B( I+1, J )*( E( I ) )
   50       CONTINUE
   65    CONTINUE
      ENDIF
*
      RETURN
*
*     End of SPTTRS
*
      END
