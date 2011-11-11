      SUBROUTINE ZDTTRF( N, DL, D, DU, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*
*     Written by Andrew J. Cleary, November 1996.
*     Modified from ZGTTRF:
*  -- LAPACK routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         D( * ), DL( * ), DU( * )
*     ..
*
*  Purpose
*  =======
*
*  ZDTTRF computes an LU factorization of a complex tridiagonal matrix A
*  using elimination without partial pivoting.
*
*  The factorization has the form
*     A = L * U
*  where L is a product of unit lower bidiagonal
*  matrices and U is upper triangular with nonzeros in only the main
*  diagonal and first superdiagonal.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  DL      (input/output) COMPLEX array, dimension (N-1)
*          On entry, DL must contain the (n-1) subdiagonal elements of
*          A.
*          On exit, DL is overwritten by the (n-1) multipliers that
*          define the matrix L from the LU factorization of A.
*
*  D       (input/output) COMPLEX array, dimension (N)
*          On entry, D must contain the diagonal elements of A.
*          On exit, D is overwritten by the n diagonal elements of the
*          upper triangular matrix U from the LU factorization of A.
*
*  DU      (input/output) COMPLEX array, dimension (N-1)
*          On entry, DU must contain the (n-1) superdiagonal elements
*          of A.
*          On exit, DU is overwritten by the (n-1) elements of the first
*          superdiagonal of U.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16         FACT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Parameters ..
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'ZDTTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      DO 20 I = 1, N - 1
         IF( DL( I ).EQ.CZERO ) THEN
*
*           Subdiagonal is zero, no elimination is required.
*
            IF( D( I ).EQ.CZERO .AND. INFO.EQ.0 )
     $         INFO = I
         ELSE
*
            FACT = DL( I ) / D( I )
            DL( I ) = FACT
            D( I+1 ) = D( I+1 ) - FACT*DU( I )
         END IF
   20 CONTINUE
      IF( D( N ).EQ.CZERO .AND. INFO.EQ.0 ) THEN
         INFO = N
         RETURN
      END IF
*
      RETURN
*
*     End of ZDTTRF
*
      END
