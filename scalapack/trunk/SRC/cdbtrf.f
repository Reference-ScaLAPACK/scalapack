      SUBROUTINE CDBTRF( M, N, KL, KU, AB, LDAB, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*
*     Written by Andrew J. Cleary, University of Tennessee.
*     August, 1996.
*     Modified from CGBTRF:
*  -- LAPACK routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     August 6, 1991
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  Cdbtrf computes an LU factorization of a real m-by-n band matrix A
*  without using partial pivoting or row interchanges.
*
*  This is the blocked version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  AB      (input/output) REAL array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows KL+1 to
*          2*KL+KU+1; rows 1 to KL of the array need not be set.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
*
*          On exit, details of the factorization: U is stored as an
*          upper triangular band matrix with KL+KU superdiagonals in
*          rows 1 to KL+KU+1, and the multipliers used during the
*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*          See below for further details.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  M = N = 6, KL = 2, KU = 1:
*
*  On entry:                       On exit:
*
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
*
*  Array elements marked * are not used by the routine.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0 )
      PARAMETER          ( ZERO = 0.0E+0 )
      COMPLEX            CONE, CZERO
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I2, I3, II, J, J2, J3, JB, JJ, JM, JP,
     $                   JU, KM, KV, NB, NW
*     ..
*     .. Local Arrays ..
      COMPLEX              WORK13( LDWORK, NBMAX ),
     $                   WORK31( LDWORK, NBMAX )
*     ..
*     .. External Functions ..
      INTEGER            ILAENV, ISAMAX
      EXTERNAL           ILAENV, ISAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CDBTF2, CGEMM, CGERU, CSCAL,
     $                   CSWAP, CTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     KV is the number of superdiagonals in the factor U
*
      KV = KU
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.MIN( MIN( KL+KV+1,M ),N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CDBTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $     RETURN
*
*     Determine the block size for this environment
*
      NB = ILAENV( 1, 'CDBTRF', ' ', M, N, KL, KU )
*
*     The block size must not exceed the limit set by the size of the
*     local arrays WORK13 and WORK31.
*
      NB = MIN( NB, NBMAX )
*
      IF( NB.LE.1 .OR. NB.GT.KL ) THEN
*
*        Use unblocked code
*
         CALL CDBTF2( M, N, KL, KU, AB, LDAB, INFO )
      ELSE
*
*        Use blocked code
*
*        Zero the superdiagonal elements of the work array WORK13
*
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
*
*        Zero the subdiagonal elements of the work array WORK31
*
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
*
*        JU is the index of the last column affected by the current
*        stage of the factorization
*
         JU = 1
*
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
*
*           The active part of the matrix is partitioned
*
*              A11   A12   A13
*              A21   A22   A23
*              A31   A32   A33
*
*           Here A11, A21 and A31 denote the current block of JB columns
*           which is about to be factorized. The number of rows in the
*           partitioning are JB, I2, I3 respectively, and the numbers
*           of columns are JB, J2, J3. The superdiagonal elements of A13
*           and the subdiagonal elements of A31 lie outside the band.
*
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
*
*           J2 and J3 are computed after JU has been updated.
*
*           Factorize the current block of JB columns
*
            DO 80 JJ = J, J + JB - 1
*
*              Find pivot and test for singularity. KM is the number of
*              subdiagonal elements in the current column.
*
               KM = MIN( KL, M-JJ )
               JP = 1
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
*
*                 Compute multipliers
*
                  CALL CSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ),
     $                 1 )
*
*                 Update trailing submatrix within the band and within
*                 the current block. JM is the index of the last column
*                 which needs to be updated.
*
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ ) THEN
                     CALL CGERU( KM, JM-JJ, -CONE, AB( KV+2, JJ ), 1,
     $                          AB( KV, JJ+1 ), LDAB-1,
     $                          AB( KV+1, JJ+1 ), LDAB-1 )
                  END IF
               END IF
*
*              Copy current column of A31 into the work array WORK31
*
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 )
     $            CALL CCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1,
     $                        WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
*
*              Apply the row interchanges to the other blocks.
*
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
*
*              Update the relevant part of the trailing submatrix
*
               IF( J2.GT.0 ) THEN
*
*                 Update A12
*
                  CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J2, CONE, AB( KV+1, J ), LDAB-1,
     $                        AB( KV+1-JB, J+JB ), LDAB-1 )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A22
*
                     CALL CGEMM( 'No transpose', 'No transpose', I2, J2,
     $                           JB, -CONE, AB( KV+1+JB, J ), LDAB-1,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, CONE,
     $                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Update A32
*
                     CALL CGEMM( 'No transpose', 'No transpose', I3, J2,
     $                           JB, -CONE, WORK31, LDWORK,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, CONE,
     $                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
*
               IF( J3.GT.0 ) THEN
*
*                 Copy the lower triangle of A13 into the work array
*                 WORK13
*
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
*
*                 Update A13 in the work array
*
                  CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J3, CONE, AB( KV+1, J ), LDAB-1,
     $                        WORK13, LDWORK )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A23
*
                     CALL CGEMM( 'No transpose', 'No transpose', I2, J3,
     $                           JB, -CONE, AB( KV+1+JB, J ), LDAB-1,
     $                           WORK13, LDWORK, CONE, AB( 1+JB, J+KV ),
     $                           LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Update A33
*
                     CALL CGEMM( 'No transpose', 'No transpose', I3, J3,
     $                         JB, -CONE, WORK31, LDWORK, WORK13,
     $                         LDWORK, CONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
*
*                 Copy the lower triangle of A13 back into place
*
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
            END IF
*
*           copy the upper triangle of A31 back into place
*
            DO 170 JJ = J + JB - 1, J, -1
*
*              Copy the current column of A31 back into place
*
               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 )
     $            CALL CCOPY( NW, WORK31( 1, JJ-J+1 ), 1,
     $                        AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
*
      RETURN
*
*     End of CDBTRF
*
      END
