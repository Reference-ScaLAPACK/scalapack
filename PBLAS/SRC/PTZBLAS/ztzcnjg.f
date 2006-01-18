      SUBROUTINE ZTZCNJG( UPLO, M, N, IOFFD, ALPHA, A, LDA )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO
      INTEGER            IOFFD, LDA, M, N
      COMPLEX*16         ALPHA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZTZCNJG  conjugates  a  two-dimensional array A and then scales it by
*  the scalar alpha.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          On entry,  UPLO  specifies  which trapezoidal part of the ar-
*          ray A is to be conjugated and scaled as follows:
*             = 'L' or 'l':          the lower trapezoid of A is scaled,
*             = 'U' or 'u':          the upper trapezoid of A is scaled,
*             = 'D' or 'd':       diagonal specified by IOFFD is scaled,
*             Otherwise:                   all of the array A is scaled.
*
*  M       (input) INTEGER
*          On entry,  M  specifies the number of rows of the array A.  M
*          must be at least zero.
*
*  N       (input) INTEGER
*          On entry,  N  specifies the number of columns of the array A.
*          N must be at least zero.
*
*  IOFFD   (input) INTEGER
*          On entry, IOFFD specifies the position of the offdiagonal de-
*          limiting the upper and lower trapezoidal part of A as follows
*          (see the notes below):
*
*             IOFFD = 0  specifies the main diagonal A( i, i ),
*                        with i = 1 ... MIN( M, N ),
*             IOFFD > 0  specifies the subdiagonal   A( i+IOFFD, i ),
*                        with i = 1 ... MIN( M-IOFFD, N ),
*             IOFFD < 0  specifies the superdiagonal A( i, i-IOFFD ),
*                        with i = 1 ... MIN( M, N+IOFFD ).
*
*  ALPHA   (input) COMPLEX*16
*          On entry,  ALPHA  specifies the scalar alpha, i.e., the value
*          by which the diagonal and offdiagonal entries of the array  A
*          as specified by UPLO and IOFFD are scaled.
*
*  A       (input/output) COMPLEX*16 array
*          On entry, A is an array of dimension  (LDA,N).  Before  entry
*          with  UPLO = 'U' or 'u', the leading m by n part of the array
*          A must contain the upper trapezoidal  part  of the matrix  as
*          specified by  IOFFD to be scaled, and the strictly lower tra-
*          pezoidal part of A is not referenced; When UPLO = 'L' or 'l',
*          the leading m by n part of the array A must contain the lower
*          trapezoidal  part  of  the matrix as specified by IOFFD to be
*          scaled,  and  the strictly upper trapezoidal part of A is not
*          referenced. On exit, the entries of the  trapezoid part of  A
*          determined by UPLO and IOFFD are conjugated and scaled.
*
*  LDA     (input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least max( 1, M ).
*
*  Notes
*  =====
*                           N                                    N
*             ----------------------------                  -----------
*            |       d                    |                |           |
*          M |         d        'U'       |                |      'U'  |
*            |  'L'     'D'               |                |d          |
*            |             d              |              M |  d        |
*             ----------------------------                 |   'D'     |
*                                                          |      d    |
*               IOFFD < 0                                  | 'L'    d  |
*                                                          |          d|
*                  N                                       |           |
*             -----------                                   -----------
*            |    d   'U'|
*            |      d    |                                   IOFFD > 0
*          M |       'D' |
*            |          d|                              N
*            |  'L'      |                 ----------------------------
*            |           |                |          'U'               |
*            |           |                |d                           |
*            |           |                | 'D'                        |
*            |           |                |    d                       |
*            |           |                |'L'   d                     |
*             -----------                  ----------------------------
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JTMP, MN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZTZPAD
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( ( M.LE.0 ).OR.( N.LE.0 ) )
     $   RETURN
*
*     Start the operations
*
      IF( ALPHA.EQ.ZERO ) THEN
*
         CALL ZTZPAD( UPLO, 'N', M, N, IOFFD, ZERO, ZERO, A, LDA )
*
      ELSE IF( ALPHA.EQ.ONE ) THEN
*
         IF( LSAME( UPLO, 'L' ) ) THEN
*
            MN = MAX( 0, -IOFFD )
            DO 20 J = 1, MIN( MN, N )
               DO 10 I = 1, M
                  A( I, J ) = DCONJG( A( I, J ) )
   10          CONTINUE
   20       CONTINUE
*
            DO 40 J = MN + 1, MIN( M - IOFFD, N )
               DO 30 I = J + IOFFD, M
                  A( I, J ) = DCONJG( A( I, J ) )
   30          CONTINUE
   40       CONTINUE
*
         ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Scales the upper triangular part of the array by ALPHA.
*
            MN = MIN( M - IOFFD, N )
            DO 60 J = MAX( 0, -IOFFD ) + 1, MN
               DO 50 I = 1, J + IOFFD
                  A( I, J ) = DCONJG( A( I, J ) )
   50          CONTINUE
   60       CONTINUE
*
            DO 80 J = MAX( 0, MN ) + 1, N
               DO 70 I = 1, M
                  A( I, J ) = DCONJG( A( I, J ) )
   70          CONTINUE
   80       CONTINUE
*
         ELSE IF( LSAME( UPLO, 'D' ) ) THEN
*
*           Scales the diagonal entries by ALPHA.
*
            DO 90 J = MAX( 0, -IOFFD ) + 1, MIN( M - IOFFD, N )
               JTMP = J + IOFFD
               A( JTMP, J ) = DCONJG( A( JTMP, J ) )
   90       CONTINUE
*
         ELSE
*
*           Scales the entire array by ALPHA.
*
            DO 110 J = 1, N
               DO 100 I = 1, M
                  A( I, J ) = DCONJG( A( I, J ) )
  100          CONTINUE
  110       CONTINUE
*
         END IF
*
      ELSE
*
         IF( LSAME( UPLO, 'L' ) ) THEN
*
*           Scales the lower triangular part of the array by ALPHA.
*
            MN = MAX( 0, -IOFFD )
            DO 130 J = 1, MIN( MN, N )
               DO 120 I = 1, M
                  A( I, J ) = ALPHA * DCONJG( A( I, J ) )
  120          CONTINUE
  130       CONTINUE
*
            DO 150 J = MN + 1, MIN( M - IOFFD, N )
               DO 140 I = J + IOFFD, M
                  A( I, J ) = ALPHA * DCONJG( A( I, J ) )
  140          CONTINUE
  150       CONTINUE
*
         ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Scales the upper triangular part of the array by ALPHA.
*
            MN = MIN( M - IOFFD, N )
            DO 170 J = MAX( 0, -IOFFD ) + 1, MN
               DO 160 I = 1, J + IOFFD
                  A( I, J ) = ALPHA * DCONJG( A( I, J ) )
  160          CONTINUE
  170       CONTINUE
*
            DO 190 J = MAX( 0, MN ) + 1, N
               DO 180 I = 1, M
                  A( I, J ) = ALPHA * DCONJG( A( I, J ) )
  180          CONTINUE
  190       CONTINUE
*
         ELSE IF( LSAME( UPLO, 'D' ) ) THEN
*
*           Scales the diagonal entries by ALPHA.
*
            DO 200 J = MAX( 0, -IOFFD ) + 1, MIN( M - IOFFD, N )
               JTMP = J + IOFFD
               A( JTMP, J ) = ALPHA * DCONJG( A( JTMP, J ) )
  200       CONTINUE
*
         ELSE
*
*           Scales the entire array by ALPHA.
*
            DO 220 J = 1, N
               DO 210 I = 1, M
                  A( I, J ) = ALPHA * DCONJG( A( I, J ) )
  210          CONTINUE
  220       CONTINUE
*
         END IF
*
      END IF
*
      RETURN
*
*     End of ZTZCNJG
*
      END
