      SUBROUTINE ZTZSCAL( UPLO, M, N, IOFFD, ALPHA, A, LDA )
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
*  ZTZSCAL scales a two-dimensional array A by the scalar alpha.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          On entry,  UPLO  specifies  which trapezoidal part of the ar-
*          ray A is to be scaled as follows:
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
*          determined by UPLO and IOFFD are scaled.
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
      INTEGER            J, JTMP, MN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZSCAL, ZTZPAD
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( ( M.LE.0 ).OR.( N.LE.0 ).OR.( ALPHA.EQ.ONE ) )
     $   RETURN
*
*     Start the operations
*
      IF( ALPHA.EQ.ZERO ) THEN
         CALL ZTZPAD( UPLO, 'N', M, N, IOFFD, ZERO, ZERO, A, LDA )
         RETURN
      END IF
*
      IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Scales the lower triangular part of the array by ALPHA.
*
         MN = MAX( 0, -IOFFD )
         DO 10 J = 1, MIN( MN, N )
            CALL ZSCAL( M, ALPHA, A( 1, J ), 1 )
   10    CONTINUE
         DO 20 J = MN + 1, MIN( M - IOFFD, N )
            JTMP = J + IOFFD
            IF( M.GE.JTMP )
     $         CALL ZSCAL( M-JTMP+1, ALPHA, A( JTMP, J ), 1 )
   20    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Scales the upper triangular part of the array by ALPHA.
*
         MN = MIN( M - IOFFD, N )
         DO 30 J = MAX( 0, -IOFFD ) + 1, MN
            CALL ZSCAL( J + IOFFD, ALPHA, A( 1, J ), 1 )
   30    CONTINUE
         DO 40 J = MAX( 0, MN ) + 1, N
            CALL ZSCAL( M, ALPHA, A( 1, J ), 1 )
   40    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'D' ) ) THEN
*
*        Scales the diagonal entries by ALPHA.
*
         DO 50 J = MAX( 0, -IOFFD ) + 1, MIN( M - IOFFD, N )
            JTMP = J + IOFFD
            A( JTMP, J ) = ALPHA * A( JTMP, J )
   50    CONTINUE
*
      ELSE
*
*        Scales the entire array by ALPHA.
*
         DO 60 J = 1, N
            CALL ZSCAL( M, ALPHA, A( 1, J ), 1 )
   60    CONTINUE
*
      END IF
*
      RETURN
*
*     End of ZTZSCAL
*
      END
