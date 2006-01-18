      SUBROUTINE DTZPAD( UPLO, HERM, M, N, IOFFD, ALPHA, BETA, A, LDA )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        HERM, UPLO
      INTEGER            IOFFD, LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DTZPAD  initializes a two-dimensional array A to beta on the diagonal
*  specified by IOFFD or zeros the imaginary part of those diagonals and
*  set the offdiagonals to alpha.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          On entry,  UPLO  specifies  which trapezoidal part of the ar-
*          ray A is to be set as follows:
*             = 'L' or 'l':   Lower triangular part is set; the strictly
*                             upper triangular part of A is not changed,
*             = 'D' or 'd':   diagonal  specified  by  IOFFD is set; the
*                             rest of the array A is unchanged,
*             = 'U' or 'u':   Upper triangular part is set; the strictly
*                             lower triangular part of A is not changed,
*             Otherwise:      All of the array A is set.
*
*  HERM    (input) CHARACTER*1
*          On entry, HERM specifies what should be done to the diagonals
*          as follows.  When UPLO is 'L', 'l', 'D', 'd', 'U' or 'u'  and
*          HERM is  'Z'  or  'z', the imaginary part of the diagonals is
*          set  to  zero. Otherwise, the diagonals are set to beta.
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
*  ALPHA   (input) DOUBLE PRECISION
*          On entry,  ALPHA  specifies the scalar alpha, i.e., the value
*          to which the offdiagonal entries of the array A determined by
*          UPLO and IOFFD are set.
*
*  BETA    (input) DOUBLE PRECISION
*          On entry, BETA  specifies the scalar beta, i.e., the value to
*          which the diagonal entries specified by IOFFD of the array  A
*          are set. BETA is not referenced when UPLO is 'L', 'l', 'U' or
*          'u' and HERM is 'Z'.
*
*  A       (input/output) DOUBLE PRECISION array
*          On entry, A is an array of dimension  (LDA,N).  Before  entry
*          with UPLO = 'U', the leading m by n part of the array  A must
*          contain the upper trapezoidal part of the matrix to be set as
*          specified by  IOFFD,  and the strictly lower trapezoidal part
*          of A is not referenced;  When  UPLO = 'L', the leading m by n
*          part of the array A must contain the lower  trapezoidal  part
*          of  the  matrix  to  be  set  as  specified by IOFFD, and the
*          strictly upper  trapezoidal  part of A is not referenced.  On
*          exit, the entries  of the  trapezoid  part of A determined by
*          UPLO, HERM and IOFFD are set.
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
*     .. Local Scalars ..
      INTEGER            I, J, JTMP, MN
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
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
*     Start the operations
*
      IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Set the diagonal to BETA or zero the imaginary part of the
*        diagonals and set the strictly lower triangular part of the
*        array to ALPHA.
*
         MN = MAX( 0, -IOFFD )
         DO 20 J = 1, MIN( MN, N )
            DO 10 I = 1, M
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
*
         IF( LSAME( HERM, 'Z' ) ) THEN
            DO 40 J = MN + 1, MIN( M - IOFFD, N )
               JTMP = J + IOFFD
               DO 30 I = JTMP + 1, M
                  A( I, J ) = ALPHA
   30          CONTINUE
   40       CONTINUE
         ELSE
            DO 60 J = MN + 1, MIN( M - IOFFD, N )
               JTMP = J + IOFFD
               A( JTMP, J ) = BETA
               DO 50 I = JTMP + 1, M
                  A( I, J ) = ALPHA
   50          CONTINUE
   60       CONTINUE
         END IF
*
      ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Set the diagonal to BETA or zero the imaginary part of the
*        diagonals and set the strictly upper triangular part of the
*        array to ALPHA.
*
         MN = MIN( M - IOFFD, N )
         IF( LSAME( HERM, 'Z' ) ) THEN
            DO 80 J = MAX( 0, -IOFFD ) + 1, MN
               JTMP = J + IOFFD
               DO 70 I = 1, JTMP - 1
                  A( I, J ) = ALPHA
   70          CONTINUE
   80       CONTINUE
         ELSE
            DO 100 J = MAX( 0, -IOFFD ) + 1, MN
               JTMP = J + IOFFD
               DO 90 I = 1, JTMP - 1
                  A( I, J ) = ALPHA
   90          CONTINUE
               A( JTMP, J ) = BETA
  100       CONTINUE
         END IF
         DO 120 J = MAX( 0, MN ) + 1, N
            DO 110 I = 1, M
               A( I, J ) = ALPHA
  110       CONTINUE
  120    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'D' ) ) THEN
*
*        Set the diagonal to BETA
*
         IF( .NOT.( LSAME( HERM, 'Z' ) ) ) THEN
            IF( ( IOFFD.LT.M ).AND.( IOFFD.GT.-N ) ) THEN
               DO 130 J = MAX( 0, -IOFFD ) + 1, MIN( M - IOFFD, N )
                  A( J + IOFFD, J ) = BETA
  130          CONTINUE
            END IF
         END IF
*
      ELSE
*
*        Set the diagonals to BETA and the offdiagonals to ALPHA.
*
         DO 150 J = 1, N
            DO 140 I = 1, M
               A( I, J ) = ALPHA
  140       CONTINUE
  150    CONTINUE
         IF( ALPHA.NE.BETA .AND. IOFFD.LT.M .AND. IOFFD.GT.-N ) THEN
            DO 160 J = MAX( 0, -IOFFD ) + 1, MIN( M - IOFFD, N )
               A( J + IOFFD, J ) = BETA
  160       CONTINUE
         END IF
*
      END IF
*
      RETURN
*
*     End of DTZPAD
*
      END
