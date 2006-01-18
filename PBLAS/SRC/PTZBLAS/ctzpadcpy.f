      SUBROUTINE CTZPADCPY( UPLO, DIAG, M, N, IOFFD, A, LDA, B, LDB )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        DIAG, UPLO
      INTEGER            IOFFD, LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  CTZPADCPY copies an array A into an array B.  The unchanged part of B
*  is padded with zeros. The diagonal of B specified by IOFFD may be set
*  to ones.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          On entry,  UPLO  specifies  which trapezoidal part of the ar-
*          ray A is to be copied as follows:
*             = 'L' or 'l':   Lower  triangular  part  is  copied;   the
*                             strictly  upper  triangular  part  of B is
*                             padded with zeros,
*             = 'U' or 'u':   Upper  triangular  part  is  copied;   the
*                             strictly  lower  triangular  part  of B is
*                             padded with zeros.
*
*  DIAG    (input) CHARACTER*1
*          On entry, DIAG specifies whether or not the diagonal of B  is
*          to be set to ones or not as follows:
*
*          DIAG = 'N' or 'n': the diagonals of A  are  copied  into  the
*          diagonals of B, otherwise the diagonals of B are set to ones.
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
*  A       (input) COMPLEX array
*          On entry, A is an array of dimension  (LDA,N).  Before  entry
*          with UPLO = 'U', the leading m by n part of the array  A must
*          contain the upper trapezoidal part of the matrix to be copied
*          as specified by IOFFD, UPLO and DIAG, and  the strictly lower
*          trapezoidal part of A is not referenced; When  UPLO = 'L',the
*          leading m by n part of the array  A  must contain  the  lower
*          trapezoidal part of  the  matrix to be copied as specified by
*          IOFFD, UPLO and DIAG and the strictly upper trapezoidal  part
*          of A is not referenced.
*
*  LDA     (input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least max( 1, M ).
*
*  B       (output) COMPLEX array
*          On entry, B  is  an array of dimension (LDB,N). On exit, this
*          array  contains  the  padded copy of A as specified by IOFFD,
*          UPLO and DIAG.
*
*  LDB     (input) INTEGER
*          On entry, LDB specifies the leading dimension of the array B.
*          LDB must be at least max( 1, M ).
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
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ITMP, J, JTMP, MN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
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
         MN = MAX( 0, -IOFFD )
         DO 20 J = 1, MIN( MN, N )
            DO 10 I = 1, M
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
*
         JTMP = MIN( M - IOFFD, N )
*
         IF( LSAME( DIAG, 'N' ) ) THEN
            DO 50 J = MN + 1, JTMP
               ITMP = J + IOFFD
               DO 30 I = 1, ITMP - 1
                  B( I, J ) = ZERO
   30          CONTINUE
               DO 40 I = ITMP, M
                  B( I, J ) = A( I, J )
   40          CONTINUE
   50       CONTINUE
         ELSE
            DO 80 J = MN + 1, JTMP
               ITMP = J + IOFFD
               DO 60 I = 1, ITMP - 1
                  B( I, J ) = ZERO
   60          CONTINUE
               B( ITMP, J ) = ONE
               DO 70 I = ITMP + 1, M
                  B( I, J ) = A( I, J )
   70          CONTINUE
   80       CONTINUE
         END IF
*
         DO 100 J = JTMP + 1, N
            DO 90 I = 1, M
               B( I, J ) = ZERO
   90       CONTINUE
  100    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
         JTMP = MAX( 0, -IOFFD )
*
         DO 120 J = 1, JTMP
            DO 110 I = 1, M
               B( I, J ) = ZERO
  110       CONTINUE
  120    CONTINUE
*
         MN = MIN( M - IOFFD, N )
*
         IF( LSAME( DIAG, 'N' ) ) THEN
            DO 150 J = JTMP + 1, MN
               ITMP = J + IOFFD
               DO 130 I = 1, ITMP
                  B( I, J ) = A( I, J )
  130          CONTINUE
               DO 140 I = ITMP + 1, M
                  B( I, J ) = ZERO
  140          CONTINUE
  150       CONTINUE
         ELSE
            DO 180 J = JTMP + 1, MN
               ITMP = J + IOFFD
               DO 160 I = 1, ITMP - 1
                  B( I, J ) = A( I, J )
  160          CONTINUE
               B( ITMP, J ) = ONE
               DO 170 I = ITMP + 1, M
                  B( I, J ) = ZERO
  170          CONTINUE
  180       CONTINUE
         END IF
*
         DO 200 J = MAX( 0, MN ) + 1, N
            DO 190 I = 1, M
               B( I, J ) = A( I, J )
  190       CONTINUE
  200    CONTINUE
*
      ELSE
*
         DO 220 J = 1, N
            DO 210 I = 1, M
               B( I, J ) = A( I, J )
  210       CONTINUE
  220    CONTINUE
*
      END IF
*
      RETURN
*
*     End of CTZPADCPY
*
      END
