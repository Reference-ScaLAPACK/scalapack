      SUBROUTINE PDFILLPAD( ICTXT, M, N, A, LDA, IPRE, IPOST, CHKVAL )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, IPOST, IPRE, LDA, M, N
      DOUBLE PRECISION   CHKVAL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( * )
*     ..
*
*  Purpose
*  =======
*
*  PDFILLPAD surrounds a two dimensional local array with a guard-
*  zone initialized to the value CHKVAL. The user may later call the
*  routine PDCHEKPAD to discover if the guardzone has been
*  violated. There are three guardzones. The first is a buffer of size
*  IPRE that is before the start of the array. The second is the buffer
*  of size IPOST which is after the end of the array to be padded.
*  Finally, there is a guard zone inside every column of the array to
*  be padded, in the elements of A(M+1:LDA, J).
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation. The context itself is global.
*
*  M       (local input) INTEGER
*          The number of rows in the local array.
*
*  N       (local input) INTEGER
*          The number of columns in the local array.
*
*  A       (local input/local output) DOUBLE PRECISION, array of
*          dimension (LDA,N). A location IPRE elements in front of
*          the matrix to be padded.
*
*  LDA     (local input) INTEGER
*          The leading Dimension of the local array to be padded.
*
*  IPRE    (local input) INTEGER
*          The size of the guard zone to put before the start of
*          padded array.
*
*  IPOST   (local input) INTEGER
*          The size of the guard zone to put after padded array.
*
*  CHKVAL  (local input) DOUBLE PRECISION
*          The value to pad matrix with.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J, K
*     ..
*     .. Executable Statements ..
*
*     Put check buffer in front of A
*
      IF( IPRE.GT.0 ) THEN
         DO 10 I = 1, IPRE
            A( I ) = CHKVAL
   10    CONTINUE
      ELSE
         WRITE( *, FMT = * ) 'WARNING no pre-guardzone in PDFILLPAD'
      END IF
*
*     Put check buffer in back of A
*
      IF( IPOST.GT.0 ) THEN
         J = IPRE+LDA*N+1
         DO 20 I = J, J+IPOST-1
            A( I ) = CHKVAL
   20    CONTINUE
      ELSE
         WRITE( *, FMT = * ) 'WARNING no post-guardzone in PDFILLPAD'
      END IF
*
*     Put check buffer in all (LDA-M) gaps
*
      IF( LDA.GT.M ) THEN
         K = IPRE + M + 1
         DO 40 J = 1, N
            DO 30 I = K, K + (LDA-M) - 1
               A( I ) = CHKVAL
   30       CONTINUE
            K = K + LDA
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of PDFILLPAD
*
      END
