      SUBROUTINE CRSHFT( M, N, OFFSET, A, LDA )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            LDA, M, N, OFFSET
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  CRSHFT shifts rows of an m by n array A by OFFSET.
*
*  Arguments
*  =========
*
*  M       (local input) INTEGER
*          On entry, M  specifies the number of rows of A to be shifted.
*          M must be at least zero.
*
*  N       (local input) INTEGER
*          On entry,  N  specifies the number of columns of A. N must be
*          at least zero.
*
*  OFFSET  (local input) INTEGER
*          On entry, OFFSET  specifies  the  offset by which the rows of
*          A should be shifted. OFFSET  can be positive or negative (see
*          below for further details). When OFFSET is positive, the rows
*          are shifted to the  bottom. When OFFSET is negative, the rows
*          of A are shifted to the top.
*
*  A       (local input/local output) COMPLEX array
*          On entry,  A  is an array of dimension ( LDA, N ). On exit, A
*          contains the shifted array.
*
*  LDA     (local input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least max( 1, M+ABS(OFFSET) ).
*
*  Further Details
*  ===============
*
*            N                 N                   N                 N
*           ---               ---                 ---               ---
*          | 1 |             | 1 |               | 1 |             | 7 |
*          | 2 |    M = 3    | 2 |               | 2 |    M = 3    | 8 |
*          | 3 |             | 3 |               | 3 |             | 9 |
*          | 4 |             | 4 |               | 4 |             | 4 |
*          | 5 |   becomes   | 5 |               | 5 |   becomes   | 5 |
*          | 6 |             | 6 |               | 6 |             | 6 |
*          | 7 |             | 1 |               | 7 |             | 7 |
*          | 8 | OFFSET =  6 | 2 |               | 8 | OFFSET = -6 | 8 |
*          | 9 |             | 3 |               | 9 |             | 9 |
*           ---               ---                 ---               ---
*                OFFSET >= 0                           OFFSET <= 0
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. Executable Statements ..
*
      IF( ( OFFSET.EQ.0 ).OR.( M.LE.0 ).OR.( N.LE.0 ) )
     $   RETURN
*
      IF( OFFSET.GT.0 ) THEN
         DO 20 J = 1, N
            DO 10 I = M, 1, -1
               A( I+OFFSET, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE
         DO 40 J = 1, N
            DO 30 I = 1, M
               A( I, J ) = A( I-OFFSET, J )
   30       CONTINUE
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of CRSHFT
*
      END
