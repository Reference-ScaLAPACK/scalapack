      SUBROUTINE DCSHFT( M, N, OFFSET, A, LDA )
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
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DCSHFT shifts columns of an m by n array A by OFFSET.
*
*  Arguments
*  =========
*
*  M       (local input) INTEGER
*          On entry,  M  specifies the number of rows of A. M must be at
*          least zero.
*
*  N       (local input) INTEGER
*          On entry,  N  specifies  the  number of columns of  A  to  be
*          shifted. N must be at least zero.
*
*  OFFSET  (local input) INTEGER
*          On entry, OFFSET specifies the offset by which the columns of
*          A should be shifted. OFFSET  can be positive or negative (see
*          below for further details). When OFFSET  is positive, the co-
*          lumns are shifted to the right. When OFFSET  is negative, the
*          columns of A are shifted to the left.
*
*  A       (local input/local output) DOUBLE PRECISION array
*          On entry, A  is an array of dimension ( LDA, N+ABS(OFFSET) ).
*          On exit, A contains the shifted array.
*
*  LDA     (local input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least max( 1, M ).
*
*  Further Details
*  ===============
*
*             N=3     OFFSET=6                      -OFFSET=6     N=3
*           -------------------                     -------------------
*          | 1 2 3 4 5 6 7 8 9 |         M         | 1 2 3 4 5 6 7 8 9 |
*           -------------------                     -------------------
*                    V                                       V
*           -------------------                     -------------------
*          | 1 2 3 4 5 6 1 2 3 |         M         | 7 8 9 4 5 6 7 8 9 |
*           -------------------                     -------------------
*               OFFSET >= 0                             OFFSET <= 0
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
         DO 20 J = N, 1, -1
            DO 10 I = 1, M
               A( I, J+OFFSET ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE
         DO 40 J = 1, N
            DO 30 I = 1, M
               A( I, J ) = A( I, J-OFFSET )
   30       CONTINUE
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of DCSHFT
*
      END
