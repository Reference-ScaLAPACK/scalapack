      SUBROUTINE SLAPST( ID, N, D, INDX, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      INTEGER            INDX( * )
      REAL               D( * )
*     ..
*
*  Purpose
*  =======
*  SLAPST is a modified version of the LAPACK routine SLASRT.
*
*  Define a permutation INDX that sorts the numbers in D
*  in increasing order (if ID = 'I') or
*  in decreasing order (if ID = 'D' ).
*
*  Use Quick Sort, reverting to Insertion sort on arrays of
*  size <= 20. Dimension of STACK limits N to about 2**32.
*
*  Arguments
*  =========
*
*  ID      (input) CHARACTER*1
*          = 'I': sort D in increasing order;
*          = 'D': sort D in decreasing order.
*
*  N       (input) INTEGER
*          The length of the array D.
*
*  D       (input)  REAL array, dimension (N)
*          The array to be sorted.
*
*  INDX    (ouput) INTEGER array, dimension (N).
*          The permutation which sorts the array D.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
*     ..
*     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, ITMP, J, START, STKPNT
      REAL               D1, D2, D3, DMNMX
*     ..
*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input paramters.
*
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAPST', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      DO 10 I = 1, N
         INDX( I ) = I
   10 CONTINUE
*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   20 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
*
*        Do Insertion sort on D( START:ENDD )
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            DO 40 I = START + 1, ENDD
               DO 30 J = I, START + 1, -1
                  IF( D( INDX( J ) ).GT.D( INDX( J-1 ) ) ) THEN
                     ITMP = INDX( J )
                     INDX( J ) = INDX( J-1 )
                     INDX( J-1 ) = ITMP
                  ELSE
                     GO TO 40
                  END IF
   30          CONTINUE
   40       CONTINUE
*
         ELSE
*
*           Sort into increasing order
*
            DO 60 I = START + 1, ENDD
               DO 50 J = I, START + 1, -1
                  IF( D( INDX( J ) ).LT.D( INDX( J-1 ) ) ) THEN
                     ITMP = INDX( J )
                     INDX( J ) = INDX( J-1 )
                     INDX( J-1 ) = ITMP
                  ELSE
                     GO TO 60
                  END IF
   50          CONTINUE
   60       CONTINUE
*
         END IF
*
      ELSE IF( ENDD-START.GT.SELECT ) THEN
*
*        Partition D( START:ENDD ) and stack parts, largest one first
*
*        Choose partition entry as median of 3
*
         D1 = D( INDX( START ) )
         D2 = D( INDX( ENDD ) )
         I = ( START+ENDD ) / 2
         D3 = D( INDX( I ) )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            I = START - 1
            J = ENDD + 1
   70       CONTINUE
   80       CONTINUE
            J = J - 1
            IF( D( INDX( J ) ).LT.DMNMX )
     $         GO TO 80
   90       CONTINUE
            I = I + 1
            IF( D( INDX( I ) ).GT.DMNMX )
     $         GO TO 90
            IF( I.LT.J ) THEN
               ITMP = INDX( I )
               INDX( I ) = INDX( J )
               INDX( J ) = ITMP
               GO TO 70
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
*
*           Sort into increasing order
*
            I = START - 1
            J = ENDD + 1
  100       CONTINUE
  110       CONTINUE
            J = J - 1
            IF( D( INDX( J ) ).GT.DMNMX )
     $         GO TO 110
  120       CONTINUE
            I = I + 1
            IF( D( INDX( I ) ).LT.DMNMX )
     $         GO TO 120
            IF( I.LT.J ) THEN
               ITMP = INDX( I )
               INDX( I ) = INDX( J )
               INDX( J ) = ITMP
               GO TO 100
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 )
     $   GO TO 20
      RETURN
*
*     End of SLAPST
*
      END
