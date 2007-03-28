      SUBROUTINE SLASORTE( S, LDS, J, OUT, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      INTEGER            INFO, J, LDS
*     ..
*     .. Array Arguments ..
      REAL               OUT( J, * ), S( LDS, * )
*     ..
*
*  Purpose
*  =======
*
*  SLASORTE sorts eigenpairs so that real eigenpairs are together and
*    complex are together.  This way one can employ 2x2 shifts easily
*    since every 2nd subdiagonal is guaranteed to be zero.
*  This routine does no parallel work.
*
*  Arguments
*  =========
*
*  S       (local input/output) REAL array, dimension LDS
*          On entry, a matrix already in Schur form.
*          On exit, the diagonal blocks of S have been rewritten to pair
*             the eigenvalues.  The resulting matrix is no longer
*             similar to the input.
*
*  LDS     (local input) INTEGER
*          On entry, the leading dimension of the local array S.
*          Unchanged on exit.
*
*  J       (local input) INTEGER
*          On entry, the order of the matrix S.
*          Unchanged on exit.
*
*  OUT     (local input/output) REAL array, dimension Jx2
*          This is the work buffer required by this routine.
*
*  INFO    (local input) INTEGER
*          This is set if the input matrix had an odd number of real
*          eigenvalues and things couldn't be paired or if the input
*           matrix S was not originally in Schur form.
*          0 indicates successful completion.
*
*  Implemented by:  G. Henry, November 17, 1996
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            BOT, I, LAST, TOP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
      LAST = J
      TOP = 1
      BOT = J
      INFO = 0
      DO 10 I = J - 1, 1, -1
         IF( S( I+1, I ).EQ.ZERO ) THEN
            IF( LAST-I.EQ.2 ) THEN
               OUT( BOT-1, 1 ) = S( I+1, I+1 )
               OUT( BOT, 2 ) = S( I+2, I+2 )
               OUT( BOT-1, 2 ) = S( I+1, I+2 )
               OUT( BOT, 1 ) = S( I+2, I+1 )
               BOT = BOT - 2
            END IF
            IF( LAST-I.EQ.1 ) THEN
               IF( MOD( TOP, 2 ).EQ.1 ) THEN
*
*                 FIRST OF A PAIR
*
                  IF( ( I.EQ.J-1 ) .OR. ( I.EQ.1 ) ) THEN
                     OUT( TOP, 1 ) = S( I+1, I+1 )
                  ELSE
                     OUT( TOP, 1 ) = S( I+1, I+1 )
                  END IF
                  OUT( TOP, 2 ) = ZERO
               ELSE
*
*                 SECOND OF A PAIR
*
                  IF( ( I.EQ.J-1 ) .OR. ( I.EQ.1 ) ) THEN
                     OUT( TOP, 2 ) = S( I+1, I+1 )
                  ELSE
                     OUT( TOP, 2 ) = S( I+1, I+1 )
                  END IF
                  OUT( TOP, 1 ) = ZERO
               END IF
               TOP = TOP + 1
            END IF
            IF( LAST-I.GT.2 ) THEN
               INFO = I
               RETURN
            END IF
            LAST = I
         END IF
   10 CONTINUE
      IF( LAST.EQ.2 ) THEN
*
*        GRAB LAST DOUBLE PAIR
*
         OUT( BOT-1, 1 ) = S( 1, 1 )
         OUT( BOT, 2 ) = S( 2, 2 )
         OUT( BOT-1, 2 ) = S( 1, 2 )
         OUT( BOT, 1 ) = S( 2, 1 )
         BOT = BOT - 2
      END IF
      IF( LAST.EQ.1 .and. mod(top, 2) .eq. 0 ) THEN
*
*        GRAB SECOND PART OF LAST PAIR
*
         OUT(TOP, 2) = s(1,1)
         OUT(TOP, 1) = zero
         TOP = TOP + 1
      END IF
      IF( TOP-1.NE.BOT ) THEN
         INFO = -BOT
         RETURN
      END IF
*
*     Overwrite the S diagonals
*
      DO 20 I = 1, J, 2
         S( I, I ) = OUT( I, 1 )
         S( I+1, I ) = OUT( I+1, 1 )
         S( I, I+1 ) = OUT( I, 2 )
         S( I+1, I+1 ) = OUT( I+1, 2 )
   20 CONTINUE
*
      RETURN
*
*     End of SLASORTE
*
      END
