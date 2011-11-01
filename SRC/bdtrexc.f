      SUBROUTINE BDTREXC( N, T, LDT, IFST, ILST, NITRAF, ITRAF,
     $                    NDTRAF, DTRAF, WORK, INFO )
      IMPLICIT NONE
*
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            IFST, ILST, INFO, LDT, N, NDTRAF, NITRAF
*     ..
*     .. Array Arguments ..
      INTEGER            ITRAF( * )
      DOUBLE PRECISION   DTRAF( * ), T( LDT, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  BDTREXC reorders the real Schur factorization of a real matrix
*  A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
*  moved to row ILST.
*
*  The real Schur form T is reordered by an orthogonal similarity
*  transformation Z**T*T*Z. In contrast to the LAPACK routine DTREXC,
*  the orthogonal matrix Z is not explicitly constructed but
*  represented by paramaters contained in the arrays ITRAF and DTRAF,
*  see further details.
*
*  T must be in Schur canonical form (as returned by DHSEQR), that is,
*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
*  2-by-2 diagonal block has its diagonal elements equal and its
*  off-diagonal elements of opposite sign.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
*          On entry, the upper quasi-triangular matrix T, in Schur
*          Schur canonical form.
*          On exit, the reordered upper quasi-triangular matrix, again
*          in Schur canonical form.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  IFST    (input/output) INTEGER
*  ILST    (input/output) INTEGER
*          Specify the reordering of the diagonal blocks of T.
*          The block with row index IFST is moved to row ILST, by a
*          sequence of transpositions between adjacent blocks.
*          On exit, if IFST pointed on entry to the second row of a
*          2-by-2 block, it is changed to point to the first row; ILST
*          always points to the first row of the block in its final
*          position (which may differ from its input value by +1 or -1).
*          1 <= IFST <= N; 1 <= ILST <= N.
*
*  NITRAF  (input/output) INTEGER
*          On entry, length of the array ITRAF.
*          As a minimum requirement, NITRAF >= max(1,|ILST-IFST|).
*          If there are 2-by-2 blocks in T then NITRAF must be larger;
*          a safe choice is NITRAF >= max(1,2*|ILST-IFST|).
*          On exit, actual length of the array ITRAF.
*
*  ITRAF   (output) INTEGER array, length NITRAF
*          List of parameters for representing the transformation
*          matrix Z, see further details.
*
*  NDTRAF  (input/output) INTEGER
*          On entry, length of the array DTRAF.
*          As a minimum requirement, NDTRAF >= max(1,2*|ILST-IFST|).
*          If there are 2-by-2 blocks in T then NDTRAF must be larger;
*          a safe choice is NDTRAF >= max(1,5*|ILST-IFST|).
*          On exit, actual length of the array DTRAF.
*
*  DTRAF   (output) DOUBLE PRECISION array, length NDTRAF
*          List of parameters for representing the transformation
*          matrix Z, see further details.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          = 1:  two adjacent blocks were too close to swap (the problem
*                is very ill-conditioned); T may have been partially
*                reordered, and ILST points to the first row of the
*                current position of the block being moved.
*          = 2:  the 2 by 2 block to be reordered split into two 1 by 1
*                blocks and the second block failed to swap with an
*                adjacent block. ILST points to the first row of the
*                current position of the whole block being moved.
*
*  Further Details
*  ===============
*
*  The orthogonal transformation matrix Z is a product of NITRAF
*  elementary orthogonal transformations. The parameters defining these
*  transformations are stored in the arrays ITRAF and DTRAF as follows:
*
*  Consider the i-th transformation acting on rows/columns POS,
*  POS+1, ... If this transformation is
*
*     (1) a Givens rotation with cosine C and sine S then
*
*           ITRAF(i) = POS,
*           DTRAF(j) = C,    DTRAF(j+1) = S;
*
*     (2) a Householder reflector H = I - tau * v * v' with
*         v = [ 1; v2; v3 ] then
*
*           ITRAF(i) = N + POS,
*           DTRAF(j) = tau,  DTRAF(j+1) = v2,  DTRAF(j+2) = v3;
*
*     (3) a Householder reflector H = I - tau * v * v' with
*         v = [ v1; v2; 1 ] then
*
*           ITRAF(i) = 2*N + POS,
*           DTRAF(j) = v1,  DTRAF(j+1) = v2,  DTRAF(j+2) = tau;
*
*  Note that the parameters in DTRAF are stored consecutively.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      INTEGER            DLNGTH(3), ILNGTH(3)
*     ..
*     .. Local Scalars ..
      INTEGER            CDTRAF, CITRAF, LDTRAF, LITRAF, HERE, I, NBF,
     $                   NBL, NBNEXT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           BDLAEXC, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     .. Data Statements ..
c      DATA ( ILNGTH(I), I = 1, 3 ) / 1, 2, 4 /
c      DATA ( DLNGTH(I), I = 1, 3 ) / 2, 5, 10 /
      DATA ILNGTH(1)/1/, ILNGTH(2)/2/, ILNGTH(3)/4/
      DATA DLNGTH(1)/2/, DLNGTH(2)/5/, DLNGTH(3)/10/
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input arguments.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( IFST.LT.1 .OR. IFST.GT.N ) THEN
         INFO = -4
      ELSE IF( ILST.LT.1 .OR. ILST.GT.N ) THEN
         INFO = -5
      ELSE IF ( NITRAF.LT.MAX( 1, ABS( ILST-IFST ) ) ) THEN
         INFO = -6
      ELSE IF ( NDTRAF.LT.MAX( 1, 2*ABS( ILST-IFST ) ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTREXC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
      CITRAF = 1
      CDTRAF = 1
*
*     Determine the first row of specified block
*     and find out it is 1 by 1 or 2 by 2.
*
      IF( IFST.GT.1 ) THEN
         IF( T( IFST, IFST-1 ).NE.ZERO )
     $      IFST = IFST - 1
      END IF
      NBF = 1
      IF( IFST.LT.N ) THEN
         IF( T( IFST+1, IFST ).NE.ZERO )
     $      NBF = 2
      END IF
*
*     Determine the first row of the final block
*     and find out it is 1 by 1 or 2 by 2.
*
      IF( ILST.GT.1 ) THEN
         IF( T( ILST, ILST-1 ).NE.ZERO )
     $      ILST = ILST - 1
      END IF
      NBL = 1
      IF( ILST.LT.N ) THEN
         IF( T( ILST+1, ILST ).NE.ZERO )
     $      NBL = 2
      END IF
*
      IF( IFST.EQ.ILST )
     $   RETURN
*
      IF( IFST.LT.ILST ) THEN
*
*        Update ILST
*
         IF( NBF.EQ.2 .AND. NBL.EQ.1 )
     $      ILST = ILST - 1
         IF( NBF.EQ.1 .AND. NBL.EQ.2 )
     $      ILST = ILST + 1
*
         HERE = IFST
*
   10    CONTINUE
*
*        Swap block with next one below
*
         IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
*
*           Current block either 1 by 1 or 2 by 2
*
            NBNEXT = 1
            IF( HERE+NBF+1.LE.N ) THEN
               IF( T( HERE+NBF+1, HERE+NBF ).NE.ZERO )
     $            NBNEXT = 2
            END IF
*
            LITRAF = ILNGTH(NBF+NBNEXT-1)
            LDTRAF = DLNGTH(NBF+NBNEXT-1)
            IF( CITRAF+LITRAF-1.GT.NITRAF ) THEN
               INFO = -6
               CALL XERBLA( 'BDTREXC', -INFO )
               RETURN
            END IF
            IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
               INFO = -8
               CALL XERBLA( 'BDTREXC', -INFO )
               RETURN
            END IF
            CALL BDLAEXC( N, T, LDT, HERE, NBF, NBNEXT, ITRAF(CITRAF),
     $                    DTRAF(CDTRAF), WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               NITRAF = CITRAF - 1
               NDTRAF = CDTRAF - 1
               RETURN
            END IF
            CITRAF = CITRAF + LITRAF
            CDTRAF = CDTRAF + LDTRAF
            HERE = HERE + NBNEXT
*
*           Test if 2 by 2 block breaks into two 1 by 1 blocks
*
            IF( NBF.EQ.2 ) THEN
               IF( T( HERE+1, HERE ).EQ.ZERO )
     $            NBF = 3
            END IF
*
         ELSE
*
*           Current block consists of two 1 by 1 blocks each of which
*           must be swapped individually
*
            NBNEXT = 1
            IF( HERE+3.LE.N ) THEN
               IF( T( HERE+3, HERE+2 ).NE.ZERO )
     $            NBNEXT = 2
            END IF
            LITRAF = ILNGTH(NBNEXT)
            LDTRAF = DLNGTH(NBNEXT)
            IF( CITRAF+LITRAF-1.GT.NITRAF ) THEN
               INFO = -6
               CALL XERBLA( 'BDTREXC', -INFO )
               RETURN
            END IF
            IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
               INFO = -8
               CALL XERBLA( 'BDTREXC', -INFO )
               RETURN
            END IF
            CALL BDLAEXC( N, T, LDT, HERE+1, 1, NBNEXT, ITRAF(CITRAF),
     $                    DTRAF(CDTRAF), WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               NITRAF = CITRAF - 1
               NDTRAF = CDTRAF - 1
               RETURN
            END IF
            CITRAF = CITRAF + LITRAF
            CDTRAF = CDTRAF + LDTRAF
*
            IF( NBNEXT.EQ.1 ) THEN
*
*              Swap two 1 by 1 blocks, no problems possible
*
               LITRAF = ILNGTH(1)
               LDTRAF = DLNGTH(1)
               IF( CITRAF+LITRAF-1.GT.NITRAF ) THEN
                  INFO = -6
                  CALL XERBLA( 'BDTREXC', -INFO )
                  RETURN
               END IF
               IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
                  INFO = -8
                  CALL XERBLA( 'BDTREXC', -INFO )
                  RETURN
               END IF
               CALL BDLAEXC( N, T, LDT, HERE, 1, NBNEXT, ITRAF(CITRAF),
     $                       DTRAF(CDTRAF), WORK, INFO )
               CITRAF = CITRAF + LITRAF
               CDTRAF = CDTRAF + LDTRAF
               HERE = HERE + 1
            ELSE
*
*              Recompute NBNEXT in case 2 by 2 split
*
               IF( T( HERE+2, HERE+1 ).EQ.ZERO )
     $            NBNEXT = 1
               IF( NBNEXT.EQ.2 ) THEN
*
*                 2 by 2 Block did not split
*
                  LITRAF = ILNGTH(2)
                  LDTRAF = DLNGTH(2)
                  IF( CITRAF+LITRAF-1.GT.NITRAF ) THEN
                     INFO = -6
                     CALL XERBLA( 'BDTREXC', -INFO )
                     RETURN
                  END IF
                  IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
                     INFO = -8
                     CALL XERBLA( 'BDTREXC', -INFO )
                     RETURN
                  END IF
                  CALL BDLAEXC( N, T, LDT, HERE, 1, NBNEXT,
     $                          ITRAF(CITRAF), DTRAF(CDTRAF), WORK,
     $                          INFO )
                  IF( INFO.NE.0 ) THEN
                     INFO = 2
                     ILST = HERE
                     NITRAF = CITRAF - 1
                     NDTRAF = CDTRAF - 1
                     RETURN
                  END IF
                  CITRAF = CITRAF + LITRAF
                  CDTRAF = CDTRAF + LDTRAF
                  HERE = HERE + 2
               ELSE
*
*                 2 by 2 Block did split
*
                  LITRAF = ILNGTH(1)
                  LDTRAF = DLNGTH(1)
                  IF( CITRAF+2*LITRAF-1.GT.NITRAF ) THEN
                     INFO = -6
                     CALL XERBLA( 'BDTREXC', -INFO )
                     RETURN
                  END IF
                  IF( CDTRAF+2*LDTRAF-1.GT.NDTRAF ) THEN
                     INFO = -8
                     CALL XERBLA( 'BDTREXC', -INFO )
                     RETURN
                  END IF
                  CALL BDLAEXC( N, T, LDT, HERE, 1, 1, ITRAF(CITRAF),
     $                          DTRAF(CDTRAF), WORK, INFO )
                  CITRAF = CITRAF + LITRAF
                  CDTRAF = CDTRAF + LDTRAF
                  CALL BDLAEXC( N, T, LDT, HERE+1, 1, 1, ITRAF(CITRAF),
     $                          DTRAF(CDTRAF), WORK, INFO )
                  CITRAF = CITRAF + LITRAF
                  CDTRAF = CDTRAF + LDTRAF
                  HERE = HERE + 2
               END IF
            END IF
         END IF
         IF( HERE.LT.ILST )
     $      GO TO 10
*
      ELSE
*
         HERE = IFST
   20    CONTINUE
*
*        Swap block with next one above
*
         IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
*
*           Current block either 1 by 1 or 2 by 2
*
            NBNEXT = 1
            IF( HERE.GE.3 ) THEN
               IF( T( HERE-1, HERE-2 ).NE.ZERO )
     $            NBNEXT = 2
            END IF
*
            LITRAF = ILNGTH(NBF+NBNEXT-1)
            LDTRAF = DLNGTH(NBF+NBNEXT-1)
            IF( CITRAF+LITRAF-1.GT.NITRAF ) THEN
               INFO = -6
               CALL XERBLA( 'BDTREXC', -INFO )
               RETURN
            END IF
            IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
               INFO = -8
               CALL XERBLA( 'BDTREXC', -INFO )
               RETURN
            END IF
            CALL BDLAEXC( N, T, LDT, HERE-NBNEXT, NBNEXT, NBF,
     $                    ITRAF(CITRAF), DTRAF(CDTRAF), WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               NITRAF = CITRAF - 1
               NDTRAF = CDTRAF - 1
               RETURN
            END IF
            CITRAF = CITRAF + LITRAF
            CDTRAF = CDTRAF + LDTRAF
            HERE = HERE - NBNEXT
*
*           Test if 2 by 2 block breaks into two 1 by 1 blocks
*
            IF( NBF.EQ.2 ) THEN
               IF( T( HERE+1, HERE ).EQ.ZERO )
     $            NBF = 3
            END IF
*
         ELSE
*
*           Current block consists of two 1 by 1 blocks each of which
*           must be swapped individually
*
            NBNEXT = 1
            IF( HERE.GE.3 ) THEN
               IF( T( HERE-1, HERE-2 ).NE.ZERO )
     $            NBNEXT = 2
            END IF
            LITRAF = ILNGTH(NBNEXT)
            LDTRAF = DLNGTH(NBNEXT)
            IF( CITRAF+LITRAF-1.GT.NITRAF ) THEN
               INFO = -6
               CALL XERBLA( 'BDTREXC', -INFO )
               RETURN
            END IF
            IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
               INFO = -8
               CALL XERBLA( 'BDTREXC', -INFO )
               RETURN
            END IF
            CALL BDLAEXC( N, T, LDT, HERE-NBNEXT, NBNEXT, 1,
     $                    ITRAF(CITRAF), DTRAF(CDTRAF), WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               NITRAF = CITRAF - 1
               NDTRAF = CDTRAF - 1
               RETURN
            END IF
            CITRAF = CITRAF + LITRAF
            CDTRAF = CDTRAF + LDTRAF
*
            IF( NBNEXT.EQ.1 ) THEN
*
*              Swap two 1 by 1 blocks, no problems possible
*
               LITRAF = ILNGTH(1)
               LDTRAF = DLNGTH(1)
               IF( CITRAF+LITRAF-1.GT.NITRAF ) THEN
                  INFO = -6
                  CALL XERBLA( 'BDTREXC', -INFO )
                  RETURN
               END IF
               IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
                  INFO = -8
                  CALL XERBLA( 'BDTREXC', -INFO )
                  RETURN
               END IF
               CALL BDLAEXC( N, T, LDT, HERE, NBNEXT, 1, ITRAF(CITRAF),
     $                       DTRAF(CDTRAF), WORK, INFO )
               CITRAF = CITRAF + LITRAF
               CDTRAF = CDTRAF + LDTRAF
               HERE = HERE - 1
            ELSE
*
*              Recompute NBNEXT in case 2 by 2 split
*
               IF( T( HERE, HERE-1 ).EQ.ZERO )
     $            NBNEXT = 1
               IF( NBNEXT.EQ.2 ) THEN
*
*                 2 by 2 Block did not split
*
                  LITRAF = ILNGTH(2)
                  LDTRAF = DLNGTH(2)
                  IF( CITRAF+LITRAF-1.GT.NITRAF ) THEN
                     INFO = -6
                     CALL XERBLA( 'BDTREXC', -INFO )
                     RETURN
                  END IF
                  IF( CDTRAF+LDTRAF-1.GT.NDTRAF ) THEN
                     INFO = -8
                     CALL XERBLA( 'BDTREXC', -INFO )
                     RETURN
                  END IF
                  CALL BDLAEXC( N, T, LDT, HERE-1, 2, 1, ITRAF(CITRAF),
     $                          DTRAF(CDTRAF), WORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     INFO = 2
                     ILST = HERE
                     NITRAF = CITRAF - 1
                     NDTRAF = CDTRAF - 1
                     RETURN
                  END IF
                  CITRAF = CITRAF + LITRAF
                  CDTRAF = CDTRAF + LDTRAF
                  HERE = HERE - 2
               ELSE
*
*                 2 by 2 Block did split
*
                  LITRAF = ILNGTH(1)
                  LDTRAF = DLNGTH(1)
                  IF( CITRAF+2*LITRAF-1.GT.NITRAF ) THEN
                     INFO = -6
                     CALL XERBLA( 'BDTREXC', -INFO )
                     RETURN
                  END IF
                  IF( CDTRAF+2*LDTRAF-1.GT.NDTRAF ) THEN
                     INFO = -8
                     CALL XERBLA( 'BDTREXC', -INFO )
                     RETURN
                  END IF
                  CALL BDLAEXC( N, T, LDT, HERE, 1, 1, ITRAF(CITRAF),
     $                          DTRAF(CDTRAF), WORK, INFO )
                  CITRAF = CITRAF + LITRAF
                  CDTRAF = CDTRAF + LDTRAF
                  CALL BDLAEXC( N, T, LDT, HERE-1, 1, 1, ITRAF(CITRAF),
     $                          DTRAF(CDTRAF), WORK, INFO )
                  CITRAF = CITRAF + LITRAF
                  CDTRAF = CDTRAF + LDTRAF
                  HERE = HERE - 2
               END IF
            END IF
         END IF
         IF( HERE.GT.ILST )
     $      GO TO 20
      END IF
      ILST = HERE
      NITRAF = CITRAF - 1
      NDTRAF = CDTRAF - 1
*
      RETURN
*
*     End of BDTREXC
*
      END
