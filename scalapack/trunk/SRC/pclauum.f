      SUBROUTINE PCLAUUM( UPLO, N, A, IA, JA, DESCA )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX            A( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLAUUM computes the product U * U' or L' * L, where the triangular
*  factor U or L is stored in the upper or lower triangular part of
*  the distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1).
*
*  If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
*  overwriting the factor U in sub( A ).
*  If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
*  overwriting the factor L in sub( A ).
*
*  This is the blocked form of the algorithm, calling Level 3 PBLAS.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the triangular factor stored in the
*          distributed matrix sub( A ) is upper or lower triangular:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the triangular factor U or L. N >= 0.
*
*  A       (local input/local output) COMPLEX pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, the local pieces of the triangular factor L or U.
*          On exit, if UPLO = 'U', the upper triangle of the distributed
*          matrix sub( A ) is overwritten with the upper triangle of the
*          product U * U'; if UPLO = 'L', the lower triangle of sub( A )
*          is overwritten with the lower triangle of the product L' * L.
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JB, JN
*     ..
*     .. External Subroutines ..
      EXTERNAL           PCGEMM, PCHERK, PCLAUU2, PCTRMM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL
      EXTERNAL           ICEIL, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      IF(  LSAME( UPLO, 'U' ) ) THEN
*
*        Compute the product U * U'.
*
*        Handle first block separately
*
         JB = JN-JA+1
         CALL PCLAUU2( 'Upper', JB, A, IA, JA, DESCA )
         IF( JB.LE.N-1 ) THEN
            CALL PCHERK( 'Upper', 'No transpose', JB, N-JB, ONE, A, IA,
     $                   JA+JB, DESCA, ONE, A, IA, JA, DESCA )
         END IF
*
*        Loop over remaining block of columns
*
         DO 10 J = JN+1, JA+N-1, DESCA( NB_ )
            JB = MIN( N-J+JA, DESCA( NB_ ) )
            I = IA + J - JA
            CALL PCTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                   'Non-unit', J-JA, JB, CONE, A, I, J, DESCA,
     $                   A, IA, J, DESCA )
            CALL PCLAUU2( 'Upper', JB, A, I, J, DESCA )
            IF( J+JB.LE.JA+N-1 ) THEN
               CALL PCGEMM( 'No transpose', 'Conjugate transpose',
     $                      J-JA, JB, N-J-JB+JA, CONE, A, IA, J+JB,
     $                      DESCA, A, I, J+JB, DESCA, CONE, A, IA,
     $                      J, DESCA )
               CALL PCHERK( 'Upper', 'No transpose', JB, N-J-JB+JA, ONE,
     $                      A, I, J+JB, DESCA, ONE, A, I, J, DESCA )
            END IF
   10    CONTINUE
      ELSE
*
*        Compute the product L' * L.
*
*        Handle first block separately
*
         JB = JN-JA+1
         CALL PCLAUU2( 'Lower', JB, A, IA, JA, DESCA )
         IF( JB.LE.N-1 ) THEN
            CALL PCHERK( 'Lower', 'Conjugate transpose', JB, N-JB, ONE,
     $                   A, IA+JB, JA, DESCA, ONE, A, IA, JA, DESCA )
         END IF
*
*        Loop over remaining block of columns
*
         DO 20 J = JN+1, JA+N-1, DESCA( NB_ )
            JB = MIN( N-J+JA, DESCA( NB_ ) )
            I = IA + J - JA
            CALL PCTRMM( 'Left', 'Lower', 'Conjugate Transpose',
     $                   'Non-unit', JB, J-JA, CONE, A, I, J, DESCA, A,
     $                   I, JA, DESCA )
            CALL PCLAUU2( 'Lower', JB, A, I, J, DESCA )
            IF( J+JB.LE.JA+N-1 ) THEN
               CALL PCGEMM( 'Conjugate transpose', 'No transpose', JB,
     $                      J-JA, N-J-JB+JA, CONE, A, I+JB, J, DESCA,
     $                      A, I+JB, JA, DESCA, CONE, A, I, JA, DESCA )
               CALL PCHERK( 'Lower', 'Conjugate transpose', JB,
     $                      N-J-JB+JA, ONE, A, I+JB, J, DESCA, ONE,
     $                      A, I, J, DESCA )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of PCLAUUM
*
      END
