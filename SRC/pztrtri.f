      SUBROUTINE PZTRTRI( UPLO, DIAG, N, A, IA, JA, DESCA, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            IA, INFO, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         A( * )
*     ..
*
*  Purpose
*  =======
*
*  PZTRTRI computes the inverse of a upper or lower triangular
*  distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1).
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
*  UPLO    (global input) CHARACTER
*          Specifies whether the distributed matrix sub( A ) is upper
*          or lower triangular:
*          = 'U':  Upper triangular,
*          = 'L':  Lower triangular.
*
*  DIAG    (global input) CHARACTER
*          Specifies whether or not the distributed matrix sub( A )
*          is unit triangular:
*          = 'N':  Non-unit triangular,
*          = 'U':  Unit triangular.
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          triangular matrix sub( A ).  If UPLO = 'U', the leading
*          N-by-N upper triangular part of the matrix sub( A ) contains
*          the upper triangular matrix to be inverted, and the strictly
*          lower triangular part of sub( A ) is not referenced.
*          If UPLO = 'L', the leading N-by-N lower triangular part of
*          the matrix sub( A ) contains the lower triangular matrix,
*          and the strictly upper triangular part of sub( A ) is not
*          referenced.
*          On exit, the (triangular) inverse of the original matrix.
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
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*            > 0:  If INFO = K, A(IA+K-1,JA+K-1) is exactly zero.  The
*                  triangular matrix sub( A ) is singular and its
*                  inverse can not be computed.
*
*  ====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            I, ICOFF, ICTXT, IROFF, ICURCOL, ICURROW,
     $                   IDUMMY, II, IOFFA, J, JB, JJ, JN, LDA, MYCOL,
     $                   MYROW, NN, NPCOL, NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 2 ), IDUM2( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, IGAMX2D, INFOG2L,
     $                   PCHK1MAT, PXERBLA, PZTRTI2, PZTRMM,
     $                   PZTRSM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL
      EXTERNAL           ICEIL, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(700+CTXT_)
      ELSE
         UPPER = LSAME( UPLO, 'U' )
         NOUNIT = LSAME( DIAG, 'N' )
*
         CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, INFO )
         IF( INFO.EQ.0 ) THEN
            IROFF = MOD( IA-1, DESCA( MB_ ) )
            ICOFF = MOD( JA-1, DESCA( NB_ ) )
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
               INFO = -2
            ELSE IF( IROFF.NE.ICOFF .OR. IROFF.NE.0 ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(700+NB_)
            END IF
         END IF
*
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 1
         IF( NOUNIT ) THEN
            IDUM1( 2 ) = ICHAR( 'N' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'U' )
         END IF
         IDUM2( 2 ) = 2
*
         CALL PCHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, 2, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZTRTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      IF( NOUNIT ) THEN
         CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 II, JJ, ICURROW, ICURCOL )
*
*        Handle first block separately
*
         JB = JN-JA+1
         LDA = DESCA( LLD_ )
         IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
            IOFFA = II+(JJ-1)*LDA
            DO 10 I = 0, JB-1
               IF( A( IOFFA ).EQ.ZERO .AND. INFO.EQ.0 )
     $            INFO = I + 1
               IOFFA = IOFFA + LDA + 1
   10       CONTINUE
         END IF
         IF( MYROW.EQ.ICURROW )
     $      II = II + JB
         IF( MYCOL.EQ.ICURCOL )
     $      JJ = JJ + JB
         ICURROW = MOD( ICURROW+1, NPROW )
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
*        Loop over remaining blocks of columns
*
         DO 30 J = JN+1, JA+N-1, DESCA( NB_ )
            JB = MIN( JA+N-J, DESCA( NB_ ) )
            IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
               IOFFA = II+(JJ-1)*LDA
               DO 20 I = 0, JB-1
                  IF( A( IOFFA ).EQ.ZERO .AND. INFO.EQ.0 )
     $               INFO = J + I - JA + 1
                  IOFFA = IOFFA + LDA + 1
   20          CONTINUE
            END IF
            IF( MYROW.EQ.ICURROW )
     $         II = II + JB
            IF( MYCOL.EQ.ICURCOL )
     $         JJ = JJ + JB
            ICURROW = MOD( ICURROW+1, NPROW )
            ICURCOL = MOD( ICURCOL+1, NPCOL )
   30    CONTINUE
         CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, IDUMMY,
     $                 IDUMMY, -1, -1, MYCOL )
         IF( INFO.NE.0 )
     $      RETURN
      END IF
*
*     Use blocked code
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix
*
         JB = JN-JA+1
*
*        Handle first block of column separately
*
         CALL PZTRTI2( UPLO, DIAG, JB, A, IA, JA, DESCA, INFO )
*
*        Loop over remaining block of columns
*
         DO 40 J = JN+1, JA+N-1, DESCA( NB_ )
            JB = MIN( DESCA( NB_ ), JA+N-J )
            I = IA + J - JA
*
*           Compute rows 1:j-1 of current block column
*
            CALL PZTRMM( 'Left', UPLO, 'No transpose', DIAG, J-JA, JB,
     $                   ONE, A, IA, JA, DESCA, A, IA, J, DESCA )
            CALL PZTRSM( 'Right', UPLO, 'No transpose', DIAG, J-JA,
     $                   JB, -ONE, A, I, J, DESCA, A, IA, J, DESCA )
*
*           Compute inverse of current diagonal block
*
            CALL PZTRTI2( UPLO, DIAG, JB, A, I, J, DESCA, INFO )
*
   40    CONTINUE
*
      ELSE
*
*        Compute inverse of lower triangular matrix
*
         NN = ( ( JA+N-2 ) / DESCA( NB_ ) )*DESCA( NB_ ) + 1
         DO 50 J = NN, JN+1, -DESCA( NB_ )
            JB = MIN( DESCA( NB_ ), JA+N-J )
            I = IA + J - JA
            IF( J+JB.LE.JA+N-1 ) THEN
*
*              Compute rows j+jb:ja+n-1 of current block column
*
               CALL PZTRMM( 'Left', UPLO, 'No transpose', DIAG,
     $                      JA+N-J-JB, JB, ONE, A, I+JB, J+JB, DESCA,
     $                      A, I+JB, J, DESCA )
               CALL PZTRSM( 'Right', UPLO, 'No transpose', DIAG,
     $                      JA+N-J-JB, JB, -ONE, A, I, J, DESCA,
     $                      A, I+JB, J, DESCA )
            END IF
*
*           Compute inverse of current diagonal block
*
            CALL PZTRTI2( UPLO, DIAG, JB, A, I, J, DESCA, INFO )
*
   50    CONTINUE
*
*        Handle the last block of columns separately
*
         JB = JN-JA+1
         IF( JA+JB.LE.JA+N-1 ) THEN
*
*           Compute rows ja+jb:ja+n-1 of current block column
*
            CALL PZTRMM( 'Left', UPLO, 'No transpose', DIAG, N-JB, JB,
     $                   ONE, A, IA+JB, JA+JB, DESCA, A, IA+JB, JA,
     $                   DESCA )
            CALL PZTRSM( 'Right', UPLO, 'No transpose', DIAG, N-JB, JB,
     $                   -ONE, A, IA, JA, DESCA, A, IA+JB, JA, DESCA )
         END IF
*
*        Compute inverse of current diagonal block
*
         CALL PZTRTI2( UPLO, DIAG, JB, A, IA, JA, DESCA, INFO )
*
      END IF
*
      RETURN
*
*     End PZTRTRI
*
      END
