      SUBROUTINE PZPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
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
*  PZPOTRF computes the Cholesky factorization of an N-by-N complex
*  hermitian positive definite distributed matrix sub( A ) denoting
*  A(IA:IA+N-1, JA:JA+N-1).
*
*  The factorization has the form
*
*            sub( A ) = U' * U ,  if UPLO = 'U', or
*
*            sub( A ) = L  * L',  if UPLO = 'L',
*
*  where U is an upper triangular matrix and L is lower triangular.
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
*  This routine requires square block decomposition ( MB_A = NB_A ).
*
*  Arguments
*  =========
*
*  UPLO    (global input) CHARACTER
*          = 'U':  Upper triangle of sub( A ) is stored;
*          = 'L':  Lower triangle of sub( A ) is stored.
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          N-by-N Hermitian distributed matrix sub( A ) to be factored.
*          If UPLO = 'U', the leading N-by-N upper triangular part of
*          sub( A ) contains the upper triangular part of the matrix,
*          and its strictly lower triangular part is not referenced.
*          If UPLO = 'L', the leading N-by-N lower triangular part of
*          sub( A ) contains the lower triangular part of the distribu-
*          ted matrix, and its strictly upper triangular part is not
*          referenced. On exit, if UPLO = 'U', the upper triangular
*          part of the distributed matrix contains the Cholesky factor
*          U, if UPLO = 'L', the lower triangular part of the distribu-
*          ted matrix contains the Cholesky factor L.
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
*          > 0:  If INFO = K, the leading minor of order K,
*                A(IA:IA+K-1,JA:JA+K-1) is not positive definite, and
*                the factorization could not be completed.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            I, ICOFF, ICTXT, IROFF, J, JB, JN, MYCOL,
     $                   MYROW, NPCOL, NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK1MAT, PB_TOPGET,
     $                   PB_TOPSET, PXERBLA, PZPOTF2, PZHERK,
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
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(600+CTXT_)
      ELSE
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
         UPPER = LSAME( UPLO, 'U' )
         IF( INFO.EQ.0 ) THEN
            IROFF = MOD( IA-1, DESCA( MB_ ) )
            ICOFF = MOD( JA-1, DESCA( NB_ ) )
            IF ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( IROFF.NE.0 ) THEN
               INFO = -4
            ELSE IF( ICOFF.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(600+NB_)
            END IF
         END IF
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 1
         CALL PCHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, 1, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZPOTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      IF( UPPER ) THEN
*
*        Split-ring topology for the communication along process
*        columns, 1-tree topology along process rows.
*
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ' ' )
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', 'S-ring' )
*
*        A is upper triangular, compute Cholesky factorization A = U'*U.
*
*        Handle the first block of columns separately
*
         JN = MIN( ICEIL( JA, DESCA( NB_ ) )*DESCA(NB_), JA+N-1 )
         JB = JN - JA + 1
*
*        Perform unblocked Cholesky factorization on JB block
*
         CALL PZPOTF2( UPLO, JB, A, IA, JA, DESCA, INFO )
         IF( INFO.NE.0 )
     $      GO TO 30
*
         IF( JB+1.LE.N ) THEN
*
*           Form the row panel of U using the triangular solver
*
            CALL PZTRSM( 'Left', UPLO, 'Conjugate transpose',
     $                   'Non-Unit', JB, N-JB, CONE, A, IA, JA, DESCA,
     $                   A, IA, JA+JB, DESCA )
*
*           Update the trailing matrix, A = A - U'*U
*
            CALL PZHERK( UPLO, 'Conjugate transpose', N-JB, JB, -ONE, A,
     $                   IA, JA+JB, DESCA, ONE, A, IA+JB, JA+JB, DESCA )
         END IF
*
*        Loop over remaining block of columns
*
         DO 10 J = JN+1, JA+N-1, DESCA( NB_ )
            JB = MIN( N-J+JA, DESCA( NB_ ) )
            I = IA + J - JA
*
*           Perform unblocked Cholesky factorization on JB block
*
            CALL PZPOTF2( UPLO, JB, A, I, J, DESCA, INFO )
            IF( INFO.NE.0 ) THEN
               INFO = INFO + J - JA
               GO TO 30
            END IF
*
            IF( J-JA+JB+1.LE.N ) THEN
*
*              Form the row panel of U using the triangular solver
*
               CALL PZTRSM( 'Left', UPLO, 'Conjugate transpose',
     $                      'Non-Unit', JB, N-J-JB+JA, CONE, A, I, J,
     $                      DESCA, A, I, J+JB, DESCA )
*
*              Update the trailing matrix, A = A - U'*U
*
               CALL PZHERK( UPLO, 'Conjugate transpose', N-J-JB+JA, JB,
     $                      -ONE, A, I, J+JB, DESCA, ONE, A, I+JB,
     $                      J+JB, DESCA )
            END IF
   10    CONTINUE
*
      ELSE
*
*        1-tree topology for the communication along process columns,
*        Split-ring topology along process rows.
*
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'S-ring' )
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
*
*        A is lower triangular, compute Cholesky factorization A = L*L'
*        (right-looking)
*
*        Handle the first block of columns separately
*
         JN = MIN( ICEIL( JA, DESCA( NB_ ) )*DESCA( NB_ ), JA+N-1 )
         JB = JN - JA + 1
*
*        Perform unblocked Cholesky factorization on JB block
*
         CALL PZPOTF2( UPLO, JB, A, IA, JA, DESCA, INFO )
         IF( INFO.NE.0 )
     $      GO TO 30
*
         IF( JB+1.LE.N ) THEN
*
*           Form the column panel of L using the triangular solver
*
            CALL PZTRSM( 'Right', UPLO, 'Conjugate transpose',
     $                   'Non-Unit', N-JB, JB, CONE, A, IA, JA, DESCA,
     $                   A, IA+JB, JA, DESCA )
*
*           Update the trailing matrix, A = A - L*L'
*
            CALL PZHERK( UPLO, 'No Transpose', N-JB, JB, -ONE, A, IA+JB,
     $                   JA, DESCA, ONE, A, IA+JB, JA+JB, DESCA )
*
         END IF
*
         DO 20 J = JN+1, JA+N-1, DESCA( NB_ )
            JB = MIN( N-J+JA, DESCA( NB_ ) )
            I = IA + J - JA
*
*           Perform unblocked Cholesky factorization on JB block
*
            CALL PZPOTF2( UPLO, JB, A, I, J, DESCA, INFO )
            IF( INFO.NE.0 ) THEN
               INFO = INFO + J - JA
               GO TO 30
            END IF
*
            IF( J-JA+JB+1.LE.N ) THEN
*
*              Form the column panel of L using the triangular solver
*
               CALL PZTRSM( 'Right', UPLO, 'Conjugate transpose',
     $                      'Non-Unit', N-J-JB+JA, JB, CONE, A, I, J,
     $                      DESCA, A, I+JB, J, DESCA )
*
*              Update the trailing matrix, A = A - L*L'
*
               CALL PZHERK( UPLO, 'No Transpose', N-J-JB+JA, JB, -ONE,
     $                      A, I+JB, J, DESCA, ONE, A, I+JB, J+JB,
     $                      DESCA )
*
            END IF
   20    CONTINUE
*
      END IF
*
   30 CONTINUE
*
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      RETURN
*
*     End of PZPOTRF
*
      END
