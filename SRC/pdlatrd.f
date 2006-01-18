      SUBROUTINE PDLATRD( UPLO, N, NB, A, IA, JA, DESCA, D, E, TAU, W,
     $                    IW, JW, DESCW, WORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, IW, JA, JW, N, NB
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCW( * )
      DOUBLE PRECISION   A( * ), D( * ), E( * ), TAU( * ), W( * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLATRD reduces NB rows and columns of a real symmetric distributed
*  matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) to symmetric tridiagonal
*  form by an orthogonal similarity transformation Q' * sub( A ) * Q,
*  and returns the matrices V and W which are needed to apply the
*  transformation to the unreduced part of sub( A ).
*
*  If UPLO = 'U', PDLATRD reduces the last NB rows and columns of a
*  matrix, of which the upper triangle is supplied;
*  if UPLO = 'L', PDLATRD reduces the first NB rows and columns of a
*  matrix, of which the lower triangle is supplied.
*
*  This is an auxiliary routine called by PDSYTRD.
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
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix sub( A ) is stored:
*          = 'U': Upper triangular
*          = 'L': Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  NB      (global input) INTEGER
*          The number of rows and columns to be reduced.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          symmetric distributed matrix sub( A ).  If UPLO = 'U', the
*          leading N-by-N upper triangular part of sub( A ) contains
*          the upper triangular part of the matrix, and its strictly
*          lower triangular part is not referenced. If UPLO = 'L', the
*          leading N-by-N lower triangular part of sub( A ) contains the
*          lower triangular part of the matrix, and its strictly upper
*          triangular part is not referenced.
*          On exit, if UPLO = 'U', the last NB columns have been reduced
*          to tridiagonal form, with the diagonal elements overwriting
*          the diagonal elements of sub( A ); the elements above the
*          diagonal with the array TAU, represent the orthogonal matrix
*          Q as a product of elementary reflectors. If UPLO = 'L', the
*          first NB columns have been reduced to tridiagonal form, with
*          the diagonal elements overwriting the diagonal elements of
*          sub( A ); the elements below the diagonal with the array TAU,
*          represent the orthogonal matrix Q as a product of elementary
*          reflectors; See Further Details.
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
*  D       (local output) DOUBLE PRECISION array, dimension LOCc(JA+N-1)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i). D is tied to the distributed matrix A.
*
*  E       (local output) DOUBLE PRECISION array, dimension LOCc(JA+N-1)
*          if UPLO = 'U', LOCc(JA+N-2) otherwise. The off-diagonal
*          elements of the tridiagonal matrix T: E(i) = A(i,i+1) if
*          UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. E is tied to the
*          distributed matrix A.
*
*  TAU     (local output) DOUBLE PRECISION array, dimension
*          LOCc(JA+N-1). This array contains the scalar factors TAU of
*          the elementary reflectors. TAU is tied to the distributed
*          matrix A.
*
*  W       (local output) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_W,NB_W), This array contains
*          the local pieces of the N-by-NB_W matrix W required to
*          update the unreduced part of sub( A ).
*
*  IW      (global input) INTEGER
*          The row index in the global array W indicating the first
*          row of sub( W ).
*
*  JW      (global input) INTEGER
*          The column index in the global array W indicating the
*          first column of sub( W ).
*
*  DESCW   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix W.
*
*  WORK    (local workspace) DOUBLE PRECISION array, dimension (NB_A)
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n) H(n-1) . . . H(n-nb+1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in
*  A(ia:ia+i-2,ja+i), and tau in TAU(ja+i-1).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(nb).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in
*  A(ia+i+1:ia+n-1,ja+i-1), and tau in TAU(ja+i-1).
*
*  The elements of the vectors v together form the N-by-NB matrix V
*  which is needed, with W, to apply the transformation to the unreduced
*  part of the matrix, using a symmetric rank-2k update of the form:
*  sub( A ) := sub( A ) - V*W' - W*V'.
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5 and nb = 2:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  a   a   a   v4  v5 )              (  d                  )
*    (      a   a   v4  v5 )              (  1   d              )
*    (          a   1   v5 )              (  v1  1   a          )
*    (              d   1  )              (  v1  v2  a   a      )
*    (                  d  )              (  v1  v2  a   a   a  )
*
*  where d denotes a diagonal element of the reduced matrix, a denotes
*  an element of the original matrix that is unchanged, and vi denotes
*  an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   HALF, ONE, ZERO
      PARAMETER          ( HALF = 0.5D+0, ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IACOL, IAROW, ICTXT, II, J, JJ, JP, JWK, K,
     $                   KW, MYCOL, MYROW, NPCOL, NPROW, NQ
      DOUBLE PRECISION   ALPHA
*     ..
*     .. Local Arrays ..
      INTEGER            DESCD( DLEN_ ), DESCE( DLEN_ ), DESCWK( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, DGEBR2D, DGEBS2D,
     $                   INFOG2L, PDAXPY, PDDOT, PDELGET,
     $                   PDELSET, PDGEMV, PDLARFG, PDSCAL,
     $                   PDSYMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      EXTERNAL           LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NQ = MAX( 1, NUMROC( JA+N-1, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $          NPCOL ) )
      CALL DESCSET( DESCD, 1, JA+N-1, 1, DESCA( NB_ ), MYROW,
     $              DESCA( CSRC_ ), DESCA( CTXT_ ), 1 )
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
         CALL INFOG2L( N+IA-NB, N+JA-NB, DESCA, NPROW, NPCOL, MYROW,
     $                 MYCOL, II, JJ, IAROW, IACOL )
         CALL DESCSET( DESCWK, 1, DESCW( NB_ ), 1, DESCW( NB_ ), IAROW,
     $                 IACOL, ICTXT, 1 )
         CALL DESCSET( DESCE, 1, JA+N-1, 1, DESCA( NB_ ), MYROW,
     $                 DESCA( CSRC_ ), DESCA( CTXT_ ), 1 )
*
*        Reduce last NB columns of upper triangle
*
         DO 10 J = JA+N-1, JA+N-NB, -1
            I = IA + J - JA
            K = J - JA + 1
            KW = MOD( K-1, DESCA( MB_ ) ) + 1
*
*           Update A(IA:I,I)
*
            CALL PDGEMV( 'No transpose', K, N-K, -ONE, A, IA, J+1,
     $                   DESCA, W, IW+K-1, JW+KW, DESCW, DESCW( M_ ),
     $                   ONE, A, IA, J, DESCA, 1 )
            CALL PDGEMV( 'No transpose', K, N-K, -ONE, W, IW, JW+KW,
     $                   DESCW, A, I, J+1, DESCA, DESCA( M_ ), ONE, A,
     $                   IA, J, DESCA, 1 )
            IF( N-K.GT.0 )
     $         CALL PDELSET( A, I, J+1, DESCA, E( JP ) )
*
*           Generate elementary reflector H(i) to annihilate
*           A(IA:I-2,I)
*
            JP = MIN( JJ+KW-1, NQ )
            CALL PDLARFG( K-1, E( JP ), I-1, J, A, IA, J, DESCA, 1,
     $                    TAU )
            CALL PDELSET( A, I-1, J, DESCA, ONE )
*
*           Compute W(IW:IW+K-2,JW+KW-1)
*
            CALL PDSYMV( 'Upper', K-1, ONE, A, IA, JA, DESCA, A, IA, J,
     $                   DESCA, 1, ZERO, W, IW, JW+KW-1, DESCW, 1 )
*
            JWK = MOD( K-1, DESCWK( NB_ ) ) + 2
            CALL PDGEMV( 'Transpose', K-1, N-K, ONE, W, IW, JW+KW,
     $                   DESCW, A, IA, J, DESCA, 1, ZERO, WORK, 1, JWK,
     $                   DESCWK, DESCWK( M_ ) )
            CALL PDGEMV( 'No transpose', K-1, N-K, -ONE, A, IA, J+1,
     $                   DESCA, WORK, 1, JWK, DESCWK, DESCWK( M_ ), ONE,
     $                   W, IW, JW+KW-1, DESCW, 1 )
            CALL PDGEMV( 'Transpose', K-1, N-K, ONE, A, IA, J+1, DESCA,
     $                   A, IA, J, DESCA, 1, ZERO, WORK, 1, JWK, DESCWK,
     $                   DESCWK( M_ ) )
            CALL PDGEMV( 'No transpose', K-1, N-K, -ONE, W, IW, JW+KW,
     $                   DESCW, WORK, 1, JWK, DESCWK, DESCWK( M_ ), ONE,
     $                   W, IW, JW+KW-1, DESCW, 1 )
            CALL PDSCAL( K-1, TAU( JP ), W, IW, JW+KW-1, DESCW, 1 )
*
            CALL PDDOT( K-1, ALPHA, W, IW, JW+KW-1, DESCW, 1, A, IA, J,
     $                  DESCA, 1 )
            IF( MYCOL.EQ.IACOL )
     $         ALPHA = -HALF*TAU( JP )*ALPHA
            CALL PDAXPY( K-1, ALPHA, A, IA, J, DESCA, 1, W, IW, JW+KW-1,
     $                   DESCW, 1 )
            IF( MYCOL.EQ.IACOL ) THEN
               CALL PDELGET( 'E', ' ', D( JP ), A, I, J, DESCA )
            END IF
*
   10    CONTINUE
*
      ELSE
*
         CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, II,
     $                 JJ, IAROW, IACOL )
         CALL DESCSET( DESCWK, 1, DESCW( NB_ ), 1, DESCW( NB_ ), IAROW,
     $                 IACOL, ICTXT, 1 )
         CALL DESCSET( DESCE, 1, JA+N-2, 1, DESCA( NB_ ), MYROW,
     $                 DESCA( CSRC_ ), DESCA( CTXT_ ), 1 )
*
*        Reduce first NB columns of lower triangle
*
         DO 20 J = JA, JA+NB-1
            I = IA + J - JA
            K = J - JA + 1
*
*           Update A(J:JA+N-1,J)
*
            CALL PDGEMV( 'No transpose', N-K+1, K-1, -ONE, A, I, JA,
     $                   DESCA, W, IW+K-1, JW, DESCW, DESCW( M_ ), ONE,
     $                   A, I, J, DESCA, 1 )
            CALL PDGEMV( 'No transpose', N-K+1, K-1, -ONE, W, IW+K-1,
     $                   JW, DESCW, A, I, JA, DESCA, DESCA( M_ ), ONE,
     $                   A, I, J, DESCA, 1 )
            IF( K.GT.1 )
     $         CALL PDELSET( A, I, J-1, DESCA, E( JP ) )
*
*
*           Generate elementary reflector H(i) to annihilate
*           A(I+2:IA+N-1,I)
*
            JP = MIN( JJ+K-1, NQ )
            CALL PDLARFG( N-K, E( JP ), I+1, J, A, I+2, J, DESCA, 1,
     $                    TAU )
            CALL PDELSET( A, I+1, J, DESCA, ONE )
*
*           Compute W(IW+K:IW+N-1,JW+K-1)
*
            CALL PDSYMV( 'Lower', N-K, ONE, A, I+1, J+1, DESCA, A, I+1,
     $                   J, DESCA, 1, ZERO, W, IW+K, JW+K-1, DESCW, 1 )
*
            CALL PDGEMV( 'Transpose', N-K, K-1, ONE, W, IW+K, JW, DESCW,
     $                   A, I+1, J, DESCA, 1, ZERO, WORK, 1, 1, DESCWK,
     $                   DESCWK( M_ ) )
            CALL PDGEMV( 'No transpose', N-K, K-1, -ONE, A, I+1, JA,
     $                  DESCA, WORK, 1, 1, DESCWK, DESCWK( M_ ), ONE, W,
     $                  IW+K, JW+K-1, DESCW, 1 )
            CALL PDGEMV( 'Transpose', N-K, K-1, ONE, A, I+1, JA, DESCA,
     $                  A, I+1, J, DESCA, 1, ZERO, WORK, 1, 1, DESCWK,
     $                  DESCWK( M_ ) )
            CALL PDGEMV( 'No transpose', N-K, K-1, -ONE, W, IW+K, JW,
     $                  DESCW, WORK, 1, 1, DESCWK, DESCWK( M_ ), ONE, W,
     $                  IW+K, JW+K-1, DESCW, 1 )
            CALL PDSCAL( N-K, TAU( JP ), W, IW+K, JW+K-1, DESCW, 1 )
            CALL PDDOT( N-K, ALPHA, W, IW+K, JW+K-1, DESCW, 1, A, I+1,
     $                  J, DESCA, 1 )
            IF( MYCOL.EQ.IACOL )
     $         ALPHA = -HALF*TAU( JP )*ALPHA
            CALL PDAXPY( N-K, ALPHA, A, I+1, J, DESCA, 1, W, IW+K,
     $                   JW+K-1, DESCW, 1 )
            IF( MYCOL.EQ.IACOL ) THEN
               CALL PDELGET( 'E', ' ', D( JP ), A, I, J, DESCA )
            END IF
*
   20    CONTINUE
*
      END IF
*
*     Broadcast columnwise the diagonal elements into D.
*
      IF( MYCOL.EQ.IACOL ) THEN
         IF( MYROW.EQ.IAROW ) THEN
            CALL DGEBS2D( ICTXT, 'Columnwise', ' ', 1, NB, D( JJ ), 1 )
         ELSE
            CALL DGEBR2D( ICTXT, 'Columnwise', ' ', 1, NB, D( JJ ), 1,
     $                    IAROW, MYCOL )
         END IF
      END IF
*
      RETURN
*
*     End of PDLATRD
*
      END
