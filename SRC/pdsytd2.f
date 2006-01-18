      SUBROUTINE PDSYTD2( UPLO, N, A, IA, JA, DESCA, D, E, TAU, WORK,
     $                    LWORK, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION    A( * ), D( * ), E( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDSYTD2 reduces a real symmetric matrix sub( A ) to symmetric
*  tridiagonal form T by an orthogonal similarity transformation:
*  Q' * sub( A ) * Q = T, where sub( A ) = A(IA:IA+N-1,JA:JA+N-1).
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
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
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
*          triangular part is not referenced. On exit, if UPLO = 'U',
*          the diagonal and first superdiagonal of sub( A ) are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements above the first superdiagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of elementary reflectors; if UPLO = 'L', the diagonal
*          and first subdiagonal of sub( A ) are overwritten by the
*          corresponding elements of the tridiagonal matrix T, and the
*          elements below the first subdiagonal, with the array TAU,
*          represent the orthogonal matrix Q as a product of elementary
*          reflectors. See Further Details.
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
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                    dimension (LWORK)
*          On exit, WORK( 1 ) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= 3*N.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  INFO    (local output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(ia:ia+i-2,ja+i), and tau in TAU(ja+i-1).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in
*  A(ia+i+1:ia+n-1,ja+i-1), and tau in TAU(ja+i-1).
*
*  The contents of sub( A ) on exit are illustrated by the following
*  examples with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  d   e   v2  v3  v4 )              (  d                  )
*    (      d   e   v3  v4 )              (  e   d              )
*    (          d   e   v4 )              (  v1  e   d          )
*    (              d   e  )              (  v1  v2  e   d      )
*    (                  d  )              (  v1  v2  v3  e   d  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrix sub( A ) must verify some alignment proper-
*  ties, namely the following expression should be true:
*  ( MB_A.EQ.NB_A .AND. IROFFA.EQ.ICOFFA ) with
*  IROFFA = MOD( IA-1, MB_A ) and ICOFFA = MOD( JA-1, NB_A ).
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
      LOGICAL            LQUERY, UPPER
      INTEGER            IACOL, IAROW, ICOFFA, ICTXT, II, IK, IROFFA, J,
     $                   JJ, JK, JN, LDA, LWMIN, MYCOL, MYROW, NPCOL,
     $                   NPROW
      DOUBLE PRECISION   ALPHA, TAUI
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GRIDINFO, CHK1MAT, DAXPY,
     $                   DGEBR2D, DGEBS2D, DLARFG,
     $                   DSYMV, DSYR2, INFOG2L, PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
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
         UPPER = LSAME( UPLO, 'U' )
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
         LWMIN = 3 * N
*
         WORK( 1 ) = DBLE( LWMIN )
         LQUERY = ( LWORK.EQ.-1 )
         IF( INFO.EQ.0 ) THEN
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( IROFFA.NE.ICOFFA ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(600+NB_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -11
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDSYTD2', -INFO )
         CALL BLACS_ABORT( ICTXT, 1 )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Compute local information
*
      LDA = DESCA( LLD_ )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, II, JJ,
     $              IAROW, IACOL )
*
      IF( UPPER ) THEN
*
*        Process(IAROW, IACOL) owns block to be reduced
*
         IF( MYCOL.EQ.IACOL ) THEN
            IF( MYROW.EQ.IAROW ) THEN
*
*              Reduce the upper triangle of sub( A )
*
               DO 10 J = N-1, 1, -1
                  IK = II + J - 1
                  JK = JJ + J - 1
*
*                 Generate elementary reflector H(i) = I - tau * v * v'
*                 to annihilate A(IA:IA+J-1,JA:JA+J-1)
*
                  CALL DLARFG( J, A( IK+JK*LDA ), A( II+JK*LDA ), 1,
     $                         TAUI )
                  E( JK+1 ) = A( IK+JK*LDA )
*
                  IF( TAUI.NE.ZERO ) THEN
*
*                    Apply H(i) from both sides to
*                    A(IA:IA+J-1,JA:JA+J-1)
*
                     A( IK+JK*LDA ) = ONE
*
*                    Compute  x := tau * A * v  storing x in TAU(1:i)
*
                     CALL DSYMV( UPLO, J, TAUI, A( II+(JJ-1)*LDA ),
     $                           LDA, A( II+JK*LDA ), 1, ZERO,
     $                           TAU( JJ ), 1 )
*
*                    Compute  w := x - 1/2 * tau * (x'*v) * v
*
                     ALPHA = -HALF*TAUI*DDOT( J, TAU( JJ ), 1,
     $                                        A( II+JK*LDA ), 1 )
                     CALL DAXPY( J, ALPHA, A( II+JK*LDA ), 1,
     $                           TAU( JJ ), 1 )
*
*                    Apply the transformation as a rank-2 update:
*                       A := A - v * w' - w * v'
*
                     CALL DSYR2( UPLO, J, -ONE, A( II+JK*LDA ), 1,
     $                           TAU( JJ ), 1, A( II+(JJ-1)*LDA ),
     $                           LDA )
                     A( IK+JK*LDA ) = E( JK+1 )
                  END IF
*
*                 Copy D, E, TAU to broadcast them columnwise.
*
                  D( JK+1 ) = A( IK+1+JK*LDA )
                  WORK( J+1 ) = D( JK+1 )
                  WORK( N+J+1 ) = E( JK+1 )
                  TAU( JK+1 ) = TAUI
                  WORK( 2*N+J+1 ) = TAU( JK+1 )
*
   10          CONTINUE
               D( JJ ) = A( II+(JJ-1)*LDA )
               WORK( 1 ) = D( JJ )
               WORK( N+1 ) = ZERO
               WORK( 2*N+1 ) = ZERO
*
               CALL DGEBS2D( ICTXT, 'Columnwise', ' ', 1, 3*N, WORK, 1 )
*
            ELSE
               CALL DGEBR2D( ICTXT, 'Columnwise', ' ', 1, 3*N, WORK, 1,
     $                       IAROW, IACOL )
               DO 20 J = 2, N
                  JN = JJ + J - 1
                  D( JN ) = WORK( J )
                  E( JN ) = WORK( N+J )
                  TAU( JN ) = WORK( 2*N+J )
   20          CONTINUE
               D( JJ ) = WORK( 1 )
            END IF
         END IF
*
      ELSE
*
*        Process (IAROW, IACOL) owns block to be factorized
*
         IF( MYCOL.EQ.IACOL ) THEN
            IF( MYROW.EQ.IAROW ) THEN
*
*              Reduce the lower triangle of sub( A )
*
               DO 30 J = 1, N - 1
                  IK = II + J - 1
                  JK = JJ + J - 1
*
*                 Generate elementary reflector H(i) = I - tau * v * v'
*                 to annihilate A(IA+J-JA+2:IA+N-1,JA+J-1)
*
                  CALL DLARFG( N-J, A( IK+1+(JK-1)*LDA ),
     $                         A( IK+2+(JK-1)*LDA ), 1, TAUI )
                  E( JK ) = A( IK+1+(JK-1)*LDA )
*
                  IF( TAUI.NE.ZERO ) THEN
*
*                    Apply H(i) from both sides to
*                    A(IA+J-JA+1:IA+N-1,JA+J+1:JA+N-1)
*
                     A( IK+1+(JK-1)*LDA ) = ONE
*
*                    Compute  x := tau * A * v  storing y in TAU(i:n-1)
*
                     CALL DSYMV( UPLO, N-J, TAUI, A( IK+1+JK*LDA ),
     $                           LDA, A( IK+1+(JK-1)*LDA ), 1,
     $                           ZERO, TAU( JK ), 1 )
*
*                    Compute  w := x - 1/2 * tau * (x'*v) * v
*
                     ALPHA = -HALF*TAUI*DDOT( N-J, TAU( JK ), 1,
     $                        A( IK+1+(JK-1)*LDA ), 1 )
                     CALL DAXPY( N-J, ALPHA, A( IK+1+(JK-1)*LDA ),
     $                           1, TAU( JK ), 1 )
*
*                    Apply the transformation as a rank-2 update:
*                       A := A - v * w' - w * v'
*
                     CALL DSYR2( UPLO, N-J, -ONE,
     $                           A( IK+1+(JK-1)*LDA ), 1,
     $                           TAU( JK ), 1, A( IK+1+JK*LDA ),
     $                           LDA )
                     A( IK+1+(JK-1)*LDA ) = E( JK )
                  END IF
*
*                 Copy D(JK), E(JK), TAU(JK) to broadcast them
*                 columnwise.
*
                  D( JK ) = A( IK+(JK-1)*LDA )
                  WORK( J ) = D( JK )
                  WORK( N+J ) = E( JK )
                  TAU( JK ) = TAUI
                  WORK( 2*N+J ) = TAU( JK )
   30          CONTINUE
               JN = JJ + N - 1
               D( JN ) = A( II+N-1+(JN-1)*LDA )
               WORK( N ) = D( JN )
               TAU( JN ) = ZERO
               WORK( 2*N ) = ZERO
*
               CALL DGEBS2D( ICTXT, 'Columnwise', ' ', 1, 3*N-1, WORK,
     $                            1 )
*
            ELSE
               CALL DGEBR2D( ICTXT, 'Columnwise', ' ', 1, 3*N-1, WORK,
     $                       1, IAROW, IACOL )
               DO 40 J = 1, N - 1
                  JN = JJ + J - 1
                  D( JN ) = WORK( J )
                  E( JN ) = WORK( N+J )
                  TAU( JN ) = WORK( 2*N+J )
   40          CONTINUE
               JN = JJ + N - 1
               D( JN ) = WORK( N )
               TAU( JN ) = ZERO
            END IF
         END IF
      END IF
*
      WORK( 1 ) = DBLE( LWMIN )
*
      RETURN
*
*     End of PDSYTD2
*
      END
