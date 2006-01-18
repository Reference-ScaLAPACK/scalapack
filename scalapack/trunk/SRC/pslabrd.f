      SUBROUTINE PSLABRD( M, N, NB, A, IA, JA, DESCA, D, E, TAUQ, TAUP,
     $                    X, IX, JX, DESCX, Y, IY, JY, DESCY, WORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER             IA, IX, IY, JA, JX, JY, M, N, NB
*     ..
*     .. Array Arguments ..
      INTEGER             DESCA( * ), DESCX( * ), DESCY( * )
      REAL                A( * ), D( * ), E( * ), TAUP( * ),
     $                    TAUQ( * ), X( * ), Y( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLABRD reduces the first NB rows and columns of a real general
*  M-by-N distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) to upper
*  or lower bidiagonal form by an orthogonal transformation Q' * A * P,
*  and returns the matrices X and Y which are needed to apply the
*  transformation to the unreduced part of sub( A ).
*
*  If M >= N, sub( A ) is reduced to upper bidiagonal form; if M < N, to
*  lower bidiagonal form.
*
*  This is an auxiliary routine called by PSGEBRD.
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
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  NB      (global input) INTEGER
*          The number of leading rows and columns of sub( A ) to be
*          reduced.
*
*  A       (local input/local output) REAL pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          general distributed matrix sub( A ) to be reduced. On exit,
*          the first NB rows and columns of the matrix are overwritten;
*          the rest of the distributed matrix sub( A ) is unchanged.
*          If m >= n, elements on and below the diagonal in the first NB
*            columns, with the array TAUQ, represent the orthogonal
*            matrix Q as a product of elementary reflectors; and
*            elements above the diagonal in the first NB rows, with the
*            array TAUP, represent the orthogonal matrix P as a product
*            of elementary reflectors.
*          If m < n, elements below the diagonal in the first NB
*            columns, with the array TAUQ, represent the orthogonal
*            matrix Q as a product of elementary reflectors, and
*            elements on and above the diagonal in the first NB rows,
*            with the array TAUP, represent the orthogonal matrix P as
*            a product of elementary reflectors.
*          See Further Details.
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
*  D       (local output) REAL array, dimension
*          LOCr(IA+MIN(M,N)-1) if M >= N; LOCc(JA+MIN(M,N)-1) otherwise.
*          The distributed diagonal elements of the bidiagonal matrix
*          B: D(i) = A(ia+i-1,ja+i-1). D is tied to the distributed
*          matrix A.
*
*  E       (local output) REAL array, dimension
*          LOCr(IA+MIN(M,N)-1) if M >= N; LOCc(JA+MIN(M,N)-2) otherwise.
*          The distributed off-diagonal elements of the bidiagonal
*          distributed matrix B:
*          if m >= n, E(i) = A(ia+i-1,ja+i) for i = 1,2,...,n-1;
*          if m < n, E(i) = A(ia+i,ja+i-1) for i = 1,2,...,m-1.
*          E is tied to the distributed matrix A.
*
*  TAUQ    (local output) REAL array dimension
*          LOCc(JA+MIN(M,N)-1). The scalar factors of the elementary
*          reflectors which represent the orthogonal matrix Q. TAUQ
*          is tied to the distributed matrix A. See Further Details.
*
*  TAUP    (local output) REAL array, dimension
*          LOCr(IA+MIN(M,N)-1). The scalar factors of the elementary
*          reflectors which represent the orthogonal matrix P. TAUP
*          is tied to the distributed matrix A. See Further Details.
*
*  X       (local output) REAL pointer into the local memory
*          to an array of dimension (LLD_X,NB). On exit, the local
*          pieces of the distributed M-by-NB matrix
*          X(IX:IX+M-1,JX:JX+NB-1) required to update the unreduced
*          part of sub( A ).
*
*  IX      (global input) INTEGER
*          The row index in the global array X indicating the first
*          row of sub( X ).
*
*  JX      (global input) INTEGER
*          The column index in the global array X indicating the
*          first column of sub( X ).
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*
*  Y       (local output) REAL pointer into the local memory
*          to an array of dimension (LLD_Y,NB).  On exit, the local
*          pieces of the distributed N-by-NB matrix
*          Y(IY:IY+N-1,JY:JY+NB-1) required to update the unreduced
*          part of sub( A ).
*
*  IY      (global input) INTEGER
*          The row index in the global array Y indicating the first
*          row of sub( Y ).
*
*  JY      (global input) INTEGER
*          The column index in the global array Y indicating the
*          first column of sub( Y ).
*
*  DESCY   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Y.
*
*  WORK    (local workspace) REAL array, dimension (LWORK)
*          LWORK >= NB_A + NQ, with
*
*          NQ = NUMROC( N+MOD( IA-1, NB_Y ), NB_Y, MYCOL, IACOL, NPCOL )
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL )
*
*          INDXG2P and NUMROC are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*  Further Details
*  ===============
*
*  The matrices Q and P are represented as products of elementary
*  reflectors:
*
*     Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are real scalars, and v and u are real vectors.
*
*  If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in
*  A(ia+i-1:ia+m-1,ja+i-1); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is
*  stored on exit in A(ia+i-1,ja+i:ja+n-1); tauq is stored in
*  TAUQ(ja+i-1) and taup in TAUP(ia+i-1).
*
*  If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in
*  A(ia+i+1:ia+m-1,ja+i-1); u(1:i-1) = 0, u(i) = 1, and u(i:n) is
*  stored on exit in A(ia+i-1,ja+i:ja+n-1); tauq is stored in
*  TAUQ(ja+i-1) and taup in TAUP(ia+i-1).
*
*  The elements of the vectors v and u together form the m-by-nb matrix
*  V and the nb-by-n matrix U' which are needed, with X and Y, to apply
*  the transformation to the unreduced part of the matrix, using a block
*  update of the form:  sub( A ) := sub( A ) - V*Y' - X*U'.
*
*  The contents of sub( A ) on exit are illustrated by the following
*  examples with nb = 2:
*
*  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
*
*    (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )
*    (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )
*    (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )
*    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
*    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
*    (  v1  v2  a   a   a  )
*
*  where a denotes an element of the original matrix which is unchanged,
*  vi denotes an element of the vector defining H(i), and ui an element
*  of the vector defining G(i).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IACOL, IAROW, ICTXT, II, IPY, IW, J, JJ,
     $                   JWY, K, MYCOL, MYROW, NPCOL, NPROW
      REAL               ALPHA, TAU
      INTEGER            DESCD( DLEN_ ), DESCE( DLEN_ ),
     $                   DESCTP( DLEN_ ), DESCTQ( DLEN_ ),
     $                   DESCW( DLEN_ ), DESCWY( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, INFOG2L, PSCOPY,
     $                   PSELGET, PSELSET, PSGEMV, PSLARFG,
     $                   PSSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, II, JJ,
     $              IAROW, IACOL )
      IPY = DESCA( MB_ ) + 1
      IW = MOD( IA-1, DESCA( NB_ ) ) + 1
      ALPHA = ZERO
*
      CALL DESCSET( DESCWY, 1, N+MOD( IA-1, DESCY( NB_ ) ), 1,
     $              DESCA( NB_ ), IAROW, IACOL, ICTXT, 1 )
      CALL DESCSET( DESCW, DESCA( MB_ ), 1, DESCA( MB_ ), 1, IAROW,
     $              IACOL, ICTXT, DESCA( MB_ ) )
      CALL DESCSET( DESCTQ, 1, JA+MIN(M,N)-1, 1, DESCA( NB_ ), IAROW,
     $              DESCA( CSRC_ ), DESCA( CTXT_ ), 1 )
      CALL DESCSET( DESCTP, IA+MIN(M,N)-1, 1, DESCA( MB_ ), 1,
     $              DESCA( RSRC_ ), IACOL, DESCA( CTXT_ ),
     $              DESCA( LLD_ ) )
*
      IF( M.GE.N ) THEN
*
*        Reduce to upper bidiagonal form
*
         CALL DESCSET( DESCD, 1, JA+MIN(M,N)-1, 1, DESCA( NB_ ), MYROW,
     $                 DESCA( CSRC_ ), DESCA( CTXT_ ), 1 )
         CALL DESCSET( DESCE, IA+MIN(M,N)-1, 1, DESCA( MB_ ), 1,
     $                 DESCA( RSRC_ ), MYCOL, DESCA( CTXT_ ),
     $                 DESCA( LLD_ ) )
         DO 10 K = 1, NB
            I = IA + K - 1
            J = JA + K - 1
            JWY = IW + K
*
*           Update A(i:ia+m-1,j)
*
            IF( K.GT.1 ) THEN
               CALL PSGEMV( 'No transpose', M-K+1, K-1, -ONE, A, I, JA,
     $                      DESCA, Y, IY, JY+K-1, DESCY, 1, ONE, A, I,
     $                      J, DESCA, 1 )
               CALL PSGEMV( 'No transpose', M-K+1, K-1, -ONE, X, IX+K-1,
     $                      JX, DESCX, A, IA, J, DESCA, 1, ONE, A, I, J,
     $                      DESCA, 1 )
               CALL PSELSET( A, I-1, J, DESCA, ALPHA )
            END IF
*
*           Generate reflection Q(i) to annihilate A(i+1:ia+m-1,j)
*
            CALL PSLARFG( M-K+1, ALPHA, I, J, A, I+1, J, DESCA, 1,
     $                    TAUQ )
            CALL PSELSET( D, 1, J, DESCD, ALPHA )
            CALL PSELSET( A, I, J, DESCA, ONE )
*
*           Compute Y(IA+I:IA+N-1,J)
*
            CALL PSGEMV( 'Transpose', M-K+1, N-K, ONE, A, I, J+1, DESCA,
     $                   A, I, J, DESCA, 1, ZERO, WORK( IPY ), 1, JWY,
     $                   DESCWY, DESCWY( M_ ) )
            CALL PSGEMV( 'Transpose', M-K+1, K-1, ONE, A, I, JA, DESCA,
     $                   A, I, J, DESCA, 1, ZERO, WORK, IW, 1, DESCW,
     $                   1 )
            CALL PSGEMV( 'Transpose', K-1, N-K, -ONE, Y, IY, JY+K,
     $                   DESCY, WORK, IW, 1, DESCW, 1, ONE, WORK( IPY ),
     $                   1, JWY, DESCWY, DESCWY( M_ ) )
            CALL PSGEMV( 'Transpose', M-K+1, K-1, ONE, X, IX+K-1, JX,
     $                   DESCX, A, I, J, DESCA, 1, ZERO, WORK, IW, 1,
     $                   DESCW, 1 )
            CALL PSGEMV( 'Transpose', K-1, N-K, -ONE, A, IA, J+1, DESCA,
     $                   WORK, IW, 1, DESCW, 1, ONE, WORK( IPY ), 1,
     $                   JWY, DESCWY, DESCWY( M_ ) )
*
            CALL PSELGET( 'Rowwise', ' ', TAU, TAUQ, 1, J, DESCTQ )
            CALL PSSCAL( N-K, TAU, WORK( IPY ), 1, JWY, DESCWY,
     $                   DESCWY( M_ ) )
            CALL PSCOPY( N-K, WORK( IPY ), 1, JWY, DESCWY, DESCWY( M_ ),
     $                   Y, IY+K-1, JY+K, DESCY, DESCY( M_ ) )
*
*           Update A(i,j+1:ja+n-1)
*
            CALL PSGEMV( 'Transpose', K, N-K, -ONE, Y, IY, JY+K, DESCY,
     $                   A, I, JA, DESCA, DESCA( M_ ), ONE, A, I, J+1,
     $                   DESCA, DESCA( M_ ) )
            CALL PSGEMV( 'Transpose', K-1, N-K, -ONE, A, IA, J+1, DESCA,
     $                   X, IX+K-1, JX, DESCX, DESCX( M_ ), ONE, A, I,
     $                   J+1, DESCA, DESCA( M_ ) )
            CALL PSELSET( A, I, J, DESCA, ALPHA )
*
*           Generate reflection P(i) to annihilate A(i,j+2:ja+n-1)
*
            CALL PSLARFG( N-K, ALPHA, I, J+1, A, I,
     $                    MIN( J+2, N+JA-1 ), DESCA, DESCA( M_ ), TAUP )
            CALL PSELSET( E, I, 1, DESCE, ALPHA )
            CALL PSELSET( A, I, J+1, DESCA, ONE )
*
*           Compute X(I+1:IA+M-1,J)
*
            CALL PSGEMV( 'No transpose', M-K, N-K, ONE, A, I+1, J+1,
     $                   DESCA, A, I, J+1, DESCA, DESCA( M_ ), ZERO, X,
     $                   IX+K, JX+K-1, DESCX, 1 )
            CALL PSGEMV( 'No transpose', K, N-K, ONE, Y, IY, JY+K,
     $                   DESCY, A, I, J+1, DESCA, DESCA( M_ ), ZERO,
     $                   WORK, IW, 1, DESCW, 1 )
            CALL PSGEMV( 'No transpose', M-K, K, -ONE, A, I+1, JA,
     $                   DESCA, WORK, IW, 1, DESCW, 1, ONE, X, IX+K,
     $                   JX+K-1, DESCX, 1 )
            CALL PSGEMV( 'No transpose', K-1, N-K, ONE, A, IA, J+1,
     $                   DESCA, A, I, J+1, DESCA, DESCA( M_ ), ZERO,
     $                   WORK, IW, 1, DESCW, 1 )
            CALL PSGEMV( 'No transpose', M-K, K-1, -ONE, X, IX+K, JX,
     $                   DESCX, WORK, IW, 1, DESCW, 1, ONE, X, IX+K,
     $                   JX+K-1, DESCX, 1 )
*
            CALL PSELGET( 'Columnwise', ' ', TAU, TAUP, I, 1, DESCTP )
            CALL PSSCAL( M-K, TAU, X, IX+K, JX+K-1, DESCX, 1 )
   10    CONTINUE
*
      ELSE
*
*        Reduce to lower bidiagonal form
*
         CALL DESCSET( DESCD, IA+MIN(M,N)-1, 1, DESCA( MB_ ), 1,
     $                 DESCA( RSRC_ ), MYCOL, DESCA( CTXT_ ),
     $                 DESCA( LLD_ ) )
         CALL DESCSET( DESCE, 1, JA+MIN(M,N)-1, 1, DESCA( NB_ ), MYROW,
     $                 DESCA( CSRC_ ), DESCA( CTXT_ ), 1 )
         DO 20 K = 1, NB
            I = IA + K - 1
            J = JA + K - 1
            JWY = IW + K
*
*           Update A(i,j:ja+n-1)
*
            IF( K.GT.1 ) THEN
               CALL PSGEMV( 'Transpose', K-1, N-K+1, -ONE, Y, IY,
     $                      JY+K-1, DESCY, A, I, JA, DESCA, DESCA( M_ ),
     $                      ONE, A, I, J, DESCA, DESCA( M_ ) )
               CALL PSGEMV( 'Transpose', K-1, N-K+1, -ONE, A, IA, J,
     $                      DESCA, X, IX+K-1, JX, DESCX, DESCX( M_ ),
     $                      ONE, A, I, J, DESCA, DESCA( M_ ) )
               CALL PSELSET( A, I, J-1, DESCA, ALPHA )
            END IF
*
*           Generate reflection P(i) to annihilate A(i,j+1:ja+n-1)
*
            CALL PSLARFG( N-K+1, ALPHA, I, J, A, I, J+1, DESCA,
     $                    DESCA( M_ ), TAUP )
            CALL PSELSET( D, I, 1, DESCD, ALPHA )
            CALL PSELSET( A, I, J, DESCA, ONE )
*
*           Compute X(i+1:ia+m-1,j)
*
            CALL PSGEMV( 'No transpose', M-K, N-K+1, ONE, A, I+1, J,
     $                   DESCA, A, I, J, DESCA, DESCA( M_ ), ZERO, X,
     $                   IX+K, JX+K-1, DESCX, 1 )
            CALL PSGEMV( 'No transpose', K-1, N-K+1, ONE, Y, IY, JY+K-1,
     $                   DESCY, A, I, J, DESCA, DESCA( M_ ), ZERO,
     $                   WORK, IW, 1, DESCW, 1 )
            CALL PSGEMV( 'No transpose', M-K, K-1, -ONE, A, I+1, JA,
     $                   DESCA, WORK, IW, 1, DESCW, 1, ONE, X, IX+K,
     $                   JX+K-1, DESCX, 1 )
            CALL PSGEMV( 'No transpose', K-1, N-K+1, ONE, A, IA, J,
     $                   DESCA, A, I, J, DESCA, DESCA( M_ ), ZERO,
     $                   WORK, IW, 1, DESCW, 1 )
            CALL PSGEMV( 'No transpose', M-K, K-1, -ONE, X, IX+K, JX,
     $                   DESCX, WORK, IW, 1, DESCW, 1, ONE, X, IX+K,
     $                   JX+K-1, DESCX, 1 )
*
            CALL PSELGET( 'Columnwise', ' ', TAU, TAUP, I, 1, DESCTP )
            CALL PSSCAL( M-K, TAU, X, IX+K, JX+K-1, DESCX, 1 )
*
*           Update A(i+1:ia+m-1,j)
*
            CALL PSGEMV( 'No transpose', M-K, K-1, -ONE, A, I+1, JA,
     $                   DESCA, Y, IY, JY+K-1, DESCY, 1, ONE, A, I+1, J,
     $                   DESCA, 1 )
            CALL PSGEMV( 'No transpose', M-K, K, -ONE, X, IX+K, JX,
     $                   DESCX, A, IA, J, DESCA, 1, ONE, A, I+1, J,
     $                   DESCA, 1 )
            CALL PSELSET( A, I, J, DESCA, ALPHA )
*
*           Generate reflection Q(i) to annihilate A(i+2:ia+m-1,j)
*
            CALL PSLARFG( M-K, ALPHA, I+1, J, A, MIN( I+2, M+IA-1 ),
     $                    J, DESCA, 1, TAUQ )
            CALL PSELSET( E, 1, J, DESCE, ALPHA )
            CALL PSELSET( A, I+1, J, DESCA, ONE )
*
*           Compute Y(ia+i:ia+n-1,j)
*
            CALL PSGEMV( 'Transpose', M-K, N-K, ONE, A, I+1, J+1, DESCA,
     $                   A, I+1, J, DESCA, 1, ZERO, WORK( IPY ), 1,
     $                   JWY, DESCWY, DESCWY( M_ ) )
            CALL PSGEMV( 'Transpose', M-K, K-1, ONE, A, I+1, JA, DESCA,
     $                   A, I+1, J, DESCA, 1, ZERO, WORK, IW, 1, DESCW,
     $                   1 )
            CALL PSGEMV( 'Transpose', K-1, N-K, -ONE, Y, IY, JY+K,
     $                   DESCY, WORK, IW, 1, DESCW, 1, ONE, WORK( IPY ),
     $                   1, JWY, DESCWY, DESCWY( M_ ) )
            CALL PSGEMV( 'Transpose', M-K, K, ONE, X, IX+K, JX, DESCX,
     $                   A, I+1, J, DESCA, 1, ZERO, WORK, IW, 1, DESCW,
     $                   1 )
            CALL PSGEMV( 'Transpose', K, N-K, -ONE, A, IA, J+1, DESCA,
     $                   WORK, IW, 1, DESCW, 1, ONE, WORK( IPY ), 1,
     $                   JWY, DESCWY, DESCWY( M_ ) )
*
            CALL PSELGET( 'Rowwise', ' ', TAU, TAUQ, 1, J, DESCTQ )
            CALL PSSCAL( N-K, TAU, WORK( IPY ), 1, JWY, DESCWY,
     $                   DESCWY( M_ ) )
            CALL PSCOPY( N-K, WORK( IPY ), 1, JWY, DESCWY, DESCWY( M_ ),
     $                   Y, IY+K-1, JY+K, DESCY, DESCY( M_ ) )
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of PSLABRD
*
      END
