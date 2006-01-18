      SUBROUTINE PZGEBRD( M, N, A, IA, JA, DESCA, D, E, TAUQ, TAUP,
     $                    WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( * ), TAUP( * ), TAUQ( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZGEBRD reduces a complex general M-by-N distributed matrix
*  sub( A ) = A(IA:IA+M-1,JA:JA+N-1) to upper or lower bidiagonal
*  form B by an unitary transformation: Q' * sub( A ) * P = B.
*
*  If M >= N, B is upper bidiagonal; if M < N, B is lower bidiagonal.
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
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          general distributed matrix sub( A ). On exit, if M >= N,
*          the diagonal and the first superdiagonal of sub( A ) are
*          overwritten with the upper bidiagonal matrix B; the elements
*          below the diagonal, with the array TAUQ, represent the
*          unitary matrix Q as a product of elementary reflectors, and
*          the elements above the first superdiagonal, with the array
*          TAUP, represent the orthogonal matrix P as a product of
*          elementary reflectors. If M < N, the diagonal and the first
*          subdiagonal are overwritten with the lower bidiagonal
*          matrix B; the elements below the first subdiagonal, with the
*          array TAUQ, represent the unitary matrix Q as a product of
*          elementary reflectors, and the elements above the diagonal,
*          with the array TAUP, represent the orthogonal matrix P as a
*          product of elementary reflectors. See Further Details.
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
*  D       (local output) DOUBLE PRECISION array, dimension
*          LOCc(JA+MIN(M,N)-1) if M >= N; LOCr(IA+MIN(M,N)-1) otherwise.
*          The distributed diagonal elements of the bidiagonal matrix
*          B: D(i) = A(i,i). D is tied to the distributed matrix A.
*
*  E       (local output) DOUBLE PRECISION array, dimension
*          LOCr(IA+MIN(M,N)-1) if M >= N; LOCc(JA+MIN(M,N)-2) otherwise.
*          The distributed off-diagonal elements of the bidiagonal
*          distributed matrix B:
*          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
*          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
*          E is tied to the distributed matrix A.
*
*  TAUQ    (local output) COMPLEX*16 array dimension
*          LOCc(JA+MIN(M,N)-1). The scalar factors of the elementary
*          reflectors which represent the unitary matrix Q. TAUQ is
*          tied to the distributed matrix A. See Further Details.
*
*  TAUP    (local output) COMPLEX*16 array, dimension
*          LOCr(IA+MIN(M,N)-1). The scalar factors of the elementary
*          reflectors which represent the unitary matrix P. TAUP is
*          tied to the distributed matrix A. See Further Details.
*
*  WORK    (local workspace/local output) COMPLEX*16 array,
*                                                 dimension (LWORK)
*          On exit, WORK( 1 ) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= NB*( MpA0 + NqA0 + 1 ) + NqA0
*
*          where NB = MB_A = NB_A,
*          IROFFA = MOD( IA-1, NB ), ICOFFA = MOD( JA-1, NB ),
*          IAROW = INDXG2P( IA, NB, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB, MYCOL, CSRC_A, NPCOL ),
*          MpA0 = NUMROC( M+IROFFA, NB, MYROW, IAROW, NPROW ),
*          NqA0 = NUMROC( N+ICOFFA, NB, MYCOL, IACOL, NPCOL ).
*
*          INDXG2P and NUMROC are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  Further Details
*  ===============
*
*  The matrices Q and P are represented as products of elementary
*  reflectors:
*
*  If m >= n,
*
*     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are complex scalars, and v and u are complex
*  vectors;
*  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in
*  A(ia+i:ia+m-1,ja+i-1);
*  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in
*  A(ia+i-1,ja+i+1:ja+n-1);
*  tauq is stored in TAUQ(ja+i-1) and taup in TAUP(ia+i-1).
*
*  If m < n,
*
*     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are complex scalars, and v and u are complex
*  vectors;
*  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in
*  A(ia+i+1:ia+m-1,ja+i-1);
*  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in
*  A(ia+i-1,ja+i:ja+n-1);
*  tauq is stored in TAUQ(ja+i-1) and taup in TAUP(ia+i-1).
*
*  The contents of sub( A ) on exit are illustrated by the following
*  examples:
*
*  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
*
*    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
*    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
*    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
*    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
*    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
*    (  v1  v2  v3  v4  v5 )
*
*  where d and e denote diagonal and off-diagonal elements of B, vi
*  denotes an element of the vector defining H(i), and ui an element of
*  the vector defining G(i).
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrix sub( A ) must verify some alignment proper-
*  ties, namely the following expressions should be true:
*  ( MB_A.EQ.NB_A .AND. IROFFA.EQ.ICOFFA )
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      CHARACTER          COLCTOP, ROWCTOP
      INTEGER            I, IACOL, IAROW, ICTXT, IINFO, IOFF, IPW, IPY,
     $                   IW, J, JB, JS, JW, K, L, LWMIN, MN, MP, MYCOL,
     $                   MYROW, NB, NPCOL, NPROW, NQ
*     ..
*     .. Local Arrays ..
      INTEGER            DESCWX( DLEN_ ), DESCWY( DLEN_ ), IDUM1( 1 ),
     $                   IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, PCHK1MAT,
     $                   PB_TOPGET, PB_TOPSET, PXERBLA, PZELSET,
     $                   PZGEBD2, PZGEMM, PZLABRD
*     ..
*     .. External Functions ..
      INTEGER            INDXG2L, INDXG2P, NUMROC
      EXTERNAL           INDXG2L, INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DBLE, MAX, MIN, MOD
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
         CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, INFO )
         IF( INFO.EQ.0 ) THEN
            NB = DESCA( MB_ )
            IOFF = MOD( IA-1, DESCA( MB_ ) )
            IAROW = INDXG2P( IA, NB, MYROW, DESCA( RSRC_ ), NPROW )
            IACOL = INDXG2P( JA, NB, MYCOL, DESCA( CSRC_ ), NPCOL )
            MP = NUMROC( M+IOFF, NB, MYROW, IAROW, NPROW )
            NQ = NUMROC( N+IOFF, NB, MYCOL, IACOL, NPCOL )
            LWMIN = NB*( MP+NQ+1 ) + NQ
*
            WORK( 1 ) = DCMPLX( DBLE( LWMIN ) )
            LQUERY = ( LWORK.EQ.-1 )
            IF( IOFF.NE.MOD( JA-1, DESCA( NB_ ) ) ) THEN
               INFO = -5
            ELSE IF( NB.NE.DESCA( NB_ ) ) THEN
               INFO = -(600+NB_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -12
            END IF
         END IF
         IF( LQUERY ) THEN
            IDUM1( 1 ) = -1
         ELSE
            IDUM1( 1 ) = 1
         END IF
         IDUM2( 1 ) = 12
         CALL PCHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, 1, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.LT.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZGEBRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      MN = MIN( M, N )
      IF( MN.EQ.0 )
     $   RETURN
*
*     Initialize parameters.
*
      CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
      CALL PB_TOPGET( ICTXT, 'Combine', 'Rowwise',    ROWCTOP )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', '1-tree' )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise',    '1-tree' )
*
      IPY = MP * NB + 1
      IPW = NQ * NB + IPY
*
      CALL DESCSET( DESCWX, M+IOFF, NB, NB, NB, IAROW, IACOL, ICTXT,
     $              MAX( 1, MP ) )
      CALL DESCSET( DESCWY, NB, N+IOFF, NB, NB, IAROW, IACOL, ICTXT,
     $              NB )
*
      MP = NUMROC( M+IA-1, NB, MYROW, DESCA( RSRC_ ), NPROW )
      NQ = NUMROC( N+JA-1, NB, MYCOL, DESCA( CSRC_ ), NPCOL )
      K  = 1
      JB = NB - IOFF
      IW = IOFF + 1
      JW = IOFF + 1
*
      DO 10 L = 1, MN+IOFF-NB, NB
         I = IA + K - 1
         J = JA + K - 1
*
*        Reduce rows and columns i:i+nb-1 to bidiagonal form and return
*        the matrices X and Y which are needed to update the unreduced
*        part of the matrix.
*
         CALL PZLABRD( M-K+1, N-K+1, JB, A, I, J, DESCA, D, E, TAUQ,
     $                 TAUP, WORK, IW, JW, DESCWX, WORK( IPY ), IW,
     $                 JW, DESCWY, WORK( IPW ) )
*
*        Update the trailing submatrix A(i+nb:ia+m-1,j+nb:ja+n-1), using
*        an update of the form  A := A - V*Y' - X*U'.
*
         CALL PZGEMM( 'No transpose', 'No transpose', M-K-JB+1,
     $                N-K-JB+1, JB, -ONE, A, I+JB, J, DESCA,
     $                WORK( IPY ), IW, JW+JB, DESCWY, ONE, A, I+JB,
     $                J+JB, DESCA )
         CALL PZGEMM( 'No transpose', 'No transpose', M-K-JB+1,
     $                N-K-JB+1, JB, -ONE, WORK, IW+JB, JW, DESCWX, A, I,
     $                J+JB, DESCA, ONE, A, I+JB, J+JB, DESCA )
*
*        Copy last off-diagonal elements of B back into sub( A ).
*
         IF( M.GE.N ) THEN
            JS = MIN( INDXG2L( I+JB-1, NB, 0, DESCA( RSRC_ ), NPROW ),
     $                MP )
            IF( JS.GT.0 )
     $         CALL PZELSET( A, I+JB-1, J+JB, DESCA, DCMPLX( E( JS ) ) )
         ELSE
            JS = MIN( INDXG2L( J+JB-1, NB, 0, DESCA( CSRC_ ), NPCOL ),
     $                NQ )
            IF( JS.GT.0 )
     $         CALL PZELSET( A, I+JB, J+JB-1, DESCA, DCMPLX( E( JS ) ) )
         END IF
*
         K = K + JB
         JB = NB
         IW = 1
         JW = 1
         DESCWX( M_ ) = DESCWX( M_ ) - JB
         DESCWX( RSRC_ ) = MOD( DESCWX( RSRC_ ) + 1, NPROW )
         DESCWX( CSRC_ ) = MOD( DESCWX( CSRC_ ) + 1, NPCOL )
         DESCWY( N_ ) = DESCWY( N_ ) - JB
         DESCWY( RSRC_ ) = MOD( DESCWY( RSRC_ ) + 1, NPROW )
         DESCWY( CSRC_ ) = MOD( DESCWY( CSRC_ ) + 1, NPCOL )
*
   10 CONTINUE
*
*     Use unblocked code to reduce the remainder of the matrix.
*
      CALL PZGEBD2( M-K+1, N-K+1, A, IA+K-1, JA+K-1, DESCA, D, E, TAUQ,
     $              TAUP, WORK, LWORK, IINFO )
*
      CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise',    ROWCTOP )
*
      WORK( 1 ) = DCMPLX( DBLE( LWMIN ) )
*
      RETURN
*
*     End of PZGEBRD
*
      END
