      SUBROUTINE PSGEBD2( M, N, A, IA, JA, DESCA, D, E, TAUQ, TAUP,
     $                    WORK, LWORK, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      REAL               A( * ), D( * ), E( * ), TAUP( * ), TAUQ( * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSGEBD2 reduces a real general M-by-N distributed matrix
*  sub( A ) = A(IA:IA+M-1,JA:JA+N-1) to upper or lower bidiagonal
*  form B by an orthogonal transformation: Q' * sub( A ) * P = B.
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
*  A       (local input/local output) REAL pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          general distributed matrix sub( A ). On exit, if M >= N,
*          the diagonal and the first superdiagonal of sub( A ) are
*          overwritten with the upper bidiagonal matrix B; the elements
*          below the diagonal, with the array TAUQ, represent the
*          orthogonal matrix Q as a product of elementary reflectors,
*          and the elements above the first superdiagonal, with the
*          array TAUP, represent the orthogonal matrix P as a product
*          of elementary reflectors. If M < N, the diagonal and the
*          first subdiagonal are overwritten with the lower bidiagonal
*          matrix B; the elements below the first subdiagonal, with the
*          array TAUQ, represent the orthogonal matrix Q as a product of
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
*  D       (local output) REAL array, dimension
*          LOCc(JA+MIN(M,N)-1) if M >= N; LOCr(IA+MIN(M,N)-1) otherwise.
*          The distributed diagonal elements of the bidiagonal matrix
*          B: D(i) = A(i,i). D is tied to the distributed matrix A.
*
*  E       (local output) REAL array, dimension
*          LOCr(IA+MIN(M,N)-1) if M >= N; LOCc(JA+MIN(M,N)-2) otherwise.
*          The distributed off-diagonal elements of the bidiagonal
*          distributed matrix B:
*          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
*          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
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
*  WORK    (local workspace/local output) REAL array,
*                                                  dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= MAX( MpA0, NqA0 )
*
*          where NB = MB_A = NB_A, IROFFA = MOD( IA-1, NB )
*          IAROW = INDXG2P( IA, NB, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB, MYCOL, CSRC_A, NPCOL ),
*          MpA0 = NUMROC( M+IROFFA, NB, MYROW, IAROW, NPROW ),
*          NqA0 = NUMROC( N+IROFFA, NB, MYCOL, IACOL, NPCOL ).
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
*  where tauq and taup are real scalars, and v and u are real vectors;
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
*  where tauq and taup are real scalars, and v and u are real vectors;
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
*                  ( MB_A.EQ.NB_A .AND. IROFFA.EQ.ICOFFA )
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
      LOGICAL            LQUERY
      INTEGER            I, IACOL, IAROW, ICOFFA, ICTXT, II, IROFFA, J,
     $                   JJ, K, LWMIN, MPA0, MYCOL, MYROW, NPCOL, NPROW,
     $                   NQA0
      REAL               ALPHA
*     ..
*     .. Local Arrays ..
      INTEGER            DESCD( DLEN_ ), DESCE( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GRIDINFO, CHK1MAT, DESCSET,
     $                   INFOG2L, PSLARF, PSLARFG, PSELSET,
     $                   PXERBLA, SGEBR2D, SGEBS2D, SLARFG
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD, REAL
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
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
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            MPA0 = NUMROC( M+IROFFA, DESCA( MB_ ), MYROW, IAROW, NPROW )
            NQA0 = NUMROC( N+ICOFFA, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
            LWMIN = MAX( MPA0, NQA0 )
*
            WORK( 1 ) = REAL( LWMIN )
            LQUERY = ( LWORK.EQ.-1 )
            IF( IROFFA.NE.ICOFFA ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(600+NB_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -12
            END IF
         END IF
      END IF
*
      IF( INFO.LT.0 ) THEN
         CALL PXERBLA( ICTXT, 'PSGEBD2', -INFO )
         CALL BLACS_ABORT( ICTXT, 1 )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, II, JJ,
     $              IAROW, IACOL )
*
      IF( M.EQ.1 .AND. N.EQ.1 ) THEN
         IF( MYCOL.EQ.IACOL ) THEN
            IF( MYROW.EQ.IAROW ) THEN
               I = II+(JJ-1)*DESCA( LLD_ )
               CALL SLARFG( 1, A( I ), A( I ), 1, TAUQ( JJ ) )
               D( JJ ) = A( I )
               CALL SGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1, D( JJ ),
     $                       1 )
               CALL SGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1, TAUQ( JJ ),
     $                       1 )
            ELSE
               CALL SGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1, D( JJ ),
     $                       1, IAROW, IACOL )
               CALL SGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1, TAUQ( JJ ),
     $                       1, IAROW, IACOL )
            END IF
         END IF
         IF( MYROW.EQ.IAROW )
     $      TAUP( II ) = ZERO
         RETURN
      END IF
*
      ALPHA = ZERO
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
         DO 10 K = 1, N
            I = IA + K - 1
            J = JA + K - 1
*
*           Generate elementary reflector H(j) to annihilate
*           A(ia+i:ia+m-1,j)
*
            CALL PSLARFG( M-K+1, ALPHA, I, J, A, MIN( I+1, M+IA-1 ),
     $                    J, DESCA, 1, TAUQ )
            CALL PSELSET( D, 1, J, DESCD, ALPHA )
            CALL PSELSET( A, I, J, DESCA, ONE )
*
*           Apply H(i) to A(i:ia+m-1,i+1:ja+n-1) from the left
*
            CALL PSLARF( 'Left', M-K+1, N-K, A, I, J, DESCA, 1, TAUQ, A,
     $                   I, J+1, DESCA, WORK )
            CALL PSELSET( A, I, J, DESCA, ALPHA )
*
            IF( K.LT.N ) THEN
*
*              Generate elementary reflector G(i) to annihilate
*              A(i,ja+j+1:ja+n-1)
*
               CALL PSLARFG( N-K, ALPHA, I, J+1, A, I,
     $                       MIN( J+2, JA+N-1 ), DESCA, DESCA( M_ ),
     $                       TAUP )
               CALL PSELSET( E, I, 1, DESCE, ALPHA )
               CALL PSELSET( A, I, J+1, DESCA, ONE )
*
*              Apply G(i) to A(i+1:ia+m-1,i+1:ja+n-1) from the right
*
               CALL PSLARF( 'Right', M-K, N-K, A, I, J+1, DESCA,
     $                      DESCA( M_ ), TAUP, A, I+1, J+1, DESCA,
     $                      WORK )
               CALL PSELSET( A, I, J+1, DESCA, ALPHA )
            ELSE
               CALL PSELSET( TAUP, I, 1, DESCE, ZERO )
            END IF
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
         DO 20 K = 1, M
            I = IA + K - 1
            J = JA + K - 1
*
*           Generate elementary reflector G(i) to annihilate
*           A(i,ja+j:ja+n-1)
*
            CALL PSLARFG( N-K+1, ALPHA, I, J, A, I,
     $                    MIN( J+1, JA+N-1 ), DESCA, DESCA( M_ ), TAUP )
            CALL PSELSET( D, I, 1, DESCD, ALPHA )
            CALL PSELSET( A, I, J, DESCA, ONE )
*
*           Apply G(i) to A(i:ia+m-1,j:ja+n-1) from the right
*
            CALL PSLARF( 'Right', M-K, N-K+1, A, I, J, DESCA,
     $                   DESCA( M_ ), TAUP, A, MIN( I+1, IA+M-1 ), J,
     $                   DESCA, WORK )
            CALL PSELSET( A, I, J, DESCA, ALPHA )
*
            IF( K.LT.M ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:ia+m-1,j)
*
               CALL PSLARFG( M-K, ALPHA, I+1, J, A,
     $                       MIN( I+2, IA+M-1 ), J, DESCA, 1, TAUQ )
               CALL PSELSET( E, 1, J, DESCE, ALPHA )
               CALL PSELSET( A, I+1, J, DESCA, ONE )
*
*              Apply H(i) to A(i+1:ia+m-1,j+1:ja+n-1) from the left
*
               CALL PSLARF( 'Left', M-K, N-K, A, I+1, J, DESCA, 1, TAUQ,
     $                      A, I+1, J+1, DESCA, WORK )
               CALL PSELSET( A, I+1, J, DESCA, ALPHA )
            ELSE
               CALL PSELSET( TAUQ, 1, J, DESCE, ZERO )
            END IF
   20    CONTINUE
      END IF
*
      WORK( 1 ) = REAL( LWMIN )
*
      RETURN
*
*     End of PSGEBD2
*
      END
