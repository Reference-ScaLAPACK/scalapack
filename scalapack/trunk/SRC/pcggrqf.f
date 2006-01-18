      SUBROUTINE PCGGRQF( M, P, N, A, IA, JA, DESCA, TAUA, B, IB, JB,
     $                    DESCB, TAUB, WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, IB, INFO, JA, JB, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX            A( * ), B( * ), TAUA( * ), TAUB( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PCGGRQF computes a generalized RQ factorization of
*  an M-by-N matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1)
*  and a P-by-N matrix sub( B ) = B(IB:IB+P-1,JB:JB+N-1):
*
*              sub( A ) = R*Q,        sub( B ) = Z*T*Q,
*
*  where Q is an N-by-N unitary matrix, Z is a P-by-P unitary matrix,
*  and R and T assume one of the forms:
*
*  if M <= N,  R = ( 0  R12 ) M,   or if M > N,  R = ( R11 ) M-N,
*                   N-M  M                           ( R21 ) N
*                                                       N
*
*  where R12 or R21 is upper triangular, and
*
*  if P >= N,  T = ( T11 ) N  ,   or if P < N,  T = ( T11  T12 ) P,
*                  (  0  ) P-N                         P   N-P
*                     N
*
*  where T11 is upper triangular.
*
*  In particular, if sub( B ) is square and nonsingular, the GRQ
*  factorization of sub( A ) and sub( B ) implicitly gives the RQ
*  factorization of sub( A )*inv( sub( B ) ):
*
*               sub( A )*inv( sub( B ) ) = (R*inv(T))*Z'
*
*  where inv( sub( B ) ) denotes the inverse of the matrix sub( B ),
*  and Z' denotes the conjugate transpose of matrix Z.
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
*          The number of rows to be operated on i.e the number of
*          rows of the distributed submatrix sub( A ).  M >= 0.
*
*  P       (global input) INTEGER
*          The number of rows to be operated on i.e the number of
*          rows of the distributed submatrix sub( B ).  P >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrices sub( A ) and sub( B ).
*          N >= 0.
*
*  A       (local input/local output) COMPLEX pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, the local pieces of the M-by-N distributed matrix
*          sub( A ) which is to be factored. On exit, if M <= N, the
*          upper triangle of A( IA:IA+M-1, JA+N-M:JA+N-1 ) contains the
*          M by M upper triangular matrix R; if M >= N, the elements on
*          and above the (M-N)-th subdiagonal contain the M by N upper
*          trapezoidal matrix R; the remaining elements, with the array
*          TAUA, represent the unitary matrix Q as a product of
*          elementary reflectors (see Further Details).
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
*  TAUA    (local output) COMPLEX, array, dimension LOCr(IA+M-1)
*          This array contains the scalar factors of the elementary
*          reflectors which represent the unitary matrix Q. TAUA is
*          tied to the distributed matrix A (see Further Details).
*
*  B       (local input/local output) COMPLEX pointer into the
*          local memory to an array of dimension (LLD_B, LOCc(JB+N-1)).
*          On entry, the local pieces of the P-by-N distributed matrix
*          sub( B ) which is to be factored.  On exit, the elements on
*          and above the diagonal of sub( B ) contain the min(P,N) by N
*          upper trapezoidal matrix T (T is upper triangular if P >= N);
*          the elements below the diagonal, with the array TAUB,
*          represent the unitary matrix Z as a product of elementary
*          reflectors (see Further Details).
*
*  IB      (global input) INTEGER
*          The row index in the global array B indicating the first
*          row of sub( B ).
*
*  JB      (global input) INTEGER
*          The column index in the global array B indicating the
*          first column of sub( B ).
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
*
*  TAUB    (local output) COMPLEX, array, dimension
*          LOCc(JB+MIN(P,N)-1). This array contains the scalar factors
*          TAUB of the elementary reflectors which represent the unitary
*          matrix Z. TAUB is tied to the distributed matrix B (see
*          Further Details).
*
*  WORK    (local workspace/local output) COMPLEX array,
*                                                   dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= MAX( MB_A * ( MpA0 + NqA0 + MB_A ),
*                        MAX( (MB_A*(MB_A-1))/2, (PpB0 + NqB0)*MB_A ) +
*                             MB_A * MB_A,
*                        NB_B * ( PpB0 + NqB0 + NB_B ) ), where
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW  = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL  = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          MpA0   = NUMROC( M+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*          NqA0   = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),
*
*          IROFFB = MOD( IB-1, MB_B ), ICOFFB = MOD( JB-1, NB_B ),
*          IBROW  = INDXG2P( IB, MB_B, MYROW, RSRC_B, NPROW ),
*          IBCOL  = INDXG2P( JB, NB_B, MYCOL, CSRC_B, NPCOL ),
*          PpB0   = NUMROC( P+IROFFB, MB_B, MYROW, IBROW, NPROW ),
*          NqB0   = NUMROC( N+ICOFFB, NB_B, MYCOL, IBCOL, NPCOL ),
*
*          and NUMROC, INDXG2P are ScaLAPACK tool functions;
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
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(ia)' H(ia+1)' . . . H(ia+k-1)', where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - taua * v * v'
*
*  where taua is a complex scalar, and v is a complex vector with
*  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; conjg(v(1:n-k+i-1)) is stored on
*  exit in A(ia+m-k+i-1,ja:ja+n-k+i-2), and taua in TAUA(ia+m-k+i-1).
*  To form Q explicitly, use ScaLAPACK subroutine PCUNGRQ.
*  To use Q to update another matrix, use ScaLAPACK subroutine PCUNMRQ.
*
*  The matrix Z is represented as a product of elementary reflectors
*
*     Z = H(jb) H(jb+1) . . . H(jb+k-1), where k = min(p,n).
*
*  Each H(i) has the form
*
*     H(i) = I - taub * v * v'
*
*  where taub is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:p) is stored on exit in
*  B(ib+i:ib+p-1,jb+i-1), and taub in TAUB(jb+i-1).
*  To form Z explicitly, use ScaLAPACK subroutine PCUNGQR.
*  To use Z to update another matrix, use ScaLAPACK subroutine PCUNMQR.
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrices sub( A ) and sub( B ) must verify some
*  alignment properties, namely the following expression should be true:
*
*  ( NB_A.EQ.NB_B .AND. ICOFFA.EQ.ICOFFB .AND. IACOL.EQ.IBCOL )
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            IACOL, IAROW, IBCOL, IBROW, ICOFFA, ICOFFB,
     $                   ICTXT, IROFFA, IROFFB, LWMIN, MPA0, MYCOL,
     $                   MYROW, NPCOL, NPROW, NQA0, NQB0, PPB0
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCGEQRF, PCGERQF,
     $                   PCHK2MAT, PCUNMRQ, PXERBLA
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, INT, MAX, MIN, MOD, REAL
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
         INFO = -707
      ELSE
         CALL CHK1MAT( M, 1, N, 3, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( P, 2, N, 3, IB, JB, DESCB, 12, INFO )
         IF( INFO.EQ.0 ) THEN
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            ICOFFB = MOD( JB-1, DESCB( NB_ ) )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $                       NPROW )
            IBCOL = INDXG2P( JB, DESCB( NB_ ), MYCOL, DESCB( CSRC_ ),
     $                       NPCOL )
            MPA0 = NUMROC( M+IROFFA, DESCA( MB_ ), MYROW, IAROW, NPROW )
            NQA0 = NUMROC( N+ICOFFA, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
            PPB0 = NUMROC( P+IROFFB, DESCB( MB_ ), MYROW, IBROW, NPROW )
            NQB0 = NUMROC( N+ICOFFB, DESCB( NB_ ), MYCOL, IBCOL, NPCOL )
            LWMIN = MAX( DESCA( MB_ ) * ( MPA0 + NQA0 + DESCA( MB_ ) ),
     $        MAX( MAX( ( DESCA( MB_ )*( DESCA( MB_ ) - 1 ) ) / 2,
     $        ( PPB0 + NQB0 ) * DESCA( MB_ ) ) +
     $          DESCA( MB_ ) * DESCA( MB_ ),
     $        DESCB( NB_ ) * ( PPB0 + NQB0 + DESCB( NB_ ) ) ) )
*
            WORK( 1 ) = CMPLX( REAL( LWMIN ) )
            LQUERY = ( LWORK.EQ.-1 )
            IF( IACOL.NE.IBCOL .OR. ICOFFA.NE.ICOFFB ) THEN
               INFO = -11
            ELSE IF( DESCA( NB_ ).NE.DESCB( NB_ ) ) THEN
               INFO = -1204
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -1207
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -15
            END IF
         END IF
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 1 ) = -1
         ELSE
            IDUM1( 1 ) = 1
         END IF
         IDUM2( 1 ) = 15
         CALL PCHK2MAT( M, 1, N, 3, IA, JA, DESCA, 7, P, 2, N, 3, IB,
     $                  JB, DESCB, 12, 1, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCGGRQF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     RQ factorization of M-by-N matrix sub( A ): sub( A ) = R*Q
*
      CALL PCGERQF( M, N, A, IA, JA, DESCA, TAUA, WORK, LWORK, INFO )
      LWMIN = INT( WORK( 1 ) )
*
*     Update sub( B ) := sub( B )*Q'
*
      CALL PCUNMRQ( 'Right', 'Conjugate Transpose', P, N, MIN( M, N ),
     $              A, MAX( IA, IA+M-N ), JA, DESCA, TAUA, B, IB, JB,
     $              DESCB, WORK, LWORK, INFO )
      LWMIN = MAX( LWMIN, INT( WORK( 1 ) ) )
*
*     QR factorization of P-by-N matrix sub( B ): sub( B ) = Z*T
*
      CALL PCGEQRF( P, N, B, IB, JB, DESCB, TAUB, WORK, LWORK, INFO )
      WORK( 1 ) = CMPLX( REAL( MAX( LWMIN, INT( WORK( 1 ) ) ) ) )
*
      RETURN
*
*     End of PCGGRQF
*
      END
