      SUBROUTINE PZGGQRF( N, M, P, A, IA, JA, DESCA, TAUA, B, IB, JB,
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
      COMPLEX*16         A( * ), B( * ), TAUA( * ), TAUB( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZGGQRF computes a generalized QR factorization of
*  an N-by-M matrix sub( A ) = A(IA:IA+N-1,JA:JA+M-1) and
*  an N-by-P matrix sub( B ) = B(IB:IB+N-1,JB:JB+P-1):
*
*              sub( A ) = Q*R,        sub( B ) = Q*T*Z,
*
*  where Q is an N-by-N unitary matrix, Z is a P-by-P unitary matrix,
*  and R and T assume one of the forms:
*
*  if N >= M,  R = ( R11 ) M  ,   or if N < M,  R = ( R11  R12 ) N,
*                  (  0  ) N-M                         N   M-N
*                     M
*
*  where R11 is upper triangular, and
*
*  if N <= P,  T = ( 0  T12 ) N,   or if N > P,  T = ( T11 ) N-P,
*                   P-N  N                           ( T21 ) P
*                                                       P
*
*  where T12 or T21 is upper triangular.
*
*  In particular, if sub( B ) is square and nonsingular, the GQR
*  factorization of sub( A ) and sub( B ) implicitly gives the QR
*  factorization of inv( sub( B ) )* sub( A ):
*
*               inv( sub( B ) )*sub( A )= Z'*(inv(T)*R)
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
*  N       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrices sub( A ) and sub( B ). N >= 0.
*
*  M       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( A ).  M >= 0.
*
*  P       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( B ).  P >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+M-1)).
*          On entry, the local pieces of the N-by-M distributed matrix
*          sub( A ) which is to be factored.  On exit, the elements on
*          and above the diagonal of sub( A ) contain the min(N,M) by M
*          upper trapezoidal matrix R (R is upper triangular if N >= M);
*          the elements below the diagonal, with the array TAUA,
*          represent the unitary matrix Q as a product of min(N,M)
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
*  TAUA    (local output) COMPLEX*16, array, dimension
*          LOCc(JA+MIN(N,M)-1). This array contains the scalar factors
*          TAUA of the elementary reflectors which represent the unitary
*          matrix Q. TAUA is tied to the distributed matrix A. (see
*          Further Details).
*
*  B       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_B, LOCc(JB+P-1)).
*          On entry, the local pieces of the N-by-P distributed matrix
*          sub( B ) which is to be factored. On exit, if N <= P, the
*          upper triangle of B(IB:IB+N-1,JB+P-N:JB+P-1) contains the
*          N by N upper triangular matrix T; if N > P, the elements on
*          and above the (N-P)-th subdiagonal contain the N by P upper
*          trapezoidal matrix T; the remaining elements, with the array
*          TAUB, represent the unitary matrix Z as a product of
*          elementary reflectors (see Further Details).
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
*  TAUB    (local output) COMPLEX*16, array, dimension LOCr(IB+N-1)
*          This array contains the scalar factors of the elementary
*          reflectors which represent the unitary matrix Z. TAUB is
*          tied to the distributed matrix B (see Further Details).
*
*  WORK    (local workspace/local output) COMPLEX*16 array,
*                                                  dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= MAX( NB_A * ( NpA0 + MqA0 + NB_A ),
*                        MAX( (NB_A*(NB_A-1))/2, (PqB0 + NpB0)*NB_A ) +
*                             NB_A * NB_A,
*                        MB_B * ( NpB0 + PqB0 + MB_B ) ), where
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW  = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL  = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          NpA0   = NUMROC( N+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*          MqA0   = NUMROC( M+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),
*
*          IROFFB = MOD( IB-1, MB_B ), ICOFFB = MOD( JB-1, NB_B ),
*          IBROW  = INDXG2P( IB, MB_B, MYROW, RSRC_B, NPROW ),
*          IBCOL  = INDXG2P( JB, NB_B, MYCOL, CSRC_B, NPCOL ),
*          NpB0   = NUMROC( N+IROFFB, MB_B, MYROW, IBROW, NPROW ),
*          PqB0   = NUMROC( P+ICOFFB, NB_B, MYCOL, IBCOL, NPCOL ),
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
*     Q = H(ja) H(ja+1) . . . H(ja+k-1), where k = min(n,m).
*
*  Each H(i) has the form
*
*     H(i) = I - taua * v * v'
*
*  where taua is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in
*  A(ia+i:ia+n-1,ja+i-1), and taua in TAUA(ja+i-1).
*  To form Q explicitly, use ScaLAPACK subroutine PZUNGQR.
*  To use Q to update another matrix, use ScaLAPACK subroutine PZUNMQR.
*
*  The matrix Z is represented as a product of elementary reflectors
*
*     Z = H(ib)' H(ib+1)' . . . H(ib+k-1)', where k = min(n,p).
*
*  Each H(i) has the form
*
*     H(i) = I - taub * v * v'
*
*  where taub is a complex scalar, and v is a complex vector with
*  v(p-k+i+1:p) = 0 and v(p-k+i) = 1; conjg(v(1:p-k+i-1)) is stored on
*  exit in B(ib+n-k+i-1,jb:jb+p-k+i-2), and taub in TAUB(ib+n-k+i-1).
*  To form Z explicitly, use ScaLAPACK subroutine PZUNGRQ.
*  To use Z to update another matrix, use ScaLAPACK subroutine PZUNMRQ.
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrices sub( A ) and sub( B ) must verify some
*  alignment properties, namely the following expression should be true:
*
*  ( MB_A.EQ.MB_B .AND. IROFFA.EQ.IROFFB .AND. IAROW.EQ.IBROW )
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            IACOL, IAROW, IBCOL, IBROW, ICOFFA, ICOFFB,
     $                   ICTXT, IROFFA, IROFFB, LWMIN, MQA0, MYCOL,
     $                   MYROW, NPA0, NPB0, NPCOL, NPROW, PQB0
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK2MAT, PXERBLA,
     $                   PZGEQRF, PZGERQF, PZUNMQR
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, INT, MAX, MIN, MOD
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
         CALL CHK1MAT( N, 1, M, 2, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 1, P, 3, IB, JB, DESCB, 12, INFO )
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
            NPA0 = NUMROC( N+IROFFA, DESCA( MB_ ), MYROW, IAROW, NPROW )
            MQA0 = NUMROC( M+ICOFFA, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
            NPB0 = NUMROC( N+IROFFB, DESCB( MB_ ), MYROW, IBROW, NPROW )
            PQB0 = NUMROC( P+ICOFFB, DESCB( NB_ ), MYCOL, IBCOL, NPCOL )
            LWMIN = MAX( DESCA( NB_ ) * ( NPA0 + MQA0 + DESCA( NB_ ) ),
     $        MAX( MAX( ( DESCA( NB_ )*( DESCA( NB_ ) - 1 ) ) / 2,
     $         ( PQB0 + NPB0 ) * DESCA( NB_ ) ) +
     $           DESCA( NB_ ) * DESCA( NB_ ),
     $         DESCB( MB_ ) * ( NPB0 + PQB0 + DESCB( MB_ ) ) ) )
*
            WORK( 1 ) = DCMPLX( DBLE( LWMIN ) )
            LQUERY = ( LWORK.EQ.-1 )
            IF( IAROW.NE.IBROW .OR. IROFFA.NE.IROFFB ) THEN
               INFO = -10
            ELSE IF( DESCA( MB_ ).NE.DESCB( MB_ ) ) THEN
               INFO = -1203
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -1207
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -15
            END IF
         END IF
         IF( LQUERY ) THEN
            IDUM1( 1 ) = -1
         ELSE
            IDUM1( 1 ) = 1
         END IF
         IDUM2( 1 ) = 15
         CALL PCHK2MAT( N, 1, M, 2, IA, JA, DESCA, 7, N, 1, P, 3, IB,
     $                  JB, DESCB, 12, 1, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZGGQRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     QR factorization of N-by-M matrix sub( A ): sub( A ) = Q*R
*
      CALL PZGEQRF( N, M, A, IA, JA, DESCA, TAUA, WORK, LWORK, INFO )
      LWMIN = INT( WORK( 1 ) )
*
*     Update sub( B ) := Q'*sub( B ).
*
      CALL PZUNMQR( 'Left', 'Conjugate Transpose', N, P, MIN( N, M ), A,
     $              IA, JA, DESCA, TAUA, B, IB, JB, DESCB, WORK, LWORK,
     $              INFO )
      LWMIN = MIN( LWMIN, INT( WORK( 1 ) ) )
*
*     RQ factorization of N-by-P matrix sub( B ): sub( B ) = T*Z.
*
      CALL PZGERQF( N, P, B, IB, JB, DESCB, TAUB, WORK, LWORK, INFO )
      WORK( 1 ) = DCMPLX( DBLE( MAX( LWMIN, INT( WORK( 1 ) ) ) ) )
*
      RETURN
*
*     End of PZGGQRF
*
      END
