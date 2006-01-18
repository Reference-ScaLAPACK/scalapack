      SUBROUTINE PDLATRZ( M, N, L, A, IA, JA, DESCA, TAU, WORK )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, JA, L, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLATRZ reduces the M-by-N ( M<=N ) real upper trapezoidal matrix
*  sub( A ) = [ A(IA:IA+M-1,JA:JA+M-1) A(IA:IA+M-1,JA+N-L:JA+N-1) ] to
*  upper triangular form by means of orthogonal transformations.
*
*  The upper trapezoidal matrix sub( A ) is factored as
*
*     sub( A ) = ( R  0 ) * Z,
*
*  where Z is an N-by-N orthogonal matrix and R is an M-by-M upper
*  triangular matrix.
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
*  L       (global input) INTEGER
*          The columns of the distributed submatrix sub( A ) containing
*          the meaningful part of the Householder reflectors. L > 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, the local pieces of the M-by-N distributed matrix
*          sub( A ) which is to be factored. On exit, the leading M-by-M
*          upper triangular part of sub( A ) contains the upper trian-
*          gular matrix R, and elements N-L+1 to N of the first M rows
*          of sub( A ), with the array TAU, represent the orthogonal
*          matrix Z as a product of M elementary reflectors.
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
*  TAU     (local output) DOUBLE PRECISION array, dimension LOCr(IA+M-1)
*          This array contains the scalar factors of the elementary
*          reflectors. TAU is tied to the distributed matrix A.
*
*  WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK)
*          LWORK >= Nq0 + MAX( 1, Mp0 ), where
*
*          IROFF = MOD( IA-1, MB_A ), ICOFF = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          Mp0   = NUMROC( M+IROFF, MB_A, MYROW, IAROW, NPROW ),
*          Nq0   = NUMROC( N+ICOFF, NB_A, MYCOL, IACOL, NPCOL ),
*
*          and NUMROC, INDXG2P are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*  Further Details
*  ===============
*
*  The  factorization is obtained by Householder's method.  The kth
*  transformation matrix, Z( k ), which is used to introduce zeros into
*  the (m - k + 1)th row of sub( A ), is given in the form
*
*     Z( k ) = ( I     0   ),
*              ( 0  T( k ) )
*
*  where
*
*     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),
*                                                 (   0    )
*                                                 ( z( k ) )
*
*  tau is a scalar and z( k ) is an ( n - m ) element vector.
*  tau and z( k ) are chosen to annihilate the elements of the kth row
*  of sub( A ).
*
*  The scalar tau is returned in the kth element of TAU and the vector
*  u( k ) in the kth row of sub( A ), such that the elements of z( k )
*  are in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned
*  in the upper triangular part of sub( A ).
*
*  Z is given by
*
*     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IAROW, ICTXT, II, J, J1, MP, MYCOL, MYROW,
     $                   NPCOL, NPROW
      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           INFOG1L, PDELSET, PDLARFG, PDLARZ
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      MP = NUMROC( IA+M-1, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $             NPROW )
*
      IF( M.EQ.N ) THEN
*
         CALL INFOG1L( IA, DESCA( MB_ ), NPROW, MYROW, DESCA( RSRC_ ),
     $                 II, IAROW )
         DO 10 I = II, MP
            TAU( I ) = ZERO
   10    CONTINUE
*
      ELSE
*
         J1 = JA + N - L
         DO 20 I = IA+M-1, IA, -1
            J = JA + I - IA
*
*           Generate elementary reflector H(i) to annihilate
*           [ A(i, j) A(i,j1:ja+n-1) ]
*
            CALL PDLARFG( L+1, AII, I, J, A, I, J1, DESCA, DESCA( M_ ),
     $                         TAU )
*
*           Apply H(i) to A(ia:i-1,j:ja+n-1) from the right
*
            CALL PDLARZ( 'Right', I-IA, JA+N-J, L, A, I, J1, DESCA,
     $                   DESCA( M_ ), TAU, A, IA, J, DESCA, WORK )
            CALL PDELSET( A, I, J, DESCA, AII )
*
   20    CONTINUE
*
      END IF
*
      RETURN
*
*     End of PDLATRZ
*
      END
