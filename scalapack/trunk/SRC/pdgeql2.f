      SUBROUTINE PDGEQL2( M, N, A, IA, JA, DESCA, TAU, WORK, LWORK,
     $                    INFO )
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
      DOUBLE PRECISION   A( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDGEQL2 computes a QL factorization of a real distributed M-by-N
*  matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) = Q * L.
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
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, the local pieces of the M-by-N distributed matrix
*          sub( A ) which is to be factored. On exit, if M >= N, the
*          lower triangle of the distributed submatrix
*          A( IA+M-N:IA+M-1, JA:JA+N-1 ) contains the N-by-N lower
*          triangular matrix L; if M <= N, the elements on and below
*          the (N-M)-th superdiagonal contain the M by N lower
*          trapezoidal matrix L; the remaining elements, with the
*          array TAU, represent the orthogonal matrix Q as a product of
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
*  TAU     (local output) DOUBLE PRECISION array, dimension LOCc(JA+N-1)
*          This array contains the scalar factors of the elementary
*          reflectors. TAU is tied to the distributed matrix A.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                   dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= Mp0 + MAX( 1, Nq0 ), where
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
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(ja+k-1) . . . H(ja+1) H(ja), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in
*  A(ia:ia+m-k+i-2,ja+n-k+i-1), and tau in TAU(ja+n-k+i-1).
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
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            I, IACOL, IAROW, ICTXT, II, J, JJ, K, LWMIN,
     $                   MP, MYCOL, MYROW, NPCOL, NPROW, NQ
      DOUBLE PRECISION   AJJ, ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GRIDINFO, CHK1MAT, DGEBR2D,
     $                   DGEBS2D, DLARFG, DSCAL, INFOG2L,
     $                   PDELSET, PDLARF, PDLARFG, PB_TOPGET,
     $                   PB_TOPSET, PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, MOD
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
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            MP = NUMROC( M+MOD( IA-1, DESCA( MB_ ) ), DESCA( MB_ ),
     $                   MYROW, IAROW, NPROW )
            NQ = NUMROC( N+MOD( JA-1, DESCA( NB_ ) ), DESCA( NB_ ),
     $                   MYCOL, IACOL, NPCOL )
            LWMIN = MP + MAX( 1, NQ )
*
            WORK( 1 ) = DBLE( LWMIN )
            LQUERY = ( LWORK.EQ.-1 )
            IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY )
     $         INFO = -9
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDGEQL2', -INFO )
         CALL BLACS_ABORT( ICTXT, 1 )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'D-ring' )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
*
      IF( DESCA( M_ ).EQ.1 ) THEN
         IF( MYCOL.EQ.IACOL )
     $      NQ = NQ - MOD( JA-1, DESCA( NB_ ) )
         CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, II,
     $                 JJ, IAROW, IACOL )
         IACOL = INDXG2P( JA+N-1, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                    NPCOL )
         IF( MYROW.EQ.IAROW ) THEN
            IF( MYCOL.EQ.IACOL ) THEN
               I = II+(JJ+NQ-2)*DESCA( LLD_ )
               AJJ = A( I )
               CALL DLARFG( 1, AJJ, A( I ), 1, TAU( JJ+NQ-1 ) )
               IF( N.GT.1 ) THEN
                  ALPHA = ONE - TAU( JJ+NQ-1 )
                  CALL DGEBS2D( ICTXT, 'Rowwise', ' ', 1, 1, ALPHA, 1 )
                  CALL DSCAL( NQ-1, ALPHA, A( II+(JJ-1)*DESCA( LLD_ ) ),
     $                        DESCA( LLD_ ) )
               END IF
               CALL DGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1,
     $                       TAU( JJ+NQ-1 ), 1 )
               A( I ) = AJJ
            ELSE
               IF( N.GT.1 ) THEN
                  CALL DGEBR2D( ICTXT, 'Rowwise', ' ', 1, 1, ALPHA,
     $                          1, IAROW, IACOL )
                  CALL DSCAL( NQ, ALPHA, A( II+(JJ-1)*DESCA( LLD_ ) ),
     $                        DESCA( LLD_ ) )
               END IF
            END IF
         ELSE IF( MYCOL.EQ.IACOL ) THEN
            CALL DGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1,
     $                         TAU( JJ+NQ-1 ), 1, IAROW, IACOL )
         END IF
*
      ELSE
*
         K = MIN( M, N )
         DO 10 J = JA+K-1, JA, -1
            I = IA + J - JA
*
*           Generate elementary reflector H(j) to annihilate
*           A(ia:i+m-k-1,j+n-k)
*
            CALL PDLARFG( M-K+I-IA+1, AJJ, M-K+I, N-K+J, A, IA,
     $                    N-K+J, DESCA, 1, TAU )
*
*           Apply H(j) to A(ia:i+m-k,ja:j+n-k-1) from the left
*
            CALL PDELSET( A, I+M-K, J+N-K, DESCA, ONE )
            CALL PDLARF( 'Left', M-K+I-IA+1, N-K+J-JA, A, IA, N-K+J,
     $                   DESCA, 1, TAU, A, IA, JA, DESCA, WORK )
            CALL PDELSET( A, I+M-K, J+N-K, DESCA, AJJ )
*
   10    CONTINUE
*
      END IF
*
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      WORK( 1 ) = DBLE( LWMIN )
*
      RETURN
*
*     End of PDGEQL2
*
      END
