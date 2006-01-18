      SUBROUTINE PDORGLQ( M, N, K, A, IA, JA, DESCA, TAU, WORK, LWORK,
     $                    INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, K, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDORGLQ generates an M-by-N real distributed matrix Q denoting
*  A(IA:IA+M-1,JA:JA+N-1) with orthonormal rows, which is defined as
*  the first M rows of a product of K elementary reflectors of order N
*
*        Q  =  H(k) . . . H(2) H(1)
*
*  as returned by PDGELQF.
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
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix Q. M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix Q. N >= M >= 0.
*
*  K       (global input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. M >= K >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, the i-th row must contain the vector which defines
*          the elementary reflector H(i), IA <= i <= IA+K-1, as
*          returned by PDGELQF in the K rows of its distributed matrix
*          argument A(IA:IA+K-1,JA:*). On exit, this array contains the
*          local pieces of the M-by-N distributed matrix Q.
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
*  TAU     (local input) DOUBLE PRECISION array, dimension LOCr(IA+K-1).
*          This array contains the scalar factors TAU(i) of the
*          elementary reflectors H(i) as returned by PDGELQF.
*          TAU is tied to the distributed matrix A.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                      dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= MB_A * ( MpA0 + NqA0 + MB_A ), where
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          MpA0 = NUMROC( M+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*          NqA0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),
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
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            I, IACOL, IAROW, IB, ICTXT, IINFO, IL, IN, IPW,
     $                   J, LWMIN, MPA0, MYCOL, MYROW, NPCOL, NPROW,
     $                   NQA0
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 2 ), IDUM2( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK1MAT, PDLARFB,
     $                   PDLARFT, PDLASET, PDORGL2, PB_TOPGET,
     $                   PB_TOPSET, PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, INDXG2P, NUMROC
      EXTERNAL           ICEIL, INDXG2P, NUMROC
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
         INFO = -(700+CTXT_)
      ELSE
         CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 7, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            MPA0 = NUMROC( M+MOD( IA-1, DESCA( MB_ ) ), DESCA( MB_ ),
     $                     MYROW, IAROW, NPROW )
            NQA0 = NUMROC( N+MOD( JA-1, DESCA( NB_ ) ), DESCA( NB_ ),
     $                     MYCOL, IACOL, NPCOL )
            LWMIN = DESCA( MB_ ) * ( MPA0 + NQA0 + DESCA( MB_ ) )
*
            WORK( 1 ) = DBLE( LWMIN )
            LQUERY = ( LWORK.EQ.-1 )
            IF( N.LT.M ) THEN
               INFO = -2
            ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
               INFO = -3
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -10
            END IF
         END IF
         IDUM1( 1 ) = K
         IDUM2( 1 ) = 3
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 10
         CALL PCHK1MAT( M, 1, N, 2, IA, JA, DESCA, 7, 2, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDORGLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 )
     $   RETURN
*
      IPW = DESCA( MB_ ) * DESCA( MB_ ) + 1
      IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+K-1 )
      IL = MAX( ( (IA+K-2) / DESCA( MB_ ) ) * DESCA( MB_ ) + 1, IA )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ' ' )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', 'D-ring' )
*
      CALL PDLASET( 'All', IA+M-IL, IL-IA, ZERO, ZERO, A, IL, JA,
     $              DESCA )
*
*     Use unblocked code for the last or only block.
*
      CALL PDORGL2( IA+M-IL, N-IL+IA, IA+K-IL, A, IL, JA+IL-IA, DESCA,
     $              TAU, WORK, LWORK, IINFO )
*
*     Is there at least one block of rows to loop over ?
*
      IF( IL.GT.IN+1 ) THEN
*
*        Use blocked code
*
         DO 10 I = IL-DESCA( MB_ ), IN+1, -DESCA( MB_ )
            IB = MIN( DESCA( MB_ ), IA+M-I )
            J = JA + I - IA
*
            IF( I+IB.LE.IA+M-1 ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL PDLARFT( 'Forward', 'Rowwise', N-I+IA, IB, A, I, J,
     $                       DESCA, TAU, WORK, WORK( IPW ) )
*
*              Apply H' to A(i+ib:ia+m-1,j:ja+n-1) from the right
*
               CALL PDLARFB( 'Right', 'Transpose', 'Forward', 'Rowwise',
     $                       M-I-IB+IA, N-I+IA, IB, A, I, J, DESCA,
     $                       WORK, A, I+IB, J, DESCA, WORK( IPW ) )
            END IF
*
*           Apply H' to columns j:ja+n-1 of current block
*
            CALL PDORGL2( IB, N-I+IA, IB, A, I, J, DESCA, TAU, WORK,
     $                    LWORK, IINFO )
*
*           Set columns ia:i-1 of current block to zero
*
            CALL PDLASET( 'All', IB, I-IA, ZERO, ZERO, A, I, JA, DESCA )
   10    CONTINUE
*
      END IF
*
*     Handle first block separately
*
      IF( IL.GT.IA ) THEN
*
         IB = IN - IA + 1
*
*        Form the triangular factor of the block reflector
*        H = H(i) H(i+1) . . . H(i+ib-1)
*
         CALL PDLARFT( 'Forward', 'Rowwise', N, IB, A, IA, JA, DESCA,
     $                 TAU, WORK, WORK( IPW ) )
*
*        Apply H' to A(ia+ib:ia+m-1,ja:ja+n-1) from the right
*
         CALL PDLARFB( 'Right', 'Transpose', 'Forward', 'Rowwise', M-IB,
     $                 N, IB, A, IA, JA, DESCA, WORK, A, IA+IB, JA,
     $                 DESCA, WORK( IPW ) )
*
*        Apply H' to columns ja:ja+n-1 of current block
*
         CALL PDORGL2( IB, N, IB, A, IA, JA, DESCA, TAU, WORK, LWORK,
     $                 IINFO )
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
*     End of PDORGLQ
*
      END
