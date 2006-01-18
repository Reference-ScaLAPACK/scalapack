      SUBROUTINE PSORMR3( SIDE, TRANS, M, N, K, L, A, IA, JA, DESCA,
     $                    TAU, C, IC, JC, DESCC, WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            IA, IC, INFO, JA, JC, K, L, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCC( * )
      REAL               A( * ), C( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSORMR3 overwrites the general real M-by-N distributed matrix
*  sub( C ) = C(IC:IC+M-1,JC:JC+N-1) with
*
*                       SIDE = 'L'          SIDE = 'R'
*  TRANS = 'N':      Q * sub( C )         sub(  C ) * Q
*  TRANS = 'T':      Q**T * sub( C )      sub( C ) * Q**T
*
*  where Q is a real orthogonal distributed matrix defined as the
*  product of K elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by PSTZRZF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
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
*  SIDE    (global input) CHARACTER
*          = 'L': apply Q or Q**T from the Left;
*          = 'R': apply Q or Q**T from the Right.
*
*  TRANS   (global input) CHARACTER
*          = 'N':  No transpose, apply Q;
*          = 'T':  Transpose, apply Q**T.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( C ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( C ). N >= 0.
*
*  K       (global input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q.  If SIDE = 'L', M >= K >= 0, if SIDE = 'R',
*          N >= K >= 0.
*
*  L       (global input) INTEGER
*          The columns of the distributed submatrix sub( A ) containing
*          the meaningful part of the Householder reflectors.
*          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
*
*  A       (local input) REAL pointer into the local memory
*          to an array of dimension (LLD_A,LOCc(JA+M-1)) if SIDE='L',
*          and (LLD_A,LOCc(JA+N-1)) if SIDE='R', where
*          LLD_A >= MAX(1,LOCr(IA+K-1)); On entry, the i-th row must
*          contain the vector which defines the elementary reflector
*          H(i), IA <= i <= IA+K-1, as returned by PSTZRZF in the
*          K rows of its distributed matrix argument A(IA:IA+K-1,JA:*).
*          A(IA:IA+K-1,JA:*) is modified by the routine but restored on
*          exit.
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
*  TAU     (local input) REAL, array, dimension LOCc(IA+K-1).
*          This array contains the scalar factors TAU(i) of the
*          elementary reflectors H(i) as returned by PSTZRZF.
*          TAU is tied to the distributed matrix A.
*
*  C       (local input/local output) REAL pointer into the
*          local memory to an array of dimension (LLD_C,LOCc(JC+N-1)).
*          On entry, the local pieces of the distributed matrix sub(C).
*          On exit, sub( C ) is overwritten by Q*sub( C ) or Q'*sub( C )
*          or sub( C )*Q' or sub( C )*Q.
*
*  IC      (global input) INTEGER
*          The row index in the global array C indicating the first
*          row of sub( C ).
*
*  JC      (global input) INTEGER
*          The column index in the global array C indicating the
*          first column of sub( C ).
*
*  DESCC   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix C.
*
*  WORK    (local workspace/local output) REAL array,
*                                                   dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          If SIDE = 'L', LWORK >= MpC0 + MAX( MAX( 1, NqC0 ), NUMROC(
*                  NUMROC( M+IROFFC,MB_A,0,0,NPROW ),MB_A,0,0,LCMP ) );
*          if SIDE = 'R', LWORK >= NqC0 + MAX( 1, MpC0 );
*
*          where LCMP = LCM / NPROW with LCM = ICLM( NPROW, NPCOL ),
*
*          IROFFC = MOD( IC-1, MB_C ), ICOFFC = MOD( JC-1, NB_C ),
*          ICROW = INDXG2P( IC, MB_C, MYROW, RSRC_C, NPROW ),
*          ICCOL = INDXG2P( JC, NB_C, MYCOL, CSRC_C, NPCOL ),
*          MpC0 = NUMROC( M+IROFFC, MB_C, MYROW, ICROW, NPROW ),
*          NqC0 = NUMROC( N+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL ),
*
*          ILCM, INDXG2P and NUMROC are ScaLAPACK tool functions;
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
*  INFO    (local output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrices A(IA:*, JA:*) and C(IC:IC+M-1,JC:JC+N-1)
*  must verify some alignment properties, namely the following
*  expressions should be true:
*
*  If SIDE = 'L',
*    ( NB_A.EQ.MB_C .AND. ICOFFA.EQ.IROFFC )
*  If SIDE = 'R',
*    ( NB_A.EQ.NB_C .AND. ICOFFA.EQ.ICOFFC .AND. IACOL.EQ.ICCOL )
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
      LOGICAL            LEFT, LQUERY, NOTRAN
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            I, I1, I2, I3, IACOL, ICC, ICCOL, ICOFFA,
     $                   ICOFFC, ICROW, ICTXT, IROFFC, JAA, JCC, LCM,
     $                   LCMP, LWMIN, MI, MPC0, MYCOL, MYROW, NI, NPCOL,
     $                   NPROW, NQ, NQC0
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GRIDINFO, CHK1MAT, PSLARZ,
     $                   PB_TOPGET, PB_TOPSET, PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILCM, INDXG2P, NUMROC
      EXTERNAL           ILCM, INDXG2P, LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MOD, REAL
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
         INFO = -(900+CTXT_)
      ELSE
         LEFT = LSAME( SIDE, 'L' )
         NOTRAN = LSAME( TRANS, 'N' )
*
*        NQ is the order of Q
*
         IF( LEFT ) THEN
            NQ = M
            CALL CHK1MAT( K, 5, M, 3, IA, JA, DESCA, 10, INFO )
         ELSE
            NQ = N
            CALL CHK1MAT( K, 5, N, 4, IA, JA, DESCA, 10, INFO )
         END IF
         CALL CHK1MAT( M, 3, N, 4, IC, JC, DESCC, 15, INFO )
         IF( INFO.EQ.0 ) THEN
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFC = MOD( IC-1, DESCC( MB_ ) )
            ICOFFC = MOD( JC-1, DESCC( NB_ ) )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            ICROW = INDXG2P( IC, DESCC( MB_ ), MYROW, DESCC( RSRC_ ),
     $                       NPROW )
            ICCOL = INDXG2P( JC, DESCC( NB_ ), MYCOL, DESCC( CSRC_ ),
     $                       NPCOL )
            MPC0 = NUMROC( M+IROFFC, DESCC( MB_ ), MYROW, ICROW, NPROW )
            NQC0 = NUMROC( N+ICOFFC, DESCC( NB_ ), MYCOL, ICCOL, NPCOL )
*
            IF( LEFT ) THEN
               LCM = ILCM( NPROW, NPCOL )
               LCMP = LCM / NPROW
               LWMIN = MPC0 + MAX( MAX( 1, NQC0 ), NUMROC( NUMROC(
     $                 M+IROFFC, DESCA( MB_ ), 0, 0, NPROW ),
     $                 DESCA( MB_ ), 0, 0, LCMP ) )
            ELSE
               LWMIN = NQC0 + MAX( 1, MPC0 )
            END IF
*
            WORK( 1 ) = REAL( LWMIN )
            LQUERY = ( LWORK.EQ.-1 )
            IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
               INFO = -2
            ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
               INFO = -5
            ELSE IF( L.LT.0 .OR. L.GT.NQ ) THEN
               INFO = -6
            ELSE IF( LEFT .AND. DESCA( NB_ ).NE.DESCC( MB_ ) ) THEN
               INFO = -(1000+NB_)
            ELSE IF( LEFT .AND. ICOFFA.NE.IROFFC ) THEN
               INFO = -13
            ELSE IF( .NOT.LEFT .AND. ICOFFA.NE.ICOFFC ) THEN
               INFO = -14
            ELSE IF( .NOT.LEFT .AND. IACOL.NE.ICCOL ) THEN
               INFO = -14
            ELSE IF( .NOT.LEFT .AND. DESCA( NB_ ).NE.DESCC( NB_ ) ) THEN
               INFO = -(1500+NB_)
            ELSE IF( ICTXT.NE.DESCC( CTXT_ ) ) THEN
               INFO = -(1500+CTXT_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -17
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PSORMR3', -INFO )
         CALL BLACS_ABORT( ICTXT, 1 )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      IF( ( LEFT .AND. .NOT.NOTRAN .OR. .NOT.LEFT .AND. NOTRAN ) ) THEN
         I1 = IA
         I2 = IA + K - 1
         I3 = 1
      ELSE
         I1 = IA + K - 1
         I2 = IA
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
         JCC = JC
         JAA = JA + M - L
      ELSE
         MI = M
         ICC = IC
         JAA = JA + N - L
         CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ' ' )
         IF( NOTRAN ) THEN
            CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', 'I-ring' )
         ELSE
            CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', 'D-ring' )
         END IF
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) or H(i)' is applied to C(ic+i-ia:icc+m-1,jc:jc+n-1)
*
            MI = M - I + IA
            ICC = IC + I - IA
         ELSE
*
*           H(i) or H(i)' is applied to C(ic:ic+m-1,jc+i-ia:jc+n-1)
*
            NI = N - I + IA
            JCC = JC + I - IA
         END IF
*
*        Apply H(i) or H(i)'
*
         CALL PSLARZ( SIDE, MI, NI, L, A, I, JAA, DESCA, DESCA( M_ ),
     $                TAU, C, ICC, JCC, DESCC, WORK )
*
   10 CONTINUE
*
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      WORK( 1 ) = REAL( LWMIN )
*
      RETURN
*
*     End of PSORMR3
*
      END
