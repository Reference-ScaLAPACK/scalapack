      SUBROUTINE PDORMTR( SIDE, UPLO, TRANS, M, N, A, IA, JA, DESCA,
     $                    TAU, C, IC, JC, DESCC, WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, UPLO
      INTEGER            IA, IC, INFO, JA, JC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCC( * )
      DOUBLE PRECISION   A( * ), C( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDORMTR overwrites the general real M-by-N distributed matrix
*  sub( C ) = C(IC:IC+M-1,JC:JC+N-1) with
*
*                       SIDE = 'L'           SIDE = 'R'
*  TRANS = 'N':      Q * sub( C )          sub( C ) * Q
*  TRANS = 'T':      Q**T * sub( C )       sub( C ) * Q**T
*
*  where Q is a real orthogonal distributed matrix of order nq, with
*  nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the
*  product of nq-1 elementary reflectors, as returned by PDSYTRD:
*
*  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
*
*  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
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
*  UPLO    (global input) CHARACTER
*          = 'U': Upper triangle of A(IA:*,JA:*) contains elementary
*                 reflectors from PDSYTRD;
*          = 'L': Lower triangle of A(IA:*,JA:*) contains elementary
*                 reflectors from PDSYTRD.
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
*  A       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_A,LOCc(JA+M-1)) if SIDE='L',
*          or (LLD_A,LOCc(JA+N-1)) if SIDE = 'R'. The vectors which
*          define the elementary reflectors, as returned by PDSYTRD.
*          If SIDE = 'L', LLD_A >= max(1,LOCr(IA+M-1));
*          if SIDE = 'R', LLD_A >= max(1,LOCr(IA+N-1)).
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
*  TAU     (local input) DOUBLE PRECISION array, dimension LTAU, where
*          if SIDE = 'L' and UPLO = 'U', LTAU = LOCc(M_A),
*          if SIDE = 'L' and UPLO = 'L', LTAU = LOCc(JA+M-2),
*          if SIDE = 'R' and UPLO = 'U', LTAU = LOCc(N_A),
*          if SIDE = 'R' and UPLO = 'L', LTAU = LOCc(JA+N-2).
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by PDSYTRD. TAU is tied to the
*          distributed matrix A.
*
*  C       (local input/local output) DOUBLE PRECISION pointer into the
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
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                     dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input)  INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*
*          If UPLO = 'U',
*            IAA = IA, JAA = JA+1, ICC = IC, JCC = JC;
*          else UPLO = 'L',
*            IAA = IA+1, JAA = JA;
*            if SIDE = 'L',
*              ICC = IC+1; JCC = JC;
*            else
*              ICC = IC; JCC = JC+1;
*            end if
*          end if
*
*          If SIDE = 'L',
*            MI = M-1; NI = N;
*            LWORK >= MAX( (NB_A*(NB_A-1))/2, (NqC0 + MpC0)*NB_A ) +
*                     NB_A * NB_A
*          else if SIDE = 'R',
*            MI = M; MI = N-1;
*            LWORK >= MAX( (NB_A*(NB_A-1))/2, ( NqC0 + MAX( NpA0 +
*                     NUMROC( NUMROC( NI+ICOFFC, NB_A, 0, 0, NPCOL ),
*                             NB_A, 0, 0, LCMQ ), MpC0 ) )*NB_A ) +
*                     NB_A * NB_A
*          end if
*
*          where LCMQ = LCM / NPCOL with LCM = ICLM( NPROW, NPCOL ),
*
*          IROFFA = MOD( IAA-1, MB_A ), ICOFFA = MOD( JAA-1, NB_A ),
*          IAROW = INDXG2P( IAA, MB_A, MYROW, RSRC_A, NPROW ),
*          NpA0 = NUMROC( NI+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*
*          IROFFC = MOD( ICC-1, MB_C ), ICOFFC = MOD( JCC-1, NB_C ),
*          ICROW = INDXG2P( ICC, MB_C, MYROW, RSRC_C, NPROW ),
*          ICCOL = INDXG2P( JCC, NB_C, MYCOL, CSRC_C, NPCOL ),
*          MpC0 = NUMROC( MI+IROFFC, MB_C, MYROW, ICROW, NPROW ),
*          NqC0 = NUMROC( NI+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL ),
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
*  INFO    (global output) INTEGER
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
*    ( MB_A.EQ.MB_C .AND. IROFFA.EQ.IROFFC .AND. IAROW.EQ.ICROW )
*  If SIDE = 'R',
*    ( MB_A.EQ.NB_C .AND. IROFFA.EQ.ICOFFC )
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
      LOGICAL            LEFT, LQUERY, NOTRAN, UPPER
      INTEGER            IAA, IAROW, ICC, ICCOL, ICOFFC, ICROW, ICTXT,
     $                   IINFO, IROFFA, IROFFC, JAA, JCC, LCM, LCMQ,
     $                   LWMIN, MI, MPC0, MYCOL, MYROW, NI, NPA0, NPCOL,
     $                   NPROW, NQ, NQC0
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 4 ), IDUM2( 4 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK2MAT, PDORMQL,
     $                   PDORMQR, PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILCM, INDXG2P, NUMROC
      EXTERNAL           ILCM, INDXG2P, LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, ICHAR, MAX, MOD
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
         UPPER = LSAME( UPLO, 'U' )
*
         IF( UPPER ) THEN
            IAA = IA
            JAA = JA+1
            ICC = IC
            JCC = JC
         ELSE
            IAA = IA+1
            JAA = JA
            IF( LEFT ) THEN
               ICC = IC + 1
               JCC = JC
            ELSE
               ICC = IC
               JCC = JC + 1
            END IF
         END IF
*
*        NQ is the order of Q
*
         IF( LEFT ) THEN
            NQ = M
            MI = M - 1
            NI = N
            CALL CHK1MAT( MI, 4, NQ-1, 4, IAA, JAA, DESCA, 9, INFO )
         ELSE
            NQ = N
            MI = M
            NI = N - 1
            CALL CHK1MAT( NI, 5, NQ-1, 5, IAA, JAA, DESCA, 9, INFO )
         END IF
         CALL CHK1MAT( MI, 4, NI, 5, ICC, JCC, DESCC, 14, INFO )
         IF( INFO.EQ.0 ) THEN
            IROFFA = MOD( IAA-1, DESCA( MB_ ) )
            IROFFC = MOD( ICC-1, DESCC( MB_ ) )
            ICOFFC = MOD( JCC-1, DESCC( NB_ ) )
            IAROW = INDXG2P( IAA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            ICROW = INDXG2P( ICC, DESCC( MB_ ), MYROW, DESCC( RSRC_ ),
     $                       NPROW )
            ICCOL = INDXG2P( JCC, DESCC( NB_ ), MYCOL, DESCC( CSRC_ ),
     $                       NPCOL )
            MPC0 = NUMROC( MI+IROFFC, DESCC( MB_ ), MYROW, ICROW,
     $                     NPROW )
            NQC0 = NUMROC( NI+ICOFFC, DESCC( NB_ ), MYCOL, ICCOL,
     $                     NPCOL )
*
            IF( LEFT ) THEN
               LWMIN = MAX( ( DESCA( NB_ ) * ( DESCA( NB_ ) - 1 ) ) / 2,
     $                      ( MPC0 + NQC0 ) * DESCA( NB_ ) ) +
     $                 DESCA( NB_ ) * DESCA( NB_ )
            ELSE
               NPA0 = NUMROC( NI+IROFFA, DESCA( MB_ ), MYROW, IAROW,
     $                        NPROW )
               LCM = ILCM( NPROW, NPCOL )
               LCMQ = LCM / NPCOL
               LWMIN =  MAX( ( DESCA( NB_ ) * ( DESCA( NB_ ) - 1 ) )
     $                  / 2, ( NQC0 + MAX( NPA0 + NUMROC( NUMROC(
     $                  NI+ICOFFC, DESCA( NB_ ), 0, 0, NPCOL ),
     $                  DESCA( NB_ ), 0, 0, LCMQ ), MPC0 ) ) *
     $                  DESCA( NB_ ) ) + DESCA( NB_ ) * DESCA( NB_ )
            END IF
*
            WORK( 1 ) = DBLE( LWMIN )
            LQUERY = ( LWORK.EQ.-1 )
            IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -2
            ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND.
     $               .NOT.LSAME( TRANS, 'T' ) ) THEN
               INFO = -3
            ELSE IF( .NOT.LEFT .AND. DESCA( MB_ ).NE.DESCC( NB_ ) ) THEN
               INFO = -(900+NB_)
            ELSE IF( LEFT .AND. IROFFA.NE.IROFFC ) THEN
               INFO = -12
            ELSE IF( LEFT .AND. IAROW.NE.ICROW ) THEN
               INFO = -12
            ELSE IF( .NOT.LEFT .AND. IROFFA.NE.ICOFFC ) THEN
               INFO = -13
            ELSE IF( LEFT .AND. DESCA( MB_ ).NE.DESCC( MB_ ) ) THEN
               INFO = -(1400+MB_)
            ELSE IF( ICTXT.NE.DESCC( CTXT_ ) ) THEN
               INFO = -(1400+CTXT_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -16
            END IF
         END IF
*
         IF( LEFT ) THEN
            IDUM1( 1 ) = ICHAR( 'L' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'R' )
         END IF
         IDUM2( 1 ) = 1
         IF( UPPER ) THEN
            IDUM1( 2 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'L' )
         END IF
         IDUM2( 2 ) = 2
         IF( NOTRAN ) THEN
            IDUM1( 3 ) = ICHAR( 'N' )
         ELSE
            IDUM1( 3 ) = ICHAR( 'T' )
         END IF
         IDUM2( 3 ) = 3
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 4 ) = -1
         ELSE
            IDUM1( 4 ) = 1
         END IF
         IDUM2( 4 ) = 16
         IF( LEFT ) THEN
            CALL PCHK2MAT( MI, 4, NQ-1, 4, IAA, JAA, DESCA, 9, MI, 4,
     $                     NI, 5, ICC, JCC, DESCC, 14, 4, IDUM1, IDUM2,
     $                     INFO )
         ELSE
            CALL PCHK2MAT( NI, 5, NQ-1, 5, IAA, JAA, DESCA, 9, MI, 4,
     $                     NI, 5, ICC, JCC, DESCC, 14, 4, IDUM1, IDUM2,
     $                     INFO )
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDORMTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. NQ.EQ.1 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Q was determined by a call to PDSYTRD with UPLO = 'U'
*
         CALL PDORMQL( SIDE, TRANS, MI, NI, NQ-1, A, IAA, JAA, DESCA,
     $                 TAU, C, ICC, JCC, DESCC, WORK, LWORK, IINFO )
*
      ELSE
*
*        Q was determined by a call to PDSYTRD with UPLO = 'L'
*
         CALL PDORMQR( SIDE, TRANS, MI, NI, NQ-1, A, IAA, JAA, DESCA,
     $                 TAU, C, ICC, JCC, DESCC, WORK, LWORK, IINFO )
*
      END IF
*
      WORK( 1 ) = DBLE( LWMIN )
*
      RETURN
*
*     End of PDORMTR
*
      END
