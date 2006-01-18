      SUBROUTINE PZUNMBR( VECT, SIDE, TRANS, M, N, K, A, IA, JA, DESCA,
     $                    TAU, C, IC, JC, DESCC, WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, VECT
      INTEGER            IA, IC, INFO, JA, JC, K, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCC( * )
      COMPLEX*16         A( * ), C( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  If VECT = 'Q', PZUNMBR overwrites the general complex distributed
*  M-by-N matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1) with
*
*                       SIDE = 'L'           SIDE = 'R'
*  TRANS = 'N':      Q * sub( C )          sub( C ) * Q
*  TRANS = 'C':      Q**H * sub( C )       sub( C ) * Q**H
*
*  If VECT = 'P', PZUNMBR overwrites sub( C ) with
*
*                       SIDE = 'L'           SIDE = 'R'
*  TRANS = 'N':      P * sub( C )          sub( C ) * P
*  TRANS = 'C':      P**H * sub( C )       sub( C ) * P**H
*
*  Here Q and P**H are the unitary distributed matrices determined by
*  PZGEBRD when reducing a complex distributed matrix A(IA:*,JA:*) to
*  bidiagonal form: A(IA:*,JA:*) = Q * B * P**H. Q and P**H are defined
*  as products of elementary reflectors H(i) and G(i) respectively.
*
*  Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
*  order of the unitary matrix Q or P**H that is applied.
*
*  If VECT = 'Q', A(IA:*,JA:*) is assumed to have been an NQ-by-K
*  matrix:
*  if nq >= k, Q = H(1) H(2) . . . H(k);
*  if nq < k, Q = H(1) H(2) . . . H(nq-1).
*
*  If VECT = 'P', A(IA:*,JA:*) is assumed to have been a K-by-NQ
*  matrix:
*  if k < nq, P = G(1) G(2) . . . G(k);
*  if k >= nq, P = G(1) G(2) . . . G(nq-1).
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
*  VECT    (global input) CHARACTER
*          = 'Q': apply Q or Q**H;
*          = 'P': apply P or P**H.
*
*  SIDE    (global input) CHARACTER
*          = 'L': apply Q, Q**H, P or P**H from the Left;
*          = 'R': apply Q, Q**H, P or P**H from the Right.
*
*  TRANS   (global input) CHARACTER
*          = 'N':  No transpose, apply Q or P;
*          = 'C':  Conjugate transpose, apply Q**H or P**H.
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
*          If VECT = 'Q', the number of columns in the original
*          distributed matrix reduced by PZGEBRD.
*          If VECT = 'P', the number of rows in the original
*          distributed matrix reduced by PZGEBRD.
*          K >= 0.
*
*  A       (local input) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_A,LOCc(JA+MIN(NQ,K)-1)) if
*          VECT='Q', and (LLD_A,LOCc(JA+NQ-1)) if VECT = 'P'. NQ = M
*          if SIDE = 'L', and NQ = N otherwise. The vectors which
*          define the elementary reflectors H(i) and G(i), whose
*          products determine the matrices Q and P, as returned by
*          PZGEBRD.
*          If VECT = 'Q', LLD_A >= max(1,LOCr(IA+NQ-1));
*          if VECT = 'P', LLD_A >= max(1,LOCr(IA+MIN(NQ,K)-1)).
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
*  TAU     (local input) COMPLEX*16 array, dimension
*          LOCc(JA+MIN(NQ,K)-1) if VECT = 'Q', LOCr(IA+MIN(NQ,K)-1) if
*          VECT = 'P', TAU(i) must contain the scalar factor of the
*          elementary  reflector H(i) or G(i), which determines Q or P,
*          as returned by PDGEBRD in its array argument TAUQ or TAUP.
*          TAU is tied to the distributed matrix A.
*
*  C       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_C,LOCc(JC+N-1)).
*          On entry, the local pieces of the distributed matrix sub(C).
*          On exit, if VECT='Q', sub( C ) is overwritten by Q*sub( C )
*          or Q'*sub( C ) or sub( C )*Q' or sub( C )*Q; if VECT='P,
*          sub( C ) is overwritten by P*sub( C ) or P'*sub( C ) or
*          sub( C )*P or sub( C )*P'.
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
*  WORK    (local workspace/local output) COMPLEX*16 array,
*                                                     dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          If SIDE = 'L',
*            NQ = M;
*            if( (VECT = 'Q' and NQ >= K) or (VECT <> 'Q' and NQ > K) ),
*               IAA=IA; JAA=JA; MI=M; NI=N; ICC=IC; JCC=JC;
*            else
*               IAA=IA+1; JAA=JA; MI=M-1; NI=N; ICC=IC+1; JCC=JC;
*            end if
*          else if SIDE = 'R',
*            NQ = N;
*            if( (VECT = 'Q' and NQ >= K) or (VECT <> 'Q' and NQ > K) ),
*               IAA=IA; JAA=JA; MI=M; NI=N; ICC=IC; JCC=JC;
*            else
*               IAA=IA; JAA=JA+1; MI=M; NI=N-1; ICC=IC; JCC=JC+1;
*            end if
*          end if
*
*          If VECT = 'Q',
*            If SIDE = 'L',
*              LWORK >= MAX( (NB_A*(NB_A-1))/2, (NqC0 + MpC0)*NB_A ) +
*                       NB_A * NB_A
*            else if SIDE = 'R',
*              LWORK >= MAX( (NB_A*(NB_A-1))/2, ( NqC0 + MAX( NpA0 +
*                       NUMROC( NUMROC( NI+ICOFFC, NB_A, 0, 0, NPCOL ),
*                               NB_A, 0, 0, LCMQ ), MpC0 ) )*NB_A ) +
*                       NB_A * NB_A
*            end if
*          else if VECT <> 'Q',
*            if SIDE = 'L',
*              LWORK >= MAX( (MB_A*(MB_A-1))/2, ( MpC0 + MAX( MqA0 +
*                       NUMROC( NUMROC( MI+IROFFC, MB_A, 0, 0, NPROW ),
*                               MB_A, 0, 0, LCMP ), NqC0 ) )*MB_A ) +
*                       MB_A * MB_A
*            else if SIDE = 'R',
*              LWORK >= MAX( (MB_A*(MB_A-1))/2, (MpC0 + NqC0)*MB_A ) +
*                       MB_A * MB_A
*            end if
*          end if
*
*          where LCMP = LCM / NPROW, LCMQ = LCM / NPCOL, with
*          LCM = ICLM( NPROW, NPCOL ),
*
*          IROFFA = MOD( IAA-1, MB_A ), ICOFFA = MOD( JAA-1, NB_A ),
*          IAROW = INDXG2P( IAA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JAA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          MqA0 = NUMROC( MI+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),
*          NpA0 = NUMROC( NI+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*
*          IROFFC = MOD( ICC-1, MB_C ), ICOFFC = MOD( JCC-1, NB_C ),
*          ICROW = INDXG2P( ICC, MB_C, MYROW, RSRC_C, NPROW ),
*          ICCOL = INDXG2P( JCC, NB_C, MYCOL, CSRC_C, NPCOL ),
*          MpC0 = NUMROC( MI+IROFFC, MB_C, MYROW, ICROW, NPROW ),
*          NqC0 = NUMROC( NI+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL ),
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
*  Alignment requirements
*  ======================
*
*  The distributed submatrices A(IA:*, JA:*) and C(IC:IC+M-1,JC:JC+N-1)
*  must verify some alignment properties, namely the following
*  expressions should be true:
*
*  If VECT = 'Q',
*    If SIDE = 'L',
*      ( MB_A.EQ.MB_C .AND. IROFFA.EQ.IROFFC .AND. IAROW.EQ.ICROW )
*     If SIDE = 'R',
*      ( MB_A.EQ.NB_C .AND. IROFFA.EQ.ICOFFC )
*  else
*     If SIDE = 'L',
*       ( MB_A.EQ.MB_C .AND. ICOFFA.EQ.IROFFC )
*     If SIDE = 'R',
*       ( NB_A.EQ.NB_C .AND. ICOFFA.EQ.ICOFFC .AND. IACOL.EQ.ICCOL )
*  end if
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
      LOGICAL            APPLYQ, LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            IAA, IACOL, IAROW, ICC, ICCOL, ICOFFA, ICOFFC,
     $                   ICROW, ICTXT, IINFO, IROFFA, IROFFC, JAA, JCC,
     $                   LCM, LCMP, LCMQ, LWMIN, MI, MPC0, MQA0, MYCOL,
     $                   MYROW, NI, NPA0, NPCOL, NPROW, NQ, NQC0
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 5 ), IDUM2( 5 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK1MAT, PXERBLA,
     $                   PZUNMLQ, PZUNMQR
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILCM, INDXG2P, NUMROC
      EXTERNAL           ILCM, INDXG2P, LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, ICHAR, MAX, MOD
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
         INFO = -(1000+CTXT_)
      ELSE
         APPLYQ = LSAME( VECT, 'Q' )
         LEFT = LSAME( SIDE, 'L' )
         NOTRAN = LSAME( TRANS, 'N' )
*
*        NQ is the order of Q or P
*
         IF( LEFT ) THEN
            NQ = M
            IF( ( APPLYQ .AND. NQ.GE.K ) .OR.
     $          ( .NOT.APPLYQ .AND. NQ.GT.K ) ) THEN
               IAA = IA
               JAA = JA
               MI = M
               NI = N
               ICC = IC
               JCC = JC
            ELSE
               IAA = IA + 1
               JAA = JA
               MI = M - 1
               NI = N
               ICC = IC + 1
               JCC = JC
            END IF
*
            IF( APPLYQ ) THEN
               CALL CHK1MAT( M, 4, K, 6, IA, JA, DESCA, 10, INFO )
            ELSE
               CALL CHK1MAT( K, 6, M, 4, IA, JA, DESCA, 10, INFO )
            END IF
         ELSE
            NQ = N
            IF( ( APPLYQ .AND. NQ.GE.K ) .OR.
     $          ( .NOT.APPLYQ .AND. NQ.GT.K ) ) THEN
               IAA = IA
               JAA = JA
               MI = M
               NI = N
               ICC = IC
               JCC = JC
            ELSE
               IAA = IA
               JAA = JA + 1
               MI = M
               NI = N - 1
               ICC = IC
               JCC = JC + 1
            END IF
*
            IF( APPLYQ ) THEN
               CALL CHK1MAT( N, 5, K, 6, IA, JA, DESCA, 10, INFO )
            ELSE
               CALL CHK1MAT( K, 6, N, 5, IA, JA, DESCA, 10, INFO )
            END IF
         END IF
         CALL CHK1MAT( M, 4, N, 5, IC, JC, DESCC, 15, INFO )
*
         IF( INFO.EQ.0 ) THEN
            IROFFA = MOD( IAA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JAA-1, DESCA( NB_ ) )
            IROFFC = MOD( ICC-1, DESCC( MB_ ) )
            ICOFFC = MOD( JCC-1, DESCC( NB_ ) )
            IACOL = INDXG2P( JAA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
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
            IF( APPLYQ ) THEN
               IF( LEFT ) THEN
                  LWMIN = MAX( ( DESCA( NB_ ) * ( DESCA( NB_ ) - 1 ) )
     $                    / 2, ( MPC0 + NQC0 ) * DESCA( NB_ ) ) +
     $                    DESCA( NB_ ) * DESCA( NB_ )
               ELSE
                  NPA0 = NUMROC( NI+IROFFA, DESCA( MB_ ), MYROW, IAROW,
     $                           NPROW )
                  LCM = ILCM( NPROW, NPCOL )
                  LCMQ = LCM / NPCOL
                  LWMIN =  MAX( ( DESCA( NB_ ) * ( DESCA( NB_ ) - 1 ) )
     $                     / 2, ( NQC0 + MAX( NPA0 + NUMROC( NUMROC(
     $                     NI+ICOFFC, DESCA( NB_ ), 0, 0, NPCOL ),
     $                     DESCA( NB_ ), 0, 0, LCMQ ), MPC0 ) ) *
     $                     DESCA( NB_ ) ) + DESCA( NB_ ) * DESCA( NB_ )
               END IF
            ELSE
*
               IF( LEFT ) THEN
                  MQA0 = NUMROC( MI+ICOFFA, DESCA( NB_ ), MYCOL, IACOL,
     $                           NPCOL )
                  LCM = ILCM( NPROW, NPCOL )
                  LCMP = LCM / NPROW
                  LWMIN =  MAX( ( DESCA( MB_ ) * ( DESCA( MB_ ) - 1 ) )
     $                     / 2, ( MPC0 + MAX( MQA0 + NUMROC( NUMROC(
     $                     MI+IROFFC, DESCA( MB_ ), 0, 0, NPROW ),
     $                     DESCA( MB_ ), 0, 0, LCMP ), NQC0 ) ) *
     $                     DESCA( MB_ ) ) + DESCA( MB_ ) * DESCA( MB_ )
               ELSE
                  LWMIN = MAX( ( DESCA( MB_ ) * ( DESCA( MB_ ) - 1 ) )
     $                    / 2, ( MPC0 + NQC0 ) * DESCA( MB_ ) ) +
     $                    DESCA( MB_ ) * DESCA( MB_ )
               END IF
*
            END IF
*
            WORK( 1 ) = DCMPLX( DBLE( LWMIN ) )
            LQUERY = ( LWORK.EQ.-1 )
            IF( .NOT.APPLYQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
               INFO = -2
            ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
               INFO = -3
            ELSE IF( K.LT.0 ) THEN
               INFO = -6
            ELSE IF( APPLYQ .AND. .NOT.LEFT .AND.
     $               DESCA( MB_ ).NE.DESCC( NB_ ) ) THEN
               INFO = -(1000+NB_)
            ELSE IF( APPLYQ .AND. LEFT .AND. IROFFA.NE.IROFFC ) THEN
               INFO = -13
            ELSE IF( APPLYQ .AND. LEFT .AND. IAROW.NE.ICROW ) THEN
               INFO = -13
            ELSE IF( .NOT.APPLYQ .AND. LEFT .AND.
     $               ICOFFA.NE.IROFFC ) THEN
               INFO = -13
            ELSE IF( .NOT.APPLYQ .AND. .NOT.LEFT .AND.
     $               IACOL.NE.ICCOL ) THEN
               INFO = -14
            ELSE IF( APPLYQ .AND. .NOT.LEFT .AND.
     $               IROFFA.NE.ICOFFC ) THEN
               INFO = -14
            ELSE IF( .NOT.APPLYQ .AND. .NOT.LEFT .AND.
     $               ICOFFA.NE.ICOFFC ) THEN
               INFO = -14
            ELSE IF( APPLYQ .AND. LEFT .AND.
     $               DESCA( MB_ ).NE.DESCC( MB_ ) ) THEN
               INFO = -(1500+MB_)
            ELSE IF( .NOT.APPLYQ .AND. LEFT .AND.
     $               DESCA( MB_ ).NE.DESCC( MB_ ) ) THEN
               INFO = -(1500+MB_)
            ELSE IF( APPLYQ .AND. .NOT.LEFT .AND.
     $               DESCA( MB_ ).NE.DESCC( NB_ ) ) THEN
               INFO = -(1500+NB_)
            ELSE IF( .NOT.APPLYQ .AND. .NOT.LEFT .AND.
     $               DESCA( NB_ ).NE.DESCC( NB_ ) ) THEN
               INFO = -(1500+NB_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -17
            END IF
         END IF
*
         IF( APPLYQ ) THEN
            IDUM1( 1 ) = ICHAR( 'Q' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'P' )
         END IF
         IDUM2( 1 ) = 1
         IF( LEFT ) THEN
            IDUM1( 2 ) = ICHAR( 'L' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'R' )
         END IF
         IDUM2( 2 ) = 2
         IF( NOTRAN ) THEN
            IDUM1( 3 ) = ICHAR( 'N' )
         ELSE
            IDUM1( 3 ) = ICHAR( 'C' )
         END IF
         IDUM2( 3 ) = 3
         IDUM1( 4 ) = K
         IDUM2( 4 ) = 6
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 5 ) = -1
         ELSE
            IDUM1( 5 ) = 1
         END IF
         IDUM2( 5 ) = 17
         IF( APPLYQ ) THEN
            IF( LEFT ) THEN
               CALL PCHK2MAT( M, 4, K, 6, IA, JA, DESCA, 10, M, 4, N,
     $                        5, IC, JC, DESCC, 15, 5, IDUM1, IDUM2,
     $                        INFO )
            ELSE
               CALL PCHK2MAT( N, 5, K, 6, IA, JA, DESCA, 10, M, 4, N,
     $                        5, IC, JC, DESCC, 15, 5, IDUM1, IDUM2,
     $                        INFO )
            END IF
         ELSE
            IF( LEFT ) THEN
               CALL PCHK2MAT( K, 6, M, 4, IA, JA, DESCA, 10, M, 4, N,
     $                        5, IC, JC, DESCC, 15, 5, IDUM1, IDUM2,
     $                        INFO )
            ELSE
               CALL PCHK2MAT( K, 6, N, 5, IA, JA, DESCA, 10, M, 4, N,
     $                        5, IC, JC, DESCC, 15, 5, IDUM1, IDUM2,
     $                        INFO )
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZUNMBR', -INFO )
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
      IF( APPLYQ ) THEN
*
*        Apply Q
*
         IF( NQ.GE.K ) THEN
*
*           Q was determined by a call to PZGEBRD with nq >= k
*
            CALL PZUNMQR( SIDE, TRANS, M, N, K, A, IA, JA, DESCA, TAU,
     $                    C, IC, JC, DESCC, WORK, LWORK, IINFO )
         ELSE IF( NQ.GT.1 ) THEN
*
*           Q was determined by a call to PZGEBRD with nq < k
*
            CALL PZUNMQR( SIDE, TRANS, MI, NI, NQ-1, A, IA+1, JA, DESCA,
     $                    TAU, C, ICC, JCC, DESCC, WORK, LWORK, IINFO )
         END IF
      ELSE
*
*        Apply P
*
         IF( NOTRAN ) THEN
            TRANST = 'C'
         ELSE
            TRANST = 'N'
         END IF
         IF( NQ.GT.K ) THEN
*
*           P was determined by a call to PZGEBRD with nq > k
*
            CALL PZUNMLQ( SIDE, TRANST, M, N, K, A, IA, JA, DESCA, TAU,
     $                    C, IC, JC, DESCC, WORK, LWORK, IINFO )
         ELSE IF( NQ.GT.1 ) THEN
*
*           P was determined by a call to PZGEBRD with nq <= k
*
            CALL PZUNMLQ( SIDE, TRANST, MI, NI, NQ-1, A, IA, JA+1,
     $                    DESCA, TAU, C, ICC, JCC, DESCC, WORK, LWORK,
     $                    IINFO )
         END IF
      END IF
*
      WORK( 1 ) = DCMPLX( DBLE( LWMIN ) )
*
      RETURN
*
*     End of PZUNMBR
*
      END
