      SUBROUTINE PCUNMHR( SIDE, TRANS, M, N, ILO, IHI, A, IA, JA, DESCA,
     $                    TAU, C, IC, JC, DESCC, WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            IA, IC, IHI, ILO, INFO, JA, JC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCC( * )
      COMPLEX            A( * ), C( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PCUNMHR overwrites the general complex M-by-N distributed matrix
*  sub( C ) = C(IC:IC+M-1,JC:JC+N-1) with
*
*                       SIDE = 'L'           SIDE = 'R'
*  TRANS = 'N':      Q * sub( C )          sub( C ) * Q
*  TRANS = 'C':      Q**H * sub( C )       sub( C ) * Q**H
*
*  where Q is a complex unitary distributed matrix of order nq, with
*  nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the
*  product of IHI-ILO elementary reflectors, as returned by PCGEHRD:
*
*  Q = H(ilo) H(ilo+1) . . . H(ihi-1).
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
*          = 'L': apply Q or Q**H from the Left;
*          = 'R': apply Q or Q**H from the Right.
*
*  TRANS   (global input) CHARACTER
*          = 'N':  No transpose, apply Q;
*          = 'C':  Conjugate transpose, apply Q**H.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( C ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( C ). N >= 0.
*
*  ILO     (global input) INTEGER
*  IHI     (global input) INTEGER
*          ILO and IHI must have the same values as in the previous call
*          of PCGEHRD. Q is equal to the unit matrix except in the
*          distributed submatrix Q(ia+ilo:ia+ihi-1,ia+ilo:ja+ihi-1).
*          If SIDE = 'L', 1 <= ILO <= IHI <= max(1,M);
*          if SIDE = 'R', 1 <= ILO <= IHI <= max(1,N);
*          ILO and IHI are relative indexes.
*
*  A       (local input) COMPLEX pointer into the local memory
*          to an array of dimension (LLD_A,LOCc(JA+M-1)) if SIDE='L',
*          and (LLD_A,LOCc(JA+N-1)) if SIDE = 'R'. The vectors which
*          define the elementary reflectors, as returned by PCGEHRD.
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
*  TAU     (local input) COMPLEX, array, dimension LOCc(JA+M-2)
*          if SIDE = 'L', and LOCc(JA+N-2) if SIDE = 'R'. This array
*          contains the scalar factors TAU(j) of the elementary
*          reflectors H(j) as returned by PCGEHRD. TAU is tied to
*          the distributed matrix A.
*
*  C       (local input/local output) COMPLEX pointer into the
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
*  WORK    (local workspace/local output) COMPLEX array,
*                                                     dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*
*          IAA = IA + ILO; JAA = JA+ILO-1;
*          If SIDE = 'L',
*            MI = IHI-ILO; NI = N; ICC = IC + ILO; JCC = JC;
*            LWORK >= MAX( (NB_A*(NB_A-1))/2, (NqC0 + MpC0)*NB_A ) +
*                     NB_A * NB_A
*          else if SIDE = 'R',
*            MI = M; NI = IHI-ILO; ICC = IC; JCC = JC + ILO;
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
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            IAA, IAROW, ICC, ICCOL, ICOFFC, ICROW, ICTXT,
     $                   IINFO, IROFFA, IROFFC, JAA, JCC, LCM, LCMQ,
     $                   LWMIN, MI, MPC0, MYCOL, MYROW, NH, NI, NPA0,
     $                   NPCOL, NPROW, NQ, NQC0
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 5 ), IDUM2( 5 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK2MAT, PCUNMQR,
     $                   PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILCM, INDXG2P, NUMROC
      EXTERNAL           ILCM, INDXG2P, LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, ICHAR, MAX, MIN, MOD, REAL
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
      NH = IHI - ILO
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(1000+CTXT_)
      ELSE
         LEFT = LSAME( SIDE, 'L' )
         NOTRAN = LSAME( TRANS, 'N' )
         IAA = IA + ILO
         JAA = JA + ILO - 1
*
*        NQ is the order of Q
*
         IF( LEFT ) THEN
            NQ = M
            MI = NH
            NI = N
            ICC = IC + ILO
            JCC = JC
            CALL CHK1MAT( M, 3, M, 3, IA, JA, DESCA, 10, INFO )
         ELSE
            NQ = N
            MI = M
            NI = NH
            ICC = IC
            JCC = JC + ILO
            CALL CHK1MAT( N, 4, N, 4, IA, JA, DESCA, 10, INFO )
         END IF
         CALL CHK1MAT( M, 3, N, 4, IC, JC, DESCC, 15, INFO )
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
            WORK( 1 ) = CMPLX( REAL( LWMIN ) )
            LQUERY = ( LWORK.EQ.-1 )
            IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
               INFO = -2
            ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, NQ ) ) THEN
               INFO = -5
            ELSE IF( IHI.LT.MIN( ILO, NQ ) .OR. IHI.GT.NQ ) THEN
               INFO = -6
            ELSE IF( .NOT.LEFT .AND. DESCA( MB_ ).NE.DESCC( NB_ ) ) THEN
               INFO = -(1000+NB_)
            ELSE IF( LEFT .AND. IROFFA.NE.IROFFC ) THEN
               INFO = -13
            ELSE IF( LEFT .AND. IAROW.NE.ICROW ) THEN
               INFO = -13
            ELSE IF( .NOT.LEFT .AND. IROFFA.NE.ICOFFC ) THEN
               INFO = -14
            ELSE IF( LEFT .AND. DESCA( MB_ ).NE.DESCC( MB_ ) ) THEN
               INFO = -(1500+MB_)
            ELSE IF( ICTXT.NE.DESCC( CTXT_ ) ) THEN
               INFO = -(1500+CTXT_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -17
            END IF
         END IF
*
         IF( LEFT ) THEN
            IDUM1( 1 ) = ICHAR( 'L' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'R' )
         END IF
         IDUM2( 1 ) = 1
         IF( NOTRAN ) THEN
            IDUM1( 2 ) = ICHAR( 'N' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'C' )
         END IF
         IDUM2( 2 ) = 2
         IDUM1( 3 ) = ILO
         IDUM2( 3 ) = 5
         IDUM1( 4 ) = IHI
         IDUM2( 4 ) = 6
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 5 ) = -1
         ELSE
            IDUM1( 5 ) = 1
         END IF
         IDUM2( 5 ) = 17
         IF( LEFT ) THEN
            CALL PCHK2MAT( M, 3, M, 3, IA, JA, DESCA, 10, M, 3, N, 4,
     $                     IC, JC, DESCC, 15, 5, IDUM1, IDUM2, INFO )
         ELSE
            CALL PCHK2MAT( N, 4, N, 4, IA, JA, DESCA, 10, M, 3, N, 4,
     $                     IC, JC, DESCC, 15, 5, IDUM1, IDUM2, INFO )
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCUNMHR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. NH.EQ.0 )
     $   RETURN
*
      CALL PCUNMQR( SIDE, TRANS, MI, NI, NH, A, IAA, JAA, DESCA, TAU,
     $              C, ICC, JCC, DESCC, WORK, LWORK, IINFO )
*
      WORK( 1 ) = CMPLX( REAL( LWMIN ) )
*
      RETURN
*
*     End of PCUNMHR
*
      END
