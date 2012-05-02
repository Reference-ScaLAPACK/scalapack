      SUBROUTINE PCLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V,
     $                    IV, JV, DESCV, T, C, IC, JC, DESCC, WORK )
*
*  -- ScaLAPACK auxiliary routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            IC, IV, JC, JV, K, L, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCC( * ), DESCV( * )
      COMPLEX            C( * ), T( * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLARZB applies a complex block reflector Q or its conjugate
*  transpose Q**H to a complex M-by-N distributed matrix sub( C )
*  denoting C(IC:IC+M-1,JC:JC+N-1), from the left or the right.
*
*  Q is a product of k elementary reflectors as returned by PCTZRZF.
*
*  Currently, only STOREV = 'R' and DIRECT = 'B' are supported.
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
*  DIRECT  (global input) CHARACTER
*          Indicates how H is formed from a product of elementary
*          reflectors
*          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (global input) CHARACTER
*          Indicates how the vectors which define the elementary
*          reflectors are stored:
*          = 'C': Columnwise                        (not supported yet)
*          = 'R': Rowwise
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
*          The order of the matrix T (= the number of elementary
*          reflectors whose product defines the block reflector).
*
*  L       (global input) INTEGER
*          The columns of the distributed submatrix sub( A ) containing
*          the meaningful part of the Householder reflectors.
*          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
*
*  V       (local input) COMPLEX pointer into the local memory
*          to an array of dimension (LLD_V, LOCc(JV+M-1)) if SIDE = 'L',
*          (LLD_V, LOCc(JV+N-1)) if SIDE = 'R'. It contains the local
*          pieces of the distributed vectors V representing the
*          Householder transformation as returned by PCTZRZF.
*          LLD_V >= LOCr(IV+K-1).
*
*  IV      (global input) INTEGER
*          The row index in the global array V indicating the first
*          row of sub( V ).
*
*  JV      (global input) INTEGER
*          The column index in the global array V indicating the
*          first column of sub( V ).
*
*  DESCV   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix V.
*
*  T       (local input) COMPLEX array, dimension MB_V by MB_V
*          The lower triangular matrix T in the representation of the
*          block reflector.
*
*  C       (local input/local output) COMPLEX pointer into the
*          local memory to an array of dimension (LLD_C,LOCc(JC+N-1)).
*          On entry, the M-by-N distributed matrix sub( C ). On exit,
*          sub( C ) is overwritten by Q*sub( C ) or Q'*sub( C ) or
*          sub( C )*Q or sub( C )*Q'.
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
*  WORK    (local workspace) COMPLEX array, dimension (LWORK)
*          If STOREV = 'C',
*            if SIDE = 'L',
*              LWORK >= ( NqC0 + MpC0 ) * K
*            else if SIDE = 'R',
*              LWORK >= ( NqC0 + MAX( NpV0 + NUMROC( NUMROC( N+ICOFFC,
*                         NB_V, 0, 0, NPCOL ), NB_V, 0, 0, LCMQ ),
*                         MpC0 ) ) * K
*            end if
*          else if STOREV = 'R',
*            if SIDE = 'L',
*              LWORK >= ( MpC0 + MAX( MqV0 + NUMROC( NUMROC( M+IROFFC,
*                         MB_V, 0, 0, NPROW ), MB_V, 0, 0, LCMP ),
*                         NqC0 ) ) * K
*            else if SIDE = 'R',
*              LWORK >= ( MpC0 + NqC0 ) * K
*            end if
*          end if
*
*          where LCMQ = LCM / NPCOL with LCM = ICLM( NPROW, NPCOL ),
*
*          IROFFV = MOD( IV-1, MB_V ), ICOFFV = MOD( JV-1, NB_V ),
*          IVROW = INDXG2P( IV, MB_V, MYROW, RSRC_V, NPROW ),
*          IVCOL = INDXG2P( JV, NB_V, MYCOL, CSRC_V, NPCOL ),
*          MqV0 = NUMROC( M+ICOFFV, NB_V, MYCOL, IVCOL, NPCOL ),
*          NpV0 = NUMROC( N+IROFFV, MB_V, MYROW, IVROW, NPROW ),
*
*          IROFFC = MOD( IC-1, MB_C ), ICOFFC = MOD( JC-1, NB_C ),
*          ICROW = INDXG2P( IC, MB_C, MYROW, RSRC_C, NPROW ),
*          ICCOL = INDXG2P( JC, NB_C, MYCOL, CSRC_C, NPCOL ),
*          MpC0 = NUMROC( M+IROFFC, MB_C, MYROW, ICROW, NPROW ),
*          NpC0 = NUMROC( N+ICOFFC, MB_C, MYROW, ICROW, NPROW ),
*          NqC0 = NUMROC( N+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL ),
*
*          ILCM, INDXG2P and NUMROC are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrices V(IV:*, JV:*) and C(IC:IC+M-1,JC:JC+N-1)
*  must verify some alignment properties, namely the following
*  expressions should be true:
*
*  If STOREV = 'Columnwise'
*    If SIDE = 'Left',
*      ( MB_V.EQ.MB_C .AND. IROFFV.EQ.IROFFC .AND. IVROW.EQ.ICROW )
*    If SIDE = 'Right',
*      ( MB_V.EQ.NB_C .AND. IROFFV.EQ.ICOFFC )
*  else if STOREV = 'Rowwise'
*    If SIDE = 'Left',
*      ( NB_V.EQ.MB_C .AND. ICOFFV.EQ.IROFFC )
*    If SIDE = 'Right',
*      ( NB_V.EQ.NB_C .AND. ICOFFV.EQ.ICOFFC .AND. IVCOL.EQ.ICCOL )
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
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                     ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT
      CHARACTER          COLBTOP, TRANST
      INTEGER            ICCOL1, ICCOL2, ICOFFC1, ICOFFC2, ICOFFV,
     $                   ICROW1, ICROW2, ICTXT, IIBEG, IIC1, IIC2,
     $                   IIEND, IINXT, IIV, ILEFT, INFO, IOFFC2, IOFFV,
     $                   IPT, IPV, IPW, IROFFC1, IROFFC2, ITOP, IVCOL,
     $                   IVROW, J, JJBEG, JJEND, JJNXT, JJC1, JJC2, JJV,
     $                   LDC, LDV, LV, LW, MBC, MBV, MPC1, MPC2, MPC20,
     $                   MQV, MQV0, MYCOL, MYDIST, MYROW, NBC, NBV,
     $                   NPCOL, NPROW, NQC1, NQC2, NQCALL, NQV
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GRIDINFO, CGEBR2D,
     $                   CGEBS2D, CGEMM, CGSUM2D, CLACGV,
     $                   CLAMOV, CLASET, CTRBR2D, CTRBS2D,
     $                   CTRMM, INFOG2L, PBCMATADD, PBCTRAN,
     $                   PB_TOPGET, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, NUMROC
      EXTERNAL           ICEIL, LSAME, NUMROC
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 .OR. K.LE.0 )
     $   RETURN
*
*     Get grid parameters
*
      ICTXT = DESCC( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Check for currently supported options
*
      INFO = 0
      IF( .NOT.LSAME( DIRECT, 'B' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( STOREV, 'R' ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCLARZB', -INFO )
         CALL BLACS_ABORT( ICTXT, 1 )
         RETURN
      END IF
*
      LEFT = LSAME( SIDE, 'L' )
      IF( LSAME( TRANS, 'N' ) ) THEN
          TRANST = 'C'
      ELSE
          TRANST = 'N'
      END IF
*
      CALL INFOG2L( IV, JV, DESCV, NPROW, NPCOL, MYROW, MYCOL, IIV, JJV,
     $              IVROW, IVCOL )
      MBV = DESCV( MB_ )
      NBV = DESCV( NB_ )
      ICOFFV = MOD( JV-1, NBV )
      NQV = NUMROC( L+ICOFFV, NBV, MYCOL, IVCOL, NPCOL )
      IF( MYCOL.EQ.IVCOL )
     $   NQV = NQV - ICOFFV
      LDV = DESCV( LLD_ )
      IIV = MIN( IIV, LDV )
      JJV = MIN( JJV, MAX( 1, NUMROC( DESCV( N_ ), NBV, MYCOL,
     $                                DESCV( CSRC_ ), NPCOL ) ) )
      IOFFV = IIV + ( JJV-1 ) * LDV
      MBC = DESCC( MB_ )
      NBC = DESCC( NB_ )
      NQCALL = NUMROC( DESCC( N_ ), NBC, MYCOL, DESCC( CSRC_ ), NPCOL )
      CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL, IIC1,
     $              JJC1, ICROW1, ICCOL1 )
      LDC = DESCC( LLD_ )
      IIC1 = MIN( IIC1, LDC )
      JJC1 = MIN( JJC1, MAX( 1, NQCALL ) )
*
      IF( LEFT ) THEN
         IROFFC1 = MOD( IC-1, MBC )
         MPC1 = NUMROC( K+IROFFC1, MBC, MYROW, ICROW1, NPROW )
         IF( MYROW.EQ.ICROW1 )
     $      MPC1 = MPC1 - IROFFC1
         ICOFFC1 = MOD( JC-1, NBC )
         NQC1 = NUMROC( N+ICOFFC1, NBC, MYCOL, ICCOL1, NPCOL )
         IF( MYCOL.EQ.ICCOL1 )
     $      NQC1 = NQC1 - ICOFFC1
         CALL INFOG2L( IC+M-L, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIC2, JJC2, ICROW2, ICCOL2 )
         IROFFC2 = MOD( IC+M-L-1, MBC )
         MPC2 = NUMROC( L+IROFFC2, MBC, MYROW, ICROW2, NPROW )
         IF( MYROW.EQ.ICROW2 )
     $      MPC2 = MPC2 - IROFFC2
         ICOFFC2 = ICOFFC1
         NQC2 = NQC1
      ELSE
         IROFFC1 = MOD( IC-1, MBC )
         MPC1 = NUMROC( M+IROFFC1, MBC, MYROW, ICROW1, NPROW )
         IF( MYROW.EQ.ICROW1 )
     $      MPC1 = MPC1 - IROFFC1
         ICOFFC1 = MOD( JC-1, NBC )
         NQC1 = NUMROC( K+ICOFFC1, NBC, MYCOL, ICCOL1, NPCOL )
         IF( MYCOL.EQ.ICCOL1 )
     $      NQC1 = NQC1 - ICOFFC1
         CALL INFOG2L( IC, JC+N-L, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIC2, JJC2, ICROW2, ICCOL2 )
         IROFFC2 = IROFFC1
         MPC2 = MPC1
         ICOFFC2 = MOD( JC+N-L-1, NBC )
         NQC2 = NUMROC( L+ICOFFC2, NBC, MYCOL, ICCOL2, NPCOL )
         IF( MYCOL.EQ.ICCOL2 )
     $      NQC2 = NQC2 - ICOFFC2
      END IF
      IIC2 = MIN( IIC2, LDC )
      JJC2 = MIN( JJC2, NQCALL )
      IOFFC2 = IIC2 + ( JJC2-1 ) * LDC
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form Q*sub( C ) or Q'*sub( C )
*
*        IROFFC2 = ICOFFV is required by the current transposition
*        routine PBCTRAN
*
         MQV0 = NUMROC( M+ICOFFV, NBV, MYCOL, IVCOL, NPCOL )
         IF( MYCOL.EQ.IVCOL ) THEN
            MQV = MQV0 - ICOFFV
         ELSE
            MQV = MQV0
         END IF
         IF( MYROW.EQ.ICROW2 ) THEN
            MPC20 = MPC2 + IROFFC2
         ELSE
            MPC20 = MPC2
         END IF
*
*        Locally V( IOFFV ) is K x MQV, C( IOFFC2 ) is MPC2 x NQC2
*        WORK( IPV ) is MPC20 x K = [ . V( IOFFV ) ]'
*        WORK( IPW ) is K x MQV0  = [ . V( IOFFV ) ]
*        WORK( IPT ) is the workspace for PBCTRAN
*
         IPV = 1
         IPW = IPV + MPC20 * K
         IPT = IPW + K * MQV0
         LV = MAX( 1, MPC20 )
         LW = MAX( 1, K )
*
         IF( MYROW.EQ.IVROW ) THEN
            IF( MYCOL.EQ.IVCOL ) THEN
               CALL CLAMOV( 'All', K, MQV, V( IOFFV ), LDV,
     $                      WORK( IPW+ICOFFV*LW ), LW )
            ELSE
               CALL CLAMOV( 'All', K, MQV, V( IOFFV ), LDV,
     $                      WORK( IPW ), LW )
            END IF
         END IF
*
*        WORK( IPV ) = WORK( IPW )' (replicated) is MPC20 x K
*
         CALL PBCTRAN( ICTXT, 'Rowwise', 'Conjugate transpose', K,
     $                 M+ICOFFV, DESCV( NB_ ), WORK( IPW ), LW, ZERO,
     $                 WORK( IPV ), LV, IVROW, IVCOL, ICROW2, -1,
     $                 WORK( IPT ) )
*
*        WORK( IPV ) = ( . V )' -> WORK( IPV ) = V' is MPC2 x K
*
         IF( MYROW.EQ.ICROW2 )
     $      IPV = IPV + IROFFC2
*
*        WORK( IPW ) becomes NQC2 x K = C( IOFFC2 )' * V'
*        WORK( IPW ) = C( IOFFC2 )' * V'  (NQC2 x MPC2 x K) -> NQC2 x K
*
         LW = MAX( 1, NQC2 )
*
         IF( MPC2.GT.0 ) THEN
            CALL CGEMM( 'Transpose', 'No transpose', NQC2, K, MPC2,
     $                  ONE, C( IOFFC2 ), LDC, WORK( IPV ), LV, ZERO,
     $                  WORK( IPW ), LW )
         ELSE
            CALL CLASET( 'All', NQC2, K, ZERO, ZERO, WORK( IPW ), LW )
         END IF
*
*        WORK( IPW ) = WORK( IPW ) + C1 ( NQC1 = NQC2 )
*
         IF( MPC1.GT.0 ) THEN
            MYDIST = MOD( MYROW-ICROW1+NPROW, NPROW )
            ITOP = MAX( 0, MYDIST * MBC - IROFFC1 )
            IIBEG = IIC1
            IIEND = IIC1 + MPC1 - 1
            IINXT = MIN( ICEIL( IIBEG, MBC ) * MBC, IIEND )
*
   10       CONTINUE
            IF( IIBEG.LE.IINXT ) THEN
               CALL PBCMATADD( ICTXT, 'Transpose', NQC2, IINXT-IIBEG+1,
     $                         ONE, C( IIBEG+(JJC1-1)*LDC ), LDC, ONE,
     $                         WORK( IPW+ITOP ), LW )
               MYDIST = MYDIST + NPROW
               ITOP = MYDIST * MBC - IROFFC1
               IIBEG = IINXT +1
               IINXT = MIN( IINXT+MBC, IIEND )
               GO TO 10
            END IF
         END IF
*
         CALL CGSUM2D( ICTXT, 'Columnwise', ' ', NQC2, K, WORK( IPW ),
     $                 LW, IVROW, MYCOL )
*
*        WORK( IPW ) = WORK( IPW ) * T' or WORK( IPW ) * T
*
         IF( MYROW.EQ.IVROW ) THEN
            IF( MYCOL.EQ.IVCOL ) THEN
*
*              Broadcast the block reflector to the other columns.
*
               CALL CTRBS2D( ICTXT, 'Rowwise', ' ', 'Lower', 'Non unit',
     $                       K, K, T, MBV )
            ELSE
               CALL CTRBR2D( ICTXT, 'Rowwise', ' ', 'Lower', 'Non unit',
     $                       K, K, T, MBV, MYROW, IVCOL )
            END IF
            CALL CTRMM( 'Right', 'Lower', TRANST, 'Non unit', NQC2, K,
     $                  ONE, T, MBV, WORK( IPW ), LW )
*
            CALL CGEBS2D( ICTXT, 'Columnwise', ' ', NQC2, K,
     $                    WORK( IPW ), LW )
         ELSE
            CALL CGEBR2D( ICTXT, 'Columnwise', ' ', NQC2, K,
     $                    WORK( IPW ), LW, IVROW, MYCOL )
         END IF
*
*        C1 = C1 - WORK( IPW )
*
         IF( MPC1.GT.0 ) THEN
            MYDIST = MOD( MYROW-ICROW1+NPROW, NPROW )
            ITOP = MAX( 0, MYDIST * MBC - IROFFC1 )
            IIBEG = IIC1
            IIEND = IIC1 + MPC1 - 1
            IINXT = MIN( ICEIL( IIBEG, MBC ) * MBC, IIEND )
*
   20       CONTINUE
            IF( IIBEG.LE.IINXT ) THEN
               CALL PBCMATADD( ICTXT, 'Transpose', IINXT-IIBEG+1, NQC2,
     $                         -ONE, WORK( IPW+ITOP ), LW, ONE,
     $                         C( IIBEG+(JJC1-1)*LDC ), LDC )
               MYDIST = MYDIST + NPROW
               ITOP = MYDIST * MBC - IROFFC1
               IIBEG = IINXT +1
               IINXT = MIN( IINXT+MBC, IIEND )
               GO TO 20
            END IF
         END IF
*
*            C2            C2      -     V'      *     W'
*        C( IOFFC2 ) = C( IOFFC2 ) - WORK( IPV ) * WORK( IPW )'
*                      MPC2 x NQC2    MPC2 x K      K x NQC2
*
         DO 30 J = 1, K
            CALL CLACGV( MPC2, WORK( IPV+(J-1)*LV ), 1 )
   30    CONTINUE
         CALL CGEMM( 'No transpose', 'Transpose', MPC2, NQC2, K, -ONE,
     $               WORK( IPV ), LV, WORK( IPW ), LW, ONE,
     $               C( IOFFC2 ), LDC )
*
      ELSE
*
*        Form sub( C ) * Q or sub( C ) * Q'
*
*        Locally V( IOFFV ) is K x NQV, C( IOFFC2 ) is MPC2 x NQC2
*        WORK( IPV ) is K x NQV = V( IOFFV ), NQV = NQC2
*        WORK( IPW ) is MPC2 x K = C( IOFFC2 ) * V( IOFFV )'
*
         IPV = 1
         IPW = IPV + K * NQC2
         LV = MAX( 1, K )
         LW = MAX( 1, MPC2 )
*
*        Broadcast V to the other process rows.
*
         CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
         IF( MYROW.EQ.IVROW ) THEN
            CALL CGEBS2D( ICTXT, 'Columnwise', COLBTOP, K, NQC2,
     $                    V( IOFFV ), LDV )
            IF( MYCOL.EQ.IVCOL )
     $         CALL CTRBS2D( ICTXT, 'Columnwise', COLBTOP, 'Lower',
     $                       'Non unit', K, K, T, MBV )
            CALL CLAMOV( 'All', K, NQC2, V( IOFFV ), LDV, WORK( IPV ),
     $                   LV )
         ELSE
            CALL CGEBR2D( ICTXT, 'Columnwise', COLBTOP, K, NQC2,
     $                    WORK( IPV ), LV, IVROW, MYCOL )
            IF( MYCOL.EQ.IVCOL )
     $         CALL CTRBR2D( ICTXT, 'Columnwise', COLBTOP, 'Lower',
     $                       'Non unit', K, K, T, MBV, IVROW, MYCOL )
         END IF
*
*        WORK( IPV ) is K x NQC2 = V = V( IOFFV )
*        WORK( IPW ) = C( IOFFC2 ) * V'  (MPC2 x NQC2 x K) -> MPC2 x K
*
         IF( NQC2.GT.0 ) THEN
            CALL CGEMM( 'No Transpose', 'Transpose', MPC2, K, NQC2,
     $                 ONE, C( IOFFC2 ), LDC, WORK( IPV ), LV, ZERO,
     $                 WORK( IPW ), LW )
         ELSE
            CALL CLASET( 'All', MPC2, K, ZERO, ZERO, WORK( IPW ), LW )
         END IF
*
*        WORK( IPW ) = WORK( IPW ) + C1 ( MPC1 = MPC2 )
*
         IF( NQC1.GT.0 ) THEN
            MYDIST = MOD( MYCOL-ICCOL1+NPCOL, NPCOL )
            ILEFT = MAX( 0, MYDIST * NBC - ICOFFC1 )
            JJBEG = JJC1
            JJEND = JJC1 + NQC1 - 1
            JJNXT = MIN( ICEIL( JJBEG, NBC ) * NBC, JJEND )
*
   40       CONTINUE
            IF( JJBEG.LE.JJNXT ) THEN
               CALL PBCMATADD( ICTXT, 'No transpose', MPC2,
     $                         JJNXT-JJBEG+1, ONE,
     $                         C( IIC1+(JJBEG-1)*LDC ), LDC, ONE,
     $                         WORK( IPW+ILEFT*LW ), LW )
               MYDIST = MYDIST + NPCOL
               ILEFT = MYDIST * NBC - ICOFFC1
               JJBEG = JJNXT +1
               JJNXT = MIN( JJNXT+NBC, JJEND )
               GO TO 40
            END IF
         END IF
*
         CALL CGSUM2D( ICTXT, 'Rowwise', ' ', MPC2, K, WORK( IPW ),
     $                 LW, MYROW, IVCOL )
*
*        WORK( IPW ) = WORK( IPW ) * T' or WORK( IPW ) * T
*
         IF( MYCOL.EQ.IVCOL ) THEN
            DO 50 J = 1, K
               CALL CLACGV( K-J+1, T( J+(J-1)*MBV ), 1 )
   50       CONTINUE
            CALL CTRMM( 'Right', 'Lower', TRANS, 'Non unit', MPC2, K,
     $                  ONE, T, MBV, WORK( IPW ), LW )
            CALL CGEBS2D( ICTXT, 'Rowwise', ' ', MPC2, K, WORK( IPW ),
     $                    LW )
            DO 60 J = 1, K
               CALL CLACGV( K-J+1, T( J+(J-1)*MBV ), 1 )
   60       CONTINUE
         ELSE
            CALL CGEBR2D( ICTXT, 'Rowwise', ' ', MPC2, K, WORK( IPW ),
     $                    LW, MYROW, IVCOL )
         END IF
*
*        C1 = C1 - WORK( IPW )
*
         IF( NQC1.GT.0 ) THEN
            MYDIST = MOD( MYCOL-ICCOL1+NPCOL, NPCOL )
            ILEFT = MAX( 0, MYDIST * NBC - ICOFFC1 )
            JJBEG = JJC1
            JJEND = JJC1 + NQC1 - 1
            JJNXT = MIN( ICEIL( JJBEG, NBC ) * NBC, JJEND )
*
   70       CONTINUE
            IF( JJBEG.LE.JJNXT ) THEN
               CALL PBCMATADD( ICTXT, 'No transpose', MPC2,
     $                         JJNXT-JJBEG+1, -ONE,
     $                         WORK( IPW+ILEFT*LW ), LW, ONE,
     $                         C( IIC1+(JJBEG-1)*LDC ), LDC )
               MYDIST = MYDIST + NPCOL
               ILEFT = MYDIST * NBC - ICOFFC1
               JJBEG = JJNXT +1
               JJNXT = MIN( JJNXT+NBC, JJEND )
               GO TO 70
            END IF
         END IF
*
*            C2           C2      -     W       *  conjg( V )
*        C( IOFFC ) = C( IOFFC )  - WORK( IPW ) * conjg( WORK( IPV ) )
*                     MPC2 x NQC2    MPC2 x K      K x NQC2
*
         DO 80 J = 1, NQC2
            CALL CLACGV( K, WORK( IPV+(J-1)*LV ), 1 )
   80    CONTINUE
         IF( IOFFC2.GT.0 )
     $      CALL CGEMM( 'No transpose', 'No transpose', MPC2, NQC2, K,
     $                  -ONE, WORK( IPW ), LW, WORK( IPV ), LV, ONE,
     $                  C( IOFFC2 ), LDC )
*
      END IF
*
      RETURN
*
*     End of PCLARZB
*
      END
