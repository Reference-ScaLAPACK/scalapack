      SUBROUTINE PCLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, IV,
     $                    JV, DESCV, T, C, IC, JC, DESCC, WORK )
*
*  -- ScaLAPACK auxiliary routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, DIRECT, STOREV
      INTEGER            IC, IV, JC, JV, K, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCC( * ), DESCV( * )
      COMPLEX            C( * ), T( * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLARFB applies a complex block reflector Q or its conjugate
*  transpose Q**H to a complex M-by-N distributed matrix sub( C )
*  denoting C(IC:IC+M-1,JC:JC+N-1), from the left or the right.
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
*          Indicates how Q is formed from a product of elementary
*          reflectors
*          = 'F': Q = H(1) H(2) . . . H(k) (Forward)
*          = 'B': Q = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (global input) CHARACTER
*          Indicates how the vectors which define the elementary
*          reflectors are stored:
*          = 'C': Columnwise
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
*  V       (local input) COMPLEX pointer into the local memory
*          to an array of dimension ( LLD_V, LOCc(JV+K-1) ) if
*          STOREV = 'C', ( LLD_V, LOCc(JV+M-1)) if STOREV = 'R' and
*          SIDE = 'L', ( LLD_V, LOCc(JV+N-1) ) if STOREV = 'R' and
*          SIDE = 'R'. It contains the local pieces of the distributed
*          vectors V representing the Householder transformation.
*          See further details.
*          If STOREV = 'C' and SIDE = 'L', LLD_V >= MAX(1,LOCr(IV+M-1));
*          if STOREV = 'C' and SIDE = 'R', LLD_V >= MAX(1,LOCr(IV+N-1));
*          if STOREV = 'R', LLD_V >= LOCr(IV+K-1).
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
*          if STOREV = 'R' and NB_V by NB_V if STOREV = 'C'. The trian-
*          gular matrix T in the representation of the block reflector.
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
      PARAMETER          ( ONE  = ( 1.0E+0, 0.0E+0 ),
     $                     ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            FORWARD
      CHARACTER          COLBTOP, ROWBTOP, TRANST, UPLO
      INTEGER            HEIGHT, IBASE, ICCOL, ICOFFC, ICOFFV, ICROW,
     $                   ICTXT, II, IIBEG, IIC, IIEND, IINXT, IIV,
     $                   ILASTCOL, ILASTROW, ILEFT, IOFF, IOFFC, IOFFV,
     $                   IPT, IPV, IPW, IPW1, IRIGHT, IROFFC, IROFFV,
     $                   ITOP, IVCOL, IVROW, JJ, JJBEG, JJC, JJEND,
     $                   JJNXT, JJV, KP, KQ, LDC, LDV, LV, LW, MBV, MPC,
     $                   MPC0, MQV, MQV0, MYCOL, MYDIST, MYROW, NBV,
     $                   NPV, NPV0, NPCOL, NPROW, NQC, NQC0, WIDE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CGEBR2D, CGEBS2D,CGEMM,
     $                   CGSUM2D, CLAMOV, CLASET, CTRBR2D,
     $                   CTRBS2D, CTRMM, INFOG1L, INFOG2L, PB_TOPGET,
     $                   PBCTRAN
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
      IF( LSAME( TRANS, 'N' ) ) THEN
          TRANST = 'C'
      ELSE
          TRANST = 'N'
      END IF
      FORWARD = LSAME( DIRECT, 'F' )
      IF( FORWARD ) THEN
         UPLO = 'U'
      ELSE
         UPLO = 'L'
      END IF
*
      CALL INFOG2L( IV, JV, DESCV, NPROW, NPCOL, MYROW, MYCOL, IIV, JJV,
     $              IVROW, IVCOL )
      CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL, IIC, JJC,
     $              ICROW, ICCOL )
      LDC = DESCC( LLD_ )
      LDV = DESCV( LLD_ )
      IIC = MIN( IIC, LDC )
      IIV = MIN( IIV, LDV )
      IROFFC = MOD( IC-1, DESCC( MB_ ) )
      ICOFFC = MOD( JC-1, DESCC( NB_ ) )
      MBV = DESCV( MB_ )
      NBV = DESCV( NB_ )
      IROFFV = MOD( IV-1, MBV )
      ICOFFV = MOD( JV-1, NBV )
      MPC = NUMROC( M+IROFFC, DESCC( MB_ ), MYROW, ICROW, NPROW )
      NQC = NUMROC( N+ICOFFC, DESCC( NB_ ), MYCOL, ICCOL, NPCOL )
      IF( MYCOL.EQ.ICCOL )
     $   NQC = NQC - ICOFFC
      IF( MYROW.EQ.ICROW )
     $   MPC = MPC - IROFFC
      JJC = MIN( JJC, MAX( 1, JJC+NQC-1 ) )
      JJV = MIN( JJV, MAX( 1, NUMROC( DESCV( N_ ), NBV, MYCOL,
     $                                DESCV( CSRC_ ), NPCOL ) ) )
      IOFFC = IIC + ( JJC-1 ) * LDC
      IOFFV = IIV + ( JJV-1 ) * LDV
*
      IF( LSAME( STOREV, 'C' ) ) THEN
*
*        V is stored columnwise
*
         IF( LSAME( SIDE, 'L' ) ) THEN
*
*           Form  Q*sub( C )  or  Q'*sub( C  )
*
*           Locally V( IOFFV ) is MPV x K, C( IOFFC ) is MPC x NQC
*           WORK( IPV ) is MPC x K = V( IOFFV ), MPC = MPV
*           WORK( IPW ) is NQC x K = C( IOFFC )' * V( IOFFV )
*
            IPV = 1
            IPW = IPV + MPC * K
            LV = MAX( 1, MPC )
            LW = MAX( 1, NQC )
*
*           Broadcast V to the other process columns.
*
            CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
            IF( MYCOL.EQ.IVCOL ) THEN
               CALL CGEBS2D( ICTXT, 'Rowwise', ROWBTOP, MPC, K,
     $                       V( IOFFV ), LDV )
               IF( MYROW.EQ.IVROW )
     $            CALL CTRBS2D( ICTXT, 'Rowwise', ROWBTOP, UPLO,
     $                          'Non unit', K, K, T, NBV )
               CALL CLAMOV( 'All', MPC, K, V( IOFFV ), LDV, WORK( IPV ),
     $                      LV )
            ELSE
               CALL CGEBR2D( ICTXT, 'Rowwise', ROWBTOP, MPC, K,
     $                       WORK( IPV ), LV, MYROW, IVCOL )
               IF( MYROW.EQ.IVROW )
     $            CALL CTRBR2D( ICTXT, 'Rowwise', ROWBTOP, UPLO,
     $                          'Non unit', K, K, T, NBV, MYROW, IVCOL )
            END IF
*
            IF( FORWARD ) THEN
*
*              WORK(IPV) = ( V1 ) where V1 is unit lower triangular,
*                          ( V2 ) zeroes upper triangular part of V1
*
               MYDIST = MOD( MYROW-IVROW+NPROW, NPROW )
               ITOP = MAX( 0, MYDIST*MBV - IROFFV )
               IIBEG = IIV
               IIEND = IIBEG + MPC - 1
               IINXT = MIN( ICEIL( IIBEG, MBV )*MBV, IIEND )
*
   10          CONTINUE
               IF( K-ITOP .GT.0 ) THEN
                  CALL CLASET( 'Upper', IINXT-IIBEG+1, K-ITOP, ZERO,
     $                         ONE, WORK( IPV+IIBEG-IIV+ITOP*LV ), LV )
                  MYDIST = MYDIST + NPROW
                  ITOP = MYDIST * MBV - IROFFV
                  IIBEG = IINXT + 1
                  IINXT = MIN( IINXT+MBV, IIEND )
                  GO TO 10
               END IF
*
            ELSE
*
*              WORK(IPV) = ( V1 ) where V2 is unit upper triangular,
*                          ( V2 ) zeroes lower triangular part of V2
*
               JJ = JJV
               IOFF = MOD( IV+M-K-1, MBV )
               CALL INFOG1L( IV+M-K, MBV, NPROW, MYROW, DESCV( RSRC_ ),
     $                       II, ILASTROW )
               KP = NUMROC( K+IOFF, MBV, MYROW, ILASTROW, NPROW )
               IF( MYROW.EQ.ILASTROW )
     $            KP = KP - IOFF
               MYDIST = MOD( MYROW-ILASTROW+NPROW, NPROW )
               ITOP = MYDIST * MBV - IOFF
               IBASE = MIN( ITOP+MBV, K )
               ITOP = MIN( MAX( 0, ITOP ), K )
*
   20          CONTINUE
               IF( JJ.LE.( JJV+K-1 ) ) THEN
                  HEIGHT = IBASE - ITOP
                  CALL CLASET( 'All', KP, ITOP-JJ+JJV, ZERO, ZERO,
     $                         WORK( IPV+II-IIV+(JJ-JJV)*LV ), LV )
                  CALL CLASET( 'Lower', KP, HEIGHT, ZERO, ONE,
     $                         WORK( IPV+II-IIV+ITOP*LV ), LV )
                  KP = MAX( 0, KP - HEIGHT )
                  II = II + HEIGHT
                  JJ = JJV + IBASE
                  MYDIST = MYDIST + NPROW
                  ITOP = MYDIST * MBV - IOFF
                  IBASE = MIN( ITOP + MBV, K )
                  ITOP = MIN( ITOP, K )
                  GO TO 20
               END IF
*
            END IF
*
*           WORK( IPW ) = C( IOFFC )' * V  (NQC x MPC x K) -> NQC x K
*
            IF( MPC.GT.0 ) THEN
               CALL CGEMM( 'Conjugate transpose', 'No transpose', NQC,
     $                    K, MPC, ONE, C( IOFFC ), LDC, WORK( IPV ), LV,
     $                    ZERO, WORK( IPW ), LW )
            ELSE
               CALL CLASET( 'All', NQC, K, ZERO, ZERO, WORK( IPW ), LW )
            END IF
*
            CALL CGSUM2D( ICTXT, 'Columnwise', ' ', NQC, K, WORK( IPW ),
     $                    LW, IVROW, MYCOL )
*
            IF( MYROW.EQ.IVROW ) THEN
*
*              WORK( IPW ) = WORK( IPW ) * T' or WORK( IPW ) * T
*
               CALL CTRMM( 'Right', UPLO, TRANST, 'Non unit', NQC, K,
     $                     ONE, T, NBV, WORK( IPW ), LW )
               CALL CGEBS2D( ICTXT, 'Columnwise', ' ', NQC, K,
     $                       WORK( IPW ), LW )
            ELSE
               CALL CGEBR2D( ICTXT, 'Columnwise', ' ', NQC, K,
     $                       WORK( IPW ), LW, IVROW, MYCOL )
            END IF
*
*               C            C      -     V       *     W'
*           C( IOFFC ) = C( IOFFC ) - WORK( IPV ) * WORK( IPW )'
*                        MPC x NQC    MPC x K         K x NQC
*
            CALL CGEMM( 'No transpose', 'Conjugate transpose', MPC, NQC,
     $                  K, -ONE, WORK( IPV ), LV, WORK( IPW ), LW, ONE,
     $                  C( IOFFC ), LDC )
*
         ELSE
*
*           Form sub( C )*Q or sub( C )*Q'
*
*           ICOFFC = IROFFV is required by the current transposition
*           routine PBCTRAN
*
            NPV0 = NUMROC( N+IROFFV, MBV, MYROW, IVROW, NPROW )
            IF( MYROW.EQ.IVROW ) THEN
               NPV = NPV0 - IROFFV
            ELSE
               NPV = NPV0
            END IF
            IF( MYCOL.EQ.ICCOL ) THEN
               NQC0 = NQC + ICOFFC
            ELSE
               NQC0 = NQC
            END IF
*
*           Locally V( IOFFV ) is NPV x K C( IOFFC ) is MPC x NQC
*           WORK( IPV ) is K x NQC0 = [ . V( IOFFV ) ]'
*           WORK( IPW ) is NPV0 x K = [ . V( IOFFV )' ]'
*           WORK( IPT ) is the workspace for PBCTRAN
*
            IPV = 1
            IPW = IPV + K * NQC0
            IPT = IPW + NPV0 * K
            LV = MAX( 1, K )
            LW = MAX( 1, NPV0 )
*
            IF( MYCOL.EQ.IVCOL ) THEN
               IF( MYROW.EQ.IVROW ) THEN
                  CALL CLASET( 'All', IROFFV, K, ZERO, ZERO,
     $                         WORK( IPW ), LW )
                  IPW1 = IPW + IROFFV
                  CALL CLAMOV( 'All', NPV, K, V( IOFFV ), LDV,
     $                         WORK( IPW1 ), LW )
               ELSE
                  IPW1 = IPW
                  CALL CLAMOV( 'All', NPV, K, V( IOFFV ), LDV,
     $                         WORK( IPW1 ), LW )
               END IF
*
               IF( FORWARD ) THEN
*
*                 WORK(IPW) = ( . V1' V2' )' where V1 is unit lower
*                 triangular, zeroes upper triangular part of V1
*
                  MYDIST = MOD( MYROW-IVROW+NPROW, NPROW )
                  ITOP = MAX( 0, MYDIST*MBV - IROFFV )
                  IIBEG = IIV
                  IIEND = IIBEG + NPV - 1
                  IINXT = MIN( ICEIL( IIBEG, MBV )*MBV, IIEND )
*
   30             CONTINUE
                  IF( ( K-ITOP ).GT.0 ) THEN
                     CALL CLASET( 'Upper', IINXT-IIBEG+1, K-ITOP, ZERO,
     $                            ONE, WORK( IPW1+IIBEG-IIV+ITOP*LW ),
     $                            LW )
                     MYDIST = MYDIST + NPROW
                     ITOP = MYDIST * MBV - IROFFV
                     IIBEG = IINXT + 1
                     IINXT = MIN( IINXT+MBV, IIEND )
                     GO TO 30
                  END IF
*
               ELSE
*
*                 WORK( IPW ) = ( . V1' V2' )' where V2 is unit upper
*                 triangular, zeroes lower triangular part of V2.
*
                  JJ = JJV
                  CALL INFOG1L( IV+N-K, MBV, NPROW, MYROW,
     $                          DESCV( RSRC_ ), II, ILASTROW )
                  IOFF = MOD( IV+N-K-1, MBV )
                  KP = NUMROC( K+IOFF, MBV, MYROW, ILASTROW, NPROW )
                  IF( MYROW.EQ.ILASTROW )
     $               KP = KP - IOFF
                  MYDIST = MOD( MYROW-ILASTROW+NPROW, NPROW )
                  ITOP = MYDIST * MBV - IOFF
                  IBASE = MIN( ITOP+MBV, K )
                  ITOP = MIN( MAX( 0, ITOP ), K )
*
   40             CONTINUE
                  IF( JJ.LE.( JJV+K-1 ) ) THEN
                     HEIGHT = IBASE - ITOP
                     CALL CLASET( 'All', KP, ITOP-JJ+JJV, ZERO, ZERO,
     $                            WORK( IPW1+II-IIV+(JJ-JJV)*LW ), LW )
                     CALL CLASET( 'Lower', KP, HEIGHT, ZERO, ONE,
     $                            WORK( IPW1+II-IIV+ITOP*LW ), LW )
                     KP = MAX( 0, KP - HEIGHT )
                     II = II + HEIGHT
                     JJ = JJV + IBASE
                     MYDIST = MYDIST + NPROW
                     ITOP = MYDIST * MBV - IOFF
                     IBASE = MIN( ITOP + MBV, K )
                     ITOP = MIN( ITOP, K )
                     GO TO 40
                  END IF
               END IF
            END IF
*
            CALL PBCTRAN( ICTXT, 'Columnwise', 'Conjugate transpose',
     $                    N+IROFFV, K, MBV, WORK( IPW ), LW, ZERO,
     $                    WORK( IPV ), LV, IVROW, IVCOL, -1, ICCOL,
     $                    WORK( IPT ) )
*
*           WORK( IPV ) = ( . V' ) -> WORK( IPV ) = V' is K x NQC
*
            IF( MYCOL.EQ.ICCOL )
     $         IPV = IPV + ICOFFC * LV
*
*           WORK( IPW ) becomes MPC x K = C( IOFFC ) * V
*           WORK( IPW ) = C( IOFFC ) * V  (MPC x NQC x K) -> MPC x K
*
            LW = MAX( 1, MPC )
*
            IF( NQC.GT.0 ) THEN
               CALL CGEMM( 'No transpose', 'Conjugate transpose', MPC,
     $                     K, NQC, ONE, C( IOFFC ), LDC, WORK( IPV ),
     $                     LV, ZERO, WORK( IPW ), LW )
            ELSE
               CALL CLASET( 'All', MPC, K, ZERO, ZERO, WORK( IPW ), LW )
            END IF
*
            CALL CGSUM2D( ICTXT, 'Rowwise', ' ', MPC, K, WORK( IPW ),
     $                    LW, MYROW, IVCOL )
*
*           WORK( IPW ) = WORK( IPW ) * T' or WORK( IPW ) * T
*
            IF( MYCOL.EQ.IVCOL ) THEN
               IF( MYROW.EQ.IVROW ) THEN
*
*                 Broadcast the block reflector to the other rows.
*
                  CALL CTRBS2D( ICTXT, 'Columnwise', ' ', UPLO,
     $                          'Non unit', K, K, T, NBV )
               ELSE
                  CALL CTRBR2D( ICTXT, 'Columnwise', ' ', UPLO,
     $                          'Non unit', K, K, T, NBV, IVROW, MYCOL )
               END IF
               CALL CTRMM( 'Right', UPLO, TRANS, 'Non unit', MPC, K,
     $                     ONE, T, NBV, WORK( IPW ), LW )
*
               CALL CGEBS2D( ICTXT, 'Rowwise', ' ', MPC, K, WORK( IPW ),
     $                       LW )
            ELSE
               CALL CGEBR2D( ICTXT, 'Rowwise', ' ', MPC, K, WORK( IPW ),
     $                       LW, MYROW, IVCOL )
            END IF
*
*               C            C      -     W       *     V'
*           C( IOFFC ) = C( IOFFC ) - WORK( IPW ) * WORK( IPV )
*                        MPC x NQC    MPC x K         K x NQC
*
            CALL CGEMM( 'No transpose', 'No transpose', MPC, NQC, K,
     $                  -ONE, WORK( IPW ), LW, WORK( IPV ), LV, ONE,
     $                  C( IOFFC ), LDC )
         END IF
*
      ELSE
*
*        V is stored rowwise
*
         IF( LSAME( SIDE, 'L' ) ) THEN
*
*           Form Q*sub( C ) or Q'*sub( C )
*
*           IROFFC = ICOFFV is required by the current transposition
*           routine PBCTRAN
*
            MQV0 = NUMROC( M+ICOFFV, NBV, MYCOL, IVCOL, NPCOL )
            IF( MYCOL.EQ.IVCOL ) THEN
               MQV = MQV0 - ICOFFV
            ELSE
               MQV = MQV0
            END IF
            IF( MYROW.EQ.ICROW ) THEN
               MPC0 = MPC + IROFFC
            ELSE
               MPC0 = MPC
            END IF
*
*           Locally V( IOFFV ) is K x MQV, C( IOFFC ) is MPC x NQC
*           WORK( IPV ) is MPC0 x K = [ . V( IOFFV ) ]'
*           WORK( IPW ) is K x MQV0 = [ . V( IOFFV ) ]
*           WORK( IPT ) is the workspace for PBCTRAN
*
            IPV = 1
            IPW = IPV + MPC0 * K
            IPT = IPW + K * MQV0
            LV = MAX( 1, MPC0 )
            LW = MAX( 1, K )
*
            IF( MYROW.EQ.IVROW ) THEN
               IF( MYCOL.EQ.IVCOL ) THEN
                  CALL CLASET( 'All', K, ICOFFV, ZERO, ZERO,
     $                         WORK( IPW ), LW )
                  IPW1 = IPW + ICOFFV * LW
                  CALL CLAMOV( 'All', K, MQV, V( IOFFV ), LDV,
     $                         WORK( IPW1 ), LW )
               ELSE
                  IPW1 = IPW
                  CALL CLAMOV( 'All', K, MQV, V( IOFFV ), LDV,
     $                         WORK( IPW1 ), LW )
               END IF
*
               IF( FORWARD ) THEN
*
*                 WORK( IPW ) = ( . V1 V2 ) where V1 is unit upper
*                 triangular, zeroes lower triangular part of V1
*
                  MYDIST = MOD( MYCOL-IVCOL+NPCOL, NPCOL )
                  ILEFT = MAX( 0, MYDIST * NBV - ICOFFV )
                  JJBEG = JJV
                  JJEND = JJV + MQV - 1
                  JJNXT = MIN( ICEIL( JJBEG, NBV ) * NBV, JJEND )
*
   50             CONTINUE
                  IF( ( K-ILEFT ).GT.0 ) THEN
                     CALL CLASET( 'Lower', K-ILEFT, JJNXT-JJBEG+1, ZERO,
     $                            ONE,
     $                            WORK( IPW1+ILEFT+(JJBEG-JJV)*LW ),
     $                            LW )
                     MYDIST = MYDIST + NPCOL
                     ILEFT = MYDIST * NBV - ICOFFV
                     JJBEG = JJNXT + 1
                     JJNXT = MIN( JJNXT+NBV, JJEND )
                     GO TO 50
                  END IF
*
               ELSE
*
*                 WORK( IPW ) = ( . V1 V2 ) where V2 is unit lower
*                 triangular, zeroes upper triangular part of V2.
*
                  II = IIV
                  CALL INFOG1L( JV+M-K, NBV, NPCOL, MYCOL,
     $                          DESCV( CSRC_ ), JJ, ILASTCOL )
                  IOFF = MOD( JV+M-K-1, NBV )
                  KQ = NUMROC( K+IOFF, NBV, MYCOL, ILASTCOL, NPCOL )
                  IF( MYCOL.EQ.ILASTCOL )
     $               KQ = KQ - IOFF
                  MYDIST = MOD( MYCOL-ILASTCOL+NPCOL, NPCOL )
                  ILEFT = MYDIST * NBV - IOFF
                  IRIGHT = MIN( ILEFT+NBV, K )
                  ILEFT = MIN( MAX( 0, ILEFT ), K )
*
   60             CONTINUE
                  IF( II.LE.( IIV+K-1 ) ) THEN
                     WIDE = IRIGHT - ILEFT
                     CALL CLASET( 'All', ILEFT-II+IIV, KQ, ZERO, ZERO,
     $                            WORK( IPW1+II-IIV+(JJ-JJV)*LW ), LW )
                     CALL CLASET( 'Upper', WIDE, KQ, ZERO, ONE,
     $                            WORK( IPW1+ILEFT+(JJ-JJV)*LW ), LW )
                     KQ = MAX( 0, KQ - WIDE )
                     II = IIV + IRIGHT
                     JJ = JJ + WIDE
                     MYDIST = MYDIST + NPCOL
                     ILEFT = MYDIST * NBV - IOFF
                     IRIGHT = MIN( ILEFT + NBV, K )
                     ILEFT = MIN( ILEFT, K )
                     GO TO 60
                  END IF
               END IF
            END IF
*
*           WORK( IPV ) = WORK( IPW )' (replicated) is MPC0 x K
*
            CALL PBCTRAN( ICTXT, 'Rowwise', 'Conjugate transpose', K,
     $                    M+ICOFFV, NBV, WORK( IPW ), LW, ZERO,
     $                    WORK( IPV ), LV, IVROW, IVCOL, ICROW, -1,
     $                    WORK( IPT ) )
*
*           WORK( IPV ) = ( . V )' -> WORK( IPV ) = V' is MPC x K
*
            IF( MYROW.EQ.ICROW )
     $         IPV = IPV + IROFFC
*
*           WORK( IPW ) becomes NQC x K = C( IOFFC )' * V'
*           WORK( IPW ) = C( IOFFC )' * V'  (NQC x MPC x K) -> NQC x K
*
            LW = MAX( 1, NQC )
*
            IF( MPC.GT.0 ) THEN
               CALL CGEMM( 'Conjugate transpose', 'No transpose', NQC,
     $                     K, MPC, ONE, C( IOFFC ), LDC, WORK( IPV ),
     $                     LV, ZERO, WORK( IPW ), LW )
            ELSE
               CALL CLASET( 'All', NQC, K, ZERO, ZERO, WORK( IPW ), LW )
            END IF
*
            CALL CGSUM2D( ICTXT, 'Columnwise', ' ', NQC, K, WORK( IPW ),
     $                    LW, IVROW, MYCOL )
*
*           WORK( IPW ) = WORK( IPW ) * T' or WORK( IPW ) * T
*
            IF( MYROW.EQ.IVROW ) THEN
               IF( MYCOL.EQ.IVCOL ) THEN
*
*                 Broadcast the block reflector to the other columns.
*
                  CALL CTRBS2D( ICTXT, 'Rowwise', ' ', UPLO, 'Non unit',
     $                          K, K, T, MBV )
               ELSE
                  CALL CTRBR2D( ICTXT, 'Rowwise', ' ', UPLO, 'Non unit',
     $                          K, K, T, MBV, MYROW, IVCOL )
               END IF
               CALL CTRMM( 'Right', UPLO, TRANST, 'Non unit', NQC, K,
     $                     ONE, T, MBV, WORK( IPW ), LW )
*
               CALL CGEBS2D( ICTXT, 'Columnwise', ' ', NQC, K,
     $                       WORK( IPW ), LW )
            ELSE
               CALL CGEBR2D( ICTXT, 'Columnwise', ' ', NQC, K,
     $                       WORK( IPW ), LW, IVROW, MYCOL )
            END IF
*
*               C            C      -     V'      *     W'
*           C( IOFFC ) = C( IOFFC ) - WORK( IPV ) * WORK( IPW )'
*                        MPC x NQC    MPC x K         K x NQC
*
            CALL CGEMM( 'No transpose', 'Conjugate transpose', MPC, NQC,
     $                  K, -ONE, WORK( IPV ), LV, WORK( IPW ), LW, ONE,
     $                  C( IOFFC ), LDC )
*
         ELSE
*
*           Form Q*sub( C ) or Q'*sub( C )
*
*           Locally V( IOFFV ) is K x NQV, C( IOFFC ) is MPC x NQC
*           WORK( IPV ) is K x NQV = V( IOFFV ), NQV = NQC
*           WORK( IPW ) is MPC x K = C( IOFFC ) * V( IOFFV )'
*
            IPV = 1
            IPW = IPV + K * NQC
            LV = MAX( 1, K )
            LW = MAX( 1, MPC )
*
*           Broadcast V to the other process rows.
*
            CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
            IF( MYROW.EQ.IVROW ) THEN
               CALL CGEBS2D( ICTXT, 'Columnwise', COLBTOP, K, NQC,
     $                       V( IOFFV ), LDV )
               IF( MYCOL.EQ.IVCOL )
     $            CALL CTRBS2D( ICTXT, 'Columnwise', COLBTOP, UPLO,
     $                          'Non unit', K, K, T, MBV )
               CALL CLAMOV( 'All', K, NQC, V( IOFFV ), LDV, WORK( IPV ),
     $                      LV )
            ELSE
               CALL CGEBR2D( ICTXT, 'Columnwise', COLBTOP, K, NQC,
     $                       WORK( IPV ), LV, IVROW, MYCOL )
               IF( MYCOL.EQ.IVCOL )
     $            CALL CTRBR2D( ICTXT, 'Columnwise', COLBTOP, UPLO,
     $                          'Non unit', K, K, T, MBV, IVROW, MYCOL )
            END IF
*
            IF( FORWARD ) THEN
*
*              WORK(IPW) = ( V1 V2 ) where V1 is unit upper
*              triangular, zeroes lower triangular part of V1
*
               MYDIST = MOD( MYCOL-IVCOL+NPCOL, NPCOL )
               ILEFT = MAX( 0, MYDIST * NBV - ICOFFV )
               JJBEG = JJV
               JJEND = JJV + NQC - 1
               JJNXT = MIN( ICEIL( JJBEG, NBV ) * NBV, JJEND )
*
   70          CONTINUE
               IF( ( K-ILEFT ).GT.0 ) THEN
                  CALL CLASET( 'Lower', K-ILEFT, JJNXT-JJBEG+1, ZERO,
     $                            ONE, WORK( IPV+ILEFT+(JJBEG-JJV)*LV ),
     $                            LV )
                  MYDIST = MYDIST + NPCOL
                  ILEFT = MYDIST * NBV - ICOFFV
                  JJBEG = JJNXT + 1
                  JJNXT = MIN( JJNXT+NBV, JJEND )
                  GO TO 70
               END IF
*
            ELSE
*
*              WORK( IPW ) = ( . V1 V2 ) where V2 is unit lower
*              triangular, zeroes upper triangular part of V2.
*
               II = IIV
               CALL INFOG1L( JV+N-K, NBV, NPCOL, MYCOL, DESCV( CSRC_ ),
     $                       JJ, ILASTCOL )
               IOFF = MOD( JV+N-K-1, NBV )
               KQ = NUMROC( K+IOFF, NBV, MYCOL, ILASTCOL, NPCOL )
               IF( MYCOL.EQ.ILASTCOL )
     $            KQ = KQ - IOFF
               MYDIST = MOD( MYCOL-ILASTCOL+NPCOL, NPCOL )
               ILEFT = MYDIST * NBV - IOFF
               IRIGHT = MIN( ILEFT+NBV, K )
               ILEFT = MIN( MAX( 0, ILEFT ), K )
*
   80          CONTINUE
               IF( II.LE.( IIV+K-1 ) ) THEN
                  WIDE = IRIGHT - ILEFT
                  CALL CLASET( 'All', ILEFT-II+IIV, KQ, ZERO, ZERO,
     $                         WORK( IPV+II-IIV+(JJ-JJV)*LV ), LV )
                  CALL CLASET( 'Upper', WIDE, KQ, ZERO, ONE,
     $                         WORK( IPV+ILEFT+(JJ-JJV)*LV ), LV )
                  KQ = MAX( 0, KQ - WIDE )
                  II = IIV + IRIGHT
                  JJ = JJ + WIDE
                  MYDIST = MYDIST + NPCOL
                  ILEFT = MYDIST * NBV - IOFF
                  IRIGHT = MIN( ILEFT + NBV, K )
                  ILEFT = MIN( ILEFT, K )
                  GO TO 80
               END IF
*
            END IF
*
*           WORK( IPV ) is K x NQC = V = V( IOFFV )
*           WORK( IPW ) = C( IOFFC ) * V'  (MPC x NQC x K) -> MPC x K
*
            IF( NQC.GT.0 ) THEN
               CALL CGEMM( 'No transpose', 'Conjugate transpose', MPC,
     $                    K, NQC, ONE, C( IOFFC ), LDC, WORK( IPV ),
     $                    LV, ZERO, WORK( IPW ), LW )
            ELSE
               CALL CLASET( 'All', MPC, K, ZERO, ZERO, WORK( IPW ), LW )
            END IF
*
            CALL CGSUM2D( ICTXT, 'Rowwise', ' ', MPC, K, WORK( IPW ),
     $                    LW, MYROW, IVCOL )
*
*           WORK( IPW ) = WORK( IPW ) * T' or WORK( IPW ) * T
*
            IF( MYCOL.EQ.IVCOL ) THEN
               CALL CTRMM( 'Right', UPLO, TRANS, 'Non unit', MPC, K,
     $                     ONE, T, MBV, WORK( IPW ), LW )
               CALL CGEBS2D( ICTXT, 'Rowwise', ' ', MPC, K, WORK( IPW ),
     $                       LW )
            ELSE
               CALL CGEBR2D( ICTXT, 'Rowwise', ' ', MPC, K, WORK( IPW ),
     $                       LW, MYROW, IVCOL )
            END IF
*
*               C            C      -     W       *     V
*           C( IOFFC ) = C( IOFFC ) - WORK( IPW ) * WORK( IPV )
*                        MPC x NQC    MPC x K         K x NQC
*
            CALL CGEMM( 'No transpose', 'No transpose', MPC, NQC, K,
     $                  -ONE, WORK( IPW ), LW, WORK( IPV ), LV, ONE,
     $                  C( IOFFC ), LDC )
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PCLARFB
*
      END
