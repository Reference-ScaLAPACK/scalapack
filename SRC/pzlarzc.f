      SUBROUTINE PZLARZC( SIDE, M, N, L, V, IV, JV, DESCV, INCV, TAU, C,
     $                    IC, JC, DESCC, WORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            IC, INCV, IV, JC, JV, L, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCC( * ), DESCV( * )
      COMPLEX*16         C( * ), TAU( * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLARZC applies a complex elementary reflector Q**H to a
*  complex M-by-N distributed matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1),
*  from either the left or the right. Q is represented in the form
*
*        Q = I - tau * v * v'
*
*  where tau is a complex scalar and v is a complex vector.
*
*  If tau = 0, then Q is taken to be the unit matrix.
*
*  Q is a product of k elementary reflectors as returned by PZTZRZF.
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
*  Because vectors may be viewed as a subclass of matrices, a
*  distributed vector is considered to be a distributed matrix.
*
*  Restrictions
*  ============
*
*  If SIDE = 'Left' and INCV = 1, then the row process having the first
*  entry V(IV,JV) must also own C(IC+M-L,JC:JC+N-1). Moreover,
*  MOD(IV-1,MB_V) must be equal to MOD(IC+N-L-1,MB_C), if INCV=M_V, only
*  the last equality must be satisfied.
*
*  If SIDE = 'Right' and INCV = M_V then the column process having the
*  first entry V(IV,JV) must also own C(IC:IC+M-1,JC+N-L) and
*  MOD(JV-1,NB_V) must be equal to MOD(JC+N-L-1,NB_C), if INCV = 1 only
*  the last equality must be satisfied.
*
*  Arguments
*  =========
*
*  SIDE    (global input) CHARACTER
*          = 'L': form  Q**H * sub( C ),
*          = 'R': form  sub( C ) * Q**H.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( C ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( C ). N >= 0.
*
*  L       (global input) INTEGER
*          The columns of the distributed submatrix sub( A ) containing
*          the meaningful part of the Householder reflectors.
*          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
*
*  V       (local input) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_V,*) containing the local
*          pieces of the distributed vectors V representing the
*          Householder transformation Q,
*             V(IV:IV+L-1,JV) if SIDE = 'L' and INCV = 1,
*             V(IV,JV:JV+L-1) if SIDE = 'L' and INCV = M_V,
*             V(IV:IV+L-1,JV) if SIDE = 'R' and INCV = 1,
*             V(IV,JV:JV+L-1) if SIDE = 'R' and INCV = M_V,
*
*          The vector v in the representation of Q. V is not used if
*          TAU = 0.
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
*  INCV    (global input) INTEGER
*          The global increment for the elements of V. Only two values
*          of INCV are supported in this version, namely 1 and M_V.
*          INCV must not be zero.
*
*  TAU     (local input) COMPLEX*16, array, dimension  LOCc(JV) if
*          INCV = 1, and LOCr(IV) otherwise. This array contains the
*          Householder scalars related to the Householder vectors.
*          TAU is tied to the distributed matrix V.
*
*  C       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_C, LOCc(JC+N-1) ),
*          containing the local pieces of sub( C ). On exit, sub( C )
*          is overwritten by the Q**H * sub( C ) if SIDE = 'L', or
*          sub( C ) * Q**H if SIDE = 'R'.
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
*  WORK    (local workspace) COMPLEX*16 array, dimension (LWORK)
*          If INCV = 1,
*            if SIDE = 'L',
*              if IVCOL = ICCOL,
*                LWORK >= NqC0
*              else
*                LWORK >= MpC0 + MAX( 1, NqC0 )
*              end if
*            else if SIDE = 'R',
*              LWORK >= NqC0 + MAX( MAX( 1, MpC0 ), NUMROC( NUMROC(
*                       N+ICOFFC,NB_V,0,0,NPCOL ),NB_V,0,0,LCMQ ) )
*            end if
*          else if INCV = M_V,
*            if SIDE = 'L',
*              LWORK >= MpC0 + MAX( MAX( 1, NqC0 ), NUMROC( NUMROC(
*                       M+IROFFC,MB_V,0,0,NPROW ),MB_V,0,0,LCMP ) )
*            else if SIDE = 'R',
*              if IVROW = ICROW,
*                LWORK >= MpC0
*              else
*                LWORK >= NqC0 + MAX( 1, MpC0 )
*              end if
*            end if
*          end if
*
*          where LCM is the least common multiple of NPROW and NPCOL and
*          LCM = ILCM( NPROW, NPCOL ), LCMP = LCM / NPROW,
*          LCMQ = LCM / NPCOL,
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
*  Alignment requirements
*  ======================
*
*  The distributed submatrices V(IV:*, JV:*) and C(IC:IC+M-1,JC:JC+N-1)
*  must verify some alignment properties, namely the following
*  expressions should be true:
*
*  MB_V = NB_V,
*
*  If INCV = 1,
*    If SIDE = 'Left',
*      ( MB_V.EQ.MB_C .AND. IROFFV.EQ.IROFFC .AND. IVROW.EQ.ICROW )
*    If SIDE = 'Right',
*      ( MB_V.EQ.NB_A .AND. MB_V.EQ.NB_C .AND. IROFFV.EQ.ICOFFC )
*  else if INCV = M_V,
*    If SIDE = 'Left',
*      ( MB_V.EQ.NB_V .AND. MB_V.EQ.MB_C .AND. ICOFFV.EQ.IROFFC )
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
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE  = ( 1.0D+0, 0.0D+0 ),
     $                     ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            CCBLCK, CRBLCK, LEFT
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            ICCOL1, ICCOL2, ICOFFC1, ICOFFC2, ICOFFV,
     $                   ICROW1, ICROW2, ICTXT, IIC1, IIC2, IIV, IOFFC1,
     $                   IOFFC2, IOFFV, IPW, IROFFC1, IROFFC2, IROFFV,
     $                   IVCOL, IVROW, JJC1, JJC2, JJV, LDC, LDV, MPC2,
     $                   MPV, MYCOL, MYROW, NCC, NCV, NPCOL, NPROW,
     $                   NQC2, NQV, RDEST
      COMPLEX*16         TAULOC( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, PB_TOPGET, PBZTRNV,
     $                   ZAXPY, ZCOPY, ZGEBR2D, ZGEBS2D,
     $                   ZGEMV, ZGERC, ZGERV2D, ZGESD2D,
     $                   ZGSUM2D, ZLASET
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      EXTERNAL           LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
*     Get grid parameters.
*
      ICTXT = DESCC( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Figure local indexes
*
      LEFT = LSAME( SIDE, 'L' )
      CALL INFOG2L( IV, JV, DESCV, NPROW, NPCOL, MYROW, MYCOL, IIV, JJV,
     $              IVROW, IVCOL )
      IROFFV = MOD( IV-1, DESCV( NB_ ) )
      MPV = NUMROC( L+IROFFV, DESCV( MB_ ), MYROW, IVROW, NPROW )
      IF( MYROW.EQ.IVROW )
     $   MPV = MPV - IROFFV
      ICOFFV = MOD( JV-1, DESCV( NB_ ) )
      NQV = NUMROC( L+ICOFFV, DESCV( NB_ ), MYCOL, IVCOL, NPCOL )
      IF( MYCOL.EQ.IVCOL )
     $   NQV = NQV - ICOFFV
      LDV = DESCV( LLD_ )
      NCV = NUMROC( DESCV( N_ ), DESCV( NB_ ), MYCOL, DESCV( CSRC_ ),
     $              NPCOL )
      LDV = DESCV( LLD_ )
      IIV = MIN( IIV, LDV )
      JJV = MIN( JJV, NCV )
      IOFFV = IIV+(JJV-1)*LDV
      NCC = NUMROC( DESCC( N_ ), DESCC( NB_ ), MYCOL, DESCC( CSRC_ ),
     $              NPCOL )
      CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $              IIC1, JJC1, ICROW1, ICCOL1 )
      IROFFC1 = MOD( IC-1, DESCC( MB_ ) )
      ICOFFC1 = MOD( JC-1, DESCC( NB_ ) )
      LDC = DESCC( LLD_ )
      IIC1 = MIN( IIC1, LDC )
      JJC1 = MIN( JJC1, MAX( 1, NCC ) )
      IOFFC1 = IIC1 + ( JJC1-1 ) * LDC
*
      IF( LEFT ) THEN
         CALL INFOG2L( IC+M-L, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIC2, JJC2, ICROW2, ICCOL2 )
         IROFFC2 = MOD( IC+M-L-1, DESCC( MB_ ) )
         ICOFFC2 = MOD( JC-1, DESCC( NB_ ) )
         NQC2 = NUMROC( N+ICOFFC2, DESCC( NB_ ), MYCOL, ICCOL2, NPCOL )
         IF( MYCOL.EQ.ICCOL2 )
     $      NQC2 = NQC2 - ICOFFC2
      ELSE
         CALL INFOG2L( IC, JC+N-L, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIC2, JJC2, ICROW2, ICCOL2 )
         IROFFC2 = MOD( IC-1, DESCC( MB_ ) )
         MPC2 = NUMROC( M+IROFFC2, DESCC( MB_ ), MYROW, ICROW2, NPROW )
         IF( MYROW.EQ.ICROW2 )
     $      MPC2 = MPC2 - IROFFC2
         ICOFFC2 = MOD( JC+N-L-1, DESCC( NB_ ) )
      END IF
      IIC2 = MIN( IIC2, LDC )
      JJC2 = MIN( JJC2, NCC )
      IOFFC2 = IIC2 + ( JJC2-1 ) * LDC
*
*     Is sub( C ) only distributed over a process row ?
*
      CRBLCK = ( M.LE.(DESCC( MB_ )-IROFFC1) )
*
*     Is sub( C ) only distributed over a process column ?
*
      CCBLCK = ( N.LE.(DESCC( NB_ )-ICOFFC1) )
*
      IF( LEFT ) THEN
*
         IF( CRBLCK ) THEN
            RDEST = ICROW2
         ELSE
            RDEST = -1
         END IF
*
         IF( CCBLCK ) THEN
*
*           sub( C ) is distributed over a process column
*
            IF( DESCV( M_ ).EQ.INCV ) THEN
*
*              Transpose row vector V (ICOFFV = IROFFC2)
*
               IPW = MPV+1
               CALL PBZTRNV( ICTXT, 'Rowwise', 'Transpose', M,
     $                       DESCV( NB_ ), IROFFC2, V( IOFFV ), LDV,
     $                       ZERO,
     $                       WORK, 1, IVROW, IVCOL, ICROW2, ICCOL2,
     $                       WORK( IPW ) )
*
*              Perform the local computation within a process column
*
               IF( MYCOL.EQ.ICCOL2 ) THEN
*
                  IF( MYROW.EQ.IVROW ) THEN
*
                     CALL ZGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1,
     $                             TAU( IIV ), 1 )
                     TAULOC( 1 ) = DCONJG( TAU( IIV ) )
*
                  ELSE
*
                     CALL ZGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1,
     $                             TAULOC, 1, IVROW, MYCOL )
                     TAULOC( 1 ) = DCONJG( TAULOC( 1 ) )
*
                  END IF
*
                  IF( TAULOC( 1 ).NE.ZERO ) THEN
*
*                    w := sub( C )' * v
*
                     IF( MPV.GT.0 ) THEN
                        CALL ZGEMV( 'Conjugate transpose', MPV, NQC2,
     $                              ONE, C( IOFFC2 ), LDC, WORK, 1,
     $                              ZERO, WORK( IPW ), 1 )
                     ELSE
                        CALL ZLASET( 'All', NQC2, 1, ZERO, ZERO,
     $                               WORK( IPW ), MAX( 1, NQC2 ) )
                     END IF
                     IF( MYROW.EQ.ICROW1 )
     $                  CALL ZAXPY( NQC2, ONE, C( IOFFC1 ), LDC,
     $                              WORK( IPW ), MAX( 1, NQC2 ) )
*
                     CALL ZGSUM2D( ICTXT, 'Columnwise', ' ', NQC2, 1,
     $                             WORK( IPW ), MAX( 1, NQC2 ), RDEST,
     $                             MYCOL )
*
*                    sub( C ) := sub( C ) - v * w'
*
                     IF( MYROW.EQ.ICROW1 )
     $                  CALL ZAXPY( NQC2, -TAULOC( 1 ), WORK( IPW ),
     $                              MAX( 1, NQC2 ), C( IOFFC1 ), LDC )
                     CALL ZGERC( MPV, NQC2, -TAULOC( 1 ), WORK, 1,
     $                           WORK( IPW ), 1, C( IOFFC2 ), LDC )
                  END IF
*
               END IF
*
            ELSE
*
*              V is a column vector
*
               IF( IVCOL.EQ.ICCOL2 ) THEN
*
*                 Perform the local computation within a process column
*
                  IF( MYCOL.EQ.ICCOL2 ) THEN
*
                     TAULOC( 1 ) = DCONJG( TAU( JJV ) )
*
                     IF( TAULOC( 1 ).NE.ZERO ) THEN
*
*                       w := sub( C )' * v
*
                        IF( MPV.GT.0 ) THEN
                           CALL ZGEMV( 'Conjugate transpose', MPV, NQC2,
     $                              ONE, C( IOFFC2 ), LDC, V( IOFFV ),
     $                              1, ZERO, WORK, 1 )
                        ELSE
                           CALL ZLASET( 'All', NQC2, 1, ZERO, ZERO,
     $                                  WORK, MAX( 1, NQC2 ) )
                        END IF
                        IF( MYROW.EQ.ICROW1 )
     $                     CALL ZAXPY( NQC2, ONE, C( IOFFC1 ), LDC,
     $                                 WORK, MAX( 1, NQC2 ) )
*
                        CALL ZGSUM2D( ICTXT, 'Columnwise', ' ', NQC2, 1,
     $                                WORK, MAX( 1, NQC2 ), RDEST,
     $                                MYCOL )
*
*                       sub( C ) := sub( C ) - v * w'
*
                        IF( MYROW.EQ.ICROW1 )
     $                     CALL ZAXPY( NQC2, -TAULOC( 1 ), WORK,
     $                                 MAX( 1, NQC2 ), C( IOFFC1 ),
     $                                 LDC )
                        CALL ZGERC( MPV, NQC2, -TAULOC( 1 ), V( IOFFV ),
     $                              1, WORK, 1, C( IOFFC2 ), LDC )
                     END IF
*
                  END IF
*
               ELSE
*
*                 Send V and TAU to the process column ICCOL2
*
                  IF( MYCOL.EQ.IVCOL ) THEN
*
                     IPW = MPV+1
                     CALL ZCOPY( MPV, V( IOFFV ), 1, WORK, 1 )
                     WORK( IPW ) = TAU( JJV )
                     CALL ZGESD2D( ICTXT, IPW, 1, WORK, IPW, MYROW,
     $                             ICCOL2 )
*
                  ELSE IF( MYCOL.EQ.ICCOL2 ) THEN
*
                     IPW = MPV+1
                     CALL ZGERV2D( ICTXT, IPW, 1, WORK, IPW, MYROW,
     $                             IVCOL )
                     TAULOC( 1 ) = DCONJG( WORK( IPW ) )
*
                     IF( TAULOC( 1 ).NE.ZERO ) THEN
*
*                       w := sub( C )' * v
*
                        IF( MPV.GT.0 ) THEN
                           CALL ZGEMV( 'Conjugate transpose', MPV, NQC2,
     $                                 ONE, C( IOFFC2 ), LDC, WORK, 1,
     $                                 ZERO, WORK( IPW ), 1 )
                        ELSE
                           CALL ZLASET( 'All', NQC2, 1, ZERO, ZERO,
     $                                  WORK( IPW ), MAX( 1, NQC2 ) )
                        END IF
                        IF( MYROW.EQ.ICROW1 )
     $                     CALL ZAXPY( NQC2, ONE, C( IOFFC1 ), LDC,
     $                                 WORK( IPW ), MAX( 1, NQC2 ) )
*
                        CALL ZGSUM2D( ICTXT, 'Columnwise', ' ', NQC2, 1,
     $                                WORK( IPW ), MAX( 1, NQC2 ),
     $                                RDEST, MYCOL )
*
*                       sub( C ) := sub( C ) - v * w'
*
                        IF( MYROW.EQ.ICROW1 )
     $                     CALL ZAXPY( NQC2, -TAULOC( 1 ), WORK( IPW ),
     $                                 MAX( 1, NQC2 ), C( IOFFC1 ),
     $                                 LDC )
                        CALL ZGERC( MPV, NQC2, -TAULOC( 1 ), WORK, 1,
     $                              WORK( IPW ), 1, C( IOFFC2 ), LDC )
                     END IF
*
                  END IF
*
               END IF
*
            END IF
*
         ELSE
*
*           sub( C ) is a proper distributed matrix
*
            IF( DESCV( M_ ).EQ.INCV ) THEN
*
*              Transpose and broadcast row vector V (ICOFFV=IROFFC2)
*
               IPW = MPV+1
               CALL PBZTRNV( ICTXT, 'Rowwise', 'Transpose', M,
     $                       DESCV( NB_ ), IROFFC2, V( IOFFV ), LDV,
     $                       ZERO,
     $                       WORK, 1, IVROW, IVCOL, ICROW2, -1,
     $                       WORK( IPW ) )
*
*              Perform the local computation within a process column
*
               IF( MYROW.EQ.IVROW ) THEN
*
                  CALL ZGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1,
     $                          TAU( IIV ), 1 )
                  TAULOC( 1 ) = DCONJG( TAU( IIV ) )
*
               ELSE
*
                  CALL ZGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1, TAULOC,
     $                          1, IVROW, MYCOL )
                  TAULOC( 1 ) = DCONJG( TAULOC( 1 ) )
*
               END IF
*
               IF( TAULOC( 1 ).NE.ZERO ) THEN
*
*                 w := sub( C )' * v
*
                  IF( MPV.GT.0 ) THEN
                     CALL ZGEMV( 'Conjugate transpose', MPV, NQC2, ONE,
     $                           C( IOFFC2 ), LDC, WORK, 1, ZERO,
     $                           WORK( IPW ), 1 )
                  ELSE
                     CALL ZLASET( 'All', NQC2, 1, ZERO, ZERO,
     $                            WORK( IPW ), MAX( 1, NQC2 ) )
                  END IF
                  IF( MYROW.EQ.ICROW1 )
     $               CALL ZAXPY( NQC2, ONE, C( IOFFC1 ), LDC,
     $                           WORK( IPW ), MAX( 1, NQC2 ) )
*
                  CALL ZGSUM2D( ICTXT, 'Columnwise', ' ', NQC2, 1,
     $                          WORK( IPW ), MAX( 1, NQC2 ), RDEST,
     $                          MYCOL )
*
*                 sub( C ) := sub( C ) - v * w'
*
                  IF( MYROW.EQ.ICROW1 )
     $               CALL ZAXPY( NQC2, -TAULOC( 1 ), WORK( IPW ),
     $                           MAX( 1, NQC2 ), C( IOFFC1 ), LDC )
                  CALL ZGERC( MPV, NQC2, -TAULOC( 1 ), WORK, 1,
     $                        WORK( IPW ), 1, C( IOFFC2 ), LDC )
               END IF
*
            ELSE
*
*              Broadcast column vector V
*
               CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
               IF( MYCOL.EQ.IVCOL ) THEN
*
                  IPW = MPV+1
                  CALL ZCOPY( MPV, V( IOFFV ), 1, WORK, 1 )
                  WORK( IPW ) = TAU( JJV )
                  CALL ZGEBS2D( ICTXT, 'Rowwise', ROWBTOP, IPW, 1,
     $                          WORK, IPW )
                  TAULOC( 1 ) = DCONJG( TAU( JJV ) )
*
               ELSE
*
                  IPW = MPV+1
                  CALL ZGEBR2D( ICTXT, 'Rowwise', ROWBTOP, IPW, 1, WORK,
     $                          IPW, MYROW, IVCOL )
                  TAULOC( 1 ) = DCONJG( WORK( IPW ) )
*
               END IF
*
               IF( TAULOC( 1 ).NE.ZERO ) THEN
*
*                 w := sub( C )' * v
*
                  IF( MPV.GT.0 ) THEN
                     CALL ZGEMV( 'Conjugate transpose', MPV, NQC2, ONE,
     $                           C( IOFFC2 ), LDC, WORK, 1, ZERO,
     $                           WORK( IPW ), 1 )
                  ELSE
                     CALL ZLASET( 'All', NQC2, 1, ZERO, ZERO,
     $                            WORK( IPW ), MAX( 1, NQC2 ) )
                  END IF
                  IF( MYROW.EQ.ICROW1 )
     $               CALL ZAXPY( NQC2, ONE, C( IOFFC1 ), LDC,
     $                           WORK( IPW ), MAX( 1, NQC2 ) )
*
                  CALL ZGSUM2D( ICTXT, 'Columnwise', ' ', NQC2, 1,
     $                          WORK( IPW ), MAX( 1, NQC2 ), RDEST,
     $                          MYCOL )
*
*                 sub( C ) := sub( C ) - v * w'
*
                  IF( MYROW.EQ.ICROW1 )
     $               CALL ZAXPY( NQC2, -TAULOC( 1 ), WORK( IPW ),
     $                           MAX( 1, NQC2 ), C( IOFFC1 ), LDC )
                  CALL ZGERC( MPV, NQC2, -TAULOC( 1 ), WORK, 1,
     $                        WORK( IPW ), 1, C( IOFFC2 ), LDC )
               END IF
*
            END IF
*
         END IF
*
      ELSE
*
         IF( CCBLCK ) THEN
            RDEST = MYROW
         ELSE
            RDEST = -1
         END IF
*
         IF( CRBLCK ) THEN
*
*           sub( C ) is distributed over a process row
*
            IF( DESCV( M_ ).EQ.INCV ) THEN
*
*              V is a row vector
*
               IF( IVROW.EQ.ICROW2 ) THEN
*
*                 Perform the local computation within a process row
*
                  IF( MYROW.EQ.ICROW2 ) THEN
*
                     TAULOC( 1 ) = DCONJG( TAU( IIV ) )
*
                     IF( TAULOC( 1 ).NE.ZERO ) THEN
*
*                       w := sub( C ) * v
*
                        IF( NQV.GT.0 ) THEN
                           CALL ZGEMV( 'No transpose', MPC2, NQV, ONE,
     $                                 C( IOFFC2 ), LDC, V( IOFFV ),
     $                                 LDV, ZERO, WORK, 1 )
                        ELSE
                           CALL ZLASET( 'All', MPC2, 1, ZERO, ZERO,
     $                                  WORK, MAX( 1, MPC2 ) )
                        END IF
                        IF( MYCOL.EQ.ICCOL1 )
     $                     CALL ZAXPY( MPC2, ONE, C( IOFFC1 ), 1,
     $                                   WORK, 1 )
*
                        CALL ZGSUM2D( ICTXT, 'Rowwise', ' ', MPC2, 1,
     $                                WORK, MAX( 1, MPC2 ), RDEST,
     $                               ICCOL2 )
*
                        IF( MYCOL.EQ.ICCOL1 )
     $                     CALL ZAXPY( MPC2, -TAULOC( 1 ), WORK, 1,
     $                                 C( IOFFC1 ), 1 )
*
*                       sub( C ) := sub( C ) - w * v'
*
                        CALL ZGERC( MPC2, NQV, -TAULOC( 1 ), WORK, 1,
     $                              V( IOFFV ), LDV, C( IOFFC2 ), LDC )
                     END IF
*
                  END IF
*
               ELSE
*
*                 Send V and TAU to the process row ICROW2
*
                  IF( MYROW.EQ.IVROW ) THEN
*
                     IPW = NQV+1
                     CALL ZCOPY( NQV, V( IOFFV ), LDV, WORK, 1 )
                     WORK( IPW ) = TAU( IIV )
                     CALL ZGESD2D( ICTXT, IPW, 1, WORK, IPW, ICROW2,
     $                             MYCOL )
*
                  ELSE IF( MYROW.EQ.ICROW2 ) THEN
*
                     IPW = NQV+1
                     CALL ZGERV2D( ICTXT, IPW, 1, WORK, IPW, IVROW,
     $                             MYCOL )
                     TAULOC( 1 ) = DCONJG( WORK( IPW ) )
*
                     IF( TAULOC( 1 ).NE.ZERO ) THEN
*
*                       w := sub( C ) * v
*
                        IF( NQV.GT.0 ) THEN
                           CALL ZGEMV( 'No transpose', MPC2, NQV, ONE,
     $                                 C( IOFFC2 ), LDC, WORK, 1, ZERO,
     $                                 WORK( IPW ), 1 )
                        ELSE
                           CALL ZLASET( 'All', MPC2, 1, ZERO, ZERO,
     $                                  WORK( IPW ), MAX( 1, MPC2 ) )
                        END IF
                        IF( MYCOL.EQ.ICCOL1 )
     $                     CALL ZAXPY( MPC2, ONE, C( IOFFC1 ), 1,
     $                                   WORK( IPW ), 1 )
                        CALL ZGSUM2D( ICTXT, 'Rowwise', ' ', MPC2, 1,
     $                                WORK( IPW ), MAX( 1, MPC2 ),
     $                                RDEST, ICCOL2 )
                        IF( MYCOL.EQ.ICCOL1 )
     $                     CALL ZAXPY( MPC2, -TAULOC( 1 ), WORK( IPW ),
     $                                 1, C( IOFFC1 ), 1 )
*
*                       sub( C ) := sub( C ) - w * v'
*
                        CALL ZGERC( MPC2, NQV, -TAULOC( 1 ),
     $                              WORK( IPW ), 1, WORK, 1,
     $                              C( IOFFC2 ), LDC )
                     END IF
*
                  END IF
*
               END IF
*
            ELSE
*
*              Transpose column vector V (IROFFV = ICOFFC2)
*
               IPW = NQV+1
               CALL PBZTRNV( ICTXT, 'Columnwise', 'Transpose', N,
     $                       DESCV( MB_ ), ICOFFC2, V( IOFFV ), 1, ZERO,
     $                       WORK, 1, IVROW, IVCOL, ICROW2, ICCOL2,
     $                       WORK( IPW ) )
*
*              Perform the local computation within a process column
*
               IF( MYROW.EQ.ICROW2 ) THEN
*
                  IF( MYCOL.EQ.IVCOL ) THEN
*
                     CALL ZGEBS2D( ICTXT, 'Rowwise', ' ', 1, 1,
     $                             TAU( JJV ), 1 )
                     TAULOC( 1 ) = DCONJG( TAU( JJV ) )
*
                  ELSE
*
                     CALL ZGEBR2D( ICTXT, 'Rowwise', ' ', 1, 1, TAULOC,
     $                             1, MYROW, IVCOL )
                     TAULOC( 1 ) = DCONJG( TAULOC( 1 ) )
*
                  END IF
*
                  IF( TAULOC( 1 ).NE.ZERO ) THEN
*
*                    w := sub( C ) * v
*
                     IF( NQV.GT.0 ) THEN
                        CALL ZGEMV( 'No transpose', MPC2, NQV, ONE,
     $                              C( IOFFC2 ), LDC, WORK, 1, ZERO,
     $                              WORK( IPW ), 1 )
                     ELSE
                        CALL ZLASET( 'All', MPC2, 1, ZERO, ZERO,
     $                               WORK( IPW ), MAX( 1, MPC2 ) )
                     END IF
                     IF( MYCOL.EQ.ICCOL1 )
     $                  CALL ZAXPY( MPC2, ONE, C( IOFFC1 ), 1,
     $                              WORK( IPW ), 1 )
                     CALL ZGSUM2D( ICTXT, 'Rowwise', ' ', MPC2, 1,
     $                             WORK( IPW ), MAX( 1, MPC2 ), RDEST,
     $                             ICCOL2 )
                     IF( MYCOL.EQ.ICCOL1 )
     $                  CALL ZAXPY( MPC2, -TAULOC( 1 ), WORK( IPW ), 1,
     $                              C( IOFFC1 ), 1 )
*
*                    sub( C ) := sub( C ) - w * v'
*
                     CALL ZGERC( MPC2, NQV, -TAULOC( 1 ), WORK( IPW ),
     $                           1, WORK, 1, C( IOFFC2 ), LDC )
                  END IF
*
               END IF
*
            END IF
*
         ELSE
*
*           sub( C ) is a proper distributed matrix
*
            IF( DESCV( M_ ).EQ.INCV ) THEN
*
*              Broadcast row vector V
*
               CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise',
     $                         COLBTOP )
               IF( MYROW.EQ.IVROW ) THEN
*
                  IPW = NQV+1
                  CALL ZCOPY( NQV, V( IOFFV ), LDV, WORK, 1 )
                  WORK( IPW ) = TAU( IIV )
                  CALL ZGEBS2D( ICTXT, 'Columnwise', COLBTOP, IPW, 1,
     $                          WORK, IPW )
                  TAULOC( 1 ) = DCONJG( TAU( IIV ) )
*
               ELSE
*
                  IPW = NQV+1
                  CALL ZGEBR2D( ICTXT, 'Columnwise', COLBTOP, IPW, 1,
     $                          WORK, IPW, IVROW, MYCOL )
                  TAULOC( 1 ) = DCONJG( WORK( IPW ) )
*
               END IF
*
               IF( TAULOC( 1 ).NE.ZERO ) THEN
*
*                 w := sub( C ) * v
*
                  IF( NQV.GT.0 ) THEN
                     CALL ZGEMV( 'No Transpose', MPC2, NQV, ONE,
     $                           C( IOFFC2 ), LDC, WORK, 1, ZERO,
     $                           WORK( IPW ), 1 )
                  ELSE
                     CALL ZLASET( 'All', MPC2, 1, ZERO, ZERO,
     $                            WORK( IPW ), MAX( 1, MPC2 ) )
                  END IF
                  IF( MYCOL.EQ.ICCOL1 )
     $               CALL ZAXPY( MPC2, ONE, C( IOFFC1 ), 1,
     $                           WORK( IPW ), 1 )
*
                  CALL ZGSUM2D( ICTXT, 'Rowwise', ' ', MPC2, 1,
     $                          WORK( IPW ), MAX( 1, MPC2 ), RDEST,
     $                          ICCOL2 )
                  IF( MYCOL.EQ.ICCOL1 )
     $               CALL ZAXPY( MPC2, -TAULOC( 1 ), WORK( IPW ), 1,
     $                           C( IOFFC1 ), 1 )
*
*                 sub( C ) := sub( C ) - w * v'
*
                  CALL ZGERC( MPC2, NQV, -TAULOC( 1 ), WORK( IPW ), 1,
     $                        WORK, 1, C( IOFFC2 ), LDC )
               END IF
*
            ELSE
*
*              Transpose and broadcast column vector V (ICOFFC2=IROFFV)
*
               IPW = NQV+1
               CALL PBZTRNV( ICTXT, 'Columnwise', 'Transpose', N,
     $                       DESCV( MB_ ), ICOFFC2, V( IOFFV ), 1, ZERO,
     $                       WORK, 1, IVROW, IVCOL, -1, ICCOL2,
     $                       WORK( IPW ) )
*
*              Perform the local computation within a process column
*
               IF( MYCOL.EQ.IVCOL ) THEN
*
                  CALL ZGEBS2D( ICTXT, 'Rowwise', ' ', 1, 1, TAU( JJV ),
     $                          1 )
                  TAULOC( 1 ) = DCONJG( TAU( JJV ) )
*
               ELSE
*
                  CALL ZGEBR2D( ICTXT, 'Rowwise', ' ', 1, 1, TAULOC, 1,
     $                          MYROW, IVCOL )
                  TAULOC( 1 ) = DCONJG( TAULOC( 1 ) )
*
               END IF
*
               IF( TAULOC( 1 ).NE.ZERO ) THEN
*
*                 w := sub( C ) * v
*
                  IF( NQV.GT.0 ) THEN
                     CALL ZGEMV( 'No transpose', MPC2, NQV, ONE,
     $                           C( IOFFC2 ), LDC, WORK, 1, ZERO,
     $                           WORK( IPW ), 1 )
                  ELSE
                     CALL ZLASET( 'All', MPC2, 1, ZERO, ZERO,
     $                            WORK( IPW ), MAX( 1, MPC2 ) )
                  END IF
                  IF( MYCOL.EQ.ICCOL1 )
     $               CALL ZAXPY( MPC2, ONE, C( IOFFC1 ), 1,
     $                           WORK( IPW ), 1 )
                  CALL ZGSUM2D( ICTXT, 'Rowwise', ' ', MPC2, 1,
     $                          WORK( IPW ), MAX( 1, MPC2 ), RDEST,
     $                          ICCOL2 )
                  IF( MYCOL.EQ.ICCOL1 )
     $               CALL ZAXPY( MPC2, -TAULOC( 1 ), WORK( IPW ), 1,
     $                           C( IOFFC1 ), 1 )
*
*                 sub( C ) := sub( C ) - w * v'
*
                  CALL ZGERC( MPC2, NQV, -TAULOC( 1 ), WORK( IPW ), 1,
     $                        WORK, 1, C( IOFFC2 ), LDC )
               END IF
*
            END IF
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PZLARZC
*
      END
