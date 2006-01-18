      SUBROUTINE PDLARF( SIDE, M, N, V, IV, JV, DESCV, INCV, TAU,
     $                   C, IC, JC, DESCC, WORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            IC, INCV, IV, JC, JV, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCC( * ), DESCV( * )
      DOUBLE PRECISION   C( * ), TAU( * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLARF applies a real elementary reflector Q (or Q**T) to a real
*  M-by-N distributed matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1), from
*  either the left or the right. Q is represented in the form
*
*        Q = I - tau * v * v'
*
*  where tau is a real scalar and v is a real vector.
*
*  If tau = 0, then Q is taken to be the unit matrix.
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
*  entry V(IV,JV) must also have the first row of sub( C ). Moreover,
*  MOD(IV-1,MB_V) must be equal to MOD(IC-1,MB_C), if INCV=M_V, only
*  the last equality must be satisfied.
*
*  If SIDE = 'Right' and INCV = M_V then the column process having the
*  first entry V(IV,JV) must also have the first column of sub( C ) and
*  MOD(JV-1,NB_V) must be equal to MOD(JC-1,NB_C), if INCV = 1 only the
*  last equality must be satisfied.
*
*  Arguments
*  =========
*
*  SIDE    (global input) CHARACTER
*          = 'L': form  Q * sub( C ),
*          = 'R': form  sub( C ) * Q, Q = Q**T.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( C ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( C ). N >= 0.
*
*  V       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_V,*) containing the local
*          pieces of the distributed vectors V representing the
*          Householder transformation Q,
*             V(IV:IV+M-1,JV) if SIDE = 'L' and INCV = 1,
*             V(IV,JV:JV+M-1) if SIDE = 'L' and INCV = M_V,
*             V(IV:IV+N-1,JV) if SIDE = 'R' and INCV = 1,
*             V(IV,JV:JV+N-1) if SIDE = 'R' and INCV = M_V,
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
*  TAU     (local input) DOUBLE PRECISION array, dimension  LOCc(JV) if
*          INCV = 1, and LOCr(IV) otherwise. This array contains the
*          Householder scalars related to the Householder vectors.
*          TAU is tied to the distributed matrix V.
*
*  C       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_C, LOCc(JC+N-1) ),
*          containing the local pieces of sub( C ). On exit, sub( C )
*          is overwritten by the Q * sub( C ) if SIDE = 'L', or
*          sub( C ) * Q if SIDE = 'R'.
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
*  WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK)
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
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE  = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            CCBLCK, CRBLCK
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            ICCOL, ICOFF, ICROW, ICTXT, IIC, IIV, IOFFC,
     $                   IOFFV, IPW, IROFF, IVCOL, IVROW, JJC, JJV, LDC,
     $                   LDV, MYCOL, MYROW, MP, NCC, NCV, NPCOL, NPROW,
     $                   NQ, RDEST
      DOUBLE PRECISION   TAULOC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DCOPY, DGEBR2D, DGEBS2D,
     $                   DGEMV, DGER, DGERV2D, DGESD2D,
     $                   DGSUM2D, DLASET, INFOG2L, PB_TOPGET,
     $                   PBDTRNV
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
      CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL, IIC, JJC,
     $              ICROW, ICCOL )
      CALL INFOG2L( IV, JV, DESCV, NPROW, NPCOL, MYROW, MYCOL, IIV, JJV,
     $              IVROW, IVCOL )
      NCC = NUMROC( DESCC( N_ ), DESCC( NB_ ), MYCOL, DESCC( CSRC_ ),
     $              NPCOL )
      NCV = NUMROC( DESCV( N_ ), DESCV( NB_ ), MYCOL, DESCV( CSRC_ ),
     $              NPCOL )
      LDC = DESCC( LLD_ )
      LDV = DESCV( LLD_ )
      IIC = MIN( IIC, LDC )
      IIV = MIN( IIV, LDV )
      JJC = MIN( JJC, NCC )
      JJV = MIN( JJV, NCV )
      IOFFC = IIC+(JJC-1)*LDC
      IOFFV = IIV+(JJV-1)*LDV
*
      IROFF = MOD( IC-1, DESCC( MB_ ) )
      ICOFF = MOD( JC-1, DESCC( NB_ ) )
      MP = NUMROC( M+IROFF, DESCC( MB_ ), MYROW, ICROW, NPROW )
      NQ = NUMROC( N+ICOFF, DESCC( NB_ ), MYCOL, ICCOL, NPCOL )
      IF( MYROW.EQ.ICROW )
     $   MP = MP - IROFF
      IF( MYCOL.EQ.ICCOL )
     $   NQ = NQ - ICOFF
*
*     Is sub( C ) only distributed over a process row ?
*
      CRBLCK = ( M.LE.(DESCC( MB_ )-IROFF) )
*
*     Is sub( C ) only distributed over a process column ?
*
      CCBLCK = ( N.LE.(DESCC( NB_ )-ICOFF) )
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
         IF( CRBLCK ) THEN
            RDEST = ICROW
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
*              Transpose row vector V
*
               IPW = MP+1
               CALL PBDTRNV( ICTXT, 'Rowwise', 'Transpose', M,
     $                       DESCV( NB_ ), IROFF, V( IOFFV ), LDV, ZERO,
     $                       WORK, 1, IVROW, IVCOL, ICROW, ICCOL,
     $                       WORK( IPW ) )
*
*              Perform the local computation within a process column
*
               IF( MYCOL.EQ.ICCOL ) THEN
*
                  IF( MYROW.EQ.IVROW ) THEN
*
                     CALL DGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1,
     $                             TAU( IIV ), 1 )
                     TAULOC = TAU( IIV )
*
                  ELSE
*
                     CALL DGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1,
     $                             TAULOC, 1, IVROW, MYCOL )
*
                  END IF
*
                  IF( TAULOC.NE.ZERO ) THEN
*
*                    w := sub( C )' * v
*
                     IF( MP.GT.0 ) THEN
                        CALL DGEMV( 'Transpose', MP, NQ, ONE,
     $                              C( IOFFC ), LDC, WORK, 1, ZERO,
     $                              WORK( IPW ), 1 )
                     ELSE
                        CALL DLASET( 'All', NQ, 1, ZERO, ZERO,
     $                               WORK( IPW ), MAX( 1, NQ ) )
                     END IF
                     CALL DGSUM2D( ICTXT, 'Columnwise', ' ', NQ, 1,
     $                             WORK( IPW ), MAX( 1, NQ ), RDEST,
     $                             MYCOL )
*
*                    sub( C ) := sub( C ) - v * w'
*
                     CALL DGER( MP, NQ, -TAULOC, WORK, 1, WORK( IPW ),
     $                          1, C( IOFFC ), LDC )
                  END IF
*
               END IF
*
            ELSE
*
*              V is a column vector
*
               IF( IVCOL.EQ.ICCOL ) THEN
*
*                 Perform the local computation within a process column
*
                  IF( MYCOL.EQ.ICCOL ) THEN
*
                     TAULOC = TAU( JJV )
*
                     IF( TAULOC.NE.ZERO ) THEN
*
*                       w := sub( C )' * v
*
                        IF( MP.GT.0 ) THEN
                           CALL DGEMV( 'Transpose', MP, NQ, ONE,
     $                                 C( IOFFC ), LDC, V( IOFFV ), 1,
     $                                 ZERO, WORK, 1 )
                        ELSE
                           CALL DLASET( 'All', NQ, 1, ZERO, ZERO,
     $                                  WORK, MAX( 1, NQ ) )
                        END IF
                        CALL DGSUM2D( ICTXT, 'Columnwise', ' ', NQ, 1,
     $                                WORK, MAX( 1, NQ ), RDEST, MYCOL )
*
*                       sub( C ) := sub( C ) - v * w'
*
                        CALL DGER( MP, NQ, -TAULOC, V( IOFFV ), 1, WORK,
     $                             1, C( IOFFC ), LDC )
                     END IF
*
                  END IF
*
               ELSE
*
*                 Send V and TAU to the process column ICCOL
*
                  IF( MYCOL.EQ.IVCOL ) THEN
*
                     IPW = MP+1
                     CALL DCOPY( MP, V( IOFFV ), 1, WORK, 1 )
                     WORK( IPW ) = TAU( JJV )
                     CALL DGESD2D( ICTXT, IPW, 1, WORK, IPW, MYROW,
     $                             ICCOL )
*
                  ELSE IF( MYCOL.EQ.ICCOL ) THEN
*
                     IPW = MP+1
                     CALL DGERV2D( ICTXT, IPW, 1, WORK, IPW, MYROW,
     $                             IVCOL )
                     TAULOC = WORK( IPW )
*
                     IF( TAULOC.NE.ZERO ) THEN
*
*                       w := sub( C )' * v
*
                        IF( MP.GT.0 ) THEN
                           CALL DGEMV( 'Transpose', MP, NQ, ONE,
     $                                 C( IOFFC ), LDC, WORK, 1, ZERO,
     $                                 WORK( IPW ), 1 )
                        ELSE
                           CALL DLASET( 'All', NQ, 1, ZERO, ZERO,
     $                                  WORK( IPW ), MAX( 1, NQ ) )
                        END IF
                        CALL DGSUM2D( ICTXT, 'Columnwise', ' ', NQ, 1,
     $                                WORK( IPW ), MAX( 1, NQ ), RDEST,
     $                                MYCOL )
*
*                       sub( C ) := sub( C ) - v * w'
*
                        CALL DGER( MP, NQ, -TAULOC, WORK, 1,
     $                             WORK( IPW ), 1, C( IOFFC ), LDC )
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
*              Transpose and broadcast row vector V
*
               IPW = MP+1
               CALL PBDTRNV( ICTXT, 'Rowwise', 'Transpose', M,
     $                       DESCV( NB_ ), IROFF, V( IOFFV ), LDV, ZERO,
     $                       WORK, 1, IVROW, IVCOL, ICROW, -1,
     $                       WORK( IPW ) )
*
*              Perform the local computation within a process column
*
               IF( MYROW.EQ.IVROW ) THEN
*
                  CALL DGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1,
     $                          TAU( IIV ), 1 )
                  TAULOC = TAU( IIV )
*
               ELSE
*
                  CALL DGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1, TAULOC,
     $                          1, IVROW, MYCOL )
*
               END IF
*
               IF( TAULOC.NE.ZERO ) THEN
*
*                 w := sub( C )' * v
*
                  IF( MP.GT.0 ) THEN
                     IF( IOFFC.GT.0 )
     $                  CALL DGEMV( 'Transpose', MP, NQ, ONE,
     $                              C( IOFFC ), LDC, WORK, 1, ZERO,
     $                              WORK( IPW ), 1 )
                  ELSE
                     CALL DLASET( 'All', NQ, 1, ZERO, ZERO,
     $                            WORK( IPW ), MAX( 1, NQ ) )
                  END IF
                  CALL DGSUM2D( ICTXT, 'Columnwise', ' ', NQ, 1,
     $                          WORK( IPW ), MAX( 1, NQ ), RDEST,
     $                          MYCOL )
*
*                 sub( C ) := sub( C ) - v * w'
*
                  IF( IOFFC.GT.0 )
     $               CALL DGER( MP, NQ, -TAULOC, WORK, 1, WORK( IPW ),
     $                          1, C( IOFFC ), LDC )
               END IF
*
            ELSE
*
*              Broadcast column vector V
*
               CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
               IF( MYCOL.EQ.IVCOL ) THEN
*
                  IPW = MP+1
                  CALL DCOPY( MP, V( IOFFV ), 1, WORK, 1 )
                  WORK(IPW) = TAU( JJV )
                  CALL DGEBS2D( ICTXT, 'Rowwise', ROWBTOP, IPW, 1,
     $                          WORK, IPW )
                  TAULOC = TAU( JJV )
*
               ELSE
*
                  IPW = MP+1
                  CALL DGEBR2D( ICTXT, 'Rowwise', ROWBTOP, IPW, 1, WORK,
     $                          IPW, MYROW, IVCOL )
                  TAULOC = WORK( IPW )
*
               END IF
*
               IF( TAULOC.NE.ZERO ) THEN
*
*                 w := sub( C )' * v
*
                  IF( MP.GT.0 ) THEN
                     IF( IOFFC.GT.0 )
     $                  CALL DGEMV( 'Transpose', MP, NQ, ONE,
     $                              C( IOFFC ), LDC, WORK, 1, ZERO,
     $                              WORK( IPW ), 1 )
                  ELSE
                     CALL DLASET( 'All', NQ, 1, ZERO, ZERO,
     $                            WORK( IPW ), MAX( 1, NQ ) )
                  END IF
                  CALL DGSUM2D( ICTXT, 'Columnwise', ' ', NQ, 1,
     $                          WORK( IPW ), MAX( 1, NQ ), RDEST,
     $                          MYCOL )
*
*                 sub( C ) := sub( C ) - v * w'
*
                  IF( IOFFC.GT.0 )
     $               CALL DGER( MP, NQ, -TAULOC, WORK, 1, WORK( IPW ),
     $                          1, C( IOFFC ), LDC )
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
               IF( IVROW.EQ.ICROW ) THEN
*
*                 Perform the local computation within a process row
*
                  IF( MYROW.EQ.ICROW ) THEN
*
                     TAULOC = TAU( IIV )
*
                     IF( TAULOC.NE.ZERO ) THEN
*
*                       w := sub( C ) * v
*
                        IF( NQ.GT.0 ) THEN
                           CALL DGEMV( 'No transpose', MP, NQ, ONE,
     $                                 C( IOFFC ), LDC, V( IOFFV ), LDV,
     $                                 ZERO, WORK, 1 )
                        ELSE
                           CALL DLASET( 'All', MP, 1, ZERO, ZERO,
     $                                  WORK, MAX( 1, MP ) )
                        END IF
                        CALL DGSUM2D( ICTXT, 'Rowwise', ' ', MP, 1,
     $                                WORK, MAX( 1, MP ), RDEST, ICCOL )
*
*                       sub( C ) := sub( C ) - w * v'
*
                        IF( IOFFV.GT.0 .AND. IOFFC.GT.0 )
     $                     CALL DGER( MP, NQ, -TAULOC, WORK, 1,
     $                                V( IOFFV ), LDV, C( IOFFC ), LDC )
                     END IF
*
                  END IF
*
               ELSE
*
*                 Send V and TAU to the process row ICROW
*
                  IF( MYROW.EQ.IVROW ) THEN
*
                     IPW = NQ+1
                     CALL DCOPY( NQ, V( IOFFV ), LDV, WORK, 1 )
                     WORK(IPW) = TAU( IIV )
                     CALL DGESD2D( ICTXT, IPW, 1, WORK, IPW, ICROW,
     $                             MYCOL )
*
                  ELSE IF( MYROW.EQ.ICROW ) THEN
*
                     IPW = NQ+1
                     CALL DGERV2D( ICTXT, IPW, 1, WORK, IPW, IVROW,
     $                             MYCOL )
                     TAULOC = WORK( IPW )
*
                     IF( TAULOC.NE.ZERO ) THEN
*
*                       w := sub( C ) * v
*
                        IF( NQ.GT.0 ) THEN
                           CALL DGEMV( 'No transpose', MP, NQ, ONE,
     $                                 C( IOFFC ), LDC, WORK, 1, ZERO,
     $                                 WORK( IPW ), 1 )
                        ELSE
                           CALL DLASET( 'All', MP, 1, ZERO, ZERO,
     $                                  WORK( IPW ), MAX( 1, MP ) )
                        END IF
                        CALL DGSUM2D( ICTXT, 'Rowwise', ' ', MP, 1,
     $                                WORK( IPW ), MAX( 1, MP ), RDEST,
     $                                ICCOL )
*
*                       sub( C ) := sub( C ) - w * v'
*
                        CALL DGER( MP, NQ, -TAULOC, WORK( IPW ), 1,
     $                             WORK, 1, C( IOFFC ), LDC )
                     END IF
*
                  END IF
*
               END IF
*
            ELSE
*
*              Transpose column vector V
*
               IPW = NQ+1
               CALL PBDTRNV( ICTXT, 'Columnwise', 'Transpose', N,
     $                       DESCV( MB_ ), ICOFF, V( IOFFV ), 1, ZERO,
     $                       WORK, 1, IVROW, IVCOL, ICROW, ICCOL,
     $                       WORK( IPW ) )
*
*              Perform the local computation within a process column
*
               IF( MYROW.EQ.ICROW ) THEN
*
                  IF( MYCOL.EQ.IVCOL ) THEN
*
                     CALL DGEBS2D( ICTXT, 'Rowwise', ' ', 1, 1,
     $                             TAU( JJV ), 1 )
                     TAULOC = TAU( JJV )
*
                  ELSE
*
                     CALL DGEBR2D( ICTXT, 'Rowwise', ' ', 1, 1, TAULOC,
     $                             1, MYROW, IVCOL )
*
                  END IF
*
                  IF( TAULOC.NE.ZERO ) THEN
*
*                    w := sub( C ) * v
*
                     IF( NQ.GT.0 ) THEN
                        CALL DGEMV( 'No transpose', MP, NQ, ONE,
     $                              C( IOFFC ), LDC, WORK, 1, ZERO,
     $                              WORK( IPW ), 1 )
                     ELSE
                        CALL DLASET( 'All', MP, 1, ZERO, ZERO,
     $                               WORK( IPW ), MAX( 1, MP ) )
                     END IF
                     CALL DGSUM2D( ICTXT, 'Rowwise', ' ', MP, 1,
     $                             WORK( IPW ), MAX( 1, MP ), RDEST,
     $                             ICCOL )
*
*                    sub( C ) := sub( C ) - w * v'
*
                     CALL DGER( MP, NQ, -TAULOC, WORK( IPW ), 1, WORK,
     $                          1, C( IOFFC ), LDC )
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
                  IPW = NQ+1
                  IF( IOFFV.GT.0 )
     $               CALL DCOPY( NQ, V( IOFFV ), LDV, WORK, 1 )
                  WORK(IPW) = TAU( IIV )
                  CALL DGEBS2D( ICTXT, 'Columnwise', COLBTOP, IPW, 1,
     $                          WORK, IPW )
                  TAULOC = TAU( IIV )
*
               ELSE
*
                  IPW = NQ+1
                  CALL DGEBR2D( ICTXT, 'Columnwise', COLBTOP, IPW, 1,
     $                          WORK, IPW, IVROW, MYCOL )
                  TAULOC = WORK( IPW )
*
               END IF
*
               IF( TAULOC.NE.ZERO ) THEN
*
*                 w := sub( C ) * v
*
                  IF( NQ.GT.0 ) THEN
                     CALL DGEMV( 'No Transpose', MP, NQ, ONE,
     $                           C( IOFFC ), LDC, WORK, 1, ZERO,
     $                           WORK( IPW ), 1 )
                  ELSE
                     CALL DLASET( 'All', MP, 1, ZERO, ZERO,
     $                            WORK( IPW ), MAX( 1, MP ) )
                  END IF
                  CALL DGSUM2D( ICTXT, 'Rowwise', ' ', MP, 1,
     $                          WORK( IPW ), MAX( 1, MP ), RDEST,
     $                          ICCOL )
*
*                 sub( C ) := sub( C ) - w * v'
*
                  IF( IOFFC.GT.0 )
     $               CALL DGER( MP, NQ, -TAULOC, WORK( IPW ), 1, WORK,
     $                          1, C( IOFFC ), LDC )
               END IF
*
            ELSE
*
*              Transpose and broadcast column vector V
*
               IPW = NQ+1
               CALL PBDTRNV( ICTXT, 'Columnwise', 'Transpose', N,
     $                       DESCV( MB_ ), ICOFF, V( IOFFV ), 1, ZERO,
     $                       WORK, 1, IVROW, IVCOL, -1, ICCOL,
     $                       WORK( IPW ) )
*
*              Perform the local computation within a process column
*
               IF( MYCOL.EQ.IVCOL ) THEN
*
                  CALL DGEBS2D( ICTXT, 'Rowwise', ' ', 1, 1, TAU( JJV ),
     $                          1 )
                  TAULOC = TAU( JJV )
*
               ELSE
*
                  CALL DGEBR2D( ICTXT, 'Rowwise', ' ', 1, 1, TAULOC, 1,
     $                          MYROW, IVCOL )
*
               END IF
*
               IF( TAULOC.NE.ZERO ) THEN
*
*                 w := sub( C ) * v
*
                  IF( NQ.GT.0 ) THEN
                     CALL DGEMV( 'No transpose', MP, NQ, ONE,
     $                           C( IOFFC ), LDC, WORK, 1, ZERO,
     $                           WORK( IPW ), 1 )
                  ELSE
                     CALL DLASET( 'All', MP, 1, ZERO, ZERO, WORK( IPW ),
     $                            MAX( 1, MP ) )
                  END IF
                  CALL DGSUM2D( ICTXT, 'Rowwise', ' ', MP, 1,
     $                          WORK( IPW ), MAX( 1, MP ), RDEST,
     $                          ICCOL )
*
*                 sub( C ) := sub( C ) - w * v'
*
                  CALL DGER( MP, NQ, -TAULOC, WORK( IPW ), 1, WORK, 1,
     $                       C( IOFFC ), LDC )
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
*     End of PDLARF
*
      END
