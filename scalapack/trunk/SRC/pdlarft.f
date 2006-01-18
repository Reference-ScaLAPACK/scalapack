      SUBROUTINE PDLARFT( DIRECT, STOREV, N, K, V, IV, JV, DESCV, TAU,
     $                    T, WORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            IV, JV, K, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCV( * )
      DOUBLE PRECISION   TAU( * ), T( * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLARFT forms the triangular factor T of a real block reflector H
*  of order n, which is defined as a product of k elementary reflectors.
*
*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*
*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*
*  If STOREV = 'C', the vector which defines the elementary reflector
*  H(i) is stored in the i-th column of the distributed matrix V, and
*
*     H  =  I - V * T * V'
*
*  If STOREV = 'R', the vector which defines the elementary reflector
*  H(i) is stored in the i-th row of the distributed matrix V, and
*
*     H  =  I - V' * T * V
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
*  DIRECT  (global input) CHARACTER*1
*          Specifies the order in which the elementary reflectors are
*          multiplied to form the block reflector:
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (global input) CHARACTER*1
*          Specifies how the vectors which define the elementary
*          reflectors are stored (see also Further Details):
*          = 'C': columnwise
*          = 'R': rowwise
*
*  N       (global input) INTEGER
*          The order of the block reflector H. N >= 0.
*
*  K       (global input) INTEGER
*          The order of the triangular factor T (= the number of
*          elementary reflectors). 1 <= K <= MB_V (= NB_V).
*
*  V       (input/output) DOUBLE PRECISION pointer into the local memory
*          to an array of local dimension (LOCr(IV+N-1),LOCc(JV+K-1))
*          if STOREV = 'C', and (LOCr(IV+K-1),LOCc(JV+N-1)) if
*          STOREV = 'R'. The distributed matrix V contains the
*          Householder vectors. See further details.
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
*  TAU     (local input) DOUBLE PRECISION array, dimension LOCr(IV+K-1)
*          if INCV = M_V, and LOCc(JV+K-1) otherwise. This array
*          contains the Householder scalars related to the Householder
*          vectors.  TAU is tied to the distributed matrix V.
*
*  T       (local output) DOUBLE PRECISION array, dimension (NB_V,NB_V)
*          if STOREV = 'Col', and (MB_V,MB_V) otherwise. It contains
*          the k-by-k triangular factor of the block reflector asso-
*          ciated with V. If DIRECT = 'F', T is upper triangular;
*          if DIRECT = 'B', T is lower triangular.
*
*  WORK    (local workspace) DOUBLE PRECISION array,
*                                          dimension (K*(K-1)/2)
*
*  Further Details
*  ===============
*
*  The shape of the matrix V and the storage of the vectors which define
*  the H(i) is best illustrated by the following example with n = 5 and
*  k = 3. The elements equal to 1 are not stored; the corresponding
*  array elements are modified but restored on exit. The rest of the
*  array is not used.
*
*  DIRECT = 'F' and STOREV = 'C':   DIRECT = 'F' and STOREV = 'R':
*
*  V( IV:IV+N-1,    (  1       )    V( IV:IV+K-1,    (  1 v1 v1 v1 v1 )
*     JV:JV+K-1 ) = ( v1  1    )       JV:JV+N-1 ) = (     1 v2 v2 v2 )
*                   ( v1 v2  1 )                     (        1 v3 v3 )
*                   ( v1 v2 v3 )
*                   ( v1 v2 v3 )
*
*  DIRECT = 'B' and STOREV = 'C':   DIRECT = 'B' and STOREV = 'R':
*
*  V( IV:IV+N-1,    ( v1 v2 v3 )    V( IV:IV+K-1,    ( v1 v1  1       )
*     JV:JV+K-1 ) = ( v1 v2 v3 )       JV:JV+N-1 ) = ( v2 v2 v2  1    )
*                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*                   (     1 v3 )
*                   (        1 )
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
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FORWARD
      INTEGER            ICOFF, ICTXT, II, IIV, IROFF, IVCOL, IVROW,
     $                   ITMP0, ITMP1, IW, JJ, JJV, LDV, MICOL, MIROW,
     $                   MYCOL, MYROW, NP, NPCOL, NPROW, NQ
      DOUBLE PRECISION   VII
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DCOPY, DGEMV, DGSUM2D,
     $                   DLASET, DTRMV, INFOG2L
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 .OR. K.LE.0 )
     $   RETURN
*
      ICTXT = DESCV( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      FORWARD = LSAME( DIRECT, 'F' )
      CALL INFOG2L( IV, JV, DESCV, NPROW, NPCOL, MYROW, MYCOL,
     $              IIV, JJV, IVROW, IVCOL )
*
      IF( LSAME( STOREV, 'C' ) .AND. MYCOL.EQ.IVCOL ) THEN
*
         IW = 1
         LDV = DESCV( LLD_ )
         IROFF = MOD( IV-1, DESCV( MB_ ) )
*
         IF( FORWARD ) THEN
*
*           DIRECT = 'Forward', STOREV = 'Columnwise'
*
            NP = NUMROC( N+IROFF, DESCV( MB_ ), MYROW, IVROW, NPROW )
            IF( MYROW.EQ.IVROW ) THEN
               NP = NP - IROFF
               II = IIV  + 1
            ELSE
               II = IIV
            END IF
            IF( IROFF+1.EQ.DESCV( MB_ ) ) THEN
               MIROW = MOD( IVROW+1, NPROW )
            ELSE
               MIROW = IVROW
            END IF
            ITMP0 = 0
*
            DO 10 JJ = JJV+1, JJV+K-1
*
               IF( MYROW.EQ.MIROW ) THEN
                  VII = V( II+(JJ-1)*LDV )
                  V( II+(JJ-1)*LDV ) = ONE
               END IF
*
*              T(1:i-1,i) = -tau( jv+i-1 ) *
*              V(iv+i-1:iv+n-1,jv:jv+i-2)' * V(iv+i-1:iv+n-1,jv+i-1)
*
               ITMP0 = ITMP0 + 1
               IF( NP-II+IIV.GT.0 ) THEN
                  CALL DGEMV( 'Transpose', NP-II+IIV, ITMP0,
     $                        -TAU( JJ ), V( II+(JJV-1)*LDV ), LDV,
     $                        V( II+(JJ-1)*LDV ), 1, ZERO,
     $                        WORK( IW ), 1 )
               ELSE
                  CALL DLASET( 'All', ITMP0, 1, ZERO, ZERO, WORK( IW ),
     $                         ITMP0 )
               END IF
*
               IW = IW + ITMP0
               IF( MYROW.EQ.MIROW ) THEN
                  V( II+(JJ-1)*LDV ) = VII
                  II = II + 1
               END IF
*
               IF( MOD( IV+ITMP0, DESCV( MB_ ) ).EQ.0 )
     $            MIROW = MOD( MIROW+1, NPROW )
*
   10       CONTINUE
*
            CALL DGSUM2D( ICTXT, 'Columnwise', ' ', IW-1, 1, WORK, IW-1,
     $                    IVROW, MYCOL )
*
            IF( MYROW.EQ.IVROW ) THEN
*
               IW = 1
               ITMP0 = 0
               ITMP1 = 1
*
               T( ITMP1 ) = TAU( JJV )
*
               DO 20 JJ = JJV+1, JJV+K-1
*
*                 T(1:j-1,j) = T(1:j-1,1:j-1) * T(1:j-1,j)
*
                  ITMP0 = ITMP0 + 1
                  ITMP1 = ITMP1 + DESCV( NB_ )
                  CALL DCOPY( ITMP0, WORK( IW ), 1, T( ITMP1 ), 1 )
                  IW = IW + ITMP0
*
                  CALL DTRMV( 'Upper', 'No transpose', 'Non-unit',
     $                        ITMP0, T, DESCV( NB_ ), T( ITMP1 ), 1 )
                  T(ITMP1+ITMP0) = TAU( JJ )
*
   20          CONTINUE
*
            END IF
*
         ELSE
*
*           DIRECT = 'Backward', STOREV = 'Columnwise'
*
            NP = NUMROC( N+IROFF-1, DESCV( MB_ ), MYROW, IVROW, NPROW )
            IF( MYROW.EQ.IVROW )
     $         NP = NP - IROFF
            MIROW = INDXG2P( IV+N-2, DESCV( MB_ ), MYROW,
     $                       DESCV( RSRC_ ), NPROW )
            II = IIV + NP - 1
            ITMP0 = 0
*
            DO 30 JJ = JJV+K-2, JJV, -1
*
               IF( MYROW.EQ.MIROW ) THEN
                  VII = V( II+(JJ-1)*LDV )
                  V( II+(JJ-1)*LDV ) = ONE
               END IF
*
*              T(1:i-1,i) = -tau( jv+i-1 ) *
*              V(iv:iv+n-k+i-1,jv+i:jv+k-1)' * V(iv:iv+n-k+i-1,jv+i-1)
*
               ITMP0 = ITMP0 + 1
               IF( II-IIV+1.GT.0 ) THEN
                  CALL DGEMV( 'Transpose', II-IIV+1, ITMP0, -TAU( JJ ),
     $                        V( IIV+JJ*LDV ), LDV,
     $                        V( IIV+(JJ-1)*LDV ), 1, ZERO,
     $                        WORK( IW ), 1 )
               ELSE
                  CALL DLASET( 'All', ITMP0, 1, ZERO, ZERO, WORK( IW ),
     $                         ITMP0 )
               END IF
*
               IW = IW + ITMP0
               IF( MYROW.EQ.MIROW ) THEN
                  V( II+(JJ-1)*LDV ) = VII
                  II = II - 1
               END IF
*
               IF( MOD( IV+N-ITMP0-2, DESCV(MB_) ).EQ.0 )
     $            MIROW = MOD( MIROW+NPROW-1, NPROW )
*
   30       CONTINUE
*
            CALL DGSUM2D( ICTXT, 'Columnwise', ' ', IW-1, 1, WORK, IW-1,
     $                    IVROW, MYCOL )
*
            IF( MYROW.EQ.IVROW ) THEN
*
               IW = 1
               ITMP0 = 0
               ITMP1 = K + 1 + (K-1) * DESCV( NB_ )
*
               T( ITMP1-1 ) = TAU( JJV+K-1 )
*
               DO 40 JJ = JJV+K-2, JJV, -1
*
*                 T(j+1:k,j) = T(j+1:k,j+1:k) * T(j+1:k,j)
*
                  ITMP0 = ITMP0 + 1
                  ITMP1 = ITMP1 - DESCV( NB_ ) - 1
                  CALL DCOPY( ITMP0, WORK( IW ), 1, T( ITMP1 ), 1 )
                  IW = IW + ITMP0
*
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit',
     $                        ITMP0, T( ITMP1+DESCV( NB_ ) ),
     $                        DESCV( NB_ ), T( ITMP1 ), 1 )
                  T( ITMP1-1 ) = TAU( JJ )
*
   40          CONTINUE
*
            END IF
*
         END IF
*
      ELSE IF( LSAME( STOREV, 'R' ) .AND. MYROW.EQ.IVROW ) THEN
*
         IW = 1
         LDV = DESCV( LLD_ )
         ICOFF = MOD( JV-1, DESCV( NB_ ) )
*
         IF( FORWARD ) THEN
*
*           DIRECT = 'Forward', STOREV = 'Rowwise'
*
            NQ = NUMROC( N+ICOFF, DESCV( NB_ ), MYCOL, IVCOL, NPCOL )
            IF( MYCOL.EQ.IVCOL ) THEN
               NQ = NQ - ICOFF
               JJ = JJV  + 1
            ELSE
               JJ = JJV
            END IF
            IF( ICOFF+1.EQ.DESCV( NB_ ) ) THEN
               MICOL = MOD( IVCOL+1, NPCOL )
            ELSE
               MICOL = IVCOL
            END IF
            ITMP0 = 0
*
            DO 50 II = IIV+1, IIV+K-1
*
               IF( MYCOL.EQ.MICOL ) THEN
                  VII = V( II+(JJ-1)*LDV )
                  V( II+(JJ-1)*LDV ) = ONE
               END IF
*
*              T(1:i-1,i) = -tau( iv+i-1 ) *
*              V(iv+i-1,jv+i-1:jv+n-1) * V(iv:iv+i-2,jv+i-1:jv+n-1)'
*
               ITMP0 = ITMP0 + 1
               IF( NQ-JJ+JJV.GT.0 ) THEN
                  CALL DGEMV( 'No transpose', ITMP0, NQ-JJ+JJV,
     $                        -TAU(II), V( IIV+(JJ-1)*LDV ), LDV,
     $                        V( II+(JJ-1)*LDV ), LDV, ZERO,
     $                        WORK( IW ), 1 )
               ELSE
                  CALL DLASET( 'All', ITMP0, 1, ZERO, ZERO,
     $                         WORK( IW ), ITMP0 )
               END IF
*
               IW = IW + ITMP0
               IF( MYCOL.EQ.MICOL ) THEN
                  V( II+(JJ-1)*LDV ) = VII
                  JJ = JJ + 1
               END IF
*
               IF( MOD( JV+ITMP0, DESCV( NB_ ) ).EQ.0 )
     $            MICOL = MOD( MICOL+1, NPCOL )
*
   50       CONTINUE
*
            CALL DGSUM2D( ICTXT, 'Rowwise', ' ', IW-1, 1, WORK, IW-1,
     $                    MYROW, IVCOL )
*
            IF( MYCOL.EQ.IVCOL ) THEN
*
               IW = 1
               ITMP0 = 0
               ITMP1 = 1
*
               T( ITMP1 ) = TAU( IIV )
*
               DO 60 II = IIV+1, IIV+K-1
*
*                 T(1:i-1,i) = T(1:i-1,1:i-1) * T(1:i-1,i)
*
                  ITMP0 = ITMP0 + 1
                  ITMP1 = ITMP1 + DESCV( MB_ )
                  CALL DCOPY( ITMP0, WORK( IW ), 1, T( ITMP1 ), 1 )
                  IW = IW + ITMP0
*
                  CALL DTRMV( 'Upper', 'No transpose', 'Non-unit',
     $                        ITMP0, T, DESCV( MB_ ), T( ITMP1 ), 1 )
                  T( ITMP1+ITMP0 ) = TAU( II )
*
   60          CONTINUE
*
            END IF
*
         ELSE
*
*           DIRECT = 'Backward', STOREV = 'Rowwise'
*
            NQ = NUMROC( N+ICOFF-1, DESCV( NB_ ), MYCOL, IVCOL, NPCOL )
            IF( MYCOL.EQ.IVCOL )
     $         NQ = NQ - ICOFF
            MICOL = INDXG2P( JV+N-2, DESCV( NB_ ), MYCOL,
     $                       DESCV( CSRC_ ), NPCOL )
            JJ = JJV + NQ - 1
            ITMP0 = 0
*
            DO 70 II = IIV+K-2, IIV, -1
*
               IF( MYCOL.EQ.MICOL ) THEN
                  VII = V( II+(JJ-1)*LDV )
                  V( II+(JJ-1)*LDV ) = ONE
               END IF
*
*              T(i+1:k,i) = -tau( iv+i-1 ) *
*              V(iv+i:iv+k-1,jv:jv+n-k+i-1)' * V(iv+i-1,jv:jv+n-k+i-1)'
*
               ITMP0 = ITMP0 + 1
               IF( JJ-JJV+1.GT.0 ) THEN
                  CALL DGEMV( 'No transpose', ITMP0, JJ-JJV+1,
     $                        -TAU( II ), V( II+1+(JJV-1)*LDV ), LDV,
     $                        V( II+(JJV-1)*LDV ), LDV, ZERO,
     $                        WORK( IW ), 1 )
               ELSE
                  CALL DLASET( 'All', ITMP0, 1, ZERO, ZERO,
     $                         WORK( IW ), ITMP0 )
               END IF
*
               IW = IW + ITMP0
               IF( MYCOL.EQ.MICOL ) THEN
                  V( II+(JJ-1)*LDV ) = VII
                  JJ = JJ - 1
               END IF
*
               IF( MOD( JV+N-ITMP0-2, DESCV( NB_ ) ).EQ.0 )
     $            MICOL = MOD( MICOL+NPCOL-1, NPCOL )
*
   70       CONTINUE
*
            CALL DGSUM2D( ICTXT, 'Rowwise', ' ', IW-1, 1, WORK, IW-1,
     $                    MYROW, IVCOL )
*
            IF( MYCOL.EQ.IVCOL ) THEN
*
               IW = 1
               ITMP0 = 0
               ITMP1 = K + 1 + (K-1) * DESCV( MB_ )
*
               T( ITMP1-1 ) = TAU( IIV+K-1 )
*
               DO 80 II = IIV+K-2, IIV, -1
*
*                 T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)
*
                  ITMP0 = ITMP0 + 1
                  ITMP1 = ITMP1 - DESCV( MB_ ) - 1
                  CALL DCOPY( ITMP0, WORK( IW ), 1, T( ITMP1 ), 1 )
                  IW = IW + ITMP0
*
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit',
     $                        ITMP0, T( ITMP1+DESCV( MB_ ) ),
     $                        DESCV( MB_ ), T( ITMP1 ), 1 )
                  T( ITMP1-1 ) = TAU( II )
*
   80          CONTINUE
*
            END IF
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PDLARFT
*
      END
