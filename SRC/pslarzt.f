      SUBROUTINE PSLARZT( DIRECT, STOREV, N, K, V, IV, JV, DESCV, TAU,
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
      REAL               TAU( * ), T( * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLARZT forms the triangular factor T of a real block reflector
*  H of order > n, which is defined as a product of k elementary
*  reflectors as returned by PSTZRZF.
*
*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*
*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*
*  If STOREV = 'C', the vector which defines the elementary reflector
*  H(i) is stored in the i-th column of the array V, and
*
*     H  =  I - V * T * V'
*
*  If STOREV = 'R', the vector which defines the elementary reflector
*  H(i) is stored in the i-th row of the array V, and
*
*     H  =  I - V' * T * V
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
*  DIRECT  (global input) CHARACTER
*          Specifies the order in which the elementary reflectors are
*          multiplied to form the block reflector:
*          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (global input) CHARACTER
*          Specifies how the vectors which define the elementary
*          reflectors are stored (see also Further Details):
*          = 'C': columnwise                        (not supported yet)
*          = 'R': rowwise
*
*  N       (global input) INTEGER
*          The number of meaningful entries of the block reflector H.
*          N >= 0.
*
*  K       (global input) INTEGER
*          The order of the triangular factor T (= the number of
*          elementary reflectors). 1 <= K <= MB_V (= NB_V).
*
*  V       (input/output) REAL pointer into the local memory
*          to an array of local dimension (LOCr(IV+K-1),LOCc(JV+N-1)).
*          The distributed matrix V contains the Householder vectors.
*          See further details.
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
*  TAU     (local input) REAL, array, dimension LOCr(IV+K-1)
*          if INCV = M_V, and LOCc(JV+K-1) otherwise. This array
*          contains the Householder scalars related to the Householder
*          vectors.  TAU is tied to the distributed matrix V.
*
*  T       (local output) REAL array, dimension (MB_V,MB_V)
*          It contains the k-by-k triangular factor of the block
*          reflector associated with V. T is lower triangular.
*
*  WORK    (local workspace) REAL array,
*                                           dimension (K*(K-1)/2)
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
*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*
*                                              ______V_____
*         ( v1 v2 v3 )                        /            \
*         ( v1 v2 v3 )                      ( v1 v1 v1 v1 v1 . . . . 1 )
*     V = ( v1 v2 v3 )                      ( v2 v2 v2 v2 v2 . . . 1   )
*         ( v1 v2 v3 )                      ( v3 v3 v3 v3 v3 . . 1     )
*         ( v1 v2 v3 )
*            .  .  .
*            .  .  .
*            1  .  .
*               1  .
*                  1
*
*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*
*                                                        ______V_____
*            1                                          /            \
*            .  1                           ( 1 . . . . v1 v1 v1 v1 v1 )
*            .  .  1                        ( . 1 . . . v2 v2 v2 v2 v2 )
*            .  .  .                        ( . . 1 . . v3 v3 v3 v3 v3 )
*            .  .  .
*         ( v1 v2 v3 )
*         ( v1 v2 v3 )
*     V = ( v1 v2 v3 )
*         ( v1 v2 v3 )
*         ( v1 v2 v3 )
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICOFF, ICTXT, II, IIV, INFO, IVCOL, IVROW,
     $                   ITMP0, ITMP1, IW, JJV, LDV, MYCOL, MYROW,
     $                   NPCOL, NPROW, NQ
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GRIDINFO, INFOG2L, PXERBLA,
     $                   SCOPY, SGEMV, SGSUM2D, SLASET,
     $                   STRMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      EXTERNAL           LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCV( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Check for currently supported options
*
      INFO = 0
      IF( .NOT.LSAME( DIRECT, 'B' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( STOREV, 'R' ) ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PSLARZT', -INFO )
         CALL BLACS_ABORT( ICTXT, 1 )
         RETURN
      END IF
*
      CALL INFOG2L( IV, JV, DESCV, NPROW, NPCOL, MYROW, MYCOL,
     $              IIV, JJV, IVROW, IVCOL )
*
      IF( MYROW.EQ.IVROW ) THEN
         IW = 1
         ITMP0 = 0
         LDV = DESCV( LLD_ )
         ICOFF = MOD( JV-1, DESCV( NB_ ) )
         NQ = NUMROC( N+ICOFF, DESCV( NB_ ), MYCOL, IVCOL, NPCOL )
         IF( MYCOL.EQ.IVCOL )
     $      NQ = NQ - ICOFF
*
         DO 10 II = IIV+K-2, IIV, -1
*
*           T(i+1:k,i) = -tau( iv+i-1 ) *
*                     V(iv+i:iv+k-1,jv:jv+n-1) * V(iv+i-1,jv:jv+n-1)'
*
            ITMP0 = ITMP0 + 1
            IF( NQ.GT.0 ) THEN
               CALL SGEMV( 'No transpose', ITMP0, NQ, -TAU( II ),
     $                     V( II+1+(JJV-1)*LDV ), LDV,
     $                     V( II+(JJV-1)*LDV ), LDV, ZERO, WORK( IW ),
     $                     1 )
            ELSE
               CALL SLASET( 'All', ITMP0, 1, ZERO, ZERO, WORK( IW ),
     $                      ITMP0 )
            END IF
            IW = IW + ITMP0
*
   10    CONTINUE
*
         CALL SGSUM2D( ICTXT, 'Rowwise', ' ', IW-1, 1, WORK, IW-1,
     $                 MYROW, IVCOL )
*
         IF( MYCOL.EQ.IVCOL ) THEN
*
            IW = 1
            ITMP0 = 0
            ITMP1 = K + 1 + (K-1) * DESCV( MB_ )
*
            T( ITMP1-1 ) = TAU( IIV+K-1 )
*
            DO 20 II = IIV+K-2, IIV, -1
*
*              T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)
*
               ITMP0 = ITMP0 + 1
               ITMP1 = ITMP1 - DESCV( MB_ ) - 1
               CALL SCOPY( ITMP0, WORK( IW ), 1, T( ITMP1 ), 1 )
               IW = IW + ITMP0
*
               CALL STRMV( 'Lower', 'No transpose', 'Non-unit', ITMP0,
     $                     T( ITMP1+DESCV( MB_ ) ), DESCV( MB_ ),
     $                     T( ITMP1 ), 1 )
               T( ITMP1-1 ) = TAU( II )
*
   20       CONTINUE
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PSLARZT
*
      END
