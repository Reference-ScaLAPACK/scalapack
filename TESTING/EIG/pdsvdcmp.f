      SUBROUTINE PDSVDCMP( M, N, JOBTYPE, S, SC, U, UC, IU, JU, DESCU,
     $                     VT, VTC, IVT, JVT, DESCVT, THRESH, RESULT,
     $                     DELTA, WORK, LWORK )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IU, IVT, JOBTYPE, JU, JVT, LWORK, M, N
      DOUBLE PRECISION   DELTA, THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            DESCU( * ), DESCVT( * ), RESULT( * )
      DOUBLE PRECISION   S( * ), SC( * ), U( * ), UC( * ), VT( * ),
     $                   VTC( * ), WORK( * )
*     ..
*
*  Purpose
*  ========
*  Testing how accurately "full" and "partial" decomposition options
*  provided by PDGESVD correspond to each other.
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
*   Arguments
*   ==========
*
*   M        (global input) INTEGER
*            Number of rows of the distributed matrix, for which
*            SVD was calculated
*
*   N        (global input) INTEGER
*            Number of columns of the distributed matrix, for which
*            SVD was calculated
*
*   JOBTYPE  (global input) INTEGER
*            Depending on the value of this parameter,
*            the following comparisons are performed:
*
*            JOBTYPE  |   COMPARISON
*            -------------------------------------------
*               2     |  | U - UC | / ( M ulp ) > THRESH,
*               3     |  | VT - VTC | / ( N ulp ) > THRESH
*
*            In addition, for JOBTYPE = 2:4 comparison
*            | S1 - S2 | / ( SIZE ulp |S| ) > THRESH
*            is performed. Positive result of any of the comparisons
*            typically indicates erroneous computations and sets
*            to one corresponding element of array RESULT
*
*    S       (global input) DOUBLE PRECISION  array of singular values
*            calculated for JOBTYPE equal to  1
*
*    SC      (global input) DOUBLE PRECISION  array of singular values
*            calculated for JOBTYPE nonequal to 1
*
*    U       (local input) DOUBLE PRECISION array of left singular
*            vectors calculated for JOBTYPE equal to 1, local
*            dimension (MP, SIZEQ), global dimension (M, SIZE)
*
*    UC      (local input) DOUBLE PRECISION  array of left singular
*            vectors calculated for JOBTYPE non equal to 1, local
*            dimension (MP, SIZEQ), global dimension (M, SIZE)
*
*    IU      (global input) INTEGER
*            The row index in the global array U indicating the first
*            row of sub( U ).
*
*    JU      (global input) INTEGER
*            The column index in the global array U indicating the
*            first column of sub( U ).
*
*    DESCU   (global input) INTEGER array of dimension DLEN_
*            The array descriptor for the distributed matrix U and UC
*
*    V       (local input) DOUBLE PRECISION array of right singular
*            vectors calculated for JOBTYPE equal to 1, local
*            dimension (SIZEP, NQ), global dimension (SIZE, N)
*
*    VC      (local input) DOUBLE PRECISION array of right singular
*            vectors calculated for JOBTYPE non equal to 1, local
*            dimension (SIZEP, NQ), global dimension (SIZE, N)
*
*    IVT     (global input) INTEGER
*            The row index in the global array VT indicating the first
*            row of sub( VT ).
*
*    JVT     (global input) INTEGER
*            The column index in the global array VT indicating the
*            first column of sub( VT ).
*
*    DESCVT  (global input) INTEGER array of dimension DLEN_
*            The array descriptor for the distributed matrix VT and
*            VTC
*
*    THRESH  (global input) DOUBLE PRECISION
*            The threshold value for the test ratios.  A result is
*            included in the output file if RESULT >= THRESH.  The test
*            ratios are scaled to be O(1), so THRESH should be a small
*            multiple of 1, e.g., 10 or 100. To have every test ratio
*            printed, use THRESH = 0.
*
*    RESULT  (global input/output) INTEGER array.
*            Every nonzero entry corresponds to erroneous computation.
*
*    DELTA   (global output) DOUBLE PRECISION
*            maximum of the available of the following three values
*            | U - UC | / ( M ulp THRESH ),
*            | VT - VT | / ( N ulp THRESH ),
*            | S1 - S2 | / ( SIZE ulp |S| THRESH )
*
*    WORK    (local workspace/output) DOUBLE PRECISION array,
*            dimension (LWORK)
*            On exit, WORK(1) returns the optimal LWORK.
*
*    LWORK   (local input) INTEGER
*            The dimension of the array WORK.
*
* ======================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            COLPTR, I, INFO, J, LWMIN, MYCOL, MYROW, NPCOL,
     $                   NPROW, NQ, RESULTS, SIZE, SIZEPOS, SIZEQ
      DOUBLE PRECISION   ACCUR, CMP, NORMDIFS, NORMDIFU, NORMDIFV,
     $                   NORMS, ULP
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      DOUBLE PRECISION   DLANGE, PDLAMCH, PDLANGE
      EXTERNAL           NUMROC, DLANGE, PDLAMCH, PDLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*     This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*DLEN_*DTYPE_*MB_*M_*N_*RSRC_.LT.0 ) 
     $   RETURN
*
      RESULTS = 0
      NORMDIFS = 0
      NORMDIFU = 0
      NORMDIFV = 0
      SIZE = MIN( M, N )
*
*     Sizepos is a number of parameters to pdsvdcmp plus one. It's used
*     for the error reporting.
*
      SIZEPOS = 17
      INFO = 0
      CALL BLACS_GRIDINFO( DESCU( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -607
      ELSE
         CALL CHK1MAT( M, 1, SIZE, SIZEPOS, 1, 1, DESCU, 8, INFO )
         CALL CHK1MAT( SIZE, SIZEPOS, N, 2, 1, 1, DESCVT, 11, INFO )
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Calculate workspace.
*
         SIZEQ = NUMROC( SIZE, DESCU( NB_ ), MYCOL, 0, NPCOL )
         NQ = NUMROC( N, DESCVT( NB_ ), MYCOL, 0, NPCOL )
         LWMIN = MAX( SIZEQ, NQ ) + 4
         WORK( 1 ) = LWMIN
         IF( LWORK.EQ.-1 )
     $      GO TO 60
         IF( LWORK.LT.LWMIN ) THEN
            INFO = -16
         ELSE IF( THRESH.LE.0 ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCU( CTXT_ ), 'PDSVDCMP', -INFO )
         RETURN
      END IF
*
      ULP = PDLAMCH( DESCU( CTXT_ ), 'P' )
*
*     Make comparison of singular values.
*
      NORMS = DLANGE( '1', SIZE, 1, S, SIZE, WORK )
      DO 10 I = 1, SIZE
         SC( I ) = S( I ) - SC( I )
   10 CONTINUE
*
      NORMDIFS = DLANGE( '1', SIZE, 1, SC, SIZE, WORK )
      ACCUR = ULP*SIZE*NORMS*THRESH
*
      IF( NORMDIFS.GT.ACCUR )
     $   RESULTS = 1
      IF( NORMDIFS.EQ.0 .AND. ACCUR.EQ.0 ) THEN
         NORMDIFS = 0
      ELSE
         NORMDIFS = NORMDIFS / ACCUR
      END IF
*
      IF( JOBTYPE.EQ.2 ) THEN
*
         RESULT( 5 ) = RESULTS
         ACCUR = ULP*M*THRESH
         DO 30 J = 1, SIZEQ
            COLPTR = DESCU( LLD_ )*( J-1 )
            DO 20 I = 1, DESCU( LLD_ )
               UC( I+COLPTR ) = U( I+COLPTR ) - UC( I+COLPTR )
   20       CONTINUE
   30    CONTINUE
*
         NORMDIFU = PDLANGE( '1', M, SIZE, UC, IU, JU, DESCU, WORK )
*
         IF( NORMDIFU.GE.ACCUR )
     $      RESULT( 6 ) = 1
         IF( NORMDIFU.EQ.0 .AND. ACCUR.EQ.0 ) THEN
            NORMDIFU = 0
         ELSE
            NORMDIFU = NORMDIFU / ACCUR
         END IF
*
      ELSE IF( JOBTYPE.EQ.3 ) THEN
*
         RESULT( 7 ) = RESULTS
         ACCUR = ULP*N*THRESH
         DO 50 J = 1, NQ
            COLPTR = DESCVT( LLD_ )*( J-1 )
            DO 40 I = 1, DESCVT( LLD_ )
               VTC( I+COLPTR ) = VT( I+COLPTR ) - VTC( I+COLPTR )
   40       CONTINUE
   50    CONTINUE
*
         NORMDIFV = PDLANGE( '1', SIZE, N, VTC, IVT, JVT, DESCVT, WORK )
*
         IF( NORMDIFV.GE.ACCUR )
     $      RESULT( 8 ) = 1
*
         IF( NORMDIFV.EQ.0 .AND. ACCUR.EQ.0 ) THEN
            NORMDIFV = 0
         ELSE
            NORMDIFV = NORMDIFV / ACCUR
         END IF
*
      ELSE IF( JOBTYPE.EQ.4 ) THEN
*
         RESULT( 9 ) = RESULTS
*
      END IF
*
      CMP = MAX( NORMDIFV, NORMDIFU )
      DELTA = MAX( CMP, NORMDIFS )
*
   60 CONTINUE
*
*     End of PDSVDCMP
*
      RETURN
      END
