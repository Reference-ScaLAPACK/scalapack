      SUBROUTINE PDSVDCHK( M, N, A, IA, JA, DESCA, U, IU, JU, DESCU, VT,
     $                     IVT, JVT, DESCVT, S, THRESH, WORK, LWORK,
     $                     RESULT, CHK, MTM )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, IU, IVT, JA, JU, JVT, LWORK, M, N
      DOUBLE PRECISION   CHK, MTM, THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCU( * ), DESCVT( * ),
     $                   RESULT( * )
      DOUBLE PRECISION   A( * ), S( * ), U( * ), VT( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  For given two-dimensional matrices A, U, VT, and one-dimensional
*  array D compute the following four tests:
*
*  (1)   | A - U*diag(S) VT | / ( |A| max(M,N) ulp )
*
*  (2)   | I - U'*U | / ( M ulp )
*
*  (3)   | I - VT*VT' | / ( N ulp ),
*
*  (4)   S contains SIZE = MIN( M, N )  nonnegative values in
*   decreasing order.
*   It then compares result of computations (1)-(3)
*   with TRESH and returns results of comparisons and test (4) in
*   RESULT(I). When the i-th test fails, value of RESULT( I ) is set
*   to 1.
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
*     MP = number of local rows in A and  U
*     NQ = number of local columns in A and VT
*     SIZEP = number of local rows in VT
*     SIZEQ = number of local columns in U
*
*  M      (global input) INTEGER
*          Matrix size.
*          The number of global rows in A and U and
*
*  N      (global input) INTEGER
*          The number of global columns in A and VT.
*
*  A       (input) block cyclic distributed DOUBLE PRECISION array,
*          global dimension (M, N), local dimension (DESCA( DLEN_ ), NQ)
*          Contains the original test matrix.
*
*  IA      (global input) INTEGER
*          The global row index of the submatrix of the distributed
*          matrix A to operate on.
*
*  JA      (global input) INTEGER
*          The global column index of the submatrix of the distributed
*          matrix A to operate on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_
*          The array descriptor for the distributed matrix A.
*
*  U       (local input) DOUBLE PRECISION array
*           global dimension (M, SIZE), local dimension
*          (DESCU( DLEN_ ), SIZEQ)
*           Contains left singular vectors of matrix A.
*
*  IU      (global input) INTEGER
*          The global row index of the submatrix of the distributed
*          matrix U to operate on.
*
*  JU      (global input) INTEGER
*          The global column index of the submatrix of the distributed
*          matrix U to operate on.
*
*  DESCU   (global and local input) INTEGER array of dimension DLEN_
*          The array descriptor for the distributed matrix U.
*
*  VT       (local input) DOUBLE PRECISION array
*           global dimension (SIZE, N), local dimension
*           (DESCVT( DLEN_ ), NQ)
*           Contains right singular vectors of matrix A.
*
*  IVT     (global input) INTEGER
*          The global row index of the submatrix of the distributed
*          matrix VT to operate on.
*
*  JVT      (global input) INTEGER
*          The global column index of the submatrix of the distributed
*          matrix VT to operate on.
*
*  DESCVT   (global and local input) INTEGER array of dimension DLEN_
*          The array descriptor for the distributed matrix VT.
*
*  S       (global input) DOUBLE PRECISION array, dimension (SIZE)
*          Contains the computed singular values
*
*  THRESH  (input) DOUBLE PRECISION
*          A test will count as "failed" if the "error", computed as
*          described below, exceeds THRESH.  Note that the error
*          is scaled to be O(1), so THRESH should be a reasonably
*          small multiple of 1, e.g., 10 or 100.  In particular,
*          it should not depend on the precision (single vs. double)
*          or the size of the matrix.  It must be at least zero.
*
*  WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (local input) INTEGER
*          The length of the array WORK.
*          LWORK >= 1 + SIZEQ*SIZEP + MAX[WORK(pdlange(size,size)),
*          WORK(pdlange(m,n))],
*          where
*          SIZEQ = NUMROC( SIZE, DESCU( NB_ ), MYCOL, 0, NPCOL ),
*          SIZEP = NUMROC( SIZE, DESCVT( MB_ ), MYROW, 0, NPROW ),
*          and worekspaces required to call pdlange are
*          WORK(pdlange(size,size)) < MAX(SIZEQ0,2) < SIZEB +2,
*          WORK(pdlange(m,n)) < MAX(NQ0,2) < SIZEB +2,
*          SIZEB = MAX(M, N)
*          Finally, upper limit on required workspace is
*          LWORK >  1 + SIZEQ*SIZEP + SIZEB + 2
*
*  RESULT  (global input/output) INTEGER array. Four first elements of
*          the array are set to 0 or 1 depending on passing four
*          respective tests ( see above in Purpose ). The elements of
*          RESULT are set to
*          0 if the test passes i.e.
*            | A - U*diag(S)*VT | / ( |A| max(M,N) ulp ) <= THRESH
*          1 if the test fails  i.e.
*            | A - U*diag(S)*VT | / ( |A| max(M,N) ulp ) >  THRESH
*
*  CHK     (global output) DOUBLE PRECISION
*           value of the | A - U*diag(S) VT | / ( |A| max(M,N) ulp )
*
*  MTM     (global output) DOUBLE PRECISION
*           maximum of the two values:
*            | I - U'*U | / ( M ulp ) and  | I - VT*VT' | / ( N ulp )
*
* ======================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE, MONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, MONE = -1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, LDR, LOCALCOL, LWMIN, MP, MX, MYCOL,
     $                   MYROW, NPCOL, NPROW, NQ, PCOL, PTRR, PTRWORK,
     $                   SIZE, SIZEP, SIZEPOS, SIZEQ
      DOUBLE PRECISION   FIRST, NORMA, NORMAI, NORMU, NORMVT, SECOND,
     $                   THRESHA, ULP
*     ..
*     .. Local Arrays ..
      INTEGER            DESCR( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            INDXG2L, INDXG2P, NUMROC
      DOUBLE PRECISION   PDLAMCH, PDLANGE
      EXTERNAL           INDXG2L, INDXG2P, NUMROC, PDLAMCH, PDLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCINIT, DSCAL,
     $                   PDELSET, PDGEMM, PDLASET, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*     This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*DTYPE_*M_*N_*RSRC_.LT.0 ) RETURN
*
*     Test the input parameters.
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
      SIZE = MIN( M, N )
*
*     Sizepos is a number of parameters to pdsvdchk plus one. It's used
*     for the error reporting.
*
      SIZEPOS = 22
      IF( NPROW.EQ.-1 ) THEN
         INFO = -607
      ELSE
         CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, INFO )
         CALL CHK1MAT( M, 1, SIZE, SIZEPOS, IU, JU, DESCU, 10, INFO )
         CALL CHK1MAT( SIZE, SIZEPOS, N, 2, IVT, JVT, DESCVT, 14, INFO )
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*     Calculate workspace
*
         MP = NUMROC( M, DESCA( MB_ ), MYROW, 0, NPROW )
         NQ = NUMROC( N, DESCA( NB_ ), MYCOL, 0, NPCOL )
         SIZEP = NUMROC( SIZE, DESCVT( MB_ ), MYROW, 0, NPROW )
         SIZEQ = NUMROC( SIZE, DESCU( NB_ ), MYCOL, 0, NPCOL )
         MX = MAX( SIZEQ, NQ )
         LWMIN = 2 + SIZEQ*SIZEP + MAX( 2, MX )
         WORK( 1 ) = LWMIN
         IF( LWORK.EQ.-1 )
     $      GO TO 40
         IF( LWORK.LT.LWMIN ) THEN
            INFO = -18
         ELSE IF( THRESH.LE.0 ) THEN
            INFO = -16
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PDSVDCHK', -INFO )
         RETURN
      END IF
*
      LDR = MAX( 1, SIZEP )
      ULP = PDLAMCH( DESCA( CTXT_ ), 'P' )
      NORMAI = PDLANGE( '1', M, N, A, IA, JA, DESCA, WORK )
*
*     Allocate array R of global dimension SIZE x SIZE for testing
*
      PTRR = 2
      PTRWORK = PTRR + SIZEQ*SIZEP
*
      CALL DESCINIT( DESCR, SIZE, SIZE, DESCVT( MB_ ), DESCU( NB_ ), 0,
     $               0, DESCA( CTXT_ ), LDR, INFO )
*
*     Test 2. Form identity matrix R  and make  check norm(U'*U - I )
*
      CALL PDLASET( 'Full', SIZE, SIZE, ZERO, ONE, WORK( PTRR ), 1, 1,
     $              DESCR )
      CALL PDGEMM( 'T', 'N', SIZE, SIZE, M, ONE, U, IU, JU, DESCU, U, 
     $             IU, JU, DESCU, MONE, WORK( PTRR ), 1, 1, DESCR )
*
      NORMU = PDLANGE( '1', SIZE, SIZE, WORK( PTRR ), 1, 1, DESCR,
     $        WORK( PTRWORK ) )
*
      NORMU = NORMU / ULP / SIZE / THRESH
      IF( NORMU.GT.1. )
     $   RESULT( 2 ) = 1
*
*     Test3. Form identity matrix R  and check norm(VT*VT' - I )
*
      CALL PDLASET( 'Full', SIZE, SIZE, ZERO, ONE, WORK( PTRR ), 1, 1,
     $              DESCR )
      CALL PDGEMM( 'N', 'T', SIZE, SIZE, N, ONE, VT, IVT, JVT, DESCVT,
     $             VT, IVT, JVT, DESCVT, MONE, WORK( PTRR ), 
     $             1, 1, DESCR )
      NORMVT = PDLANGE( '1', SIZE, SIZE, WORK( PTRR ), 1, 1, DESCR,
     $         WORK( PTRWORK ) )
*
      NORMVT = NORMVT / ULP / SIZE / THRESH
      IF( NORMVT.GT.1. )
     $   RESULT( 3 ) = 1
*
      MTM = MAX( NORMVT, NORMU )*THRESH
*
*     Test 1.
*     Initialize R = diag( S )
*
      CALL PDLASET( 'Full', SIZE, SIZE, ZERO, ZERO, WORK( PTRR ), 1, 1,
     $              DESCR )
*
      DO 10 I = 1, SIZE
         CALL PDELSET( WORK( PTRR ), I, I, DESCR, S( I ) )
   10 CONTINUE
*
*     Calculate U = U*R
*
      DO 20 I = 1, SIZE
         PCOL = INDXG2P( I, DESCU( NB_ ), 0, 0, NPCOL )
         LOCALCOL = INDXG2L( I, DESCU( NB_ ), 0, 0, NPCOL )
         IF( MYCOL.EQ.PCOL ) THEN
            CALL DSCAL( MP, S( I ), U( ( LOCALCOL-1 )*DESCU( LLD_ )+1 ),
     $                  1 )
         END IF
   20 CONTINUE
*
*     Calculate A = U*VT - A
*
      CALL PDGEMM( 'N', 'N', M, N, SIZE, ONE, U, IU, JU, DESCU, VT,
     $             IVT, JVT, DESCVT, MONE, A, IA, JA, DESCA )
*
      NORMA = PDLANGE( '1', M, N, A, IA, JA, DESCA, WORK( PTRWORK ) )
      THRESHA = NORMAI*MAX( M, N )*ULP*THRESH
*
      IF( NORMA.GT.THRESHA )
     $   RESULT( 1 ) = 1
*
      IF( THRESHA.EQ.0 ) THEN
         CHK = 0.0D0
      ELSE
         CHK = NORMA / THRESHA*THRESH
      END IF
*
*     Test 4.
*
      DO 30 I = 1, SIZE - 1
         FIRST = S( I )
         SECOND = S( I+1 )
         IF( FIRST.LT.SECOND )
     $      RESULT( 4 ) = 1
   30 CONTINUE
   40 CONTINUE
      RETURN
      END
