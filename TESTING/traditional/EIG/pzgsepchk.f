*
*
      SUBROUTINE PZGSEPCHK( IBTYPE, MS, NV, A, IA, JA, DESCA, B, IB, JB,
     $                      DESCB, THRESH, Q, IQ, JQ, DESCQ, C, IC, JC,
     $                      DESCC, W, WORK, LWORK, TSTNRM, RESULT )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 15, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, IB, IBTYPE, IC, IQ, JA, JB, JC, JQ, LWORK,
     $                   MS, NV, RESULT
      DOUBLE PRECISION   THRESH, TSTNRM
*     ..
*     .. Array Arguments ..
*
      INTEGER            DESCA( * ), DESCB( * ), DESCC( * ), DESCQ( * )
      DOUBLE PRECISION   W( * ), WORK( * )
      COMPLEX*16         A( * ), B( * ), C( * ), Q( * )
*     ..
*
*
*  Purpose
*  =======
*
*  PZGSEPCHK checks a decomposition of the form
*
*     A Q   =  B Q D or
*     A B Q =  Q D or
*     B A Q =  Q D
*
*  where A is a symmetric matrix, B is
*  symmetric positive definite, Q is orthogonal, and D is diagonal.
*
*  One of the following test ratios is computed:
*
*  IBTYPE = 1:  TSTNRM = | A Q - B Q D | / ( |A| |Q| n ulp )
*
*  IBTYPE = 2:  TSTNRM = | A B Q - Q D | / ( |A| |Q| n ulp )
*
*  IBTYPE = 3:  TSTNRM = | B A Q - Q D | / ( |A| |Q| n ulp )
*
*
*  Notes
*  =====
*
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
*
*  Arguments
*  =========
*
*     MP = number of local rows in A, B and Q
*     MQ = number of local columns in A
*     NQ = number of local columns in B and Q
*
*  IBTYPE   (input) INTEGER
*          The form of the symmetric generalized eigenproblem.
*          = 1:  A*Q = (lambda)*B*Q
*          = 2:  A*B*Q = (lambda)*Q
*          = 3:  B*A*Q = (lambda)*Q
*
*  MS      (global input) INTEGER
*          Matrix size.
*          The number of global rows in A, B, C and Q
*          Also, the number of columns in A
*
*  NV      (global input) INTEGER
*          Number of eigenvectors
*          The number of global columns in C and Q.
*
*  A       (local input) REAL pointer to an
*          array in local memory of dimension (LLD_A, LOCc(JA+N-1)).
*          This array contains the local pieces of the M-by-N
*          distributed test matrix A
*
*  IA      (global input) INTEGER
*          A's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCA   (global and local input) INTEGER array of dimension 8
*          The array descriptor for the distributed matrix A.
*
*  B       (local input) REAL pointer to an
*          array in local memory of dimension (LLD_B, LOCc(JB+N-1)).
*          This array contains the local pieces of the M-by-N
*          distributed test matrix B
*
*  IB      (global input) INTEGER
*          B's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JB      (global input) INTEGER
*          B's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCB   (global and local input) INTEGER array of dimension 8
*          The array descriptor for the distributed matrix B.
*
*  THRESH  (input) REAL
*          A test will count as "failed" if the "error", computed as
*          described below, exceeds THRESH.  Note that the error
*          is scaled to be O(1), so THRESH should be a reasonably
*          small multiple of 1, e.g., 10 or 100.  In particular,
*          it should not depend on the precision (single vs. double)
*          or the size of the matrix.  It must be at least zero.
*
*  Q       (local input) REAL array
*          global dimension (MS, NV),
*          local dimension (DESCA( DLEN_ ), NQ)
*
*          Contains the eigenvectors as computed by PSSYEVX
*
*  IQ      (global input) INTEGER
*          Q's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JQ      (global input) INTEGER
*          Q's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCQ   (global and local input) INTEGER array of dimension 8
*          The array descriptor for the distributed matrix Q.
*
*  C       (local workspace) REAL array,
*          global dimension (MS, NV),
*          local dimension (DESCA( DLEN_ ), MQ)
*
*          Accumulator for computing AQ -QL
*
*  IC      (global input) INTEGER
*          C's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JC      (global input) INTEGER
*          C's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCC   (global and local input) INTEGER array of dimension 8
*          The array descriptor for the distributed matrix C.
*
*  W       (global input) REAL array, dimension (NV)
*
*          Contains the computed eigenvalues
*
*  WORK    (local workspace) REAL array,
*                            dimension (LWORK)
*
*  LWORK   (local input) INTEGER
*          The length of the array WORK.
*          LWORK >= NUMROC( NV, DESCA( NB_ ), MYCOL, 0, NPCOL )
*
*  TSTNRM  (global output) REAL
*
*  RESULT  (global output) INTEGER
*          0 if the test passes
*          1 if the test fails
*
*     .. Local Scalars ..
*
      INTEGER            I, INFO, MYCOL, MYROW, NPCOL, NPROW, NQ
      DOUBLE PRECISION   ANORM, ULP
*     ..
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16         CONE, CNEGONE, CZERO
      PARAMETER          ( CONE = 1.0D+0, CNEGONE = -1.0D+0,
     $                   CZERO = 0.0D+0 )
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      DOUBLE PRECISION   DLAMCH, PZLANGE
      EXTERNAL           NUMROC, DLAMCH, PZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PXERBLA, PZDSCAL,
     $                   PZGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      RESULT = 0
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
      INFO = 0
      CALL CHK1MAT( MS, 1, MS, 2, IA, JA, DESCA, 7, INFO )
      CALL CHK1MAT( MS, 1, MS, 2, IB, JB, DESCB, 11, INFO )
      CALL CHK1MAT( MS, 1, NV, 2, IQ, JQ, DESCQ, 16, INFO )
      CALL CHK1MAT( MS, 1, NV, 2, IB, JB, DESCB, 20, INFO )
*
      IF( INFO.EQ.0 ) THEN
*
         NQ = NUMROC( NV, DESCA( NB_ ), MYCOL, 0, NPCOL )
*
         IF( IQ.NE.1 ) THEN
            INFO = -14
         ELSE IF( JQ.NE.1 ) THEN
            INFO = -15
         ELSE IF( IA.NE.1 ) THEN
            INFO = -5
         ELSE IF( JA.NE.1 ) THEN
            INFO = -6
         ELSE IF( IB.NE.1 ) THEN
            INFO = -9
         ELSE IF( JB.NE.1 ) THEN
            INFO = -10
         ELSE IF( LWORK.LT.NQ ) THEN
            INFO = -23
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PZGSEPCHK', -INFO )
         RETURN
      END IF
*
      RESULT = 0
      ULP = DLAMCH( 'Epsilon' )
*
*     Compute product of Max-norms of A and Q.
*
      ANORM = PZLANGE( 'M', MS, MS, A, IA, JA, DESCA, WORK )*
     $        PZLANGE( 'M', MS, NV, Q, IQ, JQ, DESCQ, WORK )
      IF( ANORM.EQ.ZERO )
     $   ANORM = ONE
*
      IF( IBTYPE.EQ.1 ) THEN
*
*        Norm of AQ - BQD
*
*        C = AQ
*
         CALL PZGEMM( 'N', 'N', MS, NV, MS, CONE, A, IA, JA, DESCA, Q,
     $                IQ, JQ, DESCQ, CZERO, C, IC, JC, DESCC )
*
*        Q = QD
*
         DO 10 I = 1, NV
            CALL PZDSCAL( MS, W( I ), Q, IQ, JQ+I-1, DESCQ, 1 )
   10    CONTINUE
*
*        C = C - BQ  (i.e. AQ-BQD)
*
         CALL PZGEMM( 'N', 'N', MS, NV, MS, CONE, B, IB, JB, DESCB, Q,
     $                IQ, JQ, DESCQ, CNEGONE, C, IC, JC, DESCC )
*
         TSTNRM = ( PZLANGE( 'M', MS, NV, C, IC, JC, DESCC, WORK ) /
     $            ANORM ) / ( MAX( MS, 1 )*ULP )
*
*
      ELSE IF( IBTYPE.EQ.2 ) THEN
*
*        Norm of ABQ - QD
*
*
*        C = BQ
*
         CALL PZGEMM( 'N', 'N', MS, NV, MS, CONE, B, IB, JB, DESCB, Q,
     $                IQ, JQ, DESCQ, CZERO, C, IC, JC, DESCC )
*
*        Q = QD
*
         DO 20 I = 1, NV
            CALL PZDSCAL( MS, W( I ), Q, IQ, JQ+I-1, DESCQ, 1 )
   20    CONTINUE
*
*        Q = AC - Q
*
         CALL PZGEMM( 'N', 'N', MS, NV, MS, CONE, A, IA, JA, DESCA, C,
     $                IC, JC, DESCC, CNEGONE, Q, IQ, JQ, DESCQ )
*
         TSTNRM = ( PZLANGE( 'M', MS, NV, Q, IQ, JQ, DESCQ, WORK ) /
     $            ANORM ) / ( MAX( MS, 1 )*ULP )
*
      ELSE IF( IBTYPE.EQ.3 ) THEN
*
*        Norm of BAQ - QD
*
*
*        C = AQ
*
         CALL PZGEMM( 'N', 'N', MS, NV, MS, CONE, A, IA, JA, DESCA, Q,
     $                IQ, JQ, DESCQ, CZERO, C, IC, JC, DESCC )
*
*        Q = QD
*
         DO 30 I = 1, NV
            CALL PZDSCAL( MS, W( I ), Q, IQ, JQ+I-1, DESCQ, 1 )
   30    CONTINUE
*
*        Q = BC - Q
*
         CALL PZGEMM( 'N', 'N', MS, NV, MS, CONE, B, IB, JB, DESCB, C,
     $                IC, JC, DESCC, CNEGONE, Q, IQ, JQ, DESCQ )
*
         TSTNRM = ( PZLANGE( 'M', MS, NV, Q, IQ, JQ, DESCQ, WORK ) /
     $            ANORM ) / ( MAX( MS, 1 )*ULP )
*
      END IF
*
      IF( TSTNRM.GT.THRESH .OR. ( TSTNRM-TSTNRM.NE.0.0D0 ) ) THEN
         RESULT = 1
      END IF
      RETURN
*
*     End of PZGSEPCHK
*
      END
