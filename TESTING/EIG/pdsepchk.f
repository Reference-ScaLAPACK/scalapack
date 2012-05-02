*
*
      SUBROUTINE PDSEPCHK( MS, NV, A, IA, JA, DESCA, EPSNORMA, THRESH,
     $                     Q, IQ, JQ, DESCQ, C, IC, JC, DESCC, W, WORK,
     $                     LWORK, TSTNRM, RESULT )
*
*  -- ScaLAPACK routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      INTEGER            IA, IC, IQ, JA, JC, JQ, LWORK, MS, NV, RESULT
      DOUBLE PRECISION   EPSNORMA, THRESH, TSTNRM
*     ..
*     .. Array Arguments ..
*
      INTEGER            DESCA( * ), DESCC( * ), DESCQ( * )
      DOUBLE PRECISION   A( * ), C( * ), Q( * ), W( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  Compute |AQ- QL| / (EPSNORMA * N)
*  where EPSNORMA = (abstol + eps)*norm(A) when called by pdsqpsubtst.
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
*     MP = number of local rows in A, C and Q
*     MQ = number of local columns in A
*     NQ = number of local columns in C and Q
*
*  MS      (global input) INTEGER
*          Matrix size.
*          The number of global rows in A, C and Q
*          Also, the number of global columns in A
*
*  NV      (global input) INTEGER
*          Number of eigenvectors
*          The number of global columns in C and Q.
*
*  A       (local input) DOUBLE PRECISION pointer to an
*          array in local memory of dimension (LLD_A, LOCc(JA+N-1)).
*          This array contains the local pieces of the MS-by-MS
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
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  EPSNORMA (input) DOUBLE PRECISION
*          abstol + eps * inf.norm(A)
*          Abstol is absolute tolerence for the eigenvalues and is set
*          in the calling routines, pdsepsubtst and pdsqpsubtst.
*
*  THRESH  (input) DOUBLE PRECISION
*          A test will count as "failed" if the "error", computed as
*          described below, exceeds THRESH.  Note that the error
*          is scaled to be O(1), so THRESH should be a reasonably
*          small multiple of 1, e.g., 10 or 100.  In particular,
*          it should not depend on the precision (single vs. double)
*          or the size of the matrix.  It must be at least zero.
*
*  Q       (local input) DOUBLE PRECISION array
*          global dimension (MS, NV), local dimension (DESCA(DLEN_), NQ)
*
*          Contains the eigenvectors as computed by PDSYEVX
*
*  IQ      (global input) INTEGER
*          Q's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JQ      (global input) INTEGER
*          Q's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Q.
*
*  C       (local workspace) DOUBLE PRECISION array,
*          global dimension (NV, NV), local dimension (DESCA(DLEN_), MQ)
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
*  DESCC   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix C.
*
*  W       (global input) DOUBLE PRECISION array, dimension (NV)
*
*          Contains the computed eigenvalues
*
*  WORK    (local workspace) DOUBLE PRECISION array,
*                            dimension (LWORK)
*
*  LWORK   (local input) INTEGER
*          The length of the array WORK.
*          LWORK >= NUMROC( NV, DESCA( NB_ ), MYCOL, 0, NPCOL )
*
*  TSTNRM  (global output) DOUBLE PRECISION
*          |AQ- QL| / ( EPSNROMA * MS )
*
*  RESULT  (global output) INTEGER
*          0 if the test passes i.e.
*            |AQ -QL| / (abstol + eps * norm(A) ) <= n* THRESH
*          1 if the test fails  i.e.
*            |AQ -QL| / (abstol + eps * norm(A) ) > n * THRESH
*
*     .. Local Scalars ..
*
      INTEGER            INFO, J, LOCALCOL, MP, MYCOL, MYROW, NPCOL,
     $                   NPROW, NQ, PCOL
      DOUBLE PRECISION   NORM
*     ..
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE, NEGONE
      PARAMETER          ( ONE = 1.0D+0, NEGONE = -1.0D+0 )
*     ..
*     .. External Functions ..
      INTEGER            INDXG2L, INDXG2P, NUMROC
      DOUBLE PRECISION   PDLANGE
      EXTERNAL           INDXG2L, INDXG2P, NUMROC, PDLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DLACPY, DSCAL, PDGEMM,
     $                   PXERBLA
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
      CALL CHK1MAT( MS, 1, MS, 1, IA, JA, DESCA, 6, INFO )
      CALL CHK1MAT( MS, 1, NV, 2, IQ, JQ, DESCQ, 12, INFO )
      CALL CHK1MAT( MS, 1, NV, 2, IC, JC, DESCC, 16, INFO )
*
      IF( INFO.EQ.0 ) THEN
*
         MP = NUMROC( MS, DESCA( MB_ ), MYROW, 0, NPROW )
         NQ = NUMROC( NV, DESCA( NB_ ), MYCOL, 0, NPCOL )
*
         IF( IQ.NE.1 ) THEN
            INFO = -10
         ELSE IF( JQ.NE.1 ) THEN
            INFO = -11
         ELSE IF( IA.NE.1 ) THEN
            INFO = -4
         ELSE IF( JA.NE.1 ) THEN
            INFO = -5
         ELSE IF( IC.NE.1 ) THEN
            INFO = -14
         ELSE IF( JC.NE.1 ) THEN
            INFO = -15
         ELSE IF( LWORK.LT.NQ ) THEN
            INFO = -19
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PDSEPCHK', -INFO )
         RETURN
      END IF
*
*     C = Q * W
*
      CALL DLACPY( 'A', MP, NQ, Q, DESCQ( LLD_ ), C, DESCC( LLD_ ) )
*
*
      DO 10 J = 1, NV
         PCOL = INDXG2P( J, DESCC( NB_ ), 0, 0, NPCOL )
         LOCALCOL = INDXG2L( J, DESCC( NB_ ), 0, 0, NPCOL )
*
         IF( MYCOL.EQ.PCOL ) THEN
            CALL DSCAL( MP, W( J ), C( ( LOCALCOL-1 )*DESCC( LLD_ )+1 ),
     $                  1 )
         END IF
   10 CONTINUE
*
*
*     C = C - A * Q
*
      CALL PDGEMM( 'N', 'N', MS, NV, MS, NEGONE, A, 1, 1, DESCA, Q, 1,
     $             1, DESCQ, ONE, C, 1, 1, DESCC )
*
*     Compute the norm of C
*
*
      NORM = PDLANGE( 'M', MS, NV, C, 1, 1, DESCC, WORK )
*
*
      TSTNRM = NORM / EPSNORMA / MAX( MS, 1 )
*
      IF( TSTNRM.GT.THRESH .OR. ( TSTNRM-TSTNRM.NE.0.0D0 ) ) THEN
         RESULT = 1
      END IF
*
*
      RETURN
*
*     End of PDSEPCHK
*
      END
