*
*
      SUBROUTINE PDSEPQTQ( MS, NV, THRESH, Q, IQ, JQ, DESCQ, C, IC, JC,
     $                     DESCC, PROCDIST, ICLUSTR, GAP, WORK, LWORK,
     $                     QTQNRM, INFO, RES )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IC, INFO, IQ, JC, JQ, LWORK, MS, NV, RES
      DOUBLE PRECISION   QTQNRM, THRESH
*     ..
*     .. Array Arguments ..
*
      INTEGER            DESCC( * ), DESCQ( * ), ICLUSTR( * ),
     $                   PROCDIST( * )
      DOUBLE PRECISION   C( * ), GAP( * ), Q( * ), WORK( * )
*     ..
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
*  Purpose
*  =======
*
*  Compute |I - QT * Q| / (ulp * n)
*
*  Arguments
*  =========
*
*     NP = number of local rows in C
*     NQ = number of local columns in C and Q
*
*  MS      (global input) INTEGER
*          Matrix size.
*          The number of global rows in Q
*
*  NV      (global input) INTEGER
*          Number of eigenvectors
*          The number of global columns in C and Q
*
*  THRESH  (global input) DOUBLE PRECISION
*          A test will count as "failed" if the "error", computed as
*          described below, exceeds THRESH.  Note that the error
*          is scaled to be O(1), so THRESH should be a reasonably
*          small multiple of 1, e.g., 10 or 100.  In particular,
*          it should not depend on the precision (single vs. double)
*          or the size of the matrix.  It must be at least zero.
*
*  Q       (local input) DOUBLE PRECISION array,
*          global dimension (MS, NV), local dimension (LDQ, NQ)
*
*          Contains the eigenvectors as computed by PDSTEIN
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
*  C       (local workspace)  DOUBLE PRECISION array,
*          global dimension (NV, NV), local dimension (DESCC(DLEN_), NQ)
*
*          Accumulator for computing I - QT * Q
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
*  W       (input) DOUBLE PRECISION array, dimension (NV)
*          All procesors have an identical copy of W()
*
*          Contains the computed eigenvalues
*
*  PROCDIST (global input) INTEGER array dimension (NPROW*NPCOL+1)
*          Identifies which eigenvectors are the last to be computed
*          by a given process
*
*  ICLUSTR (global input) INTEGER array dimension (2*P)
*          This input array contains indices of eigenvectors
*          corresponding to a cluster of eigenvalues that could not be
*          orthogonalized due to insufficient workspace.
*          This should be the output of PDSTEIN.
*
*  GAP     (global input) DOUBLE PRECISION array, dimension (P)
*          This input array contains the gap between eigenvalues whose
*          eigenvectors could not be orthogonalized.
*
*  WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (local input) INTEGER
*          The length of the array WORK.
*          LWORK >= 2 + MAX( DESCC( MB_ ), 2 )*( 2*NP0+MQ0 )
*          Where:
*          NP0 = NUMROC( NV, DESCC( MB_ ), 0, 0, NPROW )
*          MQ0 = NUMROC( NV, DESCC( NB_ ), 0, 0, NPCOL )
*
*  QTQNRM  (global output) DOUBLE PRECISION
*          |QTQ -I| / EPS
*
*  RES     (global output) INTEGER
*          0 if the test passes i.e. |I - QT * Q| / (ulp * n) <= THRESH
*          1 if the test fails  i.e. |I - QT * Q| / (ulp * n) > THRESH
*
*
*     .. Parameters ..
*
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE, NEGONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   NEGONE = -1.0D+0 )
*     ..
*     .. Intrinsic Functions ..
*
      INTRINSIC          DBLE, MAX
*     ..
*     .. Local Scalars ..
      INTEGER            CLUSTER, FIRSTP, IMAX, IMIN, JMAX, JMIN, LWMIN,
     $                   MQ0, MYCOL, MYROW, NEXTP, NP0, NPCOL, NPROW
      DOUBLE PRECISION   NORM, QTQNRM2, ULP
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      DOUBLE PRECISION   PDLAMCH, PDLANGE
      EXTERNAL           NUMROC, PDLAMCH, PDLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PDGEMM, PDLASET,
     $                   PDMATADD, PXERBLA
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*
      RES = 0
      ULP = PDLAMCH( DESCC( CTXT_ ), 'P' )
*
      CALL BLACS_GRIDINFO( DESCC( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
      CALL CHK1MAT( MS, 1, MS, 2, IQ, JQ, DESCQ, 7, INFO )
      CALL CHK1MAT( NV, 1, MS, 2, IC, JC, DESCC, 11, INFO )
*
      IF( INFO.EQ.0 ) THEN
         NP0 = NUMROC( NV, DESCC( MB_ ), 0, 0, NPROW )
         MQ0 = NUMROC( NV, DESCC( NB_ ), 0, 0, NPCOL )
*
         LWMIN = 2 + MAX( DESCC( MB_ ), 2 )*( 2*NP0+MQ0 )
*
         IF( IQ.NE.1 ) THEN
            INFO = -5
         ELSE IF( JQ.NE.1 ) THEN
            INFO = -6
         ELSE IF( IC.NE.1 ) THEN
            INFO = -9
         ELSE IF( JC.NE.1 ) THEN
            INFO = -10
         ELSE IF( LWORK.LT.LWMIN ) THEN
            INFO = -16
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCC( CTXT_ ), 'PDSEPQTQ', -INFO )
         RETURN
      END IF
*
*     C = Identity matrix
*
      CALL PDLASET( 'A', NV, NV, ZERO, ONE, C, IC, JC, DESCC )
*
*     C = C - QT * Q
*
      IF( NV*MS.GT.0 ) THEN
         CALL PDGEMM( 'Transpose', 'N', NV, NV, MS, NEGONE, Q, 1, 1,
     $                DESCQ, Q, 1, 1, DESCQ, ONE, C, 1, 1, DESCC )
      END IF
*
*     Allow for poorly orthogonalized eigenvectors for large clusters
*
      NORM = PDLANGE( '1', NV, NV, C, 1, 1, DESCC, WORK )
      QTQNRM = NORM / ( DBLE( MAX( MS, 1 ) )*ULP )
*
      CLUSTER = 1
   10 CONTINUE
      DO 20 FIRSTP = 1, NPROW*NPCOL
         IF( PROCDIST( FIRSTP ).GE.ICLUSTR( 2*( CLUSTER-1 )+1 ) )
     $      GO TO 30
   20 CONTINUE
   30 CONTINUE
*
      IMIN = ICLUSTR( 2*CLUSTER-1 )
      JMAX = ICLUSTR( 2*CLUSTER )
*
*
      IF( IMIN.EQ.0 )
     $   GO TO 60
*
      DO 40 NEXTP = FIRSTP, NPROW*NPCOL
         IMAX = PROCDIST( NEXTP )
         JMIN = IMAX + 1
*
*
         CALL PDMATADD( IMAX-IMIN+1, JMAX-JMIN+1, ZERO, C, IMIN, JMIN,
     $                  DESCC, GAP( CLUSTER ) / 0.01D+0, C, IMIN, JMIN,
     $                  DESCC )
         CALL PDMATADD( JMAX-JMIN+1, IMAX-IMIN+1, ZERO, C, JMIN, IMIN,
     $                  DESCC, GAP( CLUSTER ) / 0.01D+0, C, JMIN, IMIN,
     $                  DESCC )
         IMIN = IMAX
*
         IF( ICLUSTR( 2*CLUSTER ).LT.PROCDIST( NEXTP+1 ) )
     $      GO TO 50
   40 CONTINUE
   50 CONTINUE
*
      CLUSTER = CLUSTER + 1
      GO TO 10
   60 CONTINUE
*
*     Compute the norm of C
*
      NORM = PDLANGE( '1', NV, NV, C, 1, 1, DESCC, WORK )
*
      QTQNRM2 = NORM / ( DBLE( MAX( MS, 1 ) )*ULP )
*
      IF( QTQNRM2.GT.THRESH ) THEN
         RES = 1
         QTQNRM = QTQNRM2
      END IF
      RETURN
*
*     End of PDSEPQTQ
*
      END
