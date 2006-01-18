      SUBROUTINE PCGESVX( FACT, TRANS, N, NRHS, A, IA, JA, DESCA, AF,
     $                    IAF, JAF, DESCAF, IPIV, EQUED, R, C, B, IB,
     $                    JB, DESCB, X, IX, JX, DESCX, RCOND, FERR,
     $                    BERR, WORK, LWORK, RWORK, LRWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, FACT, TRANS
      INTEGER            IA, IAF, IB, INFO, IX, JA, JAF, JB, JX, LRWORK,
     $                   LWORK, N, NRHS
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCAF( * ), DESCB( * ),
     $                   DESCX( * ), IPIV( * )
      REAL               BERR( * ), C( * ), FERR( * ), R( * ),
     $                   RWORK( * )
      COMPLEX            A( * ), AF( * ), B( * ), WORK( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  PCGESVX uses the LU factorization to compute the solution to a
*  complex system of linear equations
*
*        A(IA:IA+N-1,JA:JA+N-1) * X = B(IB:IB+N-1,JB:JB+NRHS-1),
*
*  where A(IA:IA+N-1,JA:JA+N-1) is an N-by-N matrix and X and
*  B(IB:IB+N-1,JB:JB+NRHS-1) are N-by-NRHS matrices.
*
*  Error bounds on the solution and a condition estimate are also
*  provided.
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
*  Description
*  ===========
*
*  In the following description, A denotes A(IA:IA+N-1,JA:JA+N-1),
*  B denotes B(IB:IB+N-1,JB:JB+NRHS-1) and X denotes
*  X(IX:IX+N-1,JX:JX+NRHS-1).
*
*  The following steps are performed:
*
*  1. If FACT = 'E', real scaling factors are computed to equilibrate
*     the system:
*        TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
*        TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
*        TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
*     Whether or not the system will be equilibrated depends on the
*     scaling of the matrix A, but if equilibration is used, A is
*     overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
*     or diag(C)*B (if TRANS = 'T' or 'C').
*
*  2. If FACT = 'N' or 'E', the LU decomposition is used to factor the
*     matrix A (after equilibration if FACT = 'E') as
*        A = P * L * U,
*     where P is a permutation matrix, L is a unit lower triangular
*     matrix, and U is upper triangular.
*
*  3. The factored form of A is used to estimate the condition number
*     of the matrix A.  If the reciprocal of the condition number is
*     less than machine precision, steps 4-6 are skipped.
*
*  4. The system of equations is solved for X using the factored form
*     of A.
*
*  5. Iterative refinement is applied to improve the computed solution
*     matrix and calculate error bounds and backward error estimates
*     for it.
*
*  6. If FACT = 'E' and equilibration was used, the matrix X is
*     premultiplied by diag(C) (if TRANS = 'N') or diag(R) (if
*     TRANS = 'T' or 'C') so that it solves the original system
*     before equilibration.
*
*  Arguments
*  =========
*
*  FACT    (global input) CHARACTER
*          Specifies whether or not the factored form of the matrix
*          A(IA:IA+N-1,JA:JA+N-1) is supplied on entry, and if not,
*          whether the matrix A(IA:IA+N-1,JA:JA+N-1) should be
*          equilibrated before it is factored.
*          = 'F':  On entry, AF(IAF:IAF+N-1,JAF:JAF+N-1) and IPIV con-
*                  tain the factored form of A(IA:IA+N-1,JA:JA+N-1).
*                  If EQUED is not 'N', the matrix
*                  A(IA:IA+N-1,JA:JA+N-1) has been equilibrated with
*                  scaling factors given by R and C.
*                  A(IA:IA+N-1,JA:JA+N-1), AF(IAF:IAF+N-1,JAF:JAF+N-1),
*                  and IPIV are not modified.
*          = 'N':  The matrix A(IA:IA+N-1,JA:JA+N-1) will be copied to
*                  AF(IAF:IAF+N-1,JAF:JAF+N-1) and factored.
*          = 'E':  The matrix A(IA:IA+N-1,JA:JA+N-1) will be equili-
*                  brated if necessary, then copied to
*                  AF(IAF:IAF+N-1,JAF:JAF+N-1) and factored.
*
*  TRANS   (global input) CHARACTER
*          Specifies the form of the system of equations:
*          = 'N':  A(IA:IA+N-1,JA:JA+N-1) * X(IX:IX+N-1,JX:JX+NRHS-1)
*                = B(IB:IB+N-1,JB:JB+NRHS-1)     (No transpose)
*          = 'T':  A(IA:IA+N-1,JA:JA+N-1)**T * X(IX:IX+N-1,JX:JX+NRHS-1)
*                = B(IB:IB+N-1,JB:JB+NRHS-1)  (Transpose)
*          = 'C':  A(IA:IA+N-1,JA:JA+N-1)**H * X(IX:IX+N-1,JX:JX+NRHS-1)
*                = B(IB:IB+N-1,JB:JB+NRHS-1)  (Conjugate transpose)
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix A(IA:IA+N-1,JA:JA+N-1).
*          N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right-hand sides, i.e., the number of columns
*          of the distributed submatrices B(IB:IB+N-1,JB:JB+NRHS-1) and
*          X(IX:IX+N-1,JX:JX+NRHS-1).  NRHS >= 0.
*
*  A       (local input/local output) COMPLEX pointer into
*          the local memory to an array of local dimension
*          (LLD_A,LOCc(JA+N-1)).  On entry, the N-by-N matrix
*          A(IA:IA+N-1,JA:JA+N-1).  If FACT = 'F' and EQUED is not 'N',
*          then A(IA:IA+N-1,JA:JA+N-1) must have been equilibrated by
*          the scaling factors in R and/or C.  A(IA:IA+N-1,JA:JA+N-1) is
*          not modified if FACT = 'F' or  'N', or if FACT = 'E' and
*          EQUED = 'N' on exit.
*
*          On exit, if EQUED .ne. 'N', A(IA:IA+N-1,JA:JA+N-1) is scaled
*          as follows:
*          EQUED = 'R':  A(IA:IA+N-1,JA:JA+N-1) :=
*                                      diag(R) * A(IA:IA+N-1,JA:JA+N-1)
*          EQUED = 'C':  A(IA:IA+N-1,JA:JA+N-1) :=
*                                      A(IA:IA+N-1,JA:JA+N-1) * diag(C)
*          EQUED = 'B':  A(IA:IA+N-1,JA:JA+N-1) :=
*                            diag(R) * A(IA:IA+N-1,JA:JA+N-1) * diag(C).
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  AF      (local input or local output) COMPLEX pointer
*          into the local memory to an array of local dimension
*          (LLD_AF,LOCc(JA+N-1)).  If FACT = 'F', then
*          AF(IAF:IAF+N-1,JAF:JAF+N-1) is an input argument and on
*          entry contains the factors L and U from the factorization
*          A(IA:IA+N-1,JA:JA+N-1) = P*L*U as computed by PCGETRF.
*          If EQUED .ne. 'N', then AF is the factored form of the
*          equilibrated matrix A(IA:IA+N-1,JA:JA+N-1).
*
*          If FACT = 'N', then AF(IAF:IAF+N-1,JAF:JAF+N-1) is an output
*          argument and on exit returns the factors L and U from the
*          factorization A(IA:IA+N-1,JA:JA+N-1) = P*L*U of the original
*          matrix A(IA:IA+N-1,JA:JA+N-1).
*
*          If FACT = 'E', then AF(IAF:IAF+N-1,JAF:JAF+N-1) is an output
*          argument and on exit returns the factors L and U from the
*          factorization A(IA:IA+N-1,JA:JA+N-1) = P*L*U of the equili-
*          brated matrix A(IA:IA+N-1,JA:JA+N-1) (see the description of
*          A(IA:IA+N-1,JA:JA+N-1) for the form of the equilibrated
*          matrix).
*
*  IAF     (global input) INTEGER
*          The row index in the global array AF indicating the first
*          row of sub( AF ).
*
*  JAF     (global input) INTEGER
*          The column index in the global array AF indicating the
*          first column of sub( AF ).
*
*  DESCAF  (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix AF.
*
*  IPIV    (local input or local output) INTEGER array, dimension
*          LOCr(M_A)+MB_A. If FACT = 'F', then IPIV is an input argu-
*          ment and on entry contains the pivot indices from the fac-
*          torization A(IA:IA+N-1,JA:JA+N-1) = P*L*U as computed by
*          PCGETRF; IPIV(i) -> The global row local row i was
*          swapped with.  This array must be aligned with
*          A( IA:IA+N-1, * ).
*
*          If FACT = 'N', then IPIV is an output argument and on exit
*          contains the pivot indices from the factorization
*          A(IA:IA+N-1,JA:JA+N-1) = P*L*U of the original matrix
*          A(IA:IA+N-1,JA:JA+N-1).
*
*          If FACT = 'E', then IPIV is an output argument and on exit
*          contains the pivot indices from the factorization
*          A(IA:IA+N-1,JA:JA+N-1) = P*L*U of the equilibrated matrix
*          A(IA:IA+N-1,JA:JA+N-1).
*
*  EQUED   (global input or global output) CHARACTER
*          Specifies the form of equilibration that was done.
*          = 'N':  No equilibration (always true if FACT = 'N').
*          = 'R':  Row equilibration, i.e., A(IA:IA+N-1,JA:JA+N-1) has
*                  been premultiplied by diag(R).
*          = 'C':  Column equilibration, i.e., A(IA:IA+N-1,JA:JA+N-1)
*                  has been postmultiplied by diag(C).
*          = 'B':  Both row and column equilibration, i.e.,
*                  A(IA:IA+N-1,JA:JA+N-1) has been replaced by
*                  diag(R) * A(IA:IA+N-1,JA:JA+N-1) * diag(C).
*          EQUED is an input variable if FACT = 'F'; otherwise, it is an
*          output variable.
*
*  R       (local input or local output) REAL array,
*                                                  dimension LOCr(M_A).
*          The row scale factors for A(IA:IA+N-1,JA:JA+N-1).
*          If EQUED = 'R' or 'B', A(IA:IA+N-1,JA:JA+N-1) is multiplied
*          on the left by diag(R); if EQUED='N' or 'C', R is not acces-
*          sed.  R is an input variable if FACT = 'F'; otherwise, R is
*          an output variable.
*          If FACT = 'F' and EQUED = 'R' or 'B', each element of R must
*          be positive.
*          R is replicated in every process column, and is aligned
*          with the distributed matrix A.
*
*  C       (local input or local output) REAL array,
*                                                  dimension LOCc(N_A).
*          The column scale factors for A(IA:IA+N-1,JA:JA+N-1).
*          If EQUED = 'C' or 'B', A(IA:IA+N-1,JA:JA+N-1) is multiplied
*          on the right by diag(C); if EQUED = 'N' or 'R', C is not
*          accessed.  C is an input variable if FACT = 'F'; otherwise,
*          C is an output variable.  If FACT = 'F' and EQUED = 'C' or
*          'B', each element of C must be positive.
*          C is replicated in every process row, and is aligned with
*          the distributed matrix A.
*
*  B       (local input/local output) COMPLEX pointer
*          into the local memory to an array of local dimension
*          (LLD_B,LOCc(JB+NRHS-1) ).  On entry, the N-by-NRHS right-hand
*          side matrix B(IB:IB+N-1,JB:JB+NRHS-1). On exit, if
*          EQUED = 'N', B(IB:IB+N-1,JB:JB+NRHS-1) is not modified; if
*          TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
*          diag(R)*B(IB:IB+N-1,JB:JB+NRHS-1); if TRANS = 'T' or 'C'
*          and EQUED = 'C' or 'B', B(IB:IB+N-1,JB:JB+NRHS-1) is over-
*          written by diag(C)*B(IB:IB+N-1,JB:JB+NRHS-1).
*
*  IB      (global input) INTEGER
*          The row index in the global array B indicating the first
*          row of sub( B ).
*
*  JB      (global input) INTEGER
*          The column index in the global array B indicating the
*          first column of sub( B ).
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
*
*  X       (local input/local output) COMPLEX pointer
*          into the local memory to an array of local dimension
*          (LLD_X, LOCc(JX+NRHS-1)).  If INFO = 0, the N-by-NRHS
*          solution matrix X(IX:IX+N-1,JX:JX+NRHS-1) to the original
*          system of equations.  Note that A(IA:IA+N-1,JA:JA+N-1) and
*          B(IB:IB+N-1,JB:JB+NRHS-1) are modified on exit if
*          EQUED .ne. 'N', and the solution to the equilibrated system
*          is inv(diag(C))*X(IX:IX+N-1,JX:JX+NRHS-1) if TRANS = 'N'
*          and EQUED = 'C' or 'B', or
*          inv(diag(R))*X(IX:IX+N-1,JX:JX+NRHS-1) if TRANS = 'T' or 'C'
*          and EQUED = 'R' or 'B'.
*
*  IX      (global input) INTEGER
*          The row index in the global array X indicating the first
*          row of sub( X ).
*
*  JX      (global input) INTEGER
*          The column index in the global array X indicating the
*          first column of sub( X ).
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*
*  RCOND   (global output) REAL
*          The estimate of the reciprocal condition number of the matrix
*          A(IA:IA+N-1,JA:JA+N-1) after equilibration (if done).  If
*          RCOND is less than the machine precision (in particular, if
*          RCOND = 0), the matrix is singular to working precision.
*          This condition is indicated by a return code of INFO > 0.
*
*  FERR    (local output) REAL array, dimension LOCc(N_B)
*          The estimated forward error bounds for each solution vector
*          X(j) (the j-th column of the solution matrix
*          X(IX:IX+N-1,JX:JX+NRHS-1). If XTRUE is the true solution,
*          FERR(j) bounds the magnitude of the largest entry in
*          (X(j) - XTRUE) divided by the magnitude of the largest entry
*          in X(j).  The estimate is as reliable as the estimate for
*          RCOND, and is almost always a slight overestimate of the
*          true error.  FERR is replicated in every process row, and is
*          aligned with the matrices B and X.
*
*  BERR    (local output) REAL array, dimension LOCc(N_B).
*          The componentwise relative backward error of each solution
*          vector X(j) (i.e., the smallest relative change in
*          any entry of A(IA:IA+N-1,JA:JA+N-1) or
*          B(IB:IB+N-1,JB:JB+NRHS-1) that makes X(j) an exact solution).
*          BERR is replicated in every process row, and is aligned
*          with the matrices B and X.
*
*  WORK    (local workspace/local output) COMPLEX array,
*                                                    dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK = MAX( PCGECON( LWORK ), PCGERFS( LWORK ) )
*                  + LOCr( N_A ).
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  RWORK   (local workspace/local output) REAL array,
*                                                  dimension (LRWORK)
*          On exit, RWORK(1) returns the minimal and optimal LRWORK.
*
*  LRWORK  (local or global input) INTEGER
*          The dimension of the array RWORK.
*          LRWORK is local input and must be at least
*          LRWORK = 2*LOCc(N_A).
*
*          If LRWORK = -1, then LRWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, and i is
*                <= N:  U(IA+I-1,IA+I-1) is exactly zero.  The
*                       factorization has been completed, but the
*                       factor U is exactly singular, so the solution
*                       and error bounds could not be computed.
*                = N+1: RCOND is less than machine precision.  The
*                       factorization has been completed, but the
*                       matrix is singular to working precision, and
*                       the solution and error bounds have not been
*                       computed.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLEQU, EQUIL, LQUERY, NOFACT, NOTRAN, ROWEQU
      CHARACTER          NORM
      INTEGER            CONWRK, I, IACOL, IAROW, IAFROW, IBROW, IBCOL,
     $                   ICOFFA, ICOFFB, ICOFFX, ICTXT, IDUMM,
     $                   IIA, IIB, IIX,
     $                   INFEQU, IROFFA, IROFFAF, IROFFB,
     $                   IROFFX, IXCOL, IXROW, J, JJA, JJB, JJX,
     $                   LCM, LCMQ,
     $                   LRWMIN, LWMIN, MYCOL, MYROW, NP, NPCOL, NPROW,
     $                   NQ, NQB, NRHSQ, RFSWRK
      REAL               AMAX, ANORM, BIGNUM, COLCND, RCMAX, RCMIN,
     $                   ROWCND, SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER            CDESC( DLEN_ ), IDUM1( 5 ), IDUM2( 5 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, PCHK2MAT,
     $                   INFOG2L, PCGECON, PCGEEQU, PCGERFS,
     $                   PCGETRF, PCGETRS, PCLACPY,
     $                   PCLAQGE, PSCOPY, PXERBLA, SGEBR2D,
     $                   SGEBS2D, SGAMN2D, SGAMX2D
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, ILCM, INDXG2P, NUMROC
      REAL               PSLAMCH, PCLANGE
      EXTERNAL           ICEIL, ILCM, INDXG2P, LSAME, NUMROC, PCLANGE,
     $                   PSLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, MAX, MIN, MOD, REAL
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(800+CTXT_)
      ELSE
         CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 8, INFO )
         IF( LSAME( FACT, 'F' ) )
     $      CALL CHK1MAT( N, 3, N, 3, IAF, JAF, DESCAF, 12, INFO )
         CALL CHK1MAT( N, 3, NRHS, 4, IB, JB, DESCB, 20, INFO )
         CALL CHK1MAT( N, 3, NRHS, 4, IX, JX, DESCX, 24, INFO )
         NOFACT = LSAME( FACT, 'N' )
         EQUIL = LSAME( FACT, 'E' )
         NOTRAN = LSAME( TRANS, 'N' )
         IF( NOFACT .OR. EQUIL ) THEN
            EQUED = 'N'
            ROWEQU = .FALSE.
            COLEQU = .FALSE.
         ELSE
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
            SMLNUM = PSLAMCH( ICTXT, 'Safe minimum' )
            BIGNUM = ONE / SMLNUM
         END IF
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IAFROW = INDXG2P( IAF, DESCAF( MB_ ), MYROW,
     $                        DESCAF( RSRC_ ), NPROW )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $                       NPROW )
            IXROW = INDXG2P( IX, DESCX( MB_ ), MYROW, DESCX( RSRC_ ),
     $                       NPROW )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            IROFFAF = MOD( IAF-1, DESCAF( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            ICOFFB = MOD( JB-1, DESCB( NB_ ) )
            IROFFX = MOD( IX-1, DESCX( MB_ ) )
            ICOFFX = MOD( JX-1, DESCX( NB_ ) )
            CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW,
     $                    MYCOL, IIA, JJA, IAROW, IACOL )
            NP = NUMROC( N+IROFFA, DESCA( MB_ ), MYROW, IAROW,
     $                   NPROW )
            IF( MYROW.EQ.IAROW )
     $         NP = NP-IROFFA
            NQ = NUMROC( N+ICOFFA, DESCA( NB_ ), MYCOL, IACOL,
     $                   NPCOL )
            IF( MYCOL.EQ.IACOL )
     $         NQ = NQ-ICOFFA
            NQB = ICEIL( N+IROFFA, DESCA( NB_ )*NPCOL )
            LCM = ILCM( NPROW, NPCOL )
            LCMQ = LCM / NPCOL
            CONWRK = 2*NP + 2*NQ + MAX( 2, MAX( DESCA( NB_ )*
     $               MAX( 1, ICEIL( NPROW-1, NPCOL ) ), NQ +
     $               DESCA( NB_ )*
     $               MAX( 1, ICEIL( NPCOL-1, NPROW ) ) ) )
            RFSWRK = 3*NP
            IF( LSAME( TRANS, 'N' ) ) THEN
               RFSWRK = RFSWRK + NP + NQ +
     $                 ICEIL( NQB, LCMQ )*DESCA( NB_ )
            ELSE IF( LSAME( TRANS, 'T' ).OR.LSAME( TRANS, 'C' ) ) THEN
               RFSWRK = RFSWRK + NP + NQ
            END IF
            LWMIN = MAX( CONWRK, RFSWRK )
            LRWMIN = MAX( 2*NQ, NP )
            RWORK( 1 ) = REAL( LRWMIN )
            IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND.
     $          .NOT.LSAME( FACT, 'F' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND.
     $         .NOT. LSAME( TRANS, 'C' ) ) THEN
               INFO = -2
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( ICOFFA.NE.0 .OR. IROFFA.NE.ICOFFA ) THEN
               INFO = -7
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(800+NB_)
            ELSE IF( IAFROW.NE.IAROW ) THEN
               INFO = -10
            ELSE IF( IROFFAF.NE.0 ) THEN
               INFO = -10
            ELSE IF( ICTXT.NE.DESCAF( CTXT_ ) ) THEN
               INFO = -(1200+CTXT_)
            ELSE IF( LSAME( FACT, 'F' ) .AND. .NOT.
     $         ( ROWEQU .OR. COLEQU .OR. LSAME( EQUED, 'N' ) ) ) THEN
               INFO = -13
            ELSE
               IF( ROWEQU ) THEN
                  RCMIN = BIGNUM
                  RCMAX = ZERO
                  DO 10 J = IIA, IIA + NP - 1
                     RCMIN = MIN( RCMIN, R( J ) )
                     RCMAX = MAX( RCMAX, R( J ) )
   10             CONTINUE
                  CALL SGAMN2D( ICTXT, 'Columnwise', ' ', 1, 1, RCMIN,
     $                          1, IDUMM, IDUMM, -1, -1, MYCOL )
                  CALL SGAMX2D( ICTXT, 'Columnwise', ' ', 1, 1, RCMAX,
     $                          1, IDUMM, IDUMM, -1, -1, MYCOL )
                  IF( RCMIN.LE.ZERO ) THEN
                     INFO = -14
                  ELSE IF( N.GT.0 ) THEN
                     ROWCND = MAX( RCMIN, SMLNUM ) /
     $                        MIN( RCMAX, BIGNUM )
                  ELSE
                     ROWCND = ONE
                  END IF
               END IF
               IF( COLEQU .AND. INFO.EQ.0 ) THEN
                  RCMIN = BIGNUM
                  RCMAX = ZERO
                  DO 20 J = JJA, JJA+NQ-1
                     RCMIN = MIN( RCMIN, C( J ) )
                     RCMAX = MAX( RCMAX, C( J ) )
   20             CONTINUE
                  CALL SGAMN2D( ICTXT, 'Rowwise', ' ', 1, 1, RCMIN,
     $                          1, IDUMM, IDUMM, -1, -1, MYCOL )
                  CALL SGAMX2D( ICTXT, 'Rowwise', ' ', 1, 1, RCMAX,
     $                          1, IDUMM, IDUMM, -1, -1, MYCOL )
                  IF( RCMIN.LE.ZERO ) THEN
                     INFO = -15
                  ELSE IF( N.GT.0 ) THEN
                     COLCND = MAX( RCMIN, SMLNUM ) /
     $                        MIN( RCMAX, BIGNUM )
                  ELSE
                     COLCND = ONE
                  END IF
               END IF
            END IF
         END IF
*
         WORK( 1 ) = REAL( LWMIN )
         RWORK( 1 ) = REAL( LRWMIN )
         LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
         IF( INFO.EQ.0 ) THEN
            IF( IBROW.NE.IAROW ) THEN
               INFO = -18
            ELSE IF( IXROW.NE.IBROW ) THEN
               INFO = -22
            ELSE IF( DESCB( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(2000+NB_)
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -(2000+CTXT_)
            ELSE IF( DESCX( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(2400+NB_)
            ELSE IF( ICTXT.NE.DESCX( CTXT_ ) ) THEN
               INFO = -(2400+CTXT_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -29
            ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -31
            END IF
            IDUM1( 1 ) = ICHAR( FACT )
            IDUM2( 1 ) = 1
            IDUM1( 2 ) = ICHAR( TRANS )
            IDUM2( 2 ) = 2
            IF( LSAME( FACT, 'F' ) ) THEN
               IDUM1( 3 ) = ICHAR( EQUED )
               IDUM2( 3 ) = 14
               IF( LWORK.EQ.-1 ) THEN
                  IDUM1( 4 ) = -1
               ELSE
                  IDUM1( 4 ) = 1
               END IF
               IDUM2( 4 ) = 29
               IF( LRWORK.EQ.-1 ) THEN
                  IDUM1( 5 ) = -1
               ELSE
                  IDUM1( 5 ) = 1
               END IF
               IDUM2( 5 ) = 31
               CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 8, N, 3,
     $                        NRHS, 4, IB, JB, DESCB, 20, 5, IDUM1,
     $                        IDUM2, INFO )
            ELSE
               IF( LWORK.EQ.-1 ) THEN
                  IDUM1( 3 ) = -1
               ELSE
                  IDUM1( 3 ) = 1
               END IF
               IDUM2( 3 ) = 29
               IF( LRWORK.EQ.-1 ) THEN
                  IDUM1( 4 ) = -1
               ELSE
                  IDUM1( 4 ) = 1
               END IF
               IDUM2( 4 ) = 31
               CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 8, N, 3,
     $                        NRHS, 4, IB, JB, DESCB, 20, 4, IDUM1,
     $                        IDUM2, INFO )
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCGESVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      IF( EQUIL ) THEN
*
*        Compute row and column scalings to equilibrate the matrix A.
*
         CALL PCGEEQU( N, N, A, IA, JA, DESCA, R, C, ROWCND, COLCND,
     $                 AMAX, INFEQU )
         IF( INFEQU.EQ.0 ) THEN
*
*           Equilibrate the matrix.
*
            CALL PCLAQGE( N, N, A, IA, JA, DESCA, R, C, ROWCND, COLCND,
     $                    AMAX, EQUED )
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         END IF
      END IF
*
*     Scale the right-hand side.
*
      CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL, IIB,
     $              JJB, IBROW, IBCOL )
      NP = NUMROC( N+IROFFB, DESCB( MB_ ), MYROW, IBROW, NPROW )
      NRHSQ = NUMROC( NRHS+ICOFFB, DESCB( NB_ ), MYCOL, IBCOL, NPCOL )
      IF( MYROW.EQ.IBROW )
     $   NP = NP-IROFFB
      IF( MYCOL.EQ.IBCOL )
     $   NRHSQ = NRHSQ-ICOFFB
*
      IF( NOTRAN ) THEN
         IF( ROWEQU ) THEN
            DO 40 J = JJB, JJB+NRHSQ-1
               DO 30 I = IIB, IIB+NP-1
                  B( I+( J-1 )*DESCB( LLD_ ) ) = R( I )*
     $                                     B( I+( J-1 )*DESCB( LLD_ ) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( COLEQU ) THEN
*
*        Transpose the Column scale factors
*
         CALL DESCSET( CDESC, 1, N+ICOFFA, 1, DESCA( NB_ ), MYROW,
     $                 IACOL, ICTXT, 1 )
         CALL PSCOPY( N, C, 1, JA, CDESC, CDESC( LLD_ ), RWORK, IB, JB,
     $                DESCB, 1 )
         IF( MYCOL.EQ.IBCOL ) THEN
            CALL SGEBS2D( ICTXT, 'Rowwise', ' ', NP, 1, RWORK( IIB ),
     $                    DESCB( LLD_ ) )
         ELSE
            CALL SGEBR2D( ICTXT, 'Rowwise', ' ', NP, 1, RWORK( IIB ),
     $                    DESCB( LLD_ ), MYROW, IBCOL )
         END IF
         DO 60 J = JJB, JJB+NRHSQ-1
            DO 50 I = IIB, IIB+NP-1
               B( I+( J-1 )*DESCB( LLD_ ) ) = RWORK( I )*
     $                                     B( I+( J-1 )*DESCB( LLD_ ) )
   50       CONTINUE
   60    CONTINUE
      END IF
*
      IF( NOFACT.OR.EQUIL ) THEN
*
*        Compute the LU factorization of A.
*
         CALL PCLACPY( 'Full', N, N, A, IA, JA, DESCA, AF, IAF, JAF,
     $                 DESCAF )
         CALL PCGETRF( N, N, AF, IAF, JAF, DESCAF, IPIV, INFO )
*
*        Return if INFO is non-zero.
*
         IF( INFO.NE.0 ) THEN
            IF( INFO.GT.0 )
     $         RCOND = ZERO
            RETURN
         END IF
      END IF
*
*     Compute the norm of the matrix A.
*
      IF( NOTRAN ) THEN
         NORM = '1'
      ELSE
         NORM = 'I'
      END IF
      ANORM = PCLANGE( NORM, N, N, A, IA, JA, DESCA, RWORK )
*
*     Compute the reciprocal of the condition number of A.
*
      CALL PCGECON( NORM, N, AF, IAF, JAF, DESCAF, ANORM, RCOND, WORK,
     $              LWORK, RWORK, LRWORK, INFO )
*
*     Return if the matrix is singular to working precision.
*
      IF( RCOND.LT.PSLAMCH( ICTXT, 'Epsilon' ) ) THEN
         INFO = IA + N
         RETURN
      END IF
*
*     Compute the solution matrix X.
*
      CALL PCLACPY( 'Full', N, NRHS, B, IB, JB, DESCB, X, IX, JX,
     $              DESCX )
      CALL PCGETRS( TRANS, N, NRHS, AF, IAF, JAF, DESCAF, IPIV, X, IX,
     $              JX, DESCX, INFO )
*
*     Use iterative refinement to improve the computed solution and
*     compute error bounds and backward error estimates for it.
*
      CALL PCGERFS( TRANS, N, NRHS, A, IA, JA, DESCA, AF, IAF, JAF,
     $              DESCAF, IPIV, B, IB, JB, DESCB, X, IX, JX, DESCX,
     $              FERR, BERR, WORK, LWORK, RWORK, LRWORK, INFO )
*
*     Transform the solution matrix X to a solution of the original
*     system.
*
      CALL INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL, IIX,
     $              JJX, IXROW, IXCOL )
      NP = NUMROC( N+IROFFX, DESCX( MB_ ), MYROW, IXROW, NPROW )
      NRHSQ = NUMROC( NRHS+ICOFFX, DESCX( NB_ ), MYCOL, IXCOL, NPCOL )
      IF( MYROW.EQ.IBROW )
     $   NP = NP-IROFFX
      IF( MYCOL.EQ.IBCOL )
     $   NRHSQ = NRHSQ-ICOFFX
*
      IF( NOTRAN ) THEN
         IF( COLEQU ) THEN
*
*           Transpose the column scaling factors
*
            CALL DESCSET( CDESC, 1, N+ICOFFA, 1, DESCA( NB_ ), MYROW,
     $                    IACOL, ICTXT, 1 )
            CALL PSCOPY( N, C, 1, JA, CDESC, CDESC( LLD_ ), RWORK, IX,
     $                   JX, DESCX, 1 )
            IF( MYCOL.EQ.IBCOL ) THEN
                CALL SGEBS2D( ICTXT, 'Rowwise', ' ', NP, 1,
     $                        RWORK( IIX ), DESCX( LLD_ ) )
            ELSE
                CALL SGEBR2D( ICTXT, 'Rowwise', ' ', NP, 1,
     $                        RWORK( IIX ), DESCX( LLD_ ), MYROW,
     $                        IBCOL )
            END IF
*
            DO 80 J = JJX, JJX+NRHSQ-1
               DO 70 I = IIX, IIX+NP-1
                  X( I+( J-1 )*DESCX( LLD_ ) ) = RWORK( I )*
     $                                      X( I+( J-1 )*DESCX( LLD_ ) )
   70          CONTINUE
   80       CONTINUE
            DO 90 J = JJX, JJX+NRHSQ-1
               FERR( J ) = FERR( J ) / COLCND
   90       CONTINUE
         END IF
      ELSE IF( ROWEQU ) THEN
         DO 110 J = JJX, JJX+NRHSQ-1
            DO 100 I = IIX, IIX+NP-1
               X( I+( J-1 )*DESCX( LLD_ ) ) = R( I )*
     $                                     X( I+( J-1 )*DESCX( LLD_ ) )
  100       CONTINUE
  110    CONTINUE
         DO 120 J = JJX, JJX+NRHSQ-1
            FERR( J ) = FERR( J ) / ROWCND
  120    CONTINUE
      END IF
*
      WORK( 1 ) = REAL( LWMIN )
      RWORK( 1 ) = REAL( LRWMIN )
*
      RETURN
*
*     End of PCGESVX
*
      END
