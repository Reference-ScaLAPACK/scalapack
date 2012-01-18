      SUBROUTINE PSTRSEN( JOB, COMPQ, SELECT, PARA, N, T, IT, JT,
     $     DESCT, Q, IQ, JQ, DESCQ, WR, WI, M, S, SEP, WORK, LWORK,
     $     IWORK, LIWORK, INFO )
*
*     Contribution from the Department of Computing Science and HPC2N,
*     Umea University, Sweden
*
*  -- ScaLAPACK computational routine (version 2.0.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     Univ. of Colorado Denver and University of California, Berkeley.
*     January, 2012
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, JOB
      INTEGER            INFO, LIWORK, LWORK, M, N,
     $                   IT, JT, IQ, JQ
      REAL   S, SEP
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( N )
      INTEGER            PARA( 6 ), DESCT( * ), DESCQ( * ), IWORK( * )
      REAL               Q( * ), T( * ), WI( * ), WORK( * ), WR( * )
*     ..
*
*  Purpose
*  =======
*
*  PSTRSEN reorders the real Schur factorization of a real matrix
*  A = Q*T*Q**T, so that a selected cluster of eigenvalues appears
*  in the leading diagonal blocks of the upper quasi-triangular matrix
*  T, and the leading columns of Q form an orthonormal basis of the
*  corresponding right invariant subspace. The reordering is performed
*  by PSTRORD.
*
*  Optionally the routine computes the reciprocal condition numbers of
*  the cluster of eigenvalues and/or the invariant subspace. SCASY
*  library is needed for condition estimation.
*
*  T must be in Schur form (as returned by PSLAHQR), that is, block
*  upper triangular with 1-by-1 and 2-by-2 diagonal blocks.
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
*  JOB     (global input) CHARACTER*1
*          Specifies whether condition numbers are required for the
*          cluster of eigenvalues (S) or the invariant subspace (SEP):
*          = 'N': none;
*          = 'E': for eigenvalues only (S);
*          = 'V': for invariant subspace only (SEP);
*          = 'B': for both eigenvalues and invariant subspace (S and
*                 SEP).
*
*  COMPQ   (global input) CHARACTER*1
*          = 'V': update the matrix Q of Schur vectors;
*          = 'N': do not update Q.
*
*  SELECT  (global input) LOGICAL  array, dimension (N)
*          SELECT specifies the eigenvalues in the selected cluster. To
*          select a real eigenvalue w(j), SELECT(j) must be set to
*          .TRUE.. To select a complex conjugate pair of eigenvalues
*          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
*          either SELECT(j) or SELECT(j+1) or both must be set to
*          .TRUE.; a complex conjugate pair of eigenvalues must be
*          either both included in the cluster or both excluded.
*
*  PARA    (global input) INTEGER*6
*          Block parameters (some should be replaced by calls to
*          PILAENV and others by meaningful default values):
*          PARA(1) = maximum number of concurrent computational windows
*                    allowed in the algorithm;
*                    0 < PARA(1) <= min(NPROW,NPCOL) must hold;
*          PARA(2) = number of eigenvalues in each window;
*                    0 < PARA(2) < PARA(3) must hold;
*          PARA(3) = window size; PARA(2) < PARA(3) < DESCT(MB_)
*                    must hold;
*          PARA(4) = minimal percentage of flops required for
*                    performing matrix-matrix multiplications instead
*                    of pipelined orthogonal transformations;
*                    0 <= PARA(4) <= 100 must hold;
*          PARA(5) = width of block column slabs for row-wise
*                    application of pipelined orthogonal
*                    transformations in their factorized form;
*                    0 < PARA(5) <= DESCT(MB_) must hold.
*          PARA(6) = the maximum number of eigenvalues moved together
*                    over a process border; in practice, this will be
*                    approximately half of the cross border window size
*                    0 < PARA(6) <= PARA(2) must hold;
*
*  N       (global input) INTEGER
*          The order of the globally distributed matrix T. N >= 0.
*
*  T       (local input/output) REAL array,
*          dimension (LLD_T,LOCc(N)).
*          On entry, the local pieces of the global distributed
*          upper quasi-triangular matrix T, in Schur form. On exit, T is
*          overwritten by the local pieces of the reordered matrix T,
*          again in Schur form, with the selected eigenvalues in the
*          globally leading diagonal blocks.
*
*  IT      (global input) INTEGER
*  JT      (global input) INTEGER
*          The row and column index in the global array T indicating the
*          first column of sub( T ). IT = JT = 1 must hold.
*
*  DESCT   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the global distributed matrix T.
*
*  Q       (local input/output) REAL array,
*          dimension (LLD_Q,LOCc(N)).
*          On entry, if COMPQ = 'V', the local pieces of the global
*          distributed matrix Q of Schur vectors.
*          On exit, if COMPQ = 'V', Q has been postmultiplied by the
*          global orthogonal transformation matrix which reorders T; the
*          leading M columns of Q form an orthonormal basis for the
*          specified invariant subspace.
*          If COMPQ = 'N', Q is not referenced.
*
*  IQ      (global input) INTEGER
*  JQ      (global input) INTEGER
*          The column index in the global array Q indicating the
*          first column of sub( Q ). IQ = JQ = 1 must hold.
*
*  DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the global distributed matrix Q.
*
*  WR      (global output) REAL array, dimension (N)
*  WI      (global output) REAL array, dimension (N)
*          The real and imaginary parts, respectively, of the reordered
*          eigenvalues of T. The eigenvalues are in principle stored in
*          the same order as on the diagonal of T, with WR(i) = T(i,i)
*          and, if T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0
*          and WI(i+1) = -WI(i).
*          Note also that if a complex eigenvalue is sufficiently
*          ill-conditioned, then its value may differ significantly
*          from its value before reordering.
*
*  M       (global output) INTEGER
*          The dimension of the specified invariant subspace.
*          0 <= M <= N.
*
*  S       (global output) REAL
*          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
*          condition number for the selected cluster of eigenvalues.
*          S cannot underestimate the true reciprocal condition number
*          by more than a factor of sqrt(N). If M = 0 or N, S = 1.
*          If JOB = 'N' or 'V', S is not referenced.
*
*  SEP     (global output) REAL
*          If JOB = 'V' or 'B', SEP is the estimated reciprocal
*          condition number of the specified invariant subspace. If
*          M = 0 or N, SEP = norm(T).
*          If JOB = 'N' or 'E', SEP is not referenced.
*
*  WORK    (local workspace/output) REAL array,
*          dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (local input) INTEGER
*          The dimension of the array WORK.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by PXERBLA.
*
*  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
*
*  LIWORK  (local input) INTEGER
*          The dimension of the array IWORK.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value.
*          If the i-th argument is an array and the j-entry had
*          an illegal value, then INFO = -(i*1000+j), if the i-th
*          argument is a scalar and had an illegal value, then INFO = -i.
*          > 0: here we have several possibilites
*            *) Reordering of T failed because some eigenvalues are too
*               close to separate (the problem is very ill-conditioned);
*               T may have been partially reordered, and WR and WI
*               contain the eigenvalues in the same order as in T.
*               On exit, INFO = {the index of T where the swap failed}.
*            *) A 2-by-2 block to be reordered split into two 1-by-1
*               blocks and the second block failed to swap with an
*               adjacent block.
*               On exit, INFO = {the index of T where the swap failed}.
*            *) If INFO = N+1, there is no valid BLACS context (see the
*               BLACS documentation for details).
*            *) If INFO = N+2, the routines used in the calculation of
*               the condition numbers raised a positive warning flag
*               (see the documentation for PGESYCTD and PSYCTCON of the
*               SCASY library).
*            *) If INFO = N+3, PGESYCTD raised an input error flag;
*               please report this bug to the authors (see below).
*               If INFO = N+4, PSYCTCON raised an input error flag;
*               please report this bug to the authors (see below).
*          In a future release this subroutine may distinguish between
*          the case 1 and 2 above.
*
*  Method
*  ======
*
*  This routine performs parallel eigenvalue reordering in real Schur
*  form by parallelizing the approach proposed in [3]. The condition
*  number estimation part is performed by using techniques and code
*  from SCASY, see http://www.cs.umu.se/research/parallel/scasy.
*
*  Additional requirements
*  =======================
*
*  The following alignment requirements must hold:
*  (a) DESCT( MB_ ) = DESCT( NB_ ) = DESCQ( MB_ ) = DESCQ( NB_ )
*  (b) DESCT( RSRC_ ) = DESCQ( RSRC_ )
*  (c) DESCT( CSRC_ ) = DESCQ( CSRC_ )
*
*  All matrices must be blocked by a block factor larger than or
*  equal to two (3). This to simplify reordering across processor
*  borders in the presence of 2-by-2 blocks.
*
*  Limitations
*  ===========
*
*  This algorithm cannot work on submatrices of T and Q, i.e.,
*  IT = JT = IQ = JQ = 1 must hold. This is however no limitation
*  since PSLAHQR does not compute Schur forms of submatrices anyway.
*
*  References
*  ==========
*
*  [1] Z. Bai and J. W. Demmel; On swapping diagonal blocks in real
*      Schur form, Linear Algebra Appl., 186:73--95, 1993. Also as
*      LAPACK Working Note 54.
*
*  [2] Z. Bai, J. W. Demmel, and A. McKenney; On computing condition
*      numbers for the nonsymmetric eigenvalue problem, ACM Trans.
*      Math. Software, 19(2):202--223, 1993. Also as LAPACK Working
*      Note 13.
*
*  [3] D. Kressner; Block algorithms for reordering standard and
*      generalized Schur forms, ACM TOMS, 32(4):521-532, 2006.
*      Also LAPACK Working Note 171.
*
*  [4] R. Granat, B. Kagstrom, and D. Kressner; Parallel eigenvalue
*      reordering in real Schur form, Concurrency and Computations:
*      Practice and Experience, 21(9):1225-1250, 2009. Also as
*      LAPACK Working Note 192.
*
*  Parallel execution recommendations
*  ==================================
*
*  Use a square grid, if possible, for maximum performance. The block
*  parameters in PARA should be kept well below the data distribution
*  block size. In particular, see [3,4] for recommended settings for
*  these parameters.
*
*  In general, the parallel algorithm strives to perform as much work
*  as possible without crossing the block borders on the main block
*  diagonal.
*
*  Contributors
*  ============
*
*  Implemented by Robert Granat, Dept. of Computing Science and HPC2N,
*  Umea University, Sweden, March 2007,
*  in collaboration with Bo Kagstrom and Daniel Kressner.
*  Modified by Meiyue Shao, October 2011.
*
*  Revisions
*  =========
*
*  Please send bug-reports to granat@cs.umu.se
*
*  Keywords
*  ========
*
*  Real Schur form, eigenvalue reordering, Sylvester matrix equation
*
*  =====================================================================
*     ..
*     .. Parameters ..
      CHARACTER          TOP
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      REAL               ZERO, ONE
      PARAMETER          ( TOP = '1-Tree',
     $                     BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     ZERO = 0.0, ONE = 1.0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, WANTBH, WANTQ, WANTS, WANTSP
      INTEGER            ICOFFT12, ICTXT, IDUM1, IDUM2, IERR, ILOC1,
     $                   IPW1, ITER, ITT, JLOC1, JTT, K, LIWMIN, LLDT,
     $                   LLDQ, LWMIN, MMAX, MMIN, MYROW, MYCOL, N1, N2,
     $                   NB, NOEXSY, NPCOL, NPROCS, NPROW, SPACE,
     $                   T12ROWS, T12COLS, TCOLS, TCSRC, TROWS, TRSRC,
     $                   WRK1, IWRK1, WRK2, IWRK2, WRK3, IWRK3
      REAL               DPDUM1, ELEM, EST, SCALE, RNORM
*     .. Local Arrays ..
      INTEGER            DESCT12( DLEN_ ), MBNB2( 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      REAL               PSLANGE
      EXTERNAL           LSAME, NUMROC, PSLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCINIT,
     $                   IGAMX2D, INFOG2L, PSLACPY, PSTRORD, PXERBLA,
     $                   PCHK1MAT, PCHK2MAT
*     $                   , PGESYCTD, PSYCTCON
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCT( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW*NPCOL
*
*     Test if grid is O.K., i.e., the context is valid
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = N+1
      END IF
*
*     Check if workspace
*
      LQUERY = LWORK.EQ.-1 .OR. LIWORK.EQ.-1
*
*     Test dimensions for local sanity
*
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 5, N, 5, IT, JT, DESCT, 9, INFO )
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 5, N, 5, IQ, JQ, DESCQ, 13, INFO )
      END IF
*
*     Check the blocking sizes for alignment requirements
*
      IF( INFO.EQ.0 ) THEN
         IF( DESCT( MB_ ).NE.DESCT( NB_ ) ) INFO = -(1000*9 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCQ( MB_ ).NE.DESCQ( NB_ ) ) INFO = -(1000*13 + MB_)
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( DESCT( MB_ ).NE.DESCQ( MB_ ) ) INFO = -(1000*9 + MB_)
      END IF
*
*     Check the blocking sizes for minimum sizes
*
      IF( INFO.EQ.0 ) THEN
         IF( N.NE.DESCT( MB_ ) .AND. DESCT( MB_ ).LT.3 )
     $        INFO = -(1000*9 + MB_)
         IF( N.NE.DESCQ( MB_ ) .AND. DESCQ( MB_ ).LT.3 )
     $        INFO = -(1000*13 + MB_)
      END IF
*
*     Check parameters in PARA
*
      NB = DESCT( MB_ )
      IF( INFO.EQ.0 ) THEN
         IF( PARA(1).LT.1 .OR. PARA(1).GT.MIN(NPROW,NPCOL) )
     $        INFO = -(1000 * 4 + 1)
         IF( PARA(2).LT.1 .OR. PARA(2).GE.PARA(3) )
     $        INFO = -(1000 * 4 + 2)
         IF( PARA(3).LT.1 .OR. PARA(3).GT.NB )
     $        INFO = -(1000 * 4 + 3)
         IF( PARA(4).LT.0 .OR. PARA(4).GT.100 )
     $        INFO = -(1000 * 4 + 4)
         IF( PARA(5).LT.1 .OR. PARA(5).GT.NB )
     $        INFO = -(1000 * 4 + 5)
         IF( PARA(6).LT.1 .OR. PARA(6).GT.PARA(2) )
     $        INFO = -(1000 * 4 + 6)
      END IF
*
*     Check requirements on IT, JT, IQ and JQ
*
      IF( INFO.EQ.0 ) THEN
         IF( IT.NE.1 ) INFO = -7
         IF( JT.NE.IT ) INFO = -8
         IF( IQ.NE.1 ) INFO = -11
         IF( JQ.NE.IQ ) INFO = -12
      END IF
*
*     Test input parameters for global sanity
*
      IF( INFO.EQ.0 ) THEN
         CALL PCHK1MAT( N, 5, N, 5, IT, JT, DESCT, 9, 0, IDUM1,
     $        IDUM2, INFO )
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL PCHK1MAT( N, 5, N, 5, IQ, JQ, DESCQ, 13, 0, IDUM1,
     $        IDUM2, INFO )
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL PCHK2MAT( N, 5, N, 5, IT, JT, DESCT, 9, N, 5, N, 5,
     $        IQ, JQ, DESCQ, 13, 0, IDUM1, IDUM2, INFO )
      END IF
*
*     Decode and test the input parameters
*
      IF( INFO.EQ.0 .OR. LQUERY ) THEN
         WANTBH = LSAME( JOB, 'B' )
         WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
         WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
         WANTQ = LSAME( COMPQ, 'V' )
*
         IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.WANTS .AND. .NOT.WANTSP )
     $        THEN
            INFO = -1
         ELSEIF( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) THEN
            INFO = -2
         ELSEIF( N.LT.0 ) THEN
            INFO = -4
         ELSE
*
*           Extract local leading dimension
*
            LLDT = DESCT( LLD_ )
            LLDQ = DESCQ( LLD_ )
*
*           Check the SELECT vector for consistency and set M to the
*           dimension of the specified invariant subspace.
*
            M = 0
            DO 10 K = 1, N
*
*              IWORK(1:N) is an integer copy of SELECT.
*
               IF( SELECT(K) ) THEN
                  IWORK(K) = 1
               ELSE
                  IWORK(K) = 0
               END IF
               IF( K.LT.N ) THEN
                  CALL INFOG2L( K+1, K, DESCT, NPROW, NPCOL,
     $                 MYROW, MYCOL, ITT, JTT, TRSRC, TCSRC )
                  IF( MYROW.EQ.TRSRC .AND. MYCOL.EQ.TCSRC ) THEN
                     ELEM = T( (JTT-1)*LLDT + ITT )
                     IF( ELEM.NE.ZERO ) THEN
                        IF( SELECT(K) .AND. .NOT.SELECT(K+1) ) THEN
*                           INFO = -3
                           IWORK(K+1) = 1
                        ELSEIF( .NOT.SELECT(K) .AND. SELECT(K+1) ) THEN
*                           INFO = -3
                           IWORK(K) = 1
                        END IF
                     END IF
                  END IF
               END IF
               IF( SELECT(K) ) M = M + 1
 10         CONTINUE
            MMAX = M
            MMIN = M
            IF( NPROCS.GT.1 )
     $           CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, MMAX, 1, -1,
     $                -1, -1, -1, -1 )
            IF( NPROCS.GT.1 )
     $           CALL IGAMN2D( ICTXT, 'All', TOP, 1, 1, MMIN, 1, -1,
     $                -1, -1, -1, -1 )
            IF( MMAX.GT.MMIN ) THEN
               M = MMAX
               IF( NPROCS.GT.1 )
     $              CALL IGAMX2D( ICTXT, 'All', TOP, N, 1, IWORK, N,
     $                   -1, -1, -1, -1, -1 )
            END IF
*
*           Set parameters for deep pipelining in parallel
*           Sylvester solver.
*
            MBNB2( 1 ) = MIN( MAX( PARA( 3 ), PARA( 2 )*2 ), NB )
            MBNB2( 2 ) = MBNB2( 1 )
*
*           Compute needed workspace
*
            N1 = M
            N2 = N - M
            IF( WANTS ) THEN
c               CALL PGESYCTD( 'Solve', 'Schur', 'Schur', 'Notranspose',
c     $              'Notranspose', -1, 'Demand', N1, N2, T, 1, 1, DESCT,
c     $              T, N1+1, N1+1, DESCT, T, 1, N1+1, DESCT, MBNB2,
c     $              WORK, -1, IWORK(N+1), -1, NOEXSY, SCALE, IERR )
               WRK1 = INT(WORK(1))
               IWRK1 = IWORK(N+1)
               WRK1 = 0
               IWRK1 = 0
            ELSE
               WRK1 = 0
               IWRK1 = 0
            END IF
*
            IF( WANTSP ) THEN
c               CALL PSYCTCON( 'Notranspose', 'Notranspose', -1,
c     $              'Demand', N1, N2, T, 1, 1, DESCT, T, N1+1, N1+1,
c     $              DESCT, MBNB2, WORK, -1, IWORK(N+1), -1, EST, ITER,
c     $              IERR )
               WRK2 = INT(WORK(1))
               IWRK2 = IWORK(N+1)
               WRK2 = 0
               IWRK2 = 0
            ELSE
               WRK2 = 0
               IWRK2 = 0
            END IF
*
            TROWS = NUMROC( N, NB, MYROW, DESCT(RSRC_), NPROW )
            TCOLS = NUMROC( N, NB, MYCOL, DESCT(CSRC_), NPCOL )
            WRK3 = N + 7*NB**2 + 2*TROWS*PARA( 3 ) + TCOLS*PARA( 3 ) +
     $           MAX( TROWS*PARA( 3 ), TCOLS*PARA( 3 ) )
            IWRK3 = 5*PARA( 1 ) + PARA(2)*PARA(3) -
     $           PARA(2) * (PARA(2) + 1 ) / 2
*
            IF( WANTSP ) THEN
               LWMIN = MAX( 1, MAX( WRK2, WRK3) )
               LIWMIN = MAX( 1, MAX( IWRK2, IWRK3 ) )+N
            ELSE IF( LSAME( JOB, 'N' ) ) THEN
               LWMIN = MAX( 1, WRK3 )
               LIWMIN = IWRK3+N
            ELSE IF( LSAME( JOB, 'E' ) ) THEN
               LWMIN = MAX( 1, MAX( WRK1, WRK3) )
               LIWMIN = MAX( 1, MAX( IWRK1, IWRK3 ) )+N
            END IF
*
            IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -20
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -22
            END IF
         END IF
      END IF
*
*     Global maximum on info
*
      IF( NPROCS.GT.1 )
     $     CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, INFO, 1, -1, -1, -1,
     $          -1, -1 )
*
*     Return if some argument is incorrect
*
      IF( INFO.NE.0 .AND. .NOT.LQUERY ) THEN
         M = 0
         S = ONE
         SEP = ZERO
         CALL PXERBLA( ICTXT, 'PSTRSEN', -INFO )
         RETURN
      ELSEIF( LQUERY ) THEN
         WORK( 1 ) = FLOAT(LWMIN)
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( M.EQ.N .OR. M.EQ.0 ) THEN
         IF( WANTS )
     $        S = ONE
         IF( WANTSP )
     $        SEP = PSLANGE( '1', N, N, T, IT, JT, DESCT, WORK )
         GO TO 50
      END IF
*
*     Reorder the eigenvalues.
*
      CALL PSTRORD( COMPQ, IWORK, PARA, N, T, IT, JT,
     $     DESCT, Q, IQ, JQ, DESCQ, WR, WI, M, WORK, LWORK,
     $     IWORK(N+1), LIWORK-N, INFO )
*
      IF( WANTS ) THEN
*
*        Solve Sylvester equation T11*R - R*T2 = scale*T12 for R in
*        parallel.
*
*        Copy T12 to workspace.
*
         CALL INFOG2L( 1, N1+1, DESCT, NPROW, NPCOL, MYROW,
     $        MYCOL, ILOC1, JLOC1, TRSRC, TCSRC )
         ICOFFT12 = MOD( N1, NB )
         T12ROWS = NUMROC( N1, NB, MYROW, TRSRC, NPROW )
         T12COLS = NUMROC( N2+ICOFFT12, NB, MYCOL, TCSRC, NPCOL )
         CALL DESCINIT( DESCT12, N1, N2+ICOFFT12, NB, NB, TRSRC,
     $        TCSRC, ICTXT, MAX(1,T12ROWS), IERR )
         CALL PSLACPY( 'All', N1, N2, T, 1, N1+1, DESCT, WORK,
     $        1, 1+ICOFFT12, DESCT12 )
*
*        Solve the equation to get the solution in workspace.
*
         SPACE = DESCT12( LLD_ ) * T12COLS
         IPW1 = 1 + SPACE
c         CALL PGESYCTD( 'Solve', 'Schur', 'Schur', 'Notranspose',
c     $        'Notranspose', -1, 'Demand', N1, N2, T, 1, 1, DESCT, T,
c     $        N1+1, N1+1, DESCT, WORK, 1, 1+ICOFFT12, DESCT12, MBNB2,
c     $        WORK(IPW1), LWORK-SPACE+1, IWORK(N+1), LIWORK-N, NOEXSY,
c     $        SCALE, IERR )
         IF( IERR.LT.0 ) THEN
            INFO = N+3
         ELSE
            INFO = N+2
         END IF
*
*        Estimate the reciprocal of the condition number of the cluster
*        of eigenvalues.
*
         RNORM = PSLANGE( 'Frobenius', N1, N2, WORK, 1, 1+ICOFFT12,
     $        DESCT12, DPDUM1 )
         IF( RNORM.EQ.ZERO ) THEN
            S = ONE
         ELSE
            S = SCALE / ( SQRT( SCALE*SCALE / RNORM+RNORM )*
     $           SQRT( RNORM ) )
         END IF
      END IF
*
      IF( WANTSP ) THEN
*
*        Estimate sep(T11,T21) in parallel.
*
c         CALL PSYCTCON( 'Notranspose', 'Notranspose', -1, 'Demand', N1,
c     $        N2, T, 1, 1, DESCT, T, N1+1, N1+1, DESCT, MBNB2, WORK,
c     $        LWORK, IWORK(N+1), LIWORK-N, EST, ITER, IERR )
         EST = EST * SQRT(FLOAT(N1*N2))
         SEP = ONE / EST
         IF( IERR.LT.0 ) THEN
            INFO = N+4
         ELSE
            INFO = N+2
         END IF
      END IF
*
*     Return to calling program.
*
 50   CONTINUE
*
      RETURN
*
*     End of PSTRSEN
*
      END
*
