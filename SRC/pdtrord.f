      SUBROUTINE PDTRORD( COMPQ, SELECT, PARA, N, T, IT, JT,
     $     DESCT, Q, IQ, JQ, DESCQ, WR, WI, M, WORK, LWORK,
     $     IWORK, LIWORK, INFO )
*
*     Contribution from the Department of Computing Science and HPC2N,
*     Umea University, Sweden
*
*  -- ScaLAPACK routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            INFO, LIWORK, LWORK, M, N,
     $                   IT, JT, IQ, JQ
*     ..
*     .. Array Arguments ..
      INTEGER            SELECT( * )
      INTEGER            PARA( 6 ), DESCT( * ), DESCQ( * ), IWORK( * )
      DOUBLE PRECISION   Q( * ), T( * ), WI( * ), WORK( * ), WR( * )
*     ..
*
*  Purpose
*  =======
*
*  PDTRORD reorders the real Schur factorization of a real matrix
*  A = Q*T*Q**T, so that a selected cluster of eigenvalues appears
*  in the leading diagonal blocks of the upper quasi-triangular matrix
*  T, and the leading columns of Q form an orthonormal basis of the
*  corresponding right invariant subspace.
*
*  T must be in Schur form (as returned by PDLAHQR), that is, block
*  upper triangular with 1-by-1 and 2-by-2 diagonal blocks.
*
*  This subroutine uses a delay and accumulate procedure for performing
*  the off-diagonal updates (see references for details).
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
*
*  COMPQ   (global input) CHARACTER*1
*          = 'V': update the matrix Q of Schur vectors;
*          = 'N': do not update Q.
*
*  SELECT  (global input/output) INTEGER array, dimension (N)
*          SELECT specifies the eigenvalues in the selected cluster. To
*          select a real eigenvalue w(j), SELECT(j) must be set to 1.
*          To select a complex conjugate pair of eigenvalues
*          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
*          either SELECT(j) or SELECT(j+1) or both must be set to 1;
*          a complex conjugate pair of eigenvalues must be
*          either both included in the cluster or both excluded.
*          On output, the (partial) reordering is displayed.
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
*  T       (local input/output) DOUBLE PRECISION array,
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
*  Q       (local input/output) DOUBLE PRECISION array,
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
*  WR      (global output) DOUBLE PRECISION array, dimension (N)
*  WI      (global output) DOUBLE PRECISION array, dimension (N)
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
*  WORK    (local workspace/output) DOUBLE PRECISION array,
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
*          In a future release this subroutine may distinguish between
*          the case 1 and 2 above.
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
*  equal to two (3). This is to simplify reordering across processor
*  borders in the presence of 2-by-2 blocks.
*
*  Limitations
*  ===========
*
*  This algorithm cannot work on submatrices of T and Q, i.e.,
*  IT = JT = IQ = JQ = 1 must hold. This is however no limitation
*  since PDLAHQR does not compute Schur forms of submatrices anyway.
*
*  References
*  ==========
*
*  [1] Z. Bai and J. W. Demmel; On swapping diagonal blocks in real
*      Schur form, Linear Algebra Appl., 186:73--95, 1993. Also as
*      LAPACK Working Note 54.
*
*  [2] D. Kressner; Block algorithms for reordering standard and
*      generalized Schur forms, ACM TOMS, 32(4):521-532, 2006.
*      Also LAPACK Working Note 171.
*
*  [3] R. Granat, B. Kagstrom, and D. Kressner; Parallel eigenvalue
*      reordering in real Schur form, Concurrency and Computations:
*      Practice and Experience, 21(9):1225-1250, 2009. Also as
*      LAPACK Working Note 192.
*
*  Parallel execution recommendations
*  ==================================
*
*  Use a square grid, if possible, for maximum performance. The block
*  parameters in PARA should be kept well below the data distribution
*  block size. In particular, see [3] for recommended settings for
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
*  Real Schur form, eigenvalue reordering
*
*  =====================================================================
*     ..
*     .. Parameters ..
      CHARACTER          TOP
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( TOP = '1-Tree',
     $                     BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, PAIR, SWAP, WANTQ,
     $                   ISHH, FIRST, SKIP1CR, BORDER, LASTWAIT
      INTEGER            NPROW, NPCOL, MYROW, MYCOL, NB, NPROCS,
     $                   IERR, DIM1, INDX, LLDT, TRSRC, TCSRC, ILOC1,
     $                   JLOC1, MYIERR, ICTXT,
     $                   RSRC1, CSRC1, ILOC3, JLOC3, TRSRC3,
     $                   TCSRC3, ILOC, JLOC, TRSRC4, TCSRC4,
     $                   FLOPS, I, ILO, IHI, J, K, KK, KKS,
     $                   KS, LIWMIN, LWMIN, MMULT, N1, N2,
     $                   NCB, NDTRAF, NITRAF, NWIN, NUMWIN, PDTRAF,
     $                   PITRAF, PDW, WINEIG, WINSIZ, LLDQ,
     $                   RSRC, CSRC, ILILO, ILIHI, ILSEL, IRSRC,
     $                   ICSRC, IPIW, IPW1, IPW2, IPW3, TIHI, TILO,
     $                   LIHI, WINDOW, LILO, LSEL, BUFFER,
     $                   NMWIN2, BUFFLEN, LROWS, LCOLS, ILOC2, JLOC2,
     $                   WNEICR, WINDOW0, RSRC4, CSRC4, LIHI4, RSRC3,
     $                   CSRC3, RSRC2, CSRC2, LIHIC, LIHI1, ILEN4,
     $                   SELI4, ILEN1, DIM4, IPW4, QROWS, TROWS,
     $                   TCOLS, IPW5, IPW6, IPW7, IPW8, JLOC4,
     $                   EAST, WEST, ILOC4, SOUTH, NORTH, INDXS,
     $                   ITT, JTT, ILEN, DLEN, INDXE, TRSRC1, TCSRC1,
     $                   TRSRC2, TCSRC2, ILOS, DIR, TLIHI, TLILO, TLSEL,
     $                   ROUND, LAST, WIN0S, WIN0E, WINE, MMAX, MMIN
      DOUBLE PRECISION   ELEM, ELEM1, ELEM2, ELEM3, ELEM4, SN, CS, TMP,
     $                   ELEM5
*     ..
*     .. Local Arrays ..
      INTEGER            IBUFF( 8 ), IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC, INDXG2P, INDXG2L
      EXTERNAL           LSAME, NUMROC, INDXG2P, INDXG2L
*     ..
*     .. External Subroutines ..
      EXTERNAL           PDLACPY, PXERBLA, PCHK1MAT, PCHK2MAT,
     $                   DGEMM, DLAMOV, ILACPY, CHK1MAT,
     $                   INFOG2L, DGSUM2D, DGESD2D, DGERV2D, DGEBS2D,
     $                   DGEBR2D, IGSUM2D, BLACS_GRIDINFO, IGEBS2D,
     $                   IGEBR2D, IGAMX2D, IGAMN2D, BDLAAPP, BDTREXC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT, MIN
*     ..
*     .. Local Functions ..
      INTEGER            ICEIL
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters.
*
      ICTXT = DESCT( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW*NPCOL
*
*     Test if grid is O.K., i.e., the context is valid.
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = N+1
      END IF
*
*     Check if workspace query.
*
      LQUERY = LWORK.EQ.-1 .OR. LIWORK.EQ.-1
*
*     Test dimensions for local sanity.
*
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 5, N, 5, IT, JT, DESCT, 9, INFO )
      END IF
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 5, N, 5, IQ, JQ, DESCQ, 13, INFO )
      END IF
*
*     Check the blocking sizes for alignment requirements.
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
*     Check the blocking sizes for minimum sizes.
*
      IF( INFO.EQ.0 ) THEN
         IF( N.NE.DESCT( MB_ ) .AND. DESCT( MB_ ).LT.3 )
     $      INFO = -(1000*9 + MB_)
         IF( N.NE.DESCQ( MB_ ) .AND. DESCQ( MB_ ).LT.3 )
     $      INFO = -(1000*13 + MB_)
      END IF
*
*     Check parameters in PARA.
*
      NB = DESCT( MB_ )
      IF( INFO.EQ.0 ) THEN
         IF( PARA(1).LT.1 .OR. PARA(1).GT.MIN(NPROW,NPCOL) )
     $      INFO = -(1000 * 4 + 1)
         IF( PARA(2).LT.1 .OR. PARA(2).GE.PARA(3) )
     $      INFO = -(1000 * 4 + 2)
         IF( PARA(3).LT.1 .OR. PARA(3).GT.NB )
     $      INFO = -(1000 * 4 + 3)
         IF( PARA(4).LT.0 .OR. PARA(4).GT.100 )
     $      INFO = -(1000 * 4 + 4)
         IF( PARA(5).LT.1 .OR. PARA(5).GT.NB )
     $      INFO = -(1000 * 4 + 5)
         IF( PARA(6).LT.1 .OR. PARA(6).GT.PARA(2) )
     $      INFO = -(1000 * 4 + 6)
      END IF
*
*     Check requirements on IT, JT, IQ and JQ.
*
      IF( INFO.EQ.0 ) THEN
         IF( IT.NE.1 ) INFO = -6
         IF( JT.NE.IT ) INFO = -7
         IF( IQ.NE.1 ) INFO = -10
         IF( JQ.NE.IQ ) INFO = -11
      END IF
*
*     Test input parameters for global sanity.
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
*     Decode and test the input parameters.
*
      IF( INFO.EQ.0 .OR. LQUERY ) THEN
*
         WANTQ = LSAME( COMPQ, 'V' )
         IF( N.LT.0 ) THEN
            INFO = -4
         ELSE
*
*           Extract local leading dimension.
*
            LLDT = DESCT( LLD_ )
            LLDQ = DESCQ( LLD_ )
*
*           Check the SELECT vector for consistency and set M to the
*           dimension of the specified invariant subspace.
*
            M = 0
            DO 10 K = 1, N
               IF( K.LT.N ) THEN
                  CALL INFOG2L( K+1, K, DESCT, NPROW, NPCOL,
     $                 MYROW, MYCOL, ITT, JTT, TRSRC, TCSRC )
                  IF( MYROW.EQ.TRSRC .AND. MYCOL.EQ.TCSRC ) THEN
                     ELEM = T( (JTT-1)*LLDT + ITT )
                     IF( ELEM.NE.ZERO ) THEN
                        IF( SELECT(K).NE.0 .AND.
     $                       SELECT(K+1).EQ.0 ) THEN
*                           INFO = -2
                           SELECT(K+1) = 1
                        ELSEIF( SELECT(K).EQ.0 .AND.
     $                          SELECT(K+1).NE.0 ) THEN
*                           INFO = -2
                           SELECT(K) = 1
                        END IF
                     END IF
                  END IF
               END IF
               IF( SELECT(K).NE.0 ) M = M + 1
 10         CONTINUE
            MMAX = M
            MMIN = M
            IF( NPROCS.GT.1 )
     $         CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, MMAX, 1, -1,
     $              -1, -1, -1, -1 )
            IF( NPROCS.GT.1 )
     $         CALL IGAMN2D( ICTXT, 'All', TOP, 1, 1, MMIN, 1, -1,
     $              -1, -1, -1, -1 )
            IF( MMAX.GT.MMIN ) THEN
               M = MMAX
               IF( NPROCS.GT.1 )
     $            CALL IGAMX2D( ICTXT, 'All', TOP, N, 1, SELECT, N,
     $                 -1, -1, -1, -1, -1 )
            END IF
*
*           Compute needed workspace.
*
            N1 = M
            N2 = N - M
*
            TROWS = NUMROC( N, NB, MYROW, DESCT(RSRC_), NPROW )
            TCOLS = NUMROC( N, NB, MYCOL, DESCT(CSRC_), NPCOL )
            LWMIN = N + 7*NB**2 + 2*TROWS*PARA( 3 ) + TCOLS*PARA( 3 ) +
     $           MAX( TROWS*PARA( 3 ), TCOLS*PARA( 3 ) )
            LIWMIN = 5*PARA( 1 ) + PARA( 2 )*PARA( 3 ) -
     $           PARA( 2 ) * ( PARA( 2 ) + 1 ) / 2
*
            IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -17
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -19
            END IF
         END IF
      END IF
*
*     Global maximum on info.
*
      IF( NPROCS.GT.1 )
     $   CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, INFO, 1, -1, -1, -1,
     $        -1, -1 )
*
*     Return if some argument is incorrect.
*
      IF( INFO.NE.0 .AND. .NOT.LQUERY ) THEN
         M = 0
         CALL PXERBLA( ICTXT, 'PDTRORD', -INFO )
         RETURN
      ELSEIF( LQUERY ) THEN
         WORK( 1 ) = DBLE(LWMIN)
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( M.EQ.N .OR. M.EQ.0 ) GO TO 545
*
*     Set parameters.
*
      NUMWIN = PARA( 1 )
      WINEIG = MAX( PARA( 2 ), 2 )
      WINSIZ = MIN( MAX( PARA( 3 ), PARA( 2 )*2 ), NB )
      MMULT  = PARA( 4 )
      NCB    = PARA( 5 )
      WNEICR = PARA( 6 )
*
*     Insert some pointers into INTEGER workspace.
*
*     Information about all the active windows is stored
*     in IWORK( 1:5*NUMWIN ). Each processor has a copy.
*       LILO: start position
*       LIHI: stop position
*       LSEL: number of selected eigenvalues
*       RSRC: processor id (row)
*       CSRC: processor id (col)
*     IWORK( IPIW+ ) contain information of orthogonal transformations.
*
      ILILO = 1
      ILIHI = ILILO + NUMWIN
      ILSEL = ILIHI + NUMWIN
      IRSRC = ILSEL + NUMWIN
      ICSRC = IRSRC + NUMWIN
      IPIW  = ICSRC + NUMWIN
*
*     Insert some pointers into DOUBLE PRECISION workspace - for now we
*     only need two pointers.
*
      IPW1 = 1
      IPW2 = IPW1 + NB
*
*     Collect the selected blocks at the top-left corner of T.
*
*     Globally: ignore eigenvalues that are already in order.
*     ILO is a global variable and is kept updated to be consistent
*     throughout the process mesh.
*
      ILO = 0
 40   CONTINUE
      ILO = ILO + 1
      IF( ILO.LE.N ) THEN
         IF( SELECT(ILO).NE.0 ) GO TO 40
      END IF
*
*     Globally: start the collection at the top of the matrix. Here,
*     IHI is a global variable and is kept updated to be consistent
*     throughout the process mesh.
*
      IHI = N
*
*     Globally:  While ( ILO <= M ) do
 50   CONTINUE
*
      IF( ILO.LE.M ) THEN
*
*        Depending on the value of ILO, find the diagonal block index J,
*        such that T(1+(J-1)*NB:1+J*NB,1+(J-1)*NB:1+J*NB) contains the
*        first unsorted eigenvalue. Check that J does not point to a
*        block with only one selected eigenvalue in the last position
*        which belongs to a splitted 2-by-2 block.
*
         ILOS = ILO - 1
 52      CONTINUE
         ILOS = ILOS + 1
         IF( SELECT(ILOS).EQ.0 ) GO TO 52
         IF( ILOS.LT.N ) THEN
            IF( SELECT(ILOS+1).NE.0 .AND. MOD(ILOS,NB).EQ.0 ) THEN
               CALL PDELGET( 'All', TOP, ELEM, T, ILOS+1, ILOS, DESCT )
               IF( ELEM.NE.ZERO ) GO TO 52
            END IF
         END IF
         J = ICEIL(ILOS,NB)
*
*        Globally: Set start values of LILO and LIHI for all processes.
*        Choose also the number of selected eigenvalues at top of each
*        diagonal block such that the number of eigenvalues which remain
*        to be reordered is an integer multiple of WINEIG.
*
*        All the information is saved into the INTEGER workspace such
*        that all processors are aware of each others operations.
*
*        Compute the number of concurrent windows.
*
         NMWIN2 = (ICEIL(IHI,NB)*NB - (ILO-MOD(ILO,NB)+1)+1) / NB
         NMWIN2 = MIN( MIN( NUMWIN, NMWIN2 ), ICEIL(N,NB) - J + 1 )
*
*        For all windows, set LSEL = 0 and find a proper start value of
*        LILO such that LILO points at the first non-selected entry in
*        the corresponding diagonal block of T.
*
         DO 80 K = 1, NMWIN2
            IWORK( ILSEL+K-1) = 0
            IWORK( ILILO+K-1) = MAX( ILO, (J-1)*NB+(K-1)*NB+1 )
            LILO = IWORK( ILILO+K-1 )
 82         CONTINUE
            IF( SELECT(LILO).NE.0 .AND. LILO.LT.(J+K-1)*NB ) THEN
               LILO = LILO + 1
               IF( LILO.LE.N ) GO TO 82
            END IF
            IWORK( ILILO+K-1 ) = LILO
*
*           Fix each LILO to ensure that no 2-by-2 block is cut in top
*           of the submatrix (LILO:LIHI,LILO:LIHI).
*
            LILO = IWORK(ILILO+K-1)
            IF( LILO.GT.NB ) THEN
               CALL PDELGET( 'All', TOP, ELEM, T, LILO, LILO-1, DESCT )
               IF( ELEM.NE.ZERO ) THEN
                  IF( LILO.LT.(J+K-1)*NB ) THEN
                     IWORK(ILILO+K-1) = IWORK(ILILO+K-1) + 1
                  ELSE
                     IWORK(ILILO+K-1) = IWORK(ILILO+K-1) - 1
                  END IF
               END IF
            END IF
*
*           Set a proper LIHI value for each window. Also find the
*           processors corresponding to the corresponding windows.
*
            IWORK( ILIHI+K-1 ) =  IWORK( ILILO+K-1 )
            IWORK( IRSRC+K-1 ) = INDXG2P( IWORK(ILILO+K-1), NB, MYROW,
     $           DESCT( RSRC_ ), NPROW )
            IWORK( ICSRC+K-1 ) = INDXG2P( IWORK(ILILO+K-1), NB, MYCOL,
     $           DESCT( CSRC_ ), NPCOL )
            TILO = IWORK(ILILO+K-1)
            TIHI = MIN( N, ICEIL( TILO, NB ) * NB )
            DO 90 KK = TIHI, TILO, -1
               IF( SELECT(KK).NE.0 ) THEN
                  IWORK(ILIHI+K-1) = MAX(IWORK(ILIHI+K-1) , KK )
                  IWORK(ILSEL+K-1) = IWORK(ILSEL+K-1) + 1
                  IF( IWORK(ILSEL+K-1).GT.WINEIG ) THEN
                     IWORK(ILIHI+K-1) = KK
                     IWORK(ILSEL+K-1) = 1
                  END IF
               END IF
 90         CONTINUE
*
*           Fix each LIHI to avoid that bottom of window cuts 2-by-2
*           block. We exclude such a block if located on block (process)
*           border and on window border or if an inclusion would cause
*           violation on the maximum number of eigenvalues to reorder
*           inside each window. If only on window border, we include it.
*           The excluded block is included automatically later when a
*           subcluster is reordered into the block from South-East.
*
            LIHI = IWORK(ILIHI+K-1)
            IF( LIHI.LT.N ) THEN
               CALL PDELGET( 'All', TOP, ELEM, T, LIHI+1, LIHI, DESCT )
               IF( ELEM.NE.ZERO ) THEN
                  IF( ICEIL( LIHI, NB ) .NE. ICEIL( LIHI+1, NB ) .OR.
     $                 IWORK( ILSEL+K-1 ).EQ.WINEIG ) THEN
                     IWORK( ILIHI+K-1 ) = IWORK( ILIHI+K-1 ) - 1
                     IF( IWORK( ILSEL+K-1 ).GT.2 )
     $                  IWORK( ILSEL+K-1 ) = IWORK( ILSEL+K-1 ) - 1
                  ELSE
                     IWORK( ILIHI+K-1 ) = IWORK( ILIHI+K-1 ) + 1
                     IF( SELECT(LIHI+1).NE.0 )
     $                  IWORK( ILSEL+K-1 ) = IWORK( ILSEL+K-1 ) + 1
                  END IF
               END IF
            END IF
 80      CONTINUE
*
*        Fix the special cases of LSEL = 0 and LILO = LIHI for each
*        window by assuring that the stop-condition for local reordering
*        is fulfilled directly. Do this by setting LIHI = startposition
*        for the corresponding block and LILO = LIHI + 1.
*
         DO 85 K = 1, NMWIN2
            LILO = IWORK( ILILO + K - 1 )
            LIHI = IWORK( ILIHI + K - 1 )
            LSEL = IWORK( ILSEL + K - 1 )
            IF( LSEL.EQ.0 .OR. LILO.EQ.LIHI ) THEN
               LIHI = IWORK( ILIHI + K - 1 )
               IWORK( ILIHI + K - 1 ) = (ICEIL(LIHI,NB)-1)*NB + 1
               IWORK( ILILO + K - 1 ) = IWORK( ILIHI + K - 1 ) + 1
            END IF
 85      CONTINUE
*
*        Associate all processors with the first computational window
*        that should be activated, if possible.
*
         LILO = IHI
         LIHI = ILO
         LSEL = M
         FIRST = .TRUE.
         DO 95 WINDOW = 1, NMWIN2
            RSRC = IWORK(IRSRC+WINDOW-1)
            CSRC = IWORK(ICSRC+WINDOW-1)
            IF( MYROW.EQ.RSRC .OR. MYCOL.EQ.CSRC ) THEN
               TLILO = IWORK( ILILO + WINDOW - 1 )
               TLIHI = IWORK( ILIHI + WINDOW - 1 )
               TLSEL = IWORK( ILSEL + WINDOW - 1 )
               IF( (.NOT. ( LIHI .GE. LILO + LSEL ) ) .AND.
     $              ( (TLIHI .GE. TLILO + TLSEL) .OR. FIRST ) ) THEN
                  IF( FIRST ) FIRST = .FALSE.
                  LILO = TLILO
                  LIHI = TLIHI
                  LSEL = TLSEL
                  GO TO 97
               END IF
            END IF
 95      CONTINUE
 97      CONTINUE
*
*        Exclude all processors that are not involved in any
*        computational window right now.
*
         IERR = 0
         IF( LILO.EQ.IHI .AND. LIHI.EQ.ILO .AND. LSEL.EQ.M )
     $      GO TO 114
*
*        Make sure all processors associated with a compuational window
*        enter the local reordering the first time.
*
         FIRST = .TRUE.
*
*        Globally for all computational windows:
*        While ( LIHI >= LILO + LSEL ) do
         ROUND = 1
 130     CONTINUE
         IF( FIRST .OR. ( LIHI .GE. LILO + LSEL ) ) THEN
*
*           Perform computations in parallel: loop through all
*           compuational windows, do local reordering and accumulate
*           transformations, broadcast them in the corresponding block
*           row and columns and compute the corresponding updates.
*
            DO 110 WINDOW = 1, NMWIN2
               RSRC = IWORK(IRSRC+WINDOW-1)
               CSRC = IWORK(ICSRC+WINDOW-1)
*
*              The process on the block diagonal computes the
*              reordering.
*
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  LILO = IWORK(ILILO+WINDOW-1)
                  LIHI = IWORK(ILIHI+WINDOW-1)
                  LSEL = IWORK(ILSEL+WINDOW-1)
*
*                 Compute the local value of I -- start position.
*
                  I = MAX( LILO, LIHI - WINSIZ + 1 )
*
*                 Fix my I to avoid that top of window cuts a 2-by-2
*                 block.
*
                  IF( I.GT.LILO ) THEN
                     CALL INFOG2L( I, I-1, DESCT, NPROW, NPCOL, MYROW,
     $                    MYCOL, ILOC, JLOC, RSRC, CSRC )
                     IF( T( LLDT*(JLOC-1) + ILOC ).NE.ZERO )
     $                  I = I + 1
                  END IF
*
*                 Compute local indicies for submatrix to operate on.
*
                  CALL INFOG2L( I, I, DESCT, NPROW, NPCOL,
     $                 MYROW, MYCOL, ILOC1, JLOC1, RSRC, CSRC )
*
*                 The active window is ( I:LIHI, I:LIHI ). Reorder
*                 eigenvalues within this window and pipeline
*                 transformations.
*
                  NWIN = LIHI - I + 1
                  KS = 0
                  PITRAF = IPIW
                  PDTRAF = IPW2
*
                  PAIR = .FALSE.
                  DO 140 K = I, LIHI
                     IF( PAIR ) THEN
                        PAIR = .FALSE.
                     ELSE
                        SWAP = SELECT( K ).NE.0
                        IF( K.LT.LIHI ) THEN
                           CALL INFOG2L( K+1, K, DESCT, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC, CSRC )
                           IF( T( LLDT*(JLOC-1) + ILOC ).NE.ZERO )
     $                        PAIR = .TRUE.
                        END IF
                        IF( SWAP ) THEN
                           KS = KS + 1
*
*                       Swap the K-th block to position I+KS-1.
*
                           IERR = 0
                           KK  = K - I + 1
                           KKS = KS
                           IF( KK.NE.KS ) THEN
                              NITRAF = LIWORK - PITRAF + 1
                              NDTRAF = LWORK - PDTRAF + 1
                              CALL BDTREXC( NWIN,
     $                             T(LLDT*(JLOC1-1) + ILOC1), LLDT, KK,
     $                             KKS, NITRAF, IWORK( PITRAF ), NDTRAF,
     $                             WORK( PDTRAF ), WORK(IPW1), IERR )
                              PITRAF = PITRAF + NITRAF
                              PDTRAF = PDTRAF + NDTRAF
*
*                             Update array SELECT.
*
                              IF ( PAIR ) THEN
                                 DO 150 J = I+KK-1, I+KKS, -1
                                    SELECT(J+1) = SELECT(J-1)
 150                             CONTINUE
                                 SELECT(I+KKS-1) = 1
                                 SELECT(I+KKS) = 1
                              ELSE
                                 DO 160 J = I+KK-1, I+KKS, -1
                                    SELECT(J) = SELECT(J-1)
 160                             CONTINUE
                                 SELECT(I+KKS-1) = 1
                              END IF
*
                              IF ( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
*
*                                Some blocks are too close to swap:
*                                prepare to leave in a clean fashion. If
*                                IERR.EQ.2, we must update SELECT to
*                                account for the fact that the 2 by 2
*                                block to be reordered did split and the
*                                first part of this block is already
*                                reordered.
*
                                 IF ( IERR.EQ.2 ) THEN
                                    SELECT( I+KKS-3 ) = 1
                                    SELECT( I+KKS-1 ) = 0
                                    KKS = KKS + 1
                                 END IF
*
*                                Update off-diagonal blocks immediately.
*
                                 GO TO 170
                              END IF
                              KS = KKS
                           END IF
                           IF( PAIR )
     $                        KS = KS + 1
                        END IF
                     END IF
 140              CONTINUE
               END IF
 110        CONTINUE
 170        CONTINUE
*
*           The on-diagonal processes save their information from the
*           local reordering in the integer buffer. This buffer is
*           broadcasted to updating processors, see below.
*
            DO 175 WINDOW = 1, NMWIN2
               RSRC = IWORK(IRSRC+WINDOW-1)
               CSRC = IWORK(ICSRC+WINDOW-1)
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  IBUFF( 1 ) = I
                  IBUFF( 2 ) = NWIN
                  IBUFF( 3 ) = PITRAF
                  IBUFF( 4 ) = KS
                  IBUFF( 5 ) = PDTRAF
                  IBUFF( 6 ) = NDTRAF
                  ILEN = PITRAF - IPIW
                  DLEN = PDTRAF - IPW2
                  IBUFF( 7 ) = ILEN
                  IBUFF( 8 ) = DLEN
               END IF
 175        CONTINUE
*
*           For the updates with respect to the local reordering, we
*           organize the updates in two phases where the update
*           "direction" (controlled by the DIR variable below) is first
*           chosen to be the corresponding rows, then the corresponding
*           columns.
*
            DO 1111 DIR = 1, 2
*
*           Broadcast information about the reordering and the
*           accumulated transformations: I, NWIN, PITRAF, NITRAF,
*           PDTRAF, NDTRAF. If no broadcast is performed, use an
*           artificial value of KS to prevent updating indicies for
*           windows already finished (use KS = -1).
*
            DO 111 WINDOW = 1, NMWIN2
               RSRC = IWORK(IRSRC+WINDOW-1)
               CSRC = IWORK(ICSRC+WINDOW-1)
               IF( MYROW.EQ.RSRC .OR. MYCOL.EQ.CSRC ) THEN
                  LILO = IWORK(ILILO+WINDOW-1)
                  LIHI = IWORK(ILIHI+WINDOW-1)
                  LSEL = IWORK(ILSEL+WINDOW-1)
               END IF
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  IF( NPCOL.GT.1 .AND. DIR.EQ.1 )
     $               CALL IGEBS2D( ICTXT, 'Row', TOP, 8, 1, IBUFF, 8 )
                  IF( NPROW.GT.1 .AND. DIR.EQ.2 )
     $                 CALL IGEBS2D( ICTXT, 'Col', TOP, 8, 1, IBUFF, 8 )
               ELSEIF( MYROW.EQ.RSRC .OR. MYCOL.EQ.CSRC ) THEN
                  IF( NPCOL.GT.1 .AND. DIR.EQ.1 .AND. MYROW.EQ.RSRC )
     $                 THEN
                     IF( FIRST .OR. (LIHI .GE. LILO + LSEL) ) THEN
                        CALL IGEBR2D( ICTXT, 'Row', TOP, 8, 1, IBUFF, 8,
     $                       RSRC, CSRC )
                        I = IBUFF( 1 )
                        NWIN = IBUFF( 2 )
                        PITRAF = IBUFF( 3 )
                        KS = IBUFF( 4 )
                        PDTRAF = IBUFF( 5 )
                        NDTRAF = IBUFF( 6 )
                        ILEN = IBUFF( 7 )
                        DLEN = IBUFF( 8 )
                     ELSE
                        ILEN = 0
                        DLEN = 0
                        KS = -1
                     END IF
                  END IF
                  IF( NPROW.GT.1 .AND. DIR.EQ.2 .AND. MYCOL.EQ.CSRC )
     $                 THEN
                     IF( FIRST .OR. (LIHI .GE. LILO + LSEL) ) THEN
                        CALL IGEBR2D( ICTXT, 'Col', TOP, 8, 1, IBUFF, 8,
     $                       RSRC, CSRC )
                        I = IBUFF( 1 )
                        NWIN = IBUFF( 2 )
                        PITRAF = IBUFF( 3 )
                        KS = IBUFF( 4 )
                        PDTRAF = IBUFF( 5 )
                        NDTRAF = IBUFF( 6 )
                        ILEN = IBUFF( 7 )
                        DLEN = IBUFF( 8 )
                     ELSE
                        ILEN = 0
                        DLEN = 0
                        KS = -1
                     END IF
                  END IF
               END IF
*
*              Broadcast the accumulated transformations - copy all
*              information from IWORK(IPIW:PITRAF-1) and
*              WORK(IPW2:PDTRAF-1) to a buffer and broadcast this
*              buffer in the corresponding block row and column.  On
*              arrival, copy the information back to the correct part of
*              the workspace. This step is avoided if no computations
*              were performed at the diagonal processor, i.e.,
*              BUFFLEN = 0.
*
               IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                  BUFFER = PDTRAF
                  BUFFLEN = DLEN + ILEN
                  IF( BUFFLEN.NE.0 ) THEN
                     DO 180 INDX = 1, ILEN
                        WORK( BUFFER+INDX-1 ) =
     $                       DBLE( IWORK(IPIW+INDX-1) )
 180                 CONTINUE
                     CALL DLAMOV( 'All', DLEN, 1, WORK( IPW2 ),
     $                    DLEN, WORK(BUFFER+ILEN), DLEN )
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) THEN
                        CALL DGEBS2D( ICTXT, 'Row', TOP, BUFFLEN, 1,
     $                       WORK(BUFFER), BUFFLEN )
                     END IF
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 ) THEN
                        CALL DGEBS2D( ICTXT, 'Col', TOP, BUFFLEN, 1,
     $                       WORK(BUFFER), BUFFLEN )
                     END IF
                  END IF
               ELSEIF( MYROW.EQ.RSRC .OR. MYCOL.EQ.CSRC ) THEN
                  IF( NPCOL.GT.1 .AND. DIR.EQ.1 .AND. MYROW.EQ.RSRC )
     $                 THEN
                     BUFFER = PDTRAF
                     BUFFLEN = DLEN + ILEN
                     IF( BUFFLEN.NE.0 ) THEN
                        CALL DGEBR2D( ICTXT, 'Row', TOP, BUFFLEN, 1,
     $                       WORK(BUFFER), BUFFLEN, RSRC, CSRC )
                     END IF
                  END IF
                  IF( NPROW.GT.1 .AND. DIR.EQ.2 .AND. MYCOL.EQ.CSRC )
     $                 THEN
                     BUFFER = PDTRAF
                     BUFFLEN = DLEN + ILEN
                     IF( BUFFLEN.NE.0 ) THEN
                        CALL DGEBR2D( ICTXT, 'Col', TOP, BUFFLEN, 1,
     $                       WORK(BUFFER), BUFFLEN, RSRC, CSRC )
                     END IF
                  END IF
                  IF((NPCOL.GT.1.AND.DIR.EQ.1.AND.MYROW.EQ.RSRC).OR.
     $                 (NPROW.GT.1.AND.DIR.EQ.2.AND.MYCOL.EQ.CSRC ) )
     $                 THEN
                     IF( BUFFLEN.NE.0 ) THEN
                        DO 190 INDX = 1, ILEN
                           IWORK(IPIW+INDX-1) =
     $                          INT(WORK( BUFFER+INDX-1 ))
 190                    CONTINUE
                        CALL DLAMOV( 'All', DLEN, 1,
     $                       WORK( BUFFER+ILEN ), DLEN,
     $                       WORK( IPW2 ), DLEN )
                     END IF
                  END IF
               END IF
 111        CONTINUE
*
*           Now really perform the updates by applying the orthogonal
*           transformations to the out-of-window parts of T and Q. This
*           step is avoided if no reordering was performed by the on-
*           diagonal processor from the beginning, i.e., BUFFLEN = 0.
*
*           Count number of operations to decide whether to use
*           matrix-matrix multiplications for updating off-diagonal
*           parts or not.
*
            DO 112 WINDOW = 1, NMWIN2
               RSRC = IWORK(IRSRC+WINDOW-1)
               CSRC = IWORK(ICSRC+WINDOW-1)
*
               IF( (MYROW.EQ.RSRC .AND. DIR.EQ.1 ).OR.
     $              (MYCOL.EQ.CSRC .AND. DIR.EQ.2 ) ) THEN
                  LILO = IWORK(ILILO+WINDOW-1)
                  LIHI = IWORK(ILIHI+WINDOW-1)
                  LSEL = IWORK(ILSEL+WINDOW-1)
*
*                 Skip update part for current WINDOW if BUFFLEN = 0.
*
                  IF( BUFFLEN.EQ.0 ) GO TO 295
*
                  NITRAF = PITRAF - IPIW
                  ISHH = .FALSE.
                  FLOPS = 0
                  DO 200 K = 1, NITRAF
                     IF( IWORK( IPIW + K - 1 ).LE.NWIN ) THEN
                        FLOPS = FLOPS + 6
                     ELSE
                        FLOPS = FLOPS + 11
                        ISHH = .TRUE.
                     END IF
 200              CONTINUE
*
*                 Compute amount of work space necessary for performing
*                 matrix-matrix multiplications.
*
                  PDW = BUFFER
                  IPW3 = PDW + NWIN*NWIN
               ELSE
                  FLOPS = 0
               END IF
*
               IF( FLOPS.NE.0 .AND.
     $              ( FLOPS*100 ) / ( 2*NWIN*NWIN ) .GE. MMULT ) THEN
*
*                 Update off-diagonal blocks and Q using matrix-matrix
*                 multiplications; if there are no Householder
*                 reflectors it is preferable to take the triangular
*                 block structure of the transformation matrix into
*                 account.
*
                  CALL DLASET( 'All', NWIN, NWIN, ZERO, ONE,
     $                 WORK( PDW ), NWIN )
                  CALL BDLAAPP( 1, NWIN, NWIN, NCB, WORK( PDW ), NWIN,
     $                 NITRAF, IWORK(IPIW), WORK( IPW2 ), WORK(IPW3) )
*
                  IF( ISHH ) THEN
*
*                    Loop through the local blocks of the distributed
*                    matrices T and Q and update them according to the
*                    performed reordering.
*
*                    Update the columns of T and Q affected by the
*                    reordering.
*
                     IF( DIR.EQ.2 ) THEN
                        DO 210 INDX = 1, I-1, NB
                           CALL INFOG2L( INDX, I, DESCT, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 )
     $                          THEN
                              LROWS = MIN(NB,I-INDX)
                              CALL DGEMM( 'No transpose',
     $                             'No transpose', LROWS, NWIN, NWIN,
     $                             ONE, T((JLOC-1)*LLDT+ILOC), LLDT,
     $                             WORK( PDW ), NWIN, ZERO,
     $                             WORK(IPW3), LROWS )
                              CALL DLAMOV( 'All', LROWS, NWIN,
     $                             WORK(IPW3), LROWS,
     $                             T((JLOC-1)*LLDT+ILOC), LLDT )
                           END IF
 210                    CONTINUE
                        IF( WANTQ ) THEN
                           DO 220 INDX = 1, N, NB
                              CALL INFOG2L( INDX, I, DESCQ, NPROW,
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                             RSRC1, CSRC1 )
                              IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 )
     $                             THEN
                                 LROWS = MIN(NB,N-INDX+1)
                                 CALL DGEMM( 'No transpose',
     $                                'No transpose', LROWS, NWIN, NWIN,
     $                                ONE, Q((JLOC-1)*LLDQ+ILOC), LLDQ,
     $                                WORK( PDW ), NWIN, ZERO,
     $                                WORK(IPW3), LROWS )
                                 CALL DLAMOV( 'All', LROWS, NWIN,
     $                                WORK(IPW3), LROWS,
     $                                Q((JLOC-1)*LLDQ+ILOC), LLDQ )
                              END IF
 220                       CONTINUE
                        END IF
                     END IF
*
*                    Update the rows of T affected by the reordering
*
                     IF( DIR.EQ.1 ) THEN
                        IF( LIHI.LT.N ) THEN
                           IF( MOD(LIHI,NB).GT.0 ) THEN
                              INDX = LIHI + 1
                              CALL INFOG2L( I, INDX, DESCT, NPROW,
     $                            NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                            RSRC1, CSRC1 )
                              IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 )
     $                             THEN
                                 LCOLS = MOD( MIN( NB-MOD(LIHI,NB),
     $                                N-LIHI ), NB )
                                 CALL DGEMM( 'Transpose',
     $                                'No Transpose', NWIN, LCOLS, NWIN,
     $                                ONE, WORK( PDW ), NWIN,
     $                                T((JLOC-1)*LLDT+ILOC), LLDT, ZERO,
     $                                WORK(IPW3), NWIN )
                                 CALL DLAMOV( 'All', NWIN, LCOLS,
     $                                WORK(IPW3), NWIN,
     $                                T((JLOC-1)*LLDT+ILOC), LLDT )
                              END IF
                           END IF
                           INDXS = ICEIL(LIHI,NB)*NB + 1
                           DO 230 INDX = INDXS, N, NB
                              CALL INFOG2L( I, INDX, DESCT, NPROW,
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                             RSRC1, CSRC1 )
                              IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 )
     $                             THEN
                                 LCOLS = MIN( NB, N-INDX+1 )
                                 CALL DGEMM( 'Transpose',
     $                                'No Transpose', NWIN, LCOLS, NWIN,
     $                                ONE, WORK( PDW ), NWIN,
     $                                T((JLOC-1)*LLDT+ILOC), LLDT, ZERO,
     $                                WORK(IPW3), NWIN )
                                 CALL DLAMOV( 'All', NWIN, LCOLS,
     $                                WORK(IPW3), NWIN,
     $                                T((JLOC-1)*LLDT+ILOC), LLDT )
                              END IF
 230                       CONTINUE
                        END IF
                     END IF
                  ELSE
*
*                    The NWIN-by-NWIN matrix U containing the
*                    accumulated orthogonal transformations has the
*                    following structure:
*
*                                  [ U11  U12 ]
*                              U = [          ],
*                                  [ U21  U22 ]
*
*                    where U21 is KS-by-KS upper triangular and U12 is
*                    (NWIN-KS)-by-(NWIN-KS) lower triangular.
*
*                    Update the columns of T and Q affected by the
*                    reordering.
*
*                    Compute T2*U21 + T1*U11 in workspace.
*
                     IF( DIR.EQ.2 ) THEN
                        DO 240 INDX = 1, I-1, NB
                           CALL INFOG2L( INDX, I, DESCT, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 )
     $                          THEN
                              JLOC1 = INDXG2L( I+NWIN-KS, NB, MYCOL,
     $                             DESCT( CSRC_ ), NPCOL )
                              LROWS = MIN(NB,I-INDX)
                              CALL DLAMOV( 'All', LROWS, KS,
     $                             T((JLOC1-1)*LLDT+ILOC ), LLDT,
     $                             WORK(IPW3), LROWS )
                              CALL DTRMM( 'Right', 'Upper',
     $                              'No transpose',
     $                             'Non-unit', LROWS, KS, ONE,
     $                             WORK( PDW+NWIN-KS ), NWIN,
     $                             WORK(IPW3), LROWS )
                              CALL DGEMM( 'No transpose',
     $                             'No transpose', LROWS, KS, NWIN-KS,
     $                             ONE, T((JLOC-1)*LLDT+ILOC), LLDT,
     $                             WORK( PDW ), NWIN, ONE, WORK(IPW3),
     $                             LROWS )
*
*                             Compute T1*U12 + T2*U22 in workspace.
*
                              CALL DLAMOV( 'All', LROWS, NWIN-KS,
     $                             T((JLOC-1)*LLDT+ILOC), LLDT,
     $                             WORK( IPW3+KS*LROWS ), LROWS )
                              CALL DTRMM( 'Right', 'Lower',
     $                             'No transpose', 'Non-unit',
     $                             LROWS, NWIN-KS, ONE,
     $                             WORK( PDW+NWIN*KS ), NWIN,
     $                             WORK( IPW3+KS*LROWS ), LROWS )
                              CALL DGEMM( 'No transpose',
     $                             'No transpose', LROWS, NWIN-KS, KS,
     $                             ONE, T((JLOC1-1)*LLDT+ILOC), LLDT,
     $                             WORK( PDW+NWIN*KS+NWIN-KS ), NWIN,
     $                             ONE, WORK( IPW3+KS*LROWS ), LROWS )
*
*                             Copy workspace to T.
*
                              CALL DLAMOV( 'All', LROWS, NWIN,
     $                             WORK(IPW3), LROWS,
     $                             T((JLOC-1)*LLDT+ILOC), LLDT )
                           END IF
 240                    CONTINUE
                        IF( WANTQ ) THEN
*
*                          Compute Q2*U21 + Q1*U11 in workspace.
*
                           DO 250 INDX = 1, N, NB
                              CALL INFOG2L( INDX, I, DESCQ, NPROW,
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                             RSRC1, CSRC1 )
                              IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 )
     $                             THEN
                                 JLOC1 = INDXG2L( I+NWIN-KS, NB,
     $                                MYCOL, DESCQ( CSRC_ ), NPCOL )
                                 LROWS = MIN(NB,N-INDX+1)
                                 CALL DLAMOV( 'All', LROWS, KS,
     $                                Q((JLOC1-1)*LLDQ+ILOC ), LLDQ,
     $                                WORK(IPW3), LROWS )
                                 CALL DTRMM( 'Right', 'Upper',
     $                                'No transpose', 'Non-unit',
     $                                LROWS, KS, ONE,
     $                                WORK( PDW+NWIN-KS ), NWIN,
     $                                WORK(IPW3), LROWS )
                                 CALL DGEMM( 'No transpose',
     $                                'No transpose', LROWS, KS,
     $                                NWIN-KS, ONE,
     $                                Q((JLOC-1)*LLDQ+ILOC), LLDQ,
     $                                WORK( PDW ), NWIN, ONE,
     $                                WORK(IPW3), LROWS )
*
*                                Compute Q1*U12 + Q2*U22 in workspace.
*
                                 CALL DLAMOV( 'All', LROWS, NWIN-KS,
     $                                Q((JLOC-1)*LLDQ+ILOC), LLDQ,
     $                                WORK( IPW3+KS*LROWS ), LROWS)
                                 CALL DTRMM( 'Right', 'Lower',
     $                                'No transpose', 'Non-unit',
     $                                LROWS, NWIN-KS, ONE,
     $                                WORK( PDW+NWIN*KS ), NWIN,
     $                                WORK( IPW3+KS*LROWS ), LROWS)
                                 CALL DGEMM( 'No transpose',
     $                                'No transpose', LROWS, NWIN-KS,
     $                                KS, ONE, Q((JLOC1-1)*LLDQ+ILOC),
     $                                LLDQ, WORK(PDW+NWIN*KS+NWIN-KS),
     $                                NWIN, ONE, WORK( IPW3+KS*LROWS ),
     $                                LROWS )
*
*                                Copy workspace to Q.
*
                                 CALL DLAMOV( 'All', LROWS, NWIN,
     $                                WORK(IPW3), LROWS,
     $                                Q((JLOC-1)*LLDQ+ILOC), LLDQ )
                              END IF
 250                       CONTINUE
                        END IF
                     END IF
*
                     IF( DIR.EQ.1 ) THEN
                        IF ( LIHI.LT.N ) THEN
*
*                          Compute U21**T*T2 + U11**T*T1 in workspace.
*
                           IF( MOD(LIHI,NB).GT.0 ) THEN
                              INDX = LIHI + 1
                              CALL INFOG2L( I, INDX, DESCT, NPROW,
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                             RSRC1, CSRC1 )
                              IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 )
     $                             THEN
                                 ILOC1 = INDXG2L( I+NWIN-KS, NB, MYROW,
     $                                DESCT( RSRC_ ), NPROW )
                                 LCOLS = MOD( MIN( NB-MOD(LIHI,NB),
     $                                N-LIHI ), NB )
                                 CALL DLAMOV( 'All', KS, LCOLS,
     $                                T((JLOC-1)*LLDT+ILOC1), LLDT,
     $                                WORK(IPW3), NWIN )
                                 CALL DTRMM( 'Left', 'Upper',
     $                                'Transpose', 'Non-unit', KS,
     $                                LCOLS, ONE, WORK( PDW+NWIN-KS ),
     $                                NWIN, WORK(IPW3), NWIN )
                                 CALL DGEMM( 'Transpose',
     $                                'No transpose', KS, LCOLS,
     $                                NWIN-KS, ONE, WORK(PDW), NWIN,
     $                                T((JLOC-1)*LLDT+ILOC), LLDT, ONE,
     $                                WORK(IPW3), NWIN )
*
*                                Compute U12**T*T1 + U22**T*T2 in
*                                workspace.
*
                                 CALL DLAMOV( 'All', NWIN-KS, LCOLS,
     $                                T((JLOC-1)*LLDT+ILOC), LLDT,
     $                                WORK( IPW3+KS ), NWIN )
                                 CALL DTRMM( 'Left', 'Lower',
     $                                'Transpose', 'Non-unit',
     $                                NWIN-KS, LCOLS, ONE,
     $                                WORK( PDW+NWIN*KS ), NWIN,
     $                                WORK( IPW3+KS ), NWIN )
                                 CALL DGEMM( 'Transpose',
     $                                'No Transpose', NWIN-KS, LCOLS,
     $                                KS, ONE,
     $                                WORK( PDW+NWIN*KS+NWIN-KS ),
     $                                NWIN, T((JLOC-1)*LLDT+ILOC1),
     $                                LLDT, ONE, WORK( IPW3+KS ),
     $                                NWIN )
*
*                                Copy workspace to T.
*
                                 CALL DLAMOV( 'All', NWIN, LCOLS,
     $                                WORK(IPW3), NWIN,
     $                                T((JLOC-1)*LLDT+ILOC), LLDT )
                              END IF
                           END IF
                           INDXS = ICEIL(LIHI,NB)*NB + 1
                           DO 260 INDX = INDXS, N, NB
                              CALL INFOG2L( I, INDX, DESCT, NPROW,
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                             RSRC1, CSRC1 )
                              IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC1 )
     $                             THEN
*
*                                Compute U21**T*T2 + U11**T*T1 in
*                                workspace.
*
                                 ILOC1 = INDXG2L( I+NWIN-KS, NB,
     $                                MYROW, DESCT( RSRC_ ), NPROW )
                                 LCOLS = MIN( NB, N-INDX+1 )
                                 CALL DLAMOV( 'All', KS, LCOLS,
     $                                T((JLOC-1)*LLDT+ILOC1), LLDT,
     $                                WORK(IPW3), NWIN )
                                 CALL DTRMM( 'Left', 'Upper',
     $                                'Transpose', 'Non-unit', KS,
     $                                LCOLS, ONE,
     $                                WORK( PDW+NWIN-KS ), NWIN,
     $                                WORK(IPW3), NWIN )
                                 CALL DGEMM( 'Transpose',
     $                                'No transpose', KS, LCOLS,
     $                                NWIN-KS, ONE, WORK(PDW), NWIN,
     $                                T((JLOC-1)*LLDT+ILOC), LLDT, ONE,
     $                                WORK(IPW3), NWIN )
*
*                                Compute U12**T*T1 + U22**T*T2 in
*                                workspace.
*
                                 CALL DLAMOV( 'All', NWIN-KS, LCOLS,
     $                                T((JLOC-1)*LLDT+ILOC), LLDT,
     $                                WORK( IPW3+KS ), NWIN )
                                 CALL DTRMM( 'Left', 'Lower',
     $                                'Transpose', 'Non-unit',
     $                                NWIN-KS, LCOLS, ONE,
     $                                WORK( PDW+NWIN*KS ), NWIN,
     $                                WORK( IPW3+KS ), NWIN )
                                 CALL DGEMM( 'Transpose',
     $                                'No Transpose', NWIN-KS, LCOLS,
     $                                KS, ONE,
     $                                WORK( PDW+NWIN*KS+NWIN-KS ),
     $                                NWIN, T((JLOC-1)*LLDT+ILOC1),
     $                                LLDT, ONE, WORK(IPW3+KS), NWIN )
*
*                                Copy workspace to T.
*
                                 CALL DLAMOV( 'All', NWIN, LCOLS,
     $                                WORK(IPW3), NWIN,
     $                                T((JLOC-1)*LLDT+ILOC), LLDT )
                              END IF
 260                       CONTINUE
                        END IF
                     END IF
                  END IF
               ELSEIF( FLOPS.NE.0 ) THEN
*
*                 Update off-diagonal blocks and Q using the pipelined
*                 elementary transformations.
*
                  IF( DIR.EQ.2 ) THEN
                     DO 270 INDX = 1, I-1, NB
                        CALL INFOG2L( INDX, I, DESCT, NPROW, NPCOL,
     $                       MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                        IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                           LROWS = MIN(NB,I-INDX)
                           CALL BDLAAPP( 1, LROWS, NWIN, NCB,
     $                          T((JLOC-1)*LLDT+ILOC ), LLDT, NITRAF,
     $                          IWORK(IPIW), WORK( IPW2 ),
     $                          WORK(IPW3) )
                        END IF
 270                 CONTINUE
                     IF( WANTQ ) THEN
                        DO 280 INDX = 1, N, NB
                           CALL INFOG2L( INDX, I, DESCQ, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 )
     $                          THEN
                              LROWS = MIN(NB,N-INDX+1)
                              CALL BDLAAPP( 1, LROWS, NWIN, NCB,
     $                             Q((JLOC-1)*LLDQ+ILOC), LLDQ, NITRAF,
     $                             IWORK(IPIW), WORK( IPW2 ),
     $                             WORK(IPW3) )
                           END IF
 280                    CONTINUE
                     END IF
                  END IF
                  IF( DIR.EQ.1 ) THEN
                     IF( LIHI.LT.N ) THEN
                        IF( MOD(LIHI,NB).GT.0 ) THEN
                           INDX = LIHI + 1
                           CALL INFOG2L( I, INDX, DESCT, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 )
     $                          THEN
                              LCOLS = MOD( MIN( NB-MOD(LIHI,NB),
     $                             N-LIHI ), NB )
                              CALL BDLAAPP( 0, NWIN, LCOLS, NCB,
     $                             T((JLOC-1)*LLDT+ILOC), LLDT, NITRAF,
     $                             IWORK(IPIW), WORK( IPW2 ),
     $                             WORK(IPW3) )
                           END IF
                        END IF
                        INDXS = ICEIL(LIHI,NB)*NB + 1
                        DO 290 INDX = INDXS, N, NB
                           CALL INFOG2L( I, INDX, DESCT, NPROW, NPCOL,
     $                          MYROW, MYCOL, ILOC, JLOC, RSRC1, CSRC1 )
                           IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 )
     $                          THEN
                              LCOLS = MIN( NB, N-INDX+1 )
                              CALL BDLAAPP( 0, NWIN, LCOLS, NCB,
     $                             T((JLOC-1)*LLDT+ILOC), LLDT, NITRAF,
     $                             IWORK(IPIW), WORK( IPW2 ),
     $                             WORK(IPW3) )
                           END IF
 290                    CONTINUE
                     END IF
                  END IF
               END IF
*
*              If I was not involved in the updates for the current
*              window or the window was fully processed, I go here and
*              try again for the next window.
*
 295           CONTINUE
*
*              Update LIHI and LIHI depending on the number of
*              eigenvalues really moved - for on-diagonal processes we
*              do this update only once since each on-diagonal process
*              is only involved with one window at one time. The
*              indicies are updated in three cases:
*                1) When some reordering was really performed
*                   -- indicated by BUFFLEN > 0.
*                2) When no selected eigenvalues was found in the
*                   current window -- indicated by KS = 0.
*                3) When some selected eigenvalues was found in the
*                   current window but no one of them was moved
*                   (KS > 0 and BUFFLEN = 0)
*              False index updating is avoided by sometimes setting
*              KS = -1. This will affect processors involved in more
*              than one window and where the first one ends up with
*              KS = 0 and for the second one is done already.
*
               IF( MYROW.EQ.RSRC.AND.MYCOL.EQ.CSRC ) THEN
                  IF( DIR.EQ.2 ) THEN
                     IF( BUFFLEN.NE.0 .OR. KS.EQ.0 .OR.
     $                    ( BUFFLEN.EQ.0 .AND. KS.GT.0 ) )
     $                  LIHI = I + KS - 1
                     IWORK( ILIHI+WINDOW-1 ) = LIHI
                     IF( .NOT. LIHI.GE.LILO+LSEL ) THEN
                        LILO = LILO + LSEL
                        IWORK( ILILO+WINDOW-1 ) = LILO
                     END IF
                  END IF
               ELSEIF( MYROW.EQ.RSRC .AND. DIR.EQ.1 ) THEN
                  IF( BUFFLEN.NE.0 .OR. KS.EQ.0 .OR.
     $                 ( BUFFLEN.EQ.0 .AND. KS.GT.0 ) )
     $               LIHI = I + KS - 1
                  IWORK( ILIHI+WINDOW-1 ) = LIHI
                  IF( .NOT. LIHI.GE.LILO+LSEL ) THEN
                     LILO = LILO + LSEL
                     IWORK( ILILO+WINDOW-1 ) = LILO
                  END IF
               ELSEIF( MYCOL.EQ.CSRC .AND. DIR.EQ.2 ) THEN
                  IF( BUFFLEN.NE.0 .OR. KS.EQ.0 .OR.
     $                 ( BUFFLEN.EQ.0 .AND. KS.GT.0 ) )
     $               LIHI = I + KS - 1
                  IWORK( ILIHI+WINDOW-1 ) = LIHI
                  IF( .NOT. LIHI.GE.LILO+LSEL ) THEN
                     LILO = LILO + LSEL
                     IWORK( ILILO+WINDOW-1 ) = LILO
                  END IF
               END IF
*
 112        CONTINUE
*
*           End of direction loop for updates with respect to local
*           reordering.
*
 1111       CONTINUE
*
*           Associate each process with one of the corresponding
*           computational windows such that the test for another round
*           of local reordering is carried out properly. Since the
*           column updates were computed after the row updates, it is
*           sufficient to test for changing the association to the
*           window in the corresponding process row.
*
            DO 113 WINDOW = 1, NMWIN2
               RSRC = IWORK( IRSRC + WINDOW - 1 )
               IF( MYROW.EQ.RSRC .AND. (.NOT. LIHI.GE.LILO+LSEL ) ) THEN
                  LILO = IWORK( ILILO + WINDOW - 1 )
                  LIHI = IWORK( ILIHI + WINDOW - 1 )
                  LSEL = IWORK( ILSEL + WINDOW - 1 )
               END IF
 113        CONTINUE
*
*           End While ( LIHI >= LILO + LSEL )
            ROUND = ROUND + 1
            IF( FIRST ) FIRST = .FALSE.
            GO TO 130
         END IF
*
*        All processors excluded from the local reordering go here.
*
 114     CONTINUE
*
*        Barrier to collect the processes before proceeding.
*
         CALL BLACS_BARRIER( ICTXT, 'All' )
*
*        Compute global maximum of IERR so that we know if some process
*        experienced a failure in the reordering.
*
         MYIERR = IERR
         IF( NPROCS.GT.1 )
     $      CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, IERR, 1, -1,
     $           -1, -1, -1, -1 )
*
         IF( IERR.NE.0 ) THEN
*
*           When calling BDTREXC, the block at position I+KKS-1 failed
*           to swap.
*
            IF( MYIERR.NE.0 ) INFO = MAX(1,I+KKS-1)
            IF( NPROCS.GT.1 )
     $         CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, INFO, 1, -1,
     $              -1, -1, -1, -1 )
            GO TO 300
         END IF
*
*        Now, for each compuational window, move the selected
*        eigenvalues across the process border. Do this by forming the
*        processors into groups of four working together to bring the
*        window over the border. The processes are numbered as follows
*
*                1 | 2
*                --+--
*                3 | 4
*
*        where '|' and '-' denotes the process (and block) borders.
*        This implies that the cluster to be reordered over the border
*        is held by process 4, process 1 will receive the cluster after
*        the reordering, process 3 holds the local (2,1)th element of a
*        2-by-2 diagonal block located on the block border and process 2
*        holds the closest off-diagonal part of the window that is
*        affected by the cross-border reordering.
*
*        The active window is now ( I : LIHI[4], I : LIHI[4] ), where
*        I = MAX( ILO, LIHI - 2*MOD(LIHI,NB) ). If this active window is
*        too large compared to the value of PARA( 6 ), it will be
*        truncated in both ends such that a maximum of PARA( 6 )
*        eigenvalues is reordered across the border this time.
*
*        The active window will be collected and built in workspace at
*        process 1 and 4, which both compute the reordering and return
*        the updated parts to the corresponding processes 2-3. Next, the
*        accumulated transformations are broadcasted for updates in the
*        block rows and column that corresponds to the process rows and
*        columns where process 1 and 4 reside.
*
*        The off-diagonal blocks are updated by the processes receiving
*        from the broadcasts of the orthogonal transformations. Since
*        the active window is split over the process borders, the
*        updates of T and Q requires that stripes of block rows of
*        columns are exchanged between neighboring processes in the
*        corresponding process rows and columns.
*
*        First, form each group of processors involved in the
*        crossborder reordering. Do this in two (or three) phases:
*        1) Reorder each odd window over the border.
*        2) Reorder each even window over the border.
*        3) Reorder the last odd window over the border, if it was not
*           processed in the first phase.
*
*        When reordering the odd windows over the border, we must make
*        sure that no process row or column is involved in both the
*        first and the last window at the same time. This happens when
*        the total number of windows is odd, greater than one and equal
*        to the minumum process mesh dimension. Therefore the last odd
*        window may be reordered over the border at last.
*
         LASTWAIT = NMWIN2.GT.1 .AND. MOD(NMWIN2,2).EQ.1 .AND.
     $        NMWIN2.EQ.MIN(NPROW,NPCOL)
*
         LAST = 0
 308     CONTINUE
         IF( LASTWAIT ) THEN
            IF( LAST.EQ.0 ) THEN
               WIN0S = 1
               WIN0E = 2
               WINE = NMWIN2 - 1
            ELSE
               WIN0S = NMWIN2
               WIN0E = NMWIN2
               WINE = NMWIN2
            END IF
         ELSE
            WIN0S = 1
            WIN0E = 2
            WINE = NMWIN2
         END IF
         DO 310 WINDOW0 = WIN0S, WIN0E
            DO 320 WINDOW = WINDOW0, WINE, 2
*
*              Define the process holding the down-right part of the
*              window.
*
               RSRC4 = IWORK(IRSRC+WINDOW-1)
               CSRC4 = IWORK(ICSRC+WINDOW-1)
*
*              Define the other processes in the group of four.
*
               RSRC3 = RSRC4
               CSRC3 = MOD( CSRC4 - 1 + NPCOL, NPCOL )
               RSRC2 = MOD( RSRC4 - 1 + NPROW, NPROW )
               CSRC2 = CSRC4
               RSRC1 = RSRC2
               CSRC1 = CSRC3
               IF( ( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) .OR.
     $             ( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) .OR.
     $             ( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) .OR.
     $             ( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) ) THEN
*
*                 Compute the correct active window - for reordering
*                 into a block that has not been active at all before,
*                 we try to reorder as many of our eigenvalues over the
*                 border as possible without knowing of the situation on
*                 the other side - this may cause very few eigenvalues
*                 to be reordered over the border this time (perhaps not
*                 any) but this should be an initial problem.  Anyway,
*                 the bottom-right position of the block will be at
*                 position LIHIC.
*
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     LIHI4 = ( IWORK( ILILO + WINDOW - 1 ) +
     $                    IWORK( ILIHI + WINDOW - 1 ) ) / 2
                     LIHIC = MIN(LIHI4,(ICEIL(LIHI4,NB)-1)*NB+WNEICR)
*
*                    Fix LIHIC to avoid that bottom of window cuts
*                    2-by-2 block and make sure all processors in the
*                    group knows about the correct value.
*
                     IF( (.NOT. LIHIC.LE.NB) .AND. LIHIC.LT.N ) THEN
                        ILOC = INDXG2L( LIHIC+1, NB, MYROW,
     $                       DESCT( RSRC_ ), NPROW )
                        JLOC = INDXG2L( LIHIC, NB, MYCOL,
     $                       DESCT( CSRC_ ), NPCOL )
                        IF( T( (JLOC-1)*LLDT+ILOC ).NE.ZERO ) THEN
                           IF( MOD( LIHIC, NB ).EQ.1 .OR.
     $                          ( MOD( LIHIC, NB ).EQ.2 .AND.
     $                          SELECT(LIHIC-2).EQ.0 ) )
     $                          THEN
                              LIHIC = LIHIC + 1
                           ELSE
                              LIHIC = LIHIC - 1
                           END IF
                        END IF
                     END IF
                     IF( RSRC4.NE.RSRC1 .OR. CSRC4.NE.CSRC1 )
     $                  CALL IGESD2D( ICTXT, 1, 1, LIHIC, 1, RSRC1,
     $                       CSRC1 )
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 )
     $                  CALL IGESD2D( ICTXT, 1, 1, LIHIC, 1, RSRC2,
     $                       CSRC2 )
                     IF( RSRC4.NE.RSRC3 .OR. CSRC4.NE.CSRC3 )
     $                  CALL IGESD2D( ICTXT, 1, 1, LIHIC, 1, RSRC3,
     $                       CSRC3 )
                  END IF
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( RSRC4.NE.RSRC1 .OR. CSRC4.NE.CSRC1 )
     $                  CALL IGERV2D( ICTXT, 1, 1, LIHIC, 1, RSRC4,
     $                       CSRC4 )
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 )
     $                  CALL IGERV2D( ICTXT, 1, 1, LIHIC, 1, RSRC4,
     $                       CSRC4 )
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     IF( RSRC4.NE.RSRC3 .OR. CSRC4.NE.CSRC3 )
     $                  CALL IGERV2D( ICTXT, 1, 1, LIHIC, 1, RSRC4,
     $                       CSRC4 )
                  END IF
*
*                 Avoid going over the border with the first window if
*                 it resides in the block where the last global position
*                 T(ILO,ILO) is or ILO has been updated to point to a
*                 position right of T(LIHIC,LIHIC).
*
                  SKIP1CR = WINDOW.EQ.1 .AND.
     $                 ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
*
*                 Decide I, where to put top of window, such that top of
*                 window does not cut 2-by-2 block. Make sure that we do
*                 not end up in a situation where a 2-by-2 block
*                 splitted on the border is left in its original place
*                 -- this can cause infinite loops.
*                 Remedy: make sure that the part of the window that
*                 resides left to the border is at least of dimension
*                 two (2) in case we have 2-by-2 blocks in top of the
*                 cross border window.
*
*                 Also make sure all processors in the group knows about
*                 the correct value of I. When skipping the crossborder
*                 reordering, just set I = LIHIC.
*
                  IF( .NOT. SKIP1CR ) THEN
                     IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                        IF( WINDOW.EQ.1 ) THEN
                           LIHI1 = ILO
                        ELSE
                           LIHI1 = IWORK( ILIHI + WINDOW - 2 )
                        END IF
                        I = MAX( LIHI1,
     $                       MIN( LIHIC-2*MOD(LIHIC,NB) + 1,
     $                       (ICEIL(LIHIC,NB)-1)*NB - 1  ) )
                        ILOC = INDXG2L( I, NB, MYROW, DESCT( RSRC_ ),
     $                       NPROW )
                        JLOC = INDXG2L( I-1, NB, MYCOL, DESCT( CSRC_ ),
     $                       NPCOL )
                        IF( T( (JLOC-1)*LLDT+ILOC ).NE.ZERO )
     $                     I = I - 1
                        IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 )
     $                     CALL IGESD2D( ICTXT, 1, 1, I, 1, RSRC4,
     $                          CSRC4 )
                        IF( RSRC1.NE.RSRC2 .OR. CSRC1.NE.CSRC2 )
     $                     CALL IGESD2D( ICTXT, 1, 1, I, 1, RSRC2,
     $                          CSRC2 )
                        IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 )
     $                     CALL IGESD2D( ICTXT, 1, 1, I, 1, RSRC3,
     $                          CSRC3 )
                     END IF
                     IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                        IF( RSRC1.NE.RSRC2 .OR. CSRC1.NE.CSRC2 )
     $                     CALL IGERV2D( ICTXT, 1, 1, I, 1, RSRC1,
     $                          CSRC1 )
                     END IF
                     IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                        IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 )
     $                     CALL IGERV2D( ICTXT, 1, 1, I, 1, RSRC1,
     $                          CSRC1 )
                     END IF
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 )
     $                     CALL IGERV2D( ICTXT, 1, 1, I, 1, RSRC1,
     $                          CSRC1 )
                     END IF
                  ELSE
                     I = LIHIC
                  END IF
*
*                 Finalize computation of window size: active window is
*                 now (I:LIHIC,I:LIHIC).
*
                  NWIN = LIHIC - I + 1
                  KS = 0
*
*                 Skip rest of this part if appropriate.
*
                  IF( SKIP1CR ) GO TO 360
*
*                 Divide workspace -- put active window in
*                 WORK(IPW2:IPW2+NWIN**2-1) and orthogonal
*                 transformations in WORK(IPW3:...).
*
                  CALL DLASET( 'All', NWIN, NWIN, ZERO, ZERO,
     $                 WORK( IPW2 ), NWIN )
*
                  PITRAF = IPIW
                  IPW3 = IPW2 + NWIN*NWIN
                  PDTRAF = IPW3
*
*                 Exchange the current view of SELECT for the active
*                 window between process 1 and 4 to make sure that
*                 exactly the same job is performed for both processes.
*
                  IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) THEN
                     ILEN4 = MOD(LIHIC,NB)
                     SELI4 = ICEIL(I,NB)*NB+1
                     ILEN1 = NWIN - ILEN4
                     IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                        CALL IGESD2D( ICTXT, ILEN1, 1, SELECT(I),
     $                       ILEN1, RSRC4, CSRC4 )
                        CALL IGERV2D( ICTXT, ILEN4, 1, SELECT(SELI4),
     $                       ILEN4, RSRC4, CSRC4 )
                     END IF
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        CALL IGESD2D( ICTXT, ILEN4, 1, SELECT(SELI4),
     $                       ILEN4, RSRC1, CSRC1 )
                        CALL IGERV2D( ICTXT, ILEN1, 1, SELECT(I),
     $                       ILEN1, RSRC1, CSRC1 )
                     END IF
                  END IF
*
*                 Form the active window by a series of point-to-point
*                 sends and receives.
*
                  DIM1 = NB - MOD(I-1,NB)
                  DIM4 = NWIN - DIM1
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     ILOC = INDXG2L( I, NB, MYROW, DESCT( RSRC_ ),
     $                    NPROW )
                     JLOC = INDXG2L( I, NB, MYCOL, DESCT( CSRC_ ),
     $                    NPCOL )
                     CALL DLAMOV( 'All', DIM1, DIM1,
     $                    T((JLOC-1)*LLDT+ILOC), LLDT, WORK(IPW2),
     $                    NWIN )
                     IF( RSRC1.NE.RSRC4 .OR. CSRC1.NE.CSRC4 ) THEN
                        CALL DGESD2D( ICTXT, DIM1, DIM1,
     $                       WORK(IPW2), NWIN, RSRC4, CSRC4 )
                        CALL DGERV2D( ICTXT, DIM4, DIM4,
     $                       WORK(IPW2+DIM1*NWIN+DIM1), NWIN, RSRC4,
     $                       CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     ILOC = INDXG2L( I+DIM1, NB, MYROW, DESCT( RSRC_ ),
     $                    NPROW )
                     JLOC = INDXG2L( I+DIM1, NB, MYCOL, DESCT( CSRC_ ),
     $                    NPCOL )
                     CALL DLAMOV( 'All', DIM4, DIM4,
     $                    T((JLOC-1)*LLDT+ILOC), LLDT,
     $                    WORK(IPW2+DIM1*NWIN+DIM1), NWIN )
                     IF( RSRC4.NE.RSRC1 .OR. CSRC4.NE.CSRC1 ) THEN
                        CALL DGESD2D( ICTXT, DIM4, DIM4,
     $                       WORK(IPW2+DIM1*NWIN+DIM1), NWIN, RSRC1,
     $                       CSRC1 )
                        CALL DGERV2D( ICTXT, DIM1, DIM1,
     $                       WORK(IPW2), NWIN, RSRC1, CSRC1 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                     ILOC = INDXG2L( I, NB, MYROW, DESCT( RSRC_ ),
     $                    NPROW )
                     JLOC = INDXG2L( I+DIM1, NB, MYCOL, DESCT( CSRC_ ),
     $                    NPCOL )
                     CALL DLAMOV( 'All', DIM1, DIM4,
     $                    T((JLOC-1)*LLDT+ILOC), LLDT,
     $                    WORK(IPW2+DIM1*NWIN), NWIN )
                     IF( RSRC2.NE.RSRC1 .OR. CSRC2.NE.CSRC1 ) THEN
                        CALL DGESD2D( ICTXT, DIM1, DIM4,
     $                       WORK(IPW2+DIM1*NWIN), NWIN, RSRC1, CSRC1 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                     IF( RSRC2.NE.RSRC4 .OR. CSRC2.NE.CSRC4 ) THEN
                        CALL DGESD2D( ICTXT, DIM1, DIM4,
     $                       WORK(IPW2+DIM1*NWIN), NWIN, RSRC4, CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     ILOC = INDXG2L( I+DIM1, NB, MYROW, DESCT( RSRC_ ),
     $                    NPROW )
                     JLOC = INDXG2L( I+DIM1-1, NB, MYCOL,
     $                    DESCT( CSRC_ ), NPCOL )
                     CALL DLAMOV( 'All', 1, 1,
     $                    T((JLOC-1)*LLDT+ILOC), LLDT,
     $                    WORK(IPW2+(DIM1-1)*NWIN+DIM1), NWIN )
                     IF( RSRC3.NE.RSRC1 .OR. CSRC3.NE.CSRC1 ) THEN
                        CALL DGESD2D( ICTXT, 1, 1,
     $                       WORK(IPW2+(DIM1-1)*NWIN+DIM1), NWIN,
     $                       RSRC1, CSRC1 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     IF( RSRC3.NE.RSRC4 .OR. CSRC3.NE.CSRC4 ) THEN
                        CALL DGESD2D( ICTXT, 1, 1,
     $                       WORK(IPW2+(DIM1-1)*NWIN+DIM1), NWIN,
     $                       RSRC4, CSRC4 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( RSRC1.NE.RSRC2 .OR. CSRC1.NE.CSRC2 ) THEN
                        CALL DGERV2D( ICTXT, DIM1, DIM4,
     $                       WORK(IPW2+DIM1*NWIN), NWIN, RSRC2,
     $                       CSRC2 )
                     END IF
                     IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) THEN
                        CALL DGERV2D( ICTXT, 1, 1,
     $                       WORK(IPW2+(DIM1-1)*NWIN+DIM1), NWIN,
     $                       RSRC3, CSRC3 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 ) THEN
                        CALL DGERV2D( ICTXT, DIM1, DIM4,
     $                       WORK(IPW2+DIM1*NWIN), NWIN, RSRC2,
     $                       CSRC2 )
                     END IF
                     IF( RSRC4.NE.RSRC3 .OR. CSRC4.NE.CSRC3 ) THEN
                        CALL DGERV2D( ICTXT, 1, 1,
     $                       WORK(IPW2+(DIM1-1)*NWIN+DIM1), NWIN,
     $                       RSRC3, CSRC3 )
                     END IF
                  END IF
*
*                 Compute the reordering (just as in the total local
*                 case) and accumulate the transformations (ONLY
*                 ON-DIAGONAL PROCESSES).
*
                  IF( ( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) .OR.
     $                ( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) ) THEN
                     PAIR = .FALSE.
                     DO 330 K = I, LIHIC
                        IF( PAIR ) THEN
                           PAIR = .FALSE.
                        ELSE
                           SWAP = SELECT( K ).NE.0
                           IF( K.LT.LIHIC ) THEN
                              ELEM = WORK(IPW2+(K-I)*NWIN+K-I+1)
                              IF( ELEM.NE.ZERO )
     $                           PAIR = .TRUE.
                           END IF
                           IF( SWAP ) THEN
                              KS = KS + 1
*
*                             Swap the K-th block to position I+KS-1.
*
                              IERR = 0
                              KK  = K - I + 1
                              KKS = KS
                              IF( KK.NE.KS ) THEN
                                 NITRAF = LIWORK - PITRAF + 1
                                 NDTRAF = LWORK - PDTRAF + 1
                                 CALL BDTREXC( NWIN, WORK(IPW2), NWIN,
     $                                KK, KKS, NITRAF, IWORK( PITRAF ),
     $                                NDTRAF, WORK( PDTRAF ),
     $                                WORK(IPW1), IERR )
                                 PITRAF = PITRAF + NITRAF
                                 PDTRAF = PDTRAF + NDTRAF
*
*                                Update array SELECT.
*
                                 IF ( PAIR ) THEN
                                    DO 340 J = I+KK-1, I+KKS, -1
                                       SELECT(J+1) = SELECT(J-1)
 340                                CONTINUE
                                    SELECT(I+KKS-1) = 1
                                    SELECT(I+KKS) = 1
                                 ELSE
                                    DO 350 J = I+KK-1, I+KKS, -1
                                       SELECT(J) = SELECT(J-1)
 350                                CONTINUE
                                    SELECT(I+KKS-1) = 1
                                 END IF
*
                                 IF ( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
*
                                    IF ( IERR.EQ.2 ) THEN
                                       SELECT( I+KKS-3 ) = 1
                                       SELECT( I+KKS-1 ) = 0
                                       KKS = KKS + 1
                                    END IF
*
                                    GO TO 360
                                 END IF
                                 KS = KKS
                              END IF
                              IF( PAIR )
     $                           KS = KS + 1
                           END IF
                        END IF
 330                 CONTINUE
                  END IF
 360              CONTINUE
*
*                 Save information about the reordering.
*
                  IF( ( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) .OR.
     $                 ( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) ) THEN
                     IBUFF( 1 ) = I
                     IBUFF( 2 ) = NWIN
                     IBUFF( 3 ) = PITRAF
                     IBUFF( 4 ) = KS
                     IBUFF( 5 ) = PDTRAF
                     IBUFF( 6 ) = NDTRAF
                     ILEN = PITRAF - IPIW + 1
                     DLEN = PDTRAF - IPW3 + 1
                     IBUFF( 7 ) = ILEN
                     IBUFF( 8 ) = DLEN
*
*                    Put reordered data back into global matrix if a
*                    reordering took place.
*
                     IF( .NOT. SKIP1CR ) THEN
                        IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                           ILOC = INDXG2L( I, NB, MYROW, DESCT( RSRC_ ),
     $                          NPROW )
                           JLOC = INDXG2L( I, NB, MYCOL, DESCT( CSRC_ ),
     $                          NPCOL )
                           CALL DLAMOV( 'All', DIM1, DIM1, WORK(IPW2),
     $                          NWIN, T((JLOC-1)*LLDT+ILOC), LLDT )
                        END IF
                        IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                           ILOC = INDXG2L( I+DIM1, NB, MYROW,
     $                          DESCT( RSRC_ ), NPROW )
                           JLOC = INDXG2L( I+DIM1, NB, MYCOL,
     $                          DESCT( CSRC_ ), NPCOL )
                           CALL DLAMOV( 'All', DIM4, DIM4,
     $                          WORK(IPW2+DIM1*NWIN+DIM1), NWIN,
     $                          T((JLOC-1)*LLDT+ILOC), LLDT )
                        END IF
                     END IF
                  END IF
*
*                 Break if appropriate -- IBUFF(3:8) may now contain
*                 nonsens, but that's no problem. The processors outside
*                 the cross border group only needs to know about I and
*                 NWIN to get a correct value of SKIP1CR (see below) and
*                 to skip the cross border updates if necessary.
*
                  IF( WINDOW.EQ.1 .AND. SKIP1CR ) GO TO 325
*
*                 Return reordered data to process 2 and 3.
*
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( RSRC1.NE.RSRC3 .OR. CSRC1.NE.CSRC3 ) THEN
                        CALL DGESD2D( ICTXT, 1, 1,
     $                       WORK( IPW2+(DIM1-1)*NWIN+DIM1 ), NWIN,
     $                       RSRC3, CSRC3 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                     IF( RSRC4.NE.RSRC2 .OR. CSRC4.NE.CSRC2 ) THEN
                        CALL DGESD2D( ICTXT, DIM1, DIM4,
     $                       WORK( IPW2+DIM1*NWIN), NWIN, RSRC2,
     $                       CSRC2 )
                     END IF
                  END IF
                  IF( MYROW.EQ.RSRC2 .AND. MYCOL.EQ.CSRC2 ) THEN
                     ILOC = INDXG2L( I, NB, MYROW, DESCT( RSRC_ ),
     $                    NPROW )
                     JLOC = INDXG2L( I+DIM1, NB, MYCOL,
     $                    DESCT( CSRC_ ), NPCOL )
                     IF( RSRC2.NE.RSRC4 .OR. CSRC2.NE.CSRC4 ) THEN
                        CALL DGERV2D( ICTXT, DIM1, DIM4,
     $                       WORK(IPW2+DIM1*NWIN), NWIN, RSRC4, CSRC4 )
                     END IF
                     CALL DLAMOV( 'All', DIM1, DIM4,
     $                    WORK( IPW2+DIM1*NWIN ), NWIN,
     $                    T((JLOC-1)*LLDT+ILOC), LLDT )
                  END IF
                  IF( MYROW.EQ.RSRC3 .AND. MYCOL.EQ.CSRC3 ) THEN
                     ILOC = INDXG2L( I+DIM1, NB, MYROW,
     $                    DESCT( RSRC_ ), NPROW )
                     JLOC = INDXG2L( I+DIM1-1, NB, MYCOL,
     $                    DESCT( CSRC_ ), NPCOL )
                     IF( RSRC3.NE.RSRC1 .OR. CSRC3.NE.CSRC1 ) THEN
                        CALL DGERV2D( ICTXT, 1, 1,
     $                       WORK( IPW2+(DIM1-1)*NWIN+DIM1 ), NWIN,
     $                       RSRC1, CSRC1 )
                     END IF
                     T((JLOC-1)*LLDT+ILOC) =
     $                    WORK( IPW2+(DIM1-1)*NWIN+DIM1 )
                  END IF
               END IF
*
 325           CONTINUE
*
 320        CONTINUE
*
*           For the crossborder updates, we use the same directions as
*           in the local reordering case above.
*
            DO 2222 DIR = 1, 2
*
*              Broadcast information about the reordering.
*
               DO 321 WINDOW = WINDOW0, WINE, 2
                  RSRC4 = IWORK(IRSRC+WINDOW-1)
                  CSRC4 = IWORK(ICSRC+WINDOW-1)
                  RSRC1 = MOD( RSRC4 - 1 + NPROW, NPROW )
                  CSRC1 = MOD( CSRC4 - 1 + NPCOL, NPCOL )
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 )
     $                  CALL IGEBS2D( ICTXT, 'Row', TOP, 8, 1,
     $                       IBUFF, 8 )
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 )
     $                  CALL IGEBS2D( ICTXT, 'Col', TOP, 8, 1,
     $                       IBUFF, 8 )
                     SKIP1CR = WINDOW.EQ.1 .AND.
     $                    ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                  ELSEIF( MYROW.EQ.RSRC1 .OR. MYCOL.EQ.CSRC1 ) THEN
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 .AND.
     $                    MYROW.EQ.RSRC1 ) THEN
                        CALL IGEBR2D( ICTXT, 'Row', TOP, 8, 1,
     $                       IBUFF, 8, RSRC1, CSRC1 )
                        I = IBUFF( 1 )
                        NWIN = IBUFF( 2 )
                        PITRAF = IBUFF( 3 )
                        KS = IBUFF( 4 )
                        PDTRAF = IBUFF( 5 )
                        NDTRAF = IBUFF( 6 )
                        ILEN = IBUFF( 7 )
                        DLEN = IBUFF( 8 )
                        BUFFLEN = ILEN + DLEN
                        IPW3 = IPW2 + NWIN*NWIN
                        DIM1 = NB - MOD(I-1,NB)
                        DIM4 = NWIN - DIM1
                        LIHIC = NWIN + I - 1
                        SKIP1CR = WINDOW.EQ.1 .AND.
     $                       ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                     END IF
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 .AND.
     $                    MYCOL.EQ.CSRC1 ) THEN
                        CALL IGEBR2D( ICTXT, 'Col', TOP, 8, 1,
     $                       IBUFF, 8, RSRC1, CSRC1 )
                        I = IBUFF( 1 )
                        NWIN = IBUFF( 2 )
                        PITRAF = IBUFF( 3 )
                        KS = IBUFF( 4 )
                        PDTRAF = IBUFF( 5 )
                        NDTRAF = IBUFF( 6 )
                        ILEN = IBUFF( 7 )
                        DLEN = IBUFF( 8 )
                        BUFFLEN = ILEN + DLEN
                        IPW3 = IPW2 + NWIN*NWIN
                        DIM1 = NB - MOD(I-1,NB)
                        DIM4 = NWIN - DIM1
                        LIHIC = NWIN + I - 1
                        SKIP1CR = WINDOW.EQ.1 .AND.
     $                       ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                     END IF
                  END IF
                  IF( RSRC1.NE.RSRC4 ) THEN
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        IF( NPCOL.GT.1 .AND. DIR.EQ.1 )
     $                     CALL IGEBS2D( ICTXT, 'Row', TOP, 8, 1,
     $                          IBUFF, 8 )
                        SKIP1CR = WINDOW.EQ.1 .AND.
     $                       ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                     ELSEIF( MYROW.EQ.RSRC4 ) THEN
                        IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) THEN
                           CALL IGEBR2D( ICTXT, 'Row', TOP, 8, 1,
     $                          IBUFF, 8, RSRC4, CSRC4 )
                           I = IBUFF( 1 )
                           NWIN = IBUFF( 2 )
                           PITRAF = IBUFF( 3 )
                           KS = IBUFF( 4 )
                           PDTRAF = IBUFF( 5 )
                           NDTRAF = IBUFF( 6 )
                           ILEN = IBUFF( 7 )
                           DLEN = IBUFF( 8 )
                           BUFFLEN = ILEN + DLEN
                           IPW3 = IPW2 + NWIN*NWIN
                           DIM1 = NB - MOD(I-1,NB)
                           DIM4 = NWIN - DIM1
                           LIHIC = NWIN + I - 1
                           SKIP1CR = WINDOW.EQ.1 .AND.
     $                          ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                        END IF
                     END IF
                  END IF
                  IF( CSRC1.NE.CSRC4 ) THEN
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        IF( NPROW.GT.1 .AND. DIR.EQ.2 )
     $                     CALL IGEBS2D( ICTXT, 'Col', TOP, 8, 1,
     $                          IBUFF, 8 )
                        SKIP1CR = WINDOW.EQ.1 .AND.
     $                       ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                     ELSEIF( MYCOL.EQ.CSRC4 ) THEN
                        IF( NPROW.GT.1 .AND. DIR.EQ.2 ) THEN
                           CALL IGEBR2D( ICTXT, 'Col', TOP, 8, 1,
     $                          IBUFF, 8, RSRC4, CSRC4 )
                           I = IBUFF( 1 )
                           NWIN = IBUFF( 2 )
                           PITRAF = IBUFF( 3 )
                           KS = IBUFF( 4 )
                           PDTRAF = IBUFF( 5 )
                           NDTRAF = IBUFF( 6 )
                           ILEN = IBUFF( 7 )
                           DLEN = IBUFF( 8 )
                           BUFFLEN = ILEN + DLEN
                           IPW3 = IPW2 + NWIN*NWIN
                           DIM1 = NB - MOD(I-1,NB)
                           DIM4 = NWIN - DIM1
                           LIHIC = NWIN + I - 1
                           SKIP1CR = WINDOW.EQ.1 .AND.
     $                          ICEIL(LIHIC,NB).LE.ICEIL(ILO,NB)
                        END IF
                     END IF
                  END IF
*
*                 Skip rest of broadcasts and updates if appropriate.
*
                  IF( SKIP1CR ) GO TO 326
*
*                 Broadcast the orthogonal transformations.
*
                  IF( MYROW.EQ.RSRC1 .AND. MYCOL.EQ.CSRC1 ) THEN
                     BUFFER = PDTRAF
                     BUFFLEN = DLEN + ILEN
                     IF( (NPROW.GT.1 .AND. DIR.EQ.2) .OR.
     $                   (NPCOL.GT.1 .AND. DIR.EQ.1) ) THEN
                        DO 370 INDX = 1, ILEN
                           WORK( BUFFER+INDX-1 ) =
     $                          DBLE( IWORK(IPIW+INDX-1) )
 370                    CONTINUE
                        CALL DLAMOV( 'All', DLEN, 1, WORK( IPW3 ),
     $                       DLEN, WORK(BUFFER+ILEN), DLEN )
                     END IF
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) THEN
                        CALL DGEBS2D( ICTXT, 'Row', TOP, BUFFLEN, 1,
     $                       WORK(BUFFER), BUFFLEN )
                     END IF
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 ) THEN
                        CALL DGEBS2D( ICTXT, 'Col', TOP, BUFFLEN, 1,
     $                       WORK(BUFFER), BUFFLEN )
                     END IF
                  ELSEIF( MYROW.EQ.RSRC1 .OR. MYCOL.EQ.CSRC1 ) THEN
                     IF( NPCOL.GT.1 .AND. DIR.EQ.1 .AND.
     $                    MYROW.EQ.RSRC1 ) THEN
                        BUFFER = PDTRAF
                        BUFFLEN = DLEN + ILEN
                        CALL DGEBR2D( ICTXT, 'Row', TOP, BUFFLEN, 1,
     $                       WORK(BUFFER), BUFFLEN, RSRC1, CSRC1 )
                     END IF
                     IF( NPROW.GT.1 .AND. DIR.EQ.2 .AND.
     $                    MYCOL.EQ.CSRC1 ) THEN
                        BUFFER = PDTRAF
                        BUFFLEN = DLEN + ILEN
                        CALL DGEBR2D( ICTXT, 'Col', TOP, BUFFLEN, 1,
     $                       WORK(BUFFER), BUFFLEN, RSRC1, CSRC1 )
                     END IF
                     IF( (NPCOL.GT.1.AND.DIR.EQ.1.AND.MYROW.EQ.RSRC1)
     $                    .OR. (NPROW.GT.1.AND.DIR.EQ.2.AND.
     $                    MYCOL.EQ.CSRC1) ) THEN
                        DO 380 INDX = 1, ILEN
                           IWORK(IPIW+INDX-1) =
     $                          INT( WORK( BUFFER+INDX-1 ) )
 380                    CONTINUE
                        CALL DLAMOV( 'All', DLEN, 1,
     $                       WORK( BUFFER+ILEN ), DLEN,
     $                       WORK( IPW3 ), DLEN )
                     END IF
                  END IF
                  IF( RSRC1.NE.RSRC4 ) THEN
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        BUFFER = PDTRAF
                        BUFFLEN = DLEN + ILEN
                        IF( NPCOL.GT.1 .AND. DIR.EQ.1 ) THEN
                           DO 390 INDX = 1, ILEN
                              WORK( BUFFER+INDX-1 ) =
     $                             DBLE( IWORK(IPIW+INDX-1) )
 390                       CONTINUE
                           CALL DLAMOV( 'All', DLEN, 1, WORK( IPW3 ),
     $                          DLEN, WORK(BUFFER+ILEN), DLEN )
                           CALL DGEBS2D( ICTXT, 'Row', TOP, BUFFLEN,
     $                          1, WORK(BUFFER), BUFFLEN )
                        END IF
                     ELSEIF( MYROW.EQ.RSRC4 .AND. DIR.EQ.1 .AND.
     $                    NPCOL.GT.1 ) THEN
                        BUFFER = PDTRAF
                        BUFFLEN = DLEN + ILEN
                        CALL DGEBR2D( ICTXT, 'Row', TOP, BUFFLEN,
     $                       1, WORK(BUFFER), BUFFLEN, RSRC4, CSRC4 )
                        DO 400 INDX = 1, ILEN
                           IWORK(IPIW+INDX-1) =
     $                          INT( WORK( BUFFER+INDX-1 ) )
 400                    CONTINUE
                        CALL DLAMOV( 'All', DLEN, 1,
     $                       WORK( BUFFER+ILEN ), DLEN,
     $                       WORK( IPW3 ), DLEN )
                     END IF
                  END IF
                  IF( CSRC1.NE.CSRC4 ) THEN
                     IF( MYROW.EQ.RSRC4 .AND. MYCOL.EQ.CSRC4 ) THEN
                        BUFFER = PDTRAF
                        BUFFLEN = DLEN + ILEN
                        IF( NPROW.GT.1 .AND. DIR.EQ.2 ) THEN
                           DO 395 INDX = 1, ILEN
                              WORK( BUFFER+INDX-1 ) =
     $                             DBLE( IWORK(IPIW+INDX-1) )
 395                       CONTINUE
                           CALL DLAMOV( 'All', DLEN, 1, WORK( IPW3 ),
     $                          DLEN, WORK(BUFFER+ILEN), DLEN )
                           CALL DGEBS2D( ICTXT, 'Col', TOP, BUFFLEN,
     $                          1, WORK(BUFFER), BUFFLEN )
                        END IF
                     ELSEIF( MYCOL.EQ.CSRC4 .AND. DIR.EQ.2 .AND.
     $                    NPROW.GT.1 ) THEN
                        BUFFER = PDTRAF
                        BUFFLEN = DLEN + ILEN
                        CALL DGEBR2D( ICTXT, 'Col', TOP, BUFFLEN, 1,
     $                       WORK(BUFFER), BUFFLEN, RSRC4, CSRC4 )
                        DO 402 INDX = 1, ILEN
                           IWORK(IPIW+INDX-1) =
     $                          INT( WORK( BUFFER+INDX-1 ) )
 402                    CONTINUE
                        CALL DLAMOV( 'All', DLEN, 1,
     $                       WORK( BUFFER+ILEN ), DLEN,
     $                       WORK( IPW3 ), DLEN )
                     END IF
                  END IF
*
 326              CONTINUE
*
 321           CONTINUE
*
*              Compute crossborder updates.
*
               DO 322 WINDOW = WINDOW0, WINE, 2
                  IF( WINDOW.EQ.1 .AND. SKIP1CR ) GO TO 327
                  RSRC4 = IWORK(IRSRC+WINDOW-1)
                  CSRC4 = IWORK(ICSRC+WINDOW-1)
                  RSRC1 = MOD( RSRC4 - 1 + NPROW, NPROW )
                  CSRC1 = MOD( CSRC4 - 1 + NPCOL, NPCOL )
*
*                 Prepare workspaces for updates:
*                   IPW3 holds now the orthogonal transformations
*                   IPW4 holds the explicit orthogonal matrix, if formed
*                   IPW5 holds the crossborder block column of T
*                   IPW6 holds the crossborder block row of T
*                   IPW7 holds the crossborder block column of Q
*                        (if WANTQ=.TRUE.)
*                   IPW8 points to the leftover workspace used as lhs in
*                        matrix multiplications
*
                  IF( ((MYCOL.EQ.CSRC1.OR.MYCOL.EQ.CSRC4).AND.DIR.EQ.2)
     $                 .OR. ((MYROW.EQ.RSRC1.OR.MYROW.EQ.RSRC4).AND.
     $                 DIR.EQ.1)) THEN
                     IPW4 = BUFFER
                     IF( DIR.EQ.2 ) THEN
                        IF( WANTQ ) THEN
                           QROWS = NUMROC( N, NB, MYROW, DESCQ( RSRC_ ),
     $                          NPROW )
                        ELSE
                           QROWS = 0
                        END IF
                        TROWS = NUMROC( I-1, NB, MYROW, DESCT( RSRC_ ),
     $                       NPROW )
                     ELSE
                        QROWS = 0
                        TROWS = 0
                     END IF
                     IF( DIR.EQ.1 ) THEN
                        TCOLS = NUMROC( N - (I+DIM1-1), NB, MYCOL,
     $                       CSRC4, NPCOL )
                        IF( MYCOL.EQ.CSRC4 ) TCOLS = TCOLS - DIM4
                     ELSE
                        TCOLS = 0
                     END IF
                     IPW5 = IPW4 + NWIN*NWIN
                     IPW6 = IPW5 + TROWS * NWIN
                     IF( WANTQ ) THEN
                        IPW7 = IPW6 + NWIN * TCOLS
                        IPW8 = IPW7 + QROWS * NWIN
                     ELSE
                        IPW8 = IPW6 + NWIN * TCOLS
                     END IF
                  END IF
*
*                 Let each process row and column involved in the updates
*                 exchange data in T and Q with their neighbours.
*
                  IF( DIR.EQ.2 ) THEN
                     IF( MYCOL.EQ.CSRC1 .OR. MYCOL.EQ.CSRC4 ) THEN
                        DO 410 INDX = 1, NPROW
                           IF( MYCOL.EQ.CSRC1 ) THEN
                              CALL INFOG2L( 1+(INDX-1)*NB, I, DESCT,
     $                             NPROW, NPCOL, MYROW, MYCOL, ILOC,
     $                             JLOC1, RSRC, CSRC1 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLAMOV( 'All', TROWS, DIM1,
     $                                T((JLOC1-1)*LLDT+ILOC), LLDT,
     $                                WORK(IPW5), TROWS )
                                 IF( NPCOL.GT.1 ) THEN
                                    EAST = MOD( MYCOL + 1, NPCOL )
                                    CALL DGESD2D( ICTXT, TROWS, DIM1,
     $                                   WORK(IPW5), TROWS, RSRC,
     $                                   EAST )
                                    CALL DGERV2D( ICTXT, TROWS, DIM4,
     $                                   WORK(IPW5+TROWS*DIM1), TROWS,
     $                                   RSRC, EAST )
                                 END IF
                              END IF
                           END IF
                           IF( MYCOL.EQ.CSRC4 ) THEN
                              CALL INFOG2L( 1+(INDX-1)*NB, I+DIM1,
     $                             DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $                             ILOC, JLOC4, RSRC, CSRC4 )
                              IF( MYROW.EQ.RSRC ) THEN
                                 CALL DLAMOV( 'All', TROWS, DIM4,
     $                                T((JLOC4-1)*LLDT+ILOC), LLDT,
     $                                WORK(IPW5+TROWS*DIM1), TROWS )
                                 IF( NPCOL.GT.1 ) THEN
                                    WEST = MOD( MYCOL-1+NPCOL, NPCOL )
                                    CALL DGESD2D( ICTXT, TROWS, DIM4,
     $                                   WORK(IPW5+TROWS*DIM1), TROWS,
     $                                   RSRC, WEST )
                                    CALL DGERV2D( ICTXT, TROWS, DIM1,
     $                                   WORK(IPW5), TROWS, RSRC,
     $                                   WEST )
                                 END IF
                              END IF
                           END IF
 410                    CONTINUE
                     END IF
                  END IF
*
                  IF( DIR.EQ.1 ) THEN
                     IF( MYROW.EQ.RSRC1 .OR. MYROW.EQ.RSRC4 ) THEN
                        DO 420 INDX = 1, NPCOL
                           IF( MYROW.EQ.RSRC1 ) THEN
                              IF( INDX.EQ.1 ) THEN
                                 CALL INFOG2L( I, LIHIC+1, DESCT, NPROW,
     $                                NPCOL, MYROW, MYCOL, ILOC1, JLOC,
     $                                RSRC1, CSRC )
                              ELSE
                                 CALL INFOG2L( I,
     $                                (ICEIL(LIHIC,NB)+(INDX-2))*NB+1,
     $                                DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $                                ILOC1, JLOC, RSRC1, CSRC )
                              END IF
                              IF( MYCOL.EQ.CSRC ) THEN
                                 CALL DLAMOV( 'All', DIM1, TCOLS,
     $                                T((JLOC-1)*LLDT+ILOC1), LLDT,
     $                                WORK(IPW6), NWIN )
                                 IF( NPROW.GT.1 ) THEN
                                    SOUTH = MOD( MYROW + 1, NPROW )
                                    CALL DGESD2D( ICTXT, DIM1, TCOLS,
     $                                   WORK(IPW6), NWIN, SOUTH,
     $                                   CSRC )
                                    CALL DGERV2D( ICTXT, DIM4, TCOLS,
     $                                   WORK(IPW6+DIM1), NWIN, SOUTH,
     $                                   CSRC )
                                 END IF
                              END IF
                           END IF
                           IF( MYROW.EQ.RSRC4 ) THEN
                              IF( INDX.EQ.1 ) THEN
                                 CALL INFOG2L( I+DIM1, LIHIC+1, DESCT,
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC4,
     $                                JLOC, RSRC4, CSRC )
                              ELSE
                                 CALL INFOG2L( I+DIM1,
     $                                (ICEIL(LIHIC,NB)+(INDX-2))*NB+1,
     $                                DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $                                ILOC4, JLOC, RSRC4, CSRC )
                              END IF
                              IF( MYCOL.EQ.CSRC ) THEN
                                 CALL DLAMOV( 'All', DIM4, TCOLS,
     $                                T((JLOC-1)*LLDT+ILOC4), LLDT,
     $                                WORK(IPW6+DIM1), NWIN )
                                 IF( NPROW.GT.1 ) THEN
                                    NORTH = MOD( MYROW-1+NPROW, NPROW )
                                    CALL DGESD2D( ICTXT, DIM4, TCOLS,
     $                                   WORK(IPW6+DIM1), NWIN, NORTH,
     $                                   CSRC )
                                    CALL DGERV2D( ICTXT, DIM1, TCOLS,
     $                                   WORK(IPW6), NWIN, NORTH,
     $                                   CSRC )
                                 END IF
                              END IF
                           END IF
 420                    CONTINUE
                     END IF
                  END IF
*
                  IF( DIR.EQ.2 ) THEN
                     IF( WANTQ ) THEN
                        IF( MYCOL.EQ.CSRC1 .OR. MYCOL.EQ.CSRC4 ) THEN
                           DO 430 INDX = 1, NPROW
                              IF( MYCOL.EQ.CSRC1 ) THEN
                                 CALL INFOG2L( 1+(INDX-1)*NB, I, DESCQ,
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC,
     $                                JLOC1, RSRC, CSRC1 )
                                 IF( MYROW.EQ.RSRC ) THEN
                                    CALL DLAMOV( 'All', QROWS, DIM1,
     $                                   Q((JLOC1-1)*LLDQ+ILOC), LLDQ,
     $                                   WORK(IPW7), QROWS )
                                    IF( NPCOL.GT.1 ) THEN
                                       EAST = MOD( MYCOL + 1, NPCOL )
                                       CALL DGESD2D( ICTXT, QROWS, DIM1,
     $                                      WORK(IPW7), QROWS, RSRC,
     $                                      EAST )
                                       CALL DGERV2D( ICTXT, QROWS, DIM4,
     $                                      WORK(IPW7+QROWS*DIM1),
     $                                      QROWS, RSRC, EAST )
                                    END IF
                                 END IF
                              END IF
                              IF( MYCOL.EQ.CSRC4 ) THEN
                                 CALL INFOG2L( 1+(INDX-1)*NB, I+DIM1,
     $                                DESCQ, NPROW, NPCOL, MYROW, MYCOL,
     $                                ILOC, JLOC4, RSRC, CSRC4 )
                                 IF( MYROW.EQ.RSRC ) THEN
                                    CALL DLAMOV( 'All', QROWS, DIM4,
     $                                   Q((JLOC4-1)*LLDQ+ILOC), LLDQ,
     $                                   WORK(IPW7+QROWS*DIM1), QROWS )
                                    IF( NPCOL.GT.1 ) THEN
                                       WEST = MOD( MYCOL-1+NPCOL,
     $                                      NPCOL )
                                       CALL DGESD2D( ICTXT, QROWS, DIM4,
     $                                      WORK(IPW7+QROWS*DIM1),
     $                                      QROWS, RSRC, WEST )
                                       CALL DGERV2D( ICTXT, QROWS, DIM1,
     $                                      WORK(IPW7), QROWS, RSRC,
     $                                      WEST )
                                    END IF
                                 END IF
                              END IF
 430                       CONTINUE
                        END IF
                     END IF
                  END IF
*
 327              CONTINUE
*
 322           CONTINUE
*
               DO 323 WINDOW = WINDOW0, WINE, 2
                  RSRC4 = IWORK(IRSRC+WINDOW-1)
                  CSRC4 = IWORK(ICSRC+WINDOW-1)
                  RSRC1 = MOD( RSRC4 - 1 + NPROW, NPROW )
                  CSRC1 = MOD( CSRC4 - 1 + NPCOL, NPCOL )
                  FLOPS = 0
                  IF( ((MYCOL.EQ.CSRC1.OR.MYCOL.EQ.CSRC4).AND.DIR.EQ.2)
     $                 .OR. ((MYROW.EQ.RSRC1.OR.MYROW.EQ.RSRC4).AND.
     $                 DIR.EQ.1) ) THEN
*
*                    Skip this part of the updates if appropriate.
*
                     IF( WINDOW.EQ.1 .AND. SKIP1CR ) GO TO 328
*
*                    Count number of operations to decide whether to use
*                    matrix-matrix multiplications for updating
*                    off-diagonal parts or not.
*
                     NITRAF = PITRAF - IPIW
                     ISHH = .FALSE.
                     DO 405 K = 1, NITRAF
                        IF( IWORK( IPIW + K - 1 ).LE.NWIN ) THEN
                           FLOPS = FLOPS + 6
                        ELSE
                           FLOPS = FLOPS + 11
                           ISHH = .TRUE.
                        END IF
 405                 CONTINUE
*
*                    Perform updates in parallel.
*
                     IF( FLOPS.NE.0 .AND.
     $                    ( 2*FLOPS*100 )/( 2*NWIN*NWIN ) .GE. MMULT )
     $                    THEN
*
                        CALL DLASET( 'All', NWIN, NWIN, ZERO, ONE,
     $                       WORK( IPW4 ), NWIN )
                        WORK(IPW8) = DBLE(MYROW)
                        WORK(IPW8+1) = DBLE(MYCOL)
                        CALL BDLAAPP( 1, NWIN, NWIN, NCB, WORK( IPW4 ),
     $                       NWIN, NITRAF, IWORK(IPIW), WORK( IPW3 ),
     $                       WORK(IPW8) )
*
*                       Test if sparsity structure of orthogonal matrix
*                       can be exploited (see below).
*
                        IF( ISHH .OR. DIM1.NE.KS .OR. DIM4.NE.KS ) THEN
*
*                          Update the columns of T and Q affected by the
*                          reordering.
*
                           IF( DIR.EQ.2 ) THEN
                              DO 440 INDX = 1, MIN(I-1,1+(NPROW-1)*NB),
     $                             NB
                                 IF( MYCOL.EQ.CSRC1 ) THEN
                                    CALL INFOG2L( INDX, I, DESCT, NPROW,
     $                                   NPCOL, MYROW, MYCOL, ILOC,
     $                                   JLOC, RSRC, CSRC1 )
                                    IF( MYROW.EQ.RSRC ) THEN
                                       CALL DGEMM( 'No transpose',
     $                                      'No transpose', TROWS, DIM1,
     $                                      NWIN, ONE, WORK( IPW5 ),
     $                                      TROWS, WORK( IPW4 ), NWIN,
     $                                      ZERO, WORK(IPW8), TROWS )
                                       CALL DLAMOV( 'All', TROWS, DIM1,
     $                                      WORK(IPW8), TROWS,
     $                                      T((JLOC-1)*LLDT+ILOC),
     $                                      LLDT )
                                    END IF
                                 END IF
                                 IF( MYCOL.EQ.CSRC4 ) THEN
                                    CALL INFOG2L( INDX, I+DIM1, DESCT,
     $                                   NPROW, NPCOL, MYROW, MYCOL,
     $                                   ILOC, JLOC, RSRC, CSRC4 )
                                    IF( MYROW.EQ.RSRC ) THEN
                                       CALL DGEMM( 'No transpose',
     $                                      'No transpose', TROWS, DIM4,
     $                                      NWIN, ONE, WORK( IPW5 ),
     $                                      TROWS,
     $                                      WORK( IPW4+NWIN*DIM1 ),
     $                                      NWIN, ZERO, WORK(IPW8),
     $                                      TROWS )
                                       CALL DLAMOV( 'All', TROWS, DIM4,
     $                                      WORK(IPW8), TROWS,
     $                                      T((JLOC-1)*LLDT+ILOC),
     $                                      LLDT )
                                    END IF
                                 END IF
 440                          CONTINUE
*
                              IF( WANTQ ) THEN
                                 DO 450 INDX = 1, MIN(N,1+(NPROW-1)*NB),
     $                                NB
                                    IF( MYCOL.EQ.CSRC1 ) THEN
                                       CALL INFOG2L( INDX, I, DESCQ,
     $                                      NPROW, NPCOL, MYROW, MYCOL,
     $                                      ILOC, JLOC, RSRC, CSRC1 )
                                       IF( MYROW.EQ.RSRC ) THEN
                                          CALL DGEMM( 'No transpose',
     $                                         'No transpose', QROWS,
     $                                         DIM1, NWIN, ONE,
     $                                         WORK( IPW7 ), QROWS,
     $                                         WORK( IPW4 ), NWIN,
     $                                         ZERO, WORK(IPW8),
     $                                         QROWS )
                                          CALL DLAMOV( 'All', QROWS,
     $                                         DIM1, WORK(IPW8), QROWS,
     $                                         Q((JLOC-1)*LLDQ+ILOC),
     $                                         LLDQ )
                                       END IF
                                    END IF
                                    IF( MYCOL.EQ.CSRC4 ) THEN
                                       CALL INFOG2L( INDX, I+DIM1,
     $                                      DESCQ, NPROW, NPCOL, MYROW,
     $                                      MYCOL, ILOC, JLOC, RSRC,
     $                                      CSRC4 )
                                       IF( MYROW.EQ.RSRC ) THEN
                                          CALL DGEMM( 'No transpose',
     $                                         'No transpose', QROWS,
     $                                         DIM4, NWIN, ONE,
     $                                         WORK( IPW7 ), QROWS,
     $                                         WORK( IPW4+NWIN*DIM1 ),
     $                                         NWIN, ZERO, WORK(IPW8),
     $                                         QROWS )
                                          CALL DLAMOV( 'All', QROWS,
     $                                         DIM4, WORK(IPW8), QROWS,
     $                                         Q((JLOC-1)*LLDQ+ILOC),
     $                                         LLDQ )
                                       END IF
                                    END IF
 450                             CONTINUE
                              END IF
                           END IF
*
*                          Update the rows of T affected by the
*                          reordering.
*
                           IF( DIR.EQ.1 ) THEN
                              IF ( LIHIC.LT.N ) THEN
                                 IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC4
     $                               .AND.MOD(LIHIC,NB).NE.0 ) THEN
                                    INDX = LIHIC + 1
                                    CALL INFOG2L( I, INDX, DESCT, NPROW,
     $                                   NPCOL, MYROW, MYCOL, ILOC,
     $                                   JLOC, RSRC1, CSRC4 )
                                    CALL DGEMM( 'Transpose',
     $                                   'No Transpose', DIM1, TCOLS,
     $                                   NWIN, ONE, WORK(IPW4), NWIN,
     $                                   WORK( IPW6 ), NWIN, ZERO,
     $                                   WORK(IPW8), DIM1 )
                                    CALL DLAMOV( 'All', DIM1, TCOLS,
     $                                   WORK(IPW8), DIM1,
     $                                   T((JLOC-1)*LLDT+ILOC), LLDT )
                                 END IF
                                 IF( MYROW.EQ.RSRC4.AND.MYCOL.EQ.CSRC4
     $                               .AND.MOD(LIHIC,NB).NE.0 ) THEN
                                    INDX = LIHIC + 1
                                    CALL INFOG2L( I+DIM1, INDX, DESCT,
     $                                   NPROW, NPCOL, MYROW, MYCOL,
     $                                   ILOC, JLOC, RSRC4, CSRC4 )
                                    CALL DGEMM( 'Transpose',
     $                                  'No Transpose', DIM4, TCOLS,
     $                                   NWIN, ONE,
     $                                   WORK( IPW4+DIM1*NWIN ), NWIN,
     $                                   WORK( IPW6), NWIN, ZERO,
     $                                   WORK(IPW8), DIM4 )
                                    CALL DLAMOV( 'All', DIM4, TCOLS,
     $                                   WORK(IPW8), DIM4,
     $                                   T((JLOC-1)*LLDT+ILOC), LLDT )
                                 END IF
                                 INDXS = ICEIL(LIHIC,NB)*NB + 1
                                 INDXE = MIN(N,INDXS+(NPCOL-2)*NB)
                                 DO 460 INDX = INDXS, INDXE, NB
                                    IF( MYROW.EQ.RSRC1 ) THEN
                                       CALL INFOG2L( I, INDX, DESCT,
     $                                      NPROW, NPCOL, MYROW, MYCOL,
     $                                      ILOC, JLOC, RSRC1, CSRC )
                                       IF( MYCOL.EQ.CSRC ) THEN
                                          CALL DGEMM( 'Transpose',
     $                                         'No Transpose', DIM1,
     $                                         TCOLS, NWIN, ONE,
     $                                         WORK( IPW4 ), NWIN,
     $                                         WORK( IPW6 ), NWIN,
     $                                         ZERO, WORK(IPW8), DIM1 )
                                          CALL DLAMOV( 'All', DIM1,
     $                                         TCOLS, WORK(IPW8), DIM1,
     $                                         T((JLOC-1)*LLDT+ILOC),
     $                                         LLDT )
                                       END IF
                                    END IF
                                    IF( MYROW.EQ.RSRC4 ) THEN
                                       CALL INFOG2L( I+DIM1, INDX,
     $                                      DESCT, NPROW, NPCOL, MYROW,
     $                                      MYCOL, ILOC, JLOC, RSRC4,
     $                                      CSRC )
                                       IF( MYCOL.EQ.CSRC ) THEN
                                          CALL DGEMM( 'Transpose',
     $                                         'No Transpose', DIM4,
     $                                         TCOLS, NWIN, ONE,
     $                                         WORK( IPW4+NWIN*DIM1 ),
     $                                         NWIN, WORK( IPW6 ),
     $                                         NWIN, ZERO, WORK(IPW8),
     $                                         DIM4 )
                                          CALL DLAMOV( 'All', DIM4,
     $                                         TCOLS, WORK(IPW8), DIM4,
     $                                         T((JLOC-1)*LLDT+ILOC),
     $                                         LLDT )
                                       END IF
                                    END IF
 460                             CONTINUE
                              END IF
                           END IF
                        ELSE
*
*                          The NWIN-by-NWIN matrix U containing the
*                          accumulated orthogonal transformations has
*                          the following structure:
*
*                                        [ U11  U12 ]
*                                    U = [          ],
*                                        [ U21  U22 ]
*
*                          where U21 is KS-by-KS upper triangular and
*                          U12 is (NWIN-KS)-by-(NWIN-KS) lower
*                          triangular. For reordering over the border
*                          the structure is only exploited when the
*                          border cuts the columns of U conformally with
*                          the structure itself. This happens exactly
*                          when all eigenvalues in the subcluster was
*                          moved to the other side of the border and
*                          fits perfectly in their new positions, i.e.,
*                          the reordering stops when the last eigenvalue
*                          to cross the border is reordered to the
*                          position closest to the border. Tested by
*                          checking is KS = DIM1 = DIM4 (see above).
*                          This should hold quite often. But this branch
*                          is entered only if all involved eigenvalues
*                          are real.
*
*                          Update the columns of T and Q affected by the
*                          reordering.
*
*                          Compute T2*U21 + T1*U11 on the left side of
*                          the border.
*
                           IF( DIR.EQ.2 ) THEN
                              INDXE = MIN(I-1,1+(NPROW-1)*NB)
                              DO 470 INDX = 1, INDXE, NB
                                 IF( MYCOL.EQ.CSRC1 ) THEN
                                    CALL INFOG2L( INDX, I, DESCT, NPROW,
     $                                   NPCOL, MYROW, MYCOL, ILOC,
     $                                   JLOC, RSRC, CSRC1 )
                                    IF( MYROW.EQ.RSRC ) THEN
                                       CALL DLAMOV( 'All', TROWS, KS,
     $                                      WORK( IPW5+TROWS*DIM4),
     $                                      TROWS, WORK(IPW8), TROWS )
                                       CALL DTRMM( 'Right', 'Upper',
     $                                      'No transpose',
     $                                      'Non-unit', TROWS, KS,
     $                                      ONE, WORK( IPW4+DIM4 ),
     $                                      NWIN, WORK(IPW8), TROWS )
                                       CALL DGEMM( 'No transpose',
     $                                      'No transpose', TROWS, KS,
     $                                      DIM4, ONE, WORK( IPW5 ),
     $                                      TROWS, WORK( IPW4 ), NWIN,
     $                                      ONE, WORK(IPW8), TROWS )
                                       CALL DLAMOV( 'All', TROWS, KS,
     $                                      WORK(IPW8), TROWS,
     $                                      T((JLOC-1)*LLDT+ILOC),
     $                                      LLDT )
                                    END IF
                                 END IF
*
*                                Compute T1*U12 + T2*U22 on the right
*                                side of the border.
*
                                 IF( MYCOL.EQ.CSRC4 ) THEN
                                    CALL INFOG2L( INDX, I+DIM1, DESCT,
     $                                   NPROW, NPCOL, MYROW, MYCOL,
     $                                   ILOC, JLOC, RSRC, CSRC4 )
                                    IF( MYROW.EQ.RSRC ) THEN
                                       CALL DLAMOV( 'All', TROWS, DIM4,
     $                                      WORK(IPW5), TROWS,
     $                                      WORK( IPW8 ), TROWS )
                                       CALL DTRMM( 'Right', 'Lower',
     $                                      'No transpose',
     $                                      'Non-unit', TROWS, DIM4,
     $                                      ONE, WORK( IPW4+NWIN*KS ),
     $                                      NWIN, WORK( IPW8 ), TROWS )
                                       CALL DGEMM( 'No transpose',
     $                                      'No transpose', TROWS, DIM4,
     $                                      KS, ONE,
     $                                      WORK( IPW5+TROWS*DIM4),
     $                                      TROWS,
     $                                      WORK( IPW4+NWIN*KS+DIM4 ),
     $                                      NWIN, ONE, WORK( IPW8 ),
     $                                      TROWS )
                                       CALL DLAMOV( 'All', TROWS, DIM4,
     $                                      WORK(IPW8), TROWS,
     $                                      T((JLOC-1)*LLDT+ILOC),
     $                                      LLDT )
                                    END IF
                                 END IF
 470                          CONTINUE
                              IF( WANTQ ) THEN
*
*                                Compute Q2*U21 + Q1*U11 on the left
*                                side of border.
*
                                 INDXE = MIN(N,1+(NPROW-1)*NB)
                                 DO 480 INDX = 1, INDXE, NB
                                    IF( MYCOL.EQ.CSRC1 ) THEN
                                       CALL INFOG2L( INDX, I, DESCQ,
     $                                      NPROW, NPCOL, MYROW, MYCOL,
     $                                      ILOC, JLOC, RSRC, CSRC1 )
                                       IF( MYROW.EQ.RSRC ) THEN
                                          CALL DLAMOV( 'All', QROWS, KS,
     $                                         WORK( IPW7+QROWS*DIM4),
     $                                         QROWS, WORK(IPW8),
     $                                         QROWS )
                                          CALL DTRMM( 'Right', 'Upper',
     $                                         'No transpose',
     $                                         'Non-unit', QROWS,
     $                                         KS, ONE,
     $                                         WORK( IPW4+DIM4 ), NWIN,
     $                                         WORK(IPW8), QROWS )
                                          CALL DGEMM( 'No transpose',
     $                                         'No transpose', QROWS,
     $                                         KS, DIM4, ONE,
     $                                         WORK( IPW7 ), QROWS,
     $                                         WORK( IPW4 ), NWIN, ONE,
     $                                         WORK(IPW8), QROWS )
                                          CALL DLAMOV( 'All', QROWS, KS,
     $                                         WORK(IPW8), QROWS,
     $                                         Q((JLOC-1)*LLDQ+ILOC),
     $                                         LLDQ )
                                       END IF
                                    END IF
*
*                                   Compute Q1*U12 + Q2*U22 on the right
*                                   side of border.
*
                                    IF( MYCOL.EQ.CSRC4 ) THEN
                                       CALL INFOG2L( INDX, I+DIM1,
     $                                      DESCQ, NPROW, NPCOL, MYROW,
     $                                      MYCOL, ILOC, JLOC, RSRC,
     $                                      CSRC4 )
                                       IF( MYROW.EQ.RSRC ) THEN
                                          CALL DLAMOV( 'All', QROWS,
     $                                         DIM4, WORK(IPW7), QROWS,
     $                                         WORK( IPW8 ), QROWS )
                                          CALL DTRMM( 'Right', 'Lower',
     $                                         'No transpose',
     $                                         'Non-unit', QROWS,
     $                                         DIM4, ONE,
     $                                         WORK( IPW4+NWIN*KS ),
     $                                         NWIN, WORK( IPW8 ),
     $                                         QROWS )
                                          CALL DGEMM( 'No transpose',
     $                                         'No transpose', QROWS,
     $                                         DIM4, KS, ONE,
     $                                         WORK(IPW7+QROWS*(DIM4)),
     $                                         QROWS,
     $                                         WORK(IPW4+NWIN*KS+DIM4),
     $                                         NWIN, ONE, WORK( IPW8 ),
     $                                         QROWS )
                                          CALL DLAMOV( 'All', QROWS,
     $                                         DIM4, WORK(IPW8), QROWS,
     $                                         Q((JLOC-1)*LLDQ+ILOC),
     $                                         LLDQ )
                                       END IF
                                    END IF
 480                             CONTINUE
                              END IF
                           END IF
*
                           IF( DIR.EQ.1 ) THEN
                              IF ( LIHIC.LT.N ) THEN
*
*                                Compute U21**T*T2 + U11**T*T1 on the
*                                upper side of the border.
*
                                 IF( MYROW.EQ.RSRC1.AND.MYCOL.EQ.CSRC4
     $                               .AND.MOD(LIHIC,NB).NE.0 ) THEN
                                    INDX = LIHIC + 1
                                    CALL INFOG2L( I, INDX, DESCT, NPROW,
     $                                   NPCOL, MYROW, MYCOL, ILOC,
     $                                   JLOC, RSRC1, CSRC4 )
                                    CALL DLAMOV( 'All', KS, TCOLS,
     $                                   WORK( IPW6+DIM4 ), NWIN,
     $                                   WORK(IPW8), KS )
                                    CALL DTRMM( 'Left', 'Upper',
     $                                   'Transpose', 'Non-unit',
     $                                   KS, TCOLS, ONE,
     $                                   WORK( IPW4+DIM4 ), NWIN,
     $                                   WORK(IPW8), KS )
                                    CALL DGEMM( 'Transpose',
     $                                   'No transpose', KS, TCOLS,
     $                                   DIM4, ONE, WORK(IPW4), NWIN,
     $                                   WORK(IPW6), NWIN, ONE,
     $                                   WORK(IPW8), KS )
                                    CALL DLAMOV( 'All', KS, TCOLS,
     $                                   WORK(IPW8), KS,
     $                                   T((JLOC-1)*LLDT+ILOC), LLDT )
                                 END IF
*
*                                Compute U12**T*T1 + U22**T*T2 on the
*                                lower side of the border.
*
                                 IF( MYROW.EQ.RSRC4.AND.MYCOL.EQ.CSRC4
     $                               .AND.MOD(LIHIC,NB).NE.0 ) THEN
                                    INDX = LIHIC + 1
                                    CALL INFOG2L( I+DIM1, INDX, DESCT,
     $                                   NPROW, NPCOL, MYROW, MYCOL,
     $                                   ILOC, JLOC, RSRC4, CSRC4 )
                                    CALL DLAMOV( 'All', DIM4, TCOLS,
     $                                   WORK( IPW6 ), NWIN,
     $                                   WORK( IPW8 ), DIM4 )
                                    CALL DTRMM( 'Left', 'Lower',
     $                                   'Transpose', 'Non-unit',
     $                                   DIM4, TCOLS, ONE,
     $                                   WORK( IPW4+NWIN*KS ), NWIN,
     $                                   WORK( IPW8 ), DIM4 )
                                    CALL DGEMM( 'Transpose',
     $                                   'No Transpose', DIM4, TCOLS,
     $                                   KS, ONE,
     $                                   WORK( IPW4+NWIN*KS+DIM4 ),
     $                                   NWIN, WORK( IPW6+DIM1 ), NWIN,
     $                                   ONE, WORK( IPW8), DIM4 )
                                    CALL DLAMOV( 'All', DIM4, TCOLS,
     $                                   WORK(IPW8), DIM4,
     $                                   T((JLOC-1)*LLDT+ILOC), LLDT )
                                 END IF
*
*                                Compute U21**T*T2 + U11**T*T1 on upper
*                                side on border.
*
                                 INDXS = ICEIL(LIHIC,NB)*NB+1
                                 INDXE = MIN(N,INDXS+(NPCOL-2)*NB)
                                 DO 490 INDX = INDXS, INDXE, NB
                                    IF( MYROW.EQ.RSRC1 ) THEN
                                       CALL INFOG2L( I, INDX, DESCT,
     $                                      NPROW, NPCOL, MYROW, MYCOL,
     $                                      ILOC, JLOC, RSRC1, CSRC )
                                       IF( MYCOL.EQ.CSRC ) THEN
                                          CALL DLAMOV( 'All', KS, TCOLS,
     $                                         WORK( IPW6+DIM4 ), NWIN,
     $                                         WORK(IPW8), KS )
                                          CALL DTRMM( 'Left', 'Upper',
     $                                         'Transpose',
     $                                         'Non-unit', KS,
     $                                         TCOLS, ONE,
     $                                         WORK( IPW4+DIM4 ), NWIN,
     $                                         WORK(IPW8), KS )
                                          CALL DGEMM( 'Transpose',
     $                                         'No transpose', KS,
     $                                         TCOLS, DIM4, ONE,
     $                                         WORK(IPW4), NWIN,
     $                                         WORK(IPW6), NWIN, ONE,
     $                                         WORK(IPW8), KS )
                                          CALL DLAMOV( 'All', KS, TCOLS,
     $                                         WORK(IPW8), KS,
     $                                         T((JLOC-1)*LLDT+ILOC),
     $                                         LLDT )
                                       END IF
                                    END IF
*
*                                   Compute U12**T*T1 + U22**T*T2 on
*                                   lower side of border.
*
                                    IF( MYROW.EQ.RSRC4 ) THEN
                                       CALL INFOG2L( I+DIM1, INDX,
     $                                      DESCT, NPROW, NPCOL, MYROW,
     $                                      MYCOL, ILOC, JLOC, RSRC4,
     $                                      CSRC )
                                       IF( MYCOL.EQ.CSRC ) THEN
                                          CALL DLAMOV( 'All', DIM4,
     $                                         TCOLS, WORK( IPW6 ),
     $                                         NWIN, WORK( IPW8 ),
     $                                         DIM4 )
                                          CALL DTRMM( 'Left', 'Lower',
     $                                         'Transpose',
     $                                         'Non-unit', DIM4,
     $                                         TCOLS, ONE,
     $                                         WORK( IPW4+NWIN*KS ),
     $                                         NWIN, WORK( IPW8 ),
     $                                         DIM4 )
                                          CALL DGEMM( 'Transpose',
     $                                         'No Transpose', DIM4,
     $                                         TCOLS, KS, ONE,
     $                                         WORK(IPW4+NWIN*KS+DIM4),
     $                                         NWIN, WORK( IPW6+DIM1 ),
     $                                         NWIN, ONE, WORK( IPW8),
     $                                         DIM4 )
                                          CALL DLAMOV( 'All', DIM4,
     $                                         TCOLS, WORK(IPW8), DIM4,
     $                                         T((JLOC-1)*LLDT+ILOC),
     $                                         LLDT )
                                       END IF
                                    END IF
 490                             CONTINUE
                              END IF
                           END IF
                        END IF
                     ELSEIF( FLOPS.NE.0 ) THEN
*
*                       Update off-diagonal blocks and Q using the
*                       pipelined elementary transformations. Now we
*                       have a delicate problem: how to do this without
*                       redundant work? For now, we let the processes
*                       involved compute the whole crossborder block
*                       rows and column saving only the part belonging
*                       to the corresponding side of the border. To make
*                       this a realistic alternative, we have modified
*                       the ratio r_flops (see Reference [2] above) to
*                       give more favor to the ordinary matrix
*                       multiplication.
*
                        IF( DIR.EQ.2 ) THEN
                           INDXE =  MIN(I-1,1+(NPROW-1)*NB)
                           DO 500 INDX = 1, INDXE, NB
                              CALL INFOG2L( INDX, I, DESCT, NPROW,
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                             RSRC, CSRC )
                              IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC )
     $                             THEN
                                 CALL BDLAAPP( 1, TROWS, NWIN, NCB,
     $                                WORK(IPW5), TROWS, NITRAF,
     $                                IWORK(IPIW), WORK( IPW3 ),
     $                                WORK(IPW8) )
                                 CALL DLAMOV( 'All', TROWS, DIM1,
     $                                WORK(IPW5), TROWS,
     $                                T((JLOC-1)*LLDT+ILOC ), LLDT )
                              END IF
                              CALL INFOG2L( INDX, I+DIM1, DESCT, NPROW,
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                             RSRC, CSRC )
                              IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC )
     $                             THEN
                                 IF( NPCOL.GT.1 )
     $                                CALL BDLAAPP( 1, TROWS, NWIN, NCB,
     $                                WORK(IPW5), TROWS, NITRAF,
     $                                IWORK(IPIW), WORK( IPW3 ),
     $                                WORK(IPW8) )
                                 CALL DLAMOV( 'All', TROWS, DIM4,
     $                                WORK(IPW5+TROWS*DIM1), TROWS,
     $                                T((JLOC-1)*LLDT+ILOC ), LLDT )
                              END IF
 500                       CONTINUE
                           IF( WANTQ ) THEN
                              INDXE = MIN(N,1+(NPROW-1)*NB)
                              DO 510 INDX = 1, INDXE, NB
                                 CALL INFOG2L( INDX, I, DESCQ, NPROW,
     $                                NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                                RSRC, CSRC )
                                 IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC )
     $                                THEN
                                    CALL BDLAAPP( 1, QROWS, NWIN, NCB,
     $                                   WORK(IPW7), QROWS, NITRAF,
     $                                   IWORK(IPIW), WORK( IPW3 ),
     $                                   WORK(IPW8) )
                                    CALL DLAMOV( 'All', QROWS, DIM1,
     $                                   WORK(IPW7), QROWS,
     $                                   Q((JLOC-1)*LLDQ+ILOC ), LLDQ )
                                 END IF
                                 CALL INFOG2L( INDX, I+DIM1, DESCQ,
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC,
     $                                JLOC, RSRC, CSRC )
                                 IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC )
     $                                THEN
                                    IF( NPCOL.GT.1 )
     $                                   CALL BDLAAPP( 1, QROWS, NWIN,
     $                                   NCB, WORK(IPW7), QROWS,
     $                                   NITRAF, IWORK(IPIW),
     $                                   WORK( IPW3 ), WORK(IPW8) )
                                    CALL DLAMOV( 'All', QROWS, DIM4,
     $                                   WORK(IPW7+QROWS*DIM1), QROWS,
     $                                   Q((JLOC-1)*LLDQ+ILOC ), LLDQ )
                                 END IF
 510                          CONTINUE
                           END IF
                        END IF
*
                        IF( DIR.EQ.1 ) THEN
                           IF( LIHIC.LT.N ) THEN
                              INDX = LIHIC + 1
                              CALL INFOG2L( I, INDX, DESCT, NPROW,
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                             RSRC, CSRC )
                              IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC.AND.
     $                            MOD(LIHIC,NB).NE.0 ) THEN
                                 CALL BDLAAPP( 0, NWIN, TCOLS, NCB,
     $                                WORK( IPW6 ), NWIN, NITRAF,
     $                                IWORK(IPIW), WORK( IPW3 ),
     $                                WORK(IPW8) )
                                 CALL DLAMOV( 'All', DIM1, TCOLS,
     $                                WORK( IPW6 ), NWIN,
     $                                T((JLOC-1)*LLDT+ILOC), LLDT )
                              END IF
                              CALL INFOG2L( I+DIM1, INDX, DESCT, NPROW,
     $                             NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                             RSRC, CSRC )
                              IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC.AND.
     $                             MOD(LIHIC,NB).NE.0 ) THEN
                                 IF( NPROW.GT.1 )
     $                                CALL BDLAAPP( 0, NWIN, TCOLS, NCB,
     $                                WORK( IPW6 ), NWIN, NITRAF,
     $                                IWORK(IPIW), WORK( IPW3 ),
     $                                WORK(IPW8) )
                                 CALL DLAMOV( 'All', DIM4, TCOLS,
     $                                WORK( IPW6+DIM1 ), NWIN,
     $                                T((JLOC-1)*LLDT+ILOC), LLDT )
                              END IF
                              INDXS = ICEIL(LIHIC,NB)*NB + 1
                              INDXE = MIN(N,INDXS+(NPCOL-2)*NB)
                              DO 520 INDX = INDXS, INDXE, NB
                                 CALL INFOG2L( I, INDX, DESCT, NPROW,
     $                                NPCOL, MYROW, MYCOL, ILOC, JLOC,
     $                                RSRC, CSRC )
                                 IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC )
     $                                THEN
                                    CALL BDLAAPP( 0, NWIN, TCOLS, NCB,
     $                                   WORK(IPW6), NWIN, NITRAF,
     $                                   IWORK(IPIW), WORK( IPW3 ),
     $                                   WORK(IPW8) )
                                    CALL DLAMOV( 'All', DIM1, TCOLS,
     $                                   WORK( IPW6 ), NWIN,
     $                                   T((JLOC-1)*LLDT+ILOC), LLDT )
                                 END IF
                                 CALL INFOG2L( I+DIM1, INDX, DESCT,
     $                                NPROW, NPCOL, MYROW, MYCOL, ILOC,
     $                                JLOC, RSRC, CSRC )
                                 IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC )
     $                                THEN
                                    IF( NPROW.GT.1 )
     $                                   CALL BDLAAPP( 0, NWIN, TCOLS,
     $                                   NCB, WORK(IPW6), NWIN, NITRAF,
     $                                   IWORK(IPIW), WORK( IPW3 ),
     $                                   WORK(IPW8) )
                                    CALL DLAMOV( 'All', DIM4, TCOLS,
     $                                   WORK( IPW6+DIM1 ), NWIN,
     $                                   T((JLOC-1)*LLDT+ILOC), LLDT )
                                 END IF
 520                          CONTINUE
                           END IF
                        END IF
                     END IF
                  END IF
*
 328              CONTINUE
*
 323           CONTINUE
*
*              End of loops over directions (DIR).
*
 2222       CONTINUE
*
*           End of loops over diagonal blocks for reordering over the
*           block diagonal.
*
 310     CONTINUE
         LAST = LAST + 1
         IF( LASTWAIT .AND. LAST.LT.2 ) GO TO 308
*
*        Barrier to collect the processes before proceeding.
*
         CALL BLACS_BARRIER( ICTXT, 'All' )
*
*        Compute global maximum of IERR so that we know if some process
*        experienced a failure in the reordering.
*
         MYIERR = IERR
         IF( NPROCS.GT.1 )
     $      CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, IERR, 1, -1,
     $           -1, -1, -1, -1 )
*
         IF( IERR.NE.0 ) THEN
*
*           When calling BDTREXC, the block at position I+KKS-1 failed
*           to swap.
*
            IF( MYIERR.NE.0 ) INFO = MAX(1,I+KKS-1)
            IF( NPROCS.GT.1 )
     $         CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, INFO, 1, -1,
     $              -1, -1, -1, -1 )
            GO TO 300
         END IF
*
*        Do a global update of the SELECT vector.
*
         DO 530 K = 1, N
            RSRC = INDXG2P( K, NB, MYROW, DESCT( RSRC_ ), NPROW )
            CSRC = INDXG2P( K, NB, MYCOL, DESCT( CSRC_ ), NPCOL )
            IF( MYROW.NE.RSRC .OR. MYCOL.NE.CSRC )
     $         SELECT( K ) = 0
 530     CONTINUE
         IF( NPROCS.GT.1 )
     $      CALL IGSUM2D( ICTXT, 'All', TOP, N, 1, SELECT, N, -1, -1 )
*
*        Find the global minumum of ILO and IHI.
*
         ILO = ILO - 1
 523     CONTINUE
         ILO = ILO + 1
         IF( ILO.LE.N ) THEN
            IF( SELECT(ILO).NE.0 ) GO TO 523
         END IF
         IHI = IHI + 1
 527     CONTINUE
         IHI = IHI - 1
         IF( IHI.GE.1 ) THEN
            IF( SELECT(IHI).EQ.0 ) GO TO 527
         END IF
*
*        End While ( ILO <= M )
         GO TO 50
      END IF
*
 300  CONTINUE
*
*     In case an error occured, do an additional global update of
*     SELECT.
*
      IF( INFO.NE.0 ) THEN
         DO 540 K = 1, N
            RSRC = INDXG2P( K, NB, MYROW, DESCT( RSRC_ ), NPROW )
            CSRC = INDXG2P( K, NB, MYCOL, DESCT( CSRC_ ), NPCOL )
            IF( MYROW.NE.RSRC .OR. MYCOL.NE.CSRC )
     $           SELECT( K ) = 0
 540     CONTINUE
         IF( NPROCS.GT.1 )
     $        CALL IGSUM2D( ICTXT, 'All', TOP, N, 1, SELECT, N, -1, -1 )
      END IF
*
 545  CONTINUE
*
*     Store the output eigenvalues in WR and WI: first let all the
*     processes compute the eigenvalue inside their diagonal blocks in
*     parallel, except for the eigenvalue located next to a block
*     border. After that, compute all eigenvalues located next to the
*     block borders. Finally, do a global summation over WR and WI so
*     that all processors receive the result. Notice: real eigenvalues
*     extracted from a non-canonical 2-by-2 block are not stored in
*     any particular order.
*
      DO 550 K = 1, N
         WR( K ) = ZERO
         WI( K ) = ZERO
 550  CONTINUE
*
*     Loop 560: extract eigenvalues from the blocks which are not laid
*     out across a border of the processor mesh, except for those 1x1
*     blocks on the border.
*
      PAIR = .FALSE.
      DO 560 K = 1, N
         IF( .NOT. PAIR ) THEN
            BORDER = ( K.NE.N .AND. MOD( K, NB ).EQ.0 ) .OR.
     %           ( K.NE.1 .AND. MOD( K, NB ).EQ.1 )
            IF( .NOT. BORDER ) THEN
               CALL INFOG2L( K, K, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $              ILOC1, JLOC1, TRSRC1, TCSRC1 )
               IF( MYROW.EQ.TRSRC1 .AND. MYCOL.EQ.TCSRC1 ) THEN
                  ELEM1 = T((JLOC1-1)*LLDT+ILOC1)
                  IF( K.LT.N ) THEN
                     ELEM3 = T((JLOC1-1)*LLDT+ILOC1+1)
                  ELSE
                     ELEM3 = ZERO
                  END IF
                  IF( ELEM3.NE.ZERO ) THEN
                     ELEM2 = T((JLOC1)*LLDT+ILOC1)
                     ELEM4 = T((JLOC1)*LLDT+ILOC1+1)
                     CALL DLANV2( ELEM1, ELEM2, ELEM3, ELEM4,
     $                    WR( K ), WI( K ), WR( K+1 ), WI( K+1 ), SN,
     $                    CS )
                     PAIR = .TRUE.
                  ELSE
                     IF( K.GT.1 ) THEN
                        TMP = T((JLOC1-2)*LLDT+ILOC1)
                        IF( TMP.NE.ZERO ) THEN
                           ELEM1 = T((JLOC1-2)*LLDT+ILOC1-1)
                           ELEM2 = T((JLOC1-1)*LLDT+ILOC1-1)
                           ELEM3 = T((JLOC1-2)*LLDT+ILOC1)
                           ELEM4 = T((JLOC1-1)*LLDT+ILOC1)
                           CALL DLANV2( ELEM1, ELEM2, ELEM3, ELEM4,
     $                          WR( K-1 ), WI( K-1 ), WR( K ),
     $                          WI( K ), SN, CS )
                        ELSE
                           WR( K ) = ELEM1
                        END IF
                     ELSE
                        WR( K ) = ELEM1
                     END IF
                  END IF
               END IF
            END IF
         ELSE
            PAIR = .FALSE.
         END IF
 560  CONTINUE
*
*     Loop 570: extract eigenvalues from the blocks which are laid
*     out across a border of the processor mesh. The processors are
*     numbered as below:
*
*                1 | 2
*                --+--
*                3 | 4
*
      DO 570 K = NB, N-1, NB
         CALL INFOG2L( K, K, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $        ILOC1, JLOC1, TRSRC1, TCSRC1 )
         CALL INFOG2L( K, K+1, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $        ILOC2, JLOC2, TRSRC2, TCSRC2 )
         CALL INFOG2L( K+1, K, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $        ILOC3, JLOC3, TRSRC3, TCSRC3 )
         CALL INFOG2L( K+1, K+1, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $        ILOC4, JLOC4, TRSRC4, TCSRC4 )
         IF( MYROW.EQ.TRSRC2 .AND. MYCOL.EQ.TCSRC2 ) THEN
            ELEM2 = T((JLOC2-1)*LLDT+ILOC2)
            IF( TRSRC1.NE.TRSRC2 .OR. TCSRC1.NE.TCSRC2 )
     $         CALL DGESD2D( ICTXT, 1, 1, ELEM2, 1, TRSRC1, TCSRC1 )
         END IF
         IF( MYROW.EQ.TRSRC3 .AND. MYCOL.EQ.TCSRC3 ) THEN
            ELEM3 = T((JLOC3-1)*LLDT+ILOC3)
            IF( TRSRC1.NE.TRSRC3 .OR. TCSRC1.NE.TCSRC3 )
     $         CALL DGESD2D( ICTXT, 1, 1, ELEM3, 1, TRSRC1, TCSRC1 )
         END IF
         IF( MYROW.EQ.TRSRC4 .AND. MYCOL.EQ.TCSRC4 ) THEN
            WORK(1) = T((JLOC4-1)*LLDT+ILOC4)
            IF( K+1.LT.N ) THEN
               WORK(2) = T((JLOC4-1)*LLDT+ILOC4+1)
            ELSE
               WORK(2) = ZERO
            END IF
            IF( TRSRC1.NE.TRSRC4 .OR. TCSRC1.NE.TCSRC4 )
     $         CALL DGESD2D( ICTXT, 2, 1, WORK, 2, TRSRC1, TCSRC1 )
         END IF
         IF( MYROW.EQ.TRSRC1 .AND. MYCOL.EQ.TCSRC1 ) THEN
            ELEM1 = T((JLOC1-1)*LLDT+ILOC1)
            IF( TRSRC1.NE.TRSRC2 .OR. TCSRC1.NE.TCSRC2 )
     $         CALL DGERV2D( ICTXT, 1, 1, ELEM2, 1, TRSRC2, TCSRC2 )
            IF( TRSRC1.NE.TRSRC3 .OR. TCSRC1.NE.TCSRC3 )
     $         CALL DGERV2D( ICTXT, 1, 1, ELEM3, 1, TRSRC3, TCSRC3 )
            IF( TRSRC1.NE.TRSRC4 .OR. TCSRC1.NE.TCSRC4 )
     $         CALL DGERV2D( ICTXT, 2, 1, WORK, 2, TRSRC4, TCSRC4 )
            ELEM4 = WORK(1)
            ELEM5 = WORK(2)
            IF( ELEM5.EQ.ZERO ) THEN
               IF( WR( K ).EQ.ZERO .AND. WI( K ).EQ.ZERO ) THEN
                  CALL DLANV2( ELEM1, ELEM2, ELEM3, ELEM4, WR( K ),
     $                 WI( K ), WR( K+1 ), WI( K+1 ), SN, CS )
               ELSEIF( WR( K+1 ).EQ.ZERO .AND. WI( K+1 ).EQ.ZERO ) THEN
                  WR( K+1 ) = ELEM4
               END IF
            ELSEIF( WR( K ).EQ.ZERO .AND. WI( K ).EQ.ZERO ) THEN
               WR( K ) = ELEM1
            END IF
         END IF
 570  CONTINUE
*
      IF( NPROCS.GT.1 ) THEN
         CALL DGSUM2D( ICTXT, 'All', TOP, N, 1, WR, N, -1, -1 )
         CALL DGSUM2D( ICTXT, 'All', TOP, N, 1, WI, N, -1, -1 )
      END IF
*
*     Store storage requirements in workspaces.
*
      WORK( 1 ) = DBLE(LWMIN)
      IWORK( 1 ) = LIWMIN
*
*     Return to calling program.
*
      RETURN
*
*     End of PDTRORD
*
      END
*
