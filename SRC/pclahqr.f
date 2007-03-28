      SUBROUTINE PCLAHQR( WANTT, WANTZ, N, ILO, IHI, A, DESCA, W, ILOZ,
     $                    IHIZ, Z, DESCZ, WORK, LWORK, IWORK, ILWORK,
     $                    INFO )
*
*  -- ScaLAPACK routine (version 1.7.3) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     1.7.3: March   22, 2006
*            modification suggested by Mark Fahey and Greg Henry
*     1.7.0: July    31, 2001
*
*     .. Scalar Arguments ..
      LOGICAL            WANTT, WANTZ
      INTEGER            IHI, IHIZ, ILO, ILOZ, ILWORK, INFO, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * ), IWORK( * )
      COMPLEX            A( * ), W( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLAHQR is an auxiliary routine used to find the Schur decomposition
*    and or eigenvalues of a matrix already in Hessenberg form from
*    cols ILO to IHI.
*  If Z = I, and WANTT=WANTZ=.TRUE., H gets replaced with Z'HZ,
*    with Z'Z=I, and H in Schur form.
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
*                                 array.  LLD_A >= MAX(1,LOCp(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCp( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCq( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCp() and LOCq() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCp( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCq( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCp( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCq( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  WANTT   (global input) LOGICAL
*          = .TRUE. : the full Schur form T is required;
*          = .FALSE.: only eigenvalues are required.
*
*  WANTZ   (global input) LOGICAL
*          = .TRUE. : the matrix of Schur vectors Z is required;
*          = .FALSE.: Schur vectors are not required.
*
*  N       (global input) INTEGER
*          The order of the Hessenberg matrix A (and Z if WANTZ).
*          N >= 0.
*
*  ILO     (global input) INTEGER
*  IHI     (global input) INTEGER
*          It is assumed that A is already upper quasi-triangular in
*          rows and columns IHI+1:N, and that A(ILO,ILO-1) = 0 (unless
*          ILO = 1). PCLAHQR works primarily with the Hessenberg
*          submatrix in rows and columns ILO to IHI, but applies
*          transformations to all of H if WANTT is .TRUE..
*          1 <= ILO <= max(1,IHI); IHI <= N.
*
*  A       (global input/output) COMPLEX array, dimension
*          (DESCA(LLD_),*)
*          On entry, the upper Hessenberg matrix A.
*          On exit, if WANTT is .TRUE., A is upper triangular in rows
*          and columns ILO:IHI.  If WANTT is .FALSE., the contents of
*          A are unspecified on exit.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  W      (global replicated output) COMPLEX array, dimension (N)
*          The computed eigenvalues ILO to IHI are stored in the
*          corresponding elements of W.  If WANTT is .TRUE., the
*          eigenvalues are stored in the same order as on the diagonal
*          of the Schur form returned in A.  A may be returned with
*          larger diagonal blocks until the next release.
*
*  ILOZ    (global input) INTEGER
*  IHIZ    (global input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if WANTZ is .TRUE..
*          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
*
*  Z       (global input/output) COMPLEX array.
*          If WANTZ is .TRUE., on entry Z must contain the current
*          matrix Z of transformations accumulated by PCHSEQR, and on
*          exit Z has been updated; transformations are applied only to
*          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
*          If WANTZ is .FALSE., Z is not referenced.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*  WORK    (local output) COMPLEX array of size LWORK
*          (Unless LWORK=-1, in which case WORK must be at least size 1)
*
*  LWORK   (local input) INTEGER
*          WORK(LWORK) is a local array and LWORK is assumed big enough
*          so that LWORK >= 3*N +
*                MAX( 2*MAX(DESCZ(LLD_),DESCA(LLD_)) + 2*LOCq(N),
*                     7*Ceil(N/HBL)/LCM(NPROW,NPCOL)) +
*                MAX( 2*N, (8*LCM(NPROW,NPCOL)+2)**2 )
*          If LWORK=-1, then WORK(1) gets set to the above number and
*          the code returns immediately.
*
*  IWORK   (global and local input) INTEGER array of size ILWORK
*          This will hold some of the IBLK integer arrays.
*          This is held as a place holder for a future release.
*          Currently unreferenced.
*
*  ILWORK  (local input) INTEGER
*          This will hold the size of the IWORK array.
*          This is held as a place holder for a future release.
*          Currently unreferenced.
*
*  INFO    (global output) INTEGER
*          < 0: parameter number -INFO incorrect or inconsistent
*          = 0: successful exit
*          > 0: PCLAHQR failed to compute all the eigenvalues ILO to IHI
*               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
*               elements i+1:ihi of W contains those eigenvalues
*               which have been successfully computed.
*
*  Logic:
*       This algorithm is very similar to SLAHQR.  Unlike SLAHQR,
*       instead of sending one double shift through the largest
*       unreduced submatrix, this algorithm sends multiple double shifts
*       and spaces them apart so that there can be parallelism across
*       several processor row/columns.  Another critical difference is
*       that this algorithm aggregrates multiple transforms together in
*       order to apply them in a block fashion.
*
*  Important Local Variables:
*       IBLK = The maximum number of bulges that can be computed.
*           Currently fixed.  Future releases this won't be fixed.
*       HBL  = The square block size (HBL=DESCA(MB_)=DESCA(NB_))
*       ROTN = The number of transforms to block together
*       NBULGE = The number of bulges that will be attempted on the
*           current submatrix.
*       IBULGE = The current number of bulges started.
*       K1(*),K2(*) = The current bulge loops from K1(*) to K2(*).
*
*  Subroutines:
*       From LAPACK, this routine calls:
*           CLAHQR     -> Serial QR used to determine shifts and
*                         eigenvalues
*           CLARFG     -> Determine the Householder transforms
*
*       This ScaLAPACK, this routine calls:
*           PCLACONSB  -> To determine where to start each iteration
*           CLAMSH     -> Sends multiple shifts through a small
*                         submatrix to see how the consecutive
*                         subdiagonals change (if PCLACONSB indicates
*                         we can start a run in the middle)
*           PCLAWIL    -> Given the shift, get the transformation
*           PCLACP3    -> Parallel array to local replicated array copy
*                         & back.
*           CLAREF     -> Row/column reflector applier.  Core routine
*                         here.
*           PCLASMSUB  -> Finds negligible subdiagonal elements.
*
*  Current Notes and/or Restrictions:
*       1.) This code requires the distributed block size to be square
*           and at least six (6); unlike simpler codes like LU, this
*           algorithm is extremely sensitive to block size.  Unwise
*           choices of too small a block size can lead to bad
*           performance.
*       2.) This code requires A and Z to be distributed identically
*           and have identical contxts.  A future version may allow Z to
*           have a different contxt to 1D row map it to all nodes (so no
*           communication on Z is necessary.)
*       3.) This code does not currently block the initial transforms
*           so that none of the rows or columns for any bulge are
*           completed until all are started.  To offset pipeline
*           start-up it is recommended that at least 2*LCM(NPROW,NPCOL)
*           bulges are used (if possible)
*       4.) The maximum number of bulges currently supported is fixed at
*           32.  In future versions this will be limited only by the
*           incoming WORK and IWORK array.
*       5.) The matrix A must be in upper Hessenberg form.  If elements
*           below the subdiagonal are nonzero, the resulting transforms
*           may be nonsimilar.  This is also true with the LAPACK
*           routine CLAHQR.
*       6.) For this release, this code has only been tested for
*           RSRC_=CSRC_=0, but it has been written for the general case.
*       7.) Currently, all the eigenvalues are distributed to all the
*           nodes.  Future releases will probably distribute the
*           eigenvalues by the column partitioning.
*       8.) The internals of this routine are subject to change.
*       9.) To optimize this for your architecture, try tuning CLAREF.
*       10.) This code has only been tested for WANTZ = .TRUE. and may
*           behave unpredictably for WANTZ set to .FALSE.
*
*  Further Details
*  ===============
*
*  Contributed by Mark Fahey, June, 2000.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               RONE
      PARAMETER          ( RONE = 1.0E+0 )
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ),
     $                   ONE = ( 1.0E+0, 0.0E+0 ) )
      REAL               CONST
      PARAMETER          ( CONST = 1.50E+0 )
      INTEGER            IBLK
      PARAMETER          ( IBLK = 32 )
*     ..
*     .. Local Scalars ..
      LOGICAL            SKIP
      INTEGER            CONTXT, DOWN, HBL, I, I1, I2, IAFIRST, IBULGE,
     $                   ICBUF, ICOL, ICOL1, ICOL2, IDIA, IERR, II,
     $                   IRBUF, IROW, IROW1, IROW2, ISPEC, ISTART,
     $                   ISTARTCOL, ISTARTROW, ISTOP, ISUB, ISUP,
     $                   ITERMAX, ITMP1, ITMP2, ITN, ITS, IZBUF, J,
     $                   JAFIRST, JBLK, JJ, K, KI, L, LCMRC, LDA, LDZ,
     $                   LEFT, LIHIH, LIHIZ, LILOH, LILOZ, LOCALI1,
     $                   LOCALI2, LOCALK, LOCALM, M, MODKM1, MYCOL,
     $                   MYROW, NBULGE, NH, NODE, NPCOL, NPROW, NQ, NR,
     $                   NUM, NZ, RIGHT, ROTN, UP, VECSIDX
      REAL               CS, OVFL, S, SMLNUM, ULP, UNFL
      COMPLEX            CDUM, H10, H11, H22, H33, H43H34, H44, SN, SUM,
     $                   T1, T1COPY, T2, T3, V1SAVE, V2, V2SAVE, V3,
     $                   V3SAVE
*     ..
*     .. Local Arrays ..
      INTEGER            ICURCOL( IBLK ), ICURROW( IBLK ), K1( IBLK ),
     $                   K2( IBLK ), KCOL( IBLK ), KP2COL( IBLK ),
     $                   KP2ROW( IBLK ), KROW( IBLK )
      COMPLEX            S1( 2*IBLK, 2*IBLK ), SMALLA( 6, 6, IBLK ),
     $                   VCOPY( 3 )
*     ..
*     .. External Functions ..
      INTEGER            ILCM, NUMROC
      REAL               PSLAMCH
      EXTERNAL           ILCM, NUMROC, PSLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGAMN2D, IGEBR2D, IGEBS2D,
     $                   INFOG1L, INFOG2L, PSLABAD, PXERBLA, PCLACONSB,
     $                   PCLACP3, PCLASMSUB, PCLAWIL, PCROT, CCOPY,
     $                   CGEBR2D, CGEBS2D, CGERV2D, CGESD2D, CGSUM2D,
     $                   CLAHQR2, CLAMSH, CLANV2, CLAREF, CLARFG
*     ..
*     .. Intrinsic Functions ..
*
      INTRINSIC          ABS, REAL, CONJG, AIMAG, MAX, MIN, MOD
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
      ITERMAX = 30*( IHI-ILO+1 )
      IF( N.EQ.0 )
     $   RETURN
*
*     NODE (IAFIRST,JAFIRST) OWNS A(1,1)
*
      HBL = DESCA( MB_ )
      CONTXT = DESCA( CTXT_ )
      LDA = DESCA( LLD_ )
      IAFIRST = DESCA( RSRC_ )
      JAFIRST = DESCA( CSRC_ )
      LDZ = DESCZ( LLD_ )
      CALL BLACS_GRIDINFO( CONTXT, NPROW, NPCOL, MYROW, MYCOL )
      NODE = MYROW*NPCOL + MYCOL
      NUM = NPROW*NPCOL
      LEFT = MOD( MYCOL+NPCOL-1, NPCOL )
      RIGHT = MOD( MYCOL+1, NPCOL )
      UP = MOD( MYROW+NPROW-1, NPROW )
      DOWN = MOD( MYROW+1, NPROW )
      LCMRC = ILCM( NPROW, NPCOL )
      IF( ( NPROW.LE.3 ) .OR. ( NPCOL.LE.3 ) ) THEN
         SKIP = .TRUE.
      ELSE
         SKIP = .FALSE.
      END IF
*
*     Determine the number of columns we have so we can check workspace
*
      NQ = NUMROC( N, HBL, MYCOL, JAFIRST, NPCOL )
      JJ = N / HBL
      IF( JJ*HBL.LT.N )
     $   JJ = JJ + 1
      JJ = 7*JJ / LCMRC
      JJ = 3*N + MAX( 2*MAX( LDA, LDZ )+2*NQ, JJ )
      JJ = JJ + MAX( 2*N, ( 8*LCMRC+2 )**2 )
      IF( LWORK.EQ.-1 ) THEN
         WORK( 1 ) = JJ
         RETURN
      END IF
      IF( LWORK.LT.JJ ) THEN
         INFO = -14
      END IF
      IF( DESCZ( CTXT_ ).NE.DESCA( CTXT_ ) ) THEN
         INFO = -( 1300+CTXT_ )
      END IF
      IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
         INFO = -( 700+NB_ )
      END IF
      IF( DESCZ( MB_ ).NE.DESCZ( NB_ ) ) THEN
         INFO = -( 1300+NB_ )
      END IF
      IF( DESCA( MB_ ).NE.DESCZ( MB_ ) ) THEN
         INFO = -( 1300+MB_ )
      END IF
      IF( ( DESCA( RSRC_ ).NE.0 ) .OR. ( DESCA( CSRC_ ).NE.0 ) ) THEN
         INFO = -( 700+RSRC_ )
      END IF
      IF( ( DESCZ( RSRC_ ).NE.0 ) .OR. ( DESCZ( CSRC_ ).NE.0 ) ) THEN
         INFO = -( 1300+RSRC_ )
      END IF
      IF( ( ILO.GT.N ) .OR. ( ILO.LT.1 ) ) THEN
         INFO = -4
      END IF
      IF( ( IHI.GT.N ) .OR. ( IHI.LT.1 ) ) THEN
         INFO = -5
      END IF
      IF( HBL.LT.5 ) THEN
         INFO = -( 700+MB_ )
      END IF
      CALL IGAMN2D( CONTXT, 'ALL', ' ', 1, 1, INFO, 1, ITMP1, ITMP2, -1,
     $              -1, -1 )
      IF( INFO.LT.0 ) THEN
         CALL PXERBLA( CONTXT, 'PCLAHQR', -INFO )
         RETURN
      END IF
*
*     Set work array indices
*
      VECSIDX = 0
      IDIA = 3*N
      ISUB = 3*N
      ISUP = 3*N
      IRBUF = 3*N
      ICBUF = 3*N
      IZBUF = 5*N
*
*     Find a value for ROTN
*
      ROTN = HBL / 3
      ROTN = MIN( ROTN, HBL-2 )
      ROTN = MAX( ROTN, 1 )
*
      IF( ILO.EQ.IHI ) THEN
         CALL INFOG2L( ILO, ILO, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 IROW, ICOL, II, JJ )
         IF( ( MYROW.EQ.II ) .AND. ( MYCOL.EQ.JJ ) ) THEN
            W( ILO ) = A( ( ICOL-1 )*LDA+IROW )
         ELSE
            W( ILO ) = ZERO
         END IF
         RETURN
      END IF
*
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
*
      CALL INFOG1L( ILOZ, HBL, NPROW, MYROW, IAFIRST, LILOZ, LIHIZ )
      LIHIZ = NUMROC( IHIZ, HBL, MYROW, IAFIRST, NPROW )
*
*     Set machine-dependent constants for the stopping criterion.
*     If NORM(H) <= SQRT(OVFL), overflow should not occur.
*
      UNFL = PSLAMCH( CONTXT, 'SAFE MINIMUM' )
      OVFL = RONE / UNFL
      CALL PSLABAD( CONTXT, UNFL, OVFL )
      ULP = PSLAMCH( CONTXT, 'PRECISION' )
      SMLNUM = UNFL*( NH / ULP )
*
*     I1 and I2 are the indices of the first row and last column of H
*     to which transformations must be applied. If eigenvalues only are
*     being computed, I1 and I2 are set inside the main loop.
*
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
*
*     ITN is the total number of QR iterations allowed.
*
      ITN = ITERMAX
*
*     The main loop begins here. I is the loop index and decreases from
*     IHI to ILO in steps of our schur block size (<=2*IBLK). Each
*     iteration of the loop works  with the active submatrix in rows
*     and columns L to I.   Eigenvalues I+1 to IHI have already
*     converged. Either L = ILO or the global A(L,L-1) is negligible
*     so that the matrix splits.
*
      I = IHI
   10 CONTINUE
      L = ILO
      IF( I.LT.ILO )
     $   GO TO 570
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 or 2 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      DO 540 ITS = 0, ITN
*
*        Look for a single small subdiagonal element.
*
         CALL PCLASMSUB( A, DESCA, I, L, K, SMLNUM, WORK( IRBUF+1 ),
     $                   LWORK-IRBUF )
         L = K
*
         IF( L.GT.ILO ) THEN
*
*           H(L,L-1) is negligible
*
            CALL INFOG2L( L, L-1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IROW, ICOL, ITMP1, ITMP2 )
            IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
               A( ( ICOL-1 )*LDA+IROW ) = ZERO
            END IF
            WORK( ISUB+L-1 ) = ZERO
         END IF
*
*        Exit from loop if a submatrix of order 1 or 2 has split off.
*
         IF( WANTT ) THEN
*           For Schur form, use 2x2 blocks
            IF( L.GE.I-1 ) THEN
               GO TO 550
            END IF
         ELSE
*           If we don't want the Schur form, use bigger blocks.
            IF( L.GE.I-( 2*IBLK-1 ) ) THEN
               GO TO 550
            END IF
         END IF
*
*        Now the active submatrix is in rows and columns L to I. If
*        eigenvalues only are being computed, only the active submatrix
*        need be transformed.
*
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
*
*        Copy submatrix of size 2*JBLK and prepare to do generalized
*           Wilkinson shift or an exceptional shift
*
         JBLK = MIN( IBLK, ( ( I-L+1 ) / 2 )-1 )
         IF( JBLK.GT.LCMRC ) THEN
*
*           Make sure it's divisible by LCM (we want even workloads!)
*
            JBLK = JBLK - MOD( JBLK, LCMRC )
         END IF
         JBLK = MIN( JBLK, 2*LCMRC )
         JBLK = MAX( JBLK, 1 )
*
         CALL PCLACP3( 2*JBLK, I-2*JBLK+1, A, DESCA, S1, 2*IBLK, -1, -1,
     $                 0 )
         IF( ( ITS.EQ.20 .OR. ITS.EQ.40 ) .AND. ( JBLK.GT.1 ) ) THEN
*
*           Exceptional shift.
*
            DO 20 II = 2*JBLK, 2, -1
               S1( II, II ) = CONST*( CABS1( S1( II, II ) )+
     $                        CABS1( S1( II, II-1 ) ) )
               S1( II, II-1 ) = ZERO
               S1( II-1, II ) = ZERO
   20       CONTINUE
            S1( 1, 1 ) = CONST*CABS1( S1( 1, 1 ) )
         ELSE
            CALL CLAHQR2( .FALSE., .FALSE., 2*JBLK, 1, 2*JBLK, S1,
     $                   2*IBLK, WORK( IRBUF+1 ), 1, 2*JBLK, Z, LDZ,
     $                   IERR )
*
*           Prepare to use Wilkinson's double shift
*
            H44 = S1( 2*JBLK, 2*JBLK )
            H33 = S1( 2*JBLK-1, 2*JBLK-1 )
            H43H34 = S1( 2*JBLK-1, 2*JBLK )*S1( 2*JBLK, 2*JBLK-1 )
*
         END IF
*
*        Look for two consecutive small subdiagonal elements:
*           PCLACONSB is the routine that does this.
*
         CALL PCLACONSB( A, DESCA, I, L, M, H44, H33, H43H34,
     $                   WORK( IRBUF+1 ), LWORK-IRBUF )
*
*        Double-shift QR step
*
*        NBULGE is the number of bulges that will be attempted
*
         ISTOP = MIN( M+ROTN-1-MOD( M-( M / HBL )*HBL-1, ROTN ), I-2 )
         ISTOP = MIN( ISTOP, M+HBL-3-MOD( M-1, HBL ) )
         ISTOP = MIN( ISTOP, I2-2 )
         ISTOP = MAX( ISTOP, M )
         NBULGE = ( I-1-ISTOP ) / HBL
*
*        Do not exceed maximum determined.
*
         NBULGE = MIN( NBULGE, JBLK )
         IF( NBULGE.GT.LCMRC ) THEN
*
*           Make sure it's divisible by LCM (we want even workloads!)
*
            NBULGE = NBULGE - MOD( NBULGE, LCMRC )
         END IF
         NBULGE = MAX( NBULGE, 1 )
*
*        If we are starting in the middle because of consecutive small
*           subdiagonal elements, we need to see how many bulges we
*           can send through without breaking the consecutive small
*           subdiagonal property.
*
         IF( ( NBULGE.GT.1 ) .AND. ( M.GT.L ) ) THEN
*
*           Copy a chunk of elements from global A(M-1:,M-1:)
*
            CALL INFOG2L( M+2, M+2, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IROW1, ICOL1, ITMP1, ITMP2 )
            II = MIN( 4*NBULGE+2, N-M+2 )
            CALL PCLACP3( II, M-1, A, DESCA, WORK( IRBUF+1 ), II, ITMP1,
     $                    ITMP2, 0 )
            IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
*
*              Find a new NBULGE based on the bulges we have.
*
               CALL CLAMSH( S1, 2*IBLK, NBULGE, JBLK, WORK( IRBUF+1 ),
     $                      II, II, ULP )
               IF( NUM.GT.1 ) THEN
                  CALL IGEBS2D( CONTXT, 'ALL', ' ', 1, 1, NBULGE, 1 )
               END IF
            ELSE
*
*              Everyone needs to receive the new NBULGE
*
               CALL IGEBR2D( CONTXT, 'ALL', ' ', 1, 1, NBULGE, 1, ITMP1,
     $                       ITMP2 )
            END IF
         END IF
*
*        IBULGE is the number of bulges going so far
*
         IBULGE = 1
*
*        "A" row defs : main row transforms from LOCALK to LOCALI2
*
         CALL INFOG1L( M, HBL, NPCOL, MYCOL, JAFIRST, ITMP1, LOCALK )
         LOCALK = NQ
         CALL INFOG1L( 1, HBL, NPCOL, MYCOL, JAFIRST, ICOL1, LOCALI2 )
         LOCALI2 = NUMROC( I2, HBL, MYCOL, JAFIRST, NPCOL )
*
*        "A" col defs : main col transforms from LOCALI1 to LOCALM
*
         CALL INFOG1L( I1, HBL, NPROW, MYROW, IAFIRST, LOCALI1, ICOL1 )
         CALL INFOG1L( 1, HBL, NPROW, MYROW, IAFIRST, LOCALM, ICOL1 )
         ICOL1 = NUMROC( MIN( M+3, I ), HBL, MYROW, IAFIRST, NPROW )
*
*        Which row & column will start the bulges
*
         ISTARTROW = MOD( ( M+1 ) / HBL, NPROW ) + IAFIRST
         ISTARTCOL = MOD( ( M+1 ) / HBL, NPCOL ) + JAFIRST
*
         CALL INFOG1L( M, HBL, NPROW, MYROW, IAFIRST, II, ITMP2 )
         CALL INFOG1L( M, HBL, NPCOL, MYCOL, JAFIRST, JJ, ITMP2 )
         CALL INFOG1L( 1, HBL, NPROW, MYROW, IAFIRST, ISTOP,
     $                 KP2ROW( 1 ) )
         KP2ROW( 1 ) = NUMROC( M+2, HBL, MYROW, IAFIRST, NPROW )
         CALL INFOG1L( 1, HBL, NPCOL, MYCOL, JAFIRST, ISTOP,
     $                 KP2COL( 1 ) )
         KP2COL( 1 ) = NUMROC( M+2, HBL, MYCOL, JAFIRST, NPCOL )
*
*        Set all values for bulges.  All bulges are stored in
*          intermediate steps as loops over KI.  Their current "task"
*          over the global M to I-1 values is always K1(KI) to K2(KI).
*          However, because there are many bulges, K1(KI) & K2(KI) might
*          go past that range while later bulges (KI+1,KI+2,etc..) are
*          finishing up.  Even if ROTN=1, in order to minimize border
*          communication sometimes K1(KI)=HBL-2 & K2(KI)=HBL-1 so both
*          border messages can be handled at once.
*
*        Rules:
*              If MOD(K1(KI)-1,HBL) < HBL-2 then MOD(K2(KI)-1,HBL)<HBL-2
*              If MOD(K1(KI)-1,HBL) = HBL-1 then MOD(K2(KI)-1,HBL)=HBL-1
*              K2(KI)-K1(KI) <= ROTN
*
*        We first hit a border when MOD(K1(KI)-1,HBL)=HBL-2 and we hit
*        it again when MOD(K1(KI)-1,HBL)=HBL-1.
*
         DO 30 KI = 1, NBULGE
            K1( KI ) = M
            ISTOP = MIN( M+ROTN-1-MOD( M-( M / HBL )*HBL-1, ROTN ),
     $              I-2 )
            ISTOP = MIN( ISTOP, M+HBL-3-MOD( M-1, HBL ) )
            ISTOP = MIN( ISTOP, I2-2 )
            ISTOP = MAX( ISTOP, M )
            IF( ( MOD( M-1, HBL ).EQ.HBL-2 ) .AND.
     $          ( ISTOP.LT.MIN( I-2, I2-2 ) ) ) THEN
               ISTOP = ISTOP + 1
            END IF
            K2( KI ) = ISTOP
            ICURROW( KI ) = ISTARTROW
            ICURCOL( KI ) = ISTARTCOL
            KROW( KI ) = II
            KCOL( KI ) = JJ
            IF( KI.GT.1 )
     $         KP2ROW( KI ) = KP2ROW( 1 )
            IF( KI.GT.1 )
     $         KP2COL( KI ) = KP2COL( 1 )
   30    CONTINUE
*
*        Get first transform on node who owns M+2,M+2
*
         DO 31 ITMP1 = 1, 3
            VCOPY(ITMP1) = ZERO
   31    CONTINUE
         ITMP1 = ISTARTROW
         ITMP2 = ISTARTCOL
         CALL PCLAWIL( ITMP1, ITMP2, M, A, DESCA, H44, H33, H43H34,
     $                 VCOPY )
         V1SAVE = VCOPY( 1 )
         V2SAVE = VCOPY( 2 )
         V3SAVE = VCOPY( 3 )
*
*        The main implicit shift Francis loops over the bulges starts
*           here!
*
         IF( K2( IBULGE ).LE.I-1 ) THEN
   40       CONTINUE
            IF( ( K1( IBULGE ).GE.M+5 ) .AND. ( IBULGE.LT.NBULGE ) )
     $           THEN
               IF( ( MOD( K2( IBULGE )+2, HBL ).EQ.MOD( K2( IBULGE+1 )+
     $             2, HBL ) ) .AND. ( K1( 1 ).LE.I-1 ) ) THEN
                  H44 = S1( 2*JBLK-2*IBULGE, 2*JBLK-2*IBULGE )
                  H33 = S1( 2*JBLK-2*IBULGE-1, 2*JBLK-2*IBULGE-1 )
                  H43H34 = S1( 2*JBLK-2*IBULGE-1, 2*JBLK-2*IBULGE )*
     $                     S1( 2*JBLK-2*IBULGE, 2*JBLK-2*IBULGE-1 )
                  ITMP1 = ISTARTROW
                  ITMP2 = ISTARTCOL
                  CALL PCLAWIL( ITMP1, ITMP2, M, A, DESCA, H44, H33,
     $                          H43H34, VCOPY )
                  V1SAVE = VCOPY( 1 )
                  V2SAVE = VCOPY( 2 )
                  V3SAVE = VCOPY( 3 )
                  IBULGE = IBULGE + 1
               END IF
            END IF
*
*        When we hit a border, there are row and column transforms that
*          overlap over several processors and the code gets very
*          "congested."  As a remedy, when we first hit a border, a 6x6
*          *local* matrix is generated on one node (called SMALLA) and
*          work is done on that.  At the end of the border, the data is
*          passed back and everything stays a lot simpler.
*
            DO 120 KI = 1, IBULGE
*
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
               K = ISTART
               MODKM1 = MOD( K-1, HBL )
               IF( ( MODKM1.GE.HBL-2 ) .AND. ( K.LE.I-1 ) ) THEN
                  DO 81 ITMP1 = 1, 6
                     DO 82 ITMP2 = 1, 6
                        SMALLA(ITMP1, ITMP2, KI) = ZERO
   82                CONTINUE
   81             CONTINUE
                  IF( ( MODKM1.EQ.HBL-2 ) .AND. ( K.LT.I-1 ) ) THEN
*
*                 Copy 6 elements from global A(K-1:K+4,K-1:K+4)
*
                     ITMP1 = ICURROW( KI )
                     ITMP2 = ICURCOL( KI )
                     CALL PCLACP3( MIN( 6, N-K+2 ), K-1, A, DESCA,
     $                             SMALLA( 1, 1, KI ), 6, ITMP1, ITMP2,
     $                             0 )
                  END IF
                  IF( MODKM1.EQ.HBL-1 ) THEN
*
*                 Copy 6 elements from global A(K-2:K+3,K-2:K+3)
*
                     CALL INFOG2L( K+1, K+1, DESCA, NPROW, NPCOL, MYROW,
     $                             MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                     CALL PCLACP3( MIN( 6, N-K+3 ), K-2, A, DESCA,
     $                             SMALLA( 1, 1, KI ), 6, ITMP1, ITMP2,
     $                             0 )
                  END IF
               END IF
*
*
*           CLAHQR used to have a single row application and a single
*              column application to H.  Here we do something a little
*              more clever.  We break each transformation down into 3
*              parts:
*                  1.) The minimum amount of work it takes to determine
*                        a group of ROTN transformations (this is on
*                        the critical path.) (Loops 50-120)
*                  (the data is broadcast now: loops 180-240)
*                  2.) The small work it takes so that each of the rows
*                        and columns is at the same place.  For example,
*                        all ROTN row transforms are all complete
*                        through some column TMP.  (Loops 250-260)
*                  3.) The majority of the row and column transforms
*                        are then applied in a block fashion.
*                        (row transforms are in loops 280-380)
*                        (col transforms are in loops 400-540)
*
*           Each of these three parts are further subdivided into 3
*           parts:
*               A.) Work at the start of a border when
*                       MOD(ISTART-1,HBL) = HBL-2
*               B.) Work at the end of a border when
*                       MOD(ISTART-1,HBL) = HBL-1
*               C.) Work in the middle of the block when
*                       MOD(ISTART-1,HBL) < HBL-2
*
*           Further optimization is met with the boolean SKIP.  A border
*              communication can be broken into several parts for
*              efficient parallelism:
*                 Loop over all the bulges, just sending the data out
*                 Loop over all the bulges, just doing the work
*                 Loop over all the bulges, just sending the data back.
*
*
               IF( ( MYROW.EQ.ICURROW( KI ) ) .AND.
     $             ( MYCOL.EQ.ICURCOL( KI ) ) .AND.
     $             ( MODKM1.EQ.HBL-2 ) .AND.
     $             ( ISTART.LT.MIN( I-1, ISTOP+1 ) ) ) THEN
                  K = ISTART
                  NR = MIN( 3, I-K+1 )
                  IF( K.GT.M ) THEN
                     CALL CCOPY( NR, SMALLA( 2, 1, KI ), 1, VCOPY, 1 )
                  ELSE
                     VCOPY( 1 ) = V1SAVE
                     VCOPY( 2 ) = V2SAVE
                     VCOPY( 3 ) = V3SAVE
                  END IF
                  CALL CLARFG( NR, VCOPY( 1 ), VCOPY( 2 ), 1, T1COPY )
                  IF( K.GT.M ) THEN
                     SMALLA( 2, 1, KI ) = VCOPY( 1 )
                     SMALLA( 3, 1, KI ) = ZERO
                     IF( K.LT.I-1 )
     $                  SMALLA( 4, 1, KI ) = ZERO
                  ELSE IF( M.GT.L ) THEN
*
*                 Following differs in comparison to pslahqr.
*
                     SMALLA( 2, 1, KI ) = SMALLA( 2, 1, KI ) -
     $                                    CONJG( T1COPY )*
     $                                    SMALLA( 2, 1, KI )
                  END IF
                  V2 = VCOPY( 2 )
                  T2 = T1COPY*V2
                  WORK( VECSIDX+( K-1 )*3+1 ) = VCOPY( 2 )
                  WORK( VECSIDX+( K-1 )*3+2 ) = VCOPY( 3 )
                  WORK( VECSIDX+( K-1 )*3+3 ) = T1COPY
                  IF( NR.EQ.3 ) THEN
*
*                    Do some work so next step is ready...
*
                     T1 = T1COPY
                     V3 = VCOPY( 3 )
                     T3 = T1*V3
                     ITMP1 = MIN( 6, I2+2-K )
                     ITMP2 = MAX( I1-K+2, 1 )
                     DO 50 J = 2, ITMP1
                        SUM = CONJG( T1 )*SMALLA( 2, J, KI ) +
     $                        CONJG( T2 )*SMALLA( 3, J, KI ) +
     $                        CONJG( T3 )*SMALLA( 4, J, KI )
                        SMALLA( 2, J, KI ) = SMALLA( 2, J, KI ) - SUM
                        SMALLA( 3, J, KI ) = SMALLA( 3, J, KI ) - SUM*V2
                        SMALLA( 4, J, KI ) = SMALLA( 4, J, KI ) - SUM*V3
   50                CONTINUE
                     DO 60 J = ITMP2, 5
                        SUM = T1*SMALLA( J, 2, KI ) +
     $                        T2*SMALLA( J, 3, KI ) +
     $                        T3*SMALLA( J, 4, KI )
                        SMALLA( J, 2, KI ) = SMALLA( J, 2, KI ) - SUM
                        SMALLA( J, 3, KI ) = SMALLA( J, 3, KI ) -
     $                                       SUM*CONJG( V2 )
                        SMALLA( J, 4, KI ) = SMALLA( J, 4, KI ) -
     $                                       SUM*CONJG( V3 )
   60                CONTINUE
                  END IF
               END IF
*
               IF( ( MOD( ISTOP-1, HBL ).EQ.HBL-1 ) .AND.
     $             ( MYROW.EQ.ICURROW( KI ) ) .AND.
     $             ( MYCOL.EQ.ICURCOL( KI ) ) .AND.
     $             ( ISTART.LE.MIN( I, ISTOP ) ) ) THEN
                  K = ISTOP
                  NR = MIN( 3, I-K+1 )
                  IF( K.GT.M ) THEN
                     CALL CCOPY( NR, SMALLA( 3, 2, KI ), 1, VCOPY, 1 )
                  ELSE
                     VCOPY( 1 ) = V1SAVE
                     VCOPY( 2 ) = V2SAVE
                     VCOPY( 3 ) = V3SAVE
                  END IF
                  CALL CLARFG( NR, VCOPY( 1 ), VCOPY( 2 ), 1, T1COPY )
                  IF( K.GT.M ) THEN
                     SMALLA( 3, 2, KI ) = VCOPY( 1 )
                     SMALLA( 4, 2, KI ) = ZERO
                     IF( K.LT.I-1 )
     $                  SMALLA( 5, 2, KI ) = ZERO
*
*                    Set a subdiagonal to zero now if it's possible
*
                     IF( ( K-2.GT.M ) .AND. ( MOD( K-1, HBL ).GT.1 ) )
     $                    THEN
                        H11 = SMALLA( 1, 1, KI )
                        H10 = SMALLA( 2, 1, KI )
                        H22 = SMALLA( 2, 2, KI )
                        S = CABS1( H11 ) + CABS1( H22 )
                        IF( CABS1( H10 ).LE.MAX( ULP*S, SMLNUM ) ) THEN
                           SMALLA( 2, 1, KI ) = ZERO
                        END IF
                     END IF
                  ELSE IF( M.GT.L ) THEN
*
*                 Following differs in comparison to pslahqr.
*
                     SMALLA( 3, 2, KI ) = SMALLA( 3, 2, KI ) -
     $                                    CONJG( T1COPY )*
     $                                    SMALLA( 3, 2, KI )
                  END IF
                  V2 = VCOPY( 2 )
                  T2 = T1COPY*V2
                  WORK( VECSIDX+( K-1 )*3+1 ) = VCOPY( 2 )
                  WORK( VECSIDX+( K-1 )*3+2 ) = VCOPY( 3 )
                  WORK( VECSIDX+( K-1 )*3+3 ) = T1COPY
                  IF( NR.EQ.3 ) THEN
*
*                    Do some work so next step is ready...
*
                     T1 = T1COPY
                     V3 = VCOPY( 3 )
                     T3 = T1*V3
                     ITMP1 = MIN( 6, I2-K+3 )
                     ITMP2 = MAX( I1-K+3, 1 )
                     DO 70 J = 3, ITMP1
                        SUM = CONJG( T1 )*SMALLA( 3, J, KI ) +
     $                        CONJG( T2 )*SMALLA( 4, J, KI ) +
     $                        CONJG( T3 )*SMALLA( 5, J, KI )
                        SMALLA( 3, J, KI ) = SMALLA( 3, J, KI ) - SUM
                        SMALLA( 4, J, KI ) = SMALLA( 4, J, KI ) - SUM*V2
                        SMALLA( 5, J, KI ) = SMALLA( 5, J, KI ) - SUM*V3
   70                CONTINUE
                     DO 80 J = ITMP2, 6
                        SUM = T1*SMALLA( J, 3, KI ) +
     $                        T2*SMALLA( J, 4, KI ) +
     $                        T3*SMALLA( J, 5, KI )
                        SMALLA( J, 3, KI ) = SMALLA( J, 3, KI ) - SUM
                        SMALLA( J, 4, KI ) = SMALLA( J, 4, KI ) -
     $                                       SUM*CONJG( V2 )
                        SMALLA( J, 5, KI ) = SMALLA( J, 5, KI ) -
     $                                       SUM*CONJG( V3 )
   80                CONTINUE
                  END IF
               END IF
*
               IF( ( MODKM1.EQ.0 ) .AND. ( ISTART.LE.I-1 ) .AND.
     $             ( MYROW.EQ.ICURROW( KI ) ) .AND.
     $             ( RIGHT.EQ.ICURCOL( KI ) ) ) THEN
*
*              (IROW1,ICOL1) is (I,J)-coordinates of H(ISTART,ISTART)
*
                  IROW1 = KROW( KI )
                  ICOL1 = KCOL( KI )
*
*                 The ELSE part of this IF needs updated VCOPY, this
*                 was not necessary in PSLAHQR.
*
                  IF( ISTART.GT.M ) THEN
                     VCOPY( 1 ) = SMALLA( 4, 3, KI )
                     VCOPY( 2 ) = SMALLA( 5, 3, KI )
                     VCOPY( 3 ) = SMALLA( 6, 3, KI )
                     NR = MIN( 3, I-ISTART+1 )
                     CALL CLARFG( NR, VCOPY( 1 ), VCOPY( 2 ), 1,
     $                            T1COPY )
                     A( ( ICOL1-2 )*LDA+IROW1 ) = VCOPY( 1 )
                     A( ( ICOL1-2 )*LDA+IROW1+1 ) = ZERO
                     IF( ISTART.LT.I-1 ) THEN
                        A( ( ICOL1-2 )*LDA+IROW1+2 ) = ZERO
                     END IF
                  ELSE
*
*                    If NPCOL.NE.1 THEN we need updated VCOPY.
*
                     NR = MIN( 3, I-ISTART+1 )
                     IF( NPCOL.EQ.1 ) THEN
                        VCOPY( 1 ) = V1SAVE
                        VCOPY( 2 ) = V2SAVE
                        VCOPY( 3 ) = V3SAVE
                     ELSE
*
*                    Get updated VCOPY from RIGHT
*
                        CALL CGERV2D( CONTXT, 3, 1, VCOPY, 3, MYROW,
     $                                RIGHT )
                     END IF
                     CALL CLARFG( NR, VCOPY( 1 ), VCOPY( 2 ), 1,
     $                            T1COPY )
                     IF( M.GT.L ) THEN
*
*                    Following differs in comparison to pslahqr.
*
                        A( ( ICOL1-2 )*LDA+IROW1 ) = A( ( ICOL1-2 )*LDA+
     $                     IROW1 )*CONJG( ONE-T1COPY )
                     END IF
                  END IF
               END IF
*
               IF( ( MYROW.EQ.ICURROW( KI ) ) .AND.
     $             ( MYCOL.EQ.ICURCOL( KI ) ) .AND.
     $             ( ( ( MODKM1.EQ.HBL-2 ) .AND. ( ISTART.EQ.I-
     $             1 ) ) .OR. ( ( MODKM1.LT.HBL-2 ) .AND. ( ISTART.LE.I-
     $             1 ) ) ) ) THEN
*
*              (IROW1,ICOL1) is (I,J)-coordinates of H(ISTART,ISTART)
*
                  IROW1 = KROW( KI )
                  ICOL1 = KCOL( KI )
                  DO 110 K = ISTART, ISTOP
*
*                    Create and do these transforms
*
                     NR = MIN( 3, I-K+1 )
                     IF( K.GT.M ) THEN
                        IF( MOD( K-1, HBL ).EQ.0 ) THEN
                           VCOPY( 1 ) = SMALLA( 4, 3, KI )
                           VCOPY( 2 ) = SMALLA( 5, 3, KI )
                           VCOPY( 3 ) = SMALLA( 6, 3, KI )
                        ELSE
                           VCOPY( 1 ) = A( ( ICOL1-2 )*LDA+IROW1 )
                           VCOPY( 2 ) = A( ( ICOL1-2 )*LDA+IROW1+1 )
                           IF( NR.EQ.3 ) THEN
                              VCOPY( 3 ) = A( ( ICOL1-2 )*LDA+IROW1+2 )
                           END IF
                        END IF
                     ELSE
                        VCOPY( 1 ) = V1SAVE
                        VCOPY( 2 ) = V2SAVE
                        VCOPY( 3 ) = V3SAVE
                     END IF
*
*                    Must send uptodate copy of VCOPY to left.
*
                     IF( NPCOL.GT.1 .AND. ISTART.LE.M .AND.
     $                   MOD( K-1, HBL ).EQ.0 ) THEN
                        CALL CGESD2D( CONTXT, 3, 1, VCOPY, 3, MYROW,
     $                                LEFT )
                     END IF
                     CALL CLARFG( NR, VCOPY( 1 ), VCOPY( 2 ), 1,
     $                            T1COPY )
                     IF( K.GT.M ) THEN
                        IF( MOD( K-1, HBL ).GT.0 ) THEN
                           A( ( ICOL1-2 )*LDA+IROW1 ) = VCOPY( 1 )
                           A( ( ICOL1-2 )*LDA+IROW1+1 ) = ZERO
                           IF( K.LT.I-1 ) THEN
                              A( ( ICOL1-2 )*LDA+IROW1+2 ) = ZERO
                           END IF
*
*                       Set a subdiagonal to zero now if it's possible
*
                           IF( ( IROW1.GT.2 ) .AND. ( ICOL1.GT.2 ) .AND.
     $                         ( K-2.GT.M ) .AND. ( MOD( K-1,
     $                         HBL ).GT.1 ) ) THEN
                              H11 = A( ( ICOL1-3 )*LDA+IROW1-2 )
                              H10 = A( ( ICOL1-3 )*LDA+IROW1-1 )
                              H22 = A( ( ICOL1-2 )*LDA+IROW1-1 )
                              S = CABS1( H11 ) + CABS1( H22 )
                              IF( CABS1( H10 ).LE.MAX( ULP*S, SMLNUM ) )
     $                             THEN
                                 A( ( ICOL1-3 )*LDA+IROW1-1 ) = ZERO
                              END IF
                           END IF
                        END IF
                     ELSE IF( M.GT.L ) THEN
                        IF( MOD( K-1, HBL ).GT.0 ) THEN
*
*                       Following differs in comparison to pslahqr.
*
                           A( ( ICOL1-2 )*LDA+IROW1 ) = A( ( ICOL1-2 )*
     $                        LDA+IROW1 )*CONJG( ONE-T1COPY )
                        END IF
                     END IF
                     V2 = VCOPY( 2 )
                     T2 = T1COPY*V2
                     WORK( VECSIDX+( K-1 )*3+1 ) = VCOPY( 2 )
                     WORK( VECSIDX+( K-1 )*3+2 ) = VCOPY( 3 )
                     WORK( VECSIDX+( K-1 )*3+3 ) = T1COPY
                     T1 = T1COPY
                     IF( K.LT.ISTOP ) THEN
*
*                       Do some work so next step is ready...
*
                        V3 = VCOPY( 3 )
                        T3 = T1*V3
                        DO 90 J = ( ICOL1-1 )*LDA + IROW1,
     $                          ( MIN( K2( KI )+1, I-1 )+ICOL1-K-1 )*
     $                          LDA + IROW1, LDA
                           SUM = CONJG( T1 )*A( J ) +
     $                           CONJG( T2 )*A( J+1 ) +
     $                           CONJG( T3 )*A( J+2 )
                           A( J ) = A( J ) - SUM
                           A( J+1 ) = A( J+1 ) - SUM*V2
                           A( J+2 ) = A( J+2 ) - SUM*V3
   90                   CONTINUE
                        DO 100 J = IROW1 + 1, IROW1 + 3
                           SUM = T1*A( ( ICOL1-1 )*LDA+J ) +
     $                           T2*A( ICOL1*LDA+J ) +
     $                           T3*A( ( ICOL1+1 )*LDA+J )
                           A( ( ICOL1-1 )*LDA+J ) = A( ( ICOL1-1 )*LDA+
     $                        J ) - SUM
                           A( ICOL1*LDA+J ) = A( ICOL1*LDA+J ) -
     $                                        SUM*CONJG( V2 )
                           A( ( ICOL1+1 )*LDA+J ) = A( ( ICOL1+1 )*LDA+
     $                        J ) - SUM*CONJG( V3 )
  100                   CONTINUE
                     END IF
                     IROW1 = IROW1 + 1
                     ICOL1 = ICOL1 + 1
  110             CONTINUE
               END IF
  120       CONTINUE
*
*           First part of applying the transforms is complete.
*           Broadcasts of the Householder data is done here.
*
            DO 130 KI = 1, IBULGE
*
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
*
*              Broadcast Householder information from the block
*
               IF( ( MYROW.EQ.ICURROW( KI ) ) .AND. ( NPCOL.GT.1 ) .AND.
     $             ( ISTART.LE.ISTOP ) ) THEN
                  IF( MYCOL.NE.ICURCOL( KI ) ) THEN
                     CALL CGEBR2D( CONTXT, 'ROW', ' ',
     $                             3*( ISTOP-ISTART+1 ), 1,
     $                             WORK( VECSIDX+( ISTART-1 )*3+1 ),
     $                             3*( ISTOP-ISTART+1 ), MYROW,
     $                             ICURCOL( KI ) )
                  ELSE
                     CALL CGEBS2D( CONTXT, 'ROW', ' ',
     $                             3*( ISTOP-ISTART+1 ), 1,
     $                             WORK( VECSIDX+( ISTART-1 )*3+1 ),
     $                             3*( ISTOP-ISTART+1 ) )
                  END IF
               END IF
  130       CONTINUE
*
*           Now do column transforms and finish work
*
            DO 140 KI = 1, IBULGE
*
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
*
               IF( ( MYCOL.EQ.ICURCOL( KI ) ) .AND. ( NPROW.GT.1 ) .AND.
     $             ( ISTART.LE.ISTOP ) ) THEN
                  IF( MYROW.NE.ICURROW( KI ) ) THEN
                     CALL CGEBR2D( CONTXT, 'COL', ' ',
     $                             3*( ISTOP-ISTART+1 ), 1,
     $                             WORK( VECSIDX+( ISTART-1 )*3+1 ),
     $                             3*( ISTOP-ISTART+1 ), ICURROW( KI ),
     $                             MYCOL )
                  ELSE
                     CALL CGEBS2D( CONTXT, 'COL', ' ',
     $                             3*( ISTOP-ISTART+1 ), 1,
     $                             WORK( VECSIDX+( ISTART-1 )*3+1 ),
     $                             3*( ISTOP-ISTART+1 ) )
                  END IF
               END IF
  140       CONTINUE
*
*
*           Now do make up work to have things in block fashion
*
            DO 160 KI = 1, IBULGE
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
*
               MODKM1 = MOD( ISTART-1, HBL )
               IF( ( MYROW.EQ.ICURROW( KI ) ) .AND.
     $             ( MYCOL.EQ.ICURCOL( KI ) ) .AND.
     $             ( ( ( MODKM1.EQ.HBL-2 ) .AND. ( ISTART.EQ.I-
     $             1 ) ) .OR. ( ( MODKM1.LT.HBL-2 ) .AND. ( ISTART.LE.I-
     $             1 ) ) ) ) THEN
*
*                 (IROW1,ICOL1) is (I,J)-coordinates of H(ISTART,ISTART)
*
                  IROW1 = KROW( KI )
                  ICOL1 = KCOL( KI )
                  DO 150 K = ISTART, ISTOP
*
*              Catch up on column & border work
*
                     NR = MIN( 3, I-K+1 )
                     V2 = WORK( VECSIDX+( K-1 )*3+1 )
                     V3 = WORK( VECSIDX+( K-1 )*3+2 )
                     T1 = WORK( VECSIDX+( K-1 )*3+3 )
                     T2 = T1*V2
                     IF( K.LT.ISTOP ) THEN
*
*                 Do some work so next step is ready...
*
                        T3 = T1*V3
                        CALL CLAREF( 'Col', A, LDA, .FALSE., Z, LDZ,
     $                               .FALSE., ICOL1, ICOL1, ISTART,
     $                               ISTOP, MIN( ISTART+1, I )-K+IROW1,
     $                               IROW1, LILOZ, LIHIZ,
     $                               WORK( VECSIDX+1 ), V2, V3, T1, T2,
     $                               T3 )
                        IROW1 = IROW1 + 1
                        ICOL1 = ICOL1 + 1
                     ELSE
                        IF( ( NR.EQ.3 ) .AND. ( MOD( K-1,
     $                      HBL ).LT.HBL-2 ) ) THEN
                           T3 = T1*V3
                           CALL CLAREF( 'Row', A, LDA, .FALSE., Z, LDZ,
     $                                  .FALSE., IROW1, IROW1, ISTART,
     $                                  ISTOP, ICOL1, MIN( MIN( K2( KI )
     $                                  +1, I-1 ), I2 )-K+ICOL1, LILOZ,
     $                                  LIHIZ, WORK( VECSIDX+1 ), V2,
     $                                  V3, T1, T2, T3 )
                        END IF
                     END IF
  150             CONTINUE
               END IF
*
*           Send SMALLA back again.
*
               K = ISTART
               MODKM1 = MOD( K-1, HBL )
               IF( ( MODKM1.GE.HBL-2 ) .AND. ( K.LE.I-1 ) ) THEN
                  IF( ( MODKM1.EQ.HBL-2 ) .AND. ( K.LT.I-1 ) ) THEN
*
*                    Copy 6 elements from global A(K-1:K+4,K-1:K+4)
*
                     ITMP1 = ICURROW( KI )
                     ITMP2 = ICURCOL( KI )
                     CALL PCLACP3( MIN( 6, N-K+2 ), K-1, A, DESCA,
     $                             SMALLA( 1, 1, KI ), 6, ITMP1, ITMP2,
     $                             1 )
*
                  END IF
                  IF( MODKM1.EQ.HBL-1 ) THEN
*
*                    Copy 6 elements from global A(K-2:K+3,K-2:K+3)
*
                     ITMP1 = ICURROW( KI )
                     ITMP2 = ICURCOL( KI )
                     CALL PCLACP3( MIN( 6, N-K+3 ), K-2, A, DESCA,
     $                             SMALLA( 1, 1, KI ), 6, ITMP1, ITMP2,
     $                             1 )
                  END IF
               END IF
*
  160       CONTINUE
*
  170       CONTINUE
*
*           Now start major set of block ROW reflections
*
            DO 180 KI = 1, IBULGE
               IF( ( MYROW.NE.ICURROW( KI ) ) .AND.
     $             ( DOWN.NE.ICURROW( KI ) ) )GO TO 180
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
*
               IF( ( ISTOP.GT.ISTART ) .AND.
     $             ( MOD( ISTART-1, HBL ).LT.HBL-2 ) .AND.
     $             ( ICURROW( KI ).EQ.MYROW ) ) THEN
                  IROW1 = MIN( K2( KI )+1, I-1 ) + 1
                  CALL INFOG1L( IROW1, HBL, NPCOL, MYCOL, JAFIRST,
     $                          ITMP1, ITMP2 )
                  ITMP2 = LOCALI2
                  II = KROW( KI )
                  CALL CLAREF( 'Row', A, LDA, WANTZ, Z, LDZ, .TRUE., II,
     $                         II, ISTART, ISTOP, ITMP1, ITMP2, LILOZ,
     $                         LIHIZ, WORK( VECSIDX+1 ), V2, V3, T1, T2,
     $                         T3 )
               END IF
  180       CONTINUE
*
            DO 220 KI = 1, IBULGE
               IF( KROW( KI ).GT.KP2ROW( KI ) )
     $            GO TO 220
               IF( ( MYROW.NE.ICURROW( KI ) ) .AND.
     $             ( DOWN.NE.ICURROW( KI ) ) )GO TO 220
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
               IF( ( ISTART.EQ.ISTOP ) .OR.
     $             ( MOD( ISTART-1, HBL ).GE.HBL-2 ) .OR.
     $             ( ICURROW( KI ).NE.MYROW ) ) THEN
                  DO 210 K = ISTART, ISTOP
                     V2 = WORK( VECSIDX+( K-1 )*3+1 )
                     V3 = WORK( VECSIDX+( K-1 )*3+2 )
                     T1 = WORK( VECSIDX+( K-1 )*3+3 )
                     NR = MIN( 3, I-K+1 )
                     IF( ( NR.EQ.3 ) .AND. ( KROW( KI ).LE.
     $                   KP2ROW( KI ) ) ) THEN
                        IF( ( K.LT.ISTOP ) .AND.
     $                      ( MOD( K-1, HBL ).LT.HBL-2 ) ) THEN
                           ITMP1 = MIN( K2( KI )+1, I-1 ) + 1
                        ELSE
                           IF( MOD( K-1, HBL ).LT.HBL-2 ) THEN
                              ITMP1 = MIN( K2( KI )+1, I-1 ) + 1
                           END IF
                           IF( MOD( K-1, HBL ).EQ.HBL-2 ) THEN
                              ITMP1 = MIN( K+4, I2 ) + 1
                           END IF
                           IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                              ITMP1 = MIN( K+3, I2 ) + 1
                           END IF
                        END IF
*
*                    Find local coor of rows K through K+2
*
                        IROW1 = KROW( KI )
                        IROW2 = KP2ROW( KI )
                        IF( ( K.GT.ISTART ) .AND.
     $                      ( MOD( K-1, HBL ).GE.HBL-2 ) ) THEN
                           IF( DOWN.EQ.ICURROW( KI ) ) THEN
                              IROW1 = IROW1 + 1
                           END IF
                           IF( MYROW.EQ.ICURROW( KI ) ) THEN
                              IROW2 = IROW2 + 1
                           END IF
                        END IF
                        CALL INFOG1L( ITMP1, HBL, NPCOL, MYCOL, JAFIRST,
     $                                ICOL1, ICOL2 )
                        ICOL2 = LOCALI2
                        IF( ( MOD( K-1, HBL ).LT.HBL-2 ) .OR.
     $                      ( NPROW.EQ.1 ) ) THEN
                           T2 = T1*V2
                           T3 = T1*V3
                           CALL CLAREF( 'Row', A, LDA, WANTZ, Z, LDZ,
     $                                  .FALSE., IROW1, IROW1, ISTART,
     $                                  ISTOP, ICOL1, ICOL2, LILOZ,
     $                                  LIHIZ, WORK( VECSIDX+1 ), V2,
     $                                  V3, T1, T2, T3 )
                        END IF
                        IF( ( MOD( K-1, HBL ).EQ.HBL-2 ) .AND.
     $                      ( NPROW.GT.1 ) ) THEN
                           IF( IROW1.NE.IROW2 ) THEN
                              CALL CGESD2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                      A( ( ICOL1-1 )*LDA+IROW1 ),
     $                                      LDA, DOWN, MYCOL )
                              IF( SKIP .AND. ( ISTART.EQ.ISTOP ) ) THEN
                                 CALL CGERV2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                         A( ( ICOL1-1 )*LDA+
     $                                         IROW1 ), LDA, DOWN,
     $                                         MYCOL )
                              END IF
                           ELSE IF( SKIP ) THEN
                              CALL CGERV2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                      WORK( IRBUF+1 ), 2, UP,
     $                                      MYCOL )
                              T2 = T1*V2
                              T3 = T1*V3
                              DO 190 J = ICOL1, ICOL2
                                 SUM = CONJG( T1 )*
     $                                 WORK( IRBUF+2*( J-ICOL1 )+1 ) +
     $                                 CONJG( T2 )*WORK( IRBUF+2*
     $                                 ( J-ICOL1 )+2 ) +
     $                                 CONJG( T3 )*A( ( J-1 )*LDA+
     $                                 IROW1 )
                                 WORK( IRBUF+2*( J-ICOL1 )+1 )
     $                              = WORK( IRBUF+2*( J-ICOL1 )+1 ) -
     $                              SUM
                                 WORK( IRBUF+2*( J-ICOL1 )+2 )
     $                              = WORK( IRBUF+2*( J-ICOL1 )+2 ) -
     $                              SUM*V2
                                 A( ( J-1 )*LDA+IROW1 ) = A( ( J-1 )*
     $                              LDA+IROW1 ) - SUM*V3
  190                         CONTINUE
                              IF( ISTART.EQ.ISTOP ) THEN
                                 CALL CGESD2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                         WORK( IRBUF+1 ), 2, UP,
     $                                         MYCOL )
                              END IF
                           END IF
                        END IF
                        IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                      ( NPROW.GT.1 ) ) THEN
                           IF( IROW1.EQ.IROW2 ) THEN
                              IF( ISTART.EQ.ISTOP ) THEN
                                 CALL CGESD2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                         A( ( ICOL1-1 )*LDA+IROW1-
     $                                         1 ), LDA, DOWN, MYCOL )
                              END IF
                              IF( SKIP ) THEN
                                 CALL CGERV2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                         A( ( ICOL1-1 )*LDA+IROW1-
     $                                         1 ), LDA, DOWN, MYCOL )
                              END IF
                           ELSE IF( SKIP ) THEN
                              IF( ISTART.EQ.ISTOP ) THEN
                                 CALL CGERV2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                         WORK( IRBUF+1 ), 2, UP,
     $                                         MYCOL )
                              END IF
                              T2 = T1*V2
                              T3 = T1*V3
                              DO 200 J = ICOL1, ICOL2
                                 SUM = CONJG( T1 )*
     $                                 WORK( IRBUF+2*( J-ICOL1 )+2 ) +
     $                                 CONJG( T2 )*A( ( J-1 )*LDA+
     $                                 IROW1 ) + CONJG( T3 )*
     $                                 A( ( J-1 )*LDA+IROW1+1 )
                                 WORK( IRBUF+2*( J-ICOL1 )+2 )
     $                              = WORK( IRBUF+2*( J-ICOL1 )+2 ) -
     $                              SUM
                                 A( ( J-1 )*LDA+IROW1 ) = A( ( J-1 )*
     $                              LDA+IROW1 ) - SUM*V2
                                 A( ( J-1 )*LDA+IROW1+1 ) = A( ( J-1 )*
     $                              LDA+IROW1+1 ) - SUM*V3
  200                         CONTINUE
                              CALL CGESD2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                      WORK( IRBUF+1 ), 2, UP,
     $                                      MYCOL )
*
                           END IF
                        END IF
                     END IF
  210             CONTINUE
               END IF
  220       CONTINUE
*
            IF( SKIP )
     $         GO TO 290
*
            DO 260 KI = 1, IBULGE
               IF( KROW( KI ).GT.KP2ROW( KI ) )
     $            GO TO 260
               IF( ( MYROW.NE.ICURROW( KI ) ) .AND.
     $             ( DOWN.NE.ICURROW( KI ) ) )GO TO 260
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
               IF( ( ISTART.EQ.ISTOP ) .OR.
     $             ( MOD( ISTART-1, HBL ).GE.HBL-2 ) .OR.
     $             ( ICURROW( KI ).NE.MYROW ) ) THEN
                  DO 250 K = ISTART, ISTOP
                     V2 = WORK( VECSIDX+( K-1 )*3+1 )
                     V3 = WORK( VECSIDX+( K-1 )*3+2 )
                     T1 = WORK( VECSIDX+( K-1 )*3+3 )
                     NR = MIN( 3, I-K+1 )
                     IF( ( NR.EQ.3 ) .AND. ( KROW( KI ).LE.
     $                   KP2ROW( KI ) ) ) THEN
                        IF( ( K.LT.ISTOP ) .AND.
     $                      ( MOD( K-1, HBL ).LT.HBL-2 ) ) THEN
                           ITMP1 = MIN( K2( KI )+1, I-1 ) + 1
                        ELSE
                           IF( MOD( K-1, HBL ).LT.HBL-2 ) THEN
                              ITMP1 = MIN( K2( KI )+1, I-1 ) + 1
                           END IF
                           IF( MOD( K-1, HBL ).EQ.HBL-2 ) THEN
                              ITMP1 = MIN( K+4, I2 ) + 1
                           END IF
                           IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                              ITMP1 = MIN( K+3, I2 ) + 1
                           END IF
                        END IF
*
*                    Find local coor of rows K through K+2
*
                        IROW1 = KROW( KI )
                        IROW2 = KP2ROW( KI )
                        IF( ( K.GT.ISTART ) .AND.
     $                      ( MOD( K-1, HBL ).GE.HBL-2 ) ) THEN
                           IF( DOWN.EQ.ICURROW( KI ) ) THEN
                              IROW1 = IROW1 + 1
                           END IF
                           IF( MYROW.EQ.ICURROW( KI ) ) THEN
                              IROW2 = IROW2 + 1
                           END IF
                        END IF
                        CALL INFOG1L( ITMP1, HBL, NPCOL, MYCOL, JAFIRST,
     $                                ICOL1, ICOL2 )
                        ICOL2 = LOCALI2
                        IF( ( MOD( K-1, HBL ).EQ.HBL-2 ) .AND.
     $                      ( NPROW.GT.1 ) ) THEN
                           IF( IROW1.EQ.IROW2 ) THEN
                              CALL CGERV2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                      WORK( IRBUF+1 ), 2, UP,
     $                                      MYCOL )
                              T2 = T1*V2
                              T3 = T1*V3
                              DO 230 J = ICOL1, ICOL2
                                 SUM = CONJG( T1 )*
     $                                 WORK( IRBUF+2*( J-ICOL1 )+1 ) +
     $                                 CONJG( T2 )*WORK( IRBUF+2*
     $                                 ( J-ICOL1 )+2 ) +
     $                                 CONJG( T3 )*A( ( J-1 )*LDA+
     $                                 IROW1 )
                                 WORK( IRBUF+2*( J-ICOL1 )+1 )
     $                              = WORK( IRBUF+2*( J-ICOL1 )+1 ) -
     $                              SUM
                                 WORK( IRBUF+2*( J-ICOL1 )+2 )
     $                              = WORK( IRBUF+2*( J-ICOL1 )+2 ) -
     $                              SUM*V2
                                 A( ( J-1 )*LDA+IROW1 ) = A( ( J-1 )*
     $                              LDA+IROW1 ) - SUM*V3
  230                         CONTINUE
                              IF( ISTART.EQ.ISTOP ) THEN
                                 CALL CGESD2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                         WORK( IRBUF+1 ), 2, UP,
     $                                         MYCOL )
                              END IF
                           END IF
                        END IF
                        IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                      ( NPROW.GT.1 ) ) THEN
                           IF( IROW1.NE.IROW2 ) THEN
                              IF( ISTART.EQ.ISTOP ) THEN
                                 CALL CGERV2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                         WORK( IRBUF+1 ), 2, UP,
     $                                         MYCOL )
                              END IF
                              T2 = T1*V2
                              T3 = T1*V3
                              DO 240 J = ICOL1, ICOL2
                                 SUM = CONJG( T1 )*
     $                                 WORK( IRBUF+2*( J-ICOL1 )+2 ) +
     $                                 CONJG( T2 )*A( ( J-1 )*LDA+
     $                                 IROW1 ) + CONJG( T3 )*
     $                                 A( ( J-1 )*LDA+IROW1+1 )
                                 WORK( IRBUF+2*( J-ICOL1 )+2 )
     $                              = WORK( IRBUF+2*( J-ICOL1 )+2 ) -
     $                              SUM
                                 A( ( J-1 )*LDA+IROW1 ) = A( ( J-1 )*
     $                              LDA+IROW1 ) - SUM*V2
                                 A( ( J-1 )*LDA+IROW1+1 ) = A( ( J-1 )*
     $                              LDA+IROW1+1 ) - SUM*V3
  240                         CONTINUE
                              CALL CGESD2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                      WORK( IRBUF+1 ), 2, UP,
     $                                      MYCOL )
                           END IF
                        END IF
                     END IF
  250             CONTINUE
               END IF
  260       CONTINUE
*
            DO 280 KI = 1, IBULGE
               IF( KROW( KI ).GT.KP2ROW( KI ) )
     $            GO TO 280
               IF( ( MYROW.NE.ICURROW( KI ) ) .AND.
     $             ( DOWN.NE.ICURROW( KI ) ) )GO TO 280
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
               IF( ( ISTART.EQ.ISTOP ) .OR.
     $             ( MOD( ISTART-1, HBL ).GE.HBL-2 ) .OR.
     $             ( ICURROW( KI ).NE.MYROW ) ) THEN
                  DO 270 K = ISTART, ISTOP
                     V2 = WORK( VECSIDX+( K-1 )*3+1 )
                     V3 = WORK( VECSIDX+( K-1 )*3+2 )
                     T1 = WORK( VECSIDX+( K-1 )*3+3 )
                     NR = MIN( 3, I-K+1 )
                     IF( ( NR.EQ.3 ) .AND. ( KROW( KI ).LE.
     $                   KP2ROW( KI ) ) ) THEN
                        IF( ( K.LT.ISTOP ) .AND.
     $                      ( MOD( K-1, HBL ).LT.HBL-2 ) ) THEN
                           ITMP1 = MIN( K2( KI )+1, I-1 ) + 1
                        ELSE
                           IF( MOD( K-1, HBL ).LT.HBL-2 ) THEN
                              ITMP1 = MIN( K2( KI )+1, I-1 ) + 1
                           END IF
                           IF( MOD( K-1, HBL ).EQ.HBL-2 ) THEN
                              ITMP1 = MIN( K+4, I2 ) + 1
                           END IF
                           IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                              ITMP1 = MIN( K+3, I2 ) + 1
                           END IF
                        END IF
*
*                    Find local coor of rows K through K+2
*
                        IROW1 = KROW( KI )
                        IROW2 = KP2ROW( KI )
                        IF( ( K.GT.ISTART ) .AND.
     $                      ( MOD( K-1, HBL ).GE.HBL-2 ) ) THEN
                           IF( DOWN.EQ.ICURROW( KI ) ) THEN
                              IROW1 = IROW1 + 1
                           END IF
                           IF( MYROW.EQ.ICURROW( KI ) ) THEN
                              IROW2 = IROW2 + 1
                           END IF
                        END IF
                        CALL INFOG1L( ITMP1, HBL, NPCOL, MYCOL, JAFIRST,
     $                                ICOL1, ICOL2 )
                        ICOL2 = LOCALI2
                        IF( ( MOD( K-1, HBL ).EQ.HBL-2 ) .AND.
     $                      ( NPROW.GT.1 ) ) THEN
                           IF( IROW1.NE.IROW2 ) THEN
                              IF( ISTART.EQ.ISTOP ) THEN
                                 CALL CGERV2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                         A( ( ICOL1-1 )*LDA+
     $                                         IROW1 ), LDA, DOWN,
     $                                         MYCOL )
                              END IF
                           END IF
                        END IF
                        IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                      ( NPROW.GT.1 ) ) THEN
                           IF( IROW1.EQ.IROW2 ) THEN
                              CALL CGERV2D( CONTXT, 2, ICOL2-ICOL1+1,
     $                                      A( ( ICOL1-1 )*LDA+IROW1-
     $                                      1 ), LDA, DOWN, MYCOL )
                           END IF
                        END IF
                     END IF
  270             CONTINUE
               END IF
  280       CONTINUE
*
  290       CONTINUE
*
*           Now start major set of block COL reflections
*
            DO 300 KI = 1, IBULGE
               IF( ( MYCOL.NE.ICURCOL( KI ) ) .AND.
     $             ( RIGHT.NE.ICURCOL( KI ) ) )GO TO 300
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
*
               IF( ( ( MOD( ISTART-1, HBL ).LT.HBL-2 ) .OR. ( NPCOL.EQ.
     $             1 ) ) .AND. ( ICURCOL( KI ).EQ.MYCOL ) .AND.
     $             ( I-ISTOP+1.GE.3 ) ) THEN
                  K = ISTART
                  IF( ( K.LT.ISTOP ) .AND. ( MOD( K-1,
     $                HBL ).LT.HBL-2 ) ) THEN
                     ITMP1 = MIN( ISTART+1, I ) - 1
                  ELSE
                     IF( MOD( K-1, HBL ).LT.HBL-2 ) THEN
                        ITMP1 = MIN( K+3, I )
                     END IF
                     IF( MOD( K-1, HBL ).EQ.HBL-2 ) THEN
                        ITMP1 = MAX( I1, K-1 ) - 1
                     END IF
                     IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                        ITMP1 = MAX( I1, K-2 ) - 1
                     END IF
                  END IF
*
                  ICOL1 = KCOL( KI )
                  CALL INFOG1L( I1, HBL, NPROW, MYROW, IAFIRST, IROW1,
     $                          IROW2 )
                  IROW2 = NUMROC( ITMP1, HBL, MYROW, IAFIRST, NPROW )
                  IF( IROW1.LE.IROW2 ) THEN
                     ITMP2 = IROW2
                  ELSE
                     ITMP2 = -1
                  END IF
                  CALL CLAREF( 'Col', A, LDA, WANTZ, Z, LDZ, .TRUE.,
     $                         ICOL1, ICOL1, ISTART, ISTOP, IROW1,
     $                         IROW2, LILOZ, LIHIZ, WORK( VECSIDX+1 ),
     $                         V2, V3, T1, T2, T3 )
                  K = ISTOP
                  IF( MOD( K-1, HBL ).LT.HBL-2 ) THEN
*
*                 Do from ITMP1+1 to MIN(K+3,I)
*
                     IF( MOD( K-1, HBL ).LT.HBL-3 ) THEN
                        IROW1 = ITMP2 + 1
                        IF( MOD( ( ITMP1 / HBL ), NPROW ).EQ.MYROW )
     $                       THEN
                           IF( ITMP2.GT.0 ) THEN
                              IROW2 = ITMP2 + MIN( K+3, I ) - ITMP1
                           ELSE
                              IROW2 = IROW1 - 1
                           END IF
                        ELSE
                           IROW2 = IROW1 - 1
                        END IF
                     ELSE
                        CALL INFOG1L( ITMP1+1, HBL, NPROW, MYROW,
     $                                IAFIRST, IROW1, IROW2 )
                        IROW2 = NUMROC( MIN( K+3, I ), HBL, MYROW,
     $                          IAFIRST, NPROW )
                     END IF
                     V2 = WORK( VECSIDX+( K-1 )*3+1 )
                     V3 = WORK( VECSIDX+( K-1 )*3+2 )
                     T1 = WORK( VECSIDX+( K-1 )*3+3 )
                     T2 = T1*V2
                     T3 = T1*V3
                     ICOL1 = KCOL( KI ) + ISTOP - ISTART
                     CALL CLAREF( 'Col', A, LDA, .FALSE., Z, LDZ,
     $                            .FALSE., ICOL1, ICOL1, ISTART, ISTOP,
     $                            IROW1, IROW2, LILOZ, LIHIZ,
     $                            WORK( VECSIDX+1 ), V2, V3, T1, T2,
     $                            T3 )
                  END IF
               END IF
  300       CONTINUE
*
            DO 360 KI = 1, IBULGE
               IF( KCOL( KI ).GT.KP2COL( KI ) )
     $            GO TO 360
               IF( ( MYCOL.NE.ICURCOL( KI ) ) .AND.
     $             ( RIGHT.NE.ICURCOL( KI ) ) )GO TO 360
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
               IF( MOD( ISTART-1, HBL ).GE.HBL-2 ) THEN
*
*              INFO is found in a buffer
*
                  ISPEC = 1
               ELSE
*
*              All INFO is local
*
                  ISPEC = 0
               END IF
               DO 350 K = ISTART, ISTOP
*
                  V2 = WORK( VECSIDX+( K-1 )*3+1 )
                  V3 = WORK( VECSIDX+( K-1 )*3+2 )
                  T1 = WORK( VECSIDX+( K-1 )*3+3 )
                  NR = MIN( 3, I-K+1 )
                  IF( ( NR.EQ.3 ) .AND. ( KCOL( KI ).LE.KP2COL( KI ) ) )
     $                 THEN
*
                     IF( ( K.LT.ISTOP ) .AND.
     $                   ( MOD( K-1, HBL ).LT.HBL-2 ) ) THEN
                        ITMP1 = MIN( ISTART+1, I ) - 1
                     ELSE
                        IF( MOD( K-1, HBL ).LT.HBL-2 ) THEN
                           ITMP1 = MIN( K+3, I )
                        END IF
                        IF( MOD( K-1, HBL ).EQ.HBL-2 ) THEN
                           ITMP1 = MAX( I1, K-1 ) - 1
                        END IF
                        IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                           ITMP1 = MAX( I1, K-2 ) - 1
                        END IF
                     END IF
                     IF( MOD( K-1, HBL ).LT.HBL-2 ) THEN
                        ICOL1 = KCOL( KI ) + K - ISTART
                        ICOL2 = KP2COL( KI ) + K - ISTART
                     ELSE
                        ICOL1 = KCOL( KI )
                        ICOL2 = KP2COL( KI )
                        IF( K.GT.ISTART ) THEN
                           IF( RIGHT.EQ.ICURCOL( KI ) ) THEN
                              ICOL1 = ICOL1 + 1
                           END IF
                           IF( MYCOL.EQ.ICURCOL( KI ) ) THEN
                              ICOL2 = ICOL2 + 1
                           END IF
                        END IF
                     END IF
                     CALL INFOG1L( I1, HBL, NPROW, MYROW, IAFIRST,
     $                             IROW1, IROW2 )
                     IROW2 = NUMROC( ITMP1, HBL, MYROW, IAFIRST, NPROW )
                     IF( ( MOD( K-1, HBL ).EQ.HBL-2 ) .AND.
     $                   ( NPCOL.GT.1 ) ) THEN
                        IF( ICOL1.NE.ICOL2 ) THEN
                           CALL CGESD2D( CONTXT, IROW2-IROW1+1, 2,
     $                                   A( ( ICOL1-1 )*LDA+IROW1 ),
     $                                   LDA, MYROW, RIGHT )
                           IF( ( ISTART.EQ.ISTOP ) .AND. SKIP ) THEN
                              CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      A( ( ICOL1-1 )*LDA+IROW1 ),
     $                                      LDA, MYROW, RIGHT )
                           END IF
                        ELSE IF( SKIP ) THEN
                           T2 = T1*V2
                           T3 = T1*V3
                           CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                   WORK( ICBUF+1 ), IROW2-IROW1+1,
     $                                   MYROW, LEFT )
                           II = ICBUF - IROW1 + 1
                           JJ = ICBUF + IROW2 - 2*IROW1 + 2
                           DO 310 J = IROW1, IROW2
                              SUM = T1*WORK( II+J ) + T2*WORK( JJ+J ) +
     $                              T3*A( ( ICOL1-1 )*LDA+J )
                              WORK( II+J ) = WORK( II+J ) - SUM
                              WORK( JJ+J ) = WORK( JJ+J ) -
     $                                       SUM*CONJG( V2 )
                              A( ( ICOL1-1 )*LDA+J ) = A( ( ICOL1-1 )*
     $                           LDA+J ) - SUM*CONJG( V3 )
  310                      CONTINUE
                           IF( ISTART.EQ.ISTOP ) THEN
                              CALL CGESD2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      WORK( ICBUF+1 ),
     $                                      IROW2-IROW1+1, MYROW, LEFT )
                           END IF
                        END IF
                     END IF
                     IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                   ( NPCOL.GT.1 ) ) THEN
                        IF( ICOL1.EQ.ICOL2 ) THEN
                           IF( ISTART.EQ.ISTOP ) THEN
                              CALL CGESD2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      A( ( ICOL1-2 )*LDA+IROW1 ),
     $                                      LDA, MYROW, RIGHT )
                           END IF
                           IF( SKIP ) THEN
                              CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      A( ( ICOL1-2 )*LDA+IROW1 ),
     $                                      LDA, MYROW, RIGHT )
                           END IF
                        ELSE IF( SKIP ) THEN
                           IF( ISTART.EQ.ISTOP ) THEN
                              CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      WORK( ICBUF+1 ),
     $                                      IROW2-IROW1+1, MYROW, LEFT )
                           END IF
                           T2 = T1*V2
                           T3 = T1*V3
                           II = ICBUF + IROW2 - 2*IROW1 + 2
                           DO 320 J = IROW1, IROW2
                              SUM = T1*WORK( J+II ) +
     $                              T2*A( ( ICOL1-1 )*LDA+J ) +
     $                              T3*A( ICOL1*LDA+J )
                              WORK( J+II ) = WORK( J+II ) - SUM
                              A( ( ICOL1-1 )*LDA+J ) = A( ( ICOL1-1 )*
     $                           LDA+J ) - SUM*CONJG( V2 )
                              A( ICOL1*LDA+J ) = A( ICOL1*LDA+J ) -
     $                                           SUM*CONJG( V3 )
  320                      CONTINUE
                           CALL CGESD2D( CONTXT, IROW2-IROW1+1, 2,
     $                                   WORK( ICBUF+1 ), IROW2-IROW1+1,
     $                                   MYROW, LEFT )
                        END IF
                     END IF
*
*                    If we want Z and we haven't already done any Z
*
                     IF( ( WANTZ ) .AND. ( MOD( K-1,
     $                   HBL ).GE.HBL-2 ) .AND. ( NPCOL.GT.1 ) ) THEN
*
*                       Accumulate transformations in the matrix Z
*
                        IROW1 = LILOZ
                        IROW2 = LIHIZ
                        IF( MOD( K-1, HBL ).EQ.HBL-2 ) THEN
                           IF( ICOL1.NE.ICOL2 ) THEN
                              CALL CGESD2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      Z( ( ICOL1-1 )*LDZ+IROW1 ),
     $                                      LDZ, MYROW, RIGHT )
                              IF( ( ISTART.EQ.ISTOP ) .AND. SKIP ) THEN
                                 CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                         Z( ( ICOL1-1 )*LDZ+
     $                                         IROW1 ), LDZ, MYROW,
     $                                         RIGHT )
                              END IF
                           ELSE IF( SKIP ) THEN
                              CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      WORK( IZBUF+1 ),
     $                                      IROW2-IROW1+1, MYROW, LEFT )
                              T2 = T1*V2
                              T3 = T1*V3
                              ICOL1 = ( ICOL1-1 )*LDZ
                              II = IZBUF - IROW1 + 1
                              JJ = IZBUF + IROW2 - 2*IROW1 + 2
                              DO 330 J = IROW1, IROW2
                                 SUM = T1*WORK( II+J ) +
     $                                 T2*WORK( JJ+J ) + T3*Z( ICOL1+J )
                                 WORK( II+J ) = WORK( II+J ) - SUM
                                 WORK( JJ+J ) = WORK( JJ+J ) -
     $                                          SUM*CONJG( V2 )
                                 Z( ICOL1+J ) = Z( ICOL1+J ) -
     $                                          SUM*CONJG( V3 )
  330                         CONTINUE
                              IF( ISTART.EQ.ISTOP ) THEN
                                 CALL CGESD2D( CONTXT, IROW2-IROW1+1, 2,
     $                                         WORK( IZBUF+1 ),
     $                                         IROW2-IROW1+1, MYROW,
     $                                         LEFT )
                              END IF
                           END IF
                        END IF
                        IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                           IF( ICOL1.EQ.ICOL2 ) THEN
                              IF( ISTART.EQ.ISTOP ) THEN
                                 CALL CGESD2D( CONTXT, IROW2-IROW1+1, 2,
     $                                         Z( ( ICOL1-2 )*LDZ+
     $                                         IROW1 ), LDZ, MYROW,
     $                                         RIGHT )
                              END IF
                              IF( SKIP ) THEN
                                 CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                         Z( ( ICOL1-2 )*LDZ+
     $                                         IROW1 ), LDZ, MYROW,
     $                                         RIGHT )
                              END IF
                           ELSE IF( SKIP ) THEN
                              IF( ISTART.EQ.ISTOP ) THEN
                                 CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                         WORK( IZBUF+1 ),
     $                                         IROW2-IROW1+1, MYROW,
     $                                         LEFT )
                              END IF
                              T2 = T1*V2
                              T3 = T1*V3
                              ICOL1 = ( ICOL1-1 )*LDZ
                              II = IZBUF + IROW2 - 2*IROW1 + 2
                              DO 340 J = IROW1, IROW2
                                 SUM = T1*WORK( II+J ) +
     $                                 T2*Z( J+ICOL1 ) +
     $                                 T3*Z( J+ICOL1+LDZ )
                                 WORK( II+J ) = WORK( II+J ) - SUM
                                 Z( J+ICOL1 ) = Z( J+ICOL1 ) -
     $                                          SUM*CONJG( V2 )
                                 Z( J+ICOL1+LDZ ) = Z( J+ICOL1+LDZ ) -
     $                                              SUM*CONJG( V3 )
  340                         CONTINUE
                              CALL CGESD2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      WORK( IZBUF+1 ),
     $                                      IROW2-IROW1+1, MYROW, LEFT )
                           END IF
                        END IF
                     END IF
                  END IF
  350          CONTINUE
  360       CONTINUE
*
            IF( SKIP )
     $         GO TO 450
*
            DO 420 KI = 1, IBULGE
               IF( KCOL( KI ).GT.KP2COL( KI ) )
     $            GO TO 420
               IF( ( MYCOL.NE.ICURCOL( KI ) ) .AND.
     $             ( RIGHT.NE.ICURCOL( KI ) ) )GO TO 420
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
               IF( MOD( ISTART-1, HBL ).GE.HBL-2 ) THEN
*
*                 INFO is found in a buffer
*
                  ISPEC = 1
               ELSE
*
*                 All INFO is local
*
                  ISPEC = 0
               END IF
               DO 410 K = ISTART, ISTOP
*
                  V2 = WORK( VECSIDX+( K-1 )*3+1 )
                  V3 = WORK( VECSIDX+( K-1 )*3+2 )
                  T1 = WORK( VECSIDX+( K-1 )*3+3 )
                  NR = MIN( 3, I-K+1 )
                  IF( ( NR.EQ.3 ) .AND. ( KCOL( KI ).LE.KP2COL( KI ) ) )
     $                 THEN
*
                     IF( ( K.LT.ISTOP ) .AND.
     $                   ( MOD( K-1, HBL ).LT.HBL-2 ) ) THEN
                        ITMP1 = MIN( ISTART+1, I ) - 1
                     ELSE
                        IF( MOD( K-1, HBL ).LT.HBL-2 ) THEN
                           ITMP1 = MIN( K+3, I )
                        END IF
                        IF( MOD( K-1, HBL ).EQ.HBL-2 ) THEN
                           ITMP1 = MAX( I1, K-1 ) - 1
                        END IF
                        IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                           ITMP1 = MAX( I1, K-2 ) - 1
                        END IF
                     END IF
                     IF( MOD( K-1, HBL ).LT.HBL-2 ) THEN
                        ICOL1 = KCOL( KI ) + K - ISTART
                        ICOL2 = KP2COL( KI ) + K - ISTART
                     ELSE
                        ICOL1 = KCOL( KI )
                        ICOL2 = KP2COL( KI )
                        IF( K.GT.ISTART ) THEN
                           IF( RIGHT.EQ.ICURCOL( KI ) ) THEN
                              ICOL1 = ICOL1 + 1
                           END IF
                           IF( MYCOL.EQ.ICURCOL( KI ) ) THEN
                              ICOL2 = ICOL2 + 1
                           END IF
                        END IF
                     END IF
                     CALL INFOG1L( I1, HBL, NPROW, MYROW, IAFIRST,
     $                             IROW1, IROW2 )
                     IROW2 = NUMROC( ITMP1, HBL, MYROW, IAFIRST, NPROW )
                     IF( ( MOD( K-1, HBL ).EQ.HBL-2 ) .AND.
     $                   ( NPCOL.GT.1 ) ) THEN
                        IF( ICOL1.EQ.ICOL2 ) THEN
                           CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                   WORK( ICBUF+1 ), IROW2-IROW1+1,
     $                                   MYROW, LEFT )
                           T2 = T1*V2
                           T3 = T1*V3
                           II = ICBUF - IROW1 + 1
                           JJ = ICBUF + IROW2 - 2*IROW1 + 2
                           DO 370 J = IROW1, IROW2
                              SUM = T1*WORK( II+J ) + T2*WORK( JJ+J ) +
     $                              T3*A( ( ICOL1-1 )*LDA+J )
                              WORK( II+J ) = WORK( II+J ) - SUM
                              WORK( JJ+J ) = WORK( JJ+J ) -
     $                                       SUM*CONJG( V2 )
                              A( ( ICOL1-1 )*LDA+J ) = A( ( ICOL1-1 )*
     $                           LDA+J ) - SUM*CONJG( V3 )
  370                      CONTINUE
                           IF( ISTART.EQ.ISTOP ) THEN
                              CALL CGESD2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      WORK( ICBUF+1 ),
     $                                      IROW2-IROW1+1, MYROW, LEFT )
                           END IF
                        END IF
                     END IF
                     IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                   ( NPCOL.GT.1 ) ) THEN
                        IF( ICOL1.NE.ICOL2 ) THEN
                           IF( ISTART.EQ.ISTOP ) THEN
                              CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      WORK( ICBUF+1 ),
     $                                      IROW2-IROW1+1, MYROW, LEFT )
                           END IF
                           T2 = T1*V2
                           T3 = T1*V3
                           II = ICBUF + IROW2 - 2*IROW1 + 2
                           DO 380 J = IROW1, IROW2
                              SUM = T1*WORK( J+II ) +
     $                              T2*A( ( ICOL1-1 )*LDA+J ) +
     $                              T3*A( ICOL1*LDA+J )
                              WORK( J+II ) = WORK( J+II ) - SUM
                              A( ( ICOL1-1 )*LDA+J ) = A( ( ICOL1-1 )*
     $                           LDA+J ) - SUM*CONJG( V2 )
                              A( ICOL1*LDA+J ) = A( ICOL1*LDA+J ) -
     $                                           SUM*CONJG( V3 )
  380                      CONTINUE
                           CALL CGESD2D( CONTXT, IROW2-IROW1+1, 2,
     $                                   WORK( ICBUF+1 ), IROW2-IROW1+1,
     $                                   MYROW, LEFT )
                        END IF
                     END IF
*
*
*                 If we want Z and we haven't already done any Z
                     IF( ( WANTZ ) .AND. ( MOD( K-1,
     $                   HBL ).GE.HBL-2 ) .AND. ( NPCOL.GT.1 ) ) THEN
*
*                    Accumulate transformations in the matrix Z
*
                        IROW1 = LILOZ
                        IROW2 = LIHIZ
                        IF( MOD( K-1, HBL ).EQ.HBL-2 ) THEN
                           IF( ICOL1.EQ.ICOL2 ) THEN
                              CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      WORK( IZBUF+1 ),
     $                                      IROW2-IROW1+1, MYROW, LEFT )
                              T2 = T1*V2
                              T3 = T1*V3
                              ICOL1 = ( ICOL1-1 )*LDZ
                              II = IZBUF - IROW1 + 1
                              JJ = IZBUF + IROW2 - 2*IROW1 + 2
                              DO 390 J = IROW1, IROW2
                                 SUM = T1*WORK( II+J ) +
     $                                 T2*WORK( JJ+J ) + T3*Z( ICOL1+J )
                                 WORK( II+J ) = WORK( II+J ) - SUM
                                 WORK( JJ+J ) = WORK( JJ+J ) -
     $                                          SUM*CONJG( V2 )
                                 Z( ICOL1+J ) = Z( ICOL1+J ) -
     $                                          SUM*CONJG( V3 )
  390                         CONTINUE
                              IF( ISTART.EQ.ISTOP ) THEN
                                 CALL CGESD2D( CONTXT, IROW2-IROW1+1, 2,
     $                                         WORK( IZBUF+1 ),
     $                                         IROW2-IROW1+1, MYROW,
     $                                         LEFT )
                              END IF
                           END IF
                        END IF
                        IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                           IF( ICOL1.NE.ICOL2 ) THEN
                              IF( ISTART.EQ.ISTOP ) THEN
                                 CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                         WORK( IZBUF+1 ),
     $                                         IROW2-IROW1+1, MYROW,
     $                                         LEFT )
                              END IF
                              T2 = T1*V2
                              T3 = T1*V3
                              ICOL1 = ( ICOL1-1 )*LDZ
                              II = IZBUF + IROW2 - 2*IROW1 + 2
                              DO 400 J = IROW1, IROW2
                                 SUM = T1*WORK( II+J ) +
     $                                 T2*Z( J+ICOL1 ) +
     $                                 T3*Z( J+ICOL1+LDZ )
                                 WORK( II+J ) = WORK( II+J ) - SUM
                                 Z( J+ICOL1 ) = Z( J+ICOL1 ) -
     $                                          SUM*CONJG( V2 )
                                 Z( J+ICOL1+LDZ ) = Z( J+ICOL1+LDZ ) -
     $                                              SUM*CONJG( V3 )
  400                         CONTINUE
                              CALL CGESD2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      WORK( IZBUF+1 ),
     $                                      IROW2-IROW1+1, MYROW, LEFT )
                           END IF
                        END IF
                     END IF
                  END IF
  410          CONTINUE
  420       CONTINUE
*
            DO 440 KI = 1, IBULGE
               IF( KCOL( KI ).GT.KP2COL( KI ) )
     $            GO TO 440
               IF( ( MYCOL.NE.ICURCOL( KI ) ) .AND.
     $             ( RIGHT.NE.ICURCOL( KI ) ) )GO TO 440
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
               IF( MOD( ISTART-1, HBL ).GE.HBL-2 ) THEN
*
*              INFO is found in a buffer
*
                  ISPEC = 1
               ELSE
*
*              All INFO is local
*
                  ISPEC = 0
               END IF
               DO 430 K = ISTART, ISTOP
*
                  V2 = WORK( VECSIDX+( K-1 )*3+1 )
                  V3 = WORK( VECSIDX+( K-1 )*3+2 )
                  T1 = WORK( VECSIDX+( K-1 )*3+3 )
                  NR = MIN( 3, I-K+1 )
                  IF( ( NR.EQ.3 ) .AND. ( KCOL( KI ).LE.KP2COL( KI ) ) )
     $                 THEN
*
                     IF( ( K.LT.ISTOP ) .AND.
     $                   ( MOD( K-1, HBL ).LT.HBL-2 ) ) THEN
                        ITMP1 = MIN( ISTART+1, I ) - 1
                     ELSE
                        IF( MOD( K-1, HBL ).LT.HBL-2 ) THEN
                           ITMP1 = MIN( K+3, I )
                        END IF
                        IF( MOD( K-1, HBL ).EQ.HBL-2 ) THEN
                           ITMP1 = MAX( I1, K-1 ) - 1
                        END IF
                        IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                           ITMP1 = MAX( I1, K-2 ) - 1
                        END IF
                     END IF
                     IF( MOD( K-1, HBL ).LT.HBL-2 ) THEN
                        ICOL1 = KCOL( KI ) + K - ISTART
                        ICOL2 = KP2COL( KI ) + K - ISTART
                     ELSE
                        ICOL1 = KCOL( KI )
                        ICOL2 = KP2COL( KI )
                        IF( K.GT.ISTART ) THEN
                           IF( RIGHT.EQ.ICURCOL( KI ) ) THEN
                              ICOL1 = ICOL1 + 1
                           END IF
                           IF( MYCOL.EQ.ICURCOL( KI ) ) THEN
                              ICOL2 = ICOL2 + 1
                           END IF
                        END IF
                     END IF
                     CALL INFOG1L( I1, HBL, NPROW, MYROW, IAFIRST,
     $                             IROW1, IROW2 )
                     IROW2 = NUMROC( ITMP1, HBL, MYROW, IAFIRST, NPROW )
                     IF( ( MOD( K-1, HBL ).EQ.HBL-2 ) .AND.
     $                   ( NPCOL.GT.1 ) ) THEN
                        IF( ICOL1.NE.ICOL2 ) THEN
                           IF( ISTART.EQ.ISTOP ) THEN
                              CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      A( ( ICOL1-1 )*LDA+IROW1 ),
     $                                      LDA, MYROW, RIGHT )
                           END IF
                        END IF
                     END IF
                     IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                   ( NPCOL.GT.1 ) ) THEN
                        IF( ICOL1.EQ.ICOL2 ) THEN
                           CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                   A( ( ICOL1-2 )*LDA+IROW1 ),
     $                                   LDA, MYROW, RIGHT )
                        END IF
                     END IF
*
*                    If we want Z and we haven't already done any Z
*
                     IF( ( WANTZ ) .AND. ( MOD( K-1,
     $                   HBL ).GE.HBL-2 ) .AND. ( NPCOL.GT.1 ) ) THEN
*
*                       Accumulate transformations in the matrix Z
*
                        IROW1 = LILOZ
                        IROW2 = LIHIZ
                        IF( MOD( K-1, HBL ).EQ.HBL-2 ) THEN
                           IF( ICOL1.NE.ICOL2 ) THEN
                              IF( ISTART.EQ.ISTOP ) THEN
                                 CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                         Z( ( ICOL1-1 )*LDZ+
     $                                         IROW1 ), LDZ, MYROW,
     $                                         RIGHT )
                              END IF
                           END IF
                        END IF
                        IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                           IF( ICOL1.EQ.ICOL2 ) THEN
                              CALL CGERV2D( CONTXT, IROW2-IROW1+1, 2,
     $                                      Z( ( ICOL1-2 )*LDZ+IROW1 ),
     $                                      LDZ, MYROW, RIGHT )
                           END IF
                        END IF
                     END IF
                  END IF
  430          CONTINUE
  440       CONTINUE
*
*           Column work done
*
  450       CONTINUE
*
*           Now do NR=2 work
*
            DO 530 KI = 1, IBULGE
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
               IF( MOD( ISTART-1, HBL ).GE.HBL-2 ) THEN
*
*                 INFO is found in a buffer
*
                  ISPEC = 1
               ELSE
*
*                 All INFO is local
*
                  ISPEC = 0
               END IF
*
               DO 520 K = ISTART, ISTOP
*
                  V2 = WORK( VECSIDX+( K-1 )*3+1 )
                  V3 = WORK( VECSIDX+( K-1 )*3+2 )
                  T1 = WORK( VECSIDX+( K-1 )*3+3 )
                  NR = MIN( 3, I-K+1 )
                  IF( NR.EQ.2 ) THEN
                     IF ( ICURROW( KI ).EQ.MYROW ) THEN
                        T2 = T1*V2
                     END IF
                     IF ( ICURCOL( KI ).EQ.MYCOL ) THEN
                        T2 = T1*V2
                     END IF
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
                     CALL INFOG1L( K, HBL, NPCOL, MYCOL, JAFIRST, LILOH,
     $                             LIHIH )
                     LIHIH = LOCALI2
                     CALL INFOG1L( 1, HBL, NPROW, MYROW, IAFIRST, ITMP2,
     $                             ITMP1 )
                     ITMP1 = NUMROC( K+1, HBL, MYROW, IAFIRST, NPROW )
                     IF( ICURROW( KI ).EQ.MYROW ) THEN
                        IF( ( ISPEC.EQ.0 ) .OR. ( NPROW.EQ.1 ) .OR.
     $                      ( MOD( K-1, HBL ).EQ.HBL-2 ) ) THEN
                           ITMP1 = ITMP1 - 1
                           DO 460 J = ( LILOH-1 )*LDA,
     $                             ( LIHIH-1 )*LDA, LDA
                              SUM = CONJG( T1 )*A( ITMP1+J ) +
     $                              CONJG( T2 )*A( ITMP1+1+J )
                              A( ITMP1+J ) = A( ITMP1+J ) - SUM
                              A( ITMP1+1+J ) = A( ITMP1+1+J ) - SUM*V2
  460                      CONTINUE
                        ELSE
                           IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                              CALL CGERV2D( CONTXT, 1, LIHIH-LILOH+1,
     $                                      WORK( IRBUF+1 ), 1, UP,
     $                                      MYCOL )
                              DO 470 J = LILOH, LIHIH
                                 SUM = CONJG( T1 )*
     $                                 WORK( IRBUF+J-LILOH+1 ) +
     $                                 CONJG( T2 )*A( ( J-1 )*LDA+
     $                                 ITMP1 )
                                 WORK( IRBUF+J-LILOH+1 ) = WORK( IRBUF+
     $                              J-LILOH+1 ) - SUM
                                 A( ( J-1 )*LDA+ITMP1 ) = A( ( J-1 )*
     $                              LDA+ITMP1 ) - SUM*V2
  470                         CONTINUE
                              CALL CGESD2D( CONTXT, 1, LIHIH-LILOH+1,
     $                                      WORK( IRBUF+1 ), 1, UP,
     $                                      MYCOL )
                           END IF
                        END IF
                     ELSE
                        IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                      ( ICURROW( KI ).EQ.DOWN ) ) THEN
                           CALL CGESD2D( CONTXT, 1, LIHIH-LILOH+1,
     $                                   A( ( LILOH-1 )*LDA+ITMP1 ),
     $                                   LDA, DOWN, MYCOL )
                           CALL CGERV2D( CONTXT, 1, LIHIH-LILOH+1,
     $                                   A( ( LILOH-1 )*LDA+ITMP1 ),
     $                                   LDA, DOWN, MYCOL )
                        END IF
                     END IF
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to MIN(K+3,I).
*
                     CALL INFOG1L( I1, HBL, NPROW, MYROW, IAFIRST,
     $                             LILOH, LIHIH )
                     LIHIH = NUMROC( I, HBL, MYROW, IAFIRST, NPROW )
*
                     IF( ICURCOL( KI ).EQ.MYCOL ) THEN
*                       LOCAL A(LILOZ:LIHIZ,KCOL:KCOL+2)
                        IF( ( ISPEC.EQ.0 ) .OR. ( NPCOL.EQ.1 ) .OR.
     $                      ( MOD( K-1, HBL ).EQ.HBL-2 ) ) THEN
                           CALL INFOG1L( K, HBL, NPCOL, MYCOL, JAFIRST,
     $                                   ITMP1, ITMP2 )
                           ITMP2 = NUMROC( K+1, HBL, MYCOL, JAFIRST,
     $                             NPCOL )
                           DO 480 J = LILOH, LIHIH
                              SUM = T1*A( ( ITMP1-1 )*LDA+J ) +
     $                              T2*A( ITMP1*LDA+J )
                              A( ( ITMP1-1 )*LDA+J ) = A( ( ITMP1-1 )*
     $                           LDA+J ) - SUM
                              A( ITMP1*LDA+J ) = A( ITMP1*LDA+J ) -
     $                                           SUM*CONJG( V2 )
  480                      CONTINUE
                        ELSE
                           ITMP1 = KCOL( KI )
                           IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                              CALL CGERV2D( CONTXT, LIHIH-LILOH+1, 1,
     $                                      WORK( ICBUF+1 ),
     $                                      LIHIH-LILOH+1, MYROW, LEFT )
                              DO 490 J = LILOH, LIHIH
                                 SUM = T1*WORK( ICBUF+J ) +
     $                                 T2*A( ( ITMP1-1 )*LDA+J )
                                 WORK( ICBUF+J ) = WORK( ICBUF+J ) - SUM
                                 A( ( ITMP1-1 )*LDA+J )
     $                              = A( ( ITMP1-1 )*LDA+J ) -
     $                              SUM*CONJG( V2 )
  490                         CONTINUE
                              CALL CGESD2D( CONTXT, LIHIH-LILOH+1, 1,
     $                                      WORK( ICBUF+1 ),
     $                                      LIHIH-LILOH+1, MYROW, LEFT )
                           END IF
                        END IF
                     ELSE
                        IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                      ( ICURCOL( KI ).EQ.RIGHT ) ) THEN
                           ITMP1 = KCOL( KI )
                           CALL CGESD2D( CONTXT, LIHIH-LILOH+1, 1,
     $                                   A( ( ITMP1-1 )*LDA+LILOH ),
     $                                   LDA, MYROW, RIGHT )
                           CALL INFOG1L( K, HBL, NPCOL, MYCOL, JAFIRST,
     $                                   ITMP1, ITMP2 )
                           ITMP2 = NUMROC( K+1, HBL, MYCOL, JAFIRST,
     $                             NPCOL )
                           CALL CGERV2D( CONTXT, LIHIH-LILOH+1, 1,
     $                                   A( ( ITMP1-1 )*LDA+LILOH ),
     $                                   LDA, MYROW, RIGHT )
                        END IF
                     END IF
*
                     IF( WANTZ ) THEN
*
*                       Accumulate transformations in the matrix Z
*
                        IF( ICURCOL( KI ).EQ.MYCOL ) THEN
*                          LOCAL Z(LILOZ:LIHIZ,KCOL:KCOL+2)
                           IF( ( ISPEC.EQ.0 ) .OR. ( NPCOL.EQ.1 ) .OR.
     $                         ( MOD( K-1, HBL ).EQ.HBL-2 ) ) THEN
                              ITMP1 = KCOL( KI ) + K - ISTART
                              ITMP1 = ( ITMP1-1 )*LDZ
                              DO 500 J = LILOZ, LIHIZ
                                 SUM = T1*Z( J+ITMP1 ) +
     $                                 T2*Z( J+ITMP1+LDZ )
                                 Z( J+ITMP1 ) = Z( J+ITMP1 ) - SUM
                                 Z( J+ITMP1+LDZ ) = Z( J+ITMP1+LDZ ) -
     $                                              SUM*CONJG( V2 )
  500                         CONTINUE
                           ELSE
                              ITMP1 = KCOL( KI )
*                             IF WE ACTUALLY OWN COLUMN K
                              IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                                 CALL CGERV2D( CONTXT, LIHIZ-LILOZ+1, 1,
     $                                         WORK( IZBUF+1 ), LDZ,
     $                                         MYROW, LEFT )
                                 ITMP1 = ( ITMP1-1 )*LDZ
                                 DO 510 J = LILOZ, LIHIZ
                                    SUM = T1*WORK( IZBUF+J ) +
     $                                    T2*Z( J+ITMP1 )
                                    WORK( IZBUF+J ) = WORK( IZBUF+J ) -
     $                                 SUM
                                    Z( J+ITMP1 ) = Z( J+ITMP1 ) -
     $                                             SUM*CONJG( V2 )
  510                            CONTINUE
                                 CALL CGESD2D( CONTXT, LIHIZ-LILOZ+1, 1,
     $                                         WORK( IZBUF+1 ), LDZ,
     $                                         MYROW, LEFT )
                              END IF
                           END IF
                        ELSE
*
*                          NO WORK BUT NEED TO UPDATE ANYWAY????
*
                           IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                         ( ICURCOL( KI ).EQ.RIGHT ) ) THEN
                              ITMP1 = KCOL( KI )
                              ITMP1 = ( ITMP1-1 )*LDZ
                              CALL CGESD2D( CONTXT, LIHIZ-LILOZ+1, 1,
     $                                      Z( LILOZ+ITMP1 ), LDZ,
     $                                      MYROW, RIGHT )
                              CALL CGERV2D( CONTXT, LIHIZ-LILOZ+1, 1,
     $                                      Z( LILOZ+ITMP1 ), LDZ,
     $                                      MYROW, RIGHT )
                           END IF
                        END IF
                     END IF
                  END IF
  520          CONTINUE
*
*        Adjust local information for this bulge
*
               IF( NPROW.EQ.1 ) THEN
                  KROW( KI ) = KROW( KI ) + K2( KI ) - K1( KI ) + 1
                  KP2ROW( KI ) = KP2ROW( KI ) + K2( KI ) - K1( KI ) + 1
               END IF
               IF( ( MOD( K1( KI )-1, HBL ).LT.HBL-2 ) .AND.
     $             ( ICURROW( KI ).EQ.MYROW ) .AND. ( NPROW.GT.1 ) )
     $              THEN
                  KROW( KI ) = KROW( KI ) + K2( KI ) - K1( KI ) + 1
               END IF
               IF( ( MOD( K2( KI ), HBL ).LT.HBL-2 ) .AND.
     $             ( ICURROW( KI ).EQ.MYROW ) .AND. ( NPROW.GT.1 ) )
     $              THEN
                  KP2ROW( KI ) = KP2ROW( KI ) + K2( KI ) - K1( KI ) + 1
               END IF
               IF( ( MOD( K1( KI )-1, HBL ).GE.HBL-2 ) .AND.
     $             ( ( MYROW.EQ.ICURROW( KI ) ) .OR. ( DOWN.EQ.
     $             ICURROW( KI ) ) ) .AND. ( NPROW.GT.1 ) ) THEN
                  CALL INFOG1L( K2( KI )+1, HBL, NPROW, MYROW, IAFIRST,
     $                          KROW( KI ), ITMP2 )
               END IF
               IF( ( MOD( K2( KI ), HBL ).GE.HBL-2 ) .AND.
     $             ( ( MYROW.EQ.ICURROW( KI ) ) .OR. ( UP.EQ.
     $             ICURROW( KI ) ) ) .AND. ( NPROW.GT.1 ) ) THEN
                  KP2ROW( KI ) = NUMROC( K2( KI )+3, HBL, MYROW,
     $                           IAFIRST, NPROW )
               END IF
               IF( NPCOL.EQ.1 ) THEN
                  KCOL( KI ) = KCOL( KI ) + K2( KI ) - K1( KI ) + 1
                  KP2COL( KI ) = KP2COL( KI ) + K2( KI ) - K1( KI ) + 1
               END IF
               IF( ( MOD( K1( KI )-1, HBL ).LT.HBL-2 ) .AND.
     $             ( ICURCOL( KI ).EQ.MYCOL ) .AND. ( NPCOL.GT.1 ) )
     $              THEN
                  KCOL( KI ) = KCOL( KI ) + K2( KI ) - K1( KI ) + 1
               END IF
               IF( ( MOD( K2( KI ), HBL ).LT.HBL-2 ) .AND.
     $             ( ICURCOL( KI ).EQ.MYCOL ) .AND. ( NPCOL.GT.1 ) )
     $              THEN
                  KP2COL( KI ) = KP2COL( KI ) + K2( KI ) - K1( KI ) + 1
               END IF
               IF( ( MOD( K1( KI )-1, HBL ).GE.HBL-2 ) .AND.
     $             ( ( MYCOL.EQ.ICURCOL( KI ) ) .OR. ( RIGHT.EQ.
     $             ICURCOL( KI ) ) ) .AND. ( NPCOL.GT.1 ) ) THEN
                  CALL INFOG1L( K2( KI )+1, HBL, NPCOL, MYCOL, JAFIRST,
     $                          KCOL( KI ), ITMP2 )
               END IF
               IF( ( MOD( K2( KI ), HBL ).GE.HBL-2 ) .AND.
     $             ( ( MYCOL.EQ.ICURCOL( KI ) ) .OR. ( LEFT.EQ.
     $             ICURCOL( KI ) ) ) .AND. ( NPCOL.GT.1 ) ) THEN
                  KP2COL( KI ) = NUMROC( K2( KI )+3, HBL, MYCOL,
     $                           JAFIRST, NPCOL )
               END IF
               K1( KI ) = K2( KI ) + 1
               ISTOP = MIN( K1( KI )+ROTN-MOD( K1( KI ), ROTN ), I-2 )
               ISTOP = MIN( ISTOP, K1( KI )+HBL-3-
     $                 MOD( K1( KI )-1, HBL ) )
               ISTOP = MIN( ISTOP, I2-2 )
               ISTOP = MAX( ISTOP, K1( KI ) )
               IF( ( MOD( K1( KI )-1, HBL ).EQ.HBL-2 ) .AND.
     $             ( ISTOP.LT.MIN( I-2, I2-2 ) ) ) THEN
                  ISTOP = ISTOP + 1
               END IF
               K2( KI ) = ISTOP
               IF( K1( KI ).LE.ISTOP ) THEN
                  IF( ( MOD( K1( KI )-1, HBL ).EQ.HBL-2 ) .AND.
     $                ( I-K1( KI ).GT.1 ) ) THEN
*
*                    Next step switches rows & cols
*
                     ICURROW( KI ) = MOD( ICURROW( KI )+1, NPROW )
                     ICURCOL( KI ) = MOD( ICURCOL( KI )+1, NPCOL )
                  END IF
               END IF
  530       CONTINUE
*
            IF( K2( IBULGE ).LE.I-1 )
     $         GO TO 40
         END IF
*
  540 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      INFO = I
      RETURN
*
  550 CONTINUE
*
      IF( L.EQ.I ) THEN
*
*        H(I,I-1) is negligible: one eigenvalue has converged.
*
         CALL INFOG2L( I, I, DESCA, NPROW, NPCOL, MYROW, MYCOL, IROW,
     $                 ICOL, ITMP1, ITMP2 )
         IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
            W( I ) = A( ( ICOL-1 )*LDA+IROW )
         ELSE
            W( I ) = ZERO
         END IF
      ELSE IF( L.EQ.I-1 ) THEN
*
*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
*
         CALL PCLACP3( 2, I-1, A, DESCA, S1, 2*IBLK, -1, -1, 0 )
         CALL CLANV2( S1( 1, 1 ), S1( 1, 2 ), S1( 2, 1 ), S1( 2, 2 ),
     $                W( I-1 ), W( I ), CS, SN )
         CALL PCLACP3( 2, I-1, A, DESCA, S1, 2*IBLK, 0, 0, 1 )
*
         IF( NODE.NE.0 ) THEN
*           Erase the eigenvalues other eigenvalues
            W( I-1 ) = ZERO
            W( I ) = ZERO
         END IF
*
         IF( WANTT ) THEN
*
*           Apply the transformation to A.
*
            IF( I2.GT.I ) THEN
               CALL PCROT( I2-I, A, I-1, I+1, DESCA, N, A, I, I+1,
     $                     DESCA, N, CS, SN )
            END IF
            CALL PCROT( I-I1-1, A, I1, I-1, DESCA, 1, A, I1, I, DESCA,
     $                  1, CS, CONJG( SN ) )
         END IF
         IF( WANTZ ) THEN
*
*           Apply the transformation to Z.
*
            CALL PCROT( NZ, Z, ILOZ, I-1, DESCZ, 1, Z, ILOZ, I, DESCZ,
     $                  1, CS, CONJG( SN ) )
         END IF
*
      ELSE
*
*        Find the eigenvalues in H(L:I,L:I), L < I-1
*
         JBLK = I - L + 1
         IF( JBLK.LE.2*IBLK ) THEN
            CALL PCLACP3( I-L+1, L, A, DESCA, S1, 2*IBLK, 0, 0, 0 )
            CALL CLAHQR2( .FALSE., .FALSE., JBLK, 1, JBLK, S1, 2*IBLK,
     $                   W( L ), 1, JBLK, Z, LDZ, IERR )
            IF( NODE.NE.0 ) THEN
*
*              Erase the eigenvalues
*
               DO 560 K = L, I
                  W( K ) = ZERO
  560          CONTINUE
            END IF
         END IF
      END IF
*
*     Decrement number of remaining iterations, and return to start of
*     the main loop with new value of I.
*
      ITN = ITN - ITS
      I = L - 1
      GO TO 10
*
  570 CONTINUE
      CALL CGSUM2D( CONTXT, 'All', ' ', N, 1, W, N, -1, -1 )
      RETURN
*
*     END OF PCLAHQR
*
      END
