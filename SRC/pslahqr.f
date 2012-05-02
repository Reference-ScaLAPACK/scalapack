      SUBROUTINE PSLAHQR( WANTT, WANTZ, N, ILO, IHI, A, DESCA, WR, WI,
     $                    ILOZ, IHIZ, Z, DESCZ, WORK, LWORK, IWORK,
     $                    ILWORK, INFO )
*
*  -- ScaLAPACK routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      LOGICAL            WANTT, WANTZ
      INTEGER            IHI, IHIZ, ILO, ILOZ, ILWORK, INFO, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * ), IWORK( * )
      REAL               A( * ), WI( * ), WORK( * ), WR( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLAHQR is an auxiliary routine used to find the Schur decomposition
*    and or eigenvalues of a matrix already in Hessenberg form from
*    cols ILO to IHI.
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
*          ILO = 1). PSLAHQR works primarily with the Hessenberg
*          submatrix in rows and columns ILO to IHI, but applies
*          transformations to all of H if WANTT is .TRUE..
*          1 <= ILO <= max(1,IHI); IHI <= N.
*
*  A       (global input/output) REAL array, dimension
*          (DESCA(LLD_),*)
*          On entry, the upper Hessenberg matrix A.
*          On exit, if WANTT is .TRUE., A is upper quasi-triangular in
*          rows and columns ILO:IHI, with any 2-by-2 or larger diagonal
*          blocks not yet in standard form. If WANTT is .FALSE., the
*          contents of A are unspecified on exit.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  WR      (global replicated output) REAL array,
*                                                         dimension (N)
*  WI      (global replicated output) REAL array,
*                                                         dimension (N)
*          The real and imaginary parts, respectively, of the computed
*          eigenvalues ILO to IHI are stored in the corresponding
*          elements of WR and WI. If two eigenvalues are computed as a
*          complex conjugate pair, they are stored in consecutive
*          elements of WR and WI, say the i-th and (i+1)th, with
*          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
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
*  Z       (global input/output) REAL array.
*          If WANTZ is .TRUE., on entry Z must contain the current
*          matrix Z of transformations accumulated by PDHSEQR, and on
*          exit Z has been updated; transformations are applied only to
*          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
*          If WANTZ is .FALSE., Z is not referenced.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*  WORK    (local output) REAL array of size LWORK
*
*  LWORK   (local input) INTEGER
*          WORK(LWORK) is a local array and LWORK is assumed big enough
*          so that LWORK >= 3*N +
*                MAX( 2*MAX(DESCZ(LLD_),DESCA(LLD_)) + 2*LOCc(N),
*                     7*Ceil(N/HBL)/LCM(NPROW,NPCOL)) )
*
*  IWORK   (global and local input) INTEGER array of size ILWORK
*
*  ILWORK  (local input) INTEGER
*          This holds the some of the IBLK integer arrays.  This is held
*          as a place holder for the next release.
*
*  INFO    (global output) INTEGER
*          < 0: parameter number -INFO incorrect or inconsistent
*          = 0: successful exit
*          > 0: PSLAHQR failed to compute all the eigenvalues ILO to IHI
*               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
*               elements i+1:ihi of WR and WI contain those eigenvalues
*               which have been successfully computed.
*
*  Logic:
*       This algorithm is very similar to _LAHQR.  Unlike _LAHQR,
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
*       This routine calls:
*           PSLACONSB   -> To determine where to start each iteration
*           PSLAWIL   -> Given the shift, get the transformation
*           SLASORTE   -> Pair up eigenvalues so that reals are paired.
*           PSLACP3   -> Parallel array to local replicated array copy &
*                        back.
*           SLAREF   -> Row/column reflector applier.  Core routine
*                        here.
*           PSLASMSUB   -> Finds negligible subdiagonal elements.
*
*  Current Notes and/or Restrictions:
*       1.) This code requires the distributed block size to be square
*           and at least six (6); unlike simpler codes like LU, this
*           algorithm is extremely sensitive to block size.  Unwise
*           choices of too small a block size can lead to bad
*           performance.
*       2.) This code requires A and Z to be distributed identically
*           and have identical contxts.
*       3.) This release currently does not have a routine for
*           resolving the Schur blocks into regular 2x2 form after
*           this code is completed.  Because of this, a significant
*           performance impact is required while the deflation is done
*           by sometimes a single column of processors.
*       4.) This code does not currently block the initial transforms
*           so that none of the rows or columns for any bulge are
*           completed until all are started.  To offset pipeline
*           start-up it is recommended that at least 2*LCM(NPROW,NPCOL)
*           bulges are used (if possible)
*       5.) The maximum number of bulges currently supported is fixed at
*           32.  In future versions this will be limited only by the
*           incoming WORK array.
*       6.) The matrix A must be in upper Hessenberg form.  If elements
*           below the subdiagonal are nonzero, the resulting transforms
*           may be nonsimilar.  This is also true with the LAPACK
*           routine.
*       7.) For this release, it is assumed RSRC_=CSRC_=0
*       8.) Currently, all the eigenvalues are distributed to all the
*           nodes.  Future releases will probably distribute the
*           eigenvalues by the column partitioning.
*       9.) The internals of this routine are subject to change.
*
*  Implemented by:  G. Henry, November 17, 1996
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0, ONE = 1.0, HALF = 0.5 )
      REAL               CONST
      PARAMETER          ( CONST = 1.50 )
      INTEGER            IBLK
      PARAMETER          ( IBLK = 32 )
*     ..
*     .. Local Scalars ..
      INTEGER            CONTXT, DOWN, HBL, I, I1, I2, IAFIRST, IBULGE,
     $                   ICBUF, ICOL, ICOL1, ICOL2, IDIA, IERR, II,
     $                   IRBUF, IROW, IROW1, IROW2, ISPEC, ISTART,
     $                   ISTARTCOL, ISTARTROW, ISTOP, ISUB, ISUP,
     $                   ITERMAX, ITMP1, ITMP2, ITN, ITS, J, JAFIRST,
     $                   JBLK, JJ, K, KI, L, LCMRC, LDA, LDZ, LEFT,
     $                   LIHIH, LIHIZ, LILOH, LILOZ, LOCALI1, LOCALI2,
     $                   LOCALK, LOCALM, M, MODKM1, MYCOL, MYROW,
     $                   NBULGE, NH, NODE, NPCOL, NPROW, NR, NUM, NZ,
     $                   RIGHT, ROTN, UP, VECSIDX
      REAL               AVE, DISC, H00, H10, H11, H12, H21, H22, H33,
     $                   H43H34, H44, OVFL, S, SMLNUM, SUM, T1, T1COPY,
     $                   T2, T3, ULP, UNFL, V1SAVE, V2, V2SAVE, V3,
     $                   V3SAVE, CS, SN
*     ..
*     .. Local Arrays ..
      INTEGER            ICURCOL( IBLK ), ICURROW( IBLK ), K1( IBLK ),
     $                   K2( IBLK ), KCOL( IBLK ), KP2COL( IBLK ),
     $                   KP2ROW( IBLK ), KROW( IBLK ), LOCALK2( IBLK )
      REAL               S1( 2*IBLK, 2*IBLK ), SMALLA( 6, 6, IBLK ),
     $                   VCOPY( 3 )
*     ..
*     .. External Functions ..
      INTEGER            ILCM, NUMROC
      REAL               PSLAMCH
      EXTERNAL           ILCM, NUMROC, PSLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, SCOPY, SGEBR2D, SGEBS2D,
     $                   SGERV2D, SGESD2D, SGSUM2D, SLAHQR, SLAREF,
     $                   SLARFG, SLASORTE, IGAMN2D, INFOG1L, INFOG2L,
     $                   PSLABAD, PSLACONSB, PSLACP3, PSLASMSUB,
     $                   PSLAWIL, PXERBLA, SLANV2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
      ITERMAX = 30*( IHI-ILO+1 )
*     ITERMAX = 0
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
*
*     Determine the number of columns we have so we can check workspace
*
      LOCALK = NUMROC( N, HBL, MYCOL, JAFIRST, NPCOL )
      JJ = N / HBL
      IF( JJ*HBL.LT.N )
     $   JJ = JJ + 1
      JJ = 7*JJ / LCMRC
      IF( LWORK.LT.3*N+MAX( 2*MAX( LDA, LDZ )+2*LOCALK, JJ ) ) THEN
         INFO = -15
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
         CALL PXERBLA( CONTXT, 'PSLAHQR', -INFO )
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
*
*     Find a value for ROTN
*
      ROTN = HBL / 3
      ROTN = MAX( ROTN, HBL-2 )
      ROTN = MIN( ROTN, 1 )
*
      IF( ILO.EQ.IHI ) THEN
         CALL INFOG2L( ILO, ILO, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 IROW, ICOL, II, JJ )
         IF( ( MYROW.EQ.II ) .AND. ( MYCOL.EQ.JJ ) ) THEN
            WR( ILO ) = A( ( ICOL-1 )*LDA+IROW )
         ELSE
            WR( ILO ) = ZERO
         END IF
         WI( ILO ) = ZERO
         RETURN
      END IF
*
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
*
      CALL INFOG1L( ILOZ, HBL, NPROW, MYROW, 0, LILOZ, LIHIZ )
      LIHIZ = NUMROC( IHIZ, HBL, MYROW, 0, NPROW )
*
*     Set machine-dependent constants for the stopping criterion.
*     If NORM(H) <= SQRT(OVFL), overflow should not occur.
*
      UNFL = PSLAMCH( CONTXT, 'SAFE MINIMUM' )
      OVFL = ONE / UNFL
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
     $   GO TO 450
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 or 2 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      DO 420 ITS = 0, ITN
*
*        Look for a single small subdiagonal element.
*
         CALL PSLASMSUB( A, DESCA, I, L, K, SMLNUM, WORK( IRBUF+1 ),
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
         M = L - 10
*        IF ( L .GE. I - (2*IBLK-1) )
*         IF ( L .GE. I - MAX(2*IBLK-1,HBL) )
         IF( L.GE.I-1 )
     $      GO TO 430
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
         CALL PSLACP3( 2*JBLK, I-2*JBLK+1, A, DESCA, S1, 2*IBLK, -1, -1,
     $                 0 )
         IF( ITS.EQ.20 .OR. ITS.EQ.40 ) THEN
*
*           Exceptional shift.
*
            DO 20 II = 2*JBLK, 2, -1
               S1( II, II ) = CONST*( ABS( S1( II, II ) )+
     $                        ABS( S1( II, II-1 ) ) )
               S1( II, II-1 ) = ZERO
               S1( II-1, II ) = ZERO
   20       CONTINUE
            S1( 1, 1 ) = CONST*ABS( S1( 1, 1 ) )
         ELSE
            CALL SLAHQR( .FALSE., .FALSE., 2*JBLK, 1, 2*JBLK, S1,
     $                   2*IBLK, WORK( IRBUF+1 ), WORK( ICBUF+1 ), 1,
     $                   2*JBLK, Z, LDZ, IERR )
*
*           Prepare to use Wilkinson's double shift
*
            H44 = S1( 2*JBLK, 2*JBLK )
            H33 = S1( 2*JBLK-1, 2*JBLK-1 )
            H43H34 = S1( 2*JBLK-1, 2*JBLK )*S1( 2*JBLK, 2*JBLK-1 )
            IF( ( JBLK.GT.1 ) .AND. ( ITS.GT.30 ) ) THEN
               S = S1( 2*JBLK-1, 2*JBLK-2 )
               DISC = ( H33-H44 )*HALF
               DISC = DISC*DISC + H43H34
               IF( DISC.GT.ZERO ) THEN
*
*                 Real roots: Use Wilkinson's shift twice
*
                  DISC = SQRT( DISC )
                  AVE = HALF*( H33+H44 )
                  IF( ABS( H33 )-ABS( H44 ).GT.ZERO ) THEN
                     H33 = H33*H44 - H43H34
                     H44 = H33 / ( SIGN( DISC, AVE )+AVE )
                  ELSE
                     H44 = SIGN( DISC, AVE ) + AVE
                  END IF
                  H33 = H44
                  H43H34 = ZERO
               END IF
            END IF
         END IF
*
*        Look for two consecutive small subdiagonal elements:
*           PSLACONSB is the routine that does this.
*
c         CALL PSLACONSB( A, DESCA, I, L, M, H44, H33, H43H34,
c     $                   WORK( IRBUF+1 ), LWORK-IRBUF )
*
*        Skip small submatrices
*
*        IF ( M .GE. I - 5 )
*    $      GO TO 80
*
*        In principle PSLACONSB needs to check all shifts to decide
*        whether two consecutive small subdiagonal entries are suitable
*        as the starting position of the bulge chasing phase. It can be
*        dangerous to check the first pair of shifts only. Moreover it
*        is quite rare to obtain an M which is much larger than L. This
*        process is a bit expensive compared with the benefit.
*        Therefore it is sensible to abandon this routine. Total amount
*        of communications is saved in average.
*
         M = L
*        Double-shift QR step
*
*        NBULGE is the number of bulges that will be attempted
*
         ISTOP = MIN( M+ROTN-MOD( M, ROTN ), I-2 )
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
         IF( ( ITS.NE.20 ) .AND. ( ITS.NE.40 ) .AND. ( NBULGE.GT.1 ) )
     $        THEN
*
*           sort the eigenpairs so that they are in twos for double
*           shifts.  only call if several need sorting
*
            CALL SLASORTE( S1( 2*( JBLK-NBULGE )+1,
     $                     2*( JBLK-NBULGE )+1 ), 2*IBLK, 2*NBULGE,
     $                     WORK( IRBUF+1 ), IERR )
         END IF
*
*        IBULGE is the number of bulges going so far
*
         IBULGE = 1
*
*        "A" row defs : main row transforms from LOCALK to LOCALI2
*
         CALL INFOG1L( M, HBL, NPCOL, MYCOL, 0, ITMP1, LOCALK )
         LOCALK = NUMROC( N, HBL, MYCOL, 0, NPCOL )
         CALL INFOG1L( 1, HBL, NPCOL, MYCOL, 0, ICOL1, LOCALI2 )
         LOCALI2 = NUMROC( I2, HBL, MYCOL, 0, NPCOL )
*
*        "A" col defs : main col transforms from LOCALI1 to LOCALM
*
         CALL INFOG1L( I1, HBL, NPROW, MYROW, 0, LOCALI1, ICOL1 )
         ICOL1 = NUMROC( N, HBL, MYROW, 0, NPROW )
         CALL INFOG1L( 1, HBL, NPROW, MYROW, 0, LOCALM, ICOL1 )
         ICOL1 = NUMROC( MIN( M+3, I ), HBL, MYROW, 0, NPROW )
*
*        Which row & column will start the bulges
*
         ISTARTROW = MOD( ( M+1 ) / HBL, NPROW ) + IAFIRST
         ISTARTCOL = MOD( ( M+1 ) / HBL, NPCOL ) + JAFIRST
*
         CALL INFOG1L( M, HBL, NPROW, MYROW, 0, II, ITMP2 )
         ITMP2 = NUMROC( N, HBL, MYROW, 0, NPROW )
         CALL INFOG1L( M, HBL, NPCOL, MYCOL, 0, JJ, ITMP2 )
         ITMP2 = NUMROC( N, HBL, MYCOL, 0, NPCOL )
         CALL INFOG1L( 1, HBL, NPROW, MYROW, 0, ISTOP, KP2ROW( 1 ) )
         KP2ROW( 1 ) = NUMROC( M+2, HBL, MYROW, 0, NPROW )
         CALL INFOG1L( 1, HBL, NPCOL, MYCOL, 0, ISTOP, KP2COL( 1 ) )
         KP2COL( 1 ) = NUMROC( M+2, HBL, MYCOL, 0, NPCOL )
*
*        Set all values for bulges.  All bulges are stored in
*          intermediate steps as loops over KI.  Their current "task"
*          over the global M to I-1 values is always K1(KI) to K2(KI).
*          However, because there are many bulges, K1(KI) & K2(KI) might
*          go past that range while later bulges (KI+1,KI+2,etc..) are
*          finishing up.
*
*        Rules:
*              If MOD(K1(KI)-1,HBL) < HBL-2 then MOD(K2(KI)-1,HBL)<HBL-2
*              If MOD(K1(KI)-1,HBL) = HBL-2 then MOD(K2(KI)-1,HBL)=HBL-2
*              If MOD(K1(KI)-1,HBL) = HBL-1 then MOD(K2(KI)-1,HBL)=HBL-1
*              K2(KI)-K1(KI) <= ROTN
*
*        We first hit a border when MOD(K1(KI)-1,HBL)=HBL-2 and we hit
*        it again when MOD(K1(KI)-1,HBL)=HBL-1.
*
         DO 30 KI = 1, NBULGE
            K1( KI ) = M
            ISTOP = MIN( M+ROTN-MOD( M, ROTN ), I-2 )
            ISTOP = MIN( ISTOP, M+HBL-3-MOD( M-1, HBL ) )
            ISTOP = MIN( ISTOP, I2-2 )
            ISTOP = MAX( ISTOP, M )
            K2( KI ) = ISTOP
            ICURROW( KI ) = ISTARTROW
            ICURCOL( KI ) = ISTARTCOL
            LOCALK2( KI ) = ITMP1
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
         CALL PSLAWIL( ITMP1, ITMP2, M, A, DESCA, H44, H33, H43H34,
     $                 VCOPY )
         V1SAVE = VCOPY( 1 )
         V2SAVE = VCOPY( 2 )
         V3SAVE = VCOPY( 3 )
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
                  CALL PSLAWIL( ITMP1, ITMP2, M, A, DESCA, H44, H33,
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
            DO 80 KI = 1, IBULGE
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
                     CALL INFOG2L( K+2, K+2, DESCA, NPROW, NPCOL, MYROW,
     $                             MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                     CALL PSLACP3( MIN( 6, N-K+2 ), K-1, A, DESCA,
     $                             SMALLA( 1, 1, KI ), 6, ITMP1, ITMP2,
     $                             0 )
                  END IF
                  IF( MODKM1.EQ.HBL-1 ) THEN
*
*                 Copy 6 elements from global A(K-2:K+3,K-2:K+3)
*
                     CALL INFOG2L( K+1, K+1, DESCA, NPROW, NPCOL, MYROW,
     $                             MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                     CALL PSLACP3( MIN( 6, N-K+3 ), K-2, A, DESCA,
     $                             SMALLA( 1, 1, KI ), 6, ITMP1, ITMP2,
     $                             0 )
                  END IF
               END IF
*
*           SLAHQR used to have a single row application and a single
*              column application to H.  Here we do something a little
*              more clever.  We break each transformation down into 3
*              parts:
*                  1.) The minimum amount of work it takes to determine
*                        a group of ROTN transformations (this is on
*                        the critical path.) (Loops 130-180)
*                  2.) The small work it takes so that each of the rows
*                        and columns is at the same place.  For example,
*                        all ROTN row transforms are all complete
*                        through some column TMP.  (Loops within 190)
*                  3.) The majority of the row and column transforms
*                        are then applied in a block fashion.
*                        (Loops 290 on.)
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
               IF( ( MYROW.EQ.ICURROW( KI ) ) .AND.
     $             ( MYCOL.EQ.ICURCOL( KI ) ) .AND.
     $             ( MODKM1.EQ.HBL-2 ) .AND.
     $             ( ISTART.LT.MIN( I-1, ISTOP+1 ) ) ) THEN
                  K = ISTART
                  NR = MIN( 3, I-K+1 )
                  IF( K.GT.M ) THEN
                     CALL SCOPY( NR, SMALLA( 2, 1, KI ), 1, VCOPY, 1 )
                  ELSE
                     VCOPY( 1 ) = V1SAVE
                     VCOPY( 2 ) = V2SAVE
                     VCOPY( 3 ) = V3SAVE
                  END IF
                  CALL SLARFG( NR, VCOPY( 1 ), VCOPY( 2 ), 1, T1COPY )
                  IF( K.GT.M ) THEN
                     SMALLA( 2, 1, KI ) = VCOPY( 1 )
                     SMALLA( 3, 1, KI ) = ZERO
                     IF( K.LT.I-1 )
     $                  SMALLA( 4, 1, KI ) = ZERO
                  ELSE IF( M.GT.L ) THEN
                     SMALLA( 2, 1, KI ) = -SMALLA( 2, 1, KI )
                  END IF
                  V2 = VCOPY( 2 )
                  T2 = T1COPY*V2
                  WORK( VECSIDX+( K-1 )*3+1 ) = VCOPY( 2 )
                  WORK( VECSIDX+( K-1 )*3+2 ) = VCOPY( 3 )
                  WORK( VECSIDX+( K-1 )*3+3 ) = T1COPY
               END IF
*
               IF( ( MOD( ISTOP-1, HBL ).EQ.HBL-1 ) .AND.
     $             ( MYROW.EQ.ICURROW( KI ) ) .AND.
     $             ( MYCOL.EQ.ICURCOL( KI ) ) .AND.
     $             ( ISTART.LE.MIN( I, ISTOP ) ) ) THEN
                  K = ISTART
                  NR = MIN( 3, I-K+1 )
                  IF( K.GT.M ) THEN
                     CALL SCOPY( NR, SMALLA( 3, 2, KI ), 1, VCOPY, 1 )
                  ELSE
                     VCOPY( 1 ) = V1SAVE
                     VCOPY( 2 ) = V2SAVE
                     VCOPY( 3 ) = V3SAVE
                  END IF
                  CALL SLARFG( NR, VCOPY( 1 ), VCOPY( 2 ), 1, T1COPY )
                  IF( K.GT.M ) THEN
                     SMALLA( 3, 2, KI ) = VCOPY( 1 )
                     SMALLA( 4, 2, KI ) = ZERO
                     IF( K.LT.I-1 )
     $                  SMALLA( 5, 2, KI ) = ZERO
*
*                 Set a subdiagonal to zero now if it's possible
*
*                 H11 = SMALLA(1,1,KI)
*                 H10 = SMALLA(2,1,KI)
*                 H22 = SMALLA(2,2,KI)
*                 IF ( ABS(H10) .LE. MAX(ULP*(ABS(H11)+ABS(H22)),
*    $                                    SMLNUM) ) THEN
*                    SMALLA(2,1,KI) = ZERO
*     WORK(ISUB+K-2) = ZERO
*                 END IF
                  ELSE IF( M.GT.L ) THEN
                     SMALLA( 3, 2, KI ) = -SMALLA( 3, 2, KI )
                  END IF
                  V2 = VCOPY( 2 )
                  T2 = T1COPY*V2
                  WORK( VECSIDX+( K-1 )*3+1 ) = VCOPY( 2 )
                  WORK( VECSIDX+( K-1 )*3+2 ) = VCOPY( 3 )
                  WORK( VECSIDX+( K-1 )*3+3 ) = T1COPY
               END IF
*
               IF( ( MODKM1.EQ.0 ) .AND. ( ISTART.LE.I-1 ) .AND.
     $             ( MYROW.EQ.ICURROW( KI ) ) .AND.
     $             ( RIGHT.EQ.ICURCOL( KI ) ) ) THEN
*
*              (IROW1,ICOL1) is (I,J)-coordinates of H(ISTART,ISTART)
*
                  IROW1 = KROW( KI )
                  ICOL1 = LOCALK2( KI )
                  IF( ISTART.GT.M ) THEN
                     VCOPY( 1 ) = SMALLA( 4, 3, KI )
                     VCOPY( 2 ) = SMALLA( 5, 3, KI )
                     VCOPY( 3 ) = SMALLA( 6, 3, KI )
                     NR = MIN( 3, I-ISTART+1 )
                     CALL SLARFG( NR, VCOPY( 1 ), VCOPY( 2 ), 1,
     $                            T1COPY )
                     A( ( ICOL1-2 )*LDA+IROW1 ) = VCOPY( 1 )
                     A( ( ICOL1-2 )*LDA+IROW1+1 ) = ZERO
                     IF( ISTART.LT.I-1 ) THEN
                        A( ( ICOL1-2 )*LDA+IROW1+2 ) = ZERO
                     END IF
                  ELSE
                     IF( M.GT.L ) THEN
                        A( ( ICOL1-2 )*LDA+IROW1 ) = -A( ( ICOL1-2 )*
     $                     LDA+IROW1 )
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
*           (IROW1,ICOL1) is (I,J)-coordinates of H(ISTART,ISTART)
*
                  IROW1 = KROW( KI )
                  ICOL1 = LOCALK2( KI )
                  DO 70 K = ISTART, ISTOP
*
*              Create and do these transforms
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
                     CALL SLARFG( NR, VCOPY( 1 ), VCOPY( 2 ), 1,
     $                            T1COPY )
                     IF( K.GT.M ) THEN
                        IF( MOD( K-1, HBL ).GT.0 ) THEN
                           A( ( ICOL1-2 )*LDA+IROW1 ) = VCOPY( 1 )
                           A( ( ICOL1-2 )*LDA+IROW1+1 ) = ZERO
                           IF( K.LT.I-1 ) THEN
                              A( ( ICOL1-2 )*LDA+IROW1+2 ) = ZERO
                           END IF
*
*                    Set a subdiagonal to zero now if it's possible
*
*                    IF ( (IROW1.GT.2) .AND. (ICOL1.GT.2) .AND.
*    $                    (MOD(K-1,HBL) .GT. 1) ) THEN
*                       H11 = A((ICOL1-3)*LDA+IROW1-2)
*                       H10 = A((ICOL1-3)*LDA+IROW1-1)
*                       H22 = A((ICOL1-2)*LDA+IROW1-1)
*                       IF ( ABS(H10).LE.MAX(ULP*(ABS(H11)+ABS(H22)),
*    $                                       SMLNUM) ) THEN
*                           A((ICOL1-3)*LDA+IROW1-1) = ZERO
*                       END IF
*                    END IF
                        END IF
                     ELSE IF( M.GT.L ) THEN
                        IF( MOD( K-1, HBL ).GT.0 ) THEN
                           A( ( ICOL1-2 )*LDA+IROW1 ) = -A( ( ICOL1-2 )*
     $                        LDA+IROW1 )
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
*                 Do some work so next step is ready...
*
                        V3 = VCOPY( 3 )
                        T3 = T1*V3
                        DO 50 J = ICOL1, MIN( K2( KI )+1, I-1 ) +
     $                          ICOL1 - K
                           SUM = A( ( J-1 )*LDA+IROW1 ) +
     $                           V2*A( ( J-1 )*LDA+IROW1+1 ) +
     $                           V3*A( ( J-1 )*LDA+IROW1+2 )
                           A( ( J-1 )*LDA+IROW1 ) = A( ( J-1 )*LDA+
     $                        IROW1 ) - SUM*T1
                           A( ( J-1 )*LDA+IROW1+1 ) = A( ( J-1 )*LDA+
     $                        IROW1+1 ) - SUM*T2
                           A( ( J-1 )*LDA+IROW1+2 ) = A( ( J-1 )*LDA+
     $                        IROW1+2 ) - SUM*T3
   50                   CONTINUE
                        ITMP1 = LOCALK2( KI )
                        DO 60 J = IROW1 + 1, IROW1 + 3
                           SUM = A( ( ICOL1-1 )*LDA+J ) +
     $                           V2*A( ICOL1*LDA+J ) +
     $                           V3*A( ( ICOL1+1 )*LDA+J )
                           A( ( ICOL1-1 )*LDA+J ) = A( ( ICOL1-1 )*LDA+
     $                        J ) - SUM*T1
                           A( ICOL1*LDA+J ) = A( ICOL1*LDA+J ) - SUM*T2
                           A( ( ICOL1+1 )*LDA+J ) = A( ( ICOL1+1 )*LDA+
     $                        J ) - SUM*T3
   60                   CONTINUE
                     END IF
                     IROW1 = IROW1 + 1
                     ICOL1 = ICOL1 + 1
   70             CONTINUE
               END IF
*
               IF( MODKM1.EQ.HBL-2 ) THEN
                  IF( ( DOWN.EQ.ICURROW( KI ) ) .AND.
     $                ( RIGHT.EQ.ICURCOL( KI ) ) .AND. ( NUM.GT.1 ) )
     $                 THEN
                     CALL SGERV2D( CONTXT, 3, 1,
     $                             WORK( VECSIDX+( ISTART-1 )*3+1 ), 3,
     $                             DOWN, RIGHT )
                  END IF
                  IF( ( MYROW.EQ.ICURROW( KI ) ) .AND.
     $                ( MYCOL.EQ.ICURCOL( KI ) ) .AND. ( NUM.GT.1 ) )
     $                 THEN
                     CALL SGESD2D( CONTXT, 3, 1,
     $                             WORK( VECSIDX+( ISTART-1 )*3+1 ), 3,
     $                             UP, LEFT )
                  END IF
                  IF( ( DOWN.EQ.ICURROW( KI ) ) .AND.
     $                ( NPCOL.GT.1 ) .AND. ( ISTART.LE.ISTOP ) ) THEN
                     JJ = MOD( ICURCOL( KI )+NPCOL-1, NPCOL )
                     IF( MYCOL.NE.JJ ) THEN
                        CALL SGEBR2D( CONTXT, 'ROW', ' ',
     $                                3*( ISTOP-ISTART+1 ), 1,
     $                                WORK( VECSIDX+( ISTART-1 )*3+1 ),
     $                                3*( ISTOP-ISTART+1 ), MYROW, JJ )
                     ELSE
                        CALL SGEBS2D( CONTXT, 'ROW', ' ',
     $                                3*( ISTOP-ISTART+1 ), 1,
     $                                WORK( VECSIDX+( ISTART-1 )*3+1 ),
     $                                3*( ISTOP-ISTART+1 ) )
                     END IF
                  END IF
               END IF
*
*           Broadcast Householder information from the block
*
               IF( ( MYROW.EQ.ICURROW( KI ) ) .AND. ( NPCOL.GT.1 ) .AND.
     $             ( ISTART.LE.ISTOP ) ) THEN
                  IF( MYCOL.NE.ICURCOL( KI ) ) THEN
                     CALL SGEBR2D( CONTXT, 'ROW', ' ',
     $                             3*( ISTOP-ISTART+1 ), 1,
     $                             WORK( VECSIDX+( ISTART-1 )*3+1 ),
     $                             3*( ISTOP-ISTART+1 ), MYROW,
     $                             ICURCOL( KI ) )
                  ELSE
                     CALL SGEBS2D( CONTXT, 'ROW', ' ',
     $                             3*( ISTOP-ISTART+1 ), 1,
     $                             WORK( VECSIDX+( ISTART-1 )*3+1 ),
     $                             3*( ISTOP-ISTART+1 ) )
                  END IF
               END IF
   80       CONTINUE
*
*        Now do column transforms and finish work
*
            DO 90 KI = 1, IBULGE
*
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
*
               IF( MOD( ISTART-1, HBL ).EQ.HBL-2 ) THEN
                  IF( ( RIGHT.EQ.ICURCOL( KI ) ) .AND.
     $                ( NPROW.GT.1 ) .AND. ( ISTART.LE.ISTOP ) ) THEN
                     JJ = MOD( ICURROW( KI )+NPROW-1, NPROW )
                     IF( MYROW.NE.JJ ) THEN
                        CALL SGEBR2D( CONTXT, 'COL', ' ',
     $                                3*( ISTOP-ISTART+1 ), 1,
     $                                WORK( VECSIDX+( ISTART-1 )*3+1 ),
     $                                3*( ISTOP-ISTART+1 ), JJ, MYCOL )
                     ELSE
                        CALL SGEBS2D( CONTXT, 'COL', ' ',
     $                                3*( ISTOP-ISTART+1 ), 1,
     $                                WORK( VECSIDX+( ISTART-1 )*3+1 ),
     $                                3*( ISTOP-ISTART+1 ) )
                     END IF
                  END IF
               END IF
*
               IF( ( MYCOL.EQ.ICURCOL( KI ) ) .AND. ( NPROW.GT.1 ) .AND.
     $             ( ISTART.LE.ISTOP ) ) THEN
                  IF( MYROW.NE.ICURROW( KI ) ) THEN
                     CALL SGEBR2D( CONTXT, 'COL', ' ',
     $                             3*( ISTOP-ISTART+1 ), 1,
     $                             WORK( VECSIDX+( ISTART-1 )*3+1 ),
     $                             3*( ISTOP-ISTART+1 ), ICURROW( KI ),
     $                             MYCOL )
                  ELSE
                     CALL SGEBS2D( CONTXT, 'COL', ' ',
     $                             3*( ISTOP-ISTART+1 ), 1,
     $                             WORK( VECSIDX+( ISTART-1 )*3+1 ),
     $                             3*( ISTOP-ISTART+1 ) )
                  END IF
               END IF
   90       CONTINUE
*
*        Now do make up work to have things in block fashion
*
            DO 150 KI = 1, IBULGE
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
*
               MODKM1 = MOD( ISTART-1, HBL )
               IF( ( MYROW.EQ.ICURROW( KI ) ) .AND.
     $             ( MYCOL.EQ.ICURCOL( KI ) ) .AND.
     $             ( MODKM1.EQ.HBL-2 ) .AND. ( ISTART.LT.I-1 ) ) THEN
                  K = ISTART
*
*              Catch up on column & border work
*
                  NR = MIN( 3, I-K+1 )
                  V2 = WORK( VECSIDX+( K-1 )*3+1 )
                  V3 = WORK( VECSIDX+( K-1 )*3+2 )
                  T1 = WORK( VECSIDX+( K-1 )*3+3 )
                  IF( NR.EQ.3 ) THEN
*
*                 Do some work so next step is ready...
*
*                 V3 = VCOPY( 3 )
                     T2 = T1*V2
                     T3 = T1*V3
                     ITMP1 = MIN( 6, I2+2-K )
                     ITMP2 = MAX( I1-K+2, 1 )
                     DO 100 J = 2, ITMP1
                        SUM = SMALLA( 2, J, KI ) +
     $                        V2*SMALLA( 3, J, KI ) +
     $                        V3*SMALLA( 4, J, KI )
                        SMALLA( 2, J, KI ) = SMALLA( 2, J, KI ) - SUM*T1
                        SMALLA( 3, J, KI ) = SMALLA( 3, J, KI ) - SUM*T2
                        SMALLA( 4, J, KI ) = SMALLA( 4, J, KI ) - SUM*T3
  100                CONTINUE
                     DO 110 J = ITMP2, 5
                        SUM = SMALLA( J, 2, KI ) +
     $                        V2*SMALLA( J, 3, KI ) +
     $                        V3*SMALLA( J, 4, KI )
                        SMALLA( J, 2, KI ) = SMALLA( J, 2, KI ) - SUM*T1
                        SMALLA( J, 3, KI ) = SMALLA( J, 3, KI ) - SUM*T2
                        SMALLA( J, 4, KI ) = SMALLA( J, 4, KI ) - SUM*T3
  110                CONTINUE
                  END IF
               END IF
*
               IF( ( MOD( ISTART-1, HBL ).EQ.HBL-1 ) .AND.
     $             ( ISTART.LE.ISTOP ) .AND.
     $             ( MYROW.EQ.ICURROW( KI ) ) .AND.
     $             ( MYCOL.EQ.ICURCOL( KI ) ) ) THEN
                  K = ISTOP
*
*              Catch up on column & border work
*
                  NR = MIN( 3, I-K+1 )
                  V2 = WORK( VECSIDX+( K-1 )*3+1 )
                  V3 = WORK( VECSIDX+( K-1 )*3+2 )
                  T1 = WORK( VECSIDX+( K-1 )*3+3 )
                  IF( NR.EQ.3 ) THEN
*
*                 Do some work so next step is ready...
*
*                 V3 = VCOPY( 3 )
                     T2 = T1*V2
                     T3 = T1*V3
                     ITMP1 = MIN( 6, I2-K+3 )
                     ITMP2 = MAX( I1-K+3, 1 )
                     DO 120 J = 3, ITMP1
                        SUM = SMALLA( 3, J, KI ) +
     $                        V2*SMALLA( 4, J, KI ) +
     $                        V3*SMALLA( 5, J, KI )
                        SMALLA( 3, J, KI ) = SMALLA( 3, J, KI ) - SUM*T1
                        SMALLA( 4, J, KI ) = SMALLA( 4, J, KI ) - SUM*T2
                        SMALLA( 5, J, KI ) = SMALLA( 5, J, KI ) - SUM*T3
  120                CONTINUE
                     DO 130 J = ITMP2, 6
                        SUM = SMALLA( J, 3, KI ) +
     $                        V2*SMALLA( J, 4, KI ) +
     $                        V3*SMALLA( J, 5, KI )
                        SMALLA( J, 3, KI ) = SMALLA( J, 3, KI ) - SUM*T1
                        SMALLA( J, 4, KI ) = SMALLA( J, 4, KI ) - SUM*T2
                        SMALLA( J, 5, KI ) = SMALLA( J, 5, KI ) - SUM*T3
  130                CONTINUE
                  END IF
               END IF
*
               MODKM1 = MOD( ISTART-1, HBL )
               IF( ( MYROW.EQ.ICURROW( KI ) ) .AND.
     $             ( MYCOL.EQ.ICURCOL( KI ) ) .AND.
     $             ( ( ( MODKM1.EQ.HBL-2 ) .AND. ( ISTART.EQ.I-
     $             1 ) ) .OR. ( ( MODKM1.LT.HBL-2 ) .AND. ( ISTART.LE.I-
     $             1 ) ) ) ) THEN
*
*           (IROW1,ICOL1) is (I,J)-coordinates of H(ISTART,ISTART)
*
                  IROW1 = KROW( KI )
                  ICOL1 = LOCALK2( KI )
                  DO 140 K = ISTART, ISTOP
*
*              Catch up on column & border work
*
                     NR = MIN( 3, I-K+1 )
                     V2 = WORK( VECSIDX+( K-1 )*3+1 )
                     V3 = WORK( VECSIDX+( K-1 )*3+2 )
                     T1 = WORK( VECSIDX+( K-1 )*3+3 )
                     IF( K.LT.ISTOP ) THEN
*
*                 Do some work so next step is ready...
*
                        T2 = T1*V2
                        T3 = T1*V3
                        CALL SLAREF( 'Col', A, LDA, .FALSE., Z, LDZ,
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
                           T2 = T1*V2
                           T3 = T1*V3
                           CALL SLAREF( 'Row', A, LDA, .FALSE., Z, LDZ,
     $                                  .FALSE., IROW1, IROW1, ISTART,
     $                                  ISTOP, ICOL1, MIN( MIN( K2( KI )
     $                                  +1, I-1 ), I2 )-K+ICOL1, LILOZ,
     $                                  LIHIZ, WORK( VECSIDX+1 ), V2,
     $                                  V3, T1, T2, T3 )
                        END IF
                     END IF
  140             CONTINUE
               END IF
*
*           Send SMALLA back again.
*
               K = ISTART
               MODKM1 = MOD( K-1, HBL )
               IF( ( MODKM1.GE.HBL-2 ) .AND. ( K.LE.I-1 ) ) THEN
                  IF( ( MODKM1.EQ.HBL-2 ) .AND. ( K.LT.I-1 ) ) THEN
*
*                 Copy 6 elements from global A(K-1:K+4,K-1:K+4)
*
                     CALL INFOG2L( K+2, K+2, DESCA, NPROW, NPCOL, MYROW,
     $                             MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                     CALL PSLACP3( MIN( 6, N-K+2 ), K-1, A, DESCA,
     $                             SMALLA( 1, 1, KI ), 6, ITMP1, ITMP2,
     $                             1 )
*
                  END IF
                  IF( MODKM1.EQ.HBL-1 ) THEN
*
*                 Copy 6 elements from global A(K-2:K+3,K-2:K+3)
*
                     CALL INFOG2L( K+1, K+1, DESCA, NPROW, NPCOL, MYROW,
     $                             MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                     CALL PSLACP3( MIN( 6, N-K+3 ), K-2, A, DESCA,
     $                             SMALLA( 1, 1, KI ), 6, ITMP1, ITMP2,
     $                             1 )
                  END IF
               END IF
*
  150       CONTINUE
*
*        Now start major set of block ROW reflections
*
            DO 160 KI = 1, IBULGE
               IF( ( MYROW.NE.ICURROW( KI ) ) .AND.
     $             ( DOWN.NE.ICURROW( KI ) ) )GO TO 160
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
*
               IF( ( ISTOP.GT.ISTART ) .AND.
     $             ( MOD( ISTART-1, HBL ).LT.HBL-2 ) .AND.
     $             ( ICURROW( KI ).EQ.MYROW ) ) THEN
                  IROW1 = MIN( K2( KI )+1, I-1 ) + 1
                  CALL INFOG1L( IROW1, HBL, NPCOL, MYCOL, 0, ITMP1,
     $                          ITMP2 )
                  ITMP2 = NUMROC( I2, HBL, MYCOL, 0, NPCOL )
                  II = KROW( KI )
                  CALL SLAREF( 'Row', A, LDA, WANTZ, Z, LDZ, .TRUE., II,
     $                         II, ISTART, ISTOP, ITMP1, ITMP2, LILOZ,
     $                         LIHIZ, WORK( VECSIDX+1 ), V2, V3, T1, T2,
     $                         T3 )
               END IF
  160       CONTINUE
*
            DO 180 KI = 1, IBULGE
               IF( KROW( KI ).GT.KP2ROW( KI ) )
     $            GO TO 180
               IF( ( MYROW.NE.ICURROW( KI ) ) .AND.
     $             ( DOWN.NE.ICURROW( KI ) ) )GO TO 180
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
               IF( ( ISTART.EQ.ISTOP ) .OR.
     $             ( MOD( ISTART-1, HBL ).GE.HBL-2 ) .OR.
     $             ( ICURROW( KI ).NE.MYROW ) ) THEN
                  DO 170 K = ISTART, ISTOP
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
                        CALL INFOG1L( ITMP1, HBL, NPCOL, MYCOL, 0,
     $                                ICOL1, ICOL2 )
                        ICOL2 = NUMROC( I2, HBL, MYCOL, 0, NPCOL )
                        IF( ( MOD( K-1, HBL ).LT.HBL-2 ) .OR.
     $                      ( NPROW.EQ.1 ) ) THEN
                           T2 = T1*V2
                           T3 = T1*V3
                           CALL SLAREF( 'Row', A, LDA, WANTZ, Z, LDZ,
     $                                  .FALSE., IROW1, IROW1, ISTART,
     $                                  ISTOP, ICOL1, ICOL2, LILOZ,
     $                                  LIHIZ, WORK( VECSIDX+1 ), V2,
     $                                  V3, T1, T2, T3 )
                        END IF
                        IF( ( MOD( K-1, HBL ).EQ.HBL-2 ) .AND.
     $                      ( NPROW.GT.1 ) ) THEN
                           IF( IROW1.EQ.IROW2 ) THEN
                              CALL SGESD2D( CONTXT, 1, ICOL2-ICOL1+1,
     $                                      A( ( ICOL1-1 )*LDA+IROW2 ),
     $                                      LDA, UP, MYCOL )
                           END IF
                        END IF
                        IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                      ( NPROW.GT.1 ) ) THEN
                           IF( IROW1.EQ.IROW2 ) THEN
                              CALL SGESD2D( CONTXT, 1, ICOL2-ICOL1+1,
     $                                      A( ( ICOL1-1 )*LDA+IROW1 ),
     $                                      LDA, DOWN, MYCOL )
                           END IF
                        END IF
                     END IF
  170             CONTINUE
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
                        IROW1 = KROW( KI ) + K - ISTART
                        IROW2 = KP2ROW( KI ) + K - ISTART
                        CALL INFOG1L( ITMP1, HBL, NPCOL, MYCOL, 0,
     $                                ICOL1, ICOL2 )
                        ICOL2 = NUMROC( I2, HBL, MYCOL, 0, NPCOL )
                        IF( ( MOD( K-1, HBL ).EQ.HBL-2 ) .AND.
     $                      ( NPROW.GT.1 ) ) THEN
                           IF( IROW1.NE.IROW2 ) THEN
                              CALL SGERV2D( CONTXT, 1, ICOL2-ICOL1+1,
     $                                      WORK( IRBUF+1 ), 1, DOWN,
     $                                      MYCOL )
                              T2 = T1*V2
                              T3 = T1*V3
                              DO 190 J = ICOL1, ICOL2
                                 SUM = A( ( J-1 )*LDA+IROW1 ) +
     $                                 V2*A( ( J-1 )*LDA+IROW1+1 ) +
     $                                 V3*WORK( IRBUF+J-ICOL1+1 )
                                 A( ( J-1 )*LDA+IROW1 ) = A( ( J-1 )*
     $                              LDA+IROW1 ) - SUM*T1
                                 A( ( J-1 )*LDA+IROW1+1 ) = A( ( J-1 )*
     $                              LDA+IROW1+1 ) - SUM*T2
                                 WORK( IRBUF+J-ICOL1+1 ) = WORK( IRBUF+
     $                              J-ICOL1+1 ) - SUM*T3
  190                         CONTINUE
                              CALL SGESD2D( CONTXT, 1, ICOL2-ICOL1+1,
     $                                      WORK( IRBUF+1 ), 1, DOWN,
     $                                      MYCOL )
                           END IF
                        END IF
                        IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                      ( NPROW.GT.1 ) ) THEN
                           IF( IROW1.NE.IROW2 ) THEN
                              CALL SGERV2D( CONTXT, 1, ICOL2-ICOL1+1,
     $                                      WORK( IRBUF+1 ), 1, UP,
     $                                      MYCOL )
                              T2 = T1*V2
                              T3 = T1*V3
                              DO 200 J = ICOL1, ICOL2
                                 SUM = WORK( IRBUF+J-ICOL1+1 ) +
     $                                 V2*A( ( J-1 )*LDA+IROW1 ) +
     $                                 V3*A( ( J-1 )*LDA+IROW1+1 )
                                 WORK( IRBUF+J-ICOL1+1 ) = WORK( IRBUF+
     $                              J-ICOL1+1 ) - SUM*T1
                                 A( ( J-1 )*LDA+IROW1 ) = A( ( J-1 )*
     $                              LDA+IROW1 ) - SUM*T2
                                 A( ( J-1 )*LDA+IROW1+1 ) = A( ( J-1 )*
     $                              LDA+IROW1+1 ) - SUM*T3
  200                         CONTINUE
                              CALL SGESD2D( CONTXT, 1, ICOL2-ICOL1+1,
     $                                      WORK( IRBUF+1 ), 1, UP,
     $                                      MYCOL )
                           END IF
                        END IF
                     END IF
  210             CONTINUE
               END IF
  220       CONTINUE
*
            DO 240 KI = 1, IBULGE
               IF( KROW( KI ).GT.KP2ROW( KI ) )
     $            GO TO 240
               IF( ( MYROW.NE.ICURROW( KI ) ) .AND.
     $             ( DOWN.NE.ICURROW( KI ) ) )GO TO 240
               ISTART = MAX( K1( KI ), M )
               ISTOP = MIN( K2( KI ), I-1 )
               IF( ( ISTART.EQ.ISTOP ) .OR.
     $             ( MOD( ISTART-1, HBL ).GE.HBL-2 ) .OR.
     $             ( ICURROW( KI ).NE.MYROW ) ) THEN
                  DO 230 K = ISTART, ISTOP
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
                        IROW1 = KROW( KI ) + K - ISTART
                        IROW2 = KP2ROW( KI ) + K - ISTART
                        CALL INFOG1L( ITMP1, HBL, NPCOL, MYCOL, 0,
     $                                ICOL1, ICOL2 )
                        ICOL2 = NUMROC( I2, HBL, MYCOL, 0, NPCOL )
                        IF( ( MOD( K-1, HBL ).EQ.HBL-2 ) .AND.
     $                      ( NPROW.GT.1 ) ) THEN
                           IF( IROW1.EQ.IROW2 ) THEN
                              CALL SGERV2D( CONTXT, 1, ICOL2-ICOL1+1,
     $                                      A( ( ICOL1-1 )*LDA+IROW2 ),
     $                                      LDA, UP, MYCOL )
                           END IF
                        END IF
                        IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                      ( NPROW.GT.1 ) ) THEN
                           IF( IROW1.EQ.IROW2 ) THEN
                              CALL SGERV2D( CONTXT, 1, ICOL2-ICOL1+1,
     $                                      A( ( ICOL1-1 )*LDA+IROW1 ),
     $                                      LDA, DOWN, MYCOL )
                           END IF
                        END IF
                     END IF
  230             CONTINUE
               END IF
  240       CONTINUE
  250       CONTINUE
*
*        Now start major set of block COL reflections
*
            DO 260 KI = 1, IBULGE
               IF( ( MYCOL.NE.ICURCOL( KI ) ) .AND.
     $             ( RIGHT.NE.ICURCOL( KI ) ) )GO TO 260
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
                  CALL INFOG1L( I1, HBL, NPROW, MYROW, 0, IROW1, IROW2 )
                  IROW2 = NUMROC( ITMP1, HBL, MYROW, 0, NPROW )
                  IF( IROW1.LE.IROW2 ) THEN
                     ITMP2 = IROW2
                  ELSE
                     ITMP2 = -1
                  END IF
                  CALL SLAREF( 'Col', A, LDA, WANTZ, Z, LDZ, .TRUE.,
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
                        CALL INFOG1L( ITMP1+1, HBL, NPROW, MYROW, 0,
     $                                IROW1, IROW2 )
                        IROW2 = NUMROC( MIN( K+3, I ), HBL, MYROW, 0,
     $                          NPROW )
                     END IF
                     V2 = WORK( VECSIDX+( K-1 )*3+1 )
                     V3 = WORK( VECSIDX+( K-1 )*3+2 )
                     T1 = WORK( VECSIDX+( K-1 )*3+3 )
                     T2 = T1*V2
                     T3 = T1*V3
                     ICOL1 = KCOL( KI ) + ISTOP - ISTART
                     CALL SLAREF( 'Col', A, LDA, .FALSE., Z, LDZ,
     $                            .FALSE., ICOL1, ICOL1, ISTART, ISTOP,
     $                            IROW1, IROW2, LILOZ, LIHIZ,
     $                            WORK( VECSIDX+1 ), V2, V3, T1, T2,
     $                            T3 )
                  END IF
               END IF
  260       CONTINUE
*
            DO 320 KI = 1, IBULGE
               IF( KCOL( KI ).GT.KP2COL( KI ) )
     $            GO TO 320
               IF( ( MYCOL.NE.ICURCOL( KI ) ) .AND.
     $             ( RIGHT.NE.ICURCOL( KI ) ) )GO TO 320
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
*
               DO 310 K = ISTART, ISTOP
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
                     ICOL1 = KCOL( KI ) + K - ISTART
                     ICOL2 = KP2COL( KI ) + K - ISTART
                     CALL INFOG1L( I1, HBL, NPROW, MYROW, 0, IROW1,
     $                             IROW2 )
                     IROW2 = NUMROC( ITMP1, HBL, MYROW, 0, NPROW )
                     IF( ( MOD( K-1, HBL ).EQ.HBL-2 ) .AND.
     $                   ( NPCOL.GT.1 ) ) THEN
                        IF( ICOL1.EQ.ICOL2 ) THEN
                           CALL SGESD2D( CONTXT, IROW2-IROW1+1, 1,
     $                                   A( ( ICOL1-1 )*LDA+IROW1 ),
     $                                   LDA, MYROW, LEFT )
                           CALL SGERV2D( CONTXT, IROW2-IROW1+1, 1,
     $                                   A( ( ICOL1-1 )*LDA+IROW1 ),
     $                                   LDA, MYROW, LEFT )
                        ELSE
                           CALL SGERV2D( CONTXT, IROW2-IROW1+1, 1,
     $                                   WORK( ICBUF+1 ), IROW2-IROW1+1,
     $                                   MYROW, RIGHT )
                           T2 = T1*V2
                           T3 = T1*V3
                           DO 270 J = IROW1, IROW2
                              SUM = A( ( ICOL1-1 )*LDA+J ) +
     $                              V2*A( ICOL1*LDA+J ) +
     $                              V3*WORK( ICBUF+J-IROW1+1 )
                              A( ( ICOL1-1 )*LDA+J ) = A( ( ICOL1-1 )*
     $                           LDA+J ) - SUM*T1
                              A( ICOL1*LDA+J ) = A( ICOL1*LDA+J ) -
     $                                           SUM*T2
                              WORK( ICBUF+J-IROW1+1 ) = WORK( ICBUF+J-
     $                           IROW1+1 ) - SUM*T3
  270                      CONTINUE
                           CALL SGESD2D( CONTXT, IROW2-IROW1+1, 1,
     $                                   WORK( ICBUF+1 ), IROW2-IROW1+1,
     $                                   MYROW, RIGHT )
                        END IF
                     END IF
                     IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                   ( NPCOL.GT.1 ) ) THEN
                        IF( ICOL1.EQ.ICOL2 ) THEN
                           CALL SGESD2D( CONTXT, IROW2-IROW1+1, 1,
     $                                   A( ( ICOL1-1 )*LDA+IROW1 ),
     $                                   LDA, MYROW, RIGHT )
                           CALL SGERV2D( CONTXT, IROW2-IROW1+1, 1,
     $                                   A( ( ICOL1-1 )*LDA+IROW1 ),
     $                                   LDA, MYROW, RIGHT )
                        ELSE
                           CALL SGERV2D( CONTXT, IROW2-IROW1+1, 1,
     $                                   WORK( ICBUF+1 ), IROW2-IROW1+1,
     $                                   MYROW, LEFT )
                           T2 = T1*V2
                           T3 = T1*V3
                           DO 280 J = IROW1, IROW2
                              SUM = WORK( ICBUF+J-IROW1+1 ) +
     $                              V2*A( ( ICOL1-1 )*LDA+J ) +
     $                              V3*A( ICOL1*LDA+J )
                              WORK( ICBUF+J-IROW1+1 ) = WORK( ICBUF+J-
     $                           IROW1+1 ) - SUM*T1
                              A( ( ICOL1-1 )*LDA+J ) = A( ( ICOL1-1 )*
     $                           LDA+J ) - SUM*T2
                              A( ICOL1*LDA+J ) = A( ICOL1*LDA+J ) -
     $                                           SUM*T3
  280                      CONTINUE
                           CALL SGESD2D( CONTXT, IROW2-IROW1+1, 1,
     $                                   WORK( ICBUF+1 ), IROW2-IROW1+1,
     $                                   MYROW, LEFT )
                        END IF
                     END IF
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
                              CALL SGESD2D( CONTXT, IROW2-IROW1+1, 1,
     $                                      Z( ( ICOL1-1 )*LDZ+IROW1 ),
     $                                      LDZ, MYROW, LEFT )
                              CALL SGERV2D( CONTXT, IROW2-IROW1+1, 1,
     $                                      Z( ( ICOL1-1 )*LDZ+IROW1 ),
     $                                      LDZ, MYROW, LEFT )
                           ELSE
                              CALL SGERV2D( CONTXT, IROW2-IROW1+1, 1,
     $                                      WORK( ICBUF+1 ),
     $                                      IROW2-IROW1+1, MYROW,
     $                                      RIGHT )
                              T2 = T1*V2
                              T3 = T1*V3
                              ICOL1 = ( ICOL1-1 )*LDZ
                              DO 290 J = IROW1, IROW2
                                 SUM = Z( ICOL1+J ) +
     $                                 V2*Z( ICOL1+J+LDZ ) +
     $                                 V3*WORK( ICBUF+J-IROW1+1 )
                                 Z( J+ICOL1 ) = Z( J+ICOL1 ) - SUM*T1
                                 Z( J+ICOL1+LDZ ) = Z( J+ICOL1+LDZ ) -
     $                                              SUM*T2
                                 WORK( ICBUF+J-IROW1+1 ) = WORK( ICBUF+
     $                              J-IROW1+1 ) - SUM*T3
  290                         CONTINUE
                              CALL SGESD2D( CONTXT, IROW2-IROW1+1, 1,
     $                                      WORK( ICBUF+1 ),
     $                                      IROW2-IROW1+1, MYROW,
     $                                      RIGHT )
                           END IF
                        END IF
                        IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                           IF( ICOL1.EQ.ICOL2 ) THEN
                              CALL SGESD2D( CONTXT, IROW2-IROW1+1, 1,
     $                                      Z( ( ICOL1-1 )*LDZ+IROW1 ),
     $                                      LDZ, MYROW, RIGHT )
                              CALL SGERV2D( CONTXT, IROW2-IROW1+1, 1,
     $                                      Z( ( ICOL1-1 )*LDZ+IROW1 ),
     $                                      LDZ, MYROW, RIGHT )
                           ELSE
                              CALL SGERV2D( CONTXT, IROW2-IROW1+1, 1,
     $                                      WORK( ICBUF+1 ),
     $                                      IROW2-IROW1+1, MYROW, LEFT )
                              T2 = T1*V2
                              T3 = T1*V3
                              ICOL1 = ( ICOL1-1 )*LDZ
                              DO 300 J = IROW1, IROW2
                                 SUM = WORK( ICBUF+J-IROW1+1 ) +
     $                                 V2*Z( J+ICOL1 ) +
     $                                 V3*Z( J+ICOL1+LDZ )
                                 WORK( ICBUF+J-IROW1+1 ) = WORK( ICBUF+
     $                              J-IROW1+1 ) - SUM*T1
                                 Z( J+ICOL1 ) = Z( J+ICOL1 ) - SUM*T2
                                 Z( J+ICOL1+LDZ ) = Z( J+ICOL1+LDZ ) -
     $                                              SUM*T3
  300                         CONTINUE
                              CALL SGESD2D( CONTXT, IROW2-IROW1+1, 1,
     $                                      WORK( ICBUF+1 ),
     $                                      IROW2-IROW1+1, MYROW, LEFT )
                           END IF
                        END IF
                     END IF
                     IF( ICURCOL( KI ).EQ.MYCOL ) THEN
                        IF( ( ISPEC.EQ.0 ) .OR. ( NPCOL.EQ.1 ) ) THEN
                           LOCALK2( KI ) = LOCALK2( KI ) + 1
                        END IF
                     ELSE
                        IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                      ( ICURCOL( KI ).EQ.RIGHT ) ) THEN
                           IF( K.GT.M ) THEN
                              LOCALK2( KI ) = LOCALK2( KI ) + 2
                           ELSE
                              LOCALK2( KI ) = LOCALK2( KI ) + 1
                           END IF
                        END IF
                        IF( ( MOD( K-1, HBL ).EQ.HBL-2 ) .AND.
     $                      ( I-K.EQ.2 ) .AND. ( ICURCOL( KI ).EQ.
     $                      RIGHT ) ) THEN
                           LOCALK2( KI ) = LOCALK2( KI ) + 2
                        END IF
                     END IF
                  END IF
  310          CONTINUE
  320       CONTINUE
*
*        Column work done
*
  330       CONTINUE
*
*        Now do NR=2 work
*
            DO 410 KI = 1, IBULGE
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
*
               DO 400 K = ISTART, ISTOP
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
                     CALL INFOG1L( K, HBL, NPCOL, MYCOL, 0, LILOH,
     $                             LIHIH )
                     LIHIH = NUMROC( I2, HBL, MYCOL, 0, NPCOL )
                     CALL INFOG1L( 1, HBL, NPROW, MYROW, 0, ITMP2,
     $                             ITMP1 )
                     ITMP1 = NUMROC( K+1, HBL, MYROW, 0, NPROW )
                     IF( ICURROW( KI ).EQ.MYROW ) THEN
                        IF( ( ISPEC.EQ.0 ) .OR. ( NPROW.EQ.1 ) .OR.
     $                      ( MOD( K-1, HBL ).EQ.HBL-2 ) ) THEN
                           ITMP1 = ITMP1 - 1
                           DO 340 J = ( LILOH-1 )*LDA,
     $                             ( LIHIH-1 )*LDA, LDA
                              SUM = A( ITMP1+J ) + V2*A( ITMP1+1+J )
                              A( ITMP1+J ) = A( ITMP1+J ) - SUM*T1
                              A( ITMP1+1+J ) = A( ITMP1+1+J ) - SUM*T2
  340                      CONTINUE
                        ELSE
                           IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                              CALL SGERV2D( CONTXT, 1, LIHIH-LILOH+1,
     $                                      WORK( IRBUF+1 ), 1, UP,
     $                                      MYCOL )
                              DO 350 J = LILOH, LIHIH
                                 SUM = WORK( IRBUF+J-LILOH+1 ) +
     $                                 V2*A( ( J-1 )*LDA+ITMP1 )
                                 WORK( IRBUF+J-LILOH+1 ) = WORK( IRBUF+
     $                              J-LILOH+1 ) - SUM*T1
                                 A( ( J-1 )*LDA+ITMP1 ) = A( ( J-1 )*
     $                              LDA+ITMP1 ) - SUM*T2
  350                         CONTINUE
                              CALL SGESD2D( CONTXT, 1, LIHIH-LILOH+1,
     $                                      WORK( IRBUF+1 ), 1, UP,
     $                                      MYCOL )
                           END IF
                        END IF
                     ELSE
                        IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                      ( ICURROW( KI ).EQ.DOWN ) ) THEN
                           CALL SGESD2D( CONTXT, 1, LIHIH-LILOH+1,
     $                                   A( ( LILOH-1 )*LDA+ITMP1 ),
     $                                   LDA, DOWN, MYCOL )
                           CALL SGERV2D( CONTXT, 1, LIHIH-LILOH+1,
     $                                   A( ( LILOH-1 )*LDA+ITMP1 ),
     $                                   LDA, DOWN, MYCOL )
                        END IF
                     END IF
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to MIN(K+3,I).
*
                     CALL INFOG1L( I1, HBL, NPROW, MYROW, 0, LILOH,
     $                             LIHIH )
                     LIHIH = NUMROC( I, HBL, MYROW, 0, NPROW )
*
                     IF( ICURCOL( KI ).EQ.MYCOL ) THEN
*                 LOCAL A(LILOZ:LIHIZ,LOCALK2:LOCALK2+2)
                        IF( ( ISPEC.EQ.0 ) .OR. ( NPCOL.EQ.1 ) .OR.
     $                      ( MOD( K-1, HBL ).EQ.HBL-2 ) ) THEN
                           CALL INFOG1L( K, HBL, NPCOL, MYCOL, 0, ITMP1,
     $                                   ITMP2 )
                           ITMP2 = NUMROC( K+1, HBL, MYCOL, 0, NPCOL )
                           DO 360 J = LILOH, LIHIH
                              SUM = A( ( ITMP1-1 )*LDA+J ) +
     $                              V2*A( ITMP1*LDA+J )
                              A( ( ITMP1-1 )*LDA+J ) = A( ( ITMP1-1 )*
     $                           LDA+J ) - SUM*T1
                              A( ITMP1*LDA+J ) = A( ITMP1*LDA+J ) -
     $                                           SUM*T2
  360                      CONTINUE
                        ELSE
                           ITMP1 = LOCALK2( KI )
                           IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                              CALL SGERV2D( CONTXT, LIHIH-LILOH+1, 1,
     $                                      WORK( ICBUF+1 ),
     $                                      LIHIH-LILOH+1, MYROW, LEFT )
                              DO 370 J = LILOH, LIHIH
                                 SUM = WORK( ICBUF+J ) +
     $                                 V2*A( ( ITMP1-1 )*LDA+J )
                                 WORK( ICBUF+J ) = WORK( ICBUF+J ) -
     $                                             SUM*T1
                                 A( ( ITMP1-1 )*LDA+J )
     $                              = A( ( ITMP1-1 )*LDA+J ) - SUM*T2
  370                         CONTINUE
                              CALL SGESD2D( CONTXT, LIHIH-LILOH+1, 1,
     $                                      WORK( ICBUF+1 ),
     $                                      LIHIH-LILOH+1, MYROW, LEFT )
                           END IF
                        END IF
                     ELSE
                        IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                      ( ICURCOL( KI ).EQ.RIGHT ) ) THEN
                           ITMP1 = KCOL( KI )
                           CALL SGESD2D( CONTXT, LIHIH-LILOH+1, 1,
     $                                   A( ( ITMP1-1 )*LDA+LILOH ),
     $                                   LDA, MYROW, RIGHT )
                           CALL INFOG1L( K, HBL, NPCOL, MYCOL, 0, ITMP1,
     $                                   ITMP2 )
                           ITMP2 = NUMROC( K+1, HBL, MYCOL, 0, NPCOL )
                           CALL SGERV2D( CONTXT, LIHIH-LILOH+1, 1,
     $                                   A( ( ITMP1-1 )*LDA+LILOH ),
     $                                   LDA, MYROW, RIGHT )
                        END IF
                     END IF
*
                     IF( WANTZ ) THEN
*
*                 Accumulate transformations in the matrix Z
*
                        IF( ICURCOL( KI ).EQ.MYCOL ) THEN
*                    LOCAL Z(LILOZ:LIHIZ,LOCALK2:LOCALK2+2)
                           IF( ( ISPEC.EQ.0 ) .OR. ( NPCOL.EQ.1 ) .OR.
     $                         ( MOD( K-1, HBL ).EQ.HBL-2 ) ) THEN
                              ITMP1 = KCOL( KI ) + K - ISTART
                              ITMP1 = ( ITMP1-1 )*LDZ
                              DO 380 J = LILOZ, LIHIZ
                                 SUM = Z( J+ITMP1 ) +
     $                                 V2*Z( J+ITMP1+LDZ )
                                 Z( J+ITMP1 ) = Z( J+ITMP1 ) - SUM*T1
                                 Z( J+ITMP1+LDZ ) = Z( J+ITMP1+LDZ ) -
     $                                              SUM*T2
  380                         CONTINUE
                              LOCALK2( KI ) = LOCALK2( KI ) + 1
                           ELSE
                              ITMP1 = LOCALK2( KI )
*                       IF WE ACTUALLY OWN COLUMN K
                              IF( MOD( K-1, HBL ).EQ.HBL-1 ) THEN
                                 CALL SGERV2D( CONTXT, LIHIZ-LILOZ+1, 1,
     $                                         WORK( ICBUF+1 ), LDZ,
     $                                         MYROW, LEFT )
                                 ITMP1 = ( ITMP1-1 )*LDZ
                                 DO 390 J = LILOZ, LIHIZ
                                    SUM = WORK( ICBUF+J ) +
     $                                    V2*Z( J+ITMP1 )
                                    WORK( ICBUF+J ) = WORK( ICBUF+J ) -
     $                                 SUM*T1
                                    Z( J+ITMP1 ) = Z( J+ITMP1 ) - SUM*T2
  390                            CONTINUE
                                 CALL SGESD2D( CONTXT, LIHIZ-LILOZ+1, 1,
     $                                         WORK( ICBUF+1 ), LDZ,
     $                                         MYROW, LEFT )
                                 LOCALK2( KI ) = LOCALK2( KI ) + 1
                              END IF
                           END IF
                        ELSE
*
*                    NO WORK BUT NEED TO UPDATE ANYWAY????
*
                           IF( ( MOD( K-1, HBL ).EQ.HBL-1 ) .AND.
     $                         ( ICURCOL( KI ).EQ.RIGHT ) ) THEN
                              ITMP1 = KCOL( KI )
                              ITMP1 = ( ITMP1-1 )*LDZ
                              CALL SGESD2D( CONTXT, LIHIZ-LILOZ+1, 1,
     $                                      Z( LILOZ+ITMP1 ), LDZ,
     $                                      MYROW, RIGHT )
                              CALL SGERV2D( CONTXT, LIHIZ-LILOZ+1, 1,
     $                                      Z( LILOZ+ITMP1 ), LDZ,
     $                                      MYROW, RIGHT )
                              LOCALK2( KI ) = LOCALK2( KI ) + 1
                           END IF
                        END IF
                     END IF
                  END IF
  400          CONTINUE
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
                  CALL INFOG1L( K2( KI )+1, HBL, NPROW, MYROW, 0,
     $                          KROW( KI ), ITMP2 )
                  ITMP2 = NUMROC( N, HBL, MYROW, 0, NPROW )
               END IF
               IF( ( MOD( K2( KI ), HBL ).GE.HBL-2 ) .AND.
     $             ( ( MYROW.EQ.ICURROW( KI ) ) .OR. ( UP.EQ.
     $             ICURROW( KI ) ) ) .AND. ( NPROW.GT.1 ) ) THEN
                  CALL INFOG1L( 1, HBL, NPROW, MYROW, 0, ITMP2,
     $                          KP2ROW( KI ) )
                  KP2ROW( KI ) = NUMROC( K2( KI )+3, HBL, MYROW, 0,
     $                           NPROW )
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
                  CALL INFOG1L( K2( KI )+1, HBL, NPCOL, MYCOL, 0,
     $                          KCOL( KI ), ITMP2 )
                  ITMP2 = NUMROC( N, HBL, MYCOL, 0, NPCOL )
               END IF
               IF( ( MOD( K2( KI ), HBL ).GE.HBL-2 ) .AND.
     $             ( ( MYCOL.EQ.ICURCOL( KI ) ) .OR. ( LEFT.EQ.
     $             ICURCOL( KI ) ) ) .AND. ( NPCOL.GT.1 ) ) THEN
                  CALL INFOG1L( 1, HBL, NPCOL, MYCOL, 0, ITMP2,
     $                          KP2COL( KI ) )
                  KP2COL( KI ) = NUMROC( K2( KI )+3, HBL, MYCOL, 0,
     $                           NPCOL )
               END IF
               K1( KI ) = K2( KI ) + 1
               ISTOP = MIN( K1( KI )+ROTN-MOD( K1( KI ), ROTN ), I-2 )
               ISTOP = MIN( ISTOP, K1( KI )+HBL-3-
     $                 MOD( K1( KI )-1, HBL ) )
               ISTOP = MIN( ISTOP, I2-2 )
               ISTOP = MAX( ISTOP, K1( KI ) )
*        ISTOP = MIN( ISTOP , I-1 )
               K2( KI ) = ISTOP
               IF( K1( KI ).EQ.ISTOP ) THEN
                  IF( ( MOD( ISTOP-1, HBL ).EQ.HBL-2 ) .AND.
     $                ( I-ISTOP.GT.1 ) ) THEN
*
*              Next step switches rows & cols
*
                     ICURROW( KI ) = MOD( ICURROW( KI )+1, NPROW )
                     ICURCOL( KI ) = MOD( ICURCOL( KI )+1, NPCOL )
                  END IF
               END IF
  410       CONTINUE
            IF( K2( IBULGE ).LE.I-1 )
     $         GO TO 40
         END IF
*
  420 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      INFO = I
      RETURN
*
  430 CONTINUE
*
      IF( L.EQ.I ) THEN
*
*        H(I,I-1) is negligible: one eigenvalue has converged.
*
         CALL INFOG2L( I, I, DESCA, NPROW, NPCOL, MYROW, MYCOL, IROW,
     $                 ICOL, ITMP1, ITMP2 )
         IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
            WR( I ) = A( ( ICOL-1 )*LDA+IROW )
         ELSE
            WR( I ) = ZERO
         END IF
         WI( I ) = ZERO
      ELSE IF( L.EQ.I-1 ) THEN
*
*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
*
         CALL PSELGET( 'All', ' ', H11, A, L, L, DESCA )
         CALL PSELGET( 'All', ' ', H21, A, I, L, DESCA )
         CALL PSELGET( 'All', ' ', H12, A, L, I, DESCA )
         CALL PSELGET( 'All', ' ', H22, A, I, I, DESCA )
         CALL SLANV2( H11, H12, H21, H22, WR( L ), WI( L ), WR( I ),
     $                WI( I ), CS, SN )
         IF( NODE .NE. 0 ) THEN
            WR( L ) = ZERO
            WR( I ) = ZERO
            WI( L ) = ZERO
            WI( I ) = ZERO
         ENDIF
      ELSE
*
*        Find the eigenvalues in H(L:I,L:I), L < I-1
*
         JBLK = I - L + 1
         IF( JBLK.LE.2*IBLK ) THEN
            CALL PSLACP3( I-L+1, L, A, DESCA, S1, 2*IBLK, 0, 0, 0 )
            CALL SLAHQR( .FALSE., .FALSE., JBLK, 1, JBLK, S1, 2*IBLK,
     $                   WR( L ), WI( L ), 1, JBLK, Z, LDZ, IERR )
            IF( NODE.NE.0 ) THEN
*
*           Erase the eigenvalues
*
               DO 440 K = L, I
                  WR( K ) = ZERO
                  WI( K ) = ZERO
  440          CONTINUE
            END IF
         END IF
      END IF
*
*     Decrement number of remaining iterations, and return to start of
*     the main loop with new value of I.
*
      ITN = ITN - ITS
      IF( M.EQ.L-10 ) THEN
         I = L - 1
      ELSE
         I = M
      END IF
*     I = L - 1
      GO TO 10
*
  450 CONTINUE
      CALL SGSUM2D( CONTXT, 'All', ' ', N, 1, WR, N, -1, -1 )
      CALL SGSUM2D( CONTXT, 'All', ' ', N, 1, WI, N, -1, -1 )
      RETURN
*
*     END OF PSLAHQR
*
      END
