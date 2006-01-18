      SUBROUTINE PDSTEBZ( ICTXT, RANGE, ORDER, N, VL, VU, IL, IU,
     $                    ABSTOL, D, E, M, NSPLIT, W, IBLOCK, ISPLIT,
     $                    WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          ORDER, RANGE
      INTEGER            ICTXT, IL, INFO, IU, LIWORK, LWORK, M, N,
     $                   NSPLIT
      DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDSTEBZ computes the eigenvalues of a symmetric tridiagonal matrix in
*  parallel. The user may ask for all eigenvalues, all eigenvalues in
*  the interval [VL, VU], or the eigenvalues indexed IL through IU. A
*  static partitioning of work is done at the beginning of PDSTEBZ which
*  results in all processes finding an (almost) equal number of
*  eigenvalues.
*
*  NOTE : It is assumed that the user is on an IEEE machine. If the user
*         is not on an IEEE mchine, set the compile time flag NO_IEEE
*         to 1 (in SLmake.inc). The features of IEEE arithmetic that
*         are needed for the "fast" Sturm Count are : (a) infinity
*         arithmetic (b) the sign bit of a single precision floating
*         point number is assumed be in the 32nd bit position
*         (c) the sign of negative zero.
*
*  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
*  Matrix", Report CS41, Computer Science Dept., Stanford
*  University, July 21, 1966.
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle.
*
*  RANGE   (global input) CHARACTER
*          Specifies which eigenvalues are to be found.
*          = 'A': ("All")   all eigenvalues will be found.
*          = 'V': ("Value") all eigenvalues in the interval
*                           [VL, VU] will be found.
*          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
*                           entire matrix) will be found.
*
*  ORDER   (global input) CHARACTER
*          Specifies the order in which the eigenvalues and their block
*          numbers are stored in W and IBLOCK.
*          = 'B': ("By Block") the eigenvalues will be grouped by
*                              split-off block (see IBLOCK, ISPLIT) and
*                              ordered from smallest to largest within
*                              the block.
*          = 'E': ("Entire matrix")
*                              the eigenvalues for the entire matrix
*                              will be ordered from smallest to largest.
*
*  N       (global input) INTEGER
*          The order of the tridiagonal matrix T.  N >= 0.
*
*  VL      (global input) DOUBLE PRECISION
*          If RANGE='V', the lower bound of the interval to be searched
*          for eigenvalues.  Eigenvalues less than VL will not be
*          returned.  Not referenced if RANGE='A' or 'I'.
*
*  VU      (global input) DOUBLE PRECISION
*          If RANGE='V', the upper bound of the interval to be searched
*          for eigenvalues.  Eigenvalues greater than VU will not be
*          returned.  VU must be greater than VL.  Not referenced if
*          RANGE='A' or 'I'.
*
*  IL      (global input) INTEGER
*          If RANGE='I', the index (from smallest to largest) of the
*          smallest eigenvalue to be returned.  IL must be at least 1.
*          Not referenced if RANGE='A' or 'V'.
*
*  IU      (global input) INTEGER
*          If RANGE='I', the index (from smallest to largest) of the
*          largest eigenvalue to be returned.  IU must be at least IL
*          and no greater than N.  Not referenced if RANGE='A' or 'V'.
*
*  ABSTOL  (global input) DOUBLE PRECISION
*          The absolute tolerance for the eigenvalues.  An eigenvalue
*          (or cluster) is considered to be located if it has been
*          determined to lie in an interval whose width is ABSTOL or
*          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
*          will be used, where |T| means the 1-norm of T.
*          Eigenvalues will be computed most accurately when ABSTOL is
*          set to the underflow threshold DLAMCH('U'), not zero.
*          Note : If eigenvectors are desired later by inverse iteration
*          ( PDSTEIN ), ABSTOL should be set to 2*PDLAMCH('S').
*
*  D       (global input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the tridiagonal matrix T.  To
*          avoid overflow, the matrix must be scaled so that its largest
*          entry is no greater than overflow**(1/2) * underflow**(1/4)
*          in absolute value, and for greatest accuracy, it should not
*          be much smaller than that.
*
*  E       (global input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) off-diagonal elements of the tridiagonal matrix T.
*          To avoid overflow, the matrix must be scaled so that its
*          largest entry is no greater than overflow**(1/2) *
*          underflow**(1/4) in absolute value, and for greatest
*          accuracy, it should not be much smaller than that.
*
*  M       (global output) INTEGER
*          The actual number of eigenvalues found. 0 <= M <= N.
*          (See also the description of INFO=2)
*
*  NSPLIT  (global output) INTEGER
*          The number of diagonal blocks in the matrix T.
*          1 <= NSPLIT <= N.
*
*  W       (global output) DOUBLE PRECISION array, dimension (N)
*          On exit, the first M elements of W contain the eigenvalues
*          on all processes.
*
*  IBLOCK  (global output) INTEGER array, dimension (N)
*          At each row/column j where E(j) is zero or small, the
*          matrix T is considered to split into a block diagonal
*          matrix.  On exit IBLOCK(i) specifies which block (from 1
*          to the number of blocks) the eigenvalue W(i) belongs to.
*          NOTE:  in the (theoretically impossible) event that bisection
*          does not converge for some or all eigenvalues, INFO is set
*          to 1 and the ones for which it did not are identified by a
*          negative block number.
*
*  ISPLIT  (global output) INTEGER array, dimension (N)
*          The splitting points, at which T breaks up into submatrices.
*          The first submatrix consists of rows/columns 1 to ISPLIT(1),
*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
*          etc., and the NSPLIT-th consists of rows/columns
*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
*          (Only the first NSPLIT elements will actually be used, but
*          since the user cannot know a priori what value NSPLIT will
*          have, N words must be reserved for ISPLIT.)
*
*  WORK    (local workspace) DOUBLE PRECISION array,
*          dimension ( MAX( 5*N, 7 ) )
*
*  LWORK   (local input) INTEGER
*          size of array WORK must be >= MAX( 5*N, 7 )
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  IWORK   (local workspace) INTEGER array, dimension ( MAX( 4*N, 14 ) )
*
*  LIWORK  (local input) INTEGER
*          size of array IWORK must be >= MAX( 4*N, 14, NPROCS )
*          If LIWORK = -1, then LIWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0 :  successful exit
*          < 0 :  if INFO = -i, the i-th argument had an illegal value
*          > 0 :  some or all of the eigenvalues failed to converge or
*                 were not computed:
*              = 1 : Bisection failed to converge for some eigenvalues;
*                    these eigenvalues are flagged by a negative block
*                    number.  The effect is that the eigenvalues may not
*                    be as accurate as the absolute and relative
*                    tolerances. This is generally caused by arithmetic
*                    which is less accurate than PDLAMCH says.
*              = 2 : There is a mismatch between the number of
*                    eigenvalues output and the number desired.
*              = 3 : RANGE='i', and the Gershgorin interval initially
*                    used was incorrect. No eigenvalues were computed.
*                    Probable cause: your machine has sloppy floating
*                    point arithmetic.
*                    Cure: Increase the PARAMETER "FUDGE", recompile,
*                    and try again.
*
*  Internal Parameters
*  ===================
*
*  RELFAC  DOUBLE PRECISION, default = 2.0
*          The relative tolerance.  An interval [a,b] lies within
*          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),
*          where "ulp" is the machine precision (distance from 1 to
*          the next larger floating point number.)
*
*  FUDGE   DOUBLE PRECISION, default = 2.0
*          A "fudge factor" to widen the Gershgorin intervals.  Ideally,
*          a value of 1 should work, but on machines with sloppy
*          arithmetic, this needs to be larger.  The default for
*          publicly released versions should be large enough to handle
*          the worst machine around.  Note that this has no effect
*          on the accuracy of the solution.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, ICHAR, MAX, MIN, MOD
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            BLACS_PNUM
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           LSAME, BLACS_PNUM, PDLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_FREEBUFF, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDMAP, DGEBR2D,
     $                   DGEBS2D, DGERV2D, DGESD2D, DLASRT2, GLOBCHK,
     $                   IGEBR2D, IGEBS2D, IGERV2D, IGESD2D, IGSUM2D,
     $                   PDLAEBZ, PDLAIECTB, PDLAIECTL, PDLAPDCT,
     $                   PDLASNBT, PXERBLA
*     ..
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            BIGNUM, DESCMULT
      PARAMETER          ( BIGNUM = 10000, DESCMULT = 100 )
      DOUBLE PRECISION   ZERO, ONE, TWO, FIVE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0,
     $                   FIVE = 5.0D+0, HALF = 1.0D+0 / TWO )
      DOUBLE PRECISION   FUDGE, RELFAC
      PARAMETER          ( FUDGE = 2.0D+0, RELFAC = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            BLKNO, FOUND, I, IBEGIN, IEFLAG, IEND, IFRST,
     $                   IINFO, ILAST, ILOAD, IM, IMYLOAD, IN, INDRIW1,
     $                   INDRIW2, INDRW1, INDRW2, INXTLOAD, IOFF,
     $                   IORDER, IOUT, IRANGE, IRECV, IREM, ITMP1,
     $                   ITMP2, J, JB, K, LAST, LEXTRA, LREQ, MYCOL,
     $                   MYROW, NALPHA, NBETA, NCMP, NEIGINT, NEXT, NGL,
     $                   NGLOB, NGU, NINT, NPCOL, NPROW, OFFSET,
     $                   ONEDCONTEXT, P, PREV, REXTRA, RREQ, SELF,
     $                   TORECV
      DOUBLE PRECISION   ALPHA, ATOLI, BETA, BNORM, DRECV, DSEND, GL,
     $                   GU, INITVL, INITVU, LSAVE, MID, PIVMIN, RELTOL,
     $                   SAFEMN, TMP1, TMP2, TNORM, ULP
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM( 5, 2 )
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Set up process grid
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      INFO = 0
      M = 0
*
*     Decode RANGE
*
      IF( LSAME( RANGE, 'A' ) ) THEN
         IRANGE = 1
      ELSE IF( LSAME( RANGE, 'V' ) ) THEN
         IRANGE = 2
      ELSE IF( LSAME( RANGE, 'I' ) ) THEN
         IRANGE = 3
      ELSE
         IRANGE = 0
      END IF
*
*     Decode ORDER
*
      IF( LSAME( ORDER, 'B' ) ) THEN
         IORDER = 2
      ELSE IF( LSAME( ORDER, 'E' ) .OR. LSAME( ORDER, 'A' ) ) THEN
         IORDER = 1
      ELSE
         IORDER = 0
      END IF
*
*     Check for Errors
*
      IF( NPROW.EQ.-1 ) THEN
         INFO = -1
      ELSE
*
*     Get machine constants
*
         SAFEMN = PDLAMCH( ICTXT, 'S' )
         ULP = PDLAMCH( ICTXT, 'P' )
         RELTOL = ULP*RELFAC
         IDUM( 1, 1 ) = ICHAR( RANGE )
         IDUM( 1, 2 ) = 2
         IDUM( 2, 1 ) = ICHAR( ORDER )
         IDUM( 2, 2 ) = 3
         IDUM( 3, 1 ) = N
         IDUM( 3, 2 ) = 4
         NGLOB = 5
         IF( IRANGE.EQ.3 ) THEN
            IDUM( 4, 1 ) = IL
            IDUM( 4, 2 ) = 7
            IDUM( 5, 1 ) = IU
            IDUM( 5, 2 ) = 8
         ELSE
            IDUM( 4, 1 ) = 0
            IDUM( 4, 2 ) = 0
            IDUM( 5, 1 ) = 0
            IDUM( 5, 2 ) = 0
         END IF
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            WORK( 1 ) = ABSTOL
            IF( IRANGE.EQ.2 ) THEN
               WORK( 2 ) = VL
               WORK( 3 ) = VU
            ELSE
               WORK( 2 ) = ZERO
               WORK( 3 ) = ZERO
            END IF
            CALL DGEBS2D( ICTXT, 'ALL', ' ', 3, 1, WORK, 3 )
         ELSE
            CALL DGEBR2D( ICTXT, 'ALL', ' ', 3, 1, WORK, 3, 0, 0 )
         END IF
         LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
         IF( INFO.EQ.0 ) THEN
            IF( IRANGE.EQ.0 ) THEN
               INFO = -2
            ELSE IF( IORDER.EQ.0 ) THEN
               INFO = -3
            ELSE IF( IRANGE.EQ.2 .AND. VL.GE.VU ) THEN
               INFO = -5
            ELSE IF( IRANGE.EQ.3 .AND. ( IL.LT.1 .OR. IL.GT.MAX( 1,
     $              N ) ) ) THEN
               INFO = -6
            ELSE IF( IRANGE.EQ.3 .AND. ( IU.LT.MIN( N,
     $              IL ) .OR. IU.GT.N ) ) THEN
               INFO = -7
            ELSE IF( LWORK.LT.MAX( 5*N, 7 ) .AND. .NOT.LQUERY ) THEN
               INFO = -18
            ELSE IF( LIWORK.LT.MAX( 4*N, 14, NPROW*NPCOL ) .AND. .NOT.
     $              LQUERY ) THEN
               INFO = -20
            ELSE IF( IRANGE.EQ.2 .AND. ( ABS( WORK( 2 )-VL ).GT.FIVE*
     $              ULP*ABS( VL ) ) ) THEN
               INFO = -5
            ELSE IF( IRANGE.EQ.2 .AND. ( ABS( WORK( 3 )-VU ).GT.FIVE*
     $              ULP*ABS( VU ) ) ) THEN
               INFO = -6
            ELSE IF( ABS( WORK( 1 )-ABSTOL ).GT.FIVE*ULP*ABS( ABSTOL ) )
     $               THEN
               INFO = -9
            END IF
         END IF
         IF( INFO.EQ.0 )
     $      INFO = BIGNUM
         CALL GLOBCHK( ICTXT, NGLOB, IDUM, 5, IWORK, INFO )
         IF( INFO.EQ.BIGNUM ) THEN
            INFO = 0
         ELSE IF( MOD( INFO, DESCMULT ).EQ.0 ) THEN
            INFO = -INFO / DESCMULT
         ELSE
            INFO = -INFO
         END IF
      END IF
      WORK( 1 ) = DBLE( MAX( 5*N, 7 ) )
      IWORK( 1 ) = MAX( 4*N, 14, NPROW*NPCOL )
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDSTEBZ', -INFO )
         RETURN
      ELSE IF( LWORK.EQ.-1 .AND. LIWORK.EQ.-1 ) THEN
         RETURN
      END IF
*
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      K = 1
      DO 20 I = 0, NPROW - 1
         DO 10 J = 0, NPCOL - 1
            IWORK( K ) = BLACS_PNUM( ICTXT, I, J )
            K = K + 1
   10    CONTINUE
   20 CONTINUE
*
      P = NPROW*NPCOL
      NPROW = 1
      NPCOL = P
*
      CALL BLACS_GET( ICTXT, 10, ONEDCONTEXT )
      CALL BLACS_GRIDMAP( ONEDCONTEXT, IWORK, NPROW, NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ONEDCONTEXT, I, J, K, SELF )
*
*     Simplifications:
*
      IF( IRANGE.EQ.3 .AND. IL.EQ.1 .AND. IU.EQ.N )
     $   IRANGE = 1
*
      NEXT = MOD( SELF+1, P )
      PREV = MOD( P+SELF-1, P )
*
*     Compute squares of off-diagonals, splitting points and pivmin.
*     Interleave diagonals and off-diagonals.
*
      INDRW1 = MAX( 2*N, 4 )
      INDRW2 = INDRW1 + 2*N
      INDRIW1 = MAX( 2*N, 8 )
      NSPLIT = 1
      WORK( INDRW1+2*N ) = ZERO
      PIVMIN = ONE
*
      DO 30 I = 1, N - 1
         TMP1 = E( I )**2
         J = 2*I
         WORK( INDRW1+J-1 ) = D( I )
         IF( ABS( D( I+1 )*D( I ) )*ULP**2+SAFEMN.GT.TMP1 ) THEN
            ISPLIT( NSPLIT ) = I
            NSPLIT = NSPLIT + 1
            WORK( INDRW1+J ) = ZERO
         ELSE
            WORK( INDRW1+J ) = TMP1
            PIVMIN = MAX( PIVMIN, TMP1 )
         END IF
   30 CONTINUE
      WORK( INDRW1+2*N-1 ) = D( N )
      ISPLIT( NSPLIT ) = N
      PIVMIN = PIVMIN*SAFEMN
*
*     Compute Gershgorin interval [gl,gu] for entire matrix
*
      GU = D( 1 )
      GL = D( 1 )
      TMP1 = ZERO
*
      DO 40 I = 1, N - 1
         TMP2 = ABS( E( I ) )
         GU = MAX( GU, D( I )+TMP1+TMP2 )
         GL = MIN( GL, D( I )-TMP1-TMP2 )
         TMP1 = TMP2
   40 CONTINUE
      GU = MAX( GU, D( N )+TMP1 )
      GL = MIN( GL, D( N )-TMP1 )
      TNORM = MAX( ABS( GL ), ABS( GU ) )
      GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
      GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN
*
      IF( ABSTOL.LE.ZERO ) THEN
         ATOLI = ULP*TNORM
      ELSE
         ATOLI = ABSTOL
      END IF
*
*     Find out if on an IEEE machine, the sign bit is the
*     32nd bit (Big Endian) or the 64th bit (Little Endian)
*
      IF( IRANGE.EQ.1 .OR. NSPLIT.EQ.1 ) THEN
         CALL PDLASNBT( IEFLAG )
      ELSE
         IEFLAG = 0
      END IF
      LEXTRA = 0
      REXTRA = 0
*
*     Form Initial Interval containing desired eigenvalues
*
      IF( IRANGE.EQ.1 ) THEN
         INITVL = GL
         INITVU = GU
         WORK( 1 ) = GL
         WORK( 2 ) = GU
         IWORK( 1 ) = 0
         IWORK( 2 ) = N
         IFRST = 1
         ILAST = N
      ELSE IF( IRANGE.EQ.2 ) THEN
         IF( VL.GT.GL ) THEN
            IF( IEFLAG.EQ.0 ) THEN
               CALL PDLAPDCT( VL, N, WORK( INDRW1+1 ), PIVMIN, IFRST )
            ELSE IF( IEFLAG.EQ.1 ) THEN
               CALL PDLAIECTB( VL, N, WORK( INDRW1+1 ), IFRST )
            ELSE
               CALL PDLAIECTL( VL, N, WORK( INDRW1+1 ), IFRST )
            END IF
            IFRST = IFRST + 1
            INITVL = VL
         ELSE
            INITVL = GL
            IFRST = 1
         END IF
         IF( VU.LT.GU ) THEN
            IF( IEFLAG.EQ.0 ) THEN
               CALL PDLAPDCT( VU, N, WORK( INDRW1+1 ), PIVMIN, ILAST )
            ELSE IF( IEFLAG.EQ.1 ) THEN
               CALL PDLAIECTB( VU, N, WORK( INDRW1+1 ), ILAST )
            ELSE
               CALL PDLAIECTL( VU, N, WORK( INDRW1+1 ), ILAST )
            END IF
            INITVU = VU
         ELSE
            INITVU = GU
            ILAST = N
         END IF
         WORK( 1 ) = INITVL
         WORK( 2 ) = INITVU
         IWORK( 1 ) = IFRST - 1
         IWORK( 2 ) = ILAST
      ELSE IF( IRANGE.EQ.3 ) THEN
         WORK( 1 ) = GL
         WORK( 2 ) = GU
         IWORK( 1 ) = 0
         IWORK( 2 ) = N
         IWORK( 5 ) = IL - 1
         IWORK( 6 ) = IU
         CALL PDLAEBZ( 0, N, 2, 1, ATOLI, RELTOL, PIVMIN,
     $                 WORK( INDRW1+1 ), IWORK( 5 ), WORK, IWORK, NINT,
     $                 LSAVE, IEFLAG, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 3
            GO TO 230
         END IF
         IF( NINT.GT.1 ) THEN
            IF( IWORK( 5 ).EQ.IL-1 ) THEN
               WORK( 2 ) = WORK( 4 )
               IWORK( 2 ) = IWORK( 4 )
            ELSE
               WORK( 1 ) = WORK( 3 )
               IWORK( 1 ) = IWORK( 3 )
            END IF
            IF( IWORK( 1 ).LT.0 .OR. IWORK( 1 ).GT.IL-1 .OR.
     $          IWORK( 2 ).LE.MIN( IU-1, IWORK( 1 ) ) .OR.
     $          IWORK( 2 ).GT.N ) THEN
               INFO = 3
               GO TO 230
            END IF
         END IF
         LEXTRA = IL - 1 - IWORK( 1 )
         REXTRA = IWORK( 2 ) - IU
         INITVL = WORK( 1 )
         INITVU = WORK( 2 )
         IFRST = IL
         ILAST = IU
      END IF
*     NVL = IFRST - 1
*     NVU = ILAST
      GL = INITVL
      GU = INITVU
      NGL = IWORK( 1 )
      NGU = IWORK( 2 )
      IM = 0
      FOUND = 0
      INDRIW2 = INDRIW1 + NGU - NGL
      IEND = 0
      IF( IFRST.GT.ILAST )
     $   GO TO 100
      IF( IFRST.EQ.1 .AND. ILAST.EQ.N )
     $   IRANGE = 1
*
*     Find Eigenvalues -- Loop Over Blocks
*
      DO 90 JB = 1, NSPLIT
         IOFF = IEND
         IBEGIN = IOFF + 1
         IEND = ISPLIT( JB )
         IN = IEND - IOFF
         IF( JB.NE.1 ) THEN
            IF( IRANGE.NE.1 ) THEN
               FOUND = IM
*
*              Find total number of eigenvalues found thus far
*
               CALL IGSUM2D( ONEDCONTEXT, 'All', ' ', 1, 1, FOUND, 1,
     $                       -1, -1 )
            ELSE
               FOUND = IOFF
            END IF
         END IF
*         IF( SELF.GE.P )
*     $      GO TO 30
         IF( IN.NE.N ) THEN
*
*           Compute Gershgorin interval [gl,gu] for split matrix
*
            GU = D( IBEGIN )
            GL = D( IBEGIN )
            TMP1 = ZERO
*
            DO 50 J = IBEGIN, IEND - 1
               TMP2 = ABS( E( J ) )
               GU = MAX( GU, D( J )+TMP1+TMP2 )
               GL = MIN( GL, D( J )-TMP1-TMP2 )
               TMP1 = TMP2
   50       CONTINUE
*
            GU = MAX( GU, D( IEND )+TMP1 )
            GL = MIN( GL, D( IEND )-TMP1 )
            BNORM = MAX( ABS( GL ), ABS( GU ) )
            GL = GL - FUDGE*BNORM*ULP*IN - FUDGE*PIVMIN
            GU = GU + FUDGE*BNORM*ULP*IN + FUDGE*PIVMIN
*
*           Compute ATOLI for the current submatrix
*
            IF( ABSTOL.LE.ZERO ) THEN
               ATOLI = ULP*BNORM
            ELSE
               ATOLI = ABSTOL
            END IF
*
            IF( GL.LT.INITVL ) THEN
               GL = INITVL
               IF( IEFLAG.EQ.0 ) THEN
                  CALL PDLAPDCT( GL, IN, WORK( INDRW1+2*IOFF+1 ),
     $                           PIVMIN, NGL )
               ELSE IF( IEFLAG.EQ.1 ) THEN
                  CALL PDLAIECTB( GL, IN, WORK( INDRW1+2*IOFF+1 ), NGL )
               ELSE
                  CALL PDLAIECTL( GL, IN, WORK( INDRW1+2*IOFF+1 ), NGL )
               END IF
            ELSE
               NGL = 0
            END IF
            IF( GU.GT.INITVU ) THEN
               GU = INITVU
               IF( IEFLAG.EQ.0 ) THEN
                  CALL PDLAPDCT( GU, IN, WORK( INDRW1+2*IOFF+1 ),
     $                           PIVMIN, NGU )
               ELSE IF( IEFLAG.EQ.1 ) THEN
                  CALL PDLAIECTB( GU, IN, WORK( INDRW1+2*IOFF+1 ), NGU )
               ELSE
                  CALL PDLAIECTL( GU, IN, WORK( INDRW1+2*IOFF+1 ), NGU )
               END IF
            ELSE
               NGU = IN
            END IF
            IF( NGL.GE.NGU )
     $         GO TO 90
            WORK( 1 ) = GL
            WORK( 2 ) = GU
            IWORK( 1 ) = NGL
            IWORK( 2 ) = NGU
         END IF
         OFFSET = FOUND - NGL
         BLKNO = JB
*
*        Do a static partitioning of work so that each process
*        has to find an (almost) equal number of eigenvalues
*
         NCMP = NGU - NGL
         ILOAD = NCMP / P
         IREM = NCMP - ILOAD*P
         ITMP1 = MOD( SELF-FOUND, P )
         IF( ITMP1.LT.0 )
     $      ITMP1 = ITMP1 + P
         IF( ITMP1.LT.IREM ) THEN
            IMYLOAD = ILOAD + 1
         ELSE
            IMYLOAD = ILOAD
         END IF
         IF( IMYLOAD.EQ.0 ) THEN
            GO TO 90
         ELSE IF( IN.EQ.1 ) THEN
            WORK( INDRW2+IM+1 ) = WORK( INDRW1+2*IOFF+1 )
            IWORK( INDRIW1+IM+1 ) = BLKNO
            IWORK( INDRIW2+IM+1 ) = OFFSET + 1
            IM = IM + 1
            GO TO 90
         ELSE
            INXTLOAD = ILOAD
            ITMP2 = MOD( SELF+1-FOUND, P )
            IF( ITMP2.LT.0 )
     $         ITMP2 = ITMP2 + P
            IF( ITMP2.LT.IREM )
     $         INXTLOAD = INXTLOAD + 1
            LREQ = NGL + ITMP1*ILOAD + MIN( IREM, ITMP1 )
            RREQ = LREQ + IMYLOAD
            IWORK( 5 ) = LREQ
            IWORK( 6 ) = RREQ
            TMP1 = WORK( 1 )
            ITMP1 = IWORK( 1 )
            CALL PDLAEBZ( 1, IN, 1, 1, ATOLI, RELTOL, PIVMIN,
     $                    WORK( INDRW1+2*IOFF+1 ), IWORK( 5 ), WORK,
     $                    IWORK, NINT, LSAVE, IEFLAG, IINFO )
            ALPHA = WORK( 1 )
            BETA = WORK( 2 )
            NALPHA = IWORK( 1 )
            NBETA = IWORK( 2 )
            DSEND = BETA
            IF( NBETA.GT.RREQ+INXTLOAD ) THEN
               NBETA = RREQ
               DSEND = ALPHA
            END IF
            LAST = MOD( FOUND+MIN( NGU-NGL, P )-1, P )
            IF( LAST.LT.0 )
     $         LAST = LAST + P
            IF( SELF.NE.LAST ) THEN
               CALL DGESD2D( ONEDCONTEXT, 1, 1, DSEND, 1, 0, NEXT )
               CALL IGESD2D( ONEDCONTEXT, 1, 1, NBETA, 1, 0, NEXT )
            END IF
            IF( SELF.NE.MOD( FOUND, P ) ) THEN
               CALL DGERV2D( ONEDCONTEXT, 1, 1, DRECV, 1, 0, PREV )
               CALL IGERV2D( ONEDCONTEXT, 1, 1, IRECV, 1, 0, PREV )
            ELSE
               DRECV = TMP1
               IRECV = ITMP1
            END IF
            WORK( 1 ) = MAX( LSAVE, DRECV )
            IWORK( 1 ) = IRECV
            ALPHA = MAX( ALPHA, WORK( 1 ) )
            NALPHA = MAX( NALPHA, IRECV )
            IF( BETA-ALPHA.LE.MAX( ATOLI, RELTOL*MAX( ABS( ALPHA ),
     $          ABS( BETA ) ) ) ) THEN
               MID = HALF*( ALPHA+BETA )
               DO 60 J = OFFSET + NALPHA + 1, OFFSET + NBETA
                  WORK( INDRW2+IM+1 ) = MID
                  IWORK( INDRIW1+IM+1 ) = BLKNO
                  IWORK( INDRIW2+IM+1 ) = J
                  IM = IM + 1
   60          CONTINUE
               WORK( 2 ) = ALPHA
               IWORK( 2 ) = NALPHA
            END IF
         END IF
         NEIGINT = IWORK( 2 ) - IWORK( 1 )
         IF( NEIGINT.LE.0 )
     $      GO TO 90
*
*        Call the main computational routine
*
         CALL PDLAEBZ( 2, IN, NEIGINT, 1, ATOLI, RELTOL, PIVMIN,
     $                 WORK( INDRW1+2*IOFF+1 ), IWORK, WORK, IWORK,
     $                 IOUT, LSAVE, IEFLAG, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 1
         END IF
         DO 80 I = 1, IOUT
            MID = HALF*( WORK( 2*I-1 )+WORK( 2*I ) )
            IF( I.GT.IOUT-IINFO )
     $         BLKNO = -BLKNO
            DO 70 J = OFFSET + IWORK( 2*I-1 ) + 1,
     $              OFFSET + IWORK( 2*I )
               WORK( INDRW2+IM+1 ) = MID
               IWORK( INDRIW1+IM+1 ) = BLKNO
               IWORK( INDRIW2+IM+1 ) = J
               IM = IM + 1
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
*     Find out total number of eigenvalues computed
*
  100 CONTINUE
      M = IM
      CALL IGSUM2D( ONEDCONTEXT, 'ALL', ' ', 1, 1, M, 1, -1, -1 )
*
*     Move the eigenvalues found to their final destinations
*
      DO 130 I = 1, P
         IF( SELF.EQ.I-1 ) THEN
            CALL IGEBS2D( ONEDCONTEXT, 'ALL', ' ', 1, 1, IM, 1 )
            IF( IM.NE.0 ) THEN
               CALL IGEBS2D( ONEDCONTEXT, 'ALL', ' ', IM, 1,
     $                       IWORK( INDRIW2+1 ), IM )
               CALL DGEBS2D( ONEDCONTEXT, 'ALL', ' ', IM, 1,
     $                       WORK( INDRW2+1 ), IM )
               CALL IGEBS2D( ONEDCONTEXT, 'ALL', ' ', IM, 1,
     $                       IWORK( INDRIW1+1 ), IM )
               DO 110 J = 1, IM
                  W( IWORK( INDRIW2+J ) ) = WORK( INDRW2+J )
                  IBLOCK( IWORK( INDRIW2+J ) ) = IWORK( INDRIW1+J )
  110          CONTINUE
            END IF
         ELSE
            CALL IGEBR2D( ONEDCONTEXT, 'ALL', ' ', 1, 1, TORECV, 1, 0,
     $                    I-1 )
            IF( TORECV.NE.0 ) THEN
               CALL IGEBR2D( ONEDCONTEXT, 'ALL', ' ', TORECV, 1, IWORK,
     $                       TORECV, 0, I-1 )
               CALL DGEBR2D( ONEDCONTEXT, 'ALL', ' ', TORECV, 1, WORK,
     $                       TORECV, 0, I-1 )
               CALL IGEBR2D( ONEDCONTEXT, 'ALL', ' ', TORECV, 1,
     $                       IWORK( N+1 ), TORECV, 0, I-1 )
               DO 120 J = 1, TORECV
                  W( IWORK( J ) ) = WORK( J )
                  IBLOCK( IWORK( J ) ) = IWORK( N+J )
  120          CONTINUE
            END IF
         END IF
  130 CONTINUE
      IF( NSPLIT.GT.1 .AND. IORDER.EQ.1 ) THEN
*
*        Sort the eigenvalues
*
*
         DO 140 I = 1, M
            IWORK( M+I ) = I
  140    CONTINUE
         CALL DLASRT2( 'I', M, W, IWORK( M+1 ), IINFO )
         DO 150 I = 1, M
            IWORK( I ) = IBLOCK( I )
  150    CONTINUE
         DO 160 I = 1, M
            IBLOCK( I ) = IWORK( IWORK( M+I ) )
  160    CONTINUE
      END IF
      IF( IRANGE.EQ.3 .AND. ( LEXTRA.GT.0 .OR. REXTRA.GT.0 ) ) THEN
*
*        Discard unwanted eigenvalues (occurs only when RANGE = 'I',
*        and eigenvalues IL, and/or IU are in a cluster)
*
         DO 170 I = 1, M
            WORK( I ) = W( I )
            IWORK( I ) = I
            IWORK( M+I ) = I
  170    CONTINUE
         DO 190 I = 1, LEXTRA
            ITMP1 = I
            DO 180 J = I + 1, M
               IF( WORK( J ).LT.WORK( ITMP1 ) ) THEN
                  ITMP1 = J
               END IF
  180       CONTINUE
            TMP1 = WORK( I )
            WORK( I ) = WORK( ITMP1 )
            WORK( ITMP1 ) = TMP1
            IWORK( IWORK( M+ITMP1 ) ) = I
            IWORK( IWORK( M+I ) ) = ITMP1
            ITMP2 = IWORK( M+I )
            IWORK( M+I ) = IWORK( M+ITMP1 )
            IWORK( M+ITMP1 ) = ITMP2
  190    CONTINUE
         DO 210 I = 1, REXTRA
            ITMP1 = M - I + 1
            DO 200 J = M - I, LEXTRA + 1, -1
               IF( WORK( J ).GT.WORK( ITMP1 ) ) THEN
                  ITMP1 = J
               END IF
  200       CONTINUE
            TMP1 = WORK( M-I+1 )
            WORK( M-I+1 ) = WORK( ITMP1 )
            WORK( ITMP1 ) = TMP1
            IWORK( IWORK( M+ITMP1 ) ) = M - I + 1
            IWORK( IWORK( 2*M-I+1 ) ) = ITMP1
            ITMP2 = IWORK( 2*M-I+1 )
            IWORK( 2*M-I+1 ) = IWORK( M+ITMP1 )
            IWORK( M+ITMP1 ) = ITMP2
*           IWORK( ITMP1 ) = 1
  210    CONTINUE
         J = 0
         DO 220 I = 1, M
            IF( IWORK( I ).GT.LEXTRA .AND. IWORK( I ).LE.M-REXTRA ) THEN
               J = J + 1
               W( J ) = WORK( IWORK( I ) )
               IBLOCK( J ) = IBLOCK( I )
            END IF
  220    CONTINUE
         M = M - LEXTRA - REXTRA
      END IF
      IF( M.NE.ILAST-IFRST+1 ) THEN
         INFO = 2
      END IF
*
  230 CONTINUE
      CALL BLACS_FREEBUFF( ONEDCONTEXT, 1 )
      CALL BLACS_GRIDEXIT( ONEDCONTEXT )
      RETURN
*
*     End of PDSTEBZ
*
      END
*
      SUBROUTINE PDLAEBZ( IJOB, N, MMAX, MINP, ABSTOL, RELTOL, PIVMIN,
     $                    D, NVAL, INTVL, INTVLCT, MOUT, LSAVE, IEFLAG,
     $                    INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*
*     .. Scalar Arguments ..
      INTEGER            IEFLAG, IJOB, INFO, MINP, MMAX, MOUT, N
      DOUBLE PRECISION   ABSTOL, LSAVE, PIVMIN, RELTOL
*     ..
*     .. Array Arguments ..
      INTEGER            INTVLCT( * ), NVAL( * )
      DOUBLE PRECISION   D( * ), INTVL( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAEBZ contains the iteration loop which computes the eigenvalues
*  contained in the input intervals [ INTVL(2*j-1), INTVL(2*j) ] where
*  j = 1,...,MINP. It uses and computes the function N(w), which is
*  the count of eigenvalues of a symmetric tridiagonal matrix less than
*  or equal to its argument w.
*
*  This is a ScaLAPACK internal subroutine and arguments are not
*  checked for unreasonable values.
*
*  Arguments
*  =========
*
*  IJOB    (input) INTEGER
*          Specifies the computation done by PDLAEBZ
*          = 0 : Find an interval with desired values of N(w) at the
*                endpoints of the interval.
*          = 1 : Find a floating point number contained in the initial
*                interval with a desired value of N(w).
*          = 2 : Perform bisection iteration to find eigenvalues of T.
*
*  N       (input) INTEGER
*          The order of the tridiagonal matrix T. N >= 1.
*
*  MMAX    (input) INTEGER
*          The maximum number of intervals that may be generated. If
*          more than MMAX intervals are generated, then PDLAEBZ will
*          quit with INFO = MMAX+1.
*
*  MINP    (input) INTEGER
*          The initial number of intervals. MINP <= MMAX.
*
*  ABSTOL  (input) DOUBLE PRECISION
*          The minimum (absolute) width of an interval. When an interval
*          is narrower than ABSTOL, or than RELTOL times the larger (in
*          magnitude) endpoint, then it is considered to be sufficiently
*          small, i.e., converged.
*          This must be at least zero.
*
*  RELTOL  (input) DOUBLE PRECISION
*          The minimum relative width of an interval. When an interval
*          is narrower than ABSTOL, or than RELTOL times the larger (in
*          magnitude) endpoint, then it is considered to be sufficiently
*          small, i.e., converged.
*          Note : This should be at least radix*machine epsilon.
*
*  PIVMIN  (input) DOUBLE PRECISION
*          The minimum absolute of a "pivot" in the "paranoid"
*          implementation of the Sturm sequence loop. This must be at
*          least max_j |e(j)^2| *safe_min, and at least safe_min, where
*          safe_min is at least the smallest number that can divide 1.0
*          without overflow.
*          See PDLAPDCT for the "paranoid" implementation of the Sturm
*          sequence loop.
*
*  D       (input) DOUBLE PRECISION array, dimension (2*N - 1)
*          Contains the diagonals and the squares of the off-diagonal
*          elements of the tridiagonal matrix T. These elements are
*          assumed to be interleaved in memory for better cache
*          performance. The diagonal entries of T are in the entries
*          D(1),D(3),...,D(2*N-1), while the squares of the off-diagonal
*          entries are D(2),D(4),...,D(2*N-2). To avoid overflow, the
*          matrix must be scaled so that its largest entry is no greater
*          than overflow**(1/2) * underflow**(1/4) in absolute value,
*          and for greatest accuracy, it should not be much smaller
*          than that.
*
*  NVAL    (input/output) INTEGER array, dimension (4)
*          If IJOB = 0, the desired values of N(w) are in NVAL(1) and
*                       NVAL(2).
*          If IJOB = 1, NVAL(2) is the desired value of N(w).
*          If IJOB = 2, not referenced.
*          This array will, in general, be reordered on output.
*
*  INTVL   (input/output) DOUBLE PRECISION array, dimension (2*MMAX)
*          The endpoints of the intervals. INTVL(2*j-1) is the left
*          endpoint of the j-th interval, and INTVL(2*j) is the right
*          endpoint of the j-th interval. The input intervals will,
*          in general, be modified, split and reordered by the
*          calculation.
*          On input, INTVL contains the MINP input intervals.
*          On output, INTVL contains the converged intervals.
*
*  INTVLCT (input/output) INTEGER array, dimension (2*MMAX)
*          The counts at the endpoints of the intervals. INTVLCT(2*j-1)
*          is the count at the left endpoint of the j-th interval, i.e.,
*          the function value N(INTVL(2*j-1)), and INTVLCT(2*j) is the
*          count at the right endpoint of the j-th interval.
*          On input, INTVLCT contains the counts at the endpoints of
*          the MINP input intervals.
*          On output, INTVLCT contains the counts at the endpoints of
*          the converged intervals.
*
*  MOUT    (output) INTEGER
*          The number of intervals output.
*
*  LSAVE   (output) DOUBLE PRECISION
*          If IJOB = 0 or 2, not referenced.
*          If IJOB = 1, this is the largest floating point number
*          encountered which has count N(w) = NVAL(1).
*
*  IEFLAG  (input) INTEGER
*          A flag which indicates whether N(w) should be speeded up by
*          exploiting IEEE Arithmetic.
*
*  INFO    (output) INTEGER
*          = 0        : All intervals converged.
*          = 1 - MMAX : The last INFO intervals did not converge.
*          = MMAX + 1 : More than MMAX intervals were generated.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           PDLAECV, PDLAIECTB, PDLAIECTL, PDLAPDCT
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TWO, HALF
      PARAMETER          ( ZERO = 0.0D+0, TWO = 2.0D+0,
     $                   HALF = 1.0D+0 / TWO )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ITMAX, J, K, KF, KL, KLNEW, L, LCNT, LREQ,
     $                   NALPHA, NBETA, NMID, RCNT, RREQ
      DOUBLE PRECISION   ALPHA, BETA, MID
*     ..
*     .. Executable Statements ..
*
      KF = 1
      KL = MINP + 1
      INFO = 0
      IF( INTVL( 2 )-INTVL( 1 ).LE.ZERO ) THEN
         INFO = MINP
         MOUT = KF
         RETURN
      END IF
      IF( IJOB.EQ.0 ) THEN
*
*        Check if some input intervals have "converged"
*
         CALL PDLAECV( 0, KF, KL, INTVL, INTVLCT, NVAL,
     $                 MAX( ABSTOL, PIVMIN ), RELTOL )
         IF( KF.GE.KL )
     $      GO TO 60
*
*        Compute upper bound on number of iterations needed
*
         ITMAX = INT( ( LOG( INTVL( 2 )-INTVL( 1 )+PIVMIN )-
     $           LOG( PIVMIN ) ) / LOG( TWO ) ) + 2
*
*        Iteration Loop
*
         DO 20 I = 1, ITMAX
            KLNEW = KL
            DO 10 J = KF, KL - 1
               K = 2*J
*
*           Bisect the interval and find the count at that point
*
               MID = HALF*( INTVL( K-1 )+INTVL( K ) )
               IF( IEFLAG.EQ.0 ) THEN
                  CALL PDLAPDCT( MID, N, D, PIVMIN, NMID )
               ELSE IF( IEFLAG.EQ.1 ) THEN
                  CALL PDLAIECTB( MID, N, D, NMID )
               ELSE
                  CALL PDLAIECTL( MID, N, D, NMID )
               END IF
               LREQ = NVAL( K-1 )
               RREQ = NVAL( K )
               IF( KL.EQ.1 )
     $            NMID = MIN( INTVLCT( K ),
     $                   MAX( INTVLCT( K-1 ), NMID ) )
               IF( NMID.LE.NVAL( K-1 ) ) THEN
                  INTVL( K-1 ) = MID
                  INTVLCT( K-1 ) = NMID
               END IF
               IF( NMID.GE.NVAL( K ) ) THEN
                  INTVL( K ) = MID
                  INTVLCT( K ) = NMID
               END IF
               IF( NMID.GT.LREQ .AND. NMID.LT.RREQ ) THEN
                  L = 2*KLNEW
                  INTVL( L-1 ) = MID
                  INTVL( L ) = INTVL( K )
                  INTVLCT( L-1 ) = NVAL( K )
                  INTVLCT( L ) = INTVLCT( K )
                  INTVL( K ) = MID
                  INTVLCT( K ) = NVAL( K-1 )
                  NVAL( L-1 ) = NVAL( K )
                  NVAL( L ) = NVAL( L-1 )
                  NVAL( K ) = NVAL( K-1 )
                  KLNEW = KLNEW + 1
               END IF
   10       CONTINUE
            KL = KLNEW
            CALL PDLAECV( 0, KF, KL, INTVL, INTVLCT, NVAL,
     $                    MAX( ABSTOL, PIVMIN ), RELTOL )
            IF( KF.GE.KL )
     $         GO TO 60
   20    CONTINUE
      ELSE IF( IJOB.EQ.1 ) THEN
         ALPHA = INTVL( 1 )
         BETA = INTVL( 2 )
         NALPHA = INTVLCT( 1 )
         NBETA = INTVLCT( 2 )
         LSAVE = ALPHA
         LREQ = NVAL( 1 )
         RREQ = NVAL( 2 )
   30    CONTINUE
         IF( NBETA.NE.RREQ .AND. BETA-ALPHA.GT.
     $       MAX( ABSTOL, RELTOL*MAX( ABS( ALPHA ), ABS( BETA ) ) ) )
     $        THEN
*
*           Bisect the interval and find the count at that point
*
            MID = HALF*( ALPHA+BETA )
            IF( IEFLAG.EQ.0 ) THEN
               CALL PDLAPDCT( MID, N, D, PIVMIN, NMID )
            ELSE IF( IEFLAG.EQ.1 ) THEN
               CALL PDLAIECTB( MID, N, D, NMID )
            ELSE
               CALL PDLAIECTL( MID, N, D, NMID )
            END IF
            NMID = MIN( NBETA, MAX( NALPHA, NMID ) )
            IF( NMID.GE.RREQ ) THEN
               BETA = MID
               NBETA = NMID
            ELSE
               ALPHA = MID
               NALPHA = NMID
               IF( NMID.EQ.LREQ )
     $            LSAVE = ALPHA
            END IF
            GO TO 30
         END IF
         KL = KF
         INTVL( 1 ) = ALPHA
         INTVL( 2 ) = BETA
         INTVLCT( 1 ) = NALPHA
         INTVLCT( 2 ) = NBETA
      ELSE IF( IJOB.EQ.2 ) THEN
*
*        Check if some input intervals have "converged"
*
         CALL PDLAECV( 1, KF, KL, INTVL, INTVLCT, NVAL,
     $                 MAX( ABSTOL, PIVMIN ), RELTOL )
         IF( KF.GE.KL )
     $      GO TO 60
*
*        Compute upper bound on number of iterations needed
*
         ITMAX = INT( ( LOG( INTVL( 2 )-INTVL( 1 )+PIVMIN )-
     $           LOG( PIVMIN ) ) / LOG( TWO ) ) + 2
*
*        Iteration Loop
*
         DO 50 I = 1, ITMAX
            KLNEW = KL
            DO 40 J = KF, KL - 1
               K = 2*J
               MID = HALF*( INTVL( K-1 )+INTVL( K ) )
               IF( IEFLAG.EQ.0 ) THEN
                  CALL PDLAPDCT( MID, N, D, PIVMIN, NMID )
               ELSE IF( IEFLAG.EQ.1 ) THEN
                  CALL PDLAIECTB( MID, N, D, NMID )
               ELSE
                  CALL PDLAIECTL( MID, N, D, NMID )
               END IF
               LCNT = INTVLCT( K-1 )
               RCNT = INTVLCT( K )
               NMID = MIN( RCNT, MAX( LCNT, NMID ) )
*
*              Form New Interval(s)
*
               IF( NMID.EQ.LCNT ) THEN
                  INTVL( K-1 ) = MID
               ELSE IF( NMID.EQ.RCNT ) THEN
                  INTVL( K ) = MID
               ELSE IF( KLNEW.LT.MMAX+1 ) THEN
                  L = 2*KLNEW
                  INTVL( L-1 ) = MID
                  INTVL( L ) = INTVL( K )
                  INTVLCT( L-1 ) = NMID
                  INTVLCT( L ) = INTVLCT( K )
                  INTVL( K ) = MID
                  INTVLCT( K ) = NMID
                  KLNEW = KLNEW + 1
               ELSE
                  INFO = MMAX + 1
                  RETURN
               END IF
   40       CONTINUE
            KL = KLNEW
            CALL PDLAECV( 1, KF, KL, INTVL, INTVLCT, NVAL,
     $                    MAX( ABSTOL, PIVMIN ), RELTOL )
            IF( KF.GE.KL )
     $         GO TO 60
   50    CONTINUE
      END IF
   60 CONTINUE
      INFO = MAX( KL-KF, 0 )
      MOUT = KL - 1
      RETURN
*
*     End of PDLAEBZ
*
      END
*
*
      SUBROUTINE PDLAECV( IJOB, KF, KL, INTVL, INTVLCT, NVAL, ABSTOL,
     $                    RELTOL )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*
*     .. Scalar Arguments ..
      INTEGER            IJOB, KF, KL
      DOUBLE PRECISION   ABSTOL, RELTOL
*     ..
*     .. Array Arguments ..
      INTEGER            INTVLCT( * ), NVAL( * )
      DOUBLE PRECISION   INTVL( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAECV checks if the input intervals [ INTVL(2*i-1), INTVL(2*i) ],
*  i = KF, ... , KL-1, have "converged".
*  PDLAECV modifies KF to be the index of the last converged interval,
*  i.e., on output, all intervals [ INTVL(2*i-1), INTVL(2*i) ], i < KF,
*  have converged. Note that the input intervals may be reordered by
*  PDLAECV.
*
*  This is a SCALAPACK internal procedure and arguments are not checked
*  for unreasonable values.
*
*  Arguments
*  =========
*
*  IJOB    (input) INTEGER
*          Specifies the criterion for "convergence" of an interval.
*          = 0 : When an interval is narrower than ABSTOL, or than
*                RELTOL times the larger (in magnitude) endpoint, then
*                it is considered to have "converged".
*          = 1 : When an interval is narrower than ABSTOL, or than
*                RELTOL times the larger (in magnitude) endpoint, or if
*                the counts at the endpoints are identical to the counts
*                specified by NVAL ( see NVAL ) then the interval is
*                considered to have "converged".
*
*  KF      (input/output) INTEGER
*          On input, the index of the first input interval is 2*KF-1.
*          On output, the index of the last converged interval
*          is 2*KF-3.
*
*  KL      (input) INTEGER
*          The index of the last input interval is 2*KL-3.
*
*  INTVL   (input/output) DOUBLE PRECISION array, dimension (2*(KL-KF))
*          The endpoints of the intervals. INTVL(2*j-1) is the left
*          oendpoint f the j-th interval, and INTVL(2*j) is the right
*          endpoint of the j-th interval. The input intervals will,
*          in general, be reordered on output.
*          On input, INTVL contains the KL-KF input intervals.
*          On output, INTVL contains the converged intervals, 1 thru'
*          KF-1, and the unconverged intervals, KF thru' KL-1.
*
*  INTVLCT (input/output) INTEGER array, dimension (2*(KL-KF))
*          The counts at the endpoints of the intervals. INTVLCT(2*j-1)
*          is the count at the left endpoint of the j-th interval, i.e.,
*          the function value N(INTVL(2*j-1)), and INTVLCT(2*j) is the
*          count at the right endpoint of the j-th interval. This array
*          will, in general, be reordered on output.
*          See the comments in PDLAEBZ for more on the function N(w).
*
*  NVAL    (input/output) INTEGER array, dimension (2*(KL-KF))
*          The desired counts, N(w), at the endpoints of the
*          corresponding intervals.  This array will, in general,
*          be reordered on output.
*
*  ABSTOL  (input) DOUBLE PRECISION
*          The minimum (absolute) width of an interval. When an interval
*          is narrower than ABSTOL, or than RELTOL times the larger (in
*          magnitude) endpoint, then it is considered to be sufficiently
*          small, i.e., converged.
*          Note : This must be at least zero.
*
*  RELTOL  (input) DOUBLE PRECISION
*          The minimum relative width of an interval. When an interval
*          is narrower than ABSTOL, or than RELTOL times the larger (in
*          magnitude) endpoint, then it is considered to be sufficiently
*          small, i.e., converged.
*          Note : This should be at least radix*machine epsilon.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Local Scalars ..
      LOGICAL            CONDN
      INTEGER            I, ITMP1, ITMP2, J, K, KFNEW
      DOUBLE PRECISION   TMP1, TMP2, TMP3, TMP4
*     ..
*     .. Executable Statements ..
*
      KFNEW = KF
      DO 10 I = KF, KL - 1
         K = 2*I
         TMP3 = INTVL( K-1 )
         TMP4 = INTVL( K )
         TMP1 = ABS( TMP4-TMP3 )
         TMP2 = MAX( ABS( TMP3 ), ABS( TMP4 ) )
         CONDN = TMP1.LT.MAX( ABSTOL, RELTOL*TMP2 )
         IF( IJOB.EQ.0 )
     $      CONDN = CONDN .OR. ( ( INTVLCT( K-1 ).EQ.NVAL( K-1 ) ) .AND.
     $              INTVLCT( K ).EQ.NVAL( K ) )
         IF( CONDN ) THEN
            IF( I.GT.KFNEW ) THEN
*
*              Reorder Intervals
*
               J = 2*KFNEW
               TMP1 = INTVL( K-1 )
               TMP2 = INTVL( K )
               ITMP1 = INTVLCT( K-1 )
               ITMP2 = INTVLCT( K )
               INTVL( K-1 ) = INTVL( J-1 )
               INTVL( K ) = INTVL( J )
               INTVLCT( K-1 ) = INTVLCT( J-1 )
               INTVLCT( K ) = INTVLCT( J )
               INTVL( J-1 ) = TMP1
               INTVL( J ) = TMP2
               INTVLCT( J-1 ) = ITMP1
               INTVLCT( J ) = ITMP2
               IF( IJOB.EQ.0 ) THEN
                  ITMP1 = NVAL( K-1 )
                  NVAL( K-1 ) = NVAL( J-1 )
                  NVAL( J-1 ) = ITMP1
                  ITMP1 = NVAL( K )
                  NVAL( K ) = NVAL( J )
                  NVAL( J ) = ITMP1
               END IF
            END IF
            KFNEW = KFNEW + 1
         END IF
   10 CONTINUE
      KF = KFNEW
      RETURN
*
*     End of PDLAECV
*
      END
*
      SUBROUTINE PDLAPDCT( SIGMA, N, D, PIVMIN, COUNT )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*
*     .. Scalar Arguments ..
      INTEGER            COUNT, N
      DOUBLE PRECISION   PIVMIN, SIGMA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAPDCT counts the number of negative eigenvalues of (T - SIGMA I).
*  This implementation of the Sturm Sequence loop has conditionals in
*  the innermost loop to avoid overflow and determine the sign of a
*  floating point number. PDLAPDCT will be referred to as the "paranoid"
*  implementation of the Sturm Sequence loop.
*
*  This is a SCALAPACK internal procedure and arguments are not checked
*  for unreasonable values.
*
*  Arguments
*  =========
*
*  SIGMA   (input) DOUBLE PRECISION
*          The shift. PDLAPDCT finds the number of eigenvalues of T less
*          than or equal to SIGMA.
*
*  N       (input) INTEGER
*          The order of the tridiagonal matrix T. N >= 1.
*
*  D       (input) DOUBLE PRECISION array, dimension (2*N - 1)
*          Contains the diagonals and the squares of the off-diagonal
*          elements of the tridiagonal matrix T. These elements are
*          assumed to be interleaved in memory for better cache
*          performance. The diagonal entries of T are in the entries
*          D(1),D(3),...,D(2*N-1), while the squares of the off-diagonal
*          entries are D(2),D(4),...,D(2*N-2). To avoid overflow, the
*          matrix must be scaled so that its largest entry is no greater
*          than overflow**(1/2) * underflow**(1/4) in absolute value,
*          and for greatest accuracy, it should not be much smaller
*          than that.
*
*  PIVMIN  (input) DOUBLE PRECISION
*          The minimum absolute of a "pivot" in this "paranoid"
*          implementation of the Sturm sequence loop. This must be at
*          least max_j |e(j)^2| *safe_min, and at least safe_min, where
*          safe_min is at least the smallest number that can divide 1.0
*          without overflow.
*
*  COUNT   (output) INTEGER
*          The count of the number of eigenvalues of T less than or
*          equal to SIGMA.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   TMP
*     ..
*     .. Executable Statements ..
*
      TMP = D( 1 ) - SIGMA
      IF( ABS( TMP ).LE.PIVMIN )
     $   TMP = -PIVMIN
      COUNT = 0
      IF( TMP.LE.ZERO )
     $   COUNT = 1
      DO 10 I = 3, 2*N - 1, 2
         TMP = D( I ) - D( I-1 ) / TMP - SIGMA
         IF( ABS( TMP ).LE.PIVMIN )
     $      TMP = -PIVMIN
         IF( TMP.LE.ZERO )
     $      COUNT = COUNT + 1
   10 CONTINUE
*
      RETURN
*
*     End of PDLAPDCT
*
      END
