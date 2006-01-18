*
*
      SUBROUTINE PCLATMS( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX,
     $                    KL, KU, PACK, A, IA, JA, DESCA, ORDER, WORK,
     $                    LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          DIST, PACK, SYM
      INTEGER            IA, INFO, JA, KL, KU, LWORK, M, MODE, N, ORDER
      REAL               COND, DMAX
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), ISEED( 4 )
      REAL               D( * )
      COMPLEX            A( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*     PCLATMS generates random Hermitian matrices with specified
*     eigenvalues for testing SCALAPACK programs.
*
*     PCLATMS operates by applying the following sequence of
*     operations:
*
*       Set the diagonal to D, where D may be input or
*          computed according to MODE, COND, DMAX, and SYM
*          as described below.
*
*           Generate a dense M x N matrix by multiplying D on the left
*               and the right by random unitary matrices, then:
*
*           Reduce the bandwidth according to KL and KU, using
*           Householder transformations.
*           ### bandwidth reduction NOT SUPPORTED ###
*
*  Arguments
*  =========
*
*  M      - (global input) INTEGER
*           The number of rows of A. Not modified.
*
*  N      - (global input) INTEGER
*           The number of columns of A. Not modified.
*           ### M .ne. N unsupported
*
*  DIST   - (global input) CHARACTER*1
*           On entry, DIST specifies the type of distribution to be used
*           to generate the random eigen-/singular values.
*           'U' => UNIFORM( 0, 1 )  ( 'U' for uniform )
*           'S' => UNIFORM( -1, 1 ) ( 'S' for symmetric )
*           'N' => NORMAL( 0, 1 )   ( 'N' for normal )
*           Not modified.
*
*  ISEED  - (global input) INTEGER array, dimension ( 4 )
*           On entry ISEED specifies the seed of the random number
*           generator. They should lie between 0 and 4095 inclusive,
*           and ISEED(4) should be odd. The random number generator
*           uses a linear congruential sequence limited to small
*           integers, and so should produce machine independent
*           random numbers. The values of ISEED are changed on
*           exit, and can be used in the next call to CLATMS
*           to continue the same random number sequence.
*           Changed on exit.
*
*  SYM    - (global input) CHARACTER*1
*           If SYM='S' or 'H', the generated matrix is Hermitian, with
*             eigenvalues specified by D, COND, MODE, and DMAX; they
*             may be positive, negative, or zero.
*           If SYM='P', the generated matrix is Hermitian, with
*             eigenvalues (= singular values) specified by D, COND,
*             MODE, and DMAX; they will not be negative.
*           If SYM='N', the generated matrix is nonsymmetric, with
*             singular values specified by D, COND, MODE, and DMAX;
*             they will not be negative.
*           ### SYM = 'N' NOT SUPPORTED ###
*           Not modified.
*
*  D      - (local input/output) REAL array,
*           dimension ( MIN( M , N ) )
*           This array is used to specify the singular values or
*           eigenvalues of A (see SYM, above.)  If MODE=0, then D is
*           assumed to contain the singular/eigenvalues, otherwise
*           they will be computed according to MODE, COND, and DMAX,
*           and placed in D.
*           Modified if MODE is nonzero.
*
*  MODE   - (global input) INTEGER
*           On entry this describes how the singular/eigenvalues are to
*           be specified:
*           MODE = 0 means use D as input
*           MODE = 1 sets D(1)=1 and D(2:N)=1.0/COND
*           MODE = 2 sets D(1:N-1)=1 and D(N)=1.0/COND
*           MODE = 3 sets D(I)=COND**(-(I-1)/(N-1))
*           MODE = 4 sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND)
*           MODE = 5 sets D to random numbers in the range
*                    ( 1/COND , 1 ) such that their logarithms
*                    are uniformly distributed.
*           MODE = 6 set D to random numbers from same distribution
*                    as the rest of the matrix.
*           MODE < 0 has the same meaning as ABS(MODE), except that
*              the order of the elements of D is reversed.
*           Thus if MODE is positive, D has entries ranging from
*              1 to 1/COND, if negative, from 1/COND to 1,
*           If SYM='S' or 'H', and MODE is neither 0, 6, nor -6, then
*              the elements of D will also be multiplied by a random
*              sign (i.e., +1 or -1.)
*           Not modified.
*
*  COND   - (global input) REAL
*           On entry, this is used as described under MODE above.
*           If used, it must be >= 1. Not modified.
*
*  DMAX   - (global input) REAL
*           If MODE is neither -6, 0 nor 6, the contents of D, as
*           computed according to MODE and COND, will be scaled by
*           DMAX / max(abs(D(i))); thus, the maximum absolute eigen- or
*           singular value (which is to say the norm) will be abs(DMAX).
*           Note that DMAX need not be positive: if DMAX is negative
*           (or zero), D will be scaled by a negative number (or zero).
*           Not modified.
*
*  KL     - (global input) INTEGER
*           This specifies the lower bandwidth of the  matrix. For
*           example, KL=0 implies upper triangular, KL=1 implies upper
*           Hessenberg, and KL being at least M-1 means that the matrix
*           has full lower bandwidth.  KL must equal KU if the matrix
*           is Hermitian.
*           Not modified.
*           ### 1 <= KL < N-1 is NOT SUPPORTED ###
*
*  KU     - (global input) INTEGER
*           This specifies the upper bandwidth of the  matrix. For
*           example, KU=0 implies lower triangular, KU=1 implies lower
*           Hessenberg, and KU being at least N-1 means that the matrix
*           has full upper bandwidth.  KL must equal KU if the matrix
*           is Hermitian.
*           Not modified.
*           ### 1 <= KU < N-1 is NOT SUPPORTED ###
*
*  PACK   - (global input) CHARACTER*1
*           This specifies packing of matrix as follows:
*           'N' => no packing
*           ### PACK must be 'N'  all other options NOT SUPPORTED ###
*
*  A      - (local output) COMPLEX array
*           Global dimension (M, N), local dimension (MP, NQ)
*           On exit A is the desired test matrix.
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
*  ORDER  - (input) INTEGER
*           The number of reflectors used to define the orthogonal
*           matrix Q.  A = Q * D * Q'
*           Higher ORDER requires more computation and communication.
*
*  WORK   - (local input/output) COMPLEX array,
*           dimension (LWORK)
*
*  LWORK  - (local input) INTEGER dimension of WORK
*           LWORK >= SIZETMS as returned by PCLASIZESEP
*
*  INFO   - (global output) INTEGER
*           Error code.  On exit, INFO will be set to one of the
*           following values:
*             0 => normal return
*            -1 => M negative or unequal to N and SYM='S', 'H', or 'P'
*            -2 => N negative
*            -3 => DIST illegal string
*            -5 => SYM illegal string
*            -7 => MODE not in range -6 to 6
*            -8 => COND less than 1.0, and MODE neither -6, 0 nor 6
*           -10 => KL negative
*           -11 => KU negative, or SYM='S' or 'H' and KU not equal to KL
*           -16 => DESCA is inconsistent
*           -17 => ORDER not in the range 0 to N inclusive
*            1  => Error return from SLATM1
*            2  => Cannot scale to DMAX (max. sing. value is 0)
*            3  => Error return from PCLAGHE
*
*-----------------------------------------------------------------------
*
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IDIST, IINFO, IPACK, IRSIGN, ISYM, LLB,
     $                   MNMIN, MYCOL, MYROW, NP, NPCOL, NPROW, NQ
      REAL               ALPHA, TEMP
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      EXTERNAL           LSAME, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, CLASET, PCHK1MAT,
     $                   PCLAGHE, PXERBLA, SLATM1, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     1)      Decode and Test the input parameters.
*             Initialize flags & seed.
*
*
      INFO = 0
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      IF( ( MYROW.GE.NPROW .OR. MYROW.LT.0 ) .OR.
     $    ( MYCOL.GE.NPCOL .OR. MYCOL.LT.0 ) )RETURN
*
      NP = NUMROC( N, DESCA( MB_ ), MYROW, 0, NPROW )
      NQ = NUMROC( N, DESCA( NB_ ), MYCOL, 0, NPCOL )
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Decode DIST
*
      IF( LSAME( DIST, 'U' ) ) THEN
         IDIST = 1
      ELSE IF( LSAME( DIST, 'S' ) ) THEN
         IDIST = 2
      ELSE IF( LSAME( DIST, 'N' ) ) THEN
         IDIST = 3
      ELSE
         IDIST = -1
      END IF
*
*     Decode SYM
*
      IF( LSAME( SYM, 'N' ) ) THEN
         ISYM = 1
         IRSIGN = 0
      ELSE IF( LSAME( SYM, 'P' ) ) THEN
         ISYM = 2
         IRSIGN = 0
      ELSE IF( LSAME( SYM, 'S' ) ) THEN
         ISYM = 2
         IRSIGN = 1
      ELSE IF( LSAME( SYM, 'H' ) ) THEN
         ISYM = 2
         IRSIGN = 1
      ELSE
         ISYM = -1
      END IF
*
*     Decode PACK
*
      IF( LSAME( PACK, 'N' ) ) THEN
         IPACK = 0
      ELSE
         IPACK = 1
      END IF
*
*     Set certain internal parameters
*
      MNMIN = MIN( M, N )
      LLB = MIN( KL, M-1 )
*
      IF( ORDER.EQ.0 )
     $   ORDER = N
*
*     Set INFO if an error
*
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 1600+CTXT_ )
      ELSE
         CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 16, INFO )
         IF( INFO.EQ.0 ) THEN
            IF( M.NE.N .AND. ISYM.NE.1 ) THEN
               INFO = -2
            ELSE IF( IDIST.EQ.-1 ) THEN
               INFO = -3
            ELSE IF( ISYM.EQ.-1 ) THEN
               INFO = -5
            ELSE IF( ABS( MODE ).GT.6 ) THEN
               INFO = -7
            ELSE IF( ( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) .AND. COND.LT.
     $               ONE ) THEN
               INFO = -8
            ELSE IF( KL.LT.0 ) THEN
               INFO = -10
            ELSE IF( KU.LT.0 .OR. ( ISYM.NE.1 .AND. KL.NE.KU ) ) THEN
               INFO = -11
            ELSE IF( ( ORDER.LT.0 ) .OR. ( ORDER.GT.N ) ) THEN
               INFO = -17
            END IF
         END IF
         CALL PCHK1MAT( M, 1, N, 2, IA, JA, DESCA, 16, 0, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
*     Check for unsupported features
*
      IF( ISYM.NE.2 ) THEN
         INFO = -5
      ELSE IF( IPACK.NE.0 ) THEN
         INFO = -12
      ELSE IF( KL.GT.0 .AND. KL.LT.M-1 ) THEN
         INFO = -10
      ELSE IF( KU.GT.0 .AND. KU.LT.N-1 ) THEN
         INFO = -11
      ELSE IF( LLB.NE.0 .AND. LLB.NE.M-1 ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PCLATMS', -INFO )
         RETURN
      END IF
*
*     Initialize random number generator
*
      DO 10 I = 1, 4
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
   10 CONTINUE
*
      IF( MOD( ISEED( 4 ), 2 ).NE.1 )
     $   ISEED( 4 ) = ISEED( 4 ) + 1
*
*     2)      Set up D  if indicated.
*
*             Compute D according to COND and MODE
*
      CALL SLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, IINFO )
*
      IF( IINFO.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
*
*
      IF( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) THEN
*
*        Scale by DMAX
*
         TEMP = ABS( D( 1 ) )
         DO 20 I = 2, MNMIN
            TEMP = MAX( TEMP, ABS( D( I ) ) )
   20    CONTINUE
*
         IF( TEMP.GT.ZERO ) THEN
            ALPHA = DMAX / TEMP
         ELSE
            INFO = 2
            RETURN
         END IF
*
         CALL SSCAL( MNMIN, ALPHA, D, 1 )
*
      END IF
*
      CALL CLASET( 'A', NP, NQ, CZERO, CZERO, A, DESCA( LLD_ ) )
*
*     Hermitian -- A = U D U'
*
      CALL PCLAGHE( M, LLB, D, A, IA, JA, DESCA, ISEED, ORDER, WORK,
     $              LWORK, IINFO )
*
      RETURN
*
*     End of PCLATMS
*
      END
