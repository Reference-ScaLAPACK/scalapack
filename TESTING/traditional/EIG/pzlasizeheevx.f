*
*
      SUBROUTINE PZLASIZEHEEVX( WKNOWN, RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN, MAXSIZE, VECSIZE, VALSIZE )
*
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      LOGICAL            WKNOWN
      CHARACTER          RANGE
      INTEGER            IL, IU, MAXSIZE, N, VALSIZE, VECSIZE
      DOUBLE PRECISION   VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), ISEED( 4 )
      DOUBLE PRECISION   WIN( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLASIZEHEEVX computes the amount of memory needed by PZHEEVX
*  to ensure:
*    1)  Orthogonal Eigenvectors
*    2)  Eigenvectors
*    3)  Eigenvalues
*
*  Arguments
*  =========
*
*  WKNOWN  (global input) INTEGER
*          .FALSE.:  WIN does not contain the eigenvalues
*          .TRUE.:   WIN does contain the eigenvalues
*
*  RANGE   (global input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the interval [VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  N       (global input) INTEGER
*          Size of the matrix to be tested.  (global size)
*
*  DESCA   (global input) INTEGER array dimension ( DLEN_ )
*
*  VL      (global input/output ) DOUBLE PRECISION
*          If RANGE='V', the lower bound of the interval to be searched
*          for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
*          If VL > VU, RANGE='V' and WKNOWN = .TRUE., VL is set
*          to a random value near an entry in WIN
*
*  VU      (global input/output ) DOUBLE PRECISION
*          If RANGE='V', the upper bound of the interval to be searched
*          for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
*          If VL > VU, RANGE='V' and WKNOWN = .TRUE., VU is set
*          to a random value near an entry in WIN
*
*  IL      (global input/output ) INTEGER
*          If RANGE='I', the index (from smallest to largest) of the
*          smallest eigenvalue to be returned.  IL >= 1.
*          Not referenced if RANGE = 'A' or 'V'.
*          If IL < 0, RANGE='I' and WKNOWN = .TRUE., IL is set
*          to a random value from 1 to N
*
*  IU      (global input/output ) INTEGER
*          If RANGE='I', the index (from smallest to largest) of the
*          largest eigenvalue to be returned.  min(IL,N) <= IU <= N.
*          Not referenced if RANGE = 'A' or 'V'.
*          If IU < 0, RANGE='I' and WKNOWN = .TRUE., IU is set
*          to a random value from IL to N
*
*  ISEED   (global input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*          ISEED is not touched unless IL, IU, VL or VU are modified.
*
*  WIN     (global input) DOUBLE PRECISION array, dimension (N)
*          If WKNOWN=1, WIN contains the eigenvalues of the matrix.
*
*  MAXSIZE (global output) INTEGER
*          Workspace required to guarantee that PZHEEVX will return
*          orthogonal eigenvectors.  IF WKNOWN=0, MAXSIZE is set to a
*          a value which guarantees orthogonality no matter what the
*          spectrum is.  If WKNOWN=1, MAXSIZE is set to a value which
*          guarantees orthogonality on a matrix with eigenvalues given
*          by WIN.
*
*  VECSIZE (global output) INTEGER
*          Workspace required to guarantee that PZHEEVX
*          will compute eigenvectors.
*
*  VALSIZE (global output) INTEGER
*          Workspace required to guarantee that PZHEEVX
*          will compute eigenvalues.
*
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   TWENTY
      PARAMETER          ( TWENTY = 20.0D0 )
*     ..
*     .. Local Scalars ..
*
      INTEGER            CLUSTERSIZE, I, ILMIN, IUMAX, MAXCLUSTERSIZE,
     $                   MQ0, MYCOL, MYIL, MYIU, MYROW, NB, NEIG, NN,
     $                   NP0, NPCOL, NPROW
      DOUBLE PRECISION   ANORM, EPS, ORFAC, SAFMIN, VLMIN, VUMAX
*     ..
*     .. External Functions ..
*
*
      LOGICAL            LSAME
      INTEGER            ICEIL, NUMROC
      DOUBLE PRECISION   DLARAN, PDLAMCH
      EXTERNAL           LSAME, ICEIL, NUMROC, DLARAN, PDLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      ORFAC = 1.0D-3
*
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      EPS = PDLAMCH( DESCA( CTXT_ ), 'Precision' )
      SAFMIN = PDLAMCH( DESCA( CTXT_ ), 'Safe Minimum' )
      NB = DESCA( MB_ )
      NN = MAX( N, NB, 2 )
      NP0 = NUMROC( NN, NB, 0, 0, NPROW )
*
      VALSIZE = 5*NN + 4*N
*
      IF( WKNOWN ) THEN
         ANORM = SAFMIN / EPS
         IF( N.GE.1 )
     $      ANORM = MAX( ABS( WIN( 1 ) ), ABS( WIN( N ) ), ANORM )
*
         IF( LSAME( RANGE, 'I' ) ) THEN
            IF( IL.LT.0 )
     $         IL = INT( DLARAN( ISEED )*DBLE( N ) ) + 1
            IF( IU.LT.0 )
     $         IU = INT( DLARAN( ISEED )*DBLE( N-IL ) ) + IL
            IF( N.EQ.0 )
     $         IU = 0
         ELSE IF( LSAME( RANGE, 'V' ) ) THEN
            IF( VL.GT.VU ) THEN
               MYIL = INT( DLARAN( ISEED )*DBLE( N ) ) + 1
               MYIU = INT( DLARAN( ISEED )*DBLE( N-MYIL ) ) + MYIL
               VL = WIN( MYIL ) + TWENTY*EPS*ABS( WIN( MYIL ) )
               VU = WIN( MYIU ) + TWENTY*EPS*ABS( WIN( MYIU ) )
               VU = MAX( VU, VL+EPS*TWENTY*ABS( VL )+SAFMIN )
            END IF
         END IF
*
      END IF
      IF( LSAME( RANGE, 'V' ) ) THEN
*
*     Compute ILMIN, IUMAX (based on VL, VU and WIN)
*
         IF( WKNOWN ) THEN
            VLMIN = VL - TWENTY*EPS*ANORM
            VUMAX = VU + TWENTY*EPS*ANORM
            ILMIN = 1
            IUMAX = 0
            DO 10 I = 1, N
               IF( WIN( I ).LT.VLMIN )
     $            ILMIN = ILMIN + 1
               IF( WIN( I ).LT.VUMAX )
     $            IUMAX = IUMAX + 1
   10       CONTINUE
         ELSE
            ILMIN = 1
            IUMAX = N
         END IF
      ELSE IF( LSAME( RANGE, 'I' ) ) THEN
         ILMIN = IL
         IUMAX = IU
      ELSE IF( LSAME( RANGE, 'A' ) ) THEN
         ILMIN = 1
         IUMAX = N
      END IF
*
      NEIG = IUMAX - ILMIN + 1
*
      MQ0 = NUMROC( MAX( NEIG, NB, 2 ), NB, 0, 0, NPCOL )
      VECSIZE = 4*N + MAX( 5*NN, NP0*MQ0 ) +
     $          ICEIL( NEIG, NPROW*NPCOL )*NN
*
      IF( WKNOWN ) THEN
         CLUSTERSIZE = 1
         MAXCLUSTERSIZE = 1
         DO 20 I = ILMIN + 1, IUMAX
            IF( ( WIN( I )-WIN( I-1 ) ).LT.ORFAC*2*ANORM ) THEN
               CLUSTERSIZE = CLUSTERSIZE + 1
               IF( CLUSTERSIZE.GT.MAXCLUSTERSIZE )
     $            MAXCLUSTERSIZE = CLUSTERSIZE
            ELSE
               CLUSTERSIZE = 1
            END IF
   20    CONTINUE
         IF( CLUSTERSIZE.GT.MAXCLUSTERSIZE )
     $      MAXCLUSTERSIZE = CLUSTERSIZE
      ELSE
         MAXCLUSTERSIZE = N
      END IF
*
      MAXSIZE = VECSIZE + MAX( ( MAXCLUSTERSIZE-1 ), 0 )*N
*
*
      RETURN
*
*     End of PZLASIZEHEEVX
*
      END
