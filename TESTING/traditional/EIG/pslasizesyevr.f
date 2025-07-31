      SUBROUTINE PSLASIZESYEVR( WKNOWN, RANGE, N, DESCA, VL, VU, IL, IU,
     $                          ISEED, WIN, MAXSIZE, VECSIZE, VALSIZE )
*
*  -- ScaLAPACK routine (@(MODE)version *TBA*) --
*     University of California, Berkeley and
*     University of Tennessee, Knoxville. 
*     October 21, 2006
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      LOGICAL            WKNOWN
      CHARACTER          RANGE
      INTEGER            IL, IU, MAXSIZE, N, VALSIZE, VECSIZE
      REAL               VL, VU

*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), ISEED( 4 )
      REAL               WIN( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLASIZESYEVR computes the amount of memory needed by PSSYEVR
*  to ensure:
*    1)  Orthogonal Eigenvectors
*    2)  Eigenpairs with small residual norms
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
*  VL      (global input/output ) REAL            
*          If RANGE='V', the lower bound of the interval to be searched
*          for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
*          If VL > VU, RANGE='V' and WKNOWN = .TRUE., VL is set
*          to a random value near an entry in WIN
*
*  VU      (global input/output ) REAL            
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
*  WIN     (global input) REAL             array, dimension (N)
*          If WKNOWN=1, WIN contains the eigenvalues of the matrix.
*
*  MAXSIZE (global output) INTEGER
*          Workspace required to guarantee that PSSYEVR will return
*          orthogonal eigenvectors.  IF WKNOWN=0, MAXSIZE is set to a
*          a value which guarantees orthogonality no matter what the
*          spectrum is.  If WKNOWN=1, MAXSIZE is set to a value which
*          guarantees orthogonality on a matrix with eigenvalues given
*          by WIN.
*
*  VECSIZE (global output) INTEGER
*          Workspace required to guarantee that PSSYEVR
*          will compute eigenvectors.
*
*  VALSIZE (global output) INTEGER
*          Workspace required to guarantee that PSSYEVR
*          will compute eigenvalues.
*
*
*     .. Parameters ..
      INTEGER            CTXT_, MB_
      PARAMETER          ( CTXT_ = 2, MB_ = 5 )
      REAL               TWENTY
      PARAMETER          ( TWENTY = 20.0E0 )
*     ..
*     .. Local Scalars ..
*
      INTEGER            ILMIN, IUMAX, 
     $                   MQ0, MYCOL, MYIL, MYIU, MYROW, NB, NEIG, NN,
     $                   NP0, NPCOL, NPROW
      REAL               ANORM, EPS, SAFMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, NUMROC
      REAL               SLARAN, PSLAMCH
      EXTERNAL           LSAME, ICEIL, NUMROC, SLARAN, PSLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL, INT, MAX

*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      EPS = PSLAMCH( DESCA( CTXT_ ), 'Precision' )
      SAFMIN = PSLAMCH( DESCA( CTXT_ ), 'Safe Minimum' )
      NB = DESCA( MB_ )
      NN = MAX( N, NB, 2 )
      NP0 = NUMROC( NN, NB, 0, 0, NPROW )

      VALSIZE = 3 + 5*N + MAX( 12*NN, NB*( NP0+1 ) )

      IF( WKNOWN ) THEN
         ANORM = SAFMIN / EPS
         IF( N.GE.1 )
     $      ANORM = MAX( ABS( WIN( 1 ) ), ABS( WIN( N ) ), ANORM )
         IF( LSAME( RANGE, 'I' ) ) THEN
            IF( IL.LT.0 )
     $         IL = INT( SLARAN( ISEED )*REAL( N ) ) + 1
            IF( IU.LT.0 )
     $         IU = INT( SLARAN( ISEED )*REAL( N-IL ) ) + IL
            IF( N.EQ.0 )
     $         IU = 0
         ELSE IF( LSAME( RANGE, 'V' ) ) THEN
            IF( VL.GT.VU ) THEN
               MYIL = INT( SLARAN( ISEED )*REAL( N ) ) + 1
               MYIU = INT( SLARAN( ISEED )*REAL( N-MYIL ) ) + MYIL
               VL = WIN( MYIL ) - TWENTY*EPS*ABS( WIN( MYIL ) )
               VU = WIN( MYIU ) + TWENTY*EPS*ABS( WIN( MYIU ) )
               VU = MAX( VU, VL+EPS*TWENTY*ABS( VL )+SAFMIN )
            END IF
         END IF
*
      END IF
      IF( LSAME( RANGE, 'V' ) ) THEN
*        We do not know how many eigenvalues will be computed
         ILMIN = 1
         IUMAX = N
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
*
      VECSIZE = 3 + 5*N + MAX( 18*NN, NP0*MQ0+2*NB*NB ) + 
     $          (2 + ICEIL( NEIG, NPROW*NPCOL ))*NN

      VALSIZE = MAX(3, VALSIZE)
      VECSIZE = MAX(3, VECSIZE)
      MAXSIZE = VECSIZE
*
      RETURN
*
*     End of PSLASIZESYEVR
*
      END
