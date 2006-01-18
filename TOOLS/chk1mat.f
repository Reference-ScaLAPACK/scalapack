      SUBROUTINE CHK1MAT( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     $                    DESCAPOS0, INFO )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            DESCAPOS0, IA, INFO, JA, MA, MAPOS0, NA, NAPOS0
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
*     ..
*
*  Purpose
*  =======
*
*  CHK1MAT checks that the values associated with one distributed matrix
*  make sense from a local viewpoint
*
*  Arguments
*  =========
*
*  MA      (global input) INTEGER
*          The number or matrix rows of A being operated on.
*
*  MAPOS0  (global input) INTEGER
*          Where in the calling routine's parameter list MA appears.
*
*  NA      (global input) INTEGER
*          The number of matrix columns of A being operated on.
*
*  NAPOS0  (global input) INTEGER
*          Where in the calling routine's parameter list NA appears.
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
*  DESCAPOS0 (global input) INTEGER
*          Where in the calling routine's parameter list DESCA
*          appears.  Note that we assume IA and JA are respectively 2
*          and 1 entries behind DESCA.
*
*  INFO    (local input/local output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            DESCMULT, BIGNUM
      PARAMETER          ( DESCMULT = 100, BIGNUM = DESCMULT*DESCMULT )
*     ..
*     .. Local Scalars ..
      INTEGER            DESCAPOS, IAPOS, JAPOS, MAPOS, NAPOS, MYCOL,
     $                   MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX
*     ..
*     .. Executable Statements ..
*
*     Want to find errors with MIN( ), so if no error, set it to a big
*     number.  If there already is an error, multiply by the the des-
*     criptor multiplier
*
      IF( INFO.GE.0 ) THEN
         INFO = BIGNUM
      ELSE IF( INFO.LT.-DESCMULT ) THEN
         INFO = -INFO
      ELSE
         INFO = -INFO * DESCMULT
      END IF
*
*     Figure where in parameter list each parameter was, factoring in
*     descriptor multiplier
*
      MAPOS = MAPOS0 * DESCMULT
      NAPOS = NAPOS0 * DESCMULT
      IAPOS = (DESCAPOS0-2) * DESCMULT
      JAPOS = (DESCAPOS0-1) * DESCMULT
      DESCAPOS = DESCAPOS0 * DESCMULT
*
*     Get grid parameters
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
*     Check that matrix values make sense from local viewpoint
*
      IF( DESCA( DTYPE_ ) .NE. BLOCK_CYCLIC_2D ) THEN
         INFO = MIN( INFO, DESCAPOS+DTYPE_ )
      ELSE IF( MA.LT.0 ) THEN
         INFO = MIN( INFO, MAPOS )
      ELSE IF( NA.LT.0 ) THEN
         INFO = MIN( INFO, NAPOS )
      ELSE IF( IA.LT.1 ) THEN
         INFO = MIN( INFO, IAPOS )
      ELSE IF( JA.LT.1 ) THEN
         INFO = MIN( INFO, JAPOS )
      ELSE IF( DESCA( MB_ ).LT.1 ) THEN
         INFO = MIN( INFO, DESCAPOS+MB_ )
      ELSE IF( DESCA( NB_ ).LT.1 ) THEN
         INFO = MIN( INFO, DESCAPOS+NB_ )
      ELSE IF( DESCA( RSRC_ ).LT.0 .OR. DESCA( RSRC_ ).GE.NPROW ) THEN
         INFO = MIN( INFO, DESCAPOS+RSRC_ )
      ELSE IF( DESCA( CSRC_ ).LT.0 .OR. DESCA( CSRC_ ).GE.NPCOL ) THEN
         INFO = MIN( INFO, DESCAPOS+CSRC_ )
      ELSE IF( DESCA( LLD_ ).LT.1 ) THEN
            INFO = MIN( INFO, DESCAPOS+LLD_ )
      ELSE IF( DESCA( LLD_ ) .LT.
     $         NUMROC( DESCA( M_ ), DESCA( MB_ ), MYROW, DESCA(RSRC_),
     $                 NPROW ) ) THEN
         IF( NUMROC( DESCA( N_ ), DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $               NPCOL ) .GT. 0 )
     $      INFO = MIN( INFO, DESCAPOS+LLD_ )
      END IF
*
      IF( MA.EQ.0 .OR. NA.EQ.0 ) THEN
*
*        NULL matrix, relax some checks
*
         IF( DESCA(M_).LT.0 )
     $      INFO = MIN( INFO, DESCAPOS+M_ )
         IF( DESCA(N_).LT.0 )
     $      INFO = MIN( INFO, DESCAPOS+N_ )
*
      ELSE
*
*        more rigorous checks for non-degenerate matrices
*
         IF( DESCA( M_ ).LT.1 ) THEN
            INFO = MIN( INFO, DESCAPOS+M_ )
         ELSE IF( DESCA( N_ ).LT.1 ) THEN
            INFO = MIN( INFO, DESCAPOS+N_ )
         ELSE
            IF( IA.GT.DESCA( M_ ) ) THEN
               INFO = MIN( INFO, IAPOS )
            ELSE IF( JA.GT.DESCA( N_ ) ) THEN
               INFO = MIN( INFO, JAPOS )
            ELSE
               IF( IA+MA-1.GT.DESCA( M_ ) )
     $            INFO = MIN( INFO, MAPOS )
               IF( JA+NA-1.GT.DESCA( N_ ) )
     $            INFO = MIN( INFO, NAPOS )
            END IF
         END IF
*
      END IF
*
*     Prepare output: set info = 0 if no error, and divide by
*     DESCMULT if error is not in a descriptor entry
*
      IF( INFO.EQ.BIGNUM ) THEN
         INFO = 0
      ELSE IF( MOD( INFO, DESCMULT ).EQ.0 ) THEN
         INFO = -INFO / DESCMULT
      ELSE
         INFO = -INFO
      END IF
*
      RETURN
*
*     End CHK1MAT
*
      END
