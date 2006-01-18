      SUBROUTINE PCHK1MAT( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     $                     DESCAPOS0, NEXTRA, EX, EXPOS, INFO )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            DESCAPOS0, IA, INFO, JA, MA, MAPOS0, NA,
     $                   NAPOS0, NEXTRA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), EX( NEXTRA ), EXPOS( NEXTRA )
*     ..
*
*  Purpose
*  =======
*
*  PCHK1MAT checks that the values associated with one distributed
*  matrix are consistant across the entire process grid.
*
*  Notes
*  =====
*
*  This routine checks that all values are the same across the grid.
*  It does no local checking; it is therefore legal to abuse the
*  definitions of the non-descriptor arguments, i.e., if the routine
*  you are checking does not possess a MA value, you may pass some
*  other integer that must be global into this argument instead.
*
*  Arguments
*  =========
*
*  MA      (global input) INTEGER
*          The global number of matrix rows of A being operated on.
*
*  MAPOS0  (global input) INTEGER
*          Where in the calling routine's parameter list MA appears.
*
*  NA      (global input) INTEGER
*          The global number of matrix columns of A being operated on.
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
*  NEXTRA  (global input) INTEGER
*          The number of extra parameters (i.e., besides the ones
*          above) to check.  NEXTRA <= LDW - 11.
*
*  EX      (local input) INTEGER array of dimension (NEXTRA)
*          The values of these extra parameters
*
*  EXPOS   (local input) INTEGER array of dimension (NEXTRA)
*          The parameter list positions of these extra values.
*
*  INFO    (local input/global output) INTEGER
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
      INTEGER            BIGNUM, DESCMULT, LDW
      PARAMETER          ( DESCMULT = 100, BIGNUM = DESCMULT * DESCMULT,
     $                     LDW = 25 )
*     ..
*     .. Local Scalars ..
      INTEGER            DESCPOS, K
*     ..
*     .. Local Arrays ..
      INTEGER            IWORK( LDW, 2 ), IWORK2( LDW )
*     ..
*     .. External Subroutines ..
      EXTERNAL           GLOBCHK
*     ..
*     .. Executable Statements ..
*
*     Want to find errors with MIN( ), so if no error, set it to a big
*     number. If there already is an error, multiply by the the
*     descriptor multiplier.
*
      IF( INFO.GE.0 ) THEN
         INFO = BIGNUM
      ELSE IF( INFO.LT.-DESCMULT ) THEN
         INFO = -INFO
      ELSE
         INFO = -INFO * DESCMULT
      END IF
*
*     Pack values and their positions in the parameter list, factoring
*     in the descriptor multiplier
*
      IWORK( 1, 1 ) = MA
      IWORK( 1, 2 ) = MAPOS0 * DESCMULT
      IWORK( 2, 1 ) = NA
      IWORK( 2, 2 ) = NAPOS0 * DESCMULT
      IWORK( 3, 1 ) = IA
      IWORK( 3, 2 ) = (DESCAPOS0-2) * DESCMULT
      IWORK( 4, 1 ) = JA
      IWORK( 4, 2 ) = (DESCAPOS0-1) * DESCMULT
      DESCPOS = DESCAPOS0 * DESCMULT
*
      IWORK(  5, 1 ) = DESCA( DTYPE_ )
      IWORK(  5, 2 ) = DESCPOS + DTYPE_
      IWORK(  6, 1 ) = DESCA( M_ )
      IWORK(  6, 2 ) = DESCPOS + M_
      IWORK(  7, 1 ) = DESCA( N_ )
      IWORK(  7, 2 ) = DESCPOS + N_
      IWORK(  8, 1 ) = DESCA( MB_ )
      IWORK(  8, 2 ) = DESCPOS + MB_
      IWORK(  9, 1 ) = DESCA( NB_ )
      IWORK(  9, 2 ) = DESCPOS + NB_
      IWORK( 10, 1 ) = DESCA( RSRC_ )
      IWORK( 10, 2 ) = DESCPOS + RSRC_
      IWORK( 11, 1 ) = DESCA( CSRC_ )
      IWORK( 11, 2 ) = DESCPOS + CSRC_
*
      IF( NEXTRA.GT.0 ) THEN
         DO 10 K = 1, NEXTRA
            IWORK( 11+K, 1 ) = EX( K )
            IWORK( 11+K, 2 ) = EXPOS( K )
   10    CONTINUE
      END IF
      K = 11 + NEXTRA
*
*     Get the smallest error detected anywhere (BIGNUM if no error)
*
      CALL GLOBCHK( DESCA( CTXT_ ), K, IWORK, LDW, IWORK2, INFO )
*
*     Prepare output: set info = 0 if no error, and divide by DESCMULT if
*     error is not in a descriptor entry
*
      IF( INFO .EQ. BIGNUM ) THEN
         INFO = 0
      ELSE IF( MOD( INFO, DESCMULT ) .EQ. 0 ) THEN
         INFO = -INFO / DESCMULT
      ELSE
         INFO = -INFO
      END IF
*
      RETURN
*
*     End of PCHK1MAT
*
      END
*
      SUBROUTINE PCHK2MAT( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     $                     DESCAPOS0, MB, MBPOS0, NB, NBPOS0, IB, JB,
     $                     DESCB, DESCBPOS0, NEXTRA, EX, EXPOS, INFO )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            DESCAPOS0, DESCBPOS0, IA, IB, INFO, JA, JB, MA,
     $                   MAPOS0, MB, MBPOS0, NA, NAPOS0, NB, NBPOS0,
     $                   NEXTRA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( 8 ), EX( NEXTRA ),
     $                   EXPOS( NEXTRA )
*     ..
*
*  Purpose
*  =======
*
*  PCHK2MAT checks that the values associated with two distributed
*  matrices are consistant across the entire process grid.
*
*  Notes
*  =====
*
*  This routine checks that all values are the same across the grid.
*  It does no local checking; it is therefore legal to abuse the
*  definitions of the non-descriptor arguments, i.e., if the routine
*  you are checking does not possess a MA value, you may pass some
*  other integer that must be global into this argument instead.
*
*  Arguments
*  =========
*
*  MA      (global input) INTEGER
*          The global number of matrix rows of A being operated on.
*
*  MAPOS0  (global input) INTEGER
*          Where in the calling routine's parameter list MA appears.
*
*  NA      (global input) INTEGER
*          The global number of matrix columns of A being operated on.
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
*  MB      (global input) INTEGER
*          The global number of matrix rows of B being operated on.
*
*  MBPOS0  (global input) INTEGER
*          Where in the calling routine's parameter list MB appears.
*
*  NB      (global input) INTEGER
*          The global number of matrix columns of B being operated on.
*
*  NBPOS0  (global input) INTEGER
*          Where in the calling routine's parameter list NB appears.
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
*  DESCBPOS0 (global input) INTEGER
*          Where in the calling routine's parameter list DESCB
*          appears. Note that we assume IB and JB are respectively 2
*          and 1 entries behind DESCB.
*
*  NEXTRA  (global input) INTEGER
*          The number of extra parameters (i.e., besides the ones
*          above) to check.  NEXTRA <= LDW - 22.
*
*  EX      (local input) INTEGER array of dimension (NEXTRA)
*          The values of these extra parameters
*
*  EXPOS   (local input) INTEGER array of dimension (NEXTRA)
*          The parameter list positions of these extra values.
*
*  INFO    (local input/global output) INTEGER
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
      INTEGER            DESCMULT, BIGNUM, LDW
      PARAMETER          ( DESCMULT = 100, BIGNUM = DESCMULT * DESCMULT,
     $                     LDW = 35 )
*     ..
*     .. Local Scalars ..
      INTEGER            K, DESCPOS
*     ..
*     .. Local Arrays ..
      INTEGER            IWORK( LDW, 2 ), IWORK2( LDW )
*     ..
*     .. External Subroutines ..
      EXTERNAL           GLOBCHK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     Want to find errors with MIN( ), so if no error, set it to a big
*     number. If there already is an error, multiply by the the
*     descriptor multiplier.
*
      IF( INFO.GE.0 ) THEN
         INFO = BIGNUM
      ELSE IF( INFO.LT.-DESCMULT ) THEN
         INFO = -INFO
      ELSE
         INFO = -INFO * DESCMULT
      END IF
*
*     Pack values and their positions in the parameter list, factoring
*     in the descriptor multiplier
*
      IWORK( 1, 1 ) = MA
      IWORK( 1, 2 ) = MAPOS0 * DESCMULT
      IWORK( 2, 1 ) = NA
      IWORK( 2, 2 ) = NAPOS0 * DESCMULT
      IWORK( 3, 1 ) = IA
      IWORK( 3, 2 ) = (DESCAPOS0-2) * DESCMULT
      IWORK( 4, 1 ) = JA
      IWORK( 4, 2 ) = (DESCAPOS0-1) * DESCMULT
      DESCPOS = DESCAPOS0 * DESCMULT
*
      IWORK(  5, 1 ) = DESCA( DTYPE_ )
      IWORK(  5, 2 ) = DESCPOS + DTYPE_
      IWORK(  6, 1 ) = DESCA( M_ )
      IWORK(  6, 2 ) = DESCPOS + M_
      IWORK(  7, 1 ) = DESCA( N_ )
      IWORK(  7, 2 ) = DESCPOS + N_
      IWORK(  8, 1 ) = DESCA( MB_ )
      IWORK(  8, 2 ) = DESCPOS + MB_
      IWORK(  9, 1 ) = DESCA( NB_ )
      IWORK(  9, 2 ) = DESCPOS + NB_
      IWORK( 10, 1 ) = DESCA( RSRC_ )
      IWORK( 10, 2 ) = DESCPOS + RSRC_
      IWORK( 11, 1 ) = DESCA( CSRC_ )
      IWORK( 11, 2 ) = DESCPOS + CSRC_
*
      IWORK( 12, 1 ) = MB
      IWORK( 12, 2 ) = MBPOS0 * DESCMULT
      IWORK( 13, 1 ) = NB
      IWORK( 13, 2 ) = NBPOS0 * DESCMULT
      IWORK( 14, 1 ) = IB
      IWORK( 14, 2 ) = (DESCBPOS0-2) * DESCMULT
      IWORK( 15, 1 ) = JB
      IWORK( 15, 2 ) = (DESCBPOS0-1) * DESCMULT
      DESCPOS = DESCBPOS0 * DESCMULT
*
      IWORK( 16, 1 ) = DESCB( DTYPE_ )
      IWORK( 16, 2 ) = DESCPOS + DTYPE_
      IWORK( 17, 1 ) = DESCB( M_ )
      IWORK( 17, 2 ) = DESCPOS + M_
      IWORK( 18, 1 ) = DESCB( N_ )
      IWORK( 18, 2 ) = DESCPOS + N_
      IWORK( 19, 1 ) = DESCB( MB_ )
      IWORK( 19, 2 ) = DESCPOS + MB_
      IWORK( 20, 1 ) = DESCB( NB_ )
      IWORK( 20, 2 ) = DESCPOS + NB_
      IWORK( 21, 1 ) = DESCB( RSRC_ )
      IWORK( 21, 2 ) = DESCPOS + RSRC_
      IWORK( 22, 1 ) = DESCB( CSRC_ )
      IWORK( 22, 2 ) = DESCPOS + CSRC_
*
      IF( NEXTRA.GT.0 ) THEN
         DO 10 K = 1, NEXTRA
            IWORK( 22+K, 1 ) = EX( K )
            IWORK( 22+K, 2 ) = EXPOS( K )
   10    CONTINUE
      END IF
      K = 22 + NEXTRA
*
*     Get the smallest error detected anywhere (BIGNUM if no error)
*
      CALL GLOBCHK( DESCA( CTXT_ ), K, IWORK, LDW, IWORK2, INFO )
*
*     Prepare output: set info = 0 if no error, and divide by DESCMULT
*     if error is not in a descriptor entry.
*
      IF( INFO.EQ.BIGNUM ) THEN
         INFO = 0
      ELSE IF( MOD( INFO, DESCMULT ) .EQ. 0 ) THEN
         INFO = -INFO / DESCMULT
      ELSE
         INFO = -INFO
      END IF
*
      RETURN
*
*     End of PCHK2MAT
*
      END
*
      SUBROUTINE GLOBCHK( ICTXT, N, X, LDX, IWORK, INFO )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, INFO, LDX, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( N ), X( LDX, 2 )
*     ..
*
*  Purpose
*  =======
*
*  GLOBCHK checks that values in X(i,1) are the same on all processes
*  in the process grid indicated by ICTXT.
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle indicating the context over which
*          the values are to be the same.
*
*  N       (global input) INTEGER
*          The number of values to be compared.
*
*  X       (local input) INTEGER array, dimension (N,2)
*          The 1st column contains the values which should be the same
*          on all processes.  The 2nd column indicates where in the
*          calling routine's parameter list the corresponding value
*          from column 1 came from.
*
*  LDX     (local input) INTEGER
*          The leading dimension of the array X. LDX >= MAX(1,N).
*
*  IWORK   (local workspace) INTEGER array, dimension (N)
*          Used to receive other processes' values for comparing with X.
*
*  INFO    (local input/global output) INTEGER
*          On entry, the smallest error flag so far generated, or BIGNUM
*          for no error. On exit:
*          = BIGNUM : no error
*          < 0: if INFO = -i*100, the i-th argument had an illegal
*               value, or was different between processes.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            K, MYROW, MYCOL
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGAMN2D, IGEBR2D, IGEBS2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, IWORK, K, MYROW, MYCOL )
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         CALL IGEBS2D( ICTXT, 'All', ' ', N, 1, X, N )
      ELSE
         CALL IGEBR2D( ICTXT, 'All', ' ', N, 1, IWORK, N, 0, 0 )
         DO 10 K = 1, N
            IF( X( K, 1 ).NE.IWORK( K ) )
     $         INFO = MIN( INFO, X( K, 2 ) )
   10    CONTINUE
      END IF
*
      CALL IGAMN2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, K, K, -1, -1, 0 )
*
      RETURN
*
*     End GLOBCHK
*
      END
