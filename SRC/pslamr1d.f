      SUBROUTINE PSLAMR1D( N, A, IA, JA, DESCA, B, IB, JB, DESCB )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     October 15, 1999
*
*     .. Scalar Arguments ..
      INTEGER            IA, IB, JA, JB, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      REAL               A( * ), B( * )
*     ..
*
*  Bugs
*  ====
*
*  I am not sure that this works correctly when IB and JB are not equal
*  to 1.  Indeed, I suspect that IB should always be set to 1 or ignored
*  with 1 used in its place.
*
*  PSLAMR1D has not been tested except withint the contect of
*  PSSYPTRD, the prototype reduction to tridiagonal form code.
*
*  Purpose
*
*  =======
*
*  PSLAMR1D redistributes a one-dimensional row vector from one data
*  decomposition to another.
*
*  This is an auxiliary routine called by PSSYTRD to redistribute D, E
*  and TAU.
*
*  Notes
*  =====
*
*  Although all processes call PSGEMR2D, only the processes that own
*  the first column of A send data and only processes that own the
*  first column of B receive data.  The calls to SGEBS2D/SGEBR2D
*  spread the data down.
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The size of the matrix to be transposed.
*
*  A       (local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LOCc(JA+N-1)).
*          On output, A is replicated across all processes in
*          this processor column.
*
*  IA      (global input) INTEGER
*          A's global row index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  B       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LOCc(JB+N-1)).
*
*  IB      (global input) INTEGER
*          B's global row index,  NOT USED
*
*  JB      (global input) INTEGER
*          B's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
*
*  WORK    (local workspace) COMPLEX*16 array, dimension ( LWORK )
*
*  LWORK   (local input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= NB * NUMROC( N, 1, 0, 0, NPROW )
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICTXT, MYCOL, MYROW, NPCOL, NPROW, NQ
*     ..
*     .. Local Arrays ..
      INTEGER            DESCAA( DLEN_ ), DESCBB( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PSGEMR2D, SGEBR2D, SGEBS2D
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      DO 10 I = 1, DLEN_
         DESCAA( I ) = DESCA( I )
         DESCBB( I ) = DESCB( I )
   10 CONTINUE
*
      DESCAA( M_ ) = 1
      DESCBB( M_ ) = 1
      DESCAA( LLD_ ) = 1
      DESCBB( LLD_ ) = 1
*
      ICTXT = DESCB( CTXT_ )
      CALL PSGEMR2D( 1, N, A, IA, JA, DESCAA, B, IB, JB, DESCBB, ICTXT )
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NQ = NUMROC( N, DESCB( NB_ ), MYCOL, 0, NPCOL )
*
      IF( MYROW.EQ.0 ) THEN
         CALL SGEBS2D( ICTXT, 'C', ' ', NQ, 1, B, NQ )
      ELSE
         CALL SGEBR2D( ICTXT, 'C', ' ', NQ, 1, B, NQ, 0, MYCOL )
      END IF
*
      RETURN
*
*     End of PSLAMR1D
*
      END
