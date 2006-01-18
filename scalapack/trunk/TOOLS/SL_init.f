      SUBROUTINE SL_INIT( ICTXT, NPROW, NPCOL )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, NPCOL, NPROW
*     ..
*
*  Purpose
*  =======
*
*  SL_INIT initializes an NPROW x NPCOL process grid using a row-major
*  ordering  of  the  processes. This routine retrieves a default system
*  context  which  will  include all available processes. In addition it
*  spawns the processes if needed.
*
*  Arguments
*  =========
*
*  ICTXT   (global output) INTEGER
*          ICTXT specifies the BLACS context handle identifying the
*          created process grid.  The context itself is global.
*
*  NPROW   (global input) INTEGER
*          NPROW specifies the number of process rows in the grid
*          to be created.
*
*  NPCOL   (global input) INTEGER
*          NPCOL specifies the number of process columns in the grid
*          to be created.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            IAM, NPROCS
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GET, BLACS_GRIDINIT, BLACS_PINFO,
     $                   BLACS_SETUP
*     ..
*     .. Executable Statements ..
*
*     Get starting information
*
      CALL BLACS_PINFO( IAM, NPROCS )
*
*     If machine needs additional set up, do it now
*
      IF( NPROCS.LT.1 ) THEN
         IF( IAM.EQ.0 )
     $      NPROCS = NPROW * NPCOL
         CALL BLACS_SETUP( IAM, NPROCS )
      END IF
*
*     Define process grid
*
      CALL BLACS_GET( -1, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
*
      RETURN
*
*     End of SL_INIT
*
      END
