***********************************************************************
*
*     Auxiliary subroutine for eigenpair assignments
*
***********************************************************************
      SUBROUTINE PMPCOL( MYPROC, NPROCS, IIL, NEEDIL, NEEDIU, 
     $                   PMYILS, PMYIUS,
     $                   COLBRT, FRSTCL, LASTCL )

      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      INTEGER            FRSTCL, IIL, LASTCL, MYPROC, NEEDIL, NEEDIU,
     $                   NPROCS
      LOGICAL COLBRT
*     ..
*     .. Array Arguments ..
      INTEGER            PMYILS( * ), PMYIUS( * )
*     ..
*
*  Purpose
*  =======
*
*  Using the output from PMPIM2 and given the information on
*  eigenvalue clusters, PMPCOL finds the collaborators of MYPROC.
*
*  Arguments
*  =========
*
*  MYPROC  (input) INTEGER
*          The processor number, 0 <= MYPROC < NPROCS
*
*  NPROCS  (input) INTEGER
*          The total number of processors available
*
*  IIL     (input) INTEGER
*          The index of the leftmost eigenvalue in W
*
*  NEEDIL  (input) INTEGER
*          The leftmost position in W needed by MYPROC
*
*  NEEDIU  (input) INTEGER
*          The rightmost position in W needed by MYPROC
*
*  PMYILS  (input) INTEGER array
*          For each processor p,  PMYILS(p) is the index
*          of the first eigenvalue in W to be computed
*          PMYILS(p) equals zero if p stays idle
*
*  PMYIUS  (input) INTEGER array
*          For each processor p,  PMYIUS(p) is the index
*          of the last eigenvalue in W to be computed
*          PMYIUS(p) equals zero if p stays idle
*
*  COLBRT  (output) LOGICAL
*          TRUE if MYPROC collaborates.
*
*  FRSTCL  (output) INTEGER
*  LASTCL  FIRST and LAST collaborator of MYPROC   
*          MYPROC collaborates with
*          FRSTCL, ..., MYPROC-1, MYPROC+1, ...,LASTCL 
*          If MYPROC == FRSTCL, there are no collaborators 
*          on the left. IF MYPROC == LASTCL, there are no
*          collaborators on the right.
*          If FRSTCL == 0 and LASTCL = NPROCS-1, then
*          MYPROC collaborates with everybody
*

*     .. Local Scalars ..
      INTEGER I, NEEDIIL, NEEDIIU
*     ..
*     .. Executable Statements ..
*     Compute global eigenvalue index from position in W
      NEEDIIL = NEEDIL + IIL - 1
      NEEDIIU = NEEDIU + IIL - 1

*     Find processor responsible for NEEDIL, this is the first
*     collaborator
      DO 1 I = 1, NPROCS
         IF( PMYILS(I).GT.NEEDIIL) GOTO 2
         FRSTCL = I-1
 1    CONTINUE
 2    CONTINUE

*     Find processor responsible for NEEDIU, this is the last
*     collaborator
      DO 3 I = NPROCS,1,-1
         IF( PMYIUS(I).LT.NEEDIIU ) THEN  
*          Need to check special case: does this proc work at all?
           IF( PMYIUS(I).GT.0 )
     $        GOTO 4
         ENDIF
         LASTCL = I-1
 3    CONTINUE
 4    CONTINUE

*     Decide if there is a collaboration
      IF( (FRSTCL.LT.MYPROC).OR.(LASTCL.GT.MYPROC) ) THEN
         COLBRT = .TRUE.
      ELSE
         COLBRT = .FALSE.
      ENDIF

      RETURN
      END
