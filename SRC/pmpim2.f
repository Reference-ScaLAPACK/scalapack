***********************************************************************
*
*     Below are auxiliary subroutines for eigenpair assignments
*
***********************************************************************

      SUBROUTINE PMPIM2( IL, IU, NPROCS, PMYILS, PMYIUS )

      IMPLICIT NONE

      INTEGER IL, IU, M, NPROCS, PRCCTR
      INTEGER PMYILS( * ), PMYIUS( * )
*
*     Purpose
*     =======
*
*     PMPIM2 is the scheduling subroutine.
*     It computes for all processors the eigenpair range assignments.
*
*     Arguments
*     =========
*
*     IL, IU  (input) INTEGER
*             The range of eigenpairs to be computed
*
*     NPROCS  (input) INTEGER
*             The total number of processors available
*
*     PMYILS  (output) INTEGER array
*             For each processor p,  PMYILS(p) is the index 
*             of the first eigenvalue in W to be computed
*             PMYILS(p) equals zero if p stays idle
*
*     PMYIUS  (output) INTEGER array
*             For each processor p,  PMYIUS(p) is the index
*             of the last eigenvalue in W to be computed
*             PMYIUS(p) equals zero if p stays idle
*

      M = IU - IL + 1

      IF ( NPROCS.GT.M ) THEN
         DO 10 PRCCTR = 0, NPROCS-1
            IF ( PRCCTR.LT.M ) THEN
               PMYILS(PRCCTR+1) = PRCCTR + IL
               PMYIUS(PRCCTR+1) = PRCCTR + IL
            ELSE
               PMYILS(PRCCTR+1) = 0
               PMYIUS(PRCCTR+1) = 0
            END IF
 10      CONTINUE
      ELSE
         DO 20 PRCCTR = 0, NPROCS-1
            PMYILS(PRCCTR+1) = (PRCCTR * (M / NPROCS)) + IL
            IF (PRCCTR.LT.MOD(M, NPROCS)) THEN
               PMYILS(PRCCTR+1) = PMYILS(PRCCTR+1) + PRCCTR
               PMYIUS(PRCCTR+1) = PMYILS(PRCCTR+1) + M / NPROCS
            ELSE
               PMYILS(PRCCTR+1) = PMYILS(PRCCTR+1) + MOD(M, NPROCS)
               PMYIUS(PRCCTR+1) = PMYILS(PRCCTR+1) + M / NPROCS - 1
            END IF
 20      CONTINUE
      END IF

      RETURN
      END


