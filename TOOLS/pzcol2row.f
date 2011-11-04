      SUBROUTINE PZCOL2ROW( ICTXT, M, N, NB, VS, LDVS, VD, LDVD, RSRC,
     $                      CSRC, RDEST, CDEST, WORK)
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            CDEST, CSRC, ICTXT, LDVD, LDVS, M, N, NB,
     $                   RDEST, RSRC
*     ..
*     .. Array Arguments ..
      COMPLEX*16         VD( LDVD, * ), VS( LDVS, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  Take a block of vectors with M total rows which are distributed over
*  a column of processes, and distribute those rows over a row of
*  processes. This routine minimizes communication by sending all
*  information it has that a given process in the RDEST needs at once.
*  To do this it uses the least common multiple (LCM) concept.  This is
*  simply the realization that if I have part of a vector split over a
*  process column consisting of P processes, and I want to send all of
*  that vector that I own to a new vector distributed over Q processes
*  within a process row, that after I find the process in RDEST that
*  owns the row of the vector I'm currently looking at, he will want
*  every ( (LCM(P,Q) / P ) block of my vector (the block being of size
*  NB x N).
*
*  Arguments
*  =========
*
*  Rem:  MP, resp. NQ, denotes the number of local rows, resp. local
*  ====  columns, necessary to store a global vector of dimension M
*        across P processes, resp. N over Q processes.
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation. The context itself is global.
*
*  M       (global input) INTEGER
*          The number of global rows each vector has.
*
*  N       (global input) INTEGER
*          The number of vectors in the vector block.
*
*  NB      (global input) INTEGER
*          The blocking factor used to divide the rows of the vector
*          amongst the processes of a column.
*
*  VS      (local input) COMPLEX*16
*          Array of dimension (LDVS,N), the block of vectors stored on
*          process column CSRC to be put into memory VD, and stored
*          on process row RDEST.
*
*  LDVS    (local input) INTEGER
*          The leading dimension of VS, LDVS >= MAX( 1, MP ).
*
*  VD      (local output) COMPLEX*16
*          Array of dimension (LDVD,N), on output, the contents of VS
*          stored on process row RDEST will be here.
*
*  LDVD    (local input) INTEGER
*          The leading dimension of VD, LDVD >= MAX( 1, MQ ).
*
*  RSRC    (global input) INTEGER
*          The process row the distributed block of vectors VS begins
*          on.
*
*  CSRC    (global input) INTEGER
*          The process column VS is distributed over.
*
*  RDEST   (global input) INTEGER
*          The process row to distribute VD over.
*
*  CDEST   (global input) INTEGER
*          The process column that VD begins on.
*
*  WORK    (local workspace) COMPLEX*16
*          Array of dimension (LDW), the required size of work varies:
*          if( nprow.eq.npcol ) then
*             LDW = 0; WORK not accessed.
*          else
*             lcm = least common multiple of process rows and columns.
*             Mp  = number of rows of VS on my process.
*             nprow = number of process rows
*             CEIL = the ceiling of given operation
*             LDW = NB*N*CEIL( CEIL( Mp/NB )/(LCM/nprow) )
*          end if
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            CBLKSKIP, ICPY, ICDEST, II, IRSRC, ISTART, JB,
     $                   JJ, K, LCM, MP, MQ, MYCOL, MYDIST, MYROW,
     $                   NBLOCKS, NPCOL, NPROW, RBLKSKIP
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, ZGESD2D, ZGERV2D, ZLACPY
*     ..
*     .. External Functions ..
      INTEGER            ILCM, NUMROC
      EXTERNAL           ILCM, NUMROC
*     ..
*     .. Executable Statements ..
*
	ICPY = 0
*	
*     Get grid parameters.
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     If we are not in special case for NPROW = NPCOL where there
*     is no copying required
*
      IF( NPROW.NE.NPCOL ) THEN
         LCM = ILCM( NPROW, NPCOL )
         RBLKSKIP = LCM / NPCOL
         CBLKSKIP = LCM / NPROW
*
*        If I have part of VS, the source vector(s)
*
         IF( MYCOL.EQ.CSRC ) THEN
*
            ISTART = 1
*
*           Figure my distance from RSRC: the process in RDEST the same
*           distance from CDEST will want my first block
*
            MYDIST = MOD( NPROW+MYROW-RSRC, NPROW )
            MP = NUMROC( M, NB, MYROW, RSRC, NPROW )
            ICDEST = MOD( CDEST+MYDIST, NPCOL )
*
*           Loop over all possible destination processes
*
            DO 20 K = 1, CBLKSKIP
               JJ = 1
*
*              If I am not destination process
*
               IF( (MYCOL.NE.ICDEST).OR.(MYROW.NE.RDEST) ) THEN
*
*                 Pack all data I own that destination needs
*
                  DO 10 II = ISTART, MP, NB*CBLKSKIP
                     JB = MIN(NB, MP-II+1)
                     CALL ZLACPY( 'G', JB, N, VS(II,1), LDVS,
     $                            WORK(JJ), JB )
                     JJ = JJ + NB*N
   10             CONTINUE
*
*                 Figure how many rows are to be sent and send them if
*                 necessary (NOTE: will send extra if NB > JB)
*
                  JJ = JJ - 1
                  IF( JJ.GT.0 )
     $               CALL ZGESD2D( ICTXT, JJ, 1, WORK, JJ, RDEST,
     $                             ICDEST )
*
               ELSE
*
*                 I am both source and destination, save where to start
*                 copying from for later use.
*
                  ICPY = ISTART
               END IF
*
               ISTART = ISTART + NB
               ICDEST = MOD(ICDEST+NPROW, NPCOL)
   20       CONTINUE
         END IF
*
*        If I should receive info into VD
*
         IF( MYROW.EQ.RDEST ) THEN
*
            ISTART = 1
*
*           Figure my distance from CDEST: the process in CSRC the same
*           distance from RSRC will have my first block.
*
            MYDIST = MOD( NPCOL+MYCOL-CDEST, NPCOL )
            MQ = NUMROC( M, NB, MYCOL, CDEST, NPCOL )
            IRSRC = MOD( RSRC+MYDIST, NPROW )
            DO 50 K = 1, RBLKSKIP
*
*              If I don't already possess the required data
*
               IF( (MYCOL.NE.CSRC).OR.(MYROW.NE.IRSRC) ) THEN
*
*                 Figure how many rows to receive, and receive them
*                 NOTE: may receive to much -- NB instead of JB
*
                  NBLOCKS = (MQ - ISTART + NB) / NB
                  JJ = ((NBLOCKS+RBLKSKIP-1) / RBLKSKIP)*NB
                  IF( JJ.GT.0 )
     $               CALL ZGERV2D( ICTXT, JJ, N, WORK, JJ, IRSRC, CSRC )
*
*                 Copy data to destination vector
*
                  JJ = 1
                  DO 30 II = ISTART, MQ, NB*RBLKSKIP
                     JB = MIN( NB, MQ-II+1 )
                     CALL ZLACPY( 'G', JB, N, WORK(JJ), JB,
     $                            VD(II,1), LDVD )
                     JJ = JJ + NB*N
   30             CONTINUE
*
*                 If I am both source and destination
*
               ELSE
                  JJ = ICPY
                  DO 40 II = ISTART, MQ, NB*RBLKSKIP
                     JB = MIN( NB, MQ-II+1 )
                     CALL ZLACPY( 'G', JB, N, VS(JJ,1), LDVS,
     $                            VD(II,1), LDVD )
                     JJ = JJ + NB*CBLKSKIP
   40             CONTINUE
               END IF
               ISTART = ISTART + NB
               IRSRC = MOD( IRSRC+NPCOL, NPROW )
   50       CONTINUE
         END IF
*
*     If NPROW = NPCOL, there is a one-to-one correspondance between
*     process rows and columns, so no work space or copying required
*
      ELSE
*
         IF( MYCOL.EQ.CSRC ) THEN
*
*           Figure my distance from RSRC: the process in RDEST the same
*           distance from CDEST will want my piece of the vector.
*
            MYDIST = MOD( NPROW+MYROW-RSRC, NPROW )
            MP = NUMROC( M, NB, MYROW, RSRC, NPROW )
            ICDEST = MOD( CDEST+MYDIST, NPCOL )
*
            IF( (MYCOL.NE.ICDEST).OR.(MYROW.NE.RDEST) ) THEN
               CALL ZGESD2D( ICTXT, MP, N, VS, LDVS, RDEST, ICDEST )
            ELSE
               CALL ZLACPY( 'G', MP, N, VS, LDVS, VD, LDVD )
            END IF
         END IF
*
         IF( MYROW.EQ.RDEST ) THEN
*
*           Figure my distance from CDEST: the process in CSRC the same
*           distance from RSRC will have my piece of the vector.
*
            MYDIST = MOD( NPCOL+MYCOL-CDEST, NPCOL )
            MQ = NUMROC( M, NB, MYCOL, CDEST, NPCOL )
            IRSRC = MOD( RSRC+MYDIST, NPROW )
*
            IF( (MYROW.NE.IRSRC).OR.(MYCOL.NE.CSRC) )
     $         CALL ZGERV2D( ICTXT, MQ, N, VD, LDVD, IRSRC, CSRC )
*
          END IF
*
      END IF
*
      RETURN
*
*     End of PZCOL2ROW
*
      END
