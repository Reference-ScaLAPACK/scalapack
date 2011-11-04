      SUBROUTINE PZROW2COL( ICTXT, M, N, NB, VS, LDVS, VD, LDVD,
     $                     RSRC, CSRC, RDEST, CDEST, WORK)
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
*  a row of processes, and distribute those rows over a column of
*  processes. This routine minimizes communication by sending all
*  information it has that a given process in the CDEST needs at once.
*  To do this it uses the least common multiple (LCM) concept.  This is
*  simply the realization that if I have part of a vector split over a
*  process row consisting of Q processes, and I want to send all of that
*  vector that I own to a new vector distributed over P processes within
*  a process column, that after I find the process in RDEST that owns
*  the row of the vector I'm currently looking at, he will want every
*  ( (LCM(P,Q)/Q ) block of my vector (the block being of size NB x N).
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
*          The number of vectors in the vector block
*
*  NB      (global input) INTEGER
*          The blocking factor used to divide the rows of the vector
*          amongst the processes of a row.
*
*  VS      (local input) COMPLEX*16
*          Array of dimension (LDVS,N), the block of vectors stored on
*          process row RSRC to be put into memory VD, and stored on
*          process column CDEST.
*
*  LDVS    (local input) INTEGER
*          The leading dimension of VS.
*
*  VD      (local output) COMPLEX*16
*          Array of dimension (LDVD,N), on output, the contents of VD
*          stored on process column CDEST will be here.
*
*  LDVD    (local input) INTEGER
*          The leading dimension of VD.
*
*  RSRC    (global input) INTEGER
*          The process row VS is distributed over.
*
*  CSRC    (global input) INTEGER
*          The process column the distributed block of vectors VS
*          begins on.
*
*  RDEST   (global input) INTEGER
*          The process row that VD begins on.
*
*  CDEST   (global input) INTEGER
*          The process column to distribute VD over.
*
*  WORK    (local workspace) COMPLEX*16
*          Array, dimension (LDW). The required size of work varies:
*          if( nprow.eq.npcol ) then
*             LDW = 0; WORK not accessed.
*          else
*            lcm = least common multiple of process rows and columns.
*             Mq  = number of rows of VS on my process.
*             npcol = number of process columns
*             CEIL = the ceiling of given operation
*             LDW = NB*N*CEIL( CEIL( Mq/NB )/(LCM/npcol) )
*          end if
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            CBLKSKIP, ICPY, II, ISTART, ICSRC, IRDEST, JB,
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
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     If we are not in special case for NPROW = NPCOL where there is no
*     copying required
*
      IF( NPROW .NE. NPCOL ) THEN
         LCM = ILCM( NPROW, NPCOL )
         RBLKSKIP = LCM / NPCOL
         CBLKSKIP = LCM / NPROW
*
*        If I have part of VS, the source vector(s)
*
         IF( MYROW.EQ.RSRC ) THEN
*
            ISTART = 1
*
*           Figure my distance from CSRC: the process in CDEST the same
*           distance from RDEST will want my first block
*
            MYDIST = MOD( NPCOL+MYCOL-CSRC, NPCOL )
            MQ = NUMROC( M, NB, MYCOL, CSRC, NPCOL )
            IRDEST = MOD( RDEST+MYDIST, NPROW )
*
*           Loop over all possible destination processes
*
            DO 20 K = 1, RBLKSKIP
               JJ = 1
*
*              If I am not destination process
*
               IF( (MYROW.NE.IRDEST).OR.(MYCOL.NE.CDEST) ) THEN
*
*                 Pack all data I own that destination needs
*
                  DO 10 II = ISTART, MQ, NB*RBLKSKIP
                     JB = MIN( NB, MQ-II+1 )
                     CALL ZLACPY( 'G', JB, N, VS(II,1), LDVS,
     $                            WORK(JJ), JB )
                     JJ = JJ + NB*N
   10             CONTINUE
*
*                 Figure how many rows are to be sent and send them if
*                 necessary, NOTE: will send extra if NB > JB
*
                  JJ = JJ - 1
                  IF( JJ.GT.0 )
     $               CALL ZGESD2D( ICTXT, JJ, 1, WORK, JJ, IRDEST,
     $                             CDEST )
*
*              I am both source and destination, save where to start
*              copying from for later use
*
               ELSE
                  ICPY = ISTART
               END IF
*
               ISTART = ISTART + NB
               IRDEST = MOD( IRDEST+NPCOL, NPROW )
   20       CONTINUE
         END IF
*
*        If I should receive info into VD
*
         IF( MYCOL.EQ.CDEST ) THEN
*
            ISTART = 1
*
*           Figure my distance from CDEST: the process in CSRC the same
*           distance from RSRC will have my first block
*
            MYDIST = MOD( NPROW+MYROW-RDEST, NPROW )
            MP = NUMROC( M, NB, MYROW, RDEST, NPROW )
            ICSRC = MOD( CSRC+MYDIST, NPCOL )
*
*           Loop over all sending processes
*
            DO 50 K = 1, CBLKSKIP
*
*              If I don't already possess the required data
*
               IF( (MYROW.NE.RSRC).OR.(MYCOL.NE.ICSRC) ) THEN
*
*                 Figure how many rows to receive, and receive them
*                 NOTE: may receive to much -- NB instead of JB
*
                  NBLOCKS = (MP - ISTART + NB) / NB
                  JJ = ((NBLOCKS+CBLKSKIP-1) / CBLKSKIP)*NB
                  IF( JJ.GT.0 )
     $               CALL ZGERV2D( ICTXT, JJ, N, WORK, JJ, RSRC, ICSRC )
*
*                 Copy data to destination vector
*
                  JJ = 1
                  DO 30 II = ISTART, MP, NB*CBLKSKIP
                     JB = MIN( NB, MP-II+1 )
                     CALL ZLACPY( 'G', JB, N, WORK(JJ), JB, VD(II,1),
     $                            LDVD )
                     JJ = JJ + NB*N
   30             CONTINUE
*
*              If I am both source and destination
*
               ELSE
                  JJ = ICPY
                  DO 40 II = ISTART, MP, NB*CBLKSKIP
                     JB = MIN( NB, MP-II+1 )
                     CALL ZLACPY( 'G', JB, N, VS(JJ,1), LDVS, VD(II,1),
     $                            LDVD )
                     JJ = JJ + NB*RBLKSKIP
   40             CONTINUE
               END IF
               ISTART = ISTART + NB
               ICSRC = MOD( ICSRC+NPROW, NPCOL )
   50       CONTINUE
         END IF
*
*     if NPROW = NPCOL, there is a one-to-one correspondance between
*     process rows and columns, so no work space or copying required
*
      ELSE
*
         IF( MYROW.EQ.RSRC ) THEN
*
*           Figure my distance from CSRC: the process in CDEST the same
*           distance from RDEST will want my piece of the vector
*
            MYDIST = MOD( NPCOL+MYCOL-CSRC, NPCOL )
            MQ = NUMROC( M, NB, MYCOL, CSRC, NPCOL )
            IRDEST = MOD( RDEST+MYDIST, NPROW )
            IF( (MYROW.NE.IRDEST).OR.(MYCOL.NE.CDEST) ) THEN
               CALL ZGESD2D( ICTXT, MQ, N, VS, LDVS, IRDEST, CDEST )
            ELSE
               CALL ZLACPY( 'G', MQ, N, VS, LDVS, VD, LDVD )
            END IF
         END IF
         IF( MYCOL.EQ.CDEST ) THEN
*
*           Figure my distance from RDEST: the process in RSRC the same
*           distance from CSRC will have my piece of the vector
*
            MYDIST = MOD( NPROW+MYROW-RDEST, NPROW )
            MP = NUMROC( M, NB, MYROW, RDEST, NPROW )
            ICSRC = MOD( CSRC+MYDIST, NPCOL )
            IF( (MYCOL.NE.ICSRC).OR.(MYROW.NE. RSRC) )
     $         CALL ZGERV2D( ICTXT, MP, N, VD, LDVD, RSRC, ICSRC )
         END IF
      END IF
*
      RETURN
*
*     End of PZROW2COL
*
      END
