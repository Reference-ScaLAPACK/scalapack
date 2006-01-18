      SUBROUTINE PBCTRNV( ICONTXT, XDIST, TRANS, N, NB, NZ, X, INCX,
     $                    BETA, Y, INCY, IXROW, IXCOL, IYROW, IYCOL,
     $                    WORK )
*
*  -- PB-BLAS routine (version 2.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory.
*     April 28, 1996
*
*     Jaeyoung Choi, Oak Ridge National Laboratory
*     Jack Dongarra, University of Tennessee and Oak Ridge National Lab.
*     David Walker,  Oak Ridge National Laboratory
*
*     .. Scalar Arguments ..
      CHARACTER*1        TRANS, XDIST
      INTEGER            ICONTXT, INCX, INCY, IXCOL, IXROW, IYCOL,
     $                   IYROW, N, NB, NZ
      COMPLEX            BETA
*     ..
*     .. Array Arguments ..
      COMPLEX            WORK( * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  PBCTRNV transposes a column vector to row vector, or a row vector to
*  column vector by reallocating data distribution.
*
*     Y := X'
*
*  where X and Y are N vectors.
*
*  Parameters
*  ==========
*
*  ICONTXT (input) INTEGER
*          ICONTXT is the BLACS mechanism for partitioning communication
*          space.  A defining property of a context is that a message in
*          a context cannot be sent or received in another context.  The
*          BLACS context includes the definition of a grid, and each
*          process' coordinates in it.
*
*  XDIST   (input) CHARACTER*1
*          XDIST specifies whether X is a column vector or a row vector,
*
*            XDIST = 'C',  X is a column vector (distributed columnwise)
*            XDIST = 'R',  X is a row vector    (distributed rowwise)
*
*  TRANS   (input) CHARACTER*1
*          TRANS specifies whether the transposed format is transpose
*          or conjugate transpose.  If the vectors X and Y are real,
*          the argument is ignored.
*
*             TRANS = 'T',  transpose
*             TRANS = 'C',  conjugate transpose
*
*  N       (input) INTEGER
*          N specifies the (global) number of the vector X and the
*          vector Y.  N >= 0.
*
*  NB      (input) INTEGER
*          NB specifies the block size of vectors X and Y.  NB >= 0.
*
*  NZ      (input) INTEGER
*          NZ is the column offset to specify the column distance from
*          the beginning of the block to the first element of the
*          vector X, and the row offset to the first element of the
*          vector Y if XDIST = 'C'.
*          Otherwise, it is row offset to specify the row distance
*          from the beginning of the block to the first element of the
*          vector X, and the column offset to the first element of the
*          vector Y.  0 < NZ <= NB.
*
*  X       (input) COMPLEX array of dimension at least
*          ( 1 + (Np-1) * abs(INCX)) in IXCOL if XDIST = 'C', or
*          ( 1 + (Nq-1) * abs(INCX)) in IXROW if XDIST = 'R'.
*          The incremented array X must contain the vector X.
*
*  INCX    (input) INTEGER
*          INCX specifies the increment for the elements of X.
*          INCX <> 0.
*
*  BETA    (input) COMPLEX
*          BETA specifies scaler beta.
*
*  Y       (input/output) COMPLEX array of dimension at least
*          ( 1 + (Nq-1) * abs(INCY)) in IYROW if XDIST = 'C', or
*          ( 1 + (Np-1) * abs(INCY)) in IYCOL if XDIST = 'R', or
*          The incremented array Y must contain the vector Y.
*          Y will not be referenced if beta is zero.
*
*  INCY    (input) INTEGER
*          INCY specifies the increment for the elements of Y.
*          INCY <> 0.
*
*  IXROW   (input) INTEGER
*          IXROW specifies a row of the process template, which holds
*          the first element of the vector X. If X is a row vector and
*          all rows of processes have a copy of X, then set IXROW = -1.
*
*  IXCOL   (input) INTEGER
*          IXCOL specifies  a column of the process template,
*          which holds the first element of the vector X.  If  X is  a
*          column block and all columns of processes have a copy of X,
*          then set IXCOL = -1.
*
*  IYROW   (input) INTEGER
*          IYROW specifies the current row process which holds the
*          first element of the vector Y, which is transposed of X.
*          If X  is a column vector and the transposed  row vector Y is
*          distributed all rows of processes, set IYROW = -1.
*
*  IYCOL   (input) INTEGER
*          IYCOL specifies  the current column process  which holds
*          the first element of the vector Y, which is transposed of Y.
*          If X is a row block and the transposed column vector Y is
*          distributed all columns of processes, set IYCOL = -1.
*
*  WORK    (workspace) COMPLEX array of dimension Size(WORK).
*          It needs extra working space of x**T or x**H.
*
*  Parameters Details
*  ==================
*
*  Nx      It is a local portion  of N owned by a process, where x is
*          replaced by  either p (=NPROW) or q (=NPCOL)).  The value is
*          determined by N, NB, NZ, x, and MI, where NB is a block size,
*          NZ is a offset from the beginning of the block,  and MI is a
*          row or column position  in a process template. Nx is equal
*          to  or less than Nx0 = CEIL( N+NZ, NB*x ) * NB.
*
*  Communication Scheme
*  ====================
*
*  The communication scheme of the routine is set to '1-tree', which is
*  fan-out.  (For details, see BLACS user's guide.)
*
*  Memory Requirement of WORK
*  ==========================
*
*  NN   = N + NZ
*  Npb  = CEIL( NN, NB*NPROW )
*  Nqb  = CEIL( NN, NB*NPCOL )
*  LCMP = LCM / NPROW
*  LCMQ = LCM / NPCOL
*
*   (1) XDIST = 'C'
*     (a) IXCOL != -1
*         Size(WORK) = CEIL(Nqb,LCMQ)*NB
*     (b) IXCOL = -1
*         Size(WORK) = CEIL(Nqb,LCMQ)*NB * MIN(LCMQ,CEIL(NN,NB))
*
*   (2) XDIST = 'R'
*     (a) IXROW != -1
*         Size(WORK) = CEIL(Npb,LCMP)*NB
*     (b) IXROW = -1
*         Size(WORK) = CEIL(Npb,LCMP)*NB * MIN(LCMP,CEIL(NN,NB))
*
*  Notes
*  -----
*  More precise space can be computed as
*
*  CEIL(Npb,LCMP)*NB => NUMROC( NUMROC(NN,NB,0,0,NPROW), NB, 0, 0, LCMP)
*  CEIL(Nqb,LCMQ)*NB => NUMROC( NUMROC(NN,NB,0,0,NPCOL), NB, 0, 0, LCMQ)
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE  = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLFORM, ROWFORM
      INTEGER            I, IDEX, IGD, INFO, JDEX, JYCOL, JYROW, JZ, KZ,
     $                   LCM, LCMP, LCMQ, MCCOL, MCROW, MRCOL, MRROW,
     $                   MYCOL, MYROW, NN, NP, NP0, NP1, NPCOL, NPROW,
     $                   NQ, NQ0, NQ1
      COMPLEX            TBETA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILCM, ICEIL, NUMROC
      EXTERNAL           LSAME, ILCM, ICEIL, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CGEBR2D, CGEBS2D, CGERV2D,
     $                   CGESD2D, PBCTR2A1, PBCTR2B1, PBCTRGET,
     $                   PBCTRST1, PBCVECADD, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      IF( N.EQ.0 ) RETURN
*
      CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      COLFORM = LSAME( XDIST, 'C' )
      ROWFORM = LSAME( XDIST, 'R' )
*
*     Test the input parameters.
*
      INFO = 0
      IF( ( .NOT.COLFORM ) .AND. ( .NOT.ROWFORM ) ) THEN
         INFO = 2
      ELSE IF( N   .LT.0                          ) THEN
         INFO = 4
      ELSE IF( NB  .LT.1                          ) THEN
         INFO = 5
      ELSE IF( NZ  .LT.0 .OR. NZ.GE.NB            ) THEN
         INFO = 6
      ELSE IF( INCX.EQ.0                          ) THEN
         INFO = 8
      ELSE IF( INCY.EQ.0                          ) THEN
         INFO = 11
      ELSE IF( IXROW.LT.-1 .OR. IXROW.GE.NPROW .OR.
     $       ( IXROW.EQ.-1 .AND. COLFORM )        ) THEN
         INFO = 12
      ELSE IF( IXCOL.LT.-1 .OR. IXCOL.GE.NPCOL .OR.
     $       ( IXCOL.EQ.-1 .AND. ROWFORM )        ) THEN
         INFO = 13
      ELSE IF( IYROW.LT.-1 .OR. IYROW.GE.NPROW .OR.
     $       ( IYROW.EQ.-1 .AND. ROWFORM )        ) THEN
         INFO = 14
      ELSE IF( IYCOL.LT.-1 .OR. IYCOL.GE.NPCOL .OR.
     $       ( IYCOL.EQ.-1 .AND. COLFORM )        ) THEN
         INFO = 15
      END IF
*
   10 CONTINUE
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICONTXT, 'PBCTRNV ', INFO )
         RETURN
      END IF
*
*     Start the operations.
*
*     LCM : the least common multiple of NPROW and NPCOL
*
      LCM  = ILCM( NPROW, NPCOL )
      LCMP = LCM   / NPROW
      LCMQ = LCM   / NPCOL
      IGD  = NPCOL / LCMP
      NN   = N + NZ
*
*     When x is a column vector
*
      IF( COLFORM ) THEN
*
*       Form  y <== x'  ( x is a column vector )
*
*                                        ||
*                                        ||
*            _____________               ||
*            -----(y)-----      <==     (x)
*                                        ||
*                                        ||
*                                        ||
*
        IF(      IXROW.LT.0  .OR. IXROW.GE.NPROW ) THEN
          INFO = 12
        ELSE IF( IXCOL.LT.-1 .OR. IXCOL.GE.NPCOL ) THEN
          INFO = 13
        ELSE IF( IYROW.LT.-1 .OR. IYROW.GE.NPROW ) THEN
          INFO = 14
        ELSE IF( IYCOL.LT.0  .OR. IYCOL.GE.NPCOL ) THEN
          INFO = 15
        END IF
        IF( INFO.NE.0 ) GO TO 10
*
*       MRROW : row relative position in template from IXROW
*       MRCOL : column relative position in template from IYCOL
*
        MRROW = MOD( NPROW+MYROW-IXROW, NPROW )
        MRCOL = MOD( NPCOL+MYCOL-IYCOL, NPCOL )
        JYROW = IYROW
        IF( IYROW.EQ.-1 ) JYROW = IXROW
*
        NP  = NUMROC( NN, NB, MYROW, IXROW, NPROW )
        IF( MRROW.EQ.0 ) NP = NP - NZ
        NQ  = NUMROC( NN, NB, MYCOL, IYCOL, NPCOL )
        IF( MRCOL.EQ.0 ) NQ = NQ - NZ
        NQ0 = NUMROC( NUMROC(NN, NB, 0, 0, NPCOL), NB, 0, 0, LCMQ )
*
*       When a column process of IXCOL has a column block A,
*
        IF( IXCOL .GE. 0 ) THEN
          TBETA = ZERO
          IF( MYROW.EQ.JYROW ) TBETA = BETA
          KZ = NZ
*
          DO 20 I = 0, MIN( LCM, ICEIL(NN,NB) ) - 1
            MCROW = MOD( MOD(I, NPROW) + IXROW, NPROW )
            MCCOL = MOD( MOD(I, NPCOL) + IYCOL, NPCOL )
            IF( LCMQ.EQ.1 )  NQ0 = NUMROC( NN, NB, I, 0, NPCOL )
            JDEX  = (I/NPCOL) * NB
            IF( MRCOL.EQ.0 ) JDEX = MAX(0, JDEX-NZ)
*
*           A source node copies the blocks to WORK, and send it
*
            IF( MYROW.EQ.MCROW .AND. MYCOL.EQ.IXCOL ) THEN
*
*             The source node is a destination node
*
              IDEX = (I/NPROW) * NB
              IF( MRROW.EQ.0 ) IDEX = MAX( 0, IDEX-NZ )
              IF( MYROW.EQ.JYROW .AND. MYCOL.EQ.MCCOL ) THEN
                CALL PBCTR2B1( ICONTXT, TRANS, NP-IDEX, NB, KZ,
     $                          X(IDEX*INCX+1), INCX, TBETA,
     $                          Y(JDEX*INCY+1), INCY, LCMP, LCMQ )
*
*             The source node sends blocks to a destination node
*
              ELSE
                CALL PBCTR2B1( ICONTXT, TRANS, NP-IDEX, NB, KZ,
     $                         X(IDEX*INCX+1), INCX, ZERO, WORK, 1,
     $                         LCMP, 1 )
                CALL CGESD2D( ICONTXT, 1, NQ0-KZ, WORK, 1,
     $                        JYROW, MCCOL )
              END IF
*
*           A destination node receives the copied vector
*
            ELSE IF( MYROW.EQ.JYROW .AND. MYCOL.EQ.MCCOL ) THEN
              IF( LCMQ.EQ.1 .AND. TBETA.EQ.ZERO ) THEN
                CALL CGERV2D( ICONTXT, 1, NQ0-KZ, Y, INCY,
     $                        MCROW, IXCOL )
              ELSE
                CALL CGERV2D( ICONTXT, 1, NQ0-KZ, WORK, 1,
     $                        MCROW, IXCOL )
                CALL PBCTR2A1( ICONTXT, NQ-JDEX, NB, KZ, WORK, 1, TBETA,
     $                         Y(JDEX*INCY+1), INCY, LCMQ*NB )
              END IF
            END IF
            KZ = 0
   20     CONTINUE
*
*         Broadcast a row block of WORK in each column of template
*
          IF( IYROW.EQ.-1 ) THEN
            IF( MYROW.EQ.JYROW ) THEN
              CALL CGEBS2D( ICONTXT, 'Col', '1-tree', 1, NQ, Y, INCY )
            ELSE
              CALL CGEBR2D( ICONTXT, 'Col', '1-tree', 1, NQ, Y, INCY,
     $                     JYROW, MYCOL )
             END IF
          END IF
*
*       When all column procesors have a copy of the column block A,
*
        ELSE
          IF( LCMQ.EQ.1 ) NQ0 = NQ
*
*         Processors, which have diagonal blocks of X, copy them to
*         WORK array in transposed form
*
          KZ = 0
          IF( MRROW.EQ.0 ) KZ = NZ
          JZ = 0
          IF( MRROW.EQ.0 .AND. MYCOL.EQ.IYCOL ) JZ = NZ
*
          DO 30 I = 0, LCMP - 1
            IF( MRCOL.EQ.MOD(NPROW*I+MRROW, NPCOL) ) THEN
              IDEX = MAX( 0, I*NB-KZ )
              IF( LCMQ.EQ.1 .AND. (IYROW.EQ.-1.OR.IYROW.EQ.MYROW) ) THEN
                 CALL PBCTR2B1( ICONTXT, TRANS, NP-IDEX, NB, JZ,
     $                          X(IDEX*INCX+1), INCX, BETA, Y, INCY,
     $                          LCMP, 1 )
              ELSE
                 CALL PBCTR2B1( ICONTXT, TRANS, NP-IDEX, NB, JZ,
     $                          X(IDEX*INCX+1), INCX, ZERO, WORK, 1,
     $                          LCMP, 1 )
              END IF
            END IF
   30     CONTINUE
*
*         Get diagonal blocks of A for each column of the template
*
          MCROW = MOD( MOD(MRCOL, NPROW) + IXROW, NPROW )
          IF( LCMQ.GT.1 ) THEN
            MCCOL = MOD( NPCOL+MYCOL-IYCOL, NPCOL )
            CALL PBCTRGET( ICONTXT, 'Row', 1, NQ0, ICEIL( NN, NB ),
     $                     WORK, 1, MCROW, MCCOL, IGD, MYROW, MYCOL,
     $                     NPROW, NPCOL )
          END IF
*
*         Broadcast a row block of WORK in every row of template
*
          IF( IYROW.EQ.-1 ) THEN
            IF( MYROW.EQ.MCROW ) THEN
              IF( LCMQ.GT.1 ) THEN
                KZ = 0
                IF( MYCOL.EQ.IYCOL ) KZ = NZ
                CALL PBCTRST1( ICONTXT, 'Row', NQ, NB, KZ, WORK, 1,
     $                         BETA, Y, INCY, LCMP, LCMQ, NQ0 )
              END IF
              CALL CGEBS2D( ICONTXT, 'Col', '1-tree', 1, NQ, Y, INCY )
            ELSE
              CALL CGEBR2D( ICONTXT, 'Col', '1-tree', 1, NQ, Y, INCY,
     $                      MCROW, MYCOL )
            END IF
*
*         Send a row block of WORK to the destination row
*
          ELSE
            IF( LCMQ.EQ.1 ) THEN
              IF( MYROW.EQ.MCROW ) THEN
                IF( MYROW.NE.IYROW )
     $            CALL CGESD2D( ICONTXT, 1, NQ0, WORK, 1, IYROW, MYCOL )
              ELSE IF( MYROW.EQ.IYROW ) THEN
                IF( BETA.EQ.ZERO ) THEN
                  CALL CGERV2D( ICONTXT, 1, NQ0, Y, INCY, MCROW, MYCOL )
                ELSE
                  CALL CGERV2D( ICONTXT, 1, NQ0, WORK, 1, MCROW, MYCOL )
                  CALL PBCVECADD( ICONTXT, 'G', NQ0, ONE, WORK, 1,
     $                            BETA, Y, INCY )
                END IF
              END IF
*
            ELSE
              NQ1 = NQ0 * MIN( LCMQ, MAX( 0, ICEIL(NN,NB)-MCCOL ) )
              IF( MYROW.EQ.MCROW ) THEN
                IF( MYROW.NE.IYROW )
     $            CALL CGESD2D( ICONTXT, 1, NQ1, WORK, 1, IYROW, MYCOL )
              ELSE IF( MYROW.EQ.IYROW ) THEN
                CALL CGERV2D( ICONTXT, 1, NQ1, WORK, 1, MCROW, MYCOL )
              END IF
*
              IF( MYROW.EQ.IYROW ) THEN
                KZ = 0
                IF( MYCOL.EQ.IYCOL ) KZ = NZ
                CALL PBCTRST1( ICONTXT, 'Row', NQ, NB, KZ, WORK, 1,
     $                         BETA, Y, INCY, LCMP, LCMQ, NQ0 )
              END IF
            END IF
          END IF
        END IF
*
*     When x is a row vector
*
      ELSE
*
*       Form  y <== x'  ( x is a row block )
*
*           ||
*           ||
*           ||               _____________
*          (y)      <==      -----(x)-----
*           ||
*           ||
*           ||
*
        IF(      IXROW.LT.-1 .OR. IXROW.GE.NPROW ) THEN
          INFO = 12
        ELSE IF( IXCOL.LT.0  .OR. IXCOL.GE.NPCOL ) THEN
          INFO = 13
        ELSE IF( IYROW.LT.0  .OR. IYROW.GE.NPROW ) THEN
          INFO = 14
        ELSE IF( IYCOL.LT.-1 .OR. IYCOL.GE.NPCOL ) THEN
          INFO = 15
        END IF
        IF( INFO.NE.0 ) GO TO 10
*
*       MRROW : row relative position in template from IYROW
*       MRCOL : column relative position in template from IXCOL
*
        MRROW = MOD( NPROW+MYROW-IYROW, NPROW )
        MRCOL = MOD( NPCOL+MYCOL-IXCOL, NPCOL )
        JYCOL = IYCOL
        IF( IYCOL.EQ.-1 ) JYCOL = IXCOL
*
        NP  = NUMROC( NN, NB, MYROW, IYROW, NPROW )
        IF( MRROW.EQ.0 ) NP = NP - NZ
        NQ  = NUMROC( NN, NB, MYCOL, IXCOL, NPCOL )
        IF( MRCOL.EQ.0 ) NQ = NQ - NZ
        NP0 = NUMROC( NUMROC(NN, NB, 0, 0, NPROW), NB, 0, 0, LCMP )
*
*       When a row process of IXROW has a row block A,
*
        IF( IXROW .GE. 0 ) THEN
          TBETA = ZERO
          IF( MYCOL.EQ.JYCOL ) TBETA = BETA
          KZ = NZ
*
          DO 40 I = 0, MIN( LCM, ICEIL(NN,NB) ) - 1
            MCROW = MOD( MOD(I, NPROW) + IYROW, NPROW )
            MCCOL = MOD( MOD(I, NPCOL) + IXCOL, NPCOL )
            IF( LCMP.EQ.1 ) NP0 = NUMROC( NN, NB, I, 0, NPROW )
            JDEX  = (I/NPROW) * NB
            IF( MRROW.EQ.0 ) JDEX = MAX(0, JDEX-NZ)
*
*           A source node copies the blocks to WORK, and send it
*
            IF( MYROW.EQ.IXROW .AND. MYCOL.EQ.MCCOL ) THEN
*
*             The source node is a destination node
*
              IDEX = (I/NPCOL) * NB
              IF( MRCOL.EQ.0 ) IDEX = MAX( 0, IDEX-NZ )
              IF( MYROW.EQ.MCROW .AND. MYCOL.EQ.JYCOL ) THEN
                CALL PBCTR2B1( ICONTXT, TRANS, NQ-IDEX, NB, KZ,
     $                         X(IDEX*INCX+1), INCX, TBETA,
     $                         Y(JDEX*INCY+1), INCY, LCMQ, LCMP )
*
*             The source node sends blocks to a destination node
*
              ELSE
                CALL PBCTR2B1( ICONTXT, TRANS, NQ-IDEX, NB, KZ,
     $                         X(IDEX*INCX+1), INCX, ZERO, WORK, 1,
     $                         LCMQ, 1 )
                CALL CGESD2D( ICONTXT, 1, NP0-KZ, WORK, 1,
     $                        MCROW, JYCOL )
              END IF
*
*           A destination node receives the copied blocks
*
            ELSE IF( MYROW.EQ.MCROW .AND. MYCOL.EQ.JYCOL ) THEN
              IF( LCMP.EQ.1 .AND. TBETA.EQ.ZERO ) THEN
                CALL CGERV2D( ICONTXT, 1, NP0-KZ, Y, INCY,
     $                        IXROW, MCCOL )
              ELSE
                CALL CGERV2D( ICONTXT, 1, NP0-KZ, WORK, 1,
     $                        IXROW, MCCOL )
                CALL PBCTR2A1( ICONTXT, NP-JDEX, NB, KZ, WORK, 1, TBETA,
     $                         Y(JDEX*INCY+1), INCY, LCMP*NB )
              END IF
            END IF
            KZ = 0
   40     CONTINUE
*
*         Broadcast a column vector Y in each row of template
*
          IF( IYCOL.EQ.-1 ) THEN
            IF( MYCOL.EQ.JYCOL ) THEN
              CALL CGEBS2D( ICONTXT, 'Row', '1-tree', 1, NP, Y, INCY )
            ELSE
              CALL CGEBR2D( ICONTXT, 'Row', '1-tree', 1, NP, Y, INCY,
     $                      MYROW, JYCOL )
            END IF
          END IF
*
*       When all row procesors have a copy of the row block A,
*
        ELSE
          IF( LCMP.EQ.1 ) NP0 = NP
*
*         Processors, which have diagonal blocks of A, copy them to
*         WORK array in transposed form
*
          KZ = 0
          IF( MRCOL.EQ.0 ) KZ = NZ
          JZ = 0
          IF( MRCOL.EQ.0 .AND. MYROW.EQ.IYROW ) JZ = NZ
*
          DO 50 I = 0, LCMQ-1
            IF( MRROW.EQ.MOD(NPCOL*I+MRCOL, NPROW) ) THEN
              IDEX = MAX( 0, I*NB-KZ )
              IF( LCMP.EQ.1 .AND. (IYCOL.EQ.-1.OR.IYCOL.EQ.MYCOL) ) THEN
                CALL PBCTR2B1( ICONTXT, TRANS, NQ-IDEX, NB, JZ,
     $                          X(IDEX*INCX+1), INCX, BETA, Y, INCY,
     $                          LCMQ, 1 )
              ELSE
                CALL PBCTR2B1( ICONTXT, TRANS, NQ-IDEX, NB, JZ,
     $                         X(IDEX*INCX+1), INCX, ZERO, WORK, 1,
     $                         LCMQ, 1 )
              END IF
            END IF
   50     CONTINUE
*
*         Get diagonal blocks of A for each row of the template
*
          MCCOL = MOD( MOD(MRROW, NPCOL) + IXCOL, NPCOL )
          IF( LCMP.GT.1 ) THEN
            MCROW = MOD( NPROW+MYROW-IYROW, NPROW )
            CALL PBCTRGET( ICONTXT, 'Col', 1, NP0, ICEIL( NN, NB ),
     $                     WORK, 1, MCROW, MCCOL, IGD, MYROW, MYCOL,
     $                     NPROW, NPCOL )
          END IF
*
*         Broadcast a column block of WORK in every column of template
*
          IF( IYCOL.EQ.-1 ) THEN
            IF( MYCOL.EQ.MCCOL ) THEN
              IF( LCMP.GT.1 ) THEN
                KZ = 0
                IF( MYROW.EQ.IYROW ) KZ = NZ
                CALL PBCTRST1( ICONTXT, 'Col', NP, NB, KZ, WORK, 1,
     $                         BETA, Y, INCY, LCMP, LCMQ, NP0 )
              END IF
              CALL CGEBS2D( ICONTXT, 'Row', '1-tree', 1, NP, Y, INCY )
            ELSE
              CALL CGEBR2D( ICONTXT, 'Row', '1-tree', 1, NP, Y, INCY,
     $                      MYROW, MCCOL )
            END IF
*
*         Send a column block of WORK to the destination column
*
          ELSE
            IF( LCMP.EQ.1 ) THEN
              IF( MYCOL.EQ.MCCOL ) THEN
                IF( MYCOL.NE.IYCOL )
     $            CALL CGESD2D( ICONTXT, 1, NP, WORK, 1, MYROW, IYCOL )
              ELSE IF( MYCOL.EQ.IYCOL ) THEN
                IF( BETA.EQ.ZERO ) THEN
                  CALL CGERV2D( ICONTXT, 1, NP, Y, INCY, MYROW, MCCOL )
                ELSE
                  CALL CGERV2D( ICONTXT, 1, NP, WORK, 1, MYROW, MCCOL )
                  CALL PBCVECADD( ICONTXT, 'G', NP, ONE, WORK, 1, BETA,
     $                            Y, INCY )
                END IF
              END IF
*
            ELSE
              NP1 = NP0 * MIN( LCMP, MAX( 0, ICEIL(NN,NB)-MCROW ) )
              IF( MYCOL.EQ.MCCOL ) THEN
                IF( MYCOL.NE.IYCOL )
     $            CALL CGESD2D( ICONTXT, 1, NP1, WORK, 1, MYROW, IYCOL )
              ELSE IF( MYCOL.EQ.IYCOL ) THEN
                CALL CGERV2D( ICONTXT, 1, NP1, WORK, 1, MYROW, MCCOL )
              END IF
*
              IF( MYCOL.EQ.IYCOL ) THEN
                KZ = 0
                IF( MYROW.EQ.IYROW ) KZ = NZ
                CALL PBCTRST1( ICONTXT, 'Col', NP, NB, KZ, WORK, 1,
     $                         BETA, Y, INCY, LCMP, LCMQ, NP0 )
              END IF
            END IF
          END IF
        END IF
      END IF
*
      RETURN
*
*     End of PBCTRNV
*
      END
*
*=======================================================================
*     SUBROUTINE PBCTR2A1
*=======================================================================
*
      SUBROUTINE PBCTR2A1( ICONTXT, N, NB, NZ, X, INCX, BETA, Y, INCY,
     $                     INTV )
*
*  -- PB-BLAS routine (version 2.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory.
*     April 28, 1996
*
*     .. Scalar Arguments ..
      INTEGER              ICONTXT, N, NB, NZ, INCX, INCY, INTV
      COMPLEX              BETA
*     ..
*     .. Array Arguments ..
      COMPLEX              X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*     y <== x
*     y is a scattered vector, copied from a condensed vector x.
*
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC            MIN
*     ..
*     .. External Functions ..
      INTEGER              ICEIL
      EXTERNAL             ICEIL
*     ..
*     .. External Subroutines ..
      EXTERNAL             PBCVECADD
*     ..
*     .. Parameters ..
      COMPLEX              ONE
      PARAMETER          ( ONE  = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Variables ..
      INTEGER              IX, IY, JZ, K, ITER
*
      IX = 0
      IY = 0
      JZ = NZ
      ITER = ICEIL( N+NZ, INTV )
*
      IF( ITER.GT.1 ) THEN
         CALL PBCVECADD( ICONTXT, 'G', NB-JZ, ONE, X(IX*INCX+1), INCX,
     $                   BETA, Y(IY*INCY+1), INCY )
         IX = IX + NB   - JZ
         IY = IY + INTV - JZ
         JZ = 0
*
         DO 10 K = 2, ITER-1
            CALL PBCVECADD( ICONTXT, 'G', NB, ONE, X(IX*INCX+1), INCX,
     $                      BETA, Y(IY*INCY+1), INCY )
            IX = IX + NB
            IY = IY + INTV
   10    CONTINUE
      END IF
*
      CALL PBCVECADD( ICONTXT, 'G', MIN( N-IY, NB-JZ ), ONE,
     $                X(IX*INCX+1), INCX, BETA, Y(IY*INCY+1), INCY )
*
      RETURN
*
*     End of PBCTR2A1
*
      END
*
*=======================================================================
*     SUBROUTINE PBCTR2B1
*=======================================================================
*
      SUBROUTINE PBCTR2B1( ICONTXT, TRANS, N, NB, NZ, X, INCX, BETA, Y,
     $                     INCY, JINX, JINY )
*
*  -- PB-BLAS routine (version 2.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory.
*     April 28, 1996
*
*     .. Scalar Arguments ..
      CHARACTER*1          TRANS
      INTEGER              ICONTXT, N, NB, NZ, INCX, INCY, JINX, JINY
      COMPLEX              BETA
*     ..
*     .. Array Arguments ..
      COMPLEX              X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*     y <== x + beta * y
*     y is a condensed vector, copied from a scattered vector x
*
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC            MIN
*     ..
*     .. External Functions ..
      INTEGER              ICEIL
      EXTERNAL             ICEIL
*     ..
*     .. External Subroutines ..
      EXTERNAL             PBCVECADD
*     ..
*     .. Parameters ..
      COMPLEX              ONE
      PARAMETER          ( ONE  = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Variables ..
      INTEGER              IX, IY, JZ, K, ITER, LENX, LENY
*
      IF( JINX.EQ.1 .AND. JINY.EQ.1 ) THEN
         CALL PBCVECADD( ICONTXT, TRANS, N, ONE, X, INCX, BETA,
     $                   Y, INCY )
*
      ELSE
         IX   = 0
         IY   = 0
         JZ   = NZ
         LENX = NB * JINX
         LENY = NB * JINY
         ITER = ICEIL( N+NZ, LENX )
*
         IF( ITER.GT.1 ) THEN
            CALL PBCVECADD( ICONTXT, TRANS, NB-JZ, ONE, X(IX*INCX+1),
     $                      INCX, BETA, Y(IY*INCY+1), INCY )
            IX = IX + LENX - JZ
            IY = IY + LENY - JZ
            JZ = 0
*
            DO 10 K = 2, ITER-1
               CALL PBCVECADD( ICONTXT, TRANS, NB, ONE, X(IX*INCX+1),
     $                         INCX, BETA, Y(IY*INCY+1), INCY )
               IX = IX + LENX
               IY = IY + LENY
   10       CONTINUE
         END IF
*
         CALL PBCVECADD( ICONTXT, TRANS, MIN( N-IX, NB-JZ ), ONE,
     $                   X(IX*INCX+1), INCX, BETA, Y(IY*INCY+1), INCY )
      END IF
*
      RETURN
*
*     End of PBCTR2B1
*
      END
