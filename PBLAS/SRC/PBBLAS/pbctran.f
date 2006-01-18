      SUBROUTINE PBCTRAN( ICONTXT, ADIST, TRANS, M, N, NB, A, LDA, BETA,
     $                    C, LDC, IAROW, IACOL, ICROW, ICCOL, WORK )
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
      CHARACTER*1        ADIST, TRANS
      INTEGER            IACOL, IAROW, ICCOL, ICONTXT, ICROW, LDA, LDC,
     $                   M, N, NB
      COMPLEX            BETA
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PBCTRAN  transposes  a column block to row block, or a row block to
*  column block by reallocating data distribution.
*
*     C := A^T + beta*C, or C := A^C + beta*C
*
*  where A is an M-by-N matrix  and C is an N-by-M matrix, and the size
*  of M or N is limited to its block size NB.
*
*  The first elements  of the matrices A, and C  should  be  located  at
*  the beginnings of their first blocks. (not the middle of the blocks.)
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
*  ADIST  - (input) CHARACTER*1
*           ADIST specifies whether A is a column block or a row block.
*
*              ADIST = 'C',  A is a column block
*              ADIST = 'R',  A is a row block
*
*  TRANS  - (input) CHARACTER*1
*           TRANS specifies whether the transposed format is transpose
*           or conjugate transpose.  If the matrices A and C are real,
*           the argument is ignored.
*
*              TRANS = 'T',  transpose
*              TRANS = 'C',  conjugate transpose
*
*  M      - (input) INTEGER
*           M specifies the (global) number of rows of the matrix (block
*           column or block row) A and of columns of the matrix C.
*           M >= 0.
*
*  N      - (input) INTEGER
*           N specifies the (global) number of columns of the matrix
*           (block column or block row) A  and of columns of the matrix
*           C.  N >= 0.
*
*  NB     - (input) INTEGER
*           NB specifies  the column block size of the matrix A and the
*           row block size of the matrix C when ADIST = 'C'.  Otherwise,
*           it specifies  the row block size of the matrix A and the
*           column block size of the matrix C. NB >= 1.
*
*  A       (input) COMPLEX array of DIMENSION ( LDA, Lx ),
*          where Lx is N  when ADIST = 'C', or Nq when ADIST = 'R'.
*          Before entry with  ADIST = 'C',  the leading Mp by N part of
*          the array A must contain the matrix A, otherwise the leading
*          M by Nq part of the array A  must contain the matrix A.  See
*          parameter details for the values of Mp and Nq.
*
*  LDA     (input) INTEGER
*          LDA specifies the leading dimension of (local) A as declared
*          in the calling (sub) program.  LDA >= MAX(1,Mp) when
*          ADIST = 'C', or LDA >= MAX(1,M) otherwise.
*
*  BETA    (input) COMPLEX
*          BETA specifies scaler beta.
*
*  C       (input/output) COMPLEX array of DIMENSION ( LDC, Lx ),
*          where Lx is Mq when ADIST = 'C', or N when ADIST = 'R'.
*          If ADIST = 'C', the leading N-by-Mq part of the array C
*          contains the (local) matrix C, otherwise the leading
*          Np-by-M part of the array C must contain the (local) matrix
*          C.  C will not be referenced if beta is zero.
*
*  LDC     (input) INTEGER
*          LDC specifies the leading dimension of (local) C as declared
*          in the calling (sub) program. LDC >= MAX(1,N) when ADIST='C',
*          or LDC >= MAX(1,Np) otherwise.
*
*  IAROW   (input) INTEGER
*          IAROW specifies  a row  of the process  template,
*          which holds the first block  of the matrix  A. If A is a row
*          of blocks (ADIST = 'R') and all rows of processes have a copy
*          of A, then set IAROW = -1.
*
*  IACOL   (input) INTEGER
*          IACOL specifies  a column of the process template,
*          which holds  the first block  of the matrix A.  If  A is  a
*          column of blocks (ADIST = 'C') and all columns of processes
*          have a copy of A, then set IACOL = -1.
*
*  ICROW   (input) INTEGER
*          ICROW specifies the current row process which holds
*          the first block  of the matrix C,  which is transposed of A.
*          If C is a row of blocks (ADIST = 'C') and the transposed
*          row block C is distributed all rows of processes, set
*          ICROW = -1.
*
*  ICCOL   (input) INTEGER
*          ICCOL specifies  the current column process which holds
*          the first block of the matrix C,  which is transposed of A.
*          If C is a column of blocks (ADIST = 'R') and the transposed
*          column block C is distributed all columns of processes,
*          set ICCOL = -1.
*
*  WORK    (workspace) COMPLEX array of dimension Size(WORK).
*          It needs extra working space of A'.
*
*  Parameters Details
*  ==================
*
*  Lx      It is  a local portion  of L  owned  by  a process,  (L is
*          replaced by M, or N,  and x is replaced by either p (=NPROW)
*          or q (=NPCOL)).  The value is  determined by  L, LB, x,  and
*          MI, where  LB is  a block size  and  MI is a  row  or column
*          position  in a process template.  Lx is  equal to  or less
*          than Lx0 = CEIL( L, LB*x ) * LB.
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
*  Mqb  = CEIL( M, NB*NPCOL )
*  Npb  = CEIL( N, NB*NPROW )
*  LCMQ = LCM / NPCOL
*  LCMP = LCM / NPROW
*
*  (1) ADIST = 'C'
*   (a) IACOL != -1
*       Size(WORK) = N * CEIL(Mqb,LCMQ)*NB
*   (b) IACOL = -1
*       Size(WORK) = N * CEIL(Mqb,LCMQ)*NB * MIN(LCMQ,CEIL(M,NB))
*
*  (2) ADIST = 'R'
*   (a) IAROW != -1
*       Size(WORK) = M * CEIL(Npb,LCMP)*NB
*   (b) IAROW = -1
*       Size(WORK) = M * CEIL(Npb,LCMP)*NB * MIN(LCMP,CEIL(N,NB))
*
*  Notes
*  -----
*  More precise space can be computed as
*
*  CEIL(Mqb,LCMQ)*NB => NUMROC( NUMROC(M,NB,0,0,NPCOL), NB, 0, 0, LCMQ )
*  CEIL(Npb,LCMP)*NB => NUMROC( NUMROC(N,NB,0,0,NPROW), NB, 0, 0, LCMP )
*
*  =====================================================================
*
*     ..
*     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE  = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLFORM, ROWFORM
      INTEGER            I, IDEX, IGD, INFO, JCCOL, JCROW, JDEX, LCM,
     $                   LCMP, LCMQ, MCCOL, MCROW, ML, MP, MQ, MQ0,
     $                   MRCOL, MRROW, MYCOL, MYROW, NP, NP0, NPCOL,
     $                   NPROW, NQ
      COMPLEX            TBETA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILCM, ICEIL, NUMROC
      EXTERNAL           ILCM, ICEIL, LSAME, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CGEBR2D, CGEBS2D, CGERV2D,
     $                   CGESD2D, PBCMATADD, PBCTR2AF, PBCTR2AT,
     $                   PBCTR2BT, PBCTRGET, PBCTRSRT, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
*
      CALL BLACS_GRIDINFO( ICONTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      COLFORM = LSAME( ADIST, 'C' )
      ROWFORM = LSAME( ADIST, 'R' )
*
*     Test the input parameters.
*
      INFO = 0
      IF( ( .NOT.COLFORM ) .AND. ( .NOT.ROWFORM ) ) THEN
         INFO = 2
      ELSE IF( M .LT.0                            ) THEN
         INFO = 4
      ELSE IF( N .LT.0                            ) THEN
         INFO = 5
      ELSE IF( NB.LT.1                            ) THEN
         INFO = 6
      ELSE IF( IAROW.LT.-1 .OR. IAROW.GE.NPROW .OR.
     $       ( IAROW.EQ.-1 .AND. COLFORM )        ) THEN
         INFO = 12
      ELSE IF( IACOL.LT.-1 .OR. IACOL.GE.NPCOL .OR.
     $       ( IACOL.EQ.-1 .AND. ROWFORM )        ) THEN
         INFO = 13
      ELSE IF( ICROW.LT.-1 .OR. ICROW.GE.NPROW .OR.
     $       ( ICROW.EQ.-1 .AND. ROWFORM )        ) THEN
         INFO = 14
      ELSE IF( ICCOL.LT.-1 .OR. ICCOL.GE.NPCOL .OR.
     $       ( ICCOL.EQ.-1 .AND. COLFORM )        ) THEN
         INFO = 15
      END IF
*
   10 CONTINUE
      IF( INFO .NE. 0 ) THEN
         CALL PXERBLA( ICONTXT, 'PBCTRAN ', INFO )
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
*
*     When A is a column block
*
      IF( COLFORM ) THEN
*
*       Form  C <== A'  ( A is a column block )
*                                         _
*                                        | |
*                                        | |
*            _____________               | |
*           |______C______|     <==      |A|
*                                        | |
*                                        | |
*                                        |_|
*
*       MRROW : row relative position in template from IAROW
*       MRCOL : column relative position in template from ICCOL
*
        MRROW = MOD( NPROW+MYROW-IAROW, NPROW )
        MRCOL = MOD( NPCOL+MYCOL-ICCOL, NPCOL )
        JCROW = ICROW
        IF( ICROW.EQ.-1 ) JCROW = IAROW
*
        MP  = NUMROC( M, NB, MYROW, IAROW, NPROW )
        MQ  = NUMROC( M, NB, MYCOL, ICCOL, NPCOL )
        MQ0 = NUMROC( NUMROC(M, NB, 0, 0, NPCOL), NB, 0, 0, LCMQ )
*
        IF( LDA.LT.MP .AND.
     $         ( IACOL.EQ.MYCOL .OR. IACOL.EQ.-1 ) ) THEN
           INFO = 8
        ELSE IF( LDC.LT.N .AND.
     $         ( ICROW.EQ.MYROW .OR. ICROW.EQ.-1 ) ) THEN
           INFO = 11
        END IF
        IF( INFO.NE.0 ) GO TO 10
*
*       When a column process of IACOL has a column block A,
*
        IF( IACOL.GE.0 ) THEN
          TBETA = ZERO
          IF( MYROW.EQ.JCROW ) TBETA = BETA
*
          DO 20 I = 0, MIN( LCM, ICEIL(M,NB) ) - 1
            MCROW = MOD( MOD(I, NPROW) + IAROW, NPROW )
            MCCOL = MOD( MOD(I, NPCOL) + ICCOL, NPCOL )
            IF( LCMQ.EQ.1 )  MQ0 = NUMROC( M, NB, I, 0, NPCOL )
            JDEX = (I/NPCOL) * NB
*
*           A source node copies the blocks to WORK, and send it
*
            IF( MYROW.EQ.MCROW .AND. MYCOL.EQ.IACOL ) THEN
*
*             The source node is a destination node
*
              IDEX = (I/NPROW) * NB
              IF( MYROW.EQ.JCROW .AND. MYCOL.EQ.MCCOL ) THEN
                CALL PBCTR2AT( ICONTXT, 'Col', TRANS, MP-IDEX, N, NB,
     $                         A(IDEX+1,1), LDA, TBETA, C(1,JDEX+1),
     $                         LDC, LCMP, LCMQ )
*
*             The source node sends blocks to a destination node
*
              ELSE
                CALL PBCTR2BT( ICONTXT, 'Col', TRANS, MP-IDEX, N, NB,
     $                         A(IDEX+1,1), LDA, ZERO, WORK, N,
     $                         LCMP*NB )
                CALL CGESD2D( ICONTXT, N, MQ0, WORK, N, JCROW, MCCOL )
              END IF
*
*           A destination node receives the copied blocks
*
            ELSE IF( MYROW.EQ.JCROW .AND. MYCOL.EQ.MCCOL ) THEN
              IF( LCMQ.EQ.1 .AND. TBETA.EQ.ZERO ) THEN
                CALL CGERV2D( ICONTXT, N, MQ0, C, LDC, MCROW, IACOL )
              ELSE
                CALL CGERV2D( ICONTXT, N, MQ0, WORK, N, MCROW, IACOL )
                CALL PBCTR2AF( ICONTXT, 'Row', N, MQ-JDEX, NB, WORK, N,
     $                         TBETA, C(1,JDEX+1), LDC, LCMP, LCMQ,
     $                         MQ0 )
              END IF
            END IF
   20     CONTINUE
*
*         Broadcast a row block of C in each column of template
*
          IF( ICROW.EQ.-1 ) THEN
            IF( MYROW.EQ.JCROW ) THEN
              CALL CGEBS2D( ICONTXT, 'Col', '1-tree', N, MQ, C, LDC )
            ELSE
              CALL CGEBR2D( ICONTXT, 'Col', '1-tree', N, MQ, C, LDC,
     $                      JCROW, MYCOL )
            END IF
          END IF
*
*       When all column procesors have a copy of the column block A,
*
        ELSE
          IF( LCMQ.EQ.1 ) MQ0 = MQ
*
*         Processors, which have diagonal blocks of A, copy them to
*         WORK array in transposed form
*
          DO 30 I = 0, LCMP-1
            IF( MRCOL.EQ.MOD( NPROW*I+MRROW, NPCOL ) ) THEN
              IF( LCMQ.EQ.1.AND.(ICROW.EQ.-1.OR.ICROW.EQ.MYROW) ) THEN
                 CALL PBCTR2BT( ICONTXT, 'Col', TRANS, MP-I*NB, N, NB,
     $                          A(I*NB+1,1), LDA, BETA, C, LDC,
     $                          LCMP*NB )
              ELSE
                 CALL PBCTR2BT( ICONTXT, 'Col', TRANS, MP-I*NB, N, NB,
     $                          A(I*NB+1,1), LDA, ZERO, WORK, N,
     $                          LCMP*NB )
              END IF
            END IF
   30     CONTINUE
*
*         Get diagonal blocks of A for each column of the template
*
          MCROW = MOD( MOD(MRCOL,NPROW)+IAROW, NPROW )
          IF( LCMQ.GT.1 ) THEN
            MCCOL = MOD( NPCOL+MYCOL-ICCOL, NPCOL )
            CALL PBCTRGET( ICONTXT, 'Row', N, MQ0, ICEIL(M,NB), WORK, N,
     $                     MCROW,  MCCOL, IGD, MYROW, MYCOL, NPROW,
     $                     NPCOL )
          END IF
*
*         Broadcast a row block of WORK in every row of template
*
          IF( ICROW.EQ.-1 ) THEN
            IF( MYROW.EQ.MCROW ) THEN
              IF( LCMQ.GT.1 )
     $          CALL PBCTRSRT( ICONTXT, 'Row', N, MQ, NB, WORK, N, BETA,
     $                         C, LDC, LCMP, LCMQ, MQ0 )
              CALL CGEBS2D( ICONTXT, 'Col', '1-tree', N, MQ, C, LDC )
            ELSE
              CALL CGEBR2D( ICONTXT, 'Col', '1-tree', N, MQ, C, LDC,
     $                      MCROW, MYCOL )
            END IF
*
*         Send a row block of WORK to the destination row
*
          ELSE
            IF( LCMQ.EQ.1 ) THEN
              IF( MYROW.EQ.MCROW ) THEN
                IF( MYROW.NE.ICROW )
     $            CALL CGESD2D( ICONTXT, N, MQ, WORK, N, ICROW, MYCOL )
              ELSE IF( MYROW.EQ.ICROW ) THEN
                IF( BETA.EQ.ZERO ) THEN
                  CALL CGERV2D( ICONTXT, N, MQ, C, LDC, MCROW, MYCOL )
                ELSE
                  CALL CGERV2D( ICONTXT, N, MQ, WORK, N, MCROW, MYCOL )
                  CALL PBCMATADD( ICONTXT, 'G', N, MQ, ONE, WORK, N,
     $                            BETA, C, LDC )
                END IF
              END IF
*
            ELSE
              ML = MQ0 * MIN( LCMQ, MAX(0,ICEIL(M,NB)-MCCOL) )
              IF( MYROW.EQ.MCROW ) THEN
                IF( MYROW.NE.ICROW )
     $            CALL CGESD2D( ICONTXT, N, ML, WORK, N, ICROW, MYCOL )
              ELSE IF( MYROW.EQ.ICROW ) THEN
                CALL CGERV2D( ICONTXT, N, ML, WORK, N, MCROW, MYCOL )
              END IF
*
              IF( MYROW.EQ.ICROW )
     $          CALL PBCTRSRT( ICONTXT, 'Row', N, MQ, NB, WORK, N, BETA,
     $                         C, LDC, LCMP, LCMQ, MQ0 )
            END IF
          END IF
*
        END IF
*
*     When A is a row block
*
      ELSE
*
*        Form  C <== A'  ( A is a row block )
*            _
*           | |
*           | |
*           | |                _____________
*           |C|      <==      |______A______|
*           | |
*           | |
*           |_|
*
*        MRROW : row relative position in template from ICROW
*        MRCOL : column relative position in template from IACOL
*
         MRROW = MOD( NPROW+MYROW-ICROW, NPROW )
         MRCOL = MOD( NPCOL+MYCOL-IACOL, NPCOL )
         JCCOL = ICCOL
         IF( ICCOL.EQ.-1 ) JCCOL = IACOL
*
         NP  = NUMROC( N, NB, MYROW, ICROW, NPROW )
         NQ  = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
         NP0 = NUMROC( NUMROC(N, NB, 0, 0, NPROW), NB, 0, 0, LCMP )
*
         IF( LDA.LT.M .AND.
     $          ( IAROW.EQ.MYROW .OR. IAROW.EQ.-1 ) ) THEN
            INFO = 8
         ELSE IF( LDC.LT.NP .AND.
     $          ( ICCOL.EQ.MYCOL .OR. ICCOL.EQ.-1 ) ) THEN
            INFO = 11
         END IF
         IF( INFO.NE.0 ) GO TO 10
*
*        When a row process of IAROW has a row block A,
*
         IF( IAROW.GE.0 ) THEN
           TBETA = ZERO
           IF( MYCOL.EQ.JCCOL ) TBETA = BETA
*
           DO 40 I = 0, MIN( LCM, ICEIL(N,NB) ) - 1
             MCROW = MOD( MOD(I, NPROW) + ICROW, NPROW )
             MCCOL = MOD( MOD(I, NPCOL) + IACOL, NPCOL )
             IF( LCMP.EQ.1 )  NP0 = NUMROC( N, NB, I, 0, NPROW )
             IDEX = (I/NPROW) * NB
*
*            A source node copies the blocks to WORK, and send it
*
             IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.MCCOL ) THEN
*
*              The source node is a destination node
*
               JDEX = (I/NPCOL) * NB
               IF( MYROW.EQ.MCROW .AND. MYCOL.EQ.JCCOL ) THEN
                 CALL PBCTR2AT( ICONTXT, 'Row', TRANS, M, NQ-JDEX, NB,
     $                          A(1,JDEX+1), LDA, TBETA, C(IDEX+1,1),
     $                          LDC, LCMP, LCMQ )
*
*              The source node sends blocks to a destination node
*
               ELSE
                 CALL PBCTR2BT( ICONTXT, 'Row', TRANS, M, NQ-JDEX, NB,
     $                          A(1,JDEX+1), LDA, ZERO, WORK, NP0,
     $                          LCMQ*NB )
                 CALL CGESD2D( ICONTXT, NP0, M, WORK, NP0,
     $                         MCROW, JCCOL )
               END IF
*
*           A destination node receives the copied blocks
*
            ELSE IF( MYROW.EQ.MCROW .AND. MYCOL.EQ.JCCOL ) THEN
              IF( LCMP.EQ.1 .AND. TBETA.EQ.ZERO ) THEN
                CALL CGERV2D( ICONTXT, NP0, M, C, LDC, IAROW, MCCOL )
              ELSE
                CALL CGERV2D( ICONTXT, NP0, M, WORK, NP0, IAROW, MCCOL )
                CALL PBCTR2AF( ICONTXT, 'Col', NP-IDEX, M, NB, WORK,
     $                         NP0, TBETA, C(IDEX+1,1), LDC, LCMP, LCMQ,
     $                         NP0 )
              END IF
            END IF
   40     CONTINUE
*
*         Broadcast a column block of WORK in each row of template
*
          IF( ICCOL.EQ.-1 ) THEN
            IF( MYCOL.EQ.JCCOL ) THEN
              CALL CGEBS2D( ICONTXT, 'Row', '1-tree', NP, M, C, LDC )
            ELSE
              CALL CGEBR2D( ICONTXT, 'Row', '1-tree', NP, M, C, LDC,
     $                       MYROW, JCCOL )
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
          DO 50 I = 0, LCMQ-1
            IF( MRROW.EQ.MOD(NPCOL*I+MRCOL, NPROW) ) THEN
              IF( LCMP.EQ.1.AND.(ICCOL.EQ.-1.OR.ICCOL.EQ.MYCOL) ) THEN
                CALL PBCTR2BT( ICONTXT, 'Row', TRANS, M, NQ-I*NB, NB,
     $                         A(1,I*NB+1), LDA, BETA, C, LDC,
     $                         LCMQ*NB )
              ELSE
                CALL PBCTR2BT( ICONTXT, 'Row', TRANS, M, NQ-I*NB, NB,
     $                         A(1,I*NB+1), LDA, ZERO, WORK, NP0,
     $                         LCMQ*NB )
              END IF
            END IF
   50     CONTINUE
*
*         Get diagonal blocks of A for each row of the template
*
          MCCOL = MOD( MOD(MRROW, NPCOL)+IACOL, NPCOL )
          IF( LCMP.GT.1 ) THEN
            MCROW = MOD( NPROW+MYROW-ICROW, NPROW )
            CALL PBCTRGET( ICONTXT, 'Col', NP0, M, ICEIL(N,NB), WORK,
     $                     NP0, MCROW, MCCOL, IGD, MYROW, MYCOL, NPROW,
     $                     NPCOL )
          END IF
*
*         Broadcast a column block of WORK in every column of template
*
          IF( ICCOL.EQ.-1 ) THEN
            IF( MYCOL.EQ.MCCOL ) THEN
              IF( LCMP.GT.1 )
     $          CALL PBCTRSRT( ICONTXT, 'Col', NP, M, NB, WORK, NP0,
     $                         BETA, C, LDC, LCMP, LCMQ, NP0 )
              CALL CGEBS2D( ICONTXT, 'Row', '1-tree', NP, M, C, LDC )
            ELSE
              CALL CGEBR2D( ICONTXT, 'Row', '1-tree', NP, M, C, LDC,
     $                       MYROW, MCCOL )
            END IF
*
*         Send a column block of WORK to the destination column
*
          ELSE
            IF( LCMP.EQ.1 ) THEN
              IF( MYCOL.EQ.MCCOL ) THEN
                IF( MYCOL.NE.ICCOL )
     $            CALL CGESD2D( ICONTXT, NP, M, WORK, NP, MYROW, ICCOL )
              ELSE IF( MYCOL.EQ.ICCOL ) THEN
                IF( BETA.EQ.ZERO ) THEN
                  CALL CGERV2D( ICONTXT, NP, M, C, LDC, MYROW, MCCOL )
                ELSE
                  CALL CGERV2D( ICONTXT, NP, M, WORK, NP, MYROW, MCCOL )
                  CALL PBCMATADD( ICONTXT, 'G', NP, M, ONE, WORK, NP,
     $                            BETA, C, LDC )
                END IF
              END IF
*
            ELSE
              ML = M * MIN( LCMP, MAX( 0, ICEIL(N,NB) - MCROW ) )
              IF( MYCOL.EQ.MCCOL ) THEN
                IF( MYCOL.NE.ICCOL )
     $            CALL CGESD2D( ICONTXT, NP0, ML, WORK, NP0,
     $                          MYROW, ICCOL )
              ELSE IF( MYCOL.EQ.ICCOL ) THEN
                CALL CGERV2D( ICONTXT, NP0, ML, WORK, NP0,
     $                        MYROW, MCCOL )
              END IF
*
              IF( MYCOL.EQ.ICCOL )
     $          CALL PBCTRSRT( ICONTXT, 'Col', NP, M, NB, WORK, NP0,
     $                         BETA, C, LDC, LCMP, LCMQ, NP0 )
            END IF
          END IF
*
        END IF
      END IF
*
      RETURN
*
*     End of PBCTRAN
*
      END
*
*=======================================================================
*     SUBROUTINE PBCTR2AT
*=======================================================================
*
      SUBROUTINE PBCTR2AT( ICONTXT, ADIST, TRANS, M, N, NB, A, LDA,
     $                     BETA, B, LDB, LCMP, LCMQ )
*
*  -- PB-BLAS routine (version 2.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory.
*     April 28, 1996
*
*     .. Scalar Arguments ..
      CHARACTER*1        ADIST, TRANS
      INTEGER            ICONTXT, LCMP, LCMQ, LDA, LDB, M, N, NB
      COMPLEX            BETA
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  PBCTR2AT forms   B <== A^T + beta*B, or A^C + beta*B
*  B is a ((conjugate) transposed) scattered block row (or column),
*  copied from a scattered block column (or row) of A
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE  = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            IA, IB, K, INTV, JNTV
*     ..
*     .. External Subroutines ..
      EXTERNAL           PBCMATADD
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL
      EXTERNAL           LSAME, ICEIL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Excutable Statements ..
*
      IF( LCMP.EQ.LCMQ ) THEN
         CALL PBCMATADD( ICONTXT, TRANS, N, M, ONE, A, LDA, BETA, B,
     $                   LDB )
*
      ELSE
*
*        If A is a column block ( ADIST = 'C' ),
*
         IF( LSAME( ADIST, 'C' ) ) THEN
            INTV = LCMP * NB
            JNTV = LCMQ * NB
            IA = 1
            IB = 1
            DO 10 K = 1, ICEIL( M, INTV )
               CALL PBCMATADD( ICONTXT, TRANS, N, MIN( M-IA+1, NB ),
     $                         ONE, A(IA,1), LDA, BETA, B(1,IB), LDB )
               IA = IA + INTV
               IB = IB + JNTV
   10       CONTINUE
*
*        If A is a row block ( ADIST = 'R' ),
*
         ELSE
            INTV = LCMP * NB
            JNTV = LCMQ * NB
            IA = 1
            IB = 1
            DO 20 K = 1, ICEIL( N, JNTV )
               CALL PBCMATADD( ICONTXT, TRANS, MIN( N-IA+1, NB ), M,
     $                         ONE, A(1,IA), LDA, BETA, B(IB,1), LDB )
               IA = IA + JNTV
               IB = IB + INTV
   20       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of PBCTR2AT
*
      END
*
*=======================================================================
*     SUBROUTINE PBCTR2BT
*=======================================================================
*
      SUBROUTINE PBCTR2BT( ICONTXT, ADIST, TRANS, M, N, NB, A, LDA,
     $                     BETA, B, LDB, INTV )
*
*  -- PB-BLAS routine (version 2.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory.
*     April 28, 1996
*
*     .. Scalar Arguments ..
      CHARACTER*1        ADIST, TRANS
      INTEGER            ICONTXT, INTV, LDA, LDB, M, N, NB
      COMPLEX            BETA
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  PBCTR2BT forms T <== A^T + beta*T or A^C + beta*T, where T is a
*  ((conjugate) transposed) condensed block row (or column), copied from
*  a scattered block column (or row) of A
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE  = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            IA, IB, K
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL
      EXTERNAL           LSAME, ICEIL
*     ..
*     .. External Subroutines ..
      EXTERNAL           PBCMATADD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Excutable Statements ..
*
      IF( INTV.EQ.NB ) THEN
         CALL PBCMATADD( ICONTXT, TRANS, N, M, ONE, A, LDA, BETA, B,
     $                   LDB )
*
      ELSE
*
*        If A is a column block ( ADIST = 'C' ),
*
         IF( LSAME( ADIST, 'C' ) ) THEN
            IA = 1
            IB = 1
            DO 10 K = 1, ICEIL( M, INTV )
               CALL PBCMATADD( ICONTXT, TRANS, N, MIN( M-IA+1, NB ),
     $                         ONE, A(IA,1), LDA, BETA, B(1,IB), LDB )
               IA = IA + INTV
               IB = IB + NB
   10       CONTINUE
*
*        If A is a row block (ADIST = 'R'),
*
         ELSE
            IA = 1
            IB = 1
            DO 20 K = 1, ICEIL( N, INTV )
               CALL PBCMATADD( ICONTXT, TRANS, MIN( N-IA+1, NB ), M,
     $                         ONE, A(1,IA), LDA, BETA, B(IB,1), LDB )
               IA = IA + INTV
               IB = IB + NB
   20       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of PBCTR2BT
*
      END
*
*=======================================================================
*     SUBROUTINE PBCTR2AF
*=======================================================================
*
      SUBROUTINE PBCTR2AF( ICONTXT, ADIST, M, N, NB, A, LDA, BETA, B,
     $                     LDB, LCMP, LCMQ, NINT )
*
*  -- PB-BLAS routine (version 2.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory.
*     April 28, 1996
*
*     .. Scalar Arguments ..
      CHARACTER*1          ADIST
      INTEGER              ICONTXT, M, N, NB, LDA, LDB, LCMP, LCMQ, NINT
      COMPLEX              BETA
*     ..
*     .. Array Arguments ..
      COMPLEX              A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  PBCTR2AF forms  T <== A + BETA*T, where T is a scattered block
*  row (or column) copied from a (condensed) block column (or row) of A
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            JA, JB, K, INTV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL
      EXTERNAL           LSAME, ICEIL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( ADIST, 'R' ) ) THEN
         INTV = NB * LCMQ
         JA = 1
         JB = 1
         DO 10 K = 1, ICEIL( NINT, NB )
            CALL PBCMATADD( ICONTXT, 'G', M, MIN( N-JB+1, NB ), ONE,
     $                      A(1,JA), LDA, BETA, B(1,JB), LDB )
            JA = JA + NB
            JB = JB + INTV
   10    CONTINUE
*
*     if( LSAME( ADIST, 'C' ) ) then
*
      ELSE
         INTV = NB * LCMP
         JA = 1
         JB = 1
         DO 20 K = 1, ICEIL( NINT, NB )
            CALL PBCMATADD( ICONTXT, 'G', MIN( M-JB+1, NB ), N, ONE,
     $                      A(JA,1), LDA, BETA, B(JB,1), LDB )
            JA = JA + NB
            JB = JB + INTV
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of PBCTR2AF
*
      END
