***********************************************************************
*     Test program for ScaLAPACK-style routine PDHSEQR                *
***********************************************************************
*
*     Contributor: Robert Granat and Meiyue Shao
*     This version is of Feb 2011.
*
      PROGRAM PDHSEQRDRIVER
*
*     Declarations
*
      IMPLICIT NONE
*     ...Parameters...
      LOGICAL           BALANCE, COMPHESS, COMPRESI,
     $                  COMPORTH
      LOGICAL           DEBUG, PRN, TIMESTEPS, BARR,
     $                  UNI_LAPACK
      INTEGER           SLV_MIN, SLV_MAX
      PARAMETER         ( DEBUG = .FALSE.,
     $                    PRN = .FALSE.,
     $                    TIMESTEPS = .TRUE.,
     $                    COMPHESS = .TRUE.,
     $                    COMPRESI = .TRUE.,
     $                    COMPORTH = .TRUE.,
     $                    BALANCE = .TRUE.,
     $                    BARR = .FALSE.,
*     Solver: 1-PDLAQR1, 2-PDHSEQR.
     $                    SLV_MIN = 2, SLV_MAX = 2,
     $                    UNI_LAPACK = .TRUE. )
      INTEGER           N, NB, ARSRC, ACSRC
      PARAMETER         (
*     Problem size.
     $                    N = 500, NB = 50,
*     What processor should hold the first element in A?
     $                    ARSRC = 0, ACSRC = 0 )
      INTEGER           BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                  LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER         ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                    CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                    RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER           DPALLOC, INTALLC
      INTEGER           DPSIZ, INTSZ, NOUT, IZERO
      PARAMETER         ( DPSIZ = 8, DPALLOC = 8 000 000,
     $                    INTSZ = 4, INTALLC = 8 000 000,
     $                    NOUT = 6, IZERO = 0)
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D+00, ONE = 1.0D+00, TWO = 2.0D+00 )
*
*     ...Local Scalars...
      INTEGER           ICTXT, IAM, NPROW, NPCOL, MYROW, MYCOL,
     $                  SYS_NPROCS, NPROCS, AROWS, ACOLS, TEMP_ICTXT
      INTEGER           THREADS
      INTEGER           INFO, KTOP, KBOT, ILO, IHI, I
      INTEGER           IPA, IPACPY, IPQ, WR1, WI1, WR2, WI2, IPW1,
     $                  IPW2, IPIW
      INTEGER           TOTIT, SWEEPS, TOTNS, HESS
      DOUBLE PRECISION  EPS, THRESH
      DOUBLE PRECISION  STAMP, TOTTIME, T_BA, T_GEN, T_HS, T_SCH, T_QR,
     $                  T_RES, ITPEREIG, SWPSPEIG, NSPEIG, SPEEDUP, 
     $                  EFFICIENCY
      DOUBLE PRECISION  RNORM, ANORM, R1, ORTH, O1, O2, DPDUM, ELEM1,
     $                  ELEM2, ELEM3, EDIFF
      INTEGER           SOLVER
      CHARACTER*6       PASSED
*
*     ...Local Arrays...
      INTEGER           DESCA( DLEN_ ), DESCQ( DLEN_ ), DESCVEC( DLEN_ )
      DOUBLE PRECISION  SCALE( N )
      DOUBLE PRECISION, ALLOCATABLE :: MEM(:)
      INTEGER, ALLOCATABLE :: IMEM(:)
*
*     ...Intrinsic Functions...
      INTRINSIC         INT, DBLE, SQRT, MAX, MIN
*
*     ...External Functions...
      INTEGER           NUMROC
      DOUBLE PRECISION  PDLAMCH, PDLANGE, MPI_WTIME
      EXTERNAL          BLACS_PINFO, BLACS_GET, BLACS_GRIDINIT,
     $                  BLACS_GRIDINFO, BLACS_GRIDEXIT, BLACS_EXIT
      EXTERNAL          NUMROC, PDLAMCH, PDLASET, PDGEHRD, PDLANGE
      EXTERNAL          DGEBAL, DGEHRD
      EXTERNAL          MPI_WTIME
      EXTERNAL          PDGEBAL
      EXTERNAL          PDMATGEN2
*
*     ...Executable statements...
*
      CALL BLACS_PINFO( IAM, SYS_NPROCS )
      NPROW = INT( SQRT( DBLE(SYS_NPROCS) ) )
      NPCOL = SYS_NPROCS / NPROW
      CALL BLACS_GET( 0, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, '2D', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
c      print*, iam, ictxt, myrow, mycol
c      IF ( ( MYROW.GE.NPROW ) .OR. ( MYCOL.GE.NPCOL ) ) GO TO 777
      IF ( ICTXT.LT.0 ) GO TO 777
*
*     Read out the number of underlying threads and set stack size in
*     kilobytes.
*
	THRESH = 30.0
      TOTTIME = MPI_WTIME()
      T_GEN = 0.0D+00
      T_RES = 0.0D+00
      T_SCH = 0.0D+00
*
*     Allocate and Init memory with zeros.
*
      INFO = 0
      ALLOCATE ( MEM( DPALLOC ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate MEM. INFO = ', INFO
         GO TO 777
      END IF
      ALLOCATE ( IMEM( INTALLC ), STAT = INFO )
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) '% Could not allocate IMEM. INFO = ', INFO
         GO TO 777
      END IF
      MEM( 1:DPALLOC ) = ZERO
      IMEM( 1:INTALLC ) = IZERO
*
*     Get machine epsilon.
*
      EPS = PDLAMCH( ICTXT, 'Epsilon' )      
*
*     Print welcoming message.
*
      IF( IAM.EQ.0 ) THEN
         WRITE(*,*)
         WRITE(*,*) 'ScaLAPACK Test for PDHSEQR'
         WRITE(*,*) 
         WRITE(*,*) 'epsilon   = ', EPS
         WRITE(*,*) 'threshold = ', THRESH
         WRITE(*,*)
         WRITE(*,*) 'Residual and Orthogonality Residual computed by:'
         WRITE(*,*)
         WRITE(*,*) 'Residual      = ',
     $   ' || T - Q^T*A*Q ||_F / ( ||A||_F * eps * sqrt(N) )'
     	   WRITE(*,*)
         WRITE(*,*) 'Orthogonality = ',
     $   ' MAX( || I - Q^T*Q ||_F, || I - Q*Q^T ||_F ) / ',
     $   ' (eps * N)'
     	   WRITE(*,*) 
     	   WRITE(*,*) 
     $  'Test passes if both residuals are less then threshold'        
	   WRITE( NOUT, * )
	   WRITE( NOUT, FMT = 9995 )
	   WRITE( NOUT, FMT = 9994 )
      END IF
*
*     Loop over problem parameters.
*
      DO KTOP = 1, 1
      DO KBOT = N, N
      DO SOLVER = SLV_MAX, SLV_MIN, -1
*
*        Set INFO to zero for this run.
*
         INFO = 0
         NPROCS = NPROW*NPCOL
         TEMP_ICTXT = ICTXT
*
*        Count the number of rows and columns of current problem
*        for the current block sizes and grid properties.
*
         STAMP = MPI_WTIME()
         AROWS = NUMROC( N, NB, MYROW, 0, NPROW )
         ACOLS = NUMROC( N, NB, MYCOL, 0, NPCOL )
*
*        Set up matrix descriptors.
*
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ': Set up descriptors...'
         IF( BARR ) CALL BLACS_BARRIER(ICTXT, 'A')
         CALL DESCINIT( DESCA, N, N, NB, NB, MIN(ARSRC,NPROW-1),
     $        MIN(NPCOL-1,ACSRC), TEMP_ICTXT, MAX(1, AROWS), INFO )
         IF ( INFO.NE.0 ) THEN
            WRITE(*,*) "% DESCINIT DESCA failed, INFO =", INFO
            GO TO 999
         END IF
         CALL DESCINIT( DESCQ, N, N, NB, NB, MIN(ARSRC,NPROW-1),
     $        MIN(NPCOL-1,ACSRC), TEMP_ICTXT, MAX(1, AROWS), INFO )
         IF ( INFO.NE.0 ) THEN
            WRITE(*,*) "% DESCINIT DESCQ failed, INFO =", INFO
            GO TO 999
         END IF
         CALL DESCINIT( DESCVEC, N, 1, N, 1, MIN(ARSRC,NPROW-1),
     $        MIN(NPCOL-1,ACSRC), TEMP_ICTXT, N, INFO )
         IF ( INFO.NE.0 ) THEN
            WRITE(*,*) "% DESCINIT DESCVEC failed, INFO =", INFO
            GO TO 999
         END IF
*
*        Assign pointer for ScaLAPACK arrays - first set DP memory.
*
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ': Assign pointers...'
         IPA    = 1
         IPACPY = IPA + DESCA( LLD_ ) * ACOLS
         IPQ    = IPACPY + DESCA( LLD_ ) * ACOLS
         WR1    = IPQ + DESCQ( LLD_ ) * ACOLS
         WI1    = WR1 + N
         WR2    = WI1 + N
         WI2    = WR2 + N
         IPW1   = WI2 + N
         IPW2   = IPW1 + DESCA( LLD_ ) * ACOLS
         IF( DEBUG ) WRITE(*,*) '% (IPW2,DPALLOC):', IPW2, DPALLOC
*         PRINT*, '%', IPA, IPACPY, IPQ, WR1, WI1, WR2, WI2, IPW1, IPW2
         IF( IPW2+DESCA(LLD_)*ACOLS .GT. DPALLOC+1 ) THEN
            WRITE(*,*) '% Not enough DP memory!'
            GO TO 999
         END IF
*
*        Then set integer memory pointers.
*
         IPIW = 1
*
*        Generate testproblem.
*
         IF( BARR ) CALL BLACS_BARRIER(ICTXT, 'A')
         CALL PDLASET( 'All over', N, N, ZERO, ONE, MEM(IPQ), 1, 1,
     $        DESCQ )
         CALL PDMATGEN2( TEMP_ICTXT, 'Random', 'NoDiagDominant',
     $        N, N, NB, NB, MEM(IPA), DESCA( LLD_ ), 0, 0, 7, 0,
     $        AROWS, 0, ACOLS, MYROW, MYCOL, NPROW, NPCOL )
         IF( .NOT. COMPHESS ) THEN
            CALL PDLASET( 'Lower triangular', N-2, N-2, ZERO, ZERO,
     $           MEM(IPA), 3, 1, DESCA )
            CALL PDLASET( 'All over', N, N, ZERO, ONE, MEM(IPQ),
     $           1, 1, DESCQ )
            IF( KTOP.GT.1 )
     $           CALL PDLASET( 'Lower triangular', KTOP-1, KTOP-1,
     $           ZERO, ZERO, MEM(IPA), 2, 1, DESCQ )
            IF( KBOT.LT.N )
     $           CALL PDLASET( 'Lower triangular', N-KBOT, N-KBOT,
     $           ZERO, ZERO, MEM(IPA), KBOT+1, KBOT, DESCQ )
         END IF
*
*        Do balancing if general matrix.
*
         IF( BARR ) CALL BLACS_BARRIER(ICTXT, 'A')
         T_BA = MPI_WTIME()
         IF( COMPHESS .AND. BALANCE ) THEN
            IF( NPROCS.EQ.1 .AND. SOLVER.NE.2 .AND. UNI_LAPACK ) THEN
               IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == dgebal =='
               CALL DGEBAL( 'Both', N, MEM(IPA), DESCA(LLD_), ILO,
     $              IHI, SCALE, INFO )
               IF ( INFO.NE.0 ) THEN
                  WRITE(*,*) "% DGEBAL failed, INFO =", INFO
                  GO TO 999
               END IF
            ELSE
               IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == pdgebal =='
               CALL PDGEBAL( 'Both', N, MEM(IPA), DESCA, ILO, IHI,
     $              SCALE, INFO )
               IF ( INFO.NE.0 ) THEN
                  WRITE(*,*) "% PDGEBAL failed, INFO =", INFO
                  GO TO 999
               END IF
            END IF
         ELSEIF( COMPHESS ) THEN
            ILO = 1
            IHI = N
         ELSE
            ILO = KTOP
            IHI = KBOT
         END IF
         T_BA = MPI_WTIME() - T_BA
c         IF( TIMESTEPS.AND.IAM.EQ.0 ) WRITE(*,*)
c     $      ' %%% Balancing took in seconds:',T_BA
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ': ILO,IHI=',ILO,IHI
*
*        Make a copy of A.
*
         IF( BARR ) CALL BLACS_BARRIER(ICTXT, 'A')
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ': Copy matrix A'
         CALL PDLACPY( 'All', N, N, MEM(IPA), 1, 1, DESCA, MEM(IPACPY),
     $                 1, 1, DESCA )
*
*        Print matrices to screen in debugging mode.
*
         IF( PRN )
     $      CALL PDLAPRNT( N, N, MEM(IPACPY), 1, 1, DESCA, 0, 0,
     $           'A', NOUT, MEM(IPW1) )
         T_GEN = T_GEN + MPI_WTIME() - STAMP - T_BA
c         IF( TIMESTEPS.AND.IAM.EQ.0 ) WRITE(*,*)
c     $      ' %%% Generation took in seconds:',T_GEN
*
*        Only compute the Hessenberg form if necessary.
*
         T_HS = MPI_WTIME()
         IF( .NOT. COMPHESS ) GO TO 30
*
*        Reduce A to Hessenberg form.
*
         IF( BARR ) CALL BLACS_BARRIER(ICTXT, 'A')
         IF( DEBUG ) WRITE(*,*) '% #', IAM,
     $      ': Reduce to Hessenberg form...N=',N, ILO,IHI
*         PRINT*, '% PDGEHRD: IPW2,MEM(IPW2)', IPW2, MEM(IPW2)
         IF( NPROCS.EQ.1 .AND. SOLVER.NE.2 .AND. UNI_LAPACK ) THEN
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == dgehrd =='
            CALL DGEHRD( N, ILO, IHI, MEM(IPA), DESCA(LLD_),
     $           MEM(IPW1), MEM(IPW2), -1, INFO )
            IF (DPALLOC-IPW2.LT.MEM(IPW2)) THEN
               WRITE(*,*) "% Not enough memory for DGEHRD"
               GO TO 999
            END IF
            CALL DGEHRD( N, ILO, IHI, MEM(IPA), DESCA(LLD_),
     $           MEM(IPW1), MEM(IPW2), DPALLOC-IPW2, INFO )
            IF ( INFO.NE.0 ) THEN
               WRITE(*,*) "% DGEHRD failed, INFO =", INFO
               GO TO 999
            END IF
         ELSE
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == pdgehrd =='
            CALL PDGEHRD( N, ILO, IHI, MEM(IPA), 1, 1, DESCA, MEM(IPW1),
     $           MEM(IPW2), -1, INFO )
            IF (DPALLOC-IPW2.LT.MEM(IPW2)) THEN
               WRITE(*,*) "% Not enough memory for PDGEHRD"
               GO TO 999
            END IF
            CALL PDGEHRD( N, ILO, IHI, MEM(IPA), 1, 1, DESCA, MEM(IPW1),
     $           MEM(IPW2), DPALLOC-IPW2, INFO )
            IF ( INFO.NE.0 ) THEN
               WRITE(*,*) "% PDGEHRD failed, INFO =", INFO
               GO TO 999
            END IF
         END IF
*
*        Form Q explicitly.
*
         IF( BARR ) CALL BLACS_BARRIER(ICTXT, 'A')
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ':Form Q explicitly'
*         PRINT*, '% PDORMHR: IPW2,MEM(IPW2)', IPW2, MEM(IPW2)
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == pdormhr =='
         CALL PDORMHR( 'L', 'N', N, N, ILO, IHI, MEM(IPA), 1, 1,
     $        DESCA, MEM(IPW1), MEM(IPQ), 1, 1, DESCQ, MEM(IPW2),
     $        -1, INFO )
         IF (DPALLOC-IPW2.LT.MEM(IPW2)) THEN
            WRITE(*,*) "% Not enough memory for PDORMHR"
            GO TO 999
         END IF
         CALL PDORMHR( 'L', 'N', N, N, ILO, IHI, MEM(IPA), 1, 1,
     $        DESCA, MEM(IPW1), MEM(IPQ), 1, 1, DESCQ, MEM(IPW2),
     $        DPALLOC-IPW2, INFO )
         IF ( INFO.NE.0 ) THEN
            WRITE(*,*) "% PDORMHR failed, INFO =", INFO
            GO TO 999
         END IF
*
*        Extract the upper Hessenberg part of A.
*
         CALL PDLASET( 'Lower triangular', N-2, N-2, ZERO, ZERO,
     $        MEM(IPA), 3, 1, DESCA )
*
*        Print reduced matrix A in debugging mode.
*
         IF( PRN ) THEN
            CALL PDLAPRNT( N, N, MEM(IPA), 1, 1, DESCA, 0, 0, 'H', NOUT,
     $           MEM(IPW1) )
            CALL PDLAPRNT( N, N, MEM(IPQ), 1, 1, DESCQ, 0, 0, 'Q', NOUT,
     $           MEM(IPW1) )
         END IF
*
 30      CONTINUE
         T_HS = MPI_WTIME() - T_HS
c         IF( TIMESTEPS.AND.IAM.EQ.0 ) WRITE(*,*)
c     $      ' %%% Hessenberg took in seconds:',T_HS
*
*        Compute the real Schur form of the Hessenberg matrix A.
*
         IF( BARR ) CALL BLACS_BARRIER(ICTXT, 'A')
         T_QR = MPI_WTIME()
         IF( SOLVER.EQ.1 ) THEN
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == pdlaqr1 =='
*            PRINT*, '% PDLAQR1: IPW1,MEM(IPW1)', IPW1, MEM(IPW1)
            CALL PDLAQR1( .TRUE., .TRUE., N, ILO, IHI, MEM(IPA), DESCA,
     $           MEM(WR1), MEM(WI1), ILO, IHI, MEM(IPQ), DESCQ,
     $           MEM(IPW1), -1, IMEM, -1, INFO )
            IF (DPALLOC-IPW1.LT.MEM(IPW1)) THEN
               WRITE(*,*) "% Not enough DP memory for PDLAQR1"
               GO TO 999
            END IF
            IF (INTALLC.LT.IMEM(1)) THEN
               WRITE(*,*) "% Not enough INT memory for PDLAQR1"
               GO TO 999
            END IF
            CALL PDLAQR1( .TRUE., .TRUE., N, ILO, IHI, MEM(IPA), DESCA,
     $           MEM(WR1), MEM(WI1), ILO, IHI, MEM(IPQ), DESCQ,
     $           MEM(IPW1), DPALLOC-IPW1+1, IMEM, INTALLC, INFO )
            IF (INFO.NE.0) THEN
               WRITE(*,*) "% PDLAQR1: INFO =", INFO
            END IF
         ELSEIF( SOLVER.EQ.2 ) THEN
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': == pdhseqr =='
*            PRINT*, '% PDHSEQR: IPW1,MEM(IPW1)', IPW1, MEM(IPW1)
         IF( BARR ) CALL BLACS_BARRIER(ICTXT, 'A')
            CALL PDHSEQR( 'Schur', 'Vectors', N, ILO, IHI, MEM(IPA),
     $           DESCA, MEM(WR2), MEM(WI2), MEM(IPQ), DESCQ, MEM(IPW1),
     $           -1, IMEM, -1, INFO )
         IF( BARR ) CALL BLACS_BARRIER(ICTXT, 'A')
            IF (DPALLOC-IPW1.LT.MEM(IPW1)) THEN
               WRITE(*,*) "% Not enough DP memory for PDHSEQR"
               GO TO 999
            END IF
            IF (INTALLC.LT.IMEM(1)) THEN
               WRITE(*,*) "% Not enough INT memory for PDHSEQR"
               GO TO 999
            END IF
         IF( BARR ) CALL BLACS_BARRIER(ICTXT, 'A')
            CALL PDHSEQR( 'Schur', 'Vectors', N, ILO, IHI, MEM(IPA),
     $           DESCA, MEM(WR2), MEM(WI2), MEM(IPQ), DESCQ, MEM(IPW1),
     $           DPALLOC-IPW1+1, IMEM, INTALLC, INFO )
            IF (INFO.NE.0) THEN
               WRITE(*,*) "% PDHSEQR: INFO =", INFO
            END IF
         ELSE
             WRITE(*,*) '% ERROR: Illegal SOLVER number!'
             GO TO 999
         END IF
         T_QR = MPI_WTIME() - T_QR
c         IF( TIMESTEPS.AND.IAM.EQ.0 ) WRITE(*,*)
c     $      ' %%% QR-algorithm took in seconds:',T_QR
         T_SCH = T_SCH + T_QR + T_HS + T_BA
*         TOTIT = IMEM(1)
*         SWEEPS = IMEM(2)
*         TOTNS = IMEM(3)
         ITPEREIG = DBLE(TOTIT) / DBLE(N)
         SWPSPEIG = DBLE(SWEEPS) / DBLE(N)
         NSPEIG = DBLE(TOTNS) / DBLE(N)
*
*        Print reduced matrix A in debugging mode.
*
         IF( PRN ) THEN
            CALL PDLAPRNT( N, N, MEM(IPA), 1, 1, DESCA, 0, 0, 'T',
     $           NOUT, MEM(IPW1) )
            CALL PDLAPRNT( N, N, MEM(IPQ), 1, 1, DESCQ, 0, 0, 'Z',
     $           NOUT, MEM(IPW1) )
         END IF
*
*        Check that returned Schur form is really a quasi-triangular
*        matrix.
*
         HESS = 0
         DO I = 1, N-1
            IF( I.GT.1 ) THEN
               CALL PDELGET( 'All', '1-Tree', ELEM1, MEM(IPA), I, I-1,
     $              DESCA )
            ELSE
               ELEM1 = ZERO
            END IF
            CALL PDELGET( 'All', '1-Tree', ELEM2, MEM(IPA), I+1, I,
     $           DESCA )
            IF( I.LT.N-1 ) THEN
               CALL PDELGET( 'All', '1-Tree', ELEM3, MEM(IPA), I+2, I+1,
     $              DESCA )
            ELSE
               ELEM3 = ZERO
            END IF
            IF( ELEM2.NE.ZERO .AND. ABS(ELEM1)+ABS(ELEM2)+ABS(ELEM3).GT.
     $         ABS(ELEM2) ) HESS = HESS + 1
         END DO
*
*        Compute residual norms and other results:
*
*           1) RNORM = || T - Q'*A*Q ||_F / ||A||_F
*           2) ORTH  = MAX( || I - Q'*Q ||_F, || I - Q*Q' ||_F ) /
*                  (epsilon*N)
*
         STAMP = MPI_WTIME()
         IF( COMPRESI ) THEN
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': Compute residuals 1'
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': pdgemm 3'
            CALL PDGEMM( 'N', 'N', N, N, N, ONE, MEM(IPACPY), 1, 1,
     $           DESCA, MEM(IPQ), 1, 1, DESCQ, ZERO, MEM(IPW1), 1, 1,
     $           DESCA )
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': pdgemm 4'
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': N=',N
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': DESCA=',DESCA(1:DLEN_)
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': DESCQ=',DESCQ(1:DLEN_)
            CALL PDGEMM( 'T', 'N', N, N, N, -ONE, MEM(IPQ), 1, 1,
     $           DESCQ, MEM(IPW1), 1, 1, DESCA, ONE, MEM(IPA), 1, 1,
     $           DESCA )
            R1 = PDLANGE( 'Frobenius', N, N, MEM(IPA), 1, 1, DESCA,
     $           DPDUM )
            ANORM = PDLANGE( 'Frobenius', N, N, MEM(IPACPY), 1, 1,
     $           DESCA, DPDUM )
            IF( ANORM.GT.ZERO )THEN
               RNORM = R1 / (ANORM*EPS*SQRT(DBLE(N)))
            ELSE
               RNORM = R1
            END IF
         ELSE
            RNORM = 0.0D0
         END IF
*
         IF( COMPORTH ) THEN
            IF( DEBUG ) WRITE(*,*) '% #', IAM, ': Compute residuals 2'
            CALL PDLASET( 'All', N, N, ZERO, ONE, MEM(IPW1), 1, 1,
     $           DESCQ )
            CALL PDLACPY( 'All', N, N, MEM(IPQ), 1, 1, DESCQ, MEM(IPW2),
     $           1, 1, DESCQ )
            CALL PDGEMM( 'T', 'N', N, N, N, -ONE, MEM(IPQ), 1, 1, DESCQ,
     $           MEM(IPW2), 1, 1, DESCQ, ONE, MEM(IPW1), 1, 1, DESCQ )
            O1 = PDLANGE( 'Frobenius', N, N, MEM(IPW1), 1, 1, DESCQ,
     $           DPDUM )
            CALL PDLASET( 'All', N, N, ZERO, ONE, MEM(IPW1), 1, 1,
     $           DESCQ )
            CALL PDGEMM( 'N', 'T', N, N, N, -ONE, MEM(IPQ), 1, 1, DESCQ,
     $           MEM(IPW2), 1, 1, DESCQ, ONE, MEM(IPW1), 1, 1, DESCQ )
            O2 = PDLANGE( 'Frobenius', N, N, MEM(IPW1), 1, 1, DESCQ,
     $           DPDUM )
            ORTH = MAX(O1,O2) / (EPS*DBLE(N))
         ELSE
            ORTH = 0.0D0
         END IF
*
         T_RES = T_RES + MPI_WTIME() - STAMP
c         IF( TIMESTEPS.AND.IAM.EQ.0 ) WRITE(*,*)
c     $      ' %%% Residuals took in seconds:',T_RES
         TOTTIME = MPI_WTIME() - TOTTIME
c         IF( IAM.EQ.0 ) WRITE(*,*)
c     $      ' %%% Total execution time in seconds:', TOTTIME
*
*
*        Print results to screen.
*
	   IF( (ORTH.GT.THRESH).OR.(RNORM.GT.THRESH) ) THEN
	      PASSED = 'FAILED'
	   ELSE
	      PASSED = 'PASSED'
	   END IF
         IF( DEBUG ) WRITE(*,*) '% #', IAM, ': Print results...'
         IF( IAM.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9993 ) N, NB, NPROW, NPCOL, T_QR, PASSED
         END IF
         CALL BLACS_BARRIER( ICTXT, 'All' )
      END DO
      END DO
      END DO
 999  CONTINUE
*
*     Deallocate MEM and IMEM.
*
      DEALLOCATE( MEM, IMEM )
*
      CALL BLACS_GRIDEXIT( ICTXT )
*
 777  CONTINUE
*
      CALL BLACS_EXIT( 0 )
*
*     Format specifications.
*
 6666 FORMAT(A2,A3,A6,A4,A5,A6,A3,A3,A3,A9,A9,A9,A8,A8,A9,A9,A9,A9,A9,
     $       A9,A9,A9,A9,A9,A9,A5,A5,A8,A5,A5)
 7777 FORMAT(A2,I3,I6,I4,I5,I6,I3,I3,I3,F9.2,F9.2,F9.2,F8.2,F8.2,F9.2,
     $       F9.2,F9.2,F9.2,F9.2,F9.2,F9.2,F9.2,E9.2,E9.2,E9.2,I5,I5,
     $       F8.4,I5,I5,A2)
 9995 FORMAT( '    N  NB    P    Q  QR Time  CHECK' )
 9994 FORMAT( '----- --- ---- ---- -------- ------' )
 9993 FORMAT( I5, 1X, I3, 1X, I4, 1X, I4, 1X, F8.2, 1X, A6 )
          
*
      END
