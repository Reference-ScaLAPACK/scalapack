      FUNCTION PSCHKRESI( N, A, IA, JA, DESCA, T, IT, JT, DESCT, 
     $                   Z, IZ, JZ, DESCZ, WORK, LWORK )
      IMPLICIT NONE
      INTEGER          N, IA, JA, IT, JT, IZ, JZ, LWORK
      REAL             A(*), T(*), Z(*), WORK(LWORK)
      INTEGER          DESCA(*), DESCT(*), DESCZ(*)
      REAL             PSCHKRESI
      LOGICAL          DEBUG
      INTEGER          BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                 LLD_, MB_, M_, NB_, N_, RSRC_
      REAL             ONE, ZERO
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     ONE = 1.0, ZERO = 0.0,
     $                     DEBUG = .FALSE. )
      INTEGER          NUMROC
      REAL             DPDUM, R1, ANORM
      REAL             PSLANGE
      EXTERNAL         PSLANGE, PSGEMM, PSLACPY, PSLASET, NUMROC
      INTEGER          IPW1, IPW2, IPW3, ACOLS, TCOLS, ZCOLS, MYROW, 
     $                 MYCOL, NPROW, NPCOL, ICTXT
*
      ICTXT = DESCZ( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      TCOLS = NUMROC( DESCT(N_), DESCT(NB_), MYCOL, DESCT(CSRC_), 
     $     NPCOL )
      IPW1 = 1
      IPW2 = IPW1 + DESCT(LLD_)*TCOLS
      IPW3 = IPW2 + DESCT(LLD_)*TCOLS
      IF( IPW3-1.GT.LWORK ) THEN
         CALL PXERBLA( ICTXT, 'PSCHKRESI', -9 )
         RETURN
      END IF
*
      CALL PSLACPY( 'All', N, N, T, IT, JT, DESCT, WORK(IPW1),
     $     IT, JT, DESCT )
*
      CALL PSGEMM( 'N', 'N', N, N, N, 1D0, A, IA, JA, DESCA, Z, IZ, JZ, 
     $     DESCZ, 0D0, WORK(IPW2), IT, JT, DESCT )
      CALL PSGEMM( 'T', 'N', N, N, N, -1D0, Z, IZ, JZ, DESCZ, 
     $     WORK(IPW2), IT, JT, DESCT, 1D0, WORK(IPW1), IT, JT, DESCT )
*
      R1 = PSLANGE( 'Frobenius', N, N, WORK(IPW1), IT, JT, DESCT, 
     $     DPDUM )
      ANORM = PSLANGE( 'Frobenius', N, N, A, IA, JA, DESCA, DPDUM )
*
      IF(ANORM.LT.1.0e-7) ANORM = ONE
      PSCHKRESI = R1 / ANORM
      END

