      FUNCTION PDCHKRESI( N, A, IA, JA, DESCA, T, IT, JT, DESCT, 
     $                   Z, IZ, JZ, DESCZ, WORK, LWORK )
      IMPLICIT NONE
      INTEGER          N, IA, JA, IT, JT, IZ, JZ, LWORK
      DOUBLE PRECISION A(*), T(*), Z(*), WORK(LWORK)
      INTEGER          DESCA(*), DESCT(*), DESCZ(*)
      DOUBLE PRECISION PDCHKRESI
      LOGICAL          DEBUG
      INTEGER          BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                 LLD_, MB_, M_, NB_, N_, RSRC_
      DOUBLE PRECISION ONE, ZERO
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     ONE = 1.0D+00, ZERO = 0.0D+00,
     $                     DEBUG = .FALSE. )
      INTEGER          NUMROC
      DOUBLE PRECISION DPDUM, R1, ANORM
      DOUBLE PRECISION PDLANGE
      EXTERNAL         PDLANGE, PDGEMM, PDLACPY, PDLASET, NUMROC
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
         CALL PXERBLA( ICTXT, 'PDCHKRESI', -9 )
         RETURN
      END IF
*
      CALL PDLACPY( 'All', N, N, T, IT, JT, DESCT, WORK(IPW1),
     $     IT, JT, DESCT )
*
      CALL PDGEMM( 'N', 'N', N, N, N, 1D0, A, IA, JA, DESCA, Z, IZ, JZ, 
     $     DESCZ, 0D0, WORK(IPW2), IT, JT, DESCT )
      CALL PDGEMM( 'T', 'N', N, N, N, -1D0, Z, IZ, JZ, DESCZ, 
     $     WORK(IPW2), IT, JT, DESCT, 1D0, WORK(IPW1), IT, JT, DESCT )
*
      R1 = PDLANGE( 'Frobenius', N, N, WORK(IPW1), IT, JT, DESCT, 
     $     DPDUM )
      ANORM = PDLANGE( 'Frobenius', N, N, A, IA, JA, DESCA, DPDUM )
*
      IF(ANORM .LE. 1.0D-16) ANORM = ONE
      PDCHKRESI = R1 / ANORM
      END

