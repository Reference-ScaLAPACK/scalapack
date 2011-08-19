      SUBROUTINE PCLASIZESEPR( DESCA, IPREPAD, IPOSTPAD, SIZEMQRLEFT,
     $                         SIZEMQRRIGHT, SIZEQRF, SIZETMS, SIZEQTQ,
     $                         SIZECHK, SIZEHEEVR, RSIZEHEEVR, 
     $                         ISIZEHEEVR, SIZESUBTST, RSIZESUBTST, 
     $                         ISIZESUBTST, SIZETST, RSIZETST,
     $                         ISIZETST )
*
*  -- ScaLAPACK routine (@(MODE)version *TBA*) --
*     University of California, Berkeley and
*     University of Tennessee, Knoxville. 
*     October 21, 2006
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            IPOSTPAD, IPREPAD, ISIZEHEEVR, ISIZESUBTST,
     $                   ISIZETST, RSIZEHEEVR, RSIZESUBTST, RSIZETST,
     $                   SIZECHK, SIZEHEEVR, SIZEMQRLEFT, SIZEMQRRIGHT,
     $                   SIZEQRF, SIZEQTQ, SIZESUBTST, SIZETMS, SIZETST
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
*
*  Purpose
*  =======
*
*  PCLASIZESEPR computes the amount of memory needed by
*  various SEPR test routines, as well as PCHEEVR itself.
*
*  Arguments
*  =========
*
*  DESCA        (global input) INTEGER array dimension ( DLEN_ )
*               Array descriptor for dense matrix.
*
*  SIZEMQRLEFT  LWORK for the 1st PCUNMQR call in PCLAGHE
*
*  SIZEMQRRIGHT LWORK for the 2nd PCUNMQR call in PCLAGHE
*
*  SIZEQRF      LWORK for PCGEQRF in PCLAGHE
*
*  SIZETMS      LWORK for PCLATMS
*
*  SIZEQTQ      LWORK for PCSEPQTQ
*
*  SIZECHK      LWORK for PCSEPCHK
*
*  SIZEHEEVR    LWORK for PCHEEVR
*
*  RSIZEHEEVR   LRWORK for PCHEEVR
*
*  ISIZEHEEVR   LIWORK for PCHEEVR
*
*  SIZESUBTST   LWORK for PCSEPRSUBTST
*
*  RSIZESUBTST  LRWORK for PCSEPRSUBTST
*
*  ISIZESUBTST  LIWORK for PCSEPRSUBTST
*
*  SIZETST      LWORK for PCSEPRTST
*
*  RSIZETST     LRWORK for PCSEPRTST
*
*  ISIZETST     LIWORK for PCSEPRTST
*
*
*     .. Parameters ..
      INTEGER            CTXT_, M_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( 
     $                   CTXT_ = 2, M_ = 3, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            CSRC_A, IACOL, IAROW, ICOFFA, IROFFA, LCM,
     $                   LCMQ, LDA, MQ0, MYCOL, MYROW, N, NB, NEIG, NN,
     $                   NNP, NP, NP0, NPCOL, NPROW, NQ, RSRC_A
      INTEGER            ANB, ICTXT, NHETRD_LWOPT, NPS, SQNPC
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, ILCM, INDXG2P, NUMROC
      EXTERNAL           ICEIL, ILCM, INDXG2P, NUMROC
      INTEGER            PJLAENV
      EXTERNAL           PJLAENV
*
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      N = DESCA( M_ )
      NB = DESCA( MB_ )
      RSRC_A = DESCA( RSRC_ )
      CSRC_A = DESCA( CSRC_ )
*
      LDA = DESCA( LLD_ )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
      LCM = ILCM( NPROW, NPCOL )
      LCMQ = LCM / NPCOL
      IROFFA = 0
      ICOFFA = 0
      IAROW = INDXG2P( 1, NB, MYROW, RSRC_A, NPROW )
      IACOL = INDXG2P( 1, NB, MYCOL, CSRC_A, NPCOL )
      NP = NUMROC( N+IROFFA, NB, MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFFA, NB, MYCOL, IACOL, NPCOL )
      SIZEMQRLEFT = MAX( ( NB*( NB-1 ) ) / 2, ( NP+NQ )*NB ) + NB*NB
      SIZEMQRRIGHT = MAX( ( NB*( NB-1 ) ) / 2,
     $               ( NQ+MAX( NP+NUMROC( NUMROC( N+ICOFFA, NB, 0, 0,
     $               NPCOL ), NB, 0, 0, LCMQ ), NP ) )*NB ) + NB*NB
      SIZEQRF = NB*NP + NB*NQ + NB*NB
      SIZETMS = ( LDA+1 )*MAX( 1, NQ ) +
     $          MAX( SIZEMQRLEFT, SIZEMQRRIGHT, SIZEQRF )
*
      NP0 = NUMROC( N, DESCA( MB_ ), 0, 0, NPROW )
      MQ0 = NUMROC( N, DESCA( NB_ ), 0, 0, NPCOL )
      SIZEQTQ = 2 + MAX( DESCA( MB_ ), 2 )*( 2*NP0+MQ0 )
      SIZECHK = NUMROC( N, DESCA( NB_ ), MYCOL, 0, NPCOL )
*
      NEIG = N
      NN = MAX( N, NB, 2 ) + 1
      NP0 = NUMROC( NN, NB, 0, 0, NPROW )
      MQ0 = NUMROC( MAX( NEIG, NB, 2 ), NB, 0, 0, NPCOL )
      NNP = MAX( N, NPROW*NPCOL+1, 4 )
*
*
      SIZEHEEVR = 1+N + ( NP0+MQ0+NB )*NB
      SIZEHEEVR = MAX(3, SIZEHEEVR)
      RSIZEHEEVR = 1 + 5*N + MAX( 18*NN, NP0*MQ0+2*NB*NB ) +
     $            (2 + ICEIL( NEIG, NPROW*NPCOL ))*NN
      RSIZEHEEVR = MAX(3, RSIZEHEEVR)
*
      ISIZEHEEVR = 12*NNP + 2*N
*
      ICTXT = DESCA( CTXT_ )
      ANB = PJLAENV( ICTXT, 3, 'PCHETTRD', 'L', 0, 0, 0, 0 )
      SQNPC = INT( SQRT( DBLE( NPROW*NPCOL ) ) )
      NPS = MAX( NUMROC( N, 1, 0, 0, SQNPC ), 2*ANB )
      NHETRD_LWOPT = 2*( ANB+1 )*( 4*NPS+2 ) + ( NPS+2 )*NPS
      SIZEHEEVR = MAX( SIZEHEEVR, N + NHETRD_LWOPT )
*
      SIZESUBTST = MAX( SIZETMS,  SIZEHEEVR ) +
     $             IPREPAD + IPOSTPAD
      RSIZESUBTST = MAX( SIZEQTQ, SIZECHK, RSIZEHEEVR ) +
     $             IPREPAD + IPOSTPAD
      ISIZESUBTST = ISIZEHEEVR + IPREPAD + IPOSTPAD
*
*     Allow room for A, COPYA, Z, WORK
*
      SIZETST = 3*( LDA*NP+IPREPAD+IPOSTPAD ) + SIZESUBTST
*
*     Allow room for DIAG, WIN, WNEW, GAP, RWORK
*
      RSIZETST = 4*( N+IPREPAD+IPOSTPAD ) + RSIZESUBTST
*
*     Allow room for IFAIL, ICLUSTR, and IWORK 
*     (only needed for PCHEEVX)
*
      ISIZETST = N + 2*NPROW*NPCOL + 2*( IPREPAD+IPOSTPAD ) +
     $           ISIZESUBTST
*
*
      RETURN
      END
