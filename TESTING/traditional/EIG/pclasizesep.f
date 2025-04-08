      SUBROUTINE PCLASIZESEP( DESCA, IPREPAD, IPOSTPAD, SIZEMQRLEFT,
     $                        SIZEMQRRIGHT, SIZEQRF, SIZETMS, RSIZEQTQ,
     $                        RSIZECHK, SIZEHEEVX, RSIZEHEEVX,
     $                        ISIZEHEEVX, SIZEHEEVD, RSIZEHEEVD,
     $                        ISIZEHEEVD, SIZESUBTST, RSIZESUBTST,
     $                        ISIZESUBTST, SIZETST, RSIZETST, ISIZETST )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IPOSTPAD, IPREPAD, ISIZEHEEVD, ISIZEHEEVX,
     $                   ISIZESUBTST, ISIZETST, RSIZECHK, RSIZEHEEVD,
     $                   RSIZEHEEVX, RSIZEQTQ, RSIZESUBTST, RSIZETST,
     $                   SIZEHEEVD, SIZEHEEVX, SIZEMQRLEFT,
     $                   SIZEMQRRIGHT, SIZEQRF, SIZESUBTST, SIZETMS,
     $                   SIZETST
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLASIZESEP computes the amount of memory needed by
*    various SEP test routines, as well as HEEVX itself
*
*  Arguments
*  =========
*
*  DESCA        (global input) INTEGER array dimension ( DLEN_ )
*               Array descriptor as passed to PCHEEVX
*
*  SIZEMQRLEFT  LWORK for the 1st PCUNMQR call in PCLAGHE
*
*  SIZEMQRRIGHT LWORK for the 2nd PCUNMQR call in PCLAGHE
*
*  SIZEQRF      LWORK for PCGEQRF in PCLAGHE
*
*  SIZETMS      LWORK for PCLATMS
*
*  RSIZEQTQ      LWORK for PCSEPQTQ (nexer complex)
*
*  RSIZECHK      LWORK for PCSEPCHK
*
*  SIZEHEEVX    LWORK for PCHEEVX
*
*  RSIZEHEEVX   LRWORK for PCHEEVX
*
*  ISIZEHEEVX   LIWORK for PCHEEVX
*
*  SIZEHEEVD    LWORK for PCHEEVD
*
*  RSIZEHEEVD   LRWORK for PCHEEVD
*
*  ISIZEHEEVD   LIWORK for PCHEEVD
*
*  SIZESUBTST   LWORK for PCSUBTST
*
*  RSIZESUBTST  LRWORK for PCSUBTST
*
*  ISIZESUBTST  LIWORK for PCSUBTST
*
*  SIZETST      LWORK for PCTST
*
*  RSIZETST     LRWORK for PCTST
*
*  ISIZETST     LIWORK for PCTST
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            ANB, CSRC_A, IACOL, IAROW, ICOFFA, ICTXT,
     $                   IROFFA, LCM, LCMQ, LDA, MQ0, MYCOL, MYROW, N,
     $                   NB, NEIG, NHETRD_LWOPT, NN, NNP, NP, NP0,
     $                   NPCOL, NPROW, NPS, NQ, RSRC_A, SIZECHK,
     $                   SIZEQTQ, SQNPC
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, ILCM, INDXG2P, NUMROC, PJLAENV
      EXTERNAL           ICEIL, ILCM, INDXG2P, NUMROC, PJLAENV
*     ..
**     .. Executable Statements ..
*       This is just to keep ftnchek happy
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, REAL, SQRT
*     ..
*     .. Executable Statements ..
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
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
      SIZEQTQ = 0
      SIZECHK = 0
      RSIZEQTQ = 2 + MAX( DESCA( MB_ ), 2 )*( 2*NP0+MQ0 )
      RSIZECHK = NUMROC( N, DESCA( NB_ ), MYCOL, 0, NPCOL )
*
      NEIG = N
      NN = MAX( N, NB, 2 )
      NP0 = NUMROC( NN, NB, 0, 0, NPROW )
      MQ0 = NUMROC( MAX( NEIG, NB, 2 ), NB, 0, 0, NPCOL )
      SIZEHEEVX = N + ( NP0+MQ0+NB )*NB
      RSIZEHEEVX = 4*N + MAX( 5*NN, NP0*MQ0 ) +
     $             ICEIL( NEIG, NPROW*NPCOL )*NN
      NNP = MAX( N, NPROW*NPCOL+1, 4 )
      ISIZEHEEVX = 6*NNP
*
      ICTXT = DESCA( CTXT_ )
      ANB = PJLAENV( ICTXT, 3, 'PCHETTRD', 'L', 0, 0, 0, 0 )
      SQNPC = INT( SQRT( REAL( NPROW*NPCOL ) ) )
      NPS = MAX( NUMROC( N, 1, 0, 0, SQNPC ), 2*ANB )
      NHETRD_LWOPT = 2*( ANB+1 )*( 4*NPS+2 ) + ( NPS+2 )*NPS
*
      SIZEHEEVX = MAX( SIZEHEEVX, N+NHETRD_LWOPT )
*
      SIZEHEEVD = SIZEHEEVX
      RSIZEHEEVD = 7*N + 3*NP0*MQ0
      ISIZEHEEVD = 7*N + 8*NPCOL + 2
      SIZESUBTST = MAX( SIZETMS, SIZEQTQ, SIZECHK, SIZEHEEVX,
     $             SIZEHEEVD ) + IPREPAD + IPOSTPAD
      RSIZESUBTST = MAX( RSIZEHEEVX, RSIZEHEEVD, RSIZEQTQ, RSIZECHK ) +
     $              IPREPAD + IPOSTPAD
      ISIZESUBTST = MAX( ISIZEHEEVX, ISIZEHEEVD ) + IPREPAD + IPOSTPAD
*
*     Allow room for A, COPYA and Z and WORK
*
      SIZETST = 3*( LDA*NP+IPREPAD+IPOSTPAD ) + SIZESUBTST
*
*     Room for DIAG, WIN, WNEW, GAP and RWORK
*
      RSIZETST = 4*( N+IPREPAD+IPOSTPAD ) + RSIZESUBTST
*
*     Allow room for IFAIL, ICLUSTR, and IWORK (all in PCHEEVX)
*
      ISIZETST = N + 2*NPROW*NPCOL + 2*( IPREPAD+IPOSTPAD ) +
     $           ISIZESUBTST
*
      RETURN
      END
