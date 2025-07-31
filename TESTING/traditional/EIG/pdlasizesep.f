*
*
      SUBROUTINE PDLASIZESEP( DESCA, IPREPAD, IPOSTPAD, SIZEMQRLEFT,
     $                        SIZEMQRRIGHT, SIZEQRF, SIZETMS, SIZEQTQ,
     $                        SIZECHK, SIZESYEVX, ISIZESYEVX,
     $                        SIZESUBTST, ISIZESUBTST, SIZETST,
     $                        ISIZETST )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IPOSTPAD, IPREPAD, ISIZESUBTST, ISIZESYEVX,
     $                   ISIZETST, SIZECHK, SIZEMQRLEFT, SIZEMQRRIGHT,
     $                   SIZEQRF, SIZEQTQ, SIZESUBTST, SIZESYEVX,
     $                   SIZETMS, SIZETST
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLASIZESEP computes the amount of memory needed by
*    various SEP test routines, as well as SYEVX itself
*
*  Arguments
*  =========
*
*  DESCA        (global input) INTEGER array dimension ( DLEN_ )
*               Array descriptor as passed to PDSYEVX
*
*  SIZEMQRLEFT  LWORK for the 1st PDORMQR call in PDLAGSY
*
*  SIZEMQRRIGHT LWORK for the 2nd PDORMQR call in PDLAGSY
*
*  SIZEQRF      LWORK for PDGEQRF in PDLAGSY
*
*  SIZETMS      LWORK for PDLATMS
*
*  SIZEQTQ      LWORK for PDSEPQTQ (nexer complex)
*
*  SIZECHK      LWORK for PDSEPCHK
*
*  SIZESYEVX    LWORK for PDSYEVX
*
*  ISIZESYEVX   LIWORK for PDSYEVX
*
*  SIZESUBTST   LWORK for PDSUBTST
*
*  ISIZESUBTST  LIWORK for PDSUBTST
*
*  SIZETST      LWORK for PDTST
*
*  ISIZETST     LIWORK for PDTST
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            CSRC_A, IACOL, IAROW, ICOFFA, IROFFA, LCM,
     $                   LCMQ, LDA, MQ0, MYCOL, MYROW, N, NB, NEIG, NN,
     $                   NNP, NP, NP0, NPCOL, NPROW, NQ, RSRC_A
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, ILCM, INDXG2P, NUMROC
      EXTERNAL           ICEIL, ILCM, INDXG2P, NUMROC
*     ..
**     .. Executable Statements ..
*       This is just to keep ftnchek happy
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
      SIZEQTQ = 2 + MAX( DESCA( MB_ ), 2 )*( 2*NP0+MQ0 )
      SIZECHK = NUMROC( N, DESCA( NB_ ), MYCOL, 0, NPCOL )
*
      NEIG = N
      NN = MAX( N, NB, 2 )
      NP0 = NUMROC( NN, NB, 0, 0, NPROW )
      MQ0 = NUMROC( MAX( NEIG, NB, 2 ), NB, 0, 0, NPCOL )
      SIZESYEVX = 5*N + MAX( 5*NN, NP0*MQ0+2*NB*NB ) +
     $            ICEIL( NEIG, NPROW*NPCOL )*NN
      NNP = MAX( N, NPROW*NPCOL+1, 4 )
      ISIZESYEVX = 6*NNP
*
      SIZESUBTST = MAX( SIZETMS, SIZEQTQ, SIZECHK, SIZESYEVX ) +
     $             IPREPAD + IPOSTPAD
      ISIZESUBTST = ISIZESYEVX + IPREPAD + IPOSTPAD
*
*
*     Allow room for A, COPYA and Z and DIAG, WIN, WNEW, GAP, WORK
*
      SIZETST = 3*( LDA*NP+IPREPAD+IPOSTPAD ) +
     $          4*( N+IPREPAD+IPOSTPAD ) + SIZESUBTST
*
*     Allow room for IFAIL, ICLUSTR, and IWORK (all in PDSYEVX)
*
      ISIZETST = N + 2*NPROW*NPCOL + 2*( IPREPAD+IPOSTPAD ) +
     $           ISIZESUBTST
*
      RETURN
      END
