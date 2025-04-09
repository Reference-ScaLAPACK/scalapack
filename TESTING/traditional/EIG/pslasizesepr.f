      SUBROUTINE PSLASIZESEPR( DESCA, IPREPAD, IPOSTPAD, SIZEMQRLEFT,
     $                         SIZEMQRRIGHT, SIZEQRF, SIZETMS, SIZEQTQ,
     $                         SIZECHK, SIZESYEVR, ISIZESYEVR,
     $                         SIZESUBTST, ISIZESUBTST, SIZETST,
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
      INTEGER            IPOSTPAD, IPREPAD, ISIZESUBTST, ISIZESYEVR,
     $                   ISIZETST, SIZECHK, SIZEMQRLEFT, SIZEMQRRIGHT,
     $                   SIZEQRF, SIZEQTQ, SIZESUBTST, SIZESYEVR,
     $                   SIZETMS, SIZETST
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
*
*  Purpose
*  =======
*
*  PSLASIZESEPR computes the amount of memory needed by
*  various SEPR test routines, as well as PSSYEVR itself.
*
*  Arguments
*  =========
*
*  DESCA        (global input) INTEGER array dimension ( DLEN_ )
*               Array descriptor for dense matrix.
*
*  SIZEMQRLEFT  LWORK for the 1st PSORMQR call in PSLAGSY
*
*  SIZEMQRRIGHT LWORK for the 2nd PSORMQR call in PSLAGSY
*
*  SIZEQRF      LWORK for PSGEQRF in PSLAGSY
*
*  SIZETMS      LWORK for PSLATMS
*
*  SIZEQTQ      LWORK for PSSEPQTQ
*
*  SIZECHK      LWORK for PSSEPCHK
*
*  SIZESYEVR    LWORK for PSSYEVR
*
*  ISIZESYEVR   LIWORK for PSSYEVR
*
*  SIZESUBTST   LWORK for PSSEPRSUBTST
*
*  ISIZESUBTST  LIWORK for PSSEPRSUBTST
*
*  SIZETST      LWORK for PSSEPRTST
*
*  ISIZETST     LIWORK for PSSEPRTST
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
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, ILCM, INDXG2P, NUMROC
      EXTERNAL           ICEIL, ILCM, INDXG2P, NUMROC
*
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
      SIZESYEVR = 1 + 5*N + MAX( 18*NN, NP0*MQ0+2*NB*NB ) +
     $            (2 + ICEIL( NEIG, NPROW*NPCOL ))*NN
      SIZESYEVR = MAX(3, SIZESYEVR)
*
      ISIZESYEVR = 12*NNP + 2*N
*
      SIZESUBTST = MAX( SIZETMS, SIZEQTQ, SIZECHK, SIZESYEVR ) +
     $             IPREPAD + IPOSTPAD
      ISIZESUBTST = ISIZESYEVR + IPREPAD + IPOSTPAD
*
*     Allow room for A, COPYA and Z and DIAG, WIN, WNEW, GAP, WORK
*
      SIZETST = 3*( LDA*NP+IPREPAD+IPOSTPAD ) +
     $          4*( N+IPREPAD+IPOSTPAD ) + SIZESUBTST
*
*     Allow room for IFAIL, ICLUSTR, and IWORK 
*     (only needed for PSSYEVX)
*
      ISIZETST = N + 2*NPROW*NPCOL + 2*( IPREPAD+IPOSTPAD ) +
     $           ISIZESUBTST
*
*
      RETURN
      END
