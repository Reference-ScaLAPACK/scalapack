*
*
      SUBROUTINE PSLASIZESYEV( JOBZ, N, DESCA, MINSIZE )
*
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ
      INTEGER            MINSIZE, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLASIZESYEV computes the amount of memory needed by PSSYEV
*  to calculate:
*    1)  Eigenvectors and eigenvalues if JOBZ = 'V'
*    2)  Eigenvalues only if JOBZ = 'N'
*
*  Arguments
*  =========
*
*     NP = the number of rows local to a given process.
*     NQ = the number of columns local to a given process.
*
*  JOBZ     (global input) CHARACTER*1
*           Specifies whether or not to compute the eigenvectors:
*           = 'N':  Compute eigenvalues only.
*           = 'V':  Compute eigenvalues and eigenvectors.
*
*  N        (global input) INTEGER
*           Size of the matrix to be tested.  (global size)
*
*  DESCA    (global input) INTEGER array dimension ( DLEN_ )
*
*  MINSIZE  (global output) INTEGER
*           Workspace required for PSSYEV to:
*           1)  Eigenvectors and eigenvalues if JOBZ = 'V'
*           2)  Eigenvalues only if JOBZ = 'N'
*
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      LOGICAL            WANTZ
      INTEGER            CONTEXTC, CSRC_A, IACOL, IAROW, ICOFFA, IROFFA, 
     $                   LCM, LCMQ, LDC, MQ0, MYCOL, MYPCOLC, MYPROWC, 
     $                   MYROW, NB, NN, NP, NP0, NPCOL, NPCOLC, NPROCS, 
     $                   NPROW, NPROWC, NQ, NRC, QRMEM, RSRC_A, 
     $                   SIZEMQRLEFT, SIZEMQRRIGHT
*     ..
*     .. External Functions ..
*
*
      LOGICAL            LSAME
      INTEGER            ILCM, INDXG2P, NUMROC, SL_GRIDRESHAPE
      EXTERNAL           ILCM, INDXG2P, LSAME, NUMROC, SL_GRIDRESHAPE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, BLACS_GRIDEXIT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      NB = DESCA( MB_ )
      N = DESCA( M_ )
      RSRC_A = DESCA( RSRC_ )
      CSRC_A = DESCA( CSRC_ )
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
      NN = MAX( N, NB, 2 )
      NP0 = NUMROC( NN, NB, 0, 0, NPROW )
      MQ0 = NUMROC( NN, NB, 0, 0, NPCOL )
      NPROCS = NPROW*NPCOL
      WANTZ = LSAME( JOBZ, 'V' )
      LDC = 0
*
*     Create the new context that is used in PSSYEV
*
      IF( WANTZ ) THEN
         CONTEXTC = SL_GRIDRESHAPE( DESCA( CTXT_ ), 0, 1, 1, NPROCS, 1 )
         CALL BLACS_GRIDINFO( CONTEXTC, NPROWC, NPCOLC, MYPROWC,
     $                        MYPCOLC )
         NRC = NUMROC( N, NB, MYPROWC, 0, NPROCS)
         LDC = MAX( 1, NRC )
         CALL BLACS_GRIDEXIT( CONTEXTC )
      END IF

*
*     Compute the total amount of space needed
*
      IF( WANTZ ) THEN
        QRMEM = 5*N + MAX(  2*NP0 +MQ0 + NB*NN, 2*NN-2 ) + N*LDC 
        MINSIZE = MAX ( SIZEMQRLEFT, SIZEMQRRIGHT, QRMEM )
      ELSE
         MINSIZE = 5*N + 2*NP0 +MQ0 + NB*NN 
      END IF
*
      RETURN
*
*     End of PSLASIZESYEV
*
      END
