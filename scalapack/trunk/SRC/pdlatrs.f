      SUBROUTINE PDLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, IA,
     $                    JA, DESCA, X, IX, JX, DESCX, SCALE, CNORM,
     $                    WORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            IA, IX, JA, JX, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCX( * )
      DOUBLE PRECISION   A( * ), CNORM( * ),
     $                   X( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLATRS solves a triangular system. This routine in unfinished
*  at this time, but will be part of the next release.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*
*     .. Local Scalars ..
      INTEGER            ICTXT, IIX, IROFF, JJX, MYCOL, MYROW, NP,
     $                   NPCOL, NPROW, LDX, IXCOL, IXROW
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGEBR2D, DGEBS2D, INFOG2L,
     $                   PDTRSV
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     *****  NO SCALING ***** Call PDTRSV for all cases  *****
*
      SCALE = ONE
      CALL PDTRSV( UPLO, TRANS, DIAG, N, A, IA, JA, DESCA, X, IX, JX,
     $             DESCX, 1 )
*
      CALL INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL, IIX, JJX,
     $              IXROW, IXCOL )
      LDX = DESCX( LLD_ )
      IROFF = MOD( IX-1, DESCX(MB_) )
      NP = NUMROC( N+IROFF, DESCX( MB_ ), MYROW, IXROW, NPROW )
      IF( MYROW.EQ.IXROW )
     $   NP = NP - IROFF
      IF( MYCOL.EQ.IXCOL ) THEN
         CALL DGEBS2D( ICTXT, 'R', ' ', NP, 1, X( IIX+(JJX-1)*LDX ),
     $                 LDX )
      ELSE
         CALL DGEBR2D( ICTXT, 'R', ' ', NP, 1, X( IIX+(JJX-1)*LDX ),
     $                 LDX, MYROW, IXCOL )
      END IF
*
      RETURN
*
*     End of PDLATRS
*
      END
