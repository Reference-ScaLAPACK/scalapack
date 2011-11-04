      SUBROUTINE DESC_CONVERT( DESC_IN, DESC_OUT, INFO )
*
*
*     .. Array Arguments ..
      INTEGER DESC_IN( * ), DESC_OUT( * ), INFO
*     ..
*
*  Purpose
*  =======
*
*  Converts descriptors from one type to another if they are compatible.
*
*  Supports *ONLY* an output descriptor type of 1D_horizontal (type
*     number 501) or 1D_vertical (number 502).
*  Supports only one-dimensional 1xP input grids if descriptor_in is 2D.
*
*  Arguments
*  =========
*
*  DESC_IN: (input) input descriptor
*
*  DESC_OUT: (output) output descriptor (required to be 1D_horizontal
*            in this release).
*
*  INFO: (output) return code
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*
*     .. Local Scalars ..
      INTEGER         DESC_TYPE, DESC_TYPE_IN, ICTXT
      INTEGER         CSRC, RSRC, MB, NB, LLDA
      INTEGER         M, N, NPROW, NPCOL, IDUM1, IDUM2
*
*     .. External routines ..
*     EXTERNAL        BLACS_GRIDINFO
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
      DESC_TYPE_IN = DESC_IN( 1 )
*
*     .. Initialize Variables ..
*
	RSRC = 0
	NB = 0
	N = 0
	MB = 0
	M = 0
	LLDA = 0
	CSRC = 0
*	
      IF( DESC_TYPE_IN .EQ. BLOCK_CYCLIC_2D ) THEN
         ICTXT = DESC_IN( CTXT_ )
         RSRC = DESC_IN( RSRC_ )
         CSRC = DESC_IN( CSRC_ )
         MB = DESC_IN( MB_ )
         NB = DESC_IN( NB_ )
         LLDA = DESC_IN( LLD_ )
         M = DESC_IN( M_ )
         N = DESC_IN( N_ )
         CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, IDUM1, IDUM2 )
      ELSEIF ( DESC_TYPE_IN .EQ. 502 ) THEN
         ICTXT = DESC_IN( 2 )
         RSRC = DESC_IN( 5 )
         CSRC = 1
         MB = DESC_IN( 4 )
         NB = 1
         LLDA = DESC_IN( 6 )
         M = DESC_IN( 3 )
         N = 1
         NPROW = 0
         NPCOL = 1
      ELSEIF ( DESC_TYPE_IN .EQ. 501 ) THEN
         ICTXT = DESC_IN( 2 )
         RSRC = 1
         CSRC = DESC_IN( 5 )
         MB = 1
         NB = DESC_IN( 4 )
         LLDA = DESC_IN( 6 )
         M = 1
         N = DESC_IN( 3 )
         NPROW = 1
         NPCOL = 0
      ENDIF
*
*
      DESC_TYPE = DESC_OUT( 1 )
*
      IF( DESC_TYPE .EQ. 501 ) THEN
         IF( NPROW .NE. 1 )THEN
            INFO = -1
            RETURN
         ENDIF
         DESC_OUT( 2 ) = ICTXT
         DESC_OUT( 5 ) = CSRC
         DESC_OUT( 4 ) = NB
         DESC_OUT( 6 ) = LLDA
         DESC_OUT( 3 ) = N
      ELSEIF( DESC_TYPE .EQ. 502 ) THEN
         IF( NPCOL .NE. 1 )THEN
            INFO = -1
            RETURN
         ENDIF
         DESC_OUT( 2 ) = ICTXT
         DESC_OUT( 5 ) = RSRC
         DESC_OUT( 4 ) = MB
         DESC_OUT( 6 ) = LLDA
         DESC_OUT( 3 ) = M
      ENDIF
*
      RETURN
*
*     End of DESC_CONVERT
*
      END
