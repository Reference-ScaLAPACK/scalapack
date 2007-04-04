      SUBROUTINE PCLAREAD( FILNAM, A, DESCA, IRREAD, ICREAD, WORK )
*
*  -- ScaLAPACK tools routine (version 1.8) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
* 
*     written by Antoine Petitet, August 1995 (petitet@cs.utk.edu)
*     adapted by Julie Langou, April 2007 (julie@cs.utk.edu)
* 
*     .. Scalar Arguments ..
      INTEGER            ICREAD, IRREAD
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      FILNAM
      INTEGER            DESCA( * )
      COMPLEX            A( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLAREAD reads from a file named FILNAM a matrix and distribute
*  it to the process grid.
*
*  Only the process of coordinates {IRREAD, ICREAD} read the file.
*
*  WORK must be of size >= MB_ = DESCA( MB_ ).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NIN
      PARAMETER          ( NIN = 11 )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            H, I, IB, ICTXT, ICURCOL, ICURROW, II, J, JB,
     $                   JJ, K, LDA, M, MYCOL, MYROW, N, NPCOL, NPROW
      REAL               REAL_PART, IMAG_PART
*    ..
*     .. Local Arrays ..
      INTEGER            IWORK( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, CGERV2D, CGESD2D,
     $                   IGEBS2D, IGEBR2D
*     ..
*     .. External Functions ..
      INTEGER            ICEIL
      EXTERNAL           ICEIL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( MYROW.EQ.IRREAD .AND. MYCOL.EQ.ICREAD ) THEN
         OPEN( NIN, FILE=FILNAM, STATUS='OLD' )
         READ( NIN, FMT = * ) ( IWORK( I ), I = 1, 2 )
         CALL IGEBS2D( ICTXT, 'All', ' ', 2, 1, IWORK, 2 )
      ELSE
         CALL IGEBR2D( ICTXT, 'All', ' ', 2, 1, IWORK, 2, IRREAD,
     $                 ICREAD )
      END IF
      M = IWORK( 1 )
      N = IWORK( 2 )
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( M.GT.DESCA( M_ ).OR. N.GT.DESCA( N_ ) ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            WRITE( *, FMT = * ) 'PCLAREAD: Matrix too big to fit in'
            WRITE( *, FMT = * ) 'Abort ...'
         END IF
         CALL BLACS_ABORT( ICTXT, 0 )
      END IF
*
      II = 1
      JJ = 1
      ICURROW = DESCA( RSRC_ )
      ICURCOL = DESCA( CSRC_ )
      LDA = DESCA( LLD_ )
*
*     Loop over column blocks
*
      DO 50 J = 1, N, DESCA( NB_ )
         JB = MIN(  DESCA( NB_ ), N-J+1 )
         DO 40 H = 0, JB-1
*
*           Loop over block of rows
*
            DO 30 I = 1, M, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), M-I+1 )
               IF( ICURROW.EQ.IRREAD .AND. ICURCOL.EQ.ICREAD ) THEN
                  IF( MYROW.EQ.IRREAD .AND. MYCOL.EQ.ICREAD ) THEN
                     DO 10 K = 0, IB-1
                   READ( NIN , FMT = *) REAL_PART, IMAG_PART 
                  A( II+K+(JJ+H-1)*LDA ) = CMPLX(REAL_PART, IMAG_PART)
   10                CONTINUE
                  END IF
               ELSE
                  IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                     CALL CGERV2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                             LDA, IRREAD, ICREAD )
                   ELSE IF( MYROW.EQ.IRREAD .AND. MYCOL.EQ.ICREAD ) THEN
                     DO 20 K = 1, IB
                        READ( NIN, FMT = * ) REAL_PART, IMAG_PART
                        WORK(K)=CMPLX(REAL_PART,IMAG_PART)
   20                CONTINUE
                     CALL CGESD2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                             ICURROW, ICURCOL )
                  END IF
               END IF
               IF( MYROW.EQ.ICURROW )
     $            II = II + IB
               ICURROW = MOD( ICURROW+1, NPROW )
   30       CONTINUE
*
            II = 1
            ICURROW = DESCA( RSRC_ )
   40    CONTINUE
*
         IF( MYCOL.EQ.ICURCOL )
     $      JJ = JJ + JB
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
   50 CONTINUE
*
      IF( MYROW.EQ.IRREAD .AND. MYCOL.EQ.ICREAD ) THEN
         CLOSE( NIN )
      END IF
*
      RETURN
*
*     End of PCLAREAD
*
      END
