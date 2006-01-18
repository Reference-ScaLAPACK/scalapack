      SUBROUTINE PCLAWIL( II, JJ, M, A, DESCA, H44, H33, H43H34, V )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     July 31, 2001
*
*     .. Scalar Arguments ..
      INTEGER            II, JJ, M
      COMPLEX            H33, H43H34, H44
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX            A( * ), V( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLAWIL gets the transform given by H44,H33, & H43H34 into V
*     starting at row M.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  II      (global input) INTEGER
*          Row owner of H(M+2,M+2)
*
*  JJ      (global input) INTEGER
*          Column owner of H(M+2,M+2)
*
*  M       (global input) INTEGER
*          On entry, this is where the transform starts (row M.)
*          Unchanged on exit.
*
*  A       (global input) COMPLEX array, dimension
*          (DESCA(LLD_),*)
*          On entry, the Hessenberg matrix.
*          Unchanged on exit.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*          Unchanged on exit.
*
*  H44
*  H33
*  H43H34  (global input) COMPLEX
*          These three values are for the double shift QR iteration.
*          Unchanged on exit.
*
*  V       (global output) COMPLEX array of size 3.
*          Contains the transform on ouput.
*
*  Further Details
*  ===============
*
*  Implemented by:  M. Fahey, May 28, 1999
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            CONTXT, DOWN, HBL, ICOL, IROW, JSRC, LDA, LEFT,
     $                   MODKM1, MYCOL, MYROW, NPCOL, NPROW, NUM, RIGHT,
     $                   RSRC, UP
      REAL               S
      COMPLEX            CDUM, H11, H12, H21, H22, H33S, H44S, V1, V2,
     $                   V3
*     ..
*     .. Local Arrays ..
      COMPLEX            BUF( 4 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, CGERV2D, CGESD2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL, AIMAG, MOD
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
      HBL = DESCA( MB_ )
      CONTXT = DESCA( CTXT_ )
      LDA = DESCA( LLD_ )
      CALL BLACS_GRIDINFO( CONTXT, NPROW, NPCOL, MYROW, MYCOL )
      LEFT = MOD( MYCOL+NPCOL-1, NPCOL )
      RIGHT = MOD( MYCOL+1, NPCOL )
      UP = MOD( MYROW+NPROW-1, NPROW )
      DOWN = MOD( MYROW+1, NPROW )
      NUM = NPROW*NPCOL
*
*     On node (II,JJ) collect all DIA,SUP,SUB info from M, M+1
*
      MODKM1 = MOD( M+1, HBL )
      IF( MODKM1.EQ.0 ) THEN
         IF( ( MYROW.EQ.II ) .AND. ( RIGHT.EQ.JJ ) .AND.
     $       ( NPCOL.GT.1 ) ) THEN
            CALL INFOG2L( M+2, M+1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IROW, ICOL, RSRC, JSRC )
            BUF( 1 ) = A( ( ICOL-1 )*LDA+IROW )
            CALL CGESD2D( CONTXT, 1, 1, BUF, 1, II, JJ )
         END IF
         IF( ( DOWN.EQ.II ) .AND. ( RIGHT.EQ.JJ ) .AND. ( NUM.GT.1 ) )
     $        THEN
            CALL INFOG2L( M, M, DESCA, NPROW, NPCOL, MYROW, MYCOL, IROW,
     $                    ICOL, RSRC, JSRC )
            BUF( 1 ) = A( ( ICOL-1 )*LDA+IROW )
            BUF( 2 ) = A( ( ICOL-1 )*LDA+IROW+1 )
            BUF( 3 ) = A( ICOL*LDA+IROW )
            BUF( 4 ) = A( ICOL*LDA+IROW+1 )
            CALL CGESD2D( CONTXT, 4, 1, BUF, 4, II, JJ )
         END IF
         IF( ( MYROW.EQ.II ) .AND. ( MYCOL.EQ.JJ ) ) THEN
            CALL INFOG2L( M+2, M+2, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IROW, ICOL, RSRC, JSRC )
            IF( NPCOL.GT.1 ) THEN
               CALL CGERV2D( CONTXT, 1, 1, V3, 1, MYROW, LEFT )
            ELSE
               V3 = A( ( ICOL-2 )*LDA+IROW )
            END IF
            IF( NUM.GT.1 ) THEN
               CALL CGERV2D( CONTXT, 4, 1, BUF, 4, UP, LEFT )
               H11 = BUF( 1 )
               H21 = BUF( 2 )
               H12 = BUF( 3 )
               H22 = BUF( 4 )
            ELSE
               H11 = A( ( ICOL-3 )*LDA+IROW-2 )
               H21 = A( ( ICOL-3 )*LDA+IROW-1 )
               H12 = A( ( ICOL-2 )*LDA+IROW-2 )
               H22 = A( ( ICOL-2 )*LDA+IROW-1 )
            END IF
         END IF
      END IF
      IF( MODKM1.EQ.1 ) THEN
         IF( ( DOWN.EQ.II ) .AND. ( RIGHT.EQ.JJ ) .AND. ( NUM.GT.1 ) )
     $        THEN
            CALL INFOG2L( M, M, DESCA, NPROW, NPCOL, MYROW, MYCOL, IROW,
     $                    ICOL, RSRC, JSRC )
            CALL CGESD2D( CONTXT, 1, 1, A( ( ICOL-1 )*LDA+IROW ), 1, II,
     $                    JJ )
         END IF
         IF( ( DOWN.EQ.II ) .AND. ( MYCOL.EQ.JJ ) .AND. ( NPROW.GT.1 ) )
     $        THEN
            CALL INFOG2L( M, M+1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IROW, ICOL, RSRC, JSRC )
            CALL CGESD2D( CONTXT, 1, 1, A( ( ICOL-1 )*LDA+IROW ), 1, II,
     $                    JJ )
         END IF
         IF( ( MYROW.EQ.II ) .AND. ( RIGHT.EQ.JJ ) .AND.
     $       ( NPCOL.GT.1 ) ) THEN
            CALL INFOG2L( M+1, M, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IROW, ICOL, RSRC, JSRC )
            CALL CGESD2D( CONTXT, 1, 1, A( ( ICOL-1 )*LDA+IROW ), 1, II,
     $                    JJ )
         END IF
         IF( ( MYROW.EQ.II ) .AND. ( MYCOL.EQ.JJ ) ) THEN
            CALL INFOG2L( M+2, M+2, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IROW, ICOL, RSRC, JSRC )
            IF( NUM.GT.1 ) THEN
               CALL CGERV2D( CONTXT, 1, 1, H11, 1, UP, LEFT )
            ELSE
               H11 = A( ( ICOL-3 )*LDA+IROW-2 )
            END IF
            IF( NPROW.GT.1 ) THEN
               CALL CGERV2D( CONTXT, 1, 1, H12, 1, UP, MYCOL )
            ELSE
               H12 = A( ( ICOL-2 )*LDA+IROW-2 )
            END IF
            IF( NPCOL.GT.1 ) THEN
               CALL CGERV2D( CONTXT, 1, 1, H21, 1, MYROW, LEFT )
            ELSE
               H21 = A( ( ICOL-3 )*LDA+IROW-1 )
            END IF
            H22 = A( ( ICOL-2 )*LDA+IROW-1 )
            V3 = A( ( ICOL-2 )*LDA+IROW )
         END IF
      END IF
      IF( ( MYROW.NE.II ) .OR. ( MYCOL.NE.JJ ) )
     $   RETURN
*
      IF( MODKM1.GT.1 ) THEN
         CALL INFOG2L( M+2, M+2, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 IROW, ICOL, RSRC, JSRC )
         H11 = A( ( ICOL-3 )*LDA+IROW-2 )
         H21 = A( ( ICOL-3 )*LDA+IROW-1 )
         H12 = A( ( ICOL-2 )*LDA+IROW-2 )
         H22 = A( ( ICOL-2 )*LDA+IROW-1 )
         V3 = A( ( ICOL-2 )*LDA+IROW )
      END IF
*
      H44S = H44 - H11
      H33S = H33 - H11
      V1 = ( H33S*H44S-H43H34 ) / H21 + H12
      V2 = H22 - H11 - H33S - H44S
      S = CABS1( V1 ) + CABS1( V2 ) + CABS1( V3 )
      V1 = V1 / S
      V2 = V2 / S
      V3 = V3 / S
      V( 1 ) = V1
      V( 2 ) = V2
      V( 3 ) = V3
*
      RETURN
*
*     End of PCLAWIL
*
      END
