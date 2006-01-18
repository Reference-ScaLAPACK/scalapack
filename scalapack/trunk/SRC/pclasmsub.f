      SUBROUTINE PCLASMSUB( A, DESCA, I, L, K, SMLNUM, BUF, LWORK )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     July 31, 2001
*
*     .. Scalar Arguments ..
      INTEGER            I, K, L, LWORK
      REAL               SMLNUM
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX            A( * ), BUF( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLASMSUB looks for a small subdiagonal element from the bottom
*     of the matrix that it can safely set to zero.
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
*  A       (global input) COMPLEX array, dimension (DESCA(LLD_),*)
*          On entry, the Hessenberg matrix whose tridiagonal part is
*          being scanned.
*          Unchanged on exit.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  I       (global input) INTEGER
*          The global location of the bottom of the unreduced
*          submatrix of A.
*          Unchanged on exit.
*
*  L       (global input) INTEGER
*          The global location of the top of the unreduced submatrix
*          of A.
*          Unchanged on exit.
*
*  K       (global output) INTEGER
*          On exit, this yields the bottom portion of the unreduced
*          submatrix.  This will satisfy: L <= M  <= I-1.
*
*  SMLNUM  (global input) REAL
*          On entry, a "small number" for the given matrix.
*          Unchanged on exit.
*
*  BUF     (local output) COMPLEX array of size LWORK.
*
*  LWORK   (global input) INTEGER
*          On exit, LWORK is the size of the work buffer.
*          This must be at least 2*Ceil( Ceil( (I-L)/HBL ) /
*                                        LCM(NPROW,NPCOL) )
*          Here LCM is least common multiple, and NPROWxNPCOL is the
*          logical grid size.
*
*   Notes:
*
*     This routine does a global maximum and must be called by all
*     processes.
*
*     This code is basically a parallelization of the following snip
*     of LAPACK code from CLAHQR:
*
*        Look for a single small subdiagonal element.
*
*        DO 20 K = I, L + 1, -1
*           TST1 = CABS1( H( K-1, K-1 ) ) + CABS1( H( K, K ) )
*           IF( TST1.EQ.ZERO )
*    $         TST1 = CLANHS( '1', I-L+1, H( L, L ), LDH, WORK )
*           IF( CABS1( H( K, K-1 ) ).LE.MAX( ULP*TST1, SMLNUM ) )
*    $         GO TO 30
*  20    CONTINUE
*  30    CONTINUE
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
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            CONTXT, DOWN, HBL, IBUF1, IBUF2, ICOL1, ICOL2,
     $                   II, III, IRCV1, IRCV2, IROW1, IROW2, ISRC,
     $                   ISTR1, ISTR2, ITMP1, ITMP2, JJ, JJJ, JSRC, LDA,
     $                   LEFT, MODKM1, MYCOL, MYROW, NPCOL, NPROW, NUM,
     $                   RIGHT, UP
      REAL               TST1, ULP
      COMPLEX            CDUM, H10, H11, H22
*     ..
*     .. External Functions ..
      INTEGER            ILCM, NUMROC
      REAL               PSLAMCH
      EXTERNAL           ILCM, NUMROC, PSLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGAMX2D, INFOG1L, INFOG2L,
     $                   CGERV2D, CGESD2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL, AIMAG, MAX, MOD
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
      ULP = PSLAMCH( CONTXT, 'PRECISION' )
      CALL BLACS_GRIDINFO( CONTXT, NPROW, NPCOL, MYROW, MYCOL )
      LEFT = MOD( MYCOL+NPCOL-1, NPCOL )
      RIGHT = MOD( MYCOL+1, NPCOL )
      UP = MOD( MYROW+NPROW-1, NPROW )
      DOWN = MOD( MYROW+1, NPROW )
      NUM = NPROW*NPCOL
*
*     BUFFER1 STARTS AT BUF(ISTR1+1) AND WILL CONTAINS IBUF1 ELEMENTS
*     BUFFER2 STARTS AT BUF(ISTR2+1) AND WILL CONTAINS IBUF2 ELEMENTS
*
      ISTR1 = 0
      ISTR2 = ( ( I-L ) / HBL )
      IF( ISTR2*HBL.LT.( I-L ) )
     $   ISTR2 = ISTR2 + 1
      II = ISTR2 / ILCM( NPROW, NPCOL )
      IF( II*ILCM( NPROW, NPCOL ).LT.ISTR2 ) THEN
         ISTR2 = II + 1
      ELSE
         ISTR2 = II
      END IF
      IF( LWORK.LT.2*ISTR2 ) THEN
*
*        Error!
*
         RETURN
      END IF
      CALL INFOG2L( I, I, DESCA, NPROW, NPCOL, MYROW, MYCOL, IROW1,
     $              ICOL1, II, JJ )
      MODKM1 = MOD( I-1+HBL, HBL )
*
*     COPY OUR RELEVANT PIECES OF TRIADIAGONAL THAT WE OWE INTO
*     2 BUFFERS TO SEND TO WHOMEVER OWNS H(K,K) AS K MOVES DIAGONALLY
*     UP THE TRIDIAGONAL
*
      IBUF1 = 0
      IBUF2 = 0
      IRCV1 = 0
      IRCV2 = 0
      DO 10 K = I, L + 1, -1
         IF( ( MODKM1.EQ.0 ) .AND. ( DOWN.EQ.II ) .AND.
     $       ( RIGHT.EQ.JJ ) ) THEN
*
*           WE MUST PACK H(K-1,K-1) AND SEND IT DIAGONAL DOWN
*
            IF( ( DOWN.NE.MYROW ) .OR. ( RIGHT.NE.MYCOL ) ) THEN
               CALL INFOG2L( K-1, K-1, DESCA, NPROW, NPCOL, MYROW,
     $                       MYCOL, IROW1, ICOL1, ISRC, JSRC )
               IBUF1 = IBUF1 + 1
               BUF( ISTR1+IBUF1 ) = A( ( ICOL1-1 )*LDA+IROW1 )
            END IF
         END IF
         IF( ( MODKM1.EQ.0 ) .AND. ( MYROW.EQ.II ) .AND.
     $       ( RIGHT.EQ.JJ ) ) THEN
*
*           WE MUST PACK H(K  ,K-1) AND SEND IT RIGHT
*
            IF( NPCOL.GT.1 ) THEN
               CALL INFOG2L( K, K-1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW1, ICOL1, ISRC, JSRC )
               IBUF2 = IBUF2 + 1
               BUF( ISTR2+IBUF2 ) = A( ( ICOL1-1 )*LDA+IROW1 )
            END IF
         END IF
*
*        ADD UP THE RECEIVES
*
         IF( ( MYROW.EQ.II ) .AND. ( MYCOL.EQ.JJ ) ) THEN
            IF( ( MODKM1.EQ.0 ) .AND. ( ( NPROW.GT.1 ) .OR. ( NPCOL.GT.
     $          1 ) ) ) THEN
*
*              WE MUST RECEIVE H(K-1,K-1) FROM DIAGONAL UP
*
               IRCV1 = IRCV1 + 1
            END IF
            IF( ( MODKM1.EQ.0 ) .AND. ( NPCOL.GT.1 ) ) THEN
*
*              WE MUST RECEIVE H(K  ,K-1) FROM LEFT
*
               IRCV2 = IRCV2 + 1
            END IF
         END IF
*
*        POSSIBLY CHANGE OWNERS (OCCURS ONLY WHEN MOD(K-1,HBL) = 0)
*
         IF( MODKM1.EQ.0 ) THEN
            II = II - 1
            JJ = JJ - 1
            IF( II.LT.0 )
     $         II = NPROW - 1
            IF( JJ.LT.0 )
     $         JJ = NPCOL - 1
         END IF
         MODKM1 = MODKM1 - 1
         IF( MODKM1.LT.0 )
     $      MODKM1 = HBL - 1
   10 CONTINUE
*
*     SEND DATA ON TO THE APPROPRIATE NODE IF THERE IS ANY DATA TO SEND
*
      IF( IBUF1.GT.0 ) THEN
         CALL CGESD2D( CONTXT, IBUF1, 1, BUF( ISTR1+1 ), IBUF1, DOWN,
     $                 RIGHT )
      END IF
      IF( IBUF2.GT.0 ) THEN
         CALL CGESD2D( CONTXT, IBUF2, 1, BUF( ISTR2+1 ), IBUF2, MYROW,
     $                 RIGHT )
      END IF
*
*     RECEIVE APPROPRIATE DATA IF THERE IS ANY
*
      IF( IRCV1.GT.0 ) THEN
         CALL CGERV2D( CONTXT, IRCV1, 1, BUF( ISTR1+1 ), IRCV1, UP,
     $                 LEFT )
      END IF
      IF( IRCV2.GT.0 ) THEN
         CALL CGERV2D( CONTXT, IRCV2, 1, BUF( ISTR2+1 ), IRCV2, MYROW,
     $                 LEFT )
      END IF
*
*     START MAIN LOOP
*
      IBUF1 = 0
      IBUF2 = 0
      CALL INFOG2L( I, I, DESCA, NPROW, NPCOL, MYROW, MYCOL, IROW1,
     $              ICOL1, II, JJ )
      MODKM1 = MOD( I-1+HBL, HBL )
*
*        LOOK FOR A SINGLE SMALL SUBDIAGONAL ELEMENT.
*
*        Start loop for subdiagonal search
*
      DO 40 K = I, L + 1, -1
         IF( ( MYROW.EQ.II ) .AND. ( MYCOL.EQ.JJ ) ) THEN
            IF( MODKM1.EQ.0 ) THEN
*
*                 Grab information from WORK array
*
               IF( NUM.GT.1 ) THEN
                  IBUF1 = IBUF1 + 1
                  H11 = BUF( ISTR1+IBUF1 )
               ELSE
                  H11 = A( ( ICOL1-2 )*LDA+IROW1-1 )
               END IF
               IF( NPCOL.GT.1 ) THEN
                  IBUF2 = IBUF2 + 1
                  H10 = BUF( ISTR2+IBUF2 )
               ELSE
                  H10 = A( ( ICOL1-2 )*LDA+IROW1 )
               END IF
            ELSE
*
*                 Information is local
*
               H11 = A( ( ICOL1-2 )*LDA+IROW1-1 )
               H10 = A( ( ICOL1-2 )*LDA+IROW1 )
            END IF
            H22 = A( ( ICOL1-1 )*LDA+IROW1 )
            TST1 = CABS1( H11 ) + CABS1( H22 )
            IF( TST1.EQ.ZERO ) THEN
*
*                 FIND SOME NORM OF THE LOCAL H(L:I,L:I)
*
               CALL INFOG1L( L, HBL, NPROW, MYROW, 0, IROW1, III )
               IROW2 = NUMROC( I, HBL, MYROW, 0, NPROW )
               CALL INFOG1L( L, HBL, NPCOL, MYCOL, 0, ICOL1, III )
               ICOL2 = NUMROC( I, HBL, MYCOL, 0, NPCOL )
               DO 30 III = IROW1, IROW2
                  DO 20 JJJ = ICOL1, ICOL2
                     TST1 = TST1 + CABS1( A( ( JJJ-1 )*LDA+III ) )
   20             CONTINUE
   30          CONTINUE
            END IF
            IF( CABS1( H10 ).LE.MAX( ULP*TST1, SMLNUM ) )
     $         GO TO 50
            IROW1 = IROW1 - 1
            ICOL1 = ICOL1 - 1
         END IF
         MODKM1 = MODKM1 - 1
         IF( MODKM1.LT.0 )
     $      MODKM1 = HBL - 1
         IF( ( MODKM1.EQ.HBL-1 ) .AND. ( K.GT.2 ) ) THEN
            II = MOD( II+NPROW-1, NPROW )
            JJ = MOD( JJ+NPCOL-1, NPCOL )
            CALL INFOG2L( K-1, K-1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IROW1, ICOL1, ITMP1, ITMP2 )
         END IF
   40 CONTINUE
   50 CONTINUE
      CALL IGAMX2D( CONTXT, 'ALL', ' ', 1, 1, K, 1, ITMP1, ITMP2, -1,
     $              -1, -1 )
      RETURN
*
*     End of PCLASMSUB
*
      END
