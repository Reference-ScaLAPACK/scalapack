      SUBROUTINE PCLACONSB( A, DESCA, I, L, M, H44, H33, H43H34, BUF,
     $                      LWORK )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     July 31, 2001
*
*     .. Scalar Arguments ..
      INTEGER            I, L, LWORK, M
      COMPLEX            H33, H43H34, H44
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX            A( * ), BUF( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLACONSB looks for two consecutive small subdiagonal elements by
*     seeing the effect of starting a double shift QR iteration
*     given by H44, H33, & H43H34 and see if this would make a
*     subdiagonal negligible.
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
*  A       (global input) COMPLEX array, dimension
*          (DESCA(LLD_),*)
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
*  M       (global output) INTEGER
*          On exit, this yields the starting location of the QR double
*          shift.  This will satisfy: L <= M  <= I-2.
*
*  H44
*  H33
*  H43H34  (global input) COMPLEX
*          These three values are for the double shift QR iteration.
*
*  BUF     (local output) COMPLEX array of size LWORK.
*
*  LWORK   (global input) INTEGER
*          On exit, LWORK is the size of the work buffer.
*          This must be at least 7*Ceil( Ceil( (I-L)/HBL ) /
*                                        LCM(NPROW,NPCOL) )
*          Here LCM is least common multiple, and NPROWxNPCOL is the
*          logical grid size.
*
*  Logic:
*  ======
*
*        Two consecutive small subdiagonal elements will stall
*        convergence of a double shift if their product is small
*        relatively even if each is not very small.  Thus it is
*        necessary to scan the "tridiagonal portion of the matrix."  In
*        the LAPACK algorithm ZLAHQR, a loop of M goes from I-2 down to
*        L and examines
*        H(m,m),H(m+1,m+1),H(m+1,m),H(m,m+1),H(m-1,m-1),H(m,m-1), and
*        H(m+2,m-1).  Since these elements may be on separate
*        processors, the first major loop (10) goes over the tridiagonal
*        and has each node store whatever values of the 7 it has that
*        the node owning H(m,m) does not.  This will occur on a border
*        and can happen in no more than 3 locations per block assuming
*        square blocks.  There are 5 buffers that each node stores these
*        values:  a buffer to send diagonally down and right, a buffer
*        to send up, a buffer to send left, a buffer to send diagonally
*        up and left and a buffer to send right.  Each of these buffers
*        is actually stored in one buffer BUF where BUF(ISTR1+1) starts
*        the first buffer, BUF(ISTR2+1) starts the second, etc..  After
*        the values are stored, if there are any values that a node
*        needs, they will be sent and received.  Then the next major
*        loop passes over the data and searches for two consecutive
*        small subdiagonals.
*
*  Notes:
*
*     This routine does a global maximum and must be called by all
*     processes.
*
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
      INTEGER            CONTXT, DOWN, HBL, IBUF1, IBUF2, IBUF3, IBUF4,
     $                   IBUF5, ICOL1, II, IRCV1, IRCV2, IRCV3, IRCV4,
     $                   IRCV5, IROW1, ISRC, ISTR1, ISTR2, ISTR3, ISTR4,
     $                   ISTR5, JJ, JSRC, LDA, LEFT, MODKM1, MYCOL,
     $                   MYROW, NPCOL, NPROW, NUM, RIGHT, UP
      REAL               S, TST1, ULP
      COMPLEX            CDUM, H00, H10, H11, H12, H21, H22, H33S, H44S,
     $                   V1, V2, V3
*     ..
*     .. External Functions ..
      INTEGER            ILCM
      REAL               PSLAMCH
      EXTERNAL           ILCM, PSLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGAMX2D, INFOG2L, PXERBLA,
     $                   CGERV2D, CGESD2D
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
      ULP = PSLAMCH( CONTXT, 'PRECISION' )
      CALL BLACS_GRIDINFO( CONTXT, NPROW, NPCOL, MYROW, MYCOL )
      LEFT = MOD( MYCOL+NPCOL-1, NPCOL )
      RIGHT = MOD( MYCOL+1, NPCOL )
      UP = MOD( MYROW+NPROW-1, NPROW )
      DOWN = MOD( MYROW+1, NPROW )
      NUM = NPROW*NPCOL
*
*     BUFFER1 starts at BUF(ISTR1+1) and will contain IBUF1 elements
*     BUFFER2 starts at BUF(ISTR2+1) and will contain IBUF2 elements
*     BUFFER3 starts at BUF(ISTR3+1) and will contain IBUF3 elements
*     BUFFER4 starts at BUF(ISTR4+1) and will contain IBUF4 elements
*     BUFFER5 starts at BUF(ISTR5+1) and will contain IBUF5 elements
*
      ISTR1 = 0
      ISTR2 = ( ( I-L-1 ) / HBL )
      IF( ISTR2*HBL.LT.( I-L-1 ) )
     $   ISTR2 = ISTR2 + 1
      II = ISTR2 / ILCM( NPROW, NPCOL )
      IF( II*ILCM( NPROW, NPCOL ).LT.ISTR2 ) THEN
         ISTR2 = II + 1
      ELSE
         ISTR2 = II
      END IF
      IF( LWORK.LT.7*ISTR2 ) THEN
         CALL PXERBLA( CONTXT, 'PCLACONSB', 10 )
         RETURN
      END IF
      ISTR3 = 3*ISTR2
      ISTR4 = ISTR3 + ISTR2
      ISTR5 = ISTR3 + ISTR3
      CALL INFOG2L( I-2, I-2, DESCA, NPROW, NPCOL, MYROW, MYCOL, IROW1,
     $              ICOL1, II, JJ )
      MODKM1 = MOD( I-3+HBL, HBL )
*
*     Copy our relevant pieces of triadiagonal that we owe into
*     5 buffers to send to whomever owns H(M,M) as M moves diagonally
*     up the tridiagonal
*
      IBUF1 = 0
      IBUF2 = 0
      IBUF3 = 0
      IBUF4 = 0
      IBUF5 = 0
      IRCV1 = 0
      IRCV2 = 0
      IRCV3 = 0
      IRCV4 = 0
      IRCV5 = 0
      DO 10 M = I - 2, L, -1
         IF( ( MODKM1.EQ.0 ) .AND. ( DOWN.EQ.II ) .AND.
     $       ( RIGHT.EQ.JJ ) .AND. ( M.GT.L ) ) THEN
*
*           We must pack H(M-1,M-1) and send it diagonal down
*
            IF( ( DOWN.NE.MYROW ) .OR. ( RIGHT.NE.MYCOL ) ) THEN
               CALL INFOG2L( M-1, M-1, DESCA, NPROW, NPCOL, MYROW,
     $                       MYCOL, IROW1, ICOL1, ISRC, JSRC )
               IBUF1 = IBUF1 + 1
               BUF( ISTR1+IBUF1 ) = A( ( ICOL1-1 )*LDA+IROW1 )
            END IF
         END IF
         IF( ( MODKM1.EQ.0 ) .AND. ( MYROW.EQ.II ) .AND.
     $       ( RIGHT.EQ.JJ ) .AND. ( M.GT.L ) ) THEN
*
*           We must pack H(M  ,M-1) and send it right
*
            IF( NPCOL.GT.1 ) THEN
               CALL INFOG2L( M, M-1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW1, ICOL1, ISRC, JSRC )
               IBUF5 = IBUF5 + 1
               BUF( ISTR5+IBUF5 ) = A( ( ICOL1-1 )*LDA+IROW1 )
            END IF
         END IF
         IF( ( MODKM1.EQ.HBL-1 ) .AND. ( UP.EQ.II ) .AND.
     $       ( MYCOL.EQ.JJ ) ) THEN
*
*           We must pack H(M+1,M) and send it up
*
            IF( NPROW.GT.1 ) THEN
               CALL INFOG2L( M+1, M, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW1, ICOL1, ISRC, JSRC )
               IBUF2 = IBUF2 + 1
               BUF( ISTR2+IBUF2 ) = A( ( ICOL1-1 )*LDA+IROW1 )
            END IF
         END IF
         IF( ( MODKM1.EQ.HBL-1 ) .AND. ( MYROW.EQ.II ) .AND.
     $       ( LEFT.EQ.JJ ) ) THEN
*
*           We must pack H(M  ,M+1) and send it left
*
            IF( NPCOL.GT.1 ) THEN
               CALL INFOG2L( M, M+1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW1, ICOL1, ISRC, JSRC )
               IBUF3 = IBUF3 + 1
               BUF( ISTR3+IBUF3 ) = A( ( ICOL1-1 )*LDA+IROW1 )
            END IF
         END IF
         IF( ( MODKM1.EQ.HBL-1 ) .AND. ( UP.EQ.II ) .AND.
     $       ( LEFT.EQ.JJ ) ) THEN
*
*           We must pack H(M+1,M+1) & H(M+2,M+1) and send it
*           diagonally up
*
            IF( ( UP.NE.MYROW ) .OR. ( LEFT.NE.MYCOL ) ) THEN
               CALL INFOG2L( M+1, M+1, DESCA, NPROW, NPCOL, MYROW,
     $                       MYCOL, IROW1, ICOL1, ISRC, JSRC )
               IBUF4 = IBUF4 + 2
               BUF( ISTR4+IBUF4-1 ) = A( ( ICOL1-1 )*LDA+IROW1 )
               BUF( ISTR4+IBUF4 ) = A( ( ICOL1-1 )*LDA+IROW1+1 )
            END IF
         END IF
         IF( ( MODKM1.EQ.HBL-2 ) .AND. ( UP.EQ.II ) .AND.
     $       ( MYCOL.EQ.JJ ) ) THEN
*
*           We must pack H(M+2,M+1) and send it up
*
            IF( NPROW.GT.1 ) THEN
               CALL INFOG2L( M+2, M+1, DESCA, NPROW, NPCOL, MYROW,
     $                       MYCOL, IROW1, ICOL1, ISRC, JSRC )
               IBUF2 = IBUF2 + 1
               BUF( ISTR2+IBUF2 ) = A( ( ICOL1-1 )*LDA+IROW1 )
            END IF
         END IF
*
*        Add up the receives
*
         IF( ( MYROW.EQ.II ) .AND. ( MYCOL.EQ.JJ ) ) THEN
            IF( ( MODKM1.EQ.0 ) .AND. ( M.GT.L ) .AND.
     $          ( ( NPROW.GT.1 ) .OR. ( NPCOL.GT.1 ) ) ) THEN
*
*              We must receive H(M-1,M-1) from diagonal up
*
               IRCV1 = IRCV1 + 1
            END IF
            IF( ( MODKM1.EQ.0 ) .AND. ( NPCOL.GT.1 ) .AND. ( M.GT.L ) )
     $           THEN
*
*              We must receive H(M  ,M-1) from left
*
               IRCV5 = IRCV5 + 1
            END IF
            IF( ( MODKM1.EQ.HBL-1 ) .AND. ( NPROW.GT.1 ) ) THEN
*
*              We must receive H(M+1,M  ) from down
*
               IRCV2 = IRCV2 + 1
            END IF
            IF( ( MODKM1.EQ.HBL-1 ) .AND. ( NPCOL.GT.1 ) ) THEN
*
*              We must receive H(M  ,M+1) from right
*
               IRCV3 = IRCV3 + 1
            END IF
            IF( ( MODKM1.EQ.HBL-1 ) .AND.
     $          ( ( NPROW.GT.1 ) .OR. ( NPCOL.GT.1 ) ) ) THEN
*
*              We must receive H(M+1:M+2,M+1) from diagonal down
*
               IRCV4 = IRCV4 + 2
            END IF
            IF( ( MODKM1.EQ.HBL-2 ) .AND. ( NPROW.GT.1 ) ) THEN
*
*              We must receive H(M+2,M+1) from down
*
               IRCV2 = IRCV2 + 1
            END IF
         END IF
*
*        Possibly change owners (occurs only when MOD(M-1,HBL) = 0)
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
*
*     Send data on to the appropriate node if there is any data to send
*
      IF( IBUF1.GT.0 ) THEN
         CALL CGESD2D( CONTXT, IBUF1, 1, BUF( ISTR1+1 ), IBUF1, DOWN,
     $                 RIGHT )
      END IF
      IF( IBUF2.GT.0 ) THEN
         CALL CGESD2D( CONTXT, IBUF2, 1, BUF( ISTR2+1 ), IBUF2, UP,
     $                 MYCOL )
      END IF
      IF( IBUF3.GT.0 ) THEN
         CALL CGESD2D( CONTXT, IBUF3, 1, BUF( ISTR3+1 ), IBUF3, MYROW,
     $                 LEFT )
      END IF
      IF( IBUF4.GT.0 ) THEN
         CALL CGESD2D( CONTXT, IBUF4, 1, BUF( ISTR4+1 ), IBUF4, UP,
     $                 LEFT )
      END IF
      IF( IBUF5.GT.0 ) THEN
         CALL CGESD2D( CONTXT, IBUF5, 1, BUF( ISTR5+1 ), IBUF5, MYROW,
     $                 RIGHT )
      END IF
*
*     Receive appropriate data if there is any
*
      IF( IRCV1.GT.0 ) THEN
         CALL CGERV2D( CONTXT, IRCV1, 1, BUF( ISTR1+1 ), IRCV1, UP,
     $                 LEFT )
      END IF
      IF( IRCV2.GT.0 ) THEN
         CALL CGERV2D( CONTXT, IRCV2, 1, BUF( ISTR2+1 ), IRCV2, DOWN,
     $                 MYCOL )
      END IF
      IF( IRCV3.GT.0 ) THEN
         CALL CGERV2D( CONTXT, IRCV3, 1, BUF( ISTR3+1 ), IRCV3, MYROW,
     $                 RIGHT )
      END IF
      IF( IRCV4.GT.0 ) THEN
         CALL CGERV2D( CONTXT, IRCV4, 1, BUF( ISTR4+1 ), IRCV4, DOWN,
     $                 RIGHT )
      END IF
      IF( IRCV5.GT.0 ) THEN
         CALL CGERV2D( CONTXT, IRCV5, 1, BUF( ISTR5+1 ), IRCV5, MYROW,
     $                 LEFT )
      END IF
*
*     Start main loop
*
      IBUF1 = 0
      IBUF2 = 0
      IBUF3 = 0
      IBUF4 = 0
      IBUF5 = 0
      CALL INFOG2L( I-2, I-2, DESCA, NPROW, NPCOL, MYROW, MYCOL, IROW1,
     $              ICOL1, II, JJ )
      MODKM1 = MOD( I-3+HBL, HBL )
      IF( ( MYROW.EQ.II ) .AND. ( MYCOL.EQ.JJ ) .AND.
     $    ( MODKM1.NE.HBL-1 ) ) THEN
         CALL INFOG2L( I-2, I-1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 IROW1, ICOL1, ISRC, JSRC )
      END IF
*
*     Look for two consecutive small subdiagonal elements.
*
      DO 20 M = I - 2, L, -1
*
*        Determine the effect of starting the double-shift QR
*        iteration at row M, and see if this would make H(M,M-1)
*        negligible.
*
         IF( ( MYROW.EQ.II ) .AND. ( MYCOL.EQ.JJ ) ) THEN
            IF( MODKM1.EQ.0 ) THEN
               H22 = A( ( ICOL1-1 )*LDA+IROW1+1 )
               H11 = A( ( ICOL1-2 )*LDA+IROW1 )
               V3 = A( ( ICOL1-1 )*LDA+IROW1+2 )
               H21 = A( ( ICOL1-2 )*LDA+IROW1+1 )
               H12 = A( ( ICOL1-1 )*LDA+IROW1 )
               IF( M.GT.L ) THEN
                  IF( NUM.GT.1 ) THEN
                     IBUF1 = IBUF1 + 1
                     H00 = BUF( ISTR1+IBUF1 )
                  ELSE
                     H00 = A( ( ICOL1-3 )*LDA+IROW1-1 )
                  END IF
                  IF( NPCOL.GT.1 ) THEN
                     IBUF5 = IBUF5 + 1
                     H10 = BUF( ISTR5+IBUF5 )
                  ELSE
                     H10 = A( ( ICOL1-3 )*LDA+IROW1 )
                  END IF
               END IF
            END IF
            IF( MODKM1.EQ.HBL-1 ) THEN
               CALL INFOG2L( M, M, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW1, ICOL1, ISRC, JSRC )
               H11 = A( ( ICOL1-1 )*LDA+IROW1 )
               IF( NUM.GT.1 ) THEN
                  IBUF4 = IBUF4 + 2
                  H22 = BUF( ISTR4+IBUF4-1 )
                  V3 = BUF( ISTR4+IBUF4 )
               ELSE
                  H22 = A( ICOL1*LDA+IROW1+1 )
                  V3 = A( ( ICOL1+1 )*LDA+IROW1+1 )
               END IF
               IF( NPROW.GT.1 ) THEN
                  IBUF2 = IBUF2 + 1
                  H21 = BUF( ISTR2+IBUF2 )
               ELSE
                  H21 = A( ( ICOL1-1 )*LDA+IROW1+1 )
               END IF
               IF( NPCOL.GT.1 ) THEN
                  IBUF3 = IBUF3 + 1
                  H12 = BUF( ISTR3+IBUF3 )
               ELSE
                  H12 = A( ICOL1*LDA+IROW1 )
               END IF
               IF( M.GT.L ) THEN
                  H00 = A( ( ICOL1-2 )*LDA+IROW1-1 )
                  H10 = A( ( ICOL1-2 )*LDA+IROW1 )
               END IF
*
*              Adjust ICOL1 for next iteration where MODKM1=HBL-2
*
               ICOL1 = ICOL1 + 1
            END IF
            IF( MODKM1.EQ.HBL-2 ) THEN
               H22 = A( ( ICOL1-1 )*LDA+IROW1+1 )
               H11 = A( ( ICOL1-2 )*LDA+IROW1 )
               IF( NPROW.GT.1 ) THEN
                  IBUF2 = IBUF2 + 1
                  V3 = BUF( ISTR2+IBUF2 )
               ELSE
                  V3 = A( ( ICOL1-1 )*LDA+IROW1+2 )
               END IF
               H21 = A( ( ICOL1-2 )*LDA+IROW1+1 )
               H12 = A( ( ICOL1-1 )*LDA+IROW1 )
               IF( M.GT.L ) THEN
                  H00 = A( ( ICOL1-3 )*LDA+IROW1-1 )
                  H10 = A( ( ICOL1-3 )*LDA+IROW1 )
               END IF
            END IF
            IF( ( MODKM1.LT.HBL-2 ) .AND. ( MODKM1.GT.0 ) ) THEN
               H22 = A( ( ICOL1-1 )*LDA+IROW1+1 )
               H11 = A( ( ICOL1-2 )*LDA+IROW1 )
               V3 = A( ( ICOL1-1 )*LDA+IROW1+2 )
               H21 = A( ( ICOL1-2 )*LDA+IROW1+1 )
               H12 = A( ( ICOL1-1 )*LDA+IROW1 )
               IF( M.GT.L ) THEN
                  H00 = A( ( ICOL1-3 )*LDA+IROW1-1 )
                  H10 = A( ( ICOL1-3 )*LDA+IROW1 )
               END IF
            END IF
            H44S = H44 - H11
            H33S = H33 - H11
            V1 = ( H33S*H44S-H43H34 ) / H21 + H12
            V2 = H22 - H11 - H33S - H44S
            S = CABS1( V1 ) + CABS1( V2 ) + CABS1( V3 )
            V1 = V1 / S
            V2 = V2 / S
            V3 = V3 / S
            IF( M.EQ.L )
     $         GO TO 30
            TST1 = CABS1( V1 )*( CABS1( H00 )+CABS1( H11 )+
     $             CABS1( H22 ) )
            IF( CABS1( H10 )*( CABS1( V2 )+CABS1( V3 ) ).LE.ULP*TST1 )
     $         GO TO 30
*
*           Slide indices diagonally up one for next iteration
*
            IROW1 = IROW1 - 1
            ICOL1 = ICOL1 - 1
         END IF
         IF( M.EQ.L ) THEN
*
*           Stop regardless of which node we are
*
            GO TO 30
         END IF
*
*        Possibly change owners if on border
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
   20 CONTINUE
   30 CONTINUE
*
      CALL IGAMX2D( CONTXT, 'ALL', ' ', 1, 1, M, 1, L, L, -1, -1, -1 )
*
      RETURN
*
*     End of PCLACONSB
*
      END
