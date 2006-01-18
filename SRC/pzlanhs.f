      DOUBLE PRECISION   FUNCTION PZLANHS( NORM, N, A, IA, JA, DESCA,
     $                                     WORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            IA, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         A( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLANHS returns the value of the one norm, or the Frobenius norm,
*  or the infinity norm, or the element of largest absolute value of a
*  Hessenberg distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1).
*
*  PZLANHS returns the value
*
*     ( max(abs(A(i,j))),  NORM = 'M' or 'm' with IA <= i <= IA+N-1,
*     (                                      and  JA <= j <= JA+N-1,
*     (
*     ( norm1( sub( A ) ), NORM = '1', 'O' or 'o'
*     (
*     ( normI( sub( A ) ), NORM = 'I' or 'i'
*     (
*     ( normF( sub( A ) ), NORM = 'F', 'f', 'E' or 'e'
*
*  where norm1 denotes the  one norm of a matrix (maximum column sum),
*  normI denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
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
*  NORM    (global input) CHARACTER
*          Specifies the value to be returned in PZLANHS as described
*          above.
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on i.e the
*          number of rows and columns of the distributed submatrix
*          sub( A ). When N = 0, PZLANHS is set to zero. N >= 0.
*
*  A       (local input) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_A, LOCc(JA+N-1) ) containing
*          the local pieces of sub( A ).
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  WORK    (local workspace) DOUBLE PRECISION array dimension (LWORK)
*          LWORK >=   0 if NORM = 'M' or 'm' (not referenced),
*                   Nq0 if NORM = '1', 'O' or 'o',
*                   Mp0 if NORM = 'I' or 'i',
*                     0 if NORM = 'F', 'f', 'E' or 'e' (not referenced),
*          where
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          Np0 = NUMROC( N+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*          Nq0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),
*
*          INDXG2P and NUMROC are ScaLAPACK tool functions; MYROW,
*          MYCOL, NPROW and NPCOL can be determined by calling the
*          subroutine BLACS_GRIDINFO.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IACOL, IAROW, ICTXT, II, IIA, ICOFF, INXTROW,
     $                   IOFFA, IROFF, J, JB, JJ, JJA, JN, KK, LDA, LL,
     $                   MYCOL, MYROW, NP, NPCOL, NPROW, NQ
      DOUBLE PRECISION   SCALE, SUM, VALUE
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DCOMBSSQ, DGEBR2D,
     $                   DGEBS2D, DGAMX2D, DGSUM2D,
     $                   INFOG2L, PDTREECOMB, ZLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, IDAMAX, NUMROC
      EXTERNAL           LSAME, ICEIL, IDAMAX, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD, SQRT
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters.
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      IROFF = MOD( IA-1, DESCA( MB_ ) )
      ICOFF = MOD( JA-1, DESCA( NB_ ) )
      NP = NUMROC( N+IROFF, DESCA( MB_ ), MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFF, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IF( MYROW.EQ.IAROW )
     $   NP = NP - IROFF
      IF( MYCOL.EQ.IACOL )
     $   NQ = NQ - ICOFF
      LDA = DESCA( LLD_ )
      IOFFA = ( JJA - 1 ) * LDA
*
      IF( N.EQ.0 ) THEN
*
         VALUE = ZERO
*
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
         VALUE = ZERO
*
*        Find max(abs(A(i,j))).
*
         II = IIA
         JJ = JJA
         JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
         JB = JN-JA+1
*
*        Only one process row
*
         IF( NPROW.EQ.1 ) THEN
*
*           Handle first block of columns separately
*
            IF( MYCOL.EQ.IACOL ) THEN
               DO 20 LL = JJ, JJ+JB-1
                  DO 10 KK = IIA, MIN( II+LL-JJ+1, IIA+NP-1 )
                     VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
   10             CONTINUE
                  IOFFA = IOFFA + LDA
   20          CONTINUE
               JJ = JJ + JB
            END IF
*
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 50 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  DO 40 LL = JJ, JJ+JB-1
                     DO 30 KK = IIA, MIN( II+LL-JJ+1, IIA+NP-1 )
                        VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
   30                CONTINUE
                     IOFFA = IOFFA + LDA
   40             CONTINUE
                  JJ = JJ + JB
               END IF
*
               II = II + JB
               IACOL = MOD( IACOL+1, NPCOL )
*
   50       CONTINUE
*
         ELSE
*
*           Handle first block of columns separately
*
            INXTROW = MOD( IAROW+1, NPROW )
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  DO 70 LL = JJ, JJ + JB -1
                     DO 60 KK = IIA, MIN( II+LL-JJ+1, IIA+NP-1 )
                        VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
   60                CONTINUE
                     IOFFA = IOFFA + LDA
   70             CONTINUE
               ELSE
                  DO 90 LL = JJ, JJ+JB-1
                     DO 80 KK = IIA, MIN( II-1, IIA+NP-1 )
                        VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
   80                CONTINUE
                     IOFFA = IOFFA + LDA
   90             CONTINUE
                  IF( MYROW.EQ.INXTROW .AND. II.LE.IIA+NP-1 )
     $               VALUE = MAX( VALUE, ABS( A( II+(JJ+JB-2)*LDA ) ) )
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = INXTROW
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 140 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     DO 110 LL = JJ, JJ + JB -1
                        DO 100 KK = IIA, MIN( II+LL-JJ+1, IIA+NP-1 )
                           VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
  100                   CONTINUE
                        IOFFA = IOFFA + LDA
  110                CONTINUE
                  ELSE
                     DO 130 LL = JJ, JJ + JB -1
                        DO 120 KK = IIA, MIN( II-1, IIA+NP-1 )
                           VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
  120                   CONTINUE
                        IOFFA = IOFFA + LDA
  130                CONTINUE
                     IF( MYROW.EQ.INXTROW .AND. II.LE.IIA+NP-1 )
     $                  VALUE = MAX( VALUE,
     $                               ABS( A( II+(JJ+JB-2)*LDA ) ) )
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = INXTROW
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  140       CONTINUE
*
         END IF
*
*        Gather the intermediate results to process (0,0).
*
         CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, VALUE, 1, KK, LL, -1,
     $                 0, 0 )
*
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' ) THEN
*
         VALUE = ZERO
         II = IIA
         JJ = JJA
         JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
         JB = JN-JA+1
*
*        Only one process row
*
         IF( NPROW.EQ.1 ) THEN
*
*           Handle first block of columns separately
*
            IF( MYCOL.EQ.IACOL ) THEN
               DO 160 LL = JJ, JJ+JB-1
                  SUM = ZERO
                  DO 150 KK = IIA, MIN( II+LL-JJ+1, IIA+NP-1 )
                     SUM = SUM + ABS( A( IOFFA+KK ) )
  150             CONTINUE
                  IOFFA = IOFFA + LDA
                  WORK( LL-JJA+1 ) = SUM
  160          CONTINUE
               JJ = JJ + JB
            END IF
*
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 190 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  DO 180 LL = JJ, JJ+JB-1
                     SUM = ZERO
                     DO 170 KK = IIA, MIN( II+LL-JJ+1, IIA+NP-1 )
                        SUM = SUM + ABS( A( IOFFA+KK ) )
  170                CONTINUE
                     IOFFA = IOFFA + LDA
                     WORK( LL-JJA+1 ) = SUM
  180             CONTINUE
                  JJ = JJ + JB
               END IF
*
               II = II + JB
               IACOL = MOD( IACOL+1, NPCOL )
*
  190       CONTINUE
*
         ELSE
*
*           Handle first block of columns separately
*
            INXTROW = MOD( IAROW+1, NPROW )
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  DO 210 LL = JJ, JJ + JB -1
                     SUM = ZERO
                     DO 200 KK = IIA, MIN( II+LL-JJ+1, IIA+NP-1 )
                        SUM = SUM + ABS( A( IOFFA+KK ) )
  200                CONTINUE
                     IOFFA = IOFFA + LDA
                     WORK( LL-JJA+1 ) = SUM
  210             CONTINUE
               ELSE
                  DO 230 LL = JJ, JJ + JB -1
                     SUM = ZERO
                     DO 220 KK = IIA, MIN( II-1, IIA+NP-1 )
                        SUM = SUM + ABS( A( IOFFA+KK ) )
  220                CONTINUE
                     IOFFA = IOFFA + LDA
                     WORK( LL-JJA+1 ) = SUM
  230             CONTINUE
                  IF( MYROW.EQ.INXTROW .AND. II.LE.IIA+NP-1 )
     $               WORK( JJ+JB-JJA ) = WORK( JJ+JB-JJA ) +
     $                                   ABS( A( II+(JJ+JB-2)*LDA ) )
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = INXTROW
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 280 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     DO 250 LL = JJ, JJ + JB -1
                        SUM = ZERO
                        DO 240 KK = IIA, MIN( II+LL-JJ+1, IIA+NP-1 )
                           SUM = SUM + ABS( A( IOFFA+KK ) )
  240                   CONTINUE
                        IOFFA = IOFFA + LDA
                        WORK( LL-JJA+1 ) = SUM
  250                CONTINUE
                  ELSE
                     DO 270 LL = JJ, JJ + JB -1
                        SUM = ZERO
                        DO 260 KK = IIA, MIN( II-1, IIA+NP-1 )
                           SUM = SUM + ABS( A( IOFFA+KK ) )
  260                   CONTINUE
                        IOFFA = IOFFA + LDA
                        WORK( LL-JJA+1 ) = SUM
  270                CONTINUE
                     IF( MYROW.EQ.INXTROW .AND. II.LE.IIA+NP-1 )
     $                  WORK( JJ+JB-JJA ) = WORK( JJ+JB-JJA ) +
     $                                      ABS( A( II+(JJ+JB-2)*LDA ) )
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = INXTROW
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  280       CONTINUE
*
         END IF
*
*        Find sum of global matrix columns and store on row 0 of
*        process grid
*
         CALL DGSUM2D( ICTXT, 'Columnwise', ' ', 1, NQ, WORK, 1,
     $                 0, MYCOL )
*
*        Find maximum sum of columns for 1-norm
*
         IF( MYROW.EQ.0 ) THEN
            IF( NQ.GT.0 ) THEN
               VALUE = WORK( IDAMAX( NQ, WORK, 1 ) )
            ELSE
               VALUE = ZERO
            END IF
            CALL DGAMX2D( ICTXT, 'Rowwise', ' ', 1, 1, VALUE, 1, KK, LL,
     $                    -1, 0, 0 )
         END IF
*
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
         DO 290 KK = IIA, IIA+NP-1
            WORK( KK ) = ZERO
  290    CONTINUE
*
         II = IIA
         JJ = JJA
         JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
         JB = JN-JA+1
*
*        Only one process row
*
         IF( NPROW.EQ.1 ) THEN
*
*           Handle first block of columns separately
*
            IF( MYCOL.EQ.IACOL ) THEN
               DO 310 LL = JJ, JJ+JB-1
                  DO 300 KK = IIA, MIN( II+LL-JJ+1, IIA+NP-1 )
                     WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                  ABS( A( IOFFA+KK ) )
  300             CONTINUE
                  IOFFA = IOFFA + LDA
  310          CONTINUE
               JJ = JJ + JB
            END IF
*
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 340 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  DO 330 LL = JJ, JJ+JB-1
                     DO 320 KK = IIA, MIN( II+LL-JJ+1, IIA+NP-1 )
                        WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                     ABS( A( IOFFA+KK ) )
  320                CONTINUE
                     IOFFA = IOFFA + LDA
  330             CONTINUE
                  JJ = JJ + JB
               END IF
*
               II = II + JB
               IACOL = MOD( IACOL+1, NPCOL )
*
  340       CONTINUE
*
         ELSE
*
*           Handle first block of columns separately
*
            INXTROW = MOD( IAROW+1, NPROW )
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  DO 360 LL = JJ, JJ + JB -1
                     DO 350 KK = IIA, MIN( II+LL-JJ+1, IIA+NP-1 )
                        WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                     ABS( A( IOFFA+KK ) )
  350                CONTINUE
                     IOFFA = IOFFA + LDA
  360             CONTINUE
               ELSE
                  DO 380 LL = JJ, JJ + JB -1
                     DO 370 KK = IIA, MIN( II-1, IIA+NP-1 )
                        WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                     ABS( A( IOFFA+KK ) )
  370                CONTINUE
                     IOFFA = IOFFA + LDA
  380             CONTINUE
                  IF( MYROW.EQ.INXTROW .AND. II.LE.IIA+NP-1 )
     $               WORK( II-IIA+1 ) = WORK( II-IIA+1 ) +
     $                                  ABS( A( II+(JJ+JB-2)*LDA ) )
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = INXTROW
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 430 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     DO 400 LL = JJ, JJ + JB -1
                        DO 390 KK = IIA, MIN( II+LL-JJ+1, IIA+NP-1 )
                           WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                        ABS( A( IOFFA+KK ) )
  390                   CONTINUE
                        IOFFA = IOFFA + LDA
  400                CONTINUE
                  ELSE
                     DO 420 LL = JJ, JJ + JB -1
                        DO 410 KK = IIA, MIN( II-1, IIA+NP-1 )
                           WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                        ABS(A(IOFFA+KK))
  410                   CONTINUE
                        IOFFA = IOFFA + LDA
  420                CONTINUE
                     IF( MYROW.EQ.INXTROW .AND. II.LE.IIA+NP-1 )
     $                  WORK( II-IIA+1 ) = WORK( II-IIA+1 ) +
     $                                     ABS( A( II+(JJ+JB-2)*LDA ) )
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = INXTROW
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  430       CONTINUE
*
         END IF
*
*        Find sum of global matrix rows and store on column 0 of
*        process grid
*
         CALL DGSUM2D( ICTXT, 'Rowwise', ' ', NP, 1, WORK, MAX( 1, NP ),
     $                 MYROW, 0 )
*
*        Find maximum sum of rows for Infinity-norm
*
         IF( MYCOL.EQ.0 ) THEN
            IF( NP.GT.0 ) THEN
               VALUE = WORK( IDAMAX( NP, WORK, 1 ) )
            ELSE
               VALUE = ZERO
            END IF
            CALL DGAMX2D( ICTXT, 'Columnwise', ' ', 1, 1, VALUE, 1, KK,
     $                    LL, -1, 0, 0 )
         END IF
*
      ELSE IF( LSAME( NORM, 'F' ) .OR. LSAME( NORM, 'E' ) ) THEN
*
         SCALE = ZERO
         SUM = ONE
         II = IIA
         JJ = JJA
         JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
         JB = JN-JA+1
*
*        Only one process row
*
         IF( NPROW.EQ.1 ) THEN
*
*           Handle first block of columns separately
*
            IF( MYCOL.EQ.IACOL ) THEN
               DO 440 LL = JJ, JJ+JB-1
                  CALL ZLASSQ( MIN( II+LL-JJ+1, IIA+NP-1 )-IIA+1,
     $                         A( IIA+IOFFA ), 1, SCALE, SUM )
                  IOFFA = IOFFA + LDA
  440          CONTINUE
               JJ = JJ + JB
            END IF
*
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 460 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  DO 450 LL = JJ, JJ+JB-1
                     CALL ZLASSQ( MIN( II+LL-JJ+1, IIA+NP-1 )-IIA+1,
     $                            A( IIA+IOFFA ), 1, SCALE, SUM )
                     IOFFA = IOFFA + LDA
  450             CONTINUE
                  JJ = JJ + JB
               END IF
*
               II = II + JB
               IACOL = MOD( IACOL+1, NPCOL )
*
  460       CONTINUE
*
         ELSE
*
*           Handle first block of columns separately
*
            INXTROW = MOD( IAROW+1, NPROW )
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  DO 470 LL = JJ, JJ + JB -1
                     CALL ZLASSQ( MIN( II+LL-JJ+1, IIA+NP-1 )-IIA+1,
     $                            A( IIA+IOFFA ), 1, SCALE, SUM )
                     IOFFA = IOFFA + LDA
  470             CONTINUE
               ELSE
                  DO 480 LL = JJ, JJ + JB -1
                     CALL ZLASSQ( MIN( II-1, IIA+NP-1 )-IIA+1,
     $                            A( IIA+IOFFA ), 1, SCALE, SUM )
                     IOFFA = IOFFA + LDA
  480             CONTINUE
                  IF( MYROW.EQ.INXTROW .AND. II.LE.IIA+NP-1 )
     $               CALL ZLASSQ( 1, A( II+(JJ+JB-2)*LDA ), 1,
     $                            SCALE, SUM )
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = INXTROW
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 510 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     DO 490 LL = JJ, JJ + JB -1
                        CALL ZLASSQ( MIN( II+LL-JJ+1, IIA+NP-1 )-IIA+1,
     $                               A( IIA+IOFFA ), 1, SCALE, SUM )
                        IOFFA = IOFFA + LDA
  490                CONTINUE
                  ELSE
                     DO 500 LL = JJ, JJ + JB -1
                        CALL ZLASSQ( MIN( II-1, IIA+NP-1 )-IIA+1,
     $                               A( IIA+IOFFA ), 1, SCALE, SUM )
                        IOFFA = IOFFA + LDA
  500                CONTINUE
                     IF( MYROW.EQ.INXTROW .AND. II.LE.IIA+NP-1 )
     $                  CALL ZLASSQ( 1, A( II+(JJ+JB-2)*LDA ), 1,
     $                               SCALE, SUM )
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = INXTROW
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  510       CONTINUE
*
         END IF
*
*        Perform the global scaled sum
*
         RWORK( 1 ) = SCALE
         RWORK( 2 ) = SUM
         CALL PDTREECOMB( ICTXT, 'All', 2, RWORK, 0, 0, DCOMBSSQ )
         VALUE = RWORK( 1 ) * SQRT( RWORK( 2 ) )
*
      END IF
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         CALL DGEBS2D( ICTXT, 'All', ' ', 1, 1, VALUE, 1 )
      ELSE
         CALL DGEBR2D( ICTXT, 'All', ' ', 1, 1, VALUE, 1, 0, 0 )
      END IF
*
      PZLANHS = VALUE
*
      RETURN
*
*     End of PZLANHS
*
      END
