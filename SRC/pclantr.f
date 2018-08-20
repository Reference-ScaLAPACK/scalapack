      REAL               FUNCTION PCLANTR( NORM, UPLO, DIAG, M, N, A,
     $                                     IA, JA, DESCA, WORK )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            IA, JA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      REAL               WORK( * )
      COMPLEX            A( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLANTR returns the value of the one norm, or the Frobenius norm,
*  or the infinity norm, or the element of largest absolute value of a
*  trapezoidal or triangular distributed matrix sub( A ) denoting
*  A(IA:IA+M-1, JA:JA+N-1).
*
*  PCLANTR returns the value
*
*     ( max(abs(A(i,j))),  NORM = 'M' or 'm' with ia <= i <= ia+m-1,
*     (                                      and  ja <= j <= ja+n-1,
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
*          Specifies the value to be returned in PCLANTR as described
*          above.
*
*  UPLO    (global input) CHARACTER
*          Specifies whether the matrix sub( A ) is upper or lower
*          trapezoidal.
*          = 'U':  Upper trapezoidal
*          = 'L':  Lower trapezoidal
*          Note that sub( A ) is triangular instead of trapezoidal
*          if M = N.
*
*  DIAG    (global input) CHARACTER
*          Specifies whether or not the distributed matrix sub( A ) has
*          unit diagonal.
*          = 'N':  Non-unit diagonal
*          = 'U':  Unit diagonal
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( A ). When M = 0, PCLANTR is
*          set to zero. M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( A ). When N = 0,
*          PCLANTR is set to zero. N >= 0.
*
*  A       (local input) COMPLEX pointer into the local memory
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
*  WORK    (local workspace) REAL array dimension (LWORK)
*          LWORK >=   0 if NORM = 'M' or 'm' (not referenced),
*                   Nq0 if NORM = '1', 'O' or 'o',
*                   Mp0 if NORM = 'I' or 'i',
*                     0 if NORM = 'F', 'f', 'E' or 'e' (not referenced),
*          where
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          Mp0 = NUMROC( M+IROFFA, MB_A, MYROW, IAROW, NPROW ),
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
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UDIAG
      INTEGER            IACOL, IAROW, ICTXT, II, IIA, ICOFF, IOFFA,
     $                   IROFF, J, JB, JJ, JJA, JN, KK, LDA, LL, MP,
     $                   MYCOL, MYROW, NP, NPCOL, NPROW, NQ
      REAL               SUM, VALUE
*     ..
*     .. Local Arrays ..
      REAL               SSQ( 2 ), COLSSQ( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CLASSQ, INFOG2L, PSTREECOMB,
     $                   SCOMBSSQ, SGEBR2D, SGEBS2D,
     $                   SGAMX2D, SGSUM2D
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, ISAMAX, NUMROC
      EXTERNAL           LSAME, ICEIL, ISAMAX, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      UDIAG = LSAME( DIAG, 'U' )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      IROFF = MOD( IA-1, DESCA( MB_ ) )
      ICOFF = MOD( JA-1, DESCA( NB_ ) )
      MP = NUMROC( M+IROFF, DESCA( MB_ ), MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFF, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IF( MYROW.EQ.IAROW )
     $   MP = MP - IROFF
      IF( MYCOL.EQ.IACOL )
     $   NQ = NQ - ICOFF
      LDA = DESCA( LLD_ )
      IOFFA = ( JJA - 1 ) * LDA
*
      IF( MIN( M, N ).EQ.0 ) THEN
*
         VALUE = ZERO
*
************************************************************************
* max norm
*
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         IF( UDIAG ) THEN
            VALUE = ONE
         ELSE
            VALUE = ZERO
         END IF
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Upper triangular matrix
*
            II = IIA
            JJ = JJA
            JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
            JB = JN-JA+1
*
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  IF( UDIAG ) THEN
                     DO 20 LL = JJ, JJ + JB -1
                        DO 10 KK = IIA, MIN(II+LL-JJ-1,IIA+MP-1)
                           VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
   10                   CONTINUE
                        IOFFA = IOFFA + LDA
   20                CONTINUE
                  ELSE
                     DO 40 LL = JJ, JJ + JB -1
                        DO 30 KK = IIA, MIN( II+LL-JJ, IIA+MP-1 )
                           VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
   30                   CONTINUE
                        IOFFA = IOFFA + LDA
   40                CONTINUE
                  END IF
               ELSE
                  DO 60 LL = JJ, JJ + JB -1
                     DO 50 KK = IIA, MIN( II-1, IIA+MP-1 )
                        VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
   50                CONTINUE
                     IOFFA = IOFFA + LDA
   60             CONTINUE
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 130 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     IF( UDIAG ) THEN
                        DO 80 LL = JJ, JJ + JB -1
                           DO 70 KK = IIA, MIN( II+LL-JJ-1, IIA+MP-1 )
                              VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
   70                      CONTINUE
                           IOFFA = IOFFA + LDA
   80                   CONTINUE
                     ELSE
                        DO 100 LL = JJ, JJ + JB -1
                           DO 90 KK = IIA, MIN( II+LL-JJ, IIA+MP-1 )
                              VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
   90                      CONTINUE
                           IOFFA = IOFFA + LDA
  100                   CONTINUE
                     END IF
                  ELSE
                     DO 120 LL = JJ, JJ + JB -1
                        DO 110 KK = IIA, MIN( II-1, IIA+MP-1 )
                           VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
  110                   CONTINUE
                        IOFFA = IOFFA + LDA
  120                CONTINUE
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  130       CONTINUE
*
         ELSE
*
*           Lower triangular matrix
*
            II = IIA
            JJ = JJA
            JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
            JB = JN-JA+1
*
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  IF( UDIAG ) THEN
                     DO 150 LL = JJ, JJ + JB -1
                        DO 140 KK = II+LL-JJ+1, IIA+MP-1
                           VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
  140                   CONTINUE
                        IOFFA = IOFFA + LDA
  150                CONTINUE
                  ELSE
                     DO 170 LL = JJ, JJ + JB -1
                        DO 160 KK = II+LL-JJ, IIA+MP-1
                           VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
  160                   CONTINUE
                        IOFFA = IOFFA + LDA
  170                CONTINUE
                  END IF
               ELSE
                  DO 190 LL = JJ, JJ + JB -1
                     DO 180 KK = II, IIA+MP-1
                        VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
  180                CONTINUE
                     IOFFA = IOFFA + LDA
  190             CONTINUE
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 260 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     IF( UDIAG ) THEN
                        DO 210 LL = JJ, JJ + JB -1
                           DO 200 KK = II+LL-JJ+1, IIA+MP-1
                              VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
  200                      CONTINUE
                           IOFFA = IOFFA + LDA
  210                   CONTINUE
                     ELSE
                        DO 230 LL = JJ, JJ + JB -1
                           DO 220 KK = II+LL-JJ, IIA+MP-1
                              VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
  220                      CONTINUE
                           IOFFA = IOFFA + LDA
  230                   CONTINUE
                     END IF
                  ELSE
                     DO 250 LL = JJ, JJ + JB -1
                        DO 240 KK = II, IIA+MP-1
                           VALUE = MAX( VALUE, ABS( A( IOFFA+KK ) ) )
  240                   CONTINUE
                        IOFFA = IOFFA + LDA
  250                CONTINUE
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  260       CONTINUE
*
         END IF
*
*        Gather the intermediate results to process (0,0).
*
         CALL SGAMX2D( ICTXT, 'All', ' ', 1, 1, VALUE, 1, KK, LL, -1,
     $                 0, 0 )
*
************************************************************************
* one norm
*
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' ) THEN
*
         VALUE = ZERO
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Upper triangular matrix
*
            II = IIA
            JJ = JJA
            JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
            JB = JN-JA+1
*
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  IF( UDIAG ) THEN
                     DO 280 LL = JJ, JJ + JB -1
                        SUM = ZERO
                        DO 270 KK = IIA, MIN( II+LL-JJ-1, IIA+MP-1 )
                           SUM = SUM + ABS( A( IOFFA+KK ) )
  270                   CONTINUE
*                       Unit diagonal entry
                        KK = II+LL-JJ
                        IF (KK <= IIA+MP-1) THEN
                           SUM = SUM + ONE
                        ENDIF
                        IOFFA = IOFFA + LDA
                        WORK( LL-JJA+1 ) = SUM
  280                CONTINUE
                  ELSE
                     DO 300 LL = JJ, JJ + JB -1
                        SUM = ZERO
                        DO 290 KK = IIA, MIN( II+LL-JJ, IIA+MP-1 )
                           SUM = SUM + ABS( A( IOFFA+KK ) )
  290                   CONTINUE
                        IOFFA = IOFFA + LDA
                        WORK( LL-JJA+1 ) = SUM
  300                CONTINUE
                  END IF
               ELSE
                  DO 320 LL = JJ, JJ + JB -1
                     SUM = ZERO
                     DO 310 KK = IIA, MIN( II-1, IIA+MP-1 )
                        SUM = SUM + ABS( A( IOFFA+KK ) )
  310                CONTINUE
                     IOFFA = IOFFA + LDA
                     WORK( LL-JJA+1 ) = SUM
  320             CONTINUE
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 390 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     IF( UDIAG ) THEN
                        DO 340 LL = JJ, JJ + JB -1
                           SUM = ZERO
                           DO 330 KK = IIA, MIN( II+LL-JJ-1, IIA+MP-1 )
                              SUM = SUM + ABS( A( IOFFA+KK ) )
  330                      CONTINUE
*                          Unit diagonal entry
                           KK = II+LL-JJ
                           IF (KK <= IIA+MP-1) THEN
                              SUM = SUM + ONE
                           ENDIF
                           IOFFA = IOFFA + LDA
                           WORK( LL-JJA+1 ) = SUM
  340                   CONTINUE
                     ELSE
                        DO 360 LL = JJ, JJ + JB -1
                           SUM = ZERO
                           DO 350 KK = IIA, MIN( II+LL-JJ, IIA+MP-1 )
                              SUM = SUM + ABS( A( IOFFA+KK ) )
  350                      CONTINUE
                           IOFFA = IOFFA + LDA
                           WORK( LL-JJA+1 ) = SUM
  360                   CONTINUE
                     END IF
                  ELSE
                     DO 380 LL = JJ, JJ + JB -1
                        SUM = ZERO
                        DO 370 KK = IIA, MIN( II-1, IIA+MP-1 )
                           SUM = SUM + ABS( A( IOFFA+KK ) )
  370                   CONTINUE
                        IOFFA = IOFFA + LDA
                        WORK( LL-JJA+1 ) = SUM
  380                CONTINUE
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  390       CONTINUE
*
         ELSE
*
*           Lower triangular matrix
*
            II = IIA
            JJ = JJA
            JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
            JB = JN-JA+1
*
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  IF( UDIAG ) THEN
                     DO 410 LL = JJ, JJ + JB -1
                        SUM = ONE
                        DO 400 KK = II+LL-JJ+1, IIA+MP-1
                           SUM = SUM + ABS( A( IOFFA+KK ) )
  400                   CONTINUE
                        IOFFA = IOFFA + LDA
                        WORK( LL-JJA+1 ) = SUM
  410                CONTINUE
                  ELSE
                     DO 430 LL = JJ, JJ + JB -1
                        SUM = ZERO
                        DO 420 KK = II+LL-JJ, IIA+MP-1
                           SUM = SUM + ABS( A( IOFFA+KK ) )
  420                   CONTINUE
                        IOFFA = IOFFA + LDA
                        WORK( LL-JJA+1 ) = SUM
  430                CONTINUE
                  END IF
               ELSE
                  DO 450 LL = JJ, JJ + JB -1
                     SUM = ZERO
                     DO 440 KK = II, IIA+MP-1
                        SUM = SUM + ABS( A( IOFFA+KK ) )
  440                CONTINUE
                     IOFFA = IOFFA + LDA
                     WORK( LL-JJA+1 ) = SUM
  450             CONTINUE
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 520 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     IF( UDIAG ) THEN
                        DO 470 LL = JJ, JJ + JB -1
                           SUM = ONE
                           DO 460 KK = II+LL-JJ+1, IIA+MP-1
                              SUM = SUM + ABS( A( IOFFA+KK ) )
  460                      CONTINUE
                           IOFFA = IOFFA + LDA
                           WORK( LL-JJA+1 ) = SUM
  470                   CONTINUE
                     ELSE
                        DO 490 LL = JJ, JJ + JB -1
                           SUM = ZERO
                           DO 480 KK = II+LL-JJ, IIA+MP-1
                              SUM = SUM + ABS( A( IOFFA+KK ) )
  480                      CONTINUE
                           IOFFA = IOFFA + LDA
                           WORK( LL-JJA+1 ) = SUM
  490                   CONTINUE
                     END IF
                  ELSE
                     DO 510 LL = JJ, JJ + JB -1
                        SUM = ZERO
                        DO 500 KK = II, IIA+MP-1
                           SUM = SUM + ABS( A( IOFFA+KK ) )
  500                   CONTINUE
                        IOFFA = IOFFA + LDA
                        WORK( LL-JJA+1 ) = SUM
  510                CONTINUE
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  520       CONTINUE
*
         END IF
*
*        Find sum of global matrix columns and store on row 0 of
*        process grid
*
         CALL SGSUM2D( ICTXT, 'Columnwise', ' ', 1, NQ, WORK, 1,
     $                 0, MYCOL )
*
*        Find maximum sum of columns for 1-norm
*
         IF( MYROW.EQ.0 ) THEN
            IF( NQ.GT.0 ) THEN
               VALUE = WORK( ISAMAX( NQ, WORK, 1 ) )
            ELSE
               VALUE = ZERO
            END IF
            CALL SGAMX2D( ICTXT, 'Rowwise', ' ', 1, 1, VALUE, 1, KK, LL,
     $                    -1, 0, 0 )
         END IF
*
************************************************************************
* infinity norm
*
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
         IF( LSAME( UPLO, 'U' ) ) THEN
               DO 540 KK = IIA, IIA+MP-1
                  WORK( KK ) = ZERO
  540          CONTINUE
         ELSE
               DO 570 KK = IIA, IIA+MP-1
                  WORK( KK ) = ZERO
  570          CONTINUE
         END IF
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Upper triangular matrix
*
            II = IIA
            JJ = JJA
            JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
            JB = JN-JA+1
*
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  IF( UDIAG ) THEN
                     DO 590 LL = JJ, JJ + JB -1
                        DO 580 KK = IIA, MIN( II+LL-JJ-1, IIA+MP-1 )
                           WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                        ABS( A( IOFFA+KK ) )
  580                   CONTINUE
*                       Unit diagonal entry
                        KK = II+LL-JJ
                        IF (KK <= IIA+MP-1) THEN
                           WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) + ONE
                        ENDIF
                        IOFFA = IOFFA + LDA
  590                CONTINUE
                  ELSE
                     DO 610 LL = JJ, JJ + JB -1
                        DO 600 KK = IIA, MIN( II+LL-JJ, IIA+MP-1 )
                           WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                        ABS( A( IOFFA+KK ) )
  600                   CONTINUE
                        IOFFA = IOFFA + LDA
  610                CONTINUE
                  END IF
               ELSE
                  DO 630 LL = JJ, JJ + JB -1
                     DO 620 KK = IIA, MIN( II-1, IIA+MP-1 )
                        WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                     ABS( A( IOFFA+KK ) )
  620                CONTINUE
                     IOFFA = IOFFA + LDA
  630             CONTINUE
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 700 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     IF( UDIAG ) THEN
                        DO 650 LL = JJ, JJ + JB -1
                           DO 640 KK = IIA, MIN( II+LL-JJ-1, IIA+MP-1 )
                              WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                           ABS( A( IOFFA+KK ) )
  640                      CONTINUE
*                          Unit diagonal entry
                           KK = II+LL-JJ
                           IF (KK <= IIA+MP-1) THEN
                              WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) + ONE
                           ENDIF
                           IOFFA = IOFFA + LDA
  650                   CONTINUE
                     ELSE
                        DO 670 LL = JJ, JJ + JB -1
                           DO 660 KK = IIA, MIN( II+LL-JJ, IIA+MP-1 )
                              WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                           ABS( A( IOFFA+KK ) )
  660                      CONTINUE
                           IOFFA = IOFFA + LDA
  670                   CONTINUE
                     END IF
                  ELSE
                     DO 690 LL = JJ, JJ + JB -1
                        DO 680 KK = IIA, MIN( II-1, IIA+MP-1 )
                           WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                        ABS( A( IOFFA+KK ) )
  680                   CONTINUE
                        IOFFA = IOFFA + LDA
  690                CONTINUE
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  700       CONTINUE
*
         ELSE
*
*           Lower triangular matrix
*
            II = IIA
            JJ = JJA
            JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
            JB = JN-JA+1
*
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  IF( UDIAG ) THEN
                     DO 720 LL = JJ, JJ + JB -1
*                       Unit diagonal entry
                        KK = II+LL-JJ
                        WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) + ONE
                        DO 710 KK = II+LL-JJ+1, IIA+MP-1
                           WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                        ABS( A( IOFFA+KK ) )
  710                   CONTINUE
                        IOFFA = IOFFA + LDA
  720                CONTINUE
                  ELSE
                     DO 740 LL = JJ, JJ + JB -1
                        DO 730 KK = II+LL-JJ, IIA+MP-1
                           WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                        ABS( A( IOFFA+KK ) )
  730                   CONTINUE
                        IOFFA = IOFFA + LDA
  740                CONTINUE
                  END IF
               ELSE
                  DO 760 LL = JJ, JJ + JB -1
                     DO 750 KK = II, IIA+MP-1
                        WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                     ABS( A( IOFFA+KK ) )
  750                CONTINUE
                     IOFFA = IOFFA + LDA
  760             CONTINUE
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 830 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     IF( UDIAG ) THEN
                        DO 780 LL = JJ, JJ + JB -1
*                          Unit diagonal entry
                           KK = II+LL-JJ
                           WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) + ONE
                           DO 770 KK = II+LL-JJ+1, IIA+MP-1
                              WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                           ABS( A( IOFFA+KK ) )
  770                      CONTINUE
                           IOFFA = IOFFA + LDA
  780                   CONTINUE
                     ELSE
                        DO 800 LL = JJ, JJ + JB -1
                           DO 790 KK = II+LL-JJ, IIA+MP-1
                              WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                           ABS( A( IOFFA+KK ) )
  790                      CONTINUE
                           IOFFA = IOFFA + LDA
  800                   CONTINUE
                     END IF
                  ELSE
                     DO 820 LL = JJ, JJ + JB -1
                        DO 810 KK = II, IIA+MP-1
                           WORK( KK-IIA+1 ) = WORK( KK-IIA+1 ) +
     $                                        ABS( A( IOFFA+KK ) )
  810                   CONTINUE
                        IOFFA = IOFFA + LDA
  820                CONTINUE
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  830       CONTINUE
*
         END IF
*
*        Find sum of global matrix rows and store on column 0 of
*        process grid
*
         CALL SGSUM2D( ICTXT, 'Rowwise', ' ', MP, 1, WORK, MAX( 1, MP ),
     $                 MYROW, 0 )
*
*        Find maximum sum of rows for Infinity-norm
*
         IF( MYCOL.EQ.0 ) THEN
            IF( MP.GT.0 ) THEN
               VALUE = WORK( ISAMAX( MP, WORK, 1 ) )
            ELSE
               VALUE = ZERO
            END IF
            CALL SGAMX2D( ICTXT, 'Columnwise', ' ', 1, 1, VALUE, 1, KK,
     $                    LL, -1, 0, 0 )
         END IF
*
************************************************************************
* Frobenius norm
* SSQ(1) is scale
* SSQ(2) is sum-of-squares
*
      ELSE IF( LSAME( NORM, 'F' ) .OR. LSAME( NORM, 'E' ) ) THEN
*
         IF( UDIAG ) THEN
            SSQ(1) = ONE
            SSQ(2) = REAL( MIN( M, N ) ) / REAL( NPROW*NPCOL )
         ELSE
            SSQ(1) = ZERO
            SSQ(2) = ONE
         END IF
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           ***********************
*           Upper triangular matrix
*
            II = IIA
            JJ = JJA
            JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
            JB = JN-JA+1
*
*           First block column of sub-matrix.
*
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
*                 This process has part of current block column,
*                 including diagonal block.
                  IF( UDIAG ) THEN
                     DO 840 LL = JJ, JJ + JB -1
                        COLSSQ(1) = ZERO
                        COLSSQ(2) = ONE
                        CALL CLASSQ( MIN( II+LL-JJ-1, IIA+MP-1 )-IIA+1,
     $                               A( IIA+IOFFA ), 1,
     $                               COLSSQ(1), COLSSQ(2) )
                        CALL SCOMBSSQ( SSQ, COLSSQ )
                        IOFFA = IOFFA + LDA
  840                CONTINUE
                  ELSE
                     DO 850 LL = JJ, JJ + JB -1
                        COLSSQ(1) = ZERO
                        COLSSQ(2) = ONE
                        CALL CLASSQ( MIN( II+LL-JJ, IIA+MP-1 )-IIA+1,
     $                               A( IIA+IOFFA ), 1,
     $                               COLSSQ(1), COLSSQ(2) )
                        CALL SCOMBSSQ( SSQ, COLSSQ )
                        IOFFA = IOFFA + LDA
  850                CONTINUE
                  END IF
               ELSE
*                 This rank has part of current block column,
*                 but not diagonal block.
*                 It seems this lassq will be length 0, since ii = iia.
                  DO 860 LL = JJ, JJ + JB -1
                     COLSSQ(1) = ZERO
                     COLSSQ(2) = ONE
                     CALL CLASSQ( MIN( II-1, IIA+MP-1 )-IIA+1,
     $                            A( IIA+IOFFA ), 1,
     $                            COLSSQ(1), COLSSQ(2) )
                     CALL SCOMBSSQ( SSQ, COLSSQ )
                     IOFFA = IOFFA + LDA
  860             CONTINUE
               END IF
               JJ = JJ + JB
            END IF
*
*           If this process has part of current block row, advance ii,
*           then advance iarow, iacol to next diagonal block.
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block columns
*
            DO 900 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     IF( UDIAG ) THEN
                        DO 870 LL = JJ, JJ + JB -1
                           COLSSQ(1) = ZERO
                           COLSSQ(2) = ONE
                           CALL CLASSQ( MIN(II+LL-JJ-1, IIA+MP-1)-IIA+1,
     $                                  A( IIA+IOFFA ), 1,
     $                                  COLSSQ(1), COLSSQ(2) )
                           CALL SCOMBSSQ( SSQ, COLSSQ )
                           IOFFA = IOFFA + LDA
  870                   CONTINUE
                     ELSE
                        DO 880 LL = JJ, JJ + JB -1
                           COLSSQ(1) = ZERO
                           COLSSQ(2) = ONE
                           CALL CLASSQ( MIN( II+LL-JJ, IIA+MP-1 )-IIA+1,
     $                                  A( IIA+IOFFA ), 1,
     $                                  COLSSQ(1), COLSSQ(2) )
                           CALL SCOMBSSQ( SSQ, COLSSQ )
                           IOFFA = IOFFA + LDA
  880                   CONTINUE
                     END IF
                  ELSE
                     DO 890 LL = JJ, JJ + JB -1
                        COLSSQ(1) = ZERO
                        COLSSQ(2) = ONE
                        CALL CLASSQ( MIN( II-1, IIA+MP-1 )-IIA+1,
     $                               A( IIA+IOFFA ), 1,
     $                               COLSSQ(1), COLSSQ(2) )
                        CALL SCOMBSSQ( SSQ, COLSSQ )
                        IOFFA = IOFFA + LDA
  890                CONTINUE
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  900       CONTINUE
*
         ELSE
*
*           ***********************
*           Lower triangular matrix
*
            II = IIA
            JJ = JJA
            JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
            JB = JN-JA+1
*
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  IF( UDIAG ) THEN
                     DO 910 LL = JJ, JJ + JB -1
                        COLSSQ(1) = ZERO
                        COLSSQ(2) = ONE
                        CALL CLASSQ( IIA+MP-(II+LL-JJ+1),
     $                               A( II+LL-JJ+1+IOFFA ), 1,
     $                               COLSSQ(1), COLSSQ(2) )
                        CALL SCOMBSSQ( SSQ, COLSSQ )
                        IOFFA = IOFFA + LDA
  910                CONTINUE
                  ELSE
                     DO 920 LL = JJ, JJ + JB -1
                        COLSSQ(1) = ZERO
                        COLSSQ(2) = ONE
                        CALL CLASSQ( IIA+MP-(II+LL-JJ),
     $                               A( II+LL-JJ+IOFFA ), 1,
     $                               COLSSQ(1), COLSSQ(2) )
                        CALL SCOMBSSQ( SSQ, COLSSQ )
                        IOFFA = IOFFA + LDA
  920                CONTINUE
                  END IF
               ELSE
                  DO 930 LL = JJ, JJ + JB -1
                     COLSSQ(1) = ZERO
                     COLSSQ(2) = ONE
                     CALL CLASSQ( IIA+MP-II, A( II+IOFFA ), 1,
     $                            COLSSQ(1), COLSSQ(2) )
                     CALL SCOMBSSQ( SSQ, COLSSQ )
                     IOFFA = IOFFA + LDA
  930             CONTINUE
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 970 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     IF( UDIAG ) THEN
                        DO 940 LL = JJ, JJ + JB -1
                           COLSSQ(1) = ZERO
                           COLSSQ(2) = ONE
                           CALL CLASSQ( IIA+MP-(II+LL-JJ+1),
     $                                  A( II+LL-JJ+1+IOFFA ), 1,
     $                                  COLSSQ(1), COLSSQ(2) )
                           CALL SCOMBSSQ( SSQ, COLSSQ )
                           IOFFA = IOFFA + LDA
  940                   CONTINUE
                     ELSE
                        DO 950 LL = JJ, JJ + JB -1
                           COLSSQ(1) = ZERO
                           COLSSQ(2) = ONE
                           CALL CLASSQ( IIA+MP-(II+LL-JJ),
     $                                  A( II+LL-JJ+IOFFA ), 1,
     $                                  COLSSQ(1), COLSSQ(2) )
                           CALL SCOMBSSQ( SSQ, COLSSQ )
                           IOFFA = IOFFA + LDA
  950                   CONTINUE
                     END IF
                  ELSE
                     DO 960 LL = JJ, JJ + JB -1
                        COLSSQ(1) = ZERO
                        COLSSQ(2) = ONE
                        CALL CLASSQ( IIA+MP-II, A( II+IOFFA ), 1,
     $                               COLSSQ(1), COLSSQ(2) )
                        CALL SCOMBSSQ( SSQ, COLSSQ )
                        IOFFA = IOFFA + LDA
  960                CONTINUE
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  970       CONTINUE
*
         END IF
*
*        ***********************
*        Perform the global scaled sum
*
         CALL PSTREECOMB( ICTXT, 'All', 2, SSQ, 0, 0, SCOMBSSQ )
         VALUE = SSQ( 1 ) * SQRT( SSQ( 2 ) )
*
      END IF
*
*     Broadcast the result to every process in the grid.
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         CALL SGEBS2D( ICTXT, 'All', ' ', 1, 1, VALUE, 1 )
      ELSE
         CALL SGEBR2D( ICTXT, 'All', ' ', 1, 1, VALUE, 1, 0, 0 )
      END IF
*
      PCLANTR = VALUE
*
      RETURN
*
*     End of PCLANTR
*
      END
