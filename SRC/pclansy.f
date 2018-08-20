      REAL             FUNCTION PCLANSY( NORM, UPLO, N, A, IA, JA,
     $                                   DESCA, WORK )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            IA, JA, N
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
*  PCLANSY  returns the value of the one norm, or the Frobenius norm,
*  or the infinity norm, or the element of largest absolute value of a
*  real symmetric distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1).
*
*  PCLANSY returns the value
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
*  where norm1  denotes the  one norm of a matrix (maximum column sum),
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
*          Specifies the value to be returned in PCLANSY as described
*          above.
*
*  UPLO    (global input) CHARACTER
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix sub( A ) is to be referenced.
*          = 'U':  Upper triangular part of sub( A ) is referenced,
*          = 'L':  Lower triangular part of sub( A ) is referenced.
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on i.e the
*          number of rows and columns of the distributed submatrix
*          sub( A ). When N = 0, PCLANSY is set to zero. N >= 0.
*
*  A       (local input) COMPLEX pointer into the local memory
*          to an array of dimension (LLD_A, LOCc(JA+N-1)) containing the
*          local pieces of the symmetric distributed matrix sub( A ).
*          If UPLO = 'U', the leading N-by-N upper triangular part of
*          sub( A ) contains the upper triangular matrix which norm is
*          to be computed, and the strictly lower triangular part of
*          this matrix is not referenced.  If UPLO = 'L', the leading
*          N-by-N lower triangular part of sub( A ) contains the lower
*          triangular matrix which norm is to be computed, and the
*          strictly upper triangular part of sub( A ) is not referenced.
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
*          LWORK >= 0 if NORM = 'M' or 'm' (not referenced),
*                   2*Nq0+Np0+LDW if NORM = '1', 'O', 'o', 'I' or 'i',
*                     where LDW is given by:
*                     IF( NPROW.NE.NPCOL ) THEN
*                        LDW = MB_A*CEIL(CEIL(Np0/MB_A)/(LCM/NPROW))
*                     ELSE
*                        LDW = 0
*                     END IF
*                   0 if NORM = 'F', 'f', 'E' or 'e' (not referenced),
*
*          where LCM is the least common multiple of NPROW and NPCOL
*          LCM = ILCM( NPROW, NPCOL ) and CEIL denotes the ceiling
*          operation (ICEIL).
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          Np0 = NUMROC( N+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*          Nq0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),
*
*          ICEIL, ILCM, INDXG2P and NUMROC are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
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
      INTEGER            I, IAROW, IACOL, IB, ICOFF, ICTXT, ICURCOL,
     $                   ICURROW, II, IIA, IN, IROFF, ICSR, ICSR0,
     $                   IOFFA, IRSC, IRSC0, IRSR, IRSR0, JJ, JJA, K,
     $                   LDA, LL, MYCOL, MYROW, NP, NPCOL, NPROW, NQ
      REAL               SUM, VALUE
*     ..
*     .. Local Arrays ..
      REAL               SSQ( 2 ), COLSSQ( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CLASSQ, PSCOL2ROW,
     $                   PSTREECOMB, SAXPY, SCOMBSSQ,
     $                   SGAMX2D, SGSUM2D, SGEBR2D, SGEBS2D
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, ISAMAX, NUMROC
      EXTERNAL           ICEIL, ISAMAX, LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD, SQRT
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters and local indexes.
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $              IIA, JJA, IAROW, IACOL )
*
      IROFF = MOD( IA-1, DESCA( MB_ ) )
      ICOFF = MOD( JA-1, DESCA( NB_ ) )
      NP = NUMROC( N+IROFF, DESCA( MB_ ), MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFF, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      ICSR = 1
      IRSR = ICSR + NQ
      IRSC = IRSR + NQ
      IF( MYROW.EQ.IAROW ) THEN
         IRSC0 = IRSC + IROFF
         NP = NP - IROFF
      ELSE
         IRSC0 = IRSC
      END IF
      IF( MYCOL.EQ.IACOL ) THEN
         ICSR0 = ICSR + ICOFF
         IRSR0 = IRSR + ICOFF
         NQ = NQ - ICOFF
      ELSE
         ICSR0 = ICSR
         IRSR0 = IRSR
      END IF
      IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+N-1 )
      LDA = DESCA( LLD_ )
*
*     If the matrix is symmetric, we address only a triangular portion
*     of the matrix.  A sum of row (column) i of the complete matrix
*     can be obtained by adding along row i and column i of the the
*     triangular matrix, stopping/starting at the diagonal, which is
*     the point of reflection.  The pictures below demonstrate this.
*     In the following code, the row sums created by --- rows below are
*     refered to as ROWSUMS, and the column sums shown by | are refered
*     to as COLSUMS. Infinity-norm = 1-norm = ROWSUMS+COLSUMS.
*
*      UPLO = 'U'                        UPLO = 'L'
*      ____i______                       ___________
*     |\   |      |                     |\          |
*     | \  |      |                     | \         |
*     |  \ |      |                     |  \        |
*     |   \|------| i                  i|---\       |
*     |    \      |                     |   |\      |
*     |     \     |                     |   | \     |
*     |      \    |                     |   |  \    |
*     |       \   |                     |   |   \   |
*     |        \  |                     |   |    \  |
*     |         \ |                     |   |     \ |
*     |__________\|                     |___|______\|
*                                           i
*
*     II, JJ  : local indices into array A
*     ICURROW : process row containing diagonal block
*     ICURCOL : process column containing diagonal block
*     IRSC0   : pointer to part of work used to store the ROWSUMS while
*               they are stored along a process column
*     IRSR0   : pointer to part of work used to store the ROWSUMS after
*               they have been transposed to be along a process row
*
      II = IIA
      JJ = JJA
*
      IF( N.EQ.0 ) THEN
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
         VALUE = ZERO
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Handle first block separately
*
            IB = IN-IA+1
*
*           Find COLMAXS
*
            IF( MYCOL.EQ.IACOL ) THEN
               DO 20 K = (JJ-1)*LDA, (JJ+IB-2)*LDA, LDA
                  IF( II.GT.IIA ) THEN
                     DO 10 LL = IIA, II-1
                        VALUE = MAX( VALUE, ABS( A( LL+K ) ) )
   10                CONTINUE
                  END IF
                  IF( MYROW.EQ.IAROW )
     $               II = II + 1
   20          CONTINUE
*
*              Reset local indices so we can find ROWMAXS
*
               IF( MYROW.EQ.IAROW )
     $            II = II - IB
*
            END IF
*
*           Find ROWMAXS
*
            IF( MYROW.EQ.IAROW ) THEN
               DO 40 K = II, II+IB-1
                  IF( JJ.LE.JJA+NQ-1 ) THEN
                     DO 30 LL = (JJ-1)*LDA, (JJA+NQ-2)*LDA, LDA
                        VALUE = MAX( VALUE, ABS( A( K+LL ) ) )
   30                CONTINUE
                  END IF
                  IF( MYCOL.EQ.IACOL )
     $               JJ = JJ + 1
   40          CONTINUE
               II = II + IB
            ELSE IF( MYCOL.EQ.IACOL ) THEN
               JJ = JJ + IB
            END IF
*
            ICURROW = MOD( IAROW+1, NPROW )
            ICURCOL = MOD( IACOL+1, NPCOL )
*
*           Loop over the remaining rows/columns of the matrix.
*
            DO 90 I = IN+1, IA+N-1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IA+N-I )
*
*              Find COLMAXS
*
               IF( MYCOL.EQ.ICURCOL ) THEN
                  DO 60 K = (JJ-1)*LDA, (JJ+IB-2)*LDA, LDA
                     IF( II.GT.IIA ) THEN
                        DO 50 LL = IIA, II-1
                           VALUE = MAX( VALUE, ABS( A( LL+K ) ) )
   50                   CONTINUE
                     END IF
                     IF( MYROW.EQ.ICURROW )
     $                  II = II + 1
   60             CONTINUE
*
*                 Reset local indices so we can find ROWMAXS
*
                  IF( MYROW.EQ.ICURROW )
     $               II = II - IB
               END IF
*
*              Find ROWMAXS
*
               IF( MYROW.EQ.ICURROW ) THEN
                  DO 80 K = II, II+IB-1
                     IF( JJ.LE.JJA+NQ-1 ) THEN
                        DO 70 LL = (JJ-1)*LDA, (JJA+NQ-2)*LDA, LDA
                           VALUE = MAX( VALUE, ABS( A( K+LL ) ) )
   70                   CONTINUE
                     END IF
                     IF( MYCOL.EQ.ICURCOL )
     $                  JJ = JJ + 1
   80             CONTINUE
                  II = II + IB
               ELSE IF( MYCOL.EQ.ICURCOL ) THEN
                  JJ = JJ + IB
               END IF
               ICURROW = MOD( ICURROW+1, NPROW )
               ICURCOL = MOD( ICURCOL+1, NPCOL )
   90       CONTINUE
*
         ELSE
*
*           Handle first block separately
*
            IB = IN-IA+1
*
*           Find COLMAXS
*
            IF( MYCOL.EQ.IACOL ) THEN
               DO 110 K = (JJ-1)*LDA, (JJ+IB-2)*LDA, LDA
                  IF( II.LE.IIA+NP-1 ) THEN
                     DO 100 LL = II, IIA+NP-1
                        VALUE = MAX( VALUE, ABS( A( LL+K ) ) )
  100                CONTINUE
                  END IF
                  IF( MYROW.EQ.IAROW )
     $               II = II + 1
  110          CONTINUE
*
*              Reset local indices so we can find ROWMAXS
*
               IF( MYROW.EQ.IAROW )
     $            II = II - IB
            END IF
*
*           Find ROWMAXS
*
            IF( MYROW.EQ.IAROW ) THEN
               DO 130 K = 0, IB-1
                  IF( JJ.GT.JJA ) THEN
                     DO 120 LL = (JJA-1)*LDA, (JJ-2)*LDA, LDA
                        VALUE = MAX( VALUE, ABS( A( II+LL ) ) )
  120                CONTINUE
                  END IF
                  II = II + 1
                  IF( MYCOL.EQ.IACOL )
     $               JJ = JJ + 1
  130          CONTINUE
            ELSE IF( MYCOL.EQ.IACOL ) THEN
               JJ = JJ + IB
            END IF
*
            ICURROW = MOD( IAROW+1, NPROW )
            ICURCOL = MOD( IACOL+1, NPCOL )
*
*           Loop over rows/columns of global matrix.
*
            DO 180 I = IN+1, IA+N-1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IA+N-I )
*
*              Find COLMAXS
*
               IF( MYCOL.EQ.ICURCOL ) THEN
                  DO 150 K = (JJ-1)*LDA, (JJ+IB-2)*LDA, LDA
                     IF( II.LE.IIA+NP-1 ) THEN
                        DO 140 LL = II, IIA+NP-1
                           VALUE = MAX( VALUE, ABS( A( LL+K ) ) )
  140                   CONTINUE
                     END IF
                      IF( MYROW.EQ.ICURROW )
     $                   II = II + 1
  150             CONTINUE
*
*                 Reset local indices so we can find ROWMAXS
*
                  IF( MYROW.EQ.ICURROW )
     $               II = II - IB
               END IF
*
*              Find ROWMAXS
*
               IF( MYROW.EQ.ICURROW ) THEN
                  DO 170 K = 0, IB-1
                     IF( JJ.GT.JJA ) THEN
                        DO 160 LL = (JJA-1)*LDA, (JJ-2)*LDA, LDA
                           VALUE = MAX( VALUE, ABS( A( II+LL ) ) )
  160                   CONTINUE
                     END IF
                     II = II + 1
                     IF( MYCOL.EQ.ICURCOL )
     $                  JJ = JJ + 1
  170             CONTINUE
               ELSE IF( MYCOL.EQ.ICURCOL ) THEN
                  JJ = JJ + IB
               END IF
               ICURROW = MOD( ICURROW+1, NPROW )
               ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  180       CONTINUE
*
         END IF
*
*        Gather the result on process (IAROW,IACOL).
*
         CALL SGAMX2D( ICTXT, 'All', ' ', 1, 1, VALUE, 1, I, K, -1,
     $                 IAROW, IACOL )
*
************************************************************************
* one or inf norm
*
      ELSE IF( LSAME( NORM, 'I' ) .OR. LSAME( NORM, 'O' ) .OR.
     $         NORM.EQ.'1' ) THEN
*
*        Find normI( sub( A ) ) ( = norm1( sub( A ) ), since sub( A ) is
*        symmetric).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Handle first block separately
*
            IB = IN-IA+1
*
*           Find COLSUMS
*
            IF( MYCOL.EQ.IACOL ) THEN
               IOFFA = ( JJ - 1 ) * LDA
               DO 200 K = 0, IB-1
                  SUM = ZERO
                  IF( II.GT.IIA ) THEN
                     DO 190 LL = IIA, II-1
                        SUM = SUM + ABS( A( LL+IOFFA ) )
  190                CONTINUE
                  END IF
                  IOFFA = IOFFA + LDA
                  WORK( JJ+K-JJA+ICSR0 ) = SUM
                  IF( MYROW.EQ.IAROW )
     $               II = II + 1
  200          CONTINUE
*
*              Reset local indices so we can find ROWSUMS
*
               IF( MYROW.EQ.IAROW )
     $            II = II - IB
*
            END IF
*
*           Find ROWSUMS
*
            IF( MYROW.EQ.IAROW ) THEN
               DO 220 K = II, II+IB-1
                  SUM = ZERO
                  IF( JJA+NQ.GT.JJ ) THEN
                     DO 210 LL = (JJ-1)*LDA, (JJA+NQ-2)*LDA, LDA
                        SUM = SUM + ABS( A( K+LL ) )
  210                CONTINUE
                  END IF
                  WORK( K-IIA+IRSC0 ) = SUM
                  IF( MYCOL.EQ.IACOL )
     $               JJ = JJ + 1
  220          CONTINUE
               II = II + IB
            ELSE IF( MYCOL.EQ.IACOL ) THEN
               JJ = JJ + IB
            END IF
*
            ICURROW = MOD( IAROW+1, NPROW )
            ICURCOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining rows/columns of global matrix.
*
            DO 270 I = IN+1, IA+N-1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IA+N-I )
*
*              Find COLSUMS
*
               IF( MYCOL.EQ.ICURCOL ) THEN
                  IOFFA = ( JJ - 1 ) * LDA
                  DO 240 K = 0, IB-1
                     SUM = ZERO
                     IF( II.GT.IIA ) THEN
                        DO 230 LL = IIA, II-1
                           SUM = SUM + ABS( A( IOFFA+LL ) )
  230                   CONTINUE
                     END IF
                     IOFFA = IOFFA + LDA
                     WORK( JJ+K-JJA+ICSR0 ) = SUM
                     IF( MYROW.EQ.ICURROW )
     $                  II = II + 1
  240             CONTINUE
*
*                 Reset local indices so we can find ROWSUMS
*
                  IF( MYROW.EQ.ICURROW )
     $               II = II - IB
*
               END IF
*
*              Find ROWSUMS
*
               IF( MYROW.EQ.ICURROW ) THEN
                  DO 260 K = II, II+IB-1
                     SUM = ZERO
                     IF( JJA+NQ.GT.JJ ) THEN
                        DO 250 LL = (JJ-1)*LDA, (JJA+NQ-2)*LDA, LDA
                           SUM = SUM + ABS( A( K+LL ) )
  250                   CONTINUE
                     END IF
                     WORK( K-IIA+IRSC0 ) = SUM
                     IF( MYCOL.EQ.ICURCOL )
     $                  JJ = JJ + 1
  260             CONTINUE
                  II = II + IB
               ELSE IF( MYCOL.EQ.ICURCOL ) THEN
                  JJ = JJ + IB
               END IF
*
               ICURROW = MOD( ICURROW+1, NPROW )
               ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  270       CONTINUE
*
         ELSE
*
*           Handle first block separately
*
            IB = IN-IA+1
*
*           Find COLSUMS
*
            IF( MYCOL.EQ.IACOL ) THEN
               IOFFA = (JJ-1)*LDA
               DO 290 K = 0, IB-1
                  SUM = ZERO
                  IF( IIA+NP.GT.II ) THEN
                     DO 280 LL = II, IIA+NP-1
                        SUM = SUM + ABS( A( IOFFA+LL ) )
  280                CONTINUE
                  END IF
                  IOFFA = IOFFA + LDA
                  WORK( JJ+K-JJA+ICSR0 ) = SUM
                  IF( MYROW.EQ.IAROW )
     $               II = II + 1
  290          CONTINUE
*
*              Reset local indices so we can find ROWSUMS
*
               IF( MYROW.EQ.IAROW )
     $            II = II - IB
*
            END IF
*
*           Find ROWSUMS
*
            IF( MYROW.EQ.IAROW ) THEN
               DO 310 K = II, II+IB-1
                  SUM = ZERO
                  IF( JJ.GT.JJA ) THEN
                     DO 300 LL = (JJA-1)*LDA, (JJ-2)*LDA, LDA
                        SUM = SUM + ABS( A( K+LL ) )
  300                CONTINUE
                  END IF
                  WORK( K-IIA+IRSC0 ) = SUM
                  IF( MYCOL.EQ.IACOL )
     $               JJ = JJ + 1
  310          CONTINUE
               II = II + IB
            ELSE IF( MYCOL.EQ.IACOL ) THEN
               JJ = JJ + IB
            END IF
*
            ICURROW = MOD( IAROW+1, NPROW )
            ICURCOL = MOD( IACOL+1, NPCOL )
*
*           Loop over rows/columns of global matrix.
*
            DO 360 I = IN+1, IA+N-1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IA+N-I )
*
*              Find COLSUMS
*
               IF( MYCOL.EQ.ICURCOL ) THEN
                  IOFFA = ( JJ - 1 ) * LDA
                  DO 330 K = 0, IB-1
                     SUM = ZERO
                     IF( IIA+NP.GT.II ) THEN
                        DO 320 LL = II, IIA+NP-1
                           SUM = SUM + ABS( A( LL+IOFFA ) )
  320                   CONTINUE
                     END IF
                     IOFFA = IOFFA + LDA
                     WORK( JJ+K-JJA+ICSR0 ) = SUM
                     IF( MYROW.EQ.ICURROW )
     $                  II = II + 1
  330             CONTINUE
*
*                 Reset local indices so we can find ROWSUMS
*
                  IF( MYROW.EQ.ICURROW )
     $               II = II - IB
*
               END IF
*
*              Find ROWSUMS
*
               IF( MYROW.EQ.ICURROW ) THEN
                  DO 350 K = II, II+IB-1
                     SUM = ZERO
                     IF( JJ.GT.JJA ) THEN
                        DO 340 LL = (JJA-1)*LDA, (JJ-2)*LDA, LDA
                           SUM = SUM + ABS( A( K+LL ) )
  340                   CONTINUE
                     END IF
                     WORK(K-IIA+IRSC0) = SUM
                     IF( MYCOL.EQ.ICURCOL )
     $                  JJ = JJ + 1
  350             CONTINUE
                  II = II + IB
               ELSE IF( MYCOL.EQ.ICURCOL ) THEN
                  JJ = JJ + IB
               END IF
*
               ICURROW = MOD( ICURROW+1, NPROW )
               ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  360       CONTINUE
         END IF
*
*        After calls to SGSUM2D, process row 0 will have global
*        COLSUMS and process column 0 will have global ROWSUMS.
*        Transpose ROWSUMS and add to COLSUMS to get global row/column
*        sum, the max of which is the infinity or 1 norm.
*
         IF( MYCOL.EQ.IACOL )
     $      NQ = NQ + ICOFF
         CALL SGSUM2D( ICTXT, 'Columnwise', ' ', 1, NQ, WORK( ICSR ), 1,
     $                 IAROW, MYCOL )
         IF( MYROW.EQ.IAROW )
     $      NP = NP + IROFF
         CALL SGSUM2D( ICTXT, 'Rowwise', ' ', NP, 1, WORK( IRSC ),
     $                 MAX( 1, NP ), MYROW, IACOL )
*
         CALL PSCOL2ROW( ICTXT, N, 1, DESCA( MB_ ), WORK( IRSC ),
     $                   MAX( 1, NP ), WORK( IRSR ), MAX( 1, NQ ),
     $                   IAROW, IACOL, IAROW, IACOL, WORK( IRSC+NP ) )
*
         IF( MYROW.EQ.IAROW ) THEN
            IF( MYCOL.EQ.IACOL )
     $         NQ = NQ - ICOFF
            CALL SAXPY( NQ, ONE, WORK( IRSR0 ), 1, WORK( ICSR0 ), 1 )
            IF( NQ.LT.1 ) THEN
               VALUE = ZERO
            ELSE
               VALUE = WORK( ISAMAX( NQ, WORK( ICSR0 ), 1 ) )
            END IF
            CALL SGAMX2D( ICTXT, 'Rowwise', ' ', 1, 1, VALUE, 1, I, K,
     $                    -1, IAROW, IACOL )
         END IF
*
************************************************************************
* Frobenius norm
* SSQ(1) is scale
* SSQ(2) is sum-of-squares
*
      ELSE IF( LSAME( NORM, 'F' ) .OR. LSAME( NORM, 'E' ) ) THEN
*
*        Find normF( sub( A ) ).
*
         SSQ(1) = ZERO
         SSQ(2) = ONE
*
*        Add off-diagonal entries, first
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Handle first block separately
*
            IB = IN-IA+1
*
            IF( MYCOL.EQ.IACOL ) THEN
               DO 370 K = (JJ-1)*LDA, (JJ+IB-2)*LDA, LDA
                  COLSSQ(1) = ZERO
                  COLSSQ(2) = ONE
                  CALL CLASSQ( II-IIA, A( IIA+K ), 1,
     $                         COLSSQ(1), COLSSQ(2) )
                  IF( MYROW.EQ.IAROW )
     $               II = II + 1
                  CALL CLASSQ( II-IIA, A( IIA+K ), 1,
     $                         COLSSQ(1), COLSSQ(2) )
                  CALL SCOMBSSQ( SSQ, COLSSQ )
  370          CONTINUE
*
               JJ = JJ + IB
            ELSE IF( MYROW.EQ.IAROW ) THEN
               II = II + IB
            END IF
*
            ICURROW = MOD( IAROW+1, NPROW )
            ICURCOL = MOD( IACOL+1, NPCOL )
*
*           Loop over rows/columns of global matrix.
*
            DO 390 I = IN+1, IA+N-1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IA+N-I )
*
               IF( MYCOL.EQ.ICURCOL ) THEN
                  DO 380 K = (JJ-1)*LDA, (JJ+IB-2)*LDA, LDA
                     COLSSQ(1) = ZERO
                     COLSSQ(2) = ONE
                     CALL CLASSQ( II-IIA, A( IIA+K ), 1,
     $                            COLSSQ(1), COLSSQ(2) )
                     IF( MYROW.EQ.ICURROW )
     $                  II = II + 1
                     CALL CLASSQ( II-IIA, A (IIA+K ), 1,
     $                            COLSSQ(1), COLSSQ(2) )
                     CALL SCOMBSSQ( SSQ, COLSSQ )
  380             CONTINUE
*
                  JJ = JJ + IB
               ELSE IF( MYROW.EQ.ICURROW ) THEN
                  II = II + IB
               END IF
*
               ICURROW = MOD( ICURROW+1, NPROW )
               ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  390       CONTINUE
*
         ELSE
*
*           Handle first block separately
*
            IB = IN-IA+1
*
            IF( MYCOL.EQ.IACOL ) THEN
               DO 400 K = (JJ-1)*LDA, (JJ+IB-2)*LDA, LDA
                  COLSSQ(1) = ZERO
                  COLSSQ(2) = ONE
                  CALL CLASSQ( IIA+NP-II, A( II+K ), 1,
     $                         COLSSQ(1), COLSSQ(2) )
                  IF( MYROW.EQ.IAROW )
     $               II = II + 1
                  CALL CLASSQ( IIA+NP-II, A( II+K ), 1,
     $                         COLSSQ(1), COLSSQ(2) )
                  CALL SCOMBSSQ( SSQ, COLSSQ )
  400          CONTINUE
*
               JJ = JJ + IB
            ELSE IF( MYROW.EQ.IAROW ) THEN
               II = II + IB
            END IF
*
            ICURROW = MOD( IAROW+1, NPROW )
            ICURCOL = MOD( IACOL+1, NPCOL )
*
*           Loop over rows/columns of global matrix.
*
            DO 420 I = IN+1, IA+N-1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IA+N-I )
*
               IF( MYCOL.EQ.ICURCOL ) THEN
                  DO 410 K = (JJ-1)*LDA, (JJ+IB-2)*LDA, LDA
                     COLSSQ(1) = ZERO
                     COLSSQ(2) = ONE
                     CALL CLASSQ( IIA+NP-II, A( II+K ), 1,
     $                            COLSSQ(1), COLSSQ(2) )
                     IF( MYROW.EQ.ICURROW )
     $                  II = II + 1
                     CALL CLASSQ( IIA+NP-II, A( II+K ), 1,
     $                            COLSSQ(1), COLSSQ(2) )
                     CALL SCOMBSSQ( SSQ, COLSSQ )
  410             CONTINUE
*
                  JJ = JJ + IB
               ELSE IF( MYROW.EQ.ICURROW ) THEN
                  II = II + IB
               END IF
*
               ICURROW = MOD( ICURROW+1, NPROW )
               ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  420       CONTINUE
*
         END IF
*
*        Perform the global scaled sum
*
         CALL PSTREECOMB( ICTXT, 'All', 2, SSQ, IAROW, IACOL,
     $                    SCOMBSSQ )
         VALUE = SSQ( 1 ) * SQRT( SSQ( 2 ) )
*
      END IF
*
*     Broadcast the result to the other processes
*
      IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL ) THEN
          CALL SGEBS2D( ICTXT, 'All', ' ', 1, 1, VALUE, 1 )
      ELSE
          CALL SGEBR2D( ICTXT, 'All', ' ', 1, 1, VALUE, 1, IAROW,
     $                  IACOL )
      END IF
*
      PCLANSY = VALUE
*
      RETURN
*
*     End of PCLANSY
*
      END
