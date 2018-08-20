      REAL               FUNCTION PSLANGE( NORM, M, N, A, IA, JA, DESCA,
     $                                     WORK )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            IA, JA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      REAL               A( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLANGE returns the value of the one norm, or the Frobenius norm,
*  or the infinity norm, or the element of largest absolute value of a
*  distributed matrix sub( A ) = A(IA:IA+M-1, JA:JA+N-1).
*
*  PSLANGE returns the value
*
*     ( max(abs(A(i,j))),  NORM = 'M' or 'm' with IA <= i <= IA+M-1,
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
*          Specifies the value to be returned in PSLANGE as described
*          above.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( A ). When M = 0, PSLANGE
*          is set to zero. M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( A ). When N = 0,
*          PSLANGE is set to zero. N >= 0.
*
*  A       (local input) REAL pointer into the local memory
*          to an array of dimension (LLD_A, LOCc(JA+N-1)) containing the
*          local pieces of the distributed matrix sub( A ).
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
      INTEGER            I, IACOL, IAROW, ICTXT, II, ICOFF, IOFFA,
     $                   IROFF, J, JJ, LDA, MP, MYCOL, MYROW, NPCOL,
     $                   NPROW, NQ
      REAL               SUM, VALUE
*     ..
*     .. Local Arrays ..
      REAL               SSQ( 2 ), COLSSQ( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, PSTREECOMB,
     $                   SCOMBSSQ, SGEBR2D, SGEBS2D,
     $                   SGAMX2D, SGSUM2D, SLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ISAMAX, NUMROC
      EXTERNAL           LSAME, ISAMAX, NUMROC
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
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, II, JJ,
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
         VALUE = ZERO
         IF( NQ.GT.0 .AND. MP.GT.0 ) THEN
            IOFFA = (JJ-1)*LDA
            DO 20 J = JJ, JJ+NQ-1
               DO 10 I = II, MP+II-1
                  VALUE = MAX( VALUE, ABS( A( IOFFA+I ) ) )
   10          CONTINUE
               IOFFA = IOFFA + LDA
   20       CONTINUE
         END IF
         CALL SGAMX2D( ICTXT, 'All', ' ', 1, 1, VALUE, 1, I, J, -1,
     $                 0, 0 )
*
************************************************************************
* one norm
*
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' ) THEN
*
*        Find norm1( sub( A ) ).
*
         IF( NQ.GT.0 ) THEN
            IOFFA = ( JJ - 1 ) * LDA
            DO 40 J = JJ, JJ+NQ-1
               SUM = ZERO
               IF( MP.GT.0 ) THEN
                  DO 30 I = II, MP+II-1
                     SUM = SUM + ABS( A( IOFFA+I ) )
   30             CONTINUE
               END IF
               IOFFA = IOFFA + LDA
               WORK( J-JJ+1 ) = SUM
   40       CONTINUE
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
            CALL SGAMX2D( ICTXT, 'Rowwise', ' ', 1, 1, VALUE, 1, I, J,
     $                    -1, 0, 0 )
         END IF
*
************************************************************************
* inf norm
*
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
*        Find normI( sub( A ) ).
*
         IF( MP.GT.0 ) THEN
            IOFFA = II + ( JJ - 1 ) * LDA
            DO 60 I = II, II+MP-1
               SUM = ZERO
               IF( NQ.GT.0 ) THEN
                  DO 50 J = IOFFA, IOFFA + NQ*LDA - 1, LDA
                     SUM = SUM + ABS( A( J ) )
   50             CONTINUE
               END IF
               WORK( I-II+1 ) = SUM
               IOFFA = IOFFA + 1
   60       CONTINUE
         END IF
*
*        Find sum of global matrix rows and store on column 0 of
*        process grid
*
         CALL SGSUM2D( ICTXT, 'Rowwise', ' ', MP, 1, WORK, MAX( 1, MP ),
     $                 MYROW, 0 )
*
*        Find maximum sum of rows for supnorm
*
         IF( MYCOL.EQ.0 ) THEN
            IF( MP.GT.0 ) THEN
               VALUE = WORK( ISAMAX( MP, WORK, 1 ) )
            ELSE
               VALUE = ZERO
            END IF
            CALL SGAMX2D( ICTXT, 'Columnwise', ' ', 1, 1, VALUE, 1, I,
     $                    J, -1, 0, 0 )
         END IF
*
************************************************************************
* Frobenius norm
* SSQ(1) is scale
* SSQ(2) is sum-of-squares
*
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF( sub( A ) ).
*
         SSQ(1) = ZERO
         SSQ(2) = ONE
         IOFFA = II + ( JJ - 1 ) * LDA
         IF( NQ.GT.0 ) THEN
             DO 70 J = IOFFA, IOFFA + NQ*LDA - 1, LDA
                COLSSQ(1) = ZERO
                COLSSQ(2) = ONE
                CALL SLASSQ( MP, A( J ), 1, COLSSQ(1), COLSSQ(2) )
                CALL SCOMBSSQ( SSQ, COLSSQ )
   70        CONTINUE
         END IF
*
*        Perform the global scaled sum
*
         CALL PSTREECOMB( ICTXT, 'All', 2, SSQ, 0, 0, SCOMBSSQ )
         VALUE = SSQ( 1 ) * SQRT( SSQ( 2 ) )
*
      END IF
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         CALL SGEBS2D( ICTXT, 'All', ' ', 1, 1, VALUE, 1 )
      ELSE
         CALL SGEBR2D( ICTXT, 'All', ' ', 1, 1, VALUE, 1, 0, 0 )
      END IF
*
      PSLANGE = VALUE
*
      RETURN
*
*     End of PSLANGE
*
      END
