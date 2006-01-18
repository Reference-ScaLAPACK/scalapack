      SUBROUTINE PZLAQGE( M, N, A, IA, JA, DESCA, R, C, ROWCND, COLCND,
     $                    AMAX, EQUED )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED
      INTEGER            IA, JA, M, N
      DOUBLE PRECISION   AMAX, COLCND, ROWCND
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   C( * ), R( * )
      COMPLEX*16         A( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLAQGE equilibrates a general M-by-N distributed matrix
*  sub( A ) = A(IA:IA+M-1,JA:JA+N-1) using the row and scaling
*  factors in the vectors R and C.
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
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1))
*          containing on entry the M-by-N matrix sub( A ). On exit,
*          the equilibrated distributed matrix.  See EQUED for the
*          form of the equilibrated distributed submatrix.
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
*  R       (local input) DOUBLE PRECISION array, dimension LOCr(M_A)
*          The row scale factors for sub( A ). R is aligned with the
*          distributed matrix A, and replicated across every process
*          column. R is tied to the distributed matrix A.
*
*  C       (local input) DOUBLE PRECISION array, dimension LOCc(N_A)
*          The column scale factors of sub( A ). C is aligned with the
*          distributed matrix A, and replicated down every process
*          row. C is tied to the distributed matrix A.
*
*  ROWCND  (global input) DOUBLE PRECISION
*          The global ratio of the smallest R(i) to the largest R(i),
*          IA <= i <= IA+M-1.
*
*  COLCND  (global input) DOUBLE PRECISION
*          The global ratio of the smallest C(i) to the largest C(i),
*          JA <= j <= JA+N-1.
*
*  AMAX    (global input) DOUBLE PRECISION
*          Absolute value of largest distributed submatrix entry.
*
*  EQUED   (global output) CHARACTER
*          Specifies the form of equilibration that was done.
*          = 'N':  No equilibration
*          = 'R':  Row equilibration, i.e., sub( A ) has been pre-
*                  multiplied by diag(R(IA:IA+M-1)),
*          = 'C':  Column equilibration, i.e., sub( A ) has been post-
*                  multiplied by diag(C(JA:JA+N-1)),
*          = 'B':  Both row and column equilibration, i.e., sub( A )
*                  has been replaced by
*                  diag(R(IA:IA+M-1)) * sub( A ) * diag(C(JA:JA+N-1)).
*
*  Internal Parameters
*  ===================
*
*  THRESH is a threshold value used to decide if row or column scaling
*  should be done based on the ratio of the row or column scaling
*  factors.  If ROWCND < THRESH, row scaling is done, and if
*  COLCND < THRESH, column scaling is done.
*
*  LARGE and SMALL are threshold values used to decide if row scaling
*  should be done based on the absolute size of the largest matrix
*  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE, THRESH
      PARAMETER          ( ONE = 1.0D+0, THRESH = 0.1D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IACOL, IAROW, ICOFF, ICTXT, IIA, IOFFA,
     $                   IROFF, J, JJA, LDA, MP, MYCOL, MYROW, NPCOL,
     $                   NPROW, NQ
      DOUBLE PRECISION   CJ, LARGE, SMALL
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           NUMROC, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF
*
*     Get grid parameters and compute local indexes
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
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
*
*     Initialize LARGE and SMALL.
*
      SMALL = PDLAMCH( ICTXT, 'Safe minimum' ) /
     $        PDLAMCH( ICTXT, 'Precision' )
      LARGE = ONE / SMALL
*
      IF( ROWCND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE )
     $     THEN
*
*        No row scaling
*
         IF( COLCND.GE.THRESH ) THEN
*
*           No column scaling
*
            EQUED = 'N'
*
         ELSE
*
*           Column scaling
*
            IOFFA = (JJA-1)*LDA
            DO 20 J = JJA, JJA+NQ-1
               CJ = C( J )
               DO 10 I = IIA, IIA+MP-1
                  A( IOFFA + I ) = CJ*A( IOFFA + I )
   10          CONTINUE
               IOFFA = IOFFA + LDA
   20       CONTINUE
            EQUED = 'C'
         END IF
*
      ELSE IF( COLCND.GE.THRESH ) THEN
*
*        Row scaling, no column scaling
*
         IOFFA = (JJA-1)*LDA
         DO 40 J = JJA, JJA+NQ-1
            DO 30 I = IIA, IIA+MP-1
               A( IOFFA + I ) = R( I )*A( IOFFA + I )
   30       CONTINUE
            IOFFA = IOFFA + LDA
   40    CONTINUE
         EQUED = 'R'
*
      ELSE
*
*        Row and column scaling
*
         IOFFA = (JJA-1)*LDA
         DO 60 J = JJA, JJA+NQ-1
            CJ = C( J )
            DO 50 I = IIA, IIA+MP-1
               A( IOFFA + I ) = CJ*R( I )*A( IOFFA + I )
   50       CONTINUE
            IOFFA = IOFFA + LDA
   60    CONTINUE
         EQUED = 'B'
*
      END IF
*
      RETURN
*
*     End of PZLAQGE
*
      END
