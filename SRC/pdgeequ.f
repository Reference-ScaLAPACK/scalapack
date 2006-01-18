      SUBROUTINE PDGEEQU( M, N, A, IA, JA, DESCA, R, C, ROWCND, COLCND,
     $                    AMAX, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, M, N
      DOUBLE PRECISION   AMAX, COLCND, ROWCND
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), C( * ), R( * )
*     ..
*
*  Purpose
*  =======
*
*  PDGEEQU computes row and column scalings intended to equilibrate an
*  M-by-N distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA:JA+N-1) and
*  reduce its condition number.  R returns the row scale factors and C
*  the column scale factors, chosen to try to make the largest entry in
*  each row and column of the distributed matrix B with elements
*  B(i,j) = R(i) * A(i,j) * C(j) have absolute value 1.
*
*  R(i) and C(j) are restricted to be between SMLNUM = smallest safe
*  number and BIGNUM = largest safe number.  Use of these scaling
*  factors is not guaranteed to reduce the condition number of
*  sub( A ) but works well in practice.
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
*  A       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension ( LLD_A, LOCc(JA+N-1) ), the
*          local pieces of the M-by-N distributed matrix whose
*          equilibration factors are to be computed.
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
*  R       (local output) DOUBLE PRECISION array, dimension LOCr(M_A)
*          If INFO = 0 or INFO > IA+M-1, R(IA:IA+M-1) contains the row
*          scale factors for sub( A ). R is aligned with the distributed
*          matrix A, and replicated across every process column. R is
*          tied to the distributed matrix A.
*
*  C       (local output) DOUBLE PRECISION array, dimension LOCc(N_A)
*          If INFO = 0,  C(JA:JA+N-1) contains the column scale factors
*          for sub( A ). C is aligned with the distributed matrix A, and
*          replicated down every process row. C is tied to the distri-
*          buted matrix A.
*
*  ROWCND  (global output) DOUBLE PRECISION
*          If INFO = 0 or INFO > IA+M-1, ROWCND contains the ratio of
*          the smallest R(i) to the largest R(i) (IA <= i <= IA+M-1).
*          If ROWCND >= 0.1 and AMAX is neither too large nor too small,
*          it is not worth scaling by R(IA:IA+M-1).
*
*  COLCND  (global output) DOUBLE PRECISION
*          If INFO = 0, COLCND contains the ratio of the smallest C(j)
*          to the largest C(j) (JA <= j <= JA+N-1). If COLCND >= 0.1, it
*          is not worth scaling by C(JA:JA+N-1).
*
*  AMAX    (global output) DOUBLE PRECISION
*          Absolute value of largest distributed matrix element.  If
*          AMAX is very close to overflow or very close to underflow,
*          the matrix should be scaled.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = i,  and i is
*                <= M:  the i-th row of the distributed matrix sub( A )
*                       is exactly zero,
*                >  M:  the (i-M)-th column of the distributed
*                       matrix sub( A ) is exactly zero.
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
      CHARACTER          COLCTOP, ROWCTOP
      INTEGER            I, IACOL, IAROW, ICOFF, ICTXT, IDUMM, IIA,
     $                   IOFFA, IROFF, J, JJA, LDA, MP, MYCOL, MYROW,
     $                   NPCOL, NPROW, NQ
      DOUBLE PRECISION   BIGNUM, RCMAX, RCMIN, SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER            DESCC( DLEN_ ), DESCR( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, DGAMN2D,
     $                   DGAMX2D, IGAMX2D, INFOG2L, PCHK1MAT, PB_TOPGET,
     $                   PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            INDXL2G, NUMROC
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           INDXL2G, NUMROC, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters.
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(600+CTXT_)
      ELSE
         CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, INFO )
         CALL PCHK1MAT( M, 1, N, 2, IA, JA, DESCA, 6, 0, IDUMM, IDUMM,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDGEEQU', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         ROWCND = ONE
         COLCND = ONE
         AMAX = ZERO
         RETURN
      END IF
*
      CALL PB_TOPGET( ICTXT, 'Combine', 'Rowwise', ROWCTOP )
      CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
*
*     Get machine constants and local indexes.
*
      SMLNUM = PDLAMCH( ICTXT, 'S' )
      BIGNUM = ONE / SMLNUM
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
*     Assign descriptors for R and C arrays
*
      CALL DESCSET( DESCR, M, 1, DESCA( MB_ ), 1, 0, 0, ICTXT,
     $               MAX( 1, MP ) )
      CALL DESCSET( DESCC, 1, N, 1, DESCA( NB_ ), 0, 0, ICTXT, 1 )
*
*     Compute row scale factors.
*
      DO 10 I = IIA, IIA+MP-1
         R( I ) = ZERO
   10 CONTINUE
*
*     Find the maximum element in each row.
*
      IOFFA = (JJA-1)*LDA
      DO 30 J = JJA, JJA+NQ-1
         DO 20 I = IIA, IIA+MP-1
            R( I ) = MAX( R( I ), ABS( A( IOFFA + I ) ) )
   20    CONTINUE
         IOFFA = IOFFA + LDA
   30 CONTINUE
      CALL DGAMX2D( ICTXT, 'Rowwise', ROWCTOP, MP, 1, R( IIA ),
     $              MAX( 1, MP ), IDUMM, IDUMM, -1, -1, MYCOL )
*
*     Find the maximum and minimum scale factors.
*
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 40 I = IIA, IIA+MP-1
         RCMAX = MAX( RCMAX, R( I ) )
         RCMIN = MIN( RCMIN, R( I ) )
   40 CONTINUE
      CALL DGAMX2D( ICTXT, 'Columnwise', COLCTOP, 1, 1, RCMAX, 1, IDUMM,
     $              IDUMM, -1, -1, MYCOL )
      CALL DGAMN2D( ICTXT, 'Columnwise', COLCTOP, 1, 1, RCMIN, 1, IDUMM,
     $              IDUMM, -1, -1, MYCOL )
      AMAX = RCMAX
*
      IF( RCMIN.EQ.ZERO ) THEN
*
*        Find the first zero scale factor and return an error code.
*
         DO 50 I = IIA, IIA+MP-1
            IF( R( I ).EQ.ZERO .AND. INFO.EQ.0 )
     $         INFO = INDXL2G( I, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                NPROW ) - IA + 1
   50    CONTINUE
         CALL IGAMX2D( ICTXT, 'Columnwise', COLCTOP, 1, 1, INFO, 1,
     $                 IDUMM, IDUMM, -1, -1, MYCOL )
         IF( INFO.NE.0 )
     $      RETURN
      ELSE
*
*        Invert the scale factors.
*
         DO 60 I = IIA, IIA+MP-1
            R( I ) = ONE / MIN( MAX( R( I ), SMLNUM ), BIGNUM )
   60    CONTINUE
*
*        Compute ROWCND = min(R(I)) / max(R(I))
*
         ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
*
      END IF
*
*     Compute column scale factors
*
      DO 70 J = JJA, JJA+NQ-1
         C( J ) = ZERO
   70 CONTINUE
*
*     Find the maximum element in each column,
*     assuming the row scaling computed above.
*
      IOFFA = (JJA-1)*LDA
      DO 90 J = JJA, JJA+NQ-1
         DO 80 I = IIA, IIA+MP-1
            C( J ) = MAX( C( J ), ABS( A( IOFFA + I ) )*R( I ) )
   80    CONTINUE
         IOFFA = IOFFA + LDA
   90 CONTINUE
      CALL DGAMX2D( ICTXT, 'Columnwise', COLCTOP, 1, NQ, C( JJA ),
     $              1, IDUMM, IDUMM, -1, -1, MYCOL )
*
*     Find the maximum and minimum scale factors.
*
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 100 J = JJA, JJA+NQ-1
         RCMIN = MIN( RCMIN, C( J ) )
         RCMAX = MAX( RCMAX, C( J ) )
  100 CONTINUE
      CALL DGAMX2D( ICTXT, 'Columnwise', COLCTOP, 1, 1, RCMAX, 1, IDUMM,
     $              IDUMM, -1, -1, MYCOL )
      CALL DGAMN2D( ICTXT, 'Columnwise', COLCTOP, 1, 1, RCMIN, 1, IDUMM,
     $              IDUMM, -1, -1, MYCOL )
*
      IF( RCMIN.EQ.ZERO ) THEN
*
*        Find the first zero scale factor and return an error code.
*
         DO 110 J = JJA, JJA+NQ-1
            IF( C( J ).EQ.ZERO .AND. INFO.EQ.0 )
     $         INFO = M + INDXL2G( J, DESCA( NB_ ), MYCOL,
     $                DESCA( CSRC_ ), NPCOL ) - JA + 1
  110    CONTINUE
         CALL IGAMX2D( ICTXT, 'Columnwise', COLCTOP, 1, 1, INFO, 1,
     $                 IDUMM, IDUMM, -1, -1, MYCOL )
         IF( INFO.NE.0 )
     $      RETURN
      ELSE
*
*        Invert the scale factors.
*
         DO 120 J = JJA, JJA+NQ-1
            C( J ) = ONE / MIN( MAX( C( J ), SMLNUM ), BIGNUM )
  120    CONTINUE
*
*        Compute COLCND = min(C(J)) / max(C(J))
*
         COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
*
      END IF
*
      RETURN
*
*     End of PDGEEQU
*
      END
