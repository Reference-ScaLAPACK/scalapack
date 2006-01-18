      SUBROUTINE PCLATTRS( UPLO, TRANS, DIAG, NORMIN, N, A, IA, JA,
     $                     DESCA, X, IX, JX, DESCX, SCALE, CNORM, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     July 31, 2001
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            IA, INFO, IX, JA, JX, N
      REAL               SCALE
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCX( * )
      REAL               CNORM( * )
      COMPLEX            A( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLATTRS solves one of the triangular systems
*
*     A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
*
*  with scaling to prevent overflow.  Here A is an upper or lower
*  triangular matrix, A**T denotes the transpose of A, A**H denotes the
*  conjugate transpose of A, x and b are n-element vectors, and s is a
*  scaling factor, usually less than or equal to 1, chosen so that the
*  components of x will be less than the overflow threshold.  If the
*  unscaled problem will not cause overflow, the Level 2 PBLAS routine
*  PCTRSV is called. If the matrix A is singular (A(j,j) = 0 for some j)
*  then s is set to 0 and a non-trivial solution to A*x = 0 is returned.
*
*  This is very slow relative to PCTRSV.  This should only be used
*  when scaling is necessary to control overflow, or when it is modified
*  to scale better.
*  Notes
*
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
*  and assume that its process grid has dimension r x c.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the r processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the c processes of
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
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  TRANS   (global input) CHARACTER*1
*          Specifies the operation applied to A.
*          = 'N':  Solve A * x = s*b     (No transpose)
*          = 'T':  Solve A**T * x = s*b  (Transpose)
*          = 'C':  Solve A**H * x = s*b  (Conjugate transpose)
*
*  DIAG    (global input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  NORMIN  (global input) CHARACTER*1
*          Specifies whether CNORM has been set or not.
*          = 'Y':  CNORM contains the column norms on entry
*          = 'N':  CNORM is not set on entry.  On exit, the norms will
*                  be computed and stored in CNORM.
*
*  N       (global input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (local input) COMPLEX array, dimension (DESCA(LLD_),*)
*          The triangular matrix A.  If UPLO = 'U', the leading n by n
*          upper triangular part of the array A contains the upper
*          triangular matrix, and the strictly lower triangular part of
*          A is not referenced.  If UPLO = 'L', the leading n by n lower
*          triangular part of the array A contains the lower triangular
*          matrix, and the strictly upper triangular part of A is not
*          referenced.  If DIAG = 'U', the diagonal elements of A are
*          also not referenced and are assumed to be 1.
*
*  IA      (global input) pointer to INTEGER
*          The global row index of the submatrix of the distributed
*          matrix A to operate on.
*
*  JA      (global input) pointer to INTEGER
*          The global column index of the submatrix of the distributed
*          matrix A to operate on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  X       (local input/output) COMPLEX array,
*                                             dimension (DESCX(LLD_),*)
*          On entry, the right hand side b of the triangular system.
*          On exit, X is overwritten by the solution vector x.
*
*  IX      (global input) pointer to INTEGER
*          The global row index of the submatrix of the distributed
*          matrix X to operate on.
*
*  JX      (global input) pointer to INTEGER
*          The global column index of the submatrix of the distributed
*          matrix X to operate on.
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*
*  SCALE   (global output) REAL
*          The scaling factor s for the triangular system
*             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.
*          If SCALE = 0, the matrix A is singular or badly scaled, and
*          the vector x is an exact or approximate solution to A*x = 0.
*
*  CNORM   (global input or global output) REAL array,
*                                                       dimension (N)
*          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
*          contains the norm of the off-diagonal part of the j-th column
*          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
*          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
*          must be greater than or equal to the 1-norm.
*
*          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
*          returns the 1-norm of the offdiagonal part of the j-th column
*          of A.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -k, the k-th argument had an illegal value
*
*  Further Details
*  ======= =======
*
*  A rough bound on x is computed; if that is less than overflow, PCTRSV
*  is called, otherwise, specific code is used which checks for possible
*  overflow or divide-by-zero at every operation.
*
*  A columnwise scheme is used for solving A*x = b.  The basic algorithm
*  if A is lower triangular is
*
*       x[1:n] := b[1:n]
*       for j = 1, ..., n
*            x(j) := x(j) / A(j,j)
*            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
*       end
*
*  Define bounds on the components of x after j iterations of the loop:
*     M(j) = bound on x[1:j]
*     G(j) = bound on x[j+1:n]
*  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
*
*  Then for iteration j+1 we have
*     M(j+1) <= G(j) / | A(j+1,j+1) |
*     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
*            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
*
*  where CNORM(j+1) is greater than or equal to the infinity-norm of
*  column j+1 of A, not counting the diagonal.  Hence
*
*     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
*                  1<=i<=j
*  and
*
*     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
*                                   1<=i< j
*
*  Since |x(j)| <= M(j), we use the Level 2 PBLAS routine PCTRSV if the
*  reciprocal of the largest M(j), j=1,..,n, is larger than
*  max(underflow, 1/overflow).
*
*  The bound on x(j) is also used to determine when a step in the
*  columnwise method can be performed without fear of overflow.  If
*  the computed bound is greater than a large constant, x is scaled to
*  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
*  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
*
*  Similarly, a row-wise scheme is used to solve A**T *x = b  or
*  A**H *x = b.  The basic algorithm for A upper triangular is
*
*       for j = 1, ..., n
*            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
*       end
*
*  We simultaneously compute two bounds
*       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
*       M(j) = bound on x(i), 1<=i<=j
*
*  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
*  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
*  Then the bound on x(j) is
*
*       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
*
*            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
*                      1<=i<=j
*
*  and we can safely call PCTRSV if 1/M(n) and 1/G(n) are both greater
*  than max(underflow, 1/overflow).
*
*  Last modified by: Mark R. Fahey, August 2000
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, HALF, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0, HALF = 0.5E+0, ONE = 1.0E+0,
     $                   TWO = 2.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            CONTXT, CSRC, I, ICOL, ICOLX, IMAX, IROW,
     $                   IROWX, ITMP1, ITMP1X, ITMP2, ITMP2X, J, JFIRST,
     $                   JINC, JLAST, LDA, LDX, MB, MYCOL, MYROW, NB,
     $                   NPCOL, NPROW, RSRC
      REAL               BIGNUM, GROW, REC, SMLNUM, TJJ, TMAX, TSCAL,
     $                   XBND, XJ, XMAX
      COMPLEX            CSUMJ, TJJS, USCAL, XJTMP, ZDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ISAMAX
      REAL               PSLAMCH
      COMPLEX            CLADIV
      EXTERNAL           LSAME, ISAMAX, PSLAMCH, CLADIV
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, SGSUM2D, SSCAL, INFOG2L,
     $                   PSCASUM, PSLABAD, PXERBLA, PCAMAX, PCAXPY,
     $                   PCDOTC, PCDOTU, PCSSCAL, PCLASET, PCSCAL,
     $                   PCTRSV, CGEBR2D, CGEBS2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL, CMPLX, CONJG, AIMAG, MAX, MIN
*     ..
*     .. Statement Functions ..
      REAL               CABS1, CABS2
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      CABS2( ZDUM ) = ABS( REAL( ZDUM ) / 2.E0 ) +
     $                ABS( AIMAG( ZDUM ) / 2.E0 )
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
*
      CONTXT = DESCA( CTXT_ )
      RSRC = DESCA( RSRC_ )
      CSRC = DESCA( CSRC_ )
      MB = DESCA( MB_ )
      NB = DESCA( NB_ )
      LDA = DESCA( LLD_ )
      LDX = DESCX( LLD_ )
*
*     Test the input parameters.
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT.
     $         LSAME( NORMIN, 'N' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      END IF
*
      CALL BLACS_GRIDINFO( CONTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( CONTXT, 'PCLATTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine machine dependent parameters to control overflow.
*
      SMLNUM = PSLAMCH( CONTXT, 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      CALL PSLABAD( CONTXT, SMLNUM, BIGNUM )
      SMLNUM = SMLNUM / PSLAMCH( CONTXT, 'Precision' )
      BIGNUM = ONE / SMLNUM
      SCALE = ONE
*
*
      IF( LSAME( NORMIN, 'N' ) ) THEN
*
*        Compute the 1-norm of each column, not including the diagonal.
*
         IF( UPPER ) THEN
*
*           A is upper triangular.
*
            CNORM( 1 ) = ZERO
            DO 10 J = 2, N
               CALL PSCASUM( J-1, CNORM( J ), A, IA, JA+J-1, DESCA, 1 )
   10       CONTINUE
         ELSE
*
*           A is lower triangular.
*
            DO 20 J = 1, N - 1
               CALL PSCASUM( N-J, CNORM( J ), A, IA+J, JA+J-1, DESCA,
     $                       1 )
   20       CONTINUE
            CNORM( N ) = ZERO
         END IF
         CALL SGSUM2D( CONTXT, 'Row', ' ', N, 1, CNORM, 1, -1, -1 )
      END IF
*
*     Scale the column norms by TSCAL if the maximum element in CNORM is
*     greater than BIGNUM/2.
*
      IMAX = ISAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      IF( TMAX.LE.BIGNUM*HALF ) THEN
         TSCAL = ONE
      ELSE
         TSCAL = HALF / ( SMLNUM*TMAX )
         CALL SSCAL( N, TSCAL, CNORM, 1 )
      END IF
*
*     Compute a bound on the computed solution vector to see if the
*     Level 2 PBLAS routine PCTRSV can be used.
*
      XMAX = ZERO
      CALL PCAMAX( N, ZDUM, IMAX, X, IX, JX, DESCX, 1 )
      XMAX = CABS2( ZDUM )
      CALL SGSUM2D( CONTXT, 'Row', ' ', 1, 1, XMAX, 1, -1, -1 )
      XBND = XMAX
*
      IF( NOTRAN ) THEN
*
*        Compute the growth in A * x = b.
*
         IF( UPPER ) THEN
            JFIRST = N
            JLAST = 1
            JINC = -1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
         END IF
*
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 50
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, G(0) = max{x(i), i=1,...,n}.
*
            GROW = HALF / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 30 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 50
*
*              TJJS = A( J, J )
               CALL INFOG2L( IA+J-1, JA+J-1, DESCA, NPROW, NPCOL, MYROW,
     $                       MYCOL, IROW, ICOL, ITMP1, ITMP2 )
               IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
                  TJJS = A( ( ICOL-1 )*LDA+IROW )
                  CALL CGEBS2D( CONTXT, 'All', ' ', 1, 1, TJJS, 1 )
               ELSE
                  CALL CGEBR2D( CONTXT, 'All', ' ', 1, 1, TJJS, 1,
     $                          ITMP1, ITMP2 )
               END IF
               TJJ = CABS1( TJJS )
*
               IF( TJJ.GE.SMLNUM ) THEN
*
*                 M(j) = G(j-1) / abs(A(j,j))
*
                  XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               ELSE
*
*                 M(j) could overflow, set XBND to 0.
*
                  XBND = ZERO
               END IF
*
               IF( TJJ+CNORM( J ).GE.SMLNUM ) THEN
*
*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
*
                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               ELSE
*
*                 G(j) could overflow, set GROW to 0.
*
                  GROW = ZERO
               END IF
   30       CONTINUE
            GROW = XBND
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) )
            DO 40 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 50
*
*              G(j) = G(j-1)*( 1 + CNORM(j) )
*
               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   40       CONTINUE
         END IF
   50    CONTINUE
*
      ELSE
*
*        Compute the growth in A**T * x = b  or  A**H * x = b.
*
         IF( UPPER ) THEN
            JFIRST = 1
            JLAST = N
            JINC = 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
         END IF
*
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 80
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, M(0) = max{x(i), i=1,...,n}.
*
            GROW = HALF / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 60 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 80
*
*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
*
               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )
*
*              TJJS = A( J, J )
               CALL INFOG2L( IA+J-1, JA+J-1, DESCA, NPROW, NPCOL, MYROW,
     $                       MYCOL, IROW, ICOL, ITMP1, ITMP2 )
               IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
                  TJJS = A( ( ICOL-1 )*LDA+IROW )
                  CALL CGEBS2D( CONTXT, 'All', ' ', 1, 1, TJJS, 1 )
               ELSE
                  CALL CGEBR2D( CONTXT, 'All', ' ', 1, 1, TJJS, 1,
     $                          ITMP1, ITMP2 )
               END IF
               TJJ = CABS1( TJJS )
*
               IF( TJJ.GE.SMLNUM ) THEN
*
*                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
*
                  IF( XJ.GT.TJJ )
     $               XBND = XBND*( TJJ / XJ )
               ELSE
*
*                 M(j) could overflow, set XBND to 0.
*
                  XBND = ZERO
               END IF
   60       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( ONE, HALF / MAX( XBND, SMLNUM ) )
            DO 70 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 80
*
*              G(j) = ( 1 + CNORM(j) )*G(j-1)
*
               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
   70       CONTINUE
         END IF
   80    CONTINUE
      END IF
*
      IF( ( GROW*TSCAL ).GT.SMLNUM ) THEN
*
*        Use the Level 2 PBLAS solve if the reciprocal of the bound on
*        elements of X is not too small.
*
         CALL PCTRSV( UPLO, TRANS, DIAG, N, A, IA, JA, DESCA, X, IX, JX,
     $                DESCX, 1 )
      ELSE
*
*        Use a Level 1 PBLAS solve, scaling intermediate results.
*
         IF( XMAX.GT.BIGNUM*HALF ) THEN
*
*           Scale X so that its components are less than or equal to
*           BIGNUM in absolute value.
*
            SCALE = ( BIGNUM*HALF ) / XMAX
            CALL PCSSCAL( N, SCALE, X, IX, JX, DESCX, 1 )
            XMAX = BIGNUM
         ELSE
            XMAX = XMAX*TWO
         END IF
*
         IF( NOTRAN ) THEN
*
*           Solve A * x = b
*
            DO 100 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
*
*              XJ = CABS1( X( J ) )
               CALL INFOG2L( IX+J-1, JX, DESCX, NPROW, NPCOL, MYROW,
     $                       MYCOL, IROWX, ICOLX, ITMP1X, ITMP2X )
               IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) ) THEN
                  XJTMP = X( IROWX )
                  CALL CGEBS2D( CONTXT, 'All', ' ', 1, 1, XJTMP, 1 )
               ELSE
                  CALL CGEBR2D( CONTXT, 'All', ' ', 1, 1, XJTMP, 1,
     $                          ITMP1X, ITMP2X )
               END IF
               XJ = CABS1( XJTMP )
               IF( NOUNIT ) THEN
*                 TJJS = A( J, J )*TSCAL
                  CALL INFOG2L( IA+J-1, JA+J-1, DESCA, NPROW, NPCOL,
     $                          MYROW, MYCOL, IROW, ICOL, ITMP1, ITMP2 )
                  IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
                     TJJS = A( ( ICOL-1 )*LDA+IROW )*TSCAL
                     CALL CGEBS2D( CONTXT, 'All', ' ', 1, 1, TJJS, 1 )
                  ELSE
                     CALL CGEBR2D( CONTXT, 'All', ' ', 1, 1, TJJS, 1,
     $                             ITMP1, ITMP2 )
                  END IF
               ELSE
                  TJJS = TSCAL
                  IF( TSCAL.EQ.ONE )
     $               GO TO 90
               END IF
               TJJ = CABS1( TJJS )
               IF( TJJ.GT.SMLNUM ) THEN
*
*                    abs(A(j,j)) > SMLNUM:
*
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by 1/b(j).
*
                        REC = ONE / XJ
                        CALL PCSSCAL( N, REC, X, IX, JX, DESCX, 1 )
                        XJTMP = XJTMP*REC
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
*                 X( J ) = CLADIV( X( J ), TJJS )
*                 XJ = CABS1( X( J ) )
                  XJTMP = CLADIV( XJTMP, TJJS )
                  XJ = CABS1( XJTMP )
                  IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
     $                 THEN
                     X( IROWX ) = XJTMP
                  END IF
               ELSE IF( TJJ.GT.ZERO ) THEN
*
*                    0 < abs(A(j,j)) <= SMLNUM:
*
                  IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
*                       to avoid overflow when dividing by A(j,j).
*
                     REC = ( TJJ*BIGNUM ) / XJ
                     IF( CNORM( J ).GT.ONE ) THEN
*
*                          Scale by 1/CNORM(j) to avoid overflow when
*                          multiplying x(j) times column j.
*
                        REC = REC / CNORM( J )
                     END IF
                     CALL PCSSCAL( N, REC, X, IX, JX, DESCX, 1 )
                     XJTMP = XJTMP*REC
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
*                 X( J ) = CLADIV( X( J ), TJJS )
*                 XJ = CABS1( X( J ) )
                  XJTMP = CLADIV( XJTMP, TJJS )
                  XJ = CABS1( XJTMP )
                  IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
     $                 THEN
                     X( IROWX ) = XJTMP
                  END IF
               ELSE
*
*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                    scale = 0, and compute a solution to A*x = 0.
*
                  CALL PCLASET( ' ', N, 1, CZERO, CZERO, X, IX, JX,
     $                          DESCX )
                  IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
     $                 THEN
                     X( IROWX ) = CONE
                  END IF
                  XJTMP = CONE
                  XJ = ONE
                  SCALE = ZERO
                  XMAX = ZERO
               END IF
   90          CONTINUE
*
*              Scale x if necessary to avoid overflow when adding a
*              multiple of column j of A.
*
               IF( XJ.GT.ONE ) THEN
                  REC = ONE / XJ
                  IF( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) THEN
*
*                    Scale x by 1/(2*abs(x(j))).
*
                     REC = REC*HALF
                     CALL PCSSCAL( N, REC, X, IX, JX, DESCX, 1 )
                     XJTMP = XJTMP*REC
                     SCALE = SCALE*REC
                  END IF
               ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN
*
*                 Scale x by 1/2.
*
                  CALL PCSSCAL( N, HALF, X, IX, JX, DESCX, 1 )
                  XJTMP = XJTMP*HALF
                  SCALE = SCALE*HALF
               END IF
*
               IF( UPPER ) THEN
                  IF( J.GT.1 ) THEN
*
*                    Compute the update
*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
*
                     ZDUM = -XJTMP*TSCAL
                     CALL PCAXPY( J-1, ZDUM, A, IA, JA+J-1, DESCA, 1, X,
     $                            IX, JX, DESCX, 1 )
                     CALL PCAMAX( J-1, ZDUM, IMAX, X, IX, JX, DESCX, 1 )
                     XMAX = CABS1( ZDUM )
                     CALL SGSUM2D( CONTXT, 'Row', ' ', 1, 1, XMAX, 1,
     $                             -1, -1 )
                  END IF
               ELSE
                  IF( J.LT.N ) THEN
*
*                    Compute the update
*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
*
                     ZDUM = -XJTMP*TSCAL
                     CALL PCAXPY( N-J, ZDUM, A, IA+J, JA+J-1, DESCA, 1,
     $                            X, IX+J, JX, DESCX, 1 )
                     CALL PCAMAX( N-J, ZDUM, I, X, IX+J, JX, DESCX, 1 )
                     XMAX = CABS1( ZDUM )
                     CALL SGSUM2D( CONTXT, 'Row', ' ', 1, 1, XMAX, 1,
     $                             -1, -1 )
                  END IF
               END IF
  100       CONTINUE
*
         ELSE IF( LSAME( TRANS, 'T' ) ) THEN
*
*           Solve A**T * x = b
*
            DO 120 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
*              XJ = CABS1( X( J ) )
               CALL INFOG2L( IX+J-1, JX, DESCX, NPROW, NPCOL, MYROW,
     $                       MYCOL, IROWX, ICOLX, ITMP1X, ITMP2X )
               IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) ) THEN
                  XJTMP = X( IROWX )
                  CALL CGEBS2D( CONTXT, 'All', ' ', 1, 1, XJTMP, 1 )
               ELSE
                  CALL CGEBR2D( CONTXT, 'All', ' ', 1, 1, XJTMP, 1,
     $                          ITMP1X, ITMP2X )
               END IF
               XJ = CABS1( XJTMP )
               USCAL = CMPLX( TSCAL )
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  REC = REC*HALF
                  IF( NOUNIT ) THEN
*                    TJJS = A( J, J )*TSCAL
                     CALL INFOG2L( IA+J-1, JA+J-1, DESCA, NPROW, NPCOL,
     $                             MYROW, MYCOL, IROW, ICOL, ITMP1,
     $                             ITMP2 )
                     IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) )
     $                    THEN
                        TJJS = A( ( ICOL-1 )*LDA+IROW )*TSCAL
                        CALL CGEBS2D( CONTXT, 'All', ' ', 1, 1, TJJS,
     $                                1 )
                     ELSE
                        CALL CGEBR2D( CONTXT, 'All', ' ', 1, 1, TJJS, 1,
     $                                ITMP1, ITMP2 )
                     END IF
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.ONE ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = CLADIV( USCAL, TJJS )
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL PCSSCAL( N, REC, X, IX, JX, DESCX, 1 )
                     XJTMP = XJTMP*REC
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
*
               CSUMJ = CZERO
               IF( USCAL.EQ.CONE ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call PCDOTU to perform the dot product.
*
                  IF( UPPER ) THEN
                     CALL PCDOTU( J-1, CSUMJ, A, IA, JA+J-1, DESCA, 1,
     $                            X, IX, JX, DESCX, 1 )
                  ELSE IF( J.LT.N ) THEN
                     CALL PCDOTU( N-J, CSUMJ, A, IA+J, JA+J-1, DESCA, 1,
     $                            X, IX+J, JX, DESCX, 1 )
                  END IF
                  IF( MYCOL.EQ.ITMP2X ) THEN
                     CALL CGEBS2D( CONTXT, 'Row', ' ', 1, 1, CSUMJ, 1 )
                  ELSE
                     CALL CGEBR2D( CONTXT, 'Row', ' ', 1, 1, CSUMJ, 1,
     $                             MYROW, ITMP2X )
                  END IF
               ELSE
*
*                 Otherwise, scale column of A by USCAL before dot
*                 product.  Below is not the best way to do it.
*
                  IF( UPPER ) THEN
*                    DO 130 I = 1, J - 1
*                       CSUMJ = CSUMJ + ( A( I, J )*USCAL )*X( I )
* 130                CONTINUE
                     ZDUM = CONJG( USCAL )
                     CALL PCSCAL( J-1, ZDUM, A, IA, JA+J-1, DESCA, 1 )
                     CALL PCDOTU( J-1, CSUMJ, A, IA, JA+J-1, DESCA, 1,
     $                            X, IX, JX, DESCX, 1 )
                     ZDUM = CLADIV( ZDUM, USCAL )
                     CALL PCSCAL( J-1, ZDUM, A, IA, JA+J-1, DESCA, 1 )
                  ELSE IF( J.LT.N ) THEN
*                    DO 140 I = J + 1, N
*                       CSUMJ = CSUMJ + ( A( I, J )*USCAL )*X( I )
*  140               CONTINUE
                     ZDUM = CONJG( USCAL )
                     CALL PCSCAL( N-J, ZDUM, A, IA+J, JA+J-1, DESCA, 1 )
                     CALL PCDOTU( N-J, CSUMJ, A, IA+J, JA+J-1, DESCA, 1,
     $                            X, IX+J, JX, DESCX, 1 )
                     ZDUM = CLADIV( ZDUM, USCAL )
                     CALL PCSCAL( N-J, ZDUM, A, IA+J, JA+J-1, DESCA, 1 )
                  END IF
                  IF( MYCOL.EQ.ITMP2X ) THEN
                     CALL CGEBS2D( CONTXT, 'Row', ' ', 1, 1, CSUMJ, 1 )
                  ELSE
                     CALL CGEBR2D( CONTXT, 'Row', ' ', 1, 1, CSUMJ, 1,
     $                             MYROW, ITMP2X )
                  END IF
               END IF
*
               IF( USCAL.EQ.CMPLX( TSCAL ) ) THEN
*
*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
*                 X( J ) = X( J ) - CSUMJ
*                 XJ = CABS1( X( J ) )
                  XJTMP = XJTMP - CSUMJ
                  XJ = CABS1( XJTMP )
*                  IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
*     $               X( IROWX ) = XJTMP
                  IF( NOUNIT ) THEN
*                    TJJS = A( J, J )*TSCAL
                     CALL INFOG2L( IA+J-1, JA+J-1, DESCA, NPROW, NPCOL,
     $                             MYROW, MYCOL, IROW, ICOL, ITMP1,
     $                             ITMP2 )
                     IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) )
     $                    THEN
                        TJJS = A( ( ICOL-1 )*LDA+IROW )*TSCAL
                        CALL CGEBS2D( CONTXT, 'All', ' ', 1, 1, TJJS,
     $                                1 )
                     ELSE
                        CALL CGEBR2D( CONTXT, 'All', ' ', 1, 1, TJJS, 1,
     $                                ITMP1, ITMP2 )
                     END IF
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE )
     $                  GO TO 110
                  END IF
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           REC = ONE / XJ
                           CALL PCSSCAL( N, REC, X, IX, JX, DESCX, 1 )
                           XJTMP = XJTMP*REC
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
*                    X( J ) = CLADIV( X( J ), TJJS )
                     XJTMP = CLADIV( XJTMP, TJJS )
                     IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
     $                    THEN
                        X( IROWX ) = XJTMP
                     END IF
                  ELSE IF( TJJ.GT.ZERO ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL PCSSCAL( N, REC, X, IX, JX, DESCX, 1 )
                        XJTMP = XJTMP*REC
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
*                    X( J ) = CLADIV( X( J ), TJJS )
                     XJTMP = CLADIV( XJTMP, TJJS )
                     IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
     $                    THEN
                        X( IROWX ) = XJTMP
                     END IF
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0 and compute a solution to A**T *x = 0.
*
                     CALL PCLASET( ' ', N, 1, CZERO, CZERO, X, IX, JX,
     $                             DESCX )
                     IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
     $                    THEN
                        X( IROWX ) = CONE
                     END IF
                     XJTMP = CONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  110             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
*                 product has already been divided by 1/A(j,j).
*
*                 X( J ) = CLADIV( X( J ), TJJS ) - CSUMJ
                  XJTMP = CLADIV( XJTMP, TJJS ) - CSUMJ
                  IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
     $                 THEN
                     X( IROWX ) = XJTMP
                  END IF
               END IF
               XMAX = MAX( XMAX, CABS1( XJTMP ) )
  120       CONTINUE
*
         ELSE
*
*           Solve A**H * x = b
*
            DO 140 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               CALL INFOG2L( IX+J-1, JX, DESCX, NPROW, NPCOL, MYROW,
     $                       MYCOL, IROWX, ICOLX, ITMP1X, ITMP2X )
               IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) ) THEN
                  XJTMP = X( IROWX )
                  CALL CGEBS2D( CONTXT, 'All', ' ', 1, 1, XJTMP, 1 )
               ELSE
                  CALL CGEBR2D( CONTXT, 'All', ' ', 1, 1, XJTMP, 1,
     $                          ITMP1X, ITMP2X )
               END IF
               XJ = CABS1( XJTMP )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  REC = REC*HALF
                  IF( NOUNIT ) THEN
*                    TJJS = CONJG( A( J, J ) )*TSCAL
                     CALL INFOG2L( IA+J-1, JA+J-1, DESCA, NPROW, NPCOL,
     $                             MYROW, MYCOL, IROW, ICOL, ITMP1,
     $                             ITMP2 )
                     IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) )
     $                    THEN
                        TJJS = CONJG( A( ( ICOL-1 )*LDA+IROW ) )*TSCAL
                        CALL CGEBS2D( CONTXT, 'All', ' ', 1, 1, TJJS,
     $                                1 )
                     ELSE
                        CALL CGEBR2D( CONTXT, 'All', ' ', 1, 1, TJJS, 1,
     $                                ITMP1, ITMP2 )
                     END IF
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.ONE ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = CLADIV( USCAL, TJJS )
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL PCSSCAL( N, REC, X, IX, JX, DESCX, 1 )
                     XJTMP = XJTMP*REC
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
*
               CSUMJ = CZERO
               IF( USCAL.EQ.CONE ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call PCDOTC to perform the dot product.
*
                  IF( UPPER ) THEN
                     CALL PCDOTC( J-1, CSUMJ, A, IA, JA+J-1, DESCA, 1,
     $                            X, IX, JX, DESCX, 1 )
                  ELSE IF( J.LT.N ) THEN
                     CALL PCDOTC( N-J, CSUMJ, A, IA+J, JA+J-1, DESCA, 1,
     $                            X, IX+J, JX, DESCX, 1 )
                  END IF
                  IF( MYCOL.EQ.ITMP2X ) THEN
                     CALL CGEBS2D( CONTXT, 'Row', ' ', 1, 1, CSUMJ, 1 )
                  ELSE
                     CALL CGEBR2D( CONTXT, 'Row', ' ', 1, 1, CSUMJ, 1,
     $                             MYROW, ITMP2X )
                  END IF
               ELSE
*
*                 Otherwise, scale column of A by USCAL before dot
*                 product.  Below is not the best way to do it.
*
                  IF( UPPER ) THEN
*                    DO 180 I = 1, J - 1
*                       CSUMJ = CSUMJ + ( CONJG( A( I, J ) )*USCAL )*
*    $                          X( I )
* 180                CONTINUE
                     ZDUM = CONJG( USCAL )
                     CALL PCSCAL( J-1, ZDUM, A, IA, JA+J-1, DESCA, 1 )
                     CALL PCDOTC( J-1, CSUMJ, A, IA, JA+J-1, DESCA, 1,
     $                            X, IX, JX, DESCX, 1 )
                     ZDUM = CLADIV( CONE, ZDUM )
                     CALL PCSCAL( J-1, ZDUM, A, IA, JA+J-1, DESCA, 1 )
                  ELSE IF( J.LT.N ) THEN
*                    DO 190 I = J + 1, N
*                       CSUMJ = CSUMJ + ( CONJG( A( I, J ) )*USCAL )*
*    $                          X( I )
* 190                CONTINUE
                     ZDUM = CONJG( USCAL )
                     CALL PCSCAL( N-J, ZDUM, A, IA+J, JA+J-1, DESCA, 1 )
                     CALL PCDOTC( N-J, CSUMJ, A, IA+J, JA+J-1, DESCA, 1,
     $                            X, IX+J, JX, DESCX, 1 )
                     ZDUM = CLADIV( CONE, ZDUM )
                     CALL PCSCAL( N-J, ZDUM, A, IA+J, JA+J-1, DESCA, 1 )
                  END IF
                  IF( MYCOL.EQ.ITMP2X ) THEN
                     CALL CGEBS2D( CONTXT, 'Row', ' ', 1, 1, CSUMJ, 1 )
                  ELSE
                     CALL CGEBR2D( CONTXT, 'Row', ' ', 1, 1, CSUMJ, 1,
     $                             MYROW, ITMP2X )
                  END IF
               END IF
*
               IF( USCAL.EQ.CMPLX( TSCAL ) ) THEN
*
*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
*                 X( J ) = X( J ) - CSUMJ
*                 XJ = CABS1( X( J ) )
                  XJTMP = XJTMP - CSUMJ
                  XJ = CABS1( XJTMP )
*                  IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
*     $               X( IROWX ) = XJTMP
                  IF( NOUNIT ) THEN
*                    TJJS = CONJG( A( J, J ) )*TSCAL
                     CALL INFOG2L( IA+J-1, JA+J-1, DESCA, NPROW, NPCOL,
     $                             MYROW, MYCOL, IROW, ICOL, ITMP1,
     $                             ITMP2 )
                     IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) )
     $                    THEN
                        TJJS = CONJG( A( ( ICOL-1 )*LDA+IROW ) )*TSCAL
                        CALL CGEBS2D( CONTXT, 'All', ' ', 1, 1, TJJS,
     $                                1 )
                     ELSE
                        CALL CGEBR2D( CONTXT, 'All', ' ', 1, 1, TJJS, 1,
     $                                ITMP1, ITMP2 )
                     END IF
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE )
     $                  GO TO 130
                  END IF
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                  TJJ = CABS1( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           REC = ONE / XJ
                           CALL PCSSCAL( N, REC, X, IX, JX, DESCX, 1 )
                           XJTMP = XJTMP*REC
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
*                    X( J ) = CLADIV( X( J ), TJJS )
                     XJTMP = CLADIV( XJTMP, TJJS )
                     IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
     $                  X( IROWX ) = XJTMP
                  ELSE IF( TJJ.GT.ZERO ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL PCSSCAL( N, REC, X, IX, JX, DESCX, 1 )
                        XJTMP = XJTMP*REC
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
*                    X( J ) = CLADIV( X( J ), TJJS )
                     XJTMP = CLADIV( XJTMP, TJJS )
                     IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
     $                  X( IROWX ) = XJTMP
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0 and compute a solution to A**H *x = 0.
*
                     CALL PCLASET( ' ', N, 1, CZERO, CZERO, X, IX, JX,
     $                             DESCX )
                     IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
     $                  X( IROWX ) = CONE
                     XJTMP = CONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  130             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
*                 product has already been divided by 1/A(j,j).
*
*                 X( J ) = CLADIV( X( J ), TJJS ) - CSUMJ
                  XJTMP = CLADIV( XJTMP, TJJS ) - CSUMJ
                  IF( ( MYROW.EQ.ITMP1X ) .AND. ( MYCOL.EQ.ITMP2X ) )
     $               X( IROWX ) = XJTMP
               END IF
               XMAX = MAX( XMAX, CABS1( XJTMP ) )
  140       CONTINUE
         END IF
         SCALE = SCALE / TSCAL
      END IF
*
*     Scale the column norms by 1/TSCAL for return.
*
      IF( TSCAL.NE.ONE ) THEN
         CALL SSCAL( N, ONE / TSCAL, CNORM, 1 )
      END IF
*
      RETURN
*
*     End of PCLATTRS
*
      END
