      DOUBLE PRECISION FUNCTION PDQRT17( TRANS, IRESID, M, N, NRHS, A,
     $                                   IA, JA, DESCA, X, IX, JX,
     $                                   DESCX, B, IB, JB, DESCB, WORK,
     $                                   RWORK )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IA, IB, IRESID, IX, JA, JB, JX, M, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCX( * )
      DOUBLE PRECISION   A( * ), B( * ), WORK( * ), X( * )
      DOUBLE PRECISION   RWORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDQRT17 computes the ratio
*
*     || R'*op( sub( A ) ) ||/(||sub( A )||*alpha*max(M,N,NRHS)*eps)
*
*  where R = op( sub( A ) )*sub( X ) - B, op(A) is A or A', and
*
*     alpha = ||B|| if IRESID = 1 (zero-residual problem)
*     alpha = ||R|| if IRESID = 2 (otherwise).
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
*  TRANS   (global input) CHARACTER*1
*          Specifies whether or not the transpose of sub( A ) is used.
*          = 'N':  No transpose, op( sub( A ) ) = sub( A ).
*          = 'T':  Transpose, op( sub( A ) ) = sub( A' ).
*
*  IRESID  (global input) INTEGER
*          IRESID = 1 indicates zero-residual problem.
*          IRESID = 2 indicates non-zero residual.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*          If TRANS = 'N', the number of rows of the distributed
*          submatrix sub( B ).
*          If TRANS = 'T', the number of rows of the distributed
*          submatrix sub( X ).
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*          If TRANS = 'N', the number of rows of the distributed
*          submatrix sub( X ). Otherwise N is the number of rows of
*          the distributed submatrix sub( B ).
*
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the distributed submatrices sub( X ) and sub( B ).
*          NRHS >= 0.
*
*  A       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_A,LOCc(JA+N-1)). This array
*          contains the local pieces of the distributed M-by-N
*          submatrix sub( A ).
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
*  X       (local input) DOUBLE PRECISION pointer into the local
*          memory to an array of dimension (LLD_X,LOCc(JX+NRHS-1)).
*          If TRANS = 'N', this array contains the local pieces of the
*          N-by-NRHS distributed submatrix sub( X ). Otherwise, this
*          array contains the local pieces of the M-by-NRHS distributed
*          submatrix sub( X ).
*
*  IX      (global input) INTEGER
*          The row index in the global array X indicating the first
*          row of sub( X ).
*
*  JX      (global input) INTEGER
*          The column index in the global array X indicating the
*          first column of sub( X ).
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*
*  B       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_B,LOCc(JB+NRHS-1)).
*          If TRANS='N', this array contains the local pieces of the
*          distributed M-by-NRHS submatrix operand sub( B ). Otherwise,
*          this array contains the local pieces of the distributed
*          N-by-NRHS submatrix operand sub( B ).
*
*  IB      (global input) INTEGER
*          The row index in the global array B indicating the first
*          row of sub( B ).
*
*  JB      (global input) INTEGER
*          The column index in the global array B indicating the
*          first column of sub( B ).
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
*
*  WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK)
*          If TRANS = 'N', LWORK >= Mp0 * NRHSq0 + NRHSp0 * Nq0
*          otherwise       LWORK >= Np0 * NRHSq0 + NRHSp0 * Mq0
*
*  RWORK   (local workspace) DOUBLE PRECISION array, dimension (LRWORK)
*          LRWORK >= Nq0, if TRANS = 'N', and LRWORK >= Mp0 otherwise.
*
*          where
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          Mp0 = NUMROC( M+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*          Np0 = NUMROC( N+ICOFFA, NB_A, MYROW, IAROW, NPROW ),
*          Mq0 = NUMROC( M+IROFFA, NB_A, MYCOL, IACOL, NPCOL ),
*          Nq0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),
*
*          IROFFB = MOD( IB-1, MB_B ), ICOFFB = MOD( JB-1, NB_B ),
*          IBROW = INDXG2P( IB, MB_B, MYROW, RSRC_B, NPROW ),
*          IBCOL = INDXG2P( JB, NB_B, MYCOL, CSRC_B, NPCOL ),
*          NRHSp0 = NUMROC( NRHS+ICOFFB, NB_B, MYROW, IBROW, NPROW ),
*          NRHSq0 = NUMROC( NRHS+ICOFFB, NB_B, MYCOL, IBCOL, NPCOL ).
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
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IACOL, IBCOL, IBROW, ICOFFB, ICTXT, INFO,
     $                   IOFFA, IROFFB, ISCL, IW, IW2, JW, JW2, MYCOL,
     $                   NRHSQ, NRHSP, MYROW, NCOLS, NPCOL, NPROW,
     $                   NROWS, NROWSP
      DOUBLE PRECISION   ERR, NORMA, NORMB, NORMRS, NORMX, SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER            DESCW( DLEN_ ), DESCW2( DLEN_ )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      DOUBLE PRECISION   PDLAMCH, PDLANGE
      EXTERNAL           INDXG2P, LSAME, NUMROC, PDLAMCH, PDLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, PDGEMM, PDLACPY,
     $                   PDLASCL, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Executable Statements ..
*
      PDQRT17 = ZERO
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      INFO = 0
      IF( LSAME( TRANS, 'N' ) ) THEN
         NROWS = M
         NCOLS = N
         IOFFA = MOD( JA-1, DESCA( NB_ ) )
      ELSE IF( LSAME( TRANS, 'T' ) ) THEN
         NROWS = N
         NCOLS = M
         IOFFA = MOD( IA-1, DESCA( MB_ ) )
      ELSE
         CALL PXERBLA( ICTXT, 'PDQRT17', -1 )
         RETURN
      END IF
*
      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.LE.0 )
     $   RETURN
*
      IROFFB = MOD( IA-1, DESCA( MB_ ) )
      ICOFFB = MOD( JA-1, DESCA( NB_ ) )
      IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $                 NPROW )
      IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                 NPCOL )
      IBCOL = INDXG2P( JB, DESCB( NB_ ), MYCOL, DESCB( CSRC_ ),
     $                 NPCOL )
*
      NRHSQ = NUMROC( NRHS+ICOFFB, DESCB( NB_ ), MYCOL, IBCOL, NPCOL )
      NRHSP = NUMROC( NRHS+IROFFB, DESCB( NB_ ), MYROW, IBROW, NPROW )
      NROWSP = NUMROC( NROWS+IROFFB, DESCA( MB_ ), MYROW, IBROW, NPROW )
*
*     Assign array descriptor DESCW for workspace WORK, where DESCW
*     holds a copy of the distributed submatrix sub( B )
*
      CALL DESCSET( DESCW, NROWS+IROFFB, NRHS+ICOFFB, DESCB( MB_ ),
     $              DESCB( NB_ ), IBROW, IBCOL, ICTXT, MAX( 1,
     $              NROWSP ) )
*
*     Assign array descriptor DESCW2 for workspace WORK, where DESCW2
*     holds a copy of the distributed submatrix sub( X**T )
*
      CALL DESCSET( DESCW2, NRHS+ICOFFB, NCOLS+IOFFA, DESCX( NB_ ),
     $              DESCX( MB_ ), IBROW, IACOL, ICTXT, MAX( 1,
     $              NRHSP ) )
*
      NORMA = PDLANGE( 'One-norm', M, N, A, IA, JA, DESCA, RWORK )
      SMLNUM = PDLAMCH( ICTXT, 'Safe minimum' )
      SMLNUM = SMLNUM / PDLAMCH( ICTXT, 'Precision' )
      ISCL = 0
*
*     compute residual and scale it
*
      IW = 1 + IROFFB
      JW = 1 + ICOFFB
      CALL PDLACPY( 'All', NROWS, NRHS, B, IB, JB, DESCB, WORK, IW, JW,
     $              DESCW )
      CALL PDGEMM( TRANS, 'No transpose', NROWS, NRHS, NCOLS, -ONE, A,
     $             IA, JA, DESCA, X, IX, JX, DESCX, ONE, WORK, IW, JW,
     $             DESCW )
      NORMRS = PDLANGE( 'Max', NROWS, NRHS, WORK, IW, JW, DESCW,
     $                  RWORK )
      IF( NORMRS.GT.SMLNUM ) THEN
         ISCL = 1
         CALL PDLASCL( 'General', NORMRS, ONE, NROWS, NRHS, WORK,
     $                 IW, JW, DESCW, INFO )
      END IF
*
*     compute R'*sub( A )
*
      IW2 = 1 + ICOFFB
      JW2 = 1 + IOFFA
      CALL PDGEMM( 'Transpose', TRANS, NRHS, NCOLS, NROWS, ONE, WORK,
     $             IW, JW, DESCW, A, IA, JA, DESCA, ZERO,
     $             WORK( NROWSP*NRHSQ+1 ), IW2, JW2, DESCW2 )
*
*     compute and properly scale error
*
      ERR = PDLANGE( 'One-norm', NRHS, NCOLS, WORK( NROWSP*NRHSQ+1 ),
     $               IW2, JW2, DESCW2, RWORK )
      IF( NORMA.NE.ZERO )
     $   ERR = ERR / NORMA
*
      IF( ISCL.EQ.1 )
     $   ERR = ERR*NORMRS
*
      IF( IRESID.EQ.1 ) THEN
         NORMB = PDLANGE( 'One-norm', NROWS, NRHS, B, IB, JB, DESCB,
     $                    RWORK )
         IF( NORMB.NE.ZERO )
     $      ERR = ERR / NORMB
      ELSE
         NORMX = PDLANGE( 'One-norm', NCOLS, NRHS, X, IX, JX, DESCX,
     $                    RWORK )
         IF( NORMX.NE.ZERO )
     $      ERR = ERR / NORMX
      END IF
*
      PDQRT17 = ERR / ( PDLAMCH( ICTXT, 'Epsilon' ) *
     $                  DBLE( MAX( M, N, NRHS ) ) )
*
      RETURN
*
*     End of PDQRT17
*
      END
