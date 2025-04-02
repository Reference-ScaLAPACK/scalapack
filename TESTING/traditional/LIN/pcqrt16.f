      SUBROUTINE PCQRT16( TRANS, M, N, NRHS, A, IA, JA, DESCA, X, IX,
     $                    JX, DESCX, B, IB, JB, DESCB, RWORK, RESID )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IA, IB, IX, JA, JB, JX, M, N, NRHS
      REAL               RESID
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCX( * )
      REAL               RWORK( * )
      COMPLEX            A( * ), B( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  PCQRT16 computes the residual for a solution of a system of linear
*  equations  sub( A )*sub( X ) = B  or  sub( A' )*sub( X ) = B:
*     RESID = norm(B - sub( A )*sub( X ) ) /
*             ( max(m,n) * norm(sub( A ) ) * norm(sub( X ) ) * EPS ),
*  where EPS is the machine epsilon, sub( A ) denotes
*  A(IA:IA+N-1,JA,JA+N-1), and sub( X ) denotes
*  X(IX:IX+N-1, JX:JX+NRHS-1).
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
*          Specifies the form of the system of equations:
*          = 'N':  sub( A )*sub( X ) = sub( B )
*          = 'T':  sub( A' )*sub( X )= sub( B ), where A' is the
*                  transpose of sub( A ).
*          = 'C':  sub( A' )*sub( X )= B, where A' is the conjugate
*                  transpose of sub( A ).
*
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the distributed submatrix sub( B ). NRHS >= 0.
*
*  A       (local input) COMPLEX pointer into the local
*          memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          The original M x N matrix A.
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
*  X       (local input) COMPLEX pointer into the local
*          memory to an array of dimension (LLD_X,LOCc(JX+NRHS-1)). This
*          array contains the local pieces of the computed solution
*          distributed vectors for the system of linear equations.
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
*  B       (local input/local output) COMPLEX pointer into
*          the local memory to an array of dimension
*          (LLD_B,LOCc(JB+NRHS-1)).  On entry, this array contains the
*          local pieces of the distributes right hand side vectors for
*          the system of linear equations. On exit, sub( B ) is over-
*          written with the difference sub( B ) - sub( A )*sub( X ) or
*          sub( B ) - sub( A )'*sub( X ).
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
*  RWORK   (local workspace) REAL array, dimension (LRWORK)
*          LWORK >= Nq0 if TRANS = 'N', and LRWORK >= Mp0 otherwise.
*
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
*  RESID   (global output) REAL
*          The maximum over the number of right hand sides of
*          norm( sub( B )- sub( A )*sub( X ) ) /
*          ( max(m,n) * norm( sub( A ) ) * norm( sub( X ) ) * EPS ).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            ICTXT, IDUMM, J, MYCOL, MYROW, N1, N2, NPCOL,
     $                   NPROW
      REAL               ANORM, BNORM, EPS, XNORM
*     ..
*     .. Local Arrays ..
      REAL               TEMP( 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               PCLANGE, PSLAMCH
      EXTERNAL           LSAME, PCLANGE, PSLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PCGEMM, PSCASUM,
     $                   SGAMX2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Quick exit if M = 0 or N = 0 or NRHS = 0
*
      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.EQ.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
      IF( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) ) THEN
         ANORM = PCLANGE( 'I', M, N, A, IA, JA, DESCA, RWORK )
         N1 = N
         N2 = M
      ELSE
         ANORM = PCLANGE( '1', M, N, A, IA, JA, DESCA, RWORK )
         N1 = M
         N2 = N
      END IF
*
      EPS = PSLAMCH( ICTXT, 'Epsilon' )
*
*     Compute  B - sub( A )*sub( X )  (or  B - sub( A' )*sub( X ) ) and
*     store in B.
*
      CALL PCGEMM( TRANS, 'No transpose', N1, NRHS, N2, -CONE, A, IA,
     $            JA, DESCA, X, IX, JX, DESCX, CONE, B, IB, JB, DESCB )
*
*     Compute the maximum over the number of right hand sides of
*        norm( sub( B ) - sub( A )*sub( X ) ) /
*        ( max(m,n) * norm( sub( A ) ) * norm( sub( X ) ) * EPS ).
*
      RESID = ZERO
      DO 10 J = 1, NRHS
*
         CALL PSCASUM( N1, BNORM, B, IB, JB+J-1, DESCB, 1 )
         CALL PSCASUM( N2, XNORM, X, IX, JX+J-1, DESCX, 1 )
*
*        Only the process columns owning the vector operands will have
*        the correct result, the other will have zero.
*
         TEMP( 1 ) = BNORM
         TEMP( 2 ) = XNORM
         IDUMM = 0
         CALL SGAMX2D( ICTXT, 'All', ' ', 2, 1, TEMP, 2, IDUMM, IDUMM,
     $                 -1, -1, IDUMM )
         BNORM = TEMP( 1 )
         XNORM = TEMP( 2 )
*
*        Every processes have ANORM, BNORM and XNORM now.
*
         IF( ANORM.EQ.ZERO .AND. BNORM.EQ.ZERO ) THEN
            RESID = ZERO
         ELSE IF( ANORM.LE.ZERO .OR. XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) /
     $              ( MAX( M, N )*EPS ) )
         END IF
*
   10 CONTINUE
*
      RETURN
*
*     End of PCQRT16
*
      END
