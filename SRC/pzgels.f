      SUBROUTINE PZGELS( TRANS, M, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     $                   DESCB, WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX*16         A( * ), B( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZGELS solves overdetermined or underdetermined complex linear
*  systems involving an M-by-N matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1),
*  or its conjugate-transpose, using a QR or LQ factorization of
*  sub( A ).  It is assumed that sub( A ) has full rank.
*
*  The following options are provided:
*
*  1. If TRANS = 'N' and m >= n:  find the least squares solution of
*     an overdetermined system, i.e., solve the least squares problem
*                  minimize || sub( B ) - sub( A )*X ||.
*
*  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
*     an underdetermined system sub( A ) * X = sub( B ).
*
*  3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
*     an undetermined system sub( A )**H * X = sub( B ).
*
*  4. If TRANS = 'C' and m < n:  find the least squares solution of
*     an overdetermined system, i.e., solve the least squares problem
*                  minimize || sub( B ) - sub( A )**H * X ||.
*
*  where sub( B ) denotes B( IB:IB+M-1, JB:JB+NRHS-1 ) when TRANS = 'N'
*  and B( IB:IB+N-1, JB:JB+NRHS-1 ) otherwise. Several right hand side
*  vectors b and solution vectors x can be handled in a single call;
*  When TRANS = 'N', the solution vectors are stored as the columns of
*  the N-by-NRHS right hand side matrix sub( B ) and the M-by-NRHS
*  right hand side matrix sub( B ) otherwise.
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
*  TRANS   (global input) CHARACTER
*          = 'N': the linear system involves sub( A );
*          = 'C': the linear system involves sub( A )**H.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of
*          rows of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e. the number of columns
*          of the distributed submatrices sub( B ) and X.  NRHS >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of local dimension
*          ( LLD_A, LOCc(JA+N-1) ).  On entry, the M-by-N matrix A.
*          if M >= N, sub( A ) is overwritten by details of its QR
*            factorization as returned by PZGEQRF;
*          if M <  N, sub( A ) is overwritten by details of its LQ
*            factorization as returned by PZGELQF.
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
*  B       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of local dimension
*          (LLD_B, LOCc(JB+NRHS-1)).  On entry, this array contains the
*          local pieces of the distributed matrix B of right hand side
*          vectors, stored columnwise;
*          sub( B ) is M-by-NRHS if TRANS='N', and N-by-NRHS otherwise.
*          On exit, sub( B ) is overwritten by the solution vectors,
*          stored columnwise:  if TRANS = 'N' and M >= N, rows 1 to N
*          of sub( B ) contain the least squares solution vectors; the
*          residual sum of squares for the solution in each column is
*          given by the sum of squares of elements N+1 to M in that
*          column; if TRANS = 'N' and M < N, rows 1 to N of sub( B )
*          contain the minimum norm solution vectors; if TRANS = 'C'
*          and M >= N, rows 1 to M of sub( B ) contain the minimum norm
*          solution vectors; if TRANS = 'C' and M < N, rows 1 to M of
*          sub( B ) contain the least squares solution vectors; the
*          residual sum of squares for the solution in each column is
*          given by the sum of squares of elements M+1 to N in that
*          column.
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
*  WORK    (local workspace/local output) COMPLEX*16 array,
*                                                  dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= LTAU + MAX( LWF, LWS ) where
*          If M >= N, then
*            LTAU = NUMROC( JA+MIN(M,N)-1, NB_A, MYCOL, CSRC_A, NPCOL ),
*            LWF  = NB_A * ( MpA0 + NqA0 + NB_A )
*            LWS  = MAX( (NB_A*(NB_A-1))/2, (NRHSqB0 + MpB0)*NB_A ) +
*                   NB_A * NB_A
*          Else
*            LTAU = NUMROC( IA+MIN(M,N)-1, MB_A, MYROW, RSRC_A, NPROW ),
*            LWF  = MB_A * ( MpA0 + NqA0 + MB_A )
*            LWS  = MAX( (MB_A*(MB_A-1))/2, ( NpB0 + MAX( NqA0 +
*                   NUMROC( NUMROC( N+IROFFB, MB_A, 0, 0, NPROW ),
*                   MB_A, 0, 0, LCMP ), NRHSqB0 ) )*MB_A ) +
*                   MB_A * MB_A
*          End if
*
*          where LCMP = LCM / NPROW with LCM = ILCM( NPROW, NPCOL ),
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          MpA0 = NUMROC( M+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*          NqA0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),
*
*          IROFFB = MOD( IB-1, MB_B ), ICOFFB = MOD( JB-1, NB_B ),
*          IBROW = INDXG2P( IB, MB_B, MYROW, RSRC_B, NPROW ),
*          IBCOL = INDXG2P( JB, NB_B, MYCOL, CSRC_B, NPCOL ),
*          MpB0 = NUMROC( M+IROFFB, MB_B, MYROW, IBROW, NPROW ),
*          NpB0 = NUMROC( N+IROFFB, MB_B, MYROW, IBROW, NPROW ),
*          NRHSqB0 = NUMROC( NRHS+ICOFFB, NB_B, MYCOL, IBCOL, NPCOL ),
*
*          ILCM, INDXG2P and NUMROC are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
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
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, TPSD
      INTEGER            BROW, IACOL, IAROW, IASCL, IBCOL, IBROW, IBSCL,
     $                   ICOFFA, ICOFFB, ICTXT, IPW, IROFFA, IROFFB,
     $                   LCM, LCMP, LTAU, LWF, LWMIN, LWS, MPA0, MPB0,
     $                   MYCOL, MYROW, NPB0, NPCOL, NPROW, NQA0,
     $                   NRHSQB0, SCLLEN
      DOUBLE PRECISION   ANRM, BIGNUM, BNRM, SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 2 ), IDUM2( 2 )
      DOUBLE PRECISION   RWORK( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILCM
      INTEGER            INDXG2P, NUMROC
      DOUBLE PRECISION   PDLAMCH, PZLANGE
      EXTERNAL           ILCM, INDXG2P, LSAME, NUMROC, PDLAMCH,
     $                   PZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK2MAT, PZGELQF,
     $                   PZGEQRF, PDLABAD, PZLASCL, PZLASET,
     $                   PZTRSM, PZUNMLQ, PZUNMQR, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, ICHAR, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 800 + CTXT_ )
      ELSE
         CALL CHK1MAT( M, 2, N, 3, IA, JA, DESCA, 8, INFO )
         IF ( M .GE. N ) THEN
            CALL CHK1MAT( M, 2, NRHS, 4, IB, JB, DESCB, 12, INFO )
         ELSE
            CALL CHK1MAT( N, 3, NRHS, 4, IB, JB, DESCB, 12, INFO )
         ENDIF
         IF( INFO.EQ.0 ) THEN
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            ICOFFB = MOD( JB-1, DESCB( NB_ ) )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IACOL = INDXG2P( IA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            MPA0 = NUMROC( M+IROFFA, DESCA( MB_ ), MYROW, IAROW, NPROW )
            NQA0 = NUMROC( N+ICOFFA, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
*
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $                       NPROW )
            IBCOL = INDXG2P( IB, DESCB( NB_ ), MYCOL, DESCB( CSRC_ ),
     $                       NPCOL )
            NRHSQB0 = NUMROC( NRHS+ICOFFB, DESCB( NB_ ), MYCOL, IBCOL,
     $                        NPCOL )
            IF( M.GE.N ) THEN
               MPB0 = NUMROC( M+IROFFB, DESCB( MB_ ), MYROW, IBROW,
     $                        NPROW )
               LTAU = NUMROC( JA+MIN(M,N)-1, DESCA( NB_ ), MYCOL,
     $                        DESCA( CSRC_ ), NPCOL )
               LWF  = DESCA( NB_ ) * ( MPA0 + NQA0 + DESCA( NB_ ) )
               LWS = MAX( ( DESCA( NB_ )*( DESCA( NB_ ) - 1 ) ) / 2,
     $               ( MPB0 + NRHSQB0 ) * DESCA( NB_ ) ) +
     $               DESCA( NB_ )*DESCA( NB_ )
            ELSE
               LCM = ILCM( NPROW, NPCOL )
               LCMP = LCM / NPROW
               NPB0 = NUMROC( N+IROFFB, DESCB( MB_ ), MYROW, IBROW,
     $                        NPROW )
               LTAU = NUMROC( IA+MIN(M,N)-1, DESCA( MB_ ), MYROW,
     $                        DESCA( RSRC_ ), NPROW )
               LWF  = DESCA( MB_ ) * ( MPA0 + NQA0 + DESCA( MB_ ) )
               LWS  = MAX( ( DESCA( MB_ )*( DESCA( MB_ ) - 1 ) ) / 2,
     $                ( NPB0 + MAX( NQA0 + NUMROC( NUMROC( N+IROFFB,
     $                DESCA( MB_ ), 0, 0, NPROW ), DESCA( MB_ ), 0, 0,
     $                LCMP ), NRHSQB0 ) )*DESCA( MB_ ) ) +
     $                DESCA( MB_ ) * DESCA( MB_ )
            END IF
            LWMIN = LTAU + MAX( LWF, LWS )
            WORK( 1 ) = DCMPLX( DBLE( LWMIN ) )
            LQUERY = ( LWORK.EQ.-1 )
*
            TPSD = .TRUE.
            IF( LSAME( TRANS, 'N' ) )
     $         TPSD = .FALSE.
*
            IF( .NOT.( LSAME( TRANS, 'N' ) .OR.
     $          LSAME( TRANS, 'C' ) ) ) THEN
               INFO = -1
            ELSE IF( M.LT.0 ) THEN
               INFO = -2
            ELSE IF( N.LT.0 ) THEN
               INFO = -3
            ELSE IF( NRHS.LT.0 ) THEN
               INFO = -4
            ELSE IF( M.GE.N .AND. IROFFA.NE.IROFFB ) THEN
               INFO = -10
            ELSE IF( M.GE.N .AND. IAROW.NE.IBROW ) THEN
               INFO = -10
            ELSE IF( M.LT.N .AND. ICOFFA.NE.IROFFB ) THEN
               INFO = -10
            ELSE IF( M.GE.N .AND. DESCA( MB_ ).NE.DESCB( MB_ ) ) THEN
               INFO = -( 1200 + MB_ )
            ELSE IF( M.LT.N .AND. DESCA( NB_ ).NE.DESCB( MB_ ) ) THEN
               INFO = -( 1200 + MB_ )
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -( 1200 + CTXT_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -14
            END IF
         END IF
*
         IF( .NOT.TPSD ) THEN
            IDUM1( 1 ) = ICHAR( 'N' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'C' )
         END IF
         IDUM2( 1 ) = 1
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 14
         CALL PCHK2MAT( M, 2, N, 3, IA, JA, DESCA, 8, N, 3, NRHS, 4,
     $                  IB, JB, DESCB, 12, 2, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZGELS', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         CALL PZLASET( 'Full', MAX( M, N ), NRHS, CZERO, CZERO, B,
     $                 IB, JB, DESCB )
         RETURN
      END IF
*
*     Get machine parameters
*
      SMLNUM = PDLAMCH( ICTXT, 'S' )
      SMLNUM = SMLNUM / PDLAMCH( ICTXT, 'P' )
      BIGNUM = ONE / SMLNUM
      CALL PDLABAD( ICTXT, SMLNUM, BIGNUM )
*
*     Scale A, B if max entry outside range [SMLNUM,BIGNUM]
*
      ANRM = PZLANGE( 'M', M, N, A, IA, JA, DESCA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL PZLASCL( 'G', ANRM, SMLNUM, M, N, A, IA, JA, DESCA,
     $                 INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL PZLASCL( 'G', ANRM, BIGNUM, M, N, A, IA, JA, DESCA,
     $                 INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*
*        Matrix all zero. Return zero solution.
*
         CALL PZLASET( 'F', MAX( M, N ), NRHS, CZERO, CZERO, B, IB,
     $                 JB, DESCB )
         GO TO 10
      END IF
*
      BROW = M
      IF( TPSD )
     $   BROW = N
*
      BNRM = PZLANGE( 'M', BROW, NRHS, B, IB, JB, DESCB, RWORK )
*
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL PZLASCL( 'G', BNRM, SMLNUM, BROW, NRHS, B, IB, JB,
     $                 DESCB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL PZLASCL( 'G', BNRM, BIGNUM, BROW, NRHS, B, IB, JB,
     $                 DESCB, INFO )
         IBSCL = 2
      END IF
*
      IPW = LTAU + 1
*
      IF( M.GE.N ) THEN
*
*        compute QR factorization of A
*
         CALL PZGEQRF( M, N, A, IA, JA, DESCA, WORK, WORK( IPW ),
     $                 LWORK-LTAU, INFO )
*
*        workspace at least N, optimally N*NB
*
         IF( .NOT.TPSD ) THEN
*
*           Least-Squares Problem min || A * X - B ||
*
*           B(IB:IB+M-1,JB:JB+NRHS-1) := Q' * B(IB:IB+M-1,JB:JB+NRHS-1)
*
            CALL PZUNMQR( 'Left', 'Conjugate transpose', M, NRHS, N, A,
     $                    IA, JA, DESCA, WORK, B, IB, JB, DESCB,
     $                    WORK( IPW ), LWORK-LTAU, INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
*           B(IB:IB+N-1,JB:JB+NRHS-1) := inv(R) *
*                                        B(IB:IB+N-1,JB:JB+NRHS-1)
*
            CALL PZTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $                   NRHS, CONE, A, IA, JA, DESCA, B, IB, JB,
     $                   DESCB )
*
            SCLLEN = N
*
         ELSE
*
*           Overdetermined system of equations sub( A )' * X = sub( B )
*
*           sub( B ) := inv(R') * sub( B )
*
            CALL PZTRSM( 'Left', 'Upper', 'Conjugate transpose',
     $                   'Non-unit', N, NRHS, CONE, A, IA, JA, DESCA,
     $                   B, IB, JB, DESCB )
*
*           B(IB+N:IB+M-1,JB:JB+NRHS-1) = ZERO
*
            CALL PZLASET( 'All', M-N, NRHS, CZERO, CZERO, B, IB+N, JB,
     $                    DESCB )
*
*           B(IB:IB+M-1,JB:JB+NRHS-1) := Q(1:N,:) *
*                                        B(IB:IB+N-1,JB:JB+NRHS-1)
*
            CALL PZUNMQR( 'Left', 'No transpose', M, NRHS, N, A, IA, JA,
     $                    DESCA, WORK, B, IB, JB, DESCB, WORK( IPW ),
     $                    LWORK-LTAU, INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
            SCLLEN = M
*
         END IF
*
      ELSE
*
*        Compute LQ factorization of sub( A )
*
         CALL PZGELQF( M, N, A, IA, JA, DESCA, WORK, WORK( IPW ),
     $                 LWORK-LTAU, INFO )
*
*        workspace at least M, optimally M*NB.
*
         IF( .NOT.TPSD ) THEN
*
*           underdetermined system of equations sub( A ) * X = sub( B )
*
*           B(IB:IB+M-1,JB:JB+NRHS-1) := inv(L) *
*                                        B(IB:IB+M-1,JB:JB+NRHS-1)
*
            CALL PZTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', M,
     $                   NRHS, CONE, A, IA, JA, DESCA, B, IB, JB,
     $                   DESCB )
*
*           B(IB+M:IB+N-1,JB:JB+NRHS-1) = 0
*
            CALL PZLASET( 'All', N-M, NRHS, CZERO, CZERO, B, IB+M, JB,
     $                    DESCB )
*
*           B(IB:IB+N-1,JB:JB+NRHS-1) := Q(1:N,:)' *
*                                        B(IB:IB+M-1,JB:JB+NRHS-1)
*
            CALL PZUNMLQ( 'Left', 'Conjugate transpose', N, NRHS, M, A,
     $                    IA, JA, DESCA, WORK, B, IB, JB, DESCB,
     $                    WORK( IPW ), LWORK-LTAU, INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
            SCLLEN = N
*
         ELSE
*
*           overdetermined system min || A' * X - B ||
*
*           B(IB:IB+N-1,JB:JB+NRHS-1) := Q * B(IB:IB+N-1,JB:JB+NRHS-1)
*
            CALL PZUNMLQ( 'Left', 'No transpose', N, NRHS, M, A, IA, JA,
     $                    DESCA, WORK, B, IB, JB, DESCB, WORK( IPW ),
     $                    LWORK-LTAU, INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
*           B(IB:IB+M-1,JB:JB+NRHS-1) := inv(L') *
*                                        B(IB:IB+M-1,JB:JB+NRHS-1)
*
            CALL PZTRSM( 'Left', 'Lower', 'Conjugate transpose',
     $                   'Non-unit', M, NRHS, CONE, A, IA, JA, DESCA,
     $                   B, IB, JB, DESCB )
*
            SCLLEN = M
*
         END IF
*
      END IF
*
*     Undo scaling
*
      IF( IASCL.EQ.1 ) THEN
         CALL PZLASCL( 'G', ANRM, SMLNUM, SCLLEN, NRHS, B, IB, JB,
     $                 DESCB, INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL PZLASCL( 'G', ANRM, BIGNUM, SCLLEN, NRHS, B, IB, JB,
     $                 DESCB, INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL PZLASCL( 'G', SMLNUM, BNRM, SCLLEN, NRHS, B, IB, JB,
     $                 DESCB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL PZLASCL( 'G', BIGNUM, BNRM, SCLLEN, NRHS, B, IB, JB,
     $                 DESCB, INFO )
      END IF
*
   10 CONTINUE
*
      WORK( 1 ) = DCMPLX( DBLE( LWMIN ) )
*
      RETURN
*
*     End of PZGELS
*
      END
