      SUBROUTINE PDPOCON( UPLO, N, A, IA, JA, DESCA, ANORM, RCOND, WORK,
     $                    LWORK, IWORK, LIWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, LIWORK, LWORK, N
      DOUBLE PRECISION   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), IWORK( * )
      DOUBLE PRECISION   A( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDPOCON estimates the reciprocal of the condition number (in the
*  1-norm) of a real symmetric positive definite distributed matrix
*  using the Cholesky factorization A = U**T*U or A = L*L**T computed by
*  PDPOTRF.
*
*  An estimate is obtained for norm(inv(A(IA:IA+N-1,JA:JA+N-1))), and
*  the reciprocal of the condition number is computed as
*             RCOND = 1 / ( norm( A(IA:IA+N-1,JA:JA+N-1)      ) *
*                           norm( inv(A(IA:IA+N-1,JA:JA+N-1)) ) ).
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
*  UPLO    (global input) CHARACTER
*          Specifies whether the factor stored in
*          A(IA:IA+N-1,JA:JA+N-1) is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The order of the distributed matrix A(IA:IA+N-1,JA:JA+N-1).
*          N >= 0.
*
*  A       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension ( LLD_A, LOCc(JA+N-1) ). On entry,
*          this array contains the local pieces of the factors L or U
*          from the Cholesky factorization A(IA:IA+N-1,JA:JA+N-1) = U'*U
*          or L*L', as computed by PDPOTRF.
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
*  ANORM   (global input) DOUBLE PRECISION
*          The 1-norm (or infinity-norm) of the symmetric distributed
*          matrix A(IA:IA+N-1,JA:JA+N-1).
*
*  RCOND   (global output) DOUBLE PRECISION
*          The reciprocal of the condition number of the distributed
*          matrix A(IA:IA+N-1,JA:JA+N-1), computed as
*             RCOND = 1 / ( norm( A(IA:IA+N-1,JA:JA+N-1)      ) *
*                           norm( inv(A(IA:IA+N-1,JA:JA+N-1)) ) ).
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                   dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= 2*LOCr(N+MOD(IA-1,MB_A)) + 2*LOCc(N+MOD(JA-1,NB_A))+
*          MAX( 2, MAX(NB_A*CEIL(NPROW-1,NPCOL),LOCc(N+MOD(JA-1,NB_A)) +
*          NB_A*CEIL(NPCOL-1,NPROW)) ).
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  IWORK   (local workspace/local output) INTEGER array,
*                                                 dimension (LIWORK)
*          On exit, IWORK(1) returns the minimal and optimal LIWORK.
*
*  LIWORK  (local or global input) INTEGER
*          The dimension of the array IWORK.
*          LIWORK is local input and must be at least
*          LIWORK >= LOCr(N+MOD(IA-1,MB_A)).
*
*          If LIWORK = -1, then LIWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
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
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      CHARACTER          CBTOP, COLCTOP, NORMIN, ROWCTOP
      INTEGER            IACOL, IAROW, ICOFF, ICTXT, IIA, IPNL, IPNU,
     $                   IPV, IPW, IPX, IROFF, IV, IX, IXX, JJA, JV,
     $                   JX, KASE, LIWMIN, LWMIN, MYCOL, MYROW, NP,
     $                   NPCOL, NPROW, NPMOD, NQ, NQMOD
      DOUBLE PRECISION   AINVNM, SCALE, SL, SU, SMLNUM
      DOUBLE PRECISION   WMAX
*     ..
*     .. Local Arrays ..
      INTEGER            DESCV( DLEN_ ), DESCX( DLEN_ ), IDUM1( 3 ),
     $                   IDUM2( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, DGEBR2D,
     $                   DGEBS2D, INFOG2L, PCHK1MAT, PDAMAX,
     $                   PDLATRS, PDLACON, PDRSCL, PB_TOPGET,
     $                   PB_TOPSET, PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, INDXG2P, NUMROC
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           ICEIL, INDXG2P, LSAME, NUMROC, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, ICHAR, MAX, MOD
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
         INFO = -(600+CTXT_)
      ELSE
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
         IF( INFO.EQ.0 ) THEN
            UPPER = LSAME( UPLO, 'U' )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            NPMOD = NUMROC( N + MOD( IA-1, DESCA( MB_ ) ), DESCA( MB_ ),
     $                      MYROW, IAROW, NPROW )
            NQMOD = NUMROC( N + MOD( JA-1, DESCA( NB_ ) ), DESCA( NB_ ),
     $                      MYCOL, IACOL, NPCOL )
            LWMIN = 2*NPMOD + 2*NQMOD +
     $              MAX( 2, MAX( DESCA( NB_ )*
     $                   MAX( 1, ICEIL( NPROW-1, NPCOL ) ), NQMOD +
     $                   DESCA( NB_ )*
     $                   MAX( 1, ICEIL( NPCOL-1, NPROW ) ) ) )
            WORK( 1 ) = DBLE( LWMIN )
            LIWMIN = NPMOD
            IWORK( 1 ) = LIWMIN
            LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( ANORM.LT.ZERO ) THEN
               INFO = -7
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -10
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               IWORK( 1 ) = LIWMIN
               INFO = -12
            END IF
         END IF
*
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 1
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 10
         IF( LIWORK.EQ.-1 ) THEN
            IDUM1( 3 ) = -1
         ELSE
            IDUM1( 3 ) = 1
         END IF
         IDUM2( 3 ) = 12
         CALL PCHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, 3, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDPOCON', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         RCOND = ONE
         RETURN
      END IF
*
      CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
      CALL PB_TOPGET( ICTXT, 'Combine', 'Rowwise',    ROWCTOP )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', '1-tree' )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise',    '1-tree' )
*
      SMLNUM = PDLAMCH( ICTXT, 'Safe minimum' )
      IROFF = MOD( IA-1, DESCA( MB_ ) )
      ICOFF = MOD( JA-1, DESCA( NB_ ) )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      NP = NUMROC( N+IROFF, DESCA( MB_ ), MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFF, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IV = IROFF + 1
      IX = IV
      JV = ICOFF + 1
      JX = JV
*
      IPX  = 1
      IPV  = IPX + NP
      IPNL = IPV + NP
      IPNU = IPNL + NQ
      IPW  = IPNU + NQ
*
      CALL DESCSET( DESCV, N+IROFF, 1, DESCA( MB_ ), 1, IAROW, MYCOL,
     $              ICTXT, MAX( 1, NP ) )
      CALL DESCSET( DESCX, N+IROFF, 1, DESCA( MB_ ), 1, IAROW, MYCOL,
     $              ICTXT, MAX( 1, NP ) )
*
*     Estimate the 1-norm (or I-norm) of inv(A).
*
      AINVNM = ZERO
      KASE   = 0
      NORMIN = 'N'
*
   10 CONTINUE
      CALL PDLACON( N, WORK( IPV ), IV, JV, DESCV, WORK( IPX ), IX, JX,
     $              DESCX, IWORK, AINVNM, KASE )
      IF( KASE.NE.0 ) THEN
         IF( UPPER ) THEN
*
*           Multiply by inv(U').
*
            DESCX( CSRC_ ) = IACOL
            CALL PDLATRS( 'Upper', 'Transpose', 'Non-unit', NORMIN,
     $                    N, A, IA, JA, DESCA, WORK( IPX ), IX, JX,
     $                    DESCX, SL, WORK( IPNL ), WORK( IPW ) )
            DESCX( CSRC_ ) = MYCOL
            NORMIN = 'Y'
*
*           Multiply by inv(U).
*
            DESCX( CSRC_ ) = IACOL
            CALL PDLATRS( 'Upper', 'No transpose', 'Non-unit', NORMIN,
     $                    N, A, IA, JA, DESCA, WORK( IPX ), IX, JX,
     $                    DESCX, SU, WORK( IPNU ), WORK( IPW ) )
            DESCX( CSRC_ ) = MYCOL
         ELSE
*
*           Multiply by inv(L).
*
            DESCX( CSRC_ ) = IACOL
            CALL PDLATRS( 'Lower', 'No transpose', 'Non-unit', NORMIN,
     $                    N, A, IA, JA, DESCA, WORK( IPX ), IX, JX,
     $                    DESCX, SL, WORK( IPNL ), WORK( IPW ) )
            DESCX( CSRC_ ) = MYCOL
            NORMIN = 'Y'
*
*           Multiply by inv(L').
*
            DESCX( CSRC_ ) = IACOL
            CALL PDLATRS( 'Lower', 'Transpose', 'Non-unit', NORMIN,
     $                    N, A, IA, JA, DESCA, WORK( IPX ), IX, JX,
     $                    DESCX, SU, WORK( IPNU ), WORK( IPW ) )
            DESCX( CSRC_ ) = MYCOL
         END IF
*
*        Multiply by 1/SCALE if doing so will not cause overflow.
*
         SCALE = SL*SU
         IF( SCALE.NE.ONE ) THEN
            CALL PDAMAX( N, WMAX, IXX, WORK( IPX ), IX, JX, DESCX, 1 )
            IF( DESCX( M_ ).EQ.1 .AND. N.EQ.1 ) THEN
               CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', CBTOP )
               IF( MYROW.EQ.IAROW ) THEN
                  CALL DGEBS2D( ICTXT, 'Column', CBTOP, 1, 1, WMAX, 1 )
               ELSE
                  CALL DGEBR2D( ICTXT, 'Column', CBTOP, 1, 1, WMAX, 1,
     $                          IAROW, MYCOL )
               END IF
            END IF
            IF( SCALE.LT.ABS( WMAX )*SMLNUM .OR. SCALE.EQ.ZERO )
     $         GO TO 20
            CALL PDRSCL( N, SCALE, WORK( IPX ), IX, JX, DESCX, 1 )
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO )
     $   RCOND = ( ONE / AINVNM ) / ANORM
*
   20 CONTINUE
*
      CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise',    ROWCTOP )
*
      RETURN
*
*     End of PDPOCON
*
      END
