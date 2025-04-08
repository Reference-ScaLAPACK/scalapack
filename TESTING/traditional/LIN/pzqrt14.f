      DOUBLE PRECISION FUNCTION PZQRT14( TRANS, M, N, NRHS, A, IA, JA,
     $                                   DESCA, X, IX, JX, DESCX, WORK )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IA, IX, JA, JX, M, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCX( * )
      COMPLEX*16         A( * ), WORK( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  PZQRT14 checks whether sub( X ) is in the row space of sub( A ) or
*  sub( A )', where sub( A ) denotes A( IA:IA+M-1, JA:JA+N-1 ) and
*  sub( X ) denotes X( IX:IX+N-1, JX:JX+NRHS-1 ) if TRANS = 'N', and
*  X( IX:IX+N-1, JX:JX+NRHS-1 ) otherwise.  It does so by scaling both
*  sub( X ) and sub( A ) such that their norms are in the range
*  [sqrt(eps), 1/sqrt(eps)], then computing an LQ factorization of
*  [sub( A )',sub( X )]' (if TRANS = 'N') or a QR factorization of
*  [sub( A ),sub( X )] otherwise, and returning the norm of the trailing
*  triangle, scaled by MAX(M,N,NRHS)*eps.
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
*          = 'N':  No transpose, check for sub( X ) in the row space of
*                  sub( A ),
*          = 'C':  Conjugate transpose, check for sub( X ) in row space
*                  of sub( A )'.
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
*          of the distributed submatrix sub( X ). NRHS >= 0.
*
*  A       (local input) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_A, LOCc(JA+N-1)). This array
*          contains the local pieces of the M-by-N distributed matrix
*          sub( A ).
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
*  X       (local input) COMPLEX*16 pointer into the local
*          memory to an array of dimension (LLD_X,LOCc(JX+NRHS-1)).
*          On entry, this array contains the local pieces of the
*          N-by-NRHS distributed submatrix sub( X ) if TRANS = 'N',
*          and the M-by-NRHS distributed submatrix sub( X ) otherwise.
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
*  WORK    (local workspace) COMPLEX*16 array dimension (LWORK)
*          If TRANS='N', LWORK >= MNRHSP * NQ + LTAU + LWF and
*          LWORK >= MP * NNRHSQ + LTAU + LWF otherwise, where
*
*          IF TRANS='N', (LQ fact)
*            MNRHSP = NUMROC( M+NRHS+IROFFA, MB_A, MYROW, IAROW,
*                             NPROW )
*            LTAU   = NUMROC( IA+MIN( M+NRHS, N )-1, MB_A, MYROW,
*                             RSRC_A, NPROW )
*            LWF    = MB_A * ( MB_A + MNRHSP + NQ0 )
*          ELSE         (QR fact)
*            NNRHSQ = NUMROC( N+NRHS+ICOFFA, NB_A, MYCOL, IACOL,
*                             NPCOL )
*            LTAU   = NUMROC( JA+MIN( M, N+NRHS )-1, NB_A, MYCOL,
*                             CSRC_A, NPCOL )
*            LWF    = NB_A * ( NB_A + MP0 + NNRHSQ )
*          END IF
*
*          and,
*
*          IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
*          MP0 = NUMROC( M+IROFFA, MB_A, MYROW, IAROW, NPROW ),
*          NQ0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ).
*
*          INDXG2P and NUMROC are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
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
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            TPSD
      INTEGER            IACOL, IAROW, ICOFFA, ICTXT, IDUM, IIA, INFO,
     $                   IPTAU, IPW, IPWA, IROFFA, IWA, IWX, J, JJA,
     $                   JWA, JWX, LDW, LWORK, MPWA, MPW, MQW, MYCOL,
     $                   MYROW, NPCOL, NPROW, NPW, NQWA, NQW
      DOUBLE PRECISION   ANRM, ERR, XNRM
      COMPLEX*16         AMAX
*     ..
*     .. Local Arrays ..
      INTEGER            DESCW( DLEN_ ), IDUM1( 1 ), IDUM2( 1 )
      DOUBLE PRECISION   RWORK( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            NUMROC
      DOUBLE PRECISION   PDLAMCH, PZLANGE
      EXTERNAL           LSAME, NUMROC, PDLAMCH, PZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, DGAMX2D, INFOG2L,
     $                   PXERBLA, PZMAX1, PZCOPY, PZGELQF,
     $                   PZGEQRF, PZLACGV, PZLACPY, PZLASCL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      PZQRT14 = ZERO
*
      IPWA = 1
      IROFFA = MOD( IA-1, DESCA( MB_ ) )
      ICOFFA = MOD( JA-1, DESCA( NB_ ) )
      IWA = IROFFA + 1
      JWA = ICOFFA + 1
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA,
     $              JJA, IAROW, IACOL )
      MPWA = NUMROC( M+IROFFA, DESCA( MB_ ), MYROW, IAROW, NPROW )
      NQWA = NUMROC( N+ICOFFA, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
*
      INFO = 0
      IF( LSAME( TRANS, 'N' ) ) THEN
         IF( N.LE.0 .OR. NRHS.LE.0 )
     $      RETURN
         TPSD = .FALSE.
         MPW = NUMROC( M+NRHS+IROFFA, DESCA( MB_ ), MYROW, IAROW,
     $                 NPROW )
         NQW = NQWA
*
*        Assign descriptor DESCW for workspace WORK and pointers to
*        matrices sub( A ) and sub( X ) in workspace
*
         IWX = IWA + M
         JWX = JWA
         LDW = MAX( 1, MPW )
         CALL DESCSET( DESCW, M+NRHS+IROFFA, N+ICOFFA, DESCA( MB_ ),
     $                 DESCA( NB_ ), IAROW, IACOL, ICTXT, LDW )
*
      ELSE IF( LSAME( TRANS, 'C' ) ) THEN
         IF( M.LE.0 .OR. NRHS.LE.0 )
     $      RETURN
         TPSD = .TRUE.
         MPW = MPWA
         NQW = NUMROC( N+NRHS+ICOFFA, DESCA( NB_ ), MYCOL, IACOL,
     $                 NPCOL )
*
*        Assign descriptor DESCW for workspace WORK and pointers to
*        matrices sub( A ) and sub( X ) in workspace
*
         IWX = IWA
         JWX = JWA + N
         LDW = MAX( 1, MPW )
         CALL DESCSET( DESCW, M+IROFFA, N+NRHS+ICOFFA, DESCA( MB_ ),
     $                 DESCA( NB_ ), IAROW, IACOL, ICTXT, LDW )
      ELSE
         CALL PXERBLA( ICTXT, 'PZQRT14', -1 )
         RETURN
      END IF
*
*     Copy and scale sub( A )
*
      IPTAU = IPWA + MPW*NQW
      CALL PZLACPY( 'All', M, N, A, IA, JA, DESCA, WORK( IPWA ), IWA,
     $              JWA, DESCW )
      RWORK( 1 ) = ZERO
      ANRM = PZLANGE( 'M', M, N, WORK( IPWA ), IWA, JWA, DESCW, RWORK )
      IF( ANRM.NE.ZERO )
     $   CALL PZLASCL( 'G', ANRM, ONE, M, N, WORK( IPWA ), IWA,
     $                 JWA, DESCW, INFO )
*
*     Copy sub( X ) or sub( X )' into the right place and scale it
*
      IF( TPSD ) THEN
*
*        Copy sub( X ) into columns jwa+n:jwa+n+nrhs-1 of work
*
         DO 10 J = 1, NRHS
            CALL PZCOPY( M, X, IX, JX+J-1, DESCX, 1, WORK( IPWA ), IWX,
     $                   JWX+J-1, DESCW, 1 )
   10    CONTINUE
         XNRM = PZLANGE( 'M', M, NRHS, WORK( IPWA ), IWX, JWX, DESCW,
     $                   RWORK )
         IF( XNRM.NE.ZERO )
     $      CALL PZLASCL( 'G', XNRM, ONE, M, NRHS, WORK( IPWA ), IWX,
     $                    JWX, DESCW, INFO )
*
*        Compute QR factorization of work(iwa:iwa+m-1,jwa:jwa+n+nrhs-1)
*
         MQW = NUMROC( M+ICOFFA, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
         IPW = IPTAU + MIN( MQW, NQW )
         LWORK = DESCW( NB_ ) * ( MPW + NQW + DESCW( NB_ ) )
         CALL PZGEQRF( M, N+NRHS, WORK( IPWA ), IWA, JWA, DESCW,
     $                WORK( IPTAU ), WORK( IPW ), LWORK, INFO )
*
*        Compute largest entry in upper triangle of
*        work(iwa+n:iwa+m-1,jwa+n:jwa+n+nrhs-1)
*
         ERR = ZERO
         IF( N.LT.M ) THEN
            DO 20 J = JWX, JWA+N+NRHS-1
               CALL PZMAX1( MIN(M-N,J-JWX+1), AMAX, IDUM, WORK( IPWA ),
     $                      IWA+N, J, DESCW, 1 )
               ERR = MAX( ERR, ABS( AMAX ) )
   20       CONTINUE
         END IF
         CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERR, 1, IDUM1, IDUM2,
     $                 -1, -1, 0 )
*
      ELSE
*
*        Copy sub( X )' into rows iwa+m:iwa+m+nrhs-1 of work
*
         DO 30 J = 1, NRHS
            CALL PZCOPY( N, X, IX, JX+J-1, DESCX, 1, WORK( IPWA ),
     $                   IWX+J-1, JWX, DESCW, DESCW( M_ ) )
            CALL PZLACGV( N, WORK( IPWA ), IWX+J-1, JWX, DESCW,
     $                    DESCW( M_ ) )
   30    CONTINUE
*
         XNRM = PZLANGE( 'M', NRHS, N, WORK( IPWA ), IWX, JWX, DESCW,
     $                   RWORK )
         IF( XNRM.NE.ZERO )
     $      CALL PZLASCL( 'G', XNRM, ONE, NRHS, N, WORK( IPWA ), IWX,
     $                    JWX, DESCW, INFO )
*
*        Compute LQ factorization of work(iwa:iwa+m+nrhs-1,jwa:jwa+n-1)
*
         NPW = NUMROC( N+IROFFA, DESCA( MB_ ), MYROW, IAROW, NPROW )
         IPW = IPTAU + MIN( MPW, NPW )
         LWORK = DESCW( MB_ ) * ( MPW + NQW + DESCW( MB_ ) )
         CALL PZGELQF( M+NRHS, N, WORK( IPWA ), IWA, JWA, DESCW,
     $                 WORK( IPTAU ), WORK( IPW ), LWORK, INFO )
*
*        Compute largest entry in lower triangle in
*        work(iwa+m:iwa+m+nrhs-1,jwa+m:jwa+n-1)
*
         ERR = ZERO
         DO 40 J = JWA+M, MIN( JWA+N-1, JWA+M+NRHS-1 )
            CALL PZMAX1( JWA+M+NRHS-J, AMAX, IDUM, WORK( IPWA ),
     $                   IWX+J-JWA-M, J, DESCW, 1 )
            ERR = MAX( ERR, ABS( AMAX ) )
   40    CONTINUE
         CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERR, 1, IDUM1, IDUM2,
     $                 -1, -1, 0 )
*
      END IF
*
      PZQRT14 = ERR / ( DBLE( MAX( M, N, NRHS ) ) *
     $          PDLAMCH( ICTXT, 'Epsilon' ) )
*
      RETURN
*
*     End of PZQRT14
*
      END
