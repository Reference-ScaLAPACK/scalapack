      SUBROUTINE PSTRCON( NORM, UPLO, DIAG, N, A, IA, JA, DESCA, RCOND,
     $                    WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            IA, JA, INFO, LIWORK, LWORK, N
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), IWORK( * )
      REAL               A( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSTRCON estimates the reciprocal of the condition number of a
*  triangular distributed matrix A(IA:IA+N-1,JA:JA+N-1), in either the
*  1-norm or the infinity-norm.
*
*  The norm of A(IA:IA+N-1,JA:JA+N-1) is computed and an estimate is
*  obtained for norm(inv(A(IA:IA+N-1,JA:JA+N-1))), then the reciprocal
*  of the condition number is computed as
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
*  NORM    (global input) CHARACTER
*          Specifies whether the 1-norm condition number or the
*          infinity-norm condition number is required:
*          = '1' or 'O':  1-norm;
*          = 'I':         Infinity-norm.
*
*  UPLO    (global input) CHARACTER
*          = 'U':  A(IA:IA+N-1,JA:JA+N-1) is upper triangular;
*          = 'L':  A(IA:IA+N-1,JA:JA+N-1) is lower triangular.
*
*  DIAG    (global input) CHARACTER
*          = 'N':  A(IA:IA+N-1,JA:JA+N-1) is non-unit triangular;
*          = 'U':  A(IA:IA+N-1,JA:JA+N-1) is unit triangular.
*
*  N       (global input) INTEGER
*          The order of the distributed matrix A(IA:IA+N-1,JA:JA+N-1).
*          N >= 0.
*
*  A       (local input) REAL pointer into the local memory
*          to an array of dimension ( LLD_A, LOCc(JA+N-1) ). This array
*          contains the local pieces of the triangular distributed
*          matrix A(IA:IA+N-1,JA:JA+N-1). If UPLO = 'U', the leading
*          N-by-N upper triangular part of this distributed matrix con-
*          tains the upper triangular matrix, and its strictly lower
*          triangular part is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of this ditributed
*          matrix contains the lower triangular matrix, and the strictly
*          upper triangular part is not referenced. If DIAG = 'U', the
*          diagonal elements of A(IA:IA+N-1,JA:JA+N-1) are also not
*          referenced and are assumed to be 1.
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
*  RCOND   (global output) REAL
*          The reciprocal of the condition number of the distributed
*          matrix A(IA:IA+N-1,JA:JA+N-1), computed as
*             RCOND = 1 / ( norm( A(IA:IA+N-1,JA:JA+N-1)      ) *
*                           norm( inv(A(IA:IA+N-1,JA:JA+N-1)) ) ).
*
*  WORK    (local workspace/local output) REAL array,
*                                                   dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= 2*LOCr(N+MOD(IA-1,MB_A)) + LOCc(N+MOD(JA-1,NB_A))
*          + MAX( 2, MAX( NB_A*MAX( 1, CEIL(NPROW-1,NPCOL) ),
*                         LOCc(N+MOD(JA-1,NB_A)) +
*                         NB_A*MAX( 1, CEIL(NPCOL-1,NPROW) ) ).
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  IWORK   (local workspace/local output) INTEGER array,
*                                                   dimension (LIWORK)
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
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, NOUNIT, ONENRM, UPPER
      CHARACTER          CBTOP, COLCTOP, NORMIN, ROWCTOP
      INTEGER            IACOL, IAROW, ICOFF, ICTXT, IIA, IPN, IPV, IPW,
     $                   IPX, IROFF, IV, IX, IXX, JJA, JV, JX, KASE,
     $                   KASE1, LIWMIN, LWMIN, MYCOL, MYROW, NP, NPCOL,
     $                   NPMOD, NPROW, NQ, NQMOD
      REAL               AINVNM, ANORM, SCALE, SMLNUM
      REAL               WMAX
*     ..
*     .. Local Arrays ..
      INTEGER            DESCV( DLEN_ ), DESCX( DLEN_ ), IDUM1( 5 ),
     $                   IDUM2( 5 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, INFOG2L,
     $                   PCHK1MAT, PSAMAX, PSLATRS, PSLACON,
     $                   PSRSCL, PB_TOPGET, PB_TOPSET, PXERBLA, SGEBR2D,
     $                   SGEBS2D
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, INDXG2P, NUMROC
      REAL               PSLAMCH, PSLANTR
      EXTERNAL           ICEIL, INDXG2P, LSAME, NUMROC, PSLAMCH,
     $                   PSLANTR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, ICHAR, MAX, MOD, REAL
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
         CALL CHK1MAT( N, 4, N, 4, IA, JA, DESCA, 8, INFO )
         IF( INFO.EQ.0 ) THEN
            UPPER = LSAME( UPLO, 'U' )
            ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
            NOUNIT = LSAME( DIAG, 'N' )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IACOL = INDXG2P( JA, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            NPMOD = NUMROC( N + MOD( IA-1, DESCA( MB_ ) ), DESCA( MB_ ),
     $                      MYROW, IAROW, NPROW )
            NQMOD = NUMROC( N + MOD( JA-1, DESCA( NB_ ) ), DESCA( NB_ ),
     $                      MYCOL, IACOL, NPCOL )
            LWMIN = 2*NPMOD + NQMOD +
     $              MAX( 2, MAX( DESCA( NB_ )*
     $                   MAX( 1, ICEIL( NPROW-1, NPCOL ) ), NQMOD +
     $                   DESCA( NB_ )*
     $                   MAX( 1, ICEIL( NPCOL-1, NPROW ) ) ) )
            WORK( 1 ) = REAL( LWMIN )
            LIWMIN = NPMOD
            IWORK( 1 ) = LIWMIN
            LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
            IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -2
            ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
               INFO = -3
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -11
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -13
            END IF
         END IF
*
         IF( ONENRM ) THEN
            IDUM1( 1 ) = ICHAR( '1' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'I' )
         END IF
         IDUM2( 1 ) = 1
         IF( UPPER ) THEN
            IDUM1( 2 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'L' )
         END IF
         IDUM2( 2 ) = 2
         IF( NOUNIT ) THEN
            IDUM1( 3 ) = ICHAR( 'N' )
         ELSE
            IDUM1( 3 ) = ICHAR( 'U' )
         END IF
         IDUM2( 3 ) = 3
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 4 ) = -1
         ELSE
            IDUM1( 4 ) = 1
         END IF
         IDUM2( 4 ) = 11
         IF( LIWORK.EQ.-1 ) THEN
            IDUM1( 5 ) = -1
         ELSE
            IDUM1( 5 ) = 1
         END IF
         IDUM2( 5 ) = 13
         CALL PCHK1MAT( N, 4, N, 4, IA, JA, DESCA, 8, 5, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PSTRCON', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      END IF
*
      CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
      CALL PB_TOPGET( ICTXT, 'Combine', 'Rowwise',    ROWCTOP )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', '1-tree' )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise',    '1-tree' )
*
      RCOND = ZERO
      SMLNUM = PSLAMCH( ICTXT, 'Safe minimum' )*REAL( MAX( 1, N ) )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      IROFF = MOD( IA-1, DESCA( MB_ ) )
      ICOFF = MOD( JA-1, DESCA( NB_ ) )
      NP = NUMROC( N+IROFF, DESCA( MB_ ), MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFF, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IV = IROFF + 1
      IX = IV
      JV = ICOFF + 1
      JX = JV
*
      IPX  = 1
      IPV  = IPX + NP
      IPN = IPV + NP
      IPW  = IPN + NQ
*
      CALL DESCSET( DESCV, N+IROFF, 1, DESCA( MB_ ), 1, IAROW, MYCOL,
     $              ICTXT, MAX( 1, NP ) )
      CALL DESCSET( DESCX, N+IROFF, 1, DESCA( MB_ ), 1, IAROW, MYCOL,
     $              ICTXT, MAX( 1, NP ) )
*
*     Compute the norm of the triangular matrix A.
*
      ANORM = PSLANTR( NORM, UPLO, DIAG, N, N, A, IA, JA, DESCA, WORK )
*
*     Continue only if ANORM > 0.
*
      IF( ANORM.GT.ZERO ) THEN
*
*        Estimate the norm of the inverse of A.
*
         AINVNM = ZERO
         NORMIN = 'N'
         IF( ONENRM ) THEN
            KASE1 = 1
         ELSE
            KASE1 = 2
         END IF
         KASE = 0
   10    CONTINUE
         CALL PSLACON( N, WORK( IPV ), IV, JV, DESCV, WORK( IPX ),
     $                 IX, JX, DESCX, IWORK, AINVNM, KASE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.KASE1 ) THEN
*
*              Multiply by inv(A).
*
               DESCX( CSRC_ ) = IACOL
               CALL PSLATRS( UPLO, 'No transpose', DIAG, NORMIN,
     $                       N, A, IA, JA, DESCA, WORK( IPX ), IX, JX,
     $                       DESCX, SCALE, WORK( IPN ), WORK( IPW ) )
               DESCX( CSRC_ ) = MYCOL
            ELSE
*
*              Multiply by inv(A').
*
               DESCX( CSRC_ ) = IACOL
               CALL PSLATRS( UPLO, 'Transpose', DIAG, NORMIN,
     $                       N, A, IA, JA, DESCA, WORK( IPX ), IX, JX,
     $                       DESCX, SCALE, WORK( IPN ), WORK( IPW ) )
               DESCX( CSRC_ ) = MYCOL
            END IF
            NORMIN = 'Y'
*
*           Multiply by 1/SCALE if doing so will not cause overflow.
*
            IF( SCALE.NE.ONE ) THEN
               CALL PSAMAX( N, WMAX, IXX, WORK( IPX ), IX, JX,
     $                   DESCX, 1 )
               IF( DESCX( M_ ).EQ.1 .AND. N.EQ.1 ) THEN
                  CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise',
     $                          CBTOP )
                  IF( MYROW.EQ.IAROW ) THEN
                     CALL SGEBS2D( ICTXT, 'Column', CBTOP, 1, 1, WMAX,
     $                             1 )
                  ELSE
                     CALL SGEBR2D( ICTXT, 'Column', CBTOP, 1, 1, WMAX,
     $                             1, IAROW, MYCOL )
                  END IF
               END IF
               IF( SCALE.LT.ABS( WMAX )*SMLNUM .OR. SCALE.EQ.ZERO )
     $            GO TO 20
               CALL PSRSCL( N, SCALE, WORK( IPX ), IX, JX, DESCX, 1 )
            END IF
            GO TO 10
         END IF
*
*        Compute the estimate of the reciprocal condition number.
*
         IF( AINVNM.NE.ZERO )
     $      RCOND = ( ONE / ANORM ) / AINVNM
      END IF
*
   20 CONTINUE
*
      CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise',    ROWCTOP )
*
      RETURN
*
*     End of PSTRCON
*
      END
