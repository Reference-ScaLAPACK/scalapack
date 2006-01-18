      SUBROUTINE PZHENTRD( UPLO, N, A, IA, JA, DESCA, D, E, TAU, WORK,
     $                     LWORK, RWORK, LRWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     October 15, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, LRWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   D( * ), E( * ), RWORK( * )
      COMPLEX*16         A( * ), TAU( * ), WORK( * )
*     ..
*  Bugs
*  ====
*
*
*  Support for UPLO='U' is limited to calling the old, slow, PZHETRD
*  code.
*
*
*  Purpose
*
*  =======
*
*  PZHENTRD is a prototype version of PZHETRD which uses tailored
*  codes (either the serial, ZHETRD, or the parallel code, PZHETTRD)
*  when the workspace provided by the user is adequate.
*
*
*  PZHENTRD reduces a complex Hermitian matrix sub( A ) to Hermitian
*  tridiagonal form T by an unitary similarity transformation:
*  Q' * sub( A ) * Q = T, where sub( A ) = A(IA:IA+N-1,JA:JA+N-1).
*
*  Features
*  ========
*
*  PZHENTRD is faster than PZHETRD on almost all matrices,
*  particularly small ones (i.e. N < 500 * sqrt(P) ), provided that
*  enough workspace is available to use the tailored codes.
*
*  The tailored codes provide performance that is essentially
*  independent of the input data layout.
*
*  The tailored codes place no restrictions on IA, JA, MB or NB.
*  At present, IA, JA, MB and NB are restricted to those values allowed
*  by PZHETRD to keep the interface simple.  These restrictions are
*  documented below.  (Search for "restrictions".)
*
*  Notes
*  =====
*
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
*
*  Arguments
*  =========
*
*  UPLO    (global input) CHARACTER
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix sub( A ) is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          Hermitian distributed matrix sub( A ).  If UPLO = 'U', the
*          leading N-by-N upper triangular part of sub( A ) contains
*          the upper triangular part of the matrix, and its strictly
*          lower triangular part is not referenced. If UPLO = 'L', the
*          leading N-by-N lower triangular part of sub( A ) contains the
*          lower triangular part of the matrix, and its strictly upper
*          triangular part is not referenced. On exit, if UPLO = 'U',
*          the diagonal and first superdiagonal of sub( A ) are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements above the first superdiagonal,
*          with the array TAU, represent the unitary matrix Q as a
*          product of elementary reflectors; if UPLO = 'L', the diagonal
*          and first subdiagonal of sub( A ) are overwritten by the
*          corresponding elements of the tridiagonal matrix T, and the
*          elements below the first subdiagonal, with the array TAU,
*          represent the unitary matrix Q as a product of elementary
*          reflectors. See Further Details.
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
*  D       (local output) DOUBLE PRECISION array, dimension LOCc(JA+N-1)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i). D is tied to the distributed matrix A.
*
*  E       (local output) DOUBLE PRECISION array, dimension LOCc(JA+N-1)
*          if UPLO = 'U', LOCc(JA+N-2) otherwise. The off-diagonal
*          elements of the tridiagonal matrix T: E(i) = A(i,i+1) if
*          UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. E is tied to the
*          distributed matrix A.
*
*  TAU     (local output) COMPLEX*16, array, dimension
*          LOCc(JA+N-1). This array contains the scalar factors TAU of
*          the elementary reflectors. TAU is tied to the distributed
*          matrix A.
*
*  WORK    (local workspace/local output) COMPLEX*16 array,
*                                                  dimension (LWORK)
*          On exit, WORK( 1 ) returns the optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= MAX( NB * ( NP +1 ), 3 * NB )
*
*          For optimal performance, greater workspace is needed, i.e.
*            LWORK >= 2*( ANB+1 )*( 4*NPS+2 ) + ( NPS + 4 ) * NPS
*            ICTXT = DESCA( CTXT_ )
*            ANB = PJLAENV( ICTXT, 3, 'PZHETTRD', 'L', 0, 0, 0, 0 )
*            SQNPC = INT( SQRT( DBLE( NPROW * NPCOL ) ) )
*            NPS = MAX( NUMROC( N, 1, 0, 0, SQNPC ), 2*ANB )
*
*            NUMROC is a ScaLAPACK tool functions;
*            PJLAENV is a ScaLAPACK envionmental inquiry function
*            MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*            the subroutine BLACS_GRIDINFO.
*
*
*  RWORK   (local workspace/local output) COMPLEX*16 array,
*                                                  dimension (LRWORK)
*          On exit, RWORK( 1 ) returns the optimal LRWORK.
*
*  LRWORK  (local or global input) INTEGER
*          The dimension of the array RWORK.
*          LRWORK is local input and must be at least
*          LRWORK >= 1
*
*          For optimal performance, greater workspace is needed, i.e.
*            LRWORK >= MAX( 2 * N )
*
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(ia:ia+i-2,ja+i), and tau in TAU(ja+i-1).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in
*  A(ia+i+1:ia+n-1,ja+i-1), and tau in TAU(ja+i-1).
*
*  The contents of sub( A ) on exit are illustrated by the following
*  examples with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  d   e   v2  v3  v4 )              (  d                  )
*    (      d   e   v3  v4 )              (  e   d              )
*    (          d   e   v4 )              (  v1  e   d          )
*    (              d   e  )              (  v1  v2  e   d      )
*    (                  d  )              (  v1  v2  v3  e   d  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrix sub( A ) must verify some alignment proper-
*  ties, namely the following expression should be true:
*  ( MB_A.EQ.NB_A .AND. IROFFA.EQ.ICOFFA .AND. IROFFA.EQ.0 ) with
*  IROFFA = MOD( IA-1, MB_A ) and ICOFFA = MOD( JA-1, NB_A ).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      CHARACTER          COLCTOP, ROWCTOP
      INTEGER            ANB, CTXTB, I, IACOL, IAROW, ICOFFA, ICTXT,
     $                   IINFO, INDB, INDRD, INDRE, INDTAU, INDW, IPW,
     $                   IROFFA, J, JB, JX, K, KK, LLRWORK, LLWORK,
     $                   LRWMIN, LWMIN, MINSZ, MYCOL, MYCOLB, MYROW,
     $                   MYROWB, NB, NP, NPCOL, NPCOLB, NPROW, NPROWB,
     $                   NPS, NQ, ONEPMIN, ONEPRMIN, SQNPC, TTLRWMIN,
     $                   TTLWMIN
*     ..
*     .. Local Arrays ..
      INTEGER            DESCB( DLEN_ ), DESCW( DLEN_ ), IDUM1( 3 ),
     $                   IDUM2( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GET, BLACS_GRIDEXIT, BLACS_GRIDINFO,
     $                   BLACS_GRIDINIT, CHK1MAT, DESCSET, IGAMN2D,
     $                   PCHK1MAT, PDLAMR1D, PB_TOPGET, PB_TOPSET,
     $                   PXERBLA, PZELSET, PZHER2K, PZHETD2, PZHETTRD,
     $                   PZLAMR1D, PZLATRD, PZTRMR2D, ZHETRD
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2L, INDXG2P, NUMROC, PJLAENV
      EXTERNAL           LSAME, INDXG2L, INDXG2P, NUMROC, PJLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, ICHAR, INT, MAX, MIN, MOD, SQRT
*     ..
*     .. Executable Statements ..
*
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 600+CTXT_ )
      ELSE
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
         UPPER = LSAME( UPLO, 'U' )
         IF( INFO.EQ.0 ) THEN
            NB = DESCA( NB_ )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IAROW = INDXG2P( IA, NB, MYROW, DESCA( RSRC_ ), NPROW )
            IACOL = INDXG2P( JA, NB, MYCOL, DESCA( CSRC_ ), NPCOL )
            NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
            NQ = MAX( 1, NUMROC( N+JA-1, NB, MYCOL, DESCA( CSRC_ ),
     $           NPCOL ) )
            LWMIN = MAX( ( NP+1 )*NB, 3*NB )
            ANB = PJLAENV( ICTXT, 3, 'PZHETTRD', 'L', 0, 0, 0, 0 )
            MINSZ = PJLAENV( ICTXT, 5, 'PZHETTRD', 'L', 0, 0, 0, 0 )
            SQNPC = INT( SQRT( DBLE( NPROW*NPCOL ) ) )
            NPS = MAX( NUMROC( N, 1, 0, 0, SQNPC ), 2*ANB )
            TTLWMIN = 2*( ANB+1 )*( 4*NPS+2 ) + ( NPS+2 )*NPS
            LRWMIN = 1
            TTLRWMIN = 2*NPS
*
            WORK( 1 ) = DCMPLX( DBLE( TTLWMIN ) )
            RWORK( 1 ) = DBLE( TTLRWMIN )
            LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 )
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
*
*            The following two restrictions are not necessary provided
*            that either of the tailored codes are used.
*
            ELSE IF( IROFFA.NE.ICOFFA .OR. ICOFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 600+NB_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -11
            ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -13
            END IF
         END IF
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
         IDUM2( 2 ) = 11
         IF( LRWORK.EQ.-1 ) THEN
            IDUM1( 3 ) = -1
         ELSE
            IDUM1( 3 ) = 1
         END IF
         IDUM2( 3 ) = 13
         CALL PCHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, 3, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZHENTRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*
      ONEPMIN = N*N + 3*N + 1
      LLWORK = LWORK
      CALL IGAMN2D( ICTXT, 'A', ' ', 1, 1, LLWORK, 1, 1, -1, -1, -1,
     $              -1 )
*
      ONEPRMIN = 2*N
      LLRWORK = LRWORK
      CALL IGAMN2D( ICTXT, 'A', ' ', 1, 1, LLRWORK, 1, 1, -1, -1, -1,
     $              -1 )
*
*
*     Use the serial, LAPACK, code:  ZTRD on small matrices if we
*     we have enough space.
*
      NPROWB = 0
      IF( ( N.LT.MINSZ .OR. SQNPC.EQ.1 ) .AND. LLWORK.GE.ONEPMIN .AND.
     $    LLRWORK.GE.ONEPRMIN .AND. .NOT.UPPER ) THEN
         NPROWB = 1
         NPS = N
      ELSE
         IF( LLWORK.GE.TTLWMIN .AND. LLRWORK.GE.TTLRWMIN .AND. .NOT.
     $       UPPER ) THEN
            NPROWB = SQNPC
         END IF
      END IF
*
      IF( NPROWB.GE.1 ) THEN
         NPCOLB = NPROWB
         SQNPC = NPROWB
         INDB = 1
         INDRD = 1
         INDRE = INDRD + NPS
         INDTAU = INDB + NPS*NPS
         INDW = INDTAU + NPS
         LLWORK = LLWORK - INDW + 1
*
         CALL BLACS_GET( ICTXT, 10, CTXTB )
         CALL BLACS_GRIDINIT( CTXTB, 'Row major', SQNPC, SQNPC )
         CALL BLACS_GRIDINFO( CTXTB, NPROWB, NPCOLB, MYROWB, MYCOLB )
         CALL DESCSET( DESCB, N, N, 1, 1, 0, 0, CTXTB, NPS )
*
         CALL PZTRMR2D( UPLO, 'N', N, N, A, IA, JA, DESCA, WORK( INDB ),
     $                  1, 1, DESCB, ICTXT )
*
*
*        Only those processors in context CTXTB are needed for a while
*
         IF( NPROWB.GT.0 ) THEN
*
            IF( NPROWB.EQ.1 ) THEN
               CALL ZHETRD( UPLO, N, WORK( INDB ), NPS, RWORK( INDRD ),
     $                      RWORK( INDRE ), WORK( INDTAU ),
     $                      WORK( INDW ), LLWORK, INFO )
            ELSE
*
               CALL PZHETTRD( 'L', N, WORK( INDB ), 1, 1, DESCB,
     $                        RWORK( INDRD ), RWORK( INDRE ),
     $                        WORK( INDTAU ), WORK( INDW ), LLWORK,
     $                        INFO )
*
            END IF
         END IF
*
*           All processors participate in moving the data back to the
*           way that PZHENTRD expects it.
*
         CALL PDLAMR1D( N-1, RWORK( INDRE ), 1, 1, DESCB, E, 1, JA,
     $                  DESCA )
*
         CALL PDLAMR1D( N, RWORK( INDRD ), 1, 1, DESCB, D, 1, JA,
     $                  DESCA )
*
         CALL PZLAMR1D( N, WORK( INDTAU ), 1, 1, DESCB, TAU, 1, JA,
     $                  DESCA )
*
         CALL PZTRMR2D( UPLO, 'N', N, N, WORK( INDB ), 1, 1, DESCB, A,
     $                  IA, JA, DESCA, ICTXT )
*
         IF( MYROWB.GE.0 )
     $      CALL BLACS_GRIDEXIT( CTXTB )
*
      ELSE
*
         CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
         CALL PB_TOPGET( ICTXT, 'Combine', 'Rowwise', ROWCTOP )
         CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', '1-tree' )
         CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise', '1-tree' )
*
         IPW = NP*NB + 1
*
         IF( UPPER ) THEN
*
*        Reduce the upper triangle of sub( A ).
*
            KK = MOD( JA+N-1, NB )
            IF( KK.EQ.0 )
     $         KK = NB
            CALL DESCSET( DESCW, N, NB, NB, NB, IAROW,
     $                    INDXG2P( JA+N-KK, NB, MYCOL, DESCA( CSRC_ ),
     $                    NPCOL ), ICTXT, MAX( 1, NP ) )
*
            DO 10 K = N - KK + 1, NB + 1, -NB
               JB = MIN( N-K+1, NB )
               I = IA + K - 1
               J = JA + K - 1
*
*           Reduce columns I:I+NB-1 to tridiagonal form and form
*           the matrix W which is needed to update the unreduced part of
*           the matrix
*
               CALL PZLATRD( UPLO, K+JB-1, JB, A, IA, JA, DESCA, D, E,
     $                       TAU, WORK, 1, 1, DESCW, WORK( IPW ) )
*
*           Update the unreduced submatrix A(IA:I-1,JA:J-1), using an
*           update of the form:
*           A(IA:I-1,JA:J-1) := A(IA:I-1,JA:J-1) - V*W' - W*V'
*
               CALL PZHER2K( UPLO, 'No transpose', K-1, JB, -CONE, A,
     $                       IA, J, DESCA, WORK, 1, 1, DESCW, ONE, A,
     $                       IA, JA, DESCA )
*
*           Copy last superdiagonal element back into sub( A )
*
               JX = MIN( INDXG2L( J, NB, 0, IACOL, NPCOL ), NQ )
               CALL PZELSET( A, I-1, J, DESCA, DCMPLX( E( JX ) ) )
*
               DESCW( CSRC_ ) = MOD( DESCW( CSRC_ )+NPCOL-1, NPCOL )
*
   10       CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
            CALL PZHETD2( UPLO, MIN( N, NB ), A, IA, JA, DESCA, D, E,
     $                    TAU, WORK, LWORK, IINFO )
*
         ELSE
*
*        Reduce the lower triangle of sub( A )
*
            KK = MOD( JA+N-1, NB )
            IF( KK.EQ.0 )
     $         KK = NB
            CALL DESCSET( DESCW, N, NB, NB, NB, IAROW, IACOL, ICTXT,
     $                    MAX( 1, NP ) )
*
            DO 20 K = 1, N - NB, NB
               I = IA + K - 1
               J = JA + K - 1
*
*           Reduce columns I:I+NB-1 to tridiagonal form and form
*           the matrix W which is needed to update the unreduced part
*           of the matrix
*
               CALL PZLATRD( UPLO, N-K+1, NB, A, I, J, DESCA, D, E, TAU,
     $                       WORK, K, 1, DESCW, WORK( IPW ) )
*
*           Update the unreduced submatrix A(I+NB:IA+N-1,I+NB:IA+N-1),
*           using an update of the form: A(I+NB:IA+N-1,I+NB:IA+N-1) :=
*           A(I+NB:IA+N-1,I+NB:IA+N-1) - V*W' - W*V'
*
               CALL PZHER2K( UPLO, 'No transpose', N-K-NB+1, NB, -CONE,
     $                       A, I+NB, J, DESCA, WORK, K+NB, 1, DESCW,
     $                       ONE, A, I+NB, J+NB, DESCA )
*
*           Copy last subdiagonal element back into sub( A )
*
               JX = MIN( INDXG2L( J+NB-1, NB, 0, IACOL, NPCOL ), NQ )
               CALL PZELSET( A, I+NB, J+NB-1, DESCA, DCMPLX( E( JX ) ) )
*
               DESCW( CSRC_ ) = MOD( DESCW( CSRC_ )+1, NPCOL )
*
   20       CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
            CALL PZHETD2( UPLO, KK, A, IA+K-1, JA+K-1, DESCA, D, E, TAU,
     $                    WORK, LWORK, IINFO )
         END IF
*
         CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
         CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise', ROWCTOP )
*
      END IF
*
      WORK( 1 ) = DCMPLX( DBLE( TTLWMIN ) )
      RWORK( 1 ) = DBLE( TTLRWMIN )
*
      RETURN
*
*     End of PZHENTRD
*
      END
