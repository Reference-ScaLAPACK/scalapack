      SUBROUTINE PDGEHD2( N, ILO, IHI, A, IA, JA, DESCA, TAU, WORK,
     $                    LWORK, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER             IA, IHI, ILO, INFO, JA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER             DESCA( * )
      DOUBLE PRECISION    A( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDGEHD2 reduces a real general distributed matrix sub( A )
*  to upper Hessenberg form H by an orthogonal similarity transforma-
*  tion:  Q' * sub( A ) * Q = H, where
*  sub( A ) = A(IA+N-1:IA+N-1,JA+N-1:JA+N-1).
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
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  ILO     (global input) INTEGER
*  IHI     (global input) INTEGER
*          It is assumed that sub( A ) is already upper triangular in
*          rows IA:IA+ILO-2 and IA+IHI:IA+N-1 and columns JA:JA+JLO-2
*          and JA+JHI:JA+N-1. See Further Details. If N > 0,
*          1 <= ILO <= IHI <= N; otherwise set ILO = 1, IHI = N.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the N-by-N
*          general distributed matrix sub( A ) to be reduced. On exit,
*          the upper triangle and the first subdiagonal of sub( A ) are
*          overwritten with the upper Hessenberg matrix H, and the ele-
*          ments below the first subdiagonal, with the array TAU, repre-
*          sent the orthogonal matrix Q as a product of elementary
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
*  TAU     (local output) DOUBLE PRECISION array, dimension LOCc(JA+N-2)
*          The scalar factors of the elementary reflectors (see Further
*          Details). Elements JA:JA+ILO-2 and JA+IHI:JA+N-2 of TAU are
*          set to zero. TAU is tied to the distributed matrix A.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*                                                    dimension (LWORK)
*          On exit, WORK( 1 ) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK >= NB + MAX( NpA0, NB )
*
*          where NB = MB_A = NB_A, IROFFA = MOD( IA-1, NB ),
*          IAROW = INDXG2P( IA, NB, MYROW, RSRC_A, NPROW ),
*          NpA0 = NUMROC( IHI+IROFFA, NB, MYROW, IAROW, NPROW ),
*
*          INDXG2P and NUMROC are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  INFO    (local output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of (ihi-ilo) elementary
*  reflectors
*
*     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
*  exit in A(ia+ilo+i:ia+ihi-1,ja+ilo+i-2), and tau in TAU(ja+ilo+i-2).
*
*  The contents of A(IA:IA+N-1,JA:JA+N-1) are illustrated by the follo-
*  wing example, with n = 7, ilo = 2 and ihi = 6:
*
*  on entry                         on exit
*
*  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
*  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
*  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
*  (                         a )    (                          a )
*
*  where a denotes an element of the original matrix sub( A ), h denotes
*  a modified element of the upper Hessenberg matrix H, and vi denotes
*  an element of the vector defining H(ja+ilo+i-2).
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrix sub( A ) must verify some alignment proper-
*  ties, namely the following expression should be true:
*  ( MB_A.EQ.NB_A .AND. IROFFA.EQ.ICOFFA )
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IAROW, ICOFFA, ICTXT, IROFFA, J, K, LWMIN,
     $                   MYCOL, MYROW, NPA0, NPCOL, NPROW
      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GRIDINFO, CHK1MAT, PDELSET,
     $                   PDLARF, PDLARFG, PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, MOD
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
         INFO = -(700+CTXT_)
      ELSE
         CALL CHK1MAT( N, 1, N, 1, IA, JA, DESCA, 7, INFO )
         IF( INFO.EQ.0 ) THEN
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            NPA0 = NUMROC( IHI+IROFFA, DESCA( MB_ ), MYROW, IAROW,
     $                     NPROW )
            LWMIN = DESCA( NB_ ) + MAX( NPA0, DESCA( NB_ ) )
*
            WORK( 1 ) = DBLE( LWMIN )
            LQUERY = ( LWORK.EQ.-1 )
            IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
               INFO = -2
            ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
               INFO = -3
            ELSE IF( IROFFA.NE.ICOFFA ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(700+NB_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -10
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDGEHD2', -INFO )
         CALL BLACS_ABORT( ICTXT, 1 )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      DO 10 K = ILO, IHI-1
         I = IA + K - 1
         J = JA + K - 1
*
*        Compute elementary reflector H(j) to annihilate
*        A(i+2:ihi+ia-1,j)
*
         CALL PDLARFG( IHI-K, AII, I+1, J, A, MIN( I+2, N+IA-1 ), J,
     $                 DESCA, 1, TAU )
         CALL PDELSET( A, I+1, J, DESCA, ONE )
*
*        Apply H(k) to A(ia:ihi+ia-1,j+1:ihi+ja-1) from the right
*
         CALL PDLARF( 'Right', IHI, IHI-K, A, I+1, J, DESCA, 1, TAU, A,
     $                IA, J+1, DESCA, WORK )
*
*        Apply H(j) to A(i+1:ia+ihi-1,j+1:ja+n-1) from the left
*
         CALL PDLARF( 'Left', IHI-K, N-K, A, I+1, J, DESCA, 1, TAU, A,
     $                I+1, J+1, DESCA, WORK )
*
         CALL PDELSET( A, I+1, J, DESCA, AII )
   10 CONTINUE
*
      WORK( 1 ) = DBLE( LWMIN )
*
      RETURN
*
*     End of PDGEHD2
*
      END
