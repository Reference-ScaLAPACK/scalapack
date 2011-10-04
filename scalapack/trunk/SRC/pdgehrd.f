      SUBROUTINE PDGEHRD( N, ILO, IHI, A, IA, JA, DESCA, TAU, WORK,
     $                    LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
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
*  PDGEHRD reduces a real general distributed matrix sub( A )
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
*          rows IA:IA+ILO-2 and IA+IHI:IA+N-1 and columns JA:JA+ILO-2
*          and JA+IHI:JA+N-1. See Further Details. If N > 0,
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
*          LWORK >= NB*NB + NB*MAX( IHIP+1, IHLP+INLQ )
*
*          where NB = MB_A = NB_A, IROFFA = MOD( IA-1, NB ),
*          ICOFFA = MOD( JA-1, NB ), IOFF = MOD( IA+ILO-2, NB ),
*          IAROW = INDXG2P( IA, NB, MYROW, RSRC_A, NPROW ),
*          IHIP = NUMROC( IHI+IROFFA, NB, MYROW, IAROW, NPROW ),
*          ILROW = INDXG2P( IA+ILO-1, NB, MYROW, RSRC_A, NPROW ),
*          IHLP = NUMROC( IHI-ILO+IOFF+1, NB, MYROW, ILROW, NPROW ),
*          ILCOL = INDXG2P( JA+ILO-1, NB, MYCOL, CSRC_A, NPCOL ),
*          INLQ = NUMROC( N-ILO+IOFF+1, NB, MYCOL, ILCOL, NPCOL ),
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
*  v(1:I) = 0, v(I+1) = 1 and v(IHI+1:N) = 0; v(I+2:IHI) is stored on
*  exit in A(IA+ILO+I:IA+IHI-1,JA+ILO+I-2), and tau in TAU(JA+ILO+I-2).
*
*  The contents of A(IA:IA+N-1,JA:JA+N-1) are illustrated by the follow-
*  ing example, with N = 7, ILO = 2 and IHI = 6:
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
*  where a denotes an element of the original matrix sub( A ), H denotes
*  a modified element of the upper Hessenberg matrix H, and vi denotes
*  an element of the vector defining H(JA+ILO+I-2).
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
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      CHARACTER          COLCTOP, ROWCTOP
      INTEGER            I, IACOL, IAROW, IB, ICOFFA, ICTXT, IHIP,
     $                   IHLP, IIA, IINFO, ILCOL, ILROW, IMCOL, INLQ,
     $                   IOFF, IPT, IPW, IPY, IROFFA, J, JJ, JJA, JY,
     $                   K, L, LWMIN, MYCOL, MYROW, NB, NPCOL, NPROW,
     $                   NQ
      DOUBLE PRECISION   EI
*     ..
*     .. Local Arrays ..
      INTEGER            DESCY( DLEN_ ), IDUM1( 3 ), IDUM2( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, INFOG1L,
     $                   INFOG2L, PCHK1MAT, PDGEMM, PDGEHD2,
     $                   PDLAHRD, PDLARFB, PB_TOPGET, PB_TOPSET, PXERBLA
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
            NB = DESCA( NB_ )
            IROFFA = MOD( IA-1, NB )
            ICOFFA = MOD( JA-1, NB )
            CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIA, JJA, IAROW, IACOL )
            IHIP = NUMROC( IHI+IROFFA, NB, MYROW, IAROW, NPROW )
            IOFF = MOD( IA+ILO-2, NB )
            ILROW = INDXG2P( IA+ILO-1, NB, MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IHLP = NUMROC( IHI-ILO+IOFF+1, NB, MYROW, ILROW, NPROW )
            ILCOL = INDXG2P( JA+ILO-1, NB, MYCOL, DESCA( CSRC_ ),
     $                       NPCOL )
            INLQ = NUMROC( N-ILO+IOFF+1, NB, MYCOL, ILCOL, NPCOL )
            LWMIN = NB*( NB + MAX( IHIP+1, IHLP+INLQ ) )
*
            WORK( 1 ) = DBLE( LWMIN )
            LQUERY = ( LWORK.EQ.-1 )
            IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
               INFO = -2
            ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
               INFO = -3
C            ELSE IF( IROFFA.NE.ICOFFA .OR. IROFFA.NE.0 ) THEN
            ELSE IF( IROFFA.NE.ICOFFA ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(700+NB_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -10
            END IF
         END IF
         IDUM1( 1 ) = ILO
         IDUM2( 1 ) = 2
         IDUM1( 2 ) = IHI
         IDUM2( 2 ) = 3
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 3 ) = -1
         ELSE
            IDUM1( 3 ) = 1
         END IF
         IDUM2( 3 ) = 10
         CALL PCHK1MAT( N, 1, N, 1, IA, JA, DESCA, 7, 3, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDGEHRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Set elements JA:JA+ILO-2 and JA+JHI-1:JA+N-2 of TAU to zero.
*
      NQ = NUMROC( JA+N-2, NB, MYCOL, DESCA( CSRC_ ), NPCOL )
      CALL INFOG1L( JA+ILO-2, NB, NPCOL, MYCOL, DESCA( CSRC_ ), JJ,
     $              IMCOL )
      DO 10 J = JJA, MIN( JJ, NQ )
         TAU( J ) = ZERO
   10 CONTINUE
*
      CALL INFOG1L( JA+IHI-1, NB, NPCOL, MYCOL, DESCA( CSRC_ ), JJ,
     $              IMCOL )
      DO 20 J = JJ, NQ
         TAU( J ) = ZERO
   20 CONTINUE
*
*     Quick return if possible
*
      IF( IHI-ILO.LE.0 )
     $   RETURN
*
      CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
      CALL PB_TOPGET( ICTXT, 'Combine', 'Rowwise',    ROWCTOP )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', '1-tree' )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise',    '1-tree' )
*
      IPT = 1
      IPY = IPT + NB * NB
      IPW = IPY + IHIP * NB
      CALL DESCSET( DESCY, IHI+IROFFA, NB, NB, NB, IAROW, ILCOL, ICTXT,
     $              MAX( 1, IHIP ) )
*
      K = ILO
      IB = NB - IOFF
      JY = IOFF + 1
*
*     Loop over remaining block of columns
*
      DO 30 L = 1, IHI-ILO+IOFF-NB, NB
         I = IA + K - 1
         J = JA + K - 1
*
*        Reduce columns j:j+ib-1 to Hessenberg form, returning the
*        matrices V and T of the block reflector H = I - V*T*V'
*        which performs the reduction, and also the matrix Y = A*V*T
*
         CALL PDLAHRD( IHI, K, IB, A, IA, J, DESCA, TAU, WORK( IPT ),
     $                 WORK( IPY ), 1, JY, DESCY, WORK( IPW ) )
*
*        Apply the block reflector H to A(ia:ia+ihi-1,j+ib:ja+ihi-1)
*        from the right, computing  A := A - Y * V'.
*        V(i+ib,ib-1) must be set to 1.
*
         CALL PDELSET2( EI, A, I+IB, J+IB-1, DESCA, ONE )
         CALL PDGEMM( 'No transpose', 'Transpose', IHI, IHI-K-IB+1, IB,
     $                -ONE, WORK( IPY ), 1, JY, DESCY, A, I+IB, J,
     $                DESCA, ONE, A, IA, J+IB, DESCA )
         CALL PDELSET( A, I+IB, J+IB-1, DESCA, EI )
*
*        Apply the block reflector H to A(i+1:ia+ihi-1,j+ib:ja+n-1) from
*        the left
*
         CALL PDLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise',
     $                 IHI-K, N-K-IB+1, IB, A, I+1, J, DESCA,
     $                 WORK( IPT ), A, I+1, J+IB, DESCA, WORK( IPY ) )
*
         K = K + IB
         IB = NB
         JY = 1
         DESCY( CSRC_ ) = MOD( DESCY( CSRC_ ) + 1, NPCOL )
*
   30 CONTINUE
*
*     Use unblocked code to reduce the rest of the matrix
*
      CALL PDGEHD2( N, K, IHI, A, IA, JA, DESCA, TAU, WORK, LWORK,
     $              IINFO )
*
      CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
      CALL PB_TOPSET( ICTXT, 'Combine', 'Rowwise',    ROWCTOP )
*
      WORK( 1 ) = DBLE( LWMIN )
*
      RETURN
*
*     End of PDGEHRD
*
      END
