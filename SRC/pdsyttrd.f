      SUBROUTINE PDSYTTRD( UPLO, N, A, IA, JA, DESCA, D, E, TAU, WORK,
     $                     LWORK, INFO )
*
*  -- ScaLAPACK routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), D( * ), E( * ), TAU( * ), WORK( * )
*     ..
*
*     Purpose
*
*     =======
*
*     PDSYTTRD reduces a complex Hermitian matrix sub( A ) to Hermitian
*     tridiagonal form T by an unitary similarity transformation:
*     Q' * sub( A ) * Q = T, where sub( A ) = A(IA:IA+N-1,JA:JA+N-1).
*
*     Notes
*     =====
*
*     Each global data object is described by an associated description
*     vector.  This vector stores the information required to establish
*     the mapping between an object element and its corresponding
*     process and memory location.
*
*     Let A be a generic term for any 2D block cyclicly distributed
*     array.
*     Such a global array has an associated description vector DESCA.
*     In the following comments, the character _ should be read as
*     "of the global array".
*
*     NOTATION        STORED IN      EXPLANATION
*     --------------- -------------- -----------------------------------
*     DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*     DTYPE_A = 1.
*     CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle,
*     indicating the BLACS process grid A is distribu-
*     ted over. The context itself is glo-
*     bal, but the handle (the integer
*     value) may vary.
*     M_A    (global) DESCA( M_ )    The number of rows in the global
*     array A.
*     N_A    (global) DESCA( N_ )    The number of columns in the global
*     array A.
*     MB_A   (global) DESCA( MB_ )   The blocking factor used to
*     distribute the rows of the array.
*     NB_A   (global) DESCA( NB_ )   The blocking factor used to
*     distribute the columns of the array.
*     RSRC_A (global) DESCA( RSRC_ ) The process row over which the
*     first row of the array A is distributed.
*     CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*     first column of the array A is
*     distributed.
*     LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*     array.  LLD_A >= MAX(1,LOCp(M_A)).
*
*     Let K be the number of rows or columns of a distributed matrix,
*     and assume that its process grid has dimension p x q.
*     LOCp( K ) denotes the number of elements of K that a process
*     would receive if K were distributed over the p processes of its
*     process column.
*     Similarly, LOCq( K ) denotes the number of elements of K that a
*     process would receive if K were distributed over the q processes
*     of its process row.
*     The values of LOCp() and LOCq() may be determined via a call to
*     the ScaLAPACK tool function, NUMROC:
*     LOCp( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*     LOCq( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*     An upper bound for these quantities may be computed by:
*     LOCp( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*     LOCq( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*     Arguments
*     =========
*
*     UPLO    (global input) CHARACTER
*     Specifies whether the upper or lower triangular part of the
*     Hermitian matrix sub( A ) is stored:
*     = 'U':  Upper triangular
*     = 'L':  Lower triangular
*
*     N       (global input) INTEGER
*     The number of rows and columns to be operated on, i.e. the
*     order of the distributed submatrix sub( A ). N >= 0.
*
*     A  (local input/local output) DOUBLE PRECISION pointer into the
*     local memory to an array of dimension (LLD_A,LOCq(JA+N-1)).
*     On entry, this array contains the local pieces of the
*     Hermitian distributed matrix sub( A ).  If UPLO = 'U', the
*     leading N-by-N upper triangular part of sub( A ) contains
*     the upper triangular part of the matrix, and its strictly
*     lower triangular part is not referenced. If UPLO = 'L', the
*     leading N-by-N lower triangular part of sub( A ) contains the
*     lower triangular part of the matrix, and its strictly upper
*     triangular part is not referenced. On exit, if UPLO = 'U',
*     the diagonal and first superdiagonal of sub( A ) are over-
*     written by the corresponding elements of the tridiagonal
*     matrix T, and the elements above the first superdiagonal,
*     with the array TAU, represent the unitary matrix Q as a
*     product of elementary reflectors; if UPLO = 'L', the diagonal
*     and first subdiagonal of sub( A ) are overwritten by the
*     corresponding elements of the tridiagonal matrix T, and the
*     elements below the first subdiagonal, with the array TAU,
*     represent the unitary matrix Q as a product of elementary
*     reflectors. See Further Details.
*
*     IA      (global input) INTEGER
*     The row index in the global array A indicating the first
*     row of sub( A ).
*
*     JA      (global input) INTEGER
*     The column index in the global array A indicating the
*     first column of sub( A ).
*
*     DESCA   (global and local input) INTEGER array of dimension DLEN_.
*     The array descriptor for the distributed matrix A.
*
*     D       (local output) DOUBLE PRECISION array, dim LOCq(JA+N-1)
*     The diagonal elements of the tridiagonal matrix T:
*     D(i) = A(i,i). D is tied to the distributed matrix A.
*
*     E       (local output) DOUBLE PRECISION array, dim LOCq(JA+N-1)
*     if UPLO = 'U', LOCq(JA+N-2) otherwise. The off-diagonal
*     elements of the tridiagonal matrix T: E(i) = A(i,i+1) if
*     UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. E is tied to the
*     distributed matrix A.
*
*     TAU     (local output) DOUBLE PRECISION array, dimension
*     LOCq(JA+N-1). This array contains the scalar factors TAU of
*     the elementary reflectors. TAU is tied to the distributed
*     matrix A.
*
*     WORK  (local workspace) DOUBLE PRECISION array, dimension (LWORK)
*     On exit, WORK( 1 ) returns the minimal and optimal workspace
*
*     LWORK   (local input) INTEGER
*     The dimension of the array WORK.
*     LWORK >= 2*( ANB+1 )*( 4*NPS+2 ) + NPS
*     Where:
*         NPS = MAX( NUMROC( N, 1, 0, 0, NPROW ), 2*ANB )
*         ANB = PJLAENV( DESCA( CTXT_ ), 3, 'PDSYTTRD', 'L', 0, 0,
*           0, 0 )
*
*         NUMROC is a ScaLAPACK tool function;
*         PJLAENV is a ScaLAPACK envionmental inquiry function
*         MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*         the subroutine BLACS_GRIDINFO.
*
*     INFO    (global output) INTEGER
*     = 0:  successful exit
*     < 0:  If the i-th argument is an array and the j-entry had
*     an illegal value, then INFO = -(i*100+j), if the i-th
*     argument is a scalar and had an illegal value, then
*     INFO = -i.
*
*     Further Details
*     ===============
*
*     If UPLO = 'U', the matrix Q is represented as a product of
*     elementary reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*     Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*     where tau is a complex scalar, and v is a complex vector with
*     v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*     A(ia:ia+i-2,ja+i), and tau in TAU(ja+i-1).
*
*     If UPLO = 'L', the matrix Q is represented as a product of
*     elementary reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*     Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*     where tau is a complex scalar, and v is a complex vector with
*     v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in
*     A(ia+i+1:ia+n-1,ja+i-1), and tau in TAU(ja+i-1).
*
*     The contents of sub( A ) on exit are illustrated by the following
*     examples with n = 5:
*
*     if UPLO = 'U':                       if UPLO = 'L':
*
*     (  d   e   v2  v3  v4 )              (  d                  )
*     (      d   e   v3  v4 )              (  e   d              )
*     (          d   e   v4 )              (  v1  e   d          )
*     (              d   e  )              (  v1  v2  e   d      )
*     (                  d  )              (  v1  v2  v3  e   d  )
*
*     where d and e denote diagonal and off-diagonal elements of T, and
*     vi denotes an element of the vector defining H(i).
*
*     Data storage requirements
*     =========================
*
*     PDSYTTRD is not intended to be called directly.  All users are
*     encourage to call PDSYTRD which will then call PDHETTRD if
*     appropriate.  A must be in cyclic format (i.e. MB = NB = 1),
*     the process grid must be square ( i.e. NPROW = NPCOL ) and
*     only lower triangular storage is supported.
*
*     Local variables
*     ===============
*
*     PDSYTTRD uses five local arrays:
*       WORK ( InV ) dimension ( NP, ANB+1): array V
*       WORK ( InH ) dimension ( NP, ANB+1): array H
*       WORK ( InVT ) dimension ( NQ, ANB+1): transpose of the array V
*       WORK ( InHT ) dimension ( NQ, ANB+1): transpose of the array H
*       WORK ( InVTT ) dimension ( NQ, 1): transpose of the array VT
*
*     Arrays V and H are replicated across all processor columns.
*     Arrays V^T and H^T are replicated across all processor rows.
*
*         WORK ( InVT ), or V^T, is stored as a tall skinny
*         array ( NQ x ANB-1 ) for efficiency.  Since only the lower
*         triangular portion of A is updated, Av is computed as:
*         tril(A) * v + v^T * tril(A,-1).  This is performed as
*         two local triangular matrix-vector multiplications (both in
*         MVR2) followed by a transpose and a sum across the columns.
*         In the local computation, WORK( InVT ) is used to compute
*         tril(A) * v and WORK( InV ) is used to compute
*         v^T * tril(A,-1)
*
*     The following variables are global indices into A:
*       INDEX:  The current global row and column number.
*       MAXINDEX:  The global row and column for the first row and
*       column in the trailing block of A.
*       LIIB, LIJB:  The first row, column in
*
*     The following variables point into the arrays A, V, H, V^T, H^T:
*       BINDEX  =INDEX-MININDEX: The column index in V, H, V^T, H^T.
*       LII:  local index I:  The local row number for row INDEX
*       LIJ:  local index J:  The local column number for column INDEX
*       LIIP1:  local index I+1:  The local row number for row INDEX+1
*       LIJP1:  local index J+1:  The local col number for col INDEX+1
*       LTLI: lower triangular local index I:  The local row for the
*         upper left entry in tril( A(INDEX, INDEX) )
*       LTLIP1: lower triangular local index I+1:  The local row for the
*         upper left entry in tril( A(INDEX+1, INDEX+1) )
*
*         Details:  The distinction between LII and LTLI (and between
*         LIIP1 and LTLIP1) is subtle.  Within the current processor
*         column (i.e. MYCOL .eq. CURCOL) they are the same.  However,
*         on some processors, A( LII, LIJ ) points to an element
*         above the diagonal, on these processors, LTLI = LII+1.
*
*     The following variables give the number of rows and/or columns
*     in various matrices:
*       NP:  The number of local rows in A( 1:N, 1:N )
*       NQ:  The number of local columns in A( 1:N, 1:N )
*       NPM0:  The number of local rows in A( INDEX:N, INDEX:N )
*       NQM0:  The number of local columns in A( INDEX:N, INDEX:N )
*       NPM1:  The number of local rows in A( INDEX+1:N, INDEX:N )
*       NQM1:  The number of local columns in A( INDEX+1:N, INDEX:N )
*       LTNM0:  The number of local rows & columns in
*         tril( A( INDEX:N, INDEX:N ) )
*       LTNM1:  The number of local rows & columns in
*         tril( A( INDEX+1:N, INDEX+1:N ) )
*         NOTE:  LTNM0 == LTNM1 on all processors except the diagonal
*         processors, i.e. those where MYCOL == MYROW.
*
*         Invariants:
*           NP = NPM0 + LII - 1
*           NQ = NQM0 + LIJ - 1
*           NP = NPM1 + LIIP1 - 1
*           NQ = NQM1 + LIJP1 - 1
*           NP = LTLI + LTNM0 - 1
*           NP = LTLIP1 + LTNM1 - 1
*
*       Temporary variables.  The following variables are used within
*       a few lines after they are set and do hold state from one loop
*       iteration to the next:
*
*     The matrix A:
*       The matrix A does not hold the same values that it would
*       in an unblocked code nor the values that it would hold in
*       in a blocked code.
*
*       The value of A is confusing.  It is easiest to state the
*       difference between trueA and A at the point that MVR2 is called,
*       so we will start there.
*
*       Let trueA be the value that A would
*       have at a given point in an unblocked code and A
*       be the value that A has in this code at the same point.
*
*       At the time of the call to MVR2,
*       trueA = A + V' * H + H' * V
*       where H = H( MAXINDEX:N, 1:BINDEX ) and
*       V = V( MAXINDEX:N, 1:BINDEX ).
*
*       At the bottom of the inner loop,
*       trueA = A +  V' * H + H' * V + v' * h + h' * v
*       where H = H( MAXINDEX:N, 1:BINDEX ) and
*       V = V( MAXINDEX:N, 1:BINDEX ) and
*       v = V( liip1:N, BINDEX+1 ) and
*       h = H( liip1:N, BINDEX+1 )
*
*       At the top of the loop, BINDEX gets incremented, hence:
*       trueA = A +  V' * H + H' * V + v' * h + h' * v
*       where H = H( MAXINDEX:N, 1:BINDEX-1 ) and
*       V = V( MAXINDEX:N, 1:BINDEX-1 ) and
*       v = V( liip1:N, BINDEX ) and
*       h = H( liip1:N, BINDEX )
*
*
*       A gets updated at the bottom of the outer loop
*       After this update, trueA = A + v' * h + h' * v
*       where v = V( liip1:N, BINDEX ) and
*       h = H( liip1:N, BINDEX ) and BINDEX = 0
*       Indeed, the previous loop invariant as stated above for the
*       top of the loop still holds, but with BINDEX = 0, H and V
*       are null matrices.
*
*       After the current column of A is updated,
*         trueA( INDEX, INDEX:N ) = A( INDEX, INDEX:N )
*       the rest of A is untouched.
*
*       After the current block column of A is updated,
*       trueA = A + V' * H + H' * V
*       where H = H( MAXINDEX:N, 1:BINDEX ) and
*       V = V( MAXINDEX:N, 1:BINDEX )
*
*       This brings us back to the point at which mvr2 is called.
*
*
*     Details of the parallelization:
*
*       We delay spreading v across to all processor columns (which
*       would naturally happen at the bottom of the loop) in order to
*       combine the spread of v( : , i-1 ) with the spread of h( : , i )
*
*       In order to compute h( :, i ), we must update A( :, i )
*       which means that the processor column owning A( :, i ) must
*       have: c, tau, v( i, i ) and h( i, i ).
*
*       The traditional
*       way of computing v (and the one used in pzlatrd.f and
*       zlatrd.f) is:
*         v = tau * v
*         c = v' * h
*         alpha = - tau * c / 2
*         v = v + alpha * h
*       However, the traditional way of computing v requires that tau
*       be broadcast to all processors in the current column (to compute
*       v = tau * v) and then a sum-to-all is required (to
*       compute v' * h ).  We use the following formula instead:
*         c = v' * h
*         v = tau * ( v - c * tau' * h / 2 )
*       The above formula allows tau to be spread down in the
*       same call to DGSUM2D which performs the sum-to-all of c.
*
*       The computation of v, which could be performed in any processor
*       column (or other procesor subsets), is performed in the
*       processor column that owns A( :, i+1 ) so that A( :, i+1 )
*       can be updated prior to spreading v across.
*
*       We keep the block column of A up-to-date to minimize the
*       work required in updating the current column of A.  Updating
*       the block column of A is reasonably load balanced whereas
*       updating the current column of A is not (only the current
*       processor column is involved).
*
*     In the following overview of the steps performed, M in the
*     margin indicates message traffic and C indicates O(n^2 nb/sqrt(p))
*     or more flops per processor.
*
*     Inner loop:
*       A( index:n, index ) -= ( v * ht(bindex) + h * vt( bindex) )
*M      h = house( A(index:n, index) )
*M      Spread v, h across
*M      vt = v^T; ht = h^T
*       A( index+1:n, index+1:maxindex ) -=
*         ( v * ht(index+1:maxindex) + h *vt(index+1:maxindex) )
*C      v = tril(A) * h; vt = ht * tril(A,-1)
*MorC   v = v - H*V*h - V*H*h
*M      v = v + vt^T
*M      c = v' * h
*       v = tau * ( v - c * tau' * h / 2 )
*C    A = A - H*V - V*H
*
*
*
*     =================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   Z_ONE, Z_NEGONE, Z_ZERO
      PARAMETER          ( Z_ONE = 1.0D0, Z_NEGONE = -1.0D0,
     $                   Z_ZERO = 0.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*
*
*     .. Local Scalars ..
*
*
      LOGICAL            BALANCED, INTERLEAVE, TWOGEMMS, UPPER
      INTEGER            ANB, BINDEX, CURCOL, CURROW, I, ICTXT, INDEX,
     $                   INDEXA, INDEXINH, INDEXINV, INH, INHB, INHT,
     $                   INHTB, INTMP, INV, INVB, INVT, INVTB, J, LDA,
     $                   LDV, LDZG, LII, LIIB, LIIP1, LIJ, LIJB, LIJP1,
     $                   LTLIP1, LTNM1, LWMIN, MAXINDEX, MININDEX,
     $                   MYCOL, MYFIRSTROW, MYROW, MYSETNUM, NBZG, NP,
     $                   NPB, NPCOL, NPM0, NPM1, NPROW, NPS, NPSET, NQ,
     $                   NQB, NQM1, NUMROWS, NXTCOL, NXTROW, PBMAX,
     $                   PBMIN, PBSIZE, PNB, ROWSPERPROC
      DOUBLE PRECISION   ALPHA, BETA, C, CONJTOPH, CONJTOPV, NORM,
     $                   ONEOVERBETA, SAFMAX, SAFMIN, TOPH, TOPNV,
     $                   TOPTAU, TOPV
*     ..
*     .. Local Arrays ..
*
*
*
*
      INTEGER            IDUM1( 1 ), IDUM2( 1 )
      DOUBLE PRECISION   CC( 3 ), DTMP( 5 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DCOMBNRM2, DGEBR2D,
     $                   DGEBS2D, DGEMM, DGEMV, DGERV2D, DGESD2D,
     $                   DGSUM2D, DLAMOV, DSCAL, DTRMVT, PCHK1MAT,
     $                   PDTREECOMB, PXERBLA
*     ..
*     .. External Functions ..
*
      LOGICAL            LSAME
      INTEGER            ICEIL, NUMROC, PJLAENV
      DOUBLE PRECISION   DNRM2, PDLAMCH
      EXTERNAL           LSAME, ICEIL, NUMROC, PJLAENV, DNRM2, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, ICHAR, MAX, MIN, MOD, SIGN, SQRT
*     ..
*
*
*     .. Executable Statements ..
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*
*
*     Further details
*     ===============
*
*     At the top of the loop, v and nh have been computed but not
*     spread across.  Hence, A is out-of-date even after the
*     rank 2k update.  Furthermore, we compute the next v before
*     nh is spread across.
*
*     I claim that if we used a sum-to-all on NV, by summing CC within
*     each column, that we could compute NV locally and could avoid
*     spreading V across.  Bruce claims that sum-to-all can be made
*     to cost no more than sum-to-one on the Paragon.  If that is
*     true, this would be a win.  But,
*     the BLACS sum-to-all is just a sum-to-one followed by a broadcast,
*     and hence the present scheme is better for now.
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      SAFMAX = SQRT( PDLAMCH( ICTXT, 'O' ) ) / N
      SAFMIN = SQRT( PDLAMCH( ICTXT, 'S' ) )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 600+CTXT_ )
      ELSE
*
*     Here we set execution options for PDSYTTRD
*
         PNB = PJLAENV( ICTXT, 2, 'PDSYTTRD', 'L', 0, 0, 0, 0 )
         ANB = PJLAENV( ICTXT, 3, 'PDSYTTRD', 'L', 0, 0, 0, 0 )
*
         INTERLEAVE = ( PJLAENV( ICTXT, 4, 'PDSYTTRD', 'L', 1, 0, 0,
     $                0 ).EQ.1 )
         TWOGEMMS = ( PJLAENV( ICTXT, 4, 'PDSYTTRD', 'L', 2, 0, 0,
     $              0 ).EQ.1 )
         BALANCED = ( PJLAENV( ICTXT, 4, 'PDSYTTRD', 'L', 3, 0, 0,
     $              0 ).EQ.1 )
*
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
*
*
         UPPER = LSAME( UPLO, 'U' )
         IF( INFO.EQ.0 .AND. DESCA( NB_ ).NE.1 )
     $      INFO = 600 + NB_
         IF( INFO.EQ.0 ) THEN
*
*
*           Here is the arithmetic:
*             Let maxnpq = max( np, nq, 2 * ANB )
*             LDV = 4 * max( np, nq ) + 2
*             LWMIN = 2 * ( ANB + 1 ) * LDV + MAX( np, 2 * ANB )
*             = 2 * ( ANB + 1 ) * ( 4 * NPS + 2 ) + NPS
*
*           This overestimates memory requirements when ANB > NP/2
*           Memory requirements are lower when interleave = .false.
*           Hence, we could have two sets of memory requirements,
*           one for interleave and one for
*
*
            NPS = MAX( NUMROC( N, 1, 0, 0, NPROW ), 2*ANB )
            LWMIN = 2*( ANB+1 )*( 4*NPS+2 ) + NPS
*
            WORK( 1 ) = DBLE( LWMIN )
            IF( .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( IA.NE.1 ) THEN
               INFO = -4
            ELSE IF( JA.NE.1 ) THEN
               INFO = -5
            ELSE IF( NPROW.NE.NPCOL ) THEN
               INFO = -( 600+CTXT_ )
            ELSE IF( DESCA( DTYPE_ ).NE.1 ) THEN
               INFO = -( 600+DTYPE_ )
            ELSE IF( DESCA( MB_ ).NE.1 ) THEN
               INFO = -( 600+MB_ )
            ELSE IF( DESCA( NB_ ).NE.1 ) THEN
               INFO = -( 600+NB_ )
            ELSE IF( DESCA( RSRC_ ).NE.0 ) THEN
               INFO = -( 600+RSRC_ )
            ELSE IF( DESCA( CSRC_ ).NE.0 ) THEN
               INFO = -( 600+CSRC_ )
            ELSE IF( LWORK.LT.LWMIN ) THEN
               INFO = -11
            END IF
         END IF
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 1
*
         CALL PCHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, 1, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDSYTTRD', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*
*
*     Reduce the lower triangle of sub( A )
      NP = NUMROC( N, 1, MYROW, 0, NPROW )
      NQ = NUMROC( N, 1, MYCOL, 0, NPCOL )
*
      NXTROW = 0
      NXTCOL = 0
*
      LIIP1 = 1
      LIJP1 = 1
      NPM1 = NP
      NQM1 = NQ
*
      LDA = DESCA( LLD_ )
      ICTXT = DESCA( CTXT_ )
*
*
*
*     Miscellaneous details:
*     Put tau, D and E in the right places
*     Check signs
*     Place all the arrays in WORK, control their placement
*     in  memory.
*
*
*
*     Loop invariants
*     A(LIIP1, LIJ) points to the first element of A(I+1,J)
*     NPM1,NQM1 = the number of rows, cols in A( LII+1:N,LIJ+1:N )
*     A(LII:N,LIJ:N) is one step out of date.
*     proc( CURROW, CURCOL ) owns A(LII,LIJ)
*     proc( NXTROW, CURCOL ) owns A(LIIP1,LIJ)
*
      INH = 1
*
      IF( INTERLEAVE ) THEN
*
*        H and V are interleaved to minimize memory movement
*        LDV has to be twice as large to accomodate interleaving.
*        In addition, LDV is doubled again to allow v, h and
*        toptau to be spreaad across and transposed in a
*        single communication operation with minimum memory
*        movement.
*
*        We could reduce LDV back to 2*MAX(NPM1,NQM1)
*        by increasing the memory movement required in
*        the spread and transpose of v, h and toptau.
*        However, since the non-interleaved path already
*        provides a mear minimum memory requirement option,
*        we did not provide this additional path.
*
         LDV = 4*( MAX( NPM1, NQM1 ) ) + 2
*
         INH = 1
*
         INV = INH + LDV / 2
         INVT = INH + ( ANB+1 )*LDV
*
         INHT = INVT + LDV / 2
         INTMP = INVT + LDV*( ANB+1 )
*
      ELSE
         LDV = MAX( NPM1, NQM1 )
*
         INHT = INH + LDV*( ANB+1 )
         INV = INHT + LDV*( ANB+1 )
*
*        The code works without this +1, but only because of a
*        coincidence.  Without the +1, WORK(INVT) gets trashed, but
*        WORK(INVT) is only used once and when it is used, it is
*        multiplied by WORK( INH ) which is zero.  Hence, the fact
*        that WORK(INVT) is trashed has no effect.
*
         INVT = INV + LDV*( ANB+1 ) + 1
         INTMP = INVT + LDV*( 2*ANB )
*
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDSYTTRD', -INFO )
         WORK( 1 ) = DBLE( LWMIN )
         RETURN
      END IF
*
*
*        The satisfies the loop invariant: trueA = A - V * HT - H * VT,
*        (where V, H, VT and HT all have BINDEX+1 rows/columns)
*        the first ANB times through the loop.
*
*
*
*     Setting either ( InH and InHT ) or InV to Z_ZERO
*     is adequate except in the face of NaNs.
*
*
      DO 10 I = 1, NP
         WORK( INH+I-1 ) = Z_ZERO
         WORK( INV+I-1 ) = Z_ZERO
   10 CONTINUE
      DO 20 I = 1, NQ
         WORK( INHT+I-1 ) = Z_ZERO
   20 CONTINUE
*
*
*
      TOPNV = Z_ZERO
*
      LTLIP1 = LIJP1
      LTNM1 = NPM1
      IF( MYCOL.GT.MYROW ) THEN
         LTLIP1 = LTLIP1 + 1
         LTNM1 = LTNM1 - 1
      END IF
*
*
      DO 210 MININDEX = 1, N - 1, ANB
*
*
         MAXINDEX = MIN( MININDEX+ANB-1, N )
         LIJB = NUMROC( MAXINDEX, 1, MYCOL, 0, NPCOL ) + 1
         LIIB = NUMROC( MAXINDEX, 1, MYROW, 0, NPROW ) + 1
*
         NQB = NQ - LIJB + 1
         NPB = NP - LIIB + 1
         INHTB = INHT + LIJB - 1
         INVTB = INVT + LIJB - 1
         INHB = INH + LIIB - 1
         INVB = INV + LIIB - 1
*
*
*
*
         DO 160 INDEX = MININDEX, MIN( MAXINDEX, N-1 )
*
            BINDEX = INDEX - MININDEX
*
            CURROW = NXTROW
            CURCOL = NXTCOL
*
            NXTROW = MOD( CURROW+1, NPROW )
            NXTCOL = MOD( CURCOL+1, NPCOL )
*
            LII = LIIP1
            LIJ = LIJP1
            NPM0 = NPM1
*
            IF( MYROW.EQ.CURROW ) THEN
               NPM1 = NPM1 - 1
               LIIP1 = LIIP1 + 1
            END IF
            IF( MYCOL.EQ.CURCOL ) THEN
               NQM1 = NQM1 - 1
               LIJP1 = LIJP1 + 1
               LTLIP1 = LTLIP1 + 1
               LTNM1 = LTNM1 - 1
            END IF
*
*
*
*
*     V = NV, VT = NVT, H = NH, HT = NHT
*
*
*     Update the current column of A
*
*
            IF( MYCOL.EQ.CURCOL ) THEN
*
               INDEXA = LII + ( LIJ-1 )*LDA
               INDEXINV = INV + LII - 1 + ( BINDEX-1 )*LDV
               INDEXINH = INH + LII - 1 + ( BINDEX-1 )*LDV
               CONJTOPH = WORK( INHT+LIJ-1+BINDEX*LDV )
               CONJTOPV = TOPNV
*
               IF( INDEX.GT.1 ) THEN
                  DO 30 I = 0, NPM0 - 1
*                  A( INDEXA+I ) = A( INDEXA+I )
                     A( INDEXA+I ) = A( INDEXA+I ) -
     $                               WORK( INDEXINV+LDV+I )*CONJTOPH -
     $                               WORK( INDEXINH+LDV+I )*CONJTOPV
   30             CONTINUE
               END IF
*
*
            END IF
*
*
            IF( MYCOL.EQ.CURCOL ) THEN
*
*     Compute the householder vector
*
               IF( MYROW.EQ.CURROW ) THEN
                  DTMP( 2 ) = A( LII+( LIJ-1 )*LDA )
               ELSE
                  DTMP( 2 ) = ZERO
               END IF
               IF( MYROW.EQ.NXTROW ) THEN
                  DTMP( 3 ) = A( LIIP1+( LIJ-1 )*LDA )
                  DTMP( 4 ) = ZERO
               ELSE
                  DTMP( 3 ) = ZERO
                  DTMP( 4 ) = ZERO
               END IF
*
               NORM = DNRM2( NPM1, A( LIIP1+( LIJ-1 )*LDA ), 1 )
               DTMP( 1 ) = NORM
*
*              IF DTMP(5) = 1.0, NORM is too large and might cause
*              overflow, hence PDTREECOMB must be called.  IF DTMP(5)
*              is zero on output, DTMP(1) can be trusted.
*
               DTMP( 5 ) = ZERO
               IF( DTMP( 1 ).GE.SAFMAX .OR. DTMP( 1 ).LT.SAFMIN ) THEN
                  DTMP( 5 ) = ONE
                  DTMP( 1 ) = ZERO
               END IF
*
               DTMP( 1 ) = DTMP( 1 )*DTMP( 1 )
               CALL DGSUM2D( ICTXT, 'C', ' ', 5, 1, DTMP, 5, -1,
     $                       CURCOL )
               IF( DTMP( 5 ).EQ.ZERO ) THEN
                  DTMP( 1 ) = SQRT( DTMP( 1 ) )
               ELSE
                  DTMP( 1 ) = NORM
                  CALL PDTREECOMB( ICTXT, 'C', 1, DTMP, -1, MYCOL,
     $                             DCOMBNRM2 )
               END IF
*
               NORM = DTMP( 1 )
*
               D( LIJ ) = DTMP( 2 )
               IF( MYROW.EQ.CURROW .AND. MYCOL.EQ.CURCOL ) THEN
                  A( LII+( LIJ-1 )*LDA ) = D( LIJ )
               END IF
*
*
               ALPHA = DTMP( 3 )
*
               NORM = SIGN( NORM, ALPHA )
*
               IF( NORM.EQ.ZERO ) THEN
                  TOPTAU = ZERO
               ELSE
                  BETA = NORM + ALPHA
                  TOPTAU = BETA / NORM
                  ONEOVERBETA = 1.0D0 / BETA
*
                  CALL DSCAL( NPM1, ONEOVERBETA,
     $                        A( LIIP1+( LIJ-1 )*LDA ), 1 )
               END IF
*
               IF( MYROW.EQ.NXTROW ) THEN
                  A( LIIP1+( LIJ-1 )*LDA ) = Z_ONE
               END IF
*
               TAU( LIJ ) = TOPTAU
               E( LIJ ) = -NORM
*
            END IF
*
*
*     Spread v, nh, toptau across
*
            DO 40 I = 0, NPM1 - 1
               WORK( INV+LIIP1-1+BINDEX*LDV+NPM1+I ) = A( LIIP1+I+
     $            ( LIJ-1 )*LDA )
   40       CONTINUE
*
            IF( MYCOL.EQ.CURCOL ) THEN
               WORK( INV+LIIP1-1+BINDEX*LDV+NPM1+NPM1 ) = TOPTAU
               CALL DGEBS2D( ICTXT, 'R', ' ', NPM1+NPM1+1, 1,
     $                       WORK( INV+LIIP1-1+BINDEX*LDV ),
     $                       NPM1+NPM1+1 )
            ELSE
               CALL DGEBR2D( ICTXT, 'R', ' ', NPM1+NPM1+1, 1,
     $                       WORK( INV+LIIP1-1+BINDEX*LDV ),
     $                       NPM1+NPM1+1, MYROW, CURCOL )
               TOPTAU = WORK( INV+LIIP1-1+BINDEX*LDV+NPM1+NPM1 )
            END IF
            DO 50 I = 0, NPM1 - 1
               WORK( INH+LIIP1-1+( BINDEX+1 )*LDV+I ) = WORK( INV+LIIP1-
     $            1+BINDEX*LDV+NPM1+I )
   50       CONTINUE
*
            IF( INDEX.LT.N ) THEN
               IF( MYROW.EQ.NXTROW .AND. MYCOL.EQ.CURCOL )
     $            A( LIIP1+( LIJ-1 )*LDA ) = E( LIJ )
            END IF
*
*     Transpose v, nh
*
*
            IF( MYROW.EQ.MYCOL ) THEN
               DO 60 I = 0, NPM1 + NPM1
                  WORK( INVT+LIJP1-1+BINDEX*LDV+I ) = WORK( INV+LIIP1-1+
     $               BINDEX*LDV+I )
   60          CONTINUE
            ELSE
               CALL DGESD2D( ICTXT, NPM1+NPM1, 1,
     $                       WORK( INV+LIIP1-1+BINDEX*LDV ), NPM1+NPM1,
     $                       MYCOL, MYROW )
               CALL DGERV2D( ICTXT, NQM1+NQM1, 1,
     $                       WORK( INVT+LIJP1-1+BINDEX*LDV ), NQM1+NQM1,
     $                       MYCOL, MYROW )
            END IF
*
            DO 70 I = 0, NQM1 - 1
               WORK( INHT+LIJP1-1+( BINDEX+1 )*LDV+I ) = WORK( INVT+
     $            LIJP1-1+BINDEX*LDV+NQM1+I )
   70       CONTINUE
*
*
*           Update the current block column of A
*
            IF( INDEX.GT.1 ) THEN
               DO 90 J = LIJP1, LIJB - 1
                  DO 80 I = 0, NPM1 - 1
*
                     A( LIIP1+I+( J-1 )*LDA ) = A( LIIP1+I+( J-1 )*LDA )
     $                   - WORK( INV+LIIP1-1+BINDEX*LDV+I )*
     $                  WORK( INHT+J-1+BINDEX*LDV ) -
     $                  WORK( INH+LIIP1-1+BINDEX*LDV+I )*
     $                  WORK( INVT+J-1+BINDEX*LDV )
   80             CONTINUE
   90          CONTINUE
            END IF
*
*
*
*     Compute NV = A * NHT; NVT = A * NH
*
*           These two lines are necessary because these elements
*           are not always involved in the calls to DTRMVT
*           for two reasons:
*           1)  On diagonal processors, the call to TRMVT
*               involves only LTNM1-1 elements
*           2)  On some processes, NQM1 < LTM1 or  LIIP1 < LTLIP1
*               and when the results are combined across all processes,
*               uninitialized values may be included.
            WORK( INV+LIIP1-1+( BINDEX+1 )*LDV ) = Z_ZERO
            WORK( INVT+LIJP1-1+( BINDEX+1 )*LDV+NQM1-1 ) = Z_ZERO
*
*
            IF( MYROW.EQ.MYCOL ) THEN
               IF( LTNM1.GT.1 ) THEN
                  CALL DTRMVT( 'L', LTNM1-1,
     $                         A( LTLIP1+1+( LIJP1-1 )*LDA ), LDA,
     $                         WORK( INVT+LIJP1-1+( BINDEX+1 )*LDV ), 1,
     $                         WORK( INH+LTLIP1+1-1+( BINDEX+1 )*LDV ),
     $                         1, WORK( INV+LTLIP1+1-1+( BINDEX+1 )*
     $                         LDV ), 1, WORK( INHT+LIJP1-1+( BINDEX+
     $                         1 )*LDV ), 1 )
               END IF
               DO 100 I = 1, LTNM1
                  WORK( INVT+LIJP1+I-1-1+( BINDEX+1 )*LDV )
     $               = WORK( INVT+LIJP1+I-1-1+( BINDEX+1 )*LDV ) +
     $               A( LTLIP1+I-1+( LIJP1+I-1-1 )*LDA )*
     $               WORK( INH+LTLIP1+I-1-1+( BINDEX+1 )*LDV )
  100          CONTINUE
            ELSE
               IF( LTNM1.GT.0 )
     $            CALL DTRMVT( 'L', LTNM1, A( LTLIP1+( LIJP1-1 )*LDA ),
     $                         LDA, WORK( INVT+LIJP1-1+( BINDEX+1 )*
     $                         LDV ), 1, WORK( INH+LTLIP1-1+( BINDEX+
     $                         1 )*LDV ), 1, WORK( INV+LTLIP1-1+
     $                         ( BINDEX+1 )*LDV ), 1,
     $                         WORK( INHT+LIJP1-1+( BINDEX+1 )*LDV ),
     $                         1 )
*
            END IF
*
*
*     We take advantage of the fact that:
*     A * sum( B ) = sum ( A * B ) for matrices A,B
*
*     trueA = A + V * HT + H * VT
*     hence:  (trueA)v = Av' + V * HT * v + H * VT * v
*     VT * v = sum_p_in_NPROW ( VTp * v )
*     H * VT * v = H * sum (VTp * v) = sum ( H * VTp * v )
*
*     v = v + V * HT * h + H * VT * h
*
*
*
*     tmp = HT * nh1
            DO 110 I = 1, 2*( BINDEX+1 )
               WORK( INTMP-1+I ) = 0
  110       CONTINUE
*
            IF( BALANCED ) THEN
               NPSET = NPROW
               MYSETNUM = MYROW
               ROWSPERPROC = ICEIL( NQB, NPSET )
               MYFIRSTROW = MIN( NQB+1, 1+ROWSPERPROC*MYSETNUM )
               NUMROWS = MIN( ROWSPERPROC, NQB-MYFIRSTROW+1 )
*
*
*     tmp = HT * v
*
               CALL DGEMV( 'C', NUMROWS, BINDEX+1, Z_ONE,
     $                     WORK( INHTB+MYFIRSTROW-1 ), LDV,
     $                     WORK( INHTB+MYFIRSTROW-1+( BINDEX+1 )*LDV ),
     $                     1, Z_ZERO, WORK( INTMP ), 1 )
*     tmp2 = VT * v
               CALL DGEMV( 'C', NUMROWS, BINDEX+1, Z_ONE,
     $                     WORK( INVTB+MYFIRSTROW-1 ), LDV,
     $                     WORK( INHTB+MYFIRSTROW-1+( BINDEX+1 )*LDV ),
     $                     1, Z_ZERO, WORK( INTMP+BINDEX+1 ), 1 )
*
*
               CALL DGSUM2D( ICTXT, 'C', ' ', 2*( BINDEX+1 ), 1,
     $                       WORK( INTMP ), 2*( BINDEX+1 ), -1, -1 )
            ELSE
*     tmp = HT * v
*
               CALL DGEMV( 'C', NQB, BINDEX+1, Z_ONE, WORK( INHTB ),
     $                     LDV, WORK( INHTB+( BINDEX+1 )*LDV ), 1,
     $                     Z_ZERO, WORK( INTMP ), 1 )
*     tmp2 = VT * v
               CALL DGEMV( 'C', NQB, BINDEX+1, Z_ONE, WORK( INVTB ),
     $                     LDV, WORK( INHTB+( BINDEX+1 )*LDV ), 1,
     $                     Z_ZERO, WORK( INTMP+BINDEX+1 ), 1 )
*
            END IF
*
*
*
            IF( BALANCED ) THEN
               MYSETNUM = MYCOL
*
               ROWSPERPROC = ICEIL( NPB, NPSET )
               MYFIRSTROW = MIN( NPB+1, 1+ROWSPERPROC*MYSETNUM )
               NUMROWS = MIN( ROWSPERPROC, NPB-MYFIRSTROW+1 )
*
               CALL DGSUM2D( ICTXT, 'R', ' ', 2*( BINDEX+1 ), 1,
     $                       WORK( INTMP ), 2*( BINDEX+1 ), -1, -1 )
*
*
*     v = v + V * tmp
               IF( INDEX.GT.1. ) THEN
                  CALL DGEMV( 'N', NUMROWS, BINDEX+1, Z_NEGONE,
     $                        WORK( INVB+MYFIRSTROW-1 ), LDV,
     $                        WORK( INTMP ), 1, Z_ONE,
     $                        WORK( INVB+MYFIRSTROW-1+( BINDEX+1 )*
     $                        LDV ), 1 )
*
*     v = v + H * tmp2
                  CALL DGEMV( 'N', NUMROWS, BINDEX+1, Z_NEGONE,
     $                        WORK( INHB+MYFIRSTROW-1 ), LDV,
     $                        WORK( INTMP+BINDEX+1 ), 1, Z_ONE,
     $                        WORK( INVB+MYFIRSTROW-1+( BINDEX+1 )*
     $                        LDV ), 1 )
               END IF
*
            ELSE
*     v = v + V * tmp
               CALL DGEMV( 'N', NPB, BINDEX+1, Z_NEGONE, WORK( INVB ),
     $                     LDV, WORK( INTMP ), 1, Z_ONE,
     $                     WORK( INVB+( BINDEX+1 )*LDV ), 1 )
*
*
*     v = v + H * tmp2
               CALL DGEMV( 'N', NPB, BINDEX+1, Z_NEGONE, WORK( INHB ),
     $                     LDV, WORK( INTMP+BINDEX+1 ), 1, Z_ONE,
     $                     WORK( INVB+( BINDEX+1 )*LDV ), 1 )
*
            END IF
*
*
*     Transpose NV and add it back into NVT
*
            IF( MYROW.EQ.MYCOL ) THEN
               DO 120 I = 0, NQM1 - 1
                  WORK( INTMP+I ) = WORK( INVT+LIJP1-1+( BINDEX+1 )*LDV+
     $                              I )
  120          CONTINUE
            ELSE
               CALL DGESD2D( ICTXT, NQM1, 1,
     $                       WORK( INVT+LIJP1-1+( BINDEX+1 )*LDV ),
     $                       NQM1, MYCOL, MYROW )
               CALL DGERV2D( ICTXT, NPM1, 1, WORK( INTMP ), NPM1, MYCOL,
     $                       MYROW )
*
            END IF
            DO 130 I = 0, NPM1 - 1
               WORK( INV+LIIP1-1+( BINDEX+1 )*LDV+I ) = WORK( INV+LIIP1-
     $            1+( BINDEX+1 )*LDV+I ) + WORK( INTMP+I )
  130       CONTINUE
*
*     Sum-to-one NV rowwise (within a row)
*
            CALL DGSUM2D( ICTXT, 'R', ' ', NPM1, 1,
     $                    WORK( INV+LIIP1-1+( BINDEX+1 )*LDV ), NPM1,
     $                    MYROW, NXTCOL )
*
*
*     Dot product c = NV * NH
*     Sum-to-all c within next processor column
*
*
            IF( MYCOL.EQ.NXTCOL ) THEN
               CC( 1 ) = Z_ZERO
               DO 140 I = 0, NPM1 - 1
                  CC( 1 ) = CC( 1 ) + WORK( INV+LIIP1-1+( BINDEX+1 )*
     $                      LDV+I )*WORK( INH+LIIP1-1+( BINDEX+1 )*LDV+
     $                      I )
  140          CONTINUE
               IF( MYROW.EQ.NXTROW ) THEN
                  CC( 2 ) = WORK( INV+LIIP1-1+( BINDEX+1 )*LDV )
                  CC( 3 ) = WORK( INH+LIIP1-1+( BINDEX+1 )*LDV )
               ELSE
                  CC( 2 ) = Z_ZERO
                  CC( 3 ) = Z_ZERO
               END IF
               CALL DGSUM2D( ICTXT, 'C', ' ', 3, 1, CC, 3, -1, NXTCOL )
*
               TOPV = CC( 2 )
               C = CC( 1 )
               TOPH = CC( 3 )
*
               TOPNV = TOPTAU*( TOPV-C*TOPTAU / 2*TOPH )
*
*
*     Compute V = Tau * (V - C * Tau' / 2 * H )
*
*
               DO 150 I = 0, NPM1 - 1
                  WORK( INV+LIIP1-1+( BINDEX+1 )*LDV+I ) = TOPTAU*
     $               ( WORK( INV+LIIP1-1+( BINDEX+1 )*LDV+I )-C*TOPTAU /
     $               2*WORK( INH+LIIP1-1+( BINDEX+1 )*LDV+I ) )
  150          CONTINUE
*
            END IF
*
*
  160    CONTINUE
*
*
*     Perform the rank2k update
*
         IF( MAXINDEX.LT.N ) THEN
*
            DO 170 I = 0, NPM1 - 1
               WORK( INTMP+I ) = WORK( INH+LIIP1-1+ANB*LDV+I )
  170       CONTINUE
*
*
*
            IF( .NOT.TWOGEMMS ) THEN
               IF( INTERLEAVE ) THEN
                  LDZG = LDV / 2
               ELSE
                  CALL DLAMOV( 'A', LTNM1, ANB, WORK( INHT+LIJP1-1 ),
     $                         LDV, WORK( INVT+LIJP1-1+ANB*LDV ), LDV )
*
                  CALL DLAMOV( 'A', LTNM1, ANB, WORK( INV+LTLIP1-1 ),
     $                         LDV, WORK( INH+LTLIP1-1+ANB*LDV ), LDV )
                  LDZG = LDV
               END IF
               NBZG = ANB*2
            ELSE
               LDZG = LDV
               NBZG = ANB
            END IF
*
*
            DO 180 PBMIN = 1, LTNM1, PNB
*
               PBSIZE = MIN( PNB, LTNM1-PBMIN+1 )
               PBMAX = MIN( LTNM1, PBMIN+PNB-1 )
               CALL DGEMM( 'N', 'C', PBSIZE, PBMAX, NBZG, Z_NEGONE,
     $                     WORK( INH+LTLIP1-1+PBMIN-1 ), LDZG,
     $                     WORK( INVT+LIJP1-1 ), LDZG, Z_ONE,
     $                     A( LTLIP1+PBMIN-1+( LIJP1-1 )*LDA ), LDA )
               IF( TWOGEMMS ) THEN
                  CALL DGEMM( 'N', 'C', PBSIZE, PBMAX, ANB, Z_NEGONE,
     $                        WORK( INV+LTLIP1-1+PBMIN-1 ), LDZG,
     $                        WORK( INHT+LIJP1-1 ), LDZG, Z_ONE,
     $                        A( LTLIP1+PBMIN-1+( LIJP1-1 )*LDA ), LDA )
               END IF
  180       CONTINUE
*
*
*
            DO 190 I = 0, NPM1 - 1
               WORK( INV+LIIP1-1+I ) = WORK( INV+LIIP1-1+ANB*LDV+I )
               WORK( INH+LIIP1-1+I ) = WORK( INTMP+I )
  190       CONTINUE
            DO 200 I = 0, NQM1 - 1
               WORK( INHT+LIJP1-1+I ) = WORK( INHT+LIJP1-1+ANB*LDV+I )
  200       CONTINUE
*
*
         END IF
*
*     End of the update A code
*
  210 CONTINUE
*
      IF( MYCOL.EQ.NXTCOL ) THEN
         IF( MYROW.EQ.NXTROW ) THEN
*
            D( NQ ) = A( NP+( NQ-1 )*LDA )
*
            CALL DGEBS2D( ICTXT, 'C', ' ', 1, 1, D( NQ ), 1 )
         ELSE
            CALL DGEBR2D( ICTXT, 'C', ' ', 1, 1, D( NQ ), 1, NXTROW,
     $                    NXTCOL )
         END IF
      END IF
*
*
*
*
      WORK( 1 ) = DBLE( LWMIN )
      RETURN
*
*     End of PDSYTTRD
*
*
      END
