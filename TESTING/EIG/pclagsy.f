*
*
      SUBROUTINE PCLAGHE( N, K, D, A, IA, JA, DESCA, ISEED, ORDER, WORK,
     $                    LWORK, INFO )
*
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, K, LWORK, N, ORDER
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), ISEED( 4 )
      REAL               D( * )
      COMPLEX            A( * ), WORK( * )
*     ..
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
*  Purpose
*  =======
*
*  PCLAGHE generates a real Hermitian matrix A, by pre- and post-
*  multiplying a real diagonal matrix D with a random orthogonal matrix:
*  A = U*D*U'.
*
*  This is just a quick implementation which will be replaced in the
*  future.  The random orthogonal matrix is computed by creating a
*  random matrix and running QR on it.  This requires vastly more
*  computation than necessary, but not significantly more communication
*  than is used in the rest of this rouinte, and hence is not that much
*  slower than an efficient solution.
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The size of the matrix A.  N >= 0.
*
*  K       (global input) INTEGER
*          The number of nonzero subdiagonals within the band of A.
*          0 <= K <= N-1.
*          ### K must be 0 or N-1, 0 < K < N-1 is not supported yet.
*
*  D       (global input) COMPLEX array, dimension (N)
*          The diagonal elements of the diagonal matrix D.
*
*  A       (local output) COMPLEX array
*          Global dimension (N, N), local dimension (NP, NQ)
*          The generated n by n Hermitian matrix A (the full matrix is
*          stored).
*
*  IA      (global input) INTEGER
*          A's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  ISEED   (global input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated and will remain identical on
*          all processes in the context.
*
*  ORDER   (global input) INTEGER
*          Number of reflectors in the matrix Q
*          At present, ORDER .NE. N is not supported
*
*  WORK    (local workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (local input) INTEGER dimension of WORK
*        LWORK >= SIZETMS as returned by PCLASIZESEP
*
*
*  INFO    (local output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      COMPLEX            CZERO
      PARAMETER          ( CZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            CSRC_A, I, IACOL, IAROW, ICOFFA, II, IIROW,
     $                   INDAA, INDTAU, INDWORK, IPOSTPAD, IPREPAD,
     $                   IROFFA, ISIZEHEEVX, ISIZESUBTST, ISIZETST,
     $                   JJCOL, LDAA, LII, LIII, LJJ, LJJJ, LWMIN, MAXI,
     $                   MB_A, MYCOL, MYROW, NB_A, NP, NPCOL, NPROW, NQ,
     $                   RSIZECHK, RSIZEHEEVX, RSIZEQTQ, RSIZESUBTST,
     $                   RSIZETST, RSRC_A, SIZEHEEVX, SIZEMQRLEFT,
     $                   SIZEMQRRIGHT, SIZEQRF, SIZESUBTST, SIZETMS,
     $                   SIZETST, SIZEHEEVD, RSIZEHEEVD, ISIZEHEEVD
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, CLASET, PCGEQRF,
     $                   PCLASIZESEP, PCMATGEN, PCUNMQR, PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
*
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Initialize grid information
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
*     Check LWORK
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 700+CTXT_ )
      ELSE
         CALL CHK1MAT( N, 1, N, 1, IA, JA, DESCA, 7, INFO )
      END IF
*
      LDAA = DESCA( LLD_ )
      MB_A = DESCA( MB_ )
      NB_A = DESCA( NB_ )
      RSRC_A = DESCA( RSRC_ )
      CSRC_A = DESCA( CSRC_ )
      IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW )
      IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL )
      IROFFA = MOD( IA-1, MB_A )
      ICOFFA = MOD( JA-1, NB_A )
      NP = NUMROC( N+IROFFA, MB_A, MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL )
      IPREPAD = 0
      IPOSTPAD = 0
      CALL PCLASIZESEP( DESCA, IPREPAD, IPOSTPAD, SIZEMQRLEFT,
     $                  SIZEMQRRIGHT, SIZEQRF, SIZETMS, RSIZEQTQ,
     $                  RSIZECHK, SIZEHEEVX, RSIZEHEEVX, ISIZEHEEVX,
     $                  SIZEHEEVD, RSIZEHEEVD, ISIZEHEEVD,
     $                  SIZESUBTST, RSIZESUBTST, ISIZESUBTST, SIZETST,
     $                  RSIZETST, ISIZETST )
      LWMIN = SIZETMS
*
*     Test the input arguments
*
      IF( INFO.EQ.0 ) THEN
         IF( K.LT.0 .OR. K.GT.N-1 ) THEN
            INFO = -2
         ELSE IF( N.NE.ORDER ) THEN
            INFO = -9
         ELSE IF( LWORK.LT.LWMIN ) THEN
            INFO = -11
         END IF
      END IF
      IF( INFO.LT.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PCLAGHE', -INFO )
         RETURN
      END IF
*
      INDAA = 1
      INDTAU = INDAA + LDAA*MAX( 1, NQ )
      INDWORK = INDTAU + MAX( 1, NQ )
*
      IF( K.NE.0 ) THEN
         CALL CLASET( 'A', LDAA, NQ, CZERO, CZERO, WORK( INDAA ), LDAA )
*
*
*        Build a random matrix
*
*
         CALL PCMATGEN( DESCA( CTXT_ ), 'N', 'N', N, ORDER,
     $                  DESCA( MB_ ), DESCA( NB_ ), WORK( INDAA ),
     $                  DESCA( LLD_ ), DESCA( RSRC_ ), DESCA( CSRC_ ),
     $                  ISEED( 1 ), 0, NP, 0, NQ, MYROW, MYCOL, NPROW,
     $                  NPCOL )
         CALL PCGEQRF( N, ORDER, WORK( INDAA ), IA, JA, DESCA,
     $                 WORK( INDTAU ), WORK( INDWORK ), SIZEQRF, INFO )
*
      END IF
*
*     Build a diagonal matrix A with the eigenvalues specified in D
*
      CALL CLASET( 'A', NP, NQ, CZERO, CZERO, A, DESCA( LLD_ ) )
*
      IIROW = 0
      JJCOL = 0
      LII = 1
      LJJ = 1
*
      DO 20 II = 1, N, DESCA( MB_ )
         MAXI = MIN( N, II+DESCA( MB_ )-1 )
         IF( ( MYROW.EQ.IIROW ) .AND. ( MYCOL.EQ.JJCOL ) ) THEN
            LIII = LII
            LJJJ = LJJ
            DO 10 I = II, MAXI
               A( LIII+( LJJJ-1 )*DESCA( LLD_ ) ) = D( I )
               LIII = LIII + 1
               LJJJ = LJJJ + 1
   10       CONTINUE
         END IF
         IF( MYROW.EQ.IIROW )
     $      LII = LII + DESCA( MB_ )
         IF( MYCOL.EQ.JJCOL )
     $      LJJ = LJJ + DESCA( MB_ )
         IIROW = MOD( IIROW+1, NPROW )
         JJCOL = MOD( JJCOL+1, NPCOL )
   20 CONTINUE
*
*     A = Q * A
*
      IF( K.NE.0 ) THEN
*
         CALL PCUNMQR( 'L', 'Conjugate transpose', N, N, ORDER,
     $                 WORK( INDAA ), IA, JA, DESCA, WORK( INDTAU ), A,
     $                 IA, JA, DESCA, WORK( INDWORK ), SIZEMQRLEFT,
     $                 INFO )
*
*
*        A = A * Q'
*
*
         CALL PCUNMQR( 'R', 'N', N, N, ORDER, WORK( INDAA ), IA, JA,
     $                 DESCA, WORK( INDTAU ), A, IA, JA, DESCA,
     $                 WORK( INDWORK ), SIZEMQRRIGHT, INFO )
*
      END IF
*
*     End of PCLAGHE
*
      END
