      SUBROUTINE PDLAGGE( M, N, D, A, IA, JA, DESCA, ISEED, ORDER, WORK,
     $                    LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, LWORK, M, N, ORDER
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), ISEED( 4 )
      DOUBLE PRECISION   A( * ), D( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAGGE generates a real symmetric matrix A, by pre- and post-
*  multiplying a real diagonal matrix D with a random orthogonal
*  matrices:  A = U*D*VT.
*
*  This is just a quick implementation which will be replaced in the
*  future. The random matrix A1(m,n) is generated and random left
*  orthogonal matrix U(m,m) is obtained  by running QR on A1:
*     A1(m,n) = U(m,m)*R,
*  where U(m,m) is a product of min(m,n) Householder rotations.
*  Afterwards the space of A1 is reused for a second random matrix
*  A2(m,n), which is used to obtain the right orthogonal matrix VT(n,n)
*  by running LQ on A2:
*     A2(m,n) = L*VT(n,n).
*  This requires vastly more computation than necessary, but not 
*  significantly more communication than is used in the rest of this 
*  routine, and hence is not that much slower than an efficient 
*  solution.
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
*  M       (global input) INTEGER
*           Number of rows of the matrix A. M >= 0.
*
*  N       (global input) INTEGER
*          Number  of columns of matrix A.  N >= 0.
*
*  D       (local input) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of the diagonal matrix D.
*
*  A       (local output) DOUBLE PRECISION array
*          Global dimension (M, N), local dimension (MP, NQ)
*
*  IA      (global input) INTEGER
*          The global row index of the submatrix of the distributed
*          matrix A to operate on.
*
*  JA      (global input) INTEGER
*          The global column index of the submatrix of the distributed
*          matrix A to operate on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_
*          The array descriptor for the distributed matrix A.
*
*  ISEED   (global input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd. On exit, the seed is updated and will remain identical
*          on all processes in the context.
*
*  ORDER   (global input) INTEGER
*          Number of reflectors in the matrix Q
*          At present, ORDER .NE. N is not supported
*
*  WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (local input) INTEGER dimension of WORK
*          LWORK >= MAX( QR_WORK, LQ_WORK )
*          QR_WORK = LDAA*MAX( 1, NQ ) + 200 + MAX( 1, DTAU1 ) +
*                    MAX( SIZEMQRLEFT, SIZEQRF)
*          LQ_WORK = LDAA*MAX( 1, NQ ) + 200 + MAX( 1, DTAU2) +
*                    MAX( SIZEMLQRIGHT, SIZEQRF )
*          Where:
*          LDAA = DESCA( LLD_ )
*          MB_A = DESCA( MB_ )
*          NB_A = DESCA( NB_ )
*          RSRC_A = DESCA( RSRC_ )
*          CSRC_A = DESCA( CSRC_ )
*          LCM = ILCM( NPROW, NPCOL )
*          LCMQ = LCM / NPCOL
*          IROFFA = MOD( IA-1, MB_A )
*          ICOFFA = MOD( JA-1, NB_A )
*          IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW )
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL )
*          MP = NUMROC( N+IROFFA, MB_A, MYROW, IAROW, NPROW )
*          NQ = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL )
*          DTAU1 = NUMROC( JA + SIZE- 1, NB_A, MYCOL, IACOL, NPROW )
*          DTAU2 = NUMROC( IA + SIZE- 1, MB_A, MYROW, IAROW, NPROW )
*          SIZEMQRLEFT = MAX( (MB_A*(MB_A-1))/2, ( MP + NQ ) * MB_A )
*               + ( MP + NB_A ) * NB_A
*          SIZEMLQRIGHT =  MAX( (MB_A*(MB_A-1))/2, (MP + NQ)*MB_A ) +
*                      MB_A * MB_A
*          SIZEQRF = NB_A*NP + MB_A*NQ + NB_A*NB_A
*
*  INFO    (local output) INTEGER
*
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
* ======================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            CSRC_A, DTAU1, DTAU2, I, IACOL, IAROW, ICOFFA,
     $                   IROFFA, LCM, LCMQ, LDAA, LQ_WORK, LWMIN, MB_A,
     $                   MP, MYCOL, MYROW, NB_A, NPCOL, NPROW, NQ,
     $                   PTR2AA, PTR2TAU, PTR2WORK, QR_WORK, RSRC_A,
     $                   SIZE, SIZELQF, SIZEMLQRIGHT, SIZEMQRLEFT,
     $                   SIZEQRF
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PDELSET, PDGELQF,
     $                   PDGEQRF, PDLASET, PDMATGEN, PDORMLQ, PDORMQR,
     $                   PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILCM, INDXG2P, NUMROC
      EXTERNAL           ILCM, INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*     This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*DLEN_*DTYPE_*M_*N_.LT.0 )RETURN
*
*     Initialize grid information.
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
*     Check LWORK.
*
      INFO = 0
      SIZE = MIN( M, N )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -607
      ELSE
         CALL CHK1MAT( M, 1, N, 2, IA, JA, DESCA, 8, INFO )
      END IF
*     Calculation of a minimum workspace.
      LDAA = DESCA( LLD_ )
      MB_A = DESCA( MB_ )
      NB_A = DESCA( NB_ )
      RSRC_A = DESCA( RSRC_ )
      CSRC_A = DESCA( CSRC_ )
      LCM = ILCM( NPROW, NPCOL )
      LCMQ = LCM / NPCOL
      IROFFA = MOD( IA-1, MB_A )
      ICOFFA = MOD( JA-1, NB_A )
      IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW )
      IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL )
      DTAU1 = NUMROC( JA+SIZE-1, NB_A, MYCOL, IACOL, NPCOL )
      DTAU2 = NUMROC( IA+SIZE-1, MB_A, MYROW, IAROW, NPROW )
      MP = NUMROC( M+IROFFA, MB_A, MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL )
*
      SIZEMQRLEFT = MAX( ( MB_A*( MB_A-1 ) ) / 2, ( MP+NQ )*MB_A ) +
     $              ( MP+NB_A )*NB_A
      SIZEMLQRIGHT = MAX( ( MB_A*( MB_A-1 ) ) / 2, ( MP+NQ )*MB_A ) +
     $               MB_A*MB_A
      SIZEQRF = NB_A*MP + MB_A*NQ + NB_A*NB_A + 100
      SIZELQF = NB_A*( MP+NQ+NB_A ) + 100
*
      QR_WORK = LDAA*MAX( 1, NQ ) + 200 + MAX( 1, DTAU1 ) +
     $          MAX( SIZEMQRLEFT, SIZEQRF )
      LQ_WORK = LDAA*MAX( 1, NQ ) + 200 + MAX( 1, DTAU2 ) +
     $          MAX( SIZEMLQRIGHT, SIZELQF )
      LWMIN = MAX( QR_WORK, LQ_WORK )
      WORK( 1 ) = LWMIN
      IF( LWORK.EQ.-1 )
     $   GO TO 20
*
*     Test the input arguments.
*
      IF( INFO.EQ.0 ) THEN
         IF( SIZE.NE.ORDER ) THEN
            INFO = -9
         ELSE IF( LWORK.LT.LWMIN ) THEN
            INFO = -11
         END IF
      END IF
      IF( INFO.LT.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PDLAGGE', -INFO )
         RETURN
      END IF
*
*     Build a diagonal matrix A with the eigenvalues specified in D.
*
      CALL PDLASET( 'Full', M, N, ZERO, ZERO, A, IA, JA, DESCA )
      DO 10 I = 1, SIZE
         CALL PDELSET( A, I, I, DESCA, D( I ) )
   10 CONTINUE
*
*     Local dimension of array TAU in tis case is LOCc(JA+MIN(M,N)-1).
*
      PTR2AA = 2
      PTR2TAU = PTR2AA + LDAA*MAX( 1, NQ ) + 100
      PTR2WORK = PTR2TAU + MAX( 1, DTAU1 ) + 100
*
      CALL PDLASET( 'All', M, N, ZERO, ZERO, WORK( PTR2AA ), IA, JA,
     $              DESCA )
*
*     Build a random matrix AA1.
*
      CALL PDMATGEN( DESCA( CTXT_ ), 'N', 'N', M, N, DESCA( MB_ ),
     $               DESCA( NB_ ), WORK( PTR2AA ), DESCA( LLD_ ),
     $               DESCA( RSRC_ ), DESCA( CSRC_ ), ISEED( 1 ), 0, MP,
     $               0, NQ, MYROW, MYCOL, NPROW, NPCOL )
*
*     Produce QR decomposition AA1 -> U*R.
*
      CALL PDGEQRF( M, N, WORK( PTR2AA ), IA, JA, DESCA,
     $              WORK( PTR2TAU ), WORK( PTR2WORK ), SIZEQRF, INFO )
*
*     A = U*A.
*
      CALL PDORMQR( 'L', 'N', M, N, SIZE, WORK( PTR2AA ), IA, JA, DESCA,
     $              WORK( PTR2TAU ), A, IA, JA, DESCA, WORK( PTR2WORK ),
     $              SIZEMQRLEFT, INFO )
*
*     Reinitialize pointer to WORK array. Dimension of array TAU in
*     this case is  LOCr(IA+MIN(M,N)-1).
*
      PTR2WORK = PTR2TAU + MAX( 1, DTAU2 ) + 100
*
*     Use the same workspace to generate a random matrix AA2.
*
      CALL PDMATGEN( DESCA( CTXT_ ), 'N', 'N', M, N, DESCA( MB_ ),
     $               DESCA( NB_ ), WORK( PTR2AA ), DESCA( LLD_ ),
     $               DESCA( RSRC_ ), DESCA( CSRC_ ), ISEED( 2 ), 0, MP,
     $               0, NQ, MYROW, MYCOL, NPROW, NPCOL )
*
*     Produce LQ decomposition of random matrix  AA2 -> L*VT.
*
      CALL PDGELQF( M, N, WORK( PTR2AA ), IA, JA, DESCA,
     $              WORK( PTR2TAU ), WORK( PTR2WORK ), SIZELQF, INFO )
*
*     Calculate A = A*VT.
*
      CALL PDORMLQ( 'R', 'N', M, N, SIZE, WORK( PTR2AA ), IA, JA, DESCA,
     $              WORK( PTR2TAU ), A, IA, JA, DESCA, WORK( PTR2WORK ),
     $              SIZEMLQRIGHT, INFO )
*
*     End of PDLAGGE
*
   20 CONTINUE
      RETURN
      END
