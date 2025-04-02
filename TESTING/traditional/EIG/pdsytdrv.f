      SUBROUTINE PDSYTDRV( UPLO, N, A, IA, JA, DESCA, D, E, TAU, WORK,
     $                     INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), D( * ), E( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDSYTDRV computes sub( A ) = A(IA:IA+N-1,JA:JA+N-1) from Q, the
*  symmetric tridiagonal matrix T (or D and E), and TAU, which were
*  computed by PDSYTRD:  sub( A ) := Q * T * Q'.
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
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix sub( A ) is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          This array contains the local pieces of sub( A ). On entry,
*          if UPLO='U', the diagonal and first superdiagonal of sub( A )
*          have the corresponding elements of the tridiagonal matrix T,
*          and the elements above the first superdiagonal, with the
*          array TAU, represent the orthogonal matrix Q as a product of
*          elementary reflectors, and the strictly lower triangular part
*          of sub( A ) is not referenced. If UPLO='L', the diagonal and
*          first subdiagonal of sub( A ) have the corresponding elements
*          of the tridiagonal matrix T, and the elements below the first
*          subdiagonal, with the array TAU, represent the orthogonal
*          matrix Q as a product of elementary reflectors, and the
*          strictly upper triangular part of sub( A ) is not referenced.
*          On exit, if UPLO = 'U', the upper triangular part of the
*          distributed symmetric matrix sub( A ) is recovered.
*          If UPLO='L', the lower triangular part of the distributed
*          symmetric matrix sub( A ) is recovered.
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
*  D       (local input) DOUBLE PRECISION array, dimension LOCc(JA+N-1)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i). D is tied to the distributed matrix A.
*
*  E       (local input) DOUBLE PRECISION array, dimension LOCc(JA+N-1)
*          if UPLO = 'U', LOCc(JA+N-2) otherwise. The off-diagonal
*          elements of the tridiagonal matrix T: E(i) = A(i,i+1) if
*          UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. E is tied to the
*          distributed matrix A.
*
*  TAU     (local input) DOUBLE PRECISION, array, dimension
*          LOCc(JA+N-1). This array contains the scalar factors TAU of
*          the elementary reflectors. TAU is tied to the distributed
*          matrix A.
*
*  WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK)
*          LWORK >= 2 * NB *( NB + NP )
*
*          where NB = MB_A = NB_A,
*          NP = NUMROC( N, NB, MYROW, IAROW, NPROW ),
*          IAROW = INDXG2P( IA, NB, MYROW, RSRC_A, NPROW ).
*
*          INDXG2P and NUMROC are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*  INFO    (global output) INTEGER
*          On exit, if INFO <> 0, a discrepancy has been found between
*          the diagonal and off-diagonal elements of A and the copies
*          contained in the arrays D and E.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   EIGHT, HALF, ONE, ZERO
      PARAMETER          ( EIGHT = 8.0D+0, HALF = 0.5D+0, ONE = 1.0D+0,
     $                     ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IACOL, IAROW, ICTXT, II, IPT, IPV, IPX,
     $                   IPY, J, JB, JJ, JL, K, MYCOL, MYROW, NB, NP,
     $                   NPCOL, NPROW
      DOUBLE PRECISION   ADDBND, D1, D2, E1, E2
*     ..
*     .. Local Arrays ..
      INTEGER            DESCD( DLEN_ ), DESCE( DLEN_ ), DESCV( DLEN_ ),
     $                   DESCT( DLEN_ )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           INDXG2P, LSAME, NUMROC, PDLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, INFOG2L, IGSUM2D,
     $                   PDELGET, PDGEMM, PDLACPY,
     $                   PDLARFT, PDLASET, PDSYMM,
     $                   PDSYR2K, PDTRMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD
*     ..
*     .. Executable statements ..
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      INFO = 0
      NB = DESCA( MB_ )
      UPPER = LSAME( UPLO, 'U' )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, II, JJ,
     $              IAROW, IACOL )
      NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
*
      IPT = 1
      IPV = NB * NB + IPT
      IPX = NB * NP + IPV
      IPY = NB * NP + IPX
*
      CALL DESCSET( DESCD, 1, JA+N-1, 1, DESCA( NB_ ), MYROW,
     $              DESCA( CSRC_ ), DESCA( CTXT_ ), 1 )
*
      ADDBND = EIGHT * PDLAMCH( ICTXT, 'eps' )
*
      IF( UPPER ) THEN
*
         CALL DESCSET( DESCE, 1, JA+N-1, 1, DESCA( NB_ ), MYROW,
     $                 DESCA( CSRC_ ), DESCA( CTXT_ ), 1 )
*
         DO 10 J = 0, N-1
            D1 = ZERO
            E1 = ZERO
            D2 = ZERO
            E2 = ZERO
            CALL PDELGET( ' ', ' ', D2, D, 1, JA+J, DESCD )
            CALL PDELGET( 'Columnwise', ' ', D1, A, IA+J, JA+J, DESCA )
            IF( J.LT.(N-1) ) THEN
               CALL PDELGET( ' ', ' ', E2, E, 1, JA+J+1, DESCE )
               CALL PDELGET( 'Columnwise', ' ', E1, A, IA+J, JA+J+1,
     $                       DESCA )
            END IF
*
            IF( ( ABS( D1 - D2 ).GT.( ABS( D2 ) * ADDBND ) ) .OR.
     $          ( ABS( E1 - E2 ).GT.( ABS( E2 ) * ADDBND ) ) )
     $         INFO = INFO + 1
   10    CONTINUE
*
*        Compute the upper triangle of sub( A ).
*
         CALL DESCSET( DESCV, N, NB, NB, NB, IAROW, IACOL, ICTXT,
     $                 MAX( 1, NP ) )
         CALL DESCSET( DESCT, NB, NB, NB, NB, IAROW, IACOL, ICTXT, NB )
*
         DO 20 K = 0, N-1, NB
            JB = MIN( NB, N-K )
            I = IA + K
            J = JA + K
*
*           Compute the lower triangular matrix T.
*
            CALL PDLARFT( 'Backward', 'Columnwise', K+JB-1, JB, A, IA,
     $                    J, DESCA, TAU, WORK( IPT ), WORK( IPV ) )
*
*           Copy Householder vectors into WORK( IPV ).
*
            CALL PDLACPY( 'All', K+JB-1, JB, A, IA, J, DESCA,
     $                    WORK( IPV ), 1, 1, DESCV )
*
            IF( K.GT.0 ) THEN
               CALL PDLASET( 'Lower', JB+1, JB, ZERO, ONE, WORK( IPV ),
     $                       K, 1, DESCV )
            ELSE
               CALL PDLASET( 'Lower', JB, JB-1, ZERO, ONE, WORK( IPV ),
     $                       1, 2, DESCV )
               CALL PDLASET( 'Ge', JB, 1, ZERO, ZERO, WORK( IPV ), 1,
     $                       1, DESCV )
            END IF
*
*           Zero out the strict upper triangular part of A.
*
            IF( K.GT.0 ) THEN
               CALL PDLASET( 'Ge', K-1, JB, ZERO, ZERO, A, IA, J,
     $                       DESCA )
               CALL PDLASET( 'Upper', JB-1, JB-1, ZERO, ZERO, A, I-1,
     $                       J+1, DESCA )
            ELSE IF( JB.GT.1 ) THEN
               CALL PDLASET( 'Upper', JB-2, JB-2, ZERO, ZERO, A, IA,
     $                       J+2, DESCA )
            END IF
*
*           (1) X := A * V * T'
*
            CALL PDSYMM( 'Left', 'Upper', K+JB, JB, ONE, A, IA, JA,
     $                   DESCA, WORK( IPV ), 1, 1, DESCV, ZERO,
     $                   WORK( IPX ), 1, 1, DESCV )
            CALL PDTRMM( 'Right', 'Lower', 'Transpose', 'Non-Unit',
     $                   K+JB, JB, ONE, WORK( IPT ), 1, 1, DESCT,
     $                   WORK( IPX ), 1, 1, DESCV )
*
*           (2) X := X - 1/2 * V * (T * V' * X)
*
            CALL PDGEMM( 'Transpose', 'No transpose', JB, JB, K+JB, ONE,
     $                   WORK( IPV ), 1, 1, DESCV, WORK( IPX ), 1, 1,
     $                   DESCV, ZERO, WORK( IPY ), 1, 1, DESCT )
            CALL PDTRMM( 'Left', 'Lower', 'No transpose', 'Non-Unit',
     $                   JB, JB, ONE, WORK( IPT ), 1, 1, DESCT,
     $                   WORK( IPY ), 1, 1, DESCT )
            CALL PDGEMM( 'No tranpose', 'No transpose', K+JB, JB, JB,
     $                   -HALF, WORK( IPV ), 1, 1, DESCV, WORK( IPY ),
     $                   1, 1, DESCT, ONE, WORK( IPX ), 1, 1, DESCV )
*
*           (3) A := A - X * V' - V * X'
*
            CALL PDSYR2K( 'Upper', 'No transpose', K+JB, JB, -ONE,
     $                    WORK( IPV ), 1, 1, DESCV, WORK( IPX ), 1, 1,
     $                    DESCV, ONE, A, IA, JA, DESCA )
*
            DESCV( CSRC_ ) = MOD( DESCV( CSRC_ ) + 1, NPCOL )
            DESCT( CSRC_ ) = MOD( DESCT( CSRC_ ) + 1, NPCOL )
*
   20    CONTINUE
*
      ELSE
*
         CALL DESCSET( DESCE, 1, JA+N-2, 1, DESCA( NB_ ), MYROW,
     $                 DESCA( CSRC_ ), DESCA( CTXT_ ), 1 )
*
         DO 30 J = 0, N-1
            D1 = ZERO
            E1 = ZERO
            D2 = ZERO
            E2 = ZERO
            CALL PDELGET( ' ', ' ', D2, D, 1, JA+J, DESCD )
            CALL PDELGET( 'Columnwise', ' ', D1, A, IA+J, JA+J, DESCA )
            IF( J.LT.(N-1) ) THEN
               CALL PDELGET( ' ', ' ', E2, E, 1, JA+J, DESCE )
               CALL PDELGET( 'Columnwise', ' ', E1, A, IA+J+1, JA+J,
     $                       DESCA )
            END IF
*
            IF( ( ABS( D1 - D2 ).GT.( ABS( D2 ) * ADDBND ) ) .OR.
     $          ( ABS( E1 - E2 ).GT.( ABS( E2 ) * ADDBND ) ) )
     $         INFO = INFO + 1
   30    CONTINUE
*
*        Compute the lower triangle of sub( A ).
*
         JL = MAX( ( ( JA+N-2 ) / NB ) * NB + 1, JA )
         IACOL = INDXG2P( JL, NB, MYCOL, DESCA( CSRC_ ), NPCOL )
         CALL DESCSET( DESCV, N, NB, NB, NB, IAROW, IACOL, ICTXT,
     $                 MAX( 1, NP ) )
         CALL DESCSET( DESCT, NB, NB, NB, NB, INDXG2P( IA+JL-JA+1, NB,
     $                 MYROW, DESCA( RSRC_ ), NPROW ), IACOL, ICTXT,
     $                 NB )
*
         DO 40 J = JL, JA, -NB
            K  = J - JA + 1
            I  = IA + K - 1
            JB = MIN( N-K+1, NB )
*
*           Compute upper triangular matrix T from TAU.
*
            CALL PDLARFT( 'Forward', 'Columnwise', N-K, JB, A, I+1, J,
     $                    DESCA, TAU, WORK( IPT ), WORK( IPV ) )
*
*           Copy Householder vectors into WORK( IPV ).
*
            CALL PDLACPY( 'Lower', N-K, JB, A, I+1, J, DESCA,
     $                    WORK( IPV ), K+1, 1, DESCV )
            CALL PDLASET( 'Upper', N-K, JB, ZERO, ONE, WORK( IPV ),
     $                    K+1, 1, DESCV )
            CALL PDLASET( 'Ge', 1, JB, ZERO, ZERO, WORK( IPV ), K, 1,
     $                    DESCV )
*
*           Zero out the strict lower triangular part of A.
*
            CALL PDLASET( 'Lower', N-K-1, JB, ZERO, ZERO, A, I+2, J,
     $                    DESCA )
*
*           (1) X := A * V * T'
*
            CALL PDSYMM( 'Left', 'Lower', N-K+1, JB, ONE, A, I, J,
     $                   DESCA, WORK( IPV ), K, 1, DESCV, ZERO,
     $                   WORK( IPX ), K, 1, DESCV )
            CALL PDTRMM( 'Right', 'Upper', 'Transpose', 'Non-Unit',
     $                   N-K+1, JB, ONE, WORK( IPT ), 1, 1, DESCT,
     $                   WORK( IPX ), K, 1, DESCV )
*
*           (2) X := X - 1/2 * V * (T * V' * X)
*
            CALL PDGEMM( 'Transpose', 'No transpose', JB, JB, N-K+1,
     $                   ONE, WORK( IPV ), K, 1, DESCV, WORK( IPX ),
     $                   K, 1, DESCV, ZERO, WORK( IPY ), 1, 1, DESCT )
            CALL PDTRMM( 'Left', 'Upper', 'No transpose', 'Non-Unit',
     $                   JB, JB, ONE, WORK( IPT ), 1, 1, DESCT,
     $                   WORK( IPY ), 1, 1, DESCT )
            CALL PDGEMM( 'No transpose', 'No transpose', N-K+1, JB, JB,
     $                   -HALF, WORK( IPV ), K, 1, DESCV, WORK( IPY ),
     $                   1, 1, DESCT, ONE, WORK( IPX ), K, 1, DESCV )
*
*           (3) A := A - X * V' - V * X'
*
            CALL PDSYR2K( 'Lower', 'No tranpose', N-K+1, JB, -ONE,
     $                    WORK( IPV ), K, 1, DESCV, WORK( IPX ), K, 1,
     $                    DESCV, ONE, A, I, J, DESCA )
*
            DESCV( CSRC_ ) = MOD( DESCV( CSRC_ ) + NPCOL - 1, NPCOL )
            DESCT( RSRC_ ) = MOD( DESCT( RSRC_ ) + NPROW - 1, NPROW )
            DESCT( CSRC_ ) = MOD( DESCT( CSRC_ ) + NPCOL - 1, NPCOL )
*
   40    CONTINUE
*
      END IF
*
      CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, 0 )
*
      RETURN
*
*     End of PDSYTDRV
*
      END
