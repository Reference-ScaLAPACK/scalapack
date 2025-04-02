      SUBROUTINE PDPOTRRV( UPLO, N, A, IA, JA, DESCA, WORK )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 28, 2001
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ),  WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDPOTRRV recomputes sub( A ) = A(IA:IA+N-1,JA:JA+N-1) from L or U
*  computed by PDPOTRF. The routine performs the Cholesky factorization
*  in reverse.
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
*          symmetric distributed matrix sub( A ) is stored:
*          stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, the local pieces of the factors L or U of the
*          distributed matrix sub( A ) from the Cholesky factorization.
*          On exit, the original distributed matrix sub( A ) is
*          restored.
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
*  WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK)
*          LWORK >= MB_A*NB_A.
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
      LOGICAL            UPPER
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            IACOL, IAROW, ICTXT, IL, J, JB, JL, JN, MYCOL,
     $                   MYROW, NPCOL, NPROW
*     .. Local Arrays ..
      INTEGER            DESCW( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, PDLACPY, PDLASET,
     $                   PDSYRK, PDTRMM, PB_TOPGET, PB_TOPSET
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, INDXG2P
      EXTERNAL           ICEIL, INDXG2P, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      UPPER = LSAME( UPLO, 'U' )
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      JL = MAX( ( ( JA+N-2 ) / DESCA( NB_ ) ) * DESCA( NB_ ) + 1, JA )
      IL = MAX( ( ( IA+N-2 ) / DESCA( MB_ ) ) * DESCA( MB_ ) + 1, IA )
      IAROW = INDXG2P( IL, DESCA( MB_ ), MYROW, DESCA( RSRC_ ), NPROW )
      IACOL = INDXG2P( JL, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ), NPCOL )
*
*     Define array descriptor for working array WORK
*
      CALL DESCSET( DESCW, DESCA( MB_ ), DESCA( NB_ ), DESCA( MB_ ),
     $              DESCA( NB_ ), IAROW, IACOL, ICTXT, DESCA( MB_ ) )
*
      IF ( UPPER ) THEN
*
*       Compute A from the Cholesky factor U : A = U'*U.
*
        CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ' ' )
        CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', 'S-ring' )
*
        DO 10 J = JL, JN+1, -DESCA( NB_ )
*
           JB = MIN( JA+N-J, DESCA( NB_ ) )
*
*          Update the trailing matrix, A = A + U'*U
*
           CALL PDSYRK( 'Upper', 'Transpose', JA+N-J-JB, JB, ONE, A, IL,
     $                  J+JB, DESCA, ONE, A, IL+JB, J+JB, DESCA )
*
*          Copy current diagonal block of A into workspace
*
           CALL PDLACPY( 'All', JB, JB, A, IL, J, DESCA, WORK, 1, 1,
     $                   DESCW )
*
*          Zero strict lower triangular part of diagonal block, to make
*          it U1.
*
           CALL PDLASET( 'Lower', JB-1, JB, ZERO, ZERO, A, IL+1, J,
     $                   DESCA )
*
*          Update the row panel U with the triangular matrix
*
           CALL PDTRMM( 'Left', 'Upper', 'Transpose', 'Non-Unit', JB,
     $                  N-J+JA, ONE, WORK, 1, 1, DESCW, A, IL, J,
     $                  DESCA )
*
*          Restore the strict lower triangular part of diagonal block.
*
           CALL PDLACPY( 'Lower', JB-1, JB, WORK, 2, 1, DESCW, A,
     $                    IL+1, J, DESCA )
*
           IL = IL - DESCA( MB_ )
           DESCW( RSRC_ ) = MOD( DESCW( RSRC_ ) + NPROW - 1, NPROW )
           DESCW( CSRC_ ) = MOD( DESCW( CSRC_ ) + NPCOL - 1, NPCOL )
*
   10   CONTINUE
*
*       Handle first block separately
*
        JB = MIN( JN-JA+1, DESCA( NB_ ) )
*
*       Update the trailing matrix, A = A + U'*U
*
        CALL PDSYRK( 'Upper', 'Transpose', N-JB, JB, ONE, A, IA, JA+JB,
     $               DESCA, ONE, A, IA+JB, JA+JB, DESCA )
*
*       Copy current diagonal block of A into workspace
*
        CALL PDLACPY( 'All', JB, JB, A, IA, JA, DESCA, WORK, 1, 1,
     $                DESCW )
*
*       Zero strict lower triangular part of diagonal block, to make
*       it U1.
*
        CALL PDLASET( 'Lower', JB-1, JB, ZERO, ZERO, A, IA+1, JA,
     $                DESCA )
*
*       Update the row panel U with the triangular matrix
*
        CALL PDTRMM( 'Left', 'Upper', 'Transpose', 'Non-Unit', JB,
     $               N, ONE, WORK, 1, 1, DESCW, A, IA, JA, DESCA )
*
*       Restore the strict lower triangular part of diagonal block.
*
        CALL PDLACPY( 'Lower', JB-1, JB, WORK, 2, 1, DESCW, A, IA+1,
     $                JA, DESCA )
*
      ELSE
*
*       Compute A from the Cholesky factor L : A = L*L'.
*
        CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'S-ring' )
        CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
*
        DO 20 J = JL, JN+1, -DESCA( NB_ )
*
           JB = MIN( JA+N-J, DESCA( NB_ ) )
*
*          Update the trailing matrix, A = A + L*L'
*
           CALL PDSYRK( 'Lower', 'No transpose', IA+N-IL-JB, JB, ONE, A,
     $                  IL+JB, J, DESCA, ONE, A, IL+JB, J+JB, DESCA )
*
*          Copy current diagonal block of A into workspace
*
           CALL PDLACPY( 'All', JB, JB, A, IL, J, DESCA, WORK, 1, 1,
     $                   DESCW )
*
*          Zero strict upper triangular part of diagonal block, to make
*          it L1.
*
           CALL PDLASET( 'Upper', JB, JB-1, ZERO, ZERO, A, IL, J+1,
     $                   DESCA )
*
*          Update the column panel L with the triangular matrix
*
           CALL PDTRMM( 'Right', 'Lower', 'Transpose', 'Non-Unit',
     $                   IA+N-IL, JB, ONE, WORK, 1, 1, DESCW, A, IL,
     $                   J, DESCA )
*
*          Restore the strict upper triangular part of diagonal block.
*
           CALL PDLACPY( 'Upper', JB, JB-1, WORK, 1, 2, DESCW, A,
     $                   IL, J+1, DESCA )
*
           IL = IL - DESCA( MB_ )
           DESCW( RSRC_ ) = MOD( DESCW( RSRC_ ) + NPROW - 1, NPROW )
           DESCW( CSRC_ ) = MOD( DESCW( CSRC_ ) + NPCOL - 1, NPCOL )
*
   20   CONTINUE
*
*       Handle first block separately
*
        JB = MIN( JN-JA+1, DESCA( NB_ ) )
*
*       Update the trailing matrix, A = A + L*L'
*
        CALL PDSYRK( 'Lower', 'No transpose', N-JB, JB, ONE, A,
     $               IA+JB, JA, DESCA, ONE, A, IA+JB, JA+JB, DESCA )
*
*       Copy current diagonal block of A into workspace
*
        CALL PDLACPY( 'All', JB, JB, A, IA, JA, DESCA, WORK, 1, 1,
     $                DESCW )
*
*       Zero strict upper triangular part of diagonal block, to make
*       it L1.
*
        CALL PDLASET( 'Upper', JB, JB-1, ZERO, ZERO, A, IA, JA+1,
     $                DESCA )
*
*       Update the column panel L with the triangular matrix
*
        CALL PDTRMM( 'Right', 'Lower', 'Transpose', 'Non-Unit', N, JB,
     $               ONE, WORK, 1, 1, DESCW, A, IA, JA, DESCA )
*
*       Restore the strict upper triangular part of diagonal block.
*
        CALL PDLACPY( 'Upper', JB, JB-1, WORK, 1, 2, DESCW, A, IA,
     $                JA+1, DESCA )
*
      END IF
*
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      RETURN
*
*     End of PDPOTRRV
*
      END
