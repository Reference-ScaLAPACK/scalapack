      SUBROUTINE PDLAUU2( UPLO, N, A, IA, JA, DESCA )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAUU2 computes the product U * U' or L' * L, where the triangular
*  factor U or L is stored in the upper or lower triangular part of
*  the matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1).
*
*  If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
*  overwriting the factor U in sub( A ).
*  If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
*  overwriting the factor L in sub( A ).
*
*  This is the unblocked form of the algorithm, calling Level 2 BLAS.
*  No communication is performed by this routine, the matrix to operate
*  on should be strictly local to one process.
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
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the triangular factor stored in the matrix
*          sub( A ) is upper or lower triangular:
*          = 'U':  Upper triangular,
*          = 'L':  Lower triangular.
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the order of the triangular factor U or L.  N >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, the local pieces of the triangular factor L or U.
*          On exit, if UPLO = 'U', the upper triangle of the distributed
*          matrix sub( A ) is overwritten with the upper triangle of the
*          product U * U'; if UPLO = 'L', the lower triangle of sub( A )
*          is overwritten with the lower triangle of the product L' * L.
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
      INTEGER            IACOL, IAROW, ICURR, IDIAG, IIA, IOFFA, JJA,
     $                   LDA, MYCOL, MYROW, NA, NPCOL, NPROW
      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, DGEMV, DSCAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT, LSAME
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get grid parameters and compute local indexes
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
*
      IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL ) THEN
*
         LDA = DESCA( LLD_ )
         IDIAG = IIA + ( JJA - 1 ) * LDA
         IOFFA = IDIAG
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Compute the product U * U'.
*
            DO 10 NA = N-1, 1, -1
               AII = A( IDIAG )
               ICURR = IDIAG + LDA
               A( IDIAG ) = AII*AII + DDOT( NA, A( ICURR ), LDA,
     $                                    A( ICURR ), LDA )
               CALL DGEMV( 'No transpose', N-NA-1, NA, ONE,
     $                     A( IOFFA+LDA ), LDA, A( ICURR ), LDA, AII,
     $                     A( IOFFA ), 1 )
               IDIAG = IDIAG + LDA + 1
               IOFFA = IOFFA + LDA
   10       CONTINUE
            AII = A( IDIAG )
            CALL DSCAL( N, AII, A( IOFFA ), 1 )
*
         ELSE
*
*           Compute the product L' * L.
*
            DO 20 NA = 1, N-1
               AII = A( IDIAG )
               ICURR = IDIAG + 1
               A(IDIAG) = AII*AII + DDOT( N-NA, A( ICURR ), 1,
     $                                    A( ICURR ), 1 )
               CALL DGEMV( 'Transpose', N-NA, NA-1, ONE, A( IOFFA+1 ),
     $                     LDA, A( ICURR ), 1, AII, A( IOFFA ), LDA )
               IDIAG = IDIAG + LDA + 1
               IOFFA = IOFFA + 1
   20       CONTINUE
            AII = A( IDIAG )
            CALL DSCAL( N, AII, A( IOFFA ), LDA )
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PDLAUU2
*
      END
