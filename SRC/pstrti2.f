      SUBROUTINE PSTRTI2( UPLO, DIAG, N, A, IA, JA, DESCA, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            IA, INFO, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      REAL               A( * )
*     ..
*
*  Purpose
*  =======
*
*  PSTRTI2 computes the inverse of a real upper or lower triangular
*  block matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1). This matrix should be
*  contained in one and only one process memory space (local operation).
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
*          = 'U':  sub( A ) is upper triangular;
*          = 'L':  sub( A ) is lower triangular.
*
*  DIAG    (global input) CHARACTER*1
*          = 'N':  sub( A ) is non-unit triangular
*          = 'U':  sub( A ) is unit triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) REAL pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)),
*          this array contains the local pieces of the triangular matrix
*          sub( A ). If UPLO = 'U', the leading N-by-N upper triangular
*          part of the matrix sub( A ) contains the upper triangular
*          matrix, and the strictly lower triangular part of sub( A )
*          is not referenced.  If UPLO = 'L', the leading N-by-N lower
*          triangular part of the matrix sub( A ) contains the lower
*          triangular matrix, and the strictly upper triangular part
*          of sub( A ) is not referenced. If DIAG = 'U', the diagonal
*          elements of sub( A ) are also not referenced and are assumed
*          to be 1.  On exit, the (triangular) inverse of the original
*          matrix, in the same storage format.
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
*  INFO    (local output) INTEGER
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
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            IACOL, IAROW, ICTXT, ICURR, IDIAG, IIA, IOFFA,
     $                   JJA, LDA, MYCOL, MYROW, NA, NPCOL, NPROW
      REAL               AJJ
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GRIDINFO, CHK1MAT, INFOG2L,
     $                   PXERBLA, SSCAL, STRMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
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
         CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, INFO )
         UPPER = LSAME( UPLO, 'U' )
         NOUNIT = LSAME( DIAG, 'N' )
         IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
            INFO = -1
         ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
            INFO = -2
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PSTRTI2', -INFO )
         CALL BLACS_ABORT( ICTXT, 1 )
         RETURN
      END IF
*
*     Compute local indexes
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
*
      IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL ) THEN
*
         LDA = DESCA( LLD_ )
*
         IF( UPPER ) THEN
*
            IOFFA = IIA + ( JJA - 1 ) * LDA
            ICURR = IOFFA + LDA
*
            IF( NOUNIT ) THEN
*
*              Compute inverse of upper non-unit triangular matrix.
*
               A( IOFFA ) = ONE / A( IOFFA )
               IDIAG = ICURR + 1
               DO 10 NA = 1, N-1
                  A( IDIAG ) = ONE / A( IDIAG )
                  AJJ = -A( IDIAG )
*
*                 Compute elements 1:j-1 of j-th column.
*
                  CALL STRMV( 'Upper', 'No transpose', DIAG, NA,
     $                        A( IOFFA ), LDA, A( ICURR ), 1 )
                  CALL SSCAL( NA, AJJ, A( ICURR ), 1 )
                  IDIAG = IDIAG + LDA + 1
                  ICURR = ICURR + LDA
   10          CONTINUE
*
            ELSE
*
*              Compute inverse of upper unit triangular matrix.
*
               DO 20 NA = 1, N-1
*
*                 Compute elements 1:j-1 of j-th column.
*
                  CALL STRMV( 'Upper', 'No transpose', DIAG, NA,
     $                        A( IOFFA ), LDA, A( ICURR ), 1 )
                  CALL SSCAL( NA, -ONE, A( ICURR ), 1 )
                  ICURR = ICURR + LDA
   20          CONTINUE
*
            END IF
*
         ELSE
*
            ICURR = IIA + N - 1 + ( JJA + N - 2 ) * LDA
            IOFFA = ICURR - LDA
*
            IF( NOUNIT ) THEN
*
*              Compute inverse of lower non-unit triangular matrix.
*
               A( ICURR ) = ONE / A( ICURR )
               IDIAG = IOFFA - 1
               DO 30 NA = 1, N-1
                  A( IDIAG ) = ONE / A( IDIAG )
                  AJJ = -A( IDIAG )
*
*                 Compute elements j+1:n of j-th column.
*
                  CALL STRMV( 'Lower', 'No transpose', DIAG, NA,
     $                        A( ICURR ), LDA, A( IOFFA ), 1 )
                  CALL SSCAL( NA, AJJ, A( IOFFA ), 1 )
                  ICURR = IDIAG
                  IDIAG = IDIAG - LDA - 1
                  IOFFA = IDIAG + 1
   30          CONTINUE
*
            ELSE
*
*              Compute inverse of lower unit triangular matrix.
*
               DO 40 NA = 1, N-1
*
*                 Compute elements j+1:n of j-th column.
*
                  CALL STRMV( 'Lower', 'No transpose', DIAG, NA,
     $                     A( ICURR ), LDA, A( IOFFA ), 1 )
                  CALL SSCAL( NA, -ONE, A( IOFFA ), 1 )
                  ICURR = ICURR - LDA - 1
                  IOFFA = ICURR - LDA
   40          CONTINUE
*
            END IF
*
         END IF
*
      END IF
*
*     End of PSTRTI2
*
      END
