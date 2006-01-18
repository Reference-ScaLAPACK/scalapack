      SUBROUTINE PCPOTF2( UPLO, N, A, IA, JA, DESCA, INFO )
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
      COMPLEX            A( * )
*     ..
*
*  Purpose
*  =======
*
*  PCPOTF2 computes the Cholesky factorization of a complex hermitian
*  positive definite distributed matrix sub( A )=A(IA:IA+N-1,JA:JA+N-1).
*
*  The factorization has the form
*
*            sub( A ) = U' * U ,  if UPLO = 'U', or
*
*            sub( A ) = L  * L',  if UPLO = 'L',
*
*  where U is an upper triangular matrix and L is lower triangular.
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
*  This routine requires N <= NB_A-MOD(JA-1, NB_A) and square block
*  decomposition ( MB_A = NB_A ).
*
*  Arguments
*  =========
*
*  UPLO    (global input) CHARACTER
*          = 'U':  Upper triangle of sub( A ) is stored;
*          = 'L':  Lower triangle of sub( A ) is stored.
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) COMPLEX pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          N-by-N symmetric distributed matrix sub( A ) to be factored.
*          If UPLO = 'U', the leading N-by-N upper triangular part of
*          sub( A ) contains the upper triangular part of the matrix,
*          and its strictly lower triangular part is not referenced.
*          If UPLO = 'L', the leading N-by-N lower triangular part of
*          sub( A ) contains the lower triangular part of the distribu-
*          ted matrix, and its strictly upper triangular part is not
*          referenced.  On exit, if UPLO = 'U', the upper triangular
*          part of the distributed matrix contains the Cholesky factor
*          U, if UPLO = 'L', the lower triangular part of the distribu-
*          ted matrix contains the Cholesky factor L.
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
*          > 0:  If INFO = K, the leading minor of order K,
*                A(IA:IA+K-1,JA:JA+K-1) is not positive definite, and
*                the factorization could not be completed.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      CHARACTER          COLBTOP, ROWBTOP
      INTEGER            IACOL, IAROW, ICOFF, ICTXT, ICURR, IDIAG, IIA,
     $                   IOFFA, IROFF, J, JJA, LDA, MYCOL, MYROW,
     $                   NPCOL, NPROW
      REAL               AJJ
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GRIDINFO, CHK1MAT, CGEMV,
     $                   CLACGV, CSSCAL, IGEBR2D, IGEBS2D,
     $                   INFOG2L, PB_TOPGET, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD, REAL, SQRT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX            CDOTC
      EXTERNAL           LSAME, CDOTC
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters.
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(600+CTXT_)
      ELSE
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 6, INFO )
         IF( INFO.EQ.0 ) THEN
            UPPER = LSAME( UPLO, 'U' )
            IROFF = MOD( IA-1, DESCA( MB_ ) )
            ICOFF = MOD( JA-1, DESCA( NB_ ) )
            IF ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( N+ICOFF.GT.DESCA( NB_ ) ) THEN
               INFO = -2
            ELSE IF( IROFF.NE.0 ) THEN
               INFO = -4
            ELSE IF( ICOFF.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(600+NB_)
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCPOTF2', -INFO )
         CALL BLACS_ABORT( ICTXT, 1 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Compute local information
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
      CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
*
      IF ( UPPER ) THEN
*
*        Process (IAROW, IACOL) owns block to be factorized
*
         IF( MYROW.EQ.IAROW ) THEN
            IF( MYCOL.EQ.IACOL ) THEN
*
*              Compute the Cholesky factorization A = U'*U.
*
               LDA = DESCA( LLD_ )
               IDIAG = IIA + ( JJA - 1 ) * LDA
               IOFFA = IDIAG
*
               DO 10 J = JA, JA+N-1
*
*                 Compute U(J,J) and test for non-positive-definiteness.
*
                  AJJ = REAL( A( IDIAG ) ) -
     $                  CDOTC( J-JA, A( IOFFA ), 1, A( IOFFA ), 1 )
                  IF( AJJ.LE.ZERO ) THEN
                     A( IDIAG ) = AJJ
                     INFO = J - JA + 1
                     GO TO 20
                  END IF
                  AJJ = SQRT( AJJ )
                  A( IDIAG ) = AJJ
*
*                 Compute elements J+1:JA+N-1 of row J.
*
                  IF( J.LT.JA+N-1 ) THEN
                     ICURR = IDIAG + LDA
                     CALL CLACGV( J-JA, A( IOFFA ), 1 )
                     CALL CGEMV( 'Transpose', J-JA, JA+N-J-1, -CONE,
     $                           A( IOFFA+LDA ), LDA, A( IOFFA ), 1,
     $                           CONE, A( ICURR ), LDA )
                     CALL CLACGV( J-JA, A( IOFFA ), 1 )
                     CALL CSSCAL( JA+N-J-1, ONE / AJJ, A( ICURR ),
     $                            LDA )
                  END IF
                  IDIAG = IDIAG + LDA + 1
                  IOFFA = IOFFA + LDA
   10          CONTINUE
*
   20          CONTINUE
*
*              Broadcast INFO to all processes in my IAROW.
*
               CALL IGEBS2D( ICTXT, 'Rowwise', ROWBTOP, 1, 1, INFO, 1 )
*
            ELSE
*
               CALL IGEBR2D( ICTXT, 'Rowwise', ROWBTOP, 1, 1, INFO, 1,
     $                       MYROW, IACOL )
            END IF
*
*           IAROW bcasts along columns so that everyone has INFO
*
            CALL IGEBS2D( ICTXT, 'Columnwise', COLBTOP, 1, 1, INFO, 1 )
*
         ELSE
*
            CALL IGEBR2D( ICTXT, 'Columnwise', COLBTOP, 1, 1, INFO, 1,
     $                    IAROW, MYCOL )
*
         END IF
*
      ELSE
*
*        Process (IAROW, IACOL) owns block to be factorized
*
         IF( MYCOL.EQ.IACOL ) THEN
            IF( MYROW.EQ.IAROW ) THEN
*
*              Compute the Cholesky factorization A = L*L'.
*
               LDA = DESCA( LLD_ )
               IDIAG = IIA + ( JJA - 1 ) * LDA
               IOFFA = IDIAG
*
               DO 30 J = JA, JA+N-1
*
*                 Compute L(J,J) and test for non-positive-definiteness.
*
                  AJJ = REAL( A( IDIAG ) ) -
     $                  CDOTC( J-JA, A( IOFFA ), LDA, A( IOFFA ), LDA )
                  IF ( AJJ.LE.ZERO ) THEN
                     A( IDIAG ) = AJJ
                     INFO = J - JA + 1
                     GO TO 40
                  END IF
                  AJJ = SQRT( AJJ )
                  A( IDIAG ) = AJJ
*
*                 Compute elements J+1:JA+N-1 of column J.
*
                  IF( J.LT.JA+N-1 ) THEN
                     ICURR = IDIAG + 1
                     CALL CLACGV( J-JA, A( IOFFA ), LDA )
                     CALL CGEMV( 'No transpose', JA+N-J-1, J-JA, -CONE,
     $                           A( IOFFA+1 ), LDA, A( IOFFA ), LDA,
     $                           CONE, A( ICURR ), 1 )
                     CALL CLACGV( J-JA, A( IOFFA ), LDA )
                     CALL CSSCAL( JA+N-J-1, ONE / AJJ, A( ICURR ), 1 )
                  END IF
                  IDIAG = IDIAG + LDA + 1
                  IOFFA = IOFFA + 1
   30          CONTINUE
*
   40          CONTINUE
*
*              Broadcast INFO to everyone in IACOL
*
               CALL IGEBS2D( ICTXT, 'Columnwise', COLBTOP, 1, 1, INFO,
     $                       1 )
*
            ELSE
*
               CALL IGEBR2D( ICTXT, 'Columnwise', COLBTOP, 1, 1, INFO,
     $                       1, IAROW, MYCOL )
*
            END IF
*
*           IACOL bcasts INFO along rows so that everyone has it
*
            CALL IGEBS2D( ICTXT, 'Rowwise', ROWBTOP, 1, 1, INFO, 1 )
*
         ELSE
*
            CALL IGEBR2D( ICTXT, 'Rowwise', ROWBTOP, 1, 1, INFO, 1,
     $                    MYROW, IACOL )
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PCPOTF2
*
      END
