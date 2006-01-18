      SUBROUTINE PCLASET( UPLO, M, N, ALPHA, BETA, A, IA, JA, DESCA )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, JA, M, N
      COMPLEX            ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX            A( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLASET initializes an M-by-N distributed matrix sub( A ) denoting
*  A(IA:IA+M-1,JA:JA+N-1) to BETA on the diagonal and ALPHA on the
*  offdiagonals.
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
*          Specifies the part of the distributed matrix sub( A ) to be
*          set:
*          = 'U':      Upper triangular part is set; the strictly lower
*                      triangular part of sub( A ) is not changed;
*          = 'L':      Lower triangular part is set; the strictly upper
*                      triangular part of sub( A ) is not changed;
*          Otherwise:  All of the matrix sub( A ) is set.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  ALPHA   (global input) COMPLEX
*          The constant to which the offdiagonal elements are to be
*          set.
*
*  BETA    (global input) COMPLEX
*          The constant to which the diagonal elements are to be set.
*
*  A       (local output) COMPLEX pointer into the local memory
*          to an array of dimension (LLD_A,LOCc(JA+N-1)).  This array
*          contains the local pieces of the distributed matrix sub( A )
*          to be set.  On exit, the leading M-by-N submatrix sub( A )
*          is set as follows:
*
*          if UPLO = 'U', A(IA+i-1,JA+j-1) = ALPHA, 1<=i<=j-1, 1<=j<=N,
*          if UPLO = 'L', A(IA+i-1,JA+j-1) = ALPHA, j+1<=i<=M, 1<=j<=N,
*          otherwise,     A(IA+i-1,JA+j-1) = ALPHA, 1<=i<=M, 1<=j<=N,
*                                                   IA+i.NE.JA+j,
*          and, for all UPLO, A(IA+i-1,JA+i-1) = BETA, 1<=i<=min(M,N).
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
*     ..
*     .. Local Scalars ..
      INTEGER            I, IAA, IBLK, IN, ITMP, J, JAA, JBLK, JN, JTMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           PCLASE2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL
      EXTERNAL           ICEIL, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      IF( M.LE.( DESCA( MB_ ) - MOD( IA-1, DESCA( MB_ ) ) ) .OR.
     $    N.LE.( DESCA( NB_ ) - MOD( JA-1, DESCA( NB_ ) ) ) ) THEN
         CALL PCLASE2( UPLO, M, N, ALPHA, BETA, A, IA, JA, DESCA )
      ELSE
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+M-1 )
            CALL PCLASE2( UPLO, IN-IA+1, N, ALPHA, BETA, A, IA, JA,
     $                    DESCA )
            DO 10 I = IN+1, IA+M-1, DESCA( MB_ )
               ITMP = I-IA
               IBLK = MIN( DESCA( MB_ ), M-ITMP )
               JAA = JA + ITMP
               CALL PCLASE2( UPLO, IBLK, N-ITMP, ALPHA, BETA,
     $                       A, I, JAA, DESCA )
   10       CONTINUE
         ELSE IF( LSAME( UPLO, 'L' ) ) THEN
            JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
            CALL PCLASE2( UPLO, M, JN-JA+1, ALPHA, BETA, A, IA, JA,
     $                    DESCA )
            DO 20 J = JN+1, JA+N-1, DESCA( NB_ )
               JTMP = J-JA
               JBLK = MIN( DESCA( NB_ ), N-JTMP )
               IAA = IA + JTMP
               CALL PCLASE2( UPLO, M-JTMP, JBLK, ALPHA, BETA, A, IAA,
     $                       J, DESCA )
   20       CONTINUE
         ELSE
            IF( M.LE.N ) THEN
               IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ),
     $                   IA+M-1 )
               CALL PCLASE2( UPLO, IN-IA+1, N, ALPHA, BETA, A, IA,
     $                       JA, DESCA )
               DO 30 I = IN+1, IA+M-1, DESCA( MB_ )
                  ITMP = I-IA
                  IBLK = MIN( DESCA( MB_ ), M-ITMP )
                  CALL PCLASE2( UPLO, IBLK, I-IA, ALPHA, ALPHA, A, I,
     $                          JA, DESCA )
                  CALL PCLASE2( UPLO, IBLK, N-I+IA, ALPHA, BETA, A, I,
     $                          JA+I-IA, DESCA )
   30          CONTINUE
            ELSE
               JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ),
     $                   JA+N-1 )
               CALL PCLASE2( UPLO, M, JN-JA+1, ALPHA, BETA, A, IA,
     $                       JA, DESCA )
               DO 40 J = JN+1, JA+N-1, DESCA( NB_ )
                  JTMP = J-JA
                  JBLK = MIN( DESCA( NB_ ), N-JTMP )
                  CALL PCLASE2( UPLO, J-JA, JBLK, ALPHA, ALPHA, A, IA,
     $                          J, DESCA )
                  CALL PCLASE2( UPLO, M-J+JA, JBLK, ALPHA, BETA, A,
     $                          IA+J-JA, J, DESCA )
   40          CONTINUE
            END IF
         END IF
*
      END IF
*
      RETURN
*
*     End of PCLASET
*
      END
