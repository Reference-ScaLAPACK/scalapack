      SUBROUTINE PSLACPY( UPLO, M, N, A, IA, JA, DESCA, B, IB, JB,
     $                    DESCB )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, IB, JA, JB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      REAL               A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLACPY copies all or part of a distributed matrix A to another
*  distributed matrix B.  No communication is performed, PSLACPY
*  performs a local copy sub( A ) := sub( B ), where sub( A ) denotes
*  A(IA:IA+M-1,JA:JA+N-1) and sub( B ) denotes B(IB:IB+M-1,JB:JB+N-1).
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
*          copied:
*          = 'U':   Upper triangular part is copied; the strictly
*                   lower triangular part of sub( A ) is not referenced;
*          = 'L':   Lower triangular part is copied; the strictly
*                   upper triangular part of sub( A ) is not referenced;
*          Otherwise:  All of the matrix sub( A ) is copied.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input) REAL pointer into the local memory
*          to an array of dimension (LLD_A, LOCc(JA+N-1) ). This array
*          contains the local pieces of the distributed matrix sub( A )
*          to be copied from.
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
*  B       (local output) REAL pointer into the local memory
*          to an array of dimension (LLD_B, LOCc(JB+N-1) ). This array
*          contains on exit the local pieces of the distributed matrix
*          sub( B ) set as follows:
*
*          if UPLO = 'U', B(IB+i-1,JB+j-1) = A(IA+i-1,JA+j-1),
*                         1<=i<=j, 1<=j<=N;
*          if UPLO = 'L', B(IB+i-1,JB+j-1) = A(IA+i-1,JA+j-1),
*                         j<=i<=M, 1<=j<=N;
*          otherwise,     B(IB+i-1,JB+j-1) = A(IA+i-1,JA+j-1),
*                         1<=i<=M, 1<=j<=N.
*
*  IB      (global input) INTEGER
*          The row index in the global array B indicating the first
*          row of sub( B ).
*
*  JB      (global input) INTEGER
*          The column index in the global array B indicating the
*          first column of sub( B ).
*
*  DESCB   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix B.
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
      INTEGER            I, IAA, IBB, IBLK, IN, ITMP, J, JAA, JBB,
     $                   JBLK, JN, JTMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           PSLACP2
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
      IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+M-1 )
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
*
      IF( M.LE.( DESCA( MB_ ) - MOD( IA-1, DESCA( MB_ ) ) ) .OR.
     $    N.LE.( DESCA( NB_ ) - MOD( JA-1, DESCA( NB_ ) ) ) ) THEN
         CALL PSLACP2( UPLO, M, N, A, IA, JA, DESCA,
     $                 B, IB, JB, DESCB )
      ELSE
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            CALL PSLACP2( UPLO, IN-IA+1, N, A, IA, JA, DESCA,
     $                    B, IB, JB, DESCB )
            DO 10 I = IN+1, IA+M-1, DESCA( MB_ )
               ITMP = I-IA
               IBLK = MIN( DESCA( MB_ ), M-ITMP )
               IBB = IB + ITMP
               JBB = JB + ITMP
               JAA = JA + ITMP
               CALL PSLACP2( UPLO, IBLK, N-ITMP, A, I, JAA, DESCA,
     $                       B, IBB, JBB, DESCB )
   10       CONTINUE
         ELSE IF( LSAME( UPLO, 'L' ) ) THEN
            CALL PSLACP2( UPLO, M, JN-JA+1, A, IA, JA, DESCA,
     $                    B, IB, JB, DESCB )
            DO 20 J = JN+1, JA+N-1, DESCA( NB_ )
               JTMP = J-JA
               JBLK = MIN( DESCA( NB_ ), N-JTMP )
               IBB = IB + JTMP
               JBB = JB + JTMP
               IAA = IA + JTMP
               CALL PSLACP2( UPLO, M-JTMP, JBLK, A, IAA, J, DESCA,
     $                       B, IBB, JBB, DESCB )
   20       CONTINUE
         ELSE
            IF( M.LE.N ) THEN
               CALL PSLACP2( UPLO, IN-IA+1, N, A, IA, JA, DESCA,
     $                       B, IB, JB, DESCB )
               DO 30 I = IN+1, IA+M-1, DESCA( MB_ )
                  ITMP = I-IA
                  IBLK = MIN( DESCA( MB_ ), M-ITMP )
                  IBB = IB+ITMP
                  CALL PSLACP2( UPLO, IBLK, N, A, I, JA, DESCA,
     $                          B, IBB, JB, DESCB )
   30          CONTINUE
            ELSE
               CALL PSLACP2( UPLO, M, JN-JA+1, A, IA, JA, DESCA,
     $                       B, IB, JB, DESCB )
               DO 40 J = JN+1, JA+N-1, DESCA( NB_ )
                  JTMP = J-JA
                  JBLK = MIN( DESCA( NB_ ), N-JTMP )
                  JBB = JB+JTMP
                  CALL PSLACP2( UPLO, M, JBLK, A, IA, J, DESCA,
     $                          B, IB, JBB, DESCB )
   40          CONTINUE
            END IF
         END IF
*
      END IF
*
      RETURN
*
*     End of PSLACPY
*
      END
