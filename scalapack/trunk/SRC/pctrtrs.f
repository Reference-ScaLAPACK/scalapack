      SUBROUTINE PCTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     $                    B, IB, JB, DESCB, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX            A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PCTRTRS solves a triangular system of the form
*
*     sub( A ) * X = sub( B )  or  sub( A )**T * X = sub( B ) or
*
*     sub( A )**H * X = sub( B ),
*
*  where sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1) and is a triangular
*  distributed matrix of order N, and B(IB:IB+N-1,JB:JB+NRHS-1) is an
*  N-by-NRHS distributed matrix denoted by sub( B ). A check is made
*  to verify that sub( A ) is nonsingular.
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
*          = 'U':  sub( A ) is upper triangular;
*          = 'L':  sub( A ) is lower triangular.
*
*  TRANS   (global input) CHARACTER
*          Specifies the form of the system of equations:
*          = 'N': Solve sub( A )    * X = sub( B ) (No transpose)
*          = 'T': Solve sub( A )**T * X = sub( B ) (Transpose)
*          = 'C': Solve sub( A )**H * X = sub( B ) (Conjugate transpose)
*
*  DIAG    (global input) CHARACTER
*          = 'N':  sub( A ) is non-unit triangular;
*          = 'U':  sub( A ) is unit triangular.
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on i.e the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the distributed matrix sub( B ). NRHS >= 0.
*
*  A       (local input) COMPLEX pointer into the local memory
*          to an array of dimension (LLD_A,LOCc(JA+N-1) ). This array
*          contains the local pieces of the distributed triangular
*          matrix sub( A ).  If UPLO = 'U', the leading N-by-N upper
*          triangular part of sub( A ) contains the upper triangular
*          matrix, and the strictly lower triangular part of sub( A )
*          is not referenced.  If UPLO = 'L', the leading N-by-N lower
*          triangular part of sub( A ) contains the lower triangular
*          matrix, and the strictly upper triangular part of sub( A )
*          is not referenced.  If DIAG = 'U', the diagonal elements of
*          sub( A ) are also not referenced and are assumed to be 1.
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
*  B       (local input/local output) COMPLEX pointer into the
*          local memory to an array of dimension
*          (LLD_B,LOCc(JB+NRHS-1)).  On entry, this array contains the
*          local pieces of the right hand side distributed matrix
*          sub( B ). On exit, if INFO = 0, sub( B ) is overwritten by
*          the solution matrix X.
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
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = i, the i-th diagonal element of sub( A ) is
*                zero, indicating that the submatrix is singular and the
*                solutions X have not been computed.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IAROW, IBROW, ICOFFA, ICTXT, ICURCOL,
     $                   ICURROW, IROFFA, IROFFB, IDUM, II, IOFFA, J,
     $                   JBLK, JJ, JN, LDA, LL, MYCOL, MYROW, NPCOL,
     $                   NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 3 ), IDUM2( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, IGAMX2D, INFOG2L,
     $                   PCHK2MAT, PCTRSM, PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, INDXG2P
      EXTERNAL           ICEIL, INDXG2P, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -907
      ELSE
         UPPER = LSAME( UPLO, 'U' )
         NOUNIT = LSAME( DIAG, 'N' )
         NOTRAN = LSAME( TRANS, 'N' )
*
         CALL CHK1MAT( N, 4, N, 4, IA, JA, DESCA, 9, INFO )
         CALL CHK1MAT( N, 4, NRHS, 5, IB, JB, DESCB, 13, INFO )
         IF( INFO.EQ.0 ) THEN
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $                       NPROW )
            IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.NOTRAN .AND.  .NOT.LSAME( TRANS, 'T' ) .AND.
     $               .NOT.LSAME( TRANS, 'C' ) ) THEN
               INFO = -2
            ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
               INFO = -3
            ELSE IF( IROFFA.NE.ICOFFA .OR. IROFFA.NE.0 ) THEN
               INFO = -8
            ELSE IF( IROFFA.NE.IROFFB .OR. IAROW.NE.IBROW ) THEN
               INFO = -11
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -904
            ELSE IF( DESCB( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -1304
            END IF
         END IF
*
         IF( UPPER ) THEN
            IDUM1( 1 ) = ICHAR( 'U' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'L' )
         END IF
         IDUM2( 1 ) = 1
         IF( NOTRAN ) THEN
            IDUM1( 2 ) = ICHAR( 'N' )
         ELSE IF( LSAME( TRANS, 'T' ) ) THEN
            IDUM1( 2 ) = ICHAR( 'T' )
         ELSE IF( LSAME( TRANS, 'C' ) ) THEN
            IDUM1( 2 ) = ICHAR( 'C' )
         END IF
         IDUM2( 2 ) = 2
         IF( NOUNIT ) THEN
            IDUM1( 3 ) = ICHAR( 'N' )
         ELSE
            IDUM1( 3 ) = ICHAR( 'D' )
         END IF
         IDUM2( 3 ) = 3
         CALL PCHK2MAT( N, 4, N, 4, IA, JA, DESCA, 9, N, 4, NRHS, 5,
     $                  IB, JB, DESCB, 13, 3, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCTRTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      IF( NOUNIT ) THEN
          CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                  II, JJ, ICURROW, ICURCOL )
          JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
          LDA = DESCA( LLD_ )
          IOFFA = II + ( JJ - 1 ) * LDA
*
*         Handle first block separately
*
          JBLK = JN-JA+1
          IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
             LL = IOFFA
             DO 10 I = 0, JBLK-1
                IF( A( LL ).EQ.ZERO .AND. INFO.EQ.0 )
     $             INFO = I + 1
                LL = IOFFA + LDA + 1
   10        CONTINUE
          END IF
          IF( MYROW.EQ.ICURROW )
     $       IOFFA = IOFFA + JBLK
          IF( MYCOL.EQ.ICURCOL )
     $       IOFFA = IOFFA + JBLK*LDA
          ICURROW = MOD( ICURROW+1, NPROW )
          ICURCOL = MOD( ICURCOL+1, NPCOL )
*
*         Loop over remaining blocks of columns
*
          DO 30 J = JN+1, JA+N-1, DESCA( NB_ )
             JBLK = MIN( JA+N-J, DESCA( NB_ ) )
             IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                LL = IOFFA
                DO 20 I = 0, JBLK-1
                   IF( A( LL ).EQ.ZERO .AND. INFO.EQ.0 )
     $                INFO = J + I - JA + 1
                   LL = IOFFA + LDA + 1
   20           CONTINUE
             END IF
             IF( MYROW.EQ.ICURROW )
     $          IOFFA = IOFFA + JBLK
             IF( MYCOL.EQ.ICURCOL )
     $          IOFFA = IOFFA + JBLK*LDA
             ICURROW = MOD( ICURROW+1, NPROW )
             ICURCOL = MOD( ICURCOL+1, NPCOL )
   30     CONTINUE
          CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, IDUM, IDUM,
     $                  -1, -1, MYCOL )
          IF( INFO.NE.0 )
     $       RETURN
      END IF
*
*     Solve A * x = b,  A**T * x = b,  or  A**H * x = b.
*
      CALL PCTRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, IA, JA,
     $             DESCA, B, IB, JB, DESCB )
*
      RETURN
*
*     End of PCTRTRS
*
      END
