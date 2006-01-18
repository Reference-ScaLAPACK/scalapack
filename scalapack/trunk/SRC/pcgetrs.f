      SUBROUTINE PCGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     $                    IB, JB, DESCB, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      COMPLEX            A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PCGETRS solves a system of distributed linear equations
*
*                   op( sub( A ) ) * X = sub( B )
*
*  with a general N-by-N distributed matrix sub( A ) using the LU
*  factorization computed by PCGETRF.
*  sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1), op( A ) = A, A**T or A**H
*  and sub( B ) denotes B(IB:IB+N-1,JB:JB+NRHS-1).
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
*  This routine requires square block data decomposition ( MB_A=NB_A ).
*
*  Arguments
*  =========
*
*  TRANS   (global input) CHARACTER
*          Specifies the form of the system of equations:
*          = 'N':  sub( A )    * X = sub( B )  (No transpose)
*          = 'T':  sub( A )**T * X = sub( B )  (Transpose)
*          = 'C':  sub( A )**H * X = sub( B )  (Conjugate transpose)
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the distributed submatrix sub( B ). NRHS >= 0.
*
*  A       (local input) COMPLEX pointer into the local
*          memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the factors
*          L and U from the factorization sub( A ) = P*L*U; the unit
*          diagonal elements of L are not stored.
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
*  IPIV    (local input) INTEGER array, dimension ( LOCr(M_A)+MB_A )
*          This array contains the pivoting information.
*          IPIV(i) -> The global row local row i was swapped with.
*          This array is tied to the distributed matrix A.
*
*  B       (local input/local output) COMPLEX pointer into the
*          local memory to an array of dimension
*          (LLD_B,LOCc(JB+NRHS-1)).  On entry, the right hand sides
*          sub( B ). On exit, sub( B ) is overwritten by the solution
*          distributed matrix X.
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
*  INFO    (global output) INTEGER
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
      COMPLEX            ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
      INTEGER            IAROW, IBROW, ICOFFA, ICTXT, IROFFA, IROFFB,
     $                   MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            DESCIP( DLEN_ ), IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, PCHK2MAT,
     $                   PCLAPIV, PCTRSM, PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, MOD
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
         NOTRAN = LSAME( TRANS, 'N' )
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 2, NRHS, 3, IB, JB, DESCB, 12, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $                       NPROW )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
               INFO = -1
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(700+NB_)
            ELSE IF( IROFFB.NE.0 .OR. IBROW.NE.IAROW ) THEN
               INFO = -10
            ELSE IF( DESCB( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(1200+NB_)
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -(1200+CTXT_)
            END IF
         END IF
         IF( NOTRAN ) THEN
            IDUM1( 1 ) = ICHAR( 'N' )
         ELSE IF( LSAME( TRANS, 'T' ) ) THEN
            IDUM1( 1 ) = ICHAR( 'T' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'C' )
         END IF
         IDUM2( 1 ) = 1
         CALL PCHK2MAT( N, 2, N, 2, IA, JA, DESCA, 7, N, 2, NRHS, 3,
     $                  IB, JB, DESCB, 12, 1, IDUM1, IDUM2, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      CALL DESCSET( DESCIP, DESCA( M_ ) + DESCA( MB_ )*NPROW, 1,
     $              DESCA( MB_ ), 1, DESCA( RSRC_ ), MYCOL, ICTXT,
     $              DESCA( MB_ ) + NUMROC( DESCA( M_ ), DESCA( MB_ ),
     $              MYROW, DESCA( RSRC_ ), NPROW ) )
*
      IF( NOTRAN ) THEN
*
*        Solve sub( A ) * X = sub( B ).
*
*        Apply row interchanges to the right hand sides.
*
         CALL PCLAPIV( 'Forward', 'Row', 'Col', N, NRHS, B, IB, JB,
     $                 DESCB, IPIV, IA, 1, DESCIP, IDUM1 )
*
*        Solve L*X = sub( B ), overwriting sub( B ) with X.
*
         CALL PCTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $                ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
*
*        Solve U*X = sub( B ), overwriting sub( B ) with X.
*
         CALL PCTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $                NRHS, ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
      ELSE
*
*        Solve sub( A )' * X = sub( B ).
*
*        Solve U'*X = sub( B ), overwriting sub( B ) with X.
*
         CALL PCTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS,
     $                ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
*
*        Solve L'*X = sub( B ), overwriting sub( B ) with X.
*
         CALL PCTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE,
     $                A, IA, JA, DESCA, B, IB, JB, DESCB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL PCLAPIV( 'Backward', 'Row', 'Col', N, NRHS, B, IB, JB,
     $                 DESCB, IPIV, IA, 1, DESCIP, IDUM1 )
*
      END IF
*
      RETURN
*
*     End of PCGETRS
*
      END
