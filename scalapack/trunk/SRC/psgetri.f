      SUBROUTINE PSGETRI( N, A, IA, JA, DESCA, IPIV, WORK, LWORK,
     $                    IWORK, LIWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7.4) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     v1.7.4: May 10, 2006 
*     v1.7:   May 1,  1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), IPIV( * ), IWORK( * )
      REAL               A( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSGETRI computes the inverse of a distributed matrix using the LU
*  factorization computed by PSGETRF. This method inverts U and then
*  computes the inverse of sub( A ) = A(IA:IA+N-1,JA:JA+N-1) denoted
*  InvA by solving the system InvA*L = inv(U) for InvA.
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
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) REAL pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, the local pieces of the L and U obtained by the
*          factorization sub( A ) = P*L*U computed by PSGETRF. On
*          exit, if INFO = 0, sub( A ) contains the inverse of the
*          original distributed matrix sub( A ).
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
*  IPIV    (local input) INTEGER array, dimension LOCr(M_A)+MB_A
*          keeps track of the pivoting information. IPIV(i) is the
*          global row index the local row i was swapped with.  This
*          array is tied to the distributed matrix A.
*
*  WORK    (local workspace/local output) REAL array,
*                                                     dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK = LOCr(N+MOD(IA-1,MB_A))*NB_A. WORK is used to keep a
*          copy of at most an entire column block of sub( A ).
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  IWORK   (local workspace/local output) INTEGER array,
*                                                    dimension (LIWORK)
*          On exit, IWORK(1) returns the minimal and optimal LIWORK.
*
*  LIWORK  (local or global input) INTEGER
*          The dimension of the array IWORK used as workspace for
*          physically transposing the pivots.
*          LIWORK is local input and must be at least
*          if NPROW == NPCOL then
*            LIWORK = LOCc( N_A + MOD(JA-1, NB_A) ) + NB_A,
*          else
*            LIWORK =  LOCc( N_A + MOD(JA-1, NB_A) ) +
*                      MAX( CEIL(CEIL(LOCr(M_A)/MB_A)/(LCM/NPROW)),
*                           NB_A )
*              where LCM is the least common multiple of process
*              rows and columns (NPROW and NPCOL).
*          end if
*
*          If LIWORK = -1, then LIWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*            > 0:  If INFO = K, U(IA+K-1,IA+K-1) is exactly zero; the
*                  matrix is singular and its inverse could not be
*                  computed.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IACOL, IAROW, ICOFF, ICTXT, IROFF, IW, J,
     $                   JB, JN, LCM, LIWMIN, LWMIN, MP, MYCOL, MYROW,
     $                   NN, NP, NPCOL, NPROW, NQ
*     ..
*     .. Local Arrays ..
      INTEGER            DESCW( DLEN_ ), IDUM1( 2 ), IDUM2( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, PCHK1MAT,
     $                   PSGEMM, PSLACPY, PSLASET, PSLAPIV,
     $                   PSTRSM, PSTRTRI, PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, ILCM, INDXG2P, NUMROC
      EXTERNAL           ICEIL, ILCM, INDXG2P, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD, REAL
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
         INFO = -(500+CTXT_)
      ELSE
         CALL CHK1MAT( N, 1, N, 1, IA, JA, DESCA, 5, INFO )
         IF( INFO.EQ.0 ) THEN
            IROFF = MOD( IA-1, DESCA( MB_ ) )
            ICOFF = MOD( JA-1, DESCA( NB_ ) )
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            NP = NUMROC( N+IROFF, DESCA( MB_ ), MYROW, IAROW, NPROW )
            LWMIN = NP * DESCA( NB_ )
*
            MP = NUMROC( DESCA( M_ ), DESCA( MB_ ), MYROW,
     $                   DESCA( RSRC_ ), NPROW )
            NQ = NUMROC( DESCA( N_ ), DESCA( NB_ ), MYCOL,
     $                   DESCA( CSRC_ ), NPCOL )
            IF( NPROW.EQ.NPCOL ) THEN
               LIWMIN = NQ + DESCA( NB_ )
            ELSE
*
* Use the formula for the workspace given in PxLAPIV
* to compute the minimum size LIWORK for IWORK
*
* The formula in PxLAPIV is
*   LDW = LOCc( M_P + MOD(IP-1, MB_P) ) +
*         MB_P * CEIL( CEIL(LOCr(M_P)/MB_P) / (LCM/NPROW) )
*
* where 
*   M_P     is the global length of the pivot vector
*           MP = DESCA( M_ ) + DESCA( MB_ ) * NPROW
*   I_P     is IA
*           I_P = IA
*   MB_P    is the block size use for the block cyclic distribution of the 
*           pivot vector
*           MB_P = DESCA (MB_ )
*   LOCc ( . ) 
*           NUMROC ( . , DESCA ( NB_ ), MYCOL, DESCA ( CSRC_ ), NPCOL )
*   LOCr ( . )
*           NUMROC ( . , DESCA ( MB_ ), MYROW, DESCA ( RSRC_ ), NPROW )
*   CEIL ( X / Y )
*           ICEIL( X, Y )
*   LCM 
*           LCM = ILCM( NPROW, NPCOL )
*
               LCM = ILCM( NPROW, NPCOL )
               LIWMIN = NUMROC( DESCA( M_ ) + DESCA( MB_ ) * NPROW
     $                  + MOD ( IA - 1, DESCA( MB_ ) ), DESCA ( NB_ ),
     $                  MYCOL, DESCA( CSRC_ ), NPCOL ) +
     $                  MAX ( DESCA( MB_ ) * ICEIL ( ICEIL(
     $                  NUMROC( DESCA( M_ ) + DESCA( MB_ ) * NPROW,
     $                  DESCA( MB_ ), MYROW, DESCA( RSRC_ ), NPROW ),
     $                  DESCA( MB_ ) ), LCM / NPROW ), DESCA( NB_ ) )
*
            END IF
*
            WORK( 1 ) = REAL( LWMIN )
            IWORK( 1 ) = LIWMIN
            LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
            IF( IROFF.NE.ICOFF .OR. IROFF.NE.0 ) THEN
               INFO = -4
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(500+NB_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -8
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -10
            END IF
         END IF
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 1 ) = -1
         ELSE
            IDUM1( 1 ) = 1
         END IF
         IDUM2( 1 ) = 8
         IF( LIWORK.EQ.-1 ) THEN
            IDUM1( 2 ) = -1
         ELSE
            IDUM1( 2 ) = 1
         END IF
         IDUM2( 2 ) = 10
         CALL PCHK1MAT( N, 1, N, 1, IA, JA, DESCA, 5, 2, IDUM1, IDUM2,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PSGETRI', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form inv(U).  If INFO > 0 from PSTRTRI, then U is singular,
*     and the inverse is not computed.
*
      CALL PSTRTRI( 'Upper', 'Non-unit', N, A, IA, JA, DESCA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
*     Define array descriptor for working array WORK
*
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      NN = ( ( JA+N-2 ) / DESCA( NB_ ) ) * DESCA( NB_ ) + 1
      IACOL = INDXG2P( NN, DESCA( NB_ ), MYCOL, DESCA( CSRC_ ), NPCOL )
      CALL DESCSET( DESCW, N+IROFF, DESCA( NB_ ), DESCA( MB_ ),
     $              DESCA( NB_ ), IAROW, IACOL, ICTXT, MAX( 1, NP ) )
      IW = IROFF + 1
*
*     Solve the equation inv(A)*L=inv(U) for inv(A) using blocked code.
*
      DO 10 J = NN, JN+1, -DESCA( NB_ )
         JB = MIN( DESCA( NB_ ), JA+N-J )
         I = IA + J - JA
*
*        Copy current block column of L to WORK and replace with zeros.
*
         CALL PSLACPY( 'Lower', JA+N-1-J, JB, A, I+1, J, DESCA,
     $                 WORK, IW+J-JA+1, 1, DESCW )
         CALL PSLASET( 'Lower', JA+N-1-J, JB, ZERO, ZERO, A, I+1, J,
     $                 DESCA )
*
*        Compute current block column of inv(A).
*
         IF( J+JB.LE.JA+N-1 )
     $      CALL PSGEMM( 'No transpose', 'No transpose', N, JB,
     $                   JA+N-J-JB, -ONE, A, IA, J+JB, DESCA, WORK,
     $                   IW+J+JB-JA, 1, DESCW, ONE, A, IA, J, DESCA )
         CALL PSTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                ONE, WORK, IW+J-JA, 1, DESCW, A, IA, J, DESCA )
         DESCW( CSRC_ ) = MOD( DESCW( CSRC_ ) + NPCOL - 1, NPCOL )
*
   10 CONTINUE
*
*     Handle the last block of columns separately
*
      JB = JN-JA+1
*
*     Copy current block column of L to WORK and replace with zeros.
*
      CALL PSLACPY( 'Lower', N-1, JB, A, IA+1, JA, DESCA, WORK, IW+1,
     $              1, DESCW )
      CALL PSLASET( 'Lower', N-1, JB, ZERO, ZERO, A, IA+1, JA, DESCA )
*
*     Compute current block column of inv(A).
*
      IF( JA+JB.LE.JA+N-1 )
     $   CALL PSGEMM( 'No transpose', 'No transpose', N, JB,
     $                N-JB, -ONE, A, IA, JA+JB, DESCA, WORK, IW+JB, 1,
     $                DESCW, ONE, A, IA, JA, DESCA )
      CALL PSTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $             ONE, WORK, IW, 1, DESCW, A, IA, JA, DESCA )
*
*     Use the row pivots and apply them to the columns of the global
*     matrix.
*
      CALL DESCSET( DESCW, DESCA( M_ ) + DESCA( MB_ )*NPROW, 1,
     $              DESCA( MB_ ), 1, DESCA( RSRC_ ), MYCOL, ICTXT,
     $              MP+DESCA( MB_ ) )
      CALL PSLAPIV( 'Backward', 'Columns', 'Column', N, N, A, IA,
     $              JA, DESCA, IPIV, IA, 1, DESCW, IWORK )
*
      WORK( 1 ) = REAL( LWMIN )
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of PSGETRI
*
      END
