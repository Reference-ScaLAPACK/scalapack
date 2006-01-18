      SUBROUTINE PSPOSVX( FACT, UPLO, N, NRHS, A, IA, JA, DESCA, AF,
     $                    IAF, JAF, DESCAF, EQUED, SR, SC, B, IB, JB,
     $                    DESCB, X, IX, JX, DESCX, RCOND, FERR, BERR,
     $                    WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, FACT, UPLO
      INTEGER            IA, IAF, IB, INFO, IX, JA, JAF, JB, JX, LIWORK,
     $                   LWORK, N, NRHS
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCAF( * ), DESCB( * ),
     $                   DESCX( * ), IWORK( * )
      REAL               A( * ), AF( * ),
     $                   B( * ), BERR( * ), FERR( * ),
     $                   SC( * ), SR( * ), WORK( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  PSPOSVX uses the Cholesky factorization A = U**T*U or A = L*L**T to
*  compute the solution to a real system of linear equations
*
*        A(IA:IA+N-1,JA:JA+N-1) * X = B(IB:IB+N-1,JB:JB+NRHS-1),
*
*  where A(IA:IA+N-1,JA:JA+N-1) is an N-by-N matrix and X and
*  B(IB:IB+N-1,JB:JB+NRHS-1) are N-by-NRHS matrices.
*
*  Error bounds on the solution and a condition estimate are also
*  provided.  In the following comments Y denotes Y(IY:IY+M-1,JY:JY+K-1)
*  a M-by-K matrix where Y can be A, AF, B and X.
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
*  Description
*  ===========
*
*  The following steps are performed:
*
*  1. If FACT = 'E', real scaling factors are computed to equilibrate
*     the system:
*        diag(SR) * A * diag(SC) * inv(diag(SC)) * X = diag(SR) * B
*     Whether or not the system will be equilibrated depends on the
*     scaling of the matrix A, but if equilibration is used, A is
*     overwritten by diag(SR)*A*diag(SC) and B by diag(SR)*B.
*
*  2. If FACT = 'N' or 'E', the Cholesky decomposition is used to
*     factor the matrix A (after equilibration if FACT = 'E') as
*        A = U**T* U,  if UPLO = 'U', or
*        A = L * L**T,  if UPLO = 'L',
*     where U is an upper triangular matrix and L is a lower triangular
*     matrix.
*
*  3. The factored form of A is used to estimate the condition number
*     of the matrix A.  If the reciprocal of the condition number is
*     less than machine precision, steps 4-6 are skipped.
*
*  4. The system of equations is solved for X using the factored form
*     of A.
*
*  5. Iterative refinement is applied to improve the computed solution
*     matrix and calculate error bounds and backward error estimates
*     for it.
*
*  6. If equilibration was used, the matrix X is premultiplied by
*     diag(SR) so that it solves the original system before
*     equilibration.
*
*  Arguments
*  =========
*
*  FACT    (global input) CHARACTER
*          Specifies whether or not the factored form of the matrix A is
*          supplied on entry, and if not, whether the matrix A should be
*          equilibrated before it is factored.
*          = 'F':  On entry, AF contains the factored form of A.
*                  If EQUED = 'Y', the matrix A has been equilibrated
*                  with scaling factors given by S.  A and AF will not
*                  be modified.
*          = 'N':  The matrix A will be copied to AF and factored.
*          = 'E':  The matrix A will be equilibrated if necessary, then
*                  copied to AF and factored.
*
*  UPLO    (global input) CHARACTER
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix A(IA:IA+N-1,JA:JA+N-1).
*          N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the distributed submatrices B and X.  NRHS >= 0.
*
*  A       (local input/local output) REAL pointer into
*          the local memory to an array of local dimension
*          ( LLD_A, LOCc(JA+N-1) ).
*          On entry, the symmetric matrix A, except if FACT = 'F' and
*          EQUED = 'Y', then A must contain the equilibrated matrix
*          diag(SR)*A*diag(SC).  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.  A is not modified if
*          FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit.
*
*          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by
*          diag(SR)*A*diag(SC).
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
*  AF      (local input or local output) REAL pointer
*          into the local memory to an array of local dimension
*          ( LLD_AF, LOCc(JA+N-1)).
*          If FACT = 'F', then AF is an input argument and on entry
*          contains the triangular factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T, in the same storage
*          format as A.  If EQUED .ne. 'N', then AF is the factored form
*          of the equilibrated matrix diag(SR)*A*diag(SC).
*
*          If FACT = 'N', then AF is an output argument and on exit
*          returns the triangular factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T of the original
*          matrix A.
*
*          If FACT = 'E', then AF is an output argument and on exit
*          returns the triangular factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T of the equilibrated
*          matrix A (see the description of A for the form of the
*          equilibrated matrix).
*
*  IAF     (global input) INTEGER
*          The row index in the global array AF indicating the first
*          row of sub( AF ).
*
*  JAF     (global input) INTEGER
*          The column index in the global array AF indicating the
*          first column of sub( AF ).
*
*  DESCAF  (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix AF.
*
*  EQUED   (global input/global output) CHARACTER
*          Specifies the form of equilibration that was done.
*          = 'N':  No equilibration (always true if FACT = 'N').
*          = 'Y':  Equilibration was done, i.e., A has been replaced by
*                  diag(SR) * A * diag(SC).
*          EQUED is an input variable if FACT = 'F'; otherwise, it is an
*          output variable.
*
*  SR      (local input/local output) REAL array,
*                                                    dimension (LLD_A)
*          The scale factors for A distributed across process rows;
*          not accessed if EQUED = 'N'.  SR is an input variable if
*          FACT = 'F'; otherwise, SR is an output variable.
*          If FACT = 'F' and EQUED = 'Y', each element of SR must be
*          positive.
*
*  SC      (local input/local output) REAL array,
*                                              dimension (LOC(N_A))
*          The scale factors for A distributed across
*          process columns; not accessed if EQUED = 'N'. SC is an input
*          variable if FACT = 'F'; otherwise, SC is an output variable.
*          If FACT = 'F' and EQUED = 'Y', each element of SC must be
*          positive.
*
*  B       (local input/local output) REAL pointer into
*          the local memory to an array of local dimension
*          ( LLD_B, LOCc(JB+NRHS-1) ).
*          On entry, the N-by-NRHS right-hand side matrix B.
*          On exit, if EQUED = 'N', B is not modified; if TRANS = 'N'
*          and EQUED = 'R' or 'B', B is overwritten by diag(R)*B; if
*          TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is overwritten
*          by diag(C)*B.
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
*  X       (local input/local output) REAL pointer into
*          the local memory to an array of local dimension
*          ( LLD_X, LOCc(JX+NRHS-1) ).
*          If INFO = 0, the N-by-NRHS solution matrix X to the original
*          system of equations.  Note that A and B are modified on exit
*          if EQUED .ne. 'N', and the solution to the equilibrated
*          system is inv(diag(SC))*X if TRANS = 'N' and EQUED = 'C' or
*          'B', or inv(diag(SR))*X if TRANS = 'T' or 'C' and EQUED = 'R'
*          or 'B'.
*
*  IX      (global input) INTEGER
*          The row index in the global array X indicating the first
*          row of sub( X ).
*
*  JX      (global input) INTEGER
*          The column index in the global array X indicating the
*          first column of sub( X ).
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*
*  RCOND   (global output) REAL
*          The estimate of the reciprocal condition number of the matrix
*          A after equilibration (if done).  If RCOND is less than the
*          machine precision (in particular, if RCOND = 0), the matrix
*          is singular to working precision.  This condition is
*          indicated by a return code of INFO > 0, and the solution and
*          error bounds are not computed.
*
*  FERR    (local output) REAL array, dimension (LOC(N_B))
*          The estimated forward error bounds for each solution vector
*          X(j) (the j-th column of the solution matrix X).
*          If XTRUE is the true solution, FERR(j) bounds the magnitude
*          of the largest entry in (X(j) - XTRUE) divided by
*          the magnitude of the largest entry in X(j).  The quality of
*          the error bound depends on the quality of the estimate of
*          norm(inv(A)) computed in the code; if the estimate of
*          norm(inv(A)) is accurate, the error bound is guaranteed.
*
*  BERR    (local output) REAL array, dimension (LOC(N_B))
*          The componentwise relative backward error of each solution
*          vector X(j) (i.e., the smallest relative change in
*          any entry of A or B that makes X(j) an exact solution).
*
*  WORK    (local workspace/local output) REAL array,
*                                                dimension (LWORK)
*          On exit, WORK(1) returns the minimal and optimal LWORK.
*
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK is local input and must be at least
*          LWORK = MAX( PSPOCON( LWORK ), PSPORFS( LWORK ) )
*                  + LOCr( N_A ).
*          LWORK = 3*DESCA( LLD_ )
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  IWORK   (local workspace/local output) INTEGER array,
*                                                  dimension (LIWORK)
*          On exit, IWORK(1) returns the minimal and optimal LIWORK.
*
*  LIWORK  (local or global input) INTEGER
*          The dimension of the array IWORK.
*          LIWORK is local input and must be at least
*          LIWORK = DESCA( LLD_ )
*          LIWORK = LOCr(N_A).
*
*          If LIWORK = -1, then LIWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*
*  INFO    (global output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, and i is
*               <= N: if INFO = i, the leading minor of order i of A
*                     is not positive definite, so the factorization
*                     could not be completed, and the solution and error
*                     bounds could not be computed.
*               = N+1: RCOND is less than machine precision.  The
*                     factorization has been completed, but the matrix
*                     is singular to working precision, and the solution
*                     and error bounds have not been computed.
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
*     ..
*     .. Local Scalars ..
      LOGICAL            EQUIL, LQUERY, NOFACT, RCEQU
      INTEGER            I, IACOL, IAROW, IAFROW, IBROW, IBCOL, ICOFF,
     $                   ICOFFA, ICTXT, IDUMM, IIA, IIB, IIX, INFEQU,
     $                   IROFF, IROFFA, IROFFAF, IROFFB, IROFFX, IXCOL,
     $                   IXROW, J, JJA, JJB, JJX, LDB, LDX, LIWMIN,
     $                   LWMIN, MYCOL, MYROW, NP, NPCOL, NPROW, NRHSQ,
     $                   NQ
      REAL               AMAX, ANORM, BIGNUM, SCOND, SMAX, SMIN, SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 5 ), IDUM2( 5 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PCHK2MAT, INFOG2L,
     $                   PSPOCON, PSPOEQU, PSPORFS,
     $                   PSPOTRF, PSPOTRS,
     $                   PSLACPY, PSLAQSY, PXERBLA,
     $                   SGAMN2D, SGAMX2D
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      REAL               PSLANSY, PSLAMCH
      EXTERNAL           INDXG2P, LSAME, NUMROC, PSLANSY, PSLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, MAX, MIN, MOD
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
         INFO = -(800+CTXT_)
      ELSE
         CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 8, INFO )
         IF( LSAME( FACT, 'F' ) )
     $      CALL CHK1MAT( N, 3, N, 3, IAF, JAF, DESCAF, 12, INFO )
         CALL CHK1MAT( N, 3, NRHS, 4, IB, JB, DESCB, 20, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),
     $                       NPROW )
            IAFROW = INDXG2P( IAF, DESCAF( MB_ ), MYROW,
     $                        DESCAF( RSRC_ ), NPROW )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),
     $                       NPROW )
            IXROW = INDXG2P( IX, DESCX( MB_ ), MYROW, DESCX( RSRC_ ),
     $                       NPROW )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            IROFFAF = MOD( IAF-1, DESCAF( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            IROFFX = MOD( IX-1, DESCX( MB_ ) )
            CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW,
     $                    MYCOL, IIA, JJA, IAROW, IACOL )
            NP = NUMROC( N+IROFFA, DESCA( MB_ ), MYROW, IAROW, NPROW )
            IF( MYROW.EQ.IAROW )
     $         NP = NP-IROFFA
            NQ = NUMROC( N+ICOFFA, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
            IF( MYCOL.EQ.IACOL )
     $         NQ = NQ-ICOFFA
            LWMIN = 3*DESCA( LLD_ )
            LIWMIN = NP
            NOFACT = LSAME( FACT, 'N' )
            EQUIL = LSAME( FACT, 'E' )
            IF( NOFACT .OR. EQUIL ) THEN
               EQUED = 'N'
               RCEQU = .FALSE.
            ELSE
               RCEQU = LSAME( EQUED, 'Y' )
               SMLNUM = PSLAMCH( ICTXT, 'Safe minimum' )
               BIGNUM = ONE / SMLNUM
            END IF
            IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND.
     $          .NOT.LSAME( FACT, 'F' ) ) THEN
               INFO = -1
            ELSE IF( .NOT.LSAME( UPLO, 'U' ) .AND.
     $               .NOT.LSAME( UPLO, 'L' ) ) THEN
               INFO = -2
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( ICOFFA.NE.0 .OR. IROFFA.NE.ICOFFA ) THEN
               INFO = -7
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(800+NB_)
            ELSE IF( IAFROW.NE.IAROW ) THEN
               INFO = -10
            ELSE IF( IROFFAF.NE.0 ) THEN
               INFO = -10
            ELSE IF( ICTXT.NE.DESCAF( CTXT_ ) ) THEN
               INFO = -(1200+CTXT_)
            ELSE IF( LSAME( FACT, 'F' ) .AND. .NOT.
     $         ( RCEQU .OR. LSAME( EQUED, 'N' ) ) ) THEN
               INFO = -13
            ELSE
               IF( RCEQU ) THEN
*
                  SMIN = BIGNUM
                  SMAX = ZERO
                  DO 10 J = IIA, IIA + NP - 1
                     SMIN = MIN( SMIN, SR( J ) )
                     SMAX = MAX( SMAX, SR( J ) )
   10             CONTINUE
                  CALL SGAMN2D( ICTXT, 'Columnwise', ' ', 1, 1, SMIN,
     $                          1, IDUMM, IDUMM, -1, -1, MYCOL )
                  CALL SGAMX2D( ICTXT, 'Columnwise', ' ', 1, 1, SMAX,
     $                          1, IDUMM, IDUMM, -1, -1, MYCOL )
                  IF( SMIN.LE.ZERO ) THEN
                     INFO = -14
                  ELSE IF( N.GT.0 ) THEN
                     SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
                  ELSE
                     SCOND = ONE
                  END IF
               END IF
            END IF
         END IF
*
         WORK( 1 ) = REAL( LWMIN )
         IWORK( 1 ) = LIWMIN
         LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
         IF( INFO.EQ.0 ) THEN
            IF( IBROW.NE.IAROW ) THEN
               INFO = -18
            ELSE IF( IXROW.NE.IBROW ) THEN
               INFO = -22
            ELSE IF( DESCB( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(2000+NB_)
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -(2000+CTXT_)
            ELSE IF( ICTXT.NE.DESCX( CTXT_ ) ) THEN
               INFO = -(2400+CTXT_)
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -28
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -30
            END IF
            IDUM1( 1 ) = ICHAR( FACT )
            IDUM2( 1 ) = 1
            IDUM1( 2 ) = ICHAR( UPLO )
            IDUM2( 2 ) = 2
            IF( LSAME( FACT, 'F' ) ) THEN
               IDUM1( 3 ) = ICHAR( EQUED )
               IDUM2( 3 ) = 13
               IF( LWORK.EQ.-1 ) THEN
                  IDUM1( 4 ) = -1
               ELSE
                  IDUM1( 4 ) = 1
               END IF
               IDUM2( 4 ) = 28
               IF( LIWORK.EQ.-1 ) THEN
                  IDUM1( 5 ) = -1
               ELSE
                  IDUM1( 5 ) = 1
               END IF
               IDUM2( 5 ) = 30
               CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 8, N, 3, NRHS,
     $                        4, IB, JB, DESCB, 19, 5, IDUM1, IDUM2,
     $                        INFO )
            ELSE
               IF( LWORK.EQ.-1 ) THEN
                  IDUM1( 3 ) = -1
               ELSE
                  IDUM1( 3 ) = 1
               END IF
               IDUM2( 3 ) = 28
               IF( LIWORK.EQ.-1 ) THEN
                  IDUM1( 4 ) = -1
               ELSE
                  IDUM1( 4 ) = 1
               END IF
               IDUM2( 4 ) = 30
               CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 8, N, 3, NRHS,
     $                        4, IB, JB, DESCB, 19, 4, IDUM1, IDUM2,
     $                        INFO )
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PSPOSVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      IF( EQUIL ) THEN
*
*        Compute row and column scalings to equilibrate the matrix A.
*
         CALL PSPOEQU( N, A, IA, JA, DESCA, SR, SC, SCOND, AMAX,
     $                 INFEQU )
*
         IF( INFEQU.EQ.0 ) THEN
*
*           Equilibrate the matrix.
*
            CALL PSLAQSY( UPLO, N, A, IA, JA, DESCA, SR, SC, SCOND,
     $                    AMAX, EQUED )
            RCEQU = LSAME( EQUED, 'Y' )
         END IF
      END IF
*
*     Scale the right-hand side.
*
      CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL, IIB,
     $              JJB, IBROW, IBCOL )
      LDB =  DESCB( LLD_ )
      IROFF = MOD( IB-1, DESCB( MB_ ) )
      ICOFF = MOD( JB-1, DESCB( NB_ ) )
      NP = NUMROC( N+IROFF, DESCB( MB_ ), MYROW, IBROW, NPROW )
      NRHSQ = NUMROC( NRHS+ICOFF, DESCB( NB_ ), MYCOL, IBCOL, NPCOL )
      IF( MYROW.EQ.IBROW ) NP = NP-IROFF
      IF( MYCOL.EQ.IBCOL ) NRHSQ = NRHSQ-ICOFF
*
      IF( RCEQU ) THEN
         DO 30 J = JJB, JJB+NRHSQ-1
            DO 20 I = IIB, IIB+NP-1
               B( I + ( J-1 )*LDB ) = SR( I )*B( I + ( J-1 )*LDB )
   20       CONTINUE
   30    CONTINUE
      END IF
*
      IF( NOFACT .OR. EQUIL ) THEN
*
*        Compute the Cholesky factorization A = U'*U or A = L*L'.
*
         CALL PSLACPY( 'Full', N, N, A, IA, JA, DESCA, AF, IAF, JAF,
     $                 DESCAF )
         CALL PSPOTRF( UPLO, N, AF, IAF, JAF, DESCAF, INFO )
*
*        Return if INFO is non-zero.
*
         IF( INFO.NE.0 ) THEN
            IF( INFO.GT.0 )
     $         RCOND = ZERO
            RETURN
         END IF
      END IF
*
*     Compute the norm of the matrix A.
*
      ANORM = PSLANSY( '1', UPLO, N, A, IA, JA, DESCA, WORK )
*
*     Compute the reciprocal of the condition number of A.
*
      CALL PSPOCON( UPLO, N, AF, IAF, JAF, DESCAF, ANORM, RCOND, WORK,
     $              LWORK, IWORK, LIWORK, INFO )
*
*     Return if the matrix is singular to working precision.
*
      IF( RCOND.LT.PSLAMCH( ICTXT, 'Epsilon' ) ) THEN
         INFO = IA + N
         RETURN
      END IF
*
*     Compute the solution matrix X.
*
      CALL PSLACPY( 'Full', N, NRHS, B, IB, JB, DESCB, X, IX, JX,
     $              DESCX )
      CALL PSPOTRS( UPLO, N, NRHS, AF, IAF, JAF, DESCAF, X, IX, JX,
     $              DESCX, INFO )
*
*     Use iterative refinement to improve the computed solution and
*     compute error bounds and backward error estimates for it.
*
      CALL PSPORFS( UPLO, N, NRHS, A, IA, JA, DESCA, AF, IAF, JAF,
     $              DESCAF, B, IB, JB, DESCB, X, IX, JX, DESCX, FERR,
     $              BERR, WORK, LWORK, IWORK, LIWORK, INFO )
*
*     Transform the solution matrix X to a solution of the original
*     system.
*
      CALL INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL, IIX,
     $              JJX, IXROW, IXCOL )
      LDX = DESCX( LLD_ )
      IROFF = MOD( IX-1, DESCX( MB_ ) )
      ICOFF = MOD( JX-1, DESCX( NB_ ) )
      NP = NUMROC( N+IROFF, DESCX( MB_ ), MYROW, IXROW, NPROW )
      NRHSQ = NUMROC( NRHS+ICOFF, DESCX( NB_ ), MYCOL, IXCOL, NPCOL )
      IF( MYROW.EQ.IBROW ) NP = NP-IROFF
      IF( MYCOL.EQ.IBCOL ) NRHSQ = NRHSQ-ICOFF
*
      IF( RCEQU ) THEN
         DO 50 J = JJX, JJX+NRHSQ-1
            DO 40 I = IIX, IIX+NP-1
               X( I + ( J-1 )*LDX ) = SR( I )*X( I + ( J-1 )*LDX )
   40       CONTINUE
   50    CONTINUE
         DO 60 J = JJX, JJX+NRHSQ-1
            FERR( J ) = FERR( J ) / SCOND
   60    CONTINUE
      END IF
*
      WORK( 1 ) = REAL( LWMIN )
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of PSPOSVX
*
      END
