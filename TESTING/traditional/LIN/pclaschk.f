      SUBROUTINE PCLASCHK( SYMM, DIAG, N, NRHS, X, IX, JX, DESCX,
     $                     IASEED, IA, JA, DESCA, IBSEED, ANORM, RESID,
     $                     WORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, SYMM
      INTEGER            IA, IASEED, IBSEED, IX, JA, JX, N, NRHS
      REAL               ANORM, RESID
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCX( * )
      COMPLEX            WORK( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLASCHK computes the residual
*  || sub( A )*sub( X ) - B || / (|| sub( A ) ||*|| sub( X ) ||*eps*N)
*  to check the accuracy of the factorization and solve steps in the
*  LU and Cholesky decompositions, where sub( A ) denotes
*  A(IA:IA+N-1,JA,JA+N-1), sub( X ) denotes X(IX:IX+N-1, JX:JX+NRHS-1).
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
*  SYMM      (global input) CHARACTER
*          if SYMM = 'H', sub( A ) is a hermitian distributed matrix,
*          otherwise sub( A ) is a general distributed matrix.
*
*  DIAG    (global input) CHARACTER
*          If DIAG = 'D', sub( A ) is diagonally dominant.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right-hand-sides, i.e the number of columns
*          of the distributed matrix sub( X ). NRHS >= 0.
*
*  X       (local input) COMPLEX pointer into the local memory
*          to an array of dimension (LLD_X,LOCc(JX+NRHS-1). This array
*          contains the local pieces of the answer vector(s) sub( X ) of
*          sub( A ) sub( X ) - B, split up over a column of processes.
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
*  IASEED  (global input) INTEGER
*          The seed number to generate the original matrix Ao.
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
*  IBSEED  (global input) INTEGER
*          The seed number to generate the original matrix B.
*
*  ANORM   (global input) REAL
*          The 1-norm or infinity norm of the distributed matrix
*          sub( A ).
*
*  RESID   (global output) REAL
*          The residual error:
*          ||sub( A )*sub( X )-B|| / (||sub( A )||*||sub( X )||*eps*N).
*
*  WORK    (local workspace) COMPLEX array, dimension (LWORK)
*          LWORK >= MAX(1,Np)*NB_X + Nq*NB_X + MAX( MAX(NQ*MB_A,2*NB_X),
*          NB_X * NUMROC( NUMROC(N,MB_X,0,0,NPCOL), MB_X, 0, 0, LCMQ ) )
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
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                     ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            IACOL, IAROW, IB, ICOFF, ICTXT, ICURCOL, IDUMM,
     $                   II, IIA, IIX, IOFFX, IPA, IPB, IPW, IPX, IROFF,
     $                   IXCOL, IXROW, J, JBRHS, JJ, JJA, JJX, LDX,
     $                   MYCOL, MYROW, NP, NPCOL, NPROW, NQ
      REAL               DIVISOR, EPS, RESID1
      COMPLEX            BETA
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CGAMX2D, CGEMM, CGSUM2D,
     $                   CLASET, PBCTRAN, PCMATGEN, SGEBR2D,
     $                   SGEBS2D, SGERV2D, SGESD2D
*     ..
*     .. External Functions ..
      INTEGER            ICAMAX, NUMROC
      REAL               PSLAMCH
      EXTERNAL           ICAMAX, NUMROC, PSLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD, REAL
*     ..
*     .. Executable Statements ..
*
*     Get needed initial parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      EPS = PSLAMCH( ICTXT, 'eps' )
      RESID = 0.0E+0
      DIVISOR = ANORM * EPS * REAL( N )
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      CALL INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL, IIX, JJX,
     $              IXROW, IXCOL )
      IROFF = MOD( IA-1, DESCA( MB_ ) )
      ICOFF = MOD( JA-1, DESCA( NB_ ) )
      NP = NUMROC( N+IROFF, DESCA( MB_ ), MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFF, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
*
      LDX = MAX( 1, NP )
      IPB = 1
      IPX = IPB + NP * DESCX( NB_ )
      IPA = IPX + NQ * DESCX( NB_ )
*
      IF( MYROW.EQ.IAROW )
     $   NP = NP - IROFF
      IF( MYCOL.EQ.IACOL )
     $   NQ = NQ - ICOFF
*
      ICURCOL = IXCOL
*
*     Loop over the rhs
*
      DO 40 J = 1, NRHS, DESCX( NB_ )
         JBRHS = MIN( DESCX( NB_ ), NRHS-J+1 )
*
*        Transpose x from ICURCOL to all rows
*
         IOFFX = IIX + ( JJX - 1 ) * DESCX( LLD_ )
         CALL PBCTRAN( ICTXT, 'Column', 'Transpose', N, JBRHS,
     $              DESCX( MB_ ), X( IOFFX ), DESCX( LLD_ ), ZERO,
     $              WORK( IPX ), JBRHS, IXROW, ICURCOL, -1, IACOL,
     $              WORK( IPA ) )
*
*        Regenerate B in IXCOL
*
         IF( MYCOL.EQ.ICURCOL ) THEN
            CALL PCMATGEN( ICTXT, 'N', 'N', DESCX( M_ ), DESCX( N_ ),
     $                     DESCX( MB_ ), DESCX( NB_ ), WORK( IPB ), LDX,
     $                     IXROW, IXCOL, IBSEED, IIX-1, NP, JJX-1,
     $                     JBRHS, MYROW, MYCOL, NPROW, NPCOL )
            BETA = ONE
         ELSE
            BETA = ZERO
         END IF
*
         IF( NQ.GT.0 ) THEN
            DO 10 II = IIA, IIA+NP-1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IIA+NP-II )
*
*              Regenerate ib rows of the matrix A(IA:IA+N-1,JA:JA+N-1).
*
               CALL PCMATGEN( ICTXT, SYMM, DIAG, DESCA( M_ ),
     $                        DESCA( N_ ), DESCA( MB_ ), DESCA( NB_ ),
     $                        WORK( IPA ), IB, DESCA( RSRC_ ),
     $                        DESCA( CSRC_ ), IASEED, II-1, IB,
     $                        JJA-1, NQ, MYROW, MYCOL, NPROW, NPCOL )
*
*              Compute B <= B - A * X.
*
               CALL CGEMM( 'No transpose', 'Transpose', IB, JBRHS, NQ,
     $                     -ONE, WORK( IPA ), IB, WORK( IPX ), JBRHS,
     $                     BETA, WORK( IPB+II-IIA ), LDX )
*
   10       CONTINUE
*
         ELSE IF( MYCOL.NE.ICURCOL ) THEN
*
            CALL CLASET( 'All', NP, JBRHS, ZERO, ZERO, WORK( IPB ),
     $                   LDX )
*
         END IF
*
*        Add B rowwise to ICURCOL
*
         CALL CGSUM2D( ICTXT, 'Row', ' ', NP, JBRHS, WORK( IPB ), LDX,
     $                 MYROW, ICURCOL )
*
         IF( MYCOL.EQ.ICURCOL ) THEN
*
*           Figure || A * X - B || & || X ||
*
            IPW = IPA + JBRHS
            DO 20 JJ = 0, JBRHS - 1
               IF( NP.GT.0 ) THEN
                  II = ICAMAX( NP, WORK( IPB+JJ*LDX ), 1 )
                  WORK( IPA+JJ ) = ABS( WORK( IPB+II-1+JJ*LDX ) )
                  WORK( IPW+JJ ) = ABS( X( IOFFX + ICAMAX( NP,
     $            X( IOFFX + JJ*DESCX( LLD_ ) ), 1 )-1+JJ*
     $            DESCX( LLD_ ) ) )
               ELSE
                  WORK( IPA+JJ ) = ZERO
                  WORK( IPW+JJ ) = ZERO
               END IF
   20       CONTINUE
*
*           After CGAMX2D computation,
*              WORK(IPB) has the maximum of || Ax - b ||, and
*              WORK(IPX) has the maximum of || X ||.
*
            CALL CGAMX2D( ICTXT, 'Column', ' ', 1, 2*JBRHS,
     $                    WORK( IPA ), 1, IDUMM, IDUMM, -1, 0, ICURCOL )
*
*           Calculate residual = ||Ax-b|| / (||x||*||A||*eps*N)
*
            IF( MYROW.EQ.0 ) THEN
               DO 30 JJ = 0, JBRHS - 1
                  RESID1 = REAL( WORK( IPA+JJ ) ) /
     $                     ( REAL( WORK( IPW+JJ ) )*DIVISOR )
                  IF( RESID.LT.RESID1 )
     $               RESID = RESID1
   30          CONTINUE
               IF( MYCOL.NE.0 )
     $            CALL SGESD2D( ICTXT, 1, 1, RESID, 1, 0, 0 )
            END IF
*
         ELSE IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
*
            CALL SGERV2D( ICTXT, 1, 1, RESID1, 1, 0, ICURCOL )
            IF( RESID.LT.RESID1 )
     $         RESID = RESID1
*
         END IF
*
         IF( MYCOL.EQ.ICURCOL )
     $      JJX = JJX + JBRHS
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
   40 CONTINUE
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         CALL SGEBS2D( ICTXT, 'All', ' ', 1, 1, RESID, 1 )
      ELSE
         CALL SGEBR2D( ICTXT, 'All', ' ', 1, 1, RESID, 1, 0, 0 )
      END IF
*
      RETURN
*
*     End of PCLASCHK
*
      END
