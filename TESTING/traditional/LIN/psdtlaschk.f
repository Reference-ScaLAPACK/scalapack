      SUBROUTINE PSDTLASCHK( SYMM, UPLO, TRANS, N, BWL, BWU, NRHS, X,
     $                       IX, JX, DESCX, IASEED, A, IA, JA, DESCA,
     $                       IBSEED, ANORM, RESID, WORK, WORKSIZ )
*
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          SYMM, TRANS, UPLO
      INTEGER            BWL, BWU, IA, IASEED, IBSEED,
     $                   IX, JA, JX, N, NRHS, WORKSIZ
      REAL               ANORM, RESID
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCX( * )
      REAL               A( * ), WORK( * ), X( * )
*     .. External Functions ..
      LOGICAL            LSAME
*     ..
*
*  Purpose
*  =======
*
*  PSDTLASCHK computes the residual
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
*          if SYMM = 'S', sub( A ) is a symmetric distributed band
*          matrix, otherwise sub( A ) is a general distributed matrix.
*
*  UPLO    (global input) CHARACTER
*          if SYMM = 'S', then
*            if UPLO = 'L', the lower half of the matrix is stored
*            if UPLO = 'U', the upper half of the matrix is stored
*          if SYMM != 'S' or 'H', then
*            if UPLO = 'D', the matrix is stable during factorization
*                           without interchanges
*            if UPLO != 'D', the matrix is general
*
*  TRANS   if TRANS= 'T', A 'Transpose' is used as the
*            coefficient matrix in the solve.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  NRHS    (global input) INTEGER
*          The number of right-hand-sides, i.e the number of columns
*          of the distributed matrix sub( X ). NRHS >= 1.
*
*  X       (local input) REAL pointer into the local memory
*          to an array of dimension (LLD_X,LOCq(JX+NRHS-1). This array
*          contains the local pieces of the answer vector(s) sub( X ) of
*          sub( A ) sub( X ) - B, split up over a column of processes.
*
*  IX      (global input) INTEGER
*          The row index in the global array X indicating the first
*          row of sub( X ).
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*
*  IASEED  (global input) INTEGER
*          The seed number to generate the original matrix Ao.
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
*  WORK    (local workspace) REAL array, dimension (LWORK)
*          IF SYMM='S'
*          LWORK >= max(5,NB)+2*NB
*          IF SYMM!='S' or 'H'
*          LWORK >= max(5,NB)+2*NB
*
*  WORKSIZ (local input) size of WORK.
*
*  =====================================================================
*
*  Code Developer: Andrew J. Cleary, University of Tennessee.
*    Current address: Lawrence Livermore National Labs.
*  This version released: August, 2001.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            INT_ONE
      PARAMETER          ( INT_ONE = 1 )
*     ..
*     .. Local Scalars ..
      INTEGER            IACOL, IAROW, ICTXT,
     $                   IIA, IIX, IPB, IPW,
     $                   IXCOL, IXROW, J, JJA, JJX, LDA,
     $                   MYCOL, MYROW, NB, NP, NPCOL, NPROW, NQ
      INTEGER            I, START
      INTEGER            BW, INFO, IPPRODUCT, WORK_MIN
      REAL               DIVISOR, EPS, RESID1, NORMX
*     ..
*     .. Local Arrays ..
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PBSTRAN, PSMATGEN,
     $                   SGAMX2D, SGEBR2D, SGEBS2D, SGEMM,
     $                   SGERV2D, SGESD2D, SGSUM2D, SLASET
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX, NUMROC
      REAL               PSLAMCH
      EXTERNAL           ISAMAX, NUMROC, PSLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD, REAL
*     ..
*     .. Executable Statements ..
*
*     Get needed initial parameters
*
      ICTXT = DESCA( CTXT_ )
      NB = DESCA( NB_ )
*
      IF( LSAME( SYMM, 'S' ) ) THEN
         BW = BWL
         START = 1
         WORK_MIN = MAX(5,NB)+2*NB
      ELSE
         BW = MAX(BWL, BWU)
         IF( LSAME( UPLO, 'D' )) THEN
            START = 1
         ELSE
            START = 2
         ENDIF
         WORK_MIN = MAX(5,NB)+2*NB
      ENDIF
*
      IF ( WORKSIZ .LT. WORK_MIN ) THEN
       CALL PXERBLA( ICTXT, 'PSTLASCHK', -18 )
       RETURN
      END IF
*
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
      NP = NUMROC( (3), DESCA( MB_ ), MYROW, 0, NPROW )
      NQ = NUMROC( N, DESCA( NB_ ), MYCOL, 0, NPCOL )
*
      IPB = 1
      IPPRODUCT = 1 + DESCA( NB_ )
      IPW = 1 + 2*DESCA( NB_ )
*
      LDA = DESCA( LLD_ )
*
*     Regenerate A
*
               IF( LSAME( SYMM, 'S' )) THEN
                 CALL PSBMATGEN( ICTXT, UPLO, 'D', BW, BW, N, BW+1,
     $                           DESCA( NB_ ), A, DESCA( LLD_ ), 0, 0,
     $                           IASEED, MYROW, MYCOL, NPROW, NPCOL )
               ELSE
*
                 CALL PSBMATGEN( ICTXT, 'N', UPLO, BWL, BWU, N,
     $                           DESCA( MB_ ), DESCA( NB_ ), A,
     $                           DESCA( LLD_ ), 0, 0, IASEED, MYROW,
     $                           MYCOL, NPROW, NPCOL )
               ENDIF
*
*           Matrix formed above has the diagonals shifted from what was
*              input to the tridiagonal routine. Shift them back.
*
*              Send elements to neighboring processors
*
               IF( MYCOL.GT.0 ) THEN
                  CALL SGESD2D( ICTXT, 1, 1, A( START+2), LDA,
     $                     MYROW, MYCOL-1 )
               ENDIF
*
               IF( MYCOL.LT.NPCOL-1 ) THEN
                  CALL SGESD2D( ICTXT, 1, 1,
     $                     A( START+( DESCA( NB_ )-1 )*LDA ),
     $                     LDA, MYROW, MYCOL+1 )
               ENDIF
*
*              Shift local elements
*
               DO 220 I=0,DESCA( NB_ )-1
                  A( START+2+(I)*LDA ) = A( START+2+(I+1)*LDA )
  220          CONTINUE
*
               DO 230 I=DESCA( NB_ )-1,0,-1
                  A( START+(I+1)*LDA ) = A( START+(I)*LDA )
  230          CONTINUE
*
*              Receive elements from neighboring processors
*
               IF( MYCOL.LT.NPCOL-1 ) THEN
                  CALL SGERV2D( ICTXT, 1, 1,
     $                     A( START+2+( DESCA( NB_ )-1 )*LDA ),
     $                     LDA, MYROW, MYCOL+1 )
               ENDIF
*
               IF( MYCOL.GT.0 ) THEN
                  CALL SGERV2D( ICTXT, 1, 1, A( START), LDA,
     $                               MYROW, MYCOL-1 )
               ENDIF
*
*     Loop over the rhs
*
      RESID = 0.0
*
      DO 40 J = 1, NRHS
*
*           Multiply A * current column of X
*
*
            CALL PSGBDCMV( BWL+BWU+1, BWL, BWU, TRANS, N, A, 1, DESCA,
     $                     1, X( 1 + (J-1)*DESCX( LLD_ )), 1, DESCX,
     $                     WORK( IPPRODUCT ), WORK( IPW ),
     $                     (INT_ONE+2)*INT_ONE, INFO )
*
*
*           Regenerate column of B
*
            CALL PSMATGEN( DESCX( CTXT_ ), 'No', 'No', DESCX( M_ ),
     $                     DESCX( N_ ), DESCX( MB_ ), DESCX( NB_ ),
     $                     WORK( IPB ), DESCX( LLD_ ), DESCX( RSRC_ ),
     $                     DESCX( CSRC_ ), IBSEED, 0, NQ, J-1, 1, MYCOL,
     $                     MYROW, NPCOL, NPROW )
*
*           Figure || A * X - B || & || X ||
*
            CALL PSAXPY( N, -ONE, WORK( IPPRODUCT ), 1, 1, DESCX, 1,
     $                   WORK( IPB ), 1, 1, DESCX, 1 )
*
            CALL PSNRM2( N, NORMX,
     $           X, 1, J, DESCX, 1 )
*
            CALL PSNRM2( N, RESID1,
     $           WORK( IPB ), 1, 1, DESCX, 1 )
*
*
*           Calculate residual = ||Ax-b|| / (||x||*||A||*eps*N)
*
            RESID1 = RESID1 / ( NORMX*DIVISOR )
*
            RESID = MAX( RESID, RESID1 )
*
   40 CONTINUE
*
      RETURN
*
*     End of PSTLASCHK
*
      END
