      SUBROUTINE PDZSUM1( N, ASUM, X, IX, JX, DESCX, INCX )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER        IX, INCX, JX, N
      DOUBLE PRECISION   ASUM
*     ..
*     .. Array Arguments ..
      INTEGER       DESCX( * )
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  PDZSUM1 returns the sum of absolute values of a complex
*  distributed vector sub( X ) in ASUM,
*
*  where sub( X ) denotes X(IX:IX+N-1,JX:JX), if INCX = 1,
*                         X(IX:IX,JX:JX+N-1), if INCX = M_X.
*
*  Based on PDZASUM from the Level 1 PBLAS. The change is
*  to use the 'genuine' absolute value.
*
*  The serial version of this routine was originally contributed by
*  Nick Higham for use with ZLACON.
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
*  Because vectors may be viewed as a subclass of matrices, a
*  distributed vector is considered to be a distributed matrix.
*
*  When the result of a vector-oriented PBLAS call is a scalar, it will
*  be made available only within the scope which owns the vector(s)
*  being operated on.  Let X be a generic term for the input vector(s).
*  Then, the processes which receive the answer will be (note that if
*  an operation involves more than one vector, the processes which re-
*  ceive the result will be the union of the following calculation for
*  each vector):
*
*  If N = 1, M_X = 1 and INCX = 1, then one can't determine if a process
*  row or process column owns the vector operand, therefore only the
*  process of coordinate {RSRC_X, CSRC_X} receives the result;
*
*  If INCX = M_X, then sub( X ) is a vector distributed over a process
*  row. Each process part of this row receives the result;
*
*  If INCX = 1, then sub( X ) is a vector distributed over a process
*  column. Each process part of this column receives the result;
*
*  Parameters
*  ==========
*
*  N       (global input) pointer to INTEGER
*          The number of components of the distributed vector sub( X ).
*          N >= 0.
*
*  ASUM    (local output) pointer to DOUBLE PRECISION
*          The sum of absolute values of the distributed vector sub( X )
*          only in its scope.
*
*  X       (local input) COMPLEX*16 array containing the local
*          pieces of a distributed matrix of dimension of at least
*              ( (JX-1)*M_X + IX + ( N - 1 )*abs( INCX ) )
*          This array contains the entries of the distributed vector
*          sub( X ).
*
*  IX      (global input) pointer to INTEGER
*          The global row index of the submatrix of the distributed
*          matrix X to operate on.
*
*  JX      (global input) pointer to INTEGER
*          The global column index of the submatrix of the distributed
*          matrix X to operate on.
*
*  DESCX   (global and local input) INTEGER array of dimension 8.
*          The array descriptor of the distributed matrix X.
*
*  INCX    (global input) pointer to INTEGER
*          The global increment for the elements of X. Only two values
*          of INCX are supported in this version, namely 1 and M_X.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          CCTOP, RCTOP
      INTEGER            ICOFF, ICTXT, IIX, IROFF, IXCOL, IXROW, JJX,
     $                   LDX, MYCOL, MYROW, NP, NPCOL, NPROW, NQ
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGSUM2D, INFOG2L, PB_TOPGET
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      DOUBLE PRECISION   DZSUM1
      EXTERNAL           DZSUM1, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MOD
*     ..
*     .. Executable Statements ..
*
      ICTXT = DESCX( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Quick return if possible
*
      ASUM = ZERO
      IF( N.LE.0 )
     $   RETURN
*
      LDX = DESCX( LLD_ )
      CALL INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL, IIX, JJX,
     $              IXROW, IXCOL )
*
      IF( INCX.EQ.1 .AND. DESCX( M_ ).EQ.1 .AND. N.EQ.1 ) THEN
         IF( MYROW.EQ.IXROW .AND. MYCOL.EQ.IXCOL ) THEN
            ASUM = ABS( X( IIX+(JJX-1)*LDX ) )
         END IF
         RETURN
      END IF
*
      IF( INCX.EQ.DESCX( M_ ) ) THEN
*
*        X is distributed over a process row
*
         IF( MYROW.EQ.IXROW ) THEN
            CALL PB_TOPGET( ICTXT, 'Combine', 'Rowwise', RCTOP )
            ICOFF = MOD( JX-1, DESCX( NB_ ) )
            NQ = NUMROC( N+ICOFF, DESCX( NB_ ), MYCOL, IXCOL, NPCOL )
            IF( MYCOL.EQ.IXCOL )
     $         NQ = NQ-ICOFF
            ASUM = DZSUM1( NQ, X( IIX+(JJX-1)*LDX ), LDX )
            CALL DGSUM2D( ICTXT, 'Rowwise', RCTOP, 1, 1, ASUM, 1,
     $                    -1, MYCOL )
         END IF
*
      ELSE
*
*        X is distributed over a process column
*
         IF( MYCOL.EQ.IXCOL ) THEN
            CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', CCTOP )
            IROFF = MOD( IX-1, DESCX( MB_ ) )
            NP = NUMROC( N+IROFF, DESCX( MB_ ), MYROW, IXROW, NPROW )
            IF( MYROW.EQ.IXROW )
     $         NP = NP-IROFF
            ASUM = DZSUM1( NP, X( IIX+(JJX-1)*LDX ), 1 )
            CALL DGSUM2D( ICTXT, 'Columnwise', CCTOP, 1, 1, ASUM, 1,
     $                    -1, MYCOL )
         END IF
*
      END IF
*
      RETURN
*
*     End of PDZSUM1
*
      END
