      SUBROUTINE PZLACGV( N, X, IX, JX, DESCX, INCX )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INCX, IX, JX, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCX( * )
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLACGV conjugates a complex vector of length N, sub( X ), where
*  sub( X ) denotes X(IX,JX:JX+N-1) if INCX = DESCX( M_ ) and
*  X(IX:IX+N-1,JX) if INCX = 1, and
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
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The length of the distributed vector sub( X ).
*
*  X       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_X,*).
*          On entry the vector to be conjugated
*             x( i )  = X(IX+(JX-1)*M_X +(i-1)*INCX ), 1 <= i <= N.
*          On exit the conjugated vector.
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
*  INCX    (global input) INTEGER
*          The global increment for the elements of X. Only two values
*          of INCX are supported in this version, namely 1 and M_X.
*          INCX must not be zero.
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
      INTEGER            I, ICOFFX, ICTXT, IIX, IOFFX, IROFFX, IXCOL,
     $                   IXROW, JJX, LDX, MYCOL, MYROW, NP, NPCOL,
     $                   NPROW, NQ
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters.
*
      ICTXT = DESCX( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Figure local indexes
*
      CALL INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL,
     $              IIX, JJX, IXROW, IXCOL )
*
      LDX = DESCX( LLD_ )
      IF( INCX.EQ.DESCX( M_ ) ) THEN
*
*        sub( X ) is rowwise distributed.
*
         IF( MYROW.NE.IXROW )
     $      RETURN
         ICOFFX = MOD( JX-1, DESCX( NB_ ) )
         NQ = NUMROC( N+ICOFFX, DESCX( NB_ ), MYCOL, IXCOL, NPCOL )
         IF( MYCOL.EQ.IXCOL )
     $      NQ = NQ - ICOFFX
*
         IF( NQ.GT.0 ) THEN
            IOFFX = IIX+(JJX-1)*LDX
            DO 10 I = 1, NQ
               X( IOFFX ) = DCONJG( X( IOFFX ) )
               IOFFX = IOFFX + LDX
   10       CONTINUE
         END IF
*
      ELSE IF( INCX.EQ.1 ) THEN
*
*        sub( X ) is columnwise distributed.
*
         IF( MYCOL.NE.IXCOL )
     $      RETURN
         IROFFX = MOD( IX-1, DESCX( MB_ ) )
         NP = NUMROC( N+IROFFX, DESCX( MB_ ), MYROW, IXROW, NPROW )
         IF( MYROW.EQ.IXROW )
     $      NP = NP - IROFFX
*
         IF( NP.GT.0 ) THEN
            IOFFX = IIX+(JJX-1)*LDX
            DO 20 I = IOFFX, IOFFX+NP-1
              X( I ) = DCONJG( X( I ) )
   20       CONTINUE
         END IF
*
      END IF
*
      RETURN
*
*     End of PZLACGV
*
      END
