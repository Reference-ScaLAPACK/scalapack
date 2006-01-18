      SUBROUTINE PCMAX1( N, AMAX, INDX, X, IX, JX, DESCX, INCX )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
      INTEGER            INDX, INCX, IX, JX, N
      COMPLEX            AMAX
*     ..
*     .. Array Arguments ..
      INTEGER            DESCX( * )
      COMPLEX            X( * )
*     ..
*
*  Purpose
*  =======
*
*  PCMAX1 computes the global index of the maximum element in absolute
*  value of a distributed vector sub( X ). The global index is returned
*  in INDX and the value is returned in AMAX,
*
*  where sub( X ) denotes X(IX:IX+N-1,JX) if INCX = 1,
*                         X(IX,JX:JX+N-1) if INCX = M_X.
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
*  Based on PCAMAX from Level 1 PBLAS. The change is to use the
*  'genuine' absolute value.
*
*  The serial version was contributed to LAPACK by Nick Higham for use
*  with CLACON.
*
*  Arguments
*  =========
*
*  N       (global input) pointer to INTEGER
*          The number of components of the distributed vector sub( X ).
*          N >= 0.
*
*  AMAX    (global output) pointer to REAL
*          The absolute value of the largest entry of the distributed
*          vector sub( X ) only in the scope of sub( X ).
*
*  INDX    (global output) pointer to INTEGER
*          The global index of the element of the distributed vector
*          sub( X ) whose real part has maximum absolute value.
*
*  X       (local input) COMPLEX array containing the local
*          pieces of a distributed matrix of dimension of at least
*              ( (JX-1)*M_X + IX + ( N - 1 )*abs( INCX ) )
*          This array contains the entries of the distributed vector
*          sub( X ).
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
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          CBTOP, CCTOP, RBTOP, RCTOP
      INTEGER            ICOFF, ICTXT, IDUMM, IIX, IROFF, IXCOL, IXROW,
     $                   JJX, LCINDX, LDX, MAXPOS, MYCOL, MYROW, NP,
     $                   NPCOL, NPROW, NQ
*     ..
*     .. Local Arrays ..
      COMPLEX            WORK( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CCOMBAMAX1, CGAMX2D,
     $                   IGEBR2D, IGEBS2D, INFOG2L, PCTREECOMB,
     $                   PB_TOPGET
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICMAX1, INDXL2G, NUMROC
      EXTERNAL           ICMAX1, INDXL2G, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, MOD, NINT, REAL
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCX( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Quick return if possible.
*
      INDX = 0
      AMAX = ZERO
      IF( N.LE.0 )
     $   RETURN
*
*     Retrieve local information for vector X.
*
      LDX = DESCX( LLD_ )
      CALL INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL, IIX, JJX,
     $              IXROW, IXCOL )
*
      IF( INCX.EQ.1 .AND. DESCX( M_ ).EQ.1 .AND. N.EQ.1 ) THEN
         INDX = JX
         AMAX = X( IIX+(JJX-1)*LDX )
         RETURN
      END IF
*
*     Find the maximum value and its index
*
      IF( INCX.EQ.DESCX( M_ ) ) THEN
*
         IF( MYROW.EQ.IXROW ) THEN
*
            ICOFF = MOD( JX-1, DESCX( NB_ ) )
            NQ = NUMROC( N+ICOFF, DESCX( NB_ ), MYCOL, IXCOL, NPCOL )
            IF( MYCOL.EQ.IXCOL )
     $         NQ = NQ-ICOFF
*
            CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', RBTOP )
*
            IF( LSAME( RBTOP, ' ' ) ) THEN
*
               IF( NQ.GT.0 ) THEN
                  LCINDX = JJX-1+ICMAX1( NQ, X( IIX+(JJX-1)*LDX ), LDX )
                  WORK( 1 ) = X( IIX+(LCINDX-1)*LDX )
                  WORK( 2 ) = CMPLX( REAL( INDXL2G( LCINDX,
     $              DESCX( NB_ ), MYCOL, DESCX( CSRC_ ), NPCOL ) ) )
               ELSE
                  WORK( 1 ) = ZERO
                  WORK( 2 ) = ZERO
               END IF
*
               CALL PCTREECOMB( ICTXT, 'Row', 2, WORK, -1, MYCOL,
     $                          CCOMBAMAX1 )
*
               AMAX = WORK( 1 )
               IF( AMAX.EQ.ZERO ) THEN
                  INDX = JX
               ELSE
                  INDX = NINT( REAL( WORK( 2 ) ) )
               END IF
*
            ELSE
*
               CALL PB_TOPGET( ICTXT, 'Combine', 'Rowwise', RCTOP )
*
               IF( NQ.GT.0 ) THEN
                  LCINDX = JJX-1+ICMAX1( NQ, X( IIX+(JJX-1)*LDX ), LDX )
                  AMAX = X( IIX + (LCINDX-1)*LDX )
               ELSE
                  AMAX = ZERO
               END IF
*
*              Find the maximum value
*
               CALL CGAMX2D( ICTXT, 'Rowwise', RCTOP, 1, 1, AMAX, 1,
     $                       IDUMM, MAXPOS, 1, -1, MYROW )
*
               IF( AMAX.NE.ZERO ) THEN
*
*                 Broadcast corresponding global index
*
                  IF( MYCOL.EQ.MAXPOS ) THEN
                     INDX = INDXL2G( LCINDX, DESCX( NB_ ), MYCOL,
     $                               DESCX( CSRC_ ), NPCOL )
                     CALL IGEBS2D( ICTXT, 'Rowwise', RBTOP, 1, 1, INDX,
     $                             1 )
                  ELSE
                     CALL IGEBR2D( ICTXT, 'Rowwise', RBTOP, 1, 1, INDX,
     $                             1, MYROW, MAXPOS )
                  END IF
*
               ELSE
*
                  INDX = JX
*
               END IF
*
            END IF
*
         END IF
*
      ELSE
*
         IF( MYCOL.EQ.IXCOL ) THEN
*
            IROFF = MOD( IX-1, DESCX( MB_ ) )
            NP = NUMROC( N+IROFF, DESCX( MB_ ), MYROW, IXROW, NPROW )
            IF( MYROW.EQ.IXROW )
     $         NP = NP-IROFF
*
            CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', CBTOP )
*
            IF( LSAME( CBTOP, ' ' ) ) THEN
*
               IF( NP.GT.0 ) THEN
                  LCINDX = IIX-1+ICMAX1( NP, X( IIX+(JJX-1)*LDX ), 1 )
                  WORK( 1 ) = X( LCINDX + (JJX-1)*LDX )
                  WORK( 2 ) = CMPLX( REAL( INDXL2G( LCINDX,
     $              DESCX( MB_ ), MYROW, DESCX( RSRC_ ), NPROW ) ) )
               ELSE
                  WORK( 1 ) = ZERO
                  WORK( 2 ) = ZERO
               END IF
*
               CALL PCTREECOMB( ICTXT, 'Column', 2, WORK, -1, MYCOL,
     $                          CCOMBAMAX1 )
*
               AMAX = WORK( 1 )
               IF( AMAX.EQ.ZERO ) THEN
                  INDX = IX
               ELSE
                  INDX = NINT( REAL( WORK( 2 ) ) )
               END IF
*
            ELSE
*
               CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', CCTOP )
*
               IF( NP.GT.0 ) THEN
                  LCINDX = IIX-1+ICMAX1( NP, X( IIX+(JJX-1)*LDX ), 1 )
                  AMAX = X( LCINDX + (JJX-1)*LDX )
               ELSE
                  AMAX = ZERO
               END IF
*
*              Find the maximum value
*
               CALL CGAMX2D( ICTXT, 'Columnwise', CCTOP, 1, 1, AMAX, 1,
     $                       MAXPOS, IDUMM, 1, -1, MYCOL )
*
               IF( AMAX.NE.ZERO ) THEN
*
*                 Broadcast corresponding global index
*
                  IF( MYROW.EQ.MAXPOS ) THEN
                     INDX = INDXL2G( LCINDX, DESCX( MB_ ), MYROW,
     $                               DESCX( RSRC_ ), NPROW )
                     CALL IGEBS2D( ICTXT, 'Columnwise', CBTOP, 1, 1,
     $                             INDX, 1 )
                  ELSE
                     CALL IGEBR2D( ICTXT, 'Columnwise', CBTOP, 1, 1,
     $                             INDX, 1, MAXPOS, MYCOL )
                  END IF
*
               ELSE
*
                  INDX = IX
*
               END IF
*
            END IF
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PCMAX1
*
      END
*
      SUBROUTINE CCOMBAMAX1 ( V1, V2 )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Array Arguments ..
      COMPLEX            V1( 2 ), V2( 2 )
*     ..
*
*  Purpose
*  =======
*
*  CCOMBAMAX1 finds the element having maximum real part absolute
*  value as well as its corresponding globl index.
*
*  Arguments
*  =========
*
*  V1        (local input/local output) COMPLEX array of
*            dimension 2.  The first maximum absolute value element and
*            its global index. V1(1) = AMAX, V1(2) = INDX.
*
*  V2        (local input) COMPLEX array of dimension 2.
*            The second maximum absolute value element and its global
*            index. V2(1) = AMAX, V2(2) = INDX.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL
*     ..
*     .. Executable Statements ..
*
      IF( ABS( REAL( V1( 1 ) ) ).LT.ABS( REAL( V2( 1 ) ) ) ) THEN
         V1( 1 ) = V2( 1 )
         V1( 2 ) = V2( 2 )
      END IF
*
      RETURN
*
*     End of CCOMBAMAX1
*
      END
