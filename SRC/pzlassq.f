      SUBROUTINE PZLASSQ( N, X, IX, JX, DESCX, INCX, SCALE, SUMSQ )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IX, INCX, JX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      INTEGER            DESCX( * )
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where x( i ) = sub( X ) = abs( X( IX+(JX-1)*DESCX(M_)+(i-1)*INCX ) ).
*  The value of sumsq is assumed to be at least unity and the value of
*  ssq will then satisfy
*
*     1.0 .le. ssq .le. ( sumsq + 2*n ).
*
*  scale is assumed to be non-negative and scl returns the value
*
*     scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
*            i
*
*  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
*  SCALE and SUMSQ are overwritten by scl and ssq respectively.
*
*  The routine makes only one pass through the vector sub( X ).
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
*  The result are only available in the scope of sub( X ), i.e if
*  sub( X ) is distributed along a process row, the correct results are
*  only available in this process row of the grid. Similarly if sub( X )
*  is distributed along a process column, the correct results are only
*  available in this process column of the grid.
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The length of the distributed vector sub( X ).
*
*  X       (input) COMPLEX*16
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X(IX+(JX-1)*M_X +(i-1)*INCX ), 1 <= i <= n.
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
*  SCALE   (local input/local output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (local input/local output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
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
      INTEGER            I, ICOFF, ICTXT, IIX, IOFF, IROFF, IXCOL,
     $                   IXROW, JJX, LDX, MYCOL, MYROW, NP, NPCOL,
     $                   NPROW, NQ
      DOUBLE PRECISION   TEMP1
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   WORK( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DCOMBSSQ, INFOG2L, PDTREECOMB
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MOD
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
      CALL INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL, IIX, JJX,
     $              IXROW, IXCOL )
*
      LDX = DESCX( LLD_ )
      IF( INCX.EQ.DESCX( M_ ) ) THEN
*
*        X is rowwise distributed.
*
         IF( MYROW.NE.IXROW )
     $      RETURN
         ICOFF = MOD( JX, DESCX( NB_ ) )
         NQ = NUMROC( N+ICOFF, DESCX( NB_ ), MYCOL, IXCOL, NPCOL )
         IF( MYCOL.EQ.IXCOL )
     $      NQ = NQ - ICOFF
*
*        Code direct from LAPACK's ZLASSQ, (save subroutine call)
*
         IF( NQ.GT.0 ) THEN
            IOFF = IIX + ( JJX - 1 ) * LDX
            DO 10 I = 1, NQ
               IF( DBLE( X( IOFF ) ).NE.ZERO ) THEN
                  TEMP1 = ABS( DBLE( X( IOFF ) ) )
                  IF( SCALE.LT.TEMP1 ) THEN
                     SUMSQ = 1 + SUMSQ * ( SCALE / TEMP1 )**2
                     SCALE = TEMP1
                  ELSE
                     SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
                  END IF
               END IF
               IF( DIMAG( X( IOFF ) ).NE.ZERO ) THEN
                  TEMP1 = ABS( DIMAG( X( IOFF ) ) )
                  IF( SCALE.LT.TEMP1 ) THEN
                     SUMSQ = 1 + SUMSQ * ( SCALE / TEMP1 )**2
                     SCALE = TEMP1
                  ELSE
                     SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
                  END IF
               END IF
               IOFF = IOFF + LDX
   10       CONTINUE
         END IF
*
*        Take local result and find global
*
         WORK( 1 ) = SCALE
         WORK( 2 ) = SUMSQ
*
         CALL PDTREECOMB( ICTXT, 'Rowwise', 2, WORK, -1, IXCOL,
     $                    DCOMBSSQ )
*
         SCALE = WORK( 1 )
         SUMSQ = WORK( 2 )
*
      ELSE IF( INCX.EQ.1 ) THEN
*
*        X is columnwise distributed.
*
         IF( MYCOL.NE.IXCOL )
     $      RETURN
         IROFF = MOD( IX, DESCX( MB_ ) )
         NP = NUMROC( N+IROFF, DESCX( MB_ ), MYROW, IXROW, NPROW )
         IF( MYROW.EQ.IXROW )
     $      NP = NP - IROFF
*
*        Code direct from LAPACK's ZLASSQ, (save subroutine call)
*
         IF( NP.GT.0 ) THEN
            IOFF = IIX + ( JJX - 1 ) * LDX
            DO 20 I = 1, NP
               IF( DBLE( X( IOFF ) ).NE.ZERO ) THEN
                  TEMP1 = ABS( DBLE( X( IOFF ) ) )
                  IF( SCALE.LT.TEMP1 ) THEN
                     SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                     SCALE = TEMP1
                  ELSE
                     SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
                  END IF
               END IF
               IF( DIMAG( X( IOFF ) ).NE.ZERO ) THEN
                  TEMP1 = ABS( DIMAG( X( IOFF ) ) )
                  IF( SCALE.LT.TEMP1 ) THEN
                     SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                     SCALE = TEMP1
                  ELSE
                     SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
                  END IF
               END IF
               IOFF = IOFF + 1
   20       CONTINUE
         END IF
*
*        Take local result and find global
*
         WORK( 1 ) = SCALE
         WORK( 2 ) = SUMSQ
*
         CALL PDTREECOMB( ICTXT, 'Columnwise', 2, WORK, -1, IXCOL,
     $                    DCOMBSSQ )
*
         SCALE = WORK( 1 )
         SUMSQ = WORK( 2 )
*
      END IF
*
      RETURN
*
*     End of PZLASSQ
*
      END
