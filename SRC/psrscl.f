      SUBROUTINE PSRSCL( N, SA, SX, IX, JX, DESCX, INCX )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IX, INCX, JX, N
      REAL               SA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCX( * )
      REAL               SX( * )
*     ..
*
*  Purpose
*  =======
*
*  PSRSCL multiplies an N-element real distributed vector sub( X ) by
*  the real scalar 1/a. This is done without overflow or underflow as
*  long as the final result sub( X )/a does not overflow or underflow.
*
*  where sub( X ) denotes X(IX:IX+N-1,JX:JX), if INCX = 1,
*                         X(IX:IX,JX:JX+N-1), if INCX = M_X.
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
*  Such a global array has an associated description vector descA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DT_A   (global) descA[ DT_ ]   The descriptor type.  In this case,
*                                 DT_A = 1.
*  CTXT_A (global) descA[ CTXT_ ] The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) descA[ M_ ]    The number of rows in the global
*                                 array A.
*  N_A    (global) descA[ N_ ]    The number of columns in the global
*                                 array A.
*  MB_A   (global) descA[ MB_ ]   The blocking factor used to distribu-
*                                 te the rows of the array.
*  NB_A   (global) descA[ NB_ ]   The blocking factor used to distribu-
*                                 te the columns of the array.
*  RSRC_A (global) descA[ RSRC_ ] The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) descA[ CSRC_ ] The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  descA[ LLD_ ]  The leading dimension of the local
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
*  Because vectors may be seen as particular matrices, a distributed
*  vector is considered to be a distributed matrix.
*
*  Arguments
*  =========
*
*  N       (global input) pointer to INTEGER
*          The number of components of the distributed vector sub( X ).
*          N >= 0.
*
*  SA      (global input) REAL
*          The scalar a which is used to divide each component of
*          sub( X ).  SA must be >= 0, or the subroutine will divide by
*          zero.
*
*  SX      (local input/local output) REAL array
*          containing the local pieces of a distributed matrix of
*          dimension of at least
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
* =====================================================================
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
      LOGICAL            DONE
      INTEGER            ICTXT, MYCOL, MYROW, NPCOL, NPROW
      REAL               BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PSLABAD, PSSCAL
*     ..
*     .. External Functions ..
      REAL               PSLAMCH
      EXTERNAL           PSLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCX( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = PSLAMCH( ICTXT, 'S' )
      BIGNUM = ONE / SMLNUM
      CALL PSLABAD( ICTXT, SMLNUM, BIGNUM )
*
*     Initialize the denominator to SA and the numerator to 1.
*
      CDEN = SA
      CNUM = ONE
*
   10 CONTINUE
      CDEN1 = CDEN*SMLNUM
      CNUM1 = CNUM / BIGNUM
      IF( ABS( CDEN1 ).GT.ABS( CNUM ) .AND. CNUM.NE.ZERO ) THEN
*
*        Pre-multiply sub( X ) by SMLNUM if CDEN is large compared to
*        CNUM.
*
         MUL = SMLNUM
         DONE = .FALSE.
         CDEN = CDEN1
      ELSE IF( ABS( CNUM1 ).GT.ABS( CDEN ) ) THEN
*
*        Pre-multiply sub( X ) by BIGNUM if CDEN is small compared to
*        CNUM.
*
         MUL = BIGNUM
         DONE = .FALSE.
         CNUM = CNUM1
      ELSE
*
*        Multiply sub( X ) by CNUM / CDEN and return.
*
         MUL = CNUM / CDEN
         DONE = .TRUE.
      END IF
*
*     Scale the vector sub( X ) by MUL
*
      CALL PSSCAL( N, MUL, SX, IX, JX, DESCX, INCX )
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of PSRSCL
*
      END
