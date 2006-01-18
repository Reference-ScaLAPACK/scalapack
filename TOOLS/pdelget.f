      SUBROUTINE PDELGET( SCOPE, TOP, ALPHA, A, IA, JA, DESCA )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER*1        SCOPE, TOP
      INTEGER            IA, JA
      DOUBLE PRECISION   ALPHA
*     ..
*     .. Array arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * )
*     ..
*
*  Purpose
*  =======
*
*  PDELGET sets alpha to the distributed matrix entry A( IA, JA ).
*  The value of alpha is set according to the scope.
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
*  SCOPE   (global input) CHARACTER*1
*          The BLACS scope in which alpha is updated.
*          If SCOPE = 'R', alpha is updated only in the process row
*                          containing A( IA, JA ),
*          If SCOPE = 'C', alpha is updated only in the process column
*                          containing A( IA, JA ),
*          If SCOPE = 'A', alpha is updated in all the processes of the
*                          grid,
*          otherwise alpha is updated only in the process containing
*           A( IA, JA ).
*
*  TOP     (global input) CHARACTER*1
*          The topology to be used if broadcast is needed.
*
*  ALPHA   (global output) DOUBLE PRECISION, the scalar alpha.
*
*  A       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_A,*) containing the local
*          pieces of the distributed matrix A.
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
      INTEGER            IACOL, IAROW, ICTXT, IIA, IOFFA, JJA, MYCOL,
     $                   MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGEBR2D, DGEBS2D, INFOG2L
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters.
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
*
      ALPHA = ZERO
*
      IF( LSAME( SCOPE, 'R' ) ) THEN
         IF( MYROW.EQ.IAROW ) THEN
            IF( MYCOL.EQ.IACOL ) THEN
               IOFFA = IIA+(JJA-1)*DESCA( LLD_ )
               CALL DGEBS2D( ICTXT, SCOPE, TOP, 1, 1, A( IOFFA ), 1 )
               ALPHA = A( IOFFA )
            ELSE
               CALL DGEBR2D( ICTXT, SCOPE, TOP, 1, 1, ALPHA, 1,
     $                       IAROW, IACOL )
            END IF
         END IF
      ELSE IF( LSAME( SCOPE, 'C' ) ) THEN
         IF( MYCOL.EQ.IACOL ) THEN
            IF( MYROW.EQ.IAROW ) THEN
               IOFFA = IIA+(JJA-1)*DESCA( LLD_ )
               CALL DGEBS2D( ICTXT, SCOPE, TOP, 1, 1, A( IOFFA ), 1 )
               ALPHA = A( IOFFA )
            ELSE
               CALL DGEBR2D( ICTXT, SCOPE, TOP, 1, 1, ALPHA, 1,
     $                       IAROW, IACOL )
            END IF
         END IF
      ELSE IF( LSAME( SCOPE, 'A' ) ) THEN
         IF( ( MYROW.EQ.IAROW ).AND.( MYCOL.EQ.IACOL ) ) THEN
            IOFFA = IIA+(JJA-1)*DESCA( LLD_ )
            CALL DGEBS2D( ICTXT, SCOPE, TOP, 1, 1, A( IOFFA ), 1 )
            ALPHA = A( IOFFA )
         ELSE
            CALL DGEBR2D( ICTXT, SCOPE, TOP, 1, 1, ALPHA, 1,
     $                    IAROW, IACOL )
         END IF
      ELSE
         IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL )
     $      ALPHA = A( IIA+(JJA-1)*DESCA( LLD_ ) )
      END IF
*
      RETURN
*
*     End of PDELGET
*
      END
