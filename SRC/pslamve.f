      SUBROUTINE PSLAMVE( UPLO, M, N, A, IA, JA, DESCA, B, IB, JB,
     $                    DESCB, DWORK )
*
*     Contribution from the Department of Computing Science and HPC2N,
*     Umea University, Sweden
*
*  -- ScaLAPACK auxiliary routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, IB, JA, JB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      REAL               A( * ), B( * ), DWORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLAMVE copies all or part of a distributed matrix A to another
*  distributed matrix B. There is no alignment assumptions at all
*  except that A and B are of the same size.
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
*  UPLO    (global input) CHARACTER
*          Specifies the part of the distributed matrix sub( A ) to be
*          copied:
*          = 'U':   Upper triangular part is copied; the strictly
*                   lower triangular part of sub( A ) is not referenced;
*          = 'L':   Lower triangular part is copied; the strictly
*                   upper triangular part of sub( A ) is not referenced;
*          Otherwise:  All of the matrix sub( A ) is copied.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input) REAL             pointer into the local memory
*          to an array of dimension (LLD_A, LOCc(JA+N-1) ). This array
*          contains the local pieces of the distributed matrix sub( A )
*          to be copied from.
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
*  B       (local output) REAL             pointer into the local memory
*          to an array of dimension (LLD_B, LOCc(JB+N-1) ). This array
*          contains on exit the local pieces of the distributed matrix
*          sub( B ).
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
*  DWORK   (local workspace) REAL             array
*          If UPLO = 'U' or UPLO = 'L' and number of processors > 1,
*          the length of DWORK is at least as large as the length of B.
*          Otherwise, DWORK is not referenced.
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
      LOGICAL            UPPER, LOWER, FULL
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, MYPROC,
     $                   NPROCS, AROWS, ACOLS, K, SPROC, SRSRC, SCSRC,
     $                   RPROC, RRSRC, RCSRC, COUNT, J, I, IIA, JJA,
     $                   IIB, JJB, BRSRC, BCSRC, RAROWS, RACOLS,
     $                   INDEX, IDUM, NUMREC, NUMSND
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAMOV, INFOG2L
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, NUMROC, INDXL2G
      EXTERNAL           ICEIL, LSAME, NUMROC, INDXL2G
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Find underlying mesh properties.
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Decode input parameters.
*
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT. UPPER ) LOWER = LSAME( UPLO, 'L' )
      FULL = (.NOT. UPPER) .AND. (.NOT. LOWER)
*
*     Assign indiviual numbers based on column major ordering.
*
      NPROCS = NPROW*NPCOL
*
*     Do redistribution operation.
*
      IF( NPROCS.EQ.1 ) THEN
         CALL SLAMOV( UPLO, M, N, A((JA-1)*DESCA(LLD_)+IA),
     $        DESCA(LLD_), B((JB-1)*DESCB(LLD_)+IB),
     $        DESCB(LLD_) )
      ELSEIF( FULL ) THEN
         CALL PSGEMR2D( M, N, A, IA, JA, DESCA, B, IB, JB, DESCB,
     $        ICTXT )
      ELSE
         CALL PSGEMR2D( M, N, A, IA, JA, DESCA, DWORK, IB, JB, DESCB,
     $        ICTXT )
         CALL PSLACPY( UPLO, M, N, DWORK, IB, JB, DESCB, B, IB, JB,
     $        DESCB )
      END IF
*
      RETURN
*
*     End of PSLAMVE
*
      END
