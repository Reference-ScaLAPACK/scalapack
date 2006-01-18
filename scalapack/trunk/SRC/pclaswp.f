      SUBROUTINE PCLASWP( DIREC, ROWCOL, N, A, IA, JA, DESCA, K1, K2,
     $                    IPIV )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          DIREC, ROWCOL
      INTEGER            IA, JA, K1, K2, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), IPIV( * )
      COMPLEX            A( * )
*     ..
*
*  Purpose:
*  ========
*
*  PCLASWP performs a series of row or column interchanges on
*  the distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1).  One
*  interchange is initiated for each of rows or columns K1 trough K2 of
*  sub( A ). This routine assumes that the pivoting information has
*  already been broadcast along the process row or column.
*  Also note that this routine will only work for K1-K2 being in the
*  same MB (or NB) block.  If you want to pivot a full matrix, use
*  PCLAPIV.
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
*  DIREC   (global input) CHARACTER
*          Specifies in which order the permutation is applied:
*          = 'F' (Forward)
*          = 'B' (Backward)
*
*  ROWCOL  (global input) CHARACTER
*          Specifies if the rows or columns are permuted:
*          = 'R' (Rows)
*          = 'C' (Columns)
*
*  N       (global input) INTEGER
*          If ROWCOL = 'R', the length of the rows of the distributed
*          matrix A(*,JA:JA+N-1) to be permuted;
*          If ROWCOL = 'C', the length of the columns of the distributed
*          matrix A(IA:IA+N-1,*) to be permuted.
*
*  A       (local input/local output) COMPLEX pointer into the
*          local memory to an array of dimension (LLD_A, * ).
*          On entry, this array contains the local pieces of the distri-
*          buted matrix to which the row/columns interchanges will be
*          applied. On exit the permuted distributed matrix.
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
*  K1      (global input) INTEGER
*          The first element of IPIV for which a row or column inter-
*          change will be done.
*
*  K2      (global input) INTEGER
*          The last element of IPIV for which a row or column inter-
*          change will be done.
*
*  IPIV    (local input) INTEGER array, dimension LOCr(M_A)+MB_A for
*          row pivoting and LOCc(N_A)+NB_A for column pivoting.  This
*          array is tied to the matrix A, IPIV(K) = L implies rows
*          (or columns) K and L are to be interchanged.
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
      INTEGER            I, ICURCOL, ICURROW, IIA, IP, J, JJA, JP,
     $                   MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, PCSWAP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
      IF( LSAME( ROWCOL, 'R' ) ) THEN
         IF( LSAME( DIREC, 'F' ) ) THEN
            CALL INFOG2L( K1, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIA, JJA, ICURROW, ICURCOL )
            DO 10 I = K1, K2
               IP = IPIV( IIA+I-K1 )
               IF( IP.NE.I )
     $            CALL PCSWAP( N, A, I, JA, DESCA, DESCA( M_ ), A, IP,
     $                         JA, DESCA, DESCA( M_ ) )
   10       CONTINUE
         ELSE
            CALL INFOG2L( K2, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIA, JJA, ICURROW, ICURCOL )
            DO 20 I = K2, K1, -1
               IP = IPIV( IIA+I-K1 )
               IF( IP.NE.I )
     $            CALL PCSWAP( N, A, I, JA, DESCA, DESCA( M_ ), A, IP,
     $                         JA, DESCA, DESCA( M_ ) )
   20       CONTINUE
         END IF
      ELSE
         IF( LSAME( DIREC, 'F' ) ) THEN
            CALL INFOG2L( IA, K1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIA, JJA, ICURROW, ICURCOL )
            DO 30 J = K1, K2
               JP = IPIV( JJA+J-K1 )
               IF( JP.NE.J )
     $            CALL PCSWAP( N, A, IA, J, DESCA, 1, A, IA, JP,
     $                         DESCA, 1 )
   30       CONTINUE
         ELSE
            CALL INFOG2L( IA, K2, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIA, JJA, ICURROW, ICURCOL )
            DO 40 J = K2, K1, -1
               JP = IPIV( JJA+J-K1 )
               IF( JP.NE.J )
     $            CALL PCSWAP( N, A, IA, J, DESCA, 1, A, IA, JP,
     $                         DESCA, 1 )
   40       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End PCLASWP
*
      END
