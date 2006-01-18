      SUBROUTINE DESCSET( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT,
     $                    LLD )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            ICSRC, ICTXT, IRSRC, LLD, M, MB, N, NB
*     ..
*     .. Array Arguments ..
      INTEGER            DESC( * )
*     ..
*
*  Purpose
*  =======
*
*  DESCSET initializes a descriptor vector with the 8 input arguments
*  M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD.
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
*  DESC    (output) INTEGER array of dimension DLEN.
*          The array descriptor of a distributed matrix to be set.
*
*  M       (global input) INTEGER
*          The number of rows in the distributed matrix. M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns in the distributed matrix. N >= 0.
*
*  MB      (global input) INTEGER
*          The blocking factor used to distribute the rows of the
*          matrix. MB >= 1.
*
*  NB      (global input) INTEGER
*          The blocking factor used to distribute the columns of the
*          matrix. NB >= 1.
*
*  IRSRC   (global input) INTEGER
*          The process row over which the first row of the matrix is
*          distributed. 0 <= IRSRC < NPROW.
*
*  ICSRC   (global input) INTEGER
*          The process column over which the first column of the
*          matrix is distributed. 0 <= ICSRC < NPCOL.
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation on the matrix. The context itself is global.
*
*  LLD     (local input)  INTEGER
*          The leading dimension of the local array storing the local
*          blocks of the distributed matrix. LLD >= MAX(1,LOCr(M)).
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
*
*     .. Executable Statements ..
*
      DESC( DTYPE_ ) = BLOCK_CYCLIC_2D
      DESC( M_ ) = M
      DESC( N_ ) = N
      DESC( MB_ ) = MB
      DESC( NB_ ) = NB
      DESC( RSRC_ ) = IRSRC
      DESC( CSRC_ ) = ICSRC
      DESC( CTXT_ ) = ICTXT
      DESC( LLD_ ) = LLD
*
      RETURN
*
*     End DESCSET
*
      END
