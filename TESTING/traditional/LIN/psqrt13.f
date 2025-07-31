      SUBROUTINE PSQRT13( SCALE, M, N, A, IA, JA, DESCA, NORMA, ISEED,
     $                    WORK )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, ISEED, JA, M, N, SCALE
      REAL               NORMA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      REAL               WORK( * )
      REAL               A( * )
*     ..
*
*  Purpose
*  =======
*
*  PSQRT13 generates a full-rank matrix that may be scaled to have
*  large or small norm.
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
*  SCALE   (global input) INTEGER
*          SCALE = 1: normally scaled matrix
*          SCALE = 2: matrix scaled up
*          SCALE = 3: matrix scaled down
*
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local output) REAL pointer into the local memory
*          to an array of dimension (LLD_A,LOCc(JA+N-1)). This array
*          contains the local pieces of the distributed matrix sub( A ).
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
*  NORMA   (global output) REAL
*          The one-norm of A.
*
*  ISEED   (global input/global output) INTEGER
*          Seed for random number generator.
*
*  WORK    (local workspace) REAL               array, dimension (LWORK)
*          LWORK >= Nq0, where
*
*          ICOFFA = MOD( JA-1, NB_A ),
*          IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ), and
*          Nq0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ).
*
*          INDXG2P and NUMROC are ScaLAPACK tool functions; MYROW,
*          MYCOL, NPROW and NPCOL can be determined by calling the
*          subroutine BLACS_GRIDINFO.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IACOL, IAROW, ICOFFA, ICTXT, IIA, INFO,
     $                   IROFFA, J, JJA, MP, MYCOL, MYROW, NPCOL,
     $                   NPROW, NQ
      REAL               AJJ, ASUM, BIGNUM, SMLNUM
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      REAL               PSLAMCH, PSLANGE
      EXTERNAL           NUMROC, PSLAMCH, PSLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, PSLABAD, PSLASCL,
     $                   PSMATGEN, PSASUM, PSELGET, PSELSET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD, SIGN
*     ..
*     .. Executable Statements ..
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
*     generate the matrix
*
      IROFFA = MOD( IA-1, DESCA( MB_ ) )
      ICOFFA = MOD( JA-1, DESCA( NB_ ) )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA,
     $              JJA, IAROW, IACOL )
      MP = NUMROC( M+IROFFA, DESCA( MB_ ), MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFFA, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IF( MYROW.EQ.IAROW )
     $   MP = MP - IROFFA
      IF( MYCOL.EQ.IACOL )
     $   NQ = NQ  - ICOFFA
*
      CALL PSMATGEN( ICTXT, 'N', 'N', DESCA( M_ ), DESCA( N_ ),
     $               DESCA( MB_ ), DESCA( NB_ ), A, DESCA( LLD_ ),
     $               DESCA( RSRC_ ), DESCA( CSRC_ ), ISEED, IIA-1, MP,
     $               JJA-1, NQ, MYROW, MYCOL, NPROW, NPCOL )
*
      DO 10 J = JA, JA+N-1
         I = IA + J - JA
         IF( I.LE.IA+M-1 ) THEN
            CALL PSASUM( M, ASUM, A, IA, J, DESCA, 1 )
            CALL PSELGET( 'Column', ' ', AJJ, A, I, J, DESCA )
            AJJ = AJJ + SIGN( ASUM, AJJ )
            CALL PSELSET( A, I, J, DESCA, AJJ )
         END IF
   10 CONTINUE
*
*     scaled versions
*
      IF( SCALE.NE.1 ) THEN
*
         NORMA = PSLANGE( 'M', M, N, A, IA, JA, DESCA, WORK )
         SMLNUM = PSLAMCH( ICTXT, 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
         CALL PSLABAD( ICTXT, SMLNUM, BIGNUM )
         SMLNUM = SMLNUM / PSLAMCH( ICTXT, 'Epsilon' )
         BIGNUM = ONE / SMLNUM
*
         IF( SCALE.EQ.2 ) THEN
*
*           matrix scaled up
*
            CALL PSLASCL( 'General', NORMA, BIGNUM, M, N, A, IA,
     $                   JA, DESCA, INFO )
*
         ELSE IF( SCALE.EQ.3 ) THEN
*
*           matrix scaled down
*
            CALL PSLASCL( 'General', NORMA, SMLNUM, M, N, A, IA,
     $                   JA, DESCA, INFO )
*
         END IF
*
      END IF
*
      NORMA = PSLANGE( 'One-norm', M, N, A, IA, JA, DESCA, WORK )
*
      RETURN
*
*     End of PSQRT13
*
      END
