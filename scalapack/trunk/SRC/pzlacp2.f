      SUBROUTINE PZLACP2( UPLO, M, N, A, IA, JA, DESCA, B, IB, JB,
     $                     DESCB )
*
*  -- ScaLAPACK auxiliary routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, IB, JA, JB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX*16         A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLACP2 copies all or part of a distributed matrix A to another
*  distributed matrix B.  No communication is performed, PZLACP2
*  performs a local copy sub( A ) := sub( B ), where sub( A ) denotes
*  A(IA:IA+M-1,JA:JA+N-1) and sub( B ) denotes B(IB:IB+M-1,JB:JB+N-1).
*  PZLACP2 requires that only dimension of the matrix operands is
*  distributed.
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
*  A       (local input) COMPLEX*16 pointer into the local memory
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
*  B       (local output) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_B, LOCc(JB+N-1) ). This array
*          contains on exit the local pieces of the distributed matrix
*          sub( B ) set as follows:
*
*          if UPLO = 'U', B(IB+i-1,JB+j-1) = A(IA+i-1,JA+j-1),
*                         1<=i<=j, 1<=j<=N;
*          if UPLO = 'L', B(IB+i-1,JB+j-1) = A(IA+i-1,JA+j-1),
*                         j<=i<=M, 1<=j<=N;
*          otherwise,     B(IB+i-1,JB+j-1) = A(IA+i-1,JA+j-1),
*                         1<=i<=M, 1<=j<=N.
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
      INTEGER            HEIGHT, IACOL, IAROW, IBASE, IBCOL, IBROW,
     $                   ICOFFA, IIA, IIAA, IIB, IIBB, IIBEGA, IIBEGB,
     $                   IIENDA, IINXTA, IINXTB, ILEFT, IRIGHT, IROFFA,
     $                   ITOP, JJA, JJAA, JJB, JJBB, JJBEGA, JJBEGB,
     $                   JJENDA, JJNXTA, JJNXTB, LDA, LDB, MBA, MP,
     $                   MPAA, MYCOL, MYDIST, MYROW, NBA, NPCOL, NPROW,
     $                   NQ, NQAA, WIDE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, ZLAMOV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, NUMROC
      EXTERNAL           ICEIL, LSAME, NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Get grid parameters
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      CALL INFOG2L( IB, JB, DESCB, NPROW, NPCOL, MYROW, MYCOL, IIB, JJB,
     $              IBROW, IBCOL )
*
      MBA    = DESCA( MB_ )
      NBA    = DESCA( NB_ )
      LDA    = DESCA( LLD_ )
      IROFFA = MOD( IA-1, MBA )
      ICOFFA = MOD( JA-1, NBA )
      LDB    = DESCB( LLD_ )
*
      IF( N.LE.( NBA-ICOFFA ) ) THEN
*
*        It is assumed that the local columns JJA:JJA+N-1 of the matrix
*        A are in the same process column (IACOL).
*
*                         N
*                JJA             JJA+N-1
*         /      ---------------------    \
*   IROFFA|      |                   |    |
*         \      |...................|    |          ( IAROW )
*           IIA  |x                  |    |   MBA = DESCA( MB_ )
*                | x                 |    |
*                |--x----------------|    /
*                |   x               |
*                |    x              |        ITOP
*                |     x             |          |
*                |      x            |      /-------\
*                |-------x-----------|      |-------x-----------|
*                |        x          |      |        x          |
*                |         x         |      |         x         |
*                |          x        |      |          x        |
*                |           x       |      |           x       |
*                |------------x------|      |------------x------|
*                |             x     |      \____________/
*                |              x    |            |
*                |               x   |          IBASE
*                |                x  |
*                |-----------------x-|          Local picture
*                |                  x|
*                |                   |
*                |                   |
*                |                   |
*                |-------------------|
*                |                   |
*                .                   .
*                .                   .
*                .      (IACOL)      .
*
         IF( MYCOL.EQ.IACOL ) THEN
*
            MP = NUMROC( M+IROFFA, MBA, MYROW, IAROW, NPROW )
            IF( MP.LE.0 )
     $         RETURN
            IF( MYROW.EQ.IAROW )
     $         MP = MP - IROFFA
            MYDIST = MOD( MYROW-IAROW+NPROW, NPROW )
            ITOP   = MYDIST * MBA - IROFFA
*
            IF( LSAME( UPLO, 'U' ) ) THEN
*
               ITOP   = MAX( 0, ITOP )
               IIBEGA = IIA
               IIENDA = IIA + MP - 1
               IINXTA = MIN( ICEIL( IIBEGA, MBA ) * MBA, IIENDA )
               IIBEGB = IIB
               IINXTB = IIBEGB + IINXTA - IIBEGA
*
   10          CONTINUE
               IF( ( N-ITOP ).GT.0 ) THEN
                  CALL ZLAMOV( UPLO, IINXTA-IIBEGA+1, N-ITOP,
     $                         A( IIBEGA+(JJA+ITOP-1)*LDA ), LDA,
     $                         B( IIBEGB+(JJB+ITOP-1)*LDB ), LDB )
                  MYDIST = MYDIST + NPROW
                  ITOP   = MYDIST * MBA - IROFFA
                  IIBEGA = IINXTA + 1
                  IINXTA = MIN( IINXTA+MBA, IIENDA )
                  IIBEGB = IINXTB + 1
                  IINXTB = IIBEGB + IINXTA - IIBEGA
                  GO TO 10
               END IF
*
            ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
               MPAA  = MP
               IIAA  = IIA
               JJAA  = JJA
               IIBB  = IIB
               JJBB  = JJB
               IBASE = MIN( ITOP + MBA, N )
               ITOP  = MIN( MAX( 0, ITOP ), N )
*
   20          CONTINUE
               IF( JJAA.LE.( JJA+N-1 ) ) THEN
                  HEIGHT = IBASE - ITOP
                  CALL ZLAMOV( 'All', MPAA, ITOP-JJAA+JJA,
     $                         A( IIAA+(JJAA-1)*LDA ), LDA,
     $                         B( IIBB+(JJBB-1)*LDB ), LDB )
                  CALL ZLAMOV( UPLO, MPAA, HEIGHT,
     $                         A( IIAA+(JJA+ITOP-1)*LDA ), LDA,
     $                         B( IIBB+(JJB+ITOP-1)*LDB ), LDB )
                  MPAA   = MAX( 0, MPAA - HEIGHT )
                  IIAA   = IIAA + HEIGHT
                  JJAA   = JJA  + IBASE
                  IIBB   = IIBB + HEIGHT
                  JJBB   = JJB  + IBASE
                  MYDIST = MYDIST + NPROW
                  ITOP   = MYDIST * MBA - IROFFA
                  IBASE  = MIN( ITOP + MBA, N )
                  ITOP   = MIN( ITOP, N )
                  GO TO 20
               END IF
*
            ELSE
*
               CALL ZLAMOV( 'All', MP, N, A( IIA+(JJA-1)*LDA ),
     $                      LDA, B( IIB+(JJB-1)*LDB ), LDB )
*
            END IF
*
         END IF
*
      ELSE IF( M.LE.( MBA-IROFFA ) ) THEN
*
*        It is assumed that the local rows IIA:IIA+M-1 of the matrix A
*        are in the same process row (IAROW).
*
*            ICOFFA
*             / \JJA
*        IIA  ------------------ ....            --------
*             | .x  |    |    |                 / |    | \
*             | . x |    |    |            ILEFT| |    | |
*             | .  x     |    |                 | |    | |
*             | .   x    |    |                 \ x    | |
*             | .   |x   |    |                   |x   | | IRIGHT
*             | .   | x  |    |                   | x  | |
*    (IAROW)  | .   |  x |    |                   |  x | |
*             | .   |   x|    |                   |   x| |
*             | .   |    x    |                   |    x /
*             | .   |    |x   |                   |    |
*             | .   |    | x  |                   |    |
*             | .   |    |  x |                   |    |
*             | .   |    |   x|                   |    |
*    IIA+M-1  ------------------ ....            -------
*              NB_A
*             (IACOL)                          Local picture
*
         IF( MYROW.EQ.IAROW ) THEN
*
            NQ = NUMROC( N+ICOFFA, NBA, MYCOL, IACOL, NPCOL )
            IF( NQ.LE.0 )
     $         RETURN
            IF( MYCOL.EQ.IACOL )
     $         NQ = NQ - ICOFFA
            MYDIST = MOD( MYCOL-IACOL+NPCOL, NPCOL )
            ILEFT  = MYDIST * NBA - ICOFFA
*
            IF( LSAME( UPLO, 'L' ) ) THEN
*
               ILEFT  = MAX( 0, ILEFT )
               JJBEGA = JJA
               JJENDA = JJA + NQ - 1
               JJNXTA = MIN( ICEIL( JJBEGA, NBA ) * NBA, JJENDA )
               JJBEGB = JJB
               JJNXTB = JJBEGB + JJNXTA - JJBEGA
*
   30          CONTINUE
               IF( ( M-ILEFT ).GT.0 ) THEN
                  CALL ZLAMOV( UPLO, M-ILEFT, JJNXTA-JJBEGA+1,
     $                         A( IIA+ILEFT+(JJBEGA-1)*LDA ), LDA,
     $                         B( IIB+ILEFT+(JJBEGB-1)*LDB ), LDB )
                  MYDIST = MYDIST + NPCOL
                  ILEFT  = MYDIST * NBA - ICOFFA
                  JJBEGA = JJNXTA +1
                  JJNXTA = MIN( JJNXTA+NBA, JJENDA )
                  JJBEGB = JJNXTB +1
                  JJNXTB = JJBEGB + JJNXTA - JJBEGA
                  GO TO 30
               END IF
*
            ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
               NQAA   = NQ
               IIAA   = IIA
               JJAA   = JJA
               IIBB   = IIB
               JJBB   = JJB
               IRIGHT = MIN( ILEFT + NBA, M )
               ILEFT  = MIN( MAX( 0, ILEFT ), M )
*
   40          CONTINUE
               IF( IIAA.LE.( IIA+M-1 ) ) THEN
                  WIDE = IRIGHT - ILEFT
                  CALL ZLAMOV( 'All', ILEFT-IIAA+IIA, NQAA,
     $                         A( IIAA+(JJAA-1)*LDA ), LDA,
     $                         B( IIBB+(JJBB-1)*LDB ), LDB )
                  CALL ZLAMOV( UPLO, WIDE, NQAA,
     $                         A( IIA+ILEFT+(JJAA-1)*LDA ), LDA,
     $                         B( IIB+ILEFT+(JJBB-1)*LDB ), LDB )
                  NQAA   = MAX( 0, NQAA - WIDE )
                  IIAA   = IIA  + IRIGHT
                  JJAA   = JJAA + WIDE
                  IIBB   = IIB  + IRIGHT
                  JJBB   = JJBB + WIDE
                  MYDIST = MYDIST + NPCOL
                  ILEFT  = MYDIST * NBA - ICOFFA
                  IRIGHT = MIN( ILEFT + NBA, M )
                  ILEFT  = MIN( ILEFT, M )
                  GO TO 40
               END IF
*
            ELSE
*
               CALL ZLAMOV( 'All', M, NQ, A( IIA+(JJA-1)*LDA ),
     $                      LDA, B( IIB+(JJB-1)*LDB ), LDB )
*
            END IF
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PZLACP2
*
      END
