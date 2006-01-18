      SUBROUTINE PZLASE2( UPLO, M, N, ALPHA, BETA, A, IA, JA, DESCA )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            IA, JA, M, N
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         A( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLASE2 initializes an M-by-N distributed matrix sub( A ) denoting
*  A(IA:IA+M-1,JA:JA+N-1) to BETA on the diagonal and ALPHA on the
*  offdiagonals.  PZLASE2 requires that only dimension of the matrix
*  operand is distributed.
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
*          set:
*          = 'U':      Upper triangular part is set; the strictly lower
*                      triangular part of sub( A ) is not changed;
*          = 'L':      Lower triangular part is set; the strictly upper
*                      triangular part of sub( A ) is not changed;
*          Otherwise:  All of the matrix sub( A ) is set.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  ALPHA   (global input) COMPLEX*16
*          The constant to which the offdiagonal elements are to be
*          set.
*
*  BETA    (global input) COMPLEX*16
*          The constant to which the diagonal elements are to be set.
*
*  A       (local output) COMPLEX*16 pointer into the local memory
*          to an array of dimension (LLD_A,LOCc(JA+N-1)).  This array
*          contains the local pieces of the distributed matrix sub( A )
*          to be set.  On exit, the leading M-by-N submatrix sub( A )
*          is set as follows:
*
*          if UPLO = 'U', A(IA+i-1,JA+j-1) = ALPHA, 1<=i<=j-1, 1<=j<=N,
*          if UPLO = 'L', A(IA+i-1,JA+j-1) = ALPHA, j+1<=i<=M, 1<=j<=N,
*          otherwise,     A(IA+i-1,JA+j-1) = ALPHA, 1<=i<=M, 1<=j<=N,
*                                                   IA+i.NE.JA+j,
*          and, for all UPLO, A(IA+i-1,JA+i-1) = BETA, 1<=i<=min(M,N).
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
*     ..
*     .. Local Scalars ..
      INTEGER            HEIGHT, IACOL, IAROW, IBASE, ICOFFA, II, IIA,
     $                   IIBEG, IIEND, IINXT, ILEFT, IRIGHT, IROFFA,
     $                   ITOP, JJ, JJA, JJBEG, JJEND, JJNXT, LDA, MBA,
     $                   MP, MPA, MYCOL, MYDIST, MYROW, NBA, NPCOL,
     $                   NPROW, NQ, NQA, WIDE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, ZLASET
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
      MBA = DESCA( MB_ )
      NBA = DESCA( NB_ )
      LDA = DESCA( LLD_ )
      IROFFA = MOD( IA-1, MBA )
      ICOFFA = MOD( JA-1, NBA )
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
*         \      |...................|    |       ( IAROW )
*           IIA  |x                  |    | MB_A
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
            MPA = NUMROC( M+IROFFA, MBA, MYROW, IAROW, NPROW )
            IF( MPA.LE.0 )
     $         RETURN
            IF( MYROW.EQ.IAROW )
     $         MPA = MPA - IROFFA
            MYDIST = MOD( MYROW-IAROW+NPROW, NPROW )
            ITOP = MYDIST * MBA - IROFFA
*
            IF( LSAME( UPLO, 'U' ) ) THEN
*
               ITOP = MAX( 0, ITOP )
               IIBEG = IIA
               IIEND = IIA + MPA - 1
               IINXT = MIN( ICEIL( IIBEG, MBA ) * MBA, IIEND )
*
   10          CONTINUE
               IF( ( N-ITOP ).GT.0 ) THEN
                  CALL ZLASET( UPLO, IINXT-IIBEG+1, N-ITOP, ALPHA, BETA,
     $                         A( IIBEG+(JJA+ITOP-1)*LDA ), LDA )
                  MYDIST = MYDIST + NPROW
                  ITOP = MYDIST * MBA - IROFFA
                  IIBEG = IINXT +1
                  IINXT = MIN( IINXT+MBA, IIEND )
                  GO TO 10
               END IF
*
            ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
               II = IIA
               JJ = JJA
               MP = MPA
               IBASE = MIN( ITOP+MBA, N )
               ITOP = MIN( MAX( 0, ITOP ), N )
*
   20          CONTINUE
               IF( JJ.LE.( JJA+N-1 ) ) THEN
                  HEIGHT = IBASE - ITOP
                  CALL ZLASET( 'All', MP, ITOP-JJ+JJA, ALPHA, ALPHA,
     $                         A( II+(JJ-1)*LDA ), LDA )
                  CALL ZLASET( UPLO, MP, HEIGHT, ALPHA, BETA,
     $                         A( II+(JJA+ITOP-1)*LDA ), LDA )
                  MP = MAX( 0, MP - HEIGHT )
                  II = II + HEIGHT
                  JJ = JJA + IBASE
                  MYDIST = MYDIST + NPROW
                  ITOP = MYDIST * MBA - IROFFA
                  IBASE = MIN( ITOP + MBA, N )
                  ITOP = MIN( ITOP, N )
                  GO TO 20
               END IF
*
            ELSE
*
               II = IIA
               JJ = JJA
               MP = MPA
               IBASE = MIN( ITOP+MBA, N )
               ITOP = MIN( MAX( 0, ITOP ), N )
*
   30          CONTINUE
               IF( JJ.LE.( JJA+N-1 ) ) THEN
                  HEIGHT = IBASE - ITOP
                  CALL ZLASET( 'All', MPA, ITOP-JJ+JJA, ALPHA, ALPHA,
     $                         A( IIA+(JJ-1)*LDA ), LDA )
                  CALL ZLASET( 'All', MPA-MP, HEIGHT, ALPHA, ALPHA,
     $                         A( IIA+(JJA+ITOP-1)*LDA ), LDA )
                  CALL ZLASET( 'All', MP, HEIGHT, ALPHA, BETA,
     $                         A( II+(JJA+ITOP-1)*LDA ), LDA )
                  MP = MAX( 0, MP - HEIGHT )
                  II = II + HEIGHT
                  JJ = JJA + IBASE
                  MYDIST = MYDIST + NPROW
                  ITOP = MYDIST * MBA - IROFFA
                  IBASE = MIN( ITOP + MBA, N )
                  ITOP = MIN( ITOP, N )
                  GO TO 30
               END IF
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
            NQA = NUMROC( N+ICOFFA, NBA, MYCOL, IACOL, NPCOL )
            IF( NQA.LE.0 )
     $         RETURN
            IF( MYCOL.EQ.IACOL )
     $         NQA = NQA - ICOFFA
            MYDIST = MOD( MYCOL-IACOL+NPCOL, NPCOL )
            ILEFT = MYDIST * NBA - ICOFFA
*
            IF( LSAME( UPLO, 'L' ) ) THEN
*
               ILEFT = MAX( 0, ILEFT )
               JJBEG = JJA
               JJEND = JJA + NQA - 1
               JJNXT = MIN( ICEIL( JJBEG, NBA ) * NBA, JJEND )
*
   40          CONTINUE
               IF( ( M-ILEFT ).GT.0 ) THEN
                  CALL ZLASET( UPLO, M-ILEFT, JJNXT-JJBEG+1, ALPHA,
     $                         BETA, A( IIA+ILEFT+(JJBEG-1)*LDA ), LDA )
                  MYDIST = MYDIST + NPCOL
                  ILEFT = MYDIST * NBA - ICOFFA
                  JJBEG = JJNXT +1
                  JJNXT = MIN( JJNXT+NBA, JJEND )
                  GO TO 40
               END IF
*
            ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
               II = IIA
               JJ = JJA
               NQ = NQA
               IRIGHT = MIN( ILEFT+NBA, M )
               ILEFT = MIN( MAX( 0, ILEFT ), M )
*
   50          CONTINUE
               IF( II.LE.( IIA+M-1 ) ) THEN
                  WIDE = IRIGHT - ILEFT
                  CALL ZLASET( 'All', ILEFT-II+IIA, NQ, ALPHA, ALPHA,
     $                         A( II+(JJ-1)*LDA ), LDA )
                  CALL ZLASET( UPLO, WIDE, NQ, ALPHA, BETA,
     $                         A( IIA+ILEFT+(JJ-1)*LDA ), LDA )
                  NQ = MAX( 0, NQ - WIDE )
                  II = IIA + IRIGHT
                  JJ = JJ + WIDE
                  MYDIST = MYDIST + NPCOL
                  ILEFT = MYDIST * NBA - ICOFFA
                  IRIGHT = MIN( ILEFT + NBA, M )
                  ILEFT = MIN( ILEFT, M )
                  GO TO 50
               END IF
*
            ELSE
*
               II = IIA
               JJ = JJA
               NQ = NQA
               IRIGHT = MIN( ILEFT+NBA, M )
               ILEFT = MIN( MAX( 0, ILEFT ), M )
*
   60          CONTINUE
               IF( II.LE.( IIA+M-1 ) ) THEN
                  WIDE = IRIGHT - ILEFT
                  CALL ZLASET( 'All', ILEFT-II+IIA, NQA, ALPHA, ALPHA,
     $                         A( II+(JJA-1)*LDA ), LDA )
                  CALL ZLASET( 'All', WIDE, NQA-NQ, ALPHA, ALPHA,
     $                         A( IIA+ILEFT+(JJA-1)*LDA ), LDA )
                  CALL ZLASET( 'All', WIDE, NQ, ALPHA, BETA,
     $                         A( IIA+ILEFT+(JJ-1)*LDA ), LDA )
                  NQ = MAX( 0, NQ - WIDE )
                  II = IIA + IRIGHT
                  JJ = JJ + WIDE
                  MYDIST = MYDIST + NPCOL
                  ILEFT = MYDIST * NBA - ICOFFA
                  IRIGHT = MIN( ILEFT + NBA, M )
                  ILEFT = MIN( ILEFT, M )
                  GO TO 60
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
*     End of PZLASE2
*
      END
