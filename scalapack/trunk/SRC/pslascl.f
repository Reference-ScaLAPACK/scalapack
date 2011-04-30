      SUBROUTINE PSLASCL( TYPE, CFROM, CTO, M, N, A, IA, JA, DESCA,
     $                    INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            IA, INFO, JA, M, N
      REAL               CFROM, CTO
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      REAL               A( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLASCL multiplies the M-by-N real distributed matrix sub( A )
*  denoting A(IA:IA+M-1,JA:JA+N-1) by the real scalar CTO/CFROM.  This
*  is done without over/underflow as long as the final result
*  CTO * A(I,J) / CFROM does not over/underflow. TYPE specifies that
*  sub( A ) may be full, upper triangular, lower triangular or upper
*  Hessenberg.
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
*  TYPE    (global input) CHARACTER
*          TYPE indices the storage type of the input distributed
*          matrix.
*          = 'G':  sub( A ) is a full matrix,
*          = 'L':  sub( A ) is a lower triangular matrix,
*          = 'U':  sub( A ) is an upper triangular matrix,
*          = 'H':  sub( A ) is an upper Hessenberg matrix.
*
*  CFROM   (global input) REAL
*  CTO     (global input) REAL
*          The distributed matrix sub( A ) is multiplied by CTO/CFROM.
*          A(I,J) is computed without over/underflow if the final
*          result CTO * A(I,J) / CFROM can be represented without
*          over/underflow.  CFROM must be nonzero.
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) REAL pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          This array contains the local pieces of the distributed
*          matrix sub( A ). On exit, this array contains the local
*          pieces of the distributed matrix multiplied by CTO/CFROM.
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
*  INFO    (local output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ONE, ZERO
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            IACOL, IAROW, ICOFFA, ICTXT, ICURCOL, ICURROW,
     $                   IIA, II, INXTROW, IOFFA, IROFFA, ITYPE, J, JB,
     $                   JJA, JJ, JN, KK, LDA, LL, MYCOL, MYROW, MP,
     $                   NPCOL, NPROW, NQ
      REAL               BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, INFOG2L, PXERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, SISNAN
      INTEGER            ICEIL, NUMROC
      REAL               PSLAMCH
      EXTERNAL           SISNAN, ICEIL, LSAME, NUMROC, PSLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      IF( NPROW.EQ.-1 ) THEN
         INFO = -907
      ELSE
         INFO = 0
         CALL CHK1MAT( M, 4, N, 6, IA, JA, DESCA, 9, INFO )
         IF( INFO.EQ.0 ) THEN
            IF( LSAME( TYPE, 'G' ) ) THEN
               ITYPE = 0
            ELSE IF( LSAME( TYPE, 'L' ) ) THEN
               ITYPE = 1
            ELSE IF( LSAME( TYPE, 'U' ) ) THEN
               ITYPE = 2
            ELSE IF( LSAME( TYPE, 'H' ) ) THEN
               ITYPE = 3
            ELSE
               ITYPE = -1
            END IF
            IF( ITYPE.EQ.-1 ) THEN
               INFO = -1
            ELSE IF( CFROM.EQ.ZERO .OR. SISNAN(CFROM) ) THEN
               INFO = -4
            ELSE IF( SISNAN(CTO) ) THEN
               INFO = -5
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PSLASCL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = PSLAMCH( ICTXT, 'S' )
      BIGNUM = ONE / SMLNUM
*
      CFROMC = CFROM
      CTOC = CTO
*
*     Compute local indexes
*
      LDA = DESCA( LLD_ )
      IROFFA = MOD( IA-1, DESCA( MB_ ) )
      ICOFFA = MOD( JA-1, DESCA( NB_ ) )
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      MP = NUMROC( M+IROFFA, DESCA( MB_ ), MYROW, IAROW, NPROW )
      IF( MYROW.EQ.IAROW )
     $   MP = MP - IROFFA
      NQ = NUMROC( N+ICOFFA, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IF( MYCOL.EQ.IACOL )
     $   NQ = NQ - ICOFFA
*
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      IF( CFROM1.EQ.CFROMC ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = .TRUE.
         CTO1 = CTOC
      ELSE
         CTO1 = CTOC / BIGNUM
         IF( CTO1.EQ.CTOC ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            MUL = CTOC
            DONE = .TRUE.
            CFROMC = ONE
         ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
            MUL = SMLNUM
            DONE = .FALSE.
            CFROMC = CFROM1
         ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
            MUL = BIGNUM
            DONE = .FALSE.
            CTOC = CTO1
         ELSE
            MUL = CTOC / CFROMC
            DONE = .TRUE.
         END IF
      END IF
*
      IOFFA = ( JJA - 1 ) * LDA
      ICURROW = IAROW
      ICURCOL = IACOL
*
      IF( ITYPE.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 JJ = JJA, JJA+NQ-1
            DO 20 II = IIA, IIA+MP-1
               A( IOFFA+II ) = A( IOFFA+II ) * MUL
   20       CONTINUE
            IOFFA = IOFFA + LDA
   30    CONTINUE
*
      ELSE IF( ITYPE.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         II = IIA
         JJ = JJA
         JB = JN-JA+1
*
         IF( MYCOL.EQ.ICURCOL ) THEN
            IF( MYROW.EQ.ICURROW ) THEN
               DO 50 LL = JJ, JJ + JB -1
                  DO 40 KK = II+LL-JJ, IIA+MP-1
                     A( IOFFA+KK ) = A( IOFFA+KK ) * MUL
   40             CONTINUE
                  IOFFA = IOFFA + LDA
   50          CONTINUE
            ELSE
               DO 70 LL = JJ, JJ + JB -1
                  DO 60 KK = II, IIA+MP-1
                     A( IOFFA+KK ) = A( IOFFA+KK ) * MUL
   60             CONTINUE
                  IOFFA = IOFFA + LDA
   70          CONTINUE
            END IF
            JJ = JJ + JB
         END IF
*
         IF( MYROW.EQ.ICURROW )
     $      II = II + JB
         ICURROW = MOD( ICURROW+1, NPROW )
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
*        Loop over remaining block of columns
*
         DO 120 J = JN+1, JA+N-1, DESCA( NB_ )
            JB = MIN( JA+N-J, DESCA( NB_ ) )
*
            IF( MYCOL.EQ.ICURCOL ) THEN
               IF( MYROW.EQ.ICURROW ) THEN
                  DO 90 LL = JJ, JJ + JB -1
                     DO 80 KK = II+LL-JJ, IIA+MP-1
                        A( IOFFA+KK ) = A( IOFFA+KK ) * MUL
   80                CONTINUE
                     IOFFA = IOFFA + LDA
   90             CONTINUE
               ELSE
                  DO 110 LL = JJ, JJ + JB -1
                     DO 100 KK = II, IIA+MP-1
                        A( IOFFA+KK ) = A( IOFFA+KK ) * MUL
  100                CONTINUE
                     IOFFA = IOFFA + LDA
  110             CONTINUE
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.ICURROW )
     $         II = II + JB
            ICURROW = MOD( ICURROW+1, NPROW )
            ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  120    CONTINUE
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         II = IIA
         JJ = JJA
         JB = JN-JA+1
*
         IF( MYCOL.EQ.ICURCOL ) THEN
            IF( MYROW.EQ.ICURROW ) THEN
               DO 140 LL = JJ, JJ + JB -1
                  DO 130 KK = IIA, MIN(II+LL-JJ,IIA+MP-1)
                     A( IOFFA+KK ) = A( IOFFA+KK ) * MUL
  130             CONTINUE
                  IOFFA = IOFFA + LDA
  140          CONTINUE
            ELSE
               DO 160 LL = JJ, JJ + JB -1
                  DO 150 KK = IIA, MIN(II-1,IIA+MP-1)
                     A( IOFFA+KK ) = A( IOFFA+KK ) * MUL
  150             CONTINUE
                  IOFFA = IOFFA + LDA
  160          CONTINUE
            END IF
            JJ = JJ + JB
         END IF
*
         IF( MYROW.EQ.ICURROW )
     $      II = II + JB
         ICURROW = MOD( ICURROW+1, NPROW )
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
*        Loop over remaining block of columns
*
         DO 210 J = JN+1, JA+N-1, DESCA( NB_ )
            JB = MIN( JA+N-J, DESCA( NB_ ) )
*
            IF( MYCOL.EQ.ICURCOL ) THEN
               IF( MYROW.EQ.ICURROW ) THEN
                  DO 180 LL = JJ, JJ + JB -1
                     DO 170 KK = IIA, MIN(II+LL-JJ,IIA+MP-1)
                        A( IOFFA+KK ) = A( IOFFA+KK )*MUL
  170                CONTINUE
                     IOFFA = IOFFA + LDA
  180             CONTINUE
               ELSE
                  DO 200 LL = JJ, JJ + JB -1
                     DO 190 KK = IIA, MIN(II-1,IIA+MP-1)
                        A( IOFFA+KK ) = A( IOFFA+KK ) * MUL
  190                CONTINUE
                     IOFFA = IOFFA + LDA
  200             CONTINUE
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.ICURROW )
     $         II = II + JB
            ICURROW = MOD( ICURROW+1, NPROW )
            ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  210    CONTINUE
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         II = IIA
         JJ = JJA
         JB = JN-JA+1
*
*        Only one process row
*
         IF( NPROW.EQ.1 ) THEN
*
*           Handle first block of columns separately
*
            IF( MYCOL.EQ.ICURCOL ) THEN
               DO 230 LL = JJ, JJ+JB-1
                  DO 220 KK = IIA, MIN( II+LL-JJ+1, IIA+MP-1 )
                     A( IOFFA+KK ) = A( IOFFA+KK )*MUL
  220             CONTINUE
                  IOFFA = IOFFA + LDA
  230          CONTINUE
               JJ = JJ + JB
            END IF
*
            ICURCOL = MOD( ICURCOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 260 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.ICURCOL ) THEN
                  DO 250 LL = JJ, JJ+JB-1
                     DO 240 KK = IIA, MIN( II+LL-JJ+1, IIA+MP-1 )
                        A( IOFFA+KK ) = A( IOFFA+KK )*MUL
  240                CONTINUE
                     IOFFA = IOFFA + LDA
  250             CONTINUE
                  JJ = JJ + JB
               END IF
*
               II = II + JB
               ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  260       CONTINUE
*
         ELSE
*
*           Handle first block of columns separately
*
            INXTROW = MOD( ICURROW+1, NPROW )
            IF( MYCOL.EQ.ICURCOL ) THEN
               IF( MYROW.EQ.ICURROW ) THEN
                  DO 280 LL = JJ, JJ + JB -1
                     DO 270 KK = IIA, MIN(II+LL-JJ+1,IIA+MP-1)
                        A( IOFFA+KK ) = A( IOFFA+KK ) * MUL
  270                CONTINUE
                     IOFFA = IOFFA + LDA
  280             CONTINUE
               ELSE
                  DO 300 LL = JJ, JJ + JB -1
                     DO 290 KK = IIA, MIN(II-1,IIA+MP-1)
                        A( IOFFA+KK ) = A( IOFFA+KK ) * MUL
  290                CONTINUE
                     IOFFA = IOFFA + LDA
  300             CONTINUE
                  IF( MYROW.EQ.INXTROW .AND. II.LE.IIA+MP-1 )
     $               A( II+(JJ+JB-2)*LDA ) = A( II+(JJ+JB-2)*LDA ) * MUL
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.ICURROW )
     $         II = II + JB
            ICURROW = INXTROW
            ICURROW = MOD( ICURROW+1, NPROW )
            ICURCOL = MOD( ICURCOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 350 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.ICURCOL ) THEN
                  IF( MYROW.EQ.ICURROW ) THEN
                     DO 320 LL = JJ, JJ + JB -1
                        DO 310 KK = IIA, MIN( II+LL-JJ+1, IIA+MP-1 )
                           A( IOFFA+KK ) = A( IOFFA+KK ) * MUL
  310                   CONTINUE
                        IOFFA = IOFFA + LDA
  320                CONTINUE
                  ELSE
                     DO 340 LL = JJ, JJ + JB -1
                        DO 330 KK = IIA, MIN( II-1, IIA+MP-1 )
                           A( IOFFA+KK ) = A( IOFFA+KK ) * MUL
  330                   CONTINUE
                        IOFFA = IOFFA + LDA
  340                CONTINUE
                     IF( MYROW.EQ.INXTROW .AND. II.LE.IIA+MP-1 )
     $                  A( II+(JJ+JB-2)*LDA ) = A( II+(JJ+JB-2)*LDA ) *
     $                                          MUL
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.ICURROW )
     $            II = II + JB
               ICURROW = INXTROW
               ICURROW = MOD( ICURROW+1, NPROW )
               ICURCOL = MOD( ICURCOL+1, NPCOL )
*
  350       CONTINUE
*
         END IF
*
      END IF
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of PSLASCL
*
      END
