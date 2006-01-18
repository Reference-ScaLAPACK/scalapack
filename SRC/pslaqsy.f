      SUBROUTINE PSLAQSY( UPLO, N, A, IA, JA, DESCA, SR, SC, SCOND,
     $                    AMAX, EQUED )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, UPLO
      INTEGER            IA, JA, N
      REAL               AMAX, SCOND
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      REAL               A( * ), SC( * ), SR( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLAQSY equilibrates a symmetric distributed matrix
*  sub( A ) = A(IA:IA+N-1,JA:JA+N-1) using the scaling factors in the
*  vectors SR and SC.
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
*          Specifies whether the upper or lower triangular part of the
*          symmetric distributed matrix sub( A ) is to be referenced:
*             = 'U':  Upper triangular
*             = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (input/output) REAL pointer into the local
*          memory to an array of local dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, the local pieces of the distributed symmetric
*          matrix sub( A ). If UPLO = 'U', the leading N-by-N upper
*          triangular part of sub( A ) contains the upper triangular
*          part of the matrix, and the strictly lower triangular part
*          of sub( A ) is not referenced.  If UPLO = 'L', the leading
*          N-by-N lower triangular part of sub( A ) contains the lower
*          triangular part of the matrix, and the strictly upper trian-
*          gular part of sub( A ) is not referenced.
*          On exit, if EQUED = 'Y', the equilibrated matrix:
*              diag(SR(IA:IA+N-1)) * sub( A ) * diag(SC(JA:JA+N-1)).
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
*  SR      (local input) REAL array, dimension LOCr(M_A)
*          The scale factors for A(IA:IA+M-1,JA:JA+N-1). SR is aligned
*          with the distributed matrix A, and replicated across every
*          process column. SR is tied to the distributed matrix A.
*
*  SC      (local input) REAL array, dimension LOCc(N_A)
*          The scale factors for sub( A ). SC is aligned with the dis-
*          tributed matrix A, and replicated down every process row.
*          SC is tied to the distributed matrix A.
*
*  SCOND   (global input) REAL
*          Ratio of the smallest SR(i) (respectively SC(j)) to the
*          largest SR(i) (respectively SC(j)), with IA <= i <= IA+N-1
*          and JA <= j <= JA+N-1.
*
*  AMAX    (global input) REAL
*          Absolute value of the largest distributed submatrix entry.
*
*  EQUED   (output) CHARACTER*1
*          Specifies whether or not equilibration was done.
*          = 'N':  No equilibration.
*          = 'Y':  Equilibration was done, i.e., sub( A ) has been re-
*                  placed by:
*                  diag(SR(IA:IA+N-1)) * sub( A ) * diag(SC(JA:JA+N-1)).
*
*  Internal Parameters
*  ===================
*
*  THRESH is a threshold value used to decide if scaling should be done
*  based on the ratio of the scaling factors.  If SCOND < THRESH,
*  scaling is done.
*
*  LARGE and SMALL are threshold values used to decide if scaling should
*  be done based on the absolute size of the largest matrix element.
*  If AMAX > LARGE or AMAX < SMALL, scaling is done.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ONE, THRESH
      PARAMETER          ( ONE = 1.0E+0, THRESH = 0.1E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IACOL, IAROW, ICTXT, II, IIA, IOFFA, IROFF, J,
     $                   JB, JJ, JJA, JN, KK, LDA, LL, MYCOL, MYROW, NP,
     $                   NPCOL, NPROW
      REAL               CJ, LARGE, SMALL
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, NUMROC
      REAL               PSLAMCH
      EXTERNAL           ICEIL, LSAME, NUMROC, PSLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF
*
*     Get grid parameters and compute local indexes
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      LDA = DESCA( LLD_ )
*
*     Initialize LARGE and SMALL.
*
      SMALL = PSLAMCH( ICTXT, 'Safe minimum' ) /
     $        PSLAMCH( ICTXT, 'Precision' )
      LARGE = ONE / SMALL
*
      IF( SCOND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE ) THEN
*
*        No equilibration
*
         EQUED = 'N'
*
      ELSE
*
         II = IIA
         JJ = JJA
         JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
         JB = JN-JA+1
*
*        Replace A by diag(S) * A * diag(S).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Upper triangle of A(IA:IA+N-1,JA:JA+N-1) is stored.
*           Handle first block separately
*
            IOFFA = (JJ-1)*LDA
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  DO 20 LL = JJ, JJ + JB -1
                     CJ = SC( LL )
                     DO 10 KK = IIA, II+LL-JJ+1
                        A( IOFFA + KK ) = CJ*SR( KK )*A( IOFFA + KK )
   10                CONTINUE
                     IOFFA = IOFFA + LDA
   20             CONTINUE
               ELSE
                  IOFFA = IOFFA + JB*LDA
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 70 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     DO 40 LL = JJ, JJ + JB -1
                        CJ = SC( LL )
                        DO 30 KK = IIA, II+LL-JJ+1
                           A( IOFFA + KK ) = CJ*SR( KK )*A( IOFFA + KK )
   30                   CONTINUE
                        IOFFA = IOFFA + LDA
   40                CONTINUE
                  ELSE
                     DO 60 LL = JJ, JJ + JB -1
                        CJ = SC( LL )
                        DO 50 KK = IIA, II-1
                           A( IOFFA + KK ) = CJ*SR( KK )*A( IOFFA + KK )
   50                   CONTINUE
                        IOFFA = IOFFA + LDA
   60                CONTINUE
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
   70       CONTINUE
*
         ELSE
*
*           Lower triangle of A(IA:IA+N-1,JA:JA+N-1) is stored.
*           Handle first block separately
*
            IROFF = MOD( IA-1, DESCA( MB_ ) )
            NP = NUMROC( N+IROFF, DESCA( MB_ ), MYROW, IAROW, NPROW )
            IF( MYROW.EQ.IAROW )
     $         NP = NP-IROFF
*
            IOFFA = (JJ-1)*LDA
            IF( MYCOL.EQ.IACOL ) THEN
               IF( MYROW.EQ.IAROW ) THEN
                  DO 90 LL = JJ, JJ + JB -1
                     CJ = SC( LL )
                     DO 80 KK = II+LL-JJ, IIA+NP-1
                        A( IOFFA + KK ) = CJ*SR( KK )*A( IOFFA + KK )
   80                CONTINUE
                     IOFFA = IOFFA + LDA
   90             CONTINUE
               ELSE
                  DO 110 LL = JJ, JJ + JB -1
                     CJ = SC( LL )
                     DO 100 KK = II, IIA+NP-1
                        A( IOFFA + KK ) = CJ*SR( KK )*A( IOFFA + KK )
  100                CONTINUE
                     IOFFA = IOFFA + LDA
  110             CONTINUE
               END IF
               JJ = JJ + JB
            END IF
*
            IF( MYROW.EQ.IAROW )
     $         II = II + JB
            IAROW = MOD( IAROW+1, NPROW )
            IACOL = MOD( IACOL+1, NPCOL )
*
*           Loop over remaining block of columns
*
            DO 160 J = JN+1, JA+N-1, DESCA( NB_ )
               JB = MIN( JA+N-J, DESCA( NB_ ) )
*
               IF( MYCOL.EQ.IACOL ) THEN
                  IF( MYROW.EQ.IAROW ) THEN
                     DO 130 LL = JJ, JJ + JB -1
                        CJ = SC( LL )
                        DO 120 KK = II+LL-JJ, IIA+NP-1
                           A( IOFFA + KK ) = CJ*SR( KK )*A( IOFFA + KK )
  120                   CONTINUE
                        IOFFA = IOFFA + LDA
  130                CONTINUE
                  ELSE
                     DO 150 LL = JJ, JJ + JB -1
                        CJ = SC( LL )
                        DO 140 KK = II, IIA+NP-1
                           A( IOFFA + KK ) = CJ*SR( KK )*A( IOFFA + KK )
  140                   CONTINUE
                        IOFFA = IOFFA + LDA
  150                CONTINUE
                  END IF
                  JJ = JJ + JB
               END IF
*
               IF( MYROW.EQ.IAROW )
     $            II = II + JB
               IAROW = MOD( IAROW+1, NPROW )
               IACOL = MOD( IACOL+1, NPCOL )
*
  160       CONTINUE
*
         END IF
*
         EQUED = 'Y'
*
      END IF
*
      RETURN
*
*     End of PSLAQSY
*
      END
