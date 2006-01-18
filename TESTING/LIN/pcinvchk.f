      SUBROUTINE PCINVCHK( MATTYP, N, A, IA, JA, DESCA, IASEED, ANORM,
     $                     FRESID, RCOND, WORK )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 28, 2001
*
*     .. Scalar Arguments ..
      INTEGER            IA, IASEED, JA, N
      REAL               ANORM, FRESID, RCOND
*     ..
*     .. Array Arguments ..
      CHARACTER*3        MATTYP
      INTEGER            DESCA( * )
      COMPLEX            A( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PCINVCHK computes the scaled residual
*
*  || sub( A ) * inv( sub( A ) ) - I || / ( || sub( A ) || * N * eps ),
*
*  where sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1). to check the result
*  returned by the matrix inversion routines.
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
*  MATTYP   (global input) CHARACTER*3
*           The type of the distributed matrix to be generated:
*           if MATTYP = 'GEN' then GENeral matrix,
*           if MATTYP = 'UTR' then Upper TRiangular matrix,
*           if MATTYP = 'LTR' then Lower TRiangular matrix,
*           if MATTYP = 'UPD' then (Upper) Hermitian Positive Definite,
*           if MATTYP = 'LPD' then (Lower) Hermitian Positive Definite.
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input) COMPLEX pointer into the local memory
*          to an array of local dimension (LLD_A, LOCc(JA+N-1)).  On
*          entry, sub( A ) contains the distributed matrix inverse
*          computed by PCGETRI, PCPOTRI or PCTRTRI.
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
*  IASEED  (global input) INTEGER
*          Seed for the random generation of sub( A ).
*
*  ANORM   (global input) REAL
*          The 1-norm of the original matrix sub( A ).
*
*  FRESID  (global output) REAL
*          The inversion residual.
*
*  RCOND   (global output) REAL
*          The condition number of the original distributed matrix.
*          RCOND = || sub( A ) ||.|| sub( A )^{-1} || where ||A||
*          denotes the 1-norm of A.
*
*  WORK    (local workspace) COMPLEX array, dimension
*             MAX(2*LOCr(N_A+MOD(IA-1,MB_A))*MB_A, LDW)
*          where LDW is the workspace requirement for the norm computa-
*          tions, see PCLANGE, PCLANHE, PCLANSY and PCLANTR.
*
* =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      COMPLEX            ZERO,          ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          AFORM, DIAG, UPLO
      INTEGER            ICTXT, ICURCOL, ICURROW, II, IIA, IPW, IROFF,
     $                   IW, J, JB, JJA, JN, KK, MYCOL, MYROW, NP,
     $                   NPCOL, NPROW
      REAL               AUXNORM, EPS, NRMINVAXA, TEMP
*     ..
*     .. Local Arrays ..
      INTEGER            DESCW( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, INFOG2L, PCGEMM,
     $                   PCHEMM, PCLASET, PCMATGEN, PCTRMM
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      INTEGER            ICEIL, NUMROC
      REAL               PCLANGE, PCLANHE, PCLANTR, PSLAMCH
      EXTERNAL           ICEIL, LSAMEN, NUMROC, PCLANGE, PCLANHE,
     $                   PCLANSY, PCLANTR, PSLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
      EPS = PSLAMCH( DESCA( CTXT_ ), 'eps' )
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Compute the condition number
*
      IF( LSAMEN( 1, MATTYP( 1:1 ), 'U' ) ) THEN
         UPLO = 'U'
      ELSE
         UPLO = 'L'
      END IF
*
      IF( LSAMEN( 3, MATTYP, 'GEN' ) ) THEN
*
         AFORM = 'N'
         DIAG = 'D'
         AUXNORM = PCLANGE( '1', N, N, A, IA, JA, DESCA, WORK )
*
      ELSE IF( LSAMEN( 2, MATTYP( 2:3 ), 'TR' ) ) THEN
*
         AFORM = 'N'
         DIAG = 'D'
         AUXNORM = PCLANTR( '1', UPLO, 'Non unit', N, N, A, IA, JA,
     $                      DESCA, WORK )
      ELSE IF( LSAMEN( 2, MATTYP( 2:3 ), 'PD' ) ) THEN
*
         AFORM = 'H'
         DIAG = 'D'
         AUXNORM = PCLANHE( '1', UPLO, N, A, IA, JA, DESCA, WORK )
*
      END IF
      RCOND   = ANORM*AUXNORM
*
*     Compute inv(A)*A
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              ICURROW, ICURCOL )
*
*     Define array descriptor for working array WORK
*
      IROFF = MOD( IA-1, DESCA( MB_ ) )
      NP = NUMROC( N+IROFF, DESCA( MB_ ), MYROW, ICURROW, NPROW )
      CALL DESCSET( DESCW, N+IROFF, DESCA( NB_ ), DESCA( MB_ ),
     $              DESCA( NB_ ), ICURROW, ICURCOL, DESCA( CTXT_ ),
     $              MAX( 1, NP ) )
      IPW = DESCW( LLD_ ) * DESCW( NB_ ) + 1
*
      IF( MYROW.EQ.ICURROW ) THEN
         II = IROFF + 1
         NP = NP - IROFF
      ELSE
         II = 1
      END IF
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      JB = JN - JA + 1
*
*     Handle first block separately, regenerate a block of columns of A
*
      IW = IROFF + 1
      IF( MYCOL.EQ.ICURCOL ) THEN
         IF( LSAMEN( 2, MATTYP( 2:3 ), 'TR' ) ) THEN
            CALL PCMATGEN( ICTXT, AFORM, DIAG, DESCA( M_ ), DESCA( N_ ),
     $                     DESCW( MB_ ), DESCW( NB_ ), WORK,
     $                     DESCW( LLD_ ), DESCA( RSRC_ ),
     $                     DESCA( CSRC_ ), IASEED, IIA-1, NP,
     $                     JJA-1, JB, MYROW, MYCOL, NPROW, NPCOL )
            IF( LSAMEN( 3, MATTYP, 'UTR' ) ) THEN
               CALL PCLASET( 'Lower', N-1, JB, ZERO, ZERO, WORK, IW+1,
     $                       1, DESCW )
            ELSE
               CALL PCLASET( 'Upper', JB-1, JB-1, ZERO, ZERO, WORK, IW,
     $                       2, DESCW )
            END IF
         ELSE
            CALL PCMATGEN( ICTXT, AFORM, DIAG, DESCA( M_ ), DESCA( N_ ),
     $                     DESCW( MB_ ), DESCW( NB_ ), WORK( IPW ),
     $                     DESCW( LLD_ ), DESCA( RSRC_ ),
     $                     DESCA( CSRC_ ), IASEED,
     $                     IIA-1, NP, JJA-1, JB, MYROW, MYCOL, NPROW,
     $                     NPCOL )
         END IF
      END IF
*
*     Multiply A^{-1}*A
*
      IF( LSAMEN( 3, MATTYP, 'GEN' ) ) THEN
*
         CALL PCGEMM( 'No tranpose', 'No transpose', N, JB, N, ONE, A,
     $                IA, JA, DESCA, WORK( IPW ), IW, 1, DESCW, ZERO,
     $                WORK, IW, 1, DESCW )
*
      ELSE IF( LSAMEN( 2, MATTYP( 2:3 ), 'TR' ) ) THEN
*
         CALL PCTRMM( 'Left', UPLO, 'No tranpose', 'Non unit', N, JB,
     $                ONE, A, IA, JA, DESCA, WORK, IW, 1, DESCW )
*
      ELSE IF( LSAMEN( 2, MATTYP( 2:3 ), 'PD' ) ) THEN
*
         CALL PCHEMM( 'Left', UPLO, N, JB, ONE, A, IA, JA, DESCA,
     $                WORK(IPW), IW, 1, DESCW, ZERO, WORK, IW, 1,
     $                DESCW )
*
      END IF
*
*     subtract the identity matrix to the diagonal block of these cols
*
      IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
         DO 10 KK = 0, JB-1
            WORK( II+KK*(DESCW(LLD_)+1) ) =
     $                 WORK( II+KK*(DESCW( LLD_ )+1) )-ONE
   10    CONTINUE
      END IF
*
      NRMINVAXA = PCLANGE( '1', N, JB, WORK, IW, 1, DESCW, WORK( IPW ) )
*
      IF( MYROW.EQ.ICURROW )
     $   II = II + JB
      IF( MYCOL.EQ.ICURCOL )
     $   JJA = JJA + JB
      ICURROW = MOD( ICURROW+1, NPROW )
      ICURCOL = MOD( ICURCOL+1, NPCOL )
      DESCW( CSRC_ ) = ICURCOL
*
      DO 30 J = JN+1, JA+N-1, DESCA( NB_ )
*
         JB = MIN( N-J+JA, DESCA( NB_ ) )
*
*        regenerate a block of columns of A
*
         IF( MYCOL.EQ.ICURCOL ) THEN
            IF( LSAMEN( 2, MATTYP( 2:3 ), 'TR' ) ) THEN
               CALL PCMATGEN( ICTXT, AFORM, DIAG, DESCA( M_ ),
     $                        DESCA( N_ ), DESCW( MB_ ), DESCW( NB_ ),
     $                        WORK, DESCW( LLD_ ), DESCA( RSRC_ ),
     $                        DESCA( CSRC_ ),
     $                        IASEED, IIA-1, NP, JJA-1, JB, MYROW,
     $                        MYCOL, NPROW, NPCOL )
               IF( LSAMEN( 3, MATTYP, 'UTR' ) ) THEN
                  CALL PCLASET( 'Lower', JA+N-J-1, JB, ZERO, ZERO,
     $                       WORK, IW+J-JA+1, 1, DESCW )
               ELSE
                  CALL PCLASET( 'All', J-JA, JB, ZERO, ZERO, WORK, IW,
     $                          1, DESCW )
                  CALL PCLASET( 'Upper', JB-1, JB-1, ZERO, ZERO,
     $                          WORK, IW+J-JA, 2, DESCW )
               END IF
            ELSE
               CALL PCMATGEN( ICTXT, AFORM, DIAG, DESCA( M_ ),
     $                        DESCA( N_ ), DESCW( MB_ ), DESCW( NB_ ),
     $                        WORK( IPW ), DESCW( LLD_ ),
     $                        DESCA( RSRC_ ), DESCA( CSRC_ ), IASEED,
     $                        IIA-1, NP,
     $                        JJA-1, JB, MYROW, MYCOL, NPROW, NPCOL )
            END IF
         END IF
*
*        Multiply A^{-1}*A
*
         IF( LSAMEN( 3, MATTYP, 'GEN' ) ) THEN
*
            CALL PCGEMM( 'No tranpose', 'No transpose', N, JB, N, ONE,
     $                   A, IA, JA, DESCA, WORK( IPW ), IW, 1, DESCW,
     $                   ZERO, WORK, IW, 1, DESCW )
*
         ELSE IF( LSAMEN( 2, MATTYP(2:3), 'TR' ) ) THEN
*
            CALL PCTRMM( 'Left', UPLO, 'No tranpose', 'Non unit', N, JB,
     $                   ONE, A, IA, JA, DESCA, WORK, IW, 1, DESCW )
*
         ELSE IF( LSAMEN( 2, MATTYP( 2:3 ), 'PD' ) ) THEN
*
            CALL PCHEMM( 'Left', UPLO, N, JB, ONE, A, IA, JA, DESCA,
     $                   WORK(IPW), IW, 1, DESCW, ZERO, WORK, IW, 1,
     $                   DESCW )
*
         END IF
*
*        subtract the identity matrix to the diagonal block of these
*        cols
*
         IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
            DO 20 KK = 0, JB-1
               WORK( II+KK*(DESCW( LLD_ )+1) ) =
     $                   WORK( II+KK*(DESCW( LLD_ )+1) ) - ONE
   20       CONTINUE
         END IF
*
*        Compute the 1-norm of these JB cols
*
         TEMP = PCLANGE( '1', N, JB, WORK, IW, 1, DESCW, WORK( IPW ) )
         NRMINVAXA = MAX( TEMP, NRMINVAXA )
*
         IF( MYROW.EQ.ICURROW )
     $      II = II + JB
         IF( MYCOL.EQ.ICURCOL )
     $      JJA = JJA + JB
         ICURROW = MOD( ICURROW+1, NPROW )
         ICURCOL = MOD( ICURCOL+1, NPCOL )
         DESCW( CSRC_ ) = ICURCOL
*
   30 CONTINUE
*
*     Compute the scaled residual
*
      FRESID = NRMINVAXA / ( N * EPS * ANORM )
*
      RETURN
*
*     End of PCINVCHK
*
      END
