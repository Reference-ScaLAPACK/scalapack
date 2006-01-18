      SUBROUTINE PCPOEQU( N, A, IA, JA, DESCA, SR, SC, SCOND, AMAX,
     $                    INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, N
      REAL               AMAX, SCOND
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      REAL               SC( * ), SR( * )
      COMPLEX            A( * )
*     ..
*
*  Purpose
*  =======
*
*  PCPOEQU computes row and column scalings intended to
*  equilibrate a distributed Hermitian positive definite matrix
*  sub( A ) = A(IA:IA+N-1,JA:JA+N-1) and reduce its condition number
*  (with respect to the two-norm).  SR and SC contain the scale
*  factors, S(i) = 1/sqrt(A(i,i)), chosen so that the scaled distri-
*  buted matrix B with elements B(i,j) = S(i)*A(i,j)*S(j) has ones on
*  the  diagonal.  This choice of SR and SC puts the condition number
*  of B within a factor N of the smallest possible condition number
*  over all possible diagonal scalings.
*
*  The scaling factor are stored along process rows in SR and along
*  process columns in SC. The duplication of information simplifies
*  greatly the application of the factors.
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
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on i.e the
*          order of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input) COMPLEX pointer into the local memory to an
*          array of local dimension ( LLD_A, LOCc(JA+N-1) ), the
*          N-by-N Hermitian positive definite distributed matrix
*          sub( A ) whose scaling factors are to be computed.  Only the
*          diagonal elements of sub( A ) are referenced.
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
*  SR      (local output) REAL array, dimension LOCr(M_A)
*          If INFO = 0, SR(IA:IA+N-1) contains the row scale factors
*          for sub( A ). SR is aligned with the distributed matrix A,
*          and replicated across every process column. SR is tied to the
*          distributed matrix A.
*
*  SC      (local output) REAL array, dimension LOCc(N_A)
*          If INFO = 0, SC(JA:JA+N-1) contains the column scale factors
*          for A(IA:IA+M-1,JA:JA+N-1). SC is aligned with the distribu-
*          ted matrix A, and replicated down every process row. SC is
*          tied to the distributed matrix A.
*
*  SCOND   (global output) REAL
*          If INFO = 0, SCOND contains the ratio of the smallest SR(i)
*          (or SC(j)) to the largest SR(i) (or SC(j)), with
*          IA <= i <= IA+N-1 and JA <= j <= JA+N-1. If SCOND >= 0.1
*          and AMAX is neither too large nor too small, it is not worth
*          scaling by SR (or SC).
*
*  AMAX    (global output) REAL
*          Absolute value of largest matrix element.  If AMAX is very
*          close to overflow or very close to underflow, the matrix
*          should be scaled.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = K, the K-th diagonal entry of sub( A ) is
*                nonpositive.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          ALLCTOP, COLCTOP, ROWCTOP
      INTEGER            IACOL, IAROW, ICOFF, ICTXT, ICURCOL, ICURROW,
     $                   IDUMM, II, IIA, IOFFA, IOFFD, IROFF, J, JB, JJ,
     $                   JJA, JN, LDA, LL, MYCOL, MYROW, NP, NPCOL,
     $                   NPROW, NQ
      REAL               AII, SMIN
*     ..
*     .. Local Arrays ..
      INTEGER            DESCSC( DLEN_ ), DESCSR( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, IGAMN2D,
     $                   INFOG2L, PCHK1MAT, PB_TOPGET, PXERBLA,
     $                   SGAMN2D, SGAMX2D, SGSUM2D
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, NUMROC
      REAL               PSLAMCH
      EXTERNAL           ICEIL, NUMROC, PSLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters.
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(500+CTXT_)
      ELSE
         CALL CHK1MAT( N, 1, N, 1, IA, JA, DESCA, 5, INFO )
         CALL PCHK1MAT( N, 1, N, 1, IA, JA, DESCA, 5, 0, IDUMM, IDUMM,
     $                  INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCPOEQU', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         SCOND = ONE
         AMAX = ZERO
         RETURN
      END IF
*
      CALL PB_TOPGET( ICTXT, 'Combine', 'All', ALLCTOP )
      CALL PB_TOPGET( ICTXT, 'Combine', 'Rowwise', ROWCTOP )
      CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
*
*     Compute some local indexes
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      IROFF = MOD( IA-1, DESCA( MB_ ) )
      ICOFF = MOD( JA-1, DESCA( NB_ ) )
      NP = NUMROC( N+IROFF, DESCA( MB_ ), MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFF, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IF( MYROW.EQ.IAROW )
     $   NP = NP - IROFF
      IF( MYCOL.EQ.IACOL )
     $   NQ = NQ - ICOFF
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      LDA = DESCA( LLD_ )
*
*     Assign descriptors for SR and SC arrays
*
      CALL DESCSET( DESCSR, N, 1, DESCA( MB_ ), 1, 0, 0, ICTXT,
     $               MAX( 1, NP ) )
      CALL DESCSET( DESCSC, 1, N, 1, DESCA( NB_ ), 0, 0, ICTXT, 1 )
*
*     Initialize the scaling factors to zero.
*
      DO 10 II = IIA, IIA+NP-1
         SR( II ) = ZERO
   10 CONTINUE
*
      DO 20 JJ = JJA, JJA+NQ-1
         SC( JJ ) = ZERO
   20 CONTINUE
*
*     Find the minimum and maximum diagonal elements.
*     Handle first block separately.
*
      II = IIA
      JJ = JJA
      JB = JN-JA+1
      SMIN = ONE / PSLAMCH( ICTXT, 'S' )
      AMAX = ZERO
*
      IOFFA = II+(JJ-1)*LDA
      IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL ) THEN
         IOFFD = IOFFA
         DO 30 LL = 0, JB-1
            AII = REAL( A( IOFFD ) )
            SR( II+LL ) = AII
            SC( JJ+LL ) = AII
            SMIN = MIN( SMIN, AII )
            AMAX = MAX( AMAX, AII )
            IF( AII.LE.ZERO .AND. INFO.EQ.0 )
     $         INFO = LL + 1
            IOFFD = IOFFD + LDA + 1
   30    CONTINUE
      END IF
*
      IF( MYROW.EQ.IAROW ) THEN
         II = II + JB
         IOFFA = IOFFA + JB
      END IF
      IF( MYCOL.EQ.IACOL ) THEN
         JJ = JJ + JB
         IOFFA = IOFFA + JB*LDA
      END IF
      ICURROW = MOD( IAROW+1, NPROW )
      ICURCOL = MOD( IACOL+1, NPCOL )
*
*     Loop over remaining blocks of columns
*
      DO 50 J = JN+1, JA+N-1, DESCA( NB_ )
         JB = MIN( N-J+JA, DESCA( NB_ ) )
*
         IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
            IOFFD = IOFFA
            DO 40 LL = 0, JB-1
               AII = REAL( A( IOFFD ) )
               SR( II+LL ) = AII
               SC( JJ+LL ) = AII
               SMIN = MIN( SMIN, AII )
               AMAX = MAX( AMAX, AII )
               IF( AII.LE.ZERO .AND. INFO.EQ.0 )
     $            INFO = J + LL - JA + 1
               IOFFD = IOFFD + LDA + 1
   40       CONTINUE
         END IF
*
         IF( MYROW.EQ.ICURROW ) THEN
            II = II + JB
            IOFFA = IOFFA + JB
         END IF
         IF( MYCOL.EQ.ICURCOL ) THEN
            JJ = JJ + JB
            IOFFA = IOFFA + JB*LDA
         END IF
         ICURROW = MOD( ICURROW+1, NPROW )
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
   50 CONTINUE
*
*     Compute scaling factors
*
      CALL SGSUM2D( ICTXT, 'Columnwise', COLCTOP, 1, NQ, SC( JJA ),
     $              1, -1, MYCOL )
      CALL SGSUM2D( ICTXT, 'Rowwise', ROWCTOP, NP, 1, SR( IIA ),
     $              MAX( 1, NP ), -1, MYCOL )
*
      CALL SGAMX2D( ICTXT, 'All', ALLCTOP, 1, 1, AMAX, 1, IDUMM, IDUMM,
     $              -1, -1, MYCOL )
      CALL SGAMN2D( ICTXT, 'All', ALLCTOP, 1, 1, SMIN, 1, IDUMM, IDUMM,
     $              -1, -1, MYCOL )
*
      IF( SMIN.LE.ZERO ) THEN
*
*        Find the first non-positive diagonal element and return.
*
         CALL IGAMN2D( ICTXT, 'All', ALLCTOP, 1, 1, INFO, 1, II, JJ, -1,
     $                 -1, MYCOL )
         RETURN
*
      ELSE
*
*        Set the scale factors to the reciprocals
*        of the diagonal elements.
*
         DO 60 II = IIA, IIA+NP-1
            SR( II ) = ONE / SQRT( SR( II ) )
   60    CONTINUE
*
         DO 70 JJ = JJA, JJA+NQ-1
            SC( JJ ) = ONE / SQRT( SC( JJ ) )
   70    CONTINUE
*
*        Compute SCOND = min(S(I)) / max(S(I))
*
         SCOND = SQRT( SMIN ) / SQRT( AMAX )
*
      END IF
*
      RETURN
*
*     End of PCPOEQU
*
      END
