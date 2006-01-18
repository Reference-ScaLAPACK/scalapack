      SUBROUTINE PSSTEDC( COMPZ, N, D, E, Q, IQ, JQ, DESCQ, WORK, LWORK,
     $                    IWORK, LIWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 13, 2000
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, IQ, JQ, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCQ( * ), IWORK( * )
      REAL               D( * ), E( * ), Q( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*  PSSTEDC computes all eigenvalues and eigenvectors of a
*  symmetric tridiagonal matrix in parallel, using the divide and
*  conquer algorithm.
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.  See SLAED3 for details.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.    (NOT IMPLEMENTED YET)
*          = 'I':  Compute eigenvectors of tridiagonal matrix also.
*          = 'V':  Compute eigenvectors of original dense symmetric
*                  matrix also.  On entry, Z contains the orthogonal
*                  matrix used to reduce the original matrix to
*                  tridiagonal form.            (NOT IMPLEMENTED YET)
*
*  N       (global input) INTEGER
*          The order of the tridiagonal matrix T.  N >= 0.
*
*  D       (global input/output) REAL array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in descending order.
*
*  E       (global input/output) REAL array, dimension (N-1)
*          On entry, the subdiagonal elements of the tridiagonal matrix.
*          On exit, E has been destroyed.
*
*  Q       (local output) REAL array,
*          local dimension ( LLD_Q, LOCc(JQ+N-1))
*          Q  contains the orthonormal eigenvectors of the symmetric
*          tridiagonal matrix.
*          On output, Q is distributed across the P processes in block
*          cyclic format.
*
*  IQ      (global input) INTEGER
*          Q's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JQ      (global input) INTEGER
*          Q's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*
*  WORK    (local workspace/output) REAL array,
*          dimension (LWORK)
*          On output, WORK(1) returns the workspace needed.
*
*  LWORK   (local input/output) INTEGER,
*          the dimension of the array WORK.
*          LWORK = 6*N + 2*NP*NQ
*          NP = NUMROC( N, NB, MYROW, DESCQ( RSRC_ ), NPROW )
*          NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
*
*          If LWORK = -1, the LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          size for the WORK array.  The required workspace is returned
*          as the first element of WORK and no error message is issued
*          by PXERBLA.
*
*  IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          LIWORK = 2 + 7*N + 8*NPCOL
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  The algorithm failed to compute the INFO/(N+1) th
*                eigenvalue while working on the submatrix lying in
*                global rows and columns mod(INFO,N+1).
*
*  Further Details
*  ======= =======
*
*  Contributed by Francoise Tisseur, University of Manchester.
*
*  Reference:  F. Tisseur and J. Dongarra, "A Parallel Divide and
*              Conquer Algorithm for the Symmetric Eigenvalue Problem
*              on Distributed Memory Architectures",
*              SIAM J. Sci. Comput., 6:20 (1999), pp. 2223--2236.
*              (see also LAPACK Working Note 132)
*                http://www.netlib.org/lapack/lawns/lawn132.ps
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            ICOFFQ, IIQ, IPQ, IQCOL, IQROW, IROFFQ, JJQ,
     $                   LDQ, LIWMIN, LWMIN, MYCOL, MYROW, NB, NP,
     $                   NPCOL, NPROW, NQ
      REAL               ORGNRM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      REAL               SLANST
      EXTERNAL           INDXG2P, LSAME, NUMROC, SLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, INFOG2L, PSLAED0,
     $                   PSLASRT, PXERBLA, SLASCL, SSTEDC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD, REAL
*     ..
*     .. Executable Statements ..
*
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Test the input parameters.
*
      CALL BLACS_GRIDINFO( DESCQ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      LDQ = DESCQ( LLD_ )
      NB = DESCQ( NB_ )
      NP = NUMROC( N, NB, MYROW, DESCQ( RSRC_ ), NPROW )
      NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 600+CTXT_ )
      ELSE
         CALL CHK1MAT( N, 2, N, 2, IQ, JQ, DESCQ, 8, INFO )
         IF( INFO.EQ.0 ) THEN
            NB = DESCQ( NB_ )
            IROFFQ = MOD( IQ-1, DESCQ( MB_ ) )
            ICOFFQ = MOD( JQ-1, DESCQ( NB_ ) )
            IQROW = INDXG2P( IQ, NB, MYROW, DESCQ( RSRC_ ), NPROW )
            IQCOL = INDXG2P( JQ, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
            LWMIN = 6*N + 2*NP*NQ
            LIWMIN = 2 + 7*N + 8*NPCOL
            WORK( 1 ) = REAL( LWMIN )
            IWORK( 1 ) = LIWMIN
            LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
            IF( .NOT.LSAME( COMPZ, 'I' ) ) THEN
               INFO = -1
            ELSE IF( N.LT.0 ) THEN
               INFO = -2
            ELSE IF( IROFFQ.NE.ICOFFQ .OR. ICOFFQ.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCQ( MB_ ).NE.DESCQ( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -10
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -12
            END IF
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCQ( CTXT_ ), 'PSSTEDC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return
*
      IF( N.EQ.0 )
     $   GO TO 10
      CALL INFOG2L( IQ, JQ, DESCQ, NPROW, NPCOL, MYROW, MYCOL, IIQ, JJQ,
     $              IQROW, IQCOL )
      IF( N.EQ.1 ) THEN
         IF( MYROW.EQ.IQROW .AND. MYCOL.EQ.IQCOL )
     $      Q( 1 ) = ONE
         GO TO 10
      END IF
*
*     If N is smaller than the minimum divide size NB, then
*     solve the problem with the serial divide and conquer
*     code locally.
*
      IF( N.LE.NB ) THEN
         IF( ( MYROW.EQ.IQROW ) .AND. ( MYCOL.EQ.IQCOL ) ) THEN
            IPQ = IIQ + ( JJQ-1 )*LDQ
            CALL SSTEDC( 'I', N, D, E, Q( IPQ ), LDQ, WORK, LWORK,
     $                   IWORK, LIWORK, INFO )
            IF( INFO.NE.0 ) THEN
               INFO = ( N+1 ) + N
               GO TO 10
            END IF
         END IF
         GO TO 10
      END IF
*
*     If P=NPROW*NPCOL=1, solve the problem with SSTEDC.
*
      IF( NPCOL*NPROW.EQ.1 ) THEN
         IPQ = IIQ + ( JJQ-1 )*LDQ
         CALL SSTEDC( 'I', N, D, E, Q( IPQ ), LDQ, WORK, LWORK, IWORK,
     $                LIWORK, INFO )
         GO TO 10
      END IF
*
*     Scale matrix to allowable range, if necessary.
*
      ORGNRM = SLANST( 'M', N, D, E )
      IF( ORGNRM.NE.ZERO ) THEN
         CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, INFO )
         CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, N-1, 1, E, N-1, INFO )
      END IF
*
      CALL PSLAED0( N, D, E, Q, IQ, JQ, DESCQ, WORK, IWORK, INFO )
*
*     Sort eigenvalues and corresponding eigenvectors
*
      CALL PSLASRT( 'I', N, D, Q, IQ, JQ, DESCQ, WORK, LWORK, IWORK,
     $              LIWORK, INFO )
*
*           Scale back.
*
      IF( ORGNRM.NE.ZERO )
     $   CALL SLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, INFO )
*
   10 CONTINUE
*
      IF( LWORK.GT.0 )
     $   WORK( 1 ) = REAL( LWMIN )
      IF( LIWORK.GT.0 )
     $   IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of PSSTEDC
*
      END
