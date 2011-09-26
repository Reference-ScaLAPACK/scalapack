      SUBROUTINE PDSYEV( JOBZ, UPLO, N, A, IA, JA, DESCA, W,
     $                   Z, IZ, JZ, DESCZ, WORK, LWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            IA, INFO, IZ, JA, JZ, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * )
      DOUBLE PRECISION   A( * ), W( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PDSYEV computes all eigenvalues and, optionally, eigenvectors
*  of a real symmetric matrix A by calling the recommended sequence
*  of ScaLAPACK routines.
*
*  In its present form, PDSYEV assumes a homogeneous system and makes
*  no checks for consistency of the eigenvalues or eigenvectors across
*  the different processes.  Because of this, it is possible that a
*  heterogeneous system may return incorrect results without any error
*  messages.
*
*  Notes
*  =====
*  A description vector is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This vector stores the information required to
*  establish the mapping between a matrix entry and its corresponding
*  process and memory location.
*
*  In the following comments, the character _ should be read as
*  "of the distributed matrix".  Let A be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_) The descriptor type.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the distributed
*                                 matrix A.
*  N_A    (global) DESCA( N_ )    The number of columns in the distri-
*                                 buted matrix A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of A.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of A.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the matrix A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of A is distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array storing the local blocks of the
*                                 distributed matrix A.
*                                 LLD_A >= MAX(1,LOCr(M_A)).
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
*
*  Arguments
*  =========
*
*     NP = the number of rows local to a given process.
*     NQ = the number of columns local to a given process.
*
*  JOBZ    (global input) CHARACTER*1
*          Specifies whether or not to compute the eigenvectors:
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns of the matrix A.  N >= 0.
*
*  A       (local input/workspace) block cyclic DOUBLE PRECISION array,
*          global dimension (N, N), local dimension ( LLD_A,
*          LOCc(JA+N-1) )
*
*          On entry, the symmetric matrix A.  If UPLO = 'U', only the
*          upper triangular part of A is used to define the elements of
*          the symmetric matrix.  If UPLO = 'L', only the lower
*          triangular part of A is used to define the elements of the
*          symmetric matrix.
*
*          On exit, the lower triangle (if UPLO='L') or the upper
*          triangle (if UPLO='U') of A, including the diagonal, is
*          destroyed.
*
*  IA      (global input) INTEGER
*          A's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*          If DESCA( CTXT_ ) is incorrect, PDSYEV cannot guarantee
*          correct error reporting.
*
*  W       (global output) DOUBLE PRECISION array, dimension (N)
*          If INFO=0, the eigenvalues in ascending order.
*
*  Z       (local output) DOUBLE PRECISION array,
*          global dimension (N, N),
*          local dimension ( LLD_Z, LOCc(JZ+N-1) )
*          If JOBZ = 'V', then on normal exit the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix
*          corresponding to the selected eigenvalues.
*          If JOBZ = 'N', then Z is not referenced.
*
*  IZ      (global input) INTEGER
*          Z's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JZ      (global input) INTEGER
*          Z's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*          DESCZ( CTXT_ ) must equal DESCA( CTXT_ )
*
*  WORK    (local workspace/output) DOUBLE PRECISION array,
*          dimension (LWORK)
*          Version 1.0:  on output, WORK(1) returns the workspace
*          needed to guarantee completion.
*          If the input parameters are incorrect, WORK(1) may also be
*          incorrect.
*
*          If JOBZ='N' WORK(1) = minimal=optimal amount of workspace
*          If JOBZ='V' WORK(1) = minimal workspace required to
*             generate all the eigenvectors.
*
*
*  LWORK   (local input) INTEGER
*          See below for definitions of variables used to define LWORK.
*          If no eigenvectors are requested (JOBZ = 'N') then
*             LWORK >= 5*N + SIZESYTRD + 1
*          where
*             SIZESYTRD = The workspace requirement for PDSYTRD
*                         and is MAX( NB * ( NP +1 ), 3 * NB )
*          If eigenvectors are requested (JOBZ = 'V' ) then
*             the amount of workspace required to guarantee that all
*             eigenvectors are computed is:
*
*             QRMEM = 2*N-2
*             LWMIN = 5*N + N*LDC + MAX( SIZEMQRLEFT, QRMEM ) + 1
*
*          Variable definitions:
*             NB = DESCA( MB_ ) = DESCA( NB_ ) =
*                  DESCZ( MB_ ) = DESCZ( NB_ )
*             NN = MAX( N, NB, 2 )
*             DESCA( RSRC_ ) = DESCA( RSRC_ ) = DESCZ( RSRC_ ) =
*                              DESCZ( CSRC_ ) = 0
*             NP = NUMROC( NN, NB, 0, 0, NPROW )
*             NQ = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPCOL )
*             NRC = NUMROC( N, NB, MYPROWC, 0, NPROCS)
*             LDC = MAX( 1, NRC )
*             SIZEMQRLEFT = The workspace requirement for PDORMTR
*                           when it's SIDE argument is 'L'.
*
*          With MYPROWC defined when a new context is created as:
*             CALL BLACS_GET( DESCA( CTXT_ ), 0, CONTEXTC )
*             CALL BLACS_GRIDINIT( CONTEXTC, 'R', NPROCS, 1 )
*             CALL BLACS_GRIDINFO( CONTEXTC, NPROWC, NPCOLC, MYPROWC,
*                                  MYPCOLC )
*
*          If LWORK = -1, the LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          size for the WORK array.  The required workspace is returned
*          as the first element of WORK and no error message is issued
*          by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = 1 through N, the i(th) eigenvalue did not
*                converge in DSTEQR2 after a total of 30*N iterations.
*                If INFO = N+1, then PDSYEV has detected heterogeneity
*                by finding that eigenvalues were not identical across
*                the process grid.  In this case, the accuracy of
*                the results from PDSYEV cannot be guaranteed.
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrices A(IA:*, JA:*) and Z(IZ:IZ+M-1,JZ:JZ+N-1)
*  must verify some alignment properties, namely the following
*  expressions should be true:
*
*  ( MB_A.EQ.NB_A.EQ.MB_Z .AND. IROFFA.EQ.IROFFZ .AND. IROFFA.EQ.0 .AND.
*    IAROW.EQ.IZROW )
*  where
*  IROFFA = MOD( IA-1, MB_A ) and ICOFFA = MOD( JA-1, NB_A ).
*
*  =====================================================================
*
*  Version 1.4 limitations:
*     DESCA(MB_) = DESCA(NB_)
*     DESCA(M_) = DESCZ(M_)
*     DESCA(N_) = DESCZ(N_)
*     DESCA(MB_) = DESCZ(MB_)
*     DESCA(NB_) = DESCZ(NB_)
*     DESCA(RSRC_) = DESCZ(RSRC_)
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   FIVE, ONE, TEN, ZERO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0,
     $                     TEN = 10.0D+0, FIVE = 5.0D+0 )
      INTEGER            IERREIN, IERRCLS, IERRSPC, IERREBZ, ITHVAL
      PARAMETER          ( IERREIN = 1, IERRCLS = 2, IERRSPC = 4,
     $                   IERREBZ = 8, ITHVAL = 10 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, WANTZ
      INTEGER            CONTEXTC, CSRC_A, I, IACOL, IAROW, ICOFFA, 
     $                   IINFO, INDD, INDD2, INDE, INDE2, INDTAU, 
     $                   INDWORK, INDWORK2, IROFFA, IROFFZ, ISCALE, 
     $                   IZROW, J, K, LDC, LLWORK, LWMIN, MB_A, MB_Z, 
     $                   MYCOL, MYPCOLC, MYPROWC, MYROW, NB, NB_A, NB_Z,
     $                   NP, NPCOL, NPCOLC, NPROCS, NPROW, NPROWC, NQ, 
     $                   NRC, QRMEM, RSRC_A, RSRC_Z, SIZEMQRLEFT, 
     $                   SIZESYTRD
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, 
     $                   SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER            DESCQR( 9 ), IDUM1( 3 ), IDUM2( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC, SL_GRIDRESHAPE
      DOUBLE PRECISION   PDLAMCH, PDLANSY
      EXTERNAL           LSAME, NUMROC, PDLAMCH, PDLANSY,
     $                   SL_GRIDRESHAPE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDEXIT, BLACS_GRIDINFO, CHK1MAT, DCOPY,
     $                   DESCINIT, DSCAL, DSTEQR2, PCHK1MAT, PCHK2MAT,
     $                   PDELGET, PDGEMR2D, PDLASCL, PDLASET, PDORMTR,
     $                   PDSYTRD, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, ICHAR, MAX, MIN, MOD, SQRT, INT
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Quick return
*
      IF( N.EQ.0 ) RETURN
*
*     Test the input arguments.
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
*
      WANTZ = LSAME( JOBZ, 'V' )
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 700+CTXT_ )
      ELSE IF( WANTZ ) THEN
         IF( DESCA( CTXT_ ).NE.DESCZ( CTXT_ ) ) THEN
            INFO = -( 1200+CTXT_ )
         END IF
      END IF
      IF( INFO .EQ. 0 ) THEN
         CALL CHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, INFO )
         IF( WANTZ )
     $      CALL CHK1MAT( N, 3, N, 3, IZ, JZ, DESCZ, 12, INFO )
*
         IF( INFO.EQ.0 ) THEN
*
*           Get machine constants.
*
            SAFMIN = PDLAMCH( DESCA( CTXT_ ), 'Safe minimum' )
            EPS = PDLAMCH( DESCA( CTXT_ ), 'Precision' )
            SMLNUM = SAFMIN / EPS
            BIGNUM = ONE / SMLNUM
            RMIN = SQRT( SMLNUM )
            RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
            NPROCS = NPROW*NPCOL
            NB_A = DESCA( NB_ )
            MB_A = DESCA( MB_ )
            NB = NB_A
            LOWER = LSAME( UPLO, 'L' )
*
            RSRC_A = DESCA( RSRC_ )
            CSRC_A = DESCA( CSRC_ )
            IROFFA = MOD( IA-1, MB_A )
            ICOFFA = MOD( JA-1, NB_A )
            IAROW = INDXG2P( 1, NB_A, MYROW, RSRC_A, NPROW )
            IACOL = INDXG2P( 1, MB_A, MYCOL, CSRC_A, NPCOL )
            NP = NUMROC( N+IROFFA, NB, MYROW, IAROW, NPROW )
            NQ = NUMROC( N+ICOFFA, NB, MYCOL, IACOL, NPCOL )

            IF( WANTZ ) THEN
               NB_Z = DESCZ( NB_ )
               MB_Z = DESCZ( MB_ )
               RSRC_Z = DESCZ( RSRC_ )
               IROFFZ = MOD( IZ-1, MB_A )
               IZROW = INDXG2P( 1, NB_A, MYROW, RSRC_Z, NPROW )
               SIZEMQRLEFT = MAX( ( NB_A*( NB_A-1 ) ) / 2, ( NP+NQ )*
     $                            NB_A ) + NB_A*NB_A
            ELSE
               SIZEMQRLEFT = 0
               IROFFZ = 0
               IZROW = 0
            END IF
            SIZESYTRD = MAX( NB * ( NP +1 ), 3 * NB )
*
*           Initialize the context of the single column distributed
*           matrix required by DSTEQR2. This specific distribution
*           allows each process to do 1/pth of the work updating matrix
*           Q during DSTEQR2 and achieve some parallelization to an
*           otherwise serial subroutine.
*
            LDC = 0
            IF( WANTZ ) THEN
               CONTEXTC = SL_GRIDRESHAPE( DESCA( CTXT_ ), 0, 1, 1,
     $                                    NPROCS, 1 )
               CALL BLACS_GRIDINFO( CONTEXTC, NPROWC, NPCOLC, MYPROWC,
     $                              MYPCOLC )
               NRC = NUMROC( N, NB_A, MYPROWC, 0, NPROCS)
               LDC = MAX( 1, NRC )
               CALL DESCINIT( DESCQR, N, N, NB, NB, 0, 0, CONTEXTC,
     $                        LDC, INFO )
            END IF

*
*           Set up pointers into the WORK array
*
            INDTAU = 1
            INDE = INDTAU + N
            INDD = INDE + N
            INDD2 = INDD + N
            INDE2 = INDD2 + N
            INDWORK = INDE2 + N
            INDWORK2 = INDWORK + N*LDC
            LLWORK = LWORK - INDWORK + 1
*
*           Compute the total amount of space needed
*
            QRMEM = 2*N-2
            IF( WANTZ ) THEN
               LWMIN = 5*N + N*LDC + MAX( SIZEMQRLEFT, QRMEM ) + 1
            ELSE
               LWMIN = 5*N + SIZESYTRD + 1
            END IF
*
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
               INFO = -1
            ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
               INFO = -2
            ELSE IF( LWORK.LT.LWMIN .AND. LWORK.NE.-1 ) THEN
               INFO = -14
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 700+NB_ )
            END IF
            IF( WANTZ ) THEN
               IF( IROFFA.NE.IROFFZ ) THEN
                  INFO = -10
               ELSE IF( IAROW.NE.IZROW ) THEN
                  INFO = -10
               ELSE IF( DESCA( M_ ).NE.DESCZ( M_ ) ) THEN
                  INFO = -( 1200+M_ )
               ELSE IF( DESCA( N_ ).NE.DESCZ( N_ ) ) THEN
                  INFO = -( 1200+N_ )
               ELSE IF( DESCA( MB_ ).NE.DESCZ( MB_ ) ) THEN
                  INFO = -( 1200+MB_ )
               ELSE IF( DESCA( NB_ ).NE.DESCZ( NB_ ) ) THEN
                  INFO = -( 1200+NB_ )
               ELSE IF( DESCA( RSRC_ ).NE.DESCZ( RSRC_ ) ) THEN
                  INFO = -( 1200+RSRC_ )
               ELSE IF( DESCA( CTXT_ ).NE.DESCZ( CTXT_ ) ) THEN
                  INFO = -( 1200+CTXT_ )
               ENDIF
            END IF
         END IF
         IF( WANTZ ) THEN
            IDUM1( 1 ) = ICHAR( 'V' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'N' )
         END IF
         IDUM2( 1 ) = 1
         IF( LOWER ) THEN
            IDUM1( 2 ) = ICHAR( 'L' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'U' )
         END IF
         IDUM2( 2 ) = 2
         IF( LWORK.EQ.-1 ) THEN
            IDUM1( 3 ) = -1
         ELSE
            IDUM1( 3 ) = 1
         END IF
         IDUM2( 3 ) = 3
         IF( LSAME( JOBZ, 'V' ) ) THEN
            CALL PCHK2MAT( N, 3, N, 3, IA, JA, DESCA, 7, N, 3, N, 3,
     $                     IZ, JZ, DESCZ, 12, 3, IDUM1, IDUM2, INFO )
         ELSE
            CALL PCHK1MAT( N, 3, N, 3, IA, JA, DESCA, 7, 3, IDUM1,
     $                     IDUM2, INFO )
         END IF
*
*     Write the required workspace for lwork queries.
*
         WORK( 1 ) = DBLE( LWMIN )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PDSYEV', -INFO )
         IF( WANTZ ) CALL BLACS_GRIDEXIT( CONTEXTC )
         RETURN
      ELSE IF( LWORK .EQ. -1 ) THEN
         IF( WANTZ ) CALL BLACS_GRIDEXIT( CONTEXTC )
         RETURN
      END IF
*
*     Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
*
      ANRM = PDLANSY( 'M', UPLO, N, A, IA, JA, DESCA, WORK( INDWORK ) )
*

      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
*
      IF( ISCALE.EQ.1 ) THEN
         CALL PDLASCL( UPLO, ONE, SIGMA, N, N, A, IA, JA, DESCA, IINFO )
      END IF
*
*     Reduce symmetric matrix to tridiagonal form.
*
      CALL PDSYTRD( UPLO, N, A, IA, JA, DESCA, WORK( INDD ),
     $              WORK( INDE ), WORK( INDTAU ), WORK( INDWORK ),
     $              LLWORK, IINFO )
*
*     Copy the values of D, E to all processes.
*
      DO 10 I=1,N
         CALL PDELGET( 'A', ' ', WORK(INDD2+I-1), A,
     $                 I+IA-1, I+JA-1, DESCA )
 10   CONTINUE
      IF( LSAME( UPLO, 'U') ) THEN
          DO 20 I=1,N-1
             CALL PDELGET( 'A', ' ', WORK(INDE2+I-1), A, 
     $                     I+IA-1, I+JA, DESCA )
 20       CONTINUE
      ELSE
          DO 30 I=1,N-1
             CALL PDELGET( 'A', ' ', WORK(INDE2+I-1), A,
     $                     I+IA, I+JA-1, DESCA )
 30       CONTINUE
      ENDIF
*
      IF( WANTZ ) THEN
*
         CALL PDLASET( 'Full', N, N, ZERO, ONE, WORK( INDWORK ), 1, 1,
     $                 DESCQR )
*
*        DSTEQR2 is a modified version of LAPACK's DSTEQR.  The
*        modifications allow each process to perform partial updates
*        to matrix Q.
*
         CALL DSTEQR2( 'I', N, WORK( INDD2 ), WORK( INDE2 ),
     $                 WORK( INDWORK ), LDC, NRC, WORK( INDWORK2 ), 
     $                 INFO )
*
         CALL PDGEMR2D( N, N, WORK( INDWORK ), 1, 1, DESCQR, Z, IA, JA,
     $                  DESCZ, CONTEXTC )
*
         CALL PDORMTR( 'L', UPLO, 'N', N, N, A, IA, JA, DESCA,
     $                 WORK( INDTAU ), Z, IZ, JZ, DESCZ,
     $                 WORK( INDWORK ), LLWORK, IINFO )
*
      ELSE
*
         CALL DSTEQR2( 'N', N, WORK( INDD2 ), WORK( INDE2 ),
     $                 WORK( INDWORK ), 1, 1, WORK( INDWORK2 ),
     $                 INFO )
      ENDIF
*
*     Copy eigenvalues from workspace to output array
*
      CALL DCOPY( N, WORK( INDD2 ), 1, W, 1 )
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE .EQ. 1 ) THEN
         CALL DSCAL( N, ONE / SIGMA, W, 1 )
      END IF
*
*     Free up resources
*
      IF( WANTZ ) THEN
           CALL BLACS_GRIDEXIT( CONTEXTC )
      END IF
*
*     Compare every ith eigenvalue, or all if there are only a few,
*     across the process grid to check for heterogeneity.
*
      IF( N.LE.ITHVAL ) THEN
         J = N
         K = 1
      ELSE
         J = N/ITHVAL
         K = ITHVAL
      END IF
*
      DO 40 I = 1, J
         WORK( I+INDTAU ) = W( (I-1)*K+1 )
         WORK( I+INDE ) = W( (I-1)*K+1 )
 40   CONTINUE
*
      CALL DGAMN2D( DESCA( CTXT_ ), 'a', ' ', J, 1, WORK( 1+INDTAU ),
     $              J, 1, 1, -1, -1, 0 )
      CALL DGAMX2D( DESCA( CTXT_ ), 'a', ' ', J, 1, WORK( 1+INDE ),
     $              J, 1, 1, -1, -1, 0 )
*
      DO 50 I = 1, J
         IF( INFO.EQ.0 .AND. ( WORK( I+INDTAU )-WORK( I+INDE )
     $        .NE. ZERO ) )THEN 
            INFO = N+1
         END IF
 50   CONTINUE
*
      RETURN
*
*     End of PDSYEV
*
      END
