      SUBROUTINE PCHEEVR( JOBZ, RANGE, UPLO, N, A, IA, JA, 
     $                    DESCA, VL, VU, IL, IU, M, NZ, W, Z, IZ,
     $                    JZ, DESCZ, 
     $                    WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK,
     $                    INFO )

      IMPLICIT NONE
*
*  -- ScaLAPACK routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO

      INTEGER            IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LRWORK,
     $                   LWORK, M, N, NZ
      REAL             VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * ), IWORK( * )
      REAL               W( * ), RWORK( * )
      COMPLEX            A( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PCHEEVR computes selected eigenvalues and, optionally, eigenvectors
*  of a complex Hermitian matrix A distributed in 2D blockcyclic format
*  by calling the recommended sequence of ScaLAPACK routines.  
*
*  First, the matrix A is reduced to real symmetric tridiagonal form.
*  Then, the eigenproblem is solved using the parallel MRRR algorithm.
*  Last, if eigenvectors have been computed, a backtransformation is done.
*
*  Upon successful completion, each processor stores a copy of all computed
*  eigenvalues in W. The eigenvector matrix Z is stored in 
*  2D blockcyclic format distributed over all processors.
*
*  For constructive feedback and comments, please contact cvoemel@lbl.gov
*  C. Voemel
*
*
*  Arguments
*  =========
*
*  JOBZ    (global input) CHARACTER*1
*          Specifies whether or not to compute the eigenvectors:
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (global input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the interval [VL,VU] will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  UPLO    (global input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (global input) INTEGER
*          The number of rows and columns of the matrix A.  N >= 0
*
*  A       (local input/workspace) 2D block cyclic COMPLEX array,
*          global dimension (N, N),
*          local dimension ( LLD_A, LOCc(JA+N-1) )
*          (see Notes below for more detailed explanation of 2d arrays)  
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
*          It should be set to 1 when operating on a full matrix.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*          It should be set to 1 when operating on a full matrix.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          (The ScaLAPACK descriptor length is DLEN_ = 9.)
*          The array descriptor for the distributed matrix A.
*          The descriptor stores details about the 2D block-cyclic 
*          storage, see the notes below. 
*          If DESCA is incorrect, PCHEEVR cannot work correctly.
*          Also note the array alignment requirements specified below
*
*  VL      (global input) REAL 
*          If RANGE='V', the lower bound of the interval to be searched
*          for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
*
*  VU      (global input) REAL 
*          If RANGE='V', the upper bound of the interval to be searched
*          for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (global input) INTEGER
*          If RANGE='I', the index (from smallest to largest) of the
*          smallest eigenvalue to be returned.  IL >= 1.
*          Not referenced if RANGE = 'A'.
*
*  IU      (global input) INTEGER
*          If RANGE='I', the index (from smallest to largest) of the
*          largest eigenvalue to be returned.  min(IL,N) <= IU <= N.
*          Not referenced if RANGE = 'A'.
*
*  M       (global output) INTEGER
*          Total number of eigenvalues found.  0 <= M <= N.
*
*  NZ      (global output) INTEGER
*          Total number of eigenvectors computed.  0 <= NZ <= M.
*          The number of columns of Z that are filled.
*          If JOBZ .NE. 'V', NZ is not referenced.
*          If JOBZ .EQ. 'V', NZ = M 
*
*  W       (global output) REAL array, dimension (N)
*          On normal exit, the first M entries contain the selected
*          eigenvalues in ascending order.
*
*  Z       (local output) COMPLEX array,
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
*          It should be set to 1 when operating on a full matrix.
*
*  JZ      (global input) INTEGER
*          Z's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*          It should be set to 1 when operating on a full matrix.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*          DESCZ( CTXT_ ) must equal DESCA( CTXT_ )
*
*  WORK    (local workspace/output) COMPLEX  array,
*          dimension (LWORK)
*          WORK(1) returns workspace adequate workspace to allow
*          optimal performance.
*
*  LWORK  (local input) INTEGER
*          Size of WORK array, must be at least 3.
*          If only eigenvalues are requested:
*            LWORK >= N + MAX( NB * ( NP00 + 1 ), NB * 3 )
*          If eigenvectors are requested:
*            LWORK >= N + ( NP00 + MQ00 + NB ) * NB
*          For definitions of NP00 & MQ00, see LRWORK. 
*
*          For optimal performance, greater workspace is needed, i.e.
*            LWORK >= MAX( LWORK, NHETRD_LWORK )
*          Where LWORK is as defined above, and
*          NHETRD_LWORK = N + 2*( ANB+1 )*( 4*NPS+2 ) +
*            ( NPS + 1 ) * NPS
*
*          ICTXT = DESCA( CTXT_ )
*          ANB = PJLAENV( ICTXT, 3, 'PCHETTRD', 'L', 0, 0, 0, 0 )
*          SQNPC = SQRT( REAL( NPROW * NPCOL ) )
*          NPS = MAX( NUMROC( N, 1, 0, 0, SQNPC ), 2*ANB )
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the
*          optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*          NOTE THAT FOR OPTIMAL PERFORMANCE, LWOPT IS RETURNED
*          (THE OPTIMUM WORKSPACE) RATHER THAN THE MINIMUM NECESSARY
*          WORKSPACE LWMIN WHEN A WORKSPACE QUERY IS ISSUED.
*          FOR VERY SMALL MATRICES, LWOPT >> LWMIN.
*
*  RWORK    (local workspace/output) REAL  array,
*          dimension (LRWORK)
*          On return, RWORK(1) contains the optimal amount of
*          workspace required for efficient execution.
*          if JOBZ='N' RWORK(1) = optimal amount of workspace
*             required to compute the eigenvalues.
*          if JOBZ='V' RWORK(1) = optimal amount of workspace
*             required to compute eigenvalues and eigenvectors.
*
*  LRWORK  (local input) INTEGER
*          Size of RWORK, must be at least 3.
*          See below for definitions of variables used to define LRWORK.
*          If no eigenvectors are requested (JOBZ = 'N') then
*             LRWORK >= 2 + 5 * N + MAX( 12 * N, NB * ( NP00 + 1 ) )
*          If eigenvectors are requested (JOBZ = 'V' ) then
*             the amount of workspace required is:
*             LRWORK >= 2 + 5 * N + MAX( 18*N, NP00 * MQ00 + 2 * NB * NB ) +
*               (2 + ICEIL( NEIG, NPROW*NPCOL))*N
*
*          Variable definitions:
*             NEIG = number of eigenvectors requested
*             NB = DESCA( MB_ ) = DESCA( NB_ ) =
*                  DESCZ( MB_ ) = DESCZ( NB_ )
*             NN = MAX( N, NB, 2 )
*             DESCA( RSRC_ ) = DESCA( NB_ ) = DESCZ( RSRC_ ) =
*                              DESCZ( CSRC_ ) = 0
*             NP00 = NUMROC( NN, NB, 0, 0, NPROW )
*             MQ00 = NUMROC( MAX( NEIG, NB, 2 ), NB, 0, 0, NPCOL )
*             ICEIL( X, Y ) is a ScaLAPACK function returning
*             ceiling(X/Y)
*
*          If LRWORK = -1, then LRWORK is global input and a workspace
*          query is assumed; the routine only calculates the size
*          required for optimal performance for all work arrays. Each of
*          these values is returned in the first entry of the
*          corresponding work arrays, and no error message is issued by
*          PXERBLA.
*
*  IWORK   (local workspace) INTEGER array
*          On return, IWORK(1) contains the amount of integer workspace
*          required.
*
*  LIWORK  (local input) INTEGER
*          size of IWORK
*
*          Let  NNP = MAX( N, NPROW*NPCOL + 1, 4 ). Then:
*          LIWORK >= 12*NNP + 2*N when the eigenvectors are desired
*          LIWORK >= 10*NNP + 2*N when only the eigenvalues have to be computed
*          
*          If LIWORK = -1, then LIWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
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
*  PCHEEVR assumes IEEE 754 standard compliant arithmetic. 
*
*  Alignment requirements
*  ======================
*
*  The distributed submatrices A(IA:*, JA:*) and Z(IZ:IZ+M-1,JZ:JZ+N-1)
*  must satisfy the following alignment properties:
*
*  1.Identical (quadratic) dimension: 
*    DESCA(M_) = DESCZ(M_) = DESCA(N_) = DESCZ(N_)
*  2.Quadratic conformal blocking: 
*    DESCA(MB_) = DESCA(NB_) = DESCZ(MB_) = DESCZ(NB_)
*    DESCA(RSRC_) = DESCZ(RSRC_)
*  3.MOD( IA-1, MB_A ) = MOD( IZ-1, MB_Z ) = 0
*  4.IAROW = IZROW
*
*
*     .. Parameters ..
      INTEGER            CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_
      PARAMETER          ( CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8 )
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, COLBRT, DOBCST, FINISH, FIRST, INDEIG,
     $                   LOWER, LQUERY, VALEIG, VSTART, WANTZ
      INTEGER            ANB, DOL, DOU, DSTCOL, DSTROW, EIGCNT, FRSTCL,
     $                   I, IAROW, ICTXT, IIL, IINDERR, IINDWLC, IINFO,
     $                   IIU, IM, INDD, INDD2, INDE, INDE2, INDERR,
     $                   INDILU, INDRTAU, INDRW, INDRWORK, INDTAU,
     $                   INDWLC, INDWORK, IPIL, IPIU, IPROC, IZROW,
     $                   LASTCL, LENGTHI, LENGTHI2, LIWMIN, LLRWORK,
     $                   LLWORK, LRWMIN, LRWOPT, LWMIN, LWOPT, MAXCLS,
     $                   MQ00, MYCOL, MYIL, MYIU, MYPROC, MYROW, MZ, NB,
     $                   NDEPTH, NEEDIL, NEEDIU, NHETRD_LWOPT, NNP,
     $                   NP00, NPCOL, NPROCS, NPROW, NPS, NSPLIT,
     $                   OFFSET, PARITY, RLENGTHI, RLENGTHI2, RSTARTI,
     $                   SIZE1, SIZE2, SQNPC, SRCCOL, SRCROW, STARTI,
     $                   ZOFFSET

      REAL                        PIVMIN, SAFMIN, SCALE, VLL, VUU, WL,
     $                            WU
*
*     .. Local Arrays ..
      INTEGER            IDUM1( 4 ), IDUM2( 4 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, INDXG2P, NUMROC, PJLAENV
      REAL               PSLAMCH
      EXTERNAL            ICEIL, INDXG2P, LSAME, NUMROC, PJLAENV,
     $                    PSLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL            BLACS_GRIDINFO, CHK1MAT, IGEBR2D, IGEBS2D,
     $                    IGERV2D, IGESD2D, IGSUM2D, PCELGET, PCHENTRD,
     $                    PCHK1MAT, PCHK2MAT, PCLAEVSWP, PCUNMTR,
     $                    PSLARED1D, PXERBLA, SCOPY, SGEBR2D, SGEBS2D,
     $                    SGERV2D, SGESD2D, SLARRC, SLASRT2,
     $                    SSTEGR2A, SSTEGR2B, SSTEGR2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC           ABS, CMPLX, ICHAR, INT, MAX, MIN, MOD, REAL,
     $                    SQRT
*     ..
*     .. Executable Statements ..
*

      INFO = 0

***********************************************************************
*
*     Decode character arguments to find out what the code should do
*
***********************************************************************
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

***********************************************************************
*
*     GET MACHINE PARAMETERS
*
***********************************************************************
      ICTXT = DESCA( CTXT_ )
      SAFMIN = PSLAMCH( ICTXT, 'Safe minimum' )

***********************************************************************
*
*     Set up pointers into the (complex) WORK array
*     
***********************************************************************
      INDTAU = 1
      INDWORK = INDTAU + N
      LLWORK = LWORK - INDWORK + 1

***********************************************************************
*
*     Set up pointers into the RWORK array
*     
***********************************************************************
      INDRTAU = 1
      INDD = INDRTAU + N
      INDE = INDD + N + 1
      INDD2 = INDE + N + 1
      INDE2 = INDD2 + N
      INDRWORK = INDE2 + N
      LLRWORK = LRWORK - INDRWORK + 1

***********************************************************************
*
*     BLACS PROCESSOR GRID SETUP
*
***********************************************************************
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )


      NPROCS = NPROW * NPCOL
      MYPROC = MYROW * NPCOL + MYCOL
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 800+CTXT_ )
      ELSE IF( WANTZ ) THEN
         IF( ICTXT.NE.DESCZ( CTXT_ ) ) THEN
            INFO = -( 2100+CTXT_ )
         END IF
      END IF

***********************************************************************
*
*     COMPUTE REAL WORKSPACE
*
***********************************************************************
      IF ( ALLEIG ) THEN
         MZ = N
      ELSE IF ( INDEIG ) THEN
         MZ = IU - IL + 1
      ELSE
*        Take upper bound for VALEIG case
         MZ = N
      END IF
*     
      NB =  DESCA( NB_ )
      NP00 = NUMROC( N, NB, 0, 0, NPROW )
      MQ00 = NUMROC( MZ, NB, 0, 0, NPCOL )            
      IF ( WANTZ ) THEN
         INDRW = INDRWORK + MAX(18*N, NP00*MQ00 + 2*NB*NB)
         LRWMIN = INDRW - 1 + (ICEIL(MZ, NPROCS) + 2)*N
         LWMIN = N + MAX((NP00 + MQ00 + NB) * NB, 3 * NB)
      ELSE
         INDRW = INDRWORK + 12*N
         LRWMIN = INDRW - 1
         LWMIN = N + MAX( NB*( NP00 + 1 ), 3 * NB ) 
      END IF

*     The code that validates the input requires 3 workspace entries
      LRWMIN = MAX(3, LRWMIN)
      LRWOPT = LRWMIN
      LWMIN = MAX(3, LWMIN)
      LWOPT = LWMIN
*
      ANB = PJLAENV( ICTXT, 3, 'PCHETTRD', 'L', 0, 0, 0, 0 )
      SQNPC = INT( SQRT( REAL( NPROCS ) ) )
      NPS = MAX( NUMROC( N, 1, 0, 0, SQNPC ), 2*ANB )
      NHETRD_LWOPT = 2*( ANB+1 )*( 4*NPS+2 ) + ( NPS+4 )*NPS
      LWOPT = MAX( LWOPT, N+NHETRD_LWOPT )
*
      SIZE1 = INDRW - INDRWORK

***********************************************************************
*
*     COMPUTE INTEGER WORKSPACE
*
***********************************************************************
      NNP = MAX( N, NPROCS+1, 4 )
      IF ( WANTZ ) THEN
        LIWMIN = 12*NNP + 2*N 
      ELSE
        LIWMIN = 10*NNP + 2*N
      END IF

***********************************************************************
*
*     Set up pointers into the IWORK array
*     
***********************************************************************
*     Pointer to eigenpair distribution over processors
      INDILU = LIWMIN - 2*NPROCS + 1            
      SIZE2 = INDILU - 2*N 
	

***********************************************************************
*
*     Test the input arguments.
*
***********************************************************************
      IF( INFO.EQ.0 ) THEN
         CALL CHK1MAT( N, 4, N, 4, IA, JA, DESCA, 8, INFO )
         IF( WANTZ )
     $      CALL CHK1MAT( N, 4, N, 4, IZ, JZ, DESCZ, 21, INFO )
*
         IF( INFO.EQ.0 ) THEN
            IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
               INFO = -1
            ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
               INFO = -2
            ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
               INFO = -3
            ELSE IF( MOD( IA-1, DESCA( MB_ ) ).NE.0 ) THEN
               INFO = -6
            ELSE IF( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) THEN
               INFO = -10
            ELSE IF( INDEIG .AND. ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) )
     $                THEN
               INFO = -11
            ELSE IF( INDEIG .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ))
     $                THEN
               INFO = -12
            ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -21
            ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -23
            ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
               INFO = -25
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -( 800+NB_ )
            END IF
            IF( WANTZ ) THEN
               IAROW = INDXG2P( 1, DESCA( NB_ ), MYROW, 
     $                       DESCA( RSRC_ ), NPROW )
               IZROW = INDXG2P( 1, DESCA( NB_ ), MYROW, 
     $                          DESCZ( RSRC_ ), NPROW )
               IF( IAROW.NE.IZROW ) THEN
                  INFO = -19
               ELSE IF( MOD( IA-1, DESCA( MB_ ) ).NE.
     $             MOD( IZ-1, DESCZ( MB_ ) ) ) THEN
                  INFO = -19
               ELSE IF( DESCA( M_ ).NE.DESCZ( M_ ) ) THEN
                  INFO = -( 2100+M_ )
               ELSE IF( DESCA( N_ ).NE.DESCZ( N_ ) ) THEN
                  INFO = -( 2100+N_ )
               ELSE IF( DESCA( MB_ ).NE.DESCZ( MB_ ) ) THEN
                  INFO = -( 2100+MB_ )
               ELSE IF( DESCA( NB_ ).NE.DESCZ( NB_ ) ) THEN
                  INFO = -( 2100+NB_ )
               ELSE IF( DESCA( RSRC_ ).NE.DESCZ( RSRC_ ) ) THEN
                  INFO = -( 2100+RSRC_ )
               ELSE IF( DESCA( CSRC_ ).NE.DESCZ( CSRC_ ) ) THEN
                  INFO = -( 2100+CSRC_ )
               ELSE IF( ICTXT.NE.DESCZ( CTXT_ ) ) THEN
                  INFO = -( 2100+CTXT_ )
               END IF
            END IF
         END IF
         IDUM2( 1 ) = 1
         IF( LOWER ) THEN
            IDUM1( 2 ) = ICHAR( 'L' )
         ELSE
            IDUM1( 2 ) = ICHAR( 'U' )
         END IF
         IDUM2( 2 ) = 2
         IF( ALLEIG ) THEN
            IDUM1( 3 ) = ICHAR( 'A' )
         ELSE IF( INDEIG ) THEN
            IDUM1( 3 ) = ICHAR( 'I' )
         ELSE
            IDUM1( 3 ) = ICHAR( 'V' )
         END IF
         IDUM2( 3 ) = 3
         IF( LQUERY ) THEN
            IDUM1( 4 ) = -1
         ELSE
            IDUM1( 4 ) = 1
         END IF
         IDUM2( 4 ) = 4
         IF( WANTZ ) THEN
            IDUM1( 1 ) = ICHAR( 'V' )
            CALL PCHK2MAT( N, 4, N, 4, IA, JA, DESCA, 8, N, 4, N, 4,IZ,
     $                     JZ, DESCZ, 21, 4, IDUM1, IDUM2, INFO )
         ELSE
            IDUM1( 1 ) = ICHAR( 'N' )
            CALL PCHK1MAT( N, 4, N, 4, IA, JA, DESCA, 8, 4, IDUM1,
     $                     IDUM2, INFO )
         END IF
         WORK( 1 ) = CMPLX( LWOPT )
         RWORK( 1 ) = REAL( LRWOPT )
         IWORK( 1 ) = LIWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PCHEEVR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

***********************************************************************
*
*     Quick return if possible
*
***********************************************************************
      IF( N.EQ.0 ) THEN
         IF( WANTZ ) THEN
            NZ = 0
         END IF
         M = 0
         WORK( 1 ) = CMPLX( LWOPT )
         RWORK( 1 ) = REAL( LRWOPT )
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF

      IF( VALEIG ) THEN
         VLL = VL
         VUU = VU
      ELSE
         VLL = ZERO
         VUU = ZERO
      END IF
*
*     No scaling done here, leave this to MRRR kernel.
*     Scale tridiagonal rather than full matrix.
*
***********************************************************************
*
*     REDUCE MATRIX TO REAL SYMMETRIC TRIDIAGONAL FORM.
*
***********************************************************************


      CALL PCHENTRD( UPLO, N, A, IA, JA, DESCA, RWORK( INDD ),
     $               RWORK( INDE ), WORK( INDTAU ), WORK( INDWORK ),
     $               LLWORK, RWORK( INDRWORK ), LLRWORK,IINFO )


      IF (IINFO .NE. 0) THEN
         CALL PXERBLA( ICTXT, 'PCHENTRD', -IINFO )
         RETURN
      END IF

***********************************************************************
*
*     DISTRIBUTE TRIDIAGONAL TO ALL PROCESSORS
*
***********************************************************************
      OFFSET = 0
      IF( IA.EQ.1 .AND. JA.EQ.1 .AND. 
     $    DESCA( RSRC_ ).EQ.0 .AND. DESCA( CSRC_ ).EQ.0 )
     $   THEN
         CALL PSLARED1D( N, IA, JA, DESCA, RWORK( INDD ), 
     $                   RWORK( INDD2 ), RWORK( INDRWORK ), LLRWORK )
*
         CALL PSLARED1D( N, IA, JA, DESCA, RWORK( INDE ), 
     $                   RWORK( INDE2 ), RWORK( INDRWORK ), LLRWORK )
         IF( .NOT.LOWER )
     $      OFFSET = 1
      ELSE
         DO 10 I = 1, N
            CALL PCELGET( 'A', ' ', WORK( INDWORK ), A, 
     $                    I+IA-1, I+JA-1, DESCA )
            RWORK( INDD2+I-1 ) = REAL( WORK( INDWORK ) )
   10    CONTINUE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 I = 1, N - 1
               CALL PCELGET( 'A', ' ', WORK( INDWORK ), A, 
     $                       I+IA-1, I+JA, DESCA )
               RWORK( INDE2+I-1 ) = REAL( WORK( INDWORK ) )
   20       CONTINUE
         ELSE
            DO 30 I = 1, N - 1
               CALL PCELGET( 'A', ' ', WORK( INDWORK ), A,
     $                       I+IA, I+JA-1, DESCA )
               RWORK( INDE2+I-1 ) = REAL( WORK( INDWORK ) )
   30       CONTINUE
         END IF
      END IF




***********************************************************************
*
*     SET IIL, IIU
*
***********************************************************************
      IF ( ALLEIG ) THEN 
         IIL = 1
         IIU = N
      ELSE IF ( INDEIG ) THEN
         IIL = IL
         IIU = IU
      ELSE IF ( VALEIG ) THEN
         CALL SLARRC('T', N, VLL, VUU, RWORK( INDD2 ), 
     $    RWORK( INDE2 + OFFSET ), SAFMIN, EIGCNT, IIL, IIU, INFO)
*        Refine upper bound N that was taken 
         MZ = EIGCNT
         IIL = IIL + 1
      ENDIF

      IF(MZ.EQ.0) THEN
         M = 0
         IF( WANTZ ) THEN
            NZ = 0
         END IF
         WORK( 1 ) = REAL( LWOPT )
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF

      MYIL = 0
      MYIU = 0
      M = 0
      IM = 0

***********************************************************************
*
*     COMPUTE WORK ASSIGNMENTS
*
***********************************************************************

*
*     Each processor computes the work assignments for all processors
*
      CALL PMPIM2( IIL, IIU, NPROCS,
     $             IWORK(INDILU), IWORK(INDILU+NPROCS) )
*
*     Find local work assignment
*
      MYIL = IWORK(INDILU+MYPROC)
      MYIU = IWORK(INDILU+NPROCS+MYPROC)


      ZOFFSET = MAX(0, MYIL - IIL - 1)
      FIRST = ( MYIL .EQ. IIL )


***********************************************************************
*
*     CALLS TO MRRR KERNEL
*
***********************************************************************
      IF(.NOT.WANTZ) THEN
*
*        Compute eigenvalues only.
*
         IINFO = 0
         IF ( MYIL.GT.0 ) THEN
            DOL = 1
            DOU = MYIU - MYIL + 1
            CALL SSTEGR2( JOBZ, 'I', N,  RWORK( INDD2 ),
     $                  RWORK( INDE2+OFFSET ), VLL, VUU, MYIL, MYIU,
     $                  IM, W( 1 ), RWORK( INDRW ), N, 
     $                  MYIU - MYIL + 1,
     $                  IWORK( 1 ), RWORK( INDRWORK ), SIZE1, 
     $                  IWORK( 2*N+1 ), SIZE2, 
     $                  DOL, DOU, ZOFFSET, IINFO )
*           SSTEGR2 zeroes out the entire W array, so we can't just give
*           it the part of W we need.  So here we copy the W entries into
*           their correct location
            DO 49 I = 1, IM
              W( MYIL-IIL+I ) = W( I )
 49         CONTINUE
*           W( MYIL ) is at W( MYIL - IIL + 1 )
*           W( X ) is at W(X - IIL + 1 )
         END IF
         IF (IINFO .NE. 0) THEN
            CALL PXERBLA( ICTXT, 'SSTEGR2', -IINFO )
            RETURN
         END IF
      ELSEIF ( WANTZ .AND. NPROCS.EQ.1 ) THEN
*
*        Compute eigenvalues and -vectors, but only on one processor
*
         IINFO = 0
         IF ( MYIL.GT.0 ) THEN
            DOL = MYIL - IIL + 1
            DOU = MYIU - IIL + 1
            CALL SSTEGR2( JOBZ, 'I', N,  RWORK( INDD2 ),
     $                  RWORK( INDE2+OFFSET ), VLL, VUU, IIL, IIU,
     $                  IM, W( 1 ), RWORK( INDRW ), N, 
     $                  N,
     $                  IWORK( 1 ), RWORK( INDRWORK ), SIZE1, 
     $                  IWORK( 2*N+1 ), SIZE2, DOL, DOU,
     $                  ZOFFSET, IINFO )
         ENDIF
         IF (IINFO .NE. 0) THEN
            CALL PXERBLA( ICTXT, 'SSTEGR2', -IINFO )
            RETURN
         END IF
      ELSEIF ( WANTZ ) THEN
*        Compute representations in parallel.
*        Share eigenvalue computation for root between all processors
*        Then compute the eigenvectors. 
         IINFO = 0
*        Part 1. compute root representations and root eigenvalues
         IF ( MYIL.GT.0 ) THEN
            DOL = MYIL - IIL + 1
            DOU = MYIU - IIL + 1
            CALL SSTEGR2A( JOBZ, 'I', N,  RWORK( INDD2 ),
     $                  RWORK( INDE2+OFFSET ), VLL, VUU, IIL, IIU,
     $                  IM, W( 1 ), RWORK( INDRW ), N, 
     $                  N, RWORK( INDRWORK ), SIZE1, 
     $                  IWORK( 2*N+1 ), SIZE2, DOL, 
     $                  DOU, NEEDIL, NEEDIU,
     $                  INDERR, NSPLIT, PIVMIN, SCALE, WL, WU,
     $                  IINFO )
         ENDIF
         IF (IINFO .NE. 0) THEN
            CALL PXERBLA( ICTXT, 'SSTEGR2A', -IINFO )
            RETURN
         END IF
*
*        The second part of parallel MRRR, the representation tree
*        construction begins. Upon successful completion, the 
*        eigenvectors have been computed. This is indicated by
*        the flag FINISH.
*
         VSTART = .TRUE.
         FINISH = (MYIL.LE.0)
C        Part 2. Share eigenvalues and uncertainties between all processors
         IINDERR = INDRWORK + INDERR - 1

*


*
*        There are currently two ways to communicate eigenvalue information
*        using the BLACS.
*        1.) BROADCAST
*        2.) POINT2POINT between collaborators (those processors working
*            jointly on a cluster.
*        For efficiency, BROADCAST has been disabled.
*        At a later stage, other more efficient communication algorithms 
*        might be implemented, e. g. group or tree-based communication.

         DOBCST = .FALSE.
         IF(DOBCST) THEN
*           First gather everything on the first processor.
*           Then use BROADCAST-based communication 
            DO 45 I = 2, NPROCS
               IF (MYPROC .EQ. (I - 1)) THEN
                  DSTROW = 0
                  DSTCOL = 0
                  STARTI = DOL
                  IWORK(1) = STARTI
                  IF(MYIL.GT.0) THEN
                     LENGTHI = MYIU - MYIL + 1
                  ELSE
                     LENGTHI = 0
                  ENDIF
                  IWORK(2) = LENGTHI
                  CALL IGESD2D( ICTXT, 2, 1, IWORK, 2, 
     $                    DSTROW, DSTCOL )
                  IF (( STARTI.GE.1 ) .AND. ( LENGTHI.GE.1 )) THEN
                     LENGTHI2 = 2*LENGTHI
*                    Copy eigenvalues into communication buffer
                     CALL SCOPY(LENGTHI,W( STARTI ),1,
     $                          RWORK( INDD ), 1)                    
*                    Copy uncertainties into communication buffer
                     CALL SCOPY(LENGTHI,RWORK(IINDERR+STARTI-1),1,
     $                          RWORK( INDD+LENGTHI ), 1)                    
*                    send buffer
                     CALL SGESD2D( ICTXT, LENGTHI2, 
     $                    1, RWORK( INDD ), LENGTHI2,
     $                    DSTROW, DSTCOL )
                  END IF
               ELSE IF (MYPROC .EQ. 0) THEN
                  SRCROW = (I-1) / NPCOL
                  SRCCOL = MOD(I-1, NPCOL)
                  CALL IGERV2D( ICTXT, 2, 1, IWORK, 2, 
     $                    SRCROW, SRCCOL )
                  STARTI = IWORK(1)
                  LENGTHI = IWORK(2)
                  IF (( STARTI.GE.1 ) .AND. ( LENGTHI.GE.1 )) THEN
                     LENGTHI2 = 2*LENGTHI
*                    receive buffer
                     CALL SGERV2D( ICTXT, LENGTHI2, 1,
     $                 RWORK(INDD), LENGTHI2, SRCROW, SRCCOL )
*                    copy eigenvalues from communication buffer
                     CALL SCOPY( LENGTHI, RWORK(INDD), 1,
     $                          W( STARTI ), 1)                    
*                    copy uncertainties (errors) from communication buffer
                     CALL SCOPY(LENGTHI,RWORK(INDD+LENGTHI),1,
     $                          RWORK( IINDERR+STARTI-1 ), 1)     
                  END IF
               END IF
  45        CONTINUE
            LENGTHI = IIU - IIL + 1
            LENGTHI2 = LENGTHI * 2
            IF (MYPROC .EQ. 0) THEN
*              Broadcast eigenvalues and errors to all processors
               CALL SCOPY(LENGTHI,W ,1, RWORK( INDD ), 1)                 
               CALL SCOPY(LENGTHI,RWORK( IINDERR ),1,
     $                          RWORK( INDD+LENGTHI ), 1)                    
               CALL SGEBS2D( ICTXT, 'A', ' ', LENGTHI2, 1, 
     $              RWORK(INDD), LENGTHI2 )
            ELSE
               SRCROW = 0
               SRCCOL = 0
               CALL SGEBR2D( ICTXT, 'A', ' ', LENGTHI2, 1,
     $             RWORK(INDD), LENGTHI2, SRCROW, SRCCOL )
               CALL SCOPY( LENGTHI, RWORK(INDD), 1, W, 1)
               CALL SCOPY(LENGTHI,RWORK(INDD+LENGTHI),1,
     $                          RWORK( IINDERR ), 1)                   
            END IF
         ELSE
*           Enable point2point communication between collaborators

*           Find collaborators of MYPROC            
            IF( (NPROCS.GT.1).AND.(MYIL.GT.0) ) THEN
               CALL PMPCOL( MYPROC, NPROCS, IIL, NEEDIL, NEEDIU, 
     $                   IWORK(INDILU), IWORK(INDILU+NPROCS),
     $                   COLBRT, FRSTCL, LASTCL )
            ELSE
               COLBRT = .FALSE.
            ENDIF

            IF(COLBRT) THEN
*              If the processor collaborates with others,
*              communicate information. 
               DO 47 IPROC = FRSTCL, LASTCL
                  IF (MYPROC .EQ. IPROC) THEN
                     STARTI = DOL
                     IWORK(1) = STARTI
                     LENGTHI = MYIU - MYIL + 1
                     IWORK(2) = LENGTHI
                     
                     IF ((STARTI.GE.1) .AND. (LENGTHI.GE.1)) THEN
*                       Copy eigenvalues into communication buffer
                        CALL SCOPY(LENGTHI,W( STARTI ),1,
     $                              RWORK(INDD), 1)                    
*                       Copy uncertainties into communication buffer
                        CALL SCOPY(LENGTHI,
     $                          RWORK( IINDERR+STARTI-1 ),1,
     $                          RWORK(INDD+LENGTHI), 1)                    
                     ENDIF

                     DO 46 I = FRSTCL, LASTCL                      
                        IF(I.EQ.MYPROC) GOTO 46
                        DSTROW = I/ NPCOL
                        DSTCOL = MOD(I, NPCOL)
                        CALL IGESD2D( ICTXT, 2, 1, IWORK, 2, 
     $                             DSTROW, DSTCOL )
                        IF ((STARTI.GE.1) .AND. (LENGTHI.GE.1)) THEN
                           LENGTHI2 = 2*LENGTHI
*                          send buffer
                           CALL SGESD2D( ICTXT, LENGTHI2, 
     $                          1, RWORK(INDD), LENGTHI2,
     $                          DSTROW, DSTCOL )
                        END IF
  46                 CONTINUE
                  ELSE
                     SRCROW = IPROC / NPCOL
                     SRCCOL = MOD(IPROC, NPCOL)
                     CALL IGERV2D( ICTXT, 2, 1, IWORK, 2, 
     $                             SRCROW, SRCCOL )
                     RSTARTI = IWORK(1)
                     RLENGTHI = IWORK(2)
                     IF ((RSTARTI.GE.1 ) .AND. (RLENGTHI.GE.1 )) THEN
                        RLENGTHI2 = 2*RLENGTHI
                        CALL SGERV2D( ICTXT, RLENGTHI2, 1,
     $                      RWORK(INDE), RLENGTHI2,
     $                      SRCROW, SRCCOL )
*                       copy eigenvalues from communication buffer
                        CALL SCOPY( RLENGTHI,RWORK(INDE), 1,
     $                          W( RSTARTI ), 1)                    
*                       copy uncertainties (errors) from communication buffer
                        CALL SCOPY(RLENGTHI,RWORK(INDE+RLENGTHI),1,
     $                          RWORK( IINDERR+RSTARTI-1 ), 1)                    
                     END IF
                  END IF
  47           CONTINUE
            ENDIF
         ENDIF

*        Part 3. Compute representation tree and eigenvectors.
*                What follows is a loop in which the tree
*                is constructed in parallel from top to bottom,
*                on level at a time, until all eigenvectors
*                have been computed.
*      
 100     CONTINUE
         IF ( MYIL.GT.0 ) THEN
            CALL SSTEGR2B( JOBZ, N,  RWORK( INDD2 ),
     $                  RWORK( INDE2+OFFSET ), 
     $                  IM, W( 1 ), RWORK( INDRW ), N, N,
     $                  IWORK( 1 ), RWORK( INDRWORK ), SIZE1, 
     $                  IWORK( 2*N+1 ), SIZE2, DOL, 
     $                  DOU, NEEDIL, NEEDIU, INDWLC,
     $                  PIVMIN, SCALE, WL, WU,
     $                  VSTART, FINISH, 
     $                  MAXCLS, NDEPTH, PARITY, ZOFFSET, IINFO )
            IINDWLC = INDRWORK + INDWLC - 1
            IF(.NOT.FINISH) THEN
               IF((NEEDIL.LT.DOL).OR.(NEEDIU.GT.DOU)) THEN
                  CALL PMPCOL( MYPROC, NPROCS, IIL, NEEDIL, NEEDIU,
     $                 IWORK(INDILU), IWORK(INDILU+NPROCS),
     $                   COLBRT, FRSTCL, LASTCL )
               ELSE
                  COLBRT = .FALSE.
                  FRSTCL = MYPROC
                  LASTCL = MYPROC
               ENDIF
*
*              Check if this processor collaborates, i.e. 
*              communication is needed.
*
               IF(COLBRT) THEN
                  DO 147 IPROC = FRSTCL, LASTCL
                     IF (MYPROC .EQ. IPROC) THEN
                        STARTI = DOL
                        IWORK(1) = STARTI
                        IF(MYIL.GT.0) THEN
                           LENGTHI = MYIU - MYIL + 1
                        ELSE
                           LENGTHI = 0
                        ENDIF
                        IWORK(2) = LENGTHI
                        IF ((STARTI.GE.1).AND.(LENGTHI.GE.1)) THEN
*                          Copy eigenvalues into communication buffer
                           CALL SCOPY(LENGTHI,
     $                          RWORK( IINDWLC+STARTI-1 ),1,
     $                          RWORK(INDD), 1)                    
*                          Copy uncertainties into communication buffer
                           CALL SCOPY(LENGTHI,
     $                          RWORK( IINDERR+STARTI-1 ),1,
     $                          RWORK(INDD+LENGTHI), 1)                    
                        ENDIF
                      
                        DO 146 I = FRSTCL, LASTCL                      
                           IF(I.EQ.MYPROC) GOTO 146
                           DSTROW = I/ NPCOL
                           DSTCOL = MOD(I, NPCOL)
                           CALL IGESD2D( ICTXT, 2, 1, IWORK, 2, 
     $                             DSTROW, DSTCOL )
                           IF ((STARTI.GE.1).AND.(LENGTHI.GE.1)) THEN
                              LENGTHI2 = 2*LENGTHI
*                             send buffer
                              CALL SGESD2D( ICTXT, LENGTHI2, 
     $                             1, RWORK(INDD), LENGTHI2,
     $                             DSTROW, DSTCOL )
                           END IF
 146                    CONTINUE
                     ELSE
                        SRCROW = IPROC / NPCOL
                        SRCCOL = MOD(IPROC, NPCOL)
                        CALL IGERV2D( ICTXT, 2, 1, IWORK, 2, 
     $                             SRCROW, SRCCOL )
                        RSTARTI = IWORK(1)
                        RLENGTHI = IWORK(2)
                        IF ((RSTARTI.GE.1).AND.(RLENGTHI.GE.1)) THEN
                           RLENGTHI2 = 2*RLENGTHI
                           CALL SGERV2D( ICTXT,RLENGTHI2, 1,
     $                         RWORK(INDE),RLENGTHI2,
     $                         SRCROW, SRCCOL )
*                          copy eigenvalues from communication buffer
                           CALL SCOPY(RLENGTHI,RWORK(INDE), 1,
     $                          RWORK( IINDWLC+RSTARTI-1 ), 1)        
*                          copy uncertainties (errors) from communication buffer
                           CALL SCOPY(RLENGTHI,RWORK(INDE+RLENGTHI),
     $                          1,RWORK( IINDERR+RSTARTI-1 ), 1)            
                        END IF
                      END IF
 147              CONTINUE
               ENDIF
               GOTO 100         
            ENDIF
         ENDIF
         IF (IINFO .NE. 0) THEN
            CALL PXERBLA( ICTXT, 'SSTEGR2B', -IINFO )
            RETURN
         END IF
*
      ENDIF

*
***********************************************************************
*
*     MAIN PART ENDS HERE
*
***********************************************************************
*

***********************************************************************
*
*     ALLGATHER: EACH PROCESSOR SENDS ITS EIGENVALUES TO THE FIRST ONE,
*                THEN THE FIRST PROCESSOR BROADCASTS ALL EIGENVALUES
*
***********************************************************************

      DO 50 I = 2, NPROCS
         IF (MYPROC .EQ. (I - 1)) THEN
            DSTROW = 0
            DSTCOL = 0
            STARTI = MYIL - IIL + 1
            IWORK(1) = STARTI
            IF(MYIL.GT.0) THEN
               LENGTHI = MYIU - MYIL + 1
            ELSE
               LENGTHI = 0
            ENDIF
            IWORK(2) = LENGTHI
            CALL IGESD2D( ICTXT, 2, 1, IWORK, 2, 
     $                    DSTROW, DSTCOL )
            IF ((STARTI.GE.1).AND.(LENGTHI.GE.1)) THEN
               CALL SGESD2D( ICTXT, LENGTHI, 
     $              1, W( STARTI ), LENGTHI,
     $              DSTROW, DSTCOL )
            ENDIF
         ELSE IF (MYPROC .EQ. 0) THEN
            SRCROW = (I-1) / NPCOL
            SRCCOL = MOD(I-1, NPCOL)
            CALL IGERV2D( ICTXT, 2, 1, IWORK, 2, 
     $                    SRCROW, SRCCOL )
            STARTI = IWORK(1)
            LENGTHI = IWORK(2)
            IF ((STARTI.GE.1).AND.(LENGTHI.GE.1)) THEN
               CALL SGERV2D( ICTXT, LENGTHI, 1,
     $                 W( STARTI ), LENGTHI, SRCROW, SRCCOL )
            ENDIF
         ENDIF
   50 CONTINUE

*     Accumulate M from all processors
      M = IM
      CALL IGSUM2D( ICTXT, 'A', ' ', 1, 1, M, 1, -1, -1 )

*     Broadcast eigenvalues to all processors
      IF (MYPROC .EQ. 0) THEN
*        Send eigenvalues
         CALL SGEBS2D( ICTXT, 'A', ' ', M, 1, W, M )
      ELSE
         SRCROW = 0
         SRCCOL = 0
         CALL SGEBR2D( ICTXT, 'A', ' ', M, 1,
     $           W, M, SRCROW, SRCCOL )
      END IF
*
*     Sort the eigenvalues and keep permutation in IWORK to
*     sort the eigenvectors accordingly
*
      DO 160 I = 1, M
         IWORK( NPROCS+1+I ) = I
  160 CONTINUE
      CALL SLASRT2( 'I', M, W, IWORK( NPROCS+2 ), IINFO )
      IF (IINFO.NE.0) THEN
         CALL PXERBLA( ICTXT, 'SLASRT2', -IINFO )
         RETURN
      END IF

***********************************************************************
*
*     TRANSFORM Z FROM 1D WORKSPACE INTO 2D BLOCKCYCLIC STORAGE     
*
***********************************************************************
      IF ( WANTZ ) THEN
         DO 170 I = 1, M
            IWORK( M+NPROCS+1+IWORK( NPROCS+1+I ) ) = I
  170    CONTINUE
*        Store NVS in IWORK(1:NPROCS+1) for PCLAEVSWP
         IWORK( 1 ) = 0
         DO 180 I = 1, NPROCS
*           Find IL and IU for processor i-1
*           Has already been computed by PMPIM2 and stored
            IPIL = IWORK(INDILU+I-1)
            IPIU = IWORK(INDILU+NPROCS+I-1)
            IF (IPIL .EQ. 0) THEN
               IWORK( I + 1 ) = IWORK( I )
            ELSE
               IWORK( I + 1 ) = IWORK( I ) + IPIU - IPIL + 1
            ENDIF
  180    CONTINUE

         IF ( FIRST ) THEN
            CALL PCLAEVSWP(N, RWORK( INDRW ), N, Z, IZ, JZ, 
     $       DESCZ, IWORK( 1 ), IWORK( NPROCS+M+2 ), RWORK( INDRWORK ), 
     $       SIZE1 )
         ELSE
            CALL PCLAEVSWP(N, RWORK( INDRW + N ), N, Z, IZ, JZ, 
     $       DESCZ, IWORK( 1 ), IWORK( NPROCS+M+2 ), RWORK( INDRWORK ),
     $       SIZE1 )
         END IF
*
         NZ = M
*

***********************************************************************
*
*       Compute eigenvectors of A from eigenvectors of T
*
***********************************************************************
         IF( NZ.GT.0 ) THEN
           CALL PCUNMTR( 'L', UPLO, 'N', N, NZ, A, IA, JA, DESCA,
     $                    WORK( INDTAU ), Z, IZ, JZ, DESCZ,
     $                    WORK( INDWORK ), LLWORK, IINFO )
         END IF
         IF (IINFO.NE.0) THEN
            CALL PXERBLA( ICTXT, 'PCUNMTR', -IINFO )
            RETURN
         END IF
*

      END IF
*
      WORK( 1 ) = CMPLX( LWOPT )
      RWORK( 1 ) = REAL( LRWOPT )
      IWORK( 1 ) = LIWMIN

      RETURN
*
*     End of PCHEEVR
*
      END
