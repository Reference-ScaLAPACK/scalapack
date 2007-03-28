      SUBROUTINE PZSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, ORFAC, Z, IZ,
     $                    JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, IFAIL,
     $                    ICLUSTR, GAP, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INFO, IZ, JZ, LIWORK, LWORK, M, N
      DOUBLE PRECISION   ORFAC
*     ..
*     .. Array Arguments ..
      INTEGER            DESCZ( * ), IBLOCK( * ), ICLUSTR( * ),
     $                   IFAIL( * ), ISPLIT( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), GAP( * ), W( * ), WORK( * )
      COMPLEX*16         Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PZSTEIN computes the eigenvectors of a symmetric tridiagonal matrix
*  in parallel, using inverse iteration. The eigenvectors found
*  correspond to user specified eigenvalues. PZSTEIN does not
*  orthogonalize vectors that are on different processes. The extent
*  of orthogonalization is controlled by the input parameter LWORK.
*  Eigenvectors that are to be orthogonalized are computed by the same
*  process. PZSTEIN decides on the allocation of work among the
*  processes and then calls DSTEIN2 (modified LAPACK routine) on each
*  individual process. If insufficient workspace is allocated, the
*  expected orthogonalization may not be done.
*
*  Note : If the eigenvectors obtained are not orthogonal, increase
*         LWORK and run the code again.
*
*  Notes
*  =====
*
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
*  and assume that its process grid has dimension r x c.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the r processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the c processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*
*  Arguments
*  =========
*
*          P = NPROW * NPCOL is the total number of processes
*
*  N       (global input) INTEGER
*          The order of the tridiagonal matrix T.  N >= 0.
*
*  D       (global input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the tridiagonal matrix T.
*
*  E       (global input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) off-diagonal elements of the tridiagonal matrix T.
*
*  M       (global input) INTEGER
*          The total number of eigenvectors to be found. 0 <= M <= N.
*
*  W       (global input/global output) DOUBLE PRECISION array, dim (M)
*          On input, the first M elements of W contain all the
*          eigenvalues for which eigenvectors are to be computed. The
*          eigenvalues should be grouped by split-off block and ordered
*          from smallest to largest within the block (The output array
*          W from PDSTEBZ with ORDER='b' is expected here). This
*          array should be replicated on all processes.
*          On output, the first M elements contain the input
*          eigenvalues in ascending order.
*
*          Note : To obtain orthogonal vectors, it is best if
*          eigenvalues are computed to highest accuracy ( this can be
*          done by setting ABSTOL to the underflow threshold =
*          DLAMCH('U') --- ABSTOL is an input parameter
*          to PDSTEBZ )
*
*  IBLOCK  (global input) INTEGER array, dimension (N)
*          The submatrix indices associated with the corresponding
*          eigenvalues in W -- 1 for eigenvalues belonging to the
*          first submatrix from the top, 2 for those belonging to
*          the second submatrix, etc. (The output array IBLOCK
*          from PDSTEBZ is expected here).
*
*  ISPLIT  (global input) INTEGER array, dimension (N)
*          The splitting points, at which T breaks up into submatrices.
*          The first submatrix consists of rows/columns 1 to ISPLIT(1),
*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
*          etc., and the NSPLIT-th consists of rows/columns
*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N (The output array
*          ISPLIT from PDSTEBZ is expected here.)
*
*  ORFAC   (global input) DOUBLE PRECISION
*          ORFAC specifies which eigenvectors should be orthogonalized.
*          Eigenvectors that correspond to eigenvalues which are within
*          ORFAC*||T|| of each other are to be orthogonalized.
*          However, if the workspace is insufficient (see LWORK), this
*          tolerance may be decreased until all eigenvectors to be
*          orthogonalized can be stored in one process.
*          No orthogonalization will be done if ORFAC equals zero.
*          A default value of 10^-3 is used if ORFAC is negative.
*          ORFAC should be identical on all processes.
*
*  Z       (local output) COMPLEX*16 array,
*          dimension (DESCZ(DLEN_), N/npcol + NB)
*          Z contains the computed eigenvectors associated with the
*          specified eigenvalues. Any vector which fails to converge is
*          set to its current iterate after MAXITS iterations ( See
*          DSTEIN2 ).
*          On output, Z is distributed across the P processes in block
*          cyclic format.
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
*
*  WORK    (local workspace/global output) DOUBLE PRECISION array,
*          dimension ( LWORK )
*          On output, WORK(1) gives a lower bound on the
*          workspace ( LWORK ) that guarantees the user desired
*          orthogonalization (see ORFAC).
*          Note that this may overestimate the minimum workspace needed.
*
*  LWORK   (local input) integer
*          LWORK  controls the extent of orthogonalization which can be
*          done. The number of eigenvectors for which storage is
*          allocated on each process is
*                NVEC = floor(( LWORK- max(5*N,NP00*MQ00) )/N).
*          Eigenvectors corresponding to eigenvalue clusters of size
*          NVEC - ceil(M/P) + 1 are guaranteed to be orthogonal ( the
*          orthogonality is similar to that obtained from ZSTEIN2).
*          Note : LWORK  must be no smaller than:
*               max(5*N,NP00*MQ00) + ceil(M/P)*N,
*          and should have the same input value on all processes.
*          It is the minimum value of LWORK input on different processes
*          that is significant.
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  IWORK   (local workspace/global output) INTEGER array,
*          dimension ( 3*N+P+1 )
*          On return, IWORK(1) contains the amount of integer workspace
*          required.
*          On return, the IWORK(2) through IWORK(P+2) indicate
*          the eigenvectors computed by each process. Process I computes
*          eigenvectors indexed IWORK(I+2)+1 thru' IWORK(I+3).
*
*  LIWORK  (local input) INTEGER
*          Size of array IWORK.  Must be >= 3*N + P + 1
*
*          If LIWORK = -1, then LIWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          and optimal size for all work arrays. Each of these
*          values is returned in the first entry of the corresponding
*          work array, and no error message is issued by PXERBLA.
*
*  IFAIL   (global output) integer array, dimension (M)
*          On normal exit, all elements of IFAIL are zero.
*          If one or more eigenvectors fail to converge after MAXITS
*          iterations (as in ZSTEIN), then INFO > 0 is returned.
*          If mod(INFO,M+1)>0, then
*          for I=1 to mod(INFO,M+1), the eigenvector
*          corresponding to the eigenvalue W(IFAIL(I)) failed to
*          converge ( W refers to the array of eigenvalues on output ).
*
*  ICLUSTR (global output) integer array, dimension (2*P)
*          This output array contains indices of eigenvectors
*          corresponding to a cluster of eigenvalues that could not be
*          orthogonalized due to insufficient workspace (see LWORK,
*          ORFAC and INFO). Eigenvectors corresponding to clusters of
*          eigenvalues indexed ICLUSTR(2*I-1) to ICLUSTR(2*I), I = 1 to
*          INFO/(M+1), could not be orthogonalized due to lack of
*          workspace. Hence the eigenvectors corresponding to these
*          clusters may not be orthogonal. ICLUSTR is a zero terminated
*          array --- ( ICLUSTR(2*K).NE.0 .AND. ICLUSTR(2*K+1).EQ.0 )
*          if and only if K is the number of clusters.
*
*  GAP     (global output) DOUBLE PRECISION array, dimension (P)
*          This output array contains the gap between eigenvalues whose
*          eigenvectors could not be orthogonalized. The INFO/M output
*          values in this array correspond to the INFO/(M+1) clusters
*          indicated by the array ICLUSTR. As a result, the dot product
*          between eigenvectors corresponding to the I^th cluster may be
*          as high as ( O(n)*macheps ) / GAP(I).
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          < 0 :  if INFO = -I, the I-th argument had an illegal value
*          > 0 :  if mod(INFO,M+1) = I, then I eigenvectors failed to
*                 converge in MAXITS iterations. Their indices are
*                 stored in the array IFAIL.
*                 if INFO/(M+1) = I, then eigenvectors corresponding to
*                 I clusters of eigenvalues could not be orthogonalized
*                 due to insufficient workspace. The indices of the
*                 clusters are stored in the array ICLUSTR.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, MOD
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, NUMROC
      EXTERNAL           ICEIL, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DGEBR2D, DGEBS2D,
     $                   DLASRT2, DSTEIN2, IGAMN2D, IGEBR2D, IGEBS2D,
     $                   PCHK1MAT, PXERBLA, PZLAEVSWP
*     ..
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, NEGONE, ODM1, FIVE, ODM3, ODM18
      PARAMETER          ( ZERO = 0.0D+0, NEGONE = -1.0D+0,
     $                   ODM1 = 1.0D-1, FIVE = 5.0D+0, ODM3 = 1.0D-3,
     $                   ODM18 = 1.0D-18 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SORTED
      INTEGER            B1, BN, BNDRY, CLSIZ, COL, I, IFIRST, IINFO,
     $                   ILAST, IM, INDRW, ITMP, J, K, LGCLSIZ, LLWORK,
     $                   LOAD, LOCINFO, MAXVEC, MQ00, MYCOL, MYROW,
     $                   NBLK, NERR, NEXT, NP00, NPCOL, NPROW, NVS,
     $                   OLNBLK, P, ROW, SELF, TILL, TOTERR
      DOUBLE PRECISION   DIFF, MINGAP, ONENRM, ORGFAC, ORTOL, TMPFAC
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      CALL BLACS_GRIDINFO( DESCZ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      SELF = MYROW*NPCOL + MYCOL
*
*     Make sure that we belong to this context (before calling PCHK1MAT)
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 1200+CTXT_ )
      ELSE
*
*        Make sure that NPROW>0 and NPCOL>0 before calling NUMROC
*
         CALL CHK1MAT( N, 1, N, 1, IZ, JZ, DESCZ, 12, INFO )
         IF( INFO.EQ.0 ) THEN
*
*           Now we know that our context is good enough to
*           perform the rest of the checks
*
            NP00 = NUMROC( N, DESCZ( MB_ ), 0, 0, NPROW )
            MQ00 = NUMROC( M, DESCZ( NB_ ), 0, 0, NPCOL )
            P = NPROW*NPCOL
*
*           Compute the maximum number of vectors per process
*
            LLWORK = LWORK
            CALL IGAMN2D( DESCZ( CTXT_ ), 'A', ' ', 1, 1, LLWORK, 1, 1,
     $                    1, -1, -1, -1 )
            INDRW = MAX( 5*N, NP00*MQ00 )
            IF( N.NE.0 )
     $         MAXVEC = ( LLWORK-INDRW ) / N
            LOAD = ICEIL( M, P )
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
               TMPFAC = ORFAC
               CALL DGEBS2D( DESCZ( CTXT_ ), 'ALL', ' ', 1, 1, TMPFAC,
     $                       1 )
            ELSE
               CALL DGEBR2D( DESCZ( CTXT_ ), 'ALL', ' ', 1, 1, TMPFAC,
     $                       1, 0, 0 )
            END IF
*
            LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
            IF( M.LT.0 .OR. M.GT.N ) THEN
               INFO = -4
            ELSE IF( MAXVEC.LT.LOAD .AND. .NOT.LQUERY ) THEN
               INFO = -14
            ELSE IF( LIWORK.LT.3*N+P+1 .AND. .NOT.LQUERY ) THEN
               INFO = -16
            ELSE
               DO 10 I = 2, M
                  IF( IBLOCK( I ).LT.IBLOCK( I-1 ) ) THEN
                     INFO = -6
                     GO TO 20
                  END IF
                  IF( IBLOCK( I ).EQ.IBLOCK( I-1 ) .AND. W( I ).LT.
     $                W( I-1 ) ) THEN
                     INFO = -5
                     GO TO 20
                  END IF
   10          CONTINUE
   20          CONTINUE
               IF( INFO.EQ.0 ) THEN
                  IF( ABS( TMPFAC-ORFAC ).GT.FIVE*ABS( TMPFAC ) )
     $               INFO = -8
               END IF
            END IF
*
         END IF
         IDUM1( 1 ) = M
         IDUM2( 1 ) = 4
         CALL PCHK1MAT( N, 1, N, 1, IZ, JZ, DESCZ, 12, 1, IDUM1, IDUM2,
     $                  INFO )
         WORK( 1 ) = DBLE( MAX( 5*N, NP00*MQ00 )+ICEIL( M, P )*N )
         IWORK( 1 ) = 3*N + P + 1
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCZ( CTXT_ ), 'PZSTEIN', -INFO )
         RETURN
      ELSE IF( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 ) THEN
         RETURN
      END IF
*
      DO 30 I = 1, M
         IFAIL( I ) = 0
   30 CONTINUE
      DO 40 I = 1, P + 1
         IWORK( I ) = 0
   40 CONTINUE
      DO 50 I = 1, P
         GAP( I ) = NEGONE
         ICLUSTR( 2*I-1 ) = 0
         ICLUSTR( 2*I ) = 0
   50 CONTINUE
*
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
      IF( ORFAC.GE.ZERO ) THEN
         TMPFAC = ORFAC
      ELSE
         TMPFAC = ODM3
      END IF
      ORGFAC = TMPFAC
*
*     Allocate the work among the processes
*
      ILAST = M / LOAD
      IF( MOD( M, LOAD ).EQ.0 )
     $   ILAST = ILAST - 1
      OLNBLK = -1
      NVS = 0
      NEXT = 1
      IM = 0
      ONENRM = ZERO
      DO 100 I = 0, ILAST - 1
         NEXT = NEXT + LOAD
         J = NEXT - 1
         IF( J.GT.NVS ) THEN
            NBLK = IBLOCK( NEXT )
            IF( NBLK.EQ.IBLOCK( NEXT-1 ) .AND. NBLK.NE.OLNBLK ) THEN
*
*              Compute orthogonalization criterion
*
               IF( NBLK.EQ.1 ) THEN
                  B1 = 1
               ELSE
                  B1 = ISPLIT( NBLK-1 ) + 1
               END IF
               BN = ISPLIT( NBLK )
*
               ONENRM = ABS( D( B1 ) ) + ABS( E( B1 ) )
               ONENRM = MAX( ONENRM, ABS( D( BN ) )+ABS( E( BN-1 ) ) )
               DO 60 J = B1 + 1, BN - 1
                  ONENRM = MAX( ONENRM, ABS( D( J ) )+ABS( E( J-1 ) )+
     $                     ABS( E( J ) ) )
   60          CONTINUE
               OLNBLK = NBLK
            END IF
            TILL = NVS + MAXVEC
   70       CONTINUE
            J = NEXT - 1
            IF( TMPFAC.GT.ODM18 ) THEN
               ORTOL = TMPFAC*ONENRM
               DO 80 J = NEXT - 1, MIN( TILL, M-1 )
                  IF( IBLOCK( J+1 ).NE.IBLOCK( J ) .OR. W( J+1 )-
     $                W( J ).GE.ORTOL ) THEN
                     GO TO 90
                  END IF
   80          CONTINUE
               IF( J.EQ.M .AND. TILL.GE.M )
     $            GO TO 90
               TMPFAC = TMPFAC*ODM1
               GO TO 70
            END IF
   90       CONTINUE
            J = MIN( J, TILL )
         END IF
         IF( SELF.EQ.I )
     $      IM = MAX( 0, J-NVS )
*
         IWORK( I+1 ) = NVS
         NVS = MAX( J, NVS )
  100 CONTINUE
      IF( SELF.EQ.ILAST )
     $   IM = M - NVS
      IWORK( ILAST+1 ) = NVS
      DO 110 I = ILAST + 2, P + 1
         IWORK( I ) = M
  110 CONTINUE
*
      CLSIZ = 1
      LGCLSIZ = 1
      ILAST = 0
      NBLK = 0
      BNDRY = 2
      K = 1
      DO 140 I = 1, M
         IF( IBLOCK( I ).NE.NBLK ) THEN
            NBLK = IBLOCK( I )
            IF( NBLK.EQ.1 ) THEN
               B1 = 1
            ELSE
               B1 = ISPLIT( NBLK-1 ) + 1
            END IF
            BN = ISPLIT( NBLK )
*
            ONENRM = ABS( D( B1 ) ) + ABS( E( B1 ) )
            ONENRM = MAX( ONENRM, ABS( D( BN ) )+ABS( E( BN-1 ) ) )
            DO 120 J = B1 + 1, BN - 1
               ONENRM = MAX( ONENRM, ABS( D( J ) )+ABS( E( J-1 ) )+
     $                  ABS( E( J ) ) )
  120       CONTINUE
*
         END IF
         IF( I.GT.1 ) THEN
            DIFF = W( I ) - W( I-1 )
            IF( IBLOCK( I ).NE.IBLOCK( I-1 ) .OR. I.EQ.M .OR. DIFF.GT.
     $          ORGFAC*ONENRM ) THEN
               IFIRST = ILAST
               IF( I.EQ.M ) THEN
                  IF( IBLOCK( M ).NE.IBLOCK( M-1 ) .OR. DIFF.GT.ORGFAC*
     $                ONENRM ) THEN
                     ILAST = M - 1
                  ELSE
                     ILAST = M
                  END IF
               ELSE
                  ILAST = I - 1
               END IF
               CLSIZ = ILAST - IFIRST
               IF( CLSIZ.GT.1 ) THEN
                  IF( LGCLSIZ.LT.CLSIZ )
     $               LGCLSIZ = CLSIZ
                  MINGAP = ONENRM
  130             CONTINUE
                  IF( BNDRY.GT.P+1 )
     $               GO TO 150
                  IF( IWORK( BNDRY ).GT.IFIRST .AND. IWORK( BNDRY ).LT.
     $                ILAST ) THEN
                     MINGAP = MIN( W( IWORK( BNDRY )+1 )-
     $                        W( IWORK( BNDRY ) ), MINGAP )
                  ELSE IF( IWORK( BNDRY ).GE.ILAST ) THEN
                     IF( MINGAP.LT.ONENRM ) THEN
                        ICLUSTR( 2*K-1 ) = IFIRST + 1
                        ICLUSTR( 2*K ) = ILAST
                        GAP( K ) = MINGAP / ONENRM
                        K = K + 1
                     END IF
                     GO TO 140
                  END IF
                  BNDRY = BNDRY + 1
                  GO TO 130
               END IF
            END IF
         END IF
  140 CONTINUE
  150 CONTINUE
      INFO = ( K-1 )*( M+1 )
*
*     Call DSTEIN2 to find the eigenvectors
*
      CALL DSTEIN2( N, D, E, IM, W( IWORK( SELF+1 )+1 ),
     $              IBLOCK( IWORK( SELF+1 )+1 ), ISPLIT, ORGFAC,
     $              WORK( INDRW+1 ), N, WORK, IWORK( P+2 ),
     $              IFAIL( IWORK( SELF+1 )+1 ), LOCINFO )
*
*     Redistribute the eigenvector matrix to conform with the block
*     cyclic distribution of the input matrix
*
*
      DO 160 I = 1, M
         IWORK( P+1+I ) = I
  160 CONTINUE
*
      CALL DLASRT2( 'I', M, W, IWORK( P+2 ), IINFO )
*
      DO 170 I = 1, M
         IWORK( M+P+1+IWORK( P+1+I ) ) = I
  170 CONTINUE
*
*
      DO 180 I = 1, LOCINFO
         ITMP = IWORK( SELF+1 ) + I
         IFAIL( ITMP ) = IFAIL( ITMP ) + ITMP - I
         IFAIL( ITMP ) = IWORK( M+P+1+IFAIL( ITMP ) )
  180 CONTINUE
*
      DO 190 I = 1, K - 1
         ICLUSTR( 2*I-1 ) = IWORK( M+P+1+ICLUSTR( 2*I-1 ) )
         ICLUSTR( 2*I ) = IWORK( M+P+1+ICLUSTR( 2*I ) )
  190 CONTINUE
*
*
* Still need to apply the above permutation to IFAIL
*
*
      TOTERR = 0
      DO 210 I = 1, P
         IF( SELF.EQ.I-1 ) THEN
            CALL IGEBS2D( DESCZ( CTXT_ ), 'ALL', ' ', 1, 1, LOCINFO, 1 )
            IF( LOCINFO.NE.0 ) THEN
               CALL IGEBS2D( DESCZ( CTXT_ ), 'ALL', ' ', LOCINFO, 1,
     $                       IFAIL( IWORK( I )+1 ), LOCINFO )
               DO 200 J = 1, LOCINFO
                  IFAIL( TOTERR+J ) = IFAIL( IWORK( I )+J )
  200          CONTINUE
               TOTERR = TOTERR + LOCINFO
            END IF
         ELSE
*
            ROW = ( I-1 ) / NPCOL
            COL = MOD( I-1, NPCOL )
*
            CALL IGEBR2D( DESCZ( CTXT_ ), 'ALL', ' ', 1, 1, NERR, 1,
     $                    ROW, COL )
            IF( NERR.NE.0 ) THEN
               CALL IGEBR2D( DESCZ( CTXT_ ), 'ALL', ' ', NERR, 1,
     $                       IFAIL( TOTERR+1 ), NERR, ROW, COL )
               TOTERR = TOTERR + NERR
            END IF
         END IF
  210 CONTINUE
      INFO = INFO + TOTERR
*
*
      CALL PZLAEVSWP( N, WORK( INDRW+1 ), N, Z, IZ, JZ, DESCZ, IWORK,
     $                IWORK( M+P+2 ), WORK, INDRW )
*
      DO 220 I = 2, P
         IWORK( I ) = IWORK( M+P+1+IWORK( I ) )
  220 CONTINUE
*
*
*     Sort the IWORK array
*
*
  230 CONTINUE
      SORTED = .TRUE.
      DO 240 I = 2, P - 1
         IF( IWORK( I ).GT.IWORK( I+1 ) ) THEN
            ITMP = IWORK( I+1 )
            IWORK( I+1 ) = IWORK( I )
            IWORK( I ) = ITMP
            SORTED = .FALSE.
         END IF
  240 CONTINUE
      IF( .NOT.SORTED )
     $   GO TO 230
*
      DO 250 I = P + 1, 1, -1
         IWORK( I+1 ) = IWORK( I )
  250 CONTINUE
*
      WORK( 1 ) = ( LGCLSIZ+LOAD-1 )*N + INDRW
      IWORK( 1 ) = 3*N + P + 1
*
*     End of PZSTEIN
*
      END
