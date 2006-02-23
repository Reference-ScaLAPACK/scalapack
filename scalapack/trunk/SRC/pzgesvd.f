
      SUBROUTINE PZGESVD(JOBU,JOBVT,M,N,A,IA,JA,DESCA,S,U,IU,JU,DESCU,
     +                   VT,IVT,JVT,DESCVT,WORK,LWORK,RWORK,INFO)
*
*  -- ScaLAPACK routine (version 1.7) --
*     Univ. of Tennessee, Oak Ridge National Laboratory
*     and Univ. of California Berkeley.
*     Jan 2006

*
*     .. Scalar Arguments ..
      CHARACTER JOBU,JOBVT
      INTEGER IA,INFO,IU,IVT,JA,JU,JVT,LWORK,M,N
*     ..
*     .. Array Arguments ..
      INTEGER DESCA(*),DESCU(*),DESCVT(*)
      COMPLEX*16 A(*),U(*),VT(*),WORK(*)
      DOUBLE PRECISION S(*)
      DOUBLE PRECISION RWORK(*)
*     ..
*
*  Purpose
*  =======
*
*  PZGESVD computes the singular value decomposition (SVD) of an
*  M-by-N matrix A, optionally computing the left and/or right
*  singular vectors. The SVD is written as
*
*       A = U * SIGMA * transpose(V)
*
*  where SIGMA is an M-by-N matrix which is zero except for its
*  min(M,N) diagonal elements, U is an M-by-M orthogonal matrix, and
*  V is an N-by-N orthogonal matrix. The diagonal elements of SIGMA
*  are the singular values of A and the columns of U and V are the
*  corresponding right and left singular vectors, respectively. The
*  singular values are returned in array S in decreasing order and
*  only the first min(M,N) columns of U and rows of VT = V**T are
*  computed.
*
*  Notes
*  =====
*  Each global data object is described by an associated description
*  vector. This vector stores the information required to establish
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
*  Let K be the number of rows or columns of a distributed matrix, and
*  assume that its process grid has dimension r x c. LOCr( K ) denotes
*  the number of elements of K that a process would receive if K were
*  distributed over the r processes of its process column. Similarly,
*  LOCc( K ) denotes the number of elements of K that a process would
*  receive if K were distributed over the c processes of its process
*  row. The values of LOCr() and LOCc() may be determined via a call
*  to the ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*          MP = number of local rows in A and U
*          NQ = number of local columns in A and VT
*          SIZE = min( M, N )
*          SIZEQ = number of local columns in U
*          SIZEP = number of local rows in VT
*
*  JOBU    (global input) CHARACTER*1
*          Specifies options for computing U:
*          = 'V':  the first SIZE columns of U (the left singular
*                  vectors) are returned in the array U;
*          = 'N':  no columns of U (no left singular vectors) are
*                  computed.
*
*  JOBVT   (global input) CHARACTER*1
*          Specifies options for computing V**T:
*          = 'V':  the first SIZE rows of V**T (the right singular
*                  vectors) are returned in the array VT;
*          = 'N':  no rows of V**T (no right singular vectors) are
*                  computed.
*
*  M       (global input) INTEGER
*          The number of rows of the input matrix A.  M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns of the input matrix A.  N >= 0.
*
*  A       (local input/workspace) block cyclic COMPLEX*16
*          array,
*          global dimension (M, N), local dimension (MP, NQ)
*          On exit, the contents of A are destroyed.
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global input) INTEGER array of dimension DLEN_
*          The array descriptor for the distributed matrix A.
*
*  S       (global output) DOUBLE PRECISION   array, dimension SIZE
*          The singular values of A, sorted so that S(i) >= S(i+1).
*
*  U       (local output) COMPLEX*16         array, local dimension
*          (MP, SIZEQ), global dimension (M, SIZE)
*          if JOBU = 'V', U contains the first min(m,n) columns of U
*          if JOBU = 'N', U is not referenced.
*
*  IU      (global input) INTEGER
*          The row index in the global array U indicating the first
*          row of sub( U ).
*
*  JU      (global input) INTEGER
*          The column index in the global array U indicating the
*          first column of sub( U ).
*
*  DESCU   (global input) INTEGER array of dimension DLEN_
*          The array descriptor for the distributed matrix U.
*
*  VT      (local output) COMPLEX*16         array, local dimension
*          (SIZEP, NQ), global dimension (SIZE, N).
*          If JOBVT = 'V', VT contains the first SIZE rows of
*          V**T. If JOBVT = 'N', VT is not referenced.
*
*  IVT     (global input) INTEGER
*          The row index in the global array VT indicating the first
*          row of sub( VT ).
*
*  JVT     (global input) INTEGER
*          The column index in the global array VT indicating the
*          first column of sub( VT ).
*
*  DESCVT   (global input) INTEGER array of dimension DLEN_
*          The array descriptor for the distributed matrix VT.
*
*  WORK    (local workspace/output) COMPLEX*16         array, dimension
*          (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (local input) INTEGER
*          The dimension of the array WORK.
*
*          LWORK >= 1 + 2*SIZEB + MAX(WATOBD, WBDTOSVD),
*
*          where SIZEB = MAX(M,N), and WATOBD and WBDTOSVD refer,
*          respectively, to the workspace required to bidiagonalize
*          the matrix A and to go from the bidiagonal matrix to the
*          singular value decomposition U*S*VT.
*
*          For WATOBD, the following holds:
*
*          WATOBD = MAX(MAX(WPZLANGE,WPZGEBRD),
*                       MAX(WPZLARED2D,WP(pre)LARED1D)),
*
*          where WPZLANGE, WPZLARED1D, WPZLARED2D, WPZGEBRD are the
*          workspaces required respectively for the subprograms
*          PZLANGE, PDLARED1D, PDLARED2D, PZGEBRD. Using the
*          standard notation
*
*          MP = NUMROC( M, MB, MYROW, DESCA( CTXT_ ), NPROW),
*          NQ = NUMROC( N, NB, MYCOL, DESCA( LLD_ ), NPCOL),
*
*          the workspaces required for the above subprograms are
*
*          WPZLANGE = MP,
*          WPDLARED1D = NQ0,
*          WPDLARED2D = MP0,
*          WPZGEBRD = NB*(MP + NQ + 1) + NQ,
*
*          where NQ0 and MP0 refer, respectively, to the values obtained
*          at MYCOL = 0 and MYROW = 0. In general, the upper limit for
*          the workspace is given by a workspace required on
*          processor (0,0):
*
*          WATOBD <= NB*(MP0 + NQ0 + 1) + NQ0.
*
*          In case of a homogeneous process grid this upper limit can
*          be used as an estimate of the minimum workspace for every
*          processor.
*
*          For WBDTOSVD, the following holds:
*
*          WBDTOSVD = SIZE*(WANTU*NRU + WANTVT*NCVT) +
*                     MAX(WZBDSQR,
*                         MAX(WANTU*WPZORMBRQLN, WANTVT*WPZORMBRPRT)),
*
*          where
*
*                          1, if left(right) singular vectors are wanted
*          WANTU(WANTVT) =
*                          0, otherwise
*
*          and WZBDSQR, WPZORMBRQLN and WPZORMBRPRT refer respectively
*          to the workspace required for the subprograms ZBDSQR,
*          PZUNMBR(QLN), and PZUNMBR(PRT), where QLN and PRT are the
*          values of the arguments VECT, SIDE, and TRANS in the call
*          to PZUNMBR. NRU is equal to the local number of rows of
*          the matrix U when distributed 1-dimensional "column" of
*          processes. Analogously, NCVT is equal to the local number
*          of columns of the matrix VT when distributed across
*          1-dimensional "row" of processes. Calling the LAPACK
*          procedure ZBDSQR requires
*
*          WZBDSQR = MAX(1, 4*SIZE )
*
*          on every processor. Finally,
*
*          WPZORMBRQLN = MAX( (NB*(NB-1))/2, (SIZEQ+MP)*NB)+NB*NB,
*          WPZORMBRPRT = MAX( (MB*(MB-1))/2, (SIZEP+NQ)*MB )+MB*MB,
*
*          If LWORK = -1, then LWORK is global input and a workspace
*          query is assumed; the routine only calculates the minimum
*          size for the work array. The required workspace is returned
*          as the first element of WORK and no error message is issued
*          by PXERBLA.
*
*  RWORK   (workspace) REAL             array, dimension (1+4*SIZEB)
*          On exit, if INFO = 0, RWORK(1) returns the necessary size
*          for RWORK.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value

*          > 0:  if ZBDSQR did not converge
*                If INFO = MIN(M,N) + 1, then PZGESVD has detected
*                heterogeneity by finding that eigenvalues were not
*                identical across the process grid. In this case, the
*                accuracy of the results from PZGESVD cannot be
*                guaranteed.
*
*  =====================================================================
*
*  The results of PZGEBRD, and therefore PZGESVD, may vary slightly
*  from run to run with the same input data. If repeatability is an
*  issue, call BLACS_SET with the appropriate option after defining
*  the process grid.
*
*  Alignment requirements
*  ======================
*
*  The routine PZGESVD inherits the same alignement requirement as
*  the routine PZGEBRD, namely:
*
*  The distributed submatrix sub( A ) must verify some alignment proper-
*  ties, namely the following expressions should be true:
*  ( MB_A.EQ.NB_A .AND. IROFFA.EQ.ICOFFA )
*          where NB = MB_A = NB_A,
*          IROFFA = MOD( IA-1, NB ), ICOFFA = MOD( JA-1, NB ),
*
*  =====================================================================
*
*
*     .. Parameters ..
      INTEGER BLOCK_CYCLIC_2D,DLEN_,DTYPE_,CTXT_,M_,N_,MB_,NB_,RSRC_,
     +        CSRC_,LLD_,ITHVAL
      PARAMETER (BLOCK_CYCLIC_2D=1,DLEN_=9,DTYPE_=1,CTXT_=2,M_=3,N_=4,
     +          MB_=5,NB_=6,RSRC_=7,CSRC_=8,LLD_=9,ITHVAL=10)
      COMPLEX*16 ZERO,ONE
      PARAMETER (ZERO= ((0.0D+0,0.0D+0)),ONE= ((1.0D+0,0.0D+0)))
      DOUBLE PRECISION DZERO,DONE
      PARAMETER (DZERO=0.0D+0,DONE=1.0D+0)
*     ..
*     .. Local Scalars ..
      CHARACTER UPLO
      INTEGER CONTEXTC,CONTEXTR,I,INDD,INDD2,INDE,INDE2,INDTAUP,INDTAUQ,
     +        INDU,INDV,INDWORK,IOFFD,IOFFE,ISCALE,J,K,LDU,LDVT,LLWORK,
     +        LWMIN,MAXIM,MB,MP,MYPCOL,MYPCOLC,MYPCOLR,MYPROW,MYPROWC,
     +        MYPROWR,NB,NCVT,NPCOL,NPCOLC,NPCOLR,NPROCS,NPROW,NPROWC,
     +        NPROWR,NQ,NRU,SIZE,SIZEB,SIZEP,SIZEPOS,SIZEQ,WANTU,WANTVT,
     +        WATOBD,WBDTOSVD,WZBDSQR,WPZGEBRD,WPZLANGE,WPZORMBRPRT,
     +        WPZORMBRQLN
      DOUBLE PRECISION ANRM,BIGNUM,EPS,RMAX,RMIN,SAFMIN,SIGMA,SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER DESCTU(DLEN_),DESCTVT(DLEN_),IDUM1(3),IDUM2(3)
      DOUBLE PRECISION C(1,1)
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      INTEGER NUMROC
      DOUBLE PRECISION PDLAMCH,PZLANGE
      EXTERNAL LSAME,NUMROC,PDLAMCH,PZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GET,BLACS_GRIDEXIT,BLACS_GRIDINFO,BLACS_GRIDINIT,
     +         CHK1MAT,ZBDSQR,DESCINIT,DGAMN2D,DGAMX2D,DSCAL,IGAMX2D,
     +         IGEBR2D,IGEBS2D,PCHK1MAT,PZGEBRD,PZGEMR2D,PDLARED1D,
     +         PDLARED2D,PZLASCL,PZLASET,PZUNMBR,PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN,SQRT,DBLE
      INTRINSIC DCMPLX
*     ..
*     .. Executable Statements ..
*     This is just to keep ftnchek happy
      IF (BLOCK_CYCLIC_2D*DTYPE_*LLD_*MB_*M_*NB_*N_.LT.0) RETURN
*
      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYPROW,MYPCOL)
      ISCALE = 0
      INFO = 0
*
      IF (NPROW.EQ.-1) THEN
          INFO = - (800+CTXT_)
      ELSE
*
          SIZE = MIN(M,N)
          SIZEB = MAX(M,N)
          NPROCS = NPROW*NPCOL
          IF (M.GE.N) THEN
              IOFFD = JA - 1
              IOFFE = IA - 1
              SIZEPOS = 1
          ELSE
              IOFFD = IA - 1
              IOFFE = JA - 1
              SIZEPOS = 3
          END IF
*
          IF (LSAME(JOBU,'V')) THEN
              WANTU = 1
          ELSE
              WANTU = 0
          END IF
          IF (LSAME(JOBVT,'V')) THEN
              WANTVT = 1
          ELSE
              WANTVT = 0
          END IF
*
          CALL CHK1MAT(M,3,N,4,IA,JA,DESCA,8,INFO)
          IF (WANTU.EQ.1) THEN
              CALL CHK1MAT(M,3,SIZE,SIZEPOS,IU,JU,DESCU,13,INFO)
          END IF
          IF (WANTVT.EQ.1) THEN
              CALL CHK1MAT(SIZE,SIZEPOS,N,4,IVT,JVT,DESCVT,17,INFO)
          END IF
          CALL IGAMX2D(DESCA(CTXT_),'A',' ',1,1,INFO,1,1,1,-1,-1,0)
*
          IF (INFO.EQ.0) THEN
*
*           Set up pointers into the WORK array.
*
              INDD = 2
              INDE = INDD + SIZEB + IOFFD
              INDD2 = INDE + SIZEB + IOFFE
              INDE2 = INDD2 + SIZEB + IOFFD
*
              INDTAUQ = 2
              INDTAUP = INDTAUQ + SIZEB + JA - 1
              INDWORK = INDTAUP + SIZEB + IA - 1
              LLWORK = LWORK - INDWORK + 1
*
*           Initialize contexts for "column" and "row" process matrices.
*
              CALL BLACS_GET(DESCA(CTXT_),10,CONTEXTC)
              CALL BLACS_GRIDINIT(CONTEXTC,'R',NPROCS,1)
              CALL BLACS_GRIDINFO(CONTEXTC,NPROWC,NPCOLC,MYPROWC,
     +                            MYPCOLC)
              CALL BLACS_GET(DESCA(CTXT_),10,CONTEXTR)
              CALL BLACS_GRIDINIT(CONTEXTR,'R',1,NPROCS)
              CALL BLACS_GRIDINFO(CONTEXTR,NPROWR,NPCOLR,MYPROWR,
     +                            MYPCOLR)
*
*           Set local dimensions of matrices (this is for MB=NB=1).
*
              NRU = NUMROC(M,1,MYPROWC,0,NPROCS)
              NCVT = NUMROC(N,1,MYPCOLR,0,NPROCS)
              NB = DESCA(NB_)
              MB = DESCA(MB_)
              MP = NUMROC(M,MB,MYPROW,DESCA(RSRC_),NPROW)
              NQ = NUMROC(N,NB,MYPCOL,DESCA(CSRC_),NPCOL)
              IF (WANTVT.EQ.1) THEN
                  SIZEP = NUMROC(SIZE,DESCVT(MB_),MYPROW,DESCVT(RSRC_),
     +                    NPROW)
              ELSE
                  SIZEP = 0
              END IF
              IF (WANTU.EQ.1) THEN
                  SIZEQ = NUMROC(SIZE,DESCU(NB_),MYPCOL,DESCU(CSRC_),
     +                    NPCOL)
              ELSE
                  SIZEQ = 0
              END IF
*
*           Transmit MAX(NQ0, MP0).
*
              IF (MYPROW.EQ.0 .AND. MYPCOL.EQ.0) THEN
                  MAXIM = MAX(NQ,MP)
                  CALL IGEBS2D(DESCA(CTXT_),'All',' ',1,1,MAXIM,1)
              ELSE
                  CALL IGEBR2D(DESCA(CTXT_),'All',' ',1,1,MAXIM,1,0,0)
              END IF
*
              WPZLANGE = MP
              WPZGEBRD = NB* (MP+NQ+1) + NQ
              WATOBD = MAX(MAX(WPZLANGE,WPZGEBRD),MAXIM)
*
              WZBDSQR = MAX(1,4*SIZE)
              WPZORMBRQLN = MAX((NB* (NB-1))/2, (SIZEQ+MP)*NB) + NB*NB
              WPZORMBRPRT = MAX((MB* (MB-1))/2, (SIZEP+NQ)*MB) + MB*MB
              WBDTOSVD = SIZE* (WANTU*NRU+WANTVT*NCVT) +
     +                   MAX(WZBDSQR,MAX(WANTU*WPZORMBRQLN,
     +                   WANTVT*WPZORMBRPRT))
*
*           Finally, calculate required workspace.
*
              LWMIN = 1 + 2*SIZEB + MAX(WATOBD,WBDTOSVD)
              WORK(1) = DCMPLX(LWMIN,0D+00)
              RWORK(1) = DBLE(1+4*SIZEB)
*
              IF (WANTU.NE.1 .AND. .NOT. (LSAME(JOBU,'N'))) THEN
                  INFO = -1
              ELSE IF (WANTVT.NE.1 .AND. .NOT. (LSAME(JOBVT,'N'))) THEN
                  INFO = -2
              ELSE IF (LWORK.LT.LWMIN .AND. LWORK.NE.-1) THEN
                  INFO = -19
              END IF
*
          END IF
*
          IDUM1(1) = WANTU
          IDUM1(2) = WANTVT
          IF (LWORK.EQ.-1) THEN
              IDUM1(3) = -1
          ELSE
              IDUM1(3) = 1
          END IF
          IDUM2(1) = 1
          IDUM2(2) = 2
          IDUM2(3) = 19
          CALL PCHK1MAT(M,3,N,4,IA,JA,DESCA,8,3,IDUM1,IDUM2,INFO)
          IF (INFO.EQ.0) THEN
              IF (WANTU.EQ.1) THEN
                  CALL PCHK1MAT(M,3,SIZE,4,IU,JU,DESCU,13,0,IDUM1,IDUM2,
     +                          INFO)
              END IF
              IF (WANTVT.EQ.1) THEN
                  CALL PCHK1MAT(SIZE,3,N,4,IVT,JVT,DESCVT,17,0,IDUM1,
     +                          IDUM2,INFO)
              END IF
          END IF
*
      END IF
*
      IF (INFO.NE.0) THEN
          CALL PXERBLA(DESCA(CTXT_),'PZGESVD',-INFO)
          RETURN
      ELSE IF (LWORK.EQ.-1) THEN
          GO TO 40
      END IF
*
*     Quick return if possible.
*
      IF (M.LE.0 .OR. N.LE.0) GO TO 40
*
*     Get machine constants.
*
      SAFMIN = PDLAMCH(DESCA(CTXT_),'Safe minimum')
      EPS = PDLAMCH(DESCA(CTXT_),'Precision')
      SMLNUM = SAFMIN/EPS
      BIGNUM = DONE/SMLNUM
      RMIN = SQRT(SMLNUM)
      RMAX = MIN(SQRT(BIGNUM),DONE/SQRT(SQRT(SAFMIN)))
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = PZLANGE('1',M,N,A,IA,JA,DESCA,WORK(INDWORK))
      IF (ANRM.GT.DZERO .AND. ANRM.LT.RMIN) THEN
          ISCALE = 1
          SIGMA = RMIN/ANRM
      ELSE IF (ANRM.GT.RMAX) THEN
          ISCALE = 1
          SIGMA = RMAX/ANRM
      END IF
*
      IF (ISCALE.EQ.1) THEN
          CALL PZLASCL('G',DONE,SIGMA,M,N,A,IA,JA,DESCA,INFO)
      END IF
*
      CALL PZGEBRD(M,N,A,IA,JA,DESCA,RWORK(INDD),RWORK(INDE),
     +             WORK(INDTAUQ),WORK(INDTAUP),WORK(INDWORK),LLWORK,
     +             INFO)
*
*     Copy D and E to all processes.
*     Array D is in local array of dimension:
*     LOCc(JA+MIN(M,N)-1) if M >= N; LOCr(IA+MIN(M,N)-1) otherwise.
*     Array E is in local array of dimension
*     LOCr(IA+MIN(M,N)-1) if M >= N; LOCc(JA+MIN(M,N)-2) otherwise.
*
      IF (M.GE.N) THEN
*        Distribute D
          CALL PDLARED1D(N+IOFFD,IA,JA,DESCA,RWORK(INDD),RWORK(INDD2),
     +                   WORK(INDWORK),LLWORK)
*        Distribute E
          CALL PDLARED2D(M+IOFFE,IA,JA,DESCA,RWORK(INDE),RWORK(INDE2),
     +                   WORK(INDWORK),LLWORK)
      ELSE
*        Distribute D
          CALL PDLARED2D(M+IOFFD,IA,JA,DESCA,RWORK(INDD),RWORK(INDD2),
     +                   WORK(INDWORK),LLWORK)
*        Distribute E
          CALL PDLARED1D(N+IOFFE,IA,JA,DESCA,RWORK(INDE),RWORK(INDE2),
     +                   WORK(INDWORK),LLWORK)
      END IF
*
*     Prepare for calling PZBDSQR.
*
      IF (M.GE.N) THEN
          UPLO = 'U'
      ELSE
          UPLO = 'L'
      END IF
*
      INDU = INDWORK
      INDV = INDU + SIZE*NRU*WANTU
      INDWORK = INDV + SIZE*NCVT*WANTVT
*
      LDU = MAX(1,NRU)
      LDVT = MAX(1,SIZE)
*
      CALL DESCINIT(DESCTU,M,SIZE,1,1,0,0,CONTEXTC,LDU,INFO)
      CALL DESCINIT(DESCTVT,SIZE,N,1,1,0,0,CONTEXTR,LDVT,INFO)
*
      IF (WANTU.EQ.1) THEN
          CALL PZLASET('Full',M,SIZE,ZERO,ONE,WORK(INDU),1,1,DESCTU)
      ELSE
          NRU = 0
      END IF
*
      IF (WANTVT.EQ.1) THEN
          CALL PZLASET('Full',SIZE,N,ZERO,ONE,WORK(INDV),1,1,DESCTVT)
      ELSE
          NCVT = 0
      END IF
*
      CALL ZBDSQR(UPLO,SIZE,NCVT,NRU,0,RWORK(INDD2+IOFFD),
     +            RWORK(INDE2+IOFFE),WORK(INDV),SIZE,WORK(INDU),LDU,C,1,
     +            WORK(INDWORK),INFO)
*
*     Redistribute elements of U and VT in the block-cyclic fashion.
*
      IF (WANTU.EQ.1) CALL PZGEMR2D(M,SIZE,WORK(INDU),1,1,DESCTU,U,IU,
     +                              JU,DESCU,DESCU(CTXT_))
*
      IF (WANTVT.EQ.1) CALL PZGEMR2D(SIZE,N,WORK(INDV),1,1,DESCTVT,VT,
     +                               IVT,JVT,DESCVT,DESCVT(CTXT_))
*
*     Set to ZERO "non-square" elements of the larger matrices U, VT.
*
      IF (M.GT.N .AND. WANTU.EQ.1) THEN
          CALL PZLASET('Full',M-SIZE,SIZE,ZERO,ZERO,U,IA+SIZE,JU,DESCU)
      ELSE IF (N.GT.M .AND. WANTVT.EQ.1) THEN
          CALL PZLASET('Full',SIZE,N-SIZE,ZERO,ZERO,VT,IVT,JVT+SIZE,
     +                 DESCVT)
      END IF
*
*     Multiply Householder rotations from bidiagonalized matrix.
*
      IF (WANTU.EQ.1) CALL PZUNMBR('Q','L','N',M,SIZE,N,A,IA,JA,DESCA,
     +                             WORK(INDTAUQ),U,IU,JU,DESCU,
     +                             WORK(INDWORK),LLWORK,INFO)
*
      IF (WANTVT.EQ.1) CALL PZUNMBR('P','R','C',SIZE,N,M,A,IA,JA,DESCA,
     +                              WORK(INDTAUP),VT,IVT,JVT,DESCVT,
     +                              WORK(INDWORK),LLWORK,INFO)
*
*     Copy singular values into output array S.
*
      DO 10 I = 1,SIZE
          S(I) = RWORK(INDD2+IOFFD+I-1)
   10 CONTINUE
*
*     If matrix was scaled, then rescale singular values appropriately.
*
      IF (ISCALE.EQ.1) THEN
          CALL DSCAL(SIZE,ONE/SIGMA,S,1)
      END IF
*
*     Compare every ith eigenvalue, or all if there are only a few,
*     across the process grid to check for heterogeneity.
*
      IF (SIZE.LE.ITHVAL) THEN
          J = SIZE
          K = 1
      ELSE
          J = SIZE/ITHVAL
          K = ITHVAL
      END IF
*
      DO 20 I = 1,J
          RWORK(I+INDE) = S((I-1)*K+1)
          RWORK(I+INDD2) = S((I-1)*K+1)
   20 CONTINUE
*
      CALL DGAMN2D(DESCA(CTXT_),'a',' ',J,1,RWORK(1+INDE),J,1,1,-1,-1,0)
      CALL DGAMX2D(DESCA(CTXT_),'a',' ',J,1,RWORK(1+INDD2),J,1,1,-1,-1,
     +             0)
*
      DO 30 I = 1,J
          IF ((RWORK(I+INDE)-RWORK(I+INDD2)).NE.DZERO) THEN
              INFO = SIZE + 1
          END IF
   30 CONTINUE
*
   40 CONTINUE
*
      CALL BLACS_GRIDEXIT(CONTEXTC)
      CALL BLACS_GRIDEXIT(CONTEXTR)
*
*     End of PZGESVD
*
      RETURN
      END
