      RECURSIVE SUBROUTINE PDLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H,
     $                              DESCH, ILOZ, IHIZ, Z, DESCZ, NS, ND,
     $                              SR, SI, V, DESCV, NH, T, DESCT, NV,
     $                              WV, DESCW, WORK, LWORK, IWORK,
     $                              LIWORK, RECLEVEL )
*
*     Contribution from the Department of Computing Science and HPC2N,
*     Umea University, Sweden
*
*  -- ScaLAPACK auxiliary routine (version 2.0.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     Univ. of Colorado Denver and University of California, Berkeley.
*     January, 2012
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LWORK, N, ND, NH, NS,
     $                   NV, NW, LIWORK, RECLEVEL
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      INTEGER            DESCH( * ), DESCZ( * ), DESCT( * ), DESCV( * ),
     $                   DESCW( * ), IWORK( * )
      DOUBLE PRECISION   H( * ), SI( KBOT ), SR( KBOT ), T( * ),
     $                   V( * ), WORK( * ), WV( * ),
     $                   Z( * )
*     ..
*
*  Purpose
*  =======
*
*  Aggressive early deflation:
*
*  This subroutine accepts as input an upper Hessenberg matrix H and
*  performs an orthogonal similarity transformation designed to detect
*  and deflate fully converged eigenvalues from a trailing principal
*  submatrix.  On output H has been overwritten by a new Hessenberg
*  matrix that is a perturbation of an orthogonal similarity
*  transformation of H.  It is to be hoped that the final version of H
*  has many zero subdiagonal entries.
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
*  WANTT   (global input) LOGICAL
*          If .TRUE., then the Hessenberg matrix H is fully updated
*          so that the quasi-triangular Schur factor may be
*          computed (in cooperation with the calling subroutine).
*          If .FALSE., then only enough of H is updated to preserve
*          the eigenvalues.
*
*  WANTZ   (global input) LOGICAL
*          If .TRUE., then the orthogonal matrix Z is updated so
*          so that the orthogonal Schur factor may be computed
*          (in cooperation with the calling subroutine).
*          If .FALSE., then Z is not referenced.
*
*  N       (global input) INTEGER
*          The order of the matrix H and (if WANTZ is .TRUE.) the
*          order of the orthogonal matrix Z.
*
*  KTOP    (global input) INTEGER
*          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
*          KBOT and KTOP together determine an isolated block
*          along the diagonal of the Hessenberg matrix.
*
*  KBOT    (global input) INTEGER
*          It is assumed without a check that either
*          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
*          determine an isolated block along the diagonal of the
*          Hessenberg matrix.
*
*  NW      (global input) INTEGER
*          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
*
*  H       (local input/output) DOUBLE PRECISION array, dimension
*             (DESCH(LLD_),*)
*          On input the initial N-by-N section of H stores the
*          Hessenberg matrix undergoing aggressive early deflation.
*          On output H has been transformed by an orthogonal
*          similarity transformation, perturbed, and the returned
*          to Hessenberg form that (it is to be hoped) has some
*          zero subdiagonal entries.
*
*  DESCH   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix H.
*
*  ILOZ    (global input) INTEGER
*  IHIZ    (global input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension
*             (DESCH(LLD_),*)
*          IF WANTZ is .TRUE., then on output, the orthogonal
*          similarity transformation mentioned above has been
*          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
*          If WANTZ is .FALSE., then Z is unreferenced.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*  NS      (global output) INTEGER
*          The number of unconverged (ie approximate) eigenvalues
*          returned in SR and SI that may be used as shifts by the
*          calling subroutine.
*
*  ND      (global output) INTEGER
*          The number of converged eigenvalues uncovered by this
*          subroutine.
*
*  SR      (global output) DOUBLE PRECISION array, dimension KBOT
*  SI      (global output) DOUBLE PRECISION array, dimension KBOT
*          On output, the real and imaginary parts of approximate
*          eigenvalues that may be used for shifts are stored in
*          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
*          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
*          The real and imaginary parts of converged eigenvalues
*          are stored in SR(KBOT-ND+1) through SR(KBOT) and
*          SI(KBOT-ND+1) through SI(KBOT), respectively.
*
*  V       (global workspace) DOUBLE PRECISION array, dimension 
*             (DESCV(LLD_),*)
*          An NW-by-NW distributed work array.
*
*  DESCV   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix V.
*
*  NH      (input) INTEGER scalar
*          The number of columns of T.  NH.GE.NW.
*
*  T       (global workspace) DOUBLE PRECISION array, dimension 
*             (DESCV(LLD_),*)
*
*  DESCT   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix T.
*
*  NV      (global input) INTEGER
*          The number of rows of work array WV available for
*          workspace.  NV.GE.NW.
*
*  WV      (global workspace) DOUBLE PRECISION array, dimension 
*             (DESCW(LLD_),*)
*
*  DESCW   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix WV.
*
*  WORK    (local workspace) DOUBLE PRECISION array, dimension LWORK.
*          On exit, WORK(1) is set to an estimate of the optimal value
*          of LWORK for the given values of N, NW, KTOP and KBOT.
*
*  LWORK   (local input) INTEGER
*          The dimension of the work array WORK.  LWORK = 2*NW
*          suffices, but greater efficiency may result from larger
*          values of LWORK.
*
*          If LWORK = -1, then a workspace query is assumed; PDLAQR3
*          only estimates the optimal workspace size for the given
*          values of N, NW, KTOP and KBOT.  The estimate is returned
*          in WORK(1).  No error message related to LWORK is issued
*          by XERBLA.  Neither H nor Z are accessed.
*
*  IWORK   (local workspace) INTEGER array, dimension (LIWORK)
*
*  LIWORK  (local input) INTEGER
*          The length of the workspace array IWORK
*
*  ================================================================
*  Based on contributions by
*        Robert Granat and Meiyue Shao,
*        Department of Computing Science and HPC2N,
*        Umea University, Sweden
*
*  ================================================================
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      INTEGER            RECMAX
      LOGICAL            SORTGRAD
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9, RECMAX = 3,
     $                     SORTGRAD = .FALSE. )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S,
     $                   SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP,
     $                   ELEM, ELEM1, ELEM2, ELEM3, R1, ANORM, RNORM,
     $                   RESAED
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL,
     $                   KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2, LWK3,
     $                   LWKOPT, NMIN, LLDH, LLDZ, LLDT, LLDV, LLDWV,
     $                   ICTXT, NPROW, NMAX, NPCOL, MYROW, MYCOL, NB,
     $                   IROFFH, M, RCOLS, TAUROWS, RROWS, TAUCOLS,
     $                   ITAU, IR, IPW, NPROCS, MLOC, IROFFHH,
     $                   ICOFFHH, HHRSRC, HHCSRC, HHROWS, HHCOLS,
     $                   IROFFZZ, ICOFFZZ, ZZRSRC, ZZCSRC, ZZROWS,
     $                   ZZCOLS, IERR, TZROWS0, TZCOLS0, IERR0, IPT0,
     $                   IPZ0, IPW0, NB2, ROUND, LILST, KK, LILST0,
     $                   IWRK1, RSRC, CSRC, LWK4, LWK5, IWRK2, LWK6,
     $                   LWK7, LWK8, ILWKOPT, TZROWS, TZCOLS, NSEL,
     $                   NPMIN, ICTXT_NEW, MYROW_NEW, MYCOL_NEW
      LOGICAL            BULGE, SORTED, LQUERY
*     ..
*     .. Local Arrays ..
      INTEGER            PAR( 6 ), DESCR( DLEN_ ),
     $                   DESCTAU( DLEN_ ), DESCHH( DLEN_ ),
     $                   DESCZZ( DLEN_ ), DESCTZ0( DLEN_ ),
     $                   PMAP( 64*64 )
      DOUBLE PRECISION   DDUM( 1 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, PDLANGE
      INTEGER            PILAENVX, NUMROC, INDXG2P, ICEIL, BLACS_PNUM
      EXTERNAL           DLAMCH, PILAENVX, NUMROC, INDXG2P, PDLANGE,
     $                   MPI_WTIME, ICEIL, BLACS_PNUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           PDCOPY, PDGEHRD, PDGEMM, DLABAD, PDLACPY,
     $                   PDLAQR1, DLANV2, PDLAQR0, PDLARF, PDLARFG,
     $                   PDLASET, PDTRORD, PDELGET, PDELSET,
     $                   PDLAMVE, BLACS_GRIDINFO, BLACS_GRIDMAP,
     $                   BLACS_GRIDEXIT, PDGEMR2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
      ICTXT = DESCH( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW*NPCOL
*
*     Extract local leading dimensions, blockfactors, offset for
*     keeping the alignment requirements and size of deflation window.
*
      LLDH  = DESCH( LLD_ )
      LLDZ  = DESCZ( LLD_ )
      LLDT  = DESCT( LLD_ )
      LLDV  = DESCV( LLD_ )
      LLDWV = DESCW( LLD_ )
      NB = DESCH( MB_ )
      IROFFH = MOD( KTOP - 1, NB )
      JW = MIN( NW, KBOT-KTOP+1 )
      NSEL = NB+JW
*
*     Extract environment variables for parallel eigenvalue reordering.
*
      PAR(1) = PILAENVX(ICTXT, 17, 'PDLAQR3', 'SV', JW, NB, -1, -1)
      PAR(2) = PILAENVX(ICTXT, 18, 'PDLAQR3', 'SV', JW, NB, -1, -1)
      PAR(3) = PILAENVX(ICTXT, 19, 'PDLAQR3', 'SV', JW, NB, -1, -1)
      PAR(4) = PILAENVX(ICTXT, 20, 'PDLAQR3', 'SV', JW, NB, -1, -1)
      PAR(5) = PILAENVX(ICTXT, 21, 'PDLAQR3', 'SV', JW, NB, -1, -1)
      PAR(6) = PILAENVX(ICTXT, 22, 'PDLAQR3', 'SV', JW, NB, -1, -1)
*
*     Check if workspace query.
*
      LQUERY = LWORK.EQ.-1 .OR. LIWORK.EQ.-1
*
*     Estimate optimal workspace.
*
      IF( JW.LE.2 ) THEN
         LWKOPT = 1
      ELSE
*
*        Workspace query calls to PDGEHRD and PDORMHR.
*
         TAUROWS = NUMROC( 1, 1, MYCOL, DESCV(RSRC_), NPROW )
         TAUCOLS = NUMROC( JW+IROFFH, NB, MYCOL, DESCV(CSRC_),
     $        NPCOL )
         CALL PDGEHRD( JW, 1, JW, T, 1, 1, DESCT, WORK, WORK, -1,
     $        INFO )
         LWK1 = INT( WORK( 1 ) ) + TAUROWS*TAUCOLS
*
*        Workspace query call to PDORMHR.
*
         CALL PDORMHR( 'Right', 'No', JW, JW, 1, JW, T, 1, 1, DESCT,
     $        WORK, V, 1, 1, DESCV, WORK, -1, INFO )
         LWK2 = INT( WORK( 1 ) )
*
*        Workspace query call to PDLAQR0.
*
         NMIN = PILAENVX( ICTXT, 12, 'PDLAQR3', 'SV', JW, 1, JW, LWORK )
         NMAX = ( N-1 ) / 3
         IF( JW+IROFFH.GT.NMIN .AND. JW+IROFFH.LE.NMAX
     $        .AND. RECLEVEL.LT.RECMAX ) THEN
            CALL PDLAQR0( .TRUE., .TRUE., JW+IROFFH, 1+IROFFH,
     $           JW+IROFFH, T, DESCT, SR, SI, 1, JW, V, DESCV,
     $           WORK, -1, IWORK, LIWORK-NSEL, INFQR, 
     $           RECLEVEL+1 )
            LWK3 = INT( WORK( 1 ) )
            IWRK1 = IWORK( 1 )
         ELSE
            RSRC = DESCT( RSRC_ )
            CSRC = DESCT( CSRC_ )
            DESCT( RSRC_ ) = 0
            DESCT( CSRC_ ) = 0
            CALL PDLAQR1( .TRUE., .TRUE., JW+IROFFH, 1, JW+IROFFH, T,
     $           DESCT, SR, SI, 1, JW+IROFFH, V, DESCV, WORK, -1,
     $           IWORK, LIWORK-NSEL, INFQR )
            DESCT( RSRC_ ) = RSRC
            DESCT( CSRC_ ) = CSRC
            LWK3 = INT( WORK( 1 ) )
            IWRK1 = IWORK( 1 )
         END IF
*
*        Workspace in case of alignment problems.
*
         TZROWS0 = NUMROC( JW+IROFFH, NB, MYROW, 0, NPROW )
         TZCOLS0 = NUMROC( JW+IROFFH, NB, MYCOL, 0, NPCOL )
         LWK4 = 2 * TZROWS0*TZCOLS0
*
*        Workspace check for reordering.
*
         CALL PDTRORD( 'Vectors', IWORK, PAR, JW+IROFFH, T, 1, 1,
     $        DESCT, V, 1, 1, DESCV, DDUM, DDUM, MLOC, WORK, -1,
     $        IWORK, LIWORK-NSEL, INFO )
         LWK5 = INT( WORK( 1 ) )
         IWRK2 = IWORK( 1 )
*
*        Extra workspace for reflecting back spike
*        (workspace for PDLARF approximated for simplicity).
*
         RROWS =  NUMROC( N+IROFFH, NB, MYROW, DESCV(RSRC_), NPROW )
         RCOLS =  NUMROC( 1, 1, MYCOL, DESCV(CSRC_), NPCOL )
         LWK6 = RROWS*RCOLS + TAUROWS*TAUCOLS +
     $        2*ICEIL(ICEIL(JW+IROFFH,NB),NPROW)*NB
     $         *ICEIL(ICEIL(JW+IROFFH,NB),NPCOL)*NB
*
*        Extra workspace needed by PBLAS update calls
*        (also estimated for simplicity).
*
         LWK7 = MAX( ICEIL(ICEIL(JW,NB),NPROW)*NB *
     $               ICEIL(ICEIL(N-KBOT,NB),NPCOL)*NB,
     $               ICEIL(ICEIL(IHIZ-ILOZ+1,NB),NPROW)*NB *
     $               ICEIL(ICEIL(JW,NB),NPCOL)*NB,
     $               ICEIL(ICEIL(KBOT-JW,NB),NPROW)*NB *
     $               ICEIL(ICEIL(JW,NB),NPCOL)*NB )
*
*        Residual check workspace.
*
         LWK8 = 0
*
*        Optimal workspace.
*
         LWKOPT = MAX( LWK1, LWK2, LWK3+LWK4, LWK5, LWK6, LWK7, LWK8 )
         ILWKOPT = MAX( IWRK1, IWRK2 )
      END IF
*
*     Quick return in case of workspace query.
*
      WORK( 1 ) = DBLE( LWKOPT )
*
*     IWORK(1:NSEL) is used as the array SELECT for PDTRORD.
*
      IWORK( 1 ) = ILWKOPT + NSEL
      IF( LQUERY )
     $   RETURN
*
*     Nothing to do for an empty active block ...
      NS = 0
      ND = 0
      IF( KTOP.GT.KBOT )
     $   RETURN
*     ... nor for an empty deflation window.
*
      IF( NW.LT.1 )
     $   RETURN
*
*     Machine constants.
*
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
*
*     Setup deflation window.
*
      JW = MIN( NW, KBOT-KTOP+1 )
      KWTOP = KBOT - JW + 1
      IF( KWTOP.EQ.KTOP ) THEN
         S = ZERO
      ELSE
         CALL PDELGET( 'All', '1-Tree', S, H, KWTOP, KWTOP-1, DESCH )
      END IF
*
      IF( KBOT.EQ.KWTOP ) THEN
*
*        1-by-1 deflation window: not much to do.
*
         CALL PDELGET( 'All', '1-Tree', SR( KWTOP ), H, KWTOP, KWTOP,
     $        DESCH )
         SI( KWTOP ) = ZERO
         NS = 1
         ND = 0
         IF( ABS( S ).LE.MAX( SMLNUM, ULP*ABS( SR( KWTOP ) ) ) )
     $        THEN
            NS = 0
            ND = 1
            IF( KWTOP.GT.KTOP )
     $         CALL PDELSET( H, KWTOP, KWTOP-1 , DESCH, ZERO )
         END IF
         RETURN
      END IF
*
      IF( KWTOP.EQ.KTOP .AND. KBOT-KWTOP.EQ.1 ) THEN
*
*        2-by-2 deflation window: a little more to do.
*
         CALL PDELGET( 'All', '1-Tree', AA, H, KWTOP, KWTOP, DESCH )
         CALL PDELGET( 'All', '1-Tree', BB, H, KWTOP, KWTOP+1, DESCH )
         CALL PDELGET( 'All', '1-Tree', CC, H, KWTOP+1, KWTOP, DESCH )
         CALL PDELGET( 'All', '1-Tree', DD, H, KWTOP+1, KWTOP+1, DESCH )
         CALL DLANV2( AA, BB, CC, DD, SR(KWTOP), SI(KWTOP),
     $        SR(KWTOP+1), SI(KWTOP+1), CS, SN )
         NS = 0
         ND = 2
         IF( CC.EQ.ZERO ) THEN
            I = KWTOP
            IF( I+2.LE.N .AND. WANTT )
     $         CALL PDROT( N-I-1, H, I, I+2, DESCH, DESCH(M_), H, I+1,
     $              I+2, DESCH, DESCH(M_), CS, SN, WORK, LWORK, INFO )
            IF( I.GT.1 )
     $         CALL PDROT( I-1, H, 1, I, DESCH, 1, H, 1, I+1, DESCH, 1,
     $              CS, SN, WORK, LWORK, INFO )
            IF( WANTZ )
     $         CALL PDROT( IHIZ-ILOZ+1, Z, ILOZ, I, DESCZ, 1, Z, ILOZ,
     $              I+1, DESCZ, 1, CS, SN, WORK, LWORK, INFO )
            CALL PDELSET( H, I, I, DESCH, AA )
            CALL PDELSET( H, I, I+1, DESCH, BB )
            CALL PDELSET( H, I+1, I, DESCH, CC )
            CALL PDELSET( H, I+1, I+1, DESCH, DD )
         END IF
         WORK( 1 ) = DBLE( LWKOPT )
         RETURN
      END IF
*
*     Calculate new value for IROFFH in case deflation window
*     was adjusted.
*
      IROFFH = MOD( KWTOP - 1, NB )
*
*     Adjust number of rows and columns of T matrix descriptor
*     to prepare for call to PDBTRORD.
*
      DESCT( M_ ) = JW+IROFFH
      DESCT( N_ ) = JW+IROFFH
*
*     Convert to spike-triangular form.  (In case of a rare QR failure,
*     this routine continues to do aggressive early deflation using that
*     part of the deflation window that converged using INFQR here and
*     there to keep track.)
*
*     Copy the trailing submatrix to the working space.
*
      CALL PDLASET( 'All', IROFFH, JW+IROFFH, ZERO, ONE, T, 1, 1,
     $     DESCT )
      CALL PDLASET( 'All', JW, IROFFH, ZERO, ZERO, T, 1+IROFFH, 1,
     $     DESCT )
      CALL PDLACPY( 'All', 1, JW, H, KWTOP, KWTOP, DESCH, T, 1+IROFFH,
     $     1+IROFFH, DESCT )
      CALL PDLACPY( 'Upper', JW-1, JW-1, H, KWTOP+1, KWTOP, DESCH, T,
     $     1+IROFFH+1, 1+IROFFH, DESCT )
      IF( JW.GT.2 )
     $   CALL PDLASET( 'Lower', JW-2, JW-2, ZERO, ZERO, T, 1+IROFFH+2,
     $        1+IROFFH, DESCT )
      CALL PDLACPY( 'All', JW-1, 1, H, KWTOP+1, KWTOP+JW-1, DESCH, T,
     $     1+IROFFH+1, 1+IROFFH+JW-1, DESCT )
*
*     Initialize the working orthogonal matrix.
*
      CALL PDLASET( 'All', JW+IROFFH, JW+IROFFH, ZERO, ONE, V, 1, 1,
     $     DESCV )
*
*     Compute the Schur form of T.
*
      NPMIN = PILAENVX( ICTXT, 23, 'PDLAQR3', 'SV', JW, NB, NPROW,
     $     NPCOL )
      NMIN = PILAENVX( ICTXT, 12, 'PDLAQR3', 'SV', JW, 1, JW, LWORK )
      NMAX = ( N-1 ) / 3
      IF( MIN(NPROW, NPCOL).LE.NPMIN+1 .OR. RECLEVEL.GE.1 ) THEN
*
*        The AED window is large enough.
*        Compute the Schur decomposition with all processors.
*
         IF( JW+IROFFH.GT.NMIN .AND. JW+IROFFH.LE.NMAX
     $        .AND. RECLEVEL.LT.RECMAX ) THEN
            CALL PDLAQR0( .TRUE., .TRUE., JW+IROFFH, 1+IROFFH,
     $           JW+IROFFH, T, DESCT, SR( KWTOP-IROFFH ),
     $           SI( KWTOP-IROFFH ), 1+IROFFH, JW+IROFFH, V, DESCV,
     $           WORK, LWORK, IWORK(NSEL+1), LIWORK-NSEL, INFQR,
     $           RECLEVEL+1 )
         ELSE
            IF( DESCT(RSRC_).EQ.0 .AND. DESCT(CSRC_).EQ.0 ) THEN
               IF( JW+IROFFH.GT.DESCT( MB_ ) ) THEN
                  CALL PDLAQR1( .TRUE., .TRUE., JW+IROFFH, 1,
     $                 JW+IROFFH, T, DESCT, SR( KWTOP-IROFFH ),
     $                 SI( KWTOP-IROFFH ), 1, JW+IROFFH, V,
     $                 DESCV, WORK, LWORK, IWORK(NSEL+1), LIWORK-NSEL,
     $                 INFQR )
               ELSE
                  IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                     CALL DLAHQR( .TRUE., .TRUE., JW+IROFFH, 1+IROFFH,
     $                    JW+IROFFH, T, DESCT(LLD_),
     $                    SR( KWTOP-IROFFH ), SI( KWTOP-IROFFH ),
     $                    1+IROFFH, JW+IROFFH, V, DESCV(LLD_), INFQR )
                  ELSE
                     INFQR = 0
                  END IF
                  IF( NPROCS.GT.1 )
     $               CALL IGAMN2D( ICTXT, 'All', '1-Tree', 1, 1, INFQR,
     $                    1, -1, -1, -1, -1, -1 )
               END IF
            ELSEIF( JW+IROFFH.LE.DESCT( MB_ ) ) THEN
               IF( MYROW.EQ.DESCT(RSRC_) .AND. MYCOL.EQ.DESCT(CSRC_) )
     $              THEN
                  CALL DLAHQR( .TRUE., .TRUE., JW+IROFFH, 1+IROFFH,
     $                 JW+IROFFH, T, DESCT(LLD_),
     $                 SR( KWTOP-IROFFH ), SI( KWTOP-IROFFH ),
     $                 1+IROFFH, JW+IROFFH, V, DESCV(LLD_), INFQR )
               ELSE
                  INFQR = 0
               END IF
               IF( NPROCS.GT.1 )
     $         CALL IGAMN2D( ICTXT, 'All', '1-Tree', 1, 1, INFQR,
     $              1, -1, -1, -1, -1, -1 )
            ELSE
               TZROWS0 = NUMROC( JW+IROFFH, NB, MYROW, 0, NPROW )
               TZCOLS0 = NUMROC( JW+IROFFH, NB, MYCOL, 0, NPCOL )
               CALL DESCINIT( DESCTZ0, JW+IROFFH, JW+IROFFH, NB, NB, 0,
     $              0, ICTXT, MAX(1,TZROWS0), IERR0 )
               IPT0 = 1
               IPZ0 = IPT0 + MAX(1,TZROWS0)*TZCOLS0
               IPW0 = IPZ0 + MAX(1,TZROWS0)*TZCOLS0
               CALL PDLAMVE( 'All', JW+IROFFH, JW+IROFFH, T, 1, 1,
     $              DESCT, WORK(IPT0), 1, 1, DESCTZ0, WORK(IPW0) )
               CALL PDLASET( 'All', JW+IROFFH, JW+IROFFH, ZERO, ONE,
     $              WORK(IPZ0), 1, 1, DESCTZ0 )
               CALL PDLAQR1( .TRUE., .TRUE., JW+IROFFH, 1,
     $              JW+IROFFH, WORK(IPT0), DESCTZ0,
     $              SR( KWTOP-IROFFH ), SI( KWTOP-IROFFH ),
     $              1, JW+IROFFH, WORK(IPZ0),
     $              DESCTZ0, WORK(IPW0), LWORK-IPW0+1, IWORK(NSEL+1),
     $              LIWORK-NSEL, INFQR )
               CALL PDLAMVE( 'All', JW+IROFFH, JW+IROFFH, WORK(IPT0), 1,
     $              1, DESCTZ0, T, 1, 1, DESCT, WORK(IPW0) )
               CALL PDLAMVE( 'All', JW+IROFFH, JW+IROFFH, WORK(IPZ0), 1,
     $              1, DESCTZ0, V, 1, 1, DESCV, WORK(IPW0) )
            END IF
         END IF
      ELSE
*
*        The AED window is too small.
*        Redistribute the AED window to a subgrid
*        and do the computation on the subgrid.
*
         ICTXT_NEW = ICTXT
         DO 20 I = 0, NPMIN-1
            DO 10 J = 0, NPMIN-1
               PMAP( J+1+I*NPMIN ) = BLACS_PNUM( ICTXT, I, J )
 10         CONTINUE
 20      CONTINUE
         CALL BLACS_GRIDMAP( ICTXT_NEW, PMAP, NPMIN, NPMIN, NPMIN )
         CALL BLACS_GRIDINFO( ICTXT_NEW, NPMIN, NPMIN, MYROW_NEW,
     $        MYCOL_NEW )
         IF( MYROW.GE.NPMIN .OR. MYCOL.GE.NPMIN ) ICTXT_NEW = -1
         IF( ICTXT_NEW.GE.0 ) THEN
            TZROWS0 = NUMROC( JW, NB, MYROW_NEW, 0, NPMIN )
            TZCOLS0 = NUMROC( JW, NB, MYCOL_NEW, 0, NPMIN )
            CALL DESCINIT( DESCTZ0, JW, JW, NB, NB, 0,
     $           0, ICTXT_NEW, MAX(1,TZROWS0), IERR0 )
            IPT0 = 1
            IPZ0 = IPT0 + MAX(1,TZROWS0)*MAX(1,TZCOLS0)
            IPW0 = IPZ0 + MAX(1,TZROWS0)*MAX(1,TZCOLS0)
         ELSE
            IPT0 = 1
            IPZ0 = 2
            IPW0 = 3
            DESCTZ0( CTXT_ ) = -1
            INFQR = 0
         END IF
         CALL PDGEMR2D( JW, JW, T, 1+IROFFH, 1+IROFFH, DESCT,
     $        WORK(IPT0), 1, 1, DESCTZ0, ICTXT )
         IF( ICTXT_NEW.GE.0 ) THEN
            CALL PDLASET( 'All', JW, JW, ZERO, ONE, WORK(IPZ0), 1, 1,
     $           DESCTZ0 )
            NMIN = PILAENVX( ICTXT_NEW, 12, 'PDLAQR3', 'SV', JW, 1, JW,
     $           LWORK )
            IF( JW.GT.NMIN .AND. JW.LE.NMAX .AND. RECLEVEL.LT.1 ) THEN
               CALL PDLAQR0( .TRUE., .TRUE., JW, 1, JW, WORK(IPT0),
     $              DESCTZ0, SR( KWTOP ), SI( KWTOP ), 1, JW,
     $              WORK(IPZ0), DESCTZ0, WORK(IPW0), LWORK-IPW0+1,
     $              IWORK(NSEL+1), LIWORK-NSEL, INFQR,
     $              RECLEVEL+1 )
            ELSE
               CALL PDLAQR1( .TRUE., .TRUE., JW, 1, JW, WORK(IPT0),
     $              DESCTZ0, SR( KWTOP ), SI( KWTOP ), 1, JW,
     $              WORK(IPZ0), DESCTZ0, WORK(IPW0), LWORK-IPW0+1,
     $              IWORK(NSEL+1), LIWORK-NSEL, INFQR )
            END IF
         END IF
         CALL PDGEMR2D( JW, JW, WORK(IPT0), 1, 1, DESCTZ0, T, 1+IROFFH,
     $        1+IROFFH, DESCT, ICTXT )
         CALL PDGEMR2D( JW, JW, WORK(IPZ0), 1, 1, DESCTZ0, V, 1+IROFFH,
     $        1+IROFFH, DESCV, ICTXT )
         IF( ICTXT_NEW.GE.0 )
     $      CALL BLACS_GRIDEXIT( ICTXT_NEW )
         IF( MYROW+MYCOL.GT.0 ) THEN
            DO 40 J = 0, JW-1
               SR( KWTOP+J ) = ZERO
               SI( KWTOP+J ) = ZERO
 40         CONTINUE
         END IF
         CALL IGAMN2D( ICTXT, 'All', '1-Tree', 1, 1, INFQR, 1, -1, -1,
     $        -1, -1, -1 )
         CALL DGSUM2D( ICTXT, 'All', ' ', JW, 1, SR(KWTOP), JW, -1, -1 )
         CALL DGSUM2D( ICTXT, 'All', ' ', JW, 1, SI(KWTOP), JW, -1, -1 )
      END IF
*
*     Adjust INFQR for offset from block border in submatrices.
*
      IF( INFQR.NE.0 )
     $   INFQR = INFQR - IROFFH
*
*     PDTRORD needs a clean margin near the diagonal.
*
      DO 50 J = 1, JW - 3
         CALL PDELSET( T, J+2, J, DESCT, ZERO )
         CALL PDELSET( T, J+3, J, DESCT, ZERO )
 50   CONTINUE
      IF( JW.GT.2 )
     $   CALL PDELSET( T, JW, JW-2, DESCT, ZERO )
*
*     Check local residual for AED Schur decomposition.
*
      RESAED = 0.0D+00
*
*     Clean up the array SELECT for PDTRORD.
*
      DO 60 J = 1, NSEL
         IWORK( J ) = 0
 60   CONTINUE
*
*     Set local M counter to zero.
*
      MLOC = 0
*
*     Outer deflation detection loop (label 80).
*     In this loop a bunch of undeflatable eigenvalues
*     are moved simultaneously.
*
      DO 70 J = 1, IROFFH + INFQR
         IWORK( J ) = 1
 70   CONTINUE
*
      NS = JW
      ILST = INFQR + 1 + IROFFH
      IF( ILST.GT.1 ) THEN
         CALL PDELGET( 'All', '1-Tree', ELEM, T, ILST, ILST-1, DESCT )
         BULGE = ELEM.NE.ZERO
         IF( BULGE ) ILST = ILST+1
      END IF
*
 80   CONTINUE
      IF( ILST.LE.NS+IROFFH ) THEN
*
*        Find the top-left corner of the local window.
*
         LILST = MAX(ILST,NS+IROFFH-NB+1)
         IF( LILST.GT.1 ) THEN
            CALL PDELGET( 'All', '1-Tree', ELEM, T, LILST, LILST-1,
     $           DESCT )
            BULGE = ELEM.NE.ZERO
            IF( BULGE ) LILST = LILST+1
         END IF
*
*        Lock all eigenvalues outside the local window.
*
         DO 90 J = IROFFH+1, LILST-1
            IWORK( J ) = 1
 90      CONTINUE
         LILST0 = LILST
*
*        Inner deflation detection loop (label 100).
*        In this loop, the undeflatable eigenvalues are moved to the
*        top-left corner of the local window.
*
 100     CONTINUE
         IF( LILST.LE.NS+IROFFH ) THEN
            IF( NS.EQ.1 ) THEN
               BULGE = .FALSE.
            ELSE
               CALL PDELGET( 'All', '1-Tree', ELEM, T, NS+IROFFH,
     $              NS+IROFFH-1, DESCT )
               BULGE = ELEM.NE.ZERO
            END IF
*
*           Small spike tip test for deflation.
*
            IF( .NOT.BULGE ) THEN
*
*              Real eigenvalue.
*
               CALL PDELGET( 'All', '1-Tree', ELEM, T, NS+IROFFH,
     $              NS+IROFFH, DESCT )
               FOO = ABS( ELEM )
               IF( FOO.EQ.ZERO )
     $            FOO = ABS( S )
               CALL PDELGET( 'All', '1-Tree', ELEM, V, 1+IROFFH,
     $              NS+IROFFH, DESCV )
               IF( ABS( S*ELEM ).LE.MAX( SMLNUM, ULP*FOO ) ) THEN
*
*                 Deflatable.
*
                  NS = NS - 1
               ELSE
*
*                 Undeflatable: move it up out of the way.
*
                  IFST = NS
                  DO 110 J = LILST, JW+IROFFH
                     IWORK( J ) = 0
 110              CONTINUE
                  IWORK( IFST+IROFFH ) = 1
                  CALL PDTRORD( 'Vectors', IWORK, PAR, JW+IROFFH, T, 1,
     $                 1, DESCT, V, 1, 1, DESCV, WORK,
     $                 WORK(JW+IROFFH+1), MLOC,
     $                 WORK(2*(JW+IROFFH)+1), LWORK-2*(JW+IROFFH),
     $                 IWORK(NSEL+1), LIWORK-NSEL, INFO )
*
*                 Adjust the array SELECT explicitly so that it does not
*                 rely on the output of PDTRORD.
*
                  IWORK( IFST+IROFFH ) = 0
                  IWORK( LILST ) = 1
                  LILST = LILST + 1
*
*                 In case of a rare exchange failure, adjust the
*                 pointers ILST and LILST to the current place to avoid
*                 unexpected behaviors.
*
                  IF( INFO.NE.0 ) THEN
                     LILST = MAX(INFO, LILST)
                     ILST = MAX(INFO, ILST)
                  END IF
               END IF
            ELSE
*
*              Complex conjugate pair.
*
               CALL PDELGET( 'All', '1-Tree', ELEM1, T, NS+IROFFH,
     $              NS+IROFFH, DESCT )
               CALL PDELGET( 'All', '1-Tree', ELEM2, T, NS+IROFFH,
     $              NS+IROFFH-1, DESCT )
               CALL PDELGET( 'All', '1-Tree', ELEM3, T, NS+IROFFH-1,
     $              NS+IROFFH, DESCT )
               FOO = ABS( ELEM1 ) + SQRT( ABS( ELEM2 ) )*
     $              SQRT( ABS( ELEM3 ) )
               IF( FOO.EQ.ZERO )
     $            FOO = ABS( S )
               CALL PDELGET( 'All', '1-Tree', ELEM1, V, 1+IROFFH,
     $              NS+IROFFH, DESCV )
               CALL PDELGET( 'All', '1-Tree', ELEM2, V, 1+IROFFH,
     $              NS+IROFFH-1, DESCV )
               IF( MAX( ABS( S*ELEM1 ), ABS( S*ELEM2 ) ).LE.
     $              MAX( SMLNUM, ULP*FOO ) ) THEN
*
*                 Deflatable.
*
                  NS = NS - 2
               ELSE
*
*                 Undeflatable: move them up out of the way.
*
                  IFST = NS
                  DO 120 J = LILST, JW+IROFFH
                     IWORK( J ) = 0
 120              CONTINUE
                  IWORK( IFST+IROFFH ) = 1
                  IWORK( IFST+IROFFH-1 ) = 1
                  CALL PDTRORD( 'Vectors', IWORK, PAR, JW+IROFFH, T, 1,
     $                 1, DESCT, V, 1, 1, DESCV, WORK,
     $                 WORK(JW+IROFFH+1), MLOC,
     $                 WORK(2*(JW+IROFFH)+1), LWORK-2*(JW+IROFFH),
     $                 IWORK(NSEL+1), LIWORK-NSEL, INFO )
*
*                 Adjust the array SELECT explicitly so that it does not
*                 rely on the output of PDTRORD.
*
                  IWORK( IFST+IROFFH ) = 0
                  IWORK( IFST+IROFFH-1 ) = 0
                  IWORK( LILST ) = 1
                  IWORK( LILST+1 ) = 1
                  LILST = LILST + 2
*
*                 In case of a rare exchange failure, adjust the
*                 pointers ILST and LILST to the current place to avoid
*                 unexpected behaviors.
*
                  IF( INFO.NE.0 ) THEN
                     LILST = MAX(INFO, LILST)
                     ILST = MAX(INFO, ILST)
                  END IF
               END IF
            END IF
*
*           End of inner deflation detection loop.
*
            GO TO 100
         END IF
*
*        Unlock the eigenvalues outside the local window.
*        Then undeflatable eigenvalues are moved to the proper position.
*
         DO 130 J = ILST, LILST0-1
            IWORK( J ) = 0
 130     CONTINUE
         CALL PDTRORD( 'Vectors', IWORK, PAR, JW+IROFFH, T, 1, 1,
     $        DESCT, V, 1, 1, DESCV, WORK, WORK(JW+IROFFH+1),
     $        M, WORK(2*(JW+IROFFH)+1), LWORK-2*(JW+IROFFH),
     $        IWORK(NSEL+1), LIWORK-NSEL, INFO )
         ILST = M + 1
*
*        In case of a rare exchange failure, adjust the pointer ILST to
*        the current place to avoid unexpected behaviors.
*
         IF( INFO.NE.0 )
     $      ILST = MAX(INFO, ILST)
*
*        End of outer deflation detection loop.
*
         GO TO 80
      END IF

*
*     Post-reordering step: copy output eigenvalues to output.
*
      CALL DCOPY( JW, WORK(1+IROFFH), 1, SR( KWTOP ), 1 )
      CALL DCOPY( JW, WORK(JW+2*IROFFH+1), 1, SI( KWTOP ), 1 )
*
*     Check local residual for reordered AED Schur decomposition.
*
      RESAED = 0.0D+00
*
*     Return to Hessenberg form.
*
      IF( NS.EQ.0 )
     $   S = ZERO
*
      IF( NS.LT.JW .AND. SORTGRAD ) THEN
*
*        Sorting diagonal blocks of T improves accuracy for
*        graded matrices.  Bubble sort deals well with exchange
*        failures. Eigenvalues/shifts from T are also restored.
*
         ROUND = 0
         SORTED = .FALSE.
         I = NS + 1 + IROFFH
 140     CONTINUE
         IF( SORTED )
     $      GO TO 180
         SORTED = .TRUE.
         ROUND = ROUND + 1
*
         KEND = I - 1
         I = INFQR + 1 + IROFFH
         IF( I.EQ.NS+IROFFH ) THEN
            K = I + 1
         ELSE IF( SI( KWTOP-IROFFH + I-1 ).EQ.ZERO ) THEN
            K = I + 1
         ELSE
            K = I + 2
         END IF
 150     CONTINUE
         IF( K.LE.KEND ) THEN
            IF( K.EQ.I+1 ) THEN
               EVI = ABS( SR( KWTOP-IROFFH+I-1 ) )
            ELSE
               EVI = ABS( SR( KWTOP-IROFFH+I-1 ) ) +
     $              ABS( SI( KWTOP-IROFFH+I-1 ) )
            END IF
*
            IF( K.EQ.KEND ) THEN
               EVK = ABS( SR( KWTOP-IROFFH+K-1 ) )
            ELSEIF( SI( KWTOP-IROFFH+K-1 ).EQ.ZERO ) THEN
               EVK = ABS( SR( KWTOP-IROFFH+K-1 ) )
            ELSE
               EVK = ABS( SR( KWTOP-IROFFH+K-1 ) ) +
     $              ABS( SI( KWTOP-IROFFH+K-1 ) )
            END IF
*
            IF( EVI.GE.EVK ) THEN
               I = K
            ELSE
               MLOC = 0
               SORTED = .FALSE.
               IFST = I
               ILST = K
               DO 160 J = 1, I-1
                  IWORK( J ) = 1
                  MLOC = MLOC + 1
 160           CONTINUE
               IF( K.EQ.I+2 ) THEN
                  IWORK( I ) = 0
                  IWORK(I+1) = 0
               ELSE
                  IWORK( I ) = 0
               END IF
               IF( K.NE.KEND .AND. SI( KWTOP-IROFFH+K-1 ).NE.ZERO ) THEN
                  IWORK( K ) = 1
                  IWORK(K+1) = 1
                  MLOC = MLOC + 2
               ELSE
                  IWORK( K ) = 1
                  IF( K.LT.KEND ) IWORK(K+1) = 0
                  MLOC = MLOC + 1
               END IF
               DO 170 J = K+2, JW+IROFFH
                  IWORK( J ) = 0
 170           CONTINUE
               CALL PDTRORD( 'Vectors', IWORK, PAR, JW+IROFFH, T, 1, 1,
     $              DESCT, V, 1, 1, DESCV, WORK, WORK(JW+IROFFH+1), M,
     $              WORK(2*(JW+IROFFH)+1), LWORK-2*(JW+IROFFH),
     $              IWORK(NSEL+1), LIWORK-NSEL, IERR )
               CALL DCOPY( JW, WORK(1+IROFFH), 1, SR( KWTOP ), 1 )
               CALL DCOPY( JW, WORK(JW+2*IROFFH+1), 1, SI( KWTOP ), 1 )
               IF( IERR.EQ.0 ) THEN
                  I = ILST
               ELSE
                  I = K
               END IF
            END IF
            IF( I.EQ.KEND ) THEN
               K = I + 1
            ELSE IF( SI( KWTOP-IROFFH+I-1 ).EQ.ZERO ) THEN
               K = I + 1
            ELSE
               K = I + 2
            END IF
            GO TO 150
         END IF
         GO TO 140
 180     CONTINUE
      END IF
*
*     Restore number of rows and columns of T matrix descriptor.
*
      DESCT( M_ ) = NW+IROFFH
      DESCT( N_ ) = NH+IROFFH
*
      IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
*
*           Reflect spike back into lower triangle.
*
            RROWS = NUMROC( NS+IROFFH, NB, MYROW, DESCV(RSRC_), NPROW )
            RCOLS = NUMROC( 1, 1, MYCOL, DESCV(CSRC_), NPCOL )
            CALL DESCINIT( DESCR, NS+IROFFH, 1, NB, 1, DESCV(RSRC_),
     $           DESCV(CSRC_), ICTXT, MAX(1, RROWS), INFO )
            TAUROWS = NUMROC( 1, 1, MYCOL, DESCV(RSRC_), NPROW )
            TAUCOLS = NUMROC( JW+IROFFH, NB, MYCOL, DESCV(CSRC_),
     $           NPCOL )
            CALL DESCINIT( DESCTAU, 1, JW+IROFFH, 1, NB, DESCV(RSRC_),
     $           DESCV(CSRC_), ICTXT, MAX(1, TAUROWS), INFO )
*
            IR = 1
            ITAU = IR + DESCR( LLD_ ) * RCOLS
            IPW  = ITAU + DESCTAU( LLD_ ) * TAUCOLS
*
            CALL PDLASET( 'All', NS+IROFFH, 1, ZERO, ZERO, WORK(ITAU),
     $           1, 1, DESCTAU )
*
            CALL PDCOPY( NS, V, 1+IROFFH, 1+IROFFH, DESCV, DESCV(M_),
     $           WORK(IR), 1+IROFFH, 1, DESCR, 1 )
            CALL PDLARFG( NS, BETA, 1+IROFFH, 1, WORK(IR), 2+IROFFH, 1,
     $           DESCR, 1, WORK(ITAU+IROFFH) )
            CALL PDELSET( WORK(IR), 1+IROFFH, 1, DESCR, ONE )
*
            CALL PDLASET( 'Lower', JW-2, JW-2, ZERO, ZERO, T, 3+IROFFH,
     $           1+IROFFH, DESCT )
*
            CALL PDLARF( 'Left', NS, JW, WORK(IR), 1+IROFFH, 1, DESCR,
     $           1, WORK(ITAU+IROFFH), T, 1+IROFFH, 1+IROFFH,
     $           DESCT, WORK( IPW ) )
            CALL PDLARF( 'Right', NS, NS, WORK(IR), 1+IROFFH, 1, DESCR,
     $           1, WORK(ITAU+IROFFH), T, 1+IROFFH, 1+IROFFH,
     $           DESCT, WORK( IPW ) )
            CALL PDLARF( 'Right', JW, NS, WORK(IR), 1+IROFFH, 1, DESCR,
     $           1, WORK(ITAU+IROFFH), V, 1+IROFFH, 1+IROFFH,
     $           DESCV, WORK( IPW ) )
*
            ITAU = 1
            IPW = ITAU + DESCTAU( LLD_ ) * TAUCOLS
            CALL PDGEHRD( JW+IROFFH, 1+IROFFH, NS+IROFFH, T, 1, 1,
     $           DESCT, WORK(ITAU), WORK( IPW ), LWORK-IPW+1, INFO )
         END IF
*
*        Copy updated reduced window into place.
*
         IF( KWTOP.GT.1 ) THEN
            CALL PDELGET( 'All', '1-Tree', ELEM, V, 1+IROFFH,
     $           1+IROFFH, DESCV )
            CALL PDELSET( H, KWTOP, KWTOP-1, DESCH, S*ELEM )
         END IF
         CALL PDLACPY( 'Upper', JW-1, JW-1, T, 1+IROFFH+1, 1+IROFFH,
     $        DESCT, H, KWTOP+1, KWTOP, DESCH )
         CALL PDLACPY( 'All', 1, JW, T, 1+IROFFH, 1+IROFFH, DESCT, H,
     $        KWTOP, KWTOP, DESCH )
         CALL PDLACPY( 'All', JW-1, 1, T, 1+IROFFH+1, 1+IROFFH+JW-1,
     $        DESCT, H, KWTOP+1, KWTOP+JW-1, DESCH )
*
*        Accumulate orthogonal matrix in order to update
*        H and Z, if requested.
*
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
            CALL PDORMHR( 'Right', 'No', JW+IROFFH, NS+IROFFH, 1+IROFFH,
     $           NS+IROFFH, T, 1, 1, DESCT, WORK(ITAU), V, 1,
     $           1, DESCV, WORK( IPW ), LWORK-IPW+1, INFO )
         END IF
*
*        Update vertical slab in H.
*
         IF( WANTT ) THEN
            LTOP = 1
         ELSE
            LTOP = KTOP
         END IF
         KLN = MAX( 0, KWTOP-LTOP )
         IROFFHH = MOD( LTOP-1, NB )
         ICOFFHH = MOD( KWTOP-1, NB )
         HHRSRC = INDXG2P( LTOP, NB, MYROW, DESCH(RSRC_), NPROW )
         HHCSRC = INDXG2P( KWTOP, NB, MYCOL, DESCH(CSRC_), NPCOL )
         HHROWS = NUMROC( KLN+IROFFHH, NB, MYROW, HHRSRC, NPROW )
         HHCOLS = NUMROC( JW+ICOFFHH, NB, MYCOL, HHCSRC, NPCOL )
         CALL DESCINIT( DESCHH, KLN+IROFFHH, JW+ICOFFHH, NB, NB,
     $        HHRSRC, HHCSRC, ICTXT, MAX(1, HHROWS), IERR )
         CALL PDGEMM( 'No', 'No', KLN, JW, JW, ONE, H, LTOP,
     $        KWTOP, DESCH, V, 1+IROFFH, 1+IROFFH, DESCV, ZERO,
     $        WORK, 1+IROFFHH, 1+ICOFFHH, DESCHH )
         CALL PDLACPY( 'All', KLN, JW, WORK, 1+IROFFHH, 1+ICOFFHH,
     $        DESCHH, H, LTOP, KWTOP, DESCH )
*
*        Update horizontal slab in H.
*
         IF( WANTT ) THEN
            KLN = N-KBOT
            IROFFHH = MOD( KWTOP-1, NB )
            ICOFFHH = MOD( KBOT, NB )
            HHRSRC = INDXG2P( KWTOP, NB, MYROW, DESCH(RSRC_), NPROW )
            HHCSRC = INDXG2P( KBOT+1, NB, MYCOL, DESCH(CSRC_), NPCOL )
            HHROWS = NUMROC( JW+IROFFHH, NB, MYROW, HHRSRC, NPROW )
            HHCOLS = NUMROC( KLN+ICOFFHH, NB, MYCOL, HHCSRC, NPCOL )
            CALL DESCINIT( DESCHH, JW+IROFFHH, KLN+ICOFFHH, NB, NB,
     $           HHRSRC, HHCSRC, ICTXT, MAX(1, HHROWS), IERR )
            CALL PDGEMM( 'Tr', 'No', JW, KLN, JW, ONE, V,
     $           1+IROFFH, 1+IROFFH, DESCV, H, KWTOP, KBOT+1,
     $           DESCH, ZERO, WORK, 1+IROFFHH, 1+ICOFFHH, DESCHH )
            CALL PDLACPY( 'All', JW, KLN, WORK, 1+IROFFHH, 1+ICOFFHH,
     $           DESCHH, H, KWTOP, KBOT+1, DESCH )
         END IF
*
*        Update vertical slab in Z.
*
         IF( WANTZ ) THEN
            KLN = IHIZ-ILOZ+1
            IROFFZZ = MOD( ILOZ-1, NB )
            ICOFFZZ = MOD( KWTOP-1, NB )
            ZZRSRC = INDXG2P( ILOZ, NB, MYROW, DESCZ(RSRC_), NPROW )
            ZZCSRC = INDXG2P( KWTOP, NB, MYCOL, DESCZ(CSRC_), NPCOL )
            ZZROWS = NUMROC( KLN+IROFFZZ, NB, MYROW, ZZRSRC, NPROW )
            ZZCOLS = NUMROC( JW+ICOFFZZ, NB, MYCOL, ZZCSRC, NPCOL )
            CALL DESCINIT( DESCZZ, KLN+IROFFZZ, JW+ICOFFZZ, NB, NB,
     $           ZZRSRC, ZZCSRC, ICTXT, MAX(1, ZZROWS), IERR )
            CALL PDGEMM( 'No', 'No', KLN, JW, JW, ONE, Z, ILOZ,
     $           KWTOP, DESCZ, V, 1+IROFFH, 1+IROFFH, DESCV,
     $           ZERO, WORK, 1+IROFFZZ, 1+ICOFFZZ, DESCZZ )
            CALL PDLACPY( 'All', KLN, JW, WORK, 1+IROFFZZ, 1+ICOFFZZ,
     $           DESCZZ, Z, ILOZ, KWTOP, DESCZ )
         END IF
      END IF
*
*     Return the number of deflations (ND) and the number of shifts (NS).
*     (Subtracting INFQR from the spike length takes care of the case of
*     a rare QR failure while calculating eigenvalues of the deflation
*     window.)
*
      ND = JW - NS
      NS = NS - INFQR
*
*     Return optimal workspace.
*
      WORK( 1 ) = DBLE( LWKOPT )
      IWORK( 1 ) = ILWKOPT + NSEL
*
*     End of PDLAQR3
*
      END
