      SUBROUTINE PSLAQR4( WANTT, WANTZ, N, ILO, IHI, A, DESCA, WR, WI,
     $                    ILOZ, IHIZ, Z, DESCZ, T, LDT, V, LDV, WORK,
     $                    LWORK, INFO )
*
*     Contribution from the Department of Computing Science and HPC2N,
*     Umea University, Sweden
*
*  -- ScaLAPACK auxiliary routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      LOGICAL            WANTT, WANTZ
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDT, LDV, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * )
      REAL               A( * ), T( LDT, * ), V( LDV, * ), WI( * ),
     $                   WORK( * ), WR( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLAQR4 is an auxiliary routine used to find the Schur decomposition
*  and or eigenvalues of a matrix already in Hessenberg form from cols
*  ILO to IHI.  This routine requires that the active block is small
*  enough, i.e. IHI-ILO+1 .LE. LDT, so that it can be solved by LAPACK.
*  Normally, it is called by PSLAQR1.  All the inputs are assumed to be
*  valid without checking.
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
*          = .TRUE. : the full Schur form T is required;
*          = .FALSE.: only eigenvalues are required.
*
*  WANTZ   (global input) LOGICAL
*          = .TRUE. : the matrix of Schur vectors Z is required;
*          = .FALSE.: Schur vectors are not required.
*
*  N       (global input) INTEGER
*          The order of the Hessenberg matrix A (and Z if WANTZ).
*          N >= 0.
*
*  ILO     (global input) INTEGER
*  IHI     (global input) INTEGER
*          It is assumed that A is already upper quasi-triangular in
*          rows and columns IHI+1:N, and that A(ILO,ILO-1) = 0 (unless
*          ILO = 1). PSLAQR4 works primarily with the Hessenberg
*          submatrix in rows and columns ILO to IHI, but applies
*          transformations to all of H if WANTT is .TRUE..
*          1 <= ILO <= max(1,IHI); IHI <= N.
*
*  A       (global input/output) REAL             array, dimension
*          (DESCA(LLD_),*)
*          On entry, the upper Hessenberg matrix A.
*          On exit, if WANTT is .TRUE., A is upper quasi-triangular in
*          rows and columns ILO:IHI, with any 2-by-2 or larger diagonal
*          blocks not yet in standard form. If WANTT is .FALSE., the
*          contents of A are unspecified on exit.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  WR      (global replicated output) REAL             array,
*                                                         dimension (N)
*  WI      (global replicated output) REAL             array,
*                                                         dimension (N)
*          The real and imaginary parts, respectively, of the computed
*          eigenvalues ILO to IHI are stored in the corresponding
*          elements of WR and WI. If two eigenvalues are computed as a
*          complex conjugate pair, they are stored in consecutive
*          elements of WR and WI, say the i-th and (i+1)th, with
*          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
*          eigenvalues are stored in the same order as on the diagonal
*          of the Schur form returned in A.  A may be returned with
*          larger diagonal blocks until the next release.
*
*  ILOZ    (global input) INTEGER
*  IHIZ    (global input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if WANTZ is .TRUE..
*          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
*
*  Z       (global input/output) REAL             array.
*          If WANTZ is .TRUE., on entry Z must contain the current
*          matrix Z of transformations accumulated by PDHSEQR, and on
*          exit Z has been updated; transformations are applied only to
*          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
*          If WANTZ is .FALSE., Z is not referenced.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*  T       (local workspace) REAL             array, dimension LDT*NW.
*
*  LDT     (local input) INTEGER
*          The leading dimension of the array T.
*          LDT >= IHI-ILO+1.
*
*  V       (local workspace) REAL             array, dimension LDV*NW.
*
*  LDV     (local input) INTEGER
*          The leading dimension of the array V.
*          LDV >= IHI-ILO+1.
*
*  WORK    (local workspace) REAL             array, dimension LWORK.
*
*  LWORK   (local input) INTEGER
*          The dimension of the work array WORK.
*          LWORK >= IHI-ILO+1.
*          WORK(LWORK) is a local array and LWORK is assumed big enough.
*          Typically LWORK >= 4*LDS*LDS if this routine is called by
*          PSLAQR1. (LDS = 385, see PSLAQR1)
*
*  INFO    (global output) INTEGER
*          < 0: parameter number -INFO incorrect or inconsistent;
*          = 0: successful exit;
*          > 0: PSLAQR4 failed to compute all the eigenvalues ILO to IHI
*               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
*               elements i+1:ihi of WR and WI contain those eigenvalues
*               which have been successfully computed.
*
*  ================================================================
*  Implemented by
*        Meiyue Shao, Department of Computing Science and HPC2N,
*        Umea University, Sweden
*
*  ================================================================
*  References:
*        B. Kagstrom, D. Kressner, and M. Shao,
*        On Aggressive Early Deflation in Parallel Variants of the QR
*        Algorithm.
*        Para 2010, to appear.
*
*  ================================================================
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0, ONE = 1.0 )
*     ..
*     .. Local Scalars ..
      INTEGER            CONTXT, HBL, I, I1, I2, IAFIRST, ICOL, ICOL1,
     $                   ICOL2, II, IROW, IROW1, IROW2, ITMP1, ITMP2,
     $                   IERR, J, JAFIRST, JJ, K, L, LDA, LDZ, LLDTMP,
     $                   MYCOL, MYROW, NODE, NPCOL, NPROW, NH, NMIN, NZ,
     $                   HSTEP, VSTEP, KKROW, KKCOL, KLN, LTOP, LEFT,
     $                   RIGHT, UP, DOWN, D1, D2
*     ..
*     .. Local Arrays ..
      INTEGER            DESCT( 9 ), DESCV( 9 ), DESCWH( 9 ),
     $                   DESCWV( 9 )
*     ..
*     .. External Functions ..
      INTEGER            NUMROC, ILAENV
      EXTERNAL           NUMROC, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, SLASET,
     $                   SLAHQR, SLAQR4, DESCINIT, PSGEMM, PSGEMR2D,
     $                   SGEMM, SLAMOV, SGESD2D, SGERV2D,
     $                   SGEBS2D, SGEBR2D, IGEBS2D, IGEBR2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
      IF( N.EQ.0 .OR. NH.EQ.0 )
     $   RETURN
*
*     NODE (IAFIRST,JAFIRST) OWNS A(1,1)
*
      HBL = DESCA( MB_ )
      CONTXT = DESCA( CTXT_ )
      LDA = DESCA( LLD_ )
      IAFIRST = DESCA( RSRC_ )
      JAFIRST = DESCA( CSRC_ )
      LDZ = DESCZ( LLD_ )
      CALL BLACS_GRIDINFO( CONTXT, NPROW, NPCOL, MYROW, MYCOL )
      NODE = MYROW*NPCOL + MYCOL
      LEFT = MOD( MYCOL+NPCOL-1, NPCOL )
      RIGHT = MOD( MYCOL+1, NPCOL )
      UP = MOD( MYROW+NPROW-1, NPROW )
      DOWN = MOD( MYROW+1, NPROW )
*
*     I1 and I2 are the indices of the first row and last column of A
*     to which transformations must be applied.
*
      I = IHI
      L = ILO
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
         LTOP = 1
      ELSE
         I1 = L
         I2 = I
         LTOP = L
      END IF
*
*     Copy the diagonal block to local and call LAPACK.
*
      CALL INFOG2L( ILO, ILO, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $     IROW, ICOL, II, JJ )
      IF ( MYROW .EQ. II ) THEN
         CALL DESCINIT( DESCT, NH, NH, NH, NH, II, JJ, CONTXT,
     $        LDT, IERR )
         CALL DESCINIT( DESCV, NH, NH, NH, NH, II, JJ, CONTXT,
     $        LDV, IERR )
      ELSE
         CALL DESCINIT( DESCT, NH, NH, NH, NH, II, JJ, CONTXT,
     $        1, IERR )
         CALL DESCINIT( DESCV, NH, NH, NH, NH, II, JJ, CONTXT,
     $        1, IERR )
      END IF
      CALL PSGEMR2D( NH, NH, A, ILO, ILO, DESCA, T, 1, 1, DESCT,
     $     CONTXT )
      IF ( MYROW .EQ. II .AND. MYCOL .EQ. JJ ) THEN
         CALL SLASET( 'All', NH, NH, ZERO, ONE, V, LDV )
         NMIN = ILAENV( 12, 'SLAQR3', 'SV', NH, 1, NH, LWORK )
         IF( NH .GT. NMIN ) THEN
            CALL SLAQR4( .TRUE., .TRUE., NH, 1, NH, T, LDT, WR( ILO ),
     $           WI( ILO ), 1, NH, V, LDV, WORK, LWORK, INFO )
*           Clean up the scratch used by SLAQR4.
            CALL SLASET( 'L', NH-2, NH-2, ZERO, ZERO, T( 3, 1 ), LDT )
         ELSE
            CALL SLAHQR( .TRUE., .TRUE., NH, 1, NH, T, LDT, WR( ILO ),
     $           WI( ILO ), 1, NH, V, LDV, INFO )
         END IF
         CALL SGEBS2D( CONTXT, 'All', ' ', NH, NH, V, LDV )
         CALL IGEBS2D( CONTXT, 'All', ' ', 1, 1, INFO, 1 )
      ELSE
         CALL SGEBR2D( CONTXT, 'All', ' ', NH, NH, V, LDV, II, JJ )
         CALL IGEBR2D( CONTXT, 'All', ' ', 1, 1, INFO, 1, II, JJ )
      END IF
      IF( INFO .NE. 0 ) INFO = INFO+ILO-1
*
*     Copy the local matrix back to the diagonal block.
*
      CALL PSGEMR2D( NH, NH, T, 1, 1, DESCT, A, ILO, ILO, DESCA,
     $     CONTXT )
*
*     Update T and Z.
*
      IF( MOD( ILO-1, HBL )+NH .LE. HBL ) THEN
*
*        Simplest case: the diagonal block is located on one processor.
*        Call SGEMM directly to perform the update.
*
         HSTEP = LWORK / NH
         VSTEP = HSTEP
*
         IF( WANTT ) THEN
*
*           Update horizontal slab in A.
*
            CALL INFOG2L( ILO, I+1, DESCA, NPROW, NPCOL, MYROW,
     $           MYCOL, IROW, ICOL, II, JJ )
            IF( MYROW .EQ. II ) THEN
               ICOL1 = NUMROC( N, HBL, MYCOL, JAFIRST, NPCOL )
               DO 10 KKCOL = ICOL, ICOL1, HSTEP
                  KLN = MIN( HSTEP, ICOL1-KKCOL+1 )
                  CALL SGEMM( 'T', 'N', NH, KLN, NH, ONE, V,
     $                 LDV, A( IROW+(KKCOL-1)*LDA ), LDA, ZERO, WORK,
     $                 NH )
                  CALL SLAMOV( 'A', NH, KLN, WORK, NH,
     $                 A( IROW+(KKCOL-1)*LDA ), LDA )
   10          CONTINUE
            END IF
*
*           Update vertical slab in A.
*
            CALL INFOG2L( LTOP, ILO, DESCA, NPROW, NPCOL, MYROW,
     $           MYCOL, IROW, ICOL, II, JJ )
            IF( MYCOL .EQ. JJ ) THEN
               CALL INFOG2L( ILO-1, ILO, DESCA, NPROW, NPCOL,
     $              MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
               IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
               DO 20 KKROW = IROW, IROW1, VSTEP
                  KLN = MIN( VSTEP, IROW1-KKROW+1 )
                  CALL SGEMM( 'N', 'N', KLN, NH, NH, ONE,
     $                 A( KKROW+(ICOL-1)*LDA ), LDA, V, LDV, ZERO,
     $                 WORK, KLN )
                  CALL SLAMOV( 'A', KLN, NH, WORK, KLN,
     $                 A( KKROW+(ICOL-1)*LDA ), LDA )
   20          CONTINUE
            END IF
         END IF
*
*        Update vertical slab in Z.
*
         IF( WANTZ ) THEN
            CALL INFOG2L( ILOZ, ILO, DESCZ, NPROW, NPCOL, MYROW,
     $           MYCOL, IROW, ICOL, II, JJ )
            IF( MYCOL .EQ. JJ ) THEN
               CALL INFOG2L( IHIZ, ILO, DESCZ, NPROW, NPCOL,
     $              MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
               IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
               DO 30 KKROW = IROW, IROW1, VSTEP
                  KLN = MIN( VSTEP, IROW1-KKROW+1 )
                  CALL SGEMM( 'N', 'N', KLN, NH, NH, ONE,
     $                 Z( KKROW+(ICOL-1)*LDZ ), LDZ, V, LDV, ZERO,
     $                 WORK, KLN )
                  CALL SLAMOV( 'A', KLN, NH, WORK, KLN,
     $                 Z( KKROW+(ICOL-1)*LDZ ), LDZ )
   30          CONTINUE
            END IF
         END IF
*
      ELSE IF( MOD( ILO-1, HBL )+NH .LE. 2*HBL ) THEN
*
*        More complicated case: the diagonal block lay on a 2x2
*        processor mesh.
*        Call SGEMM locally and communicate by pair.
*
         D1 = HBL - MOD( ILO-1, HBL )
         D2 = NH - D1
         HSTEP = LWORK / NH
         VSTEP = HSTEP
*
         IF( WANTT ) THEN
*
*           Update horizontal slab in A.
*
            CALL INFOG2L( ILO, I+1, DESCA, NPROW, NPCOL, MYROW,
     $           MYCOL, IROW, ICOL, II, JJ )
            IF( MYROW .EQ. UP ) THEN
               IF( MYROW .EQ. II ) THEN
                  ICOL1 = NUMROC( N, HBL, MYCOL, JAFIRST, NPCOL )
                  DO 40 KKCOL = ICOL, ICOL1, HSTEP
                     KLN = MIN( HSTEP, ICOL1-KKCOL+1 )
                     CALL SGEMM( 'T', 'N', NH, KLN, NH, ONE, V,
     $                    NH, A( IROW+(KKCOL-1)*LDA ), LDA, ZERO,
     $                    WORK, NH )
                     CALL SLAMOV( 'A', NH, KLN, WORK, NH,
     $                    A( IROW+(KKCOL-1)*LDA ), LDA )
   40             CONTINUE
               END IF
            ELSE
               IF( MYROW .EQ. II ) THEN
                  ICOL1 = NUMROC( N, HBL, MYCOL, JAFIRST, NPCOL )
                  DO 50 KKCOL = ICOL, ICOL1, HSTEP
                     KLN = MIN( HSTEP, ICOL1-KKCOL+1 )
                     CALL SGEMM( 'T', 'N', D2, KLN, D1, ONE,
     $                    V( 1, D1+1 ), LDV, A( IROW+(KKCOL-1)*LDA ),
     $                    LDA, ZERO, WORK( D1+1 ), NH )
                     CALL SGESD2D( CONTXT, D2, KLN, WORK( D1+1 ),
     $                    NH, DOWN, MYCOL )
                     CALL SGERV2D( CONTXT, D1, KLN, WORK, NH, DOWN,
     $                    MYCOL )
                     CALL SGEMM( 'T', 'N', D1, KLN, D1, ONE,
     $                    V, LDV, A( IROW+(KKCOL-1)*LDA ), LDA, ONE,
     $                    WORK, NH )
                     CALL SLAMOV( 'A', D1, KLN, WORK, NH,
     $                    A( IROW+(KKCOL-1)*LDA ), LDA )
   50             CONTINUE
               ELSE IF( UP .EQ. II ) THEN
                  ICOL1 = NUMROC( N, HBL, MYCOL, JAFIRST, NPCOL )
                  DO 60 KKCOL = ICOL, ICOL1, HSTEP
                     KLN = MIN( HSTEP, ICOL1-KKCOL+1 )
                     CALL SGEMM( 'T', 'N', D1, KLN, D2, ONE,
     $                    V( D1+1, 1 ), LDV, A( IROW+(KKCOL-1)*LDA ),
     $                    LDA, ZERO, WORK, NH )
                     CALL SGESD2D( CONTXT, D1, KLN, WORK, NH, UP,
     $                    MYCOL )
                     CALL SGERV2D( CONTXT, D2, KLN, WORK( D1+1 ),
     $                    NH, UP, MYCOL )
                     CALL SGEMM( 'T', 'N', D2, KLN, D2, ONE,
     $                    V( D1+1, D1+1 ), LDV,
     $                    A( IROW+(KKCOL-1)*LDA ), LDA, ONE,
     $                    WORK( D1+1 ), NH )
                     CALL SLAMOV( 'A', D2, KLN, WORK( D1+1 ), NH,
     $                    A( IROW+(KKCOL-1)*LDA ), LDA )
   60             CONTINUE
               END IF
            END IF
*
*           Update vertical slab in A.
*
            CALL INFOG2L( LTOP, ILO, DESCA, NPROW, NPCOL, MYROW,
     $           MYCOL, IROW, ICOL, II, JJ )
            IF( MYCOL .EQ. LEFT ) THEN
               IF( MYCOL .EQ. JJ ) THEN
                  CALL INFOG2L( ILO-1, ILO, DESCA, NPROW, NPCOL,
     $                 MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 70 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL SGEMM( 'N', 'N', KLN, NH, NH, ONE,
     $                    A( KKROW+(ICOL-1)*LDA ), LDA, V, LDV,
     $                    ZERO, WORK, KLN )
                     CALL SLAMOV( 'A', KLN, NH, WORK, KLN,
     $                    A( KKROW+(ICOL-1)*LDA ), LDA )
   70             CONTINUE
               END IF
            ELSE
               IF( MYCOL .EQ. JJ ) THEN
                  CALL INFOG2L( ILO-1, ILO, DESCA, NPROW, NPCOL,
     $                 MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 80 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL SGEMM( 'N', 'N', KLN, D2, D1, ONE,
     $                    A( KKROW+(ICOL-1)*LDA ), LDA, V( 1, D1+1 ),
     $                    LDV, ZERO, WORK( 1+D1*KLN ), KLN )
                     CALL SGESD2D( CONTXT, KLN, D2, WORK( 1+D1*KLN ),
     $                    KLN, MYROW, RIGHT )
                     CALL SGERV2D( CONTXT, KLN, D1, WORK, KLN, MYROW,
     $                    RIGHT )
                     CALL SGEMM( 'N', 'N', KLN, D1, D1, ONE,
     $                    A( KKROW+(ICOL-1)*LDA ), LDA, V, LDV, ONE,
     $                    WORK, KLN )
                     CALL SLAMOV( 'A', KLN, D1, WORK, KLN,
     $                    A( KKROW+(ICOL-1)*LDA ), LDA )
   80             CONTINUE
               ELSE IF ( LEFT .EQ. JJ ) THEN
                  CALL INFOG2L( ILO-1, ILO, DESCA, NPROW, NPCOL,
     $                 MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 90 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL SGEMM( 'N', 'N', KLN, D1, D2, ONE,
     $                    A( KKROW+(ICOL-1)*LDA ), LDA, V( D1+1, 1 ),
     $                    LDV, ZERO, WORK, KLN )
                     CALL SGESD2D( CONTXT, KLN, D1, WORK, KLN, MYROW,
     $                    LEFT )
                     CALL SGERV2D( CONTXT, KLN, D2, WORK( 1+D1*KLN ),
     $                    KLN, MYROW, LEFT )
                     CALL SGEMM( 'N', 'N', KLN, D2, D2, ONE,
     $                    A( KKROW+(ICOL-1)*LDA ), LDA, V( D1+1, D1+1 ),
     $                    LDV, ONE, WORK( 1+D1*KLN ), KLN )
                     CALL SLAMOV( 'A', KLN, D2, WORK( 1+D1*KLN ), KLN,
     $                    A( KKROW+(ICOL-1)*LDA ), LDA )
   90             CONTINUE
               END IF
            END IF
         END IF
*
*        Update vertical slab in Z.
*
         IF( WANTZ ) THEN
            CALL INFOG2L( ILOZ, ILO, DESCZ, NPROW, NPCOL, MYROW,
     $           MYCOL, IROW, ICOL, II, JJ )
            IF( MYCOL .EQ. LEFT ) THEN
               IF( MYCOL .EQ. JJ ) THEN
                  CALL INFOG2L( IHIZ, ILO, DESCZ, NPROW, NPCOL,
     $                 MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 100 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL SGEMM( 'N', 'N', KLN, NH, NH, ONE,
     $                    Z( KKROW+(ICOL-1)*LDZ ), LDZ, V, LDV, ZERO,
     $                    WORK, KLN )
                     CALL SLAMOV( 'A', KLN, NH, WORK, KLN,
     $                    Z( KKROW+(ICOL-1)*LDZ ), LDZ )
  100             CONTINUE
               END IF
            ELSE
               IF( MYCOL .EQ. JJ ) THEN
                  CALL INFOG2L( IHIZ, ILO, DESCZ, NPROW, NPCOL,
     $                 MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 110 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL SGEMM( 'N', 'N', KLN, D2, D1, ONE,
     $                    Z( KKROW+(ICOL-1)*LDZ ), LDZ, V( 1, D1+1 ),
     $                    LDV, ZERO, WORK( 1+D1*KLN ), KLN )
                     CALL SGESD2D( CONTXT, KLN, D2, WORK( 1+D1*KLN ),
     $                    KLN, MYROW, RIGHT )
                     CALL SGERV2D( CONTXT, KLN, D1, WORK, KLN, MYROW,
     $                    RIGHT )
                     CALL SGEMM( 'N', 'N', KLN, D1, D1, ONE,
     $                    Z( KKROW+(ICOL-1)*LDZ ), LDZ, V, LDV, ONE,
     $                    WORK, KLN )
                     CALL SLAMOV( 'A', KLN, D1, WORK, KLN,
     $                    Z( KKROW+(ICOL-1)*LDZ ), LDZ )
  110             CONTINUE
               ELSE IF( LEFT .EQ. JJ ) THEN
                  CALL INFOG2L( IHIZ, ILO, DESCZ, NPROW, NPCOL,
     $                 MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 120 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL SGEMM( 'N', 'N', KLN, D1, D2, ONE,
     $                    Z( KKROW+(ICOL-1)*LDZ ), LDZ, V( D1+1, 1 ),
     $                    LDV, ZERO, WORK, KLN )
                     CALL SGESD2D( CONTXT, KLN, D1, WORK, KLN, MYROW,
     $                    LEFT )
                     CALL SGERV2D( CONTXT, KLN, D2, WORK( 1+D1*KLN ),
     $                    KLN, MYROW, LEFT )
                     CALL SGEMM( 'N', 'N', KLN, D2, D2, ONE,
     $                    Z( KKROW+(ICOL-1)*LDZ ), LDZ,
     $                    V( D1+1, D1+1 ), LDV, ONE, WORK( 1+D1*KLN ),
     $                    KLN )
                     CALL SLAMOV( 'A', KLN, D2, WORK( 1+D1*KLN ),
     $                    KLN, Z( KKROW+(ICOL-1)*LDZ ), LDZ )
  120             CONTINUE
               END IF
            END IF
         END IF
*
      ELSE
*
*        Most complicated case: the diagonal block lay across the border
*        of the processor mesh.
*        Treat V as a distributed matrix and call PSGEMM.
*
         HSTEP = LWORK / NH * NPCOL
         VSTEP = LWORK / NH * NPROW
         LLDTMP = NUMROC( NH, NH, MYROW, 0, NPROW )
         LLDTMP = MAX( 1, LLDTMP )
         CALL DESCINIT( DESCV, NH, NH, NH, NH, 0, 0, CONTXT,
     $        LLDTMP, IERR )
         CALL DESCINIT( DESCWH, NH, HSTEP, NH, LWORK / NH, 0, 0,
     $        CONTXT, LLDTMP, IERR )
*
         IF( WANTT ) THEN
*
*           Update horizontal slab in A.
*
            DO 130 KKCOL = I+1, N, HSTEP
               KLN = MIN( HSTEP, N-KKCOL+1 )
               CALL PSGEMM( 'T', 'N', NH, KLN, NH, ONE, V, 1, 1,
     $              DESCV, A, ILO, KKCOL, DESCA, ZERO, WORK, 1, 1,
     $              DESCWH )
               CALL PSGEMR2D( NH, KLN, WORK, 1, 1, DESCWH, A,
     $              ILO, KKCOL, DESCA, CONTXT )
  130       CONTINUE
*
*           Update vertical slab in A.
*
            DO 140 KKROW = LTOP, ILO-1, VSTEP
               KLN = MIN( VSTEP, ILO-KKROW )
               LLDTMP = NUMROC( KLN, LWORK / NH, MYROW, 0, NPROW )
               LLDTMP = MAX( 1, LLDTMP )
               CALL DESCINIT( DESCWV, KLN, NH, LWORK / NH, NH, 0, 0,
     $              CONTXT, LLDTMP, IERR )
               CALL PSGEMM( 'N', 'N', KLN, NH, NH, ONE, A, KKROW,
     $              ILO, DESCA, V, 1, 1, DESCV, ZERO, WORK, 1, 1,
     $              DESCWV )
               CALL PSGEMR2D( KLN, NH, WORK, 1, 1, DESCWV, A, KKROW,
     $              ILO, DESCA, CONTXT )
  140       CONTINUE
         END IF
*
*        Update vertical slab in Z.
*
         IF( WANTZ ) THEN
            DO 150 KKROW = ILOZ, IHIZ, VSTEP
               KLN = MIN( VSTEP, IHIZ-KKROW+1 )
               LLDTMP = NUMROC( KLN, LWORK / NH, MYROW, 0, NPROW )
               LLDTMP = MAX( 1, LLDTMP )
               CALL DESCINIT( DESCWV, KLN, NH, LWORK / NH, NH, 0, 0,
     $              CONTXT, LLDTMP, IERR )
               CALL PSGEMM( 'N', 'N', KLN, NH, NH, ONE, Z, KKROW,
     $              ILO, DESCZ, V, 1, 1, DESCV, ZERO, WORK, 1, 1,
     $              DESCWV )
               CALL PSGEMR2D( KLN, NH, WORK, 1, 1, DESCWV, Z,
     $              KKROW, ILO, DESCZ, CONTXT )
  150       CONTINUE
         END IF
      END IF
*
*     END OF PSLAQR4
*
      END
