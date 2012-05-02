      SUBROUTINE PSLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, A, DESCA,
     $                    ILOZ, IHIZ, Z, DESCZ, NS, ND, SR, SI, T, LDT,
     $                    V, LDV, WR, WI, WORK, LWORK )
*
*     Contribution from the Department of Computing Science and HPC2N,
*     Umea University, Sweden
*
*  -- ScaLAPACK routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDT, LDV, LWORK, N, ND,
     $                   NS, NW
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCZ( * )
      REAL               A( * ), SI( KBOT ), SR( KBOT ), T( LDT, * ),
     $                   V( LDV, * ), WORK( * ), WI( * ), WR( * ),
     $                   Z( * )
*     ..
*
*  Purpose
*  =======
*
*  Aggressive early deflation:
*
*  PSLAQR2 accepts as input an upper Hessenberg matrix A and performs an
*  orthogonal similarity transformation designed to detect and deflate
*  fully converged eigenvalues from a trailing principal submatrix.  On
*  output A has been overwritten by a new Hessenberg matrix that is a
*  perturbation of an orthogonal similarity transformation of A.  It is
*  to be hoped that the final version of H has many zero subdiagonal
*  entries.
*
*  This routine handles small deflation windows which is affordable by
*  one processor. Normally, it is called by PSLAQR1. All the inputs are
*  assumed to be valid without checking.
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
*  KBOT    (global input) INTEGER
*          It is assumed without a check that either
*          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
*          determine an isolated block along the diagonal of the
*          Hessenberg matrix. However, H(KTOP,KTOP-1)=0 is not
*          essentially necessary if WANTT is .TRUE. .
*
*  NW      (global input) INTEGER
*          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
*          Normally NW .GE. 3 if PSLAQR2 is called by PSLAQR1.
*
*  A       (local input/output) REAL             array, dimension
*          (DESCH(LLD_),*)
*          On input the initial N-by-N section of A stores the
*          Hessenberg matrix undergoing aggressive early deflation.
*          On output A has been transformed by an orthogonal
*          similarity transformation, perturbed, and the returned
*          to Hessenberg form that (it is to be hoped) has some
*          zero subdiagonal entries.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  ILOZ    (global input) INTEGER
*  IHIZ    (global input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
*
*  Z       (input/output) REAL             array, dimension
*          (DESCH(LLD_),*)
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
*  SR      (global output) REAL             array, dimension KBOT
*  SI      (global output) REAL             array, dimension KBOT
*          On output, the real and imaginary parts of approximate
*          eigenvalues that may be used for shifts are stored in
*          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
*          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
*          On proc #0, the real and imaginary parts of converged
*          eigenvalues are stored in SR(KBOT-ND+1) through SR(KBOT) and
*          SI(KBOT-ND+1) through SI(KBOT), respectively. On other
*          processors, these entries are set to zero.
*
*  T       (local workspace) REAL             array, dimension LDT*NW.
*
*  LDT     (local input) INTEGER
*          The leading dimension of the array T.
*          LDT >= NW.
*
*  V       (local workspace) REAL             array, dimension LDV*NW.
*
*  LDV     (local input) INTEGER
*          The leading dimension of the array V.
*          LDV >= NW.
*
*  WR      (local workspace) REAL             array, dimension KBOT.
*  WI      (local workspace) REAL             array, dimension KBOT.
*
*  WORK    (local workspace) REAL             array, dimension LWORK.
*
*  LWORK   (local input) INTEGER
*          WORK(LWORK) is a local array and LWORK is assumed big enough
*          so that LWORK >= NW*NW.
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
     $                   ICOL2, INFO, II, IROW, IROW1, IROW2, ITMP1,
     $                   ITMP2, J, JAFIRST, JJ, K, L, LDA, LDZ, LLDTMP,
     $                   MYCOL, MYROW, NODE, NPCOL, NPROW, DBLK,
     $                   HSTEP, VSTEP, KKROW, KKCOL, KLN, LTOP, LEFT,
     $                   RIGHT, UP, DOWN, D1, D2
*     ..
*     .. Local Arrays ..
      INTEGER            DESCT( 9 ), DESCV( 9 ), DESCWH( 9 ),
     $                   DESCWV( 9 )
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, SLASET,
     $                   SLAQR3, DESCINIT, PSGEMM, PSGEMR2D, SGEMM,
     $                   SLAMOV, SGESD2D, SGERV2D, SGEBS2D, SGEBR2D,
     $                   IGEBS2D, IGEBR2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
      IF( N.EQ.0 )
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
      I = KBOT
      L = KTOP
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
*     Begin Aggressive Early Deflation.
*
      DBLK = NW
      CALL INFOG2L( I-DBLK+1, I-DBLK+1, DESCA, NPROW, NPCOL, MYROW,
     $     MYCOL, IROW, ICOL, II, JJ )
      IF ( MYROW .EQ. II ) THEN
         CALL DESCINIT( DESCT, DBLK, DBLK, DBLK, DBLK, II, JJ, CONTXT,
     $        LDT, INFO )
         CALL DESCINIT( DESCV, DBLK, DBLK, DBLK, DBLK, II, JJ, CONTXT,
     $        LDV, INFO )
      ELSE
         CALL DESCINIT( DESCT, DBLK, DBLK, DBLK, DBLK, II, JJ, CONTXT,
     $        1, INFO )
         CALL DESCINIT( DESCV, DBLK, DBLK, DBLK, DBLK, II, JJ, CONTXT,
     $        1, INFO )
      END IF
      CALL PSGEMR2D( DBLK, DBLK, A, I-DBLK+1, I-DBLK+1, DESCA, T, 1, 1,
     $     DESCT, CONTXT )
      IF ( MYROW .EQ. II .AND. MYCOL .EQ. JJ ) THEN
         CALL SLASET( 'All', DBLK, DBLK, ZERO, ONE, V, LDV )
         CALL SLAQR3( .TRUE., .TRUE., DBLK, 1, DBLK, DBLK-1, T, LDT, 1,
     $        DBLK, V, LDV, NS, ND, WR, WI, WORK, DBLK, DBLK,
     $        WORK( DBLK*DBLK+1 ), DBLK, DBLK, WORK( 2*DBLK*DBLK+1 ),
     $        DBLK, WORK( 3*DBLK*DBLK+1 ), LWORK-3*DBLK*DBLK )
         CALL SGEBS2D( CONTXT, 'All', ' ', DBLK, DBLK, V, LDV )
         CALL IGEBS2D( CONTXT, 'All', ' ', 1, 1, ND, 1 )
      ELSE
         CALL SGEBR2D( CONTXT, 'All', ' ', DBLK, DBLK, V, LDV, II, JJ )
         CALL IGEBR2D( CONTXT, 'All', ' ', 1, 1, ND, 1, II, JJ )
      END IF
*
      IF( ND .GT. 0 ) THEN
*
*        Copy the local matrix back to the diagonal block.
*
         CALL PSGEMR2D( DBLK, DBLK, T, 1, 1, DESCT, A, I-DBLK+1,
     $        I-DBLK+1, DESCA, CONTXT )
*
*        Update T and Z.
*
         IF( MOD( I-DBLK, HBL )+DBLK .LE. HBL ) THEN
*
*           Simplest case: the deflation window is located on one
*           processor.
*           Call SGEMM directly to perform the update.
*
            HSTEP = LWORK / DBLK
            VSTEP = HSTEP
*
*           Update horizontal slab in A.
*
            IF( WANTT ) THEN
               CALL INFOG2L( I-DBLK+1, I+1, DESCA, NPROW, NPCOL, MYROW,
     $              MYCOL, IROW, ICOL, II, JJ )
               IF( MYROW .EQ. II ) THEN
                  ICOL1 = NUMROC( N, HBL, MYCOL, JAFIRST, NPCOL )
                  DO 10 KKCOL = ICOL, ICOL1, HSTEP
                     KLN = MIN( HSTEP, ICOL1-KKCOL+1 )
                     CALL SGEMM( 'T', 'N', DBLK, KLN, DBLK, ONE, V,
     $                    LDV, A( IROW+(KKCOL-1)*LDA ), LDA, ZERO, WORK,
     $                    DBLK )
                     CALL SLAMOV( 'A', DBLK, KLN, WORK, DBLK,
     $                    A( IROW+(KKCOL-1)*LDA ), LDA )
   10             CONTINUE
               END IF
            END IF
*
*           Update vertical slab in A.
*
            CALL INFOG2L( LTOP, I-DBLK+1, DESCA, NPROW, NPCOL, MYROW,
     $           MYCOL, IROW, ICOL, II, JJ )
            IF( MYCOL .EQ. JJ ) THEN
               CALL INFOG2L( I-DBLK, I-DBLK+1, DESCA, NPROW, NPCOL,
     $              MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
               IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
               DO 20 KKROW = IROW, IROW1, VSTEP
                  KLN = MIN( VSTEP, IROW1-KKROW+1 )
                  CALL SGEMM( 'N', 'N', KLN, DBLK, DBLK, ONE,
     $                 A( KKROW+(ICOL-1)*LDA ), LDA, V, LDV, ZERO, WORK,
     $                 KLN )
                  CALL SLAMOV( 'A', KLN, DBLK, WORK, KLN,
     $                 A( KKROW+(ICOL-1)*LDA ), LDA )
   20          CONTINUE
            END IF
*
*           Update vertical slab in Z.
*
            IF( WANTZ ) THEN
               CALL INFOG2L( ILOZ, I-DBLK+1, DESCZ, NPROW, NPCOL, MYROW,
     $              MYCOL, IROW, ICOL, II, JJ )
               IF( MYCOL .EQ. JJ ) THEN
                  CALL INFOG2L( IHIZ, I-DBLK+1, DESCZ, NPROW, NPCOL,
     $                 MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 30 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL SGEMM( 'N', 'N', KLN, DBLK, DBLK, ONE,
     $                    Z( KKROW+(ICOL-1)*LDZ ), LDZ, V, LDV, ZERO,
     $                    WORK, KLN )
                     CALL SLAMOV( 'A', KLN, DBLK, WORK, KLN,
     $                    Z( KKROW+(ICOL-1)*LDZ ), LDZ )
   30             CONTINUE
               END IF
            END IF
*
         ELSE IF( MOD( I-DBLK, HBL )+DBLK .LE. 2*HBL ) THEN
*
*           More complicated case: the deflation window lay on a 2x2
*           processor mesh.
*           Call SGEMM locally and communicate by pair.
*
            D1 = HBL - MOD( I-DBLK, HBL )
            D2 = DBLK - D1
            HSTEP = LWORK / DBLK
            VSTEP = HSTEP
*
*           Update horizontal slab in A.
*
            IF( WANTT ) THEN
               CALL INFOG2L( I-DBLK+1, I+1, DESCA, NPROW, NPCOL, MYROW,
     $              MYCOL, IROW, ICOL, II, JJ )
               IF( MYROW .EQ. UP ) THEN
                  IF( MYROW .EQ. II ) THEN
                     ICOL1 = NUMROC( N, HBL, MYCOL, JAFIRST, NPCOL )
                     DO 40 KKCOL = ICOL, ICOL1, HSTEP
                        KLN = MIN( HSTEP, ICOL1-KKCOL+1 )
                        CALL SGEMM( 'T', 'N', DBLK, KLN, DBLK, ONE, V,
     $                       DBLK, A( IROW+(KKCOL-1)*LDA ), LDA, ZERO,
     $                       WORK, DBLK )
                        CALL SLAMOV( 'A', DBLK, KLN, WORK, DBLK,
     $                       A( IROW+(KKCOL-1)*LDA ), LDA )
   40                CONTINUE
                  END IF
               ELSE
                  IF( MYROW .EQ. II ) THEN
                     ICOL1 = NUMROC( N, HBL, MYCOL, JAFIRST, NPCOL )
                     DO 50 KKCOL = ICOL, ICOL1, HSTEP
                        KLN = MIN( HSTEP, ICOL1-KKCOL+1 )
                        CALL SGEMM( 'T', 'N', D2, KLN, D1, ONE,
     $                       V( 1, D1+1 ), LDV, A( IROW+(KKCOL-1)*LDA ),
     $                       LDA, ZERO, WORK( D1+1 ), DBLK )
                        CALL SGESD2D( CONTXT, D2, KLN, WORK( D1+1 ),
     $                       DBLK, DOWN, MYCOL )
                        CALL SGERV2D( CONTXT, D1, KLN, WORK, DBLK, DOWN,
     $                       MYCOL )
                        CALL SGEMM( 'T', 'N', D1, KLN, D1, ONE,
     $                       V, LDV, A( IROW+(KKCOL-1)*LDA ), LDA, ONE,
     $                       WORK, DBLK )
                        CALL SLAMOV( 'A', D1, KLN, WORK, DBLK,
     $                       A( IROW+(KKCOL-1)*LDA ), LDA )
   50                CONTINUE
                  ELSE IF( UP .EQ. II ) THEN
                     ICOL1 = NUMROC( N, HBL, MYCOL, JAFIRST, NPCOL )
                     DO 60 KKCOL = ICOL, ICOL1, HSTEP
                        KLN = MIN( HSTEP, ICOL1-KKCOL+1 )
                        CALL SGEMM( 'T', 'N', D1, KLN, D2, ONE,
     $                       V( D1+1, 1 ), LDV, A( IROW+(KKCOL-1)*LDA ),
     $                       LDA, ZERO, WORK, DBLK )
                        CALL SGESD2D( CONTXT, D1, KLN, WORK, DBLK, UP,
     $                       MYCOL )
                        CALL SGERV2D( CONTXT, D2, KLN, WORK( D1+1 ),
     $                       DBLK, UP, MYCOL )
                        CALL SGEMM( 'T', 'N', D2, KLN, D2, ONE,
     $                       V( D1+1, D1+1 ), LDV,
     $                       A( IROW+(KKCOL-1)*LDA ), LDA, ONE,
     $                       WORK( D1+1 ), DBLK )
                        CALL SLAMOV( 'A', D2, KLN, WORK( D1+1 ), DBLK,
     $                       A( IROW+(KKCOL-1)*LDA ), LDA )
   60                CONTINUE
                  END IF
               END IF
            END IF
*
*           Update vertical slab in A.
*
            CALL INFOG2L( LTOP, I-DBLK+1, DESCA, NPROW, NPCOL, MYROW,
     $           MYCOL, IROW, ICOL, II, JJ )
            IF( MYCOL .EQ. LEFT ) THEN
               IF( MYCOL .EQ. JJ ) THEN
                  CALL INFOG2L( I-DBLK, I-DBLK+1, DESCA, NPROW, NPCOL,
     $                 MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 70 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL SGEMM( 'N', 'N', KLN, DBLK, DBLK, ONE,
     $                    A( KKROW+(ICOL-1)*LDA ), LDA, V, LDV, ZERO,
     $                    WORK, KLN )
                     CALL SLAMOV( 'A', KLN, DBLK, WORK, KLN,
     $                    A( KKROW+(ICOL-1)*LDA ), LDA )
   70             CONTINUE
               END IF
            ELSE
               IF( MYCOL .EQ. JJ ) THEN
                  CALL INFOG2L( I-DBLK, I-DBLK+1, DESCA, NPROW, NPCOL,
     $                 MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                  IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                  DO 80 KKROW = IROW, IROW1, VSTEP
                     KLN = MIN( VSTEP, IROW1-KKROW+1 )
                     CALL SGEMM( 'N', 'N', KLN, D2, D1, ONE,
     $                    A( KKROW+(ICOL-1)*LDA ), LDA,
     $                    V( 1, D1+1 ), LDV, ZERO, WORK( 1+D1*KLN ),
     $                    KLN )
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
                  CALL INFOG2L( I-DBLK, I-DBLK+1, DESCA, NPROW, NPCOL,
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
*
*           Update vertical slab in Z.
*
            IF( WANTZ ) THEN
               CALL INFOG2L( ILOZ, I-DBLK+1, DESCZ, NPROW, NPCOL, MYROW,
     $              MYCOL, IROW, ICOL, II, JJ )
               IF( MYCOL .EQ. LEFT ) THEN
                  IF( MYCOL .EQ. JJ ) THEN
                     CALL INFOG2L( IHIZ, I-DBLK+1, DESCZ, NPROW, NPCOL,
     $                    MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                     IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                     DO 100 KKROW = IROW, IROW1, VSTEP
                        KLN = MIN( VSTEP, IROW1-KKROW+1 )
                        CALL SGEMM( 'N', 'N', KLN, DBLK, DBLK, ONE,
     $                       Z( KKROW+(ICOL-1)*LDZ ), LDZ, V, LDV, ZERO,
     $                       WORK, KLN )
                        CALL SLAMOV( 'A', KLN, DBLK, WORK, KLN,
     $                       Z( KKROW+(ICOL-1)*LDZ ), LDZ )
  100                CONTINUE
                  END IF
               ELSE
                  IF( MYCOL .EQ. JJ ) THEN
                     CALL INFOG2L( IHIZ, I-DBLK+1, DESCZ, NPROW, NPCOL,
     $                    MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                     IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                     DO 110 KKROW = IROW, IROW1, VSTEP
                        KLN = MIN( VSTEP, IROW1-KKROW+1 )
                        CALL SGEMM( 'N', 'N', KLN, D2, D1, ONE,
     $                       Z( KKROW+(ICOL-1)*LDZ ), LDZ,
     $                       V( 1, D1+1 ), LDV, ZERO, WORK( 1+D1*KLN ),
     $                       KLN )
                        CALL SGESD2D( CONTXT, KLN, D2, WORK( 1+D1*KLN ),
     $                       KLN, MYROW, RIGHT )
                        CALL SGERV2D( CONTXT, KLN, D1, WORK, KLN, MYROW,
     $                       RIGHT )
                        CALL SGEMM( 'N', 'N', KLN, D1, D1, ONE,
     $                       Z( KKROW+(ICOL-1)*LDZ ), LDZ, V, LDV, ONE,
     $                       WORK, KLN )
                        CALL SLAMOV( 'A', KLN, D1, WORK, KLN,
     $                       Z( KKROW+(ICOL-1)*LDZ ), LDZ )
  110                CONTINUE
                  ELSE IF( LEFT .EQ. JJ ) THEN
                     CALL INFOG2L( IHIZ, I-DBLK+1, DESCZ, NPROW, NPCOL,
     $                    MYROW, MYCOL, IROW1, ICOL1, ITMP1, ITMP2 )
                     IF( MYROW .NE. ITMP1 ) IROW1 = IROW1-1
                     DO 120 KKROW = IROW, IROW1, VSTEP
                        KLN = MIN( VSTEP, IROW1-KKROW+1 )
                        CALL SGEMM( 'N', 'N', KLN, D1, D2, ONE,
     $                       Z( KKROW+(ICOL-1)*LDZ ), LDZ,
     $                       V( D1+1, 1 ), LDV, ZERO, WORK, KLN )
                        CALL SGESD2D( CONTXT, KLN, D1, WORK, KLN, MYROW,
     $                       LEFT )
                        CALL SGERV2D( CONTXT, KLN, D2, WORK( 1+D1*KLN ),
     $                       KLN, MYROW, LEFT )
                        CALL SGEMM( 'N', 'N', KLN, D2, D2, ONE,
     $                       Z( KKROW+(ICOL-1)*LDZ ), LDZ,
     $                       V( D1+1, D1+1 ), LDV, ONE,
     $                       WORK( 1+D1*KLN ), KLN )
                        CALL SLAMOV( 'A', KLN, D2, WORK( 1+D1*KLN ),
     $                       KLN, Z( KKROW+(ICOL-1)*LDZ ), LDZ )
  120                CONTINUE
                  END IF
               END IF
            END IF
*
         ELSE
*
*           Most complicated case: the deflation window lay across the
*           border of the processor mesh.
*           Treat V as a distributed matrix and call PSGEMM.
*
            HSTEP = LWORK / DBLK * NPCOL
            VSTEP = LWORK / DBLK * NPROW
            LLDTMP = NUMROC( DBLK, DBLK, MYROW, 0, NPROW )
            LLDTMP = MAX( 1, LLDTMP )
            CALL DESCINIT( DESCV, DBLK, DBLK, DBLK, DBLK, 0, 0, CONTXT,
     $           LLDTMP, INFO )
            CALL DESCINIT( DESCWH, DBLK, HSTEP, DBLK, LWORK / DBLK, 0,
     $           0, CONTXT, LLDTMP, INFO )
*
*           Update horizontal slab in A.
*
            IF( WANTT ) THEN
               DO 130 KKCOL = I+1, N, HSTEP
                  KLN = MIN( HSTEP, N-KKCOL+1 )
                  CALL PSGEMM( 'T', 'N', DBLK, KLN, DBLK, ONE, V, 1, 1,
     $                 DESCV, A, I-DBLK+1, KKCOL, DESCA, ZERO, WORK, 1,
     $                 1, DESCWH )
                  CALL PSGEMR2D( DBLK, KLN, WORK, 1, 1, DESCWH, A,
     $                 I-DBLK+1, KKCOL, DESCA, CONTXT )
  130          CONTINUE
            END IF
*
*           Update vertical slab in A.
*
            DO 140 KKROW = LTOP, I-DBLK, VSTEP
               KLN = MIN( VSTEP, I-DBLK-KKROW+1 )
               LLDTMP = NUMROC( KLN, LWORK / DBLK, MYROW, 0, NPROW )
               LLDTMP = MAX( 1, LLDTMP )
               CALL DESCINIT( DESCWV, KLN, DBLK, LWORK / DBLK, DBLK, 0,
     $              0, CONTXT, LLDTMP, INFO )
               CALL PSGEMM( 'N', 'N', KLN, DBLK, DBLK, ONE, A, KKROW,
     $              I-DBLK+1, DESCA, V, 1, 1, DESCV, ZERO, WORK, 1, 1,
     $              DESCWV )
               CALL PSGEMR2D( KLN, DBLK, WORK, 1, 1, DESCWV, A, KKROW,
     $              I-DBLK+1, DESCA, CONTXT )
  140       CONTINUE
*
*           Update vertical slab in Z.
*
            IF( WANTZ ) THEN
               DO 150 KKROW = ILOZ, IHIZ, VSTEP
                  KLN = MIN( VSTEP, IHIZ-KKROW+1 )
                  LLDTMP = NUMROC( KLN, LWORK / DBLK, MYROW, 0, NPROW )
                  LLDTMP = MAX( 1, LLDTMP )
                  CALL DESCINIT( DESCWV, KLN, DBLK, LWORK / DBLK, DBLK,
     $                 0, 0, CONTXT, LLDTMP, INFO )
                  CALL PSGEMM( 'N', 'N', KLN, DBLK, DBLK, ONE, Z, KKROW,
     $                 I-DBLK+1, DESCZ, V, 1, 1, DESCV, ZERO, WORK, 1,
     $                 1, DESCWV )
                  CALL PSGEMR2D( KLN, DBLK, WORK, 1, 1, DESCWV, Z,
     $                 KKROW, I-DBLK+1, DESCZ, CONTXT )
  150          CONTINUE
            END IF
         END IF
*
*        Extract converged eigenvalues.
*
         II = 0
  160    CONTINUE
            IF( II .EQ. ND-1 .OR. WI( DBLK-II ) .EQ. ZERO ) THEN
               IF( NODE .EQ. 0 ) THEN
                  SR( I-II ) = WR( DBLK-II )
               ELSE
                  SR( I-II ) = ZERO
               END IF
               SI( I-II ) = ZERO
               II = II + 1
            ELSE
               IF( NODE .EQ. 0 ) THEN
                  SR( I-II-1 ) = WR( DBLK-II-1 )
                  SR( I-II ) = WR( DBLK-II )
                  SI( I-II-1 ) = WI( DBLK-II-1 )
                  SI( I-II ) = WI( DBLK-II )
               ELSE
                  SR( I-II-1 ) = ZERO
                  SR( I-II ) = ZERO
                  SI( I-II-1 ) = ZERO
                  SI( I-II ) = ZERO
               END IF
               II = II + 2
            END IF
         IF( II .LT. ND ) GOTO 160
      END IF
*
*     END OF PSLAQR2
*
      END
