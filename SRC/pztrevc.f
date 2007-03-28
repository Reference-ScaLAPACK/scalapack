      SUBROUTINE PZTREVC( SIDE, HOWMNY, SELECT, N, T, DESCT, VL, DESCVL,
     $                    VR, DESCVR, MM, M, WORK, RWORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     July 31, 2001
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      INTEGER            DESCT( * ), DESCVL( * ), DESCVR( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         T( * ), VL( * ), VR( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZTREVC computes some or all of the right and/or left eigenvectors of
*  a complex upper triangular matrix T in parallel.
*
*  The right eigenvector x and the left eigenvector y of T corresponding
*  to an eigenvalue w are defined by:
*
*               T*x = w*x,     y'*T = w*y'
*
*  where y' denotes the conjugate transpose of the vector y.
*
*  If all eigenvectors are requested, the routine may either return the
*  matrices X and/or Y of right or left eigenvectors of T, or the
*  products Q*X and/or Q*Y, where Q is an input unitary
*  matrix. If T was obtained from the Schur factorization of an
*  original matrix A = Q*T*Q', then Q*X and Q*Y are the matrices of
*  right or left eigenvectors of A.
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
*  Arguments
*  =========
*
*  SIDE    (global input) CHARACTER*1
*          = 'R':  compute right eigenvectors only;
*          = 'L':  compute left eigenvectors only;
*          = 'B':  compute both right and left eigenvectors.
*
*  HOWMNY  (global input) CHARACTER*1
*          = 'A':  compute all right and/or left eigenvectors;
*          = 'B':  compute all right and/or left eigenvectors,
*                  and backtransform them using the input matrices
*                  supplied in VR and/or VL;
*          = 'S':  compute selected right and/or left eigenvectors,
*                  specified by the logical array SELECT.
*
*  SELECT  (global input) LOGICAL array, dimension (N)
*          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
*          computed.
*          If HOWMNY = 'A' or 'B', SELECT is not referenced.
*          To select the eigenvector corresponding to the j-th
*          eigenvalue, SELECT(j) must be set to .TRUE..
*
*  N       (global input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (global input/output) COMPLEX*16 array, dimension
*          (DESCT(LLD_),*)
*          The upper triangular matrix T.  T is modified, but restored
*          on exit.
*
*  DESCT   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix T.
*
*  VL      (global input/output) COMPLEX*16 array, dimension
*          (DESCVL(LLD_),MM)
*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
*          contain an N-by-N matrix Q (usually the unitary matrix Q of
*          Schur vectors returned by ZHSEQR).
*          On exit, if SIDE = 'L' or 'B', VL contains:
*          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
*          if HOWMNY = 'B', the matrix Q*Y;
*          if HOWMNY = 'S', the left eigenvectors of T specified by
*                           SELECT, stored consecutively in the columns
*                           of VL, in the same order as their
*                           eigenvalues.
*          If SIDE = 'R', VL is not referenced.
*
*  DESCVL  (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix VL.
*
*  VR      (global input/output) COMPLEX*16 array, dimension
*          (DESCVR(LLD_),MM)
*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
*          contain an N-by-N matrix Q (usually the unitary matrix Q of
*          Schur vectors returned by ZHSEQR).
*          On exit, if SIDE = 'R' or 'B', VR contains:
*          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
*          if HOWMNY = 'B', the matrix Q*X;
*          if HOWMNY = 'S', the right eigenvectors of T specified by
*                           SELECT, stored consecutively in the columns
*                           of VR, in the same order as their
*                           eigenvalues.
*          If SIDE = 'L', VR is not referenced.
*
*  DESCVR  (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix VR.
*
*  MM      (global input) INTEGER
*          The number of columns in the arrays VL and/or VR. MM >= M.
*
*  M       (global output) INTEGER
*          The number of columns in the arrays VL and/or VR actually
*          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M
*          is set to N.  Each selected eigenvector occupies one
*          column.
*
*  WORK    (local workspace) COMPLEX*16 array,
*                                         dimension ( 2*DESCT(LLD_) )
*          Additional workspace may be required if PZLATTRS is updated
*          to use WORK.
*
*  RWORK   (local workspace) DOUBLE PRECISION array,
*                                          dimension ( DESCT(LLD_) )
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The algorithm used in this program is basically backward (forward)
*  substitution.  It is the hope that scaling would be used to make the
*  the code robust against possible overflow.  But scaling has not yet
*  been implemented in PZLATTRS which is called by this routine to solve
*  the triangular systems.  PZLATTRS just calls PZTRSV.
*
*  Each eigenvector is normalized so that the element of largest
*  magnitude has magnitude 1; here the magnitude of a complex number
*  (x,y) is taken to be |x| + |y|.
*
*  Further Details
*  ===============
*
*  Implemented by Mark R. Fahey, June, 2000
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, OVER, RIGHTV, SOMEV
      INTEGER            CONTXT, CSRC, I, ICOL, II, IROW, IS, ITMP1,
     $                   ITMP2, J, K, KI, LDT, LDVL, LDVR, LDW, MB,
     $                   MYCOL, MYROW, NB, NPCOL, NPROW, RSRC
      REAL               SELF
      DOUBLE PRECISION   OVFL, REMAXD, SCALE, SMIN, SMLNUM, ULP, UNFL
      COMPLEX*16         CDUM, REMAXC, SHIFT
*     ..
*     .. Local Arrays ..
      INTEGER            DESCW( DLEN_ )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           LSAME, PDLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCINIT, DGSUM2D, IGAMN2D,
     $                   INFOG2L, PDLABAD, PDZASUM, PXERBLA, PZAMAX,
     $                   PZCOPY, PZDSCAL, PZGEMV, PZLASET, PZLATTRS,
     $                   ZGSUM2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      CONTXT = DESCT( CTXT_ )
      RSRC = DESCT( RSRC_ )
      CSRC = DESCT( CSRC_ )
      MB = DESCT( MB_ )
      NB = DESCT( NB_ )
      LDT = DESCT( LLD_ )
      LDW = LDT
      LDVR = DESCVR( LLD_ )
      LDVL = DESCVL( LLD_ )
*
      CALL BLACS_GRIDINFO( CONTXT, NPROW, NPCOL, MYROW, MYCOL )
      SELF = MYROW*NPCOL + MYCOL
*
*     Decode and test the input parameters
*
      BOTHV = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV = LSAME( SIDE, 'L' ) .OR. BOTHV
*
      ALLV = LSAME( HOWMNY, 'A' )
      OVER = LSAME( HOWMNY, 'B' ) .OR. LSAME( HOWMNY, 'O' )
      SOMEV = LSAME( HOWMNY, 'S' )
*
*     Set M to the number of columns required to store the selected
*     eigenvectors.
*
      IF( SOMEV ) THEN
         M = 0
         DO 10 J = 1, N
            IF( SELECT( J ) )
     $         M = M + 1
   10    CONTINUE
      ELSE
         M = N
      END IF
*
      INFO = 0
      IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -1
      ELSE IF( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( MM.LT.M ) THEN
         INFO = -11
      END IF
      CALL IGAMN2D( CONTXT, 'ALL', ' ', 1, 1, INFO, 1, ITMP1, ITMP2, -1,
     $              -1, -1 )
      IF( INFO.LT.0 ) THEN
         CALL PXERBLA( CONTXT, 'PZTREVC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Set the constants to control overflow.
*
      UNFL = PDLAMCH( CONTXT, 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL PDLABAD( CONTXT, UNFL, OVFL )
      ULP = PDLAMCH( CONTXT, 'Precision' )
      SMLNUM = UNFL*( N / ULP )
*
*     Store the diagonal elements of T in working array WORK( LDW+1 ).
*
      DO 20 I = 1, N
         CALL INFOG2L( I, I, DESCT, NPROW, NPCOL, MYROW, MYCOL, IROW,
     $                 ICOL, ITMP1, ITMP2 )
         IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
            WORK( LDW+IROW ) = T( ( ICOL-1 )*LDT+IROW )
         END IF
   20 CONTINUE
*
*     Compute 1-norm of each column of strictly upper triangular
*     part of T to control overflow in triangular solver.  Computed,
*     but not used.  For use in PZLATTRS.
*
      RWORK( 1 ) = ZERO
      DO 30 J = 2, N
         CALL PDZASUM( J-1, RWORK( J ), T, 1, J, DESCT, 1 )
   30 CONTINUE
*     I replicate the norms in RWORK.  Should they be distributed
*     over the process rows?
      CALL DGSUM2D( CONTXT, 'Row', ' ', N, 1, RWORK, N, -1, -1 )
*
      IF( RIGHTV ) THEN
*
*        Compute right eigenvectors.
*
*        Need to set the distribution pattern of WORK
*
         CALL DESCINIT( DESCW, N, 1, NB, 1, RSRC, CSRC, CONTXT, LDW,
     $                  INFO )
*
         IS = M
         DO 70 KI = N, 1, -1
*
            IF( SOMEV ) THEN
               IF( .NOT.SELECT( KI ) )
     $            GO TO 70
            END IF
*
            SMIN = ZERO
            SHIFT = CZERO
            CALL INFOG2L( KI, KI, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $                    IROW, ICOL, ITMP1, ITMP2 )
            IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
               SHIFT = T( ( ICOL-1 )*LDT+IROW )
               SMIN = MAX( ULP*( CABS1( SHIFT ) ), SMLNUM )
            END IF
            CALL DGSUM2D( CONTXT, 'ALL', ' ', 1, 1, SMIN, 1, -1, -1 )
            CALL ZGSUM2D( CONTXT, 'ALL', ' ', 1, 1, SHIFT, 1, -1, -1 )
*
            CALL INFOG2L( 1, 1, DESCW, NPROW, NPCOL, MYROW, MYCOL, IROW,
     $                    ICOL, ITMP1, ITMP2 )
            IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
               WORK( 1 ) = CONE
            END IF
*
*           Form right-hand side.  Distribute rhs onto first column
*           of processor grid.
*
            IF( KI.GT.1 ) THEN
               CALL PZCOPY( KI-1, T, 1, KI, DESCT, 1, WORK, 1, 1, DESCW,
     $                      1 )
            END IF
            DO 40 K = 1, KI - 1
               CALL INFOG2L( K, 1, DESCW, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW, ICOL, ITMP1, ITMP2 )
               IF( MYROW.EQ.ITMP1 .AND. MYCOL.EQ.ITMP2 ) THEN
                  WORK( IROW ) = -WORK( IROW )
               END IF
   40       CONTINUE
*
*           Solve the triangular system:
*              (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.
*
            DO 50 K = 1, KI - 1
               CALL INFOG2L( K, K, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW, ICOL, ITMP1, ITMP2 )
               IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
                  T( ( ICOL-1 )*LDT+IROW ) = T( ( ICOL-1 )*LDT+IROW ) -
     $               SHIFT
                  IF( CABS1( T( ( ICOL-1 )*LDT+IROW ) ).LT.SMIN ) THEN
                     T( ( ICOL-1 )*LDT+IROW ) = DCMPLX( SMIN )
                  END IF
               END IF
   50       CONTINUE
*
            IF( KI.GT.1 ) THEN
               CALL PZLATTRS( 'Upper', 'No transpose', 'Non-unit', 'Y',
     $                        KI-1, T, 1, 1, DESCT, WORK, 1, 1, DESCW,
     $                        SCALE, RWORK, INFO )
               CALL INFOG2L( KI, 1, DESCW, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW, ICOL, ITMP1, ITMP2 )
               IF( MYROW.EQ.ITMP1 .AND. MYCOL.EQ.ITMP2 ) THEN
                  WORK( IROW ) = DCMPLX( SCALE )
               END IF
            END IF
*
*           Copy the vector x or Q*x to VR and normalize.
*
            IF( .NOT.OVER ) THEN
               CALL PZCOPY( KI, WORK, 1, 1, DESCW, 1, VR, 1, IS, DESCVR,
     $                      1 )
*
               CALL PZAMAX( KI, REMAXC, II, VR, 1, IS, DESCVR, 1 )
               REMAXD = ONE / MAX( CABS1( REMAXC ), UNFL )
               CALL PZDSCAL( KI, REMAXD, VR, 1, IS, DESCVR, 1 )
*
               CALL PZLASET( ' ', N-KI, 1, CZERO, CZERO, VR, KI+1, IS,
     $                       DESCVR )
            ELSE
               IF( KI.GT.1 )
     $            CALL PZGEMV( 'N', N, KI-1, CONE, VR, 1, 1, DESCVR,
     $                         WORK, 1, 1, DESCW, 1, DCMPLX( SCALE ),
     $                         VR, 1, KI, DESCVR, 1 )
*
               CALL PZAMAX( N, REMAXC, II, VR, 1, KI, DESCVR, 1 )
               REMAXD = ONE / MAX( CABS1( REMAXC ), UNFL )
               CALL PZDSCAL( N, REMAXD, VR, 1, KI, DESCVR, 1 )
            END IF
*
*           Set back the original diagonal elements of T.
*
            DO 60 K = 1, KI - 1
               CALL INFOG2L( K, K, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW, ICOL, ITMP1, ITMP2 )
               IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
                  T( ( ICOL-1 )*LDT+IROW ) = WORK( LDW+IROW )
               END IF
   60       CONTINUE
*
            IS = IS - 1
   70    CONTINUE
      END IF
*
      IF( LEFTV ) THEN
*
*        Compute left eigenvectors.
*
*        Need to set the distribution pattern of WORK
*
         CALL DESCINIT( DESCW, N, 1, MB, 1, RSRC, CSRC, CONTXT, LDW,
     $                  INFO )
*
         IS = 1
         DO 110 KI = 1, N
*
            IF( SOMEV ) THEN
               IF( .NOT.SELECT( KI ) )
     $            GO TO 110
            END IF
*
            SMIN = ZERO
            SHIFT = CZERO
            CALL INFOG2L( KI, KI, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $                    IROW, ICOL, ITMP1, ITMP2 )
            IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
               SHIFT = T( ( ICOL-1 )*LDT+IROW )
               SMIN = MAX( ULP*( CABS1( SHIFT ) ), SMLNUM )
            END IF
            CALL DGSUM2D( CONTXT, 'ALL', ' ', 1, 1, SMIN, 1, -1, -1 )
            CALL ZGSUM2D( CONTXT, 'ALL', ' ', 1, 1, SHIFT, 1, -1, -1 )
*
            CALL INFOG2L( N, 1, DESCW, NPROW, NPCOL, MYROW, MYCOL, IROW,
     $                    ICOL, ITMP1, ITMP2 )
            IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
               WORK( IROW ) = CONE
            END IF
*
*           Form right-hand side.
*
            IF( KI.LT.N ) THEN
               CALL PZCOPY( N-KI, T, KI, KI+1, DESCT, N, WORK, KI+1, 1,
     $                      DESCW, 1 )
            END IF
            DO 80 K = KI + 1, N
               CALL INFOG2L( K, 1, DESCW, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW, ICOL, ITMP1, ITMP2 )
               IF( MYROW.EQ.ITMP1 .AND. MYCOL.EQ.ITMP2 ) THEN
                  WORK( IROW ) = -DCONJG( WORK( IROW ) )
               END IF
   80       CONTINUE
*
*           Solve the triangular system:
*              (T(KI+1:N,KI+1:N) - T(KI,KI))'*X = SCALE*WORK.
*
            DO 90 K = KI + 1, N
               CALL INFOG2L( K, K, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW, ICOL, ITMP1, ITMP2 )
               IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
                  T( ( ICOL-1 )*LDT+IROW ) = T( ( ICOL-1 )*LDT+IROW ) -
     $               SHIFT
                  IF( CABS1( T( ( ICOL-1 )*LDT+IROW ) ).LT.SMIN )
     $               T( ( ICOL-1 )*LDT+IROW ) = DCMPLX( SMIN )
               END IF
   90       CONTINUE
*
            IF( KI.LT.N ) THEN
               CALL PZLATTRS( 'Upper', 'Conjugate transpose', 'Nonunit',
     $                        'Y', N-KI, T, KI+1, KI+1, DESCT, WORK,
     $                        KI+1, 1, DESCW, SCALE, RWORK, INFO )
               CALL INFOG2L( KI, 1, DESCW, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW, ICOL, ITMP1, ITMP2 )
               IF( MYROW.EQ.ITMP1 .AND. MYCOL.EQ.ITMP2 ) THEN
                  WORK( IROW ) = DCMPLX( SCALE )
               END IF
            END IF
*
*           Copy the vector x or Q*x to VL and normalize.
*
            IF( .NOT.OVER ) THEN
               CALL PZCOPY( N-KI+1, WORK, KI, 1, DESCW, 1, VL, KI, IS,
     $                      DESCVL, 1 )
*
               CALL PZAMAX( N-KI+1, REMAXC, II, VL, KI, IS, DESCVL, 1 )
               REMAXD = ONE / MAX( CABS1( REMAXC ), UNFL )
               CALL PZDSCAL( N-KI+1, REMAXD, VL, KI, IS, DESCVL, 1 )
*
               CALL PZLASET( ' ', KI-1, 1, CZERO, CZERO, VL, 1, IS,
     $                       DESCVL )
            ELSE
               IF( KI.LT.N )
     $            CALL PZGEMV( 'N', N, N-KI, CONE, VL, 1, KI+1, DESCVL,
     $                         WORK, KI+1, 1, DESCW, 1, DCMPLX( SCALE ),
     $                         VL, 1, KI, DESCVL, 1 )
*
               CALL PZAMAX( N, REMAXC, II, VL, 1, KI, DESCVL, 1 )
               REMAXD = ONE / MAX( CABS1( REMAXC ), UNFL )
               CALL PZDSCAL( N, REMAXD, VL, 1, KI, DESCVL, 1 )
            END IF
*
*           Set back the original diagonal elements of T.
*
            DO 100 K = KI + 1, N
               CALL INFOG2L( K, K, DESCT, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW, ICOL, ITMP1, ITMP2 )
               IF( ( MYROW.EQ.ITMP1 ) .AND. ( MYCOL.EQ.ITMP2 ) ) THEN
                  T( ( ICOL-1 )*LDT+IROW ) = WORK( LDW+IROW )
               END IF
  100       CONTINUE
*
            IS = IS + 1
  110    CONTINUE
      END IF
*
      RETURN
*
*     End of PZTREVC
*
      END
