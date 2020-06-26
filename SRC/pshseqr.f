      SUBROUTINE PSHSEQR( JOB, COMPZ, N, ILO, IHI, H, DESCH, WR, WI, Z,
     $                    DESCZ, WORK, LWORK, IWORK, LIWORK, INFO )
*
*     Contribution from the Department of Computing Science and HPC2N,
*     Umea University, Sweden
*
*  -- ScaLAPACK driver routine (version 2.0.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     Univ. of Colorado Denver and University of California, Berkeley.
*     January, 2012
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LWORK, LIWORK, N
      CHARACTER          COMPZ, JOB
*     ..
*     .. Array Arguments ..
      INTEGER            DESCH( * ) , DESCZ( * ), IWORK( * )
      REAL               H( * ), WI( N ), WORK( * ), WR( N ), Z( * )
*     ..
*  Purpose
*  =======
*
*  PSHSEQR computes the eigenvalues of an upper Hessenberg matrix H
*  and, optionally, the matrices T and Z from the Schur decomposition
*  H = Z*T*Z**T, where T is an upper quasi-triangular matrix (the
*  Schur form), and Z is the orthogonal matrix of Schur vectors.
*
*  Optionally Z may be postmultiplied into an input orthogonal
*  matrix Q so that this routine can give the Schur factorization
*  of a matrix A which has been reduced to the Hessenberg form H
*  by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
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
*  JOB     (global input) CHARACTER*1
*          = 'E':  compute eigenvalues only;
*          = 'S':  compute eigenvalues and the Schur form T.
*
*  COMPZ   (global input) CHARACTER*1
*          = 'N':  no Schur vectors are computed;
*          = 'I':  Z is initialized to the unit matrix and the matrix Z
*                  of Schur vectors of H is returned;
*          = 'V':  Z must contain an orthogonal matrix Q on entry, and
*                  the product Q*Z is returned.
*
*  N       (global input) INTEGER
*          The order of the Hessenberg matrix H (and Z if WANTZ).
*          N >= 0.
*
*  ILO     (global input) INTEGER
*  IHI     (global input) INTEGER
*          It is assumed that H is already upper triangular in rows
*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*          set by a previous call to PSGEBAL, and then passed to PSGEHRD
*          when the matrix output by PSGEBAL is reduced to Hessenberg
*          form. Otherwise ILO and IHI should be set to 1 and N
*          respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
*          If N = 0, then ILO = 1 and IHI = 0.
*
*  H       (global input/output) REAL             array, dimension
*          (DESCH(LLD_),*)
*          On entry, the upper Hessenberg matrix H.
*          On exit, if JOB = 'S', H is upper quasi-triangular in
*          rows and columns ILO:IHI, with 1-by-1 and 2-by-2 blocks on
*          the main diagonal.  The 2-by-2 diagonal blocks (corresponding
*          to complex conjugate pairs of eigenvalues) are returned in
*          standard form, with H(i,i) = H(i+1,i+1) and
*          H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and JOB = 'E', the
*          contents of H are unspecified on exit.
*
*  DESCH   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix H.
*
*  WR      (global output) REAL             array, dimension (N)
*  WI      (global output) REAL             array, dimension (N)
*          The real and imaginary parts, respectively, of the computed
*          eigenvalues ILO to IHI are stored in the corresponding
*          elements of WR and WI. If two eigenvalues are computed as a
*          complex conjugate pair, they are stored in consecutive
*          elements of WR and WI, say the i-th and (i+1)th, with
*          WI(i) > 0 and WI(i+1) < 0. If JOB = 'S', the
*          eigenvalues are stored in the same order as on the diagonal
*          of the Schur form returned in H.
*
*  Z       (global input/output) REAL             array.
*          If COMPZ = 'V', on entry Z must contain the current
*          matrix Z of accumulated transformations from, e.g., PSGEHRD,
*          and on exit Z has been updated; transformations are applied
*          only to the submatrix Z(ILO:IHI,ILO:IHI).
*          If COMPZ = 'N', Z is not referenced.
*          If COMPZ = 'I', on entry Z need not be set and on exit,
*          if INFO = 0, Z contains the orthogonal matrix Z of the Schur
*          vectors of H.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*  WORK    (local workspace) REAL             array, dimension(LWORK)
*
*  LWORK   (local input) INTEGER
*          The length of the workspace array WORK.
*
*  IWORK   (local workspace) INTEGER array, dimension (LIWORK)
*
*  LIWORK  (local input) INTEGER
*          The length of the workspace array IWORK.
*
*  INFO    (output) INTEGER
*          =    0:  successful exit
*          .LT. 0:  if INFO = -i, the i-th argument had an illegal
*                   value (see also below for -7777 and -8888).
*          .GT. 0:  if INFO = i, PSHSEQR failed to compute all of
*                   the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
*                   and WI contain those eigenvalues which have been
*                   successfully computed.  (Failures are rare.)
*
*                If INFO .GT. 0 and JOB = 'E', then on exit, the
*                remaining unconverged eigenvalues are the eigen-
*                values of the upper Hessenberg matrix rows and
*                columns ILO through INFO of the final, output
*                value of H.
*
*                If INFO .GT. 0 and JOB   = 'S', then on exit
*
*           (*)  (initial value of H)*U  = U*(final value of H)
*
*                where U is an orthogonal matrix.  The final
*                value of H is upper Hessenberg and quasi-triangular
*                in rows and columns INFO+1 through IHI.
*
*                If INFO .GT. 0 and COMPZ = 'V', then on exit
*
*                  (final value of Z)  =  (initial value of Z)*U
*
*                where U is the orthogonal matrix in (*) (regard-
*                less of the value of JOB.)
*
*                If INFO .GT. 0 and COMPZ = 'I', then on exit
*                      (final value of Z)  = U
*                where U is the orthogonal matrix in (*) (regard-
*                less of the value of JOB.)
*
*                If INFO .GT. 0 and COMPZ = 'N', then Z is not
*                accessed.
*
*          = -7777: PSLAQR0 failed to converge and PSLAQR1 was called
*                   instead. This could happen. Mostly due to a bug.
*                   Please, send a bug report to the authors.
*          = -8888: PSLAQR1 failed to converge and PSLAQR0 was called
*                   instead. This should not happen.
*
*     ================================================================
*     Based on contributions by
*        Robert Granat, Department of Computing Science and HPC2N,
*        Umea University, Sweden.
*     ================================================================
*
*     Restrictions: The block size in H and Z must be square and larger
*     than or equal to six (6) due to restrictions in PSLAQR1, PSLAQR5
*     and SLAQR6. Moreover, H and Z need to be distributed identically
*     with the same context.
*
*     ================================================================
*     References:
*       K. Braman, R. Byers, and R. Mathias,
*       The Multi-Shift QR Algorithm Part I: Maintaining Well Focused
*       Shifts, and Level 3 Performance.
*       SIAM J. Matrix Anal. Appl., 23(4):929--947, 2002.
*
*       K. Braman, R. Byers, and R. Mathias,
*       The Multi-Shift QR Algorithm Part II: Aggressive Early
*       Deflation.
*       SIAM J. Matrix Anal. Appl., 23(4):948--973, 2002.
*
*       R. Granat, B. Kagstrom, and D. Kressner,
*       A Novel Parallel QR Algorithm for Hybrid Distributed Momory HPC
*       Systems.
*       SIAM J. Sci. Comput., 32(4):2345--2378, 2010.
*
*     ================================================================
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      LOGICAL            CRSOVER
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9,
     $                     CRSOVER = .TRUE. )
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
      INTEGER            NL
      PARAMETER          ( NL = 49 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0e0, ONE = 1.0e0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, KBOT, NMIN, LLDH, LLDZ, ICTXT, NPROW, NPCOL,
     $                   MYROW, MYCOL, HROWS, HCOLS, IPW, NH, NB,
     $                   II, JJ, HRSRC, HCSRC, NPROCS, ILOC1, JLOC1,
     $                   HRSRC1, HCSRC1, K, ILOC2, JLOC2, ILOC3, JLOC3,
     $                   ILOC4, JLOC4, HRSRC2, HCSRC2, HRSRC3, HCSRC3,
     $                   HRSRC4, HCSRC4, LIWKOPT
      LOGICAL            INITZ, LQUERY, WANTT, WANTZ, PAIR, BORDER
      REAL               TMP1, TMP2, TMP3, TMP4, DUM1, DUM2, DUM3,
     $                   DUM4, ELEM1, ELEM4,
     $                   CS, SN, ELEM5, TMP, LWKOPT
*     ..
*     .. Local Arrays ..
      INTEGER            DESCH2( DLEN_ )
      REAL               ELEM2( 1 ), ELEM3( 1 )
*     ..
*     .. External Functions ..
      INTEGER            PILAENVX, NUMROC, ICEIL
      LOGICAL            LSAME
      EXTERNAL           PILAENVX, LSAME, NUMROC, ICEIL
*     ..
*     .. External Subroutines ..
      EXTERNAL           PSLACPY, PSLAQR1, PSLAQR0, PSLASET, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          FLOAT, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Decode and check the input parameters.
*
      INFO = 0
      ICTXT = DESCH( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NPROCS = NPROW*NPCOL
      IF( NPROW.EQ.-1 ) INFO = -(600+CTXT_)
      IF( INFO.EQ.0 ) THEN
         WANTT = LSAME( JOB, 'S' )
         INITZ = LSAME( COMPZ, 'I' )
         WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
         LLDH = DESCH( LLD_ )
         LLDZ = DESCZ( LLD_ )
         NB = DESCH( MB_ )
         LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
         IF( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) THEN
            INFO = -1
         ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
            INFO = -2
         ELSE IF( N.LT.0 ) THEN
            INFO = -3
         ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
            INFO = -4
         ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
            INFO = -5
         ELSEIF( DESCZ( CTXT_ ).NE.DESCH( CTXT_ ) ) THEN
            INFO = -( 1000+CTXT_ )
         ELSEIF( DESCH( MB_ ).NE.DESCH( NB_ ) ) THEN
            INFO = -( 700+NB_ )
         ELSEIF( DESCZ( MB_ ).NE.DESCZ( NB_ ) ) THEN
            INFO = -( 1000+NB_ )
         ELSEIF( DESCH( MB_ ).NE.DESCZ( MB_ ) ) THEN
            INFO = -( 1000+MB_ )
         ELSEIF( DESCH( MB_ ).LT.6 ) THEN
            INFO = -( 700+NB_ )
         ELSEIF( DESCZ( MB_ ).LT.6 ) THEN
            INFO = -( 1000+MB_ )
         ELSE
            CALL CHK1MAT( N, 3, N, 3, 1, 1, DESCH, 7, INFO )
            IF( INFO.EQ.0 )
     $         CALL CHK1MAT( N, 3, N, 3, 1, 1, DESCZ, 11, INFO )
            IF( INFO.EQ.0 )
     $         CALL PCHK2MAT( N, 3, N, 3, 1, 1, DESCH, 7, N, 3, N, 3,
     $              1, 1, DESCZ, 11, 0, IWORK, IWORK, INFO )
         END IF
      END IF
*
*     Compute required workspace.
*
      CALL PSLAQR1( WANTT, WANTZ, N, ILO, IHI, H, DESCH, WR, WI,
     $     ILO, IHI, Z, DESCZ, WORK, -1, IWORK, -1, INFO )
      LWKOPT = WORK(1)
      LIWKOPT = IWORK(1)
      CALL PSLAQR0( WANTT, WANTZ, N, ILO, IHI, H, DESCH, WR, WI,
     $     ILO, IHI, Z, DESCZ, WORK, -1, IWORK, -1, INFO, 0 )
      IF( N.LT.NL ) THEN
         HROWS = NUMROC( NL, NB, MYROW, DESCH(RSRC_), NPROW )
         HCOLS = NUMROC( NL, NB, MYCOL, DESCH(CSRC_), NPCOL )
         WORK(1) = WORK(1) + FLOAT(2*HROWS*HCOLS)
      END IF
      LWKOPT = MAX( LWKOPT, WORK(1) )
      LIWKOPT = MAX( LIWKOPT, IWORK(1) )
      WORK(1) = LWKOPT
      IWORK(1) = LIWKOPT
*
      IF( .NOT.LQUERY .AND. LWORK.LT.INT(LWKOPT) ) THEN
         INFO = -13
      ELSEIF( .NOT.LQUERY .AND. LIWORK.LT.LIWKOPT ) THEN
         INFO = -15
      END IF
*
      IF( INFO.NE.0 ) THEN
*
*        Quick return in case of invalid argument.
*
         CALL PXERBLA( ICTXT, 'PSHSEQR', -INFO )
         RETURN
*
      ELSE IF( N.EQ.0 ) THEN
*
*        Quick return in case N = 0; nothing to do.
*
         RETURN
*
      ELSE IF( LQUERY ) THEN
*
*        Quick return in case of a workspace query.
*
         RETURN
*
      ELSE
*
*        Copy eigenvalues isolated by PSGEBAL.
*
         DO 10 I = 1, ILO - 1
            CALL INFOG2L( I, I, DESCH, NPROW, NPCOL, MYROW, MYCOL, II,
     $           JJ, HRSRC, HCSRC )
            IF( MYROW.EQ.HRSRC .AND. MYCOL.EQ.HCSRC ) THEN
               WR( I ) = H( (JJ-1)*LLDH + II )
            ELSE
               WR( I ) = ZERO
            END IF
            WI( I ) = ZERO
   10    CONTINUE
         IF( ILO.GT.1 )
     $      CALL SGSUM2D( ICTXT, 'All', '1-Tree', ILO-1, 1, WR, N, -1,
     $           -1 )
         DO 20 I = IHI + 1, N
            CALL INFOG2L( I, I, DESCH, NPROW, NPCOL, MYROW, MYCOL, II,
     $           JJ, HRSRC, HCSRC )
            IF( MYROW.EQ.HRSRC .AND. MYCOL.EQ.HCSRC ) THEN
               WR( I ) = H( (JJ-1)*LLDH + II )
            ELSE
               WR( I ) = ZERO
            END IF
            WI( I ) = ZERO
   20    CONTINUE
         IF( IHI.LT.N )
     $      CALL SGSUM2D( ICTXT, 'All', '1-Tree', N-IHI, 1, WR(IHI+1),
     $           N, -1, -1 )
*
*        Initialize Z, if requested.
*
         IF( INITZ )
     $      CALL PSLASET( 'A', N, N, ZERO, ONE, Z, 1, 1, DESCZ )
*
*        Quick return if possible.
*
         NPROCS = NPROW*NPCOL
         IF( ILO.EQ.IHI ) THEN
            CALL INFOG2L( ILO, ILO, DESCH, NPROW, NPCOL, MYROW,
     $           MYCOL, II, JJ, HRSRC, HCSRC )
            IF( MYROW.EQ.HRSRC .AND. MYCOL.EQ.HCSRC ) THEN
               WR( ILO ) = H( (JJ-1)*LLDH + II )
               IF( NPROCS.GT.1 )
     $            CALL SGEBS2D( ICTXT, 'All', '1-Tree', 1, 1, WR(ILO),
     $                 1 )
            ELSE
               CALL SGEBR2D( ICTXT, 'All', '1-Tree', 1, 1, WR(ILO),
     $              1, HRSRC, HCSRC )
            END IF
            WI( ILO ) = ZERO
            RETURN
         END IF
*
*        PSLAQR1/PSLAQR0 crossover point.
*
         NH = IHI-ILO+1
         NMIN = PILAENVX( ICTXT, 12, 'PSHSEQR',
     $        JOB( : 1 ) // COMPZ( : 1 ), N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
*
*        PSLAQR0 for big matrices; PSLAQR1 for small ones.
*
         IF( (.NOT. CRSOVER .AND. NH.GT.NTINY) .OR. NH.GT.NMIN .OR.
     $        DESCH(RSRC_).NE.0 .OR. DESCH(CSRC_).NE.0 ) THEN
            CALL PSLAQR0( WANTT, WANTZ, N, ILO, IHI, H, DESCH, WR, WI,
     $           ILO, IHI, Z, DESCZ, WORK, LWORK, IWORK, LIWORK, INFO,
     $           0 )
            IF( INFO.GT.0 .AND. ( DESCH(RSRC_).NE.0 .OR.
     $           DESCH(CSRC_).NE.0 ) ) THEN
*
*              A rare PSLAQR0 failure!  PSLAQR1 sometimes succeeds
*              when PSLAQR0 fails.
*
               KBOT = INFO
               CALL PSLAQR1( WANTT, WANTZ, N, ILO, IHI, H, DESCH, WR,
     $              WI, ILO, IHI, Z, DESCZ, WORK, LWORK, IWORK,
     $              LIWORK, INFO )
               INFO = -7777
            END IF
         ELSE
*
*           Small matrix.
*
            CALL PSLAQR1( WANTT, WANTZ, N, ILO, IHI, H, DESCH, WR, WI,
     $           ILO, IHI, Z, DESCZ, WORK, LWORK, IWORK, LIWORK, INFO )
*
            IF( INFO.GT.0 ) THEN
*
*              A rare PSLAQR1 failure!  PSLAQR0 sometimes succeeds
*              when PSLAQR1 fails.
*
               KBOT = INFO
*
               IF( N.GE.NL ) THEN
*
*                 Larger matrices have enough subdiagonal scratch
*                 space to call PSLAQR0 directly.
*
                  CALL PSLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, DESCH,
     $                 WR, WI, ILO, IHI, Z, DESCZ, WORK, LWORK,
     $                 IWORK, LIWORK, INFO, 0 )
               ELSE
*
*                 Tiny matrices don't have enough subdiagonal
*                 scratch space to benefit from PSLAQR0.  Hence,
*                 tiny matrices must be copied into a larger
*                 array before calling PSLAQR0.
*
                  HROWS = NUMROC( NL, NB, MYROW, DESCH(RSRC_), NPROW )
                  HCOLS = NUMROC( NL, NB, MYCOL, DESCH(CSRC_), NPCOL )
                  CALL DESCINIT( DESCH2, NL, NL, NB, NB, DESCH(RSRC_),
     $                 DESCH(CSRC_), ICTXT, MAX(1, HROWS), INFO )
                  CALL PSLACPY( 'All', N, N, H, 1, 1, DESCH, WORK, 1,
     $                 1, DESCH2 )
                  CALL PSELSET( WORK, N+1, N, DESCH2, ZERO )
                  CALL PSLASET( 'All', NL, NL-N, ZERO, ZERO, WORK, 1,
     $                 N+1, DESCH2 )
                  IPW = 1 + DESCH2(LLD_)*HCOLS
                  CALL PSLAQR0( WANTT, WANTZ, NL, ILO, KBOT, WORK,
     $                 DESCH2, WR, WI, ILO, IHI, Z, DESCZ,
     $                 WORK(IPW), LWORK-IPW+1, IWORK,
     $                 LIWORK, INFO, 0 )
                  IF( WANTT .OR. INFO.NE.0 )
     $               CALL PSLACPY( 'All', N, N, WORK, 1, 1, DESCH2,
     $                    H, 1, 1, DESCH )
               END IF
               INFO = -8888
            END IF
         END IF
*
*        Clear out the trash, if necessary.
*
         IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 )
     $      CALL PSLASET( 'L', N-2, N-2, ZERO, ZERO, H, 3, 1, DESCH )
*
*        Force any 2-by-2 blocks to be complex conjugate pairs of
*        eigenvalues by removing false such blocks.
*
         DO 30 I = ILO, IHI-1
            CALL PSELGET( 'All', ' ', TMP3, H, I+1, I, DESCH )
            IF( TMP3.NE.0.0E+00 ) THEN
               CALL PSELGET( 'All', ' ', TMP1, H, I, I, DESCH )
               CALL PSELGET( 'All', ' ', TMP2, H, I, I+1, DESCH )
               CALL PSELGET( 'All', ' ', TMP4, H, I+1, I+1, DESCH )
               CALL SLANV2( TMP1, TMP2, TMP3, TMP4, DUM1, DUM2, DUM3,
     $              DUM4, CS, SN )
               IF( TMP3.EQ.0.0E+00 ) THEN
                  IF( WANTT ) THEN
                     IF( I+2.LE.N )
     $                  CALL PSROT( N-I-1, H, I, I+2, DESCH,
     $                       DESCH(M_), H, I+1, I+2, DESCH, DESCH(M_),
     $                       CS, SN, WORK, LWORK, INFO )
                     CALL PSROT( I-1, H, 1, I, DESCH, 1, H, 1, I+1,
     $                    DESCH, 1, CS, SN, WORK, LWORK, INFO )
                  END IF
                  IF( WANTZ ) THEN
                     CALL PSROT( N, Z, 1, I, DESCZ, 1, Z, 1, I+1, DESCZ,
     $                    1, CS, SN, WORK, LWORK, INFO )
                  END IF
                  CALL PSELSET( H, I, I, DESCH, TMP1 )
                  CALL PSELSET( H, I, I+1, DESCH, TMP2 )
                  CALL PSELSET( H, I+1, I, DESCH, TMP3 )
                  CALL PSELSET( H, I+1, I+1, DESCH, TMP4 )
               END IF
            END IF
 30      CONTINUE
*
*        Read out eigenvalues: first let all the processes compute the
*        eigenvalue inside their diagonal blocks in parallel, except for
*        the eigenvalue located next to a block border. After that,
*        compute all eigenvalues located next to the block borders.
*        Finally, do a global summation over WR and WI so that all
*        processors receive the result.
*
         DO 40 K = ILO, IHI
            WR( K ) = ZERO
            WI( K ) = ZERO
 40      CONTINUE
         NB = DESCH( MB_ )
*
*        Loop 50: extract eigenvalues from the blocks which are not laid
*        out across a border of the processor mesh, except for those 1x1
*        blocks on the border.
*
         PAIR = .FALSE.
         DO 50 K = ILO, IHI
            IF( .NOT. PAIR ) THEN
               BORDER = MOD( K, NB ).EQ.0 .OR. ( K.NE.1 .AND.
     $              MOD( K, NB ).EQ.1 )
               IF( .NOT. BORDER ) THEN
                  CALL INFOG2L( K, K, DESCH, NPROW, NPCOL, MYROW,
     $                 MYCOL, ILOC1, JLOC1, HRSRC1, HCSRC1 )
                  IF( MYROW.EQ.HRSRC1 .AND. MYCOL.EQ.HCSRC1 ) THEN
                     ELEM1 = H((JLOC1-1)*LLDH+ILOC1)
                     IF( K.LT.N ) THEN
                        ELEM3( 1 ) = H((JLOC1-1)*LLDH+ILOC1+1)
                     ELSE
                        ELEM3( 1 ) = ZERO
                     END IF
                     IF( ELEM3( 1 ).NE.ZERO ) THEN
                        ELEM2( 1 ) = H((JLOC1)*LLDH+ILOC1)
                        ELEM4 = H((JLOC1)*LLDH+ILOC1+1)
                        CALL SLANV2( ELEM1, ELEM2( 1 ), ELEM3( 1 ),
     $                       ELEM4, WR( K ), WI( K ), WR( K+1 ),
     $                       WI( K+1 ), SN, CS )
                        PAIR = .TRUE.
                     ELSE
                        IF( K.GT.1 ) THEN
                           TMP = H((JLOC1-2)*LLDH+ILOC1)
                           IF( TMP.NE.ZERO ) THEN
                              ELEM1 = H((JLOC1-2)*LLDH+ILOC1-1)
                              ELEM2( 1 ) = H((JLOC1-1)*LLDH+ILOC1-1)
                              ELEM3( 1 ) = H((JLOC1-2)*LLDH+ILOC1)
                              ELEM4 = H((JLOC1-1)*LLDH+ILOC1)
                              CALL SLANV2( ELEM1, ELEM2( 1 ),
     $                             ELEM3( 1 ), ELEM4, WR( K-1 ),
     $                             WI( K-1 ), WR( K ), WI( K ), SN, CS )
                           ELSE
                              WR( K ) = ELEM1
                           END IF
                        ELSE
                           WR( K ) = ELEM1
                        END IF
                     END IF
                  END IF
               END IF
            ELSE
               PAIR = .FALSE.
            END IF
 50      CONTINUE
*
*        Loop 60: extract eigenvalues from the blocks which are laid
*        out across a border of the processor mesh. The processors are
*        numbered as below:
*
*                        1 | 2
*                        --+--
*                        3 | 4
*
         DO 60 K = ICEIL(ILO,NB)*NB, IHI-1, NB
            CALL INFOG2L( K, K, DESCH, NPROW, NPCOL, MYROW, MYCOL,
     $           ILOC1, JLOC1, HRSRC1, HCSRC1 )
            CALL INFOG2L( K, K+1, DESCH, NPROW, NPCOL, MYROW, MYCOL,
     $           ILOC2, JLOC2, HRSRC2, HCSRC2 )
            CALL INFOG2L( K+1, K, DESCH, NPROW, NPCOL, MYROW, MYCOL,
     $           ILOC3, JLOC3, HRSRC3, HCSRC3 )
            CALL INFOG2L( K+1, K+1, DESCH, NPROW, NPCOL, MYROW, MYCOL,
     $           ILOC4, JLOC4, HRSRC4, HCSRC4 )
            IF( MYROW.EQ.HRSRC2 .AND. MYCOL.EQ.HCSRC2 ) THEN
               ELEM2( 1 ) = H((JLOC2-1)*LLDH+ILOC2)
               IF( HRSRC1.NE.HRSRC2 .OR. HCSRC1.NE.HCSRC2 )
     $            CALL SGESD2D( ICTXT, 1, 1, ELEM2, 1, HRSRC1, HCSRC1)
            END IF
            IF( MYROW.EQ.HRSRC3 .AND. MYCOL.EQ.HCSRC3 ) THEN
               ELEM3( 1 ) = H((JLOC3-1)*LLDH+ILOC3)
               IF( HRSRC1.NE.HRSRC3 .OR. HCSRC1.NE.HCSRC3 )
     $            CALL SGESD2D( ICTXT, 1, 1, ELEM3, 1, HRSRC1, HCSRC1)
            END IF
            IF( MYROW.EQ.HRSRC4 .AND. MYCOL.EQ.HCSRC4 ) THEN
               WORK(1) = H((JLOC4-1)*LLDH+ILOC4)
               IF( K+1.LT.N ) THEN
                  WORK(2) = H((JLOC4-1)*LLDH+ILOC4+1)
               ELSE
                  WORK(2) = ZERO
               END IF
               IF( HRSRC1.NE.HRSRC4 .OR. HCSRC1.NE.HCSRC4 )
     $            CALL SGESD2D( ICTXT, 2, 1, WORK, 2, HRSRC1, HCSRC1 )
            END IF
            IF( MYROW.EQ.HRSRC1 .AND. MYCOL.EQ.HCSRC1 ) THEN
               ELEM1 = H((JLOC1-1)*LLDH+ILOC1)
               IF( HRSRC1.NE.HRSRC2 .OR. HCSRC1.NE.HCSRC2 )
     $            CALL SGERV2D( ICTXT, 1, 1, ELEM2, 1, HRSRC2, HCSRC2)
               IF( HRSRC1.NE.HRSRC3 .OR. HCSRC1.NE.HCSRC3 )
     $            CALL SGERV2D( ICTXT, 1, 1, ELEM3, 1, HRSRC3, HCSRC3)
               IF( HRSRC1.NE.HRSRC4 .OR. HCSRC1.NE.HCSRC4 )
     $            CALL SGERV2D( ICTXT, 2, 1, WORK, 2, HRSRC4, HCSRC4 )
               ELEM4 = WORK(1)
               ELEM5 = WORK(2)
               IF( ELEM5.EQ.ZERO ) THEN
                  IF( WR( K ).EQ.ZERO .AND. WI( K ).EQ.ZERO ) THEN
                     CALL SLANV2( ELEM1, ELEM2( 1 ), ELEM3( 1 ), ELEM4,
     $                    WR( K ), WI( K ), WR( K+1 ), WI( K+1 ), SN,
     $                    CS )
                  ELSEIF( WR( K+1 ).EQ.ZERO .AND. WI( K+1 ).EQ.ZERO )
     $                 THEN
                     WR( K+1 ) = ELEM4
                  END IF
               ELSEIF( WR( K ).EQ.ZERO .AND. WI( K ).EQ.ZERO )
     $              THEN
                  WR( K ) = ELEM1
               END IF
            END IF
 60      CONTINUE
*
         IF( NPROCS.GT.1 ) THEN
            CALL SGSUM2D( ICTXT, 'All', ' ', IHI-ILO+1, 1, WR(ILO), N,
     $           -1, -1 )
            CALL SGSUM2D( ICTXT, 'All', ' ', IHI-ILO+1, 1, WI(ILO), N,
     $           -1, -1 )
         END IF
*
      END IF
*
      WORK(1) = LWKOPT
      IWORK(1) = LIWKOPT
      RETURN
*
*     End of PSHSEQR
*
      END
