      SUBROUTINE PSLAED3( ICTXT, K, N, NB, D, DROW, DCOL, RHO, DLAMDA,
     $                    W, Z, U, LDU, BUF, INDX, INDCOL, INDROW,
     $                    INDXR, INDXC, CTOT, NPCOL, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      INTEGER            DCOL, DROW, ICTXT, INFO, K, LDU, N, NB, NPCOL
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            CTOT( 0: NPCOL-1, 4 ), INDCOL( * ),
     $                   INDROW( * ), INDX( * ), INDXC( * ), INDXR( * )
      REAL               BUF( * ), D( * ), DLAMDA( * ), U( LDU, * ),
     $                   W( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PSLAED3 finds the roots of the secular equation, as defined by the
*  values in D, W, and RHO, between 1 and K.  It makes the
*  appropriate calls to SLAED4
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  ICTXT  (global input) INTEGER
*         The BLACS context handle, indicating the global context of
*         the operation on the matrix. The context itself is global.
*
*  K      (output) INTEGER
*         The number of non-deflated eigenvalues, and the order of the
*         related secular equation. 0 <= K <=N.
*
*  N      (input) INTEGER
*         The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  NB      (global input) INTEGER
*          The blocking factor used to distribute the columns of the
*          matrix. NB >= 1.
*
*  D      (input/output) REAL array, dimension (N)
*         On entry, D contains the eigenvalues of the two submatrices to
*         be combined.
*         On exit, D contains the trailing (N-K) updated eigenvalues
*         (those which were deflated) sorted into increasing order.
*
*  DROW   (global input) INTEGER
*          The process row over which the first row of the matrix D is
*          distributed. 0 <= DROW < NPROW.
*
*  DCOL   (global input) INTEGER
*          The process column over which the first column of the
*          matrix D is distributed. 0 <= DCOL < NPCOL.
*
*  RHO    (global input/output) REAL
*         On entry, the off-diagonal element associated with the rank-1
*         cut which originally split the two submatrices which are now
*         being recombined.
*         On exit, RHO has been modified to the value required by
*         PSLAED3.
*
*  DLAMDA (global output) REAL array, dimension (N)
*         A copy of the first K eigenvalues which will be used by
*         SLAED4 to form the secular equation.
*
*  W      (global output) REAL array, dimension (N)
*         The first k values of the final deflation-altered z-vector
*         which will be passed to SLAED4.
*
*  Z      (global input) REAL array, dimension (N)
*         On entry, Z contains the updating vector (the last
*         row of the first sub-eigenvector matrix and the first row of
*         the second sub-eigenvector matrix).
*         On exit, the contents of Z have been destroyed by the updating
*         process.
*
*  U     (global output) REAL array
*         global dimension (N, N), local dimension (LDU, NQ).
*         (See PSLAED0 for definition of NQ.)
*         Q  contains the orthonormal eigenvectors of the symmetric
*         tridiagonal matrix.
*
*  LDU    (input) INTEGER
*         The leading dimension of the array U.
*
*  BUF    (workspace) REAL array, dimension 3*N
*
*
*  INDX   (workspace) INTEGER array, dimension (N)
*         The permutation used to sort the contents of DLAMDA into
*         ascending order.
*
*  INDCOL (workspace) INTEGER array, dimension (N)
*
*
*  INDROW (workspace) INTEGER array, dimension (N)
*
*
*  INDXR (workspace) INTEGER array, dimension (N)
*
*
*  INDXC (workspace) INTEGER array, dimension (N)
*
*  CTOT   (workspace) INTEGER array, dimension( NPCOL, 4)
*
*  NPCOL   (global input) INTEGER
*          The total number of columns over which the distributed
*           submatrix is distributed.
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  The algorithm failed to compute the ith eigenvalue.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            COL, GI, I, IINFO, IIU, IPD, IU, J, JJU, JU,
     $                   KK, KL, KLC, KLR, MYCOL, MYKL, MYKLR, MYROW,
     $                   NPROW, PDC, PDR, ROW
      REAL               AUX, TEMP
*     ..
*     .. External Functions ..
      INTEGER            INDXG2L
      REAL               SLAMC3, SNRM2
      EXTERNAL           INDXG2L, SLAMC3, SNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, SCOPY, SGEBR2D, SGEBS2D,
     $                   SGERV2D, SGESD2D, SLAED4
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      IINFO = 0
*
*     Quick return if possible
*
      IF( K.EQ.0 )
     $   RETURN
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      ROW = DROW
      COL = DCOL
      DO 20 I = 1, N, NB
         DO 10 J = 0, NB - 1
            INDROW( I+J ) = ROW
            INDCOL( I+J ) = COL
   10    CONTINUE
         ROW = MOD( ROW+1, NPROW )
         COL = MOD( COL+1, NPCOL )
   20 CONTINUE
*
      MYKL = CTOT( MYCOL, 1 ) + CTOT( MYCOL, 2 ) + CTOT( MYCOL, 3 )
      KLR = MYKL / NPROW
      IF( MYROW.EQ.DROW ) THEN
         MYKLR = KLR + MOD( MYKL, NPROW )
      ELSE
         MYKLR = KLR
      END IF
      PDC = 1
      COL = DCOL
   30 CONTINUE
      IF( MYCOL.NE.COL ) THEN
         PDC = PDC + CTOT( COL, 1 ) + CTOT( COL, 2 ) + CTOT( COL, 3 )
         COL = MOD( COL+1, NPCOL )
         GO TO 30
      END IF
      PDR = PDC
      KL = KLR + MOD( MYKL, NPROW )
      ROW = DROW
   40 CONTINUE
      IF( MYROW.NE.ROW ) THEN
         PDR = PDR + KL
         KL = KLR
         ROW = MOD( ROW+1, NPROW )
         GO TO 40
      END IF
*
      DO 50 I = 1, K
         DLAMDA( I ) = SLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
         Z( I ) = ONE
   50 CONTINUE
      IF( MYKLR.GT.0 ) THEN
         KK = PDR
         DO 80 I = 1, MYKLR
            CALL SLAED4( K, KK, DLAMDA, W, BUF, RHO, BUF( K+I ), IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = KK
            END IF
*
*     ..Compute part of z
*
            DO 60 J = 1, KK - 1
               Z( J ) = Z( J )*( BUF( J ) /
     $                  ( DLAMDA( J )-DLAMDA( KK ) ) )
   60       CONTINUE
            Z( KK ) = Z( KK )*BUF( KK )
            DO 70 J = KK + 1, K
               Z( J ) = Z( J )*( BUF( J ) /
     $                  ( DLAMDA( J )-DLAMDA( KK ) ) )
   70       CONTINUE
            KK = KK + 1
   80    CONTINUE
*
         IF( MYROW.NE.DROW ) THEN
            CALL SCOPY( K, Z, 1, BUF, 1 )
            CALL SGESD2D( ICTXT, K+MYKLR, 1, BUF, K+MYKLR, DROW, MYCOL )
         ELSE
            IPD = 2*K + 1
            CALL SCOPY( MYKLR, BUF( K+1 ), 1, BUF( IPD ), 1 )
            IF( KLR.GT.0 ) THEN
               IPD = MYKLR + IPD
               ROW = MOD( DROW+1, NPROW )
               DO 100 I = 1, NPROW - 1
                  CALL SGERV2D( ICTXT, K+KLR, 1, BUF, K+KLR, ROW,
     $                          MYCOL )
                  CALL SCOPY( KLR, BUF( K+1 ), 1, BUF( IPD ), 1 )
                  DO 90 J = 1, K
                     Z( J ) = Z( J )*BUF( J )
   90             CONTINUE
                  IPD = IPD + KLR
                  ROW = MOD( ROW+1, NPROW )
  100          CONTINUE
            END IF
         END IF
      END IF
*
      IF( MYROW.EQ.DROW ) THEN
         IF( MYCOL.NE.DCOL .AND. MYKL.NE.0 ) THEN
            CALL SCOPY( K, Z, 1, BUF, 1 )
            CALL SCOPY( MYKL, BUF( 2*K+1 ), 1, BUF( K+1 ), 1 )
            CALL SGESD2D( ICTXT, K+MYKL, 1, BUF, K+MYKL, MYROW, DCOL )
         ELSE IF( MYCOL.EQ.DCOL ) THEN
            IPD = 2*K + 1
            COL = DCOL
            KL = MYKL
            DO 120 I = 1, NPCOL - 1
               IPD = IPD + KL
               COL = MOD( COL+1, NPCOL )
               KL = CTOT( COL, 1 ) + CTOT( COL, 2 ) + CTOT( COL, 3 )
               IF( KL.NE.0 ) THEN
                  CALL SGERV2D( ICTXT, K+KL, 1, BUF, K+KL, MYROW, COL )
                  CALL SCOPY( KL, BUF( K+1 ), 1, BUF( IPD ), 1 )
                  DO 110 J = 1, K
                     Z( J ) = Z( J )*BUF( J )
  110             CONTINUE
               END IF
  120       CONTINUE
            DO 130 I = 1, K
               Z( I ) = SIGN( SQRT( -Z( I ) ), W( I ) )
  130       CONTINUE
*
         END IF
      END IF
*
*     Diffusion
*
      IF( MYROW.EQ.DROW .AND. MYCOL.EQ.DCOL ) THEN
         CALL SCOPY( K, Z, 1, BUF, 1 )
         CALL SCOPY( K, BUF( 2*K+1 ), 1, BUF( K+1 ), 1 )
         CALL SGEBS2D( ICTXT, 'All', ' ', 2*K, 1, BUF, 2*K )
      ELSE
         CALL SGEBR2D( ICTXT, 'All', ' ', 2*K, 1, BUF, 2*K, DROW, DCOL )
         CALL SCOPY( K, BUF, 1, Z, 1 )
      END IF
*
*     Copy of D at the good place
*
      KLC = 0
      KLR = 0
      DO 140 I = 1, K
         GI = INDX( I )
         D( GI ) = BUF( K+I )
         COL = INDCOL( GI )
         ROW = INDROW( GI )
         IF( COL.EQ.MYCOL ) THEN
            KLC = KLC + 1
            INDXC( KLC ) = I
         END IF
         IF( ROW.EQ.MYROW ) THEN
            KLR = KLR + 1
            INDXR( KLR ) = I
         END IF
  140 CONTINUE
*
*     Compute eigenvectors of the modified rank-1 modification.
*
      IF( MYKL.NE.0 ) THEN
         DO 180 J = 1, MYKL
            KK = INDXC( J )
            JU = INDX( KK )
            JJU = INDXG2L( JU, NB, J, J, NPCOL )
            CALL SLAED4( K, KK, DLAMDA, W, BUF, RHO, AUX, IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = KK
            END IF
            IF( K.EQ.1 .OR. K.EQ.2 ) THEN
               DO 150 I = 1, KLR
                  KK = INDXR( I )
                  IU = INDX( KK )
                  IIU = INDXG2L( IU, NB, J, J, NPROW )
                  U( IIU, JJU ) = BUF( KK )
  150          CONTINUE
               GO TO 180
            END IF
*
            DO 160 I = 1, K
               BUF( I ) = Z( I ) / BUF( I )
  160       CONTINUE
            TEMP = SNRM2( K, BUF, 1 )
            DO 170 I = 1, KLR
               KK = INDXR( I )
               IU = INDX( KK )
               IIU = INDXG2L( IU, NB, J, J, NPROW )
               U( IIU, JJU ) = BUF( KK ) / TEMP
  170       CONTINUE
  180    CONTINUE
      END IF
*
  190 CONTINUE
      RETURN
*
*     End of PSLAED3
*
      END
