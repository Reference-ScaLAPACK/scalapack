      SUBROUTINE PZGEBDRV( M, N, A, IA, JA, DESCA, D, E, TAUQ, TAUP,
     $                     WORK, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER              INFO, IA, JA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER              DESCA( * )
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( * ), TAUP( * ), TAUQ( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZGEBDRV computes sub( A ) = A(IA:IA+M-1,JA:JA+N-1) from sub( A ),
*  Q, P returned by PZGEBRD:
*
*                         sub( A ) := Q * B * P'.
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
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of sub( A )
*          as returned by PZGEBRD. On exit, the original distribu-
*          ted matrix sub( A ) is restored.
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  D       (local input) DOUBLE PRECISION array, dimension
*          LOCc(JA+MIN(M,N)-1) if M >= N; LOCr(IA+MIN(M,N)-1) otherwise.
*          The distributed diagonal elements of the bidiagonal matrix
*          B: D(i) = A(i,i). D is tied to the distributed matrix A.
*
*  E       (local input) DOUBLE PRECISION array, dimension
*          LOCr(IA+MIN(M,N)-1) if M >= N; LOCc(JA+MIN(M,N)-2) otherwise.
*          The distributed off-diagonal elements of the bidiagonal
*          distributed matrix B:
*          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
*          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
*          E is tied to the distributed matrix A.
*
*  TAUQ    (local input) COMPLEX*16 array dimension
*          LOCc(JA+MIN(M,N)-1). The scalar factors of the elementary
*          reflectors which represent the unitary matrix Q. TAUQ is
*          tied to the distributed matrix A. See Further Details.
*
*  TAUP    (local input) COMPLEX*16 array, dimension
*          LOCr(IA+MIN(M,N)-1). The scalar factors of the elementary
*          reflectors which represent the unitary matrix P. TAUP is
*          tied to the distributed matrix A. See Further Details.
*
*  WORK    (local workspace) COMPLEX*16 array, dimension (LWORK)
*          LWORK >= 2*NB*( MP + NQ + NB )
*
*          where NB = MB_A = NB_A,
*          IROFFA = MOD( IA-1, NB ), ICOFFA = MOD( JA-1, NB ),
*          IAROW = INDXG2P( IA, NB, MYROW, RSRC_A, NPROW ),
*          IACOL = INDXG2P( JA, NB, MYCOL, CSRC_A, NPCOL ),
*          MP = NUMROC( M+IROFFA, NB, MYROW, IAROW, NPROW ),
*          NQ = NUMROC( N+ICOFFA, NB, MYCOL, IACOL, NPCOL ).
*
*          INDXG2P and NUMROC are ScaLAPACK tool functions;
*          MYROW, MYCOL, NPROW and NPCOL can be determined by calling
*          the subroutine BLACS_GRIDINFO.
*
*  INFO    (global output) INTEGER
*          On exit, if INFO <> 0, a discrepancy has been found between
*          the diagonal and off-diagonal elements of A and the copies
*          contained in the arrays D and E.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   REIGHT, RZERO
      PARAMETER          ( REIGHT = 8.0D+0, RZERO = 0.0D+0 )
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IACOL, IAROW, ICTXT, IIA, IL, IPTP, IPTQ,
     $                   IPV, IPW, IPWK, IOFF, IV, J, JB, JJA, JL, JV,
     $                   K, MN, MP, MYCOL, MYROW, NB, NPCOL, NPROW, NQ
      DOUBLE PRECISION   ADDBND, D2, E2
      COMPLEX*16         D1, E1
*     ..
*     .. Local Arrays ..
      INTEGER            DESCD( DLEN_ ), DESCE( DLEN_ ), DESCV( DLEN_ ),
     $                   DESCW( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DESCSET, IGSUM2D, INFOG2L,
     $                   PDELGET, PZLACPY, PZLARFB, PZLARFT,
     $                   PZLASET, PZELGET
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P, NUMROC
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           INDXG2P, NUMROC, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DCMPLX, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      INFO = 0
      NB = DESCA( MB_ )
      IOFF = MOD( IA-1, NB )
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,
     $              IAROW, IACOL )
      MP = NUMROC( M+IOFF, NB, MYROW, IAROW, NPROW )
      NQ = NUMROC( N+IOFF, NB, MYCOL, IACOL, NPCOL )
      IPV  = 1
      IPW  = IPV + MP*NB
      IPTP = IPW + NQ*NB
      IPTQ = IPTP + NB*NB
      IPWK = IPTQ + NB*NB
*
      IV = 1
      JV = 1
      MN = MIN( M, N )
      IL = MAX( ( (IA+MN-2) / NB )*NB + 1, IA )
      JL = MAX( ( (JA+MN-2) / NB )*NB + 1, JA )
      IAROW = INDXG2P( IL, NB, MYROW, DESCA( RSRC_ ), NPROW )
      IACOL = INDXG2P( JL, NB, MYCOL, DESCA( CSRC_ ), NPCOL )
      CALL DESCSET( DESCV, IA+M-IL, NB, NB, NB, IAROW, IACOL, ICTXT,
     $              MAX( 1, MP ) )
      CALL DESCSET( DESCW, NB, JA+N-JL, NB, NB, IAROW, IACOL, ICTXT,
     $              NB )
*
      ADDBND = REIGHT * PDLAMCH( ICTXT, 'eps' )
*
*     When A is an upper bidiagonal form
*
      IF( M.GE.N ) THEN
*
         CALL DESCSET( DESCD, 1, JA+MN-1, 1, DESCA( NB_ ), MYROW,
     $                 DESCA( CSRC_ ), DESCA( CTXT_ ), 1 )
         CALL DESCSET( DESCE, IA+MN-1, 1, DESCA( MB_ ), 1,
     $                 DESCA( RSRC_ ), MYCOL, DESCA( CTXT_ ),
     $                 DESCA( LLD_ ) )
*
         DO 10 J = 0, MN-1
            D1 = ZERO
            E1 = ZERO
            D2 = RZERO
            E2 = RZERO
            CALL PDELGET( ' ', ' ', D2, D, 1, JA+J, DESCD )
            CALL PZELGET( 'Columnwise', ' ', D1, A, IA+J, JA+J, DESCA )
            IF( J.LT.(MN-1) ) THEN
               CALL PDELGET( ' ', ' ', E2, E, IA+J, 1, DESCE )
               CALL PZELGET( 'Rowwise', ' ', E1, A, IA+J, JA+J+1,
     $                       DESCA )
            END IF
*
            IF( ( ABS( D1-DCMPLX( D2 ) ).GT.( ABS( D2 )*ADDBND ) ) .OR.
     $          ( ABS( E1-DCMPLX( E2 ) ).GT.( ABS( E2 )*ADDBND ) ) )
     $         INFO = INFO + 1
   10    CONTINUE
*
         DO 20 J = JL, JA+NB-IOFF, -NB
            JB = MIN( JA+N-J, NB )
            I  = IA + J - JA
            K  = I - IA + 1
*
*           Compute upper triangular matrix TQ from TAUQ.
*
            CALL PZLARFT( 'Forward', 'Columnwise', M-K+1, JB, A, I, J,
     $                    DESCA, TAUQ, WORK( IPTQ ), WORK( IPWK ) )
*
*           Copy Householder vectors into workspace.
*
            CALL PZLACPY( 'Lower', M-K+1, JB, A, I, J, DESCA,
     $                    WORK( IPV ), IV, JV, DESCV )
            CALL PZLASET( 'Upper', M-K+1, JB, ZERO, ONE, WORK( IPV ),
     $                    IV, JV, DESCV )
*
*           Zero out the strict lower triangular part of A.
*
            CALL PZLASET( 'Lower', M-K, JB, ZERO, ZERO, A, I+1, J,
     $                    DESCA )
*
*           Compute upper triangular matrix TP from TAUP.
*
            CALL PZLARFT( 'Forward', 'Rowwise', N-K, JB, A, I, J+1,
     $                    DESCA, TAUP, WORK( IPTP ), WORK( IPWK ) )
*
*           Copy Householder vectors into workspace.
*
            CALL PZLACPY( 'Upper', JB, N-K, A, I, J+1, DESCA,
     $                    WORK( IPW ), IV, JV+1, DESCW )
            CALL PZLASET( 'Lower', JB, N-K, ZERO, ONE, WORK( IPW ), IV,
     $                    JV+1, DESCW )
*
*           Zero out the strict+1 upper triangular part of A.
*
            CALL PZLASET( 'Upper', JB, N-K-1, ZERO, ZERO, A, I, J+2,
     $                    DESCA )
*
*           Apply block Householder transformation from Left.
*
            CALL PZLARFB( 'Left', 'No transpose', 'Forward',
     $                    'Columnwise', M-K+1, N-K+1, JB, WORK( IPV ),
     $                    IV, JV, DESCV, WORK( IPTQ ), A, I, J, DESCA,
     $                    WORK( IPWK ) )
*
*           Apply block Householder transformation from Right.
*
            CALL PZLARFB( 'Right', 'Conjugate transpose', 'Forward',
     $                    'Rowwise', M-K+1, N-K, JB, WORK( IPW ), IV,
     $                    JV+1, DESCW, WORK( IPTP ), A, I, J+1, DESCA,
     $                    WORK( IPWK ) )
*
            DESCV( M_ ) = DESCV( M_ ) + NB
            DESCV( RSRC_ ) = MOD( DESCV( RSRC_ ) + NPROW - 1, NPROW )
            DESCV( CSRC_ ) = MOD( DESCV( CSRC_ ) + NPCOL - 1, NPCOL )
            DESCW( N_ ) = DESCW( N_ ) + NB
            DESCW( RSRC_ ) = DESCV( RSRC_ )
            DESCW( CSRC_ ) = DESCV( CSRC_ )
*
   20    CONTINUE
*
*        Handle first block separately
*
         JB = MIN( N, NB - IOFF )
         IV = IOFF + 1
         JV = IOFF + 1
*
*        Compute upper triangular matrix TQ from TAUQ.
*
         CALL PZLARFT( 'Forward', 'Columnwise', M, JB, A, IA, JA, DESCA,
     $                 TAUQ, WORK( IPTQ ), WORK( IPWK ) )
*
*        Copy Householder vectors into workspace.
*
         CALL PZLACPY( 'Lower', M, JB, A, IA, JA, DESCA, WORK( IPV ),
     $                 IV, JV, DESCV )
         CALL PZLASET( 'Upper', M, JB, ZERO, ONE, WORK( IPV ), IV, JV,
     $                 DESCV )
*
*        Zero out the strict lower triangular part of A.
*
         CALL PZLASET( 'Lower', M-1, JB, ZERO, ZERO, A, IA+1, JA,
     $                 DESCA )
*
*        Compute upper triangular matrix TP from TAUP.
*
         CALL PZLARFT( 'Forward', 'Rowwise', N-1, JB, A, IA, JA+1,
     $                 DESCA, TAUP, WORK( IPTP ), WORK( IPWK ) )
*
*        Copy Householder vectors into workspace.
*
         CALL PZLACPY( 'Upper', JB, N-1, A, IA, JA+1, DESCA,
     $                 WORK( IPW ), IV, JV+1, DESCW )
         CALL PZLASET( 'Lower', JB, N-1, ZERO, ONE, WORK( IPW ), IV,
     $                 JV+1, DESCW )
*
*        Zero out the strict+1 upper triangular part of A.
*
         CALL PZLASET( 'Upper', JB, N-2, ZERO, ZERO, A, IA, JA+2,
     $                 DESCA )
*
*        Apply block Householder transformation from left.
*
         CALL PZLARFB( 'Left', 'No transpose', 'Forward', 'Columnwise',
     $                 M, N, JB, WORK( IPV ), IV, JV, DESCV,
     $                 WORK( IPTQ ), A, IA, JA, DESCA, WORK( IPWK ) )
*
*        Apply block Householder transformation from right.
*
         CALL PZLARFB( 'Right', 'Conjugate transpose', 'Forward',
     $                 'Rowwise', M, N-1, JB, WORK( IPW ), IV, JV+1,
     $                 DESCW, WORK( IPTP ), A, IA, JA+1, DESCA,
     $                 WORK( IPWK ) )
*
      ELSE
*
         CALL DESCSET( DESCD, IA+MN-1, 1, DESCA( MB_ ), 1,
     $                 DESCA( RSRC_ ), MYCOL, DESCA( CTXT_ ),
     $                 DESCA( LLD_ ) )
         CALL DESCSET( DESCE, 1, JA+MN-2, 1, DESCA( NB_ ), MYROW,
     $                 DESCA( CSRC_ ), DESCA( CTXT_ ), 1 )
*
         DO 30 J = 0, MN-1
            D1 = ZERO
            E1 = ZERO
            D2 = RZERO
            E2 = RZERO
            CALL PDELGET( ' ', ' ', D2, D, IA+J, 1, DESCD )
            CALL PZELGET( 'Rowwise', ' ', D1, A, IA+J, JA+J, DESCA )
            IF( J.LT.(MN-1) ) THEN
               CALL PDELGET( ' ', ' ', E2, E, 1, JA+J, DESCE )
               CALL PZELGET( 'Columnwise', ' ', E1, A, IA+J+1, JA+J,
     $                       DESCA )
            END IF
*
            IF( ( ABS( D1-DCMPLX( D2 ) ).GT.( ABS( D2 )*ADDBND ) ) .OR.
     $          ( ABS( E1-DCMPLX( E2 ) ).GT.( ABS( E2 )*ADDBND ) ) )
     $         INFO = INFO + 1
   30    CONTINUE
*
         DO 40 I = IL, IA+NB-IOFF, -NB
            JB = MIN( IA+M-I, NB )
            J  = JA + I - IA
            K  = J - JA + 1
*
*           Compute upper triangular matrix TQ from TAUQ.
*
            CALL PZLARFT( 'Forward', 'Columnwise', M-K, JB, A, I+1, J,
     $                    DESCA, TAUQ, WORK( IPTQ ), WORK( IPWK ) )
*
*           Copy Householder vectors into workspace.
*
            CALL PZLACPY( 'Lower', M-K, JB, A, I+1, J, DESCA,
     $                    WORK( IPV ), IV+1, JV, DESCV )
            CALL PZLASET( 'Upper', M-K, JB, ZERO, ONE, WORK( IPV ),
     $                    IV+1, JV, DESCV )
*
*           Zero out the strict lower triangular part of A.
*
            CALL PZLASET( 'Lower', M-K-1, JB, ZERO, ZERO, A, I+2, J,
     $                    DESCA )
*
*           Compute upper triangular matrix TP from TAUP.
*
            CALL PZLARFT( 'Forward', 'Rowwise', N-K+1, JB, A, I, J,
     $                    DESCA, TAUP, WORK( IPTP ), WORK( IPWK ) )
*
*           Copy Householder vectors into workspace.
*
            CALL PZLACPY( 'Upper', JB, N-K+1, A, I, J, DESCA,
     $                    WORK( IPW ), IV, JV, DESCW )
            CALL PZLASET( 'Lower', JB, N-K+1, ZERO, ONE, WORK( IPW ),
     $                    IV, JV, DESCW )
*
*           Zero out the strict+1 upper triangular part of A.
*
            CALL PZLASET( 'Upper', JB, N-K, ZERO, ZERO, A, I, J+1,
     $                    DESCA )
*
*           Apply block Householder transformation from Left.
*
            CALL PZLARFB( 'Left', 'No transpose', 'Forward',
     $                    'Columnwise', M-K, N-K+1, JB, WORK( IPV ),
     $                    IV+1, JV, DESCV, WORK( IPTQ ), A, I+1, J,
     $                    DESCA, WORK( IPWK ) )
*
*           Apply block Householder transformation from Right.
*
            CALL PZLARFB( 'Right', 'Conjugate transpose', 'Forward',
     $                    'Rowwise', M-K+1, N-K+1, JB, WORK( IPW ), IV,
     $                    JV, DESCW, WORK( IPTP ), A, I, J, DESCA,
     $                    WORK( IPWK ) )
*
            DESCV( M_ ) = DESCV( M_ ) + NB
            DESCV( RSRC_ ) = MOD( DESCV( RSRC_ ) + NPROW - 1, NPROW )
            DESCV( CSRC_ ) = MOD( DESCV( CSRC_ ) + NPCOL - 1, NPCOL )
            DESCW( N_ ) = DESCW( N_ ) + NB
            DESCW( RSRC_ ) = DESCV( RSRC_ )
            DESCW( CSRC_ ) = DESCV( CSRC_ )
*
   40    CONTINUE
*
*        Handle first block separately
*
         JB = MIN( M, NB - IOFF )
         IV = IOFF + 1
         JV = IOFF + 1
*
*        Compute upper triangular matrix TQ from TAUQ.
*
         CALL PZLARFT( 'Forward', 'Columnwise', M-1, JB, A, IA+1, JA,
     $                 DESCA, TAUQ, WORK( IPTQ ), WORK( IPWK ) )
*
*        Copy Householder vectors into workspace.
*
         CALL PZLACPY( 'Lower', M-1, JB, A, IA+1, JA, DESCA,
     $                 WORK( IPV ), IV+1, JV, DESCV )
         CALL PZLASET( 'Upper', M-1, JB, ZERO, ONE, WORK( IPV ), IV+1,
     $                 JV, DESCV )
*
*        Zero out the strict lower triangular part of A.
*
         CALL PZLASET( 'Lower', M-2, JB, ZERO, ZERO, A, IA+2, JA,
     $                 DESCA )
*
*        Compute upper triangular matrix TP from TAUP.
*
         CALL PZLARFT( 'Forward', 'Rowwise', N, JB, A, IA, JA, DESCA,
     $                 TAUP, WORK( IPTP ), WORK( IPWK ) )
*
*        Copy Householder vectors into workspace.
*
         CALL PZLACPY( 'Upper', JB, N, A, IA, JA, DESCA, WORK( IPW ),
     $                 IV, JV, DESCW )
         CALL PZLASET( 'Lower', JB, N, ZERO, ONE, WORK( IPW ), IV, JV,
     $                 DESCW )
*
*        Zero out the strict+1 upper triangular part of A.
*
         CALL PZLASET( 'Upper', JB, N-1, ZERO, ZERO, A, IA, JA+1,
     $                 DESCA )
*
*        Apply block Householder transformation from left
*
         CALL PZLARFB( 'Left', 'No transpose', 'Forward', 'Columnwise',
     $                 M-1, N, JB, WORK( IPV ), IV+1, JV, DESCV,
     $                 WORK( IPTQ ), A, IA+1, JA, DESCA, WORK( IPWK ) )
*
*        Apply block Householder transformation from right
*
         CALL PZLARFB( 'Right', 'Conjugate transpose', 'Forward',
     $                 'Rowwise', M, N, JB, WORK( IPW ), IV, JV, DESCW,
     $                 WORK( IPTP ), A, IA, JA, DESCA, WORK( IPWK ) )
      END IF
*
      CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, 0 )
*
      RETURN
*
*     End of PZGEBDRV
*
      END
