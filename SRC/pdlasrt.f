      SUBROUTINE PDLASRT( ID, N, D, Q, IQ, JQ, DESCQ, WORK, LWORK, 
     $                    IWORK, LIWORK, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver
*     May 1 2012
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, IQ, JQ, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCQ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), Q( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLASRT Sort the numbers in D in increasing order and the
*  corresponding vectors in Q.
*
*  Arguments
*  =========
*
*  ID      (global input) CHARACTER*1
*          = 'I': sort D in increasing order;
*          = 'D': sort D in decreasing order. (NOT IMPLEMENTED YET)
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( Q ). N >= 0.
*
*  D       (global input/output) DOUBLE PRECISION array, dimmension (N)
*          On exit, the number in D are sorted in increasing order.
*
*  Q       (local input) DOUBLE PRECISION pointer into the local memory
*          to an array of dimension (LLD_Q, LOCc(JQ+N-1) ). This array
*          contains the local pieces of the distributed matrix sub( A )
*          to be copied from.
*
*  IQ      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( Q ).
*
*  JQ      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( Q ).
*
*  DESCQ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  WORK    (local workspace/local output) DOUBLE PRECISION array,
*          dimension (LWORK)
*  LWORK   (local or global input) INTEGER
*          The dimension of the array WORK.
*          LWORK = MAX( N, NP * ( NB + NQ ))
*          where
*          NP = NUMROC( N, NB, MYROW, IAROW, NPROW ),
*          NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
*
*  IWORK   (local workspace/local output) INTEGER array,
*                                                  dimension (LIWORK)
*
*  LIWORK (local or global input) INTEGER
*          The dimension of the array IWORK.
*          LIWORK = N + 2*NB + 2*NPCOL
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            CL, COL, DUMMY, I, ICTXT, IID, IIQ, INDCOL,
     $                   INDX, INDXC, INDXG, IPQ, IPQ2, IPW, IPWORK, J,
     $                   JJQ, K, L, LDQ, LEND, LIWMIN, LWMIN, MYCOL,
     $                   MYROW, NB, ND, NP, NPCOL, NPROW, NQ, PSQ, QCOL,
     $                   QTOT, SBUF
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2L, INDXG2P, NUMROC
      EXTERNAL           INDXG2L, INDXG2P, LSAME, NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PXERBLA, DCOPY,
     $                   DGERV2D, DGESD2D, DLAMOV, DLAPST
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      IF( N.EQ.0 )
     $   RETURN
*
      ICTXT = DESCQ( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 600+CTXT_ )
      ELSE
         CALL CHK1MAT( N, 1, N, 1, IQ, JQ, DESCQ, 6, INFO )
         IF( INFO.EQ.0 ) THEN
            NB = DESCQ( NB_ )
            LDQ = DESCQ( LLD_ )
            NP = NUMROC( N, NB, MYROW, DESCQ( RSRC_ ), NPROW )
            NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
            LWMIN = MAX( N, NP*( NB+NQ ) )
            LIWMIN = N + 2*( NB+NPCOL )
            IF( .NOT.LSAME( ID, 'I' ) ) THEN
               INFO = -1
            ELSE IF( N.LT.0 ) THEN
               INFO = -2
            ELSE IF( LWORK.LT.LWMIN ) THEN
               INFO = -9
            ELSE IF( LIWORK.LT.LIWMIN ) THEN
               INFO = -11
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDLASRT', -INFO )
         RETURN
      END IF
*
*     Set Pointers
*
      INDXC = 1
      INDX = INDXC + N
      INDXG = INDX
      INDCOL = INDXG + NB
      QTOT = INDCOL + NB
      PSQ = QTOT + NPCOL
*
      IID = 1
      IPQ2 = 1
      IPW = IPQ2 + NP*NQ
*
      DUMMY = 0
      IIQ = INDXG2L( IQ, NB, DUMMY, DUMMY, NPROW )
*
*     Sort the eigenvalues in D
*
      CALL DLAPST( 'I', N, D, IWORK( INDX ), INFO )
*
      DO 10 L = 0, N - 1
         WORK( IID+L ) = D( IWORK( INDX+L ) )
         IWORK( INDXC-1+IWORK( INDX+L ) ) = IID + L
   10 CONTINUE
      CALL DCOPY( N, WORK, 1, D, 1 )
*
      ND = 0
   20 CONTINUE
      IF( ND.LT.N ) THEN
         LEND = MIN( NB, N-ND )
         J = JQ + ND
         QCOL = INDXG2P( J, NB, DUMMY, DESCQ( CSRC_ ), NPCOL )
         K = 0
         DO 30 L = 0, LEND - 1
            I = JQ - 1 + IWORK( INDXC+ND+L )
            CL = INDXG2P( I, NB, DUMMY, DESCQ( CSRC_ ), NPCOL )
            IWORK( INDCOL+L ) = CL
            IF( MYCOL.EQ.CL ) THEN
               IWORK( INDXG+K ) = IWORK( INDXC+ND+L )
               K = K + 1
            END IF
   30    CONTINUE
*
         IF( MYCOL.EQ.QCOL ) THEN
            DO 40 CL = 0, NPCOL - 1
               IWORK( QTOT+CL ) = 0
   40       CONTINUE
            DO 50 L = 0, LEND - 1
               IWORK( QTOT+IWORK( INDCOL+L ) ) = IWORK( QTOT+
     $            IWORK( INDCOL+L ) ) + 1
   50       CONTINUE
            IWORK( PSQ ) = 1
            DO 60 CL = 1, NPCOL - 1
               IWORK( PSQ+CL ) = IWORK( PSQ+CL-1 ) + IWORK( QTOT+CL-1 )
   60       CONTINUE
            DO 70 L = 0, LEND - 1
               CL = IWORK( INDCOL+L )
               I = JQ + ND + L
               JJQ = INDXG2L( I, NB, DUMMY, DUMMY, NPCOL )
               IPQ = IIQ + ( JJQ-1 )*LDQ
               IPWORK = IPW + ( IWORK( PSQ+CL )-1 )*NP
               CALL DCOPY( NP, Q( IPQ ), 1, WORK( IPWORK ), 1 )
               IWORK( PSQ+CL ) = IWORK( PSQ+CL ) + 1
   70       CONTINUE
            IWORK( PSQ ) = 1
            DO 80 CL = 1, NPCOL - 1
               IWORK( PSQ+CL ) = IWORK( PSQ+CL-1 ) + IWORK( QTOT+CL-1 )
   80       CONTINUE
            DO 90 L = 0, K - 1
               I = IWORK( INDXG+L )
               JJQ = INDXG2L( I, NB, DUMMY, DUMMY, NPCOL )
               IPQ = IPQ2 + ( JJQ-1 )*NP
               IPWORK = IPW + ( IWORK( PSQ+MYCOL )-1 )*NP
               CALL DCOPY( NP, WORK( IPWORK ), 1, WORK( IPQ ), 1 )
               IWORK( PSQ+MYCOL ) = IWORK( PSQ+MYCOL ) + 1
   90       CONTINUE
            DO 100 CL = 1, NPCOL - 1
               COL = MOD( MYCOL+CL, NPCOL )
               SBUF = IWORK( QTOT+COL )
               IF( SBUF.NE.0 ) THEN
                  IPWORK = IPW + ( IWORK( PSQ+COL )-1 )*NP
                  CALL DGESD2D( DESCQ( CTXT_ ), NP, SBUF,
     $                          WORK( IPWORK ), NP, MYROW, COL )
               END IF
  100       CONTINUE
*
         ELSE
*
            IF( K.NE.0 ) THEN
               CALL DGERV2D( DESCQ( CTXT_ ), NP, K, WORK( IPW ), NP,
     $                       MYROW, QCOL )
               DO 110 L = 0, K - 1
                  I = JQ - 1 + IWORK( INDXG+L )
                  JJQ = INDXG2L( I, NB, DUMMY, DUMMY, NPCOL )
                  IPQ = 1 + ( JJQ-1 )*NP
                  IPWORK = IPW + L*NP
                  CALL DCOPY( NP, WORK( IPWORK ), 1, WORK( IPQ ), 1 )
  110          CONTINUE
            END IF
         END IF
         ND = ND + NB
         GO TO 20
      END IF
      CALL DLAMOV( 'Full', NP, NQ, WORK, NP, Q( IIQ ), LDQ )
*
*     End of PDLASRT
*
      END
