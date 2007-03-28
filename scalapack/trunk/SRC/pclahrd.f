      SUBROUTINE PCLAHRD( N, K, NB, A, IA, JA, DESCA, TAU, T, Y, IY, JY,
     $                    DESCY, WORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER             IA, IY, JA, JY, K, N, NB
*     ..
*     .. Array Arguments ..
      INTEGER             DESCA( * ), DESCY( * )
      COMPLEX             A( * ), T( * ), TAU( * ), WORK( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  PCLAHRD reduces the first NB columns of a complex general
*  N-by-(N-K+1) distributed matrix A(IA:IA+N-1,JA:JA+N-K) so that
*  elements below the k-th subdiagonal are zero. The reduction is
*  performed by an unitary similarity transformation Q' * A * Q. The
*  routine returns the matrices V and T which determine Q as a block
*  reflector I - V*T*V', and also the matrix Y = A * V * T.
*
*  This is an auxiliary routine called by PCGEHRD. In the following
*  comments sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1).
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The number of rows and columns to be operated on, i.e. the
*          order of the distributed submatrix sub( A ).
*          N >= 0.
*
*  K       (global input) INTEGER
*          The offset for the reduction. Elements below the k-th
*          subdiagonal in the first NB columns are reduced to zero.
*
*  NB      (global input) INTEGER
*          The number of columns to be reduced.
*
*  A       (local input/local output) COMPLEX pointer into
*          the local memory to an array of dimension (LLD_A,
*          LOCc(JA+N-K)). On entry, this array contains the the local
*          pieces of the N-by-(N-K+1) general distributed matrix
*          A(IA:IA+N-1,JA:JA+N-K). On exit, the elements on and above
*          the k-th subdiagonal in the first NB columns are overwritten
*          with the corresponding elements of the reduced distributed
*          matrix; the elements below the k-th subdiagonal, with the
*          array TAU, represent the matrix Q as a product of elementary
*          reflectors. The other columns of A(IA:IA+N-1,JA:JA+N-K) are
*          unchanged. See Further Details.
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
*  TAU     (local output) COMPLEX array, dimension LOCc(JA+N-2)
*          The scalar factors of the elementary reflectors (see Further
*          Details). TAU is tied to the distributed matrix A.
*
*  T       (local output) COMPLEX array, dimension (NB_A,NB_A)
*          The upper triangular matrix T.
*
*  Y       (local output) COMPLEX pointer into the local memory
*          to an array of dimension (LLD_Y,NB_A). On exit, this array
*          contains the local pieces of the N-by-NB distributed
*          matrix Y. LLD_Y >= LOCr(IA+N-1).
*
*  IY      (global input) INTEGER
*          The row index in the global array Y indicating the first
*          row of sub( Y ).
*
*  JY      (global input) INTEGER
*          The column index in the global array Y indicating the
*          first column of sub( Y ).
*
*  DESCY   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Y.
*
*  WORK    (local workspace) COMPLEX array, dimension (NB)
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of nb elementary reflectors
*
*     Q = H(1) H(2) . . . H(nb).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
*  A(ia+i+k:ia+n-1,ja+i-1), and tau in TAU(ja+i-1).
*
*  The elements of the vectors v together form the (n-k+1)-by-nb matrix
*  V which is needed, with T and Y, to apply the transformation to the
*  unreduced part of the matrix, using an update of the form:
*  A(ia:ia+n-1,ja:ja+n-k) := (I-V*T*V')*(A(ia:ia+n-1,ja:ja+n-k)-Y*V').
*
*  The contents of A(ia:ia+n-1,ja:ja+n-k) on exit are illustrated by the
*  following example with n = 7, k = 3 and nb = 2:
*
*     ( a   h   a   a   a )
*     ( a   h   a   a   a )
*     ( a   h   a   a   a )
*     ( h   h   a   a   a )
*     ( v1  h   a   a   a )
*     ( v1  v2  a   a   a )
*     ( v1  v2  a   a   a )
*
*  where a denotes an element of the original matrix
*  A(ia:ia+n-1,ja:ja+n-k), h denotes a modified element of the upper
*  Hessenberg matrix H, and vi denotes an element of the vector
*  defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            IPROC
      INTEGER            I, IACOL, IAROW, ICTXT, IOFF, II, J, JJ, JL,
     $                   JT, JW, L, MYROW, MYCOL, NPCOL, NPROW, NQ
      COMPLEX            EI
*     ..
*     .. Local Arrays ..
      INTEGER            DESCW( DLEN_ )
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CAXPY, CCOPY, CSCAL,
     $                   CTRMV, DESCSET, INFOG2L, PCELSET,
     $                   PCGEMV, PCLACGV, PCLARFG, PCSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IOFF = MOD( JA-1, DESCA( NB_ ) )
      CALL INFOG2L( IA+K, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, II,
     $              JJ, IAROW, IACOL )
*
      IPROC = ( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL )
      NQ = NUMROC( N+JA-1, DESCA( NB_ ), MYCOL, IACOL, NPCOL )
      IF( MYCOL.EQ.IACOL )
     $   NQ = NQ - IOFF
*
      EI = ZERO

      JW = IOFF + 1
      CALL DESCSET( DESCW, 1, DESCA( MB_ ), 1, DESCA( MB_ ), IAROW,
     $              IACOL, ICTXT, 1 )
*
      DO 10 L = 1, NB
         I = IA + K + L - 2
         J = JA + L - 1
*
         IF( L.GT.1 ) THEN
*
*           Update A(ia:ia+n-1,j)
*
*           Compute i-th column of A - Y * V'
*
            CALL PCLACGV( L-1, A, I, JA, DESCA, DESCA( M_ ) )
            CALL PCGEMV( 'No transpose', N, L-1, -ONE, Y, IY, JY, DESCY,
     $                   A, I, JA, DESCA, DESCA( M_ ), ONE, A, IA, J,
     $                   DESCA, 1 )
            CALL PCLACGV( L-1, A, I, JA, DESCA, DESCA( M_ ) )
*
*           Apply I - V * T' * V' to this column (call it b) from the
*           left, using the last column of T as workspace
*
*           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
*                    ( V2 )             ( b2 )
*
*           where V1 is unit lower triangular
*
*           w := V1' * b1
*
            IF( IPROC ) THEN
               CALL CCOPY( L-1, A( (JJ+L-2)*DESCA( LLD_ )+II ), 1,
     $                     WORK( JW ), 1 )
               CALL CTRMV( 'Lower', 'Conjugate transpose', 'Unit', L-1,
     $                     A( (JJ-1)*DESCA( LLD_ )+II ), DESCA( LLD_ ),
     $                     WORK( JW ), 1 )
            END IF
*
*           w := w + V2'*b2
*
            CALL PCGEMV( 'Conjugate transpose', N-K-L+1, L-1, ONE, A,
     $                   I+1, JA, DESCA, A, I+1, J, DESCA, 1, ONE, WORK,
     $                   1, JW, DESCW, DESCW( M_ ) )
*
*           w := T'*w
*
            IF( IPROC )
     $         CALL CTRMV( 'Upper', 'Conjugate transpose', 'Non-unit',
     $                     L-1, T, DESCA( NB_ ), WORK( JW ), 1 )
*
*           b2 := b2 - V2*w
*
            CALL PCGEMV( 'No transpose', N-K-L+1, L-1, -ONE, A, I+1, JA,
     $                   DESCA, WORK, 1, JW, DESCW, DESCW( M_ ), ONE,
     $                   A, I+1, J, DESCA, 1 )
*
*           b1 := b1 - V1*w
*
            IF( IPROC ) THEN
               CALL CTRMV( 'Lower', 'No transpose', 'Unit', L-1,
     $                     A( (JJ-1)*DESCA( LLD_ )+II ), DESCA( LLD_ ),
     $                     WORK( JW ), 1 )
               CALL CAXPY( L-1, -ONE, WORK( JW ), 1,
     $                     A( ( JJ+L-2 )*DESCA( LLD_ )+II ), 1 )
            END IF
            CALL PCELSET( A, I, J-1, DESCA, EI )
         END IF
*
*        Generate the elementary reflector H(i) to annihilate
*        A(ia+k+i:ia+n-1,j)
*
         CALL PCLARFG( N-K-L+1, EI, I+1, J, A, MIN( I+2, N+IA-1 ), J,
     $                 DESCA, 1, TAU )
         CALL PCELSET( A, I+1, J, DESCA, ONE )
*
*        Compute  Y(iy:y+n-1,jy+l-1)
*
         CALL PCGEMV( 'No transpose', N, N-K-L+1, ONE, A, IA, J+1,
     $                DESCA, A, I+1, J, DESCA, 1, ZERO, Y, IY, JY+L-1,
     $                DESCY, 1 )
         CALL PCGEMV( 'Conjugate transpose', N-K-L+1, L-1, ONE, A, I+1,
     $                JA, DESCA, A, I+1, J, DESCA, 1, ZERO, WORK, 1, JW,
     $                DESCW, DESCW( M_ ) )
         CALL PCGEMV( 'No transpose', N, L-1, -ONE, Y, IY, JY, DESCY,
     $                WORK, 1, JW, DESCW, DESCW( M_ ), ONE, Y, IY,
     $                JY+L-1, DESCY, 1 )
         JL = MIN( JJ+L-1, JA+NQ-1 )
         CALL PCSCAL( N, TAU( JL ), Y, IY, JY+L-1, DESCY, 1 )
*
*        Compute T(1:i,i)
*
         IF( IPROC ) THEN
            JT = ( L-1 ) * DESCA( NB_ )
            CALL CSCAL( L-1, -TAU( JL ), WORK( JW ), 1 )
            CALL CCOPY( L-1, WORK( JW ), 1, T( JT+1 ), 1 )
            CALL CTRMV( 'Upper', 'No transpose', 'Non-unit', L-1, T,
     $                  DESCA( NB_ ), T( JT+1 ), 1 )
            T( JT+L ) = TAU( JL )
         END IF
   10 CONTINUE
*
      CALL PCELSET( A, K+NB+IA-1, J, DESCA, EI )
*
      RETURN
*
*     End of PCLAHRD
*
      END
