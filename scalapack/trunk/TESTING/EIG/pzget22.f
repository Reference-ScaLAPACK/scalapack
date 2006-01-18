      SUBROUTINE PZGET22( TRANSA, TRANSE, TRANSW, N, A, DESCA, E, DESCE,
     $                    W, WORK, DESCW, RWORK, RESULT )
*
*  -- ScaLAPACK testing routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     August 14, 2001
*
*     .. Scalar Arguments ..
      CHARACTER          TRANSA, TRANSE, TRANSW
      INTEGER            N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCE( * ), DESCW( * )
      DOUBLE PRECISION   RESULT( 2 ), RWORK( * )
      COMPLEX*16         A( * ), E( * ), W( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PZGET22 does an eigenvector check.
*
*  The basic test is:
*
*     RESULT(1) = | A E  -  E W | / ( |A| |E| ulp )
*
*  using the 1-norm.  It also tests the normalization of E:
*
*     RESULT(2) = max | m-norm(E(j)) - 1 | / ( n ulp )
*                  j
*
*  where E(j) is the j-th eigenvector, and m-norm is the max-norm of a
*  vector.  The max-norm of a complex n-vector x in this case is the
*  maximum of |re(x(i)| + |im(x(i)| over i = 1, ..., n.
*
*  Arguments
*  ==========
*
*  TRANSA  (input) CHARACTER*1
*          Specifies whether or not A is transposed.
*          = 'N':  No transpose
*          = 'T':  Transpose
*          = 'C':  Conjugate transpose
*
*  TRANSE  (input) CHARACTER*1
*          Specifies whether or not E is transposed.
*          = 'N':  No transpose, eigenvectors are in columns of E
*          = 'T':  Transpose, eigenvectors are in rows of E
*          = 'C':  Conjugate transpose, eigenvectors are in rows of E
*
*  TRANSW  (input) CHARACTER*1
*          Specifies whether or not W is transposed.
*          = 'N':  No transpose
*          = 'T':  Transpose, same as TRANSW = 'N'
*          = 'C':  Conjugate transpose, use -WI(j) instead of WI(j)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (*)
*          The matrix whose eigenvectors are in E.
*
*  DESCA   (input) INTEGER array, dimension(*)
*
*  E       (input) COMPLEX*16 array, dimension (*)
*          The matrix of eigenvectors. If TRANSE = 'N', the eigenvectors
*          are stored in the columns of E, if TRANSE = 'T' or 'C', the
*          eigenvectors are stored in the rows of E.
*
*  DESCE   (input) INTEGER array, dimension(*)
*
*  W       (input) COMPLEX*16 array, dimension (N)
*          The eigenvalues of A.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (*)
*  DESCW   (input) INTEGER array, dimension(*)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  RESULT  (output) DOUBLE PRECISION array, dimension (2)
*          RESULT(1) = | A E  -  E W | / ( |A| |E| ulp )
*          RESULT(2) = max | m-norm(E(j)) - 1 | / ( n ulp )
*                       j
*  Further Details
*  ===============
*
*  Contributed by Mark Fahey, June, 2000
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          NORMA, NORME
      INTEGER            ICOL, II, IROW, ITRNSE, ITRNSW, J, JCOL, JJ,
     $                   JROW, JVEC, LDA, LDE, LDW, MB, MYCOL, MYROW,
     $                   NB, NPCOL, NPROW, CONTXT, CA, CSRC, RA, RSRC
      DOUBLE PRECISION   ANORM, ENORM, ENRMAX, ENRMIN, ERRNRM, TEMP1,
     $                   ULP, UNFL
      COMPLEX*16         CDUM, WTEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PDLAMCH, PZLANGE
      EXTERNAL           LSAME, PDLAMCH, PZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGAMN2D, DGAMX2D, INFOG2L,
     $                   PZAXPY, PZGEMM, PZLASET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Initialize RESULT (in case N=0)
*
      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 )
     $   RETURN
*
      CONTXT = DESCA( CTXT_ )
      RSRC = DESCA( RSRC_ )
      CSRC = DESCA( CSRC_ )
      NB = DESCA( NB_ )
      MB = DESCA( MB_ )
      LDA = DESCA( LLD_ )
      LDE = DESCE( LLD_ )
      LDW = DESCW( LLD_ )
      CALL BLACS_GRIDINFO( CONTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      UNFL = PDLAMCH( CONTXT, 'Safe minimum' )
      ULP = PDLAMCH( CONTXT, 'Precision' )
*
      ITRNSE = 0
      ITRNSW = 0
      NORMA = 'O'
      NORME = 'O'
*
      IF( LSAME( TRANSA, 'T' ) .OR. LSAME( TRANSA, 'C' ) ) THEN
         NORMA = 'I'
      END IF
*
      IF( LSAME( TRANSE, 'T' ) ) THEN
         ITRNSE = 1
         NORME = 'I'
      ELSE IF( LSAME( TRANSE, 'C' ) ) THEN
         ITRNSE = 2
         NORME = 'I'
      END IF
*
      IF( LSAME( TRANSW, 'C' ) ) THEN
         ITRNSW = 1
      END IF
*
*     Normalization of E:
*
      ENRMIN = ONE / ULP
      ENRMAX = ZERO
      IF( ITRNSE.EQ.0 ) THEN
         DO 20 JVEC = 1, N
            TEMP1 = ZERO
            DO 10 J = 1, N
               CALL INFOG2L( J, JVEC, DESCE, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW, ICOL, II, JJ )
               IF( ( MYROW.EQ.II ) .AND. ( MYCOL.EQ.JJ ) ) THEN
                  TEMP1 = MAX( TEMP1, CABS1( E( ( ICOL-1 )*LDE+
     $                    IROW ) ) )
               END IF
   10       CONTINUE
            IF( MYCOL.EQ.JJ ) THEN
               CALL DGAMX2D( CONTXT, 'Col', ' ', 1, 1, TEMP1, 1, RA, CA,
     $                       -1, -1, -1 )
               ENRMIN = MIN( ENRMIN, TEMP1 )
               ENRMAX = MAX( ENRMAX, TEMP1 )
            END IF
   20    CONTINUE
         CALL DGAMX2D( CONTXT, 'Row', ' ', 1, 1, ENRMAX, 1, RA, CA, -1,
     $                 -1, -1 )
         CALL DGAMN2D( CONTXT, 'Row', ' ', 1, 1, ENRMIN, 1, RA, CA, -1,
     $                 -1, -1 )
      ELSE
         DO 40 J = 1, N
            TEMP1 = ZERO
            DO 30 JVEC = 1, N
               CALL INFOG2L( J, JVEC, DESCE, NPROW, NPCOL, MYROW, MYCOL,
     $                       IROW, ICOL, II, JJ )
               IF( ( MYROW.EQ.II ) .AND. ( MYCOL.EQ.JJ ) ) THEN
                  TEMP1 = MAX( TEMP1, CABS1( E( ( ICOL-1 )*LDE+
     $                    IROW ) ) )
               END IF
   30       CONTINUE
            IF( MYROW.EQ.II ) THEN
               CALL DGAMX2D( CONTXT, 'Row', ' ', 1, 1, TEMP1, 1, RA, CA,
     $                       -1, -1, -1 )
               ENRMIN = MIN( ENRMIN, TEMP1 )
               ENRMAX = MAX( ENRMAX, TEMP1 )
            END IF
   40    CONTINUE
         CALL DGAMX2D( CONTXT, 'Row', ' ', 1, 1, ENRMAX, 1, RA, CA, -1,
     $                 -1, -1 )
         CALL DGAMN2D( CONTXT, 'Row', ' ', 1, 1, ENRMIN, 1, RA, CA, -1,
     $                 -1, -1 )
      END IF
*
*     Norm of A:
*
      ANORM = MAX( PZLANGE( NORMA, N, N, A, 1, 1, DESCA, RWORK ), UNFL )
*
*     Norm of E:
*
      ENORM = MAX( PZLANGE( NORME, N, N, E, 1, 1, DESCE, RWORK ), ULP )
*
*     Norm of error:
*
*     Error =  AE - EW
*
      CALL PZLASET( 'Full', N, N, CZERO, CZERO, WORK, 1, 1, DESCW )
*
      DO 60 JCOL = 1, N
         IF( ITRNSW.EQ.0 ) THEN
            WTEMP = W( JCOL )
         ELSE
            WTEMP = DCONJG( W( JCOL ) )
         END IF
*
         IF( ITRNSE.EQ.0 ) THEN
            CALL PZAXPY( N, WTEMP, E, 1, JCOL, DESCE, 1, WORK, 1, JCOL,
     $                   DESCW, 1 )
         ELSE IF( ITRNSE.EQ.1 ) THEN
            CALL PZAXPY( N, WTEMP, E, JCOL, 1, DESCE, N, WORK, 1, JCOL,
     $                   DESCW, 1 )
         ELSE
            CALL PZAXPY( N, DCONJG( WTEMP ), E, JCOL, 1, DESCE, N, WORK,
     $                   1, JCOL, DESCW, 1 )
            DO 50 JROW = 1, N
               CALL INFOG2L( JROW, JCOL, DESCW, NPROW, NPCOL, MYROW,
     $                       MYCOL, IROW, ICOL, II, JJ )
               IF( ( MYROW.EQ.II ) .AND. ( MYCOL.EQ.JJ ) ) THEN
                  WORK( ( JCOL-1 )*LDW+JROW )
     $               = DCONJG( WORK( ( JCOL-1 )*LDW+JROW ) )
               END IF
   50       CONTINUE
         END IF
   60 CONTINUE
*
      CALL PZGEMM( TRANSA, TRANSE, N, N, N, CONE, A, 1, 1, DESCA, E, 1,
     $             1, DESCE, -CONE, WORK, 1, 1, DESCW )
*
      ERRNRM = PZLANGE( 'One', N, N, WORK, 1, 1, DESCW, RWORK ) / ENORM
*
*     Compute RESULT(1) (avoiding under/overflow)
*
      IF( ANORM.GT.ERRNRM ) THEN
         RESULT( 1 ) = ( ERRNRM / ANORM ) / ULP
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT( 1 ) = ( MIN( ERRNRM, ANORM ) / ANORM ) / ULP
         ELSE
            RESULT( 1 ) = MIN( ERRNRM / ANORM, ONE ) / ULP
         END IF
      END IF
*
*     Compute RESULT(2) : the normalization error in E.
*
      RESULT( 2 ) = MAX( ABS( ENRMAX-ONE ), ABS( ENRMIN-ONE ) ) /
     $              ( DBLE( N )*ULP )
*
      RETURN
*
*     End of PZGET22
*
      END
