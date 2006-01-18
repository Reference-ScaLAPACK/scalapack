*
*
      SUBROUTINE PDLAEVSWP( N, ZIN, LDZI, Z, IZ, JZ, DESCZ, NVS, KEY,
     $                      WORK, LWORK )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 15, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IZ, JZ, LDZI, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCZ( * ), KEY( * ), NVS( * )
      DOUBLE PRECISION   WORK( * ), Z( * ), ZIN( LDZI, * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAEVSWP moves the eigenvectors (potentially unsorted) from
*  where they are computed, to a ScaLAPACK standard block cyclic
*  array, sorted so that the corresponding eigenvalues are sorted.
*
*  Notes
*  =====
*
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
*
*  Arguments
*  =========
*
*     NP = the number of rows local to a given process.
*     NQ = the number of columns local to a given process.
*
*  N       (global input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  ZIN     (local input) DOUBLE PRECISION array,
*          dimension ( LDZI, NVS(iam) )
*          The eigenvectors on input.  Each eigenvector resides entirely
*          in one process.  Each process holds a contiguous set of
*          NVS(iam) eigenvectors.  The first eigenvector which the
*          process holds is:  sum for i=[0,iam-1) of NVS(i)
*
*  LDZI    (locl input) INTEGER
*          leading dimension of the ZIN array
*
*  Z       (local output) DOUBLE PRECISION array
*          global dimension (N, N), local dimension (DESCZ(DLEN_), NQ)
*          The eigenvectors on output.  The eigenvectors are distributed
*          in a block cyclic manner in both dimensions, with a
*          block size of NB.
*
*  IZ      (global input) INTEGER
*          Z's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JZ      (global input) INTEGER
*          Z's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCZ   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix Z.
*
*  NVS     (global input) INTEGER array, dimension( nprocs+1 )
*          nvs(i) = number of processes
*          number of eigenvectors held by processes [0,i-1)
*          nvs(1) = number of eigen vectors held by [0,1-1) == 0
*          nvs(nprocs+1) = number of eigen vectors held by [0,nprocs) ==
*            total number of eigenvectors
*
*  KEY     (global input) INTEGER array, dimension( N )
*          Indicates the actual index (after sorting) for each of the
*          eigenvectors.
*
*  WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (local input) INTEGER dimension of WORK
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            CYCLIC_I, CYCLIC_J, DIST, I, IAM, II, INCII, J,
     $                   MAXI, MAXII, MINI, MINII, MYCOL, MYROW, NB,
     $                   NBUFSIZE, NPCOL, NPROCS, NPROW, PCOL, RECVCOL,
     $                   RECVFROM, RECVROW, SENDCOL, SENDROW, SENDTO
*     ..
*     .. External Functions ..
      INTEGER            INDXG2L, INDXG2P
      EXTERNAL           INDXG2L, INDXG2P
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGERV2D, DGESD2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
      CALL BLACS_GRIDINFO( DESCZ( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      IAM = MYROW + MYCOL*NPROW
      IAM = MYROW*NPCOL + MYCOL
*
      NB = DESCZ( MB_ )
*
      NPROCS = NPROW*NPCOL
*
*     If PxSTEIN operates on a sub-matrix of a global matrix, the
*     key [] that contains the indicies of the eigenvectors is refe-
*     renced to the dimensions of the sub-matrix and not the global
*     distrubited matrix. Because of this, PxLAEVSWP will incorrectly
*     map the eigenvectors to the global eigenvector matrix, Z, unless
*     the key[] elements are shifted as below.
*
      DO 10 J = DESCZ( N_ ), 1, -1
         KEY( J ) = KEY( J-JZ+1 ) + ( JZ-1 )
   10 CONTINUE
*
      DO 110 DIST = 0, NPROCS - 1
*
         SENDTO = MOD( IAM+DIST, NPROCS )
         RECVFROM = MOD( NPROCS+IAM-DIST, NPROCS )
*
         SENDROW = MOD( SENDTO, NPROW )
         SENDCOL = SENDTO / NPROW
         RECVROW = MOD( RECVFROM, NPROW )
         RECVCOL = RECVFROM / NPROW
*
         SENDROW = SENDTO / NPCOL
         SENDCOL = MOD( SENDTO, NPCOL )
         RECVROW = RECVFROM / NPCOL
         RECVCOL = MOD( RECVFROM, NPCOL )
*
*        Figure out what I have that process "sendto" wants
*
         NBUFSIZE = 0
*
*        We are looping through the eigenvectors that I presently own.
*
         DO 40 J = NVS( 1+IAM ) + JZ, NVS( 1+IAM+1 ) + JZ - 1
            PCOL = INDXG2P( KEY( J ), DESCZ( NB_ ), -1, DESCZ( CSRC_ ),
     $             NPCOL )
            IF( SENDCOL.EQ.PCOL ) THEN
               MINII = MOD( SENDROW+DESCZ( RSRC_ ), NPROW )*
     $                 DESCZ( MB_ ) + 1
               MAXII = DESCZ( M_ )
               INCII = DESCZ( MB_ )*NPROW
               DO 30 II = MINII, MAXII, INCII
                  MINI = MAX( II, IZ )
                  MAXI = MIN( II+DESCZ( MB_ )-1, N+IZ-1 )
                  DO 20 I = MINI, MAXI, 1
                     NBUFSIZE = NBUFSIZE + 1
                     WORK( NBUFSIZE ) = ZIN( I+1-IZ,
     $                                  J-NVS( 1+IAM )+1-JZ )
   20             CONTINUE
   30          CONTINUE
            END IF
   40    CONTINUE
*
*
         IF( MYROW.NE.SENDROW .OR. MYCOL.NE.SENDCOL )
     $      CALL DGESD2D( DESCZ( CTXT_ ), NBUFSIZE, 1, WORK, NBUFSIZE,
     $                    SENDROW, SENDCOL )
*
*
*        Figure out what process "recvfrom" has that I want
*
         NBUFSIZE = 0
         DO 70 J = NVS( 1+RECVFROM ) + JZ,
     $           NVS( 1+RECVFROM+1 ) + JZ - 1, 1
            PCOL = INDXG2P( KEY( J ), DESCZ( NB_ ), -1, DESCZ( CSRC_ ),
     $             NPCOL )
            IF( MYCOL.EQ.PCOL ) THEN
               MINII = MOD( MYROW+DESCZ( RSRC_ ), NPROW )*DESCZ( MB_ ) +
     $                 1
               MAXII = DESCZ( M_ )
               INCII = DESCZ( MB_ )*NPROW
               DO 60 II = MINII, MAXII, INCII
                  MINI = MAX( II, IZ )
                  MAXI = MIN( II+NB-1, N+IZ-1 )
                  DO 50 I = MINI, MAXI, 1
                     NBUFSIZE = NBUFSIZE + 1
   50             CONTINUE
   60          CONTINUE
            END IF
   70    CONTINUE
*
*
*
         IF( MYROW.NE.RECVROW .OR. MYCOL.NE.RECVCOL )
     $      CALL DGERV2D( DESCZ( CTXT_ ), 1, NBUFSIZE, WORK, 1, RECVROW,
     $                    RECVCOL )
*
         NBUFSIZE = 0
         DO 100 J = NVS( 1+RECVFROM ) + JZ,
     $           NVS( 1+RECVFROM+1 ) + JZ - 1, 1
            PCOL = INDXG2P( KEY( J ), DESCZ( NB_ ), -1, DESCZ( CSRC_ ),
     $             NPCOL )
            IF( MYCOL.EQ.PCOL ) THEN
               CYCLIC_J = INDXG2L( KEY( J ), DESCZ( MB_ ), -1, -1,
     $                    NPCOL )
               CYCLIC_I = 1
               MINII = MOD( MYROW+DESCZ( RSRC_ ), NPROW )*DESCZ( MB_ ) +
     $                 1
               MAXII = DESCZ( M_ )
               INCII = DESCZ( MB_ )*NPROW
               DO 90 II = MINII, MAXII, INCII
                  MINI = MAX( II, IZ )
                  CYCLIC_I = INDXG2L( MINI, DESCZ( MB_ ), -1, -1,
     $                       NPROW )
                  MAXI = MIN( II+NB-1, N+IZ-1 )
                  DO 80 I = MINI, MAXI, 1
                     NBUFSIZE = NBUFSIZE + 1
                     Z( CYCLIC_I+( CYCLIC_J-1 )*DESCZ( LLD_ ) )
     $                  = WORK( NBUFSIZE )
                     CYCLIC_I = CYCLIC_I + 1
   80             CONTINUE
   90          CONTINUE
            END IF
  100    CONTINUE
*
  110 CONTINUE
      RETURN
*
*     End of PDLAEVSWP
*
      END
