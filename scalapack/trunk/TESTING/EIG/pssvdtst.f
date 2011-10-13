      SUBROUTINE PSSVDTST( M, N, NPROW, NPCOL, NB, ISEED, THRESH, WORK,
     $                     RESULT, LWORK, NOUT )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997          
*
*     .. Scalar Arguments ..
      INTEGER            LWORK, M, N, NB, NOUT, NPCOL, NPROW
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 ), RESULT( 9 )
      REAL               WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PSSVDTST checks the singular value decomposition (SVD) routine
*  PSGESVD. PSGESVD factors A = U diag(S) VT, where U and VT are
*  orthogonal and diag(S) is diagonal with the entries of the array
*  S on its diagonal. The entries of S are the singular values, stored
*  in decreasing order. U and VT can be optionally not computed,
*  computed and overwritten on A, or computed partially.
*
*  A is M by N. Let SIZE = min( M, N ). S has dimension SIZE by SIZE.
*  U is M by SIZE and VT is SIZE by  N. PDGESVD optionally calculates
*  U and VT, depending on the values of its parameters JOBU and JOBVT.
*  There are four possible  combinations of "job" parameters for a call
*  to PDGESVD, that correspond to four values of internal index JOBTYPE.
*  The table below shows the mapping between "job" parameters of 
*  PDGESVD and respective values of the index JOBTYPE together 
*  with matrices computed for each type of the job.
*
*
*              |   JOBU = 'V'     |       JOBU = 'N'
*   ---------- -------------------------------------------
*   JOBVT = 'V'|  JOBTYPE = 1     |   JOBTYPE = 3
*              |  U1, S1, VT1     |   S3, VT3
*   ---------- ------------------------------------------
*   JOBVT = 'N'|  JOBTYPE = 2     |   JOBTYPE = 4
*              |  U2, S2          |   S4
*
*
*  When PSSVDTST is called, a number of matrix "types" are specified.
*  For each type of matrix, and for the minimal workspace as well as
*  for larger than minimal workspace an M x N matrix "A" with known
*  singular values is generated and used to test the SVD routines.
*  For each matrix, A will be factored as A = U diag(S) VT and the
*  following 9 tests computed:
*
*  (1)   | A - U1 diag(S1) VT1 | / ( |A| max(M,N) ulp )
*
*  (2)   | I - U1'U1 | / ( M ulp )
*
*  (3)   | I - VT1 VT1' | / ( N ulp ),
*
*  (4)   S1 contains SIZE  nonnegative values in decreasing order.
*        (Return 0 if true, 1/ULP if false.)
*
*  (5)   | S1 - S2 | / ( SIZE ulp |S| )
*
*  (6)   | U1 - U2 | / ( M ulp )
*
*  (7)   | S1 - S3 | / ( SIZE ulp |S| )
*
*  (8)   | VT1 - VT3 | / ( N ulp )
*
*  (9)   | S1 - S4 | / (  SIZE ulp |S| )
*
*  Currently, the list of possible  matrix types is:
*
*  (1)  The zero matrix.
*
*  (2)  The identity matrix.
*
*  (3)  A diagonal matrix with evenly spaced entries
*       1, ..., ULP.
*       (ULP = (first number larger than 1) - 1 )
*
*  (4)  A matrix of the form  U D VT, where U, VT are orthogonal and
*       D has evenly spaced entries 1, ..., ULP.
*
*  (5) Same as (4), but multiplied by SQRT( overflow threshold )
*
*  (6) Same as (4), but multiplied by SQRT( underflow threshold )
*
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
*  ==========
*
*  M      (global input) INTEGER dimension
*         The value of the matrix row dimension.
*
*  N      (global input) INTEGER dimension
*         The value of the matrix column dimension.
*
*  NPROW  (global input) INTEGER
*         Number of process rows
*
*  NPCOL  (global input) INTEGER
*         Number of process columns
*
*  NB     (global input) INTEGER
*         The block size of the matrix A. NB >=1.
*
*  ISEED  (global input/local output) INTEGER array, dimension (4)
*         On entry, the seed of the random number generator.  The array
*         elements should be between 0 and 4095; if not they will be
*         reduced mod 4096.  Also, ISEED(4) must be odd.
*         On exit, ISEED is changed and can be used in the next call to
*         SDRVBD to continue the same random number sequence.
*
*  THRESH (global input) REAL
*         The threshold value for the test ratios.  A result is
*         included in the output file if RESULT >= THRESH.  The test
*         ratios are scaled to be O(1), so THRESH should be a small
*         multiple of 1, e.g., 10 or 100.  To have every test ratio
*         printed, use THRESH = 0.
*
*  RESULT (global input/output) INTEGER array of dimension 9. Initially
*         RESULT( I ) = 0. On the output, RESULT ( I ) = 1 if test I
*         ( see above ) wasn't passed.
*
*  WORK   (local workspace) REAL array, dimension (LWORK)
*
*  LWORK  (local input) INTEGER
*         Dimension of the array WORK. It is defined as follows
*         LWORK = 1 + 2*LDA*NQ + 3*SIZE +
*         MAX(WPSLAGGE, LDU*SIZEQ + LDVT*NQ + MAX(LDU*SIZEQ, LDVT*NQ)
*         + WPSGESVD + MAX( WPSSVDCHK, WPSSVDCMP)),
*         where WPSLAGGE, WPSGESVD, WPSSVDCHK,  WPSSVDCMP are amounts
*         of workspace required respectively by PSLAGGE, PSGESVD,
*         PSSVDCHK, PSSVDCMP.
*         Here
*         LDA = NUMROC( M, NB, MYROW, 0, NPROW ), LDU = LDA,
*         LDVT = NUMROC( SIZE, NB, MYROW, 0, NPROW ),
*         NQ = NUMROC( N, NB, MYCOL, 0, NPCOL ),
*         SIZEQ = NUMROC( SIZE, NB, MYCOL, 0, NPCOL ).
*         Values of the variables WPSLAGGE, WPSGESVD, WPSSVDCHK,
*         WPSSVDCMP are found by "dummy" calls to
*         the respective routines. In every "dummy" call, variable
*         LWORK is set to -1, thus causing respective routine
*         immediately return required workspace in WORK(1) without
*         executing any calculations
*
*  NOUT   (local input) INTEGER
*         The unit number for output file.  Only used on node 0.
*         NOUT = 6, output to screen,
*         NOUT = 0, output to stderr.
*  =====================================================================
*
*     .. Parameters ..
      INTEGER             BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                    MB_, NB_, RSRC_, CSRC_, LLD_, NTYPES
      PARAMETER           ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                    CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                    RSRC_ = 7, CSRC_ = 8, LLD_ = 9, NTYPES = 6 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          HETERO, JOBU, JOBVT
      INTEGER            CONTEXT, DINFO, I, IA, IAM, INFO, ITYPE, IU,
     $                   IVT, JA, JOBTYPE, JU, JVT, LDA, LDU, LDVT,
     $                   LLWORK, LWMIN, MYCOL, MYROW, NNODES, NQ, PASS,
     $                   PTRA, PTRAC, PTRD, PTRS, PTRSC, PTRU, PTRUC, 
     $                   PTRVT, PTRVTC, PTRWORK, SETHET, SIZE, SIZEQ,
     $                   WPSGESVD, WPSLAGGE, WPSSVDCHK, WPSSVDCMP
      REAL               CHK, DELTA, H, MTM, OVFL, RTOVFL, RTUNFL, ULP,
     $                   UNFL
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $                   BLACS_SET,
     $                   DESCINIT, SGAMN2D, SGAMX2D, SLABAD, SSCAL,
     $                   IGAMN2D, IGAMX2D, IGEBR2D, IGEBS2D, PSELSET,
     $                   PSGESVD, PSLACPY, PSLAGGE, PSLASET, PSSVDCHK,
     $                   PSSVDCMP, PXERBLA, SLBOOT, SLCOMBINE, SLTIMER
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      REAL               PSLAMCH
      EXTERNAL           NUMROC, PSLAMCH
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCU( DLEN_ ),
     $                   DESCVT( DLEN_ ), ITMP( 2 )
      DOUBLE PRECISION   CTIME( 1 ), WTIME( 1 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*     This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*DTYPE_*LLD_*MB_*M_*NB_*N_*RSRC_.LT.0 )
     $   RETURN
*
      CALL BLACS_PINFO( IAM, NNODES )
      CALL BLACS_GET( -1, 0, CONTEXT )
      CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     If this process is not a part of the contex, bail out now.
*
      IF( ( MYROW.GE.NPROW ) .OR. ( MYROW.LT.0 ) .OR.
     $    ( MYCOL.GE.NPCOL ) .OR. ( MYCOL.LT.0 ) )GO TO 110
      CALL BLACS_SET( CONTEXT, 15, 1 )
      INFO = 0
*
*     Check input parameters.
*
      IF( M.LE.0 ) THEN
         INFO = -1
      ELSE IF( N.LE.0 ) THEN
         INFO = -2
      ELSE IF( NPROW.LE.0 ) THEN
         INFO = -3
      ELSE IF( NPCOL.LE.0 ) THEN
         INFO = -4
      ELSE IF( NB.LE.0 ) THEN
         INFO = -5
      ELSE IF( THRESH.LE.0 ) THEN
         INFO = -7
      END IF
*
      SIZE = MIN( M, N )
*
*     Initialize matrix descriptors.
*
      IA  = 1
      JA  = 1
      IU  = 1
      JU  = 1
      IVT = 1
      JVT = 1
*
      LDA = NUMROC( M, NB, MYROW, 0, NPROW )
      LDA = MAX( 1, LDA )
      NQ = NUMROC( N, NB, MYCOL, 0, NPCOL )
      LDU = LDA
      SIZEQ = NUMROC( SIZE, NB, MYCOL, 0, NPCOL )
      LDVT = NUMROC( SIZE, NB, MYROW, 0, NPROW )
      LDVT = MAX( 1, LDVT )
      CALL DESCINIT( DESCA, M, N, NB, NB, 0, 0, CONTEXT, LDA, DINFO )
      CALL DESCINIT( DESCU, M, SIZE, NB, NB, 0, 0, CONTEXT, LDU, DINFO )
      CALL DESCINIT( DESCVT, SIZE, N, NB, NB, 0, 0, CONTEXT, LDVT,
     $               DINFO )
*
*     Set some pointers to work array in order to do "dummy" calls.
*
      PTRA = 2
      PTRAC = PTRA + LDA*NQ
      PTRD = PTRAC + LDA*NQ
      PTRS = PTRD + SIZE
      PTRSC = PTRS + SIZE
      PTRWORK = PTRSC + SIZE
*
      PTRU = PTRWORK
      PTRVT = PTRWORK
      PTRUC = PTRWORK
      PTRVTC = PTRWORK
*
*     "Dummy" calls -- return required workspace in work(1) without
*     any calculation. 
*
      CALL PSLAGGE( M, N, WORK( PTRD ), WORK( PTRA ), IA, JA, DESCA,
     $              ISEED, SIZE, WORK( PTRWORK ), -1, DINFO )
      WPSLAGGE = INT( WORK( PTRWORK ) )
*
      CALL PSGESVD( 'V', 'V', M, N, WORK( PTRA ), IA, JA, DESCA, 
     $              WORK( PTRS ), WORK( PTRU ), IU, JU, DESCU, 
     $              WORK( PTRVT ), IVT, JVT, DESCVT,
     $              WORK( PTRWORK ), -1, DINFO )
      WPSGESVD = INT( WORK( PTRWORK ) )
*
      CALL PSSVDCHK( M, N, WORK( PTRAC ), IA, JA, DESCA, WORK( PTRUC ),
     $               IU, JU, DESCU, WORK( PTRVT ), IVT, JVT, DESCVT, 
     $               WORK( PTRS ), THRESH, WORK( PTRWORK ), -1, 
     $               RESULT, CHK, MTM )
      WPSSVDCHK = INT( WORK( PTRWORK ) )
*
      CALL PSSVDCMP( M, N, 1, WORK( PTRS ), WORK( PTRSC ), WORK( PTRU ),
     $               WORK( PTRUC ), IU, JU, DESCU, WORK( PTRVT ),
     $               WORK( PTRVTC ), IVT, JVT, DESCVT, THRESH, 
     $               RESULT, DELTA, WORK( PTRWORK ), -1 )
      WPSSVDCMP = INT( WORK( PTRWORK ) )
*
*     Calculation of workspace at last.
*
      LWMIN = 1 + 2*LDA*NQ + 3*SIZE +
     $        MAX( WPSLAGGE, LDU*SIZEQ+LDVT*NQ+MAX( LDU*SIZEQ,
     $        LDVT*NQ )+WPSGESVD+MAX( WPSSVDCHK, WPSSVDCMP ) )
      WORK( 1 ) = LWMIN
*
*     If this is a "dummy" call, return.
*
      IF( LWORK.EQ.-1 )
     $   GO TO 120
      IF( INFO.EQ.0 ) THEN
         IF( LWORK.LT.LWMIN ) THEN
            INFO = -10
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( DESCA( CTXT_ ), 'PSSVDTST', -INFO )
         RETURN
      END IF
*
      ULP = PSLAMCH( CONTEXT, 'P' )
      UNFL = PSLAMCH( CONTEXT, 'Safe min' )
      OVFL = ONE / UNFL
      CALL SLABAD( UNFL, OVFL )
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )
*
*     This ensures that everyone starts out with the same seed.
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         CALL IGEBS2D( CONTEXT, 'a', ' ', 4, 1, ISEED, 4 )
      ELSE
         CALL IGEBR2D( CONTEXT, 'a', ' ', 4, 1, ISEED, 4, 0, 0 )
      END IF
*
*     Loop over matrix types.
*
      DO 100 ITYPE = 1, NTYPES
*
         PASS = 0
         SETHET = 0
         PTRWORK = PTRSC + SIZE
         LLWORK = LWORK - PTRWORK + 1
*
*        Compute A.
*
         IF( ITYPE.EQ.1 ) THEN
*
*           Zero Matrix.
*
            DO 10 I = 1, SIZE
               WORK( PTRD+I-1 ) = ZERO
   10       CONTINUE
*
            CALL PSLASET( 'All', M, N, ZERO, ZERO, WORK( PTRA ), 
     $                    IA, JA, DESCA )
*
         ELSE IF( ITYPE.EQ.2 ) THEN
*
*           Identity Matrix.
*
            DO 20 I = 1, SIZE
               WORK( PTRD+I-1 ) = ONE
   20       CONTINUE
*
            CALL PSLASET( 'All', M, N, ZERO, ONE, WORK( PTRA ),
     $                    IA, JA, DESCA )
*
         ELSE IF( ITYPE.GT.2 ) THEN
*
*           Preset Singular Values.
*
            IF( SIZE.NE.1 ) THEN
               H = ( ULP-1 ) / ( SIZE-1 )
               DO 30 I = 1, SIZE
                  WORK( PTRD+I-1 ) = 1 + H*( I-1 )
   30          CONTINUE
            ELSE
               WORK( PTRD ) = 1
            END IF
*
            IF( ITYPE.EQ.3 ) THEN
*
*              Diagonal Matrix with specified singular values.
*
               CALL PSLASET( 'All', M, N, ZERO, ZERO, WORK( PTRA ),
     $                       IA, JA, DESCA )
*
               DO 40 I = 1, SIZE
                  CALL PSELSET( WORK( PTRA ), I, I, DESCA,
     $                          WORK( PTRD+I-1 ) )
   40          CONTINUE
*
            ELSE IF( ITYPE.EQ.4 ) THEN
*
*              General matrix with specified singular values.
*
               CALL PSLAGGE( M, N, WORK( PTRD ), WORK( PTRA ), IA, JA,
     $                       DESCA, ISEED, SIZE, WORK( PTRWORK ),
     $                       LLWORK, INFO )
*
            ELSE IF( ITYPE.EQ.5 ) THEN
*
*              Singular values scaled by overflow.
*
               CALL SSCAL( SIZE, RTOVFL, WORK( PTRD ), 1 )
*
               CALL PSLAGGE( M, N, WORK( PTRD ), WORK( PTRA ), IA, JA,
     $                       DESCA, ISEED, SIZE, WORK( PTRWORK ),
     $                       LLWORK, INFO )
*
            ELSE IF( ITYPE.EQ.6 ) THEN
*
*              Singular values scaled by underflow.
*
               CALL SSCAL( SIZE, RTUNFL, WORK( PTRD ), 1 )
               CALL PSLAGGE( M, N, WORK( PTRD ), WORK( PTRA ), IA, JA,
     $                       DESCA, ISEED, SIZE, WORK( PTRWORK ),
     $                       LLWORK, INFO )
*
            END IF
*
         END IF
*
*        Set mapping between JOBTYPE and calling parameters of
*        PSGESVD, reset pointers to WORK array to save space.
*
         DO 80 JOBTYPE = 1, 4
*
            IF( JOBTYPE.EQ.1 ) THEN
               JOBU = 'V'
               JOBVT = 'V'
               PTRVT = PTRU + LDU*SIZEQ
               PTRUC = PTRVT + LDVT*NQ
               PTRWORK = PTRUC + LDU*SIZEQ
               LLWORK = LWORK - PTRWORK + 1
            ELSE IF( JOBTYPE.EQ.2 ) THEN
               JOBU = 'V'
               JOBVT = 'N'
            ELSE IF( JOBTYPE.EQ.3 ) THEN
               JOBU = 'N'
               JOBVT = 'V'
               PTRVTC = PTRUC
               PTRWORK = PTRVTC + LDVT*NQ
               LLWORK = LWORK - PTRWORK + 1
            ELSE IF( JOBTYPE.EQ.4 ) THEN
               JOBU = 'N'
               JOBVT = 'N'
               PTRWORK = PTRUC
               LLWORK = LWORK - PTRWORK + 1
            END IF
*
*           Duplicate matrix A.
*
            CALL PSLACPY( 'A', M, N, WORK( PTRA ), IA, JA, DESCA,
     $                    WORK( PTRAC ), IA, JA, DESCA )
*
*           Test SVD calculation with minimum amount of workspace
*           calculated earlier.
*
            IF( JOBTYPE.EQ.1 ) THEN
*
*              Run SVD.
               CALL SLBOOT
               CALL BLACS_BARRIER( CONTEXT, 'All' )
               CALL SLTIMER( 1 )
*
               CALL PSGESVD( JOBU, JOBVT, M, N, WORK( PTRAC ), IA, JA,
     $                       DESCA, WORK( PTRS ), WORK( PTRU ), IU, JU,
     $                       DESCU, WORK( PTRVT ), IVT, JVT, DESCVT,
     $                       WORK( PTRWORK ), WPSGESVD, INFO )
*
               CALL SLTIMER( 1 )
               CALL SLCOMBINE( CONTEXT, 'All', '>', 'W', 1, 1, WTIME )
               CALL SLCOMBINE( CONTEXT, 'All', '>', 'C', 1, 1, CTIME )
*
*              Check INFO. Different INFO for  different processes mean
*              something went wrong.
*
               ITMP( 1 ) = INFO
               ITMP( 2 ) = INFO
*
               CALL IGAMN2D( DESCA( CTXT_ ), 'a', ' ', 1, 1, ITMP, 1, 1,
     $                       1, -1, -1, 0 )
               CALL IGAMX2D( DESCA( CTXT_ ), 'a', ' ', 1, 1, ITMP( 2 ),
     $                       1, 1, 1, -1, -1, 0 )
*
               IF( ITMP( 1 ).NE.ITMP( 2 ) ) THEN
                  IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
                     WRITE( NOUT, FMT = * )
     $                  'Different processes return different INFO'
                     GO TO 120
                  END IF
               END IF
*
*              If INFO is  negative PXERBLA tells you. So the only thing
*              is to check for positive INFO -- detected heterogeneous
*              system.
*
               IF( INFO.EQ.( SIZE+1 ) ) THEN
                  HETERO = 'P'
                  SETHET = 1
               END IF
*
*              If INFO was fine do more exhaustive check.
*
               IF( INFO.EQ.ZERO ) THEN
*
                  DO 50 I = 1, SIZE
                     WORK( I+PTRWORK ) = WORK( I+PTRS-1 )
                     WORK( I+SIZE+PTRWORK ) = WORK( I+PTRS-1 )
   50             CONTINUE
*
                  CALL SGAMN2D( DESCA( CTXT_ ), 'a', ' ', SIZE, 1,
     $                          WORK( 1+PTRWORK ), SIZE, 1, 1, -1, -1,
     $                          0 )
                  CALL SGAMX2D( DESCA( CTXT_ ), 'a', ' ', SIZE, 1,
     $                          WORK( 1+SIZE+PTRWORK ), SIZE, 1, 1, -1,
     $                          -1, 0 )
*
                  DO 60 I = 1, SIZE
                     IF( ABS( WORK( I+PTRWORK )-WORK( SIZE+I+
     $                   PTRWORK ) ).GT.ZERO ) THEN
                        WRITE( NOUT, FMT = * )'I= ', I, ' MIN=',
     $                     WORK( I+PTRWORK ), ' MAX=',
     $                     WORK( SIZE+I+PTRWORK )
                        HETERO = 'T'
                        SETHET = 1
                        GO TO 70
                     END IF
*
   60             CONTINUE
   70             CONTINUE
*
               END IF
*
               IF( SETHET.NE.1 )
     $            HETERO = 'N'
*
*              After PSGESVD AC got screwed up -- need to copy again.
*
               CALL PSLACPY( 'A', M, N, WORK( PTRA ), IA, JA, DESCA,
     $                       WORK( PTRAC ), IA, JA, DESCA )
*
*              PSSVDCHK screws up U. So before the call to PSSVDCHK
*              U is copied to UC and a pointer to UC is passed to
*              PSSVDCHK.
*
               CALL PSLACPY( 'A', M, SIZE, WORK( PTRU ), IU, JU, DESCU,
     $                       WORK( PTRUC ), IU, JU, DESCU )
*
*              Run tests 1 - 4.
*
               CALL PSSVDCHK( M, N, WORK( PTRAC ), IA, JA, DESCA,
     $                        WORK( PTRUC ), IU, JU, DESCU, 
     $                        WORK( PTRVT ), IVT, JVT, DESCVT, 
     $                        WORK( PTRS ), THRESH, WORK( PTRWORK ), 
     $                        LLWORK, RESULT, CHK, MTM )
*
            ELSE
*
*              Once again test PSGESVD with min workspace.
*
               CALL PSGESVD( JOBU, JOBVT, M, N, WORK( PTRAC ), IA, JA,
     $                       DESCA, WORK( PTRSC ), WORK( PTRUC ), IU,
     $                       JU, DESCU, WORK( PTRVTC ), IVT, JVT, 
     $                       DESCVT, WORK( PTRWORK ), WPSGESVD, INFO )
*
               CALL PSSVDCMP( M, N, JOBTYPE, WORK( PTRS ),
     $                        WORK( PTRSC ), WORK( PTRU ),
     $                        WORK( PTRUC ), IU, JU, DESCU,
     $                        WORK( PTRVT ), WORK( PTRVTC ), IVT, JVT,
     $                        DESCVT, THRESH, RESULT, DELTA,
     $                        WORK( PTRWORK ), LLWORK )
*
            END IF
*
   80    CONTINUE
*
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            DO 90 I = 1, 9
               IF( RESULT( I ).EQ.1 ) THEN
                  PASS = 1
                  WRITE( NOUT, FMT = * )'Test I = ', I, 'has failed'
                  WRITE( NOUT, FMT = * )' '
               END IF
   90       CONTINUE
            IF( PASS.EQ.0 ) THEN
               WRITE( NOUT, FMT = 9999 )'Passed', WTIME( 1 ),
     $            CTIME( 1 ), M, N, NPROW, NPCOL, NB, ITYPE, CHK, MTM,
     $            DELTA, HETERO
            END IF
         END IF
  100 CONTINUE
      CALL BLACS_GRIDEXIT( CONTEXT )
  110 CONTINUE
*
 9999 FORMAT( A6, 2E10.3, 2I6, 2I4, I5, I6, 3F6.2, 4X, A1 )
  120 CONTINUE
*
*     End of PSSVDTST
*
      RETURN
      END
