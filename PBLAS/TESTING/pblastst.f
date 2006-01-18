      SUBROUTINE PVDIMCHK( ICTXT, NOUT, N, MATRIX, IX, JX, DESCX, INCX,
     $                     INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      INTEGER            ICTXT, INCX, INFO, IX, JX, N, NOUT
*     ..
*     .. Array Arguments ..
      INTEGER            DESCX( * )
*     ..
*
*  Purpose
*  =======
*
*  PVDIMCHK checks the validity of the input test dimensions. In case of
*  an invalid parameter or discrepancy between the parameters, this rou-
*  tine  displays  error  messages and returns an non-zero error code in
*  INFO.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information  required to
*  establish the  mapping  between a  matrix entry and its corresponding
*  process and memory location.
*
*  In  the  following  comments,   the character _  should  be  read  as
*  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN       EXPLANATION
*  ---------------- --------------- ------------------------------------
*  DTYPE_A (global) DESCA( DTYPE_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT_  ) The BLACS context handle, indicating
*                                   the NPROW x NPCOL BLACS process grid
*                                   A  is distributed over.  The context
*                                   itself  is  global,  but  the handle
*                                   (the integer value) may vary.
*  M_A     (global) DESCA( M_     ) The  number of rows in the distribu-
*                                   ted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N_     ) The number of columns in the distri-
*                                   buted matrix A, N_A >= 0.
*  IMB_A   (global) DESCA( IMB_   ) The number of rows of the upper left
*                                   block of the matrix A, IMB_A > 0.
*  INB_A   (global) DESCA( INB_   ) The  number  of columns of the upper
*                                   left   block   of   the   matrix  A,
*                                   INB_A > 0.
*  MB_A    (global) DESCA( MB_    ) The blocking factor used to  distri-
*                                   bute the last  M_A-IMB_A rows of  A,
*                                   MB_A > 0.
*  NB_A    (global) DESCA( NB_    ) The blocking factor used to  distri-
*                                   bute the last  N_A-INB_A  columns of
*                                   A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC_  ) The process row over which the first
*                                   row of the matrix  A is distributed,
*                                   NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC_  ) The  process  column  over which the
*                                   first  column of  A  is distributed.
*                                   NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD_   ) The  leading  dimension of the local
*                                   array  storing  the  local blocks of
*                                   the distributed matrix A,
*                                   IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                   ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_NUMROC:
*  Lr( IA, K ) = PB_NUMROC( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_NUMROC( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  NOUT    (global input) INTEGER
*          On entry, NOUT specifies the unit number for the output file.
*          When NOUT is 6, output to screen,  when  NOUT is 0, output to
*          stderr. NOUT is only defined for process 0.
*
*  MATRIX  (global input) CHARACTER*1
*          On entry,  MATRIX  specifies the one character matrix identi-
*          fier.
*
*  IX      (global input) INTEGER
*          On entry, IX  specifies X's global row index, which points to
*          the beginning of the submatrix sub( X ).
*
*  JX      (global input) INTEGER
*          On entry, JX  specifies X's global column index, which points
*          to the beginning of the submatrix sub( X ).
*
*  DESCX   (global and local input) INTEGER array
*          On entry, DESCX  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix X.
*
*  INCX    (global input) INTEGER
*          On entry,  INCX   specifies  the  global  increment  for  the
*          elements of  X.  Only two values of  INCX   are  supported in
*          this version, namely 1 and M_X. INCX  must not be zero.
*
*  INFO    (global output) INTEGER
*          On exit,  when  INFO  is  zero,  no  error has been detected,
*          otherwise an error has been detected.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      INTEGER         MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL        BLACS_GRIDINFO, IGSUM2D
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( N.LT.0 ) THEN
         INFO = 1
      ELSE IF( N.EQ.0 ) THEN
         IF( DESCX( M_ ).LT.0 )
     $      INFO = 1
         IF( DESCX( N_ ).LT.0 )
     $      INFO = 1
      ELSE
         IF( INCX.EQ.DESCX( M_ ) .AND.
     $      DESCX( N_ ).LT.( JX+N-1 ) ) THEN
            INFO = 1
         ELSE IF( INCX.EQ.1 .AND. INCX.NE.DESCX( M_ ) .AND.
     $      DESCX( M_ ).LT.( IX+N-1 ) ) THEN
            INFO = 1
         ELSE
            IF( IX.GT.DESCX( M_ ) ) THEN
               INFO = 1
            ELSE IF( JX.GT.DESCX( N_ ) ) THEN
               INFO = 1
            END IF
         END IF
      END IF
*
*     Check all processes for an error
*
      CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, 0 )
*
      IF( INFO.NE.0 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9999 ) MATRIX
            WRITE( NOUT, FMT = 9998 ) N, MATRIX, IX, MATRIX, JX, MATRIX,
     $                                INCX
            WRITE( NOUT, FMT = 9997 ) MATRIX, DESCX( M_ ), MATRIX,
     $                                DESCX( N_ )
            WRITE( NOUT, FMT = * )
         END IF
      END IF
*
 9999 FORMAT( 'Incompatible arguments for matrix ', A1, ':' )
 9998 FORMAT( 'N = ', I6, ', I', A1, ' = ', I6, ', J', A1, ' = ',
     $        I6, ',INC', A1, ' = ', I6 )
 9997 FORMAT( 'DESC', A1, '( M_ ) = ', I6, ', DESC', A1, '( N_ ) = ',
     $        I6, '.' )
*
      RETURN
*
*     End of PVDIMCHK
*
      END
      SUBROUTINE PMDIMCHK( ICTXT, NOUT, M, N, MATRIX, IA, JA, DESCA,
     $                     INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      INTEGER            ICTXT, INFO, IA, JA, M, N, NOUT
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
*     ..
*
*  Purpose
*  =======
*
*  PMDIMCHK checks the validity of the input test dimensions. In case of
*  an invalid parameter or discrepancy between the parameters, this rou-
*  tine  displays  error  messages and returns an non-zero error code in
*  INFO.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information  required to
*  establish the  mapping  between a  matrix entry and its corresponding
*  process and memory location.
*
*  In  the  following  comments,   the character _  should  be  read  as
*  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN       EXPLANATION
*  ---------------- --------------- ------------------------------------
*  DTYPE_A (global) DESCA( DTYPE_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT_  ) The BLACS context handle, indicating
*                                   the NPROW x NPCOL BLACS process grid
*                                   A  is distributed over.  The context
*                                   itself  is  global,  but  the handle
*                                   (the integer value) may vary.
*  M_A     (global) DESCA( M_     ) The  number of rows in the distribu-
*                                   ted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N_     ) The number of columns in the distri-
*                                   buted matrix A, N_A >= 0.
*  IMB_A   (global) DESCA( IMB_   ) The number of rows of the upper left
*                                   block of the matrix A, IMB_A > 0.
*  INB_A   (global) DESCA( INB_   ) The  number  of columns of the upper
*                                   left   block   of   the   matrix  A,
*                                   INB_A > 0.
*  MB_A    (global) DESCA( MB_    ) The blocking factor used to  distri-
*                                   bute the last  M_A-IMB_A rows of  A,
*                                   MB_A > 0.
*  NB_A    (global) DESCA( NB_    ) The blocking factor used to  distri-
*                                   bute the last  N_A-INB_A  columns of
*                                   A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC_  ) The process row over which the first
*                                   row of the matrix  A is distributed,
*                                   NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC_  ) The  process  column  over which the
*                                   first  column of  A  is distributed.
*                                   NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD_   ) The  leading  dimension of the local
*                                   array  storing  the  local blocks of
*                                   the distributed matrix A,
*                                   IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                   ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_NUMROC:
*  Lr( IA, K ) = PB_NUMROC( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_NUMROC( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  NOUT    (global input) INTEGER
*          On entry, NOUT specifies the unit number for the output file.
*          When NOUT is 6, output to screen,  when  NOUT is 0, output to
*          stderr. NOUT is only defined for process 0.
*
*  MATRIX  (global input) CHARACTER*1
*          On entry,  MATRIX  specifies the one character matrix identi-
*          fier.
*
*  IA      (global input) INTEGER
*          On entry, IA  specifies A's global row index, which points to
*          the beginning of the submatrix sub( A ).
*
*  JA      (global input) INTEGER
*          On entry, JA  specifies A's global column index, which points
*          to the beginning of the submatrix sub( A ).
*
*  DESCA   (global and local input) INTEGER array
*          On entry, DESCA  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix A.
*
*  INFO    (global output) INTEGER
*          On exit,  when  INFO  is  zero,  no  error has been detected,
*          otherwise an error has been detected.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      INTEGER         MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL        BLACS_GRIDINFO, IGSUM2D
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( ( M.LT.0 ).OR.( N.LT.0 ) ) THEN
         INFO = 1
      ELSE IF( ( M.EQ.0 ).OR.( N.EQ.0 ) )THEN
         IF( DESCA( M_ ).LT.0 )
     $      INFO = 1
         IF( DESCA( N_ ).LT.0 )
     $      INFO = 1
      ELSE
         IF( DESCA( M_ ).LT.( IA+M-1 ) )
     $      INFO = 1
         IF( DESCA( N_ ).LT.( JA+N-1 ) )
     $      INFO = 1
      END IF
*
*     Check all processes for an error
*
      CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, 0 )
*
      IF( INFO.NE.0 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9999 ) MATRIX
            WRITE( NOUT, FMT = 9998 ) M, N, MATRIX, IA, MATRIX, JA
            WRITE( NOUT, FMT = 9997 ) MATRIX, DESCA( M_ ), MATRIX,
     $                                DESCA( N_ )
            WRITE( NOUT, FMT = * )
         END IF
      END IF
*
 9999 FORMAT( 'Incompatible arguments for matrix ', A1, ':' )
 9998 FORMAT( 'M = ', I6, ', N = ', I6, ', I', A1, ' = ', I6,
     $        ', J', A1, ' = ', I6 )
 9997 FORMAT( 'DESC', A1, '( M_ ) = ', I6, ', DESC', A1, '( N_ ) = ',
     $        I6, '.' )
*
      RETURN
*
*     End of PMDIMCHK
*
      END
      SUBROUTINE PVDESCCHK( ICTXT, NOUT, MATRIX, DESCX, DTX, MX, NX,
     $                      IMBX, INBX, MBX, NBX, RSRCX, CSRCX, INCX,
     $                      MPX, NQX, IPREX, IMIDX, IPOSTX, IGAP,
     $                      GAPMUL, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      INTEGER            CSRCX, DTX, GAPMUL, ICTXT, IGAP, IMBX, IMIDX,
     $                   INBX, INCX, INFO, IPOSTX, IPREX, MBX, MPX, MX,
     $                   NBX, NOUT, NQX, NX, RSRCX
*     ..
*     .. Array Arguments ..
      INTEGER            DESCX( * )
*     ..
*
*  Purpose
*  =======
*
*  PVDESCCHK  checks  the validity of the input test parameters and ini-
*  tializes  the  descriptor DESCX and the scalar variables MPX, NQX. In
*  case  of  an  invalid parameter, this routine displays error messages
*  and return an non-zero error code in INFO.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information  required to
*  establish the  mapping  between a  matrix entry and its corresponding
*  process and memory location.
*
*  In  the  following  comments,   the character _  should  be  read  as
*  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN       EXPLANATION
*  ---------------- --------------- ------------------------------------
*  DTYPE_A (global) DESCA( DTYPE_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT_  ) The BLACS context handle, indicating
*                                   the NPROW x NPCOL BLACS process grid
*                                   A  is distributed over.  The context
*                                   itself  is  global,  but  the handle
*                                   (the integer value) may vary.
*  M_A     (global) DESCA( M_     ) The  number of rows in the distribu-
*                                   ted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N_     ) The number of columns in the distri-
*                                   buted matrix A, N_A >= 0.
*  IMB_A   (global) DESCA( IMB_   ) The number of rows of the upper left
*                                   block of the matrix A, IMB_A > 0.
*  INB_A   (global) DESCA( INB_   ) The  number  of columns of the upper
*                                   left   block   of   the   matrix  A,
*                                   INB_A > 0.
*  MB_A    (global) DESCA( MB_    ) The blocking factor used to  distri-
*                                   bute the last  M_A-IMB_A rows of  A,
*                                   MB_A > 0.
*  NB_A    (global) DESCA( NB_    ) The blocking factor used to  distri-
*                                   bute the last  N_A-INB_A  columns of
*                                   A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC_  ) The process row over which the first
*                                   row of the matrix  A is distributed,
*                                   NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC_  ) The  process  column  over which the
*                                   first  column of  A  is distributed.
*                                   NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD_   ) The  leading  dimension of the local
*                                   array  storing  the  local blocks of
*                                   the distributed matrix A,
*                                   IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                   ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_NUMROC:
*  Lr( IA, K ) = PB_NUMROC( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_NUMROC( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  NOUT    (global input) INTEGER
*          On entry, NOUT specifies the unit number for the output file.
*          When NOUT is 6, output to screen,  when  NOUT is 0, output to
*          stderr. NOUT is only defined for process 0.
*
*  MATRIX  (global input) CHARACTER*1
*          On entry,  MATRIX  specifies the one character matrix identi-
*          fier.
*
*  DESCX   (global output) INTEGER array
*          On entry, DESCX  is an array of dimension DLEN_. DESCX is the
*          array descriptor to be set.
*
*  DTYPEX  (global input) INTEGER
*          On entry, DTYPEX  specifies the descriptor type. In this ver-
*          sion, DTYPEX must be BLOCK_CYCLIC_INB_2D.
*
*  MX      (global input) INTEGER
*          On entry, MX  specifies the number of rows in the matrix.  MX
*          must be at least zero.
*
*  NX      (global input) INTEGER
*          On  entry,  NX specifies the number of columns in the matrix.
*          NX must be at least zero.
*
*  IMBX    (global input) INTEGER
*          On entry, IMBX specifies the row blocking factor used to dis-
*          tribute  the  first  IMBX rows of the matrix. IMBX must be at
*          least one.
*
*  INBX    (global input) INTEGER
*          On entry,  INBX  specifies the column blocking factor used to
*          distribute  the  first  INBX columns of the matrix. INBX must
*          be at least one.
*
*  MBX     (global input) INTEGER
*          On entry, MBX  specifies the row blocking factor used to dis-
*          tribute the rows of the matrix. MBX must be at least one.
*
*  NBX     (global input) INTEGER
*          On entry, NBX  specifies  the  column blocking factor used to
*          distribute  the  columns  of the matrix. NBX must be at least
*          one.
*
*  RSRCX   (global input) INTEGER
*          On entry, RSRCX  specifies the process row in which the first
*          row  of  the  matrix resides. When RSRCX is -1, the matrix is
*          row replicated,  otherwise  RSCRX  must  be at least zero and
*          strictly less than NPROW.
*
*  CSRCX   (global input) INTEGER
*          On entry,  CSRCX  specifies  the  process column in which the
*          first column of the matrix resides.  When  CSRCX  is -1,  the
*          matrix is column replicated, otherwise CSCRX must be at least
*          zero and strictly less than NPCOL.
*
*  INCX    (global input) INTEGER
*          On entry,  INCX  specifies  the global vector increment. INCX
*          must be one or MX.
*
*  MPX     (local output) INTEGER
*          On exit, MPX is Lr( 1, MX ).
*
*  NQX     (local output) INTEGER
*          On exit, NQX is Lc( 1, NX ).
*
*  IPREX   (local output) INTEGER
*          On exit,  IPREX  specifies  the size of the guard zone to put
*          before the start of the local padded array.
*
*  IMIDX   (local output) INTEGER
*          On exit,  IMIDX  specifies  the  ldx-gap of the guard zone to
*          put after each column of the local padded array.
*
*  IPOSTX  (local output) INTEGER
*          On exit,  IPOSTX  specifies the size of the guard zone to put
*          after the local padded array.
*
*  IGAP    (global input) INTEGER
*          On entry, IGAP specifies the size of the ldx-gap.
*
*  GAPMUL  (global input) INTEGER
*          On entry,  GAPMUL  is  a constant factor controlling the size
*          of the pre- and post guardzone.
*
*  INFO    (global output) INTEGER
*          On exit,  when  INFO  is  zero,  no  error has been detected,
*          otherwise an error has been detected.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      INTEGER            LLDX, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGSUM2D, PB_DESCINIT2
*     ..
*     .. External Functions ..
      INTEGER            PB_NUMROC
      EXTERNAL           PB_NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Verify descriptor type DTYPE_
*
      IF( DTX.NE.BLOCK_CYCLIC_2D_INB ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9999 ) MATRIX, 'DTYPE', MATRIX, DTX,
     $                                BLOCK_CYCLIC_2D_INB
         INFO = 1
      END IF
*
*     Verify global matrix dimensions (M_,N_) are correct
*
      IF( MX.LT.0 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9998 ) MATRIX, 'M', MATRIX, MX
         INFO = 1
      ELSE IF( NX.LT.0 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9997 ) MATRIX, 'N', MATRIX, NX
         INFO = 1
      END IF
*
*     Verify if blocking factors (IMB_, INB_) are correct
*
      IF( IMBX.LT.1 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9996 ) MATRIX, 'IMB', MATRIX, IMBX
         INFO = 1
      ELSE IF( INBX.LT.1 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9995 ) MATRIX, 'INB', MATRIX, INBX
         INFO = 1
      END IF
*
*     Verify if blocking factors (MB_, NB_) are correct
*
      IF( MBX.LT.1 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9994 ) MATRIX, 'MB', MATRIX, MBX
         INFO = 1
      ELSE IF( NBX.LT.1 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9993 ) MATRIX, 'NB', MATRIX, NBX
         INFO = 1
      END IF
*
*     Verify if origin process coordinates (RSRC_, CSRC_) are valid
*
      IF( RSRCX.LT.-1 .OR. RSRCX.GE.NPROW ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9992 ) MATRIX
            WRITE( NOUT, FMT = 9990 ) 'RSRC', MATRIX, RSRCX, NPROW
         END IF
         INFO = 1
      ELSE IF( CSRCX.LT.-1 .OR. CSRCX.GE.NPCOL ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9991 ) MATRIX
            WRITE( NOUT, FMT = 9990 ) 'CSRC', MATRIX, CSRCX, NPCOL
         END IF
         INFO = 1
      END IF
*
*     Check input increment value
*
      IF( INCX.NE.1 .AND. INCX.NE.MX ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9989 ) MATRIX
            WRITE( NOUT, FMT = 9988 ) 'INC', MATRIX, INCX, MATRIX, MX
         END IF
         INFO = 1
      END IF
*
*     Check all processes for an error
*
      CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, 0 )
*
      IF( INFO.NE.0 ) THEN
*
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9987 ) MATRIX
            WRITE( NOUT, FMT = * )
         END IF
*
      ELSE
*
*        Compute local testing leading dimension
*
         MPX    = PB_NUMROC( MX, 1, IMBX, MBX, MYROW, RSRCX, NPROW )
         NQX    = PB_NUMROC( NX, 1, INBX, NBX, MYCOL, CSRCX, NPCOL )
         IPREX  = MAX( GAPMUL*NBX, MPX )
         IMIDX  = IGAP
         IPOSTX = MAX( GAPMUL*NBX, NQX )
         LLDX   = MAX( 1, MPX ) + IMIDX
*
         CALL PB_DESCINIT2( DESCX, MX, NX, IMBX, INBX, MBX, NBX, RSRCX,
     $                      CSRCX, ICTXT, LLDX, INFO )
*
*        Check all processes for an error
*
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, 0 )
*
         IF( INFO.NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
               WRITE( NOUT, FMT = 9987 ) MATRIX
               WRITE( NOUT, FMT = * )
            END IF
         END IF
*
      END IF
*
 9999 FORMAT( 2X, '>> Invalid matrix ', A1, ' descriptor type ', A5, A1,
     $        ': ', I6, ' should be ', I3, '.' )
 9998 FORMAT( 2X, '>> Invalid matrix ', A1, ' row dimension ', A1, A1,
     $        ': ', I6, ' should be at least 1.' )
 9997 FORMAT( 2X, '>> Invalid matrix ', A1, ' column dimension ', A1,
     $        A1, ': ', I6, ' should be at least 1.' )
 9996 FORMAT( 2X, '>> Invalid matrix ', A1, ' first row block size ',
     $        A3, A1, ': ', I6, ' should be at least 1.' )
 9995 FORMAT( 2X, '>> Invalid matrix ', A1, ' first column block size ',
     $        A3, A1,': ', I6, ' should be at least 1.' )
 9994 FORMAT( 2X, '>> Invalid matrix ', A1, ' row block size ', A2, A1,
     $        ': ', I6, ' should be at least 1.' )
 9993 FORMAT( 2X, '>> Invalid matrix ', A1, ' column block size ', A2,
     $        A1,': ', I6, ' should be at least 1.' )
 9992 FORMAT( 2X, '>> Invalid matrix ', A1, ' row process source:' )
 9991 FORMAT( 2X, '>> Invalid matrix ', A1, ' column process source:' )
 9990 FORMAT( 2X, '>> ', A4, A1, '= ', I6, ' should be >= -1 and < ',
     $        I6, '.' )
 9989 FORMAT( 2X, '>> Invalid vector ', A1, ' increment:' )
 9988 FORMAT( 2X, '>> ', A3, A1, '= ', I6, ' should be 1 or M', A1,
     $        ' = ', I6, '.' )
 9987 FORMAT( 2X, '>> Invalid matrix ', A1, ' descriptor: going on to ',
     $        'next test case.' )
*
      RETURN
*
*     End of PVDESCCHK
*
      END
      SUBROUTINE PMDESCCHK( ICTXT, NOUT, MATRIX, DESCA, DTA, MA, NA,
     $                      IMBA, INBA, MBA, NBA, RSRCA, CSRCA, MPA,
     $                      NQA, IPREA, IMIDA, IPOSTA, IGAP, GAPMUL,
     $                      INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      INTEGER            CSRCA, DTA, GAPMUL, ICTXT, IGAP, IMBA, IMIDA,
     $                   INBA, INFO, IPOSTA, IPREA, MA, MBA, MPA, NA,
     $                   NBA, NOUT, NQA, RSRCA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
*     ..
*
*  Purpose
*  =======
*
*  PMDESCCHK  checks  the validity of the input test parameters and ini-
*  tializes  the  descriptor DESCA and the scalar variables MPA, NQA. In
*  case  of  an  invalid parameter, this routine displays error messages
*  and return an non-zero error code in INFO.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information  required to
*  establish the  mapping  between a  matrix entry and its corresponding
*  process and memory location.
*
*  In  the  following  comments,   the character _  should  be  read  as
*  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN       EXPLANATION
*  ---------------- --------------- ------------------------------------
*  DTYPE_A (global) DESCA( DTYPE_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT_  ) The BLACS context handle, indicating
*                                   the NPROW x NPCOL BLACS process grid
*                                   A  is distributed over.  The context
*                                   itself  is  global,  but  the handle
*                                   (the integer value) may vary.
*  M_A     (global) DESCA( M_     ) The  number of rows in the distribu-
*                                   ted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N_     ) The number of columns in the distri-
*                                   buted matrix A, N_A >= 0.
*  IMB_A   (global) DESCA( IMB_   ) The number of rows of the upper left
*                                   block of the matrix A, IMB_A > 0.
*  INB_A   (global) DESCA( INB_   ) The  number  of columns of the upper
*                                   left   block   of   the   matrix  A,
*                                   INB_A > 0.
*  MB_A    (global) DESCA( MB_    ) The blocking factor used to  distri-
*                                   bute the last  M_A-IMB_A rows of  A,
*                                   MB_A > 0.
*  NB_A    (global) DESCA( NB_    ) The blocking factor used to  distri-
*                                   bute the last  N_A-INB_A  columns of
*                                   A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC_  ) The process row over which the first
*                                   row of the matrix  A is distributed,
*                                   NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC_  ) The  process  column  over which the
*                                   first  column of  A  is distributed.
*                                   NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD_   ) The  leading  dimension of the local
*                                   array  storing  the  local blocks of
*                                   the distributed matrix A,
*                                   IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                   ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_NUMROC:
*  Lr( IA, K ) = PB_NUMROC( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_NUMROC( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  NOUT    (global input) INTEGER
*          On entry, NOUT specifies the unit number for the output file.
*          When NOUT is 6, output to screen,  when  NOUT is 0, output to
*          stderr. NOUT is only defined for process 0.
*
*  MATRIX  (global input) CHARACTER*1
*          On entry,  MATRIX  specifies the one character matrix identi-
*          fier.
*
*  DESCA   (global output) INTEGER array
*          On entry, DESCA  is an array of dimension DLEN_. DESCA is the
*          array descriptor to be set.
*
*  DTYPEA  (global input) INTEGER
*          On entry, DTYPEA  specifies the descriptor type. In this ver-
*          sion, DTYPEA must be BLOCK_CYCLIC_INB_2D.
*
*  MA      (global input) INTEGER
*          On entry, MA  specifies the number of rows in the matrix.  MA
*          must be at least zero.
*
*  NA      (global input) INTEGER
*          On  entry,  NA specifies the number of columns in the matrix.
*          NA must be at least zero.
*
*  IMBA    (global input) INTEGER
*          On entry, IMBA specifies the row blocking factor used to dis-
*          tribute  the  first  IMBA rows of the matrix. IMBA must be at
*          least one.
*
*  INBA    (global input) INTEGER
*          On entry,  INBA  specifies the column blocking factor used to
*          distribute  the  first  INBA columns of the matrix. INBA must
*          be at least one.
*
*  MBA     (global input) INTEGER
*          On entry, MBA  specifies the row blocking factor used to dis-
*          tribute the rows of the matrix. MBA must be at least one.
*
*  NBA     (global input) INTEGER
*          On entry, NBA  specifies  the  column blocking factor used to
*          distribute  the  columns  of the matrix. NBA must be at least
*          one.
*
*  RSRCA   (global input) INTEGER
*          On entry, RSRCA  specifies the process row in which the first
*          row  of  the  matrix resides. When RSRCA is -1, the matrix is
*          row replicated,  otherwise  RSCRA  must  be at least zero and
*          strictly less than NPROW.
*
*  CSRCA   (global input) INTEGER
*          On entry,  CSRCA  specifies  the  process column in which the
*          first column of the matrix resides.  When  CSRCA  is -1,  the
*          matrix is column replicated, otherwise CSCRA must be at least
*          zero and strictly less than NPCOL.
*
*  MPA     (local output) INTEGER
*          On exit, MPA is Lr( 1, MA ).
*
*  NQA     (local output) INTEGER
*          On exit, NQA is Lc( 1, NA ).
*
*  IPREA   (local output) INTEGER
*          On exit,  IPREA  specifies  the size of the guard zone to put
*          before the start of the local padded array.
*
*  IMIDA   (local output) INTEGER
*          On exit,  IMIDA  specifies  the  lda-gap of the guard zone to
*          put after each column of the local padded array.
*
*  IPOSTA  (local output) INTEGER
*          On exit,  IPOSTA  specifies the size of the guard zone to put
*          after the local padded array.
*
*  IGAP    (global input) INTEGER
*          On entry, IGAP specifies the size of the lda-gap.
*
*  GAPMUL  (global input) INTEGER
*          On entry,  GAPMUL  is  a constant factor controlling the size
*          of the pre- and post guardzone.
*
*  INFO    (global output) INTEGER
*          On exit,  when  INFO  is  zero,  no  error has been detected,
*          otherwise an error has been detected.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      INTEGER            LLDA, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGSUM2D, PB_DESCINIT2
*     ..
*     .. External Functions ..
      INTEGER            PB_NUMROC
      EXTERNAL           PB_NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Verify descriptor type DTYPE_
*
      IF( DTA.NE.BLOCK_CYCLIC_2D_INB ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9999 ) MATRIX, 'DTYPE', MATRIX, DTA,
     $                                BLOCK_CYCLIC_2D_INB
         INFO = 1
      END IF
*
*     Verify global matrix dimensions (M_,N_) are correct
*
      IF( MA.LT.0 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9998 ) MATRIX, 'M', MATRIX, MA
         INFO = 1
      ELSE IF( NA.LT.0 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9997 ) MATRIX, 'N', MATRIX, NA
         INFO = 1
      END IF
*
*     Verify if blocking factors (IMB_, INB_) are correct
*
      IF( IMBA.LT.1 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9996 ) MATRIX, 'IMB', MATRIX, IMBA
         INFO = 1
      ELSE IF( INBA.LT.1 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9995 ) MATRIX, 'INB', MATRIX, INBA
         INFO = 1
      END IF
*
*     Verify if blocking factors (MB_, NB_) are correct
*
      IF( MBA.LT.1 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9994 ) MATRIX, 'MB', MATRIX, MBA
         INFO = 1
      ELSE IF( NBA.LT.1 ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 )
     $      WRITE( NOUT, FMT = 9993 ) MATRIX, 'NB', MATRIX, NBA
         INFO = 1
      END IF
*
*     Verify if origin process coordinates (RSRC_, CSRC_) are valid
*
      IF( RSRCA.LT.-1 .OR. RSRCA.GE.NPROW ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9992 ) MATRIX
            WRITE( NOUT, FMT = 9990 ) 'RSRC', MATRIX, RSRCA, NPROW
         END IF
         INFO = 1
      ELSE IF( CSRCA.LT.-1 .OR. CSRCA.GE.NPCOL ) THEN
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9991 ) MATRIX
            WRITE( NOUT, FMT = 9990 ) 'CSRC', MATRIX, CSRCA, NPCOL
         END IF
         INFO = 1
      END IF
*
*     Check all processes for an error
*
      CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, 0 )
*
      IF( INFO.NE.0 ) THEN
*
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            WRITE( NOUT, FMT = 9989 ) MATRIX
            WRITE( NOUT, FMT = * )
         END IF
*
      ELSE
*
*        Compute local testing leading dimension
*
         MPA    = PB_NUMROC( MA, 1, IMBA, MBA, MYROW, RSRCA, NPROW )
         NQA    = PB_NUMROC( NA, 1, INBA, NBA, MYCOL, CSRCA, NPCOL )
         IPREA  = MAX( GAPMUL*NBA, MPA )
         IMIDA  = IGAP
         IPOSTA = MAX( GAPMUL*NBA, NQA )
         LLDA   = MAX( 1, MPA ) + IMIDA
*
         CALL PB_DESCINIT2( DESCA, MA, NA, IMBA, INBA, MBA, NBA, RSRCA,
     $                      CSRCA, ICTXT, LLDA, INFO )
*
*        Check all processes for an error
*
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, 0 )
*
         IF( INFO.NE.0 ) THEN
            IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
               WRITE( NOUT, FMT = 9989 ) MATRIX
               WRITE( NOUT, FMT = * )
            END IF
         END IF
*
      END IF
*
 9999 FORMAT( 2X, '>> Invalid matrix ', A1, ' descriptor type ', A5, A1,
     $        ': ', I6, ' should be ', I3, '.' )
 9998 FORMAT( 2X, '>> Invalid matrix ', A1, ' row dimension ', A1, A1,
     $        ': ', I6, ' should be at least 1.' )
 9997 FORMAT( 2X, '>> Invalid matrix ', A1, ' column dimension ', A1,
     $        A1, ': ', I6, ' should be at least 1.' )
 9996 FORMAT( 2X, '>> Invalid matrix ', A1, ' first row block size ',
     $        A3, A1, ': ', I6, ' should be at least 1.' )
 9995 FORMAT( 2X, '>> Invalid matrix ', A1, ' first column block size ',
     $        A3, A1,': ', I6, ' should be at least 1.' )
 9994 FORMAT( 2X, '>> Invalid matrix ', A1, ' row block size ', A2, A1,
     $        ': ', I6, ' should be at least 1.' )
 9993 FORMAT( 2X, '>> Invalid matrix ', A1, ' column block size ', A2,
     $        A1,': ', I6, ' should be at least 1.' )
 9992 FORMAT( 2X, '>> Invalid matrix ', A1, ' row process source:' )
 9991 FORMAT( 2X, '>> Invalid matrix ', A1, ' column process source:' )
 9990 FORMAT( 2X, '>> ', A4, A1, '= ', I6, ' should be >= -1 and < ',
     $        I6, '.' )
 9989 FORMAT( 2X, '>> Invalid matrix ', A1, ' descriptor: going on to ',
     $        'next test case.' )
*
      RETURN
*
*     End of PMDESCCHK
*
      END
      SUBROUTINE PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, INFOT, NOUT
      CHARACTER*(*)      SNAME
*     ..
*
*  Purpose
*  =======
*
*  PCHKPBE  tests  whether a PBLAS routine has detected an error when it
*  should.  This routine does a global operation to ensure all processes
*  have detected this error.  If  an  error  has  been detected an error
*  message is displayed.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information  required to
*  establish the  mapping  between a  matrix entry and its corresponding
*  process and memory location.
*
*  In  the  following  comments,   the character _  should  be  read  as
*  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN       EXPLANATION
*  ---------------- --------------- ------------------------------------
*  DTYPE_A (global) DESCA( DTYPE_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT_  ) The BLACS context handle, indicating
*                                   the NPROW x NPCOL BLACS process grid
*                                   A  is distributed over.  The context
*                                   itself  is  global,  but  the handle
*                                   (the integer value) may vary.
*  M_A     (global) DESCA( M_     ) The  number of rows in the distribu-
*                                   ted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N_     ) The number of columns in the distri-
*                                   buted matrix A, N_A >= 0.
*  IMB_A   (global) DESCA( IMB_   ) The number of rows of the upper left
*                                   block of the matrix A, IMB_A > 0.
*  INB_A   (global) DESCA( INB_   ) The  number  of columns of the upper
*                                   left   block   of   the   matrix  A,
*                                   INB_A > 0.
*  MB_A    (global) DESCA( MB_    ) The blocking factor used to  distri-
*                                   bute the last  M_A-IMB_A rows of  A,
*                                   MB_A > 0.
*  NB_A    (global) DESCA( NB_    ) The blocking factor used to  distri-
*                                   bute the last  N_A-INB_A  columns of
*                                   A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC_  ) The process row over which the first
*                                   row of the matrix  A is distributed,
*                                   NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC_  ) The  process  column  over which the
*                                   first  column of  A  is distributed.
*                                   NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD_   ) The  leading  dimension of the local
*                                   array  storing  the  local blocks of
*                                   the distributed matrix A,
*                                   IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                   ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_NUMROC:
*  Lr( IA, K ) = PB_NUMROC( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_NUMROC( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  NOUT    (global input) INTEGER
*          On entry, NOUT specifies the unit number for the output file.
*          When NOUT is 6, output to screen,  when  NOUT is 0, output to
*          stderr. NOUT is only defined for process 0.
*
*  SNAME   (global input) CHARACTER*(*)
*          On entry, SNAME specifies the subroutine  name  calling  this
*          subprogram.
*
*  INFOT   (global input) INTEGER
*          On entry, INFOT specifies the position of the wrong argument.
*          If  the  PBLAS  error  handler is called, INFO will be set to
*          -INFOT.  This  routine  verifies if the error was reported by
*          all processes by doing a global sum, and assert the result to
*          be NPROW * NPCOL.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            GERR, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGSUM2D
*     ..
*     .. Common Blocks ..
      INTEGER            INFO, NBLOG
      COMMON             /INFOC/INFO, NBLOG
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      GERR = 0
      IF( INFO.NE.-INFOT )
     $   GERR = 1
*
      CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, GERR, 1, -1, 0 )
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         IF( GERR.EQ.( NPROW * NPCOL ) ) THEN
            WRITE( NOUT, FMT = 9999 ) SNAME, INFO, -INFOT
         END IF
      END IF
*
 9999 FORMAT( 1X, A7, ': *** ERROR *** ERROR CODE RETURNED = ', I6,
     $        ' SHOULD HAVE BEEN ', I6 )
*
      RETURN
*
*     End of PCHKPBE
*
      END
      REAL FUNCTION PSDIFF( X, Y )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      REAL               X, Y
*     ..
*
*  Purpose
*  =======
*
*  PSDIFF returns the scalar difference X - Y. Similarly to the
*  BLAS tester, this routine allows for the possibility of computing a
*  more accurate difference if necessary.
*
*  Arguments
*  =========
*
*  X       (input) REAL
*          The real scalar X.
*
*  Y       (input) REAL
*          The real scalar Y.
*
*  =====================================================================
*
*     .. Executable Statements ..
*
      PSDIFF = X - Y
*
      RETURN
*
*     End of PSDIFF
*
      END
*
      DOUBLE PRECISION FUNCTION PDDIFF( X, Y )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
*     ..
*
*  Purpose
*  =======
*
*  PDDIFF returns the scalar difference X - Y. Similarly to the
*  BLAS tester, this routine allows for the possibility of computing a
*  more accurate difference if necessary.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*          The real scalar X.
*
*  Y       (input) DOUBLE PRECISION
*          The real scalar Y.
*
*  =====================================================================
*
*     .. Executable Statements ..
*
      PDDIFF = X - Y
*
      RETURN
*
*     End of PDDIFF
*
      END
      SUBROUTINE PXERBLA( ICTXT, SRNAME, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, INFO
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      SRNAME
*     ..
*
*  Purpose
*  =======
*
*  PXERBLA is an error handler for the ScaLAPACK routines.  It is called
*  by a ScaLAPACK routine if an input parameter has an invalid value.  A
*  message is printed. Installers may consider modifying this routine in
*  order to call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  SRNAME  (global input) CHARACTER*(*)
*          On entry, SRNAME specifies the name of the routine which cal-
*          ling PXERBLA.
*
*  INFO    (global input) INTEGER
*          On entry, INFO  specifies the position of the invalid parame-
*          ter in the parameter list of the calling routine.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      WRITE( *, FMT = 9999 ) MYROW, MYCOL, SRNAME, INFO
*
 9999 FORMAT( '{', I5, ',', I5, '}:  On entry to ', A,
     $        ' parameter number ', I4, ' had an illegal value' )
*
      RETURN
*
*     End of PXERBLA
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 2.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      LOGICAL          FUNCTION LSAMEN( N, CA, CB )
*
*  -- LAPACK auxiliary routine (version 2.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    CA, CB
      INTEGER            N
*     ..
*
*  Purpose
*  =======
*
*  LSAMEN  tests if the first N letters of CA are the same as the
*  first N letters of CB, regardless of case.
*  LSAMEN returns .TRUE. if CA and CB are equivalent except for case
*  and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA )
*  or LEN( CB ) is less than N.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of characters in CA and CB to be compared.
*
*  CA      (input) CHARACTER*(*)
*  CB      (input) CHARACTER*(*)
*          CA and CB specify two character strings of length at least N.
*          Only the first N characters of each string will be accessed.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN
*     ..
*     .. Executable Statements ..
*
      LSAMEN = .FALSE.
      IF( LEN( CA ).LT.N .OR. LEN( CB ).LT.N )
     $   GO TO 20
*
*     Do for each character in the two strings.
*
      DO 10 I = 1, N
*
*        Test if the characters are equal using LSAME.
*
         IF( .NOT.LSAME( CA( I: I ), CB( I: I ) ) )
     $      GO TO 20
*
   10 CONTINUE
      LSAMEN = .TRUE.
*
   20 CONTINUE
      RETURN
*
*     End of LSAMEN
*
      END
      SUBROUTINE ICOPY( N, SX, INCX, SY, INCY )
*
*  -- LAPACK auxiliary test routine (version 2.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
*     ..
*     .. Array Arguments ..
      INTEGER            SX( * ), SY( * )
*     ..
*
*  Purpose
*  =======
*
*  ICOPY copies an integer vector x to an integer vector y.
*  Uses unrolled loops for increments equal to 1.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The length of the vectors SX and SY.
*
*  SX      (input) INTEGER array, dimension (1+(N-1)*abs(INCX))
*          The vector X.
*
*  INCX    (input) INTEGER
*          The spacing between consecutive elements of SX.
*
*  SY      (output) INTEGER array, dimension (1+(N-1)*abs(INCY))
*          The vector Y.
*
*  INCY    (input) INTEGER
*          The spacing between consecutive elements of SY.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IX, IY, M, MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 )
     $   GO TO 20
*
*     Code for unequal increments or equal increments not equal to 1
*
      IX = 1
      IY = 1
      IF( INCX.LT.0 )
     $   IX = ( -N+1 )*INCX + 1
      IF( INCY.LT.0 )
     $   IY = ( -N+1 )*INCY + 1
      DO 10 I = 1, N
         SY( IY ) = SX( IX )
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
*
*     Code for both increments equal to 1
*
*     Clean-up loop
*
   20 CONTINUE
      M = MOD( N, 7 )
      IF( M.EQ.0 )
     $   GO TO 40
      DO 30 I = 1, M
         SY( I ) = SX( I )
   30 CONTINUE
      IF( N.LT.7 )
     $   RETURN
   40 CONTINUE
      MP1 = M + 1
      DO 50 I = MP1, N, 7
         SY( I ) = SX( I )
         SY( I+1 ) = SX( I+1 )
         SY( I+2 ) = SX( I+2 )
         SY( I+3 ) = SX( I+3 )
         SY( I+4 ) = SX( I+4 )
         SY( I+5 ) = SX( I+5 )
         SY( I+6 ) = SX( I+6 )
   50 CONTINUE
      RETURN
*
*     End of ICOPY
*
      END
      INTEGER FUNCTION PB_NOABORT( CINFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            CINFO
*     ..
*
*  Purpose
*  =======
*
*  PB_NOABORT  transmits  the  info  parameter of a PBLAS routine to the
*  tester  and  tells the PBLAS error handler to avoid aborting on erro-
*  neous input arguments.
*
*  Notes
*  =====
*
*  This  routine  is  necessary  because of the CRAY C fortran interface
*  and  the  fact  that  the  usual PBLAS error handler routine has been
*  initially written in C.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Common Blocks ..
      INTEGER            INFO, NBLOG, NOUT
      LOGICAL            ABRTFLG
      COMMON             /INFOC/INFO, NBLOG
      COMMON             /PBERRORC/NOUT, ABRTFLG
*     ..
*     .. Executable Statements ..
*
      INFO = CINFO
      IF( ABRTFLG ) THEN
         PB_NOABORT = 0
      ELSE
         PB_NOABORT = 1
      END IF
*
      RETURN
*
*     End of PB_NOABORT
*
      END
      SUBROUTINE PB_INFOG2L( I, J, DESC, NPROW, NPCOL, MYROW, MYCOL, II,
     $                       JJ, PROW, PCOL )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            I, II, J, JJ, MYCOL, MYROW, NPCOL, NPROW, PCOL,
     $                   PROW
*     ..
*     .. Array Arguments ..
      INTEGER            DESC( * )
*     ..
*
*  Purpose
*  =======
*
*  PB_INFOG2L  computes the starting local index II, JJ corresponding to
*  the submatrix starting globally at the entry pointed by  I,  J.  This
*  routine returns the coordinates in the grid of the process owning the
*  matrix entry of global indexes I, J, namely PROW and PCOL.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information  required to
*  establish the  mapping  between a  matrix entry and its corresponding
*  process and memory location.
*
*  In  the  following  comments,   the character _  should  be  read  as
*  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN       EXPLANATION
*  ---------------- --------------- ------------------------------------
*  DTYPE_A (global) DESCA( DTYPE_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT_  ) The BLACS context handle, indicating
*                                   the NPROW x NPCOL BLACS process grid
*                                   A  is distributed over.  The context
*                                   itself  is  global,  but  the handle
*                                   (the integer value) may vary.
*  M_A     (global) DESCA( M_     ) The  number of rows in the distribu-
*                                   ted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N_     ) The number of columns in the distri-
*                                   buted matrix A, N_A >= 0.
*  IMB_A   (global) DESCA( IMB_   ) The number of rows of the upper left
*                                   block of the matrix A, IMB_A > 0.
*  INB_A   (global) DESCA( INB_   ) The  number  of columns of the upper
*                                   left   block   of   the   matrix  A,
*                                   INB_A > 0.
*  MB_A    (global) DESCA( MB_    ) The blocking factor used to  distri-
*                                   bute the last  M_A-IMB_A rows of  A,
*                                   MB_A > 0.
*  NB_A    (global) DESCA( NB_    ) The blocking factor used to  distri-
*                                   bute the last  N_A-INB_A  columns of
*                                   A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC_  ) The process row over which the first
*                                   row of the matrix  A is distributed,
*                                   NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC_  ) The  process  column  over which the
*                                   first  column of  A  is distributed.
*                                   NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD_   ) The  leading  dimension of the local
*                                   array  storing  the  local blocks of
*                                   the distributed matrix A,
*                                   IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                   ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_NUMROC:
*  Lr( IA, K ) = PB_NUMROC( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_NUMROC( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  Arguments
*  =========
*
*  I       (global input) INTEGER
*          On entry, I  specifies  the  global starting row index of the
*          submatrix. I must at least one.
*
*  J       (global input) INTEGER
*          On entry, J  specifies  the  global  starting column index of
*          the submatrix. J must at least one.
*
*  DESC    (global and local input) INTEGER array
*          On entry,  DESC is an integer array of dimension DLEN_.  This
*          is the array descriptor of the underlying matrix.
*
*  NPROW   (global input) INTEGER
*          On entry,  NPROW   specifies the total number of process rows
*          over which the matrix is distributed.  NPROW must be at least
*          one.
*
*  NPCOL   (global input) INTEGER
*          On entry, NPCOL specifies the total number of process columns
*          over which the matrix is distributed.  NPCOL must be at least
*          one.
*
*  MYROW   (local input) INTEGER
*          On entry,  MYROW  specifies the row coordinate of the process
*          whose local index  II  is determined.  MYROW must be at least
*          zero and strictly less than NPROW.
*
*  MYCOL   (local input) INTEGER
*          On entry,  MYCOL  specifies the column coordinate of the pro-
*          cess whose local index  JJ  is determined.  MYCOL  must be at
*          least zero and strictly less than NPCOL.
*
*  II      (local output) INTEGER
*          On exit, II  specifies the  local  starting  row index of the
*          submatrix. On exit, II is at least one.
*
*  JJ      (local output) INTEGER
*          On exit, JJ  specifies the local starting column index of the
*          submatrix. On exit, JJ is at least one.
*
*  PROW    (global output) INTEGER
*          On exit,  PROW  specifies  the  row coordinate of the process
*          that possesses the first row of the submatrix.  On exit, PROW
*          is -1 if DESC( RSRC_ )  is -1 on input, and,  at  least  zero
*          and strictly less than NPROW otherwise.
*
*  PCOL    (global output) INTEGER
*          On exit, PCOL  specifies the column coordinate of the process
*          that possesses the first column of the  submatrix.  On  exit,
*          PCOL is -1 if DESC( CSRC_ )  is -1 on input, and,  at   least
*          zero and strictly less than NPCOL otherwise.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      INTEGER            CSRC, I1, ILOCBLK, IMB, INB, J1, MB, MYDIST,
     $                   NB, NBLOCKS, RSRC
*     ..
*     .. Local Arrays ..
      INTEGER            DESC2( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           PB_DESCTRANS
*     ..
*     .. Executable Statements ..
*
*     Convert descriptor
*
      CALL PB_DESCTRANS( DESC, DESC2 )
*
      IMB  = DESC2( IMB_ )
      PROW = DESC2( RSRC_ )
*
*     Has every process row I ?
*
      IF( ( PROW.EQ.-1 ).OR.( NPROW.EQ.1 ) ) THEN
*
         II = I
*
      ELSE IF( I.LE.IMB ) THEN
*
*        I is in range of first block
*
         IF( MYROW.EQ.PROW ) THEN
            II = I
         ELSE
            II = 1
         END IF
*
      ELSE
*
*        I is not in first block of matrix, figure out who has it.
*
         RSRC = PROW
         MB = DESC2( MB_ )
*
         IF( MYROW.EQ.RSRC ) THEN
*
            NBLOCKS = ( I - IMB - 1 ) / MB + 1
            PROW    = PROW + NBLOCKS
            PROW    = PROW - ( PROW / NPROW ) * NPROW
*
            ILOCBLK = NBLOCKS / NPROW
*
            IF( ILOCBLK.GT.0 ) THEN
               IF( ( ILOCBLK*NPROW ).GE.NBLOCKS ) THEN
                  IF( MYROW.EQ.PROW ) THEN
                     II = I + ( ILOCBLK - NBLOCKS ) * MB
                  ELSE
                     II = IMB + ( ILOCBLK - 1 ) * MB + 1
                  END IF
               ELSE
                  II = IMB + ILOCBLK * MB + 1
               END IF
            ELSE
               II = IMB + 1
            END IF
*
         ELSE
*
            I1      = I - IMB
            NBLOCKS = ( I1 - 1 ) / MB + 1
            PROW    = PROW + NBLOCKS
            PROW    = PROW - ( PROW / NPROW ) * NPROW
*
            MYDIST  = MYROW - RSRC
            IF( MYDIST.LT.0 )
     $         MYDIST = MYDIST + NPROW
*
            ILOCBLK = NBLOCKS / NPROW
*
            IF( ILOCBLK.GT.0 ) THEN
               MYDIST = MYDIST - NBLOCKS + ILOCBLK * NPROW
               IF( MYDIST.LT.0 ) THEN
                  II = MB + ILOCBLK * MB + 1
               ELSE
                  IF( MYROW.EQ.PROW ) THEN
                     II = I1 + ( ILOCBLK - NBLOCKS + 1 ) * MB
                  ELSE
                     II = ILOCBLK * MB + 1
                  END IF
               END IF
            ELSE
               MYDIST = MYDIST - NBLOCKS
               IF( MYDIST.LT.0 ) THEN
                  II = MB + 1
               ELSE IF( MYROW.EQ.PROW ) THEN
                  II = I1 + ( 1 - NBLOCKS ) * MB
               ELSE
                  II = 1
               END IF
            END IF
         END IF
*
      END IF
*
      INB  = DESC2( INB_ )
      PCOL = DESC2( CSRC_ )
*
*     Has every process column J ?
*
      IF( ( PCOL.EQ.-1 ).OR.( NPCOL.EQ.1 ) ) THEN
*
         JJ = J
*
      ELSE IF( J.LE.INB ) THEN
*
*        J is in range of first block
*
         IF( MYCOL.EQ.PCOL ) THEN
            JJ = J
         ELSE
            JJ = 1
         END IF
*
      ELSE
*
*        J is not in first block of matrix, figure out who has it.
*
         CSRC = PCOL
         NB   = DESC2( NB_ )
*
         IF( MYCOL.EQ.CSRC ) THEN
*
            NBLOCKS = ( J - INB - 1 ) / NB + 1
            PCOL    = PCOL + NBLOCKS
            PCOL    = PCOL - ( PCOL / NPCOL ) * NPCOL
*
            ILOCBLK = NBLOCKS / NPCOL
*
            IF( ILOCBLK.GT.0 ) THEN
               IF( ( ILOCBLK*NPCOL ).GE.NBLOCKS ) THEN
                  IF( MYCOL.EQ.PCOL ) THEN
                     JJ = J + ( ILOCBLK - NBLOCKS ) * NB
                  ELSE
                     JJ = INB + ( ILOCBLK - 1 ) * NB + 1
                  END IF
               ELSE
                  JJ = INB + ILOCBLK * NB + 1
               END IF
            ELSE
               JJ = INB + 1
            END IF
*
         ELSE
*
            J1      = J - INB
            NBLOCKS = ( J1 - 1 ) / NB + 1
            PCOL    = PCOL + NBLOCKS
            PCOL    = PCOL - ( PCOL / NPCOL ) * NPCOL
*
            MYDIST  = MYCOL - CSRC
            IF( MYDIST.LT.0 )
     $         MYDIST = MYDIST + NPCOL
*
            ILOCBLK = NBLOCKS / NPCOL
*
            IF( ILOCBLK.GT.0 ) THEN
               MYDIST = MYDIST - NBLOCKS + ILOCBLK * NPCOL
               IF( MYDIST.LT.0 ) THEN
                  JJ = NB + ILOCBLK * NB + 1
               ELSE
                  IF( MYCOL.EQ.PCOL ) THEN
                     JJ = J1 + ( ILOCBLK - NBLOCKS + 1 ) * NB
                  ELSE
                     JJ = ILOCBLK * NB + 1
                  END IF
               END IF
            ELSE
               MYDIST = MYDIST - NBLOCKS
               IF( MYDIST.LT.0 ) THEN
                  JJ = NB + 1
               ELSE IF( MYCOL.EQ.PCOL ) THEN
                  JJ = J1 + ( 1 - NBLOCKS ) * NB
               ELSE
                  JJ = 1
               END IF
            END IF
         END IF
*
      END IF
*
      RETURN
*
*     End of PB_INFOG2L
*
      END
      SUBROUTINE PB_AINFOG2L( M, N, I, J, DESC, NPROW, NPCOL, MYROW,
     $                        MYCOL, IMB1, INB1, MP, NQ, II, JJ, PROW,
     $                        PCOL, RPROW, RPCOL )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            I, II, IMB1, INB1, J, JJ, M, MP, MYCOL, MYROW,
     $                   N, NPCOL, NPROW, NQ, PCOL, PROW, RPCOL, RPROW
*     ..
*     .. Array Arguments ..
      INTEGER            DESC( * )
*     ..
*
*  Purpose
*  =======
*
*  PB_AINFOG2L  computes the  starting  local row and column indexes II,
*  JJ  corresponding to  the  submatrix  starting  globally at the entry
*  pointed by I,  J. This routine returns the coordinates in the grid of
*  the  process owning  the  matrix entry of global indexes I, J, namely
*  PROW  and  PCOL. In addition, this routine computes the quantities MP
*  and  NQ,  which are respectively the local number of rows and columns
*  owned by the process of coordinate  MYROW, MYCOL corresponding to the
*  global submatrix A(I:I+M-1,J:J+N-1).  Finally, the size  of the first
*  partial block and the relative process coordinates  are also returned
*  respectively in IMB, INB and RPROW, RPCOL.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information  required to
*  establish the  mapping  between a  matrix entry and its corresponding
*  process and memory location.
*
*  In  the  following  comments,   the character _  should  be  read  as
*  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN       EXPLANATION
*  ---------------- --------------- ------------------------------------
*  DTYPE_A (global) DESCA( DTYPE_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT_  ) The BLACS context handle, indicating
*                                   the NPROW x NPCOL BLACS process grid
*                                   A  is distributed over.  The context
*                                   itself  is  global,  but  the handle
*                                   (the integer value) may vary.
*  M_A     (global) DESCA( M_     ) The  number of rows in the distribu-
*                                   ted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N_     ) The number of columns in the distri-
*                                   buted matrix A, N_A >= 0.
*  IMB_A   (global) DESCA( IMB_   ) The number of rows of the upper left
*                                   block of the matrix A, IMB_A > 0.
*  INB_A   (global) DESCA( INB_   ) The  number  of columns of the upper
*                                   left   block   of   the   matrix  A,
*                                   INB_A > 0.
*  MB_A    (global) DESCA( MB_    ) The blocking factor used to  distri-
*                                   bute the last  M_A-IMB_A rows of  A,
*                                   MB_A > 0.
*  NB_A    (global) DESCA( NB_    ) The blocking factor used to  distri-
*                                   bute the last  N_A-INB_A  columns of
*                                   A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC_  ) The process row over which the first
*                                   row of the matrix  A is distributed,
*                                   NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC_  ) The  process  column  over which the
*                                   first  column of  A  is distributed.
*                                   NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD_   ) The  leading  dimension of the local
*                                   array  storing  the  local blocks of
*                                   the distributed matrix A,
*                                   IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                   ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_NUMROC:
*  Lr( IA, K ) = PB_NUMROC( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_NUMROC( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  Arguments
*  =========
*
*  M       (global input) INTEGER
*          On entry, M specifies the global number of rows of the subma-
*          trix. M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N specifies  the  global  number  of columns of the
*          submatrix. N must be at least zero.
*
*  I       (global input) INTEGER
*          On entry, I  specifies  the  global starting row index of the
*          submatrix. I must at least one.
*
*  J       (global input) INTEGER
*          On entry, J  specifies  the global starting column  index  of
*          the submatrix. J must at least one.
*
*  DESC    (global and local input) INTEGER array
*          On entry,  DESC is an integer array of dimension DLEN_.  This
*          is the array descriptor of the underlying matrix.
*
*  NPROW   (global input) INTEGER
*          On entry,  NPROW   specifies the total number of process rows
*          over which the matrix is distributed.  NPROW must be at least
*          one.
*
*  NPCOL   (global input) INTEGER
*          On entry, NPCOL specifies the total number of process columns
*          over which the matrix is distributed.  NPCOL must be at least
*          one.
*
*  MYROW   (local input) INTEGER
*          On entry,  MYROW  specifies the row coordinate of the process
*          whose local index  II  is determined.  MYROW must be at least
*          zero and strictly less than NPROW.
*
*  MYCOL   (local input) INTEGER
*          On entry,  MYCOL  specifies the column coordinate of the pro-
*          cess whose local index  JJ  is determined.  MYCOL  must be at
*          least zero and strictly less than NPCOL.
*
*  IMB1    (global output) INTEGER
*          On exit, IMB1 specifies the number of rows of the upper  left
*          block of the submatrix. On exit,  IMB1 is less or equal  than
*          M and greater or equal than MIN( 1, M ).
*
*  INB1    (global output) INTEGER
*          On exit, INB1 specifies  the number  of  columns of the upper
*          left block of the submatrix. On exit,  INB1 is  less or equal
*          than N and greater or equal than MIN( 1, N ).
*
*  MP      (local output) INTEGER
*          On exit, MP specifies the local number of rows of the  subma-
*          trix, that the processes of row coordinate MYROW own.  MP  is
*          at least zero.
*
*  NQ      (local output) INTEGER
*          On exit, NQ specifies  the  local  number  of columns  of the
*          submatrix,  that  the processes  of column  coordinate  MYCOL
*          own. NQ is at least zero.
*
*  II      (local output) INTEGER
*          On exit, II  specifies the  local  starting  row index of the
*          submatrix. On exit, II is at least one.
*
*  JJ      (local output) INTEGER
*          On exit, JJ  specifies the  local  starting  column index  of
*          the submatrix. On exit, II is at least one.
*
*  PROW    (global output) INTEGER
*          On exit,  PROW  specifies the row coordinate of  the  process
*          that possesses the first row of the submatrix. On exit,  PROW
*          is -1 if DESC(RSRC_)  is -1 on input, and, at least zero  and
*          strictly less than NPROW otherwise.
*
*  PCOL    (global output) INTEGER
*          On exit, PCOL  specifies the column coordinate of the process
*          that possesses the first column of the  submatrix.  On  exit,
*          PCOL is -1 if DESC(CSRC_)  is -1 on input, and, at least zero
*          and strictly less than NPCOL otherwise.
*
*  RPROW   (global output) INTEGER
*          On exit, RPROW specifies  the  relative row coordinate of the
*          process that possesses the first row  I  of the submatrix. On
*          exit, RPROW is -1 if DESC(RSRC_) is  -1  on  input,  and,  at
*          least zero and strictly less than NPROW otherwise.
*
*  RPCOL   (global output) INTEGER
*          On exit, RPCOL specifies  the  relative column  coordinate of
*          the process that possesses the first column  J  of the subma-
*          trix. On exit, RPCOL is -1 if  DESC(CSRC_)  is  -1  on input,
*          and, at least zero and strictly less than NPCOL otherwise.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      INTEGER            CSRC, I1, ILOCBLK, J1, M1, MB, MYDIST, N1, NB,
     $                   NBLOCKS, RSRC
*     ..
*     .. Local Arrays ..
      INTEGER            DESC2( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           PB_DESCTRANS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Convert descriptor
*
      CALL PB_DESCTRANS( DESC, DESC2 )
*
      MB   = DESC2( MB_ )
      IMB1 = DESC2( IMB_ )
      RSRC = DESC2( RSRC_ )
*
      IF( ( RSRC.EQ.-1 ).OR.( NPROW.EQ.1 ) ) THEN
*
         II    = I
         IMB1  = IMB1 - I + 1
         IF( IMB1.LE.0 )
     $      IMB1 = ( ( -IMB1 ) / MB + 1 ) * MB + IMB1
         IMB1  = MIN( IMB1, M )
         MP    = M
         PROW  = RSRC
         RPROW = 0
*
      ELSE
*
*        Figure out PROW, II and IMB1 first
*
         IF( I.LE.IMB1 ) THEN
*
            PROW = RSRC
*
            IF( MYROW.EQ.PROW ) THEN
               II = I
            ELSE
               II = 1
            END IF
*
            IMB1 = IMB1 - I + 1
*
         ELSE
*
            I1 = I - IMB1 - 1
            NBLOCKS = I1 / MB + 1
            PROW = RSRC + NBLOCKS
            PROW = PROW - ( PROW / NPROW ) * NPROW
*
            IF( MYROW.EQ.RSRC ) THEN
*
               ILOCBLK = NBLOCKS / NPROW
*
               IF( ILOCBLK.GT.0 ) THEN
                  IF( ( ILOCBLK*NPROW ).GE.NBLOCKS ) THEN
                     IF( MYROW.EQ.PROW ) THEN
                        II = I + ( ILOCBLK - NBLOCKS ) * MB
                     ELSE
                        II = IMB1 + ( ILOCBLK - 1 ) * MB + 1
                     END IF
                  ELSE
                     II = IMB1 + ILOCBLK * MB + 1
                  END IF
               ELSE
                  II = IMB1 + 1
               END IF
*
            ELSE
*
               MYDIST = MYROW - RSRC
               IF( MYDIST.LT.0 )
     $            MYDIST = MYDIST + NPROW
*
               ILOCBLK = NBLOCKS / NPROW
*
               IF( ILOCBLK.GT.0 ) THEN
                  MYDIST = MYDIST - NBLOCKS + ILOCBLK * NPROW
                  IF( MYDIST.LT.0 ) THEN
                     II = ( ILOCBLK + 1 ) * MB + 1
                  ELSE IF( MYROW.EQ.PROW ) THEN
                     II = I1 + ( ILOCBLK - NBLOCKS + 1 ) * MB + 1
                  ELSE
                     II = ILOCBLK * MB + 1
                  END IF
               ELSE
                  MYDIST = MYDIST - NBLOCKS
                  IF( MYDIST.LT.0 ) THEN
                     II = MB + 1
                  ELSE IF( MYROW.EQ.PROW ) THEN
                     II = I1 + ( 1 - NBLOCKS ) * MB + 1
                  ELSE
                     II = 1
                  END IF
               END IF
            END IF
*
            IMB1 = NBLOCKS * MB - I1
*
         END IF
*
*        Figure out MP
*
         IF( M.LE.IMB1 ) THEN
*
            IF( MYROW.EQ.PROW ) THEN
               MP = M
            ELSE
               MP = 0
            END IF
*
         ELSE
*
            M1 = M - IMB1
            NBLOCKS = M1 / MB + 1
*
            IF( MYROW.EQ.PROW ) THEN
               ILOCBLK = NBLOCKS / NPROW
               IF( ILOCBLK.GT.0 ) THEN
                  IF( ( NBLOCKS - ILOCBLK * NPROW ).GT.0 ) THEN
                     MP = IMB1 + ILOCBLK * MB
                  ELSE
                     MP = M + MB * ( ILOCBLK - NBLOCKS )
                  END IF
               ELSE
                  MP = IMB1
               END IF
            ELSE
               MYDIST = MYROW - PROW
               IF( MYDIST.LT.0 )
     $            MYDIST = MYDIST + NPROW
               ILOCBLK = NBLOCKS / NPROW
               IF( ILOCBLK.GT.0 ) THEN
                  MYDIST = MYDIST - NBLOCKS + ILOCBLK * NPROW
                  IF( MYDIST.LT.0 ) THEN
                     MP = ( ILOCBLK + 1 ) * MB
                  ELSE IF( MYDIST.GT.0 ) THEN
                     MP = ILOCBLK * MB
                  ELSE
                     MP = M1 + MB * ( ILOCBLK - NBLOCKS + 1 )
                  END IF
               ELSE
                  MYDIST = MYDIST - NBLOCKS
                  IF( MYDIST.LT.0 ) THEN
                     MP = MB
                  ELSE IF( MYDIST.GT.0 ) THEN
                     MP = 0
                  ELSE
                     MP = M1 + MB * ( 1 - NBLOCKS )
                  END IF
               END IF
            END IF
*
         END IF
*
         IMB1 = MIN( IMB1, M )
         RPROW = MYROW - PROW
         IF( RPROW.LT.0 )
     $      RPROW = RPROW + NPROW
*
      END IF
*
      NB   = DESC2( NB_ )
      INB1 = DESC2( INB_ )
      CSRC = DESC2( CSRC_ )
*
      IF( ( CSRC.EQ.-1 ).OR.( NPCOL.EQ.1 ) ) THEN
*
         JJ    = J
         INB1  = INB1 - I + 1
         IF( INB1.LE.0 )
     $      INB1 = ( ( -INB1 ) / NB + 1 ) * NB + INB1
         INB1  = MIN( INB1, N )
         NQ    = N
         PCOL  = CSRC
         RPCOL = 0
*
      ELSE
*
*        Figure out PCOL, JJ and INB1 first
*
         IF( J.LE.INB1 ) THEN
*
            PCOL = CSRC
*
            IF( MYCOL.EQ.PCOL ) THEN
               JJ = J
            ELSE
               JJ = 1
            END IF
*
            INB1 = INB1 - J + 1
*
         ELSE
*
            J1 = J - INB1 - 1
            NBLOCKS = J1 / NB + 1
            PCOL = CSRC + NBLOCKS
            PCOL = PCOL - ( PCOL / NPCOL ) * NPCOL
*
            IF( MYCOL.EQ.CSRC ) THEN
*
               ILOCBLK = NBLOCKS / NPCOL
*
               IF( ILOCBLK.GT.0 ) THEN
                  IF( ( ILOCBLK*NPCOL ).GE.NBLOCKS ) THEN
                     IF( MYCOL.EQ.PCOL ) THEN
                        JJ = J + ( ILOCBLK - NBLOCKS ) * NB
                     ELSE
                        JJ = INB1 + ( ILOCBLK - 1 ) * NB + 1
                     END IF
                  ELSE
                     JJ = INB1 + ILOCBLK * NB + 1
                  END IF
               ELSE
                  JJ = INB1 + 1
               END IF
*
            ELSE
*
               MYDIST = MYCOL - CSRC
               IF( MYDIST.LT.0 )
     $            MYDIST = MYDIST + NPCOL
*
               ILOCBLK = NBLOCKS / NPCOL
*
               IF( ILOCBLK.GT.0 ) THEN
                  MYDIST = MYDIST - NBLOCKS + ILOCBLK * NPCOL
                  IF( MYDIST.LT.0 ) THEN
                     JJ = ( ILOCBLK + 1 ) * NB + 1
                  ELSE IF( MYCOL.EQ.PCOL ) THEN
                     JJ = J1 + ( ILOCBLK - NBLOCKS + 1 ) * NB + 1
                  ELSE
                     JJ = ILOCBLK * NB + 1
                  END IF
               ELSE
                  MYDIST = MYDIST - NBLOCKS
                  IF( MYDIST.LT.0 ) THEN
                     JJ = NB + 1
                  ELSE IF( MYCOL.EQ.PCOL ) THEN
                     JJ = J1 + ( 1 - NBLOCKS ) * NB + 1
                  ELSE
                     JJ = 1
                  END IF
               END IF
            END IF
*
            INB1 = NBLOCKS * NB - J1
*
         END IF
*
*        Figure out NQ
*
         IF( N.LE.INB1 ) THEN
*
            IF( MYCOL.EQ.PCOL ) THEN
               NQ = N
            ELSE
               NQ = 0
            END IF
*
         ELSE
*
            N1 = N - INB1
            NBLOCKS = N1 / NB + 1
*
            IF( MYCOL.EQ.PCOL ) THEN
               ILOCBLK = NBLOCKS / NPCOL
               IF( ILOCBLK.GT.0 ) THEN
                  IF( ( NBLOCKS - ILOCBLK * NPCOL ).GT.0 ) THEN
                     NQ = INB1 + ILOCBLK * NB
                  ELSE
                     NQ = N + NB * ( ILOCBLK - NBLOCKS )
                  END IF
               ELSE
                  NQ = INB1
               END IF
            ELSE
               MYDIST = MYCOL - PCOL
               IF( MYDIST.LT.0 )
     $            MYDIST = MYDIST + NPCOL
               ILOCBLK = NBLOCKS / NPCOL
               IF( ILOCBLK.GT.0 ) THEN
                  MYDIST = MYDIST - NBLOCKS + ILOCBLK * NPCOL
                  IF( MYDIST.LT.0 ) THEN
                     NQ = ( ILOCBLK + 1 ) * NB
                  ELSE IF( MYDIST.GT.0 ) THEN
                     NQ = ILOCBLK * NB
                  ELSE
                     NQ = N1 + NB * ( ILOCBLK - NBLOCKS + 1 )
                  END IF
               ELSE
                  MYDIST = MYDIST - NBLOCKS
                  IF( MYDIST.LT.0 ) THEN
                     NQ = NB
                  ELSE IF( MYDIST.GT.0 ) THEN
                     NQ = 0
                  ELSE
                     NQ = N1 + NB * ( 1 - NBLOCKS )
                  END IF
               END IF
            END IF
*
         END IF
*
         INB1 = MIN( INB1, N )
         RPCOL = MYCOL - PCOL
         IF( RPCOL.LT.0 )
     $      RPCOL = RPCOL + NPCOL
*
      END IF
*
      RETURN
*
*     End of PB_AINFOG2L
*
      END
      INTEGER FUNCTION PB_NUMROC( N, I, INB, NB, PROC, SRCPROC, NPROCS )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            I, INB, N, NB, NPROCS, PROC, SRCPROC
*     ..
*
*  Purpose
*  =======
*
*  PB_NUMROC   returns  the  local number of matrix rows/columns process
*  PROC will get  if we give out N rows/columns starting from global in-
*  dex I.
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of rows/columns being dealt
*          out. N must be at least zero.
*
*  I       (global input) INTEGER
*          On entry, I  specifies the global index of the matrix  entry.
*          I must be at least one.
*
*  INB     (global input) INTEGER
*          On entry,  INB  specifies  the size of the first block of the
*          global matrix. INB must be at least one.
*
*  NB      (global input) INTEGER
*          On entry, NB specifies the size of the blocks used to  parti-
*          tion the matrix. NB must be at least one.
*
*  PROC    (local input) INTEGER
*          On entry, PROC specifies  the coordinate of the process whose
*          local portion is determined.  PROC must be at least zero  and
*          strictly less than NPROCS.
*
*  SRCPROC (global input) INTEGER
*          On entry,  SRCPROC  specifies  the coordinate of the  process
*          that possesses the  first row or column  of the matrix.  When
*          SRCPROC = -1, the data  is not  distributed  but  replicated,
*          otherwise  SRCPROC  must be at least zero and  strictly  less
*          than NPROCS.
*
*  NPROCS  (global input) INTEGER
*          On entry,  NPROCS  specifies the total number of process rows
*          or columns over which the matrix is distributed.  NPROCS must
*          be at least one.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I1, ILOCBLK, INB1, MYDIST, N1, NBLOCKS,
     $                   SRCPROC1
*     ..
*     .. Executable Statements ..
*
      IF( ( SRCPROC.EQ.-1 ).OR.( NPROCS.EQ.1 ) ) THEN
         PB_NUMROC = N
         RETURN
      END IF
*
*     Compute coordinate of process owning I and corresponding INB
*
      IF( I.LE.INB ) THEN
*
*        I is in range of first block, i.e SRCPROC owns I.
*
         SRCPROC1 = SRCPROC
         INB1 = INB - I + 1
*
      ELSE
*
*        I is not in first block of matrix, figure out who has it
*
         I1 = I - 1 - INB
         NBLOCKS = I1 / NB + 1
         SRCPROC1 = SRCPROC + NBLOCKS
         SRCPROC1 = SRCPROC1 - ( SRCPROC1 / NPROCS ) * NPROCS
         INB1 = NBLOCKS*NB - I1
*
      END IF
*
*     Now everything is just like I=1. Search now who has N-1, Is N-1
*     in the first block ?
*
      IF( N.LE.INB1 ) THEN
         IF( PROC.EQ.SRCPROC1 ) THEN
            PB_NUMROC = N
         ELSE
            PB_NUMROC = 0
         END IF
         RETURN
      END IF
*
      N1 = N - INB1
      NBLOCKS = N1 / NB + 1
*
      IF( PROC.EQ.SRCPROC1 ) THEN
         ILOCBLK = NBLOCKS / NPROCS
         IF( ILOCBLK.GT.0 ) THEN
            IF( ( NBLOCKS - ILOCBLK * NPROCS ).GT.0 ) THEN
               PB_NUMROC = INB1 + ILOCBLK * NB
            ELSE
               PB_NUMROC = N + NB * ( ILOCBLK - NBLOCKS )
            END IF
         ELSE
            PB_NUMROC = INB1
         END IF
      ELSE
         MYDIST = PROC - SRCPROC1
         IF( MYDIST.LT.0 )
     $      MYDIST = MYDIST + NPROCS
         ILOCBLK = NBLOCKS / NPROCS
         IF( ILOCBLK.GT.0 ) THEN
            MYDIST = MYDIST - NBLOCKS + ILOCBLK * NPROCS
            IF( MYDIST.LT.0 ) THEN
               PB_NUMROC = ( ILOCBLK + 1 ) * NB
            ELSE IF( MYDIST.GT.0 ) THEN
               PB_NUMROC = ILOCBLK * NB
            ELSE
               PB_NUMROC = N1 + NB * ( ILOCBLK - NBLOCKS + 1 )
            END IF
         ELSE
            MYDIST = MYDIST - NBLOCKS
            IF( MYDIST.LT.0 ) THEN
               PB_NUMROC = NB
            ELSE IF( MYDIST.GT.0 ) THEN
               PB_NUMROC = 0
            ELSE
               PB_NUMROC = N1 + NB * ( 1 - NBLOCKS )
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of PB_NUMROC
*
      END
      INTEGER FUNCTION PB_FCEIL( NUM, DENOM )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      REAL               DENOM, NUM
*     ..
*
*  Purpose
*  =======
*
*  PB_FCEIL  returns  the  ceiling  of the division of two integers. The
*  integer operands are passed as real to avoid integer overflow.
*
*  Arguments
*  =========
*
*  NUM     (local input) REAL
*          On entry, NUM  specifies  the numerator of the fraction to be
*          evaluated.
*
*  DENOM   (local input) REAL
*          On entry, DENOM specifies  the denominator of the fraction to
*          be evaluated.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          NINT
*     ..
*     .. Executable Statements ..
*
      PB_FCEIL = NINT( ( ( NUM + DENOM - 1.0E+0 ) / DENOM ) - 0.5E+0 )
*
      RETURN
*
*     End of PB_FCEIL
*
      END
      SUBROUTINE PB_CHKMAT( ICTXT, M, MPOS0, N, NPOS0, IA, JA, DESCA,
     $                      DPOS0, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            DPOS0, IA, ICTXT, INFO, JA, M, MPOS0, N, NPOS0
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
*     ..
*
*  Purpose
*  =======
*
*  PB_CHKMAT  checks the validity of a descriptor vector  DESCA, the re-
*  lated global indexes  IA, JA from a local view point. If an inconsis-
*  tency is found among its parameters IA, JA and DESCA, the routine re-
*  turns an error code in INFO.
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies  the  number  of  rows  the submatrix
*          sub( A ).
*
*  MPOS0   (global input) INTEGER
*          On entry,  MPOS0  specifies the  position in the calling rou-
*          tine's parameter list where the formal parameter M appears.
*
*  N       (global input) INTEGER
*          On entry,  N  specifies  the  number of columns the submatrix
*          sub( A ).
*
*  NPOS0   (global input) INTEGER
*          On entry,  NPOS0  specifies the  position in the calling rou-
*          tine's parameter list where the formal parameter N appears.
*
*  IA      (global input) INTEGER
*          On entry, IA  specifies A's global row index, which points to
*          the beginning of the submatrix sub( A ).
*
*  JA      (global input) INTEGER
*          On entry, JA  specifies A's global column index, which points
*          to the beginning of the submatrix sub( A ).
*
*  DESCA   (global and local input) INTEGER array
*          On entry, DESCA  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix A.
*
*  DPOS0   (global input) INTEGER
*          On entry,  DPOS0  specifies the  position in the calling rou-
*          tine's parameter list where the formal  parameter  DESCA  ap-
*          pears.  Note that it is assumed that  IA and JA are respecti-
*          vely 2 and 1 entries behind DESCA.
*
*  INFO    (local input/local output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had an
*                illegal  value,  then  INFO = -(i*100+j),  if  the i-th
*                argument is a  scalar  and had an  illegal  value, then
*                INFO = -i.
*
*  -- Written on April 1, 1998 by
*     R. Clint Whaley, University of Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
      INTEGER            DESCMULT, BIGNUM
      PARAMETER          ( DESCMULT = 100, BIGNUM = DESCMULT*DESCMULT )
*     ..
*     .. Local Scalars ..
      INTEGER            DPOS, IAPOS, JAPOS, MP, MPOS, MYCOL, MYROW,
     $                   NPCOL, NPOS, NPROW, NQ
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA2( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PB_DESCTRANS
*     ..
*     .. External Functions ..
      INTEGER            PB_NUMROC
      EXTERNAL           PB_NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX
*     ..
*     .. Executable Statements ..
*
*     Convert descriptor
*
      CALL PB_DESCTRANS( DESCA, DESCA2 )
*
*     Want to find errors with MIN( ), so if no error, set it to a big
*     number.  If there already is an error, multiply by the the des-
*     criptor multiplier
*
      IF( INFO.GE.0 ) THEN
         INFO = BIGNUM
      ELSE IF( INFO.LT.-DESCMULT ) THEN
         INFO = -INFO
      ELSE
         INFO = -INFO * DESCMULT
      END IF
*
*     Figure where in parameter list each parameter was, factoring in
*     descriptor multiplier
*
      MPOS  = MPOS0 * DESCMULT
      NPOS  = NPOS0 * DESCMULT
      IAPOS = ( DPOS0 - 2 ) * DESCMULT
      JAPOS = ( DPOS0 - 1 ) * DESCMULT
      DPOS  = DPOS0 * DESCMULT
*
*     Get grid parameters
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Check that matrix values make sense from local viewpoint
*
      IF( M.LT.0 )
     $   INFO = MIN( INFO, MPOS )
      IF( N.LT.0 )
     $   INFO = MIN( INFO, NPOS )
      IF( IA.LT.1 )
     $   INFO = MIN( INFO, IAPOS )
      IF( JA.LT.1 )
     $   INFO = MIN( INFO, JAPOS )
      IF( DESCA2( DTYPE_ ).NE.BLOCK_CYCLIC_2D_INB )
     $   INFO = MIN( INFO, DPOS + DTYPE_ )
      IF( DESCA2( IMB_ ).LT.1 )
     $   INFO = MIN( INFO, DPOS + IMB_ )
      IF( DESCA2( INB_ ).LT.1 )
     $   INFO = MIN( INFO, DPOS + INB_ )
      IF( DESCA2( MB_ ).LT.1 )
     $   INFO = MIN( INFO, DPOS + MB_ )
      IF( DESCA2( NB_ ).LT.1 )
     $   INFO = MIN( INFO, DPOS + NB_ )
      IF( DESCA2( RSRC_ ).LT.-1 .OR. DESCA2( RSRC_ ).GE.NPROW )
     $   INFO = MIN( INFO, DPOS + RSRC_ )
      IF( DESCA2( CSRC_ ).LT.-1 .OR. DESCA2( CSRC_ ).GE.NPCOL )
     $   INFO = MIN( INFO, DPOS + CSRC_ )
      IF( DESCA2( CTXT_ ).NE.ICTXT )
     $   INFO = MIN( INFO, DPOS + CTXT_ )
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
*
*        NULL matrix, relax some checks
*
         IF( DESCA2( M_ ).LT.0 )
     $      INFO = MIN( INFO, DPOS + M_ )
         IF( DESCA2( N_ ).LT.0 )
     $      INFO = MIN( INFO, DPOS + N_ )
         IF( DESCA2( LLD_ ).LT.1 )
     $      INFO = MIN( INFO, DPOS + LLD_ )
*
      ELSE
*
*        more rigorous checks for non-degenerate matrices
*
         MP = PB_NUMROC( DESCA2( M_ ), 1, DESCA2( IMB_ ), DESCA2( MB_ ),
     $                   MYROW, DESCA2( RSRC_ ), NPROW )
*
         IF( DESCA2( M_ ).LT.1 )
     $      INFO = MIN( INFO, DPOS + M_ )
         IF( DESCA2( N_ ).LT.1 )
     $      INFO = MIN( INFO, DPOS + N_ )
         IF( IA.GT.DESCA2( M_ ) )
     $      INFO = MIN( INFO, IAPOS )
         IF( JA.GT.DESCA2( N_ ) )
     $      INFO = MIN( INFO, JAPOS )
         IF( IA+M-1.GT.DESCA2( M_ ) )
     $      INFO = MIN( INFO, MPOS )
         IF( JA+N-1.GT.DESCA2( N_ ) )
     $      INFO = MIN( INFO, NPOS )
*
         IF( DESCA2( LLD_ ).LT.MAX( 1, MP ) ) THEN
            NQ = PB_NUMROC( DESCA2( N_ ), 1, DESCA2( INB_ ),
     $                      DESCA2( NB_ ), MYCOL, DESCA2( CSRC_ ),
     $                      NPCOL )
            IF( DESCA2( LLD_ ).LT.1 ) THEN
               INFO = MIN( INFO, DPOS + LLD_ )
            ELSE IF( NQ.GT.0 ) THEN
               INFO = MIN( INFO, DPOS + LLD_ )
            END IF
         END IF
*
      END IF
*
*     Prepare output: set info = 0 if no error, and divide by
*     DESCMULT if error is not in a descriptor entry
*
      IF( INFO.EQ.BIGNUM ) THEN
         INFO = 0
      ELSE IF( MOD( INFO, DESCMULT ).EQ.0 ) THEN
         INFO = -( INFO / DESCMULT )
      ELSE
         INFO = -INFO
      END IF
*
      RETURN
*
*     End of PB_CHKMAT
*
      END
      SUBROUTINE PB_DESCTRANS( DESCIN, DESCOUT )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Array Arguments ..
      INTEGER            DESCIN( * ), DESCOUT( * )
*     ..
*
*  Purpose
*  =======
*
*  PB_DESCTRANS  converts  a  descriptor  DESCIN of type BLOCK_CYCLIC_2D
*  or   BLOCK_CYCLIC_INB_2D   into   a   descriptor   DESCOUT   of  type
*  BLOCK_CYCLIC_INB_2D.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information required  to
*  establish the  mapping between a matrix entry and  its  corresponding
*  process and memory location.
*
*  In  the  following  comments,  the  character _  should  be  read  as
*  "of the distributed  matrix".  Let  A  be a generic term for  any  2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN        EXPLANATION
*  ---------------- ---------------  -----------------------------------
*  DTYPE_A (global) DESCA( DTYPE1_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT1_  ) The BLACS context handle indicating
*                                    the   NPROW x NPCOL  BLACS  process
*                                    grid  A  is  distributed  over. The
*                                    context  itself  is global, but the
*                                    handle   (the  integer  value)  may
*                                    vary.
*  M_A     (global) DESCA( M1_     ) The  number  of rows in the distri-
*                                    buted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N1_     ) The  number  of columns in the dis-
*                                    tributed matrix A, N_A >= 0.
*  MB_A    (global) DESCA( MB1_    ) The blocking factor used to distri-
*                                    bute the rows of A, MB_A > 0.
*  NB_A    (global) DESCA( NB1_    ) The blocking factor used to distri-
*                                    bute the columns of A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC1_  ) The  process  row  over  which  the
*                                    first row of the matrix  A  is dis-
*                                    tributed, NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC1_  ) The process column  over  which the
*                                    first column of  A  is distributed.
*                                    NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD1_   ) The leading dimension  of the local
*                                    array  storing  the local blocks of
*                                    the distributed matrix A,
*                                    IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                    ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_NUMROC:
*  Lr( IA, K ) = PB_NUMROC( K, IA, MB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_NUMROC( K, JA, NB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information  required to
*  establish the  mapping  between a  matrix entry and its corresponding
*  process and memory location.
*
*  In  the  following  comments,   the character _  should  be  read  as
*  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN       EXPLANATION
*  ---------------- --------------- ------------------------------------
*  DTYPE_A (global) DESCA( DTYPE_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT_  ) The BLACS context handle, indicating
*                                   the NPROW x NPCOL BLACS process grid
*                                   A  is distributed over.  The context
*                                   itself  is  global,  but  the handle
*                                   (the integer value) may vary.
*  M_A     (global) DESCA( M_     ) The  number of rows in the distribu-
*                                   ted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N_     ) The number of columns in the distri-
*                                   buted matrix A, N_A >= 0.
*  IMB_A   (global) DESCA( IMB_   ) The number of rows of the upper left
*                                   block of the matrix A, IMB_A > 0.
*  INB_A   (global) DESCA( INB_   ) The  number  of columns of the upper
*                                   left   block   of   the   matrix  A,
*                                   INB_A > 0.
*  MB_A    (global) DESCA( MB_    ) The blocking factor used to  distri-
*                                   bute the last  M_A-IMB_A rows of  A,
*                                   MB_A > 0.
*  NB_A    (global) DESCA( NB_    ) The blocking factor used to  distri-
*                                   bute the last  N_A-INB_A  columns of
*                                   A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC_  ) The process row over which the first
*                                   row of the matrix  A is distributed,
*                                   NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC_  ) The  process  column  over which the
*                                   first  column of  A  is distributed.
*                                   NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD_   ) The  leading  dimension of the local
*                                   array  storing  the  local blocks of
*                                   the distributed matrix A,
*                                   IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                   ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_NUMROC:
*  Lr( IA, K ) = PB_NUMROC( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_NUMROC( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  Arguments
*  =========
*
*  DESCIN  (global and local input) INTEGER array
*          On entry, DESCIN  is an array of dimension DLEN1_ or DLEN_ as
*          specified by its first entry DESCIN( DTYPE_ ).  DESCIN is the
*          source  array  descriptor of type BLOCK_CYCLIC_2D  or of type
*          BLOCK_CYCLIC_2D_INB.
*
*  DESCOUT (global and local output) INTEGER array
*          On entry, DESCOUT is an array of dimension DLEN_.  DESCOUT is
*          the target array descriptor of type BLOCK_CYCLIC_2D_INB.
*
*  -- Written on April 1, 1998 by
*     R. Clint Whaley, University of Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC1_, CTXT1_, DLEN1_,
     $                   DTYPE1_, LLD1_, M1_, MB1_, N1_, NB1_, RSRC1_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN1_ = 9, DTYPE1_ = 1,
     $                   CTXT1_ = 2, M1_ = 3, N1_ = 4, MB1_ = 5,
     $                   NB1_ = 6, RSRC1_ = 7, CSRC1_ = 8, LLD1_ = 9 )
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. Executable Statements ..
*
      IF( DESCIN( DTYPE_ ).EQ.BLOCK_CYCLIC_2D ) THEN
         DESCOUT( DTYPE_ ) = BLOCK_CYCLIC_2D_INB
         DESCOUT( CTXT_  ) = DESCIN( CTXT1_ )
         DESCOUT( M_     ) = DESCIN( M1_    )
         DESCOUT( N_     ) = DESCIN( N1_    )
         DESCOUT( IMB_   ) = DESCIN( MB1_   )
         DESCOUT( INB_   ) = DESCIN( NB1_   )
         DESCOUT( MB_    ) = DESCIN( MB1_   )
         DESCOUT( NB_    ) = DESCIN( NB1_   )
         DESCOUT( RSRC_  ) = DESCIN( RSRC1_ )
         DESCOUT( CSRC_  ) = DESCIN( CSRC1_ )
         DESCOUT( LLD_   ) = DESCIN( LLD1_  )
      ELSE IF( DESCIN( DTYPE_ ).EQ.BLOCK_CYCLIC_2D_INB ) THEN
         DO 10 I = 1, DLEN_
            DESCOUT( I ) = DESCIN( I )
   10    CONTINUE
      ELSE
         DESCOUT( DTYPE_ ) = DESCIN( 1 )
         DESCOUT( CTXT_  ) = DESCIN( 2 )
         DESCOUT( M_     ) = 0
         DESCOUT( N_     ) = 0
         DESCOUT( IMB_   ) = 1
         DESCOUT( INB_   ) = 1
         DESCOUT( MB_    ) = 1
         DESCOUT( NB_    ) = 1
         DESCOUT( RSRC_  ) = 0
         DESCOUT( CSRC_  ) = 0
         DESCOUT( LLD_   ) = 1
      END IF
*
      RETURN
*
*     End of PB_DESCTRANS
*
      END
      SUBROUTINE PB_DESCSET2( DESC, M, N, IMB, INB, MB, NB, RSRC, CSRC,
     $                        CTXT, LLD )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            CSRC, CTXT, IMB, INB, LLD, M, MB, N, NB, RSRC
*     ..
*     .. Array Arguments ..
      INTEGER            DESC( * )
*     ..
*
*  Purpose
*  =======
*
*  PB_DESCSET2 uses  its  10  input  arguments  M,  N, IMB, INB, MB, NB,
*  RSRC,  CSRC,  CTXT  and LLD to initialize a descriptor vector of type
*  BLOCK_CYCLIC_2D_INB.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information required  to
*  establish the  mapping between a matrix entry and  its  corresponding
*  process and memory location.
*
*  In  the  following  comments,  the  character _  should  be  read  as
*  "of the distributed  matrix".  Let  A  be a generic term for  any  2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN        EXPLANATION
*  ---------------- ---------------  -----------------------------------
*  DTYPE_A (global) DESCA( DTYPE1_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT1_  ) The BLACS context handle indicating
*                                    the   NPROW x NPCOL  BLACS  process
*                                    grid  A  is  distributed  over. The
*                                    context  itself  is global, but the
*                                    handle   (the  integer  value)  may
*                                    vary.
*  M_A     (global) DESCA( M1_     ) The  number  of rows in the distri-
*                                    buted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N1_     ) The  number  of columns in the dis-
*                                    tributed matrix A, N_A >= 0.
*  MB_A    (global) DESCA( MB1_    ) The blocking factor used to distri-
*                                    bute the rows of A, MB_A > 0.
*  NB_A    (global) DESCA( NB1_    ) The blocking factor used to distri-
*                                    bute the columns of A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC1_  ) The  process  row  over  which  the
*                                    first row of the matrix  A  is dis-
*                                    tributed, NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC1_  ) The process column  over  which the
*                                    first column of  A  is distributed.
*                                    NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD1_   ) The leading dimension  of the local
*                                    array  storing  the local blocks of
*                                    the distributed matrix A,
*                                    IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                    ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_NUMROC:
*  Lr( IA, K ) = PB_NUMROC( K, IA, MB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_NUMROC( K, JA, NB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  Arguments
*  =========
*
*  DESC    (global and local output) INTEGER array
*          On entry, DESC is an array of  dimension  DLEN_.  DESC is the
*          array descriptor to be set.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies  the  number  of  rows of the matrix.
*          M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry,  N  specifies  the number of columns of the matrix.
*          N must be at least zero.
*
*  IMB     (global input) INTEGER
*          On entry,  IMB  specifies  the row size of the first block of
*          the global matrix distribution. IMB must be at least one.
*
*  INB     (global input) INTEGER
*          On entry,  INB  specifies  the column size of the first block
*          of the global matrix distribution. INB must be at least one.
*
*  MB      (global input) INTEGER
*          On entry,  MB  specifies  the  row size of the blocks used to
*          partition the matrix. MB must be at least one.
*
*  NB      (global input) INTEGER
*          On entry, NB  specifies the column size of the blocks used to
*          partition the matrix. NB must be at least one.
*
*  RSRC    (global input) INTEGER
*          On entry,  RSRC  specifies  the row coordinate of the process
*          that possesses the first row of the matrix.  When  RSRC = -1,
*          the data is not  distributed but replicated,  otherwise  RSRC
*          must be at least zero and strictly less than NPROW.
*
*  CSRC    (global input) INTEGER
*          On entry,  CSRC  specifies  the column coordinate of the pro-
*          cess  that  possesses  the  first column of the matrix.  When
*          CSRC = -1, the data is not distributed but replicated, other-
*          wise CSRC must be at least zero and strictly less than NPCOL.
*
*  CTXT    (local input) INTEGER
*          On entry, CTXT specifies the BLACS context handle, indicating
*          the global  communication  context.  The value of the context
*          itself is local.
*
*  LLD     (local input)  INTEGER
*          On entry, LLD  specifies  the  leading dimension of the local
*          array storing the local entries of the matrix. LLD must be at
*          least MAX( 1, Lr(1,M) ).
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Executable Statements ..
*
      DESC( DTYPE_ ) = BLOCK_CYCLIC_2D_INB
      DESC( CTXT_  ) = CTXT
      DESC( M_     ) = M
      DESC( N_     ) = N
      DESC( IMB_   ) = IMB
      DESC( INB_   ) = INB
      DESC( MB_    ) = MB
      DESC( NB_    ) = NB
      DESC( RSRC_  ) = RSRC
      DESC( CSRC_  ) = CSRC
      DESC( LLD_   ) = LLD
*
      RETURN
*
*     End of PB_DESCSET2
*
      END
      SUBROUTINE PB_DESCINIT2( DESC, M, N, IMB, INB, MB, NB, RSRC, CSRC,
     $                         CTXT, LLD, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            CSRC, CTXT, IMB, INB, INFO, LLD, M, MB, N, NB,
     $                   RSRC
*     ..
*     .. Array Arguments ..
      INTEGER            DESC( * )
*     ..
*
*  Purpose
*  =======
*
*  PB_DESCINIT2 uses  its  10  input  arguments  M, N, IMB, INB, MB, NB,
*  RSRC,  CSRC,  CTXT  and LLD to initialize a descriptor vector of type
*  BLOCK_CYCLIC_2D_INB.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information  required to
*  establish the  mapping  between a  matrix entry and its corresponding
*  process and memory location.
*
*  In  the  following  comments,   the character _  should  be  read  as
*  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN       EXPLANATION
*  ---------------- --------------- ------------------------------------
*  DTYPE_A (global) DESCA( DTYPE_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT_  ) The BLACS context handle, indicating
*                                   the NPROW x NPCOL BLACS process grid
*                                   A  is distributed over.  The context
*                                   itself  is  global,  but  the handle
*                                   (the integer value) may vary.
*  M_A     (global) DESCA( M_     ) The  number of rows in the distribu-
*                                   ted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N_     ) The number of columns in the distri-
*                                   buted matrix A, N_A >= 0.
*  IMB_A   (global) DESCA( IMB_   ) The number of rows of the upper left
*                                   block of the matrix A, IMB_A > 0.
*  INB_A   (global) DESCA( INB_   ) The  number  of columns of the upper
*                                   left   block   of   the   matrix  A,
*                                   INB_A > 0.
*  MB_A    (global) DESCA( MB_    ) The blocking factor used to  distri-
*                                   bute the last  M_A-IMB_A rows of  A,
*                                   MB_A > 0.
*  NB_A    (global) DESCA( NB_    ) The blocking factor used to  distri-
*                                   bute the last  N_A-INB_A  columns of
*                                   A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC_  ) The process row over which the first
*                                   row of the matrix  A is distributed,
*                                   NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC_  ) The  process  column  over which the
*                                   first  column of  A  is distributed.
*                                   NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD_   ) The  leading  dimension of the local
*                                   array  storing  the  local blocks of
*                                   the distributed matrix A,
*                                   IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                   ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_NUMROC:
*  Lr( IA, K ) = PB_NUMROC( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_NUMROC( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  Arguments
*  =========
*
*  DESC    (global and local output) INTEGER array
*          On entry, DESC is an array of  dimension  DLEN_.  DESC is the
*          array descriptor to be set.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies  the  number  of  rows of the matrix.
*          M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry,  N  specifies  the number of columns of the matrix.
*          N must be at least zero.
*
*  IMB     (global input) INTEGER
*          On entry,  IMB  specifies  the row size of the first block of
*          the global matrix distribution. IMB must be at least one.
*
*  INB     (global input) INTEGER
*          On entry,  INB  specifies  the column size of the first block
*          of the global matrix distribution. INB must be at least one.
*
*  MB      (global input) INTEGER
*          On entry,  MB  specifies  the  row size of the blocks used to
*          partition the matrix. MB must be at least one.
*
*  NB      (global input) INTEGER
*          On entry, NB  specifies the column size of the blocks used to
*          partition the matrix. NB must be at least one.
*
*  RSRC    (global input) INTEGER
*          On entry,  RSRC  specifies  the row coordinate of the process
*          that possesses the first row of the matrix.  When  RSRC = -1,
*          the data is not  distributed but replicated,  otherwise  RSRC
*          must be at least zero and strictly less than NPROW.
*
*  CSRC    (global input) INTEGER
*          On entry,  CSRC  specifies  the column coordinate of the pro-
*          cess  that  possesses  the  first column of the matrix.  When
*          CSRC = -1, the data is not distributed but replicated, other-
*          wise CSRC must be at least zero and strictly less than NPCOL.
*
*  CTXT    (local input) INTEGER
*          On entry, CTXT specifies the BLACS context handle, indicating
*          the global  communication  context.  The value of the context
*          itself is local.
*
*  LLD     (local input)  INTEGER
*          On entry, LLD  specifies  the  leading dimension of the local
*          array storing the local entries of the matrix. LLD must be at
*          least MAX( 1, Lr(1,M) ).
*
*  INFO    (local output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value.
*
*  Notes
*  =====
*
*  If the routine can recover from an erroneous input argument,  it will
*  return an acceptable descriptor vector.  For example,  if LLD = 0  on
*  input, DESC( LLD_ ) will  contain  the smallest leading dimension re-
*  quired to store the specified m by n matrix, INFO will however be set
*  to -11 on exit in that case.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D_INB, CSRC_, CTXT_, DLEN_,
     $                   DTYPE_, IMB_, INB_, LLD_, MB_, M_, NB_, N_,
     $                   RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D_INB = 2, DLEN_ = 11,
     $                   DTYPE_ = 1, CTXT_ = 2, M_ = 3, N_ = 4,
     $                   IMB_ = 5, INB_ = 6, MB_ = 7, NB_ = 8,
     $                   RSRC_ = 9, CSRC_ = 10, LLD_ = 11 )
*     ..
*     .. Local Scalars ..
      INTEGER            LLDMIN, MP, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PXERBLA
*     ..
*     .. External Functions ..
      INTEGER            PB_NUMROC
      EXTERNAL           PB_NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      CALL BLACS_GRIDINFO( CTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( IMB.LT.1 ) THEN
         INFO = -4
      ELSE IF( INB.LT.1 ) THEN
         INFO = -5
      ELSE IF( MB.LT.1 ) THEN
         INFO = -6
      ELSE IF( NB.LT.1 ) THEN
         INFO = -7
      ELSE IF( RSRC.LT.-1 .OR. RSRC.GE.NPROW ) THEN
         INFO = -8
      ELSE IF( CSRC.LT.-1 .OR. CSRC.GE.NPCOL ) THEN
         INFO = -9
      ELSE IF( NPROW.EQ.-1 ) THEN
         INFO = -10
      END IF
*
*     Compute minimum LLD if safe (to avoid division by 0)
*
      IF( INFO.EQ.0 ) THEN
         MP = PB_NUMROC( M, 1, IMB, MB, MYROW, RSRC, NPROW )
         IF( PB_NUMROC( N, 1, INB, NB, MYCOL, CSRC, NPCOL ).GT.0 ) THEN
            LLDMIN = MAX( 1, MP )
         ELSE
            LLDMIN = 1
         END IF
         IF( LLD.LT.LLDMIN )
     $      INFO = -11
      END IF
*
      IF( INFO.NE.0 )
     $   CALL PXERBLA( CTXT, 'PB_DESCINIT2', -INFO )
*
      DESC( DTYPE_ ) = BLOCK_CYCLIC_2D_INB
      DESC( CTXT_  ) = CTXT
      DESC( M_     ) = MAX( 0, M )
      DESC( N_     ) = MAX( 0, N )
      DESC( IMB_   ) = MAX( 1, IMB )
      DESC( INB_   ) = MAX( 1, INB )
      DESC( MB_    ) = MAX( 1, MB )
      DESC( NB_    ) = MAX( 1, NB )
      DESC( RSRC_  ) = MAX( -1, MIN( RSRC, NPROW-1 ) )
      DESC( CSRC_  ) = MAX( -1, MIN( CSRC, NPCOL-1 ) )
      DESC( LLD_   ) = MAX( LLD, LLDMIN )
*
      RETURN
*
*     End of PB_DESCINIT2
*
      END
      SUBROUTINE PB_BINFO( OFFD, M, N, IMB1, INB1, MB, NB, MRROW, MRCOL,
     $                     LCMT00, MBLKS, NBLKS, IMBLOC, INBLOC, LMBLOC,
     $                     LNBLOC, ILOW, LOW, IUPP, UPP )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ILOW, IMB1, IMBLOC, INB1, INBLOC, IUPP, LCMT00,
     $                   LMBLOC, LNBLOC, LOW, M, MB, MBLKS, MRCOL,
     $                   MRROW, N, NB, NBLKS, OFFD, UPP
*     ..
*
*  Purpose
*  =======
*
*  PB_BINFO   initializes the local information of an m by n local array
*  owned by the process of  relative  coordinates ( MRROW, MRCOL ). Note
*  that if m or n is less or equal than zero, there is no data, in which
*  case this process  does  not  need  the local information computed by
*  this routine to proceed.
*
*  Arguments
*  =========
*
*  OFFD    (global input) INTEGER
*          On entry,  OFFD  specifies the off-diagonal of the underlying
*          matrix of interest as follows:
*             OFFD = 0 specifies the main diagonal,
*             OFFD > 0 specifies lower subdiagonals, and
*             OFFD < 0 specifies upper superdiagonals.
*
*  M       (local input) INTEGER
*          On entry, M  specifies the local number of rows of the under-
*          lying matrix  owned  by the  process  of relative coordinates
*          ( MRROW, MRCOL ). M must be at least zero.
*
*  N       (local input) INTEGER
*          On entry, N  specifies the local number of columns of the un-
*          derlying matrix  owned by the process of relative coordinates
*          ( MRROW, MRCOL ). N must be at least zero.
*
*  IMB1    (global input) INTEGER
*          On input, IMB1 specifies  the global true size of  the  first
*          block of rows of the underlying global submatrix.  IMB1  must
*          be at least MIN( 1, M ).
*
*  INB1    (global input) INTEGER
*          On input, INB1 specifies  the global true size of  the  first
*          block  of  columns  of  the underlying global submatrix. INB1
*          must be at least MIN( 1, N ).
*
*  MB      (global input) INTEGER
*          On entry, MB  specifies the blocking factor used to partition
*          the rows of the matrix.  MB  must be at least one.
*
*  NB      (global input) INTEGER
*          On entry, NB  specifies the blocking factor used to partition
*          the the columns of the matrix.  NB  must be at least one.
*
*  MRROW   (local input) INTEGER
*          On entry, MRROW specifies the  relative row coordinate of the
*          process that possesses these M rows. MRROW must be least zero
*          and strictly less than NPROW.
*
*  MRCOL   (local input) INTEGER
*          On entry, MRCOL specifies  the  relative column coordinate of
*          the process that possesses these N  columns.  MRCOL  must  be
*          least zero and strictly less than NPCOL.
*
*  LCMT00  (local output) INTEGER
*          On exit, LCMT00  is the  LCM value of the left upper block of
*          this m by n local  block owned by the process of relative co-
*          ordinates ( MRROW, MRCOL ).
*
*  MBLKS   (local output) INTEGER
*          On exit, MBLKS specifies the local number of blocks  of  rows
*          corresponding to M. MBLKS must be at least zero.
*
*  NBLKS   (local output) INTEGER
*          On exit,  NBLKS  specifies  the local number of blocks of co-
*          lumns corresponding to N. NBLKS must be at least zero.
*
*  IMBLOC  (local output) INTEGER
*          On exit, IMBLOC  specifies  the  number of rows (size) of the
*          uppest blocks of this m by n local array owned by the process
*          of relative coordinates ( MRROW, MRCOL ).  IMBLOC is at least
*          MIN( 1, M ).
*
*  INBLOC  (local output) INTEGER
*          On exit, INBLOC  specifies  the  number of columns (size) of
*          the leftmost  blocks of this m by n local array owned by the
*          process of relative coordinates ( MRROW, MRCOL ).  INBLOC is
*          at least MIN( 1, N ).
*
*  LMBLOC  (local output) INTEGER
*          On exit, LMBLOC specifies the number  of  rows  (size) of the
*          lowest blocks of this m by n local array owned by the process
*          of  relative coordinates ( MRROW, MRCOL ). LMBLOC is at least
*          MIN( 1, M ).
*
*  LNBLOC  (local output) INTEGER
*          On exit, LNBLOC specifies the number of columns (size) of the
*          rightmost  blocks of this  m by n  local  array  owned by the
*          process of  relative  coordinates ( MRROW, MRCOL ). LNBLOC is
*          at least MIN( 1, N ).
*
*  ILOW    (local output) INTEGER
*          On exit, ILOW is the lower bound characterizing the first co-
*          lumn block owning offdiagonals of  this  m by n  array.  ILOW
*          must be less or equal than zero.
*
*  LOW     (global output) INTEGER
*          On exit,  LOW  is  the  lower bound characterizing the column
*          blocks with te exception of the  first  one (see ILOW) owning
*          offdiagonals of this m by n array. LOW  must be less or equal
*          than zero.
*
*  IUPP    (local output) INTEGER
*          On exit, IUPP is the upper bound characterizing the first row
*          block owning offdiagonals of this m by n array.  IUPP must be
*          greater or equal than zero.
*
*  UPP     (global output) INTEGER
*          On exit,  UPP  is  the  upper  bound  characterizing  the row
*          blocks with te exception of the  first  one (see IUPP) owning
*          offdiagonals of this m by n array. UPP  must  be  greater  or
*          equal than zero.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            TMP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Initialize LOW, ILOW, UPP, IUPP, LMBLOC, LNBLOC, IMBLOC, INBLOC,
*     MBLKS, NBLKS and LCMT00.
*
      LOW = 1 - NB
      UPP = MB - 1
*
      LCMT00 = OFFD
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
*
         IF( MRROW.GT.0 ) THEN
            IUPP = MB - 1
         ELSE
            IUPP = MAX( 0, IMB1 - 1 )
         END IF
         IMBLOC = 0
         MBLKS  = 0
         LMBLOC = 0
*
         IF( MRCOL.GT.0 ) THEN
            ILOW = 1 - NB
         ELSE
            ILOW = MIN( 0, 1 - INB1 )
         END IF
         INBLOC = 0
         NBLKS  = 0
         LNBLOC = 0
*
         LCMT00 = LCMT00 + ( LOW - ILOW + MRCOL * NB ) -
     $            ( IUPP - UPP + MRROW * MB )
*
         RETURN
*
      END IF
*
      IF( MRROW.GT.0 ) THEN
*
         IMBLOC = MIN( M, MB )
         IUPP   = MB - 1
         LCMT00 = LCMT00 - ( IMB1 - MB + MRROW * MB )
         MBLKS  = ( M - 1 ) / MB + 1
         LMBLOC = M - ( M / MB ) * MB
         IF( LMBLOC.EQ.0 )
     $      LMBLOC = MB
*
         IF( MRCOL.GT.0 ) THEN
*
            INBLOC = MIN( N, NB )
            ILOW   = 1 - NB
            LCMT00 = LCMT00 + INB1 - NB + MRCOL * NB
            NBLKS  = ( N - 1 ) / NB + 1
            LNBLOC = N - ( N / NB ) * NB
            IF( LNBLOC.EQ.0 )
     $         LNBLOC = NB
*
         ELSE
*
            INBLOC = INB1
            ILOW   = 1 - INB1
            TMP1   = N - INB1
            IF( TMP1.GT.0 ) THEN
*
*              more than one block
*
               NBLKS = ( TMP1 - 1 ) / NB + 2
               LNBLOC = TMP1 - ( TMP1 / NB ) * NB
               IF( LNBLOC.EQ.0 )
     $            LNBLOC = NB
*
            ELSE
*
               NBLKS  = 1
               LNBLOC = INB1
*
            END IF
*
         END IF
*
      ELSE
*
         IMBLOC = IMB1
         IUPP = IMB1 - 1
         TMP1 = M - IMB1
         IF( TMP1.GT.0 ) THEN
*
*           more than one block
*
            MBLKS  = ( TMP1 - 1 ) / MB + 2
            LMBLOC = TMP1 - ( TMP1 / MB ) * MB
            IF( LMBLOC.EQ.0 )
     $         LMBLOC = MB
*
         ELSE
*
            MBLKS  = 1
            LMBLOC = IMB1
*
         END IF
*
         IF( MRCOL.GT.0 ) THEN
*
            INBLOC = MIN( N, NB )
            ILOW   = 1 - NB
            LCMT00 = LCMT00 + INB1 - NB + MRCOL * NB
            NBLKS  = ( N - 1 ) / NB + 1
            LNBLOC = N - ( N / NB ) * NB
            IF( LNBLOC.EQ.0 )
     $         LNBLOC = NB
*
         ELSE
*
            INBLOC = INB1
            ILOW   = 1 - INB1
            TMP1   = N - INB1
            IF( TMP1.GT.0 ) THEN
*
*              more than one block
*
               NBLKS  = ( TMP1 - 1 ) / NB + 2
               LNBLOC = TMP1 - ( TMP1 / NB ) * NB
               IF( LNBLOC.EQ.0 )
     $            LNBLOC = NB
*
            ELSE
*
               NBLKS  = 1
               LNBLOC = INB1
*
            END IF
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PB_BINFO
*
      END
      INTEGER FUNCTION PILAENV( ICTXT, PREC )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT
      CHARACTER*1        PREC
*     ..
*
*  Purpose
*  =======
*
*  PILAENV  returns  the  logical computational block size to be used by
*  the PBLAS routines during testing and timing. This is a special  ver-
*  sion to be used only as part of the testing or timing  PBLAS programs
*  for testing different values of logical computational block sizes for
*  the PBLAS routines. It is called by the PBLAS routines to  retrieve a
*  logical computational block size value.
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  PREC    (dummy input) CHARACTER*1
*          On entry, PREC is a dummy argument.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Common Blocks ..
      INTEGER            INFO, NBLOG
      COMMON             /INFOC/INFO, NBLOG
*     ..
*     .. Executable Statements ..
*
      PILAENV = NBLOG
*
      RETURN
*
*     End of PILAENV
*
      END
      SUBROUTINE PB_LOCINFO( I, INB, NB, MYROC, SRCPROC, NPROCS,
     $                       ILOCBLK, ILOCOFF, MYDIST )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            I, ILOCBLK, ILOCOFF, INB, MYDIST, MYROC, NB,
     $                   NPROCS, SRCPROC
*     ..
*
*  Purpose
*  =======
*
*  PB_LOCINFO  computes  local information about the beginning of a sub-
*  matrix starting at the global index I.
*
*  Arguments
*  =========
*
*  I       (global input) INTEGER
*          On entry,  I  specifies  the global starting index in the ma-
*          trix. I must be at least one.
*
*  INB     (global input) INTEGER
*          On entry,  INB  specifies the size of the first block of rows
*          or columns of the matrix. INB must be at least one.
*
*  NB      (global input) INTEGER
*          On entry, NB  specifies the size of the blocks of rows or co-
*          lumns of the matrix is partitioned into.  NB must be at least
*          one.
*
*  MYROC   (local input) INTEGER
*          On entry, MYROC is the  coordinate of the process whose local
*          information  is  determined.  MYROC  is  at  least  zero  and
*          strictly less than NPROCS.
*
*  SRCPROC (global input) INTEGER
*          On entry,  SRCPROC  specifies  the coordinate of the  process
*          that possesses the  first row or column  of the matrix.  When
*          SRCPROC = -1, the data  is not  distributed  but  replicated,
*          otherwise  SRCPROC  must be at least zero and  strictly  less
*          than NPROCS.
*
*  NPROCS  (global input) INTEGER
*          On entry, NPROCS  specifies  the total number of process rows
*          or  columns  over  which the submatrix is distributed. NPROCS
*          must be at least one.
*
*  ILOCBLK (local output) INTEGER
*          On exit, ILOCBLK  specifies  the  local  row  or column block
*          coordinate  corresponding  to  the row or column I of the ma-
*          trix. ILOCBLK must be at least zero.
*
*  ILOCOFF (local output) INTEGER
*          On exit, ILOCOFF  specifies the local row offset in the block
*          of local coordinate  ILOCBLK  corresponding to the row or co-
*          lumn I of the matrix. ILOCOFF must at least zero.
*
*  MYDIST  (local output) INTEGER
*          On exit, MYDIST  specifies the relative process coordinate of
*          the process specified by MYROC to the process owning the  row
*          or column I. MYDIST  is at  least zero and strictly less than
*          NPROCS.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            ITMP, NBLOCKS, PROC
*     ..
*     .. Executable Statements ..
*
      ILOCOFF = 0
*
      IF( SRCPROC.LT.0 ) THEN
*
         MYDIST = 0
*
         IF( I.LE.INB ) THEN
*
            ILOCBLK = 0
            ILOCOFF = I - 1
*
         ELSE
*
            ITMP    = I - INB
            NBLOCKS = ( ITMP - 1 ) / NB + 1
            ILOCBLK = NBLOCKS
            ILOCOFF = ITMP - 1 - ( NBLOCKS - 1 ) * NB
*
         END IF
*
      ELSE
*
         PROC   = SRCPROC
         MYDIST = MYROC - PROC
         IF( MYDIST.LT.0 )
     $      MYDIST = MYDIST + NPROCS
*
         IF( I.LE.INB ) THEN
*
            ILOCBLK = 0
            IF( MYROC.EQ.PROC )
     $         ILOCOFF = I - 1
*
         ELSE
*
            ITMP    = I - INB
            NBLOCKS = ( ITMP - 1 ) / NB + 1
            PROC    = PROC + NBLOCKS
            PROC    = PROC - ( PROC / NPROCS ) * NPROCS
            ILOCBLK = NBLOCKS / NPROCS
*
            IF( ( ILOCBLK*NPROCS ).LT.( MYDIST-NBLOCKS ) )
     $         ILOCBLK = ILOCBLK + 1
*
            IF( MYROC.EQ.PROC )
     $         ILOCOFF = ITMP - 1 - ( NBLOCKS - 1 ) * NB
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PB_LOCINFO
*
      END
      SUBROUTINE PB_INITJMP( COLMAJ, NVIR, IMBVIR, INBVIR, IMBLOC,
     $                       INBLOC, MB, NB, RSRC, CSRC, NPROW, NPCOL,
     $                       STRIDE, JMP )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      LOGICAL            COLMAJ
      INTEGER            CSRC, IMBLOC, IMBVIR, INBLOC, INBVIR, MB, NB,
     $                   NPCOL, NPROW, NVIR, RSRC, STRIDE
*     ..
*     .. Array Arguments ..
      INTEGER            JMP( * )
*     ..
*
*  Purpose
*  =======
*
*  PB_INITJMP  initializes the jump values JMP used by the random matrix
*  generator.
*
*  Arguments
*  =========
*
*  COLMAJ  (global input) LOGICAL
*          On entry, COLMAJ specifies the ordering of the random sequen-
*          ce. When  COLMAJ is .TRUE.,  the random sequence will be used
*          for a column major ordering, and otherwise a  row-major orde-
*          ring. This impacts on the computation of the jump values.
*
*  NVIR    (global input) INTEGER
*          On entry, NVIR  specifies  the size of the underlying virtual
*          matrix. NVIR must be at least zero.
*
*  IMBVIR  (local input) INTEGER
*          On entry, IMBVIR  specifies the number of virtual rows of the
*          upper left block of the underlying virtual submatrix.  IMBVIR
*          must be at least IMBLOC.
*
*  INBVIR  (local input) INTEGER
*          On entry, INBVIR  specifies  the number of virtual columns of
*          the  upper  left  block  of the underlying virtual submatrix.
*          INBVIR must be at least INBLOC.
*
*  IMBLOC  (local input) INTEGER
*          On entry, IMBLOC specifies  the  number of rows (size) of the
*          local uppest  blocks. IMBLOC is at least zero.
*
*  INBLOC  (local input) INTEGER
*          On entry,  INBLOC  specifies the number of columns (size)  of
*          the local leftmost blocks. INBLOC is at least zero.
*
*  MB      (global input) INTEGER
*          On entry, MB specifies the size of the blocks used to  parti-
*          tion the matrix rows. MB must be at least one.
*
*  NB      (global input) INTEGER
*          On entry, NB specifies the size of the blocks used to  parti-
*          tion the matrix columns. NB must be at least one.
*
*  RSRC    (global input) INTEGER
*          On entry,  RSRC  specifies the row coordinate of the  process
*          that possesses the  first row of the matrix.  When RSRC = -1,
*          the rows are not distributed but replicated,  otherwise  RSRC
*          must be at least zero and  strictly less than NPROW.
*
*  CSRC    (global input) INTEGER
*          On entry,  CSRC  specifies  the column coordinate of the pro-
*          cess that possesses the first column of the matrix. When CSRC
*          is equal to -1,  the columns are not distributed but replica-
*          ted, otherwise  CSRC  must be at least zero and strictly less
*          than NPCOL.
*
*  NPROW   (global input) INTEGER
*          On entry,  NPROW  specifies  the total number of process rows
*          over which the matrix is distributed.  NPROW must be at least
*          one.
*
*  NPCOL   (global input) INTEGER
*          On entry,  NPCOL  specifies  the  total number of process co-
*          lumns over which the matrix is distributed.  NPCOL must be at
*          least one.
*
*  STRIDE  (global input) INTEGER
*          On entry, STRIDE specifies the number of random numbers to be
*          generated  to  compute  one  matrix  entry. In the real case,
*          STRIDE is usually 1,  where  as in the complex case STRIDE is
*          usually 2 in order to generate the real and imaginary parts.
*
*  JMP     (local output) INTEGER array
*          On entry, JMP is an array of dimension JMP_LEN. On exit, this
*          array contains  the different  jump values used by the random
*          matrix generator.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            JMP_1, JMP_COL, JMP_IMBV, JMP_INBV, JMP_LEN,
     $                   JMP_MB, JMP_NB, JMP_NPIMBLOC, JMP_NPMB,
     $                   JMP_NQINBLOC, JMP_NQNB, JMP_ROW
      PARAMETER          ( JMP_1 = 1, JMP_ROW = 2, JMP_COL = 3,
     $                   JMP_MB = 4, JMP_IMBV = 5, JMP_NPMB = 6,
     $                   JMP_NPIMBLOC = 7, JMP_NB = 8, JMP_INBV = 9,
     $                   JMP_NQNB = 10, JMP_NQINBLOC = 11,
     $                   JMP_LEN = 11 )
*     ..
*     .. Local Scalars ..
      INTEGER            NPMB, NQNB
*     ..
*     .. Executable Statements ..
*
      IF( RSRC.LT.0 ) THEN
         NPMB = MB
      ELSE
         NPMB = NPROW * MB
      END IF
      IF( CSRC.LT.0 ) THEN
         NQNB = NB
      ELSE
         NQNB = NPCOL * NB
      END IF
*
      JMP( JMP_1        ) = 1
*
      JMP( JMP_MB       ) = MB
      JMP( JMP_IMBV     ) = IMBVIR
      JMP( JMP_NPMB     ) = NPMB
      JMP( JMP_NPIMBLOC ) = IMBLOC + NPMB - MB
*
      JMP( JMP_NB       ) = NB
      JMP( JMP_INBV     ) = INBVIR
      JMP( JMP_NQNB     ) = NQNB
      JMP( JMP_NQINBLOC ) = INBLOC + NQNB - NB
*
      IF( COLMAJ ) THEN
         JMP( JMP_ROW ) = STRIDE
         JMP( JMP_COL ) = STRIDE * NVIR
      ELSE
         JMP( JMP_ROW ) = STRIDE * NVIR
         JMP( JMP_COL ) = STRIDE
      END IF
*
      RETURN
*
*     End of PB_INITJMP
*
      END
      SUBROUTINE PB_INITMULADD( MULADD0, JMP, IMULADD )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Array Arguments ..
      INTEGER            IMULADD( 4, * ), JMP( * ), MULADD0( * )
*     ..
*
*  Purpose
*  =======
*
*  PB_INITMULADD initializes the  constants a's and c's corresponding to
*  the jump values (JMP) used by the matrix generator.
*
*  Arguments
*  =========
*
*  MULADD0 (local input) INTEGER array
*          On entry,  MULADD0  is an array of dimension 4 containing the
*          encoded  initial  constants  a and c to jump from  X( n )  to
*          X( n+1 ) = a*X( n ) + c in the random sequence.  MULADD0(1:2)
*          contains respectively the 16-lower and  16-higher bits of the
*          constant  a,  and  MULADD0(3:4)  contains  the  16-lower  and
*          16-higher bits of the constant c.
*
*  JMP     (local input) INTEGER array
*          On entry, JMP is an array of dimension JMP_LEN containing the
*          different jump values used by the matrix generator.
*
*  IMULADD (local output) INTEGER array
*          On entry, IMULADD is an array of dimension ( 4, JMP_LEN ). On
*          exit,  the jth column of this array contains the encoded ini-
*          tial constants a_j and c_j to jump from X( n ) to X(n+JMP(j))
*          (= a_j*X( n ) + c_j) in the random  sequence.  IMULADD(1:2,j)
*          contains  respectively the 16-lower and 16-higher bits of the
*          constant  a_j,  and  IMULADD(3:4,j) contains the 16-lower and
*          16-higher bits of the constant c_j.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            JMP_1, JMP_COL, JMP_IMBV, JMP_INBV, JMP_LEN,
     $                   JMP_MB, JMP_NB, JMP_NPIMBLOC, JMP_NPMB,
     $                   JMP_NQINBLOC, JMP_NQNB, JMP_ROW
      PARAMETER          ( JMP_1 = 1, JMP_ROW = 2, JMP_COL = 3,
     $                   JMP_MB = 4, JMP_IMBV = 5, JMP_NPMB = 6,
     $                   JMP_NPIMBLOC = 7, JMP_NB = 8, JMP_INBV = 9,
     $                   JMP_NQNB = 10, JMP_NQINBLOC = 11,
     $                   JMP_LEN = 11 )
*     ..
*
*     .. Local Arrays ..
      INTEGER            ITMP1( 2 ), ITMP2( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           PB_JUMP
*     ..
*     .. Executable Statements ..
*
      ITMP2( 1 ) = 100
      ITMP2( 2 ) = 0
*
*     Compute IMULADD for all JMP values
*
      CALL PB_JUMP( JMP( JMP_1   ), MULADD0, ITMP2, ITMP1,
     $              IMULADD( 1, JMP_1   ) )
*
      CALL PB_JUMP( JMP( JMP_ROW ), MULADD0, ITMP1, ITMP2,
     $              IMULADD( 1, JMP_ROW ) )
      CALL PB_JUMP( JMP( JMP_COL ), MULADD0, ITMP1, ITMP2,
     $              IMULADD( 1, JMP_COL ) )
*
*     Compute constants a and c to jump JMP( * ) numbers in the
*     sequence for column- or row-major ordering of the sequence.
*
      CALL PB_JUMP( JMP( JMP_IMBV     ), IMULADD( 1, JMP_ROW ), ITMP1,
     $              ITMP2, IMULADD( 1, JMP_IMBV     ) )
      CALL PB_JUMP( JMP( JMP_MB       ), IMULADD( 1, JMP_ROW ), ITMP1,
     $              ITMP2, IMULADD( 1, JMP_MB       ) )
      CALL PB_JUMP( JMP( JMP_NPMB     ), IMULADD( 1, JMP_ROW ), ITMP1,
     $              ITMP2, IMULADD( 1, JMP_NPMB     ) )
      CALL PB_JUMP( JMP( JMP_NPIMBLOC ), IMULADD( 1, JMP_ROW ), ITMP1,
     $              ITMP2, IMULADD( 1, JMP_NPIMBLOC ) )
*
      CALL PB_JUMP( JMP( JMP_INBV     ), IMULADD( 1, JMP_COL ), ITMP1,
     $              ITMP2, IMULADD( 1, JMP_INBV     ) )
      CALL PB_JUMP( JMP( JMP_NB       ), IMULADD( 1, JMP_COL ), ITMP1,
     $              ITMP2, IMULADD( 1, JMP_NB       ) )
      CALL PB_JUMP( JMP( JMP_NQNB     ), IMULADD( 1, JMP_COL ), ITMP1,
     $              ITMP2, IMULADD( 1, JMP_NQNB     ) )
      CALL PB_JUMP( JMP( JMP_NQINBLOC ), IMULADD( 1, JMP_COL ), ITMP1,
     $              ITMP2, IMULADD( 1, JMP_NQINBLOC ) )
*
      RETURN
*
*     End of PB_INITMULADD
*
      END
      SUBROUTINE PB_SETLOCRAN( SEED, ILOCBLK, JLOCBLK, ILOCOFF, JLOCOFF,
     $                         MYRDIST, MYCDIST, NPROW, NPCOL, JMP,
     $                         IMULADD, IRAN )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ILOCBLK, ILOCOFF, JLOCBLK, JLOCOFF, MYCDIST,
     $                   MYRDIST, NPCOL, NPROW, SEED
*     ..
*     .. Array Arguments ..
      INTEGER            IMULADD( 4, * ), IRAN( * ), JMP( * )
*     ..
*
*  Purpose
*  =======
*
*  PB_SETLOCRAN locally initializes the random number generator.
*
*  Arguments
*  =========
*
*  SEED    (global input) INTEGER
*          On entry, SEED specifies a positive integer used to initiali-
*          ze the first number in the random sequence used by the matrix
*          generator. SEED must be at least zero.
*
*  ILOCBLK (local input) INTEGER
*          On entry,  ILOCBLK  specifies  the local row block coordinate
*          corresponding to the first row of the submatrix of  interest.
*          ILOCBLK must be at least zero.
*
*  ILOCOFF (local input) INTEGER
*          On entry, ILOCOFF specifies the local row offset in the block
*          of local coordinate ILOCBLK corresponding to the first row of
*          the submatrix of interest. ILOCOFF must at least zero.
*
*  JLOCBLK (local input) INTEGER
*          On entry, JLOCBLK specifies the local column block coordinate
*          corresponding to the first column of  the  submatrix of inte-
*          rest. JLOCBLK must be at least zero.
*
*  JLOCOFF (local input) INTEGER
*          On entry,  JLOCOFF  specifies  the local column offset in the
*          block of local coordinate  JLOCBLK corresponding to the first
*          column of the submatrix of interest. JLOCOFF must be at least
*          zero.
*
*  MYRDIST (local input) INTEGER
*          On entry, MYRDIST  specifies the relative row process coordi-
*          nate to the process  owning the first row of the submatrix of
*          interest. MYRDIST must be at least zero and stricly less than
*          NPROW (see the subroutine PB_LOCINFO).
*
*  MYCDIST (local input) INTEGER
*          On entry, MYCDIST specifies the relative column process coor-
*          dinate to the  process  owning the first column of the subma-
*          trix of interest.  MYCDIST  must be at least zero and stricly
*          less than NPCOL (see the subroutine PB_LOCINFO).
*
*  NPROW   (global input) INTEGER
*          On entry,  NPROW  specifies  the total number of process rows
*          over which the matrix is distributed.  NPROW must be at least
*          one.
*
*  NPCOL   (global input) INTEGER
*          On entry,  NPCOL  specifies  the  total number of process co-
*          lumns over which the matrix is distributed.  NPCOL must be at
*          least one.
*
*  JMP     (local input) INTEGER array
*          On entry, JMP is an array of dimension JMP_LEN containing the
*          different jump values used by the matrix generator.
*
*  IMULADD (local input) INTEGER array
*          On entry, IMULADD is an array of dimension (4, JMP_LEN).  The
*          jth  column  of this array contains the encoded initial cons-
*          tants a_j and c_j to jump  from  X( n ) to  X( n + JMP( j ) )
*          (= a_j * X( n ) + c_j) in the random sequence. IMULADD(1:2,j)
*          contains respectively the 16-lower and 16-higher bits of  the
*          constant a_j, and IMULADD(3:4,j)  contains  the 16-lower  and
*          16-higher bits of the constant c_j.
*
*  IRAN    (local output) INTEGER array
*          On entry, IRAN is an array of dimension 2. On exit, IRAN con-
*          tains respectively the 16-lower and 32-higher bits of the en-
*          coding of the entry of the  random sequence corresponding lo-
*          cally to the first local array entry to generate.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            JMP_1, JMP_COL, JMP_IMBV, JMP_INBV, JMP_LEN,
     $                   JMP_MB, JMP_NB, JMP_NPIMBLOC, JMP_NPMB,
     $                   JMP_NQINBLOC, JMP_NQNB, JMP_ROW
      PARAMETER          ( JMP_1 = 1, JMP_ROW = 2, JMP_COL = 3,
     $                   JMP_MB = 4, JMP_IMBV = 5, JMP_NPMB = 6,
     $                   JMP_NPIMBLOC = 7, JMP_NB = 8, JMP_INBV = 9,
     $                   JMP_NQNB = 10, JMP_NQINBLOC = 11,
     $                   JMP_LEN = 11 )
*     ..
*     .. Local Arrays ..
      INTEGER            IMULADDTMP( 4 ), ITMP( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           PB_JUMP, PB_SETRAN
*     ..
*     .. Executable Statements ..
*
*     Compute and set the value of IRAN corresponding to A( IA, JA )
*
      ITMP( 1 ) = SEED
      ITMP( 2 ) = 0
*
      CALL PB_JUMP( JMP( JMP_1 ), IMULADD( 1, JMP_1 ), ITMP, IRAN,
     $              IMULADDTMP )
*
*     Jump ILOCBLK blocks of rows + ILOCOFF rows
*
      CALL PB_JUMP( ILOCOFF, IMULADD( 1, JMP_ROW ), IRAN, ITMP,
     $              IMULADDTMP )
      IF( MYRDIST.GT.0 ) THEN
         CALL PB_JUMP( JMP( JMP_IMBV ), IMULADD( 1, JMP_ROW  ), ITMP,
     $                 IRAN, IMULADDTMP )
         CALL PB_JUMP( MYRDIST - 1,     IMULADD( 1, JMP_MB   ), IRAN,
     $                 ITMP, IMULADDTMP )
         CALL PB_JUMP( ILOCBLK,         IMULADD( 1, JMP_NPMB ), ITMP,
     $                 IRAN, IMULADDTMP )
      ELSE
         IF( ILOCBLK.GT.0 ) THEN
            CALL PB_JUMP( JMP( JMP_IMBV ), IMULADD( 1, JMP_ROW  ), ITMP,
     $                    IRAN, IMULADDTMP )
            CALL PB_JUMP( NPROW - 1,       IMULADD( 1, JMP_MB   ), IRAN,
     $                    ITMP, IMULADDTMP )
            CALL PB_JUMP( ILOCBLK - 1,     IMULADD( 1, JMP_NPMB ), ITMP,
     $                    IRAN, IMULADDTMP )
         ELSE
            CALL PB_JUMP( 0,               IMULADD( 1, JMP_1    ), ITMP,
     $                    IRAN, IMULADDTMP )
         END IF
      END IF
*
*     Jump JLOCBLK blocks of columns + JLOCOFF columns
*
      CALL PB_JUMP( JLOCOFF, IMULADD( 1, JMP_COL ), IRAN, ITMP,
     $              IMULADDTMP )
      IF( MYCDIST.GT.0 ) THEN
         CALL PB_JUMP( JMP( JMP_INBV ), IMULADD( 1, JMP_COL  ), ITMP,
     $                 IRAN, IMULADDTMP )
         CALL PB_JUMP( MYCDIST - 1,     IMULADD( 1, JMP_NB   ), IRAN,
     $                 ITMP, IMULADDTMP )
         CALL PB_JUMP( JLOCBLK,         IMULADD( 1, JMP_NQNB ), ITMP,
     $                 IRAN, IMULADDTMP )
      ELSE
         IF( JLOCBLK.GT.0 ) THEN
            CALL PB_JUMP( JMP( JMP_INBV ), IMULADD( 1, JMP_COL  ), ITMP,
     $                    IRAN, IMULADDTMP )
            CALL PB_JUMP( NPCOL - 1,       IMULADD( 1, JMP_NB   ), IRAN,
     $                    ITMP, IMULADDTMP )
            CALL PB_JUMP( JLOCBLK - 1,     IMULADD( 1, JMP_NQNB ), ITMP,
     $                    IRAN, IMULADDTMP )
         ELSE
            CALL PB_JUMP( 0,               IMULADD( 1, JMP_1    ), ITMP,
     $                    IRAN, IMULADDTMP )
         END IF
      END IF
*
      CALL PB_SETRAN( IRAN, IMULADD( 1, JMP_1 ) )
*
      RETURN
*
*     End of PB_SETLOCRAN
*
      END
      SUBROUTINE PB_LADD( J, K, I )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Array Arguments ..
      INTEGER            I( 2 ), J( 2 ), K( 2 )
*     ..
*
*  Purpose
*  =======
*
*  PB_LADD adds without carry two long positive integers K and J and put
*  the result into I.  The long integers  I, J, K are encoded on 31 bits
*  using an array of 2 integers.  The 16-lower bits  are stored  in  the
*  first entry of each array, the  15-higher  bits  in the second entry.
*  For efficiency purposes, the intrisic modulo function is inlined.
*
*  Arguments
*  =========
*
*  J       (local input) INTEGER array
*          On entry, J is an array of dimension 2 containing the encoded
*          long integer J.
*
*  K       (local input) INTEGER array
*          On entry, K is an array of dimension 2 containing the encoded
*          long integer K.
*
*  I       (local output) INTEGER array
*          On entry, I is an array of dimension 2. On exit,  this  array
*          contains the encoded long integer I.
*
*  Further Details
*  ===============
*
*            K( 2 )   K( 1 )
*          0XXXXXXX XXXXXXXX  K   I( 1 ) = MOD( K( 1 ) + J( 1 ), 2**16 )
*        +                        carry  = ( K( 1 ) + J( 1 ) ) / 2**16
*            J( 2 )   J( 1 )
*          0XXXXXXX XXXXXXXX  J   I( 2 ) = K( 2 ) + J( 2 ) + carry
*        ----------------------   I( 2 ) = MOD( I( 2 ), 2**15 )
*            I( 2 )   I( 1 )
*          0XXXXXXX XXXXXXXX  I
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            IPOW15, IPOW16
      PARAMETER          ( IPOW15 = 2**15, IPOW16 = 2**16 )
*     ..
*     .. Local Scalars ..
      INTEGER            ITMP1, ITMP2
*     ..
*     .. Executable Statements ..
*
*     I( 1 ) = MOD( K( 1 ) + J( 1 ), IPOW16 )
*
      ITMP1 = K( 1 ) + J( 1 )
      ITMP2 = ITMP1 / IPOW16
      I( 1 ) = ITMP1 - ITMP2 * IPOW16
*
*     I( 2 ) = MOD( ( K( 1 ) + J( 1 ) ) / IPOW16 + K( 2 ) + J( 2 ),
*                   IPOW15 )
*
      ITMP1 = ITMP2 + K( 2 ) + J( 2 )
      ITMP2 = ITMP1 / IPOW15
      I( 2 ) = ITMP1 - ITMP2 * IPOW15
*
      RETURN
*
*     End of PB_LADD
*
      END
      SUBROUTINE PB_LMUL( K, J, I )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Array Arguments ..
      INTEGER            I( 2 ), J( 2 ), K( 2 )
*     ..
*
*  Purpose
*  =======
*
*  PB_LMUL  multiplies  without carry two long positive integers K and J
*  and put the result into I.  The long integers  I, J, K are encoded on
*  31 bits using an array of 2 integers. The 16-lower bits are stored in
*  the first entry of each array, the 15-higher bits in the second entry
*  of each array. For efficiency purposes, the  intrisic modulo function
*  is inlined.
*
*  Arguments
*  =========
*
*  K       (local input) INTEGER array
*          On entry, K is an array of dimension 2 containing the encoded
*          long integer K.
*
*  J       (local input) INTEGER array
*          On entry, J is an array of dimension 2 containing the encoded
*          long integer J.
*
*  I       (local output) INTEGER array
*          On entry, I is an array of dimension 2. On exit,  this  array
*          contains the encoded long integer I.
*
*  Further Details
*  ===============
*
*            K( 2 )   K( 1 )
*          0XXXXXXX XXXXXXXX  K   I( 1 ) = MOD( K( 1 ) + J( 1 ), 2**16 )
*        *                        carry  = ( K( 1 ) + J( 1 ) ) / 2**16
*            J( 2 )   J( 1 )
*          0XXXXXXX XXXXXXXX  J   I( 2 ) = K( 2 ) + J( 2 ) + carry
*        ----------------------   I( 2 ) = MOD( I( 2 ), 2**15 )
*            I( 2 )   I( 1 )
*          0XXXXXXX XXXXXXXX  I
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            IPOW15, IPOW16, IPOW30
      PARAMETER          ( IPOW15 = 2**15, IPOW16 = 2**16,
     $                   IPOW30 = 2**30 )
*     ..
*     .. Local Scalars ..
      INTEGER            ITMP1, ITMP2
*     ..
*     .. Executable Statements ..
*
      ITMP1 = K( 1 ) * J( 1 )
      IF( ITMP1.LT.0 )
     $   ITMP1 = ( ITMP1 + IPOW30 ) + IPOW30
*
*     I( 1 ) = MOD( ITMP1, IPOW16 )
*
      ITMP2 = ITMP1 / IPOW16
      I( 1 ) = ITMP1 - ITMP2 * IPOW16
*
      ITMP1 = K( 1 ) * J( 2 ) + K( 2 ) * J( 1 )
      IF( ITMP1.LT.0 )
     $   ITMP1 = ( ITMP1 + IPOW30 ) + IPOW30
*
      ITMP1 = ITMP2 + ITMP1
      IF( ITMP1.LT.0 )
     $   ITMP1 = ( ITMP1 + IPOW30 ) + IPOW30
*
*     I( 2 ) = MOD( ITMP1, IPOW15 )
*
      I( 2 ) = ITMP1 - ( ITMP1 / IPOW15 ) * IPOW15
*
      RETURN
*
*     End of PB_LMUL
*
      END
      SUBROUTINE PB_JUMP( K, MULADD, IRANN, IRANM, IMA )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            K
*     ..
*     .. Array Arguments ..
      INTEGER            IMA( 4 ), IRANM( 2 ), IRANN( 2 ), MULADD( 4 )
*     ..
*
*  Purpose
*  =======
*
*  PB_JUMP  computes the constants A and C to jump K numbers in the ran-
*  dom sequence:
*
*     X( n+K ) = A * X( n ) + C.
*
*  The constants encoded in MULADD specify how to jump from entry in the
*  sequence to the next.
*
*  Arguments
*  =========
*
*  K       (local input) INTEGER
*          On entry, K specifies the number of entries  of  the sequence
*          to jump over. When K is less or equal than zero, A and C  are
*          not computed, and  IRANM  is set to  IRANN corresponding to a
*          jump of size zero.
*
*  MULADD  (local input) INTEGER array
*          On entry,  MULADD  is an  array of dimension 4 containing the
*          encoded constants a and c to  jump  from  X( n ) to  X( n+1 )
*          ( = a*X( n )+c) in the random sequence.  MULADD(1:2) contains
*          respectively the 16-lower and 16-higher bits of  the constant
*          a,  and  MULADD(3:4) contains the 16-lower and 16-higher bits
*          of the constant c.
*
*  IRANN   (local input) INTEGER array
*          On entry,  IRANN  is an array of dimension 2. This array con-
*          tains respectively the 16-lower and 16-higher bits of the en-
*          coding of X( n ).
*
*  IRANM   (local output) INTEGER array
*          On entry,  IRANM  is an  array of dimension 2.  On exit, this
*          array contains respectively the 16-lower and  16-higher  bits
*          of the encoding of X( n+K ).
*
*  IMA     (local output) INTEGER array
*          On entry, IMA is an array of dimension 4. On exit, when K is
*          greater than zero, this array contains the encoded constants
*          A and C to  jump  from X( n ) to  X( n+K ) in the random se-
*          quence.  IMA(1:2)  contains  respectively  the  16-lower and
*          16-higher bits of the constant A, and IMA(3:4)  contains the
*          16-lower  and  16-higher  bits of the constant  C. When K is
*          less or equal than zero, this array is not referenced.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. Local Arrays ..
      INTEGER            J( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           PB_LADD, PB_LMUL
*     ..
*     .. Executable Statements ..
*
      IF( K.GT.0 ) THEN
*
         IMA( 1 ) = MULADD( 1 )
         IMA( 2 ) = MULADD( 2 )
         IMA( 3 ) = MULADD( 3 )
         IMA( 4 ) = MULADD( 4 )
*
         DO 10 I = 1, K - 1
*
            CALL PB_LMUL( IMA, MULADD, J )
*
            IMA( 1 ) = J( 1 )
            IMA( 2 ) = J( 2 )
*
            CALL PB_LMUL( IMA( 3 ), MULADD, J )
            CALL PB_LADD( MULADD( 3 ), J, IMA( 3 ) )
*
   10    CONTINUE
*
         CALL PB_LMUL( IRANN, IMA, J )
         CALL PB_LADD( J, IMA( 3 ), IRANM )
*
      ELSE
*
         IRANM( 1 ) = IRANN( 1 )
         IRANM( 2 ) = IRANN( 2 )
*
      END IF
*
      RETURN
*
*     End of PB_JUMP
*
      END
      SUBROUTINE PB_SETRAN( IRAN, IAC )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Array Arguments ..
      INTEGER            IAC( 4 ), IRAN( 2 )
*     ..
*
*  Purpose
*  =======
*
*  PB_SETRAN  initializes  the random generator with the encoding of the
*  first number X( 1 ) in the sequence,  and  the constants a and c used
*  to compute the next element in the sequence:
*
*     X( n+1 ) = a * X( n ) + c.
*
*  X( 1 ), a and c are stored in the common block  RANCOM  for later use
*  (see the routines PB_SRAN or PB_DRAN).
*
*  Arguments
*  =========
*
*  IRAN    (local input) INTEGER array
*          On entry, IRAN is an array of dimension 2.  This  array  con-
*          tains respectively the 16-lower and 16-higher bits of the en-
*          coding of X( 1 ).
*
*  IAC     (local input) INTEGER array
*          On entry,  IAC  is an array of dimension 4.  IAC(1:2) contain
*          respectively the 16-lower and 16-higher bits  of the constant
*          a, and  IAC(3:4)  contain  the 16-lower and 16-higher bits of
*          the constant c.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Common Blocks ..
      INTEGER            IACS( 4 ), IRAND( 2 )
      COMMON             /RANCOM/ IRAND, IACS
*     ..
*     .. Save Statements ..
      SAVE               /RANCOM/
*     ..
*     .. Executable Statements ..
*
      IRAND( 1 ) = IRAN( 1 )
      IRAND( 2 ) = IRAN( 2 )
      IACS( 1 )  = IAC( 1 )
      IACS( 2 )  = IAC( 2 )
      IACS( 3 )  = IAC( 3 )
      IACS( 4 )  = IAC( 4 )
*
      RETURN
*
*     End of PB_SETRAN
*
      END
      SUBROUTINE PB_JUMPIT( MULADD, IRANN, IRANM )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Array Arguments ..
      INTEGER            IRANM( 2 ), IRANN( 2 ), MULADD( 4 )
*     ..
*
*  Purpose
*  =======
*
*  PB_JUMPIT  jumps  in the random sequence from the number X( n ) enco-
*  ded in IRANN to the number  X( m )  encoded in  IRANM using the cons-
*  tants A and C encoded in MULADD:
*
*     X( m ) = A * X( n ) + C.
*
*  The constants A and C obviously depend on m and n, see the subroutine
*  PB_JUMP in order to set them up.
*
*  Arguments
*  =========
*
*  MULADD  (local input) INTEGER array
*          On netry, MULADD is an array of dimension 4. MULADD(1:2) con-
*          tains  respectively  the 16-lower and 16-higher bits  of  the
*          constant  A,  and   MULADD(3:4)  contains  the  16-lower  and
*          16-higher bits of the constant C.
*
*  IRANN   (local input) INTEGER array
*          On entry,  IRANN  is an array of dimension 2. This array con-
*          tains respectively the 16-lower and 16-higher bits of the en-
*          coding of X( n ).
*
*  IRANM   (local output) INTEGER array
*          On entry,  IRANM  is an  array of dimension 2.  On exit, this
*          array contains respectively the 16-lower and  16-higher  bits
*          of the encoding of X( m ).
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Arrays ..
      INTEGER            J( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           PB_LADD, PB_LMUL
*     ..
*     .. Common Blocks ..
      INTEGER            IACS( 4 ), IRAND( 2 )
      COMMON             /RANCOM/ IRAND, IACS
*     ..
*     .. Save Statements ..
      SAVE               /RANCOM/
*     ..
*     .. Executable Statements ..
*
      CALL PB_LMUL( IRANN, MULADD, J )
      CALL PB_LADD( J, MULADD( 3 ), IRANM )
*
      IRAND( 1 ) = IRANM( 1 )
      IRAND( 2 ) = IRANM( 2 )
*
      RETURN
*
*     End of PB_JUMPIT
*
      END
