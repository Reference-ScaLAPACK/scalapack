      SUBROUTINE PZOPTEE( ICTXT, NOUT, SUBPTR, SCODE, SNAME )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, NOUT, SCODE
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      SNAME
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL           SUBPTR
*     ..
*
*  Purpose
*  =======
*
*  PZOPTEE  tests  whether  the  PBLAS respond correctly to a bad option
*  argument.
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
*  SUBPTR  (global input) SUBROUTINE
*          On entry,  SUBPTR  is  a  subroutine. SUBPTR must be declared
*          EXTERNAL in the calling subroutine.
*
*  SCODE   (global input) INTEGER
*          On entry, SCODE specifies the calling sequence code.
*
*  SNAME   (global input) CHARACTER*(*)
*          On entry,  SNAME  specifies  the subroutine name calling this
*          subprogram.
*
*  Calling sequence encodings
*  ==========================
*
*  code Formal argument list                                Examples
*
*  11   (n,      v1,v2)                                     _SWAP, _COPY
*  12   (n,s1,   v1   )                                     _SCAL, _SCAL
*  13   (n,s1,   v1,v2)                                     _AXPY, _DOT_
*  14   (n,s1,i1,v1   )                                     _AMAX
*  15   (n,u1,   v1   )                                     _ASUM, _NRM2
*
*  21   (     trans,     m,n,s1,m1,v1,s2,v2)                _GEMV
*  22   (uplo,             n,s1,m1,v1,s2,v2)                _SYMV, _HEMV
*  23   (uplo,trans,diag,  n,   m1,v1      )                _TRMV, _TRSV
*  24   (                m,n,s1,v1,v2,m1)                   _GER_
*  25   (uplo,             n,s1,v1,   m1)                   _SYR
*  26   (uplo,             n,u1,v1,   m1)                   _HER
*  27   (uplo,             n,s1,v1,v2,m1)                   _SYR2, _HER2
*
*  31   (          transa,transb,     m,n,k,s1,m1,m2,s2,m3) _GEMM
*  32   (side,uplo,                   m,n,  s1,m1,m2,s2,m3) _SYMM, _HEMM
*  33   (     uplo,trans,               n,k,s1,m1,   s2,m3) _SYRK
*  34   (     uplo,trans,               n,k,u1,m1,   u2,m3) _HERK
*  35   (     uplo,trans,               n,k,s1,m1,m2,s2,m3) _SYR2K
*  36   (     uplo,trans,               n,k,s1,m1,m2,u2,m3) _HER2K
*  37   (                             m,n,  s1,m1,   s2,m3) _TRAN_
*  38   (side,uplo,transa,       diag,m,n,  s1,m1,m2      ) _TRMM, _TRSM
*  39   (          trans,             m,n,  s1,m1,   s2,m3) _GEADD
*  40   (     uplo,trans,             m,n,  s1,m1,   s2,m3) _TRADD
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER             APOS
*     ..
*     .. External Subroutines ..
      EXTERNAL            PZCHKOPT
*     ..
*     .. Executable Statements ..
*
*     Level 2 PBLAS
*
      IF( SCODE.EQ.21 ) THEN
*
*        Check 1st (and only) option
*
         APOS = 1
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
      ELSE IF( SCODE.EQ.22 .OR. SCODE.EQ.25 .OR. SCODE.EQ.26 .OR.
     $         SCODE.EQ.27 ) THEN
*
*        Check 1st (and only) option
*
         APOS = 1
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'U', APOS )
*
      ELSE IF( SCODE.EQ.23 ) THEN
*
*        Check 1st option
*
         APOS = 1
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'U', APOS )
*
*        Check 2nd option
*
         APOS = 2
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
*        Check 3rd option
*
         APOS = 3
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'D', APOS )
*
*     Level 3 PBLAS
*
      ELSE IF( SCODE.EQ.31 ) THEN
*
*        Check 1st option
*
         APOS = 1
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
*        Check 2'nd option
*
         APOS = 2
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'B', APOS )
*
      ELSE IF( SCODE.EQ.32 ) THEN
*
*        Check 1st option
*
         APOS = 1
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'S', APOS )
*
*        Check 2nd option
*
         APOS = 2
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'U', APOS )
*
      ELSE IF( SCODE.EQ.33 .OR. SCODE.EQ.34 .OR. SCODE.EQ.35 .OR.
     $         SCODE.EQ.36 .OR. SCODE.EQ.40 ) THEN
*
*        Check 1st option
*
         APOS = 1
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'U', APOS )
*
*        Check 2'nd option
*
         APOS = 2
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
      ELSE IF( SCODE.EQ.38 ) THEN
*
*        Check 1st option
*
         APOS = 1
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'S', APOS )
*
*        Check 2nd option
*
         APOS = 2
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'U', APOS )
*
*        Check 3rd option
*
         APOS = 3
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
*        Check 4th option
*
         APOS = 4
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'D', APOS )
*
*
      ELSE IF( SCODE.EQ.39 ) THEN
*
*        Check 1st option
*
         APOS = 1
         CALL PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
      END IF
*
      RETURN
*
*     End of PZOPTEE
*
      END
      SUBROUTINE PZCHKOPT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, ARGNAM,
     $                     ARGPOS )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1         ARGNAM
      INTEGER             ARGPOS, ICTXT, NOUT, SCODE
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)       SNAME
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL            SUBPTR
*     ..
*
*  Purpose
*  =======
*
*  PZCHKOPT tests the option ARGNAM in any PBLAS routine.
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
*  SUBPTR  (global input) SUBROUTINE
*          On entry,  SUBPTR  is  a  subroutine. SUBPTR must be declared
*          EXTERNAL in the calling subroutine.
*
*  SCODE   (global input) INTEGER
*          On entry, SCODE specifies the calling sequence code.
*
*  SNAME   (global input) CHARACTER*(*)
*          On entry,  SNAME  specifies  the subroutine name calling this
*          subprogram.
*
*  ARGNAM  (global input) CHARACTER*(*)
*          On entry,  ARGNAM  specifies  the  name  of  the option to be
*          checked. ARGNAM can either be 'D', 'S', 'A', 'B', or 'U'.
*
*  ARGPOS  (global input) INTEGER
*          On entry, ARGPOS indicates the position of the option ARGNAM
*          to be tested.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            INFOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           PCHKPBE, PZCALLSUB, PZSETPBLAS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Common Blocks ..
      CHARACTER          DIAG, SIDE, TRANSA, TRANSB, UPLO
      COMMON             /PBLASC/DIAG, SIDE, TRANSA, TRANSB, UPLO
*     ..
*     .. Executable Statements ..
*
*     Reiniatilize the dummy arguments to correct values
*
      CALL PZSETPBLAS( ICTXT )
*
      IF( LSAME( ARGNAM, 'D' ) ) THEN
*
*        Generate bad DIAG option
*
         DIAG = '/'
*
      ELSE IF( LSAME( ARGNAM, 'S' ) ) THEN
*
*        Generate bad SIDE option
*
         SIDE = '/'
*
      ELSE IF( LSAME( ARGNAM, 'A' ) ) THEN
*
*        Generate bad TRANSA option
*
         TRANSA = '/'
*
      ELSE IF( LSAME( ARGNAM, 'B' ) ) THEN
*
*        Generate bad TRANSB option
*
         TRANSB = '/'
*
      ELSE IF( LSAME( ARGNAM, 'U' ) ) THEN
*
*        Generate bad UPLO option
*
         UPLO = '/'
*
      END IF
*
*     Set INFOT to the position of the bad dimension argument
*
      INFOT = ARGPOS
*
*     Call the PBLAS routine
*
      CALL PZCALLSUB( SUBPTR, SCODE )
      CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
      RETURN
*
*     End of PZCHKOPT
*
      END
      SUBROUTINE PZDIMEE( ICTXT, NOUT, SUBPTR, SCODE, SNAME )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, NOUT, SCODE
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      SNAME
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL           SUBPTR
*     ..
*
*  Purpose
*  =======
*
*  PZDIMEE  tests whether the PBLAS respond correctly to a bad dimension
*  argument.
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
*  SUBPTR  (global input) SUBROUTINE
*          On entry,  SUBPTR  is  a  subroutine. SUBPTR must be declared
*          EXTERNAL in the calling subroutine.
*
*  SCODE   (global input) INTEGER
*          On entry, SCODE specifies the calling sequence code.
*
*  SNAME   (global input) CHARACTER*(*)
*          On entry,  SNAME  specifies  the subroutine name calling this
*          subprogram.
*
*  Calling sequence encodings
*  ==========================
*
*  code Formal argument list                                Examples
*
*  11   (n,      v1,v2)                                     _SWAP, _COPY
*  12   (n,s1,   v1   )                                     _SCAL, _SCAL
*  13   (n,s1,   v1,v2)                                     _AXPY, _DOT_
*  14   (n,s1,i1,v1   )                                     _AMAX
*  15   (n,u1,   v1   )                                     _ASUM, _NRM2
*
*  21   (     trans,     m,n,s1,m1,v1,s2,v2)                _GEMV
*  22   (uplo,             n,s1,m1,v1,s2,v2)                _SYMV, _HEMV
*  23   (uplo,trans,diag,  n,   m1,v1      )                _TRMV, _TRSV
*  24   (                m,n,s1,v1,v2,m1)                   _GER_
*  25   (uplo,             n,s1,v1,   m1)                   _SYR
*  26   (uplo,             n,u1,v1,   m1)                   _HER
*  27   (uplo,             n,s1,v1,v2,m1)                   _SYR2, _HER2
*
*  31   (          transa,transb,     m,n,k,s1,m1,m2,s2,m3) _GEMM
*  32   (side,uplo,                   m,n,  s1,m1,m2,s2,m3) _SYMM, _HEMM
*  33   (     uplo,trans,               n,k,s1,m1,   s2,m3) _SYRK
*  34   (     uplo,trans,               n,k,u1,m1,   u2,m3) _HERK
*  35   (     uplo,trans,               n,k,s1,m1,m2,s2,m3) _SYR2K
*  36   (     uplo,trans,               n,k,s1,m1,m2,u2,m3) _HER2K
*  37   (                             m,n,  s1,m1,   s2,m3) _TRAN_
*  38   (side,uplo,transa,       diag,m,n,  s1,m1,m2      ) _TRMM, _TRSM
*  39   (          trans,             m,n,  s1,m1,   s2,m3) _GEADD
*  40   (     uplo,trans,             m,n,  s1,m1,   s2,m3) _TRADD
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER             APOS
*     ..
*     .. External Subroutines ..
      EXTERNAL            PZCHKDIM
*     ..
*     .. Executable Statements ..
*
*     Level 1 PBLAS
*
      IF( SCODE.EQ.11 .OR. SCODE.EQ.12 .OR. SCODE.EQ.13 .OR.
     $    SCODE.EQ.14 .OR. SCODE.EQ.15 ) THEN
*
*        Check 1st (and only) dimension
*
         APOS = 1
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'N', APOS )
*
*     Level 2 PBLAS
*
      ELSE IF( SCODE.EQ.21 ) THEN
*
*        Check 1st dimension
*
         APOS = 2
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'M', APOS )
*
*        Check 2nd dimension
*
         APOS = 3
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'N', APOS )
*
      ELSE IF( SCODE.EQ.22 .OR. SCODE.EQ.25 .OR. SCODE.EQ.26 .OR.
     $         SCODE.EQ.27 ) THEN
*
*        Check 1st (and only) dimension
*
         APOS = 2
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'N', APOS )
*
      ELSE IF( SCODE.EQ.23 ) THEN
*
*        Check 1st (and only) dimension
*
         APOS = 4
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'N', APOS )
*
      ELSE IF( SCODE.EQ.24 ) THEN
*
*        Check 1st dimension
*
         APOS = 1
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'M', APOS )
*
*        Check 2nd dimension
*
         APOS = 2
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'N', APOS )
*
*     Level 3 PBLAS
*
      ELSE IF( SCODE.EQ.31 ) THEN
*
*        Check 1st dimension
*
         APOS = 3
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'M', APOS )
*
*        Check 2nd dimension
*
         APOS = 4
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'N', APOS )
*
*        Check 3rd dimension
*
         APOS = 5
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'K', APOS )
*
      ELSE IF( SCODE.EQ.32 ) THEN
*
*        Check 1st dimension
*
         APOS = 3
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'M', APOS )
*
*        Check 2nd dimension
*
         APOS = 4
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'N', APOS )
*
      ELSE IF( SCODE.EQ.33 .OR. SCODE.EQ.34 .OR. SCODE.EQ.35 .OR.
     $         SCODE.EQ.36 ) THEN
*
*        Check 1st dimension
*
         APOS = 3
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'N', APOS )
*
*        Check 2nd dimension
*
         APOS = 4
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'K', APOS )
*
      ELSE IF( SCODE.EQ.37 ) THEN
*
*        Check 1st dimension
*
         APOS = 1
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'M', APOS )
*
*        Check 2nd dimension
*
         APOS = 2
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'N', APOS )
*
      ELSE IF( SCODE.EQ.38 ) THEN
*
*        Check 1st dimension
*
         APOS = 5
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'M', APOS )
*
*        Check 2nd dimension
*
         APOS = 6
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'N', APOS )
*
      ELSE IF( SCODE.EQ.39 ) THEN
*
*        Check 1st dimension
*
         APOS = 2
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'M', APOS )
*
*        Check 2nd dimension
*
         APOS = 3
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'N', APOS )
*
      ELSE IF( SCODE.EQ.40 ) THEN
*
*        Check 1st dimension
*
         APOS = 3
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'M', APOS )
*
*        Check 2nd dimension
*
         APOS = 4
         CALL PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'N', APOS )
*
      END IF
*
      RETURN
*
*     End of PZDIMEE
*
      END
      SUBROUTINE PZCHKDIM( ICTXT, NOUT, SUBPTR, SCODE, SNAME, ARGNAM,
     $                     ARGPOS )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1         ARGNAM
      INTEGER             ARGPOS, ICTXT, NOUT, SCODE
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)       SNAME
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL            SUBPTR
*     ..
*
*  Purpose
*  =======
*
*  PZCHKDIM tests the dimension ARGNAM in any PBLAS routine.
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
*  SUBPTR  (global input) SUBROUTINE
*          On entry,  SUBPTR  is  a  subroutine. SUBPTR must be declared
*          EXTERNAL in the calling subroutine.
*
*  SCODE   (global input) INTEGER
*          On entry, SCODE specifies the calling sequence code.
*
*  SNAME   (global input) CHARACTER*(*)
*          On entry,  SNAME  specifies  the subroutine name calling this
*          subprogram.
*
*  ARGNAM  (global input) CHARACTER*(*)
*          On entry,  ARGNAM  specifies  the name of the dimension to be
*          checked. ARGNAM can either be 'M', 'N' or 'K'.
*
*  ARGPOS  (global input) INTEGER
*          On entry, ARGPOS indicates the position of the option ARGNAM
*          to be tested.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            INFOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           PCHKPBE, PZCALLSUB, PZSETPBLAS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Common Blocks ..
      INTEGER            KDIM, MDIM, NDIM
      COMMON             /PBLASN/KDIM, MDIM, NDIM
*     ..
*     .. Executable Statements ..
*
*     Reiniatilize the dummy arguments to correct values
*
      CALL PZSETPBLAS( ICTXT )
*
      IF( LSAME( ARGNAM, 'M' ) ) THEN
*
*        Generate bad MDIM
*
         MDIM = -1
*
      ELSE IF( LSAME( ARGNAM, 'N' ) ) THEN
*
*        Generate bad NDIM
*
         NDIM = -1
*
      ELSE
*
*        Generate bad KDIM
*
         KDIM = -1
*
      END IF
*
*     Set INFOT to the position of the bad dimension argument
*
      INFOT = ARGPOS
*
*     Call the PBLAS routine
*
      CALL PZCALLSUB( SUBPTR, SCODE )
      CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
      RETURN
*
*     End of PZCHKDIM
*
      END
      SUBROUTINE PZVECEE( ICTXT, NOUT, SUBPTR, SCODE, SNAME )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER             ICTXT, NOUT, SCODE
*     ..
*     .. Array Arguments ..
      CHARACTER*7         SNAME
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL            SUBPTR
*     ..
*
*  Purpose
*  =======
*
*  PZVECEE  tests  whether  the  PBLAS respond correctly to a bad vector
*  argument.  Each  vector <vec> is described by: <vec>, I<vec>, J<vec>,
*  DESC<vec>,  INC<vec>.   Out   of  all  these,  only  I<vec>,  J<vec>,
*  DESC<vec>, and INC<vec> can be tested.
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
*  SUBPTR  (global input) SUBROUTINE
*          On entry,  SUBPTR  is  a  subroutine. SUBPTR must be declared
*          EXTERNAL in the calling subroutine.
*
*  SCODE   (global input) INTEGER
*          On entry, SCODE specifies the calling sequence code.
*
*  SNAME   (global input) CHARACTER*(*)
*          On entry,  SNAME  specifies  the subroutine name calling this
*          subprogram.
*
*  Calling sequence encodings
*  ==========================
*
*  code Formal argument list                                Examples
*
*  11   (n,      v1,v2)                                     _SWAP, _COPY
*  12   (n,s1,   v1   )                                     _SCAL, _SCAL
*  13   (n,s1,   v1,v2)                                     _AXPY, _DOT_
*  14   (n,s1,i1,v1   )                                     _AMAX
*  15   (n,u1,   v1   )                                     _ASUM, _NRM2
*
*  21   (     trans,     m,n,s1,m1,v1,s2,v2)                _GEMV
*  22   (uplo,             n,s1,m1,v1,s2,v2)                _SYMV, _HEMV
*  23   (uplo,trans,diag,  n,   m1,v1      )                _TRMV, _TRSV
*  24   (                m,n,s1,v1,v2,m1)                   _GER_
*  25   (uplo,             n,s1,v1,   m1)                   _SYR
*  26   (uplo,             n,u1,v1,   m1)                   _HER
*  27   (uplo,             n,s1,v1,v2,m1)                   _SYR2, _HER2
*
*  31   (          transa,transb,     m,n,k,s1,m1,m2,s2,m3) _GEMM
*  32   (side,uplo,                   m,n,  s1,m1,m2,s2,m3) _SYMM, _HEMM
*  33   (     uplo,trans,               n,k,s1,m1,   s2,m3) _SYRK
*  34   (     uplo,trans,               n,k,u1,m1,   u2,m3) _HERK
*  35   (     uplo,trans,               n,k,s1,m1,m2,s2,m3) _SYR2K
*  36   (     uplo,trans,               n,k,s1,m1,m2,u2,m3) _HER2K
*  37   (                             m,n,  s1,m1,   s2,m3) _TRAN_
*  38   (side,uplo,transa,       diag,m,n,  s1,m1,m2      ) _TRMM, _TRSM
*  39   (          trans,             m,n,  s1,m1,   s2,m3) _GEADD
*  40   (     uplo,trans,             m,n,  s1,m1,   s2,m3) _TRADD
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER             APOS
*     ..
*     .. External Subroutines ..
      EXTERNAL            PZCHKMAT
*     ..
*     .. Executable Statements ..
*
*     Level 1 PBLAS
*
      IF( SCODE.EQ.11 ) THEN
*
*        Check 1st vector
*
         APOS = 2
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'X', APOS )
*
*        Check 2nd vector
*
         APOS = 7
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'Y', APOS )
*
      ELSE IF( SCODE.EQ.12 .OR. SCODE.EQ.15 ) THEN
*
*        Check 1st (and only) vector
*
         APOS = 3
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'X', APOS )
*
      ELSE IF( SCODE.EQ.13 ) THEN
*
*        Check 1st vector
*
         APOS = 3
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'X', APOS )
*
*        Check 2nd vector
*
         APOS = 8
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'Y', APOS )
*
      ELSE IF( SCODE.EQ.14 ) THEN
*
*        Check 1st (and only) vector
*
         APOS = 4
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'X', APOS )
*
*     Level 2 PBLAS
*
      ELSE IF( SCODE.EQ.21 ) THEN
*
*        Check 1st vector
*
         APOS = 9
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'X', APOS )
*
*        Check 2nd vector
*
         APOS = 15
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'Y', APOS )
*
      ELSE IF( SCODE.EQ.22 ) THEN
*
*        Check 1st vector
*
         APOS = 8
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'X', APOS )
*
*        Check 2nd vector
*
         APOS = 14
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'Y', APOS )
*
      ELSE IF( SCODE.EQ.23 ) THEN
*
*        Check 1st (and only) vector
*
         APOS = 9
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'X', APOS )
*
      ELSE IF( SCODE.EQ.24 .OR. SCODE.EQ.27 ) THEN
*
*        Check 1st vector
*
         APOS = 4
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'X', APOS )
*
*        Check 2nd vector
*
         APOS = 9
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'Y', APOS )
*
      ELSE IF( SCODE.EQ.26 .OR. SCODE.EQ.27 ) THEN
*
*        Check 1'st (and only) vector
*
         APOS = 4
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'X', APOS )
*
      END IF
*
      RETURN
*
*     End of PZVECEE
*
      END
      SUBROUTINE PZMATEE( ICTXT, NOUT, SUBPTR, SCODE, SNAME )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER             ICTXT, NOUT, SCODE
*     ..
*     .. Array Arguments ..
      CHARACTER*7         SNAME
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL            SUBPTR
*     ..
*
*  Purpose
*  =======
*
*  PZMATEE  tests  whether  the  PBLAS respond correctly to a bad matrix
*  argument.  Each  matrix <mat> is described by: <mat>, I<mat>, J<mat>,
*  and DESC<mat>.  Out  of  all these, only I<vec>, J<vec> and DESC<mat>
*  can be tested.
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
*  SUBPTR  (global input) SUBROUTINE
*          On entry,  SUBPTR  is  a  subroutine. SUBPTR must be declared
*          EXTERNAL in the calling subroutine.
*
*  SCODE   (global input) INTEGER
*          On entry, SCODE specifies the calling sequence code.
*
*  SNAME   (global input) CHARACTER*(*)
*          On entry,  SNAME  specifies  the subroutine name calling this
*          subprogram.
*
*  Calling sequence encodings
*  ==========================
*
*  code Formal argument list                                Examples
*
*  11   (n,      v1,v2)                                     _SWAP, _COPY
*  12   (n,s1,   v1   )                                     _SCAL, _SCAL
*  13   (n,s1,   v1,v2)                                     _AXPY, _DOT_
*  14   (n,s1,i1,v1   )                                     _AMAX
*  15   (n,u1,   v1   )                                     _ASUM, _NRM2
*
*  21   (     trans,     m,n,s1,m1,v1,s2,v2)                _GEMV
*  22   (uplo,             n,s1,m1,v1,s2,v2)                _SYMV, _HEMV
*  23   (uplo,trans,diag,  n,   m1,v1      )                _TRMV, _TRSV
*  24   (                m,n,s1,v1,v2,m1)                   _GER_
*  25   (uplo,             n,s1,v1,   m1)                   _SYR
*  26   (uplo,             n,u1,v1,   m1)                   _HER
*  27   (uplo,             n,s1,v1,v2,m1)                   _SYR2, _HER2
*
*  31   (          transa,transb,     m,n,k,s1,m1,m2,s2,m3) _GEMM
*  32   (side,uplo,                   m,n,  s1,m1,m2,s2,m3) _SYMM, _HEMM
*  33   (     uplo,trans,               n,k,s1,m1,   s2,m3) _SYRK
*  34   (     uplo,trans,               n,k,u1,m1,   u2,m3) _HERK
*  35   (     uplo,trans,               n,k,s1,m1,m2,s2,m3) _SYR2K
*  36   (     uplo,trans,               n,k,s1,m1,m2,u2,m3) _HER2K
*  37   (                             m,n,  s1,m1,   s2,m3) _TRAN_
*  38   (side,uplo,transa,       diag,m,n,  s1,m1,m2      ) _TRMM, _TRSM
*  39   (          trans,             m,n,  s1,m1,   s2,m3) _GEADD
*  40   (     uplo,trans,             m,n,  s1,m1,   s2,m3) _TRADD
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER             APOS
*     ..
*     .. External Subroutines ..
      EXTERNAL            PZCHKMAT
*     ..
*     .. Executable Statements ..
*
*     Level 2 PBLAS
*
      IF( SCODE.EQ.21 .OR. SCODE.EQ.23 ) THEN
*
*        Check 1st (and only) matrix
*
         APOS = 5
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
      ELSE IF( SCODE.EQ.22 ) THEN
*
*        Check 1st (and only) matrix
*
         APOS = 4
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
      ELSE IF( SCODE.EQ.24 .OR. SCODE.EQ.27 ) THEN
*
*        Check 1st (and only) matrix
*
         APOS = 14
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
      ELSE IF( SCODE.EQ.25 .OR. SCODE.EQ.26 ) THEN
*
*        Check 1st (and only) matrix
*
         APOS = 9
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
*     Level 3 PBLAS
*
      ELSE IF( SCODE.EQ.31 ) THEN
*
*        Check 1st matrix
*
         APOS = 7
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
*        Check 2nd matrix
*
         APOS = 11
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'B', APOS )
*
*        Check 3nd matrix
*
         APOS = 16
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'C', APOS )
*
      ELSE IF( SCODE.EQ.32 .OR. SCODE.EQ.35 .OR. SCODE.EQ.36 ) THEN
*
*        Check 1st matrix
*
         APOS = 6
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
*        Check 2nd matrix
*
         APOS = 10
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'B', APOS )
*
*        Check 3nd matrix
*
         APOS = 15
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'C', APOS )
*
      ELSE IF( SCODE.EQ.33 .OR. SCODE.EQ.34 ) THEN
*
*        Check 1st matrix
*
         APOS = 6
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
*        Check 2nd matrix
*
         APOS = 11
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'C', APOS )
*
      ELSE IF( SCODE.EQ.37 ) THEN
*
*        Check 1st matrix
*
         APOS = 4
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
*        Check 2nd matrix
*
         APOS = 9
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'C', APOS )
*
      ELSE IF( SCODE.EQ.38 ) THEN
*
*        Check 1st matrix
*
         APOS = 8
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
*        Check 2nd matrix
*
         APOS = 12
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'B', APOS )
*
      ELSE IF( SCODE.EQ.39 ) THEN
*
*        Check 1st matrix
*
         APOS = 5
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
*        Check 2nd matrix
*
         APOS = 10
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'C', APOS )
*
      ELSE IF( SCODE.EQ.40 ) THEN
*
*        Check 1st matrix
*
         APOS = 6
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'A', APOS )
*
*        Check 2nd matrix
*
         APOS = 11
         CALL PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, 'C', APOS )
*
      END IF
*
      RETURN
*
*     End of PZMATEE
*
      END
      SUBROUTINE PZSETPBLAS( ICTXT )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT
*     ..
*
*  Purpose
*  =======
*
*  PZSETPBLAS initializes *all* the dummy arguments to correct values.
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
      DOUBLE PRECISION   RONE
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   RONE = 1.0D+0 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           PB_DESCSET2
*     ..
*     .. Common Blocks ..
      CHARACTER*1        DIAG, SIDE, TRANSA, TRANSB, UPLO
      INTEGER            IA, IB, IC, INCX, INCY, ISCLR, IX, IY, JA, JB,
     $                   JC, JX, JY, KDIM, MDIM, NDIM
      DOUBLE PRECISION   USCLR
      COMPLEX*16         SCLR
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), DESCC( DLEN_ ),
     $                   DESCX( DLEN_ ), DESCY( DLEN_ )
      COMPLEX*16         A( 2, 2 ), B( 2, 2 ), C( 2, 2 ), X( 2 ), Y( 2 )
      COMMON             /PBLASC/DIAG, SIDE, TRANSA, TRANSB, UPLO
      COMMON             /PBLASD/DESCA, DESCB, DESCC, DESCX, DESCY
      COMMON             /PBLASI/IA, IB, IC, INCX, INCY, ISCLR, IX, IY,
     $                   JA, JB, JC, JX, JY
      COMMON             /PBLASM/A, B, C
      COMMON             /PBLASN/KDIM, MDIM, NDIM
      COMMON             /PBLASS/SCLR, USCLR
      COMMON             /PBLASV/X, Y
*     ..
*     .. Executable Statements ..
*
*     Set default values for options
*
      DIAG   = 'N'
      SIDE   = 'L'
      TRANSA = 'N'
      TRANSB = 'N'
      UPLO   = 'U'
*
*     Set default values for scalars
*
      KDIM   = 1
      MDIM   = 1
      NDIM   = 1
      ISCLR  = 1
      SCLR   = ONE
      USCLR  = RONE
*
*     Set default values for distributed matrix A
*
      A( 1, 1 ) = ONE
      A( 2, 1 ) = ONE
      A( 1, 2 ) = ONE
      A( 2, 2 ) = ONE
      IA = 1
      JA = 1
      CALL PB_DESCSET2( DESCA, 2, 2, 1, 1, 1, 1, 0, 0, ICTXT, 2 )
*
*     Set default values for distributed matrix B
*
      B( 1, 1 ) = ONE
      B( 2, 1 ) = ONE
      B( 1, 2 ) = ONE
      B( 2, 2 ) = ONE
      IB = 1
      JB = 1
      CALL PB_DESCSET2( DESCB, 2, 2, 1, 1, 1, 1, 0, 0, ICTXT, 2 )
*
*     Set default values for distributed matrix C
*
      C( 1, 1 ) = ONE
      C( 2, 1 ) = ONE
      C( 1, 2 ) = ONE
      C( 2, 2 ) = ONE
      IC = 1
      JC = 1
      CALL PB_DESCSET2( DESCC, 2, 2, 1, 1, 1, 1, 0, 0, ICTXT, 2 )
*
*     Set default values for distributed matrix X
*
      X( 1 ) = ONE
      X( 2 ) = ONE
      IX = 1
      JX = 1
      CALL PB_DESCSET2( DESCX, 2, 1, 1, 1, 1, 1, 0, 0, ICTXT, 2 )
      INCX = 1
*
*     Set default values for distributed matrix Y
*
      Y( 1 ) = ONE
      Y( 2 ) = ONE
      IY = 1
      JY = 1
      CALL PB_DESCSET2( DESCY, 2, 1, 1, 1, 1, 1, 0, 0, ICTXT, 2 )
      INCY = 1
*
      RETURN
*
*     End of PZSETPBLAS
*
      END
      SUBROUTINE PZCHKMAT( ICTXT, NOUT, SUBPTR, SCODE, SNAME, ARGNAM,
     $                     ARGPOS )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1         ARGNAM
      INTEGER             ARGPOS, ICTXT, NOUT, SCODE
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)       SNAME
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL            SUBPTR
*     ..
*
*  Purpose
*  =======
*
*  PZCHKMAT tests the matrix (or vector) ARGNAM in any PBLAS routine.
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
*  SUBPTR  (global input) SUBROUTINE
*          On entry,  SUBPTR  is  a  subroutine. SUBPTR must be declared
*          EXTERNAL in the calling subroutine.
*
*  SCODE   (global input) INTEGER
*          On entry, SCODE specifies the calling sequence code.
*
*  SNAME   (global input) CHARACTER*(*)
*          On entry,  SNAME  specifies  the subroutine name calling this
*          subprogram.
*
*  ARGNAM  (global input) CHARACTER*(*)
*          On entry,  ARGNAM  specifies the name of the matrix or vector
*          to be checked.  ARGNAM can either be 'A', 'B' or 'C' when one
*          wants to check a matrix, and 'X' or 'Y' for a vector.
*
*  ARGPOS  (global input) INTEGER
*          On entry, ARGPOS indicates the position of the first argument
*          of the matrix (or vector) ARGNAM.
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
      INTEGER             DESCMULT
      PARAMETER           ( DESCMULT = 100 )
*     ..
*     .. Local Scalars ..
      INTEGER             I, INFOT, NPROW, NPCOL, MYROW, MYCOL
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PCHKPBE, PZCALLSUB, PZSETPBLAS
*     ..
*     .. External Functions ..
      LOGICAL             LSAME
      EXTERNAL            LSAME
*     ..
*     .. Common Blocks ..
      INTEGER            IA, IB, IC, INCX, INCY, ISCLR, IX, IY, JA, JB,
     $                   JC, JX, JY
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), DESCC( DLEN_ ),
     $                   DESCX( DLEN_ ), DESCY( DLEN_ )
      COMMON             /PBLASD/DESCA, DESCB, DESCC, DESCX, DESCY
      COMMON             /PBLASI/IA, IB, IC, INCX, INCY, ISCLR, IX, IY,
     $                   JA, JB, JC, JX, JY
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( LSAME( ARGNAM, 'A' ) ) THEN
*
*        Check IA. Set all other OK, bad IA
*
         CALL PZSETPBLAS( ICTXT )
         IA    = -1
         INFOT = ARGPOS + 1
         CALL PZCALLSUB( SUBPTR, SCODE )
         CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*        Check JA. Set all other OK, bad JA
*
         CALL PZSETPBLAS( ICTXT )
         JA    = -1
         INFOT = ARGPOS + 2
         CALL PZCALLSUB( SUBPTR, SCODE )
         CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*        Check DESCA. Set all other OK, bad DESCA
*
         DO 10 I = 1, DLEN_
*
*           Set I'th entry of DESCA to incorrect value, rest ok.
*
            CALL PZSETPBLAS( ICTXT )
            DESCA( I ) =  -2
            INFOT = ( ( ARGPOS + 3 ) * DESCMULT ) + I
            CALL PZCALLSUB( SUBPTR, SCODE )
            CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*           Extra tests for RSRCA, CSRCA, LDA
*
            IF( ( I.EQ.RSRC_ ) .OR. ( I.EQ.CSRC_ ) .OR.
     $          ( I.EQ.LLD_ ) ) THEN
*
               CALL PZSETPBLAS( ICTXT )
*
*              Test RSRCA >= NPROW
*
               IF( I.EQ.RSRC_ )
     $            DESCA( I ) =  NPROW
*
*              Test CSRCA >= NPCOL
*
               IF( I.EQ.CSRC_ )
     $            DESCA( I ) =  NPCOL
*
*              Test LDA >= MAX(1, PB_NUMROC(...)). Set to 1 as mat 2x2.
*
               IF( I.EQ.LLD_ ) THEN
                  IF( MYROW.EQ.0 .AND.MYCOL.EQ.0 ) THEN
                     DESCA( I ) = 1
                  ELSE
                     DESCA( I ) = 0
                  END IF
               END IF
*
               INFOT = ( ( ARGPOS + 3 ) * DESCMULT ) + I
               CALL PZCALLSUB( SUBPTR, SCODE )
               CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
            END IF
*
   10    CONTINUE
*
      ELSE IF( LSAME( ARGNAM, 'B' ) ) THEN
*
*        Check IB. Set all other OK, bad IB
*
         CALL PZSETPBLAS( ICTXT )
         IB    = -1
         INFOT = ARGPOS + 1
         CALL PZCALLSUB( SUBPTR, SCODE )
         CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*        Check JB. Set all other OK, bad JB
*
         CALL PZSETPBLAS( ICTXT )
         JB    = -1
         INFOT = ARGPOS + 2
         CALL PZCALLSUB( SUBPTR, SCODE )
         CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*        Check DESCB. Set all other OK, bad DESCB
*
         DO 20 I = 1, DLEN_
*
*           Set I'th entry of DESCB to incorrect value, rest ok.
*
            CALL PZSETPBLAS( ICTXT )
            DESCB( I ) =  -2
            INFOT = ( ( ARGPOS + 3 ) * DESCMULT ) + I
            CALL PZCALLSUB( SUBPTR, SCODE )
            CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*           Extra tests for RSRCB, CSRCB, LDB
*
            IF( ( I.EQ.RSRC_ ) .OR. ( I.EQ.CSRC_ ) .OR.
     $          ( I.EQ.LLD_ ) ) THEN
*
               CALL PZSETPBLAS( ICTXT )
*
*              Test RSRCB >= NPROW
*
               IF( I.EQ.RSRC_ )
     $            DESCB( I ) =  NPROW
*
*              Test CSRCB >= NPCOL
*
               IF( I.EQ.CSRC_ )
     $            DESCB( I ) =  NPCOL
*
*              Test LDB >= MAX(1, PB_NUMROC(...)). Set to 1 as mat 2x2.
*
               IF( I.EQ.LLD_ ) THEN
                  IF( MYROW.EQ.0 .AND.MYCOL.EQ.0 ) THEN
                     DESCB( I ) = 1
                  ELSE
                     DESCB( I ) = 0
                  END IF
               END IF
*
               INFOT = ( ( ARGPOS + 3 ) * DESCMULT ) + I
               CALL PZCALLSUB( SUBPTR, SCODE )
               CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
            END IF
*
   20    CONTINUE
*
      ELSE IF( LSAME( ARGNAM, 'C' ) ) THEN
*
*        Check IC. Set all other OK, bad IC
*
         CALL PZSETPBLAS( ICTXT )
         IC    = -1
         INFOT = ARGPOS + 1
         CALL PZCALLSUB( SUBPTR, SCODE )
         CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*        Check JC. Set all other OK, bad JC
*
         CALL PZSETPBLAS( ICTXT )
         JC    = -1
         INFOT = ARGPOS + 2
         CALL PZCALLSUB( SUBPTR, SCODE )
         CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*        Check DESCC. Set all other OK, bad DESCC
*
         DO 30 I = 1, DLEN_
*
*           Set I'th entry of DESCC to incorrect value, rest ok.
*
            CALL PZSETPBLAS( ICTXT )
            DESCC( I ) =  -2
            INFOT = ( ( ARGPOS + 3 ) * DESCMULT ) + I
            CALL PZCALLSUB( SUBPTR, SCODE )
            CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*           Extra tests for RSRCC, CSRCC, LDC
*
            IF( ( I.EQ.RSRC_ ) .OR. ( I.EQ.CSRC_ ) .OR.
     $          ( I.EQ.LLD_ ) ) THEN
*
               CALL PZSETPBLAS( ICTXT )
*
*              Test RSRCC >= NPROW
*
               IF( I.EQ.RSRC_ )
     $            DESCC( I ) =  NPROW
*
*              Test CSRCC >= NPCOL
*
               IF( I.EQ.CSRC_ )
     $            DESCC( I ) =  NPCOL
*
*              Test LDC >= MAX(1, PB_NUMROC(...)). Set to 1 as mat 2x2.
*
               IF( I.EQ.LLD_ ) THEN
                  IF( MYROW.EQ.0 .AND.MYCOL.EQ.0 ) THEN
                     DESCC( I ) = 1
                  ELSE
                     DESCC( I ) = 0
                  END IF
               END IF
*
               INFOT = ( ( ARGPOS + 3 ) * DESCMULT ) + I
               CALL PZCALLSUB( SUBPTR, SCODE )
               CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
            END IF
*
   30    CONTINUE
*
      ELSE IF( LSAME( ARGNAM, 'X' ) ) THEN
*
*        Check IX. Set all other OK, bad IX
*
         CALL PZSETPBLAS( ICTXT )
         IX    = -1
         INFOT = ARGPOS + 1
         CALL PZCALLSUB( SUBPTR, SCODE )
         CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*        Check JX. Set all other OK, bad JX
*
         CALL PZSETPBLAS( ICTXT )
         JX    = -1
         INFOT = ARGPOS + 2
         CALL PZCALLSUB( SUBPTR, SCODE )
         CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*        Check DESCX. Set all other OK, bad DESCX
*
         DO 40 I = 1, DLEN_
*
*           Set I'th entry of DESCX to incorrect value, rest ok.
*
            CALL PZSETPBLAS( ICTXT )
            DESCX( I ) =  -2
            INFOT = ( ( ARGPOS + 3 ) * DESCMULT ) + I
            CALL PZCALLSUB( SUBPTR, SCODE )
            CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*           Extra tests for RSRCX, CSRCX, LDX
*
            IF( ( I.EQ.RSRC_ ) .OR. ( I.EQ.CSRC_ ) .OR.
     $          ( I.EQ.LLD_ ) ) THEN
*
               CALL PZSETPBLAS( ICTXT )
*
*              Test RSRCX >= NPROW
*
               IF( I.EQ.RSRC_ )
     $            DESCX( I ) =  NPROW
*
*              Test CSRCX >= NPCOL
*
               IF( I.EQ.CSRC_ )
     $            DESCX( I ) =  NPCOL
*
*              Test LDX >= MAX(1, PB_NUMROC(...)). Set to 1 as mat 2x2.
*
               IF( I.EQ.LLD_ ) THEN
                  IF( MYROW.EQ.0 .AND.MYCOL.EQ.0 ) THEN
                     DESCX( I ) = 1
                  ELSE
                     DESCX( I ) = 0
                  END IF
               END IF
*
               INFOT = ( ( ARGPOS + 3 ) * DESCMULT ) + I
               CALL PZCALLSUB( SUBPTR, SCODE )
               CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
            END IF
*
   40    CONTINUE
*
*        Check INCX. Set all other OK, bad INCX
*
         CALL PZSETPBLAS( ICTXT )
         INCX  =  -1
         INFOT = ARGPOS + 4
         CALL PZCALLSUB( SUBPTR, SCODE )
         CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
      ELSE
*
*        Check IY. Set all other OK, bad IY
*
         CALL PZSETPBLAS( ICTXT )
         IY    = -1
         INFOT = ARGPOS + 1
         CALL PZCALLSUB( SUBPTR, SCODE )
         CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*        Check JY. Set all other OK, bad JY
*
         CALL PZSETPBLAS( ICTXT )
         JY    = -1
         INFOT = ARGPOS + 2
         CALL PZCALLSUB( SUBPTR, SCODE )
         CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*        Check DESCY. Set all other OK, bad DESCY
*
         DO 50 I = 1, DLEN_
*
*           Set I'th entry of DESCY to incorrect value, rest ok.
*
            CALL PZSETPBLAS( ICTXT )
            DESCY( I ) =  -2
            INFOT = ( ( ARGPOS + 3 ) * DESCMULT ) + I
            CALL PZCALLSUB( SUBPTR, SCODE )
            CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
*           Extra tests for RSRCY, CSRCY, LDY
*
            IF( ( I.EQ.RSRC_ ) .OR. ( I.EQ.CSRC_ ) .OR.
     $          ( I.EQ.LLD_ ) ) THEN
*
               CALL PZSETPBLAS( ICTXT )
*
*              Test RSRCY >= NPROW
*
               IF( I.EQ.RSRC_ )
     $            DESCY( I ) = NPROW
*
*              Test CSRCY >= NPCOL
*
               IF( I.EQ.CSRC_ )
     $            DESCY( I ) = NPCOL
*
*              Test LDY >= MAX(1, PB_NUMROC(...)). Set to 1 as mat 2x2.
*
               IF( I.EQ.LLD_ ) THEN
                  IF( MYROW.EQ.0 .AND.MYCOL.EQ.0 ) THEN
                     DESCY( I ) = 1
                  ELSE
                     DESCY( I ) = 0
                  END IF
               END IF
*
               INFOT = ( ( ARGPOS + 3 ) * DESCMULT ) + I
               CALL PZCALLSUB( SUBPTR, SCODE )
               CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
            END IF
*
   50    CONTINUE
*
*        Check INCY. Set all other OK, bad INCY
*
         CALL PZSETPBLAS( ICTXT )
         INCY =  -1
         INFOT = ARGPOS + 4
         CALL PZCALLSUB( SUBPTR, SCODE )
         CALL PCHKPBE( ICTXT, NOUT, SNAME, INFOT )
*
      END IF
*
      RETURN
*
*     End of PZCHKMAT
*
      END
      SUBROUTINE PZCALLSUB( SUBPTR, SCODE )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER             SCODE
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL            SUBPTR
*     ..
*
*  Purpose
*  =======
*
*  PZCALLSUB calls the subroutine SUBPTR with the calling sequence iden-
*  tified by SCODE.
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
*  SUBPTR  (global input) SUBROUTINE
*          On entry,  SUBPTR  is  a  subroutine. SUBPTR must be declared
*          EXTERNAL in the calling subroutine.
*
*  SCODE   (global input) INTEGER
*          On entry, SCODE specifies the calling sequence code.
*
*  Calling sequence encodings
*  ==========================
*
*  code Formal argument list                                Examples
*
*  11   (n,      v1,v2)                                     _SWAP, _COPY
*  12   (n,s1,   v1   )                                     _SCAL, _SCAL
*  13   (n,s1,   v1,v2)                                     _AXPY, _DOT_
*  14   (n,s1,i1,v1   )                                     _AMAX
*  15   (n,u1,   v1   )                                     _ASUM, _NRM2
*
*  21   (     trans,     m,n,s1,m1,v1,s2,v2)                _GEMV
*  22   (uplo,             n,s1,m1,v1,s2,v2)                _SYMV, _HEMV
*  23   (uplo,trans,diag,  n,   m1,v1      )                _TRMV, _TRSV
*  24   (                m,n,s1,v1,v2,m1)                   _GER_
*  25   (uplo,             n,s1,v1,   m1)                   _SYR
*  26   (uplo,             n,u1,v1,   m1)                   _HER
*  27   (uplo,             n,s1,v1,v2,m1)                   _SYR2, _HER2
*
*  31   (          transa,transb,     m,n,k,s1,m1,m2,s2,m3) _GEMM
*  32   (side,uplo,                   m,n,  s1,m1,m2,s2,m3) _SYMM, _HEMM
*  33   (     uplo,trans,               n,k,s1,m1,   s2,m3) _SYRK
*  34   (     uplo,trans,               n,k,u1,m1,   u2,m3) _HERK
*  35   (     uplo,trans,               n,k,s1,m1,m2,s2,m3) _SYR2K
*  36   (     uplo,trans,               n,k,s1,m1,m2,u2,m3) _HER2K
*  37   (                             m,n,  s1,m1,   s2,m3) _TRAN_
*  38   (side,uplo,transa,       diag,m,n,  s1,m1,m2      ) _TRMM, _TRSM
*  39   (          trans,             m,n,  s1,m1,   s2,m3) _GEADD
*  40   (     uplo,trans,             m,n,  s1,m1,   s2,m3) _TRADD
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
*     .. Common Blocks ..
      CHARACTER*1        DIAG, SIDE, TRANSA, TRANSB, UPLO
      INTEGER            IA, IB, IC, INCX, INCY, ISCLR, IX, IY, JA, JB,
     $                   JC, JX, JY, KDIM, MDIM, NDIM
      DOUBLE PRECISION   USCLR
      COMPLEX*16         SCLR
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), DESCC( DLEN_ ),
     $                   DESCX( DLEN_ ), DESCY( DLEN_ )
      COMPLEX*16         A( 2, 2 ), B( 2, 2 ), C( 2, 2 ), X( 2 ), Y( 2 )
      COMMON             /PBLASC/DIAG, SIDE, TRANSA, TRANSB, UPLO
      COMMON             /PBLASD/DESCA, DESCB, DESCC, DESCX, DESCY
      COMMON             /PBLASI/IA, IB, IC, INCX, INCY, ISCLR, IX, IY,
     $                   JA, JB, JC, JX, JY
      COMMON             /PBLASM/A, B, C
      COMMON             /PBLASN/KDIM, MDIM, NDIM
      COMMON             /PBLASS/SCLR, USCLR
      COMMON             /PBLASV/X, Y
*     ..
*     .. Executable Statements ..
*
*     Level 1 PBLAS
*
      IF( SCODE.EQ.11 ) THEN
*
         CALL SUBPTR( NDIM, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY,
     $                INCY )
*
      ELSE IF( SCODE.EQ.12 ) THEN
*
         CALL SUBPTR( NDIM, SCLR, X, IX, JX, DESCX, INCX )
*
      ELSE IF( SCODE.EQ.13 ) THEN
*
         CALL SUBPTR( NDIM, SCLR, X, IX, JX, DESCX, INCX, Y, IY, JY,
     $                DESCY, INCY )
*
      ELSE IF( SCODE.EQ.14 ) THEN
*
         CALL SUBPTR( NDIM, SCLR, ISCLR, X, IX, JX, DESCX, INCX )
*
      ELSE IF( SCODE.EQ.15 ) THEN
*
         CALL SUBPTR( NDIM, USCLR, X, IX, JX, DESCX, INCX )
*
*     Level 2 PBLAS
*
      ELSE IF( SCODE.EQ.21 ) THEN
*
         CALL SUBPTR( TRANSA, MDIM, NDIM, SCLR, A, IA, JA, DESCA, X, IX,
     $                JX, DESCX, INCX, SCLR, Y, IY, JY, DESCY, INCY )
*
      ELSE IF( SCODE.EQ.22 ) THEN
*
         CALL SUBPTR( UPLO, NDIM, SCLR, A, IA, JA, DESCA, X, IX, JX,
     $                DESCX, INCX, SCLR, Y, IY, JY, DESCY, INCY )
*
      ELSE IF( SCODE.EQ.23 ) THEN
*
         CALL SUBPTR( UPLO, TRANSA, DIAG, NDIM, A, IA, JA, DESCA, X, IX,
     $                JX, DESCX, INCX )
*
      ELSE IF( SCODE.EQ.24 ) THEN
*
         CALL SUBPTR( MDIM, NDIM, SCLR, X, IX, JX, DESCX, INCX, Y, IY,
     $                JY, DESCY, INCY, A, IA, JA, DESCA )
*
      ELSE IF( SCODE.EQ.25 ) THEN
*
         CALL SUBPTR( UPLO, NDIM, SCLR, X, IX, JX, DESCX, INCX, A, IA,
     $                JA, DESCA )
*
      ELSE IF( SCODE.EQ.26 ) THEN
*
         CALL SUBPTR( UPLO, NDIM, USCLR, X, IX, JX, DESCX, INCX, A, IA,
     $                JA, DESCA )
*
      ELSE IF( SCODE.EQ.27 ) THEN
*
         CALL SUBPTR( UPLO, NDIM, SCLR, X, IX, JX, DESCX, INCX, Y, IY,
     $                JY, DESCY, INCY, A, IA, JA, DESCA )
*
*     Level 3 PBLAS
*
      ELSE IF( SCODE.EQ.31 ) THEN
*
         CALL SUBPTR( TRANSA, TRANSB, MDIM, NDIM, KDIM, SCLR, A, IA, JA,
     $                DESCA, B, IB, JB, DESCB, SCLR, C, IC, JC, DESCC )
*
      ELSE IF( SCODE.EQ.32 ) THEN
*
         CALL SUBPTR( SIDE, UPLO, MDIM, NDIM, SCLR, A, IA, JA, DESCA, B,
     $                IB, JB, DESCB, SCLR, C, IC, JC, DESCC )
*
      ELSE IF( SCODE.EQ.33 ) THEN
*
         CALL SUBPTR( UPLO, TRANSA, NDIM, KDIM, SCLR, A, IA, JA, DESCA,
     $                SCLR, C, IC, JC, DESCC )
*
      ELSE IF( SCODE.EQ.34 ) THEN
*
         CALL SUBPTR( UPLO, TRANSA, NDIM, KDIM, USCLR, A, IA, JA, DESCA,
     $                USCLR, C, IC, JC, DESCC )
*
      ELSE IF( SCODE.EQ.35 ) THEN
*
         CALL SUBPTR( UPLO, TRANSA, NDIM, KDIM, SCLR, A, IA, JA, DESCA,
     $                B, IB, JB, DESCB, SCLR, C, IC, JC, DESCC )
*
      ELSE IF( SCODE.EQ.36 ) THEN
*
         CALL SUBPTR( UPLO, TRANSA, NDIM, KDIM, SCLR, A, IA, JA, DESCA,
     $                B, IB, JB, DESCB, USCLR, C, IC, JC, DESCC )
*
      ELSE IF( SCODE.EQ.37 ) THEN
*
         CALL SUBPTR( MDIM, NDIM, SCLR, A, IA, JA, DESCA, SCLR, C, IC,
     $                JC, DESCC )
*
      ELSE IF( SCODE.EQ.38 ) THEN
*
         CALL SUBPTR( SIDE, UPLO, TRANSA, DIAG, MDIM, NDIM, SCLR, A, IA,
     $                JA, DESCA, B, IB, JB, DESCB )
*
      ELSE IF( SCODE.EQ.39 ) THEN
*
         CALL SUBPTR( TRANSA, MDIM, NDIM, SCLR, A, IA, JA, DESCA, SCLR,
     $                C, IC, JC, DESCC )
*
      ELSE IF( SCODE.EQ.40 ) THEN
*
         CALL SUBPTR( UPLO, TRANSA, MDIM, NDIM, SCLR, A, IA, JA, DESCA,
     $                SCLR, C, IC, JC, DESCC )
*
      END IF
*
      RETURN
*
*     End of PZCALLSUB
*
      END
      SUBROUTINE PZERRSET( ERR, ERRMAX, XTRUE, X )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ERR, ERRMAX
      COMPLEX*16         X, XTRUE
*     ..
*
*  Purpose
*  =======
*
*  PZERRSET  computes the absolute difference ERR = |XTRUE - X| and com-
*  pares it with zero. ERRMAX accumulates the absolute error difference.
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
*  ERR     (local output) DOUBLE PRECISION
*          On exit, ERR specifies the absolute difference |XTRUE - X|.
*
*  ERRMAX  (local input/local output) DOUBLE PRECISION
*          On entry,  ERRMAX  specifies  a previously computed error. On
*          exit ERRMAX is the accumulated error MAX( ERRMAX, ERR ).
*
*  XTRUE   (local input) COMPLEX*16
*          On entry, XTRUE specifies the true value.
*
*  X       (local input) COMPLEX*16
*          On entry, X specifies the value to be compared to XTRUE.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. External Functions ..
      DOUBLE PRECISION   PDDIFF
      EXTERNAL           PDDIFF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX
*     ..
*     .. Executable Statements ..
*
      ERR = ABS( PDDIFF( DBLE( XTRUE ), DBLE( X ) ) )
      ERR = MAX( ERR, ABS( PDDIFF( DIMAG( XTRUE ), DIMAG( X ) ) ) )
*
      ERRMAX = MAX( ERRMAX, ERR )
*
      RETURN
*
*     End of PZERRSET
*
      END
      SUBROUTINE PZCHKVIN( ERRMAX, N, X, PX, IX, JX, DESCX, INCX,
     $                     INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INFO, IX, JX, N
      DOUBLE PRECISION   ERRMAX
*     ..
*     .. Array Arguments ..
      INTEGER            DESCX( * )
      COMPLEX*16         PX( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  PZCHKVIN  checks that the submatrix sub( PX ) remained unchanged. The
*  local  array  entries are compared element by element, and their dif-
*  ference  is tested against 0.0 as well as the epsilon machine. Notice
*  that  this difference should be numerically exactly the zero machine,
*  but  because of the possible fluctuation of some of the data we flag-
*  ged differently a difference less than twice the epsilon machine. The
*  largest error is also returned.
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
*  ERRMAX  (global output) DOUBLE PRECISION
*          On exit,  ERRMAX  specifies the largest absolute element-wise
*          difference between sub( X ) and sub( PX ).
*
*  N       (global input) INTEGER
*          On entry,  N  specifies  the  length of the subvector operand
*          sub( X ). N must be at least zero.
*
*  X       (local input) COMPLEX*16 array
*          On entry, X is an array of  dimension  (DESCX( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PX.
*
*  PX      (local input) COMPLEX*16 array
*          On entry, PX is an array of dimension (DESCX( LLD_ ),*). This
*          array contains the local entries of the matrix PX.
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
*          On exit, if INFO = 0, no error has been found,
*          If INFO > 0, the maximum abolute error found is in (0,eps],
*          If INFO < 0, the maximum abolute error found is in (eps,+oo).
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
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, ROWREP
      INTEGER            I, IB, ICTXT, ICURCOL, ICURROW, IIX, IN, IXCOL,
     $                   IXROW, J, JB, JJX, JN, KK, LDPX, LDX, LL,
     $                   MYCOL, MYROW, NPCOL, NPROW
      DOUBLE PRECISION   ERR, EPS
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGAMX2D, PB_INFOG2L, PZERRSET
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      ERRMAX = ZERO
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      ICTXT = DESCX( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      EPS = PDLAMCH( ICTXT, 'eps' )
*
      CALL PB_INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL, IIX,
     $                 JJX, IXROW, IXCOL )
*
      LDX    = DESCX( M_ )
      LDPX   = DESCX( LLD_ )
      ROWREP = ( IXROW.EQ.-1 )
      COLREP = ( IXCOL.EQ.-1 )
*
      IF( N.EQ.1 ) THEN
*
         IF( ( MYROW.EQ.IXROW .OR. ROWREP ) .AND.
     $       ( MYCOL.EQ.IXCOL .OR. COLREP ) )
     $      CALL PZERRSET( ERR, ERRMAX, X( IX+(JX-1)*LDX ),
     $                     PX( IIX+(JJX-1)*LDPX ) )
*
      ELSE IF( INCX.EQ.DESCX( M_ ) ) THEN
*
*        sub( X ) is a row vector
*
         JB = DESCX( INB_ ) - JX + 1
         IF( JB.LE.0 )
     $      JB = ( ( -JB ) / DESCX( NB_ ) + 1 ) * DESCX( NB_ ) + JB
         JB = MIN( JB, N )
         JN = JX + JB - 1
*
         IF( MYROW.EQ.IXROW .OR. ROWREP ) THEN
*
            ICURCOL = IXCOL
            IF( MYCOL.EQ.ICURCOL .OR. COLREP ) THEN
               DO 10 J = JX, JN
                  CALL PZERRSET( ERR, ERRMAX, X( IX+(J-1)*LDX ),
     $                           PX( IIX+(JJX-1)*LDPX ) )
                  JJX = JJX + 1
   10          CONTINUE
            END IF
            ICURCOL = MOD( ICURCOL+1, NPCOL )
*
            DO 30 J = JN+1, JX+N-1, DESCX( NB_ )
               JB = MIN( JX+N-J, DESCX( NB_ ) )
*
               IF( MYCOL.EQ.ICURCOL .OR. COLREP ) THEN
*
                  DO 20 KK = 0, JB-1
                     CALL PZERRSET( ERR, ERRMAX, X( IX+(J+KK-1)*LDX ),
     $                              PX( IIX+(JJX+KK-1)*LDPX ) )
   20             CONTINUE
*
                  JJX = JJX + JB
*
               END IF
*
               ICURCOL = MOD( ICURCOL+1, NPCOL )
*
   30       CONTINUE
*
         END IF
*
      ELSE
*
*        sub( X ) is a column vector
*
         IB = DESCX( IMB_ ) - IX + 1
         IF( IB.LE.0 )
     $      IB = ( ( -IB ) / DESCX( MB_ ) + 1 ) * DESCX( MB_ ) + IB
         IB = MIN( IB, N )
         IN = IX + IB - 1
*
         IF( MYCOL.EQ.IXCOL .OR. COLREP ) THEN
*
            ICURROW = IXROW
            IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
               DO 40 I = IX, IN
                  CALL PZERRSET( ERR, ERRMAX, X( I+(JX-1)*LDX ),
     $                           PX( IIX+(JJX-1)*LDPX ) )
                  IIX = IIX + 1
   40          CONTINUE
            END IF
            ICURROW = MOD( ICURROW+1, NPROW )
*
            DO 60 I = IN+1, IX+N-1, DESCX( MB_ )
               IB = MIN( IX+N-I, DESCX( MB_ ) )
*
               IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
*
                  DO 50 KK = 0, IB-1
                     CALL PZERRSET( ERR, ERRMAX, X( I+KK+(JX-1)*LDX ),
     $                              PX( IIX+KK+(JJX-1)*LDPX ) )
   50             CONTINUE
*
                  IIX = IIX + IB
*
               END IF
*
               ICURROW = MOD( ICURROW+1, NPROW )
*
   60       CONTINUE
*
         END IF
*
      END IF
*
      CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERRMAX, 1, KK, LL, -1,
     $              -1, -1 )
*
      IF( ERRMAX.GT.ZERO .AND. ERRMAX.LE.EPS ) THEN
         INFO = 1
      ELSE IF( ERRMAX.GT.EPS ) THEN
         INFO = -1
      END IF
*
      RETURN
*
*     End of PZCHKVIN
*
      END
      SUBROUTINE PZCHKVOUT( N, X, PX, IX, JX, DESCX, INCX, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INFO, IX, JX, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCX( * )
      COMPLEX*16         PX( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  PZCHKVOUT  checks  that the matrix PX \ sub( PX ) remained unchanged.
*  The  local array  entries  are compared element by element, and their
*  difference  is tested against 0.0 as well as the epsilon machine. No-
*  tice that this  difference should be numerically exactly the zero ma-
*  chine, but because  of  the  possible movement of some of the data we
*  flagged differently a difference less than twice the epsilon machine.
*  The largest error is reported.
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
*  N       (global input) INTEGER
*          On entry,  N  specifies  the  length of the subvector operand
*          sub( X ). N must be at least zero.
*
*  X       (local input) COMPLEX*16 array
*          On entry, X is an array of  dimension  (DESCX( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PX.
*
*  PX      (local input) COMPLEX*16 array
*          On entry, PX is an array of dimension (DESCX( LLD_ ),*). This
*          array contains the local entries of the matrix PX.
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
*          On exit, if INFO = 0, no error has been found,
*          If INFO > 0, the maximum abolute error found is in (0,eps],
*          If INFO < 0, the maximum abolute error found is in (eps,+oo).
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
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, ROWREP
      INTEGER            I, IB, ICTXT, ICURCOL, ICURROW, II, IMBX, INBX,
     $                   J, JB, JJ, KK, LDPX, LDX, LL, MBX, MPALL,
     $                   MYCOL, MYCOLDIST, MYROW, MYROWDIST, NBX, NPCOL,
     $                   NPROW, NQALL
      DOUBLE PRECISION   EPS, ERR, ERRMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGAMX2D, PZERRSET
*     ..
*     .. External Functions ..
      INTEGER            PB_NUMROC
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           PDLAMCH, PB_NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      ERRMAX = ZERO
*
*     Quick return if possible
*
      IF( ( DESCX( M_ ).LE.0 ).OR.( DESCX( N_ ).LE.0 ) )
     $   RETURN
*
*     Start the operations
*
      ICTXT = DESCX( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      EPS = PDLAMCH( ICTXT, 'eps' )
*
      MPALL   = PB_NUMROC( DESCX( M_ ), 1, DESCX( IMB_ ), DESCX( MB_ ),
     $                     MYROW, DESCX( RSRC_ ), NPROW )
      NQALL   = PB_NUMROC( DESCX( N_ ), 1, DESCX( INB_ ), DESCX( NB_ ),
     $                     MYCOL, DESCX( CSRC_ ), NPCOL )
*
      MBX     = DESCX( MB_ )
      NBX     = DESCX( NB_ )
      LDX     = DESCX( M_ )
      LDPX    = DESCX( LLD_ )
      ICURROW = DESCX( RSRC_ )
      ICURCOL = DESCX( CSRC_ )
      ROWREP  = ( ICURROW.EQ.-1 )
      COLREP  = ( ICURCOL.EQ.-1 )
      IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
         IMBX = DESCX( IMB_ )
      ELSE
         IMBX = MBX
      END IF
      IF( MYCOL.EQ.ICURCOL .OR. COLREP ) THEN
         INBX = DESCX( INB_ )
      ELSE
         INBX = NBX
      END IF
      IF( ROWREP ) THEN
         MYROWDIST = 0
      ELSE
         MYROWDIST = MOD( MYROW - ICURROW + NPROW, NPROW )
      END IF
      IF( COLREP ) THEN
         MYCOLDIST = 0
      ELSE
         MYCOLDIST = MOD( MYCOL - ICURCOL + NPCOL, NPCOL )
      END IF
      II = 1
      JJ = 1
*
      IF( INCX.EQ.DESCX( M_ ) ) THEN
*
*        sub( X ) is a row vector
*
         IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
*
            I = 1
            IF( MYCOLDIST.EQ.0 ) THEN
               J = 1
            ELSE
               J = DESCX( INB_ ) + ( MYCOLDIST - 1 ) * NBX + 1
            END IF
            JB = MIN( MAX( 0, DESCX( N_ ) - J + 1 ), INBX )
            IB = MIN( DESCX( M_ ), DESCX( IMB_ ) )
*
            DO 20 KK = 0, JB-1
               DO 10 LL = 0, IB-1
                  IF( I+LL.NE.IX .OR. J+KK.LT.JX .OR. J+KK.GT.JX+N-1 )
     $               CALL PZERRSET( ERR, ERRMAX,
     $                              X( I+LL+(J+KK-1)*LDX ),
     $                              PX( II+LL+(JJ+KK-1)*LDPX ) )
   10          CONTINUE
   20       CONTINUE
            IF( COLREP ) THEN
               J = J + INBX
            ELSE
               J = J + INBX + ( NPCOL - 1 ) * NBX
            END IF
*
            DO 50 JJ = INBX+1, NQALL, NBX
               JB = MIN( NQALL-JJ+1, NBX )
*
               DO 40 KK = 0, JB-1
                  DO 30 LL = 0, IB-1
                     IF( I+LL.NE.IX .OR. J+KK.LT.JX .OR.
     $                   J+KK.GT.JX+N-1 )
     $                  CALL PZERRSET( ERR, ERRMAX,
     $                                 X( I+LL+(J+KK-1)*LDX ),
     $                                 PX( II+LL+(JJ+KK-1)*LDPX ) )
   30             CONTINUE
   40          CONTINUE
*
               IF( COLREP ) THEN
                  J = J + NBX
               ELSE
                  J = J + NPCOL * NBX
               END IF
*
   50       CONTINUE
*
            II = II + IB
*
         END IF
*
         ICURROW = MOD( ICURROW + 1, NPROW )
*
         DO 110 I = DESCX( IMB_ ) + 1, DESCX( M_ ), MBX
            IB = MIN( DESCX( M_ ) - I + 1, MBX )
*
            IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
*
               IF( MYCOLDIST.EQ.0 ) THEN
                  J = 1
               ELSE
                  J = DESCX( INB_ ) + ( MYCOLDIST - 1 ) * NBX + 1
               END IF
*
               JJ = 1
               JB = MIN( MAX( 0, DESCX( N_ ) - J + 1 ), INBX )
               DO 70 KK = 0, JB-1
                  DO 60 LL = 0, IB-1
                     IF( I+LL.NE.IX .OR. J+KK.LT.JX .OR.
     $                   J+KK.GT.JX+N-1 )
     $                  CALL PZERRSET( ERR, ERRMAX,
     $                                 X( I+LL+(J+KK-1)*LDX ),
     $                                 PX( II+LL+(JJ+KK-1)*LDPX ) )
   60             CONTINUE
   70          CONTINUE
               IF( COLREP ) THEN
                  J = J + INBX
               ELSE
                  J = J + INBX + ( NPCOL - 1 ) * NBX
               END IF
*
               DO 100 JJ = INBX+1, NQALL, NBX
                  JB = MIN( NQALL-JJ+1, NBX )
*
                  DO 90 KK = 0, JB-1
                     DO 80 LL = 0, IB-1
                        IF( I+LL.NE.IX .OR. J+KK.LT.JX .OR.
     $                      J+KK.GT.JX+N-1 )
     $                     CALL PZERRSET( ERR, ERRMAX,
     $                                    X( I+LL+(J+KK-1)*LDX ),
     $                                    PX( II+LL+(JJ+KK-1)*LDPX ) )
   80                CONTINUE
   90             CONTINUE
*
                  IF( COLREP ) THEN
                     J = J + NBX
                  ELSE
                     J = J + NPCOL * NBX
                  END IF
*
  100          CONTINUE
*
               II = II + IB
*
            END IF
*
            ICURROW = MOD( ICURROW + 1, NPROW )
*
  110    CONTINUE
*
      ELSE
*
*        sub( X ) is a column vector
*
         IF( MYCOL.EQ.ICURCOL .OR. COLREP ) THEN
*
            J = 1
            IF( MYROWDIST.EQ.0 ) THEN
               I = 1
            ELSE
               I = DESCX( IMB_ ) + ( MYROWDIST - 1 ) * MBX + 1
            END IF
            IB = MIN( MAX( 0, DESCX( M_ ) - I + 1 ), IMBX )
            JB = MIN( DESCX( N_ ), DESCX( INB_ ) )
*
            DO 130 KK = 0, JB-1
               DO 120 LL = 0, IB-1
                  IF( J+KK.NE.JX .OR. I+LL.LT.IX .OR. I+LL.GT.IX+N-1 )
     $               CALL PZERRSET( ERR, ERRMAX,
     $                              X( I+LL+(J+KK-1)*LDX ),
     $                              PX( II+LL+(JJ+KK-1)*LDPX ) )
  120          CONTINUE
  130       CONTINUE
            IF( ROWREP ) THEN
               I = I + IMBX
            ELSE
               I = I + IMBX + ( NPROW - 1 ) * MBX
            END IF
*
            DO 160 II = IMBX+1, MPALL, MBX
               IB = MIN( MPALL-II+1, MBX )
*
               DO 150 KK = 0, JB-1
                  DO 140 LL = 0, IB-1
                     IF( J+KK.NE.JX .OR. I+LL.LT.IX .OR.
     $                   I+LL.GT.IX+N-1 )
     $                  CALL PZERRSET( ERR, ERRMAX,
     $                                 X( I+LL+(J+KK-1)*LDX ),
     $                                 PX( II+LL+(JJ+KK-1)*LDPX ) )
  140             CONTINUE
  150          CONTINUE
*
               IF( ROWREP ) THEN
                  I = I + MBX
               ELSE
                  I = I + NPROW * MBX
               END IF
*
  160       CONTINUE
*
            JJ = JJ + JB
*
         END IF
*
         ICURCOL = MOD( ICURCOL + 1, NPCOL )
*
         DO 220 J = DESCX( INB_ ) + 1, DESCX( N_ ), NBX
            JB = MIN( DESCX( N_ ) - J + 1, NBX )
*
            IF( MYCOL.EQ.ICURCOL .OR. COLREP ) THEN
*
               IF( MYROWDIST.EQ.0 ) THEN
                  I = 1
               ELSE
                  I = DESCX( IMB_ ) + ( MYROWDIST - 1 ) * MBX + 1
               END IF
*
               II = 1
               IB = MIN( MAX( 0, DESCX( M_ ) - I + 1 ), IMBX )
               DO 180 KK = 0, JB-1
                  DO 170 LL = 0, IB-1
                     IF( J+KK.NE.JX .OR. I+LL.LT.IX .OR.
     $                   I+LL.GT.IX+N-1 )
     $                  CALL PZERRSET( ERR, ERRMAX,
     $                                 X( I+LL+(J+KK-1)*LDX ),
     $                                 PX( II+LL+(JJ+KK-1)*LDPX ) )
  170             CONTINUE
  180          CONTINUE
               IF( ROWREP ) THEN
                  I = I + IMBX
               ELSE
                  I = I + IMBX + ( NPROW - 1 ) * MBX
               END IF
*
               DO 210 II = IMBX+1, MPALL, MBX
                  IB = MIN( MPALL-II+1, MBX )
*
                  DO 200 KK = 0, JB-1
                     DO 190 LL = 0, IB-1
                        IF( J+KK.NE.JX .OR. I+LL.LT.IX .OR.
     $                      I+LL.GT.IX+N-1 )
     $                     CALL PZERRSET( ERR, ERRMAX,
     $                                    X( I+LL+(J+KK-1)*LDX ),
     $                                    PX( II+LL+(JJ+KK-1)*LDPX ) )
  190                CONTINUE
  200             CONTINUE
*
                  IF( ROWREP ) THEN
                     I = I + MBX
                  ELSE
                     I = I + NPROW * MBX
                  END IF
*
  210          CONTINUE
*
               JJ = JJ + JB
*
            END IF
*
            ICURCOL = MOD( ICURCOL + 1, NPCOL )
*
  220    CONTINUE
*
      END IF
*
      CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERRMAX, 1, KK, LL, -1,
     $              -1, -1 )
*
      IF( ERRMAX.GT.ZERO .AND. ERRMAX.LE.EPS ) THEN
         INFO = 1
      ELSE IF( ERRMAX.GT.EPS ) THEN
         INFO = -1
      END IF
*
      RETURN
*
*     End of PZCHKVOUT
*
      END
      SUBROUTINE PZCHKMIN( ERRMAX, M, N, A, PA, IA, JA, DESCA, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, M, N
      DOUBLE PRECISION   ERRMAX
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         PA( * ), A( * )
*     ..
*
*  Purpose
*  =======
*
*  PZCHKMIN  checks that the submatrix sub( PA ) remained unchanged. The
*  local  array  entries are compared element by element, and their dif-
*  ference  is tested against 0.0 as well as the epsilon machine. Notice
*  that  this difference should be numerically exactly the zero machine,
*  but  because of the possible fluctuation of some of the data we flag-
*  ged differently a difference less than twice the epsilon machine. The
*  largest error is also returned.
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
*  ERRMAX  (global output) DOUBLE PRECISION
*          On exit,  ERRMAX  specifies the largest absolute element-wise
*          difference between sub( A ) and sub( PA ).
*
*  M       (global input) INTEGER
*          On entry,  M  specifies  the  number of rows of the submatrix
*          operand sub( A ). M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          operand sub( A ). N must be at least zero.
*
*  A       (local input) COMPLEX*16 array
*          On entry, A is an array of  dimension  (DESCA( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PA.
*
*  PA      (local input) COMPLEX*16 array
*          On entry, PA is an array of dimension (DESCA( LLD_ ),*). This
*          array contains the local entries of the matrix PA.
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
*          On exit, if INFO = 0, no error has been found,
*          If INFO > 0, the maximum abolute error found is in (0,eps],
*          If INFO < 0, the maximum abolute error found is in (eps,+oo).
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
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, ROWREP
      INTEGER            H, I, IACOL, IAROW, IB, ICTXT, ICURCOL,
     $                   ICURROW, II, IIA, IN, J, JB, JJ, JJA, JN, K,
     $                   KK, LDA, LDPA, LL, MYCOL, MYROW, NPCOL, NPROW
      DOUBLE PRECISION   ERR, EPS
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGAMX2D, PB_INFOG2L, PZERRSET
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
      INFO   = 0
      ERRMAX = ZERO
*
*     Quick return if posssible
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ) )
     $   RETURN
*
*     Start the operations
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      EPS = PDLAMCH( ICTXT, 'eps' )
*
      CALL PB_INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA,
     $                 JJA, IAROW, IACOL )
*
      II      = IIA
      JJ      = JJA
      LDA     = DESCA( M_ )
      LDPA    = DESCA( LLD_ )
      ICURROW = IAROW
      ICURCOL = IACOL
      ROWREP  = ( IAROW.EQ.-1 )
      COLREP  = ( IACOL.EQ.-1 )
*
*     Handle the first block of column separately
*
      JB = DESCA( INB_ ) - JA  + 1
      IF( JB.LE.0 )
     $   JB = ( ( -JB ) / DESCA( NB_ ) + 1 ) * DESCA( NB_ ) + JB
      JB = MIN( JB, N )
      JN = JA + JB - 1
*
      IF( MYCOL.EQ.ICURCOL .OR. COLREP ) THEN
*
         DO 40 H = 0, JB-1
            IB = DESCA( IMB_ ) - IA  + 1
            IF( IB.LE.0 )
     $         IB = ( ( -IB ) / DESCA( MB_ ) + 1 ) * DESCA( MB_ ) + IB
            IB = MIN( IB, M )
            IN = IA + IB - 1
            IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
               DO 10 K = 0, IB-1
                  CALL PZERRSET( ERR, ERRMAX, A( IA+K+(JA+H-1)*LDA ),
     $                           PA( II+K+(JJ+H-1)*LDPA ) )
   10          CONTINUE
               II = II + IB
            END IF
            ICURROW = MOD( ICURROW+1, NPROW )
*
*           Loop over remaining block of rows
*
            DO 30 I = IN+1, IA+M-1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IA+M-I )
               IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                  DO 20 K = 0, IB-1
                     CALL PZERRSET( ERR, ERRMAX, A( I+K+(JA+H-1)*LDA ),
     $                              PA( II+K+(JJ+H-1)*LDPA ) )
   20             CONTINUE
                  II = II + IB
               END IF
               ICURROW = MOD( ICURROW+1, NPROW )
   30       CONTINUE
*
            II = IIA
            ICURROW = IAROW
   40    CONTINUE
*
         JJ = JJ + JB
*
      END IF
*
      ICURCOL = MOD( ICURCOL+1, NPCOL )
*
*     Loop over remaining column blocks
*
      DO 90 J = JN+1, JA+N-1, DESCA( NB_ )
         JB = MIN(  DESCA( NB_ ), JA+N-J )
         IF( MYCOL.EQ.ICURCOL .OR. COLREP ) THEN
            DO 80 H = 0, JB-1
               IB = DESCA( IMB_ ) - IA  + 1
               IF( IB.LE.0 )
     $            IB = ( ( -IB ) / DESCA( MB_ ) + 1 )*DESCA( MB_ ) + IB
               IB = MIN( IB, M )
               IN = IA + IB - 1
               IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                  DO 50 K = 0, IB-1
                     CALL PZERRSET( ERR, ERRMAX, A( IA+K+(J+H-1)*LDA ),
     $                              PA( II+K+(JJ+H-1)*LDPA ) )
   50             CONTINUE
                  II = II + IB
               END IF
               ICURROW = MOD( ICURROW+1, NPROW )
*
*              Loop over remaining block of rows
*
               DO 70 I = IN+1, IA+M-1, DESCA( MB_ )
                  IB = MIN( DESCA( MB_ ), IA+M-I )
                  IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                     DO 60 K = 0, IB-1
                        CALL PZERRSET( ERR, ERRMAX,
     $                                 A( I+K+(J+H-1)*LDA ),
     $                                 PA( II+K+(JJ+H-1)*LDPA ) )
   60                CONTINUE
                     II = II + IB
                  END IF
                  ICURROW = MOD( ICURROW+1, NPROW )
   70          CONTINUE
*
               II = IIA
               ICURROW = IAROW
   80       CONTINUE
*
            JJ = JJ + JB
         END IF
*
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
   90 CONTINUE
*
      CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERRMAX, 1, KK, LL, -1,
     $              -1, -1 )
*
      IF( ERRMAX.GT.ZERO .AND. ERRMAX.LE.EPS ) THEN
         INFO = 1
      ELSE IF( ERRMAX.GT.EPS ) THEN
         INFO = -1
      END IF
*
      RETURN
*
*     End of PZCHKMIN
*
      END
      SUBROUTINE PZCHKMOUT( M, N, A, PA, IA, JA, DESCA, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         A( * ), PA( * )
*     ..
*
*  Purpose
*  =======
*
*  PZCHKMOUT  checks  that the matrix PA \ sub( PA ) remained unchanged.
*  The  local array  entries  are compared element by element, and their
*  difference  is tested against 0.0 as well as the epsilon machine. No-
*  tice that this  difference should be numerically exactly the zero ma-
*  chine, but because  of  the  possible movement of some of the data we
*  flagged differently a difference less than twice the epsilon machine.
*  The largest error is reported.
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
*          On entry,  M  specifies  the  number of rows of the submatrix
*          sub( PA ). M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N specifies the  number of columns of the submatrix
*          sub( PA ). N must be at least zero.
*
*  A       (local input) COMPLEX*16 array
*          On entry, A is an array of  dimension  (DESCA( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PA.
*
*  PA      (local input) COMPLEX*16 array
*          On entry, PA is an array of dimension (DESCA( LLD_ ),*). This
*          array contains the local entries of the matrix PA.
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
*          On exit, if INFO = 0, no error has been found,
*          If INFO > 0, the maximum abolute error found is in (0,eps],
*          If INFO < 0, the maximum abolute error found is in (eps,+oo).
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
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, ROWREP
      INTEGER            I, IB, ICTXT, ICURCOL, II, IMBA, J, JB, JJ, KK,
     $                   LDA, LDPA, LL, MPALL, MYCOL, MYROW, MYROWDIST,
     $                   NPCOL, NPROW
      DOUBLE PRECISION   EPS, ERR, ERRMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGAMX2D, PZERRSET
*     ..
*     .. External Functions ..
      INTEGER            PB_NUMROC
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           PDLAMCH, PB_NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      ERRMAX = ZERO
*
*     Quick return if possible
*
      IF( ( DESCA( M_ ).LE.0 ).OR.( DESCA( N_ ).LE.0 ) )
     $   RETURN
*
*     Start the operations
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      EPS = PDLAMCH( ICTXT, 'eps' )
*
      MPALL = PB_NUMROC( DESCA( M_ ), 1, DESCA( IMB_ ), DESCA( MB_ ),
     $                   MYROW, DESCA( RSRC_ ), NPROW )
*
      LDA    = DESCA( M_ )
      LDPA   = DESCA( LLD_ )
*
      II = 1
      JJ = 1
      ROWREP  = ( DESCA( RSRC_ ).EQ.-1 )
      COLREP  = ( DESCA( CSRC_ ).EQ.-1 )
      ICURCOL = DESCA( CSRC_ )
      IF( MYROW.EQ.DESCA( RSRC_ ) .OR. ROWREP ) THEN
         IMBA = DESCA( IMB_ )
      ELSE
         IMBA = DESCA( MB_ )
      END IF
      IF( ROWREP ) THEN
         MYROWDIST = 0
      ELSE
         MYROWDIST = MOD( MYROW - DESCA( RSRC_ ) + NPROW, NPROW )
      END IF
*
      IF( MYCOL.EQ.ICURCOL .OR. COLREP ) THEN
*
         J = 1
         IF( MYROWDIST.EQ.0 ) THEN
            I = 1
         ELSE
            I = DESCA( IMB_ ) + ( MYROWDIST - 1 ) * DESCA( MB_ ) + 1
         END IF
         IB = MIN( MAX( 0, DESCA( M_ ) - I + 1 ), IMBA )
         JB = MIN( DESCA( N_ ), DESCA( INB_ ) )
*
         DO 20 KK = 0, JB-1
            DO 10 LL = 0, IB-1
               IF( I+LL.LT.IA .OR. I+LL.GT.IA+M-1 .OR.
     $             J+KK.LT.JA .OR. J+KK.GT.JA+N-1 )
     $            CALL PZERRSET( ERR, ERRMAX, A( I+LL+(J+KK-1)*LDA ),
     $                           PA( II+LL+(JJ+KK-1)*LDPA ) )
   10       CONTINUE
   20    CONTINUE
         IF( ROWREP ) THEN
            I = I + IMBA
         ELSE
            I = I + IMBA + ( NPROW - 1 ) * DESCA( MB_ )
         END IF
*
         DO 50 II = IMBA + 1, MPALL, DESCA( MB_ )
            IB = MIN( MPALL-II+1, DESCA( MB_ ) )
*
            DO 40 KK = 0, JB-1
               DO 30 LL = 0, IB-1
                  IF( I+LL.LT.IA .OR. I+LL.GT.IA+M-1 .OR.
     $                J+KK.LT.JA .OR. J+KK.GT.JA+N-1 )
     $               CALL PZERRSET( ERR, ERRMAX,
     $                              A( I+LL+(J+KK-1)*LDA ),
     $                              PA( II+LL+(JJ+KK-1)*LDPA ) )
   30          CONTINUE
   40       CONTINUE
*
            IF( ROWREP ) THEN
               I = I + DESCA( MB_ )
            ELSE
               I = I + NPROW * DESCA( MB_ )
            END IF
*
   50    CONTINUE
*
         JJ = JJ + JB
*
      END IF
*
      ICURCOL = MOD( ICURCOL + 1, NPCOL )
*
      DO 110 J = DESCA( INB_ ) + 1, DESCA( N_ ), DESCA( NB_ )
         JB = MIN( DESCA( N_ ) - J + 1, DESCA( NB_ ) )
*
         IF( MYCOL.EQ.ICURCOL .OR. COLREP ) THEN
*
            IF( MYROWDIST.EQ.0 ) THEN
               I = 1
            ELSE
               I = DESCA( IMB_ ) + ( MYROWDIST - 1 ) * DESCA( MB_ ) + 1
            END IF
*
            II = 1
            IB = MIN( MAX( 0, DESCA( M_ ) - I + 1 ), IMBA )
            DO 70 KK = 0, JB-1
               DO 60 LL = 0, IB-1
                  IF( I+LL.LT.IA .OR. I+LL.GT.IA+M-1 .OR.
     $                J+KK.LT.JA .OR. J+KK.GT.JA+N-1 )
     $               CALL PZERRSET( ERR, ERRMAX,
     $                              A( I+LL+(J+KK-1)*LDA ),
     $                              PA( II+LL+(JJ+KK-1)*LDPA ) )
   60          CONTINUE
   70       CONTINUE
            IF( ROWREP ) THEN
               I = I + IMBA
            ELSE
               I = I + IMBA + ( NPROW - 1 ) * DESCA( MB_ )
            END IF
*
            DO 100 II = IMBA+1, MPALL, DESCA( MB_ )
               IB = MIN( MPALL-II+1, DESCA( MB_ ) )
*
               DO 90 KK = 0, JB-1
                  DO 80 LL = 0, IB-1
                     IF( I+LL.LT.IA .OR. I+LL.GT.IA+M-1 .OR.
     $                   J+KK.LT.JA .OR. J+KK.GT.JA+N-1 )
     $                  CALL PZERRSET( ERR, ERRMAX,
     $                                 A( I+LL+(J+KK-1)*LDA ),
     $                                 PA( II+LL+(JJ+KK-1)*LDPA ) )
   80             CONTINUE
   90          CONTINUE
*
               IF( ROWREP ) THEN
                  I = I + DESCA( MB_ )
               ELSE
                  I = I + NPROW * DESCA( MB_ )
               END IF
*
  100       CONTINUE
*
            JJ = JJ + JB
*
         END IF
*
         ICURCOL = MOD( ICURCOL + 1, NPCOL )
*                                                           INSERT MODE
  110 CONTINUE
*
      CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERRMAX, 1, KK, LL, -1,
     $              -1, -1 )
*
      IF( ERRMAX.GT.ZERO .AND. ERRMAX.LE.EPS ) THEN
         INFO = 1
      ELSE IF( ERRMAX.GT.EPS ) THEN
         INFO = -1
      END IF
*
      RETURN
*
*     End of PZCHKMOUT
*
      END
      SUBROUTINE PZMPRNT( ICTXT, NOUT, M, N, A, LDA, IRPRNT, ICPRNT,
     $                    CMATNM )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICPRNT, ICTXT, IRPRNT, LDA, M, N, NOUT
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      CMATNM
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  PZMPRNT prints to the standard output an array A of size m by n. Only
*  the process of coordinates ( IRPRNT, ICPRNT ) is printing.
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
*  M       (global input) INTEGER
*          On entry, M  specifies the number of rows of the matrix A.  M
*          must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the matrix A.
*          N must be at least zero.
*
*  A       (local input) COMPLEX*16 array
*          On entry,  A  is an array of dimension (LDA,N). The leading m
*          by n part of this array is printed.
*
*  LDA     (local input) INTEGER
*          On entry, LDA  specifies the leading dimension of  the  local
*          array A to be printed. LDA must be at least MAX( 1, M ).
*
*  IRPRNT  (global input) INTEGER
*          On entry, IRPRNT  specifies the process row coordinate of the
*          printing process.
*
*  ICPRNT  (global input) INTEGER
*          On entry,  ICPRNT  specifies the process column coordinate of
*          the printing process.
*
*  CMATNM  (global input) CHARACTER*(*)
*          On entry, CMATNM specifies the identifier of the matrix to be
*          printed.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DIMAG
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( ( M.LE.0 ).OR.( N.LE.0 ) )
     $   RETURN
*
*     Get grid parameters
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
*
         WRITE( NOUT, FMT = * )
         DO 20 J = 1, N
*
            DO 10 I = 1, M
*
               WRITE( NOUT, FMT = 9999 ) CMATNM, I, J,
     $                         DBLE( A( I, J ) ), DIMAG( A( I, J ) )
*
   10       CONTINUE
*
   20    CONTINUE
*
      END IF
*
 9999 FORMAT( 1X, A, '(', I6, ',', I6, ')=', D30.18, '+i*(',
     $        D30.18, ')' )
*
      RETURN
*
*     End of PZMPRNT
*
      END
      SUBROUTINE PZVPRNT( ICTXT, NOUT, N, X, INCX, IRPRNT, ICPRNT,
     $                    CVECNM )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICPRNT, ICTXT, INCX, IRPRNT, N, NOUT
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      CVECNM
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  PZVPRNT  prints  to the standard output an vector x of length n. Only
*  the process of coordinates ( IRPRNT, ICPRNT ) is printing.
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
*  N       (global input) INTEGER
*          On entry, N  specifies the length of the vector X.  N must be
*          at least zero.
*
*  X       (global input) COMPLEX*16 array
*          On   entry,   X   is   an   array   of   dimension  at  least
*          ( 1 + ( n - 1 )*abs( INCX ) ).  Before  entry,  the incremen-
*          ted array X must contain the vector x.
*
*  INCX    (global input) INTEGER.
*          On entry, INCX specifies the increment for the elements of X.
*          INCX must not be zero.
*
*  IRPRNT  (global input) INTEGER
*          On entry, IRPRNT  specifies the process row coordinate of the
*          printing process.
*
*  ICPRNT  (global input) INTEGER
*          On entry,  ICPRNT  specifies the process column coordinate of
*          the printing process.
*
*  CVECNM  (global input) CHARACTER*(*)
*          On entry, CVECNM specifies the identifier of the vector to be
*          printed.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DIMAG
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Get grid parameters
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
*
         WRITE( NOUT, FMT = * )
         DO 10 I = 1, 1 + ( N-1 )*INCX, INCX
*
            WRITE( NOUT, FMT = 9999 ) CVECNM, I, DBLE( X( I ) ),
     $                                DIMAG( X( I ) )
*
   10    CONTINUE
*
      END IF
*
 9999 FORMAT( 1X, A, '(', I6, ')=', D30.18, '+i*(', D30.18, ')' )
*
      RETURN
*
*     End of PZVPRNT
*
      END
      SUBROUTINE PZMVCH( ICTXT, TRANS, M, N, ALPHA, A, IA, JA, DESCA,
     $                   X, IX, JX, DESCX, INCX, BETA, Y, PY, IY, JY,
     $                   DESCY, INCY, G, ERR, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        TRANS
      INTEGER            IA, ICTXT, INCX, INCY, INFO, IX, IY, JA, JX,
     $                   JY, M, N
      DOUBLE PRECISION   ERR
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCX( * ), DESCY( * )
      DOUBLE PRECISION   G( * )
      COMPLEX*16         A( * ), PY( * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  PZMVCH checks the results of the computational tests.
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
*  TRANS   (global input) CHARACTER*1
*          On entry,  TRANS  specifies which matrix-vector product is to
*          be computed as follows:
*             If TRANS = 'T',
*                sub( Y ) = BETA * sub( Y ) + sub( A )**T  * sub( X ),
*             else if TRANS = 'C',
*                sub( Y ) = BETA * sub( Y ) + sub( A )**H  * sub( X ),
*             otherwise
*                sub( Y ) = BETA * sub( Y ) + sub( A )     * sub( X ).
*
*  M       (global input) INTEGER
*          On entry,  M  specifies  the  number of rows of the submatrix
*          operand matrix A. M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry,  N  specifies  the  number of columns of the subma-
*          trix operand matrix A. N must be at least zero.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (local input) COMPLEX*16 array
*          On entry, A is an array of  dimension  (DESCA( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PA.
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
*  X       (local input) COMPLEX*16 array
*          On entry, X is an array of  dimension  (DESCX( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PX.
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
*  BETA    (global input) COMPLEX*16
*          On entry, BETA specifies the scalar beta.
*
*  Y       (local input/local output) COMPLEX*16 array
*          On entry, Y is an array of  dimension  (DESCY( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PY.
*
*  PY      (local input) COMPLEX*16 array
*          On entry, PY is an array of dimension (DESCY( LLD_ ),*). This
*          array contains the local entries of the matrix PY.
*
*  IY      (global input) INTEGER
*          On entry, IY  specifies Y's global row index, which points to
*          the beginning of the submatrix sub( Y ).
*
*  JY      (global input) INTEGER
*          On entry, JY  specifies Y's global column index, which points
*          to the beginning of the submatrix sub( Y ).
*
*  DESCY   (global and local input) INTEGER array
*          On entry, DESCY  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix Y.
*
*  INCY    (global input) INTEGER
*          On entry,  INCY   specifies  the  global  increment  for  the
*          elements of  Y.  Only two values of  INCY   are  supported in
*          this version, namely 1 and M_Y. INCY  must not be zero.
*
*  G       (workspace) DOUBLE PRECISION array
*          On entry, G is an array of dimension at least MAX( M, N ).  G
*          is used to compute the gauges.
*
*  ERR     (global output) DOUBLE PRECISION
*          On exit, ERR specifies the largest error in absolute value.
*
*  INFO    (global output) INTEGER
*          On exit, if INFO <> 0, the result is less than half accurate.
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
      DOUBLE PRECISION   RZERO, RONE
      PARAMETER          ( RZERO = 0.0D+0, RONE = 1.0D+0 )
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, CTRAN, ROWREP, TRAN
      INTEGER            I, IB, ICURCOL, ICURROW, IIY, IN, IOFFA, IOFFX,
     $                   IOFFY, IYCOL, IYROW, J, JB, JJY, JN, KK, LDA,
     $                   LDPY, LDX, LDY, ML, MYCOL, MYROW, NL, NPCOL,
     $                   NPROW
      DOUBLE PRECISION   EPS, ERRI, GTMP
      COMPLEX*16         C, TBETA, YTMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGAMX2D, IGSUM2D, PB_INFOG2L
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           LSAME, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN, MOD, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
      ABS1( C ) = ABS( DBLE( C ) ) + ABS( DIMAG( C ) )
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      EPS = PDLAMCH( ICTXT, 'eps' )
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         TBETA = ONE
      ELSE
         TBETA = BETA
      END IF
*
      TRAN = LSAME( TRANS, 'T' )
      CTRAN = LSAME( TRANS, 'C' )
      IF( TRAN.OR.CTRAN ) THEN
         ML = N
         NL = M
      ELSE
         ML = M
         NL = N
      END IF
*
      LDA = MAX( 1, DESCA( M_ ) )
      LDX = MAX( 1, DESCX( M_ ) )
      LDY = MAX( 1, DESCY( M_ ) )
*
*     Compute expected result in Y using data in A, X and Y.
*     Compute gauges in G. This part of the computation is performed
*     by every process in the grid.
*
      IOFFY = IY + ( JY - 1 ) * LDY
      DO 40 I = 1, ML
         YTMP = ZERO
         GTMP = RZERO
         IOFFX = IX + ( JX - 1 ) * LDX
         IF( TRAN )THEN
            IOFFA = IA + ( JA + I - 2 ) * LDA
            DO 10 J = 1, NL
               YTMP = YTMP + A( IOFFA ) * X( IOFFX )
               GTMP = GTMP + ABS1( A( IOFFA ) ) * ABS1( X( IOFFX ) )
               IOFFA = IOFFA + 1
               IOFFX = IOFFX + INCX
   10       CONTINUE
         ELSE IF( CTRAN )THEN
            IOFFA = IA + ( JA + I - 2 ) * LDA
            DO 20 J = 1, NL
               YTMP = YTMP + DCONJG( A( IOFFA ) ) * X( IOFFX )
               GTMP = GTMP + ABS1( A( IOFFA ) ) * ABS1( X( IOFFX ) )
               IOFFA = IOFFA + 1
               IOFFX = IOFFX + INCX
   20       CONTINUE
         ELSE
            IOFFA = IA + I - 1 + ( JA - 1 ) * LDA
            DO 30 J = 1, NL
               YTMP = YTMP + A( IOFFA ) * X( IOFFX )
               GTMP = GTMP + ABS1( A( IOFFA ) ) * ABS1( X( IOFFX ) )
               IOFFA = IOFFA + LDA
               IOFFX = IOFFX + INCX
   30       CONTINUE
         END IF
         G( I ) = ABS1( ALPHA )*GTMP + ABS1( TBETA )*ABS1( Y( IOFFY ) )
         Y( IOFFY ) = ALPHA * YTMP + TBETA * Y( IOFFY )
         IOFFY = IOFFY + INCY
   40 CONTINUE
*
*     Compute the error ratio for this result.
*
      ERR  = RZERO
      INFO = 0
      LDPY = DESCY( LLD_ )
      IOFFY = IY + ( JY - 1 ) * LDY
      CALL PB_INFOG2L( IY, JY, DESCY, NPROW, NPCOL, MYROW, MYCOL, IIY,
     $                 JJY, IYROW, IYCOL )
      ICURROW = IYROW
      ICURCOL = IYCOL
      ROWREP  = ( IYROW.EQ.-1 )
      COLREP  = ( IYCOL.EQ.-1 )
*
      IF( INCY.EQ.DESCY( M_ ) ) THEN
*
*        sub( Y ) is a row vector
*
         JB = DESCY( INB_ ) - JY + 1
         IF( JB.LE.0 )
     $      JB = ( ( -JB ) / DESCY( NB_ ) + 1 ) * DESCY( NB_ ) + JB
         JB = MIN( JB, ML )
         JN = JY + JB - 1
*
         DO 50 J = JY, JN
*
            IF( ( MYROW.EQ.ICURROW .OR. ROWREP ) .AND.
     $          ( MYCOL.EQ.ICURCOL .OR. COLREP ) ) THEN
               ERRI = ABS( PY( IIY+(JJY-1)*LDPY ) - Y( IOFFY ) ) / EPS
               IF( G( J-JY+1 ).NE.RZERO )
     $            ERRI = ERRI / G( J-JY+1 )
               ERR = MAX( ERR, ERRI )
               IF( ERR*SQRT( EPS ).GE.RONE )
     $            INFO = 1
               JJY = JJY + 1
            END IF
*
            IOFFY = IOFFY + INCY
*
   50    CONTINUE
*
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
         DO 70 J = JN+1, JY+ML-1, DESCY( NB_ )
            JB = MIN( JY+ML-J, DESCY( NB_ ) )
*
            DO 60 KK = 0, JB-1
*
               IF( ( MYROW.EQ.ICURROW .OR. ROWREP ) .AND.
     $             ( MYCOL.EQ.ICURCOL .OR. COLREP ) ) THEN
                  ERRI = ABS( PY( IIY+(JJY-1)*LDPY ) - Y( IOFFY ) )/EPS
                  IF( G( J+KK-JY+1 ).NE.RZERO )
     $               ERRI = ERRI / G( J+KK-JY+1 )
                  ERR = MAX( ERR, ERRI )
                  IF( ERR*SQRT( EPS ).GE.RONE )
     $               INFO = 1
                  JJY = JJY + 1
               END IF
*
               IOFFY = IOFFY + INCY
*
   60       CONTINUE
*
            ICURCOL = MOD( ICURCOL+1, NPCOL )
*
   70    CONTINUE
*
      ELSE
*
*        sub( Y ) is a column vector
*
         IB = DESCY( IMB_ ) - IY + 1
         IF( IB.LE.0 )
     $      IB = ( ( -IB ) / DESCY( MB_ ) + 1 ) * DESCY( MB_ ) + IB
         IB = MIN( IB, ML )
         IN = IY + IB - 1
*
         DO 80 I = IY, IN
*
            IF( ( MYROW.EQ.ICURROW .OR. ROWREP ) .AND.
     $          ( MYCOL.EQ.ICURCOL .OR. COLREP ) ) THEN
               ERRI = ABS( PY( IIY+(JJY-1)*LDPY ) - Y( IOFFY ) ) / EPS
               IF( G( I-IY+1 ).NE.RZERO )
     $            ERRI = ERRI / G( I-IY+1 )
               ERR = MAX( ERR, ERRI )
               IF( ERR*SQRT( EPS ).GE.RONE )
     $            INFO = 1
               IIY = IIY + 1
            END IF
*
            IOFFY = IOFFY + INCY
*
   80    CONTINUE
*
         ICURROW = MOD( ICURROW+1, NPROW )
*
         DO 100 I = IN+1, IY+ML-1, DESCY( MB_ )
            IB = MIN( IY+ML-I, DESCY( MB_ ) )
*
            DO 90 KK = 0, IB-1
*
               IF( ( MYROW.EQ.ICURROW .OR. ROWREP ) .AND.
     $             ( MYCOL.EQ.ICURCOL .OR. COLREP ) ) THEN
                  ERRI = ABS( PY( IIY+(JJY-1)*LDPY ) - Y( IOFFY ) )/EPS
                  IF( G( I+KK-IY+1 ).NE.RZERO )
     $               ERRI = ERRI / G( I+KK-IY+1 )
                  ERR = MAX( ERR, ERRI )
                  IF( ERR*SQRT( EPS ).GE.RONE )
     $               INFO = 1
                  IIY = IIY + 1
               END IF
*
               IOFFY = IOFFY + INCY
*
   90       CONTINUE
*
            ICURROW = MOD( ICURROW+1, NPROW )
*
  100    CONTINUE
*
      END IF
*
*     If INFO = 0, all results are at least half accurate.
*
      CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, MYCOL )
      CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERR, 1, I, J, -1, -1,
     $              MYCOL )
*
      RETURN
*
*     End of PZMVCH
*
      END
      SUBROUTINE PZVMCH( ICTXT, TRANS, UPLO, M, N, ALPHA, X, IX, JX,
     $                     DESCX, INCX, Y, IY, JY, DESCY, INCY, A, PA,
     $                     IA, JA, DESCA, G, ERR, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        TRANS, UPLO
      INTEGER            IA, ICTXT, INCX, INCY, INFO, IX, IY, JA, JX,
     $                   JY, M, N
      DOUBLE PRECISION   ERR
      COMPLEX*16         ALPHA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCX( * ), DESCY( * )
      DOUBLE PRECISION   G( * )
      COMPLEX*16         A( * ), PA( * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  PZVMCH checks the results of the computational tests.
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
*  TRANS   (global input) CHARACTER*1
*          On entry,  TRANS  specifies  the operation to be performed in
*          the complex cases:
*             if TRANS = 'C',
*                sub( A ) := sub( A ) + alpha * sub( X ) * sub( Y )**H,
*             otherwise
*                sub( A ) := sub( A ) + alpha * sub( X ) * sub( Y )**T.
*
*  UPLO    (global input) CHARACTER*1
*          On entry, UPLO specifies which part of the submatrix sub( A )
*          is to be referenced as follows:
*             If UPLO = 'L', only the lower triangular part,
*             If UPLO = 'U', only the upper triangular part,
*             else the entire matrix is to be referenced.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies  the  number of rows of the submatrix
*          operand matrix A. M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry,  N  specifies  the  number of columns of the subma-
*          trix operand matrix A. N must be at least zero.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  X       (local input) COMPLEX*16 array
*          On entry, X is an array of  dimension  (DESCX( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PX.
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
*  Y       (local input) COMPLEX*16 array
*          On entry, Y is an array of  dimension  (DESCY( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PY.
*
*  IY      (global input) INTEGER
*          On entry, IY  specifies Y's global row index, which points to
*          the beginning of the submatrix sub( Y ).
*
*  JY      (global input) INTEGER
*          On entry, JY  specifies Y's global column index, which points
*          to the beginning of the submatrix sub( Y ).
*
*  DESCY   (global and local input) INTEGER array
*          On entry, DESCY  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix Y.
*
*  INCY    (global input) INTEGER
*          On entry,  INCY   specifies  the  global  increment  for  the
*          elements of  Y.  Only two values of  INCY   are  supported in
*          this version, namely 1 and M_Y. INCY  must not be zero.
*
*  A       (local input/local output) COMPLEX*16 array
*          On entry, A is an array of  dimension  (DESCA( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PA.
*
*  PA      (local input) COMPLEX*16 array
*          On entry, PA is an array of dimension (DESCA( LLD_ ),*). This
*          array contains the local entries of the matrix PA.
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
*  G       (workspace) DOUBLE PRECISION array
*          On entry, G is an array of dimension at least MAX( M, N ).  G
*          is used to compute the gauges.
*
*  ERR     (global output) DOUBLE PRECISION
*          On exit, ERR specifies the largest error in absolute value.
*
*  INFO    (global output) INTEGER
*          On exit, if INFO <> 0, the result is less than half accurate.
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
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, CTRAN, LOWER, ROWREP, UPPER
      INTEGER            I, IACOL, IAROW, IB, IBEG, ICURROW, IEND, IIA,
     $                   IN, IOFFA, IOFFX, IOFFY, J, JJA, KK, LDA, LDPA,
     $                   LDX, LDY, MYCOL, MYROW, NPCOL, NPROW
      DOUBLE PRECISION   EPS, ERRI, GTMP
      COMPLEX*16         ATMP, C
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGAMX2D, IGSUM2D, PB_INFOG2L
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           LSAME, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN, MOD, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
      ABS1( C ) = ABS( DBLE( C ) ) + ABS( DIMAG( C ) )
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      EPS = PDLAMCH( ICTXT, 'eps' )
*
      CTRAN = LSAME( TRANS, 'C' )
      UPPER = LSAME( UPLO, 'U' )
      LOWER = LSAME( UPLO, 'L' )
*
      LDA = MAX( 1, DESCA( M_ ) )
      LDX = MAX( 1, DESCX( M_ ) )
      LDY = MAX( 1, DESCY( M_ ) )
*
*     Compute expected result in A using data in A, X and Y.
*     Compute gauges in G. This part of the computation is performed
*     by every process in the grid.
*
      DO 70 J = 1, N
*
         IOFFY = IY + ( JY - 1 ) * LDY + ( J - 1 ) * INCY
*
         IF( LOWER ) THEN
            IBEG = J
            IEND = M
            DO 10 I = 1, J-1
               G( I ) = ZERO
   10       CONTINUE
         ELSE IF( UPPER ) THEN
            IBEG = 1
            IEND = J
            DO 20 I = J+1, M
               G( I ) = ZERO
   20       CONTINUE
         ELSE
            IBEG = 1
            IEND = M
         END IF
*
         DO 30 I = IBEG, IEND
*
            IOFFX = IX + ( JX - 1 ) * LDX + ( I - 1 ) * INCX
            IOFFA = IA + I - 1 + ( JA + J - 2 ) * LDA
            IF( CTRAN ) THEN
               ATMP = X( IOFFX ) * DCONJG( Y( IOFFY ) )
            ELSE
               ATMP = X( IOFFX ) * Y( IOFFY )
            END IF
            GTMP = ABS1( X( IOFFX ) ) * ABS1( Y( IOFFY ) )
            G( I ) = ABS1( ALPHA ) * GTMP + ABS1( A( IOFFA ) )
            A( IOFFA ) = ALPHA * ATMP + A( IOFFA )
*
   30    CONTINUE
*
*        Compute the error ratio for this result.
*
         INFO = 0
         ERR  = ZERO
         LDPA = DESCA( LLD_ )
         IOFFA = IA + ( JA + J - 2 ) * LDA
         CALL PB_INFOG2L( IA, JA+J-1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIA, JJA, IAROW, IACOL )
         ROWREP = ( IAROW.EQ.-1 )
         COLREP = ( IACOL.EQ.-1 )
*
         IF( MYCOL.EQ.IACOL .OR. COLREP ) THEN
*
            ICURROW = IAROW
            IB = DESCA( IMB_ ) - IA + 1
            IF( IB.LE.0 )
     $         IB = ( ( -IB ) / DESCA( MB_ ) + 1 ) * DESCA( MB_ ) + IB
            IB = MIN( IB, M )
            IN = IA + IB - 1
*
            DO 40 I = IA, IN
*
               IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                  ERRI = ABS( PA( IIA+(JJA-1)*LDPA ) - A( IOFFA ) )/EPS
                  IF( G( I-IA+1 ).NE.ZERO )
     $               ERRI = ERRI / G( I-IA+1 )
                  ERR = MAX( ERR, ERRI )
                  IF( ERR*SQRT( EPS ).GE.ONE )
     $               INFO = 1
                  IIA = IIA + 1
               END IF
*
               IOFFA = IOFFA + 1
*
   40       CONTINUE
*
            ICURROW = MOD( ICURROW+1, NPROW )
*
            DO 60 I = IN+1, IA+M-1, DESCA( MB_ )
               IB = MIN( IA+M-I, DESCA( MB_ ) )
*
               DO 50 KK = 0, IB-1
*
                  IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                     ERRI = ABS( PA( IIA+(JJA-1)*LDPA )-A( IOFFA ) )/EPS
                     IF( G( I+KK-IA+1 ).NE.ZERO )
     $                  ERRI = ERRI / G( I+KK-IA+1 )
                     ERR = MAX( ERR, ERRI )
                     IF( ERR*SQRT( EPS ).GE.ONE )
     $                  INFO = 1
                     IIA = IIA + 1
                  END IF
*
                  IOFFA = IOFFA + 1
*
   50          CONTINUE
*
               ICURROW = MOD( ICURROW+1, NPROW )
*
   60       CONTINUE
*
         END IF
*
*        If INFO = 0, all results are at least half accurate.
*
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, MYCOL )
         CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERR, 1, I, J, -1, -1,
     $                 MYCOL )
         IF( INFO.NE.0 )
     $      GO TO 80
*
   70 CONTINUE
*
   80 CONTINUE
*
      RETURN
*
*     End of PZVMCH
*
      END
      SUBROUTINE PZVMCH2( ICTXT, UPLO, M, N, ALPHA, X, IX, JX, DESCX,
     $                    INCX, Y, IY, JY, DESCY, INCY, A, PA, IA,
     $                    JA, DESCA, G, ERR, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO
      INTEGER            IA, ICTXT, INCX, INCY, INFO, IX, IY, JA, JX,
     $                   JY, M, N
      DOUBLE PRECISION   ERR
      COMPLEX*16         ALPHA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCX( * ), DESCY( * )
      DOUBLE PRECISION   G( * )
      COMPLEX*16         A( * ), PA( * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  PZVMCH2 checks the results of the computational tests.
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
*  UPLO    (global input) CHARACTER*1
*          On entry, UPLO specifies which part of the submatrix sub( A )
*          is to be referenced as follows:
*             If UPLO = 'L', only the lower triangular part,
*             If UPLO = 'U', only the upper triangular part,
*             else the entire matrix is to be referenced.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies  the  number of rows of the submatrix
*          operand matrix A. M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry,  N  specifies  the  number of columns of the subma-
*          trix operand matrix A. N must be at least zero.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  X       (local input) COMPLEX*16 array
*          On entry, X is an array of  dimension  (DESCX( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PX.
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
*  Y       (local input) COMPLEX*16 array
*          On entry, Y is an array of  dimension  (DESCY( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PY.
*
*  IY      (global input) INTEGER
*          On entry, IY  specifies Y's global row index, which points to
*          the beginning of the submatrix sub( Y ).
*
*  JY      (global input) INTEGER
*          On entry, JY  specifies Y's global column index, which points
*          to the beginning of the submatrix sub( Y ).
*
*  DESCY   (global and local input) INTEGER array
*          On entry, DESCY  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix Y.
*
*  INCY    (global input) INTEGER
*          On entry,  INCY   specifies  the  global  increment  for  the
*          elements of  Y.  Only two values of  INCY   are  supported in
*          this version, namely 1 and M_Y. INCY  must not be zero.
*
*  A       (local input/local output) COMPLEX*16 array
*          On entry, A is an array of  dimension  (DESCA( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PA.
*
*  PA      (local input) COMPLEX*16 array
*          On entry, PA is an array of dimension (DESCA( LLD_ ),*). This
*          array contains the local entries of the matrix PA.
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
*  G       (workspace) DOUBLE PRECISION array
*          On entry, G is an array of dimension at least MAX( M, N ).  G
*          is used to compute the gauges.
*
*  ERR     (global output) DOUBLE PRECISION
*          On exit, ERR specifies the largest error in absolute value.
*
*  INFO    (global output) INTEGER
*          On exit, if INFO <> 0, the result is less than half accurate.
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
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, LOWER, ROWREP, UPPER
      INTEGER            I, IACOL, IAROW, IB, IBEG, ICURROW, IEND, IIA,
     $                   IN, IOFFA, IOFFXI, IOFFXJ, IOFFYI, IOFFYJ, J,
     $                   JJA, KK, LDA, LDPA, LDX, LDY, MYCOL, MYROW,
     $                   NPCOL, NPROW
      DOUBLE PRECISION   EPS, ERRI, GTMP
      COMPLEX*16         C, ATMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGAMX2D, IGSUM2D, PB_INFOG2L
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           LSAME, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN, MOD, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
      ABS1( C ) = ABS( DBLE( C ) ) + ABS( DIMAG( C ) )
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      EPS = PDLAMCH( ICTXT, 'eps' )
*
      UPPER = LSAME( UPLO, 'U' )
      LOWER = LSAME( UPLO, 'L' )
*
      LDA = MAX( 1, DESCA( M_ ) )
      LDX = MAX( 1, DESCX( M_ ) )
      LDY = MAX( 1, DESCY( M_ ) )
*
*     Compute expected result in A using data in A, X and Y.
*     Compute gauges in G. This part of the computation is performed
*     by every process in the grid.
*
      DO 70 J = 1, N
*
         IOFFXJ = IX + ( JX - 1 ) * LDX + ( J - 1 ) * INCX
         IOFFYJ = IY + ( JY - 1 ) * LDY + ( J - 1 ) * INCY
*
         IF( LOWER ) THEN
            IBEG = J
            IEND = M
            DO 10 I = 1, J-1
               G( I ) = ZERO
   10       CONTINUE
         ELSE IF( UPPER ) THEN
            IBEG = 1
            IEND = J
            DO 20 I = J+1, M
               G( I ) = ZERO
   20       CONTINUE
         ELSE
            IBEG = 1
            IEND = M
         END IF
*
         DO 30 I = IBEG, IEND
            IOFFA = IA + I - 1 + ( JA + J - 2 ) * LDA
            IOFFXI = IX + ( JX - 1 ) * LDX + ( I - 1 ) * INCX
            IOFFYI = IY + ( JY - 1 ) * LDY + ( I - 1 ) * INCY
            ATMP = ALPHA * X( IOFFXI ) * DCONJG( Y( IOFFYJ ) )
            ATMP = ATMP + Y( IOFFYI ) * DCONJG( ALPHA * X( IOFFXJ ) )
            GTMP = ABS1( ALPHA * X( IOFFXI ) ) * ABS1( Y( IOFFYJ ) )
            GTMP = GTMP + ABS1( Y( IOFFYI ) ) *
     $                    ABS1( DCONJG( ALPHA * X( IOFFXJ ) ) )
            G( I ) = GTMP + ABS1( A( IOFFA ) )
            A( IOFFA ) = A( IOFFA ) + ATMP
*
   30    CONTINUE
*
*        Compute the error ratio for this result.
*
         INFO = 0
         ERR  = ZERO
         LDPA = DESCA( LLD_ )
         IOFFA = IA + ( JA + J - 2 ) * LDA
         CALL PB_INFOG2L( IA, JA+J-1, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIA, JJA, IAROW, IACOL )
         ROWREP = ( IAROW.EQ.-1 )
         COLREP = ( IACOL.EQ.-1 )
*
         IF( MYCOL.EQ.IACOL .OR. COLREP ) THEN
*
            ICURROW = IAROW
            IB = DESCA( IMB_ ) - IA + 1
            IF( IB.LE.0 )
     $         IB = ( ( -IB ) / DESCA( MB_ ) + 1 ) * DESCA( MB_ ) + IB
            IB = MIN( IB, M )
            IN = IA + IB - 1
*
            DO 40 I = IA, IN
*
               IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                  ERRI = ABS( PA( IIA+(JJA-1)*LDPA ) - A( IOFFA ) )/EPS
                  IF( G( I-IA+1 ).NE.ZERO )
     $               ERRI = ERRI / G( I-IA+1 )
                  ERR = MAX( ERR, ERRI )
                  IF( ERR*SQRT( EPS ).GE.ONE )
     $               INFO = 1
                  IIA = IIA + 1
               END IF
*
               IOFFA = IOFFA + 1
*
   40       CONTINUE
*
            ICURROW = MOD( ICURROW+1, NPROW )
*
            DO 60 I = IN+1, IA+M-1, DESCA( MB_ )
               IB = MIN( IA+M-I, DESCA( MB_ ) )
*
               DO 50 KK = 0, IB-1
*
                  IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                     ERRI = ABS( PA( IIA+(JJA-1)*LDPA )-A( IOFFA ) )/EPS
                     IF( G( I+KK-IA+1 ).NE.ZERO )
     $                  ERRI = ERRI / G( I+KK-IA+1 )
                     ERR = MAX( ERR, ERRI )
                     IF( ERR*SQRT( EPS ).GE.ONE )
     $                  INFO = 1
                     IIA = IIA + 1
                  END IF
*
                  IOFFA = IOFFA + 1
*
   50          CONTINUE
*
               ICURROW = MOD( ICURROW+1, NPROW )
*
   60       CONTINUE
*
         END IF
*
*        If INFO = 0, all results are at least half accurate.
*
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, MYCOL )
         CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERR, 1, I, J, -1, -1,
     $                 MYCOL )
         IF( INFO.NE.0 )
     $      GO TO 80
*
   70 CONTINUE
*
   80 CONTINUE
*
      RETURN
*
*     End of PZVMCH2
*
      END
      SUBROUTINE PZMMCH( ICTXT, TRANSA, TRANSB, M, N, K, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB, BETA, C, PC, IC,
     $                   JC, DESCC, CT, G, ERR, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            IA, IB, IC, ICTXT, INFO, JA, JB, JC, K, M, N
      DOUBLE PRECISION   ERR
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCC( * )
      DOUBLE PRECISION   G( * )
      COMPLEX*16         A( * ), B( * ), C( * ), CT( * ), PC( * )
*     ..
*
*  Purpose
*  =======
*
*  PZMMCH checks the results of the computational tests.
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
*  TRANSA  (global input) CHARACTER*1
*          On entry, TRANSA specifies if the matrix  operand  A is to be
*          transposed.
*
*  TRANSB  (global input) CHARACTER*1
*          On entry, TRANSB specifies if the matrix  operand  B is to be
*          transposed.
*
*  M       (global input) INTEGER
*          On entry, M specifies the number of rows of C.
*
*  N       (global input) INTEGER
*          On entry, N specifies the number of columns of C.
*
*  K       (global input) INTEGER
*          On entry, K specifies the number of columns (resp. rows) of A
*          when  TRANSA = 'N'  (resp. TRANSA <> 'N')  in PxGEMM, PxSYRK,
*          PxSYR2K, PxHERK and PxHER2K.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (local input) COMPLEX*16 array
*          On entry, A is an array of  dimension  (DESCA( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PA.
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
*  B       (local input) COMPLEX*16 array
*          On entry, B is an array of  dimension  (DESCB( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PB.
*
*  IB      (global input) INTEGER
*          On entry, IB  specifies B's global row index, which points to
*          the beginning of the submatrix sub( B ).
*
*  JB      (global input) INTEGER
*          On entry, JB  specifies B's global column index, which points
*          to the beginning of the submatrix sub( B ).
*
*  DESCB   (global and local input) INTEGER array
*          On entry, DESCB  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix B.
*
*  BETA    (global input) COMPLEX*16
*          On entry, BETA specifies the scalar beta.
*
*  C       (local input/local output) COMPLEX*16 array
*          On entry, C is an array of  dimension  (DESCC( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PC.
*
*  PC      (local input) COMPLEX*16 array
*          On entry, PC is an array of dimension (DESCC( LLD_ ),*). This
*          array contains the local pieces of the matrix PC.
*
*  IC      (global input) INTEGER
*          On entry, IC  specifies C's global row index, which points to
*          the beginning of the submatrix sub( C ).
*
*  JC      (global input) INTEGER
*          On entry, JC  specifies C's global column index, which points
*          to the beginning of the submatrix sub( C ).
*
*  DESCC   (global and local input) INTEGER array
*          On entry, DESCC  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix C.
*
*  CT      (workspace) COMPLEX*16 array
*          On entry, CT is an array of dimension at least MAX(M,N,K). CT
*          holds a copy of the current column of C.
*
*  G       (workspace) DOUBLE PRECISION array
*          On entry, G  is  an array of dimension at least MAX(M,N,K). G
*          is used to compute the gauges.
*
*  ERR     (global output) DOUBLE PRECISION
*          On exit, ERR specifies the largest error in absolute value.
*
*  INFO    (global output) INTEGER
*          On exit, if INFO <> 0, the result is less than half accurate.
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
      DOUBLE PRECISION   RZERO, RONE
      PARAMETER          ( RZERO = 0.0D+0, RONE = 1.0D+0 )
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, CTRANA, CTRANB, ROWREP, TRANA, TRANB
      INTEGER            I, IBB, ICCOL, ICROW, ICURROW, IIC, IN, IOFFA,
     $                   IOFFB, IOFFC, J, JJC, KK, LDA, LDB, LDC, LDPC,
     $                   MYCOL, MYROW, NPCOL, NPROW
      DOUBLE PRECISION   EPS, ERRI
      COMPLEX*16         Z
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGAMX2D, IGSUM2D, PB_INFOG2L
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           LSAME, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN, MOD, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
      ABS1( Z ) = ABS( DBLE( Z ) ) + ABS( DIMAG( Z ) )
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      EPS = PDLAMCH( ICTXT, 'eps' )
*
      TRANA = LSAME( TRANSA, 'T' ).OR.LSAME( TRANSA, 'C' )
      TRANB = LSAME( TRANSB, 'T' ).OR.LSAME( TRANSB, 'C' )
      CTRANA = LSAME( TRANSA, 'C' )
      CTRANB = LSAME( TRANSB, 'C' )
*
      LDA = MAX( 1, DESCA( M_ ) )
      LDB = MAX( 1, DESCB( M_ ) )
      LDC = MAX( 1, DESCC( M_ ) )
*
*     Compute expected result in C using data in A, B and C.
*     Compute gauges in G. This part of the computation is performed
*     by every process in the grid.
*
      DO 240 J = 1, N
*
         IOFFC = IC + ( JC + J - 2 ) * LDC
         DO 10 I = 1, M
            CT( I ) = ZERO
            G( I )  = RZERO
   10    CONTINUE
*
         IF( .NOT.TRANA .AND. .NOT.TRANB ) THEN
            DO 30 KK = 1, K
               IOFFB = IB + KK - 1 + ( JB + J - 2 ) * LDB
               DO 20 I = 1, M
                  IOFFA = IA + I - 1 + ( JA + KK - 2 ) * LDA
                  CT( I ) = CT( I ) + A( IOFFA ) * B( IOFFB )
                  G( I ) = G( I ) + ABS( A( IOFFA ) ) *
     $                     ABS( B( IOFFB ) )
   20          CONTINUE
   30       CONTINUE
         ELSE IF( TRANA .AND. .NOT.TRANB ) THEN
            IF( CTRANA ) THEN
               DO 50 KK = 1, K
                  IOFFB = IB + KK - 1 + ( JB + J - 2 ) * LDB
                  DO 40 I = 1, M
                     IOFFA = IA + KK - 1 + ( JA + I - 2 ) * LDA
                     CT( I ) = CT( I ) + DCONJG( A( IOFFA ) ) *
     $                                   B( IOFFB )
                     G( I ) = G( I ) + ABS1( A( IOFFA ) ) *
     $                        ABS1( B( IOFFB ) )
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 70 KK = 1, K
                  IOFFB = IB + KK - 1 + ( JB + J - 2 ) * LDB
                  DO 60 I = 1, M
                     IOFFA = IA + KK - 1 + ( JA + I - 2 ) * LDA
                     CT( I ) = CT( I ) + A( IOFFA ) * B( IOFFB )
                     G( I ) = G( I ) + ABS1( A( IOFFA ) ) *
     $                        ABS1( B( IOFFB ) )
   60             CONTINUE
   70          CONTINUE
            END IF
         ELSE IF( .NOT.TRANA .AND. TRANB ) THEN
            IF( CTRANB ) THEN
               DO 90 KK = 1, K
                  IOFFB = IB + J - 1 + ( JB + KK - 2 ) * LDB
                  DO 80 I = 1, M
                     IOFFA = IA + I - 1 + ( JA + KK - 2 ) * LDA
                     CT( I ) = CT( I ) + A( IOFFA ) *
     $                                   DCONJG( B( IOFFB ) )
                     G( I ) = G( I ) + ABS1( A( IOFFA ) ) *
     $                        ABS1( B( IOFFB ) )
   80             CONTINUE
   90          CONTINUE
            ELSE
               DO 110 KK = 1, K
                  IOFFB = IB + J - 1 + ( JB + KK - 2 ) * LDB
                  DO 100 I = 1, M
                     IOFFA = IA + I - 1 + ( JA + KK - 2 ) * LDA
                     CT( I ) = CT( I ) + A( IOFFA ) * B( IOFFB )
                     G( I ) = G( I ) + ABS1( A( IOFFA ) ) *
     $                        ABS1( B( IOFFB ) )
  100             CONTINUE
  110          CONTINUE
            END IF
         ELSE IF( TRANA .AND. TRANB ) THEN
            IF( CTRANA ) THEN
               IF( CTRANB ) THEN
                  DO 130 KK = 1, K
                     IOFFB = IB + J - 1 + ( JB + KK - 2 ) * LDB
                     DO 120 I = 1, M
                        IOFFA = IA + KK - 1 + ( JA + I - 2 ) * LDA
                        CT( I ) = CT( I ) + DCONJG( A( IOFFA ) ) *
     $                                      DCONJG( B( IOFFB ) )
                        G( I ) = G( I ) + ABS1( A( IOFFA ) ) *
     $                           ABS1( B( IOFFB ) )
  120                CONTINUE
  130             CONTINUE
               ELSE
                  DO 150 KK = 1, K
                     IOFFB = IB + J - 1 + ( JB + KK - 2 ) * LDB
                     DO 140 I = 1, M
                        IOFFA = IA + KK - 1 + ( JA + I - 2 ) * LDA
                        CT( I ) = CT( I ) + DCONJG( A( IOFFA ) ) *
     $                                      B( IOFFB )
                        G( I ) = G( I ) + ABS1( A( IOFFA ) ) *
     $                           ABS1( B( IOFFB ) )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
               IF( CTRANB ) THEN
                  DO 170 KK = 1, K
                     IOFFB = IB + J - 1 + ( JB + KK - 2 ) * LDB
                     DO 160 I = 1, M
                        IOFFA = IA + KK - 1 + ( JA + I - 2 ) * LDA
                        CT( I ) = CT( I ) + A( IOFFA ) *
     $                                      DCONJG( B( IOFFB ) )
                        G( I ) = G( I ) + ABS1( A( IOFFA ) ) *
     $                           ABS1( B( IOFFB ) )
  160                CONTINUE
  170             CONTINUE
               ELSE
                  DO 190 KK = 1, K
                     IOFFB = IB + J - 1 + ( JB + KK - 2 ) * LDB
                     DO 180 I = 1, M
                        IOFFA = IA + KK - 1 + ( JA + I - 2 ) * LDA
                        CT( I ) = CT( I ) + A( IOFFA ) * B( IOFFB )
                        G( I ) = G( I ) + ABS1( A( IOFFA ) ) *
     $                           ABS1( B( IOFFB ) )
  180                CONTINUE
  190             CONTINUE
               END IF
            END IF
         END IF
*
         DO 200 I = 1, M
            CT( I ) = ALPHA*CT( I ) + BETA * C( IOFFC )
            G( I ) = ABS1( ALPHA )*G( I ) +
     $               ABS1( BETA )*ABS1( C( IOFFC ) )
            C( IOFFC ) = CT( I )
            IOFFC      = IOFFC + 1
  200    CONTINUE
*
*        Compute the error ratio for this result.
*
         ERR  = RZERO
         INFO = 0
         LDPC = DESCC( LLD_ )
         IOFFC = IC + ( JC + J - 2 ) * LDC
         CALL PB_INFOG2L( IC, JC+J-1, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIC, JJC, ICROW, ICCOL )
         ICURROW = ICROW
         ROWREP  = ( ICROW.EQ.-1 )
         COLREP  = ( ICCOL.EQ.-1 )
*
         IF( MYCOL.EQ.ICCOL .OR. COLREP ) THEN
*
            IBB = DESCC( IMB_ ) - IC + 1
            IF( IBB.LE.0 )
     $         IBB = ( ( -IBB ) / DESCC( MB_ ) + 1 )*DESCC( MB_ ) + IBB
            IBB = MIN( IBB, M )
            IN = IC + IBB - 1
*
            DO 210 I = IC, IN
*
               IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                  ERRI = ABS( PC( IIC+(JJC-1)*LDPC ) -
     $                        C( IOFFC ) ) / EPS
                  IF( G( I-IC+1 ).NE.RZERO )
     $               ERRI = ERRI / G( I-IC+1 )
                  ERR = MAX( ERR, ERRI )
                  IF( ERR*SQRT( EPS ).GE.RONE )
     $               INFO = 1
                  IIC = IIC + 1
               END IF
*
               IOFFC = IOFFC + 1
*
  210       CONTINUE
*
            ICURROW = MOD( ICURROW+1, NPROW )
*
            DO 230 I = IN+1, IC+M-1, DESCC( MB_ )
               IBB = MIN( IC+M-I, DESCC( MB_ ) )
*
               DO 220 KK = 0, IBB-1
*
                  IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                     ERRI = ABS( PC( IIC+(JJC-1)*LDPC ) -
     $                           C( IOFFC ) )/EPS
                     IF( G( I+KK-IC+1 ).NE.RZERO )
     $                  ERRI = ERRI / G( I+KK-IC+1 )
                     ERR = MAX( ERR, ERRI )
                     IF( ERR*SQRT( EPS ).GE.RONE )
     $                  INFO = 1
                     IIC = IIC + 1
                  END IF
*
                  IOFFC = IOFFC + 1
*
  220          CONTINUE
*
               ICURROW = MOD( ICURROW+1, NPROW )
*
  230       CONTINUE
*
         END IF
*
*        If INFO = 0, all results are at least half accurate.
*
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, MYCOL )
         CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERR, 1, I, J, -1, -1,
     $                 MYCOL )
         IF( INFO.NE.0 )
     $      GO TO 250
*
  240 CONTINUE
*
  250 CONTINUE
*
      RETURN
*
*     End of PZMMCH
*
      END
      SUBROUTINE PZMMCH1( ICTXT, UPLO, TRANS, N, K, ALPHA, A, IA, JA,
     $                    DESCA, BETA, C, PC, IC, JC, DESCC, CT, G,
     $                    ERR, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        TRANS, UPLO
      INTEGER            IA, IC, ICTXT, INFO, JA, JC, K, N
      DOUBLE PRECISION   ERR
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCC( * )
      DOUBLE PRECISION   G( * )
      COMPLEX*16         A( * ), C( * ), CT( * ), PC( * )
*     ..
*
*  Purpose
*  =======
*
*  PZMMCH1 checks the results of the computational tests.
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
*  UPLO    (global input) CHARACTER*1
*          On entry,  UPLO  specifies which part of C should contain the
*          result.
*
*  TRANS   (global input) CHARACTER*1
*          On entry,  TRANS  specifies  whether  the  matrix A has to be
*          transposed or not before computing the matrix-matrix product.
*
*  N       (global input) INTEGER
*          On entry, N  specifies  the order  the submatrix operand C. N
*          must be at least zero.
*
*  K       (global input) INTEGER
*          On entry, K specifies the number of columns (resp. rows) of A
*          when  TRANS = 'N'  (resp. TRANS <> 'N').  K  must be at least
*          zero.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (local input) COMPLEX*16 array
*          On entry, A is an array of  dimension  (DESCA( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PA.
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
*  BETA    (global input) COMPLEX*16
*          On entry, BETA specifies the scalar beta.
*
*  C       (local input/local output) COMPLEX*16 array
*          On entry, C is an array of  dimension  (DESCC( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PC.
*
*  PC      (local input) COMPLEX*16 array
*          On entry, PC is an array of dimension (DESCC( LLD_ ),*). This
*          array contains the local pieces of the matrix PC.
*
*  IC      (global input) INTEGER
*          On entry, IC  specifies C's global row index, which points to
*          the beginning of the submatrix sub( C ).
*
*  JC      (global input) INTEGER
*          On entry, JC  specifies C's global column index, which points
*          to the beginning of the submatrix sub( C ).
*
*  DESCC   (global and local input) INTEGER array
*          On entry, DESCC  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix C.
*
*  CT      (workspace) COMPLEX*16 array
*          On entry, CT is an array of dimension at least MAX(M,N,K). CT
*          holds a copy of the current column of C.
*
*  G       (workspace) DOUBLE PRECISION array
*          On entry, G  is  an array of dimension at least MAX(M,N,K). G
*          is used to compute the gauges.
*
*  ERR     (global output) DOUBLE PRECISION
*          On exit, ERR specifies the largest error in absolute value.
*
*  INFO    (global output) INTEGER
*          On exit, if INFO <> 0, the result is less than half accurate.
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
      DOUBLE PRECISION   RZERO, RONE
      PARAMETER          ( RZERO = 0.0D+0, RONE = 1.0D+0 )
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, HTRAN, NOTRAN, ROWREP, TRAN, UPPER
      INTEGER            I, IBB, IBEG, ICCOL, ICROW, ICURROW, IEND, IIC,
     $                   IN, IOFFAK, IOFFAN, IOFFC, J, JJC, KK, LDA,
     $                   LDC, LDPC, MYCOL, MYROW, NPCOL, NPROW
      DOUBLE PRECISION   EPS, ERRI
      COMPLEX*16         Z
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGAMX2D, IGSUM2D, PB_INFOG2L
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           LSAME, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN, MOD, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
      ABS1( Z ) = ABS( DBLE( Z ) ) + ABS( DIMAG( Z ) )
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      EPS = PDLAMCH( ICTXT, 'eps' )
*
      UPPER  = LSAME( UPLO,  'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      TRAN   = LSAME( TRANS, 'T' )
      HTRAN  = LSAME( TRANS, 'H' )
*
      LDA = MAX( 1, DESCA( M_ ) )
      LDC = MAX( 1, DESCC( M_ ) )
*
*     Compute expected result in C using data in A, B and C.
*     Compute gauges in G. This part of the computation is performed
*     by every process in the grid.
*
      DO 140 J = 1, N
*
         IF( UPPER ) THEN
            IBEG = 1
            IEND = J
         ELSE
            IBEG = J
            IEND = N
         END IF
*
         DO 10 I = 1, N
            CT( I ) = ZERO
            G( I )  = RZERO
   10    CONTINUE
*
         IF( NOTRAN ) THEN
            DO 30 KK = 1, K
               IOFFAK = IA + J - 1 + ( JA + KK - 2 ) * LDA
               DO 20 I = IBEG, IEND
                  IOFFAN = IA + I - 1 + ( JA + KK - 2 ) * LDA
                  CT( I ) = CT( I ) + A( IOFFAK ) * A( IOFFAN )
                  G( I ) = G( I ) + ABS1( A( IOFFAK ) ) *
     $                     ABS1( A( IOFFAN ) )
   20          CONTINUE
   30       CONTINUE
         ELSE IF( TRAN ) THEN
            DO 50 KK = 1, K
               IOFFAK = IA + KK - 1 + ( JA + J - 2 ) * LDA
               DO 40 I = IBEG, IEND
                  IOFFAN = IA + KK - 1 + ( JA + I - 2 ) * LDA
                  CT( I ) = CT( I ) + A( IOFFAK ) * A( IOFFAN )
                  G( I ) = G( I ) + ABS1( A( IOFFAK ) ) *
     $                     ABS1( A( IOFFAN ) )
   40          CONTINUE
   50       CONTINUE
         ELSE IF( HTRAN ) THEN
            DO 70 KK = 1, K
               IOFFAK = IA + J - 1 + ( JA + KK - 2 ) * LDA
               DO 60 I = IBEG, IEND
                  IOFFAN = IA + I - 1 + ( JA + KK - 2 ) * LDA
                  CT( I ) = CT( I ) + A( IOFFAN ) *
     $                      DCONJG( A( IOFFAK ) )
                  G( I ) = G( I ) + ABS1( A( IOFFAK ) ) *
     $                     ABS1( A( IOFFAN ) )
   60          CONTINUE
   70       CONTINUE
         ELSE
            DO 90 KK = 1, K
               IOFFAK = IA + KK - 1 + ( JA + J - 2 ) * LDA
               DO 80 I = IBEG, IEND
                  IOFFAN = IA + KK - 1 + ( JA + I - 2 ) * LDA
                  CT( I ) = CT( I ) + DCONJG( A( IOFFAN ) ) *
     $                      A( IOFFAK )
                  G( I ) = G( I ) + ABS1( DCONJG( A( IOFFAN ) ) ) *
     $                     ABS1( A( IOFFAK ) )
   80          CONTINUE
   90       CONTINUE
         END IF
*
         IOFFC = IC + IBEG - 1 + ( JC + J - 2 ) * LDC
*
         DO 100 I = IBEG, IEND
            CT( I ) = ALPHA*CT( I ) + BETA * C( IOFFC )
            G( I ) = ABS1( ALPHA )*G( I ) +
     $               ABS1( BETA )*ABS1( C( IOFFC ) )
            C( IOFFC ) = CT( I )
            IOFFC = IOFFC + 1
  100    CONTINUE
*
*        Compute the error ratio for this result.
*
         ERR  = RZERO
         INFO = 0
         LDPC = DESCC( LLD_ )
         IOFFC = IC + ( JC + J - 2 ) * LDC
         CALL PB_INFOG2L( IC, JC+J-1, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIC, JJC, ICROW, ICCOL )
         ICURROW = ICROW
         ROWREP  = ( ICROW.EQ.-1 )
         COLREP  = ( ICCOL.EQ.-1 )
*
         IF( MYCOL.EQ.ICCOL .OR. COLREP ) THEN
*
            IBB = DESCC( IMB_ ) - IC + 1
            IF( IBB.LE.0 )
     $         IBB = ( ( -IBB ) / DESCC( MB_ ) + 1 )*DESCC( MB_ ) + IBB
            IBB = MIN( IBB, N )
            IN = IC + IBB - 1
*
            DO 110 I = IC, IN
*
               IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                  ERRI = ABS( PC( IIC+(JJC-1)*LDPC ) -
     $                        C( IOFFC ) ) / EPS
                  IF( G( I-IC+1 ).NE.RZERO )
     $               ERRI = ERRI / G( I-IC+1 )
                  ERR = MAX( ERR, ERRI )
                  IF( ERR*SQRT( EPS ).GE.RONE )
     $               INFO = 1
                  IIC = IIC + 1
               END IF
*
               IOFFC = IOFFC + 1
*
  110       CONTINUE
*
            ICURROW = MOD( ICURROW+1, NPROW )
*
            DO 130 I = IN+1, IC+N-1, DESCC( MB_ )
               IBB = MIN( IC+N-I, DESCC( MB_ ) )
*
               DO 120 KK = 0, IBB-1
*
                  IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                     ERRI = ABS( PC( IIC+(JJC-1)*LDPC ) -
     $                           C( IOFFC ) )/EPS
                     IF( G( I+KK-IC+1 ).NE.RZERO )
     $                  ERRI = ERRI / G( I+KK-IC+1 )
                     ERR = MAX( ERR, ERRI )
                     IF( ERR*SQRT( EPS ).GE.RONE )
     $                  INFO = 1
                     IIC = IIC + 1
                  END IF
*
                  IOFFC = IOFFC + 1
*
  120          CONTINUE
*
               ICURROW = MOD( ICURROW+1, NPROW )
*
  130       CONTINUE
*
         END IF
*
*        If INFO = 0, all results are at least half accurate.
*
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, MYCOL )
         CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERR, 1, I, J, -1, -1,
     $                 MYCOL )
         IF( INFO.NE.0 )
     $      GO TO 150
*
  140 CONTINUE
*
  150 CONTINUE
*
      RETURN
*
*     End of PZMMCH1
*
      END
      SUBROUTINE PZMMCH2( ICTXT, UPLO, TRANS, N, K, ALPHA, A, IA, JA,
     $                    DESCA, B, IB, JB, DESCB, BETA, C, PC, IC,
     $                    JC, DESCC, CT, G, ERR, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        TRANS, UPLO
      INTEGER            IA, IB, IC, ICTXT, INFO, JA, JB, JC, K, N
      DOUBLE PRECISION   ERR
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), DESCC( * )
      DOUBLE PRECISION   G( * )
      COMPLEX*16         A( * ), B( * ), C( * ), CT( * ),
     $                   PC( * )
*     ..
*
*  Purpose
*  =======
*
*  PZMMCH2 checks the results of the computational tests.
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
*  UPLO    (global input) CHARACTER*1
*          On entry,  UPLO  specifies which part of C should contain the
*          result.
*
*  TRANS   (global input) CHARACTER*1
*          On entry,  TRANS  specifies whether the matrices A and B have
*          to  be  transposed  or not before computing the matrix-matrix
*          product.
*
*  N       (global input) INTEGER
*          On entry, N  specifies  the order  the submatrix operand C. N
*          must be at least zero.
*
*  K       (global input) INTEGER
*          On entry, K specifies the number of columns (resp. rows) of A
*          and B when  TRANS = 'N' (resp. TRANS <> 'N').  K  must  be at
*          least zero.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (local input) COMPLEX*16 array
*          On entry, A is an array of  dimension  (DESCA( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PA.
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
*  B       (local input) COMPLEX*16 array
*          On entry, B is an array of  dimension  (DESCB( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PB.
*
*  IB      (global input) INTEGER
*          On entry, IB  specifies B's global row index, which points to
*          the beginning of the submatrix sub( B ).
*
*  JB      (global input) INTEGER
*          On entry, JB  specifies B's global column index, which points
*          to the beginning of the submatrix sub( B ).
*
*  DESCB   (global and local input) INTEGER array
*          On entry, DESCB  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix B.
*
*  BETA    (global input) COMPLEX*16
*          On entry, BETA specifies the scalar beta.
*
*  C       (local input/local output) COMPLEX*16 array
*          On entry, C is an array of  dimension  (DESCC( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PC.
*
*  PC      (local input) COMPLEX*16 array
*          On entry, PC is an array of dimension (DESCC( LLD_ ),*). This
*          array contains the local pieces of the matrix PC.
*
*  IC      (global input) INTEGER
*          On entry, IC  specifies C's global row index, which points to
*          the beginning of the submatrix sub( C ).
*
*  JC      (global input) INTEGER
*          On entry, JC  specifies C's global column index, which points
*          to the beginning of the submatrix sub( C ).
*
*  DESCC   (global and local input) INTEGER array
*          On entry, DESCC  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix C.
*
*  CT      (workspace) COMPLEX*16 array
*          On entry, CT is an array of dimension at least MAX(M,N,K). CT
*          holds a copy of the current column of C.
*
*  G       (workspace) DOUBLE PRECISION array
*          On entry, G  is  an array of dimension at least MAX(M,N,K). G
*          is used to compute the gauges.
*
*  ERR     (global output) DOUBLE PRECISION
*          On exit, ERR specifies the largest error in absolute value.
*
*  INFO    (global output) INTEGER
*          On exit, if INFO <> 0, the result is less than half accurate.
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
      DOUBLE PRECISION   RZERO, RONE
      PARAMETER          ( RZERO = 0.0D+0, RONE = 1.0D+0 )
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, HTRAN, NOTRAN, ROWREP, TRAN, UPPER
      INTEGER            I, IBB, IBEG, ICCOL, ICROW, ICURROW, IEND, IIC,
     $                   IN, IOFFAK, IOFFAN, IOFFBK, IOFFBN, IOFFC, J,
     $                   JJC, KK, LDA, LDB, LDC, LDPC, MYCOL, MYROW,
     $                   NPCOL, NPROW
      DOUBLE PRECISION   EPS, ERRI
      COMPLEX*16         Z
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGAMX2D, IGSUM2D, PB_INFOG2L
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           LSAME, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN, MOD, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
      ABS1( Z ) = ABS( DBLE( Z ) ) + ABS( DIMAG( Z ) )
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      EPS = PDLAMCH( ICTXT, 'eps' )
*
      UPPER = LSAME( UPLO, 'U' )
      HTRAN = LSAME( TRANS, 'H' )
      NOTRAN = LSAME( TRANS, 'N' )
      TRAN = LSAME( TRANS, 'T' )
*
      LDA = MAX( 1, DESCA( M_ ) )
      LDB = MAX( 1, DESCB( M_ ) )
      LDC = MAX( 1, DESCC( M_ ) )
*
*     Compute expected result in C using data in A, B and C.
*     Compute gauges in G. This part of the computation is performed
*     by every process in the grid.
*
      DO 140 J = 1, N
*
         IF( UPPER ) THEN
            IBEG = 1
            IEND = J
         ELSE
            IBEG = J
            IEND = N
         END IF
*
         DO 10 I = 1, N
            CT( I ) = ZERO
            G( I )  = RZERO
   10    CONTINUE
*
         IF( NOTRAN ) THEN
            DO 30 KK = 1, K
               IOFFAK = IA + J - 1 + ( JA + KK - 2 ) * LDA
               IOFFBK = IB + J - 1 + ( JB + KK - 2 ) * LDB
               DO 20 I = IBEG, IEND
                  IOFFAN = IA + I - 1 + ( JA + KK - 2 ) * LDA
                  IOFFBN = IB + I - 1 + ( JB + KK - 2 ) * LDB
                  CT( I ) = CT( I ) + ALPHA * (
     $                      A( IOFFAN ) * B( IOFFBK ) +
     $                      B( IOFFBN ) * A( IOFFAK ) )
                  G( I ) = G( I ) + ABS( ALPHA ) * (
     $                     ABS1( A( IOFFAN ) ) * ABS1( B( IOFFBK ) ) +
     $                     ABS1( B( IOFFBN ) ) * ABS1( A( IOFFAK ) ) )
   20          CONTINUE
   30       CONTINUE
         ELSE IF( TRAN ) THEN
            DO 50 KK = 1, K
               IOFFAK = IA + KK - 1 + ( JA + J - 2 ) * LDA
               IOFFBK = IB + KK - 1 + ( JB + J - 2 ) * LDB
               DO 40 I = IBEG, IEND
                  IOFFAN = IA + KK - 1 + ( JA + I - 2 ) * LDA
                  IOFFBN = IB + KK - 1 + ( JB + I - 2 ) * LDB
                  CT( I ) = CT( I ) + ALPHA * (
     $                      A( IOFFAN ) * B( IOFFBK ) +
     $                      B( IOFFBN ) * A( IOFFAK ) )
                  G( I ) = G( I ) + ABS( ALPHA ) * (
     $                     ABS1( A( IOFFAN ) ) * ABS1( B( IOFFBK ) ) +
     $                     ABS1( B( IOFFBN ) ) * ABS1( A( IOFFAK ) ) )
   40          CONTINUE
   50       CONTINUE
         ELSE IF( HTRAN ) THEN
            DO 70 KK = 1, K
               IOFFAK = IA + J - 1 + ( JA + KK - 2 ) * LDA
               IOFFBK = IB + J - 1 + ( JB + KK - 2 ) * LDB
               DO 60 I = IBEG, IEND
                  IOFFAN = IA + I - 1 + ( JA + KK - 2 ) * LDA
                  IOFFBN = IB + I - 1 + ( JB + KK - 2 ) * LDB
                  CT( I ) = CT( I ) +
     $                ALPHA * A( IOFFAN ) * DCONJG( B( IOFFBK ) ) +
     $                B( IOFFBN ) * DCONJG( ALPHA * A( IOFFAK ) )
                  G( I ) = G( I ) + ABS1( ALPHA ) * (
     $                     ABS1( A( IOFFAN ) ) * ABS1( B( IOFFBK ) ) +
     $                     ABS1( B( IOFFBN ) ) * ABS1( A( IOFFAK ) ) )
   60          CONTINUE
   70       CONTINUE
         ELSE
            DO 90 KK = 1, K
               IOFFAK = IA + KK - 1 + ( JA + J - 2 ) * LDA
               IOFFBK = IB + KK - 1 + ( JB + J - 2 ) * LDB
               DO 80 I = IBEG, IEND
                  IOFFAN = IA + KK - 1 + ( JA + I - 2 ) * LDA
                  IOFFBN = IB + KK - 1 + ( JB + I - 2 ) * LDB
                  CT( I ) = CT( I ) +
     $                   ALPHA * DCONJG( A( IOFFAN ) ) * B( IOFFBK ) +
     $                   DCONJG( ALPHA * B( IOFFBN ) ) * A( IOFFAK )
                  G( I ) = G( I ) + ABS1( ALPHA ) * (
     $                   ABS1( DCONJG( A( IOFFAN ) ) * B( IOFFBK ) ) +
     $                   ABS1( DCONJG( B( IOFFBN ) ) * A( IOFFAK ) ) )
   80          CONTINUE
   90       CONTINUE
         END IF
*
         IOFFC = IC + IBEG - 1 + ( JC + J - 2 ) * LDC
*
         DO 100 I = IBEG, IEND
            CT( I ) = CT( I ) + BETA * C( IOFFC )
            G( I ) = G( I ) + ABS1( BETA )*ABS1( C( IOFFC ) )
            C( IOFFC ) = CT( I )
            IOFFC = IOFFC + 1
  100    CONTINUE
*
*        Compute the error ratio for this result.
*
         ERR  = RZERO
         INFO = 0
         LDPC = DESCC( LLD_ )
         IOFFC = IC + ( JC + J - 2 ) * LDC
         CALL PB_INFOG2L( IC, JC+J-1, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIC, JJC, ICROW, ICCOL )
         ICURROW = ICROW
         ROWREP  = ( ICROW.EQ.-1 )
         COLREP  = ( ICCOL.EQ.-1 )
*
         IF( MYCOL.EQ.ICCOL .OR. COLREP ) THEN
*
            IBB = DESCC( IMB_ ) - IC + 1
            IF( IBB.LE.0 )
     $         IBB = ( ( -IBB ) / DESCC( MB_ ) + 1 )*DESCC( MB_ ) + IBB
            IBB = MIN( IBB, N )
            IN = IC + IBB - 1
*
            DO 110 I = IC, IN
*
               IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                  ERRI = ABS( PC( IIC+(JJC-1)*LDPC ) -
     $                        C( IOFFC ) ) / EPS
                  IF( G( I-IC+1 ).NE.RZERO )
     $               ERRI = ERRI / G( I-IC+1 )
                  ERR = MAX( ERR, ERRI )
                  IF( ERR*SQRT( EPS ).GE.RONE )
     $               INFO = 1
                  IIC = IIC + 1
               END IF
*
               IOFFC = IOFFC + 1
*
  110       CONTINUE
*
            ICURROW = MOD( ICURROW+1, NPROW )
*
            DO 130 I = IN+1, IC+N-1, DESCC( MB_ )
               IBB = MIN( IC+N-I, DESCC( MB_ ) )
*
               DO 120 KK = 0, IBB-1
*
                  IF( MYROW.EQ.ICURROW .OR. ROWREP ) THEN
                     ERRI = ABS( PC( IIC+(JJC-1)*LDPC ) -
     $                           C( IOFFC ) )/EPS
                     IF( G( I+KK-IC+1 ).NE.RZERO )
     $                  ERRI = ERRI / G( I+KK-IC+1 )
                     ERR = MAX( ERR, ERRI )
                     IF( ERR*SQRT( EPS ).GE.RONE )
     $                  INFO = 1
                     IIC = IIC + 1
                  END IF
*
                  IOFFC = IOFFC + 1
*
  120          CONTINUE
*
               ICURROW = MOD( ICURROW+1, NPROW )
*
  130       CONTINUE
*
         END IF
*
*        If INFO = 0, all results are at least half accurate.
*
         CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, MYCOL )
         CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERR, 1, I, J, -1, -1,
     $                 MYCOL )
         IF( INFO.NE.0 )
     $      GO TO 150
*
  140 CONTINUE
*
  150 CONTINUE
*
      RETURN
*
*     End of PZMMCH2
*
      END
      SUBROUTINE PZMMCH3( UPLO, TRANS, M, N, ALPHA, A, IA, JA, DESCA,
     $                    BETA, C, PC, IC, JC, DESCC, ERR, INFO )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        TRANS, UPLO
      INTEGER            IA, IC, INFO, JA, JC, M, N
      DOUBLE PRECISION   ERR
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCC( * )
      COMPLEX*16         A( * ), C( * ), PC( * )
*     ..
*
*  Purpose
*  =======
*
*  PZMMCH3 checks the results of the computational tests.
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
*  UPLO    (global input) CHARACTER*1
*          On entry,  UPLO  specifies which part of C should contain the
*          result.
*
*  TRANS   (global input) CHARACTER*1
*          On entry,  TRANS  specifies  whether  the  matrix A has to be
*          transposed  or not  before computing the  matrix-matrix addi-
*          tion.
*
*  M       (global input) INTEGER
*          On entry, M specifies the number of rows of C.
*
*  N       (global input) INTEGER
*          On entry, N specifies the number of columns of C.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (local input) COMPLEX*16 array
*          On entry, A is an array of  dimension  (DESCA( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PA.
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
*  BETA    (global input) COMPLEX*16
*          On entry, BETA specifies the scalar beta.
*
*  C       (local input/local output) COMPLEX*16 array
*          On entry, C is an array of  dimension  (DESCC( M_ ),*).  This
*          array contains a local copy of the initial entire matrix PC.
*
*  PC      (local input) COMPLEX*16 array
*          On entry, PC is an array of dimension (DESCC( LLD_ ),*). This
*          array contains the local pieces of the matrix PC.
*
*  IC      (global input) INTEGER
*          On entry, IC  specifies C's global row index, which points to
*          the beginning of the submatrix sub( C ).
*
*  JC      (global input) INTEGER
*          On entry, JC  specifies C's global column index, which points
*          to the beginning of the submatrix sub( C ).
*
*  DESCC   (global and local input) INTEGER array
*          On entry, DESCC  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix C.
*
*  ERR     (global output) DOUBLE PRECISION
*          On exit, ERR specifies the largest error in absolute value.
*
*  INFO    (global output) INTEGER
*          On exit, if INFO <> 0, the result is less than half accurate.
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
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, CTRAN, LOWER, NOTRAN, ROWREP, UPPER
      INTEGER            I, ICCOL, ICROW, ICTXT, IIC, IOFFA, IOFFC, J,
     $                   JJC, LDA, LDC, LDPC, MYCOL, MYROW, NPCOL,
     $                   NPROW
      DOUBLE PRECISION   ERR0, ERRI, PREC
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGAMX2D, IGSUM2D, PB_INFOG2L,
     $                   PZERRAXPBY
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           LSAME, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
      ICTXT = DESCC( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      PREC   = PDLAMCH( ICTXT, 'eps' )
*
      UPPER  = LSAME( UPLO,  'U' )
      LOWER  = LSAME( UPLO,  'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      CTRAN  = LSAME( TRANS, 'C' )
*
*     Compute expected result in C using data in A and C. This part of
*     the computation is performed by every process in the grid.
*
      INFO   = 0
      ERR    = ZERO
*
      LDA    = MAX( 1, DESCA( M_   ) )
      LDC    = MAX( 1, DESCC( M_   ) )
      LDPC   = MAX( 1, DESCC( LLD_ ) )
      ROWREP = ( DESCC( RSRC_ ).EQ.-1 )
      COLREP = ( DESCC( CSRC_ ).EQ.-1 )
*
      IF( NOTRAN ) THEN
*
         DO 20 J = JC, JC + N - 1
*
            IOFFC = IC + ( J  - 1          ) * LDC
            IOFFA = IA + ( JA - 1 + J - JC ) * LDA
*
            DO 10 I = IC, IC + M - 1
*
               IF( UPPER ) THEN
                  IF( ( J - JC ).GE.( I - IC ) ) THEN
                     CALL PZERRAXPBY( ERRI, ALPHA, A( IOFFA ), BETA,
     $                                C( IOFFC ), PREC )
                  ELSE
                     ERRI = ZERO
                  END IF
               ELSE IF( LOWER ) THEN
                  IF( ( J - JC ).LE.( I - IC ) ) THEN
                     CALL PZERRAXPBY( ERRI, ALPHA, A( IOFFA ), BETA,
     $                                C( IOFFC ), PREC )
                  ELSE
                     ERRI = ZERO
                  END IF
               ELSE
                  CALL PZERRAXPBY( ERRI, ALPHA, A( IOFFA ), BETA,
     $                             C( IOFFC ), PREC )
               END IF
*
               CALL PB_INFOG2L( I, J, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                          IIC, JJC, ICROW, ICCOL )
               IF( ( MYROW.EQ.ICROW .OR. ROWREP ) .AND.
     $             ( MYCOL.EQ.ICCOL .OR. COLREP ) ) THEN
                  ERR0 = ABS( PC( IIC+(JJC-1)*LDPC )-C( IOFFC ) )
                  IF( ERR0.GT.ERRI )
     $               INFO = 1
                  ERR = MAX( ERR, ERR0 )
               END IF
*
               IOFFA = IOFFA + 1
               IOFFC = IOFFC + 1
*
   10       CONTINUE
*
   20    CONTINUE
*
      ELSE IF( CTRAN ) THEN
*
         DO 40 J = JC, JC + N - 1
*
            IOFFC = IC +              ( J  - 1 ) * LDC
            IOFFA = IA + ( J - JC ) + ( JA - 1 ) * LDA
*
            DO 30 I = IC, IC + M - 1
*
               IF( UPPER ) THEN
                  IF( ( J - JC ).GE.( I - IC ) ) THEN
                     CALL PZERRAXPBY( ERRI, ALPHA, DCONJG( A( IOFFA ) ),
     $                                BETA, C( IOFFC ), PREC )
                  ELSE
                     ERRI = ZERO
                  END IF
               ELSE IF( LOWER ) THEN
                  IF( ( J - JC ).LE.( I - IC ) ) THEN
                     CALL PZERRAXPBY( ERRI, ALPHA, DCONJG( A( IOFFA ) ),
     $                                BETA, C( IOFFC ), PREC )
                  ELSE
                     ERRI = ZERO
                  END IF
               ELSE
                  CALL PZERRAXPBY( ERRI, ALPHA, DCONJG( A( IOFFA ) ),
     $                             BETA, C( IOFFC ), PREC )
               END IF
*
               CALL PB_INFOG2L( I, J, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                          IIC, JJC, ICROW, ICCOL )
               IF( ( MYROW.EQ.ICROW .OR. ROWREP ) .AND.
     $             ( MYCOL.EQ.ICCOL .OR. COLREP ) ) THEN
                  ERR0 = ABS( PC( IIC+(JJC-1)*LDPC )-C( IOFFC ) )
                  IF( ERR0.GT.ERRI )
     $               INFO = 1
                  ERR = MAX( ERR, ERR0 )
               END IF
*
               IOFFC = IOFFC + 1
               IOFFA = IOFFA + LDA
*
   30       CONTINUE
*
   40    CONTINUE
*
      ELSE
*
         DO 60 J = JC, JC + N - 1
*
            IOFFC = IC +              ( J  - 1 ) * LDC
            IOFFA = IA + ( J - JC ) + ( JA - 1 ) * LDA
*
            DO 50 I = IC, IC + M - 1
*
               IF( UPPER ) THEN
                  IF( ( J - JC ).GE.( I - IC ) ) THEN
                     CALL PZERRAXPBY( ERRI, ALPHA, A( IOFFA ), BETA,
     $                                C( IOFFC ), PREC )
                  ELSE
                     ERRI = ZERO
                  END IF
               ELSE IF( LOWER ) THEN
                  IF( ( J - JC ).LE.( I - IC ) ) THEN
                     CALL PZERRAXPBY( ERRI, ALPHA, A( IOFFA ), BETA,
     $                                C( IOFFC ), PREC )
                  ELSE
                     ERRI = ZERO
                  END IF
               ELSE
                  CALL PZERRAXPBY( ERRI, ALPHA, A( IOFFA ), BETA,
     $                             C( IOFFC ), PREC )
               END IF
*
               CALL PB_INFOG2L( I, J, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $                          IIC, JJC, ICROW, ICCOL )
               IF( ( MYROW.EQ.ICROW .OR. ROWREP ) .AND.
     $             ( MYCOL.EQ.ICCOL .OR. COLREP ) ) THEN
                  ERR0 = ABS( PC( IIC+(JJC-1)*LDPC )-C( IOFFC ) )
                  IF( ERR0.GT.ERRI )
     $               INFO = 1
                  ERR = MAX( ERR, ERR0 )
               END IF
*
               IOFFC = IOFFC + 1
               IOFFA = IOFFA + LDA
*
   50       CONTINUE
*
   60    CONTINUE
*
      END IF
*
*     If INFO = 0, all results are at least half accurate.
*
      CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, MYCOL )
      CALL DGAMX2D( ICTXT, 'All', ' ', 1, 1, ERR, 1, I, J, -1, -1,
     $              MYCOL )
*
      RETURN
*
*     End of PZMMCH3
*
      END
      SUBROUTINE PZERRAXPBY( ERRBND, ALPHA, X, BETA, Y, PREC )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ERRBND, PREC
      COMPLEX*16         ALPHA, BETA, X, Y
*     ..
*
*  Purpose
*  =======
*
*  PZERRAXPBY  serially  computes  y := beta*y + alpha * x and returns a
*  scaled relative acceptable error bound on the result.
*
*  Arguments
*  =========
*
*  ERRBND  (global output) DOUBLE PRECISION
*          On exit, ERRBND  specifies the scaled relative acceptable er-
*          ror bound.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  X       (global input) COMPLEX*16
*          On entry, X  specifies the scalar x to be scaled.
*
*  BETA    (global input) COMPLEX*16
*          On entry, BETA specifies the scalar beta.
*
*  Y       (global input/global output) COMPLEX*16
*          On entry,  Y  specifies  the scalar y to be added. On exit, Y
*          contains the resulting scalar y.
*
*  PREC    (global input) DOUBLE PRECISION
*          On entry, PREC specifies the machine precision.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO, ZERO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0,
     $                   ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   ADDBND, FACT, SUMINEG, SUMIPOS, SUMRNEG,
     $                   SUMRPOS
      COMPLEX*16         TMP
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. Executable Statements ..
*
      SUMIPOS = ZERO
      SUMINEG = ZERO
      SUMRPOS = ZERO
      SUMRNEG = ZERO
      FACT = ONE + TWO * PREC
      ADDBND = TWO * TWO * TWO * PREC
*
      TMP = ALPHA * X
      IF( DBLE( TMP ).GE.ZERO ) THEN
         SUMRPOS = SUMRPOS + DBLE( TMP ) * FACT
      ELSE
         SUMRNEG = SUMRNEG - DBLE( TMP ) * FACT
      END IF
      IF( DIMAG( TMP ).GE.ZERO ) THEN
         SUMIPOS = SUMIPOS + DIMAG( TMP ) * FACT
      ELSE
         SUMINEG = SUMINEG - DIMAG( TMP ) * FACT
      END IF
*
      TMP = BETA * Y
      IF( DBLE( TMP ).GE.ZERO ) THEN
         SUMRPOS = SUMRPOS + DBLE( TMP ) * FACT
      ELSE
         SUMRNEG = SUMRNEG - DBLE( TMP ) * FACT
      END IF
      IF( DIMAG( TMP ).GE.ZERO ) THEN
         SUMIPOS = SUMIPOS + DIMAG( TMP ) * FACT
      ELSE
         SUMINEG = SUMINEG - DIMAG( TMP ) * FACT
      END IF
*
      Y = ( BETA * Y ) + ( ALPHA * X )
*
      ERRBND = ADDBND * MAX( MAX( SUMRPOS, SUMRNEG ),
     $                       MAX( SUMIPOS, SUMINEG ) )
*
      RETURN
*
*     End of PZERRAXPBY
*
      END
      SUBROUTINE PZIPSET( TOGGLE, N, A, IA, JA, DESCA )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        TOGGLE
      INTEGER            IA, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         A( * )
*     ..
*
*  Purpose
*  =======
*
*  PZIPSET  sets the imaginary part of the diagonal entries of an n by n
*  matrix sub( A )  denoting  A( IA:IA+N-1, JA:JA+N-1 ). This is used to
*  test the  PBLAS routines  for  complex Hermitian  matrices, which are
*  either  not supposed to access or use the imaginary parts of the dia-
*  gonals, or supposed to set  them to zero. The  value  used to set the
*  imaginary part of the diagonals depends on the value of TOGGLE.
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
*  TOGGLE  (global input) CHARACTER*1
*          On entry,  TOGGLE  specifies the set-value to be used as fol-
*          lows:
*             If TOGGLE = 'Z' or 'z', the  imaginary  part of the diago-
*                                     nals are set to zero,
*             If TOGGLE = 'B' or 'b', the  imaginary  part of the diago-
*                                     nals are set to a large value.
*
*  N       (global input) INTEGER
*          On entry,  N  specifies  the  order of sub( A ). N must be at
*          least zero.
*
*  A       (local input/local output) pointer to COMPLEX*16
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix  A. On exit, the diagonals of
*          sub( A ) have been updated as specified by TOGGLE.
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
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLREP, GODOWN, GOLEFT, ROWREP
      INTEGER            I, IACOL, IAROW, ICTXT, IIA, IJOFFA, ILOW,
     $                   IMB1, IMBLOC, INB1, INBLOC, IOFFA, IOFFD, IUPP,
     $                   JJA, JOFFA, JOFFD, LCMT, LCMT00, LDA, LDAP1,
     $                   LMBLOC, LNBLOC, LOW, MB, MBLKD, MBLKS, MBLOC,
     $                   MRCOL, MRROW, MYCOL, MYROW, NB, NBLKD, NBLKS,
     $                   NBLOC, NP, NPCOL, NPROW, NQ, PMB, QNB, UPP
      DOUBLE PRECISION   ALPHA, ATMP
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA2( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PB_AINFOG2L, PB_BINFO,
     $                   PB_DESCTRANS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           LSAME, PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Convert descriptor
*
      CALL PB_DESCTRANS( DESCA, DESCA2 )
*
*     Get grid parameters
*
      ICTXT = DESCA2( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( LSAME( TOGGLE, 'Z' ) ) THEN
         ALPHA = ZERO
      ELSE IF( LSAME( TOGGLE, 'B' ) ) THEN
         ALPHA = PDLAMCH( ICTXT, 'Epsilon' )
         ALPHA = ALPHA / PDLAMCH( ICTXT, 'Safe minimum' )
      END IF
*
      CALL PB_AINFOG2L( N, N, IA, JA, DESCA2, NPROW, NPCOL, MYROW,
     $                  MYCOL, IMB1, INB1, NP, NQ, IIA, JJA, IAROW,
     $                  IACOL, MRROW, MRCOL )
*
      IF( NP.LE.0 .OR. NQ.LE.0 )
     $   RETURN
*
*     Initialize LCMT00, MBLKS, NBLKS, IMBLOC, INBLOC, LMBLOC, LNBLOC,
*     ILOW, LOW, IUPP, and UPP.
*
      MB = DESCA2( MB_ )
      NB = DESCA2( NB_ )
      CALL PB_BINFO( 0, NP, NQ, IMB1, INB1, MB, NB, MRROW, MRCOL,
     $               LCMT00, MBLKS, NBLKS, IMBLOC, INBLOC, LMBLOC,
     $               LNBLOC, ILOW, LOW, IUPP, UPP )
*
      IOFFA  = IIA - 1
      JOFFA  = JJA - 1
      ROWREP = ( DESCA2( RSRC_ ).EQ.-1 )
      COLREP = ( DESCA2( CSRC_ ).EQ.-1 )
      LDA    = DESCA2( LLD_ )
      LDAP1  = LDA + 1
*
      IF( ROWREP ) THEN
         PMB = MB
      ELSE
         PMB = NPROW * MB
      END IF
      IF( COLREP ) THEN
         QNB = NB
      ELSE
         QNB = NPCOL * NB
      END IF
*
*     Handle the first block of rows or columns separately, and update
*     LCMT00, MBLKS and NBLKS.
*
      GODOWN = ( LCMT00.GT.IUPP )
      GOLEFT = ( LCMT00.LT.ILOW )
*
      IF( .NOT.GODOWN .AND. .NOT.GOLEFT ) THEN
*
*        LCMT00 >= ILOW && LCMT00 <= IUPP
*
         IF( LCMT00.GE.0 ) THEN
            IJOFFA = IOFFA + LCMT00 + ( JOFFA - 1 ) * LDA
            DO 10 I = 1, MIN( INBLOC, MAX( 0, IMBLOC - LCMT00 ) )
               ATMP = DBLE( A( IJOFFA + I*LDAP1 ) )
               A( IJOFFA + I*LDAP1 ) = DCMPLX( ATMP, ALPHA )
   10       CONTINUE
         ELSE
            IJOFFA = IOFFA + ( JOFFA - LCMT00 - 1 ) * LDA
            DO 20 I = 1, MIN( IMBLOC, MAX( 0, INBLOC + LCMT00 ) )
               ATMP = DBLE( A( IJOFFA + I*LDAP1 ) )
               A( IJOFFA + I*LDAP1 ) = DCMPLX( ATMP, ALPHA )
   20       CONTINUE
         END IF
         GOLEFT = ( ( LCMT00 - ( IUPP - UPP + PMB ) ).LT.ILOW )
         GODOWN = .NOT.GOLEFT
*
      END IF
*
      IF( GODOWN ) THEN
*
         LCMT00 = LCMT00 - ( IUPP - UPP + PMB )
         MBLKS  = MBLKS - 1
         IOFFA  = IOFFA + IMBLOC
*
   30    CONTINUE
         IF( MBLKS.GT.0 .AND. LCMT00.GT.UPP ) THEN
            LCMT00 = LCMT00 - PMB
            MBLKS  = MBLKS - 1
            IOFFA  = IOFFA + MB
            GO TO 30
         END IF
*
         IF( MBLKS.LE.0 )
     $      RETURN
*
         LCMT  = LCMT00
         MBLKD = MBLKS
         IOFFD = IOFFA
*
         MBLOC = MB
   40    CONTINUE
         IF( MBLKD.GT.0 .AND. LCMT.GE.ILOW ) THEN
            IF( MBLKD.EQ.1 )
     $         MBLOC = LMBLOC
            IF( LCMT.GE.0 ) THEN
               IJOFFA = IOFFD + LCMT + ( JOFFA - 1 ) * LDA
               DO 50 I = 1, MIN( INBLOC, MAX( 0, MBLOC - LCMT ) )
                  ATMP = DBLE( A( IJOFFA + I*LDAP1 ) )
                  A( IJOFFA + I*LDAP1 ) = DCMPLX( ATMP, ALPHA )
   50          CONTINUE
            ELSE
               IJOFFA = IOFFD + ( JOFFA - LCMT - 1 ) * LDA
               DO 60 I = 1, MIN( MBLOC, MAX( 0, INBLOC + LCMT ) )
                  ATMP = DBLE( A( IJOFFA + I*LDAP1 ) )
                  A( IJOFFA + I*LDAP1 ) = DCMPLX( ATMP, ALPHA )
   60          CONTINUE
            END IF
            LCMT00 = LCMT
            LCMT   = LCMT - PMB
            MBLKS  = MBLKD
            MBLKD  = MBLKD - 1
            IOFFA  = IOFFD
            IOFFD  = IOFFD + MBLOC
            GO TO 40
         END IF
*
         LCMT00 = LCMT00 + LOW - ILOW + QNB
         NBLKS  = NBLKS - 1
         JOFFA  = JOFFA + INBLOC
*
      ELSE IF( GOLEFT ) THEN
*
         LCMT00 = LCMT00 + LOW - ILOW + QNB
         NBLKS  = NBLKS - 1
         JOFFA  = JOFFA + INBLOC
*
   70    CONTINUE
         IF( NBLKS.GT.0 .AND. LCMT00.LT.LOW ) THEN
            LCMT00 = LCMT00 + QNB
            NBLKS  = NBLKS - 1
            JOFFA  = JOFFA + NB
            GO TO 70
         END IF
*
         IF( NBLKS.LE.0 )
     $      RETURN
*
         LCMT  = LCMT00
         NBLKD = NBLKS
         JOFFD = JOFFA
*
         NBLOC = NB
   80    CONTINUE
         IF( NBLKD.GT.0 .AND. LCMT.LE.IUPP ) THEN
            IF( NBLKD.EQ.1 )
     $         NBLOC = LNBLOC
            IF( LCMT.GE.0 ) THEN
               IJOFFA = IOFFA + LCMT + ( JOFFD - 1 ) * LDA
               DO 90 I = 1, MIN( NBLOC, MAX( 0, IMBLOC - LCMT ) )
                  ATMP = DBLE( A( IJOFFA + I*LDAP1 ) )
                  A( IJOFFA + I*LDAP1 ) = DCMPLX( ATMP, ALPHA )
   90          CONTINUE
            ELSE
               IJOFFA = IOFFA + ( JOFFD - LCMT - 1 ) * LDA
               DO 100 I = 1, MIN( IMBLOC, MAX( 0, NBLOC + LCMT ) )
                  ATMP = DBLE( A( IJOFFA + I*LDAP1 ) )
                  A( IJOFFA + I*LDAP1 ) = DCMPLX( ATMP, ALPHA )
  100          CONTINUE
            END IF
            LCMT00 = LCMT
            LCMT   = LCMT + QNB
            NBLKS  = NBLKD
            NBLKD  = NBLKD - 1
            JOFFA  = JOFFD
            JOFFD  = JOFFD + NBLOC
            GO TO 80
         END IF
*
         LCMT00 = LCMT00 - ( IUPP - UPP + PMB )
         MBLKS  = MBLKS - 1
         IOFFA  = IOFFA + IMBLOC
*
      END IF
*
      NBLOC = NB
  110 CONTINUE
      IF( NBLKS.GT.0 ) THEN
         IF( NBLKS.EQ.1 )
     $      NBLOC = LNBLOC
  120    CONTINUE
         IF( MBLKS.GT.0 .AND. LCMT00.GT.UPP ) THEN
            LCMT00 = LCMT00 - PMB
            MBLKS  = MBLKS - 1
            IOFFA  = IOFFA + MB
            GO TO 120
         END IF
*
         IF( MBLKS.LE.0 )
     $      RETURN
*
         LCMT  = LCMT00
         MBLKD = MBLKS
         IOFFD = IOFFA
*
         MBLOC = MB
  130    CONTINUE
         IF( MBLKD.GT.0 .AND. LCMT.GE.LOW ) THEN
            IF( MBLKD.EQ.1 )
     $         MBLOC = LMBLOC
            IF( LCMT.GE.0 ) THEN
               IJOFFA = IOFFD + LCMT + ( JOFFA - 1 ) * LDA
               DO 140 I = 1, MIN( NBLOC, MAX( 0, MBLOC - LCMT ) )
                  ATMP = DBLE( A( IJOFFA + I*LDAP1 ) )
                  A( IJOFFA + I*LDAP1 ) = DCMPLX( ATMP, ALPHA )
  140          CONTINUE
            ELSE
               IJOFFA = IOFFD + ( JOFFA - LCMT - 1 ) * LDA
               DO 150 I = 1, MIN( MBLOC, MAX( 0, NBLOC + LCMT ) )
                  ATMP = DBLE( A( IJOFFA + I*LDAP1 ) )
                  A( IJOFFA + I*LDAP1 ) = DCMPLX( ATMP, ALPHA )
  150          CONTINUE
            END IF
            LCMT00 = LCMT
            LCMT   = LCMT - PMB
            MBLKS  = MBLKD
            MBLKD  = MBLKD - 1
            IOFFA  = IOFFD
            IOFFD  = IOFFD + MBLOC
            GO TO 130
         END IF
*
         LCMT00 = LCMT00 + QNB
         NBLKS  = NBLKS - 1
         JOFFA  = JOFFA + NBLOC
         GO TO 110
*
      END IF
*
      RETURN
*
*     End of PZIPSET
*
      END
      DOUBLE PRECISION   FUNCTION PDLAMCH( ICTXT, CMACH )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        CMACH
      INTEGER            ICTXT
*     ..
*
*  Purpose
*  =======
*
*
*     .. Local Scalars ..
      CHARACTER*1        TOP
      INTEGER            IDUMM
      DOUBLE PRECISION   TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGAMN2D, DGAMX2D, PB_TOPGET
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, LSAME
*     ..
*     .. Executable Statements ..
*
      TEMP = DLAMCH( CMACH )
*
      IF( LSAME( CMACH, 'E' ).OR.LSAME( CMACH, 'S' ).OR.
     $    LSAME( CMACH, 'M' ).OR.LSAME( CMACH, 'U' ) ) THEN
         CALL PB_TOPGET( ICTXT, 'Combine', 'All', TOP )
         IDUMM = 0
         CALL DGAMX2D( ICTXT, 'All', TOP, 1, 1, TEMP, 1, IDUMM,
     $                 IDUMM, -1, -1, IDUMM )
      ELSE IF( LSAME( CMACH, 'L' ).OR.LSAME( CMACH, 'O' ) ) THEN
         CALL PB_TOPGET( ICTXT, 'Combine', 'All', TOP )
         IDUMM = 0
         CALL DGAMN2D( ICTXT, 'All', TOP, 1, 1, TEMP, 1, IDUMM,
     $                 IDUMM, -1, -1, IDUMM )
      END IF
*
      PDLAMCH = TEMP
*
      RETURN
*
*     End of PDLAMCH
*
      END
      SUBROUTINE PZLASET( UPLO, M, N, ALPHA, BETA, A, IA, JA, DESCA )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO
      INTEGER            IA, JA, M, N
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         A( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLASET  initializes an m by n submatrix A(IA:IA+M-1,JA:JA+N-1) deno-
*  ted  by  sub( A )  to beta on the diagonal and alpha on the offdiago-
*  nals.
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
*  UPLO    (global input) CHARACTER*1
*          On entry, UPLO specifies the part  of  the submatrix sub( A )
*          to be set:
*             = 'L' or 'l':   Lower triangular part is set; the strictly
*                      upper triangular part of sub( A ) is not changed;
*             = 'U' or 'u':   Upper triangular part is set; the strictly
*                      lower triangular part of sub( A ) is not changed;
*             Otherwise:  All of the matrix sub( A ) is set.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( A ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( A ). N must be at least zero.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry,  ALPHA  specifies the scalar alpha, i.e., the cons-
*          tant to which the offdiagonal elements are to be set.
*
*  BETA    (global input) COMPLEX*16
*          On entry, BETA  specifies the scalar beta, i.e., the constant
*          to which the diagonal elements are to be set.
*
*  A       (local input/local output) COMPLEX*16 array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix  A  to be  set.  On exit, the
*          leading m by n submatrix sub( A ) is set as follows:
*
*          if UPLO = 'U',  A(IA+i-1,JA+j-1) = ALPHA, 1<=i<=j-1, 1<=j<=N,
*          if UPLO = 'L',  A(IA+i-1,JA+j-1) = ALPHA, j+1<=i<=M, 1<=j<=N,
*          otherwise,      A(IA+i-1,JA+j-1) = ALPHA, 1<=i<=M,   1<=j<=N,
*                                                      and IA+i.NE.JA+j,
*          and, for all UPLO,  A(IA+i-1,JA+i-1) = BETA,  1<=i<=min(M,N).
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
      LOGICAL            GODOWN, GOLEFT, ISCOLREP, ISROWREP, LOWER,
     $                   UPPER
      INTEGER            IACOL, IAROW, ICTXT, IIA, IIMAX, ILOW, IMB1,
     $                   IMBLOC, INB1, INBLOC, IOFFA, IOFFD, IUPP, JJA,
     $                   JJMAX, JOFFA, JOFFD, LCMT, LCMT00, LDA, LMBLOC,
     $                   LNBLOC, LOW, M1, MB, MBLKD, MBLKS, MBLOC, MP,
     $                   MRCOL, MRROW, MYCOL, MYROW, N1, NB, NBLKD,
     $                   NBLKS, NBLOC, NPCOL, NPROW, NQ, PMB, QNB, TMP1,
     $                   UPP
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA2( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PB_AINFOG2L, PB_BINFO,
     $                   PB_DESCTRANS, PB_ZLASET
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Convert descriptor
*
      CALL PB_DESCTRANS( DESCA, DESCA2 )
*
*     Get grid parameters
*
      ICTXT = DESCA2( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      CALL PB_AINFOG2L( M, N, IA, JA, DESCA2, NPROW, NPCOL, MYROW,
     $                  MYCOL, IMB1, INB1, MP, NQ, IIA, JJA, IAROW,
     $                  IACOL, MRROW, MRCOL )
*
      IF( MP.LE.0 .OR. NQ.LE.0 )
     $   RETURN
*
      ISROWREP = ( DESCA2( RSRC_ ).LT.0 )
      ISCOLREP = ( DESCA2( CSRC_ ).LT.0 )
      LDA      = DESCA2( LLD_ )
*
      UPPER = .NOT.( LSAME( UPLO, 'L' ) )
      LOWER = .NOT.( LSAME( UPLO, 'U' ) )
*
      IF( ( ( LOWER.AND.UPPER ).AND.( ALPHA.EQ.BETA ) ).OR.
     $    (   ISROWREP        .AND.  ISCOLREP        ) ) THEN
         IF( ( MP.GT.0 ).AND.( NQ.GT.0 ) )
     $      CALL PB_ZLASET( UPLO, MP, NQ, 0, ALPHA, BETA,
     $                      A( IIA + ( JJA - 1 ) * LDA ), LDA )
         RETURN
      END IF
*
*     Initialize LCMT00, MBLKS, NBLKS, IMBLOC, INBLOC, LMBLOC, LNBLOC,
*     ILOW, LOW, IUPP, and UPP.
*
      MB = DESCA2( MB_ )
      NB = DESCA2( NB_ )
      CALL PB_BINFO( 0, MP, NQ, IMB1, INB1, MB, NB, MRROW, MRCOL,
     $               LCMT00, MBLKS, NBLKS, IMBLOC, INBLOC, LMBLOC,
     $               LNBLOC, ILOW, LOW, IUPP, UPP )
*
      IOFFA = IIA - 1
      JOFFA = JJA - 1
      IIMAX = IOFFA + MP
      JJMAX = JOFFA + NQ
*
      IF( ISROWREP ) THEN
         PMB = MB
      ELSE
         PMB = NPROW * MB
      END IF
      IF( ISCOLREP ) THEN
         QNB = NB
      ELSE
         QNB = NPCOL * NB
      END IF
*
      M1 = MP
      N1 = NQ
*
*     Handle the first block of rows or columns separately, and update
*     LCMT00, MBLKS and NBLKS.
*
      GODOWN = ( LCMT00.GT.IUPP )
      GOLEFT = ( LCMT00.LT.ILOW )
*
      IF( .NOT.GODOWN .AND. .NOT.GOLEFT ) THEN
*
*        LCMT00 >= ILOW && LCMT00 <= IUPP
*
         GOLEFT = ( ( LCMT00 - ( IUPP - UPP + PMB ) ).LT.ILOW )
         GODOWN = .NOT.GOLEFT
*
         CALL PB_ZLASET( UPLO, IMBLOC, INBLOC, LCMT00, ALPHA, BETA,
     $                   A( IIA+JOFFA*LDA ), LDA )
         IF( GODOWN ) THEN
            IF( UPPER .AND. NQ.GT.INBLOC )
     $         CALL PB_ZLASET( 'All', IMBLOC, NQ-INBLOC, 0, ALPHA,
     $                         ALPHA, A( IIA+(JOFFA+INBLOC)*LDA ), LDA )
            IIA = IIA + IMBLOC
            M1  = M1 - IMBLOC
         ELSE
            IF( LOWER .AND. MP.GT.IMBLOC )
     $         CALL PB_ZLASET( 'All', MP-IMBLOC, INBLOC, 0, ALPHA,
     $                         ALPHA, A( IIA+IMBLOC+JOFFA*LDA ), LDA )
            JJA = JJA + INBLOC
            N1  = N1 - INBLOC
         END IF
*
      END IF
*
      IF( GODOWN ) THEN
*
         LCMT00 = LCMT00 - ( IUPP - UPP + PMB )
         MBLKS  = MBLKS - 1
         IOFFA  = IOFFA + IMBLOC
*
   10    CONTINUE
         IF( MBLKS.GT.0 .AND. LCMT00.GT.UPP ) THEN
            LCMT00 = LCMT00 - PMB
            MBLKS  = MBLKS - 1
            IOFFA  = IOFFA + MB
            GO TO 10
         END IF
*
         TMP1 = MIN( IOFFA, IIMAX ) - IIA + 1
         IF( UPPER .AND. TMP1.GT.0 ) THEN
            CALL PB_ZLASET( 'All', TMP1, N1, 0, ALPHA, ALPHA,
     $                      A( IIA+JOFFA*LDA ), LDA )
            IIA = IIA + TMP1
            M1  = M1 - TMP1
         END IF
*
         IF( MBLKS.LE.0 )
     $      RETURN
*
         LCMT  = LCMT00
         MBLKD = MBLKS
         IOFFD = IOFFA
*
         MBLOC = MB
   20    CONTINUE
         IF( MBLKD.GT.0 .AND. LCMT.GE.ILOW ) THEN
            IF( MBLKD.EQ.1 )
     $         MBLOC = LMBLOC
            CALL PB_ZLASET( UPLO, MBLOC, INBLOC, LCMT, ALPHA, BETA,
     $                      A( IOFFD+1+JOFFA*LDA ), LDA )
            LCMT00 = LCMT
            LCMT   = LCMT - PMB
            MBLKS  = MBLKD
            MBLKD  = MBLKD - 1
            IOFFA  = IOFFD
            IOFFD  = IOFFD + MBLOC
            GO TO 20
         END IF
*
         TMP1 = M1 - IOFFD + IIA - 1
         IF( LOWER .AND. TMP1.GT.0 )
     $      CALL PB_ZLASET( 'ALL', TMP1, INBLOC, 0, ALPHA, ALPHA,
     $                      A( IOFFD+1+JOFFA*LDA ), LDA )
*
         TMP1   = IOFFA - IIA + 1
         M1     = M1 - TMP1
         N1     = N1 - INBLOC
         LCMT00 = LCMT00 + LOW - ILOW + QNB
         NBLKS  = NBLKS - 1
         JOFFA  = JOFFA + INBLOC
*
         IF( UPPER .AND. TMP1.GT.0 .AND. N1.GT.0 )
     $      CALL PB_ZLASET( 'ALL', TMP1, N1, 0, ALPHA, ALPHA,
     $                      A( IIA+JOFFA*LDA ), LDA )
*
         IIA = IOFFA + 1
         JJA = JOFFA + 1
*
      ELSE IF( GOLEFT ) THEN
*
         LCMT00 = LCMT00 + LOW - ILOW + QNB
         NBLKS  = NBLKS - 1
         JOFFA  = JOFFA + INBLOC
*
   30    CONTINUE
         IF( NBLKS.GT.0 .AND. LCMT00.LT.LOW ) THEN
            LCMT00 = LCMT00 + QNB
            NBLKS  = NBLKS - 1
            JOFFA  = JOFFA + NB
            GO TO 30
         END IF
*
         TMP1 = MIN( JOFFA, JJMAX ) - JJA + 1
         IF( LOWER .AND. TMP1.GT.0 ) THEN
            CALL PB_ZLASET( 'All', M1, TMP1, 0, ALPHA, ALPHA,
     $                      A( IIA+(JJA-1)*LDA ), LDA )
            JJA = JJA + TMP1
            N1  = N1 - TMP1
         END IF
*
         IF( NBLKS.LE.0 )
     $      RETURN
*
         LCMT  = LCMT00
         NBLKD = NBLKS
         JOFFD = JOFFA
*
         NBLOC = NB
   40    CONTINUE
         IF( NBLKD.GT.0 .AND. LCMT.LE.IUPP ) THEN
            IF( NBLKD.EQ.1 )
     $         NBLOC = LNBLOC
            CALL PB_ZLASET( UPLO, IMBLOC, NBLOC, LCMT, ALPHA, BETA,
     $                      A( IIA+JOFFD*LDA ), LDA )
            LCMT00 = LCMT
            LCMT   = LCMT + QNB
            NBLKS  = NBLKD
            NBLKD  = NBLKD - 1
            JOFFA  = JOFFD
            JOFFD  = JOFFD + NBLOC
            GO TO 40
         END IF
*
         TMP1 = N1 - JOFFD + JJA - 1
         IF( UPPER .AND. TMP1.GT.0 )
     $      CALL PB_ZLASET( 'All', IMBLOC, TMP1, 0, ALPHA, ALPHA,
     $                      A( IIA+JOFFD*LDA ), LDA )
*
         TMP1   = JOFFA - JJA + 1
         M1     = M1 - IMBLOC
         N1     = N1 - TMP1
         LCMT00 = LCMT00 - ( IUPP - UPP + PMB )
         MBLKS  = MBLKS - 1
         IOFFA  = IOFFA + IMBLOC
*
         IF( LOWER .AND. M1.GT.0 .AND. TMP1.GT.0 )
     $      CALL PB_ZLASET( 'All', M1, TMP1, 0, ALPHA, ALPHA,
     $                      A( IOFFA+1+(JJA-1)*LDA ), LDA )
*
         IIA = IOFFA + 1
         JJA = JOFFA + 1
*
      END IF
*
      NBLOC = NB
   50 CONTINUE
      IF( NBLKS.GT.0 ) THEN
         IF( NBLKS.EQ.1 )
     $      NBLOC = LNBLOC
   60    CONTINUE
         IF( MBLKS.GT.0 .AND. LCMT00.GT.UPP ) THEN
            LCMT00 = LCMT00 - PMB
            MBLKS  = MBLKS - 1
            IOFFA  = IOFFA + MB
            GO TO 60
         END IF
*
         TMP1 = MIN( IOFFA, IIMAX ) - IIA + 1
         IF( UPPER .AND. TMP1.GT.0 ) THEN
            CALL PB_ZLASET( 'All', TMP1, N1, 0, ALPHA, ALPHA,
     $                      A( IIA+JOFFA*LDA ), LDA )
            IIA = IIA + TMP1
            M1  = M1 - TMP1
         END IF
*
         IF( MBLKS.LE.0 )
     $      RETURN
*
         LCMT  = LCMT00
         MBLKD = MBLKS
         IOFFD = IOFFA
*
         MBLOC = MB
   70    CONTINUE
         IF( MBLKD.GT.0 .AND. LCMT.GE.LOW ) THEN
            IF( MBLKD.EQ.1 )
     $         MBLOC = LMBLOC
            CALL PB_ZLASET( UPLO, MBLOC, NBLOC, LCMT, ALPHA, BETA,
     $                      A( IOFFD+1+JOFFA*LDA ), LDA )
            LCMT00 = LCMT
            LCMT   = LCMT - PMB
            MBLKS  = MBLKD
            MBLKD  = MBLKD - 1
            IOFFA  = IOFFD
            IOFFD  = IOFFD + MBLOC
            GO TO 70
         END IF
*
         TMP1 = M1 - IOFFD + IIA - 1
         IF( LOWER .AND. TMP1.GT.0 )
     $      CALL PB_ZLASET( 'All', TMP1, NBLOC, 0, ALPHA, ALPHA,
     $                      A( IOFFD+1+JOFFA*LDA ), LDA )
*
         TMP1   = MIN( IOFFA, IIMAX )  - IIA + 1
         M1     = M1 - TMP1
         N1     = N1 - NBLOC
         LCMT00 = LCMT00 + QNB
         NBLKS  = NBLKS - 1
         JOFFA  = JOFFA + NBLOC
*
         IF( UPPER .AND. TMP1.GT.0 .AND. N1.GT.0 )
     $      CALL PB_ZLASET( 'All', TMP1, N1, 0, ALPHA, ALPHA,
     $                      A( IIA+JOFFA*LDA ), LDA )
*
         IIA = IOFFA + 1
         JJA = JOFFA + 1
*
         GO TO 50
*
      END IF
*
      RETURN
*
*     End of PZLASET
*
      END
      SUBROUTINE PZLASCAL( TYPE, M, N, ALPHA, A, IA, JA, DESCA )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        TYPE
      INTEGER            IA, JA, M, N
      COMPLEX*16         ALPHA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         A( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLASCAL  scales the  m by n submatrix A(IA:IA+M-1,JA:JA+N-1) denoted
*  by sub( A ) by the scalar alpha. TYPE  specifies if sub( A ) is full,
*  upper triangular, lower triangular or upper Hessenberg.
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
*  TYPE    (global input) CHARACTER*1
*          On entry,  TYPE  specifies the type of the input submatrix as
*          follows:
*             = 'L' or 'l':  sub( A ) is a lower triangular matrix,
*             = 'U' or 'u':  sub( A ) is an upper triangular matrix,
*             = 'H' or 'h':  sub( A ) is an upper Hessenberg matrix,
*             otherwise sub( A ) is a  full matrix.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( A ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( A ). N  must be at least zero.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (local input/local output) COMPLEX*16 array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix  A.
*          On exit, the local entries of this array corresponding to the
*          to  the entries of the submatrix sub( A ) are  overwritten by
*          the local entries of the m by n scaled submatrix.
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
      CHARACTER*1        UPLO
      LOGICAL            GODOWN, GOLEFT, LOWER, UPPER
      INTEGER            IACOL, IAROW, ICTXT, IIA, IIMAX, ILOW, IMB1,
     $                   IMBLOC, INB1, INBLOC, IOFFA, IOFFD, ITYPE,
     $                   IUPP, JJA, JJMAX, JOFFA, JOFFD, LCMT, LCMT00,
     $                   LDA, LMBLOC, LNBLOC, LOW, M1, MB, MBLKD, MBLKS,
     $                   MBLOC, MP, MRCOL, MRROW, MYCOL, MYROW, N1, NB,
     $                   NBLKD, NBLKS, NBLOC, NPCOL, NPROW, NQ, PMB,
     $                   QNB, TMP1, UPP
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA2( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PB_AINFOG2L, PB_BINFO,
     $                   PB_DESCTRANS, PB_INFOG2L, PB_ZLASCAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            PB_NUMROC
      EXTERNAL           LSAME, PB_NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Convert descriptor
*
      CALL PB_DESCTRANS( DESCA, DESCA2 )
*
*     Get grid parameters
*
      ICTXT = DESCA2( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
         UPLO  = TYPE
         UPPER = .FALSE.
         LOWER = .TRUE.
         IOFFD = 0
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
         UPLO  = TYPE
         UPPER = .TRUE.
         LOWER = .FALSE.
         IOFFD = 0
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
         UPLO  = 'U'
         UPPER = .TRUE.
         LOWER = .FALSE.
         IOFFD = 1
      ELSE
         ITYPE = 0
         UPLO  = 'A'
         UPPER = .TRUE.
         LOWER = .TRUE.
         IOFFD = 0
      END IF
*
*     Compute local indexes
*
      IF( ITYPE.EQ.0 ) THEN
*
*        Full matrix
*
         CALL PB_INFOG2L( IA, JA, DESCA2, NPROW, NPCOL, MYROW, MYCOL,
     $                    IIA, JJA, IAROW, IACOL )
         MP = PB_NUMROC( M, IA, DESCA2( IMB_ ), DESCA2( MB_ ), MYROW,
     $                   DESCA2( RSRC_ ), NPROW )
         NQ = PB_NUMROC( N, JA, DESCA2( INB_ ), DESCA2( NB_ ), MYCOL,
     $                   DESCA2( CSRC_ ), NPCOL )
*
         IF( MP.LE.0 .OR. NQ.LE.0 )
     $      RETURN
*
         LDA   = DESCA2( LLD_ )
         IOFFA = IIA + ( JJA - 1 ) * LDA
*
         CALL PB_ZLASCAL( 'All', MP, NQ, 0, ALPHA, A( IOFFA ), LDA )
*
      ELSE
*
*        Trapezoidal matrix
*
         CALL PB_AINFOG2L( M, N, IA, JA, DESCA2, NPROW, NPCOL, MYROW,
     $                     MYCOL, IMB1, INB1, MP, NQ, IIA, JJA, IAROW,
     $                     IACOL, MRROW, MRCOL )
*
         IF( MP.LE.0 .OR. NQ.LE.0 )
     $      RETURN
*
*        Initialize LCMT00, MBLKS, NBLKS, IMBLOC, INBLOC, LMBLOC,
*        LNBLOC, ILOW, LOW, IUPP, and UPP.
*
         MB  = DESCA2( MB_ )
         NB  = DESCA2( NB_ )
         LDA = DESCA2( LLD_ )
*
         CALL PB_BINFO( IOFFD, MP, NQ, IMB1, INB1, MB, NB, MRROW,
     $                  MRCOL, LCMT00, MBLKS, NBLKS, IMBLOC, INBLOC,
     $                  LMBLOC, LNBLOC, ILOW, LOW, IUPP, UPP )
*
         M1    = MP
         N1    = NQ
         IOFFA = IIA - 1
         JOFFA = JJA - 1
         IIMAX = IOFFA + MP
         JJMAX = JOFFA + NQ
*
         IF( DESCA2( RSRC_ ).LT.0 ) THEN
            PMB = MB
         ELSE
            PMB = NPROW * MB
         END IF
         IF( DESCA2( CSRC_ ).LT.0 ) THEN
            QNB = NB
         ELSE
            QNB = NPCOL * NB
         END IF
*
*        Handle the first block of rows or columns separately, and
*        update LCMT00, MBLKS and NBLKS.
*
         GODOWN = ( LCMT00.GT.IUPP )
         GOLEFT = ( LCMT00.LT.ILOW )
*
         IF( .NOT.GODOWN .AND. .NOT.GOLEFT ) THEN
*
*           LCMT00 >= ILOW && LCMT00 <= IUPP
*
            GOLEFT = ( ( LCMT00 - ( IUPP - UPP + PMB ) ).LT.ILOW )
            GODOWN = .NOT.GOLEFT
*
            CALL PB_ZLASCAL( UPLO, IMBLOC, INBLOC, LCMT00, ALPHA,
     $                       A( IIA+JOFFA*LDA ), LDA )
            IF( GODOWN ) THEN
               IF( UPPER .AND. NQ.GT.INBLOC )
     $            CALL PB_ZLASCAL( 'All', IMBLOC, NQ-INBLOC, 0, ALPHA,
     $                             A( IIA+(JOFFA+INBLOC)*LDA ), LDA )
               IIA = IIA + IMBLOC
               M1  = M1 - IMBLOC
            ELSE
               IF( LOWER .AND. MP.GT.IMBLOC )
     $            CALL PB_ZLASCAL( 'All', MP-IMBLOC, INBLOC, 0, ALPHA,
     $                             A( IIA+IMBLOC+JOFFA*LDA ), LDA )
               JJA = JJA + INBLOC
               N1  = N1 - INBLOC
            END IF
*
         END IF
*
         IF( GODOWN ) THEN
*
            LCMT00 = LCMT00 - ( IUPP - UPP + PMB )
            MBLKS  = MBLKS - 1
            IOFFA  = IOFFA + IMBLOC
*
   10       CONTINUE
            IF( MBLKS.GT.0 .AND. LCMT00.GT.UPP ) THEN
               LCMT00 = LCMT00 - PMB
               MBLKS  = MBLKS - 1
               IOFFA  = IOFFA + MB
               GO TO 10
            END IF
*
            TMP1 = MIN( IOFFA, IIMAX ) - IIA + 1
            IF( UPPER .AND. TMP1.GT.0 ) THEN
               CALL PB_ZLASCAL( 'All', TMP1, N1, 0, ALPHA,
     $                          A( IIA+JOFFA*LDA ), LDA )
               IIA = IIA + TMP1
               M1  = M1 - TMP1
            END IF
*
            IF( MBLKS.LE.0 )
     $         RETURN
*
            LCMT  = LCMT00
            MBLKD = MBLKS
            IOFFD = IOFFA
*
            MBLOC = MB
   20       CONTINUE
            IF( MBLKD.GT.0 .AND. LCMT.GE.ILOW ) THEN
               IF( MBLKD.EQ.1 )
     $            MBLOC = LMBLOC
               CALL PB_ZLASCAL( UPLO, MBLOC, INBLOC, LCMT, ALPHA,
     $                          A( IOFFD+1+JOFFA*LDA ), LDA )
               LCMT00 = LCMT
               LCMT   = LCMT - PMB
               MBLKS  = MBLKD
               MBLKD  = MBLKD - 1
               IOFFA  = IOFFD
               IOFFD  = IOFFD + MBLOC
               GO TO 20
            END IF
*
            TMP1 = M1 - IOFFD + IIA - 1
            IF( LOWER .AND. TMP1.GT.0 )
     $         CALL PB_ZLASCAL( 'All', TMP1, INBLOC, 0, ALPHA,
     $                          A( IOFFD+1+JOFFA*LDA ), LDA )
*
            TMP1   = IOFFA - IIA + 1
            M1     = M1 - TMP1
            N1     = N1 - INBLOC
            LCMT00 = LCMT00 + LOW - ILOW + QNB
            NBLKS  = NBLKS - 1
            JOFFA  = JOFFA + INBLOC
*
            IF( UPPER .AND. TMP1.GT.0 .AND. N1.GT.0 )
     $         CALL PB_ZLASCAL( 'All', TMP1, N1, 0, ALPHA,
     $                          A( IIA+JOFFA*LDA ), LDA )
*
            IIA = IOFFA + 1
            JJA = JOFFA + 1
*
         ELSE IF( GOLEFT ) THEN
*
            LCMT00 = LCMT00 + LOW - ILOW + QNB
            NBLKS  = NBLKS - 1
            JOFFA  = JOFFA + INBLOC
*
   30       CONTINUE
            IF( NBLKS.GT.0 .AND. LCMT00.LT.LOW ) THEN
               LCMT00 = LCMT00 + QNB
               NBLKS  = NBLKS - 1
               JOFFA  = JOFFA + NB
               GO TO 30
            END IF
*
            TMP1 = MIN( JOFFA, JJMAX ) - JJA + 1
            IF( LOWER .AND. TMP1.GT.0 ) THEN
               CALL PB_ZLASCAL( 'All', M1, TMP1, 0, ALPHA,
     $                          A( IIA+(JJA-1)*LDA ), LDA )
               JJA = JJA + TMP1
               N1  = N1 - TMP1
            END IF
*
            IF( NBLKS.LE.0 )
     $         RETURN
*
            LCMT  = LCMT00
            NBLKD = NBLKS
            JOFFD = JOFFA
*
            NBLOC = NB
   40       CONTINUE
            IF( NBLKD.GT.0 .AND. LCMT.LE.IUPP ) THEN
               IF( NBLKD.EQ.1 )
     $            NBLOC = LNBLOC
               CALL PB_ZLASCAL( UPLO, IMBLOC, NBLOC, LCMT, ALPHA,
     $                          A( IIA+JOFFD*LDA ), LDA )
               LCMT00 = LCMT
               LCMT   = LCMT + QNB
               NBLKS  = NBLKD
               NBLKD  = NBLKD - 1
               JOFFA  = JOFFD
               JOFFD  = JOFFD + NBLOC
               GO TO 40
            END IF
*
            TMP1 = N1 - JOFFD + JJA - 1
            IF( UPPER .AND. TMP1.GT.0 )
     $         CALL PB_ZLASCAL( 'All', IMBLOC, TMP1, 0, ALPHA,
     $                          A( IIA+JOFFD*LDA ), LDA )
*
            TMP1   = JOFFA - JJA + 1
            M1     = M1 - IMBLOC
            N1     = N1 - TMP1
            LCMT00 = LCMT00 - ( IUPP - UPP + PMB )
            MBLKS  = MBLKS - 1
            IOFFA  = IOFFA + IMBLOC
*
            IF( LOWER .AND. M1.GT.0 .AND. TMP1.GT.0 )
     $         CALL PB_ZLASCAL( 'All', M1, TMP1, 0, ALPHA,
     $                          A( IOFFA+1+(JJA-1)*LDA ), LDA )
*
            IIA = IOFFA + 1
            JJA = JOFFA + 1
*
         END IF
*
         NBLOC = NB
   50    CONTINUE
         IF( NBLKS.GT.0 ) THEN
            IF( NBLKS.EQ.1 )
     $         NBLOC = LNBLOC
   60       CONTINUE
            IF( MBLKS.GT.0 .AND. LCMT00.GT.UPP ) THEN
               LCMT00 = LCMT00 - PMB
               MBLKS  = MBLKS - 1
               IOFFA  = IOFFA + MB
               GO TO 60
            END IF
*
            TMP1 = MIN( IOFFA, IIMAX ) - IIA + 1
            IF( UPPER .AND. TMP1.GT.0 ) THEN
               CALL PB_ZLASCAL( 'All', TMP1, N1, 0, ALPHA,
     $                          A( IIA+JOFFA*LDA ), LDA )
               IIA = IIA + TMP1
               M1  = M1 - TMP1
            END IF
*
            IF( MBLKS.LE.0 )
     $         RETURN
*
            LCMT  = LCMT00
            MBLKD = MBLKS
            IOFFD = IOFFA
*
            MBLOC = MB
   70       CONTINUE
            IF( MBLKD.GT.0 .AND. LCMT.GE.LOW ) THEN
               IF( MBLKD.EQ.1 )
     $            MBLOC = LMBLOC
               CALL PB_ZLASCAL( UPLO, MBLOC, NBLOC, LCMT, ALPHA,
     $                          A( IOFFD+1+JOFFA*LDA ), LDA )
               LCMT00 = LCMT
               LCMT   = LCMT - PMB
               MBLKS  = MBLKD
               MBLKD  = MBLKD - 1
               IOFFA  = IOFFD
               IOFFD  = IOFFD + MBLOC
               GO TO 70
            END IF
*
            TMP1 = M1 - IOFFD + IIA - 1
            IF( LOWER .AND. TMP1.GT.0 )
     $         CALL PB_ZLASCAL( 'All', TMP1, NBLOC, 0, ALPHA,
     $                          A( IOFFD+1+JOFFA*LDA ), LDA )
*
            TMP1   = MIN( IOFFA, IIMAX )  - IIA + 1
            M1     = M1 - TMP1
            N1     = N1 - NBLOC
            LCMT00 = LCMT00 + QNB
            NBLKS  = NBLKS - 1
            JOFFA  = JOFFA + NBLOC
*
            IF( UPPER .AND. TMP1.GT.0 .AND. N1.GT.0 )
     $         CALL PB_ZLASCAL( 'All', TMP1, N1, 0, ALPHA,
     $                          A( IIA+JOFFA*LDA ), LDA )
*
            IIA = IOFFA + 1
            JJA = JOFFA + 1
*
            GO TO 50
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PZLASCAL
*
      END
      SUBROUTINE PZLAGEN( INPLACE, AFORM, DIAG, OFFA, M, N, IA, JA,
     $                    DESCA, IASEED, A, LDA )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      LOGICAL            INPLACE
      CHARACTER*1        AFORM, DIAG
      INTEGER            IA, IASEED, JA, LDA, M, N, OFFA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  PZLAGEN  generates  (or regenerates)  a  submatrix  sub( A ) denoting
*  A(IA:IA+M-1,JA:JA+N-1).
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
*  INPLACE (global input) LOGICAL
*          On entry, INPLACE specifies if the matrix should be generated
*          in place or not. If INPLACE is .TRUE., the local random array
*          to be generated  will start in memory at the local memory lo-
*          cation A( 1, 1 ),  otherwise it will start at the local posi-
*          tion induced by IA and JA.
*
*  AFORM   (global input) CHARACTER*1
*          On entry, AFORM specifies the type of submatrix to be genera-
*          ted as follows:
*             AFORM = 'S', sub( A ) is a symmetric matrix,
*             AFORM = 'H', sub( A ) is a Hermitian matrix,
*             AFORM = 'T', sub( A ) is overrwritten  with  the transpose
*                          of what would normally be generated,
*             AFORM = 'C', sub( A ) is overwritten  with  the  conjugate
*                          transpose  of  what would normally be genera-
*                          ted.
*             AFORM = 'N', a random submatrix is generated.
*
*  DIAG    (global input) CHARACTER*1
*          On entry, DIAG specifies if the generated submatrix is diago-
*          nally dominant or not as follows:
*             DIAG = 'D' : sub( A ) is diagonally dominant,
*             DIAG = 'N' : sub( A ) is not diagonally dominant.
*
*  OFFA    (global input) INTEGER
*          On entry, OFFA  specifies  the  offdiagonal of the underlying
*          matrix A(1:DESCA(M_),1:DESCA(N_)) of interest when the subma-
*          trix is symmetric, Hermitian or diagonally dominant. OFFA = 0
*          specifies the main diagonal,  OFFA > 0  specifies a subdiago-
*          nal,  and OFFA < 0 specifies a superdiagonal (see further de-
*          tails).
*
*  M       (global input) INTEGER
*          On entry, M specifies the global number of matrix rows of the
*          submatrix sub( A ) to be generated. M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry,  N specifies the global number of matrix columns of
*          the  submatrix  sub( A )  to be generated. N must be at least
*          zero.
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
*  IASEED  (global input) INTEGER
*          On entry, IASEED  specifies  the  seed number to generate the
*          matrix A. IASEED must be at least zero.
*
*  A       (local output) COMPLEX*16 array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  On  exit, this array  contains the
*          local entries of the randomly generated submatrix sub( A ).
*
*  LDA     (local input) INTEGER
*          On entry,  LDA  specifies  the local leading dimension of the
*          array A. When INPLACE is .FALSE., LDA is usually DESCA(LLD_).
*          This restriction is however not enforced, and this subroutine
*          requires only that LDA >= MAX( 1, Mp ) where
*
*          Mp = PB_NUMROC( M, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW ).
*
*          PB_NUMROC  is  a ScaLAPACK tool function; MYROW, MYCOL, NPROW
*          and NPCOL  can  be determined by calling the BLACS subroutine
*          BLACS_GRIDINFO.
*
*  Further Details
*  ===============
*
*  OFFD  is  tied  to  the matrix described by  DESCA, as opposed to the
*  piece that is currently  (re)generated.  This is a global information
*  independent from the distribution  parameters.  Below are examples of
*  the meaning of OFFD for a global 7 by 5 matrix:
*
*  ---------------------------------------------------------------------
*  OFFD   |  0 -1 -2 -3 -4         0 -1 -2 -3 -4          0 -1 -2 -3 -4
*  -------|-------------------------------------------------------------
*         |     | OFFD=-1          |   OFFD=0                 OFFD=2
*         |     V                  V
*  0      |  .  d  .  .  .      -> d  .  .  .  .          .  .  .  .  .
*  1      |  .  .  d  .  .         .  d  .  .  .          .  .  .  .  .
*  2      |  .  .  .  d  .         .  .  d  .  .       -> d  .  .  .  .
*  3      |  .  .  .  .  d         .  .  .  d  .          .  d  .  .  .
*  4      |  .  .  .  .  .         .  .  .  .  d          .  .  d  .  .
*  5      |  .  .  .  .  .         .  .  .  .  .          .  .  .  d  .
*  6      |  .  .  .  .  .         .  .  .  .  .          .  .  .  .  d
*  ---------------------------------------------------------------------
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
      INTEGER            JMP_1, JMP_COL, JMP_IMBV, JMP_INBV, JMP_LEN,
     $                   JMP_MB, JMP_NB, JMP_NPIMBLOC, JMP_NPMB,
     $                   JMP_NQINBLOC, JMP_NQNB, JMP_ROW
      PARAMETER          ( JMP_1 = 1, JMP_ROW = 2, JMP_COL = 3,
     $                   JMP_MB = 4, JMP_IMBV = 5, JMP_NPMB = 6,
     $                   JMP_NPIMBLOC = 7, JMP_NB = 8, JMP_INBV = 9,
     $                   JMP_NQNB = 10, JMP_NQINBLOC = 11,
     $                   JMP_LEN = 11 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DIAGDO, SYMM, HERM, NOTRAN
      INTEGER            CSRC, I, IACOL, IAROW, ICTXT, IIA, ILOCBLK,
     $                   ILOCOFF, ILOW, IMB, IMB1, IMBLOC, IMBVIR, INB,
     $                   INB1, INBLOC, INBVIR, INFO, IOFFDA, ITMP, IUPP,
     $                   IVIR, JJA, JLOCBLK, JLOCOFF, JVIR, LCMT00,
     $                   LMBLOC, LNBLOC, LOW, MAXMN, MB, MBLKS, MP,
     $                   MRCOL, MRROW, MYCDIST, MYCOL, MYRDIST, MYROW,
     $                   NB, NBLKS, NPCOL, NPROW, NQ, NVIR, RSRC, UPP
      COMPLEX*16         ALPHA
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA2( DLEN_ ), IMULADD( 4, JMP_LEN ),
     $                   IRAN( 2 ), JMP( JMP_LEN ), MULADD0( 4 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PB_AINFOG2L, PB_BINFO,
     $                   PB_CHKMAT, PB_DESCTRANS, PB_INITJMP,
     $                   PB_INITMULADD, PB_JUMP, PB_JUMPIT, PB_LOCINFO,
     $                   PB_SETLOCRAN, PB_SETRAN, PB_ZLAGEN, PXERBLA,
     $                   PZLADOM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Data Statements ..
      DATA               ( MULADD0( I ), I = 1, 4 ) / 20077, 16838,
     $                   12345, 0 /
*     ..
*     .. Executable Statements ..
*
*     Convert descriptor
*
      CALL PB_DESCTRANS( DESCA, DESCA2 )
*
*     Test the input arguments
*
      ICTXT = DESCA2( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test the input parameters
*
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -( 1000 + CTXT_ )
      ELSE
         SYMM   = LSAME( AFORM, 'S' )
         HERM   = LSAME( AFORM, 'H' )
         NOTRAN = LSAME( AFORM, 'N' )
         DIAGDO = LSAME( DIAG, 'D' )
         IF( .NOT.( SYMM.OR.HERM.OR.NOTRAN ) .AND.
     $       .NOT.( LSAME( AFORM, 'T' )    ) .AND.
     $       .NOT.( LSAME( AFORM, 'C' )    ) ) THEN
            INFO = -2
         ELSE IF( ( .NOT.DIAGDO ) .AND.
     $            ( .NOT.LSAME( DIAG, 'N' ) ) ) THEN
            INFO = -3
         END IF
         CALL PB_CHKMAT( ICTXT, M, 5, N, 6, IA, JA, DESCA2, 10, INFO )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PZLAGEN', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( ( M.LE.0 ).OR.( N.LE.0 ) )
     $   RETURN
*
*     Start the operations
*
      MB   = DESCA2( MB_   )
      NB   = DESCA2( NB_   )
      IMB  = DESCA2( IMB_  )
      INB  = DESCA2( INB_  )
      RSRC = DESCA2( RSRC_ )
      CSRC = DESCA2( CSRC_ )
*
*     Figure out local information about the distributed matrix operand
*
      CALL PB_AINFOG2L( M, N, IA, JA, DESCA2, NPROW, NPCOL, MYROW,
     $                  MYCOL, IMB1, INB1, MP, NQ, IIA, JJA, IAROW,
     $                  IACOL, MRROW, MRCOL )
*
*     Decide where the entries shall be stored in memory
*
      IF( INPLACE ) THEN
         IIA = 1
         JJA = 1
      END IF
*
*     Initialize LCMT00, MBLKS, NBLKS, IMBLOC, INBLOC, LMBLOC, LNBLOC,
*     ILOW, LOW, IUPP, and UPP.
*
      IOFFDA = JA + OFFA - IA
      CALL PB_BINFO( IOFFDA, MP, NQ, IMB1, INB1, MB, NB, MRROW,
     $               MRCOL, LCMT00, MBLKS, NBLKS, IMBLOC, INBLOC,
     $               LMBLOC, LNBLOC, ILOW, LOW, IUPP, UPP )
*
*     Initialize ILOCBLK, ILOCOFF, MYRDIST, JLOCBLK, JLOCOFF, MYCDIST
*     This values correspond to the square virtual underlying matrix
*     of size MAX( M_ + MAX( 0, -OFFA ), N_ + MAX( 0, OFFA ) ) used
*     to set up the random sequence. For practical purposes, the size
*     of this virtual matrix is upper bounded by M_ + N_ - 1.
*
      ITMP   = MAX( 0, -OFFA )
      IVIR   = IA  + ITMP
      IMBVIR = IMB + ITMP
      NVIR   = DESCA2( M_ ) + ITMP
*
      CALL PB_LOCINFO( IVIR, IMBVIR, MB, MYROW, RSRC, NPROW, ILOCBLK,
     $                 ILOCOFF, MYRDIST )
*
      ITMP   = MAX( 0, OFFA )
      JVIR   = JA  + ITMP
      INBVIR = INB + ITMP
      NVIR   = MAX( MAX( NVIR, DESCA2( N_ ) + ITMP ),
     $              DESCA2( M_ ) + DESCA2( N_ ) - 1 )
*
      CALL PB_LOCINFO( JVIR, INBVIR, NB, MYCOL, CSRC, NPCOL, JLOCBLK,
     $                 JLOCOFF, MYCDIST )
*
      IF( SYMM .OR. HERM .OR. NOTRAN ) THEN
*
         CALL PB_INITJMP( .TRUE., NVIR, IMBVIR, INBVIR, IMBLOC, INBLOC,
     $                    MB, NB, RSRC, CSRC, NPROW, NPCOL, 2, JMP )
*
*        Compute constants to jump JMP( * ) numbers in the sequence
*
         CALL PB_INITMULADD( MULADD0, JMP, IMULADD )
*
*        Compute and set the random value corresponding to A( IA, JA )
*
         CALL PB_SETLOCRAN( IASEED, ILOCBLK, JLOCBLK, ILOCOFF, JLOCOFF,
     $                      MYRDIST, MYCDIST, NPROW, NPCOL, JMP,
     $                      IMULADD, IRAN )
*
         CALL PB_ZLAGEN( 'Lower', AFORM, A( IIA, JJA ), LDA, LCMT00,
     $                   IRAN, MBLKS, IMBLOC, MB, LMBLOC, NBLKS, INBLOC,
     $                   NB, LNBLOC, JMP, IMULADD )
*
      END IF
*
      IF( SYMM .OR. HERM .OR. ( .NOT. NOTRAN ) ) THEN
*
         CALL PB_INITJMP( .FALSE., NVIR, IMBVIR, INBVIR, IMBLOC, INBLOC,
     $                    MB, NB, RSRC, CSRC, NPROW, NPCOL, 2, JMP )
*
*        Compute constants to jump JMP( * ) numbers in the sequence
*
         CALL PB_INITMULADD( MULADD0, JMP, IMULADD )
*
*        Compute and set the random value corresponding to A( IA, JA )
*
         CALL PB_SETLOCRAN( IASEED, ILOCBLK, JLOCBLK, ILOCOFF, JLOCOFF,
     $                      MYRDIST, MYCDIST, NPROW, NPCOL, JMP,
     $                      IMULADD, IRAN )
*
         CALL PB_ZLAGEN( 'Upper', AFORM, A( IIA, JJA ), LDA, LCMT00,
     $                   IRAN, MBLKS, IMBLOC, MB, LMBLOC, NBLKS, INBLOC,
     $                   NB, LNBLOC, JMP, IMULADD )
*
      END IF
*
      IF( DIAGDO ) THEN
*
         MAXMN = MAX( DESCA2( M_ ), DESCA2( N_ ) )
         IF( HERM ) THEN
            ALPHA = DCMPLX( DBLE( 2 * MAXMN ), ZERO )
         ELSE
            ALPHA = DCMPLX( DBLE( NVIR ), DBLE( MAXMN ) )
         END IF
*
         IF( IOFFDA.GE.0 ) THEN
            CALL PZLADOM( INPLACE, MIN( MAX( 0, M-IOFFDA ), N ), ALPHA,
     $                    A, MIN( IA+IOFFDA, IA+M-1 ), JA, DESCA )
         ELSE
            CALL PZLADOM( INPLACE, MIN( M, MAX( 0, N+IOFFDA ) ), ALPHA,
     $                    A, IA, MIN( JA-IOFFDA, JA+N-1 ), DESCA )
         END IF
*
      END IF
*
      RETURN
*
*     End of PZLAGEN
*
      END
      SUBROUTINE PZLADOM( INPLACE, N, ALPHA, A, IA, JA, DESCA )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      LOGICAL            INPLACE
      INTEGER            IA, JA, N
      COMPLEX*16         ALPHA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         A( * )
*     ..
*
*  Purpose
*  =======
*
*  PZLADOM  adds alpha to the diagonal entries  of  an  n by n submatrix
*  sub( A ) denoting A( IA:IA+N-1, JA:JA+N-1 ).
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
*  INPLACE (global input) LOGICAL
*          On entry, INPLACE specifies if the matrix should be generated
*          in place or not. If INPLACE is .TRUE., the local random array
*          to be generated  will start in memory at the local memory lo-
*          cation A( 1, 1 ),  otherwise it will start at the local posi-
*          tion induced by IA and JA.
*
*  N       (global input) INTEGER
*          On entry,  N  specifies  the  global  order  of the submatrix
*          sub( A ) to be modified. N must be at least zero.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (local input/local output) COMPLEX*16 array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix A. On exit, the local entries
*          of this array corresponding to the main diagonal of  sub( A )
*          have been updated.
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
      LOGICAL            GODOWN, GOLEFT
      INTEGER            I, IACOL, IAROW, ICTXT, IIA, IJOFFA, ILOW,
     $                   IMB1, IMBLOC, INB1, INBLOC, IOFFA, IOFFD, IUPP,
     $                   JJA, JOFFA, JOFFD, LCMT, LCMT00, LDA, LDAP1,
     $                   LMBLOC, LNBLOC, LOW, MB, MBLKD, MBLKS, MBLOC,
     $                   MRCOL, MRROW, MYCOL, MYROW, NB, NBLKD, NBLKS,
     $                   NBLOC, NP, NPCOL, NPROW, NQ, PMB, QNB, UPP
      COMPLEX*16         ATMP
*     ..
*     .. Local Scalars ..
      INTEGER            DESCA2( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PB_AINFOG2L, PB_BINFO,
     $                   PB_DESCTRANS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Convert descriptor
*
      CALL PB_DESCTRANS( DESCA, DESCA2 )
*
*     Get grid parameters
*
      ICTXT = DESCA2( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( N.EQ.0 )
     $   RETURN
*
      CALL PB_AINFOG2L( N, N, IA, JA, DESCA2, NPROW, NPCOL, MYROW,
     $                  MYCOL, IMB1, INB1, NP, NQ, IIA, JJA, IAROW,
     $                  IACOL, MRROW, MRCOL )
*
*     Decide where the entries shall be stored in memory
*
      IF( INPLACE ) THEN
         IIA = 1
         JJA = 1
      END IF
*
*     Initialize LCMT00, MBLKS, NBLKS, IMBLOC, INBLOC, LMBLOC, LNBLOC,
*     ILOW, LOW, IUPP, and UPP.
*
      MB = DESCA2( MB_ )
      NB = DESCA2( NB_ )
*
      CALL PB_BINFO( 0, NP, NQ, IMB1, INB1, MB, NB, MRROW, MRCOL,
     $               LCMT00, MBLKS, NBLKS, IMBLOC, INBLOC, LMBLOC,
     $               LNBLOC, ILOW, LOW, IUPP, UPP )
*
      IOFFA  = IIA - 1
      JOFFA  = JJA - 1
      LDA    = DESCA2( LLD_ )
      LDAP1  = LDA + 1
*
      IF( DESCA2( RSRC_ ).LT.0 ) THEN
         PMB = MB
      ELSE
         PMB = NPROW * MB
      END IF
      IF( DESCA2( CSRC_ ).LT.0 ) THEN
         QNB = NB
      ELSE
         QNB = NPCOL * NB
      END IF
*
*     Handle the first block of rows or columns separately, and update
*     LCMT00, MBLKS and NBLKS.
*
      GODOWN = ( LCMT00.GT.IUPP )
      GOLEFT = ( LCMT00.LT.ILOW )
*
      IF( .NOT.GODOWN .AND. .NOT.GOLEFT ) THEN
*
*        LCMT00 >= ILOW && LCMT00 <= IUPP
*
         IF( LCMT00.GE.0 ) THEN
            IJOFFA = IOFFA+LCMT00 + ( JOFFA - 1 ) * LDA
            DO 10 I = 1, MIN( INBLOC, MAX( 0, IMBLOC - LCMT00 ) )
               ATMP = A( IJOFFA + I*LDAP1 )
               A( IJOFFA + I*LDAP1 ) = ALPHA +
     $                                 DCMPLX( ABS( DBLE(  ATMP ) ),
     $                                         ABS( DIMAG( ATMP ) ) )
   10       CONTINUE
         ELSE
            IJOFFA = IOFFA + ( JOFFA - LCMT00 - 1 ) * LDA
            DO 20 I = 1, MIN( IMBLOC, MAX( 0, INBLOC + LCMT00 ) )
               ATMP = A( IJOFFA + I*LDAP1 )
               A( IJOFFA + I*LDAP1 ) = ALPHA +
     $                                 DCMPLX( ABS( DBLE(  ATMP ) ),
     $                                         ABS( DIMAG( ATMP ) ) )
   20       CONTINUE
         END IF
         GOLEFT = ( ( LCMT00 - ( IUPP - UPP + PMB ) ).LT.ILOW )
         GODOWN = .NOT.GOLEFT
*
      END IF
*
      IF( GODOWN ) THEN
*
         LCMT00 = LCMT00 - ( IUPP - UPP + PMB )
         MBLKS  = MBLKS - 1
         IOFFA  = IOFFA + IMBLOC
*
   30    CONTINUE
         IF( MBLKS.GT.0 .AND. LCMT00.GT.UPP ) THEN
            LCMT00 = LCMT00 - PMB
            MBLKS  = MBLKS - 1
            IOFFA  = IOFFA + MB
            GO TO 30
         END IF
*
         LCMT  = LCMT00
         MBLKD = MBLKS
         IOFFD = IOFFA
*
         MBLOC = MB
   40    CONTINUE
         IF( MBLKD.GT.0 .AND. LCMT.GE.ILOW ) THEN
            IF( MBLKD.EQ.1 )
     $         MBLOC = LMBLOC
            IF( LCMT.GE.0 ) THEN
               IJOFFA = IOFFD + LCMT + ( JOFFA - 1 ) * LDA
               DO 50 I = 1, MIN( INBLOC, MAX( 0, MBLOC - LCMT ) )
                  ATMP = A( IJOFFA + I*LDAP1 )
                  A( IJOFFA + I*LDAP1 ) = ALPHA +
     $                                    DCMPLX( ABS( DBLE(  ATMP ) ),
     $                                            ABS( DIMAG( ATMP ) ) )
   50          CONTINUE
            ELSE
               IJOFFA = IOFFD + ( JOFFA - LCMT - 1 ) * LDA
               DO 60 I = 1, MIN( MBLOC, MAX( 0, INBLOC + LCMT ) )
                  ATMP = A( IJOFFA + I*LDAP1 )
                  A( IJOFFA + I*LDAP1 ) = ALPHA +
     $                                    DCMPLX( ABS( DBLE(  ATMP ) ),
     $                                            ABS( DIMAG( ATMP ) ) )
   60          CONTINUE
            END IF
            LCMT00 = LCMT
            LCMT   = LCMT - PMB
            MBLKS  = MBLKD
            MBLKD  = MBLKD - 1
            IOFFA  = IOFFD
            IOFFD  = IOFFD + MBLOC
            GO TO 40
         END IF
*
         LCMT00 = LCMT00 + LOW - ILOW + QNB
         NBLKS  = NBLKS - 1
         JOFFA  = JOFFA + INBLOC
*
      ELSE IF( GOLEFT ) THEN
*
         LCMT00 = LCMT00 + LOW - ILOW + QNB
         NBLKS  = NBLKS - 1
         JOFFA  = JOFFA + INBLOC
*
   70    CONTINUE
         IF( NBLKS.GT.0 .AND. LCMT00.LT.LOW ) THEN
            LCMT00 = LCMT00 + QNB
            NBLKS  = NBLKS - 1
            JOFFA  = JOFFA + NB
            GO TO 70
         END IF
*
         LCMT  = LCMT00
         NBLKD = NBLKS
         JOFFD = JOFFA
*
         NBLOC = NB
   80    CONTINUE
         IF( NBLKD.GT.0 .AND. LCMT.LE.IUPP ) THEN
            IF( NBLKD.EQ.1 )
     $         NBLOC = LNBLOC
            IF( LCMT.GE.0 ) THEN
               IJOFFA = IOFFA + LCMT + ( JOFFD - 1 ) * LDA
               DO 90 I = 1, MIN( NBLOC, MAX( 0, IMBLOC - LCMT ) )
                  ATMP = A( IJOFFA + I*LDAP1 )
                  A( IJOFFA + I*LDAP1 ) = ALPHA +
     $                                    DCMPLX( ABS( DBLE(  ATMP ) ),
     $                                            ABS( DIMAG( ATMP ) ) )
   90          CONTINUE
            ELSE
               IJOFFA = IOFFA + ( JOFFD - LCMT - 1 ) * LDA
               DO 100 I = 1, MIN( IMBLOC, MAX( 0, NBLOC + LCMT ) )
                  ATMP = A( IJOFFA + I*LDAP1 )
                  A( IJOFFA + I*LDAP1 ) = ALPHA +
     $                                    DCMPLX( ABS( DBLE(  ATMP ) ),
     $                                            ABS( DIMAG( ATMP ) ) )
  100          CONTINUE
            END IF
            LCMT00 = LCMT
            LCMT   = LCMT + QNB
            NBLKS  = NBLKD
            NBLKD  = NBLKD - 1
            JOFFA  = JOFFD
            JOFFD  = JOFFD + NBLOC
            GO TO 80
         END IF
*
         LCMT00 = LCMT00 - ( IUPP - UPP + PMB )
         MBLKS  = MBLKS - 1
         IOFFA  = IOFFA + IMBLOC
*
      END IF
*
      NBLOC = NB
  110 CONTINUE
      IF( NBLKS.GT.0 ) THEN
         IF( NBLKS.EQ.1 )
     $      NBLOC = LNBLOC
  120    CONTINUE
         IF( MBLKS.GT.0 .AND. LCMT00.GT.UPP ) THEN
            LCMT00 = LCMT00 - PMB
            MBLKS  = MBLKS - 1
            IOFFA  = IOFFA + MB
            GO TO 120
         END IF
*
         LCMT  = LCMT00
         MBLKD = MBLKS
         IOFFD = IOFFA
*
         MBLOC = MB
  130    CONTINUE
         IF( MBLKD.GT.0 .AND. LCMT.GE.LOW ) THEN
            IF( MBLKD.EQ.1 )
     $         MBLOC = LMBLOC
            IF( LCMT.GE.0 ) THEN
               IJOFFA = IOFFD + LCMT + ( JOFFA - 1 ) * LDA
               DO 140 I = 1, MIN( NBLOC, MAX( 0, MBLOC - LCMT ) )
                  ATMP = A( IJOFFA + I*LDAP1 )
                  A( IJOFFA + I*LDAP1 ) = ALPHA +
     $                                    DCMPLX( ABS( DBLE(  ATMP ) ),
     $                                            ABS( DIMAG( ATMP ) ) )
  140          CONTINUE
            ELSE
               IJOFFA = IOFFD + ( JOFFA - LCMT - 1 ) * LDA
               DO 150 I = 1, MIN( MBLOC, MAX( 0, NBLOC + LCMT ) )
                  ATMP = A( IJOFFA + I*LDAP1 )
                  A( IJOFFA + I*LDAP1 ) = ALPHA +
     $                                    DCMPLX( ABS( DBLE(  ATMP ) ),
     $                                            ABS( DIMAG( ATMP ) ) )
  150          CONTINUE
            END IF
            LCMT00 = LCMT
            LCMT   = LCMT - PMB
            MBLKS  = MBLKD
            MBLKD  = MBLKD - 1
            IOFFA  = IOFFD
            IOFFD  = IOFFD + MBLOC
            GO TO 130
         END IF
*
         LCMT00 = LCMT00 + QNB
         NBLKS  = NBLKS - 1
         JOFFA  = JOFFA + NBLOC
         GO TO 110
*
      END IF
*
      RETURN
*
*     End of PZLADOM
*
      END
      SUBROUTINE PB_PZLAPRNT( M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT,
     $                        CMATNM, NOUT, WORK )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            IA, ICPRNT, IRPRNT, JA, M, N, NOUT
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      CMATNM
      INTEGER            DESCA( * )
      COMPLEX*16         A( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PB_PZLAPRNT  prints to the standard output a submatrix sub( A ) deno-
*  ting A(IA:IA+M-1,JA:JA+N-1). The local pieces are sent and printed by
*  the process of coordinates (IRPRNT, ICPRNT).
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
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( A ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( A ). N must be at least zero.
*
*  A       (local input) COMPLEX*16 array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix A.
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
*  IRPRNT  (global input) INTEGER
*          On entry, IRPRNT specifies the row index of the printing pro-
*          cess.
*
*  ICPRNT  (global input) INTEGER
*          On entry, ICPRNT specifies the  column  index of the printing
*          process.
*
*  CMATNM  (global input) CHARACTER*(*)
*          On entry, CMATNM is the name of the matrix to be printed.
*
*  NOUT    (global input) INTEGER
*          On entry, NOUT specifies the output unit number. When NOUT is
*          equal to 6, the submatrix is printed on the screen.
*
*  WORK    (local workspace) COMPLEX*16 array
*          On entry, WORK is a work array of dimension at least equal to
*          MAX( IMB_A, MB_A ).
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
      INTEGER            MYCOL, MYROW, NPCOL, NPROW, PCOL, PROW
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA2( DLEN_ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PB_DESCTRANS, PB_PZLAPRN2
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( ( M.LE.0 ).OR.( N.LE.0 ) )
     $   RETURN
*
*     Convert descriptor
*
      CALL PB_DESCTRANS( DESCA, DESCA2 )
*
      CALL BLACS_GRIDINFO( DESCA2( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
      IF( DESCA2( RSRC_ ).GE.0 ) THEN
         IF( DESCA2( CSRC_ ).GE.0 ) THEN
            CALL PB_PZLAPRN2( M, N, A, IA, JA, DESCA2, IRPRNT, ICPRNT,
     $                        CMATNM, NOUT, DESCA2( RSRC_ ),
     $                        DESCA2( CSRC_ ), WORK )
         ELSE
            DO 10 PCOL = 0, NPCOL - 1
               IF( ( MYROW.EQ.IRPRNT ).AND.( MYCOL.EQ.ICPRNT ) )
     $            WRITE( NOUT, * ) 'Colum-replicated array -- ' ,
     $                             'copy in process column: ', PCOL
               CALL PB_PZLAPRN2( M, N, A, IA, JA, DESCA2, IRPRNT,
     $                           ICPRNT, CMATNM, NOUT, DESCA2( RSRC_ ),
     $                           PCOL, WORK )
   10       CONTINUE
         END IF
      ELSE
         IF( DESCA2( CSRC_ ).GE.0 ) THEN
            DO 20 PROW = 0, NPROW - 1
               IF( ( MYROW.EQ.IRPRNT ).AND.( MYCOL.EQ.ICPRNT ) )
     $            WRITE( NOUT, * ) 'Row-replicated array -- ' ,
     $                             'copy in process row: ', PROW
               CALL PB_PZLAPRN2( M, N, A, IA, JA, DESCA2, IRPRNT,
     $                           ICPRNT, CMATNM, NOUT, PROW,
     $                           DESCA2( CSRC_ ), WORK )
   20       CONTINUE
         ELSE
            DO 40 PROW = 0, NPROW - 1
               DO 30 PCOL = 0, NPCOL - 1
                  IF( ( MYROW.EQ.IRPRNT ).AND.( MYCOL.EQ.ICPRNT ) )
     $               WRITE( NOUT, * ) 'Replicated array -- ' ,
     $                      'copy in process (', PROW, ',', PCOL, ')'
                  CALL PB_PZLAPRN2( M, N, A, IA, JA, DESCA2, IRPRNT,
     $                              ICPRNT, CMATNM, NOUT, PROW, PCOL,
     $                              WORK )
   30          CONTINUE
   40       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of PB_PZLAPRNT
*
      END
      SUBROUTINE PB_PZLAPRN2( M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT,
     $                        CMATNM, NOUT, PROW, PCOL, WORK )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            IA, ICPRNT, IRPRNT, JA, M, N, NOUT, PCOL, PROW
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      CMATNM
      INTEGER            DESCA( * )
      COMPLEX*16         A( * ), WORK( * )
*     ..
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
      LOGICAL            AISCOLREP, AISROWREP
      INTEGER            H, I, IACOL, IAROW, IB, ICTXT, ICURCOL,
     $                   ICURROW, II, IIA, IN, J, JB, JJ, JJA, JN, K,
     $                   LDA, LDW, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_GRIDINFO, PB_INFOG2L,
     $                   ZGERV2D, ZGESD2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DIMAG, MIN
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      CALL PB_INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIA, JJA, IAROW, IACOL )
      II = IIA
      JJ = JJA
      IF( DESCA( RSRC_ ).LT.0 ) THEN
         AISROWREP = .TRUE.
         IAROW     = PROW
         ICURROW   = PROW
      ELSE
         AISROWREP = .FALSE.
         ICURROW   = IAROW
      END IF
      IF( DESCA( CSRC_ ).LT.0 ) THEN
         AISCOLREP = .TRUE.
         IACOL     = PCOL
         ICURCOL   = PCOL
      ELSE
         AISCOLREP = .FALSE.
         ICURCOL   = IACOL
      END IF
      LDA = DESCA( LLD_ )
      LDW = MAX( DESCA( IMB_ ), DESCA( MB_ ) )
*
*     Handle the first block of column separately
*
      JB = DESCA( INB_ ) - JA + 1
      IF( JB.LE.0 )
     $   JB = ( (-JB) / DESCA( NB_ ) + 1 ) * DESCA( NB_ ) + JB
      JB = MIN( JB, N )
      JN = JA+JB-1
      DO 60 H = 0, JB-1
         IB = DESCA( IMB_ ) - IA + 1
         IF( IB.LE.0 )
     $      IB = ( (-IB) / DESCA( MB_ ) + 1 ) * DESCA( MB_ ) + IB
         IB = MIN( IB, M )
         IN = IA+IB-1
         IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
            IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
               DO 10 K = 0, IB-1
                  WRITE( NOUT, FMT = 9999 )
     $                   CMATNM, IA+K, JA+H,
     $                   DBLE( A( II+K+(JJ+H-1)*LDA ) ),
     $                   DIMAG( A( II+K+(JJ+H-1)*LDA ) )
   10          CONTINUE
            END IF
         ELSE
            IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
               CALL ZGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ), LDA,
     $                       IRPRNT, ICPRNT )
            ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
               CALL ZGERV2D( ICTXT, IB, 1, WORK, LDW, ICURROW, ICURCOL )
               DO 20 K = 1, IB
                  WRITE( NOUT, FMT = 9999 )
     $                   CMATNM, IA+K-1, JA+H, DBLE( WORK( K ) ),
     $                   DIMAG( WORK( K ) )
   20          CONTINUE
            END IF
         END IF
         IF( MYROW.EQ.ICURROW )
     $      II = II + IB
         IF( .NOT.AISROWREP )
     $      ICURROW = MOD( ICURROW+1, NPROW )
         CALL BLACS_BARRIER( ICTXT, 'All' )
*
*        Loop over remaining block of rows
*
         DO 50 I = IN+1, IA+M-1, DESCA( MB_ )
            IB = MIN( DESCA( MB_ ), IA+M-I )
            IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
               IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  DO 30 K = 0, IB-1
                     WRITE( NOUT, FMT = 9999 )
     $                      CMATNM, I+K, JA+H,
     $                      DBLE( A( II+K+(JJ+H-1)*LDA ) ),
     $                      DIMAG( A( II+K+(JJ+H-1)*LDA ) )
   30             CONTINUE
               END IF
            ELSE
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                  CALL ZGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                          LDA, IRPRNT, ICPRNT )
               ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  CALL ZGERV2D( ICTXT, IB, 1, WORK, LDW, ICURROW,
     $                          ICURCOL )
                  DO 40 K = 1, IB
                     WRITE( NOUT, FMT = 9999 )
     $                      CMATNM, I+K-1, JA+H, DBLE( WORK( K ) ),
     $                      DIMAG( WORK( K ) )
   40             CONTINUE
               END IF
            END IF
            IF( MYROW.EQ.ICURROW )
     $         II = II + IB
            IF( .NOT.AISROWREP )
     $         ICURROW = MOD( ICURROW+1, NPROW )
            CALL BLACS_BARRIER( ICTXT, 'All' )
   50    CONTINUE
*
         II = IIA
         ICURROW = IAROW
   60 CONTINUE
*
      IF( MYCOL.EQ.ICURCOL )
     $   JJ = JJ + JB
      IF( .NOT.AISCOLREP )
     $   ICURCOL = MOD( ICURCOL+1, NPCOL )
      CALL BLACS_BARRIER( ICTXT, 'All' )
*
*     Loop over remaining column blocks
*
      DO 130 J = JN+1, JA+N-1, DESCA( NB_ )
         JB = MIN(  DESCA( NB_ ), JA+N-J )
         DO 120 H = 0, JB-1
            IB = DESCA( IMB_ )-IA+1
            IF( IB.LE.0 )
     $         IB = ( (-IB) / DESCA( MB_ ) + 1 ) * DESCA( MB_ ) + IB
            IB = MIN( IB, M )
            IN = IA+IB-1
            IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
               IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  DO 70 K = 0, IB-1
                     WRITE( NOUT, FMT = 9999 )
     $                      CMATNM, IA+K, J+H,
     $                      DBLE( A( II+K+(JJ+H-1)*LDA ) ),
     $                      DIMAG( A( II+K+(JJ+H-1)*LDA ) )
   70             CONTINUE
               END IF
            ELSE
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                  CALL ZGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                          LDA, IRPRNT, ICPRNT )
               ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  CALL ZGERV2D( ICTXT, IB, 1, WORK, LDW, ICURROW,
     $                          ICURCOL )
                  DO 80 K = 1, IB
                     WRITE( NOUT, FMT = 9999 )
     $                      CMATNM, IA+K-1, J+H, DBLE( WORK( K ) ),
     $                      DIMAG( WORK( K ) )
   80             CONTINUE
               END IF
            END IF
            IF( MYROW.EQ.ICURROW )
     $         II = II + IB
            ICURROW = MOD( ICURROW+1, NPROW )
            CALL BLACS_BARRIER( ICTXT, 'All' )
*
*           Loop over remaining block of rows
*
            DO 110 I = IN+1, IA+M-1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IA+M-I )
               IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
                  IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                     DO 90 K = 0, IB-1
                        WRITE( NOUT, FMT = 9999 )
     $                         CMATNM, I+K, J+H,
     $                         DBLE( A( II+K+(JJ+H-1)*LDA ) ),
     $                         DIMAG( A( II+K+(JJ+H-1)*LDA ) )
   90                CONTINUE
                  END IF
               ELSE
                  IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                     CALL ZGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                             LDA, IRPRNT, ICPRNT )
                   ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                     CALL ZGERV2D( ICTXT, IB, 1, WORK, LDW, ICURROW,
     $                             ICURCOL )
                     DO 100 K = 1, IB
                        WRITE( NOUT, FMT = 9999 )
     $                         CMATNM, I+K-1, J+H, DBLE( WORK( K ) ),
     $                         DIMAG( WORK( K ) )
  100                CONTINUE
                  END IF
               END IF
               IF( MYROW.EQ.ICURROW )
     $            II = II + IB
               IF( .NOT.AISROWREP )
     $            ICURROW = MOD( ICURROW+1, NPROW )
               CALL BLACS_BARRIER( ICTXT, 'All' )
  110       CONTINUE
*
            II = IIA
            ICURROW = IAROW
  120    CONTINUE
*
         IF( MYCOL.EQ.ICURCOL )
     $      JJ = JJ + JB
         IF( .NOT.AISCOLREP )
     $      ICURCOL = MOD( ICURCOL+1, NPCOL )
         CALL BLACS_BARRIER( ICTXT, 'All' )
*
  130 CONTINUE
*
 9999 FORMAT( 1X, A, '(', I6, ',', I6, ')=', D30.18, '+i*(',
     $        D30.18, ')' )
*
      RETURN
*
*     End of PB_PZLAPRN2
*
      END
      SUBROUTINE PB_ZFILLPAD( ICTXT, M, N, A, LDA, IPRE, IPOST, CHKVAL )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, IPOST, IPRE, LDA, M, N
      COMPLEX*16         CHKVAL
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( * )
*     ..
*
*  Purpose
*  =======
*
*  PB_ZFILLPAD surrounds a two dimensional local array with a guard-zone
*  initialized to the value CHKVAL. The user may later call the  routine
*  PB_ZCHEKPAD to discover if the guardzone has been violated. There are
*  three guardzones. The first is a buffer of size  IPRE  that is before
*  the start of the array. The second is the buffer of size IPOST  which
*  is after the end of the array to be padded. Finally, there is a guard
*  zone inside every column of the array to be padded, in  the  elements
*  of A(M+1:LDA, J).
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  M       (local input) INTEGER
*          On entry, M  specifies the number of rows in the local  array
*          A.  M must be at least zero.
*
*  N       (local input) INTEGER
*          On entry, N  specifies the number of columns in the local ar-
*          ray A. N must be at least zero.
*
*  A       (local input/local output) COMPLEX*16 array
*          On entry,  A  is an array of dimension (LDA,N). On exit, this
*          array is the padded array.
*
*  LDA     (local input) INTEGER
*          On entry,  LDA  specifies  the leading dimension of the local
*          array to be padded. LDA must be at least MAX( 1, M ).
*
*  IPRE    (local input) INTEGER
*          On entry, IPRE specifies the size of  the  guard zone  to put
*          before the start of the padded array.
*
*  IPOST   (local input) INTEGER
*          On entry, IPOST specifies the size of the  guard zone  to put
*          after the end of the padded array.
*
*  CHKVAL  (local input) COMPLEX*16
*          On entry, CHKVAL specifies the value to pad the array with.
*
*  -- Written on April 1, 1998 by
*     R. Clint Whaley, University of Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J, K
*     ..
*     .. Executable Statements ..
*
*     Put check buffer in front of A
*
      IF( IPRE.GT.0 ) THEN
         DO 10 I = 1, IPRE
            A( I ) = CHKVAL
   10    CONTINUE
      ELSE
         WRITE( *, FMT = '(A)' )
     $          'WARNING no pre-guardzone in PB_ZFILLPAD'
      END IF
*
*     Put check buffer in back of A
*
      IF( IPOST.GT.0 ) THEN
         J = IPRE+LDA*N+1
         DO 20 I = J, J+IPOST-1
            A( I ) = CHKVAL
   20    CONTINUE
      ELSE
         WRITE( *, FMT = '(A)' )
     $          'WARNING no post-guardzone in PB_ZFILLPAD'
      END IF
*
*     Put check buffer in all (LDA-M) gaps
*
      IF( LDA.GT.M ) THEN
         K = IPRE + M + 1
         DO 40 J = 1, N
            DO 30 I = K, K + ( LDA - M ) - 1
               A( I ) = CHKVAL
   30       CONTINUE
            K = K + LDA
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of PB_ZFILLPAD
*
      END
      SUBROUTINE PB_ZCHEKPAD( ICTXT, MESS, M, N, A, LDA, IPRE, IPOST,
     $                        CHKVAL )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, IPOST, IPRE, LDA, M, N
      COMPLEX*16         CHKVAL
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      MESS
      COMPLEX*16         A( * )
*     ..
*
*  Purpose
*  =======
*
*  PB_ZCHEKPAD checks that the padding around a local array has not been
*  overwritten since the call to PB_ZFILLPAD.  Three types of errors are
*  reported:
*
*  1) Overwrite in pre-guardzone.  This indicates a memory overwrite has
*  occurred in the  first  IPRE  elements which form a buffer before the
*  beginning of A. Therefore, the error message:
*     'Overwrite in  pre-guardzone: loc(  5) =         18.00000'
*  tells that the 5th element of the IPRE long buffer has been overwrit-
*  ten with the value 18, where it should still have the value CHKVAL.
*
*  2) Overwrite in post-guardzone. This indicates a memory overwrite has
*  occurred in the last IPOST elements which form a buffer after the end
*  of A. Error reports are refered from the end of A.  Therefore,
*     'Overwrite in post-guardzone: loc( 19) =         24.00000'
*  tells  that the  19th element after the end of A was overwritten with
*  the value 24, where it should still have the value of CHKVAL.
*
*  3) Overwrite in lda-m gap.  Tells you elements between M and LDA were
*  overwritten.  So,
*     'Overwrite in lda-m gap: A( 12,  3) =         22.00000'
*  tells  that the element at the 12th row and 3rd column of A was over-
*  written with the value of 22, where it should still have the value of
*  CHKVAL.
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  MESS    (local input) CHARACTER*(*)
*          On entry, MESS is a ttring containing a user-defined message.
*
*  M       (local input) INTEGER
*          On entry, M  specifies the number of rows in the local  array
*          A.  M must be at least zero.
*
*  N       (local input) INTEGER
*          On entry, N  specifies the number of columns in the local ar-
*          ray A. N must be at least zero.
*
*  A       (local input) COMPLEX*16 array
*          On entry,  A  is an array of dimension (LDA,N).
*
*  LDA     (local input) INTEGER
*          On entry,  LDA  specifies  the leading dimension of the local
*          array to be padded. LDA must be at least MAX( 1, M ).
*
*  IPRE    (local input) INTEGER
*          On entry, IPRE specifies the size of  the  guard zone  to put
*          before the start of the padded array.
*
*  IPOST   (local input) INTEGER
*          On entry, IPOST specifies the size of the  guard zone  to put
*          after the end of the padded array.
*
*  CHKVAL  (local input) COMPLEX*16
*          On entry, CHKVAL specifies the value to pad the array with.
*
*
*  -- Written on April 1, 1998 by
*     R. Clint Whaley, University of Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      CHARACTER*1        TOP
      INTEGER            I, IAM, IDUMM, INFO, J, K, MYCOL, MYROW, NPCOL,
     $                   NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGAMX2D, PB_TOPGET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DIMAG
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      IAM  = MYROW*NPCOL + MYCOL
      INFO = -1
*
*     Check buffer in front of A
*
      IF( IPRE.GT.0 ) THEN
         DO 10 I = 1, IPRE
            IF( A( I ).NE.CHKVAL ) THEN
               WRITE( *, FMT = 9998 ) MYROW, MYCOL, MESS, ' pre', I,
     $                                DBLE( A( I ) ), DIMAG( A( I ) )
               INFO = IAM
            END IF
   10    CONTINUE
      ELSE
         WRITE( *, FMT = * ) 'WARNING no pre-guardzone in PB_ZCHEKPAD'
      END IF
*
*     Check buffer after A
*
      IF( IPOST.GT.0 ) THEN
         J = IPRE+LDA*N+1
         DO 20 I = J, J+IPOST-1
            IF( A( I ).NE.CHKVAL ) THEN
               WRITE( *, FMT = 9998 ) MYROW, MYCOL, MESS, 'post',
     $                                I-J+1, DBLE( A( I ) ),
     $                                DIMAG( A( I ) )
               INFO = IAM
            END IF
   20    CONTINUE
      ELSE
         WRITE( *, FMT = * )
     $          'WARNING no post-guardzone buffer in PB_ZCHEKPAD'
      END IF
*
*     Check all (LDA-M) gaps
*
      IF( LDA.GT.M ) THEN
         K = IPRE + M + 1
         DO 40 J = 1, N
            DO 30 I = K, K + (LDA-M) - 1
               IF( A( I ).NE.CHKVAL ) THEN
                  WRITE( *, FMT = 9997 ) MYROW, MYCOL, MESS,
     $               I-IPRE-LDA*(J-1), J, DBLE( A( I ) ),
     $               DIMAG( A( I ) )
                  INFO = IAM
               END IF
   30       CONTINUE
            K = K + LDA
   40    CONTINUE
      END IF
*
      CALL PB_TOPGET( ICTXT, 'Combine', 'All', TOP )
      CALL IGAMX2D( ICTXT, 'All', TOP, 1, 1, INFO, 1, IDUMM, IDUMM, -1,
     $              0, 0 )
      IF( IAM.EQ.0 .AND. INFO.GE.0 ) THEN
         WRITE( *, FMT = 9999 ) INFO / NPCOL, MOD( INFO, NPCOL ), MESS
      END IF
*
 9999 FORMAT( '{', I5, ',', I5, '}:  Memory overwrite in ', A )
 9998 FORMAT( '{', I5, ',', I5, '}:  ', A, ' memory overwrite in ',
     $        A4, '-guardzone: loc(', I3, ') = ', G20.7, '+ i*',
     $        G20.7 )
 9997 FORMAT( '{', I5, ',', I5, '}: ', A, ' memory overwrite in ',
     $        'lda-m gap: loc(', I3, ',', I3, ') = ', G20.7,
     $        '+ i*', G20.7 )
*
      RETURN
*
*     End of PB_ZCHEKPAD
*
      END
      SUBROUTINE PB_ZLASET( UPLO, M, N, IOFFD, ALPHA, BETA, A, LDA )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO
      INTEGER            IOFFD, LDA, M, N
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  PB_ZLASET initializes a two-dimensional array A to beta on the diago-
*  nal specified by IOFFD and alpha on the offdiagonals.
*
*  Arguments
*  =========
*
*  UPLO    (global input) CHARACTER*1
*          On entry,  UPLO  specifies  which trapezoidal part of the ar-
*          ray A is to be set as follows:
*             = 'L' or 'l':   Lower triangular part is set; the strictly
*                             upper triangular part of A is not changed,
*             = 'U' or 'u':   Upper triangular part is set; the strictly
*                             lower triangular part of A is not changed,
*             = 'D' or 'd'    Only the diagonal of A is set,
*             Otherwise:      All of the array A is set.
*
*  M       (input) INTEGER
*          On entry,  M  specifies the number of rows of the array A.  M
*          must be at least zero.
*
*  N       (input) INTEGER
*          On entry,  N  specifies the number of columns of the array A.
*          N must be at least zero.
*
*  IOFFD   (input) INTEGER
*          On entry, IOFFD specifies the position of the offdiagonal de-
*          limiting the upper and lower trapezoidal part of A as follows
*          (see the notes below):
*
*             IOFFD = 0  specifies the main diagonal A( i, i ),
*                        with i = 1 ... MIN( M, N ),
*             IOFFD > 0  specifies the subdiagonal   A( i+IOFFD, i ),
*                        with i = 1 ... MIN( M-IOFFD, N ),
*             IOFFD < 0  specifies the superdiagonal A( i, i-IOFFD ),
*                        with i = 1 ... MIN( M, N+IOFFD ).
*
*  ALPHA   (input) COMPLEX*16
*          On entry,  ALPHA specifies the value to which the offdiagonal
*          array elements are set to.
*
*  BETA    (input) COMPLEX*16
*          On entry, BETA  specifies the value to which the diagonal ar-
*          ray elements are set to.
*
*  A       (input/output) COMPLEX*16 array
*          On entry, A is an array of dimension  (LDA,N).  Before  entry
*          with UPLO = 'U' or 'u', the leading m by n part of the  array
*          A  must  contain  the upper trapezoidal part of the matrix as
*          specified by IOFFD to be set, and  the  strictly lower trape-
*          zoidal  part of A is not referenced; When IUPLO = 'L' or 'l',
*          the leading m by n part of  the  array  A  must  contain  the
*          lower trapezoidal part of the matrix as specified by IOFFD to
*          be set,  and  the  strictly  upper  trapezoidal part of  A is
*          not referenced.
*
*  LDA     (input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least max( 1, M ).
*
*  Notes
*  =====
*                           N                                    N
*             ----------------------------                  -----------
*            |       d                    |                |           |
*          M |         d        'U'       |                |      'U'  |
*            |  'L'     'D'               |                |d          |
*            |             d              |              M |  d        |
*             ----------------------------                 |   'D'     |
*                                                          |      d    |
*               IOFFD < 0                                  | 'L'    d  |
*                                                          |          d|
*                  N                                       |           |
*             -----------                                   -----------
*            |    d   'U'|
*            |      d    |                                   IOFFD > 0
*          M |       'D' |
*            |          d|                              N
*            |  'L'      |                 ----------------------------
*            |           |                |          'U'               |
*            |           |                |d                           |
*            |           |                | 'D'                        |
*            |           |                |    d                       |
*            |           |                |'L'   d                     |
*             -----------                  ----------------------------
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J, JTMP, MN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
*     Start the operations
*
      IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Set the diagonal to BETA and the strictly lower triangular
*        part of the array to ALPHA.
*
         MN = MAX( 0, -IOFFD )
         DO 20 J = 1, MIN( MN, N )
            DO 10 I = 1, M
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
         DO 40 J = MN + 1, MIN( M - IOFFD, N )
            JTMP = J + IOFFD
            A( JTMP, J ) = BETA
            DO 30 I = JTMP + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Set the diagonal to BETA and the strictly upper triangular
*        part of the array to ALPHA.
*
         MN = MIN( M - IOFFD, N )
         DO 60 J = MAX( 0, -IOFFD ) + 1, MN
            JTMP = J + IOFFD
            DO 50 I = 1, JTMP - 1
               A( I, J ) = ALPHA
   50       CONTINUE
            A( JTMP, J ) = BETA
   60    CONTINUE
         DO 80 J = MAX( 0, MN ) + 1, N
            DO 70 I = 1, M
               A( I, J ) = ALPHA
   70       CONTINUE
   80    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'D' ) ) THEN
*
*        Set the array to BETA on the diagonal.
*
         DO 90 J = MAX( 0, -IOFFD ) + 1, MIN( M - IOFFD, N )
            A( J + IOFFD, J ) = BETA
   90    CONTINUE
*
      ELSE
*
*        Set the array to BETA on the diagonal and ALPHA on the
*        offdiagonal.
*
         DO 110 J = 1, N
            DO 100 I = 1, M
               A( I, J ) = ALPHA
  100       CONTINUE
  110    CONTINUE
         IF( ALPHA.NE.BETA .AND. IOFFD.LT.M .AND. IOFFD.GT.-N ) THEN
            DO 120 J = MAX( 0, -IOFFD ) + 1, MIN( M - IOFFD, N )
               A( J + IOFFD, J ) = BETA
  120       CONTINUE
         END IF
*
      END IF
*
      RETURN
*
*     End of PB_ZLASET
*
      END
      SUBROUTINE PB_ZLASCAL( UPLO, M, N, IOFFD, ALPHA, A, LDA )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO
      INTEGER            IOFFD, LDA, M, N
      COMPLEX*16         ALPHA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  PB_ZLASCAL scales a two-dimensional array A by the scalar alpha.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          On entry,  UPLO  specifies  which trapezoidal part of the ar-
*          ray A is to be scaled as follows:
*             = 'L' or 'l':          the lower trapezoid of A is scaled,
*             = 'U' or 'u':          the upper trapezoid of A is scaled,
*             = 'D' or 'd':       diagonal specified by IOFFD is scaled,
*             Otherwise:                   all of the array A is scaled.
*
*  M       (input) INTEGER
*          On entry,  M  specifies the number of rows of the array A.  M
*          must be at least zero.
*
*  N       (input) INTEGER
*          On entry,  N  specifies the number of columns of the array A.
*          N must be at least zero.
*
*  IOFFD   (input) INTEGER
*          On entry, IOFFD specifies the position of the offdiagonal de-
*          limiting the upper and lower trapezoidal part of A as follows
*          (see the notes below):
*
*             IOFFD = 0  specifies the main diagonal A( i, i ),
*                        with i = 1 ... MIN( M, N ),
*             IOFFD > 0  specifies the subdiagonal   A( i+IOFFD, i ),
*                        with i = 1 ... MIN( M-IOFFD, N ),
*             IOFFD < 0  specifies the superdiagonal A( i, i-IOFFD ),
*                        with i = 1 ... MIN( M, N+IOFFD ).
*
*  ALPHA   (input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (input/output) COMPLEX*16 array
*          On entry, A is an array of dimension  (LDA,N).  Before  entry
*          with  UPLO = 'U' or 'u', the leading m by n part of the array
*          A must contain the upper trapezoidal  part  of the matrix  as
*          specified by  IOFFD to be scaled, and the strictly lower tra-
*          pezoidal part of A is not referenced; When UPLO = 'L' or 'l',
*          the leading m by n part of the array A must contain the lower
*          trapezoidal  part  of  the matrix as specified by IOFFD to be
*          scaled,  and  the strictly upper trapezoidal part of A is not
*          referenced. On exit, the entries of the  trapezoid part of  A
*          determined by UPLO and IOFFD are scaled.
*
*  LDA     (input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least max( 1, M ).
*
*  Notes
*  =====
*                           N                                    N
*             ----------------------------                  -----------
*            |       d                    |                |           |
*          M |         d        'U'       |                |      'U'  |
*            |  'L'     'D'               |                |d          |
*            |             d              |              M |  d        |
*             ----------------------------                 |   'D'     |
*                                                          |      d    |
*              IOFFD < 0                                   | 'L'    d  |
*                                                          |          d|
*                  N                                       |           |
*             -----------                                   -----------
*            |    d   'U'|
*            |      d    |                                   IOFFD > 0
*          M |       'D' |
*            |          d|                              N
*            |  'L'      |                 ----------------------------
*            |           |                |          'U'               |
*            |           |                |d                           |
*            |           |                | 'D'                        |
*            |           |                |    d                       |
*            |           |                |'L'   d                     |
*             -----------                  ----------------------------
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J, JTMP, MN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
*     Start the operations
*
      IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Scales the lower triangular part of the array by ALPHA.
*
         MN = MAX( 0, -IOFFD )
         DO 20 J = 1, MIN( MN, N )
            DO 10 I = 1, M
               A( I, J ) = ALPHA * A( I, J )
   10       CONTINUE
   20    CONTINUE
         DO 40 J = MN + 1, MIN( M - IOFFD, N )
            DO 30 I = J + IOFFD, M
               A( I, J ) = ALPHA * A( I, J )
   30       CONTINUE
   40    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Scales the upper triangular part of the array by ALPHA.
*
         MN = MIN( M - IOFFD, N )
         DO 60 J = MAX( 0, -IOFFD ) + 1, MN
            DO 50 I = 1, J + IOFFD
               A( I, J ) = ALPHA * A( I, J )
   50       CONTINUE
   60    CONTINUE
         DO 80 J = MAX( 0, MN ) + 1, N
            DO 70 I = 1, M
               A( I, J ) = ALPHA * A( I, J )
   70       CONTINUE
   80    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'D' ) ) THEN
*
*        Scales the diagonal entries by ALPHA.
*
         DO 90 J = MAX( 0, -IOFFD ) + 1, MIN( M - IOFFD, N )
            JTMP = J + IOFFD
            A( JTMP, J ) = ALPHA * A( JTMP, J )
   90    CONTINUE
*
      ELSE
*
*        Scales the entire array by ALPHA.
*
         DO 110 J = 1, N
            DO 100 I = 1, M
               A( I, J ) = ALPHA * A( I, J )
  100       CONTINUE
  110    CONTINUE
*
      END IF
*
      RETURN
*
*     End of PB_ZLASCAL
*
      END
      SUBROUTINE PB_ZLAGEN( UPLO, AFORM, A, LDA, LCMT00, IRAN, MBLKS,
     $                      IMBLOC, MB, LMBLOC, NBLKS, INBLOC, NB,
     $                      LNBLOC, JMP, IMULADD )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO, AFORM
      INTEGER            IMBLOC, INBLOC, LCMT00, LDA, LMBLOC, LNBLOC,
     $                   MB, MBLKS, NB, NBLKS
*     ..
*     .. Array Arguments ..
      INTEGER            IMULADD( 4, * ), IRAN( * ), JMP( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  PB_ZLAGEN locally initializes an array A.
*
*  Arguments
*  =========
*
*  UPLO    (global input) CHARACTER*1
*          On entry, UPLO  specifies whether the lower (UPLO='L') trape-
*          zoidal part or the upper (UPLO='U') trapezoidal part is to be
*          generated  when  the  matrix  to be generated is symmetric or
*          Hermitian. For  all  the  other values of AFORM, the value of
*          this input argument is ignored.
*
*  AFORM   (global input) CHARACTER*1
*          On entry, AFORM specifies the type of submatrix to be genera-
*          ted as follows:
*             AFORM = 'S', sub( A ) is a symmetric matrix,
*             AFORM = 'H', sub( A ) is a Hermitian matrix,
*             AFORM = 'T', sub( A ) is overrwritten  with  the transpose
*                          of what would normally be generated,
*             AFORM = 'C', sub( A ) is overwritten  with  the  conjugate
*                          transpose  of  what would normally be genera-
*                          ted.
*             AFORM = 'N', a random submatrix is generated.
*
*  A       (local output) COMPLEX*16 array
*          On entry,  A  is  an array of dimension (LLD_A, *).  On exit,
*          this array contains the local entries of the randomly genera-
*          ted submatrix sub( A ).
*
*  LDA     (local input) INTEGER
*          On entry,  LDA  specifies  the local leading dimension of the
*          array A. LDA must be at least one.
*
*  LCMT00  (global input) INTEGER
*          On entry, LCMT00 is the LCM value specifying the off-diagonal
*          of the underlying matrix of interest. LCMT00=0 specifies  the
*          main diagonal, LCMT00 > 0 specifies a subdiagonal, LCMT00 < 0
*          specifies superdiagonals.
*
*  IRAN    (local input) INTEGER array
*          On entry, IRAN  is an array of dimension 2 containing respec-
*          tively the 16-lower and 16-higher bits of the encoding of the
*          entry of  the  random  sequence  corresponding locally to the
*          first local array entry to generate. Usually,  this  array is
*          computed by PB_SETLOCRAN.
*
*  MBLKS   (local input) INTEGER
*          On entry, MBLKS specifies the local number of blocks of rows.
*          MBLKS is at least zero.
*
*  IMBLOC  (local input) INTEGER
*          On entry, IMBLOC specifies  the  number of rows (size) of the
*          local uppest  blocks. IMBLOC is at least zero.
*
*  MB      (global input) INTEGER
*          On entry, MB  specifies the blocking factor used to partition
*          the rows of the matrix.  MB  must be at least one.
*
*  LMBLOC  (local input) INTEGER
*          On entry, LMBLOC specifies the number of  rows  (size) of the
*          local lowest blocks. LMBLOC is at least zero.
*
*  NBLKS   (local input) INTEGER
*          On entry,  NBLKS  specifies the local number of blocks of co-
*          lumns. NBLKS is at least zero.
*
*  INBLOC  (local input) INTEGER
*          On entry,  INBLOC  specifies the number of columns (size)  of
*          the local leftmost blocks. INBLOC is at least zero.
*
*  NB      (global input) INTEGER
*          On entry, NB  specifies the blocking factor used to partition
*          the the columns of the matrix.  NB  must be at least one.
*
*  LNBLOC  (local input) INTEGER
*          On entry,  LNBLOC  specifies  the number of columns (size) of
*          the local rightmost blocks. LNBLOC is at least zero.
*
*  JMP     (local input) INTEGER array
*          On entry, JMP is an array of dimension JMP_LEN containing the
*          different jump values used by the random matrix generator.
*
*  IMULADD (local input) INTEGER array
*          On entry, IMULADD is an array of dimension (4, JMP_LEN).  The
*          jth  column  of this array contains the encoded initial cons-
*          tants a_j and c_j to  jump  from X( n ) to  X( n + JMP( j ) )
*          (= a_j * X( n ) + c_j) in the random sequence. IMULADD(1:2,j)
*          contains respectively the 16-lower and 16-higher bits of  the
*          constant a_j, and IMULADD(3:4,j)  contains  the 16-lower  and
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
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IB, IBLK, II, IK, ITMP, JB, JBLK, JJ, JK,
     $                   JTMP, LCMTC, LCMTR, LOW, MNB, UPP
      COMPLEX*16         DUMMY
*     ..
*     .. Local Arrays ..
      INTEGER            IB0( 2 ), IB1( 2 ), IB2( 2 ), IB3( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           PB_JUMPIT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   PB_DRAND
      EXTERNAL           LSAME, PB_DRAND
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      DO 10 I = 1, 2
         IB1( I ) = IRAN( I )
         IB2( I ) = IRAN( I )
         IB3( I ) = IRAN( I )
   10 CONTINUE
*
      IF( LSAME( AFORM, 'N' ) ) THEN
*
*        Generate random matrix
*
         JJ = 1
*
         DO 50 JBLK = 1, NBLKS
*
            IF( JBLK.EQ.1 ) THEN
               JB = INBLOC
            ELSE IF( JBLK.EQ.NBLKS ) THEN
               JB = LNBLOC
            ELSE
               JB = NB
            END IF
*
            DO 40 JK = JJ, JJ + JB - 1
*
               II = 1
*
               DO 30 IBLK = 1, MBLKS
*
                  IF( IBLK.EQ.1 ) THEN
                     IB = IMBLOC
                  ELSE IF( IBLK.EQ.MBLKS ) THEN
                     IB = LMBLOC
                  ELSE
                     IB = MB
                  END IF
*
*                 Blocks are IB by JB
*
                  DO 20 IK = II, II + IB - 1
                     A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                     PB_DRAND( 0 ) )
   20             CONTINUE
*
                  II = II + IB
*
                  IF( IBLK.EQ.1 ) THEN
*
*                    Jump IMBLOC + ( NPROW - 1 ) * MB rows
*
                     CALL PB_JUMPIT( IMULADD( 1, JMP_NPIMBLOC ), IB1,
     $                               IB0 )
*
                  ELSE
*
*                    Jump NPROW * MB rows
*
                     CALL PB_JUMPIT( IMULADD( 1, JMP_NPMB ), IB1, IB0 )
*
                  END IF
*
                  IB1( 1 ) = IB0( 1 )
                  IB1( 2 ) = IB0( 2 )
*
   30          CONTINUE
*
*              Jump one column
*
               CALL PB_JUMPIT( IMULADD( 1, JMP_COL ), IB2, IB0 )
*
               IB1( 1 ) = IB0( 1 )
               IB1( 2 ) = IB0( 2 )
               IB2( 1 ) = IB0( 1 )
               IB2( 2 ) = IB0( 2 )
*
   40       CONTINUE
*
            JJ = JJ + JB
*
            IF( JBLK.EQ.1 ) THEN
*
*              Jump INBLOC + ( NPCOL - 1 ) * NB columns
*
               CALL PB_JUMPIT( IMULADD( 1, JMP_NQINBLOC ), IB3, IB0 )
*
            ELSE
*
*              Jump NPCOL * NB columns
*
               CALL PB_JUMPIT( IMULADD( 1, JMP_NQNB ), IB3, IB0 )
*
            END IF
*
            IB1( 1 ) = IB0( 1 )
            IB1( 2 ) = IB0( 2 )
            IB2( 1 ) = IB0( 1 )
            IB2( 2 ) = IB0( 2 )
            IB3( 1 ) = IB0( 1 )
            IB3( 2 ) = IB0( 2 )
*
   50    CONTINUE
*
      ELSE IF( LSAME( AFORM, 'T' ) ) THEN
*
*        Generate the transpose of the matrix that would be normally
*        generated.
*
         II = 1
*
         DO 90 IBLK = 1, MBLKS
*
            IF( IBLK.EQ.1 ) THEN
               IB = IMBLOC
            ELSE IF( IBLK.EQ.MBLKS ) THEN
               IB = LMBLOC
            ELSE
               IB = MB
            END IF
*
            DO 80 IK = II, II + IB - 1
*
               JJ = 1
*
               DO 70 JBLK = 1, NBLKS
*
                  IF( JBLK.EQ.1 ) THEN
                     JB = INBLOC
                  ELSE IF( JBLK.EQ.NBLKS ) THEN
                     JB = LNBLOC
                  ELSE
                     JB = NB
                  END IF
*
*                 Blocks are IB by JB
*
                  DO 60 JK = JJ, JJ + JB - 1
                     A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                     PB_DRAND( 0 ) )
   60             CONTINUE
*
                  JJ = JJ + JB
*
                  IF( JBLK.EQ.1 ) THEN
*
*                    Jump INBLOC + ( NPCOL - 1 ) * NB columns
*
                     CALL PB_JUMPIT( IMULADD( 1, JMP_NQINBLOC ), IB1,
     $                               IB0 )
*
                  ELSE
*
*                    Jump NPCOL * NB columns
*
                     CALL PB_JUMPIT( IMULADD( 1, JMP_NQNB ), IB1, IB0 )
*
                  END IF
*
                  IB1( 1 ) = IB0( 1 )
                  IB1( 2 ) = IB0( 2 )
*
   70          CONTINUE
*
*              Jump one row
*
               CALL PB_JUMPIT( IMULADD( 1, JMP_ROW ), IB2, IB0 )
*
               IB1( 1 ) = IB0( 1 )
               IB1( 2 ) = IB0( 2 )
               IB2( 1 ) = IB0( 1 )
               IB2( 2 ) = IB0( 2 )
*
   80       CONTINUE
*
            II = II + IB
*
            IF( IBLK.EQ.1 ) THEN
*
*              Jump IMBLOC + ( NPROW - 1 ) * MB rows
*
               CALL PB_JUMPIT( IMULADD( 1, JMP_NPIMBLOC ), IB3, IB0 )
*
            ELSE
*
*              Jump NPROW * MB rows
*
               CALL PB_JUMPIT( IMULADD( 1, JMP_NPMB ), IB3, IB0 )
*
            END IF
*
            IB1( 1 ) = IB0( 1 )
            IB1( 2 ) = IB0( 2 )
            IB2( 1 ) = IB0( 1 )
            IB2( 2 ) = IB0( 2 )
            IB3( 1 ) = IB0( 1 )
            IB3( 2 ) = IB0( 2 )
*
   90    CONTINUE
*
      ELSE IF( LSAME( AFORM, 'S' ) ) THEN
*
*        Generate a symmetric matrix
*
         IF( LSAME( UPLO, 'L' ) ) THEN
*
*           generate lower trapezoidal part
*
            JJ = 1
            LCMTC = LCMT00
*
            DO 170 JBLK = 1, NBLKS
*
               IF( JBLK.EQ.1 ) THEN
                  JB  = INBLOC
                  LOW = 1 - INBLOC
               ELSE IF( JBLK.EQ.NBLKS ) THEN
                  JB = LNBLOC
                  LOW = 1 - NB
               ELSE
                  JB  = NB
                  LOW = 1 - NB
               END IF
*
               DO 160 JK = JJ, JJ + JB - 1
*
                  II = 1
                  LCMTR = LCMTC
*
                  DO 150 IBLK = 1, MBLKS
*
                     IF( IBLK.EQ.1 ) THEN
                        IB  = IMBLOC
                        UPP = IMBLOC - 1
                     ELSE IF( IBLK.EQ.MBLKS ) THEN
                        IB  = LMBLOC
                        UPP = MB - 1
                     ELSE
                        IB  = MB
                        UPP = MB - 1
                     END IF
*
*                    Blocks are IB by JB
*
                     IF( LCMTR.GT.UPP ) THEN
*
                        DO 100 IK = II, II + IB - 1
                           DUMMY = DCMPLX( PB_DRAND( 0 ),
     $                                     PB_DRAND( 0 ) )
  100                   CONTINUE
*
                     ELSE IF( LCMTR.GE.LOW ) THEN
*
                        JTMP = JK - JJ + 1
                        MNB  = MAX( 0, -LCMTR )
*
                        IF( JTMP.LE.MIN( MNB, JB ) ) THEN
*
                           DO 110 IK = II, II + IB - 1
                              A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                              PB_DRAND( 0 ) )
  110                      CONTINUE
*
                        ELSE IF( ( JTMP.GE.( MNB + 1 )         ) .AND.
     $                           ( JTMP.LE.MIN( IB-LCMTR, JB ) ) ) THEN
*
                           ITMP = II + JTMP + LCMTR - 1
*
                           DO 120 IK = II, ITMP - 1
                              DUMMY = DCMPLX( PB_DRAND( 0 ),
     $                                        PB_DRAND( 0 ) )
  120                      CONTINUE
*
                           DO 130 IK = ITMP, II + IB - 1
                              A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                              PB_DRAND( 0 ) )
  130                      CONTINUE
*
                        END IF
*
                     ELSE
*
                        DO 140 IK = II, II + IB - 1
                           A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                           PB_DRAND( 0 ) )
  140                   CONTINUE
*
                     END IF
*
                     II = II + IB
*
                     IF( IBLK.EQ.1 ) THEN
*
*                       Jump IMBLOC + ( NPROW - 1 ) * MB rows
*
                        LCMTR = LCMTR - JMP( JMP_NPIMBLOC )
                        CALL PB_JUMPIT( IMULADD( 1, JMP_NPIMBLOC ), IB1,
     $                                  IB0 )
*
                     ELSE
*
*                       Jump NPROW * MB rows
*
                        LCMTR = LCMTR - JMP( JMP_NPMB )
                        CALL PB_JUMPIT( IMULADD( 1, JMP_NPMB ), IB1,
     $                                  IB0 )
*
                     END IF
*
                     IB1( 1 ) = IB0( 1 )
                     IB1( 2 ) = IB0( 2 )
*
  150             CONTINUE
*
*                 Jump one column
*
                  CALL PB_JUMPIT( IMULADD( 1, JMP_COL ), IB2, IB0 )
*
                  IB1( 1 ) = IB0( 1 )
                  IB1( 2 ) = IB0( 2 )
                  IB2( 1 ) = IB0( 1 )
                  IB2( 2 ) = IB0( 2 )
*
  160          CONTINUE
*
               JJ = JJ + JB
*
               IF( JBLK.EQ.1 ) THEN
*
*                 Jump INBLOC + ( NPCOL - 1 ) * NB columns
*
                  LCMTC = LCMTC + JMP( JMP_NQINBLOC )
                  CALL PB_JUMPIT( IMULADD( 1, JMP_NQINBLOC ), IB3, IB0 )
*
               ELSE
*
*                 Jump NPCOL * NB columns
*
                  LCMTC = LCMTC + JMP( JMP_NQNB )
                  CALL PB_JUMPIT( IMULADD( 1, JMP_NQNB ), IB3, IB0 )
*
               END IF
*
               IB1( 1 ) = IB0( 1 )
               IB1( 2 ) = IB0( 2 )
               IB2( 1 ) = IB0( 1 )
               IB2( 2 ) = IB0( 2 )
               IB3( 1 ) = IB0( 1 )
               IB3( 2 ) = IB0( 2 )
*
  170       CONTINUE
*
         ELSE
*
*           generate upper trapezoidal part
*
            II = 1
            LCMTR = LCMT00
*
            DO 250 IBLK = 1, MBLKS
*
               IF( IBLK.EQ.1 ) THEN
                  IB  = IMBLOC
                  UPP = IMBLOC - 1
               ELSE IF( IBLK.EQ.MBLKS ) THEN
                  IB  = LMBLOC
                  UPP = MB - 1
               ELSE
                  IB  = MB
                  UPP = MB - 1
               END IF
*
               DO 240 IK = II, II + IB - 1
*
                  JJ = 1
                  LCMTC = LCMTR
*
                  DO 230 JBLK = 1, NBLKS
*
                     IF( JBLK.EQ.1 ) THEN
                        JB  = INBLOC
                        LOW = 1 - INBLOC
                     ELSE IF( JBLK.EQ.NBLKS ) THEN
                        JB  = LNBLOC
                        LOW = 1 - NB
                     ELSE
                        JB  = NB
                        LOW = 1 - NB
                     END IF
*
*                    Blocks are IB by JB
*
                     IF( LCMTC.LT.LOW ) THEN
*
                        DO 180 JK = JJ, JJ + JB - 1
                           DUMMY = DCMPLX( PB_DRAND( 0 ),
     $                                     PB_DRAND( 0 ) )
  180                   CONTINUE
*
                     ELSE IF( LCMTC.LE.UPP ) THEN
*
                        ITMP = IK - II + 1
                        MNB  = MAX( 0, LCMTC )
*
                        IF( ITMP.LE.MIN( MNB, IB ) ) THEN
*
                           DO 190 JK = JJ, JJ + JB - 1
                              A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                              PB_DRAND( 0 ) )
  190                      CONTINUE
*
                        ELSE IF( ( ITMP.GE.( MNB + 1 )         ) .AND.
     $                           ( ITMP.LE.MIN( JB+LCMTC, IB ) ) ) THEN
*
                           JTMP = JJ + ITMP - LCMTC - 1
*
                           DO 200 JK = JJ, JTMP - 1
                              DUMMY = DCMPLX( PB_DRAND( 0 ),
     $                                        PB_DRAND( 0 ) )
  200                      CONTINUE
*
                           DO 210 JK = JTMP, JJ + JB - 1
                              A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                              PB_DRAND( 0 ) )
  210                      CONTINUE
*
                        END IF
*
                     ELSE
*
                        DO 220 JK = JJ, JJ + JB - 1
                           A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                           PB_DRAND( 0 ) )
  220                   CONTINUE
*
                     END IF
*
                     JJ = JJ + JB
*
                     IF( JBLK.EQ.1 ) THEN
*
*                       Jump INBLOC + ( NPCOL - 1 ) * NB columns
*
                        LCMTC = LCMTC + JMP( JMP_NQINBLOC )
                        CALL PB_JUMPIT( IMULADD( 1, JMP_NQINBLOC ), IB1,
     $                                  IB0 )
*
                     ELSE
*
*                       Jump NPCOL * NB columns
*
                        LCMTC = LCMTC + JMP( JMP_NQNB )
                        CALL PB_JUMPIT( IMULADD( 1, JMP_NQNB ), IB1,
     $                                  IB0 )
*
                     END IF
*
                     IB1( 1 ) = IB0( 1 )
                     IB1( 2 ) = IB0( 2 )
*
  230             CONTINUE
*
*                 Jump one row
*
                  CALL PB_JUMPIT( IMULADD( 1, JMP_ROW ), IB2, IB0 )
*
                  IB1( 1 ) = IB0( 1 )
                  IB1( 2 ) = IB0( 2 )
                  IB2( 1 ) = IB0( 1 )
                  IB2( 2 ) = IB0( 2 )
*
  240          CONTINUE
*
               II = II + IB
*
               IF( IBLK.EQ.1 ) THEN
*
*                 Jump IMBLOC + ( NPROW - 1 ) * MB rows
*
                  LCMTR = LCMTR - JMP( JMP_NPIMBLOC )
                  CALL PB_JUMPIT( IMULADD( 1, JMP_NPIMBLOC ), IB3, IB0 )
*
               ELSE
*
*                 Jump NPROW * MB rows
*
                  LCMTR = LCMTR - JMP( JMP_NPMB )
                  CALL PB_JUMPIT( IMULADD( 1, JMP_NPMB ), IB3, IB0 )
*
               END IF
*
               IB1( 1 ) = IB0( 1 )
               IB1( 2 ) = IB0( 2 )
               IB2( 1 ) = IB0( 1 )
               IB2( 2 ) = IB0( 2 )
               IB3( 1 ) = IB0( 1 )
               IB3( 2 ) = IB0( 2 )
*
  250       CONTINUE
*
         END IF
*
      ELSE IF( LSAME( AFORM, 'C' ) ) THEN
*
*        Generate the conjugate transpose of the matrix that would be
*        normally generated.
*
         II = 1
*
         DO 290 IBLK = 1, MBLKS
*
            IF( IBLK.EQ.1 ) THEN
               IB = IMBLOC
            ELSE IF( IBLK.EQ.MBLKS ) THEN
               IB = LMBLOC
            ELSE
               IB = MB
            END IF
*
            DO 280 IK = II, II + IB - 1
*
               JJ = 1
*
               DO 270 JBLK = 1, NBLKS
*
                  IF( JBLK.EQ.1 ) THEN
                     JB = INBLOC
                  ELSE IF( JBLK.EQ.NBLKS ) THEN
                     JB = LNBLOC
                  ELSE
                     JB = NB
                  END IF
*
*                 Blocks are IB by JB
*
                  DO 260 JK = JJ, JJ + JB - 1
                     A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                    -PB_DRAND( 0 ) )
  260             CONTINUE
*
                  JJ = JJ + JB
*
                  IF( JBLK.EQ.1 ) THEN
*
*                    Jump INBLOC + ( NPCOL - 1 ) * NB columns
*
                     CALL PB_JUMPIT( IMULADD( 1, JMP_NQINBLOC ), IB1,
     $                               IB0 )
*
                  ELSE
*
*                    Jump NPCOL * NB columns
*
                     CALL PB_JUMPIT( IMULADD( 1, JMP_NQNB ), IB1,
     $                               IB0 )
*
                  END IF
*
                  IB1( 1 ) = IB0( 1 )
                  IB1( 2 ) = IB0( 2 )
*
  270          CONTINUE
*
*              Jump one row
*
               CALL PB_JUMPIT( IMULADD( 1, JMP_ROW ), IB2, IB0 )
*
               IB1( 1 ) = IB0( 1 )
               IB1( 2 ) = IB0( 2 )
               IB2( 1 ) = IB0( 1 )
               IB2( 2 ) = IB0( 2 )
*
  280       CONTINUE
*
            II = II + IB
*
            IF( IBLK.EQ.1 ) THEN
*
*              Jump IMBLOC + ( NPROW - 1 ) * MB rows
*
               CALL PB_JUMPIT( IMULADD( 1, JMP_NPIMBLOC ), IB3, IB0 )
*
            ELSE
*
*              Jump NPROW * MB rows
*
               CALL PB_JUMPIT( IMULADD( 1, JMP_NPMB ), IB3, IB0 )
*
            END IF
*
            IB1( 1 ) = IB0( 1 )
            IB1( 2 ) = IB0( 2 )
            IB2( 1 ) = IB0( 1 )
            IB2( 2 ) = IB0( 2 )
            IB3( 1 ) = IB0( 1 )
            IB3( 2 ) = IB0( 2 )
*
  290    CONTINUE
*
      ELSE IF( LSAME( AFORM, 'H' ) ) THEN
*
*        Generate a Hermitian matrix
*
         IF( LSAME( UPLO, 'L' ) ) THEN
*
*           generate lower trapezoidal part
*
            JJ = 1
            LCMTC = LCMT00
*
            DO 370 JBLK = 1, NBLKS
*
               IF( JBLK.EQ.1 ) THEN
                  JB  = INBLOC
                  LOW = 1 - INBLOC
               ELSE IF( JBLK.EQ.NBLKS ) THEN
                  JB = LNBLOC
                  LOW = 1 - NB
               ELSE
                  JB  = NB
                  LOW = 1 - NB
               END IF
*
               DO 360 JK = JJ, JJ + JB - 1
*
                  II = 1
                  LCMTR = LCMTC
*
                  DO 350 IBLK = 1, MBLKS
*
                     IF( IBLK.EQ.1 ) THEN
                        IB  = IMBLOC
                        UPP = IMBLOC - 1
                     ELSE IF( IBLK.EQ.MBLKS ) THEN
                        IB  = LMBLOC
                        UPP = MB - 1
                     ELSE
                        IB  = MB
                        UPP = MB - 1
                     END IF
*
*                    Blocks are IB by JB
*
                     IF( LCMTR.GT.UPP ) THEN
*
                        DO 300 IK = II, II + IB - 1
                           DUMMY = DCMPLX( PB_DRAND( 0 ),
     $                                     PB_DRAND( 0 ) )
  300                   CONTINUE
*
                     ELSE IF( LCMTR.GE.LOW ) THEN
*
                        JTMP = JK - JJ + 1
                        MNB  = MAX( 0, -LCMTR )
*
                        IF( JTMP.LE.MIN( MNB, JB ) ) THEN
*
                           DO 310 IK = II, II + IB - 1
                              A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                              PB_DRAND( 0 ) )
  310                      CONTINUE
*
                        ELSE IF( ( JTMP.GE.( MNB + 1 )         ) .AND.
     $                           ( JTMP.LE.MIN( IB-LCMTR, JB ) ) ) THEN
*
                           ITMP = II + JTMP + LCMTR - 1
*
                           DO 320 IK = II, ITMP - 1
                              DUMMY = DCMPLX( PB_DRAND( 0 ),
     $                                        PB_DRAND( 0 ) )
  320                      CONTINUE
*
                           IF( ITMP.LE.( II + IB - 1 ) ) THEN
                              DUMMY = DCMPLX( PB_DRAND( 0 ),
     $                                       -PB_DRAND( 0 ) )
                              A( ITMP, JK ) = DCMPLX( DBLE( DUMMY ),
     $                                                ZERO )
                           END IF
*
                           DO 330 IK = ITMP + 1, II + IB - 1
                              A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                              PB_DRAND( 0 ) )
  330                      CONTINUE
*
                        END IF
*
                     ELSE
*
                        DO 340 IK = II, II + IB - 1
                           A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                           PB_DRAND( 0 ) )
  340                   CONTINUE
*
                     END IF
*
                     II = II + IB
*
                     IF( IBLK.EQ.1 ) THEN
*
*                       Jump IMBLOC + ( NPROW - 1 ) * MB rows
*
                        LCMTR = LCMTR - JMP( JMP_NPIMBLOC )
                        CALL PB_JUMPIT( IMULADD( 1, JMP_NPIMBLOC ), IB1,
     $                                  IB0 )
*
                     ELSE
*
*                       Jump NPROW * MB rows
*
                        LCMTR = LCMTR - JMP( JMP_NPMB )
                        CALL PB_JUMPIT( IMULADD( 1, JMP_NPMB ), IB1,
     $                                  IB0 )
*
                     END IF
*
                     IB1( 1 ) = IB0( 1 )
                     IB1( 2 ) = IB0( 2 )
*
  350             CONTINUE
*
*                 Jump one column
*
                  CALL PB_JUMPIT( IMULADD( 1, JMP_COL ), IB2, IB0 )
*
                  IB1( 1 ) = IB0( 1 )
                  IB1( 2 ) = IB0( 2 )
                  IB2( 1 ) = IB0( 1 )
                  IB2( 2 ) = IB0( 2 )
*
  360          CONTINUE
*
               JJ = JJ + JB
*
               IF( JBLK.EQ.1 ) THEN
*
*                 Jump INBLOC + ( NPCOL - 1 ) * NB columns
*
                  LCMTC = LCMTC + JMP( JMP_NQINBLOC )
                  CALL PB_JUMPIT( IMULADD( 1, JMP_NQINBLOC ), IB3, IB0 )
*
               ELSE
*
*                 Jump NPCOL * NB columns
*
                  LCMTC = LCMTC + JMP( JMP_NQNB )
                  CALL PB_JUMPIT( IMULADD( 1, JMP_NQNB ), IB3, IB0 )
*
               END IF
*
               IB1( 1 ) = IB0( 1 )
               IB1( 2 ) = IB0( 2 )
               IB2( 1 ) = IB0( 1 )
               IB2( 2 ) = IB0( 2 )
               IB3( 1 ) = IB0( 1 )
               IB3( 2 ) = IB0( 2 )
*
  370       CONTINUE
*
         ELSE
*
*           generate upper trapezoidal part
*
            II = 1
            LCMTR = LCMT00
*
            DO 450 IBLK = 1, MBLKS
*
               IF( IBLK.EQ.1 ) THEN
                  IB  = IMBLOC
                  UPP = IMBLOC - 1
               ELSE IF( IBLK.EQ.MBLKS ) THEN
                  IB  = LMBLOC
                  UPP = MB - 1
               ELSE
                  IB  = MB
                  UPP = MB - 1
               END IF
*
               DO 440 IK = II, II + IB - 1
*
                  JJ = 1
                  LCMTC = LCMTR
*
                  DO 430 JBLK = 1, NBLKS
*
                     IF( JBLK.EQ.1 ) THEN
                        JB  = INBLOC
                        LOW = 1 - INBLOC
                     ELSE IF( JBLK.EQ.NBLKS ) THEN
                        JB  = LNBLOC
                        LOW = 1 - NB
                     ELSE
                        JB  = NB
                        LOW = 1 - NB
                     END IF
*
*                    Blocks are IB by JB
*
                     IF( LCMTC.LT.LOW ) THEN
*
                        DO 380 JK = JJ, JJ + JB - 1
                           DUMMY = DCMPLX( PB_DRAND( 0 ),
     $                                    -PB_DRAND( 0 ) )
  380                   CONTINUE
*
                     ELSE IF( LCMTC.LE.UPP ) THEN
*
                        ITMP = IK - II + 1
                        MNB  = MAX( 0, LCMTC )
*
                        IF( ITMP.LE.MIN( MNB, IB ) ) THEN
*
                           DO 390 JK = JJ, JJ + JB - 1
                              A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                             -PB_DRAND( 0 ) )
  390                      CONTINUE
*
                        ELSE IF( ( ITMP.GE.( MNB + 1 )         ) .AND.
     $                           ( ITMP.LE.MIN( JB+LCMTC, IB ) ) ) THEN
*
                           JTMP = JJ + ITMP - LCMTC - 1
*
                           DO 400 JK = JJ, JTMP - 1
                              DUMMY = DCMPLX( PB_DRAND( 0 ),
     $                                       -PB_DRAND( 0 ) )
  400                      CONTINUE
*
                           IF( JTMP.LE.( JJ + JB - 1 ) ) THEN
                              DUMMY = DCMPLX( PB_DRAND( 0 ),
     $                                       -PB_DRAND( 0 ) )
                              A( IK, JTMP ) = DCMPLX( DBLE( DUMMY ),
     $                                                ZERO )
                           END IF
*
                           DO 410 JK = JTMP + 1, JJ + JB - 1
                              A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                             -PB_DRAND( 0 ) )
  410                      CONTINUE
*
                        END IF
*
                     ELSE
*
                        DO 420 JK = JJ, JJ + JB - 1
                           A( IK, JK ) = DCMPLX( PB_DRAND( 0 ),
     $                                          -PB_DRAND( 0 ) )
  420                   CONTINUE
*
                     END IF
*
                     JJ = JJ + JB
*
                     IF( JBLK.EQ.1 ) THEN
*
*                       Jump INBLOC + ( NPCOL - 1 ) * NB columns
*
                        LCMTC = LCMTC + JMP( JMP_NQINBLOC )
                        CALL PB_JUMPIT( IMULADD( 1, JMP_NQINBLOC ), IB1,
     $                                  IB0 )
*
                     ELSE
*
*                       Jump NPCOL * NB columns
*
                        LCMTC = LCMTC + JMP( JMP_NQNB )
                        CALL PB_JUMPIT( IMULADD( 1, JMP_NQNB ), IB1,
     $                                  IB0 )
*
                     END IF
*
                     IB1( 1 ) = IB0( 1 )
                     IB1( 2 ) = IB0( 2 )
*
  430             CONTINUE
*
*                 Jump one row
*
                  CALL PB_JUMPIT( IMULADD( 1, JMP_ROW ), IB2, IB0 )
*
                  IB1( 1 ) = IB0( 1 )
                  IB1( 2 ) = IB0( 2 )
                  IB2( 1 ) = IB0( 1 )
                  IB2( 2 ) = IB0( 2 )
*
  440          CONTINUE
*
               II = II + IB
*
               IF( IBLK.EQ.1 ) THEN
*
*                 Jump IMBLOC + ( NPROW - 1 ) * MB rows
*
                  LCMTR = LCMTR - JMP( JMP_NPIMBLOC )
                  CALL PB_JUMPIT( IMULADD( 1, JMP_NPIMBLOC ), IB3, IB0 )
*
               ELSE
*
*                 Jump NPROW * MB rows
*
                  LCMTR = LCMTR - JMP( JMP_NPMB )
                  CALL PB_JUMPIT( IMULADD( 1, JMP_NPMB ), IB3, IB0 )
*
               END IF
*
               IB1( 1 ) = IB0( 1 )
               IB1( 2 ) = IB0( 2 )
               IB2( 1 ) = IB0( 1 )
               IB2( 2 ) = IB0( 2 )
               IB3( 1 ) = IB0( 1 )
               IB3( 2 ) = IB0( 2 )
*
  450       CONTINUE
*
         END IF
*
      END IF
*
      RETURN
*
*     End of PB_ZLAGEN
*
      END
      DOUBLE PRECISION   FUNCTION PB_DRAND( IDUMM )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            IDUMM
*     ..
*
*  Purpose
*  =======
*
*  PB_DRAND generates the next number in the random sequence. This func-
*  tion ensures that this number will be in the interval ( -1.0, 1.0 ).
*
*  Arguments
*  =========
*
*  IDUMM   (local input) INTEGER
*          This argument is ignored, but necessary to a FORTRAN 77 func-
*          tion.
*
*  Further Details
*  ===============
*
*  On entry, the array IRAND stored in the common block  RANCOM contains
*  the information (2 integers)  required to generate the next number in
*  the sequence X( n ). This number is computed as
*
*     X( n ) = ( 2^16 * IRAND( 2 ) + IRAND( 1 ) ) / d,
*
*  where the constant d is the  largest  32 bit  positive  integer.  The
*  array  IRAND  is  then  updated for the generation of the next number
*  X( n+1 ) in the random sequence as follows X( n+1 ) = a * X( n ) + c.
*  The constants  a  and c  should have been preliminarily stored in the
*  array  IACS  as  2 pairs of integers. The initial set up of IRAND and
*  IACS is performed by the routine PB_SETRAN.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   PB_DRAN
      EXTERNAL           PB_DRAN
*     ..
*     .. Executable Statements ..
*
      PB_DRAND = ONE - TWO * PB_DRAN( IDUMM )
*
      RETURN
*
*     End of PB_DRAND
*
      END
      DOUBLE PRECISION   FUNCTION PB_DRAN( IDUMM )
*
*  -- PBLAS test routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            IDUMM
*     ..
*
*  Purpose
*  =======
*
*  PB_DRAN generates the next number in the random sequence.
*
*  Arguments
*  =========
*
*  IDUMM   (local input) INTEGER
*          This argument is ignored, but necessary to a FORTRAN 77 func-
*          tion.
*
*  Further Details
*  ===============
*
*  On entry, the array IRAND stored in the common block  RANCOM contains
*  the information (2 integers)  required to generate the next number in
*  the sequence X( n ). This number is computed as
*
*     X( n ) = ( 2^16 * IRAND( 2 ) + IRAND( 1 ) ) / d,
*
*  where the constant d is the  largest  32 bit  positive  integer.  The
*  array  IRAND  is  then  updated for the generation of the next number
*  X( n+1 ) in the random sequence as follows X( n+1 ) = a * X( n ) + c.
*  The constants  a  and c  should have been preliminarily stored in the
*  array  IACS  as  2 pairs of integers. The initial set up of IRAND and
*  IACS is performed by the routine PB_SETRAN.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   DIVFAC, POW16
      PARAMETER          ( DIVFAC = 2.147483648D+9,
     $                   POW16 = 6.5536D+4 )
*     ..
*     .. Local Arrays ..
      INTEGER            J( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           PB_LADD, PB_LMUL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
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
      PB_DRAN = ( DBLE( IRAND( 1 ) ) + POW16 * DBLE( IRAND( 2 ) ) ) /
     $            DIVFAC
*
      CALL PB_LMUL( IRAND, IACS, J )
      CALL PB_LADD( J, IACS( 3 ), IRAND )
*
      RETURN
*
*     End of PB_DRAN
*
      END
