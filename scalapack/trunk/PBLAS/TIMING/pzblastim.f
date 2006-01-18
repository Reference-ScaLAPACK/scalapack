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
