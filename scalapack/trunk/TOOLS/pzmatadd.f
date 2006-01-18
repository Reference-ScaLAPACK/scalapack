      SUBROUTINE PZMATADD( M, N, ALPHA, A, IA, JA, DESCA, BETA, C, IC,
     $                     JC, DESCC )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, IC, JA, JC, M, N
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCC( * )
      COMPLEX*16         A( * ), C( * )
*     ..
*
*  Purpose
*  =======
*
*  PZMATADD performs a distributed matrix-matrix addition
*
*            sub( C ) := alpha * sub( A ) + beta * sub( C ),
*
*  where sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1) and sub( A ) denotes
*  A(IA:IA+M-1,JA:JA+N-1). No communications are performed in this
*  routine, the arrays are supposed to be aligned.
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
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrices sub( A ) and sub( C ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrices  sub( A ) and
*          sub( C ). N >= 0.
*
*  ALPHA   (global input) COMPLEX*16
*          The scalar ALPHA.
*
*  A       (local input) COMPLEX*16 pointer into the local memory
*          to a local array of dimension (LLD_A, LOCc(JA+N-1) ). This
*          array contains the local pieces of the distributed matrix
*          sub( A ).
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
*  BETA    (global input) COMPLEX*16
*          The scalar BETA.
*
*  C       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_C,LOCc(JC+N-1)).
*          This array contains the local pieces of the distributed
*          matrix sub( C ).  On exit, this array contains the local
*          pieces of the resulting distributed matrix.
*
*  IC      (global input) INTEGER
*          The row index in the global array C indicating the first
*          row of sub( C ).
*
*  JC      (global input) INTEGER
*          The column index in the global array C indicating the
*          first column of sub( C ).
*
*  DESCC   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix C.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                     ONE  = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IACOL, IAROW, ICCOL, ICOFF, ICROW, IIA,
     $                   IIC, IOFFA, IOFFC, IROFF, J, JJA, JJC, LDA,
     $                   LDC, MP, MYCOL, MYROW, NPCOL, NPROW, NQ
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters.
*
      CALL BLACS_GRIDINFO( DESCA(CTXT_), NPROW, NPCOL, MYROW, MYCOL )
*
*     Quick return if possible.
*
      IF( (M.EQ.0).OR.(N.EQ.0).OR.
     $    ((ALPHA.EQ.ZERO).AND.(BETA.EQ.ONE)) )
     $   RETURN
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $              IIA, JJA, IAROW, IACOL )
      CALL INFOG2L( IC, JC, DESCC, NPROW, NPCOL, MYROW, MYCOL,
     $              IIC, JJC, ICROW, ICCOL )
*
      IROFF = MOD( IA-1, DESCA(MB_) )
      ICOFF = MOD( JA-1, DESCA(NB_) )
      MP = NUMROC( M+IROFF, DESCA(MB_), MYROW, IAROW, NPROW )
      NQ = NUMROC( N+ICOFF, DESCA(NB_), MYCOL, IACOL, NPCOL )
      IF( MYROW.EQ.IAROW )
     $   MP = MP-IROFF
      IF( MYCOL.EQ.IACOL )
     $   NQ = NQ-ICOFF
      LDA = DESCA(LLD_)
      LDC = DESCC(LLD_)
*
      IF( NQ.EQ.1 ) THEN
         IF( BETA.EQ.ZERO ) THEN
            IF( ALPHA.EQ.ZERO ) THEN
               IOFFC = IIC + (JJC-1)*LDC
               DO 10 I = IOFFC, IOFFC+MP-1
                  C( I ) = ZERO
   10          CONTINUE
            ELSE
               IOFFA = IIA + (JJA-1)*LDA
               IOFFC = IIC + (JJC-1)*LDC
               DO 20 I = IOFFC, IOFFC+MP-1
                  C( I ) = ALPHA * A( IOFFA )
                  IOFFA = IOFFA + 1
   20          CONTINUE
            END IF
         ELSE
            IF( ALPHA.EQ.ONE ) THEN
               IF( BETA.EQ.ONE ) THEN
                  IOFFA = IIA + (JJA-1)*LDA
                  IOFFC = IIC + (JJC-1)*LDC
                  DO 30 I = IOFFC, IOFFC+MP-1
                     C( I ) = C( I ) + A( IOFFA )
                     IOFFA = IOFFA + 1
   30             CONTINUE
               ELSE
                  IOFFA = IIA + (JJA-1)*LDA
                  IOFFC = IIC + (JJC-1)*LDC
                  DO 40 I = IOFFC, IOFFC+MP-1
                     C( I ) = BETA * C( I ) + A( IOFFA )
                     IOFFA = IOFFA + 1
   40             CONTINUE
               END IF
            ELSE IF( BETA.EQ.ONE ) THEN
               IOFFA = IIA + (JJA-1)*LDA
               IOFFC = IIC + (JJC-1)*LDC
               DO 50 I = IOFFC, IOFFC+MP-1
                  C( I ) = C( I ) + ALPHA * A( IOFFA )
                  IOFFA = IOFFA + 1
   50          CONTINUE
            ELSE
               IOFFA = IIA + (JJA-1)*LDA
               IOFFC = IIC + (JJC-1)*LDC
               DO 60 I = IOFFC, IOFFC+MP-1
                  C( I ) = BETA * C( I ) + ALPHA * A( IOFFA )
                  IOFFA = IOFFA + 1
   60          CONTINUE
            END IF
         END IF
      ELSE
         IF( BETA.EQ.ZERO ) THEN
            IF( ALPHA.EQ.ZERO ) THEN
               IOFFC = IIC+(JJC-1)*LDC
               DO 80 J = 1, NQ
                  DO 70 I = IOFFC, IOFFC+MP-1
                     C( I ) = ZERO
   70             CONTINUE
                  IOFFC = IOFFC + LDC
   80          CONTINUE
            ELSE
               IOFFA = IIA+(JJA-1)*LDA
               IOFFC = IIC+(JJC-1)*LDC
               DO 100 J = 1, NQ
                  DO 90 I = IOFFC, IOFFC+MP-1
                     C( I ) = ALPHA * A( IOFFA )
                     IOFFA = IOFFA + 1
   90             CONTINUE
                  IOFFA = IOFFA + LDA - MP
                  IOFFC = IOFFC + LDC
  100          CONTINUE
            END IF
         ELSE
            IF( ALPHA.EQ.ONE ) THEN
               IF( BETA.EQ.ONE ) THEN
                  IOFFA = IIA+(JJA-1)*LDA
                  IOFFC = IIC+(JJC-1)*LDC
                  DO 120 J = 1, NQ
                     DO 110 I = IOFFC, IOFFC+MP-1
                        C( I ) = C( I ) + A( IOFFA )
                        IOFFA = IOFFA + 1
  110                CONTINUE
                     IOFFA = IOFFA + LDA - MP
                     IOFFC = IOFFC + LDC
  120             CONTINUE
               ELSE
                  IOFFA = IIA+(JJA-1)*LDA
                  IOFFC = IIC+(JJC-1)*LDC
                  DO 140 J = 1, NQ
                     DO 130 I = IOFFC, IOFFC+MP-1
                        C( I ) = BETA * C( I ) + A( IOFFA )
                        IOFFA = IOFFA + 1
  130                CONTINUE
                     IOFFA = IOFFA + LDA - MP
                     IOFFC = IOFFC + LDC
  140             CONTINUE
               END IF
            ELSE IF( BETA.EQ.ONE ) THEN
               IOFFA = IIA+(JJA-1)*LDA
               IOFFC = IIC+(JJC-1)*LDC
               DO 160 J = 1, NQ
                  DO 150 I = IOFFC, IOFFC+MP-1
                     C( I ) = C( I ) + ALPHA * A( IOFFA )
                     IOFFA = IOFFA + 1
  150             CONTINUE
                  IOFFA = IOFFA + LDA - MP
                  IOFFC = IOFFC + LDC
  160          CONTINUE
            ELSE
               IOFFA = IIA+(JJA-1)*LDA
               IOFFC = IIC+(JJC-1)*LDC
               DO 180 J = 1, NQ
                  DO 170 I = IOFFC, IOFFC+MP-1
                     C( I ) = BETA * C( I ) + ALPHA * A( IOFFA )
                     IOFFA = IOFFA + 1
  170             CONTINUE
                  IOFFA = IOFFA + LDA - MP
                  IOFFC = IOFFC + LDC
  180          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of PZMATADD
*
      END
