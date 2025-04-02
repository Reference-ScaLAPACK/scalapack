      SUBROUTINE PCLATRAN( N, NB, A, IA, JA, DESCA, WORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     October 15, 1999
*
*     .. Scalar Arguments ..
      INTEGER            IA, JA, N, NB
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX            A( * ), WORK( * )
*     ..
*
*  Purpose
*
*  =======
*
*  PCLATRAN transpose a lower triangular matrix on to the upper
*  triangular portion of the same matrix.
*
*  This is an auxiliary routine called by PCHETRD.
*
*  Notes
*  =====
*
*  IA must equal 1
*  JA must equal 1
*  DESCA( MB_ ) must equal 1
*  DESCA( NB_ ) must equal 1
*  DESCA( RSRC_ ) must equal 1
*  DESCA( CSRC_ ) must equal 1
*
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The size of the matrix to be transposed.
*
*  NB      (global input) INTEGER
*          The number of rows and columns to be transposed with each
*          message sent.  NB has no impact on the result, it is striclty
*          a performance tuning parameter.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the
*          Hermitian distributed matrix sub( A ).  On entry, the
*          leading N-by-N upper triangular part of sub( A ) contains
*          the upper triangular part of the matrix. On exit, the
*          leading N-by-N lower triangular part of sub( A ) contains the
*          lower triangular part of the matrix, and its strictly upper
*          triangular part is undefined (and may have been modified).
*
*  IA      (global input) INTEGER
*          A's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*          Must be equal to 1.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*          Must be equal to 1.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*          DESCA( MB_ ) must equal 1
*          DESCA( NB_ ) must equal 1
*          DESCA( ICTXT_ ) must point to a square process grid
*            i.e. one where NPROW is equal to NPCOL
*
*  WORK    (local workspace) COMPLEX*16 array, dimension ( LWORK )
*
*          Where:
*          LWORK >= NB * NUMROC( N, 1, 0, 0, NPROW )
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
      INTEGER            I, ICTXT, IRECV, ISEND, J, JJ, JRECV, JSEND,
     $                   LDA, MAXIRECV, MAXISEND, MAXJRECV, MAXJSEND,
     $                   MINIRECV, MINISEND, MINJRECV, MINJSEND, MYCOL,
     $                   MYROW, NP, NPCOL, NPROW, NQ, RECVNB, SENDNB,
     $                   STARTCOL, STARTROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CTRRV2D, CTRSD2D
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, MIN
*     ..
*     .. Executable Statements ..
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
*     Further details
*
*     Because the processor grid is square each process needs only send
*     data to its transpose process.  (Likewsie it need only receive
*     data from its transpose process.)  Because the data decomposition
*     is cyclic, the local portion of the array is triangular.
*
*     This routine requires that the data be buffered (i.e. copied)
*     on the sending process (because of the triangular shape) and
*     unbuffered on the receiving process.  Hence, two local memory to
*     memory copies are performed within the communications routines
*     followed by a memory to memory copy outside of the communications
*     routines.  It would be nice to avoid having back to back memory
*     to memory copies (as we do presently on the receiving processor).
*     This could be done by packaging the data ourselves in the sender
*     and then unpacking it directly into the matrix.  However, this
*     code seems cleaner and so since this routine is not a significant
*     performance bottleneck we have left it this way.
*
*
*
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      ICTXT = DESCA( CTXT_ )
      LDA = DESCA( LLD_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*
      NP = NUMROC( N, 1, MYROW, 0, NPROW )
      NQ = NUMROC( N, 1, MYCOL, 0, NPCOL )
*
*
      IF( MYROW.EQ.MYCOL ) THEN
*
         DO 20 J = 1, NP
            DO 10 I = J + 1, NQ
               A( J+( I-1 )*LDA ) = CONJG( A( I+( J-1 )*LDA ) )
   10       CONTINUE
   20    CONTINUE
*
      ELSE
         IF( MYROW.GT.MYCOL ) THEN
            STARTROW = 1
            STARTCOL = 2
         ELSE
            IF( MYROW.EQ.MYCOL ) THEN
               STARTROW = 2
               STARTCOL = 2
            ELSE
               STARTROW = 2
               STARTCOL = 1
            END IF
         END IF
*
         DO 50 JJ = 1, MAX( NP, NQ ), NB
            MINJSEND = STARTCOL + JJ - 1
            MINJRECV = STARTROW + JJ - 1
            MAXJSEND = MIN( MINJSEND+NB-1, NQ )
            MAXJRECV = MIN( MINJRECV+NB-1, NP )
*
            SENDNB = MAXJSEND - MINJSEND + 1
            RECVNB = MAXJRECV - MINJRECV + 1
*
            MINISEND = 1
            MINIRECV = 1
            MAXISEND = MIN( NP, JJ+SENDNB-1 )
            MAXIRECV = MIN( NQ, JJ+RECVNB-1 )
*
            ISEND = MAXISEND - MINISEND + 1
            IRECV = MAXIRECV - MINIRECV + 1
            JSEND = MAXJSEND - MINJSEND + 1
            JRECV = MAXJRECV - MINJRECV + 1
*
*
*
            DO 40 J = MINJRECV, MAXJRECV
               DO 30 I = MINIRECV, MAXIRECV + J - MAXJRECV
                  WORK( I+( J-MINJRECV )*IRECV )
     $               = CONJG( A( J+( I-1 )*LDA ) )
   30          CONTINUE
   40       CONTINUE
*
            IF( IRECV.GT.0 .AND. JRECV.GT.0 )
     $         CALL CTRSD2D( ICTXT, 'U', 'N', IRECV, JRECV, WORK, IRECV,
     $                       MYCOL, MYROW )
*
            IF( ISEND.GT.0 .AND. JSEND.GT.0 )
     $         CALL CTRRV2D( ICTXT, 'U', 'N', ISEND, JSEND,
     $                       A( MINISEND+( MINJSEND-1 )*LDA ), LDA,
     $                       MYCOL, MYROW )
*
*
   50    CONTINUE
*
      END IF
*
      RETURN
*
*     End of PCLATRD
*
      END
