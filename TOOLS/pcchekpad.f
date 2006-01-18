      SUBROUTINE PCCHEKPAD( ICTXT, MESS, M, N, A, LDA, IPRE, IPOST,
     $                     CHKVAL )
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, IPOST, IPRE, LDA, M, N
      COMPLEX            CHKVAL
*     ..
*     .. Array Arguments ..
      CHARACTER          MESS*(*)
      COMPLEX            A( * )
*     ..
*
*  Purpose
*  =======
*
*  PCCHEKPAD checks that the padding around a local array has not
*  been overwritten since the call to PCFILLPAD.  3 types of errors
*  are reported:
*
*  1) Overwrite in pre-guardzone. This indicates a memory overwrite has
*  occurred in the first IPRE elements which form a buffer before the
*  beginning of A.  Therefore, the error message:
*     'Overwrite in  pre-guardzone: loc(  5) =         18.00000'
*  tells you that the 5th element of the IPRE long buffer has been
*  overwritten with the value 18, where it should still have the value
*  of CHKVAL.
*
*  2) Overwrite in post-guardzone. This indicates a memory overwrite has
*  occurred in the last IPOST elements which form a buffer after the end
*  of A. Error reports are refered from the end of A.  Therefore,
*     'Overwrite in post-guardzone: loc( 19) =         24.00000'
*  tells you that the 19th element after the end of A was overwritten
*  with the value 24, where it should still have the value of CHKVAL.
*
*  3) Overwrite in lda-m gap.  Tells you elements between M and LDA were
*  overwritten.  So,
*     'Overwrite in lda-m gap: A( 12,  3) =         22.00000'
*  tells you that the element at the 12th row and 3rd column of A was
*  overwritten with the value of 22, where it should still have the
*  value of CHKVAL.
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation. The context itself is global.
*
*  MESS    (local input) CHARACTER*(*)
*          String containing a user-defined message.
*
*  M       (local input) INTEGER
*          The number of rows in the local array A.
*
*  N       (input) INTEGER
*          The number of columns in the local array A.
*
*  A       (local input) COMPLEX array of dimension (LDA,N).
*          A location IPRE elements in front of the array to be checked.
*
*  LDA     (local input) INTEGER
*          The leading Dimension of the local array to be checked.
*
*  IPRE    (local input) INTEGER
*          The size of the guard zone before the start of padded array.
*
*  IPOST   (local input) INTEGER
*          The size of guard zone after the padded array.
*
*  CHKVAL  (local input) COMPLEX
*          The value the local array was padded with.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IAM, IDUMM, INFO, J, K, MYCOL, MYROW,
     $                   NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, IGAMX2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, AIMAG
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      IAM = MYROW*NPCOL + MYCOL
      INFO = -1
*
*     Check buffer in front of A
*
      IF( IPRE.GT.0 ) THEN
         DO 10 I = 1, IPRE
            IF( A( I ).NE.CHKVAL ) THEN
               WRITE( *, FMT = 9998 ) MYROW, MYCOL, MESS, ' pre', I,
     $                                REAL( A( I ) ), AIMAG( A( I ) )
               INFO = IAM
            END IF
   10    CONTINUE
      ELSE
         WRITE( *, FMT = * ) 'WARNING no pre-guardzone in PCCHEKPAD'
      END IF
*
*     Check buffer after A
*
      IF( IPOST.GT.0 ) THEN
         J = IPRE+LDA*N+1
         DO 20 I = J, J+IPOST-1
            IF( A( I ).NE.CHKVAL ) THEN
               WRITE( *, FMT = 9998 ) MYROW, MYCOL, MESS, 'post',
     $                                I-J+1, REAL( A( I ) ),
     $                                AIMAG( A( I ) )
               INFO = IAM
            END IF
   20    CONTINUE
      ELSE
         WRITE( *, FMT = * )
     $          'WARNING no post-guardzone buffer in PCCHEKPAD'
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
     $               I-IPRE-LDA*(J-1), J, REAL( A( I ) ),
     $               AIMAG( A( I ) )
                  INFO = IAM
               END IF
   30       CONTINUE
            K = K + LDA
   40    CONTINUE
      END IF
*
      CALL IGAMX2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, IDUMM, IDUMM, -1,
     $              0, 0 )
      IF( IAM.EQ.0 .AND. INFO.GE.0 ) THEN
         WRITE( *, FMT = 9999 ) INFO / NPCOL, MOD( INFO, NPCOL ), MESS
      END IF
*
 9999 FORMAT( '{', I5, ',', I5, '}:  Memory overwrite in ', A )
 9998 FORMAT( '{', I5, ',', I5, '}:  ', A, ' memory overwrite in ',
     $        A4, '-guardzone: loc(', I3, ') = ', G11.4, '+ i*',
     $        G11.4 )
 9997 FORMAT( '{', I5, ',', I5, '}: ', A, ' memory overwrite in ',
     $        'lda-m gap: loc(', I3, ',', I3, ') = ', G11.4,
     $        '+ i*', G11.4 )
*
      RETURN
*
*     End of PCCHEKPAD
*
      END
