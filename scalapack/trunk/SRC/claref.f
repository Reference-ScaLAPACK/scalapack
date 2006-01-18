      SUBROUTINE CLAREF( TYPE, A, LDA, WANTZ, Z, LDZ, BLOCK, IROW1,
     $                   ICOL1, ISTART, ISTOP, ITMP1, ITMP2, LILOZ,
     $                   LIHIZ, VECS, V2, V3, T1, T2, T3 )
*
*  -- ScaLAPACK routine (version 1.7) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     May 28, 1999
*
*     .. Scalar Arguments ..
      LOGICAL            BLOCK, WANTZ
      CHARACTER          TYPE
      INTEGER            ICOL1, IROW1, ISTART, ISTOP, ITMP1, ITMP2, LDA,
     $                   LDZ, LIHIZ, LILOZ
      COMPLEX            T1, T2, T3, V2, V3
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), VECS( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CLAREF applies one or several Householder reflectors of size 3
*     to one or two matrices (if column is specified) on either their
*     rows or columns.
*
*  Arguments
*  =========
*
*  TYPE    (global input) CHARACTER*1
*          If 'R': Apply reflectors to the rows of the matrix
*              (apply from left)
*          Otherwise: Apply reflectors to the columns of the matrix
*          Unchanged on exit.
*
*  A       (global input/output) COMPLEX array, (LDA,*)
*          On entry, the matrix to receive the reflections.
*          The updated matrix on exit.
*
*  LDA     (local input) INTEGER
*          On entry, the leading dimension of A.  Unchanged on exit.
*
*  WANTZ   (global input) LOGICAL
*          If .TRUE., then apply any column reflections to Z as well.
*          If .FALSE., then do no additional work on Z.
*
*  Z       (global input/output) COMPLEX array, (LDZ,*)
*          On entry, the second matrix to receive column reflections.
*          This is changed only if WANTZ is set.
*
*  LDZ     (local input) INTEGER
*          On entry, the leading dimension of Z.  Unchanged on exit.
*
*  BLOCK   (global input) LOGICAL
*          If .TRUE., then apply several reflectors at once and read
*             their data from the VECS array.
*          If .FALSE., apply the single reflector given by V2, V3,
*             T1, T2, and T3.
*
*  IROW1   (local input/output) INTEGER
*          On entry, the local row element of A.
*          Undefined on output.
*
*
*  ICOL1   (local input/output) INTEGER
*          On entry, the local column element of A.
*          Undefined on output.
*
*  ISTART  (global input) INTEGER
*          Specifies the "number" of the first reflector.  This is
*              used as an index into VECS if BLOCK is set.
*              ISTART is ignored if BLOCK is .FALSE..
*
*  ISTOP   (global input) INTEGER
*          Specifies the "number" of the last reflector.  This is
*              used as an index into VECS if BLOCK is set.
*              ISTOP is ignored if BLOCK is .FALSE..
*
*  ITMP1   (local input) INTEGER
*          Starting range into A.  For rows, this is the local
*              first column.  For columns, this is the local first row.
*
*  ITMP2   (local input) INTEGER
*          Ending range into A.  For rows, this is the local last
*              column.  For columns, this is the local last row.
*
*  LILOZ
*  LIHIZ   (local input) INTEGER
*          These serve the same purpose as ITMP1,ITMP2 but for Z
*              when WANTZ is set.
*
*  VECS    (global input) COMPLEX array of size 3*N (matrix size)
*          This holds the size 3 reflectors one after another and this
*              is only accessed when BLOCK is .TRUE.
*
*  V2
*  V3
*  T1
*  T2
*  T3      (global input/output) COMPLEX
*          This holds information on a single size 3 Householder
*              reflector and is read when BLOCK is .FALSE., and
*              overwritten when BLOCK is .TRUE.
*
*  Further Details
*  ===============
*
*  Implemented by:  M. Fahey, May 28, 1999
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            J, K
      COMPLEX            A1, A11, A2, A22, A3, A4, A5, B1, B2, B3, B4,
     $                   B5, H11, H22, SUM, SUM1, SUM2, SUM3, T12, T13,
     $                   T22, T23, T32, T33, TMP1, TMP2, TMP3, V22, V23,
     $                   V32, V33
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MOD
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( TYPE, 'R' ) ) THEN
         IF( BLOCK ) THEN
            DO 30 K = ISTART, ISTOP - MOD( ISTOP-ISTART+1, 3 ), 3
               V2 = VECS( ( K-1 )*3+1 )
               V3 = VECS( ( K-1 )*3+2 )
               T1 = VECS( ( K-1 )*3+3 )
               V22 = VECS( ( K-1 )*3+4 )
               V32 = VECS( ( K-1 )*3+5 )
               T12 = VECS( ( K-1 )*3+6 )
               V23 = VECS( ( K-1 )*3+7 )
               V33 = VECS( ( K-1 )*3+8 )
               T13 = VECS( ( K-1 )*3+9 )
               T2 = T1*V2
               T3 = T1*V3
               T22 = T12*V22
               T32 = T12*V32
               T23 = T13*V23
               T33 = T13*V33
               DO 10 J = ITMP1, ITMP2 - MOD( ITMP2-ITMP1+1, 2 ), 2
                  A1 = A( IROW1, J )
                  A2 = A( IROW1+1, J )
                  A3 = A( IROW1+2, J )
                  A4 = A( IROW1+3, J )
                  A5 = A( IROW1+4, J )
                  B1 = A( IROW1, J+1 )
                  B2 = A( IROW1+1, J+1 )
                  B3 = A( IROW1+2, J+1 )
                  B4 = A( IROW1+3, J+1 )
                  B5 = A( IROW1+4, J+1 )
                  SUM1 = CONJG( T1 )*A1 + CONJG( T2 )*A2 +
     $                   CONJG( T3 )*A3
                  A( IROW1, J ) = A1 - SUM1
                  H11 = A2 - SUM1*V2
                  H22 = A3 - SUM1*V3
                  TMP1 = CONJG( T1 )*B1 + CONJG( T2 )*B2 +
     $                   CONJG( T3 )*B3
                  A( IROW1, J+1 ) = B1 - TMP1
                  A11 = B2 - TMP1*V2
                  A22 = B3 - TMP1*V3
                  SUM2 = CONJG( T12 )*H11 + CONJG( T22 )*H22 +
     $                   CONJG( T32 )*A4
                  A( IROW1+1, J ) = H11 - SUM2
                  H11 = H22 - SUM2*V22
                  H22 = A4 - SUM2*V32
                  TMP2 = CONJG( T12 )*A11 + CONJG( T22 )*A22 +
     $                   CONJG( T32 )*B4
                  A( IROW1+1, J+1 ) = A11 - TMP2
                  A11 = A22 - TMP2*V22
                  A22 = B4 - TMP2*V32
                  SUM3 = CONJG( T13 )*H11 + CONJG( T23 )*H22 +
     $                   CONJG( T33 )*A5
                  A( IROW1+2, J ) = H11 - SUM3
                  A( IROW1+3, J ) = H22 - SUM3*V23
                  A( IROW1+4, J ) = A5 - SUM3*V33
                  TMP3 = CONJG( T13 )*A11 + CONJG( T23 )*A22 +
     $                   CONJG( T33 )*B5
                  A( IROW1+2, J+1 ) = A11 - TMP3
                  A( IROW1+3, J+1 ) = A22 - TMP3*V23
                  A( IROW1+4, J+1 ) = B5 - TMP3*V33
   10          CONTINUE
               DO 20 J = ITMP2 - MOD( ITMP2-ITMP1+1, 2 ) + 1, ITMP2
                  SUM = CONJG( T1 )*A( IROW1, J ) +
     $                  CONJG( T2 )*A( IROW1+1, J ) +
     $                  CONJG( T3 )*A( IROW1+2, J )
                  A( IROW1, J ) = A( IROW1, J ) - SUM
                  H11 = A( IROW1+1, J ) - SUM*V2
                  H22 = A( IROW1+2, J ) - SUM*V3
                  SUM = CONJG( T12 )*H11 + CONJG( T22 )*H22 +
     $                  CONJG( T32 )*A( IROW1+3, J )
                  A( IROW1+1, J ) = H11 - SUM
                  H11 = H22 - SUM*V22
                  H22 = A( IROW1+3, J ) - SUM*V32
                  SUM = CONJG( T13 )*H11 + CONJG( T23 )*H22 +
     $                  CONJG( T33 )*A( IROW1+4, J )
                  A( IROW1+2, J ) = H11 - SUM
                  A( IROW1+3, J ) = H22 - SUM*V23
                  A( IROW1+4, J ) = A( IROW1+4, J ) - SUM*V33
   20          CONTINUE
               IROW1 = IROW1 + 3
   30       CONTINUE
            DO 50 K = ISTOP - MOD( ISTOP-ISTART+1, 3 ) + 1, ISTOP
               V2 = VECS( ( K-1 )*3+1 )
               V3 = VECS( ( K-1 )*3+2 )
               T1 = VECS( ( K-1 )*3+3 )
               T2 = T1*V2
               T3 = T1*V3
               DO 40 J = ITMP1, ITMP2
                  SUM = CONJG( T1 )*A( IROW1, J ) +
     $                  CONJG( T2 )*A( IROW1+1, J ) +
     $                  CONJG( T3 )*A( IROW1+2, J )
                  A( IROW1, J ) = A( IROW1, J ) - SUM
                  A( IROW1+1, J ) = A( IROW1+1, J ) - SUM*V2
                  A( IROW1+2, J ) = A( IROW1+2, J ) - SUM*V3
   40          CONTINUE
               IROW1 = IROW1 + 1
   50       CONTINUE
         ELSE
            DO 60 J = ITMP1, ITMP2
               SUM = CONJG( T1 )*A( IROW1, J ) +
     $               CONJG( T2 )*A( IROW1+1, J ) +
     $               CONJG( T3 )*A( IROW1+2, J )
               A( IROW1, J ) = A( IROW1, J ) - SUM
               A( IROW1+1, J ) = A( IROW1+1, J ) - SUM*V2
               A( IROW1+2, J ) = A( IROW1+2, J ) - SUM*V3
   60       CONTINUE
         END IF
      ELSE
*
*        Do column transforms
*
         IF( BLOCK ) THEN
            DO 90 K = ISTART, ISTOP - MOD( ISTOP-ISTART+1, 3 ), 3
               V2 = VECS( ( K-1 )*3+1 )
               V3 = VECS( ( K-1 )*3+2 )
               T1 = VECS( ( K-1 )*3+3 )
               V22 = VECS( ( K-1 )*3+4 )
               V32 = VECS( ( K-1 )*3+5 )
               T12 = VECS( ( K-1 )*3+6 )
               V23 = VECS( ( K-1 )*3+7 )
               V33 = VECS( ( K-1 )*3+8 )
               T13 = VECS( ( K-1 )*3+9 )
               T2 = T1*V2
               T3 = T1*V3
               T22 = T12*V22
               T32 = T12*V32
               T23 = T13*V23
               T33 = T13*V33
               DO 70 J = ITMP1, ITMP2
                  SUM = T1*A( J, ICOL1 ) + T2*A( J, ICOL1+1 ) +
     $                  T3*A( J, ICOL1+2 )
                  A( J, ICOL1 ) = A( J, ICOL1 ) - SUM
                  H11 = A( J, ICOL1+1 ) - SUM*CONJG( V2 )
                  H22 = A( J, ICOL1+2 ) - SUM*CONJG( V3 )
                  SUM = T12*H11 + T22*H22 + T32*A( J, ICOL1+3 )
                  A( J, ICOL1+1 ) = H11 - SUM
                  H11 = H22 - SUM*CONJG( V22 )
                  H22 = A( J, ICOL1+3 ) - SUM*CONJG( V32 )
                  SUM = T13*H11 + T23*H22 + T33*A( J, ICOL1+4 )
                  A( J, ICOL1+2 ) = H11 - SUM
                  A( J, ICOL1+3 ) = H22 - SUM*CONJG( V23 )
                  A( J, ICOL1+4 ) = A( J, ICOL1+4 ) - SUM*CONJG( V33 )
   70          CONTINUE
               IF( WANTZ ) THEN
                  DO 80 J = LILOZ, LIHIZ
                     SUM = T1*Z( J, ICOL1 ) + T2*Z( J, ICOL1+1 ) +
     $                     T3*Z( J, ICOL1+2 )
                     Z( J, ICOL1 ) = Z( J, ICOL1 ) - SUM
                     H11 = Z( J, ICOL1+1 ) - SUM*CONJG( V2 )
                     H22 = Z( J, ICOL1+2 ) - SUM*CONJG( V3 )
                     SUM = T12*H11 + T22*H22 + T32*Z( J, ICOL1+3 )
                     Z( J, ICOL1+1 ) = H11 - SUM
                     H11 = H22 - SUM*CONJG( V22 )
                     H22 = Z( J, ICOL1+3 ) - SUM*CONJG( V32 )
                     SUM = T13*H11 + T23*H22 + T33*Z( J, ICOL1+4 )
                     Z( J, ICOL1+2 ) = H11 - SUM
                     Z( J, ICOL1+3 ) = H22 - SUM*CONJG( V23 )
                     Z( J, ICOL1+4 ) = Z( J, ICOL1+4 ) -
     $                                 SUM*CONJG( V33 )
   80             CONTINUE
               END IF
               ICOL1 = ICOL1 + 3
   90       CONTINUE
            DO 120 K = ISTOP - MOD( ISTOP-ISTART+1, 3 ) + 1, ISTOP
               V2 = VECS( ( K-1 )*3+1 )
               V3 = VECS( ( K-1 )*3+2 )
               T1 = VECS( ( K-1 )*3+3 )
               T2 = T1*V2
               T3 = T1*V3
               DO 100 J = ITMP1, ITMP2
                  SUM = T1*A( J, ICOL1 ) + T2*A( J, ICOL1+1 ) +
     $                  T3*A( J, ICOL1+2 )
                  A( J, ICOL1 ) = A( J, ICOL1 ) - SUM
                  A( J, ICOL1+1 ) = A( J, ICOL1+1 ) - SUM*CONJG( V2 )
                  A( J, ICOL1+2 ) = A( J, ICOL1+2 ) - SUM*CONJG( V3 )
  100          CONTINUE
               IF( WANTZ ) THEN
                  DO 110 J = LILOZ, LIHIZ
                     SUM = T1*Z( J, ICOL1 ) + T2*Z( J, ICOL1+1 ) +
     $                     T3*Z( J, ICOL1+2 )
                     Z( J, ICOL1 ) = Z( J, ICOL1 ) - SUM
                     Z( J, ICOL1+1 ) = Z( J, ICOL1+1 ) -
     $                                 SUM*CONJG( V2 )
                     Z( J, ICOL1+2 ) = Z( J, ICOL1+2 ) -
     $                                 SUM*CONJG( V3 )
  110             CONTINUE
               END IF
               ICOL1 = ICOL1 + 1
  120       CONTINUE
         ELSE
            DO 130 J = ITMP1, ITMP2
               SUM = T1*A( J, ICOL1 ) + T2*A( J, ICOL1+1 ) +
     $               T3*A( J, ICOL1+2 )
               A( J, ICOL1 ) = A( J, ICOL1 ) - SUM
               A( J, ICOL1+1 ) = A( J, ICOL1+1 ) - SUM*CONJG( V2 )
               A( J, ICOL1+2 ) = A( J, ICOL1+2 ) - SUM*CONJG( V3 )
  130       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of CLAREF
*
      END
