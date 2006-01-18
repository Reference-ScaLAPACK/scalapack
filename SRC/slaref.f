      SUBROUTINE SLAREF( TYPE, A, LDA, WANTZ, Z, LDZ, BLOCK, IROW1,
     $                    ICOL1, ISTART, ISTOP, ITMP1, ITMP2, LILOZ,
     $                    LIHIZ, VECS, V2, V3, T1, T2, T3 )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      LOGICAL            BLOCK, WANTZ
      CHARACTER          TYPE
      INTEGER            ICOL1, IROW1, ISTART, ISTOP, ITMP1, ITMP2, LDA,
     $                   LDZ, LIHIZ, LILOZ
      REAL               T1, T2, T3, V2, V3
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), VECS( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  SLAREF applies one or several Householder reflectors of size 3
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
*  A       (global input/output) REAL array, (LDA,*)
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
*  Z       (global input/output) REAL array, (LDZ,*)
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
*  VECS    (global input) REAL array of size 3*N (matrix
*                                                 size)
*          This holds the size 3 reflectors one after another and this
*              is only accessed when BLOCK is .TRUE.
*
*  V2
*  V3
*  T1
*  T2
*  T3      (global input/output) REAL
*          This holds information on a single size 3 Householder
*              reflector and is read when BLOCK is .FALSE., and
*              overwritten when BLOCK is .TRUE.
*
*  Implemented by:  G. Henry, November 17, 1996
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            J, K
      REAL               H11, H22, SUM, T12, T13, T22, T23, T32, T33,
     $                   V22, V23, V32, V33
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( TYPE, 'R' ) ) THEN
         IF( BLOCK ) THEN
            DO 20 K = ISTART, ISTOP - MOD( ISTOP-ISTART+1, 3 ), 3
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
               DO 10 J = ITMP1, ITMP2
                  SUM = A( IROW1, J ) + V2*A( IROW1+1, J ) +
     $                  V3*A( IROW1+2, J )
                  A( IROW1, J ) = A( IROW1, J ) - SUM*T1
                  H11 = A( IROW1+1, J ) - SUM*T2
                  H22 = A( IROW1+2, J ) - SUM*T3
                  SUM = H11 + V22*H22 + V32*A( IROW1+3, J )
                  A( IROW1+1, J ) = H11 - SUM*T12
                  H11 = H22 - SUM*T22
                  H22 = A( IROW1+3, J ) - SUM*T32
                  SUM = H11 + V23*H22 + V33*A( IROW1+4, J )
                  A( IROW1+2, J ) = H11 - SUM*T13
                  A( IROW1+3, J ) = H22 - SUM*T23
                  A( IROW1+4, J ) = A( IROW1+4, J ) - SUM*T33
   10          CONTINUE
               IROW1 = IROW1 + 3
   20       CONTINUE
            DO 40 K = ISTOP - MOD( ISTOP-ISTART+1, 3 ) + 1, ISTOP
               V2 = VECS( ( K-1 )*3+1 )
               V3 = VECS( ( K-1 )*3+2 )
               T1 = VECS( ( K-1 )*3+3 )
               T2 = T1*V2
               T3 = T1*V3
               DO 30 J = ITMP1, ITMP2
                  SUM = A( IROW1, J ) + V2*A( IROW1+1, J ) +
     $                  V3*A( IROW1+2, J )
                  A( IROW1, J ) = A( IROW1, J ) - SUM*T1
                  A( IROW1+1, J ) = A( IROW1+1, J ) - SUM*T2
                  A( IROW1+2, J ) = A( IROW1+2, J ) - SUM*T3
   30          CONTINUE
               IROW1 = IROW1 + 1
   40       CONTINUE
         ELSE
            DO 50 J = ITMP1, ITMP2
               SUM = A( IROW1, J ) + V2*A( IROW1+1, J ) +
     $               V3*A( IROW1+2, J )
               A( IROW1, J ) = A( IROW1, J ) - SUM*T1
               A( IROW1+1, J ) = A( IROW1+1, J ) - SUM*T2
               A( IROW1+2, J ) = A( IROW1+2, J ) - SUM*T3
   50       CONTINUE
         END IF
      ELSE
*
*        Do column transforms
*
         IF( BLOCK ) THEN
            DO 80 K = ISTART, ISTOP - MOD( ISTOP-ISTART+1, 3 ), 3
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
               DO 60 J = ITMP1, ITMP2
                  SUM = A( J, ICOL1 ) + V2*A( J, ICOL1+1 ) +
     $                  V3*A( J, ICOL1+2 )
                  A( J, ICOL1 ) = A( J, ICOL1 ) - SUM*T1
                  H11 = A( J, ICOL1+1 ) - SUM*T2
                  H22 = A( J, ICOL1+2 ) - SUM*T3
                  SUM = H11 + V22*H22 + V32*A( J, ICOL1+3 )
                  A( J, ICOL1+1 ) = H11 - SUM*T12
                  H11 = H22 - SUM*T22
                  H22 = A( J, ICOL1+3 ) - SUM*T32
                  SUM = H11 + V23*H22 + V33*A( J, ICOL1+4 )
                  A( J, ICOL1+2 ) = H11 - SUM*T13
                  A( J, ICOL1+3 ) = H22 - SUM*T23
                  A( J, ICOL1+4 ) = A( J, ICOL1+4 ) - SUM*T33
   60          CONTINUE
               IF( WANTZ ) THEN
                  DO 70 J = LILOZ, LIHIZ
                     SUM = Z( J, ICOL1 ) + V2*Z( J, ICOL1+1 ) +
     $                     V3*Z( J, ICOL1+2 )
                     Z( J, ICOL1 ) = Z( J, ICOL1 ) - SUM*T1
                     H11 = Z( J, ICOL1+1 ) - SUM*T2
                     H22 = Z( J, ICOL1+2 ) - SUM*T3
                     SUM = H11 + V22*H22 + V32*Z( J, ICOL1+3 )
                     Z( J, ICOL1+1 ) = H11 - SUM*T12
                     H11 = H22 - SUM*T22
                     H22 = Z( J, ICOL1+3 ) - SUM*T32
                     SUM = H11 + V23*H22 + V33*Z( J, ICOL1+4 )
                     Z( J, ICOL1+2 ) = H11 - SUM*T13
                     Z( J, ICOL1+3 ) = H22 - SUM*T23
                     Z( J, ICOL1+4 ) = Z( J, ICOL1+4 ) - SUM*T33
   70             CONTINUE
               END IF
               ICOL1 = ICOL1 + 3
   80       CONTINUE
            DO 110 K = ISTOP - MOD( ISTOP-ISTART+1, 3 ) + 1, ISTOP
               V2 = VECS( ( K-1 )*3+1 )
               V3 = VECS( ( K-1 )*3+2 )
               T1 = VECS( ( K-1 )*3+3 )
               T2 = T1*V2
               T3 = T1*V3
               DO 90 J = ITMP1, ITMP2
                  SUM = A( J, ICOL1 ) + V2*A( J, ICOL1+1 ) +
     $                  V3*A( J, ICOL1+2 )
                  A( J, ICOL1 ) = A( J, ICOL1 ) - SUM*T1
                  A( J, ICOL1+1 ) = A( J, ICOL1+1 ) - SUM*T2
                  A( J, ICOL1+2 ) = A( J, ICOL1+2 ) - SUM*T3
   90          CONTINUE
               IF( WANTZ ) THEN
                  DO 100 J = LILOZ, LIHIZ
                     SUM = Z( J, ICOL1 ) + V2*Z( J, ICOL1+1 ) +
     $                     V3*Z( J, ICOL1+2 )
                     Z( J, ICOL1 ) = Z( J, ICOL1 ) - SUM*T1
                     Z( J, ICOL1+1 ) = Z( J, ICOL1+1 ) - SUM*T2
                     Z( J, ICOL1+2 ) = Z( J, ICOL1+2 ) - SUM*T3
  100             CONTINUE
               END IF
               ICOL1 = ICOL1 + 1
  110       CONTINUE
         ELSE
            DO 120 J = ITMP1, ITMP2
               SUM = A( J, ICOL1 ) + V2*A( J, ICOL1+1 ) +
     $               V3*A( J, ICOL1+2 )
               A( J, ICOL1 ) = A( J, ICOL1 ) - SUM*T1
               A( J, ICOL1+1 ) = A( J, ICOL1+1 ) - SUM*T2
               A( J, ICOL1+2 ) = A( J, ICOL1+2 ) - SUM*T3
  120       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of SLAREF
*
      END
