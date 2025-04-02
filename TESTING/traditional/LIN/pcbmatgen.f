      SUBROUTINE PCBMATGEN( ICTXT, AFORM, AFORM2, BWL, BWU, N,
     $                     MB, NB, A,
     $                     LDA, IAROW, IACOL, ISEED,
     $                     MYROW, MYCOL, NPROW, NPCOL )
*
*
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 15, 1997
*
*     .. Scalar Arguments ..
*     .. Scalar Arguments ..
      CHARACTER*1        AFORM, AFORM2
      INTEGER            IACOL, IAROW, ICTXT,
     $                   ISEED, LDA, MB, MYCOL, MYROW, N,
     $                   NB, NPCOL, NPROW, BWL, BWU
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  PCBMATGEN : Parallel Complex Single precision Band MATrix GENerator.
*  (Re)Generate a distributed Band matrix A (or sub-matrix of A).
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation. The context itself is global.
*
*  AFORM   (global input) CHARACTER*1
*          if AFORM = 'L' : A is returned as a hermitian lower
*            triangular matrix, and is diagonally dominant.
*          if AFORM = 'U' : A is returned as a hermitian upper
*            triangular matrix, and is diagonally dominant.
*          if AFORM = 'G' : A is returned as a general matrix.
*          if AFORM = 'T' : A is returned as a general matrix in
*            tridiagonal-compatible form.
*
*  AFORM2  (global input) CHARACTER*1
*          if the matrix is general:
*            if AFORM2 = 'D' : A is returned diagonally dominant.
*            if AFORM2 != 'D' : A is not returned diagonally dominant.
*          if the matrix is symmetric or hermitian:
*            if AFORM2 = 'T' : A is returned in tridiagonally-compatible
*              form (a transpose form).
*            if AFORM2 != 'T' : A is returned in banded-compatible form.
*
*  M       (global input) INTEGER
*          The number of nonzero rows in the generated distributed
*           band matrix.
*
*  N       (global input) INTEGER
*          The number of columns in the generated distributed
*          matrix.
*
*  MB      (global input) INTEGER
*          The row blocking factor of the distributed matrix A.
*
*  NB      (global input) INTEGER
*          The column blocking factor of the distributed matrix A.
*
*  A       (local output) COMPLEX, pointer into the local memory to
*          an array of dimension ( LDA, * ) containing the local
*          pieces of the distributed matrix.
*
*  LDA     (local input) INTEGER
*          The leading dimension of the array containing the local
*          pieces of the distributed matrix A.
*
*  IAROW   (global input) INTEGER
*          The row processor coordinate which holds the first block
*          of the distributed matrix A.
*            A( DIAG_INDEX, I ) = A( DIAG_INDEX, I ) + BWL+BWU
*
*  IACOL   (global input) INTEGER
*          The column processor coordinate which holds the first
*          block of the distributed matrix A.
*
*  ISEED   (global input) INTEGER
*          The seed number to generate the distributed matrix A.
*
*  MYROW   (local input) INTEGER
*          The row process coordinate of the calling process.
*
*  MYCOL   (local input) INTEGER
*          The column process coordinate of the calling process.
*
*  NPROW   (global input) INTEGER
*          The number of process rows in the grid.
*
*  NPCOL   (global input) INTEGER
*          The number of process columns in the grid.
*
*  Notes
*  =====
*
*  This code is a simple wrapper around PCMATGEN, for band matrices.
*
*  =====================================================================
*
*  Code Developer: Andrew J. Cleary, University of Tennessee.
*    Current address: Lawrence Livermore National Labs.
*  This version released: August, 2001.
*
*  =====================================================================
*
*     ..
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0 )
      PARAMETER          ( ZERO = 0.0E+0 )
      COMPLEX            CONE, CZERO
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER           DIAG_INDEX, I, J, M_MATGEN, NQ, N_MATGEN,
     $                  START_INDEX
*     ..
*     .. External Subroutines ..
      EXTERNAL           PCMATGEN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, NUMROC
      EXTERNAL           ICEIL, NUMROC, LSAME
*     ..
*     .. Executable Statements ..
*
*
      IF( LSAME( AFORM, 'L' ).OR.LSAME( AFORM, 'U' ) ) THEN
         M_MATGEN = BWL + 1
         N_MATGEN = N
         START_INDEX = 1
         IF( LSAME( AFORM, 'L' ) ) THEN
            DIAG_INDEX = 1
         ELSE
            DIAG_INDEX = BWL + 1
         ENDIF
      ELSE
         M_MATGEN = BWL + BWU + 1
         N_MATGEN = N
         DIAG_INDEX = BWU + 1
         START_INDEX = 1
      ENDIF
*
      NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
*
*
*     Generate a random matrix initially
*
      IF( LSAME( AFORM, 'T' ) .OR.
     $  ( LSAME( AFORM2, 'T' ) ) ) THEN
*
          CALL PCMATGEN( ICTXT, 'T', 'N',
     $                        N_MATGEN, M_MATGEN,
     $                        NB, M_MATGEN, A( START_INDEX, 1 ),
     $                        LDA, IAROW, IACOL,
     $                        ISEED, 0, NQ, 0, M_MATGEN,
     $                        MYCOL, MYROW, NPCOL, NPROW )
*
      ELSE
*
          CALL PCMATGEN( ICTXT, 'N', 'N',
     $                        M_MATGEN, N_MATGEN,
     $                        M_MATGEN, NB, A( START_INDEX, 1 ),
     $                        LDA, IAROW, IACOL,
     $                        ISEED, 0, M_MATGEN, 0, NQ,
     $                        MYROW, MYCOL, NPROW, NPCOL )
*
*        Zero out padding at tops of columns
*
         DO 1000 J=1,NB
*
            DO 2000 I=1, LDA-M_MATGEN
*
*              Indexing goes negative; BMATGEN assumes that space
*              has been preallocated above the first column as it
*              has to be if the matrix is to be input to
*              Scalapack's band solvers.
*
               A( I-LDA+M_MATGEN, J ) = CZERO
*
 2000       CONTINUE
*
 1000    CONTINUE
*
      ENDIF
*
      IF( LSAME( AFORM2, 'D' ).OR.
     $  ( LSAME( AFORM, 'L' ).OR.LSAME( AFORM, 'U' ) ) ) THEN
*
*       Loop over diagonal elements stored on this processor.
*
*
       DO 330 I=1, NQ
         IF( LSAME( AFORM, 'T' ) .OR.
     $     ( LSAME( AFORM2, 'T' ) ) ) THEN
             IF( NPROW .EQ. 1 ) THEN
                A( I, DIAG_INDEX ) = CMPLX( REAL( A( I, DIAG_INDEX ) )
     $                               + REAL( 2*( BWL+BWU+1 ) ) )
             ENDIF
          ELSE
             IF( NPROW .EQ. 1 ) THEN
                A( DIAG_INDEX, I ) = CMPLX( REAL( A( DIAG_INDEX, I ) )
     $                               + REAL( 2*( BWL+BWU+1 ) ) )
             ENDIF
          END IF
  330  CONTINUE
*
*
      ELSE
*
*       Must add elements to keep condition of matrix in check
*
        DO 380 I=1, NQ
*
          IF( NPROW .EQ. 1 ) THEN
*
            IF( MOD(I+MYCOL*NB,2) .EQ. 1 ) THEN
                A( DIAG_INDEX+1, I ) =
     $                         CMPLX( REAL( A( DIAG_INDEX+1, I ) )
     $                         + REAL( 2*( BWL+BWU+1 ) ) )
*
            ELSE
*
                A( DIAG_INDEX-1, I ) =
     $                          CMPLX( REAL( A( DIAG_INDEX-1, I ) )
     $                          + REAL( 2*( BWL+BWU+1 ) ) )
            ENDIF
*
          ENDIF
*
  380   CONTINUE
*
      END IF
*
      RETURN
*
*     End of PCBMATGEN
*
      END
