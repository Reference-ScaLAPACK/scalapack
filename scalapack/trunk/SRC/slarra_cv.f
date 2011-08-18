      SUBROUTINE SLARRA_CV( N, D, E, E2, SPLTOL, TNRM, 
     $                    NSPLIT, ISPLIT, INFO )
*
*  -- LAPACK auxiliary routine (version *TBA*) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     July 4, 2010
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N, NSPLIT
      REAL                SPLTOL, TNRM
*     ..
*     .. Array Arguments ..
      INTEGER            ISPLIT( * )
      REAL               D( * ), E( * ), E2( * )
*
*  Purpose
*  =======
*
*  Compute the splitting points with threshold SPLTOL.
*  SLARRA sets any "small" off-diagonal elements to zero.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix. N > 0.
*
*  D       (input) REAL             array, dimension (N)
*          On entry, the N diagonal elements of the tridiagonal
*          matrix T.
*
*  E       (input/output) REAL             array, dimension (N)
*          On entry, the first (N-1) entries contain the subdiagonal
*          elements of the tridiagonal matrix T; E(N) need not be set.
*          On exit, the entries E( ISPLIT( I ) ), 1 <= I <= NSPLIT, 
*          are set to zero, the other entries of E are untouched.
*
*  E2      (input/output) REAL             array, dimension (N)
*          On entry, the first (N-1) entries contain the SQUARES of the 
*          subdiagonal elements of the tridiagonal matrix T; 
*          E2(N) need not be set.
*          On exit, the entries E2( ISPLIT( I ) ),
*          1 <= I <= NSPLIT, have been set to zero
*
*  SPLTOL (input) REAL            
*          The threshold for splitting. Two criteria can be used:
*          SPLTOL<0 : criterion based on absolute off-diagonal value
*          SPLTOL>0 : criterion that preserves relative accuracy
*
*  TNRM (input) REAL            
*          The norm of the matrix.
*
*  NSPLIT  (output) INTEGER
*          The number of blocks T splits into. 1 <= NSPLIT <= N.
*
*  ISPLIT  (output) INTEGER array, dimension (N)
*          The splitting points, at which T breaks up into blocks.
*          The first block consists of rows/columns 1 to ISPLIT(1),
*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
*          etc., and the NSPLIT-th consists of rows/columns
*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
*
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Beresford Parlett, University of California, Berkeley, USA
*     Jim Demmel, University of California, Berkeley, USA
*     Inderjit Dhillon, University of Texas, Austin, USA
*     Osni Marques, LBNL/NERSC, USA
*     Christof Voemel, University of California, Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      REAL               EABS, STOL, TMP1

*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      INFO = 0

*     Compute splitting points
      NSPLIT = 1
      IF(SPLTOL.LT.ZERO) THEN
*        Criterion based on absolute off-diagonal value
         TMP1 = ABS(SPLTOL)* TNRM
         DO 9 I = 1, N-1
            EABS = ABS( E(I) )
            IF( EABS .LE. TMP1) THEN
               E(I) = ZERO
               E2(I) = ZERO
               ISPLIT( NSPLIT ) = I
               NSPLIT = NSPLIT + 1
            END IF
 9       CONTINUE
      ELSE
*        Criterion that guarantees relative accuracy
         STOL = ABS(SPLTOL)
         DO 10 I = 1, N-1
            EABS = ABS( E(I) )
            IF( EABS .LE. STOL * SQRT(ABS(D(I)))*SQRT(ABS(D(I+1))) ) 
     $      THEN
               E(I) = ZERO
               E2(I) = ZERO
               ISPLIT( NSPLIT ) = I         
               NSPLIT = NSPLIT + 1 
            END IF                 
 10      CONTINUE      
      ENDIF
      ISPLIT( NSPLIT ) = N

      RETURN
*     
*     End of SLARRA
*
      END
