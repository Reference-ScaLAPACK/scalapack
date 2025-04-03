*  =====================================================================
*     SUBROUTINE LADD
*  =====================================================================
*
      SUBROUTINE LADD( J, K, I )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Array Arguments ..
      INTEGER            I(2), J(2), K(2)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            IPOW16, IPOW15
      PARAMETER        ( IPOW16=2**16, IPOW15=2**15 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
      I(1) = MOD( K(1)+J(1), IPOW16 )
      I(2) = MOD( (K(1)+J(1)) / IPOW16+K(2)+J(2), IPOW15 )
*
      RETURN
*
*     End of LADD
*
      END
*
*  =====================================================================
*     SUBROUTINE LMUL
*  =====================================================================
*
      SUBROUTINE LMUL( K, J, I )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Array Arguments ..
      INTEGER            I(2), J(2), K(2)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            IPOW15, IPOW16, IPOW30
      PARAMETER        ( IPOW15=2**15, IPOW16=2**16, IPOW30=2**30 )
*     ..
*     .. Local Scalars ..
      INTEGER            KT, LT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
      KT   = K(1)*J(1)
      IF( KT.LT.0 ) KT = (KT+IPOW30) + IPOW30
      I(1) = MOD(KT,IPOW16)
      LT   = K(1)*J(2) + K(2)*J(1)
      IF( LT.LT.0 ) LT = (LT+IPOW30) + IPOW30
      KT   = KT/IPOW16 + LT
      IF( KT.LT.0 ) KT = (KT+IPOW30) + IPOW30
      I(2) = MOD( KT, IPOW15 )
*
      RETURN
*
*     End of LMUL
*
      END
*
*  =====================================================================
*     SUBROUTINE XJUMPM
*  =====================================================================
*
      SUBROUTINE XJUMPM( JUMPM, MULT, IADD, IRANN, IRANM, IAM, ICM )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            JUMPM
*     ..
*     .. Array Arguments ..
      INTEGER            IADD(2), IAM(2), ICM(2), IRANM(2), IRANN(2)
      INTEGER            MULT(2)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. Local Arrays ..
      INTEGER            J(2)
*     ..
*     .. External Subroutines ..
      EXTERNAL           LADD, LMUL
*     ..
*     .. Executable Statements ..
*
      IF( JUMPM.GT.0 ) THEN
         DO 10 I = 1, 2
            IAM(I) = MULT(I)
            ICM(I) = IADD(I)
   10    CONTINUE
         DO 20 I = 1, JUMPM-1
            CALL LMUL( IAM, MULT, J )
            IAM(1) = J(1)
            IAM(2) = J(2)
            CALL LMUL( ICM, MULT, J )
            CALL LADD( IADD, J, ICM )
   20    CONTINUE
         CALL LMUL( IRANN, IAM, J )
         CALL LADD( J, ICM, IRANM )
      ELSE
         IRANM(1) = IRANN(1)
         IRANM(2) = IRANN(2)
      END IF
*
      RETURN
*
*     End of XJUMPM
*
      END
*
*  =====================================================================
*     SUBROUTINE SETRAN
*  =====================================================================
*
      SUBROUTINE SETRAN( IRAN, IA, IC )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Array Arguments ..
      INTEGER            IA(2),  IC(2), IRAN(2)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. Local Arrays ..
      INTEGER            IAS(2),  ICS(2), IRAND(2)
*     ..
*     .. Common Blocks ..
      COMMON /RANCOM/    IRAND, IAS, ICS
      SAVE   /RANCOM/
*     ..
*     .. Executable Statements ..
*
      DO 10 I = 1, 2
         IRAND(I) = IRAN(I)
         IAS(I)   = IA(I)
         ICS(I)   = IC(I)
   10 CONTINUE
*
      RETURN
*
*     End of SETRAN
*
      END
*
*  =====================================================================
*     SUBROUTINE JUMPIT
*  =====================================================================
*
      SUBROUTINE JUMPIT( MULT, IADD, IRANN, IRANM )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Array Arguments ..
      INTEGER            IADD(2), IRANM(2), IRANN(2), MULT(2)
*     ..
*
*  =====================================================================
*
*     .. Local Arrays ..
      INTEGER            IAS(2), ICS(2), IRAND(2), J(2)
*     ..
*     .. External Subroutines ..
      EXTERNAL           LADD, LMUL
*     ..
*     .. Common Blocks ..
      COMMON /RANCOM/    IRAND, IAS, ICS
      SAVE   /RANCOM/
*     ..
*     .. Executable Statements ..
*
      CALL LMUL( IRANN, MULT, J )
      CALL LADD( J, IADD, IRANM )
*
      IRAND(1) = IRANM(1)
      IRAND(2) = IRANM(2)
*
      RETURN
*
*     End of JUMPIT
*
      END
*
*  =====================================================================
*     REAL FUNCTION PSRAND
*  =====================================================================
*
      REAL FUNCTION PSRAND( IDUMM )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IDUMM
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               DIVFAC, POW16
      PARAMETER          ( DIVFAC=2.147483648E+9, POW16=6.5536E+4 )
*     ..
*     .. Local Arrays ..
      INTEGER            J( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           LADD, LMUL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Common Blocks ..
      INTEGER            IAS(2), ICS(2), IRAND(2)
      COMMON /RANCOM/    IRAND, IAS, ICS
      SAVE   /RANCOM/
*     ..
*     .. Executable Statements ..
*
      PSRAND = ( REAL(IRAND(1)) + POW16 * REAL(IRAND(2)) ) / DIVFAC
*
      CALL LMUL( IRAND, IAS, J )
      CALL LADD( J, ICS, IRAND )
*
      RETURN
*
*     End of PSRAND
*
      END
*
*  =====================================================================
*     DOUBLE PRECISION FUNCTION PDRAND
*  =====================================================================
*
      DOUBLE PRECISION FUNCTION PDRAND( IDUMM )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IDUMM
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   DIVFAC, POW16
      PARAMETER          ( DIVFAC=2.147483648D+9, POW16=6.5536D+4 )
*     ..
*     .. Local Arrays ..
      INTEGER            J(2)
*     ..
*     .. External Subroutines ..
      EXTERNAL           LADD, LMUL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Common Blocks ..
      INTEGER            IAS(2), ICS(2), IRAND(2)
      COMMON /RANCOM/    IRAND, IAS, ICS
      SAVE   /RANCOM/
*     ..
*     .. Executable Statements ..
*
      PDRAND = ( DBLE(IRAND(1)) + POW16 * DBLE(IRAND(2)) ) / DIVFAC
*
      CALL LMUL( IRAND, IAS, J )
      CALL LADD( J, ICS, IRAND )
*
      RETURN
*
*     End of PDRAND
*
      END
