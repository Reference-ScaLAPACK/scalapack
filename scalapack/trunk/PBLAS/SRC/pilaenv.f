      INTEGER FUNCTION PILAENV( ICTXT, PREC )
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT
      CHARACTER*1        PREC
*     ..
*
*  Purpose
*  =======
*
*  PILAENV  returns  the  positive integer value of the logical blocking
*  size. This value is machine and precision specific. This version pro-
*  vides a logical blocking size which should give good but  not optimal
*  performance on many  of the  currently  available  distributed-memory
*  concurrent computers.  Users are encouraged to modify this subroutine
*  to set this tuning parameter for their particular machine.
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  PREC    (global input) CHARACTER*1
*          On input, PREC specifies the precision for which the  logical
*          block size should be returned as follows:
*             PREC = 'S' or 's' single precision real,
*             PREC = 'D' or 'd' double precision real,
*             PREC = 'C' or 'c' single precision complex,
*             PREC = 'Z' or 'z' double precision complex,
*             PREC = 'I' or 'i' integer.
*
*  Notes
*  =====
*
*  Before modifying this routine to tune the library performance on your
*  system, be aware of the following:
*
*  1) The value this function returns must be STRICTLY LARGER THAN ZERO,
*
*  2) If you are planning to link your program with different  instances
*  of the library, (for example on a heterogeneous  machine),  you  MUST
*  compile each instance of the library with the  EXACT SAME  version of
*  this routine for obvious inter-operability reasons.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( PREC, 'S' ) ) THEN
*
*        Single precision real logical block size
*
         PILAENV = 32
*
      ELSE IF( LSAME( PREC, 'D' ) ) THEN
*
*        Double precision real logical block size
*
         PILAENV = 32
*
      ELSE IF( LSAME( PREC, 'C' ) ) THEN
*
*        Single precision complex logical block size
*
         PILAENV = 32
*
      ELSE IF( LSAME( PREC, 'Z' ) ) THEN
*
*        Double precision complex logical block size
*
         PILAENV = 32
*
      ELSE IF( LSAME( PREC, 'I' ) ) THEN
*
*        Integer logical block size
*
         PILAENV = 32
*
      ELSE
*
*        Probably unused
*
         PILAENV = 32
*
      END IF
*
      RETURN
*
*     End of PILAENV
*
      END
