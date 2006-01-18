      SUBROUTINE PXERBLA( ICTXT, SRNAME, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, INFO
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      SRNAME
*     ..
*
*  Purpose
*  =======
*
*  PXERBLA is an error handler for the ScaLAPACK routines.  It is called
*  by a ScaLAPACK routine if an input parameter has an invalid value.  A
*  message is printed. Installers may consider modifying this routine in
*  order to call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  SRNAME  (global input) CHARACTER*(*)
*          On entry, SRNAME specifies the name of the routine which cal-
*          ling PXERBLA.
*
*  INFO    (global input) INTEGER
*          On entry, INFO  specifies the position of the invalid parame-
*          ter in the parameter list of the calling routine.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      WRITE( *, FMT = 9999 ) MYROW, MYCOL, SRNAME, INFO
*
 9999 FORMAT( '{', I5, ',', I5, '}:  On entry to ', A,
     $        ' parameter number ', I4, ' had an illegal value' )
*
      RETURN
*
*     End of PXERBLA
*
      END
