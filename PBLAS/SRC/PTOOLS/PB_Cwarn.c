/* ---------------------------------------------------------------------
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*  ---------------------------------------------------------------------
*/
/*
*  Include files
*/
#ifdef TestingPblas
#include "../SRC/pblas.h"
#include "../SRC/PBpblas.h"
#include "../SRC/PBtools.h"
#include "../SRC/PBblacs.h"
#include "../SRC/PBblas.h"
#else
#include "../pblas.h"
#include "../PBpblas.h"
#include "../PBtools.h"
#include "../PBblacs.h"
#include "../PBblas.h"
#endif

/*
*  ---------------------------------------------------------------------
*  FORTRAN <-> C interface
*  ---------------------------------------------------------------------
*
*  These macros identifies how the PBLAS will be called as follows:
*
*  _F2C_ADD_: the FORTRAN compiler expects the name of C functions to be
*  in all lower case and to have an underscore postfixed it (Suns, Intel
*  compilers expect this).
*
*  _F2C_NOCHANGE: the FORTRAN compiler expects the name of  C  functions
*  to be in all lower case (IBM RS6K compilers do this).
*
*  _F2C_UPCASE: the  FORTRAN  compiler expects the name of  C  functions
*  to be in all upcase. (Cray compilers expect this).
*
*  _F2C_F77ISF2C: the  FORTRAN  compiler in use is f2c, a  FORTRAN  to C
*  converter.
*/
#if (_F2C_CALL_ == _F2C_ADD_ )
#define PB_NoAbort pb_noabort_
#endif
#if (_F2C_CALL_ == _F2C_UPCASE )
#define PB_NoAbort PB_NOABORT
#endif
#if (_F2C_CALL_ == _F2C_NOCHANGE )
#define PB_NoAbort pb_noabort
#endif
#if (_F2C_CALL_ == _F2C_F77ISF2C )
#define PB_NoAbort pb_noabort__
#endif

#ifdef __STDC__
void PB_Cwarn( Int ICTXT, Int LINE, char * ROUT, char * FORM, ... )
#else
void PB_Cwarn( va_alist )
va_dcl
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cwarn  is  an error handler for the PBLAS routines.  This  routine
*  displays an error message on stderr.
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  LINE    (local input) INTEGER
*          On entry,  LINE  specifies the line  number in the file where
*          the error has occured. When  LINE is not a valid line number,
*
*  ROUT    (global input) pointer to CHAR
*          On entry, ROUT specifies the name of the routine calling this
*          error handler.
*
*  FORM    (local input) pointer to CHAR
*          On entry,  FORM  is a  control  string  specifying the format
*          conversion of its following arguments.
*
*  ...     (local input)
*          On entry,  FORM  is a  control  string  specifying the format
*          On entry,  the expressions that are to be  evaluated and con-
*          verted  according  to the formats in the control string  FORM
*          and then placed in the output stream.
*
*  -- Written on April 1, 1998 by
*     R. Clint Whaley, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
   va_list        argptr;
   Int            iam, mycol, myrow, npcol, nprow;
   char           cline[100];
/* ..
*  .. External Functions ..
*/
#ifdef TestingPblas
#ifdef __STDC__
   Int            PB_NoAbort( Int * );
#else
   Int            PB_NoAbort();
#endif
#endif

#ifdef __STDC__
   va_start( argptr, FORM );
#else
   char           * ROUT, * FORM;
   Int            ICTXT, LINE;
/* ..
*  .. Executable Statements ..
*
*/
   va_start( argptr );
   ICTXT = va_arg( argptr, Int );
   LINE  = va_arg( argptr, Int );
   ROUT  = va_arg( argptr, char * );
   FORM  = va_arg( argptr, char * );
#endif

#ifdef TestingPblas
/*
*  For testing purpose only, the error is reported, but the program execution
*  is not terminated
*/
   if( PB_NoAbort( &ICTXT ) ) return;
#endif
   vsprintf( cline, FORM, argptr );
   va_end( argptr );

   Cblacs_gridinfo( ICTXT, &nprow, &npcol, &myrow, &mycol );

   if( nprow != -1 ) iam = Cblacs_pnum( ICTXT, myrow, mycol );
   else              iam = -1;
/*
*  Display an error message
*/
   if( LINE <= 0 )
      (void) fprintf( stderr, "%s'%s'\n%s{%d,%d}, %s%d, %s%d%s'%s'.\n\n",
                      "PBLAS ERROR ", cline, "from ", myrow, mycol, "pnum=",
                      iam, "Contxt=", ICTXT, ", in routine ", ROUT );
   else
      (void) fprintf( stderr, "%s'%s'\n%s{%d,%d}, %s%d, %s%d%s%d%s'%s'.\n\n",
                      "PBLAS ERROR ", cline, "from ", myrow, mycol, "pnum=",
                      iam, "Contxt=", ICTXT, ", on line ", LINE,
                      " of routine ", ROUT );
/*
*  End of PB_Cwarn
*/
}
