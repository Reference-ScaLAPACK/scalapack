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
void PB_Cabort( Int ICTXT, char * ROUT, Int INFO )
#else
void PB_Cabort( ICTXT, ROUT, INFO )
/*
*  .. Scalar Arguments ..
*/
   Int            ICTXT, INFO;
/*
*  .. Array Arguments ..
*/
   char           * ROUT;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cabort is an error handler for the PBLAS  routines.  This  routine
*  displays an error message on  stderr  by calling  PB_Cwarn, and halts
*  execution by calling Cblacs_abort().
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  ROUT    (global input) pointer to CHAR
*          On entry, ROUT specifies the name of the routine calling this
*          error handler.
*
*  INFO    (local input) INTEGER
*          The error code computed by the calling PBLAS routine.
*          = 0:  no error found
*          < 0:  If the  i-th  argument is an array and the j-entry  had
*                an illegal value, then  INFO = -(i*100+j),  if the i-th
*                argument  is  a  scalar  and had an illegal value, then
*                INFO = -i.
*
*  -- Written on April 1, 1998 by
*     R. Clint Whaley, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   Int            mycol, myrow, npcol, nprow;
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
/* ..
*  .. Executable Statements ..
*
*/
   Cblacs_gridinfo( ICTXT, &nprow, &npcol, &myrow, &mycol );
#ifdef TestingPblas
/*
*  For testing purpose only, the error is reported, but the program execution
*  is not terminated
*/
   if( PB_NoAbort( &INFO ) ) return;
#endif
   if( INFO < 0 )
   {
/*
*  Display an error message
*/
      if( INFO < DESCMULT )
         PB_Cwarn( ICTXT, -1, ROUT,
                   "Parameter number %d had an illegal value", -INFO );
      else
         PB_Cwarn( ICTXT, -1, ROUT,
                   "Parameter number %d, entry number %d had an illegal value",
                   (-INFO) / DESCMULT, (-INFO) % DESCMULT );
   }
   else
   {
/*
*  Error code is incorrect, it should be negative
*/
      PB_Cwarn( ICTXT, -1, ROUT,
                "Positive error code %d returned by %s!!!", INFO );
   }
   Cblacs_abort( ICTXT, INFO );
/*
*  End of PB_Cabort
*/
}
