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
#include "../pblas.h"
#include "../PBpblas.h"
#include "../PBtools.h"
#include "../PBblacs.h"
#include "../PBblas.h"

#ifdef __STDC__
char * PB_Cmalloc( Int LENGTH )
#else
char * PB_Cmalloc( LENGTH )
/*
*  .. Scalar Arguments ..
*/
   Int            LENGTH;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cmalloc allocates a dynamic memory buffer. In case of failure, the
*  program is stopped by calling Cblacs_abort.
*
*  Arguments
*  =========
*
*  LENGTH  (local input) INTEGER
*          On entry, LENGTH  specifies the length in bytes of the buffer
*          to be allocated.  If LENGTH is less or equal than zero,  this
*          function returns NULL.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           * bufptr = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   if( LENGTH > 0 )
   {
      if( !( bufptr = (char *) malloc( (unsigned)LENGTH ) ) )
      {
         (void) fprintf( stderr, "Not enough memory on line %d of file %s!!\n",
                         __LINE__, __FILE__ );
         Cblacs_abort( -1, -1 );
      }
   }
   return( bufptr );
/*
*  End of PB_Cmalloc
*/
}
