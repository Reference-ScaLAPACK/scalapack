#include "tools.h"

char * getpbbuf( char *mess, Int length )
{
/*
*  Purpose
*  =======
*
*  getpbbuf returns a pointer to a working buffer of size length alloca-
*  ted for the PBLAS routines.
*
* ======================================================================
*
*  .. Local Scalars ..
*/
   static char * pblasbuf = NULL;
   static Int  pbbuflen = 0, mone = -1;
/* ..
*  .. External Functions ..
*/
   void        blacs_abort_();
/* ..
*  .. Executable Statements ..
*/
   if( length >= 0 )
   {
      if( length > pbbuflen )
      {
         if( pblasbuf )
            free( pblasbuf );
         pblasbuf = (char *) malloc((unsigned)length);
         if( !pblasbuf )
         {
            fprintf( stderr,
                     "PBLAS %s ERROR: Memory allocation failed\n",
                     mess );
            blacs_abort_( &mone, &mone );
         }
         pbbuflen = length;
      }
   }
   else if( pblasbuf )
   {
      free( pblasbuf );
      pblasbuf = NULL;
      pbbuflen = 0;
   }
   return( pblasbuf );
}
