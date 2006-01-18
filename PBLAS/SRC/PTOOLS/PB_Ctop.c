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
char * PB_Ctop( int * ICTXT, char * OP, char * SCOPE, char * TOP )
#else
char * PB_Ctop( ICTXT, OP, SCOPE, TOP )
/*
*  .. Scalar Arguments ..
*/
   int            * ICTXT;
/*
*  .. Array Arguments ..
*/
   char           * OP, * SCOPE, * TOP;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Ctop  returns or initializes the row-, column- or all-  broadcast
*  or combine topologies.
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  OP      (global input) pointer to CHAR
*          On input,  OP  specifies  the BLACS operation defined as fol-
*          lows:
*             OP = 'B' or 'b', BLACS broadcast operation,
*             OP = 'C' or 'c', BLACS combine operation.
*
*  SCOPE   (global input) pointer to CHAR
*          On entry, SCOPE specifies the scope of the BLACS operation as
*          follows:
*             SCOPE = 'R' or 'r', rowwise broadcast or combine,
*             SCOPE = 'C' or 'c', column broadcast or combine,
*             SCOPE = 'A' or 'a', all broadcast or combine.
*
*  TOP     (global input) pointer to CHAR
*          On entry, TOP  is a character string specifying the BLACS to-
*          pology to be used i.e. to be set for the given operation spe-
*          cified by OP and SCOPE. If TOP = TOP_GET, the routine instead
*          returns  the  current topology in use for the given operation
*          specified by OP and SCOPE.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   static char    rbtop = CTOP_DEFAULT;
   static char    cbtop = CTOP_DEFAULT;
   static char    abtop = CTOP_DEFAULT;
   static char    rctop = CTOP_DEFAULT;
   static char    cctop = CTOP_DEFAULT;
   static char    actop = CTOP_DEFAULT;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  This BLACS topology information should be cached within a BLACS context.
*  This will be corrected in the near future. Sorry.
*/
   if( *OP == CBCAST )
   {
/*
*  BLACS broadcast operations
*/
      if( *TOP == CTOP_GET )
      {
/*
*  retrieve the current topology in SCOPE
*/
         if( *SCOPE == CROW )         { return( &rbtop ); }
         else if( *SCOPE == CCOLUMN ) { return( &cbtop ); }
         else                         { return( &abtop ); }
      }
      else
      {
/*
*  set the topology to be used from now on in SCOPE
*/
         if( *SCOPE == CROW )         { rbtop = *TOP; return( &rbtop ); }
         else if( *SCOPE == CCOLUMN ) { cbtop = *TOP; return( &cbtop ); }
         else                         { abtop = *TOP; return( &abtop ); }
      }
   }
   else
   {
/*
*  BLACS combine operations
*/
      if( *TOP == CTOP_GET )
      {
/*
*  retrieve the current topology in SCOPE
*/
         if( *SCOPE == CROW )         { return( &rctop ); }
         else if( *SCOPE == CCOLUMN ) { return( &cctop ); }
         else                         { return( &actop ); }
      }
      else
      {
/*
*  set the topology to be used from now on in SCOPE
*/
         if( *SCOPE == CROW )         { rctop = *TOP; return( &rctop ); }
         else if( *SCOPE == CCOLUMN ) { cctop = *TOP; return( &cctop ); }
         else                         { actop = *TOP; return( &actop ); }
      }
   }
/*
*  End of PB_Ctop
*/
}
