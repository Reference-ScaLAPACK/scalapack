/* ---------------------------------------------------------------------
*
*  -- PBLAS routine (version 2.0) --
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
void PB_topget_( Int * ICTXT, F_CHAR_T OP, F_CHAR_T SCOPE, F_CHAR_T TOP )
#else
void PB_topget_( ICTXT, OP, SCOPE, TOP )
/*
*  .. Scalar Arguments ..
*/
   Int            * ICTXT;
/*
*  .. Array Arguments ..
*/
   F_CHAR_T       OP, SCOPE, TOP;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_topget_ returns the row, column or all broadcast and combine topo-
*  logies.
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
*  TOP     (global output) pointer to CHAR
*          On exit, TOP  is a character string specifying the  BLACS to-
*          pology current in use for the given operation specified by OP
*          and SCOPE.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           * topc;
/* ..
*  .. Executable Statements ..
*
*/
   topc  = F2C_CHAR( TOP );
   *topc = *PB_Ctop( ICTXT, F2C_CHAR( OP ), F2C_CHAR( SCOPE ), TOP_GET );
   TOP   = C2F_CHAR( topc );
/*
*  End of PB_topget_
*/
}
