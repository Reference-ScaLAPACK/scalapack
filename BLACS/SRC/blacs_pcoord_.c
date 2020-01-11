#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_pcoord(Int ConTxt, Int nodenum, Int *prow, Int *pcol)
#else
F_VOID_FUNC blacs_pcoord_(Int *ConTxt, Int *nodenum, Int *prow, Int *pcol)
#endif
{
   BLACSCONTEXT *ctxt;

   MGetConTxt(Mpval(ConTxt), ctxt);
   if ( (Mpval(nodenum) >= 0) && (Mpval(nodenum) < ctxt->ascp.Np) )
   {
      Mpcoord(ctxt, Mpval(nodenum), *prow, *pcol);
   }
   else *prow = *pcol = -1;
}
