#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_pcoord(int ConTxt, int nodenum, int *prow, int *pcol)
#else
F_VOID_FUNC blacs_pcoord_(int *ConTxt, int *nodenum, int *prow, int *pcol)
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
