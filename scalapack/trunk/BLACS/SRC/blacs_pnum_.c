#include "Bdef.h"

#if (INTFACE == C_CALL)
int Cblacs_pnum(int ConTxt, int prow, int pcol)
#else
F_INT_FUNC blacs_pnum_(int *ConTxt, int *prow, int *pcol)
#endif
{
   BLACSCONTEXT *ctxt;

   MGetConTxt(Mpval(ConTxt), ctxt);
   if ( (Mpval(prow) >= 0) && (Mpval(prow) < ctxt->cscp.Np) &&
        (Mpval(pcol) >= 0) && (Mpval(pcol) < ctxt->rscp.Np) )
      return( Mkpnum(ctxt, Mpval(prow), Mpval(pcol)) );
   else return(-1);
}
