#include "Bdef.h"
int BI_HypBS(BLACSCONTEXT *ctxt, BLACBUFF *bp, SDRVPTR send)
{
   int bit, Np, Iam, msgid;

   Np = ctxt->scp->Np;
   if (Np < 2) return(NORV);
   Iam = ctxt->scp->Iam;
   msgid = Mscopeid(ctxt);

   for (bit=2; bit < Np; bit <<= 1);
   if (bit^Np) return(NPOW2);  /* not a power of 2 */

   for(bit=1; (bit^Np); bit <<= 1)
      send(ctxt, (Iam^bit), msgid, bp);

   return(0);  /* error-free return */
}
