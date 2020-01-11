#include "Bdef.h"

void BI_IdringBR(BLACSCONTEXT *ctxt, BLACBUFF *bp, SDRVPTR send, Int src, Int step)
{
   void BI_Srecv(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   Int Np, Iam, msgid, dest;

   Np = ctxt->scp->Np;
   Iam = ctxt->scp->Iam;
   dest = (Np + Iam + step) % Np;
   msgid = Mscopeid(ctxt);
   BI_Srecv(ctxt, BANYNODE, msgid, bp);
   if (dest != src) send(ctxt, dest, msgid, bp);
}
