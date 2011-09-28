#include "Bdef.h"

void BI_IdringBR(BLACSCONTEXT *ctxt, BLACBUFF *bp, SDRVPTR send, int src, int step)
{
   void BI_Srecv(BLACSCONTEXT *, int, int, BLACBUFF *);
   int Np, Iam, msgid, dest;

   Np = ctxt->scp->Np;
   Iam = ctxt->scp->Iam;
   dest = (Np + Iam + step) % Np;
   msgid = Mscopeid(ctxt);
   BI_Srecv(ctxt, BANYNODE, msgid, bp);
   if (dest != src) send(ctxt, dest, msgid, bp);
}
