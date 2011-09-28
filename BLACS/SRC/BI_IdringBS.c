#include "Bdef.h"

void BI_IdringBS(BLACSCONTEXT *ctxt, BLACBUFF *bp, SDRVPTR send, int step)
{
   int Np, Iam, msgid;

   Np = ctxt->scp->Np;
   if (Np < 2) return;
   Iam = ctxt->scp->Iam;
   msgid = Mscopeid(ctxt);

   send(ctxt, (Np+Iam+step)%Np, msgid, bp);
}
