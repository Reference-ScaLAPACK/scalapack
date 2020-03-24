#include "Bdef.h"

void BI_SringBS(BLACSCONTEXT *ctxt, BLACBUFF *bp, SDRVPTR send)
{
   Int Np, Iam, msgid;

   Np = ctxt->scp->Np;
   if (Np < 2) return;
   Iam = ctxt->scp->Iam;
   msgid = Mscopeid(ctxt);
   send(ctxt, (Iam + 1)%Np, msgid, bp);
   if (Np > 2) send(ctxt, (Np + Iam - 1)%Np, msgid, bp);
}
