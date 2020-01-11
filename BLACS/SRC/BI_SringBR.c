#include "Bdef.h"

void BI_SringBR(BLACSCONTEXT *ctxt, BLACBUFF *bp, SDRVPTR send, Int src)
{
   void BI_Srecv(BLACSCONTEXT *, Int, Int, BLACBUFF *);

   Int mydist;  	/* my distance from source */
   Int Np, Iam, msgid, rightedge;

   Np = ctxt->scp->Np;
   Iam = ctxt->scp->Iam;
   msgid = Mscopeid(ctxt);

   mydist = (Np + Iam - src) % Np;
   rightedge = Np/2;
   BI_Srecv(ctxt, BANYNODE, msgid, bp);

/*
 * If I'm between source & right edge of split ring, send to right
 */
   if (mydist < rightedge)
      send(ctxt, (Iam+1)%Np, msgid, bp);
/*
 * If I'm between source and left edge of split ring, send to left
 */
   else if (mydist > rightedge+1)
      send(ctxt, (Np+Iam-1)%Np, msgid, bp);
}
