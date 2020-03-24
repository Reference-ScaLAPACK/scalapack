#include "Bdef.h"

void BI_TreeBR(BLACSCONTEXT *ctxt, BLACBUFF *bp, SDRVPTR send, Int src, Int nbranches)
{
   void BI_Srecv(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   Int Np, Iam, msgid, i, j;
   Int mydist;          /* my distance from src */
   Int destdist;	/* the distance of the destination node */

   Np = ctxt->scp->Np;
   if (Np < 2) return;
   Iam = ctxt->scp->Iam;
   msgid = Mscopeid(ctxt);
   mydist = (Np + Iam - src) % Np;

/*
 * Go up to first step of tree where I send data to other nodes
 */
   for (i=nbranches; i < Np; i *= nbranches);
   for (i /= nbranches; (mydist%i); i /= nbranches);
   BI_Srecv(ctxt, BANYNODE, msgid, bp);

/*
 * While I need to send data to others
 */
   while ( (i > 1) && !(mydist%i) )
   {
      i /= nbranches;
      j = 1;
      do
      {
	 destdist = mydist + j*i;
	 if (destdist < Np)
            send(ctxt, (src+destdist)%Np, msgid, bp);
      }
      while(++j < nbranches);
   }

} /* end BI_TreeBR */
