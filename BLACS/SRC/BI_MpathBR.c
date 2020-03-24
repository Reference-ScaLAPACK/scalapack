#include "Bdef.h"

void BI_MpathBR(BLACSCONTEXT *ctxt, BLACBUFF *bp, SDRVPTR send, Int src, Int npaths)
{
   void BI_Arecv(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   Int BI_BuffIsFree(BLACBUFF *, Int);

   Int pathlen;		/* the minimal length of each path */
   Int mydist;		/* my distance from src */
   Int faredge;		/* node at far end of path */
   Int lastlong;	/* distance to node on end of last path with extra node */
   Int Np, Iam, msgid, Np_1, dest;

   msgid = Mscopeid(ctxt);
   BI_Arecv(ctxt, BANYNODE, msgid, bp);
   Np = ctxt->scp->Np;
   Iam = ctxt->scp->Iam;
   Np_1 = Np - 1;
   if (npaths == FULLCON) npaths = Np_1;

   if (npaths > 0)
   {
      dest = (Iam+1) % Np;
      mydist = (Np + Iam - src) % Np;
   }
   else
   {
      dest = (Np_1+Iam) % Np;
      mydist = (Np + src - Iam) % Np;
      npaths = -npaths;
   }
/*
 * Make sure npaths is cool
 */
   if (npaths > Np_1) npaths = Np_1;

   pathlen = Np_1 / npaths;
   lastlong = (Np_1%npaths) * (pathlen+1);
   if (lastlong)
   {
      if (mydist <= lastlong) faredge = ((mydist-1)/(pathlen+1)+1)*(pathlen+1);
      else faredge = ((lastlong-1)/(pathlen+1)+1) * (pathlen+1)
		     + ((mydist-lastlong-1)/pathlen + 1) * pathlen;
   }
   else faredge = ((mydist-1)/pathlen + 1) * pathlen;

   BI_BuffIsFree(bp, 1);	/* wait for recv to complete */
   if (mydist < faredge) send(ctxt, dest, msgid, bp);
}
