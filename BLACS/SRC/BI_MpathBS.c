#include "Bdef.h"

void BI_MpathBS(BLACSCONTEXT *ctxt, BLACBUFF *bp, SDRVPTR send, Int npaths)
{
   Int pathlen;		/* the length of each path */
   Int dist;	        /* the distance to the node closest to src on each path */
   Int pdest;           /* part of dest calculation -- saves unneeded ops */
   Int lastlong;	/* number of paths with extra node */
   Int Np, Iam, msgid, Np_1, dir;

   Np = ctxt->scp->Np;
   if (Np < 2) return;
   Iam = ctxt->scp->Iam;
   msgid = Mscopeid(ctxt);
   Np_1 = Np - 1;
   if (npaths == FULLCON) npaths = Np_1;

   if (npaths > 0)  /* paths are increasing rings */
   {
      pdest = Iam;
      dir = 1;
   }
   else             /* paths are decreasing rings */
   {
      pdest = Np + Iam;
      dir = -1;
      npaths = -npaths;
   }
/*
 * Ensure npaths is correct
 */
   if (npaths > Np_1) npaths = Np_1;
   pathlen = Np_1 / npaths;

/*
 * Loop over all long paths (paths with an extra node), if there are any
 */
   lastlong = (Np_1 % npaths) * (pathlen+1);  /* last node in long ring */
   for (dist=1; dist < lastlong; dist += pathlen+1)
      send(ctxt, (pdest+dir*dist)%Np, msgid, bp);

/*
 * Loop over all normal length paths
 */
   while (dist < Np)
   {
      send(ctxt, (pdest+dir*dist)%Np, msgid, bp);
      dist += pathlen;
   }
}
