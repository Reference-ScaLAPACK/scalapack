#include "Bdef.h"
void BI_MringComb(BLACSCONTEXT *ctxt, BLACBUFF *bp, BLACBUFF *bp2,
                  Int N, VVFUNPTR Xvvop, Int dest, Int nrings)
{
   void BI_Ssend(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   void BI_Srecv(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   void BI_MpathBS(BLACSCONTEXT *, BLACBUFF *, SDRVPTR, Int);
   void BI_MpathBR(BLACSCONTEXT *, BLACBUFF *, SDRVPTR, Int, Int);

   Int Np, Iam, msgid, i, inc, mysrc, mydest, Np_1;
   Int mydist, ringlen, myring;
   Int nearedge, faredge;  /* edge closest and farthest from dest */
   Int REBS;               /* Is result leave-on-all? */

   Np = ctxt->scp->Np;
   if (Np < 2) return;
   Iam = ctxt->scp->Iam;
   msgid = Mscopeid(ctxt);
   if (REBS = (dest == -1)) dest = 0;

   if (nrings > 0)
   {
      mydist = (Np + dest - Iam) % Np;
      inc = 1;
   }
   else
   {
      mydist = (Np + Iam - dest) % Np;
      inc = -1;
      nrings = -nrings;
   }
   Np_1 = Np - 1;
   if (nrings > Np_1) nrings = Np_1;

/*
 * If I'm not the destination
 */
   if (Iam != dest)
   {
      ringlen = Np_1 / nrings;
      myring = (mydist-1) / ringlen;
      if (myring >= nrings) myring = nrings - 1;
      nearedge = (myring*ringlen) + 1;
      faredge = nearedge + ringlen - 1;
      if (myring == nrings-1) faredge += Np_1 % nrings;
      if (mydist == nearedge) mydest = dest;
      else mydest = (Np + Iam + inc) % Np;
      if (mydist != faredge)
      {
         BI_Srecv(ctxt, (Np + Iam - inc) % Np, msgid, bp2);
	 Xvvop(N, bp->Buff, bp2->Buff);
      }
      BI_Ssend(ctxt, mydest, msgid, bp);
      if (REBS) BI_MpathBR(ctxt, bp, BI_Ssend, dest, nrings);
   }
/*
 * If I'm the destination process
 */
   else
   {
      if (!ctxt->TopsRepeat)
      {
         for(i=nrings; i; i--)
         {
            BI_Srecv(ctxt, BANYNODE, msgid, bp2);
	    Xvvop(N, bp->Buff, bp2->Buff);
         }
      }
      else
      {
         ringlen = Np_1 / nrings;
         if (inc == 1) mysrc = (Np + Iam - 1) % Np;
         else mysrc = (Iam + 1) % Np;
         for(i=nrings; i; i--)
         {
            BI_Srecv(ctxt, mysrc, msgid, bp2);
	    Xvvop(N, bp->Buff, bp2->Buff);
            if (inc == 1) mysrc = (Np + mysrc - ringlen) % Np;
            else mysrc = (mysrc + ringlen) % Np;
         }
      }
      if (REBS) BI_MpathBS(ctxt, bp, BI_Ssend, nrings);
   }
}  /* end BI_MringComb */
