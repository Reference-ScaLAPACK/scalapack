#include "Bdef.h"

void BI_TransDist(BLACSCONTEXT *ctxt, char scope, int m, int n, int *rA,
                  int *cA, int ldrc, BI_DistType *dist, int rdest, int cdest)
/*
 *  This routine translates distances (offsets from the destination node),
 *  stored in location dist, into row and column coordinates.
 */
{
   int i, j, k, dest;
   int Ng, nprow, npcol, myrow, mycol;

   Mgridinfo(ctxt, Ng, nprow, npcol, myrow, mycol);
   if (rdest == -1) rdest = cdest = 0;

   switch (scope)
   {
   case 'r':
      for (j=0; j < n; j++)
      {
         for (i=0; i < m; i++)
         {
            rA[i] = myrow;
            cA[i] = (int) (cdest + dist[i]) % npcol;
         }
         rA += ldrc;
         cA += ldrc;
         dist += m;
       }
       break;
   case 'c':
      for (j=0; j < n; j++)
      {
         for (i=0; i < m; i++)
         {
            rA[i] = (int) (rdest + dist[i]) % nprow;
            cA[i] = mycol;
         }
         rA += ldrc;
         cA += ldrc;
         dist += m;
      }
      break;
   case 'a':
      dest = Mvkpnum(ctxt, rdest, cdest);
      for (j=0; j < n; j++)
      {
         for (i=0; i < m; i++)
         {
            k = (int) (dest + dist[i]) % Ng;   /* figure node number */
            Mvpcoord(ctxt, k, rA[i], cA[i]);   /* figure node coordinates */
         }
         rA += ldrc;
         cA += ldrc;
         dist += m;
      }
   }
}
