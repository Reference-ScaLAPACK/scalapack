#include "Bdef.h"
void BI_zvvamx(Int N, char *vec1, char *vec2)
{
   DCOMPLEX *v1=(DCOMPLEX*)vec1, *v2=(DCOMPLEX*)vec2;
   double diff;
   BI_DistType *dist1, *dist2;
   Int i, k;

   k = N * sizeof(DCOMPLEX);
   i = k % sizeof(BI_DistType);
   if (i) k += sizeof(BI_DistType) - i;
   dist1 = (BI_DistType *) &vec1[k];
   dist2 = (BI_DistType *) &vec2[k];

   for (k=0; k < N; k++)
   {
      diff = Cabs(v1[k]) - Cabs(v2[k]);
      if (diff < 0)
      {
         v1[k].r = v2[k].r;
         v1[k].i = v2[k].i;
         dist1[k] = dist2[k];
      }
      else if (diff == 0)
      {
         if (dist1[k] > dist2[k])
         {
            v1[k].r = v2[k].r;
            v1[k].i = v2[k].i;
            dist1[k] = dist2[k];
         }
      }
   }
}
