#include "Bdef.h"
void BI_cvvamx(int N, char *vec1, char *vec2)
{
   SCOMPLEX *v1=(SCOMPLEX*)vec1, *v2=(SCOMPLEX*)vec2;
   float diff;
   BI_DistType *dist1, *dist2;
   int i, k;

   k = N * sizeof(SCOMPLEX);
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
