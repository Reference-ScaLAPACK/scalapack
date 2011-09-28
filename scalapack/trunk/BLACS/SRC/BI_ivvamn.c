#include "Bdef.h"
void BI_ivvamn(int N, char *vec1, char *vec2)
{
   int *v1=(int*)vec1, *v2=(int*)vec2;
   int diff;
   BI_DistType *dist1, *dist2;
   int i, k;

   k = N * sizeof(int);
   i = k % sizeof(BI_DistType);
   if (i) k += sizeof(BI_DistType) - i;
   dist1 = (BI_DistType *) &vec1[k];
   dist2 = (BI_DistType *) &vec2[k];

   for (k=0; k < N; k++)
   {
      diff = Rabs(v1[k]) - Rabs(v2[k]);
      if (diff > 0)
      {
         v1[k] = v2[k];
         dist1[k] = dist2[k];
      }
      else if (diff == 0)
      {
         if (dist1[k] > dist2[k])
         {
            v1[k] = v2[k];
            dist1[k] = dist2[k];
         }
      }
   }
}
