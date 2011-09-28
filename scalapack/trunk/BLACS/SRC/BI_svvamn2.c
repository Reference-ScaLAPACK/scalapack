#include "Bdef.h"
void BI_svvamn2(int N, char *vec1, char *vec2)
{
   int k;
   float *v1=(float*)vec1, *v2=(float*)vec2;
   float diff;

   for (k=0; k != N; k++)
   {
      diff = Rabs(v1[k]) - Rabs(v2[k]);
      if (diff > 0) v1[k] = v2[k];
      else if (diff == 0) if (v1[k] < v2[k]) v1[k] = v2[k];
   }
}
