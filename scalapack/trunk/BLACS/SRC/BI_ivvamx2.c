#include "Bdef.h"
void BI_ivvamx2(int N, char *vec1, char *vec2)
{
   int k;
   int *v1=(int*)vec1, *v2=(int*)vec2;
   int diff;

   for (k=0; k != N; k++)
   {
      diff = Rabs(v1[k]) - Rabs(v2[k]);
      if (diff < 0) v1[k] = v2[k];
      else if (diff == 0) if (v1[k] < v2[k]) v1[k] = v2[k];
   }
}
