#include "Bdef.h"

void BI_ivvamn2(Int N, char *vec1, char *vec2)
{
   Int k;
   Int *v1=(Int*)vec1, *v2=(Int*)vec2;
   Int diff;

   for (k=0; k != N; k++)
   {
      diff = Rabs(v1[k]) - Rabs(v2[k]);
      if (diff > 0) v1[k] = v2[k];
      else if (diff == 0) if (v1[k] < v2[k]) v1[k] = v2[k];
   }
}
