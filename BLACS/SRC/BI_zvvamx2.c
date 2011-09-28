#include "Bdef.h"
void BI_zvvamx2(int N, char *vec1, char *vec2)
{
   int r, i;
   double *v1=(double*)vec1, *v2=(double*)vec2;
   double diff;

   N *= 2;
   for (r=0, i=1; r != N; r += 2, i += 2)
   {
      diff = (Rabs(v1[r]) + Rabs(v1[i])) - (Rabs(v2[r]) + Rabs(v2[i]));
      if (diff < 0)
      {
         v1[r] = v2[r];
         v1[i] = v2[i];
      }
      else if (diff == 0)
      {
         if (v1[r] != v2[r])
         {
            if (v1[r] < v2[r])
            {
               v1[r] = v2[r];
               v1[i] = v2[i];
            }
         }
         else
         {
            if (v1[i] < v2[i])
            {
               v1[r] = v2[r];
               v1[i] = v2[i];
            }
         }
      }
   }
}
