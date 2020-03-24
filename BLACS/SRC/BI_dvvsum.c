#include "Bdef.h"
void BI_dvvsum(Int N, char *vec1, char *vec2)
{
   double *v1=(double*)vec1, *v2=(double*)vec2;
   Int k;
   for (k=0; k < N; k++) v1[k] += v2[k];
}
