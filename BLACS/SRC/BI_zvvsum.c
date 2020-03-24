#include "Bdef.h"
void BI_zvvsum(Int N, char *vec1, char *vec2)
{
   double *v1=(double*)vec1, *v2=(double*)vec2;
   Int k;
   N *=2;
   for (k=0; k < N; k++) v1[k] += v2[k];
}
