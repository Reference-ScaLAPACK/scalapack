#include "Bdef.h"
void BI_cvvsum(int N, char *vec1, char *vec2)
{
   float *v1=(float*)vec1, *v2=(float*)vec2;
   int k;
   N *=2;
   for (k=0; k < N; k++) v1[k] += v2[k];
}
