#include "Bdef.h"
void BI_ivvsum(int N, char *vec1, char *vec2)
{
   int *v1=(int*)vec1, *v2=(int*)vec2;
   int k;
   for (k=0; k < N; k++) v1[k] += v2[k];
}
