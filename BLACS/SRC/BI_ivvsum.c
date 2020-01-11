#include "Bdef.h"
void BI_ivvsum(Int N, char *vec1, char *vec2)
{
   Int *v1=(Int*)vec1, *v2=(Int*)vec2;
   Int k;
   for (k=0; k < N; k++) v1[k] += v2[k];
}
