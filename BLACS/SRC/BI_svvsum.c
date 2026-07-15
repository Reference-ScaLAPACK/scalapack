#include "Bdef.h"
void BI_svvsum(MpiInt N, char *vec1, char *vec2)
{
   float *v1=(float*)vec1, *v2=(float*)vec2;
   MpiInt k; // size_t might be more appropriate, as in BI_cvvsum.c BI_dvvsum.c BI_ivvsum.c BI_svvsum.c BI_zvvsum.c
   for (k=0; k < N; k++) v1[k] += v2[k];
}
