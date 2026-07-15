#include "Bdef.h"
void BI_dvvsum(MpiInt N, char *vec1, char *vec2)
{
   double *v1=(double*)vec1, *v2=(double*)vec2;
   MpiInt k; // size_t might be more appropriate, as in BI_cvvsum.c BI_dvvsum.c BI_ivvsum.c BI_svvsum.c BI_zvvsum.c
   for (k=0; k < N; k++) v1[k] += v2[k];
}
