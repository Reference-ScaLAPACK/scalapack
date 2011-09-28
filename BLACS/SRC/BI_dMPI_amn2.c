#include "Bdef.h"
void BI_dMPI_amn2(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_dvvamn2(int, char *, char *);
   BI_dvvamn2(*N, inout, in);
}
