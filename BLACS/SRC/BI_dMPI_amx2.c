#include "Bdef.h"
void BI_dMPI_amx2(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_dvvamx2(int, char *, char *);
   BI_dvvamx2(*N, inout, in);
}
