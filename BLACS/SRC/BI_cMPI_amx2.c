#include "Bdef.h"
void BI_cMPI_amx2(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_cvvamx2(int, char *, char *);
   BI_cvvamx2(*N, inout, in);
}
