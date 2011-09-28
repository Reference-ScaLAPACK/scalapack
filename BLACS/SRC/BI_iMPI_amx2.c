#include "Bdef.h"
void BI_iMPI_amx2(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_ivvamx2(int, char *, char *);
   BI_ivvamx2(*N, inout, in);
}
