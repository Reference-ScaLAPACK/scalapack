#include "Bdef.h"
void BI_sMPI_amx2(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_svvamx2(int, char *, char *);
   BI_svvamx2(*N, inout, in);
}
