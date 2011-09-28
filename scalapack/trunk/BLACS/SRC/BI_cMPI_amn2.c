#include "Bdef.h"
void BI_cMPI_amn2(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_cvvamn2(int, char *, char *);
   BI_cvvamn2(*N, inout, in);
}
