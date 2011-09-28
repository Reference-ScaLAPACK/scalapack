#include "Bdef.h"
void BI_iMPI_amn2(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_ivvamn2(int, char *, char *);
   BI_ivvamn2(*N, inout, in);
}
