#include "Bdef.h"
void BI_zMPI_sum(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_zvvsum(int, char *, char *);
   BI_zvvsum(*N, inout, in);
}
