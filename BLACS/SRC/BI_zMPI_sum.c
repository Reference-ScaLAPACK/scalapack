#include "Bdef.h"
void BI_zMPI_sum(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_zvvsum(Int, char *, char *);
   BI_zvvsum(*N, inout, in);
}
