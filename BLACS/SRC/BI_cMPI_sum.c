#include "Bdef.h"
void BI_cMPI_sum(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_cvvsum(Int, char *, char *);
   BI_cvvsum(*N, inout, in);
}
