#include "Bdef.h"
void BI_cMPI_amx2(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_cvvamx2(Int, char *, char *);
   BI_cvvamx2(*N, inout, in);
}
