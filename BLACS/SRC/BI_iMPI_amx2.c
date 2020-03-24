#include "Bdef.h"
void BI_iMPI_amx2(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_ivvamx2(Int, char *, char *);
   BI_ivvamx2(*N, inout, in);
}
