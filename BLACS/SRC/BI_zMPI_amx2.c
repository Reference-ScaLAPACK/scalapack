#include "Bdef.h"
void BI_zMPI_amx2(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_zvvamx2(Int, char *, char *);
   BI_zvvamx2(*N, inout, in);
}
