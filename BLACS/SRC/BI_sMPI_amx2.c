#include "Bdef.h"
void BI_sMPI_amx2(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_svvamx2(Int, char *, char *);
   BI_svvamx2(*N, inout, in);
}
