#include "Bdef.h"
void BI_cMPI_amn2(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_cvvamn2(Int, char *, char *);
   BI_cvvamn2(*N, inout, in);
}
