#include "Bdef.h"
void BI_iMPI_amn2(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_ivvamn2(Int, char *, char *);
   BI_ivvamn2(*N, inout, in);
}
