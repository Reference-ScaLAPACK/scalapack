#include "Bdef.h"
void BI_dMPI_amn2(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_dvvamn2(Int, char *, char *);
   BI_dvvamn2(*N, inout, in);
}
