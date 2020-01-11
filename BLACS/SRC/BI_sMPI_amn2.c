#include "Bdef.h"
void BI_sMPI_amn2(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_svvamn2(Int, char *, char *);
   BI_svvamn2(*N, inout, in);
}
