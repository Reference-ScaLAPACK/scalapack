#include "Bdef.h"
void BI_sMPI_amn2(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_svvamn2(int, char *, char *);
   BI_svvamn2(*N, inout, in);
}
