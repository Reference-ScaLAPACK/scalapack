#include "Bdef.h"
void BI_zMPI_amn2(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_zvvamn2(int, char *, char *);
   BI_zvvamn2(*N, inout, in);
}
