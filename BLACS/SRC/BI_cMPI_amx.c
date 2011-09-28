#include "Bdef.h"

void BI_cMPI_amx(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_cvvamx(int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_cvvamx(BI_AuxBuff.Len, inout, in);
}
