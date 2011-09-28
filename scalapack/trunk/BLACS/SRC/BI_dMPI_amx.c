#include "Bdef.h"

void BI_dMPI_amx(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_dvvamx(int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_dvvamx(BI_AuxBuff.Len, inout, in);
}
