#include "Bdef.h"

void BI_dMPI_amn(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_dvvamn(int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_dvvamn(BI_AuxBuff.Len, inout, in);
}
