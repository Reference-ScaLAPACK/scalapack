#include "Bdef.h"

void BI_dMPI_amx(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_dvvamx(Int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_dvvamx(BI_AuxBuff.Len, inout, in);
}
