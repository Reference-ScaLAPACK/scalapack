#include "Bdef.h"

void BI_dMPI_amn(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_dvvamn(Int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_dvvamn(BI_AuxBuff.Len, inout, in);
}
