#include "Bdef.h"

void BI_iMPI_amn(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_ivvamn(Int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_ivvamn(BI_AuxBuff.Len, inout, in);
}
