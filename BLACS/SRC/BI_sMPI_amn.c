#include "Bdef.h"

void BI_sMPI_amn(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_svvamn(Int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_svvamn(BI_AuxBuff.Len, inout, in);
}
