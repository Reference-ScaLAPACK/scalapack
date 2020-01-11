#include "Bdef.h"

void BI_sMPI_amx(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_svvamx(Int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_svvamx(BI_AuxBuff.Len, inout, in);
}
