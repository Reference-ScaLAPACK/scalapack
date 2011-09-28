#include "Bdef.h"

void BI_sMPI_amx(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_svvamx(int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_svvamx(BI_AuxBuff.Len, inout, in);
}
