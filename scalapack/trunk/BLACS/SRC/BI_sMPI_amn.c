#include "Bdef.h"

void BI_sMPI_amn(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_svvamn(int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_svvamn(BI_AuxBuff.Len, inout, in);
}
