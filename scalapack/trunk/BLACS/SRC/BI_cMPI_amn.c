#include "Bdef.h"

void BI_cMPI_amn(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_cvvamn(int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_cvvamn(BI_AuxBuff.Len, inout, in);
}
