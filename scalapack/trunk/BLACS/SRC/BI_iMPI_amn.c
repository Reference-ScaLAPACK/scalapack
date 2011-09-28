#include "Bdef.h"

void BI_iMPI_amn(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_ivvamn(int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_ivvamn(BI_AuxBuff.Len, inout, in);
}
