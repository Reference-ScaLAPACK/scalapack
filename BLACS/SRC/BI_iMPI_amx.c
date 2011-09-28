#include "Bdef.h"

void BI_iMPI_amx(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_ivvamx(int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_ivvamx(BI_AuxBuff.Len, inout, in);
}
