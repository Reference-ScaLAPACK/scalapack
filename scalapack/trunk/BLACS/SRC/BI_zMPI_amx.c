#include "Bdef.h"

void BI_zMPI_amx(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_zvvamx(int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_zvvamx(BI_AuxBuff.Len, inout, in);
}
