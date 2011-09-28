#include "Bdef.h"

void BI_zMPI_amn(void *in, void *inout, int *N, MPI_Datatype *dtype)
{
   void BI_zvvamn(int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_zvvamn(BI_AuxBuff.Len, inout, in);
}
