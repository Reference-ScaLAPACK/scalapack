#include "Bdef.h"

void BI_zMPI_amx(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_zvvamx(Int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_zvvamx(BI_AuxBuff.Len, inout, in);
}
