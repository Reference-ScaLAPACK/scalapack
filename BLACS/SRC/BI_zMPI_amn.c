#include "Bdef.h"

void BI_zMPI_amn(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_zvvamn(Int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_zvvamn(BI_AuxBuff.Len, inout, in);
}
