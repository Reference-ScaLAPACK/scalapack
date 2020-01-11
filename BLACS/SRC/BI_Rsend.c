#include "Bdef.h"


void BI_Rsend(BLACSCONTEXT *ctxt, Int dest, Int msgid, BLACBUFF *bp)
{
   Int info;

   info=MPI_Rsend(bp->Buff, bp->N, bp->dtype, dest, msgid, ctxt->scp->comm);
}
