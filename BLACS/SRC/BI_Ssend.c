#include "Bdef.h"

void BI_Ssend(BLACSCONTEXT *ctxt, Int dest, Int msgid, BLACBUFF *bp)
{
   Int info;
   info=MPI_Send(bp->Buff, bp->N, bp->dtype, dest, msgid, ctxt->scp->comm);
}
