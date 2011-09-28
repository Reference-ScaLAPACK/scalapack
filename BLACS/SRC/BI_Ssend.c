#include "Bdef.h"

void BI_Ssend(BLACSCONTEXT *ctxt, int dest, int msgid, BLACBUFF *bp)
{
   int info;
   info=MPI_Send(bp->Buff, bp->N, bp->dtype, dest, msgid, ctxt->scp->comm);
}
