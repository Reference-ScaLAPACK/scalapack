#include "Bdef.h"

void BI_Srecv(BLACSCONTEXT *ctxt, int src, int msgid, BLACBUFF *bp)
{
   int i, info;
   extern MPI_Status *BI_Stats;

   info=MPI_Recv(bp->Buff, bp->N, bp->dtype, src, msgid, ctxt->scp->comm,BI_Stats);
/*
 * If we are doing our own buffering, need to determine the true length of
 * the message just received
 */
#ifndef MpiBuffGood
   if (bp->dtype == MPI_PACKED)
   {
      info=MPI_Get_count(BI_Stats, MPI_PACKED, &i);
      if (i != MPI_UNDEFINED) bp->N = i;
      else BI_BlacsWarn(BI_ContxtNum(ctxt), __LINE__, __FILE__,
                        "MPI_Get_count returned MPI_UNDEFINED.\n");
   }
#endif
}
