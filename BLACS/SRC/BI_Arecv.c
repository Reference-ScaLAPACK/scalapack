#include "Bdef.h"

void BI_Arecv(BLACSCONTEXT *ctxt, Int src, Int msgid, BLACBUFF *bp)
{
   Int i, info;
   MpiInt errclass;

   info=MPI_Irecv(bp->Buff, bp->N, bp->dtype, src, msgid, ctxt->scp->comm,
                &bp->Aops[bp->nAops]);
   while(info != MPI_SUCCESS)
   {
      i=MPI_Error_class(info, &errclass);
      if ( (errclass != MPI_ERR_UNKNOWN) && (errclass != MPI_ERR_OTHER) &&
           (errclass != MPI_ERR_INTERN) )
      {
	  Mmpierror(info, "MPI_Irecv", ctxt, __LINE__, __FILE__);
	  BI_BlacsErr(BI_ContxtNum(ctxt), __LINE__, __FILE__,
		      "MPI error %d on call to MPI_Irecv", info);
      }
#if (BlacsDebugLvl > 0)
      else BI_BlacsWarn(BI_ContxtNum(ctxt), __LINE__, __FILE__,
"MPI error %d assumed to mean out of non-blocking resources on call to MPI_Irecv",
                     info);
#endif
      info=MPI_Irecv(bp->Buff, bp->N, bp->dtype, src, msgid, ctxt->scp->comm,
                   &bp->Aops[bp->nAops]);
   }
   bp->nAops++;
/*
 * Signal if we need to use status to figure out true length of received message
 * We only need do this if we are doing our own buffering
 */
#ifndef MpiBuffGood
   if (bp->dtype == MPI_PACKED) bp->N = -bp->nAops;
#endif
}
