#include "Bdef.h"

void BI_Asend(BLACSCONTEXT *ctxt, int dest, int msgid, BLACBUFF *bp)
{
   int i, info, errclass;

   info=MPI_Isend(bp->Buff, bp->N, bp->dtype, dest, msgid, ctxt->scp->comm,
                &bp->Aops[bp->nAops]);
   while(info != MPI_SUCCESS)
   {
      i=MPI_Error_class(info, &errclass);
      if ( (errclass != MPI_ERR_UNKNOWN) && (errclass != MPI_ERR_OTHER) &&
           (errclass != MPI_ERR_INTERN) )
      {
	  Mmpierror(info, "MPI_Isend", ctxt, __LINE__, __FILE__);
	  BI_BlacsErr(BI_ContxtNum(ctxt), __LINE__, __FILE__,
		      "MPI error %d on call to MPI_Isend", info);
      }
#if (BlacsDebugLvl > 0)
      else BI_BlacsWarn(BI_ContxtNum(ctxt), __LINE__, __FILE__,
"MPI error %d assumed to mean out of non-blocking resources on call to MPI_Isend",
                        info);
#endif
      info=MPI_Isend(bp->Buff, bp->N, bp->dtype, dest, msgid, ctxt->scp->comm,
                   &bp->Aops[bp->nAops]);
   }
   bp->nAops++;
}
