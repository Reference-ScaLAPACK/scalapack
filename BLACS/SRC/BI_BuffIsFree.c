#include "Bdef.h"

Int BI_BuffIsFree(BLACBUFF *bp, Int Wait)
/*
 *  Check to see if buff is finished with async. operations.  If Wait != 0,
 *  wait for all async. operations to complete.
 */
{
   MpiInt i, info;
   extern MPI_Status *BI_Stats;
   extern Int BI_Np;


   if (!Wait)
   {
      info=MPI_Testall(bp->nAops, bp->Aops, &i, BI_Stats);
      if (!i)
      {
/*
 *       If we are doing our own Packing, need to check true length of receive
 */
#ifndef MpiBuffGood
/*
 *       If we have an outstanding receive, make sure that when it
 *       completes we correctly set bp->N, if required
 */
         if (bp->N < 0)
         {
            if (bp->Aops[-bp->N-1] == MPI_REQUEST_NULL)
            {
               info=MPI_Get_count(&BI_Stats[(-bp->N-1)*sizeof(MPI_Status)],MPI_PACKED, &i);
               if (i != MPI_UNDEFINED) bp->N = i;
               else BI_BlacsWarn(-1, __LINE__, __FILE__,
                                 "MPI_Get_count returned MPI_UNDEFINED.\n");
            }
         }
#endif
	 return(0);
      }
   }
   else
   {
      info=MPI_Waitall(bp->nAops, bp->Aops, BI_Stats);
   }

   bp->nAops = 0;
/*
 * If we are doing our own packing, need to check true length of receive
 */
#ifndef MpiBuffGood
/*
 * If we had an outstanding receive, make sure that we correctly set bp->N,
 * if required
 */
   if (bp->N < 0)
   {
      info=MPI_Get_count(&BI_Stats[(-bp->N-1)*sizeof(MPI_Status)],MPI_PACKED, &i);
      if (i != MPI_UNDEFINED) bp->N = i;
      else BI_BlacsWarn(-1, __LINE__, __FILE__,
                        "MPI_Get_count returned MPI_UNDEFINED.\n");
   }
#endif
   return(1);
}
