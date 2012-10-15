#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_exit(int NotDone)
#else
F_VOID_FUNC blacs_exit_(int *NotDone)
#endif
{
   void BI_UpdateBuffs(BLACBUFF *);
   BLACBUFF *BI_GetBuff(int);
   int BI_BuffIsFree(BLACBUFF *, int);
   BLACBUFF *bp;
   extern BLACBUFF *BI_ReadyB, *BI_ActiveQ, BI_AuxBuff;
   int i;
   extern int BI_MaxNCtxt, BI_Np;
   extern BLACSCONTEXT **BI_MyContxts;
/*
 * Destroy all contexts
 */
   for (i=0; i < BI_MaxNCtxt; i++) if (BI_MyContxts[i]) Cblacs_gridexit(i);
   free(BI_MyContxts);

   if (BI_ReadyB) free(BI_ReadyB);
   while (BI_ActiveQ != NULL)
   {
      bp = BI_ActiveQ;
      BI_BuffIsFree(bp, 1);  /* wait for async sends to complete */
      BI_ActiveQ = bp->next;
      free(bp);
   }
   free (BI_AuxBuff.Aops);

/*
 * Reset parameters to initial values
 */
   BI_MaxNCtxt = 0;
   BI_MyContxts = NULL;
   BI_Np = -1;
   if (!Mpval(NotDone))
   {
      MPI_Finalize();
   }
   BI_ReadyB = NULL;
}
