#include "Bdef.h"

void BI_UpdateBuffs(BLACBUFF *Newbp)
{
   int BI_BuffIsFree(BLACBUFF *, int);
   BLACBUFF *bp, *bp2;
   extern BLACBUFF *BI_ReadyB, *BI_ActiveQ;

   if (Newbp)
   {
      if (BI_ActiveQ == NULL) BI_ActiveQ = Newbp->prev = Newbp;
      else
      {
	 BI_ActiveQ->prev->next = Newbp;
	 Newbp->prev = BI_ActiveQ->prev;
	 BI_ActiveQ->prev = Newbp;
      }
      Newbp->next = NULL;
      if (Newbp == BI_ReadyB) BI_ReadyB = NULL;
   }
/*
 * See if any active buffers are ready for reuse.
 */
   for (bp=BI_ActiveQ; bp != NULL; bp = bp2)
   {
      bp2 = bp->next;
      if ( BI_BuffIsFree(bp, 0) )  /* if all of buff's Aops are done */
      {
/*
 *       Remove bp from BI_ActiveQ -- update pointers
 */
	 if (bp->next) bp->next->prev = bp->prev;
	 else BI_ActiveQ->prev = bp->prev;
	 if (bp != BI_ActiveQ) bp->prev->next = bp->next;
	 else BI_ActiveQ = BI_ActiveQ->next;

/*
 *       If no ready buffer, inactive buff becomes ready
 */
	 if (BI_ReadyB == NULL) BI_ReadyB = bp;
/*
 *       If inactive buff bigger than present ready buff, release ready,
 *       and inactive buff becomes ready
 */
	 else if (BI_ReadyB->Len < bp->Len)
	 {
	    free(BI_ReadyB);
	    BI_ReadyB = bp;
	 }
/*
 *       If ready buffer exists and is bigger than inactive buff,
 *       free inactive buff
 */
	 else free(bp);
      }
   }
}  /* end BI_UpdateBuffs */
