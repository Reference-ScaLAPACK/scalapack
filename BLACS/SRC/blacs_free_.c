#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_freebuff(Int ConTxt, Int Wait)
#else
F_VOID_FUNC blacs_freebuff_(Int *ConTxt, Int *Wait)
#endif
{
   void BI_UpdateBuffs(BLACBUFF *);
   Int BI_BuffIsFree(BLACBUFF *, Int);
   extern BLACBUFF *BI_ReadyB, *BI_ActiveQ;

   if (Mpval(Wait))  /* wait for all buffers to be done */
   {
      while (BI_ActiveQ != NULL) BI_UpdateBuffs(NULL);
   }
   else BI_UpdateBuffs(NULL);

   if (BI_ReadyB)
   {
      free(BI_ReadyB);
      BI_ReadyB = NULL;
   }
}
