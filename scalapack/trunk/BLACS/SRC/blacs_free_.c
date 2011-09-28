#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_freebuff(int ConTxt, int Wait)
#else
F_VOID_FUNC blacs_freebuff_(int *ConTxt, int *Wait)
#endif
{
   void BI_UpdateBuffs(BLACBUFF *);
   int BI_BuffIsFree(BLACBUFF *, int);
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
