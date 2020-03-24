#include "Bdef.h"

/***************************************************************************
 *  If there is insufficient space to allocate a needed buffer, this       *
 *  routine is called.  It moniters active buffers for the time defined by *
 *  the user-changeable macro value BUFWAIT.  If in that time no active    *
 *  buffer becomes inactive, a hang is assumed, and the grid is killed.    *
 ***************************************************************************/
void BI_EmergencyBuff(Int length)
{
   void BI_UpdateBuffs(BLACBUFF *);

   char *cptr;
   Int i, j;
   double Mwalltime(void);
   double t1;
   extern Int BI_Np;
   extern BLACBUFF *BI_ReadyB, *BI_ActiveQ;

   j = sizeof(BLACBUFF);
   if (j % sizeof(MPI_Request))
      j += sizeof(MPI_Request) - j % sizeof(MPI_Request);
   i = j + BI_Np*sizeof(MPI_Request);
   if (i % BUFFALIGN) i += BUFFALIGN - i % BUFFALIGN;
   t1 =  Mwalltime();
   while ( (BI_ActiveQ) && (Mwalltime() - t1 < BUFWAIT) && !(BI_ReadyB) )
   {
      BI_UpdateBuffs(NULL);
      if (BI_ReadyB)
      {
         if (BI_ReadyB->Len < length)
         {
	    free(BI_ReadyB);
            cptr = malloc(length + i);
            BI_ReadyB = (BLACBUFF *) cptr;
            if (BI_ReadyB)
            {
               BI_ReadyB->nAops = 0;
               BI_ReadyB->Aops = (MPI_Request *) &cptr[j];
               BI_ReadyB->Buff = &cptr[i];
               BI_ReadyB->Len = length;
            }
         }
      }
   }
   if (BI_ReadyB == NULL)
   {
      BI_BlacsErr(-1, __LINE__, __FILE__, "BLACS out of buffer space");
   }
}
