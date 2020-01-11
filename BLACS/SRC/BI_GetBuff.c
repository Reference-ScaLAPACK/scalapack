#include "Bdef.h"
/***************************************************************************
 *  The mpi implements globally blocking sends.  I.e., a send blocks until *
 *  the dest. node issues a recv.  The BLACS assume locally-blocking sends.*
 *  Therefore, the BLACS must fake locally-blocking sends.  To do this     *
 *  requires an indeterminate number of buffers and the use of             *
 *  non-blocking sends.  However, it is very important that even though I  *
 *  provide a dynamic number of buffers, that getting these buffers does   *
 *  not take too long in the critical part of a send operation.            *
 *  Therefore, the buffer management is broken into two routines.          *
 *									   *
 *  Inside the BLACS there are two states a buffer may be in.  If the buff *
 *  is currently being used (for instance, an asynchronous send is coming  *
 *  from it), it is classified as an ACTIVE buffer, and is on the active   *
 *  buffer queue.  Otherwise, a buffer is READY: it is not being used      *
 *  and is available for the next buffer operation.                        *
 *  In order to avoid buffer proliferation, only one ready buffer is kept, *
 *  and as active buffers become inactive they either become the ready     *
 *  buffer, or are freed.						   *
 *  									   *
 *  The first routine, BI_GetBuff, checks if the ready buffer is big enough   *
 *  to fulfill the buffer request.  If not, the present ready buffer is    *
 *  is freed, and a new buffer of the required length is allocated.  If    *
 *  the buffer is of sufficent size already, no action is taken.           *
 *  This routine is purposely very short, as it is called at the beginning *
 *  of each broadcast/send operation.  All participating nodes             *
 *  are waiting on the source node, so this routine must be very cheap.	   *
 *									   *
 *  The second routine, BI_UpdateBuffs, moves the ready buffer to the active  *
 *  buffer queue (if needed).  It also checks the entire active buffer     *
 *  queue to see if any have finished their operations.  If so, they are   *
 *  are either moved to the ready buff, or freed.  This routine is called  *
 *  AFTER the send/broadcast has been started, and thus I am free to make  *
 *  it a little more complex.						   *
 ***************************************************************************/

BLACBUFF *BI_GetBuff(Int length)
{
   void BI_EmergencyBuff(Int length);

   char *cptr;
   Int i, j;
   extern Int BI_Np;
   extern BLACBUFF *BI_ReadyB;

/*
 * If ready buffer already exists, and is big enough, return it.  Otherwise,
 * free the buffer (if it exists) and get one of correct size
 */
   if (BI_ReadyB)
   {
      if (BI_ReadyB->Len >= length) return(BI_ReadyB);
      else free(BI_ReadyB);
   }
/*
 * Make sure all buffers aligned correctly
 */
   j = sizeof(BLACBUFF);
   if (j % sizeof(MPI_Request))
      j += sizeof(MPI_Request) - j % sizeof(MPI_Request);
   i = j + BI_Np*sizeof(MPI_Request);
   if (i % BUFFALIGN) i += BUFFALIGN - i % BUFFALIGN;
   cptr = malloc(i + length);
   BI_ReadyB = (BLACBUFF *) cptr;

   if (BI_ReadyB != NULL)
   {
      BI_ReadyB->nAops = 0;
      BI_ReadyB->Aops = (MPI_Request *) &cptr[j];
      BI_ReadyB->Buff = &cptr[i];
      BI_ReadyB->Len = length;
   }
   else BI_EmergencyBuff(length);

   return(BI_ReadyB);
}
