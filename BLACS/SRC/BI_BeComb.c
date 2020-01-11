#include "Bdef.h"

/*
 *  The bidirectional exchange topology (BE) is specialized for dealing with
 *  case where all nodes participating in the operation need to
 *  receive the answer.  It works best when # of nodes is some even
 *  power of two.  This topology is based on an algorithm presented by
 *  Robert van de Geijn, et al.
 */
void BI_BeComb(BLACSCONTEXT *ctxt, BLACBUFF *bp, BLACBUFF *bp2,
               Int N, VVFUNPTR Xvvop)
/*
 *  -- V1.1ALPHA (test version) BLACS routine --
 *  University of Tennessee, October 1, 1995
 *  Written by Clint Whaley.
 *
 *  Purpose
 *  =======
 *  Perform a element-by-element combine on vectors.
 *  The answer will be left on all participating processes.  Since this method
 *  uses a hypercube communication pattern, the number of nodes participating
 *  in the operation must be a power of 2 for it to perform efficiently.
 *
 *  Arguments
 *  =========
 *  CTXT    (input) pointer to BLACSCONTEXT
 *          The BLACS context where operation is taking place.
 *
 *  BP      (input/output) pointer to BLACBUFF.
 *          BLACBUFF is a special data type used by the BLACS to control
 *          buffers and the asynchronous operations coming out of them.
 *          This BLACBUFF should have a buffer who's first N elements
 *          contain the data to be combined. Additional space may be
 *          required, depending upon what combine is being performed.
 *
 *  BP2     (workspace) pointer to BLACBUFF.
 *          This BLACBUFF is used to receive information for combining with
 *          this process's information.
 *
 *  N       (input) int
 *          The number of elements in the vector to be combined.
 *
 *  Xvvop   (input) pointer to typed operation function
 *          Points to a typed function which performs the required operation
 *          (e.g. summation) on the two N-element vectors.
 *
 * ------------------------------------------------------------------------
 */
{
   void BI_Ssend(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   void BI_Srecv(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   void BI_Rsend(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   void BI_Arecv(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   void BI_UpdateBuffs(BLACBUFF *);
   BLACBUFF *BI_GetBuff(Int);
   Int BI_BuffIsFree(BLACBUFF *, Int);

   Int Np, Iam, dest, msgid, Rmsgid, np2, bit, ierr;
   extern MPI_Status *BI_Stats;

   Np = ctxt->scp->Np;
   if (Np < 2) return;
   Iam = ctxt->scp->Iam;
   msgid = Mscopeid(ctxt);
   Rmsgid = Mscopeid(ctxt);

   for (np2=4; np2 < Np; np2 <<= 1);
   if (np2 > Np) np2 >>= 1;

   if (np2 != Np)
   {
      dest = (Iam ^ np2);
      if (Iam >= np2)  /* I'm node beyond power of 2 */
      {
	 BI_Arecv(ctxt, dest, Rmsgid, bp);
	 BI_Ssend(ctxt, dest, msgid, bp);
	 BI_BuffIsFree(bp, 1);
      }
      else if (Iam < (Np^np2))  /* need to fan in contents of */
      {                         /* non-power of 2 nodes */
         BI_Srecv(ctxt, dest, msgid, bp2);
	 Xvvop(N, bp->Buff, bp2->Buff);
      }
   }

   if (Iam < np2)
   {
      for (bit=1; (bit ^ np2); bit <<= 1)
      {
         dest = Iam ^ bit;
         ierr=MPI_Sendrecv(bp->Buff, bp->N, bp->dtype, dest, msgid, bp2->Buff,
                         bp2->N, bp2->dtype, dest, msgid, ctxt->scp->comm,
                         BI_Stats);
	 Xvvop(N, bp->Buff, bp2->Buff);
      }
/*
 *  For nodes that are not part of the hypercube proper, we must
 *  send data back.
 */
      if (Iam < (Np^np2)) BI_Rsend(ctxt, (Iam ^ np2), Rmsgid, bp);
   }  /* end if (nodes inside power of 2) */
}
