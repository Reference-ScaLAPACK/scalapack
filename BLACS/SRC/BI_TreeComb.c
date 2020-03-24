#include "Bdef.h"

/*
 *  This topology supports trees with arbitrary numbers of branches at
 *  each step.  The following pictures show the tree that should be visualized
 *  when examining the algorithm.
 *
 *    TREE GLOBAL OP, NBRANCHES = 2     *    TREE GLOBAL OP, NBRANCHES = 3
 *                                      *
 * i=4   &______________                *
 *       |              \               *
 * i=2   &______         &______        * i=3     &______________________
 *       |      \        |      \       *         |          \           \
 * i=1   &__     &__     &__     &__    * i=1     &______     &______     &__
 *       |  \    |  \    |  \    |  \   *         |  \   \    |  \   \    |  \
 *       0   1   2   3   4   5   6   7  *         0   1   2   3   4   5   6   7
 */

void BI_TreeComb(BLACSCONTEXT *ctxt, BLACBUFF *bp, BLACBUFF *bp2,
                 Int N, VVFUNPTR Xvvop, Int dest, Int nbranches)
/*
 *  -- V1.1ALPHA (test version) BLACS routine --
 *  University of Tennessee, October 1, 1995
 *  Written by Clint Whaley.
 *
 *  Purpose
 *  =======
 *  Perform a element-by-element combine on vectors.
 *  If rdest1 = -1, the answer will be left on all participating processes.
 *  Otherwise, only the process at grid coordinates {rdest1, cdest1} will
 *  have the final answer.  Other Processes will have intermediate (useless)
 *  values.
 *
 *  Arguments
 *  =========
 *  CTXT      (input) pointer to BLACSCONTEXT
 *            The BLACS context where operation is taking place.
 *
 *  BP        (input/output) pointer to BLACBUFF.
 *            BLACBUFF is a special data type used by the BLACS to control
 *            buffers and the asynchronous operations coming out of them.
 *            This BLACBUFF should have a buffer who's first N elements
 *            contain the data to be combined. Additional space may be
 *            required, depending upon what combine is being performed.
 *
 *  BP2       (workspace) pointer to BLACBUFF.
 *            This BLACBUFF is used to receive information for combining with
 *            this process's information.
 *
 *  DEST      (input) Int
 *            Node to receive answer.  If DEST == -1, all nodes in receive
 *            the answer.
 *
 *  N         (input) Int
 *            The number of elements in the vector.  N >= 0.
 *
 *  Xvvop     (input) pointer to typed operation function
 *            Points to a typed function which performs the required operation
 *            (e.g. summation) on the two N-element vectors.
 *
 *  NBRANCHES (input) Int
 *            Indicates the degree of the tree to use (see picture above).
 *
 * ------------------------------------------------------------------------
 */
{
   void BI_UpdateBuffs(BLACBUFF *);
   BLACBUFF *BI_GetBuff(Int);
   Int BI_BuffIsFree(BLACBUFF *, Int);
   void BI_Ssend(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   void BI_Srecv(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   void BI_Rsend(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   void BI_Arecv(BLACSCONTEXT *, Int, Int, BLACBUFF *);

   Int Np, Iam, msgid, Rmsgid, i, j;
   Int nrcvs=0;	  /* Number of ReCeiVeS to do */
   Int REBS;	  /* should info be RE-BroadcaSt? */
   Int rightedge; /* right-most receiving node */
   Int mydist;    /* my distance from destination node */
   Int dist;
   Int src;       /* Used if we must force repeatability */

   Np = ctxt->scp->Np;
   if (Np < 2) return;
   Iam = ctxt->scp->Iam;
   msgid = Mscopeid(ctxt);
   Rmsgid = Mscopeid(ctxt);
   if (REBS = (dest == -1)) dest = 0;

   mydist = (Np + Iam - dest) % Np;
   if (REBS)
   {
      dist = mydist;
      if (mydist != 0) BI_Arecv(ctxt, BANYNODE, Rmsgid, bp);
   }

   if (nbranches == FULLCON) nbranches = Np;
   rightedge = Np - 1 - (Np-1)%nbranches;

   for (i=1; (i < Np); i *= nbranches)
   {
      if (mydist%nbranches)	/* nodes that send to other nodes */
      {
	 BI_Ssend(ctxt, (dest + (mydist-mydist%nbranches)*i)%Np, msgid, bp);
	 break;		/* I'm done */
      }
      else
      {
         if (mydist != rightedge) nrcvs = nbranches - 1;
         else nrcvs = (Np + i - 1) / i - rightedge - 1;
         mydist /= nbranches;
         rightedge /= nbranches;
         rightedge -= (rightedge % nbranches);

         if (!ctxt->TopsRepeat)
         {
            for (j=nrcvs; j; j--)
            {
               BI_Srecv(ctxt, BANYNODE, msgid, bp2);
	       Xvvop(N, bp->Buff, bp2->Buff);
            }
         }
         else
         {
            src = (Iam + i) % Np;
            for (j=nrcvs; j; j--)
            {
               BI_Srecv(ctxt, src, msgid, bp2);
	       Xvvop(N, bp->Buff, bp2->Buff);
               src = (src + i) % Np;
            }
         }
      }
   }

/*
 * Broadcast answer to everyone if RDEST == -1
 */
   if (REBS)
   {
      mydist = dist;
      for (i=2; i < Np; i <<= 1);
      if (mydist > 0) BI_BuffIsFree(bp, 1);

      while (i > 1)
      {
         if ( !(mydist%i) )
         {
            i >>= 1;
            dist = mydist + i;
	    if (dist < Np) BI_Rsend(ctxt, dist, Rmsgid, bp);
         }
         else i >>= 1;
      }
   }
} /* end BI_TreeComb */
