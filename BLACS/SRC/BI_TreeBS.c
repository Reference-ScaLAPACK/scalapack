#include "Bdef.h"

/*
 *  Tree_bs/br is a algorithm that does a broadcast send/recv such that the
 *  communication pattern is a tree with an arbitrary number of branches.
 *  The following two pairs of graphs give different ways of viewing the same
 *  algorithm.  The first pair shows the trees as they should be visualized
 *  when examining the algorithm.  The second pair are isomorphic graphs of
 *  of the first, which show the actual pattern of data movement.

 *  Note that a tree broadcast with NBRANCHES = 2 is isomorphic with a
 *  hypercube broadcast (however, it does not require the nodes be a
 *  power of two to work).
 *
 *    TREE BROADCAST, NBRANCHES = 2     *    TREE BROADCAST, NBRANCHES = 3
 *                                      *
 * i=4   &______________                *
 *       |              \               *
 * i=2   &______         &______        * i=3     &______________________
 *       |      \        |      \       *         |          \           \
 * i=1   &__     &__     &__     &__    * i=1     &______     &______     &__
 *       |  \    |  \    |  \    |  \   *         |  \   \    |  \   \    |  \
 *       0   1   2   3   4   5   6   7  *         0   1   2   3   4   5   6   7
 *
 *
 *          ISOMORPHIC GRAPHS OF ABOVE, SHOWN IN MORE FAMILIAR TERMS:
 *
 *                0                                           0
 *       _________|_________                       ___________|____________
 *      /         |         \                     /           |      |     \
 *     4          2          1                   6            3      2      1
 *    / \         |                              |           / \
 *   6   5        3                              7          4   5
 *   |
 *   7
 */
void BI_TreeBS(BLACSCONTEXT *ctxt, BLACBUFF *bp, SDRVPTR send, int nbranches)
{
   int Np, Iam, msgid, i, j;
   int destdist;	/* the distance of the destination node */

   Np = ctxt->scp->Np;
   if (Np < 2) return;
   Iam = ctxt->scp->Iam;
   msgid = Mscopeid(ctxt);

   for (i=nbranches; i < Np; i*=nbranches);
   for (i /= nbranches; i > 0; i /= nbranches)
   {
      j = 1;
      do
      {
	 destdist = i*j;
	 if (destdist < Np)
            send(ctxt, (destdist+Iam)%Np, msgid, bp);
      }
      while(++j < nbranches);
   }
} /* end BI_TreeBS */
