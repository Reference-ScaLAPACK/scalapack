#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_gridinit(Int *ConTxt, char *order, Int nprow, Int npcol)
#else
F_VOID_FUNC blacs_gridinit_(Int *ConTxt, F_CHAR order, Int *nprow, Int *npcol)
#endif
{
   void Cblacs_gridmap(Int *, Int *, Int, Int, Int);
   Int *tmpgrid, *iptr;
   Int i, j;

/*
 * Grid can be row- or column-major natural ordering when blacs_gridinit is
 * called.  Define a tmpgrid to reflect this, and call blacs_gridmap to
 * set it up
 */
   iptr = tmpgrid = (Int*) malloc( Mpval(nprow)*Mpval(npcol)*sizeof(*tmpgrid) );
   if (Mlowcase(F2C_CharTrans(order)) == 'c')
   {
      i = Mpval(npcol) * Mpval(nprow);
      for (j=0; j < i; j++) iptr[j] = j;
   }
   else
   {
      for (j=0; j < Mpval(npcol); j++)
      {
         for (i=0; i < Mpval(nprow); i++) iptr[i] = i * Mpval(npcol) + j;
         iptr += Mpval(nprow);
      }
   }
#if (INTFACE == C_CALL)
   Cblacs_gridmap(ConTxt, tmpgrid, nprow, nprow, npcol);
#else
   blacs_gridmap_(ConTxt, tmpgrid, nprow, nprow, npcol);
#endif
   free(tmpgrid);
}
