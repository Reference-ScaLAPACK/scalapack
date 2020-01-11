#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_barrier(Int ConTxt, char *scope)
#else
F_VOID_FUNC blacs_barrier_(Int *ConTxt, F_CHAR scope)
#endif
{
   char tscope;
   Int ierr;
   BLACSCONTEXT *ctxt;

   tscope = F2C_CharTrans(scope);
   tscope = Mlowcase(tscope);
   MGetConTxt(Mpval(ConTxt), ctxt);
   switch(tscope)
   {
   case 'r':
      MPI_Barrier(ctxt->rscp.comm);
      break;
   case 'c':
      MPI_Barrier(ctxt->cscp.comm);
      break;
   case 'a':
      MPI_Barrier(ctxt->ascp.comm);
      break;
   }
}
