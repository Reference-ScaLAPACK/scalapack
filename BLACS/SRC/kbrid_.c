#include "Bdef.h"

#if (INTFACE == C_CALL)
Int Ckbrid(Int ConTxt, char *scope, Int rsrc, Int csrc)
#else
F_INT_FUNC kbrid_(Int *ConTxt, F_CHAR scope, Int *rsrc, Int *csrc)
#endif
{
   Int msgid;
   char tmpscope;
   BLACSCONTEXT *ctxt;

   MGetConTxt(Mpval(ConTxt), ctxt);
   tmpscope = Mlowcase(F2C_CharTrans(scope));
   switch(tmpscope)
   {
   case 'c' :
      ctxt->scp = &ctxt->cscp;
      break;
   case 'r' :
      ctxt->scp = &ctxt->cscp;
      break;
   case 'a' :
      ctxt->scp = &ctxt->cscp;
      break;
   }
   msgid = Mscopeid(ctxt);
   return (msgid);
}
