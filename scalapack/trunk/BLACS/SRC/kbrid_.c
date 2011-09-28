#include "Bdef.h"

#if (INTFACE == C_CALL)
int Ckbrid(int ConTxt, char *scope, int rsrc, int csrc)
#else
F_INT_FUNC kbrid_(int *ConTxt, F_CHAR scope, int *rsrc, int *csrc)
#endif
{
   int msgid;
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
