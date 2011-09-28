#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_set(int ConTxt, int what, int *val)
#else
F_VOID_FUNC blacs_set_(int *ConTxt, int *what, int *val)
#endif
{
   BLACSCONTEXT *ctxt;

   switch( Mpval(what) )
   {
   case SGET_SYSCONTXT:
      BI_BlacsWarn(Mpval(ConTxt), __LINE__, __FILE__,
                   "Cannot set BLACS system context, can only BLACS_GET");
      break;
   case SGET_MSGIDS:
      BI_BlacsWarn(Mpval(ConTxt), __LINE__, __FILE__,
                   "No need to set message ID range due to MPI communicator.");
      break;
   case SGET_DEBUGLVL:
      BI_BlacsWarn(Mpval(ConTxt), __LINE__, __FILE__,
                   "Cannot set BLACS debug level; must recompile to change");
      break;
   case SGET_BLACSCONTXT:
      BI_BlacsWarn(Mpval(ConTxt), __LINE__, __FILE__,
                   "Cannot set BLACS context, can only BLACS_GET");
      break;
   case SGET_NR_BS:
      MGetConTxt(Mpval(ConTxt), ctxt);
      if (*val) ctxt->Nr_bs = *val;
      else BI_BlacsWarn(Mpval(ConTxt), __LINE__, __FILE__,
                        "BSBR nrings cannot be set to zero");
      break;
   case SGET_NB_BS:
      MGetConTxt(Mpval(ConTxt), ctxt);
      if (*val > 0) ctxt->Nb_bs = *val + 1;
      else BI_BlacsWarn(Mpval(ConTxt), __LINE__, __FILE__,
                       "Illegal BSBR nbranches (%d); must be strictly positive",
                        *val);
      break;
   case SGET_NR_CO:
      MGetConTxt(Mpval(ConTxt), ctxt);
      if (*val) ctxt->Nr_co = *val;
      else BI_BlacsWarn(Mpval(ConTxt), __LINE__, __FILE__,
                        "COMB nrings cannot be set to zero");
      break;
   case SGET_NB_CO:
      MGetConTxt(Mpval(ConTxt), ctxt);
      if (*val > 0) ctxt->Nb_co = *val + 1;
      else BI_BlacsWarn(Mpval(ConTxt), __LINE__, __FILE__,
                       "Illegal COMB nbranches (%d); must be strictly positive",
                        *val);
      break;
   case SGET_TOPSREPEAT:
      MGetConTxt(Mpval(ConTxt), ctxt);
      ctxt->TopsRepeat = *val;
      break;
   case SGET_TOPSCOHRNT:
      MGetConTxt(Mpval(ConTxt), ctxt);
      ctxt->TopsCohrnt = *val;
      break;
   default:
      BI_BlacsWarn(Mpval(ConTxt), __LINE__, __FILE__, "Unknown WHAT (%d)",
                   Mpval(what));
   }
}
