#include "Bdef.h"
#if (INTFACE == C_CALL)
void Cblacs_get(Int ConTxt, Int what, Int *val)
#else
F_VOID_FUNC blacs_get_(Int *ConTxt, Int *what, Int *val)
#endif
{
   Int Csys2blacs_handle(MPI_Comm);
   Int ierr, *iptr;
   MpiInt flag;
   Int comm;
   BLACSCONTEXT *ctxt;

   switch( Mpval(what) )
   {
   case SGET_SYSCONTXT:
      if (BI_COMM_WORLD == NULL) Cblacs_pinfo(val, &ierr);
#if (INTFACE == C_CALL)
      *val = Csys2blacs_handle(MPI_COMM_WORLD);
#else
      *val = *BI_COMM_WORLD;
#endif
      break;
   case SGET_MSGIDS:
      if (BI_COMM_WORLD == NULL) Cblacs_pinfo(val, &val[1]);
      iptr = &val[1];
      ierr=MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, (BVOID **) &iptr,&flag);
      val[0] = 0;
      val[1] = *iptr;
      break;
   case SGET_DEBUGLVL:
      *val = BlacsDebugLvl;
      break;
   case SGET_BLACSCONTXT:
      MGetConTxt(Mpval(ConTxt), ctxt);
#if (INTFACE == C_CALL)
      *val = Csys2blacs_handle(ctxt->pscp.comm);
#else  /* if user called the fortran interface to the BLACS */
      *val = MPI_Comm_c2f(ctxt->pscp.comm);
#endif
      break;
   case SGET_NR_BS:
      MGetConTxt(Mpval(ConTxt), ctxt);
      *val = ctxt->Nr_bs;
      break;
   case SGET_NB_BS:
      MGetConTxt(Mpval(ConTxt), ctxt);
      *val = ctxt->Nb_bs - 1;
      break;
   case SGET_NR_CO:
      MGetConTxt(Mpval(ConTxt), ctxt);
      *val = ctxt->Nr_co;
      break;
   case SGET_NB_CO:
      MGetConTxt(Mpval(ConTxt), ctxt);
      *val = ctxt->Nb_co - 1;
      break;
   case SGET_TOPSREPEAT:
      MGetConTxt(Mpval(ConTxt), ctxt);
      *val = ctxt->TopsRepeat;
      break;
   case SGET_TOPSCOHRNT:
      MGetConTxt(Mpval(ConTxt), ctxt);
      *val = ctxt->TopsCohrnt;
      break;
   default:
      BI_BlacsWarn(Mpval(ConTxt), __LINE__, __FILE__, "Unknown WHAT (%d)",
                Mpval(what));
   }
}
