#include "Bdef.h"

int BI_ContxtNum(BLACSCONTEXT *ctxt)
/*
 *  Returns the integer ID of ctxt
 */
{
   int i;
   extern int BI_MaxNCtxt;
   extern BLACSCONTEXT **BI_MyContxts;

   if (ctxt == NULL) return(-1);
   for (i=0; i < BI_MaxNCtxt; i++) if (BI_MyContxts[i] == ctxt) break;
   if (i == BI_MaxNCtxt)
      BI_BlacsErr(-1, -1, "BLACS INTERNAL ROUTINE", "illegal context");
   return(i);
}
