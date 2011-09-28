#include "Bdef.h"

void BI_BlacsErr(int ConTxt, int line, char *file, char *form, ...)
{
#ifdef __STDC__
   void BI_BlacsAbort(int ErrNo);
#else
   void BI_BlacsAbort();
#endif
   extern int BI_Iam;
   int myrow, mycol;
   va_list argptr;
   char cline[100];
   BLACSCONTEXT *ctxt;

   va_start(argptr, form);
   vsprintf(cline, form, argptr);
   va_end(argptr);

   if (ConTxt > -1)
   {
      MGetConTxt(ConTxt, ctxt);
      myrow = ctxt->cscp.Iam;
      mycol = ctxt->rscp.Iam;
   }
   else myrow = mycol = -1;

   fprintf(stderr,
"BLACS ERROR '%s'\nfrom {%d,%d}, pnum=%d, Contxt=%d, on line %d of file '%s'.\n\n",
           cline, myrow, mycol, BI_Iam, ConTxt, line, file);

   BI_BlacsAbort(1);
}
