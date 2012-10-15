#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_abort(int ConTxt, int ErrNo)
#else
F_VOID_FUNC blacs_abort_(int *ConTxt, int *ErrNo)
#endif
{
   void BI_BlacsAbort(int ErrNo);
   int nprow, npcol, myrow, mycol;
   extern int BI_Iam;

   Cblacs_gridinfo(Mpval(ConTxt), &nprow, &npcol, &myrow, &mycol);
   fprintf(stderr,
"{%d,%d}, pnum=%d, Contxt=%d, killed other procs, exiting with error #%d.\n\n",
           myrow, mycol, BI_Iam, Mpval(ConTxt), Mpval(ErrNo));

   BI_BlacsAbort(Mpval(ErrNo));
}
