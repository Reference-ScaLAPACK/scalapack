#include "Bdef.h"


void BI_ArgCheck(Int ConTxt, Int RoutType, char *routine, char scope,
                 char uplo, char diag, Int m, Int n, Int lda, Int nprocs,
                 Int *prows, Int *pcols)
{
#if (BlacsDebugLvl > 0)
   char *srcdest;
   Int i=1, prow, pcol, Ng, nprow, npcol, myrow, mycol;
   BLACSCONTEXT *ctxt;

   MGetConTxt(ConTxt, ctxt);
   Mgridinfo(ctxt, Ng, nprow, npcol, myrow, mycol);

   if ( (scope != 'r') && (scope != 'c') && (scope != 'a') )
      BI_BlacsErr(ConTxt, -1, routine, "Unknown scope, scope=%c", scope);
   if ( (uplo != 'u') && (uplo != 'l') )
   {
      if (RoutType != RT_COMB)
         BI_BlacsWarn(ConTxt, -1, routine,
                      "UPLO=%c, will be assumed to mean LOWER", uplo);
      else i = 0;  /* combine aux, for rect. matrix */
   }
   if ( (diag != 'u') && (diag != 'n') )
   {
      if (i) BI_BlacsWarn(ConTxt, -1, routine,
                          "DIAG=%c, will be assumed to mean NON-UNIT", diag);
   }
   if (m * n != 0)
   {
      if (m < 0)
         BI_BlacsErr(ConTxt, -1, routine, "Illegal number of rows, M=%d", m);
      if (n < 0)
         BI_BlacsErr(ConTxt, -1, routine, "Illegal number of columns, N=%d", n);
      if (lda < m)
         BI_BlacsWarn(ConTxt, -1, routine,
                      "Illegal LDA, LDA=%d, M=%d; LDA assumed to be %d",
                      lda, m, m);
   }

   if ( (RoutType == RT_RV) || (RoutType == RT_BR) ) srcdest = "SRC";
   else srcdest = "DEST";

   if (RoutType == RT_SD)
   {
      if ( (nprocs > Ng) || (nprocs < 0) )
         BI_BlacsErr(ConTxt, -1, routine,
                     "Trying to send to %d procs, but only %d in grid",
                     nprocs, Ng);
   }

   for (i=0; i < nprocs; i++)
   {
      prow = prows[i];
      pcol = pcols[i];

      if ( (prow < 0) || (prow >= nprow) )
      {
         if ( !((RoutType == RT_COMB) && (prow == -1)) )
            BI_BlacsErr(ConTxt, -1, routine,
                        "R%s out of range; R%s=%d, NPROW=%d",
                        srcdest, srcdest, prow, nprow);
      }
      if ( (pcol < 0) || (pcol >= npcol) )
      {
         if ( !((RoutType == RT_COMB) && (prow == -1)) )
            BI_BlacsErr(ConTxt, -1, routine,
                        "C%s out of range; C%s=%d, NPCOL=%d",
                        srcdest, srcdest, pcol, npcol);
      }
      if (RoutType == RT_SD)  /* point to point send */
      {
         if ( (prow == myrow) && (pcol == mycol) )
            BI_BlacsWarn(ConTxt, -1, routine, "Node sending message to itself");
      }
      else if (RoutType == RT_RV)  /* point to point send */
      {
         if ( (prow == myrow) && (pcol == mycol) )
            BI_BlacsWarn(ConTxt, -1, routine,
                         "Node recving message from itself");
      }
      else if (RoutType == RT_BR) /* broadcast/recv */
      {
         if ( (prow == myrow) && (pcol == mycol) )
            BI_BlacsErr(ConTxt, -1, routine,
                        "Node tries to recv its own broadcast");

         if (scope == 'r')
         {
            if (myrow != prow)
               BI_BlacsWarn(ConTxt, -1, routine,
                            "Row broadcast: MYROW=%d, but RSRC=%d",
                            myrow, prow);
         }
         else if (scope == 'c')
         {
            if (mycol != pcol)
            {
               BI_BlacsErr(ConTxt, -1, routine,
                           "Column broadcast: MYCOL=%d, but CSRC=%d",
                           mycol, pcol);
            }
         }
      }
   }
#endif
}
