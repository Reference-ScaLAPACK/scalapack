#include "Bdef.h"
#if (INTFACE == C_CALL)
MPI_Comm Cblacs2sys_handle(Int BlacsCtxt)
#else
Int blacs2sys_handle_(Int *BlacsCtxt)
#endif
{
#if (INTFACE == C_CALL)
   Int i[2];
   extern Int BI_MaxNSysCtxt;
   extern MPI_Comm *BI_SysContxts;

   if (BI_COMM_WORLD == NULL) Cblacs_pinfo(i, &i[1]);
   if ( (BlacsCtxt >= BI_MaxNSysCtxt) || (BlacsCtxt < 0) )
   {
      BI_BlacsErr(-1, __LINE__, __FILE__,
        "No system context corresponding to BLACS system context handle %d\n",
                  BlacsCtxt);
   }
   else if (BI_SysContxts[BlacsCtxt] == MPI_COMM_NULL)
   {
      BI_BlacsErr(-1, __LINE__, __FILE__,
        "No system context corresponding to BLACS system context handle %d\n",
                  BlacsCtxt);
   }
   return(BI_SysContxts[BlacsCtxt]);
#else
   return(*BlacsCtxt);
#endif
}
