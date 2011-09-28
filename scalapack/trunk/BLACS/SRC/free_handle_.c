#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cfree_blacs_system_handle(int ISysCtxt)
#else
void free_blacs_system_handle_(int *ISysCxt)
#endif
{
#if (INTFACE == C_CALL)
   int i, j, DEF_WORLD;
   MPI_Comm *tSysCtxt;
   extern int BI_MaxNSysCtxt;
   extern MPI_Comm *BI_SysContxts;


   if ( (ISysCtxt < BI_MaxNSysCtxt) && (ISysCtxt > 0) )
   {
      if (BI_SysContxts[ISysCtxt] != MPI_COMM_NULL)
         BI_SysContxts[ISysCtxt] = MPI_COMM_NULL;
      else BI_BlacsWarn(-1, __LINE__, __FILE__,
          "Trying to free non-existent system context handle %d", ISysCtxt);
   }
   else if (ISysCtxt == 0) return;  /* never free MPI_COMM_WORLD */
   else BI_BlacsWarn(-1, __LINE__, __FILE__,
        "Trying to free non-existent system context handle %d", ISysCtxt);

/*
 * See if we have freed enough space to decrease the size of our table
 */
   for (i=j=0; i < BI_MaxNSysCtxt; i++)
      if (BI_SysContxts[i] == MPI_COMM_NULL) j++;
/*
 * If needed, get a smaller system context array
 */
   if (j > 2*MAXNSYSCTXT)
   {
      j = BI_MaxNSysCtxt - MAXNSYSCTXT;
      tSysCtxt = (MPI_Comm *) malloc(j * sizeof(MPI_Comm));
      for (i=j=0; i < BI_MaxNSysCtxt; i++)
      {
         if (BI_SysContxts[i] != MPI_COMM_NULL)
            tSysCtxt[j++] = BI_SysContxts[i];
      }
      BI_MaxNSysCtxt -= MAXNSYSCTXT;
      for(; j < BI_MaxNSysCtxt; j++) tSysCtxt[j] = MPI_COMM_NULL;
      free(BI_SysContxts);
      BI_SysContxts = tSysCtxt;
   }
#endif
}
