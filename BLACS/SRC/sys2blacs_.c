#include "Bdef.h"

#if (INTFACE == C_CALL)
Int Csys2blacs_handle(MPI_Comm SysCtxt)
#else
Int sys2blacs_handle_(Int *SysCtxt)
#endif
{
#if (INTFACE == C_CALL)
   Int i, j, DEF_WORLD;
   MPI_Comm *tSysCtxt;
   extern Int BI_MaxNSysCtxt;
   extern MPI_Comm *BI_SysContxts;

   if (BI_COMM_WORLD == NULL) 
      Cblacs_pinfo(&i, &j);
   if (SysCtxt == MPI_COMM_NULL)
      BI_BlacsErr(-1, __LINE__, __FILE__,
                  "Cannot define a BLACS system handle based on MPI_COMM_NULL");
/*
 * See if we already have this system handle stored
 */
   for (i=0; i < BI_MaxNSysCtxt; i++)
      if (BI_SysContxts[i] == SysCtxt) return(i);
/*
 * The first time in this routine, we need to define MPI_COMM_WORLD, if it isn't
 * what is already being defined.
 */
   DEF_WORLD = ( (!BI_SysContxts) && (SysCtxt != MPI_COMM_WORLD) );
/*
 * Find free slot in system context array
 */
   for (i=0; i < BI_MaxNSysCtxt; i++)
      if (BI_SysContxts[i] == MPI_COMM_NULL) break;
/*
 * If needed, get a bigger system context array
 */
   if (i == BI_MaxNSysCtxt)
   {
      j = BI_MaxNSysCtxt + MAXNSYSCTXT;
      if ( (MAXNSYSCTXT == 1) && (DEF_WORLD) ) j++;
      tSysCtxt = (MPI_Comm *) malloc(j * sizeof(MPI_Comm));
      for (i=0; i < BI_MaxNSysCtxt; i++) tSysCtxt[i] = BI_SysContxts[i];
      BI_MaxNSysCtxt = j;
      for (j=i; j < BI_MaxNSysCtxt; j++) tSysCtxt[j] = MPI_COMM_NULL;
      if (BI_SysContxts) free(BI_SysContxts);
      BI_SysContxts = tSysCtxt;
   }
   if (DEF_WORLD) BI_SysContxts[i++] = MPI_COMM_WORLD;
   BI_SysContxts[i] = SysCtxt;
   return(i);
#else
   return(*SysCtxt);
#endif
}
