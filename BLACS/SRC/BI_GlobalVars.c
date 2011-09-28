#include "Bdef.h"
/*
 * Define global variables
 */
int BI_MaxNCtxt=0;		  /* Number of context pointers allocated */
int BI_MaxNSysCtxt=0;             /* Number of system ctxt ptrs allocated */
int BI_Iam, BI_Np=(-1);	          /* My pnum, and # of procs in system */
BLACBUFF *BI_ReadyB=NULL;         /* buffer that is ready for use */
BLACBUFF *BI_ActiveQ=NULL;        /* pointer to start of active buffer queue */
BLACBUFF BI_AuxBuff;
BLACSCONTEXT **BI_MyContxts=NULL; /* Array of pointers to my contexts */
MPI_Comm *BI_SysContxts=NULL;
int *BI_COMM_WORLD=NULL;
MPI_Status *BI_Stats=NULL;
