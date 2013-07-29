#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_gridmap(int *ConTxt, int *usermap, int ldup, int nprow0, int npcol0)
#else
F_VOID_FUNC blacs_gridmap_(int *ConTxt, int *usermap, int *ldup, int *nprow0,
                           int *npcol0)
#endif
{

   MPI_Comm Cblacs2sys_handle(int BlacsCtxt);
   MPI_Comm BI_TransUserComm(int, int, int *);

   int info, i, j, Iam, *iptr;
   int myrow, mycol, nprow, npcol, Ng;
   BLACSCONTEXT *ctxt, **tCTxts;
   MPI_Comm comm, tcomm;
   MPI_Group grp, tgrp;

   extern BLACSCONTEXT **BI_MyContxts;
   extern BLACBUFF BI_AuxBuff;
   extern int BI_Iam, BI_Np, BI_MaxNCtxt;
   extern MPI_Status *BI_Stats;

/*
 * If first call to blacs_gridmap
 */
   if (BI_MaxNCtxt == 0)
   {
      Cblacs_pinfo(&BI_Iam, &BI_Np);
      BI_AuxBuff.nAops = 0;
      BI_AuxBuff.Aops = (MPI_Request*)malloc(BI_Np*sizeof(*BI_AuxBuff.Aops));
      BI_Stats = (MPI_Status *) malloc(BI_Np * sizeof(MPI_Status));
   }

   nprow = Mpval(nprow0);
   npcol = Mpval(npcol0);
   Ng = nprow * npcol;
   if ( (Ng > BI_Np) || (nprow < 1) || (npcol < 1) )
      BI_BlacsErr(-1, -1, "BLACS_GRIDINIT/BLACS_GRIDMAP",
                  "Illegal grid (%d x %d), #procs=%d", nprow, npcol, BI_Np);
/*
 * Form MPI communicator for scope = 'all'
 */
   if (Ng > 2) i = Ng;
   else i = 2;
   iptr = (int *) malloc(i*sizeof(int));
   for (j=0; j < npcol; j++)
   {
      for (i=0; i < nprow; i++) iptr[i*npcol+j] = usermap[j*Mpval(ldup)+i];
   }
#if (INTFACE == C_CALL)
   tcomm = Cblacs2sys_handle(*ConTxt);
   MPI_Comm_group(tcomm, &grp);           /* find input comm's group */
   MPI_Group_incl(grp, Ng, iptr, &tgrp);  /* form new group */
   MPI_Comm_create(tcomm, tgrp, &comm);   /* create new comm */
   MPI_Group_free(&tgrp);
   MPI_Group_free(&grp);
#else  /* gridmap called from fortran */
   comm = BI_TransUserComm(*ConTxt, Ng, iptr);
#endif

/*
 * Weed out callers who are not participating in present grid
 */
   if (comm == MPI_COMM_NULL)
   {
      *ConTxt = NOTINCONTEXT;
      free(iptr);
      return;
   }

/*
 * ==================================================
 * Get new context and add it to my array of contexts
 * ==================================================
 */
   ctxt = (BLACSCONTEXT *) malloc(sizeof(BLACSCONTEXT));
/*
 * Find free slot in my context array
 */
   for (i=0; i < BI_MaxNCtxt; i++) if (BI_MyContxts[i] == NULL) break;
/*
 * Get bigger context pointer array, if needed
 */
   if (i == BI_MaxNCtxt)
   {
      j = BI_MaxNCtxt + MAXNCTXT;
      tCTxts = (BLACSCONTEXT **) malloc(j * sizeof(*tCTxts));
      for (i=0; i < BI_MaxNCtxt; i++) tCTxts[i] = BI_MyContxts[i];
      BI_MaxNCtxt = j;
      for(j=i; j < BI_MaxNCtxt; j++) tCTxts[j] = NULL;
      if (BI_MyContxts) free(BI_MyContxts);
      BI_MyContxts = tCTxts;
   }
   BI_MyContxts[i] = ctxt;
   *ConTxt = i;

   ctxt->ascp.comm = comm;
   MPI_Comm_dup(comm, &ctxt->pscp.comm); /* copy acomm for pcomm */
   MPI_Comm_rank(comm, &Iam);            /* find my rank in new comm */
   myrow = Iam / npcol;
   mycol = Iam % npcol;

/*
 * Form MPI communicators for scope = 'row'
 */
   MPI_Comm_split(comm, myrow, mycol, &ctxt->rscp.comm);
/*
 * Form MPI communicators for scope = 'Column'
 */
   MPI_Comm_split(comm, mycol, myrow, &ctxt->cscp.comm);

   ctxt->rscp.Np = npcol;
   ctxt->rscp.Iam = mycol;
   ctxt->cscp.Np = nprow;
   ctxt->cscp.Iam = myrow;
   ctxt->pscp.Np = ctxt->ascp.Np = Ng;
   ctxt->pscp.Iam = ctxt->ascp.Iam = Iam;
   ctxt->Nr_bs = ctxt->Nr_co = 1;
   ctxt->Nb_bs = ctxt->Nb_co = 2;
   ctxt->TopsRepeat = ctxt->TopsCohrnt = 0;

/*
 * ===========================
 * Set up the message id stuff
 * ===========================
 */
   Cblacs_get(-1, 1, iptr);
   ctxt->pscp.MinId = ctxt->rscp.MinId = ctxt->cscp.MinId =
   ctxt->ascp.MinId = ctxt->pscp.ScpId = ctxt->rscp.ScpId =
   ctxt->cscp.ScpId = ctxt->ascp.ScpId = iptr[0];
   ctxt->pscp.MaxId = ctxt->rscp.MaxId = ctxt->cscp.MaxId =
   ctxt->ascp.MaxId = iptr[1];
   free(iptr);

}
