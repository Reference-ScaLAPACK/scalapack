#include "Bdef.h"

MPI_Comm BI_TransUserComm(Int Ucomm, int Np, Int *pmap)
{
   MPI_Comm bcomm, ucomm;
   MPI_Group bgrp, ugrp;
   Int i;

   int *mpmap = (int *)malloc(Np * sizeof(int));
   for (i=0; i<Np; i++) mpmap[i] = pmap[i];

   ucomm = MPI_Comm_f2c(Ucomm);
   i=MPI_Comm_group(ucomm, &ugrp);
   i=MPI_Group_incl(ugrp, Np, mpmap, &bgrp);
   i=MPI_Comm_create(ucomm, bgrp, &bcomm);
   i=MPI_Group_free(&ugrp);
   i=MPI_Group_free(&bgrp);

   free(mpmap);

   return(bcomm);
}
