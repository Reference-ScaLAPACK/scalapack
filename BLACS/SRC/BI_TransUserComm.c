#include "Bdef.h"

MPI_Comm BI_TransUserComm(int Ucomm, int Np, int *pmap)
{
   MPI_Comm bcomm, ucomm;
   MPI_Group bgrp, ugrp;
   int i;
   ucomm = MPI_Comm_f2c(Ucomm);
   i=MPI_Comm_group(ucomm, &ugrp);
   i=MPI_Group_incl(ugrp, Np, pmap, &bgrp);
   i=MPI_Comm_create(ucomm, bgrp, &bcomm);
   i=MPI_Group_free(&ugrp);
   i=MPI_Group_free(&bgrp);

   return(bcomm);
}
