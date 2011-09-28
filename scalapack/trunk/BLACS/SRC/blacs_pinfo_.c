#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_pinfo(int *mypnum, int *nprocs)
#else
F_VOID_FUNC blacs_pinfo_(int *mypnum, int *nprocs)
#endif
{
   int ierr;
   extern int BI_Iam, BI_Np;
   int argc=0;
   char **argv=NULL;
   if (BI_COMM_WORLD == NULL)
   {
      MPI_Initialized(nprocs);

      if (!(*nprocs)) 
         ierr = MPI_Init(&argc,&argv);  // call Init and ignore argc and argv

      BI_COMM_WORLD = (int *) malloc(sizeof(int));
      *BI_COMM_WORLD = MPI_Comm_c2f(MPI_COMM_WORLD);
      MPI_Comm_size(MPI_COMM_WORLD, &BI_Np);
      MPI_Comm_rank(MPI_COMM_WORLD, &BI_Iam);
   }
   *mypnum = BI_Iam;
   *nprocs = BI_Np;
}
