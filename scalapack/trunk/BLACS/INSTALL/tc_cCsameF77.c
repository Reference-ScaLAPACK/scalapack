#include <mpi.h>

int Ccommcheck(int F77World, int f77comm)
{
   int Np, Iam, i, OK=1;

   if (sizeof(int) != sizeof(MPI_Comm)) OK=0;
   else if ((MPI_Comm) F77World != MPI_COMM_WORLD) OK=0;
   else
   {
      MPI_Comm_rank(MPI_COMM_WORLD, &Iam);
      if (Iam > 1) OK = ((MPI_Comm) f77comm == MPI_COMM_NULL);
      else
      {
         i = MPI_Comm_size((MPI_Comm) f77comm, &Np);
	 if (i != MPI_SUCCESS) OK = 0;
	 else if (Np != 2) OK = 0;
      }
   }
   MPI_Allreduce(&OK, &i, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   return(i);
}

/*
 * Fortran interfaces
 */
int CCOMMCHECK(int *F77World, int *f77comm)
{
   return(Ccommcheck(*F77World, *f77comm));
}
int ccommcheck_(int *F77World, int *f77comm)
{
   return(Ccommcheck(*F77World, *f77comm));
}
int ccommcheck(int *F77World, int *f77comm)
{
   return(Ccommcheck(*F77World, *f77comm));
}
