#include <stdio.h>
#include <mpi.h>

main(int nargs, char **args)
{
   MPI_Datatype Dtype, Dt;
   int i, j, ierr;

   MPI_Init(&nargs, &args);
   printf( "If this routine does not complete, you should set SYSERRORS = -DZeroByteTypeBug.\n");

   i = 0;
   j = 1;
   ierr = MPI_Type_indexed(1, &i, &j, MPI_INT, &Dtype);
   if (ierr == MPI_SUCCESS)
   {
      MPI_Type_commit(&Dtype);
      ierr = MPI_Type_vector(0, 1, 1, MPI_INT, &Dt);
      if (ierr != MPI_SUCCESS)
         printf("MPI_Type_vector returned %d, set SYSERRORS = -DZeroByteTypeBug\n", ierr);
      else MPI_Type_commit(&Dt);
   }
   else printf("MPI_Type_commit returned %d, set SYSERRORS = -DZeroByteTypeBug\n", ierr);
   if (ierr == MPI_SUCCESS) printf("Leave SYSERRORS blank for this system.\n");

   MPI_Finalize();
}
