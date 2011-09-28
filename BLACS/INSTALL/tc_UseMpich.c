#include <stdio.h>
#include <mpi.h>
main()
{
   MPI_Comm ccomm;
   int fcomm;
   extern void *MPIR_ToPointer();
   extern int   MPIR_FromPointer();
   extern void *MPIR_RmPointer();

   if (sizeof(int) < sizeof(int*))
   {
      fcomm = MPIR_FromPointer(MPI_COMM_WORLD);
      ccomm = (MPI_Comm) MPIR_ToPointer(fcomm);
      if (ccomm == MPI_COMM_WORLD)
         printf("Set TRANSCOMM = -DUseMpich -DPOINTER_64_BITS=1\n");
      else
         printf("Do _NOT_ set TRANSCOMM = -DUseMpich -DPOINTER_64_BITS=1\n");
   }
   else
   {
      printf("Compile and run xtc_CsameF77 for correct TRANSCOMM setting.\n");
      printf("If xtc_CsameF77 fails, leave TRANSCOMM blank.\n");
   }
}
