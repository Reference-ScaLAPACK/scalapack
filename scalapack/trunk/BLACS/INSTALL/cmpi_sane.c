#include <stdio.h>
#include "mpi.h"
/*
 * Increase/decrease this value to test if a process of a particular size can
 * be spawned to a particular machine
 */
#define WASTE_SIZE 100
#define NPROC 4
main(int narg, char **args)
/*
 * This program checks to make sure that you can run a basic program on your
 * machine using MPI.  Can increase WASTE_SIZE if you think size of executable
 * may be causing launching problems.
 */
{
   int i, Iam, Np;
   int irank[NPROC];
   double WasteOfSpace[WASTE_SIZE];
   MPI_Comm  mcom;
   MPI_Group wgrp, mgrp;
   MPI_Status stat;

   MPI_Init(&narg, &args);
   MPI_Comm_size(MPI_COMM_WORLD, &Np);
   if (Np < NPROC)
   {
      fprintf(stderr, "Not enough processes to run sanity check; need %d, but I've only got %d\n", NPROC, Np);
      MPI_Abort(MPI_COMM_WORLD, -1);
   }

   for (i=0; i != WASTE_SIZE; i++) WasteOfSpace[i] = 0.0;  /* page in Waste */
/*
 * Form context with NPROC members
 */
   for (i=0; i != NPROC; i++) irank[i] = i;
   MPI_Comm_group(MPI_COMM_WORLD, &wgrp);
   MPI_Group_incl(wgrp, NPROC, irank, &mgrp);
   MPI_Comm_create(MPI_COMM_WORLD, mgrp, &mcom);
   MPI_Group_free(&mgrp);
/*
 * Everyone in new communicator sends a message to his neighbor
 */
   if (mcom != MPI_COMM_NULL)
   {
      MPI_Comm_rank(mcom, &Iam);
/*
 *    Odd nodes receive first, so we don't hang if MPI_Send is globally blocking
 */
      if (Iam % 2)
      {
         MPI_Recv(&i, 1, MPI_INT, (NPROC+Iam-1)%NPROC, 0, mcom, &stat);
         MPI_Send(&Iam, 1, MPI_INT, (Iam+1)%NPROC, 0, mcom);
      }
      else
      {
         MPI_Send(&Iam, 1, MPI_INT, (Iam+1)%NPROC, 0, mcom);
         MPI_Recv(&i, 1, MPI_INT, (NPROC+Iam-1)%NPROC, 0, mcom, &stat);
      }
/*
 *    Make sure we've received the right information
 */
      if (i != (NPROC+Iam-1)%NPROC)
      {
         fprintf(stderr, "Communication does not seem to work properly!!\n");
         MPI_Abort(MPI_COMM_WORLD, -1);
      }
   }
   fprintf(stdout, "%d: C MPI sanity test passed\n", Iam);
   MPI_Finalize();
   exit(0);
}
