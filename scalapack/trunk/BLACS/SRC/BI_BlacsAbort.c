#include "Bdef.h"

void BI_BlacsAbort(int ErrNo)
{
   int ierr;
   fflush(stderr);
   fflush(stdout);
   ierr=MPI_Abort(MPI_COMM_WORLD, ErrNo);
}
