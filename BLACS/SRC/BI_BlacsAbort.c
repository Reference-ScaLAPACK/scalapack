#include "Bdef.h"

void BI_BlacsAbort(Int ErrNo)
{
   Int ierr;
   fflush(stderr);
   fflush(stdout);
   ierr=MPI_Abort(MPI_COMM_WORLD, ErrNo);
}
