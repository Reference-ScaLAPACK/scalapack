#include "Bdef.h"


#if (INTFACE == C_CALL)
void Cblacs_setup(int *mypnum, int *nprocs)
#else
F_VOID_FUNC blacs_setup_(int *mypnum, int *nprocs)
#endif
{
/*
 * blacs_setup same as blacs_pinfo for non-PVM versions of the BLACS
 */
   void Cblacs_pinfo(int *, int *);
   Cblacs_pinfo(mypnum, nprocs);
}
