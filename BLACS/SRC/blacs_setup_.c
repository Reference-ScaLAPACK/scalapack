#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cblacs_setup(Int *mypnum, Int *nprocs)
#else
F_VOID_FUNC blacs_setup_(Int *mypnum, Int *nprocs)
#endif
{
/*
 * blacs_setup same as blacs_pinfo for non-PVM versions of the BLACS
 */
   void Cblacs_pinfo(Int *, Int *);
   Cblacs_pinfo(mypnum, nprocs);
}
