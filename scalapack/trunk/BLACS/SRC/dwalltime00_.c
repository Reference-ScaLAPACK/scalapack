#include "Bdef.h"

#if (INTFACE == C_CALL)
double Cdwalltime00(void)
#else
F_DOUBLE_FUNC dwalltime00_(void)
#endif
{
   return(MPI_Wtime());
}
