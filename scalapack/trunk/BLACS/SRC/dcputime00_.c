#include "Bdef.h"

#if (INTFACE == C_CALL)
double Cdcputime00(void)
#else
F_DOUBLE_FUNC dcputime00_(void)
#endif
{
   return(-1.0);
}
