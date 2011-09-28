#include "Bdef.h"

#if (INTFACE == C_CALL)
int Cksendid(int ConTxt, int rdest, int cdest)
#else
F_INT_FUNC ksendid_(int *ConTxt, int *rdest, int *cdest)
#endif
{
   return(PT2PTID+1);
}  /* end ksendid */
