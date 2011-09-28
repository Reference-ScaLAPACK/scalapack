#include "Bdef.h"

#if (INTFACE == C_CALL)
int Ckrecvid(int ConTxt, int rsrc, int csrc)
#else
F_INT_FUNC krecvid_(int *ConTxt, int *rsrc, int *csrc)
#endif
{
   return(PT2PTID+1);
}  /* end krecvid */
