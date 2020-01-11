#include "Bdef.h"

#if (INTFACE == C_CALL)
Int Cksendid(Int ConTxt, Int rdest, Int cdest)
#else
F_INT_FUNC ksendid_(Int *ConTxt, Int *rdest, Int *cdest)
#endif
{
   return(PT2PTID+1);
}  /* end ksendid */
