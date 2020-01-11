#include "Bdef.h"

#if (INTFACE == C_CALL)
Int Ckrecvid(Int ConTxt, Int rsrc, Int csrc)
#else
F_INT_FUNC krecvid_(Int *ConTxt, Int *rsrc, Int *csrc)
#endif
{
   return(PT2PTID+1);
}  /* end krecvid */
