#include "Bdef.h"

void BI_Unpack(BLACSCONTEXT *ctxt, BVOID *A, BLACBUFF *bp, MPI_Datatype Dtype)
{
   MpiInt i=0, info, one=1;

/*
 * Some versions of mpich and its derivitives cannot handle 0 byte typedefs,
 * so we have set MPI_BYTE as a flag for a 0 byte message
 */
#ifdef ZeroByteTypeBug
   if (Dtype == MPI_BYTE) return;
#endif
   info=MPI_Unpack(bp->Buff, bp->Len, &i, A, one, Dtype, ctxt->scp->comm);
   info=MPI_Type_free(&Dtype);
}
