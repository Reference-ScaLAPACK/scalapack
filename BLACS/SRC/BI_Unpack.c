#include "Bdef.h"

void BI_Unpack(BLACSCONTEXT *ctxt, BVOID *A, BLACBUFF *bp, MPI_Datatype Dtype)
{
   int i=0, info, one=1;

/*
 * Some versions of mpich and its derivitives cannot handle 0 byte typedefs,
 * so we have set MPI_BYTE as a flag for a 0 byte message
 */
#ifdef ZeroByteTypeBug
   if (Dtype == MPI_BYTE) return;
#endif
   // TODO
   // MPI_Unpack_c exists but no need to invoke it for a count of one==1
   // however not clear how to handle bp->Len which may be 64-bit, see BLACS/SRC/Bdef.h:54
   info=MPI_Unpack(bp->Buff, bp->Len, &i, A, one, Dtype, ctxt->scp->comm);
   info=MPI_Type_free(&Dtype);
}
