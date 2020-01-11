#include "Bdef.h"
MPI_Datatype BI_GetMpiGeType(BLACSCONTEXT *ctxt, Int m, Int n, Int lda,
                                MPI_Datatype Dtype, Int *N)
{
   Int info;
   MPI_Datatype GeType;

/*
 * Some versions of mpich and its derivitives cannot handle 0 byte typedefs,
 * so we set type MPI_BYTE as a flag for a 0 byte message
 */
#ifdef ZeroByteTypeBug
   if ( (m < 1) || (n < 1) )
   {
      *N = 0;
      return (MPI_BYTE);
   }
#endif
   *N = 1;
   info=MPI_Type_vector(n, m, lda, Dtype, &GeType);
   info=MPI_Type_commit(&GeType);

   return(GeType);
}
