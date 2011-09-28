#include "Bdef.h"

void BI_ivmcopy(int m, int n, int *A, int lda, int *buff)
/*
 *  performs an vector to matrix copy (unpack) for the data type int
 */
{
   int i, j;

   if ( (m == lda) || (n == 1) )
   {
      m = n * m;
      for (i=0; i < m; i++) A[i] = buff[i];
   }
   else if (m == 1)
   {
      for (j=0; j < n; j++) A[j*lda] = buff[j];
   }
   else
   {
      for (j=0; j< n; j++)
      {
         for (i=0; i < m; i++) A[i] = buff[i];
         A += lda;
         buff += m;
      }
   }
}
