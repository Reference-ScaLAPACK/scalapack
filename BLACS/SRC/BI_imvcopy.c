#include "Bdef.h"
void BI_imvcopy(int m, int n, int *A, int lda, int *buff)
/*
 * Performs a matrix to vector copy (pack) for the data type int
 */
{
   int i, j;

   if ( (m == lda) || (n == 1) )
   {
      m = n * m;
      for (i=0; i < m; i++) buff[i] = A[i];
   }
   else if (m == 1)
   {
      for (j=0; j < n; j++) buff[j] = A[j*lda];
   }
   else
   {
      for (j=0; j < n; j++)
      {
         for (i=0; i < m; i++) buff[i] = A[i];
         A += lda;
         buff += m;
      }
   }
}
