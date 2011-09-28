#include "Bdef.h"


MPI_Datatype BI_GetMpiTrType(BLACSCONTEXT *ctxt, char uplo, char diag,
                                int m, int n, int lda, MPI_Datatype Dtype,
                                int *N)
{
   BLACBUFF *BI_GetBuff(int);
   MPI_Datatype TrType;
   int info, start, i, k;
   int *len, *disp;
   BLACBUFF *bp;

   if (diag == 'u') start = 1;
   else start = 0;

/*
 * Some versions of mpich and its derivitives cannot handle 0 byte typedefs,
 * so we set type MPI_BYTE as a flag for a 0 byte message
 */
#ifdef ZeroByteTypeBug
   if (m > n) i = n * (m-n) + (n*n) - (n*n)/2 + n/2 - n*start;
   else i = m * (n-m) + (m*m) - (m*m)/2 + m/2 - m*start;
   if (i < 1)
   {
      *N = 0;
      return (MPI_BYTE);
   }
#endif
   *N = 1;

/*
 * Get space to hold the length and displacement values
 */
   bp = BI_GetBuff( 2 * n * sizeof(int) );
   len = (int *) bp->Buff;
   disp = (int *) &bp->Buff[n*sizeof(int)];

   if (m > n)
   {
      if (uplo == 'u')
      {
         k = m - n + 1 - start;
         for (i=0; i < n; i++)
         {
            len[i] = k + i;
            disp[i] = i*lda;
         }
      }
      else  /* uplo = 'l' and m > n */
      {
         k = m - start;
         lda++;
         len[0] = k;
         disp[0] = start;
         for (i=1; i < n; i++)
         {
            len[i] = k - i;
            disp[i] = disp[i-1] + lda;
         }
      }
   }
   else /* m <= n */
   {
      if (uplo == 'u')
      {
         k = 1 - start;
         for (i=0; i < m; i++)
         {
            len[i] = i + k;
            disp[i] = i*lda;
         }
         for (; i < n; i++)
         {
            len[i] = m;
            disp[i] = i*lda;
         }
      }
      else  /* uplo = 'l' and m <= n */
      {
         k = n - m;
         for (i=0; i < k; i++)
         {
            len[i] = m;
            disp[i] = i*lda;
         }
         if (i < n)
         {
            k = n - start;
            len[i] = k - i;
            disp[i] = i*lda + start;
            lda++;
            for (i++; i < n; i++)
            {
               len[i] = k - i;
               disp[i] = disp[i-1] + lda;
            }
         }
      }
   }
#ifdef T3ETrError
/*
 * Get rid of 0-length segments to keep T3E happy
 */
   for (i=0; i < n; i++)
   {
      if (len[i] == 0)
      {
         for (k=i+1; k < n; k++)
         {
            len[k-1] = len[k];
            disp[k-1] = disp[k];
         }
         if (n > 0) n--;
         i--;   /* check new entry for 0-byte */
      }
   }
#endif

   i=MPI_Type_indexed(n, len, disp, Dtype, &TrType);
   i=MPI_Type_commit(&TrType);
   return(TrType);
}
