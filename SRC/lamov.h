//
//  lamov.h
//
//  Written by Lee Killough 04/19/2012
//  

#include "pblas.h"
#include <ctype.h>

extern void xerbla_(const char *, const F_INTG_FCT *, size_t);

void LACPY(const char *UPLO,
           const F_INTG_FCT *M,
           const F_INTG_FCT *N,
           const TYPE *A,
           const F_INTG_FCT *LDA,
           TYPE *B,
           const F_INTG_FCT *LDB);

void LAMOV(const char *UPLO,
           const F_INTG_FCT *M,
           const F_INTG_FCT *N,
           const TYPE *A,
           const F_INTG_FCT *LDA,
           TYPE *B,
           const F_INTG_FCT *LDB)
{
   const F_INTG_FCT m = *M;
   const F_INTG_FCT n = *N;
   const F_INTG_FCT lda = *LDA;
   const F_INTG_FCT ldb = *LDB;

   if (B + m-1 + ldb*(n-1) < A || A + m-1 + lda*(n-1) < B)
     {
       LACPY(UPLO, M, N, A, LDA, B, LDB);
     }
   else if (lda != ldb)
     {
       TYPE *tmp = malloc(sizeof(*A) * m * n);
       if (!tmp)
         {
           F_INTG_FCT info = -1;
           const char func[] = FUNC;
           xerbla_(func, &info, sizeof func);
         }
       else
         {
           LACPY(UPLO, M, N,   A, LDA, tmp,  &m);
           LACPY(UPLO, M, N, tmp,  &m,   B, LDB);
           free(tmp);
         }
     }
   else
     {
       F_INTG_FCT i, j;
       switch (toupper(*UPLO))
         {
         case 'U':
           if (A > B)
             {
               for (j=0; j<n; j++)
                 for (i=0; i<j && i<m; i++)
                   B[i+ldb*j] = A[i+lda*j];
             }
           else
             {
               for (j=n; --j>=0;)
                 for (i=j<m ? j : m; --i>=0;)
                   B[i+ldb*j] = A[i+lda*j];
             }
           break;
         
         case 'L':
           if (A > B)
             {
               for (j=0; j<n; j++)
                 for (i=j; i<m; i++)
                   B[i+ldb*j] = A[i+lda*j];
             }
           else
             {
               for (j=m<n ? m : n; --j>=0;)
                 for (i=m; --i>=j;)
                   B[i+ldb*j] = A[i+lda*j];
             }
           break;
         
         default:
           if (A > B)
             {
               for (j=0; j<n; j++)
                 for (i=0; i<m; i++)
                   B[i+ldb*j] = A[i+lda*j];
             }
           else
             {
               for (j=n; --j>=0;)
                 for (i=m; --i>=0;)
                   B[i+ldb*j] = A[i+lda*j];
             }
           break;
         }
     }
}
