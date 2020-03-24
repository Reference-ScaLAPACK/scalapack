/* ---------------------------------------------------------------------
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*  ---------------------------------------------------------------------
*/
/*
*  Include files
*/
#include "../pblas.h"
#include "../PBpblas.h"
#include "../PBtools.h"
#include "../PBblacs.h"
#include "../PBblas.h"

#ifdef __STDC__
void PB_Ctzsyr2k( PBTYP_T * TYPE, char * UPLO, Int M, Int N, Int K,
                 Int IOFFD, char * ALPHA, char * AC, Int LDAC,
                 char * BC, Int LDBC, char * AR, Int LDAR, char * BR,
                 Int LDBR, char * C, Int LDC )
#else
void PB_Ctzsyr2k( TYPE, UPLO, M, N, K, IOFFD, ALPHA, AC, LDAC, BC, LDBC,
                AR, LDAR, BR, LDBR, C, LDC )
/*
*  .. Scalar Arguments ..
*/
   char           * UPLO;
   Int            IOFFD, K, LDAC, LDAR, LDBC, LDBR, LDC, M, N;
   char           * ALPHA;
/*
*  .. Array Arguments ..
*/
   char           * AC, * AR, * BC, * BR, * C;
   PBTYP_T        * TYPE;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Ctzsyr2k  performs  the  trapezoidal symmetric o Hermitian rank 2k
*  operation:
*
*     C := alpha * AC * BR + alpha * BC * AR + C,
*
*  where alpha is a scalar, AC and BC are m by k matrices, AR and BR are
*  n by k matrices and C is an m by n trapezoidal symmetric or Hermitian
*  matrix.
*
*  Arguments
*  =========
*
*  TYPE    (local input) pointer to a PBTYP_T structure
*          On entry,  TYPE  is a pointer to a structure of type PBTYP_T,
*          that contains type information (See pblas.h).
*
*  UPLO    (input) pointer to CHAR
*          On entry, UPLO  specifies which part of the matrix C is to be
*          referenced as follows:
*
*             UPLO = 'L' or 'l' the lower trapezoid of C is referenced,
*
*             UPLO = 'U' or 'u' the upper trapezoid of C is referenced,
*
*             otherwise         all of the matrix C is referenced.
*
*  M       (input) INTEGER
*          On entry,  M  specifies the number of rows of the matrix C. M
*          must be at least zero.
*
*  N       (input) INTEGER
*          On entry, N  specifies the number of columns of the matrix C.
*          N must be at least zero.
*
*  K       (input) INTEGER
*          On entry,  K  specifies the number of columns of the matrices
*          AC and BC, and the number of rows of the matrices AR and  BR.
*          K must be at least zero.
*
*  IOFFD   (input) INTEGER
*          On entry, IOFFD specifies the position of the offdiagonal de-
*          limiting the upper and lower trapezoidal part of C as follows
*          (see the notes below):
*
*             IOFFD = 0  specifies the main diagonal C( i, i ),
*                        with i = 1 ... MIN( M, N ),
*             IOFFD > 0  specifies the subdiagonal   C( i+IOFFD, i ),
*                        with i = 1 ... MIN( M-IOFFD, N ),
*             IOFFD < 0  specifies the superdiagonal C( i, i-IOFFD ),
*                        with i = 1 ... MIN( M, N+IOFFD ).
*
*  ALPHA   (input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.
*
*  AC      (input) pointer to CHAR
*          On entry, AC is an array of dimension (LDAC,K) containing the
*          m by k matrix AC.
*
*  LDAC    (input) INTEGER
*          On entry,  LDAC  specifies the leading dimension of the array
*          AC. LDAC must be at least max( 1, M ).
*
*  BC      (input) pointer to CHAR
*          On entry, BC is an array of dimension (LDBC,K) containing the
*          m by k matrix BC.
*
*  LDBC    (input) INTEGER
*          On entry,  LDBC  specifies the leading dimension of the array
*          BC. LDBC must be at least max( 1, M ).
*
*  AR      (input) pointer to CHAR
*          On entry, AR is an array of dimension (LDAR,N) containing the
*          k by n matrix AR.
*
*  LDAR    (input) INTEGER
*          On entry,  LDAR  specifies the leading dimension of the array
*          AR. LDAR must be at least K.
*
*  BR      (input) pointer to CHAR
*          On entry, BR is an array of dimension (LDBR,N) containing the
*          k by n matrix BR.
*
*  LDBR    (input) INTEGER
*          On entry,  LDBR  specifies the leading dimension of the array
*          BR. LDBR must be at least K.
*
*  C       (input/output) pointer to CHAR
*          On entry, C is an array of dimension (LDC,N) containing the m
*          by n matrix C. Only the trapezoidal part of  C  determined by
*          UPLO and IOFFD is updated.
*
*  LDC     (input) INTEGER
*          On entry, LDC specifies the leading dimension of the array C.
*          LDC must be at least max( 1, M ).
*
*  Notes
*  =====
*                           N                                    N
*             ----------------------------                  -----------
*            |       d                    |                |           |
*          M |         d         Upper    |                |    Upper  |
*            | Lower     d                |                |d          |
*            |             d              |              M |  d        |
*             ----------------------------                 |    d      |
*                                                          |      d    |
*               IOFFD < 0                                  | Lower  d  |
*                                                          |          d|
*                  N                                       |           |
*             -----------                                   -----------
*            |    d Upper|
*            |      d    |                                   IOFFD > 0
*          M |        d  |
*            |          d|                              N
*            |  Lower    |                 ----------------------------
*            |           |                |          Upper             |
*            |           |                |d                           |
*            |           |                |  d                         |
*            |           |                |    d                       |
*            |           |                |Lower d                     |
*             -----------                  ----------------------------
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           * one;
   Int            i1, j1, m1, mn, n1, size;
   GEMM_T         gemm;
/* ..
*  .. Executable Statements ..
*
*/
   if( ( M <= 0 ) || ( N <= 0 ) ) return;

   if( Mupcase( UPLO[0] ) == CLOWER )
   {
      size = TYPE->size; one = TYPE->one; gemm = TYPE->Fgemm;
      mn   = MAX( 0, -IOFFD );
      if( ( n1 = MIN( mn, N ) ) > 0 )
      {
         gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &M, &n1, &K, ALPHA, AC,
               &LDAC, BR, &LDBR, one, C, &LDC );
         gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &M, &n1, &K, ALPHA, BC,
               &LDBC, AR, &LDAR, one, C, &LDC );
      }
      n1 = M - IOFFD;
      if( ( n1 = MIN( n1, N ) - mn ) > 0 )
      {
         i1 = ( j1 = mn ) + IOFFD;
         TYPE->Fsyr2k( C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ), &n1, &K, ALPHA,
                       Mptr( AC, i1, 0, LDAC, size ), &LDAC, Mptr( BC, i1, 0,
                       LDBC, size ), &LDBC, one, Mptr( C, i1, j1, LDC, size ),
                       &LDC );
         if( ( m1 = M - mn - n1 - IOFFD ) > 0 )
         {
            i1 += n1;
            gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &m1, &n1, &K,
                         ALPHA, Mptr( AC, i1, 0, LDAC, size ), &LDAC, Mptr( BR,
                         0, j1, LDBR, size ), &LDBR, one, Mptr( C, i1, j1, LDC,
                         size ), &LDC );
            gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &m1, &n1, &K,
                         ALPHA, Mptr( BC, i1, 0, LDBC, size ), &LDBC, Mptr( AR,
                         0, j1, LDAR, size ), &LDAR, one, Mptr( C, i1, j1, LDC,
                         size ), &LDC );
         }
      }
   }
   else if( Mupcase( UPLO[0] ) == CUPPER )
   {
      size = TYPE->size; one = TYPE->one; gemm = TYPE->Fgemm;
      mn   = M - IOFFD;  mn  = MIN( mn, N );
      if( ( n1 = mn - MAX( 0, -IOFFD ) ) > 0 )
      {
         j1 = mn - n1;
         if( ( m1 = MAX( 0, IOFFD ) ) > 0 )
         {
            gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &m1, &n1, &K, ALPHA,
                  AC, &LDAC, BR, &LDBR, one, C, &LDC );
            gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &m1, &n1, &K, ALPHA,
                  BC, &LDBC, AR, &LDAR, one, C, &LDC );
         }
         TYPE->Fsyr2k( C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ), &n1, &K, ALPHA,
                       Mptr( AC, m1, 0, LDAC, size ), &LDAC, Mptr( BC, m1, 0,
                       LDBC, size ), &LDBC, one, Mptr( C, m1, j1, LDC, size ),
                       &LDC );
      }
      if( ( n1 = N - MAX( 0, mn ) ) > 0 )
      {
         j1 = N - n1;
         gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &M, &n1, &K, ALPHA, AC,
               &LDAC, Mptr( BR, 0, j1, LDBR, size ), &LDBR, one, Mptr( C, 0, j1,
               LDC, size ), &LDC );
         gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &M, &n1, &K, ALPHA, BC,
               &LDBC, Mptr( AR, 0, j1, LDAR, size ), &LDAR, one, Mptr( C, 0, j1,
               LDC, size ), &LDC );
      }
   }
   else
   {
      one = TYPE->one; gemm = TYPE->Fgemm;
      gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &M, &N, &K, ALPHA, AC,
            &LDAC, BR, &LDBR, one, C, &LDC );
      gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &M, &N, &K, ALPHA, BC,
            &LDBC, AR, &LDAR, one, C, &LDC );
   }
/*
*  End of PB_Ctzsyr2k
*/
}
