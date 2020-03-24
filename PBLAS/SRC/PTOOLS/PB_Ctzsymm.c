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
void PB_Ctzsymm( PBTYP_T * TYPE, char * SIDE, char * UPLO, Int M, Int N,
                 Int K, Int IOFFD, char * ALPHA, char * A, Int LDA,
                 char * BC, Int LDBC, char * BR, Int LDBR, char * CC,
                 Int LDCC, char * CR, Int LDCR )
#else
void PB_Ctzsymm( TYPE, SIDE, UPLO, M, N, K, IOFFD, ALPHA, A, LDA, BC,
                 LDBC, BR, LDBR, CC, LDCC, CR, LDCR )
/*
*  .. Scalar Arguments ..
*/
   char           * SIDE, * UPLO;
   Int            IOFFD, K, LDA, LDBC, LDBR, LDCC, LDCR, M, N;
   char           * ALPHA;
/*
*  .. Array Arguments ..
*/
   PBTYP_T        * TYPE;
   char           * A, * BC, * BR, * CC, * CR;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Ctzsymm  performs the matrix-matrix  operation
*
*     C := alpha * A * B + C,
*
*     or
*
*     C := alpha * B * A + C,
*
*  where alpha is a scalar, B and C are m by k and k by n matrices and A
*  is an m by n trapezoidal symmetric or Hermitian matrix.
*
*  Arguments
*  =========
*
*  TYPE    (local input) pointer to a PBTYP_T structure
*          On entry,  TYPE  is a pointer to a structure of type PBTYP_T,
*          that contains type information (See pblas.h).
*
*  SIDE    (input) pointer to CHAR
*          On entry,  SIDE  specifies whether  op( A ) multiplies B from
*          the left or right as follows:
*
*          SIDE = 'L' or 'l'  C := alpha * A * B + C,
*
*          SIDE = 'R' or 'r'  C := alpha * B * A + C.
*
*  UPLO    (input) pointer to CHAR
*          On entry, UPLO  specifies which part of the matrix A is to be
*          referenced as follows:
*
*             UPLO = 'L' or 'l' the lower trapezoid of A is referenced,
*
*             UPLO = 'U' or 'u' the upper trapezoid of A is referenced,
*
*             otherwise         all of the matrix A is referenced.
*
*  M       (input) INTEGER
*          On entry,  M  specifies the number of rows of the matrix A. M
*          must be at least zero.
*
*  N       (input) INTEGER
*          On entry, N  specifies the number of columns of the matrix A.
*          N must be at least zero.
*
*  K       (input) INTEGER
*          On entry,  K  specifies the number of rows of the matrices BR
*          and CR and the number of columns of the matrices BC and CC. K
*          must be at least zero.
*
*  IOFFD   (input) INTEGER
*          On entry, IOFFD specifies the position of the offdiagonal de-
*          limiting the upper and lower trapezoidal part of A as follows
*          (see the notes below):
*
*             IOFFD = 0  specifies the main diagonal A( i, i ),
*                        with i = 1 ... MIN( M, N ),
*             IOFFD > 0  specifies the subdiagonal   A( i+IOFFD, i ),
*                        with i = 1 ... MIN( M-IOFFD, N ),
*             IOFFD < 0  specifies the superdiagonal A( i, i-IOFFD ),
*                        with i = 1 ... MIN( M, N+IOFFD ).
*
*  ALPHA   (input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (input) pointer to CHAR
*          On entry, A is an array of dimension (LDA,N) containing the m
*          by n matrix A. Only the trapezoidal part of  A  determined by
*          UPLO and IOFFD is referenced.
*
*  LDA     (input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least max( 1, M ).
*
*  BC      (input) pointer to CHAR
*          On entry, BC is an array of dimension (LDBC,K) containing the
*          m by k matrix BC.
*
*  LDBC    (input) INTEGER
*          On entry,  LDBC  specifies the leading dimension of the array
*          BC. LDBC must be at least max( 1, M ).
*
*  BR      (input) pointer to CHAR
*          On entry, BR is an array of dimension (LDBR,N) containing the
*          k by n matrix BR.
*
*  LDBR    (input) INTEGER
*          On entry,  LDBR  specifies the leading dimension of the array
*          BR. LDBR must be at least K.
*
*  CC      (input/output) pointer to CHAR
*          On entry, CC is an array of dimension (LDCC,K) containing the
*          m by k matrix CC. On exit, CC is overwritten by the partially
*          updated matric CC.
*
*  LDCC    (input) INTEGER
*          On entry,  LDCC  specifies the leading dimension of the array
*          CC. LDCC must be at least max( 1, M ).
*
*  CR      (input/output) pointer to CHAR
*          On entry, CR is an array of dimension (LDCR,N) containing the
*          k by n matrix CR. On exit, CR is overwritten by the partially
*          updated matrix CR.
*
*  LDCR    (input) INTEGER
*          On entry,  LDCR  specifies the leading dimension of the array
*          CR. LDCR must be at least K.
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

   if( Mupcase( SIDE[0] ) == CLEFT )
   {
      if( Mupcase( UPLO[0] ) == CLOWER )
      {
         size = TYPE->size; one = TYPE->one; gemm = TYPE->Fgemm;
         mn   = MAX( 0, -IOFFD );

         if( ( n1 = MIN( mn, N ) ) > 0 )
         {
            gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &M, &K, &n1, ALPHA, A,
                  &LDA, BR, &LDBR, one, CC, &LDCC );
            gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &K, &n1, &M, ALPHA, BC,
                  &LDBC, A, &LDA, one, CR, &LDCR );
         }
         n1 = M - IOFFD;
         if( ( n1 = MIN( n1, N ) - mn ) > 0 )
         {
            i1 = ( j1 = mn ) + IOFFD;
            TYPE->Fsymm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), &n1, &K, ALPHA,
                         Mptr( A, i1, j1,  LDA, size ), &LDA, Mptr( BC, i1, 0,
                         LDBC, size ), &LDBC, one, Mptr( CC, i1,  0, LDCC,
                         size ), &LDCC );
            if( ( m1 = M - mn - n1 - IOFFD ) > 0 )
            {
               i1 += n1;
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &m1, &K, &n1, ALPHA,
                     Mptr( A, i1, j1, LDA, size ), &LDA, Mptr( BR, 0, j1, LDBR,
                     size ), &LDBR, one, Mptr( CC, i1, 0, LDCC, size ), &LDCC );
               gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &K, &n1, &m1, ALPHA,
                     Mptr( BC, i1, 0, LDBC, size ), &LDBC, Mptr( A, i1, j1, LDA,
                     size ), &LDA, one, Mptr( CR, 0, j1, LDCR, size ), &LDCR );
            }
         }
      }
      else if( Mupcase( UPLO[0] ) == CUPPER )
      {
         size = TYPE->size; one = TYPE->one; gemm = TYPE->Fgemm;
         mn   = MIN( M - IOFFD, N );

         if( ( n1 = mn - MAX( 0, -IOFFD ) ) > 0 )
         {
            j1 = mn - n1;
            if( ( m1 = MAX( 0, IOFFD ) ) > 0 )
            {
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &m1, &K, &n1,
                            ALPHA, A, &LDA, BR, &LDBR, one, CC, &LDCC );
               gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &K, &n1, &m1,
                            ALPHA, BC, &LDBC, A, &LDA, one, CR, &LDCR );
            }
            TYPE->Fsymm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), &n1, &K, ALPHA,
                         Mptr( A,  m1, j1,  LDA, size ), &LDA,
                         Mptr( BC, m1,  0, LDBC, size ), &LDBC, one,
                         Mptr( CC, m1,  0, LDCC, size ), &LDCC );
         }
         if( ( n1 = N - MAX( 0, mn ) ) > 0 )
         {
            j1 = N - n1;
            gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &M, &K, &n1,
                         ALPHA, Mptr( A, 0, j1, LDA, size ), &LDA, Mptr( BR, 0,
                         j1, LDBR, size ), &LDBR, one, CC, &LDCC );
            gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &K, &n1, &M,
                         ALPHA, BC, &LDBC, Mptr( A, 0, j1, LDA, size ), &LDA,
                         one, Mptr( CR, 0, j1, LDCR, size ), &LDCR );
         }
      }
      else
      {
         one = TYPE->one; gemm = TYPE->Fgemm;
         gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &M, &K, &N, ALPHA, A, &LDA,
               BR, &LDBR, one, CC, &LDCC );
         gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &K, &N, &M, ALPHA, BC,
               &LDBC, A, &LDA, one, CR, &LDCR );
      }
   }
   else
   {
      if( Mupcase( UPLO[0] ) == CLOWER )
      {
         size = TYPE->size; one = TYPE->one; gemm = TYPE->Fgemm;
         mn   = MAX( 0, -IOFFD );
         if( ( n1 = MIN( mn, N ) ) > 0 )
         {
            gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &M, &K, &n1, ALPHA, A,
                  &LDA, BR, &LDBR, one, CC, &LDCC );
            gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &K, &n1, &M, ALPHA, BC,
                  &LDBC, A, &LDA, one, CR, &LDCR );
         }
         n1 = M - IOFFD;
         if( ( n1 = MIN( n1, N ) - mn ) > 0 )
         {
            i1 = ( j1 = mn ) + IOFFD;
            TYPE->Fsymm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), &K, &n1, ALPHA,
                         Mptr( A,  i1, j1,  LDA, size ), &LDA,
                         Mptr( BR,  0, j1, LDBR, size ), &LDBR, one,
                         Mptr( CR,  0, j1, LDCR, size ), &LDCR );
            if( ( m1 = M - mn - n1 - IOFFD ) > 0 )
            {
               i1 += n1;
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &m1, &K, &n1, ALPHA,
                     Mptr( A, i1, j1, LDA, size ), &LDA, Mptr( BR, 0, j1, LDBR,
                     size ), &LDBR, one, Mptr( CC, i1, 0, LDCC, size ), &LDCC );
               gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &K, &n1, &m1, ALPHA,
                     Mptr( BC, i1, 0, LDBC, size ), &LDBC, Mptr( A, i1, j1, LDA,
                     size ), &LDA, one, Mptr( CR, 0, j1, LDCR, size ), &LDCR );
            }
         }
      }
      else if( Mupcase( UPLO[0] ) == CUPPER )
      {
         size = TYPE->size; one = TYPE->one; gemm = TYPE->Fgemm;
         mn   = MIN( M - IOFFD, N );
         if( ( n1 = mn - MAX( 0, -IOFFD ) ) > 0 )
         {
            j1 = mn - n1;
            if( ( m1 = MAX( 0, IOFFD ) ) > 0 )
            {
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &m1, &K, &n1, ALPHA,
                     A, &LDA, BR, &LDBR, one, CC, &LDCC );
               gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &K, &n1, &m1, ALPHA,
                     BC, &LDBC, A, &LDA, one, CR, &LDCR );
            }
            TYPE->Fsymm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), &K, &n1, ALPHA,
                         Mptr( A, m1, j1,  LDA, size ), &LDA, Mptr( BR, 0, j1,
                         LDBR, size ), &LDBR, one, Mptr( CR, 0, j1, LDCR,
                         size ), &LDCR );
         }
         if( ( n1 = N - MAX( 0, mn ) ) > 0 )
         {
            j1 = N - n1;
            gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &M, &K, &n1, ALPHA,
                  Mptr( A, 0, j1, LDA, size ), &LDA, Mptr( BR, 0, j1, LDBR,
                  size ), &LDBR, one, CC, &LDCC );
            gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &K, &n1, &M, ALPHA, BC,
                  &LDBC, Mptr( A, 0, j1, LDA, size ), &LDA, one, Mptr( CR, 0,
                  j1, LDCR, size ), &LDCR );
         }
      }
      else
      {
         one = TYPE->one; gemm = TYPE->Fgemm;
         gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &M, &K, &N, ALPHA, A, &LDA,
               BR, &LDBR, one, CC, &LDCC );
         gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &K, &N, &M, ALPHA, BC,
               &LDBC, A, &LDA, one, CR, &LDCR );
      }
   }
/*
*  End of PB_Ctzsymm
*/
}
