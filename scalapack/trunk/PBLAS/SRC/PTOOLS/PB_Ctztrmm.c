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
void PB_Ctztrmm( PBTYP_T * TYPE, char * SIDE, char * UPLO, char * TRANS,
                 char * DIAG, int M, int N, int K, int IOFFD,
                 char * ALPHA, char * A, int LDA, char * B, int LDB,
                 char * C, int LDC )
#else
void PB_Ctztrmm( TYPE, SIDE, UPLO, TRANS, DIAG, M, N, K, IOFFD, ALPHA,
                 A, LDA, B, LDB, C, LDC )
/*
*  .. Scalar Arguments ..
*/
   char               * SIDE, * UPLO, * TRANS, * DIAG;
   int                IOFFD, K, LDA, LDB, LDC, M, N;
   char               * ALPHA;
/*
*  .. Array Arguments ..
*/
   PBTYP_T            * TYPE;
   char               * A, * B, * C;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Ctztrmm  performs the matrix-matrix  operation
*
*     C := alpha * op( A ) * B,
*
*     or
*
*     C := alpha * B * op( A ),
*
*  where  alpha  is a scalar, A  is an m by n trapezoidal triangular ma-
*  trix, B  is an k by n  matrix  when TRANS is 'N' or 'n' and an m by k
*  matrix otherwise, and op( A ) is one of
*
*     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
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
*          SIDE = 'L' or 'l'  C := alpha * op( A ) * B,
*
*          SIDE = 'R' or 'r'  C := alpha * B * op( A ).
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
*  TRANS   (input) pointer to CHAR
*          On entry,  TRANS  specifies the form of op( A ) to be used as
*          follows:
*
*             TRANS = 'N' or 'n':  op( A ) = A,
*
*             TRANS = 'T' or 't':  op( A ) = A',
*
*             TRANS = 'C' or 'c':  op( A ) = A' or conjg( A' ).
*
*  DIAG    (input) pointer to CHAR
*          On entry, DIAG  specifies whether or not A is unit triangular
*          as follows:
*
*             DIAG = 'U' or 'u'  A is assumed to be unit triangular.
*
*             DIAG = 'N' or 'n'  A is not assumed to be unit triangular.
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
*          On entry, K specifies the number of rows of the matrix B when
*          TRANS is 'N' or 'n', and the number  of columns of the matrix
*          B otherwise. K must be at least zero.
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
*          UPLO  and  IOFFD is referenced.  When  DIAG = 'U' or 'u', the
*          diagonal elements of  A  are  not  referenced either, but are
*          assumed to be unity.
*
*  LDA     (input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least max( 1, M ).
*
*  B       (input) pointer to CHAR
*          On entry, B is an array of dimension (LDB,Kb).  Before entry,
*          with  TRANS = 'N' or 'n', the array B must contain the k by n
*          matrix B corresponding to the columns of  A.  Otherwise,  the
*          array B must contain the m by k matrix B corresponding to the
*          rows of A.  When  TRANS is 'N' or 'n', LDB is at least K, and
*          Kb is at least N. Otherwise, LDB is at least max(1,M), and Kb
*          is at least K.
*
*  LDB     (input) INTEGER
*          On entry, LDB specifies the leading dimension of the array B.
*          LDB  must  be  at  least  K  when  TRANS  is  'N' or 'n'  and
*          max( 1, M ) otherwise.
*
*  C       (input/output) pointer to CHAR
*          On entry, C is an array of dimension (LDC,Kc). On exit,  with
*          TRANS = 'N' or 'n', the array C contains the m by k matrix  C
*          corresponding  to  the rows of A. Otherwise, the array C con-
*          tains the k by n matrix C corresponding to the  columns of A.
*          When TRANS is 'N' or 'n', LDC is at least max( 1, M ), and Kc
*          is at least K. Otherwise, LDC is at least K,  and  Kc  is  at
*          least N. On exit,  C  is  overwritten  by the partial updated
*          matrix C.
*
*  LDC     (input) INTEGER
*          On entry, LDC specifies the leading dimension of the array C.
*          LDC  must  be  at  least  max( 1, M ) when  TRANS  is  'N' or
*          'n' and 1 otherwise.
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
   char               * Aptr = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   if( ( M <= 0 ) || ( N <= 0 ) ) return;

   if( ( Mupcase( UPLO[0] ) == CLOWER ) || ( Mupcase( UPLO[0] ) == CUPPER ) )
   {
      Aptr = PB_Cmalloc( M * N * TYPE->size );
      TYPE->Ftzpadcpy( C2F_CHAR( UPLO ), C2F_CHAR( DIAG ), &M, &N, &IOFFD,
                       A, &LDA, Aptr, &M );
      if( Mupcase( SIDE[0] ) == CLEFT )
      {
         if( Mupcase( TRANS[0] ) == CNOTRAN )
         {
            TYPE->Fgemm( C2F_CHAR( TRANS ), C2F_CHAR( TRAN   ), &M, &K, &N,
                         ALPHA, Aptr, &M, B, &LDB, TYPE->one, C, &LDC );
         }
         else
         {
            TYPE->Fgemm( C2F_CHAR( TRANS ), C2F_CHAR( NOTRAN ), &K, &N, &M,
                         ALPHA, B, &LDB, Aptr, &M, TYPE->one, C, &LDC );
         }
      }
      else
      {
         if( Mupcase( TRANS[0] ) == CNOTRAN )
         {
            TYPE->Fgemm( C2F_CHAR( TRAN   ), C2F_CHAR( TRANS ), &K, &N, &M,
                         ALPHA, B, &LDB, Aptr, &M, TYPE->one, C, &LDC );
         }
         else
         {
            TYPE->Fgemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRANS ), &M, &K, &N,
                         ALPHA, Aptr, &M, B, &LDB, TYPE->one, C, &LDC );
         }
      }
      if( Aptr ) free( Aptr );
   }
   else
   {
      if( Mupcase( SIDE[0] ) == CLEFT )
      {
         if( Mupcase( TRANS[0] ) == CNOTRAN )
         {
            TYPE->Fgemm( C2F_CHAR( TRANS ), C2F_CHAR( TRAN   ), &M, &K, &N,
                         ALPHA, A, &LDA, B, &LDB, TYPE->one, C, &LDC );
         }
         else
         {
            TYPE->Fgemm( C2F_CHAR( TRANS ), C2F_CHAR( NOTRAN ), &K, &N, &M,
                         ALPHA, B, &LDB, A, &LDA, TYPE->one, C, &LDC );
         }
      }
      else
      {
         if( Mupcase( TRANS[0] ) == CNOTRAN )
         {
            TYPE->Fgemm( C2F_CHAR( TRAN   ), C2F_CHAR( TRANS ), &K, &N, &M,
                         ALPHA, B, &LDB, A, &LDA, TYPE->one, C, &LDC );
         }
         else
         {
            TYPE->Fgemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRANS ), &M, &K, &N,
                         ALPHA, A, &LDA, B, &LDB, TYPE->one, C, &LDC );
         }
      }
   }
/*
*  End of PB_Ctztrmm
*/
}
