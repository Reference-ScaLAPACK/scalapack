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
void PB_Ctztrmv( PBTYP_T * TYPE, char * SIDE, char * UPLO, char * TRANS,
                 char * DIAG, Int M, Int N, Int K, Int IOFFD,
                 char * ALPHA, char * A, Int LDA, char * X, Int LDX,
                 char * Y, Int LDY )
#else
void PB_Ctztrmv( TYPE, SIDE, UPLO, TRANS, DIAG, M, N, K, IOFFD, ALPHA,
                 A, LDA, X, LDX, Y, LDY )
/*
*  .. Scalar Arguments ..
*/
   char               * SIDE, * UPLO, * TRANS, * DIAG;
   Int                IOFFD, K, LDA, LDX, LDY, M, N;
   char               * ALPHA;
/*
*  .. Array Arguments ..
*/
   char               * A, * X, * Y;
   PBTYP_T            * TYPE;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Ctztrmv  performs the matrix-vector  operation
*
*     y := A * x,  or  y := A' * x,  or  y := conjg( A' ) * x,
*
*  where  alpha  and beta  are scalars, x and y are vectors, and A is an
*  m by n trapezoidal triangular matrix.
*
*  Arguments
*  =========
*
*  TYPE    (local input) pointer to a PBTYP_T structure
*          On entry,  TYPE  is a pointer to a structure of type PBTYP_T,
*          that contains type information (See pblas.h).
*
*  SIDE    (dummy) pointer to CHAR
*          In this routine, SIDE is a dummy (unused) argument.
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
*          On entry,  TRANS  specifies  the operation to be performed as
*          follows:
*
*             TRANS = 'N' or 'n':  y := A*x,
*
*             TRANS = 'T' or 't':  y := A'*x,
*
*             TRANS = 'C' or 'c':  y := A'*x  or  y := conjg( A' )*x.
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
*  K       (dummy) INTEGER
*          In this routine, K is a dummy (unused) argument.
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
*  ALPHA   (dummy) pointer to CHAR
*          In this routine, ALPHA is a dummy (unused) argument.
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
*  X       (input) pointer to CHAR
*          On entry, X is an array of dimension (LDX,Kx).  Before entry,
*          with  TRANS = 'N' or 'n', the array X must contain the n ele-
*          ment vector x corresponding to the columns of  A.  Otherwise,
*          the array X must contain the m element vector x corresponding
*          to the rows of A.  When  TRANS is 'N' or 'n', LDX is at least
*          1, and Kx is at least N. Otherwise, LDX is at least max(1,M),
*          and Kx is at least 1.
*
*  LDX     (input) INTEGER
*          On entry, LDX specifies the leading dimension of the array X.
*          LDX  must  be  at  least  1  when  TRANS  is  'N' or 'n'  and
*          max( 1, M ) otherwise.
*
*  Y       (input/output) pointer to CHAR
*          On entry, Y is an array of dimension (LDY,Ky). On exit,  with
*          TRANS = 'N' or 'n', the array Y contains the m element vector
*          y corresponding to the rows of A. Otherwise, the array Y con-
*          tains the n element vector y corresponding to the  columns of
*          A. When TRANS is 'N' or 'n', LDY is at least max( 1, M ), and
*          Ky is at least 1. Otherwise, LDY is at least 1, and Ky is  at
*          least N. On exit,  Y  is  overwritten  by the partial updated
*          vector y.
*
*  LDY     (input) INTEGER
*          On entry, LDY specifies the leading dimension of the array Y.
*          LDY  must  be  at  least  max( 1, M ) when  TRANS  is  'N' or
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
   Int                ione = 1;
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
      if( Mupcase( TRANS[0] ) == CNOTRAN )
      {
         TYPE->Fgemv( C2F_CHAR( TRANS ), &M, &N, ALPHA, Aptr, &M, X, &LDX,
                      TYPE->one, Y, &ione );
      }
      else
      {
         TYPE->Fgemv( C2F_CHAR( TRANS ), &M, &N, ALPHA, Aptr, &M, X, &ione,
                      TYPE->one, Y, &LDY );
      }
      if( Aptr ) free( Aptr );
   }
   else
   {
      if( Mupcase( TRANS[0] ) == CNOTRAN )
      {
         TYPE->Fgemv( C2F_CHAR( TRANS ), &M, &N, ALPHA, A, &LDA, X, &LDX,
                      TYPE->one, Y, &ione );
      }
      else
      {
         TYPE->Fgemv( C2F_CHAR( TRANS ), &M, &N, ALPHA, A, &LDA, X, &ione,
                      TYPE->one, Y, &LDY );
      }
   }
/*
*  End of PB_Ctztrmv
*/
}
