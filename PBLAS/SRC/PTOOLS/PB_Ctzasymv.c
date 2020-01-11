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
void PB_Ctzasymv( PBTYP_T * TYPE, char * SIDE, char * UPLO, Int M, Int N,
                  Int K, Int IOFFD, char * ALPHA, char * A, Int LDA,
                  char * XC, Int LDXC, char * XR, Int LDXR, char * YC,
                  Int LDYC, char * YR, Int LDYR )
#else
void PB_Ctzasymv( TYPE, SIDE, UPLO, M, N, K, IOFFD, ALPHA, A, LDA, XC,
                  LDXC, XR, LDXR, YC, LDYC, YR, LDYR )
/*
*  .. Scalar Arguments ..
*/
   char           * SIDE, * UPLO;
   Int            IOFFD, K, LDA, LDXC, LDXR, LDYC, LDYR, M, N;
   char           * ALPHA;
/*
*  .. Array Arguments ..
*/
   PBTYP_T        * TYPE;
   char           * A, * XC, * XR, * YC, * YR;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Ctzasymv  performs the matrix-vector  operation
*
*     y := abs( alpha )*abs( A )*abs( x )+ abs( y ),
*
*  where alpha is a real scalar, y is a real vector, x is a vector and A
*  is an m by n trapezoidal symmetric or Hermitian matrix.
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
*  XC      (input) pointer to CHAR
*          On entry, XC is an array of dimension (LDXC,1) containing the
*          m by 1 vector XC.
*
*  LDXC    (input) INTEGER
*          On entry,  LDXC  specifies the leading dimension of the array
*          XC. LDXC must be at least max( 1, M ).
*
*  XR      (input) pointer to CHAR
*          On entry, XR is an array of dimension (LDXR,N) containing the
*          1 by n vector XR.
*
*  LDXR    (input) INTEGER
*          On entry,  LDXR  specifies the leading dimension of the array
*          XR. LDXR must be at least 1.
*
*  YC      (input/output) pointer to CHAR
*          On entry, YC is an array of dimension (LDYC,1) containing the
*          m by 1 vector YC. On exit, YC is overwritten by the partially
*          updated vector y.
*
*  LDYC    (input) INTEGER
*          On entry,  LDYC  specifies the leading dimension of the array
*          YC. LDYC must be at least max( 1, M ).
*
*  YR      (input/output) pointer to CHAR
*          On entry, YR is an array of dimension (LDYR,N) containing the
*          1 by n vector YR. On exit, YR is overwritten by the partially
*          updated vector y.
*
*  LDYR    (input) INTEGER
*          On entry,  LDYR  specifies the leading dimension of the array
*          YR. LDYR must be at least 1.
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
   Int            i1, ione=1, j1, m1, mn, n1, size, usiz;
   AGEMV_T        agemv;
/* ..
*  .. Executable Statements ..
*
*/
   if( ( M <= 0 ) || ( N <= 0 ) ) return;

   if( Mupcase( UPLO[0] ) == CLOWER )
   {
      size = TYPE->size; usiz  = TYPE->usiz;
      one  = TYPE->one;  agemv = TYPE->Fagemv;
      mn   = MAX( 0, -IOFFD );
      if( ( n1 = MIN( mn, N ) ) > 0 )
      {
         agemv( C2F_CHAR( NOTRAN ), &M, &n1, ALPHA, A, &LDA, XR, &LDXR, one, YC,
                &ione );
         agemv( C2F_CHAR( TRAN   ), &M, &n1, ALPHA, A, &LDA, XC, &ione, one, YR,
                &LDYR );
      }
      n1 = M - IOFFD;
      if( ( n1 = MIN( n1, N ) - mn ) > 0 )
      {
         i1 = ( j1 = mn ) + IOFFD;
         TYPE->Fasymv( C2F_CHAR( UPLO ), &n1, ALPHA, Mptr( A, i1, j1, LDA,
                       size ), &LDA, Mptr( XC, i1, 0, LDXC, size ), &ione, one,
                       Mptr( YC, i1, 0, LDYC, usiz ), &ione );
         if( ( m1 = M - mn - n1 - IOFFD ) > 0 )
         {
            i1 += n1;
            agemv( C2F_CHAR( NOTRAN ), &m1, &n1, ALPHA, Mptr( A, i1, j1, LDA,
                   size ), &LDA, Mptr( XR, 0, j1, LDXR, size ), &LDXR, one,
                   Mptr( YC, i1, 0, LDYC, usiz ), &ione );
            agemv( C2F_CHAR( TRAN   ), &m1, &n1, ALPHA, Mptr( A, i1, j1, LDA,
                   size ), &LDA, Mptr( XC, i1, 0, LDXC, size ), &ione, one,
                   Mptr( YR, 0, j1, LDYR, usiz ), &LDYR );
         }
      }
   }
   else if( Mupcase( UPLO[0] ) == CUPPER )
   {
      size = TYPE->size; usiz  = TYPE->usiz;
      one  = TYPE->one;  agemv = TYPE->Fagemv;
      mn   = M - IOFFD;  mn    = MIN( mn, N );
      if( ( n1 = mn - MAX( 0, -IOFFD ) ) > 0 )
      {
         j1 = mn - n1;
         if( ( m1 = MAX( 0, IOFFD ) ) > 0 )
         {
            agemv( C2F_CHAR( NOTRAN ), &m1, &n1, ALPHA, A, &LDA, XR, &LDXR, one,
                   YC, &ione );
            agemv( C2F_CHAR( TRAN   ), &m1, &n1, ALPHA, A, &LDA, XC, &ione, one,
                   YR, &LDYR );
         }
         TYPE->Fasymv( C2F_CHAR( UPLO ), &n1, ALPHA, Mptr( A, m1, j1, LDA,
                       size ), &LDA, Mptr( XC, m1, 0, LDXC, size ), &ione, one,
                       Mptr( YC, m1, 0, LDYC, usiz ), &ione );
      }
      if( ( n1 = N - MAX( 0, mn ) ) > 0 )
      {
         j1 = N - n1;
         agemv( C2F_CHAR( NOTRAN ), &M, &n1, ALPHA, Mptr( A, 0, j1, LDA, size ),
                &LDA, Mptr( XR, 0, j1, LDXR, size ), &LDXR, one, YC, &ione );
         agemv( C2F_CHAR( TRAN   ), &M, &n1, ALPHA, Mptr( A, 0, j1, LDA, size ),
                &LDA, XC, &ione, one, Mptr( YR, 0, j1, LDYR, usiz ), &LDYR );
      }
   }
   else
   {
      one  = TYPE->one;  agemv = TYPE->Fagemv;
      agemv( C2F_CHAR( NOTRAN ), &M, &N, ALPHA, A, &LDA, XR, &LDXR, one, YC,
             &ione );
      agemv( C2F_CHAR( TRAN   ), &M, &N, ALPHA, A, &LDA, XC, &ione, one, YR,
             &LDYR );
   }
/*
*  End of PB_Ctzasymv
*/
}
