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
void PB_Ctzher( PBTYP_T * TYPE, char * UPLO, Int M, Int N, Int K,
                Int IOFFD, char * ALPHA, char * XC, Int LDXC, char * XR,
                Int LDXR, char * A, Int LDA )
#else
void PB_Ctzher( TYPE, UPLO, M, N, K, IOFFD, ALPHA, XC, LDXC, XR, LDXR,
               A, LDA )
/*
*  .. Scalar Arguments ..
*/
   char           * UPLO;
   Int            IOFFD, K, LDA, LDXC, LDXR, M, N;
   char           * ALPHA;
/*
*  .. Array Arguments ..
*/
   PBTYP_T        * TYPE;
   char           * A, * XC, * XR;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Ctzher performs the trapezoidal symmetric or Hermitian rank 1 ope-
*  ration:
*
*     A := alpha * XC * XR + A  or  A := alpha * XC * conjg( XR ) + A,
*
*  where  alpha  is a scalar, XC is an m element vector, XR is an n ele-
*  ment vector and A is an m by n trapezoidal symmetric or Hermitian ma-
*  trix.
*
*  Arguments
*  =========
*
*  TYPE    (local input) pointer to a PBTYP_T structure
*          On entry,  TYPE  is a pointer to a structure of type PBTYP_T,
*          that contains type information (see pblas.h).
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
*  A       (input/output) pointer to CHAR
*          On entry, A is an array of dimension (LDA,N) containing the m
*          by n matrix A. Only the trapezoidal part of  A  determined by
*          UPLO and IOFFD is updated.
*
*  LDA     (input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least max( 1, M ).
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
   Int            i1, ione=1, j1, m1, mn, n1, size;
   GERC_T         gerc;
/* ..
*  .. Executable Statements ..
*
*/
   if( ( M <= 0 ) || ( N <= 0 ) ) return;

   if( Mupcase( UPLO[0] ) == CLOWER )
   {
      size = TYPE->size; gerc = TYPE->Fgerc;
      mn   = MAX( 0, -IOFFD );
      if( ( n1 = MIN( mn, N ) ) > 0 )
         gerc( &M, &n1, ALPHA, XC, &ione, XR, &LDXR, A, &LDA );
      n1 = M - IOFFD;
      if( ( n1 = MIN( n1, N ) - mn ) > 0 )
      {
         i1 = ( j1 = mn ) + IOFFD;
         TYPE->Fher( C2F_CHAR( UPLO ), &n1, ALPHA, Mptr( XC, i1, 0, LDXC,
                     size ), &ione, Mptr( A, i1, j1, LDA, size ), &LDA );
         if( ( m1 = M - mn - n1 - IOFFD ) > 0 )
         {
            i1 += n1;
            gerc( &m1, &n1, ALPHA, Mptr( XC, i1, 0, LDXC, size ), &ione,
                  Mptr( XR, 0, j1, LDXR, size ), &LDXR, Mptr( A, i1, j1, LDA,
                  size ), &LDA );
         }
      }
   }
   else if( Mupcase( UPLO[0] ) == CUPPER )
   {
      size = TYPE->size; gerc = TYPE->Fgerc;
      mn   = M - IOFFD; mn = MIN( mn, N );
      if( ( n1 = mn - MAX( 0, -IOFFD ) ) > 0 )
      {
         j1 = mn - n1;
         if( ( m1 = MAX( 0, IOFFD ) ) > 0 )
            gerc( &m1, &n1, ALPHA, XC, &ione, XR, &LDXR, A, &LDA );
         TYPE->Fher( C2F_CHAR( UPLO ), &n1, ALPHA, Mptr( XC, m1, 0, LDXC,
                     size ), &ione, Mptr( A, m1, j1, LDA, size ), &LDA );
      }
      if( ( n1 = N - MAX( 0, mn ) ) > 0 )
      {
         j1 = N - n1;
         gerc( &M, &n1, ALPHA, XC, &ione, Mptr( XR, 0, j1, LDXR, size ), &LDXR,
               Mptr( A, 0, j1, LDA, size ), &LDA );
      }
   }
   else
   {
      TYPE->Fgerc( &M, &N, ALPHA, XC, &ione, XR, &LDXR, A, &LDA );
   }
/*
*  End of PB_Ctzher
*/
}
