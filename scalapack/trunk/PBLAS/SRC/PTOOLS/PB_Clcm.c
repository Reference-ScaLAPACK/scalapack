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
int PB_Clcm( int M, int N )
#else
int PB_Clcm( M, N )
/*
*  .. Scalar Arguments ..
*/
   int            M, N;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Clcm computes and returns the Least Common Multiple  (LCM)  of two
*  positive integers M and N. In fact, the routine computes the Greatest
*  Common Divisor (GCD) and use the property that M*N = GCD*LCM.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          On entry, M must be at least zero.
*
*  N       (input) INTEGER
*          On entry, N must be at least zero.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   int            gcd=1, m_val, n_val, t;
/* ..
*  .. Executable Statements ..
*
*/
   if( M > N ) { m_val = N; n_val = M; }
   else        { m_val = M; n_val = N; }

   while( m_val > 0 )
   {
      while( !( m_val & 1 ) )
      {
/*
*  m is even
*/
         m_val >>= 1;
/*
*  if n is odd, gcd( m, n ) = gcd( m / 2, n )
*/
         if( !( n_val & 1 ) )
         {
/*
*  otherwise gcd( m, n ) = 2 * gcd( m / 2, n / 2 )
*/
            n_val >>= 1;
            gcd   <<= 1;
         }
      }
/*
*  m is odd now. If n is odd, gcd( m, n ) = gcd( m, ( m - n ) / 2 ).
*  Otherwise, gcd( m, n ) = gcd ( m, n / 2 ).
*/
      n_val -= ( n_val & 1 ) ? m_val : 0;
      n_val >>= 1;
      while( n_val >= m_val )
      {
/*
*  If n is odd, gcd( m, n ) = gcd( m, ( m - n ) / 2 ).
*  Otherwise, gcd( m, n ) = gcd ( m, n / 2 )
*/
        n_val -= ( n_val & 1 ) ? m_val : 0;
        n_val >>= 1;
      }
/*
*  n < m, gcd( m, n ) = gcd( n, m )
*/
      t     = n_val;
      n_val = m_val;
      m_val = t;
   }
   return ( ( M * N ) / ( n_val * gcd ) );
/*
*  End of PB_Clcm
*/
}
