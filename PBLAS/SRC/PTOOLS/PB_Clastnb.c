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
int PB_Clastnb( int N, int I, int INB, int NB )
#else
int PB_Clastnb( N, I, INB, NB )
/*
*  .. Scalar Arguments ..
*/
   int            I, INB, N, NB;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Clastnb  returns  the global number of matrix rows or columns of the
*  last block, if N rows or columns are given out starting from the glo-
*  bal index I. Note that if N is equal 0, this routine returns 0.
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of rows/columns being dealt
*          out. N must be at least zero.
*
*  I       (global input) INTEGER
*          On entry, I  specifies the global index of the matrix  entry.
*          I must be at least zero.
*
*  INB     (global input) INTEGER
*          On entry,  INB  specifies  the size of the first block of the
*          global matrix distribution. INB must be at least one.
*
*  NB      (global input) INTEGER
*          On entry, NB specifies the size of the blocks used to  parti-
*          tion the matrix. NB must be at least one.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   int            lnbt;
/* ..
*  .. Executable Statements ..
*
*/
   if( ( lnbt = I + N - INB ) > 0 )
   {
      lnbt = lnbt - NB * ( ( NB + lnbt - 1 ) / NB - 1 );
      return( MIN( lnbt, N ) );
   }
   else
   {
      return( N );
   }
/*
*  End of PB_Clastnb
*/
}
