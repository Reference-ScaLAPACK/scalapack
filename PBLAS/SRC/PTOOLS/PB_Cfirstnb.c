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
Int PB_Cfirstnb( Int N, Int I, Int INB, Int NB )
#else
Int PB_Cfirstnb( N, I, INB, NB )
/*
*  .. Scalar Arguments ..
*/
   Int            I, INB, N, NB;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cfirstnb returns the global number of matrix rows or columns of the
*  first block, if  N  rows or columns  are given out  starting from the
*  global index I. Note that if N is equal 0, this routine returns 0.
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
   Int            inbt;
/* ..
*  .. Executable Statements ..
*
*/
   inbt = ( ( INB -= I ) <= 0 ? ( (-INB) / NB + 1 ) * NB + INB : INB );
   return( MIN( inbt, N ) );
/*
*  End of PB_Cfirstnb
*/
}
