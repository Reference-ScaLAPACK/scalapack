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
int PB_Cnpreroc( int N, int I, int INB, int NB, int PROC, int SRCPROC,
                 int NPROCS )
#else
int PB_Cnpreroc( N, I, INB, NB, PROC, SRCPROC, NPROCS )
/*
*  .. Scalar Arguments ..
*/
   int            I, INB, N, NB, NPROCS, PROC, SRCPROC;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cnpreroc computes  the  number of preceeding rows or columns of a
*  submatrix  that  are  possessed by processes closer to SRCPROC1 than
*  PROC where SRCPROC1 is the process owning the row or column globally
*  indexed by I. The submatrix is defined by giving out N  rows/columns
*  starting from global index I.  Therefore, if  SRCPROC=0 and  PROC=4,
*  then PB_Cnpreroc returns the number of matrix rows or columns  owned
*  by processes 0, 1, 2, and 3.
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
*  PROC    (local input) INTEGER
*          On entry, PROC specifies  the coordinate of the process whose
*          local portion is determined.  PROC must be at least zero  and
*          strictly less than NPROCS.
*
*  SRCPROC (global input) INTEGER
*          On entry,  SRCPROC  specifies  the coordinate of the  process
*          that possesses the  first row or column  of the matrix.  When
*          SRCPROC = -1, the data  is not  distributed  but  replicated,
*          otherwise  SRCPROC  must be at least zero and  strictly  less
*          than NPROCS.
*
*  NPROCS  (global input) INTEGER
*          On entry,  NPROCS  specifies the total number of process rows
*          or columns over which the matrix is distributed.  NPROCS must
*          be at least one.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   int            ilocblk, mydist, nblocks;
/* ..
*  .. Executable Statements ..
*
*/
   if( ( SRCPROC == -1 ) || ( NPROCS == 1 ) )
/*
*  The data is not distributed, or there is just one process in this dimension
*  of the grid.
*/
      return( 0 );
/*
*  Compute coordinate of process owning I and corresponding INB
*/
   if( ( INB -= I ) <= 0 )
   {
/*
*  I is not in first block, find out which process has it and update size of
*  first block
*/
      nblocks  = ( -INB ) / NB + 1;
      SRCPROC += nblocks;
      SRCPROC -= ( SRCPROC / NPROCS ) * NPROCS;
      INB     += nblocks * NB;
   }
/*
*  Now everything is just like N, I=0, INB, NB, SRCPROC, NPROCS. If I am the
*  source process, nothing preceeds me ...
*/
   if( PROC == SRCPROC ) return( 0 );
/*
*  If SRCPROC owns the N rows or columns, then return N since I cannot be the
*  source process anymore.
*/
   if( N <= INB ) return( N );
/*
*  Find out how many full blocks are globally (nblocks) and locally (ilocblk)
*  in those N entries.
*/
   nblocks = ( N - INB ) / NB + 1;
/*
*  Compute my distance from the source process so that within this process
*  coordinate system, the source process is the process such that mydist=0.
*/
   if( ( mydist = PROC - SRCPROC ) < 0 ) mydist += NPROCS;
/*
*  When mydist < nblocks - ilocblk * NPROCS, I own ilocblk + 1 full blocks,
*  when mydist > nblocks - ilocblk * NPROCS, I own ilocblk     full blocks,
*  when mydist = nblocks - ilocblk * NPROCS, either the last block is not full
*  and I own it, or the last block is full and I am the first process owning
*  only ilocblk full blocks.
*
*  Therefore, when 0 < mydist <= nblocks - ilocblk * NPROCS, the number of rows
*  or columns preceeding me is INB + ilocblk*NB + (mydist-1)*(ilocblk+1)*NB,
*  i.e. INB - NB + ( ilocblk+1 ) * NB * mydist. Otherwise, there are exactly
*  NB * ilocblk * ( NPROCS - mydist ) rows or columns after me including mine,
*  i.e N + NB * ilocblk * ( mydist - NPROCS ) rows or columns preceeding me.
*/
   if( nblocks < NPROCS )
      return( ( ( mydist <= nblocks ) ? INB + NB * ( mydist - 1 ) : N ) );

   ilocblk = nblocks / NPROCS;
   return( ( ( mydist <= ( nblocks - ilocblk * NPROCS ) ) ?
             INB - NB + ( ilocblk + 1 ) * NB * mydist :
             N + NB * ilocblk * ( mydist - NPROCS ) ) );
/*
*  End of PB_Cnpreroc
*/
}
