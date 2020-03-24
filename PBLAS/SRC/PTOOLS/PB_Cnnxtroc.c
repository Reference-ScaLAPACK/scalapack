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
Int PB_Cnnxtroc( Int N, Int I, Int INB, Int NB, Int PROC, Int SRCPROC,
                 Int NPROCS )
#else
Int PB_Cnnxtroc( N, I, INB, NB, PROC, SRCPROC, NPROCS )
/*
*  .. Scalar Arguments ..
*/
   Int            I, INB, N, NB, NPROCS, PROC, SRCPROC;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cnnxtroc computes  the number of next rows or columns of a subma-
*  trix that are possessed by processes closer to  SRCPROC1  than  PROC
*  where  SRCPROC1 is the process owning the row or column globally in-
*  dexed by I. The submatrix is defined by giving out N rows or columns
*  starting  from  global index I.  Therefore, if SRCPROC=0 and PROC=1,
*  then PB_Cnnxtroc returns the number of matrix rows or columns  owned
*  by processes 2, 3 ... NPROCS-1.
*
*  In fact, if the same exact parameters N,  I,  INB,  NB,  SRCPROC and
*  NPROCS are passed to PB_Cnpreroc, PB_Cnumroc and PB_Cnnxtroc  produ-
*  cing respectively npre, np and nnxt,  then  npre + np + nnxt = N  in
*  every process PROC.
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
   Int            ilocblk, mydist, nblocks;
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
*  Now everything is just like N, I=0, INB, NB, SRCPROC, NPROCS. If the source
*  process owns the N rows or columns, nothing follows me ...
*/
   if( N <= INB ) return( 0 );
/*
*  The discussion goes as follows: compute my distance from the source process
*  so that within this process coordinate system, the source process is the
*  process such that mydist = 0, or equivalently PROC == SRCPROC.
*
*  Find out how many full blocks are globally (nblocks) and locally (ilocblk)
*  in those N entries. Then remark that
*
*  when mydist < nblocks - ilocblk * NPROCS, I own ilocblk + 1 full blocks,
*  when mydist > nblocks - ilocblk * NPROCS, I own ilocblk     full blocks,
*  when mydist = nblocks - ilocblk * NPROCS, either the last block is not full
*  and I own it, or the last block is full and I am the first process owning
*  only ilocblk full blocks.
*/
   nblocks = ( N - INB ) / NB + 1;

   if( PROC == SRCPROC )
   {
/*
*  First note that I cannot be the source and the last process because mydist=0
*  and NPROCS > 1. Since mydist = 0 and nblocks - ilocblk * NPROCS >= 0, there
*  are only two possible cases:
*
*    1) When mydist = nblocks - ilocblk * NPROCS = 0, that is NPROCS divides
*       the global number of full blocks, then the source process SRCPROC owns
*       one more block than the other processes; Thus, N can be rewritten as
*       N  = INB + (nblocks-1) * NB + LNB with LNB >= 0 size of the last block.
*       Similarly, the local value Np corresponding to the local number of rows
*       and columns owned by the source process is INB + (ilocblk-1)*NB + LNB,
*       that is N + ( ilocblk-1 - (nblocks-1) )*NB. Therefore, there must be
*       ( nblocks - ilocblk ) * NB rows or columns following me. Note that this
*       case cannot happen when ilocblk is zero, since nblocks is at least one.
*
*    2) mydist = 0 < nblocks - ilocblk * NPROCS, the source process only owns
*       full blocks, and therefore locally INB + ilocblk * NB rows or columns.
*       Thus, N - INB - ilocblk * NB rows or columns follow me. Note that when
*       ilocblk is zero, this becomes simply N - INB.
*/
      if( nblocks < NPROCS ) return( N - INB );

      ilocblk = nblocks / NPROCS;
      return( ( ( nblocks - ilocblk * NPROCS ) ? N - INB - ilocblk * NB :
                ( nblocks - ilocblk ) * NB ) );
   }
   else
   {
/*
*  I am not the source process. Compute my distance from the source process.
*/
      if( ( mydist = PROC - SRCPROC ) < 0 ) mydist += NPROCS;
/*
*  If I am the last process i.e. mydist = NPROCS - 1, nothing follows me.
*/
      if( mydist == NPROCS - 1 ) return( 0 );
/*
*  Otherwise, when mydist >= nblocks - ilocblk * NPROCS, there are exactly
*  NB * ilocblk * ( NPROCS - mydist ) rows or columns after me including mine,
*  i.e NB * ilocblk * ( NPROCS - 1 - mydist ) rows or columns following me.
*  Finally, when 0 < mydist < nblocks - ilocblk * NPROCS, the number of rows
*  or columns preceeding me is INB + ilocblk * NB + mydist*( ilocblk+1 )*NB
*  including mine, therefore there are N-INB-NB*( ilocblk+mydist*(ilocblk+1) )
*  rows or columns following me.
*/
      if( nblocks < NPROCS )
         return( ( ( mydist < nblocks ) ? N - mydist * NB - INB : 0 ) );

      ilocblk = nblocks / NPROCS;
      return( ( ( mydist >= ( nblocks - ilocblk * NPROCS ) ) ?
                ( NPROCS - 1 - mydist ) * ilocblk * NB :
                N - INB - ( ilocblk * mydist + ilocblk + mydist )*NB ) );
   }
/*
*  End of PB_Cnnxtroc
*/
}
