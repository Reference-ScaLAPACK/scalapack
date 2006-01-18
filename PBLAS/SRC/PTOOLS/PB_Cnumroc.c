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
int PB_Cnumroc( int N, int I, int INB, int NB, int PROC, int SRCPROC,
                int NPROCS )
#else
int PB_Cnumroc( N, I, INB, NB, PROC, SRCPROC, NPROCS )
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
*  PB_Cnumroc  returns  the  local number of matrix rows/columns process
*  PROC will get  if we give out N rows/columns starting from global in-
*  dex I.
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
*          global matrix. INB must be at least one.
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
      return( N );
/*
*  Compute coordinate of process owning I and corresponding INB
*/
   if( ( INB -= I ) <= 0 )
   {
/*
*  I is not in the first block, find out which process has it and update the
*  size of first block
*/
      nblocks  = (-INB) / NB + 1;
      SRCPROC += nblocks;
      SRCPROC -= ( SRCPROC / NPROCS ) * NPROCS;
      INB     += nblocks * NB;
   }
/*
*  Now everything is just like N, I=0, INB, NB, SRCPROC, NPROCS. The discussion
*  goes as follows: compute my distance from the source process so that within
*  this process coordinate system, the source process is the process such that
*  mydist = 0, or equivalently PROC == SRCPROC.
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
   if( PROC == SRCPROC )
   {
/*
*  I am the source process, i.e. I own I (mydist = 0). When N <= INB, the
*  answer is simply N.
*/
      if( N <= INB ) return( N );
/*
*  Find out how many full blocks are globally (nblocks) and locally (ilocblk)
*  in those N entries.
*/
      nblocks = ( N - INB ) / NB + 1;
/*
*  Since mydist = 0 and nblocks - ilocblk * NPROCS >= 0, there are only two
*  possible cases:
*
*    1) When mydist = nblocks - ilocblk * NPROCS = 0, that is NPROCS divides
*       the global number of full blocks, then the source process SRCPROC owns
*       one more block than the other processes; and N can be rewritten as
*       N  = INB + (nblocks-1) * NB + LNB with LNB >= 0 size of the last block.
*       Similarly, the local value Np corresponding to N can be written as
*       Np = INB + (ilocblk-1) * NB + LNB = N + ( ilocblk-1 - (nblocks-1) )*NB.
*       Note that this case cannot happen when ilocblk is zero, since nblocks
*       is at least one.
*
*    2) mydist = 0 < nblocks - ilocblk * NPROCS, the source process only owns
*       full blocks, and therefore Np = INB + ilocblk * NB. Note that when
*       ilocblk is zero, Np is just INB.
*/
      if( nblocks < NPROCS ) return( INB );

      ilocblk = nblocks / NPROCS;
      return( ( nblocks - ilocblk * NPROCS ) ? INB + ilocblk * NB :
              N + ( ilocblk - nblocks ) * NB );
   }
   else
   {
/*
*  I am not the source process. When N <= INB, the answer is simply 0.
*/
      if( N <= INB ) return( 0 );
/*
*  Find out how many full blocks are globally (nblocks) and locally (ilocblk)
*  in those N entries
*/
      nblocks = ( N - INB ) / NB + 1;
/*
*  Compute my distance from the source process so that within this process
*  coordinate system, the source process is the process such that mydist=0.
*/
      if( ( mydist = PROC - SRCPROC ) < 0 ) mydist += NPROCS;
/*
*  When mydist < nblocks - ilocblk * NPROCS, I own ilocblk + 1 full blocks of
*  size NB since I am not the source process,
*
*  when mydist > nblocks - ilocblk * NPROCS, I own ilocblk     full blocks of
*  size NB since I am not the source process,
*
*  when mydist = nblocks - ilocblk * NPROCS,
*     either the last block is not full and I own it, in which case
*        N = INB + (nblocks - 1)*NB + LNB with LNB the size of the last block
*        such that NB > LNB > 0; the local value Np corresponding to N is given
*        by Np = ilocblk * NB + LNB = N - INB + ( ilocblk - nblocks + 1 ) * NB;
*     or the last block is full and I am the first process owning only ilocblk
*        full blocks of size NB, that is N = INB + ( nblocks - 1 ) * NB and
*        Np = ilocblk * NB = N - INB + ( ilocblk - nblocks + 1 ) * NB.
*/
      if( nblocks < NPROCS )
         return( ( mydist < nblocks ) ? NB : ( ( mydist > nblocks ) ? 0 :
                 N - INB + NB * ( 1 - nblocks ) ) );

      ilocblk = nblocks / NPROCS;
      mydist -= nblocks - ilocblk * NPROCS;
      return( ( mydist < 0 ) ? ( ilocblk + 1 ) * NB :
              ( ( mydist > 0 ) ? ilocblk * NB :
                N - INB + NB * ( ilocblk - nblocks + 1 ) ) );
   }
/*
*  End of PB_Cnumroc
*/
}
