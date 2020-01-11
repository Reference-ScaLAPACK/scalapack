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
Int PB_Cg2lrem( Int IG, Int INB, Int NB, Int MYPROC, Int SRCPROC, Int NPROCS )
#else
Int PB_Cg2lrem( IG, INB, NB, MYPROC, SRCPROC, NPROCS )
/*
*  .. Scalar Arguments ..
*/
   Int            IG, INB, NB, NPROCS, MYPROC, SRCPROC;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cg2lrem computes the local index of a matrix entry pointed  to  by
*  the global index IG. Note that when MYPROC is not the process  owning
*  this entry, this routine  returns the closest larger local index cor-
*  responding to IG just like the routine PB_Cinfog2l.
*
*  Arguments
*  =========
*
*  IG      (global input) INTEGER
*          On entry, IG specifies the global index of the matrix  entry.
*          IG must be at least zero.
*
*  INB     (global input) INTEGER
*          On entry,  INB  specifies  the size of the first block of the
*          global matrix. INB must be at least one.
*
*  NB      (global input) INTEGER
*          On entry, NB specifies the size of the blocks used to  parti-
*          tion the matrix. NB must be at least one.
*
*  MYPROC  (global input) INTEGER
*          On entry,  MYPROC  specifies  the process number in which the
*          value of the local index is to be computed. MYPROC must be at
*          least zero and strictly less than NPROCS.
*
*  SRCPROC (global input) INTEGER
*          On entry, if SRCPROC = -1, the data  is not  distributed  but
*          replicated,  in  which  case  this  routine returns IG in all
*          processes. Otherwise, the value of SRCPROC is ignored.
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
   Int            ilocblk, mydist, nblocks, proc;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  The data is not distributed, or there is just one process in this dimension
*  of the grid.
*/
   if( ( SRCPROC == -1 ) || ( NPROCS == 1 ) ) return( IG );
/*
*  IG refers to an entry in the first block
*/
   if( IG < INB ) return( ( MYPROC == SRCPROC ? IG : 0 ) );
/*
*  The discussion goes as follows: compute my distance from the source process
*  so that within this process coordinate system, the source process is the
*  process such that mydist = 0, or equivalently MYROC == SRCPROC.
*
*  Find out the global coordinate of the block IG belongs to (nblocks), as well
*  as the minimum local number of blocks that every process has.
*
*  when mydist < nblocks - ilocblk * NPROCS, I own ilocblk + 1 full blocks,
*  when mydist > nblocks - ilocblk * NPROCS, I own ilocblk     full blocks,
*  when mydist = nblocks - ilocblk * NPROCS, I own ilocblk     full blocks
*  but not IG, or I own ilocblk + 1 blocks and the entry IG refers to.
*/
   if( MYPROC == SRCPROC )
   {
/*
*  If I am the source process and there are less than NPROCS blocks, then
*  the local index in that process is INB.
*/
      nblocks = ( IG - INB ) / NB + 1;
      if( nblocks < NPROCS ) return( INB );
/*
*  IG refers to an entry that is not in the first block, find out which process
*  has it.
*/
      proc    = SRCPROC + nblocks;
      proc   -= ( proc / NPROCS ) * NPROCS;
/*
*  Since mydist = 0 and nblocks - ilocblk * NPROCS >= 0, there are only three
*  possible cases:
*
*    1) When 0 = mydist = nblocks - ilocblk * NPROCS = 0 and I don't own IG, in
*       which case II = INB + ( ilocblk - 1 ) * NB. Note that this case cannot
*       happen when ilocblk is zero, since nblocks is at least one.
*
*    2) When 0 = mydist = nblocks - ilocblk * NPROCS = 0 and I own IG, in which
*       case IG and II can respectively be written as INB + (nblocks-1)*NB + IL,
*       INB + (ilocblk-1) * NB + IL. That is II = IG + ( ilocblk - nblocks )*NB.
*       Note that this case cannot happen when ilocblk is zero, since nblocks
*       is at least one.
*
*    3) mydist = 0 < nblocks - ilocblk * NPROCS, the source process owns
*       ilocblk+1 full blocks, and therefore II = INB + ilocblk * NB. Note
*       that when ilocblk is zero, II is just INB.
*/
      ilocblk = nblocks / NPROCS;
      if( ilocblk * NPROCS >= nblocks )
         return( ( ( MYPROC == proc ) ? IG + ( ilocblk - nblocks ) * NB :
                   INB + ( ilocblk - 1 ) * NB ) );
      else
         return( INB + ilocblk * NB );
   }
   else
   {
/*
*  IG refers to an entry that is not in the first block, find out which process
*  has it.
*/
      nblocks = ( IG -= INB ) / NB + 1;
      proc    = SRCPROC + nblocks;
      proc   -= ( proc / NPROCS ) * NPROCS;
/*
*  Compute my distance from the source process so that within this process
*  coordinate system, the source process is the process such that mydist=0.
*/
      if( ( mydist = MYPROC - SRCPROC ) < 0 ) mydist += NPROCS;
/*
*  When mydist <  nblocks - ilocblk * NPROCS, I own ilocblk + 1 full blocks of
*  size NB since I am not the source process, i.e. II = ( ilocblk + 1 ) * NB.
*  When mydist >= nblocks - ilocblk * NPROCS and I don't own IG, I own ilocblk
*  full blocks of size NB, i.e. II = ilocblk * NB, otherwise I own ilocblk
*  blocks and IG, in which case IG can be written as INB + (nblocks-1)*NB + IL
*  and II = ilocblk*NB + IL = IG - INB + ( ilocblk - nblocks + 1 )*NB.
*/
      if( nblocks < NPROCS )
      {
         mydist -= nblocks;
         return( ( ( mydist < 0 ) ? NB :
                   ( ( MYPROC == proc ) ? IG + ( 1 - nblocks ) * NB : 0 ) ) );
      }
      else
      {
         ilocblk = nblocks / NPROCS;
         mydist -= nblocks - ilocblk * NPROCS;
         return( ( ( mydist < 0 ) ? ( ilocblk + 1 ) * NB :
                   ( ( MYPROC == proc ) ?
                   ( ilocblk - nblocks + 1 ) * NB + IG : ilocblk * NB ) ) );
      }
   }
/*
*  End of PB_Cg2lrem
*/
}
