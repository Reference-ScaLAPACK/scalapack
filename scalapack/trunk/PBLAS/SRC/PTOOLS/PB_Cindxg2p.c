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
int PB_Cindxg2p( int IG, int INB, int NB, int PROC, int SRCPROC, int NPROCS )
#else
int PB_Cindxg2p( IG, INB, NB, PROC, SRCPROC, NPROCS )
/*
*  .. Scalar Arguments ..
*/
   int            IG, INB, NB, NPROCS, PROC, SRCPROC;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cindxg2p computes the process coordinate  which posseses the entry
*  of a matrix specified by a global index IG.
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
*  PROC    (local dummy) INTEGER
*          On entry, PROC is a dummy argument in this case in  order  to
*          unify the calling sequence of the tool-routines.
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
/* ..
*  .. Executable Statements ..
*
*/
   if( ( IG < INB ) || ( SRCPROC == -1 ) || ( NPROCS == 1 ) )
/*
*  IG belongs to the first block, or the data is not distributed, or there is
*  just one process in this dimension of the grid.
*/
      return( SRCPROC );
/*
*  Otherwise, IG is in block 1 + ( IG - INB ) / NB. Add this to SRCPROC and
*  take the NPROCS modulo (definition of the block-cyclic data distribution).
*/
   PROC = SRCPROC + 1 + ( IG - INB ) / NB;
   return( MPosMod( PROC, NPROCS ) );
/*
*  End of PB_Cindxg2p
*/
}
