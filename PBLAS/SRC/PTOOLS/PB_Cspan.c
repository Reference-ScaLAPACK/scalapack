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
int PB_Cspan( int N, int I, int INB, int NB, int SRCPROC, int NPROCS )
#else
int PB_Cspan( N, I, INB, NB, SRCPROC, NPROCS )
/*
*  .. Scalar Arguments ..
*/
   int            I, INB, N, NB, NPROCS, SRCPROC;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cspan  returns  1  if  the rows (resp. columns) I:I+N-1 spans more
*  than one process row (resp. column) and 0 otherwise.
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
*  SRCPROC (global input) INTEGER
*          On entry, if SRCPROC = -1, the data  is not  distributed  but
*          replicated,  in  which  case  this  routine return 0.
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
/*
*  If the data is replicated or if there is only one process in this dimension
*  of the process grid, the data does not span multiple processes. Finally, the
*  case where I belongs to the first block is handled separately.
*/
   return( ( SRCPROC >= 0 ) &&
             ( ( NPROCS > 1 ) &&
               ( ( I < INB ) ? ( I + N > INB ) :
                 ( I + N > INB + ( ( I - INB ) / NB + 1 ) * NB ) ) ) );
/*
*  End of PB_Cspan
*/
}
