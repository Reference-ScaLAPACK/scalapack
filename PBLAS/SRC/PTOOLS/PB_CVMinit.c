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
void PB_CVMinit( PB_VM_T * VM, Int OFFD, Int M, Int N, Int IMB1, Int INB1,
                 Int MB, Int NB, Int MRROW, Int MRCOL, Int NPROW,
                 Int NPCOL, Int LCMB )
#else
void PB_CVMinit( VM, OFFD, M, N, IMB1, INB1, MB, NB, MRROW, MRCOL, NPROW,
                 NPCOL, LCMB )
/*
*  .. Scalar Arguments ..
*/
   Int            IMB1, INB1, LCMB, M, MB, MRCOL, MRROW, N, NB, NPCOL,
                  NPROW, OFFD;
/*
*  .. Array Arguments ..
*/
   PB_VM_T        * VM;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CVMinit  initializes a virtual matrix with the information of an m
*  by  n  local  array  owned  by  the  process  of relative coordinates
*  ( MRROW, MRCOL ).
*
*  Arguments
*  =========
*
*  VM      (local output) pointer to a PB_VM_T structure
*          On entry,  VM  is  a pointer to a structure of type  PB_VM_T.
*          On exit,  VM  points  to the initialized structure containing
*          the virtual matrix information (see pblas.h).
*
*  OFFD    (global input) INTEGER
*          On entry,  OFFD  specifies the off-diagonal of the underlying
*          matrix of interest as follows:
*             OFFD = 0 specifies the main diagonal,
*             OFFD > 0 specifies lower subdiagonals, and
*             OFFD < 0 specifies upper superdiagonals.
*
*  M       (local input) INTEGER
*          On entry, M  specifies the local number of rows of the under-
*          lying matrix  owned  by the  process  of relative coordinates
*          ( MRROW, MRCOL ). M must be at least zero.
*
*  N       (local input) INTEGER
*          On entry, N  specifies the local number of columns of the un-
*          derlying matrix  owned by the process of relative coordinates
*          ( MRROW, MRCOL ). N must be at least zero.
*
*  IMB1    (global input) INTEGER
*          On input, IMB1 specifies  the global true size of  the  first
*          block of rows of the underlying global submatrix.  IMB1  must
*          be at least MIN( 1, M ).
*
*  INB1    (global input) INTEGER
*          On input, INB1 specifies  the global true size of  the  first
*          block  of  columns  of  the underlying global submatrix. INB1
*          must be at least MIN( 1, N ).
*
*  MB      (global input) INTEGER
*          On entry, MB  specifies the blocking factor used to partition
*          the rows of the matrix.  MB  must be at least one.
*
*  NB      (global input) INTEGER
*          On entry, NB  specifies the blocking factor used to partition
*          the the columns of the matrix.  NB  must be at least one.
*
*  MRROW   (local input) INTEGER
*          On entry, MRROW specifies the  relative row coordinate of the
*          process that possesses these M rows. MRROW must be least zero
*          and strictly less than NPROW.
*
*  MRCOL   (local input) INTEGER
*          On entry, MRCOL specifies  the  relative column coordinate of
*          the process that possesses these N  columns.  MRCOL  must  be
*          least zero and strictly less than NPCOL.
*
*  NPROW   (global input) INTEGER
*          On entry,  NPROW  specifies the total number of process rows
*          over which the matrix is distributed. NPROW must be at least
*          one.
*
*  NPCOL   (global input) INTEGER
*          On entry,  NPCOL  specifies the total number of process col-
*          umns over which the matrix is distributed.  NPCOL must be at
*          least one.
*
*  LCMB    (global input) INTEGER
*          On entry,  LCMB  specifies  the  least  common  multiple  of
*          NPROW * MB and NPCOL * NB.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   Int            tmp1;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Initialize the fields of the VM structure
*/
   VM->offd   = OFFD;
   VM->lcmt00 = OFFD;
   VM->mp     = M;     VM->imb1  = IMB1;  VM->mb = MB; VM->upp = MB - 1;
   VM->prow   = MRROW; VM->nprow = NPROW;
   VM->nq     = N;     VM->inb1   = INB1; VM->nb = NB; VM->low = 1 - NB;
   VM->pcol   = MRCOL; VM->npcol  = NPCOL;
   VM->lcmb   = LCMB;

   if( ( M <= 0 ) || ( N <= 0 ) )
   {
/*
*  If the local virtual array is empty, then simplify the remaining of the
*  initialization.
*/
      VM->imbloc  = 0; VM->lmbloc = 0; VM->mblks  = 0;
      VM->iupp    = ( MRROW ? MB - 1 : ( IMB1 > 0 ? IMB1 - 1 : 0 ) );
      VM->inbloc  = 0; VM->lnbloc = 0; VM->nblks  = 0;
      VM->ilow    = ( MRCOL ? 1 - NB : ( INB1 > 0 ? 1 - INB1 : 0 ) );
      VM->lcmt00 += ( VM->low - VM->ilow + MRCOL * NB ) -
                    ( VM->iupp - VM->upp + MRROW * MB );
      return;
   }

   if( MRROW )
   {
/*
*  I am not in the first relative process row. Use the first local row block
*  size MB to initialize the VM structure.
*/
      VM->lcmt00 -= IMB1 - MB + MRROW * MB;
      VM->imbloc  = MIN( M, MB );
      VM->mblks   = ( M - 1 ) / MB + 1;
      VM->iupp    = MB - 1;
      VM->lmbloc  = M - ( M / MB ) * MB;
      if( !( VM->lmbloc ) ) VM->lmbloc = MB;

      if( MRCOL )
      {
/*
*  I am not in the first relative process column. Use the first local column
*  block size NB to initialize the VM structure.
*/
         VM->inbloc  = MIN( N, NB );
         VM->ilow    = 1 - NB;
         VM->lcmt00 += INB1 - NB + MRCOL * NB;
         VM->nblks   = ( N - 1 ) / NB + 1;
         VM->lnbloc  = N - ( N / NB ) * NB;
         if( !( VM->lnbloc ) ) VM->lnbloc = NB;
      }
      else
      {
/*
*  I am in the first relative process column. Use the first column block size
*  INB1 to initialize the VM structure.
*/
         VM->inbloc = INB1;
         VM->ilow   = 1 - INB1;
         tmp1       = N - INB1;
         if( tmp1 )
         {
/*
*  There is more than one column block. Compute the number of local column
*  blocks and the size of the last one.
*/
            VM->nblks  = ( tmp1 - 1 ) / NB + 2;
            VM->lnbloc = tmp1 - ( tmp1 / NB ) * NB;
            if( !( VM->lnbloc ) ) VM->lnbloc = NB;
         }
         else
         {
/*
*  There is only one column block.
*/
            VM->nblks  = 1;
            VM->lnbloc = INB1;
         }
      }
   }
   else
   {
/*
*  I am in the first relative process row. Use the first row block size IMB1 to
*  initialize the VM structure.
*/
      VM->imbloc = IMB1;
      VM->iupp   = IMB1 - 1;
      tmp1       = M - IMB1;
      if( tmp1 )
      {
/*
*  There is more than one row block. Compute the number of local row blocks and
*  the size of the last one.
*/
         VM->mblks  = ( tmp1 - 1 ) / MB + 2;
         VM->lmbloc = tmp1 - ( tmp1 / MB ) * MB;
         if( !( VM->lmbloc ) ) VM->lmbloc = MB;
      }
      else
      {
/*
*  There is only one row block.
*/
         VM->mblks  = 1;
         VM->lmbloc = IMB1;
      }
      if( MRCOL )
      {
/*
*  I am not in the first relative process column. Use the first local column
*  block size NB to initialize the VM structure.
*/
         VM->inbloc  = MIN( N, NB );
         VM->ilow    = 1 - NB;
         VM->lcmt00 += INB1 - NB + MRCOL * NB;
         VM->nblks   = ( N - 1 ) / NB + 1;
         VM->lnbloc  = N - ( N / NB ) * NB;
         if( !( VM->lnbloc ) ) VM->lnbloc = NB;
      }
      else
      {
/*
*  I am in the first relative process column. Use the first column block size
*  INB1 to initialize the VM structure.
*/
         VM->inbloc = INB1;
         VM->ilow   = 1 - INB1;
         tmp1       = N - INB1;
         if( tmp1 )
         {
/*
*  There is more than one column block. Compute the number of local column
*  blocks and the size of the last one.
*/
            VM->nblks  = ( tmp1 - 1 ) / NB + 2;
            VM->lnbloc = tmp1 - ( tmp1 / NB ) * NB;
            if( !( VM->lnbloc ) ) VM->lnbloc = NB;
         }
         else
         {
/*
*  There is only one column block.
*/
            VM->nblks  = 1;
            VM->lnbloc = INB1;
         }
      }
   }
/*
*  End of PB_CVMinit
*/
}
