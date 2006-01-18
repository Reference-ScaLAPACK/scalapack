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
void PB_Cinfog2l( int I, int J, int * DESC, int NPROW, int NPCOL,
                  int MYROW, int MYCOL, int * II, int * JJ,
                  int * PROW, int * PCOL )
#else
void PB_Cinfog2l( I, J, DESC, NPROW, NPCOL, MYROW, MYCOL, II, JJ,
                  PROW, PCOL )
   int            I, * II, J, * JJ, MYCOL, MYROW, NPCOL, NPROW, * PCOL,
                  * PROW;
/*
*  .. Scalar Arguments ..
*/
/*
*  .. Array Arguments ..
*/
   int            * DESC;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cinfog2l computes the starting local index II, JJ corresponding to
*  the submatrix starting globally at the entry pointed by  I,  J.  This
*  routine returns the coordinates in the grid of the process owning the
*  matrix entry of global indexes I, J, namely PROW and PCOL.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information  required to
*  establish the  mapping  between a  matrix entry and its corresponding
*  process and memory location.
*
*  In  the  following  comments,   the character _  should  be  read  as
*  "of  the  distributed  matrix".  Let  A  be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESC_A:
*
*  NOTATION         STORED IN       EXPLANATION
*  ---------------- --------------- ------------------------------------
*  DTYPE_A (global) DESCA[ DTYPE_ ] The descriptor type.
*  CTXT_A  (global) DESCA[ CTXT_  ] The BLACS context handle, indicating
*                                   the NPROW x NPCOL BLACS process grid
*                                   A  is  distributed over. The context
*                                   itself  is  global,  but  the handle
*                                   (the integer value) may vary.
*  M_A     (global) DESCA[ M_     ] The  number of rows in the distribu-
*                                   ted matrix A, M_A >= 0.
*  N_A     (global) DESCA[ N_     ] The number of columns in the distri-
*                                   buted matrix A, N_A >= 0.
*  IMB_A   (global) DESCA[ IMB_   ] The number of rows of the upper left
*                                   block of the matrix A, IMB_A > 0.
*  INB_A   (global) DESCA[ INB_   ] The  number  of columns of the upper
*                                   left   block   of   the  matrix   A,
*                                   INB_A > 0.
*  MB_A    (global) DESCA[ MB_    ] The blocking factor used to  distri-
*                                   bute the last  M_A-IMB_A  rows of A,
*                                   MB_A > 0.
*  NB_A    (global) DESCA[ NB_    ] The blocking factor used to  distri-
*                                   bute the last  N_A-INB_A  columns of
*                                   A, NB_A > 0.
*  RSRC_A  (global) DESCA[ RSRC_  ] The process row over which the first
*                                   row of the matrix  A is distributed,
*                                   NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA[ CSRC_  ] The  process column  over  which the
*                                   first column of  A  is  distributed.
*                                   NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA[ LLD_   ] The  leading dimension  of the local
*                                   array  storing  the  local blocks of
*                                   the distributed matrix A,
*                                   IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                   ELSE
*                                      LLD_A >= 1.
*
*  Let K be the number of  rows of a matrix A starting at the global in-
*  dex IA,i.e, A( IA:IA+K-1, : ). Lr( IA, K ) denotes the number of rows
*  that the process of row coordinate MYROW ( 0 <= MYROW < NPROW ) would
*  receive if these K rows were distributed over NPROW processes.  If  K
*  is the number of columns of a matrix  A  starting at the global index
*  JA, i.e, A( :, JA:JA+K-1, : ), Lc( JA, K ) denotes the number  of co-
*  lumns that the process MYCOL ( 0 <= MYCOL < NPCOL ) would  receive if
*  these K columns were distributed over NPCOL processes.
*
*  The values of Lr() and Lc() may be determined via a call to the func-
*  tion PB_Cnumroc:
*  Lr( IA, K ) = PB_Cnumroc( K, IA, IMB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_Cnumroc( K, JA, INB_A, NB_A, MYCOL, CSRC_A, NPCOL )
*
*  Arguments
*  =========
*
*  I       (global input) INTEGER
*          On entry, I  specifies  the  global starting row index of the
*          submatrix. I must at least zero.
*
*  J       (global input) INTEGER
*          On entry, J  specifies  the  global  starting column index of
*          the submatrix. J must at least zero.
*
*  DESC    (global and local input) INTEGER array
*          On entry,  DESC is an integer array of dimension DLEN_.  This
*          is the array descriptor of the underlying matrix.
*
*  NPROW   (global input) INTEGER
*          On entry,  NPROW   specifies the total number of process rows
*          over which the matrix is distributed.  NPROW must be at least
*          one.
*
*  NPCOL   (global input) INTEGER
*          On entry, NPCOL specifies the total number of process columns
*          over which the matrix is distributed.  NPCOL must be at least
*          one.
*
*  MYROW   (local input) INTEGER
*          On entry,  MYROW  specifies the row coordinate of the process
*          whose local index  II  is determined.  MYROW must be at least
*          zero and strictly less than NPROW.
*
*  MYCOL   (local input) INTEGER
*          On entry,  MYCOL  specifies the column coordinate of the pro-
*          cess whose local index  JJ  is determined.  MYCOL  must be at
*          least zero and strictly less than NPCOL.
*
*  II      (local output) INTEGER
*          On exit, II  specifies the  local  starting  row index of the
*          submatrix. On exit, II is at least zero.
*
*  JJ      (local output) INTEGER
*          On exit, JJ  specifies the local starting column index of the
*          submatrix. On exit, JJ is at least zero.
*
*  PROW    (global output) INTEGER
*          On exit,  PROW  specifies  the  row coordinate of the process
*          that possesses the first row of the submatrix.  On exit, PROW
*          is -1 if DESC( RSRC_ )  is -1 on input, and,  at  least  zero
*          and strictly less than NPROW otherwise.
*
*  PCOL    (global output) INTEGER
*          On exit, PCOL  specifies the column coordinate of the process
*          that possesses the first column of the  submatrix.  On  exit,
*          PCOL is -1 if DESC( CSRC_ )  is -1 on input, and,  at   least
*          zero and strictly less than NPCOL otherwise.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   int            ilocblk, imb, inb, mb, mydist, nb, nblocks, csrc, rsrc;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Retrieve the row distribution parameters
*/
   imb   = DESC[IMB_ ];
   *PROW = DESC[RSRC_];

   if( ( *PROW == -1 ) || ( NPROW == 1 ) )
   {
/*
*  The data is not distributed, or there is just one process row in the grid.
*/
     *II = I;
   }
   else if( I < imb )
   {
/*
*  I refers to an entry in the first block of rows
*/
     *II = ( MYROW == *PROW ? I : 0 );
   }
   else
   {
      mb   = DESC[MB_];
      rsrc = *PROW;
/*
*  The discussion goes as follows: compute my distance from the source process
*  so that within this process coordinate system, the source process is the
*  process such that mydist = 0, or equivalently MYROW == rsrc.
*
*  Find out the global coordinate of the block I belongs to (nblocks), as well
*  as the minimum local number of blocks that every process has.
*
*  when mydist < nblocks - ilocblk * NPROCS, I own ilocblk + 1 full blocks,
*  when mydist > nblocks - ilocblk * NPROCS, I own ilocblk     full blocks,
*  when mydist = nblocks - ilocblk * NPROCS, I own ilocblk     full blocks
*  but not I, or I own ilocblk + 1 blocks and the entry I refers to.
*/
      if( MYROW == rsrc )
      {
/*
*  I refers to an entry that is not in the first block, find out which process
*  has it.
*/
         nblocks = ( I - imb ) / mb + 1;
         *PROW  += nblocks;
         *PROW  -= ( *PROW / NPROW ) * NPROW;
/*
*  Since mydist = 0 and nblocks - ilocblk * NPROW >= 0, there are only three
*  possible cases:
*
*    1) When 0 = mydist = nblocks - ilocblk * NPROW = 0 and I don't own I, in
*       which case II = IMB + ( ilocblk - 1 ) * MB. Note that this case cannot
*       happen when ilocblk is zero, since nblocks is at least one.
*
*    2) When 0 = mydist = nblocks - ilocblk * NPROW = 0 and I own I, in which
*       case I and II can respectively be written as IMB + (nblocks-1)*NB + IL
*       and IMB + (ilocblk-1) * MB + IL. That is II = I + (ilocblk-nblocks)*MB.
*       Note that this case cannot happen when ilocblk is zero, since nblocks
*       is at least one.
*
*    3) mydist = 0 < nblocks - ilocblk * NPROW, the source process owns
*       ilocblk+1 full blocks, and therefore II = IMB + ilocblk * MB. Note
*       that when ilocblk is zero, II is just IMB.
*/
         if( nblocks < NPROW )
         {
            *II = imb;
         }
         else
         {
            ilocblk = nblocks / NPROW;
            if( ilocblk * NPROW >= nblocks )
            {
               *II = ( ( MYROW == *PROW ) ? I + ( ilocblk - nblocks ) * mb :
                       imb + ( ilocblk - 1 ) * mb );
            }
            else
            {
               *II =  imb + ilocblk * mb;
            }
         }
      }
      else
      {
/*
*  I refers to an entry that is not in the first block, find out which process
*  has it.
*/
         nblocks = ( I -= imb ) / mb + 1;
         *PROW  += nblocks;
         *PROW  -= ( *PROW / NPROW ) * NPROW;
/*
*  Compute my distance from the source process so that within this process
*  coordinate system, the source process is the process such that mydist=0.
*/
         if( ( mydist  = MYROW - rsrc ) < 0 ) mydist += NPROW;
/*
*  When mydist <  nblocks - ilocblk * NPROW, I own ilocblk + 1 full blocks of
*  size MB since I am not the source process, i.e. II = ( ilocblk + 1 ) * MB.
*  When mydist >= nblocks - ilocblk * NPROW and I don't own I, I own ilocblk
*  full blocks of size MB, i.e. II = ilocblk * MB, otherwise I own ilocblk
*  blocks and I, in which case I can be written as IMB + (nblocks-1)*MB + IL
*  and II = ilocblk*MB + IL = I - IMB + ( ilocblk - nblocks + 1 )*MB.
*/
         if( nblocks < NPROW )
         {
            mydist -= nblocks;
            *II     = ( ( mydist < 0 ) ? mb :
                        ( ( MYROW == *PROW ) ? I + ( 1 - nblocks ) * mb : 0 ) );
         }
         else
         {
            ilocblk = nblocks / NPROW;
            mydist -= nblocks - ilocblk * NPROW;
            *II     = ( ( mydist < 0 ) ? ( ilocblk + 1 ) * mb :
                        ( ( MYROW == *PROW ) ?
                          ( ilocblk - nblocks + 1 ) * mb + I : ilocblk * mb ) );
         }
      }
   }
/*
*  Idem for the columns
*/
   inb   = DESC[INB_ ];
   *PCOL = DESC[CSRC_];

   if( ( *PCOL == -1 ) || ( NPCOL == 1 ) )
   {
      *JJ = J;
   }
   else if( J < inb )
   {
      *JJ = ( MYCOL == *PCOL ? J : 0 );
   }
   else
   {
      nb   = DESC[NB_];
      csrc = *PCOL;

      if( MYCOL == csrc )
      {
         nblocks = ( J - inb ) / nb + 1;
         *PCOL  += nblocks;
         *PCOL  -= ( *PCOL / NPCOL ) * NPCOL;

         if( nblocks < NPCOL )
         {
            *JJ = inb;
         }
         else
         {
            ilocblk = nblocks / NPCOL;
            if( ilocblk * NPCOL >= nblocks )
            {
               *JJ = ( ( MYCOL == *PCOL ) ? J + ( ilocblk - nblocks ) * nb :
                       inb + ( ilocblk - 1 ) * nb );
            }
            else
            {
               *JJ = inb + ilocblk * nb;
            }
         }
      }
      else
      {
         nblocks = ( J -= inb ) / nb + 1;
         *PCOL  += nblocks;
         *PCOL  -= ( *PCOL / NPCOL ) * NPCOL;

         if( ( mydist = MYCOL - csrc ) < 0 ) mydist += NPCOL;

         if( nblocks < NPCOL )
         {
            mydist -= nblocks;
            *JJ     = ( ( mydist < 0 ) ? nb : ( ( MYCOL == *PCOL ) ?
                        J + ( 1 - nblocks )*nb : 0 ) );
         }
         else
         {
            ilocblk = nblocks / NPCOL;
            mydist -= nblocks - ilocblk * NPCOL;
            *JJ     = ( ( mydist < 0 ) ? ( ilocblk + 1 ) * nb :
                        ( ( MYCOL == *PCOL ) ?
                          ( ilocblk - nblocks + 1 ) * nb + J : ilocblk * nb ) );
         }
      }
   }
/*
*  End of PB_Cinfog2l
*/
}
