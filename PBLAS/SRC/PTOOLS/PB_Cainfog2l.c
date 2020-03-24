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
void PB_Cainfog2l( Int M, Int N, Int I, Int J, Int * DESC, Int NPROW,
                   Int NPCOL, Int MYROW, Int MYCOL, Int * IMB1,
                   Int * INB1, Int * MP, Int * NQ, Int * II, Int * JJ,
                   Int * PROW, Int * PCOL, Int * RPROW, Int * RPCOL )
#else
void PB_Cainfog2l( M, N, I, J, DESC, NPROW, NPCOL, MYROW, MYCOL, IMB1,
                   INB1, MP, NQ, II, JJ, PROW, PCOL, RPROW, RPCOL )
/*
*  .. Scalar Arguments ..
*/
   Int            I, * II, * IMB1, * INB1, J, * JJ, M, * MP, MYCOL,
                  MYROW, N, NPCOL, NPROW, * NQ, * PCOL, * PROW, * RPCOL,
                  * RPROW;
/*
*  .. Array Arguments ..
*/
   Int            * DESC;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cainfog2l computes the  starting  local row and column indexes II,
*  JJ  corresponding to  the  submatrix  starting  globally at the entry
*  pointed by I,  J. This routine returns the coordinates in the grid of
*  the  process owning  the  matrix entry of global indexes I, J, namely
*  PROW  and  PCOL. In addition, this routine computes the quantities MP
*  and  NQ,  which are respectively the local number of rows and columns
*  owned by the process of coordinate  MYROW, MYCOL corresponding to the
*  global submatrix A(I:I+M-1,J:J+N-1).  Finally, the size  of the first
*  partial block and the relative process coordinates  are also returned
*  respectively in IMB, INB and RPROW, RPCOL.
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
*  M       (global input) INTEGER
*          On entry, M specifies the global number of rows of the subma-
*          trix. M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N specifies  the  global  number  of columns of the
*          submatrix. N must be at least zero.
*
*  I       (global input) INTEGER
*          On entry, I  specifies  the  global starting row index of the
*          submatrix. I must at least zero.
*
*  J       (global input) INTEGER
*          On entry, J  specifies  the global starting column  index  of
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
*  IMB1    (global output) INTEGER
*          On exit, IMB1 specifies the number of rows of the upper  left
*          block of the submatrix. On exit,  IMB1 is less or equal  than
*          M and greater or equal than MIN( 1, M ).
*
*  INB1    (global output) INTEGER
*          On exit, INB1 specifies  the number  of  columns of the upper
*          left block of the submatrix. On exit,  INB1 is  less or equal
*          than N and greater or equal than MIN( 1, N ).
*
*  MP      (local output) INTEGER
*          On exit, MP specifies the local number of rows of the  subma-
*          trix, that the processes of row coordinate MYROW own.  MP  is
*          at least zero.
*
*  NQ      (local output) INTEGER
*          On exit, NQ specifies  the  local  number  of columns  of the
*          submatrix,  that  the processes  of column  coordinate  MYCOL
*          own. NQ is at least zero.
*
*  II      (local output) INTEGER
*          On exit, II  specifies the  local  starting  row index of the
*          submatrix. On exit, II is at least zero.
*
*  JJ      (local output) INTEGER
*          On exit, JJ  specifies the  local  starting  column index  of
*          the submatrix. On exit, II is at least zero.
*
*  PROW    (global output) INTEGER
*          On exit,  PROW  specifies the row coordinate of  the  process
*          that possesses the first row of the submatrix. On exit,  PROW
*          is -1 if DESC(RSRC_)  is -1 on input, and, at least zero  and
*          strictly less than NPROW otherwise.
*
*  PCOL    (global output) INTEGER
*          On exit, PCOL  specifies the column coordinate of the process
*          that possesses the first column of the  submatrix.  On  exit,
*          PCOL is -1 if DESC(CSRC_)  is -1 on input, and, at least zero
*          and strictly less than NPCOL otherwise.
*
*  RPROW   (global output) INTEGER
*          On exit, RPROW specifies  the  relative row coordinate of the
*          process that possesses the first row  I  of the submatrix. On
*          exit, RPROW is -1 if DESC(RSRC_) is  -1  on  input,  and,  at
*          least zero and strictly less than NPROW otherwise.
*
*  RPCOL   (global output) INTEGER
*          On exit, RPCOL specifies  the  relative column  coordinate of
*          the process that possesses the first column  J  of the subma-
*          trix. On exit, RPCOL is -1 if  DESC(CSRC_)  is  -1  on input,
*          and, at least zero and strictly less than NPCOL otherwise.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   Int            i1, ilocblk, j1, mb, mydist, nb, nblocks, csrc, rsrc;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Retrieve the row distribution parameters
*/
   mb   = DESC[ MB_   ];
   rsrc = DESC[ RSRC_ ];

   if( ( rsrc == -1 ) || ( NPROW == 1 ) )
   {
/*
*  The rows are not distributed, or there is just one process row in the grid.
*  Therefore, the local and global indexes are the same, as well as the local
*  and global number of rows. Finally, the relative row process coordinate is
*  zero, since every process owns all rows. Note that the size of the first
*  row block can be zero only if M is zero.
*/
      *II    = I;
      if( ( *IMB1 = DESC[IMB_] - I ) <= 0 )
         *IMB1 += ( ( -(*IMB1) ) / mb + 1 ) * mb;
      *IMB1  = MIN( *IMB1, M );
      *MP    = M;
      *PROW  = rsrc;
      *RPROW = 0;
   }
   else
   {
/*
*  Figure out PROW, II and IMB1 first.
*/
      *IMB1 = DESC[IMB_];
      if( I < *IMB1 )                  /* Is I in first block range ? */
      {
/*
*  If I is in the first block of rows, then PROW is simply rsrc, II is I in
*  this process and zero elsewhere, and the size of the first block is the
*  IMB complement.
*/
         *PROW  = rsrc;
         *II    = ( ( MYROW == *PROW ) ? I : 0 );
         *IMB1 -= I;
      }
      else
      {
/*
*  The discussion goes as follows: compute my distance from the source process
*  so that within this process coordinate system, the source row process is the
*  process such that mydist=0, or equivalently MYROW == rsrc.
*
*  Find out the global coordinate of the block of rows I belongs to (nblocks),
*  as well as the minimum local number of row blocks that every process has.
*
*  when mydist < nblocks - ilocblk * NPROW, I own ilocblk + 1 full blocks,
*  when mydist > nblocks - ilocblk * NPROW, I own ilocblk     full blocks,
*  when mydist = nblocks - ilocblk * NPROW, I own ilocblk     full blocks
*  but not I, or I own ilocblk + 1 blocks and the entry I refers to.
*/
         i1 = I - *IMB1;
         if( MYROW == rsrc )
         {
/*
*  I refers to an entry that is not in the first block, find out which process
*  has it.
*/
            nblocks = i1 / mb + 1;
            *PROW   = rsrc + nblocks;
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
*       case I and II can respectively be written as IMB + (nblocks-1)*MB + IL
*       and IMB+(ilocblk-1) * MB + IL. That is II = I + (ilocblk - nblocks)*MB.
*       Note that this case cannot happen when ilocblk is zero, since nblocks
*       is at least one.
*
*    3) mydist = 0 < nblocks - ilocblk * NPROW, the source process owns
*       ilocblk+1 full blocks, and therefore II = IMB + ilocblk * MB. Note
*       that when ilocblk is zero, II is just IMB.
*/
            if( nblocks < NPROW )
            {
               *II = *IMB1;
            }
            else
            {
               ilocblk = nblocks / NPROW;
               if( ilocblk * NPROW >= nblocks )
               {
                  *II = ( ( MYROW == *PROW ) ? I + ( ilocblk - nblocks ) * mb :
                          *IMB1 + ( ilocblk - 1 ) * mb );
               }
               else
               {
                  *II = *IMB1 + ilocblk * mb;
               }
            }
         }
         else
         {
/*
*  I is not in the first block, find out which process has it.
*/
            nblocks = i1 / mb + 1;
            *PROW   = rsrc + nblocks;
            *PROW  -= ( *PROW / NPROW ) * NPROW;
/*
*  Compute my distance from the source process so that within this process
*  coordinate system, the source process is the process such that mydist=0.
*/
            if( ( mydist = MYROW - rsrc ) < 0 ) mydist += NPROW;
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
               *II     = ( ( mydist < 0 ) ? mb : ( ( MYROW == *PROW ) ?
                           i1 + ( 1 - nblocks ) * mb : 0 ) );
            }
            else
            {
               ilocblk = nblocks / NPROW;
               mydist -= nblocks - ilocblk * NPROW;
               *II     = ( ( mydist < 0 ) ? ( ilocblk + 1 ) * mb :
                             ( ( MYROW == *PROW ) ?
                               ( ilocblk - nblocks + 1 ) * mb + i1 :
                               ilocblk * mb ) );
            }
         }
/*
*  Update the size of first block
*/
         *IMB1 = nblocks * mb - i1;
      }
/*
*  Now everything is just like M, I=0, IMB1, MB, PROW, NPROW. The discussion
*  goes as follows: compute my distance from the source process PROW so that
*  within this process coordinate system, the source process is the process
*  such that mydist = 0. Figure out MP.
*/
      if( M <= *IMB1 )
      {
/*
*  M <= IMB1: if I am the source process, i.e. I own I (mydist = 0), MP is M
*  and 0 otherwise.
*/
         *MP = ( ( MYROW == *PROW ) ? M : 0 );
      }
      else
      {
/*
*  Find out how many full blocks are globally (nblocks) and locally (ilocblk)
*  in those M entries
*/
         nblocks = ( M - *IMB1 ) / mb + 1;

         if( MYROW == *PROW )
         {
/*
*  Since mydist = 0 and nblocks - ilocblk * NPROW >= 0, there are only two
*  possible cases:
*
*    1) When mydist = nblocks - ilocblk * NPROW = 0, that is NPROW divides
*       the global number of full blocks, then the source process PROW owns
*       one more block than the other processes; and M can be rewritten as
*       M  = IMB1 + (nblocks-1) * NB + LNB with LNB >= 0 size of the last block.
*       Similarly, the local value MP corresponding to M can be written as
*       MP = IMB1 + (ilocblk-1) * MB + LMB = M + ( ilocblk-1 - (nblocks-1) )*MB.
*       Note that this case cannot happen when ilocblk is zero, since nblocks
*       is at least one.
*
*    2) mydist = 0 < nblocks - ilocblk * NPROW, the source process only owns
*       full blocks, and therefore MP = IMB1 + ilocblk * MB. Note that when
*       ilocblk is zero, MP is just IMB1.
*/
            if( nblocks < NPROW )
            {
               *MP = *IMB1;
            }
            else
            {
               ilocblk = nblocks / NPROW;
               *MP     = ( ( nblocks - ilocblk * NPROW ) ?
                           *IMB1 + ilocblk * mb :
                           M + ( ilocblk - nblocks ) * mb );
            }
         }
         else
         {
/*
*  Compute my distance from the source process so that within this process
*  coordinate system, the source process is the process such that mydist=0.
*/
            if( ( mydist = MYROW - *PROW ) < 0 ) mydist += NPROW;
/*
*  When mydist < nblocks - ilocblk * NPROW, I own ilocblk + 1 full blocks of
*  size MB since I am not the source process,
*
*  when mydist > nblocks - ilocblk * NPROW, I own ilocblk     full blocks of
*  size MB since I am not the source process,
*
*  when mydist = nblocks - ilocblk * NPROW,
*     either the last block is not full and I own it, in which case
*        M = IMB1 + (nblocks - 1)*MB + LMB with LNB the size of the last block
*        such that MB > LMB > 0; the local value MP corresponding to M is given
*        by MP = ilocblk * MB + LMB = M - IMB1 + ( ilocblk - nblocks + 1 ) * MB;
*     or the last block is full and I am the first process owning only ilocblk
*        full blocks of size MB, that is M = IMB + ( nblocks - 1 ) * MB and
*        MP = ilocblk * MB = M - IMB + ( ilocblk - nblocks + 1 ) * MB.
*/
            if( nblocks < NPROW )
            {
               mydist -= nblocks;
               *MP     = ( ( mydist < 0 ) ? mb : ( ( mydist > 0 ) ? 0 :
                           M - *IMB1 + mb * ( 1 - nblocks ) ) );
            }
            else
            {
               ilocblk = nblocks / NPROW;
               mydist -= nblocks - ilocblk * NPROW;
               *MP     = ( ( mydist < 0 ) ? ( ilocblk + 1 ) * mb :
                           ( ( mydist > 0 ) ? ilocblk * mb :
                           M - *IMB1 + mb * ( ilocblk - nblocks + 1 ) ) );
            }
         }
      }
/*
*  Finally figure out IMB1 and RPROW. Note that IMB1 can be zero when M = 0.
*/
      *IMB1  = MIN( *IMB1, M );
      if( ( *RPROW = MYROW - *PROW ) < 0 ) *RPROW += NPROW;
   }
/*
*  Idem for the columns
*/
   nb   = DESC[ NB_   ];
   csrc = DESC[ CSRC_ ];

   if( ( csrc == -1 ) || ( NPCOL == 1 ) )
   {
      *JJ    = J;
      if( ( *INB1 = DESC[INB_] - J ) <= 0 )
         *INB1 += ( ( -(*INB1) ) / nb + 1 ) * nb;
      *INB1  = MIN( *INB1, N );
      *NQ    = N;
      *PCOL  = csrc;
      *RPCOL = 0;
   }
   else
   {
      *INB1 = DESC[INB_];
      if( J < *INB1 )
      {
         *PCOL  = csrc;
         *JJ    = ( ( MYCOL == *PCOL ) ? J : 0 );
         *INB1 -= J;
      }
      else
      {
         j1 = J - *INB1;
         if( MYCOL == csrc )
         {
            nblocks = j1 / nb + 1;
            *PCOL   = csrc + nblocks;
            *PCOL  -= ( *PCOL / NPCOL ) * NPCOL;

            if( nblocks < NPCOL )
            {
               *JJ = *INB1;
            }
            else
            {
               ilocblk = nblocks / NPCOL;
               if( ilocblk * NPCOL >= nblocks )
               {
                  *JJ = ( ( MYCOL == *PCOL ) ? J + ( ilocblk - nblocks ) * nb :
                          *INB1 + ( ilocblk - 1 ) * nb );
               }
               else
               {
                  *JJ = *INB1 + ilocblk * nb;
               }
            }
         }
         else
         {
            nblocks = j1 / nb + 1;
            *PCOL   = csrc + nblocks;
            *PCOL  -= ( *PCOL / NPCOL ) * NPCOL;

            if( ( mydist  = MYCOL - csrc ) < 0 ) mydist += NPCOL;

            if( nblocks < NPCOL )
            {
               mydist -= nblocks;
               *JJ     = ( ( mydist < 0 ) ? nb : ( ( MYCOL == *PCOL ) ?
                           j1 + ( 1 - nblocks ) * nb : 0 ) );
            }
            else
            {
               ilocblk = nblocks / NPCOL;
               mydist -= nblocks - ilocblk * NPCOL;
               *JJ     = ( ( mydist < 0 ) ? ( ilocblk + 1 ) * nb :
                           ( ( MYCOL == *PCOL ) ?
                             ( ilocblk - nblocks + 1 ) * nb + j1 :
                             ilocblk * nb ) );
            }
         }
         *INB1 = nblocks * nb - j1;
      }

      if( N <= *INB1 )
      {
         *NQ = ( ( MYCOL == *PCOL ) ? N : 0 );
      }
      else
      {
         nblocks = ( N - *INB1 ) / nb + 1;

         if( MYCOL == *PCOL )
         {
            if( nblocks < NPCOL )
            {
               *NQ = *INB1;
            }
            else
            {
               ilocblk = nblocks / NPCOL;
               *NQ     = ( ( nblocks - ilocblk * NPCOL ) ?
                           *INB1 + ilocblk * nb :
                           N + ( ilocblk - nblocks ) * nb );
            }
         }
         else
         {
            if( ( mydist  = MYCOL - *PCOL ) < 0 ) mydist += NPCOL;

            if( nblocks < NPCOL )
            {
               mydist -= nblocks;
               *NQ     = ( ( mydist < 0 ) ? nb : ( ( mydist > 0 ) ? 0 :
                           N - *INB1 + nb * ( 1 - nblocks ) ) );
            }
            else
            {
               ilocblk = nblocks / NPCOL;
               mydist -= nblocks - ilocblk * NPCOL;
               *NQ     = ( ( mydist < 0 ) ? ( ilocblk + 1 ) * nb :
                           ( ( mydist > 0 ) ? ilocblk * nb :
                           N - *INB1 + nb * ( ilocblk - nblocks + 1 ) ) );
            }
         }
      }
      *INB1  = MIN( *INB1, N );
      if( ( *RPCOL = MYCOL - *PCOL ) < 0 ) *RPCOL += NPCOL;
   }
/*
*  End of PB_Cainfog2l
*/
}
