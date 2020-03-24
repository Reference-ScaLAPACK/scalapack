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
Int PB_CVMnpq( PB_VM_T * VM )
#else
Int PB_CVMnpq( VM )
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
*  PB_CVMnpq  computes the number of diagonal entries in the virtual ma-
*  specified by VM.
*
*  Arguments
*  =========
*
*  VM      (local input) pointer to a PB_VM_T structure
*          On entry,  VM  is  a pointer to a structure of type  PB_VM_T,
*          that contains the virtual matrix information (see pblas.h).
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   Int            GoEast, GoSouth, Pmb, Qnb, gcdb, ilow, imbloc, inbloc, iupp,
                  kmax, kmin, k1, k2, k3, lcmb, lcmp, lcmq, lcmt, lcmt00,
                  lmbloc, lnbloc, low, l1, l2, l3, m, mb, mblkd, mblks, mbloc,
                  n, nb, nblkd, nblks, nbloc, nlcmblks, npcol, npq=0, nprow,
                  tmp1, tmp2, upp;
/* ..
*  .. Executable Statements ..
*
*/
   m = VM->mp; n = VM->nq;
/*
*  Quick return if I don't own any data.
*/
   if( ( m == 0 ) || ( n == 0 ) ) return( 0 );
/*
*  The only valuable shortcut is when the virtual grid and the blocks are
*  square, and the offset is zero or the grid is 1x1.
*/
   mb = VM->mb; nprow = VM->nprow;
   nb = VM->nb; npcol = VM->npcol;

   if( ( ( VM->offd == 0 ) && ( VM->imb1 == VM->inb1 ) && ( mb == nb ) &&
       ( nprow == npcol ) ) || ( ( nprow == 1 ) && ( npcol == 1 ) ) )
   {
      if( VM->prow == VM->pcol ) return( MIN( m, n ) );
      else                       return( 0 );
   }
/*
*  Retrieve the contents of VM structure fields
*/
   lcmt00 = VM->lcmt00;
   mblks  = VM->mblks; imbloc = VM->imbloc; lmbloc = VM->lmbloc;
   iupp   = VM->iupp;  upp    = VM->upp;
   nblks  = VM->nblks; inbloc = VM->inbloc; lnbloc = VM->lnbloc;
   ilow   = VM->ilow;  low    = VM->low;
   lcmb   = VM->lcmb;
   Pmb    = nprow * mb;
   Qnb    = npcol * nb;
/*
*  Handle separately the first row and/or column of the LCM table. Update the
*  LCM value of the curent block lcmt00, as well as the number of rows and
*  columns mblks and nblks remaining in the LCM table.
*/
   GoSouth  = ( lcmt00 > iupp );
   GoEast   = ( lcmt00 < ilow );
/*
*  Go through the table looking for blocks owning diagonal entries.
*/
   if( !( GoSouth ) && !( GoEast ) )
   {
/*
*  The upper left block owns diagonal entries lcmt00 >= ilow && lcmt00 <= iupp
*/
      npq += ( lcmt00 >= 0 ?
          ( ( tmp2 = ( tmp1 = imbloc - lcmt00 ) > 0 ? tmp1 : 0 ) < inbloc ?
          tmp2 : inbloc ) :
          ( ( tmp2 = ( tmp1 = inbloc + lcmt00 ) > 0 ? tmp1 : 0 ) > imbloc ?
          imbloc : tmp2 ) );
/*
*  Decide whether one should go south or east in the table: Go east if
*  the block below the current one only owns lower entries. If this block,
*  however, owns diagonals, then go south.
*/
      GoSouth = !( GoEast = ( ( lcmt00 - ( iupp - upp + Pmb ) ) < ilow ) );
   }

   if( GoSouth )
   {
/*
*  Go one step south in the LCM table. Adjust the current LCM value.
*/
      lcmt00 -= iupp - upp + Pmb; mblks--;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value accordingly.
*/
      while( mblks && ( lcmt00 > upp ) ) { lcmt00 -= Pmb; mblks--; }
/*
*  Return if no more row in the LCM table.
*/
      if( mblks <= 0 ) return( npq );
/*
*  lcmt00 <= upp. The current block owns either diagonals or lower entries.
*  Save the current position in the LCM table. After this column has been
*  completely taken care of, re-start from this row and the next column of
*  the LCM table.
*/
      lcmt = lcmt00; mblkd = mblks; mbloc = mb;

      while( mblkd && ( lcmt >= ilow ) )
      {
/*
*  A block owning diagonals lcmt00 >= ilow && lcmt00 <= upp has been found.
*/
         if( mblkd == 1 ) mbloc = lmbloc;
         npq += ( lcmt >= 0 ?
             ( ( tmp2 = ( tmp1 = mbloc - lcmt ) > 0 ? tmp1 : 0 ) < inbloc ?
             tmp2 : inbloc ) :
             ( ( tmp2 = ( tmp1 = inbloc + lcmt ) > 0 ? tmp1 : 0 ) > mbloc ?
             mbloc : tmp2 ) );
/*
*  Keep going south until there are no more blocks owning diagonals
*/
         lcmt00 = lcmt; lcmt -= Pmb; mblks = mblkd--;
      }
/*
*  I am done with the first column of the LCM table. Go to the next column.
*/
      lcmt00 += low - ilow + Qnb; nblks--;
   }
   else if( GoEast )
   {
/*
*  Go one step east in the LCM table. Adjust the current LCM value.
*/
      lcmt00 += low - ilow + Qnb; nblks--;
/*
*  While there are blocks remaining that own lower entries, keep going east
*  in the LCM table. Adjust the current LCM value accordingly.
*/
      while( nblks && ( lcmt00 < low ) ) { lcmt00 += Qnb; nblks--; }
/*
*  Return if no more column in the LCM table.
*/
      if( nblks <= 0 ) return( npq );
/*
*  lcmt00 >= low. The current block owns either diagonals or upper entries. Save
*  the current position in the LCM table. After this row has been completely
*  taken care of, re-start from this column and the next row of the LCM table.
*/
      lcmt = lcmt00; nblkd = nblks; nbloc = nb;

      while( nblkd && ( lcmt <= iupp ) )
      {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= iupp has been found.
*/
         if( nblkd == 1 ) nbloc = lnbloc;
         npq += ( lcmt >= 0 ?
             ( ( tmp2 = ( tmp1 = imbloc - lcmt ) > 0 ? tmp1 : 0 ) < nbloc ?
             tmp2 : nbloc ) :
             ( ( tmp2 = ( tmp1 = nbloc + lcmt ) > 0 ? tmp1 : 0 ) > imbloc ?
             imbloc : tmp2 ) );
/*
*  Keep going east until there are no more blocks owning diagonals.
*/
         lcmt00 = lcmt; lcmt += Qnb; nblks = nblkd--;
      }
/*
*  I am done with the first row of the LCM table. Go to the next row.
*/
      lcmt00 -= iupp - upp + Pmb; mblks--;
   }
/*
*  If the current block does not have diagonal elements, find the closest one in
*  the LCM table having some.
*/
   if( lcmt00 < low || lcmt00 > upp )
   {
      while( mblks && nblks )
      {
         while( mblks && ( lcmt00 > upp ) ) { lcmt00 -= Pmb; mblks--; }
         if( lcmt00 >= low ) break;
         while( nblks && ( lcmt00 < low ) ) { lcmt00 += Qnb; nblks--; }
         if( lcmt00 <= upp ) break;
      }
   }
   if( !( mblks ) || !( nblks ) ) return( npq );
/*
*  Figure out how many "full" lcm blocks are remaining.
*/
   gcdb = ( Pmb * Qnb ) / lcmb;

   if( lcmt00 > 0 )
   {
      kmin = - ( lcmb / gcdb );
      kmax = ( lcmb - Qnb ) / gcdb;
      tmp1 = ( mblks - 1 ) / ( lcmp = lcmb / Pmb );
      tmp2 = nblks / ( lcmq = lcmb / Qnb );
   }
   else if( lcmt00 < 0 )
   {
      kmin = - ( ( lcmb - Pmb ) / gcdb );
      kmax = lcmb / gcdb;
      tmp1 = mblks / ( lcmp = lcmb / Pmb );
      tmp2 = ( nblks - 1 ) / ( lcmq = lcmb / Qnb );
   }
   else
   {
      kmin = - ( ( lcmb - Pmb ) / gcdb );
      kmax = ( lcmb - Qnb ) / gcdb;
      tmp1 = mblks / ( lcmp = lcmb / Pmb );
      tmp2 = nblks / ( lcmq = lcmb / Qnb );
   }
/*
*  The last block, even if it is an lcm block will be handled separately
*/
   nlcmblks = MIN( tmp1, tmp2 );
   if( nlcmblks ) nlcmblks--;
/*
*  Compute the lcm block part, update mblks and nblks
*/
   if( nlcmblks )
   {
      tmp2 = 0;

      k1 = -lcmt00; k1 = ICEIL( k1, gcdb ); l1 = k1 - 1;
      l1 = MIN( l1, kmax ); k1 = MAX( k1, kmin );

      k3 = upp - lcmt00; k3 = FLOOR( k3, gcdb ); k3 = MIN( k3, kmax );
      l3 = low - lcmt00; l3 = ICEIL( l3, gcdb ); l3 = MAX( l3, kmin );

      if( k1 <= k3 )
      {
         k2 = mb - nb - lcmt00;
         k2 = ICEIL( k2, gcdb );

         if( k2 < k1 )
         {
/*
*  k2 < k1
*/
            tmp1  = k3 - k1 + 1;
            tmp2  = tmp1 * ( mb - lcmt00 );
            tmp1 *= ( k3 + k1 )*gcdb;
            tmp2 += ( tmp1 > 0 ? -( tmp1 / 2 ) : (-tmp1) / 2 );
         }
         else if( k2 > k3 )
         {
/*
*  k2 = k3 + 1
*/
            tmp2  = ( k3 - k1 + 1 ) * nb;
         }
         else
         {
/*
*  k1 <= k2 <= k3
*/
            tmp1  = k3 - k2 + 1;
            tmp2  = ( k2 - k1 ) * nb + tmp1 * ( mb - lcmt00 );
            tmp1 *= ( k3 + k2 ) * gcdb;
            tmp2 += ( tmp1 > 0 ? -( tmp1 / 2 ) : (-tmp1) / 2 );
         }
      }

      if( l3 <= l1 )
      {
         l2 = mb - nb - lcmt00;
         l2 = FLOOR( l2, gcdb );

         if( l2 > l1 )
         {
/*
*  l2 > l1
*/
            tmp1  = l1 - l3 + 1;
            tmp2 += tmp1 * ( nb + lcmt00 );
            tmp1 *= ( l3 + l1 ) * gcdb;
            tmp2 += ( tmp1 > 0 ? ( tmp1 / 2 ) : -( (-tmp1) / 2 ) );
         }
         else if( l2 < l3 )
         {
/*
*  l2 = l3 - 1
*/
            tmp2 += ( l1 - l3 + 1 ) * mb;
         }
         else
         {
/*
*  l3 <= l2 <= l1
*/
            tmp1  = l2 - l3 + 1;
            tmp2 += ( l1 - l2 ) * mb + tmp1 * ( nb + lcmt00 );
            tmp1 *= ( l3 + l2 ) * gcdb;
            tmp2 += ( tmp1 > 0 ? ( tmp1 / 2 ) : -( (-tmp1) / 2 ) );
         }
      }
      npq   += nlcmblks * tmp2;

      mblks -= nlcmblks * lcmp;
      nblks -= nlcmblks * lcmq;
   }
/*
*  Handle last partial (lcm) block separately
*/
   nbloc = nb;
   while( nblks )
   {
/*
*  The current block owns diagonals. Save the current position in the LCM table.
*  After this column has been completely taken care of, re-start from this row
*  and the next column in the LCM table.
*/
      if( nblks == 1 ) nbloc = lnbloc;
      while( mblks && lcmt00 > upp ) { lcmt00 -= Pmb; mblks--; }

      if( mblks <= 0 ) return( npq );

      lcmt = lcmt00; mblkd = mblks; mbloc = mb;

      while( mblkd && ( lcmt >= low ) )
      {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
*/
         if( mblkd == 1 ) mbloc = lmbloc;
         npq += ( lcmt >= 0 ?
             ( ( tmp2 = ( tmp1 = mbloc - lcmt ) > 0 ? tmp1 : 0 ) < nbloc ?
             tmp2 : nbloc ) :
             ( ( tmp2 = ( tmp1 = nbloc + lcmt ) > 0 ? tmp1 : 0 ) > mbloc ?
             mbloc : tmp2 ) );
/*
*  Keep going south until there are no more blocks owning diagonals
*/
         lcmt00 = lcmt; lcmt -= Pmb; mblks = mblkd--;
      }
/*
*  I am done with this column of the LCM table. Go to the next column ...
*/
      lcmt00 += Qnb; nblks--;
/*
*  ... until there are no more columns.
*/
   }
/*
*  Return the number of diagonals found.
*/
   return( npq );
/*
*  End of PB_CVMnpq
*/
}
