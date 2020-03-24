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
Int PB_CVMswp( PBTYP_T * TYPE, PB_VM_T * VM, char * VROCS, char * ROCS,
               char * TRANS, Int MN, char * X, Int INCX, char * Y,
               Int INCY )
#else
Int PB_CVMswp( TYPE, VM, VROCS, ROCS, TRANS, MN, X, INCX, Y, INCY )
/*
*  .. Scalar Arguments ..
*/
   Int            INCX, INCY, MN;
/*
*  .. Array Arguments ..
*/
   char           * VROCS, * ROCS, * TRANS;
   PBTYP_T        * TYPE;
   PB_VM_T        * VM;
   char           * X, * Y;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CVMswp  swaps  a one-dimensional distributed vector X with another
*  one-dimensional distributed vector Y. This operation is triggered  by
*  a virtual distributed array.
*
*  Arguments
*  =========
*
*  TYPE    (local input) pointer to a PBTYP_T structure
*          On entry,  TYPE  is a pointer to a structure of type PBTYP_T,
*          that contains type information (See pblas.h).
*
*  VM      (local input) pointer to a PB_VM_T structure
*          On entry,  VM  is  a  pointer to a structure of type PB_VM_T,
*          that contains the virtual matrix information (see pblas.h).
*
*  VROCS   (local input) pointer to CHAR
*          On entry,  VROCS specifies if the rows or columns of the vir-
*          tual distributed array grid should be used  for  the swapping
*          operation as follows:
*             VROCS = 'R' or 'r', the rows should be used,
*             VROCS = 'C' or 'c', the columns should be used.
*
*  ROCS    (local input) pointer to CHAR
*          On entry,  ROCS  specifies if rows or columns should be swap-
*          ped as follows:
*             ROCS = 'R' or 'r', rows should be swapped,
*             ROCS = 'C' or 'c', columns should be swapped.
*
*  TRANS   (local input) pointer to CHAR
*          On entry,  TRANS  specifies if transposition should occur du-
*          ring the  swapping operation as follows:
*             TRANS = 'N' or 'n', natural swapping,
*             otherwise,          transposed swapping.
*
*  MN      (local input) INTEGER
*          On entry,  MN  specifies  the number of rows or columns to be
*          swapped. MN must be at least zero.
*
*  X       (local input/local output) pointer to CHAR
*          On  entry,  X  points  to  an  array  of  dimension  at least
*          ( 1 + ( n - 1 )*abs( INCX ) ) where n is IMBLOC+(MBLKS-2)*MB+
*          LMB when VROCS is 'R'  or  'r',  and  INBLOC+(NBLKS-2)*NB+LNB
*          otherwise. Before entry, the incremented array X must contain
*          the vector x. On exit, the entries of the incremented array X
*          are exchanged with the entries of the incremented array Y.
*
*  INCX    (local input) INTEGER
*          On entry, INCX specifies the increment for the elements of X.
*          INCX must not be zero.
*
*  Y       (local input/local output) pointer to CHAR
*          On  entry,  Y  points  to  an  array  of  dimension  at least
*          ( 1 + ( n - 1 )*abs( INCY ) ) where n is IMBLOC+(MBLKS-2)*MB+
*          LMB when VROCS is 'C'  or  'c',  and  INBLOC+(NBLKS-2)*NB+LNB
*          otherwise. Before entry, the incremented array Y must contain
*          the vector y. On exit, the entries of the incremented array Y
*          are exchanged with the entries of the incremented array X.
*
*  INCY    (local input) INTEGER
*          On entry, INCY specifies the increment for the elements of Y.
*          INCY must not be zero.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   Int            GoEast, GoSouth, Xinc, Yinc, ilow, imbloc, inbloc, iupp, kb,
                  lcmt, lcmt00, lmbloc, lnbloc, low, mb, mblkd, mblks, mbloc,
                  nb, nblkd, nblks, nbloc, notran, npcol, npq=0, nprow, pmb,
                  qnb, rows, size, tmp1, tmp2, upp;
   char           * Xptrd, * Yptrd;
/* ..
*  .. Executable Statements ..
*
*/
   mblks = VM->mblks; nblks = VM->nblks;
/*
*  Quick return if I don't own any blocks.
*/
   if( ( mblks == 0 ) || ( nblks == 0 ) ) return( 0 );
/*
*  Retrieve the contents of VM structure fields
*/
   lcmt00 = VM->lcmt00;
   imbloc = VM->imbloc; mb    = VM->mb; lmbloc = VM->lmbloc; upp = VM->upp;
   iupp   = VM->iupp;   nprow = VM->nprow;
   inbloc = VM->inbloc; nb    = VM->nb; lnbloc = VM->lnbloc; low = VM->low;
   ilow   = VM->ilow;   npcol = VM->npcol;

   notran = ( Mupcase( TRANS[0] ) == CNOTRAN );

   size   = TYPE->size;
   rows   = ( Mupcase( ROCS[0] ) == CROW );

   if( Mupcase( VROCS[0] ) == CROW )
   {
/*
*  (un)packing using rows of virtual matrix
*/
      if( rows )
      {
/*
*  (un)packing rows of mn by k array A.
*/
         Xinc = size;
         Yinc = ( notran ? size : INCY * size );
      }
      else
      {
/*
*  (un)packing columns of k by mn array A
*/
         Xinc = INCX * size;
         Yinc = ( notran ? INCY * size : size );
      }
      kb  = MN;
/*
*  From the (un)packing point of view the only valuable shortcut is when the
*  virtual grid and the blocks are square, and the offset is zero or the grid
*  is 1x1.
*/
      if( ( ( lcmt00 == 0 ) && ( VM->imb1 == VM->inb1 ) && ( mb == nb ) &&
            ( nprow == npcol ) ) || ( ( nprow == 1 ) && ( npcol == 1 ) ) )
      {
         if(  VM->prow == VM->pcol )
         {
            npq = ( ( mblks <  2 ) ? imbloc :
                    imbloc + ( mblks - 2 ) * mb + lmbloc );
            npq = MIN( npq, kb );
            if( rows ) TYPE->Fswap( &npq, X, &INCX, Y, &INCY );
            else       TYPE->Fswap( &npq, X, &INCX, Y, &INCY );
         }
         return( npq );
      }
      pmb = nprow * mb;
      qnb = npcol * nb;
/*
*  Handle separately the first row and/or column of the LCM table. Update the
*  LCM value of the curent block lcmt00, as well as the number of rows and
*  columns mblks and nblks remaining in the LCM table.
*/
      GoSouth = ( lcmt00 > iupp );
      GoEast  = ( lcmt00 < ilow );

      if( !( GoSouth ) && !( GoEast ) )
      {
/*
*  The upper left block owns diagonal entries lcmt00 >= ilow && lcmt00 <= iupp
*/
         if( lcmt00 >= 0 )
         {
            tmp1 = imbloc - lcmt00; tmp1 = MAX( 0, tmp1 );
            tmp2 = MIN( tmp1, inbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
            TYPE->Fswap( &tmp2, X+lcmt00*Xinc, &INCX, Y, &INCY );
         }
         else
         {
            tmp1 = inbloc + lcmt00; tmp1 = MAX( 0, tmp1 );
            tmp2 = MIN( tmp1, imbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
            TYPE->Fswap( &tmp2, X, &INCX, Y-lcmt00*Yinc, &INCY );
         }
         if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Decide whether one should go south or east in the table: Go east if
*  the block below the current one only owns lower entries. If this block,
*  however, owns diagonals, then go south.
*/
         GoSouth = !( GoEast = ( ( lcmt00 - ( iupp - upp + pmb ) ) < ilow ) );
      }

      if( GoSouth )
      {
/*
*  Go one step south in the LCM table. Adjust the current LCM value as well as
*  the pointer to X. The pointer to Y remains unchanged.
*/
         lcmt00 -= iupp - upp + pmb; mblks--; X += imbloc * Xinc;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value as well as the pointer to X accordingly.
*/
         while( mblks && ( lcmt00 > upp ) )
         { lcmt00 -= pmb; mblks--; X += mb * Xinc; }
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
         lcmt = lcmt00; mblkd = mblks; Xptrd = X;

         while( mblkd && ( lcmt >= ilow ) )
         {
/*
*  A block owning diagonals lcmt00 >= ilow && lcmt00 <= upp has been found.
*/
            mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
            if( lcmt >= 0 )
            {
               tmp1 = mbloc - lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, inbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               TYPE->Fswap( &tmp2, Xptrd+lcmt*Xinc, &INCX, Y, &INCY );
            }
            else
            {
               tmp1 = inbloc + lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, mbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               TYPE->Fswap( &tmp2, Xptrd, &INCX, Y-lcmt*Yinc, &INCY );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Keep going south until there are no more blocks owning diagonals
*/
            lcmt -= pmb; mblkd--; Xptrd += mbloc * Xinc;
         }
/*
*  I am done with the first column of the LCM table. Go to the next column.
*/
         lcmt00 += low - ilow + qnb; nblks--; Y += inbloc * Yinc;
      }
      else if( GoEast )
      {
/*
*  Go one step east in the LCM table. Adjust the current LCM value as
*  well as the pointer to Y. The pointer to X remains unchanged.
*/
         lcmt00 += low - ilow + qnb; nblks--; Y += inbloc * Yinc;
/*
*  While there are blocks remaining that own lower entries, keep going east
*  in the LCM table. Adjust the current LCM value as well as the pointer to
*  Y accordingly.
*/
         while( nblks && ( lcmt00 < low ) )
         { lcmt00 += qnb; nblks--; Y += nb * Yinc; }
/*
*  Return if no more column in the LCM table.
*/
         if( nblks <= 0 ) return( npq );
/*
*  lcmt00 >= low. The current block owns either diagonals or upper entries. Save
*  the current position in the LCM table. After this row has been completely
*  taken care of, re-start from this column and the next row of the LCM table.
*/
         lcmt = lcmt00; nblkd = nblks; Yptrd = Y;

         while( nblkd && ( lcmt <= iupp ) )
         {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= iupp has been found.
*/
            nbloc = ( ( nblkd == 1 ) ? lnbloc : nb );
            if( lcmt >= 0 )
            {
               tmp1 = imbloc - lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, nbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               TYPE->Fswap( &tmp2, X+lcmt*Xinc, &INCX, Yptrd, &INCY );
            }
            else
            {
               tmp1 = nbloc + lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, imbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               TYPE->Fswap( &tmp2, X, &INCX, Yptrd-lcmt*Yinc, &INCY );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Keep going east until there are no more blocks owning diagonals.
*/
            lcmt += qnb; nblkd--; Yptrd += nbloc * Yinc;
         }
/*
*  I am done with the first row of the LCM table. Go to the next row.
*/
         lcmt00 -= iupp - upp + pmb; mblks--; X += imbloc * Xinc;
      }
/*
*  Loop over the remaining columns of the LCM table.
*/
      do
      {
/*
*  If the current block does not have diagonal elements, find the closest one in
*  the LCM table having some.
*/
         if( ( lcmt00 < low ) || ( lcmt00 > upp ) )
         {
            while( mblks && nblks )
            {
               while( mblks && ( lcmt00 > upp ) )
               { lcmt00 -= pmb; mblks--; X += mb * Xinc; }
               if( lcmt00 >= low ) break;
               while( nblks && ( lcmt00 < low ) )
               { lcmt00 += qnb; nblks--; Y += nb * Yinc; }
               if( lcmt00 <= upp ) break;
            }
         }
         if( !mblks || !nblks ) return( npq );
/*
*  The current block owns diagonals. Save the current position in the LCM table.
*  After this column has been completely taken care of, re-start from this row
*  and the next column in the LCM table.
*/
         nbloc = ( ( nblks == 1 ) ? lnbloc : nb );
         lcmt = lcmt00; mblkd = mblks; Xptrd = X;

         while( mblkd && lcmt >= low )
         {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
*/
            mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
            if( lcmt >= 0 )
            {
               tmp1 = mbloc - lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, nbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               TYPE->Fswap( &tmp2, Xptrd+lcmt*Xinc, &INCX, Y, &INCY );
            }
            else
            {
               tmp1 = nbloc + lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, mbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               TYPE->Fswap( &tmp2, Xptrd, &INCX, Y-lcmt*Yinc, &INCY );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Keep going south until there are no more blocks owning diagonals
*/
            lcmt -= pmb; mblkd--; Xptrd += mbloc * Xinc;
         }
/*
*  I am done with this column of the LCM table. Go to the next column ...
*/
         lcmt00 += qnb; nblks--; Y += nbloc * Yinc;
/*
*  ... until there are no more columns.
*/
      } while( nblks > 0 );
/*
*  Return the number of diagonals found.
*/
      return( npq );
   }
   else
   {
/*
*  (un)packing using columns of virtual matrix
*/
      if( rows )
      {
/*
*  (un)packing rows of mn by k array A
*/
         Xinc = size;
         Yinc = ( notran ? size : INCY * size );
      }
      else
      {
/*
*  (un)packing columns of k by mn array A
*/
         Xinc = INCX * size;
         Yinc = ( notran ? INCY * size : size );
      }
      kb = MN;
/*
*  From the (un)packing point of view the only valuable shortcut is when the
*  virtual grid and the blocks are square, and the offset is zero or the grid
*  is 1x1.
*/
      if( ( ( lcmt00 == 0 ) && ( VM->imb1 == VM->inb1 ) && ( mb == nb ) &&
            ( nprow == npcol ) ) || ( ( nprow == 1 ) && ( npcol == 1 ) ) )
      {
         if(  VM->prow == VM->pcol )
         {
            npq = ( ( nblks <  2 ) ? inbloc :
                    inbloc + ( nblks - 2 ) * nb + lnbloc );
            npq = MIN( npq, kb );
            if( rows ) TYPE->Fswap( &npq, X, &INCX, Y, &INCY );
            else       TYPE->Fswap( &npq, X, &INCX, Y, &INCY );
         }
         return( npq );
      }
      pmb = nprow * mb;
      qnb = npcol * nb;
/*
*  Handle separately the first row and/or column of the LCM table. Update the
*  LCM value of the curent block lcmt00, as well as the number of rows and
*  columns mblks and nblks remaining in the LCM table.
*/
      GoSouth = ( lcmt00 > iupp );
      GoEast  = ( lcmt00 < ilow );

      if( !( GoSouth ) && !( GoEast ) )
      {
/*
*  The upper left block owns diagonal entries lcmt00 >= ilow && lcmt00 <= iupp
*/
         if( lcmt00 >= 0 )
         {
            tmp1 = imbloc - lcmt00; tmp1 = MAX( 0, tmp1 );
            tmp2 = MIN( tmp1, inbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
            TYPE->Fswap( &tmp2, X, &INCX, Y+lcmt00*Yinc, &INCY );
         }
         else
         {
            tmp1 = inbloc + lcmt00; tmp1 = MAX( 0, tmp1 );
            tmp2 = MIN( tmp1, imbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
            TYPE->Fswap( &tmp2, X-lcmt00*Xinc, &INCX, Y, &INCY );
         }
         if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Decide whether one should go south or east in the table: Go east if
*  the block below the current one only owns lower entries. If this block,
*  however, owns diagonals, then go south.
*/
         GoSouth = !( GoEast = ( ( lcmt00 - ( iupp - upp + pmb ) ) < ilow ) );
      }

      if( GoSouth )
      {
/*
*  Go one step south in the LCM table. Adjust the current LCM value as well as
*  the pointer to Y. The pointer to X remains unchanged.
*/
         lcmt00 -= iupp - upp + pmb; mblks--; Y += imbloc * Yinc;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value as well as the pointer to Y accordingly.
*/
         while( mblks && ( lcmt00 > upp ) )
         { lcmt00 -= pmb; mblks--; Y += mb * Yinc; }
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
         lcmt  = lcmt00; mblkd = mblks; Yptrd = Y;

         while( mblkd && ( lcmt >= ilow ) )
         {
/*
*  A block owning diagonals lcmt00 >= ilow && lcmt00 <= upp has been found.
*/
            mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
            if( lcmt >= 0 )
            {
               tmp1 = mbloc - lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, inbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               TYPE->Fswap( &tmp2, X, &INCX, Yptrd+lcmt*Yinc, &INCY );
            }
            else
            {
               tmp1 = inbloc + lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, mbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               TYPE->Fswap( &tmp2, X-lcmt*Xinc, &INCX, Yptrd, &INCY );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Keep going south until there are no more blocks owning diagonals
*/
            lcmt -= pmb; mblkd--; Yptrd += mbloc * Yinc;
         }
/*
*  I am done with the first column of the LCM table. Go to the next column.
*/
         lcmt00 += low - ilow + qnb; nblks--; X += inbloc * Xinc;
      }
      else if( GoEast )
      {
/*
*  Go one step east in the LCM table. Adjust the current LCM value as
*  well as the pointer to X. The pointer to Y remains unchanged.
*/
         lcmt00 += low - ilow + qnb; nblks--; X += inbloc * Xinc;
/*
*  While there are blocks remaining that own lower entries, keep going east
*  in the LCM table. Adjust the current LCM value as well as the pointer to
*  X accordingly.
*/
         while( nblks && ( lcmt00 < low ) )
         { lcmt00 += qnb; nblks--; X += nb * Xinc; }
/*
*  Return if no more column in the LCM table.
*/
         if( nblks <= 0 ) return( npq );
/*
*  lcmt00 >= low. The current block owns either diagonals or upper entries. Save
*  the current position in the LCM table. After this row has been completely
*  taken care of, re-start from this column and the next row of the LCM table.
*/
         lcmt  = lcmt00; nblkd = nblks; Xptrd = X;

         while( nblkd && ( lcmt <= iupp ) )
         {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= iupp has been found.
*/
            nbloc = ( ( nblkd == 1 ) ? lnbloc : nb );
            if( lcmt >= 0 )
            {
               tmp1 = imbloc - lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, nbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               TYPE->Fswap( &tmp2, Xptrd, &INCX, Y+lcmt*Yinc, &INCY );
            }
            else
            {
               tmp1 = nbloc + lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, imbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               TYPE->Fswap( &tmp2, Xptrd-lcmt*Xinc, &INCX, Y, &INCY );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Keep going east until there are no more blocks owning diagonals.
*/
            lcmt += qnb; nblkd--; Xptrd += nbloc * Xinc;
         }
/*
*  I am done with the first row of the LCM table. Go to the next row.
*/
         lcmt00 -= iupp - upp + pmb; mblks--; Y += imbloc * Yinc;
      }
/*
*  Loop over the remaining columns of the LCM table.
*/
      do
      {
/*
*  If the current block does not have diagonal elements, find the closest one in
*  the LCM table having some.
*/
         if( ( lcmt00 < low ) || ( lcmt00 > upp ) )
         {
            while( mblks && nblks )
            {
               while( mblks && ( lcmt00 > upp ) )
               { lcmt00 -= pmb; mblks--; Y += mb * Yinc; }
               if( lcmt00 >= low ) break;
               while( nblks && ( lcmt00 < low ) )
               { lcmt00 += qnb; nblks--; X += nb * Xinc; }
               if( lcmt00 <= upp ) break;
            }
         }
         if( !( mblks ) || !( nblks ) ) return( npq );
/*
*  The current block owns diagonals. Save the current position in the LCM table.
*  After this column has been completely taken care of, re-start from this row
*  and the next column in the LCM table.
*/
         nbloc = ( ( nblks == 1 ) ? lnbloc : nb );
         lcmt = lcmt00; mblkd = mblks; Yptrd = Y;
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
*/
         while( mblkd && lcmt >= low )
         {
            mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
            if( lcmt >= 0 )
            {
               tmp1 = mbloc - lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, nbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               TYPE->Fswap( &tmp2, X, &INCX, Yptrd+lcmt*Yinc, &INCY );
            }
            else
            {
               tmp1 = nbloc + lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, mbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               TYPE->Fswap( &tmp2, X-lcmt*Xinc, &INCX, Yptrd, &INCY );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Keep going south until there are no more blocks owning diagonals
*/
            lcmt -= pmb; mblkd--; Yptrd += mbloc * Yinc;
         }
/*
*  I am done with this column of the LCM table. Go to the next column ...
*/
         lcmt00 += qnb; nblks--; X += nbloc * Xinc;
/*
*  ... until there are no more columns.
*/
      } while( nblks > 0 );
/*
*  Return the number of diagonals found.
*/
      return( npq );
   }
/*
*  End of PB_CVMswp
*/
}
