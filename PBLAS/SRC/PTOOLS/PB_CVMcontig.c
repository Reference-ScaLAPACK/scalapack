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
void PB_CVMcontig( PB_VM_T * VM, int * NRPQ, int * NCPQ, int * IOFF,
                   int * JOFF )
#else
void PB_CVMcontig( VM, NRPQ, NCPQ, IOFF, JOFF )
/*
*  .. Scalar Arguments ..
*/
   int            * IOFF, * JOFF, * NCPQ, * NRPQ;
   PB_VM_T        * VM;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CVMcontig  computes  the  maximum  number  of  contiguous rows and
*  columns  corresponding  to  the  first diagonals of the local virtual
*  matrix VM. This routine also returns the row and column offset of the
*  first diagonal entry.
*
*  Arguments
*  =========
*
*  VM      (local input) pointer to a PB_VM_T structure
*          On entry,  VM  is  a  pointer to a structure of type PB_VM_T,
*          that contains the virtual matrix information (see pblas.h).
*
*  NRPQ    (local output) INTEGER
*          On exit, NRPQ specifies the number of contiguous rows corres-
*          ponding  to  the  first diagonals of the local virtual matrix
*          VM. On exit, NRPQ is at least zero.
*
*  NCPQ    (local output) INTEGER
*          On exit, NCPQ specifies the number of contiguous columns cor-
*          responding to the first diagonals of the local virtual matrix
*          VM. On exit, NRPQ is at least zero.
*
*  IOFF    (local output) INTEGER
*          On exit, IOFF is  the  local row offset of the first row cor-
*          responding to a  diagonal  entry of the Virtual matrix VM. If
*          no diagonals are found, the  value zero is returned. On exit,
*          IOFF is at least zero.
*
*  JOFF    (local output) INTEGER
*          On exit, JOFF is the local column offset of the first  column
*          corresponding  to  a diagonal entry of the Virtual matrix VM.
*          If no  diagonals  are found,  the  value zero is returned. On
*          exit, JOFF is at least zero.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   int            ColCont=1, FirstD=0, GoSouth, GoEast, RowCont=1, ilow, imbloc,
                  inbloc, iupp, lcmt, lcmtnn=0, lcmt00, lmbloc, lnbloc, low, mb,
                  mblks, mbloc, mcur=0, mcurd, md=0, nb, nblks, nbloc, ncur=0,
                  ncurd, nd=0, npq=0, pmb, qnb, tmp1, tmp2, upp;
/* ..
*  .. Executable Statements ..
*
*/
   *NRPQ = 0; *NCPQ = 0; *IOFF = 0; *JOFF = 0;
   mblks = VM->mblks; nblks = VM->nblks;
/*
*  Quick return if I don't own any blocks.
*/
   if( ( mblks == 0 ) || ( nblks == 0 ) ) return;
/*
*  Retrieve the contents of VM structure fields
*/
   lcmt00 = VM->lcmt00;
   imbloc = VM->imbloc; mb  = VM->mb;  lmbloc = VM->lmbloc;
   iupp   = VM->iupp;   upp = VM->upp; pmb    = VM->nprow * mb;
   inbloc = VM->inbloc; nb  = VM->nb;  lnbloc = VM->lnbloc;
   ilow   = VM->ilow;   low = VM->low; qnb    = VM->npcol * nb;
/*
*  Handle separately the first row and/or column of the LCM table. Update the
*  LCM value of the curent block lcmt00, as well as the coordinates of the
*  current entry in the LCM table (mcur,ncur).
*/
   GoSouth = ( lcmt00 > iupp );
   GoEast  = ( lcmt00 < ilow );
/*
*  Go through the table looking for blocks owning diagonal entries.
*/
   if( !( GoSouth ) && !( GoEast ) )
   {
/*
*  The upper left block owns diagonal entries lcmt00 >= ilow && lcmt00 <= iupp
*  Compute the number of diagonals in this block as well as lcm value (lcntnn)
*  that its neighbor should have to preserve continuity.
*/
      if( lcmt00 >= 0 )
      {
         tmp2 = ( ( tmp1 = imbloc - lcmt00 ) > 0 ? tmp1 : 0 );
         if(       tmp2 <  inbloc ) { npq = tmp2;   lcmtnn = -npq;         }
         else if ( tmp2 == inbloc ) { npq = inbloc; lcmtnn = 0;            }
         else                       { npq = inbloc; lcmtnn = lcmt00 + npq; }
         *IOFF += lcmt00;
      }
      else
      {
         tmp2 = ( ( tmp1 = inbloc + lcmt00 ) > 0 ? tmp1 : 0 );
         if(       tmp2 <  imbloc ) { npq = tmp2;   lcmtnn = npq;          }
         else if ( tmp2 == imbloc ) { npq = tmp2;   lcmtnn = 0;            }
         else                       { npq = imbloc; lcmtnn = lcmt00 - npq; }
         *JOFF -= lcmt00;
      }
/*
*  Save coordinates of last block owning diagonals. Set FirstD to one, since
*  a block owning diagonals has been found.
*/
      md = 0; nd = 0; FirstD = 1;
/*
*  Those rows and columns are obviously contiguous
*/
      *NRPQ = *NCPQ = npq;
/*
*  Decide whether one should go south or east in the table: Go east if the
*  block below the current one only owns lower entries. If this block, however,
*  owns diagonals, then go south.
*/
      GoSouth = !( GoEast = ( lcmt00 - iupp + upp - pmb < ilow ) );
   }

   if( GoSouth )
   {
/*
*  Go one step south in the LCM table. Adjust the current LCM value.
*/
      lcmt00 -= iupp - upp + pmb; mcur++;
      if( !FirstD ) *IOFF += imbloc;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value accordingly.
*/
      while( ( mcur < mblks ) && ( lcmt00 > upp ) )
      {
         lcmt00 -= pmb;
         mcur++;
         if( !FirstD ) *IOFF += mb;
      }
/*
*  Return if no more row in the LCM table.
*/
      if( mcur >= mblks ) goto l_end;
/*
*  lcmt00 <= upp. The current block owns either diagonals or lower entries.
*  Save the current position in the LCM table. After this column has been
*  completely taken care of, re-start from this row and the next column of
*  the LCM table.
*/
      lcmt = lcmt00; mbloc = mb; mcurd = mcur;

      while( ( mcurd < mblks ) && ( lcmt >= ilow ) )
      {
         if( mcurd == mblks-1 ) mbloc = lmbloc;
/*
*  A block owning diagonals lcmt00 >= ilow && lcmt00 <= upp has been found.
*  If this is not the first one, update the booleans telling if the rows
*  and/or columns are contiguous.
*/
         if( FirstD )
         {
            RowCont = RowCont &&
            ( ( ( mcurd == md+1 ) &&  ( lcmtnn  <= 0 )  &&  ( lcmt <= 0 ) ) ||
              ( ( mcurd == md ) && ( ncur == nd+1 ) && ( lcmtnn == lcmt ) ) );
            ColCont = ColCont &&
            ( ( ( ncur == nd+1 )  &&  ( lcmtnn  >= 0 )  &&  ( lcmt >= 0 ) ) ||
              ( ( ncur == nd ) && ( mcurd == md+1 ) && ( lcmtnn == lcmt ) ) );
         }
/*
*  Compute the number of diagonals in this block as well as lcm value (lcntnn)
*  that its neighbor should have to preserve continuity.
*/
         if( lcmt >= 0 )
         {
            tmp2 = ( ( tmp1 = mbloc - lcmt ) > 0 ? tmp1 : 0 );
            if(       tmp2 <  inbloc ) { npq = tmp2;   lcmtnn = -npq;       }
            else if ( tmp2 == inbloc ) { npq = inbloc; lcmtnn = 0;          }
            else                       { npq = inbloc; lcmtnn = lcmt + npq; }
            if( !FirstD ) *IOFF += lcmt;
         }
         else
         {
            tmp2 = ( ( tmp1 = inbloc + lcmt ) > 0 ? tmp1 : 0 );
            if(       tmp2 <  mbloc ) { npq = tmp2;  lcmtnn = npq;        }
            else if ( tmp2 == mbloc ) { npq = tmp2;  lcmtnn = 0;          }
            else                      { npq = mbloc; lcmtnn = lcmt - npq; }
            if( !FirstD ) *JOFF -= lcmt;
         }
/*
*  Save coordinates of last block owning diagonals. Set FirstD to one, since
*  a block owning diagonals has been found.
*/
         md = mcurd; nd = ncur; FirstD = 1;
/*
*  If rows (resp columns) are still contiguous, add those npq rows (resp.
*  columns).
*/
         if( RowCont ) *NRPQ += npq;
         if( ColCont ) *NCPQ += npq;
/*
*  Keep going south until there are no more blocks owning diagonals
*/
         lcmt00  = lcmt;
         lcmt   -= pmb;
         mcur    = mcurd++;
      }
/*
*  I am done with the first column of the LCM table. Go to the next column.
*/
      lcmt00 += low - ilow + qnb; ncur++;
      if( !FirstD ) *JOFF += inbloc;
   }
   else if( GoEast )
   {
/*
*  Go one step east in the LCM table. Adjust the current LCM value.
*/
      lcmt00 += low - ilow + qnb; ncur++;
      if( !FirstD ) *JOFF += inbloc;
/*
*  While there are blocks remaining that own lower entries, keep going east
*  in the LCM table. Adjust the current LCM value.
*/
      while( ( ncur < nblks ) && ( lcmt00 < low ) )
      {
         lcmt00 += qnb;
         ncur++;
         if( !FirstD ) *JOFF += nb;
      }
/*
*  Return if no more column in the LCM table.
*/
      if( ncur >= nblks ) goto l_end;
/*
*  lcmt00 >= low. The current block owns either diagonals or upper entries. Save
*  the current position in the LCM table. After this row has been completely
*  taken care of, re-start from this column and the next row of the LCM table.
*/
      lcmt = lcmt00; nbloc = nb; ncurd = ncur;

      while( ( ncurd < nblks ) && ( lcmt <= iupp ) )
      {
         if( ncurd == nblks-1 ) nbloc = lnbloc;
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= iupp has been found.
*  If this is not the first one, update the booleans telling if the rows
*  and/or columns are contiguous.
*/
         if( FirstD )
         {
            RowCont = RowCont &&
            ( ( ( mcur == md+1 )  &&   ( lcmtnn <= 0 )  &&  ( lcmt <= 0 ) ) ||
              ( ( mcur == md ) && ( ncurd == nd+1 ) && ( lcmtnn == lcmt ) ) );
            ColCont = ColCont &&
            ( ( ( ncurd == nd+1 ) &&   ( lcmtnn >= 0 )  &&  ( lcmt >= 0 ) ) ||
              ( ( ncurd == nd ) && ( mcur == md+1 ) && ( lcmtnn == lcmt ) ) );
         }
/*
*  Compute the number of diagonals in this block as well as lcm value (lcntnn)
*  that its neighbor should have to preserve continuity.
*/
         if( lcmt >= 0 )
         {
            tmp2 = ( ( tmp1 = imbloc - lcmt ) > 0 ? tmp1 : 0 );
            if(       tmp2 <  nbloc ) { npq = tmp2;  lcmtnn = -npq;       }
            else if ( tmp2 == nbloc ) { npq = nbloc; lcmtnn = 0;          }
            else                      { npq = nbloc; lcmtnn = lcmt + npq; }
            if( !FirstD ) *IOFF += lcmt;
         }
         else
         {
            tmp2 = ( ( tmp1 = nbloc + lcmt ) > 0 ? tmp1 : 0 );
            if(       tmp2 <  imbloc ) { npq = tmp2;   lcmtnn = npq;        }
            else if ( tmp2 == imbloc ) { npq = tmp2;   lcmtnn = 0;          }
            else                       { npq = imbloc; lcmtnn = lcmt - npq; }
            if( !FirstD ) *JOFF -= lcmt;
         }
/*
*  Save coordinates of last block owning diagonals. Set FirstD to one, since
*  a block owning diagonals has been found.
*/
         md = mcur; nd = ncurd; FirstD = 1;
/*
*  If rows (resp columns) are still contiguous, add those npq rows (resp.
*  columns).
*/
         if( RowCont ) *NRPQ += npq;
         if( ColCont ) *NCPQ += npq;
/*
*  Keep going east until there are no more blocks owning diagonals.
*/
         lcmt00 = lcmt; lcmt += qnb; ncur = ncurd++;
      }
/*
*  I am done with the first row of the LCM table. Go to the next row.
*/
      lcmt00 -= iupp - upp + pmb; mcur++;
      if( !FirstD ) *IOFF += imbloc;
   }
/*
*  Loop over the remaining columns of the LCM table.
*/
   nbloc = nb;

   while( ( RowCont || ColCont ) && ( ncur  < nblks ) )
   {
      if( ncur == nblks-1 ) nbloc = lnbloc;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value accordingly.
*/
      while( ( mcur  < mblks ) && ( lcmt00 > upp ) )
      {
         lcmt00 -= pmb;
         mcur++;
         if( !FirstD ) *IOFF += mb;
      }
/*
*  Return if no more row in the LCM table.
*/
      if( mcur >= mblks ) goto l_end;
/*
*  lcmt00 <= upp. The current block owns either diagonals or lower entries.
*  Save the current position in the LCM table. After this column has been
*  completely taken care of, re-start from this row and the next column of
*  the LCM table.
*/
      lcmt = lcmt00; mbloc = mb; mcurd = mcur;

      while( ( mcurd < mblks ) && ( lcmt >= low ) )
      {
         if( mcurd == mblks-1 ) mbloc = lmbloc;
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
*  If this is not the first one, update the booleans telling if the rows
*  and/or columns are contiguous.
*/
         if( FirstD )
         {
            RowCont = RowCont &&
            ( ( ( mcurd == md+1 )  &&  ( lcmtnn <= 0 )  &&  ( lcmt <= 0 ) ) ||
              ( ( mcurd == md ) && ( ncur == nd+1 ) && ( lcmtnn == lcmt ) ) );
            ColCont = ColCont &&
            ( ( ( ncur == nd+1 )   &&  ( lcmtnn >= 0 )  &&  ( lcmt >= 0 ) ) ||
              ( ( ncur == nd ) && ( mcurd == md+1 ) && ( lcmtnn == lcmt ) ) );
         }
/*
*  Compute the number of diagonals in this block as well as lcm value (lcntnn)
*  that its neighbor should have to preserve continuity.
*/
         if( lcmt >= 0 )
         {
            tmp2 = ( ( tmp1 = mbloc - lcmt ) > 0 ? tmp1 : 0 );
            if(       tmp2 <  nbloc ) { npq = tmp2;  lcmtnn = -npq;       }
            else if ( tmp2 == nbloc ) { npq = nbloc; lcmtnn = 0;          }
            else                      { npq = nbloc; lcmtnn = lcmt + npq; }
            if( !FirstD ) *IOFF += lcmt;
         }
         else
         {
            tmp2 = ( ( tmp1 = nbloc + lcmt ) > 0 ? tmp1 : 0 );
            if(       tmp2 <  mbloc ) { npq = tmp2;  lcmtnn = npq;        }
            else if ( tmp2 == mbloc ) { npq = tmp2;  lcmtnn = 0;          }
            else                      { npq = mbloc; lcmtnn = lcmt - npq; }
            if( !FirstD ) *JOFF -= lcmt;
         }
/*
*  Save coordinates of last block owning diagonals. Set FirstD to one, since
*  a block owning diagonals has been found.
*/
         md = mcurd; nd = ncur; FirstD = 1;
/*
*  If rows (resp columns) are still contiguous, add those npq rows (resp.
*  columns).
*/
         if( RowCont ) *NRPQ += npq;
         if( ColCont ) *NCPQ += npq;
/*
*  Keep going south until there are no more blocks owning diagonals
*/
         lcmt00  = lcmt;
         lcmt   -= pmb;
         mcur    = mcurd++;
      }
/*
*  I am done with this column of the LCM table. Go to the next column until
*  there are no more column in the table.
*/
      lcmt00 += qnb; ncur++;
      if( !FirstD ) *JOFF += nb;
   }

l_end:
/*
*  If no diagonals were found, reset IOFF and JOFF to zero.
*/
   if( !FirstD ) { *IOFF = 0; *JOFF = 0; }
/*
*  End of PB_CVMcontig
*/
}
