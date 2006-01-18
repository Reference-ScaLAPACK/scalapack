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
void PB_CVMupdate( PB_VM_T * VM, int K, int * II, int * JJ )
#else
void PB_CVMupdate( VM, K, II, JJ )
/*
*  .. Scalar Arguments ..
*/
   int            * II, * JJ, K;
   PB_VM_T        * VM;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CVMupdate  updates  the  local information of an m by n local array
*  owned by the process of  relative  coordinates ( MRROW, MRCOL ). Note
*  that if m or n is less or equal than zero, there is no data, in which
*  case this process  does  not  need  the local information computed by
*  this routine to proceed.
*
*  Arguments
*  =========
*
*  VM      (local input) pointer to a PB_VM_T structure
*          On entry,  VM  is  a  pointer to a structure of type PB_VM_T,
*          that contains the virtual matrix information (see pblas.h).
*
*  K       (global input) INTEGER
*          On entry, K  specifies  the  number of diagonal elements that
*          have been used so far. K must be at least zero.
*
*  II      (local input/local output) INTEGER
*          On entry, II  specifies the local row index to be updated. On
*          exit, II points to the local row owning the K+1 diagonal   of
*          this local block. On entry and on exit, II is at least zero.
*
*  JJ      (local input/local output) INTEGER
*          On entry, JJ  specifies the local column index to be updated.
*          On exit,  JJ  points  to  the  local  column  owning  the K+1
*          diagonal of this local block.  On entry and on exit, JJ is at
*          least zero.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   int            GoEast, GoSouth, ilow, imbloc, inbloc, ioff, ioffd, iupp,
                  joff, joffd, lcmt, lcmt00, lmbloc, lnbloc, low, mb, mblkd,
                  mblks, mbloc, nb, nblkd, nblks, nbloc, npq=0, pmb, qnb,
                  tmp1, tmp2, upp;
/* ..
*  .. Executable Statements ..
*
*/
   mblks = VM->mblks; nblks = VM->nblks;
/*
*  Quick return if I don't own any blocks or if no diagonals were found.
*/
   if( ( K <= 0 ) || ( mblks == 0 ) || ( nblks == 0 ) ) return;
/*
*  Handle the first block of rows or columns separately
*/
   ioff    = *II;
   joff    = *JJ;
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
*  LCM value of the curent block lcmt00, as well as the number of rows and
*  columns mblks and nblks remaining in the LCM table.
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
*/
      if( lcmt00 >= 0 )
      {
         npq = ( ( tmp2 = ( tmp1 = imbloc - lcmt00 ) > 0 ? tmp1 : 0 ) <
                 inbloc ? tmp2 : inbloc );
         if( K < npq )
         {
            tmp1   = lcmt00 + K;
            *II   += tmp1;
            iupp   = ( imbloc -= tmp1 ) - 1;
            if( mblks == 1 ) lmbloc = imbloc;
            *JJ   += K;
            ilow   = 1 - ( inbloc -= K );
            if( nblks == 1 ) lnbloc = inbloc;
            lcmt00 = 0;
            goto l_end;
         }
      }
      else
      {
         npq = ( ( tmp2 = ( tmp1 = inbloc + lcmt00 ) > 0 ? tmp1 : 0 ) >
                 imbloc ? imbloc : tmp2 );
         if( K < npq )
         {
            tmp1   = lcmt00 - K;
            *JJ   -= tmp1;
            ilow   = 1 - ( inbloc += tmp1 );
            if( nblks == 1 ) lnbloc = inbloc;
            *II   += K;
            iupp   = ( imbloc -= K ) - 1;
            if( mblks == 1 ) lmbloc = imbloc;
            lcmt00 = 0;
            goto l_end;
         }
      }
      K -= npq;
/*
*  Decide whether one should go south or east in the table: Go east if the
*  block below the current one only owns lower entries. If this block,
*  however, owns diagonals, then go south.
*/
      GoSouth = !( GoEast = ( ( lcmt00 - ( iupp-upp+pmb ) ) < ilow ) );
/*
*  Update the local indexes II and JJ
*/
      if( GoSouth ) *II += imbloc;
      else          *JJ += inbloc;
   }

   if( GoSouth )
   {
/*
*  Go one step south in the LCM table. Adjust the current LCM value.
*/
      lcmt00 -= iupp - upp + pmb; mblks--; ioff += imbloc;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value accordingly.
*/
      while( mblks && ( lcmt00 > upp ) ) { lcmt00 -= pmb; mblks--; ioff += mb; }
/*
*  Return if no more row in the LCM table.
*/
      if( mblks <= 0 ) goto l_end;
/*
*  lcmt00 <= upp. The current block owns either diagonals or lower entries.
*  Save the current position in the LCM table. After this column has been
*  completely taken care of, re-start from this row and the next column of
*  the LCM table.
*/
      lcmt  = lcmt00; mblkd = mblks; mbloc = mb; ioffd = ioff;

      while( mblkd && ( lcmt >= ilow ) )
      {
/*
*  A block owning diagonals lcmt00 >= ilow && lcmt00 <= upp has been found.
*/
         if( mblkd == 1 ) mbloc = lmbloc;

         if( lcmt >= 0 )
         {
            npq = ( ( tmp2 = ( tmp1 = mbloc - lcmt ) > 0 ? tmp1 : 0 ) < inbloc ?
                    tmp2 : inbloc );
            if( K < npq )
            {
               tmp1   = lcmt + K;
               *II    = ioffd + tmp1;
               iupp   = ( imbloc = mbloc - tmp1 ) - 1;
               if( mblks == 1 ) lmbloc = imbloc;
               *JJ    = joff + K;
               ilow   = 1 - ( inbloc -= K );
               if( nblks == 1 ) lnbloc = inbloc;
               lcmt00 = 0;
               goto l_end;
            }
         }
         else
         {
            npq = ( ( tmp2 = ( tmp1 = inbloc + lcmt ) > 0 ? tmp1 : 0 ) > mbloc ?
                    mbloc : tmp2 );
            if( K < npq )
            {
               tmp1   = lcmt - K;
               *JJ    = joff - tmp1;
               ilow   = 1 - ( inbloc += tmp1 );
               if( nblks == 1 ) lnbloc = inbloc;
               *II    = ioffd + K;
               iupp   = ( imbloc = mbloc - K ) - 1;
               if( mblks == 1 ) lmbloc = imbloc;
               lcmt00 = 0;
               goto l_end;
            }
         }
/*
*  Keep going south until there are no more blocks owning diagonals
*/
         K      -= npq;
         lcmt00  = lcmt;
         lcmt   -= pmb;
         mblks   = mblkd--;
         ioff    = ioffd;
         ioffd  += mbloc;
      }
/*
*  I am done with the first column of the LCM table. Go to the next column.
*/
      lcmt00 += low - ilow + qnb; nblks--; joff += inbloc;
/*
*  Update the local indexes II and JJ
*/
      *II = ioff;
      *JJ = joff;
   }
   else if( GoEast )
   {
/*
*  Go one step east in the LCM table. Adjust the current LCM value.
*/
      lcmt00 += low - ilow + qnb; nblks--; joff += inbloc;
/*
*  While there are blocks remaining that own lower entries, keep going east
*  in the LCM table. Adjust the current LCM value.
*/
      while( nblks && ( lcmt00 < low ) ) { lcmt00 += qnb; nblks--; joff += nb; }
/*
*  Return if no more column in the LCM table.
*/
      if( nblks <= 0 ) goto l_end;
/*
*  lcmt00 >= low. The current block owns either diagonals or upper entries. Save
*  the current position in the LCM table. After this row has been completely
*  taken care of, re-start from this column and the next row of the LCM table.
*/
      lcmt  = lcmt00; nblkd = nblks; nbloc = nb; joffd = joff;

      while( nblkd && ( lcmt <= iupp ) )
      {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= iupp has been found.
*/
         if( nblkd == 1 ) nbloc = lnbloc;

         if( lcmt >= 0 )
         {
            npq = ( ( tmp2 = ( tmp1 = imbloc - lcmt ) > 0 ? tmp1 : 0 ) < nbloc ?
                    tmp2 : nbloc );
            if( K < npq )
            {
               tmp1   = lcmt + K;
               *II    = ioff + tmp1;
               iupp   = ( imbloc -= tmp1 ) - 1;
               if( mblks == 1 ) lmbloc = imbloc;
               *JJ    = joffd + K;
               ilow   = 1 - ( inbloc = nbloc - K );
               if( nblks == 1 ) lnbloc = inbloc;
               lcmt00 = 0;
               goto l_end;
            }
         }
         else
         {
            npq = ( ( tmp2 = ( tmp1 = nbloc + lcmt ) > 0 ? tmp1 : 0 ) > imbloc ?
                    imbloc : tmp2 );
            if( K < npq )
            {
               tmp1   = lcmt - K;
               *JJ    = joffd - tmp1;
               ilow   = 1 - ( inbloc = nbloc + tmp1 );
               if( nblks == 1 ) lnbloc = inbloc;
               *II    = ioff + K;
               iupp   = ( imbloc -= K ) - 1;
               if( mblks == 1 ) lmbloc = imbloc;
               lcmt00 = 0;
               goto l_end;
            }
         }
/*
*  Keep going east until there are no more blocks owning diagonals.
*/
         K      -= npq;
         lcmt00  = lcmt;
         lcmt   += qnb;
         nblks   = nblkd--;
         joff    = joffd;
         joffd  += nbloc;
      }
/*
*  I am done with the first row of the LCM table. Go to the next row.
*/
      lcmt00 -= iupp - upp + pmb; mblks--; ioff += imbloc;
/*
*  Update the local indexes II and JJ
*/
      *II     = ioff;
      *JJ     = joff;
   }
/*
*  Loop over the remaining columns of the LCM table.
*/
   nbloc = nb;

   while( nblks )
   {
      if( nblks == 1 ) nbloc = lnbloc;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value accordingly.
*/
      while( mblks && ( lcmt00 > upp ) ) { lcmt00 -= pmb; mblks--; ioff += mb; }
/*
*  Return if no more row in the LCM table.
*/
      if( mblks <= 0 ) goto l_end;
/*
*  lcmt00 <= upp. The current block owns either diagonals or lower entries.
*  Save the current position in the LCM table. After this column has been
*  completely taken care of, re-start from this row and the next column of
*  the LCM table.
*/
      lcmt  = lcmt00; mblkd = mblks; mbloc = mb; ioffd = ioff;

      while( mblkd && ( lcmt >= low ) )
      {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
*/
         if( mblkd == 1 ) mbloc = lmbloc;

         if( lcmt >= 0 )
         {
            npq = ( ( tmp2 = ( tmp1 = mbloc - lcmt ) > 0 ? tmp1 : 0 ) < nbloc ?
                    tmp2 : nbloc );
            if( K < npq )
            {
               tmp1   = lcmt + K;
               *II    = ioffd + tmp1;
               iupp   = ( imbloc = mbloc - tmp1 ) - 1;
               if( mblks == 1 ) lmbloc = imbloc;
               *JJ    = joff + K;
               ilow   = 1 - ( inbloc = nbloc - K );
               if( nblks == 1 ) lnbloc = inbloc;
               lcmt00 = 0;
               goto l_end;
            }
         }
         else
         {
            npq = ( ( tmp2 = ( tmp1 = nbloc + lcmt ) > 0 ? tmp1 : 0 ) > mbloc ?
                    mbloc : tmp2 );
            if( K < npq )
            {
               tmp1   = lcmt - K;
               *JJ    = joff - tmp1;
               ilow   = 1 - ( inbloc = nbloc + tmp1 );
               if( nblks == 1 ) lnbloc = inbloc;
               *II    = ioffd + K;
               iupp   = ( imbloc = mbloc - K ) - 1;
               if( mblks == 1 ) lmbloc = imbloc;
               lcmt00 = 0;
               goto l_end;
            }
         }
/*
*  Keep going south until there are no more blocks owning diagonals
*/
         K      -= npq;
         lcmt00  = lcmt;
         lcmt   -= pmb;
         mblks   = mblkd--;
         ioff    = ioffd;
         ioffd  += mbloc;
      }
/*
*  I am done with this column of the LCM table. Go to the next column until
*  there are no more column in the table.
*/
      lcmt00 += qnb; nblks--; joff += nbloc;
/*
*  Update the local indexes II and JJ
*/
      *II     = ioff;
      *JJ     = joff;
   }

l_end:
/*
*  Update the fields of the VM structure
*/
   VM->lcmt00 = lcmt00;
   VM->mp     = ( mblks >= 2 ? imbloc + ( mblks - 2 ) * mb + lmbloc :
                               ( mblks == 1 ? imbloc : 0 ) );
   VM->imbloc = imbloc; VM->lmbloc = lmbloc; VM->mblks = mblks; VM->iupp = iupp;
   VM->nq     = ( nblks >= 2 ? inbloc + ( nblks - 2 ) * nb + lnbloc :
                               ( nblks == 1 ? inbloc : 0 ) );
   VM->inbloc = inbloc; VM->lnbloc = lnbloc; VM->nblks = nblks; VM->ilow = ilow;
/*
*  End of PB_CVMupdate
*/
}
