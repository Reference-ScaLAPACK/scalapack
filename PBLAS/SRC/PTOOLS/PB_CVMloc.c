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
Int PB_CVMloc( PBTYP_T * TYPE, PB_VM_T * VM,
               char * VROCS, char * ROCS, char * UNPA, char * TRANS,
               Int MN, Int K, char * ALPHA, char * A, Int LDA,
               char * BETA,  char * B, Int LDB )
#else
Int PB_CVMloc( TYPE, VM, VROCS, ROCS, UNPA, TRANS, MN, K, ALPHA, A,
               LDA, BETA,  B, LDB )
/*
*  .. Scalar Arguments ..
*/
   Int            K, LDA, LDB, MN;
   char           * ALPHA, * BETA;
/*
*  .. Array Arguments ..
*/
   char           * VROCS, * ROCS, * UNPA, * TRANS;
   PBTYP_T        * TYPE;
   PB_VM_T        * VM;
   char           * A, * B;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CVMloc  packs  a  one-dimensional distributed array A into another
*  one-dimensional distributed array  B,  or  unpacks  a one-dimensional
*  distributed array B into a one-dimensional distributed array A.  This
*  operation is triggered by a virtual distributed array.
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
*          tual distributed array grid should be used for the packing or
*          unpacking operation as follows:
*             VROCS = 'R' or 'r', the rows should be used,
*             VROCS = 'C' or 'c', the columns should be used.
*
*  ROCS    (local input) pointer to CHAR
*          On entry, ROCS specifies if rows or columns should be  packed
*          or unpacked as follows:
*             ROCS = 'R' or 'r', rows should be (un)packed,
*             ROCS = 'C' or 'c', columns should be (un)packed.
*
*  UNPA    (local input) pointer to CHAR
*          On entry,  UNPA specifies if the data should be packed or un-
*          packed as follows:
*             UNPA = 'P' or 'p', packing (A into B),
*             UNPA = 'U' or 'u', unpacking (B into A).
*
*  TRANS   (local input) pointer to CHAR
*          On entry,  TRANS  specifies if conjugation,  transposition or
*          conjugate transposition  should occur during the  (un)packing
*          operation as follows:
*             TRANS = 'N' or 'n', natural (un)packing,
*             TRANS = 'Z' or 'z', conjugated (un)packing,
*             TRANS = 'T' or 'T', transposed (un)packing,
*             TRANS = 'C' or 'c', conjugate transposed (un)packing.
*
*  MN      (local input) INTEGER
*          On entry,  MN  specifies  the number of rows or columns to be
*          (un)packed. MN must be at least zero.
*
*  K       (local input) INTEGER
*          On entry,  K  specifies the length of the non-distributed di-
*          mension to be (un)packed. K must be at least zero.
*
*  ALPHA   (local input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (local input/local output) pointer to CHAR
*          On entry, A points to an array of  dimension (LDA, Ka), where
*          Ka is K when ROCS is 'R' or 'r' and  when ROCS is 'C' or 'c',
*          Ka is IMBLOC+(MBLKS-2)*MB+LMB  when  VROCS  is 'R' or 'r' and
*          when VROCS is 'C' or 'c', Ka is INBLOC+(NBLKS-2)*NB+LNB. This
*          array contains unpacked data.
*
*  LDA     (local input) INTEGER
*          On entry, LDA specifies the leading dimension of the array A.
*          LDA must be at least MAX( 1, K ) when ROCS = 'C' or  'c'  and
*          MAX( 1, IMBLOC+(MBLKS-2)*MB+LMB ) when ROCS is 'R' or 'r' and
*          VROCS  is  'R'  or 'r', and MAX( 1, INBLOC+(NBLKS-2)*NB+LNB )
*          when ROCS is 'R' or 'r' and VROCS is 'C' or 'c'.
*
*  BETA    (local input) pointer to CHAR
*          On entry, BETA specifies the scalar beta.
*
*  B       (local input/local output) pointer to CHAR
*          On entry, B points to an array of  dimension (LDB, Kb).  When
*          TRANS  is  'N', 'n', 'Z' or 'z',  Kb  is K when ROCS is 'R'or
*          'r', and when ROCS is 'C' or 'c', Kb  is IMBLOC+(MBLKS-2)*MB+
*          LMB  when  VROCS  is 'C' or 'c' and when VROCS is 'R', Kb  is
*          INBLOC+(NBLKS-2)*NB+LNB.  When TRANS is 'T', 't', 'C' or 'c',
*          Kb is K when ROCS is 'C' or 'c' and  when ROCS is 'R' or 'r',
*          Kb  is  IMBLOC+(MBLKS-2)*MB+LMB when  VROCS is 'C' or 'c' and
*          when VROCS is 'R' or 'r', Kb is INBLOC+(NBLKS-2)*NB+LNB. This
*          array contains unpacked data.
*
*  LDB     (local input) INTEGER
*          On entry, LDB specifies the leading dimension of the array B.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   Int            GoEast, GoSouth, ilow, imbloc, inbloc, inca, incb, iupp, kb,
                  lcmt, lcmt00, lmbloc, lnbloc, low, mb, mblkd, mblks, mbloc,
                  * m, * n, nb, nblkd, nblks, nbloc, notran, npcol, npq=0,
                  nprow, pmb, qnb, rows, size, tmp1, tmp2, upp;
   MMADD_T        add;
   char           * aptrd, * bptrd;
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

   if( Mupcase( UNPA[0] ) == CPACKING )
   {
/*
*  B is the distributed target, A is the distributed source
*/
      if( Mupcase( TRANS[0] ) == CNOTRAN )
      {
/*
*  Add A to B
*/
         notran = 1;
         add    = TYPE->Fmmadd;
      }
      else if( Mupcase( TRANS[0] ) == CCONJG )
      {
/*
*  Add the conjugate of A to B
*/
         notran = 1;
         add    = TYPE->Fmmcadd;
      }
      else if( Mupcase( TRANS[0] ) == CTRAN )
      {
/*
*  Add the tranpose of A to B
*/
         notran = 0;
         add    = TYPE->Fmmtadd;
      }
      else
      {
/*
*  Add the conjugate tranpose of A to B
*/
         notran = 0;
         add    = TYPE->Fmmtcadd;
      }
   }
   else
   {
/*
*  A is the distributed target, B is the distributed source
*/
      if( Mupcase( TRANS[0] ) == CNOTRAN )
      {
/*
*  Add B to A
*/
         notran = 1;
         add    = TYPE->Fmmdda;
      }
      else if( Mupcase( TRANS[0] ) == CCONJG )
      {
/*
*  Add the conjugate of B to A
*/
         notran = 1;
         add    = TYPE->Fmmddac;
      }
      else if( Mupcase( TRANS[0] ) == CTRAN )
      {
/*
*  Add the tranpose of B to A
*/
         notran = 0;
         add    = TYPE->Fmmddat;
      }
      else
      {
/*
*  Add the conjugate tranpose of B to A
*/
         notran = 0;
         add    = TYPE->Fmmddact;
      }
   }

   size = TYPE->size;
   rows = ( Mupcase( ROCS[0] ) == CROW );

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
         inca = size;
         incb = ( notran ? size : LDB * size );
         m    = &tmp2;
         n    = &K;
      }
      else
      {
/*
*  (un)packing columns of k by mn array A
*/
         inca = LDA * size;
         incb = ( notran ? LDB * size : size );
         m    = &K;
         n    = &tmp2;
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
            if( rows ) add( &npq, &K, ALPHA, A, &LDA, BETA, B, &LDB );
            else       add( &K, &npq, ALPHA, A, &LDA, BETA, B, &LDB );
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
            tmp1 = imbloc - lcmt00; tmp1 = MAX( 0, tmp1 );
            tmp2 = MIN( tmp1, inbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
            add( m, n, ALPHA, A+lcmt00*inca, &LDA, BETA, B, &LDB );
         }
         else
         {
            tmp1 = inbloc + lcmt00; tmp1 = MAX( 0, tmp1 );
            tmp2 = MIN( tmp1, imbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
            add( m, n, ALPHA, A, &LDA, BETA, B-lcmt00*incb, &LDB );
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
*  the pointer to A. The pointer to B remains unchanged.
*/
         lcmt00 -= iupp - upp + pmb; mblks--; A += imbloc * inca;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value as well as the pointer to A accordingly.
*/
         while( mblks && ( lcmt00 > upp ) )
         { lcmt00 -= pmb; mblks--; A += mb * inca; }
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
         lcmt = lcmt00; mblkd = mblks; aptrd = A;

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
               add( m, n, ALPHA, aptrd+lcmt*inca, &LDA, BETA, B, &LDB );
            }
            else
            {
               tmp1 = inbloc + lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, mbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               add( m, n, ALPHA, aptrd, &LDA, BETA, B-lcmt*incb, &LDB );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Keep going south until there are no more blocks owning diagonals
*/
            lcmt -= pmb; mblkd--; aptrd += mbloc * inca;
         }
/*
*  I am done with the first column of the LCM table. Go to the next column.
*/
         lcmt00 += low - ilow + qnb; nblks--; B += inbloc * incb;
      }
      else if( GoEast )
      {
/*
*  Go one step east in the LCM table. Adjust the current LCM value as
*  well as the pointer to B. The pointer to A remains unchanged.
*/
         lcmt00 += low - ilow + qnb; nblks--; B += inbloc * incb;
/*
*  While there are blocks remaining that own lower entries, keep going east
*  in the LCM table. Adjust the current LCM value as well as the pointer to
*  B accordingly.
*/
         while( nblks && ( lcmt00 < low ) )
         { lcmt00 += qnb; nblks--; B += nb * incb; }
/*
*  Return if no more column in the LCM table.
*/
         if( nblks <= 0 ) return( npq );
/*
*  lcmt00 >= low. The current block owns either diagonals or upper entries. Save
*  the current position in the LCM table. After this row has been completely
*  taken care of, re-start from this column and the next row of the LCM table.
*/
         lcmt = lcmt00; nblkd = nblks; bptrd = B;

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
               add( m, n, ALPHA, A+lcmt*inca, &LDA, BETA, bptrd, &LDB );
            }
            else
            {
               tmp1 = nbloc + lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, imbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               add( m, n, ALPHA, A, &LDA, BETA, bptrd-lcmt*incb, &LDB );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Keep going east until there are no more blocks owning diagonals.
*/
            lcmt += qnb; nblkd--; bptrd += nbloc * incb;
         }
/*
*  I am done with the first row of the LCM table. Go to the next row.
*/
         lcmt00 -= iupp - upp + pmb; mblks--; A += imbloc * inca;
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
               { lcmt00 -= pmb; mblks--; A += mb * inca; }
               if( lcmt00 >= low ) break;
               while( nblks && ( lcmt00 < low ) )
               { lcmt00 += qnb; nblks--; B += nb * incb; }
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
         lcmt = lcmt00; mblkd = mblks; aptrd = A;

         while( mblkd && ( lcmt >= low ) )
         {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
*/
            mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
            if( lcmt >= 0 )
            {
               tmp1 = mbloc - lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, nbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               add( m, n, ALPHA, aptrd+lcmt*inca, &LDA, BETA, B, &LDB );
            }
            else
            {
               tmp1 = nbloc + lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, mbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               add( m, n, ALPHA, aptrd, &LDA, BETA, B-lcmt*incb, &LDB );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Keep going south until there are no more blocks owning diagonals
*/
            lcmt -= pmb; mblkd--; aptrd += mbloc * inca;
         }
/*
*  I am done with this column of the LCM table. Go to the next column ...
*/
         lcmt00 += qnb; nblks--; B += nbloc * incb;
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
         inca = size;
         incb = ( notran ? size : LDB * size );
         m    = &tmp2;
         n    = &K;
      }
      else
      {
/*
*  (un)packing columns of k by mn array A
*/
         inca = LDA * size;
         incb = ( notran ? LDB * size : size );
         m    = &K;
         n    = &tmp2;
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
            npq = ( ( nblks <  2 ) ? inbloc :
                    inbloc + ( nblks - 2 ) * nb + lnbloc );
            npq = MIN( npq, kb );
            if( rows ) add( &npq, &K, ALPHA, A, &LDA, BETA, B, &LDB );
            else       add( &K, &npq, ALPHA, A, &LDA, BETA, B, &LDB );
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
            add( m, n, ALPHA, A, &LDA, BETA, B+lcmt00*incb, &LDB );
         }
         else
         {
            tmp1 = inbloc + lcmt00; tmp1 = MAX( 0, tmp1 );
            tmp2 = MIN( tmp1, imbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
            add( m, n, ALPHA, A-lcmt00*inca, &LDA, BETA, B, &LDB );
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
*  the pointer to B. The pointer to A remains unchanged.
*/
         lcmt00 -= iupp - upp + pmb; mblks--; B += imbloc * incb;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value as well as the pointer to B accordingly.
*/
         while( mblks && ( lcmt00 > upp ) )
         { lcmt00 -= pmb; mblks--; B += mb * incb; }
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
         lcmt  = lcmt00; mblkd = mblks; bptrd = B;

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
               add( m, n, ALPHA, A, &LDA, BETA, bptrd+lcmt*incb, &LDB );
            }
            else
            {
               tmp1 = inbloc + lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, mbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               add( m, n, ALPHA, A-lcmt*inca, &LDA, BETA, bptrd, &LDB );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Keep going south until there are no more blocks owning diagonals
*/
            lcmt -= pmb; mblkd--; bptrd += mbloc * incb;
         }
/*
*  I am done with the first column of the LCM table. Go to the next column.
*/
         lcmt00 += low - ilow + qnb; nblks--; A += inbloc * inca;
      }
      else if( GoEast )
      {
/*
*  Go one step east in the LCM table. Adjust the current LCM value as
*  well as the pointer to A. The pointer to B remains unchanged.
*/
         lcmt00 += low - ilow + qnb; nblks--; A += inbloc * inca;
/*
*  While there are blocks remaining that own lower entries, keep going east
*  in the LCM table. Adjust the current LCM value as well as the pointer to
*  A accordingly.
*/
         while( nblks && ( lcmt00 < low ) )
         { lcmt00 += qnb; nblks--; A += nb * inca; }
/*
*  Return if no more column in the LCM table.
*/
         if( nblks <= 0 ) return( npq );
/*
*  lcmt00 >= low. The current block owns either diagonals or upper entries. Save
*  the current position in the LCM table. After this row has been completely
*  taken care of, re-start from this column and the next row of the LCM table.
*/
         lcmt  = lcmt00; nblkd = nblks; aptrd = A;

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
               add( m, n, ALPHA, aptrd, &LDA, BETA, B+lcmt*incb, &LDB );
            }
            else
            {
               tmp1 = nbloc + lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, imbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               add( m, n, ALPHA, aptrd-lcmt*inca, &LDA, BETA, B, &LDB );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Keep going east until there are no more blocks owning diagonals.
*/
            lcmt += qnb; nblkd--; aptrd += nbloc * inca;
         }
/*
*  I am done with the first row of the LCM table. Go to the next row.
*/
         lcmt00 -= iupp - upp + pmb; mblks--; B += imbloc * incb;
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
               { lcmt00 -= pmb; mblks--; B += mb * incb; }
               if( lcmt00 >= low ) break;
               while( nblks && ( lcmt00 < low ) )
               { lcmt00 += qnb; nblks--; A += nb * inca; }
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
         lcmt = lcmt00; mblkd = mblks; bptrd = B;

         while( mblkd && ( lcmt >= low ) )
         {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
*/
            mbloc = ( ( mblkd == 1 ) ? lmbloc : mb );
            if( lcmt >= 0 )
            {
               tmp1 = mbloc - lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, nbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               add( m, n, ALPHA, A, &LDA, BETA, bptrd+lcmt*incb, &LDB );
            }
            else
            {
               tmp1 = nbloc + lcmt; tmp1 = MAX( 0, tmp1 );
               tmp2 = MIN( tmp1, mbloc ); npq += ( tmp2 = MIN( tmp2, kb ) );
               add( m, n, ALPHA, A-lcmt*inca, &LDA, BETA, bptrd, &LDB );
            }
            if( ( kb -= tmp2 ) == 0 ) return( npq );
/*
*  Keep going south until there are no more blocks owning diagonals
*/
            lcmt -= pmb; mblkd--; bptrd += mbloc * incb;
         }
/*
*  I am done with this column of the LCM table. Go to the next column ...
*/
         lcmt00 += qnb; nblks--; A += nbloc * inca;
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
*  End of PB_CVMloc
*/
}
