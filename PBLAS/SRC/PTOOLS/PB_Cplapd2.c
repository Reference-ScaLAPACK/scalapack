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
void PB_Cplapd2( PBTYP_T * TYPE, char * UPLO, char * CONJUG, Int M,
                 Int N, char * ALPHA, char * BETA, char * A, Int IA,
                 Int JA, Int * DESCA )
#else
void PB_Cplapd2( TYPE, UPLO, CONJUG, M, N, ALPHA, BETA, A, IA, JA,
                 DESCA )
/*
*  .. Scalar Arguments ..
*/
   char           * CONJUG, * UPLO;
   Int            IA, JA, M, N;
   char           * ALPHA, * BETA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA;
   char           * A;
#endif
{
/*
*  .. Local Scalars ..
*/
   char           UploA, herm;
   Int            Acol, Aii, Aimb1, Ainb1, Aoffi, Ajj, Ald, Amb, Amp, Anb, Anq,
                  Aoffj, Arcol, Arow, Arrow, GoEast, GoSouth, iimax, ilow,
                  imbloc, inbloc, ioffd, iupp, izero=0, jjmax, joffd, lcmt,
                  lcmt00, lmbloc, lnbloc, low, lower, m1, mbloc, mblkd, mblks,
                  mycol, myrow, n1, nbloc, nblkd, nblks, npcol, nprow, pmb, qnb,
                  size, tmp1, upp, upper;
   TZPAD_T        pad;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Quick return if possible
*/
   if( ( M <= 0 ) || ( N <= 0 ) ) return;
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( DESCA[CTXT_], &nprow, &npcol, &myrow, &mycol );
/*
*  Retrieve sub( A )'s local information: Aii, Ajj, Arow, Acol ...
*/
   PB_Cainfog2l( M, N, IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aimb1,
                 &Ainb1, &Amp, &Anq, &Aii, &Ajj, &Arow, &Acol, &Arrow, &Arcol );
/*
*  Quick return if I don't own any of sub( A ).
*/
   if( ( Amp <= 0 ) || ( Anq <= 0 ) ) return;
/*
*  Initialize lcmt00, mblks, nblks, imbloc, inbloc, lmbloc, lnbloc, ilow, low,
*  iupp, and upp.
*/
   Amb   = DESCA[MB_ ]; Anb   = DESCA[NB_ ]; Ald   = DESCA[LLD_];
   PB_Cbinfo( 0, Amp, Anq, Aimb1, Ainb1, Amb, Anb, Arrow, Arcol, &lcmt00,
              &mblks, &nblks, &imbloc, &inbloc, &lmbloc, &lnbloc, &ilow, &low,
              &iupp, &upp );

   iimax = ( Aoffi = Aii - 1 ) + ( m1 = Amp );
   jjmax = ( Aoffj = Ajj - 1 ) + ( n1 = Anq );
   pmb = ( ( ( Arow < 0 ) || ( nprow == 1 ) ) ? Amb : nprow * Amb );
   qnb = ( ( ( Acol < 0 ) || ( npcol == 1 ) ) ? Anb : npcol * Anb );

   size  = TYPE->size; pad   = TYPE->Ftzpad;
   UploA = Mupcase( UPLO[0] );
   herm  = ( UploA == CALL ? CNOCONJG : Mupcase( CONJUG[0] ) );
   upper = ( UploA != CLOWER );
   lower = ( UploA != CUPPER );
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
   if( ( !( GoSouth ) ) && ( !( GoEast ) ) )
   {
/*
*  The upper left block owns diagonal entries lcmt00 >= ilow && lcmt00 <= iupp
*/
      pad( C2F_CHAR( UPLO ), C2F_CHAR( &herm ), &imbloc, &inbloc, &lcmt00,
           ALPHA, BETA, Mptr( A, Aii, Ajj, Ald, size ), &Ald );
/*
*  Decide whether one should go south or east in the table: Go east if
*  the block below the current one only owns lower entries. If this block,
*  however, owns diagonals, then go south.
*/
      GoSouth = !( GoEast = ( ( lcmt00 - ( iupp - upp + pmb ) ) < ilow ) );

      if( GoSouth )
      {
/*
*  When the upper triangular part of sub( A ) should be set and one is
*  planning to go south in the table, it is neccessary to take care of the
*  remaining columns of these imbloc rows immediately.
*/
         if( upper && ( Anq > inbloc ) )
         {
            tmp1 = Anq - inbloc;
            pad( C2F_CHAR( ALL ), C2F_CHAR( &herm ), &imbloc, &tmp1, &izero,
                 ALPHA, ALPHA, Mptr( A, Aii, Ajj+inbloc, Ald, size ), &Ald );
         }
         Aii += imbloc;
         m1  -= imbloc;
      }
      else
      {
/*
*  When the lower triangular part of sub( A ) should be set and one is
*  planning to go east in the table, it is neccessary to take care of the
*  remaining rows of these inbloc columns immediately.
*/
         if( lower && ( Amp > imbloc ) )
         {
            tmp1 = Amp - imbloc;
            pad( C2F_CHAR( ALL ), C2F_CHAR( &herm ), &tmp1, &inbloc, &izero,
                 ALPHA, ALPHA, Mptr( A, Aii+imbloc, Ajj, Ald, size ), &Ald );
         }
         Ajj += inbloc;
         n1  -= inbloc;
      }
   }

   if( GoSouth )
   {
/*
*  Go one step south in the LCM table. Adjust the current LCM value as well as
*  the local row index in A.
*/
      lcmt00 -= ( iupp - upp + pmb ); mblks--; Aoffi += imbloc;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value as well as the local row index in A.
*/
      while( ( mblks > 0 ) && ( lcmt00 > upp ) )
      { lcmt00 -= pmb; mblks--; Aoffi += Amb; }
/*
*  Set the upper triangular part of sub( A ) we just skipped when necessary.
*/
      tmp1 = MIN( Aoffi, iimax ) - Aii + 1;
      if( upper && ( tmp1 > 0 ) )
      {
         pad( C2F_CHAR( ALL ), C2F_CHAR( &herm ), &tmp1, &n1, &izero, ALPHA,
              ALPHA, Mptr( A, Aii, Aoffj+1, Ald, size ), &Ald );
         Aii += tmp1;
         m1  -= tmp1;
      }
/*
*  Return if no more row in the LCM table.
*/
      if( mblks <= 0 ) return;
/*
*  lcmt00 <= upp. The current block owns either diagonals or lower entries.
*  Save the current position in the LCM table. After this column has been
*  completely taken care of, re-start from this row and the next column of
*  the LCM table.
*/
      lcmt  = lcmt00; mblkd = mblks; ioffd = Aoffi;

      mbloc = Amb;
      while( ( mblkd > 0 ) && ( lcmt >= ilow ) )
      {
/*
*  A block owning diagonals lcmt00 >= ilow && lcmt00 <= upp has been found.
*/
         if( mblkd == 1 ) mbloc = lmbloc;
         pad( C2F_CHAR( UPLO ), C2F_CHAR( &herm ), &mbloc, &inbloc, &lcmt,
              ALPHA, BETA, Mptr( A, ioffd+1, Aoffj+1, Ald, size ), &Ald );
         lcmt00  = lcmt;
         lcmt   -= pmb;
         mblks   = mblkd;
         mblkd--;
         Aoffi   = ioffd;
         ioffd  += mbloc;
      }
/*
*  Set the lower triangular part of sub( A ) when necessary.
*/
      tmp1 = m1 - ioffd + Aii - 1;
      if( lower && ( tmp1 > 0 ) )
         pad( C2F_CHAR( ALL ), C2F_CHAR( &herm ), &tmp1, &inbloc, &izero, ALPHA,
              ALPHA, Mptr( A, ioffd+1, Aoffj+1, Ald, size ), &Ald );

      tmp1    = Aoffi - Aii + 1;
      m1     -= tmp1;
      n1     -= inbloc;
      lcmt00 += low - ilow + qnb;
      nblks--;
      Aoffj  += inbloc;
/*
*  When the upper triangular part of sub( A ) should be set, take care of the
*  n1 remaining columns of these tmp1 rows immediately.
*/
      if( upper && ( tmp1 > 0 ) && ( n1 > 0 ) )
         pad( C2F_CHAR( ALL ), C2F_CHAR( &herm ), &tmp1, &n1, &izero, ALPHA,
              ALPHA, Mptr( A, Aii, Aoffj+1, Ald, size ), &Ald );
      Aii = Aoffi + 1;
      Ajj = Aoffj + 1;
   }
   else if( GoEast )
   {
/*
*  Go one step east in the LCM table. Adjust the current LCM value as well as
*  the local column index in A.
*/
      lcmt00 += low - ilow + qnb; nblks--; Aoffj += inbloc;
/*
*  While there are blocks remaining that own lower entries, keep going east.
*  Adjust the current LCM value as well as the local column index in A.
*/
      while( ( nblks > 0 ) && ( lcmt00 < low ) )
      { lcmt00 += qnb; nblks--; Aoffj += Anb; }
/*
*  Set the lower triangular part of sub( A ) we just skipped when necessary.
*/
      tmp1 = MIN( Aoffj, jjmax ) - Ajj + 1;
      if( lower && ( tmp1 > 0 ) )
      {
         pad( C2F_CHAR( ALL ), C2F_CHAR( &herm ), &m1, &tmp1, &izero, ALPHA,
              ALPHA, Mptr( A, Aii, Ajj, Ald, size ), &Ald );
         Ajj += tmp1;
         n1  -= tmp1;
      }
/*
*  Return if no more column in the LCM table.
*/
      if( nblks <= 0 ) return;
/*
*  lcmt00 >= low. The current block owns either diagonals or upper entries.
*  Save the current position in the LCM table. After this row has been
*  completely taken care of, re-start from this column and the next row of
*  the LCM table.
*/
      lcmt  = lcmt00; nblkd = nblks; joffd = Aoffj;

      nbloc = Anb;
      while( ( nblkd > 0 ) && ( lcmt <= iupp ) )
      {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= iupp has been found.
*/
         if( nblkd == 1 ) nbloc = lnbloc;
         pad( C2F_CHAR( UPLO ), C2F_CHAR( &herm ), &imbloc, &nbloc, &lcmt,
              ALPHA, BETA, Mptr( A, Aii, joffd+1, Ald, size ), &Ald );
         lcmt00  = lcmt;
         lcmt   += qnb;
         nblks   = nblkd;
         nblkd--;
         Aoffj   = joffd;
         joffd  += nbloc;
      }
/*
*  Set the upper triangular part of sub( A ) when necessary.
*/
      tmp1 = n1 - joffd + Ajj - 1;
      if( upper && ( tmp1 > 0 ) )
         pad( C2F_CHAR( ALL ), C2F_CHAR( &herm ), &imbloc, &tmp1, &izero, ALPHA,
              ALPHA, Mptr( A, Aii, joffd+1, Ald, size ), &Ald );

      tmp1    = Aoffj - Ajj + 1;
      m1     -= imbloc;
      n1     -= tmp1;
      lcmt00 -= ( iupp - upp + pmb );
      mblks--;
      Aoffi  += imbloc;
/*
*  When the lower triangular part of sub( A ) should be set, take care of the
*  m1 remaining rows of these tmp1 columns immediately.
*/
      if( lower && ( m1 > 0 ) && ( tmp1 > 0 ) )
         pad( C2F_CHAR( ALL ), C2F_CHAR( &herm ), &m1, &tmp1, &izero, ALPHA,
              ALPHA, Mptr( A, Aoffi+1, Ajj, Ald, size ), &Ald );
      Aii = Aoffi + 1;
      Ajj = Aoffj + 1;
   }
/*
*  Loop over the remaining columns of the LCM table.
*/
   nbloc = Anb;
   while( nblks > 0 )
   {
      if( nblks == 1 ) nbloc = lnbloc;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value as well as the local row index in A.
*/
      while( ( mblks > 0 ) && ( lcmt00 > upp ) )
      { lcmt00 -= pmb; mblks--; Aoffi  += Amb; }
/*
*  Set the upper triangular part of sub( A ) we just skipped when necessary.
*/
      tmp1 = MIN( Aoffi, iimax ) - Aii + 1;
      if( upper && ( tmp1 > 0 ) )
      {
         pad( C2F_CHAR( ALL ), C2F_CHAR( &herm ), &tmp1, &n1, &izero, ALPHA,
              ALPHA, Mptr( A, Aii, Aoffj+1, Ald, size ), &Ald );
         Aii += tmp1;
         m1  -= tmp1;
      }
/*
*  Return if no more row in the LCM table.
*/
      if( mblks <= 0 ) return;
/*
*  lcmt00 <= upp. The current block owns either diagonals or lower entries.
*  Save the current position in the LCM table. After this column has been
*  completely taken care of, re-start from this row and the next column of
*  the LCM table.
*/
      lcmt  = lcmt00; mblkd = mblks; ioffd = Aoffi;

      mbloc = Amb;
      while( ( mblkd > 0 ) && ( lcmt >= low ) )
      {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
*/
         if( mblkd == 1 ) mbloc = lmbloc;
         pad( C2F_CHAR( UPLO ), C2F_CHAR( &herm ), &mbloc, &nbloc, &lcmt, ALPHA,
              BETA, Mptr( A, ioffd+1, Aoffj+1, Ald, size ), &Ald );
         lcmt00  = lcmt;
         lcmt   -= pmb;
         mblks   = mblkd;
         mblkd--;
         Aoffi   = ioffd;
         ioffd  += mbloc;
      }
/*
*  Set the lower triangular part of sub( A ) when necessary.
*/
      tmp1 = m1 - ioffd + Aii - 1;
      if( lower && ( tmp1 > 0 ) )
         pad( C2F_CHAR( ALL ), C2F_CHAR( &herm ), &tmp1, &nbloc, &izero, ALPHA,
              ALPHA, Mptr( A, ioffd+1, Aoffj+1, Ald, size ), &Ald );

      tmp1    = MIN( Aoffi, iimax ) - Aii + 1;
      m1     -= tmp1;
      n1     -= nbloc;
      lcmt00 += qnb;
      nblks--;
      Aoffj  += nbloc;
/*
*  When the upper triangular part of sub( A ) should be set, take care of the
*  n1 remaining columns of these tmp1 rows immediately.
*/
      if( upper && ( tmp1 > 0 ) && ( n1 > 0 ) )
         pad( C2F_CHAR( ALL ), C2F_CHAR( &herm ), &tmp1, &n1, &izero, ALPHA,
              ALPHA, Mptr( A, Aii, Aoffj+1, Ald, size ), &Ald );
      Aii = Aoffi + 1;
      Ajj = Aoffj + 1;
   }
/*
*  End of PB_Cplapd2
*/
}
