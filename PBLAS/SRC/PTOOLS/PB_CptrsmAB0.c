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
void PB_CptrsmAB0( PBTYP_T * TYPE, char * SIDE, char * UPLO, char * DIAG,
                   Int M, Int N, char * ALPHA, char * A, Int IA, Int JA,
                   Int * DESCA, char * B, Int IB, Int JB, Int * DESCB,
                   char * * C, Int * DESCC, Int * CFREE )
#else
void PB_CptrsmAB0( TYPE, SIDE, UPLO, DIAG, M, N, ALPHA, A, IA, JA, DESCA,
                   B, IB, JB, DESCB, C, DESCC, CFREE )
/*
*  .. Scalar Arguments ..
*/
   char           * DIAG, * SIDE, * UPLO;
   Int            * CFREE, IA, IB, JA, JB, M, N;
   char           * ALPHA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA, * DESCB, * DESCC;
   char           * A, * B, * * C;
#endif
{
/*
*  .. Local Scalars ..
*/
   char           btop, * negone, * one, * talpha, * zero;
   Int            Acol, Acurcol, Acurrow, Aii, Aimb, Aimb1, Ainb, Ainb1, Ajj,
                  Ald, Almb1, Alnb1, Amb, Amp0, Anq0, An, Anb, Arow, Bcol, Bii,
                  Bimb, Bimb1, Binb, Binb1, Bjj, Bld, Bmb, Bmp0, Bnb, Bnq0,
                  Brow, Cld, ctxt, k=1, kb, kblks, kbprev, ktmp, lside, mycol,
                  myrow, npcol, nprow, size, upper;
   char           * Aptr  = NULL, * Aptr0 = NULL, * Bptr = NULL, * Bptr0 = NULL,
                  * Cptr  = NULL;
   MMADD_T        mmadd;
   GEBR2D_T       brecv;
   GEBS2D_T       bsend;
   GEMM_T         gemm;
   TRSM_T         trsm;
/* ..
*  .. Executable Statements ..
*
*/
   size  = TYPE->size;
   lside = ( Mupcase( SIDE[0] ) == CLEFT );
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Retrieve sub( A )'s local information: Aii, Ajj, Arow, Acol, Ald
*/
   Ald = DESCA[LLD_];
   PB_Cinfog2l( IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aii, &Ajj,
                &Arow, &Acol );
/*
*  Retrieve sub( B )'s local information: Bii, Bjj, Brow, Bcol, Bld ...
*/
   Bimb = DESCB[IMB_]; Binb = DESCB[INB_];
   Bmb  = DESCB[MB_ ]; Bnb  = DESCB[NB_ ]; Bld = DESCB[LLD_];
   PB_Cinfog2l( IB, JB, DESCB, nprow, npcol, myrow, mycol, &Bii, &Bjj,
                &Brow, &Bcol );
/*
*  Shorcuts when sub( B ) spans only one process row or column
*/
   if( lside )
   {
      if( !( PB_Cspan( M, IB, Bimb, Bmb, DESCB[RSRC_], nprow ) ) )
      {
         *CFREE = 0;
         Binb1 = PB_Cfirstnb( N, JB, Binb, Bnb );
         PB_Cdescset( DESCC, M, N, M, Binb1, Bmb, Bnb, Brow, Bcol, ctxt, Bld );
         Bnq0  = PB_Cnumroc( N, 0, Binb1, Bnb, mycol, Bcol, npcol );

         if( ( Bnq0 > 0 ) &&
             ( ( ( Brow >= 0 ) && ( myrow == Brow ) ) || ( Brow < 0 ) ) )
         {
            *C = Mptr( B, Bii, Bjj, Bld, size );
            TYPE->Ftrsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                         C2F_CHAR( DIAG ), &M, &Bnq0, ALPHA, Mptr( A, Aii, Ajj,
                         Ald, size ), &Ald, *C, &Bld );
         }
         return;
      }
   }
   else
   {
      if( !( PB_Cspan( N, JB, Binb, Bnb, DESCB[CSRC_], npcol ) ) )
      {
         *CFREE = 0;
         Bimb1 = PB_Cfirstnb( M, IB, Bimb, Bmb );
         PB_Cdescset( DESCC, M, N, Bimb1, N, Bmb, Bnb, Brow, Bcol, ctxt, Bld );
         Bmp0  = PB_Cnumroc( M, 0, Bimb1, Bmb, myrow, Brow, nprow );

         if( ( Bmp0 > 0 ) &&
             ( ( ( Bcol >= 0 ) && ( mycol == Bcol ) ) || ( Bcol < 0 ) ) )
         {
            *C = Mptr( B, Bii, Bjj, Bld, size );
            TYPE->Ftrsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                         C2F_CHAR( DIAG ), &Bmp0, &N, ALPHA, Mptr( A, Aii, Ajj,
                         Ald, size ), &Ald, *C, &Bld );
         }
         return;
      }
   }
/*
*  Handle the general case now
*/
   An     = ( lside ? M : N );
   upper  = ( Mupcase( UPLO[0] ) == CUPPER  );
   talpha = ALPHA;
   negone = TYPE->negone;  one   = TYPE->one;     zero  = TYPE->zero;
   brecv  = TYPE->Cgebr2d; bsend = TYPE->Cgebs2d; mmadd = TYPE->Fmmadd;
   gemm   = TYPE->Fgemm;   trsm  = TYPE->Ftrsm;
/*
*  Compute more local information for sub( A ) and sub( B )
*/
   Aimb  = DESCA[IMB_]; Ainb = DESCA[INB_];
   Amb   = DESCA[MB_ ]; Anb  = DESCA[NB_ ];
   Aimb1 = PB_Cfirstnb( An, IA, Aimb, Amb );
   Almb1 = PB_Clastnb ( An, IA, Aimb, Amb );
   Amp0  = PB_Cnumroc( An, 0, Aimb1, Amb, myrow, Arow, nprow );
   Ainb1 = PB_Cfirstnb( An, JA, Ainb, Anb );
   Alnb1 = PB_Clastnb ( An, JA, Ainb, Anb );
   Anq0  = PB_Cnumroc( An, 0, Ainb1, Anb, mycol, Acol, npcol );
   if( ( Amp0 > 0 ) && ( Anq0 > 0 ) ) Aptr0 = Mptr( A, Aii, Ajj, Ald, size );

   Bimb1 = PB_Cfirstnb( M, IB, Bimb, Bmb );
   Bmp0  = PB_Cnumroc( M, 0, Bimb1, Bmb, myrow, Brow, nprow );
   Binb1 = PB_Cfirstnb( N, JB, Binb, Bnb );
   Bnq0  = PB_Cnumroc( N, 0, Binb1, Bnb, mycol, Bcol, npcol );
   if( ( Bmp0 > 0 ) && ( Bnq0 > 0 ) ) Bptr0 = Mptr( B, Bii, Bjj, Bld, size );

   if( lside )
   {
      Cld = M;
      PB_Cdescset( DESCC, M, N, M, Binb1, Bmb, Bnb, -1, Bcol, ctxt, Cld );
      if( Bnq0 > 0 ) { Cptr = *C = PB_Cmalloc( M * Bnq0 * size ); *CFREE = 1; }
      else           { *C = NULL; *CFREE = 0; return;                         }

      kblks = ( An > Aimb1 ? ( An - Aimb1 - 1 ) / Amb + 2 : 1 );
      btop  = *PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );

      if( upper )
      {
         Acurrow = PB_Cindxg2p( An-1, Aimb1, Amb, Arow, Arow, nprow );
         kb      = Almb1;
         Bptr    = Mptr( Bptr0, Bmp0 - kb, 0, Bld, size );
         Cptr    = Mptr( *C,    An   - kb, 0, Cld, size );
/*
*  Solve last block of rows of sub( B ) and broadcast it vertically to update
*  the rest of sub( B )
*/
         if( myrow == Acurrow )
         {
            trsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                  C2F_CHAR( DIAG ), &kb, &Bnq0, ALPHA, Mptr( Aptr0, Amp0-kb,
                  Anq0-kb, Ald, size ), &Ald, Bptr, &Bld );
            bsend( ctxt, COLUMN, &btop, kb, Bnq0, Bptr, Bld );
            mmadd( &kb, &Bnq0, one, Bptr, &Bld, zero, Cptr, &Cld );
            Amp0 -= kb;
            Bmp0 -= kb;
         }
         else
         {
            brecv( ctxt, COLUMN, &btop, kb, Bnq0, Cptr, Cld, Acurrow, mycol );
         }
         Acurrow = MModSub1( Acurrow, nprow );
         An     -= ( kbprev = kb );
         Anq0   -= kb;
         kblks  -= 1;
/*
*  Lookahead
*/
         while( kblks > 0 )
         {
            kb = ( kblks == 1 ? Aimb1 : Amb );

            Aptr = Mptr( Aptr0, 0,         Anq0, Ald, size );
            Bptr = Mptr( Bptr0, Bmp0 - kb, 0,    Bld, size );
            Cptr = Mptr( *C,    An,        0,    Cld, size );

            if( myrow == Acurrow )
            {
/*
*  Update the current block of rows of sub( B ) with block of rows of sub( B )
*  of previous step
*/
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &kb, &Bnq0,
                     &kbprev, negone, Mptr( Aptr, Amp0-kb, 0, Ald, size ),
                     &Ald, Cptr, &Cld, talpha, Bptr, &Bld );
/*
*  Solve the current block of rows of sub( B )
*/
               trsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                     C2F_CHAR( DIAG ), &kb, &Bnq0, one, Mptr( Aptr, Amp0-kb,
                     -kb, Ald, size ), &Ald, Bptr, &Bld );
/*
*  Broadcast the current block of rows of sub( B ) for next update
*/
               bsend( ctxt, COLUMN, &btop, kb, Bnq0, Bptr, Bld );
               mmadd( &kb, &Bnq0, one, Bptr, &Bld, zero, Mptr( Cptr, -kb, 0,
                      Cld, size ), &Cld );
/*
*  Finish update of the remaining blocks of rows of sub( B ) with block of rows
*  of sub( B ) of previous step
*/
               if( ( ktmp = Amp0 - kb ) > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &ktmp, &Bnq0,
                        &kbprev, negone, Aptr, &Ald, Cptr, &Cld, talpha, Bptr0,
                        &Bld );
               Amp0 -= kb;
               Bmp0 -= kb;
            }
            else
            {
/*
*  Update the remaining rows of sub( B ) with block of rows of sub( B ) of
*  previous step
*/
               if( Amp0 > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Amp0, &Bnq0,
                        &kbprev, negone, Aptr, &Ald, Cptr, &Cld, talpha, Bptr0,
                        &Bld );
/*
*  Receive the current block of rows of sub( B ) for next update
*/
               brecv( ctxt, COLUMN, &btop, kb, Bnq0, Mptr( Cptr, -kb, 0, Cld,
                      size ), Cld, Acurrow, mycol );
            }

            Acurrow = MModSub1( Acurrow, nprow );
            An     -= ( kbprev = kb );
            Anq0   -= kb;
            talpha  = one;
            kblks  -= 1;
         }
      }
      else
      {
         Acurrow = Arow;
         kb      = Aimb1;
/*
*  Solve first block of rows of sub( B ) and broadcast it vertically to update
*  the rest of sub( B )
*/
         if( myrow == Acurrow )
         {
            trsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                  C2F_CHAR( DIAG ), &kb, &Bnq0, ALPHA, Aptr0, &Ald, Bptr0,
                  &Bld );
            bsend( ctxt, COLUMN, &btop, kb, Bnq0, Bptr0, Bld );
            mmadd( &kb, &Bnq0, one, Bptr0, &Bld, zero, Cptr, &Cld );
            Amp0 -= kb;
            Aptr0 = Mptr( Aptr0, kb, 0, Ald, size );
            Bptr0 = Mptr( Bptr0, kb, 0, Bld, size );
         }
         else
         {
            brecv( ctxt, COLUMN, &btop, kb, Bnq0, Cptr, Cld, Acurrow, mycol );
         }
         Acurrow = MModAdd1( Acurrow, nprow );
         kbprev  = kb;
         Cptr    = Mptr( Cptr,  kb, 0, Cld, size );
         Aptr0   = Mptr( Aptr0, 0, kb, Ald, size );
         k      += 1;
/*
*  Lookahead
*/
         while( k <= kblks )
         {
            kb = ( k == kblks ? Almb1 : Amb );

            if( myrow == Acurrow )
            {
/*
*  Update the current block of rows of sub( B ) with block of rows of sub( B )
*  of previous step
*/
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &kb, &Bnq0,
                     &kbprev, negone, Mptr( Aptr0, 0, -kbprev, Ald, size ),
                     &Ald, Mptr( Cptr, -kbprev, 0, Cld, size ), &Cld, talpha,
                     Bptr0, &Bld );
/*
*  Solve the current block of rows of sub( B )
*/
               trsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                     C2F_CHAR( DIAG ), &kb, &Bnq0, one, Aptr0, &Ald, Bptr0,
                     &Bld );
/*
*  Broadcast the current block of rows of sub( B ) for next update
*/
               bsend( ctxt, COLUMN, &btop, kb, Bnq0, Bptr0, Bld );
               mmadd( &kb, &Bnq0, one, Bptr0, &Bld, zero, Cptr, &Cld );
/*
*  Finish update of the remaining blocks of rows of sub( B ) with block of rows
*  of sub( B ) of previous step
*/
               if( ( ktmp = Amp0 - kb ) > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &ktmp, &Bnq0,
                        &kbprev, negone, Mptr( Aptr0, kb, -kbprev, Ald, size ),
                        &Ald, Mptr( Cptr, -kbprev, 0, Cld, size ), &Cld, talpha,
                        Mptr( Bptr0, kb, 0, Bld, size ), &Bld );
               Amp0 -= kb;
               Aptr0 = Mptr( Aptr0, kb, 0, Ald, size );
               Bptr0 = Mptr( Bptr0, kb, 0, Bld, size );
            }
            else
            {
/*
*  Update the remaining rows of sub( B ) with block of rows of sub( B ) of
*  previous step
*/
               if( Amp0 > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Amp0, &Bnq0,
                        &kbprev, negone, Mptr( Aptr0, 0, -kbprev, Ald, size ),
                        &Ald, Mptr( Cptr, -kbprev, 0, Cld, size ), &Cld, talpha,
                        Bptr0, &Bld );
/*
*  Receive the current block of rows of sub( B ) for next update
*/
               brecv( ctxt, COLUMN, &btop, kb, Bnq0, Cptr, Cld, Acurrow,
                      mycol );
            }

            Acurrow = MModAdd1( Acurrow, nprow );
            kbprev  = kb;
            Cptr    = Mptr( Cptr,  kb, 0, Cld, size );
            Aptr0   = Mptr( Aptr0, 0, kb, Ald, size );
            talpha  = one;
            k      += 1;
         }
      }
   }
   else
   {
      Cld = MAX( 1, Bmp0 );
      PB_Cdescset( DESCC, M, N, Bimb1, N, Bmb, Bnb, Brow, -1, ctxt, Cld );
      if( Bmp0 > 0 ) { Cptr = *C = PB_Cmalloc( Bmp0 * N * size ); *CFREE = 1; }
      else           { *C = NULL; *CFREE = 0; return;                         }

      kblks = ( An > Ainb1 ? ( An - Ainb1 - 1 ) / Anb + 2 : 1 );
      btop  = *PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );

      if( upper )
      {
         Acurcol = Acol;
         kb      = Ainb1;
/*
*  Solve first block of columns of sub( B ) and broadcast it horizontally to
*  update the rest of sub( B )
*/
         if( mycol == Acurcol )
         {
            trsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                  C2F_CHAR( DIAG ), &Bmp0, &kb, ALPHA, Aptr0, &Ald, Bptr0,
                  &Bld );
            bsend( ctxt, ROW, &btop, Bmp0, kb, Bptr0, Bld );
            mmadd( &Bmp0, &kb, one, Bptr0, &Bld, zero, Cptr, &Cld );
            Anq0 -= kb;
            Aptr0 = Mptr( Aptr0, 0, kb, Ald, size );
            Bptr0 = Mptr( Bptr0, 0, kb, Bld, size );
         }
         else
         {
            brecv( ctxt, ROW, &btop, Bmp0, kb, Cptr, Cld, myrow, Acurcol );
         }
         Acurcol = MModAdd1( Acurcol, npcol );
         kbprev  = kb;
         k      += 1;
         Cptr    = Mptr( Cptr,   0, kb, Cld, size );
         Aptr0   = Mptr( Aptr0, kb,  0, Ald, size );
/*
*  Lookahead
*/
         while( k <= kblks )
         {
            kb = ( k == kblks ? Alnb1 : Anb );

            if( mycol == Acurcol )
            {
/*
*  Update the current block of columns of sub( B ) with block of columns of
*  sub( B ) of previous step
*/
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp0, &kb,
                     &kbprev, negone, Mptr( Cptr, 0, -kbprev, Cld, size ), &Cld,
                     Mptr( Aptr0, -kbprev, 0, Ald, size ), &Ald, talpha, Bptr0,
                     &Bld );
/*
*  Solve the current block of columns of sub( B )
*/
               trsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                     C2F_CHAR( DIAG ), &Bmp0, &kb, one, Aptr0, &Ald, Bptr0,
                     &Bld );
/*
*  Broadcast the current block of columns of sub( B ) for next update
*/
               bsend( ctxt, ROW, &btop, Bmp0, kb, Bptr0, Bld );
               mmadd( &Bmp0, &kb, one, Bptr0, &Bld, zero, Cptr, &Cld );
/*
*  Finish update of the remaining blocks of columns of sub( B ) with block of
*  columns of sub( B ) of previous step
*/
               if( ( ktmp = Anq0 - kb ) > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp0, &ktmp,
                        &kbprev, negone, Mptr( Cptr, 0, -kbprev, Cld, size ),
                        &Cld, Mptr( Aptr0, -kbprev, kb, Ald, size ), &Ald,
                        talpha, Mptr( Bptr0, 0, kb, Bld, size ), &Bld );
               Anq0 -= kb;
               Aptr0 = Mptr( Aptr0, 0, kb, Ald, size );
               Bptr0 = Mptr( Bptr0, 0, kb, Bld, size );
            }
            else
            {
/*
*  Update the remaining columns of sub( B ) with block of columns of sub( B )
*  of previous step
*/
               if( Anq0 > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp0, &Anq0,
                        &kbprev, negone, Mptr( Cptr, 0, -kbprev, Cld, size ),
                        &Cld, Mptr( Aptr0, -kbprev, 0, Ald, size ), &Ald,
                        talpha, Bptr0, &Bld );
/*
*  Receive the current block of columns of sub( B ) for next update
*/
               brecv( ctxt, ROW, &btop, Bmp0, kb, Cptr, Cld, myrow, Acurcol );
            }

            Acurcol = MModAdd1( Acurcol, npcol );
            kbprev  = kb;
            Cptr    = Mptr( Cptr,   0, kb, Cld, size );
            Aptr0   = Mptr( Aptr0, kb,  0, Ald, size );
            talpha  = one;
            k      += 1;
         }
      }
      else
      {
         Acurcol = PB_Cindxg2p( An-1, Ainb1, Anb, Acol, Acol, npcol );
         kb      = Alnb1;
         Bptr    = Mptr( Bptr0, 0, Bnq0 - kb, Bld, size );
         Cptr    = Mptr( *C,    0, An   - kb, Cld, size );
/*
*  Solve last block of columns of sub( B ) and broadcast it horizontally to
*  update the rest of sub( B )
*/
         if( mycol == Acurcol )
         {
            trsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                  C2F_CHAR( DIAG ), &Bmp0, &kb, ALPHA, Mptr( Aptr0, Amp0-kb,
                  Anq0-kb, Ald, size ), &Ald, Bptr, &Bld );
            bsend( ctxt, ROW, &btop, Bmp0, kb, Bptr, Bld );
            mmadd( &Bmp0, &kb, one, Bptr, &Bld, zero, Cptr, &Cld );
            Anq0 -= kb;
            Bnq0 -= kb;
         }
         else
         {
            brecv( ctxt, ROW, &btop, Bmp0, kb, Cptr, Cld, myrow, Acurcol );
         }
         Acurcol = MModSub1( Acurcol, npcol );
         An     -= ( kbprev = kb );
         Amp0   -= kb;
         kblks  -= 1;
/*
*  Lookahead
*/
         while( kblks > 0 )
         {
            kb = ( kblks == 1 ? Ainb1 : Anb );

            Aptr = Mptr( Aptr0, Amp0, 0,         Ald, size );
            Bptr = Mptr( Bptr0, 0,    Bnq0 - kb, Bld, size );
            Cptr = Mptr( *C,    0,    An,        Cld, size );

            if( mycol == Acurcol )
            {
/*
*  Update the current block of columns of sub( B ) with block of columns of
*  sub( B ) of previous step
*/
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp0, &kb,
                     &kbprev, negone, Cptr, &Cld, Mptr( Aptr, 0, Anq0-kb, Ald,
                     size ), &Ald, talpha, Bptr, &Bld );
/*
*  Solve the current block of columns of sub( B )
*/
               trsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                     C2F_CHAR( DIAG ), &Bmp0, &kb, one, Mptr( Aptr, -kb,
                     Anq0-kb, Ald, size ), &Ald, Bptr, &Bld );
/*
*  Broadcast the current block of columns of sub( B ) for next update
*/
               bsend( ctxt, ROW, &btop, Bmp0, kb, Bptr, Bld );
               mmadd( &Bmp0, &kb, one, Bptr, &Bld, zero, Mptr( Cptr, 0, -kb,
                      Cld, size ), &Cld );
/*
*  Finish update of the remaining blocks of columns of sub( B ) with block of
*  columns of sub( B ) of previous step
*/
               if( ( ktmp = Anq0 - kb ) > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp0, &ktmp,
                        &kbprev, negone, Cptr, &Cld, Aptr, &Ald, talpha, Bptr0,
                        &Bld );
               Anq0 -= kb;
               Bnq0 -= kb;
            }
            else
            {
/*
*  Update the remaining columns of sub( B ) with block of columns of sub( B )
*  of previous step
*/
               if( Anq0 > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp0, &Anq0,
                        &kbprev, negone, Cptr, &Cld, Aptr, &Ald, talpha, Bptr0,
                        &Bld );
/*
*  Receive the current block of columns of sub( B ) for next update
*/
               brecv( ctxt, ROW, &btop, Bmp0, kb, Mptr( Cptr, 0, -kb, Cld,
                      size ), Cld, myrow, Acurcol );
            }

            Acurcol = MModSub1( Acurcol, npcol );
            An     -= ( kbprev = kb );
            Amp0   -= kb;
            talpha  = one;
            kblks  -= 1;
         }
      }
   }
/*
*  End of PB_CptrsmAB0
*/
}
