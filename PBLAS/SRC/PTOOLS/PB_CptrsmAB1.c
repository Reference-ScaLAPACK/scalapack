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
void PB_CptrsmAB1( PBTYP_T * TYPE, char * SIDE, char * UPLO, char * TRANSA,
                   char * DIAG, int M, int N, char * ALPHA, char * A, int IA,
                   int JA, int * DESCA,  char * B, int IB, int JB, int * DESCB,
                   char * C, int * DESCC )
#else
void PB_CptrsmAB1( TYPE, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA, JA,
                   DESCA, B, IB, JB, DESCB, C, DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * DIAG, * SIDE, * TRANSA, * UPLO;
   int            IA, IB, JA, JB, M, N;
   char           * ALPHA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCB, * DESCC;
   char           * A, * B, * C;
#endif
{
/*
*  .. Local Scalars ..
*/
   char           * negone, * one;
   int            Acol, Acurcol, Acurrow, Aii, Aimb, Aimb1, Ainb, Ainb1, Ajj,
                  Ald, Almb1, Alnb1, Amb, Amp0, An, Anb, Anq0, Anxtrow, Anxtcol,
                  Arow, Bcol, Bii, Bimb, Binb, Bjj, Bld, Bmb, Bmp0, Bnb, Bnq0,
                  Brow, Cld, ctxt, k=1, kb, kblks, lside, mycol, myrow, npcol,
                  nprow, size, upper;
   MMADD_T        mmadd;
   GERV2D_T       recv;
   GESD2D_T       send;
   GEMM_T         gemm;
   TRSM_T         trsm;
/*
*  .. Local Arrays ..
*/
   char           * Aptr  = NULL, * Aptr0 = NULL, * Bptr = NULL, * Bptr0 = NULL,
                  * Cptr  = NULL;
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
   Ald  = DESCA[LLD_];
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
      Bnq0 = PB_Cnumroc( N, JB, Binb, Bnb, mycol, DESCB[CSRC_], npcol );
      if( Bnq0 <= 0 ) return;

      Bmp0 = PB_Cnumroc( M, IB, Bimb, Bmb, myrow, DESCB[RSRC_], nprow );

      if( !( PB_Cspan( M, IB, Bimb, Bmb, DESCB[RSRC_], nprow ) ) )
      {
         if( Bmp0 > 0 )
         {
            Bptr0 = Mptr( B, Bii, Bjj, Bld, size );
            TYPE->Fmmadd( &M, &Bnq0, TYPE->negone, C, &DESCC[LLD_], ALPHA,
                          Bptr0, &Bld );
            TYPE->Ftrsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( TRANSA ),
                         C2F_CHAR( DIAG ), &M, &Bnq0, TYPE->one, Mptr( A, Aii,
                         Ajj, Ald, size ), &Ald, Bptr0, &Bld );
         }
         return;
      }
      if( Bmp0 > 0 ) Bptr0 = Mptr( B, Bii, Bjj, Bld, size );
   }
   else
   {
      Bmp0 = PB_Cnumroc( M, IB, Bimb, Bmb, myrow, DESCB[RSRC_], nprow );
      if( Bmp0 <= 0 ) return;

      Bnq0 = PB_Cnumroc( N, JB, Binb, Bnb, mycol, DESCB[CSRC_], npcol );

      if( !( PB_Cspan( N, JB, Binb, Bnb, DESCB[CSRC_], npcol ) ) )
      {
         if( Bnq0 > 0 )
         {
            Bptr0 = Mptr( B, Bii, Bjj, Bld, size );
            TYPE->Fmmadd( &Bmp0, &N, TYPE->negone, C, &DESCC[LLD_], ALPHA,
                          Bptr0, &Bld );
            TYPE->Ftrsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( TRANSA ),
                         C2F_CHAR( DIAG ), &Bmp0, &N, TYPE->one, Mptr( A, Aii,
                         Ajj, Ald, size ), &Ald, Bptr0, &Bld );
         }
         return;
      }
      if( Bnq0 > 0 ) Bptr0 = Mptr( B, Bii, Bjj, Bld, size );
   }
/*
*  Handle the general case now
*/
   An     = ( lside ? M : N );
   upper  = ( Mupcase( UPLO[0] ) == CUPPER  );
   negone = TYPE->negone;  one  = TYPE->one;
   recv   = TYPE->Cgerv2d; send = TYPE->Cgesd2d;
   mmadd  = TYPE->Fmmadd;  gemm = TYPE->Fgemm;   trsm = TYPE->Ftrsm;
/*
*  Compute more local information for sub( A )
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

   Cld = DESCC[LLD_];

   if( lside )
   {
      kblks = ( An > Aimb1 ? ( An - Aimb1 - 1 ) / Amb + 2 : 1 );

      if( upper )
      {
         Acurrow = Arow;
         Anxtrow = MModAdd1( Acurrow, nprow );
         Aptr    = Aptr0;
         Bptr    = Bptr0;
         Cptr    = C;

         while( k <= kblks )
         {
            kb  = ( k == 1 ? Aimb1 : ( k == kblks ? Almb1 : Amb ) );
            An -= kb;

            if( myrow == Acurrow )
            {
/*
*  Add contribution of previous blocks of rows of sub( B ) to part of the
*  current block of rows of sub( B )
*/
               mmadd( &kb, &Bnq0, negone, Cptr, &Cld, ALPHA, Bptr, &Bld );
/*
*  Solve updated and current part of block of rows of sub( B )
*/
               trsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( TRANSA ),
                     C2F_CHAR( DIAG ), &kb, &Bnq0, one, Aptr, &Ald, Bptr,
                     &Bld );
/*
*  Add contribution of part of the current block of rows of sub( B ) to the
*  remaining of the contribution of previous blocks of rows of sub( B ). Send
*  this remaining part to next process row.
*/
               if( An > 0 )
               {
                  gemm( C2F_CHAR( TRANSA ), C2F_CHAR( NOTRAN ), &An, &Bnq0, &kb,
                        one, Mptr( Aptr, 0, kb, Ald, size ), &Ald, Bptr, &Bld,
                        one, Mptr( Cptr, kb, 0, Cld, size ), &Cld );
                  send( ctxt, An, Bnq0, Mptr( Cptr, kb, 0, Cld, size ), Cld,
                        Anxtrow, mycol );
               }
               Aptr = Mptr( Aptr, kb, 0, Ald, size );
               Bptr = Mptr( Bptr, kb, 0, Bld, size );
               Cptr = C;
            }
            else if( myrow == Anxtrow )
            {
/*
*  Receive contribution of previous blocks of rows of sub( B ) to be added to
*  next block of rows of sub( B )
*/
               if( An > 0 ) recv( ctxt, An, Bnq0, Cptr, Cld, Acurrow, mycol );
            }

            Aptr     = Mptr( Aptr, 0, kb, Ald, size );
            Acurrow  = Anxtrow;
            Anxtrow  = MModAdd1( Acurrow, nprow );
            k       += 1;
         }
      }
      else
      {
         k       = kblks;
         Acurrow = PB_Cindxg2p( An-1, Aimb1, Amb, Arow, Arow, nprow );
         Anxtrow = MModSub1( Acurrow, nprow );

         while( k > 0 )
         {
            kb  = ( k == 1 ? Aimb1 : ( k == kblks ? Almb1 : Amb ) );
            An -= kb;

            if( myrow == Acurrow )
            {
               Aptr = Mptr( Aptr0, Amp0 - kb, 0, Ald, size );
               Bptr = Mptr( Bptr0, Bmp0 - kb, 0, Bld, size );
               Cptr = Mptr( C,     An,        0, Cld, size );
/*
*  Add contribution of previous blocks of rows of sub( B ) to part of the
*  current block of rows of sub( B )
*/
               mmadd( &kb, &Bnq0, negone, Cptr, &Cld, ALPHA, Bptr, &Bld );
/*
*  Solve updated and current part of block of rows of sub( B )
*/
               trsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( TRANSA ),
                     C2F_CHAR( DIAG ), &kb, &Bnq0, one, Mptr( Aptr, 0, Anq0-kb,
                     Ald, size ), &Ald, Bptr, &Bld );
/*
*  Add contribution of part of the current block of rows of sub( B ) to the
*  remaining of the contribution of previous blocks of rows of sub( B ). Send
*  this remaining part to next process row.
*/
               if( An > 0 )
               {
                  gemm( C2F_CHAR( TRANSA ), C2F_CHAR( NOTRAN ), &An, &Bnq0, &kb,
                        one, Aptr, &Ald, Bptr, &Bld, one, C, &Cld );
                  send( ctxt, An, Bnq0, C, Cld, Anxtrow, mycol );
               }
               Amp0 -= kb;
               Bmp0 -= kb;
            }
            else if( myrow == Anxtrow )
            {
/*
*  Receive contribution of previous blocks of rows of sub( B ) to be added to
*  next block of rows of sub( B )
*/
               if( An > 0 ) recv( ctxt, An, Bnq0, C, Cld, Acurrow, mycol );
            }

            Anq0    -= kb;
            Acurrow  = Anxtrow;
            Anxtrow  = MModSub1( Acurrow, nprow );
            k       -= 1;
         }
      }
   }
   else
   {
      kblks = ( An > Ainb1 ? ( An - Ainb1 - 1 ) / Anb + 2 : 1 );

      if( upper )
      {
         k       = kblks;
         Acurcol = PB_Cindxg2p( An-1, Ainb1, Anb, Acol, Acol, npcol );
         Anxtcol = MModSub1( Acurcol, npcol );

         while( k > 0 )
         {
            kb  = ( k == 1 ? Ainb1 : ( k == kblks ? Alnb1 : Anb ) );
            An -= kb;

            if( mycol == Acurcol )
            {
               Aptr = Mptr( Aptr0, 0, Anq0 - kb, Ald, size );
               Bptr = Mptr( Bptr0, 0, Bnq0 - kb, Bld, size );
               Cptr = Mptr( C,     0, An,        Cld, size );
/*
*  Add contribution of previous blocks of columns of sub( B ) to part of the
*  current block of columns of sub( B )
*/
               mmadd( &Bmp0, &kb, negone, Cptr, &Cld, ALPHA, Bptr, &Bld );
/*
*  Solve updated and current part of block of columns of sub( B )
*/
               trsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( TRANSA ),
                     C2F_CHAR( DIAG ), &Bmp0, &kb, one, Mptr( Aptr, Amp0-kb, 0,
                     Ald, size ), &Ald, Bptr, &Bld );
/*
*  Add contribution of part of the current block of columns of sub( B ) to the
*  remaining of the contribution of previous blocks of columns of sub( B ).
*  Send this remaining part to next process column.
*/
               if( An > 0 )
               {
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRANSA ), &Bmp0, &An, &kb,
                        one, Bptr, &Bld, Aptr, &Ald, one, C, &Cld );
                  send( ctxt, Bmp0, An, C, Cld, myrow, Anxtcol );
               }
               Anq0 -= kb;
               Bnq0 -= kb;
            }
            else if( mycol == Anxtcol )
            {
/*
*  Receive contribution of previous blocks of columns of sub( B ) to be added
*  to next block of columns of sub( B )
*/
               if( An > 0 ) recv( ctxt, Bmp0, An, C, Cld, myrow, Acurcol );
            }

            Amp0   -= kb;
            Acurcol = Anxtcol;
            Anxtcol = MModSub1( Acurcol, npcol );
            k      -= 1;
         }
      }
      else
      {
         Acurcol = Acol;
         Anxtcol = MModAdd1( Acurcol, npcol );
         Aptr    = Aptr0;
         Bptr    = Bptr0;
         Cptr    = C;

         while( k <= kblks )
         {
            kb  = ( k == 1 ? Ainb1 : ( k == kblks ? Alnb1 : Anb ) );
            An -= kb;

            if( mycol == Acurcol )
            {
/*
*  Add contribution of previous blocks of columns of sub( B ) to part of the
*  current block of columns of sub( B )
*/
               mmadd( &Bmp0, &kb, negone, Cptr, &Cld, ALPHA, Bptr, &Bld );
/*
*              Solve updated and current part of block of columns of sub( B )
*/
               trsm( C2F_CHAR( SIDE ), C2F_CHAR( UPLO ), C2F_CHAR( TRANSA ),
                     C2F_CHAR( DIAG ), &Bmp0, &kb, one, Aptr, &Ald, Bptr,
                     &Bld );
/*
*  Add contribution of part of the current block of columns of sub( B ) to the
*  remaining of the contribution of previous blocks of columns of sub( B ).
*  Send this remaining part to next process column.
*/
               if( An > 0 )
               {
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRANSA ), &Bmp0, &An, &kb,
                        one, Bptr, &Bld, Mptr( Aptr, kb, 0, Ald, size ), &Ald,
                        one, Mptr( Cptr, 0, kb, Cld, size ), &Cld );
                  send( ctxt, Bmp0, An, Mptr( Cptr, 0, kb, Cld, size ), Cld,
                        myrow, Anxtcol );
               }
               Aptr = Mptr( Aptr, 0, kb, Ald, size );
               Bptr = Mptr( Bptr, 0, kb, Bld, size );
               Cptr = C;
            }
            else if( mycol == Anxtcol )
            {
/*
*  Receive contribution of previous blocks of columns of sub( B ) to be added
*  to next block of columns of sub( B ).
*/
               if( An > 0 ) recv( ctxt, Bmp0, An, Cptr, Cld, myrow, Acurcol );
            }

            Aptr     = Mptr( Aptr, kb, 0, Ald, size );
            Acurcol  = Anxtcol;
            Anxtcol  = MModAdd1( Acurcol, npcol );
            k       += 1;
         }
      }
   }
/*
*  End of PB_CptrsmAB1
*/
}
