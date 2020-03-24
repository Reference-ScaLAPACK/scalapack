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
void PB_CScatterV( PBTYP_T * TYPE, char * DIRECA, Int M, Int N,
                   char * A, Int IA, Int JA, Int * DESCA, char * AROC,
                   char * ALPHA, char * B, Int IB, Int JB, Int * DESCB,
                   char * BROC )
#else
void PB_CScatterV( TYPE, DIRECA, M, N, A, IA, JA, DESCA, AROC,
                   ALPHA, B, IB, JB, DESCB, BROC )
/*
*  .. Scalar Arguments ..
*/
   char           * ALPHA, * AROC, * BROC, * DIRECA;
   Int            IA, IB, JA, JB, M, N;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA, * DESCB;
   char           * A, * B;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CScatterV disaggregates the one-dimensional submatrix sub( A ) de-
*  noting  A( IA:IA+M-1, JA:JA+N-1 )  into a  two-dimensional  submatrix
*  sub( B ) denoting B( IB:IB+M-1, JB:JB+N-1 )  when  AROC  is  equal to
*  BROC and B( IB:IB+N-1, JB:JB+M-1 ) otherwise:
*
*     sub( B ) := alpha * sub( B ) + sub( A ).
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
*  TYPE    (local input) pointer to a PBTYP_T structure
*          On entry,  TYPE  is a pointer to a structure of type PBTYP_T,
*          that contains type information (See pblas.h).
*
*  DIRECA  (global input) pointer to CHAR
*          On entry,  DIRECA  specifies  the direction in which the rows
*          or columns of sub( A ) should be disaggregated as follows:
*             DIRECA = 'F' or 'f'   forward  or increasing,
*             DIRECA = 'B' or 'b'   backward or decreasing.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( A ). M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N specifies the number of columns of the  submatrix
*          sub( A ). N must be at least zero.
*
*  A       (local input) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where LLD_A
*          is   DESCA[LLD_], i.e. at least  MAX( 1, Lr( M, IA ) ),  and,
*          Ka is at least Lc( N, JA ). Before entry, this array contains
*          the local entries of the matrix A.
*
*  IA      (global input) INTEGER
*          On entry, IA  specifies A's global row index, which points to
*          the beginning of the submatrix sub( A ).
*
*  JA      (global input) INTEGER
*          On entry, JA  specifies A's global column index, which points
*          to the beginning of the submatrix sub( A ).
*
*  DESCA   (global and local input) INTEGER array
*          On entry, DESCA  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix A.
*
*  AROC    (global input) pointer to CHAR
*          On entry,  AROC  specifies  the  orientation of the submatrix
*          sub( A ). When AROC is 'R' or 'r',  sub( A ) is a row matrix,
*          and a column matrix otherwise.
*
*  ALPHA   (local input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.
*
*  B       (local output) pointer to CHAR
*          On entry, A is an array of dimension (LLD_B, Kb), where LLD_B
*          is DESCB[LLD_], i.e. at least MAX( 1, Lr( M, IB ) ) when AROC
*          and BROC are equal, and MAX( 1, Lr( N, IB ) ) otherwise, and,
*          Kb is at least Lc( N, JB )  when AROC and BROC are equal, and
*          Lc( M, JB ) otherwise. On exit, this array contains the local
*          entries of the disaggregated submatrix sub( A ).
*
*  IB      (global input) INTEGER
*          On entry, IB  specifies B's global row index, which points to
*          the beginning of the submatrix sub( B ).
*
*  JB      (global input) INTEGER
*          On entry, JB  specifies B's global column index, which points
*          to the beginning of the submatrix sub( B ).
*
*  DESCB   (global and local input) INTEGER array
*          On entry, DESCB  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix B.
*
*  BROC    (global input) pointer to CHAR
*          On entry,  BROC  specifies  the  orientation of the submatrix
*          sub( B ). When BROC is 'R' or 'r',  sub( B ) is a row matrix,
*          and a column matrix otherwise.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           * one;
   Int            Afwd, Bbufld, Bcol, Bcurcol, Bcurrow, Bii, Bimb, Bimb1, Binb,
                  Binb1, BisRow, Bjj, Bld, Bm, Bmb, Bmp, Bn, Bnb, Bnnxt, BnnxtL,
                  Bnpre, Bnq, Brow, WAfr, ctxt, kb, mycol, mydist, mydistnb,
                  myrow, nlen, npcol, nprow, offset, size, srcdist, stride,
                  stridenb, tmp;
   MMADD_T        add;
   MMSHFT_T       shft;
/*
*  .. Local Arrays ..
*/
   Int            Bd0[DLEN_], WAd[DLEN_];
   char           * Bptr = NULL, * Bbuf = NULL, * Bbufptr = NULL, * WA = NULL;
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
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   Afwd   = ( Mupcase( DIRECA[0] ) == CFORWARD );
   BisRow = ( Mupcase( BROC  [0] ) ==     CROW );
   if( Mupcase( AROC[0] ) == Mupcase( BROC[0] ) ) { Bm = M; Bn = N; }
   else                                           { Bm = N; Bn = M; }
/*
*  Retrieve sub( B )'s local information: Bii, Bjj, Brow, Bcol ...
*/
   PB_Cinfog2l( IB, JB, DESCB, nprow, npcol, myrow, mycol, &Bii, &Bjj,
                &Brow, &Bcol );
   Bimb  = DESCB[IMB_]; Binb = DESCB[INB_]; Bmb = DESCB[MB_]; Bnb  = DESCB[NB_];
   Bimb1 = PB_Cfirstnb( Bm, IB, Bimb, Bmb );
   Bmp   = PB_Cnumroc( Bm, 0, Bimb1, Bmb, myrow, Brow, nprow );
   Binb1 = PB_Cfirstnb( Bn, JB, Binb, Bnb );
   Bnq   = PB_Cnumroc( Bn, 0, Binb1, Bnb, mycol, Bcol, npcol );
   Bld  = DESCB[LLD_]; size = TYPE->size; one = TYPE->one;
   if( ( Bmp > 0 ) && ( Bnq > 0 ) ) Bptr = Mptr( B, Bii, Bjj, Bld, size );

   if( BisRow )
   {
/*
*  Compute descriptor Bd0 for sub( B ).
*/
      if( Afwd ) { Bcurrow = Brow; }
      else { Bcurrow = PB_Cindxg2p( Bm-1, Bimb1, Bmb, Brow, Brow, nprow ); }
      PB_Cdescset( Bd0, Bm, Bn, Bm, Binb1, Bmb, Bnb, Bcurrow, Bcol, ctxt, Bld );
/*
*  Align sub( A ) with sub( B )
*/
      PB_CInV( TYPE, NOCONJG, BROC, Bm, Bn, Bd0, Bm, A, IA, JA, DESCA, AROC,
               &WA, WAd, &WAfr );
/*
*  Disaggregate WA = sub( A )
*/
      if( ( Brow == -1 ) || ( nprow == 1 ) )
      {
/*
*  sub( B ) is replicated
*/
         if( Bnq > 0 )
            TYPE->Fmmadd( &Bm, &Bnq, one, WA, &WAd[LLD_], ALPHA, Bptr, &Bld );
         if( WAfr ) free( WA );
         return;
      }

      if( !( PB_Cspan( Bm, 0, Bimb1, Bmb, Brow, nprow ) ) )
      {
/*
*  sub( B ) spans only one process row
*/
         if( ( myrow == Brow ) && ( Bnq > 0 ) )
            TYPE->Fmmadd( &Bm, &Bnq, one, WA, &WAd[LLD_], ALPHA, Bptr, &Bld );
         if( WAfr ) free( WA );
         return;
      }
/*
*  sub( B ) spans more than one process row
*/
      if( Afwd )
      {
/*
*  sub( B ) is not replicated and spans more than one process row. Forward row
*  dissagregation starts in the process row where the global row IB resides.
*/
         if( ( Bmp > 0 ) && ( Bnq > 0 ) )
         {
/*
*  Compute how may rows are before and after me (Bnpre and Bnnxt).
*/
            Bnpre = PB_Cnpreroc( Bm, 0, Bimb1, Bmb, myrow, Brow, nprow );
            Bnnxt = PB_Cnnxtroc( Bm, 0, Bimb1, Bmb, myrow, Brow, nprow );
            nlen  = Bmp + Bnnxt;

            if( Bnpre > 0 )
            {
/*
*  If I don't own the row IB, then allocate and receive a buffer of length
*  ( Bmp + Bnnxt ) * Bnq from the previous process row.
*/
               Bbufptr = Bbuf = PB_Cmalloc( nlen * Bnq * size );
               Bbufld  = nlen;
               TYPE->Cgerv2d( ctxt, nlen, Bnq, Bbuf, Bbufld, MModSub1( myrow,
                              nprow ), mycol );
               kb      = Bmb;
            }
            else
            {
/*
*  Otherwise, reuse WA.
*/
               Bbufptr = Bbuf = WA;
               Bbufld  = WAd[LLD_];
               kb      = Bimb1;
            }
/*
*  Unpack the received data
*/
            if( Bnnxt > 0 )
            {
/*
*  If some rows reside in the process row following mine, then unpack my piece,
*  sort the buffer and send those Bnnxt rows to the next process row.
*/
               add      = TYPE->Fmmadd; shft = TYPE->Frshft;
               mydistnb = ( nprow - MModSub( myrow, Brow, nprow ) - 1 );
               stride   = ( mydistnb *= Bmb ) * size;

               do
               {
                  kb       = MIN( kb, nlen );
                  add( &kb, &Bnq, one, Bbufptr, &Bbufld, ALPHA, Bptr, &Bld );
                  nlen    -= kb;
                  offset   = -kb;
                  shft( &nlen, &Bnq, &offset, Bbufptr, &Bbufld );
                  Bptr    += kb*size;
                  Bbufptr += stride;
                  nlen    -= mydistnb;
                  kb       = Bmb;
               } while( nlen > 0 );
/*
*  send buffer of length Bnnxt * Bnq to the next process row.
*/
               TYPE->Cgesd2d( ctxt, Bnnxt, Bnq, Bbuf, Bbufld, MModAdd1( myrow,
                              nprow ), mycol );
            }
            else
            {
/*
*  Otherwise, I must be the last process involved in the operation, so no
*  unpacking is necessary.
*/
               TYPE->Fmmadd( &Bmp, &Bnq, one, Bbufptr, &Bbufld, ALPHA, Bptr,
                             &Bld );
            }
/*
*  If I don't own the row IB, then release the dynamically allocated buffer.
*/
            if( Bnpre > 0 ) free( Bbuf );
         }
         if( WAfr ) free( WA );
      }
      else
      {
         if( ( Bmp > 0 ) && ( Bnq > 0 ) )
         {
/*
*  Compute how may rows are before and after me (Bnpre, Bnnxt).
*/
            Bnnxt  = PB_Cnnxtroc( Bm, 0, Bimb1, Bmb, myrow,   Brow, nprow );
            BnnxtL = PB_Cnnxtroc( Bm, 0, Bimb1, Bmb, Bcurrow, Brow, nprow );
            Bnnxt  = MModSub( Bnnxt, BnnxtL, Bm );
            Bnpre  = ( nlen = Bm - Bnnxt ) - Bmp;

            if( Bnnxt > 0 )
            {
/*
*  If I don't own the row IB+M-1, then allocate and receive a buffer of length
*  ( Bm - Bnnxt ) * Bnq from the next process row.
*/
               Bbufptr = Bbuf = PB_Cmalloc( nlen * Bnq * size );
               Bbufld  = nlen;
               TYPE->Cgerv2d( ctxt, nlen, Bnq, Bbuf, Bbufld, MModAdd1( myrow,
                              nprow ), mycol );
            }
            else
            {
/*
*  Otherwise, reuse WA.
*/
               Bbufptr = Bbuf = WA;
               Bbufld = WAd[LLD_];
            }
/*
*  Unpack the received data
*/
            if( Bnpre > 0 )
            {
/*
*  If some rows reside in the process row preceeding mine, then unpack my piece,
*  sort the buffer and send those Bnpre rows to the previous process row.
*/
               add      = TYPE->Fmmadd; shft = TYPE->Frshft;
               mydist   = MModSub( Bcurrow, myrow, nprow );
               srcdist  = MModSub( Bcurrow, Brow,  nprow );
               stridenb = ( nprow - mydist - 1 ) * Bmb;

               if( mydist < srcdist )
               {
                  tmp      = ( Bimb1 + ( srcdist - mydist - 1 ) * Bmb );
                  Bbufptr += tmp * size;
                  nlen    -= tmp;
                  kb       = Bmb;
               }
               else if( mydist == srcdist )
               {
                  kb       = Bimb1;
               }
               else
               {
                  Bbufptr += stridenb * size;
                  nlen    -= stridenb;
                  kb       = Bmb;
               }

               do
               {
                  kb       = MIN( kb, nlen );
                  add( &kb, &Bnq, one, Bbufptr, &Bbufld, ALPHA, Bptr, &Bld );
                  nlen    -= kb;
                  offset   = -kb;
                  shft( &nlen, &Bnq, &offset, Bbufptr, &Bbufld );
                  Bptr    += kb*size;
                  Bbufptr += stridenb*size;
                  nlen    -= stridenb;
                  kb       = Bmb;
               } while( nlen > 0 );
/*
*  send buffer of length Bnpre * Bnq to the previous process row.
*/
               TYPE->Cgesd2d( ctxt, Bnpre, Bnq, Bbuf, Bbufld, MModSub1( myrow,
                              nprow ), mycol );
            }
            else
            {
/*
*  Otherwise, I must be the last process involved in the operation, so no
*  unpacking is necessary.
*/
               TYPE->Fmmadd( &Bmp, &Bnq, one, Bbufptr, &Bbufld, ALPHA, Bptr,
                             &Bld );
            }
/*
*  If I don't own the row IB+M-1, then release the dynamically allocated buffer.
*/
            if( Bnnxt > 0 ) free( Bbuf );
         }
         if( WAfr ) free( WA );
      }
   }
   else
   {
/*
*  Compute descriptor Bd0 for sub( B ).
*/
      if( Afwd ) { Bcurcol = Bcol; }
      else { Bcurcol = PB_Cindxg2p( Bn-1, Binb1, Bnb, Bcol, Bcol, npcol ); }
      PB_Cdescset( Bd0, Bm, Bn, Bimb1, Bn, Bmb, Bnb, Brow, Bcurcol, ctxt, Bld );
/*
*  Align sub( A ) with sub( B )
*/
      PB_CInV( TYPE, NOCONJG, BROC, Bm, Bn, Bd0, Bn, A, IA, JA, DESCA, AROC,
               &WA, WAd, &WAfr );
/*
*     Disaggregate WA = sub( A )
*/
      if( ( Bcol == -1 ) || ( npcol == 1 ) )
      {
/*
*  sub( B ) is replicated
*/
         if( Bmp > 0 )
            TYPE->Fmmadd( &Bmp, &Bn, one, WA, &WAd[LLD_], ALPHA, Bptr, &Bld );
         if( WAfr ) free( WA );
         return;
      }

      if( !( PB_Cspan( Bn, 0, Binb1, Bnb, Bcol, npcol ) ) )
      {
/*
*  sub( B ) spans only one process column
*/
         if( ( mycol == Bcol ) && ( Bmp > 0 ) )
            TYPE->Fmmadd( &Bmp, &Bn, one, WA, &WAd[LLD_], ALPHA, Bptr, &Bld );
         if( WAfr ) free( WA );
         return;
      }
/*
*  sub( B ) spans more than one process column
*/
      if( Afwd )
      {
/*
*  sub( B ) is not replicated and spans more than one process column. Forward
*  column dissagregation starts in the process column where the global column
*  JB resides.
*/
         if( ( Bmp > 0 ) && ( Bnq > 0 ) )
         {
/*
*  Compute how may columns are before and after me (Bnpre and Bnnxt).
*/
            Bnpre = PB_Cnpreroc( Bn, 0, Binb1, Bnb, mycol, Bcol, npcol );
            Bnnxt = PB_Cnnxtroc( Bn, 0, Binb1, Bnb, mycol, Bcol, npcol );
            nlen  = Bnq + Bnnxt;

            if( Bnpre > 0 )
            {
/*
*  If I don't own the column JB, then allocate and receive a buffer of length
*  Bmp * ( Bnq + Bnnxt ) from the previous process column.
*/
               Bbufptr = Bbuf = PB_Cmalloc( Bmp * nlen * size );
               Bbufld  = Bmp;
               TYPE->Cgerv2d( ctxt, Bmp, nlen, Bbuf, Bbufld, myrow,
                              MModSub1( mycol, npcol ) );
               kb = Bnb;
            }
            else
            {
/*
*  Otherwise, reuse WA.
*/
               Bbufptr = Bbuf = WA;
               Bbufld  = WAd[LLD_];
               kb      = Binb1;
            }
/*
*  Unpack the received data
*/
            if( Bnnxt > 0 )
            {
/*
*  If some columns reside in the process column following mine, then unpack my
*  piece, sort the buffer and send those Bnnxt columns to the next process
*  column.
*/
               add      = TYPE->Fmmadd; shft = TYPE->Fcshft;
               mydistnb = ( npcol - MModSub( mycol, Bcol, npcol ) - 1 );
               stride   = ( mydistnb *= Bnb ) * Bbufld * size;

               do
               {
                  kb       = MIN( kb, nlen );
                  add( &Bmp, &kb, one, Bbufptr, &Bbufld, ALPHA, Bptr, &Bld );
                  nlen    -= kb;
                  offset   = -kb;
                  shft( &Bmp, &nlen, &offset, Bbufptr, &Bbufld );
                  Bptr    += kb*Bld*size;
                  Bbufptr += stride;
                  nlen    -= mydistnb;
                  kb       = Bnb;
               } while( nlen > 0 );
/*
*  send buffer of length Bmp * Bnnxt to the next process column.
*/
               TYPE->Cgesd2d( ctxt, Bmp, Bnnxt, Bbuf, Bbufld, myrow,
                              MModAdd1( mycol, npcol ) );
            }
            else
            {
/*
*  Otherwise, I must be the last process involved in the operation, so no
*  unpacking is necessary.
*/
               TYPE->Fmmadd( &Bmp, &Bnq, one, Bbufptr, &Bbufld, ALPHA, Bptr,
                             &Bld );
            }
/*
*  If I don't own the column JB, then release the dynamically allocated buffer.
*/
            if( Bnpre > 0 ) free( Bbuf );
         }
         if( WAfr ) free( WA );
      }
      else
      {
         if( ( Bmp > 0 ) && ( Bnq > 0 ) )
         {
/*
*  Compute how may rows are before and after me (Bnpre, Bnnxt).
*/
            Bnnxt  = PB_Cnnxtroc( Bn, 0, Binb1, Bnb, mycol,   Bcol, npcol );
            BnnxtL = PB_Cnnxtroc( Bn, 0, Binb1, Bnb, Bcurcol, Bcol, npcol );
            Bnnxt  = MModSub( Bnnxt, BnnxtL, Bn );
            Bnpre  = ( nlen = Bn - Bnnxt ) - Bnq;

            if( Bnnxt > 0 )
            {
/*
*  If I don't own the column JB+N-1, then allocate and receive a buffer of
*  length Bmp * ( Bn - Bnnxt ) from the next process column.
*/
               Bbufptr = Bbuf = PB_Cmalloc( Bmp * nlen * size );
               Bbufld  = Bmp;
               TYPE->Cgerv2d( ctxt, Bmp, nlen, Bbuf, Bbufld, myrow,
                              MModAdd1( mycol, npcol ) );
            }
            else
            {
/*
*  Otherwise, reuse WA.
*/
               Bbufptr = Bbuf = WA;
               Bbufld  = WAd[LLD_];
            }
/*
*  Unpack the received data
*/
            if( Bnpre > 0 )
            {
/*
*  If some columns reside in the process column preceeding mine, then unpack my
*  piece, sort the buffer and send those Bnpre columns to the previous process
*  column.
*/
               add      = TYPE->Fmmadd; shft = TYPE->Fcshft;
               mydist   = MModSub( Bcurcol, mycol, npcol );
               srcdist  = MModSub( Bcurcol, Bcol,  npcol );
               stridenb = ( npcol - mydist - 1 ) * Bnb;

               if( mydist < srcdist )
               {
                  tmp      = ( Binb1 + ( srcdist - mydist - 1 ) * Bnb );
                  Bbufptr += tmp * Bbufld * size;
                  nlen    -= tmp;
                  kb       = Bnb;
               }
               else if( mydist == srcdist )
               {
                  kb       = Binb1;
               }
               else
               {
                  Bbufptr += stridenb * Bbufld * size;
                  nlen    -= stridenb;
                  kb       = Bnb;
               }

               do
               {
                  kb = MIN( kb, nlen );
                  add( &Bmp, &kb, one, Bbufptr, &Bbufld, ALPHA, Bptr, &Bld );
                  nlen    -= kb;
                  offset   = -kb;
                  shft( &Bmp, &nlen, &offset, Bbufptr, &Bbufld );
                  Bptr    += kb * Bld * size;
                  Bbufptr += stridenb * Bbufld * size;
                  nlen    -= stridenb;
                  kb       = Bnb;
               } while( nlen > 0 );
/*
*  send buffer of length Bmp * Bnpre to the previous process column.
*/
               TYPE->Cgesd2d( ctxt, Bmp, Bnpre, Bbuf, Bbufld, myrow,
                              MModSub1( mycol, npcol ) );
            }
            else
            {
/*
*  Otherwise, I must be the last process involved in the operation, so no
*  unpacking is necessary.
*/
               TYPE->Fmmadd( &Bmp, &Bnq, one, Bbufptr, &Bbufld, ALPHA, Bptr,
                             &Bld );
            }
/*
*  If I don't own the column JB+N-1, then release the dynamically allocated
*  buffer.
*/
            if( Bnnxt > 0 ) free( Bbuf );
         }
         if( WAfr ) free( WA );
      }
   }
/*
*  End of PB_CScatterV
*/
}
