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
void PB_CGatherV( PBTYP_T * TYPE, char * ALLOC, char * DIRECA, int M, int N,
                  char * A, int IA, int JA, int * DESCA, char * AROC,
                  char * * B, int * DESCB, int * BFREE )
#else
void PB_CGatherV( TYPE, ALLOC, DIRECA, M, N, A, IA, JA, DESCA, AROC, B,
                  DESCB, BFREE )
/*
*  .. Scalar Arguments ..
*/
   char           * ALLOC, * AROC, * DIRECA;
   int            * BFREE, IA, JA, M, N;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCB;
   char           * A, * * B;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CGatherV aggregates a submatrix sub( A ) = A(IA:IA+M-1, JA:JA+N-1)
*  into a one-dimensional multivector B. The submatrix sub( A ) is  spe-
*  cified on input to  this routine that is reused whenever possible. On
*  return, the one-dimensional multivector is specified by a pointer  to
*  some data, a descriptor array describing its layout and a logical va-
*  lue indicating if this local piece of data has been dynamically allo-
*  cated by this function.
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
*  ALLOC   (global input) pointer to CHAR
*          On entry, ALLOC specifies if data  should  be  allocated even
*          when unnecessary as follows:
*             ALLOC = 'A' or 'a'    data allocation is enforced,
*             ALLOC = 'R' or 'r'    data is reused when possible.
*
*  DIRECA  (global input) pointer to CHAR
*          On entry,  DIRECA  specifies  the direction in which the rows
*          or columns of sub( A ) should be aggregated as follows:
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
*  B       (local output) pointer to pointer to CHAR
*          On exit, * B is an array containing the aggregated  submatrix
*          sub( A ).
*
*  DESCB   (global and local output) INTEGER array
*          On exit, DESCB is a descriptor array of dimension DLEN_  des-
*          cribing the data layout of the data pointed to by * B.
*
*  BFREE   (local output) INTEGER
*          On exit,  BFREE  specifies  if  it has been possible to reuse
*          the submatrix sub( A ), i.e., if some dynamic  memory was al-
*          located for the data pointed to by *B or not. When  BFREE  is
*          zero, no dynamic memory was allocated.  Otherwise, some dyna-
*          mic  memory  was allocated by this function that one MUST re-
*          lease as soon as possible.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           * one, * zero;
   int            Afwd, AggRow, AiiD, AiiR, Ainb1D, Ainb1R, Ald, AmyprocD,
                  AmyprocR, AnR, AnbD, AnbR, AnnxtL, AnnxtR, AnpD, AnpR, AnpreR,
                  AnprocsR, ArocR, AsrcD, AsrcR, Bld, Bsrc_, ctxt, k, kb, kblks,
                  kn, ktmp, mycol, mydist, mydistnb, myrow, nlen, npcol, nprow,
                  offset, size, srcdist;
   MMADD_T        add;
   MMSHFT_T       shft;
/*
*  .. Local Arrays ..
*/
   char           * Aptr = NULL, * Bptr = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Initialize the output parameters to a default value
*/
   *BFREE = 0;
   *B     = NULL;
/*
*  Quick return if possible
*/
   if( ( M <= 0 ) || ( N <= 0 ) )
   {
      PB_Cdescset( DESCB, M, N, DESCA[IMB_], DESCA[INB_], DESCA[MB_],
                   DESCA[NB_], DESCA[RSRC_], DESCA[CSRC_], DESCA[CTXT_], 1 );
      return;
   }
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   if( ( AggRow = ( Mupcase( AROC[0] ) == CROW ) ) != 0 )
   {
/*
*  Accumulate rows of sub( A )
*/
      AnbR = DESCA[MB_]; AnbD = DESCA[NB_]; Ald = DESCA[LLD_];
      PB_Cinfog2l( IA, JA, DESCA, nprow, npcol, myrow, mycol, &AiiR, &AiiD,
                   &AsrcR, &AsrcD );
      Ainb1D = PB_Cfirstnb( N, JA, DESCA[INB_], AnbD );
      AnpD   = PB_Cnumroc( N, 0, Ainb1D, AnbD, mycol, AsrcD, npcol );
/*
*  If sub( A ) is either replicated or spans only one process row, no data needs
*  to be exchanged by the processes, the operation is purely local.
*/
      if( !( PB_Cspan( M, IA, DESCA[IMB_], AnbR, AsrcR, nprow ) ) )
      {
         if( Mupcase( ALLOC[0] ) == CREUSE )
         {
/*
*  sub( A ) can be reused
*/
            if( ( ( myrow == AsrcR ) || ( AsrcR < 0 ) ) && ( AnpD > 0 ) )
            {
/*
*  If I own some entries of sub( A ), set *B
*/
               Bld = Ald;
               *B  = Mptr( A, AiiR, AiiD, Ald, TYPE->size );
            }
            else { Bld = 1; }
         }
         else
         {
/*
*  sub( A ) cannot be reused, make a copy of it.
*/
            if( ( ( myrow == AsrcR ) || ( AsrcR < 0 ) ) && ( AnpD > 0 ) )
            {
/*
*  If I own some entries of sub( A ), allocate space for the copy, and copy the
*  data.
*/
               Bld = M;
               if( AnpD > 0 )
               {
                  size   = TYPE->size;
                  *B     = PB_Cmalloc( AnpD * M * size );
                  *BFREE = 1;
                  TYPE->Fmmadd( &M, &AnpD, TYPE->one, Mptr( A, AiiR, AiiD, Ald,
                                size ), &Ald, TYPE->zero, *B, &Bld );
               }
            }
            else { Bld = 1; }
         }
/*
*  Describe the resulting operand
*/
         PB_Cdescset( DESCB, M, N, M, Ainb1D, AnbR, AnbD, AsrcR, AsrcD, ctxt,
                      Bld );
         return;
      }

      AnR      = M;     Bsrc_ = RSRC_;
      AmyprocR = myrow; AmyprocD = mycol; AnprocsR = nprow;
      Ainb1R   = PB_Cfirstnb( M, IA, DESCA[IMB_], AnbR );
      AnpR     = PB_Cnumroc( M, 0, Ainb1R, AnbR, myrow, AsrcR, nprow );
   }
   else
   {
/*
*  Accumulate columns of sub( A )
*/
      AnbD = DESCA[MB_ ]; AnbR = DESCA[NB_ ]; Ald = DESCA[LLD_];
      PB_Cinfog2l( IA, JA, DESCA, nprow, npcol, myrow, mycol, &AiiD, &AiiR,
                   &AsrcD, &AsrcR );
      Ainb1D = PB_Cfirstnb( M, IA, DESCA[IMB_], AnbD );
      AnpD   = PB_Cnumroc( M, 0, Ainb1D, AnbD, myrow, AsrcD, nprow );
/*
*  If sub( A ) is either replicated or spans only one process column, no data
*  needs to be exchanged by the processes, the operation is purely local.
*/
      if( !( PB_Cspan( N, JA, DESCA[INB_], AnbR, AsrcR, npcol ) ) )
      {
         if( Mupcase( ALLOC[0] ) == CREUSE )
         {
/*
*  sub( A ) can be reused
*/
            Bld = Ald;
            if( ( ( mycol == AsrcR ) || ( AsrcR < 0 ) ) && ( AnpD > 0 ) )
/*
*  If I own some entries of sub( A ), set *B
*/
               *B = Mptr( A, AiiD, AiiR, Ald, TYPE->size );
         }
         else
         {
/*
*  sub( A ) cannot be reused, make a copy of it.
*/
            Bld = MAX( 1, AnpD );
            if( ( ( mycol == AsrcR ) || ( AsrcR < 0 ) ) && ( AnpD > 0 ) )
            {
/*
*  If I own some entries of sub( A ), allocate space for the copy, and copy the
*  data.
*/
               if( AnpD > 0 )
               {
                  size   = TYPE->size;
                  *B     = PB_Cmalloc( AnpD * N * size );
                  *BFREE = 1;
                  TYPE->Fmmadd( &AnpD, &N, TYPE->one, Mptr( A, AiiD, AiiR, Ald,
                                size ), &Ald, TYPE->zero, *B, &Bld );
               }
            }
         }
/*
*  Describe the resulting operand
*/
         PB_Cdescset( DESCB, M, N, Ainb1D, N, AnbD, AnbR, AsrcD, AsrcR, ctxt,
                      Bld );
         return;
      }

      AnR      = N;     Bsrc_    = CSRC_;
      AmyprocR = mycol; AmyprocD = myrow; AnprocsR = npcol;
      Ainb1R   = PB_Cfirstnb( N, JA, DESCA[INB_], AnbR );
      AnpR     = PB_Cnumroc( N, 0, Ainb1R, AnbR, mycol, AsrcR, npcol );
   }
/*
*  sub( A ) is not replicated and spans more than one process row or column.
*  Forward row (resp. column) accumulation will leave the resulting operand in
*  the process(es) where the global row IA+M-1 (resp. global column JA+N-1)
*  resides.
*/
   if( ( Afwd = ( Mupcase( DIRECA[0] ) == CFORWARD ) ) != 0 )
   {
      if( ( AnpD > 0 ) && ( AnpR > 0 ) )
      {
/*
*  Compute how may rows or columns are before me -> AnpreR
*/
         AnpreR = PB_Cnpreroc( AnR, 0, Ainb1R, AnbR, AmyprocR, AsrcR,
                               AnprocsR );

         if( AnpreR == 0 )
         {
/*
*  If zero rows or columns are before me, I must be the source, so send my piece
*  to the process after me in the grid.
*/
            if( AggRow )
            {
               TYPE->Cgesd2d( ctxt, AnpR, AnpD, Mptr( A, AiiR, AiiD, Ald,
                              TYPE->size ), Ald, MModAdd1( AmyprocR, AnprocsR ),
                              AmyprocD );
            }
            else
            {
               TYPE->Cgesd2d( ctxt, AnpD, AnpR, Mptr( A, AiiD, AiiR,
                              Ald, TYPE->size ), Ald, AmyprocD,
                              MModAdd1( AmyprocR, AnprocsR ) );
            }
         }
         else if( AnpreR > 0 )
         {
/*
*  Otherwise, allocate some space for the rows or columns I have and the ones
*  globally preceeding the ones I have, that I am about to receive.
*/
            size     = TYPE->size; one = TYPE->one; zero = TYPE->zero;
            add      = TYPE->Fmmadd;
            *B       = Bptr = PB_Cmalloc( ( AnpreR + AnpR ) * AnpD * size );
            nlen     = AnpreR;
            mydistnb = MModSub( AmyprocR, AsrcR, AnprocsR ) * AnbR;
            kblks    = ( ( ( ktmp = AnR - Ainb1R - 1 ) >= 0 ) ?
                          ( ( ktmp / AnbR ) + 1 ) / AnprocsR : 0 );
            offset   = kblks * AnbR;
            kn       = Ainb1R + mydistnb - AnbR;
            kn       = MIN( kn, AnpreR ) +
                       ( MAX( 1, kblks ) - 1 ) * mydistnb;
            if( AggRow )
            {
               shft = TYPE->Frshft;
               Aptr = Mptr( A, AiiR, AiiD, Ald, size );
               Bld  = AnpreR + AnpR;
/*
*  Receive the rows globally preceeding the ones I have
*/
               TYPE->Cgerv2d( ctxt, AnpreR, AnpD, *B, Bld, MModSub1( AmyprocR,
                              AnprocsR ), AmyprocD );
/*
*  Sort the received buffer and insert at the correct place the rows of sub( A )
*  I own (from bottom to top).
*/
               if( ( ( AnpR - 1 ) / AnbR ) == kblks )
               {
                  kb      = AnpR - offset;
                  add( &kb, &AnpD, one, Mptr( Aptr, offset, 0, Ald, size ),
                       &Ald, zero, Mptr( Bptr, nlen+offset, 0, Bld, size ),
                       &Bld );
               }

               for( k = kblks; k >= 1; k-- )
               {
                  kb      = nlen - kn;
                  shft( &kb, &AnpD, &offset, Mptr( Bptr, kn, 0, Bld, size ),
                        &Bld );
                  offset -= AnbR;
                  add( &AnbR, &AnpD, one, Mptr( Aptr, offset, 0, Ald, size ),
                       &Ald, zero, Mptr( Bptr, kn+offset, 0, Bld, size ),
                       &Bld );
                  kn     -= mydistnb;
                  nlen   -= kb;
               }

               if( AnpreR + AnpR != AnR )
               {
/*
*  If I am not the last process, i.e I am not supposed to own all of the AnR
*  rows by the end of the operation, then send the sorted buffer to the next
*  process and release the dynamically allocated buffer.
*/
                  TYPE->Cgesd2d( ctxt, AnpreR+AnpR, AnpD, *B, Bld,
                                 MModAdd1( AmyprocR, AnprocsR ), AmyprocD );
                  if( *B ) free( *B );
               }
            }
            else
            {
               shft = TYPE->Fcshft;
               Aptr = Mptr( A, AiiD, AiiR, Ald, size );
               Bld  = MAX( 1, AnpD );
/*
*  Receive the columns globally preceeding the ones I have
*/
               TYPE->Cgerv2d( ctxt, AnpD, AnpreR, *B, Bld, AmyprocD,
                              MModSub1( AmyprocR, AnprocsR ) );
/*
*  Sort the received buffer and insert at the correct place the columns of
*  sub( A ) I own (from right to left).
*/
               if( ( ( AnpR - 1 ) / AnbR ) == kblks )
               {
                  kb      = AnpR - offset;
                  add( &AnpD, &kb, one, Mptr( Aptr, 0, offset, Ald, size ),
                       &Ald, zero, Mptr( Bptr, 0, nlen+offset, Bld, size ),
                       &Bld );
               }

               for( k = kblks; k >= 1; k-- )
               {
                  kb      = nlen - kn;
                  shft( &AnpD, &kb, &offset, Mptr( Bptr, 0, kn, Bld, size ),
                        &Bld );
                  offset -= AnbR;
                  add( &AnpD, &AnbR, one, Mptr( Aptr, 0, offset, Ald, size ),
                       &Ald, zero, Mptr( Bptr, 0, kn + offset, Bld, size ),
                       &Bld );
                  kn     -= mydistnb;
                  nlen   -= kb;
               }

               if( AnpreR + AnpR != AnR )
               {
/*
*  If I am not the last process, i.e I am not supposed to own all of the AnR
*  columns by the end of the operation, then send the sorted buffer to the next
*  process and release the dynamically allocated buffer.
*/
                  TYPE->Cgesd2d( ctxt, AnpD, AnpreR+AnpR, *B, Bld, AmyprocD,
                                 MModAdd1( AmyprocR, AnprocsR ) );
                  if( *B ) free( *B );
               }
            }
         }
      }
   }
   else
   {
/*
*  Backward accumulation, compute the process row or column coordinate ArocR,
*  that is going to have the resulting operand.
*/
      ArocR = PB_Cindxg2p( AnR-1, Ainb1R, AnbR, AsrcR, AsrcR, AnprocsR );

      if( ( AnpD > 0 ) && ( AnpR > 0 ) )
      {
/*
*  Compute how may rows or columns are after me -> AnnxtR
*/
         AnnxtR = PB_Cnnxtroc( AnR, 0, Ainb1R, AnbR, AmyprocR, AsrcR,
                               AnprocsR );
         AnnxtL = PB_Cnnxtroc( AnR, 0, Ainb1R, AnbR, ArocR,    AsrcR,
                               AnprocsR );

         if( ( AnnxtR = MModSub( AnnxtR, AnnxtL, AnR ) ) == 0 )
         {
/*
*  If zero rows or columns are after me, I must be the source, so send my piece
*  to the process before me in the grid.
*/
            if( AggRow )
            {
               TYPE->Cgesd2d( ctxt, AnpR, AnpD, Mptr( A, AiiR, AiiD, Ald,
                              TYPE->size ), Ald, MModSub1( AmyprocR, AnprocsR ),
                              AmyprocD );
            }
            else
            {
               TYPE->Cgesd2d( ctxt, AnpD, AnpR, Mptr( A, AiiD, AiiR, Ald,
                              TYPE->size ), Ald, AmyprocD, MModSub1( AmyprocR,
                              AnprocsR ) );
            }
         }
         else if( AnnxtR > 0 )
         {
/*
*  Otherwise, allocate some space for the rows or columns I have and the ones
*  globally following the ones I have, that I am about to receive.
*/
            size     = TYPE->size; one = TYPE->one; zero = TYPE->zero;
            add      = TYPE->Fmmadd;
            *B       = Bptr = PB_Cmalloc( ( AnnxtR + AnpR ) * AnpD * size );
            kblks    = ( ( ( ktmp = AnR - Ainb1R - 1 ) >= 0 ) ?
                         ( ( ktmp / AnbR ) + 1 ) / AnprocsR : 0 );
            mydist   = MModSub( ArocR, AmyprocR, AnprocsR );
            mydistnb = mydist * AnbR;
            srcdist  = MModSub( ArocR, AsrcR,    AnprocsR );

            if( AggRow )
            {
               shft = TYPE->Frshft;
               Aptr = Mptr( A, AiiR, AiiD, Ald, size );
               Bld  = AnnxtR + AnpR;
/*
*  Receive the rows globally following the ones I have
*/
               TYPE->Cgerv2d( ctxt, AnnxtR, AnpD, Mptr( *B, AnpR, 0, Bld,
                              size ), Bld, MModAdd1( AmyprocR, AnprocsR ),
                              AmyprocD );
/*
*  Sort the received buffer and insert at the correct place the rows of sub( A )
*  I own (from top to bottom).
*/
               if( mydist > srcdist )
               {
                  offset  = -AnpR;
                  kb      = Ainb1R + srcdist*AnbR;
               }
               else if( mydist == srcdist )
               {
                  add( &Ainb1R, &AnpD, one, Aptr, &Ald, zero, Bptr, &Bld );
                  Aptr    = Mptr( Aptr, Ainb1R, 0, Ald, size );
                  Bptr    = Mptr( Bptr, Ainb1R, 0, Ald, size );
                  offset  = Ainb1R - AnpR;
                  kb      = mydistnb;
               }
               else
               {
                  add( &AnbR, &AnpD, one, Aptr, &Ald, zero, Bptr, &Bld );
                  Aptr    = Mptr( Aptr, AnbR, 0, Ald, size );
                  Bptr    = Mptr( Bptr, AnbR, 0, Ald, size );
                  offset  = AnbR - AnpR;
                  kb      = mydistnb;
               }

               for( k = kblks; k >= 1; k-- )
               {
                  shft( &kb, &AnpD, &offset, Bptr, &Bld );
                  Bptr    = Mptr( Bptr, kb, 0, Bld, size );
                  add( &AnbR, &AnpD, one, Aptr, &Ald, zero, Bptr, &Bld );
                  Aptr    = Mptr( Aptr, AnbR, 0, Ald, size );
                  Bptr    = Mptr( Bptr, AnbR, 0, Ald, size );
                  offset += AnbR;
                  kb      = mydistnb;
               }

               if( AnnxtR + AnpR != AnR )
               {
/*
*  If I am not the last process, i.e I am not supposed to own all of the AnR
*  rows by the end of the operation, then send the sorted buffer to the previous
*  process and release the dynamically allocated buffer.
*/
                  TYPE->Cgesd2d( ctxt, AnnxtR+AnpR, AnpD, *B, Bld,
                                 MModSub1( AmyprocR, AnprocsR ), AmyprocD );
                  if( *B ) free( *B );
               }
            }
            else
            {
               shft = TYPE->Fcshft;
               Aptr = Mptr( A, AiiD, AiiR, Ald, size );
               Bld  = MAX( 1, AnpD );
/*
*  Receive the columns globally following the ones I have
*/
               TYPE->Cgerv2d( ctxt, AnpD, AnnxtR, Mptr( *B, 0, AnpR, Bld,
                              size ), Bld, AmyprocD, MModAdd1( AmyprocR,
                              AnprocsR ) );
/*
*  Sort the received buffer and insert at the correct place the columns of
*  sub( A ) I own (from left to right).
*/
               if( mydist > srcdist )
               {
                  offset  = -AnpR;
                  kb      = Ainb1R + srcdist*AnbR;
               }
               else if( mydist == srcdist )
               {
                  add( &AnpD, &Ainb1R, one, Aptr, &Ald, zero, Bptr, &Bld );
                  Aptr    = Mptr( Aptr, 0, Ainb1R, Ald, size );
                  Bptr    = Mptr( Bptr, 0, Ainb1R, Bld, size );
                  offset  = Ainb1R - AnpR;
                  kb      = mydistnb;
               }
               else
               {
                  add( &AnpD, &AnbR, one, Aptr, &Ald, zero, Bptr, &Bld );
                  Aptr    = Mptr( Aptr, 0, AnbR, Ald, size );
                  Bptr    = Mptr( Bptr, 0, AnbR, Bld, size );
                  offset  = AnbR - AnpR;
                  kb      = mydistnb;
               }

               for( k = kblks; k >= 1; k-- )
               {
                  shft( &AnpD, &kb, &offset, Bptr, &Bld );
                  Bptr    = Mptr( Bptr, 0, kb, Bld, size );
                  add( &AnpD, &AnbR, one, Aptr, &Ald, zero, Bptr, &Bld );
                  Aptr    = Mptr( Aptr, 0, AnbR, Ald, size );
                  Bptr    = Mptr( Bptr, 0, AnbR, Bld, size );
                  offset += AnbR;
                  kb      = mydistnb;
               }

               if( AnnxtR + AnpR != AnR )
               {
/*
*  If I am not the last process, i.e I am not supposed to own all of the AnR
*  columns by the end of the operation, then send the sorted buffer to the
*  previous process and release the dynamically allocated buffer.
*/
                  TYPE->Cgesd2d( ctxt, AnpD, AnnxtR+AnpR, *B, Bld, AmyprocD,
                                 MModSub1( AmyprocR, AnprocsR ) );
                  if( *B ) free( *B );
               }
            }
         }
      }
   }
/*
*  Describe the resulting operand
*/
   if( AggRow )
   {
      PB_Cdescset( DESCB, M, N, M, Ainb1D, AnbR, AnbD, AsrcR, AsrcD, ctxt, M );
   }
   else
   {
      PB_Cdescset( DESCB, M, N, Ainb1D, N, AnbD, AnbR, AsrcD, AsrcR, ctxt,
                   MAX( 1, AnpD ) );
   }
/*
*  Compute globally in which process row or column the resulting operand is
*  residing and set *BFREE accordingly.
*/
   if( Afwd )
   {
      if( AnR + AnbR > Ainb1R + ( AnprocsR - 1 ) * AnbR )
      {
/*
*  If sub( A ) is spanning all process rows or columns of the grid, the result
*  must be in the process row or column preceeding the one owning IA or JA,
*  don't you think ?
*/
         DESCB[Bsrc_] = MModSub1( AsrcR, AnprocsR );
      }
      else
      {
/*
*  Otherwise, the result is in the process row or column where the row IA+M-1
*  or column JA+N-1 of sub( A ) resides.
*/
         DESCB[Bsrc_] = PB_Cindxg2p( AnR-1, Ainb1R, AnbR, AsrcR, AsrcR,
                                     AnprocsR );
      }
      if( ( AnpD > 0 ) && ( AnpR > 0 ) && ( AmyprocR == DESCB[Bsrc_] ) )
         *BFREE = 1;
   }
   else
   {
      if( AnR + AnbR > Ainb1R + ( AnprocsR - 1 ) * AnbR )
      {
/*
*  If sub( A ) is spanning all process rows or columns of the grid, the result
*  must be in the process row or column following the one owning IA+M-1 or
*  JA+N-1, don't you think ?
*/
         DESCB[Bsrc_] = MModAdd1( ArocR, AnprocsR );
      }
      else
      {
/*
*  Otherwise, the result is in the process row or column where the row IA or
*  column JA of sub( A ) resides.
*/
         DESCB[Bsrc_] = AsrcR;
      }
      if( ( AnpD > 0 ) && ( AnpR > 0 ) && ( AmyprocR == DESCB[Bsrc_] ) )
         *BFREE = 1;
   }
/*
*  End of PB_CGatherV
*/
}
