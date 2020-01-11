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
void PB_CpaxpbyNN( PBTYP_T * TYPE, char * CONJUG, Int M, Int N,
                   char * ALPHA,
                   char * A, Int IA, Int JA, Int * DESCA, char * AROC,
                   char * BETA,
                   char * B, Int IB, Int JB, Int * DESCB, char * BROC )
#else
void PB_CpaxpbyNN( TYPE, CONJUG, M, N, ALPHA, A, IA, JA, DESCA, AROC,
                   BETA, B, IB, JB, DESCB, BROC )
/*
*  .. Scalar Arguments ..
*/
   char           * AROC, * BROC, * CONJUG;
   Int            IA, IB, JA, JB, M, N;
   char           * ALPHA, * BETA;
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
*  PB_CpaxpbyNN adds one submatrix to another,
*
*     sub( B ) := beta * sub( B ) + alpha * sub( A ), or,
*
*     sub( B ) := beta * sub( B ) + alpha * conjg( sub( A ) ),
*
*  where both submatrices are not distributed;  sub( A ) always  denotes
*  A(IA:IA+M-1,JA:JA+N-1).  When AROC is 'R' or 'r'  sub( A ) resides in
*  a  process row, otherwise sub( A ) resides in a process column.  When
*  sub( A )  resides  in  a  process  row  and  BROC  is  'R' or  'r' or
*  sub( A ) resides  in  a  process  column and BROC is 'C' or 'c', then
*  sub( B )  denotes B(IB:IB+M-1,JB:JB+N-1), and  B(IB:IB+N-1,JB:JB+M-1)
*  otherwise.
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
*  CONJUG  (global input) pointer to CHAR
*          On  entry,  CONJUG  specifies  whether  conjg( sub( A ) )  or
*          sub( A ) should be added to sub( B ) as follows:
*             CONJUG = 'N' or 'n':
*                sub( B ) := beta*sub( B ) + alpha*sub( A ),
*             otherwise
*                sub( B ) := beta*sub( B ) + alpha*conjg( sub( A ) ).
*
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( A ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( A ). N must be at least zero.
*
*  ALPHA   (global input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied as zero then the local entries of the array  A  cor-
*          responding to the entries of the submatrix sub( A ) need  not
*          be set on input.
*
*  A       (local input) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where LLD_A
*          is at least MAX( 1, Lr( 1, IA+M-1 ) ),  and,  Ka is  at least
*          Lc( 1, JA+N-1 ). Before entry, this array contains the  local
*          entries of the matrix A.
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
*          On entry,  AROC  specifies  the  orientation of the subvector
*          sub( A ). When AROC is 'R' or 'r',  sub( A ) is a row vector,
*          and a column vector otherwise.
*
*  BETA    (global input) pointer to CHAR
*          On entry, BETA specifies the scalar beta. When  BETA  is sup-
*          plied as zero then the local entries of the array  B  corres-
*          ponding to the entries of the submatrix sub( B ) need not  be
*          set on input.
*
*  B       (local input/local output) pointer to CHAR
*          On entry, B is an array of dimension (LLD_B, Kb), where LLD_B
*          is  at  least  MAX( 1, Lr( 1, IB+M-1 ) )  when  sub( A )  and
*          sub( B )  are  both  distributed  along a process column or a
*          process row. In that case, Kb is  at  least  Lc( 1, JB+N-1 ).
*          Otherwise,  LLD_B  is at least  MAX( 1, Lr( 1, IB+N-1 ) ) and
*          Kb  is at least Lc( 1, JB+M-1 ).  Before  entry,  this  array
*          contains the local entries of the matrix B. On exit, sub( B )
*          is overwritten with the updated submatrix.
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
*          On entry,  BROC  specifies  the  orientation of the subvector
*          sub( B ). When BROC is 'R' or 'r',  sub( B ) is a row vector,
*          and a column vector otherwise.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           scope, * top;
   Int            Acol, Aii, AisR, AisRow, Ajj, Ald, AmyprocD, AmyprocR,
                  AnprocsD, AnprocsR, AprocR, Arow, Bcol, Bii, BisR, BisRow,
                  Bjj, Bld, BmyprocD, BmyprocR, BnprocsD, BnprocsR, BprocR,
                  Brow, RRorCC, csrc, ctxt, iroca, mycol, myrow, npcol, nprow,
                  p, rsrc, size;
   MMADD_T        add;
/*
*  .. Local Arrays ..
*/
   char           * buf = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Retrieve sub( A )'s local information: Aii, Ajj, Arow, Acol ...
*/
   PB_Cinfog2l( IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aii, &Ajj,
                &Arow, &Acol );
   if( ( AisRow = ( Mupcase( AROC[0] ) == CROW ) ) != 0 )
   {
      Ald    = DESCA[LLD_]; AmyprocD = mycol; AnprocsD = npcol;
      AprocR = Arow;        AmyprocR = myrow; AnprocsR = nprow;
      AisR   = ( ( Arow == -1 ) || ( AnprocsR == 1 ) );
   }
   else
   {
      Ald    = DESCA[LLD_]; AmyprocD = myrow; AnprocsD = nprow;
      AprocR = Acol;        AmyprocR = mycol; AnprocsR = npcol;
      AisR   = ( ( Acol == -1 ) || ( AnprocsR == 1 ) );
   }
/*
*  Retrieve sub( B )'s local information: Bii, Bjj, Brow, Bcol ...
*/
   PB_Cinfog2l( IB, JB, DESCB, nprow, npcol, myrow, mycol, &Bii, &Bjj,
                &Brow, &Bcol );
   if( ( BisRow = ( Mupcase( BROC[0] ) == CROW ) ) != 0 )
   {
      Bld    = DESCB[LLD_]; BmyprocD = mycol; BnprocsD = npcol;
      BprocR = Brow;        BmyprocR = myrow; BnprocsR = nprow;
      BisR   = ( ( Brow == -1 ) || ( BnprocsR == 1 ) );
   }
   else
   {
      Bld    = DESCB[LLD_]; BmyprocD = myrow; BnprocsD = nprow;
      BprocR = Bcol;        BmyprocR = mycol; BnprocsR = npcol;
      BisR   = ( ( Bcol == -1 ) || ( BnprocsR == 1 ) );
   }
/*
*  Are sub( A ) and sub( B ) both row or column vectors ?
*/
   RRorCC = ( ( AisRow && BisRow ) || ( !( AisRow ) && !( BisRow ) ) );
/*
*  Neither sub( A ) nor sub( B ) are distributed
*/
   if( !AisR )
   {
/*
*  sub( A ) is not replicated
*/
      if( !( BisR ) )
      {
/*
*  sub( B ) is not replicated
*/
         if( ( AmyprocR != AprocR ) && ( BmyprocR != BprocR ) )
/*
*  If I am not in AprocR or BprocR, then return immediately
*/
            return;

         size = TYPE->size;

         if( RRorCC )
         {
/*
*  sub( A ) and sub( B ) are both row or column vectors
*/
            if( Mupcase( CONJUG[0] ) != CNOCONJG ) add = TYPE->Fmmcadd;
            else                                   add = TYPE->Fmmadd;

            if( AprocR == BprocR )
            {
               add( &M, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ), &Ald, BETA,
                    Mptr( B, Bii, Bjj, Bld, size ), &Bld );
            }
            else
            {
/*
*  sub( A ) and sub( B ) are in a different process row or column
*/
               if( AmyprocR == AprocR )
               {
/*
*  Send sub( A ) to where sub( B ) resides.
*/
                  if( AisRow )
                     TYPE->Cgesd2d( ctxt, M, N, Mptr( A, Aii, Ajj, Ald, size ),
                                    Ald, BprocR, AmyprocD );
                  else
                     TYPE->Cgesd2d( ctxt, M, N, Mptr( A, Aii, Ajj, Ald, size ),
                                    Ald, AmyprocD, BprocR );
               }
/*
*  receive sub( A ) and add it to sub( B )
*/
               if( BmyprocR == BprocR )
               {
                  buf = PB_Cmalloc( M * N * size );
                  if( BisRow )
                     TYPE->Cgerv2d( ctxt, M, N, buf, M, AprocR, BmyprocD );
                  else
                     TYPE->Cgerv2d( ctxt, M, N, buf, M, BmyprocD, AprocR );
                  add( &M, &N, ALPHA, buf, &M, BETA, Mptr( B, Bii, Bjj, Bld,
                       size ), &Bld );
                  if( buf ) free( buf );
               }
            }
         }
         else
         {
/*
*  sub( A ) and sub( B ) are not both row or column vectors
*/
            if( Mupcase( CONJUG[0] ) != CNOCONJG ) add = TYPE->Fmmtcadd;
            else                                   add = TYPE->Fmmtadd;

            iroca = 0;
            for( p = 0; p < BnprocsD; p++ )
            {
               if( ( AprocR == p ) && ( BprocR == iroca ) )
               {
                  if( ( AmyprocR == p ) && ( AmyprocD == iroca ) )
                  {
                     add( &M, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ), &Ald,
                          BETA, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
                  }
               }
               else
               {
                  if( ( AmyprocR == AprocR ) && ( AmyprocD == iroca ) )
                  {
                     if( AisRow )
                        TYPE->Cgesd2d( ctxt, M, N, Mptr( A, Aii, Ajj, Ald,
                                       size ), Ald, p, BprocR );
                     else
                        TYPE->Cgesd2d( ctxt, M, N, Mptr( A, Aii, Ajj, Ald,
                                       size ), Ald, BprocR, p );
                  }
                  if( ( BmyprocR == BprocR ) && ( BmyprocD == p ) )
                  {
                     buf = PB_Cmalloc( M * N * size );
                     if( AisRow )
                        TYPE->Cgerv2d( ctxt, M, N, buf, M, AprocR, iroca );
                     else
                        TYPE->Cgerv2d( ctxt, M, N, buf, M, iroca, AprocR );
                     add( &M, &N, ALPHA, buf, &M, BETA, Mptr( B, Bii, Bjj, Bld,
                          size ), &Bld );
                     if( buf ) free( buf );
                  }
               }
               iroca = MModAdd1( iroca, AnprocsD );
            }
         }
      }
      else
      {
/*
*  sub( B ) is replicated
*/
         size = TYPE->size;

         if( AmyprocR == AprocR )
         {
            if( RRorCC )
            {
               if( Mupcase( CONJUG[0] ) != CNOCONJG )
                  TYPE->Fmmcadd( &M, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ),
                                 &Ald, BETA, Mptr( B, Bii, Bjj, Bld, size ),
                                 &Bld );
               else
                  TYPE->Fmmadd( &M, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ),
                                &Ald, BETA, Mptr( B, Bii, Bjj, Bld, size ),
                                &Bld );
            }
            else
            {
               if( Mupcase( CONJUG[0] ) != CNOCONJG )
                  TYPE->Fmmtcadd( &M, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ),
                                  &Ald, BETA, Mptr( B, Bii, Bjj, Bld, size ),
                                  &Bld );
               else
                  TYPE->Fmmtadd( &M, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ),
                                 &Ald, BETA, Mptr( B, Bii, Bjj, Bld, size ),
                                 &Bld );
            }
            if( AisRow ) { scope = CCOLUMN; } else { scope = CROW; }
            top = PB_Ctop( &ctxt, BCAST, &scope, TOP_GET );
            if( RRorCC )
               TYPE->Cgebs2d( ctxt, &scope, top, M, N, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld );
            else
               TYPE->Cgebs2d( ctxt, &scope, top, N, M, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld );
         }
         else
         {
            if( AisRow ) { scope = CCOLUMN; rsrc = AprocR; csrc = AmyprocD; }
            else         { scope = CROW;    rsrc = AmyprocD; csrc = AprocR; }
            top = PB_Ctop( &ctxt, BCAST, &scope, TOP_GET );
            if( RRorCC )
               TYPE->Cgebr2d( ctxt, &scope, top, M, N, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld, rsrc, csrc );
            else
               TYPE->Cgebr2d( ctxt, &scope, top, N, M, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld, rsrc, csrc );
         }
      }
   }
   else
   {
/*
*  sub( A ) is replicated
*/
      if( BisR || ( BmyprocR == BprocR ) )
      {
/*
*  If I own a piece of sub( B ), then add sub( A ) to it
*/
         size = TYPE->size;
         if( RRorCC )
         {
            if( Mupcase( CONJUG[0] ) != CNOCONJG )
               TYPE->Fmmcadd( &M, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ),
                              &Ald, BETA, Mptr( B, Bii, Bjj, Bld, size ),
                              &Bld );
            else
               TYPE->Fmmadd( &M, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ),
                             &Ald, BETA, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
         }
         else
         {
            if( Mupcase( CONJUG[0] ) != CNOCONJG )
               TYPE->Fmmtcadd( &M, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ),
                               &Ald, BETA, Mptr( B, Bii, Bjj, Bld, size ),
                               &Bld );
            else
               TYPE->Fmmtadd( &M, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ),
                              &Ald, BETA, Mptr( B, Bii, Bjj, Bld, size ),
                              &Bld );
         }
      }
   }
/*
*  End of PB_CpaxpbyNN
*/
}
