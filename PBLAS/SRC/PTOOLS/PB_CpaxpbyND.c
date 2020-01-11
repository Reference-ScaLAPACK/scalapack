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
void PB_CpaxpbyND( PBTYP_T * TYPE, char * CONJUG, Int M, Int N,
                   char * ALPHA,
                   char * A, Int IA, Int JA, Int * DESCA, char * AROC,
                   char * BETA,
                   char * B, Int IB, Int JB, Int * DESCB, char * BROC )
#else
void PB_CpaxpbyND( TYPE, CONJUG, M, N, ALPHA, A, IA, JA, DESCA, AROC,
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
*  PB_CpaxpbyND adds one submatrix to another,
*
*     sub( B ) := beta * sub( B ) + alpha * sub( A ), or,
*
*     sub( B ) := beta * sub( B ) + alpha * conjg( sub( A ) ),
*
*  where sub( A ) is not distributed and sub( B ) is distributed.
*
*  sub( A ) always  denotes A(IA:IA+M-1,JA:JA+N-1).  When AROC is 'R' or
*  'r'  sub( A ) resides in a process row, otherwise sub( A ) resides in
*  a process column. When sub( A ) resides in a  process row and BROC is
*  'R' or 'r' or sub( A ) resides in a process column and BROC is 'C' or
*  'c', then sub( B ) denotes  B( IB:IB+M-1, JB:JB+N-1 ), and  otherwise
*  sub( B ) denotes B(IB:IB+N-1,JB:JB+M-1).
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
   char           * one, * top, * zero;
   Int            Acol, Aii, AisR, AisRow, Ajj, Ald, AmyprocD, AmyprocR,
                  AnprocsD, AprocR, Aroc, Arow, Bcol, Bii, Binb1D, BisR, BisRow,
                  Bjj, Bld, BmyprocD, BmyprocR, BnD, BnbD, BnpD, BnprocsD,
                  BprocD,  BprocR, Broc, Brow, RRorCC, ctxt, k, kbb, kk, kn,
                  ktmp, mycol, mydist, myproc, myrow, npcol, nprow, p, size;
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
      BnD      = N;     Ald      = DESCA[LLD_];
      AmyprocD = mycol; AnprocsD = npcol; AmyprocR = myrow; AprocR = Arow;
      AisR     = ( ( Arow == -1 ) || ( nprow == 1 ) );
   }
   else
   {
      BnD      = M;     Ald      = DESCA[LLD_];
      AmyprocD = myrow; AnprocsD = nprow; AmyprocR = mycol; AprocR = Acol;
      AisR     = ( ( Acol == -1 ) || ( npcol == 1 ) );
   }
/*
*  Retrieve sub( B )'s local information: Bii, Bjj, Brow, Bcol ...
*/
   PB_Cinfog2l( IB, JB, DESCB, nprow, npcol, myrow, mycol, &Bii, &Bjj,
                &Brow, &Bcol );
   if( ( BisRow = ( Mupcase( BROC[0] ) == CROW ) ) != 0 )
   {
      BnbD     = DESCB[NB_]; Bld      = DESCB[LLD_];
      BprocD   = Bcol;       BmyprocD = mycol;       BnprocsD = npcol;
      BprocR   = Brow;       BmyprocR = myrow;
      BisR     = ( ( BprocR == -1 ) || ( nprow == 1 ) );
      Binb1D   = PB_Cfirstnb( BnD, JB, DESCB[INB_], BnbD );
   }
   else
   {
      BnbD     = DESCB[MB_]; Bld      = DESCB[LLD_];
      BprocD   = Brow;       BmyprocD = myrow;       BnprocsD = nprow;
      BprocR   = Bcol;       BmyprocR = mycol;
      BisR     = ( ( BprocR == -1 ) || ( npcol == 1 ) );
      Binb1D   = PB_Cfirstnb( BnD, IB, DESCB[IMB_], BnbD );
   }
/*
*  Are sub( A ) and sub( B ) both row or column vectors ?
*/
   RRorCC = ( ( AisRow && BisRow ) || ( !( AisRow ) && !( BisRow ) ) );
/*
*  sub( A ) is not distributed and sub( B ) is distributed
*/
   if( !( AisR ) )
   {
/*
*  sub( A ) is not replicated. Since this operation is local if sub( A ) and
*  sub( B ) are both row or column vectors, choose BprocR = AprocR when RRorCC,
*  and BprocR = 0 otherwise.
*/
      if( BisR ) { BprocR = ( ( RRorCC ) ? AprocR : 0 ); }
/*
*  Now, it is just like sub( B ) is not replicated, this information however is
*  kept in BisR for later use.
*/
      size = TYPE->size;

      if( ( AmyprocR == AprocR ) || ( BmyprocR == BprocR ) )
      {
         zero = TYPE->zero; one = TYPE->one;
/*
*  sub( A ) and sub( B ) are both row or column vectors
*/
         if( RRorCC )
         {
            if( Mupcase( CONJUG[0] ) != CNOCONJG ) add = TYPE->Fmmcadd;
            else                                   add = TYPE->Fmmadd;

            BnpD = PB_Cnumroc( BnD, 0, Binb1D, BnbD, BmyprocD, BprocD,
                               BnprocsD );
/*
*  sub( A ) and sub( B ) are in the same process row or column
*/
            if( AprocR == BprocR )
            {
/*
*  In each process, the non distributed part of sub( A ) is added to sub( B ).
*/
               if( BnpD > 0 )
               {
                  Broc = BprocD;
                  if( AisRow ) { kk = Bjj; ktmp = JA + N; kn = JA + Binb1D; }
                  else         { kk = Bii; ktmp = IA + M; kn = IA + Binb1D; }

                  if( BmyprocD == Broc )
                  {
                     if( AisRow )
                        add( &M, &Binb1D, ALPHA, Mptr( A, Aii, Ajj, Ald, size ),
                             &Ald, BETA, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
                     else
                        add( &Binb1D, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ),
                             &Ald, BETA, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
                     kk += Binb1D;
                  }
                  Broc = MModAdd1( Broc, BnprocsD );

                  for( k = kn; k < ktmp; k += BnbD )
                  {
                     kbb = ktmp - k; kbb = MIN( kbb, BnbD );

                     if( BmyprocD == Broc )
                     {
                        if( AisRow )
                           add( &M, &kbb, ALPHA, Mptr( A, Aii, k, Ald, size ),
                                &Ald, BETA, Mptr( B, Bii, kk, Bld, size ),
                                &Bld );
                        else
                           add( &kbb, &N, ALPHA, Mptr( A, k, Ajj, Ald, size ),
                                &Ald, BETA, Mptr( B, kk, Bjj, Bld, size ),
                                &Bld );
                        kk += kbb;
                     }
                     Broc = MModAdd1( Broc, BnprocsD );
                  }
               }
            }
            else
            {
/*
*  sub( A ) and sub( B ) are in a different process row or column
*/
               if( BmyprocR == BprocR )
               {
/*
*  If I own a piece of sub( B ), then receive the relevant piece of sub( A )
*  from the corresponding process row or column where it resides.
*/
                  if( BnpD > 0 )
                  {
                     if( BisRow )
                     {
                        buf = PB_Cmalloc( M * BnpD * size );
                        TYPE->Cgerv2d( ctxt, M, BnpD, buf, M, AprocR,
                                       BmyprocD );
                        add( &M, &BnpD, ALPHA, buf, &M, BETA, Mptr( B, Bii, Bjj,
                             Bld, size ), &Bld );
                        if( buf ) free( buf );
                     }
                     else
                     {
                        buf = PB_Cmalloc( BnpD * N * size );
                        TYPE->Cgerv2d( ctxt, BnpD, N, buf, BnpD, BmyprocD,
                                       AprocR );
                        add( &BnpD, &N, ALPHA, buf, &BnpD, BETA, Mptr( B, Bii,
                             Bjj, Bld, size ), &Bld );
                        if( buf ) free( buf );
                     }
                  }
               }

               if( AmyprocR == AprocR )
               {
/*
*  If I own sub( A ), then pack and send the distributed part that should be
*  added to the distributed part of sub( B ) that resides in my row or column.
*/
                  if( BnpD > 0 )
                  {
                     if( AisRow )
                     {
                        ktmp = JA + N;
                        kn   = JA + Binb1D;
                        buf = PB_Cmalloc( M * BnpD * size );
                     }
                     else
                     {
                        ktmp = IA + M;
                        kn   = IA + Binb1D;
                        buf = PB_Cmalloc( BnpD * N * size );
                     }
                     Broc = BprocD;
                     kk   = 0;

                     if( BmyprocD == Broc )
                     {
                        if( AisRow )
                           add( &M, &Binb1D, one, Mptr( A, Aii, Ajj, Ald,
                                size ), &Ald, zero, buf, &M );
                        else
                           add( &Binb1D, &N, one, Mptr( A, Aii, Ajj, Ald,
                                size ), &Ald, zero, buf, &BnpD );
                        kk += Binb1D;
                     }

                     Broc = MModAdd1( Broc, BnprocsD );
                     for( k = kn; k < ktmp; k += BnbD )
                     {
                        kbb = ktmp - k; kbb = MIN( kbb, BnbD );

                        if( BmyprocD == Broc )
                        {
                           if( AisRow )
                              add( &M, &kbb, one, Mptr( A, Aii, k, Ald, size ),
                                   &Ald, zero, Mptr( buf, 0, kk, M, size ),
                                   &M );
                           else
                              add( &kbb, &N, one, Mptr( A, k, Ajj, Ald, size ),
                                   &Ald, zero, Mptr( buf, kk, 0, BnpD, size ),
                                   &BnpD );
                           kk += kbb;
                        }
                        Broc = MModAdd1( Broc, BnprocsD );
                     }

                     if( AisRow )
                        TYPE->Cgesd2d( ctxt, M, BnpD, buf,    M, BprocR,
                                       AmyprocD );
                     else
                        TYPE->Cgesd2d( ctxt, BnpD, N, buf, BnpD, AmyprocD,
                                       BprocR );
                     if( buf ) free( buf );
                  }
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

            Aroc = 0;
            if( AisRow ) { ktmp = JA + N; kn = JA + Binb1D; }
            else         { ktmp = IA + M; kn = IA + Binb1D; }
/*
*  Loop over the processes in which sub( B ) resides, for each process find the
*  next process Xroc. Exchange and add the data.
*/
            for( p = 0; p < BnprocsD; p++ )
            {
               mydist = MModSub(      p, BprocD, BnprocsD );
               myproc = MModAdd( BprocD, mydist, BnprocsD );

               if( ( AprocR == p ) && ( BprocR == Aroc ) )
               {
                  if( ( AmyprocR == p ) && ( AmyprocD == Aroc ) )
                  {
/*
*  local add at the intersection of the process cross
*/
                     BnpD = PB_Cnumroc( BnD, 0, Binb1D, BnbD, p, BprocD,
                                        BnprocsD );
                     if( BnpD > 0 )
                     {
                        Broc = BprocD;
                        kk   = ( AisRow ? Bii : Bjj );

                        if( myproc == Broc )
                        {
                           if( AisRow )
                              add( &M, &Binb1D, ALPHA, Mptr( A, Aii, Ajj, Ald,
                                   size ), &Ald, BETA, Mptr( B, Bii, Bjj, Bld,
                                   size ), &Bld );
                           else
                              add( &Binb1D, &N, ALPHA, Mptr( A, Aii, Ajj, Ald,
                                   size ), &Ald, BETA, Mptr( B, Bii, Bjj, Bld,
                                   size ), &Bld );
                           kk += Binb1D;
                        }
                        Broc = MModAdd1( Broc, BnprocsD );

                        for( k = kn; k < ktmp; k += BnbD )
                        {
                           kbb = ktmp - k; kbb = MIN( kbb, BnbD );
                           if( myproc == Broc )
                           {
                              if( AisRow )
                                 add( &M, &kbb, ALPHA, Mptr( A, Aii, k, Ald,
                                      size ), &Ald, BETA, Mptr( B, kk, Bjj, Bld,
                                      size ), &Bld );
                              else
                                 add( &kbb, &N, ALPHA, Mptr( A, k, Ajj, Ald,
                                      size ), &Ald, BETA, Mptr( B, Bii, kk, Bld,
                                      size ), &Bld );
                              kk += kbb;
                           }
                           Broc = MModAdd1( Broc, BnprocsD );
                        }
                     }
                  }
               }
               else
               {
/*
*  Message exchange
*/
                  if( ( BmyprocR == BprocR ) && ( BmyprocD == p ) )
                  {
                     BnpD = PB_Cnumroc( BnD, 0, Binb1D, BnbD, p, BprocD,
                                        BnprocsD );
                     if( BnpD > 0 )
                     {
                        if( AisRow )
                        {
                           buf = PB_Cmalloc( M * BnpD * size );
                           TYPE->Cgerv2d( ctxt, BnpD, M, buf, BnpD, AprocR,
                                          Aroc );
                           TYPE->Fmmadd( &BnpD, &M, ALPHA, buf, &BnpD, BETA,
                                         Mptr( B, Bii, Bjj, Bld, size ), &Bld );
                        }
                        else
                        {
                           buf = PB_Cmalloc( BnpD * N * size );
                           TYPE->Cgerv2d( ctxt, N, BnpD, buf, N, Aroc, AprocR );
                           TYPE->Fmmadd( &N, &BnpD, ALPHA, buf, &N, BETA,
                                         Mptr( B, Bii, Bjj, Bld, size ), &Bld );
                        }
                        if( buf ) free( buf );
                     }
                  }

                  if( ( AmyprocR == AprocR ) && ( AmyprocD == Aroc ) )
                  {
                     BnpD = PB_Cnumroc( BnD, 0, Binb1D, BnbD, p, BprocD,
                                        BnprocsD );
                     if( BnpD > 0 )
                     {
                        if( AisRow )
                           buf = PB_Cmalloc( M * BnpD * size );
                        else
                           buf = PB_Cmalloc( BnpD * N * size );
                        Broc = BprocD;
                        kk   = 0;

                        if( myproc == Broc )
                        {
                           if( AisRow )
                              add( &M, &Binb1D, one, Mptr( A, Aii, Ajj, Ald,
                                   size ), &Ald, zero, buf, &BnpD );
                           else
                              add( &Binb1D, &N, one, Mptr( A, Aii, Ajj, Ald,
                                   size ), &Ald, zero, buf, &N );
                           kk += Binb1D;
                        }
                        Broc = MModAdd1( Broc, BnprocsD );

                        for( k = kn; k < ktmp; k += BnbD )
                        {
                           kbb = ktmp - k; kbb = MIN( kbb, BnbD );
                           if( myproc == Broc )
                           {
                              if( AisRow )
                                 add( &M, &kbb, one, Mptr( A, Aii, k, Ald,
                                      size ), &Ald, zero, Mptr( buf, kk, 0,
                                      BnpD, size ), &BnpD );
                              else
                                 add( &kbb, &N, one, Mptr( A, k, Ajj, Ald,
                                      size ), &Ald, zero, Mptr( buf, 0, kk,
                                      N, size ), &N );
                              kk += kbb;
                           }
                           Broc = MModAdd1( Broc, BnprocsD );
                        }
                        if( AisRow )
                           TYPE->Cgesd2d( ctxt, BnpD, M, buf, BnpD, p, BprocR );
                        else
                           TYPE->Cgesd2d( ctxt, N, BnpD, buf,    N, BprocR, p );
                        if( buf ) free( buf );
                     }
                  }
               }
               Aroc = MModAdd1( Aroc, AnprocsD );
            }
         }
      }

      if( BisR )
      {
/*
*  Replicate sub( B )
*/
         BnpD = PB_Cnumroc( BnD, 0, Binb1D, BnbD, BmyprocD, BprocD, BnprocsD );
         if( BnpD > 0 )
         {
            if( BisRow )
            {
               top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
               if( BmyprocR == BprocR )
                  TYPE->Cgebs2d( ctxt, COLUMN, top, ( AisRow ? M : N ), BnpD,
                                 Mptr( B, Bii, Bjj, Bld, size ), Bld );
               else
                  TYPE->Cgebr2d( ctxt, COLUMN, top, ( AisRow ? M : N ), BnpD,
                                 Mptr( B, Bii, Bjj, Bld, size ), Bld, BprocR,
                                 BmyprocD );
            }
            else
            {
               top = PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
               if( BmyprocR == BprocR )
                  TYPE->Cgebs2d( ctxt, ROW,    top, BnpD, ( AisRow ? M : N ),
                                 Mptr( B, Bii, Bjj, Bld, size ), Bld );
               else
                  TYPE->Cgebr2d( ctxt, ROW,    top, BnpD, ( AisRow ? M : N ),
                                 Mptr( B, Bii, Bjj, Bld, size ), Bld, BmyprocD,
                                 BprocR );
            }
         }
      }
   }
   else
   {
/*
*  sub( A ) is replicated in every process. Add the data in process row or
*  column BprocR when sub( B ) is not replicated and in every process otherwise.
*/
      if( !( BisR ) && ( BmyprocR != BprocR ) ) return;

      size = TYPE->size;

      if( RRorCC )
      {
         if( Mupcase( CONJUG[0] ) != CNOCONJG ) add = TYPE->Fmmcadd;
         else                                   add = TYPE->Fmmadd;
      }
      else
      {
         if( Mupcase( CONJUG[0] ) != CNOCONJG ) add = TYPE->Fmmtcadd;
         else                                   add = TYPE->Fmmtadd;
      }

      Broc = BprocD;
      kk   = ( BisRow ? Bjj : Bii );
      if( AisRow ) { ktmp = JA + N; kn = JA + Binb1D; }
      else         { ktmp = IA + M; kn = IA + Binb1D; }

      if( BmyprocD == Broc )
      {
         if( AisRow )
            add( &M, &Binb1D, ALPHA, Mptr( A, Aii, Ajj, Ald, size ), &Ald, BETA,
                 Mptr( B, Bii, Bjj, Bld, size ), &Bld );
         else
            add( &Binb1D, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ), &Ald, BETA,
                 Mptr( B, Bii, Bjj, Bld, size ), &Bld );
         kk += Binb1D;
      }
      Broc = MModAdd1( Broc, BnprocsD );

      for( k = kn; k < ktmp; k += BnbD )
      {
         kbb = ktmp - k; kbb = MIN( kbb, BnbD );
         if( BmyprocD == Broc )
         {
            if( BisRow ) { buf = Mptr( B, Bii, kk, Bld, size ); }
            else         { buf = Mptr( B, kk, Bjj, Bld, size ); }
            if( AisRow )
               add( &M, &kbb, ALPHA, Mptr( A, Aii, k, Ald, size ), &Ald, BETA,
                    buf, &Bld );
            else
               add( &kbb, &N, ALPHA, Mptr( A, k, Ajj, Ald, size ), &Ald, BETA,
                    buf, &Bld );
            kk += kbb;
         }
         Broc = MModAdd1( Broc, BnprocsD );
      }
   }
/*
*  End of PB_CpaxpbyND
*/
}
