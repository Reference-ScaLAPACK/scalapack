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
void PB_CpaxpbyDN( PBTYP_T * TYPE, char * CONJUG, int M, int N,
                   char * ALPHA,
                   char * A, int IA, int JA, int * DESCA, char * AROC,
                   char * BETA,
                   char * B, int IB, int JB, int * DESCB, char * BROC )
#else
void PB_CpaxpbyDN( TYPE, CONJUG, M, N, ALPHA, A, IA, JA, DESCA, AROC,
                   BETA, B, IB, JB, DESCB, BROC )
/*
*  .. Scalar Arguments ..
*/
   char           * AROC, * BROC, * CONJUG;
   int            IA, IB, JA, JB, M, N;
   char           * ALPHA, * BETA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCB;
   char           * A, * B;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CpaxpbyDN adds one submatrix to another,
*
*     sub( B ) := beta * sub( B ) + alpha * sub( A ), or,
*
*     sub( B ) := beta * sub( B ) + alpha * conjg( sub( A ) ),
*
*  where sub( A ) is distributed and sub( B ) is not distributed.
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
   char           scope, * top, * zero;
   int            Acol, Aii, Ainb1D, AisR, AisRow, Ajj, Ald, AmyprocD, AmyprocR,
                  AnD, AnbD, AnpD, AnprocsD, AprocD, AprocR, Aroc, Arow, Bcol,
                  Bii, BisR, BisRow, Bjj, Bld, Bm, BmyprocD, BmyprocR, Bn,
                  BnprocsD, BprocR, Broc, Brow, RRorCC, ctxt, izero=0, k, kbb,
                  kk, kn, ktmp, mycol, mydist, myproc, myrow, npcol, nprow, p,
                  size;
   MMADD_T        add;
   TZPAD_T        pad;
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
      AnD      = N;     AnbD     = DESCA[NB_]; Ald      = DESCA[LLD_];
      AprocD   = Acol;  AprocR   = Arow;
      AmyprocD = mycol; AmyprocR = myrow;      AnprocsD = npcol;
      AisR     = ( ( Arow == -1 ) || ( nprow == 1 ) );
      Ainb1D   = PB_Cfirstnb( AnD, JA, DESCA[INB_], AnbD );
   }
   else
   {
      AnD      = M;     AnbD     = DESCA[MB_]; Ald      = DESCA[LLD_];
      AprocD   = Arow;  AprocR   = Acol;
      AmyprocD = myrow; AmyprocR = mycol;      AnprocsD = nprow;
      AisR     = ( ( Acol == -1 ) || ( npcol == 1 ) );
      Ainb1D   = PB_Cfirstnb( AnD, IA, DESCA[IMB_], AnbD );
   }
/*
*  Retrieve sub( B )'s local information: Bii, Bjj, Brow, Bcol ...
*/
   PB_Cinfog2l( IB, JB, DESCB, nprow, npcol, myrow, mycol, &Bii, &Bjj,
                &Brow, &Bcol );
   if( ( BisRow = ( Mupcase( BROC[0] ) == CROW ) ) != 0 )
   {
      Bld      = DESCB[LLD_];
      BmyprocD = mycol;       BnprocsD = npcol;
      BprocR   = Brow;        BmyprocR = myrow;
      BisR     = ( ( BprocR == -1 ) || ( nprow == 1 ) );
   }
   else
   {
      Bld      = DESCB[LLD_];
      BmyprocD = myrow;       BnprocsD = nprow;
      BprocR   = Bcol;        BmyprocR = mycol;
      BisR     = ( ( BprocR == -1 ) || ( npcol == 1 ) );
   }
/*
*  Are sub( A ) and sub( B ) both row or column vectors ?
*/
   RRorCC = ( ( AisRow && BisRow ) || ( !( AisRow ) && !( BisRow ) ) );
/*
*  Select the local add routine accordingly
*/
   size = TYPE->size;
/*
*  sub( A ) is distributed and sub( B ) is not distributed
*/
   if( !( BisR ) )
   {
/*
*  sub( B ) is not replicated. Since this operation is local if sub( B ) and
*  sub( A ) are both row or column vectors, choose AprocR = BprocR when RRorCC,
*  and AprocR = 0 otherwise.
*/
      if( AisR ) { AprocR = ( ( RRorCC ) ? BprocR : 0 ); }
/*
*  Now, it is just like sub( A ) is not replicated, this information however is
*  kept in AisR for later use.
*/
      if( ( AmyprocR == AprocR ) || ( BmyprocR == BprocR ) )
      {
         if( RRorCC )
         {
/*
*  sub( A ) and sub( B ) are both row or column vectors
*/
            zero = TYPE->zero;
            if( Mupcase( CONJUG[0] ) != CNOCONJG ) add = TYPE->Fmmcadd;
            else                                   add = TYPE->Fmmadd;
            pad = TYPE->Ftzpad;

            AnpD = PB_Cnumroc( AnD, 0, Ainb1D, AnbD, AmyprocD, AprocD,
                               AnprocsD );
/*
*  sub( A ) and sub( B ) are in the same process row or column
*/
            if( AprocR == BprocR )
            {
/*
*  In each process, the distributed part of sub( A ) is added to sub( B ). In
*  the other processes, this replicated of sub( B ) is set to zero for later
*  reduction.
*/
               if( AnpD > 0 )
               {
                  Aroc = AprocD;
                  if( BisRow ) { kk = Ajj; ktmp = JB + N; kn = JB + Ainb1D; }
                  else         { kk = Aii; ktmp = IB + M; kn = IB + Ainb1D; }

                  if( AmyprocD == Aroc )
                  {
                     if( BisRow )
                        add( &M, &Ainb1D, ALPHA, Mptr( A, Aii, Ajj, Ald, size ),
                             &Ald, BETA, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
                     else
                        add( &Ainb1D, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ),
                             &Ald, BETA, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
                     kk += Ainb1D;
                  }
                  else
                  {
                     if( BisRow )
                        pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &M, &Ainb1D,
                             &izero, zero, zero, Mptr( B, Bii, Bjj, Bld, size ),
                             &Bld );
                     else
                        pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Ainb1D, &N,
                             &izero, zero, zero, Mptr( B, Bii, Bjj, Bld, size ),
                             &Bld );
                  }
                  Aroc = MModAdd1( Aroc, AnprocsD );

                  for( k = kn; k < ktmp; k += AnbD )
                  {
                     kbb = ktmp - k; kbb = MIN( kbb, AnbD );

                     if( AmyprocD == Aroc )
                     {
                        if( BisRow )
                           add( &M, &kbb, ALPHA, Mptr( A, Aii, kk, Ald, size ),
                                &Ald, BETA, Mptr( B, Bii, k, Bld, size ),
                                &Bld );
                        else
                           add( &kbb, &N, ALPHA, Mptr( A, kk, Ajj, Ald, size ),
                                &Ald, BETA, Mptr( B, k, Bjj, Bld, size ),
                                &Bld );
                        kk += kbb;
                     }
                     else
                     {
                        if( BisRow )
                           pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &M, &kbb,
                                &izero, zero, zero, Mptr( B, Bii, k, Bld,
                                size ), &Bld );
                        else
                           pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &kbb, &N,
                                &izero, zero, zero, Mptr( B, k, Bjj, Bld,
                                size ), &Bld );
                     }
                     Aroc = MModAdd1( Aroc, AnprocsD );
                  }
               }
               else
               {
/*
*  If I don't own any entries of sub( A ), then zero the entire sub( B )
*  residing in this process.
*/
                  pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &M, &N, &izero,
                       zero, zero, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
               }
/*
*  Replicate locally scattered sub( B ) by reducing it
*/
               scope = ( BisRow ? CROW : CCOLUMN );
               top   = PB_Ctop( &ctxt, COMBINE, &scope, TOP_GET );
               TYPE->Cgsum2d( ctxt, &scope, top, M, N, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld, -1, 0 );
            }
            else
            {
/*
*  sub( A ) and sub( B ) are in a different process row or column
*/
               if( AmyprocR == AprocR )
               {
/*
*  If I own a piece of sub( A ), then send it to the corresponding process row
*  or column where sub( B ) resides.
*/
                  if( AnpD > 0 )
                  {
                     if( AisRow )
                        TYPE->Cgesd2d( ctxt, M, AnpD, Mptr( A, Aii, Ajj, Ald,
                                       size ), Ald, BprocR, BmyprocD );
                     else
                        TYPE->Cgesd2d( ctxt, AnpD, N, Mptr( A, Aii, Ajj, Ald,
                                       size ), Ald, BmyprocD, BprocR );
                  }
               }

               if( BmyprocR == BprocR )
               {
/*
*  If I own sub( B ), then receive and unpack distributed part of sub( A ) that
*  should be added to sub( B ). Combine the results.
*/
                  if( AnpD > 0 )
                  {
                     if( BisRow )
                     {
                        ktmp = JB + N;
                        kn   = JB + Ainb1D;
                        buf = PB_Cmalloc( M * AnpD * size );
                        TYPE->Cgerv2d( ctxt, M, AnpD, buf,    M, AprocR,
                                       AmyprocD );
                     }
                     else
                     {
                        ktmp = IB + M;
                        kn   = IB + Ainb1D;
                        buf = PB_Cmalloc( AnpD * N * size );
                        TYPE->Cgerv2d( ctxt, AnpD, N, buf, AnpD, AmyprocD,
                                       AprocR );
                     }
                     Aroc = AprocD;
                     kk   = 0;

                     if( AmyprocD == Aroc )
                     {
                        if( BisRow )
                           add( &M, &Ainb1D, ALPHA, buf,    &M, BETA, Mptr( B,
                                Bii, Bjj, Bld, size ), &Bld );
                        else
                           add( &Ainb1D, &N, ALPHA, buf, &AnpD, BETA, Mptr( B,
                                Bii, Bjj, Bld, size ), &Bld );
                        kk += Ainb1D;
                     }
                     else
                     {
                        if( BisRow )
                           pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &M,
                                &Ainb1D, &izero, zero, zero, Mptr( B, Bii, Bjj,
                                Bld, size ), &Bld );
                        else
                           pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Ainb1D,
                                &N, &izero, zero, zero, Mptr( B, Bii, Bjj, Bld,
                                size ), &Bld );
                     }
                     Aroc = MModAdd1( Aroc, AnprocsD );

                     for( k = kn; k < ktmp; k += AnbD )
                     {
                        kbb = ktmp - k; kbb = MIN( kbb, AnbD );

                        if( AmyprocD == Aroc )
                        {
                           if( BisRow )
                              add( &M, &kbb, ALPHA, Mptr( buf, 0, kk, M, size ),
                                   &M, BETA, Mptr( B, Bii, k, Bld, size ),
                                   &Bld );
                           else
                              add( &kbb, &N, ALPHA, Mptr( buf, kk, 0, AnpD,
                                   size ), &AnpD, BETA, Mptr( B, k, Bjj, Bld,
                                   size ), &Bld );
                           kk += kbb;
                        }
                        else
                        {
                           if( BisRow )
                              pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &M,
                                   &kbb, &izero, zero, zero, Mptr( B, Bii, k,
                                   Bld, size ), &Bld );
                           else
                              pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &kbb,
                                   &N, &izero, zero, zero, Mptr( B, k, Bjj, Bld,
                                   size ), &Bld );
                        }
                        Aroc = MModAdd1( Aroc, AnprocsD );
                     }
                     if( buf ) free( buf );
                  }
                  else
                  {
/*
*  If I don't own any entries of sub( A ), then zero the entire sub( B )
*  residing in this process.
*/
                     pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &M, &N, &izero,
                          zero, zero, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
                  }
/*
*  Replicate locally scattered sub( B ) by reducing it
*/
                  scope = ( BisRow ? CROW : CCOLUMN );
                  top   = PB_Ctop( &ctxt, COMBINE, &scope, TOP_GET );
                  TYPE->Cgsum2d( ctxt, &scope, top, M, N, Mptr( B, Bii, Bjj,
                                 Bld, size ), Bld, -1, 0 );
               }
            }
         }
         else
         {
/*
*  sub( A ) and sub( B ) are not both row or column vectors
*/
            zero = TYPE->zero;
            if( Mupcase( CONJUG[0] ) != CNOCONJG ) add = TYPE->Fmmtcadd;
            else                                   add = TYPE->Fmmtadd;
            pad = TYPE->Ftzpad;

            Broc = 0;
            if( BisRow ) { ktmp = JB + M; kn = JB + Ainb1D; }
            else         { ktmp = IB + N; kn = IB + Ainb1D; }
/*
*  Loop over the processes in which sub( A ) resides, for each process find the
*  next process Xroc. Exchange and add the data.
*/
            for( p = 0; p < AnprocsD; p++ )
            {
               mydist = MModSub(      p, AprocD, AnprocsD );
               myproc = MModAdd( AprocD, mydist, AnprocsD );

               if( ( BprocR == p ) && ( AprocR == Broc ) )
               {
                  if( BmyprocR == p )
                  {
/*
*  local add at the intersection of the process cross
*/
                     AnpD = PB_Cnumroc( AnD, 0, Ainb1D, AnbD, p, AprocD,
                                        AnprocsD );
                     if( AnpD > 0 )
                     {
                        Aroc = AprocD;
                        kk   = ( BisRow ? Aii : Ajj );

                        if( myproc == Aroc )
                        {
                           if( BmyprocD == Broc )
                           {
                              if( AisRow )
                                 add( &M, &Ainb1D, ALPHA, Mptr( A, Aii, Ajj,
                                      Ald, size ), &Ald, BETA, Mptr( B, Bii,
                                      Bjj, Bld, size ), &Bld );
                              else
                                 add( &Ainb1D, &N, ALPHA, Mptr( A, Aii, Ajj,
                                      Ald, size ), &Ald, BETA, Mptr( B, Bii,
                                      Bjj, Bld, size ), &Bld );
                              kk += Ainb1D;
                           }
                           else
                           {
                              if( BisRow )
                                 pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &N,
                                      &Ainb1D, &izero, zero, zero, Mptr( B, Bii,
                                      Bjj, Bld, size ), &Bld );
                              else
                                 pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ),
                                      &Ainb1D, &M, &izero, zero, zero, Mptr( B,
                                      Bii, Bjj, Bld, size ), &Bld );
                           }
                        }
                        Aroc = MModAdd1( Aroc, AnprocsD );

                        for( k = kn; k < ktmp; k += AnbD )
                        {
                           kbb = ktmp - k; kbb = MIN( kbb, AnbD );
                           if( myproc == Aroc )
                           {
                              if( BmyprocD == Broc )
                              {
                                 if( AisRow )
                                    add( &M, &kbb, ALPHA, Mptr( A, Aii, kk, Ald,
                                         size ), &Ald, BETA, Mptr( B, k, Bjj,
                                         Bld, size ), &Bld );
                                 else
                                    add( &kbb, &N, ALPHA, Mptr( A, kk, Ajj, Ald,
                                         size ), &Ald, BETA, Mptr( B, Bii, k,
                                         Bld, size ), &Bld );
                                 kk += kbb;
                              }
                              else
                              {
                                 if( BisRow )
                                    pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ),
                                         &N, &kbb, &izero, zero, zero, Mptr( B,
                                         Bii, k, Bld, size ), &Bld );
                                 else
                                    pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ),
                                         &kbb, &M, &izero, zero, zero, Mptr( B,
                                         k, Bjj, Bld, size ), &Bld );
                              }
                           }
                           Aroc = MModAdd1( Aroc, AnprocsD );
                        }
                     }
                  }
               }
               else
               {
/*
*  Message exchange
*/
                  if( ( AmyprocR == AprocR ) && ( AmyprocD == p ) )
                  {
                     AnpD = PB_Cnumroc( AnD, 0, Ainb1D, AnbD, p, AprocD,
                                        AnprocsD );
                     if( AnpD > 0 )
                     {
                        if( AisRow )
                           TYPE->Cgesd2d( ctxt, M, AnpD, Mptr( A, Aii, Ajj, Ald,
                                          size ), Ald, Broc, BprocR );
                        else
                           TYPE->Cgesd2d( ctxt, AnpD, N, Mptr( A, Aii, Ajj, Ald,
                                          size ), Ald, BprocR, Broc );
                     }
                  }

                  if( BmyprocR == BprocR )
                  {
                     AnpD = PB_Cnumroc( AnD, 0, Ainb1D, AnbD, p, AprocD,
                                        AnprocsD );
                     if( AnpD > 0 )
                     {
                        Aroc = AprocD;
                        kk   = 0;

                        if( BmyprocD == Broc )
                        {
                           if( AisRow )
                           {
                              buf = PB_Cmalloc( M * AnpD * size );
                              TYPE->Cgerv2d( ctxt, M, AnpD, buf, M, AprocR, p );
                           }
                           else
                           {
                              buf = PB_Cmalloc( AnpD * N * size );
                              TYPE->Cgerv2d( ctxt, AnpD, N, buf, AnpD, p,
                                             AprocR );
                           }
                        }

                        if( myproc == Aroc )
                        {
                           if( BmyprocD == Broc )
                           {
                              if( AisRow )
                                 add( &M, &Ainb1D, ALPHA, buf,    &M, BETA,
                                      Mptr( B, Bii, Bjj, Bld, size ), &Bld );
                              else
                                 add( &Ainb1D, &N, ALPHA, buf, &AnpD, BETA,
                                      Mptr( B, Bii, Bjj, Bld, size ), &Bld );
                              kk += Ainb1D;
                           }
                           else
                           {
                              if( BisRow )
                                 pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &N,
                                      &Ainb1D, &izero, zero, zero, Mptr( B, Bii,
                                      Bjj, Bld, size ), &Bld );
                              else
                                 pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ),
                                      &Ainb1D, &M, &izero, zero, zero, Mptr( B,
                                      Bii, Bjj, Bld, size ), &Bld );
                           }
                        }
                        Aroc = MModAdd1( Aroc, AnprocsD );

                        for( k = kn; k < ktmp; k += AnbD )
                        {
                           kbb = ktmp - k; kbb = MIN( kbb, AnbD );
                           if( myproc == Aroc )
                           {
                              if( BmyprocD == Broc )
                              {
                                 if( AisRow )
                                    add( &M, &kbb, ALPHA, Mptr( buf, 0, kk, M,
                                         size ),    &M, BETA, Mptr( B, k, Bjj,
                                         Bld, size ), &Bld );
                                 else
                                    add( &kbb, &N, ALPHA, Mptr( buf, kk, 0,
                                         AnpD, size ), &AnpD, BETA, Mptr( B,
                                         Bii, k, Bld, size ), &Bld );
                                 kk += kbb;
                              }
                              else
                              {
                                 if( BisRow )
                                    pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ),
                                         &N, &kbb, &izero, zero, zero, Mptr( B,
                                         Bii, k, Bld, size ), &Bld );
                                 else
                                    pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ),
                                         &kbb, &M, &izero, zero, zero, Mptr( B,
                                         k, Bjj, Bld, size ), &Bld );
                              }
                           }
                           Aroc = MModAdd1( Aroc, AnprocsD );
                        }
                        if( ( BmyprocD == Broc ) && ( buf ) ) free( buf );
                     }
                  }
               }
               Broc = MModAdd1( Broc, BnprocsD );
            }

            if( BmyprocR == BprocR )
            {
/*
*  Replicate locally scattered sub( B ) by reducing it
*/
               scope = ( BisRow ? CROW : CCOLUMN );
               top   = PB_Ctop( &ctxt, COMBINE, &scope, TOP_GET );
               TYPE->Cgsum2d( ctxt, &scope, top, N, M, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld, -1, 0 );
            }
         }
      }

      if( BisR )
      {
/*
*  Replicate sub( B )
*/
         if( BisRow )
         {
            if( AisRow ) { Bm = M; Bn = N; }
            else         { Bm = N; Bn = M; }
            top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
            if( BmyprocR == BprocR )
               TYPE->Cgebs2d( ctxt, COLUMN, top, Bm, Bn, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld );
            else
               TYPE->Cgebr2d( ctxt, COLUMN, top, Bm, Bn, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld, BprocR, BmyprocD );
         }
         else
         {
            if( AisRow ) { Bm = N; Bn = M; }
            else         { Bm = M; Bn = N; }
            top = PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
            if( BmyprocR == BprocR )
               TYPE->Cgebs2d( ctxt, ROW,    top, Bm, Bn, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld );
            else
               TYPE->Cgebr2d( ctxt, ROW,    top, Bm, Bn, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld, BmyprocD, BprocR );
         }
      }
   }
   else
   {
/*
*  sub( B ) is replicated in every process. Add the data in process row or
*  column AprocR when sub( A ) is not replicated and in every process otherwise.
*/
      if( AisR || ( AmyprocR == AprocR ) )
      {
         zero = TYPE->zero;
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
         pad = TYPE->Ftzpad;

         AnpD = PB_Cnumroc( AnD, 0, Ainb1D, AnbD, AmyprocD, AprocD, AnprocsD );
         if( AnpD > 0 )
         {
            Aroc = AprocD;
            kk   = ( AisRow ? Ajj : Aii );

            if( BisRow ) { ktmp = JB + ( RRorCC ? N : M ); kn = JB + Ainb1D; }
            else         { ktmp = IB + ( RRorCC ? M : N ); kn = IB + Ainb1D; }

            if( AmyprocD == Aroc )
            {
               if( AisRow )
                  add( &M, &Ainb1D, ALPHA, Mptr( A, Aii, Ajj, Ald, size ), &Ald,
                       BETA, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
               else
                  add( &Ainb1D, &N, ALPHA, Mptr( A, Aii, Ajj, Ald, size ), &Ald,
                       BETA, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
               kk += Ainb1D;
            }
            else
            {
               if( RRorCC )
               {
                  if( AisRow ) { Bm = M; Bn = Ainb1D; }
                  else         { Bm = Ainb1D; Bn = N; }
               }
               else
               {
                  if( AisRow ) { Bm = Ainb1D; Bn = M; }
                  else         { Bm = N; Bn = Ainb1D; }
               }
               pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Bm, &Bn, &izero,
                    zero, zero, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
            }
            Aroc = MModAdd1( Aroc, AnprocsD );

            for( k = kn; k < ktmp; k += AnbD )
            {
               kbb = ktmp - k; kbb = MIN( kbb, AnbD );

               if( BisRow ) { buf = Mptr( B, Bii, k, Bld, size ); }
               else         { buf = Mptr( B, k, Bjj, Bld, size ); }

               if( AmyprocD == Aroc )
               {
                  if( AisRow )
                     add( &M, &kbb, ALPHA, Mptr( A, Aii, kk, Ald, size ), &Ald,
                          BETA, buf, &Bld );
                  else
                     add( &kbb, &N, ALPHA, Mptr( A, kk, Ajj, Ald, size ), &Ald,
                          BETA, buf, &Bld );
                  kk += kbb;
               }
               else
               {
                  if( RRorCC )
                  {
                     if( AisRow ) { Bm = M; Bn = kbb; }
                     else         { Bm = kbb; Bn = N; }
                  }
                  else
                  {
                     if( AisRow ) { Bm = kbb; Bn = M; }
                     else         { Bm = N; Bn = kbb; }
                  }
                  pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Bm, &Bn, &izero,
                       zero, zero, buf, &Bld );
               }
               Aroc = MModAdd1( Aroc, AnprocsD );
            }
         }
         else
         {
            if( RRorCC )
               pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &M, &N, &izero, zero,
                    zero, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
            else
               pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &N, &M, &izero, zero,
                    zero, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
         }
/*
*  Replicate locally scattered sub( B ) by reducing it in the process scope of
*  sub( A )
*/
         scope = ( AisRow ? CROW : CCOLUMN );
         top   = PB_Ctop( &ctxt, COMBINE, &scope, TOP_GET );
         if( RRorCC )
            TYPE->Cgsum2d( ctxt, &scope, top, M, N, Mptr( B, Bii, Bjj, Bld,
                           size ), Bld, -1, 0 );
         else
            TYPE->Cgsum2d( ctxt, &scope, top, N, M, Mptr( B, Bii, Bjj, Bld,
                           size ), Bld, -1, 0 );
      }

      if( !AisR )
      {
/*
*  If sub( A ) is not replicated, then broadcast the result to the other pro-
*  cesses that own a piece of sub( B ), but were not involved in the above
*  addition operation.
*/
         if( RRorCC ) { Bm = M; Bn = N; }
         else         { Bm = N; Bn = M; }

         if( AisRow )
         {
            top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
            if( AmyprocR == AprocR )
               TYPE->Cgebs2d( ctxt, COLUMN, top, Bm, Bn, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld );
            else
               TYPE->Cgebr2d( ctxt, COLUMN, top, Bm, Bn, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld, AprocR, AmyprocD );
         }
         else
         {
            top = PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
            if( AmyprocR == AprocR )
               TYPE->Cgebs2d( ctxt, ROW,    top, Bm, Bn, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld );
            else
               TYPE->Cgebr2d( ctxt, ROW,    top, Bm, Bn, Mptr( B, Bii, Bjj, Bld,
                              size ), Bld, AmyprocD, AprocR );
         }
      }
   }
/*
*  End of PB_CpaxpbyDN
*/
}
