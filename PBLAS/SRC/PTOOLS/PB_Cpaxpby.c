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
void PB_Cpaxpby( PBTYP_T * TYPE, char * CONJUG, Int M, Int N,
                 char * ALPHA,
                 char * A, Int IA, Int JA, Int * DESCA, char * AROC,
                 char * BETA,
                 char * B, Int IB, Int JB, Int * DESCB, char * BROC )
#else
void PB_Cpaxpby( TYPE, CONJUG, M, N, ALPHA, A, IA, JA, DESCA, AROC,
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
*  PB_Cpaxpby adds one submatrix to another,
*
*     sub( B ) := beta * sub( B ) + alpha * sub( A ), or,
*
*     sub( B ) := beta * sub( B ) + alpha * conjg( sub( A ) ),
*
*  where both submatrices are distributed along one dimension;  sub( A )
*  always  denotes  A(IA:IA+M-1,JA:JA+N-1).  When  AROC  is  'R'  or 'r'
*  sub( A ) is  distributed  along  a  process  row, otherwise  sub( A )
*  is distributed along a process column.  When sub( A )  is distributed
*  along a process row and BROC is 'R' or 'r' or sub( A ) is distributed
*  along a process column and BROC is 'C' or 'c', then sub( B )  denotes
*  B(IB:IB+M-1,JB:JB+N-1), and B(IB:IB+N-1,JB:JB+M-1) otherwise.
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
   char           ascope, bscope, * buf = NULL, * one, * top, tran, * zero;
   Int            Acol, Aii, AinbD, Ainb1D, AisD, AisR, AisRow, AiD, Ajj, Ald,
                  AmyprocD, AmyprocR, AnbD, AnD, AnR, AnpD, AnprocsD, AnprocsR,
                  AprocD, AprocR, Aroc, Arow, Bcol, Bii, BinbD, Binb1D, BisD,
                  BisR, BisRow, BiD, Bjj, Bld, BmyprocD, BmyprocR, BnbD, BnD,
                  BnR, BnpD, BnprocsD, BnprocsR, BprocD, BprocR, Broc, Brow,
                  BsrcD, OneBlock, OneDgrid, RRorCC, Square, cdst, csrc, ctxt,
                  dst, gcdPQ, k, l, lcmPQ, lcmb, ma, mb, mycol, myrow, na, nb,
                  npcol, npq, nprow, p, q, rdst, rsrc, size, src;
   PB_VM_T        VM;
   MMADD_T        add;
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
   Cblacs_gridinfo( ( ctxt = DESCA[ CTXT_ ] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Determine if sub( A ) is distributed or not
*/
   if( ( AisRow = ( Mupcase( AROC[0] ) == CROW ) ) != 0 )
      AisD = ( ( DESCA[CSRC_] >= 0 ) && ( ( AnprocsD = npcol ) > 1 ) );
   else
      AisD = ( ( DESCA[RSRC_] >= 0 ) && ( ( AnprocsD = nprow ) > 1 ) );
/*
*  Determine if sub( B ) is distributed or not
*/
   if( ( BisRow = ( Mupcase( BROC[0] ) == CROW ) ) != 0 )
      BisD = ( ( DESCB[CSRC_] >= 0 ) && ( ( BnprocsD = npcol ) > 1 ) );
   else
      BisD = ( ( DESCB[RSRC_] >= 0 ) && ( ( BnprocsD = nprow ) > 1 ) );
/*
*  AisD && BisD <=> both operands are indeed distributed
*/
   if( AisD && BisD )
   {
/*
*  Retrieve sub( A )'s local information: Aii, Ajj, Arow, Acol ...
*/
      PB_Cinfog2l( IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aii, &Ajj, &Arow,
                   &Acol );
      if( AisRow )
      {
         AinbD  = DESCA[INB_]; AnbD = DESCA[NB_]; Ald = DESCA[LLD_];
         AiD    = JA;    AnD    =  N;      AnR    =  M;
         AprocD = Acol; AmyprocD = mycol;
         AprocR = Arow; AmyprocR = myrow; AnprocsR = nprow;
         AisR   = ( ( DESCA[ RSRC_ ] == -1 ) || ( AnprocsR == 1 ) );
      }
      else
      {
         AinbD  = DESCA[IMB_]; AnbD = DESCA[MB_]; Ald = DESCA[LLD_];
         AiD    = IA;   AnD    =  M;      AnR    =  N;
         AprocD = Arow; AmyprocD = myrow;
         AprocR = Acol; AmyprocR = mycol; AnprocsR = npcol;
         AisR   = ( ( DESCA[ CSRC_ ] == -1 ) || ( AnprocsR == 1 ) );
      }
      Ainb1D = PB_Cfirstnb( AnD, AiD, AinbD, AnbD );
/*
*  Retrieve sub( B )'s local information: Bii, Bjj, Brow, Bcol ...
*/
      PB_Cinfog2l( IB, JB, DESCB, nprow, npcol, myrow, mycol, &Bii, &Bjj, &Brow,
                   &Bcol );
      if( BisRow )
      {
         BinbD  = DESCB[ INB_  ]; BnbD   = DESCB[ NB_   ];
         BsrcD  = DESCB[ CSRC_ ]; Bld    = DESCB[ LLD_  ];
         BiD    = JB;
         if( AisRow ) { BnD = N; BnR = M; } else { BnD = M; BnR = N; }
         BprocD = Bcol; BmyprocD = mycol;
         BprocR = Brow; BmyprocR = myrow; BnprocsR = nprow;
         BisR   = ( ( DESCB[ RSRC_ ] == -1 ) || ( BnprocsR == 1 ) );
      }
      else
      {
         BinbD  = DESCB[ IMB_  ]; BnbD   = DESCB[ MB_   ];
         BsrcD  = DESCB[ RSRC_ ]; Bld    = DESCB[ LLD_  ];
         BiD    = IB;
         if( AisRow ) { BnD = N; BnR = M; } else { BnD = M; BnR = N; }
         BprocD = Brow; BmyprocD = myrow;
         BprocR = Bcol; BmyprocR = mycol; BnprocsR = npcol;
         BisR   = ( ( DESCB[ CSRC_ ] == -1 ) || ( BnprocsR == 1 ) );
      }
      Binb1D = PB_Cfirstnb( BnD, BiD, BinbD, BnbD );
/*
*  Are sub( A ) and sub( B ) both row or column vectors ?
*/
      RRorCC = ( ( AisRow && BisRow ) || ( !( AisRow ) && !( BisRow ) ) );
/*
*  Do sub( A ) and sub( B ) span more than one process ?
*/
      OneDgrid = ( ( AnprocsD ==   1 ) && ( BnprocsD ==   1 ) );
      OneBlock = ( ( Ainb1D   >= AnD ) && ( Binb1D   >= BnD ) );
/*
*  Are sub( A ) and sub( B ) distributed in the same manner ?
*/
      Square   = ( ( Ainb1D   ==   Binb1D ) && ( AnbD == BnbD ) &&
                   ( AnprocsD == BnprocsD ) );

      if( !( AisR ) )
      {
/*
*  sub( A ) is distributed but not replicated
*/
         if( BisR )
         {
/*
*  If sub( A ) is not replicated, but sub( B ) is, a process row or column
*  BprocR need to be selected. It will contain the non-replicated vector to
*  add sub( A ) to.
*/
            if( RRorCC )
            {
/*
*  sub( A ) and sub( B ) are both row or column vectors
*/
               if( ( OneDgrid || OneBlock || Square ) && ( AprocD == BprocD ) )
               {
/*
*  sub( A ) and sub( B ) start in the same process row or column AprocD=BprocD.
*  Enforce a purely local operation by choosing BprocR to be equal to AprocR.
*/
                  BprocR = AprocR;
               }
               else
               {
/*
*  Otherwise, communication has to occur, so choose the next process row or
*  column for BprocR to maximize the number of links, i.e reduce contention.
*/
                  BprocR = MModAdd1( AprocR, AnprocsR );
               }
            }
            else
            {
/*
*  sub( A ) and sub( B ) are distributed in orthogonal directions, what is
*  chosen for YprocR does not really matter. Select the process origin.
*/
               BprocR = AprocD;
            }
         }
         else
         {
/*
*  Neither sub( A ) nor sub( B ) are replicated. If I am not in process row or
*  column AprocR and not in process row or column BprocR, then quick return.
*/
            if( ( AmyprocR != AprocR ) && ( BmyprocR != BprocR ) )
               return;
         }
      }
      else
      {
/*
*  sub( A ) is distributed and replicated (so no quick return possible)
*/
         if( BisR )
         {
/*
*  sub( B ) is distributed and replicated as well
*/
            if( RRorCC )
            {
/*
*  sub( A ) and sub( B ) are both row or column vectors
*/
               if( ( OneDgrid || OneBlock || Square ) && ( AprocD == BprocD ) )
               {
/*
*  sub( A ) and sub( B ) start in the same process row or column AprocD=BprocD.
*  Enforce a purely local operation by choosing AprocR and BprocR to be equal
*  to zero.
*/
                  AprocR = BprocR = 0;
               }
               else
               {
/*
*  Otherwise, communication has to occur, so select BprocR to be zero and the
*  next process row or column for AprocR in order to maximize the number of
*  used links, i.e reduce contention.
*/
                  BprocR = 0;
                  AprocR = MModAdd1( BprocR, BnprocsR );
               }
            }
            else
            {
/*
*  sub( A ) and sub( B ) are distributed in orthogonal directions, select the
*  origin processes.
*/
               AprocR = BprocD;
               BprocR = AprocD;
            }
         }
         else
         {
/*
*  sub( B ) is distributed, but not replicated
*/
            if( RRorCC )
            {
/*
*  sub( A ) and sub( B ) are both row or column vectors
*/
               if( ( OneDgrid || OneBlock || Square ) && ( AprocD == BprocD ) )
               {
/*
*  sub( A ) and sub( B ) start in the same process row or column AprocD=BprocD.
*  Enforce a purely local operation by choosing AprocR to be equal to BprocR.
*/
                  AprocR = BprocR;
                  if( ( AmyprocR != AprocR ) && ( BmyprocR != BprocR ) ) return;
               }
               else
               {
/*
*  Otherwise, communication has to occur, so choose the next process row or
*  column for AprocR to maximize the number of links, i.e reduce contention.
*/
                  AprocR = MModAdd1( BprocR, BnprocsR );
               }
            }
            else
            {
/*
*  sub( A ) and sub( B ) are distributed in orthogonal directions, what is
*  chosen for AprocR does not really matter. Select the origin process.
*/
               AprocR = BprocD;
               if( ( OneDgrid || OneBlock || Square ) &&
                   ( AmyprocR != AprocR ) && ( BmyprocR != BprocR ) ) return;
            }
         }
      }
/*
*  Even if sub( A ) and/or sub( B ) are replicated, only two process row or
*  column are active, namely AprocR and BprocR. If any of those operands is
*  replicated, broadcast will occur (unless there is an easy way out).
*/
      size = TYPE->size;
/*
*  A purely local operation occurs iff the operands start in the same process
*  and, if either the grid is mono-dimensional or there is a single local block
*  to be added or if both operands are aligned.
*/
      if( ( (    RRorCC   && ( AprocD == BprocD ) &&
                        ( AisR || BisR || ( AprocR == BprocR ) ) ) ||
            ( !( RRorCC ) && ( BisR || ( AprocD == BprocR ) ) &&
                             ( AisR || ( AprocR == BprocD ) ) ) ) &&
          ( OneDgrid || OneBlock || ( RRorCC && Square ) ) )
      {
         if( ( !AisR && ( AmyprocR == AprocR ) ) ||
             ( AisR && ( BisR || BmyprocR == BprocR ) ) )
         {
            AnpD = PB_Cnumroc( AnD, 0, Ainb1D, AnbD, AmyprocD, AprocD,
                               AnprocsD );
            BnpD = PB_Cnumroc( BnD, 0, Binb1D, BnbD, BmyprocD, BprocD,
                               BnprocsD );
            if( ( AnpD > 0 ) && ( BnpD > 0 ) )
            {
/*
*  Select the local add routine accordingly to RRorCC
*/
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
/*
*  Local addition
*/
               if( AisRow )
                  add( &AnR, &AnpD, ALPHA, Mptr( A, Aii, Ajj, Ald, size ), &Ald,
                       BETA, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
               else
                  add( &AnpD, &AnR, ALPHA, Mptr( A, Aii, Ajj, Ald, size ), &Ald,
                       BETA, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
            }
         }
         if( RRorCC && AisR && BisR ) return;
      }
      else if( ( RRorCC && OneDgrid ) || OneBlock || Square )
      {
/*
*  Otherwise, it may be possible to add the distributed vectors in a single
*  message exchange iff the grid is mono-dimensional and the operands are
*  distributed in the same direction, or there is just one block to be exchanged
*  or if both operands are similarly distributed in their respective direction.
*/
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

         if( ( AisR && BisR ) || ( AmyprocR == AprocR ) )
         {

            AnpD = PB_Cnumroc( AnD, 0, Ainb1D, AnbD, AmyprocD, AprocD,
                               AnprocsD );
            if( AnpD > 0 )
            {
               dst = BprocD + MModSub( AmyprocD, AprocD, AnprocsD );
               dst = MPosMod( dst, BnprocsD );
               if( AisRow ) { ma = AnR; na = AnpD; }
               else         { ma = AnpD; na = AnR; }
               if( !( AisR && BisR ) )
               {
                  if( BisRow ) { rdst = BprocR; cdst = dst; }
                  else         { rdst = dst; cdst = BprocR; }
               }
               else
               {
                  if( BisRow )
                  {
                     if( !AisRow ) { rdst = AmyprocR; }
                     else { rdst = MModAdd1( BmyprocR, BnprocsR ); }
                     cdst = dst;
                  }
                  else
                  {
                     rdst = dst;
                     if( AisRow ) { cdst = AmyprocR; }
                     else { cdst = MModAdd1( BmyprocR, BnprocsR ); }
                  }
               }
               if( ( myrow == rdst ) && ( mycol == cdst ) )
               {
                  add( &ma, &na, ALPHA, Mptr( A, Aii, Ajj, Ald, size ), &Ald,
                       BETA, Mptr( B, Bii, Bjj, Bld, size ), &Bld );
               }
               else
               {
                  TYPE->Cgesd2d( ctxt, ma, na, Mptr( A, Aii, Ajj, Ald, size ),
                                 Ald, rdst, cdst );
               }
            }
         }

         if( ( AisR && BisR ) || ( BmyprocR == BprocR ) )
         {
            BnpD = PB_Cnumroc( BnD, 0, Binb1D, BnbD, BmyprocD, BprocD,
                               BnprocsD );
            if( BnpD > 0 )
            {
               src = AprocD + MModSub( BmyprocD, BprocD, BnprocsD );
               src = MPosMod( src, AnprocsD );
               if( AisRow ) { ma = BnR; na = BnpD; }
               else         { ma = BnpD; na = BnR; }
               if( !( AisR && BisR ) )
               {
                  if( AisRow ) { rsrc = AprocR; csrc = src; }
                  else         { rsrc = src; csrc = AprocR; }
               }
               else
               {
                  if( AisRow )
                  {
                     if( !BisRow ) { rsrc = BmyprocR; }
                     else { rsrc = MModSub1( AmyprocR, AnprocsR ); }
                     csrc = src;
                  }
                  else
                  {
                     rsrc = src;
                     if( BisRow ) { csrc = BmyprocR; }
                     else { csrc = MModSub1( AmyprocR, AnprocsR ); }
                  }
               }
               if( ( myrow != rsrc ) || ( mycol != csrc ) )
               {
                  buf = PB_Cmalloc( BnpD * BnR * size );
                  TYPE->Cgerv2d( ctxt, ma, na, buf, ma, rsrc, csrc );
                  add( &ma, &na, ALPHA, buf, &ma, BETA, Mptr( B, Bii, Bjj, Bld,
                       size ), &Bld );
                  if( buf ) free( buf );
               }
            }
         }
         if( AisR && BisR ) return;
      }
      else
      {
/*
*  General case
*/
         if( RRorCC )
         {
            if( Mupcase( CONJUG[0] ) != CNOCONJG ) tran = CCONJG;
            else                                   tran = CNOTRAN;
         }
         else
         {
            if( Mupcase( CONJUG[0] ) != CNOCONJG ) tran = CCOTRAN;
            else                                   tran = CTRAN;
         }

         if( AisRow ) { ascope = CCOLUMN; ma = AnR; }
         else         { ascope = CROW;    na = AnR; }
         bscope = ( BisRow ? CCOLUMN : CROW );
         lcmb   = PB_Clcm( AnprocsD * AnbD, BnprocsD * BnbD );
         one    = TYPE->one; zero = TYPE->zero;
         gcdPQ  = PB_Cgcd( AnprocsD, BnprocsD );
         lcmPQ  = ( AnprocsD / gcdPQ ) * BnprocsD;

         for( k = 0; k < gcdPQ; k++ )
         {
            p = 0; q = k;

            for( l = 0; l < lcmPQ; l++ )
            {
               Aroc = MModAdd( AprocD, p, AnprocsD );
               Broc = MModAdd( BprocD, q, BnprocsD );

               if( ( AmyprocD == Aroc ) || ( BmyprocD == Broc ) )
               {
                  AnpD = PB_Cnumroc( AnD, 0, Ainb1D, AnbD, Aroc, AprocD,
                                     AnprocsD );
                  BnpD = PB_Cnumroc( BnD, 0, Binb1D, BnbD, Broc, BprocD,
                                     BnprocsD );
                  PB_CVMinit( &VM, 0, AnpD, BnpD, Ainb1D, Binb1D, AnbD, BnbD,
                              p, q, AnprocsD, BnprocsD, lcmb );
                  if( npq = PB_CVMnpq( &VM ) )
                  {
                     if( ( RRorCC && ( Aroc == Broc ) &&
                           ( AisR || ( AprocR == BprocR ) ) ) ||
                         ( !( RRorCC ) && ( Aroc == BprocR ) &&
                           ( AisR || ( AprocR == Broc ) ) ) )
                     {
                        if( ( BmyprocD ==  Broc ) && ( BmyprocR == BprocR ) )
                        {
                           PB_CVMloc( TYPE, &VM, ROW, &ascope, PACKING, &tran,
                                      npq, AnR, ALPHA, Mptr( A, Aii, Ajj, Ald,
                                      size ), Ald, BETA, Mptr( B, Bii, Bjj, Bld,
                                      size ), Bld );
                        }
                     }
                     else
                     {
                        if( ( AmyprocR == AprocR ) && ( AmyprocD == Aroc  ) )
                        {
                           if( AisRow ) { na = npq; } else { ma = npq; }
                           buf = PB_Cmalloc( ma * na * size );
                           PB_CVMpack( TYPE, &VM, ROW, &ascope, PACKING, NOTRAN,
                                       npq, AnR, one, Mptr( A, Aii, Ajj, Ald,
                                       size ), Ald, zero, buf, ma );
                           if( BisRow ) { rdst = BprocR; cdst = Broc; }
                           else         { rdst = Broc; cdst = BprocR; }
                           TYPE->Cgesd2d( ctxt, ma, na, buf, ma, rdst, cdst );
                           if( buf ) free ( buf );
                        }
                        if( ( BmyprocR == BprocR ) && ( BmyprocD == Broc ) )
                        {
                           if( AisRow )
                           { na = npq; rsrc = AprocR; csrc = Aroc; }
                           else
                           { ma = npq; rsrc = Aroc; csrc = AprocR; }
                           buf = PB_Cmalloc( ma * na * size );
                           TYPE->Cgerv2d( ctxt, ma, na, buf, ma, rsrc, csrc );
                           PB_CVMpack( TYPE, &VM, COLUMN, &bscope, UNPACKING,
                                       &tran, npq, AnR, BETA, Mptr( B, Bii, Bjj,
                                       Bld, size ), Bld, ALPHA, buf, ma );
                           if( buf ) free ( buf );
                        }
                     }
                  }
               }
               p = MModAdd1( p, AnprocsD );
               q = MModAdd1( q, BnprocsD );
            }
            if( AisR ) AprocR = MModAdd1( AprocR, AnprocsR );
         }
      }

      if( BisR )
      {
/*
*  Replicate sub( B )
*/
         BnpD = PB_Cnumroc( BnD, BiD, BinbD, BnbD, BmyprocD, BsrcD, BnprocsD );
         if( BnpD > 0 )
         {
            if( BisRow )
            {
               bscope = CCOLUMN;  mb   = BnR;      nb = BnpD;
               rsrc   = BprocR;   csrc = BmyprocD;
            }
            else
            {
               bscope = CROW;     mb   = BnpD;     nb = BnR;
               rsrc   = BmyprocD; csrc = BprocR;
            }
            top = PB_Ctop( &ctxt, BCAST, &bscope, TOP_GET );
            if( BmyprocR == BprocR )
            {
               TYPE->Cgebs2d( ctxt, &bscope, top, mb, nb, Mptr( B, Bii, Bjj,
                              Bld, size ), Bld );
            }
            else
            {
               TYPE->Cgebr2d( ctxt, &bscope, top, mb, nb, Mptr( B, Bii, Bjj,
                              Bld, size ), Bld, rsrc, csrc );
            }
         }
      }
   }
   else if( !( AisD ) && BisD )
   {
/*
*  sub( A ) is not distributed and sub( B ) is distributed.
*/
      PB_CpaxpbyND( TYPE, CONJUG, M, N, ALPHA, A, IA, JA, DESCA, AROC, BETA, B,
                    IB, JB, DESCB, BROC );
   }
   else if( AisD && !( BisD ) )
   {
/*
*  sub( A ) is distributed and sub( B ) is not distributed.
*/
      PB_CpaxpbyDN( TYPE, CONJUG, M, N, ALPHA, A, IA, JA, DESCA, AROC, BETA, B,
                    IB, JB, DESCB, BROC );
   }
   else
   {
/*
*  Neither sub( A ) nor sub( B ) are distributed.
*/
      PB_CpaxpbyNN( TYPE, CONJUG, M, N, ALPHA, A, IA, JA, DESCA, AROC, BETA, B,
                    IB, JB, DESCB, BROC );
   }
/*
*  End of PB_Cpaxpby
*/
}
