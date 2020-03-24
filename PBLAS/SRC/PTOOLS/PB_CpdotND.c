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
void PB_CpdotND( PBTYP_T * TYPE, Int N, char * DOT,
                 char * X, Int IX, Int JX, Int * DESCX, Int INCX,
                 char * Y, Int IY, Int JY, Int * DESCY, Int INCY,
                 VVDOT_T FDOT )
#else
void PB_CpdotND( TYPE, N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY,
                 INCY, FDOT )
/*
*  .. Scalar Arguments ..
*/
   Int            INCX, INCY, IX, IY, JX, JY, N;
   char           * DOT;
   PBTYP_T        * TYPE;
   VVDOT_T        FDOT;
/*
*  .. Array Arguments ..
*/
   Int            * DESCX, * DESCY;
   char           * X, * Y;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CpdotND forms the dot product of two subvectors,
*
*     DOT := sub( X )**T * sub( Y )  or   DOT := sub( X )**H * sub( Y ),
*
*  where
*
*     sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
*                      X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X, and,
*
*     sub( Y ) denotes Y(IY,JY:JY+N-1) if INCY = M_Y,
*                      Y(IY:IY+N-1,JY) if INCY = 1 and INCY <> M_Y.
*
*  sub( X ) is assumed to be not distributed, and sub( Y ) is assumed to
*  be distributed.
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
*  N       (global input) INTEGER
*          On entry,  N  specifies the  length of the  subvectors to  be
*          multiplied. N must be at least zero.
*
*  DOT     (local output) pointer to CHAR
*          On exit, DOT  specifies the dot product of the two subvectors
*          sub( X ) and sub( Y ) only in their scope (See below for fur-
*          ther details).
*
*  X       (local input) pointer to CHAR
*          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
*          is   at  least  MAX( 1, Lr( 1, IX ) )  when  INCX = M_X   and
*          MAX( 1, Lr( 1, IX+N-1 ) )  otherwise,  and,  Kx  is  at least
*          Lc( 1, JX+N-1 )  when  INCX = M_X  and Lc( 1, JX ) otherwise.
*          Before  entry,  this  array contains the local entries of the
*          matrix X.
*
*  IX      (global input) INTEGER
*          On entry, IX  specifies X's global row index, which points to
*          the beginning of the submatrix sub( X ).
*
*  JX      (global input) INTEGER
*          On entry, JX  specifies X's global column index, which points
*          to the beginning of the submatrix sub( X ).
*
*  DESCX   (global and local input) INTEGER array
*          On entry, DESCX  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix X.
*
*  INCX    (global input) INTEGER
*          On entry,  INCX   specifies  the  global  increment  for  the
*          elements of  X.  Only two values of  INCX   are  supported in
*          this version, namely 1 and M_X. INCX  must not be zero.
*
*  Y       (local input) pointer to CHAR
*          On entry, Y is an array of dimension (LLD_Y, Ky), where LLD_Y
*          is   at  least  MAX( 1, Lr( 1, IY ) )  when  INCY = M_Y   and
*          MAX( 1, Lr( 1, IY+N-1 ) )  otherwise,  and,  Ky  is  at least
*          Lc( 1, JY+N-1 )  when  INCY = M_Y  and Lc( 1, JY ) otherwise.
*          Before  entry,  this array  contains the local entries of the
*          matrix Y.
*
*  IY      (global input) INTEGER
*          On entry, IY  specifies Y's global row index, which points to
*          the beginning of the submatrix sub( Y ).
*
*  JY      (global input) INTEGER
*          On entry, JY  specifies Y's global column index, which points
*          to the beginning of the submatrix sub( Y ).
*
*  DESCY   (global and local input) INTEGER array
*          On entry, DESCY  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix Y.
*
*  INCY    (global input) INTEGER
*          On entry,  INCY   specifies  the  global  increment  for  the
*          elements of  Y.  Only two values of  INCY   are  supported in
*          this version, namely 1 and M_Y. INCY  must not be zero.
*
*  FDOT    (local input) pointer to a function of type VVDOT
*          On entry, FDOT points to a subroutine that computes the local
*          dot product of two vectors.
*
*  Further Details
*  ===============
*
*  When  the  result  of  a vector-oriented PBLAS call is a scalar, this
*  scalar  is set only within the process scope which owns the vector(s)
*  being operated on. Let sub( X ) be a generic term for the input  vec-
*  tor(s). Then, the processes owning the correct the answer is determi-
*  ned as follows:  if  an  operation involves more than one vector, the
*  processes receiving the result will be the union of the following set
*  of processes for each vector:
*
*  If N = 1, M_X = 1 and INCX = 1,  then  one cannot determine if a pro-
*  cess  row  or  process column owns the vector operand, therefore only
*  the process owning sub( X ) receives the correct result;
*
*  If  INCX = M_X, then sub( X )  is a vector distributed over a process
*  row. Each process in this row receives the result;
*
*  If  INCX = 1, then  sub( X )  is  a vector distributed over a process
*  column. Each process in this column receives the result;
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           * top;
   Int            RRorCC, Xcol, Xii, XisR, XisRow, Xjj, Xld, Xlinc, XmyprocD,
                  XmyprocR, XnprocsD, XnprocsR, XprocR, Xroc, Xrow, Ycol, Yii,
                  Yinb1D, YisR, YisRow, Yjj, Yld, Ylinc, YmyprocD, YmyprocR,
                  YnbD, YnpD, YnprocsD, YnprocsR, YprocD, YprocR, Yroc, Yrow,
                  ctxt, ione=1, k, kbb, kk, kn, ktmp, mycol, mydist, myproc,
                  myrow, npcol, nprow, p, size;
/*
*  .. Local Arrays ..
*/
   char           * Xptr = NULL, * Yptr = NULL, * buf = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCX[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Retrieve sub( X )'s local information: Xii, Xjj, Xrow, Xcol ...
*/
   PB_Cinfog2l( IX, JX, DESCX, nprow, npcol, myrow, mycol, &Xii, &Xjj,
                &Xrow, &Xcol );
   if( ( XisRow = ( INCX == DESCX[M_] ) ) != 0 )
   {
      Xld      = DESCX[LLD_];            Xlinc    = Xld;
      XmyprocD = mycol;                  XnprocsD = npcol;
      XprocR   = Xrow; XmyprocR = myrow; XnprocsR = nprow;
      XisR     = ( ( Xrow == -1 ) || ( XnprocsR == 1 ) );
   }
   else
   {
      Xld      = DESCX[LLD_];            Xlinc    = 1;
      XmyprocD = myrow;                  XnprocsD = nprow;
      XprocR   = Xcol; XmyprocR = mycol; XnprocsR = npcol;
      XisR     = ( ( Xcol == -1 ) || ( XnprocsR == 1 ) );
   }
/*
*  Retrieve sub( Y )'s local information: Yii, Yjj, Yrow, Ycol ...
*/
   PB_Cinfog2l( IY, JY, DESCY, nprow, npcol, myrow, mycol, &Yii, &Yjj,
                &Yrow, &Ycol );
   if( ( YisRow = ( INCY == DESCY[M_] ) ) != 0 )
   {
      YnbD     = DESCY[NB_]; Yld      = DESCY[LLD_]; Ylinc    = Yld;
      YprocR   = Yrow;       YmyprocR = myrow;       YnprocsR = nprow;
      YprocD   = Ycol;       YmyprocD = mycol;       YnprocsD = npcol;
      Yinb1D   = PB_Cfirstnb( N, JY, DESCY[INB_], YnbD );
   }
   else
   {
      YnbD     = DESCY[MB_]; Yld      = DESCY[LLD_]; Ylinc    = 1;
      YprocR   = Ycol;       YmyprocR = mycol;       YnprocsR = npcol;
      YprocD   = Yrow;       YmyprocD = myrow;       YnprocsD = nprow;
      Yinb1D   = PB_Cfirstnb( N, IY, DESCY[IMB_], YnbD );
   }
   YisR = ( ( YprocR == -1 ) || ( YnprocsR == 1 ) );
/*
*  Are sub( X ) and sub( Y ) both row or column vectors ?
*/
   RRorCC = ( ( XisRow && YisRow ) || ( !( XisRow ) && !( YisRow ) ) );
/*
*  sub( X ) is not distributed and sub( Y ) is distributed
*/
   if( !( XisR ) )
   {
/*
*  sub( X ) is not replicated. Since this operation is local if sub( X ) and
*  sub( Y ) are both row or column vectors, choose YprocR = XprocR when RRorCC,
*  and YprocR = 0 otherwise.
*/
      if( YisR ) { YprocR = ( ( RRorCC ) ? XprocR : 0 ); }
/*
*  Now, it is just like sub( Y ) is not replicated, this information however is
*  kept in YisR for later use.
*/
      if( RRorCC )
      {
/*
*  sub( X ) and sub( Y ) are both row or column vectors
*/
         if( XprocR == YprocR )
         {
/*
*  sub( X ) and sub( Y ) are in the same process row or column
*/
            if( ( XmyprocR == XprocR ) || ( YmyprocR == YprocR ) )
            {
               size = TYPE->size;
               YnpD = PB_Cnumroc( N, 0, Yinb1D, YnbD, YmyprocD, YprocD,
                                  YnprocsD );
/*
*  In a given process, the dot product is computed with sub( Y ) and the cor-
*  responding non distributed part of sub( X ). In the other processes, this
*  part of sub( X ) is simply ignored.
*/
               if( YnpD > 0 )
               {
                  Yroc = YprocD;
                  if( XisRow ) { kk = Yjj; ktmp = JX + N; kn = JX + Yinb1D; }
                  else         { kk = Yii; ktmp = IX + N; kn = IX + Yinb1D; }

                  if( YmyprocD == Yroc )
                  {
                     FDOT( &Yinb1D, DOT, Mptr( X, Xii, Xjj, Xld, size ), &Xlinc,
                           Mptr( Y, Yii, Yjj, Yld, size ), &Ylinc );
                     kk += Yinb1D;
                  }
                  Yroc = MModAdd1( Yroc, YnprocsD );

                  for( k = kn; k < ktmp; k += YnbD )
                  {
                     kbb = ktmp - k; kbb = MIN( kbb, YnbD );
                     if( YmyprocD == Yroc )
                     {
                        if( XisRow )
                           FDOT( &kbb, DOT, Mptr( X, Xii, k, Xld, size ),
                                 &Xlinc, Mptr( Y, Yii, kk, Yld, size ),
                                 &Ylinc );
                        else
                           FDOT( &kbb, DOT, Mptr( X, k, Xjj, Xld, size ),
                                 &Xlinc, Mptr( Y, kk, Yjj, Yld, size ),
                                 &Ylinc );
                        kk += kbb;
                     }
                     Yroc = MModAdd1( Yroc, YnprocsD );
                  }
               }
/*
*  Replicate locally scattered dot product by reducing it
*/
               if( XisRow )
               {
                  top = PB_Ctop( &ctxt, COMBINE, ROW, TOP_GET );
                  TYPE->Cgsum2d( ctxt, ROW,    top, 1, 1, DOT, 1, -1, 0 );
               }
               else
               {
                  top = PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
                  TYPE->Cgsum2d( ctxt, COLUMN, top, 1, 1, DOT, 1, -1, 0 );
               }
            }
         }
         else
         {
/*
*  sub( X ) and sub( Y ) are in a different process row or column
*/
            if( YmyprocR == YprocR )
            {
               size = TYPE->size;
               YnpD = PB_Cnumroc( N, 0, Yinb1D, YnbD, YmyprocD, YprocD,
                                  YnprocsD );
/*
*  If I own a piece of sub( Y ), then send it to the process row or column where
*  sub( X ) resides and receive the dot product when sub( Y ) is not replicated.
*/
               if( YisRow )
               {
                  if( YnpD > 0 )
                     TYPE->Cgesd2d( ctxt, 1, YnpD, Mptr( Y, Yii, Yjj, Yld,
                                    size ), Yld, XprocR, YmyprocD );
                  TYPE->Cgerv2d( ctxt, 1, 1, DOT, 1, XprocR, XmyprocD );
               }
               else
               {
                  if( YnpD > 0 )
                     TYPE->Cgesd2d( ctxt, YnpD, 1, Mptr( Y, Yii, Yjj, Yld,
                                    size ), Yld, YmyprocD, XprocR );
                  TYPE->Cgerv2d( ctxt, 1, 1, DOT, 1, XmyprocD, XprocR );
               }
            }

            if( XmyprocR == XprocR )
            {
               size = TYPE->size;
               YnpD = PB_Cnumroc( N, 0, Yinb1D, YnbD, YmyprocD, YprocD,
                                  YnprocsD );
/*
*  If I own sub( X ), then receive the distributed part of sub( Y ) owned by
*  the process where sub( Y ) resides in my row or column. Compute the partial
*  dot product as if sub( Y ) would reside in the same process row or column as
*  sub( X ). Combine the local results.
*/
               if( YnpD > 0 )
               {
                  buf = PB_Cmalloc( YnpD * size );
                  if( YisRow )
                     TYPE->Cgerv2d( ctxt, 1, YnpD, buf,    1, YprocR,
                                    XmyprocD );
                  else
                     TYPE->Cgerv2d( ctxt, YnpD, 1, buf, YnpD, XmyprocD,
                                    YprocR );

                  Yroc = YprocD;
                  kk   = 0;
                  if( XisRow ) { ktmp = JX + N; kn = JX + Yinb1D; }
                  else         { ktmp = IX + N; kn = IX + Yinb1D; }

                  if( YmyprocD == Yroc )
                  {
                     FDOT( &Yinb1D, DOT, Mptr( X, Xii, Xjj, Xld, size ),
                           &Xlinc, buf, &ione );
                     kk += Yinb1D;
                  }
                  Yroc = MModAdd1( Yroc, YnprocsD );

                  for( k = kn; k < ktmp; k += YnbD )
                  {
                     kbb = ktmp - k; kbb = MIN( kbb, YnbD );
                     if( YmyprocD == Yroc )
                     {
                        if( XisRow )
                           FDOT( &kbb, DOT, Mptr( X, Xii, k, Xld, size ),
                                 &Xlinc, buf+kk*size, &ione );
                        else
                           FDOT( &kbb, DOT, Mptr( X, k, Xjj, Xld, size ),
                                 &Xlinc, buf+kk*size, &ione );
                        kk += kbb;
                     }
                     Yroc = MModAdd1( Yroc, YnprocsD );
                  }
                  if( buf ) free( buf );
               }
/*
*  Combine the local results within the process row or column XprocR and
*  send the result to the process row or column YprocR when sub( Y ) is not
*  replicated.
*/
               if( XisRow )
               {
                  top = PB_Ctop( &ctxt, COMBINE, ROW,    TOP_GET );
                  TYPE->Cgsum2d( ctxt, ROW,    top, 1, 1, DOT, 1, -1, 0 );
                  if( !YisR )
                     TYPE->Cgesd2d( ctxt, 1, 1, DOT, 1, YprocR, YmyprocD );
               }
               else
               {
                  top = PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
                  TYPE->Cgsum2d( ctxt, COLUMN, top, 1, 1, DOT, 1, -1, 0 );
                  if( !YisR )
                     TYPE->Cgesd2d( ctxt, 1, 1, DOT, 1, YmyprocD, YprocR );
               }
            }
         }

         if( YisR )
         {
/*
*  If sub( Y ) is replicated, then bcast the result
*/
            if( XisRow )
            {
               top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
               if( XmyprocR == XprocR )
                  TYPE->Cgebs2d( ctxt, COLUMN, top, 1, 1, DOT, 1 );
               else
                  TYPE->Cgebr2d( ctxt, COLUMN, top, 1, 1, DOT, 1, XprocR,
                                 XmyprocD );
            }
            else
            {
               top = PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
               if( XmyprocR == XprocR )
                  TYPE->Cgebs2d( ctxt, ROW,    top, 1, 1, DOT, 1 );
               else
                  TYPE->Cgebr2d( ctxt, ROW,    top, 1, 1, DOT, 1, XmyprocD,
                                 XprocR );
            }
         }
      }
      else
      {
/*
*  sub( X ) and sub( Y ) are not both row or column vectors
*/
         if( ( XmyprocR == XprocR ) || ( YmyprocR == YprocR ) )
         {
            size = TYPE->size;
            Xroc = 0;
            if( XisRow ) { ktmp = JX + N; kn = JX + Yinb1D; }
            else         { ktmp = IX + N; kn = IX + Yinb1D; }
/*
*  Loop over the processes in which sub( Y ) resides, for each process find the
*  next process Xroc and compute the dot product. After this, it will be needed
*  to reduce the local dot produsts as above.
*/
            for( p = 0; p < YnprocsD; p++ )
            {
               mydist = MModSub( p, YprocD, YnprocsD );
               myproc = MModAdd( YprocD, mydist, YnprocsD );

               if( ( XprocR == p ) && ( YprocR == Xroc ) )
               {
/*
*  Compute locally the partial dot product at the intersection of the process
*  cross.
*/
                  if( ( XmyprocR == p ) && ( XmyprocD == Xroc ) )
                  {
                     YnpD = PB_Cnumroc( N, 0, Yinb1D, YnbD, p, YprocD,
                                        YnprocsD );
                     if( YnpD > 0 )
                     {
                        Yroc = YprocD;
                        kk   = ( XisRow ? Yii : Yjj );

                        if( myproc == Yroc )
                        {
                           FDOT( &Yinb1D, DOT, Mptr( X, Xii, Xjj, Xld, size ),
                                 &Xlinc, Mptr( Y, Yii, Yjj, Yld, size ),
                                 &Ylinc );
                           kk += Yinb1D;
                        }
                        Yroc = MModAdd1( Yroc, YnprocsD );

                        for( k = kn; k < ktmp; k += YnbD )
                        {
                           kbb = ktmp - k; kbb = MIN( kbb, YnbD );
                           if( myproc == Yroc )
                           {
                              if( XisRow )
                                 FDOT( &kbb, DOT, Mptr( X, Xii, k, Xld, size ),
                                       &Xlinc, Mptr( Y, kk, Yjj, Yld, size ),
                                       &Ylinc );
                              else
                                 FDOT( &kbb, DOT, Mptr( X, k, Xjj, Xld, size ),
                                       &Xlinc, Mptr( Y, Yii, kk, Yld, size ),
                                       &Ylinc );
                              kk += kbb;
                           }
                           Yroc = MModAdd1( Yroc, YnprocsD );
                        }
                     }
                  }
               }
               else
               {
/*
*  Message exchange
*/
                  if( ( YmyprocR == YprocR ) && ( YmyprocD == p ) )
                  {
                     YnpD = PB_Cnumroc( N, 0, Yinb1D, YnbD, p, YprocD,
                                        YnprocsD );
                     if( YnpD > 0 )
                     {
                        if( XisRow )
                           TYPE->Cgesd2d( ctxt, YnpD, 1, Mptr( Y, Yii, Yjj, Yld,
                                          size ), Yld, XprocR, Xroc );
                        else
                           TYPE->Cgesd2d( ctxt, 1, YnpD, Mptr( Y, Yii, Yjj, Yld,
                                          size ), Yld, Xroc, XprocR );
                     }
                  }

                  if( ( XmyprocR == XprocR ) && ( XmyprocD == Xroc ) )
                  {
                     YnpD = PB_Cnumroc( N, 0, Yinb1D, YnbD, p, YprocD,
                                        YnprocsD );
                     if( YnpD > 0 )
                     {
                        buf  = PB_Cmalloc( YnpD * size );
                        Yroc = YprocD;
                        kk   = 0;
/*
*  Receive the piece of sub( Y ) that I should handle
*/
                        if( XisRow )
                           TYPE->Cgerv2d( ctxt, YnpD, 1, buf, YnpD, p, YprocR );
                        else
                           TYPE->Cgerv2d( ctxt, 1, YnpD, buf,    1, YprocR, p );

                        if( myproc == Yroc )
                        {
                           FDOT( &Yinb1D, DOT, Mptr( X, Xii, Xjj, Xld, size ),
                                 &Xlinc, buf, &ione );
                           kk += Yinb1D;
                        }
                        Yroc = MModAdd1( Yroc, YnprocsD );

                        for( k = kn; k < ktmp; k += YnbD )
                        {
                           kbb = ktmp - k; kbb = MIN( kbb, YnbD );
                           if( myproc == Yroc )
                           {
                              if( XisRow )
                                 FDOT( &kbb, DOT, Mptr( X, Xii, k, Xld, size ),
                                       &Xlinc, buf+kk*size, &ione );
                              else
                                 FDOT( &kbb, DOT, Mptr( X, k, Xjj, Xld, size ),
                                       &Xlinc, buf+kk*size, &ione );
                              kk += kbb;
                           }
                           Yroc = MModAdd1( Yroc, YnprocsD );
                        }
                        if( buf ) free( buf );
                     }
                  }
               }
               Xroc = MModAdd1( Xroc, XnprocsD );
            }
/*
*  Combine the local results in sub( X )'s scope
*/
            if( XmyprocR == XprocR )
            {
               if( XisRow )
               {
                  top = PB_Ctop( &ctxt, COMBINE, ROW,    TOP_GET );
                  TYPE->Cgsum2d( ctxt, ROW,    top, 1, 1, DOT, 1, -1, 0 );
               }
               else
               {
                  top = PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
                  TYPE->Cgsum2d( ctxt, COLUMN, top, 1, 1, DOT, 1, -1, 0 );
               }
            }
         }
/*
*  Broadcast the result in sub( Y )'s scope
*/
         if( YisR || ( YmyprocR == YprocR ) )
         {
            if( YisRow )
            {
               top = PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
               if( XmyprocR == XprocR )
                  TYPE->Cgebs2d( ctxt, ROW,    top, 1, 1, DOT, 1 );
               else
                  TYPE->Cgebr2d( ctxt, ROW,    top, 1, 1, DOT, 1, YmyprocR,
                                 XprocR );
            }
            else
            {
               top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
               if( XmyprocR == XprocR )
                  TYPE->Cgebs2d( ctxt, COLUMN, top, 1, 1, DOT, 1 );
               else
                  TYPE->Cgebr2d( ctxt, COLUMN, top, 1, 1, DOT, 1, XprocR,
                                 YmyprocR );
            }
         }
      }
   }
   else
   {
/*
*  sub( X ) is replicated in every process. Compute the local dot product in
*  process row or column YprocR when sub( Y ) is not replicated and in every
*  process otherwise.
*/
      if( YisR || ( YmyprocR == YprocR ) )
      {
         size = TYPE->size;
         Yroc = YprocD;
         kk   = ( YisRow ? Yjj : Yii );

         if( XisRow ) { ktmp = JX + N; kn = JX + Yinb1D; }
         else         { ktmp = IX + N; kn = IX + Yinb1D; }

         if( YmyprocD == Yroc )
         {
            FDOT( &Yinb1D, DOT, Mptr( X, Xii, Xjj, Xld, size ), &Xlinc, Mptr( Y,
                  Yii, Yjj, Yld, size ), &Ylinc );
            kk += Yinb1D;
         }
         Yroc = MModAdd1( Yroc, YnprocsD );

         for( k = kn; k < ktmp; k += YnbD )
         {
            kbb = ktmp - k; kbb = MIN( kbb, YnbD );
            if( YmyprocD == Yroc )
            {
               if( XisRow ) { Xptr = Mptr( X, Xii,   k, Xld, size ); }
               else         { Xptr = Mptr( X,   k, Xjj, Xld, size ); }
               if( YisRow ) { Yptr = Mptr( Y, Yii,  kk, Yld, size ); }
               else         { Yptr = Mptr( Y,  kk, Yjj, Yld, size ); }
               FDOT( &kbb, DOT, Xptr, &Xlinc, Yptr, &Ylinc );
               kk += kbb;
            }
            Yroc = MModAdd1( Yroc, YnprocsD );
         }
      }

     if( YisR )
     {
/*
*  sub( Y ) is replicated, combine the results in each process row or column.
*/
         if( YisRow )
         {
            top = PB_Ctop( &ctxt, COMBINE, ROW,    TOP_GET );
            TYPE->Cgsum2d( ctxt, ROW,    top, 1, 1, DOT, 1, -1, 0 );
         }
         else
         {
            top = PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
            TYPE->Cgsum2d( ctxt, COLUMN, top, 1, 1, DOT, 1, -1, 0 );
         }
      }
      else
      {
/*
*  sub( Y ) is not replicated, combine the results in the entire grid at once.
*/
         top = PB_Ctop( &ctxt, COMBINE, ALL, TOP_GET );
         TYPE->Cgsum2d( ctxt, ALL, top, 1, 1, DOT, 1, -1, 0 );
      }
   }
/*
*  End of PB_CpdotND
*/
}
