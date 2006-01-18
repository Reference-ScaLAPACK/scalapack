/* ---------------------------------------------------------------------
*
*  -- PBLAS routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*  ---------------------------------------------------------------------
*/
/*
*  Include files
*/
#include "pblas.h"
#include "PBpblas.h"
#include "PBtools.h"
#include "PBblacs.h"
#include "PBblas.h"

#ifdef __STDC__
void psswap_( int * N,
              float * X, int * IX, int * JX, int * DESCX, int * INCX,
              float * Y, int * IY, int * JY, int * DESCY, int * INCY )
#else
void psswap_( N, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
/*
*  .. Scalar Arguments ..
*/
   int            * INCX, * INCY, * IX, * IY, * JX, * JY, * N;
/*
*  .. Array Arguments ..
*/
   int            * DESCX, * DESCY;
   float          * X, * Y;
#endif
{
/*
*  Purpose
*  =======
*
*  PSSWAP  swaps two subvectors,
*
*     sub( Y ) := sub( X ) and sub( X ) := sub( Y )
*
*  where sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
*                         X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X,
*
*        sub( Y ) denotes Y(IY,JY:JY+N-1) if INCY = M_Y,
*                         Y(IY:IY+N-1,JY) if INCY = 1 and INCY <> M_Y.
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
*  N       (global input) INTEGER
*          On entry,  N  specifies the  length of the  subvectors to  be
*          swapped. N must be at least zero.
*
*  X       (local input/local output) REAL array
*          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
*          is   at  least  MAX( 1, Lr( 1, IX ) )  when  INCX = M_X   and
*          MAX( 1, Lr( 1, IX+N-1 ) )  otherwise,  and,  Kx  is  at least
*          Lc( 1, JX+N-1 )  when  INCX = M_X  and Lc( 1, JX ) otherwise.
*          Before  entry,  this array  contains the local entries of the
*          matrix X. On exit, sub( X ) is overwritten with sub( Y ).
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
*  Y       (local input/local output) REAL array
*          On entry, Y is an array of dimension (LLD_Y, Ky), where LLD_Y
*          is   at  least  MAX( 1, Lr( 1, IY ) )  when  INCY = M_Y   and
*          MAX( 1, Lr( 1, IY+N-1 ) )  otherwise,  and,  Ky  is  at least
*          Lc( 1, JY+N-1 )  when  INCY = M_Y  and Lc( 1, JY ) otherwise.
*          Before  entry,  this array  contains the local entries of the
*          matrix Y. On exit, sub( Y ) is overwritten with sub( X ).
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
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           Xscope, Yscope, * one, * top, tran, * zero;
   int            OneBlock, OneDgrid, RRorCC, Square, Xcol, Xi, XisD, XisR,
                  Xinb1D, XinbD, XisRow, Xii, Xj, Xjj, Xld, Xlinc, Xm, XmyprocD,
                  XmyprocR, Xn, XnbD, XnpD, XnprocsD, XnprocsR, XprocD, XprocR,
                  Xroc, Xrow, Ycol, Yi, Yii, Yinb1D, YinbD, YisD, YisR, YisRow,
                  Yj, Yjj, Yld, Ylinc, Ym, YmyprocD, YmyprocR, Yn, YnbD, YnpD,
                  YnprocsD, YnprocsR, YprocD, YprocR, Yroc, Yrow, cdst, csrc,
                  ctxt, dst, gcdPQ, info, ione=1, k, l, lcmPQ, lcmb, mycol,
                  myrow, npcol, npq, nprow, p, q, rdst, rsrc, src, size;
   PBTYP_T        * type;
   PB_VM_T        VM;
/*
*  .. Local Arrays ..
*/
   int            Xd[DLEN_], Yd[DLEN_];
   char           * buf = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   PB_CargFtoC( *IX, *JX, DESCX, &Xi, &Xj, Xd );
   PB_CargFtoC( *IY, *JY, DESCY, &Yi, &Yj, Yd );
#ifndef NO_ARGCHK
/*
*  Test the input parameters
*/
   Cblacs_gridinfo( ( ctxt = Xd[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
   if( !( info = ( ( nprow == -1 ) ? -( 501 + CTXT_ ) : 0 ) ) )
   {
      PB_Cchkvec( ctxt, "PSSWAP", "X", *N, 1, Xi, Xj, Xd, *INCX,  5, &info );
      PB_Cchkvec( ctxt, "PSSWAP", "Y", *N, 1, Yi, Yj, Yd, *INCY, 10, &info );
   }
   if( info ) { PB_Cabort( ctxt, "PSSWAP", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( *N == 0 ) return;
/*
*  Retrieve process grid information
*/
#ifdef NO_ARGCHK
   Cblacs_gridinfo( ( ctxt = Xd[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
#endif
/*
*  Determine if sub( X ) is distributed or not
*/
   if( ( XisRow = ( *INCX == Xd[M_] ) ) != 0 )
      XisD = ( ( Xd[CSRC_] >= 0 ) && ( ( XnprocsD = npcol ) > 1 ) );
   else
      XisD = ( ( Xd[RSRC_] >= 0 ) && ( ( XnprocsD = nprow ) > 1 ) );
/*
*  Determine if sub( Y ) is distributed or not
*/
   if( ( YisRow = ( *INCY == Yd[M_] ) ) != 0 )
      YisD = ( ( Yd[CSRC_] >= 0 ) && ( ( YnprocsD = npcol ) > 1 ) );
   else
      YisD = ( ( Yd[RSRC_] >= 0 ) && ( ( YnprocsD = nprow ) > 1 ) );
/*
*  Are sub( X ) and sub( Y ) both row or column vectors ?
*/
   RRorCC = ( ( XisRow && YisRow ) || ( !( XisRow ) && !( YisRow ) ) );
/*
*  XisD && YisD <=> both vector operands are indeed distributed
*/
   if( XisD && YisD )
   {
/*
*  Retrieve sub( X )'s local information: Xii, Xjj, Xrow, Xcol ...
*/
      PB_Cinfog2l( Xi, Xj, Xd, nprow, npcol, myrow, mycol, &Xii, &Xjj,
                   &Xrow, &Xcol );
      if( XisRow )
      {
         XinbD  = Xd[INB_ ]; XnbD     = Xd[NB_ ]; Xld      = Xd[LLD_];
         Xlinc  = Xld;
         XprocD = Xcol;      XmyprocD = mycol;
         XprocR = Xrow;      XmyprocR = myrow;    XnprocsR = nprow;
         XisR   = ( ( Xrow == -1 ) || ( XnprocsR == 1 ) );
         Mfirstnb( Xinb1D, *N, Xj, XinbD, XnbD );
      }
      else
      {
         XinbD  = Xd[IMB_ ]; XnbD     = Xd[MB_ ]; Xld      = Xd[LLD_];
         Xlinc  = 1;
         XprocD = Xrow;      XmyprocD = myrow;
         XprocR = Xcol;      XmyprocR = mycol;    XnprocsR = npcol;
         XisR   = ( ( Xcol == -1 ) || ( XnprocsR == 1 ) );
         Mfirstnb( Xinb1D, *N, Xi, XinbD, XnbD );
      }
/*
*  Retrieve sub( Y )'s local information: Yii, Yjj, Yrow, Ycol ...
*/
      PB_Cinfog2l( Yi, Yj, Yd, nprow, npcol, myrow, mycol, &Yii, &Yjj,
                   &Yrow, &Ycol );
      if( YisRow )
      {
         YinbD  = Yd[INB_ ]; YnbD     = Yd[NB_ ]; Yld      = Yd[LLD_];
         Ylinc  = Yld;
         YprocD = Ycol;      YmyprocD = mycol;
         YprocR = Yrow;      YmyprocR = myrow;    YnprocsR = nprow;
         YisR   = ( ( Yrow == -1 ) || ( YnprocsR == 1 ) );
         Mfirstnb( Yinb1D, *N, Yj, YinbD, YnbD );
      }
      else
      {
         YinbD  = Yd[IMB_ ]; YnbD     = Yd[MB_ ]; Yld      = Yd[LLD_];
         Ylinc  = 1;
         YprocD = Yrow;      YmyprocD = myrow;
         YprocR = Ycol;      YmyprocR = mycol;    YnprocsR = npcol;
         YisR   = ( ( Ycol == -1 ) || ( YnprocsR == 1 ) );
         Mfirstnb( Yinb1D, *N, Yi, YinbD, YnbD );
      }
/*
*  Do sub( X ) and sub( Y ) span more than one process ?
*/
      OneDgrid = ( ( XnprocsD ==  1 ) && ( YnprocsD ==  1 ) );
      OneBlock = ( ( Xinb1D   >= *N ) && ( Yinb1D   >= *N ) );
/*
*  Are sub( X ) and sub( Y ) distributed in the same manner ?
*/
      Square   = ( ( Xinb1D   ==   Yinb1D ) && ( XnbD == YnbD ) &&
                   ( XnprocsD == YnprocsD ) );

      if( !( XisR ) )
      {
/*
*  sub( X ) is distributed but not replicated
*/
         if( YisR )
         {
/*
*  If sub( X ) is not replicated, but sub( Y ) is, a process row or column
*  YprocR need to be selected. It will contain the non-replicated vector to
*  swap sub( X ) with.
*/
            if( RRorCC )
            {
/*
*  sub( X ) and sub( Y ) are both row or column vectors
*/
               if( ( OneDgrid || OneBlock || Square ) && ( XprocD == YprocD ) )
               {
/*
*  sub( X ) and sub( Y ) start in the same process row or column XprocD=YprocD.
*  Enforce a purely local operation by choosing YprocR to be equal to XprocR.
*/
                  YprocR = XprocR;
               }
               else
               {
/*
*  Otherwise, communication has to occur, so choose the next process row or
*  column for YprocR to maximize the number of links, i.e reduce contention.
*/
                  YprocR = MModAdd1( XprocR, XnprocsR );
               }
            }
            else
            {
/*
*  sub( X ) and sub( Y ) are distributed in orthogonal directions, what is
*  chosen for YprocR does not really matter. Select the process origin.
*/
               YprocR = XprocD;
            }
         }
         else
         {
/*
*  Neither sub( X ) nor sub( Y ) are replicated. If I am not in process row or
*  column XprocR and not in process row or column YprocR, then quick return.
*/
            if( ( XmyprocR != XprocR ) && ( YmyprocR != YprocR ) )
               return;
         }
      }
      else
      {
/*
*  sub( X ) is distributed and replicated (so no quick return possible)
*/
         if( YisR )
         {
/*
*  sub( Y ) is distributed and replicated as well
*/
            if( RRorCC )
            {
/*
*  sub( X ) and sub( Y ) are both row or column vectors
*/
               if( ( OneDgrid || OneBlock || Square ) && ( XprocD == YprocD ) )
               {
/*
*  sub( X ) and sub( Y ) start in the same process row or column XprocD=YprocD.
*  Enforce a purely local operation by choosing XprocR and YprocR to be equal
*  to zero.
*/
                  XprocR = YprocR = 0;
               }
               else
               {
/*
*  Otherwise, communication has to occur, so select YprocR to be zero and the
*  next process row or column for XprocR in order to maximize the number of
*  used links, i.e reduce contention.
*/
                  YprocR = 0;
                  XprocR = MModAdd1( YprocR, YnprocsR );
               }
            }
            else
            {
/*
*  sub( X ) and sub( Y ) are distributed in orthogonal directions, select the
*  origin processes.
*/
               XprocR = YprocD;
               YprocR = XprocD;
            }
         }
         else
         {
/*
*  sub( Y ) is distributed, but not replicated
*/
            if( RRorCC )
            {
/*
*  sub( X ) and sub( Y ) are both row or column vectors
*/
               if( ( OneDgrid || OneBlock || Square ) && ( XprocD == YprocD ) )
               {
/*
*  sub( X ) and sub( Y ) start in the same process row or column XprocD=YprocD.
*  Enforce a purely local operation by choosing XprocR to be equal to YprocR.
*/
                  XprocR = YprocR;
               }
               else
               {
/*
*  Otherwise, communication has to occur, so choose the next process row or
*  column for XprocR to maximize the number of links, i.e reduce contention.
*/
                  XprocR = MModAdd1( YprocR, YnprocsR );
               }
            }
            else
            {
/*
*  sub( X ) and sub( Y ) are distributed in orthogonal directions, what is
*  chosen for XprocR does not really matter. Select the origin process.
*/
               XprocR = YprocD;
            }
         }
      }
/*
*  Even if sub( X ) and/or sub( Y ) are replicated, only two process row or
*  column are active, namely XprocR and YprocR. If any of those operands is
*  replicated, broadcast will occur (unless there is an easy way out).
*/
      type = PB_Cstypeset(); size = type->size;
/*
*  A purely local operation occurs iff the operands start in the same process
*  and, if either the grid is mono-dimensional or there is a single local block
*  to be swapped or if both operands are aligned.
*/
      if( ( (    RRorCC   && ( XprocD == YprocD ) && ( XprocR == YprocR ) ) ||
            ( !( RRorCC ) && ( XprocD == YprocR ) && ( XprocR == YprocD ) ) ) &&
          ( OneDgrid || OneBlock || ( RRorCC && Square ) ) )
      {
         if( ( !XisR &&         ( XmyprocR == XprocR ) &&
               !YisR &&         ( YmyprocR == YprocR ) ) ||
             ( !XisR && YisR && ( YmyprocR == YprocR ) ) ||
             ( !YisR && XisR && ( XmyprocR == XprocR ) ) ||
             (  XisR && YisR                           ) )
         {
            XnpD = PB_Cnumroc( *N, 0, Xinb1D, XnbD, XmyprocD, XprocD,
                               XnprocsD );
            YnpD = PB_Cnumroc( *N, 0, Yinb1D, YnbD, YmyprocD, YprocD,
                               YnprocsD );
            if( ( XnpD > 0 ) && ( YnpD > 0 ) )
            {
               sswap_( &XnpD,
                       Mptr( ((char *) X), Xii, Xjj, Xld, size ), &Xlinc,
                       Mptr( ((char *) Y), Yii, Yjj, Yld, size ), &Ylinc );
            }
            if( RRorCC && XisR && YisR ) return;
         }
      }
      else if( ( RRorCC && OneDgrid ) || OneBlock || Square )
      {
/*
*  Otherwise, it may be possible to swap the distributed vectors in a single
*  message exchange iff the grid is mono-dimensional and the operands are
*  distributed in the same direction, or there is just one block to be exchanged
*  or if both operands are similarly distributed in their respective direction.
*/
         if( RRorCC && ( XprocR != YprocR ) )
         {
/*
*  Both operands are distributed in the same direction, but reside in different
*  process rows or columns.
*/
            if( XmyprocR == XprocR )
            {
               XnpD = PB_Cnumroc( *N, 0, Xinb1D, XnbD, XmyprocD, XprocD,
                                  XnprocsD );
               if( XnpD > 0 )
               {
                  dst = YprocD + MModSub( XmyprocD, XprocD, XnprocsD );
                  dst = MPosMod( dst, YnprocsD );
                  if( XisRow )
                  {
                     Csgesd2d( ctxt, 1, XnpD, Mptr( ((char*) X), Xii, Xjj,
                               Xld, size ), Xld, YprocR, dst );
                     Csgerv2d( ctxt, 1, XnpD, Mptr( ((char*) X), Xii, Xjj,
                               Xld, size ), Xld, YprocR, dst );
                  }
                  else
                  {

                     Csgesd2d( ctxt, XnpD, 1, Mptr( ((char*) X), Xii, Xjj,
                               Xld, size ), Xld, dst, YprocR );
                     Csgerv2d( ctxt, XnpD, 1, Mptr( ((char*) X), Xii, Xjj,
                               Xld, size ), Xld, dst, YprocR );
                  }
               }
            }
            if( YmyprocR == YprocR )
            {
               YnpD = PB_Cnumroc( *N, 0, Yinb1D, YnbD, YmyprocD, YprocD,
                                  YnprocsD );
               if( YnpD > 0 )
               {
                  dst = XprocD + MModSub( YmyprocD, YprocD, YnprocsD );
                  dst = MPosMod( dst, XnprocsD );
                  if( YisRow )
                  {
                     Csgesd2d( ctxt, 1, YnpD, Mptr( ((char*) Y), Yii, Yjj,
                               Yld, size ), Yld, XprocR, dst );
                     Csgerv2d( ctxt, 1, YnpD, Mptr( ((char*) Y), Yii, Yjj,
                               Yld, size ), Yld, XprocR, dst );
                  }
                  else
                  {
                     Csgesd2d( ctxt, YnpD, 1, Mptr( ((char*) Y), Yii, Yjj,
                               Yld, size ), Yld, dst, XprocR );
                     Csgerv2d( ctxt, YnpD, 1, Mptr( ((char*) Y), Yii, Yjj,
                               Yld, size ), Yld, dst, XprocR );
                  }
               }
            }
         }
         else
         {
/*
*  General case when just one message needs to be exchanged
*/
            if( XmyprocR == XprocR )
            {
/*
*  The processes owning a piece of sub( X ) send it to the corresponding
*  process owning s piece of sub ( Y ).
*/
               XnpD = PB_Cnumroc( *N, 0, Xinb1D, XnbD, XmyprocD, XprocD,
                                  XnprocsD );
               if( XnpD > 0 )
               {
                  dst = YprocD + MModSub( XmyprocD, XprocD, XnprocsD );
                  dst = MPosMod( dst, YnprocsD );
                  if( YisRow ) { rdst = YprocR; cdst = dst; }
                  else         { rdst = dst; cdst = YprocR; }

                  if( ( myrow == rdst ) && ( mycol == cdst ) )
                  {
                     sswap_( &XnpD,  Mptr( ((char *) X), Xii, Xjj, Xld,
                             size ), &Xlinc, Mptr( ((char *) Y), Yii, Yjj, Yld,
                             size ), &Ylinc );
                  }
                  else
                  {
                     if( XisRow )
                        Csgesd2d( ctxt, 1, XnpD, Mptr( ((char *) X), Xii,
                                  Xjj, Xld, size ), Xld, rdst, cdst );
                     else
                        Csgesd2d( ctxt, XnpD, 1, Mptr( ((char *) X), Xii,
                                  Xjj, Xld, size ), Xld, rdst, cdst );
                  }
               }
            }
            if( YmyprocR == YprocR )
            {
/*
*  The processes owning a piece of sub( Y ) receive the corresponding piece
*  of sub( X ) and send the piece of sub( Y ) they own to the same process.
*/
               YnpD = PB_Cnumroc( *N, 0, Yinb1D, YnbD, YmyprocD, YprocD,
                                  YnprocsD );
               if( YnpD > 0 )
               {
                  src = XprocD + MModSub( YmyprocD, YprocD, YnprocsD );
                  src = MPosMod( src, XnprocsD );
                  if( XisRow ) { rsrc = XprocR; csrc = src; }
                  else         { rsrc = src; csrc = XprocR; }

                  if( ( myrow != rsrc ) || ( mycol != csrc ) )
                  {
                     buf = PB_Cmalloc( YnpD * size );
                     if( XisRow )
                        Csgerv2d( ctxt, 1, YnpD, buf,    1, rsrc, csrc );
                     else
                        Csgerv2d( ctxt, YnpD, 1, buf, YnpD, rsrc, csrc );
                     if( YisRow )
                        Csgesd2d( ctxt, 1, YnpD, Mptr( ((char *) Y), Yii,
                                  Yjj, Yld, size ), Yld, rsrc, csrc );
                     else
                        Csgesd2d( ctxt, YnpD, 1, Mptr( ((char *) Y), Yii,
                                  Yjj, Yld, size ), Yld, rsrc, csrc );
                     scopy_( &YnpD, buf, &ione, Mptr( ((char *) Y), Yii,
                             Yjj, Yld, size ), &Ylinc );
                     if( buf ) free( buf );
                  }
               }
            }
            if( XmyprocR == XprocR )
            {
/*
*  The processes owning a piece of sub( X ) receive the corresponding piece
*  of sub( Y ).
*/
               if( XnpD > 0 )
               {
                  if( ( myrow != rdst ) || ( mycol != cdst ) )
                  {
                     buf = PB_Cmalloc( XnpD * size );
                     if( YisRow )
                        Csgerv2d( ctxt, 1, XnpD, buf,    1, rdst, cdst );
                     else
                        Csgerv2d( ctxt, XnpD, 1, buf, XnpD, rdst, cdst );
                     scopy_( &XnpD, buf, &ione, Mptr( ((char *) X), Xii,
                             Xjj, Xld, size ), &Xlinc );
                     if( buf ) free( buf );
                  }
               }
            }
         }
      }
      else if( ( XmyprocR == XprocR ) || ( YmyprocR == YprocR ) )
      {
/*
*  General case
*/
         tran   = ( RRorCC ? CNOTRAN : CTRAN );
         if( XisRow ) { Xscope = CCOLUMN; Xm = 1; rsrc = XprocR; }
         else         { Xscope = CROW;    Xn = 1; csrc = XprocR; }
         if( YisRow ) { Yscope = CCOLUMN; Ym = 1; rdst = YprocR; }
         else         { Yscope = CROW;    Yn = 1; cdst = YprocR; }
         lcmb   = PB_Clcm( XnprocsD * XnbD, YnprocsD * YnbD );
         one    = type->one; zero = type->zero;
         gcdPQ  = PB_Cgcd( XnprocsD, YnprocsD );
         lcmPQ  = ( XnprocsD / gcdPQ ) * YnprocsD;

         for( k = 0; k < gcdPQ; k++ )
         {
            p = 0; q = k;

            for( l = 0; l < lcmPQ; l++ )
            {
               Xroc = MModAdd( XprocD, p, XnprocsD );
               Yroc = MModAdd( YprocD, q, YnprocsD );

               if( ( XmyprocD == Xroc ) || ( YmyprocD == Yroc ) )
               {
                  XnpD = PB_Cnumroc( *N, 0, Xinb1D, XnbD, Xroc, XprocD,
                                     XnprocsD );
                  YnpD = PB_Cnumroc( *N, 0, Yinb1D, YnbD, Yroc, YprocD,
                                     YnprocsD );
                  PB_CVMinit( &VM, 0, XnpD, YnpD, Xinb1D, Yinb1D, XnbD, YnbD,
                              p, q, XnprocsD, YnprocsD, lcmb );
                  if( npq = PB_CVMnpq( &VM ) )
                  {
                     if( (    RRorCC   && ( Xroc ==   Yroc ) &&
                           ( XprocR == YprocR ) ) ||
                         ( !( RRorCC ) && ( Xroc == YprocR ) &&
                           ( XprocR == Yroc   ) ) )
                     {
/*
*  If I am at the intersection of the process cross, or simply common to the
*  process rows or columns owning sub( X ) and sub( Y )
*/
                        if( ( YmyprocD == Yroc ) && ( YmyprocR == YprocR ) )
                        {
                           PB_CVMswp( type, &VM, ROW, &Xscope, &tran, npq,
                                      Mptr( ((char *) X), Xii, Xjj, Xld, size ),
                                      Xlinc, Mptr( ((char *) Y), Yii, Yjj, Yld,
                                      size ), Ylinc );
                        }
                     }
                     else
                     {
/*
*  Perform the message exchange: pack the data I own, send it, receive the
*  remote data, and unpack it.
*/
                        if( ( XmyprocR == XprocR ) && ( XmyprocD == Xroc ) )
                        {
                           if( XisRow ) { Xn = npq; }
                           else         { Xm = npq; }
                           if( YisRow ) { Yn = npq; cdst = Yroc; }
                           else         { Ym = npq; rdst = Yroc; }
                           buf = PB_Cmalloc( npq * size );
                           PB_CVMpack( type, &VM, ROW, &Xscope, PACKING, NOTRAN,
                                       npq, 1, one, Mptr( ((char *) X), Xii,
                                       Xjj, Xld, size ), Xld, zero, buf, Xm );
                           Csgesd2d( ctxt, Xm, Xn, buf, Xm, rdst, cdst );
                           Csgerv2d( ctxt, Ym, Yn, buf, Ym, rdst, cdst );
                           PB_CVMpack( type, &VM, ROW, &Xscope, UNPACKING,
                                       &tran, npq, 1, zero, Mptr( ((char *) X),
                                       Xii, Xjj, Xld, size ), Xld, one, buf,
                                       Ym );
                           if( buf ) free ( buf );
                        }
                        if( ( YmyprocR == YprocR ) && ( YmyprocD == Yroc  ) )
                        {
                           if( XisRow ) { Xn = npq; csrc = Xroc; }
                           else         { Xm = npq; rsrc = Xroc; }
                           if( YisRow ) { Yn = npq; }
                           else         { Ym = npq; }
                           buf = PB_Cmalloc( npq * size );
                           PB_CVMpack( type, &VM, COLUMN, &Yscope, PACKING,
                                       NOTRAN, npq, 1, one, Mptr( ((char *) Y),
                                       Yii, Yjj, Yld, size ), Yld, zero, buf,
                                       Ym );
                           Csgesd2d( ctxt, Ym, Yn, buf, Ym, rsrc, csrc );
                           Csgerv2d( ctxt, Xm, Xn, buf, Xm, rsrc, csrc );
                           PB_CVMpack( type, &VM, COLUMN, &Yscope, UNPACKING,
                                       &tran, npq, 1, zero, Mptr( ((char *) Y),
                                       Yii, Yjj, Yld, size ), Yld, one, buf,
                                       Xm );
                           if( buf ) free ( buf );
                        }
                     }
                  }
               }
               p = MModAdd1( p, XnprocsD );
               q = MModAdd1( q, YnprocsD );
            }
         }
      }

      if( XisR )
      {
/*
*  Replicate sub( X ) when necessary
*/
         XnpD = PB_Cnumroc( *N, 0, Xinb1D, XnbD, XmyprocD, XprocD, XnprocsD );
         if( XnpD > 0 )
         {
            if( XisRow )
            {
               top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
               if( XmyprocR == XprocR )
                  Csgebs2d( ctxt, COLUMN, top, 1, XnpD, Mptr( ((char *) X),
                            Xii, Xjj, Xld, size ), Xld );
               else
                  Csgebr2d( ctxt, COLUMN, top, 1, XnpD, Mptr( ((char *) X),
                            Xii, Xjj, Xld, size ), Xld, XprocR, XmyprocD );
            }
            else
            {
               top = PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
               if( XmyprocR == XprocR )
                  Csgebs2d( ctxt, ROW,    top, XnpD, 1, Mptr( ((char *) X),
                            Xii, Xjj, Xld, size ), Xld );
               else
                  Csgebr2d( ctxt, ROW,    top, XnpD, 1, Mptr( ((char *) X),
                            Xii, Xjj, Xld, size ), Xld, XmyprocD, XprocR );
            }
         }
      }

      if( YisR )
      {
/*
*  Replicate sub( Y ) when necessary
*/
         YnpD = PB_Cnumroc( *N, 0, Yinb1D, YnbD, YmyprocD, YprocD, YnprocsD );
         if( YnpD > 0 )
         {
            if( YisRow )
            {
               top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
               if( YmyprocR == YprocR )
                  Csgebs2d( ctxt, COLUMN, top, 1, YnpD, Mptr( ((char *) Y),
                            Yii, Yjj, Yld, size ), Yld );
               else
                  Csgebr2d( ctxt, COLUMN, top, 1, YnpD, Mptr( ((char *) Y),
                            Yii, Yjj, Yld, size ), Yld, YprocR, YmyprocD );
            }
            else
            {
               top = PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
               if( YmyprocR == YprocR )
                  Csgebs2d( ctxt, ROW,    top, YnpD, 1, Mptr( ((char *) Y),
                            Yii, Yjj, Yld, size ), Yld );
               else
                  Csgebr2d( ctxt, ROW,    top, YnpD, 1, Mptr( ((char *) Y),
                            Yii, Yjj, Yld, size ), Yld, YmyprocD, YprocR );
            }
         }
      }
   }
   else if( !( XisD ) && YisD )
   {
/*
*  sub( X ) is not distributed and sub( Y ) is distributed.
*/
      PB_CpswapND( PB_Cstypeset(), *N, ((char *) X), Xi, Xj, Xd, *INCX,
                   ((char *) Y), Yi, Yj, Yd, *INCY );
   }
   else if( XisD && !( YisD ) )
   {
/*
*  sub( X ) is distributed and sub( Y ) is not distributed.
*/
      PB_CpswapND( PB_Cstypeset(), *N, ((char *) Y), Yi, Yj, Yd, *INCY,
                   ((char *) X), Xi, Xj, Xd, *INCX );
   }
   else
   {
/*
*  Neither sub( X ) nor sub( Y ) are distributed.
*/
      PB_CpswapNN( PB_Cstypeset(), *N, ((char *) X), Xi, Xj, Xd, *INCX,
                   ((char *) Y), Yi, Yj, Yd, *INCY );
   }
/*
*  End of PSSWAP
*/
}
