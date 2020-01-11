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
void pzdotu_( Int * N,
              double * DOT,
              double * X, Int * IX, Int * JX, Int * DESCX, Int * INCX,
              double * Y, Int * IY, Int * JY, Int * DESCY, Int * INCY )
#else
void pzdotu_( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
/*
*  .. Scalar Arguments ..
*/
   Int            * INCX, * INCY, * IX, * IY, * JX, * JY, * N;
   double         * DOT;
/*
*  .. Array Arguments ..
*/
   Int            * DESCX, * DESCY;
   double         * X, * Y;
#endif
{
/*
*  Purpose
*  =======
*
*  PZDOTU  forms the dot product of two subvectors,
*
*     DOT := sub( X )**T * sub( Y ),
*
*  where
*
*     sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
*                      X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X, and,
*
*     sub( Y ) denotes Y(IY,JY:JY+N-1) if INCY = M_Y,
*                      Y(IY:IY+N-1,JY) if INCY = 1 and INCY <> M_Y.
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
*          multiplied. N must be at least zero.
*
*  DOT     (local output) COMPLEX*16 array
*          On exit, DOT  specifies the dot product of the two subvectors
*          sub( X ) and sub( Y ) only in their scope (See below for fur-
*          ther details).
*
*  X       (local input) COMPLEX*16 array
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
*  Y       (local input) COMPLEX*16 array
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
   char           scope, * top;
   Int            OneBlock, OneDgrid, RRorCC, Square, Xcol, Xi, Xii, XinbD,
                  Xinb1D, XisD, XisR, XisRow, Xj, Xjj, Xld, Xlinc, XmyprocD,
                  XmyprocR, XnbD, XnpD, XnprocsD, XnprocsR, XprocD, XprocR,
                  Xrow, Ycol, Yi, Yii, YinbD, Yinb1D, YisD, YisR, YisRow, Yj,
                  Yjj, Yld, Ylinc, YmyprocD, YmyprocR, YnbD, YnpD, YnprocsD,
                  YnprocsR, YprocD, YprocR, Yrow, cdst, csrc, ctxt, dst, info,
                  ione=1, mycol, myrow, npcol, nprow, rdst, rsrc, size, src;
   PBTYP_T        * type;
   VVDOT_T        dot;
/*
*  .. Local Arrays ..
*/
   char           * buf = NULL;
   Int            Xd[DLEN_], Yd[DLEN_], dbuf[ DLEN_ ];
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
   if( !( info = ( ( nprow == -1 ) ? -( 601 + CTXT_ ) : 0 ) ) )
   {
      PB_Cchkvec( ctxt, "PZDOTU", "X", *N, 1, Xi, Xj, Xd, *INCX,  6, &info );
      PB_Cchkvec( ctxt, "PZDOTU", "Y", *N, 1, Yi, Yj, Yd, *INCY, 11, &info );
   }
   if( info ) { PB_Cabort( ctxt, "PZDOTU", info ); return; }
#endif
   DOT[REAL_PART] = ZERO;
   DOT[IMAG_PART] = ZERO;
/*
*  Quick return if possible
*/
   if( *N == 0 ) return;
/*
*  Handle degenerate case
*/
   if( ( *N == 1 ) && ( ( Xd[ M_ ] == 1 ) || ( Yd[ M_ ] == 1 ) ) )
   {
      type = PB_Cztypeset();
      PB_Cpdot11( type, *N, ((char *) DOT), ((char *) X), Xi, Xj, Xd, *INCX,
                  ((char *) Y), Yi, Yj, Yd, *INCY, type->Fvvdotu );
      return;
   }
/*
*  Start the operations
*/
#ifdef NO_ARGCHK
   Cblacs_gridinfo( ( ctxt = Xd[ CTXT_ ] ), &nprow, &npcol, &myrow, &mycol );
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
*  Retrieve sub( X )'s local information: Xii, Xjj, Xrow, Xcol
*/
      PB_Cinfog2l( Xi, Xj, Xd, nprow, npcol, myrow, mycol, &Xii, &Xjj,
                   &Xrow, &Xcol );
      if( XisRow )
      {
         XinbD  = Xd[INB_]; XnbD     = Xd[NB_];
         Xld    = Xd[LLD_]; Xlinc    = Xld;
         XprocD = Xcol;     XmyprocD = mycol;
         XprocR = Xrow;     XmyprocR = myrow;   XnprocsR = nprow;
         XisR   = ( ( Xrow == -1 ) || ( XnprocsR == 1 ) );
         Mfirstnb( Xinb1D, *N, Xj, XinbD, XnbD );
      }
      else
      {
         XinbD  = Xd[IMB_]; XnbD     = Xd[MB_];
         Xld    = Xd[LLD_]; Xlinc    = 1;
         XprocD = Xrow;     XmyprocD = myrow;
         XprocR = Xcol;     XmyprocR = mycol;   XnprocsR = npcol;
         XisR   = ( ( Xcol == -1 ) || ( XnprocsR == 1 ) );
         Mfirstnb( Xinb1D, *N, Xi, XinbD, XnbD );
      }
/*
*  Retrieve sub( Y )'s local information: Yii, Yjj, Yrow, Ycol
*/
      PB_Cinfog2l( Yi, Yj, Yd, nprow, npcol, myrow, mycol, &Yii, &Yjj,
                   &Yrow, &Ycol );
      if( YisRow )
      {
         YinbD  = Yd[INB_]; YnbD     = Yd[NB_];
         Yld    = Yd[LLD_]; Ylinc    = Yld;
         YprocD = Ycol;     YmyprocD = mycol;
         YprocR = Yrow;     YmyprocR = myrow;   YnprocsR = nprow;
         YisR   = ( ( Yrow == -1 ) || ( YnprocsR == 1 ) );
         Mfirstnb( Yinb1D, *N, Yj, YinbD, YnbD );
      }
      else
      {
         YinbD  = Yd[IMB_]; YnbD     = Yd[MB_];
         Yld    = Yd[LLD_]; Ylinc    = 1;
         YprocD = Yrow;     YmyprocD = myrow;
         YprocR = Ycol;     YmyprocR = mycol;   YnprocsR = npcol;
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
*  sub( X ) is not replicated
*/
         if( YisR )
         {
/*
*  If sub( X ) is not replicated, but sub( Y ) is, a process row or column
*  YprocR need to be selected. It will contain the non-replicated vector used
*  to perform the dot product computation.
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
      type = PB_Cztypeset(); size = type->size; dot = type->Fvvdotu;
/*
*  A purely operation occurs iff the operands start in the same process and if
*  either the grid is mono-dimensional or there is a single local block to be
*  operated with or if both operands are aligned.
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
               dot( &XnpD, ((char *) DOT), Mptr( ((char *) X), Xii, Xjj, Xld,
                    size ), &Xlinc, Mptr( ((char *) Y), Yii, Yjj, Yld, size ),
                    &Ylinc );
            }
         }
/*
*  Combine the local results in sub( X )'s scope
*/
         if( ( XisR && YisR ) || ( XmyprocR == XprocR ) )
         {
            scope = ( XisRow ? CROW : CCOLUMN );
            top = PB_Ctop( &ctxt, COMBINE, &scope, TOP_GET );
            Czgsum2d( ctxt, &scope, top, 1, 1, ((char *) DOT), 1, -1, 0 );
         }
         if( RRorCC && XisR && YisR ) return;
      }
      else if( ( RRorCC && OneDgrid ) || OneBlock || Square )
      {
/*
*  Otherwise, it may be possible to compute the desired dot-product in a single
*  message exchange iff the grid is mono-dimensional and the operands are
*  distributed in the same direction, or there is just one block to be exchanged
*  or if both operands are similarly distributed in their respective direction.
*/
         if( ( YmyprocR == YprocR ) )
         {
/*
*  The processes owning a piece of sub( Y ) send it to the corresponding
*  process owning s piece of sub ( X ).
*/
            YnpD = PB_Cnumroc( *N, 0, Yinb1D, YnbD, YmyprocD, YprocD,
                               YnprocsD );
            if( YnpD > 0 )
            {
               dst = XprocD + MModSub( YmyprocD, YprocD, YnprocsD );
               dst = MPosMod( dst, XnprocsD );
               if( XisRow ) { rdst = XprocR; cdst = dst; }
               else         { rdst = dst; cdst = XprocR; }

               if( ( myrow == rdst ) && ( mycol == cdst ) )
               {
                  dot( &YnpD, ((char *) DOT), Mptr( ((char *) X), Xii, Xjj, Xld,
                       size ), &Xlinc, Mptr( ((char *) Y), Yii, Yjj, Yld,
                       size ), &Ylinc );
               }
               else
               {
                  if( YisRow )
                     Czgesd2d( ctxt, 1, YnpD, Mptr( ((char *) Y), Yii, Yjj,
                               Yld, size ), Yld, rdst, cdst );
                  else
                     Czgesd2d( ctxt, YnpD, 1, Mptr( ((char *) Y), Yii, Yjj,
                               Yld, size ), Yld, rdst, cdst );
               }
            }
         }
         if( XmyprocR == XprocR )
         {
/*
*  The processes owning a piece of sub( X ) receive the corresponding local
*  piece of sub( Y ), compute the local dot product and combine the results
*  within sub( X )'s scope.
*/
            XnpD = PB_Cnumroc( *N, 0, Xinb1D, XnbD, XmyprocD, XprocD,
                               XnprocsD );
            if( XnpD > 0 )
            {
               src = YprocD + MModSub( XmyprocD, XprocD, XnprocsD );
               src = MPosMod( src, YnprocsD );
               if( YisRow ) { rsrc = YprocR; csrc = src; }
               else         { rsrc = src; csrc = YprocR; }
               if( ( myrow != rsrc ) || ( mycol != csrc ) )
               {
                  buf = PB_Cmalloc( XnpD * size );
                  if( YisRow )
                     Czgerv2d( ctxt, 1, XnpD, buf,    1, rsrc, csrc );
                  else
                     Czgerv2d( ctxt, XnpD, 1, buf, XnpD, rsrc, csrc );
                  dot( &XnpD, ((char *) DOT), Mptr( ((char *) X), Xii, Xjj, Xld,
                       size ), &Xlinc, buf, &ione );
                  if( buf ) free( buf );
               }
            }
            if( XisRow )
            {
               top = PB_Ctop( &ctxt, COMBINE, ROW, TOP_GET );
               Czgsum2d( ctxt, ROW,    top, 1, 1, ((char*)DOT), 1, -1, 0 );
            }
            else
            {
               top = PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
               Czgsum2d( ctxt, COLUMN, top, 1, 1, ((char*)DOT), 1, -1, 0 );
            }
         }
      }
      else
      {
/*
*  General case, copy sub( Y ) within sub( X )'s scope, compute the local
*  results and combine them within sub( X )'s scope.
*/
         XnpD = PB_Cnumroc( *N, 0, Xinb1D, XnbD, XmyprocD, XprocD, XnprocsD );

         if( XisRow )
         {
            PB_Cdescset( dbuf, 1, *N, 1, Xinb1D, 1, XnbD, XprocR, XprocD, ctxt,
                         1 );
         }
         else
         {
            PB_Cdescset( dbuf, *N, 1, Xinb1D, 1, XnbD, 1, XprocD, XprocR, ctxt,
                         MAX( 1, XnpD ) );
         }
         if( ( XmyprocR == XprocR ) && ( XnpD > 0 ) )
            buf = PB_Cmalloc( XnpD * size );

         if( YisRow )
         {
            PB_Cpaxpby( type, NOCONJG, 1, *N, type->one, ((char *) Y), Yi, Yj,
                        Yd, ROW, type->zero, buf, 0, 0, dbuf, ( XisRow ? ROW :
                        COLUMN ) );
         }
         else
         {
            PB_Cpaxpby( type, NOCONJG, *N, 1, type->one, ((char *) Y), Yi, Yj,
                        Yd, COLUMN, type->zero, buf, 0, 0, dbuf, ( XisRow ?
                        ROW : COLUMN ) );
         }

         if( XmyprocR == XprocR )
         {
            if( XnpD > 0 )
            {
               dot( &XnpD, ((char *) DOT), Mptr( ((char *) X), Xii, Xjj, Xld,
                    size ), &Xlinc, buf, &ione );
               if( buf ) free( buf );
            }
            if( XisRow )
            {
               top = PB_Ctop( &ctxt, COMBINE, ROW,    TOP_GET );
               Czgsum2d( ctxt, ROW,    top, 1, 1, ((char*)DOT), 1, -1, 0 );
            }
            else
            {
               top = PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
               Czgsum2d( ctxt, COLUMN, top, 1, 1, ((char*)DOT), 1, -1, 0 );
            }
         }
      }
/*
*  Send the DOT product result within sub( Y )'s scope
*/
      if( XisR || YisR )
      {
/*
*  Either sub( X ) or sub( Y ) are replicated, so that every process should have
*  the result -> broadcast it orthogonally from sub( X )'s direction.
*/
         if( XisRow )
         {
           top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
           if( XmyprocR == XprocR )
              Czgebs2d( ctxt, COLUMN, top, 1, 1, ((char*)DOT), 1 );
           else
              Czgebr2d( ctxt, COLUMN, top, 1, 1, ((char*)DOT), 1, XprocR,
                        XmyprocD );
         }
         else
         {
           top = PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
           if( XmyprocR == XprocR )
              Czgebs2d( ctxt, ROW,    top, 1, 1, ((char*)DOT), 1 );
           else
              Czgebr2d( ctxt, ROW,    top, 1, 1, ((char*)DOT), 1, XmyprocD,
                        XprocR );
         }
      }
      else
      {
/*
*  Neither sub( X ) nor sub( Y ) are replicated
*/
         if( RRorCC )
         {
/*
*  Both sub( X ) are distributed in the same direction -> the process row or
*  column XprocR sends the result to the process row or column YprocR.
*/
            if( XprocR != YprocR )
            {
               if( XmyprocR == XprocR )
               {
                  if( XisRow )
                     Czgesd2d( ctxt, 1, 1, ((char *) DOT), 1, YprocR,
                               YmyprocD );
                  else
                     Czgesd2d( ctxt, 1, 1, ((char *) DOT), 1, YmyprocD,
                               YprocR );
               }
               else if( YmyprocR == YprocR )
               {
                  if( XisRow )
                     Czgerv2d( ctxt, 1, 1, ((char *) DOT), 1, XprocR,
                               XmyprocD );
                  else
                     Czgerv2d( ctxt, 1, 1, ((char *) DOT), 1, XmyprocD,
                               XprocR );
               }
            }
         }
         else
         {
/*
*  Otherwise, the process at the intersection of sub( X )'s and sub( Y )'s
*  scope, broadcast the result within sub( Y )'s scope.
*/
            if( YmyprocR == YprocR )
            {
               if( YisRow )
               {
                  top = PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
                  if( YmyprocD == XprocR )
                     Czgebs2d( ctxt, ROW,    top, 1, 1, ((char*)DOT), 1 );
                  else
                     Czgebr2d( ctxt, ROW,    top, 1, 1, ((char*)DOT), 1,
                               YprocR, XprocR );
               }
               else
               {
                  top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
                  if( YmyprocD == XprocR )
                     Czgebs2d( ctxt, COLUMN, top, 1, 1, ((char*)DOT), 1 );
                  else
                     Czgebr2d( ctxt, COLUMN, top, 1, 1, ((char*)DOT), 1,
                               XprocR, YprocR );
               }
            }
         }
      }
   }
   else if( !( XisD ) && YisD )
   {
/*
*  sub( X ) is not distributed and sub( Y ) is distributed.
*/
      type = PB_Cztypeset();
      PB_CpdotND( type, *N, ((char *) DOT), ((char *) X), Xi, Xj, Xd, *INCX,
                  ((char *) Y), Yi, Yj, Yd, *INCY, type->Fvvdotu );
   }
   else if( XisD && !( YisD ) )
   {
/*
*  sub( X ) is distributed and sub( Y ) is not distributed.
*/
      type = PB_Cztypeset();
      PB_CpdotND( type, *N, ((char *) DOT), ((char *) Y), Yi, Yj, Yd, *INCY,
                  ((char *) X), Xi, Xj, Xd, *INCX, type->Fvvdotu );
   }
   else
   {
/*
*     Neither sub( X ) nor sub( Y ) are distributed
*/
      type = PB_Cztypeset();
      PB_CpdotNN( type, *N, ((char *) DOT), ((char *) X), Xi, Xj, Xd, *INCX,
                  ((char *) Y), Yi, Yj, Yd, *INCY, type->Fvvdotu );
   }
/*
*  End of PZDOTU
*/
}
