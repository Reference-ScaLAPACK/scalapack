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
void PB_CpswapNN( PBTYP_T * TYPE, Int N,
                  char * X, Int IX, Int JX, Int * DESCX, Int INCX,
                  char * Y, Int IY, Int JY, Int * DESCY, Int INCY )
#else
void PB_CpswapNN( TYPE, N, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
/*
*  .. Scalar Arguments ..
*/
   Int            INCX, INCY, IX, IY, JX, JY, N;
   PBTYP_T        * TYPE;
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
*  PB_CpswapNN swaps two subvectors,
*
*     sub( Y ) := sub( X ) and sub( X ) := sub( Y )
*
*  where sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
*                         X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X,
*
*        sub( Y ) denotes Y(IY,JY:JY+N-1) if INCY = M_Y,
*                         Y(IY:IY+N-1,JY) if INCY = 1 and INCY <> M_Y.
*
*  Both subvectors are assumed to be not distributed.
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
*          swapped. N must be at least zero.
*
*  X       (local input/local output) pointer to CHAR
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
*  Y       (local input/local output) pointer to CHAR
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
   char           Xscope, Yscope, * top;
   Int            RRorCC, XYm, XYn, Xcol, Xii, XisR, XisRow, Xjj, Xld, Xlinc,
                  XmyprocD, XmyprocR, XnprocsR, XprocR, Xrow, Ycol, Yii, YisR,
                  YisRow, Yjj, Yld, Ylinc, YmyprocD, YmyprocR, YnprocsR, YprocR,
                  Yrow, csrc, ctxt, mycol, myrow, npcol, nprow, rsrc, size;
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
      Xld      = DESCX[LLD_];          Xlinc    = Xld;
      XmyprocD = mycol; XprocR = Xrow; XmyprocR = myrow; XnprocsR = nprow;
      XisR     = ( ( Xrow == -1 ) || ( XnprocsR == 1 ) );
   }
   else
   {
      Xld      = DESCX[LLD_];          Xlinc    = 1;
      XmyprocD = myrow; XprocR = Xcol; XmyprocR = mycol; XnprocsR = npcol;
      XisR     = ( ( Xcol == -1 ) || ( XnprocsR == 1 ) );
   }
/*
*  Retrieve sub( Y )'s local information: Yii, Yjj, Yrow, Ycol ...
*/
   PB_Cinfog2l( IY, JY, DESCY, nprow, npcol, myrow, mycol, &Yii, &Yjj,
                &Yrow, &Ycol );
   if( ( YisRow = ( INCY == DESCY[M_] ) ) != 0 )
   {
      Yld      = DESCY[LLD_];          Ylinc    = Yld;
      YmyprocD = mycol; YprocR = Yrow; YmyprocR = myrow; YnprocsR = nprow;
      YisR     = ( ( Yrow == -1 ) || ( YnprocsR == 1 ) );
   }
   else
   {
      Yld      = DESCY[LLD_];          Ylinc    = 1;
      YmyprocD = myrow; YprocR = Ycol; YmyprocR = mycol; YnprocsR = npcol;
      YisR     = ( ( Ycol == -1 ) || ( YnprocsR == 1 ) );
   }
/*
*  Are sub( X ) and sub( Y ) both row or column vectors ?
*/
   RRorCC = ( ( XisRow && YisRow ) || ( !( XisRow ) && !( YisRow ) ) );
/*
*  Neither sub( X ) nor sub( Y ) are distributed
*/
   if( !XisR )
   {
/*
*  sub( X ) is not replicated
*/
      if( !( YisR ) )
      {
/*
*  sub( Y ) is not replicated
*/
         if( ( XmyprocR != XprocR ) && ( YmyprocR != YprocR ) )
/*
*  If I am not in XprocR or YprocR, then return immediately
*/
            return;

         size = TYPE->size;

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
               TYPE->Fswap( &N, Mptr( X, Xii, Xjj, Xld, size ), &Xlinc, Mptr( Y,
                            Yii, Yjj, Yld, size ), &Ylinc );
            }
            else
            {
/*
*  sub( X ) and sub( Y ) are in a different process row or column
*/
               if( XmyprocR == XprocR )
               {
/*
*  Send sub( X ) to where sub( Y ) resides, and receive sub( Y ) from the same
*  location.
*/
                  if( XisRow )
                  {
                     TYPE->Cgesd2d( ctxt, 1, N, Mptr( X, Xii, Xjj, Xld, size ),
                                    Xld, YprocR, XmyprocD );
                     TYPE->Cgerv2d( ctxt, 1, N, Mptr( X, Xii, Xjj, Xld, size ),
                                    Xld, YprocR, XmyprocD );
                  }
                  else
                  {
                     TYPE->Cgesd2d( ctxt, N, 1, Mptr( X, Xii, Xjj, Xld, size ),
                                    Xld, XmyprocD, YprocR );
                     TYPE->Cgerv2d( ctxt, N, 1, Mptr( X, Xii, Xjj, Xld, size ),
                                    Xld, XmyprocD, YprocR );
                  }
               }

               if( YmyprocR == YprocR )
               {
/*
*  Send sub( Y ) to where sub( X ) resides, and receive sub( X ) from the same
*  location.
*/
                  if( YisRow )
                  {
                     TYPE->Cgesd2d( ctxt, 1, N, Mptr( Y, Yii, Yjj, Yld, size ),
                                    Yld, XprocR, YmyprocD );
                     TYPE->Cgerv2d( ctxt, 1, N, Mptr( Y, Yii, Yjj, Yld, size ),
                                    Yld, XprocR, YmyprocD );
                  }
                  else
                  {
                     TYPE->Cgesd2d( ctxt, N, 1, Mptr( Y, Yii, Yjj, Yld, size ),
                                    Yld, YmyprocD, XprocR );
                     TYPE->Cgerv2d( ctxt, N, 1, Mptr( Y, Yii, Yjj, Yld, size ),
                                    Yld, YmyprocD, XprocR );
                  }
               }
            }
         }
         else
         {
/*
*  sub( X ) and sub( Y ) are not both row or column vectors
*/
            if( XisRow )
            {
              XYm    = 1;        XYn    = N;
              Xscope = CROW;     Yscope = CCOLUMN;
              rsrc   = XprocR;   csrc   = YprocR;
            }
            else
            {
               XYm    = N;       XYn    = 1;
               Xscope = CCOLUMN; Yscope = CROW;
               rsrc   = YprocR;  csrc   = XprocR;
            }

            if( ( XmyprocR == XprocR ) && ( YmyprocR == YprocR ) )
            {
/*
*  If I am at the intersection of the process row and column, then swap and
*  broadcast sub( X ) and sub( Y ) in their respective process scope.
*/
               TYPE->Fswap( &N, Mptr( X, Xii, Xjj, Xld, size ), &Xlinc,
                            Mptr( Y, Yii, Yjj, Yld, size ), &Ylinc );
               top = PB_Ctop( &ctxt, BCAST, &Xscope, TOP_GET );
               TYPE->Cgebs2d( ctxt, &Xscope, top, XYm, XYn, Mptr( X, Xii, Xjj,
                              Xld, size ), Xld );
               top = PB_Ctop( &ctxt, BCAST, &Yscope, TOP_GET );
               TYPE->Cgebs2d( ctxt, &Yscope, top, XYn, XYm, Mptr( Y, Yii, Yjj,
                              Yld, size ), Yld );
            }
            else if( XmyprocR == XprocR )
            {
               top = PB_Ctop( &ctxt, BCAST, &Xscope, TOP_GET );
               TYPE->Cgebr2d( ctxt, &Xscope, top, XYm, XYn, Mptr( X, Xii, Xjj,
                              Xld, size ), Xld, rsrc, csrc );
            }
            else if( YmyprocR == YprocR )
            {
               top = PB_Ctop( &ctxt, BCAST, &Yscope, TOP_GET );
               TYPE->Cgebr2d( ctxt, &Yscope, top, XYn, XYm, Mptr( Y, Yii, Yjj,
                              Yld, size ), Yld, rsrc, csrc );
            }
         }
      }
      else
      {
/*
*  sub( Y ) is replicated
*/
         size = TYPE->size;

         if( YisRow ) { XYm = 1; XYn = N; }
         else         { XYm = N; XYn = 1; }

         if( XmyprocR == XprocR )
         {
/*
*  If I am in the process row (resp. column) owning sub( X ), then swap and
*  broadcast sub( Y ) in my column (resp. row).
*/
            TYPE->Fswap( &N, Mptr( X, Xii, Xjj, Xld, size ), &Xlinc, Mptr( Y,
                         Yii, Yjj, Yld, size ), &Ylinc );

            if( XisRow )
            {
               top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
               TYPE->Cgebs2d( ctxt, COLUMN, top, XYm, XYn, Mptr( Y, Yii, Yjj,
                              Yld, size ), Yld );
            }
            else
            {
               top = PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
               TYPE->Cgebs2d( ctxt, ROW,    top, XYm, XYn, Mptr( Y, Yii, Yjj,
                              Yld, size ), Yld );
            }
         }
         else
         {
/*
*  Otherwise, receive sub( Y )
*/
            if( XisRow )
            {
               top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
               TYPE->Cgebr2d( ctxt, COLUMN, top, XYm, XYn, Mptr( Y, Yii, Yjj,
                              Yld, size ), Yld, XprocR, XmyprocD );
            }
            else
            {
               top = PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
               TYPE->Cgebr2d( ctxt, ROW,    top, XYm, XYn, Mptr( Y, Yii, Yjj,
                              Yld, size ), Yld, XmyprocD, XprocR );
            }
         }
      }
   }
   else
   {
/*
*  sub( X ) is replicated
*/
      size = TYPE->size;

      if( YisR || ( YmyprocR == YprocR ) )
      {
/*
*  If I own a piece of sub( Y ), then swap
*/
         TYPE->Fswap( &N, Mptr( X, Xii, Xjj, Xld, size ), &Xlinc, Mptr( Y, Yii,
                      Yjj, Yld, size ), &Ylinc );
      }

      if( !YisR )
      {
/*
*  If sub( Y ) is not replicated, then broadcast the result to the other
*  processes that own a piece of sub( X ), but were not involved in the
*  above swap operation.
*/
         if( XisRow ) { XYm = 1; XYn = N; }
         else         { XYm = N; XYn = 1; }

         if( YisRow )
         {
            top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
            if( YmyprocR == YprocR )
               TYPE->Cgebs2d( ctxt, COLUMN, top, XYm, XYn, Mptr( X, Xii, Xjj,
                              Xld, size ), Xld );
            else
               TYPE->Cgebr2d( ctxt, COLUMN, top, XYm, XYn, Mptr( X, Xii, Xjj,
                              Xld, size ), Xld, YprocR, YmyprocD );
         }
         else
         {
            top = PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
            if( YmyprocR == YprocR )
               TYPE->Cgebs2d( ctxt, ROW,    top, XYm, XYn, Mptr( X, Xii, Xjj,
                              Xld, size ), Xld );
            else
               TYPE->Cgebr2d( ctxt, ROW,    top, XYm, XYn, Mptr( X, Xii, Xjj,
                              Xld, size ), Xld, YmyprocD, YprocR );
         }
      }
   }
/*
*  End of PB_CpswapNN
*/
}
