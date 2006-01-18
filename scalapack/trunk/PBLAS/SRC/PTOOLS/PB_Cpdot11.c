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
void PB_Cpdot11( PBTYP_T * TYPE, int N, char * DOT,
                 char * X, int IX, int JX, int * DESCX, int INCX,
                 char * Y, int IY, int JY, int * DESCY, int INCY,
                 VVDOT_T FDOT )
#else
void PB_Cpdot11( TYPE, N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY,
                 INCY, FDOT )
/*
*  .. Scalar Arguments ..
*/
   int            INCX, INCY, IX, IY, JX, JY, N;
   char           * DOT;
   PBTYP_T        * TYPE;
   VVDOT_T        FDOT;
/*
*  .. Array Arguments ..
*/
   int            * DESCX, * DESCY;
   char           * X, * Y;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cpdot11 forms the dot product of two subvectors,
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
*  One subvector at least is assumed to be degenerated.
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
   char           Xscope, Yscope, * top;
   int            RRorCC, Xcol, Xii, XisD, XisOne, XisR, XisRow, Xjj, Xld,
                  Xlinc, XmyprocD, XmyprocR, XprocD, XprocR, Xrow, Ycol, Yii,
                  YisD, YisOne, YisR, YisRow, Yjj, YmyprocD, YmyprocR, YprocD,
                  YprocR, Yrow, cdst, ctxt, ione=1, mycol, myrow, npcol, nprow,
                  rdst;
/*
*  .. Local Arrays ..
*/
   int            dbuf[DLEN_];
   char           * buf = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCX[ CTXT_ ] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Retrieve sub( X )'s local information: Xii, Xjj, Xrow, Xcol ...
*/
   PB_Cinfog2l( IX, JX, DESCX, nprow, npcol, myrow, mycol, &Xii, &Xjj,
                &Xrow, &Xcol );
   if( ( XisRow = ( INCX == DESCX[ M_ ] ) ) != 0 )
   {
      Xld    = DESCX[ LLD_  ];         Xlinc  = Xld;
      XprocD = Xcol; XmyprocD = mycol; XprocR = Xrow; XmyprocR = myrow;
      XisR   = ( ( Xrow == -1 ) || ( nprow == 1 ) );
      XisD   = ( ( Xcol >=  0 ) && ( npcol >  1 ) );
   }
   else
   {
      Xld    = DESCX[ LLD_  ];         Xlinc  = 1;
      XprocD = Xrow; XmyprocD = myrow; XprocR = Xcol; XmyprocR = mycol;
      XisR   = ( ( Xcol == -1 ) || ( npcol == 1 ) );
      XisD   = ( ( Xrow >=  0 ) && ( nprow >  1 ) );
   }
   XisOne = ( ( N == 1 ) && ( DESCX[ M_ ] == 1 ) );
/*
*  Retrieve sub( Y )'s local information: Yii, Yjj, Yrow, Ycol ...
*/
   PB_Cinfog2l( IY, JY, DESCY, nprow, npcol, myrow, mycol, &Yii, &Yjj,
                &Yrow, &Ycol );
   if( ( YisRow = ( INCY == DESCY[ M_ ] ) ) != 0 )
   {
      YprocD = Ycol; YmyprocD = mycol; YprocR = Yrow; YmyprocR = myrow;
      YisR   = ( ( Yrow == -1 ) || ( nprow == 1 ) );
      YisD   = ( ( Ycol >=  0 ) && ( npcol >  1 ) );
   }
   else
   {
      YprocD = Yrow; YmyprocD = myrow; YprocR = Ycol; YmyprocR = mycol;
      YisR   = ( ( Ycol == -1 ) || ( npcol == 1 ) );
      YisD   = ( ( Yrow >=  0 ) && ( nprow >  1 ) );
   }
   YisOne = ( ( N == 1 ) && ( DESCY[ M_ ] == 1 ) );
/*
*  Are sub( X ) and sub( Y ) both row or column vectors ?
*/
   RRorCC = ( ( XisRow && YisRow ) || ( !( XisRow ) && !( YisRow ) ) );
/*
*  Copy sub( Y ) in sub( X )'s scope
*/
   PB_Cdescset( dbuf, 1, 1, 1, 1, 1, 1, Xrow, Xcol, ctxt, 1 );
   buf = PB_Cmalloc( TYPE->size );
   PB_Cpaxpby( TYPE, NOCONJG, 1, 1, TYPE->one, Y, IY, JY, DESCY, ( YisRow ?
               ROW : COLUMN ), TYPE->zero, buf, 0, 0, dbuf, ( XisRow ? ROW :
               COLUMN ) );
/*
*  Compute the dot product in sub( X )'s scope
*/
   if( XisR || ( XmyprocR == XprocR ) )
   {
      if( ( XisD && ( XmyprocD == XprocD ) ) || ( !XisD ) )
         FDOT( &ione, DOT, Mptr( X, Xii, Xjj, Xld, TYPE->size ), &Xlinc, buf,
               &ione );
      if( XisD && !XisOne )
      {
         Xscope = ( XisRow ? CROW : CCOLUMN );
         top = PB_Ctop( &ctxt, COMBINE, &Xscope, TOP_GET );
         TYPE->Cgsum2d( ctxt, &Xscope, top, 1, 1, DOT, 1, -1, 0 );
      }
   }
   if( buf ) free( buf );
/*
*  sub( X ) or sub( Y ) is a degenerated vector
*/
   if( XisD && XisOne )
   {
/*
*  Since XisOne, sub( X ) must be a row vector
*/
      if( XisR )
      {
/*
*  sub( X ) resides in one process column ( *, XprocD )
*/
         if( RRorCC )
         {
/*
*  sub( Y ) is a row vector as well
*/
            if( YisR || YmyprocR == YprocR )
            {
/*
*  I am a process row owning sub( Y )
*/
               if( YisD && YisOne )
               {
/*
*  sub( Y ) resides in a process column ( *, YprocD )
*/
                  if( XprocD != YprocD )
                  {
                     if( XmyprocD == XprocD )
                        TYPE->Cgesd2d( ctxt, 1, 1, DOT, 1, XmyprocR, YprocD );
                     else if( YmyprocD == YprocD )
                        TYPE->Cgerv2d( ctxt, 1, 1, DOT, 1, XmyprocR, XprocD );
                  }
               }
               else
               {
/*
*  Every process in those rows needs the answer
*/
                  top = PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );
                  if( XmyprocD == XprocD )
                     TYPE->Cgebs2d( ctxt, ROW, top, 1, 1, DOT, 1 );
                  else
                     TYPE->Cgebr2d( ctxt, ROW, top, 1, 1, DOT, 1, XmyprocR,
                                    XprocD );
               }
            }
         }
         else
         {
/*
*  sub( Y ) is a column vector
*/
            if( YisR )
            {
/*
*              sub( Y ) resides in every process column
*/
               top = PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );
               if( XmyprocD == XprocD )
                  TYPE->Cgebs2d( ctxt, ROW, top, 1, 1, DOT, 1 );
               else
                  TYPE->Cgebr2d( ctxt, ROW, top, 1, 1, DOT, 1, XmyprocR,
                                 XprocD );
            }
            else if( XprocD != YprocR )
            {
/*
*  sub( Y ) resides in process column YprocR
*/
               if( XmyprocD == XprocD )
                  TYPE->Cgesd2d( ctxt, 1, 1, DOT, 1, XmyprocR, YprocR );
               if( YmyprocR == YprocR )
                  TYPE->Cgerv2d( ctxt, 1, 1, DOT, 1, XmyprocR, XprocD );
            }
         }
      }
      else
      {
/*
*  sub( X ) resides in one process ( XprocR, XprocD )
*/
         if( YisD && YisOne )
         {
/*
*  sub( Y ) resides in one process ( YprocR, YprocD ) if it is not replicated,
*  and in one process column ( *, YprocD ) otherwise
*/
            if( ( XprocD != YprocD ) || ( !YisR && ( XprocR != YprocR ) ) )
            {
/*
*  ( XprocR, XprocD ) sends DOT to ( YprocR, YprocD ) if sub( Y ) is not repli-
*  cated, and to ( XprocR, YprocD ) otherwise
*/
               rdst = ( YisR ? XprocR : YprocR );
               if( ( XmyprocR == XprocR ) && ( XmyprocD == XprocD ) )
                  TYPE->Cgesd2d( ctxt, 1, 1, DOT, 1, rdst, YprocD );
               if( ( YmyprocR == rdst ) && ( YmyprocD == YprocD ) )
                  TYPE->Cgerv2d( ctxt, 1, 1, DOT, 1, XprocR, XprocD );
            }

            if( YisR && ( YmyprocD == YprocD ) )
            {
/*
*  Broadcast DOT within process column owning sub( Y )
*/
               top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
               if( XmyprocR == XprocR )
                  TYPE->Cgebs2d( ctxt, COLUMN, top, 1, 1, DOT, 1 );
               else
                  TYPE->Cgebr2d( ctxt, COLUMN, top, 1, 1, DOT, 1, XprocR,
                                 YprocD );
            }
         }
         else if( !YisR )
         {
/*
*  sub( Y ) resides in one process row or column
*/
            if( YisRow ) { Yscope = CROW; rdst = YprocR; cdst = XprocD; }
            else      { Yscope = CCOLUMN; rdst = XprocR; cdst = YprocR; }
/*
*  ( XprocR, XprocD ) sends DOT to ( YprocR, XprocD ) if sub( Y ) is a row
*  vector and to ( XprocR, YprocR ) otherwise. If RRorCC, then YisRow and the
*  send occurs iff XprocR != YprocR; Otherwise !YisRow, and the send occurs
*  iff XprocD is not YprocR.
*/
            if( (    RRorCC   && ( XprocR != YprocR ) ) ||
                ( !( RRorCC ) && ( XprocD != YprocR ) ) )
            {
               if( ( XmyprocR == XprocR ) && ( XmyprocD == XprocD ) )
                  TYPE->Cgesd2d( ctxt, 1, 1, DOT, 1, rdst, cdst );
               if( ( myrow == rdst ) && ( mycol == cdst ) )
                  TYPE->Cgerv2d( ctxt, 1, 1, DOT, 1, XprocR, XprocD );
            }
/*
*  Broadcast the result in sub( Y )'s scope
*/
            if( ( myrow == rdst ) && ( mycol == cdst ) )
            {
               top = PB_Ctop( &ctxt, BCAST, &Yscope, TOP_GET );
               TYPE->Cgebs2d( ctxt, &Yscope, top, 1, 1, DOT, 1 );
            }
            else if( (    YisRow   && ( myrow == rdst ) ) ||
                     ( !( YisRow ) && ( mycol == cdst ) ) )
            {
               top = PB_Ctop( &ctxt, BCAST, &Yscope, TOP_GET );
               TYPE->Cgebr2d( ctxt, &Yscope, top, 1, 1, DOT, 1, rdst, cdst );
            }
         }
         else
         {
/*
*  Every process in the grid needs the answer
*/
            top = PB_Ctop( &ctxt, BCAST, ALL, TOP_GET );
            if( ( XmyprocR == XprocR ) && ( XmyprocD == XprocD ) )
            {
               TYPE->Cgebs2d( ctxt, ALL, top, 1, 1, DOT, 1 );
            }
            else
            {
               TYPE->Cgebr2d( ctxt, ALL, top, 1, 1, DOT, 1, XprocR, XprocD );
            }
         }
      }
   }
   else
   {
/*
*  If XisR, then the result has already been sent in every process of the grid
*/
      if( XisR ) return;

      if( RRorCC )
      {
/*
*  If YisD && YisOne => YisRow => XisRow, communication orthogonal to sub( X )'s
*  direction: only process column YprocD is involved.
*/
         if( YisD && YisOne && ( YmyprocD != YprocD ) ) return;

         if( YisR )
         {
/*
*  YisR and sub( Y ) is // to sub( X ) => bcast orthogonal to sub( X ) direction
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
               top = PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );
               if( XmyprocR == XprocR )
                  TYPE->Cgebs2d( ctxt, ROW, top, 1, 1, DOT, 1 );
               else
                  TYPE->Cgebr2d( ctxt, ROW, top, 1, 1, DOT, 1, XmyprocD,
                                 XprocR );
            }
         }
         else if( XprocR != YprocR )
         {
/*
*  Send from one column/row to another if they differ
*/
            if( XisRow )
            {
               if( XmyprocR == XprocR )
                  TYPE->Cgesd2d( ctxt, 1, 1, DOT, 1, YprocR, YmyprocD );
               if( YmyprocR == YprocR )
                  TYPE->Cgerv2d( ctxt, 1, 1, DOT, 1, XprocR, XmyprocD );
            }
            else
            {
               if( XmyprocR == XprocR )
                  TYPE->Cgesd2d( ctxt, 1, 1, DOT, 1, YmyprocD, YprocR );
               if( YmyprocR == YprocR )
                  TYPE->Cgerv2d( ctxt, 1, 1, DOT, 1, XmyprocD, XprocR );
            }
         }
      }
      else
      {
/*
*  If XisRow then !YisRow and thus bcast result in all rows if YisR or in
*  process row YprocR otherwise. If !YisD || ( YisD && !YisOne ), then result
*  should be sent in the same processes because they span a row or a column of
*  the grid.
*/
         if( XisRow || !( YisD ) || ( YisD && !( YisOne ) ) )
         {
            if( YisR || YmyprocR == YprocR )
            {
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
                  top = PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );
                  if( XmyprocR == XprocR )
                     TYPE->Cgebs2d( ctxt, ROW, top, 1, 1, DOT, 1 );
                  else
                     TYPE->Cgebr2d( ctxt, ROW, top, 1, 1, DOT, 1, XmyprocD,
                                    XprocR );
               }
            }
         }
         else if( XprocR != YprocD )
         {
/*
*  YisD && YisOne => YisRow => !XisRow, so the column of processes owning
*  sub( X ) send the result to the column YprocD. The process rows involved
*  in the operation depend on YisR.
*/
            if( YisR || YmyprocR == YprocR )
            {
               if( XmyprocR == XprocR )
               {
                  TYPE->Cgesd2d( ctxt, 1, 1, DOT, 1, YmyprocR, YprocD );
               }
               if( YmyprocD == YprocD )
               {
                  TYPE->Cgerv2d( ctxt, 1, 1, DOT, 1, YmyprocR, XprocR );
               }
            }
         }
      }
   }
/*
*  End of PB_Cpdot11
*/
}
