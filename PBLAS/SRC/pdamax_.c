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
void pdamax_( int * N, double * AMAX, int * INDX,
              double * X, int * IX, int * JX, int * DESCX, int * INCX )
#else
void pdamax_( N, AMAX, INDX, X, IX, JX, DESCX, INCX )
/*
*  .. Scalar Arguments ..
*/
   int            * INCX, * INDX, * IX, * JX, * N;
   double         * AMAX;
/*
*  .. Array Arguments ..
*/
   int            * DESCX;
   double         * X;
#endif
{
/*
*  Purpose
*  =======
*
*  PDAMAX  computes the global index of the maximum element in  absolute
*  value of a subvector sub( X ).  The global index is returned in  INDX
*  and the value of that element is returned in AMAX,
*
*  where
*
*     sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
*                      X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X.
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
*          On entry,  N  specifies the length of the subvector sub( X ).
*          N must be at least zero.
*
*  AMAX    (global output) DOUBLE PRECISION array
*          On exit,  AMAX  specifies the largest entry in absolute value
*          of the  subvector  sub( X )  only in its scope (See below for
*          further details).
*
*  INDX    (global output) INTEGER
*          On exit, INDX  specifies the global index of the maximum ele-
*          ment in absolute  value of the subvector sub( X ) only in its
*          scope (See below for further details).
*
*  X       (local input) DOUBLE PRECISION array
*          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
*          is   at  least  MAX( 1, Lr( 1, IX ) )  when  INCX = M_X   and
*          MAX( 1, Lr( 1, IX+N-1 ) )  otherwise,  and,  Kx  is  at least
*          Lc( 1, JX+N-1 )  when  INCX = M_X  and Lc( 1, JX ) otherwise.
*          Before  entry,  this array  contains the local entries of the
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
   char           cbtop, cctop, rbtop, rctop;
   int            Xcol, Xgindx, Xi, Xii, Ximb, Xinb, Xj, Xjj, Xlindx, Xld, Xmb,
                  Xnb, Xnp, Xnq, Xrow, Xsrc, ctxt, dist, dst, idumm, info, k,
                  maxpos, mycol, mydist, myrow, npcol, nprow, src;
/*
*  .. Local Arrays ..
*/
   int            Xd[DLEN_];
   double         work[4];
/* ..
*  .. Executable Statements ..
*
*/
   PB_CargFtoC( *IX, *JX, DESCX, &Xi, &Xj, Xd );
#ifndef NO_ARGCHK
/*
*  Test the input parameters
*/
   Cblacs_gridinfo( ( ctxt = Xd[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
   if( !( info = ( ( nprow == -1 ) ? -( 701 + CTXT_ ) : 0 ) ) )
      PB_Cchkvec( ctxt, "PDAMAX", "X", *N, 1, Xi, Xj, Xd, *INCX, 7, &info );
   if( info ) { PB_Cabort( ctxt, "PDAMAX", info ); return; }
#endif
/*
*  Initialize INDX and AMAX
*/
   *INDX = 0; *AMAX = ZERO;
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
*  Retrieve sub( X )'s local information: Xii, Xjj, Xrow, Xcol
*/
   PB_Cinfog2l( Xi, Xj, Xd, nprow, npcol, myrow, mycol, &Xii, &Xjj,
                &Xrow, &Xcol );
/*
*  Handle degenerate case separately, sub( X )'s scope is just one process
*/
   if( ( *INCX == 1 ) && ( Xd[M_] == 1 ) && ( *N == 1 ) )
   {
/*
*  Make sure I own some data and compute INDX and AMAX
*/
      if( ( ( myrow == Xrow ) || ( Xrow < 0 ) ) &&
          ( ( mycol == Xcol ) || ( Xcol < 0 ) ) )
      {
         *INDX = *JX; *AMAX = *Mptr(X,Xii,Xjj,Xd[LLD_],1);
      }
      return;
   }
   else if( *INCX == Xd[M_] )
   {
/*
*  sub( X ) resides in (a) process row(s)
*/
      if( ( myrow == Xrow ) || ( Xrow < 0 ) )
      {
         rctop = *PB_Ctop( &ctxt, COMBINE, ROW, TOP_GET );

         if( ( rctop == CTOP_DEFAULT ) || ( rctop == CTOP_TREE1 ) )
         {
/*
*  Inline the 1-tree combine for communication savings
*/
            Xinb = Xd[INB_ ]; Xnb = Xd[NB_ ]; Xsrc = Xd[CSRC_];
            Xnq = PB_Cnumroc( *N, Xj, Xinb, Xnb, mycol, Xsrc, npcol );
/*
*  Make sure I own some data and compute local INDX and AMAX
*/
            if( Xnq > 0 )
            {
               Xld = Xd[LLD_];
               Xlindx = Xjj - 1 +
                        idamax_( &Xnq, ((char*)Mptr(X,Xii,Xjj,Xld,1)), &Xld );
               Mindxl2g( Xgindx, Xlindx, Xinb, Xnb, mycol, Xsrc, npcol );
               work[0] = *Mptr(X,Xii,Xlindx,Xld,1);
               work[1] = ((double)( Xgindx+1 ));
            }
            else
            {
               work[0] = ZERO;
               work[1] = ZERO;
            }
/*
*  Combine the local results using a 1-tree topology within process column 0
*  if npcol > 1 or Xcol >= 0, i.e sub( X ) is distributed.
*/
            if( ( npcol >= 2 ) && ( Xcol >= 0 ) )
            {
               mydist = mycol;
               k      = 1;
l_10:
               if( mydist & 1 )
               {
                  dist = k * ( mydist - 1 );
                  dst  = MPosMod( dist, npcol );
                  Cdgesd2d( ctxt, 2, 1, ((char*)work), 2, myrow, dst );
                  goto l_20;
               }
               else
               {
                  dist = mycol + k;
                  src  = MPosMod( dist, npcol );

                  if( mycol < src )
                  {
                     Cdgerv2d( ctxt, 2, 1, ((char*) &work[2]), 2, myrow,
                               src );
                     if( ABS( work[0] ) < ABS( work[2] ) )
                     { work[0] = work[2]; work[1] = work[3]; }
                  }
                  mydist >>= 1;
               }
               k <<= 1;

               if( k < npcol ) goto l_10;
l_20:
/*
*  Process column 0 broadcasts the combined values of INDX and AMAX within
*  their process row.
*/
               rbtop = *PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );
               if( mycol == 0 )
               {
                  Cdgebs2d( ctxt, ROW, &rbtop, 2, 1, ((char*)work), 2 );
               }
               else
               {
                  Cdgebr2d( ctxt, ROW, &rbtop, 2, 1, ((char*)work), 2,
                            myrow, 0 );
               }
            }
/*
*  Set INDX and AMAX to the replicated answers contained in work. If AMAX is
*  zero, then select a coherent INDX.
*/
            *AMAX = work[0];
            *INDX = ( ( *AMAX == ZERO ) ? ( *JX ) : ( (int)(work[1]) ) );
         }
         else
         {
/*
*  Otherwise use the current topology settings to combine the results
*/
            Xinb = Xd[INB_ ]; Xnb = Xd[NB_ ]; Xsrc = Xd[CSRC_];
            Xnq = PB_Cnumroc( *N, Xj, Xinb, Xnb, mycol, Xsrc, npcol );
/*
*  Make sure I own some data and compute local INDX and AMAX
*/
            if( Xnq > 0 )
            {
/*
*  Compute the local maximum and its corresponding local index
*/
               Xld = Xd[LLD_];
               Xlindx = Xjj - 1 +
                        idamax_( &Xnq, ((char*)Mptr(X,Xii,Xjj,Xld,1)), &Xld );
               *AMAX = *Mptr(X,Xii,Xlindx,Xld,1);
            }
            else
            {
               *AMAX = ZERO;
            }

            if( Xcol >= 0 )
            {
/*
*  Combine leave on all the local maximum if Xcol >= 0, i.e sub( X ) is
*  distributed
*/
               Cdgamx2d( ctxt, ROW, &rctop, 1, 1, ((char*)AMAX), 1,
                         &idumm, &maxpos, 1, -1, mycol );
/*
*  Broadcast the corresponding global index
*/
               if( *AMAX != ZERO )
               {
                  rbtop = *PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );
                  if( mycol == maxpos )
                  {
                     Mindxl2g( Xgindx, Xlindx, Xinb, Xnb, mycol, Xsrc, npcol );
                     *INDX = Xgindx + 1;
                     Cigebs2d( ctxt, ROW, &rbtop, 1, 1, ((char*)INDX), 1 );
                  }
                  else
                  {
                     Cigebr2d( ctxt, ROW, &rbtop, 1, 1, ((char*)INDX), 1,
                               myrow, maxpos );
                  }
               }
               else
               {
/*
*  If AMAX is zero, then select a coherent INDX.
*/
                  *INDX = *JX;
               }
            }
            else
            {
/*
*  sub( X ) is not distributed. If AMAX is zero, then select a coherent INDX.
*/
               *INDX = ( ( *AMAX == ZERO ) ? ( *JX ) : Xlindx + 1 );
            }
         }
      }
      return;
   }
   else
   {
/*
*  sub( X ) resides in (a) process column(s)
*/
      if( ( mycol == Xcol ) || ( Xcol < 0 ) )
      {
         cctop = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );

         if( ( cctop == CTOP_DEFAULT ) || ( cctop == CTOP_TREE1 ) )
         {
/*
*  Inline the 1-tree combine for communication savings
*/
            Ximb = Xd[IMB_ ]; Xmb = Xd[MB_ ]; Xsrc = Xd[RSRC_];
            Xnp = PB_Cnumroc( *N, Xi, Ximb, Xmb, myrow, Xsrc, nprow );
/*
*  Make sure I own some data and compute local INDX and AMAX
*/
            if( Xnp > 0 )
            {
               Xld     = Xd[LLD_];
               Xlindx  = Xii - 1 +
                         idamax_( &Xnp, ((char*)Mptr(X,Xii,Xjj,Xld,1)), INCX );
               Mindxl2g( Xgindx, Xlindx, Ximb, Xmb, myrow, Xsrc, nprow );
               work[0] = *Mptr(X,Xlindx,Xjj,Xld,1);
               work[1] = ((double)( Xgindx+1 ));
            }
            else
            {
               work[0] = ZERO;
               work[1] = ZERO;
            }
/*
*  Combine the local results using a 1-tree topology within process row 0
*  if nprow > 1 or Xrow >= 0, i.e sub( X ) is distributed.
*/
            if( ( nprow >= 2 ) && ( Xrow >= 0 ) )
            {
               mydist = myrow;
               k      = 1;
l_30:
               if( mydist & 1 )
               {
                  dist = k * ( mydist - 1 );
                  dst  = MPosMod( dist, nprow );
                  Cdgesd2d( ctxt, 2, 1, ((char*)work), 2, dst, mycol );
                  goto l_40;
               }
               else
               {
                  dist = myrow + k;
                  src  = MPosMod( dist, nprow );

                  if( myrow < src )
                  {
                     Cdgerv2d( ctxt, 2, 1, ((char*) &work[2]), 2,
                               src, mycol );
                     if( ABS( work[0] ) < ABS( work[2] ) )
                     { work[0] = work[2]; work[1] = work[3]; }
                  }
                  mydist >>= 1;
               }
               k <<= 1;

               if( k < nprow ) goto l_30;
l_40:
/*
*  Process row 0 broadcasts the combined values of INDX and AMAX within their
*  process column.
*/
               cbtop = *PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
               if( myrow == 0 )
               {
                  Cdgebs2d( ctxt, COLUMN, &cbtop, 2, 1, ((char*)work), 2 );
               }
               else
               {
                  Cdgebr2d( ctxt, COLUMN, &cbtop, 2, 1, ((char*)work), 2,
                            0, mycol );
               }
            }
/*
*  Set INDX and AMAX to the replicated answers contained in work. If AMAX is
*  zero, then select a coherent INDX.
*/
            *AMAX = work[0];
            *INDX = ( ( *AMAX == ZERO ) ? ( *IX ) : ( (int)(work[1]) ) );
         }
         else
         {
/*
*  Otherwise use the current topology settings to combine the results
*/
            Ximb = Xd[IMB_ ]; Xmb = Xd[MB_ ]; Xsrc = Xd[RSRC_];
            Xnp = PB_Cnumroc( *N, Xi, Ximb, Xmb, myrow, Xsrc, nprow );
/*
*  Make sure I own some data and compute local INDX and AMAX
*/

            if( Xnp > 0 )
            {
/*
*  Compute the local maximum and its corresponding local index
*/
               Xld = Xd[LLD_];
               Xlindx = Xii - 1 +
                        idamax_( &Xnp, ((char*)Mptr(X,Xii,Xjj,Xld,1)), INCX );
               *AMAX = *Mptr(X,Xlindx,Xjj,Xld,1);
            }
            else
            {
               *AMAX = ZERO;
            }

            if( Xrow >= 0 )
            {
/*
*  Combine leave on all the local maximum if Xrow >= 0, i.e sub( X ) is
*  distributed.
*/
               Cdgamx2d( ctxt, COLUMN, &cctop, 1, 1, ((char*)AMAX), 1,
                         &maxpos, &idumm, 1, -1, mycol );
/*
*  Broadcast the corresponding global index
*/
               if( *AMAX != ZERO )
               {
                  cbtop = *PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
                  if( myrow == maxpos )
                  {
                     Mindxl2g( Xgindx, Xlindx, Ximb, Xmb, myrow, Xsrc, nprow );
                     *INDX = Xgindx + 1;
                     Cigebs2d( ctxt, COLUMN, &cbtop, 1, 1, ((char*)INDX), 1 );
                  }
                  else
                  {
                     Cigebr2d( ctxt, COLUMN, &cbtop, 1, 1, ((char*)INDX), 1,
                               maxpos, mycol );
                  }
               }
               else
               {
/*
*  If AMAX is zero, then select a coherent INDX.
*/
                  *INDX = *IX;
               }
            }
            else
            {
/*
*  sub( X ) is not distributed. If AMAX is zero, then select a coherent INDX.
*/
               *INDX = ( ( *AMAX == ZERO ) ? ( *IX ) : Xlindx + 1 );
            }
         }
      }
      return;
   }
/*
*  End of PDAMAX
*/
}
