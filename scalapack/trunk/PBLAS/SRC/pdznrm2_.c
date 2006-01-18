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
void pdznrm2_( int * N, double * NORM2,
               double * X, int * IX, int * JX, int * DESCX, int * INCX )
#else
void pdznrm2_( N, NORM2, X, IX, JX, DESCX, INCX )
/*
*  .. Scalar Arguments ..
*/
   int            * INCX, * IX, * JX, * N;
   double         * NORM2;
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
*  PDZNRM2  computes the 2-norm of a subvector sub( X ),
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
*  NORM2   (local output) DOUBLE PRECISION
*          On exit, NORM2 specifies the 2-norm of the subvector sub( X )
*          only in its scope (See below for further details).
*
*  X       (local input) COMPLEX*16 array
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
   char           * Xptr = NULL, top;
   int            Xcol, Xi, Xii, Xj, Xjj, Xld, Xnp, Xnq, Xrow, ctxt, dst, dist,
                  info, k, mycol, mydist, myrow, npcol, nprow, src, size;
   double         Xtmp, scale, ssq, temp1, temp2;
   PBTYP_T        * type;
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
   if( !( info = ( ( nprow == -1 ) ? -( 601 + CTXT_ ) : 0 ) ) )
      PB_Cchkvec( ctxt, "PDZNRM2", "X", *N, 1, Xi, Xj, Xd, *INCX, 6, &info );
   if( info ) { PB_Cabort( ctxt, "PDZNRM2", info ); return; }
#endif
/*
*  Initialize NORM2
*/
   *NORM2 = ZERO;
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
   if( ( *N == 1 ) && ( *INCX == 1 ) && ( Xd[M_] == 1 ) )
   {
/*
*  Make sure I own some data and compute NORM2
*/
      if( ( ( myrow == Xrow ) || ( Xrow < 0 ) ) &&
          ( ( mycol == Xcol ) || ( Xcol < 0 ) ) )
      {
         scale = ZERO;
         ssq   = ONE;
         type  = PB_Cztypeset();
         Xptr  = Mptr( ((char *) X), Xii, Xjj, Xd[LLD_], type->size );
         Xtmp  = ((double *) Xptr)[REAL_PART];
         if( Xtmp != ZERO )
         {
            temp1 = ABS( Xtmp );
            if( scale < temp1 )
            {
               temp2 = scale / temp1;
               ssq   = ONE + ssq * ( temp2 * temp2 );
               scale = temp1;
            }
            else
            {
               temp2 = temp1 / scale;
               ssq   = ssq +       ( temp2 * temp2 );
            }
         }
         Xtmp  = ((double *) Xptr)[IMAG_PART];
         if( Xtmp != ZERO )
         {
            temp1 = ABS( Xtmp );
            if( scale < temp1 )
            {
               temp2 = scale / temp1;
               ssq   = ONE + ssq * ( temp2 * temp2 );
               scale = temp1;
            }
            else
            {
               temp2 = temp1 / scale;
               ssq   = ssq +       ( temp2 * temp2 );
            }
         }
/*
*  Compute NORM2 = SCALE * SQRT( SSQ )
*/
         dasqrtb_( &scale, &ssq, NORM2 );
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
/*
*        Initialize SCALE and SSQ
*/
         scale = ZERO;
         ssq   = ONE;
/*
*  Make sure I own some data and compute local sum of squares
*/
         Xnq = PB_Cnumroc( *N, Xj, Xd[INB_], Xd[NB_], mycol, Xd[CSRC_], npcol );
         if( Xnq > 0 )
         {
            Xld  = Xd[LLD_];
            type = PB_Cztypeset(); size = type->size;
            Xptr = Mptr( ((char *) X), Xii, Xjj, Xld, size );

            for( k = 0; k < Xnq; k++ )
            {
               Xtmp = ((double *) Xptr)[REAL_PART];
               if( Xtmp != ZERO )
               {
                  temp1 = ABS( Xtmp );
                  if( scale < temp1 )
                  {
                     temp2 = scale / temp1;
                     ssq   = ONE + ssq * ( temp2 * temp2 );
                     scale = temp1;
                  }
                  else
                  {
                     temp2 = temp1 / scale;
                     ssq   = ssq +       ( temp2 * temp2 );
                  }
               }
               Xtmp = ((double *) Xptr)[IMAG_PART];
               if( Xtmp != ZERO )
               {
                  temp1 = ABS( Xtmp );
                  if( scale < temp1 )
                  {
                     temp2 = scale / temp1;
                     ssq   = ONE + ssq * ( temp2 * temp2 );
                     scale = temp1;
                  }
                  else
                  {
                     temp2 = temp1 / scale;
                     ssq   = ssq +       ( temp2 * temp2 );
                  }
               }
               Xptr += Xld * size;
            }
         }
/*
*  If Xnq <= 0, SCALE is zero and SSQ is one (see initialization above)
*/
         if( ( npcol >= 2 ) && ( Xcol >= 0 ) )
         {
/*
*  Combine the local sum of squares using a 1-tree topology within process row
*  0 if npcol > 1 and Xcol >= 0, i.e sub( X ) is distributed.
*/
            work[0] = scale;
            work[1] = ssq;

            mydist  = mycol;
            k       = 1;
l_10:
            if( mydist & 1 )
            {
               dist = k * ( mydist - 1 );
               dst  = MPosMod( dist, npcol );
               Cdgesd2d( ctxt, 2, 1, ((char*) work), 2, myrow, dst );
               goto l_20;
            }
            else
            {
               dist = mycol + k;
               src  = MPosMod( dist, npcol );

               if( mycol < src )
               {
                  Cdgerv2d( ctxt, 2, 1, ((char*)&work[2]), 2, myrow, src );
                  if( work[0] >= work[2] )
                  {
                     if( work[0] != ZERO )
                     {
                        temp1   = work[2] / work[0];
                        work[1] = work[1] + ( temp1 * temp1 ) * work[3];
                     }
                  }
                  else
                  {
                     temp1   = work[0] / work[2];
                     work[1] = work[3] + ( temp1 * temp1 ) * work[1];
                     work[0] = work[2];
                  }
               }
               mydist >>= 1;
            }
            k <<= 1;

            if( k < npcol ) goto l_10;
l_20:
/*
*  Process column 0 broadcasts the combined values of SCALE and SSQ within their
*  process row.
*/
            top = *PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );
            if( mycol == 0 )
            {
               Cdgebs2d( ctxt, ROW, &top, 2, 1, ((char*)work), 2 );
            }
            else
            {
               Cdgebr2d( ctxt, ROW, &top, 2, 1, ((char*)work), 2,
                         myrow, 0 );
            }
/*
*  Compute NORM2 redundantly NORM2  = WORK( 1 ) * SQRT( WORK( 2 ) )
*/
            dasqrtb_( &work[0], &work[1], NORM2 );
         }
         else
         {
/*
*  Compute NORM2 redundantly ( sub( X ) is not distributed )
*/
            dasqrtb_( &scale, &ssq, NORM2 );
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
/*
*  Initialize SCALE and SSQ
*/
         scale = ZERO;
         ssq   = ONE;
/*
*  Make sure I own some data and compute local sum of squares
*/
         Xnp = PB_Cnumroc( *N, Xi, Xd[IMB_], Xd[MB_], myrow, Xd[RSRC_], nprow );
         if( Xnp > 0 )
         {
            type = PB_Cztypeset(); size = type->size;
            Xptr = Mptr( ((char *) X), Xii, Xjj, Xd[LLD_], size );

            for( k = 0; k < Xnp; k++ )
            {
               Xtmp = ((double *) Xptr)[REAL_PART];
               if( Xtmp != ZERO )
               {
                  temp1 = ABS( Xtmp );
                  if( scale < temp1 )
                  {
                     temp2 = scale / temp1;
                     ssq   = ONE + ssq * ( temp2 * temp2 );
                     scale = temp1;
                  }
                  else
                  {
                     temp2 = temp1 / scale;
                     ssq   = ssq +       ( temp2 * temp2 );
                  }
               }
               Xtmp = ((double *) Xptr)[IMAG_PART];
               if( Xtmp != ZERO )
               {
                  temp1 = ABS( Xtmp );
                  if( scale < temp1 )
                  {
                     temp2 = scale / temp1;
                     ssq   = ONE + ssq * ( temp2 * temp2 );
                     scale = temp1;
                  }
                  else
                  {
                     temp2 = temp1 / scale;
                     ssq   = ssq +       ( temp2 * temp2 );
                  }
               }
               Xptr += size;
            }
         }
/*
*  If Xnp <= 0, SCALE is zero and SSQ is one (see initialization above)
*/
         if( ( nprow >= 2 ) && ( Xrow >= 0 ) )
         {
/*
*  Combine the local sum of squares using a 1-tree topology within process
*  column 0 if nprow > 1 and Xrow >= 0, i.e sub( X ) is distributed.
*/
            work[0] = scale;
            work[1] = ssq;

            mydist  = myrow;
            k       = 1;
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
                  Cdgerv2d( ctxt, 2, 1, ((char*)&work[2]), 2, src, mycol );
                  if( work[0] >= work[2] )
                  {
                     if( work[0] != ZERO )
                     {
                        temp1   = work[2] / work[0];
                        work[1] = work[1] + ( temp1 * temp1 ) * work[3];
                     }
                  }
                  else
                  {
                     temp1   = work[0] / work[2];
                     work[1] = work[3] + ( temp1 * temp1 ) * work[1];
                     work[0] = work[2];
                  }
               }
               mydist >>= 1;
            }
            k <<= 1;

            if( k < nprow ) goto l_30;
l_40:
/*
*  Process column 0 broadcasts the combined values of SCALE and SSQ within their
*  process column
*/
            top = *PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
            if( myrow == 0 )
            {
               Cdgebs2d( ctxt, COLUMN, &top, 2, 1, ((char*)work), 2 );
            }
            else
            {
               Cdgebr2d( ctxt, COLUMN, &top, 2, 1, ((char*)work), 2,
                         0, mycol );
            }
/*
*  Compute NORM2 redundantly NORM2 = WORK[0] * SQRT( WORK[1] )
*/
            dasqrtb_( &work[0], &work[1], NORM2 );
         }
         else
         {
/*
*  Compute NORM2 redundantly ( sub( X ) is not distributed )
*/
            dasqrtb_( &scale, &ssq, NORM2 );
         }
      }
      return;
   }
/*
*  End of PDZNRM2
*/
}
