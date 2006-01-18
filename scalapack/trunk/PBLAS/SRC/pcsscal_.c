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
void pcsscal_( int * N,
               float * ALPHA,
               float * X, int * IX, int * JX, int * DESCX, int * INCX )
#else
void pcsscal_( N, ALPHA, X, IX, JX, DESCX, INCX )
/*
*  .. Scalar Arguments ..
*/
   int            * INCX, * IX, * JX, * N;
   float          * ALPHA;
/*
*  .. Array Arguments ..
*/
   int            * DESCX;
   float          * X;
#endif
{
/*
*  Purpose
*  =======
*
*  PCSSCAL multiplies an n element subvector sub( X ) by the real scalar
*  alpha,
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
*  ALPHA   (global input) REAL
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied as zero then the local entries of the array  X  cor-
*          responding to the entries of the subvector sub( X ) need  not
*          be set on input.
*
*  X       (local input/local output) COMPLEX array
*          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
*          is   at  least  MAX( 1, Lr( 1, IX ) )  when  INCX = M_X   and
*          MAX( 1, Lr( 1, IX+N-1 ) )  otherwise,  and,  Kx  is  at least
*          Lc( 1, JX+N-1 )  when  INCX = M_X  and Lc( 1, JX ) otherwise.
*          Before  entry,  this array  contains the local entries of the
*          matrix  X.  On exit, sub( X ) is overwritten with the  scaled
*          subvector.
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
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   int            Xcol, Xi, Xii, Xj, Xjj, Xld, Xnp, Xnq, Xrow, ctxt, info,
                  mycol, myrow, npcol, nprow;
   PBTYP_T        * type;
/*
*  .. Local Arrays ..
*/
   int            Xd[DLEN_];
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
      PB_Cchkvec( ctxt, "PCSSCAL", "X", *N, 1, Xi, Xj, Xd, *INCX, 6, &info );
   if( info ) { PB_Cabort( ctxt, "PCSSCAL", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( ( *N == 0 ) || ( ALPHA[REAL_PART] == ONE ) ) return;
/*
*  Retrieve process grid information
*/
#ifdef NO_ARGCHK
   Cblacs_gridinfo( Xd[CTXT_], &nprow, &npcol, &myrow, &mycol );
#endif
/*
*  Retrieve sub( X )'s local information: Xii, Xjj, Xrow, Xcol
*/
   PB_Cinfog2l( Xi, Xj, Xd, nprow, npcol, myrow, mycol, &Xii, &Xjj, &Xrow,
                &Xcol );
/*
*  Start the operations
*/
   if( *INCX == Xd[M_] )
   {
/*
*  sub( X ) resides in (a) process row(s)
*/
      if( ( myrow == Xrow ) || ( Xrow < 0 ) )
      {
/*
*  Make sure I own some data and scale sub( X )
*/
         Xnq = PB_Cnumroc( *N, Xj, Xd[INB_], Xd[NB_], mycol, Xd[CSRC_], npcol );
         if( Xnq > 0 )
         {
            Xld = Xd[LLD_];
            type = PB_Cctypeset();
            if( ALPHA[REAL_PART] == ZERO )
            {
               cset_( &Xnq, type->zero, Mptr( ((char *) X), Xii, Xjj, Xld,
                      type->size ), &Xld );
            }
            else
            {
               csscal_( &Xnq, ((char *) ALPHA), Mptr( ((char *) X),
                        Xii, Xjj, Xld, type->size ), &Xld );
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
/*
*  Make sure I own some data and scale sub( X )
*/
         Xnp = PB_Cnumroc( *N, Xi, Xd[IMB_], Xd[MB_], myrow, Xd[RSRC_], nprow );
         if( Xnp > 0 )
         {
            type = PB_Cctypeset();
            if( ALPHA[REAL_PART] == ZERO )
            {
               cset_( &Xnp, type->zero, Mptr( ((char *) X), Xii, Xjj,
                      Xd[LLD_], type->size ), INCX );
            }
            else
            {
               csscal_( &Xnp, ((char *) ALPHA), Mptr( ((char *) X),
                        Xii, Xjj, Xd[LLD_], type->size ), INCX );
            }
         }
      }
      return;
   }
/*
*  End of PCSSCAL
*/
}
