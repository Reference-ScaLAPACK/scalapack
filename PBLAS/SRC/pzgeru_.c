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
void pzgeru_( int * M, int * N, double * ALPHA,
              double * X, int * IX, int * JX, int * DESCX, int * INCX,
              double * Y, int * IY, int * JY, int * DESCY, int * INCY,
              double * A, int * IA, int * JA, int * DESCA )
#else
void pzgeru_( M, N, ALPHA, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY,
              INCY, A, IA, JA, DESCA )
/*
*  .. Scalar Arguments ..
*/
   int            * IA, * INCX, * INCY, * IX, * IY, * JA, * JX, * JY,
                  * M, * N;
   double         * ALPHA;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCX, * DESCY;
   double         * A, * X, * Y;
#endif
{
/*
*  Purpose
*  =======
*
*  PZGERU  performs the rank 1 operation
*
*     sub( A ) := alpha*sub( X )*sub( Y )' + sub( A ),
*
*  where
*
*     sub( A ) denotes A(IA:IA+M-1,JA:JA+N-1),
*
*     sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
*                      X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X, and,
*
*     sub( Y ) denotes Y(IY,JY:JY+N-1) if INCY = M_Y,
*                      Y(IY:IY+N-1,JY) if INCY = 1 and INCY <> M_Y.
*
*  Alpha  is a scalar, sub( X )  is an m element subvector,  sub( Y ) is
*  an n element subvector and sub( A ) is an m by n submatrix.
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
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( A ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( A ). N  must be at least zero.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of the arrays  X
*          and Y corresponding to the entries of the subvectors sub( X )
*          and sub( Y ) respectively need not be set on input.
*
*  X       (local input) COMPLEX*16 array
*          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
*          is   at  least  MAX( 1, Lr( 1, IX ) )  when  INCX = M_X   and
*          MAX( 1, Lr( 1, IX+M-1 ) )  otherwise,  and,  Kx  is  at least
*          Lc( 1, JX+M-1 )  when  INCX = M_X  and Lc( 1, JX ) otherwise.
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
*  A       (local input/local output) COMPLEX*16 array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix A.
*          On exit, the entries of this array corresponding to the local
*          entries  of  the  submatrix  sub( A )  are overwritten by the
*          local entries of the m by n updated submatrix.
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
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   int            Acol, Ai, Aii, Aimb1, Ainb1, Aj, Ajj, Ald, Amb, Amp, Anb,
                  Anq, Arow, XAfr, Xi, Xj, YAfr, Yi, Yj, ctxt, info, ione=1,
                  mycol, myrow, npcol, nprow;
   PBTYP_T        * type;
/*
*  .. Local Arrays ..
*/
   int            Ad[DLEN_], Ad0[DLEN_], XAd[DLEN_], Xd[DLEN_], YAd[DLEN_],
                  Yd[DLEN_];
   char           * XA = NULL, * YA = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   PB_CargFtoC( *IA, *JA, DESCA, &Ai, &Aj, Ad );
   PB_CargFtoC( *IX, *JX, DESCX, &Xi, &Xj, Xd );
   PB_CargFtoC( *IY, *JY, DESCY, &Yi, &Yj, Yd );
#ifndef NO_ARGCHK
/*
*  Test the input parameters
*/
   Cblacs_gridinfo( ( ctxt = Xd[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
   if( !( info = ( ( nprow == -1 ) ? -( 701 + CTXT_ ) : 0 ) ) )
   {
      PB_Cchkvec( ctxt, "PZGERU", "X", *M, 1, Xi, Xj, Xd, *INCX,  7, &info );
      PB_Cchkvec( ctxt, "PZGERU", "Y", *N, 2, Yi, Yj, Yd, *INCY, 12, &info );
      PB_Cchkmat( ctxt, "PZGERU", "A", *M, 1, *N, 2, Ai, Aj, Ad, 17, &info );
   }
   if( info ) { PB_Cabort( ctxt, "PZGERU", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( ( *M == 0 ) || ( *N == 0 ) ||
       ( ( ALPHA[REAL_PART] == ZERO ) && ( ALPHA[IMAG_PART] == ZERO ) ) )
      return;
/*
*  Retrieve process grid information
*/
#ifdef NO_ARGCHK
   Cblacs_gridinfo( ( ctxt = Ad[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
#endif
/*
*  Get type structure
*/
   type = PB_Cztypeset();
/*
*  Compute descriptor Ad0 for sub( A )
*/
   PB_Cdescribe( *M, *N, Ai, Aj, Ad, nprow, npcol, myrow, mycol, &Aii, &Ajj,
                 &Ald, &Aimb1, &Ainb1, &Amb, &Anb, &Arow, &Acol, Ad0 );

/*
*  Replicate sub( X ) in process columns spanned by sub( A ) -> XA
*/
   PB_CInV( type, NOCONJG, COLUMN, *M, *N, Ad0, 1, ((char *) X), Xi, Xj, Xd,
            ( *INCX == Xd[M_] ? ROW : COLUMN ), &XA, XAd, &XAfr );
/*
*  Replicate sub( Y ) in process rows spanned by sub( A ) -> YA
*/
   PB_CInV( type, NOCONJG, ROW,    *M, *N, Ad0, 1, ((char *) Y), Yi, Yj, Yd,
            ( *INCY == Yd[M_] ? ROW : COLUMN ), &YA, YAd, &YAfr );
/*
*  Local rank-1 update iff I own some data
*/
   Amp = PB_Cnumroc( *M, 0, Aimb1, Amb, myrow, Arow, nprow );
   Anq = PB_Cnumroc( *N, 0, Ainb1, Anb, mycol, Acol, npcol );

   if( ( Amp > 0 ) && ( Anq > 0 ) )
   {
      zgeru_( &Amp, &Anq, ((char *) ALPHA), XA, &ione, YA, &YAd[LLD_],
              Mptr( ((char *) A), Aii, Ajj, Ald, type->size ), &Ald );
   }
   if( XAfr ) free( XA );
   if( YAfr ) free( YA );
/*
*  End of PZGERU
*/
}
