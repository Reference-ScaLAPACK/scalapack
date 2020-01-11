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
void pdagemv_( F_CHAR_T TRANS, Int * M, Int * N, double * ALPHA,
               double * A, Int * IA, Int * JA, Int * DESCA,
               double * X, Int * IX, Int * JX, Int * DESCX, Int * INCX,
               double * BETA,
               double * Y, Int * IY, Int * JY, Int * DESCY, Int * INCY )
#else
void pdagemv_( TRANS, M, N, ALPHA, A, IA, JA, DESCA, X, IX, JX, DESCX,
               INCX, BETA, Y, IY, JY, DESCY, INCY )
/*
*  .. Scalar Arguments ..
*/
   F_CHAR_T       TRANS;
   Int            * IA, * INCX, * INCY, * IX, * IY, * JA, * JX, * JY,
                  * M, * N;
   double         * ALPHA, * BETA;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA, * DESCX, * DESCY;
   double         * A, * X, * Y;
#endif
{
/*
*  Purpose
*  =======
*
*  PDAGEMV  performs one of the matrix-vector operations
*
*     sub( Y ) := abs( alpha )*abs( sub( A ) )*abs( sub( X ) ) +
*                 abs( beta*sub( Y ) ),
*     or
*
*     sub( Y ) := abs( alpha )*abs( sub( A )' )*abs( sub( X ) ) +
*                 abs( beta*sub( Y ) ),
*
*  where
*
*     sub( A ) denotes A(IA:IA+M-1,JA:JA+N-1).
*
*  When TRANS = 'N',
*
*     sub( X ) denotes X(IX:IX,JX:JX+N-1), if INCX = M_X,
*                      X(IX:IX+N-1,JX:JX), if INCX = 1 and INCX <> M_X,
*     and,
*
*     sub( Y ) denotes Y(IY:IY,JY:JY+M-1), if INCY = M_Y,
*                      Y(IY:IY+M-1,JY:JY), if INCY = 1 and INCY <> M_Y,
*  and, otherwise
*
*     sub( X ) denotes X(IX:IX,JX:JX+M-1), if INCX = M_X,
*                      X(IX:IX+M-1,JX:JX), if INCX = 1 and INCX <> M_X,
*     and,
*
*     sub( Y ) denotes Y(IY:IY,JY:JY+N-1), if INCY = M_Y,
*                      Y(IY:IY+N-1,JY:JY), if INCY = 1 and INCY <> M_Y.
*
*  Alpha  and  beta  are  real  scalars,  sub( Y )  is a real subvector,
*  sub( X ) is a subvector and sub( A ) is an m by n submatrix.
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
*  TRANS   (global input) CHARACTER*1
*          On entry,  TRANS  specifies the  operation to be performed as
*          follows:
*
*             TRANS = 'N' or 'n'
*                sub( Y ) := |alpha|*|sub( A ) | * |sub( X )| +
*                                                       |beta*sub( Y )|,
*
*             TRANS = 'T' or 't',
*                sub( Y ) := |alpha|*|sub( A )'| * |sub( X )| +
*                                                       |beta*sub( Y )|,
*
*             TRANS = 'C' or 'c',
*                sub( Y ) := |alpha|*|sub( A )'| * |sub( X )| +
*                                                       |beta*sub( Y )|.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( A ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( A ). N  must be at least zero.
*
*  ALPHA   (global input) DOUBLE PRECISION
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of the arrays  A
*          and X corresponding to the entries of the submatrix  sub( A )
*          and the subvector sub( X ) need not be set on input.
*
*  A       (local input) DOUBLE PRECISION array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ). Before  entry,  this array contains
*          the local entries of the matrix A.
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
*  X       (local input) DOUBLE PRECISION array
*          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
*          is  at  least  MAX( 1, Lr( 1, IX ) )   when  INCX = M_X   and
*          MAX( 1, Lr( 1, IX+Lx-1 ) )  otherwise,  and,  Kx  is at least
*          Lc( 1, JX+Lx-1 ) when  INCX = M_X  and Lc( 1, JX ) otherwise.
*          Lx is N when TRANS = 'N' or 'n' and  M  otherwise. Before en-
*          try, this array  contains the local entries of the matrix X.
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
*  BETA    (global input) DOUBLE PRECISION
*          On entry,  BETA  specifies the scalar  beta.   When  BETA  is
*          supplied  as  zero  then  the  local entries of  the array  Y
*          corresponding to the entries of the subvector  sub( Y )  need
*          not be set on input.
*
*  Y       (local input/local output) DOUBLE PRECISION array
*          On entry, Y is an array of dimension (LLD_Y, Ky), where LLD_Y
*          is   at  least  MAX( 1, Lr( 1, IY ) )  when  INCY = M_Y   and
*          MAX( 1, Lr( 1, IY+Ly-1 ) )  otherwise, and,  Ky  is  at least
*          Lc( 1, JY+Ly-1 ) when  INCY = M_Y  and Lc( 1, JY ) otherwise.
*          Ly is M when TRANS = 'N' or 'n' and  N  otherwise. Before en-
*          try, this  array  contains the local entries of the matrix Y.
*          On exit, sub( Y ) is overwritten by the updated subvector.
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
   char           TrA, Yroc, * one, * tbeta, top;
   Int            Acol, Ai, Aii, Aimb1, Ainb1, Aj, Ajj, Ald, Amb, Amp, Anb,
                  Anq, Arow, XAfr, Xi, Xj, YAfr, YApbY, YAsum, Ycol, Yi, Yii,
                  Yj, Yjj, Yld, Ynp, Ynq, Yrow, ctxt, info, ione=1, mycol,
                  myrow, nota, npcol, nprow;
   PBTYP_T        * type, * utyp;
/*
*  .. Local Arrays ..
*/
   Int            Ad [DLEN_], Ad0[DLEN_], XAd[DLEN_], Xd[DLEN_], YAd[DLEN_],
                  Yd [DLEN_];
   char           * XA = NULL, * YA = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   nota = ( ( TrA = Mupcase( F2C_CHAR( TRANS )[0] ) ) == CNOTRAN );
   PB_CargFtoC( *IA, *JA, DESCA, &Ai, &Aj, Ad );
   PB_CargFtoC( *IX, *JX, DESCX, &Xi, &Xj, Xd );
   PB_CargFtoC( *IY, *JY, DESCY, &Yi, &Yj, Yd );
#ifndef NO_ARGCHK
/*
*  Test the input parameters
*/
   Cblacs_gridinfo( ( ctxt = Ad[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
   if( !( info = ( ( nprow == -1 ) ? -( 801 + CTXT_ ) : 0 ) ) )
   {
      if( ( !nota ) && ( TrA != CTRAN ) && ( TrA != CCOTRAN ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PDAGEMV", "Illegal TRANS=%c\n", TrA );
         info = -1;
      }
      PB_Cchkmat(    ctxt, "PDAGEMV", "A", *M, 2, *N, 3, Ai, Aj, Ad,  8,
                     &info );
      if( nota )
      {
         PB_Cchkvec( ctxt, "PDAGEMV", "X", *N, 3, Xi, Xj, Xd, *INCX, 12,
                     &info );
         PB_Cchkvec( ctxt, "PDAGEMV", "Y", *M, 2, Yi, Yj, Yd, *INCY, 18,
                     &info );
      }
      else
      {
         PB_Cchkvec( ctxt, "PDAGEMV", "X", *M, 2, Xi, Xj, Xd, *INCX, 12,
                     &info );
         PB_Cchkvec( ctxt, "PDAGEMV", "Y", *N, 3, Yi, Yj, Yd, *INCY, 18,
                     &info );
      }
   }
   if( info ) { PB_Cabort( ctxt, "PDAGEMV", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( ( *M == 0 ) || ( *N == 0 ) ||
       ( ( ALPHA[REAL_PART] == ZERO ) && ( BETA[REAL_PART] == ONE ) ) )
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
   type = utyp = PB_Cdtypeset();
/*
*  When alpha is zero
*/
   if( ALPHA[REAL_PART] == ZERO )
   {
/*
*  Retrieve sub( Y )'s local information: Yii, Yjj, Yrow, Ycol
*/
      PB_Cinfog2l( Yi, Yj, Yd, nprow, npcol, myrow, mycol, &Yii, &Yjj,
                   &Yrow, &Ycol );

      if( *INCY == Yd[M_] )
      {
/*
*  sub( Y ) resides in (a) process row(s)
*/
         if( ( myrow == Yrow ) || ( Yrow < 0 ) )
         {
/*
*  Make sure I own some data and scale sub( Y )
*/
            Ynq = PB_Cnumroc( ( nota ? *M : *N ), Yj, Yd[INB_], Yd[NB_], mycol,
                              Yd[CSRC_], npcol );
            if( Ynq > 0 )
            {
               Yld = Yd[LLD_];
               dascal_( &Ynq, ((char *) BETA), Mptr( ((char *) Y), Yii,
                        Yjj, Yld, utyp->size ), &Yld );
            }
         }
      }
      else
      {
/*
*  sub( Y ) resides in (a) process column(s)
*/
         if( ( mycol == Ycol ) || ( Ycol < 0 ) )
         {
/*
*  Make sure I own some data and scale sub( Y )
*/
            Ynp = PB_Cnumroc( ( nota ? *M : *N ), Yi, Yd[IMB_], Yd[MB_], myrow,
                              Yd[RSRC_], nprow );
            if( Ynp > 0 )
            {
               dascal_( &Ynp, ((char *) BETA), Mptr( ((char *) Y), Yii,
                        Yjj, Yd[LLD_], utyp->size ), INCY );
            }
         }
      }
      return;
   }
/*
*  Compute descriptor Ad0 for sub( A )
*/
   PB_Cdescribe( *M, *N, Ai, Aj, Ad, nprow, npcol, myrow, mycol, &Aii, &Ajj,
                 &Ald, &Aimb1, &Ainb1, &Amb, &Anb, &Arow, &Acol, Ad0 );

   Yroc = ( *INCY == Yd[M_] ? CROW : CCOLUMN );

   if( nota )
   {
/*
*  Reuse sub( Y ) and/or create vector YA in process columns spanned by sub( A )
*/
      PB_CInOutV( utyp, COLUMN, *M, *N, Ad0, 1, ((char *) BETA), ((char *) Y),
                  Yi, Yj, Yd, &Yroc, &tbeta, &YA, YAd, &YAfr, &YAsum, &YApbY );
/*
*  Replicate sub( X ) in process rows spanned by sub( A ) -> XA
*/
      PB_CInV( type, NOCONJG, ROW,    *M, *N, Ad0, 1, ((char *) X), Xi, Xj, Xd,
               ( *INCX == Xd[M_] ? ROW : COLUMN ), &XA, XAd, &XAfr );
/*
*  Local matrix-vector multiply iff I own some data
*/
      Amp = PB_Cnumroc( *M, 0, Ad0[IMB_], Ad0[MB_], myrow, Ad0[RSRC_], nprow );
      Anq = PB_Cnumroc( *N, 0, Ad0[INB_], Ad0[NB_], mycol, Ad0[CSRC_], npcol );
      if( ( Amp > 0 ) && ( Anq > 0 ) )
      {
         dagemv_( TRANS, &Amp, &Anq, ((char *) ALPHA), Mptr( ((char *) A),
                  Aii, Ajj, Ald, type->size), &Ald, XA, &XAd[LLD_], tbeta,
                  YA, &ione );
      }
      if( XAfr ) free( XA );
/*
*  Combine the partial column results into YA
*/
      if( YAsum && ( Amp > 0 ) )
      {
         top = *PB_Ctop( &ctxt, COMBINE, ROW, TOP_GET );
         Cdgsum2d( ctxt, ROW, &top, Amp, 1, YA, YAd[LLD_], myrow,
                   YAd[CSRC_] );
      }
   }
   else
   {
/*
*  Reuse sub( Y ) and/or create vector YA in process rows spanned by sub( A )
*/
      PB_CInOutV( utyp, ROW, *M, *N, Ad0, 1, ((char *) BETA), ((char *) Y), Yi,
                  Yj, Yd, &Yroc, &tbeta, &YA, YAd, &YAfr, &YAsum, &YApbY );
/*
*  Replicate sub( X ) in process columns spanned by sub( A ) -> XA
*/
      PB_CInV( type, NOCONJG, COLUMN, *M, *N, Ad0, 1, ((char *) X), Xi, Xj, Xd,
               ( *INCX == Xd[M_] ? ROW : COLUMN ), &XA, XAd, &XAfr );
/*
*  Local matrix-vector multiply iff I own some data
*/
      Amp = PB_Cnumroc( *M, 0, Ad0[IMB_], Ad0[MB_], myrow, Ad0[RSRC_], nprow );
      Anq = PB_Cnumroc( *N, 0, Ad0[INB_], Ad0[NB_], mycol, Ad0[CSRC_], npcol );
      if( ( Amp > 0 ) && ( Anq > 0 ) )
      {
         dagemv_( TRANS, &Amp, &Anq, ((char *) ALPHA), Mptr( ((char *) A),
                  Aii, Ajj, Ald, type->size ), &Ald, XA, &ione, tbeta, YA,
                  &YAd[LLD_] );
      }
      if( XAfr ) free( XA );
/*
*  Combine the partial row results into YA
*/
      if( YAsum && ( Anq > 0 ) )
      {
         top = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
         Cdgsum2d( ctxt, COLUMN, &top, 1, Anq, YA, YAd[LLD_], YAd[RSRC_],
                   mycol );
      }
   }
/*
*  sub( Y ) := beta * sub( Y ) + YA (if necessary)
*/
   if( YApbY )
   {
/*
*  Retrieve sub( Y )'s local information: Yii, Yjj, Yrow, Ycol
*/
      PB_Cinfog2l( Yi, Yj, Yd, nprow, npcol, myrow, mycol, &Yii, &Yjj, &Yrow,
                   &Ycol );

      if( *INCY == Yd[M_] )
      {
/*
*  sub( Y ) resides in (a) process row(s)
*/
         if( ( myrow == Yrow ) || ( Yrow < 0 ) )
         {
/*
*  Make sure I own some data and scale sub( Y )
*/
            Ynq = PB_Cnumroc( ( nota ? *M : *N ), Yj, Yd[INB_], Yd[NB_], mycol,
                              Yd[CSRC_], npcol );
            if( Ynq > 0 )
            {
               Yld = Yd[LLD_];
               dascal_( &Ynq, ((char *) BETA), Mptr( ((char *) Y), Yii,
                        Yjj, Yld, utyp->size ), &Yld );
            }
         }
      }
      else
      {
/*
*  sub( Y ) resides in (a) process column(s)
*/
         if( ( mycol == Ycol ) || ( Ycol < 0 ) )
         {
/*
*  Make sure I own some data and scale sub( Y )
*/
            Ynp = PB_Cnumroc( ( nota ? *M : *N ), Yi, Yd[IMB_], Yd[MB_], myrow,
                              Yd[RSRC_], nprow );
            if( Ynp > 0 )
            {
               dascal_( &Ynp, ((char *) BETA), Mptr( ((char *) Y), Yii,
                        Yjj, Yd[LLD_], utyp->size ), INCY );
            }
         }
      }

      one = utyp->one;

      if( nota )
      {
         PB_Cpaxpby( utyp, NOCONJG, *M, 1, one, YA, 0, 0, YAd, COLUMN, one,
                     ((char *) Y), Yi, Yj, Yd, &Yroc );
      }
      else
      {
         PB_Cpaxpby( utyp, NOCONJG, 1, *N, one, YA, 0, 0, YAd, ROW,    one,
                     ((char *) Y), Yi, Yj, Yd, &Yroc );
      }
   }
   if( YAfr ) free( YA );
/*
*  End of PDAGEMV
*/
}
