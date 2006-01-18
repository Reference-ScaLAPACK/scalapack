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
void pdsymv_( F_CHAR_T UPLO, int * N, double * ALPHA,
              double * A, int * IA, int * JA, int * DESCA,
              double * X, int * IX, int * JX, int * DESCX, int * INCX,
              double * BETA,
              double * Y, int * IY, int * JY, int * DESCY, int * INCY )
#else
void pdsymv_( UPLO, N, ALPHA, A, IA, JA, DESCA, X, IX, JX, DESCX,
              INCX, BETA, Y, IY, JY, DESCY, INCY )
/*
*  .. Scalar Arguments ..
*/
   F_CHAR_T       UPLO;
   int            * IA, * INCX, * INCY, * IX, * IY, * JA, * JX, * JY,
                  * N;
   double         * ALPHA, * BETA;
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
*  PDSYMV  performs the matrix-vector operation
*
*     sub( Y ) := alpha*sub( A )*sub( X ) + beta*sub( Y ),
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
*  Alpha and beta are scalars, sub( X ) and sub( Y ) are n element  sub-
*  vectors and sub( A ) is an n by n symmetric submatrix.
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
*  UPLO    (global input) CHARACTER*1
*          On  entry,   UPLO  specifies  whether  the  local  pieces  of
*          the array  A  containing the  upper or lower triangular  part
*          of the symmetric submatrix  sub( A )  are to be referenced as
*          follows:
*
*             UPLO = 'U' or 'u'   Only the local pieces corresponding to
*                                 the   upper  triangular  part  of  the
*                                 symmetric submatrix sub( A ) are to be
*                                 referenced,
*
*             UPLO = 'L' or 'l'   Only the local pieces corresponding to
*                                 the   lower  triangular  part  of  the
*                                 symmetric submatrix sub( A ) are to be
*                                 referenced.
*
*  N       (global input) INTEGER
*          On entry,  N specifies the order of the  submatrix  sub( A ).
*          N must be at least zero.
*
*  ALPHA   (global input) DOUBLE PRECISION
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of the arrays  A
*          and X corresponding to the entries of the submatrix  sub( A )
*          and the subvector sub( X ) need not be set on input.
*
*  A       (local input) DOUBLE PRECISION array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix A.
*          Before  entry  with  UPLO = 'U' or 'u', this  array  contains
*          the  local  entries  of  the  upper  triangular  part  of the
*          symmetric  submatrix  sub( A ), and  the local entries of the
*          strictly lower triangular of sub( A ) are not referenced.
*          Before  entry  with  UPLO = 'L' or 'l', this  array  contains
*          the  local  entries  of  the  lower  triangular  part  of the
*          symmetric  submatrix  sub( A ), and  the local entries of the
*          strictly upper triangular of sub( A ) are not referenced.
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
*  BETA    (global input) DOUBLE PRECISION
*          On entry,  BETA  specifies the scalar  beta.   When  BETA  is
*          supplied  as  zero  then  the  local entries of  the array  Y
*          corresponding to the entries of the subvector  sub( Y )  need
*          not be set on input.
*
*  Y       (local input/local output) DOUBLE PRECISION array
*          On entry, Y is an array of dimension (LLD_Y, Ky), where LLD_Y
*          is   at  least  MAX( 1, Lr( 1, IY ) )  when  INCY = M_Y   and
*          MAX( 1, Lr( 1, IY+N-1 ) )  otherwise,  and,  Ky  is  at least
*          Lc( 1, JY+N-1 )  when  INCY = M_Y  and Lc( 1, JY ) otherwise.
*          Before  entry,  this  array contains the local entries of the
*          matrix Y.  On  exit,  sub( Y )  is overwritten by the updated
*          subvector.
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
   char           UploA, * one, top;
   int            Acol, Ai, Aii, Aimb1, Ainb1, Aj, Ajj, Akp, Akq, Ald, Amb,
                  Amp, Amp0, Anb, Anq, Anq0, Arow, XCfr, XCld, XRfr, XRld,
                  Xi, Xj, YCfr, YCpbY, YCsum, YCld, YRfr, YRpbY, YRsum,
                  YRld, Ycol, Yi, Yii, YisRow, Yj, Yjj, Yld, Ynp, Ynq, Yrow,
                  ctxt, info, ione=1, k, kb, ktmp, mycol, myrow, nb, npcol,
                  nprow, size, upper;
   double         * tbeta;
   PBTYP_T        * type;
/*
*  .. Local Arrays ..
*/
   int            Ad [DLEN_], Ad0[DLEN_], XCd[DLEN_], XRd[DLEN_], Xd[DLEN_],
                  YCd[DLEN_], YRd[DLEN_], Yd [DLEN_];
   char           * Aptr = NULL, * XC = NULL, * XR = NULL, * YC = NULL,
                  * YR   = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   upper = ( ( UploA  = Mupcase( F2C_CHAR( UPLO )[0] ) ) == CUPPER );
   PB_CargFtoC( *IA, *JA, DESCA, &Ai, &Aj, Ad );
   PB_CargFtoC( *IX, *JX, DESCX, &Xi, &Xj, Xd );
   PB_CargFtoC( *IY, *JY, DESCY, &Yi, &Yj, Yd );
#ifndef NO_ARGCHK
   Cblacs_gridinfo( ( ctxt = Ad[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Test the input parameters
*/
   if( !( info = ( ( nprow == -1 ) ? -( 701 + CTXT_ ) : 0 ) ) )
   {
      if( ( !upper ) && ( UploA != CLOWER ) )
      {
         PB_Cwarn( ctxt, __LINE__, __FILE__, "Illegal UPLO = %c\n", UploA );
         info = -1;
      }
      PB_Cchkmat( ctxt, "PDSYMV", "A", *N, 2, *N, 2, Ai, Aj, Ad,  7, &info );
      PB_Cchkvec( ctxt, "PDSYMV", "X", *N, 2, Xi, Xj, Xd, *INCX, 11, &info );
      PB_Cchkvec( ctxt, "PDSYMV", "Y", *N, 2, Yi, Yj, Yd, *INCY, 17, &info );
   }
   if( info ) { PB_Cabort( ctxt, "PDSYMV", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( ( *N == 0 ) ||
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
   type = PB_Cdtypeset();
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
            Ynq = PB_Cnumroc( *N, Yj, Yd[INB_], Yd[NB_], mycol, Yd[CSRC_],
                              npcol );
            if( Ynq > 0 )
            {
               Yld = Yd[LLD_];
               if( BETA[REAL_PART] == ZERO )
               {
                  dset_( &Ynq, ((char *) BETA), Mptr( ((char *) Y), Yii,
                         Yjj, Yld, type->size ), &Yld );
               }
               else
               {
                  dscal_( &Ynq, ((char *) BETA), Mptr( ((char *) Y), Yii,
                          Yjj, Yld, type->size ), &Yld );
               }
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
            Ynp = PB_Cnumroc( *N, Yi, Yd[IMB_], Yd[MB_], myrow, Yd[RSRC_],
                              nprow );
            if( Ynp > 0 )
            {
               if( BETA[REAL_PART] == ZERO )
               {
                  dset_( &Ynp, ((char *) BETA), Mptr( ((char *) Y), Yii,
                         Yjj, Yd[LLD_], type->size ), INCY );
               }
               else
               {
                  dscal_( &Ynp, ((char *) BETA), Mptr( ((char *) Y), Yii,
                          Yjj, Yd[LLD_], type->size ), INCY );
               }
            }
         }
      }
      return;
   }
/*
*  Compute descriptor Ad0 for sub( A )
*/
   PB_Cdescribe( *N, *N, Ai, Aj, Ad, nprow, npcol, myrow, mycol, &Aii, &Ajj,
                 &Ald, &Aimb1, &Ainb1, &Amb, &Anb, &Arow, &Acol, Ad0 );
/*
*  Reuse sub( Y ) and/or create vectors YR in process rows and YC in process
*  columns spanned by sub( A )
*/
   if( ( YisRow = ( *INCY == Yd[M_] ) ) != 0 )
   {
      PB_CInOutV( type, ROW,    *N, *N, Ad0, 1, ((char *)BETA), ((char *) Y),
                  Yi, Yj, Yd, ROW,    ((char**)(&tbeta)), &YR, YRd, &YRfr,
                  &YRsum, &YRpbY );
      PB_COutV( type, COLUMN, INIT, *N, *N, Ad0, 1, &YC, YCd, &YCfr, &YCsum );
   }
   else
   {
      PB_CInOutV( type, COLUMN, *N, *N, Ad0, 1, ((char *)BETA), ((char *) Y),
                  Yi, Yj, Yd, COLUMN, ((char**)(&tbeta)), &YC, YCd, &YCfr,
                  &YCsum, &YCpbY );
      PB_COutV( type, ROW,    INIT, *N, *N, Ad0, 1, &YR, YRd, &YRfr, &YRsum );
   }
/*
*  Replicate sub( X ) in process rows (XR) and process columns (XC) spanned by
*  sub( A )
*/
   if( *INCX == Xd[M_] )
   {
      PB_CInV( type, NOCONJG, ROW,    *N, *N, Ad0, 1, ((char *) X), Xi, Xj,  Xd,
               ROW, &XR, XRd, &XRfr );
      PB_CInV( type, NOCONJG, COLUMN, *N, *N, Ad0, 1, XR,            0,  0, XRd,
               ROW, &XC, XCd, &XCfr );
   }
   else
   {
      PB_CInV( type, NOCONJG, COLUMN, *N, *N, Ad0, 1, ((char *) X), Xi, Xj,  Xd,
               COLUMN, &XC, XCd, &XCfr );
      PB_CInV( type, NOCONJG, ROW,    *N, *N, Ad0, 1, XC,            0,  0, XCd,
               COLUMN, &XR, XRd, &XRfr );
   }

   one = type->one;
/*
*  Local matrix-vector multiply iff I own some data
*/
   Aimb1 = Ad0[IMB_ ]; Ainb1 = Ad0[INB_ ]; Amb = Ad0[MB_]; Anb = Ad0[NB_];
   Acol  = Ad0[CSRC_]; Arow  = Ad0[RSRC_];
   Amp   = PB_Cnumroc( *N, 0, Aimb1, Amb, myrow, Arow, nprow );
   Anq   = PB_Cnumroc( *N, 0, Ainb1, Anb, mycol, Acol, npcol );

   if( ( Amp > 0 ) && ( Anq > 0 ) )
   {
      size = type->size;
      Aptr = Mptr( ((char *) A), Aii, Ajj, Ald, size );

      XCld = XCd[LLD_]; XRld = XRd[LLD_]; YCld = YCd[LLD_]; YRld = YRd[LLD_];
/*
*  Scale YR or YC in the case sub( Y ) has been reused
*/
      if( YisRow )
      {
/*
*  YR resides in (a) process row(s)
*/
         if( !YRpbY )
         {
            if( ( myrow == YRd[RSRC_] ) || ( YRd[RSRC_] < 0 ) )
            {
/*
*  Make sure I own some data and scale YR
*/
               if( Anq > 0 )
               {
                  if( tbeta[REAL_PART] == ZERO )
                  {
                     dset_( &Anq, ((char *) tbeta), YR, &YRld );
                  }
                  else
                  {
                     dscal_( &Anq, ((char *) tbeta), YR, &YRld );
                  }
               }
            }
         }
      }
      else
      {
/*
*  YC resides in (a) process column(s)
*/
         if( !YCpbY )
         {
            if( ( mycol == YCd[CSRC_] ) || ( YCd[CSRC_] < 0 ) )
            {
/*
*  Make sure I own some data and scale YC
*/
               if( Amp > 0 )
               {
                  if( tbeta[REAL_PART] == ZERO )
                  {
                     dset_( &Amp, ((char *) tbeta), YC, &ione );
                  }
                  else
                  {
                     dscal_( &Amp, ((char *) tbeta), YC, &ione );
                  }
               }
            }
         }
      }
/*
*  Computational partitioning size is computed as the product of the logical
*  value returned by pilaenv_ and 2 * lcm( nprow, npcol ).
*/
      nb = 2 * pilaenv_( &ctxt, C2F_CHAR( &type->type ) ) *
           PB_Clcm( ( Arow >= 0 ? nprow : 1 ), ( Acol >= 0 ? npcol : 1 ) );

      if( upper )
      {
         for( k = 0; k < *N; k += nb )
         {
            kb   = *N - k; kb = MIN( kb, nb );
            Akp  = PB_Cnumroc( k,  0, Aimb1, Amb, myrow, Arow, nprow );
            Akq  = PB_Cnumroc( k,  0, Ainb1, Anb, mycol, Acol, npcol );
            Anq0 = PB_Cnumroc( kb, k, Ainb1, Anb, mycol, Acol, npcol );
            if( Akp > 0 && Anq0 > 0 )
            {
               dgemv_( C2F_CHAR( NOTRAN ), &Akp, &Anq0, ((char *)ALPHA),
                       Mptr( Aptr, 0, Akq, Ald, size ), &Ald, Mptr( XR, 0, Akq,
                       XRld, size ), &XRld, one, YC, &ione );
               dgemv_( C2F_CHAR( TRAN   ), &Akp, &Anq0, ((char *)ALPHA),
                       Mptr( Aptr, 0, Akq, Ald, size ), &Ald, XC, &ione, one,
                       Mptr( YR, 0, Akq, YRld, size ), &YRld );
            }
            PB_Cpsym( type, type, LEFT, UPPER, kb, 1, ((char *) ALPHA),
                      Aptr, k, k, Ad0, Mptr( XC, Akp, 0, XCld, size ), XCld,
                      Mptr( XR, 0, Akq, XRld, size ), XRld, Mptr( YC, Akp, 0,
                      YCld, size ), YCld, Mptr( YR, 0, Akq, YRld, size ), YRld,
                      PB_Ctzsymv );
         }
      }
      else
      {
         for( k = 0; k < *N; k += nb )
         {
            kb  = *N - k; ktmp = k + ( kb = MIN( kb, nb ) );
            Akp = PB_Cnumroc( k, 0, Aimb1, Amb, myrow, Arow, nprow );
            Akq = PB_Cnumroc( k, 0, Ainb1, Anb, mycol, Acol, npcol );
            PB_Cpsym( type, type, LEFT, LOWER, kb, 1, ((char *) ALPHA),
                      Aptr, k, k, Ad0, Mptr( XC, Akp, 0, XCld, size ), XCld,
                      Mptr( XR, 0, Akq, XRld, size ), XRld, Mptr( YC, Akp, 0,
                      YCld, size ), YCld, Mptr( YR, 0, Akq, YRld, size ), YRld,
                      PB_Ctzsymv );
            Akp  = PB_Cnumroc( ktmp, 0, Aimb1, Amb, myrow, Arow, nprow );
            Amp0 = Amp - Akp;
            Anq0 = PB_Cnumroc( kb,   k, Ainb1, Anb, mycol, Acol, npcol );
            if( Amp0 > 0 && Anq0 > 0 )
            {
               dgemv_( C2F_CHAR( NOTRAN ), &Amp0, &Anq0, ((char *) ALPHA),
                       Mptr( Aptr, Akp, Akq, Ald, size ), &Ald, Mptr( XR, 0,
                       Akq, XRld, size ), &XRld, one, Mptr( YC, Akp, 0, YCld,
                       size ), &ione );
               dgemv_( C2F_CHAR( TRAN   ), &Amp0, &Anq0, ((char *) ALPHA),
                       Mptr( Aptr, Akp, Akq, Ald, size ), &Ald, Mptr( XC, Akp,
                       0, XCld, size ), &ione, one, Mptr( YR, 0, Akq, YRld,
                       size ), &YRld );
            }
         }
      }
   }
   if( XCfr ) free( XC );
   if( XRfr ) free( XR );

   if( YisRow )
   {
/*
*  Combine the partial column results into YC
*/
      if( YCsum )
      {
         YCd[CSRC_] = 0;
         if( Amp > 0 )
         {
            top        = *PB_Ctop( &ctxt, COMBINE, ROW, TOP_GET );
            Cdgsum2d( ctxt, ROW, &top, Amp, 1, YC, YCd[LLD_], myrow, 0 );
         }
      }
/*
*  Combine the partial row results into YR
*/
      if( YRsum && ( Anq > 0 ) )
      {
         top = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
         Cdgsum2d( ctxt, COLUMN, &top, 1, Anq, YR, YRd[LLD_], YRd[RSRC_],
                   mycol );
      }

/*
*  YR := YR + YC
*/
      PB_Cpaxpby( type, NOCONJG, *N, 1, one, YC, 0, 0, YCd, COLUMN, one,
                  YR, 0, 0, YRd, ROW );
/*
*  sub( Y ) := beta * sub( Y ) + YR (if necessary)
*/
      if( YRpbY )
      {
         PB_Cpaxpby( type, NOCONJG, 1, *N, one, YR, 0, 0, YRd, ROW,
                     ((char *)BETA), ((char *) Y), Yi, Yj, Yd, ROW );
      }
   }
   else
   {
/*
*  Combine the partial row results into YR
*/
      if( YRsum )
      {
         YRd[RSRC_] = 0;
         if( Anq > 0 )
         {
            top = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
            Cdgsum2d( ctxt, COLUMN, &top, 1, Anq, YR, YRd[LLD_], 0,
                      mycol );
         }
      }
/*
*  Combine the partial column results into YC
*/
      if( YCsum && ( Amp > 0 ) )
      {
         top = *PB_Ctop( &ctxt, COMBINE, ROW, TOP_GET );
         Cdgsum2d( ctxt, ROW, &top, Amp, 1, YC, YCd[LLD_], myrow,
                   YCd[CSRC_] );
      }
/*
*  YC := YR + YC
*/
      PB_Cpaxpby( type, NOCONJG, 1, *N, one, YR, 0, 0, YRd, ROW, one,
                  YC, 0, 0, YCd, COLUMN );
/*
*  sub( Y ) := beta * sub( Y ) + YC (if necessary)
*/
      if( YCpbY )
      {
         PB_Cpaxpby( type, NOCONJG, *N, 1, one, YC, 0, 0, YCd,
                     COLUMN, ((char *)BETA), ((char *) Y), Yi, Yj, Yd,
                     COLUMN );
      }
   }
   if( YCfr ) free( YC );
   if( YRfr ) free( YR );
/*
*  End of PDSYMV
*/
}
