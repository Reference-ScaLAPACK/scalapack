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
void pssyr_( F_CHAR_T UPLO, int * N, float * ALPHA,
             float * X, int * IX, int * JX, int * DESCX, int * INCX,
             float * A, int * IA, int * JA, int * DESCA )
#else
void pssyr_( UPLO, N, ALPHA, X, IX, JX, DESCX, INCX, A, IA, JA, DESCA )
/*
*  .. Scalar Arguments ..
*/
   F_CHAR_T       UPLO;
   int            * IA, * INCX, * IX, * JA, * JX, * N;
   float          * ALPHA;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCX;
   float          * A, * X;
#endif
{
/*
*  Purpose
*  =======
*
*  PSSYR  performs the symmetric rank 1 operation
*
*     sub( A ) := alpha*sub( X )*sub( X )' + sub( A ),
*
*  where
*
*     sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1), and,
*
*     sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
*                      X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X.
*
*  Alpha is a scalar, sub( X ) is an n element subvector and sub( A ) is
*  an n by n symmetric submatrix.
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
*  ALPHA   (global input) REAL
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of the array   X
*          corresponding to the entries of the subvector  sub( X )  need
*          not be set on input.
*
*  X       (local input) REAL array
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
*  A       (local input/local output) REAL array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix A.
*          Before  entry  with  UPLO = 'U' or 'u', this  array  contains
*          the local entries corresponding to the upper triangular  part
*          of  the  symmetric submatrix  sub( A ), and the local entries
*          corresponding to the  strictly lower triangular  of  sub( A )
*          are not referenced.  On exit,  the upper  triangular  part of
*          sub( A ) is overwritten by the  upper triangular part  of the
*          updated submatrix.
*          Before  entry  with  UPLO = 'L' or 'l', this  array  contains
*          the local entries corresponding to the lower triangular  part
*          of  the  symmetric submatrix  sub( A ), and the local entries
*          corresponding to the  strictly upper triangular  of  sub( A )
*          are not referenced.  On exit,  the lower  triangular  part of
*          sub( A ) is overwritten by the  lower triangular part  of the
*          updated submatrix.
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
   char           UploA;
   int            Acol, Ai, Aii, Aimb1, Ainb1, Aj, Ajj, Akp, Akq, Ald, Amb,
                  Amp, Amp0, Anb, Anq, Anq0, Arow, XCfr, XCld, XRfr, XRld,
                  Xi, Xj, ctxt, info, ione=1, k, kb, ktmp, mycol, myrow, nb,
                  npcol, nprow, size, upper;
   PBTYP_T        * type;
/*
*  .. Local Arrays ..
*/
   int            Ad[DLEN_], Ad0[DLEN_], XCd0[DLEN_], XRd0[DLEN_], Xd[DLEN_];
   char           * Aptr = NULL, * XC = NULL, * XR = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   upper = ( ( UploA  = Mupcase( F2C_CHAR( UPLO )[0] ) ) == CUPPER );
   PB_CargFtoC( *IA, *JA, DESCA, &Ai, &Aj, Ad );
   PB_CargFtoC( *IX, *JX, DESCX, &Xi, &Xj, Xd );
#ifndef NO_ARGCHK
/*
*  Test the input parameters
*/
   Cblacs_gridinfo( ( ctxt = Xd[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
   if( !( info = ( ( nprow == -1 ) ? -( 701 + CTXT_ ) : 0 ) ) )
   {
      if( ( !upper ) && ( UploA != CLOWER ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PSSYR", "Illegal UPLO = %c\n", UploA );
         info = -1;
      }
      PB_Cchkvec( ctxt, "PSSYR", "X", *N, 2, Xi, Xj, Xd, *INCX,  7, &info );
      PB_Cchkmat( ctxt, "PSSYR", "A", *N, 2, *N, 2, Ai, Aj, Ad, 12, &info );
   }
   if( info ) { PB_Cabort( ctxt, "PSSYR", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( ( *N == 0 ) || ( ALPHA[REAL_PART] == ZERO ) )
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
   type = PB_Cstypeset();
/*
*  Compute descriptor Ad0 for sub( A )
*/
   PB_Cdescribe( *N, *N, Ai, Aj, Ad, nprow, npcol, myrow, mycol, &Aii, &Ajj,
                 &Ald, &Aimb1, &Ainb1, &Amb, &Anb, &Arow, &Acol, Ad0 );
/*
*  Replicate sub( X ) in process rows (XR) and process columns (XC) spanned by
*  sub( A )
*/
   if( *INCX == Xd[M_] )
   {
      PB_CInV( type, NOCONJG, ROW,    *N, *N, Ad0, 1, ((char *) X), Xi, Xj,
               Xd,   ROW,    &XR, XRd0, &XRfr );
      PB_CInV( type, NOCONJG, COLUMN, *N, *N, Ad0, 1, XR,            0,  0,
               XRd0, ROW,    &XC, XCd0, &XCfr );
   }
   else
   {
      PB_CInV( type, NOCONJG, COLUMN, *N, *N, Ad0, 1, ((char *) X), Xi, Xj,
               Xd,   COLUMN, &XC, XCd0, &XCfr );
      PB_CInV( type, NOCONJG, ROW,    *N, *N, Ad0, 1, XC,            0, 0,
               XCd0, COLUMN, &XR, XRd0, &XRfr );
   }
/*
*  Local rank-1 update if I own some data
*/
   Amp = PB_Cnumroc( *N, 0, Aimb1, Amb, myrow, Arow, nprow );
   Anq = PB_Cnumroc( *N, 0, Ainb1, Anb, mycol, Acol, npcol );

   if( ( Amp > 0 ) && ( Anq > 0 ) )
   {
      size = type->size;
      Aptr = Mptr( ((char *) A), Aii, Ajj, Ald, size );
/*
*  Computational partitioning size is computed as the product of the logical
*  value returned by pilaenv_ and 2 * lcm( nprow, npcol ).
*/
      nb   = 2 * pilaenv_( &ctxt, C2F_CHAR( &type->type ) ) *
             PB_Clcm( ( Arow >= 0 ? nprow : 1 ), ( Acol >= 0 ? npcol : 1 ) );

      XCld = XCd0[LLD_]; XRld = XRd0[LLD_];

      if( upper )
      {
         for( k = 0; k < *N; k += nb )
         {
            kb   = *N - k; kb = MIN( kb, nb );
            Akp  = PB_Cnumroc( k,  0, Aimb1, Amb, myrow, Arow, nprow );
            Akq  = PB_Cnumroc( k,  0, Ainb1, Anb, mycol, Acol, npcol );
            Anq0 = PB_Cnumroc( kb, k, Ainb1, Anb, mycol, Acol, npcol );
            if( Akp > 0 && Anq0 > 0 )
               sger_( &Akp, &Anq0, ((char *) ALPHA), XC, &ione,
                      Mptr( XR, 0, Akq, XRld, size ), &XRld, Mptr( Aptr, 0, Akq,
                      Ald, size ), &Ald );
            PB_Cpsyr( type, UPPER, kb, 1, ((char *) ALPHA), Mptr( XC, Akp, 0,
                      XCld, size ), XCld, Mptr( XR, 0, Akq, XRld, size ), XRld,
                      Aptr, k, k, Ad0, PB_Ctzsyr );
         }
      }
      else
      {
         for( k = 0; k < *N; k += nb )
         {
            kb = *N - k; ktmp = k + ( kb = MIN( kb, nb ) );
            Akp = PB_Cnumroc( k, 0, Aimb1, Amb, myrow, Arow, nprow );
            Akq = PB_Cnumroc( k, 0, Ainb1, Anb, mycol, Acol, npcol );
            PB_Cpsyr( type, LOWER, kb, 1, ((char *) ALPHA), Mptr( XC, Akp, 0,
                      XCld, size ), XCld, Mptr( XR, 0, Akq, XRld, size ), XRld,
                      Aptr, k, k, Ad0, PB_Ctzsyr );
            Akp  = PB_Cnumroc( ktmp, 0, Aimb1, Amb, myrow, Arow, nprow );
            Amp0 = Amp - Akp;
            Anq0 = PB_Cnumroc( kb,   k, Ainb1, Anb, mycol, Acol, npcol );
            if( Amp0 > 0 && Anq0 > 0 )
               sger_( &Amp0, &Anq0, ((char *) ALPHA), Mptr( XC, Akp,
                      0, XCld, size ), &ione, Mptr( XR, 0, Akq, XRld, size ),
                      &XRld, Mptr( Aptr, Akp, Akq, Ald, size ), &Ald );
         }
      }
   }
   if( XRfr ) free( XR );
   if( XCfr ) free( XC );
/*
*  End of PSSYR
*/
}
