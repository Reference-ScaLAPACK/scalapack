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
void PB_CptrsmAB( PBTYP_T * TYPE, char * VARIANT, char * SIDE, char * UPLO,
                  char * TRANSA, char * DIAG, Int M, Int N, char * ALPHA,
                  char * A, Int IA, Int JA, Int * DESCA, char * B, Int IB,
                  Int JB, Int * DESCB )
#else
void PB_CptrsmAB( TYPE, VARIANT, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A,
                  IA, JA, DESCA, B, IB, JB, DESCB )
/*
*  .. Scalar Arguments ..
*/
   char           * DIAG, * SIDE, * TRANSA, * UPLO, * VARIANT;
   Int            IA, IB, JA, JB, M, N;
   char           * ALPHA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA, * DESCB;
   char           * A, * B;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CptrsmAB  solves one of the matrix equations
*
*     op( sub( A ) )*X = alpha*sub( B ),   or
*
*     X*op( sub( A ) ) = alpha*sub( B ),
*
*  where
*
*     sub( A ) denotes A(IA:IA+M-1,JA:JA+M-1)  if SIDE = 'L',
*                      A(IA:IA+N-1,JA:JA+N-1)  if SIDE = 'R', and,
*
*     sub( B ) denotes B(IB:IB+M-1,JB:JB+N-1).
*
*  Alpha is a scalar, X and sub( B ) are m by n submatrices, sub( A ) is
*  a unit, or non-unit, upper or lower  triangular submatrix and op( Y )
*  is one of
*
*     op( Y ) = Y   or   op( Y ) = Y'   or   op( Y ) = conjg( Y' ).
*
*  The submatrix X is overwritten on sub( B ).
*
*  This  is  the  outer-product  algorithm using the logical aggregation
*  blocking technique.
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
*  VARIANT (global input) pointer to CHAR
*          On entry, VARIANT specifies whether the left-  or  right-loo-
*          king variant of the algorithm should be used for the transpo-
*          se cases only, that is TRANSA is not 'N' or 'n'. When VARIANT
*          is 'L' or 'l', the left-looking variant  is  used,  otherwise
*          the right-looking algorithm is selected.
*
*  SIDE    (global input) pointer to CHAR
*          On entry,  SIDE  specifies  whether op( sub( A ) ) appears on
*          the left or right of X as follows:
*
*             SIDE = 'L' or 'l'   op( sub( A ) )*X = alpha*sub( B ),
*
*             SIDE = 'R' or 'r'   X*op( sub( A ) ) = alpha*sub( B ).
*
*  UPLO    (global input) pointer to CHAR
*          On entry,  UPLO  specifies whether the submatrix  sub( A ) is
*          an upper or lower triangular submatrix as follows:
*
*             UPLO = 'U' or 'u'   sub( A ) is an upper triangular
*                                 submatrix,
*
*             UPLO = 'L' or 'l'   sub( A ) is a  lower triangular
*                                 submatrix.
*
*  TRANSA  (global input) pointer to CHAR
*          On entry,  TRANSA  specifies the form of op( sub( A ) ) to be
*          used in the matrix multiplication as follows:
*
*             TRANSA = 'N' or 'n'   op( sub( A ) ) = sub( A ),
*
*             TRANSA = 'T' or 't'   op( sub( A ) ) = sub( A )',
*
*             TRANSA = 'C' or 'c'   op( sub( A ) ) = conjg( sub( A )' ).
*
*  DIAG    (global input) pointer to CHAR
*          On entry,  DIAG  specifies  whether or not  sub( A )  is unit
*          triangular as follows:
*
*             DIAG = 'U' or 'u'  sub( A )  is  assumed to be unit trian-
*                                gular,
*
*             DIAG = 'N' or 'n'  sub( A ) is not assumed to be unit tri-
*                                angular.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( B ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( B ). N  must be at least zero.
*
*  ALPHA   (global input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of  the array  B
*          corresponding to the entries of the submatrix  sub( B )  need
*          not be set on input.
*
*  A       (local input) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+M-1 ) when  SIDE = 'L' or 'l'   and  is at
*          least Lc( 1, JA+N-1 ) otherwise.  Before  entry,  this  array
*          contains the local entries of the matrix A.
*          Before entry with  UPLO = 'U' or 'u', this array contains the
*          local entries corresponding to  the entries of the upper tri-
*          angular submatrix  sub( A ), and the local entries correspon-
*          ding to the entries of the strictly lower triangular part  of
*          the submatrix  sub( A )  are not referenced.
*          Before entry with  UPLO = 'L' or 'l', this array contains the
*          local entries corresponding to  the entries of the lower tri-
*          angular submatrix  sub( A ), and the local entries correspon-
*          ding to the entries of the strictly upper triangular part  of
*          the submatrix  sub( A )  are not referenced.
*          Note  that  when DIAG = 'U' or 'u', the local entries corres-
*          ponding to the  diagonal elements  of the submatrix  sub( A )
*          are not referenced either, but are assumed to be unity.
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
*  B       (local input/local output) pointer to CHAR
*          On entry, B is an array of dimension (LLD_B, Kb), where Kb is
*          at least Lc( 1, JB+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix  B.
*          On exit, the local entries of this array corresponding to the
*          to  the entries of the submatrix sub( B ) are  overwritten by
*          the local entries of the m by n  solution submatrix.
*
*  IB      (global input) INTEGER
*          On entry, IB  specifies B's global row index, which points to
*          the beginning of the submatrix sub( B ).
*
*  JB      (global input) INTEGER
*          On entry, JB  specifies B's global column index, which points
*          to the beginning of the submatrix sub( B ).
*
*  DESCB   (global and local input) INTEGER array
*          On entry, DESCB  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix B.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           conjg, * negone, * one, * talph, top, * zero;
   Int            Afr, Bcol, Bcurcol, Bcurimb1, Bcurinb1, Bcurrow, Bfr, Bii,
                  Bimb, Bimb1, Binb, Binb1, Bjj, Bld, Bmb, Bmp, Bmp0, Bnb, Bnq,
                  Bnq0, Brow, WAfr, WAoff, WBfr, WBsum, ctxt, k, kb, kbb, kmax,
                  ktmp, lside, mn, mycol, myrow, notran, npcol, nprow, size,
                  upper;
   GEMM_T         gemm;
   GSUM2D_T       gsum2d;
/*
*  .. Local Arrays ..
*/
   Int            Bd0[DLEN_], DBUFA[DLEN_], DBUFB[DLEN_], WAd[DLEN_],
                  WBd[DLEN_];
   char           * Aptr = NULL, * Bptr = NULL, * Bptr0 = NULL, * WA = NULL,
                  * WB   = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   lside  = ( Mupcase( SIDE  [0] ) == CLEFT   );
   upper  = ( Mupcase( UPLO  [0] ) == CUPPER  );
   notran = ( Mupcase( TRANSA[0] ) == CNOTRAN );
   size   = TYPE->size; negone = TYPE->negone;  one  = TYPE->one;
   zero   = TYPE->zero; gsum2d = TYPE->Cgsum2d; gemm = TYPE->Fgemm;
   talph  = ALPHA;
   kb     = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Retrieve sub( B )'s local information: Bii, Bjj, Brow, Bcol, Bld ...
*/
   Bimb = DESCB[IMB_]; Binb = DESCB[INB_];
   Bmb  = DESCB[MB_ ]; Bnb  = DESCB[NB_ ]; Bld = DESCB[LLD_];
   PB_Cinfog2l( IB, JB, DESCB, nprow, npcol, myrow, mycol, &Bii, &Bjj, &Brow,
                &Bcol );

   Bimb1 = PB_Cfirstnb( M, IB, Bimb, Bmb );
   Bmp0  = PB_Cnumroc( M, 0, Bimb1, Bmb, myrow, Brow, nprow );
   Binb1 = PB_Cfirstnb( N, JB, Binb, Bnb );
   Bnq0  = PB_Cnumroc( N, 0, Binb1, Bnb, mycol, Bcol, npcol );
   if( ( Bmp0 > 0 ) && ( Bnq0 > 0 ) ) Bptr0 = Mptr( B, Bii, Bjj, Bld, size );

   if( notran )
   {
      if( lside )
      {
         if( upper )
         {
            kmax = ( ( M - 1 ) / kb ) * kb;

            for( k = kmax; k >= 0; k -= kb )
            {
               kbb = M - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA:IA+ktmp-1, JA+k:JA+k+kbb-1 )
*/
               PB_CGatherV( TYPE, REUSE, BACKWARD, ktmp, kbb, A, IA, JA+k,
                            DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Replicate A(IA:IA+ktmp-1, JA+k:JA+k+kbb-1) over B(IB:IB+ktmp-1, JB:JB+N-1)
*/
               PB_Cdescset( Bd0, ktmp, N, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                            ctxt, Bld );
               PB_CInV( TYPE, NOCONJG, COLUMN, ktmp, N, Bd0, kbb, Aptr, 0, 0,
                        DBUFA, COLUMN, &WA, WAd, &WAfr );
/*
*  Solve B( IB+k:IB+ktmp-1, JB:JB+N-1 ) with talph
*/
               PB_CptrsmAB0( TYPE, SIDE, UPLO, DIAG, kbb, N, talph, WA, k, 0,
                             WAd, B, IB+k, JB, DESCB, &Bptr, DBUFB, &Bfr );
/*
*  Update B( IB:IB+k-1, JB:JB+N-1 )
*/
               if( k > 0 )
               {
/*
*  Replicate B( IB+k:IB+ktmp-1, JB:JB+N-1 ) over B( IB:IB+k-1, JB:JB+N-1 )
*/
                  PB_Cdescset( Bd0, k, N, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                               ctxt, Bld );
                  PB_CInV( TYPE, NOCONJG, ROW, k, N, Bd0, kbb, Bptr, 0, 0,
                           DBUFB, ROW, &WB, WBd, &WBfr );
/*
*  Local update
*/
                  Bmp = PB_Cnumroc( k, 0, Bimb1, Bmb, myrow, Brow, nprow );
                  if( ( Bmp > 0 ) && ( Bnq0 > 0 ) )
                     gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp, &Bnq0,
                           &kbb, negone, WA, &WAd[LLD_], WB, &WBd[LLD_], talph,
                           Bptr0, &Bld );
                  if( WBfr ) free( WB   );
                  talph = one;
               }
               if( WAfr ) free( WA   );
               if( Bfr  ) free( Bptr );
               if( Afr  ) free( Aptr );
            }
         }
         else
         {
            for( k = 0; k < M; k += kb )
            {
               ktmp = M - k; kbb = MIN( ktmp, kb );
/*
*  Accumulate A( IA+k:IA+M-1, JA+k:JA+k+kbb-1 )
*/
               PB_CGatherV( TYPE, REUSE, FORWARD, ktmp, kbb, A, IA+k, JA+k,
                            DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+M-1, JA+k:JA+k+kbb-1 ) over B( IB+k:IB+M-1, JB:JB+N-1 )
*/
               Bcurimb1 = PB_Cfirstnb( ktmp, IB+k, Bimb, Bmb );
               Bcurrow  = PB_Cindxg2p( k, Bimb1, Bmb, Brow, Brow, nprow );
               PB_Cdescset( Bd0, ktmp, N, Bcurimb1, Binb1, Bmb, Bnb, Bcurrow,
                            Bcol, ctxt, Bld );
               PB_CInV( TYPE, NOCONJG, COLUMN, ktmp, N, Bd0, kbb, Aptr, 0, 0,
                        DBUFA, COLUMN, &WA, WAd, &WAfr );
/*
*  Solve B( IB+k:IB+k+kbb-1, JB:JB+N-1 ) with talph
*/
               PB_CptrsmAB0( TYPE, SIDE, UPLO, DIAG, kbb, N, talph, WA, 0, 0,
                             WAd, B, IB+k, JB, DESCB, &Bptr, DBUFB, &Bfr );
/*
*  Update B( IB+k+kbb:IB+M-1, JB:JB+N-1 )
*/
               if( ( ktmp = ktmp - kbb ) > 0 )
               {
/*
*  Replicate B(IB+k:IB+k+kbb-1, JB:JB+N-1) over B(IB+k+kbb:IB+M-1, JB:JB+N-1)
*/
                  Bcurimb1 = PB_Cfirstnb( ktmp, IB+k+kbb, Bimb, Bmb );
                  Bcurrow  = PB_Cindxg2p( k+kbb, Bimb1, Bmb, Brow, Brow,
                                          nprow );
                  PB_Cdescset( Bd0, ktmp, N, Bcurimb1, Binb1, Bmb, Bnb, Bcurrow,
                               Bcol, ctxt, Bld );
                  PB_CInV( TYPE, NOCONJG, ROW, ktmp, N, Bd0, kbb, Bptr, 0, 0,
                           DBUFB, ROW, &WB, WBd, &WBfr );
/*
*  Local update
*/
                  Bmp = PB_Cnumroc( ktmp, k+kbb, Bimb1, Bmb, myrow, Brow,
                                    nprow );
                  if( ( Bmp > 0 ) && ( Bnq0 > 0 ) )
                  {
                     WAoff = PB_Cnumroc( kbb, 0, WAd[IMB_], WAd[MB_], myrow,
                                         WAd[RSRC_], nprow );
                     gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp, &Bnq0,
                           &kbb, negone, Mptr( WA, WAoff, 0, WAd[LLD_], size ),
                           &WAd[LLD_], WB, &WBd[LLD_], talph, Mptr( Bptr0,
                           Bmp0-Bmp, 0, Bld, size ), &Bld );
                  }
                  if( WBfr ) free( WB   );
                  talph = one;
               }
               if( WAfr ) free( WA   );
               if( Bfr  ) free( Bptr );
               if( Afr  ) free( Aptr );
            }
         }
      }
      else
      {
         if( upper )
         {
            for( k = 0; k < N; k += kb )
            {
               ktmp = N - k; kbb = MIN( ktmp, kb );
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 )
*/
               PB_CGatherV( TYPE, REUSE, FORWARD, kbb, ktmp, A, IA+k, JA+k,
                            DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 ) over B( IB:IB+M-1, JB+k:JB+N-1 )
*/
               Bcurinb1 = PB_Cfirstnb( ktmp, JB+k, Binb, Bnb );
               Bcurcol  = PB_Cindxg2p( k, Binb1, Bnb, Bcol, Bcol, npcol );
               PB_Cdescset( Bd0, M, ktmp, Bimb1, Bcurinb1, Bmb, Bnb, Brow,
                            Bcurcol, ctxt, Bld );
               PB_CInV( TYPE, NOCONJG, ROW, M, ktmp, Bd0, kbb, Aptr, 0, 0,
                        DBUFA, ROW, &WA, WAd, &WAfr );
/*
*  Solve B( IB:IB+M-1, JB+k:JB+k+kbb-1 ) with talph
*/
               PB_CptrsmAB0( TYPE, SIDE, UPLO, DIAG, M, kbb, talph, WA, 0, 0,
                             WAd, B, IB, JB+k, DESCB, &Bptr, DBUFB, &Bfr );
/*
*  Update B( IB:IB+M-1, JB+k+kbb:JB+N-1 )
*/
               if( ( ktmp = ktmp - kbb ) > 0 )
               {
/*
*  Replicate B(IB:IB+M-1, JB+k:JB+k+kbb-1) over B(IB:IB+M-1, JB+k+kbb:JB+N-1)
*/
                  Bcurinb1 = PB_Cfirstnb( ktmp, JB+k+kbb, Binb, Bnb );
                  Bcurcol  = PB_Cindxg2p( k+kbb, Binb1, Bnb, Bcol, Bcol,
                                          npcol );
                  PB_Cdescset( Bd0, M, ktmp, Bimb1, Bcurinb1, Bmb, Bnb, Brow,
                               Bcurcol, ctxt, Bld );
                  PB_CInV( TYPE, NOCONJG, COLUMN, M, ktmp, Bd0, kbb, Bptr, 0, 0,
                           DBUFB, COLUMN, &WB, WBd, &WBfr );
/*
*  Local update
*/
                  Bnq = PB_Cnumroc( ktmp, k+kbb, Binb1, Bnb, mycol, Bcol,
                                    npcol );
                  if( ( Bmp0 > 0 ) && ( Bnq > 0 ) )
                  {
                     WAoff = PB_Cnumroc( kbb, 0, WAd[INB_], WAd[NB_], mycol,
                                         WAd[CSRC_], npcol );
                     gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp0, &Bnq,
                           &kbb, negone, WB, &WBd[LLD_], Mptr( WA, 0, WAoff,
                           WAd[LLD_], size ), &WAd[LLD_], talph, Mptr( Bptr0,
                           0, Bnq0-Bnq, Bld, size ), &Bld );
                  }
                  if( WBfr ) free( WB   );
                  talph = one;
               }
               if( WAfr ) free( WA   );
               if( Bfr  ) free( Bptr );
               if( Afr  ) free( Aptr );
            }
         }
         else
         {
            kmax = ( ( N - 1 ) / kb ) * kb;

            for( k = kmax; k >= 0; k -= kb )
            {
               kbb = N - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA:JA+ktmp-1 )
*/
               PB_CGatherV( TYPE, REUSE, BACKWARD, kbb, ktmp, A, IA+k, JA,
                            DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+k+kbb-1, JA:JA+ktmp-1 ) over B(IB:IB+M-1, JB:JB+ktmp-1)
*/
               PB_Cdescset( Bd0, M, ktmp, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                            ctxt, Bld );
               PB_CInV( TYPE, NOCONJG, ROW, M, ktmp, Bd0, kbb, Aptr, 0, 0,
                        DBUFA, ROW, &WA, WAd, &WAfr );
/*
*  Solve B( IB:IB+M-1, JB+k:JB+ktmp-1 ) with talph
*/
               PB_CptrsmAB0( TYPE, SIDE, UPLO, DIAG, M, kbb, talph, WA, 0, k,
                             WAd, B, IB, JB+k, DESCB, &Bptr, DBUFB, &Bfr );
/*
*  Update B( IB:IB+M-1, JB:JB+k-1 )
*/
               if( k > 0 )
               {
/*
*  Replicate B( IB:IB+M-1, JB+k:JB+ktmp-1 ) over B( IB:IB+M-1, JB:JB+k-1 )
*/
                  PB_Cdescset( Bd0, M, k, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                               ctxt, Bld );
                  PB_CInV( TYPE, NOCONJG, COLUMN, M, k, Bd0, kbb, Bptr, 0, 0,
                           DBUFB, COLUMN, &WB, WBd, &WBfr );
/*
*  Local update
*/
                  Bnq = PB_Cnumroc( k, 0, Binb1, Bnb, mycol, Bcol, npcol );
                  if( ( Bmp0 > 0 ) && ( Bnq > 0 ) )
                     gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp0, &Bnq,
                           &kbb, negone, WB, &WBd[LLD_], WA, &WAd[LLD_], talph,
                           Bptr0, &Bld );
                  if( WBfr ) free( WB   );
                  talph = one;
               }
               if( WAfr ) free( WA   );
               if( Bfr  ) free( Bptr );
               if( Afr  ) free( Aptr );
            }
         }
      }
   }
   else
   {
      if( Mupcase( VARIANT[0] ) == CRIGHT )
      {
/*
*  Right looking variant for the transpose cases
*/
         conjg = ( ( Mupcase( TRANSA[0] ) == CCOTRAN ) ? CCONJG : CNOCONJG );

         if( lside )
         {
            if( !upper )
            {
/*
*  Left Lower (Conjugate) Transpose
*/
               kmax = ( ( M - 1 ) / kb ) * kb;

               for( k = kmax; k >= 0; k -= kb )
               {
                  kbb = M - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA:JA+ktmp-1 )
*/
                  PB_CGatherV( TYPE, REUSE, BACKWARD, kbb, ktmp, A, IA+k, JA,
                               DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+k+kbb-1, JA:JA+ktmp-1 )' over B(IB:IB+ktmp-1, JB:JB+N-1)
*/
                  PB_Cdescset( Bd0, ktmp, N, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                               ctxt, Bld );
                  PB_CInV( TYPE, &conjg, COLUMN, ktmp, N, Bd0, kbb, Aptr, 0, 0,
                           DBUFA, ROW, &WA, WAd, &WAfr );
/*
*  Solve B( IB+k:IB+ktmp-1, JB:JB+N-1 ) with talph
*/
                  PB_CptrsmAB0( TYPE, SIDE, UPPER, DIAG, kbb, N, talph, WA, k,
                                0, WAd, B, IB+k, JB, DESCB, &Bptr, DBUFB,
                                &Bfr );
/*
*  Update B( IB:IB+k-1, JB:JB+N-1 )
*/
                  if( k > 0 )
                  {
/*
*  Replicate B( IB+k:IB+ktmp-1, JB:JB+N-1 ) over B( IB:IB+k-1, JB:JB+N-1 )
*/
                     PB_Cdescset( Bd0, k, N, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                                  ctxt, Bld );
                     PB_CInV( TYPE, NOCONJG, ROW, k, N, Bd0, kbb, Bptr, 0, 0,
                              DBUFB, ROW, &WB, WBd, &WBfr );
/*
*  Local update
*/
                     Bmp = PB_Cnumroc( k, 0, Bimb1, Bmb, myrow, Brow, nprow );
                     if( ( Bmp > 0 ) && ( Bnq0 > 0 ) )
                        gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp,
                              &Bnq0, &kbb, negone, WA, &WAd[LLD_], WB,
                              &WBd[LLD_], talph, Bptr0, &Bld );
                     if( WBfr ) free( WB   );
                     talph = one;
                  }
                  if( WAfr ) free( WA   );
                  if( Bfr  ) free( Bptr );
                  if( Afr  ) free( Aptr );
               }
            }
            else
            {
/*
*  Left Upper (Conjugate) Transpose
*/
               for( k = 0; k < M; k += kb )
               {
                  ktmp = M - k; kbb = MIN( ktmp, kb );
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA+k:JA+M-1 )
*/
                  PB_CGatherV( TYPE, REUSE, FORWARD, kbb, ktmp, A, IA+k, JA+k,
                               DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+k+kbb-1, JA+k:JA+M-1 )' over B( IB+k:IB+M-1, JB:JB+N-1 )
*/
                  Bcurimb1 = PB_Cfirstnb( ktmp, IB+k, Bimb, Bmb );
                  Bcurrow  = PB_Cindxg2p( k, Bimb1, Bmb, Brow, Brow, nprow );
                  PB_Cdescset( Bd0, ktmp, N, Bcurimb1, Binb1, Bmb, Bnb, Bcurrow,
                               Bcol, ctxt, Bld );
                  PB_CInV( TYPE, &conjg, COLUMN, ktmp, N, Bd0, kbb, Aptr, 0, 0,
                           DBUFA, ROW, &WA, WAd, &WAfr );
/*
*  Solve B( IB+k:IB+k+kbb-1, JB:JB+N-1 ) with talph
*/
                  PB_CptrsmAB0( TYPE, SIDE, LOWER, DIAG, kbb, N, talph, WA, 0,
                                0, WAd, B, IB+k, JB, DESCB, &Bptr, DBUFB,
                                &Bfr );
/*
*  Update B( IB+k+kbb:IB+M-1, JB:JB+N-1 )
*/
                  if( ( ktmp = ktmp - kbb ) > 0 )
                  {
/*
*  Replicate B(IB+k:IB+k+kbb-1, JB:JB+N-1) over B(IB+k+kbb:IB+M-1, JB:JB+N-1)
*/
                     Bcurimb1 = PB_Cfirstnb( ktmp, IB+k+kbb, Bimb, Bmb );
                     Bcurrow  = PB_Cindxg2p( k+kbb, Bimb1, Bmb, Brow, Brow,
                                             nprow );
                     PB_Cdescset( Bd0, ktmp, N, Bcurimb1, Binb1, Bmb, Bnb,
                                  Bcurrow, Bcol, ctxt, Bld );
                     PB_CInV( TYPE, NOCONJG, ROW, ktmp, N, Bd0, kbb, Bptr, 0, 0,
                              DBUFB, ROW, &WB, WBd, &WBfr );
/*
*  Local update
*/
                     Bmp = PB_Cnumroc( ktmp, k+kbb, Bimb1, Bmb, myrow, Brow,
                                       nprow );
                     if( ( Bmp > 0 ) && ( Bnq0 > 0 ) )
                     {
                        WAoff = PB_Cnumroc( kbb, 0, WAd[IMB_], WAd[MB_], myrow,
                                            WAd[RSRC_], nprow );
                        gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp,
                              &Bnq0, &kbb, negone, Mptr( WA, WAoff, 0,
                              WAd[LLD_], size ), &WAd[LLD_], WB, &WBd[LLD_],
                              talph, Mptr( Bptr0, Bmp0-Bmp, 0, Bld, size ),
                              &Bld );
                     }
                     if( WBfr ) free( WB   );
                     talph = one;
                  }
                  if( WAfr ) free( WA   );
                  if( Bfr  ) free( Bptr );
                  if( Afr  ) free( Aptr );
               }
            }
         }
         else
         {
            if( !upper )
            {
/*
*  Right Lower (Conjugate) Transpose
*/
               for( k = 0; k < N; k += kb )
               {
                  ktmp = N - k; kbb = MIN( ktmp, kb );
/*
*  Accumulate A( IA+k:IA+N-1, JA+k:JA+k+kbb-1 )
*/
                  PB_CGatherV( TYPE, REUSE, FORWARD, ktmp, kbb, A, IA+k, JA+k,
                               DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+N-1, JA+k:JA+k+kbb-1 )' over B( IB:IB+M-1, JB+k:JB+N-1 )
*/
                  Bcurinb1 = PB_Cfirstnb( ktmp, JB+k, Binb, Bnb );
                  Bcurcol  = PB_Cindxg2p( k, Binb1, Bnb, Bcol, Bcol, npcol );
                  PB_Cdescset( Bd0, M, ktmp, Bimb1, Bcurinb1, Bmb, Bnb, Brow,
                               Bcurcol, ctxt, Bld );
                  PB_CInV( TYPE, &conjg, ROW, M, ktmp, Bd0, kbb, Aptr, 0, 0,
                           DBUFA, COLUMN, &WA, WAd, &WAfr );
/*
*  Solve B( IB:IB+M-1, JB+k:JB+k+kbb-1 ) with talph
*/
                  PB_CptrsmAB0( TYPE, SIDE, UPPER, DIAG, M, kbb, talph, WA, 0,
                                0, WAd, B, IB, JB+k, DESCB, &Bptr, DBUFB,
                                &Bfr );
/*
*  Update B( IB:IB+M-1, JB+k+kbb:JB+N-1 )
*/
                  if( ( ktmp = ktmp - kbb ) > 0 )
                  {
/*
*  Replicate B(IB:IB+M-1, JB+k:JB+k+kbb-1) over B(IB:IB+M-1, JB+k+kbb:JB+N-1)
*/
                     Bcurinb1 = PB_Cfirstnb( ktmp, JB+k+kbb, Binb, Bnb );
                     Bcurcol  = PB_Cindxg2p( k+kbb, Binb1, Bnb, Bcol, Bcol,
                                             npcol );
                     PB_Cdescset( Bd0, M, ktmp, Bimb1, Bcurinb1, Bmb, Bnb, Brow,
                                  Bcurcol, ctxt, Bld );
                     PB_CInV( TYPE, NOCONJG, COLUMN, M, ktmp, Bd0, kbb, Bptr,
                              0, 0, DBUFB, COLUMN, &WB, WBd, &WBfr );
/*
*  Local update
*/
                     Bnq = PB_Cnumroc( ktmp, k+kbb, Binb1, Bnb, mycol, Bcol,
                                       npcol );
                     if( ( Bmp0 > 0 ) && ( Bnq > 0 ) )
                     {
                        WAoff = PB_Cnumroc( kbb, 0, WAd[INB_], WAd[NB_], mycol,
                                            WAd[CSRC_], npcol );
                        gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp0,
                              &Bnq, &kbb, negone, WB, &WBd[LLD_], Mptr( WA, 0,
                              WAoff, WAd[LLD_], size ), &WAd[LLD_], talph,
                              Mptr( Bptr0, 0, Bnq0-Bnq, Bld, size ), &Bld );
                     }
                     if( WBfr ) free( WB   );
                     talph = one;
                  }
                  if( WAfr ) free( WA   );
                  if( Bfr  ) free( Bptr );
                  if( Afr  ) free( Aptr );
               }
            }
            else
            {
/*
*  Right Upper (Conjugate) Transpose
*/
               kmax = ( ( N - 1 ) / kb ) * kb;

               for( k = kmax; k >= 0; k -= kb )
               {
                  kbb = N - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA:IA+ktmp-1, JA+k:JA+k+kbb-1 )
*/
                  PB_CGatherV( TYPE, REUSE, BACKWARD, ktmp, kbb, A, IA, JA+k,
                               DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA:IA+ktmp-1, JA+k:JA+k+kbb-1 )' over B(IB:IB+M-1, JB:JB+ktmp-1)
*/
                  PB_Cdescset( Bd0, M, ktmp, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                               ctxt, Bld );
                  PB_CInV( TYPE, &conjg, ROW, M, ktmp, Bd0, kbb, Aptr, 0, 0,
                           DBUFA, COLUMN, &WA, WAd, &WAfr );
/*
*  Solve B( IB:IB+M-1, JB+k:JB+ktmp-1 ) with talph
*/
                  PB_CptrsmAB0( TYPE, SIDE, LOWER, DIAG, M, kbb, talph, WA, 0,
                                k, WAd, B, IB, JB+k, DESCB, &Bptr, DBUFB,
                                &Bfr );
/*
*  Update B( IB:IB+M-1, JB:JB+k-1 )
*/
                  if( k > 0 )
                  {
/*
*  Replicate B( IB:IB+M-1, JB+k:JB+ktmp-1 ) over B( IB:IB+M-1, JB:JB+k-1 )
*/
                     PB_Cdescset( Bd0, M, k, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                                  ctxt, Bld );
                     PB_CInV( TYPE, NOCONJG, COLUMN, M, k, Bd0, kbb, Bptr, 0, 0,
                              DBUFB, COLUMN, &WB, WBd, &WBfr );
/*
*  Local update
*/
                     Bnq = PB_Cnumroc( k, 0, Binb1, Bnb, mycol, Bcol, npcol );
                     if( ( Bmp0 > 0 ) && ( Bnq > 0 ) )
                        gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Bmp0,
                              &Bnq, &kbb, negone, WB, &WBd[LLD_], WA,
                              &WAd[LLD_], talph, Bptr0, &Bld );
                     if( WBfr ) free( WB   );
                     talph = one;
                  }
                  if( WAfr ) free( WA   );
                  if( Bfr  ) free( Bptr );
                  if( Afr  ) free( Aptr );
               }
            }
         }
      }
      else
      {
/*
*  Left looking variant for the transpose cases
*/
         if( lside )
         {
            top  = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );

            if( upper )
            {
/*
*  Accumulate A( IA:IA+Bimb1-1, JA:JA+Bimb1-1 )
*/
               PB_CGatherV( TYPE, REUSE, FORWARD, Bimb1, Bimb1, A, IA, JA,
                            DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA:IA+Bimb1-1, JA:JA+Bimb1-1 ) over B(IB:IB+Bimb1-1, JB:JB+N-1)
*/
               PB_Cdescset( Bd0, Bimb1, N, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                            ctxt, Bld );
               PB_CInV( TYPE, NOCONJG, COLUMN, Bimb1, N, Bd0, Bimb1, Aptr, 0, 0,
                        DBUFA, COLUMN, &WA, WAd, &WAfr );
/*
*  Solve B( IB:IB+Bimb1-1, JB:JB+N-1 )
*/
               if( ( ( Brow < 0 ) || ( myrow == Brow ) ) && ( Bnq0 > 0 ) )
                  TYPE->Ftrsm( C2F_CHAR( SIDE   ), C2F_CHAR( UPLO ),
                               C2F_CHAR( TRANSA ), C2F_CHAR( DIAG ),
                               &Bimb1, &Bnq0, ALPHA, WA, &WAd[LLD_],
                               Bptr0, &Bld );
               if( WAfr ) free( WA   );
               if( Afr  ) free( Aptr );
/*
*  Update and solve remaining rows of sub( B )
*/
               for( k = Bimb1; k < M; k += kb )
               {
                  kbb = M - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA:IA+ktmp-1, JA+k:JA+ktmp-1 )
*/
                  PB_CGatherV( TYPE, REUSE, FORWARD, ktmp, kbb, A, IA, JA+k,
                               DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA:IA+ktmp-1, JA+k:JA+ktmp-1 ) over B(IB:IB+ktmp-1, JB:JB+N-1)
*/
                  PB_Cdescset( Bd0, ktmp, N, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                               ctxt, Bld );
                  PB_CInV( TYPE, NOCONJG, COLUMN, ktmp, N, Bd0, kbb, Aptr, 0, 0,
                           DBUFA, COLUMN, &WA, WAd, &WAfr );
/*
*  WB := A( IA:IA+k-1, JA+k:JA+ktmp-1 )' * B( IB:IB+k-1, JB:JB+N-1 )
*/
                  PB_COutV( TYPE, ROW, INIT, ktmp, N, Bd0, kbb, &WB, WBd, &WBfr,
                            &WBsum );
                  Bmp = PB_Cnumroc( k, 0, Bimb1, Bmb, myrow, Brow, nprow );
                  if( ( Bnq0 > 0 ) && ( Bmp > 0 ) )
                     gemm( C2F_CHAR( TRANSA ), C2F_CHAR( NOTRAN ), &kbb, &Bnq0,
                           &Bmp, one, WA, &WAd[LLD_], Bptr0, &Bld, zero, WB,
                           &WBd[LLD_] );
                  if( WBsum )
                  {
                     WBd[RSRC_] = PB_Cindxg2p( k, Bimb1, Bmb, Brow, Brow,
                                               nprow );
                     if( Bnq0 > 0 )
                        gsum2d( ctxt, COLUMN, &top, kbb, Bnq0, WB, WBd[LLD_],
                                WBd[RSRC_], mycol );
                  }
/*
*  Add WB to B( IB+k:IB+ktmp-1, JB:JB+N-1 ) and solve it with
*  A( IA+k:IA+ktmp-1, JA+k:JA+ktmp-1 )
*/
                  PB_CptrsmAB1( TYPE, SIDE, UPLO, TRANSA, DIAG, kbb, N, ALPHA,
                                WA, k, 0, WAd, B, IB+k, JB, DESCB, WB, WBd );
                  if( WBfr ) free( WB   );
                  if( WAfr ) free( WA   );
                  if( Afr  ) free( Aptr );
               }
            }
            else
            {
/*
*  Solve last block of rows of sub( B )
*/
               Bcurimb1 = PB_Clastnb( M, IB, Bimb, Bmb );
               k        = M - Bcurimb1;
/*
*  Accumulate A( IA+k:IA+M-1, JA+k:JA+M-1 )
*/
               PB_CGatherV( TYPE, REUSE, BACKWARD, Bcurimb1, Bcurimb1, A, IA+k,
                            JA+k, DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+M-1, JA+k:JA+M-1 ) over B( IB+k:IB+M-1, JB:JB+N-1 )
*/
               Bcurrow = PB_Cindxg2p( k, Bimb1, Bmb, Brow, Brow, nprow );
               PB_Cdescset( Bd0, Bcurimb1, N, Bcurimb1, Binb1, Bmb, Bnb,
                            Bcurrow, Bcol, ctxt, Bld );
               PB_CInV( TYPE, NOCONJG, COLUMN, Bcurimb1, N, Bd0, Bcurimb1, Aptr,
                        0, 0, DBUFA, COLUMN, &WA, WAd, &WAfr );
/*
*  Solve B( IB+k:IB+M-1, JB:JB+N-1 )
*/
               if( ( ( Brow < 0 ) || ( myrow == Bcurrow ) ) && ( Bnq0 > 0 ) )
                  TYPE->Ftrsm( C2F_CHAR( SIDE   ), C2F_CHAR( UPLO ),
                               C2F_CHAR( TRANSA ), C2F_CHAR( DIAG ),
                               &Bcurimb1, &Bnq0, ALPHA, WA, &WAd[LLD_],
                               Mptr( Bptr0, Bmp0-Bcurimb1, 0, Bld, size ),
                               &Bld );
               if( WAfr ) free( WA   );
               if( Afr  ) free( Aptr );
               if( ( mn = M - Bcurimb1 ) <= 0 ) return;
/*
*  Update and solve remaining rows of sub( B )
*/
               kmax = ( ( mn - 1 ) / kb ) * kb;

               for( k = kmax; k >= 0; k -= kb )
               {
                  ktmp = M - k; kbb = mn - k; kbb = MIN( kbb, kb );
/*
*  Accumulate A( IA+k:IA+M-1, JA+k:JA+k+kbb-1 )
*/
                  PB_CGatherV( TYPE, REUSE, BACKWARD, ktmp, kbb, A, IA+k, JA+k,
                               DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+M-1, JA+k:JA+k+kbb-1 ) over B( IB+k:IB+M-1, JB:JB+N-1 )
*/
                  Bcurimb1 = PB_Cfirstnb( ktmp, IB+k, Bimb, Bmb );
                  Bcurrow  = PB_Cindxg2p( k, Bimb1, Bmb, Brow, Brow, nprow );
                  PB_Cdescset( Bd0, ktmp, N, Bcurimb1, Binb1, Bmb, Bnb, Bcurrow,
                               Bcol, ctxt, Bld );
                  PB_CInV( TYPE, NOCONJG, COLUMN, ktmp, N, Bd0, kbb, Aptr, 0, 0,
                           DBUFA, COLUMN, &WA, WAd, &WAfr );
/*
*  WB := A( IA+k+kbb:IA+M-1, JA+k:JA+k+kbb-1 )'* B( IB+k+kbb:IB+M-1, JB:JB+N-1 )
*/
                  PB_COutV( TYPE, ROW, INIT, ktmp, N, Bd0, kbb, &WB, WBd, &WBfr,
                            &WBsum );
                  Bmp = PB_Cnumroc( ktmp-kbb, k+kbb, Bimb1, Bmb, myrow, Brow,
                                    nprow );
                  if( ( Bnq0 > 0 ) && ( Bmp > 0 ) )
                  {
                     WAoff = PB_Cnumroc( kbb, 0, WAd[IMB_], WAd[MB_], myrow,
                                         WAd[RSRC_], nprow );
                     gemm( C2F_CHAR( TRANSA ), C2F_CHAR( NOTRAN ), &kbb, &Bnq0,
                           &Bmp, one, Mptr( WA, WAoff, 0, WAd[LLD_], size ),
                           &WAd[LLD_], Mptr( Bptr0, Bmp0-Bmp, 0, Bld, size ),
                           &Bld, zero, WB, &WBd[LLD_] );
                  }
                  if( WBsum )
                  {
                     WBd[RSRC_] = PB_Cindxg2p( k + kbb - 1, Bimb1, Bmb, Brow,
                                               Brow, nprow );
                     if( Bnq0 > 0 )
                        gsum2d( ctxt, COLUMN, &top, kbb, Bnq0, WB, WBd[LLD_],
                                WBd[RSRC_], mycol );
                  }
/*
*  Add WB to B( IB+k:IB+k+kbb-1, JB:JB+N-1 ) and solve it with
*  A( IA+k:IA+k+kbb-1, JA+k:JA+k+kbb-1 )
*/
                  PB_CptrsmAB1( TYPE, SIDE, UPLO, TRANSA, DIAG, kbb, N, ALPHA,
                                WA, 0, 0, WAd, B, IB+k, JB, DESCB, WB, WBd );
                  if( WBfr ) free( WB   );
                  if( WAfr ) free( WA   );
                  if( Afr  ) free( Aptr );
               }
            }
         }
         else
         {
            top  = *PB_Ctop( &ctxt, COMBINE, ROW, TOP_GET );

            if( upper )
            {
/*
*  Solve last block of columns of sub( B )
*/
               Bcurinb1 = PB_Clastnb( N, JB, Binb, Bnb );
               k        = N - Bcurinb1;
/*
*  Accumulate A( IA+k:IA+N-1, JA+k:JA+N-1 )
*/
               PB_CGatherV( TYPE, REUSE, BACKWARD, Bcurinb1, Bcurinb1, A, IA+k,
                            JA+k, DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+N-1, JA+k:JA+N-1 ) over B( IB:IB+M-1, JB+k:JB+N-1 )
*/
               Bcurcol = PB_Cindxg2p( k, Binb1, Bnb, Bcol, Bcol, npcol );
               PB_Cdescset( Bd0, M, Bcurinb1, Bimb1, Bcurinb1, Bmb, Bnb, Brow,
                            Bcurcol, ctxt, Bld );
               PB_CInV( TYPE, NOCONJG, ROW, M, Bcurinb1, Bd0, Bcurinb1, Aptr,
                        0, 0, DBUFA, ROW, &WA, WAd, &WAfr );
/*
*  Solve B( IB:IB+M-1, JB+k:JB+N-1 )
*/
               if( ( ( Bcol < 0 ) || ( mycol == Bcurcol ) ) && ( Bmp0 > 0 ) )
                  TYPE->Ftrsm( C2F_CHAR( SIDE   ), C2F_CHAR( UPLO ),
                               C2F_CHAR( TRANSA ), C2F_CHAR( DIAG ),
                               &Bmp0, &Bcurinb1, ALPHA, WA, &WAd[LLD_],
                               Mptr( Bptr0, 0, Bnq0-Bcurinb1, Bld, size ),
                               &Bld );
               if( WAfr ) free( WA   );
               if( Afr  ) free( Aptr );
               if( ( mn = N - Bcurinb1 ) <= 0 ) return;
/*
*  Update and solve remaining columns of sub( B )
*/
               kmax = ( ( mn - 1 ) / kb ) * kb;

               for( k = kmax; k >= 0; k -= kb )
               {
                  ktmp = N - k; kbb = mn - k; kbb = MIN( kbb, kb );
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 )
*/
                  PB_CGatherV( TYPE, REUSE, BACKWARD, kbb, ktmp, A, IA+k, JA+k,
                               DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 ) over B( IB:IB+M-1, JB+k:JB+N-1 )
*/
                  Bcurinb1 = PB_Cfirstnb( ktmp, JB+k, Binb, Bnb );
                  Bcurcol  = PB_Cindxg2p( k, Binb1, Bnb, Bcol, Bcol, npcol );
                  PB_Cdescset( Bd0, M, ktmp, Bimb1, Bcurinb1, Bmb, Bnb, Brow,
                               Bcurcol, ctxt, Bld );
                  PB_CInV( TYPE, NOCONJG, ROW, M, ktmp, Bd0, kbb, Aptr, 0, 0,
                           DBUFA, ROW, &WA, WAd, &WAfr );
/*
*  WB := B( IB:IB+M-1, JB+k+kbb:JB+N-1 ) * A(IA+k:IA+k+kbb-1, JA+k+kbb:JA+N-1)'
*/
                  PB_COutV( TYPE, COLUMN, INIT, M, ktmp, Bd0, kbb, &WB, WBd,
                            &WBfr, &WBsum );
                  Bnq = PB_Cnumroc( ktmp-kbb, k+kbb, Binb1, Bnb, mycol, Bcol,
                                    npcol );
                  if( ( Bmp0 > 0 ) && ( Bnq > 0 ) )
                  {
                     WAoff = PB_Cnumroc( kbb, 0, WAd[INB_], WAd[NB_], mycol,
                                         WAd[CSRC_], npcol );
                     gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRANSA ), &Bmp0, &kbb,
                           &Bnq, one, Mptr( Bptr0, 0, Bnq0-Bnq, Bld, size ),
                           &Bld, Mptr( WA, 0, WAoff, WAd[LLD_], size ),
                           &WAd[LLD_], zero, WB, &WBd[LLD_] );
                  }
                  if( WBsum )
                  {
                     WBd[CSRC_] = PB_Cindxg2p( k + kbb - 1, Binb1, Bnb, Bcol,
                                               Bcol, npcol );
                     if( Bmp0 > 0 )
                        gsum2d( ctxt, ROW, &top, Bmp0, kbb, WB, WBd[LLD_],
                                myrow, WBd[CSRC_] );
                  }
/*
*  Add WB to  B( IB:IB+M-1, JB+k:JB+k+kbb-1 ) and solve it with
*  A( IA+k:IA+k+kbb-1, JA+k:JA+k+kbb-1 )
*/
                  PB_CptrsmAB1( TYPE, SIDE, UPLO, TRANSA, DIAG, M, kbb, ALPHA,
                                WA, 0, 0, WAd, B, IB, JB+k, DESCB, WB, WBd );
                  if( WBfr ) free( WB   );
                  if( WAfr ) free( WA   );
                  if( Afr  ) free( Aptr );
               }
            }
            else
            {
/*
*  Accumulate A( IA:IA+Binb1-1, JA:JA+Binb1-1 )
*/
               PB_CGatherV( TYPE, REUSE, FORWARD, Binb1, Binb1, A, IA, JA,
                            DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA:IA+Binb1-1, JA:JA+Binb1-1 ) over B(IB:IB+M-1, JB:JB+Binb1-1)
*/
               PB_Cdescset( Bd0, M, Binb1, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                            ctxt, Bld );
               PB_CInV( TYPE, NOCONJG, ROW, M, Binb1, Bd0, Binb1, Aptr, 0, 0,
                        DBUFA, ROW, &WA, WAd, &WAfr );
/*
*  Solve B( IB:IB+M-1, JB:JB+Binb1-1 )
*/
               if( ( ( Bcol < 0 ) || ( mycol == Bcol ) ) && ( Bmp0 > 0 ) )
                  TYPE->Ftrsm( C2F_CHAR( SIDE   ), C2F_CHAR( UPLO ),
                               C2F_CHAR( TRANSA ), C2F_CHAR( DIAG ),
                               &Bmp0, &Binb1, ALPHA, WA, &WAd[LLD_],
                               Bptr0, &Bld );
               if( WAfr ) free( WA   );
               if( Afr  ) free( Aptr );
/*
*  Update and solve remaining columns of sub( B )
*/
               for( k = Binb1; k < N; k += kb )
               {
                  kbb = N - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA+k:IA+ktmp-1, JA:JA+ktmp-1 )
*/
                  PB_CGatherV( TYPE, REUSE, FORWARD, kbb, ktmp, A, IA+k, JA,
                               DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+ktmp-1, JA:JA+ktmp-1 ) over B( IB:IB+M-1, JB:JB+ktmp-1 )
*/
                  PB_Cdescset( Bd0, M, ktmp, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                               ctxt, Bld );
                  PB_CInV( TYPE, NOCONJG, ROW, M, ktmp, Bd0, kbb, Aptr, 0, 0,
                           DBUFA, ROW, &WA, WAd, &WAfr );
/*
*  WB := B( IB:IB+M-1, JB:JB+k-1 ) * A( IA+k:IA+ktmp-1, JA:JA+k-1 )'
*/
                  PB_COutV( TYPE, COLUMN, INIT, M, ktmp, Bd0, kbb, &WB, WBd,
                            &WBfr, &WBsum );
                  Bnq = PB_Cnumroc( k, 0, Binb1, Bnb, mycol, Bcol, npcol );
                  if( ( Bmp0 > 0 ) && ( Bnq > 0 ) )
                     gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRANSA ), &Bmp0, &kbb,
                           &Bnq, one, Bptr0, &Bld, WA, &WAd[LLD_], zero, WB,
                           &WBd[LLD_] );
                  if( WBsum )
                  {
                     WBd[CSRC_] = PB_Cindxg2p( k, Binb1, Bnb, Bcol, Bcol,
                                               npcol );
                     if( Bmp0 > 0 )
                        gsum2d( ctxt, ROW, &top, Bmp0, kbb, WB, WBd[LLD_],
                                myrow, WBd[CSRC_] );
                  }
/*
*  Add WB to B( IB:IB+M-1, JB+k:JB+ktmp-1 ) and solve it with
*  A( IA+k:IA+ktmp-1, JA+k:JA+ktmp-1 )
*/
                  PB_CptrsmAB1( TYPE, SIDE, UPLO, TRANSA, DIAG, M, kbb, ALPHA,
                                WA, 0, k, WAd, B, IB, JB+k, DESCB, WB, WBd );
                  if( WAfr ) free( WA   );
                  if( Afr  ) free( Aptr );
                  if( WBfr ) free( WB   );
               }
            }
         }
      }
   }
/*
*  End of PB_CptrsmAB
*/
}
