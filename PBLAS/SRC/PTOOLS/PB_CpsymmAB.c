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
void PB_CpsymmAB( PBTYP_T * TYPE, char * DIRECAB, char * CONJUG,
                  char * SIDE, char * UPLO, Int M, Int N, char * ALPHA,
                  char * A, Int IA, Int JA, Int * DESCA, char * B,
                  Int IB, Int JB, Int * DESCB, char * BETA, char * C,
                  Int IC, Int JC, Int * DESCC )
#else
void PB_CpsymmAB( TYPE, DIRECAB, CONJUG, SIDE, UPLO, M, N, ALPHA, A, IA,
                  JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * CONJUG, * DIRECAB, * SIDE, * UPLO;
   Int            IA, IB, IC, JA, JB, JC, M, N;
   char           * ALPHA, * BETA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA, * DESCB, * DESCC;
   char           * A, * B, * C;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CpsymmAB  performs one of the matrix-matrix operations
*
*     sub( C ) := alpha*sub( A )*sub( B ) + beta*sub( C ),
*
*  or
*
*     sub( C ) := alpha*sub( B )*sub( A ) + beta*sub( C ),
*
*  where
*
*     sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),
*
*     sub( A ) denotes A(IA:IA+M-1,JA:JA+M-1)  if SIDE = 'L',
*                      A(IA:IA+N-1,JA:JA+N-1)  if SIDE = 'R', and,
*
*     sub( B ) denotes B(IB:IB+M-1,JB:JB+N-1).
*
*  Alpha  and  beta  are scalars,  sub( A )  is a symmetric or Hermitian
*  submatrix and sub( B ) and sub( C ) are m by n submatrices.
*
*  This  is  the  outer-product  algorithm using the logical aggregation
*  blocking technique. The submatrix operand sub( C ) stays in place.
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
*  DIRECAB (global input) pointer to CHAR
*          On entry, DIRECAB specifies  the direction in which the  rows
*          or columns of sub( A ) and sub( B ) should be looped  over as
*          follows:
*             DIRECAB = 'F' or 'f'   forward  or increasing,
*             DIRECAB = 'B' or 'b'   backward or decreasing.
*
*  CONJUG  (global input) pointer to CHAR
*          On entry, CONJUG specifies whether sub( A ) is a symmetric or
*          Hermitian submatrix operand as follows:
*             CONJUG = 'N' or 'n'    sub( A ) is symmetric,
*             CONJUG = 'Z' or 'z'    sub( A ) is Hermitian.
*
*  SIDE    (global input) pointer to CHAR
*          On entry, SIDE  specifies  whether the symmetric or Hermitian
*          submatrix sub( A ) appears on the left or right in the opera-
*          tion as follows:
*
*             SIDE = 'L' or 'l'
*                   sub( C ) := alpha*sub( A )*sub( B ) + beta*sub( C ),
*
*             SIDE = 'R' or 'r'
*                   sub( C ) := alpha*sub( B )*sub( A ) + beta*sub( C ).
*
*  UPLO    (global input) pointer to CHAR
*          On  entry,   UPLO  specifies  whether  the  local  pieces  of
*          the array  A  containing the  upper or lower triangular  part
*          of the submatrix  sub( A )  are to be referenced as follows:
*             UPLO = 'U' or 'u'   Only the local pieces corresponding to
*                                 the   upper  triangular  part  of  the
*                                 submatrix sub( A ) are referenced,
*             UPLO = 'L' or 'l'   Only the local pieces corresponding to
*                                 the   lower  triangular  part  of  the
*                                 submatrix sub( A ) are referenced.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( C ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( C ). N  must be at least zero.
*
*  ALPHA   (global input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as zero then the local entries of the arrays  A and
*          B corresponding to the entries of  the  submatrices  sub( A )
*          and sub( B ) respectively need not be set on input.
*
*  A       (local input) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least  Lc( 1, JA+M-1 )  when  SIDE = 'L' or 'l'  and is at
*          at least Lc( 1, JA+N-1 ) otherwise. Before  entry, this array
*          contains the local entries of the matrix A.
*          Before  entry  with  SIDE = 'L' or 'l', this  array  contains
*          the local entries corresponding to the entries of the  m by m
*          symmetric or Hermitian submatrix  sub( A ),  such  that  when
*          UPLO = 'U' or 'u', this  array contains the local entries  of
*          the upper triangular part of the submatrix  sub( A ), and the
*          local entries  of  the strictly lower triangular of  sub( A )
*          are not referenced, and when  UPLO = 'L' or 'l',  this  array
*          contains  the local entries of the  lower triangular part  of
*          the symmetric or Hermitian submatrix sub( A ), and  the local
*          entries of the strictly upper triangular of sub( A ) are  not
*          referenced.
*          Before  entry  with  SIDE = 'R' or 'r', this  array  contains
*          the local entries corresponding to the entries of the  n by n
*          symmetric or Hermitian submatrix  sub( A ),  such  that  when
*          UPLO = 'U' or 'u', this  array contains the local entries  of
*          the upper triangular part of the submatrix sub( A ), and  the
*          local entries  of  the strictly lower triangular of  sub( A )
*          are not referenced, and when  UPLO = 'L' or 'l',  this  array
*          contains  the local entries of the  lower triangular part  of
*          the  symmetric or Hermitian submatrix sub( A ), and the local
*          entries of the strictly upper triangular of sub( A ) are  not
*          referenced.
*          Note that the  imaginary parts  of the local entries  corres-
*          ponding to the  diagonal elements  of  sub( A )  need not  be
*          set and assumed to be zero.
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
*  B       (local input) pointer to CHAR
*          On entry, B is an array of dimension (LLD_B, Kb), where Kb is
*          at least Lc( 1, JB+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix B.
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
*  BETA    (global input) pointer to CHAR
*          On entry,  BETA  specifies the scalar  beta.   When  BETA  is
*          supplied  as  zero  then  the  local entries of  the array  C
*          corresponding to  the  entries of the submatrix sub( C ) need
*          not be set on input.
*
*  C       (local input/local output) pointer to CHAR
*          On entry, C is an array of dimension (LLD_C, Kc), where Kc is
*          at least Lc( 1, JC+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix C.
*          On exit, the entries of this array corresponding to the local
*          entries  of the submatrix  sub( C )  are  overwritten by  the
*          local entries of the m by n updated submatrix.
*
*  IC      (global input) INTEGER
*          On entry, IC  specifies C's global row index, which points to
*          the beginning of the submatrix sub( C ).
*
*  JC      (global input) INTEGER
*          On entry, JC  specifies C's global column index, which points
*          to the beginning of the submatrix sub( C ).
*
*  DESCC   (global and local input) INTEGER array
*          On entry, DESCC  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix C.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           GatherDir, ScatterDir, * one, top, tran, * zero;
   Int            Afr, An, Bcol, Bcurcol, Bcurimb1, Bcurinb1, Bcurrow, Bfr, Bii,
                  Bimb, Bimb1, Binb, Binb1, Bjj, Bld, Bmb, Bmp, Bmp0, Bnb, Bnq,
                  Bnq0, Brow, Ccol, Ccurcol, Ccurimb1, Ccurinb1, Ccurrow, Cii,
                  Cimb, Cimb1, Cinb, Cinb1, Cjj, Cld, Cmb, Cmp, Cmp0, Cnb, Cnq,
                  Cnq0, Crow, WABfr, WACfr, WBCfr, WBCsum, conjg, ctxt, fwd, k,
                  kb, kbb, kend, kstart, kstep, ktmp, lside, mycol, myrow,
                  npcol, nprow, size, upper;
   GEMM_T         gemm;
   GSUM2D_T       gsum2d;
/*
*  .. Local Arrays ..
*/
   Int            Bd0 [DLEN_], Cd0 [DLEN_], DBUFA[DLEN_], DBUFB[DLEN_],
                  WABd[DLEN_], WACd[DLEN_], WBCd [DLEN_];
   char           * Aptr = NULL, * Bptr = NULL, * Bptr0 = NULL, * Cptr0 = NULL,
                  * WAB  = NULL, * WAC  = NULL, * WBC   = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  sub( C ) = beta * sub( C )
*/
   PB_Cplascal( TYPE, ALL, NOCONJG, M, N, BETA, C, IC, JC, DESCC );
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   An    = ( ( lside = ( Mupcase( SIDE[0] ) == CLEFT ) ) ? M : N );
   upper = ( Mupcase( UPLO[0] ) == CUPPER );
   tran  = ( ( conjg = ( Mupcase( CONJUG[0] ) == CCONJG ) ) ? CCOTRAN : CTRAN );

   size  = TYPE->size;  one    = TYPE->one;     zero   = TYPE->zero;
   gemm  = TYPE->Fgemm; gsum2d = TYPE->Cgsum2d;
/*
*  Figure out the loop bounds accordingly to DIRECAB
*/
   kb    = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );
   if( ( fwd = ( Mupcase( DIRECAB[0] ) == CFORWARD ) ) != 0 )
   {
      kstart = 0; kend = ( ( An - 1 ) / kb + 1 ) * kb; kstep =  kb;
      GatherDir = ScatterDir = CFORWARD;
   }
   else
   {
      kstart = ( ( An - 1 ) / kb ) * kb; kend = kstep = -kb;
      GatherDir = ScatterDir = CBACKWARD;
   }
/*
*  Compute local information for sub( B ) and sub( C )
*/
   PB_Cinfog2l( IB, JB, DESCB, nprow, npcol, myrow, mycol, &Bii, &Bjj,
                &Brow, &Bcol );
   Bimb  = DESCB[IMB_]; Binb = DESCB[INB_];
   Bmb   = DESCB[MB_ ]; Bnb  = DESCB[NB_ ]; Bld = DESCB[LLD_];
   Bimb1 = PB_Cfirstnb( M, IB, Bimb, Bmb );
   Bmp0  = PB_Cnumroc( M, 0, Bimb1, Bmb, myrow, Brow, nprow );
   Binb1 = PB_Cfirstnb( N, JB, Binb, Bnb );
   Bnq0  = PB_Cnumroc( N, 0, Binb1, Bnb, mycol, Bcol, npcol );
   if( ( Bmp0 > 0 ) && ( Bnq0 > 0 ) ) Bptr0 = Mptr( B, Bii, Bjj, Bld, size );

   PB_Cinfog2l( IC, JC, DESCC, nprow, npcol, myrow, mycol, &Cii, &Cjj,
                &Crow, &Ccol );
   Cimb  = DESCC[IMB_]; Cinb = DESCC[INB_];
   Cmb   = DESCC[MB_ ]; Cnb  = DESCC[NB_ ]; Cld = DESCC[LLD_];
   Cimb1 = PB_Cfirstnb( M, IC, Cimb, Cmb );
   Cmp0  = PB_Cnumroc( M, 0, Cimb1, Cmb, myrow, Crow, nprow );
   Cinb1 = PB_Cfirstnb( N, JC, Cinb, Cnb );
   Cnq0  = PB_Cnumroc( N, 0, Cinb1, Cnb, mycol, Ccol, npcol );
   if( ( Cmp0 > 0 ) && ( Cnq0 > 0 ) ) Cptr0 = Mptr( C, Cii, Cjj, Cld, size );

   if( lside )
   {
      top = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );

      if( upper )
      {
         for( k = kstart; k != kend; k += kstep )
         {
            kbb = An - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA:IA+k+kbb-1, JA+k:JA+k+kbb-1 )
*/
            PB_CGatherV( TYPE, ALLOCATE, &GatherDir, ktmp, kbb, A, IA, JA+k,
                         DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA:IA+k+kbb-1, JA+k:JA+k+kbb-1 ) over
*  C( IC:IC+k+kbb-1, JC:JC+N-1 ) -> WAC
*/
            PB_Cdescset( Cd0, ktmp, N, Cimb1, Cinb1, Cmb, Cnb, Crow, Ccol,
                         ctxt, Cld );
            PB_CInV( TYPE, NOCONJG, COLUMN, ktmp, N, Cd0, kbb, Aptr, 0, 0,
                     DBUFA, COLUMN, &WAC, WACd, &WACfr );
/*
*  Zero lower triangle of WAC( k:k+kbb-1, 0:kbb-1 )
*/
            if( conjg )
               PB_Cplapad( TYPE, LOWER, CONJUG,  kbb,   kbb,   zero,
                           zero, WAC, k,   0, WACd );
            else if( kbb > 1 )
               PB_Cplapad( TYPE, LOWER, NOCONJG, kbb-1, kbb-1, zero,
                           zero, WAC, k+1, 0, WACd );
/*
*  Accumulate B( IB+k:IB+k+kbb-1, JB:JB+N-1 )
*/
            PB_CGatherV( TYPE, REUSE, &GatherDir, kbb, N, B, IB+k, JB, DESCB,
                         ROW, &Bptr, DBUFB, &Bfr );
/*
*  Replicate B( IB+k:IB+k+kbb-1, JB:JB+N-1 ) over C( IC:IC+k+kbb-1, JC:JC+N-1 )
*/
            PB_CInV( TYPE, NOCONJG, ROW, ktmp, N, Cd0, kbb, Bptr, 0, 0, DBUFB,
                     ROW, &WBC, WBCd, &WBCfr );
/*
*  C( IC:IC+k+kbb-1, JC:JC+N-1 ) += ALPHA * WAC * WBC
*/
            Cmp = PB_Cnumroc( ktmp, 0, Cimb1, Cmb, myrow, Crow, nprow );
            if( ( Cmp > 0 ) && ( Cnq0 > 0 ) )
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp, &Cnq0, &kbb,
                     ALPHA, WAC, &WACd[LLD_], WBC, &WBCd[LLD_], one, Cptr0,
                     &Cld );
            if( WBCfr ) free( WBC  );
            if( Bfr   ) free( Bptr );
/*
*  Replicate WAC = A( IA:IA+k+kbb-1, JA+k:JA+k+kbb-1 ) over
*  B( IB:IB+k+kbb-1, JB:JB+N-1 ) -> WAB
*/
            PB_Cdescset( Bd0, ktmp, N, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol,
                         ctxt, Bld );
            PB_CInV( TYPE, NOCONJG, COLUMN, ktmp, N, Bd0, kbb, WAC,  0, 0, WACd,
                     COLUMN, &WAB, WABd, &WABfr );
/*
*  Zero lower triangle of WAB( k:k+kbb-1, 0:kbb-1 )
*/
            PB_Cplapad( TYPE, LOWER, NOCONJG, kbb, kbb, zero, zero, WAB, k, 0,
                        WABd );
/*
*  WBC := ALPHA*A(IA:IA+k+kbb-1, JA+k:JA+k+kbb-1)'*B( IB:IB+k+kbb-1, JB:JB+N-1 )
*/
            PB_COutV( TYPE, ROW, INIT, ktmp, N, Bd0, kbb, &WBC, WBCd, &WBCfr,
                      &WBCsum );
            Bmp = PB_Cnumroc( ktmp, 0, Bimb1, Bmb, myrow, Brow, nprow );
            if( ( Bnq0 > 0 ) && ( Bmp > 0 ) )
               gemm( C2F_CHAR( &tran ), C2F_CHAR( NOTRAN ), &kbb, &Bnq0, &Bmp,
                     ALPHA, WAB, &WABd[LLD_], Bptr0, &Bld, zero, WBC,
                     &WBCd[LLD_] );
            if( WABfr ) free( WAB  );
            if( WACfr ) free( WAC  );
            if( Afr   ) free( Aptr );

            if( WBCsum )
            {
               WBCd[RSRC_] = PB_Cindxg2p( ( fwd ? k : k + kbb - 1 ), Cimb1,
                                          Cmb, Crow, Crow, nprow );
               if( Bnq0 > 0 )
                  gsum2d( ctxt, COLUMN, &top, kbb, Bnq0, WBC, WBCd[LLD_],
                          WBCd[RSRC_], mycol );
            }
/*
*  C( IC+k:IC+k+kbb-1, JC:JC+N-1 ) := C( IC+k:IC+k+kbb-1, JC:JC+N-1 ) + WBC
*/
            PB_CScatterV( TYPE, &ScatterDir, kbb, N, WBC, 0, 0, WBCd, ROW, one,
                          C, IC+k, JC, DESCC, ROW );
            if( WBCfr ) free( WBC  );
         }
      }
      else
      {
         for( k = kstart; k != kend; k += kstep )
         {
            ktmp = An - k; kbb = MIN( ktmp, kb );
/*
*  Accumulate A( IA+k:IA+M-1, JA+k:JA+k+kbb-1 )
*/
            PB_CGatherV( TYPE, ALLOCATE, &GatherDir, ktmp, kbb, A, IA+k, JA+k,
                         DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+M-1, JA+k:JA+k+kbb-1 ) over C( IC+k:IC+M-1, JC:JC+N-1 )
*/
            Ccurimb1 = PB_Cfirstnb( ktmp, IC+k, Cimb, Cmb );
            Ccurrow  = PB_Cindxg2p( k, Cimb1, Cmb, Crow, Crow, nprow );
            PB_Cdescset( Cd0, ktmp, N, Ccurimb1, Cinb1, Cmb, Cnb, Ccurrow,
                         Ccol, ctxt, Cld );
            PB_CInV( TYPE, NOCONJG, COLUMN, ktmp, N, Cd0, kbb, Aptr, 0, 0,
                     DBUFA, COLUMN, &WAC, WACd, &WACfr );
/*
*  Zero upper triangle of WAC( 0:kbb-1, 0:kbb-1 )
*/
            if( conjg )
               PB_Cplapad( TYPE, UPPER, CONJUG,  kbb,   kbb,   zero, zero, WAC,
                           0, 0, WACd );
            else if( kbb > 1 )
               PB_Cplapad( TYPE, UPPER, NOCONJG, kbb-1, kbb-1, zero, zero, WAC,
                           0, 1, WACd );
/*
*  Accumulate B( IB+k:IB+k+kbb-1, JB:JB+N-1 )
*/
            PB_CGatherV( TYPE, REUSE, &GatherDir, kbb, N, B, IB+k, JB, DESCB,
                         ROW, &Bptr, DBUFB, &Bfr );
/*
*  Replicate B( IB+k:IB+k+kbb-1, JB:JB+N-1 ) over C( IC+k:IC+M-1, JC:JC+N-1 )
*/
            PB_CInV( TYPE, NOCONJG, ROW, ktmp, N, Cd0, kbb, Bptr, 0, 0, DBUFB,
                     ROW, &WBC, WBCd, &WBCfr );
/*
*  C( IC+k:IC+M-1, JC:JC+N-1 ) += ALPHA * WAC * WBC
*/
            Cmp = PB_Cnumroc( ktmp, k, Cimb1, Cmb, myrow, Crow, nprow );
            if( ( Cmp > 0 ) && ( Cnq0 > 0 ) )
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp, &Cnq0, &kbb,
                     ALPHA, WAC, &WACd[LLD_], WBC, &WBCd[LLD_], one,
                     Mptr( Cptr0, Cmp0-Cmp, 0, Cld, size ), &Cld );
            if( WBCfr ) free( WBC  );
            if( Bfr   ) free( Bptr );
/*
*  Replicate WAC = A( IA+k:IA+M-1, JA+k:JA+k+kbb-1 ) over
*  B( IB+k:IB+M-1, JB:JB+N-1 ) -> WAB
*/
            Bcurimb1 = PB_Cfirstnb( ktmp, IB+k, Bimb, Bmb );
            Bcurrow  = PB_Cindxg2p( k, Bimb1, Bmb, Brow, Brow, nprow );
            PB_Cdescset( Bd0, ktmp, N, Bcurimb1, Binb1, Bmb, Bnb, Bcurrow,
                         Bcol, ctxt, Bld );
            PB_CInV( TYPE, NOCONJG, COLUMN, ktmp, N, Bd0, kbb, WAC, 0, 0, WACd,
                     COLUMN, &WAB, WABd, &WABfr );
/*
*  Zero upper triangle of WAB( 0:kbb-1, 0:kbb-1 )
*/
            PB_Cplapad( TYPE, UPPER, NOCONJG, kbb, kbb, zero, zero, WAB, 0, 0,
                        WABd );
/*
*  WBC := ALPHA*A( IA+k:IA+M-1, JA+k:JA+k+kbb-1 )'*B( IB+k:IB+M-1, JB:JB+N-1 )
*/
            PB_COutV( TYPE, ROW, INIT, ktmp, N, Bd0, kbb, &WBC, WBCd, &WBCfr,
                      &WBCsum );
            Bmp = PB_Cnumroc( ktmp, k, Bimb1, Bmb, myrow, Brow, nprow );
            if( ( Bnq0 > 0 ) && ( Bmp > 0 ) )
               gemm( C2F_CHAR( &tran ), C2F_CHAR( NOTRAN ), &kbb, &Bnq0, &Bmp,
                     ALPHA, WAB, &WABd[LLD_], Mptr( Bptr0, Bmp0-Bmp, 0, Bld,
                     size ), &Bld, zero, WBC, &WBCd[LLD_] );
            if( WABfr ) free( WAB  );
            if( WACfr ) free( WAC  );
            if( Afr   ) free( Aptr );

            if( WBCsum )
            {
               WBCd[RSRC_] = PB_Cindxg2p( ( fwd ? k : k + kbb - 1 ), Cimb1,
                                          Cmb, Crow, Crow, nprow );
               if( Bnq0 > 0 )
                  gsum2d( ctxt, COLUMN, &top, kbb, Bnq0, WBC, WBCd[LLD_],
                          WBCd[RSRC_], mycol );
            }
/*
*  C( IC+k:IC+k+kbb-1, JC:JC+N-1 ) := C( IC+k:IC+k+kbb-1, JC:JC+N-1 ) + WBC
*/
            PB_CScatterV( TYPE, &ScatterDir, kbb, N, WBC, 0, 0, WBCd, ROW, one,
                          C, IC+k, JC, DESCC, ROW );
            if( WBCfr ) free( WBC  );
         }
      }
   }
   else
   {
      top = *PB_Ctop( &ctxt, COMBINE, ROW, TOP_GET );

      if( upper )
      {
         for( k = kstart; k != kend; k += kstep )
         {
            ktmp = An - k; kbb = MIN( ktmp, kb );
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 )
*/
            PB_CGatherV( TYPE, ALLOCATE, &GatherDir, kbb, ktmp, A, IA+k, JA+k,
                         DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 ) over C( IC:IC+M-1, JC+k:JC+N-1 )
*/
            Ccurinb1 = PB_Cfirstnb( ktmp, JC+k, Cinb, Cnb );
            Ccurcol  = PB_Cindxg2p( k, Cinb1, Cnb, Ccol, Ccol, npcol );
            PB_Cdescset( Cd0, M, ktmp, Cimb1, Ccurinb1, Cmb, Cnb, Crow, Ccurcol,
                         ctxt, Cld );
            PB_CInV( TYPE, NOCONJG, ROW, M, ktmp, Cd0, kbb, Aptr, 0, 0, DBUFA,
                     ROW, &WAC, WACd, &WACfr );
/*
*  Zero lower triangle of WAC( 0:kbb-1, 0:kbb-1 )
*/
            if( conjg )
               PB_Cplapad( TYPE, LOWER, CONJUG,  kbb,   kbb,   zero,
                           zero, WAC, 0, 0, WACd );
            else if( kbb > 1 )
               PB_Cplapad( TYPE, LOWER, NOCONJG, kbb-1, kbb-1, zero,
                           zero, WAC, 1, 0, WACd );
/*
*  Accumulate B( IB:IB+M-1, JB+k:JB+k+kbb-1 )
*/
            PB_CGatherV( TYPE, REUSE, &GatherDir, M, kbb, B, IB, JB+k, DESCB,
                         COLUMN, &Bptr, DBUFB, &Bfr );
/*
*  Replicate B( IB:IB+M-1, JB+k:JB+k+kbb-1 ) over C( IC:IC+M-1, JC+k:JC+N-1 )
*/
            PB_CInV( TYPE, NOCONJG, COLUMN, M, ktmp, Cd0, kbb, Bptr, 0, 0,
                     DBUFB, COLUMN, &WBC, WBCd, &WBCfr );
/*
*  C( IC:IC+M-1, JC+k:JC+N-1 ) += ALPHA * WBC * WAC
*/
            Cnq = PB_Cnumroc( ktmp, k, Cinb1, Cnb, mycol, Ccol, npcol );
            if( ( Cmp0 > 0 ) && ( Cnq > 0 ) )
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp0, &Cnq, &kbb,
                     ALPHA, WBC, &WBCd[LLD_], WAC, &WACd[LLD_], one,
                     Mptr( Cptr0, 0, Cnq0-Cnq, Cld, size ), &Cld );
            if( WBCfr ) free( WBC  );
            if( Bfr   ) free( Bptr );
/*
*  Replicate WAC = A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 ) over
*  B( IB:IB+M-1, JB+k:JB+N-1 ) -> WAB
*/
            Bcurinb1 = PB_Cfirstnb( ktmp, JB+k, Binb, Bnb );
            Bcurcol  = PB_Cindxg2p( k, Binb1, Bnb, Bcol, Bcol, npcol );
            PB_Cdescset( Bd0, M, ktmp, Bimb1, Bcurinb1, Bmb, Bnb, Brow, Bcurcol,
                         ctxt, Bld );
            PB_CInV( TYPE, NOCONJG, ROW, M, ktmp, Bd0, kbb, WAC, 0, 0, WACd,
                     ROW, &WAB, WABd, &WABfr );
/*
*  Zero lower triangle of WAB( 0:kbb-1, 0:kbb-1 )
*/
            PB_Cplapad( TYPE, LOWER, NOCONJG, kbb, kbb, zero, zero, WAB, 0, 0,
                        WABd );
/*
*  WBC := ALPHA*B( IB:IB+M-1, JB+k:JB+N-1 )*A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 )'
*/
            PB_COutV( TYPE, COLUMN, INIT, M, ktmp, Bd0, kbb, &WBC, WBCd,
                      &WBCfr, &WBCsum );
            Bnq = PB_Cnumroc( ktmp, k, Binb1, Bnb, mycol, Bcol, npcol );
            if( ( Bmp0 > 0 ) && ( Bnq > 0 ) )
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( &tran ), &Bmp0, &kbb, &Bnq,
                     ALPHA, Mptr( Bptr0, 0, Bnq0-Bnq, Bld, size ), &Bld, WAB,
                     &WABd[LLD_], zero, WBC, &WBCd[LLD_] );
            if( WABfr ) free( WAB  );
            if( WACfr ) free( WAC  );
            if( Afr   ) free( Aptr );

            if( WBCsum )
            {
               WBCd[CSRC_] = PB_Cindxg2p( ( fwd ? k : k + kbb - 1 ), Cinb1,
                                          Cnb, Ccol, Ccol, npcol );
               if( Bmp0 > 0 )
                  gsum2d( ctxt, ROW, &top, Bmp0, kbb, WBC, WBCd[LLD_], myrow,
                          WBCd[CSRC_] );
            }
/*
*  C( IC:IC+M-1, JC+k:JC+k+kbb-1 ) := C( IC:IC+M-1, JC+k:JC+k+kbb-1 ) + WBC
*/
            PB_CScatterV( TYPE, &ScatterDir, M, kbb, WBC, 0, 0, WBCd, COLUMN,
                          one, C, IC, JC+k, DESCC, COLUMN );
            if( WBCfr ) free( WBC  );
         }
      }
      else
      {
         for( k = kstart; k != kend; k += kstep )
         {
            kbb = An - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 )
*/
            PB_CGatherV( TYPE, ALLOCATE, &GatherDir, kbb, ktmp, A, IA+k, JA,
                         DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 ) over
*  C( IC:IC+M-1, JC:JC+k+kbb-1 ) -> WAC
*/
            PB_Cdescset( Cd0, M, ktmp, Cimb1, Cinb1, Cmb, Cnb, Crow, Ccol, ctxt,
                         Cld );
            PB_CInV( TYPE, NOCONJG, ROW, M, ktmp, Cd0, kbb, Aptr, 0, 0, DBUFA,
                     ROW, &WAC, WACd, &WACfr );
/*
*  Zero upper triangle of WAC( 0:kbb-1, k:k+kbb-1 )
*/
            if( conjg )
               PB_Cplapad( TYPE, UPPER, CONJUG,  kbb,   kbb,   zero, zero, WAC,
                           0, k,  WACd );
            else if( kbb > 1 )
               PB_Cplapad( TYPE, UPPER, NOCONJG, kbb-1, kbb-1, zero, zero, WAC,
                           0, k+1, WACd );
/*
*  Accumulate B( IB:IB+M-1, JB+k:JB+k+kbb-1 )
*/
            PB_CGatherV( TYPE, REUSE, &GatherDir, M, kbb, B, IB, JB+k, DESCB,
                         COLUMN, &Bptr, DBUFB, &Bfr );
/*
*  Replicate B( IB:IB+M-1, JB+k:JB+k+kbb-1 ) over C( IC:IC+M-1, JC:JC+k+kbb-1 )
*/
            PB_CInV( TYPE, NOCONJG, COLUMN, M, ktmp, Cd0, kbb, Bptr, 0, 0,
                     DBUFB, COLUMN, &WBC, WBCd, &WBCfr );
/*
*  C( IC:IC+M-1, JC:JC+k+kbb-1 ) += ALPHA * WBC * WAC
*/
            Cnq = PB_Cnumroc( ktmp, 0, Cinb1, Cnb, mycol, Ccol, npcol );
            if( ( Cmp0 > 0 ) && ( Cnq > 0 ) )
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp0, &Cnq, &kbb,
                     ALPHA, WBC, &WBCd[LLD_], WAC, &WACd[LLD_], one, Cptr0,
                     &Cld );
            if( WBCfr ) free( WBC  );
            if( Bfr   ) free( Bptr );
/*
*  Replicate WAC = A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 ) over
*  B( IB:IB+M-1, JB:JB+k+kbb-1 ) -> WAB
*/
            PB_Cdescset( Bd0, M, ktmp, Bimb1, Binb1, Bmb, Bnb, Brow, Bcol, ctxt,
                         Bld );
            PB_CInV( TYPE, NOCONJG, ROW, M, ktmp, Bd0, kbb, WAC, 0, 0, WACd,
                     ROW, &WAB, WABd, &WABfr );
/*
*  Zero upper triangle of WAB( 0:kbb-1, k:k+kbb-1 )
*/
            PB_Cplapad( TYPE, UPPER, NOCONJG, kbb, kbb, zero, zero, WAB, 0, k,
                        WABd );
/*
*  WBC := ALPHA*B( IB:IB+M-1, JB:JB+k+kbb-1 )*A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 )'
*/
            PB_COutV( TYPE, COLUMN, INIT, M, ktmp, Bd0, kbb, &WBC, WBCd, &WBCfr,
                      &WBCsum );
            Bnq = PB_Cnumroc( ktmp, 0, Binb1, Bnb, mycol, Bcol, npcol );
            if( ( Bmp0 > 0 ) && ( Bnq > 0 ) )
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( &tran ), &Bmp0, &kbb, &Bnq,
                     ALPHA, Bptr0, &Bld, WAB, &WABd[LLD_], zero, WBC,
                     &WBCd[LLD_] );
            if( WABfr ) free( WAB  );
            if( WACfr ) free( WAC  );
            if( Afr   ) free( Aptr );

            if( WBCsum )
            {
               WBCd[CSRC_] = PB_Cindxg2p( ( fwd ? k : k + kbb - 1 ), Cinb1,
                                          Cnb, Ccol, Ccol, npcol );
               if( Bmp0 > 0 )
                  gsum2d( ctxt, ROW, &top, Bmp0, kbb, WBC, WBCd[LLD_], myrow,
                          WBCd[CSRC_] );
            }
/*
*  C( IC:IC+M-1, JC+k:JC+k+kbb-1 ) := C( IC:IC+M-1, JC+k:JC+k+kbb-1 ) + WBC
*/
            PB_CScatterV( TYPE, &ScatterDir, M, kbb, WBC, 0, 0, WBCd, COLUMN,
                          one, C, IC, JC+k, DESCC, COLUMN );
            if( WBCfr ) free( WBC  );
         }
      }
   }
/*
*  End of PB_CpsymmAB
*/
}
