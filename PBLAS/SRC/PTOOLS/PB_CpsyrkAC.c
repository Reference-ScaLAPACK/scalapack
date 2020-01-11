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
void PB_CpsyrkAC( PBTYP_T * TYPE, char * DIRECA, char * CONJUG,
                  char * UPLO, char * TRANS, Int N, Int K, char * ALPHA,
                  char * A, Int IA, Int JA, Int * DESCA, char * BETA,
                  char * C, Int IC, Int JC, Int * DESCC )
#else
void PB_CpsyrkAC( TYPE, DIRECA, CONJUG, UPLO, TRANS, N, K, ALPHA, A, IA,
                 JA, DESCA, BETA, C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * CONJUG, * DIRECA, * TRANS, * UPLO;
   Int            IA, IC, JA, JC, K, N;
   char           * ALPHA, * BETA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA, * DESCC;
   char           * A, * C;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CpsyrkAC performs one of the following symmetric or Hermitian rank
*  k operations
*
*     sub( C ) := alpha*sub( A )*sub( A )' + beta*sub( C ),
*  or
*     sub( C ) := alpha*sub( A )*conjg( sub( A )' ) + beta*sub( C ),
*  or
*     sub( C ) := alpha*sub( A )'*sub( A ) + beta*sub( C ),
*  or
*     sub( C ) := alpha*conjg( sub( A )' )*sub( A ) + beta*sub( C ),
*
*  where
*
*     sub( C ) denotes C(IC:IC+N-1,JC:JC+N-1), and,
*
*     sub( A ) denotes A(IA:IA+N-1,JA:JA+K-1)  if TRANS = 'N',
*                      A(IA:IA+K-1,JA:JA+N-1)  otherwise.
*
*  Alpha  and   beta  are  scalars,  sub( C )  is  an  n by n  symmetric
*  or Hermitian submatrix and  sub( A )  is  an  n by k submatrix in the
*  first case and a k by n submatrix in the second case.
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
*  DIRECA  (global input) pointer to CHAR
*          On entry, DIRECA  specifies  the direction in which the  rows
*          or columns of sub( A ) and sub( C ) should be looped  over as
*          follows:
*             DIRECA = 'F' or 'f'   forward  or increasing,
*             DIRECA = 'B' or 'b'   backward or decreasing.
*
*  CONJUG  (global input) pointer to CHAR
*          On entry, CONJUG specifies whether sub( C ) is a symmetric or
*          Hermitian submatrix operand as follows:
*             CONJUG = 'N' or 'n'    sub( C ) is symmetric,
*             CONJUG = 'Z' or 'z'    sub( C ) is Hermitian.
*
*  UPLO    (global input) pointer to CHAR
*          On  entry,   UPLO  specifies  whether  the  local  pieces  of
*          the array  C  containing the  upper or lower triangular  part
*          of the submatrix  sub( C )  are to be referenced as follows:
*             UPLO = 'U' or 'u'   Only the local pieces corresponding to
*                                 the   upper  triangular  part  of  the
*                                 submatrix sub( C ) are referenced,
*             UPLO = 'L' or 'l'   Only the local pieces corresponding to
*                                 the   lower  triangular  part  of  the
*                                 submatrix sub( C ) are referenced.
*
*  TRANS   (global input) pointer to CHAR
*          On entry,  TRANS  specifies the  operation to be performed as
*          follows:
*
*             TRANS = 'N' or 'n'
*                  sub( C ) := alpha*sub( A )*sub( A )' + beta*sub( C ),
*             or
*                  sub( C ) := alpha*sub( A )*sub( A )' + beta*sub( C ),
*             or
*                  sub( C ) := alpha*sub( A )*conjg( sub( A )' ) +
*                              beta*sub( C ),
*
*             TRANS = 'T' or 't'
*                  sub( C ) := alpha*sub( A )'*sub( A ) + beta*sub( C ),
*             or
*                  sub( C ) := alpha*sub( A )'*sub( A ) + beta*sub( C ),
*
*             TRANS = 'C' or 'c'
*                  sub( C ) := alpha*sub( A )'*sub( A ) + beta*sub( C ),
*             or
*                  sub( C ) := alpha*conjg( sub( A )' )*sub( A ) +
*                              beta*sub( C ).
*
*  N       (global input) INTEGER
*          On entry,  N specifies the order of the  submatrix  sub( C ).
*          N must be at least zero.
*
*  K       (global input) INTEGER
*          On entry, with TRANS = 'N' or 'n',  K specifies the number of
*          columns  of  the submatrix  sub( A ), and with TRANS = 'T' or
*          't' or 'C' or 'c', K specifies the number of rows of the sub-
*          matrix sub( A ). K  must  be at least zero.
*
*  ALPHA   (global input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of  the array  A
*          corresponding to the entries of the submatrix  sub( A )  need
*          not be set on input.
*
*  A       (local input) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least  Lc( 1, JA+K-1 ) when  TRANS = 'N' or 'n', and is at
*          least Lc( 1, JA+N-1 ) otherwise.  Before  entry,  this  array
*          contains the local entries of the matrix A.
*          Before entry with TRANS = 'N' or 'n', this array contains the
*          local entries corresponding to the entries of the n by k sub-
*          matrix sub( A ), otherwise the local entries corresponding to
*          the entries of the k by n submatrix sub( A ).
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
*  BETA    (global input) pointer to CHAR
*          On entry,  BETA  specifies the scalar  beta.   When  BETA  is
*          supplied  as  zero  then  the  local entries of  the array  C
*          corresponding to the entries of the submatrix  sub( C )  need
*          not be set on input.
*
*  C       (local input/local output) pointer to CHAR
*          On entry, C is an array of dimension (LLD_C, Kc), where Kc is
*          at least Lc( 1, JC+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix C.
*          Before  entry  with  UPLO = 'U' or 'u', this  array  contains
*          the local entries corresponding to the upper triangular  part
*          of the  symmetric or Hermitian submatrix  sub( C ),  and  the
*          local entries corresponding to the  strictly lower triangular
*          of sub( C ) are not referenced. On exit, the upper triangular
*          part  of sub( C ) is overwritten by the upper triangular part
*          of the updated submatrix.
*          Before  entry  with  UPLO = 'L' or 'l', this  array  contains
*          the local entries corresponding to the lower triangular  part
*          of the  symmetric or Hermitian submatrix  sub( C ),  and  the
*          local entries corresponding to the  strictly upper triangular
*          of sub( C ) are not referenced. On exit, the lower triangular
*          part of sub( C ) is overwritten by the  lower triangular part
*          of the updated submatrix.
*          Note that the  imaginary parts  of the local entries  corres-
*          ponding to the  diagonal elements  of  sub( C )  need not  be
*          set,  they are assumed to be zero,  and on exit they are  set
*          to zero.
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
   Int            Acol, Acurcol, Acurimb1, Acurinb1, Acurrow, Afr, Aii, Aimb,
                  Aimb1, Ainb, Ainb1, Ajj, Ald, Am, Amb, Amp, Amp0, An, Anb,
                  Anq, Anq0, Arow, Ccsrc, Cimb, Cinb, Cmb, Cnb, Crsrc, WAfr,
                  WCfr, WCsum, conjg, ctxt, fwd, k, kb, kbb, kend, kstart,
                  kstep, ktmp, mycol, myrow, notran, npcol, nprow, size, upper;
   GEMM_T         gemm;
   GSUM2D_T       gsum2d;
/*
*  .. Local Arrays ..
*/
   Int            Ad0[DLEN_], DBUFA[DLEN_], WAd[DLEN_], WCd[DLEN_];
   char           * Aptr = NULL, * Aptr0 = NULL, * WA = NULL, * WC = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  sub( C ) = beta * sub( C )
*/
   PB_Cplascal( TYPE, UPLO, CONJUG, N, N, BETA, C, IC, JC, DESCC );
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   fwd    = ( Mupcase( DIRECA[0] ) == CFORWARD );
   conjg  = ( Mupcase( CONJUG[0] ) ==   CCONJG );
   upper  = ( Mupcase( UPLO  [0] ) ==   CUPPER );
   notran = ( Mupcase( TRANS [0] ) ==  CNOTRAN );
   tran   = ( conjg ? CCOTRAN : CTRAN );

   size   = TYPE->size;    one  = TYPE->one;   zero = TYPE->zero;
   gsum2d = TYPE->Cgsum2d; gemm = TYPE->Fgemm;
/*
*  Figure out the loop bounds accordingly to DIRECA
*/
   kb = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );
   if( fwd )
   {
      kstart = 0; kend = ( ( N - 1 ) / kb + 1 ) * kb; kstep = kb;
      GatherDir = CFORWARD; ScatterDir = CBACKWARD;
   }
   else
   {
      kstart = ( ( N - 1 ) / kb ) * kb; kend = kstep = -kb;
      GatherDir = CBACKWARD; ScatterDir = CFORWARD;
   }
/*
*  Compute local information for A
*/
   if( notran ) { Am = N; An = K; } else { Am = K; An = N; }
   PB_Cinfog2l( IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aii, &Ajj,
                &Arow, &Acol );
   Aimb = DESCA[IMB_]; Ainb = DESCA[INB_];
   Amb  = DESCA[MB_ ]; Anb  = DESCA[NB_ ]; Ald = DESCA[LLD_];

   Aimb1 = PB_Cfirstnb( Am, IA, Aimb, Amb );
   Amp0  = PB_Cnumroc( Am, 0, Aimb1, Amb, myrow, Arow, nprow );
   Ainb1 = PB_Cfirstnb( An, JA, Ainb, Anb );
   Anq0  = PB_Cnumroc( An, 0, Ainb1, Anb, mycol, Acol, npcol );
   if( ( Amp0 > 0 ) && ( Anq0 > 0 ) ) Aptr0 = Mptr( A, Aii, Ajj, Ald, size );

   if( notran )
   {
      top  = *PB_Ctop( &ctxt, COMBINE, ROW, TOP_GET );
      Cinb = DESCC[INB_]; Cnb = DESCC[NB_]; Ccsrc = DESCC[CSRC_];

      if( upper )
      {
         for( k = kstart; k != kend; k += kstep )
         {
            kbb = N - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA:JA+K-1 )
*/
            PB_CGatherV( TYPE, REUSE, &GatherDir, kbb, An, A, IA+k, JA, DESCA,
                         ROW, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+k+kbb-1, JA:JA+K-1 ) over A( IA:IA+k+kbb-1, JA:JA+K-1 )
*/
            PB_Cdescset( Ad0, ktmp, An, Aimb1, Ainb1, Amb, Anb, Arow, Acol,
                         ctxt, Ald );
            PB_CInV( TYPE, NOCONJG, ROW, ktmp, An, Ad0, kbb, Aptr, 0, 0, DBUFA,
                     ROW, &WA, WAd, &WAfr );
/*
*  WC := A( IA:IA+k+kbb-1, JA:JA+K-1 ) * A( IA+k:IA+k+kbb-1, JA:JA+K-1 )'
*/
            PB_COutV( TYPE, COLUMN, INIT, ktmp, An, Ad0, kbb, &WC, WCd, &WCfr,
                      &WCsum );
            Amp = PB_Cnumroc( ktmp, 0, Aimb1, Amb, myrow, Arow, nprow );
            if( ( Amp > 0 ) && ( Anq0 > 0 ) )
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( &tran ), &Amp, &kbb, &Anq0,
                     ALPHA, Aptr0, &Ald, WA, &WAd[LLD_], zero, WC, &WCd[LLD_] );
            if( WAfr ) free( WA   );
            if( Afr  ) free( Aptr );

            if( WCsum )
            {
               WCd[CSRC_] = PB_Cindxg2p( JC + ( fwd ? k : ktmp - 1 ), Cinb,
                                         Cnb, Ccsrc, Ccsrc, npcol );
               if( Amp > 0 )
                  gsum2d( ctxt, ROW, &top, Amp, kbb, WC, WCd[LLD_], myrow,
                          WCd[CSRC_] );
            }
/*
*  Zero lower triangle of WC( k:k+kbb-1, 0:kbb-1 )
*/
            if( conjg )
               PB_Cplapad( TYPE, LOWER, CONJG,   kbb,   kbb,   zero, zero, WC,
                           k,   0, WCd );
            else if( kbb > 1 )
               PB_Cplapad( TYPE, LOWER, NOCONJG, kbb-1, kbb-1, zero, zero, WC,
                           k+1, 0, WCd );
/*
*  Add WC to C( IC:IC+k+kbb-1, JC+k:JC+k+kbb-1 )
*/
            PB_CScatterV( TYPE, &ScatterDir, ktmp, kbb, WC, 0, 0, WCd, COLUMN,
                          one, C, IC, JC+k, DESCC, COLUMN );
            if( WCfr ) free( WC   );
         }
      }
      else
      {
         for( k = kstart; k != kend; k += kstep )
         {
            ktmp = N - k; kbb = MIN( ktmp, kb );
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA:JA+K-1 )
*/
            PB_CGatherV( TYPE, REUSE, &GatherDir, kbb, An, A, IA+k, JA, DESCA,
                         ROW, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA+k:IA+k+kbb-1, JA:JA+K-1 ) over A( IA+k:IA+N-1, JA:JA+K-1 )
*/
            Acurimb1 = PB_Cfirstnb( ktmp, IA+k, Aimb, Amb );
            Acurrow  = PB_Cindxg2p( k, Aimb1, Amb, Arow, Arow, nprow );
            PB_Cdescset( Ad0, ktmp, An, Acurimb1, Ainb1, Amb, Anb, Acurrow,
                         Acol, ctxt, Ald );
            PB_CInV( TYPE, NOCONJG, ROW, ktmp, An, Ad0, kbb, Aptr, 0, 0, DBUFA,
                     ROW, &WA, WAd, &WAfr );
/*
*  WC := A( IA+k:IA+N-1, JA:JA+K-1 ) * A( IA+k:IA+k+kbb-1, JA:JA+K-1 )'
*/
            PB_COutV( TYPE, COLUMN, INIT, ktmp, An, Ad0, kbb, &WC, WCd, &WCfr,
                      &WCsum );
            Amp = PB_Cnumroc( ktmp, k, Aimb1, Amb, myrow, Arow, nprow );
            if( ( Amp > 0 ) && ( Anq0 > 0 ) )
               gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( &tran ), &Amp, &kbb, &Anq0,
                     ALPHA, Mptr( Aptr0, Amp0-Amp, 0, Ald, size ), &Ald, WA,
                     &WAd[LLD_], zero, WC, &WCd[LLD_] );
            if( WAfr ) free( WA   );
            if( Afr  ) free( Aptr );

            if( WCsum )
            {
               WCd[CSRC_] = PB_Cindxg2p( JC + ( fwd ? k : k + kbb - 1 ), Cinb,
                                         Cnb, Ccsrc, Ccsrc, npcol );
               if( Amp > 0 )
                  gsum2d( ctxt, ROW, &top, Amp, kbb, WC, WCd[LLD_], myrow,
                          WCd[CSRC_] );
            }
/*
*  Zero upper triangle of WC
*/
            if( conjg )
               PB_Cplapad( TYPE, UPPER, CONJG,   kbb,   kbb,   zero, zero, WC,
                           0, 0, WCd );
            else if( kbb > 1 )
               PB_Cplapad( TYPE, UPPER, NOCONJG, kbb-1, kbb-1, zero, zero, WC,
                           0, 1, WCd );
/*
*  Add WC to C( IC+k:IC+N-1, JC+k:JC+k+kbb-1 )
*/
            PB_CScatterV( TYPE, &ScatterDir, ktmp, kbb, WC, 0, 0, WCd, COLUMN,
                          one, C, IC+k, JC+k, DESCC, COLUMN );
            if( WCfr ) free( WC   );
         }
      }
   }
   else
   {
      top  = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
      Cimb = DESCC[IMB_]; Cmb = DESCC[MB_]; Crsrc = DESCC[RSRC_];

      if( upper )
      {
         for( k = kstart; k != kend; k += kstep )
         {
            ktmp = N - k; kbb = MIN( ktmp, kb );
/*
*  Accumulate A( IA:IA+K-1, JA+k:JA+k+kbb-1 )
*/
            PB_CGatherV( TYPE, REUSE, &GatherDir, Am, kbb, A, IA, JA+k, DESCA,
                         COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA:IA+K-1, JA+k:JA+k+kbb-1 ) over A( IA:IA+K-1, JA+k:JA+N-1 )
*/
            Acurinb1 = PB_Cfirstnb( ktmp, JA+k, Ainb, Anb );
            Acurcol  = PB_Cindxg2p( k, Ainb1, Anb, Acol, Acol, npcol );
            PB_Cdescset( Ad0, Am, ktmp, Aimb1, Acurinb1, Amb, Anb, Arow,
                         Acurcol, ctxt, Ald );
            PB_CInV( TYPE, NOCONJG, COLUMN, Am, ktmp, Ad0, kbb, Aptr, 0, 0,
                     DBUFA, COLUMN, &WA, WAd, &WAfr );
/*
*  WC := A( IA:IA+K-1, JA+k:JA+k+kbb-1 )' * A( IA:IA+K-1,JA+k:JA+N-1 )
*/
            PB_COutV( TYPE, ROW, INIT, Am, ktmp, Ad0, kbb, &WC, WCd, &WCfr,
                      &WCsum );
            Anq = PB_Cnumroc( ktmp, k, Ainb1, Anb, mycol, Acol, npcol );
            if( ( Anq > 0 ) && ( Amp0 > 0 ) )
               gemm( C2F_CHAR( &tran ), C2F_CHAR( NOTRAN ), &kbb, &Anq, &Amp0,
                     ALPHA, WA, &WAd[LLD_], Mptr( Aptr0, 0, Anq0-Anq, Ald,
                     size ), &Ald, zero, WC, &WCd[LLD_] );
            if( WAfr ) free( WA   );
            if( Afr  ) free( Aptr );

            if( WCsum )
            {
               WCd[RSRC_] = PB_Cindxg2p( IC + ( fwd ? k : k + kbb - 1 ), Cimb,
                                         Cmb, Crsrc, Crsrc, nprow );
               if( Anq > 0 )
                  gsum2d( ctxt, COLUMN, &top, kbb, Anq, WC, WCd[LLD_],
                          WCd[RSRC_], mycol );
            }
/*
*  Zero lower triangle of WC
*/
            if( conjg )
               PB_Cplapad( TYPE, LOWER, CONJG,   kbb,   kbb,   zero, zero, WC,
                           0, 0, WCd );
            else if( kbb > 1 )
               PB_Cplapad( TYPE, LOWER, NOCONJG, kbb-1, kbb-1, zero, zero, WC,
                           1, 0, WCd );
/*
*  Add WC to C( IC+k:IC+k+kbb-1, JC+k:JC+N-1 )
*/
            PB_CScatterV( TYPE, &ScatterDir, kbb, ktmp, WC, 0, 0, WCd, ROW, one,
                          C, IC+k, JC+k, DESCC, ROW );
            if( WCfr ) free( WC   );
         }
      }
      else
      {
         for( k = kstart; k != kend; k += kstep )
         {
            kbb = N - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA:IA+K-1, JA+k:JA+k+kbb-1 )
*/
            PB_CGatherV( TYPE, REUSE, &GatherDir, Am, kbb, A, IA, JA+k, DESCA,
                         COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Replicate A( IA:IA+K-1, JA+k:JA+k+kbb-1 ) over A( IA:IA+K-1, JA:JA+k+kbb-1 )
*/
            PB_Cdescset( Ad0, Am, ktmp, Aimb1, Ainb1, Amb, Anb, Arow, Acol,
                         ctxt, Ald );
            PB_CInV( TYPE, NOCONJG, COLUMN, Am, ktmp, Ad0, kbb, Aptr, 0, 0,
                     DBUFA, COLUMN, &WA, WAd, &WAfr );
/*
*  WC := A( IA:IA+K-1, JA+k:JA+k+kbb-1 )' * A( IA:IA+K-1, JA:JA+k+kbb-1 )
*/
            PB_COutV( TYPE, ROW, INIT, Am, ktmp, Ad0, kbb, &WC, WCd, &WCfr,
                      &WCsum );
            Anq = PB_Cnumroc( ktmp, 0, Ainb1, Anb, mycol, Acol, npcol );
            if( ( Anq > 0 ) && ( Amp0 > 0 ) )
               gemm( C2F_CHAR( &tran ), C2F_CHAR( NOTRAN ), &kbb, &Anq, &Amp0,
                     ALPHA, WA, &WAd[LLD_], Aptr0, &Ald, zero, WC, &WCd[LLD_] );
            if( WAfr ) free( WA   );
            if( Afr  ) free( Aptr );

            if( WCsum )
            {
               WCd[RSRC_] = PB_Cindxg2p( IC + ( fwd ? k : ktmp - 1 ), Cimb,
                                         Cmb, Crsrc, Crsrc, nprow );
               if( Anq > 0 )
                  gsum2d( ctxt, COLUMN, &top, kbb, Anq, WC, WCd[LLD_],
                          WCd[RSRC_], mycol );
            }
/*
*  Zero upper triangle of WC( 0:kbb-1, k:k+kbb-1 )
*/
            if( conjg )
               PB_Cplapad( TYPE, UPPER, CONJG,   kbb,   kbb,   zero, zero, WC,
                           0, k,   WCd );
            else if( kbb > 1 )
               PB_Cplapad( TYPE, UPPER, NOCONJG, kbb-1, kbb-1, zero, zero, WC,
                           0, k+1, WCd );
/*
*  Add WC to C( IC+k:IC+k+kbb-1, JC:JC+k+kbb-1 )
*/
            PB_CScatterV( TYPE, &ScatterDir, kbb, ktmp, WC, 0, 0, WCd, ROW, one,
                          C, IC+k, JC, DESCC, ROW );
            if( WCfr ) free( WC   );
         }
      }
   }
/*
*  End of PB_CpsyrkAC
*/
}
