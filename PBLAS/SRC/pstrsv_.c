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
void pstrsv_( F_CHAR_T UPLO, F_CHAR_T TRANS, F_CHAR_T DIAG, Int * N,
              float * A, Int * IA, Int * JA, Int * DESCA,
              float * X, Int * IX, Int * JX, Int * DESCX,
              Int * INCX )
#else
void pstrsv_( UPLO, TRANS, DIAG, N, A, IA, JA, DESCA, X, IX, JX,
              DESCX, INCX )
/*
*  .. Scalar Arguments ..
*/
   F_CHAR_T       DIAG, TRANS, UPLO;
   Int            * IA, * INCX, * IX, * JA, * JX, * N;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA, * DESCX;
   float          * A, * X;
#endif
{
/*
*  Purpose
*  =======
*
*  PSTRSV  solves one of the systems of equations
*
*    sub( A )*sub( X ) = b,   or   sub( A )'*sub( X ) = b,
*
*  where
*
*     sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1), and,
*
*     sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
*                      X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X.
*
*  b and sub( X ) are  n element subvectors and  sub( A )  is an  n by n
*  unit, or non-unit, upper or lower triangular submatrix.
*
*  No test for  singularity  or  near-singularity  is included  in  this
*  routine. Such tests must be performed before calling this routine.
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
*          On entry,  UPLO  specifies whether the submatrix  sub( A ) is
*          an upper or lower triangular submatrix as follows:
*
*             UPLO = 'U' or 'u'   sub( A ) is an upper triangular
*                                 submatrix,
*
*             UPLO = 'L' or 'l'   sub( A ) is a  lower triangular
*                                 submatrix.
*
*  TRANS   (global input) CHARACTER*1
*          On entry,  TRANS  specifies the  operation to be performed as
*          follows:
*
*             TRANS = 'N' or 'n'   sub( A )  * sub( X ) = b.
*
*             TRANS = 'T' or 't'   sub( A )' * sub( X ) = b.
*
*             TRANS = 'C' or 'c'   sub( A )' * sub( X ) = b.
*
*  DIAG    (global input) CHARACTER*1
*          On entry,  DIAG  specifies  whether or not  sub( A )  is unit
*          triangular as follows:
*
*             DIAG = 'U' or 'u'  sub( A )  is  assumed to be unit trian-
*                                gular,
*
*             DIAG = 'N' or 'n'  sub( A ) is not assumed to be unit tri-
*                                angular.
*
*  N       (global input) INTEGER
*          On entry,  N specifies the order of the  submatrix  sub( A ).
*          N must be at least zero.
*
*  A       (local input) REAL array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix A.
*          Before entry with  UPLO = 'U' or 'u', this array contains the
*          local entries corresponding to  the entries of the upper tri-
*          angular submatrix  sub( A ), and the local entries correspon-
*          ding to the entries of the  strictly lower triangular part of
*          the submatrix  sub( A )  are not referenced.
*          Before entry with  UPLO = 'L' or 'l', this array contains the
*          local entries corresponding to  the entries of the lower tri-
*          angular submatrix  sub( A ), and the local entries correspon-
*          ding to the entries of the  strictly upper triangular part of
*          the submatrix  sub( A )  are not referenced.
*          Note  that  when DIAG = 'U' or 'u', the local entries corres-
*          ponding to  the diagonal elements  of the submatrix  sub( A )
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
*  X       (local input/local output) REAL array
*          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
*          is   at  least  MAX( 1, Lr( 1, IX ) )  when  INCX = M_X   and
*          MAX( 1, Lr( 1, IX+N-1 ) )  otherwise,  and,  Kx  is  at least
*          Lc( 1, JX+N-1 )  when  INCX = M_X  and Lc( 1, JX ) otherwise.
*          Before entry, this array  contains  the local  entries of the
*          matrix X. On entry, sub( X ) is the n element right-hand side
*          b. On exit, sub( X ) is overwritten with the solution subvec-
*          tor.
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
   char           DiagA, TranOp, UploA, Xroc, btop, ctop, * negone, * one,
                  * zero;
   Int            Acol, Ai, Aii, Aimb1, Ainb1, Aj, Ajj, Akp, Akq, Ald, Amb, Anb,
                  Anp, Anp0, Anq, Anq0, Arow, Asrc, XACapbX, XACfr, XACld,
                  XACsum, XARapbX, XARfr, XARld, XARsum, Xi, Xj, ctxt, info,
                  ione=1, k, kb, kbnext, kbprev, ktmp, mycol, myrow, nb, notran,
                  nounit, npcol, nprow, size, upper;
   PBTYP_T        * type;
/*
*  .. Local Arrays ..
*/
   Int            Ad[DLEN_], Ad0[DLEN_], XACd[DLEN_], XARd[DLEN_], Xd[DLEN_];
   char           * Aptr = NULL, * XAC = NULL, * XAR = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   upper  = ( ( UploA  = Mupcase( F2C_CHAR( UPLO  )[0] ) ) ==  CUPPER );
   notran = ( ( TranOp = Mupcase( F2C_CHAR( TRANS )[0] ) ) == CNOTRAN );
   nounit = ( ( DiagA  = Mupcase( F2C_CHAR( DIAG  )[0] ) ) == CNOUNIT );
   PB_CargFtoC( *IA, *JA, DESCA, &Ai, &Aj, Ad );
   PB_CargFtoC( *IX, *JX, DESCX, &Xi, &Xj, Xd );
#ifndef NO_ARGCHK
/*
*  Test the input parameters
*/
   Cblacs_gridinfo( ( ctxt = Ad[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
   if( !( info = ( ( nprow == -1 ) ? -( 801 + CTXT_ ) : 0 ) ) )
   {
      if( ( !upper ) && ( UploA != CLOWER ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PSTRSV", "Illegal UPLO = %c\n", UploA );
         info = -1;
      }
      else if( ( !notran ) && ( TranOp != CTRAN ) && ( TranOp != CCOTRAN ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PSTRSV", "Illegal TRANS = %c\n", TranOp );
         info = -2;
      }
      else if( ( !nounit ) && ( DiagA != CUNIT ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PSTRSV", "Illegal DIAG = %c\n", DiagA );
         info = -3;
      }
      PB_Cchkmat( ctxt, "PSTRSV", "A", *N, 4, *N, 4, Ai, Aj, Ad,  8, &info );
      PB_Cchkvec( ctxt, "PSTRSV", "X", *N, 4, Xi, Xj, Xd, *INCX, 12, &info );
   }
   if( info ) { PB_Cabort( ctxt, "PSTRSV", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( *N == 0 ) return;
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
   size = type->size; one = type->one; zero = type->zero; negone = type->negone;
/*
*  Compute descriptor Ad0 for sub( A )
*/
   PB_Cdescribe( *N, *N, Ai, Aj, Ad, nprow, npcol, myrow, mycol, &Aii, &Ajj,
                 &Ald, &Aimb1, &Ainb1, &Amb, &Anb, &Arow, &Acol, Ad0 );
/*
*  Computational partitioning size is computed as the product of the logical
*  value returned by pilaenv_ and 2 * lcm( nprow, npcol )
*/
   nb = 2 * pilaenv_( &ctxt, C2F_CHAR( &type->type ) ) *
        PB_Clcm( ( Arow >= 0 ? nprow : 1 ), ( Acol >= 0 ? npcol : 1 ) );

   Xroc = ( *INCX == Xd[M_] ? CROW : CCOLUMN );

   if( notran )
   {
      if( upper )
      {
/*
*  Save current and enforce ring BLACS topologies
*/
         btop = *PB_Ctop( &ctxt, BCAST,   COLUMN, TOP_GET   );
         ctop = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_GET   );
         (void)  PB_Ctop( &ctxt, BCAST,   COLUMN, TOP_DRING );
         (void)  PB_Ctop( &ctxt, COMBINE, ROW,    TOP_DRING );
/*
*  Remove next line when BLACS combine operations support ring topologies.
*/
         (void)  PB_Ctop( &ctxt, COMBINE, ROW,    TOP_DEFAULT );
/*
*  Reuse sub( X ) and/or create vector XAC in process column owning the last
*  column of sub( A )
*/
         PB_CInOutV2( type, NOCONJG, COLUMN, *N, *N, *N-1, Ad0, 1,
                      ((char *) X), Xi, Xj, Xd, &Xroc, &XAC, XACd,
                      &XACfr, &XACsum, &XACapbX );
/*
*  Create vector XAR in process rows spanned by sub( A )
*/
         PB_COutV( type, ROW,    INIT, *N, *N, Ad0, 1, &XAR, XARd, &XARfr,
                   &XARsum );
/*
*  Retrieve local quantities related to Ad0 -> sub( A )
*/
         Aimb1 = Ad0[IMB_ ]; Ainb1 = Ad0[INB_ ];
         Amb   = Ad0[MB_  ]; Anb   = Ad0[NB_  ];
         Arow  = Ad0[RSRC_]; Acol  = Ad0[CSRC_]; Ald = Ad0[LLD_];
         Anp   = PB_Cnumroc( *N, 0, Aimb1, Amb, myrow, Arow, nprow );
         Anq   = PB_Cnumroc( *N, 0, Ainb1, Anb, mycol, Acol, npcol );
         if( ( Anp > 0 ) && ( Anq > 0 ) )
            Aptr = Mptr( ((char *) A), Aii, Ajj, Ald, size );

         XACld = XACd[LLD_]; XARld = XARd[LLD_];

         for( k = ( ( *N - 1 ) / nb ) * nb; k >= 0; k -= nb )
         {
            ktmp = *N - k; kb = MIN( ktmp, nb );
/*
*  Solve logical diagonal block, XAC contains the solution scattered in multiple
*  process columns and XAR contains the solution replicated in the process rows.
*/
            Akp = PB_Cnumroc( k, 0, Aimb1, Amb, myrow, Arow, nprow );
            Akq = PB_Cnumroc( k, 0, Ainb1, Anb, mycol, Acol, npcol );
            PB_Cptrsv( type, XARsum, &UploA, &TranOp, &DiagA, kb, Aptr, k, k,
                       Ad0, Mptr( XAC, Akp, 0, XACld, size ), 1, Mptr( XAR, 0,
                       Akq, XARld, size ), XARld );
/*
*  Update: only the part of sub( X ) to be solved at the next step is locally
*  updated and combined, the remaining part of the vector to be solved later is
*  only locally updated.
*/
            if( Akp > 0 )
            {
               Anq0 = PB_Cnumroc( kb, k, Ainb1, Anb, mycol, Acol, npcol );
               if( XACsum )
               {
                  kbprev = MIN( k, nb );
                  ktmp   = PB_Cnumroc( kbprev, k-kbprev, Aimb1, Amb, myrow,
                                       Arow, nprow );
                  Akp   -= ktmp;

                  if( ktmp > 0 )
                  {
                     if( Anq0 > 0 )
                        sgemv_( TRANS, &ktmp, &Anq0, negone,
                           Mptr( Aptr, Akp, Akq,   Ald, size ), &Ald,
                           Mptr( XAR,    0, Akq, XARld, size ), &XARld, one,
                           Mptr( XAC,  Akp,   0, XACld, size ), &ione );
                     Asrc = PB_Cindxg2p( k-1, Ainb1, Anb, Acol, Acol, npcol );
                     Csgsum2d( ctxt, ROW, &ctop, ktmp, 1, Mptr( XAC, Akp,
                               0, XACld, size ), XACld, myrow, Asrc );
                     if( mycol != Asrc )
                        sset_( &ktmp, zero, Mptr( XAC, Akp, 0, XACld,
                               size ), &ione );
                  }
                  if( Akp > 0 && Anq0 > 0 )
                     sgemv_( TRANS, &Akp, &Anq0, negone,
                             Mptr( Aptr, 0, Akq,   Ald, size ),   &Ald,
                             Mptr( XAR,  0, Akq, XARld, size ), &XARld, one,
                             XAC, &ione );
               }
               else
               {
                  if( Anq0 > 0 )
                     sgemv_( TRANS, &Akp, &Anq0, negone,
                             Mptr( Aptr, 0, Akq,   Ald, size ),   &Ald,
                             Mptr( XAR,  0, Akq, XARld, size ), &XARld, one,
                             XAC, &ione );
               }
            }
         }
/*
*  Combine the scattered resulting vector XAC
*/
         if( XACsum && ( Anp > 0 ) )
         {
            Csgsum2d( ctxt, ROW, &ctop, Anp, 1, XAC, XACld, myrow,
                      XACd[CSRC_] );
         }
/*
*  sub( X ) := XAC (if necessary)
*/
         if( XACapbX )
         {
            PB_Cpaxpby( type, NOCONJG, *N, 1, one, XAC, 0, 0, XACd, COLUMN,
                        zero, ((char *) X), Xi, Xj, Xd, &Xroc );
         }
/*
*  Restore BLACS topologies
*/
         (void) PB_Ctop( &ctxt, BCAST,   COLUMN, &btop );
         (void) PB_Ctop( &ctxt, COMBINE, ROW,    &ctop );
      }
      else
      {
/*
*  Save current and enforce ring BLACS topologies
*/
         btop = *PB_Ctop( &ctxt, BCAST,   COLUMN, TOP_GET   );
         ctop = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_GET   );
         (void)  PB_Ctop( &ctxt, BCAST,   COLUMN, TOP_IRING );
         (void)  PB_Ctop( &ctxt, COMBINE, ROW,    TOP_IRING );
/*
*  Remove next line when BLACS combine operations support ring topologies.
*/
         (void)  PB_Ctop( &ctxt, COMBINE, ROW,    TOP_DEFAULT );
/*
*  Reuse sub( X ) and/or create vector XAC in process column owning the first
*  column of sub( A )
*/
         PB_CInOutV2( type, NOCONJG, COLUMN, *N, *N, 0, Ad0, 1,
                      ((char *) X), Xi, Xj, Xd, &Xroc, &XAC, XACd,
                      &XACfr, &XACsum, &XACapbX );
/*
*  Create vector XAR in process rows spanned by sub( A )
*/
         PB_COutV( type, ROW, INIT, *N, *N, Ad0, 1, &XAR, XARd, &XARfr,
                   &XARsum );
/*
*  Retrieve local quantities related to Ad0 -> sub( A )
*/
         Aimb1 = Ad0[IMB_ ]; Ainb1 = Ad0[INB_ ];
         Amb   = Ad0[MB_  ]; Anb   = Ad0[NB_  ];
         Arow  = Ad0[RSRC_]; Acol  = Ad0[CSRC_]; Ald = Ad0[LLD_];
         Anp   = PB_Cnumroc( *N, 0, Aimb1, Amb, myrow, Arow, nprow );
         Anq   = PB_Cnumroc( *N, 0, Ainb1, Anb, mycol, Acol, npcol );
         if( ( Anp > 0 ) && ( Anq > 0 ) )
            Aptr = Mptr( ((char *) A), Aii, Ajj, Ald, size );

         XACld = XACd[LLD_]; XARld = XARd[LLD_];

         for( k = 0; k < *N; k += nb )
         {
            ktmp = *N - k; kb = MIN( ktmp, nb );
/*
*  Solve logical diagonal block, XAC contains the solution scattered in multiple
*  process columns and XAR contains the solution replicated in the process rows.
*/
            Akp = PB_Cnumroc( k, 0, Aimb1, Amb, myrow, Arow, nprow );
            Akq = PB_Cnumroc( k, 0, Ainb1, Anb, mycol, Acol, npcol );
            PB_Cptrsv( type, XARsum, &UploA, &TranOp, &DiagA, kb, Aptr, k, k,
                       Ad0, Mptr( XAC, Akp, 0, XACld, size ), 1, Mptr( XAR, 0,
                       Akq, XARld, size ), XARld );
/*
*  Update: only the part of sub( X ) to be solved at the next step is locally
*  updated and combined, the remaining part of the vector to be solved later is
*  only locally updated.
*/
            Akp = PB_Cnumroc( k+kb, 0, Aimb1, Amb, myrow, Arow, nprow );
            if( ( Anp0 = Anp - Akp ) > 0 )
            {
               Anq0 = PB_Cnumroc( kb, k, Ainb1, Anb, mycol, Acol, npcol );
               if( XACsum )
               {
                  kbnext = ktmp - kb; kbnext = MIN( kbnext, nb );
                  ktmp   = PB_Cnumroc( kbnext, k+kb, Aimb1, Amb, myrow, Arow,
                                       nprow );
                  Anp0  -= ktmp;

                  if( ktmp > 0 )
                  {
                     if( Anq0 > 0 )
                        sgemv_( TRANS, &ktmp, &Anq0, negone,
                          Mptr( Aptr, Akp, Akq,   Ald, size ), &Ald,
                          Mptr( XAR,    0, Akq, XARld, size ), &XARld, one,
                          Mptr( XAC,  Akp,   0, XACld, size ), &ione );
                     Asrc = PB_Cindxg2p( k+kb, Ainb1, Anb, Acol, Acol, npcol );
                     Csgsum2d( ctxt, ROW, &ctop, ktmp, 1, Mptr( XAC, Akp,
                                 0, XACld, size ), XACld, myrow, Asrc );
                     if( mycol != Asrc )
                        sset_( &ktmp, zero, Mptr( XAC, Akp, 0, XACld,
                               size ), &ione );
                  }
                  if( Anp0 > 0 && Anq0 > 0 )
                     sgemv_( TRANS, &Anp0, &Anq0, negone,
                      Mptr( Aptr, Akp+ktmp, Akq,   Ald, size ), &Ald,
                      Mptr( XAR,         0, Akq, XARld, size ), &XARld,
                      one,
                      Mptr( XAC,  Akp+ktmp,   0, XACld, size ), &ione );
               }
               else
               {
                  if( Anq0 > 0 )
                     sgemv_( TRANS, &Anp0, &Anq0, negone,
                          Mptr( Aptr, Akp, Akq,   Ald, size ), &Ald,
                          Mptr( XAR,    0, Akq, XARld, size ), &XARld,
                          one,
                          Mptr( XAC,  Akp,   0, XACld, size ), &ione );
               }
            }
         }
/*
*  Combine the scattered resulting vector XAC
*/
         if( XACsum && ( Anp > 0 ) )
         {
            Csgsum2d( ctxt, ROW, &ctop, Anp, 1, XAC, XACld, myrow,
                      XACd[CSRC_] );
         }
/*
*  sub( X ) := XAC (if necessary)
*/
         if( XACapbX )
         {
            PB_Cpaxpby( type, NOCONJG, *N, 1, one, XAC, 0, 0, XACd,
                        COLUMN, zero, ((char *) X), Xi, Xj, Xd, &Xroc );
         }
/*
*  Restore BLACS topologies
*/
         (void) PB_Ctop( &ctxt, BCAST,   COLUMN, &btop );
         (void) PB_Ctop( &ctxt, COMBINE, ROW,    &ctop );
      }
   }
   else
   {
      if( upper )
      {
/*
*  Save current and enforce ring BLACS topologies
*/
         btop = *PB_Ctop( &ctxt, BCAST,   ROW,    TOP_GET   );
         ctop = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET   );
         (void)  PB_Ctop( &ctxt, BCAST,   ROW,    TOP_IRING );
         (void)  PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_IRING );
/*
*  Remove next line when BLACS combine operations support ring topologies.
*/
         (void)  PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_DEFAULT );
/*
*  Reuse sub( X ) and/or create vector XAR in process row owning the first row
*  of sub( A )
*/
         PB_CInOutV2( type, NOCONJG, ROW, *N, *N, 0, Ad0, 1,
                      ((char *) X), Xi, Xj, Xd, &Xroc, &XAR, XARd,
                      &XARfr, &XARsum, &XARapbX );
/*
*  Create vector XAC in process columns spanned by sub( A )
*/
         PB_COutV( type, COLUMN, INIT, *N, *N, Ad0, 1, &XAC, XACd, &XACfr,
                   &XACsum );
/*
*  Retrieve local quantities related to Ad0 -> sub( A )
*/
         Aimb1 = Ad0[IMB_ ]; Ainb1 = Ad0[INB_ ];
         Amb   = Ad0[MB_  ]; Anb   = Ad0[NB_  ];
         Arow  = Ad0[RSRC_]; Acol  = Ad0[CSRC_]; Ald = Ad0[LLD_];
         Anp   = PB_Cnumroc( *N, 0, Aimb1, Amb, myrow, Arow, nprow );
         Anq   = PB_Cnumroc( *N, 0, Ainb1, Anb, mycol, Acol, npcol );
         if( ( Anp > 0 ) && ( Anq > 0 ) )
            Aptr = Mptr( ((char *) A), Aii, Ajj, Ald, size );

         XACld = XACd[LLD_]; XARld = XARd[LLD_];

         for( k = 0; k < *N; k += nb )
         {
            ktmp = *N - k; kb = MIN( ktmp, nb );
/*
*  Solve logical diagonal block, XAR contains the solution scattered in multiple
*  process rows and XAC contains the solution replicated in the process columns.
*/
            Akp = PB_Cnumroc( k, 0, Aimb1, Amb, myrow, Arow, nprow );
            Akq = PB_Cnumroc( k, 0, Ainb1, Anb, mycol, Acol, npcol );
            PB_Cptrsv( type, XACsum, &UploA, &TranOp, &DiagA, kb, Aptr, k, k,
                       Ad0, Mptr( XAC, Akp, 0, XACld, size ), 1, Mptr( XAR, 0,
                       Akq, XARld, size ), XARld );
/*
*  Update: only the part of sub( X ) to be solved at the next step is locally
*  updated and combined, the remaining part of the vector to be solved later is
*  only locally updated.
*/
            Akq = PB_Cnumroc( k+kb, 0, Ainb1, Anb, mycol, Acol, npcol );
            if( ( Anq0 = Anq - Akq ) > 0 )
            {
               Anp0 = PB_Cnumroc( kb, k, Aimb1, Amb, myrow, Arow, nprow );
               if( XARsum )
               {
                  kbnext = ktmp - kb; kbnext = MIN( kbnext, nb );
                  ktmp   = PB_Cnumroc( kbnext, k+kb, Ainb1, Anb, mycol, Acol,
                                       npcol );
                  Anq0  -= ktmp;

                  if( ktmp > 0 )
                  {
                     if( Anp0 > 0 )
                        sgemv_( TRANS, &Anp0, &ktmp, negone,
                         Mptr( Aptr, Akp, Akq,   Ald, size ), &Ald,
                         Mptr( XAC,  Akp,   0, XACld, size ), &ione, one,
                         Mptr( XAR,    0, Akq, XARld, size ), &XARld );
                     Asrc = PB_Cindxg2p( k+kb, Aimb1, Amb, Arow, Arow, nprow );
                     Csgsum2d( ctxt, COLUMN, &ctop, 1, ktmp, Mptr( XAR, 0,
                               Akq, XARld, size ), XARld, Asrc, mycol );
                     if( myrow != Asrc )
                        sset_( &ktmp, zero, Mptr( XAR, 0, Akq, XARld,
                               size ), &XARld );
                  }
                  if( Anp0 > 0 && Anq0 > 0 )
                     sgemv_( TRANS, &Anp0, &Anq0, negone,
                     Mptr( Aptr, Akp, Akq+ktmp,   Ald, size ),   &Ald,
                     Mptr( XAC,  Akp,        0, XACld, size ),  &ione, one,
                     Mptr( XAR,    0, Akq+ktmp, XARld, size ), &XARld );
               }
               else
               {
                  if( Anp0 > 0 )
                     sgemv_( TRANS, &Anp0, &Anq0, negone,
                        Mptr( Aptr, Akp, Akq,   Ald, size ),   &Ald,
                        Mptr( XAC,  Akp,   0, XACld, size ),  &ione, one,
                        Mptr( XAR,    0, Akq, XARld, size ), &XARld );
               }
            }
         }
/*
*  Combine the scattered resulting vector XAR
*/
         if( XARsum && ( Anq > 0 ) )
         {
            Csgsum2d( ctxt, COLUMN, &ctop, 1, Anq, XAR, XARld, XARd[RSRC_],
                      mycol );
         }
/*
*  sub( X ) := XAR (if necessary)
*/
         if( XARapbX )
         {
            PB_Cpaxpby( type, NOCONJG, 1, *N, one, XAR, 0, 0, XARd, ROW, zero,
                        ((char *) X), Xi, Xj, Xd, &Xroc );
         }
/*
*  Restore BLACS topologies
*/
         (void) PB_Ctop( &ctxt, BCAST,   ROW,    &btop );
         (void) PB_Ctop( &ctxt, COMBINE, COLUMN, &ctop );
      }
      else
      {
/*
*  Save current and enforce ring BLACS topologies
*/
         btop = *PB_Ctop( &ctxt, BCAST,   ROW,    TOP_GET   );
         ctop = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET   );
         (void)  PB_Ctop( &ctxt, BCAST,   ROW,    TOP_DRING );
         (void)  PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_DRING );
/*
*  Remove next line when BLACS combine operations support ring topologies.
*/
         (void)  PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_DEFAULT );
/*
*  Reuse sub( X ) and/or create vector XAC in process row owning the last row
*  of sub( A )
*/
         PB_CInOutV2( type, NOCONJG, ROW, *N, *N, *N-1, Ad0, 1,
                      ((char *) X), Xi, Xj, Xd, &Xroc, &XAR, XARd,
                      &XARfr, &XARsum, &XARapbX );
/*
*  Create vector XAC in process columns spanned by sub( A )
*/
         PB_COutV( type, COLUMN, INIT, *N, *N, Ad0, 1, &XAC, XACd, &XACfr,
                   &XACsum );
/*
*  Retrieve local quantities related to Ad0 -> sub( A )
*/
         Aimb1 = Ad0[IMB_ ]; Ainb1 = Ad0[INB_ ];
         Amb   = Ad0[MB_  ]; Anb   = Ad0[NB_  ];
         Arow  = Ad0[RSRC_]; Acol  = Ad0[CSRC_]; Ald = Ad0[LLD_];
         Anp   = PB_Cnumroc( *N, 0, Aimb1, Amb, myrow, Arow, nprow );
         Anq   = PB_Cnumroc( *N, 0, Ainb1, Anb, mycol, Acol, npcol );
         if( ( Anp > 0 ) && ( Anq > 0 ) )
            Aptr = Mptr( ((char *) A), Aii, Ajj, Ald, size );

         XACld = XACd[LLD_]; XARld = XARd[LLD_];

         for( k = ( ( *N - 1 ) / nb ) * nb; k >= 0; k -= nb )
         {
            ktmp = *N - k; kb = MIN( ktmp, nb );
/*
*  Solve logical diagonal block, XAR contains the solution scattered in multiple
*  process rows and XAC contains the solution replicated in the process columns.
*/
            Akp = PB_Cnumroc( k,  0, Aimb1, Amb, myrow, Arow, nprow );
            Akq = PB_Cnumroc( k,  0, Ainb1, Anb, mycol, Acol, npcol );
            PB_Cptrsv( type, XACsum, &UploA, &TranOp, &DiagA, kb, Aptr, k, k,
                       Ad0, Mptr( XAC, Akp, 0, XACld, size ), 1, Mptr( XAR, 0,
                       Akq, XARld, size ), XARld );
/*
*  Update: only the part of sub( X ) to be solved at the next step is locally
*  updated and combined, the remaining part of the vector to be solved later
*  is only locally updated.
*/
            if( Akq > 0 )
            {
               Anp0 = PB_Cnumroc( kb, k, Aimb1, Amb, myrow, Arow, nprow );
               if( XARsum )
               {
                  kbprev = MIN( k, nb );
                  ktmp   = PB_Cnumroc( kbprev, k-kbprev, Ainb1, Anb, mycol,
                                       Acol, npcol );
                  Akq   -= ktmp;

                  if( ktmp > 0 )
                  {
                     if( Anp0 > 0 )
                        sgemv_( TRANS, &Anp0, &ktmp, negone,
                                Mptr( Aptr, Akp, Akq,   Ald, size ), &Ald,
                                Mptr( XAC,  Akp,   0, XACld, size ), &ione, one,
                                Mptr( XAR,    0, Akq, XARld, size ), &XARld );
                     Asrc = PB_Cindxg2p( k-1, Aimb1, Amb, Arow, Arow, nprow );
                     Csgsum2d( ctxt, COLUMN, &ctop, 1, ktmp, Mptr( XAR, 0,
                               Akq, XARld, size ), XARld, Asrc, mycol );
                     if( myrow != Asrc )
                        sset_( &ktmp, zero, Mptr( XAR, 0, Akq, XARld,
                               size ), &XARld );
                  }
                  if( Anp0 > 0 && Akq > 0 )
                     sgemv_( TRANS, &Anp0, &Akq, negone,
                             Mptr( Aptr, Akp, 0,   Ald, size ), &Ald,
                             Mptr( XAC,  Akp, 0, XACld, size ), &ione,
                             one, XAR, &XARld );
               }
               else
               {
                  if( Anp0 > 0 )
                     sgemv_( TRANS, &Anp0, &Akq, negone,
                             Mptr( Aptr, Akp, 0,   Ald, size ), &Ald,
                             Mptr( XAC,  Akp, 0, XACld, size ), &ione,
                             one, XAR, &XARld );
               }
            }
         }
/*
*  Combine the scattered resulting vector XAR
*/
         if( XARsum && ( Anq > 0 ) )
         {
            Csgsum2d( ctxt, COLUMN, &ctop, 1, Anq, XAR, XARld, XARd[RSRC_],
                      mycol );
         }
/*
*  sub( X ) := XAR (if necessary)
*/
         if( XARapbX )
         {
            PB_Cpaxpby( type, NOCONJG, 1, *N, one, XAR, 0, 0, XARd, ROW, zero,
                        ((char *) X), Xi, Xj, Xd, &Xroc );
         }
/*
*  Restore BLACS topologies
*/
         (void) PB_Ctop( &ctxt, BCAST,   ROW,    &btop );
         (void) PB_Ctop( &ctxt, COMBINE, COLUMN, &ctop );
      }
   }
   if( XACfr ) free( XAC );
   if( XARfr ) free( XAR );
/*
*  End of PSTRSV
*/
}
