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
void PB_CptrsmB( PBTYP_T * TYPE, char * DIRECB, char * SIDE,
                 char * UPLO, char * TRANSA, char * DIAG, int M, int N,
                 char * ALPHA, char * A, int IA, int JA, int * DESCA,
                 char * B, int IB, int JB, int * DESCB )
#else
void PB_CptrsmB( TYPE, DIRECB, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A,
                 IA, JA, DESCA, B, IB, JB, DESCB )
/*
*  .. Scalar Arguments ..
*/
   char           * DIAG, * DIRECB, * SIDE, * TRANSA, * UPLO;
   int            IA, IB, JA, JB, M, N;
   char           * ALPHA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCB;
   char           * A, * B;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CptrsmB  solves one of the matrix equations
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
*  This is the inner-product algorithm  using  the  logical  LCM  hybrid
*  and static blocking techniques. The submatrix operand  sub( A ) stays
*  in place.
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
*  DIRECB  (global input) pointer to CHAR
*          On entry, DIRECB  specifies  the direction in which the  rows
*          or columns of sub( B ) should be looped  over as follows:
*             DIRECB = 'F' or 'f'   forward  or increasing,
*             DIRECB = 'B' or 'b'   backward or decreasing.
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
   char           Broc, TranOp, conjg, * negone, * one, * talpha, * talph0, top,
                  * zero;
   int            Acol, Aii, Aimb1, Ainb1, Ajj, Akp, Akq, Alcmb, Ald, Amb, An,
                  Anb, Anp, Anp0, Anq, Anq0, Arow, Asrc, Astart, BcurrocR, Bfwd,
                  BiiD, BiiR, Binb1D, Binb1R, BisR, Bld, BmyprocD, BmyprocR,
                  BnD, BnR, BnbD, BnbR, BnpR, BnprocsD, BnprocsR, BrocD, BrocR,
                  BsrcR, LNorRT, WBCfr, WBCld, WBCapbX, WBCsum, WBRfr, WBRld,
                  WBRapbX, WBRsum, ctxt, izero=0, k, kb, kbnext, kbprev, ktmp,
                  lside, mycol, myrow, n, nb, nbb, notran, npcol, nprow, p=0,
                  size, tmp, upper;
   TZPAD_T        pad;
   GEMM_T         gemm;
   GSUM2D_T       gsum2d;
/*
*  .. Local Arrays ..
*/
   int            Ad0[DLEN_], DBUFB[DLEN_], WBCd[DLEN_], WBRd[DLEN_];
   char           * Aptr = NULL, * Bptr = NULL, * WBC = NULL, * WBR = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   Bfwd   = ( Mupcase( DIRECB[0] ) == CFORWARD );
   lside  = ( Mupcase( SIDE  [0] ) == CLEFT    );
   upper  = ( Mupcase( UPLO  [0] ) == CUPPER   );
   notran = ( ( TranOp = Mupcase( TRANSA[0] ) ) == CNOTRAN  );
   LNorRT = ( lside && notran ) || ( !( lside ) && !( notran ) );

   size   = TYPE->size;
   one    = TYPE->one;    zero = TYPE->zero;  negone = TYPE->negone;
   pad    = TYPE->Ftzpad; gemm = TYPE->Fgemm; gsum2d = TYPE->Cgsum2d;
   nb     = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );
/*
*  Compute local information for sub( A ) and sub( B )
*/
   if( lside )
   {
      BnD = An = M;            BnR      = N;           Broc  = CCOLUMN;
      BmyprocD = myrow;        BnprocsD = nprow;
      BmyprocR = mycol;        BnprocsR = npcol;
      BnbD     = DESCB[MB_  ]; BnbR     = DESCB[NB_ ];
      BsrcR    = DESCB[CSRC_]; Bld      = DESCB[LLD_];
      PB_Cinfog2l( IB, JB, DESCB, BnprocsD, BnprocsR, BmyprocD, BmyprocR,
                   &BiiD, &BiiR, &BrocD, &BrocR );
      Binb1D   = PB_Cfirstnb( BnD, IB, DESCB[IMB_], BnbD );
      Binb1R   = PB_Cfirstnb( BnR, JB, DESCB[INB_], BnbR );
   }
   else
   {
      BnD = An = N;            BnR      = M;           Broc  = CROW;
      BmyprocD = mycol;        BnprocsD = npcol;
      BmyprocR = myrow;        BnprocsR = nprow;
      BnbR     = DESCB[MB_  ]; BnbD     = DESCB[NB_ ];
      BsrcR    = DESCB[RSRC_]; Bld      = DESCB[LLD_];
      PB_Cinfog2l( IB, JB, DESCB, BnprocsR, BnprocsD, BmyprocR, BmyprocD,
                   &BiiR, &BiiD, &BrocR, &BrocD );
      Binb1D   = PB_Cfirstnb( BnD, JB, DESCB[INB_], BnbD );
      Binb1R   = PB_Cfirstnb( BnR, IB, DESCB[IMB_], BnbR );
   }
/*
*  Compute descriptor Ad0 for sub( A )
*/
   PB_Cdescribe( An, An, IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aii, &Ajj,
                 &Ald, &Aimb1, &Ainb1, &Amb, &Anb, &Arow, &Acol, Ad0 );
/*
*  Compute conjugate of alpha for the conjugate transpose cases
*/
   if( TranOp == CCOTRAN )
   {
      conjg = CCONJG; talpha = PB_Cmalloc( size );
      PB_Cconjg( TYPE, ALPHA, talpha );
   }
   else { conjg = CNOCONJG; talpha = ALPHA; }
/*
*  Retrieve BLACS combine topology, select backward ot forward substitution.
*/
   if( LNorRT )
   {
      top    = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_GET );
      Astart = ( upper ? An - 1 : 0 );
   }
   else
   {
      top    = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
      Astart = ( upper ? 0 : An - 1 );
   }
/*
*  Computational partitioning size is computed as the product of the logical
*  value returned by pilaenv_ and 2 * lcm( nprow, npcol ).
*/
   Alcmb  = 2 * nb * PB_Clcm( ( Arow >= 0 ? nprow : 1 ),
                              ( Acol >= 0 ? npcol : 1 ) );
/*
*  When sub( B ) is not replicated and backward pass on sub( B ), find the
*  virtual process p owning the last row or column of sub( B ).
*/
   if( !( BisR = ( ( BsrcR < 0 ) || ( BnprocsR == 1 ) ) ) && !( Bfwd ) )
   {
      tmp = PB_Cindxg2p( BnR - 1, Binb1R, BnbR, BrocR, BrocR, BnprocsR );
      p   = MModSub( tmp, BrocR, BnprocsR );
   }
/*
*  Loop over the processes rows or columns owning the BnR rows or columns of
*  sub( B ) to be processed.
*/
   n = BnR;

   while( n > 0 )
   {
/*
*  Find out who is the active process row or column as well as the number of
*  rows or columns of sub( B ) it owns.
*/
      BcurrocR = ( BisR ? -1 : MModAdd( BrocR, p, BnprocsR ) );
      BnpR     = PB_Cnumroc( BnR, 0, Binb1R, BnbR, BcurrocR, BrocR, BnprocsR );

      n       -= BnpR;
/*
*  Re-adjust the number of rows or columns to be handled at each step, in order
*  to average the message sizes and the computational granularity.
*/
      if( BnpR ) nbb = BnpR / ( ( BnpR - 1 ) / nb + 1 );

      while( BnpR )
      {
         nbb = MIN( nbb, BnpR );
/*
*  Describe the local contiguous panel of sub( B )
*/
         if( lside )
         {
            PB_Cdescset( DBUFB, BnD, nbb, Binb1D, nbb, BnbD, BnbR, BrocD,
                         BcurrocR, ctxt, Bld );
            if( BisR || ( BmyprocR == BcurrocR ) )
               Bptr = Mptr( B, BiiD, BiiR, Bld, size );
         }
         else
         {
            PB_Cdescset( DBUFB, nbb, BnD, nbb, Binb1D, BnbR, BnbD, BcurrocR,
                         BrocD, ctxt, Bld );
            if( BisR || ( BmyprocR == BcurrocR ) )
               Bptr = Mptr( B, BiiR, BiiD, Bld, size );
         }

         talph0 = talpha;

         if( LNorRT )
         {
/*
*  Reuse sub( B ) and/or create vector WBC in process column owning the first
*  or last column of sub( A )
*/
            PB_CInOutV2( TYPE, &conjg, COLUMN, An, An, Astart, Ad0, nbb, Bptr,
                         0, 0, DBUFB, &Broc, &WBC, WBCd, &WBCfr, &WBCsum,
                         &WBCapbX );
/*
*  Create WBR in process rows spanned by sub( A )
*/
            PB_COutV( TYPE, ROW, INIT, An, An, Ad0, nbb, &WBR, WBRd, &WBRfr,
                      &WBRsum );
/*
*  Retrieve local quantities related to sub( A ) -> Ad0
*/
            Aimb1 = Ad0[IMB_ ]; Ainb1 = Ad0[INB_ ];
            Amb   = Ad0[MB_  ]; Anb   = Ad0[NB_  ];
            Arow  = Ad0[RSRC_]; Acol  = Ad0[CSRC_]; Ald = Ad0[LLD_];

            Anp   = PB_Cnumroc( An, 0, Aimb1, Amb, myrow, Arow, nprow );
            Anq   = PB_Cnumroc( An, 0, Ainb1, Anb, mycol, Acol, npcol );
            if( ( Anp > 0 ) && ( Anq > 0 ) )
               Aptr = Mptr( A, Aii, Ajj, Ald, size );

            WBCld = WBCd[LLD_]; WBRld = WBRd[LLD_];

            if( upper )
            {
/*
*  sub( A ) is upper triangular
*/
               for( k = ( Astart / Alcmb ) * Alcmb; k >= 0; k -= Alcmb )
               {
                  ktmp = An - k; kb = MIN( ktmp, Alcmb );
/*
*  Solve logical diagonal block, WBC contains the solution scattered in multiple
*  process columns and WBR contains the solution replicated in the process rows.
*/
                  Akp = PB_Cnumroc( k, 0, Aimb1, Amb, myrow, Arow, nprow );
                  Akq = PB_Cnumroc( k, 0, Ainb1, Anb, mycol, Acol, npcol );
                  PB_Cptrsm( TYPE, WBRsum, SIDE, UPLO, TRANSA, DIAG, kb, nbb,
                             talph0, Aptr, k, k, Ad0, Mptr( WBC, Akp, 0, WBCld,
                             size ), WBCld, Mptr( WBR, 0, Akq, WBRld, size ),
                             WBRld );
/*
*  Update: only the part of sub( B ) to be solved at the next step is locally
*  updated and combined, the remaining part of the matrix to be solved later
*  is only locally updated.
*/
                  if( Akp > 0 )
                  {
                     Anq0 = PB_Cnumroc( kb, k, Ainb1, Anb, mycol, Acol, npcol );
                     if( WBCsum )
                     {
                        kbprev = MIN( k, Alcmb );
                        ktmp   = PB_Cnumroc( kbprev, k-kbprev, Aimb1, Amb,
                                             myrow, Arow, nprow );
                        Akp   -= ktmp;

                        if( ktmp > 0 )
                        {
                           if( Anq0 > 0 )
                              gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &ktmp,
                                    &nbb, &Anq0, negone, Mptr( Aptr, Akp, Akq,
                                    Ald, size ), &Ald, Mptr( WBR, 0, Akq, WBRld,
                                    size ), &WBRld, talph0, Mptr( WBC, Akp, 0,
                                    WBCld, size ), &WBCld );
                           Asrc = PB_Cindxg2p( k-1, Ainb1, Anb, Acol, Acol,
                                               npcol );
                           gsum2d( ctxt, ROW, &top, ktmp, nbb, Mptr( WBC, Akp,
                                   0, WBCld, size ), WBCld, myrow, Asrc );
                           if( mycol != Asrc )
                              pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &ktmp,
                                   &nbb, &izero, zero, zero, Mptr( WBC, Akp, 0,
                                   WBCld, size ), &WBCld );
                        }
                        if( ( Akp > 0 ) && ( Anq0 > 0 ) )
                           gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &Akp,
                                 &nbb, &Anq0, negone, Mptr( Aptr, 0, Akq, Ald,
                                 size ), &Ald, Mptr( WBR, 0, Akq, WBRld, size ),
                                 &WBRld, talph0, WBC, &WBCld );
                     }
                     else
                     {
                        if( Anq0 > 0 )
                           gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &Akp,
                                 &nbb, &Anq0, negone, Mptr( Aptr, 0, Akq, Ald,
                                 size ), &Ald, Mptr( WBR, 0, Akq, WBRld, size ),
                                 &WBRld, talph0, WBC, &WBCld );
                     }
                  }
                  talph0 = one;
               }
            }
            else
            {
/*
*  sub( A ) is lower triangular
*/
               for( k = 0; k < An; k += Alcmb )
               {
                  ktmp = An - k; kb = MIN( ktmp, Alcmb );
/*
*  Solve logical diagonal block, WBC contains the solution scattered in multiple
*  process columns and WBR contains the solution replicated in the process rows.
*/
                  Akp = PB_Cnumroc( k, 0, Aimb1, Amb, myrow, Arow, nprow );
                  Akq = PB_Cnumroc( k, 0, Ainb1, Anb, mycol, Acol, npcol );
                  PB_Cptrsm( TYPE, WBRsum, SIDE, UPLO, TRANSA, DIAG, kb, nbb,
                             talph0, Aptr, k, k, Ad0, Mptr( WBC, Akp, 0, WBCld,
                             size ), WBCld, Mptr( WBR, 0, Akq, WBRld, size ),
                             WBRld );
/*
*  Update: only the part of sub( B ) to be solved at the next step is locally
*  updated and combined, the remaining part of the matrix to be solved later is
*  only locally updated.
*/
                  Akp = PB_Cnumroc( k+kb, 0, Aimb1, Amb, myrow, Arow, nprow );
                  if( ( Anp0 = Anp - Akp ) > 0 )
                  {
                     Anq0 = PB_Cnumroc( kb, k, Ainb1, Anb, mycol, Acol, npcol );
                     if( WBCsum )
                     {
                        kbnext = ktmp - kb;
                        kbnext = MIN( kbnext, Alcmb );
                        ktmp   = PB_Cnumroc( kbnext, k+kb, Aimb1, Amb, myrow,
                                             Arow, nprow );
                        Anp0 -= ktmp;

                        if( ktmp > 0 )
                        {
                           if( Anq0 > 0 )
                              gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &ktmp,
                                    &nbb, &Anq0, negone, Mptr( Aptr, Akp, Akq,
                                    Ald, size ), &Ald, Mptr( WBR, 0, Akq, WBRld,
                                    size ), &WBRld, talph0, Mptr( WBC, Akp, 0,
                                    WBCld, size ), &WBCld );
                           Asrc = PB_Cindxg2p( k+kb, Ainb1, Anb, Acol, Acol,
                                               npcol );
                           gsum2d( ctxt, ROW, &top, ktmp, nbb, Mptr( WBC, Akp,
                                   0, WBCld, size ), WBCld, myrow, Asrc );
                           if( mycol != Asrc )
                              pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &ktmp,
                                   &nbb, &izero, zero, zero, Mptr( WBC, Akp, 0,
                                   WBCld, size ), &WBCld );
                        }
                        if( ( Anp0 > 0 ) && ( Anq0 > 0 ) )
                           gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &Anp0,
                                 &nbb, &Anq0, negone, Mptr( Aptr, Akp+ktmp, Akq,
                                 Ald, size ), &Ald, Mptr( WBR, 0, Akq, WBRld,
                                 size ), &WBRld, talph0, Mptr( WBC, Akp+ktmp, 0,
                                 WBCld, size ), &WBCld );
                     }
                     else
                     {
                        if( Anq0 > 0 )
                           gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &Anp0,
                                 &nbb, &Anq0, negone, Mptr( Aptr, Akp, Akq, Ald,
                                 size ), &Ald, Mptr( WBR, 0, Akq, WBRld, size ),
                                 &WBRld, talph0, Mptr( WBC, Akp, 0, WBCld,
                                 size ), &WBCld );
                     }
                  }
                  talph0 = one;
               }
            }
/*
*  Combine the scattered resulting matrix WBC
*/
            if( WBCsum && ( Anp > 0 ) )
               gsum2d( ctxt, ROW, &top, Anp, nbb, WBC, WBCld, myrow,
                       WBCd[CSRC_] );
/*
*  sub( B ) := WBC (if necessary)
*/
            if( WBCapbX )
               PB_Cpaxpby( TYPE, &conjg, An, nbb, one, WBC, 0, 0, WBCd, COLUMN,
                           zero, Bptr, 0, 0, DBUFB, &Broc );
         }
         else
         {
/*
*  Reuse sub( B ) and/or create vector WBR in process row owning the first or
*  last row of sub( A )
*/
            PB_CInOutV2( TYPE, &conjg, ROW, An, An, Astart, Ad0, nbb, Bptr,
                         0, 0, DBUFB, &Broc, &WBR, WBRd, &WBRfr, &WBRsum,
                         &WBRapbX );
/*
*  Create WBC in process columns spanned by sub( A )
*/
            PB_COutV( TYPE, COLUMN, INIT, An, An, Ad0, nbb, &WBC, WBCd, &WBCfr,
                      &WBCsum );
/*
*  Retrieve local quantities related to Ad0 -> sub( A )
*/
            Aimb1 = Ad0[IMB_ ]; Ainb1 = Ad0[INB_ ];
            Amb   = Ad0[MB_  ]; Anb   = Ad0[NB_  ];
            Arow  = Ad0[RSRC_]; Acol  = Ad0[CSRC_]; Ald = Ad0[LLD_];

            Anp   = PB_Cnumroc( An, 0, Aimb1, Amb, myrow, Arow, nprow );
            Anq   = PB_Cnumroc( An, 0, Ainb1, Anb, mycol, Acol, npcol );
            if( ( Anp > 0 ) && ( Anq > 0 ) )
               Aptr = Mptr( A, Aii, Ajj, Ald, size );

            WBCld = WBCd[LLD_]; WBRld = WBRd[LLD_];

            if( upper )
            {
/*
*  sub( A ) is upper triangular
*/
               for( k = 0; k < An; k += Alcmb )
               {
                  ktmp = An - k; kb = MIN( ktmp, Alcmb );
/*
*  Solve logical diagonal block, WBR contains the solution scattered in multiple
*  process rows and WBC contains the solution replicated in the process columns.
*/
                  Akp = PB_Cnumroc( k, 0, Aimb1, Amb, myrow, Arow, nprow );
                  Akq = PB_Cnumroc( k, 0, Ainb1, Anb, mycol, Acol, npcol );
                  PB_Cptrsm( TYPE, WBCsum, SIDE, UPLO, TRANSA, DIAG, nbb, kb,
                             talph0, Aptr, k, k, Ad0, Mptr( WBC, Akp, 0, WBCld,
                             size ), WBCld, Mptr( WBR, 0, Akq, WBRld, size ),
                             WBRld );
/*
*  Update: only the part of sub( B ) to be solved at the next step is locally
*  updated and combined, the remaining part of the matrix to be solved later is
*  only locally updated.
*/
                  Akq = PB_Cnumroc( k+kb, 0, Ainb1, Anb, mycol, Acol, npcol );
                  if( ( Anq0 = Anq - Akq ) > 0 )
                  {
                     Anp0 = PB_Cnumroc( kb, k, Aimb1, Amb, myrow, Arow, nprow );
                     if( WBRsum )
                     {
                        kbnext = ktmp - kb;
                        kbnext = MIN( kbnext, Alcmb );
                        ktmp   = PB_Cnumroc( kbnext, k+kb, Ainb1, Anb, mycol,
                                             Acol, npcol );
                        Anq0  -= ktmp;

                        if( ktmp > 0 )
                        {
                           if( Anp0 > 0 )
                              gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &nbb,
                                    &ktmp, &Anp0, negone, Mptr( WBC, Akp, 0,
                                    WBCld, size ), &WBCld, Mptr( Aptr, Akp, Akq,
                                    Ald, size ), &Ald, talph0, Mptr( WBR, 0,
                                    Akq, WBRld, size ), &WBRld );
                           Asrc = PB_Cindxg2p( k+kb, Aimb1, Amb, Arow, Arow,
                                               nprow );
                           gsum2d( ctxt, COLUMN, &top, nbb, ktmp, Mptr( WBR, 0,
                                   Akq, WBRld, size ), WBRld, Asrc, mycol );
                           if( myrow != Asrc )
                              pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &nbb,
                                   &ktmp, &izero, zero, zero, Mptr( WBR, 0, Akq,
                                   WBRld, size ), &WBRld );
                        }
                        if( ( Anp0 > 0 ) && ( Anq0 > 0 ) )
                           gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &nbb,
                                 &Anq0, &Anp0, negone, Mptr( WBC, Akp, 0, WBCld,
                                 size ), &WBCld, Mptr( Aptr, Akp, Akq+ktmp, Ald,
                                 size ), &Ald, talph0, Mptr( WBR, 0, Akq+ktmp,
                                 WBRld, size ), &WBRld );
                     }
                     else
                     {
                        if( Anp0 > 0 )
                           gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &nbb,
                                 &Anq0, &Anp0, negone, Mptr( WBC, Akp, 0, WBCld,
                                 size ), &WBCld, Mptr( Aptr, Akp, Akq, Ald,
                                 size ), &Ald, talph0, Mptr( WBR, 0, Akq, WBRld,
                                 size ), &WBRld );
                     }
                  }
                  talph0 = one;
               }
            }
            else
            {
/*
*  sub( A ) is lower triangular
*/
               for( k = ( Astart / Alcmb ) * Alcmb; k >= 0; k -= Alcmb )
               {
                  ktmp = An - k; kb = MIN( ktmp, Alcmb );
/*
*  Solve logical diagonal block, WBR contains the solution scattered in multiple
*  process rows and WBC contains the solution replicated in the process columns.
*/
                  Akp = PB_Cnumroc( k, 0, Aimb1, Amb, myrow, Arow, nprow );
                  Akq = PB_Cnumroc( k, 0, Ainb1, Anb, mycol, Acol, npcol );
                  PB_Cptrsm( TYPE, WBCsum, SIDE, UPLO, TRANSA, DIAG, nbb, kb,
                             talph0, Aptr, k, k, Ad0, Mptr( WBC, Akp, 0, WBCld,
                             size ), WBCld, Mptr( WBR, 0, Akq, WBRld, size ),
                             WBRld );
/*
*  Update: only the part of sub( B ) to be solved at the next step is locally
*  updated and combined, the remaining part of the matrix to be solved later
*  is only locally updated.
*/
                  if( Akq > 0 )
                  {
                     Anp0 = PB_Cnumroc( kb, k, Aimb1, Amb, myrow, Arow, nprow );
                     if( WBRsum )
                     {
                        kbprev = MIN( k, Alcmb );
                        ktmp   = PB_Cnumroc( kbprev, k-kbprev, Ainb1, Anb,
                                             mycol, Acol, npcol );
                        Akq   -= ktmp;

                        if( ktmp > 0 )
                        {
                           if( Anp0 > 0 )
                              gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &nbb,
                                    &ktmp, &Anp0, negone, Mptr( WBC, Akp, 0,
                                    WBCld, size ), &WBCld, Mptr( Aptr, Akp, Akq,
                                    Ald, size ), &Ald, talph0, Mptr( WBR, 0,
                                    Akq, WBRld, size ), &WBRld );
                           Asrc = PB_Cindxg2p( k-1, Aimb1, Amb, Arow, Arow,
                                               nprow );
                           gsum2d( ctxt, COLUMN, &top, nbb, ktmp, Mptr( WBR, 0,
                                   Akq, WBRld, size ), WBRld, Asrc, mycol );
                           if( myrow != Asrc )
                              pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &nbb,
                                   &ktmp, &izero, zero, zero, Mptr( WBR, 0, Akq,
                                   WBRld, size ), &WBRld );
                        }
                        if( ( Anp0 > 0 ) && ( Akq > 0 ) )
                           gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &nbb,
                                 &Akq, &Anp0, negone, Mptr( WBC, Akp, 0, WBCld,
                                 size ), &WBCld, Mptr( Aptr, Akp, 0, Ald,
                                 size ), &Ald, talph0, WBR, &WBRld );
                     }
                     else
                     {
                        if( Anp0 > 0 )
                           gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &nbb,
                                 &Akq, &Anp0, negone, Mptr( WBC, Akp, 0, WBCld,
                                 size ), &WBCld, Mptr( Aptr, Akp, 0, Ald,
                                 size ), &Ald, talph0, WBR, &WBRld );
                     }
                  }
                  talph0 = one;
               }
            }
/*
*  Combine the scattered resulting matrix WBR
*/
            if( WBRsum && ( Anq > 0 ) )
               gsum2d( ctxt, COLUMN, &top, nbb, Anq, WBR, WBRld, WBRd[RSRC_],
                       mycol );
/*
*  sub( B ) := WBR (if necessary)
*/
            if( WBRapbX )
               PB_Cpaxpby( TYPE, &conjg, nbb, An, one, WBR, 0, 0, WBRd, ROW,
                           zero, Bptr, 0, 0, DBUFB, &Broc );
         }

         if( WBCfr ) free( WBC );
         if( WBRfr ) free( WBR );
/*
*  Go to the next contiguous panel if any residing in this process row or column
*/
         BnpR -= nbb;

         if( BisR || ( BmyprocR == BcurrocR ) ) BiiR += nbb;
      }
/*
*  Go to next or previous process row or column owning some of sub( B )
*/
      if( !( BisR ) )
         p = ( Bfwd ? MModAdd1( p, BnprocsR ) : MModSub1( p, BnprocsR ) );
   }

   if( TranOp == CCOTRAN ) free( talpha );
/*
*  End of PB_CptrsmB
*/
}
