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
void PB_CptrmmB( PBTYP_T * TYPE, char * DIRECB, char * SIDE,
                 char * UPLO, char * TRANSA, char * DIAG, int M, int N,
                 char * ALPHA, char * A, int IA, int JA, int * DESCA,
                 char * B, int IB, int JB, int * DESCB )
#else
void PB_CptrmmB( TYPE, DIRECB, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A,
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
*  PB_CptrmmB  performs one of the matrix-matrix operations
*
*     sub( B ) := alpha * op( sub( A ) ) * sub( B ),
*
*     or
*
*     sub( B ) := alpha * sub( B ) * op( sub( A ) ),
*
*  where
*
*     sub( A ) denotes A(IA:IA+M-1,JA:JA+M-1)  if SIDE = 'L',
*                      A(IA:IA+N-1,JA:JA+N-1)  if SIDE = 'R', and,
*
*     sub( B ) denotes B(IB:IB+M-1,JB:JB+N-1).
*
*  Alpha  is a scalar,  sub( B )  is an m by n submatrix,  sub( A ) is a
*  unit, or non-unit, upper or lower triangular submatrix and op( X ) is
*  one of
*
*     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ).
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
*          On entry,  SIDE  specifies whether  op( sub( A ) ) multiplies
*          sub( B ) from the left or right as follows:
*
*          SIDE = 'L' or 'l'  sub( B ) := alpha*op( sub( A ) )*sub( B ),
*
*          SIDE = 'R' or 'r'  sub( B ) := alpha*sub( B )*op( sub( A ) ).
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
*          the local entries of the m by n transformed submatrix.
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
   char           Broc, GemmTa, GemmTb, TranOp, WBroc, WCroc, conjg, * one,
                  * talpha, * tbeta, top, * zero;
   int            Acol, Aii, Aimb1, Ainb1, Ajj, Alcmb, Ald, Alp, Alp0, Alq,
                  Alq0, Amb, Amp, An, Anb, Anq, Arow, BcurrocR, Bfwd, BiiD,
                  BiiR, Binb1D, Binb1R, BisR, Bld, BmyprocD, BmyprocR, BnD,
                  BnR, BnbD, BnbR, BnpR, BnprocsD, BnprocsR, BrocD, BrocR,
                  BsrcR, LNorRT, WBfr, WBld, WCfr, WCld, WCpbY, WCsum, ctxt,
                  l, lb, lside, ltmp, mycol, myrow, n, nb, nbb, notran, npcol,
                  nprow, p=0, size, tmp, upper;
   GEMM_T         gemm;
   GSUM2D_T       gsum2d;
/*
*  .. Local Arrays ..
*/
   int            Ad0[DLEN_], DBUFB[DLEN_], WCd[DLEN_], WBd[DLEN_];
   char           * Aptr = NULL, * Bptr = NULL, * WB = NULL, * WC = NULL;
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

   size   = TYPE->size;  one    = TYPE->one;     zero   = TYPE->zero;
   gemm   = TYPE->Fgemm; gsum2d = TYPE->Cgsum2d;
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

   Amp = PB_Cnumroc( An, 0, Aimb1, Amb, myrow, Arow, nprow );
   Anq = PB_Cnumroc( An, 0, Ainb1, Anb, mycol, Acol, npcol );
   if( ( Amp > 0 ) && ( Anq > 0 ) ) Aptr = Mptr( A, Aii, Ajj, Ald, size );
/*
*  Compute conjugate of alpha for the conjugate transpose cases
*/
   if( TranOp == CCOTRAN )
   {
      conjg = CCONJG; talpha = PB_Cmalloc( size );
      PB_Cconjg( TYPE, ALPHA, talpha );
   }
   else { conjg  = CNOCONJG; talpha = ALPHA; }
/*
*  Retrieve BLACS combine topology, set the transpose parameters to be passed
*  to the BLAS matrix multiply routine and finally describe the form of the
*  input and output operands.
*/
   if( LNorRT )
   {
      top    = *PB_Ctop( &ctxt, COMBINE, ROW, TOP_GET );
      GemmTa = CNOTRAN; GemmTb = ( lside ? CTRAN : TranOp );
      WCroc  = CCOLUMN; WBroc  = CROW;
   }
   else
   {
      top    = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
      GemmTb = CNOTRAN; GemmTa = ( lside ? TranOp : CTRAN );
      WCroc  = CROW;    WBroc = CCOLUMN;
   }
/*
*  Computational partitioning size is computed as the product of the logical
*  value returned by pilaenv_ and 2 * lcm( nprow, npcol ).
*/
   Alcmb = 2 * nb * PB_Clcm( ( Arow >= 0 ? nprow : 1 ),
                             ( Acol >= 0 ? npcol : 1 ) );
/*
*  When sub( B ) is not replicated and backward pass on sub( B ), find the
*  virtual process p owning the last row or column of sub( B ).
*/
   if( !( BisR = ( ( BsrcR < 0 ) || ( BnprocsR == 1 ) ) ) && !Bfwd )
   {
      tmp = PB_Cindxg2p( BnR-1, Binb1R, BnbR, BrocR, BrocR, BnprocsR );
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
/*
*  Replicate this panel in the process rows or columns spanned by sub( A ): WB
*/
         PB_CInV( TYPE, NOCONJG, &WBroc, An, An, Ad0, nbb, Bptr, 0, 0, DBUFB,
                  &Broc, &WB, WBd, &WBfr );
/*
*  Reuse sub( B ) and/or create vector WC in process columns or rows spanned by
*  sub( A )
*/
         PB_CInOutV( TYPE, &WCroc, An, An, Ad0, nbb, one, Bptr, 0, 0, DBUFB,
                     &Broc, &tbeta, &WC, WCd, &WCfr, &WCsum, &WCpbY );
/*
*  When the input data is first transposed, zero it now for later accumulation
*/
         if( notran )
            PB_Cplapad( TYPE, ALL, NOCONJG, DBUFB[M_], DBUFB[N_], zero, zero,
                        Bptr, 0, 0, DBUFB );
/*
*  Local matrix-matrix multiply iff I own some data
*/
         Aimb1 = Ad0[IMB_ ]; Ainb1 = Ad0[INB_ ]; Amb = Ad0[MB_]; Anb = Ad0[NB_];
         Arow  = Ad0[RSRC_]; Acol  = Ad0[CSRC_];
         Amp   = PB_Cnumroc( An, 0, Aimb1, Amb, myrow, Arow, nprow );
         Anq   = PB_Cnumroc( An, 0, Ainb1, Anb, mycol, Acol, npcol );

         WCld = WCd[LLD_];

         if( ( Amp > 0 ) && ( Anq > 0 ) )
         {
            WBld = WBd[LLD_];

            if( upper )
            {
/*
*  sub( A ) is upper triangular
*/
               if( LNorRT )
               {
                  for( l = 0; l < An; l += Alcmb )
                  {
                     lb  = An - l; lb = MIN( lb, Alcmb );
                     Alp = PB_Cnumroc( l, 0, Aimb1, Amb, myrow, Arow, nprow );
                     Alq = PB_Cnumroc( l, 0, Ainb1, Anb, mycol, Acol, npcol );
                     if( Alp > 0 )
                     {
                        Alq0 = PB_Cnumroc( lb, l, Ainb1, Anb, mycol, Acol,
                                           npcol );
                        gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( &GemmTb ), &Alp,
                              &nbb, &Alq0, talpha, Mptr( Aptr, 0, Alq, Ald,
                              size ), &Ald, Mptr( WB, 0, Alq, WBld, size ),
                              &WBld, one, WC, &WCld );
                     }
                     PB_Cptrm( TYPE, TYPE, SIDE, UPLO, TRANSA, DIAG, lb, nbb,
                               talpha, Aptr, l, l, Ad0, Mptr( WB, 0, Alq, WBld,
                               size ), WBld, Mptr( WC, Alp, 0, WCld, size ),
                               WCld, PB_Ctztrmm );
                  }
               }
               else
               {
                  for( l = 0; l < An; l += Alcmb )
                  {
                     lb   = An - l; lb = MIN( lb, Alcmb );
                     Alp  = PB_Cnumroc( l,  0, Aimb1, Amb, myrow, Arow, nprow );
                     Alq  = PB_Cnumroc( l,  0, Ainb1, Anb, mycol, Acol, npcol );
                     Alq0 = PB_Cnumroc( lb, l, Ainb1, Anb, mycol, Acol, npcol );
                     if( Alq0 > 0 )
                        gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( &GemmTb ), &nbb,
                              &Alq0, &Alp, talpha, WB, &WBld, Mptr( Aptr, 0,
                              Alq, Ald, size ), &Ald, one, Mptr( WC, 0, Alq,
                              WCld, size ), &WCld );
                     PB_Cptrm( TYPE, TYPE, SIDE, UPLO, TRANSA, DIAG, lb, nbb,
                               talpha, Aptr, l, l, Ad0, Mptr( WB, Alp, 0, WBld,
                               size ), WBld, Mptr( WC, 0, Alq, WCld, size ),
                               WCld, PB_Ctztrmm );
                  }
               }
            }
            else
            {
/*
*  sub( A ) is lower triangular
*/
               if( LNorRT )
               {
                  for( l = 0; l < An; l += Alcmb )
                  {
                     lb   = An - l; ltmp = l + ( lb = MIN( lb, Alcmb ) );
                     Alp  = PB_Cnumroc( l, 0, Aimb1, Amb, myrow, Arow, nprow );
                     Alq  = PB_Cnumroc( l, 0, Ainb1, Anb, mycol, Acol, npcol );
                     PB_Cptrm( TYPE, TYPE, SIDE, UPLO, TRANSA, DIAG, lb, nbb,
                               talpha, Aptr, l, l, Ad0, Mptr( WB, 0, Alq, WBld,
                               size ), WBld, Mptr( WC, Alp, 0, WCld, size ),
                               WCld, PB_Ctztrmm );
                     Alp  = PB_Cnumroc( ltmp, 0, Aimb1, Amb, myrow, Arow,
                                        nprow );
                     Alp0 = Amp - Alp;
                     Alq0 = PB_Cnumroc( lb,   l, Ainb1, Anb, mycol, Acol,
                                        npcol );
                     if( Alp0 > 0 )
                        gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( &GemmTb ), &Alp0,
                              &nbb, &Alq0, talpha, Mptr( Aptr, Alp, Alq, Ald,
                              size ), &Ald, Mptr( WB, 0, Alq, WBld, size ),
                              &WBld, one, Mptr( WC, Alp, 0, WCld, size ),
                              &WCld );
                  }
               }
               else
               {
                  for( l = 0; l < An; l += Alcmb )
                  {
                     lb   = An - l; ltmp = l + ( lb = MIN( lb, Alcmb ) );
                     Alp  = PB_Cnumroc( l, 0, Aimb1, Amb, myrow, Arow, nprow );
                     Alq  = PB_Cnumroc( l, 0, Ainb1, Anb, mycol, Acol, npcol );
                     PB_Cptrm( TYPE, TYPE, SIDE, UPLO, TRANSA, DIAG, lb, nbb,
                               talpha, Aptr, l, l, Ad0, Mptr( WB, Alp, 0, WBld,
                               size ), WBld, Mptr( WC, 0, Alq, WCld, size ),
                               WCld, PB_Ctztrmm );
                     Alp  = PB_Cnumroc( ltmp, 0, Aimb1, Amb, myrow, Arow,
                                        nprow );
                     Alp0 = Amp - Alp;
                     Alq0 = PB_Cnumroc( lb,   l, Ainb1, Anb, mycol, Acol,
                                        npcol );
                     if( Alq0 > 0 )
                        gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( &GemmTb ), &nbb,
                              &Alq0, &Alp0, talpha, Mptr( WB, Alp, 0, WBld,
                              size ), &WBld, Mptr( Aptr, Alp, Alq, Ald, size ),
                              &Ald, one, Mptr( WC, 0, Alq, WCld, size ),
                              &WCld );
                  }
               }
            }
         }
         if( WBfr ) free( WB );

         if( LNorRT )
         {
/*
*  Combine the partial column results into WC
*/
            if( WCsum && ( Amp > 0 ) )
               gsum2d( ctxt, ROW, &top, Amp, nbb, WC, WCld, myrow, WCd[CSRC_] );
/*
*  sub( B ) := WC (if necessary)
*/
            if( WCpbY )
               PB_Cpaxpby( TYPE, &conjg, An, nbb, one, WC, 0, 0, WCd, &WCroc,
                           zero, Bptr, 0, 0, DBUFB, &Broc );
         }
         else
         {
/*
*  Combine the partial row results into WC
*/
            if( WCsum && ( Anq > 0 ) )
               gsum2d( ctxt, COLUMN, &top, nbb, Anq, WC, WCld, WCd[RSRC_],
                       mycol );
/*
*  sub( B ) := WC (if necessary)
*/
            if( WCpbY )
               PB_Cpaxpby( TYPE, &conjg, nbb, An, one, WC, 0, 0, WCd, &WCroc,
                           zero, Bptr, 0, 0, DBUFB, &Broc );
         }
         if( WCfr ) free( WC );
/*
*  Go to the next contiguous panel if any residing in this process row or column
*/
         BnpR -= nbb;

         if( BisR || ( BmyprocR == BcurrocR ) ) BiiR += nbb;
      }
/*
*  Go to next or previous process row or column owning some of sub( B )
*/
      if( !BisR )
         p = ( Bfwd ? MModAdd1( p, BnprocsR ) : MModSub1( p, BnprocsR ) );
   }

   if( TranOp == CCOTRAN ) free( talpha );
/*
*  End of PB_CptrmmB
*/
}
