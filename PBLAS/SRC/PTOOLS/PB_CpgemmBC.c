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
void PB_CpgemmBC( PBTYP_T * TYPE, char * DIRECB, char * DIRECC,
                  char * TRANSA, char * TRANSB, Int M, Int N, Int K,
                  char * ALPHA, char * A, Int IA, Int JA, Int * DESCA,
                  char * B, Int IB, Int JB, Int * DESCB, char * BETA,
                  char * C, Int IC, Int JC, Int * DESCC )
#else
void PB_CpgemmBC( TYPE, DIRECB, DIRECC, TRANSA, TRANSB, M, N, K, ALPHA,
                  A, IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC,
                  DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * DIRECB, * DIRECC, * TRANSA, * TRANSB;
   Int            IA, IB, IC, JA, JB, JC, K, M, N;
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
*  PB_CpgemmBC  performs one of the matrix-matrix operations
*
*     sub( C ) := alpha*op( sub( A ) )*op( sub( B ) ) + beta*sub( C ),
*
*  where
*
*     sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),  and, op( X )  is one  of
*     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ).
*
*  Thus, op( sub( A ) ) denotes A(IA:IA+M-1,JA:JA+K-1)  if TRANSA = 'N',
*                               A(IA:IA+K-1,JA:JA+M-1)' if TRANSA = 'T',
*                        conjg(A(IA:IA+K-1,JA:JA+M-1)') if TRANSA = 'C',
*
*  and,  op( sub( B ) ) denotes B(IB:IB+K-1,JB:JB+N-1)  if TRANSB = 'N',
*                               B(IB:IB+N-1,JB:JB+K-1)' if TRANSB = 'T',
*                        conjg(B(IB:IB+N-1,JB:JB+K-1)') if TRANSB = 'C'.
*
*  Alpha and beta are scalars.  A, B and C are matrices;  op( sub( A ) )
*  is an  m by k submatrix,  op( sub( B ) )  is an  k by n submatrix and
*  sub( C ) is an m by n submatrix.
*
*  This is the inner-product algorithm using the logical LCM algorithmic
*  blocking technique. The submatrix operand sub( A ) stays in place.
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
*          On entry,  DIRECB  specifies  the direction in which the rows
*          or columns of sub( B ) should be looped over as follows:
*             DIRECB = 'F' or 'f'   forward  or increasing,
*             DIRECB = 'B' or 'b'   backward or decreasing.
*
*  DIRECC  (global input) pointer to CHAR
*          On entry,  DIRECC  specifies  the direction in which the rows
*          or columns of sub( C ) should be looped over as follows:
*             DIRECC = 'F' or 'f'   forward  or increasing,
*             DIRECC = 'B' or 'b'   backward or decreasing.
*
*  TRANSA  (global input) pointer to CHAR
*          On entry,  TRANSA  specifies the form of op( sub( A ) ) to be
*          used in the matrix multiplication as follows:
*
*             TRANSA = 'N' or 'n'   op( sub( A ) ) = sub( A ),
*             TRANSA = 'T' or 't'   op( sub( A ) ) = sub( A )',
*             TRANSA = 'C' or 'c'   op( sub( A ) ) = conjg( sub( A )' ).
*
*  TRANSB  (global input) pointer to CHAR
*          On entry,  TRANSB  specifies the form of op( sub( B ) ) to be
*          used in the matrix multiplication as follows:
*
*             TRANSB = 'N' or 'n'   op( sub( B ) ) = sub( B ),
*             TRANSB = 'T' or 't'   op( sub( B ) ) = sub( B )',
*             TRANSB = 'C' or 'c'   op( sub( B ) ) = conjg( sub( B )' ).
*
*  M       (global input) INTEGER
*          On entry,  M  specifies  the number of rows of the  submatrix
*          op( sub( A ) ) and of the submatrix sub( C ). M  must  be  at
*          least  zero.
*
*  N       (global input) INTEGER
*          On entry, N specifies the number of columns of the  submatrix
*          op( sub( B ) )  and  the  number of columns of the  submatrix
*          sub( C ). N must be at least zero.
*
*  K       (global input) INTEGER
*          On entry, K specifies the number of columns of the  submatrix
*          op( sub( A ) )  and  the  number of rows   of  the  submatrix
*          op( sub( B ) ). K must be at least  zero.
*
*  ALPHA   (global input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as zero then the local entries of the arrays  A and
*          B corresponding to the entries of  the  submatrices  sub( A )
*          and sub( B ) respectively need not be set on input.
*
*  A       (local input) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+K-1 ) when  TRANSA = 'N' or 'n', and is at
*          least  Lc( 1, JA+M-1 )  otherwise.  Before  entry, this array
*          contains the local entries of the matrix A.
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
*          at least Lc( 1, JB+N-1 ) when  TRANSB = 'N' or 'n', and is at
*          least Lc( 1, JB+K-1 )  otherwise.  Before  entry,  this array
*          contains the local entries of the matrix B.
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
*          the local entries of the matrix  C.
*          On exit, the entries of this array corresponding to the local
*          entries of the  submatrix  sub( C )  are  overwritten  by the
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
   char           Broc, GemmTa, GemmTb, TrA, TrB, * one, * talpha, * tbeta,
                  top, * zero;
   Int            Acol, Aii, Aimb1, Ainb1, Ajj, Ald, Am, Amb, Amp, An, Anb,
                  Anq, Arow, Bbufld, BcurrocR, Bfr, Bfwd, BiD, BiR, BiiD, BiiR,
                  BinbD, BinbR, Binb1D, Binb1R, BisR, Bkk, Bld, BmyprocD,
                  BmyprocR, BnbD, BnbR, BnpD, BnpR, BnprocsD, BnprocsR, Boff,
                  BrocD, BrocR, BsrcR, Bsrc_, Cbufld, Ccol, Ccurcol, Cfr, Cfwd,
                  Cii, Cimb, Cimb1, Cinb, Cinb1, CisR, Cjj, Ckk, Cld, Cmb, Cmp,
                  Cnb, Cnq, Coff, Crow, Csrc, WBfr, WCfr, WCsum, ctxt, lcmb,
                  maxp, maxpm1, maxq, mycol, myrow, n, nb, nbb, ncpq, nota,
                  notb, npcol, npq=0, nprow, nrpq, p=0, q=0, size, tmp;
   GEMM_T         gemm;
   GSUM2D_T       gsum2d;
/*
*  .. Local Arrays ..
*/
   Int            Ad0[DLEN_], DBUFB[DLEN_], DBUFC[DLEN_], WBd[DLEN_],
                  WCd[DLEN_];
   PB_VM_T        VM;
   char           * Aptr = NULL, * Bbuf = NULL, * Cbuf = NULL, * WB   = NULL,
                  * WC   = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   Cblacs_gridinfo( ( ctxt = DESCC[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   Bfwd = ( Mupcase( DIRECB[0] ) == CFORWARD );
   Cfwd = ( Mupcase( DIRECC[0] ) == CFORWARD );
   nota = ( ( TrA = Mupcase( TRANSA[0] ) ) == CNOTRAN );
   notb = ( ( TrB = Mupcase( TRANSB[0] ) ) == CNOTRAN );

   size = TYPE->size;  one    = TYPE->one;     zero = TYPE->zero;
   gemm = TYPE->Fgemm; gsum2d = TYPE->Cgsum2d;
   nb   = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );
/*
*  Compute local information for sub( A ), sub( B ) and sub( C )
*/
   if( notb )
   {
      BiD      = IB;           BiR      = JB;
      Bsrc_    = CSRC_;        Broc     = CCOLUMN;
      BinbD    = DESCB[IMB_ ]; BinbR    = DESCB[INB_];
      BnbD     = DESCB[MB_  ]; BnbR     = DESCB[NB_ ];
      BsrcR    = DESCB[Bsrc_]; Bld      = DESCB[LLD_];
      BmyprocD = myrow;        BnprocsD = nprow;
      BmyprocR = mycol;        BnprocsR = npcol;
      PB_Cinfog2l( IB, JB, DESCB, BnprocsD, BnprocsR, BmyprocD, BmyprocR,
                   &BiiD, &BiiR, &BrocD, &BrocR );
   }
   else
   {
      BiD      = JB;           BiR      = IB;
      Bsrc_    = RSRC_;        Broc     = CROW;
      BinbR    = DESCB[IMB_ ]; BinbD    = DESCB[INB_];
      BnbR     = DESCB[MB_  ]; BnbD     = DESCB[NB_ ];
      BsrcR    = DESCB[Bsrc_]; Bld      = DESCB[LLD_];
      BmyprocD = mycol;        BnprocsD = npcol;
      BmyprocR = myrow;        BnprocsR = nprow;
      PB_Cinfog2l( IB, JB, DESCB, BnprocsR, BnprocsD, BmyprocR, BmyprocD,
                   &BiiR, &BiiD, &BrocR, &BrocD );
   }
   Binb1D = PB_Cfirstnb( K, BiD, BinbD, BnbD );
   BnpD   = PB_Cnumroc( K, 0, Binb1D, BnbD, BmyprocD, BrocD, BnprocsD );
   Binb1R = PB_Cfirstnb( N, BiR, BinbR, BnbR );

   Cimb = DESCC[IMB_ ]; Cinb = DESCC[INB_];
   Cmb  = DESCC[MB_  ]; Cnb  = DESCC[NB_ ];
   Csrc = DESCC[CSRC_]; Cld  = DESCC[LLD_];
   PB_Cinfog2l( IC, JC, DESCC, nprow, npcol, myrow, mycol, &Cii, &Cjj,
                &Crow, &Ccol );
   Cimb1 = PB_Cfirstnb( M, IC, Cimb, Cmb );
   Cmp   = PB_Cnumroc( M, 0, Cimb1, Cmb, myrow, Crow, nprow );
   Cinb1 = PB_Cfirstnb( N, JC, Cinb, Cnb );
/*
*  Retrieve the BLACS combine topology, compute conjugate of alpha for the
*  conjugate transpose case and set the transpose parameters to be passed to
*  the BLAS matrix multiply routine.
*/
   if( nota )
   {
      Am     = M; An     = K;
      top    = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_GET );
      talpha = ALPHA; GemmTa = CNOTRAN; GemmTb = ( notb ? CTRAN : TrB );
   }
   else
   {
      Am     = K; An     = M;
      top    = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
      if( TrA == CCOTRAN )
      {
         talpha = PB_Cmalloc( size ); PB_Cconjg( TYPE, ALPHA, talpha );
         GemmTa = ( ( TrB == CCOTRAN ) ? CTRAN : CCOTRAN );
      }
      else
      {
         talpha = ALPHA;
         GemmTa = ( ( TrB == CCOTRAN ) ? CCOTRAN : CTRAN );
      }
      GemmTb = CNOTRAN;
   }
/*
*  Compute descriptor Ad0 for sub( A )
*/
   PB_Cdescribe( Am, An, IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aii, &Ajj,
                 &Ald, &Aimb1, &Ainb1, &Amb, &Anb, &Arow, &Acol, Ad0 );

   Amp = PB_Cnumroc( Am, 0, Aimb1, Amb, myrow, Arow, nprow );
   Anq = PB_Cnumroc( An, 0, Ainb1, Anb, mycol, Acol, npcol );
   if( ( Amp > 0 ) && ( Anq > 0 ) ) { Aptr = Mptr( A, Aii, Ajj, Ald, size ); }
/*
*  When sub( B ) is not replicated and backward pass on sub( B ), find the
*  virtual process q owning the last row or column of sub( B ).
*/
   if( !( BisR = ( ( BsrcR < 0 ) || ( BnprocsR == 1 ) ) ) && !Bfwd )
   {
      tmp = PB_Cindxg2p( N - 1, Binb1R, BnbR, BrocR, BrocR, BnprocsR );
      q   = MModSub( tmp, BrocR, BnprocsR );
   }
/*
*  When sub( C ) is not replicated and backward pass on sub( C ), find the
*  virtual process p owning the last row or column of sub( C ).
*/
   if( !( CisR = ( ( Ccol < 0 ) || ( npcol == 1 ) ) ) && !Cfwd )
   {
      tmp = PB_Cindxg2p( N - 1, Cinb1, Cnb, Ccol, Ccol, npcol );
      p   = MModSub( tmp, Ccol, npcol );
   }
/*
*  Loop over the virtual process grid induced by the rows or columns of
*  sub( B ) and sub( C ).
*/
   lcmb   = PB_Clcm( ( maxp = ( CisR ? 1 : npcol    ) ) * Cnb,
                     ( maxq = ( BisR ? 1 : BnprocsR ) ) * BnbR );
   n      = N;
   maxpm1 = maxp - 1;

   while( n > 0 )
   {
/*
*  Initialize local virtual matrix in process (p,q)
*/
      BcurrocR = ( BisR ? -1 : MModAdd( BrocR, q, BnprocsR ) );
      Bkk      = PB_Cg2lrem( BiR, BinbR, BnbR, BcurrocR, BsrcR, BnprocsR );
      BnpR     = PB_Cnumroc( N, 0, Binb1R, BnbR, BcurrocR, BrocR, BnprocsR );

      Ccurcol  = ( CisR ? -1 : MModAdd( Ccol,  p, npcol    ) );
      Ckk      = PB_Cg2lrem( JC, Cinb, Cnb, Ccurcol, Csrc, npcol );
      Cnq      = PB_Cnumroc( N, 0, Cinb1, Cnb, Ccurcol, Ccol, npcol );

      PB_CVMinit( &VM, 0, Cnq, BnpR, Cinb1, Binb1R, Cnb, BnbR, p, q,
                  maxp, maxq, lcmb );
/*
*  Find how many diagonals in this virtual process
*/
      npq = PB_CVMnpq( &VM );

      n -= npq;
/*
*  Re-adjust the number of rows or columns to be (un)packed, in order to
*  average the message sizes.
*/
      if( npq ) nbb = npq / ( ( npq - 1 ) / nb + 1 );

      while( npq )
      {
         nbb = MIN( nbb, npq );
/*
*  Find out how many rows or columns of sub( B ) and sub( C ) are contiguous
*/
         PB_CVMcontig( &VM, &nrpq, &ncpq, &Coff, &Boff );

         if( notb )
         {
/*
*  Compute the descriptor DBUFB for the buffer that will contained the packed
*  columns of sub( B ).
*/
            if( ( Bfr = ( ncpq < nbb ) ) != 0 )
            {
/*
*  If columns of sub( B ) are not contiguous, then allocate the buffer and
*  pack the nbb columns of sub( B ).
*/
               Bbufld = MAX( 1, BnpD );
               if( BisR || ( BmyprocR == BcurrocR ) )
               {
                  Bbuf   = PB_Cmalloc( BnpD * nbb * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, &Broc, PACKING, NOTRAN, nbb,
                              BnpD, one, Mptr( B, BiiD, Bkk,  Bld, size ), Bld,
                              zero, Bbuf, Bbufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( B ) directly.
*/
               Bbufld = Bld;
               if( BisR || ( BmyprocR == BcurrocR ) )
                  Bbuf = Mptr( B, BiiD, Bkk+Boff, Bld, size );
            }
            PB_Cdescset( DBUFB, K, nbb, Binb1D, nbb, BnbD, nbb, BrocD,
                         BcurrocR, ctxt, Bbufld );
         }
         else
         {
/*
*  Compute the descriptor DBUFB for the buffer that will contained the packed
*  rows of sub( B ).
*/
            if( ( Bfr = ( ncpq < nbb ) ) != 0 )
            {
/*
*  If rows of sub( B ) are not contiguous, then allocate the buffer and pack
*  the nbb rows of sub( B ).
*/
               Bbufld = nbb;
               if( BisR || ( BmyprocR == BcurrocR ) )
               {
                  Bbuf   = PB_Cmalloc( BnpD * nbb * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, &Broc, PACKING, NOTRAN, nbb,
                              BnpD, one, Mptr( B, Bkk,  BiiD, Bld, size ), Bld,
                              zero, Bbuf, Bbufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( B ) directly.
*/
               Bbufld = Bld;
               if( BisR || ( BmyprocR == BcurrocR ) )
                  Bbuf = Mptr( B, Bkk+Boff, BiiD, Bld, size );
            }
            PB_Cdescset( DBUFB, nbb, K, nbb, Binb1D, nbb, BnbD, BcurrocR,
                         BrocD, ctxt, Bbufld );
         }

         if( nota )
         {
/*
*  Replicate this panel of rows or columns of sub( B ) over sub( A ) -> WB
*/
            PB_CInV( TYPE, NOCONJG, ROW,    Am, An, Ad0, nbb, Bbuf, 0, 0,
                     DBUFB, &Broc, &WB, WBd, &WBfr );
/*
*  Allocate space for temporary results in scope of sub( A ) -> WC
*/
            PB_COutV( TYPE, COLUMN, INIT, Am, An, Ad0, nbb, &WC, WCd, &WCfr,
                      &WCsum );
/*
*  Local matrix-matrix multiply iff I own some data
*/
            if( Amp > 0 && Anq > 0 )
               gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( &GemmTb ), &Amp, &nbb,
                               &Anq, talpha, Aptr, &Ald, WB, &WBd[LLD_], zero,
                               WC, &WCd[LLD_] );
            if( WBfr ) free( WB );
            if( Bfr && ( BisR || ( BmyprocR == BcurrocR ) ) )
               if( Bbuf ) free( Bbuf );
/*
*  Accumulate the intermediate results in WC
*/
            if( WCsum )
            {
               WCd[CSRC_] = Ccurcol;
               if( Amp > 0 )
                  gsum2d( ctxt, ROW,    &top, Amp, nbb, WC, WCd[LLD_], myrow,
                          WCd[CSRC_] );
            }
/*
*  Compute the descriptor DBUFC for the buffer that will contained the packed
*  columns of sub( C ). Allocate it.
*/
            if( ( Cfr = ( nrpq < nbb ) ) != 0 )
            {
/*
*  If columns of sub( C ) are not contiguous, then allocate the buffer
*/
               Cbufld = MAX( 1, Cmp );  tbeta = zero;
               if( CisR || ( mycol == Ccurcol ) )
                  Cbuf = PB_Cmalloc( Cmp * nbb * size );
            }
            else
            {
/*
*  Otherwise re-use sub( C )
*/
               Cbufld = Cld;             tbeta = BETA;
               if( CisR || ( mycol == Ccurcol ) )
                  Cbuf = Mptr( C, Cii, Ckk+Coff, Cld, size );
            }
            PB_Cdescset( DBUFC, M, nbb, Cimb1, nbb, Cmb, nbb, Crow, Ccurcol,
                         ctxt, Cbufld );
/*
*  Cbuf := Cbuf + WC
*/
            PB_Cpaxpby( TYPE, NOCONJG, M, nbb, one, WC, 0, 0, WCd, COLUMN,
                        tbeta, Cbuf, 0, 0, DBUFC, COLUMN );
/*
*  Unpack the nbb columns of sub( C ) and release the buffer containing them.
*/
            if( Cfr && ( CisR || ( mycol == Ccurcol ) ) )
            {
               PB_CVMpack( TYPE, &VM, ROW, COLUMN, UNPACKING, NOTRAN, nbb, Cmp,
                           BETA, Mptr( C, Cii, Ckk, Cld, size ), Cld, one, Cbuf,
                           Cbufld );
               if( Cbuf ) free( Cbuf );
            }
            if( WCfr ) free( WC );
         }
         else
         {
/*
*  Replicate this panel of rows or columns of sub( B ) over sub( A ) -> WB
*/
            PB_CInV( TYPE, NOCONJG, COLUMN, Am, An, Ad0, nbb, Bbuf, 0, 0,
                     DBUFB, &Broc, &WB, WBd, &WBfr );
/*
*  Allocate space for temporary results in scope of sub( A ) -> WC
*/
            PB_COutV( TYPE, ROW,    INIT, Am, An, Ad0, nbb, &WC, WCd, &WCfr,
                      &WCsum );
/*
*  Local matrix-matrix multiply iff I own some data
*/
            if( Amp > 0 && Anq > 0 )
               gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( &GemmTb ), &nbb, &Anq,
                     &Amp, talpha, WB, &WBd[LLD_], Aptr, &Ald, zero, WC,
                     &WCd[LLD_] );
            if( WBfr ) free( WB );
            if( Bfr && ( BisR || ( BmyprocR == BcurrocR ) ) )
               if( Bbuf ) free( Bbuf );
/*
*  Accumulate the intermediate results in WC
*/
            if( WCsum )
            {
               WCd[RSRC_] = 0;
               if( Anq > 0 )
                  gsum2d( ctxt, COLUMN, &top, nbb, Anq, WC, WCd[LLD_],
                          WCd[RSRC_], mycol );
            }
/*
*  Compute the descriptor DBUFC for the buffer that will contained the packed
*  columns of sub( C ). Allocate it.
*/
            if( ( Cfr = ( nrpq < nbb ) ) != 0 )
            {
/*
*  If columns of sub( C ) are not contiguous, then allocate the buffer
*/
               Cbufld = MAX( 1, Cmp );  tbeta = zero;
               if( CisR || ( mycol == Ccurcol ) )
                  Cbuf = PB_Cmalloc( Cmp * nbb * size );
            }
            else
            {
/*
*  Otherwise re-use sub( C )
*/
               Cbufld = Cld;             tbeta = BETA;
               if( CisR || ( mycol == Ccurcol ) )
                  Cbuf = Mptr( C, Cii, Ckk+Coff, Cld, size );
            }
            PB_Cdescset( DBUFC, M, nbb, Cimb1, nbb, Cmb, nbb, Crow, Ccurcol,
                         ctxt, Cbufld );
/*
*  Cbuf := Cbuf + WC'
*/
            PB_Cpaxpby( TYPE, ( TrA == CCOTRAN ? CONJG : NOCONJG ), nbb, M,
                        one, WC, 0, 0, WCd, ROW,    tbeta, Cbuf, 0, 0, DBUFC,
                        COLUMN );
/*
*  Unpack the nbb columns of sub( C ) and release the buffer containing them.
*/
            if( Cfr && ( CisR || ( mycol == Ccurcol ) ) )
            {
               PB_CVMpack( TYPE, &VM, ROW, COLUMN, UNPACKING, NOTRAN, nbb, Cmp,
                           BETA, Mptr( C, Cii, Ckk, Cld, size ), Cld, one, Cbuf,
                           Cbufld );
               if( Cbuf ) free( Cbuf );
            }
            if( WCfr ) free( WC );
         }
/*
*  Update the local indexes of sub( B ) and sub( C )
*/
         PB_CVMupdate( &VM, nbb, &Ckk, &Bkk );

         npq -= nbb;
      }
/*
*  Go to next or previous virtual process row or column
*/
      if( ( Cfwd      && ( p == maxpm1 ) ) ||
          ( !( Cfwd ) && ( p == 0      ) ) )
         q = ( Bfwd ? MModAdd1( q, maxq ) : MModSub1( q, maxq ) );
      p = ( Cfwd ? MModAdd1( p, maxp ) : MModSub1( p, maxp ) );
   }

   if( TrA == CCOTRAN ) free( talpha );
/*
*  End of PB_CpgemmBC
*/
}
