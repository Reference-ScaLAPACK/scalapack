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
void PB_CpgemmAB( PBTYP_T * TYPE, char * DIRECA, char * DIRECB,
                  char * TRANSA, char * TRANSB, Int M, Int N, Int K,
                  char * ALPHA, char * A, Int IA, Int JA, Int * DESCA,
                  char * B, Int IB, Int JB, Int * DESCB, char * BETA,
                  char * C, Int IC, Int JC, Int * DESCC )
#else
void PB_CpgemmAB( TYPE, DIRECA, DIRECB, TRANSA, TRANSB, M, N, K, ALPHA,
                  A, IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC,
                  DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * DIRECA, * DIRECB, * TRANSA, * TRANSB;
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
*  PB_CpgemmAB  performs one of the matrix-matrix operations
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
*  This is the  outer-product  algorithm using  the  logical  LCM hybrid
*  algorithmic blocking technique.  The submatrix operand sub( C ) stays
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
*  DIRECA  (global input) pointer to CHAR
*          On entry,  DIRECA  specifies  the direction in which the rows
*          or columns of sub( A ) should be looped over as follows:
*             DIRECA = 'F' or 'f'   forward  or increasing,
*             DIRECA = 'B' or 'b'   backward or decreasing.
*
*  DIRECB  (global input) pointer to CHAR
*          On entry,  DIRECB  specifies  the direction in which the rows
*          or columns of sub( B ) should be looped over as follows:
*             DIRECB = 'F' or 'f'   forward  or increasing,
*             DIRECB = 'B' or 'b'   backward or decreasing.
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
   char           Aroc, Broc, TrA, TrB, * one, * tbeta, * zero;
   Int            ABrocs, Abufld, AcurrocR, Afr, Afwd, AiD, AiR, AiiD, AiiR,
                  AinbD, AinbR, Ainb1D, Ainb1R, AisR, AkkR, Ald, AmyprocD,
                  AmyprocR, AnbD, AnbR, AnpD, AnpR, AnprocsD, AnprocsR, Aoff,
                  ArocD, ArocR, AsrcR, Bbufld, BcurrocR, Bfr, Bfwd, BiD, BiR,
                  BiiD, BiiR, BinbD, BinbR, Binb1D, Binb1R, BisR, BkkR, Bld,
                  BmyprocD, BmyprocR, BnbD, BnbR, BnpD, BnpR, BnprocsD,
                  BnprocsR, Boff, BrocD, BrocR, BsrcR, Ccol, Cii, Cimb1, Cinb1,
                  Cjj, Cld, Cmb, Cmp, Cnb, Cnq, Crow, WAfr, WAsum, WBfr, WBsum,
                  Wkbb=0, ctxt, k, kb, kbb, lcmb, maxp, maxpm1, maxq, mycol,
                  myrow, ncpq, nota, notb, npcol, npq=0, nprow, nrpq, p=0, q=0,
                  size, tmp;
   GEMM_T         gemm;
/*
*  .. Local Arrays ..
*/
   PB_VM_T        VM;
   Int            Cd0[DLEN_], DBUFA[DLEN_], DBUFB[DLEN_], WAd0[DLEN_],
                  WBd0[DLEN_];
   char           * Abuf = NULL, * Bbuf = NULL, * Cptr = NULL, * WA = NULL,
                  * WB   = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCC[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   nota = ( ( TrA = Mupcase( TRANSA[0] ) ) == CNOTRAN );
   notb = ( ( TrB = Mupcase( TRANSB[0] ) ) == CNOTRAN );
   TrA  = ( ( TrA == CCOTRAN ) ? CCONJG : CNOCONJG );
   TrB  = ( ( TrB == CCOTRAN ) ? CCONJG : CNOCONJG );

   size = TYPE->size;
/*
*  Retrieve local information for sub( A ), sub( B ) and sub( C )
*/
   if( nota )
   {
      AiR   = JA;          Aroc = CCOLUMN;     AnprocsR = npcol;
      AinbR = DESCA[INB_]; AnbR = DESCA[NB_ ]; AsrcR    = DESCA[CSRC_];
   }
   else
   {
      AiR   = IA;          Aroc = CROW;        AnprocsR = nprow;
      AinbR = DESCA[IMB_]; AnbR = DESCA[MB_ ]; AsrcR    = DESCA[RSRC_];
   }

   if( notb )
   {
      BiR   = IB;          Broc = CROW;        BnprocsR = nprow;
      BinbR = DESCB[IMB_]; BnbR = DESCB[MB_ ]; BsrcR    = DESCB[RSRC_];
   }
   else
   {
      BiR   = JB;          Broc = CCOLUMN;     BnprocsR = npcol;
      BinbR = DESCB[INB_]; BnbR = DESCB[NB_ ]; BsrcR    = DESCB[CSRC_];
   }
/*
*  Retrieve sub( C )'s local information: Aii, Ajj, Arow, Acol ...
*/
   PB_Cdescribe( M, N, IC, JC, DESCC, nprow, npcol, myrow, mycol, &Cii, &Cjj,
                 &Cld, &Cimb1, &Cinb1, &Cmb, &Cnb, &Crow, &Ccol, Cd0 );

   Cmp = PB_Cnumroc( M, 0, Cimb1, Cmb, myrow, Crow, nprow );
   Cnq = PB_Cnumroc( N, 0, Cinb1, Cnb, mycol, Ccol, npcol );
/*
*  When sub( A ) and sub( B ) do not span more than one process row or column,
*  there is no need to pack the data.
*/
   if( !( PB_Cspan( K, AiR, AinbR, AnbR, AsrcR, AnprocsR ) ) &&
       !( PB_Cspan( K, BiR, BinbR, BnbR, BsrcR, BnprocsR ) ) )
   {
      PB_CInV( TYPE, &TrA, COLUMN, M, N, Cd0, K, A, IA, JA, DESCA, &Aroc, &WA,
               WAd0, &WAfr );
      PB_CInV( TYPE, &TrB, ROW,    M, N, Cd0, K, B, IB, JB, DESCB, &Broc, &WB,
               WBd0, &WBfr );
      if( ( Cmp > 0 ) && ( Cnq > 0 ) )
      {
/*
*  Perform the local update if I own some of sub( C )
*/
         TYPE->Fgemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp, &Cnq, &K,
                      ALPHA, WA, &WAd0[LLD_], WB, &WBd0[LLD_], BETA, Mptr( C,
                      Cii, Cjj, Cld, size ), &Cld );
      }
      if( WAfr ) free( WA );
      if( WBfr ) free( WB );
      return;
   }
/*
*  sub( A ) and sub( B ) span more than one process row or column.
*/
   Afwd  = ( Mupcase( DIRECA[0] ) == CFORWARD );
   Bfwd  = ( Mupcase( DIRECB[0] ) == CFORWARD );

   one   = TYPE->one; zero = TYPE->zero; tbeta = BETA; gemm = TYPE->Fgemm;
   kb    = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );
/*
*  Compute local information for sub( A ) and sub( B )
*/
   if( nota )
   {
      AiD      = IA;          AinbD    = DESCA[IMB_]; AnbD     = DESCA[MB_];
      Ald      = DESCA[LLD_]; AmyprocD = myrow;       AmyprocR = mycol;
      AnprocsD = nprow;
      PB_Cinfog2l( IA, JA, DESCA, AnprocsD, AnprocsR, AmyprocD, AmyprocR,
                   &AiiD, &AiiR, &ArocD, &ArocR );
   }
   else
   {
      AiD      = JA;          AinbD    = DESCA[INB_]; AnbD     = DESCA[NB_];
      Ald      = DESCA[LLD_]; AmyprocD = mycol;       AmyprocR = myrow;
      AnprocsD = npcol;
      PB_Cinfog2l( IA, JA, DESCA, AnprocsR, AnprocsD, AmyprocR, AmyprocD,
                   &AiiR, &AiiD, &ArocR, &ArocD );
   }
   Ainb1D = PB_Cfirstnb( M, AiD, AinbD, AnbD );
   AnpD   = PB_Cnumroc( M, 0, Ainb1D, AnbD, AmyprocD, ArocD, AnprocsD );
   Ainb1R = PB_Cfirstnb( K, AiR, AinbR, AnbR );
   AisR   = ( ( AsrcR < 0 ) || ( AnprocsR == 1 ) );

   if( notb )
   {
      BiD      = JB;          BinbD    = DESCB[INB_]; BnbD     = DESCB[NB_];
      Bld      = DESCB[LLD_]; BmyprocD = mycol;       BmyprocR = myrow;
      BnprocsD = npcol;
      PB_Cinfog2l( IB, JB, DESCB, BnprocsR, BnprocsD, BmyprocR, BmyprocD,
                   &BiiR, &BiiD, &BrocR, &BrocD );
   }
   else
   {
      BiD      = IB;          BinbD    = DESCB[IMB_]; BnbD     = DESCB[MB_];
      Bld      = DESCB[LLD_]; BmyprocD = myrow;       BmyprocR = mycol;
      BnprocsD = nprow;
      PB_Cinfog2l( IB, JB, DESCB, BnprocsD, BnprocsR, BmyprocD, BmyprocR,
                   &BiiD, &BiiR, &BrocD, &BrocR );
   }
   Binb1D = PB_Cfirstnb( N, BiD, BinbD, BnbD );
   BnpD   = PB_Cnumroc( N, 0, Binb1D, BnbD, BmyprocD, BrocD, BnprocsD );
   Binb1R = PB_Cfirstnb( K, BiR, BinbR, BnbR );
   BisR   = ( ( BsrcR < 0 ) || ( BnprocsR == 1 ) );
/*
*  When sub( A ) is not replicated and backward pass on sub( A ), find the
*  virtual process q owning the last row or column of sub( A ).
*/
   if( !( AisR ) && !( Afwd ) )
   {
      tmp = PB_Cindxg2p( K - 1, Ainb1R, AnbR, ArocR, ArocR, AnprocsR );
      q   = MModSub( tmp, ArocR, AnprocsR );
   }
/*
*  When sub( B ) is not replicated and backward pass on sub( B ), find the
*  virtual process p owning the last row or column of sub( B ).
*/
   if( !( BisR ) && !( Bfwd ) )
   {
      tmp = PB_Cindxg2p( K - 1, Binb1R, BnbR, BrocR, BrocR, BnprocsR );
      p   = MModSub( tmp, BrocR, BnprocsR );
   }

   if( Cmp > 0 && Cnq > 0 ) Cptr = Mptr( C, Cii, Cjj, Cld, size );
/*
*  Allocate work space in process rows and columns spanned by sub( C )
*/
   PB_COutV( TYPE, COLUMN, NOINIT, M, N, Cd0, kb, &WA, WAd0, &WAfr, &WAsum );
   PB_COutV( TYPE, ROW,    NOINIT, M, N, Cd0, kb, &WB, WBd0, &WBfr, &WBsum );
/*
*  Loop over the virtual process grid induced by the sub( A ) and sub( B )
*/
   lcmb   = PB_Clcm( ( maxp = ( BisR ? 1 : BnprocsR ) ) * BnbR,
                     ( maxq = ( AisR ? 1 : AnprocsR ) ) * AnbR );
   maxpm1 = maxp - 1;
/*
*  Find out process coordinates corresponding to first virtual process (p,q)
*/
   AcurrocR = ( AisR ? -1 : MModAdd( ArocR, q, AnprocsR ) );
   AkkR     = PB_Cg2lrem( AiR, AinbR, AnbR, AcurrocR, AsrcR, AnprocsR );
   AnpR     = PB_Cnumroc( K, 0, Ainb1R, AnbR, AcurrocR, ArocR, AnprocsR );

   BcurrocR = ( BisR ? -1 : MModAdd( BrocR, p, BnprocsR ) );
   BkkR     = PB_Cg2lrem( BiR, BinbR, BnbR, BcurrocR, BsrcR, BnprocsR );
   BnpR     = PB_Cnumroc( K, 0, Binb1R, BnbR, BcurrocR, BrocR, BnprocsR );
/*
*  Find out how many diagonals this virtual process (p,q) has
*/
   PB_CVMinit( &VM, 0, BnpR, AnpR, Binb1R, Ainb1R, BnbR, AnbR, p, q,
               maxp, maxq, lcmb );
   npq = PB_CVMnpq( &VM );

   for( k = 0; k < K; k += kb )
   {
      kbb = K - k; kbb = MIN( kbb, kb );

      while( Wkbb != kbb )
      {
/*
*  Ensure that the current virtual process (p,q) has something to contribute
*  to the replicated buffers WA and WB.
*/
         while( npq == 0 )
         {
            if( ( Bfwd      && ( p == maxpm1 ) ) ||
                ( !( Bfwd ) && ( p == 0      ) ) )
               q = ( Afwd ? MModAdd1( q, maxq ) : MModSub1( q, maxq ) );
            p = ( Bfwd ? MModAdd1( p, maxp ) : MModSub1( p, maxp ) );

            AcurrocR = ( AisR ? -1 : MModAdd( ArocR, q, AnprocsR ) );
            AkkR     = PB_Cg2lrem(  AiR, AinbR,  AnbR, AcurrocR, AsrcR,
                                   AnprocsR );
            AnpR     = PB_Cnumroc( K, 0, Ainb1R, AnbR, AcurrocR, ArocR,
                                   AnprocsR );

            BcurrocR = ( BisR ? -1 : MModAdd( BrocR, p, BnprocsR ) );
            BkkR     = PB_Cg2lrem(  BiR, BinbR,  BnbR, BcurrocR, BsrcR,
                                   BnprocsR );
            BnpR     = PB_Cnumroc( K, 0, Binb1R, BnbR, BcurrocR, BrocR,
                                   BnprocsR );

            PB_CVMinit( &VM, 0, BnpR, AnpR, Binb1R, Ainb1R, BnbR, AnbR,
                        p, q, maxp, maxq, lcmb );
            npq = PB_CVMnpq( &VM );
         }
/*
*  Current virtual process (p,q) has something, find out how many rows or
*  columns could be used: ABrocs.
*/
         if( Wkbb == 0 ) { ABrocs = ( npq < kbb ? npq : kbb ); }
         else            { ABrocs = kbb - Wkbb; ABrocs = MIN( ABrocs, npq ); }
/*
*  Find out how many rows or columns of sub( A ) and sub( B ) are contiguous
*/
         PB_CVMcontig( &VM, &nrpq, &ncpq, &Boff, &Aoff );

         if( nota )
         {
/*
*  Compute the descriptor DBUFA for the buffer that will contained the packed
*  columns of sub( A ).
*/
            if( ( Afr = ( ncpq < ABrocs ) ) != 0 )
            {
/*
*  If columns of sub( A ) are not contiguous, then allocate the buffer and
*  pack the ABrocs columns of sub( A ).
*/
               Abufld = MAX( 1, AnpD );
               if( AisR || ( AmyprocR == AcurrocR ) )
               {
                  Abuf = PB_Cmalloc( AnpD * ABrocs * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, &Aroc, PACKING, NOTRAN,
                              ABrocs, AnpD, one, Mptr( A, AiiD, AkkR, Ald,
                              size ), Ald, zero, Abuf, Abufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( A ) directly.
*/
               Abufld = Ald;
               if( AisR || ( AmyprocR == AcurrocR ) )
                     Abuf = Mptr( A, AiiD, AkkR + Aoff, Ald, size );
            }
            PB_Cdescset( DBUFA, M, ABrocs, Ainb1D, ABrocs, AnbD, ABrocs,
                         ArocD, AcurrocR, ctxt, Abufld );
         }
         else
         {
/*
*  Compute the descriptor DBUFA for the buffer that will contained the packed
*  rows of sub( A ).
*/
            if( ( Afr = ( ncpq < ABrocs ) ) != 0 )
            {
/*
*  If rows of sub( A ) are not contiguous, then allocate the buffer and
*  pack the ABrocs rows of sub( A ).
*/
               Abufld = ABrocs;
               if( AisR || ( AmyprocR == AcurrocR ) )
               {
                  Abuf = PB_Cmalloc( AnpD * ABrocs * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, &Aroc, PACKING, NOTRAN,
                              ABrocs, AnpD, one, Mptr( A, AkkR, AiiD, Ald,
                              size ), Ald, zero, Abuf, Abufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( A ) directly.
*/
               Abufld = Ald;
               if( AisR || ( AmyprocR == AcurrocR ) )
                  Abuf = Mptr( A, AkkR + Aoff, AiiD, Ald, size );
            }
            PB_Cdescset( DBUFA, ABrocs, M, ABrocs, Ainb1D, ABrocs, AnbD,
                         AcurrocR, ArocD, ctxt, Abufld );
         }

         if( notb )
         {
/*
*  Compute the descriptor DBUFB for the buffer that will contained the packed
*  rows of sub( B ).
*/
            if( ( Bfr = ( nrpq < ABrocs ) ) != 0 )
            {
/*
*  If rows of sub( B ) are not contiguous, then allocate the buffer and
*  pack the ABrocs rows of sub( B ).
*/
               Bbufld = ABrocs;
               if( BisR || ( BmyprocR == BcurrocR ) )
               {
                  Bbuf = PB_Cmalloc( BnpD * ABrocs * size );
                  PB_CVMpack( TYPE, &VM, ROW,    &Broc, PACKING, NOTRAN,
                              ABrocs, BnpD, one, Mptr( B, BkkR, BiiD, Bld,
                              size ), Bld, zero, Bbuf, Bbufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( B ) directly.
*/
               Bbufld = Bld;
               if( BisR || ( BmyprocR == BcurrocR ) )
                  Bbuf = Mptr( B, BkkR + Boff, BiiD, Bld, size );
            }
            PB_Cdescset( DBUFB, ABrocs, N, ABrocs, Binb1D, ABrocs, BnbD,
                         BcurrocR, BrocD, ctxt, Bbufld );
         }
         else
         {
/*
*  Compute the descriptor DBUFB for the buffer that will contained the packed
*  columns of sub( B ).
*/
            if( ( Bfr = ( nrpq < ABrocs ) ) != 0 )
            {
/*
*  If columns of sub( B ) are not contiguous, then allocate the buffer and
*  pack the ABrocs columns of sub( B ).
*/
               Bbufld = MAX( 1, BnpD );
               if( BisR || ( BmyprocR == BcurrocR ) )
               {
                  Bbuf = PB_Cmalloc( BnpD * ABrocs * size );
                  PB_CVMpack( TYPE, &VM, ROW,    &Broc, PACKING, NOTRAN,
                              ABrocs, BnpD, one, Mptr( B, BiiD, BkkR, Bld,
                              size ), Bld, zero, Bbuf, Bbufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( B ) directly.
*/
               Bbufld = Bld;
               if( BisR || ( BmyprocR == BcurrocR ) )
                  Bbuf = Mptr( B, BiiD, BkkR + Boff, Bld, size );
            }
            PB_Cdescset( DBUFB, N, ABrocs, Binb1D, ABrocs, BnbD, ABrocs,
                         BrocD, BcurrocR, ctxt, Bbufld );
         }
/*
*  Update the local indexes of sub( A ) and sub( B )
*/
         PB_CVMupdate( &VM, ABrocs, &BkkR, &AkkR );
/*
*  Replicate panels of rows or columns of sub( A ) and sub( B ) over sub( C )
*  -> WA, WB
*/
         PB_CInV2( TYPE, &TrA, COLUMN, M, N, Cd0, ABrocs, Abuf, 0, 0,
                   DBUFA, &Aroc, WA, Wkbb, WAd0 );
         PB_CInV2( TYPE, &TrB, ROW,    M, N, Cd0, ABrocs, Bbuf, 0, 0,
                   DBUFB, &Broc, WB, Wkbb, WBd0 );

         if( Afr & ( AisR || ( AmyprocR == AcurrocR ) ) )
            if( Abuf ) free( Abuf );
         if( Bfr & ( BisR || ( BmyprocR == BcurrocR ) ) )
            if( Bbuf ) free( Bbuf );
/*
*  ABrocs rows or columns of sub( A ) and sub( B ) have been replicated,
*  update the number of diagonals in this virtual process as well as the
*  number of rows or columns of sub( A ) and sub( B ) that are in WA, WB.
*/
         npq  -= ABrocs;
         Wkbb += ABrocs;
      }
/*
*  Perform local update
*/
      if( Cmp > 0 && Cnq > 0 )
      {
         gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp, &Cnq, &kbb,
               ALPHA, WA, &WAd0[LLD_], WB, &WBd0[LLD_], tbeta, Cptr, &Cld );
         tbeta = one;
      }

      Wkbb = 0;
   }

   if( WAfr ) free( WA );
   if( WBfr ) free( WB );
/*
*  End of PB_CpgemmAB
*/
}
