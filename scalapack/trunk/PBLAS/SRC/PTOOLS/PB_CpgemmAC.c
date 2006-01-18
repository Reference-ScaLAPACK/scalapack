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
void PB_CpgemmAC( PBTYP_T * TYPE, char * DIRECA, char * DIRECC,
                  char * TRANSA, char * TRANSB, int M, int N, int K,
                  char * ALPHA, char * A, int IA, int JA, int * DESCA,
                  char * B, int IB, int JB, int * DESCB, char * BETA,
                  char * C, int IC, int JC, int * DESCC )
#else
void PB_CpgemmAC( TYPE, DIRECA, DIRECC, TRANSA, TRANSB, M, N, K, ALPHA,
                  A, IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC,
                  DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * DIRECA, * DIRECC, * TRANSA, * TRANSB;
   int            IA, IB, IC, JA, JB, JC, K, M, N;
   char           * ALPHA, * BETA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCB, * DESCC;
   char           * A, * B, * C;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CpgemmAC  performs one of the matrix-matrix operations
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
*  blocking technique. The submatrix operand sub( B ) stays in place.
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
   char           Aroc, GemmTa, GemmTb, TrA, TrB, * one, * talpha, * tbeta,
                  top, * zero;
   int            Abufld, AcurrocR, Afr, Afwd, AiD, AiR, AiiD, AiiR, AinbD,
                  AinbR, Ainb1D, Ainb1R, AisR, Akk, Ald, AmyprocD, AmyprocR,
                  AnbD, AnbR, AnpD, AnpR, AnprocsD, AnprocsR, Aoff, ArocD,
                  ArocR, AsrcR, Asrc_, Bcol, Bii, Bimb1, Binb1, Bjj, Bld, Bm,
                  Bmb, Bmp, Bn, Bnb, Bnq, Brow, Cbufld, Ccol, Ccurrow, Cfr,
                  Cfwd, Cii, Cimb, Cimb1, Cinb, Cinb1, CisR, Cjj, Ckk, Cld,
                  Cmb, Cmp, Cnb, Cnq, Coff, Crow, Csrc, WAfr, WCfr, WCsum,
                  ctxt, lcmb, m, maxp, maxpm1, maxq, mb, mbb, mycol, myrow,
                  ncpq, nota, notb, npcol, npq=0, nprow, nrpq, p=0, q=0, size,
                  tmp;
   GEMM_T         gemm;
   GSUM2D_T       gsum2d;
/*
*  .. Local Arrays ..
*/
   PB_VM_T        VM;
   int            Bd0[DLEN_], DBUFA[DLEN_], DBUFC[DLEN_], WAd[DLEN_],
                  WCd[DLEN_];
   char           * Abuf = NULL, * Bptr = NULL, * Cbuf = NULL, * WA = NULL,
                  * WC   = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   Cblacs_gridinfo( ( ctxt = DESCC[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   Afwd = ( Mupcase( DIRECA[0] ) == CFORWARD );
   Cfwd = ( Mupcase( DIRECC[0] ) == CFORWARD );
   nota = ( ( TrA = Mupcase( TRANSA[0] ) ) == CNOTRAN );
   notb = ( ( TrB = Mupcase( TRANSB[0] ) ) == CNOTRAN );

   size = TYPE->size; one  = TYPE->one; zero = TYPE->zero;
   gemm = TYPE->Fgemm; gsum2d = TYPE->Cgsum2d;
   mb   = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );
/*
*  Compute local information for sub( A ), sub( B ) and sub( C )
*/
   if( nota )
   {
      AiD      = JA;           AiR      = IA;
      Asrc_    = RSRC_;        Aroc     = CROW;
      AinbR    = DESCA[IMB_ ]; AinbD    = DESCA[INB_];
      AnbR     = DESCA[MB_  ]; AnbD     = DESCA[NB_ ];
      AsrcR    = DESCA[Asrc_]; Ald      = DESCA[LLD_];
      AmyprocD = mycol;        AnprocsD = npcol;
      AmyprocR = myrow;        AnprocsR = nprow;
      PB_Cinfog2l( IA, JA, DESCA, AnprocsR, AnprocsD, AmyprocR, AmyprocD,
                   &AiiR, &AiiD, &ArocR, &ArocD );
   }
   else
   {
      AiD      = IA;           AiR      = JA;
      Asrc_    = CSRC_;        Aroc     = CCOLUMN;
      AinbD    = DESCA[IMB_ ]; AinbR    = DESCA[INB_];
      AnbD     = DESCA[MB_  ]; AnbR     = DESCA[NB_ ];
      AsrcR    = DESCA[Asrc_]; Ald      = DESCA[LLD_];
      AmyprocD = myrow;        AnprocsD = nprow;
      AmyprocR = mycol;        AnprocsR = npcol;
      PB_Cinfog2l( IA, JA, DESCA, AnprocsD, AnprocsR, AmyprocD, AmyprocR,
                   &AiiD, &AiiR, &ArocD, &ArocR );
   }
   Ainb1D = PB_Cfirstnb( K, AiD, AinbD, AnbD );
   AnpD   = PB_Cnumroc( K, 0, Ainb1D, AnbD, AmyprocD, ArocD, AnprocsD );
   Ainb1R = PB_Cfirstnb( M, AiR, AinbR, AnbR );

   Cimb   = DESCC[IMB_ ]; Cinb = DESCC[INB_];
   Cmb    = DESCC[MB_  ]; Cnb  = DESCC[NB_ ];
   Csrc   = DESCC[RSRC_]; Cld  = DESCC[LLD_];
   PB_Cinfog2l( IC, JC, DESCC, nprow, npcol, myrow, mycol, &Cii, &Cjj,
                &Crow, &Ccol );
   Cimb1 = PB_Cfirstnb( M, IC, Cimb, Cmb );
   Cinb1 = PB_Cfirstnb( N, JC, Cinb, Cnb );
   Cnq   = PB_Cnumroc( N, 0, Cinb1, Cnb, mycol, Ccol, npcol );
/*
*  Retrieve the BLACS combine topology, compute conjugate of alpha for the
*  conjugate transpose case and set the transpose parameters to be passed to
*  the BLAS matrix multiply routine.
*/
   if( notb )
   {
      Bm     = K; Bn     = N;
      top    = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );
      talpha = ALPHA; GemmTa = ( nota ? CTRAN : TrA ); GemmTb = CNOTRAN;
   }
   else
   {
      Bm     = N; Bn     = K;
      top    = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_GET );
      if( TrB == CCOTRAN )
      {
         talpha = PB_Cmalloc( size ); PB_Cconjg( TYPE, ALPHA, talpha );
         GemmTb = ( ( TrA == CCOTRAN ) ? CTRAN : CCOTRAN );
      }
      else
      {
         talpha = ALPHA;
         GemmTb = ( ( TrA == CCOTRAN ) ? CCOTRAN : CTRAN );
      }
      GemmTa = CNOTRAN;
   }
/*
*  Compute descriptor Bd0 for sub( B )
*/
   PB_Cdescribe( Bm, Bn, IB, JB, DESCB, nprow, npcol, myrow, mycol, &Bii, &Bjj,
                 &Bld, &Bimb1, &Binb1, &Bmb, &Bnb, &Brow, &Bcol, Bd0 );

   Bmp = PB_Cnumroc( Bm, 0, Bimb1, Bmb, myrow, Brow, nprow );
   Bnq = PB_Cnumroc( Bn, 0, Binb1, Bnb, mycol, Bcol, npcol );
   if( ( Bmp > 0 ) && ( Bnq > 0 ) ) Bptr = Mptr( B, Bii, Bjj, Bld, size );
/*
*  When sub( A ) is not replicated and backward pass on sub( A ), find the
*  virtual process q owning the last row or column of sub( A ).
*/
   if( !( AisR = ( ( AsrcR < 0 ) || ( AnprocsR == 1 ) ) ) && !Afwd )
   {
      tmp = PB_Cindxg2p( M - 1, Ainb1R, AnbR, ArocR, ArocR, AnprocsR );
      q   = MModSub( tmp, ArocR, AnprocsR );
   }
/*
*  When sub( C ) is not replicated and backward pass on sub( C ), find the
*  virtual process p owning the last row or column of sub( C ).
*/
   if( !( CisR = ( ( Crow < 0 ) || ( nprow == 1 ) ) ) && !Cfwd )
   {
      tmp = PB_Cindxg2p( M - 1, Cimb1, Cmb, Crow, Crow, nprow );
      p   = MModSub( tmp, Crow, nprow );
   }
/*
*  Loop over the virtual process grid induced by the rows or columns of
*  sub( A ) and sub( C ).
*/
   lcmb   = PB_Clcm( ( maxp = ( CisR ? 1 : nprow    ) ) * Cmb,
                     ( maxq = ( AisR ? 1 : AnprocsR ) ) * AnbR );
   m      = M;
   maxpm1 = maxp - 1;

   while( m > 0 )
   {
/*
*  Initialize local virtual matrix in process (p,q)
*/
      AcurrocR = ( AisR ? -1 : MModAdd( ArocR, q, AnprocsR ) );
      Akk      = PB_Cg2lrem( AiR, AinbR, AnbR, AcurrocR, AsrcR, AnprocsR );
      AnpR     = PB_Cnumroc( M, 0, Ainb1R, AnbR, AcurrocR, ArocR, AnprocsR );

      Ccurrow  = ( CisR ? -1 : MModAdd( Crow,  p, nprow    ) );
      Ckk      = PB_Cg2lrem( IC, Cimb, Cmb, Ccurrow, Csrc, nprow );
      Cmp      = PB_Cnumroc( M, 0, Cimb1, Cmb, Ccurrow, Crow, nprow );

      PB_CVMinit( &VM, 0, Cmp, AnpR, Cimb1, Ainb1R, Cmb, AnbR, p, q,
                  maxp, maxq, lcmb );
/*
*  Find how many diagonals in this virtual process
*/
      npq = PB_CVMnpq( &VM );

      m  -= npq;
/*
*  Re-adjust the number of rows or columns to be (un)packed, in order to
*  average the message sizes.
*/
      if( npq ) mbb = npq / ( ( npq - 1 ) / mb + 1 );

      while( npq )
      {
         mbb = MIN( mbb, npq );
/*
*  Find out how many rows or columns of sub( A ) and sub( C ) are contiguous
*/
         PB_CVMcontig( &VM, &nrpq, &ncpq, &Coff, &Aoff );

         if( nota )
         {
/*
*  Compute the descriptor DBUFA for the buffer that will contained the packed
*  columns of sub( A ).
*/
            if( ( Afr = ( ncpq < mbb ) ) != 0 )
            {
/*
*  If rows of sub( A ) are not contiguous, then allocate the buffer and
*  pack the mbb rows of sub( A ).
*/
               Abufld = mbb;
               if( AisR || ( AmyprocR == AcurrocR ) )
               {
                  Abuf   = PB_Cmalloc( AnpD * mbb * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, &Aroc, PACKING, NOTRAN, mbb,
                              AnpD, one, Mptr( A, Akk, AiiD, Ald, size ), Ald,
                              zero, Abuf, Abufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( B ) directly.
*/
               Abufld = Ald;
               if( AisR || ( AmyprocR == AcurrocR ) )
                  Abuf = Mptr( A, Akk+Aoff, AiiD, Ald, size );
            }
            PB_Cdescset( DBUFA, mbb, K, mbb, Ainb1D, mbb, AnbD, AcurrocR,
                         ArocD, ctxt, Abufld );
         }
         else
         {
/*
*  Compute the descriptor DBUFA for the buffer that will contained the packed
*  columns of sub( A ).
*/
            if( ( Afr = ( ncpq < mbb ) ) != 0 )
            {
/*
*  If columns of sub( A ) are not contiguous, then allocate the buffer and pack
*  the mbb columns of sub( A ).
*/
               Abufld = MAX( 1, AnpD );
               if( AisR || ( AmyprocR == AcurrocR ) )
               {
                  Abuf   = PB_Cmalloc( AnpD * mbb * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, &Aroc, PACKING, NOTRAN, mbb,
                              AnpD, one, Mptr( A, AiiD, Akk, Ald, size ), Ald,
                              zero, Abuf, Abufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( A ) directly.
*/
               Abufld = Ald;
               if( AisR || ( AmyprocR == AcurrocR ) )
                  Abuf = Mptr( A, AiiD, Akk+Aoff, Ald, size );
            }
            PB_Cdescset( DBUFA, K, mbb, Ainb1D, mbb, AnbD, mbb, ArocD,
                         AcurrocR, ctxt, Abufld );
         }

         if( notb )
         {
/*
*  Replicate this panel of rows or columns of sub( A ) over sub( B ) -> WA
*/
            PB_CInV( TYPE, NOCONJG, COLUMN, Bm, Bn, Bd0, mbb, Abuf, 0, 0,
                     DBUFA, &Aroc, &WA, WAd, &WAfr );
/*
*  Allocate space for temporary results in scope of sub( B ) -> WC
*/
            PB_COutV( TYPE, ROW,    INIT, Bm, Bn, Bd0, mbb, &WC, WCd, &WCfr,
                      &WCsum );
/*
*  Local matrix-matrix multiply iff I own some data
*/
            if( Bmp > 0 && Bnq > 0 )
               gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( &GemmTb ), &mbb, &Bnq, &Bmp,
                     talpha, WA, &WAd[LLD_], Bptr, &Bld, zero, WC, &WCd[LLD_] );
            if( WAfr ) free( WA );
            if( Afr && ( AisR || ( AmyprocR == AcurrocR ) ) )
               if( Abuf ) free( Abuf );
/*
*  Accumulate the intermediate results in WC
*/
            if( WCsum )
            {
               WCd[RSRC_] = Ccurrow;
               if( Bnq > 0 )
                  gsum2d( ctxt, COLUMN, &top, mbb, Bnq, WC, WCd[LLD_],
                          WCd[RSRC_], mycol );
            }
/*
*  Compute the descriptor DBUFC for the buffer that will contained the packed
*  rows of sub( C ). Allocate it.
*/
            if( ( Cfr = ( nrpq < mbb ) ) != 0 )
            {
/*
*  If rows of sub( C ) are not contiguous, then allocate the buffer
*/
               Cbufld = mbb; tbeta = zero;
               if( CisR || ( myrow == Ccurrow ) )
                  Cbuf = PB_Cmalloc( Cnq * mbb * size );
            }
            else
            {
/*
*  Otherwise re-use sub( C )
*/
               Cbufld = Cld; tbeta = BETA;
               if( CisR || ( myrow == Ccurrow ) )
                  Cbuf = Mptr( C, Ckk+Coff, Cjj, Cld, size );
            }
            PB_Cdescset( DBUFC, mbb, N, mbb, Cinb1, mbb, Cnb, Ccurrow, Ccol,
                         ctxt, Cbufld );
/*
*  Cbuf := Cbuf + WC
*/
            PB_Cpaxpby( TYPE, NOCONJG, mbb, N, one, WC, 0, 0, WCd, ROW, tbeta,
                        Cbuf, 0, 0, DBUFC, ROW );
/*
*  Unpack the mbb rows of sub( C ) and release the buffer containing them.
*/
            if( Cfr && ( CisR || ( myrow == Ccurrow ) ) )
            {
               PB_CVMpack( TYPE, &VM, ROW, ROW,    UNPACKING, NOTRAN, mbb, Cnq,
                           BETA, Mptr( C, Ckk, Cjj, Cld, size ), Cld, one, Cbuf,
                           Cbufld );
               if( Cbuf ) free( Cbuf );
            }
            if( WCfr ) free( WC );
         }
         else
         {
/*
*  Replicate this panel of rows or columns of sub( A ) over sub( B ) -> WA
*/
            PB_CInV( TYPE, NOCONJG, ROW,    Bm, Bn, Bd0, mbb, Abuf, 0, 0,
                     DBUFA, &Aroc, &WA, WAd, &WAfr );
/*
*  Allocate space for temporary results in scope of sub( A ) -> WC
*/
            PB_COutV( TYPE, COLUMN, INIT, Bm, Bn, Bd0, mbb, &WC, WCd, &WCfr,
                      &WCsum );
/*
*  Local matrix-matrix multiply iff I own some data
*/
            if( Bmp > 0 && Bnq > 0 )
               gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( &GemmTb ), &Bmp, &mbb, &Bnq,
                     talpha, Bptr, &Bld, WA, &WAd[LLD_], zero, WC, &WCd[LLD_] );
            if( WAfr ) free( WA );
            if( Afr && ( AisR || ( AmyprocR == AcurrocR ) ) )
               if( Abuf ) free( Abuf );
/*
*  Accumulate the intermediate results in WC
*/
            if( WCsum )
            {
               WCd[CSRC_] = 0;
               if( Bmp > 0 )
                  gsum2d( ctxt, ROW,    &top, Bmp, mbb, WC, WCd[LLD_], myrow,
                          WCd[CSRC_] );
            }
/*
*  Compute the descriptor DBUFC for the buffer that will contained the packed
*  rows of sub( C ). Allocate it.
*/
            if( ( Cfr = ( nrpq < mbb ) ) != 0 )
            {
/*
*  If rows of sub( C ) are not contiguous, then allocate the buffer
*/
               Cbufld = mbb; tbeta = zero;
               if( CisR || ( myrow == Ccurrow ) )
                  Cbuf = PB_Cmalloc( Cnq * mbb * size );
            }
            else
            {
/*
*  Otherwise re-use sub( C )
*/
               Cbufld = Cld; tbeta = BETA;
               if( CisR || ( myrow == Ccurrow ) )
                  Cbuf = Mptr( C, Ckk+Coff, Cjj, Cld, size );
            }
            PB_Cdescset( DBUFC, mbb, N, mbb, Cinb1, mbb, Cnb, Ccurrow, Ccol,
                         ctxt, Cbufld );
/*
*  Cbuf := Cbuf + WC'
*/
            PB_Cpaxpby( TYPE, ( TrB == CCOTRAN ? CONJG : NOCONJG ), N, mbb,
                        one, WC, 0, 0, WCd, COLUMN, tbeta, Cbuf, 0, 0, DBUFC,
                        ROW    );
/*
*  Unpack the mbb rows of sub( C ) and release the buffer containing them.
*/
            if( Cfr && ( CisR || ( myrow == Ccurrow ) ) )
            {
               PB_CVMpack( TYPE, &VM, ROW, ROW,    UNPACKING, NOTRAN, mbb, Cnq,
                           BETA, Mptr( C, Ckk, Cjj, Cld, size ), Cld, one, Cbuf,
                           Cbufld );
               if( Cbuf ) free( Cbuf );
            }
            if( WCfr ) free( WC );
         }
/*
*  Update the local indexes of sub( B ) and sub( C )
*/
         PB_CVMupdate( &VM, mbb, &Ckk, &Akk );

         npq -= mbb;
      }
/*
*  Go to next or previous virtual process row or column
*/
      if( ( Cfwd      && ( p == maxpm1 ) ) ||
          ( !( Cfwd ) && ( p == 0      ) ) )
         q = ( Afwd ? MModAdd1( q, maxq ) : MModSub1( q, maxq ) );
      p = ( Cfwd ? MModAdd1( p, maxp ) : MModSub1( p, maxp ) );
   }

   if( TrB == CCOTRAN ) free( talpha );
/*
*  End of PB_CpgemmAC
*/
}
