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
void PB_CpsymmBC( PBTYP_T * TYPE, char * DIRECAB, char * CONJUG,
                  char * SIDE, char * UPLO, int M, int N, char * ALPHA,
                  char * A, int IA, int JA, int * DESCA, char * B,
                  int IB, int JB, int * DESCB, char * BETA, char * C,
                  int IC, int JC, int * DESCC )
#else
void PB_CpsymmBC( TYPE, DIRECAB, CONJUG, SIDE, UPLO, M, N, ALPHA, A, IA,
                  JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * CONJUG, * DIRECAB, * SIDE, * UPLO;
   int            IA, IB, IC, JA, JB, JC, M, N;
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
*  PB_CpsymmBC  performs one of the matrix-matrix operations
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
*  DIRECAB (global input) pointer to CHAR
*          On entry, DIRECAB specifies  the direction in which the  rows
*          or columns of sub( B ) should be looped over as follows:
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
   char           GemmTa, GemmTb, cctop, * one, rctop, * talphaCR, * talphaRC,
                  * tbeta, * zero;
   int            Acol, Aii, Aimb1, Ainb1, Ajj, Alcmb, Ald, Alp, Alp0, Alq,
                  Alq0, Amb, Amp, An, Anb, Anq, Arow, BCfwd, BCmyprocD,
                  BCmyprocR, BCnD, BCnR, BCnprocsD, BCnprocsR, Bbufld, BcurrocR,
                  Bfr, BiD, BiR, BiiD, BiiR, BinbD, BinbR, Binb1D, Binb1R, BisR,
                  Bkk, Bld, BnbD, BnbR, BnpD, BnpR, Boff, BrocD, BrocR, BsrcR,
                  Cfr, CiD, CiR, CiiD, CiiR, CinbD, CinbR, Cinb1D, Cinb1R, Ckk,
                  CnbD, CnbR, CnpD, CnpR, Coff, CrocD, CrocR, CsrcR, Cbufld,
                  CcurrocR, CisR, Cld, WBCfr, WBCld, WBRfr, WBRld, WCCfr, WCCld,
                  WCCsum, WCRfr, WCRld, WCRsum, conjg, ctxt, l, lb, lcmb, lside,
                  ltmp, maxp, maxpm1, maxq, mycol, myrow, n, nb, nbb, ncpq,
                  npcol, npq=0, nprow, nrpq, p=0, q=0, size, tmp, upper;
   TZSYM_T        tzsymm;
   GEMM_T         gemm;
   GSUM2D_T       gsum2d;
/*
*  .. Local Arrays ..
*/
   PB_VM_T        VM;
   int            Ad0 [DLEN_], DBUFB[DLEN_], DBUFC[DLEN_], WBCd[DLEN_],
                  WBRd[DLEN_], WCCd [DLEN_], WCRd [DLEN_];
   char           * Aptr = NULL, * Bbuf = NULL, * Cbuf = NULL, * WBC  = NULL,
                  * WBR  = NULL, * WCC  = NULL, * WCR  = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   Cblacs_gridinfo( ( ctxt = DESCC[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   BCfwd = ( Mupcase( DIRECAB[0] ) == CFORWARD );
   conjg = ( Mupcase( CONJUG [0] ) ==   CCONJG );
   lside = ( Mupcase( SIDE   [0] ) ==    CLEFT );
   upper = ( Mupcase( UPLO   [0] ) ==   CUPPER );

   size  = TYPE->size;  one    = TYPE->one;     zero  = TYPE->zero;
   gemm  = TYPE->Fgemm; gsum2d = TYPE->Cgsum2d;
   nb    = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );
/*
*  Compute local information for sub( A ), sub( B ) and sub( C )
*/
   if( lside )
   {
      BCnD = An = M;            BCnR      = N;
      BCmyprocD = myrow;        BCnprocsD = nprow;
      BCmyprocR = mycol;        BCnprocsR = npcol;

      BiD       = IB;           BiR       = JB;
      BinbD     = DESCB[IMB_ ]; BinbR     = DESCB[INB_];
      BnbD      = DESCB[MB_  ]; BnbR      = DESCB[NB_ ];
      BsrcR     = DESCB[CSRC_]; Bld       = DESCB[LLD_];
      PB_Cinfog2l( IB, JB, DESCB, BCnprocsD, BCnprocsR, BCmyprocD, BCmyprocR,
                   &BiiD, &BiiR, &BrocD, &BrocR );

      CiD       = IC;           CiR       = JC;
      CinbD     = DESCC[IMB_ ]; CinbR     = DESCC[INB_];
      CnbD      = DESCC[MB_  ]; CnbR      = DESCC[NB_ ];
      CsrcR     = DESCC[CSRC_]; Cld       = DESCC[LLD_];
      PB_Cinfog2l( IC, JC, DESCC, BCnprocsD, BCnprocsR, BCmyprocD, BCmyprocR,
                   &CiiD, &CiiR, &CrocD, &CrocR );
   }
   else
   {
      BCnD = An = N;            BCnR      = M;
      BCmyprocD = mycol;        BCnprocsD = npcol;
      BCmyprocR = myrow;        BCnprocsR = nprow;

      BiD       = JB;           BiR       = IB;
      BinbR     = DESCB[IMB_ ]; BinbD     = DESCB[INB_];
      BnbR      = DESCB[MB_  ]; BnbD      = DESCB[NB_ ];
      BsrcR     = DESCB[RSRC_]; Bld       = DESCB[LLD_];
      PB_Cinfog2l( IB, JB, DESCB, BCnprocsR, BCnprocsD, BCmyprocR, BCmyprocD,
                   &BiiR, &BiiD, &BrocR, &BrocD );

      CiD       = JC;           CiR       = IC;
      CinbR     = DESCC[IMB_ ]; CinbD     = DESCC[INB_];
      CnbR      = DESCC[MB_  ]; CnbD      = DESCC[NB_ ];
      CsrcR     = DESCC[RSRC_]; Cld       = DESCC[LLD_];
      PB_Cinfog2l( IC, JC, DESCC, BCnprocsR, BCnprocsD, BCmyprocR, BCmyprocD,
                   &CiiR, &CiiD, &CrocR, &CrocD );
   }

   Binb1D = PB_Cfirstnb( BCnD, BiD, BinbD, BnbD );
   BnpD   = PB_Cnumroc( BCnD, 0, Binb1D, BnbD, BCmyprocD, BrocD, BCnprocsD );
   Binb1R = PB_Cfirstnb( BCnR, BiR, BinbR, BnbR );
   BisR   = ( ( BsrcR < 0 ) || ( BCnprocsR == 1 ) );

   Cinb1D = PB_Cfirstnb( BCnD, CiD, CinbD, CnbD );
   CnpD   = PB_Cnumroc( BCnD, 0, Cinb1D, CnbD, BCmyprocD, CrocD, BCnprocsD );
   Cinb1R = PB_Cfirstnb( BCnR, CiR, CinbR, CnbR );
   CisR   = ( ( CsrcR < 0 ) || ( BCnprocsR == 1 ) );
/*
*  Compute descriptor Ad0 for sub( A )
*/
   PB_Cdescribe( An, An, IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aii, &Ajj,
                 &Ald, &Aimb1, &Ainb1, &Amb, &Anb, &Arow, &Acol, Ad0 );

   Amp = PB_Cnumroc( An, 0, Aimb1, Amb, myrow, Arow, nprow );
   Anq = PB_Cnumroc( An, 0, Ainb1, Anb, mycol, Acol, npcol );
   if( ( Amp > 0 ) && ( Anq > 0 ) ) Aptr = Mptr( A, Aii, Ajj, Ald, size );
/*
*  Retrieve the BLACS combine topologies, compute conjugate of alpha for the
*  Hermitian case and set the transpose parameters to be passed to the BLAS
*  matrix multiply routine.
*/
   rctop = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_GET );
   cctop = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );

   if( conjg )
   {
      tzsymm = PB_Ctzhemm;
      if( lside )
      {
         talphaRC = ALPHA; GemmTa = CCOTRAN; GemmTb = CTRAN;
         talphaCR = PB_Cmalloc( size );
         PB_Cconjg( TYPE, ALPHA, talphaCR );
      }
      else
      {
         talphaCR = ALPHA; GemmTa = CTRAN; GemmTb = CCOTRAN;
         talphaRC = PB_Cmalloc( size );
         PB_Cconjg( TYPE, ALPHA, talphaRC );
      }
   }
   else
   {
      tzsymm = PB_Ctzsymm; talphaCR = talphaRC = ALPHA;
      GemmTa = CTRAN; GemmTb = CTRAN;
   }
/*
*  Computational partitioning size is computed as the product of the logical
*  value returned by pilaenv_ and 2 * lcm( nprow, npcol ).
*/
   Alcmb  = 2 * nb * PB_Clcm( ( Arow >= 0 ? nprow : 1 ),
                              ( Acol >= 0 ? npcol : 1 ) );
/*
*  When sub( B ) is not replicated and backward pass on sub( B ), find the
*  virtual process q owning the last row or column of sub( B ).
*/
   if( !( BisR ) && !( BCfwd ) )
   {
      tmp = PB_Cindxg2p( BCnR - 1, Binb1R, BnbR, BrocR, BrocR, BCnprocsR );
      q   = MModSub( tmp, BrocR, BCnprocsR );
   }
/*
*  When sub( C ) is not replicated and backward pass on sub( C ), find the
*  virtual process p owning the last row or column of sub( C ).
*/
   if( !( CisR ) && !( BCfwd ) )
   {
      tmp = PB_Cindxg2p( BCnR - 1, Cinb1R, CnbR, CrocR, CrocR, BCnprocsR );
      p   = MModSub( tmp, CrocR, BCnprocsR );
   }
/*
*  Loop over the virtual process grid induced by the rows or columns of
*  sub( B ) and sub( C ).
*/
   lcmb   = PB_Clcm( ( maxp = ( CisR ? 1 : BCnprocsR ) ) * CnbR,
                     ( maxq = ( BisR ? 1 : BCnprocsR ) ) * BnbR );
   n      = BCnR;
   maxpm1 = maxp - 1;

   while( n > 0 )
   {
/*
*  Initialize local virtual matrix in process (p,q)
*/
      BcurrocR = ( BisR ? -1 : MModAdd( BrocR, q, BCnprocsR ) );
      Bkk  = PB_Cg2lrem( BiR, BinbR, BnbR, BcurrocR, BsrcR, BCnprocsR );
      BnpR = PB_Cnumroc( BCnR, 0, Binb1R, BnbR, BcurrocR, BrocR, BCnprocsR );

      CcurrocR = ( CisR ? -1 : MModAdd( CrocR, p, BCnprocsR ) );
      Ckk  = PB_Cg2lrem( CiR, CinbR, CnbR, CcurrocR, CsrcR, BCnprocsR );
      CnpR = PB_Cnumroc( BCnR, 0, Cinb1R, CnbR, CcurrocR, CrocR, BCnprocsR );

      PB_CVMinit( &VM, 0, CnpR, BnpR, Cinb1R, Binb1R, CnbR, BnbR, p, q,
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

         if( lside )
         {
/*
*  Compute the descriptor DBUFB for the buffer that will contained the packed
*  columns of sub( B ).
*/
            if( ( Bfr = ( ncpq < nbb ) ) != 0 )
            {
/*
*  If columns of sub( B ) are not contiguous, then allocate the buffer and
*  pack the kbb columns of sub( B ).
*/
               Bbufld = MAX( 1, BnpD );
               if( BisR || ( BCmyprocR == BcurrocR ) )
               {
                  Bbuf = PB_Cmalloc( BnpD * nbb * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, COLUMN, PACKING, NOTRAN, nbb,
                              BnpD, one, Mptr( B, BiiD, Bkk, Bld, size ), Bld,
                              zero, Bbuf, Bbufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( B ) directly.
*/
               Bbufld = Bld;
               if( BisR || ( BCmyprocR == BcurrocR ) )
                  Bbuf = Mptr( B, BiiD, Bkk+Boff, Bld, size );
            }
            PB_Cdescset( DBUFB, BCnD, nbb, Binb1D, nbb, BnbD, nbb, BrocD,
                         BcurrocR, ctxt, Bbufld );
/*
*  Replicate this panel of columns of sub( B ) as well as its transposed
*  over sub( A ) -> WBC, WBR
*/
            PB_CInV( TYPE, NOCONJG, COLUMN, An, An, Ad0, nbb, Bbuf, 0, 0,
                     DBUFB, COLUMN, &WBC, WBCd, &WBCfr );
            PB_CInV( TYPE, NOCONJG, ROW,    An, An, Ad0, nbb, WBC,  0, 0,
                     WBCd,  COLUMN, &WBR, WBRd, &WBRfr );
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
*  the kbb rows of sub( B ).
*/
               Bbufld = nbb;
               if( BisR || ( BCmyprocR == BcurrocR ) )
               {
                  Bbuf = PB_Cmalloc( BnpD * nbb * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, ROW,    PACKING, NOTRAN, nbb,
                              BnpD, one, Mptr( B, Bkk, BiiD, Bld, size ), Bld,
                              zero, Bbuf, Bbufld );
               }
            }
            else
            {
/*
*  Otherwise, re-use sub( B ) directly.
*/
               Bbufld = Bld;
               if( BisR || ( BCmyprocR == BcurrocR ) )
                  Bbuf = Mptr( B, Bkk+Boff, BiiD, Bld, size );
            }
            PB_Cdescset( DBUFB, nbb, BCnD, nbb, Binb1D, nbb, BnbD, BcurrocR,
                         BrocD, ctxt, Bbufld );
/*
*  Replicate this panel of rows of sub( B ) as well as its transposed
*  over sub( A ) -> WBR, WBC
*/
            PB_CInV( TYPE, NOCONJG, ROW,    An, An, Ad0, nbb, Bbuf, 0, 0,
                     DBUFB, ROW,    &WBR, WBRd, &WBRfr );
            PB_CInV( TYPE, NOCONJG, COLUMN, An, An, Ad0, nbb, WBR,  0, 0,
                     WBRd,  ROW,    &WBC, WBCd, &WBCfr );
         }
/*
*  Allocate space for temporary results in scope of sub( A ) -> WCC, WCR
*/
         PB_COutV( TYPE, COLUMN, INIT, An, An, Ad0, nbb, &WCC, WCCd, &WCCfr,
                   &WCCsum );
         PB_COutV( TYPE, ROW,    INIT, An, An, Ad0, nbb, &WCR, WCRd, &WCRfr,
                   &WCRsum );
/*
*  Local matrix-matrix multiply iff I own some data
*/
         WBCld = WBCd[LLD_]; WBRld = WBRd[LLD_];
         WCCld = WCCd[LLD_]; WCRld = WCRd[LLD_];

         if( ( Amp > 0 ) && ( Anq > 0 ) )
         {
            if( upper )
            {
/*
*  sub( A ) is upper triangular
*/
               for( l = 0; l < An; l += Alcmb )
               {
                  lb   = An - l; lb = MIN( lb, Alcmb );
                  Alp  = PB_Cnumroc( l,  0, Aimb1, Amb, myrow, Arow, nprow );
                  Alq  = PB_Cnumroc( l,  0, Ainb1, Anb, mycol, Acol, npcol );
                  Alq0 = PB_Cnumroc( lb, l, Ainb1, Anb, mycol, Acol, npcol );
                  if( Alp > 0 && Alq0 > 0 )
                  {
                     gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( &GemmTb ), &Alp, &nbb,
                           &Alq0, talphaRC, Mptr( Aptr, 0, Alq, Ald, size ),
                           &Ald, Mptr( WBR, 0, Alq, WBRld, size ), &WBRld, one,
                           WCC, &WCCld );
                     gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( NOTRAN ), &nbb, &Alq0,
                           &Alp, talphaCR, WBC, &WBCld, Mptr( Aptr, 0, Alq, Ald,
                           size ), &Ald, one, Mptr( WCR, 0, Alq, WCRld, size ),
                           &WCRld );
                  }
                  PB_Cpsym( TYPE, TYPE, SIDE, UPPER, lb, nbb, ALPHA, Aptr, l, l,
                            Ad0, Mptr( WBC, Alp, 0, WBCld, size ), WBCld,
                            Mptr( WBR, 0, Alq, WBRld, size ), WBRld, Mptr( WCC,
                            Alp, 0, WCCld, size ), WCCld, Mptr( WCR, 0, Alq,
                            WCRld, size ), WCRld, tzsymm );
               }
            }
            else
            {
/*
*  sub( A ) is lower triangular
*/
               for( l = 0; l < An; l += Alcmb )
               {
                  lb  = An - l; ltmp = l + ( lb = MIN( lb, Alcmb ) );
                  Alp = PB_Cnumroc( l, 0, Aimb1, Amb, myrow, Arow, nprow );
                  Alq = PB_Cnumroc( l, 0, Ainb1, Anb, mycol, Acol, npcol );
                  PB_Cpsym( TYPE, TYPE, SIDE, LOWER, lb, nbb, ALPHA, Aptr, l, l,
                            Ad0, Mptr( WBC, Alp, 0, WBCld, size ), WBCld,
                            Mptr( WBR, 0, Alq, WBRld, size ), WBRld, Mptr( WCC,
                            Alp, 0, WCCld, size ), WCCld, Mptr( WCR, 0, Alq,
                            WCRld, size ), WCRld, tzsymm );
                  Alp  = PB_Cnumroc( ltmp, 0, Aimb1, Amb, myrow, Arow, nprow );
                  Alp0 = Amp - Alp;
                  Alq0 = PB_Cnumroc( lb,   l, Ainb1, Anb, mycol, Acol, npcol );
                  if( Alp0 > 0 && Alq0 > 0 )
                  {
                     gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( &GemmTb ), &Alp0, &nbb,
                           &Alq0, talphaRC, Mptr( Aptr, Alp, Alq, Ald, size ),
                           &Ald, Mptr( WBR, 0, Alq, WBRld, size ), &WBRld, one,
                           Mptr( WCC, Alp, 0, WCCld, size ), &WCCld );
                     gemm( C2F_CHAR( &GemmTa ), C2F_CHAR( NOTRAN ), &nbb, &Alq0,
                           &Alp0, talphaCR, Mptr( WBC, Alp, 0, WBCld, size ),
                           &WBCld, Mptr( Aptr, Alp, Alq, Ald, size ), &Ald, one,
                           Mptr( WCR, 0, Alq, WCRld, size ), &WCRld );
                  }
               }
            }
         }
         if( WBCfr ) free( WBC );
         if( WBRfr ) free( WBR );

         if( Bfr && ( BisR || ( BCmyprocR == BcurrocR ) ) )
            if( Bbuf ) free( Bbuf );

         if( lside )
         {
/*
*  Accumulate the intermediate results in WCC and WCR
*/
            if( WCCsum )
            {
               WCCd[CSRC_] = CcurrocR;
               if( Amp > 0 )
                  gsum2d( ctxt, ROW,    &rctop, Amp, nbb, WCC, WCCld, myrow,
                          WCCd[CSRC_] );
            }

            if( WCRsum )
            {
               WCRd[RSRC_] = 0;
               if( Anq > 0 )
                  gsum2d( ctxt, COLUMN, &cctop, nbb, Anq, WCR, WCRld,
                          WCRd[RSRC_], mycol );
            }
/*
*  WCC := WCC + WCR'
*/
            PB_Cpaxpby( TYPE, CONJUG, nbb, An, one, WCR, 0, 0, WCRd, ROW, one,
                        WCC, 0, 0, WCCd, COLUMN );
            if( WCRfr ) free( WCR );
/*
*  Compute the descriptor DBUFC for the buffer that will contained the packed
*  columns of sub( C ). Allocate it.
*/
            if( ( Cfr = ( nrpq < nbb ) ) != 0 )
            {
/*
*  If columns of sub( C ) are not contiguous, then allocate the buffer
*/
               Cbufld = MAX( 1, CnpD );  tbeta = zero;
               if( CisR || ( BCmyprocR == CcurrocR ) )
                  Cbuf = PB_Cmalloc( CnpD * nbb * size );
            }
            else
            {
/*
*  Otherwise re-use sub( C )
*/
               Cbufld = Cld;             tbeta = BETA;
               if( CisR || ( BCmyprocR == CcurrocR ) )
                  Cbuf = Mptr( C, CiiD, Ckk+Coff, Cld, size );
            }
            PB_Cdescset( DBUFC, BCnD, nbb, Cinb1D, nbb, CnbD, nbb, CrocD,
                         CcurrocR, ctxt, Cbufld );
/*
*  sub( C ) := beta * sub( C ) + WCC
*/
            PB_Cpaxpby( TYPE, NOCONJG, An, nbb, one, WCC, 0, 0, WCCd, COLUMN,
                        tbeta, Cbuf, 0, 0, DBUFC, COLUMN );
            if( WCCfr ) free( WCC );
/*
*  Unpack the kbb columns of sub( C ) and release the buffer containing them.
*/
            if( Cfr && ( CisR || ( BCmyprocR == CcurrocR ) ) )
            {
               PB_CVMpack( TYPE, &VM, ROW,   COLUMN, UNPACKING, NOTRAN, nbb,
                           CnpD, BETA, Mptr( C, CiiD, Ckk, Cld, size ), Cld,
                           one, Cbuf, Cbufld );
               if( Cbuf ) free( Cbuf );
            }
         }
         else
         {
/*
*  Accumulate the intermediate results in WCC and WCR
*/
            if( WCCsum )
            {
               WCCd[CSRC_] = 0;
               if( Amp > 0 )
                  gsum2d( ctxt, ROW,    &rctop, Amp, nbb, WCC, WCCld, myrow,
                          WCCd[CSRC_] );
            }

            if( WCRsum )
            {
               WCRd[RSRC_] = CcurrocR;
               if( Anq > 0 )
                  gsum2d( ctxt, COLUMN, &cctop, nbb, Anq, WCR, WCRld,
                          WCRd[RSRC_], mycol );
            }
/*
*  WCR := WCR + WCC'
*/
            PB_Cpaxpby( TYPE, CONJUG, An, nbb, one, WCC, 0, 0, WCCd, COLUMN,
                        one, WCR, 0, 0, WCRd, ROW );
            if( WCCfr ) free( WCC );
/*
*  Compute the descriptor DBUFC for the buffer that will contained the packed
*  rows of sub( C ). Allocate it.
*/
            if( ( Cfr = ( nrpq < nbb ) ) != 0 )
            {
/*
*  If rows of sub( C ) are not contiguous, then allocate receiving buffer.
*/
               Cbufld = nbb; tbeta = zero;
               if( CisR || ( BCmyprocR == CcurrocR ) )
                  Cbuf = PB_Cmalloc( CnpD * nbb * size );
            }
            else
            {
/*
*  Otherwise re-use sub( C )
*/
               Cbufld = Cld; tbeta = BETA;
               if( CisR || ( BCmyprocR == CcurrocR ) )
                  Cbuf = Mptr( C, Ckk+Coff, CiiD, Cld, size );
            }
            PB_Cdescset( DBUFC, nbb, BCnD, nbb, Cinb1D, nbb, CnbD, CcurrocR,
                         CrocD, ctxt, Cbufld );
/*
*  sub( C ) := beta * sub( C ) + WCR
*/
            PB_Cpaxpby( TYPE, NOCONJG, nbb, An, one, WCR, 0, 0, WCRd, ROW,
                        tbeta, Cbuf, 0, 0, DBUFC, ROW );

            if( WCRfr ) free( WCR );
/*
*  Unpack the kbb rows of sub( C ) and release the buffer containing them.
*/
            if( Cfr && ( CisR || ( BCmyprocR == CcurrocR ) ) )
            {
               PB_CVMpack( TYPE, &VM, ROW,   ROW,    UNPACKING, NOTRAN, nbb,
                           CnpD, BETA, Mptr( C, Ckk, CiiD, Cld, size ), Cld,
                           one, Cbuf, Cbufld );
               if( Cbuf ) free( Cbuf );
            }
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
      if( ( BCfwd      && ( p == maxpm1 ) ) ||
          ( !( BCfwd ) && ( p == 0      ) ) )
         q = ( BCfwd ? MModAdd1( q, maxq ) : MModSub1( q, maxq ) );
      p = ( BCfwd ? MModAdd1( p, maxp ) : MModSub1( p, maxp ) );
   }

   if( conjg ) free( ( lside ? talphaCR : talphaRC ) );
/*
*  End of PB_CpsymmBC
*/
}
