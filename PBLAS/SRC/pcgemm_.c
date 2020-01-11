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
void pcgemm_( F_CHAR_T TRANSA, F_CHAR_T TRANSB,
              Int * M, Int * N, Int * K,
              float * ALPHA,
              float * A, Int * IA, Int * JA, Int * DESCA,
              float * B, Int * IB, Int * JB, Int * DESCB,
              float * BETA,
              float * C, Int * IC, Int * JC, Int * DESCC )
#else
void pcgemm_( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA,
              B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   F_CHAR_T       TRANSA, TRANSB;
   Int            * IA, * IB, * IC, * JA, * JB, * JC, * K, * M, * N;
   float          * ALPHA, * BETA;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA, * DESCB, * DESCC;
   float          * A, * B, * C;
#endif
{
/*
*  Purpose
*  =======
*
*  PCGEMM  performs one of the matrix-matrix operations
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
*  TRANSA  (global input) CHARACTER*1
*          On entry,  TRANSA  specifies the form of op( sub( A ) ) to be
*          used in the matrix multiplication as follows:
*
*             TRANSA = 'N' or 'n'   op( sub( A ) ) = sub( A ),
*
*             TRANSA = 'T' or 't'   op( sub( A ) ) = sub( A )',
*
*             TRANSA = 'C' or 'c'   op( sub( A ) ) = conjg( sub( A )' ).
*
*  TRANSB  (global input) CHARACTER*1
*          On entry,  TRANSB  specifies the form of op( sub( B ) ) to be
*          used in the matrix multiplication as follows:
*
*             TRANSB = 'N' or 'n'   op( sub( B ) ) = sub( B ),
*
*             TRANSB = 'T' or 't'   op( sub( B ) ) = sub( B )',
*
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
*  ALPHA   (global input) COMPLEX
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as zero then the local entries of the arrays  A and
*          B corresponding to the entries of  the  submatrices  sub( A )
*          and sub( B ) respectively need not be set on input.
*
*  A       (local input) COMPLEX array
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
*  B       (local input) COMPLEX array
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
*  BETA    (global input) COMPLEX
*          On entry,  BETA  specifies the scalar  beta.   When  BETA  is
*          supplied  as  zero  then  the  local entries of  the array  C
*          corresponding to  the  entries of the submatrix sub( C ) need
*          not be set on input.
*
*  C       (local input/local output) COMPLEX array
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
   char           DirAB, DirBC, OpC, OpR, TrA, TrB, ctop, ctopsave, rtop,
                  rtopsave;
   Int            Ai, Aj, Bi, Bj, ChooseAB, ChooseBC, Ci, Cj, ForceTop, ctxt,
                  info, mycol, myrow, nb, nota, notb, npcol, nprow;
   double         ABest, ACest, BCest, tmp1, tmp2, tmp3, tmp4;
   PBTYP_T        * type;
/*
*  .. Local Arrays ..
*/
   Int            Ad[DLEN_], Bd[DLEN_], Cd[DLEN_];
/* ..
*  .. Executable Statements ..
*
*/
   nota = ( ( TrA = Mupcase( F2C_CHAR( TRANSA )[0] ) ) == CNOTRAN );
   notb = ( ( TrB = Mupcase( F2C_CHAR( TRANSB )[0] ) ) == CNOTRAN );
   PB_CargFtoC( *IA, *JA, DESCA, &Ai, &Aj, Ad );
   PB_CargFtoC( *IB, *JB, DESCB, &Bi, &Bj, Bd );
   PB_CargFtoC( *IC, *JC, DESCC, &Ci, &Cj, Cd );
#ifndef NO_ARGCHK
/*
*  Test the input parameters
*/
   Cblacs_gridinfo( ( ctxt = Ad[CTXT_] ), &nprow, &npcol, &myrow,
                    &mycol );
   if( !( info = ( ( nprow == -1 ) ? -( 1001 + CTXT_ ) : 0 ) ) )
   {
      if( ( !nota ) && ( TrA != CTRAN ) && ( TrA != CCOTRAN ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PCGEMM", "Illegal TRANSA = %c\n", TrA );
         info = -1;
      }
      else if( ( !notb ) && ( TrB != CTRAN ) && ( TrB != CCOTRAN ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PCGEMM", "Illegal TRANSB = %c\n", TrB );
         info = -2;
      }
      if( nota )
         PB_Cchkmat( ctxt, "PCGEMM", "A", *M, 3, *K, 5, Ai, Aj, Ad, 10,
                     &info );
      else
         PB_Cchkmat( ctxt, "PCGEMM", "A", *K, 5, *M, 3, Ai, Aj, Ad, 10,
                     &info );
      if( notb )
         PB_Cchkmat( ctxt, "PCGEMM", "B", *K, 5, *N, 4, Bi, Bj, Bd, 14,
                     &info );
      else
         PB_Cchkmat( ctxt, "PCGEMM", "B", *N, 4, *K, 5, Bi, Bj, Bd, 14,
                     &info );
      PB_Cchkmat(    ctxt, "PCGEMM", "C", *M, 3, *N, 4, Ci, Cj, Cd, 19,
                     &info );
   }
   if( info ) { PB_Cabort( ctxt, "PCGEMM", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( ( *M == 0 ) || ( *N == 0 ) ||
       ( ( ( ALPHA[REAL_PART] == ZERO && ALPHA[IMAG_PART] == ZERO ) ||
           ( *K == 0                                              ) ) &&
           ( BETA [REAL_PART] == ONE  && BETA [IMAG_PART] == ZERO ) ) )
      return;
/*
*  Get type structure
*/
   type = PB_Cctypeset();
/*
*  If alpha or K is zero, sub( C ) := beta * sub( C ).
*/
   if( ( ( ALPHA[REAL_PART] == ZERO ) &&
         ( ALPHA[IMAG_PART] == ZERO ) ) || ( *K == 0 ) )
   {
      if( ( BETA[REAL_PART] == ZERO ) && ( BETA[IMAG_PART] == ZERO ) )
      {
         PB_Cplapad( type, ALL, NOCONJG, *M, *N, type->zero, type->zero,
                     ((char * ) C), Ci, Cj, Cd );
      }
      else if( !( ( BETA[REAL_PART] ==  ONE ) &&
                  ( BETA[IMAG_PART] == ZERO ) ) )
      {
         PB_Cplascal( type, ALL, NOCONJG, *M, *N, ((char *) BETA),
                      ((char * ) C), Ci, Cj, Cd );
      }
      return;
   }
/*
*  Start the operations
*/
#ifdef NO_ARGCHK
   Cblacs_gridinfo( ( ctxt = Ad[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
#endif
/*
*  Algorithm selection is based on approximation of the communication volume
*  for distributed and aligned operands.
*
*  ABest: both operands sub( A ) and sub( B ) are communicated (M, N >> K)
*  ACest: both operands sub( A ) and sub( C ) are communicated (K, N >> M)
*  BCest: both operands sub( B ) and sub( C ) are communicated (M, K >> N)
*/
   ABest = (double)(*K);
   ACest = (double)(*M);
   BCest = (double)(*N);

   if( notb )
   {
      if( nota )
      {
         tmp1 = DNROC( *M, Cd[MB_], nprow ); tmp2 = DNROC( *N, Cd[NB_], npcol );
         ABest *= ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp1 ) +
                  ( ( ( Bd[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp2 );

         tmp1 = DNROC( *K, Bd[MB_], nprow ); tmp2 = DNROC( *N, Bd[NB_], npcol );
         tmp3 = DNROC( *K, Ad[NB_], npcol );
         ACest *= ( npcol == 1 ? ZERO : tmp1 ) + MAX( tmp1, tmp3 ) +
                  CBRATIO * ( nprow == 1 ? ZERO : tmp2 );

         tmp1 = DNROC( *M, Ad[MB_], nprow ); tmp2 = DNROC( *K, Ad[NB_], npcol );
                                             tmp4 = DNROC( *K, Bd[MB_], nprow );
         BCest *= CBRATIO * ( npcol == 1 ? ZERO : tmp1 ) +
                  ( nprow == 1 ? ZERO : tmp2 ) + MAX( tmp2, tmp4 );
      }
      else
      {
         tmp1 = DNROC( *M, Cd[MB_], nprow ); tmp2 = DNROC( *N, Cd[NB_], npcol );
         tmp3 = DNROC( *M, Ad[NB_], npcol );
         ABest *= ( npcol == 1 ? ZERO : tmp1 ) + MAX( tmp1, tmp3 ) +
                  ( nprow == 1 ? ZERO : tmp2 );

         tmp1 = DNROC( *K, Bd[MB_], nprow ); tmp2 = DNROC( *N, Bd[NB_], npcol );
         ACest *= ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp1 ) +
                  CBRATIO *
                  ( ( ( Bd[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp2 );

         tmp1 = DNROC( *K, Ad[MB_], nprow ); tmp2 = DNROC( *M, Bd[NB_], npcol );
                                             tmp4 = DNROC( *M, Cd[MB_], nprow );
         BCest *= ( ( ( Bd[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp1 ) +
                  CBRATIO * ( nprow == 1 ? ZERO : tmp2 ) + MAX( tmp2, tmp4 );
      }
   }
   else
   {
      if( nota )
      {
         tmp1 = DNROC( *M, Cd[MB_], nprow ); tmp2 = DNROC( *N, Cd[NB_], npcol );
                                             tmp4 = DNROC( *N, Bd[MB_], nprow );
         ABest *= ( npcol == 1 ? ZERO : tmp1 ) +
                  ( nprow == 1 ? ZERO : tmp2 ) + MAX( tmp2, tmp4 );

         tmp1 = DNROC( *N, Bd[MB_], nprow ); tmp2 = DNROC( *K, Bd[NB_], npcol );
         tmp3 = DNROC( *N, Cd[NB_], npcol );
         ACest *= CBRATIO * ( npcol == 1 ? ZERO : tmp1 ) + MAX( tmp1, tmp3 ) +
                  ( ( ( Ad[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp2 );

         tmp1 = DNROC( *M, Ad[MB_], nprow ); tmp2 = DNROC( *K, Ad[NB_], npcol );
         BCest *= CBRATIO *
                  ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp1 ) +
                  ( ( ( Bd[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp2 );
      }
      else
      {
         tmp1 = DNROC( *M, Cd[MB_], nprow ); tmp2 = DNROC( *N, Cd[NB_], npcol );
         tmp3 = DNROC( *M, Ad[NB_], npcol ); tmp4 = DNROC( *N, Bd[MB_], nprow );
         ABest *= ( npcol == 1 ? ZERO : tmp1 ) + MAX( tmp1, tmp3 ) +
                  ( nprow == 1 ? ZERO : tmp2 ) + MAX( tmp2, tmp4 );

         tmp1 = DNROC( *N, Bd[MB_], nprow ); tmp2 = DNROC( *K, Bd[NB_], npcol );
         tmp3 = DNROC( *N, Cd[NB_], npcol ); tmp4 = DNROC( *K, Ad[MB_], nprow );
         ACest *= CBRATIO * ( npcol == 1 ? ZERO : tmp1 ) + MAX( tmp1, tmp3 ) +
                  ( nprow == 1 ? ZERO : tmp2 ) + MAX( tmp2, tmp4 );

         tmp1 = DNROC( *K, Ad[MB_], nprow ); tmp2 = DNROC( *M, Ad[NB_], npcol );
         tmp3 = DNROC( *K, Bd[NB_], npcol ); tmp4 = DNROC( *M, Cd[MB_], nprow );
         BCest *= ( npcol == 1 ? ZERO : tmp1 ) + MAX( tmp1, tmp3 ) +
                  CBRATIO * ( nprow == 1 ? ZERO : tmp2 ) + MAX( tmp2, tmp4 );
      }
   }
   ChooseAB = ( ( ABest <= ( 1.3 * BCest ) ) && ( ABest <= ( 1.3 * ACest ) ) );
   ChooseBC = ( ( BCest <= ACest           ) && ( ( 1.3 * BCest ) <= ABest ) );
/*
*  BLACS topologies are enforced iff M, N and K are strictly greater than the
*  logical block size returned by pilaenv_. Otherwise, it is assumed that the
*  routine calling this routine has already selected an adequate topology.
*/
   nb       = pilaenv_( &ctxt, C2F_CHAR( &type->type ) );
   ForceTop = ( ( *M > nb ) && ( *N > nb ) && ( *K > nb ) );

   if( ChooseAB )
   {
      OpR = CBCAST;
      OpC = CBCAST;
   }
   else if( ChooseBC )
   {
      if( nota ) { OpR = CCOMBINE; OpC = CBCAST; }
      else       { OpR = CBCAST; OpC = CCOMBINE; }
   }
   else
   {
      if( notb ) { OpR = CBCAST; OpC = CCOMBINE; }
      else       { OpR = CCOMBINE; OpC = CBCAST; }
   }

   rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_GET );
   ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_GET );

   if( ForceTop )
   {
      rtopsave = rtop;
      ctopsave = ctop;
/*
*  No clear winner for the ring topologies, so that if a ring topology is
*  already selected, keep it.
*/
      if( ( rtop != CTOP_DRING ) && ( rtop != CTOP_IRING ) &&
          ( rtop != CTOP_SRING ) )
         rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_IRING );
      if( ( ctop != CTOP_DRING ) && ( ctop != CTOP_IRING ) &&
          ( ctop != CTOP_SRING ) )
         ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_IRING );
/*
*  Remove the next 4 lines when the BLACS combine operations support ring
*  topologies
*/
      if( OpR == CCOMBINE )
         rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_DEFAULT );
      if( OpC == CCOMBINE )
         ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_DEFAULT );
   }

   DirAB = ( rtop == CTOP_DRING ? CBACKWARD : CFORWARD );
   DirBC = ( ctop == CTOP_DRING ? CBACKWARD : CFORWARD );

   if( ChooseAB )
   {
      PB_CpgemmAB( type, &DirAB, &DirBC, ( nota ? NOTRAN :
                   ( ( TrA == CCOTRAN ) ? COTRAN : TRAN ) ), ( notb ? NOTRAN :
                   ( ( TrB == CCOTRAN ) ? COTRAN : TRAN ) ), *M, *N, *K,
                   ((char *)ALPHA), ((char *)A), Ai, Aj, Ad, ((char *)B), Bi,
                   Bj, Bd, ((char *)BETA), ((char *)C), Ci, Cj, Cd );
   }
   else if( ChooseBC )
   {
      PB_CpgemmBC( type, &DirAB, &DirBC, ( nota ? NOTRAN :
                   ( ( TrA == CCOTRAN ) ? COTRAN : TRAN ) ), ( notb ? NOTRAN :
                   ( ( TrB == CCOTRAN ) ? COTRAN : TRAN ) ), *M, *N, *K,
                   ((char *)ALPHA), ((char *)A), Ai, Aj, Ad, ((char *)B), Bi,
                   Bj, Bd, ((char *)BETA), ((char *)C), Ci, Cj, Cd );
   }
   else
   {
      PB_CpgemmAC( type, &DirAB, &DirBC, ( nota ? NOTRAN :
                   ( ( TrA == CCOTRAN ) ? COTRAN : TRAN ) ), ( notb ? NOTRAN :
                   ( ( TrB == CCOTRAN ) ? COTRAN : TRAN ) ), *M, *N, *K,
                   ((char *)ALPHA), ((char *)A), Ai, Aj, Ad, ((char *)B), Bi,
                   Bj, Bd, ((char *)BETA), ((char *)C), Ci, Cj, Cd );
   }
/*
*  Restore the BLACS topologies when necessary.
*/
   if( ForceTop )
   {
      rtopsave = *PB_Ctop( &ctxt, &OpR, ROW,    &rtopsave );
      ctopsave = *PB_Ctop( &ctxt, &OpC, COLUMN, &ctopsave );
   }
/*
*  End of PCGEMM
*/
}
