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
void pzhemm_( F_CHAR_T SIDE, F_CHAR_T UPLO, int * M, int * N,
              double * ALPHA,
              double * A, int * IA, int * JA, int * DESCA,
              double * B, int * IB, int * JB, int * DESCB,
              double * BETA,
              double * C, int * IC, int * JC, int * DESCC )
#else
void pzhemm_( SIDE, UPLO, M, N, ALPHA, A, IA, JA, DESCA,
              B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   F_CHAR_T       SIDE, UPLO;
   int            * IA, * IB, * IC, * JA, * JB, * JC, * M, * N;
   double         * ALPHA, * BETA;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCB, * DESCC;
   double         * A, * B, * C;
#endif
{
/*
*  Purpose
*  =======
*
*  PZHEMM  performs one of the matrix-matrix operations
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
*  Alpha  and  beta  are scalars,  sub( A ) is a Hermitian submatrix and
*  sub( B ) and sub( C ) are m by n submatrices.
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
*  SIDE    (global input) CHARACTER*1
*          On entry,  SIDE  specifies  whether the  Hermitian  submatrix
*          sub( A )  appears  on  the left or right in the operation  as
*          follows:
*
*             SIDE = 'L' or 'l'
*                   sub( C ) := alpha*sub( A )*sub( B ) + beta*sub( C ),
*
*             SIDE = 'R' or 'r'
*                   sub( C ) := alpha*sub( B )*sub( A ) + beta*sub( C ).
*
*  UPLO    (global input) CHARACTER*1
*          On  entry,   UPLO  specifies  whether  the  local  pieces  of
*          the array  A  containing the  upper or lower triangular  part
*          of the Hermitian submatrix  sub( A )  are to be referenced as
*          follows:
*
*             UPLO = 'U' or 'u'   Only the local pieces corresponding to
*                                 the   upper  triangular  part  of  the
*                                 Hermitian submatrix sub( A ) are to be
*                                 referenced,
*
*             UPLO = 'L' or 'l'   Only the local pieces corresponding to
*                                 the   lower  triangular  part  of  the
*                                 Hermitian submatrix sub( A ) are to be
*                                 referenced.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( C ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( C ). N  must be at least zero.
*
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as zero then the local entries of the arrays  A and
*          B corresponding to the entries of  the  submatrices  sub( A )
*          and sub( B ) respectively need not be set on input.
*
*  A       (local input) COMPLEX*16 array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least  Lc( 1, JA+M-1 )  when  SIDE = 'L' or 'l'  and is at
*          at least Lc( 1, JA+N-1 ) otherwise. Before  entry, this array
*          contains the local entries of the matrix A.
*          Before  entry  with  SIDE = 'L' or 'l', this  array  contains
*          the local entries corresponding to the entries of the  m by m
*          Hermitian submatrix  sub( A ), such  that  when UPLO = 'U' or
*          'u', this  array contains the local entries of the upper tri-
*          angular part of the  Hermitian submatrix  sub( A ),  and  the
*          local entries  of  the strictly lower triangular of  sub( A )
*          are not referenced, and when  UPLO = 'L' or 'l',  this  array
*          contains  the local entries of the  lower triangular part  of
*          the  Hermitian  submatrix sub( A ), and  the local entries of
*          the strictly upper triangular of sub( A ) are not referenced.
*          Before  entry  with  SIDE = 'R' or 'r', this  array  contains
*          the local entries corresponding to the entries of the  n by n
*          Hermitian submatrix  sub( A ), such  that  when UPLO = 'U' or
*          'u', this  array contains the local entries of the upper tri-
*          angular part of the  Hermitian submatrix  sub( A ),  and  the
*          local entries  of  the strictly lower triangular of  sub( A )
*          are not referenced, and when  UPLO = 'L' or 'l',  this  array
*          contains  the local entries of the  lower triangular part  of
*          the  Hermitian  submatrix sub( A ), and  the local entries of
*          the strictly upper triangular of sub( A ) are not referenced.
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
*  B       (local input) COMPLEX*16 array
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
*  BETA    (global input) COMPLEX*16
*          On entry,  BETA  specifies the scalar  beta.   When  BETA  is
*          supplied  as  zero  then  the  local entries of  the array  C
*          corresponding to  the  entries of the submatrix sub( C ) need
*          not be set on input.
*
*  C       (local input/local output) COMPLEX*16 array
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
   char           DirAB, SideOp, UploA, cbtop, cbtopsave, cctop, cctopsave,
                  rbtop, rbtopsave, rctop, rctopsave;
   int            Ai, Aj, Bi, Bj, ChooseABC, Ci, Cj, ForceTop, ctxt, info,
                  lside, mycol, myrow, nb, npcol, nprow, upper;
   double         ABCest, BCest, tmp1, tmp2, tmp3, tmp4;
   PBTYP_T        * type;
/*
*  .. Local Arrays ..
*/
   int            Ad[DLEN_], Bd[DLEN_], Cd[DLEN_];
/* ..
*  .. Executable Statements ..
*
*/
   lside = ( ( SideOp = Mupcase( F2C_CHAR( SIDE )[0] ) ) ==  CLEFT );
   upper = ( ( UploA  = Mupcase( F2C_CHAR( UPLO )[0] ) ) == CUPPER );
   PB_CargFtoC( *IA, *JA, DESCA, &Ai, &Aj, Ad );
   PB_CargFtoC( *IB, *JB, DESCB, &Bi, &Bj, Bd );
   PB_CargFtoC( *IC, *JC, DESCC, &Ci, &Cj, Cd );
#ifndef NO_ARGCHK
/*
*  Test the input parameters
*/
   Cblacs_gridinfo( ( ctxt = Ad[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
   if( !( info = ( ( nprow == -1 ) ? -( 901 + CTXT_ ) : 0 ) ) )
   {
      if( ( !lside ) && ( SideOp != CRIGHT ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PZHEMM", "Illegal SIDE = %c\n", SideOp );
         info = -1;
      }
      else if( ( !upper ) && ( UploA != CLOWER ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PZHEMM", "Illegal UPLO = %c\n", UploA );
         info = -2;
      }
      if( lside )
      {
         PB_Cchkmat( ctxt, "PZHEMM", "A", *M, 3, *M, 3, Ai, Aj, Ad,  9,
                     &info );
         PB_Cchkmat( ctxt, "PZHEMM", "B", *M, 3, *N, 4, Bi, Bj, Bd, 13,
                     &info );
      }
      else
      {
         PB_Cchkmat( ctxt, "PZHEMM", "A", *N, 4, *N, 4, Ai, Aj, Ad,  9,
                     &info );
         PB_Cchkmat( ctxt, "PZHEMM", "B", *M, 3, *N, 4, Bi, Bj, Bd, 13,
                     &info );
      }
      PB_Cchkmat(    ctxt, "PZHEMM", "C", *M, 3, *N, 4, Ci, Cj, Cd, 18,
                     &info );
   }
   if( info ) { PB_Cabort( ctxt, "PZHEMM", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( ( *M == 0 ) || ( *N == 0 ) ||
       ( ( ( ALPHA[REAL_PART] == ZERO ) && ( ALPHA[IMAG_PART] == ZERO ) ) &&
         ( ( BETA [REAL_PART] ==  ONE ) && ( BETA [IMAG_PART] == ZERO ) ) ) )
      return;
/*
*  Get type structure
*/
   type = PB_Cztypeset();
/*
*  If alpha is zero, sub( C ) := beta * sub( C ).
*/
   if( ( ALPHA[REAL_PART] == ZERO ) && ( ALPHA[IMAG_PART] == ZERO ) )
   {
      if( ( BETA[REAL_PART] == ZERO ) && ( BETA[IMAG_PART] == ZERO ) )
      {
         PB_Cplapad( type, ALL, NOCONJG, *M, *N, type->zero, type->zero,
                     ((char *) C), Ci, Cj, Cd );
      }
      else if( !( ( BETA[REAL_PART] ==  ONE ) && ( BETA[IMAG_PART] == ZERO ) ) )
      {
         PB_Cplascal( type, ALL, NOCONJG, *M, *N, ((char *) BETA), ((char *) C),
                      Ci, Cj, Cd );
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
*  ABCest: operands sub( A ), sub( B ) and sub( C ) are communicated (N >> M)
*  BCest : Both operands sub( B ) and sub( C ) are communicated      (M >> N)
*/
   if( lside )
   {
      tmp1   = DNROC( *M, Ad[MB_], nprow ); tmp2 = DNROC( *N, Bd[NB_], npcol );
      ABCest = (double)(*M) *
         ( ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp1 / TWO ) +
           ( ( ( Bd[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO :
             tmp2 + tmp2 * CBRATIO ) );
      tmp1  = DNROC( *M, Ad[MB_], nprow ); tmp2 = DNROC( *M, Ad[NB_], npcol );
      tmp3  = DNROC( *M, Bd[MB_], nprow ); tmp4 = DNROC( *M, Cd[MB_], nprow );
      BCest = (double)(*N) *
              ( CBRATIO * ( npcol == 1 ? ZERO : tmp1 ) +
                ( nprow == 1 ? ZERO : tmp2 ) + MAX( tmp2, tmp3 ) +
                ( ( ( Bd[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp1 ) +
                CBRATIO * ( nprow == 1 ? ZERO : tmp2 ) + MAX( tmp2, tmp4 ) );
   }
   else
   {
      tmp1   = DNROC( *N, Ad[NB_], npcol ); tmp2 = DNROC( *M, Bd[MB_], nprow );
      ABCest = (double)(*N) *
         ( ( ( ( Ad[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp1 / TWO ) +
           ( ( ( Bd[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO :
             tmp2 + tmp2 * CBRATIO ) );
      tmp1  = DNROC( *N, Ad[MB_], nprow ); tmp2 = DNROC( *N, Ad[NB_], npcol );
      tmp3  = DNROC( *N, Bd[NB_], npcol ); tmp4 = DNROC( *N, Cd[NB_], npcol );
      BCest = (double)(*M) *
              ( ( npcol == 1 ? ZERO : tmp1 ) + MAX( tmp1, tmp3 ) +
                CBRATIO * ( nprow == 1 ? ZERO : tmp2 ) +
                ( ( ( Bd[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp2 ) +
                CBRATIO * ( npcol == 1 ? ZERO : tmp1 ) + MAX( tmp1, tmp4 ) );
   }
/*
*  Shift a little the cross-over point between both algorithms.
*/
   ChooseABC = ( ( 1.5 * ABCest ) <= BCest );
/*
*  BLACS topologies are enforced iff M and N are strictly greater than the
*  logical block size returned by pilaenv_. Otherwise, it is assumed that the
*  routine calling this routine has already selected an adequate topology.
*/
   nb       = pilaenv_( &ctxt, C2F_CHAR( &type->type ) );
   ForceTop = ( ( *M > nb ) && ( *N > nb ) );

   rbtop = *PB_Ctop( &ctxt, BCAST,   ROW,    TOP_GET );
   rctop = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_GET );
   cbtop = *PB_Ctop( &ctxt, BCAST,   COLUMN, TOP_GET );
   cctop = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_GET );

   if( ChooseABC )
   {
      if( ForceTop )
      {
         rbtopsave = rbtop; rctopsave = rctop;
         cbtopsave = cbtop; cctopsave = cctop;

         if( lside )
         {
/*
*  No clear winner for the ring topologies, so that if a ring topology is
*  already selected, keep it.
*/
            if( ( rbtop != CTOP_DRING ) && ( rbtop != CTOP_IRING ) &&
                ( rbtop != CTOP_SRING ) )
               rbtop = *PB_Ctop( &ctxt, BCAST,   ROW,    TOP_IRING );
            if( ( ( cbtop != CTOP_DRING ) && ( cbtop != CTOP_IRING ) &&
                  ( cbtop != CTOP_SRING ) ) || ( cbtop != cctop ) )
            {
               cbtop = *PB_Ctop( &ctxt, BCAST,   COLUMN, TOP_IRING   );
               cctop = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_IRING   );
/*
*  Remove the next 2 lines when the BLACS combine operations support ring
*  topologies
*/
               rctop = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_DEFAULT );
               cctop = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_DEFAULT );
            }
         }
         else
         {
/*
*  No clear winner for the ring topologies, so that if a ring topology is
*  already selected, keep it.
*/
            if( ( cbtop != CTOP_DRING ) && ( cbtop != CTOP_IRING ) &&
                ( cbtop != CTOP_SRING ) )
               cbtop = *PB_Ctop( &ctxt, BCAST,   COLUMN, TOP_IRING );
            if( ( ( rbtop != CTOP_DRING ) && ( rbtop != CTOP_IRING ) &&
                  ( rbtop != CTOP_SRING ) ) || ( rbtop != rctop ) )
            {
               rbtop = *PB_Ctop( &ctxt, BCAST,   ROW,    TOP_IRING   );
               rctop = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_IRING   );
/*
*  Remove the next 2 lines when the BLACS combine operations support ring
*  topologies
*/
               rctop = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_DEFAULT );
               cctop = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_DEFAULT );
            }
         }
      }
      if( lside )
         DirAB = ( rbtop == CTOP_DRING ? CBACKWARD : CFORWARD );
      else
         DirAB = ( cbtop == CTOP_DRING ? CBACKWARD : CFORWARD );

      PB_CpsymmAB( type, &DirAB, CONJG, &SideOp, &UploA, *M, *N,
                   ((char *)ALPHA), ((char *)A), Ai, Aj, Ad, ((char *)B), Bi,
                   Bj, Bd, ((char *)BETA), ((char *)C), Ci, Cj, Cd );
   }
   else
   {
      if( ForceTop )
      {
         rbtopsave = rbtop; rctopsave = rctop;
         cbtopsave = cbtop; cctopsave = cctop;

         if( lside )
         {
/*
*  No clear winner for the ring topologies, so that if a ring topology is
*  already selected, keep it.
*/
            if( ( ( rbtop != CTOP_DRING ) && ( rbtop != CTOP_IRING ) &&
                  ( rbtop != CTOP_SRING ) ) || ( rbtop != rctop ) )
            {
               rbtop = *PB_Ctop( &ctxt, BCAST,   ROW,    TOP_IRING   );
               rctop = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_IRING   );
/*
*  Remove the next 2 lines when the BLACS combine operations support ring
*  topologies
*/
               rctop = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_DEFAULT );
               cctop = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_DEFAULT );
            }
            cbtop = *PB_Ctop( &ctxt, BCAST,   COLUMN, TOP_DEFAULT );
            cctop = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_DEFAULT );
         }
         else
         {
/*
*  No clear winner for the ring topologies, so that if a ring topology is
*  already selected, keep it.
*/
            if( ( ( cbtop != CTOP_DRING ) && ( cbtop != CTOP_IRING ) &&
                  ( cbtop != CTOP_SRING ) ) || ( cbtop != cctop ) )
            {
               cbtop = *PB_Ctop( &ctxt, BCAST,   COLUMN, TOP_IRING   );
               cctop = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_IRING   );
/*
*  Remove the next 2 lines when the BLACS combine operations support ring
*  topologies
*/
               rctop = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_DEFAULT );
               cctop = *PB_Ctop( &ctxt, COMBINE, COLUMN, TOP_DEFAULT );
            }
            rbtop = *PB_Ctop( &ctxt, BCAST,   ROW,    TOP_DEFAULT );
            rctop = *PB_Ctop( &ctxt, COMBINE, ROW,    TOP_DEFAULT );
         }
      }
      if( lside )
         DirAB = ( ( rbtop == CTOP_DRING || rctop == CTOP_DRING ) ?
                   CBACKWARD : CFORWARD );
      else
         DirAB = ( ( cbtop == CTOP_DRING || cctop == CTOP_DRING ) ?
                   CBACKWARD : CFORWARD );

      PB_CpsymmBC( type, &DirAB, CONJG, &SideOp, &UploA, *M, *N,
                   ((char *)ALPHA), ((char *)A), Ai, Aj, Ad, ((char *)B), Bi,
                   Bj, Bd, ((char *)BETA), ((char *)C), Ci, Cj, Cd );
   }
/*
*  Restore the BLACS topologies when necessary.
*/
   if( ForceTop )
   {
      rbtopsave = *PB_Ctop( &ctxt, BCAST,   ROW,    &rbtopsave );
      rctopsave = *PB_Ctop( &ctxt, COMBINE, ROW,    &rctopsave );
      cbtopsave = *PB_Ctop( &ctxt, BCAST,   COLUMN, &cbtopsave );
      cctopsave = *PB_Ctop( &ctxt, COMBINE, COLUMN, &cctopsave );
   }
/*
*  End of PZHEMM
*/
}
