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
void PB_Cpsyr2kA( PBTYP_T * TYPE, char * DIRECAB, char * CONJUG,
                  char * UPLO, char * TRANS, Int N, Int K, char * ALPHA,
                  char * A, Int IA, Int JA, Int * DESCA, char * B,
                  Int IB, Int JB, Int * DESCB, char * BETA, char * C,
                  Int IC, Int JC, Int * DESCC )
#else
void PB_Cpsyr2kA( TYPE, DIRECAB, CONJUG, UPLO, TRANS, N, K, ALPHA, A, IA,
                  JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * CONJUG, * DIRECAB, * TRANS, * UPLO;
   Int            IA, IB, IC, JA, JB, JC, K, N;
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
*  PB_Cpsyr2kA performs one of the following symmetric or Hermitian rank
*  2k operations
*
*     sub( C ) := alpha*sub( A )*sub( B )' + alpha*sub( B )*sub( A )' +
*                 beta*sub( C ),
*  or
*     sub( C ) := alpha*sub( A )*conjg( sub( B ) )' +
*                 conjg( alpha )*sub( B )*conjg( sub( A ) )' +
*                 beta*sub( C ),
*  or
*     sub( C ) := alpha*sub( A )'*sub( B ) + alpha*sub( B )'*sub( A ) +
*                 beta*sub( C ),
*  or
*     sub( C ) := alpha*conjg( sub( A )' )*sub( B ) +
*                 conjg( alpha )*conjg( sub( B )' )*sub( A ) +
*                 beta*sub( C ),
*
*  where
*
*     sub( C ) denotes C(IC:IC+N-1,JC:JC+N-1),
*
*     sub( A ) denotes A(IA:IA+N-1,JA:JA+K-1)  if TRANS = 'N',
*                      A(IA:IA+K-1,JA:JA+N-1)  otherwise, and,
*
*     sub( B ) denotes B(IB:IB+N-1,JB:JB+K-1)  if TRANS = 'N',
*                      B(IB:IB+K-1,JB:JB+N-1)  otherwise.
*
*  Alpha  and  beta  are  scalars,  sub( C )  is an n by n  symmetric or
*  Hermitian submatrix and sub( A ) and sub( B ) are n by k  submatrices
*  in the first case and k by n submatrices in the second case.
*
*  This is the outer-product algorithm  using  the  logical  LCM  hybrid
*  and static blocking technique. The submatrix operand  sub( C )  stays
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
*          or columns of sub( A ) and sub( B ) should  be looped over as
*          follows:
*             DIRECAB = 'F' or 'f'   forward  or increasing,
*             DIRECAB = 'B' or 'b'   backward or decreasing.
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
*               sub( C ) := alpha*sub( A )*sub( B )' +
*                           alpha*sub( B )*sub( A )' +
*                           beta*sub( C ),
*             or
*               sub( C ) := alpha*sub( A )*sub( B )' +
*                           alpha*sub( B )*sub( A )' +
*                           beta*sub( C ),
*             or
*               sub( C ) := alpha*sub( A )*conjg( sub( B )' ) +
*                           conjg( alpha )*sub( B )*conjg( sub( A )' ) +
*                           beta*sub( C ),
*
*             TRANS = 'T' or 't'
*               sub( C ) := alpha*sub( B )'*sub( A ) +
*                           alpha*sub( A )'*sub( B ) +
*                           beta*sub( C ),
*             or
*               sub( C ) := alpha*sub( B )'*sub( A ) +
*                           alpha*sub( A )'*sub( B ) +
*                           beta*sub( C ),
*
*             TRANS = 'C' or 'c'
*               sub( C ) := alpha*sub( B )'*sub( A ) +
*                           alpha*sub( A )'*sub( B ) +
*                           beta*sub( C ),
*             or
*               sub( C ) := alpha*conjg( sub( A )' )*sub( B ) +
*                           conjg( alpha )*conjg( sub( B )' )*sub( A ) +
*                           beta*sub( C ).
*
*  N       (global input) INTEGER
*          On entry,  N specifies the order of the  submatrix  sub( C ).
*          N must be at least zero.
*
*  K       (global input) INTEGER
*          On entry with  TRANS = 'N' or 'n',  K specifies the number of
*          columns of  the  submatrices  sub( A )  and  sub( B ), and on
*          entry with TRANS = 'T' or 't' or 'C' or 'c', K  specifies the
*          number of rows  of the  submatrices  sub( A )  and  sub( B ).
*          K  must  be at least zero.
*
*  ALPHA   (global input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of the arrays  A
*          and  B  corresponding  to  the  entries  of  the  submatrices
*          sub( A ) and sub( B ) respectively need not be set  on input.
*
*  A       (local input) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+K-1 ) when  TRANS = 'N' or 'n', and  is at
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
*  B       (local input) pointer to CHAR
*          On entry, B is an array of dimension (LLD_B, Kb), where Kb is
*          at least Lc( 1, JB+K-1 ) when  TRANS = 'N' or 'n', and  is at
*          least Lc( 1, JB+N-1 ) otherwise.  Before  entry,  this  array
*          contains the local entries of the matrix B.
*          Before entry with TRANS = 'N' or 'n', this array contains the
*          local entries corresponding to the entries of the n by k sub-
*          matrix sub( B ), otherwise the local entries corresponding to
*          the entries of the k by n submatrix sub( B ).
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
   char           * one, * talpha, * zero;
   Int            ABfwd, ABmyprocD, ABmyprocR, ABnprocsD, ABnprocsR, ABrocs,
                  Abufld, AcurrocR, Afr, AiD, AiR, AiiD, AiiR, AinbD, AinbR,
                  Ainb1D, Ainb1R, AisR, AkkR, Ald, AnbD, AnbR, AnpD, AnpR,
                  Aoff, ArocD, ArocR, AsrcR, Bbufld, BcurrocR, Bfr, BiD, BiR,
                  BiiD, BiiR, BinbD, BinbR, Binb1D, Binb1R, BisR, BkkR, Bld,
                  BnbD, BnbR, BnpD,
                  BnpR, Boff, BrocD, BrocR, BsrcR, Ccol, Cii, Cimb1,
                  Cinb1, Cjj, Clcmb, Cld, Clp, Clq, Cnq0, Cmb, Cmp, Cmp0,
                  Cnb, Cnq, Crow,
                  WACfr, WACld, WACsum, WARfr, WARld, WARsum,
                  WBCfr, WBCld, WBCsum, WBRfr, WBRld, WBRsum, Wkbb=0,
                  conjg, ctxt, k, kb, kbb, l, lb,
                  lcmb, ltmp, maxp, maxpm1, maxq, mycol, myrow, ncpq, notran,
                  npcol, npq, nprow, nrpq, p=0, q=0, size, tmp, upper;
   GEMM_T         gemm;
   TZSYR2_T       tzsyr2k;
/*
*  .. Local Arrays ..
*/
   PB_VM_T        VM;
   Int            Cd0  [DLEN_], DBUFA[DLEN_], DBUFB[DLEN_], WACd0[DLEN_],
                  WARd0[DLEN_], WBCd0[DLEN_], WBRd0[DLEN_];
   char           * Abuf = NULL, * Bbuf = NULL, * Cptr = NULL, * WAC = NULL,
                  * WAR  = NULL, * WBC  = NULL, * WBR  = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  sub( C ) := beta * sub( C )
*/
   PB_Cplascal( TYPE, UPLO, CONJUG, N, N, BETA, C, IC, JC, DESCC );
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCC[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   conjg = ( Mupcase( CONJUG [0] ) == CCONJG );

   size = TYPE->size; one = TYPE->one; zero = TYPE->zero; gemm = TYPE->Fgemm;
   kb   = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );
/*
*  Compute descriptor Cd0 for sub( C )
*/
   PB_Cdescribe( N, N, IC, JC, DESCC, nprow, npcol, myrow, mycol, &Cii, &Cjj,
                 &Cld, &Cimb1, &Cinb1, &Cmb, &Cnb, &Crow, &Ccol, Cd0 );

   Cmp = PB_Cnumroc( N, 0, Cimb1, Cmb, myrow, Crow, nprow );
   Cnq = PB_Cnumroc( N, 0, Cinb1, Cnb, mycol, Ccol, npcol );

   if( ( Cmp > 0 ) && ( Cnq > 0 ) )
   {
      Cptr = Mptr( C, Cii, Cjj, Cld, size );

      if( conjg )
      {
         talpha  = PB_Cmalloc( size ); PB_Cconjg( TYPE, ALPHA, talpha );
         tzsyr2k = PB_Ctzher2k;
      }
      else { talpha  = ALPHA; tzsyr2k = PB_Ctzsyr2k; }
/*
*  Computational partitioning size is computed as the product of the logical
*  value returned by pilaenv_ and 2 * lcm( nprow, npcol ).
*/
      Clcmb = 2 * kb * PB_Clcm( ( Crow >= 0 ? nprow : 1 ),
                                ( Ccol >= 0 ? npcol : 1 ) );
   }
/*
*  Retrieve local information for sub( A ) and sub( B )
*/
   if( ( notran = ( Mupcase( TRANS[0] ) == CNOTRAN  ) ) != 0 )
   {
      ABnprocsR = npcol;
      AiR = JA; AinbR = DESCA[INB_]; AnbR = DESCA[NB_]; AsrcR = DESCA[CSRC_];
      BiR = JB; BinbR = DESCB[INB_]; BnbR = DESCB[NB_]; BsrcR = DESCB[CSRC_];
   }
   else
   {
      ABnprocsR = nprow;
      AiR = IA; AinbR = DESCA[IMB_]; AnbR = DESCA[MB_]; AsrcR = DESCA[RSRC_];
      BiR = IB; BinbR = DESCB[IMB_]; BnbR = DESCB[MB_]; BsrcR = DESCB[RSRC_];
   }
/*
*  If sub( A ) and sub( B ) only spans one process row or column, then there is
*  no need to pack the data.
*/
   if( !( PB_Cspan( K, AiR, AinbR, AnbR, AsrcR, ABnprocsR ) ) &&
       !( PB_Cspan( K, BiR, BinbR, BnbR, BsrcR, ABnprocsR ) ) )
   {
/*
*  Replicate sub( A ) in process rows and columns spanned by sub( C ): WAR, WAC
*  Replicate sub( B ) in process rows and columns spanned by sub( C ): WBR, WBC
*/
      if( notran )
      {
         PB_CInV( TYPE, NOCONJG, COLUMN, N, N, Cd0, K, A, IA, JA, DESCA, COLUMN,
                  &WAC, WACd0, &WACfr );
         PB_CInV( TYPE, CONJUG,  ROW,    N, N, Cd0, K, WAC, 0, 0, WACd0, COLUMN,
                  &WAR, WARd0, &WARfr );

         PB_CInV( TYPE, NOCONJG, COLUMN, N, N, Cd0, K, B, IB, JB, DESCB, COLUMN,
                  &WBC, WBCd0, &WBCfr );
         PB_CInV( TYPE, CONJUG,  ROW,    N, N, Cd0, K, WBC, 0, 0, WBCd0, COLUMN,
                  &WBR, WBRd0, &WBRfr );
      }
      else
      {
         PB_CInV( TYPE, NOCONJG, ROW,    N, N, Cd0, K, A, IA, JA, DESCA, ROW,
                  &WAR, WARd0, &WARfr );
         PB_CInV( TYPE, CONJUG,  COLUMN, N, N, Cd0, K, WAR, 0, 0, WARd0, ROW,
                  &WAC, WACd0, &WACfr );

         PB_CInV( TYPE, NOCONJG, ROW,    N, N, Cd0, K, B, IB, JB, DESCB, ROW,
                  &WBR, WBRd0, &WBRfr );
         PB_CInV( TYPE, CONJUG,  COLUMN, N, N, Cd0, K, WBR, 0, 0, WBRd0, ROW,
                  &WBC, WBCd0, &WBCfr );
      }
/*
*  Perform the local update if I own some data
*/
      if( ( Cmp > 0 ) && ( Cnq > 0 ) )
      {
         WACld = WACd0[LLD_]; WBCld = WBCd0[LLD_];
         WARld = WARd0[LLD_]; WBRld = WBRd0[LLD_];

         if( Mupcase( UPLO[0] ) == CUPPER )
         {
            for( l = 0; l < N; l += Clcmb )
            {
               lb   = N - l; lb = MIN( lb, Clcmb );
               Clp  = PB_Cnumroc( l,  0, Cimb1, Cmb, myrow, Crow, nprow );
               Clq  = PB_Cnumroc( l,  0, Cinb1, Cnb, mycol, Ccol, npcol );
               Cnq0 = PB_Cnumroc( lb, l, Cinb1, Cnb, mycol, Ccol, npcol );
               if( Clp > 0 && Cnq0 > 0 )
               {
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Clp, &Cnq0, &K,
                        ALPHA,  WAC, &WACld, Mptr( WBR, 0, Clq, WBRld, size ),
                        &WBRld, one, Mptr( Cptr, 0, Clq, Cld, size ), &Cld );
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Clp, &Cnq0, &K,
                        talpha, WBC, &WBCld, Mptr( WAR, 0, Clq, WARld, size ),
                        &WARld, one, Mptr( Cptr, 0, Clq, Cld, size ), &Cld );
               }
               PB_Cpsyr2( TYPE, UPPER, lb, K, ALPHA, Mptr( WAC, Clp, 0, WACld,
                          size ), WACld, Mptr( WAR, 0, Clq, WARld, size ),
                          WARld, Mptr( WBC, Clp, 0, WBCld, size ), WBCld,
                          Mptr( WBR, 0, Clq, WBRld, size ), WBRld, Cptr, l, l,
                          Cd0, tzsyr2k );
            }
         }
         else
         {
            for( l = 0; l < N; l += Clcmb )
            {
               lb  = N - l; ltmp = l + ( lb = MIN( lb, Clcmb ) );
               Clp = PB_Cnumroc( l, 0, Cimb1, Cmb, myrow, Crow, nprow );
               Clq = PB_Cnumroc( l, 0, Cinb1, Cnb, mycol, Ccol, npcol );
               PB_Cpsyr2( TYPE, LOWER, lb, K, ALPHA, Mptr( WAC, Clp, 0, WACld,
                          size ), WACld, Mptr( WAR, 0, Clq, WARld, size ),
                          WARld, Mptr( WBC, Clp, 0, WBCld, size ), WBCld,
                          Mptr( WBR, 0, Clq, WBRld, size ), WBRld, Cptr, l, l,
                          Cd0, tzsyr2k );
               Clp  = PB_Cnumroc( ltmp, 0, Cimb1, Cmb, myrow, Crow, nprow );
               Cmp0 = Cmp - Clp;
               Cnq0 = PB_Cnumroc( lb,   l, Cinb1, Cnb, mycol, Ccol, npcol );
               if( Cmp0 > 0 && Cnq0 > 0 )
               {
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp0, &Cnq0,
                        &K, ALPHA, Mptr( WAC, Clp, 0, WACld, size ), &WACld,
                        Mptr( WBR, 0, Clq, WBRld, size ), &WBRld, one,
                        Mptr( Cptr, Clp, Clq, Cld, size ), &Cld );
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp0, &Cnq0,
                        &K, talpha, Mptr( WBC, Clp, 0, WBCld, size ), &WBCld,
                        Mptr( WAR, 0, Clq, WARld, size ), &WARld, one,
                        Mptr( Cptr, Clp, Clq, Cld, size ), &Cld );
               }
            }
         }
         if( conjg ) free( talpha );
      }

      if( WACfr ) free( WAC );
      if( WARfr ) free( WAR );
      if( WBCfr ) free( WBC );
      if( WBRfr ) free( WBR );

      return;
   }
/*
*  Otherwise sub( A ) and sub( B ) spans more than one process row or columns
*/
   ABfwd = ( Mupcase( DIRECAB[0] ) == CFORWARD );
   upper = ( Mupcase( UPLO   [0] ) ==   CUPPER );

   if( notran )
   {
      ABmyprocD = myrow; ABmyprocR = mycol; ABnprocsD = nprow;
      AiD = IA; AinbD = DESCA[IMB_]; AnbD = DESCA[MB_]; Ald = DESCA[LLD_];
      BiD = IB; BinbD = DESCB[IMB_]; BnbD = DESCB[MB_]; Bld = DESCB[LLD_];
      PB_Cinfog2l( IA, JA, DESCA, ABnprocsD, ABnprocsR, ABmyprocD, ABmyprocR,
                   &AiiD, &AiiR, &ArocD, &ArocR );
      PB_Cinfog2l( IB, JB, DESCB, ABnprocsD, ABnprocsR, ABmyprocD, ABmyprocR,
                   &BiiD, &BiiR, &BrocD, &BrocR );
   }
   else
   {
      ABmyprocD = mycol; ABmyprocR = myrow; ABnprocsD = npcol;
      AiD = JA; AinbD = DESCA[INB_]; AnbD = DESCA[NB_]; Ald = DESCA[LLD_];
      BiD = JB; BinbD = DESCB[INB_]; BnbD = DESCB[NB_]; Bld = DESCB[LLD_];
      PB_Cinfog2l( IA, JA, DESCA, ABnprocsR, ABnprocsD, ABmyprocR, ABmyprocD,
                   &AiiR, &AiiD, &ArocR, &ArocD );
      PB_Cinfog2l( IB, JB, DESCB, ABnprocsR, ABnprocsD, ABmyprocR, ABmyprocD,
                   &BiiR, &BiiD, &BrocR, &BrocD );
   }
   Ainb1D = PB_Cfirstnb( N, AiD, AinbD, AnbD );
   AnpD   = PB_Cnumroc( N, 0, Ainb1D, AnbD, ABmyprocD, ArocD, ABnprocsD );
   Ainb1R = PB_Cfirstnb( K, AiR, AinbR, AnbR );
   AisR = ( ( AsrcR < 0 ) || ( ABnprocsR == 1 ) );

   Binb1D = PB_Cfirstnb( N, BiD, BinbD, BnbD );
   BnpD   = PB_Cnumroc( N, 0, Binb1D, BnbD, ABmyprocD, BrocD, ABnprocsD );
   Binb1R = PB_Cfirstnb( K, BiR, BinbR, BnbR );
   BisR = ( ( BsrcR < 0 ) || ( ABnprocsR == 1 ) );
/*
*  When sub( A ) is not replicated and backward pass on sub( A ), find the
*  virtual process q owning the last row or column of sub( A ).
*/
   if( !( AisR ) && !( ABfwd ) )
   {
      tmp = PB_Cindxg2p( K - 1, Ainb1R, AnbR, ArocR, ArocR, ABnprocsR );
      q   = MModSub( tmp, ArocR, ABnprocsR );
   }
/*
*  When sub( B ) is not replicated and backward pass on sub( B ), find the
*  virtual process p owning the last row or column of sub( B ).
*/
   if( !( BisR ) && !( ABfwd ) )
   {
      tmp = PB_Cindxg2p( K - 1, Binb1R, BnbR, BrocR, BrocR, ABnprocsR );
      p   = MModSub( tmp, BrocR, ABnprocsR );
   }
/*
*  Allocate work space in process rows and columns spanned by sub( C )
*/
   PB_COutV( TYPE, COLUMN, NOINIT, N, N, Cd0, kb, &WAC, WACd0, &WACfr,
             &WACsum );
   PB_COutV( TYPE, ROW,    NOINIT, N, N, Cd0, kb, &WAR, WARd0, &WARfr,
             &WARsum );
   PB_COutV( TYPE, COLUMN, NOINIT, N, N, Cd0, kb, &WBC, WBCd0, &WBCfr,
             &WBCsum );
   PB_COutV( TYPE, ROW,    NOINIT, N, N, Cd0, kb, &WBR, WBRd0, &WBRfr,
             &WBRsum );
/*
*  Loop over the virtual process grid induced by the rows or columns of
*  sub( A ) and sub( B )
*/
   lcmb = PB_Clcm( ( maxp = ( BisR ? 1 : ABnprocsR ) ) * BnbR,
                   ( maxq = ( AisR ? 1 : ABnprocsR ) ) * AnbR );
   maxpm1  = maxp - 1;
/*
*  Find out process coordinates corresponding to first virtual process (p,q)
*/
   AcurrocR = ( AisR ? -1 : MModAdd( ArocR, q, ABnprocsR ) );
   AkkR     = PB_Cg2lrem( AiR, AinbR, AnbR, AcurrocR, AsrcR, ABnprocsR );
   AnpR     = PB_Cnumroc( K, 0, Ainb1R, AnbR, AcurrocR, ArocR, ABnprocsR );

   BcurrocR = ( BisR ? -1 : MModAdd( BrocR, p, ABnprocsR ) );
   BkkR     = PB_Cg2lrem( BiR, BinbR, BnbR, BcurrocR, BsrcR, ABnprocsR );
   BnpR     = PB_Cnumroc( K, 0, Binb1R, BnbR, BcurrocR, BrocR, ABnprocsR );
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
            if( ( ABfwd      && ( p == maxpm1 ) ) ||
                ( !( ABfwd ) && ( p == 0      ) ) )
               q = ( ABfwd ? MModAdd1( q, maxq ) : MModSub1( q, maxq ) );
            p = ( ABfwd ? MModAdd1( p, maxp ) : MModSub1( p, maxp ) );

            AcurrocR = ( AisR ? -1 : MModAdd( ArocR, q, ABnprocsR ) );
            AkkR     = PB_Cg2lrem(  AiR, AinbR,  AnbR, AcurrocR, AsrcR,
                                   ABnprocsR );
            AnpR     = PB_Cnumroc( K, 0, Ainb1R, AnbR, AcurrocR, ArocR,
                                   ABnprocsR );

            BcurrocR = ( BisR ? -1 : MModAdd( BrocR, p, ABnprocsR ) );
            BkkR     = PB_Cg2lrem(  BiR, BinbR,  BnbR, BcurrocR, BsrcR,
                                   ABnprocsR );
            BnpR     = PB_Cnumroc( K, 0, Binb1R, BnbR, BcurrocR, BrocR,
                                   ABnprocsR );

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

         if( notran )
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
               if( AisR || ( ABmyprocR == AcurrocR ) )
               {
                  Abuf = PB_Cmalloc( AnpD * ABrocs * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, COLUMN, PACKING, NOTRAN,
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
               if( AisR || ( ABmyprocR == AcurrocR ) )
                     Abuf = Mptr( A, AiiD, AkkR + Aoff, Ald, size );
            }
            PB_Cdescset( DBUFA, N, ABrocs, Ainb1D, ABrocs, AnbD, ABrocs,
                         ArocD, AcurrocR, ctxt, Abufld );
/*
*  Replicate panels of columns of sub( A ) over sub( C )
*/
            PB_CInV2( TYPE, NOCONJG, COLUMN, N, N, Cd0, ABrocs, Abuf, 0, 0,
                      DBUFA, COLUMN, WAC, Wkbb, WACd0 );
            if( Afr & ( AisR || ( ABmyprocR == AcurrocR ) ) )
               if( Abuf ) free( Abuf );

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
               if( BisR || ( ABmyprocR == BcurrocR ) )
               {
                  Bbuf = PB_Cmalloc( BnpD * ABrocs * size );
                  PB_CVMpack( TYPE, &VM, ROW,    COLUMN, PACKING, NOTRAN,
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
               if( BisR || ( ABmyprocR == BcurrocR ) )
                  Bbuf = Mptr( B, BiiD, BkkR + Boff, Bld, size );
            }
            PB_Cdescset( DBUFB, N, ABrocs, Binb1D, ABrocs, BnbD, ABrocs,
                         BrocD, BcurrocR, ctxt, Bbufld );
/*
*  Replicate panels of columns of sub( A ) over sub( C )
*/
            PB_CInV2( TYPE, NOCONJG, COLUMN, N, N, Cd0, ABrocs, Bbuf, 0, 0,
                      DBUFB, COLUMN, WBC, Wkbb, WBCd0 );
            if( Bfr & ( BisR || ( ABmyprocR == BcurrocR ) ) )
               if( Bbuf ) free( Bbuf );
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
               if( AisR || ( ABmyprocR == AcurrocR ) )
               {
                  Abuf = PB_Cmalloc( AnpD * ABrocs * size );
                  PB_CVMpack( TYPE, &VM, COLUMN, ROW,    PACKING, NOTRAN,
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
               if( AisR || ( ABmyprocR == AcurrocR ) )
                  Abuf = Mptr( A, AkkR + Aoff, AiiD, Ald, size );
            }
            PB_Cdescset( DBUFA, ABrocs, N, ABrocs, Ainb1D, ABrocs, AnbD,
                         AcurrocR, ArocD, ctxt, Abufld );
/*
*  Replicate panels of rows of sub( A ) over sub( C )
*/
            PB_CInV2( TYPE, NOCONJG, ROW,    N, N, Cd0, ABrocs, Abuf, 0, 0,
                      DBUFA, ROW,    WAR, Wkbb, WARd0 );
            if( Afr & ( AisR || ( ABmyprocR == AcurrocR ) ) )
               if( Abuf ) free( Abuf );
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
               if( BisR || ( ABmyprocR == BcurrocR ) )
               {
                  Bbuf = PB_Cmalloc( BnpD * ABrocs * size );
                  PB_CVMpack( TYPE, &VM, ROW,    ROW,    PACKING, NOTRAN,
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
               if( BisR || ( ABmyprocR == BcurrocR ) )
                  Bbuf = Mptr( B, BkkR + Boff, BiiD, Bld, size );
            }
            PB_Cdescset( DBUFB, ABrocs, N, ABrocs, Binb1D, ABrocs, BnbD,
                         BcurrocR, BrocD, ctxt, Bbufld );
/*
*  Replicate panels of rows of sub( B ) over sub( C )
*/
            PB_CInV2( TYPE, NOCONJG, ROW,    N, N, Cd0, ABrocs, Bbuf, 0, 0,
                      DBUFB, ROW,    WBR, Wkbb, WBRd0 );
            if( Bfr & ( BisR || ( ABmyprocR == BcurrocR ) ) )
               if( Bbuf ) free( Bbuf );
         }
/*
*  Update the local indexes of sub( A ) and sub( B )
*/
         PB_CVMupdate( &VM, ABrocs, &BkkR, &AkkR );
/*
*  ABrocs rows or columns of sub( A ) and sub( B ) have been replicated,
*  update the number of diagonals in this virtual process as well as the
*  number of rows or columns of sub( A ) and sub( B ) that are in WA, WB.
*/
         npq  -= ABrocs;
         Wkbb += ABrocs;
      }

      if( notran )
      {
/*
*  WAR := WAC'
*/
         PB_CInV2( TYPE, CONJUG,  ROW,    N, N, Cd0, kbb, WAC,  0, 0, WACd0,
                   COLUMN, WAR, 0, WARd0 );
/*
*  WBR := WBC'
*/
         PB_CInV2( TYPE, CONJUG,  ROW,    N, N, Cd0, kbb, WBC,  0, 0, WBCd0,
                   COLUMN, WBR, 0, WBRd0 );
      }
      else
      {
/*
*  WAC := WAR'
*/
         PB_CInV2( TYPE, CONJUG,  COLUMN, N, N, Cd0, kbb, WAR,  0, 0, WARd0,
                   ROW,    WAC, 0, WACd0 );
/*
*  WBC := WBR'
*/
         PB_CInV2( TYPE, CONJUG,  COLUMN, N, N, Cd0, kbb, WBR,  0, 0, WBRd0,
                   ROW,    WBC, 0, WBCd0 );
      }
/*
*  Perform the local update if I own some data
*/
      if( ( Cmp > 0 ) && ( Cnq > 0 ) )
      {
         WACld = WACd0[LLD_]; WBCld = WBCd0[LLD_];
         WARld = WARd0[LLD_]; WBRld = WBRd0[LLD_];

         if( upper )
         {
            for( l = 0; l < N; l += Clcmb )
            {
               lb   = N - l; lb = MIN( lb, Clcmb );
               Clp  = PB_Cnumroc( l,  0, Cimb1, Cmb, myrow, Crow, nprow );
               Clq  = PB_Cnumroc( l,  0, Cinb1, Cnb, mycol, Ccol, npcol );
               Cnq0 = PB_Cnumroc( lb, l, Cinb1, Cnb, mycol, Ccol, npcol );
               if( Clp > 0 && Cnq0 > 0 )
               {
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Clp, &Cnq0,
                        &kbb, ALPHA,  WAC, &WACld, Mptr( WBR, 0, Clq, WBRld,
                        size ), &WBRld, one, Mptr( Cptr, 0, Clq, Cld, size ),
                        &Cld );
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Clp, &Cnq0,
                        &kbb, talpha, WBC, &WBCld, Mptr( WAR, 0, Clq, WARld,
                        size ), &WARld, one, Mptr( Cptr, 0, Clq, Cld, size ),
                        &Cld );
               }
               PB_Cpsyr2( TYPE, UPPER, lb, kbb, ALPHA, Mptr( WAC, Clp, 0, WACld,
                          size ), WACld, Mptr( WAR, 0, Clq, WARld, size ),
                          WARld, Mptr( WBC, Clp, 0, WBCld, size ), WBCld,
                          Mptr( WBR, 0, Clq, WBRld, size ), WBRld, Cptr, l, l,
                          Cd0, tzsyr2k );
            }
         }
         else
         {
            for( l = 0; l < N; l += Clcmb )
            {
               lb  = N - l; ltmp = l + ( lb = MIN( lb, Clcmb ) );
               Clp = PB_Cnumroc( l, 0, Cimb1, Cmb, myrow, Crow, nprow );
               Clq = PB_Cnumroc( l, 0, Cinb1, Cnb, mycol, Ccol, npcol );

               PB_Cpsyr2( TYPE, LOWER, lb, kbb, ALPHA, Mptr( WAC, Clp, 0, WACld,
                          size ), WACld, Mptr( WAR, 0, Clq, WARld, size ),
                          WARld, Mptr( WBC, Clp, 0, WBCld, size ), WBCld,
                          Mptr( WBR, 0, Clq, WBRld, size ), WBRld, Cptr, l, l,
                          Cd0, tzsyr2k );
               Clp  = PB_Cnumroc( ltmp, 0, Cimb1, Cmb, myrow, Crow, nprow );
               Cmp0 = Cmp - Clp;
               Cnq0 = PB_Cnumroc( lb,   l, Cinb1, Cnb, mycol, Ccol, npcol );
               if( Cmp0 > 0 && Cnq0 > 0 )
               {
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp0, &Cnq0,
                        &kbb, ALPHA, Mptr( WAC, Clp, 0, WACld, size ), &WACld,
                        Mptr( WBR, 0, Clq, WBRld, size ), &WBRld, one,
                        Mptr( Cptr, Clp, Clq, Cld, size ), &Cld );
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp0, &Cnq0,
                        &kbb, talpha, Mptr( WBC, Clp, 0, WBCld, size ), &WBCld,
                        Mptr( WAR, 0, Clq, WARld, size ), &WARld, one,
                        Mptr( Cptr, Clp, Clq, Cld, size ), &Cld );
               }
            }
         }
      }

      Wkbb = 0;
   }

   if( WACfr ) free( WAC );
   if( WARfr ) free( WAR );
   if( WBCfr ) free( WBC );
   if( WBRfr ) free( WBR );

   if( conjg && ( Cmp > 0 ) && ( Cnq > 0 ) ) free( talpha );
/*
*  End of PB_Cpsyr2kA
*/
}
