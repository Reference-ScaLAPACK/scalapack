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
void pctrsm_( F_CHAR_T SIDE, F_CHAR_T UPLO, F_CHAR_T TRANS, F_CHAR_T DIAG,
              int * M, int * N, float * ALPHA,
              float * A, int * IA, int * JA, int * DESCA,
              float * B, int * IB, int * JB, int * DESCB )
#else
void pctrsm_( SIDE, UPLO, TRANS, DIAG, M, N, ALPHA,
              A, IA, JA, DESCA, B, IB, JB, DESCB )
/*
*  .. Scalar Arguments ..
*/
   F_CHAR_T       DIAG, SIDE, TRANS, UPLO;
   int            * IA, * IB, * JA, * JB, * M, * N;
   float          * ALPHA;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCB;
   float          * A, * B;
#endif
{
/*
*  Purpose
*  =======
*
*  PCTRSM  solves one of the matrix equations
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
*          On entry,  SIDE  specifies  whether op( sub( A ) ) appears on
*          the left or right of X as follows:
*
*             SIDE = 'L' or 'l'   op( sub( A ) )*X = alpha*sub( B ),
*
*             SIDE = 'R' or 'r'   X*op( sub( A ) ) = alpha*sub( B ).
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
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( B ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( B ). N  must be at least zero.
*
*  ALPHA   (global input) COMPLEX
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of  the array  B
*          corresponding to the entries of the submatrix  sub( B )  need
*          not be set on input.
*
*  A       (local input) COMPLEX array
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
*  B       (local input/local output) COMPLEX array
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
   char           DiagA, DirB, OpC, OpR, SideOp, TopC, TopR, TranOp, UploA,
                  Var, ctop, ctopsave, rtop, rtopsave;
   int            Ai, Aj, Bi, Bj, ChooseAB, ForceTop, ctxt, info, itmp, lside,
                  mycol, myrow, nb, notran, nounit, npcol, nprow, upper;
   double         ABestL, ABestR, Best, tmp1, tmp2, tmp3, tmp4;
   PBTYP_T        * type;
/*
*  .. Local Arrays ..
*/
   int            Ad[DLEN_], Bd[DLEN_];
/* ..
*  .. Executable Statements ..
*
*/
   lside  = ( ( SideOp = Mupcase( F2C_CHAR( SIDE  )[0] ) ) ==   CLEFT );
   upper  = ( ( UploA  = Mupcase( F2C_CHAR( UPLO  )[0] ) ) ==  CUPPER );
   notran = ( ( TranOp = Mupcase( F2C_CHAR( TRANS )[0] ) ) == CNOTRAN );
   nounit = ( ( DiagA  = Mupcase( F2C_CHAR( DIAG  )[0] ) ) == CNOUNIT );
   PB_CargFtoC( *IA, *JA, DESCA, &Ai, &Aj, Ad );
   PB_CargFtoC( *IB, *JB, DESCB, &Bi, &Bj, Bd );
#ifndef NO_ARGCHK
/*
*  Test the input parameters
*/
   Cblacs_gridinfo( ( ctxt = Ad[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
   if( !( info = ( ( nprow == -1 ) ? -( 1101 + CTXT_ ) : 0 ) ) )
   {
      if( ( !lside ) && ( SideOp != CRIGHT ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PCTRSM", "Illegal SIDE = %c\n", SideOp );
         info = -1;
      }
      else if( ( !upper ) && ( UploA != CLOWER ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PCTRSM", "Illegal UPLO = %c\n", UploA );
         info = -2;
      }
      else if( ( !notran ) && ( TranOp != CTRAN ) && ( TranOp != CCOTRAN ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PCTRSM", "Illegal TRANS = %c\n", TranOp );
         info = -3;
      }
      else if( ( !nounit ) && ( DiagA != CUNIT ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PCTRSM", "Illegal DIAG = %c\n", DiagA );
         info = -4;
      }
      if( lside )
         PB_Cchkmat( ctxt, "PCTRSM", "A", *M, 5, *M, 5, Ai, Aj, Ad, 11,
                     &info );
      else
         PB_Cchkmat( ctxt, "PCTRSM", "A", *N, 6, *N, 6, Ai, Aj, Ad, 11,
                     &info );
      PB_Cchkmat(    ctxt, "PCTRSM", "B", *M, 5, *N, 6, Bi, Bj, Bd, 15,
                     &info );
   }
   if( info ) { PB_Cabort( ctxt, "PCTRSM", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( *M == 0 || *N == 0 ) return;
/*
*  Get type structure
*/
   type = PB_Cctypeset();
/*
*  And when alpha is zero
*/
   if( ( ALPHA[REAL_PART] == ZERO ) && ( ALPHA[IMAG_PART] == ZERO ) )
   {
      PB_Cplapad( type, ALL, NOCONJG, *M, *N, type->zero, type->zero,
                  ((char *) B), Bi, Bj, Bd );
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
*/
   nb = pilaenv_( &ctxt, C2F_CHAR( &type->type ) );
/*
*  ABestR, ABestL : both operands sub( A ) and sub( B ) are communicated
*                   ( N >> M when SIDE is left and M >> N otherwise )
*  Best           : only sub( B ) is communicated
*                   ( M >> N when SIDE is left and N >> M otherwise )
*/
   if( lside )
   {
      if( notran )
      {
         tmp1 = DNROC( *M, Ad[MB_], nprow ); tmp4 = DNROC( *N, Bd[NB_], npcol );
         ABestR = (double)(*M) *
           ( ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp1 / TWO ) +
             ( ( ( Bd[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp4 ) );
         itmp = MIN( Ad[MB_], Ad[NB_] );
         Best = (double)(*N) *
          ( (double)(CEIL( *M, itmp )) * (double)(itmp) *
            ( ( ( Ad[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : ONE ) +
            ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : ONE ) );
         ChooseAB = ( ABestR <= ( 2.0 * Best ) );
      }
      else
      {
         tmp1 = DNROC( *M, Ad[MB_], nprow ); tmp2 = DNROC( *M, Ad[NB_], npcol );
         tmp4 = DNROC( *N, Bd[NB_], npcol );
         ABestL = (double)(*M) *
           ( ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp1 / TWO ) +
             CBRATIO *
             ( ( ( Bd[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp4 ) );
         ABestR = (double)(*M) *
           ( ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp1 / TWO ) +
             ( ( ( Bd[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp4       ) +
             MAX( tmp2, tmp1 ) / TWO );
         itmp = MIN( Ad[MB_], Ad[NB_] );
         tmp2 = DNROC( *M, Ad[NB_], npcol ); tmp3 = DNROC( *M, Bd[MB_], nprow );
         Best = (double)(*N) *
          ( (double)(CEIL( *M, itmp )) * (double)(itmp) *
            ( ( ( ( Ad[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : ONE ) +
              ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : ONE ) ) +
            MAX( tmp2, tmp3 ) );
         ChooseAB = ( ( ABestL <= ( 2.0 * Best ) ) ||
                      ( ABestR <= ( 2.0 * Best ) ) );
      }
   }
   else
   {
      if( notran )
      {
         tmp2 = DNROC( *N, Ad[NB_], npcol ); tmp3 = DNROC( *M, Bd[MB_], nprow );
         ABestR = (double)(*N) *
           ( ( ( ( Ad[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp2 / TWO ) +
             ( ( ( Bd[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp3 ) );
         itmp = MIN( Ad[MB_], Ad[NB_] );
         Best = (double)(*M) *
          ( (double)(CEIL( *N, itmp )) * (double)(itmp) *
            ( ( ( Ad[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : ONE ) +
            ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : ONE ) );
         ChooseAB = ( ABestR <= ( 2.0 * Best ) );
      }
      else
      {
         tmp1 = DNROC( *N, Ad[MB_], nprow ); tmp2 = DNROC( *N, Ad[NB_], npcol );
         tmp3 = DNROC( *M, Bd[MB_], nprow );
         ABestL = (double)(*N) *
           ( ( ( ( Ad[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp2 / TWO ) +
             CBRATIO *
             ( ( ( Bd[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp3 ) );
         ABestR = (double)(*N) *
           ( ( ( ( Ad[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp2 / TWO ) +
             ( ( ( Bd[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp3       ) +
             MAX( tmp2, tmp1 ) / TWO );
         itmp = MIN( Ad[MB_], Ad[NB_] );
         tmp1 = DNROC( *N, Ad[MB_], nprow ); tmp4 = DNROC( *N, Bd[NB_], npcol );
         Best = (double)(*M) *
          ( (double)(CEIL( *N, itmp )) * (double)(itmp) *
            ( ( ( ( Ad[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : ONE ) +
              ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : ONE ) ) +
            MAX( tmp1, tmp4 ) );
         ChooseAB = ( ( ABestL <= ( 2.0 * Best ) ) ||
                      ( ABestR <= ( 2.0 * Best ) ) );
      }
   }

/*
* Var can remain uninitialized but is nevertheless used in PB_CptrsmAB.c
*  provide a default here. TODO: does this make sense ?
*==19891==    at 0x44F81B: PB_CptrsmAB (PB_CptrsmAB.c:538)
*==19891==    by 0x427BE7: pdtrsm_ (pdtrsm_.c:488)
*==19891==    by 0x405E46: MAIN_ (pdblas3tim.f:727)
*/
   Var = CRIGHT;

   if( ChooseAB )
   {
/*
*  BLACS topologies are enforced iff M and N are strictly greater than the
*  logical block size returned by pilaenv_. Otherwise, it is assumed that the
*  routine calling this routine has already selected an adequate topology.
*/
      ForceTop = ( ( *M > nb ) && ( *N > nb ) );

      if( ForceTop )
      {
         if( lside )
         {
            if( notran )
            {
               OpR = CBCAST; OpC = CBCAST; Var = CRIGHT;
               if( upper ) { TopR = TopC = CTOP_DRING; }
               else        { TopR = TopC = CTOP_IRING; }
            }
            else
            {
               if( ABestL <= ABestR )
               {
                  OpR = CBCAST; OpC = CCOMBINE; Var = CLEFT;
                  if( upper ) { TopR = TopC = CTOP_IRING; }
                  else        { TopR = TopC = CTOP_DRING; }
               }
               else
               {
                  OpR = CBCAST; OpC = CBCAST;  Var = CRIGHT;
                  if( upper ) { TopR = TopC = CTOP_IRING; }
                  else        { TopR = TopC = CTOP_DRING; }
               }
            }
         }
         else
         {
            if( notran )
            {
               OpR = CBCAST; OpC = CBCAST; Var = CRIGHT;
               if( upper ) { TopR = TopC = CTOP_IRING; }
               else        { TopR = TopC = CTOP_DRING; }
            }
            else
            {
               if( ABestL <= ABestR )
               {
                  OpR = CCOMBINE; OpC = CBCAST; Var = CLEFT;
                  if( upper ) { TopR = TopC = CTOP_DRING; }
                  else        { TopR = TopC = CTOP_IRING; }
               }
               else
               {
                  OpR = CBCAST; OpC = CBCAST;  Var = CRIGHT;
                  if( upper ) { TopR = TopC = CTOP_DRING; }
                  else        { TopR = TopC = CTOP_IRING; }
               }
            }
         }

         rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_GET );
         ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_GET );

         if( ( rtopsave = rtop ) != TopR )
            rtop = *PB_Ctop( &ctxt, &OpR, ROW,    &TopR );
         if( ( ctopsave = ctop ) != TopC )
            ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, &TopC );
/*
*  Remove the next 4 lines when the BLACS combine operations support ring
*  topologies
*/
         if( OpR == CCOMBINE )
            rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_DEFAULT );
         if( OpC == CCOMBINE )
            ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_DEFAULT );
      }

      PB_CptrsmAB( type, &Var, &SideOp, &UploA, &TranOp, &DiagA, *M, *N,
                   ((char *)ALPHA), ((char *)A), Ai, Aj, Ad, ((char *)B),
                   Bi, Bj, Bd );
/*
*  Restore the BLACS topologies when necessary.
*/
      if( ForceTop )
      {
         rtopsave = *PB_Ctop( &ctxt, &OpR, ROW,    &rtopsave );
         ctopsave = *PB_Ctop( &ctxt, &OpC, COLUMN, &ctopsave );
      }
   }
   else
   {
/*
*  BLACS topologies are always enforced.
*/
      if( ( lside && notran ) || ( !lside && !notran ) )
      {
         OpR = CCOMBINE; OpC = CBCAST;
         if( upper ) { TopR = TopC = CTOP_DRING; }
         else        { TopR = TopC = CTOP_IRING; }
/*
*  Remove the next line when the BLACS combine operations support ring
*  topologies
*/
         TopR = CTOP_DEFAULT;
      }
      else
      {
         OpR = CBCAST; OpC = CCOMBINE;
         if( upper ) { TopR = TopC = CTOP_IRING; }
         else        { TopR = TopC = CTOP_DRING; }
/*
*  Remove the next line when the BLACS combine operations support ring
*  topologies
*/
         TopC = CTOP_DEFAULT;
      }

      rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_GET );
      ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_GET );

      if( ( rtopsave = rtop ) != TopR )
         rtop = *PB_Ctop( &ctxt, &OpR, ROW,    &TopR );
      if( ( ctopsave = ctop ) != TopC )
         ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, &TopC );

      if( lside ) DirB = ( rtop == CTOP_DRING ? CBACKWARD : CFORWARD );
      else        DirB = ( ctop == CTOP_DRING ? CBACKWARD : CFORWARD );

      PB_CptrsmB( type, &DirB, &SideOp, &UploA, &TranOp, &DiagA, *M, *N,
                  ((char *)ALPHA), ((char *)A), Ai, Aj, Ad, ((char *)B),
                  Bi, Bj, Bd );
/*
*  Restore the BLACS topologies.
*/
      rtopsave = *PB_Ctop( &ctxt, &OpR, ROW,    &rtopsave );
      ctopsave = *PB_Ctop( &ctxt, &OpC, COLUMN, &ctopsave );
   }
/*
*  End of PCTRSM
*/
}
