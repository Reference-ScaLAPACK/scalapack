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
void pcherk_( F_CHAR_T UPLO, F_CHAR_T TRANS, int * N, int * K,
              float * ALPHA,
              float * A, int * IA, int * JA, int * DESCA,
              float * BETA,
              float * C, int * IC, int * JC, int * DESCC )
#else
void pcherk_( UPLO, TRANS, N, K, ALPHA, A, IA, JA, DESCA, BETA,
              C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   F_CHAR_T       TRANS, UPLO;
   int            * IA, * IC, * JA, * JC, * K, * N;
   float          * ALPHA, * BETA;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCC;
   float          * A, * C;
#endif
{
/*
*  Purpose
*  =======
*
*  PCHERK  performs one of the Hermitian rank k operations
*
*     sub( C ) := alpha*sub( A )*conjg( sub( A )' ) + beta*sub( C ),
*
*  or
*
*     sub( C ) := alpha*conjg( sub( A )' )*sub( A ) + beta*sub( C ),
*
*  where
*
*     sub( C ) denotes C(IC:IC+N-1,JC:JC+N-1), and,
*
*     sub( A ) denotes A(IA:IA+N-1,JA:JA+K-1)  if TRANS = 'N',
*                      A(IA:IA+K-1,JA:JA+N-1)  otherwise.
*
*  Alpha and beta are real scalars,  sub( C )  is  an  n by n  Hermitian
*  submatrix and sub( A ) is an n by k submatrix in the first case and a
*  k by n submatrix in the second case.
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
*  UPLO    (global input) CHARACTER*1
*          On  entry,   UPLO  specifies  whether  the  local  pieces  of
*          the array  C  containing the  upper or lower triangular  part
*          of the Hermitian submatrix  sub( C )  are to be referenced as
*          follows:
*
*             UPLO = 'U' or 'u'   Only the local pieces corresponding to
*                                 the   upper  triangular  part  of  the
*                                 Hermitian submatrix sub( C ) are to be
*                                 referenced,
*
*             UPLO = 'L' or 'l'   Only the local pieces corresponding to
*                                 the   lower  triangular  part  of  the
*                                 Hermitian submatrix sub( C ) are to be
*                                 referenced.
*
*  TRANS   (global input) CHARACTER*1
*          On entry,  TRANS  specifies the  operation to be performed as
*          follows:
*
*             TRANS = 'N' or 'n'
*                  sub( C ) := alpha*sub( A )*conjg( sub( A )' ) +
*                              beta*sub( C ),
*
*             TRANS = 'C' or 'c'
*                  sub( C ) := alpha*conjg( sub( A )' )*sub( A ) +
*                              beta*sub( C ).
*
*  N       (global input) INTEGER
*          On entry,  N specifies the order of the  submatrix  sub( C ).
*          N must be at least zero.
*
*  K       (global input) INTEGER
*          On entry, with TRANS = 'N' or 'n',  K specifies the number of
*          columns  of  the submatrix  sub( A ), and with TRANS = 'C' or
*          'c',  K  specifies  the  number  of  rows  of  the  submatrix
*          sub( A ). K  must be at least zero.
*
*  ALPHA   (global input) REAL
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of  the array  A
*          corresponding to the entries of the submatrix  sub( A )  need
*          not be set on input.
*
*  A       (local input) COMPLEX array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least  Lc( 1, JA+K-1 ) when  TRANS = 'N' or 'n', and is at
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
*  BETA    (global input) REAL
*          On entry,  BETA  specifies the scalar  beta.   When  BETA  is
*          supplied  as  zero  then  the  local entries of  the array  C
*          corresponding to the entries of the submatrix  sub( C )  need
*          not be set on input.
*
*  C       (local input/local output) COMPLEX array
*          On entry, C is an array of dimension (LLD_C, Kc), where Kc is
*          at least Lc( 1, JC+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix C.
*          Before  entry  with  UPLO = 'U' or 'u', this  array  contains
*          the local entries corresponding to the upper triangular  part
*          of the  Hermitian  submatrix  sub( C ), and the local entries
*          corresponding to the  strictly lower triangular  of  sub( C )
*          are not  referenced.  On exit,  the upper triangular part  of
*          sub( C ) is overwritten by the  upper triangular part  of the
*          updated submatrix.
*          Before  entry  with  UPLO = 'L' or 'l', this  array  contains
*          the local entries corresponding to the lower triangular  part
*          of the  Hermitian  submatrix  sub( C ), and the local entries
*          corresponding to the  strictly upper triangular  of  sub( C )
*          are not  referenced.  On exit,  the lower triangular part  of
*          sub( C ) is overwritten by the  lower triangular part  of the
*          updated submatrix.
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
   char           DirA, OpC, OpR, TopC, TopR, TranOp, UploC, ctop, ctopsave,
                  rtop, rtopsave;
   int            Ai, Aj, ChooseAC, Ci, Cj, ForceTop, ctxt, info, mycol,
                  myrow, nb, notran, npcol, nprow, upper;
   double         Aest, ACest, tmp1, tmp2, tmp3, tmp4;
   cmplx          Calph;
   PBTYP_T        * type;
/*
*  .. Local Arrays ..
*/
   int            Ad[DLEN_], Cd[DLEN_];
/* ..
*  .. Executable Statements ..
*
*/
   upper  = ( ( UploC  = Mupcase( F2C_CHAR( UPLO  )[0] ) ) ==  CUPPER );
   notran = ( ( TranOp = Mupcase( F2C_CHAR( TRANS )[0] ) ) == CNOTRAN );
   PB_CargFtoC( *IA, *JA, DESCA, &Ai, &Aj, Ad );
   PB_CargFtoC( *IC, *JC, DESCC, &Ci, &Cj, Cd );
#ifndef NO_ARGCHK
/*
*  Test the input parameters
*/
   Cblacs_gridinfo( ( ctxt = Ad[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
   if( !( info = ( ( nprow == -1 ) ? -( 901 + CTXT_ ) : 0 ) ) )
   {
      if( ( !upper ) && ( UploC != CLOWER ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PCHERK", "Illegal UPLO = %c\n", UploC );
         info = -1;
      }
      else if( ( !notran ) && ( TranOp != CCOTRAN ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PCHERK", "Illegal TRANS = %c\n", TranOp );
         info = -2;
      }
      if( notran )
         PB_Cchkmat( ctxt, "PCHERK", "A", *N, 3, *K, 4, Ai, Aj, Ad,  9,
                     &info );
      else
         PB_Cchkmat( ctxt, "PCHERK", "A", *K, 4, *N, 3, Ai, Aj, Ad,  9,
                     &info );
      PB_Cchkmat(    ctxt, "PCHERK", "C", *N, 3, *N, 3, Ci, Cj, Cd, 14,
                     &info );
   }
   if( info ) { PB_Cabort( ctxt, "PCHERK", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( ( *N == 0 ) ||
       ( ( ( ALPHA[REAL_PART] == ZERO ) || ( *K == 0 ) ) &&
         ( BETA[REAL_PART] == ONE                      ) ) )
      return;
/*
*  Get type structure
*/
   type = PB_Cctypeset();
/*
*  And when alpha or K is zero
*/
   if( ( ALPHA[REAL_PART] == ZERO ) || ( *K == 0 ) )
   {
      if( BETA[REAL_PART] == ZERO )
      {
         PB_Cplapad( type, &UploC, NOCONJG, *N, *N, type->zero, type->zero,
                     ((char *) C), Ci, Cj, Cd );
      }
      else
      {
         PB_Cplascal( type, &UploC, CONJG,   *N, *N, ((char *) BETA),
                      ((char *) C), Ci, Cj, Cd );
      }
      return;
   }
/*
*  Start the operations
*/
#ifdef NO_ARGCHK
   Cblacs_gridinfo( ( ctxt = Ad[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
#endif
   Calph[REAL_PART] = ALPHA[REAL_PART];
   Calph[IMAG_PART] = ZERO;
/*
*  Algorithm selection is based on approximation of the communication volume
*  for distributed and aligned operands.
*
*  ACest: both operands sub( A ) and sub( C ) are communicated (K >> N)
*  Aest : only sub( A ) is communicated                        (N >> K)
*/
   if( notran )
   {
      tmp1  = DNROC( *N, Cd[MB_], nprow ); tmp3 = DNROC( *K, Ad[NB_], npcol );
      ACest = (double)(*N) *
        ( ( ( ( Ad[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp3 ) +
          ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO :
            CBRATIO * tmp1 / TWO ) );
      tmp1  = DNROC( *N, Cd[MB_], nprow ); tmp2 = DNROC( *N, Cd[NB_], npcol );
                                           tmp4 = DNROC( *N, Ad[MB_], nprow );
      Aest  = (double)(*K) *
              ( ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp1 ) +
                ( nprow == 1 ? ZERO : tmp2 ) + MAX( tmp2, tmp4 ) );
   }
   else
   {
      tmp2  = DNROC( *N, Cd[NB_], npcol ); tmp4 = DNROC( *K, Ad[MB_], nprow );
      ACest = (double)(*N) *
        ( ( ( ( Ad[CSRC_] == -1 ) || ( npcol == 1 ) ) ? ZERO : tmp4 ) +
          ( ( ( Ad[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO :
            CBRATIO * tmp2 / TWO ) );
      tmp1  = DNROC( *N, Cd[MB_], nprow ); tmp2 = DNROC( *N, Cd[NB_], npcol );
      tmp3  = DNROC( *N, Ad[NB_], npcol );
      Aest  = (double)(*K) *
              ( ( ( ( Ad[RSRC_] == -1 ) || ( nprow == 1 ) ) ? ZERO : tmp2 ) +
                ( npcol == 1 ? ZERO : tmp1 ) + MAX( tmp1, tmp3 ) );
   }
/*
*  Shift a little the cross-over point between both algorithms.
*/
   ChooseAC = ( ( 1.3 * ACest ) <= Aest );
/*
*  BLACS topologies are enforced iff N and K are strictly greater than the
*  logical block size returned by pilaenv_. Otherwise, it is assumed that the
*  routine calling this routine has already selected an adequate topology.
*/
   nb       = pilaenv_( &ctxt, C2F_CHAR( &type->type ) );
   ForceTop = ( ( *N > nb ) && ( *K > nb ) );

   if( ChooseAC )
   {
      if( notran )
      {
         OpC  = CBCAST;
         ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_GET );

         if( ForceTop )
         {
            OpR  = CCOMBINE;
            rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_GET );

            rtopsave = rtop;
            ctopsave = ctop;

            if( upper ) { TopR = CTOP_IRING; TopC = CTOP_DRING; }
            else        { TopR = CTOP_DRING; TopC = CTOP_IRING; }

            ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, &TopC );
            rtop = *PB_Ctop( &ctxt, &OpR, ROW,    &TopR );
/*
*  Remove the next line when the BLACS combine operations support ring
*  topologies
*/
            rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_DEFAULT );
         }

         DirA = ( ctop == CTOP_DRING ? CBACKWARD : CFORWARD );
      }
      else
      {
         OpR  = CBCAST;
         rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_GET );

         if( ForceTop )
         {
            OpC  = CCOMBINE;
            ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_GET );

            rtopsave = rtop;
            ctopsave = ctop;

            if( upper ) { TopR = CTOP_IRING; TopC = CTOP_DRING; }
            else        { TopR = CTOP_DRING; TopC = CTOP_IRING; }

            rtop = *PB_Ctop( &ctxt, &OpR, ROW,    &TopR );
            ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, &TopC );
/*
*  Remove the next line when the BLACS combine operations support ring
*  topologies
*/
            ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_DEFAULT );
         }
         DirA = ( rtop == CTOP_DRING ? CBACKWARD : CFORWARD );
      }

      PB_CpsyrkAC( type, &DirA, CONJG,   &UploC, ( notran ? NOTRAN : COTRAN ),
                   *N, *K, ((char *)Calph), ((char *)A), Ai, Aj, Ad,
                   ((char *)BETA), ((char *)C), Ci, Cj, Cd );
   }
   else
   {
      if( notran )
      {
         OpR  = CBCAST;
         rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_GET );

         if( ForceTop )
         {
            OpC  = CBCAST;
            ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_GET );

            rtopsave = rtop;
            ctopsave = ctop;
/*
*  No clear winner for the ring topologies, so that if a ring topology is
*  already selected, keep it.
*/
            if( ( rtop != CTOP_DRING ) && ( rtop != CTOP_IRING ) &&
                ( rtop != CTOP_SRING ) )
               rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_SRING );
            if( ( ctop != CTOP_DRING ) && ( ctop != CTOP_IRING ) &&
                ( ctop != CTOP_SRING ) )
               ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_SRING );
         }

         DirA = ( rtop == CTOP_DRING ? CBACKWARD : CFORWARD );
      }
      else
      {
         OpC  = CBCAST;
         ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_GET );

         if( ForceTop )
         {
            OpR  = CBCAST;
            rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_GET );

            rtopsave = rtop;
            ctopsave = ctop;
/*
*  No clear winner for the ring topologies, so that if a ring topology is
*  already selected, keep it.
*/
            if( ( rtop != CTOP_DRING ) && ( rtop != CTOP_IRING ) &&
                ( rtop != CTOP_SRING ) )
               rtop = *PB_Ctop( &ctxt, &OpR, ROW,    TOP_SRING );
            if( ( ctop != CTOP_DRING ) && ( ctop != CTOP_IRING ) &&
                ( ctop != CTOP_SRING ) )
               ctop = *PB_Ctop( &ctxt, &OpC, COLUMN, TOP_SRING );
         }

         DirA = ( ctop == CTOP_DRING ? CBACKWARD : CFORWARD );
      }

      PB_CpsyrkA( type, &DirA, CONJG,   &UploC, ( notran ? NOTRAN : COTRAN ),
                  *N, *K, ((char *)Calph), ((char *)A), Ai, Aj, Ad,
                  ((char *)BETA), ((char *)C), Ci, Cj, Cd );
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
*  End of PCHERK
*/
}
