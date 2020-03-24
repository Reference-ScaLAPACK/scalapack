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
void pstradd_( F_CHAR_T UPLO, F_CHAR_T TRANS, Int * M, Int * N,
               float * ALPHA,
               float * A, Int * IA, Int * JA, Int * DESCA,
               float * BETA,
               float * C, Int * IC, Int * JC, Int * DESCC )
#else
void pstradd_( UPLO, TRANS, M, N, ALPHA, A, IA, JA, DESCA, BETA,
               C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   F_CHAR_T       TRANS, UPLO;
   Int            * IA, * IC, * JA, * JC, * M, * N;
   float          * ALPHA, * BETA;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA, * DESCC;
   float          * A, * C;
#endif
{
/*
*  Purpose
*  =======
*
*  PSTRADD  adds a trapezoidal matrix to another
*
*     sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )
*
*  where
*
*     sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),  and, op( X )  is one  of
*
*     op( X ) = X   or   op( X ) = X'.
*
*  Thus, op( sub( A ) ) denotes A(IA:IA+M-1,JA:JA+N-1)   if TRANS = 'N',
*                               A(IA:IA+N-1,JA:JA+M-1)'  if TRANS = 'T',
*                               A(IA:IA+N-1,JA:JA+M-1)'  if TRANS = 'C',
*
*  Alpha  and  beta  are scalars, sub( C ) and op( sub( A ) ) are m by n
*  upper or lower trapezoidal submatrices.
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
*          On  entry, UPLO specifies whether  the  local  pieces  of the
*          array C containing the upper or lower triangular part  of the
*          triangular submatrix sub( C ) is to be referenced as follows:
*
*             UPLO = 'U' or 'u'   Only the local pieces corresponding to
*                                 the   upper  triangular  part  of  the
*                                 triangular submatrix sub( C ) is to be
*                                 referenced,
*
*             UPLO = 'L' or 'l'   Only the local pieces corresponding to
*                                 the   lower  triangular  part  of  the
*                                 triangular submatrix sub( C ) is to be
*                                 referenced.
*
*  TRANS   (global input) CHARACTER*1
*          On entry,  TRANS   specifies the form of op( sub( A ) ) to be
*          used in the matrix addition as follows:
*
*             TRANS = 'N' or 'n'   op( sub( A ) ) = sub( A ),
*
*             TRANS = 'T' or 't'   op( sub( A ) ) = sub( A )',
*
*             TRANS = 'C' or 'c'   op( sub( A ) ) = sub( A )'.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( C ) and the number of columns of the submatrix sub( A ).
*          M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( C ) and the number of rows of the submatrix sub( A ).  N
*          must be at least zero.
*
*  ALPHA   (global input) REAL
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of  the array  A
*          corresponding to the entries of the submatrix  sub( A )  need
*          not be set on input.
*
*  A       (local input) REAL array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ) when  TRANS = 'N' or 'n'  and  is at
*          least Lc( 1, JA+M-1 ) otherwise.  Before  entry,  this  array
*          contains the local entries of the matrix A.
*          Before entry with UPLO = 'U' or 'u' and TRANS = 'N' or 'n' or
*          UPLO = 'L' or 'l' and  TRANS = 'T', 'C', 't' or 'c', this ar-
*          ray contains the local entries corresponding to  the  entries
*          of the upper triangular submatrix sub( A ), and the local en-
*          tries corresponding to the entries of the strictly lower tri-
*          angular part  of the submatrix  sub( A )  are not referenced.
*          Before entry with UPLO = 'L' or 'l' and TRANS = 'N' or 'n' or
*          UPLO = 'U' or 'u' and  TRANS = 'T', 'C', 't' or 'c', this ar-
*          ray contains the local entries corresponding  to the  entries
*          of the lower triangular submatrix sub( A ), and the local en-
*          tries corresponding to the entries of the strictly upper tri-
*          angular part  of the submatrix  sub( A )  are not referenced.
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
*  C       (local input/local output) REAL array
*          On entry, C is an array of dimension (LLD_C, Kc), where Kc is
*          at least Lc( 1, JC+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix C.
*          Before  entry  with  UPLO = 'U' or 'u', this  array  contains
*          the local entries corresponding to the upper triangular  part
*          of the  triangular submatrix  sub( C ), and the local entries
*          corresponding to the  strictly lower triangular  of  sub( C )
*          are not  referenced.  On exit,  the upper triangular part  of
*          sub( C ) is overwritten by the  upper triangular part  of the
*          updated submatrix.
*          Before  entry  with  UPLO = 'L' or 'l', this  array  contains
*          the local entries corresponding to the lower triangular  part
*          of the  triangular submatrix  sub( C ), and the local entries
*          corresponding to the  strictly upper triangular  of  sub( C )
*          are not  referenced.  On exit,  the lower triangular part  of
*          sub( C ) is overwritten by the  lower triangular part  of the
*          updated submatrix.
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
   char           DirAC, TranOp, UploC, ctop, rtop;
   Int            Ai, Aj, Ci, Cj, ctxt, info, mycol, myrow, notran, npcol,
                  nprow, upper;
/*
*  .. Local Arrays ..
*/
   Int            Ad[DLEN_], Cd[DLEN_];
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
         PB_Cwarn( ctxt, __LINE__, "PSTRADD", "Illegal UPLO = %c\n", UploC );
         info = -1;
      }
      else if( ( !notran ) && ( TranOp != CTRAN ) && ( TranOp != CCOTRAN ) )
      {
         PB_Cwarn( ctxt, __LINE__, "PSTRADD", "Illegal TRANS = %c\n", TranOp );
         info = -2;
      }
      if( notran )
         PB_Cchkmat( ctxt, "PSTRADD", "A", *M, 3, *N, 4, Ai, Aj, Ad,  9,
                     &info );
      else
         PB_Cchkmat( ctxt, "PSTRADD", "A", *N, 4, *M, 3, Ai, Aj, Ad,  9,
                     &info );
      PB_Cchkmat(    ctxt, "PSTRADD", "C", *M, 3, *N, 4, Ci, Cj, Cd, 14,
                     &info );
   }
   if( info ) { PB_Cabort( ctxt, "PSTRADD", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( ( *M == 0 ) || ( *N == 0 ) ||
       ( ( ALPHA[REAL_PART] == ZERO ) && ( BETA[REAL_PART] == ONE ) ) )
      return;
/*
*  And when alpha is zero
*/
   if( ALPHA[REAL_PART] == ZERO )
   {
      if( BETA[REAL_PART] == ZERO )
      {
         PB_Cplapad( PB_Cstypeset(), &UploC, NOCONJG, *M, *N,
                     ((char *)BETA), ((char *)BETA), ((char *) C), Ci, Cj, Cd );
      }
      else
      {
         PB_Cplascal( PB_Cstypeset(), &UploC, NOCONJG, *M, *N,
                      ((char *)BETA), ((char * )C), Ci, Cj, Cd );
      }
      return;
   }
/*
*  Start the operations
*/
/*
*  This operation mainly involves point-to-point send and receive communication.
*  There is therefore no particular BLACS topology to recommend. Still, one can
*  choose the main loop direction in which the operands will be added, but not
*  transposed. This selection is based on the current setting for the BLACS
*  broadcast operations.
*/
   rtop = *PB_Ctop( &ctxt, BCAST, ROW,    TOP_GET );
   ctop = *PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );

   if( *M <= *N )
      DirAC = ( rtop == CTOP_DRING ? CBACKWARD : CFORWARD );
   else
      DirAC = ( ctop == CTOP_DRING ? CBACKWARD : CFORWARD );
   PB_Cptradd( PB_Cstypeset(), &DirAC, &UploC, ( notran ? NOTRAN :
               TRAN ), *M, *N, ((char *) ALPHA), ((char *) A), Ai, Aj, Ad,
               ((char *)  BETA), ((char *) C), Ci, Cj, Cd );
/*
*  End of PSTRADD
*/
}
