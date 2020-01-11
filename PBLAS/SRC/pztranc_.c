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
void pztranc_( Int * M, Int * N,
               double * ALPHA,
               double * A, Int * IA, Int * JA, Int * DESCA,
               double * BETA,
               double * C, Int * IC, Int * JC, Int * DESCC )
#else
void pztranc_( M, N, ALPHA, A, IA, JA, DESCA, BETA, C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   Int            * IA, * IC, * JA, * JC, * M, * N;
   double         * ALPHA, * BETA;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA, * DESCC;
   double         * A, * C;
#endif
{
/*
*  Purpose
*  =======
*
*  PZTRANC  transposes a matrix
*
*     sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )
*
*  where
*
*     sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),
*
*     sub( A ) denotes A(IA:IA+N-1,JA:JA+M-1), and,
*
*     op( X ) = conjg( X )'.
*
*  Thus, op( sub( A ) ) denotes conjg( A(IA:IA+N-1,JA:JA+M-1)' ).
*
*  Beta is a scalar, sub( C ) is an m by n submatrix, and sub( A ) is an
*  n by m submatrix.
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
*  ALPHA   (global input) COMPLEX*16
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of  the array  A
*          corresponding to the entries of the submatrix  sub( A )  need
*          not be set on input.
*
*  A       (local input) COMPLEX*16 array
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+M-1 ).  Before  entry, this array contains
*          the local entries of the matrix A.
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
*  BETA    (global input) COMPLEX*16
*          On entry,  BETA  specifies the scalar  beta.   When  BETA  is
*          supplied  as  zero  then  the  local entries of  the array  C
*          corresponding to the entries of the submatrix  sub( C )  need
*          not be set on input.
*
*  C       (local input/local output) COMPLEX*16 array
*          On entry, C is an array of dimension (LLD_C, Kc), where Kc is
*          at least Lc( 1, JC+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix C.
*          On exit, the entries of this array corresponding to the local
*          entries of the submatrix  sub( C )  are  overwritten  by  the
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
   Int            Ai, Aj, Ci, Cj, ctxt, info, mycol, myrow, npcol, nprow;
/*
*  .. Local Arrays ..
*/
   Int            Ad[DLEN_], Cd[DLEN_];
/* ..
*  .. Executable Statements ..
*
*/
   PB_CargFtoC( *IA, *JA, DESCA, &Ai, &Aj, Ad );
   PB_CargFtoC( *IC, *JC, DESCC, &Ci, &Cj, Cd );
#ifndef NO_ARGCHK
/*
*  Test the input parameters
*/
   Cblacs_gridinfo( ( ctxt = Ad[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
   if( !( info = ( ( nprow == -1 ) ? -( 701 + CTXT_ ) : 0 ) ) )
   {
      PB_Cchkmat( ctxt, "PZTRANC", "A", *N, 2, *M, 1, Ai, Aj, Ad,  7, &info );
      PB_Cchkmat( ctxt, "PZTRANC", "C", *M, 1, *N, 2, Ci, Cj, Cd, 12, &info );
   }
   if( info ) { PB_Cabort( ctxt, "PZTRANC", info ); return; }
#endif
/*
*  Quick return if possible
*/
   if( ( *M == 0 ) || ( *N == 0 ) ||
       ( ( ALPHA[REAL_PART] == ZERO && ALPHA[IMAG_PART] == ZERO ) &&
         ( BETA [REAL_PART] ==  ONE && BETA [IMAG_PART] == ZERO ) ) )
      return;
/*
*  And when alpha is zero
*/
   if( ( ALPHA[REAL_PART] == ZERO ) && ( ALPHA[IMAG_PART] == ZERO ) )
   {
      if( ( BETA[REAL_PART] == ZERO ) && ( BETA[IMAG_PART] == ZERO ) )
      {
         PB_Cplapad( PB_Cztypeset(), ALL, NOCONJG, *M, *N, ((char *)BETA),
                     ((char *)BETA), ((char *) C), Ci, Cj, Cd );
      }
      else
      {
         PB_Cplascal( PB_Cztypeset(), ALL, NOCONJG, *M, *N, ((char *)BETA),
                      ((char * )C), Ci, Cj, Cd );
      }
      return;
   }
/*
*  Start the operations
*/
   PB_Cptran( PB_Cztypeset(), CONJG,   *M, *N, ((char *) ALPHA),
              ((char *) A), Ai, Aj, Ad, ((char *)  BETA), ((char *) C),
              Ci, Cj, Cd );
/*
*  End of PZTRANC
*/
}
