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
void PB_Cptradd( PBTYP_T * TYPE, char * DIRECAC, char * UPLO, char * TRANS,
                 int M, int N, char * ALPHA, char * A, int IA, int JA,
                 int * DESCA, char * BETA, char * C, int IC, int JC,
                 int * DESCC )
#else
void PB_Cptradd( TYPE, DIRECAC, UPLO, TRANS, M, N, ALPHA, A, IA, JA, DESCA,
                 BETA, C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * DIRECAC, * TRANS, * UPLO;
   int            IA, IC, JA, JC, M, N;
   char           * ALPHA, * BETA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCC;
   char           * A, * C;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cptradd  adds a trapezoidal matrix to another
*
*     sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )
*
*  where
*
*     sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),  and, op( X )  is one  of
*
*     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ).
*
*  Thus, op( sub( A ) ) denotes A(IA:IA+M-1,JA:JA+N-1)   if TRANS = 'N',
*                               A(IA:IA+N-1,JA:JA+M-1)'  if TRANS = 'T',
*                        conjg(A(IA:IA+N-1,JA:JA+M-1)')  if TRANS = 'C',
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
*  TYPE    (local input) pointer to a PBTYP_T structure
*          On entry,  TYPE  is a pointer to a structure of type PBTYP_T,
*          that contains type information (See pblas.h).
*
*  DIRECAC (global input) pointer to CHAR
*          On entry,  DIRECAC specifies  the direction in which the rows
*          or columns of sub( A ) and sub( C ) should be looped over  as
*          follows:
*             DIRECA = 'F' or 'f'   forward  or increasing,
*             DIRECA = 'B' or 'b'   backward or decreasing.
*
*  UPLO    (global input) pointer to CHAR
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
*  TRANS   (global input) pointer to CHAR
*          On entry,  TRANS   specifies the form of op( sub( A ) ) to be
*          used in the matrix addition as follows:
*
*             TRANS = 'N' or 'n'   op( sub( A ) ) = sub( A ),
*
*             TRANS = 'T' or 't'   op( sub( A ) ) = sub( A )',
*
*             TRANS = 'C' or 'c'   op( sub( A ) ) = conjg( sub( A )' ).
*
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of the submatrices
*          sub( A ) and sub( C ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatri-
*          ces sub( A ) and sub( C ). N must be at least zero.
*
*  ALPHA   (global input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of  the array  A
*          corresponding to the entries of the submatrix  sub( A )  need
*          not be set on input.
*
*  A       (local input) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
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
   char           Dir, * one, * zero;
   int            Afr, conjg, k, kb, kbb, kend, kstart, kstep, ktmp;
/*
*  .. Local Arrays ..
*/
   int            DBUFA[DLEN_];
   char           * Aptr = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  sub( C ) := beta * sub( C )
*/
   PB_Cplascal( TYPE, UPLO, NOCONJG, M, N, BETA, C, IC, JC, DESCC );

   one  = TYPE->one; zero  = TYPE->zero;
   kb   = pilaenv_( &DESCC[CTXT_], C2F_CHAR( &TYPE->type ) );

   if( Mupcase( DIRECAC[0] ) == CFORWARD )
   {
      Dir = CFORWARD;
      kstart = 0; kend = ( ( MIN( M, N ) - 1 ) / kb + 1 ) * kb; kstep = kb;
   }
   else
   {
      Dir = CBACKWARD;
      kstart = ( ( MIN( M, N ) - 1 ) / kb ) * kb;  kend = kstep = -kb;
   }

   if( Mupcase( TRANS[0] ) == CNOTRAN )
   {
      if( Mupcase( UPLO [0] ) ==  CUPPER )
      {
         if( M >= N )
         {
            for( k = kstart; k != kend; k += kstep )
            {
               kbb = N - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA:IA+k+kbb-1, JA+k:JA+k+kbb-1 )
*/
               PB_CGatherV( TYPE, ALLOCATE, &Dir, ktmp, kbb, A, IA, JA+k, DESCA,
                            COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Scale A( IA:IA+k+kbb-1, JA+k:JA+k+kbb-1 ) by ALPHA
*/
               PB_Cplascal( TYPE, ALL, NOCONJG, ktmp, kbb, ALPHA, Aptr, 0, 0,
                            DBUFA );
/*
*  Zero lower triangle of A( IA:IA+k+kbb-1, JA+k:JA+k+kbb-1 )
*/
               if( kbb > 1 )
                  PB_Cplapad( TYPE, LOWER, NOCONJG, kbb-1, kbb-1, zero, zero,
                              Aptr, k+1, 0, DBUFA );
/*
*  C( IC:IC+k+kbb-1, JC+k:JC+k+kbb-1 ) += A( IA:IA+k+kbb-1, JA+k:JA+k+kbb-1 )
*/
               PB_CScatterV( TYPE, &Dir, ktmp, kbb, Aptr, 0, 0, DBUFA, COLUMN,
                             one, C, IC, JC+k, DESCC, COLUMN );

               if( Afr ) free( Aptr );
            }
         }
         else
         {
            for( k = kstart; k != kend; k += kstep )
            {
               kbb = M - k; kbb = MIN( kbb, kb ); ktmp = N - k;
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 )
*/
               PB_CGatherV( TYPE, ALLOCATE, &Dir, kbb, ktmp, A, IA+k, JA+k,
                            DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Scale A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 ) by ALPHA
*/
               PB_Cplascal( TYPE, ALL, NOCONJG, kbb, ktmp, ALPHA, Aptr, 0, 0,
                            DBUFA );
/*
*  Zero lower triangle of A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 )
*/
               if( kbb > 1 )
                  PB_Cplapad( TYPE, LOWER, NOCONJG, kbb-1, kbb-1, zero, zero,
                              Aptr, 1, 0, DBUFA );
/*
*  C( IC+k:IC+k+kbb-1, JC+k:JC+N-1 ) += A( IA+k:IA+k+kbb-1, JA+k:JA+N-1 )
*/
               PB_CScatterV( TYPE, &Dir, kbb, ktmp, Aptr, 0, 0, DBUFA, ROW,
                             one, C, IC+k, JC+k, DESCC, ROW );

               if( Afr ) free( Aptr );
            }
         }
      }
      else
      {
         if( M >= N )
         {
            for( k = kstart; k != kend; k += kstep )
            {
               kbb = N - k; kbb = MIN( kbb, kb ); ktmp = M - k;
/*
*  Accumulate A( IA+k:IA+M-1, JA+k:JA+k+kbb-1 )
*/
               PB_CGatherV( TYPE, ALLOCATE, &Dir, ktmp, kbb, A, IA+k, JA+k,
                            DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Scale A( IA:IA+k+kbb-1, JA+k:JA+k+kbb-1 ) by ALPHA
*/
               PB_Cplascal( TYPE, ALL, NOCONJG, ktmp, kbb, ALPHA, Aptr, 0, 0,
                            DBUFA );
/*
*  Zero upper triangle of A( IA+k:IA+M-1, JA+k:JA+k+kbb-1 )
*/
               if( kbb > 1 )
                  PB_Cplapad( TYPE, UPPER, NOCONJG, kbb-1, kbb-1, zero, zero,
                              Aptr, 0, 1, DBUFA );
/*
*  C( IC:IC+k+kbb-1, JC+k:JC+k+kbb-1 ) += A( IA:IA+k+kbb-1, JA+k:JA+k+kbb-1 )
*/
               PB_CScatterV( TYPE, &Dir, ktmp, kbb, Aptr, 0, 0, DBUFA, COLUMN,
                             one, C, IC+k, JC+k, DESCC, COLUMN );

               if( Afr ) free( Aptr );
            }
         }
         else
         {
            for( k = kstart; k != kend; k += kstep )
            {
               kbb = M - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA:JA+k+kbb-1 )
*/
               PB_CGatherV( TYPE, ALLOCATE, &Dir, kbb, ktmp, A, IA+k, JA,
                            DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Scale A( IA+k:IA+k+kbb-1, JA:JA+k+kbb-1 ) by ALPHA
*/
               PB_Cplascal( TYPE, ALL, NOCONJG, kbb, ktmp, ALPHA, Aptr, 0, 0,
                            DBUFA );
/*
*  Zero upper triangle of A( IA+k:IA+k+kbb-1, JA+k:JA:JA+k+kbb-1 )
*/
               if( kbb > 1 )
                  PB_Cplapad( TYPE, UPPER, NOCONJG, kbb-1, kbb-1, zero, zero,
                              Aptr, 0, k+1, DBUFA );
/*
*  C( IC+k:IC+k+kbb-1, JC:JC+k+kbb-1 ) += A( IA+k:IA+k+kbb-1, JA:JA+k+kbb-1 )
*/
               PB_CScatterV( TYPE, &Dir, kbb, ktmp, Aptr, 0, 0, DBUFA, ROW,
                             one, C, IC+k, JC, DESCC, ROW );

               if( Afr ) free( Aptr );
            }
         }
      }
   }
   else
   {
      conjg = ( Mupcase( TRANS[0] ) == CCOTRAN );

      if( Mupcase( UPLO [0] ) ==  CUPPER )
      {
         if( M >= N )
         {
            for( k = kstart; k != kend; k += kstep )
            {
               kbb = N - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA:JA+k+kbb-1 )
*/
               PB_CGatherV( TYPE, ALLOCATE, &Dir, kbb, ktmp, A, IA+k, JA, DESCA,
                            ROW, &Aptr, DBUFA, &Afr );
/*
*  Scale A( IA+k:IA+k+kbb-1, JA:JA+k+kbb-1 ) by ALPHA
*/
               if( conjg )
                  PB_Cplacnjg( TYPE, kbb, ktmp, ALPHA, Aptr, 0, 0, DBUFA );
               else
                  PB_Cplascal( TYPE, ALL, NOCONJG, kbb, ktmp, ALPHA, Aptr,
                               0, 0, DBUFA );
/*
*  Zero upper triangle of A( IA+k:IA+k+kbb-1, JA+k:JA+k+kbb-1 )
*/
               if( kbb > 1 )
                  PB_Cplapad( TYPE, UPPER, NOCONJG, kbb-1, kbb-1, zero, zero,
                              Aptr, 0, k+1, DBUFA );
/*
*  C( IC:IC+k+kbb-1, JC+k:JC+k+kbb-1 ) += A( IA+k:IA+k+kbb-1, JA:JA+k+kbb-1 )'
*/
               PB_CScatterV( TYPE, &Dir, kbb, ktmp, Aptr, 0, 0, DBUFA, ROW,
                             one, C, IC, JC+k, DESCC, COLUMN );

               if( Afr ) free( Aptr );
            }
         }
         else
         {
            for( k = kstart; k != kend; k += kstep )
            {
               kbb = M - k; kbb = MIN( kbb, kb ); ktmp = N - k;
/*
*  Accumulate A( IA+k:IA+N-1, JA+k:JA+k+kbb-1 )
*/
               PB_CGatherV( TYPE, ALLOCATE, &Dir, ktmp, kbb, A, IA+k, JA+k,
                            DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Scale A( IA+k:IA+N-1, JA+k:JA+k+kbb-1 ) by ALPHA
*/
               if( conjg )
                  PB_Cplacnjg( TYPE, ktmp, kbb, ALPHA, Aptr, 0, 0, DBUFA );
               else
                  PB_Cplascal( TYPE, ALL, NOCONJG, ktmp, kbb, ALPHA, Aptr,
                               0, 0, DBUFA );
/*
*  Zero upper triangle of A( IA+k:IA+N-1, JA+k:JA+k+kbb-1 )
*/
               if( kbb > 1 )
                  PB_Cplapad( TYPE, UPPER, NOCONJG, kbb-1, kbb-1, zero, zero,
                              Aptr, 0, 1, DBUFA );
/*
*  C( IC+k:IC+k+kbb-1, JC+k:JC+N-1 ) += A( IA+k:IA+N-1, JA+k:JA+k+kbb-1 )'
*/
               PB_CScatterV( TYPE, &Dir, ktmp, kbb, Aptr, 0, 0, DBUFA, COLUMN,
                             one, C, IC+k, JC+k, DESCC, ROW );

               if( Afr ) free( Aptr );
            }
         }
      }
      else
      {
         if( M >= N )
         {
            for( k = kstart; k != kend; k += kstep )
            {
               kbb = N - k; kbb = MIN( kbb, kb ); ktmp = M - k;
/*
*  Accumulate A( IA+k:IA+k+kbb-1, JA+k:JA+M-1 )
*/
               PB_CGatherV( TYPE, ALLOCATE, &Dir, kbb, ktmp, A, IA+k, JA+k,
                            DESCA, ROW, &Aptr, DBUFA, &Afr );
/*
*  Scale A( IA+k:IA+k+kbb-1, JA+k:JA+M-1 ) by ALPHA
*/
               if( conjg )
                  PB_Cplacnjg( TYPE, kbb, ktmp, ALPHA, Aptr, 0, 0, DBUFA );
               else
                  PB_Cplascal( TYPE, ALL, NOCONJG, kbb, ktmp, ALPHA, Aptr,
                               0, 0, DBUFA );
/*
*  Zero lower triangle of A( IA+k:IA+k+kbb-1, JA+k:JA+k+kbb-1 )
*/
               if( kbb > 1 )
                  PB_Cplapad( TYPE, LOWER, NOCONJG, kbb-1, kbb-1, zero, zero,
                              Aptr, 1, 0, DBUFA );
/*
*  C( IC:IC+k+kbb-1, JC+k:JC+k+kbb-1 ) += A( IA+k:IA+k+kbb-1, JA+k:JA+M-1 )'
*/
               PB_CScatterV( TYPE, &Dir, kbb, ktmp, Aptr, 0, 0, DBUFA, ROW,
                             one, C, IC+k, JC+k, DESCC, COLUMN );

               if( Afr ) free( Aptr );
            }
         }
         else
         {
            for( k = kstart; k != kend; k += kstep )
            {
               kbb = M - k; kbb = MIN( kbb, kb ); ktmp = k + kbb;
/*
*  Accumulate A( IA:IA+k+kbb-1, JA+k:JA+k+kbb-1 )
*/
               PB_CGatherV( TYPE, ALLOCATE, &Dir, ktmp, kbb, A, IA, JA+k,
                            DESCA, COLUMN, &Aptr, DBUFA, &Afr );
/*
*  Scale A( IA:IA+k+kbb-1, JA+k:JA+k+kbb-1 ) by ALPHA
*/
               if( conjg )
                  PB_Cplacnjg( TYPE, ktmp, kbb, ALPHA, Aptr, 0, 0, DBUFA );
               else
                  PB_Cplascal( TYPE, ALL, NOCONJG, ktmp, kbb, ALPHA, Aptr,
                               0, 0, DBUFA );
/*
*  Zero lower triangle of A( IA+k:IA+k+kbb-1, JA+k:JA:JA+k+kbb-1 )
*/
               if( kbb > 1 )
                  PB_Cplapad( TYPE, LOWER, NOCONJG, kbb-1, kbb-1, zero, zero,
                              Aptr, k+1, 0, DBUFA );
/*
*  C( IC+k:IC+k+kbb-1, JC:JC+k+kbb-1 ) += A( IA:IA+k+kbb-1, JA+k:JA+k+kbb-1 )'
*/
               PB_CScatterV( TYPE, &Dir, ktmp, kbb, Aptr, 0, 0, DBUFA, COLUMN,
                             one, C, IC+k, JC, DESCC, ROW );

               if( Afr ) free( Aptr );
            }
         }
      }
   }
/*
*  End of PB_Cptradd
*/
}
