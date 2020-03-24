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
void PB_Cplapad( PBTYP_T * TYPE, char * UPLO, char * CONJUG, Int M,
                 Int N, char * ALPHA, char * BETA, char * A, Int IA,
                 Int JA, Int * DESCA )
#else
void PB_Cplapad( TYPE, UPLO, CONJUG, M, N, ALPHA, BETA, A, IA, JA,
                 DESCA )
/*
*  .. Scalar Arguments ..
*/
   char           * CONJUG, * UPLO;
   Int            IA, JA, M, N;
   char           * ALPHA, * BETA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA;
   char           * A;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cplapad  initializes  an m by n  submatrix  A(IA:IA+M-1,JA:JA+N-1)
*  denoted by  sub( A )  to beta on the diagonal or  zeros the imaginary
*  part of those diagonals and set the offdiagonals to alpha.
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
*  UPLO    (global input) pointer to CHAR
*          On entry, UPLO specifies the part  of  the submatrix sub( A )
*          to be set:
*             = 'L' or 'l': Lower triangular  part is set; the  strictly
*                      upper triangular part of sub( A ) is not changed;
*             = 'U' or 'u': Upper triangular  part is set; the  strictly
*                      lower triangular part of sub( A ) is not changed;
*             Otherwise:  All of the matrix sub( A ) is set.
*
*  CONJUG  (global input) pointer to CHAR
*          On entry,  CONJUG specifies what should be done to the diago-
*          nals as follows. When UPLO is 'L', 'l', 'U' or 'u' and CONJUG
*          is 'Z' or 'z', the imaginary part of the diagonals is set  to
*          zero. Otherwise, the diagonals are set to beta.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( A ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( A ). N must be at least zero.
*
*  ALPHA   (global input) pointer to CHAR
*          On entry,  ALPHA  specifies the scalar alpha, i.e., the cons-
*          tant to which the offdiagonal elements are to be set.
*
*  BETA    (global input) pointer to CHAR
*          On entry, BETA  specifies the scalar beta, i.e., the constant
*          to which the diagonal elements are to be set. BETA is not re-
*          ferenced when UPLO is 'L', 'l', 'U' or 'u' and CONJUG is 'Z'.
*
*  A       (local input/local output) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
*          the local entries of the matrix  A  to be  set.  On exit, the
*          leading m by n submatrix sub( A ) is set as follows:
*
*          UPLO = 'L' or 'l', A(IA+i-1,JA+j-1)=ALPHA, j+1<=i<=M, 1<=j<=N
*          UPLO = 'U' or 'u', A(IA+i-1,JA+j-1)=ALPHA, 1<=i<=j-1, 1<=j<=N
*          otherwise,         A(IA+i-1,JA+j-1)=ALPHA, 1<=i<=M,   1<=j<=N
*                                                      and IA+i.NE.JA+j,
*          and, for all UPLO,  A(IA+i-1,JA+i-1) = BETA,  1<=i<=min(M,N).
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
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           type;
   Int            Acol, Aii, Aimb1, Ainb1, Ajj, Akp, Akq, Ald, Amb, Amp, Amp0,
                  Anb, Anq, Anq0, Arow, ctxt, izero=0, k, kb, ktmp, mn, mycol,
                  myrow, nb, npcol, nprow, size;
   TZPAD_T        pad;
/*
*  .. Local Arrays ..
*/
   Int            Ad0[DLEN_];
   char           * Aptr = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Quick return if possible
*/
   if( ( M <= 0 ) || ( N <= 0 ) ) return;
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Compute descriptor Ad0 for sub( A )
*/
   PB_Cdescribe( M, N, IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aii, &Ajj,
                 &Ald, &Aimb1, &Ainb1, &Amb, &Anb, &Arow, &Acol, Ad0 );
/*
*  Quick return if I don't own any of sub( A ).
*/
   Amp = PB_Cnumroc( M, 0, Aimb1, Amb, myrow, Arow, nprow );
   Anq = PB_Cnumroc( N, 0, Ainb1, Anb, mycol, Acol, npcol );
   if( ( Amp <= 0 ) || ( Anq <= 0 ) ) return;

   size = TYPE->size; type = TYPE->type; pad  = TYPE->Ftzpad;
   Aptr = Mptr( A, Aii, Ajj, Ald, size );
/*
*  When the entire sub( A ) needs to be padded and alpha is equal to beta, or
*  sub( A ) is replicated in all processes, just call the local routine.
*/
   if( type == SREAL )
   {
      if( ( ( Mupcase( UPLO[0] ) == CALL ) &&
            ( ((float*)(ALPHA))[REAL_PART] ==
              ((float*)(BETA ))[REAL_PART] ) ) ||
          ( ( ( Arow < 0 ) || ( nprow == 1 ) ) &&
            ( ( Acol < 0 ) || ( npcol == 1 ) ) ) )
      {
         pad( C2F_CHAR( UPLO ), C2F_CHAR( CONJUG ), &Amp, &Anq, &izero, ALPHA,
              BETA, Aptr, &Ald );
         return;
      }
   }
   else if( type == DREAL )
   {
      if( ( ( Mupcase( UPLO[0] ) == CALL ) &&
            ( ((double*)(ALPHA))[REAL_PART] ==
              ((double*)(BETA ))[REAL_PART] ) ) ||
          ( ( ( Arow < 0 ) || ( nprow == 1 ) ) &&
            ( ( Acol < 0 ) || ( npcol == 1 ) ) ) )
      {
         pad( C2F_CHAR( UPLO ), C2F_CHAR( CONJUG ), &Amp, &Anq, &izero, ALPHA,
              BETA, Aptr, &Ald );
         return;
      }
   }
   else if( type == SCPLX )
   {
      if( ( ( Mupcase( UPLO[0] ) == CALL ) &&
            ( ((float*)(ALPHA))[REAL_PART] ==
              ((float*)(BETA ))[REAL_PART] ) &&
            ( ((float*)(ALPHA))[IMAG_PART] ==
              ((float*)(BETA ))[IMAG_PART] ) ) ||
          ( ( ( Arow < 0 ) || ( nprow == 1 ) ) &&
            ( ( Acol < 0 ) || ( npcol == 1 ) ) ) )
      {
         pad( C2F_CHAR( UPLO ), C2F_CHAR( CONJUG ), &Amp, &Anq, &izero, ALPHA,
              BETA, Aptr, &Ald );
         return;
      }
   }
   else if( type == DCPLX )
   {
      if( ( ( Mupcase( UPLO[0] ) == CALL ) &&
            ( ((double*)(ALPHA))[REAL_PART] ==
              ((double*)(BETA ))[REAL_PART] ) &&
            ( ((double*)(ALPHA))[IMAG_PART] ==
              ((double*)(BETA ))[IMAG_PART] ) ) ||
          ( ( ( Arow < 0 ) || ( nprow == 1 ) ) &&
            ( ( Acol < 0 ) || ( npcol == 1 ) ) ) )
      {
         pad( C2F_CHAR( UPLO ), C2F_CHAR( CONJUG ), &Amp, &Anq, &izero, ALPHA,
              BETA, Aptr, &Ald );
         return;
      }
   }
/*
*  Computational partitioning size is computed as the product of the logical
*  value returned by pilaenv_ and two times the least common multiple of nprow
*  and npcol.
*/
   nb = 2 * pilaenv_( &ctxt, C2F_CHAR( &type ) ) *
        PB_Clcm( ( Arow >= 0 ? nprow : 1 ), ( Acol >= 0 ? npcol : 1 ) );

   mn = MIN( M, N );

   if( Mupcase( UPLO[0] ) == CLOWER )
   {
/*
*  Lower triangle of sub( A ): proceed by block of columns. For each block of
*  columns, operate on the logical diagonal block first and then the remaining
*  rows of that block of columns.
*/
      for( k = 0; k < mn; k += nb )
      {
         kb   = mn - k; ktmp = k + ( kb = MIN( kb, nb ) );
         PB_Cplapd2( TYPE, UPLO, CONJUG, kb, kb, ALPHA, BETA, Aptr, k, k, Ad0 );
         Akp  = PB_Cnumroc( ktmp, 0, Aimb1, Amb, myrow, Arow, nprow );
         Akq  = PB_Cnumroc( k,    0, Ainb1, Anb, mycol, Acol, npcol );
         Anq0 = PB_Cnumroc( kb,   k, Ainb1, Anb, mycol, Acol, npcol );
         if( ( Amp0 = Amp - Akp ) > 0 )
            pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Amp0, &Anq0, &izero,
                 ALPHA, ALPHA, Mptr( Aptr, Akp, Akq, Ald, size ), &Ald );
      }
   }
   else if( Mupcase( UPLO[0] ) == CUPPER )
   {
/*
*  Upper triangle of sub( A ): proceed by block of columns. For each block of
*  columns, operate on the trailing rows and then the logical diagonal block
*  of that block of columns. When M < N, the last columns of sub( A ) are
*  handled together.
*/
      for( k = 0; k < mn; k += nb )
      {
         kb   = mn - k; kb = MIN( kb, nb );
         Akp  = PB_Cnumroc( k,  0, Aimb1, Amb, myrow, Arow, nprow );
         Akq  = PB_Cnumroc( k,  0, Ainb1, Anb, mycol, Acol, npcol );
         Anq0 = PB_Cnumroc( kb, k, Ainb1, Anb, mycol, Acol, npcol );
         if( Akp > 0 )
            pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Akp, &Anq0, &izero,
                 ALPHA, ALPHA, Mptr( Aptr, 0, Akq, Ald, size ), &Ald );
         PB_Cplapd2( TYPE, UPLO, CONJUG, kb, kb, ALPHA, BETA, Aptr, k, k, Ad0 );
      }
      if( ( Anq -= ( Akq += Anq0 ) ) > 0 )
         pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Amp, &Anq, &izero, ALPHA,
              ALPHA, Mptr( Aptr, 0, Akq, Ald, size ), &Ald );
   }
   else
   {
/*
*  All of sub( A ): proceed by block of columns. For each block of columns,
*  operate on the trailing rows, then the logical diagonal block, and finally
*  the remaining rows of that block of columns. When M < N, the last columns
*  of sub( A ) are handled together.
*/
      for( k = 0; k < mn; k += nb )
      {
         kb   = mn - k; kb = MIN( kb, nb );
         Akp  = PB_Cnumroc( k,  0, Aimb1, Amb, myrow, Arow, nprow );
         Akq  = PB_Cnumroc( k,  0, Ainb1, Anb, mycol, Acol, npcol );
         Anq0 = PB_Cnumroc( kb, k, Ainb1, Anb, mycol, Acol, npcol );
         if( Akp > 0 )
            pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Akp, &Anq0, &izero,
                 ALPHA, ALPHA, Mptr( Aptr, 0, Akq, Ald, size ), &Ald );
         PB_Cplapd2( TYPE, UPLO, NOCONJG, kb, kb, ALPHA, BETA, Aptr, k, k,
                     Ad0 );
         Akp = PB_Cnumroc( k+kb, 0, Aimb1, Amb, myrow, Arow, nprow );
         if( ( Amp0 = Amp - Akp ) > 0 )
            pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Amp0, &Anq0, &izero,
                 ALPHA, ALPHA, Mptr( Aptr, Akp, Akq, Ald, size ), &Ald );
      }
      if( ( Anq -= ( Akq += Anq0 ) ) > 0 )
         pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Amp, &Anq, &izero, ALPHA,
              ALPHA, Mptr( Aptr, 0, Akq, Ald, size ), &Ald );
   }
/*
*  End of PB_Cplapad
*/
}
