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
void PB_Cplaprnt( PBTYP_T * TYPE, Int M, Int N,
                  char * A, Int IA, Int JA, Int * DESCA,
                  Int IRPRNT, Int ICPRNT, char * CMATNM )
#else
void PB_Cplaprnt( TYPE, M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT, CMATNM )
/*
*  .. Scalar Arguments ..
*/
   Int            IA, ICPRNT, IRPRNT, JA, M, N;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA;
   char           * A, * CMATNM;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cplaprnt  prints to the standard output the submatrix sub( A ) de-
*  noting A(IA:IA+M-1,JA:JA+N-1).  The local pieces of sub( A ) are sent
*  and printed by the process of coordinates (IRPRNT, ICPRNT).
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
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( A ). M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( A ). N must be at least zero.
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
*  IRPRNT  (global input) INTEGER
*          On entry, IRPRNT specifies the row index of the printing pro-
*          cess.
*
*  ICPRNT  (global input) INTEGER
*          On entry, ICPRNT specifies the  column  index of the printing
*          process.
*
*  CMATNM  (global input) pointer to CHAR
*          On entry, CMATNM is the name of the matrix to be printed.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   Int            mycol, myrow, npcol, nprow, pcol, prow;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( DESCA[CTXT_], &nprow, &npcol, &myrow, &mycol );
/*
*  When sub( A ) is replicated, each copy is printed for debugging purposes.
*/
   if( DESCA[ RSRC_ ] >= 0 )
   {
/*
*  sub( A ) is distributed onto the process rows of the grid
*/
      if( DESCA[ CSRC_ ] >= 0 )
      {
/*
*  sub( A ) is distributed onto the process columns of the grid
*/
         PB_Cplaprn2( TYPE, M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT, CMATNM,
                      DESCA[ RSRC_ ], DESCA[ CSRC_ ] );
      }
      else
      {
/*
*  sub( A ) is replicated in every process column of the grid
*/
         for( pcol = 0; pcol < npcol; pcol++ )
         {
            if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) )
               (void) fprintf( stdout,
               "Colum-replicated array -- copy in process column: %d\n", pcol );
            PB_Cplaprn2( TYPE, M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT, CMATNM,
                         DESCA[ RSRC_ ], pcol );
         }
      }
   }
   else
   {
/*
*  sub( A ) is replicated in every process row of the grid
*/
      if( DESCA[ CSRC_ ] >= 0 )
      {
/*
*  sub( A ) is distributed onto the process columns of the grid
*/
         for( prow = 0; prow < nprow; prow++ )
         {
            if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) )
               (void) fprintf( stdout,
               "Row-replicated array -- copy in process row: %d\n", prow );
            PB_Cplaprn2( TYPE, M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT, CMATNM,
                         prow, DESCA[ CSRC_ ] );
         }
      }
      else
      {
/*
*  sub( A ) is replicated in every process column of the grid
*/
         for( prow = 0; prow < nprow; prow++ )
         {
            for( pcol = 0; pcol < npcol; pcol++ )
            {
               if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) )
                  (void) fprintf( stdout,
               "Replicated array -- copy in process (%d,%d)\n", prow, pcol );
               PB_Cplaprn2( TYPE, M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT,
                            CMATNM, prow, pcol );
            }
         }
      }
   }
/*
*  End of PB_Cplaprnt
*/
}

#ifdef __STDC__
void PB_Cplaprn2( PBTYP_T * TYPE, Int M, Int N, char * A, Int IA,
                  Int JA, Int * DESCA, Int IRPRNT, Int ICPRNT,
                  char * CMATNM, Int PROW, Int PCOL )
#else
void PB_Cplaprn2( TYPE, M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT, CMATNM,
                  PROW, PCOL )
/*
*  .. Scalar Arguments ..
*/
   Int            IA, ICPRNT, IRPRNT, JA, M, N, PCOL, PROW;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA;
   char           * A, * CMATNM;
#endif
{
/*
*  .. Local Scalars ..
*/
   char           type;
   Int            Acol, Aii, AisColRep, AisRowRep, Ajj, Ald, Arow, ctxt, h, i,
                  ib, icurcol, icurrow, ii, in, j, jb, jj, jn, ldw, mycol,
                  myrow, npcol, nprow, size, usiz;
/*
*  .. Local Arrays ..
*/
   char           * buf = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Retrieve sub( A )'s local information: Aii, Ajj, Arow, Acol ...
*/
   Ald  = DESCA[LLD_];
   PB_Cinfog2l( IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aii, &Ajj, &Arow,
                &Acol );
/*
*  Save the local first index of each row and column sub( A )
*/
   ii   = Aii;
   jj   = Ajj;
/*
*  When sub( A ) is row-replicated, print the copy in process row PROW.
*  Otherwise, print the distributed matrix rows starting in process row Arow.
*/
   if( Arow < 0 ) { AisRowRep = 1; icurrow = Arow = PROW; }
   else           { AisRowRep = 0; icurrow = Arow;        }
/*
*  When sub( A ) is column-replicated, print the copy in process column PCOL.
*  Otherwise, print the distributed matrix columns starting in process column
*  Acol.
*/
   if( Acol < 0 ) { AisColRep = 1; icurcol = Acol = PCOL; }
   else           { AisColRep = 0; icurcol = Acol;        }

   type = TYPE->type; usiz = TYPE->usiz;  size = TYPE->size;
/*
*  Allocate buffer in printing process
*/
   ldw  = MAX( DESCA[ IMB_ ], DESCA[ MB_ ] );
   if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) )
      buf = PB_Cmalloc( ldw * size );
/*
*  Handle the first block of column separately
*/
   jb = PB_Cfirstnb( N, JA, DESCA[INB_], DESCA[NB_] );
   jn = JA + jb - 1;

   for( h = 0; h < jb; h++ )
   {
      ib = PB_Cfirstnb( M, IA, DESCA[IMB_], DESCA[MB_] );
      in = IA + ib - 1;

      if( ( icurrow == IRPRNT ) && ( icurcol == ICPRNT ) )
      {
         if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) )
         {
            PB_Cprnt( type, size, usiz, ib, Mptr( A, ii, jj+h, Ald, size ),
                      IA+1, JA+h+1, CMATNM );
         }
      }
      else
      {
         if( ( myrow == icurrow ) && ( mycol == icurcol ) )
         {
            TYPE->Cgesd2d( ctxt, ib, 1, Mptr( A, ii, jj+h, Ald, size ), Ald,
                           IRPRNT, ICPRNT );
         }
         else if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) )
         {
            TYPE->Cgerv2d( ctxt, ib, 1, buf, ldw, icurrow, icurcol );
            PB_Cprnt( type, size, usiz, ib, buf, IA+1, JA+h+1, CMATNM );
         }
      }
/*
*  Go to next block of rows
*/
      if( myrow == icurrow ) ii += ib;
      if( !( AisRowRep ) ) icurrow = MModAdd1( icurrow, nprow );

      Cblacs_barrier( ctxt, ALL );
/*
*  Loop over remaining block of rows
*/
      for( i = in+1; i <= IA+M-1; i += DESCA[MB_] )
      {
         ib = MIN( DESCA[MB_], IA+M-i );
         if( ( icurrow == IRPRNT ) && ( icurcol == ICPRNT ) )
         {
            if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) )
            {
               PB_Cprnt( type, size, usiz, ib, Mptr( A, ii, jj+h, Ald, size ),
                         i+1, JA+h+1, CMATNM );
            }
         }
         else
         {
            if( ( myrow == icurrow ) && ( mycol == icurcol ) )
            {
               TYPE->Cgesd2d( ctxt, ib, 1, Mptr( A, ii, jj+h, Ald, size ), Ald,
                              IRPRNT, ICPRNT );
            }
            else if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) )
            {
               TYPE->Cgerv2d( ctxt, ib, 1, buf, ldw, icurrow, icurcol );
               PB_Cprnt( type, size, usiz, ib, buf, i+1, JA+h+1, CMATNM);
            }
         }
/*
*  Go to next block of rows
*/
         if( myrow == icurrow ) ii += ib;
         if( !( AisRowRep ) ) icurrow = MModAdd1( icurrow, nprow );

         Cblacs_barrier( ctxt, ALL );
      }
/*
*  Restart at the first row to be printed
*/
      ii = Aii;
      icurrow = Arow;
   }
/*
*  Go to next block of columns
*/
   if( mycol == icurcol ) jj += jb;
   if( !( AisColRep ) ) icurcol = MModAdd1( icurcol, npcol );

   Cblacs_barrier( ctxt, ALL );
/*
*  Loop over remaining column blocks
*/
   for( j = jn+1; j <= JA+N-1; j += DESCA[NB_] )
   {
      jb = MIN( DESCA[NB_], JA+N-j );
      for( h = 0; h < jb; h++ )
      {
         ib = PB_Cfirstnb( M, IA, DESCA[IMB_], DESCA[MB_] );
         in = IA + ib - 1;

         if( ( icurrow == IRPRNT ) && ( icurcol == ICPRNT ) )
         {
            if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) )
            {
               PB_Cprnt( type, size, usiz, ib, Mptr( A, ii, jj+h, Ald, size ),
                         IA+1, j+h+1, CMATNM );
            }
         }
         else
         {
            if( ( myrow == icurrow ) && ( mycol == icurcol ) )
            {
               TYPE->Cgesd2d( ctxt, ib, 1, Mptr( A, ii, jj+h, Ald, size ), Ald,
                              IRPRNT, ICPRNT );
            }
            else if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) )
            {
               TYPE->Cgerv2d( ctxt, ib, 1, buf, ldw, icurrow, icurcol );
               PB_Cprnt( type, size, usiz, ib, buf, IA+1, j+h+1, CMATNM );
            }
         }
/*
*  Go to next block of rows
*/
         if( myrow == icurrow ) ii += ib;
         if( !( AisRowRep ) ) icurrow = MModAdd1( icurrow, nprow );

         Cblacs_barrier( ctxt, ALL );
/*
*  Loop over remaining block of rows
*/
         for( i = in+1; i <= IA+M-1; i += DESCA[MB_] )
         {
            ib = MIN( DESCA[MB_], IA+M-i );
            if( ( icurrow == IRPRNT ) && ( icurcol == ICPRNT ) )
            {
               if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) )
               {
                  PB_Cprnt( type, size, usiz, ib, Mptr( A, ii, jj+h, Ald,
                            size ), i+1, j+h+1, CMATNM );
               }
            }
            else
            {
               if( ( myrow == icurrow ) && ( mycol == icurcol ) )
               {
                  TYPE->Cgesd2d( ctxt, ib, 1, Mptr( A, ii, jj+h, Ald, size ),
                                 Ald, IRPRNT, ICPRNT );
               }
               else if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) )
               {
                  TYPE->Cgerv2d( ctxt, ib, 1, buf, ldw, icurrow, icurcol );
                  PB_Cprnt( type, size, usiz, ib, buf, i+1, j+h+1, CMATNM );
               }
            }
/*
*  Go to next block of rows
*/
            if( myrow == icurrow ) ii += ib;
            if( !( AisRowRep ) ) icurrow = MModAdd1( icurrow, nprow );

            Cblacs_barrier( ctxt, ALL );
         }
/*
*  Restart at the first row to be printed
*/
         ii      = Aii;
         icurrow = Arow;
      }
/*
*  Go to next block of columns
*/
      if( mycol == icurcol ) jj += jb;
      if( !( AisColRep ) ) icurcol = MModAdd1( icurcol, npcol );

      Cblacs_barrier( ctxt, ALL );
   }

   if( ( myrow == IRPRNT ) && ( mycol == ICPRNT ) && ( buf ) ) free( buf );
/*
*  End of PB_Cplaprn2
*/
}

#ifdef __STDC__
void PB_Cprnt( char TYPE, Int SIZE, Int USIZ, Int N, char * A, Int IA,
               Int JA, char * CMATNM )
#else
void PB_Cprnt( TYPE, SIZE, USIZ, N, A, IA, JA, CMATNM )
/*
*  .. Scalar Arguments ..
*/
   Int            IA, JA, N, SIZE, TYPE, USIZ;
/*
*  .. Array Arguments ..
*/
   char           * A, * CMATNM;
#endif
{
/*
*  .. Local Scalars ..
*/
   Int            k;
/* ..
*  .. Executable Statements ..
*
*/
   if( TYPE == INT )
      for( k = 0; k < N; k++ )
         (void) fprintf( stdout, "%s(%6d,%6d)=%8d\n",      CMATNM, IA+k, JA,
                         *((Int *)(&A[k*SIZE])) );
   else if( TYPE == SREAL )
      for( k = 0; k < N; k++ )
         (void) fprintf( stdout, "%s(%6d,%6d)=%16.8f\n",   CMATNM, IA+k, JA,
                         *((float *)(&A[k*SIZE])) );
   else if( TYPE == DREAL )
      for( k = 0; k < N; k++ )
         (void) fprintf( stdout, "%s(%6d,%6d)=%30.18f\n",  CMATNM, IA+k, JA,
                         *((double *)(&A[k*SIZE])) );
   else if( TYPE == SCPLX )
      for( k = 0; k < N; k++ )
         (void) fprintf( stdout, "%s(%6d,%6d)=%16.8f+i*(%16.8f)\n",   CMATNM,
                         IA+k, JA, *((float *)(&A[k*SIZE])),
                         *((float *)(&A[k*SIZE+USIZ])) );
   else if( TYPE == DCPLX )
      for( k = 0; k < N; k++ )
         (void) fprintf( stdout, "%s(%6d,%6d)=%30.18f+i*(%30.18f)\n", CMATNM,
                         IA+k, JA, *((double *)(&A[k*SIZE])),
                         *((double *)(&A[k*SIZE+USIZ])) );
/*
*  End of PB_Cprnt
*/
}
