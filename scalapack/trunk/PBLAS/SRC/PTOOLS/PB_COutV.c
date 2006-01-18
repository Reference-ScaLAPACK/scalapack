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
void PB_COutV( PBTYP_T * TYPE, char * ROWCOL, char * ZEROIT, int M,
               int N, int * DESCA, int K, char * * YAPTR, int * DYA,
               int * YAFREE, int * YASUM )
#else
void PB_COutV( TYPE, ROWCOL, ZEROIT, M, N, DESCA, K, YAPTR, DYA,
               YAFREE, YASUM )
/*
*  .. Scalar Arguments ..
*/
   char           * ROWCOL, * ZEROIT;
   int            * YAFREE, K, M, N, * YASUM;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DYA;
   char           * * YAPTR;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_COutV  returns a pointer to an array that contains a one-dimensio-
*  nal ouput zero subvector which is replicated over the rows or columns
*  of a submatrix described by DESCA. On return, the subvector is speci-
*  fied  by  a  pointer  to some data, a descriptor array describing its
*  layout, a  logical  value  indicating if this local piece of data has
*  been  dynamically  allocated by this function, a logical value speci-
*  fying  if  sum  reduction  should occur. This routine is specifically
*  designed for traditional Level 2 and 3 PBLAS operations using an out-
*  put only vector such as PxTRMV, or PxTRMM.
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
*  ROWCOL  (global input) pointer to CHAR
*          On entry, ROWCOL  specifies  if  this routine should return a
*          row or column subvector replicated over the underlying subma-
*          trix as follows:
*             = 'R' or 'r':                 A row subvector is returned,
*             = 'C' or 'c':              A column subvector is returned.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies the number of rows of  the underlying
*          submatrix described by DESCA. M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N specifies the number of columns of the underlying
*          submatrix described by DESCA. N must be at least zero.
*
*  DESCA   (global and local input) INTEGER array
*          On entry, DESCA  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix A.
*
*  K       (global input) INTEGER
*          On entry,  K  specifies the length of the non-distributed di-
*          mension of the subvector sub( Y ). K must be at least zero.
*
*  YAPTR   (local output) pointer to pointer to CHAR
*          On exit, * YAPTR  is an array containing the same data as the
*          subvector  sub( Y )  which is replicated over the rows or co-
*          lumns of the  underlying  matrix  as  specified by ROWCOL and
*          DESCA.
*
*  DYA     (global and local output) INTEGER array
*          On exit, DYA is a descriptor array of dimension DLEN_ descri-
*          bing the data layout of the data pointed to by * YAPTR.
*
*  YAFREE  (local output) INTEGER
*          On exit, YAFREE  specifies  if  it  was possible to reuse the
*          subvector sub( Y ),  i.e., if some dynamic memory was alloca-
*          ted for the data pointed to by * YAPTR or not. When YAFREE is
*          zero, no  dynamic memory was allocated. Otherwise, some dyna-
*          mic memory  was  allocated by this function that one MUST re-
*          lease as soon as possible.
*
*  YASUM   (global output) INTEGER
*          On exit, YASUM  specifies if a global sum reduction should be
*          performed to obtain the correct sub( Y ). When YASUM is zero,
*          no reduction  is  to be performed, otherwise reduction should
*          occur.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   int            Acol, Aimb, Ainb, Amb, Amp, Anb, Anq, Arow, Yld, ctxt,
                  izero=0, nprow, myrow, npcol, mycol;
   char           * zero;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Initialize the output parameters to a default value
*/
   *YAFREE = 0;
   *YASUM  = 0;
   *YAPTR  = NULL;
/*
*  Quick return if possible
*/
   if( ( M <= 0 ) || ( N <= 0 ) || ( K <= 0 ) )
   {
      if( Mupcase( ROWCOL[0] ) == CROW )
      {
         PB_Cdescset( DYA, K, N, 1, DESCA[INB_], 1, DESCA[NB_], DESCA[RSRC_],
                      DESCA[CSRC_], DESCA[CTXT_], 1 );
      }
      else
      {
         PB_Cdescset( DYA, M, K, DESCA[IMB_], 1, DESCA[MB_], 1, DESCA[RSRC_],
                      DESCA[CSRC_], DESCA[CTXT_], DESCA[LLD_] );
      }
      return;
   }
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   Arow = DESCA[RSRC_]; Acol = DESCA[CSRC_];

   if( Mupcase( ROWCOL[0] ) == CROW )
   {
/*
*  Want a row vector
*/
      Ainb = DESCA[INB_]; Anb = DESCA[NB_];
      Anq  = PB_Cnumroc( N, 0, Ainb, Anb, mycol, Acol, npcol );
      Yld  = MAX( 1, K );

      if( ( Arow <  0 ) || ( nprow == 1 ) ||
          ( PB_Cspan( M, 0, DESCA[IMB_], DESCA[MB_], Arow, nprow ) ) )
      {
/*
*  A spans all process rows. Y should be reduced iff A is not replicated and
*  there is more than just one process row in the process grid.
*/
         *YASUM = ( ( Arow >= 0 ) && ( nprow > 1 ) );
/*
*  Allocate the space for Y in the processes owning at least one column of A,
*  and initialize it to zero if requested.
*/
         if( Anq > 0 )
         {
            *YAPTR  = PB_Cmalloc( K * Anq * TYPE->size );
            *YAFREE = 1;
            if( Mupcase( ZEROIT[0] ) == CINIT )
            {
               zero = TYPE->zero;
               TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &K, &Anq,
                             &izero, zero, zero, *YAPTR, &Yld );
            }
         }
/*
*  Describe the newly created operand
*/
         PB_Cdescset( DYA, K, N, K, Ainb, 1, Anb, -1, Acol, ctxt, Yld );
      }
      else
      {
/*
*  A spans only one process row. There is no need to reduce Y or even to
*  allocate some space for it outside this process row.
*/
         *YASUM = 0;
         if( ( myrow == Arow ) && ( Anq > 0 ) )
         {
            *YAPTR  = PB_Cmalloc( K * Anq * TYPE->size );
            *YAFREE = 1;
            if( Mupcase( ZEROIT[0] ) == CINIT )
            {
               zero = TYPE->zero;
               TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &K, &Anq,
                             &izero, zero, zero, *YAPTR, &Yld );
            }
         }
/*
*  Describe the newly created operand
*/
         PB_Cdescset( DYA, K, N, K, Ainb, 1, Anb, Arow, Acol, ctxt, Yld );
      }
   }
   else
   {
/*
*  Want a column vector
*/
      Aimb = DESCA[ IMB_  ]; Amb  = DESCA[ MB_   ];
      Amp  = PB_Cnumroc( M, 0, Aimb, Amb, myrow, Arow, nprow );
      Yld  = MAX( 1, Amp );

      if( ( Acol <  0 ) || ( npcol == 1 ) ||
          ( PB_Cspan( N, 0, DESCA[INB_], DESCA[NB_], Acol, npcol ) ) )
      {
/*
*  A spans all process columns. Y should be reduced iff A is not replicated and
*  there is more than just one process column in the process grid.
*/
         *YASUM = ( ( Acol >= 0 ) && ( npcol > 1 ) );
/*
*  Allocate the space for Y in the processes owning at least one row of A, and
*  initialize it to zero if requested.
*/
         if( Amp > 0 )
         {
            *YAPTR  = PB_Cmalloc( Amp * K * TYPE->size );
            *YAFREE = 1;
            if( Mupcase( ZEROIT[0] ) == CINIT )
            {
               zero = TYPE->zero;
               TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Amp, &K,
                             &izero, zero, zero, *YAPTR, &Yld );
            }
         }
/*
*  Describe the newly created operand
*/
         PB_Cdescset( DYA, M, K, Aimb, K, Amb, 1, Arow, -1, ctxt, Yld );
      }
      else
      {
/*
*  A spans only one process column. There is no need to reduce Y or even to
*  allocate some space for it outside this process column.
*/
         *YASUM = 0;
         if( ( mycol == Acol ) && ( Amp > 0 ) )
         {
            *YAPTR  = PB_Cmalloc( Amp * K * TYPE->size );
            *YAFREE = 1;
            if( Mupcase( ZEROIT[0] ) == CINIT )
            {
               zero = TYPE->zero;
               TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Amp, &K,
                             &izero, zero, zero, *YAPTR, &Yld );
            }
         }
/*
*  Describe the newly created operand
*/
         PB_Cdescset( DYA, M, K, Aimb, K, Amb, 1, Arow, Acol, ctxt, Yld );
      }
   }
/*
*  End of PB_COutV
*/
}
