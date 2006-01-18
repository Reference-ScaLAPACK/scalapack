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
void PB_CInV( PBTYP_T * TYPE, char * CONJUG, char * ROWCOL, int M,
              int N, int * DESCA, int K, char * X, int IX, int JX,
              int * DESCX, char * XROC, char * * XAPTR, int * DXA,
              int * XAFREE )
#else
void PB_CInV( TYPE, CONJUG, ROWCOL, M, N, DESCA, K, X, IX, JX, DESCX,
              XROC, XAPTR, DXA, XAFREE )
/*
*  .. Scalar Arguments ..
*/
   char           * CONJUG, * ROWCOL, * XROC;
   int            * XAFREE, IX, JX, K, M, N;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   int            * DESCA, * DESCX, * DXA;
   char           * X, * * XAPTR;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CInV returns a pointer to an array that contains a one-dimensional
*  input only subvector which is replicated over the rows or columns  of
*  a submatrix described by DESCA. A subvector is specified on  input to
*  this routine that is reused whenever possible. On return, the subvec-
*  tor is specified by a pointer to some data,  a  descriptor array des-
*  cribing its layout and a logical value indicating if this local piece
*  of data has been dynamically allocated by this function. This routine
*  is  specifically  designed  for traditional Level 2 like PBLAS opera-
*  tions using an input only vector such as PxGER, PxSYR ...
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
*  CONJUG  (global input) pointer to CHAR
*          On  entry,  CONJUG  specifies  if  this routine should return
*          the conjugate subvector as follows:
*             = 'N' or 'n':           The initial subvector is returned,
*             = 'Z' or 'z':         The conjugate subvector is returned.
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
*          mension of the subvector sub( X ). K must be at least zero.
*
*  X       (local input) pointer to CHAR
*          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
*          is   at  least  MAX( 1, Lr( K, IX ) ) when XROC is 'R' or 'r'
*          and MAX( 1, Lr( 1, IX+Lx-1 ) ) otherwise, and, Kx is at least
*          Lc( 1, JX+Lx-1 )  when INCX = M_X  and Lc( K, JX ) otherwise.
*          Lx is N when ROWCOL = 'R' or 'r' and M  otherwise. Before en-
*          try, this array  contains the local entries of the  matrix X.
*
*  IX      (global input) INTEGER
*          On entry, IX  specifies X's global row index, which points to
*          the beginning of the submatrix sub( X ).
*
*  JX      (global input) INTEGER
*          On entry, JX  specifies X's global column index, which points
*          to the beginning of the submatrix sub( X ).
*
*  DESCX   (global and local input) INTEGER array
*          On entry, DESCX  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix X.
*
*  XROC    (global input) pointer to CHAR
*          On entry,  XROC  specifies  the  orientation of the subvector
*          sub( X ). When XROC is 'R' or 'r',  sub( X ) is a row vector,
*          and a column vector otherwise.
*
*  XAPTR   (local output) pointer to pointer to CHAR
*          On exit, * XAPTR is  an array containing the same data as the
*          subvector  sub( X )  which is replicated over the rows or co-
*          lumns of the  underlying  matrix  as  specified by ROWCOL and
*          DESCA.
*
*  DXA     (global and local output) INTEGER array
*          On exit, DXA is a descriptor array of dimension DLEN_ descri-
*          bing the data layout of the data pointed to by * XAPTR.
*
*  XAFREE  (local output) INTEGER
*          On exit,  XAFREE  specifies  if it has been possible to reuse
*          the subvector sub( X ), i.e., if some dynamic  memory was al-
*          located for the data pointed to by *XAPTR or not. When XAFREE
*          is zero, no dynamic memory was allocated. Otherwise, some dy-
*          namic memory was allocated by this function that one MUST re-
*          lease as soon as possible.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           * top;
   int            AColSpan, ARowSpan, Acol, Aimb, Ainb, AisD, Amb, Amp, Anb,
                  Anq, Arow, Xcol, Xii, Ximb, Ximb1, Xinb, Xinb1, XisD, XisR,
                  XisRow, Xjj, Xld=1, Xmb, Xmp, Xnb, Xnq, Xrow, ctxt, mycol,
                  myrow, npcol, nprow;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Initialize the output parameters to a default value
*/
   *XAFREE = 0;
   *XAPTR  = NULL;
/*
*  Quick return if possible
*/
   if( ( M <= 0 ) || ( N <= 0 ) || ( K <= 0 ) )
   {
      if( Mupcase( ROWCOL[0] ) == CROW )
      {
         PB_Cdescset( DXA, K, N, 1, DESCA[INB_], 1, DESCA[NB_], DESCA[RSRC_],
                      DESCA[CSRC_], DESCA[CTXT_], 1 );
      }
      else
      {
         PB_Cdescset( DXA, M, K, DESCA[IMB_], 1, DESCA[MB_], 1, DESCA[RSRC_],
                      DESCA[CSRC_], DESCA[CTXT_], DESCA[LLD_] );
      }
      return;
   }
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCX[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Retrieve sub( X )'s local information: Xii, Xjj, Xrow, Ycol
*/
   Minfog2l( IX, JX, DESCX, nprow, npcol, myrow, mycol, Xii, Xjj, Xrow, Xcol );
/*
*  Is sub( X ) distributed or not, replicated or not ?
*/
   if( ( XisRow = ( Mupcase( XROC[0] ) == CROW ) ) != 0 )
   {
      XisD = ( ( Xcol >=  0 ) && ( npcol >  1 ) );
      XisR = ( ( Xrow == -1 ) || ( nprow == 1 ) );
   }
   else
   {
      XisD = ( ( Xrow >=  0 ) && ( nprow >  1 ) );
      XisR = ( ( Xcol == -1 ) || ( npcol == 1 ) );
   }

   Arow = DESCA[ RSRC_ ]; Acol = DESCA[ CSRC_ ];

   if( Mupcase( ROWCOL[0] ) == CROW )
   {
/*
*  Want a row vector
*/
      Ainb = DESCA[ INB_  ]; Anb  = DESCA[ NB_   ];
      Mnumroc( Anq, N, 0, Ainb, Anb, mycol, Acol, npcol );
/*
*  Does A spans multiples process rows ? It does if Arow < 0.
*/
      ARowSpan = ( Arow < 0 ) ||
                 Mspan( M, 0, DESCA[IMB_], DESCA[MB_], Arow, nprow );

      if( XisRow && ( Mupcase( CONJUG[0] ) == CNOCONJG ) )
      {
/*
*  It is possible to reuse sub( X ) iff sub( X ) is already a row vector and
*  the data does not need to be conjugated.
*/
         AisD = ( ( Acol >= 0 ) && ( npcol > 1 ) );

         Xinb = DESCX[ INB_ ]; Xnb  = DESCX[ NB_  ];
         Mfirstnb( Xinb1, N, JX, Xinb, Xnb );
/*
*  sub( X ) is aligned with A (reuse condition) iff both operands are not
*  distributed, or both of them are distributed and start in the same process
*  column and either N is smaller than the first blocksize of sub( X ) and A,
*  or their column blocking factors match.
*/
         if( ( !AisD && !XisD ) ||
             ( ( AisD && XisD )  &&
               ( ( Acol == Xcol ) &&
                 ( ( ( Ainb >= N     ) && ( Xinb1 >= N ) ) ||
                   ( ( Ainb == Xinb1 ) && ( Anb == Xnb ) ) ) ) ) )
         {
/*
*  sub( X ) is aligned with A
*/
            Ximb = DESCX[ IMB_  ]; Xmb  = DESCX[ MB_   ];
            Mfirstnb( Ximb1, K, IX, Ximb, Xmb );

            if( XisR || ( !ARowSpan && ( Arow == Xrow ) ) )
            {
/*
*  If sub( X ) is replicated, or, A spans only one process row and either
*  sub( X ) is replicated or resides in the same process row than A, then
*  sub( X ) is already at the correct place.
*/
               if( Anq > 0 )
               {
                  Xld = DESCX[ LLD_ ];
                  if( ARowSpan || ( myrow == Arow ) )
                     *XAPTR = Mptr( X, Xii, Xjj, Xld, TYPE->size );
               }
               else
               {
                  Xld = 1;
               }
               MDescSet( DXA, K, N, K, Xinb1, 1, Xnb, ( ARowSpan ? -1 : Arow ),
                         Xcol, ctxt, Xld );
            }
            else if( ARowSpan )
            {
/*
*  Otherwise, we know that sub( X ) cannot be replicated, let suppose in
*  addition that A spans all process rows. sub( X ) need simply to be broadcast
*  over A.
*/
               if( myrow == Xrow )
               {
                  Xld = DESCX[ LLD_ ];
                  if( Anq > 0 )
                  {
                     *XAPTR = Mptr( X, Xii, Xjj, Xld, TYPE->size );
                     top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
                     TYPE->Cgebs2d( ctxt, COLUMN, top, K, Anq, *XAPTR, Xld );
                  }
               }
               else
               {
                  Xld = MAX( 1, K );
                  if( Anq > 0 )
                  {
                     *XAPTR  = PB_Cmalloc( K * Anq * TYPE->size );
                     *XAFREE = 1;
                     top = PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
                     TYPE->Cgebr2d( ctxt, COLUMN, top, K, Anq, *XAPTR, Xld,
                                    Xrow, mycol );
                  }
               }
               PB_Cdescset( DXA, K, N, K, Xinb1, 1, Xnb, -1, Xcol, ctxt, Xld );
            }
            else
            {
/*
*  Finally, sub( X ) is not replicated and A spans only one process row. There
*  is no need to broadcast, a send/recv is sufficient.
*/
               if( myrow == Xrow )
               {
                  Xld = DESCX[ LLD_ ];
                  if( Anq > 0 )
                  {
                     *XAPTR = Mptr( X, Xii, Xjj, Xld, TYPE->size );
                     TYPE->Cgesd2d( ctxt, K, Anq, *XAPTR, Xld, Arow, mycol );
                  }
               }
               else if( myrow == Arow )
               {
                  Xld = MAX( 1, K );
                  if( Anq > 0 )
                  {
                     *XAPTR = PB_Cmalloc( K * Anq * TYPE->size );
                     *XAFREE = 1;
                     TYPE->Cgerv2d( ctxt, K, Anq, *XAPTR, Xld, Xrow, mycol );
                  }
               }
               PB_Cdescset( DXA, K, N, K, Xinb1, 1, Xnb, Arow, Xcol, ctxt,
                            Xld );
            }
            return;
         }
      }
/*
*  sub( X ) cannot be reused, too bad ... redistribute
*/
      PB_Cdescset( DXA, K, N, K, Ainb, 1, Anb, ( ARowSpan ? -1 : Arow ), Acol,
                   ctxt, K );
      Xmp = ( ARowSpan ? K : ( ( myrow == Arow ) ? K : 0 ) );
      if( Xmp > 0 &&  Anq > 0 )
      {
         *XAPTR  = PB_Cmalloc( Anq * Xmp * TYPE->size );
         *XAFREE = 1;
      }
      if( XisRow )
      {
         PB_Cpaxpby( TYPE, CONJUG, K, N, TYPE->one, X, IX, JX, DESCX, XROC,
                     TYPE->zero, *XAPTR, 0, 0, DXA, ROW );
      }
      else
      {
         PB_Cpaxpby( TYPE, CONJUG, N, K, TYPE->one, X, IX, JX, DESCX, XROC,
                     TYPE->zero, *XAPTR, 0, 0, DXA, ROW );
      }
   }
   else
   {
/*
*  Want a column vector
*/
      Aimb = DESCA[IMB_]; Amb  = DESCA[MB_];
      Mnumroc( Amp, M, 0, Aimb, Amb, myrow, Arow, nprow );
/*
*  Does A spans multiples process columns ? It does if Acol < 0.
*/
      AColSpan = ( Acol < 0 ) ||
                 Mspan( N, 0, DESCA[INB_], DESCA[NB_], Acol, npcol );

      if( !( XisRow ) && ( Mupcase( CONJUG[0] ) == CNOCONJG ) )
      {
/*
*  It is possible to reuse sub( X ) iff sub( X ) is already a column vector and
*  the data does not need to be conjugated
*/
         AisD = ( ( Arow >=  0 ) && ( nprow > 1 ) );

         Ximb = DESCX[ IMB_ ]; Xmb  = DESCX[ MB_  ];
         Mfirstnb( Ximb1, M, IX, Ximb, Xmb );
/*
*  sub( X ) is aligned with A (reuse condition) iff both operands are not
*  distributed, or both of them are distributed and start in the same process
*  row and either M is smaller than the first blocksize of sub( X ) and A, or
*  their row blocking factors match.
*/
         if( ( !AisD && !XisD ) ||
             ( ( AisD && XisD )  &&
               ( ( Arow == Xrow ) &&
                 ( ( ( Aimb >= M     ) && ( Ximb1 >= M ) ) ||
                   ( ( Aimb == Ximb1 ) && ( Amb == Xmb ) ) ) ) ) )
         {
/*
*  sub( X ) is aligned with A
*/
            Xinb = DESCX[ INB_  ]; Xnb  = DESCX[ NB_   ];
            Mfirstnb( Xinb1, K, JX, Xinb, Xnb );

            if( XisR || ( !AColSpan && ( Acol == Xcol ) ) )
            {
/*
*  If sub( X ) is replicated, or, A spans only one process column and either
*  sub( X ) is replicated or resides in the same process columns than A, then
*  sub( X ) is already at the correct place.
*/
               if( Amp > 0 )
               {
                  Xld = DESCX[ LLD_ ];
                  if( AColSpan || ( mycol == Acol ) )
                     *XAPTR = Mptr( X, Xii, Xjj, Xld, TYPE->size );
               }
               else
               {
                  Xld = 1;
               }
               MDescSet( DXA, M, K, Ximb1, K, Xmb, 1, Xrow,
                         ( AColSpan ? -1 : Acol ), ctxt, Xld );
            }
            else if( AColSpan )
            {
/*
*  Otherwise, we know that sub( X ) is not be replicated, let suppose in
*  addition that A spans all process columns. sub( X ) need simply to be
*  broadcast over A.
*/
               if( mycol == Xcol )
               {
                  Xld = DESCX[ LLD_ ];
                  if( Amp > 0 )
                  {
                     *XAPTR = Mptr( X, Xii, Xjj, Xld, TYPE->size );
                     top = PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );
                     TYPE->Cgebs2d( ctxt, ROW, top, Amp, K, *XAPTR, Xld );
                  }
               }
               else
               {
                  Xld = MAX( 1, Amp );
                  if( Amp > 0 )
                  {
                     *XAPTR = PB_Cmalloc( Amp * K * TYPE->size );
                     *XAFREE = 1;
                     top = PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );
                     TYPE->Cgebr2d( ctxt, ROW, top, Amp, K, *XAPTR, Xld, myrow,
                                    Xcol );
                  }
               }
               PB_Cdescset( DXA, M, K, Ximb1, K, Xmb, 1, Xrow, -1, ctxt, Xld );
            }
            else
            {
/*
*  Finally, sub( X ) is not replicated and A spans only one process column.
*  There is no need to broadcast, a send/recv is sufficient.
*/
               if( mycol == Xcol )
               {
                  Xld = DESCX[ LLD_ ];
                  if( Amp > 0 )
                  {
                     *XAPTR = Mptr( X, Xii, Xjj, Xld, TYPE->size );
                     TYPE->Cgesd2d( ctxt, Amp, K, *XAPTR, Xld, myrow, Acol );
                  }
               }
               else if( mycol == Acol )
               {
                  Xld = MAX( 1, Amp );
                  if( Amp > 0 )
                  {
                     *XAPTR = PB_Cmalloc( Amp * K * TYPE->size );
                     *XAFREE = 1;
                     TYPE->Cgerv2d( ctxt, Amp, K, *XAPTR, Xld, myrow, Xcol );
                  }
               }
               PB_Cdescset( DXA, M, K, Ximb1, K, Xmb, 1, Xrow, Acol, ctxt,
                            Xld );
            }
            return;
         }
      }
/*
*  sub( X ) cannot be reused, too bad ... redistribute
*/
      PB_Cdescset( DXA, M, K, Aimb, K, Amb, 1, Arow, ( AColSpan ? -1 : Acol ),
                   ctxt, MAX( 1, Amp ) );
      Xnq = ( AColSpan ? K : ( ( mycol == Acol ) ? K : 0 ) );
      if( Xnq > 0 &&  Amp > 0 )
      {
         *XAPTR  = PB_Cmalloc( Amp * Xnq * TYPE->size );
         *XAFREE = 1;
      }
      if( XisRow )
      {
         PB_Cpaxpby( TYPE, CONJUG, K, M, TYPE->one, X, IX, JX, DESCX, XROC,
                     TYPE->zero, *XAPTR, 0, 0, DXA, COLUMN );
      }
      else
      {
         PB_Cpaxpby( TYPE, CONJUG, M, K, TYPE->one, X, IX, JX, DESCX, XROC,
                     TYPE->zero, *XAPTR, 0, 0, DXA, COLUMN );
      }
   }
/*
*  End of PB_CInV
*/
}
