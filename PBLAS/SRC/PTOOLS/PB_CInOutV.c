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
void PB_CInOutV( PBTYP_T * TYPE, char * ROWCOL, Int M, Int N, Int * DESCA,
                 Int K, char * BETA,
                 char * Y, Int IY, Int JY, Int * DESCY, char * YROC,
                 char * * TBETA, char * * YAPTR, Int * DYA,
                 Int * YAFREE, Int * YASUM, Int * YAPBY )
#else
void PB_CInOutV( TYPE, ROWCOL, M, N, DESCA, K,
                 BETA, Y, IY, JY, DESCY, YROC,
                 TBETA, YAPTR, DYA, YAFREE, YASUM, YAPBY )
/*
*  .. Scalar Arguments ..
*/
   char           * BETA, * ROWCOL, * * TBETA, * YROC;
   Int            * YAPBY, * YAFREE, IY, JY, K, M, N, * YASUM;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA, * DESCY, * DYA;
   char           * Y, * * YAPTR;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CInOutV  returns  a pointer to an array that contains a one-dimen-
*  sional  input/output subvector  which  is replicated over the rows or
*  columns of a submatrix described by DESCA. A  subvector is  specified
*  on input to this routine that is reused whenever possible. On return,
*  the subvector is specified  by  a  pointer to some data, a descriptor
*  array describing its layout, a logical value indicating if this local
*  piece of data has been dynamically allocated by this function, a  lo-
*  gical  value  specifying if sum reduction should occur, and finally a
*  logical  value specifying if it is necessary to copy back the alloca-
*  ted data to the original data. This routine is specifically  designed
*  for  traditional Level 2 like PBLAS  operations using an input/output
*  vector such as PxGEMV, PxSYMV ...
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
*  DESCA   (global and local input/output) INTEGER array
*          On entry, DESCA is an integer array of dimension DLEN_.  This
*          is the array descriptor for the matrix A. EXCEPTIONALLY, THIS
*          INTERNAL  ROUTINE  MAY  MODIFY DESCA IN ORDER TO MINIMIZE THE
*          AMOUNT OF DATA TO BE MOVED FOR THE VECTOR Y. SEE  PxGEMV  FOR
*          AN EXAMPLE.
*
*  K       (global input) INTEGER
*          On entry,  K  specifies the length of the non-distributed di-
*          mension of the subvector sub( Y ). K must be at least zero.
*
*  BETA    (global input) pointer to CHAR
*          On entry, BETA  is a scalar the input subvector sub( Y ) must
*          be scaled by.
*
*  Y       (local input) pointer to CHAR
*          On entry, Y is an array of dimension (LLD_Y, Ky), where LLD_Y
*          is   at  least  MAX( 1, Lr( K, IY ) ) when YROC is 'R' or 'r'
*          and  MAX( 1, Lr( 1, IY+Ly-1 ) )  otherwise,  and,  Ky  is  at
*          least   Lc( 1, JY+Ly-1 )   when   YROC  is  'R'  or  'r'  and
*          Lc( K, JY ) otherwise.  Ly is N when ROWCOL is 'R' or 'r' and
*          M  otherwise.  Before  entry,  this array  contains the local
*          entries of the  matrix Y.
*
*  IY      (global input) INTEGER
*          On entry, IY  specifies Y's global row index, which points to
*          the beginning of the submatrix sub( Y ).
*
*  JY      (global input) INTEGER
*          On entry, JY  specifies Y's global column index, which points
*          to the beginning of the submatrix sub( Y ).
*
*  DESCY   (global and local input) INTEGER array
*          On entry, DESCY  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix Y.
*
*  YROC    (global input) pointer to CHAR
*          On entry,  YROC  specifies  the  orientation of the subvector
*          sub( Y ). When YROC is 'R' or 'r',  sub( Y ) is a row vector,
*          and a column vector otherwise.
*
*  TBETA   (local output) pointer to pointer to CHAR
*          On exit, * TBETA  is a scalar to be used locally to scale the
*          data pointed to by * YAPTR,  in  order  to obtain the correct
*          result in the original data sub( Y ).
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
*  YAPBY   (global output) INTEGER
*          On exit, YAPBY  specifies  if  the data pointed to by * YAPTR
*          must be move back onto sub( Y ) to obtain the correct result.
*          When YAPBY is zero, no supplementary data movement is  neces-
*          sary, otherwise a data redistribution should occur.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   Int            Acol, Aimb, Ainb, AisD, AisR, Amb, Amp, Anb, Anq, Arow, Ycol,
                  Yii, Yimb, Yimb1, Yinb, Yinb1, YisD, YisR, YisRow, Yjj, Yld,
                  Ymb, Ymp, Ynb, Ynq, Yrow, ctxt, izero=0, nprow, myrow, npcol,
                  mycol;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Initialize the output parameters to a default value
*/
   *YAFREE = 0;
   *YASUM  = 0;
   *YAPBY  = 0;
   *YAPTR  = NULL;
   *TBETA  = BETA;
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
   Cblacs_gridinfo( ( ctxt = DESCY[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Retrieve sub( Y )'s local information: Yii, Yjj, Yrow, Ycol
*/
   Minfog2l( IY, JY, DESCY, nprow, npcol, myrow, mycol, Yii, Yjj, Yrow, Ycol );
/*
*  Is sub( Y ) distributed or not, replicated or not ?
*/
   if( ( YisRow = ( Mupcase( YROC[0] ) == CROW ) ) != 0 )
   {
      YisD = ( ( Ycol >=  0 ) && ( npcol >  1 ) );
      YisR = ( ( Yrow == -1 ) || ( nprow == 1 ) );
   }
   else
   {
      YisD = ( ( Yrow >=  0 ) && ( nprow >  1 ) );
      YisR = ( ( Ycol == -1 ) || ( npcol == 1 ) );
   }

   Aimb = DESCA[IMB_ ]; Ainb = DESCA[INB_ ];
   Amb  = DESCA[MB_  ]; Anb  = DESCA[NB_  ];
   Arow = DESCA[RSRC_]; Acol = DESCA[CSRC_];

   if( Mupcase( ROWCOL[0] ) == CROW )
   {
/*
*  Want a row vector
*/
      AisR = ( ( Arow < 0 ) || ( nprow == 1 ) );

      if( YisRow )
      {
/*
*  It is possible to reuse sub( Y ) iff sub( Y ) is already a row vector.
*/
         AisD = ( ( Acol >= 0 ) && ( npcol > 1 ) );

         Yinb = DESCY[INB_]; Ynb = DESCY[NB_];
         Mfirstnb( Yinb1, N, JY, Yinb, Ynb );
/*
*  sub( Y ) is aligned with A (reuse condition) iff both operands are not
*  distributed, or both of them are distributed and start in the same process
*  column and either N is smaller than the first blocksize of sub( Y ) and A,
*  or their column blocking factors match.
*/
         if( ( !AisD && !YisD ) ||
             ( ( AisD && YisD )  &&
               ( ( Acol == Ycol ) &&
                 ( ( ( Ainb >= N     ) && ( Yinb1 >= N ) ) ||
                   ( ( Ainb == Yinb1 ) && ( Anb == Ynb ) ) ) ) ) )
         {
            Mnumroc( Ynq, N, 0, Yinb1, Ynb, mycol, Ycol, npcol );
            Ymp = ( YisR ? K : ( ( myrow == Yrow ) ? K : 0 ) );

            if( YisR )
            {
/*
*  If sub( Y ) is replicated, there is no need to move sub( Y ) after the
*  operation (*YAPBY = 0), and it can be reused.
*/
               *YAPBY = 0;
               Yld    = DESCY[ LLD_ ];
               if( Ynq > 0 )
                  *YAPTR = Mptr( Y, Yii, Yjj, Yld, TYPE->size );

               if( AisR )
               {
/*
*  If A is replicated as well, use BETA in every process row, and do not combine
*  the local results.
*/
                  *TBETA = BETA;
                  *YASUM = 0;
               }
               else
               {
/*
*  Otherwise, use BETA in process row Arow and zero elsewhere. Reduce the local
*  result if there is more than one row in the process grid.
*/
                  *TBETA = ( ( myrow == Arow ) ? BETA : TYPE->zero );
                  *YASUM = ( nprow > 1 );
/*
*  If some process rows do not own any entries of A, better set sub( Y ) to zero
*  in those processes.
*/
                  Mnumroc( Amp, M, 0, Aimb, Amb, myrow, Arow, nprow );
                  if( Amp <= 0 )
                     TYPE->Ftzscal( C2F_CHAR( ALL ), &K, &Ynq, &izero, *TBETA,
                                    *YAPTR, &Yld );
               }
            }
            else
            {
/*
*  sub( Y ) is not replicated, the descriptor of A may need to be modified ...
*/
               if( AisR )
               {
/*
*  If A is replicated, use only the copy in the process row where sub( Y )
*  resides -> modify DESCA !!!
*/
                  *TBETA = BETA;
                  *YASUM = 0;
                  *YAPBY = 0;
                  Yld    = DESCY[ LLD_ ];
                  DESCA[ IMB_  ] = M;
                  DESCA[ RSRC_ ] = Yrow;
                  if( ( Ynq > 0 ) && ( Ymp > 0 ) )
                     *YAPTR = Mptr( Y, Yii, Yjj, Yld, TYPE->size );
               }
               else
               {
                  if( Mspan( M, 0, Aimb, Amb, Arow, nprow ) )
                  {
/*
*  Otherwise, A is not replicated, let assume in addition that it spans more
*  than one process row
*/
                     *YASUM = ( nprow > 1 );
                     *YAPBY = 0;

                     if( myrow == Yrow )
                     {
/*
*  Reuse sub( Y ). If there is no entries of A in the process row where sub( Y )
*  resides, better scale it by BETA immediately.
*/
                        *TBETA = BETA;
                        Yld    = DESCY[ LLD_ ];
                        if( Ynq > 0 )
                        {
                           *YAPTR = Mptr( Y, Yii, Yjj, Yld, TYPE->size );
                           Mnumroc( Amp, M, 0, Aimb, Amb, myrow, Arow, nprow );
                           if( Amp <= 0 )
                              TYPE->Ftzscal( C2F_CHAR( ALL ), &K, &Ynq, &izero,
                                             *TBETA, *YAPTR, &Yld );
                        }
                     }
                     else
                     {
/*
*  Allocate space in the other process rows and initialize to zero.
*/
                        *TBETA = TYPE->zero;
                        Yld    = MAX( 1, K );
                        if( Ynq > 0 )
                        {
                           *YAPTR  = PB_Cmalloc( K * Ynq * TYPE->size );
                           *YAFREE = 1;
                           TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ),
                                         &K, &Ynq, &izero, *TBETA, *TBETA,
                                         *YAPTR, &Yld );
                        }
                     }
                  }
                  else
                  {
/*
*  A spans only one process row
*/
                     if( Yrow == Arow )
                     {
/*
*  A and sub( Y ) resides in the same process row
*/
                        *TBETA = BETA;
                        *YASUM = 0;
                        *YAPBY = 0;
                        Yld    = DESCY[ LLD_ ];
                        if( ( myrow == Yrow ) && ( Ynq > 0 ) )
                           *YAPTR = Mptr( Y, Yii, Yjj, Yld, TYPE->size );
                     }
                     else
                     {
/*
*  If sub( Y ) resides in another process row, then allocate zero-data in
*  process row where A resides, and set *YAPBY to 1, so that this data will be
*  added (moved) after the local operation has been performed.
*/
                        *TBETA = TYPE->zero;
                        *YASUM = 0;
                        *YAPBY = 1;
                        Yrow   = Arow;
                        Yld    = MAX( 1, K );
                        if( myrow == Arow )
                        {
                           if( Ynq > 0 )
                           {
                              *YAPTR  = PB_Cmalloc( K * Ynq * TYPE->size );
                              *YAFREE = 1;
                              TYPE->Ftzpad( C2F_CHAR( ALL ),
                                            C2F_CHAR( NOCONJG ), &K, &Ynq,
                                            &izero, *TBETA, *TBETA, *YAPTR,
                                            &Yld );
                           }
                        }
                     }
                  }
               }
            }
/*
*  Describe the resulting operand. Note that when reduction should occur, Yrow
*  contains the destination row. Assuming every process row needs the result,
*  Yrow is then -1.
*/
            MDescSet( DYA, K, N, K, Yinb1, 1, Ynb, Yrow, Ycol, ctxt, Yld );
            return;
         }
      }
/*
*  sub( Y ) cannot be reused, set TBETA to zero for the local operation, and
*  force YAPBY to 1 for the later update of sub( Y ).
*/
      *TBETA = TYPE->zero;
      *YAPBY = 1;
      Mnumroc( Anq, N, 0, Ainb, Anb, mycol, Acol, npcol );
      Yld    = MAX( 1, K );

      if( YisR )
      {
/*
*  If sub( Y ) is replicated, allocate space in every process row owning some
*  columns of A and initialize it to zero. There may be some wasted space
*  (suppose A was residing in just one row), however, it is hoped that moving
*  back this data to sub( Y ) will then be cheaper ...
*/
         *YASUM = ( AisR ? 0 : ( nprow > 1 ) );
         Arow   = -1;
         if( Anq > 0 )
         {
            *YAPTR  = PB_Cmalloc( K * Anq * TYPE->size );
            *YAFREE = 1;
            TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &K, &Anq,
                          &izero, *TBETA, *TBETA, *YAPTR, &Yld );
         }
      }
      else
      {
/*
*  sub( Y ) resides within only one process row
*/
         if( AisR )
         {
/*
*  If A is replicated, then modify sub( A ) so that only one process row will
*  compute the result before moving it back to sub( Y ).
*/
            *YASUM         = 0;
            DESCA[ IMB_  ] = M;

            if( YisRow )
            {
/*
*  Choose different process row than Yrow for better performance (more links)
*  of the later move-back phase.
*/
               DESCA[RSRC_] = Arow = MModSub1( Yrow, nprow );
            }
            else
            {
               DESCA[RSRC_] = Arow = 0;
            }
            if( ( myrow == Arow ) && ( Anq > 0 ) )
            {
               *YAPTR  = PB_Cmalloc( K * Anq * TYPE->size );
               *YAFREE = 1;
               TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &K, &Anq,
                             &izero, *TBETA, *TBETA, *YAPTR, &Yld );
            }
         }
         else
         {
            if( Mspan( M, 0, Aimb, Amb, Arow, nprow ) )
            {
/*
*  If A is not replicated, and spans more than just one process row, then
*  allocate space in every process row and zero it.
*/
               *YASUM = ( nprow > 1 );
               if( Anq > 0 )
               {
                  *YAPTR  = PB_Cmalloc( K * Anq * TYPE->size );
                  *YAFREE = 1;
                  TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &K, &Anq,
                                &izero, *TBETA, *TBETA, *YAPTR, &K );
               }
            }
            else
            {
/*
*  If A is not replicated, and spans only one process row, then allocate space
*  within that process row and zero it.
*/
               *YASUM = 0;
               if( ( myrow == Arow ) && ( Anq > 0 ) )
               {
                  *YAPTR  = PB_Cmalloc( K * Anq * TYPE->size );
                  *YAFREE = 1;
                  TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &K, &Anq,
                                &izero, *TBETA, *TBETA, *YAPTR, &K );
               }
            }
         }
      }
/*
*  Describe the resulting operand. Note that when reduction should occur, Arow
*  contains the destination row. Assuming every process row needs the result,
*  Arow is then -1.
*/
      MDescSet( DYA, K, N, K, Ainb, 1, Anb, Arow, Acol, ctxt, Yld );
   }
   else
   {
/*
*  Want a column vector
*/
      AisR = ( ( Acol <  0 ) || ( npcol == 1 ) );

      if( !YisRow )
      {
/*
*  It is possible to reuse sub( Y ) iff sub( Y ) is already a column vector.
*/
         AisD = ( ( Arow >= 0 ) && ( nprow >  1 ) );

         Yimb = DESCY[IMB_]; Ymb = DESCY[MB_];
         Mfirstnb( Yimb1, M, IY, Yimb, Ymb );
/*
*  sub( Y ) is aligned with A (reuse condition) iff both operands are not
*  distributed, or both of them are distributed and start in the same process
*  row and either M is smaller than the first blocksize of sub( Y ) and A, or
*  their row blocking factors match.
*/
         if( ( !AisD && !YisD ) ||
             ( ( AisD && YisD )  &&
               ( ( Arow == Yrow ) &&
                 ( ( ( Aimb >= M     ) && ( Yimb1 >= M ) ) ||
                   ( ( Aimb == Yimb1 ) && ( Amb == Ymb ) ) ) ) ) )
         {
            Mnumroc( Ymp, M, 0, Yimb1, Ymb, myrow, Yrow, nprow );
            Ynq = ( YisR ? K : ( ( mycol == Ycol ) ? K : 0 ) );

            if( YisR )
            {
/*
*  If sub( Y ) is replicated, there is no need to move sub( Y ) after the
*  operation (*YAPBY = 0), and it can be reused.
*/
               *YAPBY = 0;
               Yld    = DESCY[ LLD_ ];
               if( Ymp > 0 )
                  *YAPTR = Mptr( Y, Yii, Yjj, Yld, TYPE->size );

               if( AisR )
               {
/*
*  If A is replicated as well, use BETA in every process column, and do not
*  combine the local results.
*/
                  *TBETA = BETA;
                  *YASUM = 0;
               }
               else
               {
/*
*  Otherwise, use BETA in process column Acol and zero elsewhere. Reduce the
*  local result if there is more than one column in the process grid.
*/
                  *TBETA = ( ( mycol == Acol ) ? BETA : TYPE->zero );
                  *YASUM = ( npcol > 1 );
/*
*  If some process columns do not own any entries of A, better set sub( Y ) to
*  zero in those processes.
*/
                  Mnumroc( Anq, N, 0, Ainb, Anb, mycol, Acol, npcol );
                  if( Anq <= 0 )
                     TYPE->Ftzscal( C2F_CHAR( ALL ), &Ymp, &K, &izero, *TBETA,
                                    *YAPTR, &Yld );
               }
            }
            else
            {
/*
*  sub( Y ) is not replicated, the descriptor of A may need to be modified ...
*/
               if( AisR )
               {
/*
*  If A is replicated, use only the copy in the process column where sub( Y )
*  resides -> modify DESCA !!!
*/
                  *TBETA = BETA;
                  *YASUM = 0;
                  *YAPBY = 0;
                  Yld    = DESCY[ LLD_ ];
                  DESCA[ INB_  ] = N;
                  DESCA[ CSRC_ ] = Ycol;
                  if( ( Ymp > 0 ) && ( Ynq > 0 ) )
                     *YAPTR = Mptr( Y, Yii, Yjj, Yld, TYPE->size );
               }
               else
               {
                  if( Mspan( N, 0, Ainb, Anb, Acol, npcol ) )
                  {
/*
*  Otherwise, A is not replicated, let assume in addition that it spans more
*  than one process column
*/
                     *YASUM = ( npcol > 1 );
                     *YAPBY = 0;

                     if( mycol == Ycol )
                     {
/*
*  Reuse sub( Y ). If there is no entries of A in the process column where
*  sub( Y ) resides, better scale it by BETA immediately.
*/
                        *TBETA = BETA;
                        Yld    = DESCY[ LLD_ ];
                        if( Ymp > 0 )
                        {
                           *YAPTR = Mptr( Y, Yii, Yjj, Yld, TYPE->size );
                           Mnumroc( Anq, N, 0, Ainb, Anb, mycol, Acol, npcol );
                           if( Anq <= 0 )
                              TYPE->Ftzscal( C2F_CHAR( ALL ), &Ymp, &K, &izero,
                                             *TBETA, *YAPTR, &Yld );
                        }
                     }
                     else
                     {
/*
*  Allocate space in the other process columns and initialize to zero.
*/
                        *TBETA = TYPE->zero;
                        Yld    = MAX( 1, Ymp );
                        if( Ymp > 0 )
                        {
                           *YAPTR  = PB_Cmalloc( Ymp * K * TYPE->size );
                           *YAFREE = 1;
                           TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ),
                                         &Ymp, &K, &izero, *TBETA, *TBETA,
                                         *YAPTR, &Yld );
                        }
                     }
                  }
                  else
                  {
/*
*  A spans only one process column
*/
                     if( Ycol == Acol )
                     {
/*
*  A and sub( Y ) resides in the same process column
*/
                        *TBETA = BETA;
                        *YASUM = 0;
                        *YAPBY = 0;
                        Yld    = DESCY[ LLD_ ];
                        if( ( mycol == Ycol ) && ( Ymp > 0 ) )
                           *YAPTR = Mptr( Y, Yii, Yjj, Yld, TYPE->size );
                     }
                     else
                     {
/*
*  If sub( Y ) resides in another process column, then allocate zero-data in
*  process column where A resides, and set *YAPBY to 1, so that this data will
*  be added (moved) after the local operation has been performed.
*/
                        *TBETA = TYPE->zero;
                        *YASUM = 0;
                        *YAPBY = 1;
                        Ycol   = Acol;
                        Yld    = MAX( 1, Ymp ) ;
                        if( mycol == Acol )
                        {
                           if( Ymp > 0 )
                           {
                              *YAPTR  = PB_Cmalloc( Ymp * K * TYPE->size );
                              *YAFREE = 1;
                              TYPE->Ftzpad( C2F_CHAR( ALL ),
                                            C2F_CHAR( NOCONJG ), &Ymp, &K,
                                            &izero, *TBETA, *TBETA, *YAPTR,
                                            &Yld );
                           }
                        }
                     }
                  }
               }
            }
/*
*  Describe the resulting operand. Note that when reduction should occur, Ycol
*  contains the destination column. Assuming every process column needs the
*  result, Ycol is then -1.
*/
            MDescSet( DYA, M, K, Yimb1, K, Ymb, 1, Yrow, Ycol, ctxt, Yld );
            return;
         }
      }
/*
*  sub( Y ) cannot be reused, set TBETA to zero for the local operation, and
*  force YAPBY to 1 for the later update of sub( Y ).
*/
      *TBETA = TYPE->zero;
      *YAPBY = 1;
      Mnumroc( Amp, M, 0, Aimb, Amb, myrow, Arow, nprow );
      Yld    = MAX( 1, Amp );

      if( YisR )
      {
/*
*  If sub( Y ) is replicated, allocate space in every process column owning some
*  rows of A and initialize it to zero. There may be some wasted space (suppose
*  A was residing in just one column), however, it is hoped that moving back
*  this data to sub( Y ) will then be cheaper ...
*/
         *YASUM = ( AisR ? 0 : ( npcol > 1 ) );
         Acol   = -1;
         if( Amp > 0 )
         {
            *YAPTR  = PB_Cmalloc( Amp * K * TYPE->size );
            *YAFREE = 1;
            TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Amp, &K,
                          &izero, *TBETA, *TBETA, *YAPTR, &Yld );
         }
      }
      else
      {
/*
*  sub( Y ) resides within only one process column
*/
         if( AisR )
         {
/*
*  If A is replicated, then modify sub( A ) so that only one process column will
*  compute the result before moving it back to sub( Y ).
*/
            *YASUM         = 0;
            DESCA[ INB_  ] = N;

            if( YisRow )
            {
               DESCA[ CSRC_ ] = Acol = 0;
            }
            else
            {
/*
*  Choose different process column than Ycol for better performance (more links)
*  of the later move-back phase.
*/
               DESCA[ CSRC_ ] = Acol = MModSub1( Ycol, npcol );
            }
            if( ( mycol == Acol ) && ( Amp > 0 ) )
            {
               *YAPTR  = PB_Cmalloc( Amp * K * TYPE->size );
               *YAFREE = 1;
               TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Amp, &K,
                             &izero, *TBETA, *TBETA, *YAPTR, &Yld );
            }
         }
         else
         {
            if( Mspan( N, 0, Ainb, Anb, Acol, npcol ) )
            {
/*
*  If A is not replicated, and spans more than just one process column, then
*  allocate space in every process column and zero it.
*/
               *YASUM = ( npcol > 1 );
               if( Amp > 0 )
               {
                  *YAPTR  = PB_Cmalloc( Amp * K * TYPE->size );
                  *YAFREE = 1;
                  TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Amp, &K,
                                &izero, *TBETA, *TBETA, *YAPTR, &Yld );
               }
            }
            else
            {
/*
*  If A is not replicated, and spans only one process column, then allocate
*  space within that process column and zero it.
*/
               *YASUM = 0;
               if( ( mycol == Acol ) && ( Amp > 0 ) )
               {
                  *YAPTR  = PB_Cmalloc( Amp * K * TYPE->size );
                  *YAFREE = 1;
                  TYPE->Ftzpad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Amp, &K,
                                &izero, *TBETA, *TBETA, *YAPTR, &Yld );
               }
            }
         }
      }
/*
*  Describe the resulting operand. Note that when reduction should occur, Acol
*  contains the destination column. Assuming every process column needs the
*  result, Acol is then -1.
*/
      MDescSet( DYA, M, K, Aimb, K, Amb, 1, Arow, Acol, ctxt, Yld );
   }
/*
*  End of PB_CInOutV
*/
}
