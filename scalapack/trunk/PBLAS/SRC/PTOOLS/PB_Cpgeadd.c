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
void PB_Cpgeadd( PBTYP_T * TYPE, char * DIRECA, char * DIRECC,
                 char * CONJUG, int M, int N, char * ALPHA, char * A,
                 int IA, int JA, int * DESCA, char * BETA,  char * C,
                 int IC, int JC, int * DESCC )
#else
void PB_Cpgeadd( TYPE, DIRECA, DIRECC, CONJUG, M, N, ALPHA, A, IA, JA,
                 DESCA, BETA, C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * CONJUG, * DIRECA, * DIRECC;
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
*  PB_Cpgeadd  adds a matrix to another
*
*     sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )
*
*  where
*
*     sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),  and, op( X )  is one  of
*
*     op( X ) = X   or   op( X ) = conjg( X ).
*
*  Thus, op( sub( A ) ) denotes A(IA:IA+M-1,JA:JA+N-1)  if CONJUG = 'N',
*                        conjg(A(IA:IA+N-1,JA:JA+M-1))  if CONJUG = 'C'.
*
*  Alpha  and  beta  are scalars, sub( C ) and op( sub( A ) ) are m by n
*  submatrices.
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
*  DIRECA  (global input) pointer to CHAR
*          On entry,  DIRECA  specifies  the direction in which the rows
*          or columns of sub( A ) should be looped over as follows:
*             DIRECA = 'F' or 'f'   forward  or increasing,
*             DIRECA = 'B' or 'b'   backward or decreasing.
*
*  DIRECC  (global input) pointer to CHAR
*          On entry,  DIRECC  specifies  the direction in which the rows
*          or columns of sub( C ) should be looped over as follows:
*             DIRECC = 'F' or 'f'   forward  or increasing,
*             DIRECC = 'B' or 'b'   backward or decreasing.
*
*  CONJUG  (global input) pointer to CHAR
*          On  entry,  CONJUG  specifies  whether  conjg( sub( A ) )  or
*          sub( A ) should be added to sub( C ) as follows:
*             CONJUG = 'N' or 'n':
*                sub( C ) := beta*sub( C ) + alpha*sub( A )'
*             otherwise
*                sub( C ) := beta*sub( C ) + alpha*conjg( sub( A ) )'.
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
   char           ACroc, * one, * talpha, * tbeta, * zero;
   int            ACmyprocD, ACmyprocR, ACnD, ACnR, ACnprocsD, ACnprocsR,
                  Abufld, AcurrocR, Afr, Afwd, AiD, AiR, AiiD, AiiR, AinbD,
                  AinbR, Ainb1D, Ainb1R, AisR, Akk, Ald, AnbD, AnbR, AnpD,
                  AnpR, Aoff, ArocD, ArocR, AsrcR, Cbufld, CcurrocR, Cfr,
                  Cfwd, CiD, CiR, CiiD, CiiR, CinbD, CinbR, Cinb1D, Cinb1R,
                  CisR, Ckk, Cld, CnbD, CnbR, CnpD, CnpR, Coff, CrocD, CrocR,
                  CsrcR, ctxt, k, kb, kbb, lcmb, maxp, maxpm1, maxpq, maxq,
                  mycol, myrow, npcol, npq, nprow, ncpq, nrpq, p=0, q=0,
                  row2row, size, tmp;
   PB_VM_T        VM;
/*
*  .. Local Arrays ..
*/
   int            DBUFA[DLEN_], DBUFC[DLEN_];
   char           * Abuf = NULL, * Cbuf = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCC[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Loop over the rows of sub( C ) when M <= N, and the columns of sub( C )
*  otherwise.
*/
   row2row = ( ( M <= N ) || ( npcol == 1 ) || ( DESCA[CSRC_] == -1 ) );

   if( row2row )
   {
      AinbR = DESCA[IMB_]; AnbR = DESCA[MB_]; AsrcR = DESCA[RSRC_];
      CinbR = DESCC[IMB_]; CnbR = DESCC[MB_]; CsrcR = DESCC[RSRC_];
/*
*  If sub( A ) and sub( C ) span only one process row, then there is no need
*  to pack the data.
*/
      if( !( PB_Cspan( M, IA, AinbR, AnbR, AsrcR, nprow ) ) &&
          !( PB_Cspan( M, IC, CinbR, CnbR, CsrcR, nprow ) ) )
      {
         PB_Cpaxpby( TYPE, CONJUG, M, N, ALPHA, A, IA, JA, DESCA, ROW, BETA,
                     C, IC, JC, DESCC, ROW );
         return;
      }
/*
*  Compute local information for sub( A ) and sub( C )
*/
      ACnR      = M;           ACnD      = N;
      ACmyprocR = myrow;       ACnprocsR = nprow;
      ACmyprocD = mycol;       ACnprocsD = npcol;      ACroc = CROW;
      AiR       = IA;          AiD       = JA;
      AinbD     = DESCA[INB_]; AnbD      = DESCA[NB_]; Ald   = DESCA[LLD_];
      PB_Cinfog2l( IA, JA, DESCA, ACnprocsR, ACnprocsD, ACmyprocR, ACmyprocD,
                   &AiiR, &AiiD, &ArocR, &ArocD );
      CiR       = IC;          CiD       = JC;
      CinbD     = DESCC[INB_]; CnbD      = DESCC[NB_]; Cld = DESCC[LLD_];
      PB_Cinfog2l( IC, JC, DESCC, ACnprocsR, ACnprocsD, ACmyprocR, ACmyprocD,
                   &CiiR, &CiiD, &CrocR, &CrocD );
   }
   else
   {
      AinbR = DESCA[INB_]; AnbR = DESCA[NB_]; AsrcR = DESCA[CSRC_];
      CinbR = DESCC[INB_]; CnbR = DESCC[NB_]; CsrcR = DESCC[CSRC_];
/*
*  If sub( A ) and sub( C ) span only one process column, then there is no need
*  to pack the data.
*/
      if( !( PB_Cspan( N, JA, AinbR, AnbR, AsrcR, npcol ) ) &&
          !( PB_Cspan( N, JC, CinbR, CnbR, CsrcR, npcol ) ) )
      {
         PB_Cpaxpby( TYPE, CONJUG, M, N, ALPHA, A, IA, JA, DESCA, COLUMN, BETA,
                     C, IC, JC, DESCC, COLUMN );
         return;
      }
/*
*  Compute local information for sub( A ) and sub( C )
*/
      ACnR      = N;           ACnD      = M;
      ACmyprocR = mycol;       ACnprocsR = npcol;
      ACmyprocD = myrow;       ACnprocsD = nprow;      ACroc = CCOLUMN;
      AiR       = JA;          AiD       = IA;
      AinbD     = DESCA[IMB_]; AnbD      = DESCA[MB_]; Ald   = DESCA[LLD_];
      PB_Cinfog2l( IA, JA, DESCA, ACnprocsD, ACnprocsR, ACmyprocD, ACmyprocR,
                   &AiiD, &AiiR, &ArocD, &ArocR );
      CiR       = JC;          CiD       = IC;
      CinbD     = DESCC[IMB_]; CnbD      = DESCC[MB_]; Cld = DESCC[LLD_];
      PB_Cinfog2l( IC, JC, DESCC, ACnprocsD, ACnprocsR, ACmyprocD, ACmyprocR,
                   &CiiD, &CiiR, &CrocD, &CrocR );
   }

   size   = TYPE->size; one = TYPE->one; zero = TYPE->zero;
   kb     = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );

   Ainb1D = PB_Cfirstnb( ACnD, AiD, AinbD, AnbD );
   AnpD   = PB_Cnumroc( ACnD, 0, Ainb1D, AnbD, ACmyprocD, ArocD, ACnprocsD );
   Ainb1R = PB_Cfirstnb( ACnR, AiR, AinbR, AnbR );
   AisR   = ( ( AsrcR < 0 ) || ( ACnprocsR == 1 ) );

   Cinb1D = PB_Cfirstnb( ACnD, CiD, CinbD, CnbD );
   CnpD   = PB_Cnumroc( ACnD, 0, Cinb1D, CnbD, ACmyprocD, CrocD, ACnprocsD );
   Cinb1R = PB_Cfirstnb( ACnR, CiR, CinbR, CnbR );
   CisR   = ( ( CsrcR < 0 ) || ( ACnprocsR == 1 ) );

   lcmb   = PB_Clcm( ( maxp = ( CisR ? 1 : ACnprocsR ) ) * CnbR,
                     ( maxq = ( AisR ? 1 : ACnprocsR ) ) * AnbR );

   Afwd   = ( Mupcase( DIRECA[0] ) == CFORWARD );
   Cfwd   = ( Mupcase( DIRECC[0] ) == CFORWARD );
/*
*  When sub( A ) is not replicated and backward pass on sub( A ), find the
*  virtual process q owning the last row or column of sub( A ).
*/
   if( !( AisR ) && !( Afwd ) )
   {
      tmp = PB_Cindxg2p( ACnR-1, Ainb1R, AnbR, ArocR, ArocR, ACnprocsR );
      q   = MModSub( tmp, ArocR, ACnprocsR );
   }
/*
*  When sub( C ) is not replicated and backward pass on sub( C ), find the
*  virtual process p owning the last row or column of sub( C ).
*/
   if( !( CisR ) && !( Cfwd ) )
   {
      tmp = PB_Cindxg2p( ACnR-1, Cinb1R, CnbR, CrocR, CrocR, ACnprocsR );
      p   = MModSub( tmp, CrocR, ACnprocsR );
   }
/*
*  Loop over the processes of the virtual grid
*/
   maxpm1 = maxp - 1; maxpq = maxp * maxq;

   for( k = 0; k < maxpq; k++ )
   {
      AcurrocR = ( AisR ? -1 : MModAdd( ArocR, q, ACnprocsR ) );
      CcurrocR = ( CisR ? -1 : MModAdd( CrocR, p, ACnprocsR ) );

      if( ( AisR || ( ACmyprocR == AcurrocR ) ) ||
          ( CisR || ( ACmyprocR == CcurrocR ) ) )
      {
         Ckk = CiiR; Akk = AiiR;
/*
*  Initialize local virtual matrix in process (p,q)
*/
         AnpR = PB_Cnumroc( ACnR, 0, Ainb1R, AnbR, AcurrocR, ArocR, ACnprocsR );
         CnpR = PB_Cnumroc( ACnR, 0, Cinb1R, CnbR, CcurrocR, CrocR, ACnprocsR );
         PB_CVMinit( &VM, 0, CnpR, AnpR, Cinb1R, Ainb1R, CnbR, AnbR, p, q,
                     maxp, maxq, lcmb );
/*
*  Figure out how many diagonal entries in this new virtual process (npq).
*/
         npq = PB_CVMnpq( &VM );
/*
*  Re-adjust the number of rows or columns to be (un)packed, in order to average
*  the message sizes.
*/
         if( npq ) kbb = npq / ( ( npq - 1 ) / kb + 1 );

         if( row2row )
         {
            while( npq )
            {
               kbb = MIN( kbb, npq );
/*
*  Find out how many rows of sub( A ) and sub( C ) are contiguous
*/
               PB_CVMcontig( &VM, &nrpq, &ncpq, &Coff, &Aoff );
/*
*  Compute the descriptor DBUFA for the buffer that will contained the packed
*  rows of sub( A ).
*/
               if( ( Afr = ( ncpq < kbb ) ) != 0 )
               {
/*
*  If rows of sub( A ) are not contiguous, then allocate the buffer and pack
*  the kbb rows of sub( A ).
*/
                  Abufld = kbb;
                  if( AisR || ( ACmyprocR == AcurrocR ) )
                  {
                     Abuf = PB_Cmalloc( AnpD * kbb * size );
                     PB_CVMpack( TYPE, &VM, COLUMN, &ACroc, PACKING, NOTRAN,
                                 kbb, AnpD, one, Mptr( A, Akk, AiiD, Ald,
                                 size ), Ald, zero,  Abuf, Abufld );
                  }
               }
               else
               {
/*
*  Otherwise, re-use sub( A ) directly.
*/
                  Abufld = Ald;
                  if( AisR || ( ACmyprocR == AcurrocR ) )
                     Abuf = Mptr( A, Akk+Aoff, AiiD, Ald, size );
               }
               PB_Cdescset( DBUFA, kbb, ACnD, kbb, Ainb1D, kbb, AnbD, AcurrocR,
                            ArocD, ctxt, Abufld );
/*
*  Compute the descriptor DBUFC for the buffer that will contained the packed
*  rows of sub( C ). Allocate it.
*/
               if( ( Cfr = ( nrpq < kbb ) ) != 0 )
               {
/*
*  If rows of sub( C ) are not contiguous, then allocate receiving buffer.
*/
                  Cbufld = kbb; talpha = one;   tbeta = zero;
                  if( CisR || ( ACmyprocR == CcurrocR ) )
                     Cbuf = PB_Cmalloc( CnpD * kbb * size );
               }
               else
               {
/*
*  Otherwise, re-use sub( C ) directly.
*/
                  Cbufld = Cld; talpha = ALPHA; tbeta = BETA;
                  if( CisR || ( ACmyprocR == CcurrocR ) )
                     Cbuf = Mptr( C, Ckk+Coff, CiiD, Cld, size );
               }
               PB_Cdescset( DBUFC, kbb, ACnD, kbb, Cinb1D, kbb, CnbD, CcurrocR,
                            CrocD, ctxt, Cbufld );
/*
*  Add the one-dimensional buffer Abuf into Cbuf.
*/
               PB_Cpaxpby( TYPE, CONJUG, kbb, ACnD, talpha, Abuf, 0, 0, DBUFA,
                           &ACroc, tbeta, Cbuf, 0, 0, DBUFC, &ACroc );
/*
*  Release the buffer containing the packed rows of sub( A )
*/
               if( Afr && ( AisR || ( ACmyprocR == AcurrocR ) ) )
                  if( Abuf ) free( Abuf );
/*
*  Unpack the kbb rows of sub( C ) and release the buffer containing them.
*/
               if( Cfr && ( CisR || ( ACmyprocR == CcurrocR ) ) )
               {
                  PB_CVMpack( TYPE, &VM, ROW, &ACroc, UNPACKING, NOTRAN, kbb,
                              CnpD, BETA, Mptr( C, Ckk, CiiD, Cld, size ), Cld,
                              ALPHA, Cbuf, Cbufld );
                  if( Cbuf ) free( Cbuf );
               }
/*
*  Update the local row indexes of sub( A ) and sub( C )
*/
               PB_CVMupdate( &VM, kbb, &Ckk, &Akk );
               npq -= kbb;
            }
         }
         else
         {
            while( npq )
            {
               kbb = MIN( kbb, npq );
/*
*  Find out how many columns of sub( A ) and sub( C ) are contiguous
*/
               PB_CVMcontig( &VM, &nrpq, &ncpq, &Coff, &Aoff );
/*
*  Compute the descriptor DBUFA for the buffer that will contained the packed
*  columns of sub( A ).
*/
               if( ( Afr = ( ncpq < kbb ) ) != 0 )
               {
/*
*  If columns of sub( A ) are not contiguous, then allocate the buffer and
*  pack the kbb columns of sub( A ).
*/
                  Abufld = MAX( 1, AnpD );
                  if( AisR || ( ACmyprocR == AcurrocR ) )
                  {
                     Abuf = PB_Cmalloc( AnpD * kbb * size );
                     PB_CVMpack( TYPE, &VM, COLUMN, &ACroc, PACKING, NOTRAN,
                                 kbb, AnpD, one, Mptr( A, AiiD, Akk, Ald,
                                 size ), Ald, zero,  Abuf, Abufld );
                  }
               }
               else
               {
/*
*  Otherwise, re-use sub( A ) directly.
*/
                  Abufld = Ald;
                  if( AisR || ( ACmyprocR == AcurrocR ) )
                     Abuf = Mptr( A, AiiD, Akk+Aoff, Ald, size );
               }
               PB_Cdescset( DBUFA, ACnD, kbb, Ainb1D, kbb, AnbD, kbb, ArocD,
                            AcurrocR, ctxt, Abufld );
/*
*  Compute the descriptor DBUFC for the buffer that will contained the packed
*  columns of sub( C ). Allocate it.
*/
               if( ( Cfr = ( nrpq < kbb ) ) != 0 )
               {
/*
*  If columns of sub( C ) are not contiguous, then allocate receiving buffer.
*/
                  Cbufld = MAX( 1, CnpD ); talpha = one;   tbeta = zero;
                  if( CisR || ( ACmyprocR == CcurrocR ) )
                     Cbuf = PB_Cmalloc( CnpD * kbb * size );
               }
               else
               {
                  Cbufld = Cld;            talpha = ALPHA; tbeta = BETA;
                  if( CisR || ( ACmyprocR == CcurrocR ) )
                     Cbuf = Mptr( C, CiiD, Ckk+Coff, Cld, size );
               }
               PB_Cdescset( DBUFC, ACnD, kbb, Cinb1D, kbb, CnbD, kbb, CrocD,
                            CcurrocR, ctxt, Cbufld );
/*
*  Add the one-dimensional buffer Abuf into Cbuf.
*/
               PB_Cpaxpby( TYPE, CONJUG, ACnD, kbb, talpha, Abuf, 0, 0, DBUFA,
                           &ACroc, tbeta, Cbuf, 0, 0, DBUFC, &ACroc );
/*
*  Release the buffer containing the packed columns of sub( A )
*/
               if( Afr && ( AisR || ( ACmyprocR == AcurrocR ) ) )
                  if( Abuf ) free( Abuf );
/*
*  Unpack the kbb columns of sub( C ) and release the buffer containing them.
*/
               if( Cfr && ( CisR || ( ACmyprocR == CcurrocR ) ) )
               {
                  PB_CVMpack( TYPE, &VM, ROW, &ACroc, UNPACKING, NOTRAN, kbb,
                              CnpD, BETA, Mptr( C, CiiD, Ckk, Cld, size ), Cld,
                              ALPHA, Cbuf, Cbufld );
                  if( Cbuf ) free( Cbuf );
               }
/*
*  Update the local row index of sub( A ) and the local column index of sub( C )
*/
               PB_CVMupdate( &VM, kbb, &Ckk, &Akk );
               npq -= kbb;
            }
         }
      }
/*
*  Go to the next virtual process (p,q)
*/
      if( ( Cfwd && ( p == maxpm1 ) ) || ( !( Cfwd ) && ( p == 0 ) ) )
         q = ( Afwd ? MModAdd1( q, maxq ) : MModSub1( q, maxq ) );
      p = ( Cfwd ? MModAdd1( p, maxp ) : MModSub1( p, maxp ) );
   }
/*
*  End of PB_Cpgeadd
*/
}
