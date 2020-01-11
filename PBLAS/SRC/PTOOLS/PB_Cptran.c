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
void PB_Cptran( PBTYP_T * TYPE, char * CONJUG, Int M, Int N,
                char * ALPHA, char * A, Int IA, Int JA, Int * DESCA,
                char * BETA,  char * C, Int IC, Int JC, Int * DESCC )
#else
void PB_Cptran( TYPE, CONJUG, M, N, ALPHA, A, IA, JA, DESCA, BETA,
                C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * CONJUG;
   Int            IA, IC, JA, JC, M, N;
   char           * ALPHA, * BETA;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA, * DESCC;
   char           * A, * C;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cptran  transposes a matrix
*
*     sub( C ) := beta*sub( C ) + alpha*op( sub( A ) )
*
*  where
*
*     sub( C ) denotes C(IC:IC+M-1,JC:JC+N-1),
*
*     sub( A ) denotes A(IA:IA+N-1,JA:JA+M-1), and,
*
*     op( X ) = X' or op( X ) = conjg( X )'.
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
*  TYPE    (local input) pointer to a PBTYP_T structure
*          On entry,  TYPE  is a pointer to a structure of type PBTYP_T,
*          that contains type information (See pblas.h).
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
*          On entry,  M  specifies the number of rows of  the  submatrix
*          sub( C ) and the number of columns of the submatrix sub( A ).
*          M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          sub( C ) and the number of rows of the submatrix sub( A ).  N
*          must be at least zero.
*
*  ALPHA   (global input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of  the array  A
*          corresponding to the entries of the submatrix  sub( A )  need
*          not be set on input.
*
*  A       (local input) pointer to CHAR
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
   char           Aroc, Croc, * one, * talpha, * tbeta, * zero;
   Int            ACnD, ACnR, Abufld, AcurrocR, Afr, AiD, AiR, AiiD, AiiR,
                  AinbD, AinbR, Ainb1D, Ainb1R, AisR, Akk, Ald, AmyprocD,
                  AmyprocR, AnbD, AnbR, AnpD, AnpR, AnprocsD, AnprocsR, Aoff,
                  ArocD, ArocR, AsrcR, Cbufld, CcurrocR, Cfr, CiD, CiR, CiiD,
                  CiiR, CinbD, CinbR, Cinb1D, Cinb1R, CisR, Ckk, Cld, CmyprocD,
                  CmyprocR, CnbD, CnbR, CnpD, CnpR, CnprocsD, CnprocsR, Coff,
                  CrocD, CrocR, CsrcR, ctxt, col2row, gcdPQ, k, kb, kbb, l,
                  lcmPQ, lcmb, maxp, maxq, mycol, myrow, ncpq, npcol, npq,
                  nprow, nrpq, p, q, size;
   PB_VM_T        VM;
/*
*  .. Local Arrays ..
*/
   Int            DBUFA[DLEN_], DBUFC[DLEN_];
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
   col2row = ( ( M <= N ) || ( nprow == 1 ) || ( DESCA[RSRC_] == -1 ) );

   if( col2row )
   {
      AinbR = DESCA[INB_]; AnbR = DESCA[NB_]; AsrcR = DESCA[CSRC_];
      CinbR = DESCC[IMB_]; CnbR = DESCC[MB_]; CsrcR = DESCC[RSRC_];
/*
*  If sub( A ) only spans one process column and sub( C ) spans only one process
*  row, then there is no need to pack the data.
*/
      if( !( PB_Cspan( M, JA, AinbR, AnbR, AsrcR, npcol ) ) &&
          !( PB_Cspan( M, IC, CinbR, CnbR, CsrcR, nprow ) ) )
      {
         PB_Cpaxpby( TYPE, CONJUG, N, M, ALPHA, A, IA, JA, DESCA, COLUMN, BETA,
                     C, IC, JC, DESCC, ROW );
         return;
      }
/*
*  Compute local information for sub( A ) and sub( C )
*/
      ACnR     = M;                ACnD     = N;
      AmyprocD = CmyprocR = myrow; AnprocsD = CnprocsR = nprow;
      AmyprocR = CmyprocD = mycol; CnprocsD = AnprocsR = npcol;
      AiD   = IA;          AiR  = JA;         Aroc = CCOLUMN;
      AinbD = DESCA[IMB_]; AnbD = DESCA[MB_]; Ald  = DESCA[LLD_];
      PB_Cinfog2l( IA, JA, DESCA, AnprocsD, AnprocsR, AmyprocD, AmyprocR,
                   &AiiD, &AiiR, &ArocD, &ArocR );
      CiD   = JC;          CiR   = IC;        Croc = CROW;
      CinbD = DESCC[INB_]; CnbD = DESCC[NB_]; Cld = DESCC[LLD_];
      PB_Cinfog2l( IC, JC, DESCC, CnprocsR, CnprocsD, CmyprocR, CmyprocD,
                   &CiiR, &CiiD, &CrocR, &CrocD );
   }
   else
   {
      AinbR = DESCA[IMB_]; AnbR = DESCA[MB_]; AsrcR = DESCA[RSRC_];
      CinbR = DESCC[INB_]; CnbR = DESCC[NB_]; CsrcR = DESCC[CSRC_];
/*
*  If sub( A ) only spans one process row and sub( C ) spans only one process
*  column, then there is no need to pack the data.
*/
      if( !( PB_Cspan( N, IA, AinbR, AnbR, AsrcR, nprow ) ) &&
          !( PB_Cspan( N, JC, CinbR, CnbR, CsrcR, npcol ) ) )
      {
         PB_Cpaxpby( TYPE, CONJUG, N, M, ALPHA, A, IA, JA, DESCA, ROW, BETA, C,
                     IC, JC, DESCC, COLUMN );
         return;
      }
/*
*  Compute local information for sub( A ) and sub( C )
*/
      ACnD     = M;                ACnR = N;
      AmyprocR = CmyprocD = myrow; AnprocsR = CnprocsD = nprow;
      AmyprocD = CmyprocR = mycol; AnprocsD = CnprocsR = npcol;

      AiD   = JA;          AiR  = IA;         Aroc = CROW;
      AinbD = DESCA[INB_]; AnbD = DESCA[NB_]; Ald  = DESCA[LLD_];
      PB_Cinfog2l( IA, JA, DESCA, AnprocsR, AnprocsD, AmyprocR, AmyprocD,
                   &AiiR, &AiiD, &ArocR, &ArocD );
      CiD   = IC;          CiR  = JC;         Croc = CCOLUMN;
      CinbD = DESCC[IMB_]; CnbD = DESCC[MB_]; Cld  = DESCC[LLD_];
      PB_Cinfog2l( IC, JC, DESCC, CnprocsD, CnprocsR, CmyprocD, CmyprocR,
                   &CiiD, &CiiR, &CrocD, &CrocR );
   }

   size   = TYPE->size; one = TYPE->one; zero = TYPE->zero;
   kb     = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );

   Ainb1D = PB_Cfirstnb( ACnD, AiD, AinbD, AnbD );
   AnpD   = PB_Cnumroc( ACnD, 0, Ainb1D, AnbD, AmyprocD, ArocD, AnprocsD );
   Ainb1R = PB_Cfirstnb( ACnR, AiR, AinbR, AnbR );
   AisR   = ( ( AsrcR < 0 ) || ( AnprocsR == 1 ) );

   Cinb1D = PB_Cfirstnb( ACnD, CiD, CinbD, CnbD );
   CnpD   = PB_Cnumroc( ACnD, 0, Cinb1D, CnbD, CmyprocD, CrocD, CnprocsD );
   Cinb1R = PB_Cfirstnb( ACnR, CiR, CinbR, CnbR );
   CisR   = ( ( CsrcR < 0 ) || ( CnprocsR == 1 ) );

   lcmb   = PB_Clcm( ( maxp = ( CisR ? 1 : CnprocsR ) ) * CnbR,
                     ( maxq = ( AisR ? 1 : AnprocsR ) ) * AnbR );
   gcdPQ  = PB_Cgcd( maxp, maxq );
   lcmPQ  = ( maxp / gcdPQ ) * maxq;
/*
*  Loop over the processes of the virtual grid
*/
   for( k = 0; k < gcdPQ; k++ )
   {
      p = 0; q = k;

      for( l = 0; l < lcmPQ; l++ )
      {
         AcurrocR = ( AisR ? -1 : MModAdd( ArocR, q, AnprocsR ) );
         CcurrocR = ( CisR ? -1 : MModAdd( CrocR, p, CnprocsR ) );

         if( ( AisR || ( AmyprocR == AcurrocR ) ) ||
             ( CisR || ( CmyprocR == CcurrocR ) ) )
         {
            Ckk = CiiR; Akk = AiiR;
/*
*  Initialize local virtual matrix in process (p,q)
*/
            CnpR = PB_Cnumroc( ACnR, 0, Cinb1R, CnbR, CcurrocR, CrocR,
                               CnprocsR );
            AnpR = PB_Cnumroc( ACnR, 0, Ainb1R, AnbR, AcurrocR, ArocR,
                               AnprocsR );
            PB_CVMinit( &VM, 0, CnpR, AnpR, Cinb1R, Ainb1R, CnbR, AnbR, p, q,
                        maxp, maxq, lcmb );
/*
*  Find how many diagonals in this virtual process
*/
            npq = PB_CVMnpq( &VM );
/*
*  Re-adjust the number of rows or columns to be (un)packed, in order to
*  average the message sizes.
*/
            if( npq ) kbb = npq / ( ( npq - 1 ) / kb + 1 );

            if( col2row )
            {
               while( npq )
               {
                  kbb = MIN( kbb, npq );
/*
*  Find out how many columns of sub( A ) and rows of sub( C ) are contiguous
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
                     if( AisR || ( AmyprocR == AcurrocR ) )
                     {
                        Abuf = PB_Cmalloc( AnpD * kbb * size );
                        PB_CVMpack( TYPE, &VM, COLUMN, &Aroc, PACKING, NOTRAN,
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
                     if( AisR || ( AmyprocR == AcurrocR ) )
                        Abuf = Mptr( A, AiiD, Akk+Aoff, Ald, size );
                  }
                  PB_Cdescset( DBUFA, ACnD, kbb, Ainb1D, kbb, AnbD, kbb, ArocD,
                               AcurrocR, ctxt, Abufld );
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
                     if( CisR || ( CmyprocR == CcurrocR ) )
                        Cbuf = PB_Cmalloc( CnpD * kbb * size );
                  }
                  else
                  {
/*
*  Otherwise, re-use sub( C ) directly.
*/
                     Cbufld = Cld; talpha = ALPHA; tbeta = BETA;
                     if( CisR || ( CmyprocR == CcurrocR ) )
                        Cbuf = Mptr( C, Ckk+Coff, CiiD, Cld, size );
                  }
                  PB_Cdescset( DBUFC, kbb, ACnD, kbb, Cinb1D, kbb, CnbD,
                               CcurrocR, CrocD, ctxt, Cbufld );
/*
*  Transpose the one-dimensional buffer Abuf into Cbuf.
*/
                  PB_Cpaxpby( TYPE, CONJUG, ACnD, kbb, talpha, Abuf, 0, 0,
                              DBUFA, &Aroc, tbeta, Cbuf, 0, 0, DBUFC, &Croc );
/*
*  Release the buffer containing the packed columns of sub( A )
*/
                  if( Afr && ( AisR || ( AmyprocR == AcurrocR ) ) )
                     if( Abuf ) free( Abuf );
/*
*  Unpack the kbb rows of sub( C ) and release the buffer containing them.
*/
                  if( Cfr && ( CisR || ( CmyprocR == CcurrocR ) ) )
                  {
                     PB_CVMpack( TYPE, &VM, ROW,    &Croc, UNPACKING, NOTRAN,
                                 kbb, CnpD, BETA, Mptr( C, Ckk, CiiD, Cld,
                                 size ), Cld, ALPHA, Cbuf, Cbufld );
                     if( Cbuf ) free( Cbuf );
                  }
/*
*  Update the local column index of sub( A ) and the local row index of sub( C )
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
*  Find out how many rows of sub( A ) and columns of sub( C ) are contiguous
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
                     if( AisR || ( AmyprocR == AcurrocR ) )
                     {
                        Abuf = PB_Cmalloc( AnpD * kbb * size );
                        PB_CVMpack( TYPE, &VM, COLUMN, &Aroc, PACKING, NOTRAN,
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
                     if( AisR || ( AmyprocR == AcurrocR ) )
                        Abuf = Mptr( A, Akk+Aoff, AiiD, Ald, size );
                  }
                  PB_Cdescset( DBUFA, kbb, ACnD, kbb, Ainb1D, kbb, AnbD,
                               AcurrocR, ArocD, ctxt, Abufld );
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
                     if( CisR || ( CmyprocR == CcurrocR ) )
                        Cbuf = PB_Cmalloc( CnpD * kbb * size );
                  }
                  else
                  {
                     Cbufld = Cld;            talpha = ALPHA; tbeta = BETA;
                     if( CisR || ( CmyprocR == CcurrocR ) )
                        Cbuf = Mptr( C, CiiD, Ckk+Coff, Cld, size );
                  }
                  PB_Cdescset( DBUFC, ACnD, kbb, Cinb1D, kbb, CnbD, kbb, CrocD,
                               CcurrocR, ctxt, Cbufld );
/*
*  Transpose the one-dimensional buffer Abuf into Cbuf.
*/
                  PB_Cpaxpby( TYPE, CONJUG, kbb, ACnD, talpha, Abuf, 0, 0,
                              DBUFA, &Aroc, tbeta, Cbuf, 0, 0, DBUFC, &Croc );
/*
*  Release the buffer containing the packed rows of sub( A )
*/
                  if( Afr && ( AisR || ( AmyprocR == AcurrocR ) ) )
                     if( Abuf ) free( Abuf );
/*
*  Unpack the kbb columns of sub( C ) and release the buffer containing them.
*/
                  if( Cfr && ( CisR || ( CmyprocR == CcurrocR ) ) )
                  {
                     PB_CVMpack( TYPE, &VM, ROW,    &Croc, UNPACKING, NOTRAN,
                                 kbb, CnpD, BETA, Mptr( C, CiiD, Ckk, Cld,
                                 size ), Cld, ALPHA, Cbuf, Cbufld );
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
         p = MModAdd1( p, maxp );
         q = MModAdd1( q, maxq );
      }
   }
/*
*  End of PB_Cptran
*/
}
