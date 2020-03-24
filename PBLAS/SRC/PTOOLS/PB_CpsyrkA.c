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
void PB_CpsyrkA( PBTYP_T * TYPE, char * DIRECA, char * CONJUG,
                 char * UPLO, char * TRANS, Int N, Int K, char * ALPHA,
                 char * A, Int IA, Int JA, Int * DESCA, char * BETA,
                 char * C, Int IC, Int JC, Int * DESCC )
#else
void PB_CpsyrkA( TYPE, DIRECA, CONJUG, UPLO, TRANS, N, K, ALPHA, A, IA,
                 JA, DESCA, BETA, C, IC, JC, DESCC )
/*
*  .. Scalar Arguments ..
*/
   char           * CONJUG, * DIRECA, * TRANS, * UPLO;
   Int            IA, IC, JA, JC, K, N;
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
*  PB_CpsyrkA  performs one of the following symmetric or Hermitian rank
*  k operations
*
*     sub( C ) := alpha*sub( A )*sub( A )' + beta*sub( C ),
*  or
*     sub( C ) := alpha*sub( A )*conjg( sub( A )' ) + beta*sub( C ),
*  or
*     sub( C ) := alpha*sub( A )'*sub( A ) + beta*sub( C ),
*  or
*     sub( C ) := alpha*conjg( sub( A )' )*sub( A ) + beta*sub( C ),
*
*  where
*
*     sub( C ) denotes C(IC:IC+N-1,JC:JC+N-1), and,
*
*     sub( A ) denotes A(IA:IA+N-1,JA:JA+K-1)  if TRANS = 'N',
*                      A(IA:IA+K-1,JA:JA+N-1)  otherwise.
*
*  Alpha  and   beta  are  scalars,  sub( C )  is  an  n by n  symmetric
*  or Hermitian submatrix and  sub( A )  is  an  n by k submatrix in the
*  first case and a k by n submatrix in the second case.
*
*  This is the outer-product algorithm  using  the  logical  LCM  hybrid
*  and static blocking techniques. The submatrix operand sub( C )  stays
*  in place.
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
*          On entry, DIRECA  specifies  the direction in which the  rows
*          or columns of sub( A ) should be looped over as follows:
*             DIRECA = 'F' or 'f'   forward  or increasing,
*             DIRECA = 'B' or 'b'   backward or decreasing.
*
*  CONJUG  (global input) pointer to CHAR
*          On entry, CONJUG specifies whether sub( C ) is a symmetric or
*          Hermitian submatrix operand as follows:
*             CONJUG = 'N' or 'n'    sub( C ) is symmetric,
*             CONJUG = 'Z' or 'z'    sub( C ) is Hermitian.
*
*  UPLO    (global input) pointer to CHAR
*          On  entry,   UPLO  specifies  whether  the  local  pieces  of
*          the array  C  containing the  upper or lower triangular  part
*          of the submatrix  sub( C )  are to be referenced as follows:
*             UPLO = 'U' or 'u'   Only the local pieces corresponding to
*                                 the   upper  triangular  part  of  the
*                                 submatrix sub( C ) are referenced,
*             UPLO = 'L' or 'l'   Only the local pieces corresponding to
*                                 the   lower  triangular  part  of  the
*                                 submatrix sub( C ) are referenced.
*
*  TRANS   (global input) pointer to CHAR
*          On entry,  TRANS  specifies the  operation to be performed as
*          follows:
*
*             TRANS = 'N' or 'n'
*                  sub( C ) := alpha*sub( A )*sub( A )' + beta*sub( C ),
*             or
*                  sub( C ) := alpha*sub( A )*sub( A )' + beta*sub( C ),
*             or
*                  sub( C ) := alpha*sub( A )*conjg( sub( A )' ) +
*                              beta*sub( C ),
*
*             TRANS = 'T' or 't'
*                  sub( C ) := alpha*sub( A )'*sub( A ) + beta*sub( C ),
*             or
*                  sub( C ) := alpha*sub( A )'*sub( A ) + beta*sub( C ),
*
*             TRANS = 'C' or 'c'
*                  sub( C ) := alpha*sub( A )'*sub( A ) + beta*sub( C ),
*             or
*                  sub( C ) := alpha*conjg( sub( A )' )*sub( A ) +
*                              beta*sub( C ).
*
*  N       (global input) INTEGER
*          On entry,  N specifies the order of the  submatrix  sub( C ).
*          N must be at least zero.
*
*  K       (global input) INTEGER
*          On entry, with TRANS = 'N' or 'n',  K specifies the number of
*          columns  of  the submatrix  sub( A ), and with TRANS = 'T' or
*          't' or 'C' or 'c', K specifies the number of rows of the sub-
*          matrix sub( A ). K  must  be at least zero.
*
*  ALPHA   (global input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
*          supplied  as  zero  then  the  local entries of  the array  A
*          corresponding to the entries of the submatrix  sub( A )  need
*          not be set on input.
*
*  A       (local input) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least  Lc( 1, JA+K-1 ) when  TRANS = 'N' or 'n', and is at
*          least Lc( 1, JA+N-1 ) otherwise.  Before  entry,  this  array
*          contains the local entries of the matrix A.
*          Before entry with TRANS = 'N' or 'n', this array contains the
*          local entries corresponding to the entries of the n by k sub-
*          matrix sub( A ), otherwise the local entries corresponding to
*          the entries of the k by n submatrix sub( A ).
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
*          Before  entry  with  UPLO = 'U' or 'u', this  array  contains
*          the local entries corresponding to the upper triangular  part
*          of the  symmetric or Hermitian submatrix  sub( C ),  and  the
*          local entries corresponding to the  strictly lower triangular
*          of sub( C ) are not referenced. On exit, the upper triangular
*          part  of sub( C ) is overwritten by the upper triangular part
*          of the updated submatrix.
*          Before  entry  with  UPLO = 'L' or 'l', this  array  contains
*          the local entries corresponding to the lower triangular  part
*          of the  symmetric or Hermitian submatrix  sub( C ),  and  the
*          local entries corresponding to the  strictly upper triangular
*          of sub( C ) are not referenced. On exit, the lower triangular
*          part of sub( C ) is overwritten by the  lower triangular part
*          of the updated submatrix.
*          Note that the  imaginary parts  of the local entries  corres-
*          ponding to the  diagonal elements  of  sub( C )  need not  be
*          set,  they are assumed to be zero,  and on exit they are  set
*          to zero.
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
   char           * one;
   Int            AcurrocR, Afwd, AiD, AiR, AiiD, AiiR, AinbD, AinbR, Ainb1D,
                  Ainb1R, AisR, Ald, AmyprocD, AmyprocR, AnbD, AnbR, AnpR,
                  AnprocsD, AnprocsR, ArocD, ArocR, Arocs, AsrcR, Ccol, Cii,
                  Cimb1, Cinb1, Cjj, Clcmb, Cld, Clp, Clq, Cnq0, Cmb, Cmp,
                  Cmp0, Cnb, Cnq, Crow, WACfr, WACld, WACsum, WARfr, WARld,
                  WARsum, Wkbb=0, ctxt, k, kb, kbb, l, lb, ltmp, maxp, mycol,
                  myrow, notran, npcol, nprow, p=0, size, tmp, upper;
   GEMM_T         gemm;
   TZSYR_T        tzsyrk;
/*
*  .. Local Arrays ..
*/
   Int            Cd0[DLEN_], DBUFA[DLEN_], WACd0[DLEN_], WARd0[DLEN_];
   char           * Aptr = NULL, * Cptr = NULL, * WAC = NULL, * WAR  = NULL;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  sub( C ) := beta * sub( C )
*/
   PB_Cplascal( TYPE, UPLO, CONJUG, N, N, BETA, C, IC, JC, DESCC );
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCC[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   size = TYPE->size; one = TYPE->one; gemm  = TYPE->Fgemm;
   kb   = pilaenv_( &ctxt, C2F_CHAR( &TYPE->type ) );
/*
*  Compute descriptor Cd0 for sub( C )
*/
   PB_Cdescribe( N, N, IC, JC, DESCC, nprow, npcol, myrow, mycol, &Cii, &Cjj,
                 &Cld, &Cimb1, &Cinb1, &Cmb, &Cnb, &Crow, &Ccol, Cd0 );

   Cmp  = PB_Cnumroc( N, 0, Cimb1, Cmb, myrow, Crow, nprow );
   Cnq  = PB_Cnumroc( N, 0, Cinb1, Cnb, mycol, Ccol, npcol );

   if( ( Cmp > 0 ) && ( Cnq > 0 ) )
   {
      Cptr   = Mptr( C, Cii, Cjj, Cld, size );
      tzsyrk = ( ( Mupcase( CONJUG[0] ) == CNOCONJG ) ? PB_Ctzsyrk :
                                                        PB_Ctzherk );
/*
*  Computational partitioning size is computed as the product of the logical
*  value returned by pilaenv_ and 2 * lcm( nprow, npcol ).
*/
      Clcmb  = 2 * kb * PB_Clcm( ( Crow >= 0 ? nprow : 1 ),
                                 ( Ccol >= 0 ? npcol : 1 ) );
   }
/*
*  Retrieve local information for sub( A )
*/
   if( ( notran = ( Mupcase( TRANS[0] ) == CNOTRAN  ) ) != 0 )
   {
      AiR   = JA; AnprocsR = npcol; AinbR = DESCA[INB_]; AnbR = DESCA[NB_];
      AsrcR = DESCA[CSRC_];
   }
   else
   {
      AiR   = IA; AnprocsR = nprow; AinbR = DESCA[IMB_]; AnbR = DESCA[MB_];
      AsrcR = DESCA[RSRC_];
   }
/*
*  If sub( A ) only spans one process row or column, then there is no need to
*  pack the data.
*/
   if( !( PB_Cspan( K, AiR, AinbR, AnbR, AsrcR, AnprocsR ) ) )
   {
/*
*  Replicate sub( A ) in process rows and columns spanned by sub( C ): WAC, WAR
*/
      if( notran )
      {
         PB_CInV( TYPE, NOCONJG, COLUMN, N, N, Cd0, K, A, IA, JA, DESCA,
                  COLUMN, &WAC, WACd0, &WACfr );
         PB_CInV( TYPE, CONJUG,  ROW,    N, N, Cd0, K, WAC, 0, 0, WACd0,
                  COLUMN, &WAR, WARd0, &WARfr );
      }
      else
      {
         PB_CInV( TYPE, NOCONJG, ROW,    N, N, Cd0, K, A, IA, JA, DESCA,
                  ROW,    &WAR, WARd0, &WARfr );
         PB_CInV( TYPE, CONJUG,  COLUMN, N, N, Cd0, K, WAR, 0, 0, WARd0,
                  ROW,    &WAC, WACd0, &WACfr );
      }
/*
*  Perform the local update if I own some data
*/
      if( ( Cmp > 0 ) && ( Cnq > 0 ) )
      {
         WACld = WACd0[LLD_]; WARld = WARd0[LLD_];

         if( Mupcase( UPLO[0] ) == CUPPER )
         {
            for( l = 0; l < N; l += Clcmb )
            {
               lb   = N - l; lb = MIN( lb, Clcmb );
               Clp  = PB_Cnumroc( l,  0, Cimb1, Cmb, myrow, Crow, nprow );
               Clq  = PB_Cnumroc( l,  0, Cinb1, Cnb, mycol, Ccol, npcol );
               Cnq0 = PB_Cnumroc( lb, l, Cinb1, Cnb, mycol, Ccol, npcol );
               if( Clp > 0 && Cnq0 > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Clp, &Cnq0, &K,
                        ALPHA, WAC, &WACld, Mptr( WAR, 0, Clq, WARld, size ),
                        &WARld, one, Mptr( Cptr, 0, Clq, Cld, size ), &Cld );
               PB_Cpsyr( TYPE, UPPER, lb, K, ALPHA, Mptr( WAC, Clp, 0, WACld,
                         size ), WACld, Mptr( WAR, 0, Clq, WARld, size ), WARld,
                         Cptr, l, l, Cd0, tzsyrk );
            }
         }
         else
         {
            for( l = 0; l < N; l += Clcmb )
            {
               lb   = N - l; ltmp = l + ( lb = MIN( lb, Clcmb ) );
               Clp  = PB_Cnumroc( l, 0, Cimb1, Cmb, myrow, Crow, nprow );
               Clq  = PB_Cnumroc( l, 0, Cinb1, Cnb, mycol, Ccol, npcol );
               PB_Cpsyr( TYPE, LOWER, lb, K, ALPHA, Mptr( WAC, Clp, 0, WACld,
                         size ), WACld, Mptr( WAR, 0, Clq, WARld, size ), WARld,
                         Cptr, l, l, Cd0, tzsyrk );
               Clp  = PB_Cnumroc( ltmp, 0, Cimb1, Cmb, myrow, Crow, nprow );
               Cmp0 = Cmp - Clp;
               Cnq0 = PB_Cnumroc( lb,   l, Cinb1, Cnb, mycol, Ccol, npcol );
               if( Cmp0 > 0 && Cnq0 > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp0, &Cnq0,
                        &K, ALPHA, Mptr( WAC, Clp, 0, WACld, size ), &WACld,
                        Mptr( WAR, 0, Clq, WARld, size ), &WARld, one,
                        Mptr( Cptr, Clp, Clq, Cld, size ), &Cld );
            }
         }
      }

      if( WACfr ) free( WAC );
      if( WARfr ) free( WAR );

      return;
   }
/*
*  Otherwise sub( A ) spans more than one process row or columns -> LCM hybrid
*/
   Afwd  = ( Mupcase( DIRECA[0] ) == CFORWARD );
   upper = ( Mupcase( UPLO  [0] ) == CUPPER   );

   if( notran )
   {
      AiD = IA; AinbD = DESCA[IMB_]; AnbD = DESCA[MB_]; Ald = DESCA[LLD_];
      AmyprocD = myrow; AmyprocR = mycol; AnprocsD = nprow;
      PB_Cinfog2l( IA, JA, DESCA, AnprocsD, AnprocsR, AmyprocD, AmyprocR,
                   &AiiD, &AiiR, &ArocD, &ArocR );
   }
   else
   {
      AiD = JA; AinbD = DESCA[INB_]; AnbD = DESCA[NB_]; Ald = DESCA[LLD_];
      AmyprocD = mycol; AmyprocR = myrow; AnprocsD = npcol;
      PB_Cinfog2l( IA, JA, DESCA, AnprocsR, AnprocsD, AmyprocR, AmyprocD,
                   &AiiR, &AiiD, &ArocR, &ArocD );
   }
   Ainb1D = PB_Cfirstnb( N, AiD, AinbD, AnbD );
   Ainb1R = PB_Cfirstnb( K, AiR, AinbR, AnbR );
   AisR   = ( ( AsrcR < 0 ) || ( AnprocsR == 1 ) );
/*
*  When sub( A ) is not replicated and backward pass on sub( A ), find the
*  virtual process (p,p) owning the last row or column of sub( A ).
*/
   if( !( AisR ) && !( Afwd ) )
   {
      tmp = PB_Cindxg2p( K - 1, Ainb1R, AnbR, ArocR, ArocR, AnprocsR );
      p   = MModSub( tmp, ArocR, AnprocsR );
   }
/*
*  Allocate work space in process rows and columns spanned by sub( C )
*/
   PB_COutV( TYPE, COLUMN, NOINIT, N, N, Cd0, kb, &WAC, WACd0, &WACfr,
             &WACsum );
   PB_COutV( TYPE, ROW,    NOINIT, N, N, Cd0, kb, &WAR, WARd0, &WARfr,
             &WARsum );
/*
*  Loop over the virtual process grid induced by the rows or columns of sub( A )
*/
   maxp     = ( AisR ? 1 : AnprocsR );
   AcurrocR = ( AisR ? -1 : MModAdd( ArocR, p, AnprocsR ) );
   AnpR     = PB_Cnumroc( K, 0, Ainb1R, AnbR, AcurrocR, ArocR, AnprocsR );

   for( k = 0; k < K; k += kb )
   {
      kbb = K - k; kbb = MIN( kbb, kb );

      while( Wkbb != kbb )
      {
/*
*  Ensure that the current virtual process (p,p) has something to contribute to
*  the replicated buffers WAC and WAR.
*/
         while( AnpR == 0 )
         {
            p        = ( Afwd ? MModAdd1( p, maxp ) : MModSub1( p, maxp ) );
            AcurrocR = ( AisR ? -1 : MModAdd( ArocR, p, AnprocsR ) );
            AnpR     = PB_Cnumroc( K, 0, Ainb1R, AnbR, AcurrocR, ArocR,
                                   AnprocsR );
         }
/*
*  Current virtual process (p,p) has something, find out how many rows or
*  columns could be used: Arocs.
*/
         if( Wkbb == 0 ) { Arocs = ( AnpR < kbb ? AnpR : kbb ); }
         else            { Arocs = kbb - Wkbb; Arocs = MIN( Arocs, AnpR ); }
/*
*  The current virtual process (p,p) has Arocs rows or columns of sub( A )
*  to contribute, replicate the data over sub( C ).
*/
         if( notran )
         {
            if( AisR || ( AmyprocR == AcurrocR ) )
            { Aptr = Mptr( A, AiiD, AiiR, Ald, size ); AiiR += Arocs; }
            PB_Cdescset( DBUFA, N, Arocs, Ainb1D, Arocs, AnbD, Arocs,
                         ArocD, AcurrocR, ctxt, Ald );
/*
*  Replicate Arocs columns of sub( A ) in process columns spanned by sub( C )
*/
            PB_CInV2( TYPE, NOCONJG, COLUMN, N, N, Cd0, Arocs, Aptr, 0, 0,
                      DBUFA, COLUMN, WAC, Wkbb, WACd0 );
         }
         else
         {
            if( AisR || ( AmyprocR == AcurrocR ) )
            { Aptr = Mptr( A, AiiR, AiiD, Ald, size ); AiiR += Arocs; }
            PB_Cdescset( DBUFA, Arocs, N, Arocs, Ainb1D, Arocs, AnbD,
                         AcurrocR, ArocD, ctxt, Ald );
/*
*  Replicate Arocs rows of sub( A ) in process rows spanned by sub( C )
*/
            PB_CInV2( TYPE, NOCONJG, ROW,    N, N, Cd0, Arocs, Aptr, 0, 0,
                      DBUFA, ROW,    WAR, Wkbb, WARd0 );
         }
/*
*  Arocs rows or columns of sub( A ) have been replicated over sub( C ),
*  update the number of diagonals in this virtual process as well as the
*  number of rows or columns of sub( A ) that are in WAR or WAC.
*/
         AnpR -= Arocs;
         Wkbb += Arocs;
      }

      if( notran )
      {
/*
*  WAR := WAC'
*/
         PB_CInV2( TYPE, CONJUG,  ROW,    N, N, Cd0, kbb, WAC,  0, 0, WACd0,
                   COLUMN, WAR, 0, WARd0 );
      }
      else
      {
/*
*  WAC := WAR'
*/
         PB_CInV2( TYPE, CONJUG,  COLUMN, N, N, Cd0, kbb, WAR,  0, 0, WARd0,
                   ROW,    WAC, 0, WACd0 );
      }
/*
*  Perform the local update if I own some data
*/
      if( ( Cmp > 0 ) && ( Cnq > 0 ) )
      {
         WACld = WACd0[LLD_]; WARld = WARd0[LLD_];

         if( upper )
         {
            for( l = 0; l < N; l += Clcmb )
            {
               lb   = N - l; lb = MIN( lb, Clcmb );
               Clp  = PB_Cnumroc( l,  0, Cimb1, Cmb, myrow, Crow, nprow );
               Clq  = PB_Cnumroc( l,  0, Cinb1, Cnb, mycol, Ccol, npcol );
               Cnq0 = PB_Cnumroc( lb, l, Cinb1, Cnb, mycol, Ccol, npcol );
               if( Clp > 0 && Cnq0 > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Clp, &Cnq0,
                        &kbb, ALPHA, WAC, &WACld, Mptr( WAR, 0, Clq, WARld,
                        size ), &WARld, one, Mptr( Cptr, 0, Clq, Cld, size ),
                        &Cld );
               PB_Cpsyr( TYPE, UPPER, lb, kbb, ALPHA, Mptr( WAC, Clp, 0, WACld,
                         size ), WACld, Mptr( WAR, 0, Clq, WARld, size ), WARld,
                         Cptr, l, l, Cd0, tzsyrk );
            }
         }
         else
         {
            for( l = 0; l < N; l += Clcmb )
            {
               lb   = N - l; ltmp = l + ( lb = MIN( lb, Clcmb ) );
               Clp  = PB_Cnumroc( l, 0, Cimb1, Cmb, myrow, Crow, nprow );
               Clq  = PB_Cnumroc( l, 0, Cinb1, Cnb, mycol, Ccol, npcol );
               PB_Cpsyr( TYPE, LOWER, lb, kbb, ALPHA, Mptr( WAC, Clp, 0, WACld,
                         size ), WACld, Mptr( WAR, 0, Clq, WARld, size ), WARld,
                         Cptr, l, l, Cd0, tzsyrk );
               Clp  = PB_Cnumroc( ltmp, 0, Cimb1, Cmb, myrow, Crow, nprow );
               Cmp0 = Cmp - Clp;
               Cnq0 = PB_Cnumroc( lb,   l, Cinb1, Cnb, mycol, Ccol, npcol );
               if( Cmp0 > 0 && Cnq0 > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( NOTRAN ), &Cmp0, &Cnq0,
                        &kbb, ALPHA, Mptr( WAC, Clp, 0, WACld, size ), &WACld,
                        Mptr( WAR, 0, Clq, WARld, size ), &WARld, one,
                        Mptr( Cptr, Clp, Clq, Cld, size ), &Cld );
            }
         }
      }

      Wkbb = 0;
   }

   if( WACfr ) free( WAC );
   if( WARfr ) free( WAR );
/*
*  End of PB_CpsyrkA
*/
}
