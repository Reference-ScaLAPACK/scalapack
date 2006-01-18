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
void PB_Cptrsv( PBTYP_T * TYPE, int FBCAST, char * UPLO, char * TRANS,
                char * DIAG, int N, char * A, int IA, int JA,
                int * DESCA, char * XC, int INCXC, char * XR,
                int INCXR )
#else
void PB_Cptrsv( TYPE, FBCAST, UPLO, TRANS, DIAG, N, A, IA, JA, DESCA,
                XC, INCXC, XR, INCXR )
/*
*  .. Scalar Arguments ..
*/
   char           * DIAG, * TRANS, * UPLO;
   int            FBCAST, IA, INCXC, INCXR, JA, N;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   int            * DESCA;
   char           * A, * XC, * XR;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cptrsv  solves one of the systems of equations
*
*    sub( A )*X = b,  or sub( A )'*X = b,  or  conjg( sub( A )' )*X = b,
*
*  where sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1).
*
*  b and X are n element subvectors and  sub( A ) is an n by n  unit, or
*  non-unit, upper or lower triangular submatrix.
*
*  No test for  singularity  or  near-singularity  is included  in  this
*  routine. Such tests must be performed before calling this routine.
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
*  FBCAST  (global input) INTEGER
*          On entry, FBCAST specifies whether the transposed of the vec-
*          tor solution should be broadcast or not when there is a  pos-
*          sible ambiguity, i.e. when sub( A ) is just one  block.  When
*          FBCAST is zero, the solution vector is not broadcast, and the
*          the solution vector is broadcast otherwise.
*
*  UPLO    (global input) pointer to CHAR
*          On entry,  UPLO  specifies whether the submatrix  sub( A ) is
*          an upper or lower triangular submatrix as follows:
*
*             UPLO = 'U' or 'u'   sub( A ) is an upper triangular
*                                 submatrix,
*
*             UPLO = 'L' or 'l'   sub( A ) is a  lower triangular
*                                 submatrix.
*
*  TRANS   (global input) pointer to CHAR
*          On entry,  TRANS  specifies the  operation to be performed as
*          follows:
*
*             TRANS = 'N' or 'n'   sub( A )  * X = b,
*
*             TRANS = 'T' or 't'   sub( A )' * X = b,
*
*             TRANS = 'C' or 'c'   conjg( sub( A )' ) * X = b.
*
*  DIAG    (global input) pointer to CHAR
*          On entry,  DIAG  specifies  whether or not  sub( A )  is unit
*          triangular as follows:
*
*             DIAG = 'U' or 'u'  sub( A )  is  assumed to be unit trian-
*                                gular,
*
*             DIAG = 'N' or 'n'  sub( A ) is not assumed to be unit tri-
*                                angular.
*
*  N       (global input) INTEGER
*          On entry,  N specifies the order of the  submatrix  sub( A ).
*          N must be at least zero.
*
*  A       (local input) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least  Lc( 0, JA+N-1 ). Before  entry, this array contains
*          the local entries of the matrix A.
*          Before  entry  with  UPLO = 'U' or 'u', this  array  contains
*          the local entries corresponding to the upper triangular  part
*          of the  triangular submatrix  sub( A ), and the local entries
*          corresponding to the  strictly lower triangular  of  sub( A )
*          are not  referenced.
*          Before  entry  with  UPLO = 'L' or 'l', this  array  contains
*          the local entries corresponding to the lower triangular  part
*          of the  triangular submatrix  sub( A ), and the local entries
*          corresponding to the  strictly upper triangular  of  sub( A )
*          are not  referenced.
*          Note  that  when DIAG = 'U' or 'u', the local entries corres-
*          ponding to the  diagonal elements  of the submatrix  sub( A )
*          are not referenced either, but are assumed to be unity.
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
*  XC      (local input/local output) pointer to CHAR
*          On entry, XC is an array of dimension (LLD_X,Kx), where Kx is
*          at least 1 and LLD_X is at least Lr( IA, N ).  Before  entry,
*          when  TRANS is 'N' or 'n' this array  contains the local  en-
*          tries  of the right-hand-side vector b. When TRANS is not 'N'
*          or 'n', the entries of XC should be zero. On exit, this array
*          contains the partial solution vector x.
*
*  INCXC   (local input) INTEGER
*          On entry,  INCXC  specifies the local increment of the vector
*          XC.
*
*  XR      (local input/local output) pointer to CHAR
*          On entry, XR is an array of dimension (LLD_X,Kx), where Kx is
*          least Lc( JA, N ) and LLD_X at  least 1. Before  entry,  when
*          TRANS is 'N' or 'n'  the entries of XR should be zero. Other-
*          wise this array contains the local entries of the right-hand-
*          side vector b.  On exit, this array contains the partial  so-
*          lution vector x.
*
*  INCXR   (local input) INTEGER
*          On entry,  INCXR  specifies the local increment of the vector
*          XR.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           btop, * negone, * one, * zero;
   int            Acol, Aii, Aimb1, Ainb1, Ais1Col, Ais1Row, AisColRep,
                  AisRowRep, Ajj, Alcol, Ald, Alrow, Amb, Anpprev, Anb, Anp,
                  Anq, Arow, Asrc, ChangeRoc=0, bcst, ctxt, ione=1, k=0, kb,
                  kbprev=0, kbsize, mb1, mycol, myrow, n1, n1last, n1p,
                  n1pprev=0, nb1, nlast, npcol, nprow, rocprev, size, tmp1,
                  tmp2;
   AXPY_T         axpy;
   COPY_T         copy;
   VVSET_T        set;
   GEMV_T         gemv;
   TRSV_T         trsv;
   GESD2D_T       send;
   GERV2D_T       recv;
   GEBS2D_T       bsend;
   GEBR2D_T       brecv;
/*
*  .. Local Arrays ..
*/
   char           * Aprev = NULL, * Xd    = NULL, * Xdprev = NULL,
                  * Xprev = NULL, * work  = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   if( N <= 0 ) return;
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );
/*
*  Retrieve sub( A )'s local information: Aii, Ajj, Arow, Acol ...
*/
   PB_Cinfog2l( IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aii, &Ajj, &Arow,
                &Acol );
/*
*  Determine if sub( A ) spans more than one process row, and/or more than one
*  process column.
*/
   Amb     = DESCA[MB_]; Anb = DESCA[NB_]; Ald = DESCA[LLD_ ];
   Aimb1   = PB_Cfirstnb( N, IA, DESCA[IMB_], Amb );
   Anp     = PB_Cnumroc( N, 0, Aimb1, Amb, myrow, Arow, nprow );
   Ais1Row = !( PB_Cspan( N, 0, Aimb1, Amb, Arow, nprow ) );
   Ainb1   = PB_Cfirstnb( N, JA, DESCA[INB_], Anb );
   Anq     = PB_Cnumroc( N, 0, Ainb1, Anb, mycol, Acol, npcol );
   Ais1Col = !( PB_Cspan( N, 0, Ainb1, Anb, Acol, npcol ) );
/*
*  When sub( A ) spans only one process, solve the system locally and return.
*/
   if( Ais1Row && Ais1Col )
   {
      if( Mupcase( TRANS[0] ) == CNOTRAN )
      {
         if( Anq > 0 )
         {
            if( Anp > 0 )
            {
               TYPE->Ftrsv( C2F_CHAR( UPLO ), C2F_CHAR( TRANS ),
                            C2F_CHAR( DIAG ), &N, Mptr( A, Aii, Ajj, Ald,
                            TYPE->size ), &Ald, XC, &INCXC );
               TYPE->Fcopy( &Anp, XC, &INCXC, XR, &INCXR );
            }
            if( ( Arow >= 0 ) && FBCAST )
            {
               btop = *PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
               if( myrow == Arow )
                  TYPE->Cgebs2d( ctxt, COLUMN, &btop, 1, Anq, XR, INCXR );
               else
                  TYPE->Cgebr2d( ctxt, COLUMN, &btop, 1, Anq, XR, INCXR, Arow,
                                 mycol );
            }
         }
      }
      else
      {
         if( Anp > 0 )
         {
            if( Anq > 0 )
            {
               TYPE->Ftrsv( C2F_CHAR( UPLO ), C2F_CHAR( TRANS ),
                            C2F_CHAR( DIAG ), &N, Mptr( A, Aii, Ajj, Ald,
                            TYPE->size ), &Ald, XR, &INCXR );
               TYPE->Fcopy( &Anq, XR, &INCXR, XC, &INCXC );
            }
            if( Acol >= 0 && FBCAST )
            {
               btop = *PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );
               if( mycol == Acol )
                  TYPE->Cgebs2d( ctxt, ROW, &btop, Anp, 1, XC, Anp );
               else
                  TYPE->Cgebr2d( ctxt, ROW, &btop, Anp, 1, XC, Anp, myrow,
                                 Acol );
            }
         }
      }
      return;
   }
/*
*  Retrieve from TYPE structure useful BLAS and BLACS functions.
*/
   size   = TYPE->size;
   negone = TYPE->negone;  one   = TYPE->one;     zero = TYPE->zero;
   axpy   = TYPE->Faxpy;   copy  = TYPE->Fcopy;   set  = TYPE->Fset;
   gemv   = TYPE->Fgemv;   trsv  = TYPE->Ftrsv;
   send   = TYPE->Cgesd2d; recv  = TYPE->Cgerv2d;
   bsend  = TYPE->Cgebs2d; brecv = TYPE->Cgebr2d;

   if( ( Anp > 0 ) && ( Anq > 0 ) ) A = Mptr( A, Aii, Ajj, Ald, size );

   if( Mupcase( TRANS[0] ) == CNOTRAN )
   {
      if( ( Anq <= 0 ) || ( Ais1Row && ( ( Arow >= 0 ) && !( FBCAST ) &&
                                         ( myrow != Arow ) ) ) ) return;
      btop = *PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
      bcst = ( ( !Ais1Row ) || ( Ais1Row && ( Arow >= 0 ) && FBCAST ) );
      AisRowRep = ( ( Arow < 0 ) || ( nprow == 1 ) );

      if( Mupcase( UPLO[0] ) ==  CUPPER )
      {
/*
*  Initiate lookahead
*/
         nlast   = ( npcol - 1 ) * Anb;
         n1      = MAX( nlast, Anb );
         nlast  += Ainb1;
         n1last  = n1 - Anb + MAX( Ainb1, Anb );
         work    = PB_Cmalloc( MIN( n1last, Anp ) * size );
         tmp1    = N-1;
         Alrow   = PB_Cindxg2p( tmp1, Aimb1, Amb, Arow, Arow, nprow );
         Alcol   = PB_Cindxg2p( tmp1, Ainb1, Anb, Acol, Acol, npcol );
         rocprev = Alcol; Anpprev = Anp; Xprev = XC; Xdprev = XR;
         Aprev   = A = Mptr( A, 0, Anq, Ald, size );
         mb1     = PB_Clastnb( N, 0, Aimb1, Amb );
         nb1     = PB_Clastnb( N, 0, Ainb1, Anb );
         tmp1    = N - ( kb = MIN( mb1, nb1 ) );
         n1      = ( ( Ais1Col || ( N - nb1 < nlast ) ) ? n1last : n1 );
         tmp2    = n1 + nb1 - kb; tmp1 -= ( tmp2 = MIN( tmp1, tmp2 ) );
         Asrc    = Arow;
         n1p     = PB_Cnumroc( tmp2, MAX( 0, tmp1 ), Aimb1, Amb, myrow, Asrc,
                               nprow );
         while( N > 0 )
         {
            kbsize = kb * size;

            if( Ais1Col || ( mycol == Alcol ) )
            {
               A   -= Ald * kbsize;
               Anq -= kb;
               Xd   = Mptr( XR, 0, Anq, INCXR, size );
            }
            if( ( Arow < 0 ) || ( myrow == Alrow ) ) { Anp -= kb; }
/*
*  Partial update of previous block
*/
            if( n1pprev > 0 )
            {
               if( ( Ais1Col || ( mycol == rocprev ) ) && ( kbprev > 0 ) )
               {
                  tmp1 = ( Anpprev - n1pprev ) * size;
                  gemv( C2F_CHAR( TRANS ), &n1pprev, &kbprev, negone,
                        Aprev+tmp1, &Ald, Xdprev, &INCXR, one, Xprev+tmp1,
                        &INCXC );
               }
/*
*  Send partial updated result to current column
*/
               if( !( Ais1Col ) && ChangeRoc )
               {
                  if( mycol == rocprev )
                  {
                     send( ctxt, n1pprev, 1, Xprev+(Anpprev-n1pprev)*size,
                           n1pprev, myrow, Alcol );
                  }
                  else if( mycol == Alcol )
                  {
                     recv( ctxt, n1pprev, 1, work, n1pprev, myrow, rocprev );
                     axpy( &n1pprev, one, work, &ione, Mptr( Xprev,
                           Anpprev-n1pprev, 0, INCXC, size ), &INCXC );
                  }
               }
            }
/*
*  Solve current diagonal block
*/
            if( Ais1Col || ( mycol == Alcol ) )
            {
               if( AisRowRep || ( myrow == Alrow ) )
               {
                  trsv( C2F_CHAR( UPLO ), C2F_CHAR( TRANS ), C2F_CHAR( DIAG ),
                        &kb, Mptr( A, Anp, 0, Ald, size ), &Ald, Mptr( XC, Anp,
                        0, INCXC, size ), &INCXC );
                  copy( &kb, Mptr( XC, Anp, 0, INCXC, size ), &INCXC, Mptr( XR,
                        0, Anq, INCXR, size ), &INCXR );
               }
               if( bcst )
               {
                  if( myrow == Alrow )
                     bsend( ctxt, COLUMN, &btop, 1, kb, Mptr( XR, 0, Anq, INCXR,
                            size ), INCXR );
                  else
                     brecv( ctxt, COLUMN, &btop, 1, kb, Mptr( XR, 0, Anq, INCXR,
                            size ), INCXR, Alrow, mycol );
               }
            }
            else
            {
               if( !( Ais1Col ) && ( AisRowRep || ( myrow == Alrow ) ) )
                  set( &kb, zero, Mptr( XC, Anp, 0, INCXC, size ), &ione );
            }
/*
*  Finish previous update
*/
            if( ( Ais1Col || ( mycol == rocprev ) ) && ( kbprev > 0 ) &&
                ( ( tmp1 = Anpprev - n1pprev ) > 0 ) )
               gemv( C2F_CHAR( TRANS ), &tmp1, &kbprev, negone, Aprev, &Ald,
                     Xdprev, &INCXR, one, Xprev, &INCXC );
/*
*  Save info of current step and update info for the next step
*/
            if( Ais1Col   || ( mycol == Alcol ) ) { Xdprev   = Xd; Aprev = A; }
            if( AisRowRep || ( myrow == Alrow ) ) { Anpprev -= kb; }

            n1pprev = n1p;
            rocprev = Alcol;
            kbprev  = kb;
            k      += kb;
            N      -= kb;

            mb1    -= kb;
            if( mb1 == 0 )
            {
               if( !( Ais1Row ) && ( Alrow >= 0 ) )
                  Alrow = MModSub1( Alrow, nprow );
               mb1 = ( N > Aimb1 ? Amb : Aimb1 );
            }

            nb1      -= kb;
            ChangeRoc = ( nb1 == 0 );

            if( ChangeRoc )
            {
               if( !( Ais1Col ) && ( Alcol >= 0 ) )
                  Alcol = MModSub1( Alcol, npcol );
               nb1 = ( N > Ainb1 ? Anb : Ainb1 );
            }
            tmp1 = N - ( kb = MIN( mb1, nb1 ) );
            n1   = ( ( Ais1Col || ( N - nb1 < nlast ) ) ? n1last : n1 );
            tmp2 = n1 + nb1 - kb; tmp1 -= ( tmp2 = MIN( tmp1, tmp2 ) );
            n1p  = PB_Cnumroc( tmp2, MAX( 0, tmp1 ), Aimb1, Amb, myrow, Asrc,
                               nprow );
         }
      }
      else
      {
/*
*  Initiate lookahead
*/
         n1    = ( MAX( npcol, 2 ) - 1 ) * Anb;
         work  = PB_Cmalloc( MIN( n1, Anp ) * size );
         Aprev = A; Xprev = XC; Xdprev = XR; Anpprev = Anp;
         mb1   = Aimb1; nb1 = Ainb1; rocprev = Acol;
         tmp1  = N - ( kb = MIN( mb1, nb1 ) ); tmp2 = n1 + nb1 - kb;
         Asrc  = Arow;
         n1p   = PB_Cnumroc( MIN( tmp1, tmp2 ), kb, Aimb1, Amb, myrow, Asrc,
                           nprow );
         while( kb > 0 )
         {
            kbsize = kb * size;
/*
*  Partial update of previous block
*/
            if( n1pprev > 0 )
            {
               if( ( Ais1Col || ( mycol == rocprev ) ) && ( kbprev > 0 ) )
                  gemv( C2F_CHAR( TRANS ), &n1pprev, &kbprev, negone, Aprev,
                        &Ald, Xdprev, &INCXR, one, Xprev, &INCXC );
/*
*  Send partial updated result to current column
*/
               if( !( Ais1Col ) && ChangeRoc )
               {
                  if( mycol == rocprev )
                  {
                     send( ctxt, n1pprev, 1, Xprev, n1pprev, myrow, Acol );
                  }
                  else if( mycol == Acol )
                  {
                     recv( ctxt, n1pprev, 1, work, n1pprev, myrow, rocprev );
                     axpy( &n1pprev, one, work, &ione, Xprev, &INCXC );
                  }
               }
            }
/*
*  Solve current diagonal block
*/
            if( Ais1Col || ( mycol == Acol ) )
            {
               if( AisRowRep || ( myrow == Arow ) )
               {
                  trsv( C2F_CHAR( UPLO ), C2F_CHAR( TRANS ), C2F_CHAR( DIAG ),
                        &kb, A, &Ald, XC, &INCXC );
                  copy( &kb, XC, &INCXC, XR, &INCXR );
               }
               if( bcst )
               {
                  if( myrow == Arow )
                     bsend( ctxt, COLUMN, &btop, 1, kb, XR, INCXR );
                  else
                     brecv( ctxt, COLUMN, &btop, 1, kb, XR, INCXR, Arow,
                            mycol );
               }
            }
            else
            {
               if( !( Ais1Col ) && ( AisRowRep || ( myrow == Arow ) ) )
                  set( &kb, zero, XC, &INCXC );
            }
/*
*  Finish previous update
*/
            if( ( Ais1Col || ( mycol == rocprev ) ) && ( kbprev > 0 ) )
            {
               if( ( tmp1 = Anpprev - n1pprev ) > 0 )
               {
                  tmp2 = n1pprev * size;
                  gemv( C2F_CHAR( TRANS ), &tmp1, &kbprev, negone, Aprev+tmp2,
                        &Ald, Xdprev, &INCXR, one, Xprev+tmp2, &INCXC );
               }
               Aprev += Ald * kbprev * size;
            }
/*
*  Save info of current step and update info for the next step
*/
            if( Ais1Col || ( mycol == Acol ) )
            { A += Ald*kbsize; Xdprev = Xd = XR; XR += INCXR*kbsize; }
            if( AisRowRep || ( myrow == Arow ) )
            {
               Xprev   = ( XC += kbsize );
               A      += kbsize;
               Aprev  += kbsize;
               Anpprev = ( Anp -= kb );
            }
            n1pprev = n1p;
            rocprev = Acol;
            kbprev  = kb;
            k      += kb;
            N      -= kb;

            mb1    -= kb;
            if( mb1 == 0 )
            {
               if( !( Ais1Row ) && ( Arow >= 0 ) )
                  Arow = MModAdd1( Arow, nprow );
               mb1 = MIN( Amb, N );
            }

            nb1      -= kb;
            ChangeRoc = ( nb1 == 0 );

            if( ChangeRoc )
            {
               if( !( Ais1Col ) && ( Acol >= 0 ) )
                  Acol = MModAdd1( Acol, npcol );
               nb1 = MIN( Anb, N );
            }
            tmp1 = N - ( kb = MIN( mb1, nb1 ) ); tmp2 = n1 + nb1 - kb;
            n1p  = PB_Cnumroc( MIN( tmp2, tmp1 ), k+kb, Aimb1, Amb, myrow, Asrc,
                               nprow );
         }
      }
   }
   else
   {
      if( ( Anp <= 0 ) || ( Ais1Col && ( ( Acol >= 0 ) && !( FBCAST ) &&
                                         ( mycol != Acol ) ) ) ) return;
      btop = *PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );
      bcst = ( ( !Ais1Col ) || ( Ais1Col && ( Acol >= 0 ) && FBCAST ) );
      AisColRep = ( ( Acol < 0 ) || ( npcol == 1 ) );

      if( Mupcase( UPLO[0] ) == CUPPER )
      {
/*
*  Initiate lookahead
*/
         n1    = ( MAX( nprow, 2 ) - 1 ) * Amb;
         work  = PB_Cmalloc( MIN( n1, Anq ) * size );
         Aprev = A; Xprev = XR; Xdprev = XC; Anpprev = Anq;
         mb1   = Aimb1; nb1 = Ainb1; rocprev = Arow;
         tmp1  = N - ( kb = MIN( mb1, nb1 ) ); tmp2 = n1 + mb1 - kb;
         Asrc  = Acol;
         n1p   = PB_Cnumroc( MIN( tmp1, tmp2 ), kb, Ainb1, Anb, mycol, Asrc,
                             npcol );
         while( kb > 0 )
         {
            kbsize = kb * size;
/*
*  Partial update of previous block
*/
            if( n1pprev > 0 )
            {
               if( ( Ais1Row || ( myrow == rocprev ) ) && ( kbprev > 0 ) )
                  gemv( C2F_CHAR( TRANS ), &kbprev, &n1pprev, negone, Aprev,
                        &Ald, Xdprev, &INCXC, one, Xprev, &INCXR );
/*
*  Send partial updated result to current row
*/
               if( !( Ais1Row ) && ChangeRoc )
               {
                  if( myrow == rocprev )
                  {
                     send( ctxt, 1, n1pprev, Xprev, INCXR, Arow, mycol );
                  }
                  else if( myrow == Arow )
                  {
                     recv( ctxt, 1, n1pprev, work, 1, rocprev, mycol );
                     axpy( &n1pprev, one, work, &ione, Xprev, &INCXR );
                  }
               }
            }
/*
*  Solve current diagonal block
*/
            if( Ais1Row || ( myrow == Arow ) )
            {
               if( AisColRep || ( mycol == Acol ) )
               {
                  trsv( C2F_CHAR( UPLO ), C2F_CHAR( TRANS ), C2F_CHAR( DIAG ),
                        &kb, A, &Ald, XR, &INCXR );
                  copy( &kb, XR, &INCXR, XC, &INCXC );
               }
               if( bcst )
               {
                  if( mycol == Acol )
                     bsend( ctxt, ROW, &btop, kb, 1, XC, kb );
                  else
                     brecv( ctxt, ROW, &btop, kb, 1, XC, kb, myrow, Acol );
               }
            }
            else
            {
               if( !( Ais1Row ) && ( AisColRep || ( mycol == Acol ) ) )
                  set( &kb, zero, XR, &INCXR );
            }
/*
*  Finish previous update
*/
            if( ( Ais1Row || ( myrow == rocprev ) ) && ( kbprev > 0 ) )
            {
               if( ( tmp1 = Anpprev - n1pprev ) > 0  )
               {
                  tmp2 = n1pprev * size;
                  gemv( C2F_CHAR( TRANS ), &kbprev, &tmp1, negone,
                        Aprev+Ald*tmp2, &Ald, Xdprev, &INCXC, one,
                        Xprev+INCXR*tmp2, &INCXR );
               }
               Aprev += kbprev * size;
            }
/*
*  Save info of current step and update info for the next step
*/
            if( Ais1Row || ( myrow == Arow ) )
            { A += kbsize; Xdprev = Xd = XC; XC += kbsize; }
            if( AisColRep || ( mycol == Acol ) )
            {
               Xprev   = ( XR += INCXR * kbsize );
               A      += Ald * kbsize;
               Anpprev = ( Anq -= kb );
               Aprev  += Ald * kbsize;
            }
            n1pprev = n1p;
            rocprev = Arow;
            kbprev  = kb;
            k      += kb;
            N      -= kb;

            nb1    -= kb;
            if( nb1 == 0 )
            {
               if( !( Ais1Col ) && ( Acol >= 0 ) )
                  Acol = MModAdd1( Acol, npcol );
               nb1 = MIN( Anb, N );
            }

            mb1      -= kb;
            ChangeRoc = ( mb1 == 0 );

            if( ChangeRoc )
            {
               if( !( Ais1Row ) && ( Arow >= 0 ) )
                  Arow = MModAdd1( Arow, nprow );
               mb1 = MIN( Amb, N );
            }
            tmp1 = N - ( kb = MIN( mb1, nb1 ) ); tmp2 = n1 + mb1 - kb;
            n1p = PB_Cnumroc( MIN( tmp2, tmp1 ), k+kb, Ainb1, Anb, mycol, Asrc,
                              npcol );
         }
      }
      else
      {
/*
*  Initiate lookahead
*/
         nlast   = ( nprow - 1 ) * Amb;
         n1      = MAX( nlast, Amb );
         nlast  += Aimb1;
         n1last  = n1 - Amb + MAX( Aimb1, Amb );
         work    = PB_Cmalloc( MIN( n1last, Anq ) * size );
         tmp1    = N-1;
         Alrow   = PB_Cindxg2p( tmp1, Aimb1, Amb, Arow, Arow, nprow );
         Alcol   = PB_Cindxg2p( tmp1, Ainb1, Anb, Acol, Acol, npcol );
         rocprev = Alrow; Anpprev = Anq; Xprev = XR; Xdprev = XC;
         Aprev   = A = Mptr( A, Anp, 0, Ald, size );
         mb1     = PB_Clastnb( N, 0, Aimb1, Amb );
         nb1     = PB_Clastnb( N, 0, Ainb1, Anb );
         tmp1    = N - ( kb = MIN( mb1, nb1 ) );
         n1      = ( ( Ais1Row || ( N - mb1 < nlast ) ) ? n1last : n1 );
         tmp2    = n1 + mb1 - kb; tmp1 -= ( tmp2 = MIN( tmp1, tmp2 ) );
         Asrc    = Acol;
         n1p     = PB_Cnumroc( tmp2, MAX( 0, tmp1 ), Ainb1, Anb, mycol, Asrc,
                               npcol );
         while( N > 0 )
         {
            kbsize = kb * size;

            if( Ais1Row || ( myrow == Alrow ) )
            {
               A   -= kbsize;
               Anp -= kb;
               Xd   = Mptr( XC, Anp, 0, INCXC, size );
            }
            if( ( Acol < 0 ) || ( mycol == Alcol ) ) { Anq -= kb; }
/*
*  Partial update of previous block
*/
            if( n1pprev > 0 )
            {
               if( ( Ais1Row || ( myrow == rocprev ) ) && ( kbprev > 0 ) )
               {
                  tmp1 = ( Anpprev - n1pprev ) * size;
                  gemv( C2F_CHAR( TRANS ), &kbprev, &n1pprev, negone,
                        Aprev+Ald*tmp1, &Ald, Xdprev, &INCXC, one,
                        Xprev+INCXR*tmp1, &INCXR );
               }
/*
*  Send partial updated result to current row
*/
               if( !( Ais1Row ) && ChangeRoc )
               {
                  if( myrow == rocprev )
                  {
                     send( ctxt, 1, n1pprev, Mptr( Xprev, 0, Anpprev-n1pprev,
                           INCXR, size ), INCXR, Alrow, mycol );
                  }
                  else if( myrow == Alrow )
                  {
                     recv( ctxt, 1, n1pprev, work, 1, rocprev, mycol );
                     axpy( &n1pprev, one, work, &ione, Mptr( Xprev, 0,
                           Anpprev-n1pprev, INCXR, size ), &INCXR );
                  }
               }
            }
/*
*  Solve current diagonal block
*/
            if( Ais1Row || ( myrow == Alrow ) )
            {
               if( AisColRep || ( mycol == Alcol ) )
               {
                  trsv( C2F_CHAR( UPLO ), C2F_CHAR( TRANS ), C2F_CHAR( DIAG ),
                        &kb, Mptr( A, 0, Anq, Ald, size ), &Ald, Mptr( XR, 0,
                        Anq, INCXR, size ), &INCXR );
                  copy( &kb, Mptr( XR, 0, Anq, INCXR, size ), &INCXR, Mptr( XC,
                        0, Anp, INCXC, size ), &INCXC );
               }
               if( bcst )
               {
                  if( mycol == Alcol )
                     bsend( ctxt, ROW, &btop, kb, 1, Mptr( XC, 0, Anp, INCXC,
                            size ), kb );
                  else
                     brecv( ctxt, ROW, &btop, kb, 1, Mptr( XC, 0, Anp, INCXC,
                            size ), kb, myrow, Alcol );
               }
            }
            else
            {
               if( !( Ais1Row ) && ( AisColRep || ( mycol == Alcol ) ) )
                  set( &kb, zero, Mptr( XR, 0, Anq, INCXR, size ), &INCXR );
            }
/*
*  Finish previous update
*/
            if( ( Ais1Row || ( myrow == rocprev ) ) && ( kbprev > 0 ) &&
                ( ( tmp1 = Anpprev - n1pprev ) > 0 ) )
               gemv( C2F_CHAR( TRANS ), &kbprev, &tmp1, negone, Aprev, &Ald,
                     Xdprev, &INCXC, one, Xprev, &INCXR );
/*
*  Save info of current step and update info for the next step
*/
            if( Ais1Row   || ( myrow == Alrow ) ) { Xdprev   = Xd; Aprev = A; }
            if( AisColRep || ( mycol == Alcol ) ) { Anpprev -= kb; }

            n1pprev = n1p;
            rocprev = Alrow;
            kbprev  = kb;
            k      += kb;
            N      -= kb;

            nb1    -= kb;
            if( nb1 == 0 )
            {
               if( !( Ais1Col ) && ( Alcol >= 0 ) )
                  Alcol = MModSub1( Alcol, npcol );
               nb1 = ( N > Ainb1 ? Anb : Ainb1 );
            }

            mb1      -= kb;
            ChangeRoc = ( mb1 == 0 );

            if( ChangeRoc )
            {
               if( !( Ais1Row ) && ( Alrow >= 0 ) )
                  Alrow = MModSub1( Alrow, nprow );
               mb1 = ( N > Aimb1 ? Amb : Aimb1 );
            }
            tmp1 = N - ( kb = MIN( mb1, nb1 ) );
            n1   = ( ( Ais1Row || ( N - mb1 < nlast ) ) ? n1last : n1 );
            tmp2 = n1 + mb1 - kb; tmp1 -= ( tmp2 = MIN( tmp1, tmp2 ) );
            n1p = PB_Cnumroc( tmp2, MAX( 0, tmp1 ), Ainb1, Anb, mycol, Asrc,
                              npcol );
         }
      }
   }
   if( work ) free( work );
/*
*  End of PB_Cptrsv
*/
}
