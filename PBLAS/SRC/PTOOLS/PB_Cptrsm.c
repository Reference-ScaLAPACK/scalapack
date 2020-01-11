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
void PB_Cptrsm( PBTYP_T * TYPE, Int FBCAST, char * SIDE, char * UPLO,
                char * TRANS, char * DIAG, Int M, Int N, char * ALPHA,
                char * A, Int IA, Int JA, Int * DESCA, char * BC,
                Int LDBC, char * BR, Int LDBR )
#else
void PB_Cptrsm( TYPE, FBCAST, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA,
                A, IA, JA, DESCA, BC, LDBC, BR, LDBR )
/*
*  .. Scalar Arguments ..
*/
   char           * ALPHA, * DIAG, * SIDE, * TRANS, * UPLO;
   Int            FBCAST, IA, JA, LDBC, LDBR, M, N;
   PBTYP_T        * TYPE;
/*
*  .. Array Arguments ..
*/
   Int            * DESCA;
   char           * A, * BC, * BR;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cptrsm  solves one of the matrix equations
*
*     op( sub( A ) ) * X = B,   or    X * op( sub( A ) ) = alpha * B,
*
*  where
*
*     sub( A ) denotes   A(IA:IA+M-1,JA:JA+M-1)  if SIDE = 'L',
*                        A(IA:IA+N-1,JA:JA+N-1)  if SIDE = 'R'.
*
*  X and B are m by n submatrices, sub( A ) is a unit, or non-unit,
*  upper or lower triangular submatrix and op( Y ) is one of
*
*     op( Y ) = Y   or   op( Y ) = Y'   or   op( Y ) = conjg( Y' ).
*
*  The submatrix X is overwritten on B.
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
*  SIDE    (global input) pointer to CHAR
*          On entry,  SIDE  specifies  whether op( sub( A ) ) appears on
*          the left or right of X as follows:
*
*             SIDE = 'L' or 'l'   op( sub( A ) ) * X = B,
*
*             SIDE = 'R' or 'r'   X * op( sub( A ) ) = B.
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
*             TRANS = 'N' or 'n'   sub( A )  * X = B,
*
*             TRANS = 'T' or 't'   sub( A )' * X = B,
*
*             TRANS = 'C' or 'c'   conjg( sub( A )' ) * X = B.
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
*  M       (global input) INTEGER
*          On entry, M  specifies the number of rows of the submatrix B.
*          M  must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns of the submatrix
*          B. N  must be at least zero.
*
*  A       (local input) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at  least  Lc( 0, JA+M-1 )  when  SIDE = 'L' or 'l' and is at
*          least  Lc( 0, JA+N-1 ) otherwise.  Before  entry, this  array
*          contains the local entries of the matrix A.
*          Before entry with  UPLO = 'U' or 'u', this array contains the
*          local entries corresponding to  the entries of the upper tri-
*          angular submatrix  sub( A ), and the local entries correspon-
*          ding to the entries of the strictly lower triangular part  of
*          the submatrix  sub( A )  are not referenced.
*          Before entry with  UPLO = 'L' or 'l', this array contains the
*          local entries corresponding to  the entries of the lower tri-
*          angular submatrix  sub( A ), and the local entries correspon-
*          ding to the entries of the strictly upper triangular part  of
*          the submatrix  sub( A )  are not referenced.
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
*  BC      (local input/local output) pointer to CHAR
*          On entry, BC is  an array of dimension (LDBC,Kbc), where  Kbc
*          is at least N when  SIDE is 'L' or 'l' and at least M  other-
*          wise. Before entry, when SIDE is 'L' or  'l' and TRANS is 'N'
*          or  'n'  or  SIDE is 'R' or 'r' and TRANS  is not 'N' or 'n',
*          this  array contains the local entries of the right-hand-side
*          matrix B. Otherwise, the entries of BC  should  be  zero.  On
*          exit, this  array contains the partial  solution matrix X.
*
*  LDBC    (local input) INTEGER
*          On entry,  LDBC  specifies the leading dimension of the array
*          BC. LDBC must  be  at  least MAX( 1, Lr( IA, M ) ) when  SIDE
*          is 'L' or 'l' and at  least MAX( 1, Lr( IA, N ) ) otherwise.
*
*  BR      (local input/local output) pointer to CHAR
*          On entry, BR is  an array of dimension (LDBR,Kbr), where  Kbr
*          is at least Lc( JA, M ) when  SIDE is 'L' or 'l' and at least
*          Lc( JA, N ) otherwise. Before entry, when SIDE is 'L' or  'l'
*          and  TRANS  is 'N' or 'n' or SIDE is 'R' or 'r' and TRANS  is
*          not 'N' or 'n', the entries of BR should be zero.  Otherwise,
*          this  array contains the local entries of the right-hand-side
*          matrix B. On exit, this  array contains the partial  solution
*          matrix X.
*
*  LDBR    (local input) INTEGER
*          On entry,  LDBR  specifies the leading dimension of the array
*          BR. LDBR must be at least MAX( 1, N ) when SIDE is 'L' or 'l'
*          and at least MAX( 1, M ) otherwise.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   char           btop, * negone, * one, * talpha1, * talpha2, * zero;
   Int            Acol, Aii, Aimb1, Ainb1, Ais1Col, Ais1Row, AisColRep,
                  AisRowRep, Ajj, Alcol, Ald, Alrow, Amb, Anpprev, Anb, Anp,
                  Anq, Arow, Asrc, ChangeRoc=0, LNorRT, Na, Nb, bcst, ctxt,
                  izero=0, k=0, kb, kbprev=0, kbsize, lside, mb1, mycol, myrow,
                  n1, n1last, n1p, n1pprev=0, nb1, nlast, notran, npcol, nprow,
                  rocprev, size, tmp1, tmp2;
   MMADD_T        add, tadd;
   TZPAD_T        pad;
   GEMM_T         gemm;
   TRSM_T         trsm;
   GESD2D_T       send;
   GERV2D_T       recv;
   GEBS2D_T       bsend;
   GEBR2D_T       brecv;
/*
*  .. Local Arrays ..
*/
   char           * Aprev = NULL, * Bd    = NULL, * Bdprev = NULL,
                  * Bprev = NULL, * work  = NULL;
/* ..
*  .. Executable Statements ..
*
*/
   if( ( M <= 0 ) || ( N <= 0 ) ) return;
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( ( ctxt = DESCA[CTXT_] ), &nprow, &npcol, &myrow, &mycol );

   lside  = ( Mupcase( SIDE [0] ) ==   CLEFT );
   notran = ( Mupcase( TRANS[0] ) == CNOTRAN );
   LNorRT = ( lside && notran ) || ( !( lside ) && !( notran ) );
   if( LNorRT ) { Na = M; Nb = N; } else { Na = N; Nb = M; }
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
   Aimb1   = PB_Cfirstnb( Na, IA, DESCA[IMB_], Amb );
   Anp     = PB_Cnumroc( Na, 0, Aimb1, Amb, myrow, Arow, nprow );
   Ais1Row = !( PB_Cspan( Na, 0, Aimb1, Amb, Arow, nprow ) );
   Ainb1   = PB_Cfirstnb( Na, JA, DESCA[INB_], Anb );
   Anq     = PB_Cnumroc( Na, 0, Ainb1, Anb, mycol, Acol, npcol );
   Ais1Col = !( PB_Cspan( Na, 0, Ainb1, Anb, Acol, npcol ) );
/*
*  When sub( A ) spans only one process, solve the system locally and return.
*/
   if( Ais1Row && Ais1Col )
   {
      if( LNorRT )
      {
         if( Anq > 0 )
         {
            if( Anp > 0 )
            {
               TYPE->Ftrsm( C2F_CHAR( ( notran ? SIDE : ( lside ? RIGHT :
                            LEFT ) ) ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                            C2F_CHAR( DIAG ), &M, &N, ALPHA, Mptr( A, Aii, Ajj,
                            Ald, TYPE->size ), &Ald, BC, &LDBC );
               TYPE->Fmmtadd( &M, &N, TYPE->one, BC, &LDBC, TYPE->zero, BR,
                              &LDBR );
            }
            if( ( Arow >= 0 ) && FBCAST )
            {
               btop = *PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
               if( myrow == Arow )
                  TYPE->Cgebs2d( ctxt, COLUMN, &btop, N, M, BR, LDBR );
               else
                  TYPE->Cgebr2d( ctxt, COLUMN, &btop, N, M, BR, LDBR, Arow,
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
               TYPE->Ftrsm( C2F_CHAR( ( notran ? SIDE : ( lside ? RIGHT :
                            LEFT ) ) ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                            C2F_CHAR( DIAG ), &M, &N, ALPHA, Mptr( A, Aii, Ajj,
                            Ald, TYPE->size ), &Ald, BR, &LDBR );
               TYPE->Fmmtadd( &M, &N, TYPE->one, BR, &LDBR, TYPE->zero, BC,
                              &LDBC );
            }
            if( ( Acol >= 0 ) && FBCAST )
            {
               btop = *PB_Ctop( &ctxt, BCAST, ROW, TOP_GET );
               if( mycol == Acol )
                  TYPE->Cgebs2d( ctxt, ROW, &btop, N, M, BC, LDBC );
               else
                  TYPE->Cgebr2d( ctxt, ROW, &btop, N, M, BC, LDBC, myrow,
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
   add    = TYPE->Fmmadd;  tadd  = TYPE->Fmmtadd; pad  = TYPE->Ftzpad;
   gemm   = TYPE->Fgemm;   trsm  = TYPE->Ftrsm;
   send   = TYPE->Cgesd2d; recv  = TYPE->Cgerv2d;
   bsend  = TYPE->Cgebs2d; brecv = TYPE->Cgebr2d;

   if( ( Anp > 0 ) && ( Anq > 0 ) ) A = Mptr( A, Aii, Ajj, Ald, size );

   if( LNorRT )
   {
/*
*  Left - No tran  or  Right - (co)Trans
*/
      if( ( Anq <= 0 ) || ( Ais1Row && ( ( Arow >= 0 ) && !( FBCAST ) &&
                                        ( myrow != Arow ) ) ) ) return;
      btop = *PB_Ctop( &ctxt, BCAST, COLUMN, TOP_GET );
      bcst = ( ( !Ais1Row ) || ( Ais1Row && ( Arow >= 0 ) && FBCAST ) );
      AisRowRep = ( ( Arow < 0 ) || ( nprow == 1 ) );

      if( Mupcase( UPLO[0] ) == CUPPER )
      {
/*
*  Initiate lookahead
*/
         nlast   = ( npcol - 1 ) * Anb;
         n1      = MAX( nlast, Anb );
         nlast  += Ainb1;
         n1last  = n1 - Anb + MAX( Ainb1, Anb );
         work    = PB_Cmalloc( Nb * MIN( n1last, Anp ) * size );
         tmp1    = Na-1;
         Alrow   = PB_Cindxg2p( tmp1, Aimb1, Amb, Arow, Arow, nprow );
         Alcol   = PB_Cindxg2p( tmp1, Ainb1, Anb, Acol, Acol, npcol );
         rocprev = Alcol; Anpprev = Anp; Bprev = BC; Bdprev = BR;
         Aprev   = A = Mptr( A, 0, Anq, Ald, size );
         mb1     = PB_Clastnb( Na, 0, Aimb1, Amb );
         nb1     = PB_Clastnb( Na, 0, Ainb1, Anb );
         tmp1    = Na - ( kb = MIN( mb1, nb1 ) );
         n1      = ( ( Ais1Col || ( Na - nb1 < nlast ) ) ? n1last : n1 );
         tmp2    = n1 + nb1 - kb; tmp1 -= ( tmp2 = MIN( tmp1, tmp2 ) );
         Asrc    = Arow;
         n1p     = PB_Cnumroc( tmp2, MAX( 0, tmp1 ), Aimb1, Amb, myrow, Asrc,
                               nprow );
         talpha1 = talpha2 = ( ( Ais1Col || ( mycol == Alcol ) ) ?
                               ALPHA : one );
         while( Na > 0 )
         {
            kbsize = kb * size;

            if( Ais1Col || ( mycol == Alcol ) )
            { A -= Ald*kbsize; Anq -= kb; Bd = Mptr( BR, 0, Anq, LDBR, size ); }
            if( ( Arow < 0 ) || ( myrow == Alrow ) ) { Anp -= kb; }
/*
*  Partial update of previous block
*/
            if( n1pprev > 0 )
            {
               if( ( Ais1Col || ( mycol == rocprev ) ) && ( kbprev > 0 ) )
               {
                  tmp1 = ( Anpprev - n1pprev ) * size;
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &n1pprev, &Nb,
                        &kbprev, negone, Aprev+tmp1, &Ald, Bdprev, &LDBR,
                        talpha1, Bprev+tmp1, &LDBC );
               }
/*
*  Send partial updated result to current column
*/
               if( !( Ais1Col ) && ChangeRoc )
               {
                  if( mycol == rocprev )
                  {
                     send( ctxt, n1pprev, Nb, Mptr( Bprev, Anpprev-n1pprev, 0,
                           LDBC, size ), LDBC, myrow, Alcol );
                  }
                  else if( mycol == Alcol )
                  {
                     recv( ctxt, n1pprev, Nb, work, n1pprev, myrow, rocprev );
                     add( &n1pprev, &Nb, one, work, &n1pprev, one, Mptr( Bprev,
                          Anpprev-n1pprev, 0, LDBC, size ), &LDBC );
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
                  trsm( C2F_CHAR( LEFT ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                        C2F_CHAR( DIAG ), &kb, &Nb, talpha2,  Mptr( A, Anp, 0,
                        Ald, size ), &Ald, Mptr( BC, Anp, 0, LDBC, size ),
                        &LDBC );
                  tadd( &kb, &Nb, one, Mptr( BC, Anp, 0, LDBC, size ), &LDBC,
                        zero, Mptr( BR, 0, Anq, LDBR, size ), &LDBR );
               }
               if( bcst )
               {
                  if( myrow == Alrow )
                     bsend( ctxt, COLUMN, &btop, Nb, kb, Mptr( BR, 0, Anq, LDBR,
                            size ), LDBR );
                  else
                     brecv( ctxt, COLUMN, &btop, Nb, kb, Mptr( BR, 0, Anq, LDBR,
                            size ), LDBR, Alrow, mycol );
               }
               talpha2 = one;
            }
            else
            {
               if( !( Ais1Col ) && ( AisRowRep || ( myrow == Alrow ) ) )
                  pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &kb, &Nb, &izero,
                       zero, zero, Mptr( BC, Anp, 0, LDBC, size ), &LDBC );
            }
/*
*  Finish previous update
*/
            if( ( Ais1Col || ( mycol == rocprev ) ) && ( kbprev > 0 ) )
            {
               if( ( tmp1 = Anpprev - n1pprev ) > 0 )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &tmp1, &Nb,
                        &kbprev, negone, Aprev, &Ald, Bdprev, &LDBR, talpha1,
                        Bprev, &LDBC );
               talpha1 = one;
            }
/*
*  Save info of current step and update info for the next step
*/
            if( Ais1Col   || ( mycol == Alcol ) ) { Bdprev   = Bd; Aprev = A; }
            if( AisRowRep || ( myrow == Alrow ) ) { Anpprev -= kb; }

            n1pprev = n1p;
            rocprev = Alcol;
            kbprev  = kb;
            k      += kb;
            Na     -= kb;

            mb1    -= kb;
            if( mb1 == 0 )
            {
               if( !( Ais1Row ) && ( Alrow >= 0 ) )
                  Alrow = MModSub1( Alrow, nprow );
               mb1 = ( Na > Aimb1 ? Amb : Aimb1 );
            }

            nb1      -= kb;
            ChangeRoc = ( nb1 == 0 );

            if( ChangeRoc )
            {
               if( !( Ais1Col ) && ( Alcol >= 0 ) )
                  Alcol = MModSub1( Alcol, npcol );
               nb1 = ( Na > Ainb1 ? Anb : Ainb1 );
            }
            tmp1 = Na - ( kb = MIN( mb1, nb1 ) );
            n1   = ( ( Ais1Col || ( Na-nb1 < nlast ) ) ? n1last : n1 );
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
         work  = PB_Cmalloc( Nb*MIN( n1, Anp )*size );
         Aprev = A; Bprev = BC, Bdprev = BR; Anpprev = Anp;
         mb1   = Aimb1; nb1 = Ainb1; rocprev = Acol;
         tmp1  = Na - ( kb = MIN( mb1, nb1 ) ); tmp2 = n1 + nb1 - kb;
         Asrc  = Arow;
         n1p   = PB_Cnumroc( MIN( tmp1, tmp2 ), kb, Aimb1, Amb, myrow, Asrc,
                             nprow );
         talpha1 = talpha2 = ( ( Ais1Col || ( mycol == Acol ) ) ?
                               ALPHA : one );
         while( kb > 0 )
         {
            kbsize = kb * size;
/*
*  Partial update of previous block
*/
            if( n1pprev > 0 )
            {
               if( ( Ais1Col || ( mycol == rocprev ) ) && ( kbprev > 0 ) )
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &n1pprev, &Nb,
                        &kbprev, negone, Aprev, &Ald, Bdprev, &LDBR, talpha1,
                        Bprev, &LDBC );
/*
*  Send partial updated result to current column
*/
               if( !( Ais1Col ) && ChangeRoc )
               {
                  if( mycol == rocprev )
                  {
                     send( ctxt, n1pprev, Nb, Bprev, LDBC, myrow, Acol );
                  }
                  else if( mycol == Acol )
                  {
                     recv( ctxt, n1pprev, Nb, work, n1pprev, myrow, rocprev );
                     add( &n1pprev, &Nb, one, work, &n1pprev, one, Bprev,
                          &LDBC );
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
                  trsm( C2F_CHAR( LEFT ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                        C2F_CHAR( DIAG ), &kb, &Nb, talpha2, A, &Ald, BC,
                        &LDBC );
                  tadd( &kb, &Nb, one, BC, &LDBC, zero, BR, &LDBR );
               }
               if( bcst )
               {
                  if( myrow == Arow )
                     bsend( ctxt, COLUMN, &btop, Nb, kb, BR, LDBR );
                  else
                     brecv( ctxt, COLUMN, &btop, Nb, kb, BR, LDBR, Arow,
                            mycol );
               }
               talpha2 = one;
            }
            else
            {
               if( !( Ais1Col ) && ( AisRowRep || ( myrow == Arow ) ) )
                  pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &kb, &Nb, &izero,
                       zero, zero, BC, &LDBC );
            }
/*
*  Finish previous update
*/
            if( ( Ais1Col || ( mycol == rocprev ) ) && ( kbprev > 0 ) )
            {
               if( ( tmp1 = Anpprev - n1pprev ) > 0 )
               {
                  tmp2 = n1pprev * size;
                  gemm( C2F_CHAR( NOTRAN ), C2F_CHAR( TRAN ), &tmp1, &Nb,
                        &kbprev, negone, Aprev+tmp2, &Ald, Bdprev, &LDBR,
                        talpha1, Bprev+tmp2, &LDBC );
               }
               Aprev += Ald * kbprev * size; talpha1 = one;
            }
/*
*  Save info of current step and update info for the next step
*/
            if( Ais1Col || ( mycol == Acol ) )
            { A += Ald*kbsize; Bdprev = Bd = BR; BR += LDBR*kbsize; }
            if( AisRowRep || ( myrow == Arow ) )
            {
               Bprev   = ( BC += kbsize );
               A      += kbsize;
               Aprev  += kbsize;
               Anpprev = ( Anp -= kb );
            }
            n1pprev = n1p;
            rocprev = Acol;
            kbprev  = kb;
            k      += kb;
            Na     -= kb;

            mb1    -= kb;
            if( mb1 == 0 )
            {
               if( !( Ais1Row ) && ( Arow >= 0 ) )
                  Arow = MModAdd1( Arow, nprow );
               mb1 = MIN( Amb, Na );
            }

            nb1      -= kb;
            ChangeRoc = ( nb1 == 0 );

            if( ChangeRoc )
            {
               if( !( Ais1Col ) && ( Acol >= 0 ) )
                  Acol = MModAdd1( Acol, npcol );
               nb1 = MIN( Anb, Na );
            }
            tmp1 = Na - ( kb = MIN( mb1, nb1 ) ); tmp2 = n1 + nb1 - kb;
            n1p  = PB_Cnumroc( MIN( tmp2, tmp1 ), k + kb, Aimb1, Amb, myrow,
                               Asrc, nprow );
         }
      }
   }
   else
   {
/*
*  Right - No tran  or  Left - (co)Trans
*/
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
         work  = PB_Cmalloc( Nb*MIN( n1, Anq )*size );
         Aprev = A; Bprev = BR, Bdprev = BC; Anpprev = Anq;
         mb1   = Aimb1; nb1 = Ainb1; rocprev = Arow;
         tmp1  = Na - ( kb = MIN( mb1, nb1 ) ); tmp2 = n1 + mb1 - kb;
         Asrc  = Acol;
         n1p   = PB_Cnumroc( MIN( tmp1, tmp2 ), kb, Ainb1, Anb, mycol, Asrc,
                             npcol );
         talpha1 = talpha2 = ( ( Ais1Row || ( myrow == Arow ) ) ?
                               ALPHA : one );
         while( kb > 0 )
         {
            kbsize = kb * size;
/*
*  Partial update of previous block
*/
            if( n1pprev > 0 )
            {
               if( ( Ais1Row || ( myrow == rocprev ) ) && ( kbprev > 0 ) )
                  gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &Nb, &n1pprev,
                        &kbprev, negone, Bdprev, &LDBC, Aprev, &Ald, talpha1,
                        Bprev, &LDBR );
/*
*  Send partial updated result to current row
*/
               if( !( Ais1Row ) && ChangeRoc )
               {
                  if( myrow == rocprev )
                  {
                     send( ctxt, Nb, n1pprev, Bprev, LDBR, Arow, mycol );
                  }
                  else if( myrow == Arow )
                  {
                     recv( ctxt, Nb, n1pprev, work, Nb, rocprev, mycol );
                     add( &Nb, &n1pprev, one, work, &Nb, one, Bprev, &LDBR );
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
                  trsm( C2F_CHAR( RIGHT ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                        C2F_CHAR( DIAG ), &Nb, &kb, talpha2, A, &Ald, BR,
                        &LDBR );
                  tadd( &Nb, &kb, one, BR, &LDBR, zero, BC, &LDBC );
               }
               if( bcst )
               {
                  if( mycol == Acol )
                     bsend( ctxt, ROW, &btop, kb, Nb, BC, LDBC );
                  else
                     brecv( ctxt, ROW, &btop, kb, Nb, BC, LDBC, myrow, Acol );
               }
               talpha2 = one;
            }
            else
            {
               if( !( Ais1Row ) && ( AisColRep || ( mycol == Acol ) ) )
                  pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Nb, &kb, &izero,
                       zero, zero, BR, &LDBR );
            }
/*
*  Finish previous update
*/
            if( ( Ais1Row || ( myrow == rocprev ) ) && ( kbprev > 0 ) )
            {
               if( ( tmp1 = Anpprev - n1pprev ) > 0  )
               {
                  tmp2 = n1pprev * size;
                  gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &Nb, &tmp1,
                        &kbprev, negone, Bdprev, &LDBC, Aprev+Ald*tmp2, &Ald,
                        talpha1, Bprev+LDBR*tmp2, &LDBR );
               }
               Aprev  += kbprev * size; talpha1 = one;
            }
/*
*  Save info of current step and update info for the next step
*/
            if( Ais1Row || ( myrow == Arow ) )
            { A += kbsize; Bdprev = Bd = BC; BC += kbsize; }
            if( AisColRep || ( mycol == Acol ) )
            {
               Bprev   = ( BR += LDBR * kbsize );
               A      += Ald * kbsize;
               Anpprev = ( Anq -= kb );
               Aprev  += Ald * kbsize;
            }
            n1pprev = n1p;
            rocprev = Arow;
            kbprev  = kb;
            k      += kb;
            Na     -= kb;

            nb1    -= kb;
            if( nb1 == 0 )
            {
               if( !( Ais1Col ) && ( Acol >= 0 ) )
                  Acol = MModAdd1( Acol, npcol );
               nb1 = MIN( Anb, Na );
            }

            mb1      -= kb;
            ChangeRoc = ( mb1 == 0 );

            if( ChangeRoc )
            {
               if( !( Ais1Row ) && ( Arow >= 0 ) )
                  Arow = MModAdd1( Arow, nprow );
               mb1 = MIN( Amb, Na );
            }
            tmp1 = Na - ( kb = MIN( mb1, nb1 ) ); tmp2 = n1 + mb1 - kb;
            n1p  = PB_Cnumroc( MIN( tmp2, tmp1 ), k + kb, Ainb1, Anb, mycol,
                               Asrc, npcol );
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
         work    = PB_Cmalloc( Nb * MIN( n1last, Anq ) * size );
         tmp1    = Na-1;
         Alrow   = PB_Cindxg2p( tmp1, Aimb1, Amb, Arow, Arow, nprow );
         Alcol   = PB_Cindxg2p( tmp1, Ainb1, Anb, Acol, Acol, npcol );
         rocprev = Alrow; Anpprev = Anq; Bprev = BR; Bdprev = BC;
         Aprev   = A = Mptr( A, Anp, 0, Ald, size );
         mb1     = PB_Clastnb( Na, 0, Aimb1, Amb );
         nb1     = PB_Clastnb( Na, 0, Ainb1, Anb );
         tmp1    = Na - ( kb = MIN( mb1, nb1 ) );
         n1      = ( ( Ais1Row || ( Na-mb1 < nlast ) ) ? n1last : n1 );
         tmp2    = n1 + mb1 - kb; tmp1 -= ( tmp2 = MIN( tmp1, tmp2 ) );
         Asrc    = Acol;
         n1p     = PB_Cnumroc( tmp2, MAX( 0, tmp1 ), Ainb1, Anb, mycol, Asrc,
                               npcol );
         talpha1 = talpha2 = ( ( Ais1Row || ( myrow == Alrow ) ) ?
                               ALPHA : one );
         while( Na > 0 )
         {
            kbsize = kb * size;

            if( Ais1Row || ( myrow == Alrow ) )
            { A -= kbsize; Anp -= kb; Bd = Mptr( BC, Anp, 0, LDBC, size ); }
            if( ( Acol < 0 ) || ( mycol == Alcol ) ) { Anq -= kb; }
/*
*  Partial update of previous block
*/
            if( n1pprev > 0 )
            {
               if( ( Ais1Row || ( myrow == rocprev ) ) && ( kbprev > 0 ) )
               {
                  tmp1 = ( Anpprev - n1pprev ) * size;
                  TYPE->Fgemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ),
                               &Nb, &n1pprev, &kbprev, negone, Bdprev,
                               &LDBC, Aprev+Ald*tmp1, &Ald, talpha1,
                               Bprev+LDBR*tmp1, &LDBR );
               }
/*
*  Send partial updated result to current row
*/
               if( !( Ais1Row ) && ChangeRoc )
               {
                  if( myrow == rocprev )
                  {
                     send( ctxt, Nb, n1pprev, Mptr( Bprev, 0, Anpprev-n1pprev,
                           LDBR, size ), LDBR, Alrow, mycol );
                  }
                  else if( myrow == Alrow )
                  {
                     recv( ctxt, Nb, n1pprev, work, Nb, rocprev, mycol );
                     add( &Nb, &n1pprev, one, work, &Nb, one, Mptr( Bprev, 0,
                          Anpprev-n1pprev, LDBR, size ), &LDBR );
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
                  trsm( C2F_CHAR( RIGHT ), C2F_CHAR( UPLO ), C2F_CHAR( NOTRAN ),
                        C2F_CHAR( DIAG ), &Nb, &kb, talpha2, Mptr( A, 0, Anq,
                        Ald, size ), &Ald, Mptr( BR, 0, Anq, LDBR, size ),
                        &LDBR );
                  tadd( &Nb, &kb, one, Mptr( BR, 0, Anq, LDBR, size ), &LDBR,
                        zero, Mptr( BC, Anp, 0, LDBC, size ), &LDBC );
               }
               if( bcst )
               {
                  if( mycol == Alcol )
                     bsend( ctxt, ROW, &btop, kb, Nb, Mptr( BC, Anp, 0, LDBC,
                            size ), LDBC );
                  else
                     brecv( ctxt, ROW, &btop, kb, Nb, Mptr( BC, Anp, 0, LDBC,
                            size ), LDBC, myrow, Alcol );
               }
               talpha2 = one;
            }
            else
            {
               if( !( Ais1Row ) && ( AisColRep || ( mycol == Alcol ) ) )
                  pad( C2F_CHAR( ALL ), C2F_CHAR( NOCONJG ), &Nb, &kb, &izero,
                       zero, zero, Mptr( BR, 0, Anq, LDBR, size ), &LDBR );
            }
/*
*  Finish previous update
*/
            if( ( Ais1Row || ( myrow == rocprev ) ) && ( kbprev > 0 ) )
            {
               if( ( tmp1 = Anpprev - n1pprev ) > 0 )
                  gemm( C2F_CHAR( TRAN ), C2F_CHAR( NOTRAN ), &Nb, &tmp1,
                        &kbprev, negone, Bdprev, &LDBC, Aprev, &Ald, talpha1,
                        Bprev, &LDBR );
               talpha1 = one;
            }
/*
*  Save info of current step and update info for the next step
*/
            if(  Ais1Row  || ( myrow == Alrow ) ) { Bdprev = Bd; Aprev = A; }
            if( AisColRep || ( mycol == Alcol ) ) { Anpprev -= kb; }

            n1pprev = n1p;
            rocprev = Alrow;
            kbprev  = kb;
            k      += kb;
            Na     -= kb;

            nb1    -= kb;
            if( nb1 == 0 )
            {
               if( !( Ais1Col ) && ( Alcol >= 0 ) )
                  Alcol = MModSub1( Alcol, npcol );
               nb1 = ( Na > Ainb1 ? Anb : Ainb1 );
            }

            mb1      -= kb;
            ChangeRoc = ( mb1 == 0 );

            if( ChangeRoc )
            {
               if( !( Ais1Row ) && ( Alrow >= 0 ) )
                  Alrow = MModSub1( Alrow, nprow );
               mb1 = ( Na > Aimb1 ? Amb : Aimb1 );
            }
            tmp1 = Na - ( kb = MIN( mb1, nb1 ) );
            n1   = ( ( Ais1Row || ( Na-mb1 < nlast ) ) ? n1last : n1 );
            tmp2 = n1 + mb1 - kb; tmp1 -= ( tmp2 = MIN( tmp1, tmp2 ) );
            n1p  = PB_Cnumroc( tmp2, MAX( 0, tmp1 ), Ainb1, Anb, mycol, Asrc,
                               npcol );
         }
      }
   }
   if( work ) free( work );
/*
*  End of PB_Cptrsm
*/
}
