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
void PB_Cptrm( PBTYP_T * TYPE, PBTYP_T * UTYP, char * SIDE, char * UPLO,
               char * TRANS, char * DIAG, int N, int K, char * ALPHA,
               char * A, int IA, int JA, int * DESCA,
               char * X, int LDX, char * Y, int LDY,
               TZTRM_T TRM )
#else
void PB_Cptrm( TYPE, UTYP, SIDE, UPLO, TRANS, DIAG, N, K, ALPHA, A,
               IA, JA, DESCA, X, LDX, Y, LDY, TRM )
/*
*  .. Scalar Arguments ..
*/
   char           * DIAG, * SIDE, * TRANS, * UPLO;
   int            IA, JA, K, LDX, LDY, N;
   char           * ALPHA;
   PBTYP_T        * TYPE, * UTYP;
   TZTRM_T        TRM;
/*
*  .. Array Arguments ..
*/
   int            * DESCA;
   char           * A, * X, * Y;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cptrm  performs a triangular matrix-matrix or matrix-vector multi-
*  plication.  In the following, sub( A )  denotes the triangular subma-
*  trix operand A( IA:IA+N-1, JA:JA+N-1 ).
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
*  UTYP    (local input) pointer to a PBTYP_T structure
*          On entry,  UTYP  is a pointer to a structure of type PBTYP_T,
*          that contains type information for the Y's (See pblas.h).
*
*  SIDE    (global input) pointer to CHAR
*          On entry,  SIDE  specifies whether  op( sub( A ) ) multiplies
*          its operand X from the left or right as follows:
*
*          SIDE = 'L' or 'l'  Y := alpha*op( sub( A ) )*X + Y,
*
*          SIDE = 'R' or 'r'  Y := alpha*X*op( sub( A ) ) + Y.
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
*             TRANS = 'N' or 'n'  Y := alpha * sub( A )  * X + Y,
*
*             TRANS = 'T' or 't'  Y := alpha * sub( A )' * X + Y,
*
*             TRANS = 'C' or 'c'  Y := alpha * sub( A )' * X + Y, or
*                                 Y := alpha * conjg(sub( A )') * X + Y.
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
*  K       (global input) INTEGER
*          On entry, K  specifies the local number of columns of the lo-
*          cal array X and the local number of rows of  the  local array
*          Y when SIDE is 'L' or 'l' and TRANS is 'N' or 'n', or SIDE is
*          'R' or 'r' and  TRANS  is 'T', 't', 'C' or 'c'.  Otherwise, K
*          specifies  the  local number of rows of the local array X and
*          the local number of columns of the local array Y. K mut be at
*          least zero.
*
*  ALPHA   (global input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.
*
*  A       (local input) pointer to CHAR
*          On entry, A is an array of dimension (LLD_A, Ka), where Ka is
*          at least Lc( 1, JA+N-1 ).  Before  entry, this array contains
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
*  X       (local input) pointer to CHAR
*          On entry,  X  is  an array of dimension (LDX,Kx), where Kx is
*          at least Lc( JA, N ) when SIDE is 'L' or 'l' and TRANS is 'N'
*          or 'n', or SIDE is 'R' or 'r' and  TRANS  is 'T', 't', 'C' or
*          'c', and K otherwise. Before  entry, this array contains  the
*          local entries of the matrix X.
*
*  LDX     (local input) INTEGER
*          On entry, LDX specifies the leading dimension of the array X.
*          LDX must be at least K when SIDE is 'L' or 'l' and  TRANS  is
*          'N' or 'n', or SIDE is 'R' or 'r' and TRANS  is 'T', 't', 'C'
*          or 'c', and max( 1, Lp( IA, N ) ) otherwise.
*
*  Y       (local input/local output) pointer to CHAR
*          On entry, Y is an array of dimension ( LDY, Ky ), where Ky is
*          at least max( 1, K ) when SIDE is 'L' or 'l' and TRANS is 'N'
*          or 'n', or SIDE is 'R' or 'r' and  TRANS  is 'T', 't', 'C' or
*          'c', and max( 1, Lc( JA, N ) ) otherwise. Before  entry, this
*          array contains the local entries of the matrix  Y.  On  exit,
*          this array contains the updated vector Y.
*
*  LDY     (local input) INTEGER
*          On entry, LDY specifies the leading dimension of the array Y.
*          LDY  must  be at least max( 1, Lr( IA, N ) ) when SIDE is 'L'
*          or 'l' and  TRANS  is 'N' or 'n', or  SIDE is 'R' or 'r' and
*          TRANS  is 'T', 't', 'C' or 'c', and max( 1, K ) otherwise.
*
*  TRM     (local input) pointer to function of type TZTRM_T
*          On entry, TRM  specifies the function performing the triangu-
*          lar multiplication of a single block.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   int            Acol, Arow, Aii, Aimb1, Ainb1, Ajj, Ald, Amp, Amb, Anb, Anq,
                  Aoffi, Aoffj, Arcol, Arrow, GoEast, GoSouth, IsColRepl,
                  IsRowRepl, Xinc, Yinc, XYii=0, XYjj=0, XYoffi=-1, XYoffj=-1,
                  XisRow, iimax, ilow, imbloc, inbloc, ioffd, ioffxy, iupp,
                  jjmax, joffd, joffxy, lcmt, lcmt00, lmbloc, lnbloc, low,
                  lower, m1, mbloc, mblkd, mblks, mycol, myrow, n1, nbloc,
                  nblkd, nblks, npcol, nprow, pmb, qnb, size, tmp1, upp, upper;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Quick return if possible
*/
   if( N <= 0 ) return;
/*
*  Retrieve process grid information
*/
   Cblacs_gridinfo( DESCA[CTXT_], &nprow, &npcol, &myrow, &mycol );
/*
*  Retrieve sub( A )'s local information: Aii, Ajj, Arow, Acol ...
*/
   PB_Cainfog2l( N, N, IA, JA, DESCA, nprow, npcol, myrow, mycol, &Aimb1,
                 &Ainb1, &Amp, &Anq, &Aii, &Ajj, &Arow, &Acol, &Arrow, &Arcol );
/*
*  Quick return if I don't own any of sub( A ) or if sub( A ) is replicated in
*  all processes.
*/
   if( ( Amp <= 0 ) || ( Anq <= 0 ) ) return;

   IsRowRepl = ( ( Arow < 0 ) || ( nprow == 1 ) );
   IsColRepl = ( ( Acol < 0 ) || ( npcol == 1 ) );
   Amb  = DESCA[ MB_ ]; Anb  = DESCA[ NB_ ]; Ald = DESCA[ LLD_ ];
   size = TYPE->size;

   if( IsRowRepl && IsColRepl )
   {
      TRM( TYPE, SIDE, UPLO, TRANS, DIAG, Amp, Anq, K, 0, ALPHA,
           Mptr( A, Aii, Ajj, Ald, size ), Ald, X, LDX, Y, LDY );
      return;
   }

   if( Mupcase( SIDE[0] ) == CLEFT )
   {
      if( Mupcase( TRANS[0] ) == CNOTRAN )
           { XisRow = 1; Xinc = LDX * size; Yinc = UTYP->size; }
      else { XisRow = 0; Xinc = size; Yinc = LDY * UTYP->size; }
   }
   else
   {
      if( Mupcase( TRANS[0] ) == CNOTRAN )
           { XisRow = 0; Xinc = size; Yinc = LDY * UTYP->size; }
      else { XisRow = 1; Xinc = LDX * size; Yinc = UTYP->size; }
   }
   upper = ( Mupcase( UPLO[0] ) == CUPPER );
   lower = ( Mupcase( UPLO[0] ) == CLOWER );
/*
*  Initialize lcmt00, mblks, nblks, imbloc, inbloc, lmbloc, lnbloc, ilow, low,
*  iupp, and upp.
*/
   PB_Cbinfo( 0, Amp, Anq, Aimb1, Ainb1, Amb, Anb, Arrow, Arcol, &lcmt00,
              &mblks, &nblks, &imbloc, &inbloc, &lmbloc, &lnbloc, &ilow, &low,
              &iupp, &upp );

   iimax = ( Aoffi = Aii - 1 ) + ( m1 = Amp );
   jjmax = ( Aoffj = Ajj - 1 ) + ( n1 = Anq );
   pmb   = ( IsRowRepl ? Amb : nprow * Amb );
   qnb   = ( IsColRepl ? Anb : npcol * Anb );
/*
*  Handle separately the first row and/or column of the LCM table. Update the
*  LCM value of the curent block lcmt00, as well as the number of rows and
*  columns mblks and nblks remaining in the LCM table.
*/
   GoSouth = ( lcmt00 > iupp );
   GoEast  = ( lcmt00 < ilow );

   if( XisRow )
   {
/*
*  Go through the table looking for blocks owning diagonal entries.
*/
      if( ( !( GoSouth ) ) && ( !( GoEast ) ) )
      {
/*
*  The upper left block owns diagonal entries lcmt00 >= ilow && lcmt00 <= iupp
*/
         TRM( TYPE, SIDE, UPLO, TRANS, DIAG, imbloc, inbloc, K, lcmt00, ALPHA,
              Mptr( A, Aii, Ajj, Ald, size ), Ald, X+XYjj*Xinc, LDX,
              Y+XYii*Yinc, LDY );
/*
*  Decide whether one should go south or east in the table: Go east if
*  the block below the current one only owns lower entries. If this block,
*  however, owns diagonals, then go south.
*/
         GoSouth = !( GoEast = ( ( lcmt00 - ( iupp - upp + pmb ) ) < ilow ) );

         if( GoSouth )
         {
/*
*  When the upper triangular part of sub( A ) should be operated with and
*  one is planning to go south in the table, it is neccessary to take care
*  of the remaining columns of these imbloc rows immediately.
*/
            if( upper && ( Anq > inbloc ) )
            {
               tmp1 = Anq - inbloc;
               TRM( TYPE, SIDE, ALL, TRANS, DIAG, imbloc, tmp1, K, 0, ALPHA,
                    Mptr( A, Aii, Ajj+inbloc, Ald, size ), Ald,
                    X+(XYjj+inbloc)*Xinc, LDX, Y+XYii*Yinc, LDY );
            }
            Aii += imbloc; XYii += imbloc; m1 -= imbloc;
         }
         else
         {
/*
*  When the lower triangular part of sub( A ) should be operated with and
*  one is planning to go east in the table, it is neccessary to take care
*  of the remaining rows of these inbloc columns immediately.
*/
            if( lower && ( Amp > imbloc ) )
            {
               tmp1 = Amp - imbloc;
               TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, inbloc, K, 0, ALPHA,
                    Mptr( A, Aii+imbloc, Ajj, Ald, size ), Ald, X+XYjj*Xinc,
                    LDX, Y+(XYii+imbloc)*Yinc, LDY );
            }
            Ajj += inbloc; XYjj += inbloc; n1 -= inbloc;
         }
      }

      if( GoSouth )
      {
/*
*  Go one step south in the LCM table. Adjust the current LCM value as well as
*  the local row indexes in A and XC.
*/
         lcmt00 -= ( iupp - upp + pmb ); mblks--;
         Aoffi  += imbloc; XYoffi += imbloc;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value as well as the local row indexes in A and XC.
*/
         while( ( mblks > 0 ) && ( lcmt00 > upp ) )
         { lcmt00 -= pmb; mblks--; Aoffi += Amb; XYoffi += Amb; }
/*
*  Operate with the upper triangular part of sub( A ) we just skipped when
*  necessary.
*/
         tmp1 = MIN( Aoffi, iimax ) - Aii + 1;
         if( upper && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, n1, K, 0, ALPHA,
                 Mptr( A, Aii, Aoffj+1, Ald, size ), Ald, X+(XYoffj+1)*Xinc,
                 LDX, Y+XYii*Yinc, LDY );
            Aii += tmp1; XYii += tmp1; m1  -= tmp1;
         }
/*
*  Return if no more row in the LCM table.
*/
         if( mblks <= 0 ) return;
/*
*  lcmt00 <= upp. The current block owns either diagonals or lower entries.
*  Save the current position in the LCM table. After this column has been
*  completely taken care of, re-start from this row and the next column of
*  the LCM table.
*/
         lcmt = lcmt00; mblkd = mblks; ioffd = Aoffi; ioffxy = XYoffi;

         mbloc = Amb;
         while( ( mblkd > 0 ) && ( lcmt >= ilow ) )
         {
/*
*  A block owning diagonals lcmt00 >= ilow && lcmt00 <= upp has been found.
*/
            if( mblkd == 1 ) mbloc = lmbloc;
            TRM( TYPE, SIDE, UPLO, TRANS, DIAG, mbloc, inbloc, K, lcmt,
                 ALPHA, Mptr( A, ioffd+1, Aoffj+1, Ald, size ), Ald,
                 X+(XYoffj+1)*Xinc, LDX, Y+(ioffxy+1)*Yinc, LDY );
            lcmt00 = lcmt;  lcmt   -= pmb;
            mblks  = mblkd; mblkd--;
            Aoffi  = ioffd; XYoffi  = ioffxy;
            ioffd += mbloc; ioffxy += mbloc;
         }
/*
*  Operate with the lower triangular part of sub( A ).
*/
         tmp1 = m1 - ioffd + Aii - 1;
         if( lower && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, inbloc, K, 0, ALPHA,
                 Mptr( A, ioffd+1, Aoffj+1, Ald, size ), Ald,
                 X+(XYoffj+1)*Xinc, LDX, Y+(ioffxy+1)*Yinc, LDY );
         }
         tmp1    = Aoffi - Aii + 1;
         m1     -= tmp1;
         n1     -= inbloc;
         lcmt00 += low - ilow + qnb;
         nblks--;
         Aoffj  += inbloc;
         XYoffj += inbloc;
/*
*  Operate with the upper triangular part of sub( A ).
*/
         if( upper && ( tmp1 > 0 ) && ( n1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, n1, K, 0, ALPHA,
                 Mptr( A, Aii, Aoffj+1, Ald, size ), Ald, X+(XYoffj+1)*Xinc,
                 LDX, Y+XYii*Yinc, LDY );
         }
         Aii  = Aoffi  + 1; Ajj  = Aoffj  + 1;
         XYii = XYoffi + 1; XYjj = XYoffj + 1;
      }
      else if( GoEast )
      {
/*
*  Go one step east in the LCM table. Adjust the current LCM value as well as
*  the local column index in A and XR.
*/
         lcmt00 += low - ilow + qnb; nblks--;
         Aoffj  += inbloc; XYoffj += inbloc;
/*
*  While there are blocks remaining that own lower entries, keep going east.
*  Adjust the current LCM value as well as the local column index in A and XR.
*/
         while( ( nblks > 0 ) && ( lcmt00 < low ) )
         { lcmt00 += qnb; nblks--; Aoffj += Anb; XYoffj += Anb; }
/*
*  Operate with the lower triangular part of sub( A ).
*/
         tmp1 = MIN( Aoffj, jjmax ) - Ajj + 1;
         if( lower && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, m1, tmp1, K, 0, ALPHA,
                 Mptr( A, Aii, Ajj, Ald, size ), Ald, X+XYjj*Xinc, LDX,
                 Y+XYii*Yinc, LDY );
            Ajj += tmp1; XYjj += tmp1; n1  -= tmp1;
         }
/*
*  Return if no more column in the LCM table.
*/
         if( nblks <= 0 ) return;
/*
*  lcmt00 >= low. The current block owns either diagonals or upper entries.
*  Save the current position in the LCM table. After this row has been
*  completely taken care of, re-start from this column and the next row of
*  the LCM table.
*/
         lcmt = lcmt00; nblkd = nblks; joffd = Aoffj; joffxy = XYoffj;

         nbloc = Anb;
         while( ( nblkd > 0 ) && ( lcmt <= iupp ) )
         {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= iupp has been found.
*/
            if( nblkd == 1 ) nbloc = lnbloc;
            TRM( TYPE, SIDE, UPLO, TRANS, DIAG, imbloc, nbloc, K, lcmt,
                 ALPHA, Mptr( A, Aii, joffd+1, Ald, size ), Ald,
                 X+(joffxy+1)*Xinc, LDX, Y+XYii*Yinc, LDY );
            lcmt00 = lcmt;  lcmt   += qnb;
            nblks  = nblkd; nblkd--;
            Aoffj  = joffd; XYoffj  = joffxy;
            joffd += nbloc; joffxy += nbloc;
         }
/*
*  Operate with the upper triangular part of sub( A ).
*/
         tmp1 = n1 - joffd + Ajj - 1;
         if( upper && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, imbloc, tmp1, K, 0, ALPHA,
                 Mptr( A, Aii, joffd+1, Ald, size ), Ald, X+(joffxy+1)*Xinc,
                 LDX, Y+XYii*Yinc, LDY );
         }
         tmp1    = Aoffj - Ajj + 1;
         m1     -= imbloc;
         n1     -= tmp1;
         lcmt00 -= ( iupp - upp + pmb );
         mblks--;
         Aoffi  += imbloc;
         XYoffi += imbloc;
/*
*  Operate with the lower triangular part of sub( A ).
*/
         if( lower && ( m1 > 0 ) && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, m1, tmp1, K, 0, ALPHA,
                 Mptr( A, Aoffi+1, Ajj, Ald, size ), Ald, X+XYjj*Xinc, LDX,
                 Y+(XYoffi+1)*Yinc, LDY );
         }
         Aii  = Aoffi  + 1; Ajj  = Aoffj  + 1;
         XYii = XYoffi + 1; XYjj = XYoffj + 1;
      }
/*
*  Loop over the remaining columns of the LCM table.
*/
      nbloc = Anb;
      while( nblks > 0 )
      {
         if( nblks == 1 ) nbloc = lnbloc;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value as well as the local row index in A and XC.
*/
         while( ( mblks > 0 ) && ( lcmt00 > upp ) )
         { lcmt00 -= pmb; mblks--; Aoffi += Amb; XYoffi += Amb; }
/*
*  Operate with the upper triangular part of sub( A ).
*/
         tmp1 = MIN( Aoffi, iimax ) - Aii + 1;
         if( upper && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, n1, K, 0, ALPHA,
                 Mptr( A, Aii, Aoffj+1, Ald, size ), Ald, X+(XYoffj+1)*Xinc,
                 LDX, Y+XYii*Yinc, LDY );
            Aii  += tmp1;
            XYii += tmp1;
            m1   -= tmp1;
         }
/*
*  Return if no more row in the LCM table.
*/
         if( mblks <= 0 ) return;
/*
*  lcmt00 <= upp. The current block owns either diagonals or lower entries.
*  Save the current position in the LCM table. After this column has been
*  completely taken care of, re-start from this row and the next column of
*  the LCM table.
*/
         lcmt = lcmt00; mblkd = mblks; ioffd = Aoffi; ioffxy = XYoffi;

         mbloc = Amb;
         while( ( mblkd > 0 ) && ( lcmt >= low ) )
         {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
*/
            if( mblkd == 1 ) mbloc = lmbloc;
            TRM( TYPE, SIDE, UPLO, TRANS, DIAG, mbloc, nbloc, K, lcmt,
                 ALPHA, Mptr( A, ioffd+1, Aoffj+1, Ald, size ), Ald,
                 X+(XYoffj+1)*Xinc, LDX, Y+(ioffxy+1)*Yinc, LDY );
            lcmt00 = lcmt;  lcmt   -= pmb;
            mblks  = mblkd; mblkd--;
            Aoffi  = ioffd; XYoffi = ioffxy;
            ioffd += mbloc; ioffxy += mbloc;
         }
/*
*  Operate with the lower triangular part of sub( A ).
*/
         tmp1 = m1 - ioffd + Aii - 1;
         if( lower && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, nbloc, K, 0, ALPHA,
                 Mptr( A, ioffd+1, Aoffj+1, Ald, size ), Ald,
                 X+(XYoffj+1)*Xinc, LDX, Y+(ioffxy+1)*Yinc, LDY );
         }

         tmp1    = MIN( Aoffi, iimax ) - Aii + 1;
         m1     -= tmp1;
         n1     -= nbloc;
         lcmt00 += qnb;
         nblks--;
         Aoffj  += nbloc;
         XYoffj += nbloc;
/*
*  Operate with the upper triangular part of sub( A ).
*/
         if( upper && ( tmp1 > 0 ) && ( n1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, n1, K, 0, ALPHA,
                 Mptr( A, Aii, Aoffj+1, Ald, size ), Ald, X+(XYoffj+1)*Xinc,
                 LDX, Y+XYii*Yinc, LDY );
         }
         Aii  = Aoffi  + 1;  Ajj = Aoffj  + 1;
         XYii = XYoffi + 1; XYjj = XYoffj + 1;
      }
   }
   else
   {
/*
*  Go through the table looking for blocks owning diagonal entries.
*/
      if( ( !( GoSouth ) ) && ( !( GoEast ) ) )
      {
/*
*  The upper left block owns diagonal entries lcmt00 >= ilow && lcmt00 <= iupp
*/
         TRM( TYPE, SIDE, UPLO, TRANS, DIAG, imbloc, inbloc, K, lcmt00, ALPHA,
              Mptr( A, Aii, Ajj, Ald, size ), Ald, X+XYii*Xinc, LDX,
              Y+XYjj*Yinc, LDY );
/*
*  Decide whether one should go south or east in the table: Go east if
*  the block below the current one only owns lower entries. If this block,
*  however, owns diagonals, then go south.
*/
         GoSouth = !( GoEast = ( ( lcmt00 - ( iupp - upp + pmb ) ) < ilow ) );

         if( GoSouth )
         {
/*
*  When the upper triangular part of sub( A ) should be operated with and
*  one is planning to go south in the table, it is neccessary to take care
*  of the remaining columns of these imbloc rows immediately.
*/
            if( upper && ( Anq > inbloc ) )
            {
               tmp1 = Anq - inbloc;
               TRM( TYPE, SIDE, ALL, TRANS, DIAG, imbloc, tmp1, K, 0, ALPHA,
                    Mptr( A, Aii, Ajj+inbloc, Ald, size ), Ald, X+XYii*Xinc,
                    LDX, Y+(XYjj+inbloc)*Yinc, LDY );
            }
            Aii += imbloc; XYii += imbloc; m1 -= imbloc;
         }
         else
         {
/*
*  When the lower triangular part of sub( A ) should be operated with and
*  one is planning to go east in the table, it is neccessary to take care
*  of the remaining rows of these inbloc columns immediately.
*/
            if( lower && ( Amp > imbloc ) )
            {
               tmp1 = Amp - imbloc;
               TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, inbloc, K, 0, ALPHA,
                    Mptr( A, Aii+imbloc, Ajj, Ald, size ), Ald,
                    X+(XYii+imbloc)*Xinc, LDX, Y+XYjj*Yinc, LDY );
            }
            Ajj += inbloc; XYjj += inbloc; n1 -= inbloc;
         }
      }

      if( GoSouth )
      {
/*
*  Go one step south in the LCM table. Adjust the current LCM value as well as
*  the local row indexes in A and XC.
*/
         lcmt00 -= ( iupp - upp + pmb ); mblks--;
         Aoffi  += imbloc; XYoffi  += imbloc;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value as well as the local row indexes in A and XC.
*/
         while( ( mblks > 0 ) && ( lcmt00 > upp ) )
         { lcmt00 -= pmb; mblks--; Aoffi += Amb; XYoffi += Amb; }
/*
*  Operate with the upper triangular part of sub( A ) we just skipped when
*  necessary.
*/
         tmp1 = MIN( Aoffi, iimax ) - Aii + 1;
         if( upper && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, n1, K, 0, ALPHA,
                 Mptr( A, Aii, Aoffj+1, Ald, size ), Ald, X+XYii*Xinc, LDX,
                 Y+(XYoffj+1)*Yinc, LDY );
            Aii += tmp1; XYii += tmp1; m1  -= tmp1;
         }
/*
*  Return if no more row in the LCM table.
*/
         if( mblks <= 0 ) return;
/*
*  lcmt00 <= upp. The current block owns either diagonals or lower entries.
*  Save the current position in the LCM table. After this column has been
*  completely taken care of, re-start from this row and the next column of
*  the LCM table.
*/
         lcmt = lcmt00; mblkd = mblks; ioffd = Aoffi; ioffxy = XYoffi;

         mbloc = Amb;
         while( ( mblkd > 0 ) && ( lcmt >= ilow ) )
         {
/*
*  A block owning diagonals lcmt00 >= ilow && lcmt00 <= upp has been found.
*/
            if( mblkd == 1 ) mbloc = lmbloc;
            TRM( TYPE, SIDE, UPLO, TRANS, DIAG, mbloc, inbloc, K, lcmt,
                 ALPHA, Mptr( A, ioffd+1, Aoffj+1, Ald, size ), Ald,
                 X+(ioffxy+1)*Xinc, LDX, Y+(XYoffj+1)*Yinc, LDY );
            lcmt00 = lcmt;  lcmt   -= pmb;
            mblks  = mblkd; mblkd--;
            Aoffi  = ioffd; XYoffi  = ioffxy;
            ioffd += mbloc; ioffxy += mbloc;
         }
/*
*  Operate with the lower triangular part of sub( A ).
*/
         tmp1 = m1 - ioffd + Aii - 1;
         if( lower && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, inbloc, K, 0, ALPHA,
                 Mptr( A, ioffd+1, Aoffj+1, Ald, size ), Ald,
                 X+(ioffxy+1)*Xinc, LDX, Y+(XYoffj+1)*Yinc, LDY );
         }
         tmp1    = Aoffi - Aii + 1;
         m1     -= tmp1;
         n1     -= inbloc;
         lcmt00 += low - ilow + qnb;
         nblks--;
         Aoffj  += inbloc;
         XYoffj += inbloc;
/*
*  Operate with the upper triangular part of sub( A ).
*/
         if( upper && ( tmp1 > 0 ) && ( n1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, n1, K, 0, ALPHA,
                 Mptr( A, Aii, Aoffj+1, Ald, size ), Ald, X+XYii*Xinc, LDX,
                 Y+(XYoffj+1)*Yinc, LDY );
         }
         Aii  = Aoffi  + 1; Ajj  = Aoffj  + 1;
         XYii = XYoffi + 1; XYjj = XYoffj + 1;
      }
      else if( GoEast )
      {
/*
*  Go one step east in the LCM table. Adjust the current LCM value as well as
*  the local column index in A and XR.
*/
         lcmt00 += low - ilow + qnb; nblks--;
         Aoffj  += inbloc; XYoffj += inbloc;
/*
*  While there are blocks remaining that own lower entries, keep going east.
*  Adjust the current LCM value as well as the local column index in A and XR.
*/
         while( ( nblks > 0 ) && ( lcmt00 < low ) )
         { lcmt00 += qnb; nblks--; Aoffj += Anb; XYoffj += Anb; }
/*
*  Operate with the lower triangular part of sub( A ).
*/
         tmp1 = MIN( Aoffj, jjmax ) - Ajj + 1;
         if( lower && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, m1, tmp1, K, 0, ALPHA,
                 Mptr( A, Aii, Ajj, Ald, size ), Ald, X+XYii*Xinc, LDX,
                 Y+XYjj*Yinc, LDY );
            Ajj += tmp1; XYjj += tmp1; n1  -= tmp1;
         }
/*
*  Return if no more column in the LCM table.
*/
         if( nblks <= 0 ) return;
/*
*  lcmt00 >= low. The current block owns either diagonals or upper entries.
*  Save the current position in the LCM table. After this row has been
*  completely taken care of, re-start from this column and the next row of
*  the LCM table.
*/
         lcmt = lcmt00; nblkd = nblks; joffd = Aoffj; joffxy = XYoffj;

         nbloc = Anb;
         while( ( nblkd > 0 ) && ( lcmt <= iupp ) )
         {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= iupp has been found.
*/
            if( nblkd == 1 ) nbloc = lnbloc;
            TRM( TYPE, SIDE, UPLO, TRANS, DIAG, imbloc, nbloc, K, lcmt,
                 ALPHA, Mptr( A, Aii, joffd+1, Ald, size ), Ald, X+XYii*Xinc,
                 LDX, Y+(joffxy+1)*Yinc, LDY );
            lcmt00 = lcmt;  lcmt   += qnb;
            nblks  = nblkd; nblkd--;
            Aoffj  = joffd; XYoffj  = joffxy;
            joffd += nbloc; joffxy += nbloc;
         }
/*
*  Operate with the upper triangular part of sub( A ).
*/
         tmp1 = n1 - joffd + Ajj - 1;
         if( upper && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, imbloc, tmp1, K, 0, ALPHA,
                 Mptr( A, Aii, joffd+1, Ald, size ), Ald, X+XYii*Xinc, LDX,
                 Y+(joffxy+1)*Yinc, LDY );
         }
         tmp1    = Aoffj - Ajj + 1;
         m1     -= imbloc;
         n1     -= tmp1;
         lcmt00 -= ( iupp - upp + pmb );
         mblks--;
         Aoffi  += imbloc;
         XYoffi += imbloc;
/*
*  Operate with the lower triangular part of sub( A ).
*/
         if( lower && ( m1 > 0 ) && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, m1, tmp1, K, 0, ALPHA,
                 Mptr( A, Aoffi+1, Ajj, Ald, size ), Ald, X+(XYoffi+1)*Xinc,
                 LDX, Y+XYjj*Yinc, LDY );
         }
         Aii  = Aoffi  + 1; Ajj  = Aoffj  + 1;
         XYii = XYoffi + 1; XYjj = XYoffj + 1;
      }
/*
*  Loop over the remaining columns of the LCM table.
*/
      nbloc = Anb;
      while( nblks > 0 )
      {
         if( nblks == 1 ) nbloc = lnbloc;
/*
*  While there are blocks remaining that own upper entries, keep going south.
*  Adjust the current LCM value as well as the local row index in A and XC.
*/
         while( ( mblks > 0 ) && ( lcmt00 > upp ) )
         { lcmt00 -= pmb; mblks--; Aoffi += Amb; XYoffi += Amb; }
/*
*  Operate with the upper triangular part of sub( A ).
*/
         tmp1 = MIN( Aoffi, iimax ) - Aii + 1;
         if( upper && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, n1, K, 0, ALPHA,
                 Mptr( A, Aii, Aoffj+1, Ald, size ), Ald, X+XYii*Xinc, LDX,
                 Y+(XYoffj+1)*Yinc, LDY );
            Aii  += tmp1;
            XYii += tmp1;
            m1   -= tmp1;
         }
/*
*  Return if no more row in the LCM table.
*/
         if( mblks <= 0 ) return;
/*
*  lcmt00 <= upp. The current block owns either diagonals or lower entries.
*  Save the current position in the LCM table. After this column has been
*  completely taken care of, re-start from this row and the next column of
*  the LCM table.
*/
         lcmt = lcmt00; mblkd = mblks; ioffd = Aoffi; ioffxy = XYoffi;

         mbloc = Amb;
         while( ( mblkd > 0 ) && ( lcmt >= low ) )
         {
/*
*  A block owning diagonals lcmt00 >= low && lcmt00 <= upp has been found.
*/
            if( mblkd == 1 ) mbloc = lmbloc;
            TRM( TYPE, SIDE, UPLO, TRANS, DIAG, mbloc, nbloc, K, lcmt,
                 ALPHA, Mptr( A, ioffd+1, Aoffj+1, Ald, size ), Ald,
                 X+(ioffxy+1)*Xinc, LDX, Y+(XYoffj+1)*Yinc, LDY );
            lcmt00 = lcmt;  lcmt   -= pmb;
            mblks  = mblkd; mblkd--;
            Aoffi  = ioffd; XYoffi = ioffxy;
            ioffd += mbloc; ioffxy += mbloc;
         }
/*
*  Operate with the lower triangular part of sub( A ).
*/
         tmp1 = m1 - ioffd + Aii - 1;
         if( lower && ( tmp1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, nbloc, K, 0, ALPHA,
                 Mptr( A, ioffd+1, Aoffj+1, Ald, size ), Ald,
                 X+(ioffxy+1)*Xinc, LDX, Y+(XYoffj+1)*Yinc, LDY );
         }

         tmp1    = MIN( Aoffi, iimax ) - Aii + 1;
         m1     -= tmp1;
         n1     -= nbloc;
         lcmt00 += qnb;
         nblks--;
         Aoffj  += nbloc;
         XYoffj += nbloc;
/*
*  Operate with the upper triangular part of sub( A ).
*/
         if( upper && ( tmp1 > 0 ) && ( n1 > 0 ) )
         {
            TRM( TYPE, SIDE, ALL, TRANS, DIAG, tmp1, n1, K, 0, ALPHA,
                 Mptr( A, Aii, Aoffj+1, Ald, size ), Ald, X+XYii*Xinc, LDX,
                 Y+(XYoffj+1)*Yinc, LDY );
         }
         Aii  = Aoffi  + 1;  Ajj = Aoffj  + 1;
         XYii = XYoffi + 1; XYjj = XYoffj + 1;
      }
   }
/*
*  End of PB_Cptrm
*/
}
