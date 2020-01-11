/* ---------------------------------------------------------------------
*
*  -- PBLAS routine (version 1.5) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 17, 1995
*
*  ---------------------------------------------------------------------
*/
/*
*  Include files
*/
#include "tools.h"

void pbchkvect( n, npos0, ix, jx, desc_X, incx, dpos0, iix, jjx, ixrow,
                ixcol, nprow, npcol, myrow, mycol, info )
/*
*  .. Scalar Arguments ..
*/
   Int         dpos0, * iix, incx, * info, ix, * ixcol, * ixrow, * jjx,
               jx, myrow, mycol, npcol, nprow, n, npos0;
/*
*  .. Array Arguments ..
*/
   Int         desc_X[];
{
/*
*
*  Purpose
*  =======
*
*  pbchkvect checks the validity of a descriptor vector DESCX, the
*  related global indexes IX, JX and the global increment INCX. It also
*  computes the starting local indexes (IIX,JJX) corresponding to the
*  submatrix starting globally at the entry pointed by (IX,JX).
*  Moreover, this routine returns the coordinates in the grid of the
*  process owning the global matrix entry of indexes (IX,JX), namely
*  (IXROW,IXCOL). The routine prevents out-of-bound memory access
*  by performing the appropriate MIN operation on iix and JJX.  Finally,
*  if an inconsistency is found among its parameters IX, JX, DESCX and
*  INCX, the routine returns an error code in info.
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The length of the vector X being operated on.
*
*  NPOS0   (global input) INTEGER
*          Where in the calling routine's parameter list N appears.
*
*  IX      (global input) INTEGER
*          X's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JX      (global input) INTEGER
*          X's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*
*  INCX    (global input) INTEGER
*          The global increment for the elements of X. Only two values
*          of INCX are supported in this version, namely 1 and M_X.
*          INCX must not be zero.
*
*  DPOS0   (global input) INTEGER
*          Where in the calling routine's parameter list DESCX
*          appears.  Note that we assume IX and JX are respectively 2
*          and 1 entries behind DESCX, and INCX is 1 entry after DESCX.
*
*  IIX     (local output) pointer to INTEGER
*          The local rows starting index of the submatrix.
*
*  JJX     (local output) pointer to INTEGER
*          The local columns starting index of the submatrix.
*
*  IXROW   (global output) pointer to INTEGER
*          The row coordinate of the process that possesses the first
*          row and column of the submatrix.
*
*  IXCOL   (global output) pointer to INTEGER
*          The column coordinate of the process that possesses the
*          first row and column of the submatrix.
*
*  NPROW   (global input) INTEGER
*          The total number of process rows over which the distributed
*          matrix is distributed.
*
*  NPCOL   (global input) INTEGER
*          The total number of process columns over which the
*          distributed matrix is distributed.
*
*  MYROW   (local input) INTEGER
*          The row coordinate of the process calling this routine.
*
*  MYCOL   (local input) INTEGER
*          The column coordinate of the process calling this routine.
*
*  INFO    (local input/local output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*
*  =====================================================================
*
*  .. Parameters ..
*/
#define DESCMULT      100
#define BIGNUM      10000
/* ..
*  .. Local Scalars ..
*/
   Int         descpos, ExtraColBlock, ExtraRowBlock, icpos, ixpos,
               jxpos, MyColBlock, MyColDist, MyRowBlock, MyRowDist,
               NColBlock, np, npos, nq, NRowBlock;
/* ..
*  .. External Functions ..
*/
   F_INTG_FCT  numroc_();
/*
*  .. Executable Statements ..
*/
   if( *info >= 0 )
      *info = BIGNUM;
   else if( *info < -DESCMULT )
      *info = -(*info);
   else
      *info = -(*info) * DESCMULT;
/*
*  Figure where in parameter list each parameter was, factoring in
*  descriptor multiplier
*/
   npos = npos0 * DESCMULT;
   ixpos = ( dpos0 - 2 ) * DESCMULT;
   jxpos = ( dpos0 - 1 ) * DESCMULT;
   icpos = ( dpos0 + 1 ) * DESCMULT;
   descpos = dpos0 * DESCMULT + 1;
/*
 * Check that we have a legal descriptor type
 */
   if(desc_X[DT_] != BLOCK_CYCLIC_2D) *info = MIN( *info, descpos + DT_ );
/*
*  Check that matrix values make sense from local viewpoint
*/
   if( n < 0 )
      *info = MIN( *info, npos );
   else if( ix < 1 )
      *info = MIN( *info, ixpos );
   else if( jx < 1 )
      *info = MIN( *info, jxpos );
   else if( desc_X[MB_] < 1 )
      *info = MIN( *info, descpos + MB_ );
   else if( desc_X[NB_] < 1 )
      *info = MIN( *info, descpos + NB_ );
   else if( ( desc_X[RSRC_] < 0 ) || ( desc_X[RSRC_] >= nprow ) )
      *info = MIN( *info, descpos + RSRC_ );
   else if( ( desc_X[CSRC_] < 0 ) || ( desc_X[CSRC_] >= npcol ) )
      *info = MIN( *info, descpos + CSRC_ );
   else if( incx != 1 && incx != desc_X[M_] )
      *info = MIN( *info, icpos );
   else if( desc_X[LLD_] < 1 )
      *info = MIN( *info, descpos + LLD_ );

   if( n == 0 )
   {
/*
*     NULL matrix, relax some checks
*/
      if( desc_X[M_] < 0 )
         *info = MIN( *info, descpos + M_ );
      if( desc_X[N_] < 0 )
         *info = MIN( *info, descpos + N_ );
   }
   else
   {
/*
*     more rigorous checks for non-degenerate matrices
*/
      if( desc_X[M_] < 1 )
         *info = MIN( *info, descpos + M_ );
      else if( desc_X[N_] < 1 )
         *info = MIN( *info, descpos + N_ );
      else if( ( incx == desc_X[M_] ) && ( jx+n-1 > desc_X[N_] ) )
         *info = MIN( *info, jxpos );
      else if( ( incx == 1 ) && ( incx != desc_X[M_] ) &&
               ( ix+n-1 > desc_X[M_] ) )
         *info = MIN( *info, ixpos );
      else
      {
         if( ix > desc_X[M_] )
            *info = MIN( *info, ixpos );
         else if( jx > desc_X[N_] )
            *info = MIN( *info, jxpos );
      }
   }
/*
*  Retrieve local information for vector X, and prepare output:
*  set info = 0 if no error, and divide by DESCMULT if error is not
*  in a descriptor entry.
*/
   if( *info == BIGNUM )
   {
      MyRowDist = ( myrow + nprow - desc_X[RSRC_] ) % nprow;
      MyColDist = ( mycol + npcol - desc_X[CSRC_] ) % npcol;
      NRowBlock = desc_X[M_] / desc_X[MB_];
      NColBlock = desc_X[N_] / desc_X[NB_];
      np = ( NRowBlock / nprow ) * desc_X[MB_];
      nq = ( NColBlock / npcol ) * desc_X[NB_];
      ExtraRowBlock = NRowBlock % nprow;
      ExtraColBlock = NColBlock % npcol;

      ix--;
      jx--;
      MyRowBlock = ix / desc_X[MB_];
      MyColBlock = jx / desc_X[NB_];
      *ixrow = ( MyRowBlock + desc_X[RSRC_] ) % nprow;
      *ixcol = ( MyColBlock + desc_X[CSRC_] ) % npcol;

      *iix = ( MyRowBlock / nprow + 1 ) * desc_X[MB_] + 1;
      *jjx = ( MyColBlock / npcol + 1 ) * desc_X[NB_] + 1;

      if( MyRowDist >= ( MyRowBlock % nprow ) )
      {
         if( myrow == *ixrow )
            *iix += ix % desc_X[MB_];
         *iix -= desc_X[MB_];
      }
      if( MyRowDist  < ExtraRowBlock )
         np += desc_X[MB_];
      else if( MyRowDist == ExtraRowBlock )
         np += ( desc_X[M_] % desc_X[MB_] );
      np = MAX( 1, np );

      if( MyColDist >= ( MyColBlock % npcol ) )
      {
         if( mycol == *ixcol )
            *jjx += jx % desc_X[NB_];
         *jjx -= desc_X[NB_];
      }
      if( MyColDist < ExtraColBlock )
         nq += desc_X[NB_];
      else if( MyColDist == ExtraColBlock )
         nq += ( desc_X[N_] % desc_X[NB_] );
      nq = MAX( 1, nq );

      *iix = MIN( *iix, np );
      *jjx = MIN( *jjx, nq );

      if( desc_X[LLD_] < np )
      {
         if( numroc_(&desc_X[N_], &desc_X[NB_], &mycol, &desc_X[CSRC_], &npcol) )
            *info = -( descpos + LLD_ );
         else *info = 0;
      }
      else *info = 0;
   }
   else if( *info % DESCMULT == 0 )
   {
      *info = -(*info) / DESCMULT;
   }
   else
   {
      *info = -(*info);
   }
}
