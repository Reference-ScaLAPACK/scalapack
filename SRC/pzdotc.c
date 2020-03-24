/* ---------------------------------------------------------------------
*
*  Mark R. Fahey
*  August 2000
*  This is a slightly modified version of pzaxpy_ from ScaLAPACK 1.0
*  which fixes a bug in the incx=1 and incy=1 case.
*
*  ---------------------------------------------------------------------
*/
/*
*  Include files
*/
#include "pblas.h"

void pzdotc_( n, dotc, X, ix, jx, desc_X, incx, Y, iy, jy, desc_Y,
              incy )
/*
*  .. Scalar Arguments ..
*/
   Int         * incx, * incy, * ix, * iy, * jx, * jy, * n;
   complex16   * dotc;
/* ..
*  .. Array Arguments ..
*/
   Int         desc_X[], desc_Y[];
   complex16   X[], Y[];
{
/*
*  Purpose
*  =======
*
*  PZDOTC forms the dot product of two distributed vectors,
*
*     dotc := sub( X )**H * sub( Y )
*
*  where sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
*                         X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X,
*
*        sub( Y ) denotes Y(IY,JY:JY+N-1) if INCY = M_Y,
*                         Y(IY:IY+N-1,JY) if INCY = 1 and INCY <> M_Y.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector descA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DT_A   (global) descA[ DT_ ]   The descriptor type.  In this case,
*                                 DT_A = 1.
*  CTXT_A (global) descA[ CTXT_ ] The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) descA[ M_ ]    The number of rows in the global
*                                 array A.
*  N_A    (global) descA[ N_ ]    The number of columns in the global
*                                 array A.
*  MB_A   (global) descA[ MB_ ]   The blocking factor used to distribu-
*                                 te the rows of the array.
*  NB_A   (global) descA[ NB_ ]   The blocking factor used to distribu-
*                                 te the columns of the array.
*  RSRC_A (global) descA[ RSRC_ ] The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) descA[ CSRC_ ] The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  descA[ LLD_ ]  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Because vectors may be seen as particular matrices, a distributed
*  vector is considered to be a distributed matrix.
*
*  If INCX = M_X and INCY = M_Y, NB_X must be equal to NB_Y, and the
*  process column having the first entries of sub( Y ) must also contain
*  the first entries of sub( X ). Moreover, the quantity
*  MOD( JX-1, NB_X ) must be equal to MOD( JY-1, NB_Y ).
*
*  If INCX = M_X, INCY = 1 and INCY <> M_Y, NB_X must be equal to MB_Y.
*  Moreover, the quantity MOD( JX-1, NB_X ) must be equal to
*  MOD( IY-1, MB_Y ).
*
*  If INCX = 1, INCX <> M_X and INCY = M_Y, MB_X must be equal to NB_Y.
*  Moreover, the quantity MOD( IX-1, MB_X ) must be equal to
*  MOD( JY-1, NB_Y ).
*
*  If INCX = 1, INCX <> M_X, INCY = 1 and INCY <> M_Y, MB_X must be
*  equal to MB_Y, and the process row having the first entries of
*  sub( Y ) must also contain the first entries of sub( X ). Moreover,
*  the quantity MOD( IX-1, MB_X ) must be equal to MOD( IY-1, MB_Y ).
*
*
*  Parameters
*  ==========
*
*  N       (global input) pointer to INTEGER
*          The length of the distributed vectors to be multiplied.
*          N >= 0.
*
*  DOTC    (local output) pointer to COMPLEX*16
*          The dot product of sub( X ) and sub( Y ) only in their scope.
*
*  X       (local input) COMPLEX*16 array containing the local
*          pieces of a distributed matrix of dimension of at least
*              ( (JX-1)*M_X + IX + ( N - 1 )*abs( INCX ) )
*          This array contains the entries of the distributed vector
*          sub( X ).
*
*  IX      (global input) pointer to INTEGER
*          The global row index of the submatrix of the distributed
*          matrix X to operate on.
*
*  JX      (global input) pointer to INTEGER
*          The global column index of the submatrix of the distributed
*          matrix X to operate on.
*
*  DESCX   (global and local input) INTEGER array of dimension 8.
*          The array descriptor of the distributed matrix X.
*
*  INCX    (global input) pointer to INTEGER
*          The global increment for the elements of X. Only two values
*          of INCX are supported in this version, namely 1 and M_X.
*
*  Y       (local input) COMPLEX*16 array containing the local
*          pieces of a distributed matrix of dimension of at least
*              ( (JY-1)*M_Y + IY + ( N - 1 )*abs( INCY ) )
*          This array contains the entries of the distributed vector
*          sub( Y ).
*
*  IY      (global input) pointer to INTEGER
*          The global row index of the submatrix of the distributed
*          matrix Y to operate on.
*
*  JY      (global input) pointer to INTEGER
*          The global column index of the submatrix of the distributed
*          matrix Y to operate on.
*
*  DESCY   (global and local input) INTEGER array of dimension 8.
*          The array descriptor of the distributed matrix Y.
*
*  INCY    (global input) pointer to INTEGER
*          The global increment for the elements of Y. Only two values
*          of INCY are supported in this version, namely 1 and M_Y.
*
*  =====================================================================
*
*  .. Local Scalars ..
*/
   char        * cbtop, * cctop, * rbtop, * rctop;
   Int         ictxt, iix, iiy, info, ixcol, ixrow, iycol, iyrow, jjx,
               jjy, lcm, lcmp, mone=-1, mycol, myrow, nn, np, np0,
               nprow, npcol, nq, nz, ione=1, tmp1, wksz;
   complex16   xwork[1], ywork[1], zero;
/* ..
*  .. PBLAS Buffer ..
*/
   complex16   * buff;
/* ..
*  .. External Functions ..
*/
   void        blacs_gridinfo_();
   void        zgebr2d_();
   void        zgebs2d_();
   void        zgerv2d_();
   void        zgesd2d_();
   void        zgsum2d_();
   void        pbchkvect();
   void        pberror_();
   char        * getpbbuf();
   char        * ptop();
   F_VOID_FCT  pbztrnv_();
   F_VOID_FCT  zzdotc_();
   F_INTG_FCT  ilcm_();
/* ..
*  .. Executable Statements ..
*
*  Get grid parameters
*/
   ictxt = desc_X[CTXT_];
   blacs_gridinfo_( &ictxt, &nprow, &npcol, &myrow, &mycol );
/*
*  Test the input parameters
*/
   info = 0;
   if( nprow == -1 )
      info = -(600+CTXT_+1);
   else
   {
      pbchkvect( *n, 1, *ix, *jx, desc_X, *incx, 6, &iix, &jjx,
                 &ixrow, &ixcol, nprow, npcol, myrow, mycol, &info );
      pbchkvect( *n, 1, *iy, *jy, desc_Y, *incy, 11, &iiy, &jjy,
                 &iyrow, &iycol, nprow, npcol, myrow, mycol, &info );

      if( info == 0 )
      {
         if( *n != 1 )
         {
            if( *incx == desc_X[M_] )
            {                 /* X is distributed along a process row */
               if( *incy == desc_Y[M_] )
               {               /* Y is distributed over a process row */
                  if( ( ixcol != iycol ) ||
                      ( ( (*jx-1) % desc_X[NB_] ) !=
                        ( (*jy-1) % desc_Y[NB_] ) ) )
                     info = -10;
                  else if( desc_Y[NB_] != desc_X[NB_] )
                     info = -(1100+NB_+1);
               }
               else if( ( *incy == 1 ) && ( *incy != desc_Y[M_] ) )
               {            /* Y is distributed over a process column */
                  if( ( (*jx-1) % desc_X[NB_] ) != ( (*iy-1) % desc_Y[MB_] ) )
                     info = -9;
                  else if( desc_Y[MB_] != desc_X[NB_] )
                     info = -(1100+MB_+1);
               }
               else
               {
                  info = -12;
               }
            }
            else if( ( *incx == 1 ) && ( *incx != desc_X[M_] ) )
            {              /* X is distributed along a process column */
               if( *incy == desc_Y[M_] )
               {                  /* Y is distributed over a process row */
                  if( ( (*ix-1) % desc_X[MB_] ) != ( (*jy-1) % desc_Y[NB_] ) )
                     info = -10;
                  else if( desc_Y[NB_] != desc_X[MB_] )
                     info = -(1100+NB_+1);
               }
               else if( ( *incy == 1 ) && ( *incy != desc_Y[M_] ) )
               {            /* Y is distributed over a process column */
                  if( ( ixrow != iyrow ) ||
                      ( ( (*ix-1) % desc_X[MB_] ) !=
                        ( (*iy-1) % desc_Y[MB_] ) ) )
                     info = -9;
                  else if( desc_Y[MB_] != desc_X[MB_] )
                     info = -(1100+MB_+1);
               }
               else
               {
                  info = -12;
               }
            }
            else
            {
               info = -7;
            }
         }
         if( ictxt != desc_Y[CTXT_] )
            info = -(1100+CTXT_+1);
      }
   }
   if( info )
   {
      pberror_( &ictxt, "PZDOTC", &info );
      return;
   }
/*
*  Quick return if possible.
*/
   dotc->re = ZERO;
   dotc->im = ZERO;
   zero.re  = ZERO;
   zero.im  = ZERO;
   if( *n == 0 ) return;
/*
*  dot <- x^{h} * y
*/
   if( *n == 1 )
   {
      if( ( myrow == ixrow ) && ( mycol == ixcol ) )
      {
         buff = &X[iix-1+(jjx-1)*desc_X[LLD_]];
         if( ( myrow != iyrow ) || ( mycol != iycol ) )
         {
            zgesd2d_( &ictxt, n, n, buff, n, &iyrow, &iycol );
            zgerv2d_( &ictxt, n, n, ywork, n, &iyrow, &iycol );
         }
         else
            *ywork = Y[iiy-1+(jjy-1)*desc_Y[LLD_]];
         zzdotc_( n, dotc, buff, n, ywork, n );
      }
      else if( ( myrow == iyrow ) && ( mycol == iycol ) )
      {
         zgesd2d_( &ictxt, n, n, &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], n,
                   &ixrow, &ixcol );
         zgerv2d_( &ictxt, n, n, xwork, n, &ixrow, &ixcol );
         zzdotc_( n, dotc, xwork, n,
                  &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], n );
      }

      if( ( *incx == desc_X[M_] ) && ( desc_X[M_] != 1 ) )
      {
         if( myrow == ixrow )
         {
            rbtop = ptop( BROADCAST, ROW, TOPGET );
            if( mycol == ixcol )
            {
               zgebs2d_( &ictxt, C2F_CHAR( ROW ), C2F_CHAR( rbtop ),
                         &ione, &ione, dotc, &ione );
            }
            else
            {
               zgebr2d_( &ictxt, C2F_CHAR( ROW ), C2F_CHAR( rbtop ),
                         &ione, &ione, dotc, &ione, &myrow, &ixcol );
            }
         }
      }
      else if( ( *incx == 1 ) && ( desc_X[M_] != 1 ) )
      {
         if( mycol == ixcol )
         {
            cbtop = ptop( BROADCAST, COLUMN, TOPGET );
            if( myrow == ixrow )
            {
               zgebs2d_( &ictxt, C2F_CHAR( COLUMN ), C2F_CHAR( cbtop ),
                         &ione, &ione, dotc, &ione );
            }
            else
            {
               zgebr2d_( &ictxt, C2F_CHAR( COLUMN ), C2F_CHAR( cbtop ),
                         &ione, &ione, dotc, &ione, &ixrow, &mycol );
            }
         }
      }

      if( ( *incy == desc_Y[M_] ) && ( desc_Y[M_] != 1 ) )
      {
         if( myrow == iyrow )
         {
            rbtop = ptop( BROADCAST, ROW, TOPGET );
            if( mycol == iycol )
            {
               zgebs2d_( &ictxt, C2F_CHAR( ROW ), C2F_CHAR( rbtop ),
                         &ione, &ione, dotc, &ione );
            }
            else
            {
               zgebr2d_( &ictxt, C2F_CHAR( ROW ), C2F_CHAR( rbtop ),
                         &ione, &ione, dotc, &ione, &myrow, &iycol );
            }
         }
      }
      else if( ( *incy == 1 ) && ( desc_Y[M_] != 1 ) )
      {
         if( mycol == iycol )
         {
            cbtop = ptop( BROADCAST, COLUMN, TOPGET );
            if( myrow == iyrow )
            {
               zgebs2d_( &ictxt, C2F_CHAR( COLUMN ), C2F_CHAR( cbtop ),
                         &ione, &ione, dotc, &ione );
            }
            else
            {
               zgebr2d_( &ictxt, C2F_CHAR( COLUMN ), C2F_CHAR( cbtop ),
                         &ione, &ione, dotc, &ione, &iyrow, &mycol );
            }
         }
      }
      return;
   }

   if( ( *incx == desc_X[M_] ) && ( *incy == desc_Y[M_] ) )
   {               /* X and Y are both distributed over a process row */
      nz = (*jx-1) % desc_Y[NB_];
      nn = *n + nz;
      nq = numroc_( &nn, &desc_X[NB_], &mycol, &ixcol, &npcol );
      if( mycol == ixcol )
         nq -= nz;
      if( ixrow == iyrow )
      {
         if( myrow == ixrow )
         {
            rctop = ptop( COMBINE, ROW, TOPGET );
            zzdotc_( &nq, dotc,
                     &X[iix-1+(jjx-1)*desc_X[LLD_]], &desc_X[LLD_],
                     &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], &desc_Y[LLD_] );
            zgsum2d_( &ictxt, C2F_CHAR( ROW ), C2F_CHAR( rctop ), &ione,
                      &ione, dotc, &ione, &mone, &mycol );
         }
      }
      else
      {
         if( myrow == ixrow )
         {
            rctop = ptop( COMBINE, ROW, TOPGET );
            zgesd2d_( &ictxt, &ione, &nq,
                      &X[iix-1+(jjx-1)*desc_X[LLD_]], &desc_X[LLD_],
                      &iyrow, &mycol );
            buff = (complex16 *)getpbbuf( "PZDOTC", nq*sizeof(complex16) );
            zgerv2d_( &ictxt, &nq, &ione, buff, &ione,
                      &ixrow, &mycol );
            zzdotc_( &nq, dotc, &X[iix-1+(jjx-1)*desc_X[LLD_]],
                          &desc_X[LLD_], buff, &ione );
            zgsum2d_( &ictxt, C2F_CHAR( ROW ), C2F_CHAR( rctop ), &ione,
                      &ione, dotc, &ione, &mone, &mycol );
         }
         else if( myrow == iyrow )
         {
            rctop = ptop( COMBINE, ROW, TOPGET );
            zgesd2d_( &ictxt, &ione, &nq,
                      &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], &desc_Y[LLD_],
                      &ixrow, &mycol );
            buff = (complex16 *)getpbbuf( "PZDOTC", nq*sizeof(complex16) );
            zgerv2d_( &ictxt, &nq, &ione, buff, &ione, &ixrow,
                      &mycol );
            zzdotc_( &nq, dotc,
                     buff, &ione,
                     &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], &desc_Y[LLD_] );
            zgsum2d_( &ictxt, C2F_CHAR( ROW ), C2F_CHAR( rctop ), &ione,
                      &ione, dotc, &ione, &mone, &mycol );
         }
      }
   }
   else if( ( *incx == 1 ) && ( *incx != desc_X[M_] ) &&
            ( *incy == 1 ) && ( *incy != desc_Y[M_] ) )
   {            /* X and Y are both distributed over a process column */
      nz = (*ix-1) % desc_X[MB_];
      nn = *n + nz;
      np = numroc_( &nn, &desc_X[MB_], &myrow, &ixrow, &nprow );
      if( myrow == ixrow )
         np -= nz;
      if( ixcol == iycol )
      {
         if( mycol == ixcol )
         {
            cctop = ptop( COMBINE, COLUMN, TOPGET );
            zzdotc_( &np, dotc,
                     &X[iix-1+(jjx-1)*desc_X[LLD_]], incx,
                     &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], incy );
            zgsum2d_( &ictxt, C2F_CHAR( COLUMN ), C2F_CHAR( cctop ),
                      &ione, &ione, dotc, &ione, &mone, &mycol );
         }
      }
      else
      {
         if( mycol == ixcol )
         {
            cctop = ptop( COMBINE, COLUMN, TOPGET );
            zgesd2d_( &ictxt, &np, &ione,
                      &X[iix-1+(jjx-1)*desc_X[LLD_]], &desc_X[LLD_],
                      &myrow, &iycol );
            buff = (complex16 *)getpbbuf( "PZDOTC", np*sizeof(complex16) );
            zgerv2d_( &ictxt, &np, &ione, buff, &ione,
                      &myrow, &iycol );
            zzdotc_( &np, dotc,
                     &X[iix-1+(jjx-1)*desc_X[LLD_]], incx,
                     buff, &ione );
            zgsum2d_( &ictxt, C2F_CHAR( COLUMN ), C2F_CHAR( cctop ),
                      &ione, &ione, dotc, &ione, &mone, &mycol );
         }
         else if( mycol == iycol )
         {
            cctop = ptop( COMBINE, COLUMN, TOPGET );
            buff = (complex16 *)getpbbuf( "PZDOTC", np*sizeof(complex16) );
            zgerv2d_( &ictxt, &np, &ione, buff, &ione,
                      &myrow, &ixcol );
            zgesd2d_( &ictxt, &np, &ione,
                      &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], &desc_Y[LLD_],
                      &myrow, &ixcol );
            zzdotc_( &np, dotc,
                     buff, &ione,
                     &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], incy );
            zgsum2d_( &ictxt, C2F_CHAR( COLUMN ), C2F_CHAR( cctop ),
                      &ione, &ione, dotc, &ione, &mone, &mycol );
         }
      }
   }
   else       /* X and Y are not distributed along the same direction */
   {
      lcm = ilcm_( &nprow, &npcol );
      if( ( *incx == 1 ) && ( *incx != desc_X[M_] ) )
      {                     /* X is distributed over a process column */
         lcmp = lcm / nprow;
         nz = (*jy-1) % desc_Y[NB_];
         nn = *n + nz;
         tmp1 = nn / desc_Y[MB_];
         np = numroc_( &nn, &desc_X[MB_], &myrow, &ixrow, &nprow );
         np0 = MYROC0( tmp1, nn, desc_X[MB_], nprow );
         tmp1 = np0 / desc_X[MB_];
         wksz = MYROC0( tmp1, np0, desc_X[MB_], lcmp );
         wksz = np + wksz;

         buff = (complex16 *)getpbbuf( "PZDOTC", wksz*sizeof(complex16) );

         if( mycol == iycol )
            jjy -= nz;
         if( myrow == ixrow )
            np -= nz;
         pbztrnv_( &ictxt, C2F_CHAR( "R" ), C2F_CHAR( "T" ), n,
                   &desc_Y[NB_], &nz, &Y[iiy-1+(jjy-1)*desc_Y[LLD_]],
                   &desc_Y[LLD_], &zero, buff, &ione, &iyrow, &iycol,
                   &ixrow, &ixcol, buff+np );
         if( mycol == ixcol )
         {
            cctop = ptop( COMBINE, COLUMN, TOPGET );
            zzdotc_( &np, dotc, &X[iix-1+(jjx-1)*desc_X[LLD_]],
                     incx, buff, &ione );
            zgsum2d_( &ictxt, C2F_CHAR( COLUMN ), C2F_CHAR( cctop ),
                      &ione, &ione, dotc, &ione, &mone, &mycol );
         }
         if( myrow == iyrow )
         {
            rbtop = ptop( BROADCAST, ROW, TOPGET );
            if( mycol == ixcol )
               zgebs2d_( &ictxt, C2F_CHAR( ROW ), C2F_CHAR( rbtop ),
                        &ione, &ione, dotc, &ione );
            else
               zgebr2d_( &ictxt, C2F_CHAR( ROW ), C2F_CHAR( rbtop ),
                        &ione, &ione, dotc, &ione, &myrow, &ixcol );
         }
      }
      else                  /* Y is distributed over a process column */
      {
         lcmp = lcm / nprow;
         nz = (*jx-1) % desc_X[NB_];
         nn = *n + nz;
         tmp1 = nn / desc_X[MB_];
         np = numroc_( &nn, desc_Y+MB_, &myrow, &iyrow, &nprow );
         np0 = MYROC0( tmp1, nn, desc_Y[MB_], nprow );
         tmp1 = np0 / desc_Y[MB_];
         wksz = MYROC0( tmp1, np0, desc_Y[MB_], lcmp );
         wksz = np + wksz;

         buff = (complex16 *)getpbbuf( "PZDOTC", wksz*sizeof(complex16) );

         if( myrow == iyrow )
            np -= nz;
         pbztrnv_( &ictxt, C2F_CHAR( "R" ), C2F_CHAR( "T" ), n,
                   &desc_X[NB_], &nz, &X[iix-1+(jjx-1)*desc_X[LLD_]],
                   &desc_X[LLD_], &zero, buff, &ione, &ixrow, &ixcol,
                   &iyrow, &iycol, buff+np );
         if( mycol == iycol )
         {
            cctop = ptop( COMBINE, COLUMN, TOPGET );
            zzdotc_( &np, dotc, buff, &ione,
                     &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], incy );
            zgsum2d_( &ictxt, C2F_CHAR( COLUMN ), C2F_CHAR( cctop ),
                      &ione, &ione, dotc, &ione, &mone, &mycol );
         }
         if( myrow == ixrow )
         {
            rbtop = ptop( BROADCAST, ROW, TOPGET );
            if( mycol == iycol )
               zgebs2d_( &ictxt, C2F_CHAR( ROW ), C2F_CHAR( rbtop ),
                         &ione, &ione, dotc, &ione );
            else
               zgebr2d_( &ictxt, C2F_CHAR( ROW ), C2F_CHAR( rbtop ),
                         &ione, &ione, dotc, &ione, &myrow, &iycol );
         }
      }
   }
}
