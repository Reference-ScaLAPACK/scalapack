/* ---------------------------------------------------------------------
*
*  -- Mark R. Fahey
*     June 28, 2000
*
*  ---------------------------------------------------------------------
*/
/*
*  Include files
*/
#include "pblas.h"

void pzrot_( Int *n, complex16 X[], Int *ix, Int *jx, Int desc_X[], Int *incx, complex16 Y[], Int *iy, Int *jy, Int desc_Y[], Int *incy, double *c, complex16 *s )
/*
*  Mark Fahey
*  June 22, 2000
*/
{
/*
*  Purpose
*  =======
*
*  PZROT applies a plane rotation, where the cos (C) is real and the
*  sin (S) is complex, and the vectors CX and CY are complex, i.e.,
*
*     [ sub( X ) ] := [  C          S  ] [ sub( X ) ]
*     [ sub( Y ) ] := [ -conjg(S)   C  ] [ sub( Y ) ]
*
*  where sub( X ) denotes X(IX,JX:JX+N-1) if INCX = M_X,
*                         X(IX:IX+N-1,JX) if INCX = 1 and INCX <> M_X,
*
*        sub( Y ) denotes Y(IY,JY:JY+N-1) if INCY = M_Y,
*                         Y(IY:IY+N-1,JY) if INCY = 1 and INCY <> M_Y,
*
*  and where C*C + S*CONJG(S) = 1.0.
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
*  Such a global array has an associated description vector DESCA.
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
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements in the vectors CX and CY.
*
*  X       (local input) COMPLEX array containing the local
*          pieces of a distributed matrix of dimension of at least
*              ( (JX-1)*M_X + IX + ( N - 1 )*abs( INCX ) )
*          This array contains the entries of the distributed vector
*          sub( X ).
*          On output, CX is overwritten with C*X + S*Y.
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
*  Y       (local input) COMPLEX array containing the local
*          pieces of a distributed matrix of dimension of at least
*              ( (JY-1)*M_Y + IY + ( N - 1 )*abs( INCY ) )
*          This array contains the entries of the distributed vector
*          sub( Y ).
*          On output, CY is overwritten with -CONJG(S)*X + C*Y.
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
*  C       (input) pointer to DOUBLE 
*  S       (input) pointer COMPLEX
*          C and S define a rotation
*             [  C          S  ]
*             [ -conjg(S)   C  ]
*          where C*C + S*CONJG(S) = 1.0.
*
* =====================================================================
*
*  .. Local Scalars ..
*/
   Int         ictxt, iix, iiy, info, ixcol, ixrow, iycol, iyrow, jjx,
               jjy, lcm, lcmp, mycol, myrow, nn, np, np0,
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
   void        zgerv2d_();
   void        zgesd2d_();
   void        pbchkvect();
   void        PB_Cabort();
   char        * getpbbuf();
   F_INTG_FCT  pbztrnv_();
   F_INTG_FCT  zrot_();
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
      info = -(500+CTXT_+1);
   else
   {
      pbchkvect( *n, 1, *ix, *jx, desc_X, *incx, 5, &iix, &jjx,
                 &ixrow, &ixcol, nprow, npcol, myrow, mycol, &info );
      pbchkvect( *n, 1, *iy, *jy, desc_Y, *incy, 10, &iiy, &jjy,
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
                     info = -9;
                  else if( desc_Y[NB_] != desc_X[NB_] )
                     info = -(1000+NB_+1);
               }
               else if( ( *incy == 1 ) && ( *incy != desc_Y[M_] ) )
               {            /* Y is distributed over a process column */
                  if( ( (*jx-1) % desc_X[NB_] ) != ( (*iy-1) % desc_Y[MB_] ) )
                     info = -8;
                  else if( desc_Y[MB_] != desc_X[NB_] )
                     info = -(1000+MB_+1);
               }
               else
               {
                  info = -11;
               }
            }
            else if( ( *incx == 1 ) && ( *incx != desc_X[M_] ) )
            {              /* X is distributed along a process column */
               if( *incy == desc_Y[M_] )
               {                  /* Y is distributed over a process row */
                  if( ( (*ix-1) % desc_X[MB_] ) != ( (*jy-1) % desc_Y[NB_] ) )
                     info = -9;
                  else if( desc_Y[NB_] != desc_X[MB_] )
                     info = -(1000+NB_+1);
               }
               else if( ( *incy == 1 ) && ( *incy != desc_Y[M_] ) )
               {            /* Y is distributed over a process column */
                  if( ( ixrow != iyrow ) ||
                      ( ( (*ix-1) % desc_X[MB_] ) !=
                        ( (*iy-1) % desc_Y[MB_] ) ) )
                     info = -8;
                  else if( desc_Y[MB_] != desc_X[MB_] )
                     info = -(1000+MB_+1);
               }
               else
               {
                  info = -11;
               }
            }
            else
            {
               info = -6;
            }
         }
         if( ictxt != desc_Y[CTXT_] )
            info = -(1000+CTXT_+1);
      }
   }
   if( info ) { PB_Cabort( ictxt, "PZROT", info ); return; }
/*
   if( info )
   {
      pberror_( &ictxt, "PZROT", &info );
      return;
   }
*/

/*
*  Quick return if possible.
*/
   zero.re  = ZERO;
   zero.im  = ZERO;
   if( *n == 0 ) return;
/*
*  rotation
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
         zrot_( n, buff, n, ywork, n, c, s );
         X[iix-1+(jjx-1)*desc_X[LLD_]] = *buff;
         if( ( myrow == iyrow ) && ( mycol == iycol ) )
            Y[iiy-1+(jjy-1)*desc_Y[LLD_]] = *ywork;
      }
      else if( ( myrow == iyrow ) && ( mycol == iycol ) )
      {
         zgesd2d_( &ictxt, n, n, &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], n,
                   &ixrow, &ixcol );
         zgerv2d_( &ictxt, n, n, xwork, n, &ixrow, &ixcol );
         zrot_( n, xwork, n, &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], n, c, s );
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
            zrot_( &nq, &X[iix-1+(jjx-1)*desc_X[LLD_]], &desc_X[LLD_],
                        &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], &desc_Y[LLD_], c, s );
         }
      }
      else
      {
         if( myrow == ixrow )
         {
            zgesd2d_( &ictxt, &ione, &nq,
                      &X[iix-1+(jjx-1)*desc_X[LLD_]], &desc_X[LLD_],
                      &iyrow, &mycol );
            buff = (complex16 *)getpbbuf( "PZROT", nq*sizeof(complex16) );
            zgerv2d_( &ictxt, &nq, &ione, buff, &nq, &iyrow, &mycol );
            zrot_( &nq, &X[iix-1+(jjx-1)*desc_X[LLD_]], &desc_X[LLD_],
                        buff, &ione, c, s );
         }
         else if( myrow == iyrow )
         {
            zgesd2d_( &ictxt, &ione, &nq,
                      &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], &desc_Y[LLD_],
                      &ixrow, &mycol );
            buff = (complex16 *)getpbbuf( "PZROT", nq*sizeof(complex16) );
            zgerv2d_( &ictxt, &nq, &ione, buff, &nq, &ixrow, &mycol );
            zrot_( &nq, buff, &ione,
                        &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], &desc_Y[LLD_], c, s );
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
            zrot_( &np, &X[iix-1+(jjx-1)*desc_X[LLD_]], incx,
                        &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], incy, c, s );
         }
      }
      else
      {
         if( mycol == ixcol )
         {
            zgesd2d_( &ictxt, &np, &ione,
                      &X[iix-1+(jjx-1)*desc_X[LLD_]], &desc_X[LLD_],
                      &myrow, &iycol );
            buff = (complex16 *)getpbbuf( "PZROT", np*sizeof(complex16) );
            zgerv2d_( &ictxt, &np, &ione, buff, &np, &myrow, &iycol );
            zrot_( &np, &X[iix-1+(jjx-1)*desc_X[LLD_]], incx,
                        buff, &ione, c, s );
         }
         else if( mycol == iycol )
         {
            zgesd2d_( &ictxt, &np, &ione,
                      &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], &desc_Y[LLD_],
                      &myrow, &ixcol );
            buff = (complex16 *)getpbbuf( "PZROT", np*sizeof(complex16) );
            zgerv2d_( &ictxt, &np, &ione, buff, &np, &myrow, &ixcol );
            zrot_( &np, buff, &ione,
                        &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], incy, c, s );
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
 
         buff = (complex16 *)getpbbuf( "PZROT", wksz*sizeof(complex16) );
 
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
            zrot_( &np, &X[iix-1+(jjx-1)*desc_X[LLD_]],
                     incx, buff, &ione, c, s );
         }
         pbztrnv_( &ictxt, C2F_CHAR( "R" ), C2F_CHAR( "T" ), n,
                   &desc_Y[NB_], &nz, buff, &ione, &zero, 
                   &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], &desc_Y[LLD_], 
                   &ixrow, &ixcol,  &iyrow, &iycol, buff+np );
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
 
         buff = (complex16 *)getpbbuf( "PZROT", wksz*sizeof(complex16) );
 
         if( myrow == iyrow )
            np -= nz;
         pbztrnv_( &ictxt, C2F_CHAR( "R" ), C2F_CHAR( "T" ), n,
                   &desc_X[NB_], &nz, &X[iix-1+(jjx-1)*desc_X[LLD_]],
                   &desc_X[LLD_], &zero, buff, &ione, &ixrow, &ixcol,
                   &iyrow, &iycol, buff+np );
         if( mycol == iycol )
         {
            zrot_( &np, buff, &ione,
                     &Y[iiy-1+(jjy-1)*desc_Y[LLD_]], incy, c, s );
         }
         pbztrnv_( &ictxt, C2F_CHAR( "R" ), C2F_CHAR( "T" ), n,
                   &desc_X[NB_], &nz, buff, &ione, &zero,
                   &X[iix-1+(jjx-1)*desc_X[LLD_]], &desc_X[LLD_], 
                   &iyrow, &iycol, &ixrow, &ixcol, buff+np );
      }
   }
}
