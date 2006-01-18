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
void PB_Cchkvec( int ICTXT, char * ROUT, char * VNAME, int N, int NPOS0,
                 int IX, int JX, int * DESCX, int INCX, int DPOS0,
                 int * INFO )
#else
void PB_Cchkvec( ICTXT, ROUT, VNAME, N, NPOS0, IX, JX, DESCX, INCX,
                 DPOS0, INFO )
/*
*  .. Scalar Arguments ..
*/
   int            DPOS0, ICTXT, IX, * INFO, INCX, JX, N, NPOS0;
/*
*  .. Array Arguments ..
*/
   char           * ROUT, * VNAME;
   int            * DESCX;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cchkvec  checks the  validity of a  descriptor vector  DESCX,  the
*  related global indexes  IX,  JX  and the global increment INCX. If an
*  inconsistency is found among its parameters  IX,  JX, DESCX and INCX,
*  the routine returns an error code in INFO.
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  ROUT    (global input) pointer to CHAR
*          On entry, ROUT specifies the name of the routine calling this
*          input error checking routine.
*
*  VNAME   (global input) pointer to CHAR
*          On entry,  VNAME specifies the name of the formal array argu-
*          ment in the calling routine.
*
*  N       (global input) INTEGER
*          On entry,  N  specifies the length of the subvector sub( X ).
*
*  NPOS0   (global input) INTEGER
*          On entry,  NPOS0  specifies the  position in the calling rou-
*          tine's parameter list where the formal parameter N appears.
*
*  IX      (global input) INTEGER
*          On entry, IX  specifies X's global row index, which points to
*          the beginning of the submatrix sub( X ).
*
*  JX      (global input) INTEGER
*          On entry, JX  specifies X's global column index, which points
*          to the beginning of the submatrix sub( X ).
*
*  DESCX   (global and local input) INTEGER array
*          On entry, DESCX  is an integer array of dimension DLEN_. This
*          is the array descriptor for the matrix X.
*
*  INCX    (global input) INTEGER
*          On entry,  INCX   specifies  the  global  increment  for  the
*          elements of  X.  Only two values of  INCX   are  supported in
*          this version, namely 1 and M_X. INCX  must not be zero.
*
*  DPOS0   (global input) INTEGER
*          On entry,  DPOS0  specifies the  position in the calling rou-
*          tine's parameter list where the formal  parameter  DESCX  ap-
*          pears.  Note that it is assumed that  IX and JX are respecti-
*          vely 2 and 1 entries behind DESCX, and  INCX is 1 entry after
*          DESCX.
*
*  INFO    (local input/local output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had an
*                illegal  value,  then  INFO = -(i*100+j),  if  the i-th
*                argument is a  scalar  and had an  illegal  value, then
*                INFO = -i.
*
*  -- Written on April 1, 1998 by
*     R. Clint Whaley, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   int            dpos, icpos, ixpos, jxpos, mycol, myrow, np, npcol, npos,
                  nprow, nq;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Want to find errors with MIN(), so if no error, set it to a big number. If
*  there already is an error, multiply by the the descriptor multiplier.
*/
   if( *INFO >= 0 )             *INFO = BIGNUM;
   else if( *INFO < -DESCMULT ) *INFO = -(*INFO);
   else                         *INFO = -(*INFO) * DESCMULT;
/*
*  Figure where in parameter list each parameter was, factoring in descriptor
*  multiplier
*/
   npos  = NPOS0 * DESCMULT;
   ixpos = ( DPOS0 - 2 ) * DESCMULT;
   jxpos = ( DPOS0 - 1 ) * DESCMULT;
   icpos = ( DPOS0 + 1 ) * DESCMULT;
   dpos  = DPOS0 * DESCMULT + 1;
/*
*  Get process grid information
*/
   Cblacs_gridinfo( ICTXT, &nprow, &npcol, &myrow, &mycol );
/*
*  Are N, IX, JX, DESCX and INCX legal inputs ?
*/
   if( N < 0 )
   {
/*
*  N must be at least zero
*/
      *INFO = MIN( *INFO, npos );
      PB_Cwarn( ICTXT, -1, ROUT, "%s sub( %s ) = %d, it must be at least 0",
                "Illegal length of", VNAME, N );
   }

   if( IX < 0 )
   {
/*
*  IX must be at least zero
*/
      *INFO = MIN( *INFO, ixpos );
      PB_Cwarn( ICTXT, -1, ROUT, "Illegal I%s = %d, I%s must be at least 1",
                VNAME, IX+1, VNAME );
   }
   if( JX < 0 )
   {
/*
*  JX must be at least zero
*/
      *INFO = MIN( *INFO, jxpos );
      PB_Cwarn( ICTXT, -1, ROUT, "Illegal J%s = %d, J%s must be at least 1",
                VNAME, JX+1, VNAME );
   }

   if( DESCX[DTYPE_] != BLOCK_CYCLIC_2D_INB )
   {
/*
*  Internally, only the descriptor type BLOCK_CYCLIC_2D_INB is supported.
*/
      *INFO = MIN( *INFO, dpos + DTYPE_ );
      PB_Cwarn( ICTXT, -1, ROUT, "%s %d for matrix %s. PBLAS accepts: %d or %d",
                "Illegal descriptor type", DESCX[DTYPE_], VNAME,
                BLOCK_CYCLIC_2D, BLOCK_CYCLIC_2D_INB );
      if( *INFO % DESCMULT == 0 ) *INFO = -( (*INFO) / DESCMULT );
      else                        *INFO = -(*INFO);
/*
*  No need to go any further ...
*/
      return;
   }

   if( DESCX[CTXT_] != ICTXT )
   {
/*
*  Check if the context of X match the other contexts. Only intra-context
*  operations are supported.
*/
      *INFO = MIN( *INFO, dpos + CTXT_ );
      PB_Cwarn( ICTXT, -1, ROUT, "DESC%s[CTXT_] = %d %s= %d", VNAME,
                DESCX[CTXT_], "does not match other operand's context ",
                ICTXT );
      if( *INFO % DESCMULT == 0 ) *INFO = -( (*INFO) / DESCMULT );
      else                        *INFO = -(*INFO);
/*
*  No need to go any further ...
*/
      return;
   }

   if( DESCX[IMB_] < 1 )
   {
/*
*  DESCX[IMB_] must be at least one
*/
      *INFO = MIN( *INFO, dpos + IMB_ );
      PB_Cwarn( ICTXT, -1, ROUT, "Illegal DESC%s[IMB_] = %d, DESC%s[IMB_] %s",
                VNAME, DESCX[IMB_], VNAME, "must be at least 1" );
   }
   if( DESCX[INB_] < 1 )
   {
/*
*  DESCX[INB_] must be at least one
*/
      *INFO = MIN( *INFO, dpos + INB_ );
      PB_Cwarn( ICTXT, -1, ROUT, "Illegal DESC%s[INB_] = %d, DESC%s[INB_] %s",
                VNAME, DESCX[INB_], VNAME, "must be at least 1" );
   }
   if( DESCX[MB_] < 1 )
   {
/*
*  DESCX[MB_] must be at least one
*/
      *INFO = MIN( *INFO, dpos + MB_ );
      PB_Cwarn( ICTXT, -1, ROUT, "Illegal DESC%s[MB_] = %d, DESC%s[MB_] %s",
                VNAME, DESCX[MB_], VNAME, "must be at least 1" );
   }
   if( DESCX[NB_] < 1 )
   {
/*
*  DESCX[NB_] must be at least one
*/
      *INFO = MIN( *INFO, dpos + NB_ );
      PB_Cwarn( ICTXT, -1, ROUT, "Illegal DESC%s[NB_] = %d, DESC%s[NB_] %s",
                VNAME, DESCX[NB_], VNAME, "must be at least 1" );
   }

   if( ( DESCX[RSRC_] < -1 ) || ( DESCX[RSRC_] >= nprow ) )
   {
/*
*  DESCX[RSRC_] must be either -1 (replication) or in the interval [0 .. nprow)
*/
      *INFO = MIN( *INFO, dpos + RSRC_ );
      PB_Cwarn( ICTXT, -1, ROUT,
                "Illegal DESC%s[RSRC_] = %d, DESC%s[RSRC_] %s%d", VNAME,
                DESCX[RSRC_], VNAME, "must be either -1, or >= 0 and < ",
                nprow );
   }
   if( ( DESCX[CSRC_] < -1 ) || ( DESCX[CSRC_] >= npcol ) )
   {
/*
*  DESCX[CSRC_] must be either -1 (replication) or in the interval [0 .. npcol)
*/
      *INFO = MIN( *INFO, dpos + CSRC_ );
      PB_Cwarn( ICTXT, -1, ROUT,
                "Illegal DESC%s[CSRC_] = %d, DESC%s[CSRC_] %s%d", VNAME,
                DESCX[CSRC_], VNAME, "must be either -1, or >= 0 and < ",
                npcol );
   }

   if( INCX != 1 && INCX != DESCX[M_] )
   {
/*
*  INCX must be either 1 or DESCX[M_]
*/
      *INFO = MIN( *INFO, icpos );
      PB_Cwarn( ICTXT, -1, ROUT,
                "Illegal INC%s = %d, INC%s should be either 1 or %d", VNAME,
                DESCX[M_], VNAME );
   }

   if( N == 0 )
   {
/*
*  NULL vector, relax some checks
*/
      if( DESCX[M_] < 0 )
      {
/*
*  DESCX[M_] must be at least 0
*/
         *INFO = MIN( *INFO, dpos + M_ );
         PB_Cwarn( ICTXT, -1, ROUT, "DESC%s[M_] = %d, it must be at least 0",
                   VNAME, DESCX[M_] );

      }
      if( DESCX[N_] < 0 )
      {
/*
*  DESCX[N_] must be at least 0
*/
         *INFO = MIN( *INFO, dpos + N_ );
         PB_Cwarn( ICTXT, -1, ROUT, "DESC%s[N_] = %d, it must be at least 0",
                   VNAME, DESCX[N_] );
      }

      if( DESCX[LLD_] < 1 )
      {
/*
*  DESCX[LLD_] must be at least 1
*/
         *INFO = MIN( *INFO, dpos + LLD_ );
         PB_Cwarn( ICTXT, -1, ROUT, "DESC%s[LLD_] = %d, it must be at least 1",
                   VNAME, DESCX[LLD_] );
      }
   }
   else
   {
/*
*  more rigorous checks for non-degenerate vector
*/
      if( DESCX[M_] < 1 )
      {
/*
*  DESCX[M_] must be at least 1
*/
         *INFO = MIN( *INFO, dpos + M_ );
         PB_Cwarn( ICTXT, -1, ROUT,
                   "Illegal DESC%s[M_] = %d, it must be at least 1", VNAME,
                   DESCX[M_]);

      }
      if( DESCX[N_] < 1 )
      {
/*
*  DESCX[N_] must be at least 1
*/
         *INFO = MIN( *INFO, dpos + N_ );
         PB_Cwarn( ICTXT, -1, ROUT,
                   "Illegal DESC%s[N_] = %d, it must be at least 1", VNAME,
                   DESCX[N_]);
      }

      if( ( DESCX[M_] >= 1 ) && ( DESCX[N_] >= 1 ) )
      {
         if( INCX == DESCX[M_] )
         {
/*
*  sub( X ) resides in (a) process row(s)
*/
            if( IX >= DESCX[M_] )
            {
/*
*  IX must be in [ 0 ... DESCX[M_]-1 ]
*/
               *INFO = MIN( *INFO, ixpos );
               PB_Cwarn( ICTXT, -1, ROUT, "%s I%s = %d, DESC%s[M_] = %d",
                         "Array subscript out of bounds:", VNAME, IX+1, VNAME,
                         DESCX[M_]);
            }
            if( JX+N > DESCX[N_] )
            {
/*
*  JX + N must be in [ 0 ... DESCX[N_]-1 ]
*/
               *INFO = MIN( *INFO, jxpos );
               PB_Cwarn( ICTXT, -1, ROUT,
                         "%s N = %d, J%s = %d, DESC%s[N_] = %d",
                         "Operation out of bounds:", N, VNAME, JX+1, VNAME,
                         DESCX[N_]);
            }
         }
         else
         {
/*
*  sub( X ) resides in (a) process column(s)
*/
            if( JX >= DESCX[N_] )
            {
/*
*  JX must be in [ 0 ... DESCX[N_] ]
*/
               *INFO = MIN( *INFO, jxpos );
               PB_Cwarn( ICTXT, -1, ROUT, "%s J%s = %d, DESC%s[N_] = %d",
                         "Array subscript out of bounds:", VNAME, JX+1, VNAME,
                         DESCX[N_]);
            }
            if( IX+N > DESCX[M_] )
            {
/*
*  IX + N must be in [ 0 ... DESCX[M_] ]
*/
               *INFO = MIN( *INFO, ixpos );
               PB_Cwarn( ICTXT, -1, ROUT,
                         "%s N = %d, I%s = %d, DESC%s[M_] = %d",
                         "Operation out of bounds:", N, VNAME, IX+1, VNAME,
                         DESCX[M_]);
            }
         }
      }
/*
*  *INFO == BIGNUM => No errors have been found so far
*/
      if( *INFO == BIGNUM )
      {
         Mnumroc( np, DESCX[M_], 0, DESCX[IMB_], DESCX[MB_], myrow,
                  DESCX[RSRC_], nprow );
         if( DESCX[LLD_] < MAX( 1, np ) )
         {
            Mnumroc( nq, DESCX[N_], 0, DESCX[INB_], DESCX[NB_], mycol,
                     DESCX[CSRC_], npcol );
/*
*  DESCX[LLD_] must be at least 1 in order to be legal and this is enough if no
*  columns of X reside in this process
*/
            if( DESCX[LLD_] < 1 )
            {
               *INFO = MIN( *INFO, dpos + LLD_ );
               PB_Cwarn( ICTXT, -1, ROUT,
                         "DESC%s[LLD_] = %d, it must be at least 1", VNAME,
                         DESCX[LLD_] );
            }
            else if( nq > 0 )
            {
/*
*  Some columns of X reside in this process, DESCX[LLD_] must be at least
*  MAX( 1, np )
*/
               *INFO = MIN( *INFO, dpos + LLD_ );
               PB_Cwarn( ICTXT, -1, ROUT,
                         "DESC%s[LLD_] = %d, it must be at least %d", VNAME,
                         DESCX[LLD_], MAX( 1, np ) );
            }
         }
      }
   }
/*
*  Prepare output: set INFO = 0 if no error, and divide by DESCMULT if error is
*  not in a descriptor entry.
*/
   if( *INFO == BIGNUM )            *INFO = 0;
   else if( *INFO % DESCMULT == 0 ) *INFO = -( (*INFO) / DESCMULT );
   else                             *INFO = -(*INFO);
/*
*  End of PB_Cchkvec
*/
}
