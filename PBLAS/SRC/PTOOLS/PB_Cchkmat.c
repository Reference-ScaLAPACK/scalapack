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
void PB_Cchkmat( int ICTXT, char * ROUT, char * MNAME, int M, int MPOS0,
                 int N, int NPOS0, int IA, int JA, int * DESCA, int DPOS0,
                 int * INFO )
#else
void PB_Cchkmat( ICTXT, ROUT, MNAME, M, MPOS0, N, NPOS0, IA, JA, DESCA,
                 DPOS0, INFO )
/*
*  .. Scalar Arguments ..
*/
   int            DPOS0, IA, ICTXT, * INFO, JA, M, MPOS0, N, NPOS0;
/*
*  .. Array Arguments ..
*/
   char           * MNAME, * ROUT;
   int            * DESCA;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cchkmat  checks the  validity of a  descriptor vector  DESCA,  the
*  related global indexes  IA, JA  from a local view point. If an incon-
*  sistency is found among its parameters IA, JA and DESCA, the  routine
*  returns an error code in INFO.
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
*  MNAME   (global input) pointer to CHAR
*          On entry,  MNAME specifies the name of the formal array argu-
*          ment in the calling routine.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies  the  number  of  rows  the submatrix
*          sub( A ).
*
*  MPOS0   (global input) INTEGER
*          On entry,  MPOS0  specifies the  position in the calling rou-
*          tine's parameter list where the formal parameter M appears.
*
*  N       (global input) INTEGER
*          On entry,  N  specifies  the  number of columns the submatrix
*          sub( A ).
*
*  NPOS0   (global input) INTEGER
*          On entry,  NPOS0  specifies the  position in the calling rou-
*          tine's parameter list where the formal parameter N appears.
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
*  DPOS0   (global input) INTEGER
*          On entry,  DPOS0  specifies the  position in the calling rou-
*          tine's parameter list where the formal  parameter  DESCA  ap-
*          pears.  Note that it is assumed that  IA and JA are respecti-
*          vely 2 and 1 entries behind DESCA.
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
   int            dpos, iapos, japos, mpos, mycol, myrow, np, npcol, nprow,
                  npos, nq;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Want to find errors with MIN( ), so if no error, set it to a big number. If
*  there already is an error, multiply by the the descriptor multiplier.
*/
   if( *INFO >= 0 )             *INFO = BIGNUM;
   else if( *INFO < -DESCMULT ) *INFO = -(*INFO);
   else                         *INFO = -(*INFO) * DESCMULT;
/*
*  Figure where in parameter list each parameter was, factoring in descriptor
*  multiplier
*/
   mpos  = MPOS0 * DESCMULT;
   npos  = NPOS0 * DESCMULT;
   iapos = ( DPOS0 - 2 ) * DESCMULT;
   japos = ( DPOS0 - 1 ) * DESCMULT;
   dpos  = DPOS0 * DESCMULT + 1;
/*
*  Get process grid information
*/
   Cblacs_gridinfo( ICTXT, &nprow, &npcol, &myrow, &mycol );
/*
*  Are M, N, IA, JA and DESCA legal inputs ?
*/
   if( M < 0 )
   {
/*
*  M must be at least zero
*/
      *INFO = MIN( *INFO, mpos );
      PB_Cwarn( ICTXT, -1, ROUT, "%s sub( %s ) = %d, it must be at least 0",
                "Illegal number of rows of", MNAME, M );
   }
   if( N < 0 )
   {
/*
*  N must be at least zero
*/
      *INFO = MIN( *INFO, npos );
      PB_Cwarn( ICTXT, -1, ROUT, "%s sub( %s ) = %d, it must be at least 0",
                "Illegal number of columns of", MNAME, N );
   }

   if( IA < 0 )
   {
/*
*  IA must be at least zero
*/
      *INFO = MIN( *INFO, iapos );
      PB_Cwarn( ICTXT, -1, ROUT, "Illegal I%s = %d, I%s must be at least 1",
                MNAME, IA+1, MNAME );
   }
   if( JA < 0 )
   {
/*
*  JA must be at least zero
*/
      *INFO = MIN( *INFO, japos );
      PB_Cwarn( ICTXT, -1, ROUT, "Illegal J%s = %d, I%s must be at least 1",
                MNAME, IA+1, MNAME );
   }

   if( DESCA[DTYPE_] != BLOCK_CYCLIC_2D_INB )
   {
/*
*  Internally, only the descriptor type BLOCK_CYCLIC_2D_INB is supported
*/
      *INFO = MIN( *INFO, dpos + DTYPE_ );
      PB_Cwarn( ICTXT, -1, ROUT, "%s %d for matrix %s. PBLAS accepts: %d or %d",
                "Illegal descriptor type", DESCA[DTYPE_], MNAME,
                BLOCK_CYCLIC_2D, BLOCK_CYCLIC_2D_INB );
      if( *INFO % DESCMULT == 0 ) *INFO = -( (*INFO) / DESCMULT );
      else                        *INFO = -(*INFO);
/*
*  No need to go any further ...
*/
      return;
   }

   if( DESCA[CTXT_] != ICTXT )
   {
/*
*  Check if the context of X match the other contexts. Only intra-context
*  operations are supported.
*/
      *INFO = MIN( *INFO, dpos + CTXT_ );
      PB_Cwarn( ICTXT, -1, ROUT, "DESC%s[CTXT_] = %d %s= %d", MNAME,
                DESCA[CTXT_], "does not match other operand's context ",
                ICTXT );
      if( *INFO % DESCMULT == 0 ) *INFO = -( (*INFO) / DESCMULT );
      else                        *INFO = -(*INFO);
/*
*  No need to go any further ...
*/
      return;
   }

   if( DESCA[IMB_] < 1 )
   {
/*
*  DESCA[IMB_] must be at least one
*/
      *INFO = MIN( *INFO, dpos + IMB_ );
      PB_Cwarn( ICTXT, -1, ROUT, "Illegal DESC%s[IMB_] = %d, DESC%s[IMB_] %s",
                MNAME, DESCA[IMB_], MNAME, "must be at least 1" );
   }
   if( DESCA[INB_] < 1 )
   {
/*
*  DESCA[INB_] must be at least one
*/
      *INFO = MIN( *INFO, dpos + INB_ );
      PB_Cwarn( ICTXT, -1, ROUT, "Illegal DESC%s[INB_] = %d, DESC%s[INB_] %s",
                MNAME, DESCA[INB_], MNAME, "must be at least 1" );
   }
   if( DESCA[MB_] < 1 )
   {
/*
*  DESCA[MB_] must be at least one
*/
      *INFO = MIN( *INFO, dpos + MB_ );
      PB_Cwarn( ICTXT, -1, ROUT, "Illegal DESC%s[MB_] = %d, DESC%s[MB_] %s",
                MNAME, DESCA[MB_], MNAME, "must be at least 1" );
   }
   if( DESCA[NB_] < 1 )
   {
/*
*  DESCA[NB_] must be at least one
*/
      *INFO = MIN( *INFO, dpos + NB_ );
      PB_Cwarn( ICTXT, -1, ROUT, "Illegal DESC%s[NB_] = %d, DESC%s[NB_] %s",
                MNAME, DESCA[NB_], MNAME, "must be at least 1" );
   }

   if( ( DESCA[RSRC_] < -1 ) || ( DESCA[RSRC_] >= nprow ) )
   {
/*
*  DESCA[RSRC_] must be either -1 (replication) or in the interval [0 .. nprow)
*/
      *INFO = MIN( *INFO, dpos + RSRC_ );
      PB_Cwarn( ICTXT, -1, ROUT,
                "Illegal DESC%s[RSRC_] = %d, DESC%s[RSRC_] %s%d", MNAME,
                DESCA[RSRC_], MNAME, "must be either -1, or >= 0 and < ",
                nprow );
   }
   if( ( DESCA[CSRC_] < -1 ) || ( DESCA[CSRC_] >= npcol ) )
   {
/*
*  DESCX[CSRC_] must be either -1 (replication) or in the interval [0 .. npcol)
*/
      *INFO = MIN( *INFO, dpos + CSRC_ );
      PB_Cwarn( ICTXT, -1, ROUT,
                "Illegal DESC%s[CSRC_] = %d, DESC%s[CSRC_] %s%d", MNAME,
                DESCA[CSRC_], MNAME, "must be either -1, or >= 0 and < ",
                npcol );
   }

   if( M == 0 || N == 0 )
   {
/*
*  NULL matrix, relax some checks
*/
      if( DESCA[M_] < 0 )
      {
/*
*  DESCX[M_] must be at least 0
*/
         *INFO = MIN( *INFO, dpos + M_ );
         PB_Cwarn( ICTXT, -1, ROUT, "DESC%s[M_] = %d, it must be at least 0",
                   MNAME, DESCA[M_] );
      }
      if( DESCA[N_] < 0 )
      {
/*
*  DESCX[N_] must be at least 0
*/
         *INFO = MIN( *INFO, dpos + N_ );
         PB_Cwarn( ICTXT, -1, ROUT, "DESC%s[N_] = %d, it must be at least 0",
                   MNAME, DESCA[N_] );
      }

      if( DESCA[LLD_] < 1 )
      {
/*
*  DESCA[LLD_] must be at least 1
*/
         *INFO = MIN( *INFO, dpos + LLD_ );
         PB_Cwarn( ICTXT, -1, ROUT, "DESC%s[LLD_] = %d, it must be at least 1",
                   MNAME, DESCA[LLD_] );
      }
   }
   else
   {
/*
*  more rigorous checks for non-degenerate matrix
*/
      if( DESCA[M_] < 1 )
      {
/*
*  DESCA[M_] must be at least 1
*/
         *INFO = MIN( *INFO, dpos + M_ );
         PB_Cwarn( ICTXT, -1, ROUT,
                   "Illegal DESC%s[M_] = %d, it must be at least 1", MNAME,
                   DESCA[M_]);
      }
      if( DESCA[N_] < 1 )
      {
/*
*  DESCA[N_] must be at least 1
*/
         *INFO = MIN( *INFO, dpos + N_ );
         PB_Cwarn( ICTXT, -1, ROUT,
                   "Illegal DESC%s[N_] = %d, it must be at least 1", MNAME,
                   DESCA[N_]);
      }

      if( ( DESCA[M_] >= 1 ) && ( DESCA[N_] >= 1 ) )
      {
         if( IA+M > DESCA[M_] )
         {
/*
*  IA + M must be in [ 0 ... DESCA[M_] ]
*/
            *INFO = MIN( *INFO, iapos );
            PB_Cwarn( ICTXT, -1, ROUT, "%s M = %d, I%s = %d, DESC%s[M_] = %d",
                      "Operation out of bounds:", M, MNAME, IA+1, MNAME,
                      DESCA[M_]);
         }
         if( JA+N > DESCA[N_] )
         {
/*
*  JA + N must be in [ 0 ... DESCA[N_] ]
*/
            *INFO = MIN( *INFO, japos );
            PB_Cwarn( ICTXT, -1, ROUT, "%s N = %d, J%s = %d, DESC%s[N_] = %d",
                      "Operation out of bounds:", N, MNAME, JA+1, MNAME,
                      DESCA[N_]);
         }
      }
/*
*  *INFO == BIGNUM => No errors have been found so far
*/
      if( *INFO == BIGNUM )
      {
         Mnumroc( np, DESCA[M_], 0, DESCA[IMB_], DESCA[MB_], myrow,
                  DESCA[RSRC_], nprow );
         if( DESCA[LLD_] < MAX( 1, np ) )
         {
            Mnumroc( nq, DESCA[N_], 0, DESCA[INB_], DESCA[NB_], mycol,
                     DESCA[CSRC_], npcol );
/*
*  DESCA[LLD_] must be at least 1 in order to be legal and this is enough if no
*  columns of A reside in this process.
*/
            if( DESCA[LLD_] < 1 )
            {
               *INFO = MIN( *INFO, dpos + LLD_ );
               PB_Cwarn( ICTXT, -1, ROUT,
                         "DESC%s[LLD_] = %d, it must be at least 1", MNAME,
                         DESCA[LLD_] );
            }
            else if( nq > 0 )
            {
/*
*  Some columns of A reside in this process, DESCA[LLD_] must be at least
*  MAX( 1, np ).
*/
               *INFO = MIN( *INFO, dpos + LLD_ );
               PB_Cwarn( ICTXT, -1, ROUT,
                         "DESC%s[LLD_] = %d, it must be at least %d", MNAME,
                         DESCA[LLD_], MAX( 1, np ) );
            }
         }
      }
   }
/*
*  Prepare output: set info = 0 if no error, and divide by DESCMULT if error is
*  not in a descriptor entry.
*/
   if( *INFO == BIGNUM )            *INFO = 0;
   else if( *INFO % DESCMULT == 0 ) *INFO = -( (*INFO) / DESCMULT );
   else                             *INFO = -(*INFO);
/*
*  End of PB_Cchkmat
*/
}
