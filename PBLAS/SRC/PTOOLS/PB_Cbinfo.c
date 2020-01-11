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
void PB_Cbinfo( Int OFFD, Int M, Int N, Int IMB1, Int INB1, Int MB,
                Int NB, Int MRROW, Int MRCOL, Int * LCMT00, Int * MBLKS,
                Int * NBLKS, Int * IMBLOC, Int * INBLOC, Int * LMBLOC,
                Int * LNBLOC, Int * ILOW, Int * LOW, Int * IUPP, Int * UPP )
#else
void PB_Cbinfo( OFFD, M, N, IMB1, INB1, MB, NB, MRROW, MRCOL, LCMT00,
                MBLKS, NBLKS, IMBLOC, INBLOC, LMBLOC, LNBLOC, ILOW, LOW,
                IUPP, UPP )
/*
*  .. Scalar Arguments ..
*/
   Int            * ILOW, IMB1, * IMBLOC, INB1, * INBLOC, * IUPP,
                  * LCMT00, * LMBLOC, * LNBLOC, * LOW, M, MB, * MBLKS,
                  MRCOL, MRROW, N, NB, * NBLKS, OFFD, * UPP;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cbinfo  initializes the local information of an m by n local array
*  owned by the process of  relative  coordinates ( MRROW, MRCOL ). Note
*  that if m or n is less or equal than zero, there is no data, in which
*  case this process  does  not  need  the local information computed by
*  this routine to proceed.
*
*  Arguments
*  =========
*
*  OFFD    (global input) INTEGER
*          On entry,  OFFD  specifies the off-diagonal of the underlying
*          matrix of interest as follows:
*             OFFD = 0 specifies the main diagonal,
*             OFFD > 0 specifies lower subdiagonals, and
*             OFFD < 0 specifies upper superdiagonals.
*
*  M       (local input) INTEGER
*          On entry, M  specifies the local number of rows of the under-
*          lying matrix  owned  by the  process  of relative coordinates
*          ( MRROW, MRCOL ). M must be at least zero.
*
*  N       (local input) INTEGER
*          On entry, N  specifies the local number of columns of the un-
*          derlying matrix  owned by the process of relative coordinates
*          ( MRROW, MRCOL ). N must be at least zero.
*
*  IMB1    (global input) INTEGER
*          On input, IMB1 specifies  the global true size of  the  first
*          block of rows of the underlying global submatrix.  IMB1  must
*          be at least MIN( 1, M ).
*
*  INB1    (global input) INTEGER
*          On input, INB1 specifies  the global true size of  the  first
*          block  of  columns  of  the underlying global submatrix. INB1
*          must be at least MIN( 1, N ).
*
*  MB      (global input) INTEGER
*          On entry, MB  specifies the blocking factor used to partition
*          the rows of the matrix.  MB  must be at least one.
*
*  NB      (global input) INTEGER
*          On entry, NB  specifies the blocking factor used to partition
*          the the columns of the matrix.  NB  must be at least one.
*
*  MRROW   (local input) INTEGER
*          On entry, MRROW specifies the  relative row coordinate of the
*          process that possesses these M rows. MRROW must be least zero
*          and strictly less than NPROW.
*
*  MRCOL   (local input) INTEGER
*          On entry, MRCOL specifies  the  relative column coordinate of
*          the process that possesses these N  columns.  MRCOL  must  be
*          least zero and strictly less than NPCOL.
*
*  LCMT00  (local output) INTEGER
*          On exit, LCMT00  is the  LCM value of the left upper block of
*          this m by n local  block owned by the process of relative co-
*          ordinates ( MRROW, MRCOL ).
*
*  MBLKS   (local output) INTEGER
*          On exit, MBLKS specifies the local number of blocks  of  rows
*          corresponding to M. MBLKS must be at least zero.
*
*  NBLKS   (local output) INTEGER
*          On exit,  NBLKS  specifies  the local number of blocks of co-
*          lumns corresponding to N. NBLKS must be at least zero.
*
*  IMBLOC  (local output) INTEGER
*          On exit, IMBLOC  specifies  the  number of rows (size) of the
*          uppest blocks of this m by n local array owned by the process
*          of relative coordinates ( MRROW, MRCOL ).  IMBLOC is at least
*          MIN( 1, M ).
*
*  INBLOC  (local output) INTEGER
*          On exit, INBLOC  specifies  the  number of columns (size) of
*          the leftmost  blocks of this m by n local array owned by the
*          process of relative coordinates ( MRROW, MRCOL ).  INBLOC is
*          at least MIN( 1, N ).
*
*  LMBLOC  (local output) INTEGER
*          On exit, LMBLOC specifies the number  of  rows  (size) of the
*          lowest blocks of this m by n local array owned by the process
*          of  relative coordinates ( MRROW, MRCOL ). LMBLOC is at least
*          MIN( 1, M ).
*
*  LNBLOC  (local output) INTEGER
*          On exit, LNBLOC specifies the number of columns (size) of the
*          rightmost  blocks of this  m by n  local  array  owned by the
*          process of  relative  coordinates ( MRROW, MRCOL ). LNBLOC is
*          at least MIN( 1, N ).
*
*  ILOW    (local output) INTEGER
*          On exit, ILOW is the lower bound characterizing the first co-
*          lumn block owning offdiagonals of  this  m by n  array.  ILOW
*          must be less or equal than zero.
*
*  LOW     (global output) INTEGER
*          On exit,  LOW  is  the  lower bound characterizing the column
*          blocks with te exception of the  first  one (see ILOW) owning
*          offdiagonals of this m by n array. LOW  must be less or equal
*          than zero.
*
*  IUPP    (local output) INTEGER
*          On exit, IUPP is the upper bound characterizing the first row
*          block owning offdiagonals of this m by n array.  IUPP must be
*          greater or equal than zero.
*
*  UPP     (global output) INTEGER
*          On exit,  UPP  is  the  upper  bound  characterizing  the row
*          blocks with te exception of the  first  one (see IUPP) owning
*          offdiagonals of this m by n array. UPP  must  be  greater  or
*          equal than zero.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   Int            tmp1;
/* ..
*  .. Executable Statements ..
*
*/
/*
*  Initialize LOW, ILOW, UPP, IUPP, LMBLOC, LNBLOC, IMBLOC, INBLOC, MBLKS,
*  NBLKS and LCMT00.
*/
   *LOW = 1 - NB;
   *UPP = MB - 1;

   *LCMT00 = OFFD;

   if( ( M <= 0 ) || ( N <= 0 ) )
   {
/*
*  If the local virtual array is empty, then simplify the remaining of the
*  initialization.
*/
      *IUPP   = ( MRROW ? MB - 1 : ( IMB1 > 0 ? IMB1 - 1 : 0 ) );
      *IMBLOC = 0;
      *MBLKS  = 0;
      *LMBLOC = 0;

      *ILOW   = ( MRCOL ? 1 - NB : ( INB1 > 0 ? 1 - INB1 : 0 ) );
      *INBLOC = 0;
      *NBLKS  = 0;
      *LNBLOC = 0;

      *LCMT00 += ( *LOW - *ILOW + MRCOL * NB ) - ( *IUPP - *UPP + MRROW * MB );

      return;
   }

   if( MRROW )
   {
/*
*  I am not in the first relative process row. Use the first local row block
*  size MB to initialize the VM structure.
*/
      *IMBLOC  = MIN( M, MB );
      *IUPP    = MB - 1;
      *LCMT00 -= IMB1 - MB + MRROW * MB;
      *MBLKS   = ( M - 1 ) / MB + 1;
      *LMBLOC  = M - ( M / MB ) * MB;
      if( !( *LMBLOC ) ) *LMBLOC = MB;

      if( MRCOL )
      {
/*
*  I am not in the first relative process column. Use the first local column
*  block size NB to initialize the VM structure.
*/
         *INBLOC  = MIN( N, NB );
         *ILOW    = 1 - NB;
         *LCMT00 += INB1 - NB + MRCOL * NB;
         *NBLKS   = ( N - 1 ) / NB + 1;
         *LNBLOC  = N - ( N / NB ) * NB;
         if( !( *LNBLOC ) ) *LNBLOC = NB;
      }
      else
      {
/*
*  I am in the first relative process column. Use the first column block size
*  INB1 to initialize the VM structure.
*/
         *INBLOC = INB1;
         *ILOW  = 1 - INB1;
         tmp1   = N - INB1;
         if( tmp1 )
         {
/*
*  There is more than one column block. Compute the number of local column
*  blocks and the size of the last one.
*/
            *NBLKS  = ( tmp1 - 1 ) / NB + 2;
            *LNBLOC = tmp1 - ( tmp1 / NB ) * NB;
            if( !( *LNBLOC ) ) *LNBLOC = NB;
         }
         else
         {
/*
*  There is only one column block.
*/
            *NBLKS  = 1;
            *LNBLOC = INB1;
         }
      }
   }
   else
   {
/*
*  I am in the first relative process row. Use the first row block size IMB1 to
*  initialize the VM structure.
*/
      *IMBLOC = IMB1;
      *IUPP   = IMB1 - 1;
      tmp1    = M - IMB1;
      if( tmp1 )
      {
/*
*  There is more than one row block. Compute the number of local row blocks and
*  the size of the last one.
*/
         *MBLKS  = ( tmp1 - 1 ) / MB + 2;
         *LMBLOC = tmp1 - ( tmp1 / MB ) * MB;
         if( !( *LMBLOC ) ) *LMBLOC = MB;
      }
      else
      {
/*
*  There is only one row block.
*/
         *MBLKS  = 1;
         *LMBLOC = IMB1;
      }

      if( MRCOL )
      {
/*
*  I am not in the first relative process column. Use the first local column
*  block size NB to initialize the VM structure.
*/
         *INBLOC  = MIN( N, NB );
         *ILOW    = 1 - NB;
         *LCMT00 += INB1 - NB + MRCOL * NB;
         *NBLKS   = ( N - 1 ) / NB + 1;
         *LNBLOC  = N - ( N / NB ) * NB;
         if( !( *LNBLOC ) ) *LNBLOC = NB;
      }
      else
      {
/*
*  I am in the first relative process column. Use the first column block size
*  INB1 to initialize the VM structure.
*/
         *INBLOC = INB1;
         *ILOW   = 1 - INB1;
         tmp1    = N - INB1;
         if( tmp1 )
         {
/*
*  There is more than one column block. Compute the number of local column
*  blocks and the size of the last one.
*/
            *NBLKS  = ( tmp1 - 1 ) / NB + 2;
            *LNBLOC = tmp1 - ( tmp1 / NB ) * NB;
            if( !( *LNBLOC ) ) *LNBLOC = NB;
         }
         else
         {
/*
*  There is only one column block.
*/
            *NBLKS  = 1;
            *LNBLOC = INB1;
         }
      }
   }
/*
*  End of PB_Cbinfo
*/
}
