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
void PB_Cdescribe( Int M, Int N, Int IA, Int JA, Int * DA, Int NPROW,
                   Int NPCOL, Int MYROW, Int MYCOL, Int * II, Int * JJ,
                   Int * LDA, Int * IMB, Int * INB, Int * MB, Int * NB,
                   Int * PROW, Int * PCOL, Int * DA0 )
#else
void PB_Cdescribe( M, N, IA, JA, DA, NPROW, NPCOL, MYROW, MYCOL, II, JJ,
                   LDA, IMB, INB, MB, NB, PROW, PCOL, DA0 )
/*
*  .. Scalar Arguments ..
*/
   Int            IA, *II, * IMB, * INB, JA, * JJ, * LDA, M, * MB, MYCOL,
                  MYROW, N, * NB, NPCOL, NPROW, * PCOL, * PROW;
/*
*  .. Array Arguments ..
*/
   Int            * DA, * DA0;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cdescribe  returns the global descriptor of a submatrix. This rou-
*  tine also computes the starting local index II, JJ  corresponding  to
*  the submatrix starting globally at the entry pointed by IA, JA.  This
*  routine returns the coordinates in the grid of the process owning the
*  matrix entry of global indexes I, J, namely PROW and PCOL.  The  true
*  global block sizes IMB, INB, MB and NB are also returned.
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
*  M       (global input) INTEGER
*          On entry, M  specifies the number of  rows  being  dealt  out
*          starting from global index IA. M  is  also the number of rows
*          of the submatrix of interest. M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry, N  specifies the number of columns being dealt  out
*          starting from global index JA. N  is  also the number of col-
*          umns of the submatrix of interest. M must be at least zero.
*
*  IA      (global input) INTEGER
*          On entry, IA specifies  the  global starting row index of the
*          submatrix. IA must at least zero.
*
*  JA      (global input) INTEGER
*          On entry, JA specifies  the  global  starting column index of
*          the submatrix. JA must at least zero.
*
*  DA      (global and local input) INTEGER array
*          On entry, DA  is an integer array of dimension DLEN_. This is
*          the array descriptor of the underlying matrix.
*
*  NPROW   (global input) INTEGER
*          On entry,  NPROW   specifies the total number of process rows
*          over which the matrix is distributed.  NPROW must be at least
*          one.
*
*  NPCOL   (global input) INTEGER
*          On entry, NPCOL specifies the total number of process columns
*          over which the matrix is distributed.  NPCOL must be at least
*          one.
*
*  MYROW   (local input) INTEGER
*          On entry,  MYROW  specifies the row coordinate of the process
*          whose local index  II  is determined.  MYROW must be at least
*          zero and strictly less than NPROW.
*
*  II      (local output) INTEGER
*          On exit, II  specifies the  local  starting  row index of the
*          submatrix. On exit, II is at least zero.
*
*  JJ      (local output) INTEGER
*          On exit, JJ  specifies the local starting column index of the
*          submatrix. On exit, JJ is at least zero.
*
*  LDA     (local output) INTEGER
*          On exit, LDA specifies the local leading dimension of the lo-
*          cal array containing the distributed matrix entries. LDA must
*          be at least one.
*
*  IMB     (global output) INTEGER
*          On exit, IMB  specifies the true global number of matrix rows
*          of the first block, if M rows are given out starting from the
*          global index IA. If M is equal zero, IMB is set to zero.
*
*  INB     (global output) INTEGER
*          On exit, INB  specifies the true global number of matrix col-
*          umns of the first block, if N columns are given out  starting
*          from the global index JA. If N is equal zero, INB  is  set to
*          zero.
*
*  MB      (global output) INTEGER
*          On exit,  MB specifies the size of the blocks used to  parti-
*          tion the matrix rows. MB is at least one.
*
*  NB      (global output) INTEGER
*          On exit,  NB specifies the size of the blocks used to  parti-
*          tion the matrix columns. NB is at least one.
*
*  PROW    (global output) INTEGER
*          On exit,  PROW  specifies  the  row coordinate of the process
*          that possesses the first row of the submatrix.  On exit, PROW
*          is -1 if DESC( RSRC_ )  is -1 on input, and,  at  least  zero
*          and strictly less than NPROW otherwise.
*
*  PCOL    (global output) INTEGER
*          On exit, PCOL  specifies the column coordinate of the process
*          that possesses the first column of the  submatrix.  On  exit,
*          PCOL is -1 if DESC( CSRC_ )  is -1 on input, and,  at   least
*          zero and strictly less than NPCOL otherwise.
*
*  DA0     (global and local output) INTEGER array
*          On exit, DA0 is an integer array of dimension DLEN_.  This is
*          the array descriptor of the submatrix A(IA:IA+M-1,JA:JA+N-1).
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/* ..
*  .. Executable Statements ..
*
*/
   *MB = DA[MB_]; Mfirstnb( *IMB, M, IA, DA[IMB_], *MB );
   *NB = DA[NB_]; Mfirstnb( *INB, N, JA, DA[INB_], *NB ); *LDA = DA[LLD_];
   Minfog2l( IA, JA, DA, NPROW, NPCOL, MYROW, MYCOL, *II, *JJ, *PROW, *PCOL );
   MDescSet( DA0, M, N, *IMB, *INB, *MB, *NB, *PROW, *PCOL, DA[CTXT_], *LDA );
/*
*  End of PB_Cdescribe
*/
}
