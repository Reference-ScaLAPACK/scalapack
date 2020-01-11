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
void PB_Cdescset( Int * DESC, Int M, Int N, Int IMB, Int INB, Int MB, Int NB,
                  Int RSRC, Int CSRC, Int CTXT, Int LLD )
#else
void PB_Cdescset( DESC, M, N, IMB, INB, MB, NB, RSRC, CSRC, CTXT, LLD )
/*
*  .. Scalar Arguments ..
*/
   Int            CSRC, CTXT, IMB, INB, LLD, M, MB, N, NB, RSRC;
/*
*  .. Array Arguments ..
*/
   Int            * DESC;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cdescset uses  its  10  input  arguments  M,  N, IMB, INB, MB, NB,
*  RSRC,  CSRC,  CTXT  and LLD to initialize a descriptor vector of type
*  BLOCK_CYCLIC_2D_INB.
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
*  DESC    (global and local output) INTEGER array
*          On entry, DESC is an array of  dimension  DLEN_.  DESC is the
*          array descriptor to be set.
*
*  M       (global input) INTEGER
*          On entry,  M  specifies  the  number  of  rows of the matrix.
*          M must be at least zero.
*
*  N       (global input) INTEGER
*          On entry,  N  specifies  the number of columns of the matrix.
*          N must be at least zero.
*
*  IMB     (global input) INTEGER
*          On entry,  IMB  specifies  the row size of the first block of
*          the global matrix distribution. IMB must be at least one.
*
*  INB     (global input) INTEGER
*          On entry,  INB  specifies  the column size of the first block
*          of the global matrix distribution. INB must be at least one.
*
*  MB      (global input) INTEGER
*          On entry,  MB  specifies  the  row size of the blocks used to
*          partition the matrix. MB must be at least one.
*
*  NB      (global input) INTEGER
*          On entry, NB  specifies the column size of the blocks used to
*          partition the matrix. NB must be at least one.
*
*  RSRC    (global input) INTEGER
*          On entry,  RSRC  specifies  the row coordinate of the process
*          that possesses the first row of the matrix.  When  RSRC = -1,
*          the data is not  distributed but replicated,  otherwise  RSRC
*          must be at least zero and strictly less than NPROW.
*
*  CSRC    (global input) INTEGER
*          On entry,  CSRC  specifies  the column coordinate of the pro-
*          cess  that  possesses  the  first column of the matrix.  When
*          CSRC = -1, the data is not distributed but replicated, other-
*          wise CSRC must be at least zero and strictly less than NPCOL.
*
*  CTXT    (local input) INTEGER
*          On entry, CTXT specifies the BLACS context handle, indicating
*          the global  communication  context.  The value of the context
*          itself is local.
*
*  LLD     (local input)  INTEGER
*          On entry, LLD  specifies  the  leading dimension of the local
*          array storing the local entries of the matrix. LLD must be at
*          least MAX( 1, Lr(1,M) ).
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
   DESC[DTYPE_] = BLOCK_CYCLIC_2D_INB;
   DESC[CTXT_ ] = CTXT;
   DESC[M_    ] = M;
   DESC[N_    ] = N;
   DESC[IMB_  ] = IMB;
   DESC[INB_  ] = INB;
   DESC[MB_   ] = MB;
   DESC[NB_   ] = NB;
   DESC[RSRC_ ] = RSRC;
   DESC[CSRC_ ] = CSRC;
   DESC[LLD_  ] = LLD;
/*
*  End of PB_Cdescset
*/
}
