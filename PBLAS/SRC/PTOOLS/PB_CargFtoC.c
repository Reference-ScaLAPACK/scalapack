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
void PB_CargFtoC( Int IF, Int JF, Int * DESCIN, Int * IC, Int * JC,
                    Int * DESCOUT )
#else
void PB_CargFtoC( IF, JF, DESCIN, IC, JC, DESCOUT )
/*
*  .. Scalar Arguments ..
*/
   Int            IF, JF, * IC, * JC;
/*
*  .. Array Arguments ..
*/
   Int            * DESCIN, * DESCOUT;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_CargFtoC converts  a  descriptor  DESCIN  of  type BLOCK_CYCLIC_2D
*  or   BLOCK_CYCLIC_INB_2D   into   a   descriptor   DESCOUT   of  type
*  BLOCK_CYCLIC_INB_2D.
*
*  Notes
*  =====
*
*  A description  vector  is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This  vector  stores  the  information required  to
*  establish the  mapping between a matrix entry and  its  corresponding
*  process and memory location.
*
*  In  the  following  comments,  the  character _  should  be  read  as
*  "of the distributed  matrix".  Let  A  be a generic term for  any  2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION         STORED IN        EXPLANATION
*  ---------------- ---------------  -----------------------------------
*  DTYPE_A (global) DESCA( DTYPE1_ ) The descriptor type.
*  CTXT_A  (global) DESCA( CTXT1_  ) The BLACS context handle indicating
*                                    the   NPROW x NPCOL  BLACS  process
*                                    grid  A  is  distributed  over. The
*                                    context  itself  is global, but the
*                                    handle   (the  integer  value)  may
*                                    vary.
*  M_A     (global) DESCA( M1_     ) The  number  of rows in the distri-
*                                    buted matrix A, M_A >= 0.
*  N_A     (global) DESCA( N1_     ) The  number  of columns in the dis-
*                                    tributed matrix A, N_A >= 0.
*  MB_A    (global) DESCA( MB1_    ) The blocking factor used to distri-
*                                    bute the rows of A, MB_A > 0.
*  NB_A    (global) DESCA( NB1_    ) The blocking factor used to distri-
*                                    bute the columns of A, NB_A > 0.
*  RSRC_A  (global) DESCA( RSRC1_  ) The  process  row  over  which  the
*                                    first row of the matrix  A  is dis-
*                                    tributed, NPROW > RSRC_A >= 0.
*  CSRC_A  (global) DESCA( CSRC1_  ) The process column  over  which the
*                                    first column of  A  is distributed.
*                                    NPCOL > CSRC_A >= 0.
*  LLD_A   (local)  DESCA( LLD1_   ) The leading dimension  of the local
*                                    array  storing  the local blocks of
*                                    the distributed matrix A,
*                                    IF( Lc( 1, N_A ) > 0 )
*                                      LLD_A >= MAX( 1, Lr( 1, M_A ) )
*                                    ELSE
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
*  Lr( IA, K ) = PB_Cnumroc( K, IA, MB_A, MB_A, MYROW, RSRC_A, NPROW )
*  Lc( JA, K ) = PB_Cnumroc( K, JA, NB_A, NB_A, MYCOL, CSRC_A, NPCOL )
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
*  IF      (global input) INTEGER
*          On entry,  IF  specifies  the global row Fortran index of the
*          distributed subarray described by DESCIN. IF must be at least
*          one.
*
*  JF      (global input) INTEGER
*          On entry,  JF  specifies  the  global column Fortran index of
*          the distributed subarray described by DESCIN. JF  must  be at
*          least one.
*
*  DESCIN  (global and local input) INTEGER array
*          On entry, DESCIN  is an array of dimension DLEN1_ or DLEN_ as
*          specified by its first entry DESCIN( DTYPE_ ).  DESCIN is the
*          source  array  descriptor of type BLOCK_CYCLIC_2D  or of type
*          BLOCK_CYCLIC_2D_INB.
*
*  IC      (global output) INTEGER
*          On exit, IC specifies the global row C index of the distribu-
*          ted subarray described by DESCOUT. IC = IF - 1, i.e IC  is at
*          least zero.
*
*  JC      (global output) INTEGER
*          On entry,  JC  specifies  the  global column Fortran index of
*          the distributed subarray described  by  DESCOUT. JC = JF - 1,
*          i.e JC is at least zero.
*
*  DESCOUT (global and local output) INTEGER array
*          On entry, DESCOUT is an array of dimension DLEN_.  DESCOUT is
*          the target array descriptor of type BLOCK_CYCLIC_2D_INB.
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
   *IC = IF - 1;
   *JC = JF - 1;

   if ( DESCIN[DTYPE_] == BLOCK_CYCLIC_2D )
   {
      DESCOUT[DTYPE_] = BLOCK_CYCLIC_2D_INB;
      DESCOUT[M_    ] = DESCIN[M1_    ];
      DESCOUT[N_    ] = DESCIN[N1_    ];
      DESCOUT[IMB_  ] = DESCIN[MB1_   ];
      DESCOUT[INB_  ] = DESCIN[NB1_   ];
      DESCOUT[MB_   ] = DESCIN[MB1_   ];
      DESCOUT[NB_   ] = DESCIN[NB1_   ];
      DESCOUT[RSRC_ ] = DESCIN[RSRC1_ ];
      DESCOUT[CSRC_ ] = DESCIN[CSRC1_ ];
      DESCOUT[CTXT_ ] = DESCIN[CTXT1_ ];
      DESCOUT[LLD_  ] = DESCIN[LLD1_  ];
   }
   else if ( DESCIN[DTYPE_] == BLOCK_CYCLIC_2D_INB )
   {
      DESCOUT[DTYPE_] = BLOCK_CYCLIC_2D_INB;
      DESCOUT[M_    ] = DESCIN[M_    ];
      DESCOUT[N_    ] = DESCIN[N_    ];
      DESCOUT[IMB_  ] = DESCIN[IMB_  ];
      DESCOUT[INB_  ] = DESCIN[INB_  ];
      DESCOUT[MB_   ] = DESCIN[MB_   ];
      DESCOUT[NB_   ] = DESCIN[NB_   ];
      DESCOUT[RSRC_ ] = DESCIN[RSRC_ ];
      DESCOUT[CSRC_ ] = DESCIN[CSRC_ ];
      DESCOUT[CTXT_ ] = DESCIN[CTXT_ ];
      DESCOUT[LLD_  ] = DESCIN[LLD_  ];
   }
   else
   {
      DESCOUT[DTYPE_] = DESCIN[0];
      DESCOUT[CTXT_ ] = DESCIN[1];
      DESCOUT[M_    ] = 0;
      DESCOUT[N_    ] = 0;
      DESCOUT[IMB_  ] = 1;
      DESCOUT[INB_  ] = 1;
      DESCOUT[MB_   ] = 1;
      DESCOUT[NB_   ] = 1;
      DESCOUT[RSRC_ ] = 0;
      DESCOUT[CSRC_ ] = 0;
      DESCOUT[LLD_  ] = 1;
   }
/*
*  End of PB_CargFtoC
*/
}
