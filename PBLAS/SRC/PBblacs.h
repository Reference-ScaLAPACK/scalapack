/* ---------------------------------------------------------------------
*
*  -- PBLAS routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*  ---------------------------------------------------------------------
*/
/*
*  This file includes  BLACS  function type definitions,  define macros,
*  and function prototypes. All PBLAS routines include this file.
*
* ----------------------------------------------------------------------
*  #define macro constants
*  ---------------------------------------------------------------------
*/
                                       /* BLACS scopes and topologies */
/* #define    CALL                'A'               (already defined) */
#define    CCOLUMN             'C'
#define    CROW                'R'

#define    CBCAST              'B'
#define    CCOMBINE            'C'
#define    CTOP_GET            '!'
#define    CTOP_DEFAULT        ' '
#define    CTOP_IRING          'I'
#define    CTOP_DRING          'D'
#define    CTOP_SRING          'S'
#define    CTOP_HYPER          'H'
#define    CTOP_FULL           'F'
#define    CTOP_MRING          'M'
#define    CTOP_TTREE          'T'
#define    CTOP_TREE1          '1'
#define    CTOP_TREE2          '2'
#define    CTOP_TREE3          '3'
#define    CTOP_TREE4          '4'
#define    CTOP_TREE5          '5'
#define    CTOP_TREE6          '6'
#define    CTOP_TREE7          '7'
#define    CTOP_TREE8          '8'
#define    CTOP_TREE9          '9'

/* #define    ALL                 "A"               (already defined) */
#define    COLUMN              "C"
#define    ROW                 "R"

#define    BCAST               "B"
#define    COMBINE             "C"
#define    TOP_GET             "!"
#define    TOP_DEFAULT         " "
#define    TOP_IRING           "I"
#define    TOP_DRING           "D"
#define    TOP_SRING           "S"
#define    TOP_HYPER           "H"
#define    TOP_FULL            "F"
#define    TOP_MRING           "M"
#define    TOP_TTREE           "T"
#define    TOP_TREE1           "1"
#define    TOP_TREE2           "2"
#define    TOP_TREE3           "3"
#define    TOP_TREE4           "4"
#define    TOP_TREE5           "5"
#define    TOP_TREE6           "6"
#define    TOP_TREE7           "7"
#define    TOP_TREE8           "8"
#define    TOP_TREE9           "9"

/*
*  ---------------------------------------------------------------------
*  Function prototypes
*  ---------------------------------------------------------------------
*/
#ifdef __STDC__
                                              /* BLACS Initialization */
void           Cblacs_pinfo    ( Int *,     Int * );
void           Cblacs_setup    ( Int *,     Int * );
void           Cblacs_get      ( Int,       Int,       Int * );
void           Cblacs_set      ( Int,       Int,       Int * );
void           Cblacs_gridinit ( Int *,     char *,    Int,
                                 Int );
void           Cblacs_gridmap  ( Int *,     Int *,     Int,
                                 Int,       Int );

                                                 /* BLACS Destruction */
void           Cblacs_freebuff ( Int,       Int );
void           Cblacs_gridexit ( Int );
void           Cblacs_abort    ( Int,       Int );
void           Cblacs_exit     ( Int );

                             /* BLACS Informational and Miscellaneous */
void           Cblacs_gridinfo ( Int,       Int *,     Int *,
                                 Int *,     Int * );
Int            Cblacs_pnum     ( Int,       Int,       Int );
void           Cblacs_pcoord   ( Int,       Int,       Int *,
                                 Int * );
void           Cblacs_barrier  ( Int,       char * );

                                                     /* BLACS Sending */
void           Cigesd2d        ( Int,       Int,       Int,
                                 char *,    Int,       Int,
                                 Int );
void           Csgesd2d        ( Int,       Int,       Int,
                                 char *,    Int,       Int,
                                 Int );
void           Cdgesd2d        ( Int,       Int,       Int,
                                 char *,    Int,       Int,
                                 Int );
void           Ccgesd2d        ( Int,       Int,       Int,
                                 char *,    Int,       Int,
                                 Int );
void           Czgesd2d        ( Int,       Int,       Int,
                                 char *,    Int,       Int,
                                 Int );

void           Citrsd2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Cstrsd2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Cdtrsd2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Cctrsd2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Cztrsd2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );

void           Cigebs2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int );
void           Csgebs2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int );
void           Cdgebs2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int );
void           Ccgebs2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int );
void           Czgebs2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int );

void           Citrbs2d        ( Int,       char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    Int );
void           Cstrbs2d        ( Int,       char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    Int );
void           Cdtrbs2d        ( Int,       char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    Int );
void           Cctrbs2d        ( Int,       char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    Int );
void           Cztrbs2d        ( Int,       char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    Int );

                                                   /* BLACS Receiving */
void           Cigerv2d        ( Int,       Int,       Int,
                                 char *,    Int,       Int,
                                 Int );
void           Csgerv2d        ( Int,       Int,       Int,
                                 char *,    Int,       Int,
                                 Int );
void           Cdgerv2d        ( Int,       Int,       Int,
                                 char *,    Int,       Int,
                                 Int );
void           Ccgerv2d        ( Int,       Int,       Int,
                                 char *,    Int,       Int,
                                 Int );
void           Czgerv2d        ( Int,       Int,       Int,
                                 char *,    Int,       Int,
                                 Int );

void           Citrrv2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Cstrrv2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Cdtrrv2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Cctrrv2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Cztrrv2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );

void           Cigebr2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Csgebr2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Cdgebr2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Ccgebr2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Czgebr2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );

void           Citrbr2d        ( Int,       char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    Int,
                                 Int,       Int );
void           Cstrbr2d        ( Int,       char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    Int,
                                 Int,       Int );
void           Cdtrbr2d        ( Int,       char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    Int,
                                 Int,       Int );
void           Cctrbr2d        ( Int,       char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    Int,
                                 Int,       Int );
void           Cztrbr2d        ( Int,       char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    Int,
                                 Int,       Int );

                                          /* BLACS Combine Operations */
void           Cigamx2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int *,     Int *,
                                 Int,       Int,       Int );
void           Csgamx2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int *,     Int *,
                                 Int,       Int,       Int );
void           Cdgamx2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int *,     Int *,
                                 Int,       Int,       Int );
void           Ccgamx2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int *,     Int *,
                                 Int,       Int,       Int );
void           Czgamx2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int *,     Int *,
                                 Int,       Int,       Int );

void           Cigamn2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int *,     Int *,
                                 Int,       Int,       Int );
void           Csgamn2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int *,     Int *,
                                 Int,       Int,       Int );
void           Cdgamn2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int *,     Int *,
                                 Int,       Int,       Int );
void           Ccgamn2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int *,     Int *,
                                 Int,       Int,       Int );
void           Czgamn2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int *,     Int *,
                                 Int,       Int,       Int );

void           Cigsum2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Csgsum2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Cdgsum2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Ccgsum2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );
void           Czgsum2d        ( Int,       char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int );

#else
                                              /* BLACS Initialization */
void           Cblacs_pinfo    ();
void           Cblacs_setup    ();
void           Cblacs_get      ();
void           Cblacs_set      ();
void           Cblacs_gridinit ();
void           Cblacs_gridmap  ();

                                                 /* BLACS Destruction */
void           Cblacs_freebuff ();
void           Cblacs_gridexit ();
void           Cblacs_abort    ();
void           Cblacs_exit     ();

                             /* BLACS Informational and Miscellaneous */
void           Cblacs_gridinfo ();
Int            Cblacs_pnum     ();
void           Cblacs_pcoord   ();
void           Cblacs_barrier  ();

                                                     /* BLACS Sending */
void           Cigesd2d        ();
void           Csgesd2d        ();
void           Cdgesd2d        ();
void           Ccgesd2d        ();
void           Czgesd2d        ();

void           Citrsd2d        ();
void           Cstrsd2d        ();
void           Cdtrsd2d        ();
void           Cctrsd2d        ();
void           Cztrsd2d        ();

void           Cigebs2d        ();
void           Csgebs2d        ();
void           Cdgebs2d        ();
void           Ccgebs2d        ();
void           Czgebs2d        ();

void           Citrbs2d        ();
void           Cstrbs2d        ();
void           Cdtrbs2d        ();
void           Cctrbs2d        ();
void           Cztrbs2d        ();

                                                   /* BLACS Receiving */
void           Cigerv2d        ();
void           Csgerv2d        ();
void           Cdgerv2d        ();
void           Ccgerv2d        ();
void           Czgerv2d        ();

void           Citrrv2d        ();
void           Cstrrv2d        ();
void           Cdtrrv2d        ();
void           Cctrrv2d        ();
void           Cztrrv2d        ();

void           Cigebr2d        ();
void           Csgebr2d        ();
void           Cdgebr2d        ();
void           Ccgebr2d        ();
void           Czgebr2d        ();

void           Citrbr2d        ();
void           Cstrbr2d        ();
void           Cdtrbr2d        ();
void           Cctrbr2d        ();
void           Cztrbr2d        ();

                                          /* BLACS Combine Operations */
void           Cigamx2d        ();
void           Csgamx2d        ();
void           Cdgamx2d        ();
void           Ccgamx2d        ();
void           Czgamx2d        ();

void           Cigamn2d        ();
void           Csgamn2d        ();
void           Cdgamn2d        ();
void           Ccgamn2d        ();
void           Czgamn2d        ();

void           Cigsum2d        ();
void           Csgsum2d        ();
void           Cdgsum2d        ();
void           Ccgsum2d        ();
void           Czgsum2d        ();

#endif
