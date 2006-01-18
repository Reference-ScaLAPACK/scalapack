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
void           Cblacs_pinfo    ( int *,     int * );
void           Cblacs_setup    ( int *,     int * );
void           Cblacs_get      ( int,       int,       int * );
void           Cblacs_set      ( int,       int,       int * );
void           Cblacs_gridinit ( int *,     char *,    int,
                                 int );
void           Cblacs_gridmap  ( int *,     int *,     int,
                                 int,       int );

                                                 /* BLACS Destruction */
void           Cblacs_freebuff ( int,       int );
void           Cblacs_gridexit ( int );
void           Cblacs_abort    ( int,       int );
void           Cblacs_exit     ( int );

                             /* BLACS Informational and Miscellaneous */
void           Cblacs_gridinfo ( int,       int *,     int *,
                                 int *,     int * );
int            Cblacs_pnum     ( int,       int,       int );
void           Cblacs_pcoord   ( int,       int,       int *,
                                 int * );
void           Cblacs_barrier  ( int,       char * );

                                                     /* BLACS Sending */
void           Cigesd2d        ( int,       int,       int,
                                 char *,    int,       int,
                                 int );
void           Csgesd2d        ( int,       int,       int,
                                 char *,    int,       int,
                                 int );
void           Cdgesd2d        ( int,       int,       int,
                                 char *,    int,       int,
                                 int );
void           Ccgesd2d        ( int,       int,       int,
                                 char *,    int,       int,
                                 int );
void           Czgesd2d        ( int,       int,       int,
                                 char *,    int,       int,
                                 int );

void           Citrsd2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Cstrsd2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Cdtrsd2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Cctrsd2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Cztrsd2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );

void           Cigebs2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int );
void           Csgebs2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int );
void           Cdgebs2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int );
void           Ccgebs2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int );
void           Czgebs2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int );

void           Citrbs2d        ( int,       char *,    char *,
                                 char *,    char *,    int,
                                 int,       char *,    int );
void           Cstrbs2d        ( int,       char *,    char *,
                                 char *,    char *,    int,
                                 int,       char *,    int );
void           Cdtrbs2d        ( int,       char *,    char *,
                                 char *,    char *,    int,
                                 int,       char *,    int );
void           Cctrbs2d        ( int,       char *,    char *,
                                 char *,    char *,    int,
                                 int,       char *,    int );
void           Cztrbs2d        ( int,       char *,    char *,
                                 char *,    char *,    int,
                                 int,       char *,    int );

                                                   /* BLACS Receiving */
void           Cigerv2d        ( int,       int,       int,
                                 char *,    int,       int,
                                 int );
void           Csgerv2d        ( int,       int,       int,
                                 char *,    int,       int,
                                 int );
void           Cdgerv2d        ( int,       int,       int,
                                 char *,    int,       int,
                                 int );
void           Ccgerv2d        ( int,       int,       int,
                                 char *,    int,       int,
                                 int );
void           Czgerv2d        ( int,       int,       int,
                                 char *,    int,       int,
                                 int );

void           Citrrv2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Cstrrv2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Cdtrrv2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Cctrrv2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Cztrrv2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );

void           Cigebr2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Csgebr2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Cdgebr2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Ccgebr2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Czgebr2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );

void           Citrbr2d        ( int,       char *,    char *,
                                 char *,    char *,    int,
                                 int,       char *,    int,
                                 int,       int );
void           Cstrbr2d        ( int,       char *,    char *,
                                 char *,    char *,    int,
                                 int,       char *,    int,
                                 int,       int );
void           Cdtrbr2d        ( int,       char *,    char *,
                                 char *,    char *,    int,
                                 int,       char *,    int,
                                 int,       int );
void           Cctrbr2d        ( int,       char *,    char *,
                                 char *,    char *,    int,
                                 int,       char *,    int,
                                 int,       int );
void           Cztrbr2d        ( int,       char *,    char *,
                                 char *,    char *,    int,
                                 int,       char *,    int,
                                 int,       int );

                                          /* BLACS Combine Operations */
void           Cigamx2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int *,     int *,
                                 int,       int,       int );
void           Csgamx2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int *,     int *,
                                 int,       int,       int );
void           Cdgamx2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int *,     int *,
                                 int,       int,       int );
void           Ccgamx2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int *,     int *,
                                 int,       int,       int );
void           Czgamx2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int *,     int *,
                                 int,       int,       int );

void           Cigamn2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int *,     int *,
                                 int,       int,       int );
void           Csgamn2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int *,     int *,
                                 int,       int,       int );
void           Cdgamn2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int *,     int *,
                                 int,       int,       int );
void           Ccgamn2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int *,     int *,
                                 int,       int,       int );
void           Czgamn2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int *,     int *,
                                 int,       int,       int );

void           Cigsum2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Csgsum2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Cdgsum2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Ccgsum2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );
void           Czgsum2d        ( int,       char *,    char *,
                                 int,       int,       char *,
                                 int,       int,       int );

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
int            Cblacs_pnum     ();
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
