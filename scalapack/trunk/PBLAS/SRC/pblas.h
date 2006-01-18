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
*  This file includes the standard C libraries, as well as system depen-
*  dent include files. All PBLAS routines include this file.
*
*  ---------------------------------------------------------------------
*  Machine Specific PBLAS macros
*  ---------------------------------------------------------------------
*/
#define    _HAL_               0
#define    _T3D_               1
#define    _T3E_               2

#ifdef T3D
#define    _MACH_              _T3D_
#endif
#ifdef T3E
#define    _MACH_              _T3E_
#endif
#ifndef _MACH_
#define    _MACH_              _HAL_
#endif
/*
*  CBRATIO is the ratio of the transfer cost per element for the combine
*  sum to one process and the broadcast operation.  This  value  is used
*  within the Level 3 PBLAS routines to decide on which algorithm to se-
*  lect.
*/
#define    CBRATIO             1.3
/*
*  ---------------------------------------------------------------------
*  Include files
*  ---------------------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>

#ifdef __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif

#if( ( _MACH_ == _T3D_ ) || ( _MACH_ == _T3E_ ) )
#include <fortran.h>
#endif
/*
*  ---------------------------------------------------------------------
*  FORTRAN <-> C interface
*  ---------------------------------------------------------------------
*
*  These macros identifies how the PBLAS will be called as follows:
*
*  _F2C_ADD_: the FORTRAN compiler expects the name of C functions to be
*  in all lower case and to have an underscore postfixed it (Suns, Intel
*  compilers expect this).
*
*  _F2C_NOCHANGE: the FORTRAN compiler expects the name of  C  functions
*  to be in all lower case (IBM RS6K compilers do this).
*
*  _F2C_UPCASE: the  FORTRAN  compiler expects the name of  C  functions
*  to be in all upcase. (Cray compilers expect this).
*
*  _F2C_F77ISF2C: the  FORTRAN  compiler in use is f2c, a  FORTRAN  to C
*  converter.
*/
#define    _F2C_ADD_           0
#define    _F2C_NOCHANGE       1
#define    _F2C_UPCASE         2
#define    _F2C_F77ISF2C       3

#ifdef UpCase
#define    _F2C_CALL_          _F2C_UPCASE
#endif

#ifdef NoChange
#define    _F2C_CALL_          _F2C_NOCHANGE
#endif

#ifdef Add_
#define    _F2C_CALL_          _F2C_ADD_
#endif

#ifdef f77IsF2C
#define    _F2C_CALL_          _F2C_F77ISF2C
#endif

#ifndef _F2C_CALL_
#define    _F2C_CALL_          _F2C_ADD_
#endif
/*
*  ---------------------------------------------------------------------
*  TYPE DEFINITIONS AND CONVERSION UTILITIES
*  ---------------------------------------------------------------------
*/
#if( ( _MACH_ == _T3D_ ) || ( _MACH_ == _T3E_ ) )

#define    float               double
                      /* Type of character argument in a FORTRAN call */
#define    F_CHAR_T            _fcd
                                    /* Character conversion utilities */
#define    F2C_CHAR(a)         ( _fcdtocp( (a) ) )
#define    C2F_CHAR(a)         ( _cptofcd( (a), 1 ) )
                                         /* Type of FORTRAN functions */
#define    F_VOID_FCT          void   fortran           /* Subroutine */
#define    F_INTG_FCT          int    fortran     /* INTEGER function */

#else                 /* Type of character argument in a FORTRAN call */

typedef    char *              F_CHAR_T;
                                    /* Character conversion utilities */
#define    F2C_CHAR(a)            (a)
#define    C2F_CHAR(a)            (a)
                                         /* Type of FORTRAN functions */
#define    F_VOID_FCT             void                  /* Subroutine */
#define    F_INTG_FCT             int             /* INTEGER function */

#endif
/*
* ----------------------------------------------------------------------
*  #typedef definitions
*  ---------------------------------------------------------------------
*/
typedef    float               cmplx  [2];
typedef    double              cmplx16[2];

#define    REAL_PART           0
#define    IMAG_PART           1

#ifdef __STDC__

typedef void           (*GESD2D_T)   ( int,       int,       int,
                                       char *,    int,       int,
                                       int );
typedef void           (*GERV2D_T)   ( int,       int,       int,
                                       char *,    int,       int,
                                       int );
typedef void           (*GEBS2D_T)   ( int,       char *,    char *,
                                       int,       int,       char *,
                                       int );
typedef void           (*GEBR2D_T)   ( int,       char *,    char *,
                                       int,       int,       char *,
                                       int,       int,       int );
typedef void           (*GSUM2D_T)   ( int,       char *,    char *,
                                       int,       int,       char *,
                                       int,       int,       int );

typedef F_VOID_FCT     (*MMADD_T)    ( int  *,    int  *,    char *,
                                       char *,    int  *,    char *,
                                       char *,    int  * );
typedef F_VOID_FCT     (*MMSHFT_T)   ( int  *,    int  *,    int *,
                                       char *,    int  * );
typedef F_VOID_FCT     (*VVDOT_T)    ( int  *,    char *,    char *,
                                       int  *,    char *,    int  * );
typedef F_VOID_FCT     (*VVSET_T)    ( int  *,    char *,    char *,
                                       int  * );
typedef F_VOID_FCT     (*TZPAD_T)    ( F_CHAR_T,  F_CHAR_T,  int  *,
                                       int  *,    int  *,    char *,
                                       char *,    char *,    int  * );
typedef F_VOID_FCT     (*TZPADCPY_T) ( F_CHAR_T,  F_CHAR_T,  int  *,
                                       int  *,    int  *,    char *,
                                       int *,     char *,    int  * );
typedef F_VOID_FCT     (*TZSET_T)    ( F_CHAR_T,  int  *,    int  *,
                                       int  *,    char *,    char *,
                                       char *,    int  * );
typedef F_VOID_FCT     (*TZSCAL_T)   ( F_CHAR_T,  int *,     int  *,
                                       int  *,    char *,    char *,
                                       int  * );

typedef F_VOID_FCT     (*AXPY_T)     ( int *,     char *,    char *,
                                       int *,     char *,    int * );
typedef F_VOID_FCT     (*COPY_T)     ( int *,     char *,    int *,
                                       char *,    int * );
typedef F_VOID_FCT     (*SWAP_T)     ( int *,     char *,    int *,
                                       char *,    int * );

typedef F_VOID_FCT     (*GEMV_T)     ( F_CHAR_T,  int *,     int *,
                                       char *,    char *,    int *,
                                       char *,    int *,     char *,
                                       char *,    int * );
typedef F_VOID_FCT     (*AGEMV_T)    ( F_CHAR_T,  int *,     int *,
                                       char *,    char *,    int *,
                                       char *,    int *,     char *,
                                       char *,    int * );
typedef F_VOID_FCT     (*SYMV_T)     ( F_CHAR_T,  int *,     char *,
                                       char *,    int *,     char *,
                                       int *,     char *,    char *,
                                       int * );
typedef F_VOID_FCT     (*ASYMV_T)    ( F_CHAR_T,  int *,     char *,
                                       char *,    int *,     char *,
                                       int *,     char *,    char *,
                                       int * );
typedef F_VOID_FCT     (*HEMV_T)     ( F_CHAR_T,  int *,     char *,
                                       char *,    int *,     char *,
                                       int *,     char *,    char *,
                                       int * );
typedef F_VOID_FCT     (*AHEMV_T)    ( F_CHAR_T,  int *,     char *,
                                       char *,    int *,     char *,
                                       int *,     char *,    char *,
                                       int * );
typedef F_VOID_FCT     (*TRMV_T)     ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                       int *,     char *,    int *,
                                       char *,    int * );
typedef F_VOID_FCT     (*ATRMV_T)    ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                       int *,     char *,    char *,
                                       int *,     char *,    int *,
                                       char *,    char *,    int * );
typedef F_VOID_FCT     (*TRSV_T)     ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                       int *,     char *,    int *,
                                       char *,    int * );
typedef F_VOID_FCT     (*GERC_T)     ( int *,     int *,     char *,
                                       char *,    int *,     char *,
                                       int *,     char *,    int * );
typedef F_VOID_FCT     (*GERU_T)     ( int *,     int *,     char *,
                                       char *,    int *,     char *,
                                       int *,     char *,    int * );
typedef F_VOID_FCT     (*SYR_T)      ( F_CHAR_T,  int *,     char *,
                                       char *,    int *,     char *,
                                       int * );
typedef F_VOID_FCT     (*HER_T)      ( F_CHAR_T,  int *,     char *,
                                       char *,    int *,     char *,
                                       int * );
typedef F_VOID_FCT     (*SYR2_T)     ( F_CHAR_T,  int *,     char *,
                                       char *,    int *,     char *,
                                       int *,     char *,    int * );
typedef F_VOID_FCT     (*HER2_T)     ( F_CHAR_T,  int *,     char *,
                                       char *,    int *,     char *,
                                       int *,     char *,    int * );

typedef F_VOID_FCT     (*GEMM_T)     ( F_CHAR_T,  F_CHAR_T,  int *,
                                       int *,     int *,     char *,
                                       char *,    int *,     char *,
                                       int *,     char *,    char *,
                                       int * );
typedef F_VOID_FCT     (*SYMM_T)     ( F_CHAR_T,  F_CHAR_T,  int *,
                                       int *,     char *,    char *,
                                       int *,     char *,    int *,
                                       char *,    char *,    int * );
typedef F_VOID_FCT     (*HEMM_T)     ( F_CHAR_T,  F_CHAR_T,  int *,
                                       int *,     char *,    char *,
                                       int *,     char *,    int *,
                                       char *,    char *,    int * );
typedef F_VOID_FCT     (*SYRK_T)     ( F_CHAR_T,  F_CHAR_T,  int *,
                                       int *,     char *,    char *,
                                       int *,     char *,    char *,
                                       int * );
typedef F_VOID_FCT     (*HERK_T)     ( F_CHAR_T,  F_CHAR_T,  int *,
                                       int *,     char *,    char *,
                                       int *,     char *,    char *,
                                       int * );
typedef F_VOID_FCT     (*SYR2K_T)    ( F_CHAR_T,  F_CHAR_T,  int *,
                                       int *,     char *,    char *,
                                       int *,     char *,    int *,
                                       char *,    char *,    int * );
typedef F_VOID_FCT     (*HER2K_T)    ( F_CHAR_T,  F_CHAR_T,  int *,
                                       int *,     char *,    char *,
                                       int *,     char *,    int *,
                                       char *,    char *,    int * );
typedef F_VOID_FCT     (*TRMM_T)     ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                       F_CHAR_T,  int *,     int *,
                                       char *,    char *,    int *,
                                       char *,    int * );
typedef F_VOID_FCT     (*TRSM_T)     ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                       F_CHAR_T,  int *,     int *,
                                       char *,    char *,    int *,
                                       char *,    int * );

#else

typedef void           (*GESD2D_T)   ();
typedef void           (*GERV2D_T)   ();
typedef void           (*GEBS2D_T)   ();
typedef void           (*GEBR2D_T)   ();
typedef void           (*GSUM2D_T)   ();

typedef F_VOID_FCT     (*MMADD_T)    ();
typedef F_VOID_FCT     (*MMSHFT_T)   ();
typedef F_VOID_FCT     (*VVDOT_T)    ();
typedef F_VOID_FCT     (*VVSET_T)    ();
typedef F_VOID_FCT     (*TZPAD_T)    ();
typedef F_VOID_FCT     (*TZPADCPY_T) ();
typedef F_VOID_FCT     (*TZSET_T)    ();
typedef F_VOID_FCT     (*TZSCAL_T)   ();

typedef F_VOID_FCT     (*AXPY_T)     ();
typedef F_VOID_FCT     (*COPY_T)     ();
typedef F_VOID_FCT     (*SWAP_T)     ();

typedef F_VOID_FCT     (*GEMV_T)     ();
typedef F_VOID_FCT     (*AGEMV_T)    ();
typedef F_VOID_FCT     (*SYMV_T)     ();
typedef F_VOID_FCT     (*ASYMV_T)    ();
typedef F_VOID_FCT     (*HEMV_T)     ();
typedef F_VOID_FCT     (*AHEMV_T)    ();
typedef F_VOID_FCT     (*TRMV_T)     ();
typedef F_VOID_FCT     (*ATRMV_T)    ();
typedef F_VOID_FCT     (*TRSV_T)     ();
typedef F_VOID_FCT     (*GERC_T)     ();
typedef F_VOID_FCT     (*GERU_T)     ();
typedef F_VOID_FCT     (*SYR_T)      ();
typedef F_VOID_FCT     (*HER_T)      ();
typedef F_VOID_FCT     (*SYR2_T)     ();
typedef F_VOID_FCT     (*HER2_T)     ();

typedef F_VOID_FCT     (*GEMM_T)     ();
typedef F_VOID_FCT     (*SYMM_T)     ();
typedef F_VOID_FCT     (*HEMM_T)     ();
typedef F_VOID_FCT     (*SYRK_T)     ();
typedef F_VOID_FCT     (*HERK_T)     ();
typedef F_VOID_FCT     (*SYR2K_T)    ();
typedef F_VOID_FCT     (*HER2K_T)    ();
typedef F_VOID_FCT     (*TRMM_T)     ();
typedef F_VOID_FCT     (*TRSM_T)     ();

#endif

typedef struct
{
   char           type;                  /* Encoding of the data type */
   int            usiz;    /* length in bytes of elementary data type */
   int            size;               /* length in bytes of data type */

   char           * zero,
                  * one,
                  * negone;   /* pointers to contants of correct type */

   GESD2D_T       Cgesd2d;                         /* BLACS functions */
   GERV2D_T       Cgerv2d;
   GEBS2D_T       Cgebs2d;
   GEBR2D_T       Cgebr2d;
   GSUM2D_T       Cgsum2d;

   MMADD_T        Fmmadd;                       /* Addition functions */
   MMADD_T        Fmmcadd;
   MMADD_T        Fmmtadd;
   MMADD_T        Fmmtcadd;
   MMADD_T        Fmmdda;
   MMADD_T        Fmmddac;
   MMADD_T        Fmmddat;
   MMADD_T        Fmmddact;

   MMSHFT_T       Fcshft;                          /* Shift functions */
   MMSHFT_T       Frshft;

   VVDOT_T        Fvvdotu;                           /* Dot functions */
   VVDOT_T        Fvvdotc;

   TZPAD_T        Ftzpad;                       /* Array pad function */
   TZPADCPY_T     Ftzpadcpy;
   VVSET_T        Fset;

   TZSCAL_T       Ftzscal;                       /* Scaling functions */
   TZSCAL_T       Fhescal;
   TZSCAL_T       Ftzcnjg;

   AXPY_T         Faxpy;                              /* Level 1 BLAS */
   COPY_T         Fcopy;
   SWAP_T         Fswap;

   GEMV_T         Fgemv;                              /* Level 2 BLAS */
   SYMV_T         Fsymv;
   HEMV_T         Fhemv;
   TRMV_T         Ftrmv;
   TRSV_T         Ftrsv;

   AGEMV_T        Fagemv;
   ASYMV_T        Fasymv;
   AHEMV_T        Fahemv;
   ATRMV_T        Fatrmv;

   GERC_T         Fgerc;
   GERU_T         Fgeru;
   SYR_T          Fsyr;
   HER_T          Fher;
   SYR2_T         Fsyr2;
   HER2_T         Fher2;

   GEMM_T         Fgemm;                              /* Level 3 BLAS */
   SYMM_T         Fsymm;
   HEMM_T         Fhemm;
   SYRK_T         Fsyrk;
   HERK_T         Fherk;
   SYR2K_T        Fsyr2k;
   HER2K_T        Fher2k;
   TRMM_T         Ftrmm;
   TRSM_T         Ftrsm;

} PBTYP_T;

#ifdef __STDC__

typedef void           (*TZSYR_T)    ( PBTYP_T *, char *,    int,
                                       int,       int,       int,
                                       char *,    char *,    int,
                                       char *,    int,       char *,
                                       int );
typedef void           (*TZSYR2_T)   ( PBTYP_T *, char *,    int,
                                       int,       int,       int,
                                       char *,    char *,    int,
                                       char *,    int,       char *,
                                       int,       char *,    int,
                                       char *,    int );
typedef void           (*TZTRM_T)    ( PBTYP_T *, char *,    char *,
                                       char *,    char *,    int,
                                       int,       int,       int,
                                       char *,    char *,    int,
                                       char *,    int,       char *,
                                       int );
typedef void           (*TZSYM_T)    ( PBTYP_T *, char *,    char *,
                                       int,       int,       int,
                                       int,       char *,    char *,
                                       int,       char *,    int,
                                       char *,    int,       char *,
                                       int,       char *,    int );
#else

typedef void           (*TZSYR_T)    ();
typedef void           (*TZSYR2_T)   ();
typedef void           (*TZTRM_T)    ();
typedef void           (*TZSYM_T)    ();

#endif

typedef struct
{
   int offd;                                /* Global diagonal offset */
   int lcmt00;                            /* LCM value of first block */

   int mp;                                    /* Local number of rows */
   int imb1;                      /* Size of first row block (global) */
   int imbloc;                       /* Size of first local row block */
   int mb;                                          /* Row block size */
   int lmbloc;                        /* Size of last local row block */
   int mblks;                           /* Number of local row blocks */
   int iupp;                /* LCM row bound for first diagonal block */
   int upp;                       /* LCM row bound for diagonal block */
   int prow;                       /* Relative row process coordinate */
   int nprow;                               /* Number of process rows */

   int nq;                                 /* Local number of columns */
   int inb1;                   /* Size of first column block (global) */
   int inbloc;                    /* Size of first local column block */
   int nb;                                       /* Column block size */
   int lnbloc;                     /* Size of last local column block */
   int nblks;                        /* Number of local column blocks */
   int ilow;             /* LCM column bound for first diagonal block */
   int low;                    /* LCM column bound for diagonal block */
   int pcol;                    /* Relative column process coordinate */
   int npcol;                            /* Number of process columns */

   int lcmb;    /* Least common multiple of nprow * mb and npcol * nb */

} PB_VM_T;

/*
*  ---------------------------------------------------------------------
*  #define macro constants
*  ---------------------------------------------------------------------
*/
#define    INT                 'I'                /* type identifiers */
#define    SREAL               'S'
#define    DREAL               'D'
#define    SCPLX               'C'
#define    DCPLX               'Z'

#define crot_ CROT
