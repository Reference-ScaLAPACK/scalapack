/* ---------------------------------------------------------------------
*
*  -- ScaLAPACK routine (version 1.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 17, 1996
*
*  ---------------------------------------------------------------------
*/

/*
* This file includes the standard C libraries, as well as system
* dependent include files.  All PBLAS routines include this file.
*/

/*
* ========================================================================
* Machine Specific PBLAS macros
* ========================================================================
*/

#define _HAL_           0
#define _T3D_           1

#ifdef T3D
#define _MACH_          _T3D_
#endif

#ifndef _MACH_
#define _MACH_          _HAL_
#endif

/*
* ========================================================================
* Include files
* ========================================================================
*/
#include <stdio.h>
#include <stdlib.h>

#if( _MACH_ == _T3D_ )
#include <fortran.h>
#endif

/*
* ========================================================================
* FORTRAN <-> C interface
* ========================================================================
*
* These macros define how the PBLAS will be called. _F2C_ADD_ assumes 
* that they will be called by FORTRAN, which expects C routines to have
* an underscore postfixed to the name (Suns, and Intel machines expect
* this). _F2C_NOCHANGE indicates that FORTRAN will be calling, and that
* it expects the name called by FORTRAN to be identical to that compiled
* by the C (RS6K's do this).  _F2C_UPCASE says it expects C routines
* called by FORTRAN to be in all upcase (CRAY wants this).   
* _F2C_F77ISF2C indicates that the fortran "compiler" in use is
* actually f2c, a FORTRAN to C converter.
*/

#define _F2C_ADD_       0
#define _F2C_NOCHANGE   1
#define _F2C_UPCASE     2
#define _F2C_F77ISF2C   3

#ifdef UpCase
#define _F2C_CALL_      _F2C_UPCASE
#endif

#ifdef NoChange
#define _F2C_CALL_      _F2C_NOCHANGE
#endif

#ifdef Add_
#define _F2C_CALL_      _F2C_ADD_
#endif

#ifdef f77IsF2C
#define _F2C_CALL_      _F2C_F77ISF2C
#endif

#ifndef _F2C_CALL_
#define _F2C_CALL_      _F2C_ADD_
#endif

/*
* ========================================================================
* TYPE DEFINITIONS AND CONVERSION UTILITIES
* ========================================================================
*/

typedef struct { float  re, im; } complex;
typedef struct { double re, im; } complex16;

#if( _MACH_ == _T3D_ )
#define float  double
                       /* Type of character argument in a FORTRAN call */
#define F_CHAR          _fcd
                                     /* Character conversion utilities */
#define F2C_CHAR(a)     ( _fcdtocp( (a) ) )
#define C2F_CHAR(a)     ( _cptofcd( (a), 1 ) )
                                          /* Type of FORTRAN functions */
#define F_VOID_FCT      void   fortran                   /* Subroutine */
#define F_INTG_FCT      int    fortran             /* INTEGER function */
#define F_DBLE_FCT      double fortran    /* DOUBLE PRECISION function */

#else
                       /* Type of character argument in a FORTRAN call */
typedef char *          F_CHAR;
                                     /* Character conversion utilities */
#define F2C_CHAR(a)     (a)
#define C2F_CHAR(a)     (a)
                                          /* Type of FORTRAN functions */
#define F_VOID_FCT      void                             /* Subroutine */
#define F_INTG_FCT      int                        /* INTEGER function */
#define F_DBLE_FCT      double            /* DOUBLE PRECISION function */

#endif

/*
* ========================================================================
* #DEFINE MACRO CONSTANTS
* ========================================================================
*/
#define    DLEN_        9                     /* Length of a descriptor */
#define    DT_          0                     /* Descriptor Type        */
#define    CTXT_        1                              /* BLACS context */
#define    M_           2                      /* Global Number of Rows */
#define    N_           3                   /* Global Number of Columns */
#define    MB_          4                          /* Row Blocking Size */
#define    NB_          5                       /* Column Blocking Size */
#define    RSRC_        6                     /* Starting Processor Row */
#define    CSRC_        7                  /* Starting Processor Column */
#define    LLD_         8                    /* Local Leading Dimension */

/*
 * Descriptor types
 */
#define    BLOCK_CYCLIC_2D                1
#define    BLOCK_CYCLIC_INB_2D            2

#define    BROADCAST    "B"              /* Blacs operation definitions */
#define    COMBINE      "C"

#define    ALL          "A"                        /* Scope definitions */
#define    COLUMN       "C"
#define    ROW          "R"

#define    TOPDEF       " " /* Default BLACS topology, PB-BLAS routines */
#define    CTOPDEF      ' '
#define    TOPGET       "!"

#define    YES          "Y"
#define    NO           "N"

#define    MULLENFAC    2

#define    ONE          1.0
#define    ZERO         0.0

/*
* ========================================================================
* PREPROCESSOR MACRO FUNCTIONS USED FOR OPTIMIZATION & CONVENIENCE
* ========================================================================
*/

#define ABS(a)   (((a) < 0) ? -(a) : (a))

#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#define CEIL(a,b) ( ((a)+(b)-1) / (b) )

#define Mlowcase(C) ( ((C) > 64 && (C) < 91) ? (C) | 32 : (C) )

#define Mupcase(C) ( ((C) > 96 && (C) < 123) ? (C) & 0xDF : (C) )

#define INDXG2L( iglob, nb, iproc, isrcproc, nprocs )\
    ( (nb) * ( ( (iglob)-1) / ( (nb) * (nprocs) ) ) +\
      ( ( (iglob) - 1 ) % (nb) ) + 1 )

#define INDXL2G( iloc, nb, iproc, isrcproc, nprocs )\
    ( (nprocs) * (nb) * ( ( (iloc) - 1 ) / (nb) ) +\
      ( ( (iloc) - 1 ) % (nb) ) +\
      ( ( (nprocs) + (iproc) - (isrcproc) ) % (nprocs) ) * (nb) + 1 )

#define INDXG2P( iglob, nb, iproc, isrcproc, nprocs ) \
    ( ( (isrcproc) + ( (iglob) - 1 ) / (nb) ) % (nprocs) )

#define MYROC0( nblocks, n, nb, nprocs )\
  ( ( (nblocks) % (nprocs) ) ? ( ( (nblocks) / (nprocs) ) * (nb) + (nb) )\
                   : ( ( (nblocks) / (nprocs) )* (nb) + ( (n) % (nb) ) ) )

#if( _F2C_CALL_ == _F2C_ADD_ )
/*
* These defines set up the naming scheme required to have a FORTRAN
* routine call a C routine (which is what the PBLAS are written in).
* No redefinition necessary to have following FORTRAN to C interface:
*           FORTRAN CALL               C DECLARATION
*           call pdgemm(...)           void pdgemm_(...)
*
* This is the default.
*/

#endif

#if( _F2C_CALL_ == _F2C_UPCASE )
/*
* These defines set up the naming scheme required to have a FORTRAN
* routine call a C routine (which is what the PBLAS are written in)
* following FORTRAN to C interface:
*           FORTRAN CALL               C DECLARATION
*           call pdgemm(...)           void PDGEMM(...)
*/
                                                            /* TOOLS */
#define ilcm_             ILCM
#define infog2l_          INFOG2L
#define numroc_           NUMROC
#define pstreecomb_       PSTREECOMB
#define pdtreecomb_       PDTREECOMB
#define pctreecomb_       PCTREECOMB
#define pztreecomb_       PZTREECOMB
#define scombamax_        SCOMBAMAX
#define dcombamax_        DCOMBAMAX
#define ccombamax_        CCOMBAMAX
#define zcombamax_        ZCOMBAMAX
#define scombnrm2_        SCOMBNRM2
#define dcombnrm2_        DCOMBNRM2

#define dlamov_           DLAMOV
#define slamov_           SLAMOV
#define clamov_           CLAMOV
#define zlamov_           ZLAMOV
#define dlacpy_           DLACPY
#define slacpy_           SLACPY
#define clacpy_           CLACPY
#define zlacpy_           ZLACPY
#define xerbla_           XERBLA
                                                            /* BLACS */
#define blacs_abort_      BLACS_ABORT
#define blacs_gridinfo_   BLACS_GRIDINFO

#define igesd2d_          IGESD2D
#define igebs2d_          IGEBS2D
#define itrsd2d_          ITRSD2D
#define itrbs2d_          ITRBS2D
#define igerv2d_          IGERV2D
#define igebr2d_          IGEBR2D
#define itrrv2d_          ITRRV2D
#define itrbr2d_          ITRBR2D
#define igamx2d_          IGAMX2D
#define igamn2d_          IGAMN2D
#define igsum2d_          IGSUM2D

#define sgesd2d_          SGESD2D
#define sgebs2d_          SGEBS2D
#define strsd2d_          STRSD2D
#define strbs2d_          STRBS2D
#define sgerv2d_          SGERV2D
#define sgebr2d_          SGEBR2D
#define strrv2d_          STRRV2D
#define strbr2d_          STRBR2D
#define sgamx2d_          SGAMX2D
#define sgamn2d_          SGAMN2D
#define sgsum2d_          SGSUM2D

#define dgesd2d_          DGESD2D
#define dgebs2d_          DGEBS2D
#define dtrsd2d_          DTRSD2D
#define dtrbs2d_          DTRBS2D
#define dgerv2d_          DGERV2D
#define dgebr2d_          DGEBR2D
#define dtrrv2d_          DTRRV2D
#define dtrbr2d_          DTRBR2D
#define dgamx2d_          DGAMX2D
#define dgamn2d_          DGAMN2D
#define dgsum2d_          DGSUM2D

#define cgesd2d_          CGESD2D
#define cgebs2d_          CGEBS2D
#define ctrsd2d_          CTRSD2D
#define ctrbs2d_          CTRBS2D
#define cgerv2d_          CGERV2D
#define cgebr2d_          CGEBR2D
#define ctrrv2d_          CTRRV2D
#define ctrbr2d_          CTRBR2D
#define cgamx2d_          CGAMX2D
#define cgamn2d_          CGAMN2D
#define cgsum2d_          CGSUM2D

#define zgesd2d_          ZGESD2D
#define zgebs2d_          ZGEBS2D
#define ztrsd2d_          ZTRSD2D
#define ztrbs2d_          ZTRBS2D
#define zgerv2d_          ZGERV2D
#define zgebr2d_          ZGEBR2D
#define ztrrv2d_          ZTRRV2D
#define ztrbr2d_          ZTRBR2D
#define zgamx2d_          ZGAMX2D
#define zgamn2d_          ZGAMN2D
#define zgsum2d_          ZGSUM2D
                                                     /* Level-1 BLAS */
#define srotg_            SROTG
#define srotmg_           SROTMG
#define srot_             SROT
#define srotm_            SROTM
#define sswap_            SSWAP
#define sscal_            SSCAL
#define scopy_            SCOPY
#define saxpy_            SAXPY
#define ssdot_            SSDOT
#define isamax_           ISAMAX

#define drotg_            DROTG
#define drotmg_           DROTMG
#define drot_             DROT
#define drotm_            DROTM
#define dswap_            DSWAP
#define dscal_            DSCAL
#define dcopy_            DCOPY
#define daxpy_            DAXPY
#define dddot_            DDDOT
#define dnrm2_            DNRM2
#define dsnrm2_           DSNRM2
#define dasum_            DASUM
#define dsasum_           DSASUM
#define idamax_           IDAMAX

#define cswap_            CSWAP
#define cscal_            CSCAL
#define csscal_           CSSCAL
#define ccopy_            CCOPY
#define caxpy_            CAXPY
#define ccdotu_           CCDOTU
#define ccdotc_           CCDOTC
#define icamax_           ICAMAX

#define zswap_            ZSWAP
#define zscal_            ZSCAL
#define zdscal_           ZDSCAL
#define zcopy_            ZCOPY
#define zaxpy_            ZAXPY
#define zzdotu_           ZZDOTU
#define zzdotc_           ZZDOTC
#define dscnrm2_          DSCNRM2
#define dznrm2_           DZNRM2
#define dscasum_          DSCASUM
#define dzasum_           DZASUM
#define izamax_           IZAMAX
                                                     /* Level-2 BLAS */
#define sgemv_            SGEMV
#define ssymv_            SSYMV
#define strmv_            STRMV
#define strsv_            STRSV
#define sger_             SGER
#define ssyr_             SSYR
#define ssyr2_            SSYR2

#define dgemv_            DGEMV
#define dsymv_            DSYMV
#define dtrmv_            DTRMV
#define dtrsv_            DTRSV
#define dger_             DGER
#define dsyr_             DSYR
#define dsyr2_            DSYR2

#define cgemv_            CGEMV
#define chemv_            CHEMV
#define ctrmv_            CTRMV
#define ctrsv_            CTRSV
#define cgeru_            CGERU
#define cgerc_            CGERC
#define cher_             CHER
#define cher2_            CHER2

#define zgemv_            ZGEMV
#define zhemv_            ZHEMV
#define ztrmv_            ZTRMV
#define ztrsv_            ZTRSV
#define zgeru_            ZGERU
#define zgerc_            ZGERC
#define zher_             ZHER
#define zher2_            ZHER2
                                                     /* Level-3 BLAS */
#define sgemm_            SGEMM
#define ssymm_            SSYMM
#define ssyrk_            SSYRK
#define ssyr2k_           SSYR2K
#define strmm_            STRMM
#define strsm_            STRSM

#define dgemm_            DGEMM
#define dsymm_            DSYMM
#define dsyrk_            DSYRK
#define dsyr2k_           DSYR2K
#define dtrmm_            DTRMM
#define dtrsm_            DTRSM

#define cgemm_            CGEMM
#define chemm_            CHEMM
#define csymm_            CSYMM
#define csyrk_            CSYRK
#define cherk_            CHERK
#define csyr2k_           CSYR2K
#define cher2k_           CHER2K
#define ctrmm_            CTRMM
#define ctrsm_            CTRSM

#define zgemm_            ZGEMM
#define zhemm_            ZHEMM
#define zsymm_            ZSYMM
#define zsyrk_            ZSYRK
#define zherk_            ZHERK
#define zsyr2k_           ZSYR2K
#define zher2k_           ZHER2K
#define ztrmm_            ZTRMM
#define ztrsm_            ZTRSM
                                   /* absolute value auxiliary PBLAS */
#define psatrmv_          PSATRMV
#define pdatrmv_          PDATRMV
#define pcatrmv_          PCATRMV
#define pzatrmv_          PZATRMV
#define psagemv_          PSAGEMV
#define pdagemv_          PDAGEMV
#define pcagemv_          PCAGEMV
#define pzagemv_          PZAGEMV
#define psasymv_          PSASYMV
#define pdasymv_          PDASYMV
#define pcahemv_          PCAHEMV
#define pzahemv_          PZAHEMV
                                                /* Auxiliary PB-BLAS */
#define pbcmatadd_        PBCMATADD
#define pbdmatadd_        PBDMATADD
#define pbsmatadd_        PBSMATADD
#define pbzmatadd_        PBZMATADD
                                                   /* Level-2 PBBLAS */
#define pbcgemv_          PBCGEMV
#define pbcgeru_          PBCGERU
#define pbcgerc_          PBCGERC
#define pbchemv_          PBCHEMV
#define pbcher_           PBCHER
#define pbcher2_          PBCHER2
#define pbctrmv_          PBCTRMV
#define pbctrnv_          PBCTRNV
#define pbctrsv_          PBCTRSV

#define pbdgemv_          PBDGEMV
#define pbdger_           PBDGER
#define pbdsymv_          PBDSYMV
#define pbdsyr_           PBDSYR
#define pbdsyr2_          PBDSYR2
#define pbdtrmv_          PBDTRMV
#define pbdtrnv_          PBDTRNV
#define pbdtrsv_          PBDTRSV

#define pbsgemv_          PBSGEMV
#define pbsger_           PBSGER
#define pbssymv_          PBSSYMV
#define pbssyr_           PBSSYR
#define pbssyr2_          PBSSYR2
#define pbstrmv_          PBSTRMV
#define pbstrnv_          PBSTRNV
#define pbstrsv_          PBSTRSV

#define pbzgemv_          PBZGEMV
#define pbzgeru_          PBZGERU
#define pbzgerc_          PBZGERC
#define pbzhemv_          PBZHEMV
#define pbzher_           PBZHER
#define pbzher2_          PBZHER2
#define pbztrmv_          PBZTRMV
#define pbztrnv_          PBZTRNV
#define pbztrsv_          PBZTRSV
                                                   /* Level-3 PBBLAS */
#define pbcgemm_          PBCGEMM
#define pbchemm_          PBCHEMM
#define pbcher2k_         PBCHER2K
#define pbcherk_          PBCHERK
#define pbcsymm_          PBCSYMM
#define pbcsyr2k_         PBCSYR2K
#define pbcsyrk_          PBCSYRK
#define pbctrmm_          PBCTRMM
#define pbctrsm_          PBCTRSM
#define pbctran_          PBCTRAN

#define pbdgemm_          PBDGEMM
#define pbdsymm_          PBDSYMM
#define pbdsyr2k_         PBDSYR2K
#define pbdsyrk_          PBDSYRK
#define pbdtrmm_          PBDTRMM
#define pbdtrsm_          PBDTRSM
#define pbdtran_          PBDTRAN

#define pbsgemm_          PBSGEMM
#define pbssymm_          PBSSYMM
#define pbssyr2k_         PBSSYR2K
#define pbssyrk_          PBSSYRK
#define pbstrmm_          PBSTRMM
#define pbstrsm_          PBSTRSM
#define pbstran_          PBSTRAN

#define pbzgemm_          PBZGEMM
#define pbzhemm_          PBZHEMM
#define pbzher2k_         PBZHER2K
#define pbzherk_          PBZHERK
#define pbzsymm_          PBZSYMM
#define pbzsyr2k_         PBZSYR2K
#define pbzsyrk_          PBZSYRK
#define pbztrmm_          PBZTRMM
#define pbztrsm_          PBZTRSM
#define pbztran_          PBZTRAN
                                                 /* Auxilliary PBLAS */
#define pberror_          PBERROR
#define pb_freebuf_       PB_FREEBUF
#define pb_topget_        PB_TOPGET
#define pb_topset_        PB_TOPSET
                                                    /* Level-1 PBLAS */
#define psrotg_           PSROTG
#define psrotmg_          PSROTMG
#define psrot_            PSROT
#define psrotm_           PSROTM
#define psswap_           PSSWAP
#define psscal_           PSSCAL
#define pscopy_           PSCOPY
#define psaxpy_           PSAXPY
#define psdot_            PSDOT
#define psnrm2_           PSNRM2
#define psasum_           PSASUM
#define psamax_           PSAMAX

#define pdrotg_           PDROTG
#define pdrotmg_          PDROTMG
#define pdrot_            PDROT
#define pdrotm_           PDROTM
#define pdswap_           PDSWAP
#define pdscal_           PDSCAL
#define pdcopy_           PDCOPY
#define pdaxpy_           PDAXPY
#define pddot_            PDDOT
#define pdnrm2_           PDNRM2
#define pdasum_           PDASUM
#define pdamax_           PDAMAX

#define pcswap_           PCSWAP
#define pcscal_           PCSCAL
#define pcsscal_          PCSSCAL
#define pccopy_           PCCOPY
#define pcaxpy_           PCAXPY
#define pcdotu_           PCDOTU
#define pcdotc_           PCDOTC
#define pscnrm2_          PSCNRM2
#define pscasum_          PSCASUM
#define pcamax_           PCAMAX
#define pcrot_            PCROT
#define crot_             CROT

#define pzswap_           PZSWAP
#define pzscal_           PZSCAL
#define pzdscal_          PZDSCAL
#define pzcopy_           PZCOPY
#define pzaxpy_           PZAXPY
#define pzdotu_           PZDOTU
#define pzdotc_           PZDOTC
#define pdznrm2_          PDZNRM2
#define pdzasum_          PDZASUM
#define pzamax_           PZAMAX
#define pzrot_            PZROT
#define zrot_             ZROT
                                                    /* Level-2 PBLAS */
#define pcgemv_           PCGEMV
#define pcgeru_           PCGERU
#define pcgerc_           PCGERC
#define pchemv_           PCHEMV
#define pcher_            PCHER
#define pcher2_           PCHER2
#define pctrmv_           PCTRMV
#define pctrsv_           PCTRSV

#define pdgemv_           PDGEMV
#define pdger_            PDGER
#define pdsymv_           PDSYMV
#define pdsyr_            PDSYR
#define pdsyr2_           PDSYR2
#define pdtrmv_           PDTRMV
#define pdtrsv_           PDTRSV

#define psgemv_           PSGEMV
#define psger_            PSGER
#define pssymv_           PSSYMV
#define pssyr_            PSSYR
#define pssyr2_           PSSYR2
#define pstrmv_           PSTRMV
#define pstrsv_           PSTRSV

#define pzgemv_           PZGEMV
#define pzgeru_           PZGERU
#define pzgerc_           PZGERC
#define pzhemv_           PZHEMV
#define pzher_            PZHER
#define pzher2_           PZHER2
#define pztrmv_           PZTRMV
#define pztrsv_           PZTRSV
                                                    /* Level-3 PBLAS */
#define pcgemm_           PCGEMM
#define pchemm_           PCHEMM
#define pcher2k_          PCHER2K
#define pcherk_           PCHERK
#define pcsymm_           PCSYMM
#define pcsyr2k_          PCSYR2K
#define pcsyrk_           PCSYRK
#define pctrmm_           PCTRMM
#define pctrsm_           PCTRSM
#define pctranu_          PCTRANU
#define pctranc_          PCTRANC

#define pdgemm_           PDGEMM
#define pdsymm_           PDSYMM
#define pdsyr2k_          PDSYR2K
#define pdsyrk_           PDSYRK
#define pdtrmm_           PDTRMM
#define pdtrsm_           PDTRSM
#define pdtran_           PDTRAN

#define psgemm_           PSGEMM
#define pssymm_           PSSYMM
#define pssyr2k_          PSSYR2K
#define pssyrk_           PSSYRK
#define pstrmm_           PSTRMM
#define pstrsm_           PSTRSM
#define pstran_           PSTRAN

#define pzgemm_           PZGEMM
#define pzhemm_           PZHEMM
#define pzher2k_          PZHER2K
#define pzherk_           PZHERK
#define pzsymm_           PZSYMM
#define pzsyr2k_          PZSYR2K
#define pzsyrk_           PZSYRK
#define pztrmm_           PZTRMM
#define pztrsm_           PZTRSM
#define pztranu_          PZTRANU
#define pztranc_          PZTRANC

#endif

#if( _F2C_CALL_ == _F2C_NOCHANGE )
/*
* These defines set up the naming scheme required to have a FORTRAN
* routine call a C routine (which is what the PBLAS are written in)
* for following FORTRAN to C interface:
*           FORTRAN CALL               C DECLARATION
*           call pdgemm(...)           void pdgemm(...)
*/
                                                            /* TOOLS */
#define ilcm_             ilcm
#define infog2l_          infog2l
#define numroc_           numroc
#define pstreecomb_       pstreecomb
#define pdtreecomb_       pdtreecomb
#define pctreecomb_       pctreecomb
#define pztreecomb_       pztreecomb
#define scombamax_        scombamax
#define dcombamax_        dcombamax
#define ccombamax_        ccombamax
#define zcombamax_        zcombamax
#define scombnrm2_        scombnrm2
#define dcombnrm2_        dcombnrm2

#define dlamov_           dlamov
#define slamov_           slamov
#define clamov_           clamov
#define zlamov_           zlamov
#define dlacpy_           dlacpy
#define slacpy_           slacpy
#define clacpy_           clacpy
#define zlacpy_           zlacpy
#define xerbla_           xerbla
                                                            /* BLACS */
#define blacs_abort_      blacs_abort
#define blacs_gridinfo_   blacs_gridinfo

#define igesd2d_          igesd2d
#define igebs2d_          igebs2d
#define itrsd2d_          itrsd2d
#define itrbs2d_          itrbs2d
#define igerv2d_          igerv2d
#define igebr2d_          igebr2d
#define itrrv2d_          itrrv2d
#define itrbr2d_          itrbr2d
#define igamx2d_          igamx2d
#define igamn2d_          igamn2d
#define igsum2d_          igsum2d

#define sgesd2d_          sgesd2d
#define sgebs2d_          sgebs2d
#define strsd2d_          strsd2d
#define strbs2d_          strbs2d
#define sgerv2d_          sgerv2d
#define sgebr2d_          sgebr2d
#define strrv2d_          strrv2d
#define strbr2d_          strbr2d
#define sgamx2d_          sgamx2d
#define sgamn2d_          sgamn2d
#define sgsum2d_          sgsum2d

#define dgesd2d_          dgesd2d
#define dgebs2d_          dgebs2d
#define dtrsd2d_          dtrsd2d
#define dtrbs2d_          dtrbs2d
#define dgerv2d_          dgerv2d
#define dgebr2d_          dgebr2d
#define dtrrv2d_          dtrrv2d
#define dtrbr2d_          dtrbr2d
#define dgamx2d_          dgamx2d
#define dgamn2d_          dgamn2d
#define dgsum2d_          dgsum2d

#define cgesd2d_          cgesd2d
#define cgebs2d_          cgebs2d
#define ctrsd2d_          ctrsd2d
#define ctrbs2d_          ctrbs2d
#define cgerv2d_          cgerv2d
#define cgebr2d_          cgebr2d
#define ctrrv2d_          ctrrv2d
#define ctrbr2d_          ctrbr2d
#define cgamx2d_          cgamx2d
#define cgamn2d_          cgamn2d
#define cgsum2d_          cgsum2d

#define zgesd2d_          zgesd2d
#define zgebs2d_          zgebs2d
#define ztrsd2d_          ztrsd2d
#define ztrbs2d_          ztrbs2d
#define zgerv2d_          zgerv2d
#define zgebr2d_          zgebr2d
#define ztrrv2d_          ztrrv2d
#define ztrbr2d_          ztrbr2d
#define zgamx2d_          zgamx2d
#define zgamn2d_          zgamn2d
#define zgsum2d_          zgsum2d
                                                     /* Level-1 BLAS */
#define srotg_            srotg
#define srotmg_           srotmg
#define srot_             srot
#define srotm_            srotm
#define sswap_            sswap
#define sscal_            sscal
#define scopy_            scopy
#define saxpy_            saxpy
#define ssdot_            ssdot
#define isamax_           isamax

#define drotg_            drotg
#define drotmg_           drotmg
#define drot_             drot
#define drotm_            drotm
#define dswap_            dswap
#define dscal_            dscal
#define dcopy_            dcopy
#define daxpy_            daxpy
#define dddot_            dddot
#define dnrm2_            dnrm2
#define dsnrm2_           dsnrm2
#define dasum_            dasum
#define dsasum_           dsasum
#define idamax_           idamax

#define cswap_            cswap
#define cscal_            cscal
#define csscal_           csscal
#define ccopy_            ccopy
#define caxpy_            caxpy
#define ccdotu_           ccdotu
#define ccdotc_           ccdotc
#define icamax_           icamax

#define zswap_            zswap
#define zscal_            zscal
#define zdscal_           zdscal
#define zcopy_            zcopy
#define zaxpy_            zaxpy
#define zzdotu_           zzdotu
#define zzdotc_           zzdotc
#define dscnrm2_          dscnrm2
#define dznrm2_           dznrm2
#define dscasum_          dscasum
#define dzasum_           dzasum
#define izamax_           izamax
                                                     /* Level-2 BLAS */
#define sgemv_            sgemv
#define ssymv_            ssymv
#define strmv_            strmv
#define strsv_            strsv
#define sger_             sger
#define ssyr_             ssyr
#define ssyr2_            ssyr2

#define dgemv_            dgemv
#define dsymv_            dsymv
#define dtrmv_            dtrmv
#define dtrsv_            dtrsv
#define dger_             dger
#define dsyr_             dsyr
#define dsyr2_            dsyr2

#define cgemv_            cgemv
#define chemv_            chemv
#define ctrmv_            ctrmv
#define ctrsv_            ctrsv
#define cgeru_            cgeru
#define cgerc_            cgerc
#define cher_             cher
#define cher2_            cher2

#define zgemv_            zgemv
#define zhemv_            zhemv
#define ztrmv_            ztrmv
#define ztrsv_            ztrsv
#define zgeru_            zgeru
#define zgerc_            zgerc
#define zher_             zher
#define zher2_            zher2
                                                     /* Level-3 BLAS */
#define sgemm_            sgemm
#define ssymm_            ssymm
#define ssyrk_            ssyrk
#define ssyr2k_           ssyr2k
#define strmm_            strmm
#define strsm_            strsm

#define dgemm_            dgemm
#define dsymm_            dsymm
#define dsyrk_            dsyrk
#define dsyr2k_           dsyr2k
#define dtrmm_            dtrmm
#define dtrsm_            dtrsm

#define cgemm_            cgemm
#define chemm_            chemm
#define csymm_            csymm
#define csyrk_            csyrk
#define cherk_            cherk
#define csyr2k_           csyr2k
#define cher2k_           cher2k
#define ctrmm_            ctrmm
#define ctrsm_            ctrsm

#define zgemm_            zgemm
#define zhemm_            zhemm
#define zsymm_            zsymm
#define zsyrk_            zsyrk
#define zherk_            zherk
#define zsyr2k_           zsyr2k
#define zher2k_           zher2k
#define ztrmm_            ztrmm
#define ztrsm_            ztrsm
                                   /* absolute value auxiliary PBLAS */
#define psatrmv_          psatrmv
#define pdatrmv_          pdatrmv
#define pcatrmv_          pcatrmv
#define pzatrmv_          pzatrmv
#define psagemv_          psagemv
#define pdagemv_          pdagemv
#define pcagemv_          pcagemv
#define pzagemv_          pzagemv
#define psasymv_          psasymv
#define pdasymv_          pdasymv
#define pcahemv_          pcahemv
#define pzahemv_          pzahemv
                                                /* Auxiliary PB-BLAS */
#define pbcmatadd_        pbcmatadd
#define pbdmatadd_        pbdmatadd
#define pbsmatadd_        pbsmatadd
#define pbzmatadd_        pbzmatadd
                                                   /* Level-2 PBBLAS */
#define pbcgemv_          pbcgemv
#define pbcgeru_          pbcgeru
#define pbcgerc_          pbcgerc
#define pbchemv_          pbchemv
#define pbcher_           pbcher
#define pbcher2_          pbcher2
#define pbctrmv_          pbctrmv
#define pbctrnv_          pbctrnv
#define pbctrsv_          pbctrsv

#define pbdgemv_          pbdgemv
#define pbdger_           pbdger
#define pbdsymv_          pbdsymv
#define pbdsyr_           pbdsyr
#define pbdsyr2_          pbdsyr2
#define pbdtrmv_          pbdtrmv
#define pbdtrnv_          pbdtrnv
#define pbdtrsv_          pbdtrsv

#define pbsgemv_          pbsgemv
#define pbsger_           pbsger
#define pbssymv_          pbssymv
#define pbssyr_           pbssyr
#define pbssyr2_          pbssyr2
#define pbstrmv_          pbstrmv
#define pbstrnv_          pbstrnv
#define pbstrsv_          pbstrsv

#define pbzgemv_          pbzgemv
#define pbzgeru_          pbzgeru
#define pbzgerc_          pbzgerc
#define pbzhemv_          pbzhemv
#define pbzher_           pbzher
#define pbzher2_          pbzher2
#define pbztrmv_          pbztrmv
#define pbztrnv_          pbztrnv
#define pbztrsv_          pbztrsv
                                                   /* Level-3 PBBLAS */
#define pbcgemm_          pbcgemm
#define pbchemm_          pbchemm
#define pbcher2k_         pbcher2k
#define pbcherk_          pbcherk
#define pbcsymm_          pbcsymm
#define pbcsyr2k_         pbcsyr2k
#define pbcsyrk_          pbcsyrk
#define pbctrmm_          pbctrmm
#define pbctrsm_          pbctrsm
#define pbctran_          pbctran

#define pbdgemm_          pbdgemm
#define pbdsymm_          pbdsymm
#define pbdsyr2k_         pbdsyr2k
#define pbdsyrk_          pbdsyrk
#define pbdtrmm_          pbdtrmm
#define pbdtrsm_          pbdtrsm
#define pbdtran_          pbdtran

#define pbsgemm_          pbsgemm
#define pbssymm_          pbssymm
#define pbssyr2k_         pbssyr2k
#define pbssyrk_          pbssyrk
#define pbstrmm_          pbstrmm
#define pbstrsm_          pbstrsm
#define pbstran_          pbstran

#define pbzgemm_          pbzgemm
#define pbzhemm_          pbzhemm
#define pbzher2k_         pbzher2k
#define pbzherk_          pbzherk
#define pbzsymm_          pbzsymm
#define pbzsyr2k_         pbzsyr2k
#define pbzsyrk_          pbzsyrk
#define pbztrmm_          pbztrmm
#define pbztrsm_          pbztrsm
#define pbztran_          pbztran
                                                 /* Auxilliary PBLAS */
#define pberror_          pberror
#define pb_freebuf_       pb_freebuf
#define pb_topget_        pb_topget
#define pb_topset_        pb_topset
                                                    /* Level-1 PBLAS */
#define psrotg_           psrotg
#define psrotmg_          psrotmg
#define psrot_            psrot
#define psrotm_           psrotm
#define psswap_           psswap
#define psscal_           psscal
#define pscopy_           pscopy
#define psaxpy_           psaxpy
#define psdot_            psdot
#define psnrm2_           psnrm2
#define psasum_           psasum
#define psamax_           psamax

#define pdrotg_           pdrotg
#define pdrotmg_          pdrotmg
#define pdrot_            pdrot
#define pdrotm_           pdrotm
#define pdswap_           pdswap
#define pdscal_           pdscal
#define pdcopy_           pdcopy
#define pdaxpy_           pdaxpy
#define pddot_            pddot
#define pdnrm2_           pdnrm2
#define pdasum_           pdasum
#define pdamax_           pdamax

#define pcswap_           pcswap
#define pcscal_           pcscal
#define pcsscal_          pcsscal
#define pccopy_           pccopy
#define pcaxpy_           pcaxpy
#define pcdotu_           pcdotu
#define pcdotc_           pcdotc
#define pscnrm2_          pscnrm2
#define pscasum_          pscasum
#define pcamax_           pcamax
#define pcrot_            pcrot
#define crot_             crot

#define pzswap_           pzswap
#define pzscal_           pzscal
#define pzdscal_          pzdscal
#define pzcopy_           pzcopy
#define pzaxpy_           pzaxpy
#define pzdotu_           pzdotu
#define pzdotc_           pzdotc
#define pdznrm2_          pdznrm2
#define pdzasum_          pdzasum
#define pzamax_           pzamax
#define pzrot_            pzrot
#define zrot_             zrot
                                                    /* Level-2 PBLAS */
#define pcgemv_           pcgemv
#define pcgeru_           pcgeru
#define pcgerc_           pcgerc
#define pchemv_           pchemv
#define pcher_            pcher
#define pcher2_           pcher2
#define pctrmv_           pctrmv
#define pctrsv_           pctrsv

#define pdgemv_           pdgemv
#define pdger_            pdger
#define pdsymv_           pdsymv
#define pdsyr_            pdsyr
#define pdsyr2_           pdsyr2
#define pdtrmv_           pdtrmv
#define pdtrsv_           pdtrsv

#define psgemv_           psgemv
#define psger_            psger
#define pssymv_           pssymv
#define pssyr_            pssyr
#define pssyr2_           pssyr2
#define pstrmv_           pstrmv
#define pstrsv_           pstrsv

#define pzgemv_           pzgemv
#define pzgeru_           pzgeru
#define pzgerc_           pzgerc
#define pzhemv_           pzhemv
#define pzher_            pzher
#define pzher2_           pzher2
#define pztrmv_           pztrmv
#define pztrsv_           pztrsv
                                                    /* Level-3 PBLAS */
#define pcgemm_           pcgemm
#define pchemm_           pchemm
#define pcher2k_          pcher2k
#define pcherk_           pcherk
#define pcsymm_           pcsymm
#define pcsyr2k_          pcsyr2k
#define pcsyrk_           pcsyrk
#define pctrmm_           pctrmm
#define pctrsm_           pctrsm
#define pctranu_          pctranu
#define pctranc_          pctranc

#define pdgemm_           pdgemm
#define pdsymm_           pdsymm
#define pdsyr2k_          pdsyr2k
#define pdsyrk_           pdsyrk
#define pdtrmm_           pdtrmm
#define pdtrsm_           pdtrsm
#define pdtran_           pdtran

#define psgemm_           psgemm
#define pssymm_           pssymm
#define pssyr2k_          pssyr2k
#define pssyrk_           pssyrk
#define pstrmm_           pstrmm
#define pstrsm_           pstrsm
#define pstran_           pstran

#define pzgemm_           pzgemm
#define pzhemm_           pzhemm
#define pzher2k_          pzher2k
#define pzherk_           pzherk
#define pzsymm_           pzsymm
#define pzsyr2k_          pzsyr2k
#define pzsyrk_           pzsyrk
#define pztrmm_           pztrmm
#define pztrsm_           pztrsm
#define pztranu_          pztranu
#define pztranc_          pztranc

#endif

#if( _F2C_CALL_ == _F2C_F77ISF2C )
/*
* These defines set up the naming scheme required to have a FORTRAN
* routine call a C routine (which is what the PBLAS are written in)
* for systems where the fortran "compiler" is actually f2c (a Fortran
* to C conversion utility).
*/

/*
* Initialization routines
*/
#define blacs_pinfo_    blacs_pinfo__
#define blacs_setup_    blacs_setup__
#define blacs_set_      blacs_set__
#define blacs_get_      blacs_get__
#define blacs_gridinit_ blacs_gridinit__
#define blacs_gridmap_  blacs_gridmap__
/*
* Destruction routines
*/
#define blacs_freebuff_ blacs_freebuff__
#define blacs_gridexit_ blacs_gridexit__
#define blacs_abort_    blacs_abort__
#define blacs_exit_     blacs_exit__
/*
* Informational & misc.
*/
#define blacs_gridinfo_ blacs_gridinfo__
#define blacs_pnum_     blacs_pnum__
#define blacs_pcoord_   blacs_pcoord__
#define blacs_barrier_  blacs_barrier__
#endif
