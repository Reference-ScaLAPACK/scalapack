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
*  This file includes PBLAS definitions. All PBLAS routines include this
*  file.
*
*  ---------------------------------------------------------------------
*  #define macro constants
*  ---------------------------------------------------------------------
*/
#if( _F2C_CALL_ == _F2C_ADD_ )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine call a C routine. No redefinition is necessary  to  have  the
*  following FORTRAN to C interface:
*
*           FORTRAN CALL                   C DECLARATION
*           CALL PDGEMM(...)               void pdgemm_(...)
*
*  This is the PBLAS default.
*/
#define    PB_freebuf_         PB_freebuf_
#define    PB_topget_          pb_topget_
#define    PB_topset_          pb_topset_

#endif

#if( _F2C_CALL_ == _F2C_UPCASE )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine call a C routine. No redefinition is necessary  to  have  the
*  following FORTRAN to C interface:
*
*           FORTRAN CALL                   C DECLARATION
*           CALL PDGEMM(...)               void PDGEMM(...)
*/
#define    pilaenv_            PILAENV
#define    PB_freebuf_         PB_FREEBUF
#define    PB_topget_          PB_TOPGET
#define    PB_topset_          PB_TOPSET
                                                     /* Level-1 PBLAS */
#define    picopy_             PICOPY
#define    pscopy_             PSCOPY
#define    pdcopy_             PDCOPY
#define    pccopy_             PCCOPY
#define    pzcopy_             PZCOPY

#define    psswap_             PSSWAP
#define    pdswap_             PDSWAP
#define    pcswap_             PCSWAP
#define    pzswap_             PZSWAP

#define    psaxpy_             PSAXPY
#define    pdaxpy_             PDAXPY
#define    pcaxpy_             PCAXPY
#define    pzaxpy_             PZAXPY

#define    psscal_             PSSCAL
#define    pdscal_             PDSCAL
#define    pcscal_             PCSCAL
#define    pzscal_             PZSCAL
#define    pcsscal_            PCSSCAL
#define    pzdscal_            PZDSCAL

#define    psasum_             PSASUM
#define    pdasum_             PDASUM
#define    pscasum_            PSCASUM
#define    pdzasum_            PDZASUM

#define    psnrm2_             PSNRM2
#define    pdnrm2_             PDNRM2
#define    pscnrm2_            PSCNRM2
#define    pdznrm2_            PDZNRM2

#define    psdot_              PSDOT
#define    pddot_              PDDOT
#define    pcdotu_             PCDOTU
#define    pzdotu_             PZDOTU
#define    pcdotc_             PCDOTC
#define    pzdotc_             PZDOTC

#define    psamax_             PSAMAX
#define    pdamax_             PDAMAX
#define    pcamax_             PCAMAX
#define    pzamax_             PZAMAX

#define    psgemv_             PSGEMV
#define    pdgemv_             PDGEMV
#define    pcgemv_             PCGEMV
#define    pzgemv_             PZGEMV

#define    psagemv_            PSAGEMV
#define    pdagemv_            PDAGEMV
#define    pcagemv_            PCAGEMV
#define    pzagemv_            PZAGEMV

#define    pssymv_             PSSYMV
#define    pdsymv_             PDSYMV
#define    pchemv_             PCHEMV
#define    pzhemv_             PZHEMV

#define    psasymv_            PSASYMV
#define    pdasymv_            PDASYMV
#define    pcahemv_            PCAHEMV
#define    pzahemv_            PZAHEMV

#define    pstrmv_             PSTRMV
#define    pdtrmv_             PDTRMV
#define    pctrmv_             PCTRMV
#define    pztrmv_             PZTRMV

#define    psatrmv_            PSATRMV
#define    pdatrmv_            PDATRMV
#define    pcatrmv_            PCATRMV
#define    pzatrmv_            PZATRMV

#define    pstrsv_             PSTRSV
#define    pdtrsv_             PDTRSV
#define    pctrsv_             PCTRSV
#define    pztrsv_             PZTRSV

#define    psger_              PSGER
#define    pdger_              PDGER
#define    pcgeru_             PCGERU
#define    pzgeru_             PZGERU
#define    pcgerc_             PCGERC
#define    pzgerc_             PZGERC

#define    pssyr_              PSSYR
#define    pdsyr_              PDSYR
#define    pcher_              PCHER
#define    pzher_              PZHER

#define    pssyr2_             PSSYR2
#define    pdsyr2_             PDSYR2
#define    pcher2_             PCHER2
#define    pzher2_             PZHER2

#define    psgemm_             PSGEMM
#define    pdgemm_             PDGEMM
#define    pcgemm_             PCGEMM
#define    pzgemm_             PZGEMM

#define    psgeadd_            PSGEADD
#define    pdgeadd_            PDGEADD
#define    pcgeadd_            PCGEADD
#define    pzgeadd_            PZGEADD

#define    pssymm_             PSSYMM
#define    pdsymm_             PDSYMM
#define    pcsymm_             PCSYMM
#define    pchemm_             PCHEMM
#define    pzsymm_             PZSYMM
#define    pzhemm_             PZHEMM

#define    pstrmm_             PSTRMM
#define    pdtrmm_             PDTRMM
#define    pctrmm_             PCTRMM
#define    pztrmm_             PZTRMM

#define    pstrsm_             PSTRSM
#define    pdtrsm_             PDTRSM
#define    pctrsm_             PCTRSM
#define    pztrsm_             PZTRSM

#define    pssyrk_             PSSYRK
#define    pdsyrk_             PDSYRK
#define    pcsyrk_             PCSYRK
#define    pcherk_             PCHERK
#define    pzsyrk_             PZSYRK
#define    pzherk_             PZHERK

#define    pssyr2k_            PSSYR2K
#define    pdsyr2k_            PDSYR2K
#define    pcsyr2k_            PCSYR2K
#define    pcher2k_            PCHER2K
#define    pzsyr2k_            PZSYR2K
#define    pzher2k_            PZHER2K

#define    pstradd_            PSTRADD
#define    pdtradd_            PDTRADD
#define    pctradd_            PCTRADD
#define    pztradd_            PZTRADD

#define    pstran_             PSTRAN
#define    pdtran_             PDTRAN
#define    pctranu_            PCTRANU
#define    pztranu_            PZTRANU
#define    pctranc_            PCTRANC
#define    pztranc_            PZTRANC

#endif

#if( _F2C_CALL_ == _F2C_NOCHANGE )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine call a C routine with the following  FORTRAN to C interface:
*
*           FORTRAN CALL                   C DECLARATION
*           CALLL PDGEMM(...)              void pdgemm(...)
*/
#define    pilaenv_            pilaenv
#define    PB_freebuf_         PB_freebuf
#define    PB_topget_          pb_topget
#define    PB_topset_          pb_topset

#define    picopy_             picopy
#define    pscopy_             pscopy
#define    pdcopy_             pdcopy
#define    pccopy_             pccopy
#define    pzcopy_             pzcopy

#define    psswap_             psswap
#define    pdswap_             pdswap
#define    pcswap_             pcswap
#define    pzswap_             pzswap

#define    psaxpy_             psaxpy
#define    pdaxpy_             pdaxpy
#define    pcaxpy_             pcaxpy
#define    pzaxpy_             pzaxpy

#define    psscal_             psscal
#define    pdscal_             pdscal
#define    pcscal_             pcscal
#define    pzscal_             pzscal
#define    pcsscal_            pcsscal
#define    pzdscal_            pzdscal

#define    psasum_             psasum
#define    pdasum_             pdasum
#define    pscasum_            pscasum
#define    pdzasum_            pdzasum

#define    psnrm2_             psnrm2
#define    pdnrm2_             pdnrm2
#define    pscnrm2_            pscnrm2
#define    pdznrm2_            pdznrm2

#define    psdot_              psdot
#define    pddot_              pddot
#define    pcdotu_             pcdotu
#define    pzdotu_             pzdotu
#define    pcdotc_             pcdotc
#define    pzdotc_             pzdotc

#define    psamax_             psamax
#define    pdamax_             pdamax
#define    pcamax_             pcamax
#define    pzamax_             pzamax

#define    psgemv_             psgemv
#define    pdgemv_             pdgemv
#define    pcgemv_             pcgemv
#define    pzgemv_             pzgemv

#define    psagemv_            psagemv
#define    pdagemv_            pdagemv
#define    pcagemv_            pcagemv
#define    pzagemv_            pzagemv

#define    pssymv_             pssymv
#define    pdsymv_             pdsymv
#define    pchemv_             pchemv
#define    pzhemv_             pzhemv

#define    psasymv_            psasymv
#define    pdasymv_            pdasymv
#define    pcahemv_            pcahemv
#define    pzahemv_            pzahemv

#define    pstrmv_             pstrmv
#define    pdtrmv_             pdtrmv
#define    pctrmv_             pctrmv
#define    pztrmv_             pztrmv

#define    psatrmv_            psatrmv
#define    pdatrmv_            pdatrmv
#define    pcatrmv_            pcatrmv
#define    pzatrmv_            pzatrmv

#define    pstrsv_             pstrsv
#define    pdtrsv_             pdtrsv
#define    pctrsv_             pctrsv
#define    pztrsv_             pztrsv

#define    psger_              psger
#define    pdger_              pdger
#define    pcgeru_             pcgeru
#define    pzgeru_             pzgeru
#define    pcgerc_             pcgerc
#define    pzgerc_             pzgerc

#define    pssyr_              pssyr
#define    pdsyr_              pdsyr
#define    pcher_              pcher
#define    pzher_              pzher

#define    pssyr2_             pssyr2
#define    pdsyr2_             pdsyr2
#define    pcher2_             pcher2
#define    pzher2_             pzher2

#define    psgeadd_            psgeadd
#define    pdgeadd_            pdgeadd
#define    pcgeadd_            pcgeadd
#define    pzgeadd_            pzgeadd

#define    psgemm_             psgemm
#define    pdgemm_             pdgemm
#define    pcgemm_             pcgemm
#define    pzgemm_             pzgemm

#define    pssymm_             pssymm
#define    pdsymm_             pdsymm
#define    pcsymm_             pcsymm
#define    pchemm_             pchemm
#define    pzsymm_             pzsymm
#define    pzhemm_             pzhemm

#define    pstrmm_             pstrmm
#define    pdtrmm_             pdtrmm
#define    pctrmm_             pctrmm
#define    pztrmm_             pztrmm

#define    pstrsm_             pstrsm
#define    pdtrsm_             pdtrsm
#define    pctrsm_             pctrsm
#define    pztrsm_             pztrsm

#define    pssyrk_             pssyrk
#define    pdsyrk_             pdsyrk
#define    pcsyrk_             pcsyrk
#define    pcherk_             pcherk
#define    pzsyrk_             pzsyrk
#define    pzherk_             pzherk

#define    pssyr2k_            pssyr2k
#define    pdsyr2k_            pdsyr2k
#define    pcsyr2k_            pcsyr2k
#define    pcher2k_            pcher2k
#define    pzsyr2k_            pzsyr2k
#define    pzher2k_            pzher2k

#define    pstradd_            pstradd
#define    pdtradd_            pdtradd
#define    pctradd_            pctradd
#define    pztradd_            pztradd

#define    pstran_             pstran
#define    pdtran_             pdtran
#define    pctranu_            pctranu
#define    pztranu_            pztranu
#define    pctranc_            pctranc
#define    pztranc_            pztranc

#endif

#if( _F2C_CALL_ == _F2C_F77ISF2C )

#define    PB_freebuf_         PB_freebuf__
#define    PB_topget_          pb_topget__
#define    PB_topset_          pb_topset__

#endif
/*
*  ---------------------------------------------------------------------
*  Function prototypes
*  ---------------------------------------------------------------------
*/
#ifdef __STDC__

void           PB_freebuf_     ( void );

void           PB_topget_      ( int *,     F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T );

void           PB_topset_      ( int *,     F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T );

void           picopy_         ( int *,     int *,     int *,
                                 int *,     int *,     int *,
                                 int *,     int *,     int *,
                                 int *,     int * );
void           pscopy_         ( int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pdcopy_         ( int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     int * );
void           pccopy_         ( int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pzcopy_         ( int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     int * );

void           psswap_         ( int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pdswap_         ( int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     int * );
void           pcswap_         ( int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pzswap_         ( int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     int * );

void           psaxpy_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int * );
void           pdaxpy_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int * );
void           pcaxpy_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int * );
void           pzaxpy_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int * );

void           psscal_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdscal_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );
void           pcscal_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pcsscal_        ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pzscal_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );
void           pzdscal_        ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );

void           psasum_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdasum_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );
void           pscasum_        ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdzasum_        ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );

void           psnrm2_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdnrm2_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );
void           pscnrm2_        ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdznrm2_        ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );

void           psdot_          ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int * );
void           pddot_          ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int * );
void           pcdotc_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int * );
void           pcdotu_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int * );
void           pzdotc_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int * );
void           pzdotu_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int * );

void           psamax_         ( int *,     float *,   int *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pdamax_         ( int *,     double *,  int *,
                                 double *,  int *,     int *,
                                 int *,     int * );
void           pcamax_         ( int *,     float *,   int *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pzamax_         ( int *,     double *,  int *,
                                 double *,  int *,     int *,
                                 int *,     int * );

void           psgemv_         ( F_CHAR_T,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdgemv_         ( F_CHAR_T,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );
void           pcgemv_         ( F_CHAR_T,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pzgemv_         ( F_CHAR_T,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );

void           psagemv_        ( F_CHAR_T,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdagemv_        ( F_CHAR_T,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );
void           pcagemv_        ( F_CHAR_T,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pzagemv_        ( F_CHAR_T,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );

void           psger_          ( int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int * );
void           pdger_          ( int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int * );
void           pcgerc_         ( int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int * );
void           pcgeru_         ( int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int * );
void           pzgerc_         ( int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int * );
void           pzgeru_         ( int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int * );

void           pssymv_         ( F_CHAR_T,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     int * );
void           pdsymv_         ( F_CHAR_T,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     int * );
void           pchemv_         ( F_CHAR_T,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     int * );
void           pzhemv_         ( F_CHAR_T,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     int * );

void           psasymv_        ( F_CHAR_T,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     int * );
void           pdasymv_        ( F_CHAR_T,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     int * );
void           pcahemv_        ( F_CHAR_T,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     int * );
void           pzahemv_        ( F_CHAR_T,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     int * );

void           pssyr_          ( F_CHAR_T,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int * );
void           pdsyr_          ( F_CHAR_T,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int * );
void           pcher_          ( F_CHAR_T,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int * );
void           pzher_          ( F_CHAR_T,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int * );

void           pssyr2_         ( F_CHAR_T,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int * );
void           pdsyr2_         ( F_CHAR_T,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int * );
void           pcher2_         ( F_CHAR_T,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int * );
void           pzher2_         ( F_CHAR_T,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int * );

void           pstrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 int *,     float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdtrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 int *,     double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int * );
void           pctrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 int *,     float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int * );
void           pztrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 int *,     double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int * );

void           psatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pdatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int * );
void           pcatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pzatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int * );

void           pstrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 int *,     float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdtrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 int *,     double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int * );
void           pctrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 int *,     float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int * );
void           pztrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 int *,     double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int * );

void           psgeadd_        ( F_CHAR_T,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int * );
void           pdgeadd_        ( F_CHAR_T,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int * );
void           pcgeadd_        ( F_CHAR_T,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int * );
void           pzgeadd_        ( F_CHAR_T,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int * );

void           psgemm_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int * );
void           pdgemm_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int * );
void           pcgemm_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int * );
void           pzgemm_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int * );

void           pssymm_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pdsymm_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );
void           pcsymm_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pzsymm_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );
void           pchemm_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pzhemm_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );

void           pssyr2k_        ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pdsyr2k_        ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );
void           pcsyr2k_        ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pzsyr2k_        ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );
void           pcher2k_        ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pzher2k_        ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );

void           pssyrk_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int * );
void           pdsyrk_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int * );
void           pcsyrk_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int * );
void           pzsyrk_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int * );
void           pcherk_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int * );
void           pzherk_         ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int * );

void           pstradd_        ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int * );
void           pdtradd_        ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int * );
void           pctradd_        ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int * );
void           pztradd_        ( F_CHAR_T,  F_CHAR_T,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int * );

void           pstran_         ( int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pdtran_         ( int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );
void           pctranc_        ( int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pztranc_        ( int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );
void           pctranu_        ( int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pztranu_        ( int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );

void           pstrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int * );
void           pdtrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int * );
void           pctrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int * );
void           pztrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int * );

void           pstrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int * );
void           pdtrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int * );
void           pctrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int * );
void           pztrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int * );
#else

void           PB_freebuf_     ();
void           PB_topget_      ();
void           PB_topset_      ();

void           picopy_         ();
void           pscopy_         ();
void           pdcopy_         ();
void           pccopy_         ();
void           pzcopy_         ();

void           psswap_         ();
void           pdswap_         ();
void           pcswap_         ();
void           pzswap_         ();

void           psaxpy_         ();
void           pdaxpy_         ();
void           pcaxpy_         ();
void           pzaxpy_         ();

void           psscal_         ();
void           pdscal_         ();
void           pcscal_         ();
void           pcsscal_        ();
void           pzscal_         ();
void           pzdscal_        ();

void           psasum_         ();
void           pdasum_         ();
void           pscasum_        ();
void           pdzasum_        ();

void           psnrm2_         ();
void           pdnrm2_         ();
void           pscnrm2_        ();
void           pdznrm2_        ();

void           psdot_          ();
void           pddot_          ();
void           pcdotc_         ();
void           pcdotu_         ();
void           pzdotc_         ();
void           pzdotu_         ();

void           psamax_         ();
void           pdamax_         ();
void           pcamax_         ();
void           pzamax_         ();

void           psgemv_         ();
void           pdgemv_         ();
void           pcgemv_         ();
void           pzgemv_         ();

void           psagemv_        ();
void           pdagemv_        ();
void           pcagemv_        ();
void           pzagemv_        ();

void           psger_          ();
void           pdger_          ();
void           pcgerc_         ();
void           pcgeru_         ();
void           pzgerc_         ();
void           pzgeru_         ();

void           pssymv_         ();
void           pdsymv_         ();
void           pchemv_         ();
void           pzhemv_         ();

void           psasymv_        ();
void           pdasymv_        ();
void           pcahemv_        ();
void           pzahemv_        ();

void           pssyr_          ();
void           pdsyr_          ();
void           pcher_          ();
void           pzher_          ();

void           pssyr2_         ();
void           pdsyr2_         ();
void           pcher2_         ();
void           pzher2_         ();

void           pstrmv_         ();
void           pdtrmv_         ();
void           pctrmv_         ();
void           pztrmv_         ();

void           psatrmv_        ();
void           pdatrmv_        ();
void           pcatrmv_        ();
void           pzatrmv_        ();

void           pstrsv_         ();
void           pdtrsv_         ();
void           pctrsv_         ();
void           pztrsv_         ();

void           psgeadd_        ();
void           pdgeadd_        ();
void           pcgeadd_        ();
void           pzgeadd_        ();

void           psgemm_         ();
void           pdgemm_         ();
void           pcgemm_         ();
void           pzgemm_         ();

void           pssymm_         ();
void           pdsymm_         ();
void           pcsymm_         ();
void           pchemm_         ();
void           pzsymm_         ();
void           pzhemm_         ();

void           pssyr2k_        ();
void           pdsyr2k_        ();
void           pcsyr2k_        ();
void           pcher2k_        ();
void           pzsyr2k_        ();
void           pzher2k_        ();

void           pssyrk_         ();
void           pdsyrk_         ();
void           pcsyrk_         ();
void           pcherk_         ();
void           pzsyrk_         ();
void           pzherk_         ();

void           pstradd_        ();
void           pdtradd_        ();
void           pctradd_        ();
void           pztradd_        ();

void           pstran_         ();
void           pdtran_         ();
void           pctranc_        ();
void           pctranu_        ();
void           pztranc_        ();
void           pztranu_        ();

void           pstrmm_         ();
void           pdtrmm_         ();
void           pctrmm_         ();
void           pztrmm_         ();

void           pstrsm_         ();
void           pdtrsm_         ();
void           pctrsm_         ();
void           pztrsm_         ();

#endif
