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

void           PB_topget_      ( Int *,     F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T );

void           PB_topset_      ( Int *,     F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T );

void           picopy_         ( Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int * );
void           pscopy_         ( Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pdcopy_         ( Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );
void           pccopy_         ( Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pzcopy_         ( Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );

void           psswap_         ( Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pdswap_         ( Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );
void           pcswap_         ( Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pzswap_         ( Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );

void           psaxpy_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pdaxpy_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int * );
void           pcaxpy_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pzaxpy_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int * );

void           psscal_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdscal_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pcscal_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pcsscal_        ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pzscal_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pzdscal_        ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psasum_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdasum_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pscasum_        ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdzasum_        ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psnrm2_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdnrm2_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pscnrm2_        ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdznrm2_        ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psdot_          ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pddot_          ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int * );
void           pcdotc_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pcdotu_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pzdotc_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int * );
void           pzdotu_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int * );

void           psamax_         ( Int *,     float *,   Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pdamax_         ( Int *,     double *,  Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );
void           pcamax_         ( Int *,     float *,   Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pzamax_         ( Int *,     double *,  Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );

void           psgemv_         ( F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdgemv_         ( F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pcgemv_         ( F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pzgemv_         ( F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psagemv_        ( F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdagemv_        ( F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pcagemv_        ( F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pzagemv_        ( F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psger_          ( Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int * );
void           pdger_          ( Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int * );
void           pcgerc_         ( Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int * );
void           pcgeru_         ( Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int * );
void           pzgerc_         ( Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int * );
void           pzgeru_         ( Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int * );

void           pssymv_         ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pdsymv_         ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     Int * );
void           pchemv_         ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pzhemv_         ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     Int * );

void           psasymv_        ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pdasymv_        ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     Int * );
void           pcahemv_        ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pzahemv_        ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     Int * );

void           pssyr_          ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int * );
void           pdsyr_          ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int * );
void           pcher_          ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int * );
void           pzher_          ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int * );

void           pssyr2_         ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int * );
void           pdsyr2_         ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int * );
void           pcher2_         ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int * );
void           pzher2_         ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int * );

void           pstrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdtrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pctrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pztrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pdatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );
void           pcatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pzatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );

void           pstrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdtrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pctrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pztrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psgeadd_        ( F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int * );
void           pdgeadd_        ( F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int * );
void           pcgeadd_        ( F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int * );
void           pzgeadd_        ( F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int * );

void           psgemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int * );
void           pdgemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int * );
void           pcgemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int * );
void           pzgemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int * );

void           pssymm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pdsymm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pcsymm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pzsymm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pchemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pzhemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );

void           pssyr2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pdsyr2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pcsyr2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pzsyr2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pcher2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pzher2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );

void           pssyrk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int * );
void           pdsyrk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int * );
void           pcsyrk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int * );
void           pzsyrk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int * );
void           pcherk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int * );
void           pzherk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int * );

void           pstradd_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int * );
void           pdtradd_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int * );
void           pctradd_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int * );
void           pztradd_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int * );

void           pstran_         ( Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pdtran_         ( Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pctranc_        ( Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pztranc_        ( Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pctranu_        ( Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pztranu_        ( Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );

void           pstrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int * );
void           pdtrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int * );
void           pctrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int * );
void           pztrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int * );

void           pstrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int * );
void           pdtrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int * );
void           pctrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int * );
void           pztrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int * );
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
