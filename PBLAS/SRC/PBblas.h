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
*  This file includes F77 BLAS definitions. All  PBLAS  routines include
*  this file.
*
*  ---------------------------------------------------------------------
*  #define macro constants
*  ---------------------------------------------------------------------
*/
#define    CNOTRAN             'N'
#define    CNOCONJG            'N'
#define    CTRAN               'T'
#define    CCONJG              'Z'
#define    CCOTRAN             'C'

#define    CALL                'A'
#define    CLOWER              'L'
#define    CUPPER              'U'
#define    CDIAGONAL           'D'

#define    CLEFT               'L'
#define    CRIGHT              'R'

#define    CUNIT               'U'
#define    CNOUNIT             'N'

#define    CINIT               'I'
#define    CNOINIT             'N'

#define    CFORWARD            'F'
#define    CBACKWARD           'B'

#define    CREUSE              'R'
#define    CALLOCATE           'A'

#define    NOTRAN              "N"
#define    NOCONJG             "N"
#define    TRAN                "T"
#define    CONJG               "Z"
#define    COTRAN              "C"

#define    ALL                 "A"
#define    LOWER               "L"
#define    UPPER               "U"
#define    DIAGONAL            "D"

#define    LEFT                "L"
#define    RIGHT               "R"

#define    UNIT                "U"
#define    NOUNIT              "N"

#define    INIT                "I"
#define    NOINIT              "N"

#define    FORWARD             "F"
#define    BACKWARD            "B"

#define    REUSE               "R"
#define    ALLOCATE            "A"

#if( _F2C_CALL_ == _F2C_ADD_ )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine. No redefinition is necessary  to  have
*  the following FORTRAN to C interface:
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE DGEMM(...)          dgemm_(...)
*
*  This is the PBLAS default.
*/
#endif

#if( _F2C_CALL_ == _F2C_UPCASE )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine with the following  FORTRAN to C inter-
*  face:
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE DGEMM(...)          DGEMM(...)
*/
#define    srot_               SROT
#define    drot_               DROT

#define    sswap_              SSWAP
#define    dswap_              DSWAP
#define    cswap_              CSWAP
#define    zswap_              ZSWAP

#define    scopy_              SCOPY
#define    dcopy_              DCOPY
#define    ccopy_              CCOPY
#define    zcopy_              ZCOPY

#define    saxpy_              SAXPY
#define    daxpy_              DAXPY
#define    caxpy_              CAXPY
#define    zaxpy_              ZAXPY

#define    sscal_              SSCAL
#define    dscal_              DSCAL
#define    cscal_              CSCAL
#define    zscal_              ZSCAL
#define    csscal_             CSSCAL
#define    zdscal_             ZDSCAL

#define    sasum_              SASUM
#define    dasum_              DASUM
#define    scasum_             SCASUM
#define    dzasum_             DZASUM

#define    snrm2_              SNRM2
#define    dnrm2_              DNRM2
#define    scnrm2_             SCNRM2
#define    dznrm2_             DZNRM2

#define    sdot_               SDOT
#define    ddot_               DDOT
#define    cdotu_              CDOTU
#define    zdotu_              ZDOTU
#define    cdotc_              CDOTC
#define    zdotc_              ZDOTC

#define    isamax_             ISAMAX
#define    idamax_             IDAMAX
#define    icamax_             ICAMAX
#define    izamax_             IZAMAX

#define    sgemv_              SGEMV
#define    dgemv_              DGEMV
#define    cgemv_              CGEMV
#define    zgemv_              ZGEMV

#define    ssymv_              SSYMV
#define    dsymv_              DSYMV
#define    chemv_              CHEMV
#define    zhemv_              ZHEMV

#define    strmv_              STRMV
#define    dtrmv_              DTRMV
#define    ctrmv_              CTRMV
#define    ztrmv_              ZTRMV

#define    strsv_              STRSV
#define    dtrsv_              DTRSV
#define    ctrsv_              CTRSV
#define    ztrsv_              ZTRSV

#define    sger_               SGER
#define    dger_               DGER
#define    cgeru_              CGERU
#define    zgeru_              ZGERU
#define    cgerc_              CGERC
#define    zgerc_              ZGERC

#define    ssyr_               SSYR
#define    dsyr_               DSYR
#define    cher_               CHER
#define    zher_               ZHER

#define    ssyr2_              SSYR2
#define    dsyr2_              DSYR2
#define    cher2_              CHER2
#define    zher2_              ZHER2

#define    sgemm_              SGEMM
#define    dgemm_              DGEMM
#define    cgemm_              CGEMM
#define    zgemm_              ZGEMM

#define    ssymm_              SSYMM
#define    dsymm_              DSYMM
#define    csymm_              CSYMM
#define    chemm_              CHEMM
#define    zsymm_              ZSYMM
#define    zhemm_              ZHEMM

#define    strmm_              STRMM
#define    dtrmm_              DTRMM
#define    ctrmm_              CTRMM
#define    ztrmm_              ZTRMM

#define    strsm_              STRSM
#define    dtrsm_              DTRSM
#define    ctrsm_              CTRSM
#define    ztrsm_              ZTRSM

#define    ssyrk_              SSYRK
#define    dsyrk_              DSYRK
#define    csyrk_              CSYRK
#define    cherk_              CHERK
#define    zsyrk_              ZSYRK
#define    zherk_              ZHERK

#define    ssyr2k_             SSYR2K
#define    dsyr2k_             DSYR2K
#define    csyr2k_             CSYR2K
#define    cher2k_             CHER2K
#define    zsyr2k_             ZSYR2K
#define    zher2k_             ZHER2K

#endif

#if( _F2C_CALL_ == _F2C_NOCHANGE )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine with the following  FORTRAN to C inter-
*  face:
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE DGEMM(...)          dgemm(...)
*/
#define    srot_               srot
#define    drot_               drot

#define    sswap_              sswap
#define    dswap_              dswap
#define    cswap_              cswap
#define    zswap_              zswap

#define    scopy_              scopy
#define    dcopy_              dcopy
#define    ccopy_              ccopy
#define    zcopy_              zcopy

#define    saxpy_              saxpy
#define    daxpy_              daxpy
#define    caxpy_              caxpy
#define    zaxpy_              zaxpy

#define    sscal_              sscal
#define    dscal_              dscal
#define    cscal_              cscal
#define    zscal_              zscal
#define    csscal_             csscal
#define    zdscal_             zdscal

#define    sasum_              sasum
#define    dasum_              dasum
#define    scasum_             scasum
#define    dzasum_             dzasum

#define    snrm2_              snrm2
#define    dnrm2_              dnrm2
#define    scnrm2_             scnrm2
#define    dznrm2_             dznrm2

#define    sdot_               sdot
#define    ddot_               ddot
#define    cdotu_              cdotu
#define    zdotu_              zdotu
#define    cdotc_              cdotc
#define    zdotc_              zdotc

#define    isamax_             isamax
#define    idamax_             idamax
#define    icamax_             icamax
#define    izamax_             izamax

#define    sgemv_              sgemv
#define    dgemv_              dgemv
#define    cgemv_              cgemv
#define    zgemv_              zgemv

#define    ssymv_              ssymv
#define    dsymv_              dsymv
#define    chemv_              chemv
#define    zhemv_              zhemv

#define    strmv_              strmv
#define    dtrmv_              dtrmv
#define    ctrmv_              ctrmv
#define    ztrmv_              ztrmv

#define    strsv_              strsv
#define    dtrsv_              dtrsv
#define    ctrsv_              ctrsv
#define    ztrsv_              ztrsv

#define    sger_               sger
#define    dger_               dger
#define    cgeru_              cgeru
#define    zgeru_              zgeru
#define    cgerc_              cgerc
#define    zgerc_              zgerc

#define    ssyr_               ssyr
#define    dsyr_               dsyr
#define    cher_               cher
#define    zher_               zher

#define    ssyr2_              ssyr2
#define    dsyr2_              dsyr2
#define    cher2_              cher2
#define    zher2_              zher2

#define    sgemm_              sgemm
#define    dgemm_              dgemm
#define    cgemm_              cgemm
#define    zgemm_              zgemm

#define    ssymm_              ssymm
#define    dsymm_              dsymm
#define    csymm_              csymm
#define    chemm_              chemm
#define    zsymm_              zsymm
#define    zhemm_              zhemm

#define    strmm_              strmm
#define    dtrmm_              dtrmm
#define    ctrmm_              ctrmm
#define    ztrmm_              ztrmm

#define    strsm_              strsm
#define    dtrsm_              dtrsm
#define    ctrsm_              ctrsm
#define    ztrsm_              ztrsm

#define    ssyrk_              ssyrk
#define    dsyrk_              dsyrk
#define    csyrk_              csyrk
#define    cherk_              cherk
#define    zsyrk_              zsyrk
#define    zherk_              zherk

#define    ssyr2k_             ssyr2k
#define    dsyr2k_             dsyr2k
#define    csyr2k_             csyr2k
#define    cher2k_             cher2k
#define    zsyr2k_             zsyr2k
#define    zher2k_             zher2k

#endif
/*
*  ---------------------------------------------------------------------
*  Function prototypes
*  ---------------------------------------------------------------------
*/
#ifdef __STDC__

Int            isamax_         ( Int *,     char *,    Int * );
Int            idamax_         ( Int *,     char *,    Int * );
Int            icamax_         ( Int *,     char *,    Int * );
Int            izamax_         ( Int *,     char *,    Int * );

F_VOID_FCT     saxpy_          ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     daxpy_          ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     caxpy_          ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zaxpy_          ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );

F_VOID_FCT     scopy_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     dcopy_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ccopy_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     zcopy_          ( Int *,     char *,    Int *,
                                 char *,    Int * );

F_VOID_FCT     sscal_          ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dscal_          ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     cscal_          ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     csscal_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zdscal_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zscal_          ( Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     sswap_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     dswap_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     cswap_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     zswap_          ( Int *,     char *,    Int *,
                                 char *,    Int * );

F_VOID_FCT     sgemv_          ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dgemv_          ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cgemv_          ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zgemv_          ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     ssymv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dsymv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     chemv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zhemv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     strmv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     dtrmv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ctrmv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ztrmv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );

F_VOID_FCT     strsv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     dtrsv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ctrsv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ztrsv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );

F_VOID_FCT     sger_           ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     dger_           ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     cgerc_          ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     cgeru_          ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zgerc_          ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zgeru_          ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );

F_VOID_FCT     ssyr_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     dsyr_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     cher_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     zher_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );

F_VOID_FCT     ssyr2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     dsyr2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     cher2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zher2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );

F_VOID_FCT     sgemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dgemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     cgemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zgemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     ssymm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     dsymm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     csymm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     zsymm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     chemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     zhemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );

F_VOID_FCT     ssyrk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dsyrk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     csyrk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zsyrk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     cherk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zherk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     ssyr2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     dsyr2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     csyr2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     zsyr2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     cher2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     zher2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );

F_VOID_FCT     strmm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     dtrmm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ctrmm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ztrmm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );

F_VOID_FCT     strsm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     dtrsm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ctrsm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ztrsm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );

#else

Int            isamax_         ();
Int            idamax_         ();
Int            icamax_         ();
Int            izamax_         ();

F_VOID_FCT     saxpy_          ();
F_VOID_FCT     daxpy_          ();
F_VOID_FCT     caxpy_          ();
F_VOID_FCT     zaxpy_          ();

F_VOID_FCT     scopy_          ();
F_VOID_FCT     dcopy_          ();
F_VOID_FCT     ccopy_          ();
F_VOID_FCT     zcopy_          ();

F_VOID_FCT     sscal_          ();
F_VOID_FCT     dscal_          ();
F_VOID_FCT     cscal_          ();
F_VOID_FCT     csscal_         ();
F_VOID_FCT     zscal_          ();
F_VOID_FCT     zdscal_         ();

F_VOID_FCT     sswap_          ();
F_VOID_FCT     dswap_          ();
F_VOID_FCT     cswap_          ();
F_VOID_FCT     zswap_          ();

F_VOID_FCT     sgemv_          ();
F_VOID_FCT     dgemv_          ();
F_VOID_FCT     cgemv_          ();
F_VOID_FCT     zgemv_          ();

F_VOID_FCT     ssymv_          ();
F_VOID_FCT     dsymv_          ();
F_VOID_FCT     chemv_          ();
F_VOID_FCT     zhemv_          ();

F_VOID_FCT     strmv_          ();
F_VOID_FCT     dtrmv_          ();
F_VOID_FCT     ctrmv_          ();
F_VOID_FCT     ztrmv_          ();

F_VOID_FCT     strsv_          ();
F_VOID_FCT     dtrsv_          ();
F_VOID_FCT     ctrsv_          ();
F_VOID_FCT     ztrsv_          ();

F_VOID_FCT     sger_           ();
F_VOID_FCT     dger_           ();
F_VOID_FCT     cgerc_          ();
F_VOID_FCT     cgeru_          ();
F_VOID_FCT     zgerc_          ();
F_VOID_FCT     zgeru_          ();

F_VOID_FCT     ssyr_           ();
F_VOID_FCT     dsyr_           ();
F_VOID_FCT     cher_           ();
F_VOID_FCT     zher_           ();

F_VOID_FCT     ssyr2_          ();
F_VOID_FCT     dsyr2_          ();
F_VOID_FCT     cher2_          ();
F_VOID_FCT     zher2_          ();

F_VOID_FCT     sgemm_          ();
F_VOID_FCT     dgemm_          ();
F_VOID_FCT     cgemm_          ();
F_VOID_FCT     zgemm_          ();

F_VOID_FCT     ssymm_          ();
F_VOID_FCT     dsymm_          ();
F_VOID_FCT     csymm_          ();
F_VOID_FCT     zsymm_          ();
F_VOID_FCT     chemm_          ();
F_VOID_FCT     zhemm_          ();

F_VOID_FCT     ssyrk_          ();
F_VOID_FCT     dsyrk_          ();
F_VOID_FCT     csyrk_          ();
F_VOID_FCT     zsyrk_          ();
F_VOID_FCT     cherk_          ();
F_VOID_FCT     zherk_          ();

F_VOID_FCT     ssyr2k_         ();
F_VOID_FCT     dsyr2k_         ();
F_VOID_FCT     csyr2k_         ();
F_VOID_FCT     zsyr2k_         ();
F_VOID_FCT     cher2k_         ();
F_VOID_FCT     zher2k_         ();

F_VOID_FCT     strmm_          ();
F_VOID_FCT     dtrmm_          ();
F_VOID_FCT     ctrmm_          ();
F_VOID_FCT     ztrmm_          ();

F_VOID_FCT     strsm_          ();
F_VOID_FCT     dtrsm_          ();
F_VOID_FCT     ctrsm_          ();
F_VOID_FCT     ztrsm_          ();

#endif
