add_subdirectory(PBBLAS)
add_subdirectory(PTZBLAS)
add_subdirectory(PTOOLS)

set (APPBLAS  pilaenv.f)

#---------------------------------------------------------------------------
#  Level 1 PBLAS.
#---------------------------------------------------------------------------

set (PIBLAS1  picopy_.c)

set (PSBLAS1  psswap_.c psscal_.c  pscopy_.c  psaxpy_.c psdot_.c  psnrm2_.c 
          psasum_.c psamax_.c)

set (PCBLAS1  pcswap_.c pcscal_.c  pcsscal_.c pccopy_.c pcaxpy_.c pcdotu_.c 
          pcdotc_.c pscnrm2_.c pscasum_.c pcamax_.c)

set (PDBLAS1  pdswap_.c pdscal_.c  pdcopy_.c  pdaxpy_.c pddot_.c  pdnrm2_.c 
          pdasum_.c pdamax_.c)

set (PZBLAS1  pzswap_.c pzscal_.c  pzdscal_.c pzcopy_.c pzaxpy_.c pzdotu_.c 
          pzdotc_.c pdznrm2_.c pdzasum_.c pzamax_.c)

#---------------------------------------------------------------------------
#  Level 2 PBLAS.
#---------------------------------------------------------------------------

set (PSBLAS2  psgemv_.c  psger_.c   pssymv_.c  pssyr_.c   pssyr2_.c  pstrmv_.c 
          pstrsv_.c  psagemv_.c psasymv_.c psatrmv_.c)

set (PCBLAS2  pcgemv_.c  pcgerc_.c  pcgeru_.c  pchemv_.c  pcher_.c   pcher2_.c 
          pctrmv_.c  pctrsv_.c  pcagemv_.c pcahemv_.c pcatrmv_.c)

set (PDBLAS2  pdgemv_.c  pdger_.c   pdsymv_.c  pdsyr_.c   pdsyr2_.c  pdtrmv_.c 
          pdtrsv_.c  pdagemv_.c pdasymv_.c pdatrmv_.c)

set (PZBLAS2  pzgemv_.c  pzgerc_.c  pzgeru_.c  pzhemv_.c  pzher_.c   pzher2_.c 
          pztrmv_.c  pztrsv_.c  pzagemv_.c pzahemv_.c pzatrmv_.c)

#---------------------------------------------------------------------------
#  Level 3 PBLAS.
#---------------------------------------------------------------------------

set (PSBLAS3  psgeadd_.c psgemm_.c pssymm_.c  pssyr2k_.c pssyrk_.c  pstradd_.c 
          pstran_.c  pstrmm_.c pstrsm_.c)

set (PCBLAS3  pcgeadd_.c pcgemm_.c pchemm_.c  pcher2k_.c pcherk_.c  pcsymm_.c  
          pcsyr2k_.c pcsyrk_.c pctradd_.c pctranc_.c pctranu_.c pctrmm_.c  
          pctrsm_.c)

set (PDBLAS3  pdgeadd_.c pdgemm_.c pdsymm_.c  pdsyr2k_.c pdsyrk_.c  pdtradd_.c 
          pdtran_.c  pdtrmm_.c pdtrsm_.c)

set (PZBLAS3  pzgeadd_.c pzgemm_.c pzhemm_.c  pzher2k_.c pzherk_.c  pzsymm_.c  
          pzsyr2k_.c pzsyrk_.c pztradd_.c pztranc_.c pztranu_.c pztrmm_.c  
          pztrsm_.c)

#---------------------------------------------------------------------------

set (pblas 
      ${PIBLAS1}
      ${PSBLAS1} ${PCBLAS1} ${PDBLAS1} ${PZBLAS1}
      ${PSBLAS2} ${PCBLAS2} ${PDBLAS2} ${PZBLAS2}
      ${PSBLAS3} ${PCBLAS3} ${PDBLAS3} ${PZBLAS3})

set (pblas-F ${APPBLAS} )