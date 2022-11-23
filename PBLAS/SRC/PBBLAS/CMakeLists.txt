set (PBSBLASAUX  pbstran.f pbsmatadd.f pbstrsrt.f pbstrget.f
             pbstrnv.f pbsvecadd.f pbstrst1.f)

set (PBCBLASAUX  pbctran.f pbcmatadd.f pbctrsrt.f pbctrget.f
             pbctrnv.f pbcvecadd.f pbctrst1.f)

set (PBDBLASAUX  pbdtran.f pbdmatadd.f pbdtrsrt.f pbdtrget.f
             pbdtrnv.f pbdvecadd.f pbdtrst1.f)

set (PBZBLASAUX  pbztran.f pbzmatadd.f pbztrsrt.f pbztrget.f
             pbztrnv.f pbzvecadd.f pbztrst1.f)

set(pbblas)
if(BUILD_SINGLE)
  list(APPEND pbblas ${PBSBLASAUX})
endif()
if(BUILD_DOUBLE)
  list(APPEND pbblas ${PBDBLASAUX})
endif()
if(BUILD_COMPLEX)
  list(APPEND pbblas ${PBCBLASAUX})
endif()
if(BUILD_COMPLEX16)
  list(APPEND pbblas ${PBZBLASAUX})
endif()
