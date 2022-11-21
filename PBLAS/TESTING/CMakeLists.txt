file(COPY ../SRC/PTOOLS/PB_Cwarn.c DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
file(COPY ../SRC/PTOOLS/PB_Cabort.c DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set (PblasErrorHandler  PB_Cwarn.c PB_Cabort.c)
set (pbtcom pblastst.f ${PblasErrorHandler})
set (spbtcom psblastst.f slamch.f ${pbtcom})
set (dpbtcom pdblastst.f dlamch.f ${pbtcom})
set (cpbtcom pcblastst.f slamch.f ${pbtcom})
set (zpbtcom pzblastst.f dlamch.f ${pbtcom})

set_property(
   SOURCE ${PblasErrorHandler}
   APPEND PROPERTY COMPILE_DEFINITIONS TestingPblas
   )

function(_pblas_test name src)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)
  add_test(NAME ${name}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${name}>
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endfunction()

if(BUILD_SINGLE)
  _pblas_test(spb1tst "psblas1tst.f;${spbtcom}")
  _pblas_test(spb2tst "psblas2tst.f;${spbtcom}")
  _pblas_test(spb3tst "psblas3tst.f;${spbtcom}")
endif()
if(BUILD_SINGLE AND BUILD_DOUBLE)
  _pblas_test(dpb1tst "pdblas1tst.f;${dpbtcom}")
  _pblas_test(dpb2tst "pdblas2tst.f;${dpbtcom}")
  _pblas_test(dpb3tst "pdblas3tst.f;${dpbtcom}")
endif()
if(BUILD_COMPLEX)
  _pblas_test(cpb1tst "pcblas1tst.f;${cpbtcom}")
  _pblas_test(cpb2tst "pcblas2tst.f;${cpbtcom}")
  _pblas_test(cpb3tst "pcblas3tst.f;${cpbtcom}")
endif()
if(BUILD_SINGLE AND BUILD_COMPLEX16)
  _pblas_test(zpb1tst "pzblas1tst.f;${zpbtcom}")
  _pblas_test(zpb2tst "pzblas2tst.f;${zpbtcom}")
  _pblas_test(zpb3tst "pzblas3tst.f;${zpbtcom}")
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-std=legacy>)  # local to this directory
endif()
