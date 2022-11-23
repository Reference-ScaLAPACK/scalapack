if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch;-std=legacy>")
endif()

file(COPY ../SRC/PTOOLS/PB_Cwarn.c DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
file(COPY ../SRC/PTOOLS/PB_Cabort.c DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set (PblasErrorHandler  PB_Cwarn.c PB_Cabort.c)
set (pbtcom pblastim.f ${PblasErrorHandler})
set (spbtcom psblastim.f ${pbtcom})
set (dpbtcom pdblastim.f ${pbtcom})
set (cpbtcom pcblastim.f ${pbtcom})
set (zpbtcom pzblastim.f ${pbtcom})

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
  _pblas_test(spb1tim "psblas1tim.f;${spbtcom}")
  _pblas_test(spb2tim "psblas2tim.f;${spbtcom}")
  _pblas_test(spb3tim "psblas3tim.f;${spbtcom}")
endif()
if(BUILD_DOUBLE)
  _pblas_test(dpb1tim "pdblas1tim.f;${dpbtcom}")
  _pblas_test(dpb2tim "pdblas2tim.f;${dpbtcom}")
  _pblas_test(dpb3tim "pdblas3tim.f;${dpbtcom}")
endif()
if(BUILD_COMPLEX)
  _pblas_test(cpb1tim "pcblas1tim.f;${cpbtcom}")
  _pblas_test(cpb2tim "pcblas2tim.f;${cpbtcom}")
  _pblas_test(cpb3tim "pcblas3tim.f;${cpbtcom}")
endif()
if(BUILD_SINGLE AND BUILD_COMPLEX16)
  _pblas_test(zpb1tim "pzblas1tim.f;${zpbtcom}")
  _pblas_test(zpb2tim "pzblas2tim.f;${zpbtcom}")
  _pblas_test(zpb3tim "pzblas3tim.f;${zpbtcom}")
endif()
