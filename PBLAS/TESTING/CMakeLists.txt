file(COPY ../SRC/PTOOLS/PB_Cwarn.c DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
file(COPY ../SRC/PTOOLS/PB_Cabort.c DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})

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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SCALAPACK_BINARY_DIR}/PBLAS/TESTING)

file(COPY PCBLAS1TST.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PCBLAS2TST.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PCBLAS3TST.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PDBLAS1TST.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PDBLAS2TST.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PDBLAS3TST.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PSBLAS1TST.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PSBLAS2TST.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PSBLAS3TST.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PZBLAS1TST.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PZBLAS2TST.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY PZBLAS3TST.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})


add_executable(spb1tst psblas1tst.f ${spbtcom})
add_executable(dpb1tst pdblas1tst.f ${dpbtcom})
add_executable(cpb1tst pcblas1tst.f ${cpbtcom})
add_executable(zpb1tst pzblas1tst.f ${zpbtcom})

add_executable(spb2tst psblas2tst.f ${spbtcom})
add_executable(dpb2tst pdblas2tst.f ${dpbtcom})
add_executable(cpb2tst pcblas2tst.f ${cpbtcom})
add_executable(zpb2tst pzblas2tst.f ${zpbtcom})

add_executable(spb3tst psblas3tst.f ${spbtcom})
add_executable(dpb3tst pdblas3tst.f ${dpbtcom})
add_executable(cpb3tst pcblas3tst.f ${cpbtcom})
add_executable(zpb3tst pzblas3tst.f ${zpbtcom})

target_link_libraries(spb1tst scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)
target_link_libraries(dpb1tst scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)
target_link_libraries(cpb1tst scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)
target_link_libraries(zpb1tst scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)

target_link_libraries(spb2tst scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)
target_link_libraries(dpb2tst scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)
target_link_libraries(cpb2tst scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)
target_link_libraries(zpb2tst scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)

target_link_libraries(spb3tst scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)
target_link_libraries(dpb3tst scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)
target_link_libraries(cpb3tst scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)
target_link_libraries(zpb3tst scalapack ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} MPI::MPI_Fortran)

add_test(spb1tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./spb1tst)
add_test(dpb1tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./dpb1tst)
add_test(cpb1tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./cpb1tst)
add_test(zpb1tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./zpb1tst)

add_test(spb2tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./spb2tst)
add_test(dpb2tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./dpb2tst)
add_test(cpb2tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./cpb2tst)
add_test(zpb2tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./zpb2tst)

add_test(spb3tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./spb3tst)
add_test(dpb3tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./dpb3tst)
add_test(cpb3tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./cpb3tst)
add_test(zpb3tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./zpb3tst)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy" )  # local to this directory
endif()
