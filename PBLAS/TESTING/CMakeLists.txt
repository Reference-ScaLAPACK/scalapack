set (pbtcom pblastst.f)
set (spbtcom psblastst.f slamch.f ${pbtcom})
set (dpbtcom pdblastst.f dlamch.f ${pbtcom})
set (cpbtcom pcblastst.f slamch.f ${pbtcom})
set (zpbtcom pzblastst.f dlamch.f ${pbtcom})

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
add_executable(spb2tst psblas2tst.f ${spbtcom})
add_executable(spb3tst psblas3tst.f ${spbtcom})

target_link_libraries(spb1tst pblas-test-support scalapack LAPACK::LAPACK MPI::MPI_Fortran)
target_link_libraries(spb2tst pblas-test-support scalapack LAPACK::LAPACK MPI::MPI_Fortran)
target_link_libraries(spb3tst pblas-test-support scalapack LAPACK::LAPACK MPI::MPI_Fortran)

add_test(spb1tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./spb1tst)
add_test(spb2tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./spb2tst)
add_test(spb3tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./spb3tst)

if(SCALAPACK_BUILD_DOUBLE)
    add_executable(dpb1tst pdblas1tst.f ${dpbtcom})
    add_executable(dpb2tst pdblas2tst.f ${dpbtcom})
    add_executable(dpb3tst pdblas3tst.f ${dpbtcom})

    target_link_libraries(dpb1tst pblas-test-support scalapack LAPACK::LAPACK MPI::MPI_Fortran)
    target_link_libraries(dpb2tst pblas-test-support scalapack LAPACK::LAPACK MPI::MPI_Fortran)
    target_link_libraries(dpb3tst pblas-test-support scalapack LAPACK::LAPACK MPI::MPI_Fortran)

    add_test(dpb1tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./dpb1tst)
    add_test(dpb2tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./dpb2tst)
    add_test(dpb3tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./dpb3tst)
endif()

if(SCALAPACK_BUILD_COMPLEX)
    add_executable(cpb1tst pcblas1tst.f ${cpbtcom})
    add_executable(cpb2tst pcblas2tst.f ${cpbtcom})
    add_executable(cpb3tst pcblas3tst.f ${cpbtcom})

    target_link_libraries(cpb1tst pblas-test-support scalapack LAPACK::LAPACK MPI::MPI_Fortran)
    target_link_libraries(cpb2tst pblas-test-support scalapack LAPACK::LAPACK MPI::MPI_Fortran)
    target_link_libraries(cpb3tst pblas-test-support scalapack LAPACK::LAPACK MPI::MPI_Fortran)

    add_test(cpb1tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./cpb1tst)
    add_test(cpb2tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./cpb2tst)
    add_test(cpb3tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./cpb3tst)
endif()

if(SCALAPACK_BUILD_COMPLEX16)
    add_executable(zpb1tst pzblas1tst.f ${zpbtcom})
    add_executable(zpb2tst pzblas2tst.f ${zpbtcom})
    add_executable(zpb3tst pzblas3tst.f ${zpbtcom})

    target_link_libraries(zpb1tst pblas-test-support scalapack LAPACK::LAPACK MPI::MPI_Fortran)
    target_link_libraries(zpb2tst pblas-test-support scalapack LAPACK::LAPACK MPI::MPI_Fortran)
    target_link_libraries(zpb3tst pblas-test-support scalapack LAPACK::LAPACK MPI::MPI_Fortran)

    add_test(zpb1tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./zpb1tst)
    add_test(zpb2tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./zpb2tst)
    add_test(zpb3tst ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./zpb3tst)
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy" )  # local to this directory
endif()
