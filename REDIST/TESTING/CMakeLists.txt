set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SCALAPACK_BINARY_DIR}/scalapack/REDIST/TESTING)

file(COPY GEMR2D.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(COPY TRMR2D.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

add_executable(xigemr pigemrdrv.c)
add_executable(xsgemr psgemrdrv.c)
add_executable(xdgemr pdgemrdrv.c)
add_executable(xcgemr pcgemrdrv.c)
add_executable(xzgemr pzgemrdrv.c)

add_executable(xitrmr pitrmrdrv.c)
add_executable(xstrmr pstrmrdrv.c)
add_executable(xdtrmr pdtrmrdrv.c)
add_executable(xctrmr pctrmrdrv.c)
add_executable(xztrmr pztrmrdrv.c)

target_link_libraries(xigemr scalapack MPI::MPI_Fortran )
target_link_libraries(xsgemr scalapack MPI::MPI_Fortran )
target_link_libraries(xdgemr scalapack MPI::MPI_Fortran )
target_link_libraries(xcgemr scalapack MPI::MPI_Fortran )
target_link_libraries(xzgemr scalapack MPI::MPI_Fortran )

target_link_libraries(xitrmr scalapack MPI::MPI_Fortran )
target_link_libraries(xstrmr scalapack MPI::MPI_Fortran )
target_link_libraries(xdtrmr scalapack MPI::MPI_Fortran )
target_link_libraries(xctrmr scalapack MPI::MPI_Fortran )
target_link_libraries(xztrmr scalapack MPI::MPI_Fortran )

#add_test(xigemr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xigemr)
#add_test(xsgemr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xsgemr)
#add_test(xdgemr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xdgemr)
#add_test(xcgemr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xcgemr)
#add_test(xzgemr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xzgemr)

#add_test(xitrmr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xitrmr)
#add_test(xstrmr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xstrmr)
#add_test(xdtrmr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xdtrmr)
#add_test(xctrmr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xctrmr)
#add_test(xztrmr ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ./xztrmr)

